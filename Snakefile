ref_gtf='ref/gencode_annotation_v37.gtf'
ref_genome='ref/gencode_genome_v37.fa'
chr_names= '/data/swamyvs/pacbio_testing/ref/chr_fasta_entries.txt'
variants= '/data/swamyvs/pacbio_testing/ref/encode_variants.vcf.gz'

rule all:
    input: 
        expand('data/talon_results/{sample}/{sample}_talon_observedOnly.gtf',sample =['ONT_RNA_RPE_D42_all', 'pacbio_RNA_RPE_D42_all','bonito-ont_RNA_RPE_D42_all']), 
        expand('data/fasta_lengths/{sample}.tsv',sample =['ONT_RNA_RPE_D42_all', 'pacbio_RNA_RPE_D42_all','bonito-ont_RNA_RPE_D42_all']),
        'data/combined_gtfs/ipscRPE-pacbio-ONT_talon_vs_st.combined.gtf'
        
       
##Note: this was not run through Snakemake
rule bonito_basecall:
    output:'fastq/bonito-ont_RNA_RPE_D42_all.fastq.gz'
    shell:
        '''
        module load bonito/0.3.6
        bonito basecaller --fastq dna_r9.4.1 /data/OGVFB_BG/iPSC_RPE_ONT_RNA-seq/fast5_pass/ | gzip -c - > fastq/bonito-ont_RNA_RPE_D42_all.fastq.gz
        '''

#### build transcriptomes 

rule get_annotation:
    output:
        gtf = ref_gtf,
        genome = ref_genome
    shell:
        '''
        wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz | gunzip -c - > {output.gtf}


        wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz | gunzip -c - > {output.genome}
        
        '''



#####long read stuff
rule clean_genome:
    input: ref_genome
    output: 'ref/gencode_genome_clean.fa'
    shell:
        '''
        python3 scripts/filterFasta.py {input} {chr_names} {output}
        '''



rule build_minimap_index:
    input: 
        'ref/gencode_genome_clean.fa'
    output: 
        'data/index/minimap2_index.mmi'
    shell:
        '''
        module load minimap2/2.18
        minimap2 -d {output} {input}
        '''


rule build_talon_db:
    input: 
        gtf=ref_gtf, 
        genome = 'ref/gencode_genome_clean.fa'
    params: 
        anno_name= lambda wildcards: f'db_{wildcards.sample}', 
        prefix= lambda wildcards: f'{wildcards.sample}', 
        outdir= lambda wildcards: f'data/talon_db/{wildcards.sample}'
    output: 
        'data/talon_db/{sample}.db'
    shell:
        '''

        talon_initialize_database \
            --f {input.gtf} --g {input.genome} \
            --a {params.anno_name} \
            --5p 100 \
            --3p 100 \
            --idprefix {params.prefix} \
            --o {params.outdir} 

        '''

rule calc_fasta_lengths:
    input:
        fa='fastq/{sample}.fastq.gz'
    output:
        tab='data/fasta_lengths/{sample}.tsv'
    shell:
        '''
        zcat {input.fa} > /tmp/{wildcards.sample}.fq
        python3 scripts/get_fasta_lengths.py --inFasta /tmp/{wildcards.sample}.fq --outFile {output.tab}
        '''


rule align_minimap2:
    input:
        fa='fastq/{sample}.fastq.gz' ,
        idx='data/index/minimap2_index.mmi', 
        genome = 'ref/gencode_genome_clean.fa'
    output:
        sam='data/sams/{sample}.sam'
    shell:
        '''
        module load minimap2/2.18
        module load samtools/1.11
        minimap2 -t 8 -a -x splice -u f --MD  {input.genome} {input.fa} | samtools sort -O sam - > {output.sam}
        
        '''
        
'''
JK build is the genome
'''

rule deprime_sams:
    input:
        sam = 'data/sams/{sample}.sam'
    output:
        sam = 'data/sams/{sample}_labeled.sam'
    params:
        prefix = lambda wildcards: f'data/sams/{wildcards.sample}'
    shell:
        '''
        talon_label_reads --f {input.sam} --g {ref_genome} --o {params.prefix}
        '''

# rule make_spliceJN_file:
#     input:
#         gtf=ref_gtf, 
#         genome = 'ref/gencode_genome_clean.fa'
#     output:
#         'data/SJ_tab.txt'
#     shell:
#         '''
#         module load samtools/1.11
#         module load  bedtools/2.30.0 
#         python TranscriptClean/accessory_scripts/get_SJs_from_gtf.py --f {input.gtf} --g {input.genome} --o {output}

#         '''


# rule TranscriptClean: # pull from github repo
#     input: 
#         sam = 'data/sams/{sample}_labeled.sam', 
#         genome = 'ref/gencode_genome_clean.fa'
#     params:
#         out_pref=lambda wildcards: f'data/sams/{wildcards.sample}',
#         tmp_dir=lambda wildcards: f'data/tmp/{wildcards.sample}/'
#     output:
#         sam='data/sams/{sample}_clean.sam',
#         bam='data/bams/{sample}.bam'
#     shell:
#         '''
#         rm -rf {params.tmp_dir}
#         module load samtools/1.11
#         module load bedtools/2.30.0 
#         /data/swamyvs/anaconda3/envs/ipsc_lr/bin/python TranscriptClean/TranscriptClean.py --sam {input.sam} \
#             --genome {input.genome}  \
#             --variants {variants}\
#             --maxLenIndel 5 \
#             --maxSJOffset 5 \
#             --correctMismatches true \
#             --correctIndels true \
#             --primaryOnly \
#             --outprefix {params.out_pref} \
#             --tmpDir {params.tmp_dir} \
#             --deleteTmp \
#             --threads 32 
#         samtools view -S -b {output.sam} > {output.bam}
#         '''





rule run_talon:
    input: 
        sam = 'data/sams/{sample}_labeled.sam',
        db = 'data/talon_db/{sample}.db', 
        genome = 'ref/gencode_genome_clean.fa'
    params:
        talon_config= lambda wildcards: f'{wildcards.sample},dummy1,dummy2,/data/swamyvs/ipsc_RPE_longread_analysis/data/sams/{wildcards.sample}_labeled.sam',
        build='gencode', 
        res_outdir=lambda wildcards: f'data/talon_results/{wildcards.sample}/',
        prefix= lambda wildcards: f'db_{wildcards.sample}', 
        outdir_pf=lambda wildcards: f'data/talon_results/{wildcards.sample}/{wildcards.sample}'
    output: 
        anno = 'data/talon_results/{sample}/_talon_read_annot.tsv', 
        qc = 'data/talon_results/{sample}/talon_QC.log',
        gtf='data/talon_results/{sample}/{sample}_talon_observedOnly.gtf',
        abundance = 'data/talon_results/{sample}/{sample}_talon_abundance.tsv'
    shell:
        '''
        tlncfg=/tmp/config.{wildcards.sample}.csv
        echo "{params.talon_config}" >  $tlncfg
        fp=$PWD
        mkdir -p {params.res_outdir}
        cd {params.res_outdir}
        talon \
            --f $tlncfg \
            --db $fp/{input.db} \
            --build {input.genome} \
            --threads 16 \
            --o $fp/{params.res_outdir}
        
        talon_create_GTF \
            --db $fp/{input.db} \
            --build {input.genome} \
            --annot {params.prefix} \
            --observed \
            --o $fp/{params.outdir_pf}
        
        talon_abundance \
            --db $fp/{input.db} \
            --build {input.genome} \
            --annot {params.prefix} \
            --o $fp/{params.outdir_pf}

        '''


rule run_stringtie:
    input: 
        gtf = ref_gtf,
        bam='data/bams/{sample}.bam'
    output: 
        'data/stringtie/{sample}.gtf'
    params:
        tag = lambda wildcards: wildcards.sample.split('_')[0]
    shell:
        '''
        module load stringtie/2.1.5
        stringtie -L -G {input.gtf} -p 8 -l {params.tag} -o {output} {input.bam}
        '''

rule merge_all_gtfs:
    input:  
        talon_gtfs = expand('data/talon_results/{sample}/{sample}_talon_observedOnly.gtf',sample =[ 'pacbio_RNA_RPE_D42_all', 'bonito-ont_RNA_RPE_D42_all'] ),
        string_tie_gtfs = expand('data/stringtie/{sample}.gtf',sample =['pacbio_RNA_RPE_D42_all', 'bonito-ont_RNA_RPE_D42_all'] )
    params:
        merge_prefix =lambda wildcards: f'data/combined_gtfs/ipscRPE-pacbio-ONT_talon_vs_st',
        label_prefix = 'ipscRPE'
    output: 
        all_merge_gtf = 'data/combined_gtfs/ipscRPE-pacbio-ONT_talon_vs_st.combined.gtf'
    shell:
        '''
        module load gffcompare/0.11.8
        gffcompare --strict-match -e 0 -d 0  -r {ref_gtf} -p {params.label_prefix} -o  {params.merge_prefix} {input.talon_gtfs} {input.string_tie_gtfs}
        '''

'''
Run this to merge this long read data with old short read data
module load gffcompare/0.11.8
 k=`ls ../ocular_transcriptomes_longread_analysis/data/stringtie/H*.gtf`
 gffcompare --strict-match -e 0 -d 0  -r ref/gencode_annotation_v37.gtf  -p ipscRPE -o  data/combined_gtfs/lr_sr_merged  \
  'data/talon_results/pacbio_RNA_RPE_D42_all/pacbio_RNA_RPE_D42_all_talon_observedOnly.gtf' \
  'data/talon_results/bonito-ont_RNA_RPE_D42_all/bonito-ont_RNA_RPE_D42_all_talon_observedOnly.gtf' \
  'data/stringtie/pacbio_RNA_RPE_D42_all.gtf' \
  'data/stringtie/bonito-ont_RNA_RPE_D42_all.gtf' $k
'''