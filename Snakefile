ref_gtf='ref/gencode_annotation_v37.gtf'
ref_genome='ref/gencode_genome_v37.fa'
chr_names= '/data/swamyvs/pacbio_testing/ref/chr_fasta_entries.txt'


rule all:
    input: 
        'data/combined_gtfs/ipscRPE-pacbio-ONT_talon_vs_st.combined.gtf',
        expand('data/fasta_lengths/{sample}.tsv',sample =['ONT_RNA_RPE_D42_all', 'pacbio_RNA_RPE_D42_all'])
        
       

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
        sam='data/sams/{sample}.sam',
        bam='data/bams/{sample}.bam'
    shell:
        '''
        module load minimap2/2.18
        module load samtools/1.11
        minimap2 -t 8 -a -x splice -u f --MD  {input.genome} {input.fa} | samtools sort -O sam - > {output.sam}
        samtools view -S -b {output.sam} > {output.bam}
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
            --threads 8 \
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
    shell:
        '''
        module load stringtie/2.1.5
        stringtie -L -G {input.gtf} -o {output} {input.bam}
        '''

rule merge_all_gtfs:
    input:  
        talon_gtfs = expand('data/talon_results/{sample}/{sample}_talon_observedOnly.gtf',sample =['ONT_RNA_RPE_D42_all', 'pacbio_RNA_RPE_D42_all'] ),
        string_tie_gtfs = expand('data/stringtie/{sample}.gtf',sample =['ONT_RNA_RPE_D42_all', 'pacbio_RNA_RPE_D42_all'] )
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
