process FASTQC {
    //conda = 'bioconda::fastqc'
    conda "/export/home/agletdinov/mambaforge/envs/multiqc"

    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*")

    script:
    """
    fastqc ${reads}
    multiqc . --filename ${sample_id}
    """
}

process TRIM_ADAPT {
    conda = 'bioconda::fastp=0.23.4'

    tag "Fastp on ${sample_id}"
    publishDir "${params.outdir}/fastp_adapt", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*R*fastq.gz'), emit: paired
    tuple val(sample_id), path('*u.fastq.gz'), emit: unpaired
    
    cpus 20
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}-t_a-R1.fastq.gz -O  ${sample_id}-t_a-R2.fastq.gz --unpaired1 ${sample_id}-t_a-u.fastq.gz --unpaired2 ${sample_id}-t_a-u.fastq.gz --thread ${task.cpus} --detect_adapter_for_pe 
    """
}

process TRIM_2_NUCL {
    conda "/export/home/agletdinov/mambaforge/envs/cutadapt/"

    tag "Cutadapt on ${sample_id}"
    //publishDir "${params.outdir}/cutadapt_trim/${mode}_n", mode: "copy"
    publishDir "${params.outdir}/cutadapt_trim", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads_p)
    tuple val(sample_id), path(reads_up)
    //each mode

    output:
    tuple val("${sample_id}"), path('*R*fastq.gz'), emit: paired
    tuple val("${sample_id}"), path('*u.fastq.gz'), emit: unpaired
    cpus 20

    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5

    script:
    """
    cutadapt -u 2 -u -2 -U 2 -U -2 -o ${sample_id}-t_2-R1.fastq.gz -p ${sample_id}-t_2-R2.fastq.gz \
        ${reads_p[0]} ${reads_p[1]} \
        -j ${task.cpus} && \
    cutadapt -u 2 -u -2 -o ${sample_id}-t_2-u.fastq.gz ${reads_up} 
    """
}

process VSEARCH {
    conda = "/export/home/agletdinov/mambaforge/envs/vsearch"

    tag "VSEARCH on ${sample_id}"
    publishDir "${params.outdir}/vsearch", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*-unmerged*fastq.gz'), emit: unmerged
    tuple val(sample_id), path('*-merged*fastq.gz'), emit: merged
    
    cpus 20
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    script:
    """
    vsearch --fastq_mergepairs ${reads[0]} --reverse ${reads[1]} --fastq_minovlen 15 --fastq_maxdiffs 3 \
    --fastqout >(gzip -c > ${sample_id}-merged.fastq.gz) --fastqout_notmerged_fwd >(gzip -c > ${sample_id}-unmerged-R1.fastq.gz) --fastqout_notmerged_rev >(gzip -c > ${sample_id}-unmerged-R2.fastq.gz) --threads ${task.cpus}
    """
}

process MERGE_UNPAIRED {
    tag "cat on ${sample_id}"
    publishDir "${params.outdir}/cat", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads_v), path(reads_t2)
    
    output:
    tuple val(sample_id), path('*-merged-u.fastq.gz')
    
    cpus 20
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    script:
    """
    cat ${reads_v} ${reads_t2} > ${sample_id}-merged-u.fastq.gz
    """
}

process KRAKEN2 {
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    cpus 40
    maxForks 4

    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads_um), path(reads_m)
    
    output:
    path('*_combined_report.kraken'), emit: combined_report
    tuple val(sample_id), path('*'), emit: full_report
    tuple val(sample_id), path('*_paired_report.kraken'), emit: id_paired_report
    tuple val(sample_id), path('*_unpaired_report.kraken'), emit: id_unpaired_report
    tuple val(sample_id), path('*_combined_report.kraken'), emit: id_combined_report
    
    script:
    """
    # Run Kraken2 for paired reads
    kraken2 --db ${params.kraken2db} --threads ${task.cpus} --gzip-compressed --report ${sample_id}_paired_report.kraken --paired ${reads_um[0]} ${reads_um[1]} > ${sample_id}_paired_kraken.txt
    
    # Run Kraken2 for unpaired read
    kraken2 --db ${params.kraken2db} --threads ${task.cpus} --gzip-compressed --report ${sample_id}_unpaired_report.kraken ${reads_m} > ${sample_id}_unpaired_kraken.txt
    
    python3 /export/home/agletdinov/mambaforge/envs/kraken2/bin/combine_kreports.py -r ${sample_id}_paired_report.kraken ${sample_id}_unpaired_report.kraken -o ${sample_id}_combined_report.kraken --only-combined
    """
}


process BRACKEN {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/bracken"
    //maxForks 1
    cpus 20

    tag "Bracken on ${sample_id}"
    publishDir "${params.outdir}/bracken", mode: "copy"
    
    input:
    tuple val(sample_id), path(kraken_report)
    
    output:
    path('*.tsv')
    
    script:
    """
    bracken -l S -t ${task.cpus} -d ${params.kraken2db} -i ${kraken_report} -o ${sample_id}_bracken_S_mqc.tsv
    """
}


process BWA_INDEX {
    conda = '/export/home/agletdinov/mambaforge/envs/bwa'

    tag "Bwa index for ${params.genome}"
    publishDir "${params.outdir}/bwa_index", mode:'copy'
    
    input:
    path(ref_fasta)
    val(prefix)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(prefix), path("${prefix}.{ann,amb,sa,bwt,pac}")

    """
    bwa index \\
        -p "${prefix}" \\
        "${ref_fasta}"
    """
}


process BWA_MEM_BAM_SORT {
    conda '/export/home/agletdinov/mambaforge/envs/bwa'
    
    tag "Bwa mem on ${sample_id}"
    publishDir "${params.outdir}/bam", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)
    tuple val(idxbase), path("bwa_index/*")
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("${sample_id}.aln.sorted.bam")
   
    cpus 20
    //maxForks params.maxForks
    
    script:
    """
    bwa mem -t ${task.cpus} "bwa_index/${idxbase}" ${reads[0]} ${reads[1]}|
    samtools view --threads ${task.cpus} -1|
    samtools sort --threads ${task.cpus} -o "${sample_id}.aln.sorted.bam" 
    """
}


process QUALIMAP {
    conda "/export/home/agletdinov/mambaforge/envs/qualimap"

    tag "Qualimap on ${sample_id}"
    publishDir "${params.outdir}/qualimap", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam_file)

    output:
    path("${sample_id}_qualimap_report")
    
    cpus 20
    script:
    """
    qualimap bamqc -bam ${bam_file} -outdir ${sample_id} --java-mem-size=20G
    """
}


process MEGAHIT {
    //conda 'bioconda::megahit'
    conda "/export/home/agletdinov/mambaforge/envs/megahit"
    //memory 500.GB
    //maxForks 2
    cpus 40

    tag "Megahit on ${sample_id}"
    publishDir "${params.outdir}/megahit", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path('*'), emit: report
    //tuple val(sample_id), path('*'), emit: id_contigs
    tuple val("${sample_id}_megahit"), path("${sample_id}/*fa"), emit: id_contigs

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \
        -o ${sample_id} --out-prefix ${sample_id} -t ${task.cpus}
    """
}

process METAPHLAN {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/metaphlan"
    //memory = '1 MB'
    //maxForks 2
    cpus 40
    tag "Metaphlan on ${sample_id}"
    publishDir "${params.outdir}/metaphlan/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path('*.txt'), emit: txt_report
    tuple val(sample_id), path('*'), emit: id_report
    path('*.txt'), emit: report

    script:
    """
    metaphlan ${contigs} --bowtie2db ${params.metaphlandb} --bowtie2out metagenome_${sample_id}.bowtie2.bz2 --nproc ${task.cpus} --input_type fasta -o profiled_metagenome_${sample_id}.txt
    """
}

process DIAMOND {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/diamond"
    //memory = '1 MB'
    //maxForks 2
    cpus 20
    tag "Diamond on ${sample_id}"
    publishDir "${params.outdir}/diamond/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)

    output:
    path('*.tsv')
 
    script:
    """
    diamond blastx \
        --very-sensitive \
        -d /export/home/public/tools/database/virus_nr_diamond.dmnd \
        --outfmt 6 qseqid sseqid pident evalue \
        -p ${task.cpus} \
        -q ${contigs} \
        -o ${sample_id}_matches.tsv \
        --max-target-seqs 10 \
        --evalue 1e-05 \
        -v \
        --log
    """
}