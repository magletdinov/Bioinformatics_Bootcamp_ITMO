#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */

params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/data/agletdinov/hackathon_itmo/WGS_project"
params.reads = "${params.results_project}/raw_data/*/*R{1,2}*.fastq.gz"
params.genome = "${params.shared}/genomes/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
//params.kraken2db = "/export/data/roev/st"
params.nodes = "${projectDir}/ncbi_taxonomy/nodes.dmp"
params.outdir = "${params.results_project}/results"

params.script1 = "${projectDir}/bin/bracken_parse.py"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    script1        : ${params.script1}
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/total_seq/total_modules.nf'
include { BRACKEN_PARSE } from './modules/total_seq/bracken_parse.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow taxonomy_analysis{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    genome = params.genome
    script1 = params.script1
    TAXONOMY_ANALYSIS(read_pairs_ch, genome)
    BRACKEN_PARSE(script1, TAXONOMY_ANALYSIS.out)
    MULTIQC(TAXONOMY_ANALYSIS.out | concat(BRACKEN_PARSE.out) | collect)
}
//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//

//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//