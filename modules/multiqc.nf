process MULTIQC {
    //conda 'multiqc'
    conda "/export/home/agletdinov/mambaforge/envs/multiqc"
    
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode:'copy'

    input:
    path("*")

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}