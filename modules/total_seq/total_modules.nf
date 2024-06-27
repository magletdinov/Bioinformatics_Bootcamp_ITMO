include { FASTQC as FQ1 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_2_NUCL; VSEARCH; MERGE_UNPAIRED; KRAKEN2; BRACKEN; BWA_INDEX; BWA_MEM_BAM_SORT; QUALIMAP; MEGAHIT; QUAST } from '../modules.nf'


workflow TAXONOMY_ANALYSIS {
  take:
    read_pairs_ch
    genome
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_2_NUCL(TRIM_ADAPT.out.paired, TRIM_ADAPT.out.unpaired)
    VSEARCH(TRIM_2_NUCL.out.paired)
    MERGE_UNPAIRED(VSEARCH.out.merged.join(TRIM_2_NUCL.out.unpaired))
    KRAKEN2(VSEARCH.out.unmerged.join(MERGE_UNPAIRED.out))
    BRACKEN(KRAKEN2.out.id_combined_report)
    MEGAHIT(VSEARCH.out.unmerged.join(MERGE_UNPAIRED.out))
    QUAST(MEGAHIT.out.id_contigs)
    //BWA_INDEX(genome, "homo")
    //BWA_MEM_BAM_SORT(MERGE_UNPAIRED.out, BWA_INDEX.out)
    //QUALIMAP(BWA_MEM_BAM_SORT.out)
    
  emit: 
     //FQ1.out | concat(KRAKEN2.out.combined_report) | collect
     FQ1.out | concat(BRACKEN.out) | collect
}