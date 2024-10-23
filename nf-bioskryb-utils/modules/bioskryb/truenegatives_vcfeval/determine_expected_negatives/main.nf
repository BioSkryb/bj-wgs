nextflow.enable.dsl=2
params.timestamp = ""

process DETERMINE_EXPECTED_NEGATIVES {
    tag "DETERMINE_EXPECTED_NEGATIVES"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path(baseline_vcf)
    path(baseline_vcf_index)
    path(baseline_regions)
    path(bed_regions)
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("expected_0_0.bed"), emit: expected_inert

    script:
    """
    #! /bin/bash
    set +u
    
    bedtools subtract -a ${baseline_regions} -b ${baseline_vcf}  > temp.bed
    
    bedtools intersect -a temp.bed -b ${bed_regions} > expected_0_0.bed 


    """
    
}