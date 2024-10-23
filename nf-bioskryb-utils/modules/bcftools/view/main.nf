nextflow.enable.dsl=2
params.timestamp = ""

process BCFTOOLS_VIEW {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(input_vcf)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_filtered_wt.vcf.gz*"), emit: filtered_wt_vcf
    
    script:
    """
    bcftools view -e 'GT="0/0"' -Oz -o ${sample_name}_filtered_wt.vcf.gz ${input_vcf[0]}

    bcftools index -t ${sample_name}_filtered_wt.vcf.gz

    """
        
}
