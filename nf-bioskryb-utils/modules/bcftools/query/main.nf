nextflow.enable.dsl=2
params.timestamp = ""

process BCFTOOLS_QUERY_L {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(input_vcf)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path(input_vcf), emit: vcf
    tuple val(sample_name), path("*id"), emit: list_names
    path("bcftools_version.yml"), emit: version
    
    script:
    """
    
    bcftools query -l ${input_vcf[0]} > list_names_${sample_name}.txt
    
    cat list_names_${sample_name}.txt | while read id;
    
    do
    
        echo -e \${id}  > \${id}.id
        
    done
    
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo  BCFtools: \$BCFTOOLS_VER > bcftools_version.yml
    
    """
    
}