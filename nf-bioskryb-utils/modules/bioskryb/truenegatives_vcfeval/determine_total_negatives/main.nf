nextflow.enable.dsl=2
params.timestamp = ""

process DETERMINE_TOTAL_NEGATIVES {
    tag "DETERMINE_TOTAL_NEGATIVES_${sample_name}_${id_sample}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(input_vcf), path(id_sample),  path(expected_inert)
    val(publish_dir)
    val(enable_publish)


  
    
    output:
    path("df_eval_negatives_${sample_name}_${id_sample}.tsv"), emit: df_eval_neg

    script:
    """
    #! /bin/bash
    set +u
    
    bcftools view -f 'PASS,.' --types snps  ${input_vcf[0]}  | \
    bedtools intersect -a ${expected_inert} -b - -wb > res.txt
    
    total_bases=`cat expected_0_0.bed | awk '{print \$3-\$2}'| awk '{sum+=\$1}END{print sum}'`;                                                     
    missed_neg=`wc -l res.txt | awk '{print \$1}'`;
    
    echo -e "${sample_name}\\t${id_sample}\\t\${total_bases}\\t\${missed_neg}" > df_eval_negatives_${sample_name}_${id_sample}.tsv
    
    
    """
    
}