nextflow.enable.dsl=2
params.timestamp = ""

process VCFEVAL {
    tag "VCFEVAL_${sample_name}_${id_sample}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/vcfeval_results/${sample_name}", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(input_vcf), path(id_sample)
    path(baseline_vcf)
    path(baseline_vcf_index)
    path(baseline_regions)
    path(bed_regions)
    path(reference_sdf_ref)
    val(output_mode)
    val(reference_name)
    val(score_metric)
    val(other_options)
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("res_vcfeval_*"), emit: results
    path("*_snp_roc.tsv.gz"), emit: metrics
    tuple val(sample_name), path("*_output.vcf.gz"), emit: df_class
    path("vcfeval_version.yml"), emit: version
    
    script:
    """
    #! /bin/bash
    set +u

    NAME=\$(echo -e ${id_sample} | sed 's/.id//')
    SAMPLE_STRING=\$(echo -e "${reference_name},\${NAME}")
    OPTIONS=\$(echo -e "${other_options} --vcf-score-field ${score_metric}")
    
    #export RTG_JAVA_OPTS='-Xms500m'
    #export RTG_MEM='${task.memory.toGiga()}g'
        
    rtg RTG_MEM=${task.memory.toGiga()}g vcfeval \
        --baseline ${baseline_vcf} \
        --threads $task.cpus \
        --calls ${input_vcf[0]} \
        --evaluation-regions ${baseline_regions} \
        --bed-regions ${bed_regions} \
        --output  "vcfeval_results_${score_metric}_\${NAME}/" \
        --template ${reference_sdf_ref} \
        \${OPTIONS} \
        --sample=\${SAMPLE_STRING} \
        --output-mode ${output_mode}

    cd "vcfeval_results_${score_metric}_\${NAME}/"
    
    for file in *;do mv \${file} "res_vcfeval_\${NAME}_\${file}";done
    
    mv res_vcfeval_* ../
    
    cd ../

    export VCFEVAL_VER=\$(rtg -v 2>&1 | head -n1 | sed 's/^.*Product: RTG Tools //; s/ .*\$//')
    echo VCFeval: \$VCFEVAL_VER > vcfeval_version.yml

    """
    
}