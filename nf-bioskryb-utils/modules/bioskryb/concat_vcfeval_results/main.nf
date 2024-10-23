nextflow.enable.dsl=2
params.timestamp = ""

process CONCAT_VCFEVAL_RESULTS {
    tag "CONCAT_VCFEVAL_RESULTS"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/vcfeval_results", enabled:"$enable_publish"
    
    
    input:
    path(res_score)
    val(publish_dir)
    val(enable_publish)

    output:
    path("vcfeval_results_mqc.tsv")
    
    script:
    """
    
    echo -e "sample_name\\tscore\\ttrue_positives_baseline\\tfalse_positives\\ttrue_positives_call\\tfalse_negatives\\tprecision\\tsensitivity\\tf_measure" > vcfeval_results_mqc.tsv
    
    find . -name "*" | cut -d "/" -f2 | grep "snp_roc"  | grep -v "_non_" | while read file;
    
    do
    
        # Extract SAMPLE_NAME from the file name using the 'sed' command
        SAMPLE_NAME=\$(echo "\${file}" | sed -E 's/res_vcfeval_(.+)_snp_roc\\.tsv\\.gz/\\1/')

        zcat \${file} | tail -n 1 | sed "s|^|\${SAMPLE_NAME}\\t|" >> vcfeval_results_mqc.tsv
        
    done
    
    """
    
}