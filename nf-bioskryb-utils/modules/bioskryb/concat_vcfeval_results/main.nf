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
    path("vcfeval_results_mqc*.tsv")
    
    script:
    """
    
    echo -e "sample_name\\tscore\\ttrue_positives_baseline\\tfalse_positives\\ttrue_positives_call\\tfalse_negatives\\tprecision\\tsensitivity\\tf_measure" > vcfeval_results_mqc_snps.tsv
    
    find . -name "*snp_roc*" | cut -d "/" -f2 | grep "snp_roc"  |  while read file;
    
    do
    
        # Extract SAMPLE_NAME from the file name using the 'sed' command
        SAMPLE_NAME=\$(echo "\${file}" | sed -E 's/res_vcfeval_(.+)_snp_roc\\.tsv\\.gz/\\1/')

        echo -e "SNP\\t\${SAMPLE_NAME}";

        zcat \${file} | tail -n 1 | sed 's/#//' | sed "s|^|\${SAMPLE_NAME}\\t|" >> vcfeval_results_mqc_snps.tsv
        
    done

    echo -e "sample_name\\tscore\\ttrue_positives_baseline\\tfalse_positives\\ttrue_positives_call\\tfalse_negatives\\tprecision\\tsensitivity\\tf_measure" > vcfeval_results_mqc_indels.tsv
    
    find . -name "*indel_roc*" | cut -d "/" -f2 | grep "indel_roc"  |  while read file;
    
    do
    
        # Extract SAMPLE_NAME from the file name using the 'sed' command
        SAMPLE_NAME=\$(echo "\${file}" | sed -E 's/res_vcfeval_(.+)_indel_roc\\.tsv\\.gz/\\1/')
        
        echo -e "Indels\\t\${SAMPLE_NAME}";
        zcat \${file} | tail -n 1 | sed 's/#//' | sed "s|^|\${SAMPLE_NAME}\\t|" >> vcfeval_results_mqc_indels.tsv
        
    done
    
    """
    
}