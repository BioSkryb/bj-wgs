nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_VCFEVAL_OVERALL_METRICS {
    tag "CREATE_VCFEVAL_OVERALL_METRICS"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/vcfeval_results", enabled:"$enable_publish"
    
    
    input:
    path(concat_vcfeval_results)
    path(files_negatives)
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("vcfeval_results_mqc.tsv")

    script:
    """
    cat vcfeval_results_mqc.tsv | grep -v "#score" > temp
    mv temp vcfeval_results_mqc.tsv

    Rscript /usr/local/bin/rscript_vcfeval_all_metrics.R
    
    mv vcfeval_results_mqc_manual.tsv vcfeval_results_mqc.tsv
   
    """
    
}
