params.publish_dir = ""

params.timestamp = ""

process CUSTOM_REPORT {
    tag 'report'
    label 'process_low'
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics", enabled:"$enable_publish"

    
    input:
    path(files)
    path(read_counts_csv)
    val(min_reads)
    path(ado_full_summary)
    val(subsample_array)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*metrics_mqc.txt", emit: mqc
    path("custom_report_version.yml"), emit: version

    
    script:
    """
    python3 /scripts/custom_report.py -s '*sentieonmetrics.txt' -m '*MAPD' -f '*.fastp.json' -i '$read_counts_csv' -r '$min_reads' -a '$ado_full_summary' -p 'wgs' -o 'Summary' -n '$subsample_array'

    echo custom_report: v0.7 > custom_report_version.yml
    
    """
    
}

workflow CUSTOM_REPORT_WF{
    take:
        ch_metrics
        ch_read_counts_csv
        ch_min_reads
        ch_ado_full_summary
        ch_subsample_array
        ch_publish_dir
        ch_enable_publish
    main:
        CUSTOM_REPORT ( 
                        ch_metrics,
                        ch_read_counts_csv,
                        ch_min_reads,
                        ch_ado_full_summary,
                        ch_subsample_array,
                        ch_publish_dir,
                        ch_enable_publish
                      )
    emit:
        mqc = CUSTOM_REPORT.out.mqc
        version = CUSTOM_REPORT.out.version
}

workflow{
    
}
