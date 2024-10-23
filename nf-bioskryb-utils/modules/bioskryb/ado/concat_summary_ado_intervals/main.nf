nextflow.enable.dsl=2
params.timestamp = ""

process CONCAT_SUMMARY_ADO_INTERVALS {
    tag "summary_ado_files"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path(df_files)
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("merged_ADO.tsv"), emit: merged_ADO
    path("ADO_plot_summary.png"), emit: plot_ADO
    path("merged_ADO_summary.tsv"), emit: summary_ADO
    path("ado_version.yml"), emit: version

    script:
    def containerName = task.container
    """
    echo -e "File_Interval\\tFreq\\tProp" > merged_ADO.tsv
    
    cat res_ADO_* | cut -f1  | sed 's/\\.tsv//' > a.txt
    
    cat res_ADO_*  | sed 's/\\[//' | sed 's/)//' | sed 's/]//'  | cut -f2 > b.txt
    
    cat res_ADO_*  | cut -f3 > c.txt
    
    cat res_ADO_*  | cut -f4 > d.txt
    
    paste -d "_" a.txt b.txt > end.txt
    
    paste -d "\\t" end.txt c.txt d.txt >> merged_ADO.tsv

    Rscript /usr/local/bin/plot_summary_ado_intervals.R merged_ADO.tsv

    export CONTAINER_NAME=${containerName}
    export ADO_VER=\$(echo \$CONTAINER_NAME | awk -F'/:' '{print \$NF}' | awk -F'_' '{print \$NF}')
    echo ADO: \$ADO_VER > ado_version.yml
    """
    
}
