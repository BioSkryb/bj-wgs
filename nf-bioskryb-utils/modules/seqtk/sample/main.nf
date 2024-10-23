nextflow.enable.dsl=2
params.timestamp = ""

process SEQTK_SAMPLE {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/primary_analyses/output/${sample_name}", enabled:"$enable_publish"
    
    input:
    tuple val(sample_name), path(reads), val(n_reads)
    val(is_fasterq)
    val(trim_reads_over_x)
    val(seqtk_sample_seed)
    val(publish_dir)
    val(enable_publish)

    
    output:
    tuple val(sample_name), path("${sample_name}_sampled_R*.fastq.gz"), emit: reads
    tuple val("${sample_name}_${n_reads}"), path("${sample_name}_sampled_R*.fastq.gz"), emit: reads_sampledName
    path("${sample_name}_sample_metadata.json"), emit: metadata
    path("seqtk_version.yml"), emit: version
    
    script:
    if(is_fasterq) {
        r1 = "${reads[0].baseName}.fastq.gz"
        if(reads.size() == 2) {
            r2 = "${reads[1].baseName}.fastq.gz"
        } else {
            r2 = ""
        }
    } else {
        r1 = "${reads[0]}"
        if(reads.size() == 2) {
            r2 = "${reads[1]}"
        } else {
            r2 = ""
        }
    }
    reads_count_paired_cmd = reads.size() == 2 ? "zcat ${r2} | wc -l | awk '{print \$1 / 4}' > read2.txt &" : " "
    reads2_cat = reads.size() == 2 ? "cat ${r2} > ${sample_name}_sampled_R2.fastq.gz &" : " "
    reads_paired_cmd = reads.size() == 2 ? "seqtk sample -s123 ${r2} \$PROPORTION | \
                                            seqtk trimfq -L $trim_reads_over_x - | \
                                            gzip -n > ${sample_name}_sampled_R2.fastq.gz &" : " "
                                            
    """
    zcat ${r1} | wc -l | awk '{print \$1 / 4}' > read1.txt &
    ${reads_count_paired_cmd}
    wait
    
    TOTAL_READS=\$(cat read*.txt|awk '{ sum+=\$0 } END { print sum}')
    PROPORTION=\$(awk -v total_reads=\$TOTAL_READS 'BEGIN { print ${n_reads}/total_reads}')
    PROPORTION=\$(echo \$PROPORTION 1 | awk '{if (\$1 > \$2) print \$2; else print \$1}')
    echo "Total Reads: \$TOTAL_READS, Sample: ${n_reads}, Proportion: \$PROPORTION"
    
    
    if [ \$(awk -v p="\$PROPORTION" 'BEGIN{if (p == 1) print 1; else print 0}') -eq 1 ]; then
        echo "cat the output"
        cat ${r1} > ${sample_name}_sampled_R1.fastq.gz &
        ${reads2_cat}
        wait
    else
        echo "seqtk sample running"
        seqtk sample -s${seqtk_sample_seed} ${r1} \$PROPORTION | seqtk trimfq -L $trim_reads_over_x - | gzip -n > ${sample_name}_sampled_R1.fastq.gz &
        ${reads_paired_cmd}
        wait
    fi
    
    touch ${sample_name}_sample_metadata.json
    echo "{\\"total_reads\\":\$TOTAL_READS," >> ${sample_name}_sample_metadata.json
    echo "\\"subsample_proportions\\":\$PROPORTION}" >> ${sample_name}_sample_metadata.json

    
    #echo "Total Subsample Reads: \$(zcat *_sampled*.fastq.gz | wc -l | awk '{print \$1 / 4}')"
    
    export SEQTK_VER=\$(seqtk 2>&1 | grep -Eo '[0-9].[0-9][-]*[a-z]*[0-9][0-9][0-9][-]*[a-z]*')
    echo Seqtk: \$SEQTK_VER > seqtk_version.yml
    """
}


workflow SEQTK_WF{
    take:
        ch_reads
        ch_is_fasterq
        ch_read_length
        ch_seqtk_sample_seed
        ch_publish_dir
        ch_enable_publish
    main:
        SEQTK_SAMPLE (
                        ch_reads,
                        ch_is_fasterq,
                        ch_read_length,
                        ch_seqtk_sample_seed,
                        ch_publish_dir,
                        ch_enable_publish
                     )
                     
    emit:
        version = SEQTK_SAMPLE.out.version
        reads = SEQTK_SAMPLE.out.reads
        reads_sampledName = SEQTK_SAMPLE.out.reads_sampledName
        metadata = SEQTK_SAMPLE.out.metadata
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , checkExists: true )
                            // .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
    ch_reads.view()
    SEQTK_WF (
                        ch_reads,
                        params.is_fasterq,
                        params.n_reads,
                        params.read_length,
                        params.seqtk_sample_seed,
                        params.publish_dir,
                        params.enable_publish,
                     )
}