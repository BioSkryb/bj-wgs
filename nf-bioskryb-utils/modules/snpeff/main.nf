nextflow.enable.dsl=2
params.timestamp = ""

process SNPEFF_ANNOTATION {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/snpeff", enabled:"$enable_publish"
    
    input:
    tuple val(sample_name), path(vcf)
    path(db)
    val(genome_name)
    val(hgvs_old)
    val(publish_dir)
    val(enable_publish)
    
    output:
    tuple val(sample_name), path("${sample_name}_snpEff.ann.vcf"), emit: vcf
    path("${sample_name}_snpEff.csv"), emit: report
    path("snpEff_version.yml"), emit: version
    
    script:
    def avail_mem = 6
    if(task.memory){
        avail_mem = task.memory.giga
    }
    opts = []
    if(hgvs_old){ opts.add("-hgvsOld") }
    
    def data_dir = "-dataDir \${PWD}/${db}"
    
    """
    snpEff ann\\
        -Xmx${avail_mem}g \\
        -canon -v \\
        -hgvs -hgvs1LetterAa \\
        -nodownload \\
        ${data_dir} \\
        ${genome_name} \\
        -csvStats ${sample_name}_snpEff.csv \\
        ${vcf[0]} \\
        > ${sample_name}_snpEff.ann.vcf
        
    export SNPEFF_VER=\$(snpEff -version 2>&1 | cut -f 2 -d ' ')
    echo SnpEff: \$SNPEFF_VER > snpEff_version.yml
        
    """
}

workflow SNPEFF_ANNOTATION_WF {
    take:
        ch_vcf
        ch_db
        ch_genome_name
        ch_hgvs_old
        ch_publish_dir
        ch_enable_publish
    
    main:
        SNPEFF_ANNOTATION (
                            ch_vcf,
                            ch_db,
                            ch_genome_name,
                            ch_hgvs_old,
                            ch_publish_dir,
                            ch_enable_publish
                          )
    emit:
        vcf = SNPEFF_ANNOTATION.out.vcf
        report = SNPEFF_ANNOTATION.out.report
        version = SNPEFF_ANNOTATION.out.version
}


workflow  {
   
       
       ch_vcf = Channel.fromFilePairs( params.vcf + "/*{.vcf.gz,.vcf.gz.tbi}", size: -1 )
        
        SNPEFF_ANNOTATION_WF (
                            ch_vcf,
                            params.snpEff_db,
                            params.snpEff_genome_name,
                            params.hgvs_old,
                            params.publishDir,
                            params.enable_publish
                          )
   
}