/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf'
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf'
include { SENTIEON_SECONDARY_WF } from '../nf-bioskryb-utils/subworkflows/sentieon_secondary_wf/main.nf'
include { SENTIEON_DRIVER_METRICS_WF } from '../nf-bioskryb-utils/modules/sentieon/driver/metrics/main.nf'
include { SENTIEON_DRIVER_COVERAGEMETRICS } from '../nf-bioskryb-utils/modules/sentieon/driver/coveragemetrics/main.nf'
include { CUSTOM_DATA_PROCESSING_WF } from '../modules/local/custom_data_processing/main.nf'
include { CUSTOM_REPORT_WF } from '../modules/local/custom_report/main.nf'
include { SCORE_VARIANTS_VCFEVAL_WF } from '../nf-bioskryb-utils/subworkflows/score_variants_vcfeval/main.nf'
include { PREPROCESS_VCF } from '../nf-bioskryb-utils/modules/bcftools/preprocess_vcf/main.nf'
include { ADO_WF } from '../nf-bioskryb-utils/subworkflows/ado/main.nf'
include { BCFTOOLS_VIEW } from '../nf-bioskryb-utils/modules/bcftools/view/main.nf'
include { SNPEFF_ANNOTATION_WF } from '../nf-bioskryb-utils/modules/snpeff/main.nf'
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf'
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf'
include { COUNT_READS_FASTQ_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf'
include { SIGPROFILERGENERATEMATRIX_WF} from '../nf-bioskryb-utils/modules/sigprofilermatrixgenerator/main.nf'


workflow WGS_WF {

    take:
        ch_input_csv
        ch_reads
        ch_dummy_file
        ch_min_reads
        ch_genomes
        ch_genome
        ch_platform
        ch_dnascope_model_selection
        ch_mode
        ch_giab_reference_name
        ch_exome_panel
        ch_variant_caller
        ch_skip_subsampling
        ch_subsample_array
        ch_read_length
        ch_seqtk_sample_seed
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_mills
        ch_mills_index
        ch_onekg_omni
        ch_onekg_omni_index
        ch_ploidy
        ch_pcrfree
        ch_tmp_dir
        ch_timestamp
        ch_skip_ado
        ch_ado_infer_germ_het
        ch_ado_slice_het_sites
        ch_ado_jg_file
        ch_ado_min_samples
        ch_ado_variant_type
        ch_ado_min_af
        ch_ado_max_af
        ch_ado_sample_prop
        ch_ado_cov_cutoff
        ch_base_metrics_intervals
        ch_skip_gene_coverage
        ch_roi_intervals
        ch_refseq
        ch_skip_variant_annotation
        ch_skip_vcfeval
        ch_vcfeval_output_mode
        ch_vcfeval_score_metric
        ch_vcfeval_other_options
        ch_snpeff_db
        ch_snpeff_genome_name
        ch_hgvs_old
        ch_multiqc_config
        ch_sigprofilermatrixgenerator_reference
        ch_report_s3_dir
        ch_skip_sigprofile
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish


    main:

        dnascope_model                      = ""
        if ( ch_genome == 'GRCh38' ) {
            dnascope_model               = ch_genomes [ ch_genome ] [ ch_platform ] [ ch_dnascope_model_selection ] [ ch_mode ] [ 'dnascope_model' ]
            vcfeval_baseline_vcf         = ch_genomes [ ch_genome ] [ ch_giab_reference_name ] [ 'vcfeval_baseline_vcf' ]
            vcfeval_baseline_vcf_index   = ch_genomes [ ch_genome ] [ ch_giab_reference_name ] [ 'vcfeval_baseline_vcf_index' ]
            vcfeval_baseline_regions     = ch_genomes [ ch_genome ] [ ch_giab_reference_name ] [ 'vcfeval_baseline_regions' ]
            giab_genome_name             = ch_genomes [ ch_genome ] [ ch_giab_reference_name ] [ 'giab_genome_name' ]
            vcfeval_reference_sdf_ref    = ch_genomes [ ch_genome ] [ 'vcfeval_reference_sdf_ref']
        }

        if ( ch_mode == 'wgs' ) {
            wgs_or_target_intervals      = ch_genomes [ ch_genome ] [ 'wgs_or_target_intervals' ]
            calling_intervals_filename   = ch_genomes [ ch_genome ] [ 'calling_intervals_filename' ]
            sigprofilermatrixgenerator_interval = ch_genomes [ ch_genome ] [ 'sigprofilermatrixgenerator_interval' ]
            if ( ch_genome == 'GRCh38' ) {
                vcfeval_bed_regions          = ch_genomes [ ch_genome ] [ ch_giab_reference_name ] [ 'vcfeval_bed_regions' ]
            }
        }

        if ( ch_mode == 'exome' ) {
            wgs_or_target_intervals      = ch_genomes [ ch_genome ] [ ch_exome_panel ] [ 'wgs_or_target_intervals' ]
            calling_intervals_filename   = ch_genomes [ ch_genome ] [ ch_exome_panel ] [ 'calling_intervals_filename' ]
            sigprofilermatrixgenerator_interval = ch_genomes [ ch_genome ] [ ch_exome_panel ] [ 'vcfeval_bed_regions' ]
            if ( ch_genome == 'GRCh38' ) {
                vcfeval_bed_regions          = ch_genomes [ ch_genome ] [ ch_exome_panel ] [ 'vcfeval_bed_regions' ]
            }
        }

        // Override variant caller for GRCm39 genome (make sure to use the below variable in the downstream instead of ch_variant_caller)
        variant_caller_modified = (ch_genome == 'GRCm39') ? 'haplotyper' : ch_variant_caller
        println "INFO: Variant caller is set to ${variant_caller_modified} for ${ch_genome} genome"

        /*
        ========================================================================================
            PROCESS PRIMARY DATA
        ========================================================================================
        */
    
        if ( ch_input_csv ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )
        } else {
            ch_input_csv = Channel.fromPath("$projectDir/assets/input_csv_dummy_file.txt", checkIfExists: true).collect()
        }

        if (!ch_skip_subsampling) {
            if (!ch_subsample_array) {
                error "ERROR: subsample_array is empty"
            } else {
                ch_subsample_array = Channel.from(ch_subsample_array.toString().split(',').collect { item -> item.trim().toInteger() })
                combined_channel = ch_reads.combine(ch_subsample_array)

                SEQTK_WF (
                    combined_channel,
                    false,
                    ch_read_length,
                    ch_seqtk_sample_seed,
                    ch_publish_dir,
                    ch_enable_publish
                )
                ch_reads = SEQTK_WF.out.reads_sampledName
            }
        }

        COUNT_READS_FASTQ_WF (
                                ch_reads,
                                ch_publish_dir,
                                ch_enable_publish
                            )
        COUNT_READS_FASTQ_WF.out.read_counts
        .map { sample_id, files, read_count_file -> 
            def read_count = read_count_file.text.trim().toLong()
            [sample_id, files, read_count]
        }
        .branch { read ->
            small: read[2] < ch_min_reads
            large: read[2] >= ch_min_reads
        }
        .set { branched_reads }

        ch_fastqs = branched_reads.large.map { sample_id, files, _read_count ->
            tuple(sample_id, files)
        }

        /*
        ========================================================================================
            SENTIEON SECONDARY WORKFLOW
        ========================================================================================
        */
        SENTIEON_SECONDARY_WF ( 
                                ch_genome,
                                ch_fastqs,
                                ch_reference,
                                ch_dbsnp,
                                ch_dbsnp_index,
                                ch_mills,
                                ch_mills_index,
                                ch_onekg_omni,
                                ch_onekg_omni_index,
                                calling_intervals_filename,
                                ch_ploidy,
                                dnascope_model,
                                ch_pcrfree,
                                ch_platform,
                                variant_caller_modified,
                                ch_publish_dir,
                                ch_enable_publish
                             )
        
        ch_bam = SENTIEON_SECONDARY_WF.out.bam
        ch_dnascope_vcf = SENTIEON_SECONDARY_WF.out.dnascope_vcf
        ch_haplotyper_vcf = SENTIEON_SECONDARY_WF.out.haplotyper_vcf
        ch_sentieon_version = SENTIEON_SECONDARY_WF.out.version
        ch_mqc_dedup_metrics = SENTIEON_SECONDARY_WF.out.mqc_dedup_metrics
        // when haplotyper is not executed recal_table will be empty. Need to combine with dummy file to avoid error
        if (SENTIEON_SECONDARY_WF.out.recal_table.count().map { count -> count == 0 }) {
            ch_bam_recal_input = ch_bam.combine(ch_dummy_file)
        } else {
            ch_bam_recal_input = ch_bam.join(SENTIEON_SECONDARY_WF.out.recal_table)
        }

        ch_bam.ifEmpty{ error "ERROR: cannot generate any alignment files.\nPlease check whether the pipeline parameters/profiles were set correctly or whether the specified files dont have low reads in them." }

        ch_bam
            .collectFile( name: "bam_files.txt", newLine: true, sort: { sample -> sample[0] }, storeDir: "${ch_tmp_dir}" )
                { sample -> sample[0] + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/alignment/" + sample[1].getName() }

        if ( (ch_variant_caller == 'dnascope' ) && (ch_genome == "GRCh37" || ch_genome == "GRCh38") ) {
    
            ch_dnascope_vcf
                .collectFile( name: "vcf_files.txt", newLine: true, sort: { item -> item[0] }, storeDir: "${ch_tmp_dir}" )
                    { item -> item[0] + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/variant_calls_dnascope/" + item[1][0].getName() }
        }            
        if (ch_variant_caller == 'haplotyper') {
            ch_haplotyper_vcf
                .collectFile( name: "vcf_files.txt", newLine: true, sort: { item -> item[0] }, storeDir: "${ch_tmp_dir}" )
                    { item -> item[0] + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/variant_calls_haplotyper/" + item[1][0].getName() }
        }

        /*
        ========================================================================================
            MUTATION SIGNATURE
        ========================================================================================
        */
        if ( !ch_skip_sigprofile ) {
            // Selecting the non empty vcf channel from dnascope or haplotyper
            ch_vcf_out = ch_dnascope_vcf.mix(ch_haplotyper_vcf).ifEmpty { null }
            SIGPROFILERGENERATEMATRIX_WF (
                ch_vcf_out.map { sample_name, vcf_files -> [sample_name, vcf_files[0], vcf_files[1]] },
                sigprofilermatrixgenerator_interval,
                ch_reference,
                ch_genome,
                ch_sigprofilermatrixgenerator_reference,
                ch_report_s3_dir,
                ch_publish_dir,
                ch_enable_publish,
                ch_disable_publish
            )
        }

        /*
        ========================================================================================
            Alelic balance (ADO) EVALUATION
        ========================================================================================
        */

        ch_ado_report = Channel.empty()
        ch_ado_plot = Channel.empty()
        ch_ado_summary_table = Channel.empty()
        ch_ado_version = Channel.empty()
        if ( ch_genome == "GRCh38" ){
            
            if ( !ch_skip_ado ) {
    
                                    ADO_WF(
                                            ch_ado_infer_germ_het,
                                            ch_ado_slice_het_sites,
                                            ch_ado_jg_file,
                                            ch_dbsnp,
                                            ch_ado_min_samples,
                                            ch_ado_variant_type,
                                            ch_ado_min_af,
                                            ch_ado_max_af,
                                            vcfeval_baseline_vcf,
                                            ch_dnascope_vcf,
                                            vcfeval_baseline_regions,
                                            ch_reference,
                                            ch_ado_sample_prop,
                                            ch_ado_cov_cutoff,
                                            ch_publish_dir,
                                            ch_disable_publish
                                          )
                ch_ado_report = ADO_WF.out.ado_summary
                ch_ado_plot = ADO_WF.out.ado_summary_plot
                ch_ado_summary_table = ADO_WF.out.ado_summary_table
                ch_ado_version = ADO_WF.out.ado_version
            } else {
                ch_ado_summary_table = Channel.fromPath("$projectDir/assets/ado_dummy_file.txt", checkIfExists: true).collect()
            }

        } else {
            ch_ado_summary_table = Channel.fromPath("$projectDir/assets/ado_dummy_file.txt", checkIfExists: true).collect()
        }
        
        /*
        ========================================================================================
            METRICS
        ========================================================================================
        */

        SENTIEON_DRIVER_METRICS_WF (
                                      ch_bam_recal_input,
                                      ch_reference,
                                      ch_base_metrics_intervals,
                                      wgs_or_target_intervals,
                                      ch_mode,
                                      "dedup", // dedup or non_dedup
                                      ch_publish_dir,
                                      ch_enable_publish
                                  )
        
        ch_sentieon_metrics_version = SENTIEON_DRIVER_METRICS_WF.out.version
        ch_sentieon_metrics = SENTIEON_DRIVER_METRICS_WF.out.metrics
                          
        if ( ch_genome == "GRCh38" ){
    
            if ( !ch_skip_gene_coverage ) {
        
                SENTIEON_DRIVER_COVERAGEMETRICS (
                                                    ch_bam_recal_input,
                                                    ch_reference,
                                                    ch_roi_intervals,
                                                    ch_refseq,
                                                    ch_publish_dir,
                                                    ch_enable_publish
                                            )
            }
        }

        ch_custom_data_processing_input = SENTIEON_SECONDARY_WF.out.dedup_metrics
                                                .combine (
                                                            ch_bam
                                                                .join( SENTIEON_DRIVER_METRICS_WF.out.metrics_tuple ),
                                                            by: 0
                                                         )


        CUSTOM_DATA_PROCESSING_WF ( 
                                    ch_custom_data_processing_input,
                                    ch_mode,
                                    ch_publish_dir,
                                    ch_disable_publish
                                  )

        ch_merge_metrics_report = CUSTOM_DATA_PROCESSING_WF.out.metrics


        CUSTOM_REPORT_WF ( 
                            ch_merge_metrics_report.collect(),
                            COUNT_READS_FASTQ_WF.out.combined_read_counts,
                            ch_min_reads,
                            ch_ado_summary_table,
                            ch_subsample_array,
                            ch_publish_dir,
                            ch_enable_publish
                         )

        ch_custom_report = CUSTOM_REPORT_WF.out.mqc


    
        /*
        ========================================================================================
            VARIANT FILTER WILDTYPES
        ========================================================================================
        */

        if ( ch_genome == "GRCh38" ) {
            if ( !ch_skip_variant_annotation || !ch_skip_vcfeval ) {

                BCFTOOLS_VIEW (
                                ch_dnascope_vcf,
                                ch_publish_dir,
                                ch_disable_publish
                )
            }
        } else if ( ch_genome == 'GRCm38') {
            if ( !ch_skip_variant_annotation ) {

                BCFTOOLS_VIEW (
                                ch_haplotyper_vcf,
                                ch_publish_dir,
                                ch_disable_publish
                )
            }
        }
    
        /*
        ========================================================================================
            VARIANT EVALUATION
        ========================================================================================
        */
        ch_vcfeval_report = Channel.empty()
        ch_bcftools_version = Channel.empty()
        ch_vcfeval_version = Channel.empty()

        if ( ch_genome == "GRCh38" ){
            
            if ( !ch_skip_vcfeval ) {

                PREPROCESS_VCF (
                                            BCFTOOLS_VIEW.out.filtered_wt_vcf,
                                            ch_reference,
                                            ch_publish_dir,
                                            ch_enable_publish
                                          )
    
                SCORE_VARIANTS_VCFEVAL_WF (
                                            PREPROCESS_VCF.out.vcf,
                                            vcfeval_baseline_vcf,
                                            vcfeval_baseline_vcf_index,
                                            vcfeval_baseline_regions,
                                            vcfeval_bed_regions,
                                            vcfeval_reference_sdf_ref,
                                            ch_vcfeval_output_mode,
                                            giab_genome_name,
                                            ch_vcfeval_score_metric,
                                            ch_vcfeval_other_options,
                                            ch_publish_dir,
                                            "vcf",
                                            ch_reference,
                                            ch_enable_publish,
                                            ch_disable_publish
                                          )
                ch_vcfeval_report = SCORE_VARIANTS_VCFEVAL_WF.out.report
                ch_bcftools_version = SCORE_VARIANTS_VCFEVAL_WF.out.bcftools_version
                ch_vcfeval_version = SCORE_VARIANTS_VCFEVAL_WF.out.vcfeval_version
                
            }
        }

      
        /*
        ========================================================================================
            VARIANT ANNOTATION
        ========================================================================================
        */
        ch_snpeff_version = Channel.empty()
        ch_snpeff_report = Channel.empty()        

        if ( ch_genome == "GRCh37" || ch_genome == "GRCh38" ){
            if ( !ch_skip_variant_annotation ) {

                SNPEFF_ANNOTATION_WF (
                                    BCFTOOLS_VIEW.out.filtered_wt_vcf,
                                    ch_snpeff_db,
                                    ch_snpeff_genome_name,
                                    ch_hgvs_old,
                                    ch_publish_dir,
                                    ch_enable_publish
                                  )
                                  
                ch_snpeff_report = SNPEFF_ANNOTATION_WF.out.report
                ch_snpeff_version = SNPEFF_ANNOTATION_WF.out.version
    
    
            }
        }

        /*
        ========================================================================================
            REPORTING
        ========================================================================================
        */

        ch_tool_versions = ch_sentieon_version.take(1).ifEmpty([])
                                    .combine(ch_sentieon_metrics_version.take(1).ifEmpty([]))
                                    .combine(ch_bcftools_version.take(1).ifEmpty([]))
                                    .combine(ch_vcfeval_version.take(1).ifEmpty([]))
                                    .combine(ch_snpeff_version.take(1).ifEmpty([]))
                                    .combine((ch_ado_version ?: Channel.empty()).ifEmpty([]))
                                    

        REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            ch_publish_dir,
                            ch_enable_publish
                          )


        collect_mqc = ch_sentieon_metrics.collect().ifEmpty([])
                        .combine( ch_custom_report.collect().ifEmpty([]))
                        .combine( ch_vcfeval_report.collect().ifEmpty([]))
                        .combine( ch_snpeff_report.collect().ifEmpty([]))
                        .combine( ch_ado_report.collect().ifEmpty([]))
                        .combine( ch_ado_summary_table.collect().ifEmpty([]))
                        .combine( ch_ado_plot.collect().ifEmpty([]))
                        .combine( ch_mqc_dedup_metrics.collect().ifEmpty([]))
                        .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))

        params_meta = [
            session_id: workflow.sessionId,
            mode: ch_mode,
            genome: ch_genome,
            variant_caller: variant_caller_modified,
            skip_Variant_Annotation: ch_skip_variant_annotation,
            skip_Evaluate_Variant_Calling: ch_skip_vcfeval,
            skip_ADO_Benchmarking: ch_skip_ado
        ]

        if (variant_caller_modified == 'dnascope') {
            params_meta['dnascope_model'] = ch_dnascope_model_selection
        }

        MULTIQC_WF ( collect_mqc,
                    params_meta,
                    ch_multiqc_config,
                    ch_publish_dir,
                    ch_enable_publish
                  )

        MULTIQC_WF.out.report.ifEmpty{ error "ERROR: cannot generate any MULTIQC Report." }


}



workflow.onComplete {
    output                              = [:]
    output["pipeline_run_name"]         = workflow.runName
    output["pipeline_name"]             = workflow.manifest.name
    output["pipeline_version"]          = workflow.manifest.version
    output["pipeline_session_id"]       = workflow.sessionId
    output["output"]                    = [:]
    output["output"]["bam"]             = [:]
    output["output"]["vcf"]             = [:]


    bam_outfile = file("$params.tmp_dir/bam_files.txt")
    bam_outfile_lines = bam_outfile.readLines()
    for ( bam_line : bam_outfile_lines ) {
        def (sample_name, bam_path) = bam_line.split('\t')
        output["output"]["bam"][sample_name] = [:]
        output["output"]["bam"][sample_name]["bam"] = bam_path
    }
    
    vcf_outfile = file("$params.tmp_dir/vcf_files.txt")
    vcf_outfile_lines = vcf_outfile.readLines()
    for ( vcf_line : vcf_outfile_lines ) {
        def (sample_name, vcf_path) = vcf_line.split('\t')
        output["output"]["vcf"][sample_name] = [:]
        output["output"]["vcf"][sample_name]["vcf"] = vcf_path
    }
    
    def output_json = groovy.json.JsonOutput.toJson(output)
    def output_json_pretty = groovy.json.JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
}


workflow.onError {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]

    
    

    def subject = """\
        [nf-wgs-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\

        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}



