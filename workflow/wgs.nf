nextflow.enable.dsl=2
import groovy.json.JsonOutput
/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf' addParams(timestamp: params.timestamp)
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_SECONDARY_WF } from '../nf-bioskryb-utils/subworkflows/sentieon_secondary_wf/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_METRICS_WF } from '../nf-bioskryb-utils/modules/sentieon/driver/metrics/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_COVERAGEMETRICS } from '../nf-bioskryb-utils/modules/sentieon/driver/coveragemetrics/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_DATA_PROCESSING_WF } from '../modules/local/custom_data_processing/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_REPORT_WF } from '../modules/local/custom_report/main.nf' addParams(timestamp: params.timestamp)
include { SCORE_VARIANTS_VCFEVAL_WF } from '../nf-bioskryb-utils/subworkflows/score_variants_vcfeval/main.nf' addParams(timestamp: params.timestamp)
include { ADO_WF } from '../nf-bioskryb-utils/subworkflows/ado/main.nf' addParams(timestamp: params.timestamp)
include { BCFTOOLS_VIEW } from '../nf-bioskryb-utils/modules/bcftools/view/main.nf' addParams(timestamp: params.timestamp)
include { SNPEFF_ANNOTATION_WF } from '../nf-bioskryb-utils/modules/snpeff/main.nf' addParams(timestamp: params.timestamp)
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf' addParams(timestamp: params.timestamp)
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf' addParams(timestamp: params.timestamp)

params.reference                    = params.genomes [ params.genome ] [ 'reference' ]
params.base_metrics_intervals       = params.genomes [ params.genome ] [ 'base_metrics_intervals' ]
params.dbsnp                        = params.genomes [ params.genome ] [ 'dbsnp' ]
params.dbsnp_index                  = params.genomes [ params.genome ] [ 'dbsnp_index' ]
params.clinvar                      = params.genomes [ params.genome ] [ 'clinvar' ]
params.clinvar_index                = params.genomes [ params.genome ] [ 'clinvar_index' ]
params.mills                        = params.genomes [ params.genome ] [ 'mills' ]
params.mills_index                  = params.genomes [ params.genome ] [ 'mills_index' ]
params.onekg_omni                   = params.genomes [ params.genome ] [ 'onekg_omni' ]
params.onekg_omni_index             = params.genomes [ params.genome ] [ 'onekg_omni_index' ]
params.dbnsfp                       = params.genomes [ params.genome ] [ 'dbnsfp' ]
params.dbnsfp_index                 = params.genomes [ params.genome ] [ 'dbnsfp_index' ]
params.snpeff_db                    = params.genomes [ params.genome ] [ 'snpeff_db' ]
params.snpeff_genome_name           = params.genomes [ params.genome ] [ 'snpeff_genome_name' ]
params.roi_intervals                = params.genomes [ params.genome ] [ 'roi_intervals' ]
params.refseq                       = params.genomes [ params.genome ] [ 'refseq' ]
params.dnascope_model               = params.genomes [ params.genome ] [ params.platform ] [ params.dnascope_model_selection ] [ params.mode ] [ 'dnascope_model' ]
params.vcfeval_baseline_vcf         = params.genomes [ params.genome ] [ params.giab_reference_name ] [ 'vcfeval_baseline_vcf' ]
params.vcfeval_baseline_vcf_index   = params.genomes [ params.genome ] [ params.giab_reference_name ] [ 'vcfeval_baseline_vcf_index' ]
params.vcfeval_baseline_regions     = params.genomes [ params.genome ] [ params.giab_reference_name ] [ 'vcfeval_baseline_regions' ]
params.giab_genome_name             = params.genomes [ params.genome ] [ params.giab_reference_name ] [ 'giab_genome_name' ]
params.vcfeval_reference_sdf_ref    = params.genomes [ params.genome ] [ 'vcfeval_reference_sdf_ref']

if ( params.mode == 'wgs' ) {
    params.wgs_or_target_intervals      = params.genomes [ params.genome ] [ 'wgs_or_target_intervals' ]
    params.calling_intervals_filename   = params.genomes [ params.genome ] [ 'calling_intervals_filename' ]
    params.vcfeval_bed_regions          = params.genomes [ params.genome ] [ params.giab_reference_name ] [ 'vcfeval_bed_regions' ]
}

if ( params.mode == 'exome' ) {
    params.wgs_or_target_intervals      = params.genomes [ params.genome ] [ params.exome_panel ] [ 'wgs_or_target_intervals' ]
    params.calling_intervals_filename   = params.genomes [ params.genome ] [ params.exome_panel ] [ 'calling_intervals_filename' ]
    params.vcfeval_bed_regions          = params.genomes [ params.genome ] [ params.exome_panel ] [ 'vcfeval_bed_regions' ]
}

workflow WGS_WF {

    take:
        ch_input_csv
        ch_reads
        ch_dummy_file

    main:

        /*
        ========================================================================================
            PROCESS PRIMARY DATA
        ========================================================================================
        */
    
        if ( params.input_csv ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        params.publish_dir,
                                        params.enable_publish
                                    )
        } else {
            ch_input_csv = Channel.fromPath("$projectDir/assets/input_csv_dummy_file.txt", checkIfExists: true).collect()
        }

        if (!params.skip_subsampling) {
            if (!params.subsample_array) {
                log.error "ERROR: subsample_array is empty"
                exit 1
            } else {
                ch_subsample_array = Channel.from(params.subsample_array.toString().split(',').collect { it.trim().toInteger() })
                combined_channel = ch_reads.combine(ch_subsample_array)

                SEQTK_WF (
                    combined_channel,
                    false,
                    params.read_length,
                    params.seqtk_sample_seed,
                    params.publish_dir,
                    params.enable_publish
                )
                ch_reads = SEQTK_WF.out.reads_sampledName
            }
        }
        
        /*
        ========================================================================================
            SENTIEON SECONDARY WORKFLOW
        ========================================================================================
        */
        SENTIEON_SECONDARY_WF ( 
                                params.genome,
                                ch_reads,
                                params.reference,
                                params.dbsnp,
                                params.dbsnp_index,
                                params.mills,
                                params.mills_index,
                                params.onekg_omni,
                                params.onekg_omni_index,
                                params.calling_intervals_filename,
                                params.ploidy,
                                params.dnascope_model,
                                params.pcrfree,
                                params.platform,
                                params.variant_caller,
                                params.publish_dir,
                                params.enable_publish
                             )
        
        ch_bam = SENTIEON_SECONDARY_WF.out.bam
        ch_dnascope_vcf = SENTIEON_SECONDARY_WF.out.dnascope_vcf
        ch_haplotyper_vcf = SENTIEON_SECONDARY_WF.out.haplotyper_vcf
        ch_sentieon_version = SENTIEON_SECONDARY_WF.out.version
        ch_mqc_dedup_metrics = SENTIEON_SECONDARY_WF.out.mqc_dedup_metrics
        // when haplotyper is not executed recal_table will be empty. Need to combine with dummy file to avoid error
        if (SENTIEON_SECONDARY_WF.out.recal_table.count().map { it == 0 }) {
            ch_bam_recal_input = ch_bam.combine(ch_dummy_file)
        } else {
            ch_bam_recal_input = ch_bam.join(SENTIEON_SECONDARY_WF.out.recal_table)
        }

        ch_bam.ifEmpty{ exit 1, "ERROR: cannot generate any alignment files.\nPlease check whether the pipeline parameters/profiles were set correctly or whether the specified files dont have low reads in them." }

        ch_bam
            .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
                { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment/" + it[1].getName() }

        if ( params.genome == "GRCh37" || params.genome == "GRCh38" ){
    
            ch_dnascope_vcf
                .collectFile( name: "vcf_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
                    { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/variant_calls_dnascope/" + it[1][0].getName() }
        }else{
            ch_haplotyper_vcf
                .collectFile( name: "vcf_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
                    { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/variant_calls_haplotyper/" + it[1][0].getName() }
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
        if ( params.genome == "GRCh38" ){
            
            if ( !params.skip_ado ) {
    
                                    ADO_WF(
                                            params.ado_infer_germ_het,
                                            params.ado_slice_het_sites,
                                            params.ado_jg_file,
                                            params.dbsnp,
                                            params.ado_min_samples,
                                            params.ado_variant_type,
                                            params.ado_min_af,
                                            params.ado_max_af,
                                            params.vcfeval_baseline_vcf,
                                            ch_dnascope_vcf,
                                            params.vcfeval_baseline_regions,
                                            params.reference,
                                            params.ado_sample_prop,
                                            params.ado_cov_cutoff,
                                            params.publish_dir,
                                            params.disable_publish
                                          )
                ch_ado_report = ADO_WF.out.ado_summary
                ch_ado_plot = ADO_WF.out.ado_summary_plot
                ch_ado_summary_table = ADO_WF.out.ado_summary_table
                ch_ado_version = ADO_WF.out.ado_version
            } else {
                ch_ado_summary_table = Channel.fromPath("$projectDir/assets/ado_dummy_file.txt", checkIfExists: true).collect()
            }

        }
        
        /*
        ========================================================================================
            METRICS
        ========================================================================================
        */

        SENTIEON_DRIVER_METRICS_WF (
                                      ch_bam_recal_input,
                                      params.reference,
                                      params.base_metrics_intervals,
                                      params.wgs_or_target_intervals,
                                      params.mode,
                                      params.publish_dir,
                                      params.enable_publish
                                  )
        
        ch_sentieon_metrics_version = SENTIEON_DRIVER_METRICS_WF.out.version
        ch_sentieon_metrics = SENTIEON_DRIVER_METRICS_WF.out.metrics
                          
        if ( params.genome == "GRCh38" ){
    
            if ( !params.skip_gene_coverage ) {
        
                SENTIEON_DRIVER_COVERAGEMETRICS (
                                                    ch_bam_recal_input,
                                                    params.reference,
                                                    params.roi_intervals,
                                                    params.refseq,
                                                    params.publish_dir,
                                                    params.enable_publish
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
                                    params.mode,
                                    params.publish_dir,
                                    params.disable_publish
                                  )

        ch_merge_metrics_report = CUSTOM_DATA_PROCESSING_WF.out.metrics
        ch_merge_metrics_version = CUSTOM_DATA_PROCESSING_WF.out.version


        CUSTOM_REPORT_WF ( 
                            ch_merge_metrics_report.collect(),
                            ch_input_csv,
                            ch_ado_summary_table,
                            params.subsample_array,
                            params.publish_dir,
                            params.enable_publish
                         )

        ch_custom_report = CUSTOM_REPORT_WF.out.mqc
        ch_custom_report_version = CUSTOM_REPORT_WF.out.version


    
        /*
        ========================================================================================
            VARIANT FILTER WILDTYPES
        ========================================================================================
        */

        if ( params.genome == "GRCh38" ) {
            if ( !params.skip_variant_annotation || !params.skip_vcfeval ) {

                BCFTOOLS_VIEW (
                                ch_dnascope_vcf,
                                params.publish_dir,
                                params.disable_publish
                )
            }
        } else if ( params.genome == 'GRCm38') {
            if ( !params.skip_variant_annotation ) {

                BCFTOOLS_VIEW (
                                ch_haplotyper_vcf,
                                params.publish_dir,
                                params.disable_publish
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

        if ( params.genome == "GRCh38" ){
            
            if ( !params.skip_vcfeval ) {
    
                SCORE_VARIANTS_VCFEVAL_WF (
                                            BCFTOOLS_VIEW.out.filtered_wt_vcf,
                                            params.vcfeval_baseline_vcf,
                                            params.vcfeval_baseline_vcf_index,
                                            params.vcfeval_baseline_regions,
                                            params.vcfeval_bed_regions,
                                            params.vcfeval_reference_sdf_ref,
                                            params.vcfeval_output_mode,
                                            params.giab_genome_name,
                                            params.vcfeval_score_metric,
                                            params.vcfeval_other_options,
                                            params.publish_dir,
                                            "vcf",
                                            params.reference,
                                            params.enable_publish,
                                            params.disable_publish
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

        if ( params.genome == "GRCh37" || params.genome == "GRCh38" ){
            if ( !params.skip_variant_annotation ) {

                SNPEFF_ANNOTATION_WF (
                                    BCFTOOLS_VIEW.out.filtered_wt_vcf,
                                    params.snpeff_db,
                                    params.snpeff_genome_name,
                                    params.hgvs_old,
                                    params.publish_dir,
                                    params.enable_publish
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
                            params.publish_dir,
                            params.enable_publish
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
            mode: params.mode,
            genome: params.genome,
            variant_caller: params.variant_caller,
            dnascope_model: params.dnascope_model_selection,
            skip_Variant_Annotation: params.skip_variant_annotation,
            skip_Evaluate_Variant_Calling: params.skip_vcfeval,
            skip_ADO_Benchmarking: params.skip_ado
        ]

        MULTIQC_WF ( collect_mqc,
                    params_meta,
                    params.multiqc_config,
                    params.publish_dir,
                    params.enable_publish
                  )

        ch_multiqc_version = MULTIQC_WF.out.version
        MULTIQC_WF.out.report.ifEmpty{ exit 1, "ERROR: cannot generate any MULTIQC Report." }


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
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
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



