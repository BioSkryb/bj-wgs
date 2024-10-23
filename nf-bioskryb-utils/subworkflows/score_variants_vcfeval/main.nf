nextflow.enable.dsl=2

// IMPORT MODULES

include { PREPROCESS_VCF } from '../../modules/bcftools/preprocess_vcf/main.nf' addParams( timestamp: params.timestamp )
include { BCFTOOLS_QUERY_L } from '../../modules/bcftools/query/main.nf' addParams( timestamp: params.timestamp )
include { VCFEVAL } from '../../modules/vcfeval/main.nf' addParams( timestamp: params.timestamp )
include { ANNOT_VCF_WITH_VCFEVAL } from '../../modules/bioskryb/annot_vcf_with_vcfeval/main.nf' addParams( timestamp: params.timestamp )
include { CONCAT_VCFEVAL_RESULTS } from '../../modules/bioskryb/concat_vcfeval_results/main.nf' addParams( timestamp: params.timestamp )
include { DETERMINE_EXPECTED_NEGATIVES } from '../../modules/bioskryb/truenegatives_vcfeval/determine_expected_negatives/main.nf' addParams( timestamp: params.timestamp )
include { DETERMINE_TOTAL_NEGATIVES } from '../../modules/bioskryb/truenegatives_vcfeval/determine_total_negatives/main.nf' addParams( timestamp: params.timestamp )
include { CREATE_VCFEVAL_OVERALL_METRICS } from '../../modules/bioskryb/create_vcfeval_overall_metrics/main.nf' addParams( timestamp: params.timestamp )


workflow SCORE_VARIANTS_VCFEVAL_WF {
    take:
        ch_vcf
        ch_baseline_vcf
        ch_baseline_vcf_index
        ch_baseline_regions
        ch_bed_regions
        ch_reference_sdf
        ch_output_mode
        ch_reference_name
        ch_score_metric
        ch_other_options
        ch_publish_dir
        ch_vcfeval_mode
        ch_reference
        ch_enable_publish
        ch_disable_publish
        

    main:
    
        PREPROCESS_VCF (
                            ch_vcf,
                            ch_reference,
                            ch_publish_dir,
                            ch_disable_publish
        )


        BCFTOOLS_QUERY_L ( 
                            PREPROCESS_VCF.out.vcf,
                            ch_publish_dir,
                            ch_disable_publish
                         )
                         
        ch_bcftools_version = BCFTOOLS_QUERY_L.out.version
        ch_list = BCFTOOLS_QUERY_L.out.list_names
        
        // ch_list.view()
        
        ch_list_new = ch_list.transpose()
        
        // ch_list_new.view()

        ch_vcfeval_input = ch_vcf.combine(ch_list_new, by: 0)
        
        ch_vcfeval_input.view()
        
        
        // ch_sublist = ch_list
        // .flatMap { item ->
        // samplename = item[0];
        // files  = item[1];
        // return [ [samplename, files ]] }
        

        // ch_vcfeval_input = ch_vcf.combine(ch_sublist, by: 0)
        // ch_vcfeval_input.view()

        VCFEVAL ( 
                    ch_vcfeval_input, 
                    ch_baseline_vcf,
                    ch_baseline_vcf_index,
                    ch_baseline_regions,
                    ch_bed_regions,
                    ch_reference_sdf, 
                    ch_output_mode, 
                    ch_reference_name, 
                    ch_score_metric, 
                    ch_other_options,
                    ch_publish_dir,
                    ch_enable_publish
                )
                
        // VCFEVAL.out.results.collect().view()
        
        ch_input_annot_vcf = PREPROCESS_VCF.out.vcf
        .combine(VCFEVAL.out.df_class,by:0)
        
        ANNOT_VCF_WITH_VCFEVAL (
        
                    ch_input_annot_vcf,
                    ch_publish_dir,
                    ch_enable_publish
        
        )
        
        DETERMINE_EXPECTED_NEGATIVES (
        
                                    ch_baseline_vcf,
                                    ch_baseline_vcf_index,
                                    ch_baseline_regions,
                                    ch_bed_regions,
                                    ch_publish_dir,
                                    ch_disable_publish
                                )
                                
        ch_vcfeval_input_neg = ch_vcfeval_input.
        combine(DETERMINE_EXPECTED_NEGATIVES.out.expected_inert)
        

        DETERMINE_TOTAL_NEGATIVES (
                                    ch_vcfeval_input_neg, 
                                    ch_publish_dir,
                                    ch_disable_publish
                                )                        
                                
                                

        CONCAT_VCFEVAL_RESULTS ( 
                                    VCFEVAL.out.results.collect(),
                                    ch_publish_dir,
                                    ch_disable_publish
                                )
                                
                                
        CREATE_VCFEVAL_OVERALL_METRICS (
                                    CONCAT_VCFEVAL_RESULTS.out,
                                    DETERMINE_TOTAL_NEGATIVES.out.df_eval_neg.collect(),
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
                                
        
        
    emit:
        report = CREATE_VCFEVAL_OVERALL_METRICS.out
        bcftools_version = ch_bcftools_version
        vcfeval_version = VCFEVAL.out.version
}


workflow{
    
    
    ch_vcf = Channel.fromFilePairs( params.vcf )
    ch_vcf.view()

    log.info "vcfeval_baseline_vcf: ${params.vcfeval_baseline_vcf} \n vcfeval_baseline_vcf_index: ${params.vcfeval_baseline_vcf_index} \n vcfeval_baseline_regions: ${params.vcfeval_baseline_regions} \n vcfeval_bed_regions: ${params.vcfeval_bed_regions} \n vcfeval_reference_sdf_ref: ${params.vcfeval_reference_sdf_ref} \n vcfeval_output_mode: ${params.vcfeval_output_mode} \n giab_reference_name: ${params.giab_reference_name} \n vcfeval_score_metric: ${params.vcfeval_score_metric} \n vcfeval_other_options: ${params.vcfeval_other_options} \n publish_dir: ${params.publish_dir} \n enable_publish: ${params.enable_publish} \n disable_publish: ${params.disable_publish}"
    
    // Channel.fromList( params.vcfeval_baseline_vcf.split(',').toList() )
    //       .map { it -> tuple(file("$it"), file("${it}.tbi")) }
    //       .set { ch_vcfeval_baseline_vcf_raw }
    // ch_vcfeval_baseline_vcf_raw
    //       .map{ it -> it.flatten().collect() }
    //       .set{ ch_vcfeval_baseline_vcf }
    
    SCORE_VARIANTS_VCFEVAL_WF ( ch_vcf,
                                params.vcfeval_baseline_vcf,
                                params.vcfeval_baseline_vcf_index,
                                params.vcfeval_baseline_regions,
                                params.vcfeval_bed_regions,
                                params.vcfeval_reference_sdf_ref,
                                params.vcfeval_output_mode,
                                params.giab_reference_name,
                                params.vcfeval_score_metric,
                                params.vcfeval_other_options,
                                params.publish_dir,
                                "vcf",
                                params.reference,
                                params.enable_publish,
                                params.disable_publish
                              )
    
}