nextflow.enable.dsl=2

// IMPORT MODULES

include { IDENTIFY_GERMLINE_HETS } from '../../modules/bioskryb/ado/identify_germline_hets/main.nf'
include { EXTRACT_HET_SITES } from '../../modules/bioskryb/ado/extract_het_sites/main.nf'
include { SLICE_HET_INTERVALS } from '../../modules/bioskryb/ado/slice_het_byintervals/main.nf'
include { CREATE_ADO_TABLE } from '../../modules/bioskryb/ado/create_ado_table/main.nf'
include { SUMMARIZE_ADO_INTERVALS } from '../../modules/bioskryb/ado/summarize_ado_intervals_r/main.nf'
include { CONCAT_SUMMARY_ADO_INTERVALS } from '../../modules/bioskryb/ado/concat_summary_ado_intervals/main.nf'

workflow ADO_WF{
    take:
        ch_infer_germ_het
        ch_slice_het_sites
        ch_jg
        ch_dbnsp
        ch_min_samples
        ch_variant_type
        ch_min_af
        ch_max_af
        ch_baseline_vcf
        ch_vcf
        ch_bed_regions
        ch_ref_fasta
        ch_sample_prop
        ch_cov_cutoff
        ch_publish_dir
        ch_enable_publish
    main:
        
        // ch_baseline_vcf.view()
        if(ch_infer_germ_het){
            IDENTIFY_GERMLINE_HETS (
                                        ch_jg,
                                        ch_dbnsp,
                                        ch_min_samples,
                                        ch_variant_type,
                                        ch_min_af,
                                        ch_max_af,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )
            ch_baseline_vcf = IDENTIFY_GERMLINE_HETS.out.vcf
            // ch_baseline_vcf.view()
        }
        
        
        EXTRACT_HET_SITES ( 
                            ch_baseline_vcf,
                            ch_variant_type,
                            ch_max_af,
                            ch_publish_dir,
                            ch_enable_publish
                          )
        ch_sites = EXTRACT_HET_SITES.out.het_vcf
        if(ch_slice_het_sites){
            SLICE_HET_INTERVALS ( 
                                    ch_sites,
                                    ch_bed_regions,
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
            ch_sites = SLICE_HET_INTERVALS.out.het_filtered_vcf
        }
        
        
        CREATE_ADO_TABLE ( 
                            ch_vcf,
                            ch_sites.collect(),
                            ch_ref_fasta,
                            ch_sample_prop,
                            ch_publish_dir,
                            ch_enable_publish
                         )
        
        SUMMARIZE_ADO_INTERVALS (
                                    CREATE_ADO_TABLE.out.ado_table,
                                    ch_cov_cutoff,
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
        all_df_sum = SUMMARIZE_ADO_INTERVALS.out.df_sum.collect()
        CONCAT_SUMMARY_ADO_INTERVALS (
                                        all_df_sum,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )
    emit:
        ado_summary = CONCAT_SUMMARY_ADO_INTERVALS.out.merged_ADO
        ado_summary_plot = CONCAT_SUMMARY_ADO_INTERVALS.out.plot_ADO
        ado_summary_table = CONCAT_SUMMARY_ADO_INTERVALS.out.summary_ADO
        ado_version = CONCAT_SUMMARY_ADO_INTERVALS.out.version

}


workflow  {
    
    
    Channel.fromList( params.jg_file.split(',').toList() )
           .map { it -> tuple(file("$it"), file("${it}.tbi")) }
           .set { ch_jg_raw }
    ch_jg_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_jg }
              
    
    Channel.fromList( params.dbsnp.split(',').toList() )
          .map { it -> tuple(file("$it"), file("${it}.tbi")) }
          .set { ch_dbsnp_raw }
    ch_dbsnp_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_dbsnp }

    Channel.fromList( params.baseline_vcf.split(',').toList() )
          .map { it -> tuple(file("$it"), file("${it}.tbi")) }
          .set { ch_baseline_vcf_raw }
    ch_baseline_vcf_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_baseline_vcf }
    
    // ch_bed_regions = Channel.fromPath(params.baseline_regions)
    
    // ch_ref_fasta =  Channel.fromPath(params.reference)
    
    
    if(params.vcf != ""){
        ch_vcf = Channel.fromFilePairs( params.vcf, size: -1)
    }
    else if(params.input_csv != ""){
        ch_vcf = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.sampleId, [ row.vcf, row.vcf + ".tbi" ] ] }
    }
       
    
    ADO_WF(
            params.infer_germ_het,
            params.slice_het_sites,
            ch_jg,
            ch_dbsnp,
            params.min_samples,
            params.variant_type,
            params.min_af,
            params.max_af,
            ch_baseline_vcf,
            ch_vcf,
            params.baseline_regions,
            params.reference,
            params.sample_prop,
            params.cov_cutoff,
            params.publish_dir,
            params.enable_publish
          )
        
        
}