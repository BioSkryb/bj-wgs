nextflow.enable.dsl=2
params.timestamp = ""


include { SENTIEON_DNASCOPE } from '../../modules/sentieon/driver/dnascope/main.nf'
include { SENTIEON_HAPLOTYPER } from '../../modules/sentieon/driver/haplotyper/main.nf'
include { SENTIEON_ALIGNMENT } from '../../modules/sentieon/driver/alignment/main.nf'

workflow SENTIEON_SECONDARY_WF {
    take:
        ch_genome
        ch_primary_reads
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_mills
        ch_mills_index
        ch_onekg_omni
        ch_onekg_omni_index
        ch_calling_intervals_filename
        ch_ploidy
        ch_dnascope_model
        ch_pcrfree
        ch_platform
        ch_variant_caller
        ch_publish_dir
        ch_enable_publish
        
    main:
        // Check for incompatible choices
        if (ch_variant_caller == 'dnascope' && !(ch_genome == "GRCh37" || ch_genome == "GRCh38")) {
            error "ERROR: dnascope is not a valid option for genomes other than GRCh37 or GRCh38."
        }

        SENTIEON_ALIGNMENT ( 
                            ch_genome,
                            ch_primary_reads,
                            ch_reference,
                            ch_dbsnp,
                            ch_dbsnp_index,
                            ch_mills,
                            ch_mills_index,
                            ch_onekg_omni,
                            ch_onekg_omni_index,
                            ch_platform,
                            ch_publish_dir,
                            ch_enable_publish
                         )
        recal_table = SENTIEON_ALIGNMENT.out.bam_recal_table.map { it[3] }

        if (ch_variant_caller == 'haplotyper' || ch_variant_caller == 'all') {
            SENTIEON_HAPLOTYPER ( 
                                    SENTIEON_ALIGNMENT.out.bam_recal_table,
                                    ch_reference,
                                    ch_calling_intervals_filename,
                                    ch_ploidy,
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
            haplotyper_vcf = SENTIEON_HAPLOTYPER.out.vcf
        }else{

            haplotyper_vcf = Channel.empty()
        
        }
                         
        if ( (ch_variant_caller == 'dnascope' || ch_variant_caller == 'all') && (ch_genome == "GRCh37" || ch_genome == "GRCh38") ) {
        
        
            SENTIEON_DNASCOPE ( 
                                 SENTIEON_ALIGNMENT.out.bam,
                                 ch_reference,
                                 ch_calling_intervals_filename,
                                 ch_dbsnp,
                                 ch_dbsnp_index,
                                 ch_dnascope_model,
                                 ch_pcrfree,
                                 ch_ploidy,
                                 "gvcf",
                                 ch_publish_dir,
                                 ch_enable_publish
                               )
            
            dnascope_vcf = SENTIEON_DNASCOPE.out.vcf
        }else{
        
        
            dnascope_vcf = Channel.empty()
            
        }
                         
    emit:
        bam = SENTIEON_ALIGNMENT.out.bam
        dedup_metrics = SENTIEON_ALIGNMENT.out.dedup_metrics
        mqc_dedup_metrics = SENTIEON_ALIGNMENT.out.mqc_dedup_metrics
        recal_table = recal_table
        haplotyper_vcf = haplotyper_vcf
        dnascope_vcf = dnascope_vcf
        version = SENTIEON_ALIGNMENT.out.version
    
}
