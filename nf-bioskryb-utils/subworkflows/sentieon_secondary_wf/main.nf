nextflow.enable.dsl=2
params.timestamp = ""


include { SENTIEON_DNASCOPE } from '../../modules/sentieon/driver/dnascope/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_HAPLOTYPER } from '../../modules/sentieon/driver/haplotyper/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_ALIGNMENT } from '../../modules/sentieon/driver/alignment/main.nf' addParams( timestamp: params.timestamp )

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

workflow {
    
    params.reference                    = WorkflowWGS.getGenomeAttribute( params, 'reference' )
    params.intervals                    = WorkflowWGS.getGenomeAttribute( params, 'base_metrics_intervals' )
    params.dbsnp                        = WorkflowWGS.getGenomeAttribute( params, 'dbsnp' )
    params.dbsnp_index                  = WorkflowWGS.getGenomeAttribute( params, 'dbsnp_index' )
    params.mills                        = WorkflowWGS.getGenomeAttribute( params, 'mills' )
    params.mills_index                  = WorkflowWGS.getGenomeAttribute( params, 'mills_index' )
    params.onekg_omni                   = WorkflowWGS.getGenomeAttribute( params, 'onekg_omni' )
    params.onekg_omni_index             = WorkflowWGS.getGenomeAttribute( params, 'onekg_omni_index' )
    params.dnascope_model               = WorkflowWGS.getGenomeAttribute( params, 'dnascope_model' )

    if ( params.mode == 'wgs' ) {
        params.calling_intervals_filename   = WorkflowWGS.getGenomeAttribute( params, 'calling_intervals_filename' )
    }
    
    if ( params.mode == 'exome' ) {
        params.calling_intervals_filename   = WorkflowWGS.getExomeAttribute( params, 'calling_intervals_filename' )

        log.info """\
        interval       : ${ params.wgs_or_target_intervals }
        \n
        """
    }

    
    ch_primary_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
    
    
    
    
    SENTIEON_SECONDARY_WF ( 
                                params.genome,
                                ch_primary_reads,
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
}