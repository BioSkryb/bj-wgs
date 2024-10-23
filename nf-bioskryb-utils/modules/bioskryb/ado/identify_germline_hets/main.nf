nextflow.enable.dsl=2
params.timestamp = ""

process IDENTIFY_GERMLINE_HETS {
    tag "identify_germline_hets"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path( jg_file )
    path( dbSNP )
    val( min_samples ) 
    val( variant_type ) 
    val( min_af )
    val( max_af )
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("germline.het.common.vcf.gz*"), emit: vcf
    path("bcftools_version.yml"), emit: version

    script:
    """
    bcftools view -Oz -i '(COUNT(GT="alt")>=${min_samples}) && (INFO/AF>${min_af} & INFO/AF < ${max_af})' -o germline.het.vcf.gz ${jg_file[0]}
    
    bcftools index -t germline.het.vcf.gz
    
    bcftools isec -p isec/ germline.het.vcf.gz ${dbSNP[0]}
    
    bcftools view -i 'type="${variant_type}"' -Oz isec/0002.vcf  > germline.het.common.vcf.gz 
    
    bcftools index -t germline.het.common.vcf.gz
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools: \$BCFTOOLS_VER > bcftools_version.yml
    """
    
}
