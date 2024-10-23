nextflow.enable.dsl=2
params.timestamp = ""

process EXTRACT_HET_SITES {
    tag "extract_het_sites"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path( baseline_vcf )
    val( variant_type )
    val( max_af )
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("*.het_1nalt.vcf.gz*"), emit: het_vcf
    path("bcftools_version.yml"), emit: version

    script:
    """
    outfile=`echo ${baseline_vcf[0]} | sed 's/.vcf.gz/.het_1nalt.vcf.gz/'`;
    
    bcftools view -Oz  -i 'TYPE="${variant_type}" && FORMAT/GT="het" && N_ALT==1' ${baseline_vcf[0]} > \${outfile}
    
    bcftools index -t \${outfile}
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools: \$BCFTOOLS_VER > bcftools_version.yml
    """
    
}
