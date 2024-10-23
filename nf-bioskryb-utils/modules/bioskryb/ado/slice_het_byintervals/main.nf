nextflow.enable.dsl=2
params.timestamp = ""

process SLICE_HET_INTERVALS {
    tag "slice_het_intervals"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path( het_vcf_file )
    path( bed_regions )
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("*.filtered_by_intervals.vcf.gz*"), emit: het_filtered_vcf
    path("bcftools_version.yml"), emit: version

    script:
    """
    outfile=`echo ${het_vcf_file[0]} | sed 's/.vcf.gz/.filtered_by_intervals.vcf.gz/'`;
    echo "Removing bed headers and sorting"
    
    cat ${bed_regions} | grep -v "^#" | grep -v "^browser" |  grep -v "^track" | sort -V -k1,1 -k2,2 > clean_intervals
    echo "Filtering hets file"
    # see https://samtools.github.io/bcftools/bcftools.html#common_options
    # is regions (-R) or targets (-T) preferable?
    
    bcftools view -R clean_intervals  ${het_vcf_file[0]}  -O z -o \${outfile}
    bcftools  index -t \${outfile}
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools: \$BCFTOOLS_VER > bcftools_version.yml
    """
    
}
