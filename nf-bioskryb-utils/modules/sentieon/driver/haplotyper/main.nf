nextflow.enable.dsl=2
params.timestamp = ""


process SENTIEON_HAPLOTYPER {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/variant_calls_haplotyper", enabled:"$enable_publish", saveAs: { it == "${sample_name}.bam*" ? "alignment/$it" : "haplotyper/$it" }

    input:
    tuple val(sample_name), path(bam), path(bai), path(recal_table)
    path fasta_ref
    path interval
    val(ploidy)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("*_haplotyper.vcf.gz*"), emit: vcf
    path("sentieon_haplotyper_version.yml"), emit: version
    
    script:
    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi

    
    sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa \
             -i ${bam} -q ${recal_table} \
             --interval ${interval} \
             --algo Haplotyper \
             --ploidy ${ploidy} \
             --emit_mode gvcf  ${sample_name}_haplotyper.vcf.gz
     

    
    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_haplotyper_version.yml
    """
}

