process PREPROCESS_VCF {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/preprocessing", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_vcf)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref_norm_qualinfo.vcf.gz*"), emit: vcf
    
    script:
    """
    # check if the VCF is for DNAscope or TNscope and preprocess accordingly

    if zgrep -q 'ID=DNAscope' ${input_vcf[0]}; then
        echo "Running DNAscope VCF preprocess"
        # filter wildtypes and normalize
        bcftools view --threads ${task.cpus} -i 'GT!="0/0"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    elif zgrep -q '##DeepVariant' ${input_vcf[0]}; then
        echo "Running DeepVariant VCF preprocess"
        # filter wildtypes and normalize
        bcftools view --threads ${task.cpus} -i 'GT!="0/0"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    elif zgrep -q 'ID=TNscope' ${input_vcf[0]}; then
        echo "Running TNscope VCF preprocess"
        normal_sample_name=\$(bcftools query -l ${input_vcf[0]} | tail -1)
        echo "Normal sample name: \${normal_sample_name}"
        bcftools view --threads ${task.cpus} -s ^\${normal_sample_name} -f PASS -e 'SVTYPE="BND" || SVTYPE="INS"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    elif zgrep -q 'ID=TNhaplotyper2' ${input_vcf[0]}; then
        echo "Running TNhaplotyper2 VCF preprocess"
        normal_sample_name=\$(bcftools query -l ${input_vcf[0]} | tail -1)
        echo "Normal sample name: \${normal_sample_name}"
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM ${input_vcf[0]} | bcftools view --threads ${task.cpus} -s ^\${normal_sample_name} -f PASS | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    fi

    bcftools index -t _temp_${sample_name}.vcf.gz

    # add quality information
    bcftools query -f'%CHROM\t%POS\t%QUAL\n' _temp_${sample_name}.vcf.gz | bgzip -c > ${sample_name}_qual.txt.gz
    tabix -s1 -b2 -e2 ${sample_name}_qual.txt.gz

    echo '##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Per-sample QUAL">' > hdr.txt
    bcftools annotate --threads ${task.cpus} -a ${sample_name}_qual.txt.gz -c CHROM,POS,FORMAT/QUAL -h hdr.txt -Oz -o ${sample_name}_noref_norm_qualinfo.vcf.gz _temp_${sample_name}.vcf.gz

    bcftools index -t ${sample_name}_noref_norm_qualinfo.vcf.gz


    """
}
