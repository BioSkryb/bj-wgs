nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_ADO_TABLE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val( sample_name ), path( input_vcf )
    path( het_vcf_file )
    path( ref_fasta_file )
    val( sample_prop )
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(sample_name), path("df_ADO_*"), emit: ado_table
    path("bcftools_version.yml"), emit: version

    script:
    """ 
    #Create table
    
    bcftools query -f '%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\n' ${het_vcf_file[0]} | awk 'rand()<${sample_prop}' > ${sample_name}.annos.tsv
    bgzip ${sample_name}.annos.tsv
    tabix -f -c "#" -b 2 -e 3 ${sample_name}.annos.tsv.gz
    
    #Add header to VCF file
    echo '##INFO=<ID=TRUTH_REF,Number=1,Type=String,Description="TRUTH ref allele">' > header_info.txt
    echo '##INFO=<ID=TRUTH_ALT,Number=1,Type=String,Description="TRUTH alt allele"' >> header_info.txt
    
    echo Created annotation file of heterozygous sites
    
    bcftools view --threads $task.cpus -H ${input_vcf[0]} | awk 'BEGIN {OFS="\\t"} \$8~/END/ {ENDVAL=\$8; sub("END=","", ENDVAL); print(\$1,\$2-1,ENDVAL,\$0)} \$8!~/END/ {print(\$1,\$2-1,\$2,\$0)}' | bedtools intersect -wb -a - -b ${het_vcf_file[0]} | awk 'BEGIN {OFS="\\t"} {print(\$4,\$2+1,\$16,\$17,\$8,\$9,\$10,\$11,\$12,\$13)}' > ${sample_name}.het_truth.g.vcf
    bcftools view --threads $task.cpus -h ${input_vcf[0]} | cat - ${sample_name}.het_truth.g.vcf | bcftools view --threads $task.cpus -Oz > ${sample_name}.het_truth.g.vcf.gz
    
    echo Converted target VCF to list sites that are known heterozygous
    
    bcftools convert --threads $task.cpus --gvcf2vcf -f ${ref_fasta_file}/genome.fa ${sample_name}.het_truth.g.vcf.gz | bcftools annotate --threads $task.cpus -h header_info.txt -c 'CHROM,FROM,TO,TRUTH_REF,TRUTH_ALT' -a ${sample_name}.annos.tsv.gz - | bcftools view --threads $task.cpus -Oz -i 'INFO/TRUTH_ALT!="."' > ${sample_name}.het_truth.annotated.g.vcf.gz
    bcftools query -i 'ALT[0]=="<NON_REF>"' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/TRUTH_ALT\\t.\\t[%AD{1}]\\t[%DP]\n' ${sample_name}.het_truth.annotated.g.vcf.gz > ${sample_name}.non_covered_allele_counts.tsv
    
    echo Identified uncovered sites
    
    for j in 0 1 2 3 4;
    do
      bcftools query -i 'INFO/TRUTH_ALT == ALT['\${j}']' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/TRUTH_ALT\\t'\$((\$j + 1))'\\t[%AD{'\$((\$j + 1))'}]\\t[%DP]\\n' ${sample_name}.het_truth.annotated.g.vcf.gz
    done > ${sample_name}.matching_allele_counts.tsv
    
    echo Identified allele frequencies at covered sites
    
    cat ${sample_name}.matching_allele_counts.tsv ${sample_name}.non_covered_allele_counts.tsv | awk 'BEGIN {OFS="\\t"} \$7!="." && \$8!=0 {FREQ=\$7/\$8; print(\$0, FREQ)} \$7=="." || \$8==0 {print(\$0, 0)}' > df_ADO_${sample_name}.tsv
    
    echo  Merged covered and uncovered sites
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools: \$BCFTOOLS_VER > bcftools_version.yml
    """
    
}
