nextflow.enable.dsl=2
params.timestamp = ""

process ANNOT_VCF_WITH_VCFEVAL {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/vcfeval_results/annotated_vcf_withvcfeval", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(input_vcf), path(df_class)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val("${sample_name}"), path("${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz*"),emit:all
    tuple val("${sample_name}"), path("${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_tp.vcf.gz*"),emit:tp
    tuple val("${sample_name}"), path("${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fp.vcf.gz*"),emit:fp
    tuple val("${sample_name}"), path("${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fn.vcf.gz*"),emit:fn

    script:
    """
 
    zcat ${df_class}  | grep -v "^#"  | cut -f1,2,4,5,8 |  sed -E 's|SYNC=[0-9]+;||' | sed -E 's|SYNC=[0-9]+,[0-9]+;||' | sed -E 's|CALL_WEIGHT=[0-9]+\\.[0-9]+;|CALL_WEIGHT;|' | sed 's|;|_|g' | bgzip -c > verdict.txt.gz


    tabix -s1 -b2 -e2 verdict.txt.gz


    echo '##INFO=<ID=FlagVCFEval,Number=1,Type=String,Description="VCFEval verdict flag">' > hdr.txt


    bcftools annotate --threads $task.cpus -a verdict.txt.gz -c CHROM,POS,REF,ALT,INFO/FlagVCFEval -h hdr.txt -Oz -o  ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz ${input_vcf[0]}

    bcftools index -t ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz

    
    bcftools filter -i 'INFO/FlagVCFEval="CALL=FP"' -Oz -o ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fp.vcf.gz ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz
    
    bcftools filter -i 'INFO/FlagVCFEval="CALL=TP"' -Oz -o ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_tp.vcf.gz ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz
    
    bcftools filter -i 'INFO/FlagVCFEval="CALL=FN"' -Oz -o ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fn.vcf.gz ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_all.vcf.gz

    bcftools index -t ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fp.vcf.gz

    bcftools index -t ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_fn.vcf.gz

    bcftools index -t ${sample_name}_noref_norm_qualinfo_annotatedwithvcfeval_tp.vcf.gz

    """
    
}