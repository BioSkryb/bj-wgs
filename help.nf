nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  BJ-WGS   P I P E L I N E
  ===================================
  fastq files    : ${ params.reads }
  genome         : ${ params.genome }
  timestamp      : ${ params.timestamp }
  publish_dir    : ${ params.publish_dir }
  \n
  """

}

def helpMessage() {

  def yellow = "\033[0;33m"
  def blue = "\033[0;34m"
  def white = "\033[0m"
  def red = "\033[0;31m"

  log.info """\

    ${blue}
    BJ-WGS

    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        ${red}
        [required]

        --input_csv                 FILE    Path to input csv file

        --publish_dir               DIR     Path of the output directory

        --genome                    STR     Reference genome to use. Available options - GRCh38, GRCm39
                                    DEFAULT: ${params.genome}

        ${yellow}
        [optional]

        --multiqc_config            DIR     Path to multiqc config
                                            DEFAULT: ${params.multiqc_config}
                                            
        --skip_fastq_merge          BOOL    Skip merging of fastq files
                                            DEFAULT: ${params.skip_fastq_merge}

        --skip_gene_coverage        BOOL    Whether to skip gene coverage
                                            DEFAULT: ${params.skip_gene_coverage}

        --skip_vcfeval              BOOL    Whether to skip vcfeval
                                            DEFAULT: ${params.skip_vcfeval}

        --skip_ado                  BOOL    Whether to skip ADO
                                            DEFAULT: ${params.skip_ado}

        --skip_variant_annotation   BOOL    Whether to skip variant annotation
                                            DEFAULT: ${params.skip_annotation}

        --platform                  STR     Select the sequencing platform. 
                                            Available options - Illumina, Ultima, Element
                                            DEFAULT: Illumina

        --dnascope_model_selection  STR     Select the dnascope model for Illumina runs.
                                            DEFAULT: Sentieon

        --subsample_array           LIST    Specifies a list of read counts for subsampling with seqtk. 
                                            To enable subsampling, ensure that the --skip_subsampling parameter is set to false. 
                                            Provide the read counts as comma-separated values.

        --min_reads                 VAL     Minimum number of reads required for analysis. Samples with fewer reads will be flagged.
                                            DEFAULT: 1000
      
        --help                      BOOL    Display help message

        

    ${yellow}

    """.stripIndent()


}

workflow{
  printHeader()
  helpMessage()
}
