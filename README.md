
# BaseJumper (BJ-WGS)

BJ-WGS pipeline is a scalable and reproducible bioinformatics pipeline to process single-cell sequencing data from ResolveDNA Whole Genome Amplification or any single-cell or bulk sequencing data. The pipeline currently supports human and mouse sequencing data but can certainly be extended to other model systems. It supports sequencing data from Illumina, Ultima, and Element. The pipeline takes raw sequencing data in form of fastq/cram files and performs alignment, removes duplicate reads, base calibrates the reads, and performs variant calling with haplotype caller and DNAScope caller. For Illumina sequencing data, users have option to use primary template amplication (PTA) corrected DNAScope model to do variant calling.

# Pipeline Overview
Following are the steps and tools that pipeline uses to perform the analyses:

- Map reads to reference genome using SENTIEON BWA MEM
- Remove duplicate reads using SENTIEON DRIVER LOCUSCOLLECTOR and SENTIEON DRIVER DEDUP
- Perform base quality score recalibration (BQSR) using SENTIEON DRIVER BQSR
- Perform variant calling with DNAScope (defualt) or HAPLOTYPER caller
- Perform variant annotation with SNPEFF, ClinVar, and dbSNP databases
- Evaluate metrics using SENTIEON DRIVER METRICS which includes Alignment, GC Bias, Insert Size, and Coverage metrics
- Evaluate variants with VCFEval to assess analytical performance (Only supported for GIAB reference samples: HG001-HG007 samples)
- Aggregate the metrics across biosamples and tools to create overall pipeline statistics summary using MULTIQC

# Running Locally

Following are instructions for running BJ-WGS in a local Ubuntu server

## Install Java

```
sudo apt-get install default-jdk

java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Sentieon License Setup

The Sentieon license is a "localhost" license that starts a lightweight license server on the localhost. This type of license is very easy to use and get started with. However, because it can be used anywhere, we restrict this license to short-term testing/evaluation only. To use this type of license, you need to set the environment variable SENTIEON_LICENSE to point to the license file on the compute nodes:
```
export SENTIEON_LICENSE=</path/to/sentieon_eval.lic>
```
The license file should be saved at the base directory of the pipeline eg: `bj-wgs/sentieon_eval.lic`
All users will need to  [ submit helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156) to get an evaluation/full pass-through BioSkryb's Sentieon license.

## Resources Required

For running the pipeline, a typical dataset requires 8 CPU cores and 50 GB of memory. For larger datasets, you may need to increase the resources to 64 CPU cores and 120 GB of memory. You can specify these resources in the command as follows:
```
--max_cpus 8 --max_memory 50.GB
```
## Test Pipeline Execution

All pipeline resources are publically available at `s3://bioskryb-public-data/pipeline_resources` users need not have to download this, and will be downloaded during nextflow run.

**Command**

example-

** csv input **

```
git clone https://github.com/BioSkryb/bj-wgs.git
cd bj-wgs
nextflow run main.nf --input_csv $PWD/tests/data/inputs/input.csv --publish_dir results/bj-wgs --max_cpus 8 --max_memory 50.GB
```

**Input Options**

The input for the pipeline can be passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 4 columns: `biosampleName`, `reads`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, `reads` have the number of reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,reads,read1,read2
chr22_testsample1,1000000,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R2_001.fastq.gz
```

- **CSV Metadata Input for Ultima**: The CSV file should have 4 columns: `biosampleName`, `reads`, `cram` and `crai`. 
`cram` and `crai` has the path to the input cram and the index files respectievely. For example:

```
biosampleName,reads,cram,crai
402468-Z0016-Kit4Pl1-DNA-SC04,1000,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/402468-Z0016-Kit4Pl1-DNA-SC04.cram,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/402468-Z0016-Kit4Pl1-DNA-SC04.cram.crai
```

- **Optional Parameters**: Users have the option to provide an additional parameter, 

- `--platform`. This parameter is used to select the sequencing platform. The default platform is Illumina, but users can also choose Ultima or Element.
- `--dnascope_model_selection`. This parameter is used to select the dnascope model for Illumina runs. The default is sentieon, but users can also choose bioskryb129 model for Illumina.
- `--subsample_array`: Specifies a list of read counts for subsampling with seqtk. To enable subsampling, ensure that the `--skip_subsampling` parameter is set to `false`. Provide the read counts as comma-separated values.


**Optional Modules**


This pipeline includes several optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_vcfeval`: Set this to `true` to exclude the vcfeval module. By default, it is set to `true`.
- `--skip_variant_annotation`: Set this to `true` to exclude the Variant Annotation module. By default, it is set to `false`.
- `--skip_gene_coverage`: Set this to `true` to exclude the gene coverage module. By default, it is set to `false`.
- `--skip_ado`: Set this to `true` to exclude the ADO module. By default, it is set to `true`.
- `--skip_subsampling`: Set this to `false` to enable the seqtk subsampling module. By default, it is set to `true`, which excludes the subsampling process.

**Outputs**


The pipeline saves its output files in the designated "publish_dir" directory. The bam files after alignment are stored in the "secondary_analyses/alignment/" subdirectory and the metrics files are saved in the "secondary_analyses/metrics/<sample_name>/" subdirectory.
For details: [BJ-WGS outputs](https://docs.basejumper.bioskryb.com/pipelines/secondary/bj-wgs/1.7.0/docs/#output-files)

**command options**

```
BJ-WGS

Usage:
    nextflow run main.nf [options]

Script Options: see nextflow.config

    
    [required]

    --input_csv                 FILE    Path to input csv file

    --publish_dir               DIR     Path of the output directory

    --genome                    STR     Reference genome to use. Available options - GRCh38, GRCm39
                                DEFAULT: GRCh38

    
    [optional]

    --multiqc_config            DIR     Path to multiqc config
                                        DEFAULT: /home/ubuntu/data/nf-wgs-pipeline/assets/multiqc

    --skip_fastq_merge          BOOL    Skip merging of fastq files
                                        DEFAULT: null

    --skip_gene_coverage        BOOL    Whether to skip gene coverage
                                        DEFAULT: false

    --skip_vcfeval              BOOL    Whether to skip vcfeval
                                        DEFAULT: true

    --skip_ado                  BOOL    Whether to skip ADO
                                        DEFAULT: true

    --skip_variant_annotation   BOOL    Whether to skip variant annotation
                                        DEFAULT: null

    --platform                  STR     Select the sequencing platform. 
                                        Available options - Illumina, Ultima, Element
                                        DEFAULT: Illumina

    --dnascope_model_selection  STR     Select the dnascope model for Illumina runs.
                                        DEFAULT: sentieon

    --subsample_array           LIST    Specifies a list of read counts for subsampling with seqtk. 
                                        To enable subsampling, ensure that the --skip_subsampling parameter is set to false. 
                                        Provide the read counts as comma-separated values.

    --help                      BOOL    Display help message

```

**Tool versions**

- `Seqtk: 1.3-r106`
- `Sentieon: 202308.01`
- `SnpEff: SnpEff 5.1d 2022-04-19`
- `VCFeval: 3.12.1`
- `BCFtools: 1.14`

**nf-test**


The BioSkryb BJ-WGS nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
sudo mv nf-test /usr/local/bin/
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("WGS_test") {

        when {
            params {
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                {assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()},

                // Check for a match in the all metrics MQC text file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/metrics/nf-wgs-pipeline_all_metrics_mqc.txt")).match("all_metrics_mqc")},

                // Check for a match in the selected metrics MQC text file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/metrics/nf-wgs-pipeline_selected_metrics_mqc.txt")).match("selected_metrics_mqc")},

                // Check for a match in the sentieonmetrics text file
                // {assert snapshot (path("${outputDir}/results_test/secondary_analyses/alignment/chr22_testsample1.dedup_sentieonmetrics.txt")).match("sentieonmetrics")},

                // Verify existence of the bam file
                {assert new File("${outputDir}/results_test/secondary_analyses/alignment/chr22_testsample1.bam.bai").exists()},

                // Verify existence of the dnascope vcf file
                {assert new File("${outputDir}/results_test/secondary_analyses/variant_calls_dnascope/chr22_testsample1_dnascope.vcf.gz.tbi").exists()},

                // Verify existence of the snpEff file
                {assert new File("${outputDir}/results_test/tertiary_analyses/variant_annotation/snpeff/chr22_testsample1_snpEff.csv").exists()},

                // Verify existence of the snpEff version file
                {assert new File("${outputDir}/results_test/tertiary_analyses/variant_annotation/snpeff/snpEff_version.yml").exists()}

            )
        }

    }
```

# Need Help?
If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References
For more information, you can refer to the following publications:

- Chung, C., Yang, X., Hevner, R. F., Kennedy, K., Vong, K. I., Liu, Y., Patel, A., Nedunuri, R., 
Barton, S. T., Noel, G., Barrows, C., Stanley, V., Mittal, S., Breuss, M. W., Schlachetzki,
J. C. M., Kingsmore, S. F., & Gleeson, J. G. (2024). Cell-type-resolved mosaicism 
reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392. [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V, Shao, D. D., 
Walsh, C. A., Yang, H. W., Johnson, M. D., & Park, P. J. (2024).
High-resolution 
detection of copy number alterations in single cells with HiScanner. BioRxiv, 
2024.04.26.587806. [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

- Zawistowski, J. S., Salas-González, I., Morozova, T. V, Blackinton, J. G., Tate, T., Arvapalli, 
D., Velivela, S., Harton, G. L., Marks, J. R., Hwang, E. S., Weigman, V. J., & West, J. A. 
A. (n.d.). Unifying genomics and transcriptomics in single cells with ResolveOME amplification
chemistry to illuminate oncogenic and drug resistance mechanisms. [https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full](https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of the standard 
quality control processes implemented through ResolveServices<sup>SM</sup>. While these 
pipelines may not be explicitly cited, they are integral to the methodologies 
described.

