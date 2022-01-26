# Puppets

Hello, we are Puppets!

Collection of pipelines for various high-throughput sequencing data.



## BAMtoFASTQ
Generating FASTQ files from BAM files by steps as follows:
- dealignment ([Samtools](http://www.htslib.org))
- splitting into lanes 

Please note that original base quality cannot be restore if base recalibration was involved in processing of the BAM file.



## EpiCons
Generating of epigenomic consensus (from chromatin state data) for specific cell type (subset of provided bed files) by steps as follows:
- filtering bed files for specific chromatin state(s)
- merging and segmenting filtered bed files ([bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/multiinter.html))
- filtering segments by minimal number of samples
- sorting ([bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/sort.html)) and joining segments ([bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html))
- filtering joined segments by minimal length (or top 5%)

Pipeline generated with BLUEPRINT data from [Carrillo-de-Santa-Pau et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28934481/) and used for genome-wide epigenomic consensus in [Mikulasova et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34933939/)



## FASTQtoBAM
Pre-processing of paired-end read DNA sequencing data, from FASTQ to BAM file by steps as follows:
- merging FASTQ files for R1 and R2
- mapping to the reference ([bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml))
- alignments sorting ([Picard](https://broadinstitute.github.io/picard/))
- duplicates marking ([Picard](https://broadinstitute.github.io/picard/))
- indel realignment ([GATK v3.8](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk))
- BQSR ([GATK v4](https://github.com/broadinstitute/gatk/releases))
- duplicates removing ([Samtools](http://www.htslib.org); muted)
- QC: [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MULTIQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)



## InSilicoWGS
Simulation of DNA abnormalities in whole-genome sequencing data, from VCF to FASTQ by steps as follows:
- filtering reference file for specific contigs (required to handle gonozomes correctly, to make 46,XX or 46,XY and not 48,XXYY)
- introducing abnormalities in reference file ([SimuG](https://github.com/yjx1217/simuG))
- generating FASTQ files ([ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm))

This particular example is to generate paired-end read sequencing from HiSeq 2000 and with read length 100bp, depth 30x and fragments size 500pb for:
- SAMPLE001: male with heterozygous losses on autosomes and hemizygous on gonozomes
- SAMPLE002: male with heterozygous gains on autosomes and hemizygous on gonozomes
- CONTROL001: wild type male
- CONTROL011: wild type female

For other abnormalities to be generated and VCF format see [SimuG](https://github.com/yjx1217/simuG) and for simulation with different parameters like sequencer, depth etc. see [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). Examples of VCF files were preparing by filtering [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) database.




## SomaticSNV
Somatic SNV calling using ([GATK v4](https://github.com/broadinstitute/gatk/releases)) by steps as follows:
- Mutect2 collection for controls
- control panel preparation
- Mutect2 collection for samples
- filtering variant based on contamination and sequencing artefact



## SRAtoPEAKS
Getting data from [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) database and processing ChiPseq of histone marks by steps as follows:
- downloading publicly available data from SRA database and converting into FASTQ files ([sra-tools](https://github.com/ncbi/sra-tools))
- data pre-processing (FASTQ to BAM)
  - mapping to the reference ([bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml)) - single-end sequencing or end-paired read sequencing (muted)
  - alignments sorting ([Picard](https://broadinstitute.github.io/picard/))
  - duplicates marking ([Picard](https://broadinstitute.github.io/picard/))
  - filtering unmapped and low quality alignments ([Samtools](http://www.htslib.org))
  - duplicates removing ([Samtools](http://www.htslib.org))
  - QC: [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MULTIQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- peak calling
  - modelling fragments size ([PhantomPeakQualTools](https://github.com/crazyhottommy/phantompeakqualtools))
  - narrow/broad peak calling ([MACS2](https://pypi.org/project/MACS2/))
  - converting bedGraph to BigWig ([bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/)

Pipeline was develop to process ChIPseq data based on [BLUEPRINT](http://dcc.blueprint-epigenome.eu/#/md/methods) method. It processes data from [GSE54379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54379) and [GSE65687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65687) dataset for T-ALL cell lines KOPT-K1, DND-41 and Jurkat. It was used to process T-ALL data in [Mikulasova et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34933939/)



