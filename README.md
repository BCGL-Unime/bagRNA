# bagRNA
The bagRNA pipeline is an exhaustive utility developed for the genome annotation of eukaryotic organisms which does not require extensive expertise in computational biology from the user.
It is ncRNA-aware, meaning that it retains all the information of non-coding sequences present in sequencing reads, integrating them in the final gff3 annotation file and renames the ncRNA
feature by its predicted biological function (lncRNA, tRNA, rRNA etc.). This pipeline integrates _ab initio_ prediction, transcript alignment and homology-based approaches. 
The modular workflow requires minimal file input - only a genome assembly, protein evidence and RNA-seq reads - and integrates essential bioinformatics processes such as transcriptome assembly,
gene model selection and correction, and structural annotation, and functional annotation all in one go.

# Pipeline
Briefly, the pipeline consists of 5 main modules: 
1) Data preparation
2) Transcriptome assembly
3) Gene models selection and correction
4) Structural annotation
5) Functional annotation

# Installation
First, clone the repository
$ git clone https://github.com/BCGL-Unime/bagRNA.git your_path

Create a conda environment and install mamba (much faster)

$ conda create -n bagRNA mamba -y
$ conda activate bagRNA

Install packages using the config file bagRNA.yaml provided in the repository

$ mamba install -f bagRNA.yaml

Run the script download_docker_images.sh to download the docker images of the tools needed to run the pipeline

$ chmod +x download_docker_images.sh
$ bash download_docker_images.sh

# Database configuration
Databases need to be configured for functional annotation by running the script download_databases_bagRNA.sh (run after download_docker_images.sh)

$ chmod +x download_databases_bagRNA.sh
$ bash download_databases_bagRNA.sh

# User-required installations (not mandatory, but recommended)
Some tools or models like GeneMark, phobius and signalp6 models cannot be freely distributed because of licencing restrictions, thus if the user wishes to use them they need to be installed separately

https://github.com/gatech-genemark/GeneMark-ETP

https://software.sbc.su.se/phobius.html

https://services.healthtech.dtu.dk/services/SignalP-6.0/  (Download the fast or slow-sequential models)


# Running the pipeline

bagRNA version 1.0.0

Usage:         bagRNA <arguments>

Description:   The bagRNA annotation pipeline takes as input a genome, proteins and RNA-seq short reads (from the STAR tsv manifest file)
               to generate a full annotation of protein-coding and non-coding genes together with their functions.
               The user is required to provide configuration files specified in the mikado documentation and a submission template file 
               (examples are present in the bagRNA GitHub page https://github.com/BCGL-Unime/bagRNA)
               Written and developed by Gabriele Rigano - Bioinformatics and Computational Genomics LAB - University of Messina - gabrielerigano99@gmail.com 

  -h, --help             Display this help message and exit

Mandatory inputs:

  --input_fasta <file>         Genome FASTA file (better if softmasked)
  --prot_evidence <file>       Protein evidence in FASTA format
  --busco_lineage <string>     BUSCO lineage (e.g. sordariomycetes)
  --Conditions <tsv file/s>    STAR manifest TSV files for conditions, an example is provided in the bagRNA github repository (https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf)  
  --mikado_config <file>       Specify the mikado configuration file, an example is provided in the bagRNA github repository (https://mikado.readthedocs.io/en/stable/Tutorial/)
  --scoring <string>           Scoring file for mikado (e.g. scerevisiae.yaml) (https://mikado.readthedocs.io/en/stable/Tutorial/Scoring_tutorial/)
  --species <string>           Species name, must be written in quotes with underscores instead of spaces (e.g. "Arabidopsis_thaliana")
  --submission_template <file> Sbt submission template .sbt, an example is provided in the bagRNA github repository

Choose one flag among:

  --helixer_lineage <string>   Lineage for helixer prediction (requires GPU and NVIDIA Container Toolkit configuration) (choose among: fungi, land_plant, vertebrate, invertebrate)
  --Helixer_gff <file>         Input GFF file from precomputed Helixer run (recommended) (https://www.plabipd.de/helixer_main.html)

Additional inputs:

  --lifted_annotation <file>   Input GFF lifted annotation from liftoff/lifton
  --GeneMark_PATH <path>       Path to GeneMarK-ET/ETP (where the gmes_petap.pl executable is)
  --databases <path>           Path to databases for functional annotation (databases can be downloaded with download_databases_bagRNA.sh script present in the bagRNA github repository

Performance and Miscellaneous options:

  --threads, -t <int>          Number of threads (default: 20)
  --jaccard_clip               Use this flag if you are expecting high gene density with UTR overlap (recommended for fungi) (default: uses --jaccard-clip)
  --max_gene_length <int>      Maximum length a gene model can be (default: 30000 bp)
  --RAM_limit_Trinity <int>    RAM limit for Trinity (default: 45G)
  --limitBAMsortRAM <int>      STAR BAM sort RAM limit
  --orientation <string>       Reads orientation (e.g. FR, RF)
  --strandedness <string>      Strandedness type (e.g. secondstrand) (https://chipster.csc.fi/manual/library-type-summary.html)
  --max_intron_length <int>    Max intron length (default: 3000)
  --codon_table <int>          Codon table (default: 1)
  --strain <string>            Strain/Isolate name (default: strain)
  --locus_tag <tag>            Locus tag (default: bagRNA)
  --no_functional_anno         Use this flag if you do NOT wish to run functional annotation (default: runs functional annotation)
  --no_antismash               Use this flag if you do NOT wish to run AntiSMASH (default: runs AntiSMASH)
  --no_tmbed                   Use this flag if you do NOT wish to run Tmbed (default: runs Tmbed)
| Argument                | Description                                             |
| ----------------------- | ------------------------------------------------------- |
| `--input_fasta`         | Genome FASTA file (preferably softmasked)               |
| `--prot_evidence`       | Protein evidence in FASTA format                        |
| `--busco_lineage`       | BUSCO lineage (e.g., `sordariomycetes`)                 |
| `--Conditions`          | STAR manifest TSV file(s)                               |
| `--mikado_config`       | Mikado configuration YAML file                          |
| `--scoring`             | Mikado scoring config (e.g., `scerevisiae.yaml`)        |
| `--species`             | Species name in quotes (e.g., `"Arabidopsis_thaliana"`) |
| `--submission_template` | `.sbt` file for GenBank submission                      |


