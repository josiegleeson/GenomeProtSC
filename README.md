# GenomeProtSC: an integrated proteogenomics analysis platform for long-read single-cell data


## Contents

- [Installation](#installation)
- [General usage](#general-usage)
  - [Database generation](#1-database-generation)
  - [Proteomics](#2-proteomics)
  - [Integration](#3-integration)
  - [Visualisation](#4-visualisation)
- [Detailed input and output descriptions](#detailed-input-and-output-descriptions)
  - [Database generation](#1-database-generation-1)
  - [Proteomics](#2-proteomics-1)
  - [Integration](#3-integration-1)
  - [Visualisation](#4-visualisation-1)

## Installation
Note: Option 1 and 2 are under development and not yet available.

### Option 1 (recommended): Run the shiny application with Docker
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```
docker run --rm -p 3838:3838 josieg/genomeprotsc:v1
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it’s stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.

### Option 2 (recommmended for downstream analysis only): Access GenomeProtSC online
https://genomeprotSC.researchsoftware.unimelb.edu.au/

### Option 3: Locally install the shiny application

The application has substantial dependencies. 

Clone this repository:
```
git clone https://github.com/josiegleeson/GenomeProtSC.git
```

Unzip the uniprot+openprot reference file in the GenomeProtSC/data directory.
```
cd GenomeProtSC/data
unzip openprot_uniprotDb_hs.txt.zip
unzip *
```

Place the reference genome FASTA in the GenomeProtSC/refs directory and re-name (human.fa or mouse.fa).
```
cd GenomeProtSC/refs
cp human.reference.genome.fa .
mv human.reference.genome.fa human.fa
```

## General usage 

GenomeProtSC is an integrated proteogenomics platform with four modules: 1) database generation, 2) proteomics (under development), 3) integration, and 4) visualisation.

### 1. Generate database

The first module processes long-read single-cell data with FLAMES and generates a custom proteome database to perform proteomics searches. The module accepts RNA sequencing FASTQ files from 10X single-cell long-read sequencing. The main output from this module is the results from FLAMES (single-cell transcriptomics), a FASTA file with candidate protein sequences and a metadata file with details of each candidate protein. 

GenomeProtSC currently supports open reading frame (ORF) identification and database generation for humand and mouse data. Users can specify an option to include short upstream ORF (uORF) and downstream ORF (dORF) protein sequences >10 amino acids (AA). Protein sequences are generated based on a user defined minimum length set to >30 AA by default.

#### Inputs:

- FASTQ file(s) (one per sample, can be gzipped)
- Reference GENCODE annotation GTF file
- Reference genome FASTA file

#### Outputs:

- FLAMES output
- Seurat output and gene/transcript counts
- FASTA file with protein sequences
- Metadata file detailing each protein entry

### 2. Proteomics

Currently under development. Please use FragPipe to process proteomics data with the custom database.

### 3. Integration

This module integrates proteomics and transcriptomics data. Peptides are associated back to transcript isoforms and mapped to spliced genomic coordinates for downstream visualisation. This generates BED12 file of transcripts, ORFs and peptides for visualisation in the UCSC genome browser and a combined GTF file for visualisation within the app. An html report is also created that provides a summary of identified known and novel transcripts, uniquely mapping peptides, and known and novel ORFs.

#### Inputs:

- Database files generated in Module 1
- Proteomics peptide data

#### Outputs:

- Peptide-to-transcript mappings with spliced genomic coordinates (CSV)
- BED12 files for visualisation in UCSC genome browser
- HTML report summarising identified transcripts, peptides and ORFs

### 4. Visualisation

The visualisation module uses ggtranscript to generate peptide mapping plots along transcript isoforms with quantitative peptide intensities and transcript expression data (median per cell cluster). This allows users to visualise transcript and peptide abundance across different experimental conditions. This module requires the combined GTF file generated in the integration step, and optionally inputs single-cell transcript counts from Module 1 and peptide intensities from external proteomics analysis. We have included gene filtering options to quickly search for features of interest, such as genes with unknown functions.

#### Inputs:

- GTF from Module 3
- Single-cell transcript counts and peptide intensities (optional)

#### Output:

- Peptide mappings along transcript isoforms (with optional quantitative heatmaps)
- Allows export of plots as PDFs


## Detailed input and output descriptions

### 1. Generate database

| Input  | File Type | Required? | Description     |
|------------------------------------|-----------|-----------|-----------------------------------------------------------------------------------------------|
| Sequencing data  | FASTQ(s)  | Yes  | Mass spec files  | Long-read single-cell sequencing files (one per sample)
| Reference annotations	 | GTF | Yes  | ENSEMBL or Gencode annotation  |


| Output  | File   | File Type | Description      |
|-------------------------------------|------------------------------------|------------|---------------------------------------------------------------------------------------------|
| Database   | proteome_database.fasta  | FASTA | Amino acid sequences of all ORFs in the data  |
| Database metadata | proteome_database_metadata.txt | TXT  | Information on each ORF in the data    |
| Database transcripts | proteome_database_transcripts.gtf | GTF  | Annotations of transcripts used to generate the database   |
| Single-cell gene counts | gene_counts.txt  | TXT  | Gene count file with cells as columns and transcripts as rows  |
| Single-cell transcript counts | transcript_counts.txt  | TXT  | Transcript count file with cells as columns and transcripts as rows  |
| Seurat object	| seurat_object.rds	| RDS	| Seurat object that can be loaded into R with ‘readRDS’ |
| Seurat cell clusters	| sample_cellbarcode_cellcluster.txt	| TXT	| Metadata file with clusters per cell per sample |
| UMAP plot	| UMAP.pdf	| PDF	| UMAP plot coloured by cell clusters |

Proteome FASTA examples and header formats: 
```
>protein_accession|CO=genomic_coordinates GA=gene_accession GN=gene_name TA=transcript_accession
MCGNNMSAPMPAVVPAARKATAAVIFLHGLGDTGHGWAEAFAGIKSPHIKYICPHAPVMPVTLNMNMAMPSWFDIVGLSPDSQEDESGIKQAAETVKALIDQEVKNGIPSNRIILGGFSQGPINSANRDISVLQCHGDCDPLVPLMFGSLTVERLKALINPANVTFKIYEGMMHSSCQQEMMDVKHFIDKLLPPID
>P10711|CO=chr1:4928137-4966584 GA=ENSMUSG00000033813.16 GN=Tcea1 TA=ENSMUST00000081551.14
MEDEVVRIAKKMDKMVQKKNAAGALDLLKELKNIPMTLELLQSTRIGMSVNALRKQSTDEEVTSLAKSLIKSWKKLLDGPSTDKDPEEKKKEPAISSQNSPEAREESSSSSNVSSRKDETNARDTYVSSFPRAPSTSDSVRLKCREMLAAALRTGDDYVAIGADEEELGSQIEEAIYQEIRNTDMKYKNRVRSRISNLKDAKNPNLRKNVLCGNIPPDLFARMTAEEMASDELKEMRKNLTKEAIREHQMAKTGGTQTDLFTCGKCKKKNCTYTQVQTRSADEPMTTFVVCNECGNRWKFC
>ORF_3|CO=chr2:53029193-53081430 GA=ENSMUSG00000061136.17 GN=Prpf40a TA=ENSMUST00000209364.3
MQATPSEAGGESPQSCLSVSRSDWTVGKPVSLLAPLIPPRSSGQPLPFGPGGRQPLRSLLVGMCSGSGRRRSSLSPTMRPGTGAERGGLMMGHPGMHYAPMGMHPMGQRANMPPVPHGMMPQMMPPMGG
```
**Note:** Unannotated ORFs are denoted by "ORF_" and variant proteins by “mORF_” followed by a unique number. UniProt or RefSeq accessions are retained for annotated proteins.

#### Open reading frame (ORF) category definitions:

| Type  | Definition    |
|----------------|-----------------------------------------------------------------------------------------------------|
| CDS  | Annotated in UniProt or RefSeq    |
| 5UTR  | Coordinates are within the 5’ UTR region of an mRNA transcript    |
| 3UTR  | Coordinates are within the 3’ UTR region of an mRNA transcript    |
| 5UTR:CDS  | Start site is within the 5’ UTR region and stop site is within the CDS region of an mRNA transcript |
| gene_overlap  | Encoded by a transcript that overlaps a region with annotated protein-coding genes   |
| intergenic | Encoded by a transcript that does not overlap a region with annotated protein-coding genes   |


### 2. Proteomics

We recommend installing and running [FragPipe](https://github.com/Nesvilab/FragPipe) for analysing mass spectrometry-based proteomics data.

| Input  | File Type | Required? | Description     |
|------------------------------------|-----------|-----------|-----------------------------------------------------------------------------------------------|
| Mass spec data  | mzML, RAW  | Yes  | Mass spec files  |
| Database (proteome_database.fasta) | FASTA | Yes  | Generated in Module 1. Amino acid sequences of all ORFs in the data  |


The output file generated is typically `peptides.txt` or `report.pr_matrix.tsv`.


### 3. Integration

| Input  | File Type | Required? | Description     |
|------------------------------------|-----------|-----------|-----------------------------------------------------------------------------------------------|
| Proteomics peptide data  | TSV/TXT  | Yes  | Peptide results. Typically, 'peptides.txt', ‘peptide.tsv’ or ‘report.pr_matrix.tsv’  |
| Database (proteome_database.fasta) | FASTA | Yes  | Generated in Module 1. Amino acid sequences of all ORFs in the data  |
| Database metadata (proteome_database_metadata.txt) | TXT  | Yes  | Generated in Module 1. Information on each ORF in the data   |
| Database transcripts (proteome_database_transcripts.gtf) | GTF  | Yes  | Generated in Module 1. Annotations of transcripts used to generate the database   |


| Output | File  | File Type | Description    |
|------------------------------|--------------------------|-----------|----------------------------------------------------------------|
| Peptide information   | peptide_info.csv   | CSV  | Main results file with peptide mapping data |
| Report | summary_report.html | HTML | Summary report  |
| Combined annotation data | combined_annotations.gtf | GTF  | Annotations of peptides, ORFs, and transcripts for visualisation |
| Peptide coordinates   | peptides.bed12  | BED12 | Peptide spliced genomic coordinates   |
| ORF coordinates   | ORFs.bed12 | BED12 | ORF spliced genomic coordinates   |
| Transcript coordinates  | transcripts.bed12  | BED12 | Transcript spliced genomic coordinates |


#### Description of `peptides_info.csv` output:

| Column Name | Description     | Class  |
|--------------------------------|----------------------------------------------------------------------------------------|--------------|
| peptide | Peptide sequence   | character  |
| accession  | Protein accession    | character  |
| PID   | protein_accession\|CO=genomic_coordinates (included for compatibility with FASTA header) | character  |
| transcript_id   | ENSEMBL or novel transcript ID     | character  |
| gene_id | ENSEMBL gene ID   | character  |
| gene_name  | Gene name/symbol   | character  |
| strand  | Strand (+ or -)   | character  |
| number_exons | Number of exons spanned by the peptide   | integer |
| transcript_length   | Transcript length (nt)     | integer |
| transcript_biotype   | Transcript biotype from Gencode    | character  |
| simplified_biotype   | Simplified transcript biotype    | character  |
| protein_length   | Protein length (AA)   | integer |
| orf_genomic_coordinates  | Genomic coordinates of ORF    | numeric |
| orf_type | Annotated, unannotated, or variant protein   | character  |
| localisation | ORF location in the genome (see ORF category definitions)  | character  |
| uniprot_status   | Review status in UniProt (reviewed/unreviewed)    | character  |
| openprot_id | OpenProt ID if present     | character  |
| molecular_weight(kDA)   | Molecular weight of protein (KDa)   | numeric |
| isoelectric_point   | Isoelectric point of ORF calculated using pKa scale EMBOSS (Rice et al., 2000)   | numeric |
| hydrophobicity   | Hydrophobicity profile of ORF calculated using Kyte-Doolittle scale (Kyte et al., 1982) | numeric |
| aliphatic_index | Aliphatic index of ORF (Ikai 1980)  | numeric |
| longest_orf_in_transcript | Longest ORF in the transcript (longest within CDS regions for known proteins)   | true/false  |
| peptide_ids_gene | Is peptide uniquely mapped to gene?    | true/false  |
| peptide_ids_orf | Is peptide uniquely mapped to ORF?    | true/false  |
| peptide_ids_transcript   | Is peptide uniquely mapped to transcript?   | true/false  |
| shared_novel_protein_peptide  | Is peptide shared with other novel proteins?  | true/false  |
| orf_identified   | Is ORF identified with unique peptide evidence?  | true/false  |
| gene_identified | Is gene identified with unique peptide evidence?  | true/false  |
| transcript_identified   | Is transcript identified with unique peptide evidence?    | true/false  |


### 4. Visualisation

| Input  | File Type | Required? | Description   |
|--------------------------------|-----------|-----------|---------------------------------------------------------|
| Combined annotations (combined_annotations.gtf)  | GTF  | Yes  | Generated in Module 3, annotations of peptides, ORFs, and transcripts    |
| Single-cell transcript counts (transcript_counts.txt)   | TXT/CSV  | No  | Generated in Module 1, transcript counts per cell barcode    |
| Peptide intensities  | TXT  | No  | Peptide intensity data ‘report.pr_matrix.tsv’  |

**Note:** There is an option to download plots as a PDF.

