
source("global.R")

# define options
option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="custom user GTF", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="reference user GTF", metavar="character"),
  make_option(c("-o", "--organism"), type="character", default=NULL,
              help="Organism", metavar="character"),
  make_option(c("-l", "--length"), type="integer", default=NULL,
              help="Minimum ORF length", metavar="integer"),
  make_option(c("-u", "--uorfs"), type="logical", default=NULL,
              help="Find uORFs", metavar="TRUE/FALSE"),
  make_option(c("-d", "--dorfs"), type="logical", default=NULL,
              help="Find dORFs", metavar="TRUE/FALSE"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# define inputs
gtf_path <- opt$gtf
reference_gtf <- opt$reference
organism <- opt$organism
min_orf_length <- opt$length
find_5_orfs <- opt$uorfs
find_3_orfs <- opt$dorfs
output_directory <- opt$savepath

# define functions

# import gtf and filter for minimum transcript counts
filter_custom_gtf <- function(customgtf, organism, outdir) {
  
  # import bambu gtf
  bambu_data <- rtracklayer::import(customgtf)
  
  # remove scaffolds and weird chromosomes
  bambu_data <- bambu_data[grep("chr", seqnames(bambu_data))]
  
  # remove rows with _ or . in the chromosome
  keep_rows <- !grepl("[._]", seqnames(bambu_data))
  
  # subset the GRanges
  bambu_data <- bambu_data[keep_rows]
  
  # filter based on strand
  okstrand <- c("+", "-")
  bambu_data <- bambu_data[strand(bambu_data) %in% okstrand]

  # convert to tibble
  bambu_df <- bambu_data %>% as_tibble()
  
  # remove version numbers for search
  bambu_df <- bambu_df %>% separate(gene_id, into="ensg_id", sep="\\.", remove = FALSE)
  
  # use mygene to search for gene names
  gene_query <- queryMany(unique(bambu_df$ensg_id), scopes="ensembl.gene", fields="symbol", species=as.character(organism),  returnall=TRUE)
  
  # make df
  gene_df <- as.data.frame(gene_query[["response"]])
  
  # if there was no name found, use original ID
  gene_df_formatted <- gene_df %>% 
    dplyr::mutate(gene_name = case_when(
      is.na(symbol) ~ query,
      !is.na(symbol) ~ symbol
    )) %>% 
    dplyr::select(query, gene_name)
  
  # merge results 
  bambu_merged <- merge(bambu_df, gene_df_formatted, by.x="ensg_id", by.y="query", all.x=T, all.y=F)
  bambu_merged$ensg_id <- NULL
  
  # make GRanges including new names
  bambu_data_gr <- makeGRangesFromDataFrame(bambu_merged,
                                            keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                            seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand",
                                            starts.in.df.are.0based=FALSE)
  
  # separate for sorting
  bambu_exons <- bambu_data_gr[bambu_data_gr$type == "exon"]
  bambu_transcripts <- bambu_data_gr[bambu_data_gr$type == "transcript"]
  
  # sort by chr and locations
  bambu_exons <- sortSeqlevels(bambu_exons)
  bambu_transcripts <- sortSeqlevels(bambu_transcripts)
  
  bambu_exons <- sort(bambu_exons)
  bambu_transcripts <- sort(bambu_transcripts)
  
  # recombine for export
  bambu_export <- c(bambu_transcripts, bambu_exons)
  
  # export filtered gtf
  export(bambu_export, paste0(outdir, "/proteome_database_transcripts.gtf"), format="gtf")
  
  print("Exported filtered GTF")
  
}

# generate FASTA of transcript sequences
get_transcript_orfs <- function(filteredgtf, organism, orf_len=30, find_UTR_5_orfs=FALSE, find_UTR_3_orfs=FALSE, referencegtf, outdir) {
  
  # set organism
  if (organism == "human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genomedb <- BSgenome.Hsapiens.UCSC.hg38
  } else if (organism == "mouse") {
    library(BSgenome.Mmusculus.UCSC.mm39)
    genomedb <- BSgenome.Mmusculus.UCSC.mm39
  }
  
  # required for UTR regions
  ref_txdb <- makeTxDbFromGFF(referencegtf)
  
  # import filtered gtf as a txdb
  txdb <- makeTxDbFromGFF(filteredgtf)
  txs <- exonsBy(txdb, by=c("tx", "gene"), use.names=TRUE)
  
  # convert txdb to GRangesList
  txs_grl <- GRangesList(txs)
  
  # ORFik function to find ORFs in GRangesList
  apply_orfik <- function(grl, orfik_min_length, orfik_max_length, orfik_type) {
    
    # extract transcript nt sequence
    tx_seqs <- extractTranscriptSeqs(genomedb, grl)
    
    # ORFik
    ORFs <- findMapORFs(grl,
                        tx_seqs, 
                        groupByTx = FALSE,
                        longestORF = orfik_type, 
                        minimumLength = as.numeric(orfik_min_length), 
                        startCodon = "ATG",
                        stopCodon = stopDefinition(1))
    
    # unlist GRL
    ORFs_unlisted <- unlist(ORFs) %>% as_tibble()
    
    # add width column, filter for ORFs of defined length
    orf_genome_coordinates <- ORFs_unlisted %>% 
      rowwise() %>% 
      dplyr::mutate(width = end - start) %>% 
      group_by(names) %>% 
      summarise(chr = seqnames[1],
                start = min(start),
                end = max(end),
                length = sum(width),
                strand = strand[1]) %>% 
      ungroup() %>% 
      dplyr::filter(length < (as.numeric(orfik_max_length)*3)-3) %>% # length ORFs < 30 AA
      dplyr::select(-length)
    
    # remove any ORFs from original ORF object if they were filtered out due to length settings above
    ORFs <- ORFs[names(ORFs) %in% orf_genome_coordinates$names]
    
    # rename
    orf_genome_coordinates$ORF_id <- orf_genome_coordinates$names
    orf_genome_coordinates$names <- NULL
    
    # convert these ORF coordinates into nucleotide sequences
    orf_seqs <- GenomicFeatures::extractTranscriptSeqs(genomedb, ORFs)
    
    # convert the nucleotide sequences to amino acid sequences
    orf_aa_seq <- Biostrings::translate(orf_seqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
    
    # create data frame of all possible ORFs
    orf_aa_seq_df <- data.frame(ORF_id = orf_aa_seq@ranges@NAMES,ORF_sequence = orf_aa_seq, row.names=NULL) 
    
    # merge with coordinates
    orf_aa_seq_df_genomic_coordinates <- left_join(orf_aa_seq_df, orf_genome_coordinates, by = "ORF_id")
    
    return(orf_aa_seq_df_genomic_coordinates)
    
  }
  
  # ORF discovery
  
  # set ORF max length to large number (to disable)
  # set longestORF to FALSE to ensure we identify CDS
  all_ORFs <- apply_orfik(txs_grl, as.numeric(orf_len), 100000, FALSE)
  
  # create tmp copy
  tmp <- all_ORFs
  
  # extract 5UTR ORFs
  if (find_UTR_5_orfs == TRUE) {
    
    # extract 5' UTR regions from ref gtf
    utrs5 <- fiveUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs5_filtered <- utrs5[names(utrs5) %in% names(txs)]
    
    # ORF max length is now the main ORF min length
    # only keep longest UTR ORFs
    utr5_ORFs <- apply_orfik(utrs5_filtered, 10, as.numeric(orf_len), TRUE)
    
    # update tmp
    tmp <- rbind(tmp, utr5_ORFs)
    
  }
  
  # extract 3UTR ORFs
  if (find_UTR_3_orfs == TRUE) {
    
    # extract 3' UTR regions from ref gtf
    utrs3 <- threeUTRsByTranscript(ref_txdb, use.names = TRUE)
    utrs3_filtered <- utrs3[names(utrs3) %in% names(txs)]
    
    # ORF max length is now the main ORF min length
    # only keep longest UTR ORFs
    utr3_ORFs <- apply_orfik(utrs3_filtered, 10, as.numeric(orf_len), TRUE)
    
    # update tmp
    tmp <- rbind(tmp, utr3_ORFs)
    
  }
  
  # rename final ORFs with new numerical IDs
  combined <- tmp %>% 
    separate(ORF_id, into=c("transcript_id"), sep=c("\\_"), remove=F) %>% 
    group_by(transcript_id) %>% 
    mutate(tx_id_number = row_number()) %>% 
    ungroup()
  
  # apply new names
  combined$ORF_id <- paste0(combined$transcript_id, "_", combined$tx_id_number)
  combined$transcript_id <- NULL
  combined$tx_id_number <- NULL
  
  # export protein seqs for python script
  write_tsv(combined, paste0(outdir, "/ORFome_aa.txt"))
  
  print("Exported ORFik data")
  
  return(combined)
  
}

# run functions
# filter and process gtf
filtered_gtf <- filter_custom_gtf(customgtf=gtf_path, organism=organism, outdir=output_directory)

# run ORFik to get ORFs
get_transcript_orfs(filteredgtf=paste0(output_directory, "/proteome_database_transcripts.gtf"), organism=organism, orf_len=min_orf_length, find_UTR_5_orfs=find_5_orfs, find_UTR_3_orfs=find_3_orfs, referencegtf=reference_gtf, outdir=output_directory)


 
 

 


