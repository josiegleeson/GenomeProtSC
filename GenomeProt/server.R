
library(shiny)
library(shinyjs)

# internal server functions
flames_server <- function(input, output, session) {
  
  # to start with fastqs a user must upload a genome, reference gtf, fastq files
  req(input$user_reference_gtf$datapath, input$user_fastq_files$datapath)  # required
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_flames <- paste0(session_id, "/singlecell_output")
  # create output dir
  system(paste0("mkdir ", outdir_flames))
  
  # set genome
  if (organism == "human") {
    genome_file <- "refs/human.fa"
  } else if (organism == "mouse") {
    genome_file <- "refs/mouse.fa"
  }
  
  # FLAMES
  command_run_flames <- paste0(
    "Rscript bin/database_module/run-flames.R",
    " -g ", genome_file,
    " -a ", input$user_reference_gtf$datapath,
    " -i ", input$user_fastq_files$datapath, # test
    " -o ", outdir_flames, "/flames"
  )
  
  # post-process with Seurat
  command_run_flames <- paste0(
    "Rscript bin/database_module/flames-postprocess.R",
    " -i ", outdir_flames, "/flames",
    " -o ", outdir_flames, "/seurat"
  )
  

  
}

database_server <- function(input, output, session) {
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_db <- paste0(session_id, "/database_output")
  # create output dir
  system(paste0("mkdir ", outdir_db))

  db_gtf_file <- paste0(session_id, "/singlecell_output/flames/")
  db_counts_file <- paste0(session_id, "/singlecell_output/flames/")
  
  flames_files <- c(paste0(session_id, "/singlecell_output/flames/bambu_transcript_annotations.gtf"), paste0(session_id, "/singlecell_output/flames/transcript_counts.txt"), paste0(session_id, "/singlecell_output/flames/novel_transcript_classes.csv"), paste0(session_id, "/singlecell_output/flames/gffcompare.tmap.txt"))
  
  # construct the command
  command_generate_proteome <- paste0(
    "Rscript bin/database_module/generate_proteome.R",
    " -g ", db_gtf_file,
    " -r ", input$user_reference_gtf$datapath,
    " -c ", db_counts_file, 
    " -m ", input$minimum_tx_count,
    " -o ", input$organism,
    " -l ", input$min_orf_length,
    " -u ", input$user_find_utr_5_orfs,
    " -d ", input$user_find_utr_3_orfs,
    " -s ", outdir_db
  )
  
  # print command
  print(command_generate_proteome)
  
  # run command
  system(command_generate_proteome)
  
  print("Generated ORFs")
  
  # set reference protein database per organism 
  if (input$organism == "human") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "mouse") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  # } else if (input$organism == "celegans") {
  #   ref_proteome <- "data/openprot_uniprotDb_c_elegans.txt"
  # } else if (input$organism == "drosophila") {
  #   ref_proteome <- "data/openprot_uniprotDb_drosophila.txt"
  # } else if (input$organism == "rat") {
  #   ref_proteome <- "data/openprot_uniprotDb_rat.txt"
  # } else if (input$organism == "zebrafish") {
  #   ref_proteome <- "data/openprot_uniprotDb_zebrafish.txt"
  }
  
  # run python script to create proteome fasta
  command_annotate_proteome <- paste0("python bin/database_module/annotate_proteome.py ", input$user_reference_gtf$datapath, " ", 
                                      ref_proteome, " ", outdir_db, "/ORFome_aa.txt ", outdir_db, "/proteome_database_transcripts.gtf ", 
                                      outdir_db, " all ", input$min_orf_length, " None")
  print(command_annotate_proteome)
  system(command_annotate_proteome)
  
  print("Annotated proteome")
  
  # get top level directory
  top_level_dir <- getwd()
  
  # zip all results files depending on input types
  if (file.exists(paste0(outdir_db, "/proteome_database.fasta")) && file.exists(paste0(outdir_db, "/proteome_database_transcripts.gtf")) && !file.exists(paste0(outdir_db, "/orf_temp.txt"))) {
    
    file.copy(bambu_files, outdir_db)
    files_to_zip_db <- c(bam_files_zip, "bambu_transcript_annotations.gtf", "transcript_counts.txt", "novel_transcript_classes.csv", "gffcompare.tmap.txt", "proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
    
    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_db <- file.path("../database_results.zip")
    
    # temp change the working dir to outdir_db
    tmp_wd <- setwd(outdir_db)
    
    # zip files
    zip(zipfile = zipfile_path_db, files = files_to_zip_db)
    
    # change back to starting wd
    setwd(top_level_dir)
    
  }
  
}

proteomics_server <- function(input, output, session) {
  
  # req(input$user_mm_data, input$user_mm_fasta)
  # 
  # # get directory path
  # dir_path <- dirname(input$user_mm_data$datapath[1])
  # print(dir_path)
  # 
  # mass_spec_names <- c(input$user_mm_data$name)
  # print(mass_spec_names)
  # 
  # # new file paths
  # new_file_paths <- file.path(dir_path, mass_spec_names)
  # 
  # # rename files
  # file.rename(input$user_mm_data$datapath, new_file_paths)
  # 
  # # copy config files
  # # replace all lines:
  # # 'MaxThreadsToUsePerFile = 3'
  # # with user set threads
  # 
  # #system(paste0("source activate mm_env; metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  # system(paste0("metamorpheus -t data/mm_configs/Task2-CalibrateTaskconfig.toml data/mm_configs/Task4-GPTMDTaskconfig.toml data/mm_configs/Task5-SearchTaskconfig.toml -s ", dir_path, " -v 'minimal' -d ", input$user_mm_fasta$datapath, " -o proteomics_output"))
  # 
  # # check files exist
  # if (file.exists("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv") && file.exists("proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")) {
  #   # create a zip file with results
  #   files_to_zip <- c("proteomics_output/Task3SearchTask/AllQuantifiedPeptides.tsv", "proteomics_output/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
  #   zipfile_path <- "proteomics_output/proteomics_results.zip"
  #   zip(zipfile = zipfile_path, files = files_to_zip)
  # }
}

integration_server <- function(input, output, session) {
  
  req(input$user_proteomics_file, input$user_post_gtf_file, input$user_fasta_file, input$user_metadata_file)  # GTF is required
  
  # store session ID
  session_id <- session$token
  # set output dir
  outdir_integ <- paste0(session_id, "/integ_output")
  # create output dir
  system(paste0("mkdir ", outdir_integ))
  
  # run Rscript
  system(paste0("Rscript bin/integration_module/map_peptides_generate_outputs.R -p ", input$user_proteomics_file$datapath, " -f ", input$user_fasta_file$datapath, " -m ", input$user_metadata_file$datapath, " -g ", input$user_post_gtf_file$datapath, " -s ", outdir_integ))
  
  # get the top level dir
  top_level_dir <- getwd()
  
  system(paste0("mkdir ", outdir_integ, "/report_images"))
  
  # create report
  rmarkdown::render(input = paste0(top_level_dir, "/bin/integration_module/integration_summary_report.Rmd"),
                    output_file = paste0(top_level_dir, "/", outdir_integ, "/summary_report.html"),
                    output_format = "html_document",
                    params = list(
                      directory = paste0(top_level_dir, "/", outdir_integ),
                      file = "peptide_info.csv"
                    ))
  
  # zip all results files
  if (file.exists(paste0(outdir_integ, "/peptide_info.csv")) && file.exists(paste0(outdir_integ, "/summary_report.html"))) {
    
    # create a zip file with results
    files_to_zip_int <- c("summary_report.html", "peptide_info.csv", "report_images/",
                          "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                          "peptides.bed12", "ORFs.bed12", "transcripts.bed12")
    
    # set the path to the ZIP file (in the session_id directory)
    zipfile_path_int <- file.path("../integration_results.zip")
    
    # temp change the working dir to outdir_integ
    tmp_wd <- setwd(outdir_integ)
    
    # zip files
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
    
    # go back to starting dir
    setwd(top_level_dir)
    
  }
  
}

# main shiny app server
server <- function(input, output, session) {
  
  # store session ID
  # create session id tmp directory each time app is run
  session_id <- session$token
  print(paste0("Session: ", session_id))
  # create the dir
  system(paste0("mkdir ", session_id))
  
  # DATABASE MODULE
  
  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)
  
  # run database function when submit is pressed
  observeEvent(input$db_submit_button, {
    
    # ensure download button remains greyed out (if submit is re-pressed)
    shinyjs::disable("db_download_button")
    shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#d3d3d3';")
    # disable submit button after it is pressed
    session$sendCustomMessage("disableButton", list(id = "db_submit_button", spinnerId = "db-loading-container"))
    
    fastq_server(input, output, session)
    bambu_server(input, output, session)
    database_server(input, output, session)
   
    # check if the zip file is created
    if (file.exists(paste0(session_id, "/database_results.zip"))) {
      file_available_db(TRUE)
    }
    
  })
  
  # enable download once files are available
  observe({
    if (file_available_db()) {
      shinyjs::enable("db_download_button")
      shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "db_submit_button", spinnerId = "db-loading-container")) # re-enable submit button
    }
  })
  
  # download handler for the database results.zip file
  output$db_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_database_results.zip")
    },
    content = function(file) {
      file.copy(paste0(session_id, "/database_results.zip"), file)
    }
  )
  
  # END DATABASE MODULE
  
  
  # PROTEOMICS MODULE
  
  # # create reactive value for the database zip
  # file_available_mm <- reactiveVal(FALSE)
  # 
  # 
  # # run database function when submit is pressed
  # observeEvent(input$proteomics_submit_button, { 
  #   session$sendCustomMessage("disableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # disable submit button
  #   proteomics_server(input, output, session)
  #   
  #   # check if the zip file is created
  #   if (file.exists("proteomics_output/proteomics_results.zip")) {
  #     file_available_mm(TRUE)
  #   }
  # })
  # 
  # # enable download once files are available
  # observe({
  #   if (file_available_mm()) {
  #     shinyjs::enable("proteomics_download_button")
  #     shinyjs::runjs("document.getElementById('proteomics_download_button').style.backgroundColor = '#4CAF50';")
  #     session$sendCustomMessage("enableButton", list(id = "proteomics_submit_button", spinnerId = "proteomics-loading-container")) # re-enable submit button
  #   }
  # })
  # 
  # # download handler for the database results.zip file
  # output$proteomics_download_button <- downloadHandler(
  #   filename = function() {
  #     paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_proteomics_results.zip")
  #   },
  #   content = function(file) {
  #     file.copy("proteomics_output/proteomics_results.zip", file)
  #   }
  # )
  
  # END PROTEOMICS MODULE
  
  
  # INTEGRATION MODULE
  
  # create reactive value for the database zip
  file_available_integ <- reactiveVal(FALSE)
  
  # run integration function when submit is pressed
  observeEvent(input$integ_submit_button, { 
    
    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button
    
    # run integration server
    integration_server(input, output, session)
    
    # check if the zip file is created
    if (file.exists(paste0(session_id, "/integration_results.zip"))) {
      file_available_integ(TRUE)
    }
    
  })
  
  # enable download once files are available
  observe({
    if (file_available_integ()) {
      shinyjs::enable("integ_download_button")
      shinyjs::runjs("document.getElementById('integ_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # re-enable submit button
    }
  })
  
  # download handler 
  output$integ_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_integration_results.zip")
    },
    content = function(file) {
      file.copy(paste0(session_id, "/integration_results.zip"), file)
    }
  )
  
  # END INTEGRATION MODULE
  
  
  # VISUALISATION MODULE
  
  data_storage <- reactiveValues()
  
  observeEvent(input$toggle_filters, {
    # check status
    filters_visible <- shinyjs::toggle(id = "filters_container")
    
    # update the toggle button text
    if (input$toggle_filters %% 2 == 1) { # odd number clicks
      updateActionLink(session, "toggle_filters", label = "▲ Hide filtering options")
    } else { # even number clicks
      updateActionLink(session, "toggle_filters", label = "▼ Gene list filtering options")
    }
  })
  
  # function to import and process combined gtf
  import_and_preprocess_gtf <- function(gtf_path) {
    # import gtf and remove gene_id version number
    gtf_import <- rtracklayer::import(gtf_path, format = "gtf") %>%
      as_tibble() %>%
      separate(gene_id, into = c("gene_id"), sep = "\\.")
    # separate gtf into multiple features
    list(
      gtf_import = gtf_import,
      res_tx_import = gtf_import %>% dplyr::filter(group_id == "transcripts"),
      res_ORF_import = gtf_import %>% dplyr::filter(group_id == "ORFs"),
      res_pep_import = gtf_import %>% dplyr::filter(group_id == "peptides")
    )
  }
  
  # function to update gene list
  update_gene_list <- function(res_tx_import, res_pep_import, uniq_map_peptides = FALSE, lncRNA_peptides = FALSE, novel_txs = FALSE, novel_txs_distinguished = FALSE, unann_orfs = FALSE, uorf_5 = FALSE, dorf_3 = FALSE) {
    
    # if check box is ticked
    if (uniq_map_peptides) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (lncRNA_peptides) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "lncRNA" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (novel_txs) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "novel" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (novel_txs_distinguished) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(transcript_biotype == "novel" & peptide_ids_orf == TRUE & peptide_ids_transcript == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (unann_orfs) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(orf_type == "unannotated" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (uorf_5) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(localisation == "5UTR" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    if (dorf_3) {
      filtered_peptides <- res_pep_import %>% dplyr::filter(localisation == "3UTR" & peptide_ids_orf == TRUE)
      res_tx_import <- res_tx_import %>% dplyr::filter(gene_id %in% filtered_peptides$gene_id)
    }
    
    # ensure gene_id overlaps
    ensembl_ids <- intersect(res_pep_import$gene_id, res_tx_import$gene_id)
    # filter for list
    genes_list <- res_tx_import %>% dplyr::filter(gene_id %in% ensembl_ids)
    
    # use gene name instead of ensembl
    if ("gene_name" %in% colnames(genes_list)) {
      unique(genes_list$gene_name)
    } else {
      unique(genes_list$gene_id)
    }
  }
  
  # function to update gene selection
  update_gene_selector <- function(session, genes_available) {
    updateSelectInput(session, "gene_selector", choices = genes_available)
  }
  
  # observe integration module file availability
  observe({
    # if the integration module has been run, use the output from this to display visualisation automatically
    if (file_available_integ()) {
      gtf_data <- import_and_preprocess_gtf(paste0(session_id, "/integ_output/combined_annotations.gtf"))
      data_storage$gtf_import <- gtf_data$gtf_import
      data_storage$res_tx_import <- gtf_data$res_tx_import
      data_storage$res_ORF_import <- gtf_data$res_ORF_import
      data_storage$res_pep_import <- gtf_data$res_pep_import
      
      genes_available <- update_gene_list(data_storage$res_tx_import, data_storage$res_pep_import, input$uniq_map_peptides, input$lncRNA_peptides, input$novel_txs, input$novel_txs_distinguished, input$unann_orfs, input$uorf_5, input$dorf_3)
      update_gene_selector(session, genes_available)
    }
  })
  
  # observe visualisation submit button
  observeEvent(input$vis_submit_button, {
    # if the submit button is pressed, import files
    session$sendCustomMessage("disableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container"))
    
    req(input$user_vis_gtf_file)
    
    gtf_data <- import_and_preprocess_gtf(input$user_vis_gtf_file$datapath)
    data_storage$gtf_import <- gtf_data$gtf_import
    data_storage$res_tx_import <- gtf_data$res_tx_import
    data_storage$res_ORF_import <- gtf_data$res_ORF_import
    data_storage$res_pep_import <- gtf_data$res_pep_import
    
    # if there are transcript counts and peptide intensities uploaded
    if (!is.null(input$user_vis_tx_count_file) & !is.null(input$user_pep_count_file)) {
      print("Counts detected")
      
      # import counts files
      data_storage$countst <- fread(input$user_vis_tx_count_file$datapath)
      data_storage$countsp <- fread(input$user_pep_count_file$datapath)
      
      # rename as per bambu counts output
      if ("TXNAME" %in% colnames(data_storage$countst) & "GENEID" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
        data_storage$countst$GENEID <- NULL
      } else if ("TXNAME" %in% colnames(data_storage$countst)) {
        data_storage$countst$transcript_id <- data_storage$countst$TXNAME
      }
      
      # filter GTF transcripts for those with counts
      data_storage$res_tx_import <- data_storage$res_tx_import %>% 
        dplyr::filter(transcript_id %in% data_storage$countst$transcript_id)
      
      # filter GTF ORFs for those with counts
      data_storage$res_ORF_import <- data_storage$res_ORF_import %>% 
        dplyr::filter(transcript_id %in% data_storage$countst$transcript_id)
      
      # match samples in both counts files
      sample_names <- intersect(colnames(data_storage$countsp), colnames(data_storage$countst))
      
      print("Samples with peptide intensities and transcript counts:")
      print(sample_names)
      
      # order sample names
      sample_names <- sample_names[order(match(sample_names, colnames(data_storage$countsp)))]
      
      # rename as per bambu counts output
      if ("Stripped.Sequence" %in% colnames(data_storage$countsp)) {
        data_storage$countsp$Peptide <- data_storage$countsp$Stripped.Sequence
        data_storage$countsp$Stripped.Sequence <- NULL
      } else if ("peptide" %in% colnames(data_storage$countsp)) {
        data_storage$countsp$Peptide <- data_storage$countsp$peptide
        data_storage$countsp$peptide <- NULL
      }
      
      # rarely, a peptide is in the data twice, so take max count value and get unique peptide IDs
      data_storage$countsp <- data_storage$countsp %>% 
        dplyr::select(Peptide, sample_names) %>% 
        dplyr::mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% 
        dplyr::group_by(Peptide) %>% 
        slice_max(sum) %>% dplyr::ungroup() %>% dplyr::select(-sum)
      
      # convert to matrix
      countsp_matrix <- as.matrix(data_storage$countsp[,-1])
      
      # set rownames as peptide IDs
      rownames(countsp_matrix) <- data_storage$countsp$Peptide
      
      # run VSN for normalisation of peptide intensities
      vsnp <- as.data.frame(justvsn(countsp_matrix))
      
      # create column of peptide IDs
      vsnp$peptide <- rownames(vsnp)
      
      # melt peptide intensities for plotting
      data_storage$countspm <- reshape2::melt(vsnp, id.vars = c("peptide"),
                                              variable.name = "sample_id", value.name = "count")
      
      # set levels for plotting
      data_storage$countspm$sample_id <- factor(as.character(data_storage$countspm$sample_id), level =  sample_names)
      
      # melt transcript counts for plotting
      data_storage$countstm <- reshape2::melt(data_storage$countst, id.vars = c("transcript_id"),
                                              variable.name = "sample_id", value.name = "count")
      # set levels for plotting
      data_storage$countstm$sample_id <- factor(as.character(data_storage$countstm$sample_id), level =  sample_names)
    }
    
    genes_available <- update_gene_list(data_storage$res_tx_import, data_storage$res_pep_import, input$uniq_map_peptides, input$lncRNA_peptides, input$novel_txs, input$unann_orfs, input$uorf_5, input$dorf_3)
    update_gene_selector(session, genes_available)
    
    # re-enable submit button after data is processed
    session$sendCustomMessage("enableButton", list(id = "vis_submit_button", spinnerId = "vis-loading-container"))
  })
  
  # make reactive list of checkbox filters
  checkbox_filters <- reactive({
    list(
      uniq_map_peptides = input$uniq_map_peptides,
      lncRNA_peptides = input$lncRNA_peptides,
      novel_txs = input$novel_txs,
      novel_txs_distinguished = input$novel_txs_distinguished,
      unann_orfs = input$unann_orfs,
      uorf_5 = input$uorf_5,
      dorf_3 = input$dorf_3
    )
  })
  
  # observe changes in any checkbox
  observe({
    # ensure GTFs are available
    req(data_storage$res_tx_import, data_storage$res_pep_import)
    
    # update gene list based on current checkbox states
    genes_available <- update_gene_list(
      data_storage$res_tx_import, 
      data_storage$res_pep_import, 
      uniq_map_peptides = checkbox_filters()$uniq_map_peptides,
      lncRNA_peptides = checkbox_filters()$lncRNA_peptides,
      novel_txs = checkbox_filters()$novel_txs,
      novel_txs_distinguished = checkbox_filters()$novel_txs_distinguished,
      unann_orfs = checkbox_filters()$unann_orfs,
      uorf_5 = checkbox_filters()$uorf_5,
      dorf_3 = checkbox_filters()$dorf_3
    )
    
    # update the gene selector with the new list
    update_gene_selector(session, genes_available)
  })
  
  # gene selection in drop down menu
  observeEvent(input$gene_selector, {
    
    req(input$gene_selector)
    
    # set gene name
    data_storage$gene_to_plot <- input$gene_selector
    print(data_storage$gene_to_plot)
    
    if (!is.null(input$user_vis_tx_count_file) & !is.null(input$user_pep_count_file)) { # if counts files are provided 
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import, data_storage$countstm, data_storage$countspm, min_intron_len=1000)
    } else { # if no counts
      data_storage$plot_obj <- plot_gene(data_storage$gene_to_plot, data_storage$res_tx_import, data_storage$res_pep_import, data_storage$res_ORF_import)
    }
    
    # print the plot
    output$plot <- renderPlot({
      suppressWarnings(print(data_storage$plot_obj))
    })
    
    # enable download button
    shinyjs::enable("vis_download_button")
    
  })
  
  # download handler for the plot
  output$vis_download_button <- downloadHandler(
    filename = function() {
      paste0(data_storage$gene_to_plot, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, h=10, w=20)
      print(data_storage$plot_obj)
      dev.off()
    }
  )
  
  # END VISUALISATION MODULE
  
  # remove session id tmp directory created each time app is run
  session$onSessionEnded(function() {
    if (dir.exists(session_id)) {
      unlink(session_id, recursive = TRUE)
    }
  })
  
}
