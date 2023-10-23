#!/usr/bin/env Rscript

library(optparse)
library(stringr)

source("load_parameters.R")
source("checkDirHierarchy.R")

args <- commandArgs(trailingOnly = TRUE)
# args <- "dea"

# Main help message
main_help <- "
Usage:
Rscript 00_pipeline_launcher.R <pipelineName> [--flag <flag_arg>]

Pipeline Order:
1) qc
2) process
3) filters
4) ctrl++
5) process (again)
6) dea
7) combine
8) process (last time)
9) dea (last time aswell)

Required Argument:
    <pipelineName>      This argument is required and must be one of the
                        following: qc, process, filters, ctrl++, dea, combine

Flag Argument:
    --flag <flag_arg>   There are optional arguments according to the required
                        argument

Example:
# For more details on each pipeline step
Rscript pipeline_launcher.R <pipelineName> --help

"

# make the call to the main help working:
# Rscript pipeline_launcher.R -h         # -h takes the place of args[1] -_-'
if (!is.na(args[1]) && grepl("^-", args[1])) {
  args[1] <- NA
}

# handle the multiple pipelineName input call main help:
# Rscript pipeline_launcher.R qc process
if (length(args) > 1 && !grepl("^-", args[1]) && !grepl("^-", args[2])) {
  cat(main_help)
  quit(status = 1)
}

# call main help message when user call unknown pipelineName:
if (!args[1] %in% c("NA", "qc", "process", "filters", "dea", "ctrl++", "combine")) {
  cat(main_help)
  quit(status = 1)
}

# defines the parameters to show in the displayed help menu
switch(args[1],
  "NA" = {
    cat(main_help)
    quit(status = 1)
  },
  "qc" = {
    option_list1 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = "Name of the dataset to process
				(required)."
        ),
        make_option(
        c("-m", "--mito_high"),
        action = "store",
        default = 10,
        type = "numeric",
        help = "Threshold for high mitochondrial gene expression,
			(default: 10).",
		metavar = "[1:100]"
      ),
      make_option(
          c("-n", "--mito_low"),
          action = "store",
          default = 1,
          type = "numeric",
          help = "Threshold for low mitochondrial gene expression
		  	(default: 1).",
		metavar = "[1:100]"
      ),
      make_option(
        c("-q", "--ribo_low"),
        action = "store",
        default = 25,
        type = "numeric",
        help = "Threshold for low ribosomal gene expression
			(default: 25).",
		metavar = "[1:100]"
      ),
      make_option(
        c("-f", "--min_feat"),
        action = "store",
        default = 200,
        type = "integer",
        help = "Minimum number of features detected in every cell
			(default: 200)."
      ),
      make_option(
        c("-c", "--min_cells"),
        action = "store",
        default = 3,
        type = "integer",
		help = "Minimum number in which a feature must be detected
			(default: 3)."
      ),
      make_option(
        c("--min_counts"),
        action = "store",
        default = 1500,
        type = "integer",
        help = "Minimum number of UMI counts per cell
			(default: 1500)."
      ),
      make_option(
        c("--max_counts"),
        action = "store",
        default = 150000,
        type = "integer",
        help = "Maximum number of UMI counts per cell
			(default: 150000)."
      )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog qc [--flag <flag_arg>]",
      option_list = option_list1,
      add_help_option = TRUE,
      description = "\nPerform quality control for single-cell RNA sequencing data.",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "process" = {
    option_list2 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = "Name of the dataset to process
				(required)."
        ),
        make_option(
            c("-f", "--filter"),
            action = "store",
            default = "complete",
            type = "character",
            help = "Dataset type, 'filtered' or 'complete'
				(default: 'complete')."
        ),
      make_option(
        c("-n", "--norm_method"),
        action = "store",
        default = "LogNormalize",
        type = "character",
        help = "Normalization method: 'NormalizeData' or 'SCTransform'
			(default: 'LogNormalize')."
      ),
      make_option(
        c("-v", "--hvg_method"),
        action = "store",
        default = "mvp",
        type = "character",
        help = "FindVariableFeatures method: 'vst', 'mvp', or 'disp'
			'vst' and 'disp' must be used with the option --hvg_number
			(default: 'mvp')."
      ),
      make_option(
        c("-w", "--hvg_number"),
        action = "store",
        default = FALSE,
        type = "integer",
        help = "Number of highly variable genes
			required with the option --hvg_method 'vst' or 'disp'
			(default: FALSE)."
      ),
      make_option(
        c("-d", "--do_scaling"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Compute scaling
			(default: FALSE)."
      ),
      make_option(
        c("-p", "--pca_npcs"),
        action = "store",
        default = 50,
        type = "integer",
        help = "Number of PCs to compute for RunPCA
			(default: 50)."
      ),
      make_option(
        c("-t", "--top_pcs"),
        action = "store",
        default = 30,
        type = "integer",
        help = "Number of PCs to select for downstream analysis
			(default: 30)."
      ),
      make_option(
        c("--pca_print"),
        action = "store",
        default = 20,
        type = "integer",
        help = "Number of features to print from the top of each PC
			(default: 20)."
      ),
      make_option(
        c("-s", "--selected_resolution"),
        action = "store",
        default = 1,
        type = "numeric",
        help = "Value of the resolution for FindClusters
			(default: 1)."
      ),
      make_option(
        c("-a", "--algo_clustering"),
        action = "store",
        default = 4,
        type = "integer",
        help = "Select the clustering algorithm for FindClusters
                1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement;
                3 = SLM algorithm; 4 = Leiden algorithm (recommended)
				(default: 4)."
      ),
      make_option(
          c("--good_quality"),
          action = "store_true",
          default = FALSE,
          type = "logical",
          help = "Process the dataset only with good quality cells, after the
		  ctrl++ pipeline
		  	(default: FALSE)."
      ),
          make_option(
              c("--combinedData"),
              action = "store_true",
              default = FALSE,
              type = "logical",
              help = "Process the dataset on combined datasets, after the
        combine pipeline
        (default: FALSE)"
      )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog process [--flag <flag_arg>]",
      option_list = option_list2,
      add_help_option = TRUE,
      description = "\nProcess and normalize single-cell RNA sequencing data.",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "filters" = {
      option_list3 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = "Name of the dataset to process
				(required)."
        ),
          make_option(
              c("--manual"),
              action = "store_true",
              default = FALSE,
              type = "logical",
              help = "Generate plots automatically if 'FALSE'. Otherwise, the
			  user needs to interact with 03_1_manualControls_<INPUT_DATASET>_pipeline.Rmd
			  	(default: FALSE)"
          ),
        make_option(
          c("-m", "--mito_thresholds"),
          action = "store",
          # default = 0,
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for mitochondrial expression, as a comma-separated list."
        ),
        make_option(
          c("-r", "--ribo_thresholds"),
          action = "store",
          # default = 0,
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for ribosomal expression, as a comma-separated list."
        ),
        make_option(
          c("-u", "--umi_thresholds"),
          action = "store",
          # default = 0,
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for read counts, as a comma-separated list."
        ),
        make_option(
          c("-f", "--feature_thresholds"),
          action = "store",
          # default = 0,
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for the number of features, as a comma-separated list."
        )
      )
      parsed <- OptionParser(
          usage = "Usage: \n\t%prog filters [--flag <flag_arg>]",
          option_list = option_list3,
          add_help_option = TRUE,
          description = "\nVisualize filtered cells on the UMAP plot of the complete dataset.",
          epilogue = "Add some details, examples, ...",
          formatter = IndentedHelpFormatter # TitleHelpFormatter
      )
  },
  "dea" = {
    option_list4 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = "Name of the dataset to process
				(required)."
        ),
      make_option(
        c("-k", "--markers_number"),
        action = "store",
        default = 20,
        type = "integer",
        help = "Reduce the list of markers to the top [-k] for each cluster at
		the selected resolution [-s]
			default: 20)."
      ),
      make_option(
          c("-S", "--scenario"),
          action = "store",
          default = 1,
          type = "integer",
          help = "Select the scenario to run:
                1 - no regression on cell cycle
                2 - global cell cycle regression, all phases are regressed
                3 - cycling cell cycle regression, G2M and S phases are regressed
				(default: 1)."
      ),
      make_option(
          c("-t", "--top_pcs"),
          action = "store",
          default = 30,
          type = "integer",
          help = "Number of PCs to select for downstream analysis
			(default: 30)"
      ),
      make_option(
          c("-f", "--filter"),
          action = "store",
          default = "filtered",
          type = "character",
          help = "Dataset type, 'filtered' or 'complete'
				(default: 'filtered')."
      ),
      make_option(
        c("-s", "--selected_resolution"),
        action = "store",
        default = 1,
        type = "numeric",
        help = "Value of the resolution for FindClusters
			(default: 1)."
      ),
      make_option(
        c("-a", "--algo_clustering"),
        action = "store",
        default = 4,
        type = "integer",
        help = "Select the clustering algorithm for FindClusters
                1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement;
                3 = SLM algorithm; 4 = Leiden algorithm (recommended)
				(default: 4)."
      ),
          make_option(
              c("--combinedData"),
              action = "store_true",
              default = FALSE,
              type = "logical",
              help = "Process the dataset on combined datasets, after the
        combine pipeline
        (default: FALSE)"
          )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog dea [--flag <flag_arg>]",
      option_list = option_list4,
      add_help_option = TRUE,
      description = "\nPerform differential expression analysis (DEA) on single-cell RNA sequencing data.",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "ctrl++" = {
    option_list5 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = "Name of the dataset to process
				(required)."
        ),
          make_option(
              c("-d", "--doublets_rate"),
              action = "store",
              default = 8,
              type = "integer",
              help = "Rate of doublets formation. Information available in the
			  sequencing kit manual
			  	(default: 8)."
          ),
          make_option(
            c("-p", "--pca_npcs"),
            action = "store",
            default = 50,
            type = "integer",
            help = "Number of PCs to compute for RunPCA
				(default: 50)."
          ),
          make_option(
              c("-t", "--top_pcs"),
              action = "store",
              default = 30,
              type = "integer",
              help = "Number of PCs to select for downstream analysis
			(default: 30)"
          ),
          make_option(
              c("-s", "--selected_resolution"),
              action = "store",
              default = 1,
              type = "numeric",
              help = "Value of the resolution for FindClusters
			  	(default: 1)."
          ),
          make_option(
              c("-a", "--algo_clustering"),
              action = "store",
              default = 4,
              type = "integer",
              help = "Select the clustering algorithm for FindClusters
                1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement;
                3 = SLM algorithm; 4 = Leiden algorithm (recommended)
				(default: 4)."
          ),
          make_option(
              c("-S", "--s_phase"),
              action = "store",
              # default = NA,
              type = "numeric",
              help = "Threshold to adjust cell cycle scoring for S phase
			  	(default: 0)."
          ),
          make_option(
              c("-G", "--g2m_phase"),
              action = "store",
              # default = NA,
              type = "numeric",
              help = "Threshold to adjust cell cycle scoring for G2M phase
			  	(default: 0)."
          ),
          make_option(
            c("-f", "--filter_features"),
            action = "store",
                        type = "character",
            help = "The list of features used to filter out unwanted cells. 
            For example, filtering out the HBA1 gene will remove all red blood cells, as it is a gene specifically expressed by such cells.
            The list should be in quotes and each feature must be separated by a colon."
          )
      )
      parsed <- OptionParser(
          usage = "Usage: \n\t%prog ctrl++ [--flag <flag_arg>]",
          option_list = option_list5,
          add_help_option = TRUE,
          description = "\nControl and adjust single-cell RNA sequencing data for doublets and cell cycle effects.",
          epilogue = "Add some details, examples, ...",
          formatter = IndentedHelpFormatter # TitleHelpFormatter
      )
  },
  "combine" = {
      option_list6 <- list(
          make_option(
              c("-I", "--input_list"),
              action = "store",
              default = NA,
              type = "character",
              help = "Comma separated list of input datasets
				(required)."
          ),
          make_option(
              c("-S", "--scenario"),
              action = "store",
              default = 1,
              type = "integer",
              help = "Select the scenario to run:
                1 - no regression on cell cycle
                2 - global cell cycle regression, all phases are regressed
                3 - cycling cell cycle regression, G2M and S phases are regressed
				(default: 1)."
          ),
          make_option(
              c("-c", "--combineMethod"),
              action = "store",
              default = NA,
              type = "character",
              help = "Select how to combine the datasets:
              - 'merge' will simply merge the datasets
              - 'blkS' performs the block-wise CCA integration from Seurat
              - 'seqS' performs the sequential CCA integration from Seurat
            (required)."
          )
      )
      parsed <- OptionParser(
          usage = "Usage: \n\t%prog combine [--flag <flag_arg>]",
          option_list = option_list6,
          add_help_option = TRUE,
          description = "\nCombine multiple datasets. Choose between merge or one
          of the two integration mode, blkS or seqS",
          epilogue = "Add some details, examples, ...",
          formatter = IndentedHelpFormatter # TitleHelpFormatter
      )
  }
)

opt <- parse_args(parsed, positional_arguments = TRUE)

# opt$options$input_dataset <- "cardioKO"
#opt$options$input_list <- "highDiet,normalDiet"
# opt$options$filter <- "filtered"
# opt$options$good_quality <- TRUE
PATH_REQUIREMENTS <- "../01_requirements/"
if (is.null(opt$options$input_list)) {
    load_parameters(paste0(PATH_REQUIREMENTS, "globalParameters_", opt$options$input_dataset,".param"))
} else {
	# 1) Split input_list based on comma(s)
	input_datasets <- strsplit(opt$options$input_list, ",")[[1]]
	
	# 2) Get the list of files present in the PATH_REQUIREMENTS directory
	file_list <- list.files(PATH_REQUIREMENTS)
	
	# 3) Match multiple patterns using grepl and &
	matching_element <- Reduce(`&`, lapply(input_datasets, grepl, file_list))
	
	# 4) Get the file that contains all the names in input_datasets
	matching_file <- file_list[matching_element]
    
	load_parameters(paste0(PATH_REQUIREMENTS, matching_file))
	# DATASET <- paste0(input_datasets, collapse = "_") # don't know if I need to keep it ???
}

if (TRUE){
    general_seed = as.numeric(general_seed)
    # qc pipeline variables ----------------------
    mito_high = if (exists("MITO_HIGH")) as.numeric(MITO_HIGH) else opt$options$mito_high
    mito_low = if (exists("MITO_LOW")) as.numeric(MITO_LOW) else opt$options$mito_low
    ribo_low = if (exists("RIBO_LOW")) as.numeric(RIBO_LOW) else opt$options$ribo_low
    min_feat = if (exists("MIN_FEAT")) as.numeric(MIN_FEAT) else opt$options$min_feat
    min_cells = if (exists("MIN_CELLS")) as.numeric(MIN_CELLS) else opt$options$min_cells
    min_counts = if (exists("MIN_COUNTS")) as.numeric(MIN_COUNTS) else opt$options$min_counts
    max_counts = if (exists("MAX_COUNTS")) as.numeric(MAX_COUNTS) else opt$options$max_counts
    # process pipeline variables -----------------
    FILTER = if (exists("FILTER")) FILTER else opt$options$filter
    norm_meth = if (exists("NORM_METH")) NORM_METH else opt$options$norm_method
    hvg_meth = if (exists("HVG_METH")) HVG_METH else opt$options$hvg_method
    hvg_num = if (exists("HVG_NUM")) as.integer(HVG_NUM) else opt$options$hvg_number
    do_scale = if (exists("DO_SCALE")) as.logical(DO_SCALE) else opt$options$do_scaling
    pca_npcs = if (exists("PCA_NPCS")) as.integer(PCA_NPCS) else opt$options$pca_npcs
    pca_print = if (exists("PCA_PRINT")) as.integer(PCA_PRINT) else opt$options$pca_print
    top_pcs = if (exists("TOP_PCS")) as.integer(TOP_PCS) else opt$options$top_pcs
    clust_res = if (exists("CLUST_RES")) as.numeric(CLUST_RES) else opt$options$selected_resolution
    clust_meth = if (exists("CLUST_METH")) as.integer(CLUST_METH) else opt$options$algo_clustering
    # advanced filter pipeline variables  --------
    mito_thresholds = if (!is.null(if (exists("MITO_THRESHOLDS")) MITO_THRESHOLDS else opt$options$mito_thresholds)) as.numeric(unlist(strsplit(if (exists("MITO_THRESHOLDS")) MITO_THRESHOLDS else opt$options$mito_thresholds, ",")))
    ribo_thresholds = if (!is.null(if (exists("RIBO_THRESHOLDS")) RIBO_THRESHOLDS else opt$options$ribo_thresholds)) as.numeric(unlist(strsplit(if (exists("RIBO_THRESHOLDS")) RIBO_THRESHOLDS else opt$options$ribo_thresholds, ",")))
    umi_thresholds = if (!is.null(if (exists("UMI_THRESHOLDS")) UMI_THRESHOLDS else opt$options$umi_thresholds)) as.numeric(unlist(strsplit(if (exists("UMI_THRESHOLDS")) UMI_THRESHOLDS else opt$options$umi_thresholds, ",")))
    feature_thresholds = if (!is.null(if (exists("FEATURE_THRESHOLDS")) FEATURE_THRESHOLDS else opt$options$feature_thresholds)) as.numeric(unlist(strsplit(if (exists("FEATURE_THRESHOLDS")) FEATURE_THRESHOLDS else opt$options$feature_thresholds, ",")))
    # deg pipeline variables ---------------------
    top_markers = opt$options$markers_number
    genes_of_interest = if (exists("GENES_OF_INTEREST")) GENES_OF_INTEREST else NULL
    # doublets removal ---------------------------
    doublets_rate = if (exists("DOUBLETS_RATE")) as.numeric(DOUBLETS_RATE) else opt$options$doublets_rate
    # cell cycle regression ----------------------
    s_thresh = if (exists("S_PHASE")) as.numeric(S_PHASE) else opt$options$s_phase
    g2m_thresh = if (exists("G2M_PHASE")) as.numeric(G2M_PHASE) else opt$options$g2m_phase
    scenario = if (exists("REGRESSION_SCENARIO")) REGRESSION_SCENARIO else opt$options$scenario         # Maybe define it into globalParams => param to be eval in step ctrl++ if set to 1 (ie not interested in CC regression)
    # combining multiple datasets
    combine_meth = if (exists("COMB_METH")) COMB_METH else opt$options$combineMethod
    # unclassified variables ---------------------
    manual = opt$options$manual
    combinedD = if(!is.null(opt$options$combinedData)) opt$options$combinedData else FALSE
    goodQ = opt$options$good_quality
    gn_col = if (exists("GENE_NAME_COLUMN")) as.numeric(GENE_NAME_COLUMN) else opt$options$gene_name_col

    ff_list = if (!is.null(if (exists("FILTER_FEATURES")) FILTER_FEATURES else opt$options$filter_features)) strsplit(if (exists("FILTER_FEATURES")) FILTER_FEATURES else opt$options$filter_features, ",")
}

# combinedD <- TRUE

if (combinedD) {
	FILTER = "filtered"
}


checkDirHierarchy()

switch(args[1],
       "qc" = {
           rmarkdown::render(
               "01_qc_pipeline.Rmd",
               output_file = paste0(PATH_OUT_HTML, "01_qc_", DATASET, "_", Sys.Date(), ".html")
           )
       },
       "process" = {
           if (goodQ) {
               rmarkdown::render(
                   "02_process_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "05_process_", DATASET, "_goodQualityCells_", Sys.Date(), ".html")
               )
           } else if (combinedD) {
               rmarkdown::render(
                   "02_process_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "08_process_", DATASET, "_combinedData_scenario_", opt$options$scenario, Sys.Date(), ".html")
               )
           } else {
               rmarkdown::render(
                   "02_process_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "02_process_", DATASET, "_", FILTER, "_", Sys.Date(), ".html")
               )
           }
       },
       "filters" = {
           if (opt$options$manual) {
               rmarkdown::render(
                   "03_1_manualControls_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "03_1_manualControls_", DATASET, "_", Sys.Date(), ".html")
               )
           } else {
               rmarkdown::render(
                   "03_filtersControls_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "03_filtersControls_", DATASET, "_", Sys.Date(), ".html")
               )   
           }
       },
       "dea" = {
           if (FILTER == "filtered") {
               rmarkdown::render(
                   "05_annotDEAviz_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "06_annotDEAviz_", DATASET, "_scenario_", opt$options$scenario, "_", Sys.Date(), ".html")
               )
           } else if (FILTER == "complete") {
               rmarkdown::render(
                   "04_1_deg_pipeline.Rmd",
                   output_file = paste0(PATH_OUT_HTML, "04_1_deg_", DATASET, "_", Sys.Date(), ".html")
               )   
           }
       },
       "ctrl++" = {
           rmarkdown::render(
               "04_additionalControls.Rmd",
               output_file = paste0(PATH_OUT_HTML, "04_additionalControls_", DATASET, "_", Sys.Date(), ".html")
           )
       },
       "combine" = {
           rmarkdown::render(
               "06_combineData_pipeline.Rmd",
               output_file = paste0(PATH_OUT_HTML, "07_combineData_", DATASET, "_", Sys.Date(), ".html")
           )
       }
)

