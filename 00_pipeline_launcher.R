#!/usr/bin/env Rscript

library(optparse)
library(stringr)
library(RcppTOML)

source("load_parameters.R")
source("checkDirHierarchy.R")
source("data_management.R")

args <- commandArgs(trailingOnly = TRUE)

#args <- ""

# Main help message
main_help <- "
Usage:
Rscript 00_pipeline_launcher.R <pipelineName> [--flag <flag_arg>]

Pipeline Order:
1) qc
2) process
3) filters
4) ctrl
5) process (again)
6) dea
7) combine
8) process (last time)
9) dea (last time aswell)
10) da
11) deg

Required Argument:
    <pipelineName>      This argument is required and must be one of the
                        following: qc, process, filters, ctrl, dea, combine, da,
                        deg

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
if (!args[1] %in% c("NA", "qc", "process", "filters", "dea", "ctrl", "combine", "da", "deg")) {
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
		  ctrl pipeline
		  	(default: FALSE)."
      ),
      make_option(
        c("--rm_clust"),
        action = "store",
        default = NA,
        type = "character",
        help = "Clusters to remove before the new clustering. All cells in the
        clusters will be removed. This information should be based on the
        observed clusters of the ctrl step.
        This parameter is only active with --good_quality, and should be a
        comma-separated list of cluster numbers."
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
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for mitochondrial expression, as a comma-separated list."
        ),
        make_option(
          c("-r", "--ribo_thresholds"),
          action = "store",
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for ribosomal expression, as a comma-separated list."
        ),
        make_option(
          c("-u", "--umi_thresholds"),
          action = "store",
          type = "character",
          help = "ONLY USED WITH --manual
          The list of the thresholds for read counts, as a comma-separated list."
        ),
        make_option(
          c("-f", "--feature_thresholds"),
          action = "store",
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
          help = "Select the scenario(s) to run:
          1 - no regression on cell cycle
          2 - global cell cycle regression, all phases are regressed
          3 - cycling cell cycle regression, G2M and S phases are regressed
      (default: 1).
      If you select multiple scenarios (separeted with a comma), they will run sequentially."
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
  "ctrl" = {
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
              type = "numeric",
              help = "Threshold to adjust cell cycle scoring for S phase
			  	(default: 0)."
          ),
          make_option(
              c("-G", "--g2m_phase"),
              action = "store",
              type = "numeric",
              help = "Threshold to adjust cell cycle scoring for G2M phase
			  	(default: 0)."
          )
      )
      parsed <- OptionParser(
          usage = "Usage: \n\t%prog ctrl [--flag <flag_arg>]",
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
              help = "Select the scenario(s) to run:
          1 - no regression on cell cycle
          2 - global cell cycle regression, all phases are regressed
          3 - cycling cell cycle regression, G2M and S phases are regressed
      (default: 1).
      If you select multiple scenarios (separeted with a comma), they will run sequentially."
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
  },
  "da" = {
    option_list7 <- list(
      make_option(
        c("-i", "--input_dataset"),
        action = "store",
        default = NA,
        type = "character",
        help = "Name of the dataset to process
          (required)."
      ),
      make_option(
        c("-m", "--daMethod"),
        action = "store",
        default = "all",
        type = "character",
        help = "Comma separated list of the methods to use for differential abundance analysis:
        - 'meld' performs the MELD method
            (required)."
      ),
      make_option(
        c("-S", "--scenario"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Select the scenario(s) to run:
          1 - no regression on cell cycle
          2 - global cell cycle regression, all phases are regressed
          3 - cycling cell cycle regression, G2M and S phases are regressed
      (default: 1).
      If you select multiple scenarios (separeted with a comma), they will run sequentially."
        )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog da [--flag <flag_arg>]",
      option_list = option_list7,
      add_help_option = TRUE,
      description = "\nPerform differential abundance analysis on the dataset.",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "deg" = {
    option_list8 = list(
        make_option(
        c("-i", "--input_dataset"),
        action = "store",
        default = NA,
        type = "character",
        help = "Name of the dataset to process
          (required)."
      ),
        make_option(
        c("-s", "--scenario"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Select the scenario(s) to run:
          1 - no regression on cell cycle
          2 - global cell cycle regression, all phases are regressed
          3 - cycling cell cycle regression, G2M and S phases are regressed
      (default: 1).
      If you select multiple scenarios (separeted with a comma), they will run sequentially."
        )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog deg [--flag <flag_arg>]",
      option_list = option_list8,
      add_help_option = TRUE,
      description = "\nPerform differential abundance analysis on the dataset.",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  }
)

opt <- parse_args(parsed, positional_arguments = TRUE)

# opt$options$input_dataset <- ""
# opt$options$input_list <- ""
# opt$options$filter <- "filtered"
# opt$options$good_quality <- TRUE
# opt$options$rm_clust <- ""
PATH_REQUIREMENTS <- "../01_requirements/"
if (is.null(opt$options$input_list)) {
    load_parameters(paste0(PATH_REQUIREMENTS, "globalParameters_", opt$options$input_dataset,".toml"))
} else {
	# 1) Split input_list based on comma(s)
	input_datasets <- strsplit(opt$options$input_list, ",")[[1]]
	
	# 2) Get the list of files present in the PATH_REQUIREMENTS directory
	file_list <- list.files(PATH_REQUIREMENTS)
	
	# 3) Match multiple patterns using grepl and &
	matching_element <- Reduce(`&`, lapply(c(input_datasets, "toml"), grepl, file_list))
	
	# 4) Get the file that contains all the names in input_datasets
	matching_file <- file_list[matching_element]
    
	load_parameters(paste0(PATH_REQUIREMENTS, matching_file))

  # 5) Merge the feature.tsv files
  merge_features(input_datasets)
	# DATASET <- paste0(input_datasets, collapse = "_") # don't know if I need to keep it ???
}

# Check for spaces after DATASET or inside DATASET
if (grepl("\\s+", DATASET)) {
    stop("Error: DATASET contains spaces. Please ensure there are no spaces in the DATASET.")
}


if(!file.exists(PATH_ATLAS_FILE) && args[1] %in% c("dea", "ctrl", "da")){
  print(paste("The provided path of the atlas file:", PATH_ATLAS_FILE, "is not valid. Some parts of the report will not be generated."))
  atlas <- FALSE
} else {
  atlas <- TRUE
}

if(!exists("PATH_MANUAL_ANNOTATION") && args[1] %in% c("da")){
    print(paste("The provided path for the manual annotation table:", PATH_MANUAL_ANNOTATION, "is not valid. Some parts of the report will not be generated."))
    manAnn <- FALSE
} else {
    manAnn <- TRUE
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
    mito_thresholds = if (!is.null(if (exists("MITO_THRESHOLDS")) MITO_THRESHOLDS else opt$options$mito_thresholds)) if (exists("MITO_THRESHOLDS")) MITO_THRESHOLDS else as.numeric(unlist(strsplit(opt$options$mito_thresholds)))
    ribo_thresholds = if (!is.null(if (exists("RIBO_THRESHOLDS")) RIBO_THRESHOLDS else opt$options$ribo_thresholds)) if (exists("RIBO_THRESHOLDS")) RIBO_THRESHOLDS else as.numeric(unlist(strsplit(opt$options$ribo_thresholds, ",")))
    umi_thresholds = if (!is.null(if (exists("UMI_THRESHOLDS")) UMI_THRESHOLDS else opt$options$umi_thresholds)) if (exists("UMI_THRESHOLDS")) UMI_THRESHOLDS else as.integer(unlist(strsplit(opt$options$umi_thresholds, ",")))
    feature_thresholds = if (!is.null(if (exists("FEATURE_THRESHOLDS")) FEATURE_THRESHOLDS else opt$options$feature_thresholds)) if (exists("FEATURE_THRESHOLDS")) FEATURE_THRESHOLDS else as.integer(unlist(strsplit(opt$options$feature_thresholds, ",")))
    # deg pipeline variables ---------------------
    top_markers = opt$options$markers_number
    genes_of_interest = if (exists("PATH_GENES_OF_INTEREST")) PATH_GENES_OF_INTEREST else NULL
    # doublets removal ---------------------------
    doublets_rate = if (exists("DOUBLETS_RATE")) as.numeric(DOUBLETS_RATE) else opt$options$doublets_rate
    # cell cycle regression ----------------------
    s_thresh = if (exists("S_PHASE")) as.numeric(S_PHASE) else opt$options$s_phase
    g2m_thresh = if (exists("G2M_PHASE")) as.numeric(G2M_PHASE) else opt$options$g2m_phase
    scenarios = if (exists("REGRESSION_SCENARIO")) REGRESSION_SCENARIO else opt$options$scenario         # Maybe define it into globalParams => param to be eval in step ctrl if set to 1 (ie not interested in CC regression)
    # combining multiple datasets
    combine_meth = if (exists("COMB_METH")) COMB_METH else opt$options$combineMethod
    # da variables --------------------------------
    da_meth = if (exists("DA_METH")) DA_METH else opt$options$daMethod
    if(!is.null(da_meth) && da_meth == "all") {
      da_meth = list("meld")
    }
    # unclassified variables ---------------------
    manual = opt$options$manual
    combinedD = if (exists("COMBINED")) COMBINED else if(!is.null(opt$options$combinedData)) opt$options$combinedData else FALSE
    goodQ = opt$options$good_quality
    gn_col = if (exists("GENE_NAME_COLUMN")) as.numeric(GENE_NAME_COLUMN) else opt$options$gene_name_col
    obs_list = if (!is.null(if (exists("OBSERVE_FEATURES")) OBSERVE_FEATURES else opt$options$observe_features)) if (exists("OBSERVE_FEATURES")) OBSERVE_FEATURES else unlist(strsplit(opt$options$observe_features, ","))
    rm_clust = if (!is.null(goodQ) && goodQ) unlist(strsplit(opt$options$rm_clust, ","))
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
               output_file = file.path(PATH_OUT_HTML, paste0("01_qc_", DATASET, "_", Sys.Date(), ".html"))
           )
       },
       "process" = {
           if (goodQ) {
              lapply(scenarios, function(scenario){
               rmarkdown::render(
                   "02_process_pipeline.Rmd",
                   output_file = file.path(PATH_OUT_HTML, paste0("05_process_", DATASET, "_goodQualityCells_", Sys.Date(), ".html"))
               )
              })
           } else if (combinedD) {
              lapply(scenarios, function(scenario){
                print(paste0("Processing scenario ", scenario, " on combined data."))
                rmarkdown::render(
                    "02_process_pipeline.Rmd",
                    output_file = file.path(PATH_OUT_HTML, paste0("08_process_", DATASET, "_combinedData_scenario_", scenario, "_", Sys.Date(), ".html"))
                )
               })
               
           } else {
              scenario = "1" # rubberband fix
              rmarkdown::render(
                "02_process_pipeline.Rmd",
                output_file = file.path(PATH_OUT_HTML, paste0("02_process_", DATASET, "_", FILTER, "_", Sys.Date(), ".html"))
               )
           }
       },
       "filters" = {
           if (opt$options$manual) {
               rmarkdown::render(
                   "03_1_manualControls_pipeline.Rmd",
                   output_file = file.path(PATH_OUT_HTML, paste0("03_1_manualControls_", DATASET, "_", Sys.Date(), ".html"))
               )
           } else {
               rmarkdown::render(
                   "03_filtersControls_pipeline.Rmd",
                   output_file = file.path(PATH_OUT_HTML, paste0("03_filtersControls_", DATASET, "_", Sys.Date(), ".html"))
               )   
           }
       },
       "dea" = {
           if (combinedD) {
              lapply(scenarios, function(scenario){
                print(paste0("Performing DEA analysis on scenario ", scenario))
                rmarkdown::render(
                    "05_annotDEAviz_pipeline.Rmd",
                    output_file = file.path(PATH_OUT_HTML, paste0("09_annotDEAviz_", DATASET, "_scenario_", scenario, "_", Sys.Date(), ".html"))
                )
              })
           } else {
              lapply(scenarios, function(scenario){
                print(paste0("Performing DEA analysis on scenario ", scenario))
                rmarkdown::render(
                    "05_annotDEAviz_pipeline.Rmd",
                    output_file = file.path(PATH_OUT_HTML, paste0("06_annotDEAviz_", DATASET, "_scenario_", scenario, "_", Sys.Date(), ".html"))
                )
              })
           }
       },
       "ctrl" = {
             print(paste0("Performing additionnal controls for scenarios ", scenarios))
             rmarkdown::render(
             "04_additionalControls.Rmd",
               output_file = file.path(PATH_OUT_HTML, paste0("04_additionalControls_", DATASET, "_scenarios_", paste(scenarios, sep="-"), "_", Sys.Date(), ".html"))
           )
       },
       "combine" = {
          lapply(scenarios, function(scenario){
            print(paste0("Combining datasets ", toString(input_datasets), " for scenario ", scenario))
            rmarkdown::render(
                "06_combineData_pipeline.Rmd",
                output_file = file.path(PATH_OUT_HTML, paste0("07_combineData_", DATASET, "_scenario_", scenario, "_",  Sys.Date(), ".html"))
            )
          })
       },
       "da" = {
        lapply(scenarios, function(scenario){
            lapply(da_meth, function(meth){
              print(paste0("Performing DA analysis using ", meth, " method on ", DATASET, " dataset for scenario ", scenario))
              if(file.exists(paste0("07_",meth,"_pipeline.Rmd"))){
                rmarkdown::render(
                  paste0("07_",meth,"_pipeline.Rmd"),
                  output_file = file.path(PATH_OUT_HTML, paste0("10_DAanalysis_", DATASET, "_meth_", meth, "_scenario_", scenario, "_", Sys.Date(), ".html"))
                )
              }
        })
        })
       },
       "deg" = {
        lapply(scenarios, function(scenario){
              print(paste0("Performing DEG analysis on ", DATASET, " dataset for scenario ", scenario))
              if(file.exists(paste0("08_DEG_analysis.Rmd"))){
                rmarkdown::render(
                  paste0("08_DEG_analysis.Rmd"),
                  output_file = file.path(PATH_OUT_HTML, paste0("11_DEGanalysis_", DATASET, "_scenario_", scenario, "_", Sys.Date(), ".html"))
                )
              }
        })
       }
)