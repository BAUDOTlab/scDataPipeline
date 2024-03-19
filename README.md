# scDataPipeline

This is a standard single-cell RNA sequencing data analysis pipeline in R and Python. The data must be first preprocessed with tools like the [CellRanger](https://github.com/10XGenomics/cellranger) pipeline from 10X Genomics to be in the form of count matrices, with the corresponding features and barcodes tsv files. The pipeline uses mostly the [Seurat version 4](https://satijalab.org/seurat/) R package thoughout the analysis, and the [scanpy](https://scanpy.readthedocs.io/en/stable/) python package for the differential composition analysis.

At each step, the state of the seurat object is saved in a RDS file and all required metadatas are stored in the analysis directory, allowing the analysis to be started from any already computed step. Each stage also produces a report in the form of an HTML file with all the relevant visualisations and information. Every figure contained in the reports is saved separately, both in PDF and PNG formats.

## Performed analyses

- Quality Control.
- Bad quality cells removal according to QC thresholds, doublet removal through [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and undesired cells flag, such as Red Blood Cells.
- Downstream analysis along 3 scenarios : no regression, global cell cycle regression, and cycling cell regression.
- Differential Expression Analysis between clusters and conditions.
- Dataset combination (via merging or integration).
- Differential Composition Analysis.
- Differential Gene Set Expression.

## Getting Started

### Requirements

- [`docker`](https://docker.com)

Install docker and the optional but useful docker compose plugin with:

```bash
sudo apt install docker.io docker-plugin-compose
```

### Installation

#### The pipeline

To set up a new project using scDataPipeline, follow this process (you must repeat it for each new project):

1. Create a folder for your project

```bash
mkdir [project_name]
cd [project_name]
```

2. Clone this repository inside the project's folder.

```bash
git clone git@github.com:BAUDOTlab/scDataPipeline.git
```

3. Create a subfolder dedicated to store the raw data, which we suggest you name `00_rawdata/`. It is where you should place your data: CellRanger output, atlas if you use one, etc... Then create a subfolder next to the previous one which **MUST** be named `01_requirements`. This is where the parameter files should be moved to and filled.

```bash
mkdir 00_rawdata 01_requirements
```

4. Copy the configuration file models into the `01_requirements` folder

```bash
cp scDataPipeline/globalParameters* 01_requirements/
```

All set! You can now put your data in the subfolder you just created and fill the configuration files (see [below](#filling-the-parameter-files)).

#### The Docker image

This pipeline runs in a docker container. The base image contains R, python and all the packages required to run the pipeline, as well as rstudio server and desktop for debugging purpose. Build the image using the following command (a pullable image might eventually be available):

```bash
docker build -t scdatapipeline .
```

Warning: the process take a little over **one hour** to complete (most of which is R package install). the `-t` option is used to name your image. You can use any name you want, as long as it is in lowercase. For more information, see the [docker build documentation](https://docs.docker.com/reference/cli/docker/image/build).

### Filling the parameter files

There are two types of parameter files: **globalParameters_MODEL_single.toml** for single datasets and **globalParameters_MODEL_combined.toml** for combined datasets. Note that in both case, the file must have the exact value of `DATASET` in its name. For example, using a dataset named `My_Dataset`, the name of the file should be `globalParameters_My_Dataset_single.toml` and the first lines should look like this (be careful of the case):

```toml
# Welcome to the config file relative to a dataset

# WARNING: the title must contain the name of the DATASET as defined in the
# file.


# ===================
# required parameters
# ===================
DATASET = "My_Dataset"
```

##### Regarding the paths

`PATH_ROOT` must contain the absolute path towards the root folder of our project. The other paths included in `path.input` and `path.output` can be either absolute paths or paths relative to `PATH_ROOT`, with the exception of `PATH_ATLAS_FILE`, which is relative to `PATH_ATLAS`.

```toml
[path]
PATH_ROOT = "XXXXX"                     # base path from which the others depends

[path.input]
PATH_INPUT_LABDATA = "XXXXX"            # path pointing to the input data (typically $PATH_ROOT/00_rawData)
PATH_ATLAS = "XXXXX"                    # path to the atlas folder
PATH_ATLAS_FILE = "XXXXX"               # /!\ if it is not absolute, it is relative to $PATH_ATLAS
# optional params
PATH_GENES_OF_INTEREST = "XXXXX"		# path to the list of specific genes of interest (text file with one gene per line), usually in $PATH_ROOT/01_requirements.

[path.output]
PATH_RDS_OBJECTS	= "50_rdsObjects"	# create path $PATH_ROOT/PATH_RDS_OBJECTS
PATH_OUT_HTML		= "51_htmlFiles"	# create path $PATH_ROOT/PATH_OUT_HTML
PATH_OUT_FIG		= "52_figures"		# create path $PATH_ROOT/PATH_OUT_FIG/DATASET
```

## How to use

We use Docker as a virtual environment. To do that, we mount a local folder to the `home` of the container. Note that by default, anyone in a docker container has **root privileges** on the host. If that is not what you want, consider running [Rootless Docker](https://docs.docker.com/engine/security/rootless/).

### Vanilla Docker

Run the following command inside the scDataPipeline folder:

```bash
docker run -it --rm --mount type=bind,source=/path/to/your/project,target=/home/[PROJECT_NAME] -e SCDP_PATH_ROOT="/home/[PROJECT_NAME]" scdatapipeline
```

This will open you an interactive bash shell in your container. The `SCDP_PATH_ROOT` is an environment variable used to tell the pipeline that it is inside a container and should not be set on the host (or have the same value as the `PATH_ROOT` field of the configuration file). Note that when you exit the container, any content it has inside (additionnal packages installed, files created outside of `/home`...) will be deleted. If you want to be able to reopen a container at a later date, remove the `--rm` option.

### Docker Compose

Docker compose offer a more complete and easier configuration of the container than vanilla docker. To use it, edit the [docker compose](docker-compose.yml) file. Then, you can start the container with this simple command:

```bash
docker compose up pipeline
```

And to stop and delete the container:

```bash
docker compose down pipeline
```

(Replace `down` with `stop` if you only with to pause the container without removing it altogether.)

### Running the Pipeline 

Once you are in the container, there is only one script you need to launch each time. Go to the `scDataPipeline` folder and run:

```
Rscript 00_pipeline_launcher.R -h
```

To see a list of the different steps, then to see the arguments of the step you want run, for example `qc`:

```
Rscript 00_pipeline_launcher.R qc -h
```

When running long steps such as `ctrl` or `dea`, it might be a good idea to use `nohup` and `&` to be able to still use the shell while the program is running like so:

```
nohup Rscript 00_pipeline_launcher.R ctrl -i [dataset_name] &
```

And then you can access the output of the script with

```
tail -f nohup.out
```

### Example workflow

If you want to stick to the pipeline's defaults or have already written all your parameters in the parameter files, you can run the following lines for the complete workflow. Note that you should run all the steps one after the other, and always check the HTML report before going on. Some steps, like `qc` and `filters --manual` are designed to be rerunned multiple times to adjust the thresholds before the effective filtering of the data.

```bash
# Run the quality control on your dataset. It does not remove any cell but allow to visualize the thresholds.
Rscript 00_pipeline_launcher.R qc -i YOUR_DATASET

# Run the initial preprocessing step, both on the complete data and the filtered dataset without cells that fail the QC.
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter complete
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter filtered

# Visualize the filtered cells on the complete UMAP and allow the evaluation of the thresholds
Rscript 00_pipeline_launcher.R filters -i YOUR_DATASET
Rscript 00_pipeline_launcher.R filters -i YOUR_DATASET --manual               # Optional step

# Run doublets removal and cell cycle regression, as well as flagging of cell expressing specific genes.
# Outputs 3 RDS files, one for each scenario : no CC regression, full regression and regression of cycling cells only.
Rscript 00_pipeline_launcher.R ctrl -i YOUR_DATASET

# Run the second preprocessing on the cleaned dataset
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter filtered --good_quality

# Perform pseudobulk differential expression analysis between clusters, using the "each versus all" method.
Rscript 00_pipeline_launcher.R dea -i YOUR_DATASET

# Run the integration of single datasets into one so that they can be easily compared
Rscript 00_pipeline_launcher.R combine -I YOUR_DATASET,YOUR_SECOND_DATASET

# Run the last preprocessing on the combined dataset
Rscript 00_pipeline_launcher.R process -i YOUR_COMBINED_DATASET

# Perform the same differential expression analysis as above, on the combined dataset
Rscript 00_pipeline_launcher.R dea -i YOUR_COMBINED_DATASET

# Run the differential composition analysis 
Rscript 00_pipeline_launcher.R da -i YOUR_COMBINED_DATASET

# Run the differential expression analysis between groups (e.g. conditions) of each clusters,
# then run gene set enrichment analysis to highlight pathways
Rscript 00_pipeline_launcher.R deg -i YOUR_COMBINED_DATASET
```

## Getting Help

Need help? Post your question on the [issue board](https://github.com/BAUDOTlab/scDataPipeline/issues) and we'll do our best to answer.

## Contributing

There are currently no clear guidelines for contributing, but anyone is free to fork this repository and create a pull request.

## License

The scDataPipeline is distributed under the **MIT license**. See [LICENSE](./LICENSE) for more details.
