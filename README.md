# scDataPipeline

> [!CAUTION]
> This tool is still a work in progress under active development. It is not complete as of now and may prove tedious to set up.

This is a standard single-cell RNA sequencing data analysis pipeline in R and Python. The data must be first preprocessed with tools like [CellRanger](https://github.com/10XGenomics/cellranger) to be in the form of count matrices, with the corresponding features and barcodes tsv files. The pipeline uses the [Seurat version 4](https://satijalab.org/seurat/) package thoughout the analysis.

At each step, the state of the seurat object is saved in a RDS file and all required metadatas are stored in the analysis directory, allowing the analysis to be started from any already computed step. Each stage also produces a report in the form of an HTML file with all the relevant visualisations and information.

## Performed analyses

- Quality Control.
- Filtering according to QC thresholds, [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and undesired cells, such as Red Blood Cells.
- Downstream analysis along 3 scenarios : no regression, global cell cycle regression, and cycling cell regression.
- Differential Expression Analysis.
- Dataset combination (via merging or integration).
- Differential Composition Analysis.
- **[TBD]** Differential Gene Set Expression.

## Getting Started

### Requirements

- [`docker`](https://docker.com)

Install docker with:

```bash
sudo apt install docker.io
```

### Installation

## The pipeline

To set up a new project using scDataPipeline, follow this process (you must repeat it for each new project):

1. Create a folder for your project
2. Clone this repository inside the project's folder.
3. Create a subfolder destined to store the raw data, which we suggest you name `00_rawdata/`. It is where to place data such as the CellRanger output or your cell atlas, either directly or through a symbolic link.
4. Create a subfolder next to the previous one which **MUST** be named `01_requirements`. This is where the parameter files should be moved to and filled.

```bash
mkdir [project_name]
cd [project_name]
git clone git@github.com:BAUDOTlab/scDataPipeline.git
mkdir 00_rawdata 01_requirements
cp scDataPipeline/globalParameters* 01_requirements/
# Fill the 00_rawdata folder with your data
```

## The Docker image

This pipeline runs in a docker container. The base image contains R, python and all the packages required to run the pipeline, as well as rstudio server and desktop for debugging purpose. Since using Rstudio desktop from inside a container requires root privileges, you may not need/be interested in it. Build the image using the following command (a pullable image might eventually be available):

```bash
nohup docker build -t scdatapipeline . &
```

The `nohup` and `&` are not part of the docker build command, but they avoid blocking the shell since the process take a little over **one hour** to complete (most of which is R package install). the `-t` option is used to name your image. You can use any name you want, as long as it is in lowercase.

### Filling the parameter files

There is two types of parameter files: **globalParameters_MODEL_single.toml** for single datasets and **globalParameters_MODEL_combined.toml** for combined datasets. Note that in both case, the file must have the exact value of `DATASET` in its name.

#### Regarding the paths

`PATH_ROOT` must contain the absolute path towards the root folder of our project. The other paths included in `path.input` and `path.output` **MUST BE** paths relative to `PATH_ROOT`, with the exception of `PATH_ATLAS_FILE`, which is relative to `PATH_ATLAS`.

## How to use

You may have noticed that the Dockerfile don't have a `COPY` command. Indeed, as we are working with omics data that can take up quite some space, we are not working directly inside the container, but on the host. It's like using docker as a virtual environment. To do that, we mount a local folder to the `home` of the container. Note that by default, anyone in a docker container has **root privileges** on the host. If that is not what you want, consider running [Rootless Docker](https://docs.docker.com/engine/security/rootless/).

```bash
docker run -it --rm --mount type=bind,source=/path/to/your/project,target=/home -e SCDP_PATH_ROOT="/home/" scdatapipeline
```

This will open you an interactive bash shell in your container. The `SCDP_PATH_ROOT` is an environment variable used to tell the pipeline that it is inside a container and should not be set on the host (or have the same value as the `PATH_ROOT` field of the configuration file). Note that when you exit the container, any content it has inside (additionnal packages installed, files created outside of `/home`...) will be deleted. If you want to be able to reopen a container at a later date, remove the `--rm` option.

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
Rscript 00_pipeline_launcher.R qc -i YOUR_DATASET
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter complete
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter filtered
Rscript 00_pipeline_launcher.R filters -i YOUR_DATASET
Rscript 00_pipeline_launcher.R filters -i YOUR_DATASET --manual               # Optional step
Rscript 00_pipeline_launcher.R ctrl -i YOUR_DATASET
Rscript 00_pipeline_launcher.R process -i YOUR_DATASET --filter filtered --good_quality
Rscript 00_pipeline_launcher.R dea -i YOUR_DATASET
Rscript 00_pipeline_launcher.R combine -I YOUR_DATASET,YOUR_SECOND_DATASET
Rscript 00_pipeline_launcher.R process -i YOUR_COMBINED_DATASET
Rscript 00_pipeline_launcher.R dea -i YOUR_COMBINED_DATASET
Rscript 00_pipeline_launcher.R da -i YOUR_COMBINED_DATASET -m meld
```

## Getting Help

Need help? Post your question on the [issue board](https://github.com/BAUDOTlab/scDataPipeline/issues) and we'll do our best to answer.

## Contributing

As this project is still under active development, it might be hard to contribute from outside the team. However, suggestions are always welcomed!

## License

The scDataPipeline is distributed under the **MIT license**. See [LICENSE](./LICENSE) for more details.
