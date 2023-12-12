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

- `conda`

Run the following command to install all the packages needed for the pipeline in a new environment:

```bash
conda env create -f environment.yml
conda activate scDataPipeline
```

This may take several minutes. Some packages must be installed manually in a R console:

- [scCustomize](https://github.com/samuel-marsh/scCustomize)
- [SeuratDisk](https://github.com/mojaveazure/seurat-disk)

```R
install.packages("scCustomize")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
```

### Installation

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

### Filling the parameter files

There is two types of parameter files: **globalParameters_MODEL_single.toml** for single datasets and **globalParameters_MODEL_combined.toml** for combined datasets. Note that in both case, the file must have the exact value of `DATASET` in its name.

#### Regarding the paths

`PATH_ROOT` must contain the absolute path towards the root folder of our project. The other paths included in `path.input` and `path.output` **MUST BE** paths relative to `PATH_ROOT`, with the exception of `PATH_ATLAS_FILE`, which is relative to `PATH_ATLAS`.

## How to use

There is only one script you need to launch each time. Go into the scDataPipeline folder, then run:

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
Rscript 00_pipeline_launcher.R da -i YOUR_DATASET -m meld
```

## Getting Help

Need help? Post your question on the [issue board](https://github.com/BAUDOTlab/scDataPipeline/issues) and we'll do our best to answer.

## Contributing

As this project is still under active development, it might be hard to contribute from outside the team. However, suggestions are always welcomed!

## License

The scDataPipeline is distributed under the **MIT license**. See [LICENSE](./LICENSE) for more details.
