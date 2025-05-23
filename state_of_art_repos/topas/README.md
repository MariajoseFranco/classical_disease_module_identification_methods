<p align="center">
  <img src="TOPAS_logo.png"
  alt="drawing" 
  style="float: right;"
  />
</p>

TOPAS (TOP-down Attachment of Seeds) is an iterative approach which aims at connecting the largest number of seed nodes in a top-down fashion through connectors which guarantee the highest flow of a Random Walk with Restart in a network of functional associations.

* Version 1.0

### Step 1: Installation

TOPAS is implemented in R (v4.1), and it was tested in a version-controlled conda environment. 

Clone the TOPAS repository and follow the instructions to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and then execute the following to generate a conda environment with all you need to run TOPAS:
```
conda env create -f environment.yml
```

### Step 2: Program Execution

Once you have created a new conda environment, make sure to activate it before running TOPAS. 
```
conda activate topas
```

You can copy this R snippet into your local computer to produce a disease module as an example. You may want to set the path to match the working directory where you cloned the TOPAS repository.

```{r}
#!/usr/bin/Rscript

setwd('~/topas') 
source('src/TOPAS.R')

## Load network
network_example = 
  readr::read_tsv(file = 'example/FunCoup5_pfc80.tsv.gz', 
                  col_names = c('V1','V2'), 
                  col_types = 'cc')
## Seed gene set
seeds_example = 
  readr::read_table(file = 'example/adrenal_gland_diseases.txt', 
                    col_names = 'V',
                    col_types = 'c')

module_example = 
  TOPAS(
    network = network_example,
    seeds = seeds_example$V,
    expansion_steps = 2,
    cores = 4
    )

readr::write_tsv(module_example, file = 'example/TOPAS_adrenal_gland_diseases.tsv')
```

**NOTE**: TOPAS is only designed to process undirected networks!


### Reproduce the benchmark

In order to reproduce all figures in the study, follow the instructions below. 

1. Unzip the file *benchmarkData.zip* and access the folder *benchmarkData*:

```
 unzip benchmarkData.zip
 cd benchmarkData
```

2. Execute the following to generate a conda environment with all you need to execute the notebooks:

```
conda env create -f environment_benchmark.yml
```

3. Once you have created a new conda environment, make sure to activate it before executing the notebooks. 

```
conda activate topas_benchmark
```

4. To execute any notebook from command line, enter the folder dedicated to one specific session of the study and render the notebook called *FiguresGenerator.Rmd*.

```
## For example
cd AssessmentOfSRR\&SCR
Rscript -e "rmarkdown::render('FiguresGenerator.Rmd')"
```

**NOTE**: All figures are pre-generated and stored in dedicated folders under *figures/* and displayed in the *html* files. However, you are welcome to reproduce the results executing the notebooks.

### Contacts ###

* Davide Buzzao (davide.buzzao@scilifelab.se)
* Erik L.L. Sonnhammer (erik.sonnhammer@scilifelab.se)