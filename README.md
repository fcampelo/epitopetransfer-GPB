# EpitopeTransfer

This repository contains data and code for the paper **EpitopeTransfer:
Improving linear B-cell epitope prediction through phylogeny-aware
transfer learning** by Lindeberg Pessoa Leite, Teófilo E. de Campos,
Francisco Pereira Lobo and Felipe Campelo.

## Table of Contents

-   [Dependencies](#Dependencies)
-   [Running models for published metrics](#Running-models-for-published-metrics)
-   [Running the complete pipeline](#Running-the-complete-pipeline)
-   [Running the analyses](#Running-the-analyses)
-   [Contact](#contact)

## Dependencies 

This project requires Python 3.10, although other versions are likely
compatible. Follow the instructions below to set up your environment.

1.  **Install Python 3.10**

Ensure that Python 3.10 is installed on your system. If not, you can install it using your system's package manager. For Ubuntu-based distributions, you can use:

``` bash
    sudo apt update
    sudo apt install software-properties-common -y
    sudo add-apt-repository ppa:deadsnakes/ppa -y
    sudo apt update
    sudo apt install python3.10 python3.10-venv python3.10-distutils -y
```

2.  **Activate epitopetransfer environment**

``` bash
  source esm2/bin/activate
  # or
  source esm1b_v1/bin/activate
  # or
  source esm1b_v2/bin/activate
```
To exit any environment, run: deactivate

## Running models for published metrics

To run the models for specific taxa or for all taxa included in the
study, use the following command format in the terminal:

``` bash
  python3.10 main.py [base_model] [taxa]
```

Replace [base_model] to esm1b or esm2, and [taxa] to the desired taxa or all (for esm2 only). Example:

``` bash
  source esm2/bin/activate
  (esm2) python3.10 main.py esm2 all # (the 'all' option is available for esm2 base model only)
```
To exit the `esm2` environment, use this command: `deactivate`

``` bash
  source esm1b_v1/bin/activate
  (esm1b_v1) python3.10 main.py esm1b bpertussis
```

``` bash
  source esm1b_v2/bin/activate
  (esm1b_v2) python3.10 main.py esm1b mononegavirales
```

Available Taxa

| **Taxa**                                        | **Taxa**                                                 |
|-----------------------------------|-------------------------------------|
| **bpertussis**: *Bordetella pertussis*          | **filoviridae**: *Filoviridae*                           |
| **corynebacterium**: *Corynebacterium*          | **ovolvulus**: *Onchocerca volvulus*                     |
| **orthopoxvirus**: *Orthopoxvirus*              | **ctrachomatis**: *Chlamydia trachomatis*                |
| **ecoli**: *Escherichia coli*                   | **human_gammaherpesvirus_4**: *Human Gammaherpesvirus 4* |
| **enterobacteriaceae**: *Enterobacteriaceae*    | **influenza_a**: *Influenza A*                           |
| **lentivirus**: *Lentivirus*                    | **cdifficile**: *Clostridioides difficile*               |
| **mtuberculosis**: *Mycobacterium tuberculosis* | **measles_morbilivirus**: *Measles morbillivirus*        |
| **paeruginosa**: *Pseudomonas aeruginosa*       | **mononegavirales**: *Mononegavirales*                   |
| **smansoni**: *Schistosoma mansoni*             |                                                          |
| **tgondii**: *Toxoplasma gondii*                |                                                          |
| **pfalciparum**: *Plasmodium falciparum*        |                                                          |


**Note:**  For the esm1b base model, activate the esm1b_v1 environment (source esm1b_v1/bin/activate) for taxa in the first column, or the esm1b_v2 environment (source esm1b_v2/bin/activate) for taxa in the second column. For the esm2 base model, activate the esm2 environment (source esm2/bin/activate) for all taxa. 


## Running the complete pipeline

All steps of the pipeline for the *B. pertussis* taxon are available in the `pipeline` folder. The same procedure applies to the other taxa as well.

Make sure to download the [ESM-2 model](https://huggingface.co/facebook/esm2_t33_650M_UR50D) or just point to this online repository, and adjust the file paths accordingly before execution.

To run the complete process, execute the following notebooks in sequence:

1. `dataset_generation_taxa.ipynb`  
2. `finetune_higher_level.ipynb`  
3. `feature_calculation_reduction.ipynb`


## Initial data generation and post-experimental analysis
The scripts to retrieve the IEDB export and pre-process, filter and consolidate
the LBCE data are available in the `./R` folder. The analysis scripts used to 
process the results of the experiments are also contained in this folder. 
Dependencies are managed using `renv`. Open project `./R/data_extraction_routines.Rproj` and run 
`renv::restore()` to install the required packages.

*****
### Contact Email:  
[Felipe Campelo](mailto:f.campelo@bristol.ac.uk) - Principal investigator  
[Lindeberg Leite](mailto:lindpessoa@gmail.com) - First author
