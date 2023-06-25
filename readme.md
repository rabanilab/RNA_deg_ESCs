# Analysis of mRNA degradation in transcriptional shut-off experiments

This Matlab package is used to analyze mRNA degradation rates within
transcriptional shut-off experiments.

calculate degradation rate based on normalized TPM values.

You can find the full method description in:

[developmental transcripts during murine embryonic stem cell differentiation via CAPRIN1-XRN2](https://www.ncbi.nlm.nih.gov/pubmed/36495875).
Viegas JO, Azad JK, Lv Y, Fishman L, Paltiel T, Pattabiraman S, Park JE, Kaganovich D, Sze SK, Rabani M, Esteban MA, Meshorer E.
Dev Cell. 2022.

## Installation

These instructions will get you a copy of the analysis software.
We provide both matlab code to use with Matlab on your local machine, 
and pre-compiled executable to use without matlab.

### Prerequisites

We developed and tested this package with Matlab R2018b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).

### Installing pre-compiled version
Determine in which directory you would like to install the package. 
This will be the installation directory.

Download the package source code from GitHub. Copy the zip file into the
installation directory, and unzip the package.

```
unzip RNA_deg_ESCs-master.zip
```

The executables are available to use.


### Installing on Matlab

Determine in which directory you would like to install the package. 
This will be the installation directory.

Download the package source code from GitHub. Copy the zip file into the
installation directory, and unzip the package.

```
unzip RNA_deg_ESCs-master.zip
```

Import the intstallation directory into matlab, by running (in Matlab):

```
addpath <path to installation directory>;
```

Now you can run in matlab all the source code that is saved in the 
installation directory (see example below).

## Example data: pre-compiled

The package that you downloaded includes an example for using the package.

To run the example via pre-compiled executables, browse to the installation directory. Then run:

```
>> make run_WT_example run_KO_example
```

In the example, we analyze two datasets: WT and KO.
This analysis will create 2 directories with the results of each of the 2 analyses:

```
example_KO
example_WT
```

The analyses directories contain the following files:

```
degradation_rates.jpg - histogram plots of degradation rates and initial expression values
degradation_rates.txt - tab delimited list of degradation rates and r-square values per gene
degradation_rates.rsq.txt - tab delimited list of degradation rates and r-square values per gene, only for genes with r-square value above 0.5
stable_genes.jpg - plot of expression levels of stable genes selected for normalization
stable_genes.txt - list of stable genes selected for normalization
```
Degradation rates are given as half-life (time until half of the initial level is reached).


## Example data: on Matlab

To run the example via matlab, browse to the installation directory. Then run:
```
>> run_example
```

In the example, we analyze two datasets: WT and KO either independently,
or relative to each other.
This analysis will create 3 directories with the results of each of the 3 analyses:

```
example_KO
example_WT
example_pairwise
```

Two directories of dataset independent analyses contain the following files:

```
degradation_rates.jpg - histogram plots of degradation rates and initial expression values
degradation_rates.txt - tab delimited list of degradation rates and r-square values per gene
degradation_rates.rsq.txt - tab delimited list of degradation rates and r-square values per gene, only for genes with r-square value above 0.5
stable_genes.jpg - plot of expression levels of stable genes selected for normalization
stable_genes.txt - list of stable genes selected for normalization
```
Degradation rates are given as half-life (time until half of the initial level is reached).


Directory of pairwise analyses contain in addition the following files:

```
degradation_rates.fasterdg.txt - list of genes with faster degradation in treatment (T) vs control (C)
degradation_rates.slowerdg.txt - list of genes with slower degradation in treatment (T) vs control (C)
degradation_rates.xy.jpg - plots of correlation between parameters in treatment (T) vs control (C)
```
## k-mer analysis

To run the k-mer analysis on the resulting degradation rates, use the following steps:
1. Download the [k-mer analysis package](https://github.com/rabanilab/cont_kmer_analysis) and follow instructions there for installation.

2. Convert the output file so that it is formatted appropriately:
```
>> make convert_file_format path=/path/to/degradation_rates.rsq.txt
```
3. Run the following k-mer analysis commands
```
>> make run_ks_test_job path_to_kmer_matrices=3utr/kmer_matrices PARAMETERS_FILE=parameters.tsv output_path=3utr/kmer_out term=half_life 
>> make ks_tests_plots ks_test_output_folder=3utr/kmer_out term=half_life output_path=3utr/kmer_out/plots
```

## License

This project is licensed under the MIT License - see source files for details.
