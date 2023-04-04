# Analysis of mRNA degradation in transcriptional shut-off experiments

This Matlab package is used to analyze mRNA degradation rates within
transcriptional shut-off experiments.

calculate degradation rate based on normalized TPM values.

You can find the full method description in:

[developmental transcripts during murine embryonic stem cell differentiation via CAPRIN1-XRN2](https://www.ncbi.nlm.nih.gov/pubmed/36495875).
Viegas JO, Azad JK, Lv Y, Fishman L, Paltiel T, Pattabiraman S, Park JE, Kaganovich D, Sze SK, Rabani M, Esteban MA, Meshorer E.
Dev Cell. 2022.

## Installation

These instructions will get you a copy of the analysis software 
to use with Matlab on your local machine. 

### Prerequisites

We developed and tested this package with Matlab R2018b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).

### Installing

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

## Example data

The package that you downloaded includes an example for using the package.

To run the example, browse to the installation directory. Then run:
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


## License

This project is licensed under the MIT License - see source files for details.
