# classifyNMIBC

This package implements a Pearson nearest-centroid classifier that assigns class labels to single samples according to the four transcriptomic UROMOL2021 classes of non-muscle-invasive bladder cancer (NMIBC): class 1, class 2a, class 2b and class 3.

The classifier code was adapted from the consensusMIBC classifier: Kamoun, A et. al. A Consensus Molecular Classification of Muscle-invasive Bladder Cancer, Eur Urol (2019), doi: https://doi.org/10.1016/j.eururo.2019.09.006
Both classifiers can be found on our online web application: http://nmibc-class.dk

A smaller, example data set is provided to run the classifier.


## Citation 
Please cite Lindskrog and Prip et al. An integrated multi-omics analysis identifies prognostic molecular subtypes of non-muscle-invasive bladder cancer. Nat Commun. 2021. PMID: 33863885. DOI: 10.1038/s41467-021-22465-w


## Install
You may install this package with devtools:
``` {r}
library(devtools)
devtools::install_github("sialindskrog/classifyNMIBC", build_vignettes = TRUE)
library(classifyNMIBC)
```

## Usage
``` {r}
classifyNMIBC(x, minCor = .2, gene_id = c("ensembl_gene_ID", "hgnc_symbol")[1])
```
'x': dataframe with unique genes in rows and samples to be classified in columns (or single named vector of gene expression values).
RNA-seq data needs to be log-transformed and micro-array data should be normalized. Gene names may be supplied as Ensembl gene IDs or HUGO gene symbols.

'minCor': a numeric value specifying a minimal threshold for best Pearson's correlation between a sample's gene expression profile and centroids profiles. A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2.

'gene_id': a character value specifying the type of gene identifiers used for the names/rownames of 'x', ensembl_gene_ID for Ensembl gene IDs or hgnc_symbol for HUGO gene symbols.


## Example

``` {r}
data(test_data)

NMIBC_class <- classifyNMIBC(test_data)

head(NMIBC_class)
#       NMIBC_class cor_pval separationLevel   Class_1  Class_2a  Class_2b   Class_3
# U0026     Class_1        0       0.4596168 0.8722377 0.7244603 0.8402362 0.7649863
# U1270     Class_1        0       0.8715953 0.9076109 0.7632527 0.7878480 0.8151044
# U0062     Class_1        0       0.7756813 0.8935798 0.7522703 0.7440653 0.8040516
# U0268     Class_1        0       0.1017993 0.8979415 0.7236878 0.8106683 0.8932611
# U1031     Class_1        0       0.8399565 0.8918952 0.7297127 0.8244154 0.8430349
# U2111     Class_3        0       0.6588400 0.6799160 0.7345242 0.6853025 0.7820521

```
The classifier returns a dataframw with 7 columns:

'NMIBC_Class': the predicted class labels of the samples.

'cor_pval': the p-value associated with the Pearson's correlation between the sample and the nearest centroid.

'separationLevel': gives a measure (ranging from 0 to 1) of how a sample is representative of its consensus class, with 0 meaning the sample is too close to the other classes to be confidently assigned to one class label, and 1 meaning the sample is very representative of its class. This separationLevel is measured as follows: (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.

The remaing four columns return the Pearson's correlation values for each sample and each centroid.
