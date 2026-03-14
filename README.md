# IDGDD
Identify Driver Genes based on Differentially expressed gene modules and Dual random walk
## Installation

Please ensure that R Studio are installed before running the code.  
Then you can install some packages that are necessary to run IDGDD code:

```{r, eval=FALSE}
library(psych)
library(igraph)
```

## Example

This is a basic example which shows you how to get a list of potential driver genes for breast cancer:

```{r, eval=FALSE}
library(psych)
library(igraph)
```
Step1: load data about breast cancer:
```{r, eval=FALSE, message=FALSE, warning = FALSE}
mul_edge_list <- read.table('edge_file.txt')
mul_gene_list <- read.table('point_file.txt')
tf_gene <- read.table('TF_trrust_rawdata.human.txt')
MIR <- read.table("mirna_MIR.txt")
PATHWAY <- read.csv("kegg_pathway.csv",F)
G_mut <- read.table('mut_brca_18846.txt',sep='\t',TRUE)
def_gene <- read.csv("brca_degs_05.csv")
mrna_exp <- read.table('Exp_brca.txt',sep='\t',TRUE)
dif_module <- read.csv("brca_nodewl60_dif_module10_1408_0001.csv")

```

Step2: Calculate the score for each gene and determine gene priority:
```{r, eval=FALSE, message=FALSE, warning = FALSE}
AdjMatrix <- Construct_PPINET(mul_edge_list,mul_gene_list,tf_gene)

brca_list <- IDGDD(AdjMatrix,G_mut,MIR,dif_module,0.4)
```
