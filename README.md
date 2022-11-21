# Enrichment_Analysis

This repo can be cloned and run locally to perform enrichment analysis across
multiple pathway databases and datasets.

**The main script `Enrichment_Analysis.R` will initialize, check databases and
check dataset files prior to launching an interactive session for user input. 
Result files will be organized under a run-specific subfolder in `output`.**

```
source("Enrichment_Analysis.R")
```
*Note: You will need to prepare a database and at least one dataset for this script to work. See below...*

## Databases
And a sample GMT database file is also provided from [WikiPathways](https://new.wikipathways.org/download.html).
You can include any GMT files that you want. Here are some suggested sources:
 * ...
 
You will need to construct an RData file containing processed versions of your GMTs. You only need to 
do this once; it will be reused by default in all subsequent runs, until you choose to update or replace it.

```
source("Build_Database.R")
```
*Note: Regularly update your GMT source files as new versions are released*

## Datasets
A sample dataset `E-GEOD-30573.csv`([Voineagu 2011](https://europepmc.org/article/MED/21614001)) is provided 
to illustrate the required format. Essentially, a `gene` column is required. Columns named `rank`, `p.value`, 
and `fold.change` are also recognized for more sophisticated analyses. All other columns are ignored.  

You can include as many dataset files as you like as well. These can be:
 * a simple list of genes (for ORA only)
 * a list of genes with a `rank` column (for GSEA only)
 * a list of genes with a `p.value` column (for GSEA and ORA)
 * a list of genes with `p.value` and `fold.change` columns (for GSEA and ORA)

*Note: For a given run, all datasets are assumed to have the same format. If you have different types of dataset files 
(per the list above), then you should run them in separate batches.*

### Identifers  
Supported gene identifiers include SYMBOL, ENSEMBL, ENTREZID, and UNIPROT, by default. But 
it is trivial to add support for [any orgDb keytypes](http://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html).
Simply edit the `supported.idTypes` list in the "Defaults and Parameters" section of
the main script `Enrichment_Analysis.R`.

### Organisms
Supported species include human, mouse and rat, by default. But it is
trivial to add support for [any orgDb species](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb).
Simply edit the `supported.org` list in the "Defaults and Parameters" section of
the main script `Enrichment_Analysis.R`.

## Output
Results are organized under a run-specific folder named by timestamp (by default) 
with subfolders for `gsea` and `ora` results, each of which includes a folder 
of `plots`. Results include TSV and XLSX tables of enriched terms and pathway, as well
as RDS versions of the complete input and output for each clusterProfiler run. Plots
include volcano, dot, cnet, heatmap, upset, and more.

**A Shiny app is available to interactively explore results and download individual 
files (see app dir). Its UI is organized in the same way as the output directory.**
