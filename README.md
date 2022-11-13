# Enrichment_Analysis

This project can be cloned and run locally to perform enrichment analysis across
multiple pathway databases and datasets.

The main script `Enrichment_Analysis.R` will initialize and check databases and
dataset files prior to launching an interactive session for user input. Once 
launched, result files will be organized under a subfolder in `output`.

A sample dataset `E-GEOD-30573.csv` is provided to illustrate the format required
for dataset files ([Voineagu 2011](https://europepmc.org/article/MED/21614001)). 
And a sample database file `WP_.gmt` is also provided from [WikiPathways](https://new.wikipathways.org/download.html).

You can include any GMT you like. Here are some suggested sources:
 * ...
 
You can include as many dataset files as you like as well. These can be:
 * a simple list of gene names (for ORA only)
 * a list of genes with a `rank` column (for GSEA only)
 * a list of genes with a `p.value` column (for GSEA or ORA)
 * a list of genes with `p.value` and `fold.change` columns (for GSEA or ORA)
  
Supported species include human, mouse and rat, by default. But it would be
trivial to add support for [any orgDb species](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb).
Simply edit the `supported.org` list in the "Defaults and Parameters" section of
the main script `Enrichment_Analysis.R`.
  
Output is organized under a run-specific folder named by timestamp (by default) 
with subfolders for `gsea` and `ora` results, each of which includes a folder 
of `plots`.