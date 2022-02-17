
<!-- README.md is generated from README.Rmd. Please edit that file -->

#  Karl's HIPPO fork

I wanted to run HIPPO on a large data set.

## What's new here
- **Added some new parameters for user** add `max_features` and `min_cluster_size` as arguments to give more control to users on these parameters that affect run time
- **Heirarchecal (rather than sequential) feature subsetting**
	- In the current implementation of HIPPO, $X$ gets subsetted to the features selected in the previous round. The sequence of clustering steps get performed on a nested sequence of features. This doesn't make sense to me... a few clustering steps may be used to dig into one cell populations, and I don't want to restrict my features to those relevant to those small cell populations when I go up the tree and try to split an earlier cluster.
	- I changed the behavior of HIPPO to restrict to the significant features at the level the cluster was last split. This way, as you go down any particular lineage, the features used to split each child are a subset of the features used to split the parent, but it does not enforce the descendants of siblings to select on the same features.
- **Sparse log plus 1** calling $log(X + 1)$ is slow because you coerce the sparse matrix into a dense one. Performing SVD on the dense matrix is also more expensive. updated code to do this step efficiently.
- **Only perform one SVD**. 
	- **NOTE:  Best I can tell, the current release of HIPPO is not centering/scaling for clustering**. The "PCs" that get clustered are just the right singular vectors output from the `irlba::irlba`, which is run without the `ceter` or `scale.` arguments.  This is all well for me though, because I can now get the "unscaled" and "scaled" PCs with only one partial SVD. 
	- The current implementation calls "unscaledpc" the PCs of $log(X + 0.1)$ and the "scaled PCs" the right singular vectors of $\log(X+1)$. I just compute (partial SVD) of $\log(X+1) = UDV^T$ and call "unscaled PCs" $VD$, and the "scaled pcs" $V$. Ignoring the difference in pseudocount, this is identical to what is currently being done. It is MUCH faster. We also avoid an expensive and unecessary transpose operation on a `dgCMatrix` (it needs to build a new index, reorganize non-zero entries)
- **Dont Stop Early** HIPPO terminates as soon as it finds one cluster that fails the split criterea (proportion of significant features AND minimum # of cells in new cluster). This can leave many other clusters untouched. Instead, if we hit an unsplittable cluster, move to the next best candidate to split. This way the program keeps trying to split clusters until all clusters are "atomic" or maximum number of clusters $K$ is reached. The user can make decisions late about how far down the tree they'd like to go.

## Other comments

In the interest of time I skip computing the dispersion statistics. Could probably make this more efficient too

In order to get HIPPO to skip clusters that already attempted splitting but hit termination criteria I set `withinss[oldk] <- 0` , if this gets recorded anywhere in the final output (and if not I think it should) it will be misleading

There are some functions like `chunkRowSum(...)` that are not being used, I was messing around with operating on representations of on-disk data e.g `DelayedArray` but abandoned that.

# HIPPO

Single cell UMI analysis tool that focuses on zero-inflation to detect
biological heterogeneity

## Getting Started

These instructions will get you a copy of the project up and running on
your local machine for development and testing purposes. See deployment
for notes on how to deploy the project on a live system.

Prerequisites 
HIPPO works on the SingleCellExperiment object. You can install the library like the following.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

## Installing

You can download the most recent version of HIPPO package as below. Please allow up to 5 minutes to
completely compile the vignette.

``` r
devtools::install_github("tk382/HIPPO", build_vignettes = TRUE)
```

## Read the data

Many single-cell data sets used in the manuscript are available in
DuoClustering2018 package
([link](https://www.bioconductor.org/packages/release/data/experiment/html/DuoClustering2018.html)).

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("DuoClustering2018")
```

For example, you can load Zhengmix4eq data set like below.

``` r
sce = DuoClustering2018::sce_full_Zhengmix4eq(metadata = FALSE)
```

For this vignette, we use a smaller toydata, a subset of Zhengmix4eq.
The toydata is included in the HIPPO package.

``` r
data(toydata)
sce = toydata
```

Alternatively, you can start from a matrix object and create
SingleCellExperiment object.

``` r
X = readRDS("../zhengmix4eq_counts.rds")
sce = SingleCellExperiment(assays = list(counts = X))
```

## Diagnostic Plot

This plot shows the zero inflation compared to the expected Poisson
line. If most genes don’t align with the black line, it shows that there
is cell heterogeneity driving the zero inflation.

``` r
hippo_diagnostic_plot(sce, 
                      show_outliers = TRUE, 
                      zvalue_thresh = 3)
```

![](README_files/figure-gfm/diagnostic-1.png)<!-- -->

## Feature Selection and Hierarchical Clustering

HIPPO assumes that the count matrix is placed in
<sce@assays@data>$counts. Some objects that we found online have the
count matrix in <sce@assays>$data$counts. In this case, HIPPO will throw
an error because it cannot found a count matrix. In this case, you have
to create another SingleCellExperiment object to assign the count matrix
in the correct slot.

Next, you can run hippo function to do the pre-processing that
simutlaneously conducts feature selection and hierarchcial clustering.
There are three arguments that help you decide the stopping criterion of
clustering procedure.

K is the maximum number of clusters that you want. HIPPO will return the
clustering results for all k = 2, 3, …, K, so you can overestimate the
number of potential clusters. The default is 10, but users are highly
recommended to adjust this.

z\_threshold is the feature selection criterion. For each round of
hierarchical clustering, hippo will find outlier genes where the z-value
of significance is greater than the threshold. For example, if you would
like to select genes with p-values less than 0.05, z\_threshold would be
1.96. The default threshold is 2, but users can use their discretion to
change this value.

outlier\_proportion is the number of outlier genes to allow. The default
is 0.01 (1%) which means the clustering procedure will automatically
stop if there are less than 1% of genes remain as important features.
With the example data set, the default choice has empirically worked
well.

``` r
set.seed(20200610)
sce = hippo(sce, 
            feature_method = "zero_inflation",
            clustering_method = "kmeans",
            K = 3, 
            outlier_proportion = 0.00001)
```

## Dimension Reduction for Each Round of HIPPO

We offer two dimension reduction methods: umap and tsne. And we offer
two separate visualization functions.

``` r
sce = hippo_dimension_reduction(sce, method="umap")
hippo_umap_plot(sce)
```

![](README_files/figure-gfm/umap-1.png)<!-- -->

``` r
sce = hippo_dimension_reduction(sce, method="tsne")
hippo_tsne_plot(sce)
```

![](README_files/figure-gfm/tsne-1.png)<!-- -->

## Visualize the selected features at each round

This function shows how the zero-inflation decreases as HIPPO proceeds
in the clustering. This function has arguments called switch\_to\_hgnc
and ref. These aim to provide the users an option to change the gene
names from ENSG IDs to HGNC symbols for ease of understanding. Many
SingleCellExperiment objects have such data embedded in rowData(sce).
Users can create a data frame with ensg and hgnc columns for the genes,
and HIPPO will automatically switch the row names of the count matrix
from ENSG IDs to HGNC symbols. The default is set to FALSE, assuming
that the row names are already HGNC symbols.

``` r
data(ensg_hgnc)
zero_proportion_plot(sce, 
                     switch_to_hgnc = TRUE, 
                     ref = ensg_hgnc)
```

![](README_files/figure-gfm/featureselection-1.png)<!-- -->

``` r
hippo_feature_heatmap(sce, k = 2, 
                      switch_to_hgnc = TRUE, 
                      ref = ensg_hgnc, 
                      top.n = 20)
```

![](README_files/figure-gfm/featureselection-2.png)<!-- -->

``` r
hippo_feature_heatmap(sce, k = 3, 
                      switch_to_hgnc = TRUE, 
                      ref = ensg_hgnc, 
                      top.n = 20)
```

![](README_files/figure-gfm/featureselection-3.png)<!-- -->

## Differential Expression Example

We also offer a differential expression analysis tool.

This function also has an option to switch the gene names to HGNC
symbols. top.n argument lets users choose how many top genes to show in
the box plot. The default is 5.

The labels of boxplots are aligned with the t-SNE or UMAP plots above.
When K is equal to 2, the color codes match with the cell groups as
separated in the dimension reduction plot.

``` r
sce = hippo_diffexp(sce, 
                  top.n = 15, 
                  switch_to_hgnc = TRUE, 
                  ref = ensg_hgnc)
```

![](README_files/figure-gfm/diffexp-1.png)<!-- -->

Each round of differential expression test results are also saved in the
list of data frames.

``` r
head(get_hippo_diffexp(sce, 1))
#>             genes  null_dev   alt_dev          pval
#> 1 ENSG00000101439 1005.8606 418.53653 9.570942e-130
#> 2 ENSG00000087086  677.0435 249.54309  5.687571e-95
#> 3 ENSG00000251562 1097.4217 709.70432  2.599070e-86
#> 4 ENSG00000011600  440.4111  87.47779  9.735996e-79
#> 5 ENSG00000227507  518.0729 308.59588  1.786349e-47
#> 6 ENSG00000137154  551.1884 344.74232  8.189475e-47
head(get_hippo_diffexp(sce, 2))
#>             genes  null_dev  alt_dev          pval
#> 1 ENSG00000019582 1701.8180 166.8158  0.000000e+00
#> 2 ENSG00000204287 1088.8127 112.4015 2.408126e-214
#> 3 ENSG00000223865  904.2276 142.0947 9.232007e-168
#> 4 ENSG00000251562 1097.4217 364.1340 1.727066e-161
#> 5 ENSG00000196126  765.7430 146.7529 1.239977e-136
#> 6 ENSG00000087086  677.0435 158.7387 9.893622e-115
```
