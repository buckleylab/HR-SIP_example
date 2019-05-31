MW-HR-SIP over multiple samples
================
Samuel Barnett
May 31, 2019

-   [Introduction](#introduction)
-   [MW-HR-SIP](#mw-hr-sip)
-   [Session Info](#session-info)

Introduction
------------

This tutorial runs through running Multiple Window High Resolution DNA-SIP (MW-HR-SIP) with a real example dataset that includes multiple treatment and control samples. The sample code in [Chapter\_Examples](Chapter_Examples.md) shows how to run a simple analysis on one treatment and control comparison. However, in a real world experiment, you will be working with more data than this, such as using multiple isotopically labled substrates, timepoints, community types. Running MW-HR-SIP using the function `HRSIP()` from package `HTSSIP` is still simple to run with multiple treatment-control comparisons, and just requires a different data management technique ... lists. Basically, instead of running `HRSIP()` on a single phyloseq object, you will make a list of phyloseq objects and then iteratively run `HRSIP()` on each element. This also allows you to use an alternative parallelization for systems with multiple processors that can make your analysis much faster. The files used in this tutorial can be found with this github site under the directory [example\_data](example_data/).

For this example the sample data is from an experiment that added glucose and cellulose to soil microcosms and harvested them 3 or 14 days after substrate addition. This means there are both a multiple substrate and timepoint components. In this experiment, both substrates are added to all microcosms. The treatment microcosms have the distinction that one of the substrates is labeled with 13C. The power of this design is that it allows you to use the same controls for both substrates within the same timepoint. The microcosms are defined as follows:

-   Treatment microcosms
    -   13C-Glu.D3: Given 13C-glucose and 12C-cellulose and harvested on day 3
    -   13C-Cel.D3: Given 13C-cellulose and 12C-glucose and harvested on day 3
    -   13C-Glu.D14: Given 13C-glucose and 12C-cellulose and harvested on day 14
    -   13C-Cel.D14: Given 13C-cellulose and 12C-glucose and harvested on day 14
-   Control microcosms
    -   12C-Con.D3: Given 12C-glucose and 12C-cellulose and harvested on day 3
    -   12C-Con.D14: Given 12C-glucose and 12C-cellulose and harvested on day 14

### R packages needed

``` r
# Packages needed for data handling
library(dplyr)
library(tidyr)
library(tibble)

# Packages needed for analysis
library(phyloseq)   # Used for handling our data format
library(HTSSIP)     # Contains the main methods used in this analysis

# Packages used to make this Rmarkdown notebook look nice
library(knitr)
library(kableExtra)
```

MW-HR-SIP
---------

This example doesn't include any prelimiary analyses, but it is recommended that some of the analyses found in [Chapter\_Examples](Chapter_Examples.md) and [addl\_prelim\_analyses](addl_prelim_analyses.md) be conducted prior to running MW-HR-SIP. These preliminary analyses can easily be run using multiple treatment-control sample pairs and examples are found in the linked pages.

#### 1. Import data

Data for this tutorial is a dataset called "example\_S2D2\_phyloseq.rds". As before this is an R object containing the data in phyloseq format.

``` r
# Import the data you using the readRDS() function
SIP.physeq <- readRDS("example_data/example_S2D2_phyloseq.rds")

# What does this phyloseq object look like?
SIP.physeq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10001 taxa and 139 samples ]
    ## sample_data() Sample Data:       [ 139 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10001 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10001 tips and 10000 internal nodes ]

#### 2. Subset data by treatment-control comparisons

Currently all of the data is in one single phyloseq object. For the analysis you need separate phyloseq objects, each containing all the fractions for a single treatment and its corresponding control. One way to do this is to just subset out a desired treatment-control pair using phyloseq function `subset_samples()`. Alternatively if you want to run MW-HR-SIP on all comparisons in one step you can convert the phyloseq object into a list of single treatment-control objects that can be itteratively run through `HRSIP()`. This latter method is shown below.

In this case you have two carbon substrates and two days. This will result in four treatment-control comparisons:

13C-Cellulose Day 3 vs. 12C-Control Day 3 13C-Glucose Day 3 vs. 12C-Control Day 3 13C-Cellulose Day 14 vs. 12C-Control Day 14 13C-Glucose Day 14 vs. 12C-Control Day 14

This means you need to subset the data by `Substrate` and by `Day`. To do this you make an expression that tells the function how to pair up samples. You will use the expression:

`(substrate=='12C-Con' & day=='${day}') | (substrate=='${substrate}' & day=='${day}')`

This expression essentially means that you group samples with the same `day` value and with either `12C-Con` or a distinct other `substrate` value.

You also need a set of the different pairs of parameters that will be used to group samples. In this case all combinations of `substrate` and `day`. This needs to only include treatment samples, so you first remove `12C-Con` from the `substrate` options.

``` r
# Set up the treatment-control pairing expression
ex <- "(substrate=='12C-Con' & day=='${day}') | (substrate=='${substrate}' & day == '${day}')"

# Get a set of subsetting parameters for the treatment samples
params <- get_treatment_params(SIP.physeq, c('substrate', 'day'), "substrate != '12C-Con'")

# Subset the data into a list of phyloseq objects, each for a different treatment-control comparison
SIP.physeq.list <- phyloseq_subset(SIP.physeq, params, ex)

# What does the resulting dataset look like?
SIP.physeq.list
```

    ## $`(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')`
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10001 taxa and 46 samples ]
    ## sample_data() Sample Data:       [ 46 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10001 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10001 tips and 10000 internal nodes ]
    ## 
    ## $`(substrate=='12C-Con' & day=='14') | (substrate=='13C-Cel' & day == '14')`
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10001 taxa and 46 samples ]
    ## sample_data() Sample Data:       [ 46 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10001 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10001 tips and 10000 internal nodes ]
    ## 
    ## $`(substrate=='12C-Con' & day=='14') | (substrate=='13C-Glu' & day == '14')`
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10001 taxa and 47 samples ]
    ## sample_data() Sample Data:       [ 47 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10001 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10001 tips and 10000 internal nodes ]
    ## 
    ## $`(substrate=='12C-Con' & day=='3') | (substrate=='13C-Glu' & day == '3')`
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10001 taxa and 46 samples ]
    ## sample_data() Sample Data:       [ 46 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10001 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10001 tips and 10000 internal nodes ]

As you can see you now have a list of 4 phyloseq objects, each containing a separate treatment and control sample set. Each list entry is named by the comparison being made. Now you can see how the expression you created above split up the phyloseq object.

#### 3. Generate list of buoyant density windows, list of sparsity threshold cutoffs, and set the p-value cutoff:

Before running the MW-HR-SIP you need to set some important parameters:

-   **Windows**: The overlaping bouyant density windows you want to analyze. You will be comparing the read counts of OTUs between the treatment and control fractions within these bouyant density windows. `density_min` refers to the minimum BD of each window. `density_max` refers to the maximum BD of each window. Make sure that there are at least 3 fractions each from the treatment and control within each window.
-   **Sparsity**: The sparsity thresholds used to remove OTUs found in very few fractions. Removing these increases your statistical power as it reduces the number of comparisons maked. You want to run this analysis with multiple sparsity threshold cutoffs and choose the one resulting in the most power or the most rejected null hypotheses.
-   **p-value cutoff**: The p-value cutoff below which an OTU is considered significantly enriched in the treatment compared to the control. OTUs with an adjusted p-value below this cutoff will be designated isotopically labeled.

``` r
# Set BD windows
windows <- data.frame(density_min=c(1.70, 1.72, 1.74), 
                     density_max=c(1.73, 1.75, 1.77))

# Set sparsity thresholds
sparsity_list <- c(0, 0.15, 0.30)

# Set pvalue_cutoff
pvalue_cutoff <- 0.05
```

#### 4. Run MW-HR-SIP:

Now you can run the MW-HR-SIP analysis on the data. Since the data has been split across a list of phyloseq objects, you can use `ldply()` to run the `HRSIP()` command iteratively over the list. Since you are comparing OTU counts between values of `substrate` (`12C-Con` controls and `13C-Glu` or `13C-Cel`), our design will be `~substrate`. This part may take some time.

One way to parallelize this if you have multiple processors is to use use set the flag `.parallel` equal to `true` within `ldply()` (i.e. `.parallel = TRUE`).

``` r
# Run MW-HR-SIP over the list of phyloseq objects
l2fc.df <- plyr::ldply(SIP.physeq.list, 
                      HRSIP, 
                      density_windows = windows,
                      design = ~substrate, 
                      padj_cutoff = pvalue_cutoff,
                      sparsity_threshold = sparsity_list)
```

    ## Sparsity threshold: 0 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold with the most rejected hypotheses: 0 
    ## Sparsity threshold: 0 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold with the most rejected hypotheses: 0.3 
    ## Sparsity threshold: 0 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold with the most rejected hypotheses: 0.15 
    ## Sparsity threshold: 0 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.7-1.73 
    ## Sparsity threshold: 0 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.72-1.75 
    ## Sparsity threshold: 0 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.15 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold: 0.3 
    ## Density window: 1.74-1.77 
    ## Sparsity threshold with the most rejected hypotheses: 0.3

Now you can view the results. You will notice that the first column `.id` contains an expression similar to the one set in step 2. This indicates the treatment-control comparison for that particular row. This column is not present when running a single phyloseq object. To filter out just the labeled OTUs you can filter your dataframe such that `padj` is less than or equal to your `pvalue_cutoff`.

``` r
# View out the first 10 results filtering to just include labeled OTUs or those with a padj <= pvalue_cutoff.
kable(head(l2fc.df[l2fc.df$padj <= pvalue_cutoff,], n=10), "html") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height="400px")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
.id
</th>
<th style="text-align:left;">
OTU
</th>
<th style="text-align:right;">
log2FoldChange
</th>
<th style="text-align:right;">
p
</th>
<th style="text-align:right;">
padj
</th>
<th style="text-align:left;">
Rank1
</th>
<th style="text-align:left;">
Rank2
</th>
<th style="text-align:left;">
Rank3
</th>
<th style="text-align:left;">
Rank4
</th>
<th style="text-align:left;">
Rank5
</th>
<th style="text-align:left;">
Rank6
</th>
<th style="text-align:left;">
Rank7
</th>
<th style="text-align:left;">
Rank8
</th>
<th style="text-align:right;">
density\_min
</th>
<th style="text-align:right;">
density\_max
</th>
<th style="text-align:right;">
sparsity\_threshold
</th>
<th style="text-align:left;">
sparsity\_apply
</th>
<th style="text-align:right;">
l2fc\_threshold
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.1136
</td>
<td style="text-align:right;">
3.110807
</td>
<td style="text-align:right;">
0.0008800
</td>
<td style="text-align:right;">
0.0442424
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
\_\_uncultured\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
171
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.995
</td>
<td style="text-align:right;">
3.607915
</td>
<td style="text-align:right;">
0.0000491
</td>
<td style="text-align:right;">
0.0054616
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
177
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.154
</td>
<td style="text-align:right;">
1.473143
</td>
<td style="text-align:right;">
0.0008285
</td>
<td style="text-align:right;">
0.0429532
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
\_\_uncultured\_Acidobacteria\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
184
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.6610
</td>
<td style="text-align:right;">
1.604973
</td>
<td style="text-align:right;">
0.0002933
</td>
<td style="text-align:right;">
0.0206198
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
\_\_uncultured\_Acidobacteria\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
223
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.91
</td>
<td style="text-align:right;">
1.385836
</td>
<td style="text-align:right;">
0.0007396
</td>
<td style="text-align:right;">
0.0402988
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
\_\_uncultured\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
253
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.167
</td>
<td style="text-align:right;">
1.628523
</td>
<td style="text-align:right;">
0.0003875
</td>
<td style="text-align:right;">
0.0257161
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Acidobacteria
</td>
<td style="text-align:left;">
\_\_DA023
</td>
<td style="text-align:left;">
\_\_uncultured\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
282
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.8567
</td>
<td style="text-align:right;">
3.081788
</td>
<td style="text-align:right;">
0.0000029
</td>
<td style="text-align:right;">
0.0006324
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Firmicutes
</td>
<td style="text-align:left;">
\_\_Bacilli
</td>
<td style="text-align:left;">
\_\_Bacillales
</td>
<td style="text-align:left;">
\_\_Bacillaceae
</td>
<td style="text-align:left;">
\_\_Bacillus
</td>
<td style="text-align:left;">
\_\_uncultured\_bacterium
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
308
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.13320
</td>
<td style="text-align:right;">
3.332413
</td>
<td style="text-align:right;">
0.0003463
</td>
<td style="text-align:right;">
0.0233572
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Firmicutes
</td>
<td style="text-align:left;">
\_\_Bacilli
</td>
<td style="text-align:left;">
\_\_Bacillales
</td>
<td style="text-align:left;">
\_\_Paenibacillaceae
</td>
<td style="text-align:left;">
\_\_Paenibacillus
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
309
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.277
</td>
<td style="text-align:right;">
2.326706
</td>
<td style="text-align:right;">
0.0001831
</td>
<td style="text-align:right;">
0.0146614
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Firmicutes
</td>
<td style="text-align:left;">
\_\_Bacilli
</td>
<td style="text-align:left;">
\_\_Bacillales
</td>
<td style="text-align:left;">
\_\_Paenibacillaceae
</td>
<td style="text-align:left;">
\_\_Paenibacillus
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
310
</td>
<td style="text-align:left;">
(substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
</td>
<td style="text-align:left;">
OTU.246
</td>
<td style="text-align:right;">
2.683650
</td>
<td style="text-align:right;">
0.0000118
</td>
<td style="text-align:right;">
0.0018075
</td>
<td style="text-align:left;">
Bacteria
</td>
<td style="text-align:left;">
\_\_Firmicutes
</td>
<td style="text-align:left;">
\_\_Bacilli
</td>
<td style="text-align:left;">
\_\_Bacillales
</td>
<td style="text-align:left;">
\_\_Paenibacillaceae
</td>
<td style="text-align:left;">
\_\_Paenibacillus
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
all
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
</tbody>
</table>

``` r
# Print how many labeled OTUs are in each sample
l2fc.df %>%
  filter(padj <= pvalue_cutoff) %>%
  group_by(.id) %>%
  summarize(n_OTUs = n()) %>%
  as.data.frame
```

    ##                                                                         .id
    ## 1 (substrate=='12C-Con' & day=='14') | (substrate=='13C-Cel' & day == '14')
    ## 2 (substrate=='12C-Con' & day=='14') | (substrate=='13C-Glu' & day == '14')
    ## 3   (substrate=='12C-Con' & day=='3') | (substrate=='13C-Cel' & day == '3')
    ## 4   (substrate=='12C-Con' & day=='3') | (substrate=='13C-Glu' & day == '3')
    ##   n_OTUs
    ## 1    116
    ## 2    298
    ## 3    201
    ## 4     59

#### 5. Save results

Don't forget to save your results.

``` r
write.table(l2fc.df, file="example_data/MWHRSIP_S2D2_output.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names=TRUE)
```

Analyses after this step are up to you and based on your own project design and question. Some example analyses can be found in [addl\_further\_analyses](addl_further_analyses.md). These examples include code for multiple treatment-control pairs like here.

Session Info
------------

``` r
sessionInfo()
```

    ## R version 3.4.4 (2018-03-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
    ## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bindrcpp_0.2.2   kableExtra_0.8.0 knitr_1.20       HTSSIP_1.4.0    
    ## [5] phyloseq_1.22.3  tibble_1.4.2     tidyr_0.8.1      dplyr_0.7.5     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] nlme_3.1-137               bitops_1.0-6              
    ##   [3] matrixStats_0.54.0         bit64_0.9-7               
    ##   [5] RColorBrewer_1.1-2         httr_1.3.1                
    ##   [7] rprojroot_1.3-2            GenomeInfoDb_1.14.0       
    ##   [9] tools_3.4.4                backports_1.1.2           
    ##  [11] R6_2.2.2                   vegan_2.5-2               
    ##  [13] rpart_4.1-13               DBI_1.0.0                 
    ##  [15] Hmisc_4.1-1                lazyeval_0.2.1            
    ##  [17] BiocGenerics_0.24.0        mgcv_1.8-24               
    ##  [19] colorspace_1.3-2           permute_0.9-4             
    ##  [21] ade4_1.7-11                nnet_7.3-12               
    ##  [23] tidyselect_0.2.4           gridExtra_2.3             
    ##  [25] DESeq2_1.18.1              bit_1.1-14                
    ##  [27] compiler_3.4.4             rvest_0.3.2               
    ##  [29] Biobase_2.38.0             htmlTable_1.12            
    ##  [31] xml2_1.2.0                 DelayedArray_0.4.1        
    ##  [33] scales_1.0.0               checkmate_1.8.5           
    ##  [35] genefilter_1.60.0          readr_1.1.1               
    ##  [37] stringr_1.3.1              digest_0.6.16             
    ##  [39] foreign_0.8-71             rmarkdown_1.10            
    ##  [41] XVector_0.18.0             base64enc_0.1-3           
    ##  [43] pkgconfig_2.0.2            htmltools_0.3.6           
    ##  [45] highr_0.7                  htmlwidgets_1.2           
    ##  [47] rlang_0.3.1                RSQLite_2.1.1             
    ##  [49] rstudioapi_0.7             bindr_0.1.1               
    ##  [51] jsonlite_1.5               BiocParallel_1.12.0       
    ##  [53] acepack_1.4.1              RCurl_1.95-4.11           
    ##  [55] magrittr_1.5               GenomeInfoDbData_1.0.0    
    ##  [57] Formula_1.2-3              biomformat_1.6.0          
    ##  [59] Matrix_1.2-14              Rcpp_1.0.1                
    ##  [61] munsell_0.5.0              S4Vectors_0.16.0          
    ##  [63] ape_5.1                    stringi_1.2.4             
    ##  [65] yaml_2.2.0                 MASS_7.3-50               
    ##  [67] SummarizedExperiment_1.8.1 zlibbioc_1.24.0           
    ##  [69] rhdf5_2.22.0               plyr_1.8.4                
    ##  [71] blob_1.1.1                 grid_3.4.4                
    ##  [73] parallel_3.4.4             lattice_0.20-35           
    ##  [75] Biostrings_2.46.0          splines_3.4.4             
    ##  [77] annotate_1.56.2            multtest_2.34.0           
    ##  [79] hms_0.4.2                  locfit_1.5-9.1            
    ##  [81] pillar_1.2.2               igraph_1.2.2              
    ##  [83] GenomicRanges_1.30.3       geneplotter_1.56.0        
    ##  [85] reshape2_1.4.3             codetools_0.2-15          
    ##  [87] stats4_3.4.4               XML_3.98-1.16             
    ##  [89] glue_1.2.0                 evaluate_0.11             
    ##  [91] latticeExtra_0.6-28        data.table_1.11.4         
    ##  [93] foreach_1.4.4              gtable_0.2.0              
    ##  [95] purrr_0.2.5                assertthat_0.2.0          
    ##  [97] ggplot2_3.1.0              xtable_1.8-3              
    ##  [99] survival_2.42-6            viridisLite_0.3.0         
    ## [101] iterators_1.0.10           memoise_1.1.0             
    ## [103] AnnotationDbi_1.40.0       IRanges_2.12.0            
    ## [105] cluster_2.0.7-1
