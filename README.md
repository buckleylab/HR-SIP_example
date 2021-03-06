# Multiple window high resolution DNA-stable isotope probing tutorials

DNA-stable isotope probing (DNA-SIP) is a method that allows us to link microbes to their function in complex environmental samples such as soil. There are many articles discussing the use of DNA-SIP
for microbial ecology studies in various environments. The tutorials found here go through some basic analyses for running multiple window high resolution DNA-SIP (MW-HR-SIP) on 16S rRNA gene amplicon sequence data.
This method allow us to identify bacterial operational taxonomic units (OTUs) that are significantly labeled with <sup>13</sup>Carbon or <sup>15</sup>Nitrogen. MW-HR-SIP and an earlier method with the same principle (high resolution DNA-SIP) are discussed more in the following papers:

* Pepe-Ranney C, Campbell AN, Koechli CN, Berthrong S and Buckley DH (2016) Unearthing the Ecology of Soil Microorganisms Using a High Resolution DNA-SIP Approach to Explore Cellulose and Xylose Metabolism in Soil. Front. Microbiol. 7:703. doi:[10.3389/fmicb.2016.00703](https://doi.org/10.3389/fmicb.2016.00703)
* Youngblut ND, Barnett SE and Buckley DH (2018) SIPSim: A Modeling Toolkit to Predict Accuracy and Aid Design of DNA-SIP Experiments. Front. Microbiol. 9:570. doi:[10.3389/fmicb.2018.00570](https://doi.org/10.3389/fmicb.2018.00570)

MW-HR-SIP can be performed in R using functions from the package [HTSSIP](https://cran.r-project.org/web/packages/HTSSIP/index.html)
available on CRAN. These tutorials will take you through some preliminary analyses, the main functions, and some basic follow-up analyses using this package and some custom code for interesting figures.

For more information on the HTSSIP package, please refer to:

Youngblut ND, Barnett SE, Buckley DH (2018) HTSSIP: An R package for analysis of high throughput sequencing data from nucleic acid stable isotope probing (SIP) experiments. PLoS ONE 13(1): e0189616. doi:[10.1371/journal.pone.0189616](https://doi.org/10.1371/journal.pone.0189616) 

You can also find more example analyses with HTSSIP in the [vignettes](https://cran.r-project.org/web/packages/HTSSIP/vignettes/HTSSIP_intro.html) on CRAN.

The [basic tutorial](Chapter_Examples.md) provided here was developed for the [MW-HR-SIP chapter](https://doi.org/10.1007/978-1-4939-9721-3_9) in the Methods in Molecular Biology book [Stable Isotope Probing](https://doi.org/10.1007/978-1-4939-9721-3):

Barnett SE, Youngblut ND, Buckley DH (2019) Data Analysis for DNA Stable Isotope Probing Experiments Using Multiple Window High-Resolution SIP. In: Dumont M., Hernández García M. (eds) Stable Isotope Probing. Methods in Molecular Biology, vol 2046. Humana, New York, NY. doi:[10.1007/978-1-4939-9721-3_9](https://doi.org/10.1007/978-1-4939-9721-3_9)

The code found in this chapter is provided in this tutorial and expanded upon.

## The tutorials

All tutorials are formatted as markdown files compatible with GitHub (.md). They were writen in Rmarkdown and original files are also included here (.Rmd).

### [Simple HTSSIP example](Chapter_Examples.md)
Runs through a simple example for analyzing a real amplicon dataset with MW-HR-SIP.
The dataset used here includes a single treatment and control pair. The code presented here is also described in the [MW-HR-SIP chapter](https://doi.org/10.1007/978-1-4939-9721-3_9) of [Stable Isotope Probing from Methods](https://doi.org/10.1007/978-1-4939-9721-3_9) in Methods in Molecular Biology.

### [Multi-sample HTSSIP example](HRSIP_multiple_samples.md)
A more complex example for analyzing an amplicon dataset that includes multiple treatments and sampling days with their corresponding controls.

### [Additional preliminary analyses](addl_prelim_analyses.md)
Examples of additional analyses that can be run prior to MW-HR-SIP:

* Beta-diversity between treatments and controls across fractions.
* Estimating buoyant density shift of the community.

### [Additional follow-up analyses](addl_further_analyses.md)
Examples of additional analyses and figures that can be done after running MW-HR-SIP:

* Examine taxonomy of labeled OTUs.
* Examine phylogeny of labeled OTUs.
 
## The data

Example data can be found in [example_data](example_data/).
All datasets are in phyloseq format. For more information on this format please refer to [this site](https://joey711.github.io/phyloseq/).

These datasets come from experiments adding <sup>13</sup>C-labeled substrates 
(<sup>13</sup>C-Amino acids, <sup>13</sup>C-Glucose, <sup>13</sup>C-Cellulose) 
to soil in microcosms. There were also some microcosms where unlabeled (<sup>12</sup>C)
substrates were added as control samples. The soils were harvested on various days. DNA was extracted
from the soil samples and either sequenced right away (unfractionated data) or added to a CsCl gradient,
run on the ultracentrifuge, and fractionated. Then the V4 region of the 16S rRNA gene was sequenced
from each fraction. All sequencing was performed with an Illumina MiSeq and processed using a standard sequence processing pipeline. 
The final data packaged into the phyloseq objects are the OTU table, taxonomy table, phylogenetic tree, and sample metadata table.

* **unfractionated_phyloseq.rds:** Includes processed sequencing data from the unfractionated samples (*i.e.* samples from the microcosms that were not fractionated but directly sequenced). These samples correspond to the fractionated samples in SIP_phyloseq.rds.  
* **SIP_phyloseq.rds:** The processed sequencing data from fractions of a single treatment and control microcosm pair.
* **example_S2D2_phyloseq.rds:** The processed sequencing data from fractions from treatments and controls harvested on two sampling days.

## Authors
Tutorials were writen by Samuel Barnett with input from Nicholas Youngblut
and Daniel Buckley. Some aspects of the tutorials were based on the HTSSIP
vignettes written by Nicholas Youngblut.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

