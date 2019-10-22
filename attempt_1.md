AD\_Plastics\_Demo
================
Anna DeVeaux
10/21/2019

# Setup

``` r
source("~/Box Sync/EEB447_Microbes_in_the_Wild/2019/source_code_all_files/Final_PAFL_Trophicstate-master/Functions_PAFL.R", local = TRUE)
source("~/Box Sync/EEB447_Microbes_in_the_Wild/2019/source_code_all_files/miseqR.R", local=TRUE)

# Install packages
library(knitr)
library(ggplot2) 
library(cowplot)
```

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

``` r
library(vegan) 
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(scales) 
library(grid) 
library(reshape2)
library(ape)
library(ade4)
library(plyr)
library(dplyr) 
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(phyloseq) 
library(magrittr)
library(geosphere)
library(matrixStats)
```

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:reshape2':
    ## 
    ##     dcast, melt

``` r
library(DT)
library(pander)
#library(DESeq2) is for RNA-seq

# Root directory
opts_knit$set(root.dir = "~/Box Sync/EEB447_Microbes_in_the_Wild/2019/Plastics_Group/")

# Checking libraries and packages for phyloseq
#library("phyloseq")
#packageVersion("phyloseq")

# Checking libraries and packages for phyloseq
#source("https://bioconductor.org/biocLite.R") 
#biocLite("phyloseq")

#If Phyloseq is giving you trouble and can't sync up with vegan
#install.packages("devtools")
#library("devtools")
#install_github("joey711/phyloseq")

# If DESeq2 and phyloseq are needed, install like this:
#source("miseqR.R")
#source('http://bioconductor.org/biocLite.R') 
#biocLite('phyloseq') 
#installation of "DESeq2"
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
```

Set our plotting theme. This will control the formating of your figures
throughout your code. But if you want to adjust, any of these can be
overriden in your code chunk.

``` r
# Set our theme for our plots down the road
theme_set(theme_bw() + theme(plot.title = element_text(face="bold", size = 12),  #Set the plot title
                                 strip.text.x = element_text(size=10, face="bold"),  #Set Y facet titles 
                                 strip.text.y = element_text(size=10, face="bold"),  #Set X facet titles 
                                 strip.background = element_blank(),  #No facet background
                                 axis.title.x = element_text(face="bold", size=12),  #Set the x-axis title
                                 axis.title.y = element_text(face="bold", size=12),  #Set the y-axis title
                                 axis.text.x = element_text(colour = "black", size=8),  #Set the x-axis labels
                                 axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
                                 legend.title = element_text(size = 8, face="bold"),  #Set the legend title 
                                 legend.text = element_text(size = 8),  #Set the legend text
                                 legend.position="right", #Default the legend position to the right
                                 plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) #top, right, bottom, left

# This will always color samples by sample type (free-living, particle-associated, plastic) in the order listed.
#labels=c("FL","PA","Plastic")
samp_colors<-c("lightseagreen","greenyellow","orangered3")

# This will always color samples by the biome definition defined in the metadata file. 
#labels=c("Basin","Non Urban","Urban","River Plume", "WWTP distal","WWTP near")
Biome_Def4_colors<-c("blue","cyan","maroon","red","orange","pink")
```

# PART II: Data import and create phyloseq object

There are three files required to create our phyloseq object: 1. the
taxonomy file from mothur, 2. the metadata file (map file) that we
created to connect the sequence ID with the sample name and 3. the
sample data file that contains all the other important information
describing each sample (group, sample date, sample location, etc.).

## Import OTU table

``` r
# Import OTU table
OTUfile <- read.table("FinalOTUs.otu.table", header=TRUE, sep="\t")
rownames(OTUfile) <- OTUfile[,1]
OTUfile = subset(OTUfile, select = -c(OTUid) )
otumat<-as.matrix(OTUfile)
```

## Import taxonomy file

``` r
# Import taxonomy file 
TAXfile <- read.table("final.FWDB.Silva.taxonomy.tsv_ed", header=FALSE, sep="\t", fill=TRUE)
rownames(TAXfile) <- TAXfile$V1
TAXfile = subset(TAXfile, select = -c(V1) )
#head(TAXfile)
```

## Change the taxonomy names

``` r
# Change the taxonomy names
taxmat<-as.matrix(TAXfile)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
   #If you want to check on the imported taxonomy file
#head(taxmat)
#dim(taxmat)
```

## Import new sample metadata; created an additional colummn of just ‘Urban vs Non-Urban’ and named it ‘BS\_sub’

``` r
map <- read.table("Erie-Plastics_New_Metadata.csv", header=TRUE, sep=",")
rownames(map) <- map$Sample_ID
mapdf<-as.data.frame(map)
#head(map)
View(map)
```

    ## Warning in system2("/usr/bin/otool", c("-L", shQuote(DSO)), stdout = TRUE):
    ## running command ''/usr/bin/otool' -L '/Library/Frameworks/R.framework/
    ## Resources/modules/R_de.so'' had status 1

## Create the phyloseq object called “plast\_phylo”

``` r
# Create the phyloseq object called "plast_phylo"
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
MAP<-sample_data(mapdf)

plast_phylo = phyloseq(OTU, TAX, MAP) #works


#add a new rank, OTU, 
colnames(tax_table(plast_phylo)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies") 
tax_table(plast_phylo) <- cbind(tax_table(plast_phylo), "OTU"= row.names(tax_table(plast_phylo)))

#head(tax_table(plast_phylo))
```

# PART III: Cleaning up phyloseq object

## Checking sequencing depth

``` r
# Check the sequencing depth of each sample
sums_plast_phylo <- data.frame(colSums(otu_table(plast_phylo)))
colnames(sums_plast_phylo) <- "Sample_TotalSeqs"
#head(sums_plast_phylo)

# ^ Gives the sequencing depth for each sample
```

## Prune out samples with fewer than 1000 sequences

We can directly prune our phyloseq object to remove samples with \<1000
sequences. Then we create a new pruned phyloseq
object.

``` r
plast_phylo_pruned <- prune_samples(sample_sums(plast_phylo) > 1000, plast_phylo)

# Takes out samples w/ less than 1000 reads
```

# PART IV: Putting this all in action

## Research Question 1: visualizing differences in MCC between sites, substrates, chemical comp etc.

Making a stacked barplot to visualize community composition (for phyla
that make up \>2% of the sample)

``` prune

# prune out phyla below 2% in each sample
plast_phylum <- plast_phylo_pruned %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Phylum)
```

Biofilm composition for hard vs. soft polysterene samples

``` r
polyst_phylo <- plast_phylum %>%
  subset_samples(Sample_type == c("hard polysterene", "soft polysterene"))

ggplot(polyst_phylo, aes(x = Sample_type, y = Abundance, fill = Phylum)) +
  facet_grid(Location~.) +
  geom_bar(stat = "identity") +
  # label axes and add a title
```

PCoA plot between MCC from water nutrient samples vs. from
plastics

``` r
# Scale reads to even depth for free-living vs. plastic-associated microbes
plast_scale <- plast_phylo_pruned %>%
  subset_samples(Sample_type == c("water", "plastic"))
  scale_reads(round = "round")

# Use colors to show diffs in free-living vs. plastic-associated microbes and shapes for each site
# Ordinate
plast_pcoa <- ordinate(
  physeq = plast_scale,
  method = "PCoA",
  distance = "bray"
)

# Plot
plot_ordination(
  physeq = plast_scale,
  ordination = plast_pcoa,
  color = "Sample_type",
  shape = "Location",
  title = "PCoA of free-living vs. plastic-associated microbes from three freshwater sites"
)
# Can add various aes as desired
```

# Current thoughts:

Need to have separate columns in metadata for free-living vs. plastic
associated and type of plastic (if NA write NA).

The stacked bar graphs were doing a weird thing when I was trying it
with the sample data… ask?

Will we be able to analyze the abundance of certain species of microbes
(human/aquatic pathogens that we want to see if are there)?

How to use ANOVA…? (specifically for measuring significance in these
plots)
