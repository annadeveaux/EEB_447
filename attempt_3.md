phyloseq\_demo
================
Anna DeVeaux
10/22/2019

``` r
taxfile <- "final.FWDB.Silva.taxonomy.tsv_ed"
metafile <- "Erie-Plastics_New_Metadata.csv"
otufile <- "FinalOTUs.otu.table"


plast_tax <- read.table(taxfile, header=F, sep="\t", fill=T)
rownames(plast_tax) <- plast_tax$V1
plast_tax <- subset(plast_tax, select = -c(V1))
taxmat <- as.matrix(plast_tax)   # converts plast_tax from a list into a matrix
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
# taxmat %>% View()

plast_meta <- read.table(metafile, header=T, sep = ",")
rownames(plast_meta) <- plast_meta$Sample_ID   # why don't we subset this out after making it row names?
plast_meta_df <- as.data.frame(plast_meta)  # convert to data frame, not matrix, b/c data =/= numerical

plast_otu <- read.table(otufile, header=T, sep="\t")
View(plast_otu)
```

    ## Warning in system2("/usr/bin/otool", c("-L", shQuote(DSO)), stdout = TRUE):
    ## running command ''/usr/bin/otool' -L '/Library/Frameworks/R.framework/
    ## Resources/modules/R_de.so'' had status 1

``` r
rownames(plast_otu) <- plast_otu$OTUid
plast_otu <- subset(plast_otu, select = -c(OTUid))
otumat <- as.matrix(plast_otu)     # again, matrix
```

``` r
plast_phylo <- phyloseq(otu_table(otumat, taxa_are_rows = T), tax_table(taxmat), sample_data(plast_meta_df))

colnames(tax_table(plast_phylo)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
tax_table(plast_phylo) <- cbind(tax_table(plast_phylo), "OTU" = row.names(tax_table(plast_phylo)))
```

# Prune out samples with less than 1000 reads and then look at sequencing depth

``` r
plast_phylo_pruned <- prune_samples(sample_sums(plast_phylo) > 1000, plast_phylo)

sample_sum_pruned_df <- data.frame(sum = sample_sums(plast_phylo_pruned))

ggplot(sample_sum_pruned_df, aes(x = sum)) +
  geom_histogram(color = "black", fill = "blue", binwidth= 2500) +
  ggtitle("Distribution of sample sequencing depth (samples with <1000 reads excluded)") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](attempt_3_files/figure-gfm/prune%20and%20look%20at%20sample%20sequencing%20depth-1.png)<!-- -->
