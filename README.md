# Amplicon-Data-Processing
DADA2 Amplicon Data Processing from Illumina Sequencing Technology
Amplicon data processing pipeline using DADA2

This repository documents an amplicon sequencing data processing pipeline in R using DADA2. It takes paired-end FASTQ files, performs quality filtering, denoising, merging, chimera removal, and assigns taxonomy with SILVA. Results are organized into a phyloseq object for downstream analysis.

**Prerequisites and installation**
- R version: 4.0 or newer
- Packages: dada2, phyloseq, ggplot2, tidyverse, BiocManager
- SILVA database: www.arb-silva.de

**Install packages if needed:**
```{r}
install.packages(c("ggplot2", "tidyverse"))
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq"))
```

**Repository structure**
- README.md: Documentation and step-by-step guide (this file)
- data/: Folder containing paired-end FASTQ files
- ref/: SILVA taxonomy training file (e.g., silva_nr99_v138.2_toSpecies_trainset.fa.gz)
- output/: Optional folder for results and figures
Tip: Your file paths in the code should point to your local data and reference files.

**Quick start**
- Place your FASTQ files in a folder (e.g., data/).
- Download SILVA training file and place it in ref/.

**Step-by-step workflow with code**
1. Load required packages
We begin by loading all the packages needed for amplicon data processing, visualization, and data handling.
```{r}
library(dada2)
library(phyloseq)
library(ggplot2)
library(tidyverse)
```

2. Set working directory and confirm FASTQ files
Set the working directory to where your FASTQ files are located. Then list the files to confirm they are present
```{r}
setwd("C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001")
path <- "C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001"
list.files(path)
```

3. Identify forward and reverse read files
Forward reads end with _R1_001.fastq and reverse reads with _R2_001.fastq. We also assign simple sample names for clarity.
```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE)) # Forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE)) # Reverse reads

sample.names <- c("Pumice", "NegativeControl")
names(fnFs) <- sample.names
names(fnRs) <- sample.names
```

4. Plot quality profiles
Inspect the quality of reads to decide truncation lengths. This step is critical for determining where to trim sequences.
```{r}
plotQualityProfile(fnFs[1])   # Forward reads for pumice sample
plotQualityProfile(fnRs[1])   # Reverse reads for pumice sample

plotQualityProfile(fnFs[2])   # Forward reads for negative control
plotQualityProfile(fnRs[2])   # Reverse reads for negative control
```

5. Decide truncation lengths
Based on the quality plots, choose truncation lengths for forward and reverse reads.
```{r}
truncLen <- c(230,190)  # First number = forward read length, second number = reverse read length
```

6. Filter and trim reads
Filter and trim reads using chosen truncation lengths. This step removes low-quality bases and ensures clean data for downstream analysis.
```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=truncLen, maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
out
```

7. Learn error rates and denoise
DADA2 learns error rates from the filtered reads and applies its algorithm to denoise sequences.
```{r}
errF <- learnErrors(filtFs, multithread=FALSE) # Forward error model
errR <- learnErrors(filtRs, multithread=FALSE) # Reverse error model

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```

8. Merge paired reads and build ASV table
Merge forward and reverse reads into full sequences, then build an ASV table. Finally, remove chimeric sequences.
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=FALSE, verbose=TRUE)
```

9. Assign taxonomy
Use the SILVA reference database to assign taxonomy to the ASVs.
```{r}
taxa <- assignTaxonomy(
  seqtab.nochim,
  "C:/Users/User/OneDrive/Desktop/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread=TRUE)
```

10. Format ASV table for phyloseq
Prepare ASV and taxonomy matrices, add sample metadata, and build a phyloseq object for downstream analysis.
```{r}
asv_mat <- t(seqtab.nochim)                          # Transpose so rows = ASVs
rownames(asv_mat) <- paste0("ASV", 1:nrow(asv_mat))  # Rename ASVs as ASV1, ASV2, etc.

tax_mat <- as.matrix(taxa)                           # Convert taxonomy to matrix
rownames(tax_mat) <- rownames(asv_mat)               # Match ASV IDs

OTU <- otu_table(asv_mat, taxa_are_rows=TRUE)        # OTU table
TAX <- tax_table(tax_mat)                            # Taxonomy table

sample_df <- data.frame(Sample=sample.names,
                        Type=c("Pumice","NegativeControl"),
                        row.names=sample.names)
SAM <- sample_data(sample_df)

sample_names(OTU) <- c("Pumice", "NegativeControl")

physeq <- phyloseq(OTU, TAX, SAM)
physeq
```

**Outputs**
- ASV table: Non-chimeric ASV counts per sample
- Taxonomy table: SILVA-based taxonomic assignments
- Phyloseq object: Integrated dataset for plotting and analysis

**Notes**
- Adjust truncation lengths: Use your quality profiles to set truncLen appropriately.
- Keep file paths portable: Consider using relative paths like data/ and ref/ for GitHub.
- Include controls: Negative controls help monitor contamination and should be documented.
