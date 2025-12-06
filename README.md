# Amplicon-Data-Processing
DADA2 Amplicon Data Processing from Illumina Sequencing Technology
Amplicon data processing pipeline using DADA2

This repository documents an amplicon sequencing data processing pipeline in R using DADA2. It takes paired-end FASTQ files, performs quality filtering, denoising, merging, chimera removal, and assigns taxonomy with SILVA. Results are organized into a phyloseq object for downstream analysis.

Prerequisites and installation
- R version: 4.0 or newer
- Packages: dada2, phyloseq, ggplot2, tidyverse, BiocManager
Install packages:
install.packages(c("ggplot2", "tidyverse"))
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq"))


Repository structure
- Assignment3.Rmd: Full, runnable pipeline with explanations
- README.md: Documentation and step-by-step guide (this file)
- data/: Folder containing paired-end FASTQ files
- ref/: SILVA taxonomy training file (e.g., silva_nr99_v138.2_toSpecies_trainset.fa.gz)
- output/: Optional folder for results and figures
Tip: Your file paths in the code should point to your local data and reference files.


Quick start
- Place your FASTQ files in a folder (e.g., data/).
- Download SILVA training file and place it in ref/.
--------------------------------------------------------------
```{r}
# Load required packages
library(dada2)
library(phyloseq)
library(ggplot2)
library(tidyverse)
```

```{r}
# Set working directory to where your FASTQ files are located
setwd("C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001")
path <- ("C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001")
# List files in the directory to confirm that the FASTQ files are present
list.files(path)

```

```{r}
# Identify forward and reverse read files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE)) # Forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE)) # Reverse reads

# Assign simple sample names manually (instead of long FASTQ names)
sample.names <- c("Pumice", "NegativeControl")
names(fnFs) <- sample.names
names(fnRs) <- sample.names
```

```{r}
# Plot quality profiles to inspect read quality
plotQualityProfile(fnFs[1])   # Forward reads for pumice sample
plotQualityProfile(fnRs[1])   # Reverse reads for pumice sample

plotQualityProfile(fnFs[2])   # Forward reads for negative control
plotQualityProfile(fnRs[2])   # Reverse reads for negative control
```

```{r}
# Decide truncation lengths based on quality plots
# First number = forward read length, second number = reverse read length
truncLen <- c(230,190)
```

```{r}
# Filter and trim reads based on chosen truncation lengths
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=truncLen, maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
out   # Shows how many reads were kept after filtering
```

```{r}
# Learn error rates from filtered reads
errF <- learnErrors(filtFs, multithread=FALSE) # Forward error model
errR <- learnErrors(filtRs, multithread=FALSE) # Reverse error model

# Apply DADA algorithm to denoise reads
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```

```{r}
# Merge paired reads into full sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Build ASV table (rows = ASVs, columns = samples)
seqtab <- makeSequenceTable(mergers)

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
```

```{r}
# Assign taxonomy using SILVA reference database
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/User/OneDrive/Desktop/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
```

```{r}
# Format ASV table for phyloseq
asv_mat <- t(seqtab.nochim)                          # Transpose so rows = ASVs
rownames(asv_mat) <- paste0("ASV", 1:nrow(asv_mat))  # Rename ASVs as ASV1, ASV2, etc.

tax_mat <- as.matrix(taxa)                           # Convert taxonomy to matrix
rownames(tax_mat) <- rownames(asv_mat)               # Match ASV IDs

OTU <- otu_table(asv_mat, taxa_are_rows=TRUE)        # OTU table
TAX <- tax_table(tax_mat)                            # Taxonomy table

# Add sample metadata
sample_df <- data.frame(Sample=sample.names,
                        Type=c("Pumice","NegativeControl"),
                        row.names=sample.names)
SAM <- sample_data(sample_df)

# Rename OTU sample names to match metadata
sample_names(OTU) <- c("Pumice", "NegativeControl")

# Build phyloseq object
physeq <- phyloseq(OTU, TAX, SAM)
physeq
```
------------------------------------------------
Step-by-step workflow with code
1. Load packages
```{r}
library(dada2)
library(phyloseq)
library(ggplot2)
library(tidyverse)
```

2. Set paths and list files
```{r}
setwd("C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001")
path <- "C:/Users/User/OneDrive/Desktop/Assignment3_ASVclass-20251106T034618Z-1-001"
list.files(path)
```

3. Identify paired-end reads and name samples
```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE)) # Forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE)) # Reverse reads

sample.names <- c("Pumice", "NegativeControl")
names(fnFs) <- sample.names
names(fnRs) <- sample.names
```

4. Inspect quality profiles
```{r}
plotQualityProfile(fnFs[1]); plotQualityProfile(fnRs[1])  # Pumice
plotQualityProfile(fnFs[2]); plotQualityProfile(fnRs[2])  # Negative control
```

5. Choose truncation lengths based on quality
```{r}
truncLen <- c(230, 190)  # c(forward, reverse)
```

6. Filter and trim reads
```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=truncLen, maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
out
```

7. Learn error rates and denoise
```{r}
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```

8. Merge reads, build ASV table, remove chimeras
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=FALSE, verbose=TRUE)
```

9. Assign taxonomy with SILVA
```{r}
taxa <- assignTaxonomy(
  seqtab.nochim,
  "C:/Users/User/OneDrive/Desktop/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread=TRUE)
```

10. Create phyloseq object
```{r}
asv_mat <- t(seqtab.nochim)
rownames(asv_mat) <- paste0("ASV", 1:nrow(asv_mat))

tax_mat <- as.matrix(taxa)
rownames(tax_mat) <- rownames(asv_mat)

OTU <- otu_table(asv_mat, taxa_are_rows=TRUE)
TAX <- tax_table(tax_mat)

sample_df <- data.frame(Sample=sample.names,
                        Type=c("Pumice","NegativeControl"),
                        row.names=sample.names)
SAM <- sample_data(sample_df)

sample_names(OTU) <- c("Pumice", "NegativeControl")

physeq <- phyloseq(OTU, TAX, SAM)
physeq
```


Outputs
- ASV table: Non-chimeric ASV counts per sample
- Taxonomy table: SILVA-based taxonomic assignments
- Phyloseq object: Integrated dataset for plotting and analysis

Notes
- Adjust truncation lengths: Use your quality profiles to set truncLen appropriately.
- Keep file paths portable: Consider using relative paths like data/ and ref/ for GitHub.
- Include controls: Negative controls help monitor contamination and should be documented.

