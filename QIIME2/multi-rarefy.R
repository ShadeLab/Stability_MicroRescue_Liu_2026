# Multi-rarefaction
# @Xipeng Liu

# Run it on a HPCC

library(readxl)
library(phyloseq)
library(BRCore)


# Load data ####

## Metadata ####

metadata <- read_excel("phyloseq/sample-metadata-dna-cdna.xlsx")

rownames(metadata) <- metadata$Sample_id


## OTU ####

# The name of first column has been manually modified as OTUID
otu <- read.table(
  "phyloseq/otu_table_dn_99.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)


## Taxonomy ####
tax_raw <- read.table("phyloseq/taxonomy.tsv", header=TRUE, sep="\t", 
                  row.names=1, fill=TRUE, stringsAsFactors=FALSE)

# Function to split and pad taxonomy
split_tax <- function(x) {
  ranks <- strsplit(as.character(x), ";")[[1]]
  # pad missing ranks with NA
  if(length(ranks) < 7) {
    ranks <- c(ranks, rep(NA, 7 - length(ranks)))
  }
  return(ranks)
}

# Apply to all rows
tax_split <- t(sapply(tax_raw$Taxon, split_tax))

# Name rows and columns
rownames(tax_split) <- rownames(tax_raw)

colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Remove QIIME2 prefixes
tax_split <- apply(tax_split, 2, function(x) gsub("^[a-z]__","", x))

tax <- tax_table(as.matrix(tax_split))


# Build phyloseq object ####

otu = otu_table(otu, taxa_are_rows = TRUE)

# Check if the names match
setdiff(sample_names(otu), rownames(metadata))

ps <- phyloseq(otu, tax, metadata)



# Multi_rarefy ####

# Perform rarefaction
otu_rare <- multi_rarefy(
  physeq_obj = ps,
  depth_level = 14714,
  num_iter = 100,
  .as = "list",
  set_seed = 2026
)

saveRDS(otu_rare, file = "phyloseq/otu_rare_depth14714_iter100.rds")

cat("Rarefaction completed! \n")



# Calculate the mean
otu_mean <- Reduce("+", otu_rare) / length(otu_rare)

saveRDS(otu_mean,"phyloseq/otu_mean_depth14714_iter100.rds")

write.table(
  otu_mean,
  "phyloseq/otu_mean_depth14714_iter100.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("Averaged OTU table saved! \n")




