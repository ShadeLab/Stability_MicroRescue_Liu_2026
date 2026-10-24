# 16S gene copy number correction
# @Xipeng Liu

# Run it on a HPCC

library(RasperGade16S)

rep_seqs <- "phyloseq/rep-seqs-dn-99-sequences.fasta"

pred.GCN <- predict_16SGCN_from_sequences(seqs=rep_seqs)

saveRDS(pred.GCN,"phyloseq/DNA_OTU_GCN_prediction.rds")

cat("GCN prediction completed\nSaved: phyloseq/DNA_OTU_GCN_prediction.rds\n")