## config
args <- commandArgs(trailingOnly=TRUE)
set_name <- args[2]
cat("set:",set_name,"\n")
source("/home/chd5n/projects/aneuploidy/scripts/config.R")

## samples
filepath <- "/scratch/chd5n/aneuploidy/raw-data/annotations/"
filename <- paste0("rna_set_",set_name,".tsv")
if (file.exists(paste0(filepath,filename))) {
  sample_set <- read_tsv(file=paste0(filepath,filename))
} else {
  cat("data not ready, run rule rank_sd first and try again.\n")
}

## gene expression
filepath <- r_dir
filename <- paste0(set_name,"_dds_sva_svA.Rdata")
if (file.exists(paste0(filepath,filename))) {
  lname <- load(file=paste0(filepath,filename))
} else {
  cat("data not ready, run rule run_de first and try again.\n")
}

## pearson correlation
expr <- counts(ddssva,normalized=TRUE)
sd_v <- setNames(sample_set$sd_tumor,sample_set$subject_id)
if (all(colnames(expr)==names(sd_v))) {
  corrs <- data.frame(ens_id=vector(),ex_sd_corr=vector())
  for (ens_id in rownames(expr)) {
    ex_v <- expr[rownames(expr)==ens_id,]
    ex_sd_corr <- cor(ex_v,sd_v,method="pearson")
    corrs <- rbind(corrs, data.frame(ens_id=ens_id,ex_sd_corr=ex_sd_corr))
  }
} else {
  cat("discordant sample names\n")
}

## rank correlations
cors_rank <- corrs[order(corrs$ex_sd_corr,decreasing=TRUE),]
cors_rank$symbol <- mapIds(org.Hs.eg.db, keys=base::substr(cors_rank$ens_id,start=1,stop=15),
                      keytype="ENSEMBL", column="SYMBOL", multiVals="first")
cors_rank$enz_id <- mapIds(org.Hs.eg.db, keys=base::substr(cors_rank$ens_id,start=1,stop=15),
                      keytype="ENSEMBL", column="ENTREZID", multiVals="first")
filepath <- table_dir
filename <- paste0(set_name,"_expr_sd_cor.tsv")
write.table(cors_rank, file=paste0(filepath,filename), sep="\t", quote=FALSE, row.names=FALSE)
cat("table save to",paste0(filepath,filename),"\n")
