library("reshape2")
library("ggplot2")
library("grid")


get_data_paths <- function(data_dir) {
  data_paths <- dir(data_dir)
  idx <- grepl("_cnts2R.tsv",data_paths)
  data_paths <- data_paths[idx]
  data_paths <- paste0(data_dir,data_paths)
  return(data_paths)
}


rank_tumors <- function(data_paths) {
  sd_tumors <- data.frame()
  for (path in data_paths) {
    data <- read.table(file=path,header=TRUE,sep='\t')
    if(nrow(data)==0) {next}
    wide_data <- dcast(data, subject_id+chrom+pos ~ t_type, value.var='maj_fract')
    subject_id <- unique(as.character(data$subject_id))
    sd_tumor <- sd(wide_data$tumor, na.rm=TRUE)
    sd_tumors <- rbind(sd_tumors, data.frame(subject_id=subject_id,sd_tumor=sd_tumor))
  }
  sd_tumors <- sd_tumors[order(sd_tumors$sd_tumor,decreasing=TRUE),]
  sd_tumors$rank <- 1:nrow(sd_tumors)
  return(sd_tumors)
}


get_pheno_data <- function(anno_dir) {
  p1_path <- paste0(anno_dir,"pheno.tsv")
  p1_data <- read.delim(p1_path)
  p2_path <- paste0(anno_dir,"coad-read.file_info")
  p2_colnames <- c("subject_id","tissue_type","file_name","uuid")
  p2_data <- read.delim(p2_path,header=FALSE,col.names=p2_colnames)
  p2_data <- p2_data[p2_data$tissue_type == "tumor",]
  keep <- p1_data$file_id %in% p2_data$uuid
  pheno_data <- p1_data[keep,]
  return(pheno_data)
}


get_rank_results <- function() {
  data_dir <- "/scratch/chd5n/aneuploidy/hetsites-data/r-cnts/"
  data_paths <- get_data_paths(data_dir)
  tumors_ranked <- rank_tumors(data_paths)
  anno_dir <- "/scratch/chd5n/aneuploidy/raw-data/annotations/"
  pheno_data <- get_pheno_data(anno_dir)
  tumors_sd_rank <- merge(tumors_ranked, pheno_data, by="subject_id", all.x=TRUE, all.y=FALSE, sort=FALSE)
  filename <- "/scratch/chd5n/aneuploidy/results/tables/tumors_sd_rank.Rdata"
  save(tumors_sd_rank,file=filename)
  filename <- "/scratch/chd5n/aneuploidy/results/tables/tumors_sd_rank.tsv"
  write.table(tumors_sd_rank,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  return(tumors_sd_rank)
}


plot_ranks <- function(tumors_sd_rank) {
  keep <- tumors_sd_rank$msi_status!="" & tumors_sd_rank$msi_status!="Indeterminate"
  plot_df <- tumors_sd_rank[keep,]
  plot_df$rank2 <- nrow(plot_df):1
  pp <- ggplot(plot_df, aes(x=rank2, y=sd_tumor, color=factor(msi_status))) +
    geom_point(size=1) +
    geom_point(aes(y=max(sd_tumor)*1.1), size=5, shape=73) +
    xlab("Aneuploidy rank") +
    ylab("Aneuploidy score") +
    labs(color="MSI status") +
    ggtitle(paste("Aneuploidy scores across",nrow(plot_df),"tumors")) +
    annotate("text", x=min(plot_df$rank2)+50, y=min(plot_df$sd_tumor)*0.9, label="less aneuploid", size=3, fontface="italic") +
    annotate("text", x=max(plot_df$rank2)-50, y=min(plot_df$sd_tumor)*0.9, label="more aneuploid", size=3, fontface="italic") +
    theme(
      plot.title=element_text(hjust=0.5),
      plot.margin=unit(c(1,1,2,1), "lines"),
      legend.position="top",
      plot.background=element_rect(fill="white"),
      panel.background=element_rect(fill="white"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.x.bottom=element_line(),
      axis.line.y.left=element_line()
    )
  filename <- "/scratch/chd5n/aneuploidy/results/plots/tumors_sd_rank.pdf"
  ggsave(filename, plot=pp, device="pdf", path=NULL,
    scale=1, width=9, height=6, units=c("in"), dpi=300, limitsize=TRUE)
  tumors_focus <- plot_df
  filename <- "/scratch/chd5n/aneuploidy/results/tables/tumors_sd_rank_focus.Rdata"
  save(tumors_focus,file=filename)
  filename <- "/scratch/chd5n/aneuploidy/results/tables/tumors_sd_rank_focus.tsv"
  write.table(tumors_focus,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  return(tumors_focus)
}


prep_rna_dge_simple <- function(tumors_focus) {
  keep <- tumors_focus$file_id_rna!=""
  tumors_focus <- tumors_focus[keep,]
  target_num <- round(nrow(tumors_focus)*(1/4))
  upper <- tumors_focus[1:target_num,]
  lower <- tumors_focus[(nrow(tumors_focus)-target_num+1):nrow(tumors_focus),]
  upper$category <- "upper"
  lower$category <- "lower"
  simple_set <- rbind(upper,lower)
  filename <- "/scratch/chd5n/aneuploidy/raw-data/annotations/rna_set_simple.tsv"
  write.table(simple_set,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  return(simple_set)
}


main <- function() {
  filename <- "/scratch/chd5n/aneuploidy/results/tables/tumors_sd_rank.Rdata"
  if (file.exists(filename)) {
    lname <- load(filename)
  } else {
    tumors_sd_rank <- get_rank_results()
  }
  tumors_focus <- plot_ranks(tumors_sd_rank)
  simple_set <- prep_rna_dge_simple(tumors_focus)
}


main()
