library("reshape2")
library("ggplot2")
library("grid")


ggp_theme_hist <- theme(
  panel.background=element_rect(fill="white", colour="black", size=1, linetype="solid"),
  panel.grid.major=element_line(color="white"),
  panel.grid.minor=element_line(color="white"),
  plot.margin=margin(t=0.25,r=0.5,b=0.25,l=0.25,unit="cm"),
  plot.title=element_text(size=12,face="bold",hjust=0.5),
  plot.subtitle=element_text(size=10,face="bold",hjust=0.5),
  axis.title.x=element_text(size=10, face="plain"),
  axis.title.y=element_text(size=10, face="plain"),
  axis.text.x=element_text(size=9, face="plain", angle=0, hjust=0.5),
  axis.text.y=element_text(size=9, face="plain"),
  legend.key=element_rect(fill="white"),
  legend.position="none",
  legend.title=element_blank(),
  strip.background=element_rect(fill="black"),
  strip.text=element_text(colour="white")
)

ggp_theme_point <- theme(
  panel.background=element_rect(fill="white"),
  panel.grid.major=element_line(color="white"),
  panel.grid.minor=element_line(color="white"),
  plot.background=element_rect(fill="white"),
  plot.margin=margin(t=1,r=1,b=1,l=1,unit="lines"),
  plot.title=element_text(size=12,face="bold",hjust=0.5),
  plot.subtitle=element_text(size=10,face="bold",hjust=0.5),
  axis.title.y=element_text(size=10, face="plain"),
  axis.text.x=element_blank(),
  axis.text.y=element_text(size=9, face="plain"),
  axis.ticks.x=element_blank(),
  axis.line.x.bottom=element_line(),
  axis.line.y.left=element_line(),
  legend.key=element_rect(fill="white"),
  legend.position="top",
  legend.title=element_blank(),
  strip.background=element_rect(fill="black"),
  strip.text=element_text(colour="white")
)


get_data_paths <- function(data_dir) {
  data_paths <- dir(data_dir)
  idx <- grepl("_cnts2R.tsv",data_paths)
  data_paths <- data_paths[idx]
  data_paths <- paste0(data_dir,data_paths)
  return(data_paths)
}


rank_tumors <- function(data_paths) {
  sd_samples <- data.frame()
  for (path in data_paths) {
    data <- read.table(file=path,header=TRUE,sep='\t')
    if(nrow(data)==0) {next}
    wide_data <- dcast(data, subject_id+chrom+pos ~ t_type, value.var='maj_fract')
    subject_id <- unique(as.character(data$subject_id))
    sd_tumor <- sd(wide_data$tumor, na.rm=TRUE)
    sd_norm <- sd(wide_data$norm, na.rm=TRUE)
    sd_samples <- rbind(sd_samples, data.frame(subject_id=subject_id,sd_tumor=sd_tumor,sd_norm=sd_norm))
  }
  sd_samples <- sd_samples[order(sd_samples$sd_tumor,decreasing=TRUE),]
  sd_samples$rank <- 1:nrow(sd_samples)
  return(sd_samples)
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
  samples_ranked <- rank_tumors(data_paths)
  anno_dir <- "/scratch/chd5n/aneuploidy/raw-data/annotations/"
  pheno_data <- get_pheno_data(anno_dir)
  samples_sd_rank <- merge(samples_ranked, pheno_data, by="subject_id", all.x=TRUE, all.y=FALSE, sort=FALSE)
  filename <- "/scratch/chd5n/aneuploidy/results/rdata/samples_sd_rank.Rdata"
  save(samples_sd_rank,file=filename)
  cat("data saved to",filename,"\n")
  filename <- "/scratch/chd5n/aneuploidy/results/tables/samples_sd_rank.tsv"
  write.table(samples_sd_rank,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  cat("table saved to",filename,"\n")
  return(samples_sd_rank)
}


plot_histogram <- function(samples_sd_rank) {
  ## make dataframe
  df1 <- data.frame(subject_id=samples_sd_rank$subject_id,
                    sd=samples_sd_rank$sd_tumor,
                    tissue_type="Tumor"
                  )
  df2 <- data.frame(subject_id=samples_sd_rank$subject_id,
                    sd=samples_sd_rank$sd_norm,
                    tissue_type="Healthy"
                  )
  df <- rbind(df1,df2)
  ## set labels
  ggp_title <- "Distribution of colorectal aneuploidy scores"
  ggp_xlab <- "Aneuploidy score"
  ggp_ylab <- "Number samples"
  ## make plot object
  ggp <- ggplot(df, aes(x=sd, fill=tissue_type)) +
          geom_histogram(
            stat="bin", position="stack",
            binwidth=0.01, alpha=0.5
          ) +
          labs(
            title=ggp_title,
            x=ggp_xlab, y=ggp_ylab
          ) +
          facet_wrap(
            ~tissue_type,
            nrow=2, ncol=1,
            scales="free_y"
          ) +
          scale_fill_manual(values=c("red","blue")) +
          ggp_theme_hist +
          xlim(NA,0.405)
  ## save plot
  filepath <- "/scratch/chd5n/aneuploidy/results/plots/"
  filename <- "sd_hist.pdf"
  ggsave(paste0(filepath,filename),plot=ggp,device="pdf",height=5,width=5,unit="in")
  cat("plot saved to",paste0(filepath,filename),"\n")
}


plot_ranks <- function(samples_sd_rank) {
  ## make dataframe
  keep <- samples_sd_rank$msi_status!="" & samples_sd_rank$msi_status!="Indeterminate"
  plot_df <- samples_sd_rank[keep,]
  plot_df$rank2 <- nrow(plot_df):1
  ## make plot object
  ggp <- ggplot(plot_df, aes(x=rank2, y=sd_tumor, color=factor(msi_status))) +
            geom_point(aes(y=sd_norm), colour="grey", size=0.25) +
            geom_point(size=1) +
            geom_point(aes(y=max(sd_tumor)*1.1), size=5, shape=73) +
            xlab("Aneuploidy rank") +
            ylab("Aneuploidy score") +
            labs(color="MSI status") +
            ggtitle(paste("Aneuploidy scores across",nrow(plot_df),"tumors")) +
            annotate("text", x=min(plot_df$rank2)+50, y=min(plot_df$sd_tumor)*0.9, label="less aneuploid", size=3, fontface="italic") +
            annotate("text", x=max(plot_df$rank2)-50, y=min(plot_df$sd_tumor)*0.9, label="more aneuploid", size=3, fontface="italic") +
            ggp_theme_point
  ## save plot
  filename <- "/scratch/chd5n/aneuploidy/results/plots/sd_rank.pdf"
  ggsave(filename, plot=ggp, device="pdf", path=NULL,
    scale=1, width=6, height=4, units=c("in"), dpi=300, limitsize=TRUE)
  cat("plot saved to",filename,"\n")
  ## write data to table
  samples_focus <- plot_df
  filename <- "/scratch/chd5n/aneuploidy/results/tables/samples_sd_rank_focus.tsv"
  write.table(samples_focus,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  cat("table saved to",filename,"\n")
  return(samples_focus)
}


prep_rna <- function(samples_focus,set_name) {
  if (set_name == "simple") {
    keep <- samples_focus$file_id_rna!=""  ## select simply tumors with RNA-seq data, regardless of MSI status
  } else if (set_name == "mss") {
    keep <- samples_focus$file_id_rna!="" & samples_focus$msi_status=="MSS"  ## select tumors with RNA-seq data and MSI-status=MSS
  } else if (set_name == "mss_msil") {
    keep <- samples_focus$file_id_rna!="" & (samples_focus$msi_status=="MSS" | samples_focus$msi_status=="MSI-L")  ## select tumors with RNA-seq data and MSI-status=MSS/MSI-L
  }
  samples_focus <- samples_focus[keep,]
  target_num <- round(nrow(samples_focus)*(1/4))
  upper <- samples_focus[1:target_num,]
  lower <- samples_focus[(nrow(samples_focus)-target_num+1):nrow(samples_focus),]
  upper$category <- "upper"
  lower$category <- "lower"
  sample_set <- rbind(upper,lower)
  filename <- paste0("/scratch/chd5n/aneuploidy/raw-data/annotations/rna_set_",set_name,".tsv")
  write.table(sample_set,file=filename,sep="\t",quote=FALSE,row.names=FALSE)
  cat("table saved to",filename,"\n")
  return(sample_set)
}


main <- function() {
  filepath <- "/scratch/chd5n/aneuploidy/results/rdata/"
  filename <- "samples_sd_rank.Rdata"
  if (file.exists(paste0(filepath,filename))) {
    lname <- load(paste0(filepath,filename))
  } else {
    samples_sd_rank <- get_rank_results()
  }
  plot_histogram(samples_sd_rank)
  samples_focus <- plot_ranks(samples_sd_rank)
  simple_set <- prep_rna(samples_focus,"simple")
  mss_set <- prep_rna(samples_focus,"mss")
  mss_msil_set <- prep_rna(samples_focus,"mss_msil")
}


main()
