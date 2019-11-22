##
## make haplotype plots across chromosomes
##

  ## requires chrom sizes file as input; file includes cumulative length


## setup
library("ggplot2")
library("scales")
library("cowplot")
library("reshape2")

ref_dir <- "/scratch/chd5n/Reference_genome/"
crunch_dir <- "/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/"
plot_dir <- "/scratch/chd5n/aneuploidy/results/plots/"

args <- commandArgs(trailingOnly=TRUE)
if (any(grepl("file_set", args))) {
  set_name <- args[grep("file_set",args)]
} else {
  subjects <- args[2:length(args)]
}

## load cumulative lengths
file <- paste0(ref_dir,"GRCh38.d1.vd1.chr1-XY.size.tsv")
hg38.offsets <- read.table(file=file,sep="\t",header=FALSE,col.names=c("chrom","len","cum_len"))
hg38.offsets$chrom <- factor(hg38.offsets$chrom,levels=unique(as.character(hg38.offsets$chrom))) ## puts chrom in order
row.names(hg38.offsets) <- levels(hg38.offsets$chrom)

## load vector of subjects to iterate over
if (any(grepl("file_set", args))) {
  file <- set_name
  file_set <- read.table(file=file, sep="\t", header=FALSE, col.names=c("subject_id","normal","tumor"))
  subjects <- as.character(file_set$subject_id)
}

## load heterozygous site allele fractions for each subject
file <- list()
for (sub in subjects) {
  file[[sub]] <- paste0(crunch_dir,sub,"_cnts2R.tsv")
}
het_data <- read.table(file=file[[1]],header=TRUE,sep='\t')
if (length(file) > 1) {
  for (i in 2:length(file)) {
    het_data <- rbind(het_data, read.table(file=file[[i]],header=TRUE,sep='\t'))
  }
}
type_levels <- c('norm','tumor')
het_data$t_type <- factor(het_data$t_type,levels=type_levels,ordered=TRUE)  ## created ordered factors for tissue type
het_data <- het_data[het_data$chrom %in% row.names(hg38.offsets),] ## filter out possible extraneous chrom

## cast data into wide format
wide_het_data <- dcast(het_data, subject_id+chrom+pos ~ t_type, value.var='maj_fract')
wide_het_data$pos_adj <- wide_het_data$pos + hg38.offsets[as.character(wide_het_data$chrom),][,"cum_len"] ## adjust positions for full karyotype

## identify outliers in normal data and tag for exclusion
tukey5 <- fivenum(wide_het_data$norm)
iqr_adj <- abs(tukey5[4] - tukey5[2]) * 1.6
iqr_min <- tukey5[3] - iqr_adj
iqr_max <- tukey5[3] + iqr_adj
wide_het_data$is_good <- wide_het_data$norm >= iqr_min & wide_het_data$norm <= iqr_max

## melt data back into long format
good_het_data <- melt(wide_het_data, id.vars=c('subject_id','chrom','pos','pos_adj','is_good'), variable.name='t_type', value.name='maj_fract')
good_het_data <- good_het_data[!is.na(good_het_data$is_good),]
good_het_data$allele_fract <- abs(good_het_data$maj_fract - 0.5) + 0.5  ## adjust fraction for .5 to 1 scale instead of 0 to 1
good_het_data <- good_het_data[!is.na(good_het_data$allele_fract),]

## set plotting parameters
x.min <- 0
x.max <- hg38.offsets[24,'len'] + hg38.offsets[24,'cum_len']
x.labels <- hg38.offsets$chrom
scale.x <- scale_x_continuous("chromosome",breaks=hg38.offsets$cum_len,labels=x.labels,limits=c(x.min,x.max),expand=expand_scale(c(0.01,0.01)))

leg.pos <- c(0.025,0.92)
leg.just <- c(0,1)
base_size <- 12

p_col <- 4
p_width <- 10
p_height <- 20

type_names <- c('Normal','Tumor')
s.color <- scale_color_manual(values=rev(hue_pal()(2)),labels=type_names)
sub_labs <- rep(type_names, nlevels(good_het_data$subject_id))
names(sub_labs) <- rep(levels(good_het_data$t_type), nlevels(good_het_data$subject_id))

if (nlevels(het_data$subject_id) <= 4) {
  p_col <- 1
  p_height <- 10
  base_size <- 16
}
if (nlevels(het_data$subject_id) <= 2) {
  p_height <- 5
}
if (nlevels(het_data$subject_id) > 9) {
  p_height <- 40
  leg.pos <- c(0.025,0.98)
  leg.just <- c(0,0)
}

## make plot
setwd(plot_dir)
if (nlevels(good_het_data$subject_id) < 2) {
  file <- paste0(levels(good_het_data$subject_id),"_hetcnts_plot.pdf")
} else {
  fn <- unlist(strsplit(unlist(strsplit(set_name,"\\."))[1],"/"))[length(unlist(strsplit(unlist(strsplit(set_name,"\\."))[1],"/")))]
  file <- paste0(fn,"_hetcnts_plot.pdf")
}
pdf(file=file,width=p_width,height=p_height)

theme_set(theme_linedraw(base_size=base_size))

theme1.leg <- theme(panel.border=element_rect(colour='black', size=1.0),
          panel.grid.major=element_line(colour='darkgrey', size=0.4, linetype='dashed'),
          plot.title=element_text(face='plain', hjust=0),
          axis.text.x=element_text(size=base_size, angle=45, vjust=1, hjust=0.8),
          axis.text.y=element_text(size=base_size),
          legend.position=leg.pos,
          legend.key=element_blank(),
          legend.background=element_rect(fill='white', color='black', linetype='solid', size=0.4),
          legend.justification=leg.just,
          legend.title=element_blank())

theme2.leg <- theme(panel.border=element_rect(colour='black', size=1.0),
          panel.grid.major=element_line(colour='darkgrey', size=0.4, linetype='dashed'),
          plot.title=element_text(face='plain', hjust=0),
          axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=0.8),
          axis.text.y=element_text(size=base_size),
          legend.position=NULL,
          legend.key=element_blank(),
          legend.background=element_rect(fill='white', color='black', linetype='solid', size=0.4),
          legend.justification=leg.just,
          legend.title=element_blank())

## density plot
p1 <- ggplot(good_het_data[good_het_data$is_good,],aes(x=maj_fract,color=t_type)) +
  stat_density(geom='line',position='identity') +
  facet_wrap(~subject_id, ncol=1) + theme1.leg + s.color +
  scale_x_continuous("allele fraction",limits=c(0.0,1.0),breaks=seq(0.0,1.0,0.25))

## karyotype-esque plot
p2 <- ggplot(good_het_data[good_het_data$is_good,],aes(x=pos_adj,y=allele_fract)) +
  geom_point(aes(color=t_type),size=0.1,show.legend=FALSE) +
  geom_vline(xintercept=c(hg38.offsets$cum_len,x.max),linetype='dashed') +
  facet_wrap(subject_id~t_type, ncol=1, labeller=labeller(t_type=sub_labs, .multi_line=FALSE)) +
  theme2.leg + s.color + scale.x + ylab("fraction homozygous")

plot_grid(p1,p2,ncol=2,rel_widths=c(0.33,0.66))

dev.off()
