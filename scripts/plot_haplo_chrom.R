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
proj_home <- "/scratch/chd5n/aneuploidy/"
data_dir <- paste0(proj_home,"raw-data/")
seq_dir <- paste0(data_dir,"sequencing/")
crunch_dir <- paste0(seq_dir,"crunch/")
plot_dir <- paste0(proj_home,"results/plots/")

## load cumulative lengths
setwd(ref_dir)
file <- "GRCh38.d1.vd1.chr1-XY.size.tsv"
hg38.offsets <- read.table(file=file,sep="\t",header=FALSE,col.names=c('chrom','len','cum_len'))
hg38.offsets$chrom <- factor(hg38.offsets$chrom,levels=unique(as.character(hg38.offsets$chrom))) ## puts chrom in order
row.names(hg38.offsets) <- levels(hg38.offsets$chrom)

## load heterozygous site allele fractions
  ## if want to do from command line, could use:
  #args = commandArgs(trailingOnly=TRUE)
  #file<-args[1]
setwd(crunch_dir)
file <- "TCGA-T9-A92H_cnts2R.tsv"
slabel <- "TCGA-T9-A92H"
het_data <- read.table(file=file,header=TRUE,sep='\t')
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
leg.just <- c(1,0)
base_size <- 12

p_col <- 4
p_width <- 10
p_height <- 12

type_names <- c('Normal','Tumor')
s.color <- scale_color_manual(values=hue_pal()(2),labels=type_names)

sub_labs <- paste(levels(good_het_data$subject_id), levels(good_het_data$t_type))
names(sub_labs) <- levels(good_het_data$t_type)

n_subjects <- nlevels(het_data$subject_id)
if (n_subjects <= 4) {
  p_col <- 1
  p_height <- 10
  leg.just <- c(0,1)
  base_size <- 16
}
if (n_subjects <= 2) {
  p_height <- 5
}

## make plot
setwd(plot_dir)
#file <- "TCGA-T9-A92H_hetcnts_plot.pdf"
file <- "sample_hetcnts_plot.pdf"
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
p1 <- ggplot(good_het_data[good_het_data$is_good,],aes(x=maj_fract,color=t_type))+stat_density(geom='line',position='identity') +
  facet_wrap(~subject_id, ncol=1) + theme1.leg + s.color + scale_x_continuous("allele fraction",limits=c(0.0,1.0),breaks=seq(0.0,1.0,0.25))

## karyotype-esque plot
p2 <- ggplot(good_het_data[good_het_data$is_good,],aes(x=pos_adj,y=allele_fract))+geom_point(aes(color=t_type),size=0.1,show.legend=FALSE) +
  geom_vline(xintercept=c(hg38.offsets$cum_len,x.max),linetype='dashed') +
  facet_wrap(~t_type, ncol=1, labeller=labeller(t_type=sub_labs)) + theme2.leg + s.color + scale.x + ylab("fraction homozygous") +
  labs(caption=slabel)

plot_grid(p1,p2,ncol=2,rel_widths=c(0.33,0.66))

dev.off()
