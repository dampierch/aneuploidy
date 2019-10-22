##
## make haplotype plots across chromosomes
##

  ## requires chrom sizes file as input; file includes cumulative length


## setup
library("ggplot2")
library("scales")
library("cowplot")

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
  #args = commandArgs(trailingOnly=TRUE)
  #file<-args[1]
setwd(crunch_dir)
file <- "TCGA-T9-A92H_cnts2R.tsv"
het_data <- read.table(file=file,header=TRUE,sep='\t')
type_levels = c('norm','tumor')
het_data$t_type <- factor(het_data$t_type,levels=type_levels,ordered=TRUE)  ## created ordered factors for tissue type

#levels(het_data$t_type) = c(paste(file,'normal'),paste(file,'tumor'))

het_data <- het_data[het_data$chrom %in% row.names(hg38.offsets),]
het_data$g_pos <- het_data$pos + hg38.offsets[as.character(het_data$chrom),]$cum_len
het_data$fract2 <- abs(het_data$fract - 0.5) + 0.5

x.min = 0
x.max = hg38.offsets[24,'len'] + hg38.offsets[24,'cum_len']

#n_samples <- nlevels(het_data$sample)

n_samples <- 1

leg.pos = c(0.025,0.92)
leg.just = c(1,0)
base_size=12

p_col <- 1
p_width <- 10
p_height <- 20

p_col <- 4
p_width <- 10
p_height <- 12
if (n_samples <= 4) {
  p_col <- 1
  p_height <- 10
  leg.just = c(0,1)
  base_size=16
}
if (n_samples <= 2) {
  p_col <- 1
  p_height <- 5.0
  base_size=16
}

type_names=c('Normal','Tumor')
s.color <- scale_color_manual(values=hue_pal()(2),labels=type_names)

#pdf(paste0("het_pdf/",out_name,"_gg_chrom2b.pdf"),width=p_width, height=p_height)

## make plot
setwd(plot_dir)
file <- "TCGA-T9-A92H_hetcnts_plot.pdf"
pdf(file=file,width=p_width,height=p_height)

theme_set(theme_linedraw(base_size=base_size))

theme1.leg <- theme(panel.border=element_rect(colour='black', size=1.0),
          panel.grid.major=element_line(colour='darkgrey', size=0.4,linetype='dashed'),
          plot.title=element_text(face='plain', hjust=0),
          axis.text.x = element_text(size=base_size, angle=45, vjust=1, hjust=0.8),
          axis.text.y = element_text(size=base_size),
          legend.position=leg.pos,
          legend.key=element_blank(),
          legend.background=element_rect(fill='white', color='black',linetype='solid',size=0.4),
          legend.justification=leg.just,
          legend.title=element_blank())

theme2.leg <- theme(panel.border=element_rect(colour='black', size=1.0),
          panel.grid.major=element_line(colour='darkgrey', size=0.4,linetype='dashed'),
          plot.title=element_text(face='plain', hjust=0),
          axis.text.x = element_text(size=10, angle=45, vjust=1, hjust=0.8),
          axis.text.y = element_text(size=base_size),
          legend.position=NULL,
          legend.key=element_blank(),
          legend.background=element_rect(fill='white', color='black',linetype='solid',size=0.4),
          legend.justification=leg.just,
          legend.title=element_blank())


x.labels = hg38.offsets$chrom

scale.x <- scale_x_continuous("chromosome",breaks=hg38.offsets$cum_len,labels=x.labels,limits=c(x.min,x.max),expand=expand_scale(c(0.01,0.01)))

p1 <- ggplot(het_data,aes(x=fract,color=t_type))+stat_density(geom='line',position='identity') +
  facet_wrap(~subject_id, ncol=1) + theme1.leg + s.color + xlab("minor allele freq.")

p2 <- ggplot(het_data,aes(x=g_pos,y=fract2))+geom_point(aes(color=t_type),size=0.1, show.legend=FALSE)+
  geom_vline(xintercept = c(hg38.offsets$cum_len,x.max),linetype='dashed') +
  facet_wrap(~t_type, ncol=1) +theme2.leg + s.color + scale.x + ylab("fraction homozygous")

# p2ab <-plot_grid(p2t,p2b,ncol=1)

plot_grid(p1,p2,ncol=2,rel_widths=c(0.33,0.66))

dev.off()
