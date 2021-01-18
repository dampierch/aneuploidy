## README
## aim: make haplotype plots across chromosomes
## method: uses chrom sizes file and parsed heterozygous site counts to create
##     density plot and karyotype-style scatter plot of major allele fractions
## notes: filters sites to exclude outlier allele fractions


library(readr)
library(ggplot2)
library(scales)
library(cowplot)
library(reshape2)


theme_default <- theme(
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour='black', size=1.0, fill=NA),
    panel.grid.major=element_line(colour='grey', size=0.4, linetype='dashed'),
    panel.grid.minor=element_line(colour="grey", size=0.2, linetype="dashed"),
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text.x=element_text(size=9, angle=45, vjust=1, hjust=0.8),
    axis.text.y=element_text(size=9),
    legend.position=c(0.025, 0.95),
    legend.key=element_blank(),
    legend.key.size=unit(0.5, "cm"),
    legend.key.width=unit(0.5, "cm"),
    legend.spacing.x=unit(0.001, "cm"),
    legend.spacing.y=unit(0.001, "cm"),
    legend.justification=c(0, 1),
    legend.title=element_blank(),
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")
)


read_chroms <- function() {
    ## load cumulative lengths
    target_dir <- "/scratch/chd5n/reference-genome/assembly/tcga/"
    target <- paste0(target_dir, "GRCh38.d1.vd1.chr1-XY.size.tsv")
    chrom_pos <- read_tsv(target, col_names=c("chrom", "len", "abs_pos0"))
    ## put chroms in order
    chrom_pos$chrom <- factor(
        chrom_pos$chrom, levels=unique(as.character(chrom_pos$chrom))
    )
    return(chrom_pos)
}


read_info <- function() {
    ## load vector of subjects to plot
    target_dir <- "/scratch/chd5n/aneuploidy/"
    target <- paste0(target_dir, "coad-read.file_info")
    df <- read_tsv(
        target, col_names=c("subject_id", "t_type", "file_name", "file_id")
    )
    subjects <- sort(unique(df$subject_id))
    return(subjects)
}


read_hetcnts <- function(subject_id, chrom_pos) {
    ## load heterozygous site maj allele fractions per subject
    target_dir <- "/scratch/chd5n/aneuploidy/hetsites-data/"
    target <- paste0(target_dir, subject_id, "_hetcnts.tsv.gz")
    data <- read_tsv(target, col_names=TRUE)
    ## factor and order tissue type
    ttype_levels <- c('norm', 'tumor')
    data$t_type <- factor(data$t_type, levels=ttype_levels, ordered=TRUE)
    ## filter out possible extraneous chroms
    data <- data[data$chrom %in% chrom_pos$chrom, ]
    return(data)
}


clean_hetcnts <- function(data, chrom_pos) {
    ## cast data into wide format
    wide_data <- reshape2::dcast(
        data,
        subject_id + chrom + pos ~ t_type,
        value.var='maj_fract'
    )
    ## set absolute genomic positions
    idx <- match(wide_data$chrom, chrom_pos$chrom)
    wide_data$pos_abs <- unlist(
        wide_data$pos + chrom_pos[idx, "abs_pos0"], use.names=FALSE
    )
    ## identify outliers in normal data and tag for exclusion
    tukey5 <- fivenum(wide_data$norm)
    r_adj <- abs(tukey5[4] - tukey5[2]) * 1.6
    r_min <- tukey5[3] - r_adj
    r_max <- tukey5[3] + r_adj
    wide_data$keep <- wide_data$norm >= r_min & wide_data$norm <= r_max
    ## melt data back into long format
    long_data <- reshape2::melt(
        wide_data,
        id.vars=c('subject_id', 'chrom', 'pos', 'pos_abs', 'keep'),
        variable.name='t_type',
        value.name='maj_fract'
    )
    long_data <- long_data[!is.na(long_data$keep), ]
    long_data <- long_data[long_data$keep, ]
    ## adjust fraction for .5 to 1 scale instead of 0 to 1
    long_data$allele_fract <- abs(long_data$maj_fract - 0.5) + 0.5
    long_data <- long_data[!is.na(long_data$allele_fract), ]
    return(long_data)
}


set_colors <- function() {
    ttype_labs <- c('Normal', 'Tumor')
    ttype_colors <- scale_color_manual(
        values=c(scales::muted("blue"), scales::muted("red")),
        labels=ttype_labs
    )
    return(ttype_colors)
}


plot_density <- function(data, ttype_colors) {
    ggp <- ggplot(data, aes(x=maj_fract, color=t_type)) +
        stat_density(geom='line', position='identity') +
        theme_default +
        theme(axis.text.x=element_text(angle=0, hjust=0.5)) +
        ttype_colors +
        scale_x_continuous(
            "Major Allele Fraction",
            limits=c(0.0, 1.0),
            breaks=seq(0.0, 1.0, 0.25)
        ) +
        ylab("Density")
    return(ggp)
}


plot_fractions <- function(data, ttype_colors, chrom_pos) {
    ## karyotype-esque plot
    strip_labs <- c(
        paste(unique(data$subject_id), "Normal",
            format(length(unique(data$pos_abs)), big.mark=","), "sites"),
        paste(unique(data$subject_id), "Tumor",
            format(length(unique(data$pos_abs)), big.mark=","), "sites")
    )
    names(strip_labs) <- levels(data$t_type)
    x_min <- 0
    x_max <- sum(
        unlist(chrom_pos[chrom_pos$chrom == "chrY", 'len'], use.names=FALSE),
        unlist(chrom_pos[chrom_pos$chrom == "chrY", 'abs_pos0'], use.names=FALSE)
    )
    x_ticks <- chrom_pos$chrom
    ggp <- ggplot(data, aes(x=pos_abs, y=allele_fract, color=t_type)) +
        geom_point(size=0.1, alpha=0.25) +
        geom_vline(
            xintercept=c(chrom_pos$abs_pos0, x_max), linetype='dashed'
        ) +
        facet_wrap(
            ~t_type, ncol=1,
            labeller=labeller(t_type=strip_labs, .multi_line=FALSE)
        ) +
        theme_default +
        theme(legend.position="none") +
        ttype_colors +
        scale_x_continuous(
            "Chromosome", breaks=chrom_pos$abs_pos0, labels=x_ticks,
            limits=c(x_min, x_max), expand=expansion(add=c(0.01, 0.01))
        ) +
        ylab("Allele Fraction") +
        scale_y_continuous(limits=c(0.5, 1.0))
    return(ggp)
}


write_plot <- function(subject_id, ggp1, ggp2) {
    target_dir <- "/scratch/chd5n/aneuploidy/plots/"
    target <- paste0(target_dir, subject_id, "_hetcnts.pdf")
    pdf(target, width=8, height=4)
    print(cowplot::plot_grid(
        ggp1, ggp2, ncol=2, rel_widths=c(0.33, 0.66), scale=c(0.95, 1)
    ))
    dev.off()
    cat("plot written to", target, "\n")
}


main <- function() {
    chrom_pos <- read_chroms()
    subjects <- read_info()
    ttype_colors <- set_colors()
    for (subject_id in subjects) {
        data <- read_hetcnts(subject_id, chrom_pos)
        data <- clean_hetcnts(data, chrom_pos)
        ggp1 <- plot_density(data, ttype_colors)
        ggp2 <- plot_fractions(data, ttype_colors, chrom_pos)
        write_plot(subject_id, ggp1, ggp2)
    }
}


main()
