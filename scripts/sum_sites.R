library(readr)
library(ggplot2)
library(reshape2)
library(scales)


ggp_theme_hist <- theme(
    panel.background=element_rect(fill="white"),
    panel.grid.major=element_line(color="grey", size=0.25, linetype="dashed"),
    panel.grid.minor=element_line(color="grey", size=0.25, linetype="dashed"),
    panel.border=element_rect(linetype="solid", fill=NA),
    plot.margin=margin(t=0.75, r=0.75, b=0.75, l=0.75, unit="lines"),
    plot.title=element_text(size=11, face="plain", hjust=0.5),
    plot.subtitle=element_text(size=10, face="plain", hjust=0.5),
    axis.title.x=element_text(size=10, face="plain"),
    axis.title.y=element_text(size=10, face="plain"),
    axis.text.x=element_text(size=9, face="plain"),
    axis.text.y=element_text(size=9, face="plain"),
    # axis.ticks.x=element_blank(),
    # axis.ticks.y=element_blank(),
    axis.line.x.bottom=element_blank(),
    axis.line.y.left=element_blank(),
    legend.key=element_rect(fill="white"),
    legend.position="top",
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")
)


read_summary <- function(count_type) {
    target_path <- "/scratch/chd5n/aneuploidy/"
    target <- paste0(target_path, "summary_", count_type, ".tsv")
    df <- read_tsv(target)
    return(df)
}


plot_hist <- function(dfwide, count_type, alpha=0.7) {
    ggp_title <- "Count Distributions"
    ggp_subtitle <- paste0(
        toupper(substr(count_type, 1, 1)),
        substr(count_type, 2, nchar(count_type))
    )
    ggp_xlab <- "Number Per Sample"
    ggp_ylab <- "Frequency"
    if (count_type == "sites") {binw <- 500} else {binw <- 1}
    dflong <- reshape2::melt(
        dfwide, id.vars="subject_id", variable.name="file_type", value.name="count"
    )
    dflong$t_type <- unlist(
        lapply(strsplit(as.character(dflong$file_type), "_"), "[", 2)
    )
    stat <- aggregate(data=dflong, count ~ file_type, median)
    colnames(stat) <- c("file_type", "median")
    df <- merge(dflong, stat, by="file_type", sort=FALSE)
    ggp <- ggplot(df, aes(x=count, group=file_type, fill=t_type)) +
        geom_histogram(binwidth=binw, alpha=alpha) +
        geom_vline(aes(xintercept=median), linetype="dashed", colour="black") +
        labs(title=ggp_title, subtitle=ggp_subtitle, x=ggp_xlab, y=ggp_ylab) +
        facet_wrap(vars(file_type), nrow=2, strip.position="top") +
        ggp_theme_hist +
        scale_fill_manual(
            values=c(scales::muted("blue"), scales::muted("red")),
            labels=c("Normal", "Tumor")
        )
    return(ggp)
}


write_plot <- function(ggp, count_type) {
    target_path <- "/scratch/chd5n/aneuploidy/plots/"
    target <- paste0(target_path, "summary_", count_type, ".pdf")
    ggsave(target, plot=ggp, device="pdf", width=5, height=5, units="in")
    cat("plot written to", target, "\n")
}


main <- function() {
    for (count_type in c("sites", "chroms")) {
        df <- read_summary(count_type)
        ggp <- plot_hist(df, count_type)
        write_plot(ggp, count_type)
    }
}


main()
