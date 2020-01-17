###
### expression analysis count prep
###


## README

  ## samples section reads pheno file from tsv into R
  ## files section sets paths to count files; paths generated
  ## from information in pheno list and hard-coded text


## parse command line
args <- commandArgs(trailingOnly=TRUE)
set_name <- args[1]
pn <- args[2]
refine <- args[3]
choose_nsv <- args[4]
cat(paste("Pass #",pn,"\n"))
cat(paste("Refine setting:",refine,"\n"))
if (refine>0) {cat(paste("Num sv chosen:",choose_nsv,"\n"))}


## source config
source("/home/chd5n/projects/aneuploidy/scripts/config.R")
numWorkers <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
cat(paste("Number workers/cores/cpus-per-task:", numWorkers, "\n"))


## set functions
get_dds <- function(set_name,count_dir,r_dir) {
  ## loads or creates deseq2 dataset from htseq counts
  cat(paste("getting dds for set",set_name),"\n\n")
  start_time <- Sys.time()
  setwd(r_dir)
  filename <- paste0(set_name,"_dds.Rdata")
  if (file.exists(filename)) {
    cat(filename,"already exists\nloading prior...\n\n")
    lname <- load(filename)
  } else {
    cat(filename,"not found\nbuilding dds now...\n\n")
    anno_fn <- paste0(anno_dir, paste0("rna_set_",set_name,".tsv"))
    anno <- read_tsv(anno_fn)
    file_paths <- data.frame(prefix=rep(count_dir,nrow(anno)), uuid=anno$file_id_rna, fn=substr(anno$file_name_rna,start=1,stop=nchar(anno$file_name_rna)-3))
    file_paths$full <- paste0(file_paths$prefix,file_paths$uuid,"/",file_paths$fn)
    sample_info <- data.frame(sampleName=anno$subject_id, fileName=file_paths$full)
    sample_table <- cbind(sample_info,data.frame(anno))
    dds <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, directory="", design=~1)
    save(dds, file=filename)
  }
  cat("done getting dds\n\n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dds)
}


pca_single <- function(vsd, pcaData, percentVar, aes_name, set_name, plab) {
  ## single frame pca plot
  color <- sym(aes_name)
  if (plab=="null_lab") {sz=2} else {sz=1}
  percentVar <- percentVar
  base_name <- set_name
  main <- paste("PCA plot for",base_name,"by",aes_name)
  color_values <- group_colors[[aes_name]]
  ggplot( pcaData, aes(PC1, PC2, color=!!color) ) +
    geom_point(size=sz) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + theme_classic() +
    theme(plot.title=element_text(hjust=0.5)) +
    ggtitle(main) +
    geom_text(aes(label=!!sym(plab)), hjust=0, vjust=0, size=3, show.legend=FALSE) #+
    if (class(colData(vsd)[,aes_name])=="character") {
      scale_color_manual(values=color_values)
    } else {
      scale_colour_gradient(low=color_values[175], high=color_values[25])
    }
}


deseq_explore <- function(dds,set_name,r_dir,plot_dir) {
  cat(paste("exploring counts with DESeq2 for set",set_name),"\n\n")
  start_time <- Sys.time()
  ## prep data
  setwd(r_dir)
  filename <- paste0(set_name,"_vsd.Rdata")
  if (file.exists(filename)) {
    cat(filename,"already exists\nloading prior...\n\n")
    lname <- load(filename)
  } else {
    cat(filename,"not found\nrunning vst now...\n\n")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- estimateSizeFactors(dds)
    vsd <- vst(dds, blind=TRUE)  ## *********************rate limiting
    save(vsd, file=filename)
  }
  ## define groups of interest
  single_factors <- c("category","disease_type","primary_site","tissue_type","stage","msi_status","vital_status","race","sex","capture_kit_name")
  intgroup <- list()
  for ( i in 1:length(single_factors) ) {
    intgroup[[i]] <- single_factors[i]
  }
  names(intgroup) <- c("aneuploidy quartile","cancer type","location","phenotype","stage","MSI","vital status","race","sex","hospital","capture kit")
  ## calculate pcs
  pcaData <- DESeq2::plotPCA(vsd, intgroup=intgroup[[1]], returnData=TRUE)  ## requires intgroup to set object
  percentVar <- round( 100 * attr(pcaData, "percentVar") )
  for (key in intgroup) {
    if (class(colData(vsd)[,key])=="character") {
      pcaData[,key] <- as.factor(colData(vsd)[,key])
    } else {
      pcaData[,key] <- colData(vsd)[,key]
    }
  }
  pcaData[,"null_lab"] <- ""
  ## make figs
  setwd(plot_dir)
  filename <- paste0(set_name,"_des_prlm_plt.pdf")
  pdf(file=filename)  ## set device
  for ( key in intgroup ) {
    print(pca_single(vsd, pcaData, percentVar, key, set_name, "null_lab"))
    print(pca_single(vsd, pcaData, percentVar, key, set_name, "name"))
  }
  dev.off()  ## close device to save plots
  cat("done exploring counts\n\n")
  end_time <- Sys.time()
  print(end_time - start_time)
}


deseq_naive <- function(dds,set_name,r_dir) {
  cat(paste("fitting naive model with DESeq2 for set",set_name),"\n\n")
  start_time <- Sys.time()
  setwd(r_dir)
  filename <- paste0(set_name,"_dds_naive.Rdata")
  if (file.exists(filename)) {
    cat(filename,"already exists\nloading prior...\n\n")
    lname <- load(filename)
  } else {
    cat(filename,"not found\nbuilding dds now...\n\n")
    colData(dds)[,"category"] <- as.factor(colData(dds)[,"category"])
    design(dds) <- ~ category
    count_lim <- 10
    sample_lim <- (1/5) * ncol(dds)
    keep <- rowSums( counts(dds) > count_lim ) >= sample_lim
    dds <- dds[keep,]
    dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(numWorkers))  ## estimate size factors and dispersions, fit negative binomial model, test
    save(dds, file=filename)
  }
  cat("done fitting naive model\n\n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(dds)
}


get_svs <- function(dds,set_name,choose_nsv=NULL) {
  cat(paste("estimating/getting svs for set",set_name),"\n\n")
  start_time <- Sys.time()
  ## check if surrogate variables already calculated; if so, then load, if not, then calculate
  if (refine<1) {
    setwd(r_dir)
    filename <- paste0(set_name,"_dds_svs_svA.Rdata")
  } else {
    setwd(r_dir)
    filename <- paste0(set_name,"_dds_svs_sv",choose_nsv,".Rdata")
  }
  if (file.exists(filename)) {
    cat(paste(filename,"already exists\nloading prior..."), "\n\n")
    lname <- load(file=filename)
    cat("done estimating/getting svs\n\n")
    end_time <- Sys.time()
    print(end_time - start_time)
    return(svs)
  } else {
    cat(paste(filename,"not found\nrunning svaseq() now..."), "\n\n")
    ## prepare data for sv calculation
    dat <- counts(dds, normalized=TRUE)
    mod <- model.matrix(~as.factor(category), data=colData(dds))   ## adjustment variables plus biological variables
    mod0 <- model.matrix(~1, data=colData(dds))      ## adjustment variables only
    if (refine<1) {
      ## estimate an appropriate number of surrogate variables to remove
      ### leek method
      n_sv <- num.sv(dat, mod, method=c("leek"))
      cat(paste("There appear to be",n_sv,"surrogate variables in the data with 'leek' method."),"\n\n")
      ### buja eyuboglu method
      n_sv <- num.sv(dat, mod, method=c("be"), B=20, seed=1)
      cat(paste("There appear to be",n_sv,"surrogate variables in the data with 'be' method."), "\n\n")
    } else {
      ## set preferred number of surrogate variables to remove
      n_sv <- choose_nsv
      cat(paste("Chosen number sv:",n_sv),"\n\n")
    }
    if (pn == 1) {
      cat("done estimating/getting svs\n\n")
      end_time <- Sys.time()
      print(end_time - start_time)
    }
    ## check pass number, only calc svs for 2nd pass
    if (pn == 2) {
      ### estimate surrogate variables
      if (refine<1) {
        cat(paste("Estimating surrogate variables with",n_sv,"surrogate variables per 'be' method."),"\n\n")
      } else {
        cat(paste("Estimating surrogate variables with",n_sv,"surrogate variables as requested."),"\n\n")
      }
      svs <- svaseq(dat, mod, mod0, n.sv=n_sv)
      cat("\n")
      ### save estimated svs
      save(svs, file=filename)
      cat("done estimating/getting svs\n\n")
      end_time <- Sys.time()
      print(end_time - start_time)
      return(svs)
    }
  }
}


deseq_sva <- function(dds,set_name,svs,r_dir) {
  ## fit informed model
  cat(paste("fitting model with svs for set",set_name), "\n\n")
  start_time <- Sys.time()
  ## check if sv-informed fit already calculated; if so, then load, if not, then calculate
  setwd(r_dir)
  if (refine<1) {
    filename <- paste0(set_name,"_dds_sva_svA.Rdata")
  } else {
    filename <- paste0(set_name,"_dds_sva_sv",svs$n.sv,".Rdata")
  }
  if (file.exists(filename)) {
    cat(paste(filename,"already exists\nloading prior..."), "\n\n")
    lname <- load(file=filename)
  } else {
    cat(paste(filename,"not found\nbuilding dds with svs now..."), "\n\n")
    ddssva <- dds  ## initiate data object
    ## add svs and set design matrix
    if (refine<1) {
      ## add num.sv svs to DESeq2 data object
      cat("adding automatic svs to model...\n\n")
      ddssva$sv1 <- svs$sv[,1]
      design(ddssva) <- ~ sv1 + category  ## simple
      if (set_name %in% c("simple")) {  ## simple needs sv2
        ddssva$sv2 <- svs$sv[,2]
        design(ddssva) <- ~ sv1 + sv2 + category
      }
      if (set_name %in% c("mss")) {  ## unknown whether any need sv3
        ddssva$sv3 <- svs$sv[,3]
        design(ddssva) <- ~ sv1 + sv2 + sv3 + category
      }
    } else {
      ## add chosen svs to DESeq2 data object
      cat("adding chosen svs to model...\n\n")
      ddssva$sv1 <- svs$sv[,1]
      ddssva$sv2 <- svs$sv[,2]
      ddssva$sv3 <- svs$sv[,3]
      design(ddssva) <- ~ sv1 + sv2 + sv3 + category
      if (svs$n.sv>3) {ddssva$sv4 <- svs$sv[,4]; design(ddssva) <- ~ sv1 + sv2 + sv3 + sv4 + category}
      if (svs$n.sv>4) {ddssva$sv5 <- svs$sv[,5]; design(ddssva) <- ~ sv1 + sv2 + sv3 + sv4 + sv5 + category}
    }
    ## fit model
    ddssva <- DESeq(ddssva, parallel=TRUE, BPPARAM=MulticoreParam(numWorkers))   ## estimate size factors and dispersions, fit negative binomial model, test
    save(ddssva, file=filename)
  }
  cat("done fitting model with svs\n\n")
  end_time <- Sys.time()
  print(end_time - start_time)
  return(ddssva)
}


## main script
dds <- get_dds(set_name,count_dir,r_dir)
deseq_explore(dds,set_name,r_dir,plot_dir)
dds <- deseq_naive(dds,set_name,r_dir)
if (pn==1) {
  get_svs(dds,set_name)  ## first pass returns no other data objects, only reports number significant svs
} else {  ## proceed only for 2nd pass
  if (refine<1) {svs <- get_svs(dds,set_name)} else {svs <- get_svs(dds,set_name,choose_nsv)}  ## get svs
  ddssva <- deseq_sva(dds,set_name,svs,r_dir)  ## fit informed model, DESeq2 + SVA
}
