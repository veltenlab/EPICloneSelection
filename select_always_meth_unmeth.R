##################################### select_always_meth_unmeth.R #####################################
#' This file computes sites that are always unmethylated/methylated in all the comparisons that
#' were conducted with RnBeads.
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

report <- '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20211020_differential/'
output <- '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/meth_unmeth/'
config_file <- '/users/lvelten/project/Methylome/src/selection_pipeline/config.yaml'
config <- yaml.load_file(config_file)

#' checkForCutSite
#'
#' This function checks for the presence of a cutsite in an amplicon using the configuration provided in the config file
#' @param input The input data frame
#' @param config The parsed configuration
#' @return A modified version of \code{input} with only those lines containing a cutsite for the specified enzyme
checkForCutSite <- function(input, config, number){
    cut.seq <- DNAString(config[['general']][['cut_seq']])
    gen.version <- config[['general']][['genome']]
    if(gen.version=='mm10'){
        suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
        genome <- BSgenome.Mmusculus.UCSC.mm10
    }else if(gen.version=='mm10'){
        suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
        genome <- BSgenome.Hsapiens.UCSC.hg19
    }else{
        stop("Genome version '", gen.version, "' not supported")
    }
    dat <- read.csv(input)
    dat <- dat[sample(row.names(dat),nrow(dat)),]
    i <- 1
    num <- 0
    vec.all <- c()
    dat$regionStart <- NA
    dat$regionEnd <- NA
    dat$cutsiteInRegion <- NA
    dat$gene <- NA
    dat$promoter <- NA
    dat$enhancer_annotation <- NA
    dat$enhancer_annotation_gene_name <- NA
    dat$GCContent <- NA
    tfs.files.all <- list.files('/users/lvelten/project/Methylome/infos/BCells/TFBS_Wilson/bed/mm10/',full.names=TRUE)
    tfs.all <- gsub('.bed', '', list.files('/users/lvelten/project/Methylome/infos/BCells/TFBS_Wilson/bed/mm10/'))
    for(tf in tfs.all){
        dat[,tf] <- NA
    }
    genes <- unlist(rnb.get.annotation('genes',assembly = 'mm10'))
    promoters <- unlist(rnb.get.annotation('promoters',assembly = 'mm10'))
    reg.elements <- read.csv('/users/lvelten/project/Methylome/infos/BCells/enhancer_catalog_lawrence/H3K27Ac.csv')
    row.names(reg.elements) <- paste0(reg.elements$Chr,'_',reg.elements$Start)
    reg.elements.tab <- read.table('/users/lvelten/project/Methylome/infos/BCells/enhancer_catalog_lawrence/H3K27Ac_mm10.bed',sep='\t')
    reg.elements.gr <- makeGRangesFromDataFrame(reg.elements.tab,seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V4')                    
    values(reg.elements.gr) <- reg.elements.tab$V5
    while(num<number && i<nrow(dat)){
        chr <- dat[i,"Chromosome"]
        start <- as.numeric(as.character(dat[i,"Start"]))
        seq <- genome[[chr]]
        region.start <- start-2
        region.end <- start+2
        sel.seq <- seq[region.start:region.end]
        res <- matchPattern(cut.seq,sel.seq)
        extended.region.start <- region.start-200
        extended.region.end <- region.end+200
        extended.region.seq <- seq[extended.region.start:extended.region.end]
        res.extended <- matchPattern(cut.seq,extended.region.seq)
        if(length(res)==1 & length(res.extended)==1){
            dat$regionStart[i] <- region.start
            dat$regionEnd[i] <- region.end
            dat$cutsiteInRegion[i] <- paste0(start(res), '-', end(res))
            freq <- alphabetFrequency(extended.region.seq)
            dat$GCContent[i] <- (freq['C']+freq['G'])/length(extended.region.seq)
            region.gr <- GRanges(paste0(chr,':',region.start,'-',region.end))
            op <- findOverlaps(region.gr,genes)
            if(length(op)>0){
                dat$gene[i] <- paste(values(genes)$symbol[subjectHits(op)],collapse=',')
            }else{
                dat$gene[i] <- NA
            }
            op <- findOverlaps(region.gr,promoters)
            if(length(op)>0){
                dat$promoter[i] <- paste(values(promoters)$symbol[subjectHits(op)],collapse=',')
            }else{
                dat$promoter[i] <- NA
            }
            op <- findOverlaps(region.gr,reg.elements.gr)
            if(length(op)>0){
                reg.element <- unname(unlist(values(reg.elements.gr))[subjectHits(op)])
                dat$enhancer_annotation[i] <- paste(reg.elements[reg.element,'Annotation'],collapse=',')
                dat$enhancer_annotation_gene_name[i] <- paste(reg.elements[reg.element,'Gene.Name'],collapse=',')
            }else{
                dat$enhancer_annotation[i] <- NA
                dat$enhancer_annotation_gene_name[i] <- NA
            }
            for(j in 1:length(tfs.all)){
                tf.gr <- makeGRangesFromDataFrame(read.table(tfs.files.all[[j]],sep = '\t'),seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
                op <- findOverlaps(region.gr,tf.gr)
                if(length(op)>0){
                    dat[i,tfs.all[j]] <- tfs.all[j]
                }else{
                    dat[i,tfs.all[j]] <- NA
                }
            }
            vec.all <- c(vec.all,i)
            num <- num+1
        }
        i <- i+1
    }
    dat <- dat[vec.all,]
    return(dat)
}

all.comparisons <- list.files(file.path(report,'differential_methylation_data'), full.names=TRUE, pattern = 'diffMethTable_site')
system(paste0('rm -rf ', output,'/always_meth.csv'))
system(paste0('rm -rf ', output,'/always_unmeth.csv'))
dmrs <- lapply(all.comparisons,function(comp){
    print(comp)
    diff.table <- as.data.frame(fread(comp))
    group.names <- colnames(diff.table)[grepl('sd.',colnames(diff.table))]
    group.names <- gsub('sd.','', group.names)
    print(paste('Processing diffMethTable for ', group.names[1], 'and', group.names[2]))
    sds.first <- diff.table[[paste0('sd.',group.names[1])]]
    sds.second <- diff.table[[paste0('sd.',group.names[2])]]
    se.sd.first <- sd(sds.first)#/sqrt(length(sds.first))
    se.sd.second <- sd(sds.second)#/sqrt(length(sds.second))
    low.sds <- (sds.first<(mean(sds.first)+2*se.sd.first)) & (sds.second<(mean(sds.second)+2*se.sd.second))
    diff.table <- diff.table[low.sds,]
    row.names(diff.table) <- paste0(diff.table$Chromosome, '_', diff.table$Start)
    diff.positive <- diff.table[diff.table[[paste0('mean.',group.names[1])]]>config[['dmcs']][['methylated']]&
                                diff.table[[paste0('mean.',group.names[2])]]>config[['dmcs']][['methylated']],]
    diff.positive$type <- 'always_meth'
    row.names(diff.positive) <- paste0(diff.positive$Chromosome, '_', diff.positive$Start)
    file.positive <- file.path(output,'always_meth.csv')
    if(file.exists(file.positive)){
        all.positive <- read.csv(file.positive,row.names = 1)
        if(nrow(all.positive)<4){
            stop('Please delete the always_meth file in the output directory')
        }
        both <- intersect(row.names(diff.positive),row.names(all.positive))
        all.positive <- data.frame(all.positive[both,],diff.positive[both,])
        write.csv(all.positive,file.positive)
    }else{
        write.csv(diff.positive,file.positive)
    }
    diff.negative <- diff.table[diff.table[[paste0('mean.',group.names[1])]]<config[['dmcs']][['unmethylated']]&
                                    diff.table[[paste0('mean.',group.names[2])]]<config[['dmcs']][['unmethylated']],]
    diff.negative$type <- 'always_unmeth'
    row.names(diff.negative) <- paste0(diff.negative$Chromosome, '_', diff.negative$Start)
    file.negative <- file.path(output,'always_unmeth.csv')
    if(file.exists(file.negative)){
        all.negative <- read.csv(file.negative, row.names=1)
        if(nrow(all.negative)<4){
            stop('Please delete the always_meth file in the output directory')
        }
        both <- intersect(row.names(diff.negative),row.names(all.negative))
        all.negative <- data.frame(all.negative[both,],diff.negative[both,])
        write.csv(all.negative,file.negative)
    }else{
        write.csv(diff.negative,file.negative)
    }
    return(list(diff.positive,diff.negative))
})

always_meth <- list.files(output,pattern='always_meth.csv', full.names=TRUE)
always_meth <- checkForCutSite(always_meth, config, 50)
write.csv(always_meth, file.path(output, 'always_meth_filtered.csv'))
always_unmeth <- list.files(output,pattern='always_unmeth.csv', full.names=TRUE)
always_unmeth <- checkForCutSite(always_unmeth, config, 50)
write.csv(always_unmeth, file.path(output, 'always_unmeth_filtered.csv'))

