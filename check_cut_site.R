suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Biostrings))
pars <- ArgumentParser()
pars$add_argument("-i","--input",type="character",
    help="The input file")
pars$add_argument("-c","--config",type="character",
    help="Path to the configuration file")
pars$add_argument("-n","--number",type="integer",
    help="Number of amplicons to be generated")
args <- pars$parse_args()
config <- read_yaml(args$config)
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
dat <- read.csv(args$input)
i <- 1
num <- 0
vec.all <- c()
dat$regionStart <- NA
dat$regionEnd <- NA
dat$cutsiteInRegion <- NA
while(num<=args$number){
    chr <- dat[i,"Chromosome"]
    start <- as.numeric(as.character(dat[i,"Start"]))
    seq <- genome[[chr]]
    region.start <- start-100
    region.end <- start+100
    sel.seq <- seq[region.start:region.end]
    res <- matchPattern(cut.seq,sel.seq)
    if(length(res)==1){
        dat$regionStart[i] <- region.start
        dat$regionEnd[i] <- region.end
        dat$cutsiteInRegion[i] <- paste0(start(res),'-',end(res))
        vec.all <- c(vec.all,i)
        num <- num+1
    }
    i <- i+1
}
dat <- dat[vec.all,]
write.csv(dat,args$input)
