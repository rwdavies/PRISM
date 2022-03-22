## load library
library("PRISM")



##
## 0 - setup
##

workdir <- "/data/smew1/rdavies/prism/"
dir.create(workdir)
setwd(workdir)
chrList <- c(1:22,"X") # list of chromosomes to analyze
nCores <- 23 # how many computer cores to use for multi-threaded options


## original set
fasta <- "human_g1k_v37.fasta.gz" # full path to fasta file
rmsk <- "repeatMasker.hg19.gz" # path to repeat masker file
simpleRepeats <- "simpleRepeat.hg19.gz" # path to simple repeats file
hotspotsFile <- "hotspots.RData" # original
hotspotsFile <- "AA_Maphotspots.RData" # one of the other ones
hotspotsFile <- "CEU_LDhotspots.RData"; chrList=1:22



## build 37 hapmap - download recombination rate files
system("wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz")

## build 37 - download fasta
system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz")

## need to manually download UCSC tracks
## to do, go to UCSC, to something like this
## https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=285525952_W7D6Vxjeh8vabp3ad5mj9rCHQGsd&clade=mammal&org=Human&db=hg19&hgta_group=allTables&hgta_track=hg19&hgta_table=simpleRepeat&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=primaryTable&hgta_outFileName=repeatMasker.hg19.gz
## and download for hg19, under group: All Tables, download simpleRepeat, and rmsk
stopifnot(file.exists(file.path(workdir, simpleRepeats)))
stopifnot(file.exists(file.path(workdir, rmsk)))

## setup fasta files. for each chromosome in chrList, make an RData object with the reference sequence for that chromosome
prepFasta(
    fasta = fasta,
    chrList = chrList,
    nCores = nCores
)


## setup repeat masker file
## for each chromosome in chrList, make an RData object with the repeat sequence for that chromosome
prepRepeat(
    rmsk = rmsk,
    simpleRepeats = simpleRepeats,
    chrList = chrList,
    nCores = nCores
)







## ## call hotspots - first, unpack some recombination rate map files
## system("tar -xzvf genetic_map_HapMapII_GRCh37.tar.gz")

## ## call hotspots - see ?callHotspots for full explanation and parameter options
## ## note - can use maxLength to filter out those that are too long
## ## or can use maxTrim to keep all hotspots but keep up to a width
## hot <- callHotspots(
##     prefix = "genetic_map_GRCh37_chr",
##     suffix = ".txt",
##     chrList = chrList,
##     nCores = nCores,
##     verbose = 2,
##     maxLength = 8000,
##     maxTrim = 2000
## )
## save(hot,file="hotspots.RData")



##
## african american hotspots
##
## system("wget http://www.well.ox.ac.uk/~anjali/AAmap/aamap.tar.gz")
## system("tar -xzvf aamap.tar.gz")

## ## fix 
## for(chr in chrList) {
##   data=read.table(paste("AAmap/AAmap.chr",chr,".txt",sep=""),header=TRUE)
##   data=cbind(data,1e6*c((data[-1,2]-data[-nrow(data),2])/(data[-1,1]-data[-nrow(data),1]),0))
##   colnames(data)[3]="rate"
##   write.table(data,file=paste("AAmap/AAmap.chr",chr,".modified.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
## }

## hot <- callHotspots(
##     prefix = "AAmap/AAmap.chr",
##     suffix = ".modified.txt",
##     chrList = chrList,
##     nCores = nCores,
##     verbose = 2,
##     maxLength = 3000,
##     positionColumn = "Physical_Position_Build36.hg18.",
##     rateColumn = "rate"
## )
## save(hot, file = "AAhotspots.RData")

##
## call broad hotspot set
##

system("wget http://www.well.ox.ac.uk/~anjali/AAmap/maps_b37.tar.gz")
system("tar -zxvf maps_b37.tar.gz")

## fix / update
for(chr in chrList) {
  data=read.table(paste("maps_b37/maps_chr.",chr,"",sep=""),header=TRUE)
  # for each column, make new rate column
  for(col in c("deCODE","COMBINED_LD","YRI_LD","CEU_LD","AA_Map","African_Enriched","Shared_Map")) {
    data=cbind(data,1e6*c((data[-1,col]-data[-nrow(data),col])/(data[-1,1]-data[-nrow(data),1]),0))
    colnames(data)[ncol(data)]=paste(col,"_rate",sep="")
  }
  write.table(data,file=paste("maps_b37/maps_chr.",chr,".modified.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}

## call hotspots in each
##for(col in c("deCODE","COMBINED_LD","YRI_LD","CEU_LD","AA_Map","African_Enriched","Shared_Map")) {
for(col in c("CEU_LD")) {    
    hot <- callHotspots(
        prefix = "maps_b37/maps_chr.",
        suffix = ".modified.txt",
        chrList = chrList,
        nCores = nCores,
        verbose = 2,
        maxLength = 3000,
        positionColumn = "Physical_Pos",
        rateColumn = paste(col,"rate",sep = "_")
    )
    save(hot,file=paste(col,"hotspots.RData",sep=""))
}

## check
for(col in c("deCODE","COMBINED_LD","YRI_LD","CEU_LD","AA_Map","African_Enriched","Shared_Map")) {
    file <- paste(col,"hotspots.RData",sep="")
    if (file.exists(file)) {
        load(file=file)
        print(c(col,nrow(hot)))
    }
}




  








##
## 1 - find some motifs
##

## find a highly enriched k-mer in the hotspots, with central enrichment, and use this to find one motif
setwd(workdir)
##hotspotsFile <- "AAhotspots.RData" # one of the other ones
hotspotsFile <- "CEU_LDhotspots.RData"
load(hotspotsFile)
use <- array(TRUE,nrow(hot))
dir <- "plots6/"
system(paste("mkdir ",dir,sep=""))
for(iMotif in 1:2) {
  hotExclude <- NULL
  plotPrefix <- paste(dir,"human.motif",iMotif,sep="")
  if(iMotif>1) hotExclude <- hot[use==FALSE,]
  hotX <- hot
  rm(hot)
  out <- getMotifs(
      hot = hotX[use,],
      hotExclude = hotExclude,
      eligible = NULL,
      plotting = "pdf",
      plotPrefix = plotPrefix,
      verbose = 2,
      nCores = nCores,
      filterSequences = TRUE,
      fasta = fasta,
      rmsk = rmsk,
      simpleRepeats = simpleRepeats,
      maxits = 10,
      initializationMethod = "simple",
      pwmValue = 0.9,
      hotspotCenterDist = 500
  )
  hot = hotX
  ## find second motif, after excluding first hotspots
  save(out, file=paste0("out.",iMotif))
  use[use==TRUE][out$regprob>0.5]=FALSE
}

quit()










##
## old scratch uncertain code
##


load(file=paste("out.",3,sep=""))
motifEnrichment=out$motifEnrichment[1:3]
lapply(motifEnrichment,head)

compX1=motifEnrichment$compX1
compX2=motifEnrichment$compX2
compX3=motifEnrichment$compX3

seqs=out$motifEnrichment[[4]]

#out=getMotifs(hot=hot,hotExclude=NULL,eligible=NULL,plotting="pdf",plotPrefix=plotPrefix,verbose=2,nCores=nCores,filterSequences=TRUE,rmsk=rmsk,simpleRepeats=simpleRepeats,maxits=10,initializationMethod="existing",pwmValue=0.9)
scorematset=out$scoremat
dimvec=out$scorematdim
pwm=out$pwm


which.max(out$regprobs>0.90)

1




hotExclude=NULL;eligible=NULL;nCores=nCores;K=8;filterSequences=TRUE;initializationMethod="simple";pwmValue=0.7;alpha=0.5;incprob=0.99999;maxits=3;plen=0.05;updatemot=1;updatealpha=1;ourprior=NULL;updateprior=1;bg=-1;plotting="pdf";plotPrefix="human.motif1";verbose=2;hotspotCenterDist=200

system("rsync -av /data/smew1/rdavies/pimp/ dense:/data/scratch/rwdavies/pimp/")


# hotspot 48, near 742 in. very GT rich
ref[5637331+ 0:9 + 742]
seqsL[[48]][0:9 + 1+742]
mask[-2:2+which.min(abs(5637331+742-mask[,"genoStart"])),]
ref[5638132:5638195] # this sequence is identified as (TG)n

rep[-2:2+which.min(abs(5637331+742-rep[,"chromStart"])),]
ref[5637331+ 0:50+742]
ref[5637331+ 0:50+742]
ref[5638133:5638161]







###
### 3 - run both motifs together 
###


### NOTE - need to write in version that either does, or does not use discovery





###
### 4 - run final iteration, no filtering for sequences
###
