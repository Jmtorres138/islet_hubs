
#Constructing functional (fgwas) credible sets from the fgwas model with the best cross validated likelihood

# Setup
args = commandArgs(trailingOnly=TRUE)
loc.id <- args[1]

"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")
library("ggplot2")
library(GenomicRanges)

serv.dir <- "/well/mccarthy/users/jason/"

work.dir <- serv.dir %&% "projects/islet_hubs/fgwas/oxford_and_imperial/"
fgwas.output.dir <- work.dir %&% "conditional/fgwas_output_files/" %&% loc.id %&% "/"
#pre <- fgwas.output.dir %&% "fgwas_run_loci-partition"


cred.set.dir <- serv.dir %&% "projects/islet_hubs/fgwas/imperial_annotations/credible_sets/"
#output_file = work_dir+"conditional/fgwas_input_files/" + loc_id + "/" +"loci_block_snps.bfs.txt.gz"

print("Making credible set file for locus: " %&% loc.id)


get_cred <- function(dframe,cname,prob=0.99){ # 99% functional credible sets
  index <- match(cname,names(dframe))
  vec <- sort(dframe[,index],decreasing=TRUE)
  count=0
  sum=0
  for (v in vec){
    count <- count + 1
    sum <- sum + v
    if (sum >= prob){
      break
    }
  }
  return(count)
}

# Determine functional credible sets (95% and 99% Credible Sets)

get_credsets <- function(probthresh,segsnps.df){
  seg.vec <- sort(unique(segsnps.df$SEGNUMBER))
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(seg.vec),style = 3)
  for (i  in 1:length(seg.vec)){
    setTxtProgressBar(pb,i)
    seg <- seg.vec[i]
    temp.df <- filter(segsnps.df,SEGNUMBER==seg) %>% arrange(desc(PPA))
    cumppa <- sum(temp.df$PPA)
    temp.df <- temp.df[1:get_cred(temp.df,"PPA",prob=probthresh*cumppa),]
    temp.df <- dplyr::select(temp.df,-pi,-chunk,-pseudologPO,-pseudoPPA,-V)
    names(temp.df)[c(2,3,4)] <- c("SNPID","CHR","POS")
    temp.df$PPA <- temp.df$PPA/cumppa # rescale ppa to reflect proportion of cummulative sum in block
    out.df <- rbind(out.df,temp.df)
  }
  return(out.df)
}

# Get nearest genes

library(GenomicFeatures)
library(org.Hs.eg.db)
library(annotate)

annot_refGene <- function(cred.df,segsnps.df){
  seg.vec <- sort(unique(segsnps.df$SEGNUMBER))
  snp.gr <- GRanges(cred.df$CHR,IRanges(cred.df$POS, cred.df$POS))

  #hg19.refseq.db <- readRDS("hg19.refseq.db.RDS")#makeTxDbFromUCSC(genome="hg19", table="refGene")
  #refseq.genes<- genes(hg19.refseq.db)
  #all.geneids <- elementMetadata(refseq.genes)$gene_id
  #all.genesymbols <- getSYMBOL(all.geneids, data='org.Hs.eg')
  #ref.df <- data.frame(gene.id=all.geneids,symbol=all.genesymbols,
  #                     stringsAsFactors = FALSE)
  #elementMetadata(refseq.genes)$gene_id <- all.genesymbols
  df <- readRDS("reference-genes.df.RDS")#as.data.frame(refseq.genes)
  sub.gr <- GRanges(seqnames=df$seqnames,
                    IRanges(start=df$start,end=df$end),
                    strand=rep("*",dim(df)[1]),
                    name=df$gene_id)
  nearestGenes <- nearest(snp.gr,sub.gr)
  res <- df$gene_id[nearestGenes]
  refseq.genes <- readRDS("refseq.genes.RDS")
  dist <- distance(snp.gr, refseq.genes[nearestGenes])
  symbol <- res
  cred.df <- cbind(symbol,cred.df)


  # Sync symbol names for each SEGNUMBER by MOST COMMON gene
  sync.df <- c()
  pb <- txtProgressBar(min=0,max=length(seg.vec),style=3)
  for (i in 1:length(seg.vec)){
    setTxtProgressBar(pb,i)
    seg <- seg.vec[i]
    temp.df <- filter(cred.df,SEGNUMBER==seg)
    top <- names(sort(table(temp.df$symbol),decreasing=TRUE))[1]
    temp.df$symbol <- rep(top,dim(temp.df)[1])
    sync.df <- rbind(sync.df,temp.df)
  }
  cred.df <- sync.df
  return(cred.df)
}

# Run


segsnps.df <- fread("cat " %&% fgwas.output.dir  %&% "loci_block_snps.bfs.txt.gz" %&% " | zmore")
cred95.df <- get_credsets(0.95,segsnps.df)
cred99.df <- get_credsets(0.99,segsnps.df)
cred95.df <- annot_refGene(cred95.df,segsnps.df)
cred99.df <- annot_refGene(cred99.df,segsnps.df)
write.table(x=cred95.df,file=cred.set.dir%&%"fgwas_credsets_95-" %&% loc.id %&% ".txt",sep="\t",
            quote=FALSE,row.names=F)
write.table(x=cred99.df,file=cred.set.dir%&%"fgwas_credsets_99-" %&% loc.id %&% ".txt",sep="\t",
            quote=FALSE,row.names=F)
