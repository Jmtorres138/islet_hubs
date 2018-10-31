
"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")
#serv.dir <- "/Users/jtorres/FUSE/"
serv.dir <- "/well/got2d/jason/"

work.dir <- serv.dir %&% "projects/islet_hubs/fgwas/oxford_and_imperial/"
fgwas.output.dir <- work.dir %&% "fgwas_output/"
cond.dir <- work.dir %&% "conditional/"
input.dir <- cond.dir %&% "fgwas_input_files/"

annot.file <- work.dir %&% "fgwas_input/" %&% "anno_input.bed"
annot.df <- fread(annot.file)

best.df <- fread(fgwas.output.dir %&% "best-joint-model.params")
annots <- unique(best.df$parameter)
annots <- gsub("_ln","",annots)

filt.annot.df <- filter(annot.df,V4 %in% annots)

write.table(x=filt.annot.df,file=input.dir%&%"best_anno_input.bed",sep="\t",quote=FALSE,
            row.names=FALSE,col.names=FALSE)
