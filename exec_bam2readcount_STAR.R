rm( list=ls(all=TRUE) ) # clean up R workspace

library(Rsubread)
library(edgeR)
args <- commandArgs(trailingOnly=T)

organism_annotation <- args[1] # mm or hg etc
samplename <- args[2] # sample ID
thread_num <- args[3] # core number
sample_dir <- args[4] # input directory

inputfile <- paste(sample_dir,"/",samplename,".Aligned.out.bam",sep="")
fc <- featureCounts(files=inputfile,annot.ext=organism_annotation,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",nthreads=thread_num)
## read count ##
setwd(sample_dir)
sfn <- paste(samplename,"_rc.txt",sep="")
write.table(fc$counts,file=sfn,quote=F,col.names=F,sep="\t")
sfn <- paste(samplename,"_rclog.txt",sep="")
write.table(fc$stat,file=sfn,quote=F,col.names=F,sep="\t")
## rpkm ##
x <- DGEList(counts=fc$counts,genes=fc$annotation[,c("GeneID","Length")])
x_rpkm <- rpkm(x,fc$genes$Length)
sfn <- paste(samplename,"_rpkm.txt",sep="")
write.table(x_rpkm,file=sfn,quote=F,col.names=F,sep="\t")

