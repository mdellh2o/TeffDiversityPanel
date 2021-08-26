options(stringsAsFactors = F)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/output"
setwd(wd)
  
#set thresholds
e.thr<-1e-20
perc.id.thr<-50
perc.sub.thr<-0

#set naming convention
newnames<-c("Et_id", "Et_protein_length", "protein_length", "id", "query_start",
            "query_end", "subject_start", "query seq", "subject seq" , "perc_identity", "evalue", "bit_score", "alignment_length", "identical", "mismatches", "gap_opens", 
            "perc_query_coverage_per_subject")

#get genes in blocks
blks<-read.table("../input/teff.genes.in.blocks.FDR005.txt", header=T,sep="\t")
blks<-blks[order(blks$attributes),]
#watch out as the file contains duplicates!
blks<-unique(blks)
#clean the attributes
attr<-blks$attributes
attr<-sub("\\;.*$", "", attr)
attr<-sub("^ID\\=", "", attr)
blks$attributes<-attr
head(blks)

blks$startstop<-paste0("lcl|",blks$seqid,":", blks$start, "-", blks$end)
dim(blks)

#get At blasts resutls
tab<-read.table("../input/blast.teff.arabidopsis.20210729.txt", header=F, comment.char="#", sep="\t")
dim(tab)
tab<-unique(tab)
dim(tab)

namesAt<-paste0(newnames, "_At")
colnames(tab)<- namesAt

tab<-data.frame(tab[,c(1:7, 10:17)])
head(tab)

#do the first merge
blstblk<-merge(blks, tab, by.x="startstop", by.y="Et_id_At", all=T)
head(blstblk)

#get Zm blasts resutls
tabZm<-read.table("../input/blast.teff.maize.20210729.txt", header=F, comment.char="#", sep="\t")
dim(tabZm)

tabZm<-unique(tabZm)
dim(tabZm)

namesZm<-paste0(newnames, "_Zm")
colnames(tabZm)<- namesZm

tabZm<-data.frame(tabZm[,c(1:7, 10:17)])
head(tabZm)

#do the second merge
blstblk2<-merge(blstblk, tabZm, by.x="startstop", by.y="Et_id_Zm", all.x=T)
head(blstblk2)

#subset the hits by chosen thresholds
e.thr<-1e-20
id.thr<-50
sub.thr<-10

#put a flag for At and for Zm
blstblk2$flagAt<-0
hitsAt<-which(blstblk2$evalue_At<e.thr & blstblk2$perc_identity_At > id.thr & blstblk2$perc_query_coverage_per_subject_At > sub.thr)
length(hitsAt)
blstblk2$flagAt[hitsAt]<-1

blstblk2$flagZm<-0
hitsZm<-which(blstblk2$evalue_Zm<e.thr & blstblk2$perc_identity_Zm > id.thr & blstblk2$perc_query_coverage_per_subject_Zm > sub.thr)
length(hitsZm)
blstblk2$flagZm[hitsZm]<-1

#now extract only significant blocks
sigblk<-read.delim("../input/sig.blocks")
head(sigblk)
dim(blstblk2)

out<-merge(blstblk2, sigblk, by.x="index", by.y="Blocks_ID", all.x=T)
head(out)
length(which(is.na(out$num_QTN)))

#get output
write.table(out, file="ALL.candidate.homology.At.Zm.txt", sep="\t", quote=F, row.names=F)

#pre-filter output
#tokeep<-unique(c(which(out$flagAt == 1), which(out$flagZm == 1)))
#outsub<-out[tokeep, ]
outsub<-out[which(out$num_QTN >= 1), ]
dim(outsub)

#get flitered output
write.table(outsub, file="SUB.candidate.homology.At.Zm.txt", sep="\t", quote=F, row.names=F)


#do some parsing
#identify prioirity blocks
blklen<-sort(table(outsub$index))
blklen
