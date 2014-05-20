suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
source("/home/rosenthal/scripts/multi/findComp.R")

#read in the maf and add in the additional columns
allmut = read.table("/home/rosenthal/projects/liver_patient/HoangAllData.txt", sep="\t", header=T, fill=T, quote="")
allmut$mutcat = as.character(allmut$SBS.on.Nontranscribed.Strand)
allmut$context = allmut$SBS.Sequence.Context.on.Nontranscribed.Strand..flanking.....10.bases.
allmut$onecontext = noquote(paste(substr(allmut$context,10,10),"x",substr(allmut$context,12,12),sep=""))

#read in the threemers table for future normalization
threemers = read.table("/home/rosenthal/values/exon_triCount.txt", sep="\t", header = T, fill=T, quote="", row.names=1)

#only look at the SNPs from the maf
onlySNPs = subset(allmut,allmut$DNA.Mutation.Type=="SBS")

#reverse complement G and T mutations (and context) to C and A
onlySNPs$revmutcat = onlySNPs$mutcat
onlySNPs$revcontext = onlySNPs$onecontext
gind = grep("G",substr(onlySNPs$mutcat,1,1))
tind = grep("T",substr(onlySNPs$mutcat,1,1))

for(i in 1:length(gind)){
  onlySNPs$revmutcat[gind[i]]=findComp_nr(onlySNPs$mutcat[gind[i]])
  onlySNPs$revcontext[gind[i]]=findComp(onlySNPs$onecontext[gind[i]])
}

for(i in 1:length(tind)){
  onlySNPs$revmutcat[tind[i]]=findComp_nr(onlySNPs$mutcat[tind[i]]) 
  onlySNPs$revcontext[tind[i]]=findComp(onlySNPs$onecontext[tind[i]])
}

exonicSNPs=onlySNPs

#table on the mutcat and context on either side, then make an array of it -- this time divided by sample
samples=unique(allmut$Sample)
lm=sort(unique(exonicSNPs$revmutcat))
lc=sort(unique(exonicSNPs$revcontext))
t = array(dim=c(6,16,length(samples)),dimnames=list(lm,lc,samples))

for(i in samples){
tm=exonicSNPs[ exonicSNPs$Sample == i, ]
x = table(factor(tm$revmutcat,levels=lm),factor(tm$revcontext,levels=lc))
y = matrix(x,length(rownames(x)),length(colnames(x)),dimnames=list(rownames(x),colnames(x)),byrow=F)
t[,,i] = y
}

r = rownames(t)
c = colnames(t)


#make matrix of each of the trinucleotide contexts
tm=list()
for(i in 1:length(r)){
for(j in 1:length(c)){
tm = c(tm,paste(substr(r[i],1,1),substr(c[j],1,1),substr(r[i],3,3),sep=""))
}}

m = matrix(tm,nrow=length(r),ncol=length(c),byrow=F)


# normalizing the raw count of mutation by the appropriate (sum for the complements) trinucleotide context
final = t

for(z in 1:length(samples)){
for(i in 1:length(r)){
for(j in 1:length(c)){
final[i,j,z]=10^6*(t[i,j,z]/(threemers[grep(m[i,j],rownames(threemers)),]))
}}}

allcol=c("yellow","springgreen2","red","darkmagenta","orange","skyblue1")
pdf("trinucleotide_contexts_perpatient.pdf")
for(p in 1:length(samples)){
test = final[,,p]
barplot(test,beside=T,cex.names=0.77,legend=rownames(test),names.arg=colnames(test),main=samples[p],ylab="Mutation Rate (per million tricontext)",col=allcol,las=2)
}
dev.off()

write.table(final, file="/home/rosenthal/projects/liver_patient/heatmap/Hoang_compare/Hoang_trinucleotide_normalized_mutations.txt", sep="\t")







