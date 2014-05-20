#!/opt/R/R-latest/bin/Rscript --no-save

suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
source("/home/rosenthal/scripts/multi/findComp.R")

#read in the maf and add in the additional columns
allmut = read.table("/taylorlab/data/oracle/tumor_oracle_Dec2013.maf", sep="\t", header=T, fill=T, quote="")
allmut$context=as.character(allmut$context)
allmut$Trns_Ctxt_Mod = as.character(allmut$Trns_Ctxt_Mod)
allmut$mutcat = noquote(paste(allmut$Reference_Allele,">",allmut$SNP_Mut_Allele,sep=""))
allmut$onecontext = noquote(paste(substr(allmut$context,3,3),"x",substr(allmut$context,nchar(allmut$context)-2,nchar(allmut$context)-2),sep=""))

#read in the threemers table for future normalization
#threemers = read.table("/home/rosenthal/scripts/multi/hg19_nuc_count_3.txt", sep="\t", header = F, fill=T, quote="", row.names=1)
threemers = read.table("/home/rosenthal/values/exon_triCount.txt", sep="\t", header = T, fill=T, quote="", row.names=1)

#only look at the SNPs from the maf
onlySNPs = subset(allmut,allmut$Variant_Type=="SNP" | allmut$Variant_Type=="SNV")

#reverse complement G and T mutations (and context) to C and A
allmut$revmutcat = noquote(substr(onlySNPs$Trns_Ctxt_Mod,3,5))
allmut$revcontext = noquote(paste(substr(allmut$Trns_Ctxt_Mod,1,1),"x",substr(allmut$Trns_Ctxt_Mod,nchar(allmut$Trns_Ctxt_Mod),nchar(allmut$Trns_Ctxt_Mod)),sep=""))
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

write.table(onlySNPs,file="/home/rosenthal/oracle_onlySNPs",sep="\t")

#only take exonic regions
#exonicSNPs = subset(onlySNPs,onlySNPs$Variant_Classification == "Silent" | onlySNPs$Variant_Classification == "Nonsense_Mutation" | onlySNPs$Variant_Classification == "Translation_Start_Site" | onlySNPs$Variant_Classification == "Missense_Mutation" | onlySNPs$Variant_Classification == "Splice_Site" | onlySNPs$Variant_Classification == "Nonstop_Mutation")

#table on the mutcat and context on either side, then make an array of it -- this time divided by sample
#samples=unique(allmut$TUMORTYPE_REV)
#lm=sort(unique(exonicSNPs$revmutcat))
#lc=sort(unique(exonicSNPs$revcontext))
#t = array(dim=c(6,16,length(samples)),dimnames=list(lm,lc,samples))

#for(i in samples){
#tm=exonicSNPs[ exonicSNPs$Sample == i, ]
#x = table(factor(tm$revmutcat,levels=lm),factor(tm$revcontext,levels=lc))
#y = matrix(x,length(rownames(x)),length(colnames(x)),dimnames=list(rownames(x),colnames(x)),byrow=F)
#t[,,i] = y
#}

#r = rownames(t)
#c = colnames(t)


#make matrix of each of the trinucleotide contexts
#tm=list()
#for(i in 1:length(r)){
#for(j in 1:length(c)){
#tm = c(tm,paste(substr(r[i],1,1),substr(c[j],1,1),substr(r[i],3,3),sep=""))
#}}

#m = matrix(tm,nrow=length(r),ncol=length(c),byrow=F)


# normalizing the raw count of mutation by the appropriate (sum for the complements) trinucleotide context
#final = t

#for(z in 1:length(samples)){
#for(i in 1:length(r)){
#for(j in 1:length(c)){
#final[i,j,z]=10^6*(t[i,j,z]/(threemers[grep(m[i,j],rownames(threemers)),]))
#}}}
