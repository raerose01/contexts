#!/opt/R/R-latest/bin/Rscript --no-save

suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
source("/home/rosenthal/scripts/multi/findComp.R")

#read in the maf and add in the additional columns
allmut = read.table("/home/sasthana/data/liver_patient/exome_06_2013/output/LP_tumor/LP_tumor.all_muts.maf", sep="\t", header=T, fill=T, quote="")
allmut$context=getSeq(Hsapiens,paste("chr",allmut$Chromosome,sep=''),allmut$Start_position-3, allmut$End_position+3,as.character=T)
allmut$mutcat = noquote(paste(allmut$Reference_Allele,">",allmut$Tumor_Seq_Allele1,sep=""))
allmut$onecontext = noquote(paste(substr(allmut$context,3,3),"x",substr(allmut$context,nchar(allmut$context)-2,nchar(allmut$context)-2),sep=""))

#read in the threemers table for future normalization
#threemers = read.table("/home/rosenthal/scripts/multi/hg19_nuc_count_3.txt", sep="\t", header = F, fill=T, quote="", row.names=1)
threemers = read.table("/home/rosenthal/values/exon_triCount.txt", sep="\t", header = T, fill=T, quote="", row.names=1)

#only look at the SNPs from the maf
onlySNPs = subset(allmut,allmut$Variant_Type=="SNP")

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
write.table(onlySNPs,"/home/rosenthal/projects/liver_patient/allmut_jcontext.maf", sep = "\t", quote=FALSE, row.names=F)

#only take exonic regions
exonicSNPs = subset(onlySNPs,onlySNPs$Variant_Classification == "Silent" | 
onlySNPs$Variant_Classification == "Nonsense_Mutation" | 
onlySNPs$Variant_Classification == "Translation_Start_Site" | 
onlySNPs$Variant_Classification == "Missense_Mutation" | 
onlySNPs$Variant_Classification == "Splice_Site" | 
onlySNPs$Variant_Classification == "Nonstop_Mutation")


#table on the mutcat and context on either side, then make a matrix of it
x = table(exonicSNPs$revmutcat,exonicSNPs$revcontext)
y = matrix(x,16,6,dimnames=list(colnames(x),rownames(x)),byrow=T)

r = rownames(y)
c = colnames(y)


#make matrix of each of the trinucleotide contexts
tm=list()
for(i in 1:length(r)){
for(j in 1:length(c)){
tm = c(tm,paste(substr(r[i],1,1),substr(c[j],1,1),substr(r[i],3,3),sep=""))
}}

#reverse complement matrix of each of the trinucleotide contexts
mt=list()
for(i in 1:length(tm)){
mt=c(mt,findComp_nr(as.character(tm[i])))
}

m = matrix(tm,nrow=length(r),ncol=length(c),byrow=T)
mr = matrix(mt,nrow=length(r),ncol=length(c),byrow=T)

# normalizing the raw count of mutation by the appropriate (sum for the complements) trinucleotide context
final = matrix("NA", nrow=length(r), ncol=length(c), dimnames=list(r,c))
test = matrix("NA", nrow=length(r), ncol=length(c), dimnames=list(r,c))
for(i in 1:length(r)){
for(j in 1:length(c)){
a = threemers[grep(m[i,j],rownames(threemers)),]
b = threemers[grep(mr[i,j],rownames(threemers)),]
test[i,j] = a+b
final[i,j]=10^6*(y[i,j]/(a+b))
}}

final
write.table(final, file="/home/rosenthal/projects/liver_patient/heatmap/trinucleotide_normalized_mutations2.txt", sep="\t")







