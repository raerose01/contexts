#!/opt/R/R-latest/bin/Rscript --no-save

suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
source("/home/rosenthal/scripts/multi/findComp.R")

add_context=function(maf_loc) {
  #add_context reads in the maf and adds 4 contexts columns
  #read in the maf and add in the additional columns
  allmut = read.table(maf_loc, sep="\t", header=T, fill=T, quote="")
  allmut$context=getSeq(Hsapiens,paste("chr",allmut$Chromosome,sep=''),allmut$Start_position-3, allmut$End_position+3,as.character=T)
  allmut$mutcat = noquote(paste(allmut$Reference_Allele,">",allmut$Tumor_Seq_Allele1,sep=""))
  allmut$onecontext = noquote(paste(substr(allmut$context,3,3),"x",substr(allmut$context,nchar(allmut$context)-2,nchar(allmut$context)-2),sep=""))

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

  return(onlySNPs)

}

exonicSNPify=function(exn) {
  allSNPs = add_context(exn)
  #only take exonic regions
  exonicSNPs = subset(allSNPs,allSNPs$Variant_Classification == "Silent" | 
  allSNPs$Variant_Classification == "Nonsense_Mutation" | 
  allSNPs$Variant_Classification == "Translation_Start_Site" | 
  allSNPs$Variant_Classification == "Missense_Mutation" | 
  allSNPs$Variant_Classification == "Splice_Site" | 
  allSNPs$Variant_Classification == "Nonstop_Mutation")

  return(exonicSNPs)

}

tricontext = function(fin) {

  #tricontext normalizes by the threemer matrix counts

  exonicSNPs = exonicSNPify(fin)

  #read in the threemers table for future normalization
  threemers = read.table("/home/rosenthal/values/exon_triCount.txt", sep="\t", header = T, fill=T, quote="", row.names=1)

  #table on the mutcat and context on either side, then make a matrix of it
  x = table(exonicSNPs$revmutcat,exonicSNPs$revcontext)
  y = matrix(x,ncol(x),nrow(x),dimnames=list(colnames(x),rownames(x)),byrow=T)

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
  #test = matrix("NA", nrow=length(r), ncol=length(c), dimnames=list(r,c))
  for(i in 1:length(r)){
  for(j in 1:length(c)){
  a = threemers[grep(m[i,j],rownames(threemers)),]
  b = threemers[grep(mr[i,j],rownames(threemers)),]
  #test[i,j] = a+b
  final[i,j]=10^6*(y[i,j]/(a+b))
  }}

  return(final)

}






