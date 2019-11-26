library(stringr)
library(dplyr)
library(pbapply)
library(parallel)
library(DECIPHER)

metax2tax <- function(metax=NULL,out=NULL){
  tab=read.table(metax,sep="\t",colClasses = "character")
  colnames(tab)=c("ID","tax")
  nc=parallel::detectCores()
  print("Converting taxon formats to tabular format, can take a while...")
  lst=pbapply::pblapply(1:length(tab$tax),function(x){
    str=stringr::str_split(tab$tax[x],"__subtri__|__super__|__sub__|;|__",simplify = T)
    return(data.frame(ID=tab$ID[x],
                      kingdom=ifelse(length(which(str=="k"))>0,str[which(str=="k")+1],NA),
                      phylum=ifelse(length(which(str=="p"))>0,str[which(str=="p")+1],NA),
                      class=ifelse(length(which(str=="c"))>0,str[which(str=="c")+1],NA),
                      order=ifelse(length(which(str=="o"))>0,str[which(str=="o")+1],NA),
                      family=ifelse(length(which(str=="f"))>0,str[which(str=="f")+1],NA),
                      genus=ifelse(length(which(str=="g"))>0,str[which(str=="g")+1],NA),
                      species=ifelse(length(which(str=="s"))>0,str[which(str=="s")+1],NA),
                      stringsAsFactors = F))
  },cl = nc-2)
  print("Concatenating dataframe.")
  lst=do.call(rbind,lst)
  write.table(lst,out,sep="\t",quote = F,row.names = F)
}

tax2metax <- function(tax=NULL,out=NULL){
  tab=read.table(tax,sep="\t",header=T,colClasses = "character")
  nc=parallel::detectCores()
    print("Converting taxon formats to metaxa format, can take a while...")
  lst=pbapply::pblapply(1:length(tab$ID),function(x){
    tmp=tab[x,!is.na(tab[x,])]
    return(data.frame(ID=tmp$ID,
                      tax=paste0(paste0(sapply(2:dim(tmp)[2],function(y) paste0(substr(colnames(tmp)[y],1,1),"__",tmp[,y])),collapse = ";"),";"),
                      stringsAsFactors = F))
  },cl=nc-2)
  lst=do.call(rbind,lst)
  write.table(lst,out,sep="\t",quote = F,row.names = F,col.names = F)
}

checkIDs <- function(metax=NULL,seq=NULL){
  tab=read.table(metax,sep="\t",colClasses = "character")
  colnames(tab)=c("ID","tax")
  dna=readDNAStringSet(seq)
  names(dna)=gsub(" ","",names(dna))
  if (length(intersect(tab$ID,names(dna)))==length(tab$ID) & length(intersect(tab$ID,names(dna)))==length(names(dna))){
  print("Perfect match between IDs. Proceed to database construction.")
  } else {
  print("IDs in taxonomy but not present in sequences:")
  print(setdiff(tab$ID,names(dna)))
  print("IDs in sequences but not present in taxonomy:")
  print(setdiff(names(dna),tab$ID))
  }
}

combine_files <- function(list=NULL,out=NULL){
  system(paste(c("cat",ls,">",out),collapse=" "))
}

make_db <- function(fasta=NULL,tax=NULL,out=NULL){

}
