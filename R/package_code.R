library(stringr)
library(dplyr)
library(pbapply)
library(parallel)
library(DECIPHER)

#' Convert from Metaxa taxonomy format to a tabular, easy-to-read format
#'
#' @param metax File name: the file you wish to convert
#' @param out File name: the output file
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


#' Convert from tabular taxonomy file to a metaxa2-format taxonomy file
#'
#' @param tax File name: the file you wish to convert
#' @param out File name: the output file
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

#' Check matching IDs between metaxa2 taxonomy and fasta files
#'
#' @param metax File name: the metaxa2 taxonomy file
#' @param seq The metaxa2 sequence file (fasta format)
checkIDs <- function(metax=NULL,seq=NULL){
  tab=read.table(metax,sep="\t",colClasses = "character")
  colnames(tab)=c("ID","tax")
  dna=Biostrings::readDNAStringSet(seq)
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

#' Combine multiple fasta or taxonomy files into one
#'
#' @param list A character vector of files that you want to combine
#' @param out the output combined file name
combine_files <- function(list=NULL,out=NULL){
  system(paste(c("cat",list,">",out),collapse=" "))
}

#' Make a metaxa2 database
#'
#' @param fasta The metaxa2 sequence data (fasta file)
#' @param tax The metaxa2 taxonomy data
#' @param outdir The directory to save the database into
#' @param name The name to give to the database
make_db <- function(fasta=NULL,tax=NULL,outdir=NULL,name=NULL){
  nc=parallel::detectCores()
  system(paste0(system.file("extdata", "metaxa2_dbb", package = "metaxa2")," -e ",fasta," -t ",tax," --cpu ",nc-2," --plus -g ",outdir,"/",name," -o ",name))
}

#' Execute this function the first time you install the package
#'
setup_package <- function(){
setwd(dirname(system.file("extdata", "metaxa2_dbb", package = "metaxa2")))
system("find *m* -type f -exec chmod +x {} +")
}
