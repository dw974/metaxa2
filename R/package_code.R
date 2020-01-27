library(DECIPHER)
library(ggplot2)
library(dplyr)

#' Convert from Metaxa taxonomy format to a tabular, easy-to-read format
#'
#' @param metax File name: the file you wish to convert
#' @param out File name: the output file. If this is set to "FALSE" no output file is written, and the tabular taxonomy format is returned.
#'
#' @return A data frame containing IDs and taxonomic classifications.
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
  },cl = 2)
  print("Concatenating dataframe.")
  lst=do.call(rbind,lst)
  if (out==F){
  } else {
    write.table(lst,out,sep="\t",quote = F,row.names = F)
  }
  return(lst)
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
  if (out==F){
  } else {
    write.table(lst,out,sep="\t",quote = F,row.names = F,col.names = F)
  }
  return(tab)
}


#'  Check matching IDs between metaxa2 taxonomy and fasta files
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
  system(paste0("mkdir ",outdir,"/",name))
  system(paste0(system.file("extdata", "metaxa2_dbb", package = "metaxa2")," -e ",fasta," -t ",tax," --cpu ",nc-2," --plus -o ",outdir,"/",name," -g ",name))

  ## PUTS DATABASE IN WEIRD LOCATION... NEED TO COPY IT OUT AGAIN.... SOMETHING LIKE
  ##system(paste0("cp -R ",system.file("extdata", "new_database", package = "metaxa2")," ",outdir))
}


#' Execute this function the first time you install the package
#'
setup_package <- function(){
  setwd(dirname(system.file("extdata", "metaxa2_dbb", package = "metaxa2")))
  system("find *m* -type f -exec chmod +x {} +")
}

#' Run metaxa2 on a dataset
#'
#' @param db location of the directory containing the metaxa2 database
#' @param input input query fasta file
#' @param out base name for output files (including folder location)
run_metaxa2<- function(db=NULL,input=NULL,out=NULL){
  nc=parallel::detectCores()
  system(paste0(system.file("extdata", "metaxa2", package = "metaxa2")," -o ",out," -f f -i ",input," --plus -g COI --cpu ",nc-2," -p ",db,"/HMMs/ -d ",db,"/blast"))
  system(paste0(system.file("extdata", "metaxa2_ttt", package = "metaxa2")," -o ",out,"_tax -i ",out,".taxonomy.txt"))
}

#' Create a fasta file containing only unique sequences from your data
#'
#' The function also outputs the file "correspondence.tab", which details the correspondence between original sequences and the defined unique sequences.
#' It also outputs the file "cross_table.tab", which summarises the number of each unique sequence in each file.
#'
#' @param list A character vector of fasta files that you want to include
#' @param outdir The folder location in which to write the results
unique_seqs <- function(list=NULL,outdir=NULL){
  #Files could be very big, so we'll start by writing a unique sequences file for each individual file.
  fls=""
  for (x in 1:length(list)){
    dna=Biostrings::readDNAStringSet(list[x])
    seqinr::write.fasta(sequences = as.list(paste(dna)),names=as.list(paste0("seq_",1:length(dna))),file.out = paste0(outdir,"/unique_",x,".fasta"))
    fls=c(fls,paste0(outdir,"/unique_",x,".fasta"))
  }
  #Then concatenate them into a single file.
  system(paste0("cat ",paste(fls,collapse = " ")," > ",outdir,"/all_sequences.fasta"))
  #delete temporary files
  system(paste0("rm ",paste(fls,collapse = " ")))
  dna=Biostrings::readDNAStringSet(paste0(outdir,"/all_sequences.fasta"))
  #Identify unique sequences and write to file
  unseq=unique(paste(dna))
  seqinr::write.fasta(sequences = as.list(unseq),names=as.list(paste0("seq_",1:length(unseq))),file.out = paste0(outdir,"/unique_sequences.fasta"))

  #find correspondence between unique sequences and original data
  cor=lapply(1:length(list),function(x){
    dna=Biostrings::readDNAStringSet(list[x])
    return(data.frame(file=list[x],ID=names(dna),seq=match(paste(dna),unseq),stringsAsFactors = F))
  })
  df=do.call(rbind,cor)
  tmp=as.matrix(table(df$file,df$seq))
  colnames(tmp)=as.list(paste0("seq_",1:length(unseq)))
  write.table(df,paste0(outdir,"/correspondence.tab"),col.names = T,row.names = F,sep="\t")
  write.table(tmp,paste0(outdir,"/cross_table.tab"),col.names=T,row.names = T,sep="\t")

}


summarise_db=function(tax=NULL){
  df=data.frame(taxonomic.level=c("species","genera","families","orders","classes","phyla","kingdoms"),
                count=c(length(unique(tax$species)),
                        length(unique(tax$genus)),
                        length(unique(tax$family)),
                        length(unique(tax$order)),
                        length(unique(tax$class)),
                        length(unique(tax$phylum)),
                        length(unique(tax$kingdom)))
  )
  df$taxonomic.level=factor(df$taxonomic.level,levels = rev(c("species","genera","families","orders","classes","phyla","kingdoms")))

  a=ggplot2::ggplot(df,aes(x=taxonomic.level,y=count,label=count,fill=taxonomic.level))+
    ggchicklet::geom_chicklet(aes(y=1),radius = grid::unit(20,"pt"))+
    geom_text(aes(y=0.5))+
    theme_void()+
    theme(axis.text.x = element_text(),axis.text.y = element_text(colour="white"),legend.position="")+
    ggtitle("Number of unique sequences assigned to different taxa")

  b=ggplot2::ggplot(df,aes(x=taxonomic.level,y=count,label=count,fill=taxonomic.level))+
    geom_col()+theme_void()+theme(axis.text = element_text(),legend.position="")

  c=ggplot(data.frame())+geom_blank()+theme_void()

  gridExtra::grid.arrange(c,c,c,c,a,c,c,b,c,c,c,c,nrow=4,ncol=3,heights=c(0.5,1,4,0.5),widths=c(0.2,2,0.2))

}

summarise_metaxa_results=function(results=NULL,correspondence=NULL){
  tab=read.table(results,header=F,stringsAsFactors = F,sep="\t")
  colnames(tab)=c("ID","tax","ident","length","RS")
  nc=parallel::detectCores()
  #Importing taxonomic results into dataframe format
  print("Importing taxonomic results into dataframe format")
  lst=lapply(1:length(tab$tax),function(x){
    str=stringr::str_split(tab$tax[x],";|__",simplify = T)
    return(data.frame(ID=tab$ID[x],
                      kingdom=ifelse(length(which(str=="k"))>0,str[which(str=="k")+1],NA),
                      phylum=ifelse(length(which(str=="p"))>0,str[which(str=="p")+1],NA),
                      class=ifelse(length(which(str=="c"))>0,str[which(str=="c")+1],NA),
                      order=ifelse(length(which(str=="o"))>0,str[which(str=="o")+1],NA),
                      family=ifelse(length(which(str=="f"))>0,str[which(str=="f")+1],NA),
                      genus=ifelse(length(which(str=="g"))>0,str[which(str=="g")+1],NA),
                      species=ifelse(length(which(str=="s"))>0,str[which(str=="s")+1],NA),
                      stringsAsFactors = F))
  })
  lst=do.call(rbind,lst)
  df=data.frame(taxonomic.level=factor(c("species","genus","family","order","class","phylum","kingdom","species","genus","family","order","class","phylum","kingdom"),levels=rev(c("species","genus","family","order","class","phylum","kingdom"))),
                count=c(sum(!is.na(lst$species)),
                        sum(!is.na(lst$genus)),
                        sum(!is.na(lst$family)),
                        sum(!is.na(lst$order)),
                        sum(!is.na(lst$class)),
                        sum(!is.na(lst$phylum)),
                        sum(!is.na(lst$kingdom)),
                        sum(is.na(lst$species)),
                        sum(is.na(lst$genus)),
                        sum(is.na(lst$family)),
                        sum(is.na(lst$order)),
                        sum(is.na(lst$class)),
                        sum(is.na(lst$phylum)),
                        sum(is.na(lst$kingdom))),
                type=factor(c("classified","classified","classified","classified","classified","classified","classified","unclassified","unclassified","unclassified","unclassified","unclassified","unclassified","unclassified")),levels=c("classified","unclassified"))
  a=ggplot2::ggplot(subset(df,df$type=="classified"),ggplot2::aes(x=taxonomic.level,y=count,label=count,fill=taxonomic.level))+
    ggchicklet::geom_chicklet(ggplot2::aes(y=1),radius = grid::unit(20,"pt"))+
    geom_text(aes(y=0.5))+
    theme_void()+
    theme(axis.text.x = element_text(),axis.text.y = element_text(colour="white"),legend.position="")+
    ggtitle("Classification summary: number and proportion of unique sequences classified to each taxonomic level.")
  b=ggplot2::ggplot(df,ggplot2::aes(x=taxonomic.level,y=count,label=count,fill=type))+
    geom_col(position="stack")+theme_void()+theme(axis.text = element_text(),legend.position="")+
    theme(legend.position="bottom")+
    scale_fill_manual(values = c("green","grey50"),name="")
  c=ggplot(data.frame())+geom_blank()+theme_void()
  gridExtra::grid.arrange(c,c,c,c,a,c,c,b,c,c,c,c,nrow=4,ncol=3,heights=c(0.5,1,4,0.5),widths=c(0.2,2,0.2))

  print("Crossing classifications with original data")
  if (!is.null(correspondence)){
    cor=as.matrix(read.table(correspondence))
    cor=cor[,!apply(cor,2,sum)==0]
    cor=as.data.frame(as.table(cor),stringsAsFactors = F)
    colnames(cor)=c("file","seq","count")
    cor=dplyr::left_join(cor,lst,by=c("seq"="ID"))
  }
  cor$taxon=sapply(1:dim(cor)[1],function(x) paste(cor$kingdom[x],cor$phylum[x],cor$class[x],cor$order[x],cor$family[x],cor$genus[x],cor$species[x],collapse = "::"))
  taxsum= cor %>% group_by(taxon) %>% dplyr::summarise(count=sum(as.numeric(count)))
  taxsum$taxon=factor(taxsum$taxon,levels=taxsum$taxon[order(taxsum$count,decreasing = F)])
  ggplot2::ggplot(taxsum[order(taxsum$count,decreasing = T)[1:50],])+geom_col(aes(x=taxon,y=count,fill=taxon))+coord_flip()+theme(legend.position = "none")

  ggplot(cor)+geom_col(aes(x=basename(file),y=count,fill=order))+coord_flip()

  xt=as.matrix(xtabs(count~taxon+file,data=cor))
  xt=t(xt)
  xt_f <- xt[,apply(xt, 2, var, na.rm=TRUE) != 0]
  pc=prcomp(xt_f,center=T,scale.=F)
  ggbiplot(pc, choices = 1:2,obs.scale = 1, var.scale = 1, varname.size = 3, labels.size=3, ellipse = TRUE, circle = TRUE,var.axes=F,labels=basename(rownames(pc$x)))

  return(xt)
}

