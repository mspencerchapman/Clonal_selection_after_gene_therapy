#SCRIPT FOR ANALYSING THE VECTOR INTEGRATION SITES 

# Here is the bash code for generating the vector coverage
# for SAMPLE in $(ls *_subset.bam|cut -c1-15);do
# samtools index ${SAMPLE}_subset.bam
# done
# for SAMPLE in $(ls *_subset.bam|cut -c1-15);do
# samtools coverage -r VECTOR:1200-3000 ${SAMPLE}_subset.bam -H
# done > coverage.txt
# ls *_subset.bam|cut -c1-15 > completed_samples.txt
# paste completed_samples.txt coverage.txt > sample_coverage.txt

library(dplyr)
library(ggplot2)
library(tidyr)
library(Rsamtools)
library(Ckmeans.1d.dp)

exp_IDs<-c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
Individual_metadata=data.frame(ID=exp_IDs,
                               new_ID=c("SCD4","SCD6","SCD6","SCD3","SCD1","SCD2"),
                               Disease=rep("SCD",6),
                               Age_at_GT=c(20,26,24,16,7,13),
                               PD_no=c("PD49229","PD53373b","PD49228","PD53374b","PD49227","PD49226"))
#Define paths for the script
root_dir<-"~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
setwd(paste0(root_dir,"/Data/vector_paired_files/"))
paired_files=list.files(path = paste0(root_dir,"/Data/vector_paired_files/"),pattern=".*vector_paired.txt")
paired_list=lapply(paired_files,function(file) {if(file.info(file)$size!=0) {read.delim(file,header=F,stringsAsFactors = F)}})
names(paired_list)<-substr(paired_files,1,15)
#paired_list<-paired_list[-which(sapply(paired_list, is.null))]
paired_list<-lapply(paired_list,function(df) if(is.null(df)){return(NULL)}else{df<-df %>% mutate_at(c("V2","V3","V7"), as.character); df<-df[!is.na(df$V4),];return(df[,1:7])})
paired_list<-Map(f=function(df,Sample_ID) if(is.null(df)){return(NULL)}else{cbind(df,data.frame(Sample=rep(Sample_ID,nrow(df))))},df=paired_list,Sample_ID=names(paired_list))
  
all_reads=dplyr::bind_rows(paired_list)
# all_reads%>%
#   filter(!(V3==4&V4>124680000 &V4<124695000))%>% #15kb region
#   filter(!(V3==4&V4>187394000 &V4<187395000))%>% #1kb region
#   filter(!(V3==11&V4>61735000 &V4<61735500))%>% #0.5kb region
#   filter(!(V3==11&V4>5240000 &V4<5350000))%>% #1Mb region
#   filter(!(V3==2&V4>33141000 &V4<33142000))%>%
#   mutate(V3=factor(V3,levels=c(1:22,"X","Y","hs37d5")))%>%
#   #filter(V3=="15")%>%
#   ggplot(aes(x=Sample,y=V4,col=V5))+
#   geom_jitter(alpha=0.5)+
#   facet_grid(cols=vars(V3),scales="free_x",space = "free")+
#   theme(axis.text.x = element_blank(),axis.text.y = element_text(size=5))+
#   coord_flip()

n_read_threshold=4
VIS_list=lapply(paired_list,function(df) {
  
  if(is.null(df)){stop(return(NULL))}
  
  #First remove PCR duplicates
  df<-df[!is.na(df$V2)&df$V2!=""&!is.na(df$V5),]
  df$V2<-as.integer(df$V2)
  df<-df[!as.logical(sapply(df$V2,bamFlagAsBitMatrix)[11,]),]
  
  #Remove reads that are from erroneous mapping
  df<-df%>%
    filter(!(V3==4&V4>124680000 &V4<124695000))%>% #15kb region
    filter(!(V3==4&V4>187394000 &V4<187395000))%>% #1kb region
    filter(!(V3==11&V4>61735000 &V4<61735500))%>% #0.5kb region
    filter(!(V3==11&V4>5240000 &V4<5350000))%>% #1Mb region
    filter(!(V3==2&V4>33141000 &V4<33142000)) #1kb region
  
  
  #Now divide the list up by chromosome
  chroms=unique(df$V3)
  chroms=chroms[!is.na(chroms)&chroms!=""]
  out_by_chrom=lapply(chroms, function(chrom) {
    df_chrom<-df%>%filter(V3==chrom)
    if(nrow(df_chrom)<n_read_threshold){
      return(NULL)
    } else {
      result <- Ckmedian.1d.dp(df_chrom$V4,k=c(1,3))
      result_short=dplyr::bind_rows(result[c(2,4)])
      
      #Merge clusters that are a single VIS
      i=1
      while(i<nrow(result_short)) {
        if((result_short$centers[i+1]-result_short$centers[i])<1000) {
          result_short$size[i]<-result_short$size[i]+result_short$size[i+1]
          result_short$centers[i]<-round(mean(result_short$centers[i:(i+1)]))
          result_short<-result_short[-(i+1),]
        } else {
          i=i+1
        }
      }

      #Only include clusters with (i) adequate reads supporting it (ii) not at sites with known homology with vector
      pass_clusters=which(result_short$size>=n_read_threshold & (chrom!=11|abs(result_short$centers-5303264)>10000) & (chrom!=4|abs(result_short$centers-124683707)>10000))
      if(length(pass_clusters)>0) {
        VIS=paste(chrom,result_short$centers[pass_clusters],sep=":")
        return(VIS)  
      } else {
        return(NULL)
      }
      
    }
  })
  out_by_chrom<-unlist(out_by_chrom)  
  return(out_by_chrom)  
  })

VIS_df<-dplyr::bind_rows(Map(VIS=VIS_list,Sample=names(VIS_list),f = function(VIS,Sample){data.frame(Sample=rep(Sample,length(VIS)),Chrom_pos=VIS)}))%>%
  separate(Chrom_pos,into=c("Chrom","Pos"),sep=":")%>%
  mutate(Pos=as.integer(Pos))

VCN_by_VIS=data.frame(Sample=names(VIS_list),VCN_by_VIS=sapply(VIS_list,length))
write.table(VIS_df,file=paste0(root_dir,"/Data/VIS_from_WGS_df.txt"),row.names = F,quote=F,sep="\t")
write.table(VCN_by_VIS,file=paste0(root_dir,"/Data/VCN_by_VIS.txt"),row.names = F,quote=F,sep = "\t")

