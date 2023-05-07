#Define custom functions for the script
reference_dir="~/Google Drive File Stream/My Drive/Reference_files/"
local_ref=paste0(reference_dir,"genome.fa")
lustre_ref="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

root_dir<-"~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
source(paste0(root_dir,"/Data/GT_functions.R"))
setwd(tree_mut_dir); source("treemut.R");setwd(root_dir)
vcf_header_path=paste0(reference_dir,"vcfHeader.txt")
plots_dir=paste0(root_dir,"/plots/")

#Manually enter the individual-level metadata data frame
exp_IDs<-c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
Individual_metadata=data.frame(ID=exp_IDs,
                               new_ID=c("SCD4","SCD6","SCD5","SCD3","SCD1","SCD2"),
                               Disease=rep("SCD",6),
                               Age_at_GT=c(20,26,24,16,7,13),
                               PD_no=c("PD49229","PD53373b","PD49228","PD53374b","PD49227","PD49226"),
                               CD34_dose=c(5.07,NA,5.15,8.26,4.86,3.55))

##------------------------------TREE PLOTTING------------------------------
sample_metadata<-read.delim(paste0(root_dir,"/Data/sample_metadata_full.tsv"),stringsAsFactors = F)
all_tree_data=readRDS(paste0(root_dir,"/Data/combined_tree_files.Rds"))
all_mut_data=readRDS(paste0(root_dir,"/Data/combined_muts_files.Rds"))

all.trees.uncorrected<-all_tree_data$all.trees.uncorrected
all.trees.cc<-all_tree_data$all.trees.cc
all.trees.cc.nodups<-all_tree_data$all.trees.cc.nodups
all.trees.ultra<-all_tree_data$all.trees.ultra

all.muts<-all_mut_data$all.muts
all.muts.nodups<-all_mut_data$all.muts.nodups

CN_change_df<-readr::read_csv(paste0(root_dir,"/Data/Copy_number_changes.csv"))
#Define colour schemes for the sample types and VCNs
scales::show_col(cell_type_colours)
cell_type_colours=c("gray75","#B3CDE3","#CCEBC5", "#9970AB")
names(cell_type_colours)<-c("Pre-GT","Pre-TDX","DP","Post-GT")

colScheme<-RColorBrewer::brewer.pal(9,"YlOrRd")[3:9]
VCN_colours=c("white",colorRampPalette(colScheme)(5))
names(VCN_colours)=c(0:4,">4")

#Create plots
pdf(paste0(plots_dir,"/Sickle_full_phylos.pdf"),width=15,6)
par(mfrow=c(1,1))
temp=Map(tree=all.trees.cc.nodups,details=all.muts.nodups,f=function(tree,details){
  print(details%>%
          dplyr::filter(coding_change_chip=="yes")%>%
          dplyr::select(mut_ref,variant_ID))
  
  #Select out just the cell types/ VCNs present in this sample set
  cell_type_colours<-cell_type_colours[which(names(cell_type_colours)%in%(sample_metadata%>%dplyr::filter(Sample%in%tree$tip.label)%>%pull(Sample_type)))]
  max_VCN=sample_metadata%>%dplyr::filter(Sample%in%tree$tip.label)%>%pull(VCN_rounded)%>%max()
  
  #if(!is.na(max_VCN)) {VCN_colours<-VCN_colours[1:(1+max_VCN)]}
  
  #Plot the tree
  tree=plot_tree(tree,cex.label=0,vspace.reserve=0.05)
  temp=plot_tree_labels(tree,
                        details,
                        type="line",
                        query.field = "Decision", #alternative is 'coding_change_chip'
                        data.frame(value=c("Oncogenic","Possible","Probable"),col="red",pch = 17,stringsAsFactors = FALSE), #if use 'coding_change_chip', value is 'Coding change mutation in driver'
                        label.field = "variant_ID",
                        cex.label = 0.8,
                        lty=2,
                        lwd=2)
  sapply(tree$tip.label,
         plot_category_tip_point,
         tree=tree,
         details=NULL,
         cat_df=sample_metadata%>%dplyr::filter(Sample%in%tree$tip.label)%>%dplyr::rename("sample"=Sample),
         cat_name="Sample_type",
         cols=cell_type_colours,
         cex=0.5,
         lwd=0.1
  )
  hm<-matrix(nrow=3,ncol=length(tree$tip.label),dimnames = list(c("Sample type","VCN","CNA"),tree$tip.label))
  for(i in 1:ncol(hm)){
    col1<-cell_type_colours[sample_metadata$Sample_type[sample_metadata$Sample==tree$tip.label[i]]]
    col2<-VCN_colours[as.character(sapply(sample_metadata$VCN_rounded[sample_metadata$Sample==tree$tip.label[i]],function(n) ifelse(n<5,n,"≥5")))]
    hm[1,i]<-col1[1]; hm[2,i]<-col2[1]
    for(i in 1:ncol(hm)){hm[3,i]<-ifelse(colnames(hm)[i]%in%CN_change_df$Sample,"purple","lightgrey")}
    hm[,"Ancestral"]<-"white"
  }
  legend("topleft", inset=c(.01,.01), title="Sample type",
         names(cell_type_colours), fill=cell_type_colours, horiz=F, cex=0.7)
  if(!is.na(max_VCN)){
    tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7)
    legend("topleft", inset=c(.01,0.2), title="VCN",
           names(VCN_colours), fill=VCN_colours, horiz=F, cex=0.7)
  } else {
    tree=add_heatmap(tree,heatmap=hm[c(1,3),,drop=F],cex.label = 0.7)
  }
})
dev.off()

##Plot the SCD3/ SCD4 trees with the post-GT samples highlighted
pdf(paste0(plots_dir,"Combined_postGT_highlighted_phylo.pdf"),width=14,6)
par(mfrow=c(1,1))
Map(tree=all.trees.cc.nodups,details=all.muts.nodups,f=function(tree,details){
  tree=plot_tree(tree,cex.label = 0,default_edge_color = "lightgray")
  add_annotation(tree=tree,
                 details=details,
                 matrices=NULL,
                 annot_function=plot_postGT_tree,
                 highlight="post",
                 sharing_cols=c("black","gray82"),
                 cat_df=sample_metadata)
  #cell_type_colours["Pre-GT"]<-"gray80"
  cat_df<-sample_metadata%>%dplyr::filter(Sample%in%tree$tip.label)%>%dplyr::rename("sample"=Sample)
  sapply(tree$tip.label,
         plot_category_tip_point,
         tree=tree,
         details=NULL,
         cat_df=cat_df,
         cat_name="Sample_type",
         cols=cell_type_colours,
         cex=0.5,
         lwd=0,
         col=NA
  )
  hm<-matrix(nrow=1,ncol=length(tree$tip.label),dimnames = list(c("Sample type"),tree$tip.label))
  for(i in 1:ncol(hm)){
    if(tree$tip.label[i]=="Ancestral"){
      col1<-"white"
    } else {
      col1<-cell_type_colours[sample_metadata$Sample_type[sample_metadata$Sample==tree$tip.label[i]]]
    }
    hm[1,i]<-col1[1]
  }
  tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7)
  legend("topleft", inset=c(.01,.01),
         title="Sample type",
         names(cell_type_colours)[names(cell_type_colours)%in%cat_df$Sample_type],
         fill=cell_type_colours[names(cell_type_colours)%in%cat_df$Sample_type],
         horiz=F,
         cex=0.7)
  
})
dev.off()

##---Plot the pre-GT trees only-----

all.trees.preGT<-lapply(all.trees.cc.nodups,function(tree){
  postGT_samples<-sample_metadata%>%filter(Sample%in%tree$tip.label & Time_point>0)%>%pull(Sample)
  if(length(postGT_samples)>0){
    tree<-drop.tip(tree,postGT_samples)
  }
  return(tree)
})

pdf(paste0(plots_dir,"Pre_GT_sickle_phylos.pdf"),width=30,height=5)
par(mfrow=c(1,6))
temp=Map(tree=all.trees.preGT,exp_ID=names(all.trees.preGT), function(tree,exp_ID) {
  tree$coords<-NULL
  tree=plot_tree(tree,cex.label=0,title = exp_ID)
})
dev.off()

lapply(all.trees.preGT,function(tree) length(tree$tip.label[!tree$tip.label=="Ancestral"]))

##---Plot the post-GT trees only-----

all.trees.postGT<-lapply(all.trees.cc.nodups,function(tree){
  preGT_samples<-sample_metadata%>%filter(Sample%in%tree$tip.label & Time_point==0)%>%pull(Sample)
  if(length(preGT_samples)>0){
    tree<-drop.tip(tree,preGT_samples)
  }
  return(tree)
})
par(mfrow=c(1,2))
temp=Map(tree=all.trees.postGT[1:2],exp_ID=names(all.trees.postGT)[1:2], function(tree,exp_ID) {
  tree$coords<-NULL
  tree=plot_tree(tree,cex.label=0,title = exp_ID)
})
lapply(all.trees.postGT,function(tree) length(tree$tip.label))

