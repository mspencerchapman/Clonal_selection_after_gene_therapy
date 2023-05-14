
#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","stringr","readr","tidyr","RColorBrewer","tibble","ape","phangorn","phytools","dichromat","seqinr","devtools","lmerTest")
bioconductor_packages=c("MutationalPatterns","BSgenome","BSgenome.Hsapiens.UCSC.hg19","TxDb.Hsapiens.UCSC.hg19.knownGene")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))


mut_sigs_theme=theme(strip.text.x = element_text(size=6,margin = margin(0.6,0,0.6,0, "mm")),
                     strip.text.y=element_text(size=6,margin = margin(0,0.6,0,0.6, "mm")),
                     #axis.text.x = element_blank(),
                     axis.text.x = element_text(size=3.5),
                     axis.text.y=element_text(size=5),
                     axis.title.x=element_text(size=7),
                     axis.title.y=element_text(size=7),
                     axis.ticks.x=element_blank(),
                     axis.ticks.y=element_line(linewidth=0.25),
                     legend.text = element_text(size=5),
                     legend.title = element_text(size=7),
                     legend.key.size=unit(2.5,"mm"),
                     strip.background = element_rect(linewidth=0.25),
                     panel.grid.major = element_line(linewidth=0.25),
                     panel.border = element_rect(linewidth=0.25))

#----------------------------------
# Set file paths and import data
#----------------------------------

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Define paths for the script
root_dir<-"~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
source(paste0(root_dir,"/Data/GT_functions.R"))
reference_dir="~/Google Drive File Stream/My Drive/Reference_files/"
vcf_header_path=paste0(reference_dir,"vcfHeader.txt")
plots_dir=paste0(root_dir,"/plots/")
local_ref=paste0(reference_dir,"genome.fa")
lustre_ref="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
vcf_dir=paste0(root_dir,"/Data/VCFs/")

#Manually enter the individual-level metadata data frame
exp_IDs<-c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
Individual_metadata=data.frame(ID=exp_IDs,
                               new_ID=c("SCD4","SCD6","SCD5","SCD3","SCD1","SCD2"),
                               Disease=rep("SCD",6),
                               Age_at_GT=c(20,26,24,16,7,13),
                               PD_no=c("PD49229","PD53373b","PD49228","PD53374b","PD49227","PD49226"))

setwd(root_dir)
sample_metadata<-read.delim(paste0(root_dir,"/Data/sample_metadata_full.tsv"),stringsAsFactors = F)
all_tree_data=readRDS(paste0(root_dir,"/Data/combined_tree_files.Rds"))
all_mut_data=readRDS(paste0(root_dir,"/Data/combined_muts_files.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(all_tree_data)) {assign(x,all_tree_data[[x]])}
for(x in names(all_mut_data)) {assign(x,all_mut_data[[x]])}

#Set the colour scheme for when colour by patient
patient_cols<-RColorBrewer::brewer.pal(7,"Set1")[c(1:5,7)]
names(patient_cols)<-paste0("SCD",1:6)

#----------------------------------
# Write VCFs of different mutation sets for reimport & analysis with MutationalPatterns
#----------------------------------

#Write vcfs of complete set of sample mutations
Map(details=all.muts.nodups,ID=names(all.trees),f=function(details,ID) if(nrow(details)>0){write.vcf(details,vcf_path = paste0(vcf_dir,ID,"_all.vcf"),vcf_header_path = vcf_header_path)})

#Now do vcfs of just the pre-GT mutations
Map(details=all.muts.nodups,tree=all.trees.cc.nodups,exp_ID=names(all.trees),f=function(details,tree,exp_ID) if(nrow(details)>0){
  pre_GT_samples=sample_metadata%>%dplyr::filter(ID==exp_ID&Time_point==0)%>%pull(Sample)
  if(length(pre_GT_samples)>0){
    pre_GT_nodes=unlist(sapply(tree$edge[,2],function(node) if(any(getTips(tree,node)%in%pre_GT_samples)){return(node)}else{return(NULL)}))
    write.vcf(details%>%dplyr::filter(node%in%pre_GT_nodes),vcf_path = paste0(vcf_dir,exp_ID,"_preGT.vcf"),vcf_header_path = vcf_header_path)
  }
})

#Now do vcfs of just developmental mutations (those with a molecular time <50)
Map(details=all.muts,tree=all.trees,exp_ID=names(all.trees),f=function(details,tree,exp_ID) if(nrow(details)>0){
  nh<-nodeHeights(tree); early_branches=tree$edge[,2][nh[,2]<50]
  write.vcf(details%>%dplyr::filter(node%in%early_branches),vcf_path = paste0(vcf_dir,exp_ID,"_early.vcf"),vcf_header_path = vcf_header_path)
})

#Look at the in vitro mutations signature - look at the private branches from confident duplicates
Map(details=all.muts,tree=all.trees,ID=names(all.trees),f=function(ID,tree,details) {
  df<-get_confident_duplicates(tree=tree,
                               sample_metadata=sample_metadata,
                               private_mut_threshold=150)
  duplicate_samples=df$Sample
  invitro_branches=which(tree$tip.label%in%unlist(duplicate_samples))
  write.vcf(details%>%filter(node%in%invitro_branches),vcf_path = paste0(vcf_dir,ID,"_invitro.vcf"),vcf_header_path = vcf_header_path)
})

#----------------------------------
# Re-import VCFs for SNV 96-mutational profile analysis
#----------------------------------

vcfs.all=MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_all.vcf",full.names = T),
                                                  sample_names = gsub("_all.vcf","",list.files(vcf_dir,pattern="_all.vcf",full.names = F)),
                                                  genome=ref_genome)
vcfs.pre=MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_preGT.vcf",full.names = T),
                                                  sample_names = gsub("_preGT.vcf","",list.files(vcf_dir,pattern="_preGT.vcf",full.names = F)),
                                                  genome=ref_genome)
vcfs.embryonic=MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_early.vcf",full.names = T),
                                                        sample_names = gsub("_early.vcf","",list.files(vcf_dir,pattern="_early.vcf",full.names = F)),
                                                        genome=ref_genome)
vcfs.invitro=MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_invitro.vcf",full.names = T),
                                                      sample_names = gsub("_invitro.vcf","",list.files(vcf_dir,pattern="_invitro.vcf",full.names = F)),
                                                      genome=ref_genome)

mut_mat.all <- MutationalPatterns::mut_matrix(vcf_list = vcfs.all, ref_genome = ref_genome)
mut_mat.pre <- MutationalPatterns::mut_matrix(vcf_list = vcfs.pre, ref_genome = ref_genome)
mut_mat.embryonic <- MutationalPatterns::mut_matrix(vcf_list = vcfs.embryonic, ref_genome = ref_genome)
mut_mat.invitro <- MutationalPatterns::mut_matrix(vcf_list = vcfs.invitro, ref_genome = ref_genome)
colnames(mut_mat.all)=colnames(mut_mat.pre)=colnames(mut_mat.embryonic)=colnames(mut_mat.invitro)<-sapply(gsub("_all.vcf","",list.files(vcf_dir,pattern="_all.vcf",full.names = F)),function(exp_ID) Individual_metadata$new_ID[Individual_metadata$ID==exp_ID]) #rename with the new IDs

SNV_MutSig_plot<-MutationalPatterns::plot_96_profile(mut_matrix = mut_mat.all[,c(5,6,4,1,3,2)],ymax=0.15,condensed = T)+ #Select just the SCD patients in order
  mut_sigs_theme+labs(x="Trinucleotide context",y="Relative contribution")
preGT_SNV_MutSig_plot<-MutationalPatterns::plot_96_profile(mut_matrix = mut_mat.pre[,c(5,6,4,1,3,2)],ymax=0.15,condensed = T)+ #Select just the SCD patients in order
  mut_sigs_theme+labs(x="Trinucleotide context",y="Relative contribution")
embryonic_SNV_MutSig_plot<-MutationalPatterns::plot_96_profile(mut_matrix = mut_mat.embryonic[,c(5,6,4,1,3,2)],ymax=0.15,condensed = T)+ #Select just the SCD patients in order
  mut_sigs_theme+labs(x="Trinucleotide context",y="Relative contribution")
invitro_SNV_MutSig_plot<-MutationalPatterns::plot_96_profile(mut_matrix = mut_mat.invitro[,c(5,6,4,1,3,2)],ymax=0.15,condensed = T)+ #Select just the SCD patients in order
  mut_sigs_theme+labs(x="Trinucleotide context",y="Relative contribution")

ggsave(filename = paste0(plots_dir,"SNV_MutSig_plot.pdf"),SNV_MutSig_plot,width=2.5,height=3)
ggsave(filename = paste0(plots_dir,"preGT_SNV_MutSig.pdf"),preGT_SNV_MutSig_plot,width=2.5,height=3)
ggsave(filename = paste0(plots_dir,"embryonic_SNV_MutSig.pdf"),embryonic_SNV_MutSig_plot,width=2.5,height=3)
ggsave(filename = paste0(plots_dir,"invitro_SNV_MutSig.pdf"),invitro_SNV_MutSig_plot,width=2.5,height=3)

#----------------------------------
# Re-import these VCFs for indel mutational profile analysis
#----------------------------------

indel_vcfs.all <- MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_all.vcf",full.names = T),
                                                           sample_names = gsub("_all.vcf","",list.files(vcf_dir,pattern="_all.vcf",full.names = F)),
                                                           genome=ref_genome,
                                                           type="indel")
indel_vcfs.pre <- MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_preGT.vcf",full.names = T),
                                                           sample_names = gsub("_preGT.vcf","",list.files(vcf_dir,pattern="_preGT.vcf",full.names = F)),
                                                           genome=ref_genome,
                                                           type="indel")
indel_vcfs.embryonic <- MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_early.vcf",full.names = T),
                                                        sample_names = gsub("_early.vcf","",list.files(vcf_dir,pattern="_early.vcf",full.names = F)),
                                                        genome=ref_genome,
                                                        type="indel")
indel_vcfs.invitro <- MutationalPatterns::read_vcfs_as_granges(vcf_files = list.files(vcf_dir,pattern="_invitro.vcf",full.names = T),
                                                      sample_names = gsub("_invitro.vcf","",list.files(vcf_dir,pattern="_invitro.vcf",full.names = F)),
                                                      genome=ref_genome,
                                                      type="indel")

indel_vcfs.all <- get_indel_context(indel_vcfs.all, ref_genome)
indel_vcfs.pre <- get_indel_context(indel_vcfs.pre, ref_genome)
indel_vcfs.embryonic <- get_indel_context(indel_vcfs.embryonic, ref_genome)
indel_vcfs.invitro <- get_indel_context(indel_vcfs.invitro, ref_genome)

indel_counts.all<- count_indel_contexts(indel_vcfs.all)
indel_counts.pre<- count_indel_contexts(indel_vcfs.pre)
indel_counts.embryonic<- count_indel_contexts(indel_vcfs.embryonic)
indel_counts.invitro<- count_indel_contexts(indel_vcfs.invitro)

colnames(indel_counts.all)=colnames(indel_counts.pre)=colnames(indel_counts.embryonic)=colnames(indel_counts.invitro)<-Individual_metadata$new_ID #rename with the new IDs

INDEL_MutSig_plot<-plot_indel_contexts(indel_counts.all[,c(5,6,4,1,3,2)], condensed = TRUE)+mut_sigs_theme
preGT_INDEL_MutSig_plot<-plot_indel_contexts(indel_counts.pre[,c(5,6,4,1,3,2)], condensed = TRUE)+mut_sigs_theme
preGT_combined_INDEL_MutSig_plot<-plot_indel_contexts(as.matrix(rowSums(indel_counts.pre)), condensed = TRUE)+mut_sigs_theme
embryonic_INDEL_MutSig_plot<-plot_indel_contexts(as.matrix(rowSums(indel_counts.embryonic)), condensed = TRUE)+mut_sigs_theme
invitro_INDEL_MutSig_plot<-plot_indel_contexts(as.matrix(rowSums(indel_counts.invitro)), condensed = TRUE)+mut_sigs_theme

ggsave(filename = paste0(plots_dir,"INDEL_MutSig_plot.pdf"),INDEL_MutSig_plot,width=7,height=3)
ggsave(filename = paste0(plots_dir,"preGT_INDEL_MutSig.pdf"),preGT_INDEL_MutSig_plot,width=7,height=3)
ggsave(filename = paste0(plots_dir,"preGT_combined_INDEL_MutSig.pdf"),preGT_combined_INDEL_MutSig_plot,width=7,height=2)
ggsave(filename = paste0(plots_dir,"embryonic_INDEL_MutSig.pdf"),embryonic_INDEL_MutSig_plot,width=7,height=2)
ggsave(filename = paste0(plots_dir,"invitro_INDEL_MutSig.pdf"),invitro_INDEL_MutSig_plot,width=7,height=2)

#----------------------------------
# Import the mutational signature extraction data (from HDP analysis)
#----------------------------------

HDP_folder=paste0(root_dir,"/Data/HDP")
mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
key_table=read.table(paste0(HDP_folder,"/key_table.txt"))

exposures_df=create_exposures_df(HDP_multi = mut_example_multi,trinuc_mut_mat = mutations,key_table=key_table)
sig_profiles=mut_mat_HDP_comp(mut_example_multi,plot=T)
extracted_signatures_real_plot<-plot_96_profile2(sig_profiles[,paste0("N",1:5)],condensed=T,ymax=0.15)+mut_sigs_theme+labs(x="Trinucleotide context")+theme(axis.text.x=element_text(size=3))
ggsave(filename = paste0(plots_dir,"extracted_signatures_real_plot.pdf"),extracted_signatures_real_plot,width=4,height=3)
extracted_signatures_artefact_plot<-plot_96_profile(sig_profiles[,paste0("N",c(0,6))],condensed=T,ymax=0.05)+mut_sigs_theme+labs(x="Trinucleotide context")
ggsave(filename = paste0(plots_dir,"extracted_signatures_artefact_plot.pdf"),extracted_signatures_artefact_plot,width=3,height=2)

#Get the signatures from individual samples - NB. this still only includes mutations on branches with >50 mutations
hdp_sigs<-dplyr::bind_rows(Map(tree=all.trees.cc,Individual=names(all.trees.uncorrected),f=function(tree,Individual) {
  sigs=get_signatures_in_samples(tree = tree,colnames(sig_profiles),exposures_df = exposures_df%>%filter(Pair==Individual))
  return(sigs)
}))

signature_cols<-RColorBrewer::brewer.pal(7,"Paired")
names(signature_cols)<-paste0("N",0:6)

MutSigsHDP_by_sample_abs_plot.preGT<-gather(hdp_sigs,-Sample,key="Signature",value="Prop")%>%
  left_join(sample_metadata%>%dplyr::select(Sample,sample_status,Time_point,new_ID,SNV_burden_tc))%>%
  mutate(Abs=Prop*SNV_burden_tc)%>%
  filter(Time_point==0 & sample_status=="PASS" & Signature%in%real_sigs)%>%
  ggplot(aes(x=Sample,y=Abs,fill=Signature))+
  geom_bar(stat="identity",col="black",size=0.05)+
  scale_fill_manual(values=signature_cols)+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme
ggsave(filename=paste0(plots_dir,"MutSigsHDP_by_sample_absolute_plot_preGT.pdf"),MutSigsHDP_by_sample_abs_plot.preGT,width=7,height=2)

MutSigsHDP_abs_contribution_per_sample<-gather(hdp_sigs,-Sample,key="Signature",value="Prop")%>%
  left_join(sample_metadata%>%dplyr::select(Sample,sample_status,Time_point,new_ID,SNV_burden_tc))%>%
  mutate(Abs=Prop*SNV_burden_tc)%>%
  filter(Time_point==0 & sample_status=="PASS")%>%
  ggplot(aes(x=new_ID,y=Abs,col=Signature))+
  geom_boxplot(col="black",outlier.shape = NA)+
  geom_jitter(alpha=0.2,size=0.2,height=0,width = 0.2)+
  scale_color_manual(values=signature_cols)+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  facet_grid(cols=vars(Signature))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  labs(col="Individual")
ggsave(filename=paste0(plots_dir,"MutSigsHDP_abs_contribution_per_sample.pdf"),MutSigsHDP_abs_contribution_per_sample,width=7,height=2)

#Alternatively re-fit signatures to mutation sets from samples using the 'fit_to_signatures' function
#Add in a 'pure embryonic' signature and a 'pure invitro' signature as these on short branches are not well represented in the HDP mutations sets
refit_signatures=cbind(sig_profiles,
                       rowSums(mut_mat.embryonic)/sum(mut_mat.embryonic),
                       rowSums(mut_mat.invitro)/sum(mut_mat.invitro))
colnames(refit_signatures)<-c(colnames(sig_profiles),"embryonic_sig","invitro_sig")
refit_signatures_plot<-plot_96_profile(refit_signatures,condensed=T,ymax=0.15)+mut_sigs_theme+labs(x="Trinucleotide context")
write.table(refit_signatures,file=paste0(root_dir,"/Data/refit_signatures.tsv"),quote=F,sep = "\t")

#Compare these to the known COSMIC signatures using the cosine similarity metric
pheatmap::pheatmap(MutationalPatterns::cos_sim_matrix(refit_signatures,get_known_signatures()),
                   cluster_rows = F,
                   cluster_cols = F)

#Now do the refitting step
sigs_refit<-Map(tree=all.trees.uncorrected, details=all.muts,exp_ID=names(all.trees.uncorrected),function (tree,details,exp_ID) {
  cat(exp_ID,sep="\n")
  tree_samples<-tree$tip.label[!tree$tip.label=="Ancestral"]
  exp_contributions<-lapply(tree_samples,function(SampleID) {
    sample_nodes=get_ancestral_nodes(node=which(tree$tip.label==SampleID),edge = tree$edge)
    sample_nodes<-sample_nodes[sample_nodes%in%details$node[details$Mut_type=="SNV"]] #Remove 'empty' nodes
    sample_node_labels=paste(sample_nodes,exp_ID,sep = "_")
    
    #Aggregate counts from these nodes in the 'mutations' matrix
    sample_trinuc_mat<-as.matrix(colSums(mutations[sample_node_labels,,drop=F]))
    sig_res<-fit_to_signatures(mut_matrix=sample_trinuc_mat,signatures = refit_signatures)
    contributions<-as.data.frame(t(sig_res$contribution/sum(sig_res$contribution)))%>%mutate(Sample=SampleID,.before="N0")
    return(contributions)
  })%>%dplyr::bind_rows()
  return(exp_contributions)
})%>%dplyr::bind_rows()

#Join the signatures data to the sample metadata df & multiply the sig contribution to the overall mutation burden to get
#the absolute contribution of each signature to each sample
sig_exposure_method="hdp"
if(sig_exposure_method=="hdp") {
  sample_metadata<-left_join(sample_metadata,hdp_sigs)
  for(sig in colnames(hdp_sigs)[-1]) {
    sample_metadata[[paste0(sig,"_abs")]]<-sample_metadata[[sig]]*sample_metadata[["SNV_burden_AR"]]
    sample_metadata[[paste0(sig,"_tc_abs")]]<-sample_metadata[[sig]]*sample_metadata[["SNV_burden_tc"]]
  }
  invitro_sigs=c("N0","N6")
  sample_metadata<-sample_metadata%>%
    mutate(SNV_burden_invitro_removed=SNV_burden_AR*(1-N0-N6),
           SNV_burden_tc_invitro_removed=SNV_burden_tc*(1-N0-N6))
  
  #Define which signatures are the invitro ones, and substract these from the mutation burden
  invitro_sigs_abs=paste(invitro_sigs,"abs",sep="_")
  real_sigs=colnames(sig_profiles)[!colnames(sig_profiles)%in%invitro_sigs]
  real_sigs_abs=paste(real_sigs,"abs",sep = "_")
  
} else if(sig_exposure_method=="refit"){
  sample_metadata<-left_join(sample_metadata,sigs_refit)
  for(sig in colnames(refit_signatures)) {
    sample_metadata[[paste0(sig,"_abs")]]<-sample_metadata[[sig]]*sample_metadata[["SNV_burden_AR"]]
    sample_metadata[[paste0(sig,"_tc_abs")]]<-sample_metadata[[sig]]*sample_metadata[["SNV_burden_tc"]]
  }
  invitro_sigs=c("invitro_sig","N0","N6")
  sample_metadata<-sample_metadata%>%
    mutate(SNV_burden_invitro_removed=SNV_burden_AR*(1-N0-N6-invitro_sig),
           SNV_burden_tc_invitro_removed=SNV_burden_tc*(1-N0-N6-invitro_sig))
  #Define which signatures are the invitro ones, and substract these from the mutation burden
  invitro_sigs_abs=paste(invitro_sigs,"abs",sep="_")
  real_sigs=colnames(refit_signatures)[!colnames(refit_signatures)%in%invitro_sigs]
  real_sigs_abs=paste(real_sigs,"abs",sep = "_")
}

invitro_sig_contribution_plot<-sample_metadata%>%
  filter(sample_status=="PASS" & !is.na(SNV_burden_AR))%>%
  dplyr::select(new_ID,Sample,all_of(invitro_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="nmuts")%>%
  ggplot(aes(x=Sample,y=nmuts,fill=Signature))+
  geom_bar(stat="identity",col="black",size=0.05)+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme

sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed) & Time_point==0)%>%
  dplyr::select(Sample,new_ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="Proportion")%>%
  dplyr::filter(Sample!="Ancestral")%>%
  ggplot(aes(x=Sample,y=Proportion,fill=factor(Signature,levels=real_sigs_abs)))+
  geom_bar(stat="identity",col="black",size=0.01,position = "fill")+
  scale_fill_manual(values = brewer.pal(8,"Paired")[-1])+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme+
  labs(fill="Signature")

#Plot signature contributions in each sample: relative & absolute
MutSigs_by_sample_relative_plot.preGT<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed) & Time_point==0)%>%
  dplyr::select(Sample,new_ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="Proportion")%>%
  dplyr::filter(Sample!="Ancestral")%>%
  ggplot(aes(x=Sample,y=Proportion,fill=factor(Signature,levels=real_sigs_abs)))+
  geom_bar(stat="identity",col="black",size=0.01,position = "fill")+
  scale_fill_manual(values = brewer.pal(8,"Paired")[-1])+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme+
  labs(fill="Signature")

ggsave(filename=paste0(plots_dir,"MutSigs_by_sample_relative_REFIT_preGT.pdf"),width=7,height=2)

MutSigs_by_sample_absolute_plot.preGT<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed)& Time_point==0)%>%
  dplyr::select(Sample,new_ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="Number_of_mutations")%>%
  mutate(Signature=gsub("_abs","",Signature))%>%
  dplyr::filter(Sample!="Ancestral")%>%
  ggplot(aes(x=Sample,y=Number_of_mutations,fill=factor(Signature,levels=real_sigs)))+
  geom_bar(stat="identity",col="black",size=0.01)+
  scale_fill_manual(values = brewer.pal(8,"Paired")[-1])+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  labs(x="Sample",y="Absolute signature contribution",fill="Signature")+
  my_theme

ggsave(filename=paste0(plots_dir,"MutSigs_by_sample_absolute_REFIT_preGT.pdf"),width=7,height=2)

#Plot signature contributions in each sample: relative & absolute
MutSigs_by_sample_relative_plot.postGT<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed) & Time_point>0)%>%
  dplyr::select(Sample,new_ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="Proportion")%>%
  dplyr::filter(Sample!="Ancestral")%>%
  ggplot(aes(x=Sample,y=Proportion,fill=factor(Signature,levels=real_sigs_abs)))+
  geom_bar(stat="identity",col="black",size=0.01,position = "fill")+
  scale_fill_manual(values = brewer.pal(8,"Paired")[-1])+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme+
  labs(fill="Signature")

ggsave(filename=paste0(plots_dir,"MutSigs_by_sample_relative_REFIT_postGT.pdf"),width=7,height=2)

MutSigs_by_sample_absolute_plot.postGT<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed)& Time_point>0)%>%
  dplyr::select(Sample,new_ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-new_ID,key = "Signature",value="Number_of_mutations")%>%
  mutate(Signature=gsub("_abs","",Signature))%>%
  dplyr::filter(Sample!="Ancestral")%>%
  ggplot(aes(x=Sample,y=Number_of_mutations,fill=factor(Signature,levels=real_sigs)))+
  geom_bar(stat="identity",col="black",size=0.01)+
  scale_fill_manual(values = brewer.pal(8,"Paired")[-1])+
  facet_grid(cols=vars(new_ID),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  labs(x="Sample",y="Absolute signature contribution",fill="Signature")+
  my_theme

ggsave(filename=paste0(plots_dir,"MutSigs_by_sample_absolute_REFIT_postGT.pdf"),width=7,height=2)

#----------------------------------
# Now look at strand bias of N2 (SBS19-like)
#----------------------------------

#The N2 signature is most prevalent in SCD2 (BCL009) & SCD3 (BCL006)
#The predominant feature is C>T mutations at CpT sites
#Therefore check for strand bias in this mutation set

mut_mat_s <- mut_matrix_stranded(vcfs.pre, ref_genome, genes_hg19)
mut_mat_s[grep("[C>T]T",rownames(mut_mat_s),fixed=T),c("BCL006","BCL009"),drop=F]
poisson.test(sum(mut_mat_s[grep("[C>T]T-tr",rownames(mut_mat_s),fixed=T),c("BCL006","BCL009")]),
             sum(mut_mat_s[grep("[C>T]T-un",rownames(mut_mat_s),fixed=T),c("BCL006","BCL009")]))

#Reveals bias for C on the transcribed strand -> G is the lesion-containing base

plot_192_profile(mut_mat_s[,4,drop=F])
strand_counts <- strand_occurrences(mut_mat_s,by=colnames(mut_mat_s))
plot_strand(strand_counts)
strand_bias <- strand_bias_test(strand_counts)
strand_bias%>%filter(type=="C>T")

#----------------------------------
# Now look at strand bias of N5 (novel signature associated with hydroxycarbamide)
#----------------------------------

#Reveals bias for T on the untranscribed strand -> T is the lesion-containing base
mut_mat_s[grep("T\\[T>[A,G]\\][A,T]",rownames(mut_mat_s)),,drop=F]
poisson.test(sum(mut_mat_s[grep("T\\[T>[A,G]\\][A,T]-tr",rownames(mut_mat_s)),c("BCL006","BCL004")]),
             sum(mut_mat_s[grep("T\\[T>[A,G]\\][A,T]-un",rownames(mut_mat_s)),c("BCL006","BCL004")]))

plot_192_profile(mut_mat_s[,4,drop=F])
strand_counts <- strand_occurrences(mut_mat_s,by=colnames(mut_mat_s))
plot_strand(strand_counts)
strand_bias <- strand_bias_test(strand_counts)
strand_bias%>%filter(type=="T>G")

#----------------------------------
# Look at correlation between N5 signature contributions and hydroxycarbamide exposure
#----------------------------------

#Read in the data relating to hydroxycarbamide treatment
HC_exposure=readxl::read_excel(path=paste0(root_dir,"/Data/HU_treatment_record.xlsx"))

#Create table of average signature contributions per sample, at each time point in each individual
average_signature_contributions<-sample_metadata%>%
  dplyr::select(new_ID,Sample,Time_point,all_of(real_sigs_abs))%>%
  gather(-new_ID,-Sample,-Time_point,key="Signature",value="Absolute_contribution")%>%
  mutate(Signature=factor(Signature,levels=real_sigs_abs))%>%
  group_by(new_ID,Time_point,Signature)%>%
  dplyr::summarise(mean_sig=mean(Absolute_contribution,na.rm=T))

N5_HDP_df<-average_signature_contributions%>%
  dplyr::filter(Time_point==0 & Signature=="N5_abs")%>%
  left_join(HC_exposure)%>%left_join(Individual_metadata)

N5_HDP_plot<-N5_HDP_df%>%
  ggplot(aes(x=years,y=mean_sig,col=new_ID))+
  geom_point(size=0.6)+
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))+
  theme_classic()+
  my_theme+
  labs(x="Hydroxycarbamide exposure\n(years of administration)",y="N5 mutation burden",col="Individual")+
  theme(legend.position = "none")
ggsave(filename=paste0(plots_dir,"N5_HDP_vs_HC_plot.pdf"),plot = N5_HDP_plot,height=2,width=1.8)

#Define the correlation
lm.N5_HDP<-lm(mean_sig~Age_at_GT+years,data=N5_HDP_df)
summary(lm.N5_HDP)
coef(lm.N5_HDP)
confint(lm.N5_HDP)

#----------------------------------
# Regression of AVERAGE numbers of distinctive N5 T>A/T>G mutations vs hydroxycarbamide exposure
#----------------------------------

#The N5 signature is likely partially contaminated by other processes (particularly N2 [SBS19]).
#Therefore a cleaner measure of the underlying process is to look only at the T>A/ T>G mutations
#that are characteristic of N5, and very rare in normal ageing.

#Get summary table of the mean mutation burden in different samples/ time points
mut_burden_summary<-sample_metadata%>%
  filter(sample_status=="PASS")%>%
  group_by(new_ID,Time_point)%>%
  dplyr::summarise(mean_SNV_burden=mean(SNV_burden_tc))

N5_muts=c("T[T>A]A","T[T>A]T","T[T>G]A","T[T>G]T") #Define the specific contexts
colSums(mut_mat.pre[N5_muts,])/colSums(mut_mat.pre) #Define the proportions of mutations in these contexts

#Get the average number of these mutations by doing:
#Overall proportion of mutations at these contexts at the pre-GT time point x average mutation burden
n_N5_muts<-(colSums(mut_mat.pre[N5_muts,])/colSums(mut_mat.pre))[c(paste0("SCD",1:6))]*mut_burden_summary%>%filter(Time_point==0)%>%pull(mean_SNV_burden)

#Convert into dataframe including the hydroxycarbamide exposure and other individual metadata
N5_df<-data.frame(new_ID=names(n_N5_muts),N5_mut_burden=n_N5_muts)%>%left_join(HC_exposure)%>%left_join(Individual_metadata)
N5_characteristic_muts_plot<-ggplot(N5_df,aes(x=years,y=N5_mut_burden,col=new_ID))+
  geom_point(size=0.6)+
  scale_color_brewer(palette = "Paired")+
  #geom_smooth(method="lm",col="black",se = F,size=0.5)+
  scale_y_continuous(limits=c(0,25))+
  theme_classic()+
  my_theme+
  labs(x="Hydroxycarbamide exposure\n(years of administration)",y="T>A / T>G at TTT / TTA mutations \n(Mean numbers per colony)",col="Individual")

ggsave(filename=paste0(plots_dir,"N5_characteristic_muts_plot.pdf"),plot = N5_characteristic_muts_plot,height=2,width=2.5)

lm.N5<-lm(N5_mut_burden~Age_at_GT+years,data=N5_df)
summary(lm.N5)
coef(lm.N5)
confint(lm.N5)

N5_df%>%
  mutate(age_component=Age_at_GT*lm.N5$coefficients[2],HC_component=HC_exposure*lm.N5$coefficients[3])

N5_sig_by_HC_plot<-N5_df%>%
  ggplot(aes(x=years,y=N5_mut_burden,col=new_ID))+
  geom_point(alpha=1)+
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(limits=c(0,25))+
  geom_smooth(method="lm",col="black",se = F,size=0.5)+
  theme_classic()+
  my_theme+
  labs(x="Hydroxycarbamide exposure\n(years of administration)",y="Mean N5 mutations per colony",col="Individual")

ggsave(filename=paste0(plots_dir,"N5_muts_by_HC_plot.pdf"),plot = N5_sig_by_HC_plot,height=2,width=2.5)

#----------------------------------
# 'PER SAMPLE' regression of numbers of distinctive N5 T>A/T>G mutations vs hydroxycarbamide exposure
#----------------------------------

#Previous regression looked at average numbers per sample
#However, alternatively can do a similar regression including each individual sample
#Need to use linear mixed effects modeling using 'individual' as a random effect

#1. Get number of N5 mutations per sample, and correct for sample sensitivity
colnames(mutations)<-rownames(mut_mat.all)
N5_muts_per_sample<-Map(tree=all.trees.uncorrected, details=all.muts,exp_ID=names(all.trees.uncorrected),function (tree,details,exp_ID) {
  cat(exp_ID,sep="\n")
  tree_samples<-tree$tip.label[!tree$tip.label=="Ancestral"]
  N5_context_contributions<-sapply(tree_samples,function(SampleID) {
    sample_nodes=get_ancestral_nodes(node=which(tree$tip.label==SampleID),edge = tree$edge)
    sample_nodes<-sample_nodes[sample_nodes%in%details$node[details$Mut_type=="SNV"]] #Remove 'empty' nodes
    sample_node_labels=paste(sample_nodes,exp_ID,sep = "_")
    #Aggregate counts from these nodes in the 'mutations' matrix
    sample_N5_muts<-sum(mutations[sample_node_labels,N5_muts,drop=F])
    return(sample_N5_muts)
  })
  df=data.frame(Sample=tree_samples,N5_muts=N5_context_contributions)
  return(df)
})%>%dplyr::bind_rows()

N5_muts_per_sample_df<-N5_muts_per_sample%>%
  left_join(sample_metadata%>%dplyr::select(ID,Sample,Age_at_GT,Time_point,somatic_SNV_sensitivity))%>%
  left_join(HC_exposure)%>%
  mutate(N5_muts_corrected=N5_muts/somatic_SNV_sensitivity)

N5_muts_per_sample_df%>%
  filter(Time_point==0)%>%
  ggplot(aes(x=years,y=N5_muts_corrected,col=new_ID))+
  geom_jitter(height=0,width=0.2,alpha=0.1,size=0.2)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme_bw()+
  my_theme

lme.res<-lme4::lmer(N5_muts_corrected~Age_at_GT+years+(1|new_ID),data=N5_muts_per_sample_df)
confint(lme.res)

#----------------------------------
# Review signature contributions with age (in pre-GT samples)
#----------------------------------

Ages=Individual_metadata$Age_at_GT;names(Ages)<-Individual_metadata$ID
sample_metadata%>%
  mutate(Age_at_GT=Ages[ID])%>%
  filter(Time_point==0)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  dplyr::select(Sample,Age_at_sampling,ID,all_of(real_sigs_abs))%>%
  gather(-Sample,-Age_at_sampling,-ID,key="Signature",value="Number_of_mutations")%>%
  ggplot(aes(x=Age_at_sampling,y=Number_of_mutations,col=ID))+
  geom_jitter(height=0,width=0.15,alpha=0.3,size=0.2)+
  facet_wrap(~Signature,scales="free")+
  scale_x_continuous(limits=c(0,25))+
  geom_smooth(aes(x=Age_at_sampling,y=Number_of_mutations),method="lm",inherit.aes = F,col="black",size=0.5)+
  theme_bw()

N1_with_age.lm<-lm(N1_abs~Age_at_sampling,data=sample_metadata%>%
     mutate(Age_at_GT=Ages[ID])%>%filter(Time_point==0)%>%
     mutate(Age_at_sampling=(Time_point/12+Age_at_GT)))
N1_with_age.lm
confint(N1_with_age.lm)

#Review which signatures are most increased post-GT
sample_metadata%>%
  filter(!is.na(SNV_burden_invitro_removed))%>%
  mutate(Age_at_GT=Ages[ID])%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  ggplot(aes(x=factor(Time_point),y=embryonic_sig_abs))+
  geom_boxplot()+
  facet_grid(~ID,scales="free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

#----------------------------------
# IN VITRO MUTATIONAL SIGNATURE
#----------------------------------

#This section attempts to:
#(1) Define an 'in vitro signature' representing the mutational profile of mutations acquired during in vitro culture
#(2) Use this as a reference to help determine whether other closely related colonies are replicates of the same colony
#    or samples deriving from different single cells that are closely related in vivo.

#The approach is to assess the mutational profile of private mutations from the two closely related samples
#If they are replicates of the same colony, these should all be mutations acquired during in vitro cloning

#Spot the potential duplicates for testing signatures
#Find sample pairs where both have fewer than 150 private mutations from MRCA
potential.dup.df=Map(tree=all.trees.uncorrected,details=all.muts,function(tree,details) {
  #Determine the sets of duplicate samples
  df<-get_duplicate_status(tree=tree,
                           sample_metadata=sample_metadata,
                           height_cut_off =150)  
  return(df)
})

#Get the trinucleotide mutation matrix for private branches of all potential duplicate samples
invitro_trinuc_mat<-Map(dup.samples=potential.dup.df,details=all.muts,tree=all.trees.uncorrected,exp_ID=names(all.trees.uncorrected),function(dup.samples,details,tree,exp_ID) {
  dup_nodes=which(tree$tip.label%in%dup.samples$Sample)
  dup_nodes<-dup_nodes[dup_nodes%in%details$node[details$Mut_type=="SNV"]]
  return(mutations[paste(dup_nodes,exp_ID,sep="_"),])
})%>%bind_rows()%>%
  as.matrix()%>%
  t()

#Rename the rows with the standard 96-profile labels
rownames(invitro_trinuc_mat)<-rownames(mut_mat.all)

#Define the "duplicate recognition signatures" - this is the standard HSPC sig ("N1") and the invitro signature
duplicate_recognition_sigs<-refit_signatures[,c("N1","invitro_sig")]
colnames(duplicate_recognition_sigs)<-c("HSPC signature","In vitro signature") #Rename these

#Review the 1st ten 96-profile signature
MutationalPatterns::plot_96_profile(mut_matrix = invitro_trinuc_mat[,1:10],ymax=0.15,condensed = T)
duplicate_recognition_sigs_plot<-MutationalPatterns::plot_96_profile(mut_matrix =duplicate_recognition_sigs,ymax=0.08)+mut_sigs_theme
ggsave(filename = paste0(plots_dir,"duplicate_recognition_sigs.pdf"),duplicate_recognition_sigs_plot,width=5,height=2.5)  

#Now fit all the private branches to the optimal combination of the HSPC & invitro signatures
sig_att<-fit_to_signatures(invitro_trinuc_mat,signatures=duplicate_recognition_sigs)

#Get a tidy table of absolute contributions of the HSPC/ invitro signatures to each sample, together with their "duplicate status"
tb <- sig_att$contribution %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Signature") %>% 
  tidyr::pivot_longer(-Signature, names_to = "Sample", 
                      values_to = "Contribution")%>%
  dplyr::mutate(ID=stringr::str_split(Sample,"_",simplify=T)[,2])%>%
  left_join(Individual_metadata,by="ID")%>%
  dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)), Signature = factor(Signature,levels = unique(Signature)))%>%
  dplyr::mutate(node=stringr::str_split(Sample,pattern="_",simplify=T)[,1])%>%
  mutate(Sample=sapply(Sample,function(x){
    node<-as.numeric(stringr::str_split(x,pattern="_",simplify=T)[,1])
    ID<-stringr::str_split(x,pattern="_",simplify=T)[,2]
    return(all.trees.uncorrected[[ID]]$tip.label[node])
    }))%>%left_join(dplyr::bind_rows(potential.dup.df))

#
present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% dplyr::pull(Signature) %>% unique()
invitro_mutation_plot<-tb%>%
  mutate(status=str_wrap(status,width=10))%>%
  ggplot(aes(x = forcats::fct_reorder(factor(Sample),dup_set), y = Contribution, fill = forcats::fct_rev(Signature))) + 
  geom_bar(stat = "identity", colour = "black",size=0) + labs(x = "Private branches with <100 mutations", y = "Absolute contribution \n (no. mutations)") + scale_fill_discrete(breaks = present_sigs) + 
  theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                     panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
                     panel.grid.major.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(cols=vars(factor(new_ID)),rows=vars(status),scales = "free",space = "free_x")+
  #facet_wrap(new_ID~status,scales = "free_x")+
  geom_hline(yintercept=15,linetype=2,size=0.5)+
  my_theme+
  labs(fill="Signature")+
  theme(legend.position="bottom")

ggsave(filename = paste0(plots_dir,"invitro_mutation_contributions.pdf"),invitro_mutation_plot,width=7,height=3)  

SCD4_plate_map<-plot_duplicate_plate_map(all.trees.uncorrected$BCL002,
                                         sample_metadata%>%
                                           mutate(plate_ID=str_replace(plate_ID,pattern="_BCL","\nBCL")),
                                         height_cut_off =150,
                                         labels="set")
ggsave(filename = paste0(plots_dir,"SCD4_plate_map.pdf"),SCD4_plate_map,width=7,height=7)  

p<-tb%>%
  left_join(sample_metadata%>%dplyr::select(Sample,Time_point))%>%
  filter(status=="Not all duplicates")%>%
  mutate(dup_set=paste(new_ID,dup_set,sep="_"))%>%
  ggplot(aes(x=Sample, y = Contribution, col=Time_point>0,fill = forcats::fct_rev(Signature))) + 
  geom_bar(stat = "identity",size=0.5) + 
  labs(x = "Private branches with <100 mutations", y = "Absolute contribution \n (no. mutations)")+
  scale_fill_discrete(breaks = present_sigs)+
  scale_color_manual(values=brewer.pal(3,"Dark2"))+
  theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                     panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
                     panel.grid.major.y = element_blank())+
  theme(axis.text.x=element_text(angle=90,size=3),axis.ticks.x=element_blank(),strip.text.x = element_text(angle=90))+
  facet_grid(cols=vars(dup_set),scales = "free",space="free")+
  #facet_wrap(~dup_set,scales = "free",nrow=4)+
  geom_hline(yintercept=15,linetype=2)+
  my_theme

potential.dup.df$BCL008%>%filter(dup_set==1)%>%View()

sig_att$contribution[,sig_att$contribution[1,]>15]
BCL002_branches=as.numeric(gsub("BCL002_","",colnames(sig_att$contribution[,sig_att$contribution[1,]>15])))

sample_metadata%>%
  filter(Sample %in% all.trees[[1]]$tip.label[BCL002_branches])

#Define mean numbers of 'in vitro muts', the number
hist(colSums(sig_att$contribution[,sig_att$contribution[1,]<=15]))
mean(colSums(sig_att$contribution[,sig_att$contribution[1,]<=15]))
invitro_mut_summary<-t(sig_att$contribution[,sig_att$contribution[1,]<=15])%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = "ID_node")%>%
  mutate(ID=stringr::str_split(ID_node,pattern="_",simplify=T)[,1],total_muts=`Bone marrow signature`+`In vitro signature`)%>%
  left_join(Individual_metadata,by="ID")%>%
  group_by(new_ID)%>%
  summarise(mean_invitro=mean(total_muts))

invitro_mut_summary%>%
  ggplot(aes(x=new_ID,y=mean_invitro))+
  geom_bar(stat="identity")+
  theme_classic()+
  my_theme+
  labs(x="Individual",y="Mean number of in vitro mutations\n(duplicate colonies)")

#Number of invitro mutations from colony duplicates
t(sig_att$contribution[,sig_att$contribution[1,]<=15])%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = "ID_node")%>%
  mutate(ID=stringr::str_split(ID_node,pattern="_",simplify=T)[,1],node=stringr::str_split(ID_node,pattern="_",simplify=T)[,2],total_muts=`Bone marrow signature`+`In vitro signature`)%>%
  mutate(Sample=mapply(FUN=function(exp_ID,tree_node){all.trees[[exp_ID]]$tip.label[tree_node]},exp_ID=ID,tree_node=as.numeric(node)))%>%
  left_join(sample_metadata,by="Sample")%>%
  ggplot(aes(x=new_ID,y=total_muts,col=Cell_type))+
  geom_boxplot(aes(group=new_ID,x=new_ID,y=total_muts),outlier.shape = NA,inherit.aes = F)+
  geom_jitter(height=0,width = 0.15,alpha=0.5,size=1)+
  #geom_point(alpha=0.5)+
  scale_y_continuous(breaks=seq(0,80,10))+
  labs(x="",y="Numbers of invitro mutations\nfrom colony duplicates",col="Colony type")+
  theme_classic()+
  my_theme

saveRDS(sample_metadata,file=paste0(root_dir,"/Data/sample_metadata_full_with_sigs.Rds"))
