## LOOK FOR ?HYDROXYCARBAMIDE SIGNATURE IN MPN PATIENTS
library(GenomicRanges)
library(IRanges)
library("Rsamtools")
library("MASS")
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(MutationalPatterns)
library(gridExtra)
library(ape)
options(stringsAsFactors = FALSE)

my_working_directory<-"/lustre/scratch119/casm/team154pc/ms56/lesion_segregation"


my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
resave_plots=F

lesion_seg_input_dir="/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data"
genome_file="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"

# lesion_seg_output_dir="~/R_work/Lesion_segregation/output/output_new"
# Phasing_MAV_file_path="~/R_work/Lesion_segregation/output/output_new/Phasing_results_MAVs_all"
# Phasing_PVV_file_path="~/R_work/Lesion_segregation/output/output_new/Phasing_results_PVVs_all"
# genome_file="~/Documents/Reference_files/genome.fa"
# output_dir="~/R_work/Lesion_segregation/"


#Set data file paths
data_sets=c("NW")
details_all_list=lapply(data_sets,function(data_set) {
  cat(data_set,sep="\n")
  samples<-readLines(paste0(lesion_seg_input_dir,"/",data_set,"_samples.txt"))
  details_list<-lapply(samples,function(SampleID) {
    cat(SampleID,sep="\n")
    info<-get_file_paths_and_project(dataset=data_set,Sample_ID=SampleID)
    load(info$filtered_muts_path)
    return(filtered_muts$COMB_mats.tree.build$mat)
  })
  names(details_list)<-samples
  return(details_list)
})
names(details_all_list)<-data_sets

Map(details=details_all_list$NW,Sample=names(details_all_list$NW),function(details,Sample) {
  write.vcf(details,vcf_path=paste0(lesion_seg_input_dir,"/NW/VCFs/",Sample,"_allmuts.vcf"),
            vcf_header_path = "/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs/mutation_vcfs/VCF_header_for_VaGrent.txt")
})

library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome,character.only=TRUE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd(paste0(lesion_seg_input_dir,"/NW/VCFs"))
vcf_files=list.files(".",pattern=".vcf")
vcfs<-MutationalPatterns::read_vcfs_as_granges(vcf_files,
                                               sample_names = stringr::str_split(vcf_files,pattern="_",simplify=T)[,1],
                                               genome=ref_genome)

mut_mat.all <- MutationalPatterns::mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

SNV_MutSig_plot<-MutationalPatterns::plot_96_profile(mut_matrix = mut_mat.all,ymax=0.08,condensed = T)
ggsave(filename="MPN_mutation_profiles.pdf",SNV_MutSig_plot,width=7,height=10)
N5_muts=c("T[T>A]A","T[T>A]T","T[T>G]A","T[T>G]T") #Define the specific contexts
colSums(mut_mat.all[N5_muts,])/colSums(mut_mat.all) #Define the proportions of mutations in these contexts



###Now do analysis of HDP output from the combined SCD/ MPN analysis
HDP_folder="~/R_work/Gene_therapy/Mutational_signature_extraction/HDP_040123/"
MPN_patients<-readLines(paste0(HDP_folder,"NW_samples.txt"))
SCD_trees=
MPN_trees=lapply(MPN_patients,function(patient) read.tree(paste0(HDP_folder,"MPN_trees/tree_",patient,".tree")))
names(MPN_trees)<-MPN_patients

mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
key_table=read.table(paste0(HDP_folder,"/key_table.txt"))

exposures_df=create_exposures_df(HDP_multi = mut_example_multi,trinuc_mut_mat = mutations,key_table=key_table)

#Get the signatures from individual samples - NB. this still only includes mutations on branches with >50 mutations
hdp_sigs<-dplyr::bind_rows(Map(tree=MPN_trees,Individual=MPN_patients,f=function(tree,Individual) {
  sigs=get_signatures_in_samples(tree = tree,signature_names=paste0("N",0:7),exposures_df = exposures_df%>%filter(Pair==Individual))
  return(sigs)
}))
cols=c(brewer.pal(12, "Paired"),"magenta","firebrick")

signature_cols<-RColorBrewer::brewer.pal(8,"Paired")
names(signature_cols)<-paste0("N",0:7)

MutSigsHDP_by_sample_abs_plot.preGT<-gather(hdp_sigs,-Sample,key="Signature",value="Prop")%>%
  mutate(ID=substr(Sample,1,6))%>%
  ggplot(aes(x=Sample,y=Prop,fill=Signature))+
  geom_bar(stat="identity",col="black",size=0.05)+
  scale_fill_manual(values=signature_cols)+
  facet_grid(cols=vars(ID),scale="free",space = "free")+
  theme(axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  my_theme

#Get the number of mutations at the specified contexts per sample
get_mut_numbers_at_contexts=function(ID,tree,contexts,trinuc_mut_mat,return_type="abs") {
  
  res<-sapply(tree$tip.label,function(Sample) {
    sample_node=which(tree$tip.label==Sample)
    nodes_included=get_ancestral_nodes(node = sample_node,edge = tree$edge,exclude_root = T)
    node_names=paste(nodes_included,ID,sep="_")
    node_names<-node_names[which(node_names%in%rownames(trinuc_mut_mat))]
    total_at_context=sum(trinuc_mut_mat[node_names,contexts])
    total=sum(trinuc_mut_mat[node_names,])
    if(return_type=="abs") {
      return(total_at_context)
    } else if(return_type=="rel") {
      return(total_at_context/total)
    }
  })
  df<-data.frame(ID=ID,Sample=tree$tip.label,n=res)
  return(df)
}

HC_exposure=readr::read_csv(paste0(HDP_folder,"HC_exposure.csv"))
N5_contexts=c("T.A.T.A","T.A.T.T","T.G.T.A","T.G.T.T")
mutations<-mutations[-which(rownames(mutations)=="162_PD5182"),]
temp=Map(tree=MPN_trees,Individual=MPN_patients,f=function(tree,Individual) {
  df<-get_mut_numbers_at_contexts(ID=Individual,tree=tree,contexts = N5_contexts,trinuc_mut_mat = mutations,return_type="abs")
})%>%dplyr::bind_rows()
head(temp)
temp%>%
  filter(Sample!="Ancestral")%>%
  left_join(HC_exposure)%>%
  ggplot(aes(x=forcats::fct_reorder(factor(ID),HC_exposure),y=n,label=paste(HC_exposure,"years\nof HC")))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(y=80,size=3)+
  geom_jitter(alpha=0.5,width=0.1,size=0.5)+
  theme_bw()+
  my_theme+
  labs(x="ID",y="Number of typical N5\nmutations per sample")

