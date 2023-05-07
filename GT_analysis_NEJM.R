#Define custom functions for the script
reference_dir="~/Google Drive/My Drive/Reference_files/"
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

#Set the patient colour scheme
patient_cols<-RColorBrewer::brewer.pal(7,"Set1")[c(1:5,7)]
names(patient_cols)<-paste0("SCD",1:6)

tree_folder=paste0(root_dir,"/Data/tree_files_with_dups/")
annotated_muts_folder=paste0(root_dir,"/Data/annot_files_with_dups/")

##1. Import the SNV trees,details matrices & sample metadata
tree_paths=list.files(tree_folder,pattern=".tree",full.names = T)
all.trees<-lapply(exp_IDs,function(exp_no){read.tree(grep(exp_no,tree_paths,value = T))})
names(all.trees)<-exp_IDs

annotated_muts_paths=list.files(annotated_muts_folder,pattern="post_mix",full.names = T)
all.muts<-lapply(exp_IDs,function(exp_no){load(grep(exp_no,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$mat)})
names(all.muts)<-exp_IDs

genotype.mats<-lapply(exp_IDs,function(exp_no){load(grep(exp_no,annotated_muts_paths,value = T));return(filtered_muts$COMB_mats.tree.build$Genotype_bin)})
names(genotype.mats)<-exp_IDs

sample_metadata<-read.delim(paste0(root_dir,"/Data/Sample_level_metadata.tsv"))

##Join the individual level metadata to the sample level metadata
sample_metadata<-left_join(sample_metadata,Individual_metadata,by="ID")
sample_metadata<-mutate(sample_metadata,Time_point=Time_point/12) #Convert the 'time point' measure to years


##2. Re-annotate drivers using the union of several clonal haem/ myeloid malignancy/ haem malignancy gene lists
myeloid_panel_bolton=read.csv(paste0(reference_dir,"Myeloid_drivers_bolton.csv"),header = F)$V1
myeloid_panel_boston=read.csv(paste0(reference_dir,"Boston_myeloid_drivers.csv"),header = F)$V1
myeloid_panel_illumina=read.csv(paste0(reference_dir,"Illumina_trusight_myeloid_panel.csv"),header = F)$V1
all_ukbb_CH=c("DNMT3A","ASXL1","PPM1D","TET2","ZNF318","SRSF2","TP53","MTA2","YLPM1","ZBTB33",
              "SF3B1","MYD88","SIK3","GNB1","BRCC3","GNAS","SRSF1","CBL","JAK2","DUSP22","CHEK2","ZNF234","SPRED2","IGLL5","BAX",
              "CCDC115","MAGEC3","SH2B3","CCL22","PHIP","KDM6A","SRCAP")
driver_list<-unique(c(myeloid_panel_bolton,myeloid_panel_boston,myeloid_panel_illumina,all_ukbb_CH))
VennDiagram::venn.diagram(x=list(myeloid_panel_bolton,myeloid_panel_boston,myeloid_panel_illumina,all_ukbb_CH),
                          category.names = c("Bolton et al.\nmyeloid panel",
                                             "Boston haem.\nmalig. panel"
                                             ,"Illumina TruSight\nmyeloid panel",
                                             "UK biobank\ndNdS analysis"),
                          filename = paste0(plots_dir,"myeloid_panel_venn.png"),
                          output=T)

all.muts<-lapply(all.muts,function(details) {
  details$coding_change <- ifelse(details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                                    "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                                    "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                                    "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                                    "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                                    "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss")|
                                                                   grepl("splice_site_variant",details$Type),
                                                                 "Coding change",
                                                                 "no")
  
  details$coding_change_chip<-ifelse(details$coding_change=="Coding change" & details$Gene%in%driver_list,"yes","no")
  details$ChromPos=paste(details$Chrom,details$Pos,sep="-")
  details$variant_ID=paste(details$Gene, details$Protein, sep = " ")
  
  return(details)
})

#Write the possible driver mutations into a spreadsheet for annotation
annotated_driver_mut_set=read_csv(paste0(root_dir,"/Data/Possible_drivers_annotated.csv"))#[,c("mut_ref","Gene","variant_ID","Decision")]

Map(details=all.muts,exp_ID=names(all.muts),f=function(details,exp_ID) details%>%dplyr::filter(coding_change_chip=="yes")%>%mutate(ID=exp_ID)%>%dplyr::select(ID,mut_ref,Gene,variant_ID))%>%
   dplyr::bind_rows()%>%
   left_join(annotated_driver_mut_set)%>%
   write_csv(file=paste0(root_dir,"/Data/Possible_drivers.csv"))
#Read back in the annotated spreadsheet

#3. Update the details matrices with the manual decisions on whether mutations are oncogenic/ possible oncogenic or VUS
all.muts<-lapply(all.muts,function(details) {
  details$Decision<-NULL
  details<-left_join(details,annotated_driver_mut_set%>%dplyr::select(mut_ref,Decision),by="mut_ref")
  return(details)
})

#Spot the duplicates for dropping

##FIRST DEFINE THE SAMPLE SETS THAT ARE TRULY IN VIVO RELATED (not technical duplicates of the same colony)
#This is decided by a fairly manual process that includes:
#(1) Examining plate positions
#(2) Looking at mutaton burdens
#(3) Seeing if the signatures are primarily the HSC signature or in vitro signature
true_invivo_related_sets=list(c("PD49227b_lo0569","PD49227b_lo0650"),
                              c("PD49226b_lo0764","PD49226b_lo0774"),
                              c("PD49226b_lo0738","PD49226b_lo0797"),
                              c("PD49229b_lo0478","PD49229b_lo0529"),
                              c("PD49229b_lo0670","PD49229b_lo0723"),
                              c("PD49229b_lo0700","PD49229b_lo0924"),
                              c("PD53373b_lo0016","PD53373b_lo0187"),
                              c("PD53373b_lo0049","PD53373b_lo0140"),
                              c("PD53373b_lo0064","PD53373b_lo0219"),
                              c("PD53373b_lo0040","PD53373b_lo0193","PD53373b_lo0239"),
                              c("PD53373b_lo0049","PD53373b_lo0140"),
                              c("PD49228b_lo0368","PD49228b_lo0445"),
                              c("PD49228b_lo0507","PD49228b_lo0646"),
                              c("PD49228b_lo0084","PD49228b_lo0290"),
                              c("PD53373b_lo0176","PD53373b_lo0537"),
                              c("PD53373b_lo0042","PD53373b_lo0488")
                              
)
drop.samples=Map(tree=all.trees,details=all.muts,function(tree,details) {
  #Determine the sets of duplicate samples
  duplicate_samples=get_duplicate_status(tree,sample_metadata,res_type="list",height_cut_off=150)
  lapply(duplicate_samples,function(x) print(paste(paste(x$Sample,collapse=" & "),"are recognised as potential duplicate sets")))
  
  #If any samples in the set are in the manually curated sample set - use this to choose any samples to drop
  #Otherwise, choose one of each set of duplicate samples to keep - the one with the lower mutation burden (as suspect higher is usually in vitro acquired)
  drop_samples_list<-lapply(duplicate_samples,function(x) {
    samples<-x$Sample
    if(any(samples %in% unlist(true_invivo_related_sets))){
      return(samples[!samples%in%unlist(true_invivo_related_sets)])
    } else {
      if(sum(samples%in%tree$tip.label)>1) {
        included_samples=samples[samples%in%tree$tip.label]
        sample_heights=sapply(included_samples,function(sample) {nodeheight(tree = tree,node=which(tree$tip.label==sample))})
        retain_sample=included_samples[sample_heights==min(sample_heights)][1]
        return(included_samples[!included_samples==retain_sample])
      }else{
        return(NULL)
      }
    }
  })

  drop_samples=unlist(drop_samples_list)
  return(drop_samples)
})

#plot_duplicate_plate_map(tree = all.trees$BCL008,sample_metadata = sample_metadata,height_cut_off = 200)

##3. These trees/ details still contain some duplicate colonies, that need to be removed manually
remove_samples_manual=list(BCL002="PD49229b_lo0521",
                           BCL003=c(""),
                           BCL004=c("PD49228b_lo0057","PD49228b_lo0523"),
                           BCL006=c(""),
                           BCL008=c("PD49227b_lo0035","PD49227b_lo0300","PD49227b_lo0104","PD49227b_lo0109","PD49227b_lo0142","PD49227b_lo0146","PD49227b_lo0149","PD49227b_lo0106","PD49227b_lo0124","PD49227b_lo0191","PD49227b_lo0126","PD49227b_lo0076","PD49227b_lo0181"),
                           BCL009=c("PD49226b_lo0206","PD49226b_lo0214","PD49226b_lo0238"))

remove_samples_list=Map(auto=drop.samples,manual=remove_samples_manual,function(auto,manual) {c(auto,manual)})


all.muts.nodups<-Map(tree=all.trees,details=all.muts,remove_samples=remove_samples_list,f=function(tree,details,remove_samples) {
  if(length(remove_samples)>0){
    output<-remove_samples_from_tree_and_update_details(remove_samples=remove_samples,tree=tree,details=details)
    return(output$details)
  } else {
    return(details)
  }
})

##4. Add sample-level data to the metadata table: SNV burden, peak VAF
sample_metadata=dplyr::bind_rows(Map(genotype_mat=genotype.mats,details=all.muts,tree=all.trees,exp_ID=names(genotype.mats),function(genotype_mat,details,tree,exp_ID){
  cat(paste0(exp_ID,"\n"))
  
  #Load the full annotated muts object - need this to be able to do the 'peak vaf' analysis
  annotated_muts_path=grep(exp_ID,annotated_muts_paths,value=T); load(annotated_muts_path)
  sample_metadata_Ind<-sample_metadata%>%dplyr::filter(ID==exp_ID)
  
  #Calculate the 'peak VAF' measure for the sample
  sample_metadata_Ind$peak_vaf=sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree$tip.label){
      return(NA)
    } else {
      max.density=function(x){
        dens<-density(x)
        return(dens$x[which.max(dens$y)])
      }
      mut_vafs=get_mut_vafs(SampleID,COMB_mats=filtered_muts$COMB_mats.tree.build,tree=tree)
      return(max.density(mut_vafs$vaf))
    }
  })
  
  #Function to calculate the mutation burden from the mutation node assignments and the tree structure
  get_uncorrected_mut_burden=function(SampleID,tree,details,Mut_type=c("SNV","INDEL","all")){
    if(Mut_type=="all"){Mut_type<-c("SNV","INDEL")}
    if(!SampleID%in%tree$tip.label){
      return(NA)
    } else {
      sample_nodes=get_ancestral_nodes(which(tree$tip.label==SampleID),tree$edge)
      return(sum(details$node[details$Mut_type%in%Mut_type]%in%sample_nodes))
    }
  }
  
  #Calculate the uncorrected SNV mutation burden for each sample
  sample_metadata_Ind$SNV_burden_u=sapply(sample_metadata_Ind$Sample,get_uncorrected_mut_burden,tree=tree,details=details,Mut_type="SNV")
  
  #Calculate the uncorrected INDEL mutation burden for each sample
  sample_metadata_Ind$INDEL_burden_u=sapply(sample_metadata_Ind$Sample,get_uncorrected_mut_burden,tree=tree,details=details,Mut_type="INDEL")
  
  #Now get the raw mutation burdens from the genotype matrix
  genotype_mat.SNV<-genotype_mat[details$mut_ref[details$Mut_type=="SNV"],]
  genotype_mat.Indel<-genotype_mat[details$mut_ref[details$Mut_type=="INDEL"],]
  raw_SNV_burdens=colSums(genotype_mat.SNV==1)
  raw_INDEL_burdens=colSums(genotype_mat.Indel==1)
  
  #Combine data
  sample_metadata_Ind<-left_join(sample_metadata_Ind,data.frame(Sample=names(raw_SNV_burdens),SNV_burden_raw=raw_SNV_burdens,INDEL_burden_raw=raw_INDEL_burdens),by="Sample")
  
  return(sample_metadata_Ind)
}))

##5. Import the library prep QC results & join to the metadata
conc_files<-list.files("~/Google Drive File Stream/My Drive/Gene_therapy_project/QC results/Library prep QC/",pattern="^6260.*MSC[_]*[Dd]",full.names = T)
Library_QC_results<-dplyr::bind_rows(lapply(conc_files,function(file){
  dat<-readxl::read_excel(file,skip=8)
  dat<-dplyr::select(dat,`SUPPLIER SAMPLE NAME`,`CONC. (ng/ul)`)%>%dplyr::rename("Sample"=`SUPPLIER SAMPLE NAME`,"Library_conc"=`CONC. (ng/ul)`)
  return(dat)
}))
sample_metadata<-Library_QC_results%>%
  dplyr::filter(!duplicated(Sample))%>%
  right_join(sample_metadata)
write.table(sample_metadata,file=paste0(root_dir,"/Data/metadata_temp.tsv"),quote = F,sep = "\t",row.names = F)

#The metadata_temp.tsv file is then fed into the 'Clonality_correct.R' script that infers the 'somatic SNV sensitivity'
# based on the peak_vaf measure (i.e. the clonality of the sample) and the germline SNV sensitivity

#Import the updated version once infer germline sensitivity/ somatic SNV sensitivity
sample_metadata<-read.delim(paste0(root_dir,"/Data/metadata_temp_updated.tsv"))
sens_df<-sample_metadata%>%
  dplyr::select(Sample,somatic_SNV_sensitivity)%>%
  dplyr::rename("SNV_sensitivity"=somatic_SNV_sensitivity)%>%
  rbind(data.frame(Sample="Ancestral",SNV_sensitivity=1))%>%
  filter(!duplicated(Sample))

sample_metadata<-Library_QC_results%>%
  dplyr::filter(!duplicated(Sample))%>%
  right_join(sample_metadata)

##6. Perform an asymptotic regression step on the INDEL and SNV mutation burdens
## This makes use of the higher coverage samples in the sample set to correct the lower coverage samples to their level
## The asymptotic regression assumes as asymptotic relationship between coverage and number of called mutations

test_groups=expand.grid(c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009"),c("Pre","Post"))

AR_results<-lapply(1:nrow(test_groups),function(i) {
  exp_ID=test_groups$Var1[i]
  print(exp_ID)
  if(test_groups$Var2[i]=="Pre") {
    test_samples=sample_metadata%>%
      dplyr::filter(ID==exp_ID&Time_point==0)%>%
      dplyr::select(Sample,ID,Time_point,Coverage,peak_vaf,SNV_burden_u,INDEL_burden_u,SNV_burden_raw,INDEL_burden_raw)
    model_samples=test_samples%>%filter(!is.na(SNV_burden_raw)&
                                          peak_vaf>0.46&
                                          INDEL_burden_raw!=0)
  } else if(test_groups$Var2[i]=="Post"){
    test_samples=sample_metadata%>%
      dplyr::filter(ID==exp_ID&Time_point>0)%>%
      dplyr::select(Sample,ID,Time_point,Coverage,peak_vaf,SNV_burden_u,INDEL_burden_u,SNV_burden_raw,INDEL_burden_raw)
    if(nrow(test_samples)==0){stop(return(NULL))}
    high_coverage_timepoint=sample_metadata%>%
      dplyr::filter(ID==exp_ID&Time_point)%>%
      group_by(Time_point,"high_cov"=Coverage>25)%>%
      dplyr::summarise(n=n())%>%
      filter(high_cov&n>5)%>%
      pull(Time_point)
    model_samples=test_samples%>%filter(Time_point==high_coverage_timepoint&
                                          !is.na(SNV_burden_raw)&
                                          peak_vaf>0.46&
                                          INDEL_burden_raw!=0)
  }
  
  #Define SNV regression
  #Model is defined with samples from a single time point - one that includes higher coverage samples
  idxs<-1:nrow(model_samples)
  model.as.SNV<-NULL
  while(length(model.as.SNV)!=3){
    #Sometimes model doesn't converge - inwhich case subsample the data and reattempt
    model.as.SNV = try(NLSstAsymptotic(sortedXyData(model_samples$Coverage[idxs],model_samples$SNV_burden_u[idxs])),silent=T)
    idxs<-sample(1:nrow(model_samples),size = ceiling(nrow(model_samples)*0.8))
  }
  mymodel.SNV = function(x){b0 = model.as.SNV[1]; b1 = model.as.SNV[2]; lrc = model.as.SNV[3]; b0 + b1*(1-exp(-exp(lrc) * x))}
  
  #Define INDEL regression
  #Model is defined with samples from a single time point - one that includes higher coverage samples
  idxs<-1:nrow(model_samples)
  model.as.INDEL<-NULL
  while(length(model.as.INDEL)!=3){
    #Sometimes model doesn't converge - inwhich case subsample the data and reattempt
    model.as.INDEL = try(NLSstAsymptotic(sortedXyData(model_samples$Coverage[idxs],model_samples$INDEL_burden_u[idxs])),silent=T)
    idxs<-sample(1:nrow(model_samples),size = ceiling(nrow(model_samples)*0.8))
  }
  mymodel.INDEL = function(x){b0 = model.as.INDEL[1]; b1 = model.as.INDEL[2]; lrc = model.as.INDEL[3]; b0 + b1*(1-exp(-exp(lrc) * x))}
  
  # model.as.SNV <- drm(sortedXyData(model_samples$Coverage,model_samples$SNV_burden_u), fct = DRC.asymReg())
  # mymodel.SNV = function(x){b0 = model.as.SNV$coefficients[1]; b1 = model.as.SNV$coefficients[2]; lrc = model.as.SNV$coefficients[3]; b0 + b1*(1-exp(-exp(lrc) * x))}
  # 
  # model.as.INDEL = drm(sortedXyData(model_samples$Coverage,model_samples$INDEL_burden_u),fct=DRC.asymReg(fixed = c(NA, -2.5, NA), names = c("init", "m", "plateau")))
  # mymodel.INDEL = function(x){b0 = model.as.INDEL$coefficients[3]; b1 = model.as.INDEL$coefficients[2]; lrc = model.as.INDEL$coefficients[1]; b0 + b1*(1-exp(-exp(lrc) * x))}
  
  
  # Plot SNV model
  myrange=1:round(max(model_samples$Coverage))
  asp.SNV = data.frame(x=myrange, y=unlist(lapply(myrange, FUN=mymodel.SNV)))
  pre_correction_plot_SNVs<-ggplot(model_samples)+ylim(1,1000)+xlim(0,50)+
    geom_point(aes(x=Coverage,y=SNV_burden_raw,col=factor(Time_point)), data=model_samples, alpha=0.25) + 
    theme_bw()+scale_fill_brewer(type="qual", palette=3)+geom_line(aes(x=x, y=y), data=asp.SNV)+
    labs(y="# Mutations",x="Coverage per colony",col="Time point")+my_theme
  
  # Plot the INDEL model
  asp.INDEL = data.frame(x=myrange, y=unlist(lapply(myrange, FUN=mymodel.INDEL)))
  pre_correction_plot_INDELs<-ggplot(model_samples)+ylim(1,40)+xlim(0,50)+
    geom_point(aes(x=Coverage,y=INDEL_burden_raw,col=factor(Time_point)), data=model_samples, alpha=0.25) + 
    theme_bw()+scale_fill_brewer(type="qual", palette=3)+geom_line(aes(x=x, y=y),data=asp.INDEL)+
    labs(y="# Mutations",x="Coverage per colony",col="Time point")+my_theme
  
  # Normalize mutation burden to a read depth (here I use 30x). Use a depth within your range of coverage.
  new_Number_mutations.SNV = function(rd, Number_mutations){
    Number_mutations + ( mymodel.SNV(30) - mymodel.SNV(rd)) ## can change depth here
  }
  new_Number_mutations.INDEL = function(rd, Number_mutations){
    Number_mutations + ( mymodel.INDEL(30) - mymodel.INDEL(rd)) ## can change depth here
  }
  
  #The model is then applied slightly more broadly - to samples from several post-GT time points
  model_samples$SNV_burden_AR = apply(model_samples, MARGIN=1, FUN=function(X)  new_Number_mutations.SNV(rd=as.numeric(X[which(colnames(model_samples)=="Coverage")]), Number_mutations=as.numeric(X[which(colnames(model_samples)=="SNV_burden_u")])))
  model_samples$INDEL_burden_AR = apply(model_samples, MARGIN=1, FUN=function(X)  new_Number_mutations.INDEL(rd=as.numeric(X[which(colnames(model_samples)=="Coverage")]), Number_mutations=as.numeric(X[which(colnames(model_samples)=="INDEL_burden_raw")])))
  
  post_correction_plot_SNVs<-ggplot(model_samples,aes(Coverage,SNV_burden_AR,col=factor(Time_point)))+ylim(0,1000)+xlim(0,50)+theme_bw()+
    geom_point(alpha=0.25)+geom_smooth(method="lm",size=0.5)+labs(y="# Mutations",x="Coverage per colony",col="Time point")+my_theme
  
  post_correction_plot_INDELs<-ggplot(model_samples,aes(Coverage,INDEL_burden_AR,col=factor(Time_point)))+ylim(0,50)+xlim(0,50)+theme_bw() +
    geom_point(alpha=0.25)+geom_smooth(method="lm",size=0.5)+labs(y="# Mutations",x="Coverage per colony",col="Time point")+my_theme
  
  comb_plot=gridExtra::arrangeGrob(pre_correction_plot_SNVs,pre_correction_plot_INDELs,post_correction_plot_SNVs,post_correction_plot_INDELs)
  plot(comb_plot)
  
  return(model_samples)
})%>%dplyr::bind_rows()

sample_metadata<-left_join(sample_metadata,AR_results)

##7. Import the vector coverage info
vector_coverage<-read.delim(paste0(root_dir,"/Data/vector_coverage_all.txt"),header=F,stringsAsFactors = F)[,c(1,8)]
VCN_by_VIS<-read.delim(paste0(root_dir,"/Data/VCN_by_VIS.txt"),header=T,stringsAsFactors = F)
colnames(vector_coverage)<-c("Sample","Vector_coverage")
sample_metadata<-sample_metadata%>%
  left_join(vector_coverage)%>%
  left_join(VCN_by_VIS)%>%
  mutate(VCN=Vector_coverage/(0.5*Coverage))%>%
  mutate(VCN_rounded=round(VCN-0.18))

##8a. Get (1) uncorrected set of trees,
#(2) trees corrected for sensitivity & clonality,
#(3) trees corrected for sensitivity & clonality - without duplicates,
#(4) Ultrametric trees

#Uncorrected trees
all.trees.uncorrected<-Map(tree=all.trees,details=all.muts,f=function(tree,details){
  tree$edge.length<-sapply(tree$edge[,2],function(node) return(sum(details$node==node & details$Mut_type=="SNV")))
  return(tree)
})

#Trees corrected for sensitivity & clonality
all.trees.cc<-Map(tree=all.trees.uncorrected,details=all.muts,remove_samples=remove_samples_list,f=function(tree,details,remove_samples) {
  tree_SNV_c<-get_corrected_tree(tree,details,sensitivity_df = sens_df,include_SNVs = T,include_indels = F,get_edge_from_tree = T)
})

#Trees corrected for sensitivity & clonality - without duplicates
all.trees.cc.nodups<-Map(tree=all.trees.uncorrected,details=all.muts,remove_samples=remove_samples_list,f=function(tree,details,remove_samples) {
  tree_SNV_c<-get_corrected_tree(tree,details,sensitivity_df = sens_df,include_SNVs = T,include_indels = F,get_edge_from_tree = T)
  tree_SNV_c<-drop.tip(tree_SNV_c,remove_samples)
})

#Trees made ultrametric from the sensitivity/ clonality-corrected tree with no duplicates
all.trees.ultra<-lapply(all.trees.cc.nodups,function(tree) {
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0
  tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(tree))
  return(tree.ultra)
})

##.Now get the SNV burdens from the tree after this correction for sensitivity based on clonality & coverage
sample_metadata=dplyr::bind_rows(Map(details=all.muts,tree=all.trees.cc,exp_ID=names(all.trees),function(details,tree,exp_ID){
  cat(paste0(exp_ID,"\n"))
  sample_metadata_Ind<-sample_metadata%>%dplyr::filter(ID==exp_ID)
  sample_metadata_Ind$SNV_burden_tc<-sapply(sample_metadata_Ind$Sample,function(SampleID) {
    if(!SampleID%in%tree$tip.label) {return(NA)} else {return(nodeheight(tree,which(tree$tip.label==SampleID)))}
  })
  return(sample_metadata_Ind)
}))

##10. Summarise the colony outcomes
sample_metadata=sample_metadata%>%
  mutate(sample_status=ifelse(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)),"PASS",
                              ifelse(Sample%in%unlist(remove_samples_list),"Duplicate",
                                     ifelse(Coverage<4,"Low coverage",
                                            ifelse(is.na(Coverage),"Not Sequenced","Non-clonal")))))
sample_metadata$sample_status=factor(sample_metadata$sample_status,levels=rev(c("PASS","Low coverage","Not Sequenced","Non-clonal","Duplicate")))

##11. Add in the individual metadata into this df
sample_metadata<-left_join(sample_metadata,Individual_metadata)

##12. Save this 'complete' metadata dataframe & the tree and mutations objects
write_delim(sample_metadata,file = paste0(root_dir,"/Data/sample_metadata_full.tsv"),delim="\t")
saveRDS(list(all.trees.uncorrected=all.trees.uncorrected,
             all.trees.cc=all.trees.cc,
             all.trees.cc.nodups=all.trees.cc.nodups,
             all.trees.ultra=all.trees.ultra),
        file=paste0(root_dir,"/Data/combined_tree_files.Rds"))
saveRDS(list(all.muts=all.muts,
             all.muts.nodups=all.muts.nodups),
        file=paste0(root_dir,"/Data/combined_muts_files.Rds"))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


##12. Check for clustering of different sample types (generally post-GT relative to pre-GT) on the tree
## Would imply differential 'transplantability' originating in embryogenesis
pdf(paste0(plots_dir,"AMOVA_cell_types_postGT_related_removed.pdf"),width=7,height=7)
par(mfrow=c(3,2))
for(tree in all.trees.ultra) {c
  
  #Don't want to include multiple samples from the postGT clones - these just reflect the transplant bottleneck, not inherited "transplantability" of embryonic clones
  remove_postGT_clone_samples<-unlist(lapply(true_invivo_related_sets,function(dup_set) {
    if(sum(dup_set%in%tree$tip.label)>1) {return(dup_set[-1])} else {return(NULL)}
  }))
  if(length(remove_postGT_clone_samples)>0) {
    cat(paste("Removing",remove_postGT_clone_samples),sep="\n")
    tree<-drop.tip(tree,remove_postGT_clone_samples)
  }
  
  
  cellkey<-sample_metadata%>%
    dplyr::filter(Sample%in%tree$tip.label)%>%
    dplyr::select(Sample,Sample_type)%>%
    dplyr::rename("Cell_type"=Sample_type)
  
  
  if(length(unique(cellkey$Cell_type))>1) {
    ID<-sample_metadata%>%
      dplyr::filter(Sample%in%tree$tip.label)%>%
      dplyr::slice(1)%>%pull(ID)
    cat(ID,sep="\n")
    
    distmat <- cophenetic(drop.tip(tree,"Ancestral"))
    
    amovapval.fn(distmat = distmat,
                 groupnames = unique(cellkey$Cell_type),
                 cell_key = cellkey,
                 iterations = 30000,
                 plottitle = paste0(ID,": Phylogenetic clustering of sample types"))
  }
}
dev.off()

##13. Review relationship between VCN and mutation burden in post-gene therapy samples
sample_metadata%>%
  filter(Sample_type=="Post-GT")%>%
  ggplot(aes(x=VCN_rounded,y=SNV_burden_AR,col=new_ID))+
  geom_jitter(width=0.2,alpha=0.5,size=0.5)+
  facet_grid(~new_ID)+
  geom_smooth(method="lm")+
  scale_color_manual(values=patient_cols)+
  theme_bw()+
  labs(x="Vector copy number",y="SNV burden")

#Now comparison of just transduced vs untransduced
transduced_vs_untransduced_df<-sample_metadata%>%
  #left_join(Individual_metadata,by="ID")%>%
  filter(Sample_type=="Post-GT" & sample_status=="PASS")%>%
  mutate(Transduced=ifelse(VCN_rounded>0.8,"Gene\nmod.","Non\nmod."),Time_point=paste(Time_point,"years"))

transduced_vs_untransduced_comparison<-lapply(c("SCD1","SCD2","SCD3","SCD4","SCD5","SCD6"),function(exp_ID) {
  time_points=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID)%>%pull(Time_point)%>%unique()
  lapply(time_points,function(time_point){
    test=t.test(x=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Gene\nmod.")%>%pull(SNV_burden_tc),
                  y=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Non\nmod.")%>%pull(SNV_burden_tc))
    return(data.frame(new_ID=exp_ID,Time_point=time_point,estimate=test$estimate[1]-test$estimate[2],p.value=test$p.value,CI_low=test$conf.int[1],CI_high=test$conf.int[2]))
  })%>%bind_rows()
})%>%bind_rows

transduced_vs_nontransduced_plot<-sample_metadata%>%
  #left_join(Individual_metadata,by="ID")%>%
  filter(Sample_type=="Post-GT" & sample_status=="PASS")%>%
  mutate(Transduced=ifelse(VCN_rounded>0.8,"Gene\nmod.","Non\nmod."),Time_point=paste(Time_point,"years"))%>%
  ggplot(aes(x=Transduced,y=SNV_burden_tc))+
  geom_boxplot(col="black",size=0.2,outlier.shape=NA)+
  geom_jitter(aes(col=new_ID),width=0.15,alpha=0.3,size=0.2)+
  scale_color_manual(values=patient_cols)+
  facet_wrap(new_ID~Time_point,nrow=2)+
  geom_text(data=transduced_vs_untransduced_comparison,aes(label=paste("p =",round(p.value,digits=2)),x=1.5,y=800),size=2,fontface="italic",inherit.aes = F)+
  geom_smooth(method="lm")+
  theme_bw()+
  labs(x="",y="SNV burden",col="ID")+
  my_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme(strip.text.x=element_text(size=4,margin = unit(c(1,0,1,0),"mm")),legend.key.size = unit(2.5,"mm"),legend.box.spacing = unit(0,"mm"))

ggsave(filename = paste0(plots_dir,"transduced_vs_nontransduced_plot.pdf"),plot=transduced_vs_nontransduced_plot,width=4.5,height=2)

##15. Review the coverage histogram across samples
coverage_summary<-sample_metadata%>%
  #left_join(Individual_metadata,by="ID")%>%
  filter(!is.na(Coverage) &Coverage >4&Disease=="SCD")%>%
  group_by(new_ID)%>%
  summarise(mean_coverage=mean(Coverage),median_coverage=median(Coverage))

coverage_histograms<-sample_metadata%>%
  filter(!is.na(Coverage) &Coverage >4&Disease=="SCD")%>%
  ggplot(aes(x=Coverage))+
  geom_histogram(fill="lightblue",col="black",bins=30)+
  scale_x_continuous(limits=c(0,40))+
  facet_wrap(~new_ID,scales="free")+
  geom_vline(data=coverage_summary,aes(xintercept=mean_coverage),col="red",linetype=2)+
  geom_text(data=coverage_summary,aes(label=paste0("Mean coverage = ",signif(mean_coverage,digits=3)," x"),x=25,y=80),size=2,inherit.aes = F)+
  theme_classic()+
  my_theme+
  labs(y="Count")

ggsave(filename=paste0(plots_dir,"Coverage_histograms.pdf"),coverage_histograms,width=3.5,height=2.5)

##16. Do plot of mutation burden for age, comparing with Emily & Nick's regression lines
HSC_prediction=function(age){
  pred=54.6+16.8*age; upper_95_CI=54.6+16.95*age; lower_95_CI=54.6+16.65*age
  return(list(mean=pred,upper_95_CI=upper_95_CI,lower_95_CI=lower_95_CI))
}

MPN_prediction=function(age){
  pred=98.3+16.5*age; upper_95_CI=98.3+18.2*age; lower_95_CI=98.3+14.8*age
  return(list(mean=pred,upper_95_CI=upper_95_CI,lower_95_CI=lower_95_CI))
}

x=c(0:26)
lm.data=data.frame(x=x,ymin=unlist(sapply(x,HSC_prediction)['lower_95_CI',]),ymax=unlist(sapply(x,HSC_prediction)['upper_95_CI',]))
lm.data.MPN=data.frame(x=x,ymin=unlist(sapply(x,MPN_prediction)['lower_95_CI',]),ymax=unlist(sapply(x,MPN_prediction)['upper_95_CI',]))


SNV_stat_to_use<-"SNV_burden_tc"
mut_burden_summary<-left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(get(SNV_stat_to_use)))%>%
  dplyr::filter(Coverage>10)%>%
  mutate(Age_at_sampling=Time_point+Age_at_GT)%>%
  group_by(new_ID,Time_point)%>%
  summarise(mean_SNV_burden=mean(get(SNV_stat_to_use)),var=var(get(SNV_stat_to_use)),Age_at_sampling=Age_at_sampling[1])

all_samples_burden_plot<-left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(get(SNV_stat_to_use)))%>%
  dplyr::filter(Coverage>10)%>%
  mutate(Age_at_sampling=Time_point+Age_at_GT)%>%
  ggplot(aes(x=Age_at_sampling,y=get(SNV_stat_to_use),col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.1,size=0.5)+
  scale_y_continuous(limits=function(x) c(0,750))+
  scale_x_continuous(limits=function(x) c(0,x[2]))+
  geom_abline(intercept = 54.6,slope = 16.8,size=0.3)+
  geom_ribbon(data=lm.data,aes(x=x,ymin=ymin, ymax=ymax),fill="darkgrey",alpha=0.7,inherit.aes = F)+
  geom_point(data=mut_burden_summary,aes(x=Age_at_sampling,y=mean_SNV_burden),shape=3,size=1,inherit.aes = F)+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_bw()+
  labs(x="Age at sampling",
       y="Single nucleotide variant burden",
       col="Individual")+
  my_theme

preGT_mut_burden_plot<-sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_stat_to_use)))%>%
  dplyr::filter(Coverage>10)%>%
  dplyr::filter(Time_point==0)%>%
  mutate(Age_at_sampling=Time_point+Age_at_GT)%>%
  ggplot(aes(x=Age_at_sampling,y=get(SNV_stat_to_use),col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.1,size=0.3)+
  scale_y_continuous(limits=function(x) c(0,750))+
  scale_x_continuous(limits=function(x) c(0,x[2]))+
  geom_abline(intercept = 54.6,slope = 16.8,size=0.3)+
  geom_ribbon(data=lm.data,aes(x=x,ymin=ymin, ymax=ymax),fill="darkgrey",alpha=0.7,inherit.aes = F)+
  geom_point(data=mut_burden_summary%>%left_join(Individual_metadata)%>%filter(Time_point==0)%>%mutate(Age_at_sampling=Time_point+Age_at_GT),aes(x=Age_at_sampling,y=mean_SNV_burden),shape=3,size=1,inherit.aes = F)+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_bw()+
  labs(x="Age at sampling",
       y="Single nucleotide variant burden",
       col="Individual")+
  my_theme

ggsave(filename = paste0(plots_dir,"preGT_mut_burden.pdf"),preGT_mut_burden_plot,width=3,height=2)

##17. Do comparisons of pre & post GT
#Perform mutation burden comparisons only using the higher coverage samples with best evidence of clonality
pre_post_comparison_df<-left_join(sample_metadata,Individual_metadata,by="ID")%>%
  dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL002","BCL004") & sample_status=="PASS"&Coverage>10 & peak_vaf>0.45 & peak_vaf<0.55)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  mutate(SNV_adjusted=SNV_burden-(16.8*Time_point/12))

pre_post_comparison_df%>%
  group_by(ID,Time_point)%>%
  summarize(mean(SNV_burden))

left_join(sample_metadata,Individual_metadata,by="ID")%>%
  dplyr::filter(!is.na(SNV_burden) & Disease=="SCD"&sample_status=="PASS"&Coverage>10 & peak_vaf>0.45 & peak_vaf<0.55)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  mutate(SNV_adjusted=SNV_burden-(16.8*Time_point/12))%>%
  ggplot(aes(x=factor(Time_point),y=SNV_burden,col=Cell_type))+
  geom_boxplot(aes(group=plate_ID))+
  geom_jitter(width = 0.1,alpha=0.2)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  #scale_color_brewer(palette="Dark2")+
  labs(x="Time point (months)", y="SNV burden")+
  theme(legend.text = element_text(size=7),legend.key.size = unit(3,"mm"))

pre_vs_post_comparison_plot<-pre_post_comparison_df%>%
  mutate(Time_point=Time_point/12)%>%
  ggplot(aes(x=Time_point,y=SNV_burden))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size = 0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.025,alpha=0.3,size=0.7,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_color_brewer(palette="Dark2")+
  scale_x_continuous(breaks=seq(0,3,1))+
  labs(x="Time point (years)", y="SNV burden")+
  my_theme

lm(SNV_burden~Time_point,data=pre_post_comparison_df%>%filter(new_ID=="SCD3")%>%mutate(Time_point=Time_point/12))
lm(SNV_burden~Time_point,data=pre_post_comparison_df%>%filter(new_ID=="SCD4")%>%mutate(Time_point=Time_point/12))

pre_vs_post_comparison_adjusted_plot<-pre_post_comparison_df%>%
  mutate(Time_point=Time_point/12)%>%
  ggplot(aes(x=Time_point,y=SNV_adjusted))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size=0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.025,alpha=0.3,size=0.7,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(250,750))+
  scale_x_continuous(breaks=seq(0,3,1))+
  labs(x="Time point (years)", y="Age-adjusted SNV burden")+
  my_theme

lm(SNV_adjusted~Time_point,data=pre_post_comparison_df%>%filter(new_ID=="SCD3")%>%mutate(Time_point=Time_point/12))
lm(SNV_adjusted~Time_point,data=pre_post_comparison_df%>%filter(new_ID=="SCD4")%>%mutate(Time_point=Time_point/12))

ggsave(paste0(plots_dir,"pre_vs_post_comparison_plot.pdf"),plot=pre_vs_post_comparison_plot,width=2.5,height=2)
ggsave(paste0(plots_dir,"pre_vs_post_comparison_adjusted_plot.pdf"),plot=pre_vs_post_comparison_adjusted_plot,width=2.5,height=2)

pre_post_comparison_df%>%
  ggplot(aes(x=factor(Time_point),y=SNV_adjusted,col=Cell_type))+
  geom_boxplot()+
  geom_jitter(width = 0.3,alpha=0.3,size=1,stroke=0)+
  facet_grid(cols=vars(ID),scales="free",space = "free")+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(250,750))+
  labs(x="Time point (months)", y="Age-adjusted SNV burden")


BCL002.test<-t.test(y=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & ID=="BCL002" & Coverage>10 &peak_vaf>0.45 &Time_point==0)%>%pull(SNV_adjusted),
       x=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & ID=="BCL002" & Coverage>10 &peak_vaf>0.45 &Time_point>0)%>%pull(SNV_adjusted),
       alternative="two.sided")

BCL004.test<-t.test(y=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & ID=="BCL004" & Coverage>10 &peak_vaf>0.45 &Time_point==0)%>%pull(SNV_adjusted),
       x=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & ID=="BCL004" & Coverage>10 &peak_vaf>0.45 &Time_point>0)%>%pull(SNV_adjusted),
       alternative="two.sided")

gt_induced_mutations_estimate_plot<-Map(test=list(BCL002.test,BCL004.test),exp_ID=c("SCD3","SCD4"),function(test,exp_ID){
  data.frame(ID=exp_ID,estimate=test$estimate[1]-test$estimate[2],CI_low=test$conf.int[1],CI_high=test$conf.int[2])
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=ID,col=ID,y=estimate,ymin=CI_low,ymax=CI_high))+
  geom_point(size=1)+
  scale_color_brewer(palette="Set1")+
  geom_errorbar(width=0.2,alpha=0.5)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-50,50))+
  theme_bw()+
  labs(x="",y="Estimate of excess mutations\nfrom gene therapy")+
  my_theme+
  theme(legend.box.spacing=unit(1,"mm"))

#Do t tests for each individual timepoint
gt_induced_mutations_by_timepoint_plot<-lapply(c("SCD3","SCD4"),function(exp_ID) {
  post_time_points<-pre_post_comparison_df%>%dplyr::filter(new_ID==exp_ID & Time_point!=0)%>%pull(Time_point)%>%unique()
  lapply(post_time_points,function(post_time_point) {
    test<-t.test(y=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & new_ID==exp_ID & Coverage>10 &peak_vaf>0.45 &Time_point==0)%>%pull(SNV_adjusted),
                 x=pre_post_comparison_df%>%filter(!is.na(SNV_burden) & new_ID==exp_ID & Coverage>10 &peak_vaf>0.45 &Time_point==post_time_point)%>%pull(SNV_adjusted),
                 alternative="two.sided")
    return(data.frame(ID=exp_ID,Time_point=post_time_point,estimate=test$estimate[1]-test$estimate[2],CI_low=test$conf.int[1],CI_high=test$conf.int[2]))
  })
})%>%bind_rows()%>%
  ggplot(aes(x=Time_point,col=ID,y=estimate,ymin=CI_low,ymax=CI_high))+
  geom_point(size=1)+
  scale_color_brewer(palette="Set1")+
  geom_errorbar(width=1,alpha=0.5)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-50,50))+
  scale_x_continuous(limits=c(0,37),breaks=seq(0,36,12))+
  facet_grid(~ID)+
  theme_bw()+
  labs(x="Time point (months)",y="Estimate of excess mutations\nfrom gene therapy")+
  my_theme

ggsave(paste0(plots_dir,"gt_induced_mutations.pdf"),plot=gt_induced_mutations_estimate_plot,width=1.8,height=2)
ggsave(paste0(plots_dir,"gt_induced_mutations_by_timepoint.pdf"),plot=gt_induced_mutations_by_timepoint_plot,width=2.5,height=2)


#Histograms of mutation burdens pre & post
left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL002","BCL004") & Coverage>10 &peak_vaf>0.45)%>%
  ggplot(aes(x=SNV_burden))+
  geom_histogram()+
  facet_grid(Time_point~ID)

#Print table of average burdens pre & post
left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL002","BCL004") & Coverage>10 &peak_vaf>0.45)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  mutate(SNV_burden_c=SNV_burden_u/somatic_SNV_sensitivity)%>%
  group_by(ID,Time_point)%>%
  summarise(median=median(SNV_burden),mean=mean(SNV_burden),var=var(SNV_burden))

##18. Print table of average burdens in pre-GT/pre-TDX & DP (where available) and do t.test comparison
left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL009") & Coverage>10)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  group_by(ID,Sample_type)%>%
  summarise(median=median(SNV_burden),mean=mean(SNV_burden),var=var(SNV_burden))

#t.test comparison of pre-TDX and donor product sample in BCL009
t.test(x=left_join(sample_metadata,Individual_metadata)%>%
         dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL009") & Coverage>10 & Sample_type=="Pre-TDX")%>%
         pull(SNV_burden),
       y=left_join(sample_metadata,Individual_metadata)%>%
         dplyr::filter(!is.na(SNV_burden) & ID%in%c("BCL009") & Coverage>10 & Sample_type=="DP")%>%
         pull(SNV_burden))

##19. Compare SNV burden distributions by Individual/ time-point
left_join(sample_metadata,Individual_metadata)%>%
  dplyr::filter(!is.na(SNV_burden) & Coverage>10)%>%
  ggplot(aes(x=SNV_burden))+
  geom_histogram(bins=20)+
  facet_grid(cols=vars(ID),rows=vars(Time_point),scales="free",space = "free")+
  theme_bw()

##20. Compare raw VCN distribution histograms for all sample types/ individuals
VCN_histos_by_ID_and_timepoint<-sample_metadata%>%
  filter(!is.na(VCN))%>%
  ggplot(aes(x=VCN-0.3))+
  geom_histogram(binwidth = 0.2)+
  #geom_density()+
  scale_x_continuous(breaks = seq(0,12,1))+
  facet_grid(Sample_type~new_ID,scales="free")+
  theme_bw()+
  my_theme+
  labs(x="Vector copy number")
ggsave(filename=paste0(plots_dir,"VCN_histos_by_ID_and_timepoint.pdf"),VCN_histos_by_ID_and_timepoint,width=6,height=4)

###--------------Plot copy number change results---------------

#Plot of copy number changes by individual
CN_changes_per_individual_plot<-CN_change_df%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(sample_metadata,by="Sample")%>%
  left_join(Individual_metadata,by="ID")%>%
  filter(sample_status=="PASS")%>%
  ggplot(aes(x=factor(Sample_type,levels=c("Pre-GT","Pre-TDX","Post-GT")),y=1,fill=factor(Copy_number_change)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15,width=0.5)+
  facet_grid(cols=vars(factor(new_ID)),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="Sample type",y="Samples",fill="Copy number abnormality")

ggsave(filename=paste0(plots_dir,"CN_changes_by_individual.pdf"),CN_changes_per_individual_plot,width=4,height=2)

#Proportion of colonies with a CNA in different sample groups - i.e. pre-GT vs post-GT (excluding LOY)
right_join(CN_change_df,sample_metadata)%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups[1:4],function(tree) tree$tip.label)))%>%
  replace_na(replace=list(Copy_number_change="Normal"))%>%
  mutate(CNA=ifelse(Copy_number_change=="Normal","Normal","Abnormal"))%>%
  group_by(CNA,Time_point>0)%>%
  summarise(n=n())%>%
  mutate(Time_point=ifelse(`Time_point > 0`,"Post-GT","Pre-GT"))%>%
  pivot_wider(id_cols=c("CNA","Time_point"),names_from="CNA",values_from="n")%>%
  replace_na(replace=list(Abnormal=0,Normal=0))%>%
  mutate(abnormal_prop=Abnormal/(Abnormal+Normal),
         CI_low=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=factor(Time_point,levels=c("Pre-GT","Post-GT")),col=factor(Time_point,levels=c("Pre-GT","Post-GT")),y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  theme_bw()+
  scale_y_continuous(limits=c(0,4))+
  theme(legend.position = "none")+
  labs(x="Donor or Recipient",y="Proportion of colonies with CNA (%)")

#Review the rounded VCN by time point for samples that contain transduced cells
VCN_histograms<-sample_metadata%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  mutate(Sample_type=factor(Sample_type,levels=c("DP","Post-GT")))%>%
  mutate(Time_point=paste("Time point =",Time_point))%>%
  ggplot(aes(x=VCN_rounded,fill=Sample_type))+
  geom_bar(stat="count",width=0.7,col=NA)+
  scale_x_continuous(breaks=seq(0,12,1))+
  facet_wrap(new_ID~Time_point,scales="free_y")+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=6),strip.text = element_text(size=6))+
  labs(x="Vector copy number",y="Count",fill="Sample type")

ggsave(filename=paste0(plots_dir,"VCN_histos_plot.pdf"),VCN_histograms,width=5,height=3.5)

vector_containing_cell_prop_plot<-sample_metadata%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  mutate(Sample_type=factor(Sample_type,levels=c("DP","Post-GT")))%>%
  group_by(ID,Time_point)%>%
  summarise(n_untransduced=sum(VCN_rounded==0),n_VCN1=sum(VCN_rounded==1),n_VCN_over1=sum(VCN_rounded>1))%>%
  mutate(untransduced=n_untransduced/(n_untransduced+n_VCN1+n_VCN_over1),
         `VCN=1`=n_VCN1/(n_untransduced+n_VCN1+n_VCN_over1),
         `VCN >1`=n_VCN_over1/(n_untransduced+n_VCN1+n_VCN_over1))%>%
  dplyr::select(-n_untransduced,-n_VCN1,-n_VCN_over1)%>%
  mutate(Transduced=1-untransduced)%>%
  gather(-ID,-Time_point,key="VCN",value="Proportion")%>%
  filter(VCN=="Transduced")%>%
  left_join(Individual_metadata,by="ID")%>%
  ggplot(aes(x=Time_point,y=Proportion,col=new_ID))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_x_continuous(limits=c(0,36),breaks=seq(0,36,6))+
  theme_classic()+
  labs(x="Time point (months)",y="Proportion of colonies\ncontaining vector",col="Individual")+
  my_theme

ggsave(filename=paste0(plots_dir,"vector_containing_cell_prop_plot.pdf"),vector_containing_cell_prop_plot,width=2.5,height=2)

##21. Compare raw vector copy number (by coverage) and the VCN by number of vector integration sites identified.
library(ggrepel)
VCN_vs_VCN_by_VIS<-sample_metadata%>%
  left_join(VCN_by_VIS)%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  ggplot(aes(x=VCN-0.3,y=VCN_by_VIS,label=Sample,col=Sample_type))+
  geom_abline(intercept = 0,slope=1,linetype=2)+
  #geom_label_repel()+
  geom_jitter(size=0.5,height=0.15,width=0,alpha=0.3) +
  scale_y_continuous(breaks = seq(0,10,1),limits=c(0,10))+
  scale_x_continuous(breaks = seq(0,10,1),limits=c(0,10))+
  facet_grid(~new_ID)+
  theme_bw()+
  my_theme+
  labs(x="Vector copy number",y="Number of vector integration\nsites identified")

ggsave(filename=paste0(plots_dir,"VCN_vs_VCN_by_VIS.pdf"),VCN_vs_VCN_by_VIS,width=6,height=2)
