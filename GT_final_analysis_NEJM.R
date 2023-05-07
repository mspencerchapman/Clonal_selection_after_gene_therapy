##------------------------------ANALYSIS-----------------------------

#This script goes through several key areas of the analysis
#1. Data summary/ technical aspects
#2. SNV Mutation burdens
#3. INDEL Mutation burdens
#4. Copy number landscape
#5. Structural variant landscape
#6. Driver mutation landscape
#7. dN/dS - for evidence of positive selection
#8. Vector copy number
root_dir="~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
setwd(root_dir)
plots_dir=paste0(root_dir,"/plots/")
#sample_metadata<-read_delim(paste0(root_dir,"/Data/sample_metadata_full.tsv"))
sample_metadata<-readRDS(paste0(root_dir,"/Data/sample_metadata_full_with_sigs.Rds"))
all_tree_data=readRDS(paste0(root_dir,"/Data/combined_tree_files.Rds"))
all_mut_data=readRDS(paste0(root_dir,"/Data/combined_muts_files.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(all_tree_data)) {assign(x,all_tree_data[[x]])}
for(x in names(all_mut_data)) {assign(x,all_mut_data[[x]])}

#Manually enter the individual-level metadata data frame
exp_IDs<-c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
Individual_metadata=data.frame(ID=exp_IDs,
                               new_ID=c("SCD4","SCD6","SCD5","SCD3","SCD1","SCD2"),
                               Disease=rep("SCD",6),
                               Age_at_GT=c(20,26,24,16,7,13),
                               PD_no=c("PD49229","PD53373b","PD49228","PD53374b","PD49227","PD49226"),
                               CD34_dose=c(5.07,NA,5.15,8.26,4.86,3.55))

patient_cols<-RColorBrewer::brewer.pal(7,"Set1")[c(1:5,7)]
names(patient_cols)<-paste0("SCD",1:6)

##########################################################################
###----------------------------DATA SUMMARY/ TECHNICAL--------------------
##########################################################################

#Review the coverage histogram across samples
coverage_summary<-sample_metadata%>%
  filter(!is.na(Coverage) &Coverage>4)%>%
  group_by(new_ID)%>%
  dplyr::summarise(n=n(),mean_coverage=mean(Coverage),median_coverage=median(Coverage))

coverage_histograms<-sample_metadata%>%
  filter(!is.na(Coverage) &Coverage>4)%>%
  ggplot(aes(x=Coverage))+
  geom_histogram(fill="lightblue",col="black",size=0.1,bins=30)+
  scale_x_continuous(limits=c(0,49))+
  facet_wrap(~new_ID,scales="free")+
  geom_vline(data=coverage_summary,aes(xintercept=mean_coverage),col="red",linetype=2)+
  geom_text(data=coverage_summary,aes(label=paste0("Mean coverage: ",signif(mean_coverage,digits=3)," x"),x=25,y=80),size=1.8,inherit.aes = F)+
  theme_classic()+
  my_theme+
  labs(y="Count")

ggsave(filename=paste0(plots_dir,"Coverage_histograms.pdf"),coverage_histograms,width=3.5,height=2.5)

sample_metadata%>%
  filter(sample_status=="PASS")%>%
  group_by(new_ID,HSPC_source,Sample_type,Time_point)%>%
  summarise(n=n())

#Review relationship between clonality ('peak VAF') and mutation burden
sample_metadata%>%
  filter(!is.na(peak_vaf) & !is.na(SNV_burden_AR))%>%
  filter(Coverage>8)%>%
  ggplot(aes(x=peak_vaf,y=SNV_burden_AR,col=ID))+
  geom_point(alpha=0.3,size=0.3)+
  scale_x_continuous(limits=c(0.3,0.7))+
  scale_y_continuous(limits=c(50,700))+
  facet_grid(Time_point~ID)+
  theme_bw()+
  geom_smooth(method="lm",col="black",size=0.2)+
  my_theme

#Summarize the outcomes of colonies by individual
colony_outcome_summary_plot<-sample_metadata%>%
  filter(sample_status!="Not Sequenced")%>%
  mutate(Time_point=ifelse(Time_point==0,"Pre-GT","Post-GT"))%>%
  ggplot(aes(x=factor(Time_point,levels=c("Pre-GT","Post-GT")),fill=sample_status))+
  geom_bar(col="black",size=0.1)+
  theme_classic()+
  facet_grid(cols=vars(new_ID))+
  scale_fill_brewer(palette="Paired")+
  scale_y_continuous(breaks=seq(0,1000,200))+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="",y="Count",fill="Outcome")

ggsave(filename = paste0(plots_dir,"colony_outcome_summary.pdf"),colony_outcome_summary_plot,width=3.5,height=2.5)

sample_metadata%>%
  mutate(sample_status=factor(sample_status,levels=rev(c("PASS","Low coverage","Not Sequenced","Non-clonal","Duplicate"))))%>%
  filter(sample_status!="Not Sequenced")%>%
  group_by(sample_status)%>%
  dplyr::summarise(n=n())

##########################################################################
###----------------------------SNV MUTATION BURDENS---------------------
##########################################################################

#First look at the pre-GT SNV burdens and compare to reference data set

#This is the regression from Mitchell et al, 2022
HSC_prediction=function(age){
  pred=54.6+16.8*age; upper_95_CI=54.6+16.95*age; lower_95_CI=54.6+16.65*age
  return(list(mean=pred,upper_95_CI=upper_95_CI,lower_95_CI=lower_95_CI))
}
lm.data=data.frame(x=0:26,ymin=unlist(sapply(0:26,HSC_prediction)['lower_95_CI',]),ymax=unlist(sapply(0:26,HSC_prediction)['upper_95_CI',]))

#Create a summary table of each individual at each time point 
#NB. include only the higher coverage samples (>10X) & those with a peak VAF > 0.45 (i.e. the most clonal samples)
Coverage_cutoff_for_mutburden_analysis=10
Peak_vaf_cutoff_for_mutburden_analysis=0.46
SNV_metric_for_mutburden_analysis="SNV_burden_AR"

mut_burden_summary<-sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)))%>%
  dplyr::filter(Sample_type!="DP")%>%
  dplyr::filter(Coverage>Coverage_cutoff_for_mutburden_analysis & peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  group_by(new_ID,Time_point)%>%
  dplyr::summarise(mean_SNV_burden=mean(get(SNV_metric_for_mutburden_analysis)),
            var=var(get(SNV_metric_for_mutburden_analysis)),
            mean_SNV_burden_invitro_removed=mean(SNV_burden_invitro_removed),
            var_invitro_removed=var(SNV_burden_invitro_removed),
            Age_at_sampling=Age_at_sampling[1],
            Disease=Disease[1])%>%
  mutate(expected_mean_SNV_burden=sapply(Age_at_sampling,function(Age) HSC_prediction(Age)$mean))%>%
  mutate(excess_mutations=mean_SNV_burden-expected_mean_SNV_burden#,
         #excess_mutations_invitro_removed=mean_SNV_burden_invitro_removed-expected_mean_SNV_burden
         )%>%
  mutate(excess_mutations_per_year=excess_mutations/Age_at_sampling,
         percentage_excess=mean_SNV_burden/expected_mean_SNV_burden#,
         #excess_mutations_per_year_invitro_removed=excess_mutations_invitro_removed/Age_at_sampling
         )

preGT_mut_burden_plot<-sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)))%>%
  dplyr::filter(Coverage>Coverage_cutoff_for_mutburden_analysis & peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  dplyr::filter(Disease=="SCD" & Time_point==0 & Sample_type!="DP")%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  ggplot(aes(x=Age_at_sampling,y=get(SNV_metric_for_mutburden_analysis),col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.1,size=0.3)+
  scale_y_continuous(limits=function(x) c(0,750))+
  scale_x_continuous(limits=function(x) c(0,x[2]))+
  scale_color_manual(values=patient_cols)+
  geom_abline(intercept = 54.6,slope = 16.8,size=0.3)+
  geom_ribbon(data=lm.data,aes(x=x,ymin=ymin, ymax=ymax),fill="darkgrey",alpha=0.7,inherit.aes = F)+
  geom_point(data=mut_burden_summary%>%filter(Disease=="SCD",Time_point==0),aes(x=Age_at_sampling,y=mean_SNV_burden),shape=3,size=1,inherit.aes = F)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme_bw()+
  labs(x="Age at sampling",
       y="Single nucleotide variant burden",
       col="Individual")+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),legend.box.spacing = unit(0,"mm"))

ggsave(filename = paste0(plots_dir,"preGT_mut_burden.pdf"),preGT_mut_burden_plot,width=2.8,height=2.5)

mut_burden_summary%>%filter(Time_point==0)

#Do a plot where all attempts have been made to remove invitro mutations
#(by refitting signatures to samples & removing all putative invitro signature mutations)
sample_metadata%>%
  ggplot(aes(x=SNV_burden_AR,y=SNV_burden_invitro_removed,col=new_ID))+
  geom_point(size=0.1,alpha=0.5)+
  scale_color_manual(values=patient_cols)+
  theme_classic()+
  labs(x="SNV burden (corrected)",
       y="SNV burden\n(in vitro signatures removed)",
       col="Individual")

preGT_mut_burden_invitro_removed_plot<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_invitro_removed))%>%
  dplyr::filter(Coverage>Coverage_cutoff_for_mutburden_analysis & peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  dplyr::filter(Disease=="SCD" & Time_point==0 & Sample_type!="DP")%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  ggplot(aes(x=Age_at_sampling,y=SNV_burden_invitro_removed,col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.1,size=0.3)+
  scale_y_continuous(limits=function(x) c(0,750))+
  scale_x_continuous(limits=function(x) c(0,x[2]))+
  scale_color_manual(values=patient_cols)+
  geom_abline(intercept = 54.6,slope = 16.8,size=0.3)+
  geom_ribbon(data=lm.data,aes(x=x,ymin=ymin, ymax=ymax),fill="darkgrey",alpha=0.7,inherit.aes = F)+
  geom_point(data=mut_burden_summary%>%filter(Disease=="SCD",Time_point==0),aes(x=Age_at_sampling,y=mean_SNV_burden_invitro_removed),shape=3,size=1,inherit.aes = F)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme_bw()+
  labs(x="Age at sampling",
       y="Single nucleotide variant burden\n(corrected for invitro mutations)",
       col="Individual")+
  my_theme

ggsave(filename = paste0(plots_dir,"preGT_mut_burden_invitro_removed.pdf"),preGT_mut_burden_invitro_removed_plot,width=3,height=2)

### REVIEW IF ANY RELATIONSHIP BETWEEN VCN and MUTATION BURDEN
sample_metadata%>%
  filter(Sample_type=="Post-GT")%>%
  ggplot(aes(x=VCN_rounded,y=get(SNV_metric_for_mutburden_analysis),col=new_ID))+
  geom_jitter(width=0.2,alpha=0.5,size=0.1)+
  scale_color_manual(values=patient_cols)+
  facet_grid(~new_ID)+
  geom_smooth(method="lm")+
  theme_bw()+
  my_theme+
  labs(x="Vector copy number",y="SNV burden",col="Individual")

##Now do comparison of just transduced vs untransduced
transduced_vs_untransduced_df<-sample_metadata%>%
  filter(Sample_type=="Post-GT" & sample_status=="PASS")%>%
  mutate(Transduced=ifelse(VCN_rounded>0.8,"Gene\nmod.","Non\nmod."),Time_point=paste(Time_point,"years"))

transduced_vs_untransduced_comparison<-lapply(Individual_metadata$new_ID,function(exp_ID) {
  time_points=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & !is.na(get(SNV_metric_for_mutburden_analysis)))%>%pull(Time_point)%>%unique()
  lapply(time_points,function(time_point){
    test=t.test(x=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Gene\nmod.")%>%pull(get(SNV_metric_for_mutburden_analysis)),
                y=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Non\nmod.")%>%pull(get(SNV_metric_for_mutburden_analysis)))
    return(data.frame(new_ID=exp_ID,Time_point=time_point,estimate=test$estimate[1]-test$estimate[2],p.value=test$p.value,CI_low=test$conf.int[1],CI_high=test$conf.int[2]))
  })%>%bind_rows()
})%>%bind_rows

transduced_vs_nontransduced_plot<-sample_metadata%>%
  filter(Sample_type=="Post-GT" & sample_status=="PASS" & !is.na(get(SNV_metric_for_mutburden_analysis)))%>%
  mutate(Transduced=ifelse(VCN_rounded>0.8,"Gene\nmod.","Non\nmod."),Time_point=paste(Time_point,"years"))%>%
  ggplot(aes(x=Transduced,y=get(SNV_metric_for_mutburden_analysis)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.15,alpha=0.2,size=0.5)+
  facet_wrap(Time_point~new_ID)+
  geom_text(data=transduced_vs_untransduced_comparison,aes(label=paste("p =",round(p.value,digits=2)),x=1.5,y=300),size=2.5,fontface="italic",inherit.aes = F)+
  geom_smooth(method="lm")+
  scale_y_continuous(limits=c(250,700))+
  theme_bw()+
  labs(x="",y="SNV burden",col="Individual")+
  my_theme

ggsave(filename = paste0(plots_dir,"SNV_burden_transduced_vs_nontransduced_plot.pdf"),plot=transduced_vs_nontransduced_plot,width=2,height=2)

### COMPARISON OF PRE- and POST-GENE THERAPY

#Perform mutation burden comparisons only using the higher coverage samples with best evidence of clonality
SNV_metric_for_mutburden_analysis="SNV_burden_tc"
pre_post_comparison_df<-sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis))&
                  sample_status=="PASS"&
                  Sample_type!="DP"&
                  Coverage>Coverage_cutoff_for_mutburden_analysis&
                  peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  mutate(SNV_adjusted=get(SNV_metric_for_mutburden_analysis)-(16.8*Time_point))
dim(pre_post_comparison_df)

pre_post_comparison_df%>%
  group_by(new_ID,Time_point)%>%
  summarize(mean_SNV_burden=mean(get(SNV_metric_for_mutburden_analysis)))

pre_post_comparison_df%>%
  filter(Sample_type!="DP")%>%
  ggplot(aes(x=factor(Time_point),y=get(SNV_metric_for_mutburden_analysis),col=Cell_type))+
  geom_boxplot(aes(group=plate_ID))+
  geom_jitter(width = 0.1,alpha=0.2)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  labs(x="Time point (months)", y="SNV burden")+
  theme(legend.text = element_text(size=7),legend.key.size = unit(3,"mm"))

pre_vs_post_comparison_plot<-pre_post_comparison_df%>%
  ggplot(aes(x=Time_point,y=get(SNV_metric_for_mutburden_analysis)))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size = 0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.05,alpha=0.15,size=0.5,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_color_brewer(palette="Dark2")+
  scale_x_continuous(breaks=seq(0,3,1))+
  labs(x="Time point (years)", y="SNV burden")+
  my_theme

all.lms<-lapply(paste0("SCD",1:6),function(this_ID) {
  cat(this_ID,sep="\n")
  res<-lm(get(SNV_metric_for_mutburden_analysis)~Time_point,data=pre_post_comparison_df%>%filter(new_ID==this_ID))
  return(res)
})

pre_vs_post_comparison_adjusted_plot<-pre_post_comparison_df%>%
  mutate(Time_point=Time_point)%>%
  ggplot(aes(x=Time_point,y=SNV_adjusted))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size=0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.025,alpha=0.15,size=0.5,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(250,750))+
  scale_x_continuous(breaks=seq(0,3,1))+
  labs(x="Time point (years)", y="Age-adjusted SNV burden")+
  my_theme

all.lms.adjusted<-lapply(paste0("SCD",1:6),function(this_ID) {
  cat(this_ID,sep="\n")
  res<-lm(SNV_adjusted~Time_point,data=pre_post_comparison_df%>%filter(new_ID==this_ID))
  return(res)
})

lapply(all.lms.adjusted,confint)

ggsave(paste0(plots_dir,"pre_vs_post_comparison_plot.pdf"),plot=pre_vs_post_comparison_plot,width=3.5,height=2)
ggsave(paste0(plots_dir,"pre_vs_post_comparison_adjusted_plot.pdf"),plot=pre_vs_post_comparison_adjusted_plot,width=3.5,height=2)

pre_post_comparison_df%>%
  ggplot(aes(x=factor(Time_point),y=SNV_adjusted,col=Cell_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.15,alpha=0.3,size=1,stroke=0)+
  facet_grid(cols=vars(ID),scales="free",space = "free")+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(250,750))+
  labs(x="Time point (months)", y="Age-adjusted SNV burden")


t.test.all<-lapply(paste0("SCD",1:6),function(this_ID) {
  t.test.res<-t.test(y=pre_post_comparison_df%>%filter(new_ID==this_ID & Time_point==0)%>%pull(SNV_adjusted),
         x=pre_post_comparison_df%>%filter(new_ID==this_ID & Time_point>0)%>%pull(SNV_adjusted),
         alternative="two.sided")
  return(t.test.res)
})

gt_induced_mutations_estimate_plot<-Map(test=t.test.all,exp_ID=paste0("SCD",1:6),function(test,exp_ID){
  data.frame(ID=exp_ID,estimate=test$estimate[1]-test$estimate[2],CI_low=test$conf.int[1],CI_high=test$conf.int[2])
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=ID,col=ID,y=estimate,ymin=CI_low,ymax=CI_high))+
  geom_point(size=0.6)+
  scale_color_manual(values=patient_cols)+
  geom_errorbar(width=0.3,alpha=0.5)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-50,50),breaks=seq(-50,50,10))+
  theme_bw()+
  labs(x="",y="Estimate of excess SNVs\nfrom gene therapy")+
  my_theme+
  theme(legend.box.spacing=unit(0,"mm"),legend.key.size = unit(3,"mm"),axis.text.x=element_text(angle=90))

ggsave(paste0(plots_dir,"gt_induced_mutations.pdf"),plot=gt_induced_mutations_estimate_plot,width=1.8,height=2)

#Histograms of mutation burdens pre & post
sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)) & new_ID%in%c("SCD3","SCD4") & Coverage>10 &peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  ggplot(aes(x=get(SNV_metric_for_mutburden_analysis)))+
  geom_histogram()+
  facet_grid(Time_point~new_ID)

#Print table of average burdens pre & post
sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)) & new_ID%in%c("SCD3","SCD4") & Coverage>10 &peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  group_by(ID,Time_point)%>%
  dplyr::summarise(median=median(get(SNV_metric_for_mutburden_analysis)),mean=mean(get(SNV_metric_for_mutburden_analysis)),var=var(get(SNV_metric_for_mutburden_analysis)))

##18. Print table of average burdens in pre-GT/pre-TDX & DP (where available) and do t.test comparison
sample_metadata%>%
  dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)) & new_ID%in%c("SCD1") & Coverage>10)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  group_by(ID,Sample_type)%>%
  dplyr::summarise(median=median(get(SNV_metric_for_mutburden_analysis)),mean=mean(get(SNV_metric_for_mutburden_analysis)),var=var(get(SNV_metric_for_mutburden_analysis)))

#t.test comparison of pre-TDX and donor product sample in BCL009
wilcox.test(x=sample_metadata%>%
         dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)) & new_ID%in%c("SCD1") & Coverage>10 & Sample_type=="Pre-GT")%>%
         pull(get(SNV_metric_for_mutburden_analysis)),
       y=sample_metadata%>%
         dplyr::filter(!is.na(get(SNV_metric_for_mutburden_analysis)) & new_ID%in%c("SCD1") & Coverage>10 & Sample_type%in%c("Pre-TDX","DP"))%>%
         pull(get(SNV_metric_for_mutburden_analysis)))

##########################################################################
###----------------------------INDEL MUTATION BURDENS---------------------
##########################################################################

sample_metadata<-sample_metadata%>%
  mutate(INDEL_burden_gc=INDEL_burden_raw/germline_INDEL_sensitivity)

#First look at the pre-GT SNV burdens and compare to reference data set
Coverage_cutoff_for_mutburden_analysis=10
Peak_vaf_cutoff_for_mutburden_analysis=0.46

INDEL_metric_for_mutburden_analysis="INDEL_burden_AR"

#This is the regression from Mitchell et al, 2022
HSC_prediction.INDEL=function(age){
  pred=3.92+0.72*age; upper_95_CI=3.92+0.7101811*age; lower_95_CI=3.92+0.7306062*age
  return(list(mean=pred,upper_95_CI=upper_95_CI,lower_95_CI=lower_95_CI))
}

lm.data.indel=data.frame(x=0:26,ymin=unlist(sapply(0:26,HSC_prediction.INDEL)['lower_95_CI',]),ymax=unlist(sapply(0:26,HSC_prediction.INDEL)['upper_95_CI',]))

#Create a summary table of each individual at each time point 
INDEL_mut_burden_summary<-sample_metadata%>%
  dplyr::filter(!is.na(get(INDEL_metric_for_mutburden_analysis)))%>%
  dplyr::filter(Sample_type!="DP"&
                  Coverage>Coverage_cutoff_for_mutburden_analysis &
                  peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  group_by(new_ID,Time_point)%>%
  dplyr::summarise(mean_INDEL_burden=mean(get(INDEL_metric_for_mutburden_analysis)),
            var=var(get(INDEL_metric_for_mutburden_analysis)),
            Age_at_sampling=Age_at_sampling[1],
            Disease=Disease[1])%>%
  mutate(expected_mean_INDEL_burden=sapply(Age_at_sampling,function(Age) HSC_prediction.INDEL(Age)$mean))%>%
  mutate(excess_mutations=mean_INDEL_burden-expected_mean_INDEL_burden)%>%
  mutate(excess_mutations_per_year=excess_mutations/Age_at_sampling)

INDEL_mut_burden_summary%>%dplyr::filter(Time_point==0)

preGT_INDEL_mut_burden_plot<-sample_metadata%>%
  dplyr::filter(!is.na(get(INDEL_metric_for_mutburden_analysis)))%>%
  dplyr::filter(Sample_type!="DP"&
                  Coverage>Coverage_cutoff_for_mutburden_analysis&
                  peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  dplyr::filter(Disease=="SCD" & Time_point==0)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  ggplot(aes(x=Age_at_sampling,y=get(INDEL_metric_for_mutburden_analysis),col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.1,size=0.3)+
  scale_y_continuous(limits=function(x) c(0,30))+
  scale_x_continuous(limits=function(x) c(0,x[2]))+
  scale_color_manual(values=patient_cols)+
  geom_abline(intercept = 54.6,slope = 16.8,size=0.3)+
  geom_ribbon(data=lm.data.indel,aes(x=x,ymin=ymin, ymax=ymax),fill="darkgrey",alpha=0.7,inherit.aes = F)+
  geom_point(data=INDEL_mut_burden_summary%>%filter(Disease=="SCD",Time_point==0),aes(x=Age_at_sampling,y=mean_INDEL_burden),shape=3,size=1,inherit.aes = F)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme_bw()+
  labs(x="Age at sampling",
       y="Indel burden",
       col="Individual")+
  my_theme

ggsave(filename = paste0(plots_dir,"preGT_INDEL_mut_burden.pdf"),preGT_INDEL_mut_burden_plot,width=3,height=2)

### REVIEW IF ANY RELATIONSHIP BETWEEN VCN and MUTATION BURDEN
INDEL_metric_for_mutburden_analysis="INDEL_burden_AR"
sample_metadata%>%
  dplyr::filter(Coverage>Coverage_cutoff_for_mutburden_analysis & peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  filter(Sample_type=="Post-GT")%>%
  ggplot(aes(x=VCN_rounded,y=get(INDEL_metric_for_mutburden_analysis),col=new_ID))+
  geom_jitter(width=0.2,alpha=0.5,size=0.1)+
  facet_grid(~new_ID)+
  geom_smooth(method="lm")+
  theme_bw()+
  my_theme+
  labs(x="Vector copy number",y="SNV burden",col="Individual")

##Now do comparison of just transduced vs untransduced
transduced_vs_untransduced_comparison.indel<-lapply(paste0("SCD",1:6),function(exp_ID) {
  cat(exp_ID,sep="\n")
  time_points=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & !is.na(get(INDEL_metric_for_mutburden_analysis)))%>%pull(Time_point)%>%unique()
  lapply(time_points,function(time_point){
    cat(time_point,sep="\n")
    test=t.test(x=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Gene\nmod.")%>%pull(get(INDEL_metric_for_mutburden_analysis)),
                y=transduced_vs_untransduced_df%>%filter(new_ID==exp_ID & Time_point==time_point & Transduced=="Non\nmod.")%>%pull(get(INDEL_metric_for_mutburden_analysis)))
    return(data.frame(new_ID=exp_ID,Time_point=time_point,estimate=test$estimate[1]-test$estimate[2],p.value=test$p.value,CI_low=test$conf.int[1],CI_high=test$conf.int[2]))
  })%>%bind_rows()
})%>%bind_rows

transduced_vs_nontransduced_plot.indel<-sample_metadata%>%
  dplyr::filter(Coverage>Coverage_cutoff_for_mutburden_analysis & peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  filter(Sample_type=="Post-GT" & sample_status=="PASS" & !is.na(get(INDEL_metric_for_mutburden_analysis)))%>%
  mutate(Transduced=ifelse(VCN_rounded>0.8,"Gene\nmod.","Non\nmod."),Time_point=paste(Time_point,"years"))%>%
  ggplot(aes(x=Transduced,y=get(INDEL_metric_for_mutburden_analysis)))+
  geom_boxplot(outlier.shape=NA,size=0.3)+
  geom_jitter(aes(col=new_ID),width=0.15,alpha=0.3,size=0.2)+
  facet_wrap(Time_point~new_ID,nrow=2)+
  geom_text(data=transduced_vs_untransduced_comparison.indel,aes(label=paste("p =",round(p.value,digits=2)),x=1.5,y=10),size=1.8,fontface="italic",inherit.aes = F)+
  geom_smooth(method="lm")+
  scale_color_manual(values=patient_cols)+
  scale_y_continuous(limits=c(5,40))+
  theme_bw()+
  labs(x="",y="Indel burden",col="Individual")+
  my_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  theme(strip.text.x=element_text(size=4,margin = unit(c(1,0,1,0),"mm")),legend.key.size = unit(2.5,"mm"),legend.box.spacing = unit(0,"mm"))

ggsave(filename = paste0(plots_dir,"INDEL_burden_transduced_vs_nontransduced_plot.pdf"),plot=transduced_vs_nontransduced_plot.indel,width=4.5,height=2)

INDEL_metric_for_mutburden_analysis="INDEL_burden_gc"
INDEL_summary_df<-sample_metadata%>%
  dplyr::filter(Sample_type!="DP"&
                  Coverage>Coverage_cutoff_for_mutburden_analysis&
                  peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  group_by(ID,Time_point)%>%
  dplyr::summarise(mean=mean(get(INDEL_metric_for_mutburden_analysis),na.rm=T),median=median(get(INDEL_metric_for_mutburden_analysis),na.rm = T),var=var(get(INDEL_metric_for_mutburden_analysis),na.rm=T),sd=sd(get(INDEL_metric_for_mutburden_analysis),na.rm=T))%>%
  left_join(Individual_metadata)

pre_post_comparison_df.INDEL<-sample_metadata%>%
  dplyr::filter(Sample_type!="DP"&
                  !is.na(get(INDEL_metric_for_mutburden_analysis))&
                  Coverage>Coverage_cutoff_for_mutburden_analysis&
                  peak_vaf>Peak_vaf_cutoff_for_mutburden_analysis)%>%
  mutate(Age_at_sampling=(Time_point/12+Age_at_GT))%>%
  mutate(INDEL_adjusted=get(INDEL_metric_for_mutburden_analysis)-(0.72*Time_point))

pre_post_comparison_df.INDEL%>%
  group_by(ID,Time_point)%>%
  summarize(mean(get(INDEL_metric_for_mutburden_analysis)))

pre_post_comparison_df.INDEL%>%
  ggplot(aes(x=factor(Time_point),y=get(INDEL_metric_for_mutburden_analysis),col=Cell_type))+
  geom_boxplot(aes(group=plate_ID))+
  geom_jitter(width = 0.1,alpha=0.2)+
  #scale_color_manual(values=patient_cols)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  labs(x="Time point (years)", y="Indel burden")+
  theme(legend.text = element_text(size=7),legend.key.size = unit(3,"mm"))

pre_vs_post_comparison_plot.INDEL<-pre_post_comparison_df.INDEL%>%
  ggplot(aes(x=Time_point,y=get(INDEL_metric_for_mutburden_analysis)))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size = 0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.025,alpha=0.3,size=0.7,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_color_brewer(palette="Dark2")+
  scale_x_continuous(breaks=seq(0,3,1))+
  scale_y_continuous(limits=c(0,100))+
  labs(x="Time point (years)", y="Indel burden")+
  my_theme

lm(get(INDEL_metric_for_mutburden_analysis)~Time_point,data=pre_post_comparison_df.INDEL%>%filter(new_ID=="SCD3"))
lm(get(INDEL_metric_for_mutburden_analysis)~Time_point,data=pre_post_comparison_df.INDEL%>%filter(new_ID=="SCD4"))

pre_vs_post_comparison_adjusted_plot.INDEL<-pre_post_comparison_df.INDEL%>%
  ggplot(aes(x=Time_point,y=INDEL_adjusted))+
  geom_boxplot(aes(group=factor(Time_point)),outlier.size=0.4,size=0.25,outlier.shape = NA)+
  geom_jitter(width = 0.025,alpha=0.3,size=0.4,stroke=0)+
  facet_grid(cols=vars(new_ID),scales="free",space = "free")+
  theme_bw()+
  geom_smooth(method="lm",col="blue",lwd=0.3)+
  scale_y_continuous(limits=c(0,75))+
  scale_x_continuous(breaks=seq(0,3,1))+
  labs(x="Time point (years)", y="Age-adjusted Indel burden")+
  my_theme

all.lms.indel.adjust<-lapply(paste0("SCD",1:6),function(this_ID) {
  if(this_ID%in%pre_post_comparison_df.INDEL$new_ID){
    return(lm(INDEL_adjusted~Time_point,data=pre_post_comparison_df.INDEL%>%filter(new_ID==this_ID)))
  } else {
    return(NULL)
  }
})

ggsave(paste0(plots_dir,"pre_vs_post_comparison_plot.INDEL.pdf"),plot=pre_vs_post_comparison_plot.INDEL,width=3.5,height=2)
ggsave(paste0(plots_dir,"pre_vs_post_comparison_adjusted_plot.INDEL.pdf"),plot=pre_vs_post_comparison_adjusted_plot.INDEL,width=3.5,height=2)


pre_post_comparison_df.INDEL%>%
  ggplot(aes(x=factor(Time_point),y=INDEL_adjusted,col=Cell_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.3,alpha=0.3,size=1,stroke=0)+
  facet_grid(cols=vars(ID),scales="free",space = "free")+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  labs(x="Time point (months)", y="Age-adjusted SNV burden")

t.test.indels<-lapply(paste0("SCD",1:6),function(this_ID) {
  if(this_ID%in%pre_post_comparison_df.INDEL$new_ID){
    res.t.test<-t.test(y=pre_post_comparison_df.INDEL%>%filter(new_ID==this_ID&Time_point==0)%>%pull(INDEL_adjusted),
           x=pre_post_comparison_df.INDEL%>%filter(new_ID==this_ID&Time_point>0)%>%pull(INDEL_adjusted),
           alternative="two.sided")
    return(res.t.test)
  } else {
    return(NULL)
  }
})

gt_induced_mutations_estimate_plot.INDEL<-Map(test=t.test.indels,exp_ID=paste0("SCD",1:6),function(test,exp_ID){
  if(!is.null(test)){
    data.frame(ID=exp_ID,estimate=test$estimate[1]-test$estimate[2],CI_low=test$conf.int[1],CI_high=test$conf.int[2])
  } else {
    NULL
  }
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=ID,col=ID,y=estimate,ymin=CI_low,ymax=CI_high))+
  geom_point(size=0.6)+
  scale_color_manual(values=patient_cols)+
  geom_errorbar(width=0.2,alpha=0.5)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-6,6),breaks=seq(-6,6,2))+
  theme_bw()+
  labs(x="",y="Estimate of excess INDELs\nfrom gene therapy")+
  my_theme+
  theme(legend.box.spacing=unit(0,"mm"),legend.key.size = unit(3,"mm"),axis.text.x=element_text(angle=90))

ggsave(paste0(plots_dir,"pre_vs_post_comparison_adjusted_plot.INDEL.pdf"),plot=pre_vs_post_comparison_adjusted_plot.INDEL,width=2.5,height=2)
ggsave(paste0(plots_dir,"gt_induced_INDEL_mutations.pdf"),plot=gt_induced_mutations_estimate_plot.INDEL,width=1.8,height=2)


##########################################################################
###----------------------------COPY NUMBER LANDSCAPE---------------------
##########################################################################

CN_change_df=read.csv(paste0(root_dir,"/Data/Copy_number_changes.csv"))
CN_change_df<-CN_change_df%>%filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))

n_CN_events_df<-data.frame(ID=names(all.trees.cc.nodups),n_CNs=sapply(names(all.trees.cc.nodups),function(ID) sum(CN_change_df$Individual==ID)),total=sapply(all.trees.cc.nodups,function(tree) length(tree$tip.label)-1))%>%
  left_join(Individual_metadata)

CN_type_vec=c("Copy neutral\nloss of heterozygosity","Deletion","Duplication")
names(CN_type_vec)=c("UPD","del","dup")
CNs_per_colony<-CN_change_df%>%
  right_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  replace_na(replace=list(Class="No_CN"))%>%
  group_by(Class,new_ID)%>%
  dplyr::summarise(n=n())%>%
  pivot_wider(names_from="Class",values_from="n")%>%
  replace_na(replace=list(del=0,UPD=0,dup=0))%>%
  mutate(total=del+dup+UPD+No_CN)%>%
  mutate(del=del/total,dup=dup/total,UPD=UPD/total,No_CN=No_CN/total)%>%
  gather(-new_ID,-total,key="type",value="Events per colony")%>%
  filter(type!="No_CN")%>%
  mutate(type=CN_type_vec[type])%>%
  ggplot(aes(x=new_ID,y=`Events per colony`,fill=factor(type)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.3,width=0.75)+
  geom_text(data=n_CN_events_df,aes(x=new_ID,y=n_CNs/total,label=n_CNs),nudge_y = +0.0005,size=2,inherit.aes = F)+
  scale_y_continuous(limits=c(0,0.007),breaks=seq(0,0.006,0.002))+
  scale_fill_brewer(palette="BuGn")+
  theme_classic()+
  labs(x="",y="Average number of CNAs per colony",fill="")+
  my_theme+
  theme(axis.text.x=element_text(angle=90))

ggsave(filename=paste0(plots_dir,"CNs_per_colony.pdf"),CNs_per_colony,width=2.8,height=2)

#Plot of copy number changes by individual
CN_changes_per_individual_plot<-CN_change_df%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  mutate(Time_point=ifelse(Time_point==0,"Pre-GT","Post-GT"))%>%
  ggplot(aes(x=factor(Time_point,levels=c("Pre-GT","Post-GT")),y=1,fill=factor(Copy_number_change)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15,width=0.5)+
  facet_grid(cols=vars(factor(new_ID)),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.key.size=unit(2.5,"mm"))+
  labs(x="Sample type",y="Samples",fill="Copy number abnormality")

ggsave(filename=paste0(plots_dir,"CN_changes_by_individual.pdf"),CN_changes_per_individual_plot,width=4,height=2)

#Proportion of colonies with a CNA in different sample groups - i.e. pre-GT vs post-GT (excluding LOY)
right_join(CN_change_df,sample_metadata)%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  replace_na(replace=list(Copy_number_change="Normal"))%>%
  mutate(CNA=ifelse(Copy_number_change=="Normal","Normal","Abnormal"))%>%
  group_by(CNA,Time_point>0)%>%
  dplyr::summarise(n=n())%>%
  mutate(Time_point=ifelse(`Time_point > 0`,"Post-GT","Pre-GT"))%>%
  pivot_wider(names_from="CNA",values_from="n")%>%
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


##########################################################################
###----------------------------STRUCTURAL VARIANT LANDSCAPE---------------
##########################################################################

SV_data<-readxl::read_excel(path=paste0(root_dir,"/Data/SV_analysis_sickle_updated050123.xlsx"),sheet=2)

#Filter out events in non-clonal samples
SV_data<-SV_data%>%
  filter(sample%in%(sample_metadata%>%dplyr::filter(sample_status!="Non-clonal")%>%pull(Sample)))%>%
  dplyr::rename("Sample"=sample)

#Filter single events that have two separate calls in the same sample (e.g. relating to the two ends of the change)
SV_change_df<-SV_data%>%
  dplyr::filter(!(SV_data%>%dplyr::select(type,Sample,gene1)%>%duplicated()))%>%
  mutate(event_length=ifelse(chro1==chro2,start2-start1,NA))

SV_tables=lapply(all.trees.uncorrected,function(tree) {
  SV_change_Ind<-SV_change_df%>%filter(Sample%in%tree$tip.label)
  SV_change_Ind$node<-sapply(SV_change_Ind$Sample,function(SampleID) which(tree$tip.label==SampleID))
  return(SV_change_Ind)
  })

SV_tables=Map(SV_table=SV_tables,tree=all.trees.cc,function(SV_table,tree){
  duplicate_sets=get_duplicate_sets(tree = tree,mut_threshold = 80)
  SV_table$duplicate_status=sapply(SV_table$Sample,function(SampleID) {
    which_set=sapply(duplicate_sets,function(set) SampleID %in% set)
    if(any(which_set)) {
      duplicate_set=duplicate_sets[which_set]
      if(all(duplicate_set[[1]]%in%SV_table$Sample)) {
        return("All_dups")
      } else {
        return("Not_all_dups")
      }
    } else {
      return("No_duplicates")
    }
  })
  return(SV_table)
})
SV_invitro=bind_rows(lapply(SV_tables,function(SV_table) SV_table%>%filter(duplicate_status=="Not_all_dups")%>%mutate(duplicate_status=as.character(duplicate_status))))
SV_tables=lapply(SV_tables,function(SV_table) SV_table%>%filter(duplicate_status!="Not_all_dups")%>%mutate(duplicate_status=as.character(duplicate_status)))

#Now reassign to node numbers from the 'nodups' tree
SV_tables=Map(SV_table=SV_tables,tree=all.trees.cc.nodups,function(SV_table,tree) {
  SV_table$node<-sapply(SV_table$Sample,function(SampleID) if(!SampleID%in%tree$tip.label) {return(NA)} else {return(which(tree$tip.label==SampleID))})
  return(SV_table%>%filter(!is.na(node))%>%arrange(node))
})

#Compressed shared SVs into a single event occurring on parental node for tree plotting
SV_tables_for_tree<-Map(df=SV_tables,tree=all.trees.cc.nodups,function(df,tree) {
  i=1
  while(i<nrow(df)) {
    if(df$chro1[i] == df$chro1[i+1] &
       df$start1[i] == df$start1[i+1] &
       (df$node[i]+1) == df$node[i+1]) {
      print(paste("Merging into single event:",df$chro1[i],df$start1[i],df$end1[i],df$type[i]))
      df$node[i]<-get_ancestor_node(node=df$node[i],tree=tree)
      df$Sample[i]<-paste(df$Sample[i],df$Sample[i+1],sep=",")
      df<-df%>%filter(!row_number() == (i+1))
    }
    i<-i+1
  }
  return(df)
})

pdf(paste0(plots_dir,"/SV_phylos.pdf"),width=7,height=2.5)
temp=Map(df=SV_tables_for_tree,tree=all.trees.cc.nodups,function(df,tree) {
  df$label=sapply(1:nrow(df),function(i) {
    if(df$type[i]=="CTX") {
      return(paste0("t(",df$chro1[i],";",df$chro2[i],") translocation"))
    } else {
      return(paste0("Chr",df$chro1[i],", ",df$type[i],", ",round(df$event_length[i]/1000)," Kb"))
    }
  })
  tree=plot_tree(tree,cex.label = 0)
  temp=plot_tree_labels(tree,
                        df,
                        type="line",
                        query.field = "type",
                        query.allowed.df=data.frame(value=c("DEL","DUP","INV","CTX"),col=brewer.pal(4,"Set1"),pch = 17,stringsAsFactors = FALSE), #if use 'coding_change_chip', value is 'Coding change mutation in driver'
                        label.field = "label",
                        cex.label = 0.4,
                        lty=1,
                        lwd=2)
})
dev.off()

SV_type_vec=c("Translocation","Inversion","Duplication","Deletion")
names(SV_type_vec)=c("CTX","INV","DUP","DEL")

SV_by_individual_plot<-bind_rows(SV_tables)%>%
  left_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  replace_na(list(gene1="Intergenic"))%>%
  mutate(type=SV_type_vec[type])%>%
  ggplot(aes(x=factor(Sample_type,levels=c("Pre-GT","Pre-TDX","Post-GT","DP")),y=1,col=type,fill=event_length))+
  geom_bar(stat="identity",position="stack",size=0.8,width=0.5)+
  facet_grid(cols=vars(factor(new_ID)),drop=F)+
  scale_color_brewer(palette="Spectral")+
  scale_fill_continuous()+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="Sample type",y="Samples",fill="SV length",col="SV type")

#Proportion of colonies with an SV in different sample groups - i.e. Individual & pre-GT vs post-GT (excluding LOY)
SVs_per_colony_pre_vs_post<-bind_rows(SV_tables)%>%
  right_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  replace_na(replace=list(type="No_SV"))%>%
  mutate(SV_status=ifelse(type=="No_SV","No_SV","SV_present"))%>%
  group_by(new_ID,SV_status,Time_point>0)%>%
  dplyr::summarise(n=n())%>%
  mutate(Time_point=ifelse(`Time_point > 0`,"Post-GT","Pre-GT"))%>%
  pivot_wider(names_from="SV_status",values_from="n")%>%
  replace_na(replace=list(SV_present=0))%>%
  mutate(abnormal_prop=SV_present/(SV_present+No_SV),
         CI_low=mapply(x=SV_present,n=(SV_present+No_SV),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=SV_present,n=(SV_present+No_SV),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=factor(Time_point,levels=c("Pre-GT","Post-GT")),col=factor(Time_point,levels=c("Pre-GT","Post-GT")),y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  theme_classic()+
  facet_grid(cols=vars(new_ID))+
  scale_y_continuous(limits=c(0,8))+
  theme(legend.position = "none",axis.text.x = element_text(angle=90))+
  labs(x="",y="Proportion of colonies with SV (%)")+
  my_theme

ggsave(filename=paste0(plots_dir,"SVs_per_colony_pre_vs_post.pdf"),SVs_per_colony_pre_vs_post,width=2.5,height=2)

#Proportion of colonies with an SV in Individuals
bind_rows(SV_tables)%>%
  right_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  replace_na(replace=list(type="No_SV"))%>%
  mutate(SV_status=ifelse(type=="No_SV","No_SV","SV_present"))%>%
  group_by(SV_status,new_ID)%>%
  dplyr::summarise(n=n())%>%
  pivot_wider(names_from="SV_status",values_from="n")%>%
  replace_na(replace=list(SV_present=0,No_SV=0))%>%
  mutate(abnormal_prop=SV_present/(SV_present+No_SV),
         CI_low=mapply(x=SV_present,n=(SV_present+No_SV),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=SV_present,n=(SV_present+No_SV),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=new_ID,y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  theme_classic()+
  scale_y_continuous(limits=c(0,7))+
  theme(legend.position = "none")+
  labs(x="",y="Proportion of colonies with SV (%)")

#Proportion of colonies with an SV in Individuals
n_events_df<-data.frame(ID=names(SV_tables),n_SVs=sapply(SV_tables,nrow),total=sapply(all.trees.cc.nodups,function(tree) length(tree$tip.label)-1))%>%
  left_join(Individual_metadata)

SVs_per_colony<-bind_rows(SV_tables)%>%
  right_join(sample_metadata,by="Sample")%>%
  filter(sample_status=="PASS")%>%
  replace_na(replace=list(type="No_SV"))%>%
  group_by(type,new_ID)%>%
  dplyr::summarise(n=n())%>%
  pivot_wider(names_from="type",values_from="n")%>%
  replace_na(replace=list(CTX=0,DEL=0,DUP=0,INV=0,No_SV=0))%>%
  mutate(total=CTX+DEL+DUP+INV+No_SV)%>%
  mutate(CTX=CTX/total,DEL=DEL/total,DUP=DUP/total,INV=INV/total,No_SV=No_SV/total)%>%
  gather(-new_ID,-total,key="type",value="Events per colony")%>%
  filter(type!="No_SV")%>%
  mutate(type=SV_type_vec[type])%>%
  ggplot(aes(x=new_ID,y=`Events per colony`,fill=factor(type)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.3,width=0.75)+
  geom_text(data=n_events_df,aes(x=new_ID,y=n_SVs/total,label=n_SVs),nudge_y = +0.002,size=2,inherit.aes = F)+
  scale_y_continuous(limits=c(0,0.04),breaks=seq(0,0.04,0.01))+
  scale_fill_brewer(palette="RdPu")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="",y="Average number of SVs per colony",fill="")+
  my_theme

ggsave(filename=paste0(plots_dir,"SVs_per_colony.pdf"),SVs_per_colony,width=2.8,height=2)


##########################################################################
###----------------------------DRIVER MUTATION LANDSCAPE---------------
##########################################################################
annotated_driver_mut_set=read_csv(paste0(root_dir,"/Data/Possible_drivers_annotated.csv"))#[,c("mut_ref","Gene","variant_ID","Decision")]

all.res_df<-Map(tree=all.trees.cc.nodups,details=all.muts.nodups,function(tree,details) {
  details$Decision<-NULL
  details<-left_join(details,annotated_driver_mut_set%>%dplyr::select(mut_ref,Decision),by="mut_ref")
  details$is.driver<-ifelse(details$Decision=="Probable"|details$Decision=="Possible",1,0)
  res_df<-drivers_per_sample_data(tree,details)
  sumstat=sapply(c("≥1 driver"=1,"≥2 drivers"=2,"≥3 drivers"=3,"≥4 drivers"=4,"≥5 drivers"=5),function(n) {sum(res_df$n_drivers>=n)/nrow(res_df)})
  return(res_df)
})

driver_proportion_by_individual_plot<-Map(df=all.res_df,Exp_ID=names(all.res_df),f=function(df,Exp_ID){
  prop_drivers_df<-df%>%
    left_join(sample_metadata%>%dplyr::select(Sample,Sample_type,Time_point,new_ID),by="Sample")%>%
    filter(Sample_type!="DP" & Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
    group_by(Time_point>0)%>%
    dplyr::summarise(new_ID=unique(new_ID),total_samples=n(),n_with_drivers=sum(n_drivers>0),n_without_drivers=sum(n_drivers==0),prop_drivers=sum(n_drivers>0)/n())
  prop_drivers_df$min.prop=sapply(1:nrow(prop_drivers_df),function(i) binom.test(prop_drivers_df$n_with_drivers[i],prop_drivers_df$total_samples[i])$conf.int[1])
  prop_drivers_df$max.prop=sapply(1:nrow(prop_drivers_df),function(i) binom.test(prop_drivers_df$n_with_drivers[i],prop_drivers_df$total_samples[i])$conf.int[2])
  prop_drivers_df
})%>%bind_rows()%>%
  mutate(Pre_or_post=factor(ifelse(`Time_point > 0`,"Post-GT","Pre-GT"),levels=c("Pre-GT","Post-GT")))%>%
  ggplot(aes(x=Pre_or_post,col=Pre_or_post,y=prop_drivers,ymin=min.prop,ymax=max.prop))+
  geom_point(size=1)+
  geom_errorbar(linetype=1,width=0.5,alpha=0.3)+
  facet_grid(cols=vars(new_ID))+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.position = "none")+
  labs(x="Time point",y="Proportion of colonies\nwith a driver mutation",col="Sample type")

ggsave(filename = paste0(plots_dir,"driver_proportion_by_individual_plot.pdf"),driver_proportion_by_individual_plot,width=3,height=2)

prop_drivers_df<-dplyr::bind_rows(all.res_df)%>%
  left_join(sample_metadata%>%dplyr::select(Sample,Time_point,new_ID),by="Sample")%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  group_by(Time_point>0)%>%
  dplyr::summarise(total_samples=n(),n_with_drivers=sum(n_drivers>0),n_without_drivers=sum(n_drivers==0),prop_drivers=sum(n_drivers>0)/n())
prop_drivers_df$min.prop=sapply(1:nrow(prop_drivers_df),function(i) binom.test(prop_drivers_df$n_with_drivers[i],prop_drivers_df$total_samples[i])$conf.int[1])
prop_drivers_df$max.prop=sapply(1:nrow(prop_drivers_df),function(i) binom.test(prop_drivers_df$n_with_drivers[i],prop_drivers_df$total_samples[i])$conf.int[2])

driver_proportion_combined_plot<-prop_drivers_df%>%
  mutate(Pre_or_post=factor(ifelse(`Time_point > 0`,"Post-GT","Pre-GT"),levels=c("Pre-GT","Post-GT")))%>%
  ggplot(aes(x=Pre_or_post,col=Pre_or_post,y=prop_drivers,ymin=min.prop,ymax=max.prop))+
  geom_point(size=1)+
  geom_errorbar(linetype=1,width=0.2,alpha=0.3)+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),legend.position = "none")+
  labs(x="Time point",y="Proportion of colonies\nwith a driver mutation",col="Sample type")

ggsave(filename = paste0(plots_dir,"driver_proportion_combined_plot.pdf"),driver_proportion_combined_plot,width=2,height=2)


bind_rows(all.res_df)%>%
  filter(!is.na(driver_ids))%>%
  left_join(sample_metadata)

#See if proportion of drivers in post-GT samples is higher than expected by chance
#Fisher's exact test
#For SCD patients - look at total numbers of different sample types
sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  group_by(ID,Sample_type)%>%
  dplyr::summarise(n=n())

sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  group_by(Time_point>0)%>%
  dplyr::summarise(n=n())

#Pool data for SCD3 and SCD4 (and driver numbers info) to make contingency table
#Shows more driver mutations in the post-GT samples than would be expected by chance
fisher.test(matrix(c(1160,1431,1,11),nrow=2))


#######################################################################
##---------------------------DNDSCV ANALYSIS---------------------------
#######################################################################

library(dndscv)
dndscv_on_details=function(details,id="this_sample",outp=1,max_muts_per_gene_per_sample = Inf,use_indel_sites=T,max_coding_muts_per_sample = Inf,...) {
  require(dndscv)
  if(!"sampleID"%in%colnames(details)) {
    details$sampleID=id
  }
  muts=details[,c("sampleID","Chrom","Pos","Ref","Alt")]
  colnames(muts)<-c("sampleID","chr","pos","ref","alt")
  
  dndscvout=dndscv(muts,outp=outp,max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,use_indel_sites=use_indel_sites,min_indels=2,max_coding_muts_per_sample = max_coding_muts_per_sample,...)
  return(dndscvout)
}
mutation_type_vec=c("All","Missense","Nonsense","Splice site","Truncating")
names(mutation_type_vec)=c("wall","wmis","wnon","wspl","wtru")

all.dndscv=Map(details=all.muts.nodups,exp_ID=names(all.muts.nodups),function(details,exp_ID) {cat(exp_ID,sep="\n");dndscv_on_details(details%>%mutate(sampleID=exp_ID))})
Individual_dNdS_plot<-Map(outp=all.dndscv,exp_ID=names(all.dndscv),function(outp,exp_ID) outp$globaldnds%>%mutate(exp_ID=exp_ID))%>%
  bind_rows()%>%
  left_join(Individual_metadata,by=c("exp_ID"="ID"))%>%
  filter(name%in%c("wall","wmis"))%>%
  mutate(name=mutation_type_vec[name])%>%
  ggplot(aes(x=new_ID,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point()+
  geom_errorbar(width=0.2)+
  theme_bw()+
  facet_grid(cols=vars(name))+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  labs(x="Indvidual",y="dN/dS ratio")+
  my_theme

ggsave(filename = paste0(plots_dir,"Individual_dNdS_plot.pdf"),Individual_dNdS_plot,width=3,height=2)

combined.dndscv=dndscv_on_details(dplyr::bind_rows(Map(details=all.muts.nodups,exp_ID=names(all.muts.nodups),function(details,exp_ID) details%>%mutate(sampleID=exp_ID))),outp=3)
combined_dNdS_plot<-combined.dndscv$globaldnds%>%
  mutate(name=mutation_type_vec[name])%>%
  ggplot(aes(x=name,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point()+
  geom_errorbar(width=0.2)+
  theme_bw()+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  labs(y="Estimated dN/dS ratio",x="Mutation type")+
  my_theme

ggsave(filename = paste0(plots_dir,"combined_dNdS_plot.pdf"),combined_dNdS_plot,width=3,height=2)

#Define function to extract the details matrix for just post or pre-GT samples
get_pre_or_postGT_details=function(tree,details,type="pre",sample_metadata){
  selected_nodes=tree$edge[,2][sapply(tree$edge[,2],function(node) {
    samples=getTips(tree,node)
    if(type=="pre"){
      return(any(sample_metadata%>%filter(Sample%in%samples)%>%pull(Time_point)==0))
    } else if(type=="post"){
      return(any(sample_metadata%>%filter(Sample%in%samples)%>%pull(Time_point)>0))
    }
  })]
  details_filtered<-details%>%filter(node%in%selected_nodes)
  return(details_filtered)
}

#Compare combined pre-GT and post-GT trees
all.muts.nodups.postGT=Map(tree=all.trees.cc.nodups,details=all.muts.nodups,function(tree,details){
  details_post=get_pre_or_postGT_details(tree=tree,details=details,type="post",sample_metadata=sample_metadata)
  return(details_post)
})

all.muts.nodups.preGT=Map(tree=all.trees.cc.nodups,details=all.muts.nodups,function(tree,details){
  details_pre=get_pre_or_postGT_details(tree=tree,details=details,type="pre",sample_metadata=sample_metadata)
  return(details_pre)
})

preGT.dndscv=dndscv_on_details(dplyr::bind_rows(all.muts.nodups.preGT),outp = 3)
head(preGT.dndscv$sel_cv)
postGT.dndscv=dndscv_on_details(dplyr::bind_rows(all.muts.nodups.postGT),outp = 3)
head(postGT.dndscv$sel_cv)

pre_vs_post_dNdS_plot<-preGT.dndscv$globaldnds%>%
  mutate(type="pre")%>%
  bind_rows(postGT.dndscv$globaldnds%>%mutate(type="post"))%>%
  mutate(type=factor(type,levels=c("pre","post")))%>%
  filter(name=="wall"|name=="wmis")%>%
  mutate(name=mutation_type_vec[name])%>%
  ggplot(aes(x=type,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point()+
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  theme_bw()+
  facet_grid(cols=vars(name))+
  labs(x="Sample timing",y="dN/dS ratio")+
  my_theme

ggsave(filename = paste0(plots_dir,"pre_vs_post_dNdS_plot.pdf"),pre_vs_post_dNdS_plot,width=3,height=2)

pre_vs_post_dNdS_plot_all_types<-preGT.dndscv$globaldnds%>%
  mutate(type="pre")%>%
  bind_rows(postGT.dndscv$globaldnds%>%mutate(type="post"))%>%
  mutate(type=factor(type,levels=c("pre","post")))%>%
  #filter(name=="wall"|name=="wmis")%>%
  mutate(name=mutation_type_vec[name])%>%
  ggplot(aes(x=type,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point(size=0.5)+
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  theme_bw()+
  facet_grid(cols=vars(name))+
  labs(x="Sample timing",y="dN/dS ratio")+
  my_theme

ggsave(filename = paste0(plots_dir,"pre_vs_post_dNdS_plot_all_types.pdf"),pre_vs_post_dNdS_plot_all_types,width=5,height=2)

#Run pre & post for each individual separately - save these as they are very memory intensive
if(file.exists(paste0(root_dir,"/Data/all.dndscv.pre.Rds"))){
  cat("Importing previously run dndscv output",sep="\n")
  all.dndscv.pre<-readRDS(paste0(root_dir,"/Data/all.dndscv.pre.Rds"))
} else {
  all.dndscv.pre=Map(details=all.muts.nodups.preGT,exp_ID=names(all.muts.nodups),function(details,exp_ID) {cat(exp_ID,sep="\n");dndscv_on_details(details%>%mutate(sampleID=exp_ID))})
  saveRDS(all.dndscv.pre,file = paste0(root_dir,"/Data/all.dndscv.pre.Rds"))
}

if(file.exists(paste0(root_dir,"/Data/all.dndscv.post.Rds"))){
  cat("Importing previously run dndscv output",sep="\n")
  all.dndscv.post<-readRDS(paste0(root_dir,"/Data/all.dndscv.post.Rds"))
} else {
  all.dndscv.post=Map(details=all.muts.nodups.postGT,exp_ID=names(all.muts.nodups),function(details,exp_ID) {cat(exp_ID,sep="\n");dndscv_on_details(details%>%mutate(sampleID=exp_ID))})
  saveRDS(all.dndscv.post,file = paste0(root_dir,"/Data/all.dndscv.post.Rds"))
}

pre_vs_post_dNdS_plot_all_types_by_individual<-bind_rows(Map(outp=all.dndscv.pre,exp_ID=names(all.dndscv.pre),function(outp,exp_ID) outp$globaldnds%>%mutate(exp_ID=exp_ID,timing="pre")))%>%
  bind_rows(bind_rows(Map(outp=all.dndscv.post,exp_ID=names(all.dndscv.post),function(outp,exp_ID) outp$globaldnds%>%mutate(exp_ID=exp_ID,timing="post"))))%>%
  mutate(timing=factor(timing,levels=c("pre","post")))%>%
  #filter(name=="wall"|name=="wmis")%>%
  mutate(name=mutation_type_vec[name])%>%
  ggplot(aes(x=timing,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point(size=0.5)+
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  theme_bw()+
  facet_grid(cols=vars(name),rows=vars(exp_ID))+
  labs(x="Sample timing",y="dN/dS ratio")+
  my_theme

ggsave(filename = paste0(plots_dir,"pre_vs_post_dNdS_plot_all_types_by_individual.pdf"),pre_vs_post_dNdS_plot_all_types_by_individual,width=5,height=5)

#Single plot of missense dN/dS
missense_dNdS_combined_categories_plot<-Map(outp=all.dndscv,exp_ID=names(all.dndscv),function(outp,exp_ID) outp$globaldnds%>%mutate(exp_ID=exp_ID))%>%
  bind_rows()%>%
  left_join(Individual_metadata,by=c("exp_ID"="ID"))%>%
  dplyr::rename("ID"="new_ID")%>%
  filter(name%in%c("wall","wmis"))%>%
  mutate(name=mutation_type_vec[name],type="By individual")%>%
  bind_rows(combined.dndscv$globaldnds%>%mutate(name=mutation_type_vec[name],ID="All",type="Combined"))%>%
  bind_rows(preGT.dndscv$globaldnds%>%mutate(ID="Pre")%>%bind_rows(postGT.dndscv$globaldnds%>%mutate(ID="Post"))%>%mutate(name=mutation_type_vec[name],type="By time point"))%>%
  dplyr::filter(name=="Missense")%>%
  dplyr::mutate(type=factor(type,levels=c("Combined","By individual","By time point")))%>%
  dplyr::mutate(ID=factor(ID,levels=c("All",paste0("SCD",1:4),"Pre","Post")))%>%
  ggplot(aes(x=ID,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point()+
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 1,col="red",linetype=2)+
  theme_bw()+
  facet_grid(cols=vars(type),scales = "free",space = "free")+
  labs(x="Sample set",y="dN/dS ratio (Missense mutations)")+
  my_theme

ggsave(filename = paste0(plots_dir,"Missense_dNdS_combined_categories.pdf"),missense_dNdS_combined_categories_plot,width=4,height=2)

#######################################################################
##---------------------------VECTOR COPY NUMBER ANALYSIS---------------
#######################################################################

##20. Compare raw VCN distribution histograms for all sample types/ individuals
sample_metadata%>%
  filter(!is.na(VCN))%>%
  ggplot(aes(x=VCN))+
  geom_histogram(binwidth = 0.2)+
  facet_grid(Sample_type~new_ID,scales="free")+
  theme_bw()+
  labs(x="Raw Vector Copy Number",y="Count")

#Review the rounded VCN by time point for samples that contain transduced cells
VCN_histograms<-sample_metadata%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  mutate(Sample_type=factor(Sample_type,levels=c("DP","Post-GT")))%>%
  mutate(Time_point=paste("Time point =",Time_point))%>%
  ggplot(aes(x=VCN_rounded,fill=Sample_type))+
  geom_bar(stat="count",width=0.7,col=NA)+
  scale_x_continuous(breaks=seq(0,12,1))+
  facet_wrap(new_ID~Time_point,scales="free_y",nrow = 3)+
  theme_bw()+
  my_theme+
  theme(strip.text.x = element_text(size=6),strip.text = element_text(size=6))+
  labs(x="Vector copy number",y="Count",fill="Sample type")

ggsave(filename=paste0(plots_dir,"VCN_histos_plot.pdf"),VCN_histograms,width=7,height=3.5)

vector_containing_cell_prop_plot<-sample_metadata%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  mutate(Sample_type=factor(Sample_type,levels=c("DP","Post-GT")))%>%
  group_by(new_ID,Time_point)%>%
  dplyr::summarise(n_untransduced=sum(VCN_rounded==0),n_VCN1=sum(VCN_rounded==1),n_VCN_over1=sum(VCN_rounded>1))%>%
  mutate(untransduced=n_untransduced/(n_untransduced+n_VCN1+n_VCN_over1),
         `VCN=1`=n_VCN1/(n_untransduced+n_VCN1+n_VCN_over1),
         `VCN >1`=n_VCN_over1/(n_untransduced+n_VCN1+n_VCN_over1))%>%
  dplyr::select(-n_untransduced,-n_VCN1,-n_VCN_over1)%>%
  mutate(Transduced=1-untransduced)%>%
  gather(-new_ID,-Time_point,key="VCN",value="Proportion")%>%
  filter(VCN=="Transduced")%>%
  left_join(Individual_metadata,by="new_ID")%>%
  ggplot(aes(x=Time_point,y=Proportion,col=new_ID))+
  geom_point(size=0.8,alpha=0.6)+
  geom_line(alpha=0.5)+
  scale_color_manual(values=patient_cols)+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_x_continuous(limits=c(0,3.2),breaks=seq(0,3,1))+
  theme_classic()+
  labs(x="Time point (years)",y="Proportion of colonies\ncontaining vector",col="Individual")+
  my_theme+
  theme(legend.box.spacing = unit(0,"mm"),legend.key.size=unit(2.5,"mm"))

ggsave(filename=paste0(plots_dir,"vector_containing_cell_prop_plot.pdf"),vector_containing_cell_prop_plot,width=2.5,height=2)

#Compare raw vector copy number (by coverage) and the VCN by number of vector integration sites identified.
sample_metadata%>%
  filter(!is.na(VCN) & Sample_type%in%c("DP","Post-GT"))%>%
  ggplot(aes(x=VCN,y=VCN_by_VIS,label=Sample,col=Sample_type))+
  geom_abline(intercept = 0,slope=1)+
  #ggrepel::geom_label_repel()+
  geom_jitter(height=0.15,width=0,alpha=0.3) +
  scale_y_continuous(breaks = seq(0,10,1),limits=c(0,10))+
  scale_x_continuous(breaks = seq(0,10,1),limits=c(0,10))+
  facet_grid(~new_ID)+
  theme_bw()+
  labs(x="VCN by coverage",y="VCN by number of\nvector integration sites",col="Sample type")


##################Other###########################

all_samples_burden_plot<-sample_metadata%>%
  dplyr::filter(!is.na(SNV_burden_AR))%>%
  dplyr::filter(Coverage>10)%>%
  mutate(Age_at_sampling=(Time_point+Age_at_GT))%>%
  ggplot(aes(x=Age_at_sampling,y=SNV_burden_AR,col=new_ID))+
  geom_jitter(height=0,width=0.15,alpha=0.3,size=0.5)+
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

##14. Review the relationship between mutation burden and peak vaf (i.e. sample clonality)
sample_metadata%>%
  filter(!is.na(peak_vaf))%>%
  filter(Coverage>8)%>%
  ggplot(aes(x=peak_vaf,y=SNV_burden_tc,col=ID))+
  geom_point(size=0.1,alpha=0.3)+
  scale_x_continuous(limits=c(0.3,0.7))+
  scale_y_continuous(limits=c(50,700))+
  facet_grid(Time_point~ID)+
  geom_smooth(method="lm",col="black",linewidth=0.5)+
  theme_bw()+
  my_theme
