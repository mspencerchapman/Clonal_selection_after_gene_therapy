#Engrafting clones inference

#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","stringr","ggridges","tidyr","RColorBrewer","ape")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
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

#----------------------------------
# Set file paths and import data
#----------------------------------

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Define paths for the script
root_dir<-"~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
source(paste0(root_dir,"/Data/GT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")

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

#Summary stats data.frame
data_stats=data.frame(ID=c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009"),
                      n_post_GT_clones=c(420,358,292,143,143,74),
                      no_of_coalescences=c(3,5,0,0,1,2))
#
get_expanded_clade_nodes(all.trees.cc.nodups$BCL009,height_cut_off = 300,min_clonal_fraction = 0,min_samples = 2)
sample_metadata%>%dplyr::filter(sample_status=="PASS")%>%group_by(ID,Time_point>0)%>%dplyr::summarise(n=n())%>%filter(`Time_point > 0`)

#----------------------------------
# Simulate the anticipated clone size distributions from different engrafting cell numbers
#----------------------------------

regenerate_clone_sizes=F
if(regenerate_clone_sizes) {
  
  clone_size_dir=paste0(root_dir,"/Data/Engrafting_cell_inference/Clone_sizes")
  
  # function to simply generate clone sizes from a birth process (much quicker than using RSimpop!!)
  generate_clone_sizes=function(n_engrafted,n_final){
    cat(paste("Generating clone size information with",n_engrafted,"engrafting clones, and a final population of",n_final,"\n"))
    clone_sizes=rep(1,n_engrafted) #all clones start with a population of 1
    for(i in 1:(n_final-n_engrafted)){i=sample(x=n_engrafted,size=1,prob=clone_sizes);clone_sizes[i]<-clone_sizes[i]+1} # increment the population by 1 u
    return(data.frame(clone=1:n_engrafted,size=clone_sizes))
  }
  all_n_engrafted=round(2^(seq(10,16.6,0.1)))
  all_n_final=c(1e5,2e5,5e5,1e6,2e6)
  clone_params_df=tidyr::expand_grid(all_n_engrafted,all_n_final)%>%dplyr::rename("n_engrafted"=all_n_engrafted,"final_population"=all_n_final)
  temp=lapply(1:nrow(clone_params_df),function(i) {
    cat(i,sep="\n")
    file_name=paste0(clone_size_dir,"/Clone_sizes_",clone_params_df$n_engrafted[i],"_",clone_params_df$final_population[i],".tsv")
    if(!file.exists(file_name)){
      clone_sizes=generate_clone_sizes(n_engrafted = clone_params_df$n_engrafted[i],n_final = clone_params_df$final_population[i])
      write.table(clone_sizes,file=file_name,quote=F,sep="\t",row.names=F)
    }
  })
}

#----------------------------------
# Import the clone sizes data frames
#----------------------------------

clone_size_dir=paste0(root_dir,"/Data/Engrafting_cell_inference/Clone_sizes")
setwd(clone_size_dir)
clone_files=list.files(pattern="Clone_sizes_")
clone_sizes_list<-lapply(clone_files,function(file){
  clone_sizes<-read.delim(file,stringsAsFactors = F)
  return(clone_sizes)
})

clone_params_df<-lapply(clone_files,function(file){
  n_engrafted<-stringr::str_split(file,pattern="_",simplify=T)[,3]
  final_population<-gsub(".tsv","",stringr::str_split(file,pattern="_",simplify=T)[,4])
  return(data.frame(n_engrafted=as.numeric(n_engrafted),final_population=as.numeric(final_population)))
})%>%dplyr::bind_rows()
names(clone_sizes_list)=apply(clone_params_df,1,paste,collapse="_")

#----------------------------------
# Run the simulations (using the clone sizes) & get posteriors
#----------------------------------

sim_res_dir=paste0(root_dir,"/Data/Engrafting_cell_inference/Simulation_results/")
for(j in 1:nrow(data_stats)){
  ID_to_run<-data_stats$ID[j]
  cat(ID_to_run,sep="\n")
  no_of_coalescences<-data_stats$no_of_coalescences[j]
  n_post_GT_clones<-data_stats$n_post_GT_clones[j]
  
  sim_res_file=paste0(sim_res_dir,"sim_res_",ID_to_run,".RDS")
  if(file.exists(sim_res_file)) {
    cat("Importing existing simulation results files")
    sim_res<-readRDS(sim_res_file)
    
    completed_params<-sim_res%>%
      dplyr::select(n_engrafted,final_population)%>%
      dplyr::filter(!duplicated(.))
    
    new_params_df<-anti_join(clone_params_df,completed_params)
    print(new_params_df)
    
    #Update the sim_res object with any simulation results from clone distributions that weren't previously complete
    if(nrow(new_params_df)>0){
      iter=1e3
      new_sim_res<-lapply(1:nrow(new_params_df),function(i) {
        cat(i,sep = "\n")
        clones_df<-clone_sizes_list[[paste(new_params_df[i,],collapse="_")]]
        n_coalescneces_vec<-sapply(1:iter,function(j) {
          pop_vec<-unlist(Map(clone=clones_df$clone,size=clones_df$size,f=function(clone,size) {rep(clone,size)}))
          sampled_clones<-sample(pop_vec,size=n_post_GT_clones,replace=F)
          #sampled_clones<-sample(clones_df$clone,size = n_post_GT_clones,replace = T,prob=clones_df$size)
          n_coalescences=length(sampled_clones)-length(unique(sampled_clones))
          return(n_coalescences)
        })
        return(data.frame(n_engrafted=new_params_df$n_engrafted[i],final_population=new_params_df$final_population[i],n_coalescences=n_coalescneces_vec))
      })%>%dplyr::bind_rows()
      
      sim_res<-bind_rows(sim_res,new_sim_res)
      saveRDS(sim_res,file=sim_res_file)
    }
    
  } else {
    cat("Running simulations")
    iter=1e3
    sim_res<-lapply(1:nrow(clone_params_df),function(i) {
      cat(i)
      clones_df<-clone_sizes_list[[paste(clone_params_df[i,],collapse="_")]]
      n_coalescneces_vec<-sapply(1:iter,function(j) {
        pop_vec<-unlist(Map(clone=clones_df$clone,size=clones_df$size,f=function(clone,size) {rep(clone,size)}))
        sampled_clones<-sample(pop_vec,size=n_post_GT_clones,replace=F)
        #sampled_clones<-sample(clones_df$clone,size = n_post_GT_clones,replace = T,prob=clones_df$size)
        n_coalescences=length(sampled_clones)-length(unique(sampled_clones))
        return(n_coalescences)
      })
      return(data.frame(n_engrafted=clone_params_df$n_engrafted[i],final_population=clone_params_df$final_population[i],n_coalescences=n_coalescneces_vec))
    })%>%dplyr::bind_rows()
    
    saveRDS(sim_res,file=sim_res_file)
    
  }
}


#----------------------------------
# Plot the results
#----------------------------------

res_df<-Map(ID_to_run=data_stats$ID,no_of_coalescences=data_stats$no_of_coalescences,f=function(ID_to_run,no_of_coalescences) {
  sim_res_file=paste0(sim_res_dir,"sim_res_",ID_to_run,".RDS")
  sim_res<-readRDS(sim_res_file)
  
  total_sims<-sim_res%>%
    #filter(final_population==5e5)%>%
    group_by(n_engrafted)%>%
    dplyr::summarise(total_sims_n=n())
  
  sim_summary<-sim_res%>%
    #filter(final_population<1e6)%>%
    mutate(result=ifelse(n_coalescences==no_of_coalescences,"PASS","FAIL"))%>%
    dplyr::count(n_engrafted,final_population,result)%>%
    tidyr::complete(n_engrafted,final_population,result,fill=list(n=0))%>%
    dplyr::mutate(prop=n/1000,ID=ID_to_run)
  
  return(sim_summary)
})%>%bind_rows()%>%
  left_join(Individual_metadata,by="ID")

#Look at likelihood of data given engrafting cell number - divide by final population
engrafting_cell_inference_by_final_pop_plot<-res_df%>%
  filter(result=="PASS")%>%
  ggplot(aes(x=n_engrafted,y=prop,col=factor(final_population)))+
  geom_point(size=0.5)+
  geom_line()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  facet_grid(rows=vars(new_ID),scales="free_y")+
  #scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_x_log10(breaks=c(1000,2000,3000,6000,10000,20000,30000,60000,100000),labels=scales::label_comma(),limits=c(1e3,1e5))+ #round(2^(seq(10,16.7,0.5)))
  labs(x="Engrafting cell number",y="Likelihood of data given engrafting cell number",fill="Result",col="Final population size\nused in simulation")+
  my_theme

ggsave(filename=paste0(plots_dir,"engrafting_cell_inference_by_final_pop_plot.pdf"),engrafting_cell_inference_by_final_pop_plot,width=5,height=6)

#Look at likelihood of data given engrafting cell number - combine results from different final population
engrafting_cell_inference_combined_plot<-res_df%>%
  group_by(new_ID,n_engrafted,result)%>%
  dplyr::summarise(n=sum(n))%>%
  tidyr::pivot_wider(names_from="result",values_from="n")%>%
  mutate(prop=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL)%>%
  ggplot(aes(x=n_engrafted,y=prop))+
  geom_point(size=0.5)+
  geom_line()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  facet_grid(rows=vars(new_ID),scales="free_y")+
  #scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_x_log10(breaks=c(1000,2000,3000,6000,10000,20000,30000,60000,100000),labels=scales::label_comma(),limits=c(1e3,1e5))+ #round(2^(seq(10,16.7,0.5)))
  labs(x="Engrafting cell number",y="Likelihood of data given engrafting cell number",fill="Result",col="Final population size\nused in simulation")+
  my_theme

ggsave(filename="engrafting_cell_inference_combined_plot.pdf",engrafting_cell_inference_combined_plot,width=3.5,height=5)

#Now display as a posterior probability by dividing by the sum of simulations matching the data in each individual
#Calculate the sum of proportions
total_props<-res_df%>%
  group_by(new_ID,n_engrafted,result)%>%
  dplyr::summarise(n=sum(n))%>%
  tidyr::pivot_wider(names_from="result",values_from="n")%>%
  mutate(prop=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL)%>%
  #pivot_wider(id_cols=c("new_ID","n_engrafted"),names_from="n_engrafted",values_from="prop")%>%
  group_by(new_ID)%>%
  dplyr::summarise(total_prop=sum(prop))

post_prob_df<-res_df%>%
  group_by(new_ID,n_engrafted,result)%>%
  dplyr::summarise(n=sum(n))%>%
  tidyr::pivot_wider(names_from="result",values_from="n")%>%
  mutate(prop=PASS/(FAIL+PASS))%>%
  dplyr::select(-FAIL)%>%
  left_join(total_props,by="new_ID")%>%
  dplyr::mutate(prob=prop/total_prop)%>%
  group_by(new_ID)%>%
  dplyr::mutate(cum_prob=cumsum(prob))

#Divide proportions by these values to get normalized 'probabilities'
engrafting_cell_inference_posterior_plot<-post_prob_df%>%
  ggplot(aes(x=n_engrafted,y=prob))+
  geom_point(size=0.5)+
  geom_line()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  facet_grid(rows=vars(new_ID),scales="free_y")+
  #scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.5))+
  scale_x_log10(breaks=c(1000,2000,3000,6000,10000,20000,30000,60000,100000),labels=scales::label_comma(),limits=c(1e3,1e5))+ #round(2^(seq(10,16.7,0.5)))
  labs(x="Engrafting cell number",y="Posterior probability")+
  my_theme

ggsave(filename="engrafting_cell_inference_posterior_plot.pdf",engrafting_cell_inference_posterior_plot,width=3.5,height=5)

##Plot the engrafting CD34 estimates against the CD34 cell dose
all_post_samples<-lapply(paste0("SCD",1:6),function(ID) {
  post_prob_this_ID<-post_prob_df%>%filter(new_ID==ID)
  post<-sample(post_prob_this_ID$n_engrafted,size = 10000,replace = T,prob=post_prob_this_ID$prob)
  return(data.frame(new_ID=ID,post=post))
})%>%dplyr::bind_rows()
post_vs_CD34<-all_post_samples%>%
  left_join(Individual_metadata)%>%
  ggplot(aes(x=CD34_dose,y=post,col=new_ID))+
  geom_jitter(height=0.01,width=0.05,size=0.05,alpha=0.01)+
  scale_y_log10(limits=c(1e2,1.1e5))+
  scale_x_continuous(limits=c(1,10))+
  scale_color_manual(values=patient_cols)+
  theme_bw()+
  my_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  labs(x="CD34+ cell dose per kg",y="Posterior",col="ID")
ggsave(filename=paste0(plots_dir,"engrafting_cell_no_vs_CD34dose_plot.pdf"),post_vs_CD34,width=3.5,height=2)

#Do LME regression of engrafting CD34 number against CD34 cell dose
lme.CD34<-lme4::lmer(post~CD34_dose+(1|new_ID),data=all_post_samples%>%
                       left_join(Individual_metadata))
summary(lme.CD34)
confint(lme.CD34)

#----------------------------------
# Compare results from our data to those obtained using vector integration site analysis
#----------------------------------

##Read in the estimates done using VIS analysis
VCN_estimates<-readxl::read_excel(paste0(root_dir,"/Data/SCDEstimatedTotalHSCpop_tidy.xlsx"))
VCN_estimates<-VCN_estimates%>%left_join(Individual_metadata,by=c("pt"="ID"))%>%
  dplyr::mutate(lower_CI=estimate-1.96*se,upper_CI=estimate+1.96*se)%>%
  dplyr::select(new_ID,stat,estimate,lower_CI,upper_CI)

#Compare estimates using different approaches
phylo_estimates_vs_VIS_estimates<-all_post_samples%>%
  group_by(new_ID)%>%
  dplyr::summarise(stat="Phylo",estimate=median(post),lower_CI=quantile(post,0.025),upper_CI=quantile(post,0.975))%>%
  bind_rows(VCN_estimates)%>%
  mutate(stat=factor(stat,levels=c("Phylo","chao","jack1","jack2","Bootstrap")))%>%
  ggplot(aes(x=stat,y=estimate,ymin=lower_CI,ymax=upper_CI,col=stat))+
  geom_point(size=0.5)+
  geom_errorbar(width=0.2)+
  facet_grid(~new_ID)+
  theme_bw()+
  my_theme+
  scale_y_log10(breaks=c(1000,3000,10000,30000,100000,300000),labels=scales::label_comma())+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Method",y="Engrafting HSC number",col="Method")
ggsave(filename=paste0(plots_dir,"phylo_estimates_vs_VIS_estimates.pdf"),phylo_estimates_vs_VIS_estimates,width=6,height=2)

#Visualize posterior as density
engrafting_cell_inference_posterior_density_plot<-all_post_samples%>%
  ggplot(aes(x=post,y=new_ID,fill=new_ID))+
  ggridges::geom_density_ridges2(scale=4)+
  scale_x_log10(breaks=c(1e3,3e3,1e4,3e4,1e5),labels=scales::label_comma())+
  theme_classic()+
  my_theme+
  scale_fill_manual(values=patient_cols)+
  labs(x="Number of engrafting HSCs",y="ID",fill="ID")

ggsave(filename=paste0(plots_dir,"engrafting_cell_inference_posterior_density_plot.pdf"),engrafting_cell_inference_posterior_density_plot,width=3.5,height=2.5)

##Other plots
sim_coalescences_plot<-sim_res%>%
  filter(final_population<1e6)%>%
  mutate(result=ifelse(n_coalescences==no_of_coalescences,"PASS","FAIL"))%>%
  ggplot(aes(x=n_coalescences,fill=n_coalescences==no_of_coalescences))+
  geom_histogram(binwidth=1)+
  facet_grid(~n_engrafted)+
  scale_fill_manual(values=c("#FDBF6F","#1F78B4"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="No of coalescences in simulation",y="Count",fill="Matches data")
ggsave(filename = paste0(ID_to_run,"_sim_coalescences.pdf"),plot=sim_coalescences_plot,height=3,width=7)

# sim_res%>%
#   filter(final_population<1e6)%>%
#   mutate(result=ifelse(n_coalescences==no_of_coalescences,"PASS","FAIL"))%>%
#   group_by(n_engrafted)%>%
#   summarise(n=n())

