#-----------------------------
# Sigs on Trees
# Tim Coorens - November 2019
#-----------------------------
library(ape)
library(ggtree)
library(RColorBrewer)
library(hdp)
library(stringr)
options(stringsAsFactors = F)

#home_dir="/lustre/scratch119/realdata/mdt1/team154/ms56/gene_therapy"
#home_dir="~/Mounts/Lustre/gene_therapy"
#home_dir="~/Mounts/Lustre/Mut_Sig_Combined/"
#setwd(paste0(home_dir,"/HDP"))
root_dir<-"~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
setwd(paste0(root_dir,"/Data/HDP/"))

mut_example_multi=readRDS("HDP_multi_chain.Rdata")

mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")
sample_remove=rownames(mutations)[rowSums(mutations)<50]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]
freq=nrow(mutations)
#freq=table(mutations)

dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
colnames(exposures)=rownames(mutations)
patients=c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
#patients=patients[patients%in%names(freq)]
sigs=rownames(exposures)
sig_profiles=mut_example_multi@comp_categ_distn$mean

for (patient in patients){
  tree <-read.tree(paste0(root_dir,"/Data/tree_files_with_dups/tree_",patient,"_m40_postMS_reduced_a_j_vaf_post_mix.tree"))
  tree_df=fortify(tree)
    cols=c(brewer.pal(12, "Paired"),"magenta","firebrick")
    #cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
    #     "turquoise1",'orange2',"chartreuse","mediumorchid3","grey20",
    #     "yellow2","mediumaquamarine","tomato4")
    names(cols)=sigs
    samples=colnames(exposures)[grepl(paste0("_",patient),colnames(exposures))]
    branches=as.numeric(stringr::str_split(samples,pattern="_",simplify=TRUE)[,1])
    
    pdf(paste0(root_dir,"/plots/",patient,"_tree_with_hdp_signatures_branch_length_newcols.pdf"),height=20,width=15)
    plot(tree,label.offset=0.01*max(tree_df$x),show.tip.label=TRUE,cex=0.3)
    for (k in 1:length(samples)){
      n=as.numeric(branches[k])
      x_end=tree_df$x[n]
      x_start=tree_df$x[tree_df$parent[n]]
      x_intv=x_end-x_start
      y=node.height(tree)[n]
      tipnum=sum(tree_df$isTip)
      for (s in sigs){
        x_end=x_start+exposures[s,samples[k]]*x_intv
        rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s])
        x_start=x_end
      }
    }
    axisPhylo(side = 1,backward=F)
    legend("topright",title="Signatures", legend=paste0("N",sigs),
           fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
    dev.off()
    
    tree_collapsed=tree
    tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
    tree_collapsed$edge.length[tree$edge.length==0]=0
    tree_collapsed_df=fortify(tree_collapsed)
    
    pdf(paste0(root_dir,"/plots/",patient,"_tree_with_hdp_signatures_equal_branch_length_newcols.pdf"),height=20,width=20)
    plot(tree_collapsed,label.offset=0.01*max(tree_collapsed_df$x),show.tip.label=TRUE,cex=0.5)
    for (k in 1:length(samples)){
      n=as.numeric(branches[k])
      x_end=tree_collapsed_df$x[n]
      x_start=tree_collapsed_df$x[tree_collapsed_df$parent[n]]
      x_intv=x_end-x_start
      y=node.height(tree_collapsed)[n]
      tipnum=sum(tree_collapsed_df$isTip)
      
      for (s in sigs){
        x_end=x_start+exposures[s,samples[k]]*x_intv
        rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s])
        x_start=x_end
      }
    }
    axisPhylo(side = 1,backward=F)
    legend("topright",title="Signatures", legend=paste0("N",sigs),
           fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
    dev.off()
}

Emily_patients=c("KX004","KX008")
for (patient in Emily_patients){
  tree_files=list.files("~/R_work/Emily_benchmarking/Emily_data/",pattern="tree_",full.names = T)
  tree <-read.tree(grep(patient,tree_files,value=T))
  tree_df=fortify(tree)
  cols=c(brewer.pal(12, "Paired"),"magenta","firebrick")
  #cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
  #     "turquoise1",'orange2',"chartreuse","mediumorchid3","grey20",
  #     "yellow2","mediumaquamarine","tomato4")
  names(cols)=sigs
  samples=colnames(exposures)[grepl(paste0("_",patient),colnames(exposures))]
  branches=as.numeric(stringr::str_split(samples,pattern="_",simplify=TRUE)[,1])
  
  pdf(paste0("~/R_work/Gene_therapy/Mutational_signature_extraction/HDP_4samps/",patient,"_tree_with_hdp_signatures_branch_length_newcols.pdf"),height=20,width=15)
  plot(tree,label.offset=0.01*max(tree_df$x),show.tip.label=TRUE,cex=0.3)
  for (k in 1:length(samples)){
    n=as.numeric(branches[k])
    x_end=tree_df$x[n]
    x_start=tree_df$x[tree_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree)[n]
    tipnum=sum(tree_df$isTip)
    for (s in sigs){
      x_end=x_start+exposures[s,samples[k]]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s])
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=paste0("N",sigs),
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
  
  tree_collapsed=tree
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  tree_collapsed$edge.length[tree$edge.length==0]=0
  tree_collapsed_df=fortify(tree_collapsed)
  
  pdf(paste0("~/R_work/Gene_therapy/Mutational_signature_extraction/HDP_4samps/",patient,"_tree_with_hdp_signatures_equal_branch_length_newcols.pdf"),height=20,width=20)
  plot(tree_collapsed,label.offset=0.01*max(tree_collapsed_df$x),show.tip.label=TRUE,cex=0.5)
  for (k in 1:length(samples)){
    n=as.numeric(branches[k])
    x_end=tree_collapsed_df$x[n]
    x_start=tree_collapsed_df$x[tree_collapsed_df$parent[n]]
    x_intv=x_end-x_start
    y=node.height(tree_collapsed)[n]
    tipnum=sum(tree_collapsed_df$isTip)
    
    for (s in sigs){
      x_end=x_start+exposures[s,samples[k]]*x_intv
      rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s])
      x_start=x_end
    }
  }
  axisPhylo(side = 1,backward=F)
  legend("topright",title="Signatures", legend=paste0("N",sigs),
         fill=cols, bty="n",cex=0.8, ncol=1, xjust=0.5)
  dev.off()
}

#Plot exposures by component & patient
as.data.frame(exposures)%>%
  tibble::rownames_to_column(var = "hdp_component")%>%
  gather(-hdp_component,key="branch",value="exposure")%>%
  mutate(patient=stringr::str_split(branch,pattern="_",simplify=T)[,2])%>%
  mutate(hdp_component=paste0("hdp_",hdp_component))%>%
  ggplot(aes(x=patient,y=exposure,col=hdp_component))+
  #geom_boxplot()+
  geom_jitter(height=0,width=0.1,alpha=0.2)+
  theme_classic()+
  facet_grid(~hdp_component)+
  theme(axis.text.x = element_text(angle = 90))

as.data.frame(exposures)%>%
  tibble::rownames_to_column(var = "hdp_component")%>%
  gather(-hdp_component,key="branch",value="exposure")%>%
  mutate(patient=stringr::str_split(branch,pattern="_",simplify=T)[,2])%>%
  mutate(hdp_component=paste0("hdp_",hdp_component))%>%
  group_by(patient,hdp_component)%>%
  summarise(median=median(exposure),mean=mean(exposure))

#Compare exposures with mutation burdens
patient="BCL002"
tree <-read.tree(paste0(home_dir,"/filtering_runs/tree_files/tree_",patient,"_m40_postMS_reduced_pval_post_mix.tree"))

as.data.frame(exposures)%>%
  tibble::rownames_to_column(var = "hdp_component")%>%
  gather(-hdp_component,key="branch",value="exposure")%>%
  mutate(patient=str_split(branch,pattern="_",simplify=T)[,2])%>%
  mutate(hdp_component=paste0("hdp_",hdp_component))%>%
  mutate(branch_no=as.numeric(str_split(branch,pattern="_",simplify=T)[,1]))%>%
  filter(patient=="BCL002")%>%
  mutate(mut_burden=sapply(branch_no,function(node) {
    if(!node%in%1:length(tree$tip.label)) {
      nodes<-get_all_node_children(node,tree)
      nodes<-nodes[nodes%in%1:length(tree$tip.label)]
    } else {
      nodes<-node
    }
    node_heights=sapply(nodes,function(node) {return(nodeheight(tree,node))})
    return(mean(node_heights))
  }))%>%
  mutate(timing=sapply(branch_no,function(node){
    if(!node%in%1:length(tree$tip.label)) {
      nodes<-get_all_node_children(node,tree)
      nodes<-nodes[nodes%in%1:length(tree$tip.label)]
    } else {
      nodes<-node
    }
    samples=sapply(nodes,function(node) {return(tree$tip.label[node])})
    Timing_tags=sapply(samples,function(sample) smry_seq$Timing_tag[smry_seq$Sample==sample])
    timings=str_split(Timing_tags,pattern = "_",simplify = T)[,1]
    if(length(unique(timings))==1){
      return(timings[1])
    }else{
      return(NA)
    }
  }))%>%
  ggplot(aes(x=mut_burden,y=exposure,col=timing))+
  #geom_point()+
  geom_boxplot()+
  facet_grid(~hdp_component)#+
  #geom_smooth(method="lm")
  

