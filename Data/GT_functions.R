get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}


get_confident_duplicates=function(tree,sample_metadata,private_mut_threshold=NA,height_cut_off=NA){
  #Default mode is the 'look back from tips' threshold - if 'private_mut_threshold' is set, this is used
  if(!is.na(private_mut_threshold)){
    duplicate_clade_samples=get_duplicate_sets(tree,mut_threshold = private_mut_threshold)
  } else if(!is.na(height_cut_off)){ #else looks for the 'height_cut_off' - and assumes
    duplicate_clades<-get_expanded_clade_nodes(tree,min_samples = 2,height_cut_off=height_cut_off,min_clonal_fraction = 0)
    duplicate_clade_samples<-lapply(duplicate_clades$nodes,function(node) getTips(tree,node=node))
  }
  if(length(duplicate_clade_samples)==0) {stop(print("No duplicates detected"))}
  
  #Define the order that wells are filled
  well_numbers<-c("01","02","03","04","05","06","07","08","09","10","11","12")
  well_filling_order=paste0(rep(LETTERS[1:8],each=12),rep(c(well_numbers),times=8))
  
  
  #Get a list of the duplicate sets including plate IDs/ well positions etc.
  dups_list<-Map(duplicate_clades=duplicate_clade_samples,dup_set=1:length(duplicate_clade_samples),f=function(duplicate_clades,dup_set){
    dup_set_metadata<-sample_metadata%>%
      dplyr::filter(Sample%in%duplicate_clades)%>%
      mutate(dup_set=dup_set)%>%
      dplyr::select(Sample,ID,new_ID,dup_set,plate_ID,position,HSPC_source,Time_point)
    return(dup_set_metadata)
  })
  
  #Test for whether these are "Confident duplicates", "Not all duplicates" or "Unclear status"
  dup_status<-lapply(dups_list,function(df) {
    if(length(unique(df$HSPC_source[df$HSPC_source!="Unknown"]))>1|length(unique(df$Time_point))>1) {
      return("Not all duplicates")
    } else {
      filling_positions<-which(well_filling_order%in%df$position)
      if(max(filling_positions)-min(filling_positions)<=7){
        return("Confident duplicates")
      } else {
        return("Unclear status") #Same plate but not nearby wells
      }
    }
  })
  
  confident_dups_df<-Map(df=dups_list,status=dup_status,f=function(df,status) {
    if(status=="Confident duplicates"){
      return(df)
    } else {
      return(NULL)
    }
  })%>%dplyr::bind_rows()
  
  return(confident_dups_df)
}

get_duplicate_status=function(tree,sample_metadata,res_type="df",private_mut_threshold=NA,height_cut_off=NA){
  #Default mode is the 'look back from tips' threshold - if 'private_mut_threshold' is set, this is used
  if(!is.na(private_mut_threshold)){
    duplicate_clade_samples=get_duplicate_sets(tree,mut_threshold = private_mut_threshold)
  } else if(!is.na(height_cut_off)){ #else looks for the 'height_cut_off' - and assumes
    duplicate_clades<-get_expanded_clade_nodes(tree,min_samples = 2,height_cut_off=height_cut_off,min_clonal_fraction = 0)
    duplicate_clade_samples<-lapply(duplicate_clades$nodes,function(node) getTips(tree,node=node))
  }
  if(length(duplicate_clade_samples)==0) {stop(print("No duplicates detected"))}
  
  #Define the order that wells are filled
  well_numbers<-c("01","02","03","04","05","06","07","08","09","10","11","12")
  well_filling_order=paste0(rep(LETTERS[1:8],each=12),rep(c(well_numbers),times=8))
  
  
  #Get a list of the duplicate sets including plate IDs/ well positions etc.
  dups_list<-Map(duplicate_clades=duplicate_clade_samples,dup_set=1:length(duplicate_clade_samples),f=function(duplicate_clades,dup_set){
    dup_set_metadata<-sample_metadata%>%
      dplyr::filter(Sample%in%duplicate_clades)%>%
      mutate(dup_set=dup_set)%>%
      dplyr::select(Sample,ID,new_ID,dup_set,plate_ID,position,HSPC_source,Time_point)
    return(dup_set_metadata)
  })
  
  #Test for whether these are "Confident duplicates", "Not all duplicates" or "Unclear status"
  dups_list<-lapply(dups_list,function(df) {
    if(length(unique(df$HSPC_source[df$HSPC_source!="Unknown"]))>1|length(unique(df$Time_point))>1) {
      status<-"Not all duplicates"
    } else {
      filling_positions<-which(well_filling_order%in%df$position)
      if(max(filling_positions)-min(filling_positions)<=7){
        status<-"Confident duplicates"
      } else {
        status<-"Unclear status" #Same plate but not nearby wells
      }
    }
    df$status<-status
    return(df)
  })
  
  if(res_type=="df") {
    return(dplyr::bind_rows(dups_list))
  } else if(res_type=="list"){
    return(dups_list)
  } else {
    warning("res_type argument must be either 'df' or 'list'")
  }
  
}


#get_duplicate_status(tree = all.trees$BCL008,sample_metadata = sample_metadata,private_mut_threshold =150)


mutlist_to_96_contexts=function(mutlist,genomeFile=lustre_ref){
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  samples=unique(mutlist$SampleID)
  trinuc_mut_mat=matrix(0,ncol=96,nrow=length(samples))
  for (n in 1:length(samples)){
    s=samples[n]
    mutations=as.data.frame(mutlist[mutlist$SampleID==s,c("Chr","Pos","Ref","Alt")])
    colnames(mutations) = c("chr","pos","ref","mut")
    mutations$pos=as.numeric(mutations$pos)
    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                       as.numeric(mutations$pos)+1))))
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    trinuc_mut_mat[n,]=freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat)=full_vec
  rownames(trinuc_mut_mat)=samples
  return(trinuc_mut_mat)
}

get_mut_vafs=function(SampleID,COMB_mats,tree) {
  sample_nodes=get_ancestral_nodes(which(tree$tip.label==SampleID),edge = tree$edge)
  mutations=COMB_mats$mat$mut_ref[COMB_mats$mat$node%in%sample_nodes]
  mutation_NR<-COMB_mats$NR[mutations,SampleID]
  mutation_NV<-COMB_mats$NV[mutations,SampleID]
  mutation_vaf<-calculate_vaf(mutation_NV,mutation_NR)
  return(data.frame(SampleID=SampleID,mut_ref=mutations,NV=mutation_NV,NR=mutation_NR,vaf=mutation_vaf))
}

remove_samples_from_tree_and_update_details=function(remove_samples,tree,details){
  if(!any(remove_samples%in%tree$tip.label)){
    cat(paste("None of the specified samples (",paste0(remove_samples,collapse=","),") are in the tree. Returning original tree and details objects."),sep = "\n")
    stop(return(list(details=details,tree=tree)))
  } else {
    remove_samples<-remove_samples[which(remove_samples%in%tree$tip.label)]
    cat(paste("Sample",remove_samples,"found in tree and will be removed"),sep="\n")
  }
  
  sample_nodes=which(tree$tip.label%in%remove_samples)
  
  #Drop samples from the tree to create the new tree
  new_tree<-drop.tip(tree,tip=remove_samples)
  
  #Match up new nodes to the old nodes - create reference list, where each new node is named with the corresponding old node
  new_tree_clades=lapply(new_tree$edge[,2],function(node) getTips(new_tree,node))
  names(new_tree_clades)<-new_tree$edge[,2]
  
  old_nodes<-tree$edge[,2]
  new_nodes<-sapply(old_nodes,function(node) {
    old_clade_samples<-getTips(tree,node)
    if(all(old_clade_samples%in%remove_samples)) {
      return(NA)
    } else {
      new_clade_samples<-old_clade_samples[!old_clade_samples%in%remove_samples]
      new_idx<-which(sapply(new_tree_clades,function(clade) {setequal(clade,new_clade_samples)}))
      new_node<-new_tree$edge[,2][new_idx]
      return(new_node)
    }
  })
  names(new_nodes)<-old_nodes
  
  details$node<-new_nodes[as.character(details$node)]
  
  #Remove muts from branches that now have all tips in remove samples i.e. (removed branches)
  cat(paste(sum(is.na(details$node)),"mutations being removed from details matrix."),sep = "\n")
  details_new<-details[!is.na(details$node),]
  
  return(list(details=details_new,tree=new_tree))
}


get_duplicate_sets=function(tree,mut_threshold){
  pseudo_terminal_nodes=sapply(tree$edge[,2][!tree$edge[,2]%in%1:length(tree$tip.label)],function(node) {
    node_height=nodeheight(tree = tree,node=node)
    samples=getTips(tree=tree,node=node)
    sample_heights=nodeHeights(tree)[tree$edge[,2]%in%which(tree$tip.label%in%samples),2]
    
    if(all((sample_heights-node_height)<mut_threshold)){ #This is the method of determining the duplicates - may need to alter the "30" for low mutation burden samples
      return(node)
    }else{
      return(NA)
    }
  })
  pseudo_terminal_nodes=pseudo_terminal_nodes[!is.na(pseudo_terminal_nodes)]
  duplicate_samples=lapply(pseudo_terminal_nodes,function(node) getTips(tree,node=node))
  return(duplicate_samples)
}

plot_duplicate_plate_map=function(tree,sample_metadata,private_mut_threshold=NA,height_cut_off=NA,labels="set"){
  my_theme<-theme(text = element_text(family="Helvetica"),
                  axis.text = element_text(size = 5),
                  axis.title = element_text(size=7),
                  legend.text = element_text(size=5),
                  legend.title = element_text(size=7),
                  strip.text = element_text(size=7),
                  legend.spacing = unit(1,"mm"),
                  legend.key.size= unit(5,"mm"))
  
  #Default mode is the 'look back from tips' threshold - if 'private_mut_threshold' is set, this is used
  if(!is.na(private_mut_threshold)){
    duplicate_clade_samples=get_duplicate_sets(tree,mut_threshold = private_mut_threshold)
  } else if(!is.na(height_cut_off)){ #else looks for the 'height_cut_off' - and assumes
    duplicate_clades<-get_expanded_clade_nodes(tree,min_samples = 2,height_cut_off=height_cut_off,min_clonal_fraction = 0)
    duplicate_clade_samples<-lapply(duplicate_clades$nodes,function(node) getTips(tree,node=node))
  }
  if(length(duplicate_clade_samples)==0) {stop(print("No duplicates detected"))}
  dups_df<-Map(duplicate_clades=duplicate_clade_samples,dup_set=1:length(duplicate_clade_samples),f=function(duplicate_clades,dup_set){
    dup_set_metadata<-sample_metadata%>%
      dplyr::filter(Sample%in%duplicate_clades)%>%
      mutate(dup_set=dup_set)
    return(dup_set_metadata)
  })%>%
    dplyr::bind_rows()
  
  plot<-dups_df%>%
    mutate(column=factor(substr(position,1,1),levels=LETTERS[1:8]),
           row=factor(substr(position,2,3),levels=c("01","02","03","04","05","06","07","08","09","10","11","12")))%>%
    mutate(Sample=stringr::str_sub(Sample,start=10))%>%
    ggplot(aes(x=column,y=row))+
    geom_point(aes(col=factor(dup_set)),size=5,alpha=0.8)+
    facet_wrap(~plate_ID,drop = F)+
    scale_color_discrete(guide="none")+
    scale_x_discrete(drop=F)+
    scale_y_discrete(drop=F)+
    theme_bw()+
    my_theme+
    labs(x="Plate column",y="Plate row",col="Duplicate set")
  if(labels=="set"){
    print(plot+geom_text(aes(label=dup_set),size=2))
  } else if(labels=="sample"){
    print(plot+geom_text(aes(label=Sample),size=2))
  }
}


create_exposures_df=function(HDP_multi,trinuc_mut_mat,key_table,minimum_branch_muts=50) {
  sample_remove=rownames(trinuc_mut_mat)[rowSums(trinuc_mut_mat)<minimum_branch_muts]
  trinuc_mut_mat=trinuc_mut_mat[!rownames(trinuc_mut_mat)%in%sample_remove,]
  key_table=key_table[!key_table$Sample%in%sample_remove,]
  freq=nrow(trinuc_mut_mat)
  
  dp_distn <- comp_dp_distn(HDP_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(trinuc_mut_mat),,drop=FALSE])
  colnames(exposures)=rownames(trinuc_mut_mat)
  rownames(exposures)<-paste0("N",rownames(exposures))
  sigs=rownames(exposures)
  sig_profiles=HDP_multi@comp_categ_distn$mean
  
  exposures_df<-as.data.frame(t(exposures),stringsAsFactors=F)%>%
    tibble::rownames_to_column("branch")%>%
    tidyr::separate(col="branch",into=c("node","Pair"),sep="_")%>%
    mutate(node=as.numeric(node))
  return(exposures_df)
}

mut_mat_HDP_comp=function(HDP_multi,ymax=0.2,plot=T){
  require(MutationalPatterns)
  sig_profiles=t(mut_example_multi@comp_categ_distn$mean)
  colnames(sig_profiles)<-paste0("N",0:(ncol(sig_profiles)-1))
  bases=c("A","C","G","T")
  subs=c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")
  rownames(sig_profiles)<-paste0(rep(rep(bases,each=4),times=6),rep(subs,each=16),rep(bases,times=24))
  if(plot){
    plot_96_profile(sig_profiles,ymax=ymax,condensed = T) 
  }
  return(sig_profiles)
}

#View signatures per sample
get_signatures_in_samples=function(tree,signature_names,exposures_df) {
  df<-dplyr::bind_cols(lapply(signature_names,function(signature) {
    sigs_in_samples=sapply(tree$tip.label,function(sample) {
      sample_node=which(tree$tip.label==sample)
      nodes_included=get_ancestral_nodes(node = sample_node,edge = tree$edge,exclude_root = T)
      branch_lengths=sapply(nodes_included,function(node) tree$edge.length[tree$edge[,2]==node])
      branch_sig_prop=sapply(nodes_included,function(node) {ifelse(node%in%exposures_df$node,sum(exposures_df[exposures_df$node==node,signature]),NA)})
      overall_contribution=weighted.mean(x=branch_sig_prop[!is.na(branch_sig_prop)],w = branch_lengths[!is.na(branch_sig_prop)])
      if(is.nan(overall_contribution)){overall_contribution<-0}
      return(overall_contribution)
    })
    return(sigs_in_samples)
  }))
  colnames(df)<-signature_names
  return(cbind(data.frame("Sample"=tree$tip.label),df))
}

drivers_per_sample_data=function(tree,details){
  require(dplyr)
  if(is.null(details$is.driver)){stop(print("Need variable 'is.driver' in the details matrix"))}
  driver_nodes=details%>%dplyr::filter(is.driver==1)%>%pull(node)
  n_drivers<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    return(sum(driver_nodes%in%tip_nodes))
  })
  driver_ids<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    driver_ids<-details%>%dplyr::filter(node %in% tip_nodes & is.driver==1)%>%pull(variant_ID)
    return(paste0(driver_ids,collapse=","))
  })
  out_df<-data.frame(Sample=tree$tip.label,n_drivers=n_drivers,driver_ids=ifelse(nchar(driver_ids)==0,NA,driver_ids))
  out_df<-out_df%>%dplyr::filter(Sample!="Ancestral")
  return(out_df)
}

