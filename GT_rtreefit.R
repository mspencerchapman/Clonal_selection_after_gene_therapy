#Script to run "rtreefit" on the farm
#This is a SLOOOW mcmc approach, so best to run on the compute farm

#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------

cran_packages=c("ape","dplyr","tidyr")
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("rtreefit", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/rtreefit")
  library("rtreefit",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Set file paths and import data
#----------------------------------

all_tree_data_path="combined_tree_files.Rds"
metadata_path="sample_metadata_full_with_sigs.Rds"

all_tree_data=readRDS(all_tree_data_path)
sample_metadata<-readRDS(metadata_path)

#----------------------------------
# Get version of tree for rtreefit
#----------------------------------

all.trees.age.nodups<-Map(tree=all_tree_data$all.trees.cc.nodups,this_ID=names(all_tree_data$all.trees.cc.nodups),function(tree,this_ID) {
  cat(this_ID,sep="\n")
  
  #Define the function to correct branch lengths for a certain number of in vitro mutations
  correct_for_invitro_muts=function(tree,n_invitro=30){
    tree_new=tree
    tree_new$coords<-NULL
    tree_new$edge.length=mapply(curr_length=tree$edge.length,node=tree$edge[,2],FUN=function(curr_length,node) {
      if(node%in%1:length(tree$tip.label)){
        return(max(curr_length-n_invitro,0))
      } else {
        return(curr_length)
      }
    })
    return(tree_new)
  }
  
  #Run the function on the tree
  tree.invitrocorrected=correct_for_invitro_muts(tree,n_invitro=10)
  
  #Create the age data frame ("agedf") required by rtreefit - the age of sampling of each colony
  tree.invitrocorrected$agedf<-data.frame(tip.label=tree.invitrocorrected$tip.label,age=sapply(tree.invitrocorrected$tip.label,function(sample){
    if(sample=="Ancestral"){
      return(1e-6) #near zero value for the ancestral branch
    } else {
      #Ages of other samples are the "age at gene therapy" + the "time point" (time post gene therapy)
      age<-sample_metadata%>%
        dplyr::filter(Sample==sample)%>%
        dplyr::select(Age_at_GT,Time_point)%>%
        rowSums(.)
      
      #Add in the 0.75 years corresponding to the 9 months of in utero growth (this is included in the model)
      total_time=age+0.75
      return(total_time)
    }
  }))
  
  #rtreefit needs integer branch lengths - therefore round these to closest integer
  tree.invitrocorrected$edge.length<-round(tree.invitrocorrected$edge.length)
  
  #Now run the rtreefit algorithm with default settings
  tree_u2<-fit_tree(tree.invitrocorrected,switch_nodes = c(),xcross = c(),niter = 10000,model = "poisson_tree")
  return(tree_u2)
})

#Saving the trees alone
cat("Saving the new tree data",sep="\n")
saveRDS(all.trees.age.nodups,file="temp.Rds")

#Add the new age trees into the all_tree_data object
cat("Combining the new tree data with existing tree data",sep="\n")
all_tree_data$all.trees.age.nodups<-lapply(age.tree=all.trees.age.nodups,function(age.tree) return(age.tree$ultratree))

#Resave the update tree data
cat("Saving the updated tree data",sep="\n")
saveRDS(all_tree_data,file = all_tree_data_path)

#Plot the trees to check structure
cat("Plotting the new trees",sep="\n")
pdf("all_age_trees.pdf",width=15,height=10)
lapply(all.trees.age.nodups,function(output) plot.phylo(output$ultratree,direction="downwards",show.tip.label=F))
dev.off()

