library(stringr)
local_ref="~/Documents/Reference_files/genome.fa"
lustre_ref="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
my_working_directory=getwd()
#Set up functions
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

#Import the mutations from the different patients
setwd("/lustre/scratch119/casm/team154pc/ms56/gene_therapy")
patients=c("BCL002","BCL003","BCL004","BCL006","BCL008","BCL009")
muts_by_patient_list=vector("list",length = length(patients))
names(muts_by_patient_list)<-patients
for(patient in patients) {
  load(paste0("filtering_runs/annotated_muts/annotated_mut_set_",patient,"_m40_postMS_reduced_a_j_vaf_post_mix"))
  details<-filtered_muts$COMB_mats.tree.build$mat
  details=details[details$Mut_type=="SNV",]
  details=details[,c("Chrom","Pos","Ref","Alt","node")]
  details$node=paste(details$node,patient,sep = "_")
  colnames(details)=c("Chr","Pos","Ref","Alt","SampleID")
  muts_by_patient_list[[patient]]<-details
}
all_muts_SCD=dplyr::bind_rows(muts_by_patient_list)

run_with_MPN_phylos=T
if(run_with_MPN_phylos) {
  #Now import the MPN data
  MPN_patients<-readLines("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/NW_samples.txt")
  MPN_muts_by_patient_list=vector("list",length = length(MPN_patients))
  names(MPN_muts_by_patient_list)<-MPN_patients
  for(MPN_patient in MPN_patients) {
    load(paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/NW/filtered_muts_",MPN_patient))
    details<-filtered_muts$COMB_mats.tree.build$mat
    details=details[details$Mut_type=="SNV",]
    details=details[,c("Chrom","Pos","Ref","Alt","node")]
    details$node=paste(details$node,MPN_patient,sep = "_")
    colnames(details)=c("Chr","Pos","Ref","Alt","SampleID")
    MPN_muts_by_patient_list[[MPN_patient]]<-details
  }
  all_muts_MPN=dplyr::bind_rows(MPN_muts_by_patient_list)
}

#Bind the mutations into a single dataframe & extract the trinuc_mut_mat using Tim's function
if(exists("all_muts_MPN")){
  all_muts<-dplyr::bind_rows(all_muts_SCD,all_muts_MPN)
} else {
  all_muts<-all_muts_SCD
}
dim(all_muts)
trinuc_mut_mat=mutlist_to_96_contexts(all_muts)

#Save the tables into the HDPx folder for further analysis
samples=rownames(trinuc_mut_mat)
key_table=data.frame(Sample=samples,Patient=str_split(samples,pattern="_",simplify=TRUE)[,2])
setwd(my_working_directory)
write.table(trinuc_mut_mat,"trinuc_mut_mat.txt")
write.table(key_table,"key_table.txt")
