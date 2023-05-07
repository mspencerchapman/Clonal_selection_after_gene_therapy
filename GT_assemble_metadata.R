#TECHNICAL SEQUENCING METADATA
#Import the coverage/ sequencing run information information
#Parse all the pipeline report xls files
root_dir="~/R_work/Gene_therapy/Gene_therapy_for_SCD_NEJM"
files=list.files(paste0(root_dir,"/Data"),pattern = "Cancer_Pipeline_Reports",full.names = T)
out_list<-lapply(files,function(file){
  print(file)
  project<-as.numeric(stringr::str_extract(tail(stringr::str_split(file,pattern="/")[[1]],n=1), "(\\d)+"))
  print(project)
  nsheet=length(readxl::excel_sheets(file))
  print("Parsing info")
  seq_info_list<-lapply(4:nsheet,function(sheet){
    sample<-readxl::read_xls(file,sheet=sheet,col_names=F,range = "B1")
    seq_run<-readxl::read_xls(file,sheet=sheet,col_names=F,range = "F10")
    lane<-readxl::read_xls(file,sheet=sheet,col_names=F,range = "G10")
    return(data.frame(Sample=sample,run=ifelse(nrow(seq_run)==0,NA,seq_run),lane=ifelse(nrow(lane)==0,NA,lane)))
  })
  seq_info_list<-lapply(seq_info_list,function(df) {colnames(df)<-c("Sample","Run","Lane");return(df)})
  print("Rearranging info")
  seq_info_df<-dplyr::bind_rows(seq_info_list)
  seq_info_df<-seq_info_df%>%
    tidyr::separate(Lane,into=c("Lane_number","Sample_number"),remove=F,sep="#")%>%
    mutate(Run_Lane=paste(Run,Lane_number,sep="_"))
  
  coverage_info<-readxl::read_xls(file,sheet=1,skip = 2)[,c(1,3)]
  colnames(coverage_info)<-c("Sample","Coverage")
  seq_info_df<-left_join(coverage_info,seq_info_df)
  seq_info_df$Project=project
  return(seq_info_df)
})

all_projects_df<-dplyr::bind_rows(out_list)

#COLONY METADATA (Colony type, time point, sample type)
#Import the sample metadata spreadsheet
metadata_path=paste0(root_dir,"/Data/Sample_metadata.xlsx")
metadata_df<-readxl::read_excel(metadata_path)%>%dplyr::select(plate_ID,position,Individual,gender,'Donor ID',sample_ID_in_COSMIC)

#Filter out the BCH data only, and sort out all of the fields. This is all a bit messy as the data is very inconsistent in its input
metadata_srt<-metadata_df%>%
  dplyr::filter(grepl("BCH",`Donor ID`))%>%
  mutate(ID=stringr::str_extract(plate_ID,pattern="BCL-([0-9]*)"))%>%
  mutate(ID=stringr::str_replace(ID,"-",""))%>%
  mutate(Plate_no=stringr::str_extract(plate_ID,"Plate(\\s*)[0-9]"))%>%
  mutate(Plate_no=stringr::str_remove(Plate_no,"\\s"))%>%
  mutate(HSPC_source=ifelse(grepl("mPB CD34+",plate_ID),"mPB_CD34+",ifelse(grepl("UMPN",plate_ID),"UMPN",ifelse(grepl("BMA",plate_ID),"BMA",ifelse(grepl("CD34+",plate_ID),"CD34+",ifelse(grepl("PBMC",plate_ID),"PB","Unknown"))))))%>%
  mutate(temp=stringr::str_extract(plate_ID,pattern="BCL-(.*)"))%>%
  mutate(temp=stringr::str_replace(temp,pattern = "BCL-([0-9]*)(.)",replacement = ""))%>%
  separate(temp,into = c("Sample_type","Cell_type"),sep="_")%>%
  mutate(Cell_type=stringr::str_trim(Cell_type))%>%
  mutate(Time_point=ifelse(grepl("12",Sample_type),12,ifelse(grepl("24|2 years|2 yeard|2YEARS",Sample_type),24,ifelse(grepl("2.5years|2.5YEARS",Sample_type),30,ifelse(grepl("3 years",Sample_type),36,ifelse(grepl("21 Months",Sample_type),21,0))))))%>%
  mutate(Sample_type=ifelse(grepl("[0-9]+",Sample_type),"Post-GT",Sample_type))%>%
  mutate(Cell_type=ifelse(grepl("MIXED",Cell_type)&grepl("NON-BFUE",Cell_type),"MIXED & NON-BFUE",Cell_type))%>%
  mutate(Cell_type=stringr::str_replace(Cell_type,pattern="#.*",replace=""))%>%
  dplyr::rename("Sample"=sample_ID_in_COSMIC)%>%
  dplyr::select(-`Donor ID`,-Individual)


#COMBINE TECHNICAL & COLONY METADATA
all_sample_level_metadata<-full_join(all_projects_df,metadata_srt)
write.table(all_sample_level_metadata,file=paste0(root_dir,"/Data/Sample_level_metadata.tsv"),quote=F,sep="\t",row.names = F)

GT_samples_already<-read_csv("~/Downloads/GT_samples_for_Hyunchul.csv")

all_sample_level_metadata%>%
  filter(!Sample%in% GT_samples_already$Sample & !is.na(Coverage) & Coverage>0)%>%
  dplyr::select(Sample,Project)%>%
  write_csv(file="~/Desktop/GT_samples_for_Hyunchul_new.csv")
