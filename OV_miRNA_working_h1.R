setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR")
load("OV_TCGA_miRNA_workingFiles.RData", ex1<- new.env())
ls.str(ex1)
# OV.miRNA.mature.working  OV.miRNA.raw.working  OV.miRNA.RPKM.working 
# OV.TCGA.clin.miRNA OV.TCGA.clin.miRNA.1  OV.TCGA.clin.working.1

dim(OV.TCGA.clin.working.1)

setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_mRNAseq_Preprocess_2016/OV_mRNAseq_preprocess_level3")
OV.mRNA.scale.est <- read.table( "OVun_mRNAseq_scaled_estimate.txt", header = TRUE, sep = "\t", row.names = 1)
sd(as.numeric(OV.mRNA.scale.est[3,]))
mean(as.numeric(OV.mRNA.scale.est[2,]))

#.....................................................................................................................
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_miRseq_isoform_level3/OV_miRseq_isoform_level3")
ov.mir.isoform <- read.table( "OV_miRseq_isoform_level3.txt", header = TRUE, sep = "\t")
# save(ov.mir.isoform, file = "mir_isoform.RData") #  2131335  7
miR.isoform.patient.ID<- (unique(ov.mir.isoform[,1])) # 461 # length(miR.isoform.patient.ID)
miR.isoform.gene.ID <- (unique(ov.mir.isoform[,2])) # 632 # length(miR.isoform.gene.ID)

sum(ov.mir.isoform[miR.isoform.patient.ID[1], 5])  # Note: miR.isoform.patient.ID[1] is treated as factor for row index here
sum(ov.mir.isoform[miR.isoform.patient.ID[2], 5])  #    ov.mir.isoform[as.numeric(miR.isoform.patient.ID[1]), 5]
sum(ov.mir.isoform[miR.isoform.patient.ID[3], 5])  # Note: factor for row index here
sum(ov.mir.isoform[miR.isoform.patient.ID[4], 5])
sum(ov.mir.isoform[miR.isoform.patient.ID[5], 5])

sum(ov.mir.isoform[as.numeric(miR.isoform.patient.ID[2]), 5]) 

#### Figure out whether sum of "reads_per_millin_miRNA_mapped" for each patient will be 1M: answer is YES!
sum <- c()
for(i in 1:length(miR.isoform.patient.ID)){
  sum[i] <- sum(ov.mir.isoform[which(ov.mir.isoform[, 1] %in% miR.isoform.patient.ID[i]), 5])
  sum <- c(sum, sum[i])
}
# or: 
sum.miR.rpm  <- function(x.id){
        # setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR")
        load("mir_isoform.RData")
        inds <- which(ov.mir.isoform[,1] %in% x.id)
        sum.x.id<- sum(ov.mir.isoform[inds, 5])
        return(sum.x.id)
    }
summary.sum <- sapply(miR.isoform.patient.ID, sum.miR.rpm)
#### end

#### Figure out the size for each patient vector (composed of all miR_isoform)
ind.finder.miRiso<- function(x.id){
  # setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR")
   load("mir_isoform.RData")
   ind.x.id<- (which(ov.mir.isoform[,1] %in% x.id))
   return(ind.x.id)
}

ind.for.patient <- sapply(miR.isoform.patient.ID, ind.finder.miRiso)
names(ind.for.patient) <- miR.isoform.patient.ID
#### end

length.for.patient <- sapply(ind.for.patient, length) # 461 length(ind.for.patient)
    summary(length.for.patient)
    X11()
    hist(length.for.patient)
    length(unique(length.for.patient))# 421

length(unique(as.character(ov.mir.isoform[,2]))) # 421 # correction: 623 unique miR

#### Here I want to make sure each patient has readings for all of 421 miRNAs
# recall: miR.isoform.patient.ID; ov.mir.isoform; ind.for.patient   
# colnames(ov.mir.isoform)
# [1] "SampleId" [2]"miRNA_ID"   [3]"isoform_coords"                
# [4] "read_count"  [5]"reads_per_million_miRNA_mapped" [6] "cross.mapped" [7] "miRNA_region"                  
patient.1.miR <- ov.mir.isoform[ind.for.patient[[1]], c(1, 2, 5)]; length(unique(as.character(patient.1.miR[,2])))
patient.2.miR <- ov.mir.isoform[ind.for.patient[[2]], c(1, 2, 5)]; length(unique(as.character(patient.2.miR[,2])))  

identical(unique(patient.1.miR[,2]), unique(patient.2.miR[,2]))
length(unique(patient.1.miR[,2])) # 391
length(unique(patient.2.miR[,2])) # 332
intersect(unique(patient.1.miR[,2]), unique(patient.2.miR[,2]))
# so, between patient.1 and patient.2, the size of miR covered is different. 

#### create collapsed miR table with all patients' information, summarize reading of all isoform for each miR.
miR.collapse.list <- list()
for (i in 1:length(miR.isoform.patient.ID)){
  each.patient.miR <- NULL; miR.list <- NULL; miR.sum <- c()
  each.patient.miR <- ov.mir.isoform[ind.for.patient[[i]], c(1, 2, 5)] # extract miR for each patient 
  miR.list <- unique(each.patient.miR[,2])
  for (j in 1:length(miR.list)){
    miR.each.sum <- sum(as.numeric(each.patient.miR[which(each.patient.miR[,2] %in% miR.list[j]), "reads_per_million_miRNA_mapped"]))
    miR.sum <- c(miR.sum, miR.each.sum)
  }
  names(miR.sum) <- miR.list
  miR.collapse.list[[i]]<- miR.sum
}
length(miR.collapse.list); names(miR.collapse.list) <- miR.isoform.patient.ID
  # test: 
    miR.collapse.list[[2]] # for loop above dese work!!!
    length(miR.collapse.list[[2]]); length(miR.collapse.list[[1]])
    names(miR.collapse.list[[2]]); names(miR.collapse.list[[1]])
    
    str(miR.collapse.list)
    
#### use sqldf package to collapse ov.mir.isoform  
    colnames(ov.mir.isoform); try1 <- ov.mir.isoform; 
                              # colnames(try1) <- c('patID', 'miRID', 'readCount', 'RPM', 'xMap', 'miRegion') 
    # [1] "SampleId"                       "miRNA_ID"                       "isoform_coords"                
    # [4] "read_count"                     "reads_per_million_miRNA_mapped" "cross.mapped"                  
    # [7] "miRNA_region"        
                               sqldf("SELECT SampleID, count(*) AS Count_Patient FROM try1 GROUP BY SampleID"); class(try1); class(ov.mir.isoform); identical(ov.mir.isoform[1,], try1[1,])                           
    
    library(sqldf)
    try1_collapse <- sqldf("SELECT SampleID, miRNA_ID, SUM(reads_per_million_miRNA_mapped) FROM 'ov.mir.isoform' 
                           GROUP BY SampleID, miRNA_ID")
    
    collasp_sql_ovmiRisoform <- sqldf("SELECT  SampleID, miRNA_ID, SUM(reads_per_million_miRNA_mapped)  FROM 'ov.mir.isoform'  GROUP BY SampleID, miRNA_ID")
    
    length(unique(collasp_sql_ovmiRisoform[,1])) # 461 unique patients
    length(unique(collasp_sql_ovmiRisoform[,2])) # 623 unique miR
    
#### use sqldf package to collapse ov.mir.isoform (~~~ END ~~~)
    
    
    
# recall: miR.isoform.patient.ID; ov.mir.isoform; ind.for.patient; miR.isoform.gene.ID; miR.collapse.list   
length(miR.isoform.patient.ID) # 461
length(unique(miR.isoform.gene.ID))  # 623
    
miR.collapse.df <- data.frame(matrix(rep(0, (623*461)), nrow=623, ncol=461)) 

gene.ID.chr <- as.character(miR.isoform.gene.ID)
rownames(miR.collapse.df)<-gene.ID.chr   

patient.ID <- as.character(miR.isoform.patient.ID)
colnames(miR.collapse.df) <- patient.ID 

identical(names(miR.collapse.list), colnames(miR.collapse.df))

# update miR.collapse.df based on miR.collapse.list
for (i in 1:ncol(miR.collapse.df)){
  valid.inds <- NULL
  valid.inds<- (which(names(miR.collapse.list[[i]]) %in% rownames(miR.collapse.df)))
  miR.collapse.df[valid.inds, i] <- miR.collapse.list[[i]] # evaluate miR reads (not NA) to each patient i
  len.invalid <- nrow(miR.collapse.df)-length(valid.inds)  # determine total number or NA miR count for each patient i
  miR.collapse.df[-valid.inds, i] <- rep(NA, len.invalid) # evaluate miR read as NA (if appropriate) for each patient
}

my.sum.rmNA <- function(my.vet){
  rm.ind <- which(is.na(my.vet)==TRUE)
  my.sum <- sum(my.vet[-rm.ind])
  return(my.sum)
}
apply(miR.collapse.df, 2, my.sum.rmNA) # 1e+06

colnames(miR.collapse.df) <- substr(colnames(miR.collapse.df), 1, 12)
save(miR.collapse.df, file="TCGA_OV_miR_complete.RData") # 623 461
 # "C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_miRseq_isoform_level3/OV_miRseq_isoform_level3"

#### !!!! Now I get miR.collapse.df with 461 patients as columes and 623 miR as rows

#### Get the working version for miR expression for 461 patients: without NA
keep.miR.ind <- c(); rm.miR.ind<-c()
for (i in 1:nrow(miR.collapse.df)){
  if (length(which(is.na(miR.collapse.df[i,])==TRUE))==0) {
    keep.miR.ind <- c(keep.miR.ind, i)
  } else {
    rm.miR.ind <- c(rm.miR.ind, i)
  }
}
keep.miR.ind
rm.miR.ind

miR.collapse.df.working <- miR.collapse.df[keep.miR.ind,]; test_col_sum_working <- sapply(miR.collapse.df.working,sum)
                                                          
miR.collapse.df.rm <- miR.collapse.df[rm.miR.ind,]

miR.collapse.df.working.percent <- miR.collapse.df.working/1000000
 save(miR.collapse.df, miR.collapse.df.working, miR.collapse.df.working.percent, miR.collapse.df.rm, 
      collasp_sql_ovmiRisoform, file="TCGA_OV_miR_complete.RData")
 #  YOU CAN FIND THIS SAVED .RDATA @ C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_miRseq_isoform_level3/OV_miRseq_isoform_level3 
 
####
save(miR.collapse.df, miR.collapse.list, 
     miR.isoform.patient.ID, ov.mir.isoform, ind.for.patient, miR.isoform.gene.ID,
     miR.collapse.df.working, miR.collapse.df.rm, miR.collapse.df.working.percent, OV.clin, OV.clin.working, miR_clin_patient_ID,
     file="iso_miR_preprocess_OV.RData"
)
#.....................................................................................................................


#.....................................................................................................................
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
OV.clin.merge <- read.delim("OV_clin_merged.txt", header = FALSE, fill = TRUE)
rownames(OV.clin.merge) <- OV.clin.merge[,1]
OV.clin.merge <- OV.clin.merge[, -1] # 3511  591
patient.ID.ov.clin <- toupper(as.character(unlist(OV.clin.merge[19, ])))

save (OV.clin.merge, patient.ID.ov.clin, file="OV_clin_merge_raw.RData")
# Note: there are 453 patient ID common in both miR.collapse.df and OV.clin.merge 
length(intersect(patient.ID.ov.clin, colnames(miR.collapse.df.working))) # 453

#----------------------------------------------------
#### working on vital-status: vital.status.collapsed # 591
#----------------------------------------------------
vital.status.ind <- grep("vital_status",rownames(OV.clin.merge))
   #  625 660 695 738 854
vital.status.mat <- OV.clin.merge[vital.status.ind, ]
colnames(vital.status.mat) <- patient.ID.ov.clin

#...collapse vital.status
vital.status.collapsed <- c()
for (i in 1:dim(vital.status.mat)[2]){
  
  if (length(which(vital.status.mat[,i]=="dead")) > 0){
    
    vital.status.collapsed <- c(vital.status.collapsed, "dead")
    
  } else { 
    vital.status.collapsed <- c(vital.status.collapsed, "alive")  
  }   
}

table(vital.status.collapsed) # alive  dead 
                              #   247   344 
names(vital.status.collapsed) <- patient.ID.ov.clin 
length(which(vital.status.collapsed=="alive"))

#----------------------------------------------------
#### working on days-to-death :  days.to.death.collapsed #
#----------------------------------------------------
days.to.death.ind <- grep("days_to_death", rownames(OV.clin.merge))
# 24 587 630 665 700
days.to.death.mat <- OV.clin.merge[days.to.death.ind, ]
colnames(days.to.death.mat) <- patient.ID.ov.clin
days.to.death.mat[1, c("TCGA-72-4231", "TCGA-72-4232", "TCGA-72-4233",
                       "TCGA-72-4234", "TCGA-72-4235", "TCGA-72-4236",
                       "TCGA-72-4237", "TCGA-72-4238", 
                       "TCGA-72-4240", "TCGA-72-4241")] <- rep(NA, 10)
View(days.to.death.mat)
#...collapse days.to.death
days.to.death.collapsed <- c()
for(i in 1:ncol(days.to.death.mat)){
  
  if (sum((is.na(days.to.death.mat[,i])==TRUE)) < nrow(days.to.death.mat)){
    m.day.death <- min(as.numeric(as.character(unlist(days.to.death.mat[,i]))), na.rm=T)
    days.to.death.collapsed <- c(days.to.death.collapsed, m.day.death)
  } else {
    days.to.death.collapsed <- c(days.to.death.collapsed, NA)
  }
}
names(days.to.death.collapsed) <- patient.ID.ov.clin
length(which(is.na(days.to.death.collapsed))) # 248

#---------------------------------------
# working on last.followUp : last.followup.collapsed; last.followup.collapsed.NUM
#---------------------------------------
last.follow.up.ind <- grep("last_followup", rownames(OV.clin.merge))
# 26 588 631 666 701
last.follow.up.mat <- OV.clin.merge[last.follow.up.ind, ]
colnames(last.follow.up.mat) <- patient.ID.ov.clin
last.follow.up.mat[1, c("TCGA-72-4231", "TCGA-72-4232", "TCGA-72-4233",
                       "TCGA-72-4234", "TCGA-72-4235", "TCGA-72-4236",
                       "TCGA-72-4237", "TCGA-72-4238", 
                       "TCGA-72-4240", "TCGA-72-4241")] <- rep(NA, 10)
rownames(last.follow.up.mat)
View(last.follow.up.mat)

last.followup.collapsed <- c()
for (i in 1:ncol(last.follow.up.mat)) {
  
  if(sum(is.na(last.follow.up.mat[, i])) < nrow(last.follow.up.mat)){
    m.fl <- max(as.numeric(as.character(last.follow.up.mat[,i])), na.rm=T)
    last.followup.collapsed <- c(last.followup.collapsed, m.fl)
  } else {
    last.followup.collapsed <- c(last.followup.collapsed, NA)
  }
}

length(which(is.na(last.followup.collapsed))) # 76
names(last.followup.collapsed) <- patient.ID.ov.clin
last.followup.collapsed.NUM <- as.numeric(last.followup.collapsed)
names(last.followup.collapsed.NUM) <- patient.ID.ov.clin

#---------------------------------------
# working on tumor_progression; tumor_recurrence
#---------------------------------------
tmr_progression.ind <- grep("tumor_progression", rownames(OV.clin.merge)) # row 27 all NA
tmr_recurrence.ind <- grep("tumor_recurrence", rownames(OV.clin.merge)) # row 28  all NA
new_tumor.ind <- grep("new_tumor",rownames(OV.clin.merge))

#---------------------------------------
# working on primary.site; grade.patient 
#---------------------------------------
site.ind <- grep("tumor_tissue_site", rownames(OV.clin.merge)) # row 852
   table(as.character(unlist(OV.clin.merge[site.ind,])))
   primary.site <- as.character(unlist(OV.clin.merge[site.ind,])); names(primary.site)<-patient.ID.ov.clin; table(primary.site) 
   primary.patient.ind <- which(primary.site=="ovary")
   primary.patient.id <- patient.ID.ov.clin[primary.patient.ind]
   
grade.ind <- grep("grade", rownames(OV.clin.merge)) # 755  
   grade.patient <- as.character(unlist(OV.clin.merge[755,]))
   names(grade.patient) <- patient.ID.ov.clin
   table(grade.patient) # g1 6; g2 67; g3 491; g4 1; gb 2; gx 10
   grade.3.patient.ID <-patient.ID.ov.clin[which(grade.patient == "g3")]

# Recall:
   grade.patient   
   primary.site
   last.followup.collapsed.NUM
   days.to.death.collapsed.NUM <- as.numeric(days.to.death.collapsed); names(days.to.death.collapsed.NUM)<-patient.ID.ov.clin
   vital.status.collapsed
 OV.clin <- data.frame(grade.patient, primary.site, 
                       vital.status.collapsed,
                       days.to.death.collapsed.NUM,
                       last.followup.collapsed.NUM
                       )
 View(OV.clin)# 591
 
 rm.ind <- which((is.na(OV.clin$days.to.death.collapsed.NUM)==TRUE) & (is.na(OV.clin$last.followup.collapsed.NUM)==TRUE))
 # 13
 dead.keep.ind<- which((OV.clin$vital.status.collapsed =="dead") & (is.na(OV.clin$days.to.death.collapsed.NUM)!=TRUE))
 # 343
 alive.keep.ind <- which((OV.clin$vital.status.collapsed =="alive") & (is.na(OV.clin$last.followup.collapsed.NUM)!=TRUE)) 
 # 234
 
## there is one patient no belonging to (rm.ind, dead.keep.ind, alive.keep.ind)
clear.ind <-c(rm.ind, dead.keep.ind, alive.keep.ind)
which(!(c(1:591) %in% clear.ind)) # 166   
(rownames(OV.clin))[166] # "TCGA-24-0968"
 
OV.clin.working <- OV.clin[c(dead.keep.ind, alive.keep.ind),] # 577 
   
save (OV.clin.merge, patient.ID.ov.clin, OV.clin, OV.clin.working, miR_clin_patient_ID, file="OV_clin_merge_raw_plus_working.RData")
# @ "C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016"

###--------------------------------------------------------
# patient IDs that show up in both miR dataset and clinical dataset:
OV.clin.working
miR.collapse.df.working

miR_clin_patient_ID <- intersect(rownames(OV.clin.working), colnames(miR.collapse.df.working))

###--------------------------------------------------------
write.csv(OV.clin.working, file="OV_clin_TCGA_working.csv")
write.csv(OV.clin, file="OV_clin_TCGA.csv")
write.csv(miR.collapse.df.working, file="miR_collapse_df_working.csv")
write.csv(miR.collapse.df.working.percent, file="miR_collapse_df_working_percent.csv")
####-------------------------------------------------------
# miR_version:

clin_for_miR <- OV.clin.working[miR_clin_patient_ID,]   # dim 451*5
miR_ovrlap_clin <- miR.collapse.df.working[, miR_clin_patient_ID] # dim  298*451

save(clin_for_miR, miR_ovrlap_clin, file="miR_clin_Assay_data.RData")
write.csv(clin_for_miR, file="clin_for_miR.csv")
write.csv(miR_ovrlap_clin, file="miR_ovrlap_clin.csv")




