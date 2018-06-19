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
save(ov.mir.isoform, file = "mir_isoform.RData")
miR.isoform.patient.ID<- (unique(ov.mir.isoform[,1])) # 461
miR.isoform.gene.ID <- (unique(ov.mir.isoform[,2])) # 632

sum(ov.mir.isoform[miR.isoform.patient.ID[1], 5])
sum(ov.mir.isoform[miR.isoform.patient.ID[4], 5])

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

#### Figure out the size for each patent vector (composed of all miR_isoform)
ind.finder.miRiso<- function(x.id){
  # setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR")
   load("mir_isoform.RData")
   ind.x.id<- (which(ov.mir.isoform[,1] %in% x.id))
   return(ind.x.id)
}

ind.for.patient <- sapply(miR.isoform.patient.ID, ind.finder.miRiso)
names(ind.for.patient) <- miR.isoform.patient.ID
#### end

length.for.patient <- sapply(ind.for.patient, length) # 461
X11()
hist(length.for.patient)
length(unique(length.for.patient))# 421

length(unique(colnames(ov.mir.isoform[2]))) # 421 

#### Here I want to make sure each patient has readings for all of 421 miRNAs
# recall: miR.isoform.patient.ID; ov.mir.isoform; ind.for.patient   
# colnames(ov.mir.isoform)
# [1] "SampleId" [2]"miRNA_ID"   [3]"isoform_coords"                
# [4] "read_count"  [5]"reads_per_million_miRNA_mapped" [6] "cross.mapped" [7] "miRNA_region"                  
patient.1.miR <- ov.mir.isoform[ind.for.patient[[1]], c(1, 2, 5)]; 
patient.2.miR <- ov.mir.isoform[ind.for.patient[[2]], c(1, 2, 5)] 

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
    length(miR.collapse.list[[2]])
    names(miR.collapse.list[[2]])
    
    str(miR.collapse.list)
    
    
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
  miR.collapse.df[valid.inds, i] <- miR.collapse.list[[i]]
  len.invalid <- nrow(miR.collapse.df)-length(valid.inds)
  miR.collapse.df[-valid.inds, i] <- rep(NA, len.invalid)
}

my.sum.rmNA <- function(my.vet){
  rm.ind <- which(is.na(my.vet)==TRUE)
  my.sum <- sum(my.vet[-rm.ind])
  return(my.sum)
}
apply(miR.collapse.df, 2, my.sum.rmNA) # 1e+06

colnames(miR.collapse.df) <- substr(colnames(miR.collapse.df), 1, 12)

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

miR.collapse.df.working <- miR.collapse.df[keep.miR.ind,]
miR.collapse.df.rm <- miR.collapse.df[rm.miR.ind,]

miR.collapse.df.working.percent <- miR.collapse.df.working/1000000

####
save(miR.collapse.df, miR.collapse.list, 
     miR.isoform.patient.ID, ov.mir.isoform, ind.for.patient, miR.isoform.gene.ID,
     miR.collapse.df.working, miR.collapse.df.rm, miR.collapse.df.working.percent, OV.clin, OV.clin.working,
     file="iso_miR_preprocess_OV_2018.RData"
)

load("iso_miR_preprocess_OV_2018.RData")

#.....................................................................................................................


#.....................................................................................................................
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
# setwd("C:/Users/we40947/Desktop/myCode/self_study/OV_miRseq/TCGA_gadc_broad_OV_clinical")
OV.clin.merge <- read.delim("OV_clin_merged.txt", header = FALSE, fill = TRUE)
rownames(OV.clin.merge) <- OV.clin.merge[,1]
OV.clin.merge <- OV.clin.merge[, -1] # 3511  591
patient.ID <- toupper(as.character(unlist(OV.clin.merge[19, ])))

# Note: there are 453 patient ID common in both miR.collapse.df and OV.clin.merge 
length(intersect(patient.ID, colnames(miR.collapse.df.working))) # 453

#----------------------------------------------------
#### working on vital-status: vital.status.collapsed # 591
#----------------------------------------------------
vital.status.ind <- grep("vital_status",rownames(OV.clin.merge))
   #  625 660 695 738 854
vital.status.mat <- OV.clin.merge[vital.status.ind, ]
colnames(vital.status.mat) <- patient.ID

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
names(vital.status.collapsed) <- patient.ID 
length(which(vital.status.collapsed=="alive"))

#----------------------------------------------------
#### working on days-to-death :  days.to.death.collapsed #
#----------------------------------------------------
days.to.death.ind <- grep("days_to_death", rownames(OV.clin.merge))
# 24 587 630 665 700
days.to.death.mat <- OV.clin.merge[days.to.death.ind, ]
colnames(days.to.death.mat) <- patient.ID
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
    days.to.death.collapsed <- c(days.to.death.collapsed, "NA")
  }
}
names(days.to.death.collapsed) <- patient.ID
length(which(days.to.death.collapsed=='NA')) # 248
days.to.death.collapsed.NUM <- as.numeric(days.to.death.collapsed); names(days.to.death.collapsed.NUM)<-patient.ID
#---------------------------------------
# working on last.followUp : last.followup.collapsed; last.followup.collapsed.NUM
#---------------------------------------
last.follow.up.ind <- grep("last_followup", rownames(OV.clin.merge))
# 26 588 631 666 701
last.follow.up.mat <- OV.clin.merge[last.follow.up.ind, ]
colnames(last.follow.up.mat) <- patient.ID
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
    last.followup.collapsed <- c(last.followup.collapsed, 'NA')
  }
}

length(which(last.followup.collapsed=='NA')) # 76
names(last.followup.collapsed) <- patient.ID
last.followup.collapsed.NUM <- as.numeric(last.followup.collapsed)
names(last.followup.collapsed.NUM) <- patient.ID

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
   primary.site <- as.character(unlist(OV.clin.merge[site.ind,])); names(primary.site)<-patient.ID; table(primary.site) 
   primary.patient.ind <- which(primary.site=="ovary")
   primary.patient.id <- patient.ID[primary.patient.ind]
   
grade.ind <- grep("grade", rownames(OV.clin.merge)) # 755  
   grade.patient <- as.character(unlist(OV.clin.merge[755,]))
   names(grade.patient) <- patient.ID
   table(grade.patient) # g1 6; g2 67; g3 491; g4 1; gb 2; gx 10
   grade.3.patient.ID <-patient.ID[which(grade.patient == "g3")]

# Recall:
   grade.patient   
   primary.site
   last.followup.collapsed.NUM
   days.to.death.collapsed.NUM <- as.numeric(days.to.death.collapsed); names(days.to.death.collapsed.NUM)<-patient.ID
   vital.status.collapsed
# End: Recall
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

OS.time.vec <- OV.clin.working$days.to.death.collapsed.NUM; names(OS.time.vec) <- rownames(OV.clin.working)
OS.time.vec[which(OV.clin.working$vital.status.collapsed=="alive")] <- OV.clin.working[which(OV.clin.working$vital.status.collapsed=="alive"), 
                                                                                       "last.followup.collapsed.NUM"] 
OV.clin.working$OS.time <- OS.time.vec

     # save(OV.clin, OV.clin.working, file="OV_clin_processed_2018.RData")

#..........................................................................................................................
# recall: miR.collapse.df.working (623x461); OV.clin.working (577x6)

mir.var <- apply(miR.collapse.df.working, 1, var)
log.mir.var <- apply(log2(miR.collapse.df.working), 1, var)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 0.7749 11.7403 16.1778 15.3213 19.2640 38.9895 

X11()
par(mfrow=c(2,2))
hist(unlist(miR.collapse.df.working), main="miR-TPM")
hist(log2(unlist(miR.collapse.df.working)), main="Log(miR-TPM)")
hist(mir.var, main="Var(miR-TPM)")
hist(log.mir.var, main="Var(Log(miR-TPM))")


common.pat.ID <- intersect(colnames(miR.collapse.df.working), rownames(OV.clin.working)) # 451

table(grade.patient);    rm.na.grade.id <- (names(grade.patient))[is.na(grade.patient)==TRUE]
grade.3.id <- (names(grade.patient))[grade.patient=="g3"]
real.grade.3.id <- grade.3.id[!(is.na(grade.3.id)==TRUE)] # 491

table(primary.site)
ov.primary <- primary.site[primary.site=="ovary"]
ov.primary.id <- names(ov.primary) # 576

## Get the overlapped id: grade3 & ov.primary.site & miR.collapse.df.working & OV.clin.working
ov.pri.g3.clin.id <- intersect(real.grade.3.id, ov.primary.id) # 487
clin.working.id <- rownames(OV.clin.working) # 577
clin.working.ov.g3.id <- intersect(ov.pri.g3.clin.id, clin.working.id) # 484

miR.working.id <- colnames(miR.collapse.df.working) # 461

miR.clin.working.id <- intersect(clin.working.ov.g3.id, miR.working.id) # 382

#### miR data and clinical data used for Overall Survival analysis
miR.exp.clin.working <- miR.collapse.df.working[, miR.clin.working.id] # 298x382
clin.miR.working <- OV.clin.working[miR.clin.working.id,] # 382x6
  save(miR.exp.clin.working, clin.miR.working, file="miR_clin_working_data_2018.RData")
  load("miR_clin_working_data_2018.RData")

identical(colnames(miR.exp.clin.working), rownames(clin.miR.working))

#..........................................................................................................................
# Recall: the data set for grade 3 OV with primary site at ovary
#        miR-seq data: miR.exp.clin.working
#        clinical data: clin.miR.working
library(survival)
os.time.g3 <- clin.miR.working$OS.time   
names(os.time.g3) <- rownames(clin.miR.working) 
os.status.g3 <- clin.miR.working$vital.status.collapsed 
names(os.status.g3) <- rownames(clin.miR.working)
os.g3.event <- (os.status.g3 == "dead")
os.g3.object <- Surv(os.time.g3/30.5, os.g3.event)
save(os.g3.object, file="OV_OS_g3_obj_2018.RData")
   load("OV_OS_g3_obj_2018.RData")

#...... To run following code, you need to laod there .RData files:
    setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/code_update_from_office")
    load("iso_miR_preprocess_OV_2018.RData")
    load("miR_clin_working_data_2018.RData")
    load("OV_OS_g3_obj_2018.RData")
   
  my.coxph.fun <- function(x){
                   load("OV_OS_g3_obj_2018.RData")
                   res.cox <- coxph(os.g3.object ~ x)
                   sum.cox <- summary(res.cox)
                   return(c(sum.cox$coefficients[c(1,5)], sum.cox$conf.int[c(1, 3, 4)]))
                   }

# test: my.coxph.fun:   try1 <- my.coxph.fun(unlist(log.miR.clin.working[1,]))
#                       try2 <- coxph(os.g3.object ~ unlist(log.miR.clin.working[1,])); try2.summary <- summary(try2)



log.miR.clin.working <- log2(miR.exp.clin.working) ; X11(); hist(unlist(log.miR.clin.working), main="OV Log(miR-TPM) N=382 miR=298")
sig.miR.screen <- t(apply(log.miR.clin.working, 1, my.coxph.fun))
sig.miR.screen.df <- data.frame(sig.miR.screen); colnames(sig.miR.screen.df)<- c("coef", "p-val","Hazard.Ratio", "lower", "upper")
sig.mature.miR.id <- (rownames(sig.miR.screen.df))[which(sig.miR.screen.df[,"p-val"]<0.05)]
   sig.cox.res <- sig.miR.screen.df[sig.mature.miR.id,] 

sig.log.miR.clin.working<- log.miR.clin.working[sig.mature.miR.id,];
target.data <- sig.log.miR.clin.working #  15 382
target.data.ttl<- log.miR.clin.working

load("mRNA_RSEM_iso_scEST.RData")
length(intersect(colnames(target.data),  patient.id.mrna.RSEM.iso.scEST))

sig.log.miR.var <- apply(sig.log.miR.clin.working, 1, var)

# hsa-mir-100 hsa-mir-101-1 hsa-mir-124-3  hsa-mir-1292   hsa-mir-140   hsa-mir-144  hsa-mir-16-1  hsa-mir-16-2  hsa-mir-1975  hsa-mir-200a 
# 34.884573      7.546704      8.638511     12.886779     17.502987     18.109952     19.183207     18.930045     17.550291     22.178821 
# hsa-mir-205 hsa-mir-219-2    hsa-mir-28   hsa-mir-410   hsa-mir-506 
# 25.016456     17.934268     20.619174     11.246111     11.463072 

log.miR.var <- apply(log.miR.clin.working, 1, var)
 # summary(log.miR.var)
 #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.725  11.707  15.853  15.287  19.187  37.747 

##### data used for clustering is based on significant miRNA-TPM: 
data.dist <- dist(t(target.data.ttl))
hc.complete <- hclust(data.dist, method="complete")
# hc.average <- hclust(data.dist, method="average")
# hc.single <- hclust(data.dist, method="single")
X11()
# par(mfrow=c(1, 3))
plot(hc.complete, main="Complete Linkage", xlab="", sub="", cex=.25)
abline(h=37, col="red", lwd=2)
# plot(hc.average, main="Average Linkage", xlab="", sub="", cex=.9)
# plot(hc.single, main="single Linkage", xlab="", sub="", cex=.9)
hc.clusters<-cutree(hc.complete, 2)
# cutree(hc.average, 2)
# cutree(hc.single, 2)
table(hc.clusters)
clust.1.TPM <- names(which(hc.clusters==1)) # 293
clust.2.TPM <- names(which(hc.clusters==2)) # 89

my.color<- hc.clusters
my.color[clust.1.TPM] <- "brown"
my.color[clust.2.TPM] <- "turquoise"

scale.target.data.ttl <- t(scale(t(target.data.ttl)))

save(log.miR.clin.working, sig.log.miR.clin.working, clust.1.TPM, clust.1.TPM, sig.mature.miR.id, 
     file="Log_g3_OVp_miR_clustering_by_miR100_miR101_2018.RData")

##### heatmap
# <https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2>
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
# lmat <- rbind(c(3,4),c(1,2))
# lhei <- c(1.5, 4); lwid <- c(4, 0.5)
lmat <- rbind(c(4,3),c(2,1)) # 4: key, 2: rowdendrogram, 3: coldendrogram, 1: heatmap
lhei <- c(1.5, 4); lwid <- c(0.5, 4)
X11()
heatmap.2(scale.target.data, main="Clustering based on miR in Hi-Grade OVCA",
          density.info="none",
          trace="none",
          margins=c(5, 12),
          col=my_palette,
          # breaks=col_breaks,
          dendrogram="column",
          Colv <- as.dendrogram(hc.complete), 
          Rowv="NA",
          ColSideColors = my.color,
          labCol = NA,
          # keysize = 0.7,           # size of color key
          # key.par = list(cex=0.5), # text size in color key
          key.title = NA,
          keysize = 0.6,           
          key.par = list(mar=c(1,2,20,2),cex=0.4), # mar(bottom, left, up, right)
          # key.par = list(mar=c(1,2,20,2),cex=0.4), # mar(bottom, left, up, right)
          cexCol = 0.3,
          cexRow = 1.5,
          symm=F,symkey=F,symbreaks=T, scale="none",
          # lmat = lmat,
          lhei=lhei,lwid=lwid
          # labRow =   ; RowSideColors = as.character(cluster)        
) # colon: tomato; ovary: green; uterus: blue; purple: vagina

par(lend = 1)
legend("topright", inset = c(0.002, 0.13),
       legend = c("clust.1.TPM", "clust.2.TPM"),
       col = c("brown", "turquoise"), 
       lty = 1,
       lwd = 5,
       bty = "n",
       cex = 0.65
)



### Kaplan-Meier plot
# Recall: clust.1.TPM; clust.2.TPM
clust.TPM <- rep(0, ncol(target.data)); names(clust.TPM) <- colnames(target.data)
clust.TPM[which(names(clust.TPM) %in% clust.1.TPM)]<- rep(1, length(clust.1.TPM)) # level 1: 293 : cluster 1 "brown"
# level 0: 89 : cluster 2 "turquoise"
clust.TPM.f <- as.factor(clust.TPM)
levels(clust.TPM.f)

os.g3.object
# sub.os.event <- (sub.update2.clin$vital.status_collapsed=="dead")     
# sur.object <- Surv(sub.update2.clin$sur.time/30.5, sub.os.event)  
fit.clust <- survfit(os.g3.object ~ clust.TPM.f, conf.type = "log-log")
log.rank.p.single<-survdiff(os.g3.object ~ clust.TPM.f,rho = 0)
p.val.single<-pchisq(log.rank.p.single$chisq, length(log.rank.p.single$n)-1, lower.tail = FALSE) 
summary(fit.clust)
X11()
# par(mar=c(6, 4.5, 4.5, 2.1), mfrow=c(1, 3))
plot(fit.clust,  mark.time = TRUE,col=c("turquoise", "brown"),lwd = 4, 
     xlab="Months after Surgery", ylab="Overall Survival Probability", main= "",
     cex.lab=1.2)

legend(75, 0.9, c("cluster.1 (N=89)", "cluster.2 (N=293)") , fill =c("turquoise", "brown"), bty = "n", cex = 1.2, text.font = 2 # 1 normal 2 bold 3 italic 4 bold and italic
) 

l<-legend(75, 0.9, c("cluster.1 (N=89)", "cluster.2 (N=293)") , fill =c("turquoise", "brown"), bty = "n", cex = 1.2, text.font = 2 # 1 normal 2 bold 3 italic 4 bold and italic
) 

legend(l$rec$left+3, l$rect$top-l$rect$h*2/3, paste("Log Rank p =", round(p.val.single, digit=4)), bty = "n", cex = 1.2, text.font = 2)

#................................................................................................................
#### Make a forest-plot like plot: 
## <https://www.r-bloggers.com/forest-plot-with-horizontal-bands/>
# Recall: sig.cox.res
#                           "coef"         "p-val"        "Hazard.Ratio" "lower"        "upper"
library(forestplot)
clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")

tabletext <- cbind( c("miRNA", "\n",rownames(sig.cox.res)),
                    c("beta", "\n", round(as.numeric(sig.cox.res[,1]), digit = 2)),
                    c("Hazard Ratio", "\n", round(as.numeric(sig.cox.res[,3]), digit = 2)),
                    c("Lower", "\n", round(as.numeric(sig.cox.res[,4]), digit = 4)),
                    c("Upper", "\n", round(as.numeric(sig.cox.res[,5]), digit = 4)),
                    c("p-value", "\n", round(as.numeric(sig.cox.res[,2]), digit = 4))
)

X11()
forestplot(labeltext=tabletext, graph.pos=3,
           mean =c(NA, NA, as.numeric(sig.cox.res[,1])),
           lower =c(NA, NA, as.numeric(sig.cox.res[,4])),
           upper =c(NA, NA, as.numeric(sig.cox.res[,5])),
           title = "Screening miRNA based on hazard ratio from Cox-PH regression fitting", 
           col=fpColors(box="royalblue", lines="darkblue", zero = "tomato"),
           xlab="Hazard Ratio",
           zero=1, boxsize=0.25, colgap=unit(9,"mm"),graphwidth=unit(100, "mm"),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=0.68),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.4)),
           lwd.ci=1.2, ci.vertices=TRUE, ci.vertices.height = 0.04
           
)
  
#................................................................................................................ 
# Recall: log.miR.clin.working, sig.log.miR.clin.working, clust.1.TPM, clust.1.TPM, sig.mature.miR.id
sig.log.miR.clin.working
corr.working.sig.miR <- t(sig.log.miR.clin.working)
corr.working.sig.miR.df <- as.data.frame(corr.working.sig.miR)
library(corrplot)
sig.miR.corr<- cor(corr.working.sig.miR, method="spearman")
sig.miR.p.mat <- cor.mtest(corr.working.sig.miR)
X11()
par(mfrow = c(1,2))
corrplot(sig.miR.corr, method="number", type = "upper", cl.cex = 0.8, pch.col = "blue")
corrplot(sig.miR.corr, p.mat = sig.miR.p.mat$p, insig = "p-value", type ="upper", sig.level = -1, pch.col = "blue", cl.cex = 0.8)
X11()
plot(corr.working.sig.miR.df, main="TCGA-OvCa miRLog(TPM)", col="lightblue") # scatter plot to check the correlation between 2 genes pairwisely. 
X11()
plot(corr.working.sig.miR.df$`hsa-mir-100`, corr.working.sig.miR.df$`hsa-mir-101-1`, pch=19, col="blue", 
     xlab="miR-100", ylab="miR-101-1")

 ## zoom in the correlation plot: length(sig.mature.miR.id), 
colnames(corr.working.sig.miR.df)
X11()
par(mfrow=c(3, 5)) 
for (i in 1:(length(sig.mature.miR.id)-1)){
  plot(corr.working.sig.miR.df[,1], corr.working.sig.miR[, (i+1)], pch=19, col = "blue", lwd=0.6,
       xlab=sig.mature.miR.id[1], ylab=sig.mature.miR.id[i+1], cex.lab= 1.65 )
} 

X11()
par(mfrow=c(3, 5)) 
plot(corr.working.sig.miR.df[,2], corr.working.sig.miR[, 1], pch=19, col = "blue", lwd=0.6,
     xlab=sig.mature.miR.id[2], ylab=sig.mature.miR.id[1], cex.lab= 1.65 )
for (i in 1:(length(sig.mature.miR.id)-2)){
                   plot(corr.working.sig.miR.df[,2], corr.working.sig.miR[, (i+2)], pch=19, col = "blue", lwd=0.6,
                        xlab=sig.mature.miR.id[2], ylab=sig.mature.miR.id[i+2], cex.lab= 1.65 )
                   } 

  ### To generalize: test all miR in log.miR.clin.working
corr.working.miR.df <- data.frame(t(log.miR.clin.working)) # 382 patient x  298 miR
colnames(corr.working.miR.df) <- rownames(log.miR.clin.working)

corr.res.rho <- c()
corr.res.pval<- c()
for (i in (1:ncol(corr.working.miR.df))){
  corr.res.rho[i] <- (cor.test(corr.working.miR.df[,"hsa-mir-100"], corr.working.miR.df[, i], method ="spearman"))$estimate
  corr.res.pval[i] <- (cor.test(corr.working.miR.df[,"hsa-mir-100"], corr.working.miR.df[, i], method ="spearman"))$p.val
}
corr.res.df.miR.100 <- data.frame(corr.res.rho, corr.res.pval)
rownames(corr.res.df.miR.100) <- colnames(corr.working.miR.df)
colnames(corr.res.df.miR.100) <- c("rho_miR100", "pval_miR100")
miR100.corr.sig.df <- corr.res.df.miR.100[which(corr.res.df.miR.100$pval_miR100<0.001), ]
miR100.sig.gene.id <- rownames(miR100.corr.sig.df)
miR100.sig.gene.m <- (t(log.miR.clin.working[miR100.sig.gene.id,])) 
miR100.sig.gene.df <- data.frame(miR100.sig.gene.m)
colnames(miR100.sig.gene.df) <- colnames(miR100.sig.gene.m)

   library(corrplot)
   sig.miR.corr<- cor(miR100.sig.gene.df, method="spearman")
   sig.miR.p.mat <- cor.mtest(miR100.sig.gene.df)
   X11()
   par(mfrow = c(1,2))
   corrplot(sig.miR.corr, method="number", type = "upper", cl.cex = 0.8, pch.col = "blue", number.cex = 0.7, diag = FALSE)
   corrplot(sig.miR.corr, p.mat = sig.miR.p.mat$p, insig = "p-value", type ="upper", sig.level = -1, pch.col = "blue", 
            cl.cex = 0.8, diag = FALSE)
   
