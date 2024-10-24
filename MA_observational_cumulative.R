# Cumulative meta analysis of observational studies: by repeatedly fitting the specified model adding one study at a time
# Load the packages:
library(readxl)
library(DescTools)
library(metafor)
library(meta)
library(dplyr)
library(rio)
library(grid)
library(tidyverse)

setwd("./syst_rev/Results")

#1.  read file
#study_tool <- read_excel("~/study_tool.xlsx",sheet = "Work_sheet_v4")
study_tool <- import("~/study_tool.xlsx",sheet = "Work_sheet_v4")
str(study_tool)

# sort by year
study_tool_sort <- study_tool[order(study_tool$Year),]

# filter case-control designs
study_tool_cohort<- dplyr::filter(study_tool_sort, Study_design == "Cohort")

# filter for BC
study_tool_bc<- dplyr::filter(study_tool_cohort, Outcome == "Breast cancer")

# filter for hormone receptors
study_hr_bc<- dplyr::filter(study_tool_bc, hormone_receptors == 'NA')

# filter out BMI strata
study_bmi_bc<- dplyr::filter(study_hr_bc, BMI == 'NA')

# filter for menopause status
study_meno_bc<- dplyr::filter(study_bmi_bc, Menopause_status == "All")

# select effect measures
study_effects<- dplyr::select(study_meno_bc,"Unique_id","first_author","Year","Name_of_study","Exposure_recoded","cases", "N","Ref_group", "Test_group_recoded",
                              "Menopause_status", "hormone_receptors", "summary_effect", "Estimate",
                              "LCI_95", "UCI_95", "Country","Adequacy_adjustment","Highest_quantile")

# set as numeric columns
cols.num <- c("cases","Estimate","LCI_95","UCI_95")
study_effects[cols.num] <- sapply(study_effects[cols.num],as.numeric)
sapply(study_effects, class)
study_effects<- as.data.frame(study_effects)
str(study_effects)

# Reverse log estimates
yi <- log(study_effects$Estimate) 
sei <- (log(study_effects$LCI_95) - log(study_effects$UCI_95)) / (2*1.96)
SE<- abs((study_effects$LCI_95 - study_effects$UCI_95) / 3.92)

# store yi and sei in data set 
study_effects$yi<- yi
study_effects$sei<- sei
study_effects$SE<- SE

# filter 0 variance in yi and sei
study_effects<- dplyr::filter(study_effects, yi != 0 & sei!= 0 & SE!= 0)
dim(study_effects) #254

# filter only highest quantile
study_effects<- dplyr::filter(study_effects, Highest_quantile == 'Yes') #study_effects[study_effects$Highest_quantile == 'Yes',]
dim(study_effects) #77

# Exposures as glycemic index
study_effects_GI <- dplyr::filter(study_effects, Exposure_recoded == "Glycemic index" )
# include
dput(study_effects_GI$Unique_id)
include_GI<- c(279, 357, 394, 467, 524, 590, 663, 710, 839, 848, 948, 
               978, 1096) #,1008,157
study_effects_GI<- study_effects_GI[study_effects_GI$Unique_id %in% include_GI,]
str(study_effects_GI) #13

# Metanalysis of GI
# fit random-effects models
res_GI <- rma(yi, sei=sei, data=study_effects_GI, method="REML", slab=paste(first_author, Year, sep=", "))
summary(res_GI)

# cumulative meta-analysis (sorted by publication year)
tmp_GI <- cumul(res_GI, order=Year)

# cumulative forest plot
pdf(file= "./syst_rev/Results/forestplot_res_GI_cum.pdf")
forest(tmp_GI, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2)),shade=TRUE, atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")
grid.text("GI cumulative and BC", .5, .9, gp=gpar(cex=2))
dev.off()

# Exposures as glycemic load
study_effects_GL <- dplyr::filter(study_effects, Exposure_recoded == "Glycemic load" )
# include
dput(study_effects_GL$Unique_id)
include_GL<- c(162, 462, 389, 715, 993, 528, 610, 1111, 362, 928, 863, 1023, 
               294, 666, 842)
study_effects_GL<- study_effects_GL[study_effects_GL$Unique_id %in% include_GL,]
str(study_effects_GL)

# Metanalysis of GL
# fit random-effects models
res_GL <- rma(yi, sei=sei, data=study_effects_GL, method="REML", slab=paste(first_author, Year, sep=", "))
summary(res_GL)

# cumulative meta-analysis (sorted by publication year)
tmp_GL <- cumul(res_GL, order=Year)

# cumulative forest plot
pdf(file= "./syst_rev/Results/forestplot_res_GL_cum.pdf")
forest(tmp_GL, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2)),shade=TRUE, atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")
grid.text("GL cumulative and BC", .5, .9, gp=gpar(cex=2))
dev.off()

# Exposures as fiber
study_effects_fiber <- dplyr::filter(study_effects, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
# 2.20 include
dput(study_effects_fiber$Unique_id)
include_fiber<- c(564, 1130, 977, 1042, 339, 693, 490, 847)
study_effects_fiber<- study_effects_fiber[study_effects_fiber$Unique_id %in% include_fiber,]
str(study_effects_fiber)

# Metanalysis of fiber intake
# fit random-effects models
res_fiber <- rma(yi, sei=sei, data=study_effects_fiber, method="REML", slab=paste(first_author, Year, sep=", "))
summary(res_fiber)

# cumulative meta-analysis (sorted by publication year)
tmp_fiber <- cumul(res_fiber, order=Year)

# cumulative forest plot
pdf(file= "./syst_rev/Results/forestplot_res_fiber_cum.pdf")
forest(tmp_fiber, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2)),shade=TRUE, atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")
grid.text("Fiber intake cumulative and BC", .5, .9, gp=gpar(cex=2))
dev.off()

# Exposures as Total sugars
study_effects_ts <- dplyr::filter(study_effects, Exposure_recoded %in% c("Total sugars"))
# 2.8 include
dput(study_effects_ts$Unique_id)
include_ts<- c(720, 923, 1071, 188)
study_effects_ts<- study_effects_ts[study_effects_ts$Unique_id %in% include_ts,]
str(study_effects_ts) #4

# Metanalysis of Total sugars
# fit random-effects models
res_ts <- rma(yi, sei=sei, data=study_effects_ts, method="REML", slab=paste(first_author, Year, sep=", "))
summary(res_ts)

# cumulative meta-analysis (sorted by publication year)
tmp_ts<- cumul(res_ts, order=Year)

# cumulative forest plot
pdf(file= "./syst_rev/Results/forestplot_res_ts_cum.pdf")
forest(tmp_ts, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2)),shade=TRUE, atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")
grid.text("Sugar intake cumulative and BC", .5, .9, gp=gpar(cex=2))
dev.off()

# Exposures as SSBs
study_effects_SSB <- dplyr::filter(study_effects, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)","Artificially sweetened beverages (ASB)", "Soft drinks", "Industrially produced juices","Sugary drinks"))
# 2.2 include
dput(study_effects_SSB$Unique_id)
include_SSB<- c(121,650,797,799) 
study_effects_SSB<- study_effects_SSB[study_effects_SSB$Unique_id %in% include_SSB,]
str(study_effects_SSB) #4

# Metanalysis of SSBs
# fit random-effects models
res_SSB<- rma(yi, sei=sei, data= study_effects_SSB, method="REML", slab=paste(first_author, Year, sep=", "))
summary(res_SSB)

# cumulative meta-analysis (sorted by publication year)
tmp_SSB<- cumul(res_SSB, order= Year)

# cumulative forest plot
pdf(file= "./syst_rev/Results/forestplot_res_SSB_cum.pdf")
forest(tmp_SSB, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2)),shade=TRUE, atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")
grid.text("Sugar-sweetened beverages (SSB) cumulative and BC", .5, .9, gp=gpar(cex=2))
dev.off()
############################## END ##############################################
