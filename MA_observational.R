# Meta analysis of observational studies (Sugar intake, added sugar, SSB,etc associated with BC)
# Load the package:
library(readxl)
library(DescTools)
library(metafor)
library(dplyr)
library(rio)
library(grid)
library(tidyverse)

setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results")

#1.  read file
#study_tool <- read_excel("~/Desktop/Skeleton_projects/nutrient_cancer/study_tool.xlsx",sheet = "Work_sheet_v4")
study_tool <- import("~/Desktop/Skeleton_projects/nutrient_cancer/study_tool.xlsx",sheet = "Work_sheet_v4")
str(study_tool)

# sort by year
#study_tool_sort <- study_tool[order(-study_tool$Year),]

# filter case-control designs
study_tool_cohort<- dplyr::filter(study_tool, Study_design == "Cohort")

# filter for BC
study_tool_bc<- dplyr::filter(study_tool_cohort, Outcome == "Breast cancer")

# filter for hormone receptors
study_hr_bc<- dplyr::filter(study_tool_bc, hormone_receptors == 'NA')

# filter out BMI strata
study_bmi_bc<- dplyr::filter(study_hr_bc, BMI == 'NA')

# filter for menopause 
study_meno_bc<- dplyr::filter(study_bmi_bc, Menopause_status == "All")

# select effect measures
study_effects<- dplyr::select(study_meno_bc,"Unique_id","Author_year","Name_of_study","Exposure_recoded","cases", "N","Ref_group", "Test_group_recoded",
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

#2. Filter according to study characteristics: BC irrespective of menopause status or hormone receptor
dput(unique(study_effects$Exposure_recoded))
# c("Sugary drinks", "Artificially sweetened beverages (ASB)", 
#   "Fruit juice", "Glycemic index", "Glycemic load", "Fiber intake", 
#   "Total sugars", "Added sugars", "Sugary foods", "Free sugar", 
#   "Fruit intake", "Fructose", "Sucrose", "Glucose", "Lactose", 
#   "Maltose", "Total fiber intake", "Sugar-sweetened beverages (SSB)", 
#   "Glycemic index: low", "Glycemic index: medium", "Glycemic index: high", 
#   "Whole grain intake", "Total fiber", "Fiber")

######## ######## Exclude redundant studies ######## ######## 
#NOTE: EPIC Italy contained in EPIC international
#NOTE: Cho E et al. 2003  NHII contained in Farvid MS et al. 2015
#NOTE: Holmes MD et al. 2004 NHI contained in Farvid MS et al. 2015

# 2.1 Exposures as SSB
study_effects_SSB <- dplyr::filter(study_effects, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)","Artificially sweetened beverages (ASB)", "Soft drinks", "Industrially produced juices","Sugary drinks"))
# 2.2 include
dput(study_effects_SSB$Unique_id)
include_SSB<- c(121,650,797,799) 
study_effects_SSB<- study_effects_SSB[study_effects_SSB$Unique_id %in% include_SSB,]
str(study_effects_SSB) #4

# 2.3 Exposures as glycemic index
study_effects_GI <- dplyr::filter(study_effects, Exposure_recoded == "Glycemic index" )
# 2.4 include
dput(study_effects_GI$Unique_id)
include_GI<- c(279, 357, 394, 467, 524, 590, 663, 710, 839, 848, 948, 
               978, 1096) #,1008,157
study_effects_GI<- study_effects_GI[study_effects_GI$Unique_id %in% include_GI,]
str(study_effects_GI) #13

# 2.5 Exposures as glycemic load
study_effects_GL <- dplyr::filter(study_effects, Exposure_recoded == "Glycemic load" )
# 2.6 include
dput(study_effects_GL$Unique_id)
include_GL<- c(294, 362, 389, 462, 528, 610, 666, 715, 842, 863, 928, 
               993, 1023, 1111) #162
study_effects_GL<- study_effects_GL[study_effects_GL$Unique_id %in% include_GL,]
str(study_effects_GL) #14

# 2.7 Exposures as Total sugars
study_effects_ts <- dplyr::filter(study_effects, Exposure_recoded %in% c("Total sugars"))
# 2.8 include
dput(study_effects_ts$Unique_id)
include_ts<- c(188, 720, 923, 1071)
study_effects_ts<- study_effects_ts[study_effects_ts$Unique_id %in% include_ts,]
str(study_effects_ts) #4

# 2.9 Exposures as sugary food
study_effects_sf <- dplyr::filter(study_effects, Exposure_recoded %in% c("Sugary foods","Sweet foods"))
# 2.10 include
dput(study_effects_sf$Unique_id)
include_sf<- c(196, 502, 640)
study_effects_sf<- study_effects_sf[study_effects_sf$Unique_id %in% include_sf,]
str(study_effects_sf) #3

# 2.11 Exposures as Added sugars
study_effects_as <- dplyr::filter(study_effects, Exposure_recoded %in% c("Added sugars"))
# 2.12 include
dput(study_effects_as$Unique_id)
include_as<- c(192, 918, 1076)
study_effects_as<- study_effects_as[study_effects_as$Unique_id %in% include_as,]
str(study_effects_as) #3

# 2.13 Exposures as common sugars: glucose, fructose, sucrose
study_effects_cs <- dplyr::filter(study_effects, Exposure_recoded %in% c("Glucose","Sucrose","Free sugar"))
# 2.14 include
dput(study_effects_cs$Unique_id)
include_cs<- c(204, 472, 634) #219 =sucrose /223=glucose
study_effects_cs<- study_effects_cs[study_effects_cs$Unique_id %in% include_cs,]
str(study_effects_cs) #3

# 2.15 Exposures as sugary drinks: SSB and fruit juice
# study_effects_sd <- dplyr::filter(study_effects, Exposure_recoded %in% c("Sugary drinks"))
# # 2.16 include
# dput(study_effects_sd$Unique_id)
# include_sd<- c(117, 200)
# study_effects_sd<- study_effects_sd[study_effects_sd$Unique_id %in% include_sd,]
# str(study_effects_sd) #2

# 2.17 Exposures as fruit juice
study_effects_fj <- dplyr::filter(study_effects, Exposure_recoded %in% c("Fruit juice"))
# 2.18 include
dput(study_effects_fj$Unique_id)
include_fj<- c(125, 269, 647)
study_effects_fj<- study_effects_fj[study_effects_fj$Unique_id %in% include_fj,]
str(study_effects_fj) #3

# 2.19 Exposures as fiber
study_effects_fiber <- dplyr::filter(study_effects, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
# 2.20 include
dput(study_effects_fiber$Unique_id)
include_fiber<- c(339, 490, 564, 693, 847, 977, 1042, 1130)#187
study_effects_fiber<- study_effects_fiber[study_effects_fiber$Unique_id %in% include_fiber,]
str(study_effects_fiber) #8

#### #### #### #### ####  MA #### #### #### #### #### 
#3. Random-effects model (using log risk ratios and variances as input)
res_SSB <- rma(yi, sei=sei, data= study_effects_SSB, method="REML",slab=Author_year,measure = "RR")
summary(res_SSB) # 0.16
reporter(res_SSB, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_GI <- rma(yi, sei=sei, data=study_effects_GI, method="REML",slab=Author_year,measure = "RR")
summary(res_GI) #* 0.024 I2: 0.03%. The more processed a food is, the higher its GI . p = 0.099 with DL estimator
#res_GI_table<- as.data.frame(coef(summary(res_GI)))
#export(res_GI_table, "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/res_GI_table.csv")
reporter(res_GI, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_GL <- rma(yi, sei=sei, data=study_effects_GL, method="REML",slab=Author_year,measure = "RR")
summary(res_GL) # 0.24
#res_GL_table<- as.data.frame(coef(summary(res_GL)))
#export(res_GL_table, "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/res_GL_table.csv")
reporter(res_GL, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_ts <- rma(yi, sei=sei, data=study_effects_ts, method="REML",slab=Author_year,measure = "RR")
summary(res_ts) # 0.5
reporter(res_ts, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_sf <- rma(yi, sei=sei, data=study_effects_sf, method="REML",slab=Author_year,measure = "RR")
summary(res_sf) # 0.4
reporter(res_sf, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_as <- rma(yi, sei=sei, data=study_effects_as, method="REML",slab=Author_year,measure = "RR")
summary(res_as) # 0.32
reporter(res_as, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_cs <- rma(yi, sei=sei, data=study_effects_cs, method="REML",slab=Author_year,measure = "RR")
summary(res_cs) # 0.35
reporter(res_cs, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

# res_sd <- rma(yi, sei=sei, data=study_effects_sd, method="REML",slab=Author_year,measure = "RR")
# summary(res_sd) # 0.0002
# reporter(res_sd, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_fj <- rma(yi, sei=sei, data=study_effects_fj, method="REML",slab=Author_year,measure = "RR")
summary(res_fj) # 0.37
reporter(res_fj, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

res_fiber <- rma(yi, sei=sei, data=study_effects_fiber, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber) # 0.0598
reporter(res_fiber, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")

# 3.1 prediction intervals for future values
pred_res_SSB<- predict(res_SSB, transf=exp, digits=2); pred_res_SSB
pred_res_GI<- predict(res_GI, transf=exp, digits=2); pred_res_GI #*
pred_res_GL<- predict(res_GL, transf=exp, digits=2); pred_res_GL
pred_res_ts<- predict(res_ts, transf=exp, digits=2); pred_res_ts
pred_res_sf<- predict(res_sf, transf=exp, digits=2); pred_res_sf
pred_res_as<- predict(res_as, transf=exp, digits=2); pred_res_as
pred_res_cs<- predict(res_cs, transf=exp, digits=2); pred_res_cs
pred_res_sd<- predict(res_sd, transf=exp, digits=2); pred_res_sd #*
pred_res_fj<- predict(res_fj, transf=exp, digits=2); pred_res_fj
# NOTE: * significant
pred_res_fiber<- predict(res_fiber, transf=exp, digits=2); pred_res_fiber

#4. Fixed-effects model 
fes_SSB <- rma(yi, sei=sei, data=study_effects_SSB, method="FE", slab=Author_year,measure = "RR") 
summary(fes_SSB) # 0.15
fes_GI <- rma(yi, sei=sei, data=study_effects_GI, method="FE", slab=Author_year,measure = "RR") 
summary(fes_GI) #* 0.0079. I2: 0%
fes_GI_table<- as.data.frame(coef(summary(fes_GI)))
#export(fes_GI_table, "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/fes_GI_table.csv")
fes_GL <- rma(yi, sei=sei, data=study_effects_GL, method="FE", slab=Author_year,measure = "RR") 
summary(fes_GL)#0.19
fes_GL_table<- as.data.frame(coef(summary(fes_GL)))
#export(fes_GL_table, "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/fes_GL_table.csv")
fes_ts <- rma(yi, sei=sei, data=study_effects_ts, method="FE", slab=Author_year,measure = "RR") 
summary(fes_ts) # 0.61
fes_sf <- rma(yi, sei=sei, data=study_effects_sf, method="FE", slab=Author_year,measure = "RR") 
summary(fes_sf) # 0.4
fes_as <- rma(yi, sei=sei, data=study_effects_as, method="FE", slab=Author_year,measure = "RR") 
summary(fes_as) # 0.29
fes_cs <- rma(yi, sei=sei, data=study_effects_cs, method="FE", slab=Author_year,measure = "RR") 
summary(fes_cs) # 0.35
fes_sd <- rma(yi, sei=sei, data=study_effects_sd, method="FE", slab=Author_year,measure = "RR") 
summary(fes_sd) #* 0.0002
fes_sd_table<- as.data.frame(coef(summary(fes_sd)))
#export(fes_sd_table, "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/fes_sd_table.csv")
fes_fj <- rma(yi, sei=sei, data=study_effects_fj, method="FE", slab=Author_year,measure = "RR") 
summary(fes_fj) # 0.37
fes_fiber <- rma(yi, sei=sei, data=study_effects_fiber, method="FE", slab=Author_year,measure = "RR") 
summary(fes_fiber) # 0.016

# 4.1 prediction intervals for future
pred_fes_SSB<- predict(fes_SSB, transf=exp, digits=2); pred_fes_SSB
pred_fes_GI<- predict(fes_GI, transf=exp, digits=2); pred_fes_GI #* 
pred_fes_GL<- predict(fes_GL, transf=exp, digits=2); pred_fes_GL
pred_fes_ts<- predict(fes_ts, transf=exp, digits=2); pred_fes_ts
pred_fes_sf<- predict(fes_sf, transf=exp, digits=2); pred_fes_sf
pred_fes_as<- predict(fes_as, transf=exp, digits=2); pred_fes_as
pred_fes_cs<- predict(fes_cs, transf=exp, digits=2); pred_fes_cs
pred_fes_sd<- predict(fes_sd, transf=exp, digits=2); pred_fes_sd # *
pred_fes_fj<- predict(fes_fj, transf=exp, digits=2); pred_fes_fj
pred_fes_fiber<- predict(fes_fiber, transf=exp, digits=2); pred_fes_fiber #*

#5. Forest plots
# helper function to add Q-test, I^2, and tau^2 estimate info into forest plot
mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

# RE-significants
# Glycemic index
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_res_GI.pdf")
forest(res_GI, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)),showweights=TRUE,
       header="First author(s) and Year",order=order(study_effects_GI$yi),xlim=c(-8,5))#mlab= mlabfun("RE Model", res_GI)
grid.text("GI and BC", .5, .9, gp=gpar(cex=2))

# funnel plot
funnel(trimfill(res_GI), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(res_GI)
plot(influence(res_GI), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(res_GI)

# Glycemic load
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_res_GL.pdf")
forest(res_GL, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)), showweights=TRUE,header="First author(s) and Year",order=order(study_effects_GL$yi),xlim=c(-8,5))#, mlab= mlabfun("RE Model", res_GI)
grid.text("GL and BC", .5, .9, gp=gpar(cex=2))
# funnel plot 
funnel(trimfill(res_GL), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(res_GL)
plot(influence(res_GL), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(res_GL)

# Fiber intake
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_res_fiber.pdf")
forest(res_fiber, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)), showweights=TRUE,header="First author(s) and Year",order=order(study_effects_fiber$yi),xlim=c(-8,5))#mlab= mlabfun("FE Model", fes_fiber)
grid.text("Fiber intake and BC", .5, .9, gp=gpar(cex=2))
# funnel plot
funnel(trimfill(res_fiber), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(res_fiber)
plot(influence(res_fiber), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(res_fiber)

# FE-significants
# Glycemic index
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_fes_GI.pdf")
forest(fes_GI, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)),showweights=TRUE,header="First author(s) and Year",order=order(study_effects_GI$yi),xlim=c(-8,5))#, mlab= mlabfun("FE Model", res_GI)
grid.text("GI and BC", .5, .9, gp=gpar(cex=2))
# funnel plot
funnel(trimfill(fes_GI), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(fes_GI)
plot(influence(fes_GI), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(fes_GI)

# Glycemic load
# pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_fes_GL.pdf")
# forest(fes_GL, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)), showweights=TRUE,header="First author(s) and Year",order=order(study_effects_GL$yi),xlim=c(-8,5))#, mlab= mlabfun("FE Model", res_GI)
# grid.text("GL and BC", .5, .9, gp=gpar(cex=2))
# # funnel plot
# funnel(trimfill(fes_GL), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# # radial plot
# radial(fes_GL)
# plot(influence(fes_GL), cex=0.8, las=1)
# dev.off()
# # regression test for funnel plot asymmetry
# regtest(fes_GL)

# Sugary drinks
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_fes_sd.pdf")
forest(fes_sd, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)), showweights=TRUE,header="First author(s) and Year",order=order(study_effects_sd$yi),xlim=c(-8,5))#, mlab= mlabfun("FE Model", res_GI)
grid.text("Sugary drinks and BC", .5, .9, gp=gpar(cex=2))
# funnel plot
funnel(trimfill(fes_sd), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(fes_sd)
plot(influence(fes_sd), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(fes_sd)

# Fiber intake
pdf(file= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results/forestplot_fes_fiber.pdf")
forest(fes_fiber, atransf=exp, at=log(c(.25, 1, 2, 2.5,4)), showweights=TRUE,header="First author(s) and Year",order=order(study_effects_fiber$yi),xlim=c(-8,5), mlab= mlabfun("FE Model", fes_fiber))
grid.text("Fiber intakes and BC", .5, .9, gp=gpar(cex=2))
# funnel plot
funnel(trimfill(fes_fiber), las=1, ylim=c(0,.8), digits=list(1L,1), legend=TRUE)
# radial plot
radial(fes_fiber)
plot(influence(fes_fiber), cex=0.8, las=1)
dev.off()
# regression test for funnel plot asymmetry
regtest(fes_fiber)

# Mixed-effects meta-regression model with moderator: Country
## Note: study_effects_GL,study_effects_sd: Not enough number of countries
# Recode country for GI
study_effects_GI$Country_recoded <- with(study_effects_GI, ifelse(Country == "USA" |  Country == "Canada", "Americas", "Europe&Asia"))
res_me_GI_country <- rma(yi, sei=sei, mods = ~ Country_recoded, data=study_effects_GI, test="knha")
res_me_GI_country
# bubble plot (with points outside of the prediction interval labeled)
regplot(res_me_GI_country,pi=TRUE, xlab="Country",label="piout", labsize=0.9, bty="l", las=1, transf=exp, refline=1, legend=TRUE)

# Mixed-effects meta-regression model with moderator:Adequacy adjustment
res_me_GI_adequ <- rma(yi, sei=sei, mods = ~ Adequacy_adjustment, data=study_effects_GI, test="knha")
res_me_GI_adequ
# bubble plot (with points outside of the prediction interval labeled)
regplot(res_me_GI_adequ,pi=TRUE, xlab="Adequacy of adjustment",label="piout", labsize=0.9, bty="l", las=1, transf=exp, refline=1, legend=TRUE)

####################### ####################### ####################### ####################### 
####################### PREMENOPUASE ########################################################
####################### ####################### ####################### ####################### 
#2.A filter for premenopause 
study_premeno_bc<- dplyr::filter(study_bmi_bc, Menopause_status == "Premenopause")

# Reverse log estimates
yi <- log(study_premeno_bc$Estimate) 
sei <- (log(study_premeno_bc$LCI_95) - log(study_premeno_bc$UCI_95)) / (2*1.96)
SE<- abs((study_premeno_bc$LCI_95 - study_premeno_bc$UCI_95) / 3.92)

# store yi and sei in data set 
study_premeno_bc$yi<- yi
study_premeno_bc$sei<- sei
study_premeno_bc$SE<- SE

# filter 0 variance in yi and sei
study_premeno_bc<- dplyr::filter(study_premeno_bc, yi != 0 & sei!= 0 & SE!= 0)
dim(study_premeno_bc) #29

# filter only highest quantile
study_premeno_bc<- dplyr::filter(study_premeno_bc, Highest_quantile == 'Yes') #study_premeno_bc[study_premeno_bc$Highest_quantile == 'Yes',]

#2A. Filter according to study characteristics: BC irrespective of menopause status
# remove
# 2A.1 Exposures as Sugary drinks: Total sugary drinks = SSB and fruit juice
study_effects_SSB_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)", "Soft drinks", "Sugary drinks", "Industrially produced juices"))
str(study_effects_SSB_A)
# 2A.2 include
dput(study_effects_SSB_A$Unique_id)
include_SSB_A<- c(129, 795, 823)
study_effects_SSB_A<- study_effects_SSB_A[study_effects_SSB_A$Unique_id %in% include_SSB_A,]
str(study_effects_SSB_A) #3

# 2A.3 Exposures as glycemic index
study_effects_GI_A <- dplyr::filter(study_premeno_bc, Exposure_recoded == "Glycemic index")
str(study_effects_GI_A)
# 2A.4 include
dput(study_effects_GI_A$Unique_id)
include_GI_A<- c(284, 409, 427, 733, 853, 983, 1013, 1101)
study_effects_GI_A<- study_effects_GI_A[study_effects_GI_A$Unique_id %in% include_GI_A,]
str(study_effects_GI_A) #8

# 2A.5 Exposures as glycemic load
study_effects_GL_A <- dplyr::filter(study_premeno_bc, Exposure_recoded == "Glycemic load")
str(study_effects_GL_A)
# 2A.6 include
dput(study_effects_GL_A$Unique_id)
include_GL_A<- c(299, 399, 432, 725, 868, 998, 1028, 1116)
study_effects_GL_A<- study_effects_GL_A[study_effects_GL_A$Unique_id %in% include_GL_A,]
str(study_effects_GL_A) #8

# 2A.7 Exposures as Total sugar
# study_effects_ts_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Total sugars","Total fructose"))
# str(study_effects_ts_A)
# # 2A.8 include
# dput(study_effects_ts_A$Unique_id)
# include_ts_A<- 
# study_effects_ts_A<- study_effects_ts_A[study_effects_ts_A$Unique_id %in% include_ts_A,]
# str(study_effects_ts_A) #

# 2A.9 Exposures as sugary food
# study_effects_sf_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Sugary foods","Sweet foods"))
# str(study_effects_sf_A)
# # 2A.10 include
# dput(study_effects_sf_A$Unique_id)
# include_sf_A<- 
# study_effects_sf_A<- study_effects_sf_A[study_effects_sf_A$Unique_id %in% include_sf_A,]
# str(study_effects_sf_A) #

# 2A.11 Exposures as Added sugar
# study_effects_as_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Added sugar"))
# str(study_effects_as_A)
# # 2A.12 include
# dput(study_effects_as_A$Unique_id)
# include_as_A<- 
# study_effects_as_A<- study_effects_as_A[study_effects_as_A$Unique_id %in% include_as_A,]
# str(study_effects_as_A) #

# 2A.13 Exposures as common sugars: glucose, fructose, sucrose
# study_effects_cs_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Glucose","Fructose","Sucrose"))
# str(study_effects_cs_A)
# # 2A.14 include
# dput(study_effects_cs_A$Unique_id)
# include_cs_A<- 
# study_effects_cs_A<- study_effects_cs_A[study_effects_cs_A$Unique_id %in% include_cs_A,]
# str(study_effects_cs_A) #

# 2.15 Exposures as sugary drinks: SSB and fruit juice
# study_effects_sd_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Sugary drinks"))
# str(study_effects_sd_A)
# # 2A.16 include
# dput(study_effects_sd_A$Unique_id)
# include_sd_A<- 
# study_effects_sd_A<- study_effects_sd_A[study_effects_sd_A$Unique_id %in% include_sd_A,]
# str(study_effects_sd_A) #1

# 2.17 Exposures as fruit juice
study_effects_fj_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Fruit juice"))
str(study_effects_fj_A)
# 2A.18 include
dput(study_effects_fj_A$Unique_id)
include_fj_A<- c(133, 259)
study_effects_fj_A<- study_effects_fj_A[study_effects_fj_A$Unique_id %in% include_fj_A,]
str(study_effects_fj_A) #2

# 2.19 Exposures as fruit juice
study_effects_fiber_A <- dplyr::filter(study_premeno_bc, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
str(study_effects_fiber_A)
# 2A.29 include
dput(study_effects_fiber_A$Unique_id)
include_fiber_A<- c(344,  697, 1135)#450,
study_effects_fiber_A<- study_effects_fiber_A[study_effects_fiber_A$Unique_id %in% include_fiber_A,]
str(study_effects_fiber_A) #2

#3A. random-effects model (using log risk ratios and variances as input)
res_SSB_A <- rma(yi, sei=sei, data=study_effects_SSB_A, method="REML",slab=Author_year,measure = "RR")
summary(res_SSB_A) # 0.46
res_GI_A <- rma(yi, sei=sei, data=study_effects_GI_A, method="REML",slab=Author_year,measure = "RR")
summary(res_GI_A) # 0.23
res_GL_A <- rma(yi, sei=sei, data=study_effects_GL_A, method="REML",slab=Author_year,measure = "RR")
summary(res_GL_A)# 0.98
# res_ts_A <- rma(yi, sei=sei, data=study_effects_ts_A, method="REML",slab=Author_year,measure = "RR")
# res_sf_A <- rma(yi, sei=sei, data=study_effects_sf_A, method="REML",slab=Author_year,measure = "RR")
# res_as_A <- rma(yi, sei=sei, data=study_effects_as_A, method="REML",slab=Author_year,measure = "RR")
# res_cs_A <- rma(yi, sei=sei, data=study_effects_cs_A, method="REML",slab=Author_year,measure = "RR")
# res_sd_A <- rma(yi, sei=sei, data=study_effects_sd_A, method="REML",slab=Author_year,measure = "RR")
res_fj_A <- rma(yi, sei=sei, data=study_effects_fj_A, method="REML",slab=Author_year,measure = "RR")
summary(res_fj_A) # 0.79
res_fiber_A <- rma(yi, sei=sei, data=study_effects_fiber_A, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_A) #0.0008 
#
reporter(res_fiber_A, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")
pred_res_fiber_A<- predict(res_fiber_A, transf=exp, digits=2); pred_res_fiber_A
#
predict(res_GI_A, transf = exp)
predict(res_GL_A, transf = exp)
predict(res_fiber_A, transf = exp)
predict(res_SSB_A, transf = exp)

####################### ####################### ####################### ####################### 
####################### POSMENOPUASE ########################################################
####################### ####################### ####################### ####################### 
#2B. filter for posmenopause 
study_posmeno_bc<- dplyr::filter(study_bmi_bc, Menopause_status == "Posmenopause")

# Reverse log estimates
yi <- log(study_posmeno_bc$Estimate) 
sei <- (log(study_posmeno_bc$LCI_95) - log(study_posmeno_bc$UCI_95)) / (2*1.96)
SE<- abs((study_posmeno_bc$LCI_95 - study_posmeno_bc$UCI_95) / 3.92)

# store yi and sei in data set 
study_posmeno_bc$yi<- yi
study_posmeno_bc$sei<- sei
study_posmeno_bc$SE<- SE

# filter 0 variance in yi and sei
study_posmeno_bc<- dplyr::filter(study_posmeno_bc, yi != 0 & sei!= 0 & SE!= 0)
dim(study_posmeno_bc) #132

# filter only highest quantile
study_posmeno_bc<- dplyr::filter(study_posmeno_bc, Highest_quantile == 'Yes') #study_posmeno_bc[study_posmeno_bc$Highest_quantile == 'Yes',]

#2B. Filter according to study characteristics: BC irrespective of menopause status
# remove
# 2B.1 Exposures as Sugary drinks: Total sugary drinks = SSB and fruit juice
study_effects_SSB_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)", "Soft drinks", "Sugary drinks", "Industrially produced juices"))
str(study_effects_SSB_B)
# 2B.2 include
dput(study_effects_SSB_B$Unique_id)
include_SSB_B<- c(141, 419, 793, 827)
study_effects_SSB_B<- study_effects_SSB_B[study_effects_SSB_B$Unique_id %in% include_SSB_B,]
str(study_effects_SSB_B) #4

# 2B.3 Exposures as glycemic index
study_effects_GI_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded == "Glycemic index")
str(study_effects_GI_B)
# 2B.4 include
dput(study_effects_GI_B$Unique_id)
include_GI_B<- c(289, 414, 437, 737, 858, 988, 1018, 1106)
study_effects_GI_B<- study_effects_GI_B[study_effects_GI_B$Unique_id %in% include_GI_B,]
str(study_effects_GI_B) #8

# 2B.5 Exposures as glycemic load
study_effects_GL_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded == "Glycemic load")
str(study_effects_GL_B)
# 2B.6 include
dput(study_effects_GL_B$Unique_id)
include_GL_B<- c(304, 404, 442, 729, 873, 1003, 1033, 1121)
study_effects_GL_B<- study_effects_GL_B[study_effects_GL_B$Unique_id %in% include_GL_B,]
str(study_effects_GL_B) #8

# 2B.7 Exposures as Total sugar
# study_effects_ts_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Total sugars","Total fructose"))
# str(study_effects_ts_B)
# # 2B.8 include
# dput(study_effects_ts_B$Unique_id)
# include_ts_B<- 
# study_effects_ts_B<- study_effects_ts_B[study_effects_ts_B$Unique_id %in% include_ts_B,]
# str(study_effects_ts_B) #

# 2B.9 Exposures as sugary food
# study_effects_sf_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Sugary foods","Sweet foods"))
# str(study_effects_sf_B)
# # 2B.10 include
# dput(study_effects_sf_B$Unique_id)
# include_sf_B<- 
# study_effects_sf_B<- study_effects_sf_B[study_effects_sf_B$Unique_id %in% include_sf_B,]
# str(study_effects_sf_B) #

# 2B.11 Exposures as Added sugar
# study_effects_as_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Added sugar"))
# str(study_effects_as_B)
# # 2B.12 include
# dput(study_effects_as_B$Unique_id)
# include_as_B<- 
# study_effects_as_B<- study_effects_as_B[study_effects_as_B$Unique_id %in% include_as_B,]
# str(study_effects_as_B) #0

# 2B.13 Exposures as common sugars: glucose, fructose, sucrose
# study_effects_cs_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Glucose","Fructose","Sucrose"))
# str(study_effects_cs_B)
# # 2B.14 include
# dput(study_effects_cs_B$Unique_id)
# include_cs_B<- 
# study_effects_cs_B<- study_effects_cs_B[study_effects_cs_B$Unique_id %in% include_cs_B,]
# str(study_effects_cs_B) #0

# 2.15 Exposures as sugary drinks: SSB and fruit juice
# study_effects_sd_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Sugary drinks"))
# str(study_effects_sd_B)
# # 2B.16 include
# dput(study_effects_sd_B$Unique_id)
# include_sd_B<- 
# study_effects_sd_B<- study_effects_sd_B[study_effects_sd_B$Unique_id %in% include_sd_B,]
# str(study_effects_sd_B) #1

# 2.17 Exposures as fruit juice
study_effects_fj_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Fruit juice"))
str(study_effects_fj_B)
# 2B.18 include
dput(study_effects_fj_B$Unique_id)
include_fj_B<- c(145, 264, 1091)
study_effects_fj_B<- study_effects_fj_B[study_effects_fj_B$Unique_id %in% include_fj_B,]
str(study_effects_fj_B) #3

# 2.19 Exposures as fruit juice
study_effects_fiber_B <- dplyr::filter(study_posmeno_bc, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
str(study_effects_fiber_B)
# 2B.20 include
dput(study_effects_fiber_B$Unique_id)
include_fiber_B<- c(349, 381, 684, 701, 766, 1055, 1140)
study_effects_fiber_B<- study_effects_fiber_B[study_effects_fiber_B$Unique_id %in% include_fiber_B,]
str(study_effects_fiber_B) #7

#3B. random-effects model (using log risk ratios and variances as input)
res_SSB_B <- rma(yi, sei=sei, data=study_effects_SSB_B, method="REML",slab=Author_year,measure = "RR")
summary(res_SSB_B) # 0.055
reporter(res_SSB_B, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")
pred_res_SSB_B<- predict(res_SSB_B, transf=exp, digits=2); pred_res_SSB_B
#
res_GI_B <- rma(yi, sei=sei, data=study_effects_GI_B, method="REML",slab=Author_year,measure = "RR")
summary(res_GI_B) # 0.0077
reporter(res_GI_B, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")
pred_res_GI_B<- predict(res_GI_B, transf=exp, digits=2); pred_res_GI_B
#
res_GL_B <- rma(yi, sei=sei, data=study_effects_GL_B, method="REML",slab=Author_year,measure = "RR")
summary(res_GL_B)# 0.06
# res_ts_B <- rma(yi, sei=sei, data=study_effects_ts_B, method="REML",slab=Author_year,measure = "RR")
# res_sf_B <- rma(yi, sei=sei, data=study_effects_sf_B, method="REML",slab=Author_year,measure = "RR")
# res_as_B <- rma(yi, sei=sei, data=study_effects_as_B, method="REML",slab=Author_year,measure = "RR")
# res_cs_B <- rma(yi, sei=sei, data=study_effects_cs_B, method="REML",slab=Author_year,measure = "RR")
# res_sd_B <- rma(yi, sei=sei, data=study_effects_sd_B, method="REML",slab=Author_year,measure = "RR")
res_fj_B <- rma(yi, sei=sei, data=study_effects_fj_B, method="REML",slab=Author_year,measure = "RR")
summary(res_fj_B) # 0.19
res_fiber_B <- rma(yi, sei=sei, data=study_effects_fiber_B, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_B) # 0.018
reporter(res_fiber_B, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")
pred_res_fiber_B<- predict(res_fiber_B, transf=exp, digits=2); pred_res_fiber_B
#
predict(res_GI_B, transf = exp)
predict(res_GL_B, transf = exp)
predict(res_fiber_B, transf = exp)
predict(res_SSB_B, transf = exp)
####################### ####################### ####################### ####################### 
####################### hormone_receptors ER+/PR+ or ER+ (only) or PR+ (only) ########################################################
####################### ####################### ####################### ####################### 
#2C. filter for posmenopause 
study_hr_bc<- dplyr::filter(study_tool_bc, hormone_receptors == "ER+/PR+"  | hormone_receptors == "ER+" | hormone_receptors == "PR+")

# Reverse log estimates
yi <- log(study_hr_bc$Estimate) 
sei <- (log(study_hr_bc$LCI_95) - log(study_hr_bc$UCI_95)) / (2*1.96)
SE<- abs((study_hr_bc$LCI_95 - study_hr_bc$UCI_95) / 3.92)

# store yi and sei in data set 
study_hr_bc$yi<- yi
study_hr_bc$sei<- sei
study_hr_bc$SE<- SE

# filter 0 variance in yi and sei
study_hr_bc<- dplyr::filter(study_hr_bc, yi != 0 & sei!= 0 & SE!= 0)
dim(study_hr_bc) #68

# filter only highest quantile
study_hr_bc<- dplyr::filter(study_hr_bc, Highest_quantile == 'Yes') #study_hr_bc[study_hr_bc$Highest_quantile == 'Yes',]

#2C. Filter according to study characteristics: BC irrespective of menopause status
# remove
# 2C.1 Exposures as Sugary drinks: Total sugary drinks = SSB and fruit juice
# study_effects_SSB_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)", "Soft drinks", "Sugary drinks", "Industrially produced juices"))
# str(study_effects_SSB_C)
# # 2C.2 include
# dput(study_effects_SSB_C$Unique_id)
# include_SSB_C<- c(141, 419, 793, 827)
# study_effects_SSB_C<- study_effects_SSB_C[study_effects_SSB_C$Unique_id %in% include_SSB_C,]
# str(study_effects_SSB_C) #1

# 2C.3 Exposures as glycemic index
study_effects_GI_C <- dplyr::filter(study_hr_bc, Exposure_recoded == "Glycemic index")
str(study_effects_GI_C)
# 2C.4 include
dput(study_effects_GI_C$Unique_id)
include_GI_C<- c(309, 532, 595, 953) #548
study_effects_GI_C<- study_effects_GI_C[study_effects_GI_C$Unique_id %in% include_GI_C,]
str(study_effects_GI_C) #4

# 2C.5 Exposures as glycemic load
study_effects_GL_C <- dplyr::filter(study_hr_bc, Exposure_recoded == "Glycemic load")
str(study_effects_GL_C)
# 2C.6 include
dput(study_effects_GL_C$Unique_id)
include_GL_C<- c(536, 615, 933) #552
study_effects_GL_C<- study_effects_GL_C[study_effects_GL_C$Unique_id %in% include_GL_C,]
str(study_effects_GL_C) #4

# 2C.7 Exposures as Total sugar
# study_effects_ts_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Total sugars","Total fructose"))
# str(study_effects_ts_C)
# # 2C.8 include
# dput(study_effects_ts_C$Unique_id)
# include_ts_C<- 
# study_effects_ts_C<- study_effects_ts_C[study_effects_ts_C$Unique_id %in% include_ts_C,]
# str(study_effects_ts_C) #

# 2C.9 Exposures as sugary food
# study_effects_sf_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugary foods","Sweet foods"))
# str(study_effects_sf_C)
# # 2C.10 include
# dput(study_effects_sf_C$Unique_id)
# include_sf_C<- 
# study_effects_sf_C<- study_effects_sf_C[study_effects_sf_C$Unique_id %in% include_sf_C,]
# str(study_effects_sf_C) #

# 2C.11 Exposures as Added sugar
# study_effects_as_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Added sugar"))
# str(study_effects_as_C)
# # 2C.12 include
# dput(study_effects_as_C$Unique_id)
# include_as_C<- 
# study_effects_as_C<- study_effects_as_C[study_effects_as_C$Unique_id %in% include_as_C,]
# str(study_effects_as_C) #0

# 2C.13 Exposures as common sugars: glucose, fructose, sucrose
# study_effects_cs_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Glucose","Fructose","Sucrose"))
# str(study_effects_cs_C)
# # 2C.14 include
# dput(study_effects_cs_C$Unique_id)
# include_cs_C<- 
# study_effects_cs_C<- study_effects_cs_C[study_effects_cs_C$Unique_id %in% include_cs_C,]
# str(study_effects_cs_C) #0

# 2.15 Exposures as sugary drinks: SSB and fruit juice
# study_effects_sd_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugary drinks"))
# str(study_effects_sd_C)
# # 2C.16 include
# dput(study_effects_sd_C$Unique_id)
# include_sd_C<- 
# study_effects_sd_C<- study_effects_sd_C[study_effects_sd_C$Unique_id %in% include_sd_C,]
# str(study_effects_sd_C) #0

# 2.17 Exposures as fruit juice
# study_effects_fj_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Fruit juice"))
# str(study_effects_fj_C)
# # 2C.18 include
# dput(study_effects_fj_C$Unique_id)
# include_fj_C<- c(145, 264, 1091)
# study_effects_fj_C<- study_effects_fj_C[study_effects_fj_C$Unique_id %in% include_fj_C,]
# str(study_effects_fj_C) #0

# 2.19 Exposures as fruit juice
study_effects_fiber_C <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
str(study_effects_fiber_C)
# 2C.20 include
dput(study_effects_fiber_C$Unique_id)
include_fiber_C<- c(495, 705, 771, 1060)
study_effects_fiber_C<- study_effects_fiber_C[study_effects_fiber_C$Unique_id %in% include_fiber_C,]
str(study_effects_fiber_C) #

#3B. random-effects model (using log risk ratios and variances as input)
#res_SSB_C <- rma(yi, sei=sei, data=study_effects_SSB_C, method="REML",slab=Author_year,measure = "RR")
res_GI_C <- rma(yi, sei=sei, data=study_effects_GI_C, method="REML",slab=Author_year,measure = "RR")
summary(res_GI_C) # 0.76
res_GL_C <- rma(yi, sei=sei, data=study_effects_GL_C, method="REML",slab=Author_year,measure = "RR")
summary(res_GL_C)# 0.07
# res_ts_C <- rma(yi, sei=sei, data=study_effects_ts_C, method="REML",slab=Author_year,measure = "RR")
# res_sf_C <- rma(yi, sei=sei, data=study_effects_sf_C, method="REML",slab=Author_year,measure = "RR")
# res_as_C <- rma(yi, sei=sei, data=study_effects_as_C, method="REML",slab=Author_year,measure = "RR")
# res_cs_C <- rma(yi, sei=sei, data=study_effects_cs_C, method="REML",slab=Author_year,measure = "RR")
# res_sd_C <- rma(yi, sei=sei, data=study_effects_sd_C, method="REML",slab=Author_year,measure = "RR")
# res_fj_C <- rma(yi, sei=sei, data=study_effects_fj_C, method="REML",slab=Author_year,measure = "RR")
res_fiber_C <- rma(yi, sei=sei, data=study_effects_fiber_C, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_C)#0.52

####################### ####################### ####################### ####################### 
####################### hormone_receptors ER-/PR- or ER- (only) or PR- (only) ########################################################
####################### ####################### ####################### ####################### 
#2C. filter for posmenopause 
study_hr_bc<- dplyr::filter(study_tool_bc, hormone_receptors == "ER-/PR-"  | hormone_receptors == "ER-" | hormone_receptors == "PR-")

# Reverse log estimates
yi <- log(study_hr_bc$Estimate) 
sei <- (log(study_hr_bc$LCI_95) - log(study_hr_bc$UCI_95)) / (2*1.96)
SE<- abs((study_hr_bc$LCI_95 - study_hr_bc$UCI_95) / 3.92)

# store yi and sei in data set 
study_hr_bc$yi<- yi
study_hr_bc$sei<- sei
study_hr_bc$SE<- SE

# filter 0 variance in yi and sei
study_hr_bc<- dplyr::filter(study_hr_bc, yi != 0 & sei!= 0 & SE!= 0)
dim(study_hr_bc) #68

# filter only highest quantile
study_hr_bc<- dplyr::filter(study_hr_bc, Highest_quantile == 'Yes') #study_hr_bc[study_hr_bc$Highest_quantile == 'Yes',]

#2D. Filter according to study characteristics: BC irrespective of menopause status
# remove
# 2D.1 Exposures as Sugary drinks: Total sugary drinks = SSB and fruit juice
# study_effects_SSB_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugar-sweetened beverages (SSB)", "Soft drinks", "Sugary drinks", "Industrially produced juices"))
# str(study_effects_SSB_D)
# # 2D.2 include
# dput(study_effects_SSB_D$Unique_id)
# include_SSB_D<- c(141, 419, 793, 827)
# study_effects_SSB_D<- study_effects_SSB_D[study_effects_SSB_D$Unique_id %in% include_SSB_D,]
# str(study_effects_SSB_D) #1

# 2D.3 Exposures as glycemic index
study_effects_GI_D <- dplyr::filter(study_hr_bc, Exposure_recoded == "Glycemic index")
str(study_effects_GI_D)
# 2D.4 include
dput(study_effects_GI_D$Unique_id)
include_GI_D<- c(314, 556, 605, 888, 963) #908, 540, 878, 898
study_effects_GI_D<- study_effects_GI_D[study_effects_GI_D$Unique_id %in% include_GI_D,]
str(study_effects_GI_D) #5

# 2D.5 Exposures as glycemic load
study_effects_GL_D <- dplyr::filter(study_hr_bc, Exposure_recoded == "Glycemic load")
str(study_effects_GL_D)
# 2D.6 include
dput(study_effects_GL_D$Unique_id)
include_GL_D<- c(544, 625, 893, 943) #903, 913, 560, 883
study_effects_GL_D<- study_effects_GL_D[study_effects_GL_D$Unique_id %in% include_GL_D,]
str(study_effects_GL_D) #4

# 2D.7 Exposures as Total sugar
# study_effects_ts_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Total sugars","Total fructose"))
# str(study_effects_ts_D)
# # 2D.8 include
# dput(study_effects_ts_D$Unique_id)
# include_ts_D<- 
# study_effects_ts_D<- study_effects_ts_D[study_effects_ts_D$Unique_id %in% include_ts_D,]
# str(study_effects_ts_D) #

# 2D.9 Exposures as sugary food
# study_effects_sf_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugary foods","Sweet foods"))
# str(study_effects_sf_D)
# # 2D.10 include
# dput(study_effects_sf_D$Unique_id)
# include_sf_D<- 
# study_effects_sf_D<- study_effects_sf_D[study_effects_sf_D$Unique_id %in% include_sf_D,]
# str(study_effects_sf_D) #

# 2D.11 Exposures as Added sugar
# study_effects_as_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Added sugar"))
# str(study_effects_as_D)
# # 2D.12 include
# dput(study_effects_as_D$Unique_id)
# include_as_D<- 
# study_effects_as_D<- study_effects_as_D[study_effects_as_D$Unique_id %in% include_as_D,]
# str(study_effects_as_D) #0

# 2D.13 Exposures as common sugars: glucose, fructose, sucrose
# study_effects_cs_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Glucose","Fructose","Sucrose"))
# str(study_effects_cs_D)
# # 2D.14 include
# dput(study_effects_cs_D$Unique_id)
# include_cs_D<- 
# study_effects_cs_D<- study_effects_cs_D[study_effects_cs_D$Unique_id %in% include_cs_D,]
# str(study_effects_cs_D) #0

# 2.15 Exposures as sugary drinks: SSB and fruit juice
# study_effects_sd_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Sugary drinks"))
# str(study_effects_sd_D)
# # 2D.16 include
# dput(study_effects_sd_D$Unique_id)
# include_sd_D<- 
# study_effects_sd_D<- study_effects_sd_D[study_effects_sd_D$Unique_id %in% include_sd_D,]
# str(study_effects_sd_D) #0

# 2.17 Exposures as fruit juice
# study_effects_fj_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Fruit juice"))
# str(study_effects_fj_D)
# # 2D.18 include
# dput(study_effects_fj_D$Unique_id)
# include_fj_D<- c(145, 264, 1091)
# study_effects_fj_D<- study_effects_fj_D[study_effects_fj_D$Unique_id %in% include_fj_D,]
# str(study_effects_fj_D) #0

# 2.19 Exposures as fruit juice
study_effects_fiber_D <- dplyr::filter(study_hr_bc, Exposure_recoded %in% c("Total fiber", "Fiber","Total fiber intake","Fiber intake"))
str(study_effects_fiber_D)
# 2D.20 include
dput(study_effects_fiber_D$Unique_id)
include_fiber_D<- c(500, 565, 709, 781, 1070)
study_effects_fiber_D<- study_effects_fiber_D[study_effects_fiber_D$Unique_id %in% include_fiber_D,]
str(study_effects_fiber_D) #5

#3B. random-effects model (using log risk ratios and variances as input)
#res_SSB_D <- rma(yi, sei=sei, data=study_effects_SSB_D, method="REML",slab=Author_year,measure = "RR")
res_GI_D <- rma(yi, sei=sei, data=study_effects_GI_D, method="REML",slab=Author_year,measure = "RR")
summary(res_GI_D) # 0.58
#
res_GL_D <- rma(yi, sei=sei, data=study_effects_GL_D, method="REML",slab=Author_year,measure = "RR")
summary(res_GL_D)# 0.0043
reporter(res_GL_D, dir= "/Users/hugopomaresmillan/Desktop/Skeleton_projects/nutrient_cancer/syst_rev/Results",format="word_document")
pred_res_GL_D<- predict(res_GL_D, transf=exp, digits=2); pred_res_GL_D
# res_ts_D <- rma(yi, sei=sei, data=study_effects_ts_D, method="REML",slab=Author_year,measure = "RR")
# res_sf_D <- rma(yi, sei=sei, data=study_effects_sf_D, method="REML",slab=Author_year,measure = "RR")
# res_as_D <- rma(yi, sei=sei, data=study_effects_as_D, method="REML",slab=Author_year,measure = "RR")
# res_cs_D <- rma(yi, sei=sei, data=study_effects_cs_D, method="REML",slab=Author_year,measure = "RR")
# res_sd_D <- rma(yi, sei=sei, data=study_effects_sd_D, method="REML",slab=Author_year,measure = "RR")
# res_fj_D <- rma(yi, sei=sei, data=study_effects_fj_D, method="REML",slab=Author_year,measure = "RR")
res_fiber_D <- rma(yi, sei=sei, data=study_effects_fiber_D, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_D)# 0.41

###############################################################################
# ONLY POSMENOPAUSE HAS RECEPTORS REPORTED for Fiber exposure
###############################################################################
#2E. filter for posmenopause 
study_posmeno_bc<- dplyr::filter(study_tool_bc, Menopause_status == "Posmenopause")
study_posmeno_hr_mebc<- dplyr::filter(study_posmeno_bc, hormone_receptors == "ER-/PR-"  | hormone_receptors == "ER-" | hormone_receptors == "PR-")

# Reverse log estimates
yi <- log(study_posmeno_hr_mebc$Estimate) 
sei <- (log(study_posmeno_hr_mebc$LCI_95) - log(study_posmeno_hr_mebc$UCI_95)) / (2*1.96)
SE<- abs((study_posmeno_hr_mebc$LCI_95 - study_posmeno_hr_mebc$UCI_95) / 3.92)

# store yi and sei in data set 
study_posmeno_hr_mebc$yi<- yi
study_posmeno_hr_mebc$sei<- sei
study_posmeno_hr_mebc$SE<- SE

# filter 0 variance in yi and sei
study_posmeno_hr_mebc<- dplyr::filter(study_posmeno_hr_mebc, yi != 0 & sei!= 0 & SE!= 0)
dim(study_posmeno_hr_mebc) #32

# filter only highest quantile
study_posmeno_hr_mebc<- dplyr::filter(study_posmeno_hr_mebc, Highest_quantile == 'Yes') #study_posmeno_hr_mebc[study_posmeno_hr_mebc$Highest_quantile == 'Yes',]

# 2E.3 Exposures as fiber
# 2E.4 include
dput(study_posmeno_hr_mebc$Unique_id)
include_fiber_E<- c(781, 1070)
study_effects_fiber_E<- study_posmeno_hr_mebc[study_posmeno_hr_mebc$Unique_id %in% include_fiber_E,]
str(study_effects_fiber_E) #2

#3.
res_fiber_E <- rma(yi, sei=sei, data=study_effects_fiber_E, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_E) # 0.14

###############################################################################
study_posmeno_hr_mebc2<- dplyr::filter(study_posmeno_bc, hormone_receptors == "ER+/PR+"  | hormone_receptors == "ER+" | hormone_receptors == "PR+")

# Reverse log estimates
yi <- log(study_posmeno_hr_mebc2$Estimate) 
sei <- (log(study_posmeno_hr_mebc2$LCI_95) - log(study_posmeno_hr_mebc2$UCI_95)) / (2*1.96)
SE<- abs((study_posmeno_hr_mebc2$LCI_95 - study_posmeno_hr_mebc2$UCI_95) / 3.92)

# store yi and sei in data set 
study_posmeno_hr_mebc2$yi<- yi
study_posmeno_hr_mebc2$sei<- sei
study_posmeno_hr_mebc2$SE<- SE

# filter 0 variance in yi and sei
study_posmeno_hr_mebc2<- dplyr::filter(study_posmeno_hr_mebc2, yi != 0 & sei!= 0 & SE!= 0)
dim(study_posmeno_hr_mebc2) #16

# filter only highest quantile
study_posmeno_hr_mebc2<- dplyr::filter(study_posmeno_hr_mebc2, Highest_quantile == 'Yes') #study_posmeno_hr_mebc2[study_posmeno_hr_mebc2$Highest_quantile == 'Yes',]

# 2E.3 Exposures as fiber
# 2E.4 include
dput(study_posmeno_hr_mebc2$Unique_id)
include_fiber_E<- c(771, 1060)
study_effects_fiber_E<- study_posmeno_hr_mebc2[study_posmeno_hr_mebc2$Unique_id %in% include_fiber_E,]
str(study_effects_fiber_E) #2

#3.
res_fiber_E <- rma(yi, sei=sei, data=study_effects_fiber_E, method="REML",slab=Author_year,measure = "RR")
summary(res_fiber_E) # 0.14


###############################################################################
### BAYESIAN #####
###############################################################################
# library(brms)
# brm_out <- brm(
#   yi | se(SE) ~ 1 + (1 | Author_year),
#   data = study_effects_fiber_E,
#   cores = 4,
#   file = "metaanalysismodel"
# )
# #
# hypothesis(brm_out, "Intercept > 0.2")
############################## END ##############################################