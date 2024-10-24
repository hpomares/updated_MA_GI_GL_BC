# Meta analysis of observational studies (Sugar intake, added sugar, SSB,etc associated with BC)
# Dose-response MA
# examples: https://rpubs.com/alecri/glst | https://rpubs.com/alecri/berlin | https://rpubs.com/alecri/qinliu
# Load the package:
library(readxl)
library(DescTools)
library(metafor)
library(dplyr)
library(rio)
library(grid)
library(dosresmeta)

setwd("/Users/hugopomaresmillan/syst_rev/Results")

#1.  read file
study_tool <- read_excel("~/study_tool.xlsx",sheet = "Work_sheet_v4")
str(study_tool)

# sort by year
#study_tool <- study_tool[order(-study_tool$Year),]

# filter case-control designs
study_tool_cohort<- dplyr::filter(study_tool, Study_design == "Cohort")
study_tool_cohort$type <- "ir"

# filter for BC
study_tool_bc<- dplyr::filter(study_tool_cohort, Outcome == "Breast cancer")

# select effect measures
study_effects<- study_tool_bc %>% select("Author_year","Name_of_study","Exposure_recoded","cases", "N","Ref_group", "Test_group","Menopause_status", "hormone_receptors", 
                                         "summary_effect", "assigned_dose", "id", "type","person_year",
                                         "Estimate","LCI_95", "UCI_95", "Country","Adequacy_adjustment")

# Reverse log estimates
yi  <- log(study_effects$Estimate) 
sei  <- (log(study_effects$LCI_95) - log(study_effects$UCI_95)) / (2*1.96)
SE<- abs((study_effects$LCI_95 - study_effects$UCI_95) / 3.92)

# store yi and sei in data set 
study_effects$yi<- yi
study_effects$sei<- sei
study_effects$SE<- SE
cols.num <- c("cases","N","assigned_dose","person_year")
study_effects[cols.num] <- sapply(study_effects[cols.num],as.numeric)
sapply(study_effects, class)
str(study_effects)

#2. Filter single study characteristics
# 2.1 Author_year
# GI: Glycemic index
# Prioritize EPIC italy rather than ORDET
# Prioritize EPIC rather than EPIC italy

study_effects_single <- dplyr::filter(study_effects, 
                                      # Author_year =="Sieri S et al. 2013" | 
                                      Author_year =="Navarro-Silvera SA et al. 2005" | 
                                        Author_year =="Romanos‐Nanclares A et al. 2021" |
                                        Author_year =="Larsson SC et al. 2009" | 
                                        Author_year =="Romieu I et al. 2012" | 
                                        Author_year =="Shikany JM et al. 2011"|
                                        Author_year =="Makarem N et al. 2017" |
                                        Author_year =="Lajous M et al. 2008" |
                                        Author_year =="Wen W et al. 2009"  |
                                        Author_year =="Holmes MD et al. 2004" |
                                        Author_year =="Higginbotham S et al. 2004" |
                                        Author_year =="Jonas CR et al. 2003" |
                                        Author_year =="Nielsen TG et al. 2005" |
                                        Author_year =="George SM et al. 2008"  # |
                                      #Author_year =="Giles GC et al. 2006"
)
# # Including case-control study design
# # Prioritize EPIC italy rather than ORDET== Sieri et al 2007
# # Nurses Health study II in Cho E et al. 2003
# study_effects_single <- dplyr::filter(study_effects, 
#                                         Author_year =="Agustin LSA et al. 2001" | 
#                                         Author_year =="Amadou A et al. 2015" | 
#                                         #Author_year =="Cho E et al. 2003" | 
#                                         Author_year =="Farvid MS et al. 2015" | 
#                                         #Author_year =="Sieri S et al. 2013" | 
#                                         Author_year =="Navarro-Silvera SA et al. 2005" | 
#                                         Author_year =="Romanos‐Nanclares A et al. 2021" |
#                                         Author_year =="Larsson SC et al. 2009" | 
#                                         Author_year =="Romieu I et al. 2012" | 
#                                         Author_year =="Shikany JM et al. 2011"|
#                                         Author_year =="Makarem N et al. 2017" |
#                                         Author_year =="Lajous M et al. 2008" |
#                                         Author_year =="Lajous M et al. 2005" |
#                                         Author_year =="Wen W et al. 2009"  |
#                                         Author_year =="Yun SH et al. 2010"  |
#                                         Author_year =="Holmes MD et al. 2004" |
#                                         Author_year =="Higginbotham S et al. 2004" |
#                                         Author_year =="Jonas CR et al. 2003" |
#                                         Author_year =="Nielsen TG et al. 2005" |
#                                         Author_year =="George SM et al. 2008"  # |
#                                       #Author_year =="Giles GC et al. 2006"
# )
# # for dose response Giles do not work

#2.2 exposure
study_effects_gi <- dplyr::filter(study_effects_single, Exposure_recoded == "Glycemic index")
dim(study_effects_gi)

#2.3 All women, not stratified by menopause or hormonal receptor status
study_effects_gi_all <- dplyr::filter(study_effects_gi, Menopause_status =="All" & hormone_receptors == "NA")
str(study_effects_gi_all) 

#2.4 fit
# fixed effects
meta_fe_gi<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                        se = SE, data = study_effects_gi_all, method = "fixed")
summary(meta_fe_gi)
predict(meta_fe_gi, delta = 1, expo = TRUE)

#2.5 fit
# Random effects
meta_re_gi<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                        se = SE, data = study_effects_gi_all, method = "reml")
summary(meta_re_gi)
predict(meta_re_gi, delta = 1, expo = TRUE)
# significant log-linear dose-response association between GI and cancer risk (p=0.04)
# identical FE ad RE

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gi, delta = 10, exp = TRUE) # associated with a significant 3%

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gi_all, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)
# 
dosex_bin <- data.frame(assigned_dose=seq(0, 100, 1))
with(predict(meta_re_gi, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model *NOT SIGNIFICANT
# square-transform the dose
quad_lin <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gi_all)
summary(quad_lin)
#plot
with(predict(quad_lin, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 3), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gi_all$assigned_dose,  probs = c(.01, .5, .9)) #c(.2, .4, .8) 
spl <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                  cases = cases, n = N, data = study_effects_gi_all)
summary(spl)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin,dataTab,exp=TRUE)
predSpl<-predict(spl,dataTab,exp=TRUE)
predSpl
round(cbind(quad_lin=predLin,spl=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl), Sigma = vcov(spl), Terms= 1:nrow(vcov(spl))) #risk of cancer may not vary in a log-linear fashion with GI
# The first of the two dose levels is excluded because it is the total raw data itself.If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3. plot RCS
xref_bin <- 0
with(predict(spl, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 10), ylab="BC relative risk", xlab= "GI, units/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 

################ PREMENOPAUSE ####################################
#2.3 All women
study_effects_gi_all_a <- dplyr::filter(study_effects_gi, Menopause_status =="Premenopause")
str(study_effects_gi_all_a) #24 

#2.4 fit
# fixed effects
meta_fe_g_a<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                         se = SE, data = study_effects_gi_all_a, method = "fixed")
summary(meta_fe_g_a)
predict(meta_fe_g_a, delta = 1, expo = TRUE)

#2.5 fit
# Random effects
meta_re_gi_a<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                          se = SE, data = study_effects_gi_all_a, method = "reml")
summary(meta_re_gi_a)
predict(meta_re_gi_a, delta = 1, expo = TRUE)
# significant log-linear dose-response association between GI and cancer risk (p=0.01)

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gi_a, delta = 10, exp = TRUE) # associated with a significant 3%

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gi_all_a, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)

# 
dosex_bin <- data.frame(assigned_dose=seq(0, 100, 1))
with(predict(meta_re_gi_a, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model *NOT SIGNIFICANT
# square-transform the dose
quad_lin_a <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gi_all_a)
summary(quad_lin_a)
#plot
with(predict(quad_lin_a, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 6), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin_a, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gi_all_a$assigned_dose,  probs = c(.35, .65, .95)) #c(.2, .4, .8) 
spl_a <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                    cases = cases, n = N, data = study_effects_gi_all_a)
summary(spl_a)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin_a,dataTab,exp=TRUE)
predSpl<-predict(spl_a,dataTab,exp=TRUE)
predSpl
round(cbind(quad_lin_a=predLin,spl_a=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl_a), Sigma = vcov(spl_a), Terms= 1:nrow(vcov(spl_a))) #risk of cancer may not vary in a log-linear fashion with GI
# The first of the two dose levels is excluded because it is the total raw data itself.If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3. plot RCS
xref_bin <- 0
with(predict(spl_a, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 5), ylab="BC relative risk", xlab= "GI, units/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 

################ POSTMENOPAUSE ####################################
#2.3 All women
study_effects_gi_all_b <- dplyr::filter(study_effects_gi, Menopause_status =="Posmenopause")
str(study_effects_gi_all_b)

# Nielsen TG et al. 2005 is in tertiles and IRR was givem
study_effects_gi_all_b <- dplyr::filter(study_effects_gi_all_b, Author_year !="Nielsen TG et al. 2005")

#2.4 fit
# fixed effects
meta_fe_g_b<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                         se = SE, data = study_effects_gi_all_b, method = "fixed")
summary(meta_fe_g_b)
predict(meta_fe_g_b, delta = 1, expo = TRUE)

#2.5 fit
# Random effects
meta_re_gi_b<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                          se = SE, data = study_effects_gi_all_b, method = "reml")
summary(meta_re_gi_b)
predict(meta_re_gi_b, delta = 1, expo = TRUE)
# significant log-linear dose-response association between GI and cancer risk (p=0.01)

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gi_b, delta = 10, exp = TRUE) # associated with a significant 3%

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gi_all_b, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)

# 
dosex_bin <- data.frame(assigned_dose=seq(0, 100, 1))
with(predict(meta_re_gi_b, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model *NOT SIGNIFICANT
# square-transform the dose
quad_lin_b <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gi_all_b)
summary(quad_lin_b)
#plot
with(predict(quad_lin_b, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 6), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin_b, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gi_all_b$assigned_dose,  probs = c(.35, .65, .95)) #c(.2, .4, .8) 
spl_b <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                    cases = cases, n = N, data = study_effects_gi_all_b)
summary(spl_b)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin_b,dataTab,exp=TRUE)
predSpl<-predict(spl_b,dataTab,exp=TRUE)
predSpl
round(cbind(quad_lin_b=predLin,spl_b=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl_b), Sigma = vcov(spl_b), Terms= 1:nrow(vcov(spl_b))) #risk of cancer may not vary in a log-linear fashion with GI
# The first of the two dose levels is excluded because it is the total raw data itself.If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3. plot RCS
xref_bin <- 0
with(predict(spl_b, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 5), ylab="BC relative risk", xlab= "GI, units/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 

####################################################
###### GLYCEMIC LOAD 
####################################################
#2. Filter single study characteristics
# 2.1 Author_year
# GL: Glycemic load
# Prioritize EPIC italy rather than ORDET
study_effects_single <- dplyr::filter(study_effects, 
                                      #Author_year =="Sieri S et al. 2013" | 
                                      Author_year =="George SM et al. 2008"  |
                                        Author_year =="Navarro-Silvera SA et al. 2005" | 
                                        Author_year =="Romanos‐Nanclares A et al. 2021" |
                                        Author_year =="Larsson SC et al. 2009" | 
                                        Author_year =="Romieu I et al. 2012" | 
                                        Author_year =="Shikany JM et al. 2011"|
                                        Author_year =="Makarem N et al. 2017" |
                                        Author_year =="Lajous M et al. 2008" |
                                        Author_year =="Wen W et al. 2009"  |
                                        Author_year =="Holmes MD et al. 2004" |
                                        Author_year =="Higginbotham S et al. 2004" |
                                        Author_year =="Jonas CR et al. 2003" |
                                        Author_year =="Nielsen TG et al. 2005" #| 
                                      #Author_year =="Giles GC et al. 2006"
)
# for dose response Giles do not work
# study_effects_single <- dplyr::filter(study_effects,
#                                         Author_year =="Agustin LSA et al. 2001" |
#                                         Author_year =="Amadou A et al. 2015" |
#                                         #Author_year =="Cho E et al. 2003" |
#                                         Author_year =="Farvid MS et al. 2015" |
#                                         #Author_year =="Sieri S et al. 2013" |
#                                         Author_year =="Navarro-Silvera SA et al. 2005" |
#                                         Author_year =="Romanos‐Nanclares A et al. 2021" |
#                                         Author_year =="Larsson SC et al. 2009" |
#                                         Author_year =="Romieu I et al. 2012" |
#                                         Author_year =="Shikany JM et al. 2011"|
#                                         Author_year =="Makarem N et al. 2017" |
#                                         Author_year =="Lajous M et al. 2008" |
#                                         Author_year =="Lajous M et al. 2005" |
#                                         Author_year =="Wen W et al. 2009"  |
#                                         Author_year =="Yun SH et al. 2010"  |
#                                         Author_year =="Holmes MD et al. 2004" |
#                                         Author_year =="Higginbotham S et al. 2004" |
#                                         Author_year =="Jonas CR et al. 2003" |
#                                        #Author_year =="Nielsen TG et al. 2005" |
#                                         Author_year =="George SM et al. 2008"  # |
#                                        #Author_year =="Giles GC et al. 2006"
# )

#2.2 exposure
study_effects_gl <- dplyr::filter(study_effects_single, Exposure_recoded =="Glycemic load")
dim(study_effects_gl)

#2.3 All women
study_effects_gl_all <- dplyr::filter(study_effects_gl, Menopause_status =="All" & hormone_receptors == "NA")
str(study_effects_gl_all)

#2.4 fit
# fixed effects
meta_fe_gl<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                        se = SE, data = study_effects_gl_all, method = "fixed")
summary(meta_fe_gl)
predict(meta_fe_gl, delta = 1, expo = TRUE)
# Not significant log-linear dose-response

#2.5 fit
# Random effects
meta_re_gl<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                        se = SE, data = study_effects_gl_all, method = "reml")
summary(meta_re_gl)
predict(meta_re_gl, delta = 1, expo = TRUE)
# Not significant log-linear dose-response

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gl, delta = 10, exp = TRUE) 

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gl_all, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)

# 
dosex_bin <- data.frame(assigned_dose=seq(0, 200, 20))
with(predict(meta_re_gl, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic load, g/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model *NOT SIGNIFICANT
# square-transform the dose
quad_lin <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gl_all)
summary(quad_lin)
#plot
with(predict(quad_lin, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 3), ylab= 'Breast cancer relative risk', xlab= 'Glycemic load, g/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gl_all$assigned_dose, c(.05, .65, .95)) #c(.2, .4, .8) 
spl <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                  cases = cases, n = N, data = study_effects_gl_all)
summary(spl)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin,dataTab,exp=TRUE)
predSpl<-predict(spl,dataTab,exp=TRUE)
predSpl
round(cbind(lin=predLin,spl=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl), Sigma = vcov(spl), Terms = 1:2) #risk of cancer may not vary in a log-linear fashion with GI
# If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3 plot RCS
xref_bin <- 0
with(predict(spl, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 10), ylab="BC relative risk", xlab= "GL, grams/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 

################ PREMENOPAUSE ####################################
#2.3 All women
study_effects_gl_all_a <- dplyr::filter(study_effects_gl, Menopause_status =="Premenopause")
str(study_effects_gl_all_a) #24 

#2.4 fit
# fixed effects
meta_fe_g_a<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                         se = SE, data = study_effects_gl_all_a, method = "fixed")
summary(meta_fe_g_a)
predict(meta_fe_g_a, delta = 1, expo = TRUE)

#2.5 fit
# Random effects
meta_re_gl_a<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                          se = SE, data = study_effects_gl_all_a, method = "reml")
summary(meta_re_gl_a)
predict(meta_re_gl_a, delta = 1, expo = TRUE)
# significant log-linear dose-response association between gl and cancer risk (p=0.01)

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gl_a, delta = 10, exp = TRUE) # associated with a significant 3%

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gl_all_a, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)

# 
dosex_bin <- data.frame(assigned_dose=seq(0, 100, 1))
with(predict(meta_re_gl_a, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model *NOT SIGNIFICANT
# square-transform the dose
quad_lin_a <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gl_all_a)
summary(quad_lin_a)
#plot
with(predict(quad_lin_a, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 6), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin_a, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gl_all_a$assigned_dose,  probs = c(.35, .65, .95)) #c(.2, .4, .8) 
spl_a <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                    cases = cases, n = N, data = study_effects_gl_all_a)
summary(spl_a)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin_a,dataTab,exp=TRUE)
predSpl<-predict(spl_a,dataTab,exp=TRUE)
predSpl
round(cbind(quad_lin_a=predLin,spl_a=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl_a), Sigma = vcov(spl_a), Terms= 1:nrow(vcov(spl_a))) #risk of cancer may not vary in a log-linear fashion with gl
# The first of the two dose levels is excluded because it is the total raw data itself.If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3. plot RCS
xref_bin <- 0
with(predict(spl_a, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 5), ylab="BC relative risk", xlab= "gl, units/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 
################ POSTMENOPAUSE ####################################
#2.3 All women
study_effects_gl_all_b <- dplyr::filter(study_effects_gl, Menopause_status =="Posmenopause")
str(study_effects_gl_all_b)

# Nielsen TG et al. 2005 is in tertiles and IRR was givem
study_effects_gl_all_b <- dplyr::filter(study_effects_gl_all_b, Author_year !="Nielsen TG et al. 2005")

#2.4 fit
# fixed effects
meta_fe_g_b<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                         se = SE, data = study_effects_gl_all_b, method = "fixed")
summary(meta_fe_g_b)
predict(meta_fe_g_b, delta = 1, expo = TRUE)

#2.5 fit
# Random effects
meta_re_gl_b<- dosresmeta(yi ~ assigned_dose, id = id, type = type, cases = cases, n = N,
                          se = SE, data = study_effects_gl_all_b, method = "reml")
summary(meta_re_gl_b)
predict(meta_re_gl_b, delta = 1, expo = TRUE)
# significant log-linear dose-response association between gl and cancer risk (p=0.01)

# The change in cancer risk associated with every 10 units/day (standard drink) can be obtained with the predict function
predict(meta_re_gl_b, delta = 10, exp = TRUE) # associated with a significant 3%

#2.6 plot linear graph
library(ggplot2)
ggplot(study_effects_gl_all_b, aes(assigned_dose, Estimate, size= SE)) + 
  geom_point (shape=1, colour= 'black') + 
  scale_size_area(max_size=20)

# 
dosex_bin <- data.frame(assigned_dose=seq(0, 100, 1))
with(predict(meta_re_gl_b, dosex_bin, order=TRUE, exp=TRUE), {plot(assigned_dose, pred, type='l', col= 'blue', ylim=c(0, 2.5), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})

#2.7 Quadratic model
# square-transform the dose
quad_lin_b <- dosresmeta(formula = yi ~ assigned_dose + I(assigned_dose^2), id = id, type = type, se = SE,cases = cases, n = N, data = study_effects_gl_all_b)
summary(quad_lin_b)
#plot
with(predict(quad_lin_b, dosex_bin, exp=TRUE), {plot(assigned_dose, pred, type='l', ylim=c(0, 6), ylab= 'Breast cancer relative risk', xlab= 'Glycemic index, units/day')  
  lines(assigned_dose, ci.lb, lty= 2)  
  lines(assigned_dose, ci.ub, lty= 2)})
points(dosex_bin$assigned_dose, predict(quad_lin_b, dosex_bin, exp=TRUE)$pred, type='l', lty=3, col='blue')

#2.8 restricted cubic splines
library(rms)
k <- quantile(study_effects_gl_all_b$assigned_dose,  probs = c(.35, .65, .95)) #c(.2, .4, .8) 
spl_b <- dosresmeta(yi ~ rcs(assigned_dose, k), id = id, se = SE, type = type, 
                    cases = cases, n = N, data = study_effects_gl_all_b)
summary(spl_b)

dataTab<-data.frame(assigned_dose=seq(0,110,10))
predLin<-predict(quad_lin_b,dataTab,exp=TRUE)
predSpl<-predict(spl_b,dataTab,exp=TRUE)
predSpl
round(cbind(quad_lin_b=predLin,spl_b=predSpl[3:5]),2)

#2.9 Testing the linearity of dose level risk (regression coefficient)
waldtest(b = coef(spl_b), Sigma = vcov(spl_b), Terms= 1:nrow(vcov(spl_b))) #risk of cancer may not vary in a log-linear fashion with gl
# The first of the two dose levels is excluded because it is the total raw data itself.If the p-value is smaller than 0.05, the null hypothesis is rejected, and the joint slope is not zero. Thus, because there is a slope or the slopes of the two dose levels have a difference, it is determined that the model has nonlinearity.

#3. plot RCS
xref_bin <- 0
with(predict(spl_b, dosex_bin, xref_bin, exp=TRUE),{plot (get("rcs(assigned_dose, k)assigned_dose"), pred, type="l", ylim=c(0.4, 5), ylab="BC relative risk", xlab= "gl, units/day", log="y", bty="l", las=1)
  matlines(get("rcs(assigned_dose, k)assigned_dose"), cbind(ci.ub, ci.lb), col=1, lty ="dashed")}) 

############################## END ##############################################
