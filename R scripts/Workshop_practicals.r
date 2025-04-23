###########################################################
############## DATA BIOKO ISLAND       ####################
############## COOK ET AL (2011)       ####################
############## PLOS ONE 6(9): e25137   ####################
###########################################################

############## LOADING THE DATA ###########################

file<-'~/Desktop/Talks_Posters/Year_2025/Workshop_malaria/data_bioko.csv'

data.bioko<-read.csv(file,header=T,stringsAsFactors=F)

colnames(data.bioko)

###########################################################
##############           PART I             ###############
###########################################################
############## ANTIBODIES AGAINST PFAMA1     ############## 
############## ESTIMATING SEROPREVALENCE     ##############
############## USING GAUSSIAN MIXTURE MODELS ##############
###########################################################

########################################################### 
############## SUMMARY STATISTICS #########################
##############  PLOTTING          ######################### 
########################################################### 

antibody<-data.bioko$Ab_pfama1_titers

summary(antibody)

hist(antibody,xlab='antibody titres',main='Antibodies against PfAMA1')

###################################################
############## MODEL ESTIMATION ##################
###################################################

library(mixtools)

res<-normalmixEM(antibody,k=2,mean.constr = c(NA,'<'))

summary(res)

plot(res,whichplots=2,density=T,las=1,main2='Model fitting',col2=c('green','red'))

legend(2500,0.0030,c('Seronegative population','Seropositive population'),lwd=2,col=c('green','red'),xjust=1)

### First component refers to seronegative population ###

threshold<-res$mu[1]+3*res$sigma[1]

threshold

abline(v=threshold,lwd=5,col='orange')

seroprevalence<-mean(antibody>threshold)

seroprevalence

###########################################################
##############           PART II            ###############
###########################################################
############## FITTING REVERSIBLE CATALYTIC  ############## 
############## MODEL                         ##############
###########################################################

file<-'~/Desktop/SERO-AID/Source_codes/sero_aid_source_code.R'

source(file)

data.bioko$Seropositive<-as.numeric(antibody>threshold)

data.seroaid.bioko<-create.seroaid.object(data.bioko$Age,data.bioko$Seropositive,data.bioko$Site)

age.prof.overall<-age.profile(data.seroaid.bioko,analysis='overall',lag=0.1)

plot.seroaid(sero.obs=age.prof.overall)

rcm.overall<-rcm.constant.analysis(data.seroaid.bioko,analysis='overall')

rcm.overall$estimates

plot.seroaid(sero.obs=age.prof.overall,sero.model=rcm.overall)

age.profile.nw<-age.profile(data.seroaid.bioko,analysis='North West',lag=0.10)

age.profile.nw

rcm.nw<-rcm.constant.analysis(data.seroaid.bioko,analysis='North West')

rcm.nw$estimates

plot.seroaid(sero.obs=age.profile.nw,sero.model=rcm.nw)

######################################################################
##############           PART III                       ##############
######################################################################
############## SAMPLE SIZE CALCULATION                  ############## 
############## FOR ESTIMATING SCR FROM NORTH WEST BIOKO ############## 
######################################################################

library(RCMsize)

?sample_s

data.bioko.nw<-data.bioko[data.bioko$Site=='North West',]

summary(data.bioko.nw$Age_in_years)

age.distribution.frequency.nw<-table(factor(data.bioko.nw$Age_in_years,levels=1:max(data.bioko.nw$Age_in_years)))

age.distribution.proportion.nw<-age.distribution.frequency.nw/sum(age.distribution.frequency.nw)

plot(age.distribution.proportion.nw,xlab='Age in Years',ylab='Proportion',las=1)

expected.scr.nw<-rcm.nw$estimates[1,'lambda.est']

expected.srr.nw<-rcm.nw$estimates[1,'rho.est']

#### Sample size calculation for estimating expected SCR with relative length of 50%

sample_s(SCR=expected.scr.nw,RL=0.50,SRR=expected.srr.nw,ages=age.distribution.proportion.nw,A_max=90,limits=c(0,1),method='asymptotic')

#### Sample size calculation for estimating expected SCR with relative length of 25%

sample_s(SCR=expected.scr.nw,RL=0.25,SRR=expected.srr.nw,ages=age.distribution.proportion.nw,A_max=90,limits=c(0,1),method='asymptotic')

#### Sample size calculation for estimating expected SCR with relative length of 25%

sample_s(SCR=expected.scr.nw,RL=0.10,SRR=expected.srr.nw,ages=age.distribution.proportion.nw,A_max=90,limits=c(0,1),method='asymptotic')


