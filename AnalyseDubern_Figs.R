rm(list=setdiff(ls(), "filepath"))
####################################################################################################
# Underlying data (Dubern 1994 Tropical Science) is reproduced in the                              #
# manuscript Donnelly and Gilligan 2023 in which results (see Fig.s 2, 3, Fig. S3.1, Table S3.2)   # 
# obtained with the present files (AnalyseDubern_Figs.R, model_1.stan, model_2.stan) are presented #
####################################################################################################
require(lme4)
library('nlme')  
library('MASS')  
library(pracma) 
library(mc2d) 
library(rstan) 
library(grid) 
library(gridExtra)
library(lmtest)
library(bridgesampling)
library(data.table)
library(rstudioapi)
library(ggplot2)
library(ggpmisc)
library(ggtext)
library(ggpubr)

options(mc.cores = parallel::detectCores()) 
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
#setwd("C:/Users/donnelr3/OneDrive - University Of Cambridge/Archive_ZC_RetentionPeriod/rFiles/RPanalysis")


#(m2_relVE_Distr_summary,'Daily transfer no.',("Vector efficiency  \n (P(acquistion)*P(inoculation))"),c(0,13*1.05),c(0,1*1.05),seq(1,13,by=1),seq(0,1,by=0.2),trString,seq(0,1,by=0.2))+
#m2_relVE_Distr_summary$mn_estimate
#seq(1,13,by=1)
#
#
tailPlot <- function(dataIn, inStrX, inStrY, xRange, yRange, xExtent, yExtent,txtIn, txtLoc){ 
  ggplot(dataIn, aes(x=xExtent,y=mn_estimate))+
    geom_errorbar(aes(ymin=sm2_5, ymax=sm97_5), width=0)+
    geom_point(aes(x=xExtent, y=sm50), size=1,color="black",show.legend = TRUE)+
    geom_hline(aes(yintercept=0),size=0.1)+
    xlab(inStrX)+
    ylab(inStrY)+
    #annotate("text", x=txtLoc[1], y=txtLoc[2], label= txtIn, size= 2.5)+
    #coord_cartesian(xlim = xRange, ylim = yRange)+
    geom_vline(xintercept = 0)+
    scale_y_continuous(expand = c(0, 0), limits = yRange, breaks = yExtent)+
    scale_x_continuous(expand = c(0, 0), limits = xRange, breaks = xExtent)
}

################################################################################################################

numWarm=400
numIter=2000
adaptVal=1-(10^-1)
numChains=4
treeDepth=15

################################################################################################################
################### 1 # Laboratory experiment fitting  Dubern 1994 #############################################
################################################################################################################

# DUBERN RETENTION DATA (see Table 1 of corresponding Donnelly and GIlligan 2023 ms)
vD1=c(0,1,2,3,4,5,6,7,8,9,10,11,12);
vD2=c(28,24,20,20,19,12,8,3,3,1,0,0,0);
vD3=c(30,30,27,23,22,19,15,12,8,7,5,3,1);
transfers=vD1
dubRData=rbind(vD1,vD2,vD3);
datasizeR<-dim(dubRData)

dat3= list(D_Ntransfers=length(vD1),
           D_Wf0=10,
           D_RepsR=30,
           D_extant1=vD3,
           D_infn=vD2[seq(1,length(vD2))])


################################################################################################################
##################  FIT FOR H. TEST PERSPECTIVE  ###############################################################
################################################################################################################
initVal=0.0001;
initf1 <- function() list(compI=array(initVal,1),
                          rhoI=array(initVal,1),
                          rhoE=array(initVal,1),
                          rhoE0=array(initVal,1))
fitB1 = stan(file = "model_1.stan",
            data = dat3,
            init=initf1,
            iter = numIter,
            warmup = numWarm,
            chains=numChains,
            control = list(max_treedepth = treeDepth,adapt_delta=adaptVal) )
stanOutB1=summary(fitB1)
stanDetailB1=stanOutB1$summary
compI_fit_B1=summary(fitB1,'compI[1]')$summary[c(1,4,6,8)]
delPeriod_fit_B1=summary(fitB1,'delPeriod[1]')$summary[c(1,4,6,8)]
compI_del_fit_B1=summary(fitB1,'compI_del[1]')$summary[c(1,4,6,8)]
loseInf_fit_B1=summary(fitB1,'rhoI[1]')$summary[c(1,4,6,8)]
pDie_fit_B1=summary(fitB1,'rhoE[1]')$summary[c(1,4,6,8)]
pSurvInitial_fit_B1=summary(fitB1,'rhoE0[1]')$summary[c(1,4,6,8)]
pLS_fit_B1=summary(fitB1,'lifespan[1]')$summary[c(1,4,6,8)]

summary(fitB1,'probCohortExtinctByRE[30]')$summary[c(1,4,6,8)]        

# Table 1 components
outs_model1=rbind(pSurvInitial_fit_B1,
  pDie_fit_B1,
  loseInf_fit_B1,
  compI_fit_B1,
  delPeriod_fit_B1,
  pLS_fit_B1,
  compI_del_fit_B1)
colnames(outs_model1)=c('mn_estimate','sm2_5','sm50','sm97_5')



trString=c('1','2','3','4','5','6','7','8','9','10','11','12','13')

##### EXTANT PREDICTION GRAPHING ##############################################    FIG 2A
pExtant_Distr_summary=numeric(0)
for (ww in 1:length(vD1)) 
  pExtant_Distr_summary=rbind(pExtant_Distr_summary,summary(fitB1,paste('pExtant_Distr[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(pExtant_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
pExtant_Distr_summary=data.frame(pExtant_Distr_summary,empirData=vD3)     # affix corresponding empirical value
pl1=tailPlot(pExtant_Distr_summary,element_blank(),'No. living cohorts     ',c(0,1.05*length(vD1)),c(0,1.05*30),seq(1,13,by=1),seq(0,30,by=5),trString,seq(0,30,by=5))+
  geom_point(aes(x=seq(1,13,by=1), y=empirData), size=1,color="blue",show.legend = TRUE)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()
##### INFECTION PREDICTION GRAPHING ###########################################    FIG 2B 
pInf_Distr_summary=numeric(0)
for (ww in 1:length(vD1)) 
  pInf_Distr_summary=rbind(pInf_Distr_summary,summary(fitB1,paste('probInfSampRE[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(pInf_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
pInf_Distr_summary=data.frame(pInf_Distr_summary,empirData=vD2)           # affix corresponding empirical value
pl2=tailPlot(pInf_Distr_summary,element_blank(),'No. infected test plants    ',c(0,1.05*length(vD1)),c(0,1.05*30),seq(1,13,by=1),seq(0,30,by=5),trString,seq(0,30,by=5))+
  geom_point(aes(x=seq(1,13,by=1), y=empirData), size=1,color="blue",show.legend = TRUE)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()
######### PROBABILITY INOCULUM WAS PRESENT GRAPHING ########################### 
anyInoculumArray=numeric(0)
for (ww in 1:(3*length(vD1)))
  anyInoculumArray=rbind(anyInoculumArray,summary(fitB1,paste('probCohortExtinctByRE[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(anyInoculumArray)=c('mn_estimate','sm2_5','sm50','sm97_5')
anyInoculumArray=data.frame(anyInoculumArray)

temp <- "alpha == "
plas=tailPlot(anyInoculumArray[length(vD1)+1:length(vD1),],'Daily transfer no.',"Probability inoculum present      ",c(0,1.05*length(vD1)),c(0,1.05*1),seq(1,13,by=1),seq(0,1,by=0.2),trString,seq(0,1,by=0.2))+
  geom_point(aes(x=seq(1,13,by=1), y=anyInoculumArray[1:length(vD1),3]), size=1,color="blue",show.legend = TRUE)+
  geom_point(aes(x=seq(1,13,by=1), y=anyInoculumArray[2*length(vD1)+1:length(vD1),3]), size=1,color="green4",show.legend = TRUE)+
  annotate("text", x=6.05, y=0.0325,fontface ="bold", label= paste(temp, 0.25), colour = "blue", size= 2.75, parse = TRUE)+
  annotate("text", x=8.2, y=0.0325,fontface ="bold", label= paste(temp, 0.5), colour = "black", size= 2.75, parse = TRUE)+
  annotate("text", x=11.35, y=0.0325,fontface ="bold", label= paste(temp, 0.75), colour = "darkgreen", size= 2.75, parse = TRUE)+
  geom_hline(yintercept=0.95, colour="red", size=0.05, linetype = "dashed")+
  geom_segment(x = 7.2, y = 0, xend = 7.2, yend = anyInoculumArray[max(which(1*(anyInoculumArray[1:10,3]>0.95)==1)),3], colour = "blue", size=0.35, linetype = "dashed")+
  geom_segment(x = 10.2, y = 0, xend = 10.2, yend = anyInoculumArray[1*length(vD1)+max(which(1*(anyInoculumArray[(1*length(vD1)+1):(2*length(vD1)),3]>0.95)==1)),3], colour = "green4", size=0.35, linetype = "dashed")+
  geom_segment(x = 9.2, y = 0, xend = 9.2, yend = anyInoculumArray[2*length(vD1)+max(which(1*(anyInoculumArray[(2*length(vD1)+1):(3*length(vD1)),3]>0.95)==1)),3], colour = "black", size=0.35, linetype = "dashed")+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()
##########################################


################################################################################################################
################################################################################################################
################################################################################################################










################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
#                   FIT 'HIGHER RESOLUTION' VERSION - INFECTIOUS PERIOD PERSPECTIVE (cf. H. TEST PERSPECTIVE)  #
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

dat3B= list(D_Ntransfers=length(vD1),
           D_Wf0=10,
           D_RepsR=30,
           D_extant1=vD3,
           D_infn=vD2[seq(1,length(vD2))])
initVal=0.0001;
initf1 <- function() list(compI=array(initVal,dim=dat3$D_Ntransfers),
                          rhoE=array(initVal,1),
                          rhoE0=array(initVal,1))
fitC1 = stan(file = "model_2.stan",
             data = dat3B,
             init=initf1,
             iter = numIter,
             warmup = numWarm,
             chains=4,
             control = list(max_treedepth = treeDepth,adapt_delta=adaptVal))
stanOutC1=summary(fitC1)
stanDetailC1=stanOutC1$summary
compI_fit_C1=summary(fitC1,'compI[1]')$summary[c(1,4,6,8)]
compI_fit_C2=summary(fitC1,'compI[2]')$summary[c(1,4,6,8)]
compI_fit_C3=summary(fitC1,'compI[3]')$summary[c(1,4,6,8)]
compI_fit_C4=summary(fitC1,'compI[4]')$summary[c(1,4,6,8)]
compI_fit_C5=summary(fitC1,'compI[5]')$summary[c(1,4,6,8)]
compI_fit_C6=summary(fitC1,'compI[6]')$summary[c(1,4,6,8)]
compI_fit_C7=summary(fitC1,'compI[7]')$summary[c(1,4,6,8)]
compI_fit_C8=summary(fitC1,'compI[8]')$summary[c(1,4,6,8)]
compI_fit_C9=summary(fitC1,'compI[9]')$summary[c(1,4,6,8)]
compI_fit_C10=summary(fitC1,'compI[10]')$summary[c(1,4,6,8)]
compI_fit_C11=summary(fitC1,'compI[11]')$summary[c(1,4,6,8)]
compI_fit_C12=summary(fitC1,'compI[12]')$summary[c(1,4,6,8)]
compI_fit_C13=summary(fitC1,'compI[13]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_1=summary(fitC1,'compI_rel[1]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_2=summary(fitC1,'compI_rel[2]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_3=summary(fitC1,'compI_rel[3]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_4=summary(fitC1,'compI_rel[4]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_5=summary(fitC1,'compI_rel[5]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_6=summary(fitC1,'compI_rel[6]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_7=summary(fitC1,'compI_rel[7]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_8=summary(fitC1,'compI_rel[8]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_9=summary(fitC1,'compI_rel[9]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_10=summary(fitC1,'compI_rel[10]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_11=summary(fitC1,'compI_rel[11]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_12=summary(fitC1,'compI_rel[12]')$summary[c(1,4,6,8)]
compI_rel_fit_C1_13=summary(fitC1,'compI_rel[13]')$summary[c(1,4,6,8)]
compImax_fit_C1=summary(fitC1,'compIMaxVal[1]')$summary[c(1,4,6,8)]
summary(fitC1,'loc_compIMaxVal[1]')$summary[c(1,4,6,8)] # for reference
pDie_fit_C1=summary(fitC1,'rhoE[1]')$summary[c(1,4,6,8)]
pSurvInitial_fit_C1=summary(fitC1,'rhoE0[1]')$summary[c(1,4,6,8)]

outs_model2=rbind(compI_fit_C1,compI_fit_C2,compI_fit_C3,compI_fit_C4,compI_fit_C5,
                  compI_fit_C6,compI_fit_C7,compI_fit_C8,compI_fit_C9,compI_fit_C10,
                  compI_fit_C11,compI_fit_C12,compI_fit_C13,compI_rel_fit_C1_1,
                  compI_rel_fit_C1_2,compI_rel_fit_C1_3,compI_rel_fit_C1_4,
                  compI_rel_fit_C1_5,compI_rel_fit_C1_6,compI_rel_fit_C1_7,
                  compI_rel_fit_C1_8,compI_rel_fit_C1_9,compI_rel_fit_C1_10,
                  compI_rel_fit_C1_10,compI_rel_fit_C1_11,compI_rel_fit_C1_12,
                  compI_rel_fit_C1_13,compImax_fit_C1,pDie_fit_C1,pSurvInitial_fit_C1)
colnames(outs_model2)=c('mn_estimate','sm2_5','sm50','sm97_5')


##### EXTANT PREDICTION GRAPHING ############################################## FIG 2C 
m2_pExtant_Distr_summary=numeric(0)
for (ww in 1:length(vD1))
  m2_pExtant_Distr_summary=rbind(m2_pExtant_Distr_summary,summary(fitC1,paste('probSampExtant[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(m2_pExtant_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
m2_pExtant_Distr_summary=data.frame(m2_pExtant_Distr_summary,empirData=vD3)
ply1=tailPlot(m2_pExtant_Distr_summary,'Daily transfer no.','No. living cohorts     ',c(0,1.05*length(vD1)),c(0,1.025*30),seq(1,13,by=1),seq(0,30,by=5),trString,seq(0,30,by=5))+
  geom_point(aes(x=1:length(vD1), y=empirData), size=1,color="blue",show.legend = TRUE)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()

##### INFECTION PREDICTION GRAPHING ############################################## FIG 2D 
m2_pInf_Distr_summary=numeric(0)
for (ww in 1:length(vD1)) 
  m2_pInf_Distr_summary=rbind(m2_pInf_Distr_summary,summary(fitC1,paste('probInfSampRE[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(m2_pInf_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
m2_pInf_Distr_summary=data.frame(m2_pInf_Distr_summary,empirData=vD2)     # affix corresponding empirical value
ply2=tailPlot(m2_pInf_Distr_summary,'Daily transfer no.','No. infected test plants    ',c(0,1.05*length(vD1)),c(0,1.025*30),seq(1,13,by=1),seq(0,30,by=5),trString,seq(0,30,by=5))+
  geom_point(aes(x=1:length(vD1), y=empirData), size=1,color="blue",show.legend = TRUE)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()


##### PERCENTAGE PEAK GRAPHING #################################################  FIG 3A
m2_propPEAK_Distr_summary=numeric(0)
for (ww in 1:length(vD1)) 
  m2_propPEAK_Distr_summary=rbind(m2_propPEAK_Distr_summary,summary(fitC1,paste('compI_rel[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(m2_propPEAK_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
m2_propPEAK_Distr_summary=data.frame(m2_propPEAK_Distr_summary,empirData=vD2)     # affix corresponding empirical value
ppPreS_p2=tailPlot(m2_propPEAK_Distr_summary*100,'Daily transfer no.','% Peak insect infectiousness       ',c(0,1.05*length(vD1)),c(0,1.025*100),seq(1,13,by=1),seq(0,100,by=20),trString,seq(0,100,by=20))+
  geom_segment(x = 7.2, y = 0, xend = 7.2, yend = 100, colour = "blue", size=0.35, linetype = "dashed")+
  geom_segment(x = 10.2, y = 0, xend = 10.2, yend = 100, colour = "green4", size=0.35, linetype = "dashed")+
  geom_segment(x = 9.2, y = 0, xend = 9.2, yend = 100, colour = "black", size=0.35, linetype = "dashed")+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()

##### IDENTIFIABILITY GRAPHING v2 ##############################################  FIG S3.1
m2_relVE_Distr_summary=numeric(0)
for (ww in 1:dat3$D_Ntransfers)
  m2_relVE_Distr_summary=rbind(m2_relVE_Distr_summary,summary(fitC1,paste('compI[',ww,']',sep=""))$summary[c(1,4,6,8)]);
colnames(m2_relVE_Distr_summary)=c('mn_estimate','sm2_5','sm50','sm97_5')
m2_relVE_Distr_summary=data.frame(m2_relVE_Distr_summary)
ppPreS_ve=tailPlot(m2_relVE_Distr_summary,'Daily transfer no.',("Vector efficiency  \n (P(acquistion)*P(inoculation))"),c(0,13*1.05),c(0,1*1.05),seq(1,13,by=1),seq(0,1,by=0.2),trString,seq(0,1,by=0.2))+
  geom_segment(x = 7.2, y = 0, xend = 7.2, yend = 100, colour = "blue", size=0.35, linetype = "dashed")+
  geom_segment(x = 10.2, y = 0, xend = 10.2, yend = 100, colour = "green4", size=0.35, linetype = "dashed")+
  geom_segment(x = 9.2, y = 0, xend = 9.2, yend = 100, colour = "black", size=0.35, linetype = "dashed")+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10))+theme_classic()

#################################################################################



######### below is for reference only, ignore #
# model 1
##### EXTANT PREDICTION GRAPHING ############################################## FIG 2 A pExtant_Distr_summary pl1
##### INFECTION PREDICTION GRAPHING ########################################### FIG 2 B pInf_Distr_summary pl2
##### WAS THERE ACTUALLY ANY INOCULUM ?! ###################################### FIG 3 B anyInoculumArray plas
# model 2
##### EXTANT PREDICTION GRAPHING ############################################## FIG 2C m2_pExtant_Distr_summary ply1
##### INFECTION PREDICTION GRAPHING ########################################### FIG 2D m2_pInf_Distr_summary ply2
##### IDENTIFIABILITY GRAPHING v2 ############################################# FIG S3.1 m2_relVE_Distr_summary ppPreS_ve

figure=ggarrange(pl1,pl2,ply1,ply2,widths = c(0.5,0.5),labels = c("A", "B", "C", "D"),ncol = 2, nrow = 2)
figure=annotate_figure(figure,
                       left = text_grob(" \n \n \n \n Model 2    ", color = "black", face = "bold",
                                        hjust = 1, x = 1),fig.lab = " \n \n   Model 1", fig.lab.face = "bold")

ggsave("RP_Figure2.pdf", width = 0.85*19, height = 0.85*13.5, units = "cm")


ggarrange(plas,ppPreS_p2,heights = c(0.5,0.5),labels = c("A", "B"),ncol = 2, nrow = 1)

ggsave("RP_Figure3.pdf", width = 0.85*19, height = 0.5*13.5, units = "cm")

ggarrange(ppPreS_ve,ncol = 1, nrow = 1)

ggsave("RP_FigureS3_1.pdf", width = 0.85*19*0.6, height = 0.85*13.5*0.6, units = "cm")
################################################################################################################
################################################################################################################

