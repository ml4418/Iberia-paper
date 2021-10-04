rm(list = ls())
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(forcats)){ install.packages("forcats");library(forcats)}
wd<-"D:"
get_MTWA<-function(MTCO,GDD0){
  if(!is.na(MTCO)&!is.na(GDD0)){
    if(MTCO<0){
      GDD0_rad<-GDD0*2*pi/365
      if(!require(rootSolve)){ install.packages("rootSolve");library(rootSolve)}
      #- GDD0 /Tmin = 2 {[(u/(1 - u)] cos-1 (-u) + ???[(1 + u)/(1 ?- u)]}.
      tryCatch({
        eq<-function(u) (2*(acos(-u)*u/(1-u)+sqrt((1+u)/(1-u))))-GDD0_rad/(-MTCO)
        #curve(eq(x), -1, 1)
        #u should be between -1 and 1
        #u<-uniroot(eq,c(-1, 1))$root
        u<-uniroot(eq,c(-1+10^(-10), 1-10^(-10)))$root
        
        #Tmax = T0 + ??T =  - Tmin (1+u)/(1 - u)
        MTWA<- (-MTCO)*(1+u)/(1-u)
      }, error=function(e){MTWA<-NA}
      )
      
    }else{
      T0<-GDD0/365
      MTWA<-2*T0-MTCO
    }
  }else{
    MTWA<-NA
  }
  
  return(MTWA)
}
########################################################################################################
#######################################    Training  ###################################################
########################################################################################################
#install.packages("fxTWAPLS")
install.packages("remotes")
remotes::install_github("special-uor/fxTWAPLS", "dev")
source(paste(wd,"/Master Project/Script/Reconstruction paper/Iberian Peninsula reconstruction functions.R",sep=""))
source(paste(wd,"/Master Project/Script/fx explorations/WAPLS2 and TWAPLS2.R",sep=""))

modern_pollen<- read.csv(paste(wd,"/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv",sep=""), row.names=1)

taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

for(i in 1:nrow(modern_pollen)){
  #i=i+1
  modern_pollen[i,"Tmax"]<-get_MTWA(MTCO=modern_pollen[i,"Tmin"],GDD0=modern_pollen[i,"gdd"])
}
summary(modern_pollen$Tmax);hist(modern_pollen$Tmax)
# Training ---------------------------------------------

#compare training results
fit_tf1_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf2_Tmin <- TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf1_Tmin_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_tf2_Tmin_pspline <- TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

fit_tf1_Tmax <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf2_Tmax <- TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf1_Tmax_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_tf2_Tmax_pspline <- TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

fit_tf1_alpha <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.002)
fit_tf2_alpha <- TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.002)
fit_tf1_alpha_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)
fit_tf2_alpha_pspline <- TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)

plot.train.sig<-function(name,fit_tf1_Tmin,fit_tf2_Tmin,nsig_tf1_Tmin,nsig_tf2_Tmin,
                         fit_tf1_Tmax,fit_tf2_Tmax,nsig_tf1_Tmax,nsig_tf2_Tmax,
                         fit_tf1_alpha,fit_tf2_alpha,nsig_tf1_alpha,nsig_tf2_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #Get the training results
  train<-cbind.data.frame(fit_tf1_Tmin[["x"]],
                          fit_tf1_Tmin[["fit"]][,nsig_tf1_Tmin],
                          fit_tf2_Tmin[["fit"]][,nsig_tf2_Tmin],
                          fit_tf1_Tmax[["x"]],
                          fit_tf1_Tmax[["fit"]][,nsig_tf1_Tmax],
                          fit_tf2_Tmax[["fit"]][,nsig_tf2_Tmax],
                          fit_tf1_alpha[["x"]],
                          fit_tf1_alpha[["fit"]][,nsig_tf1_alpha],
                          fit_tf2_alpha[["fit"]][,nsig_tf2_alpha])
  
  colnames(train)<-c("Tmin","Tmin_tf1","Tmin_tf2","Tmax","Tmax_tf1","Tmax_tf2","alpha","alpha_tf1","alpha_tf2")
  max_Tmin<-max(train[,c("Tmin","Tmin_tf1","Tmin_tf2")])
  min_Tmin<-min(train[,c("Tmin","Tmin_tf1","Tmin_tf2")])
  max_Tmax<-max(train[,c("Tmax","Tmax_tf1","Tmax_tf2")])
  min_Tmax<-min(train[,c("Tmax","Tmax_tf1","Tmax_tf2")])
  max_alpha<-max(train[,c("alpha","alpha_tf1","alpha_tf2")])
  min_alpha<-min(train[,c("alpha","alpha_tf1","alpha_tf2")])
  
  p_Tmin_tf1<-ggplot(data=train,aes(Tmin,Tmin_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmin,max_Tmin)+xlim(min_Tmin,max_Tmin)+
    annotate("text", y= max_Tmin, x = min_Tmin,label="(a)")+
    labs(y="fxTWA-PLS1 MTCO (°C)",x="MTCO (°C)")
  p_Tmin_tf2<-ggplot(data=train,aes(Tmin,Tmin_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmin,max_Tmin)+xlim(min_Tmin,max_Tmin)+
    annotate("text", y= max_Tmin, x = min_Tmin,label="(b)")+
    labs(y="fxTWA-PLS2 MTCO (°C)",x="MTCO (°C)")
  
  
  p_Tmax_tf1<-ggplot(data=train,aes(Tmax,Tmax_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmax,max_Tmax)+xlim(min_Tmax,max_Tmax)+
    annotate("text", y= max_Tmax, x = min_Tmax,label="(c)")+
    labs(y="fxTWA-PLS1 MTWA (°C)",x="MTWA (°C)")
  p_Tmax_tf2<-ggplot(data=train,aes(Tmax,Tmax_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmax,max_Tmax)+xlim(min_Tmax,max_Tmax)+
    annotate("text", y= max_Tmax, x = min_Tmax,label="(d)")+
    labs(y="fxTWA-PLS2 MTWA (°C)",x="MTWA (°C)")
  
  p_alpha_tf1<-ggplot(data=train,aes(alpha,alpha_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_alpha,max_alpha)+xlim(min_alpha,max_alpha)+
    annotate("text", y= max_alpha, x = min_alpha,label="(e)")+
    labs(y= expression("fxTWA-PLS1   "*alpha), x = expression(alpha))+
    geom_abline(slope=0, intercept=0,linetype="dashed")+
    geom_abline(slope=0, intercept=1.26,linetype="dashed")
  p_alpha_tf2<-ggplot(data=train,aes(alpha,alpha_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_alpha,max_alpha)+xlim(min_alpha,max_alpha)+
    annotate("text", y= max_alpha, x = min_alpha,label="(f)")+
    labs(y= expression("fxTWA-PLS2   "*alpha), x = expression(alpha))+
    geom_abline(slope=0, intercept=0,linetype="dashed")+
    geom_abline(slope=0, intercept=1.26,linetype="dashed")
  
  p<-ggarrange(p_Tmin_tf1,p_Tmin_tf2,p_Tmax_tf1,p_Tmax_tf2,p_alpha_tf1,p_alpha_tf2,  ncol = 2)
  ggsave(file=paste(name,"training results.jpeg"),p,width=8,height=9)
  
}

plot.train.sig("Bin version",fit_tf1_Tmin,fit_tf2_Tmin,4,3,
                         fit_tf1_Tmax,fit_tf2_Tmax,2,2,
                         fit_tf1_alpha,fit_tf2_alpha,3,4)
plot.train.sig("P-spline version",
               fit_tf1_Tmin_pspline,fit_tf2_Tmin_pspline,2,4,
               fit_tf1_Tmax_pspline,fit_tf2_Tmax_pspline,2,4,
               fit_tf1_alpha_pspline,fit_tf2_alpha_pspline,3,3)

# In step 7 of training, climate variable is regressed to the components obtained,
# MTCO
fit_tf_Tmin <- TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

# GDD0
fit_tf_gdd <- TWAPLS.w2(taxa, modern_pollen$gdd, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=20)

# alpha
fit_tf_alpha <- TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)

#MTWA
fit_tf_Tmax <- TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)


########################################################################################################
####################################  Cross validation  Table 1 #######################################
########################################################################################################
# Pseudo removed leave out cross validation
dist <- read.csv(paste(wd,"/Master Project/Data/Output data/Training fitness/Geographically and climatically close sites removed/distance.csv",sep=""), row.names = 1)
CPUS<-8
pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, modern_pollen$Tmin, cpus = CPUS)
pseudo_gdd <- fxTWAPLS::get_pseudo(dist, modern_pollen$gdd, cpus = CPUS)
pseudo_alpha <- fxTWAPLS::get_pseudo(dist, modern_pollen$alpha, cpus = CPUS)
pseudo_Tmax <- fxTWAPLS::get_pseudo(dist, modern_pollen$Tmax, cpus = CPUS)

if(!require(rlist)){install.packages("rlist");library(rlist)} 
rlist::list.save(pseudo_Tmin, 'pseudo_Tmin.rdata')
rlist::list.save(pseudo_gdd, 'pseudo_gdd.rdata')
rlist::list.save(pseudo_alpha, 'pseudo_alpha.rdata')
rlist::list.save(pseudo_Tmax, 'pseudo_Tmax.rdata')

# Pseudo removed leave out cross validation
setwd(paste(wd,"/Master Project/Data/Output data/fx explorations/Cross validation",sep=""))
pseudo_Tmin <- rlist::list.load('pseudo_Tmin.rdata')
pseudo_gdd <- rlist::list.load('pseudo_gdd.rdata')
pseudo_alpha <- rlist::list.load('pseudo_alpha.rdata')
pseudo_Tmax <- rlist::list.load('pseudo_Tmax.rdata')

#tf2 pspline
if(!require(foreach)){install.packages("foreach");library(foreach)}
cv_tf2_Tmin_pspline <- cv.pr.w(taxa,
                       modern_pollen$Tmin,
                       nPLS = 5,
                       TWAPLS.w2,
                       TWAPLS.predict.w,
                       pseudo_Tmin,
                       usefx = TRUE,
                       fx_method = "pspline",
                       bin = 0.02,
                       cpus = 8,
                       test_mode = F)   
write.csv(cv_tf2_Tmin_pspline, "cv_tf2_Tmin_pspline.csv")

cv_tf2_gdd_pspline <- cv.pr.w(taxa,
                       modern_pollen$gdd,
                       nPLS = 5,
                       TWAPLS.w2,
                       TWAPLS.predict.w,
                       pseudo_gdd,
                       usefx = TRUE,
                       fx_method = "pspline",
                       bin = 20,
                       cpus = 8,
                       test_mode = F)   
write.csv(cv_tf2_gdd_pspline, "cv_tf2_gdd_pspline.csv")

cv_tf2_alpha_pspline <- cv.pr.w(taxa,
                       modern_pollen$alpha,
                       nPLS = 5,
                       TWAPLS.w2,
                       TWAPLS.predict.w,
                       pseudo_alpha,
                       usefx = TRUE,
                       fx_method = "pspline",
                       bin = 0.002,
                       cpus = 8,
                       test_mode = F)   
write.csv(cv_tf2_alpha_pspline, "cv_tf2_alpha_pspline.csv")

cv_tf2_Tmax_pspline <- cv.pr.w(taxa,
                       modern_pollen$Tmax,
                       nPLS = 5,
                       TWAPLS.w2,
                       TWAPLS.predict.w,
                       pseudo_Tmax,
                       usefx = TRUE,
                       fx_method = "pspline",
                       bin = 0.02,
                       cpus = 8,
                       test_mode = F)   
write.csv(cv_tf2_Tmax_pspline, "cv_tf2_Tmax_pspline.csv")

rand_tf2_Tmin_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmin_pspline, n.perm = 999) 
rand_tf2_alpha_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_alpha_pspline, n.perm = 999) 
rand_tf2_Tmax_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmax_pspline, n.perm = 999) 
rand_tf2_gdd_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_gdd_pspline, n.perm = 999) 

rand_tf2_pspline<-rbind.data.frame(rand_tf2_Tmin_pspline,rand_tf2_Tmax_pspline,rand_tf2_alpha_pspline,rand_tf2_gdd_pspline)
write.csv(rand_tf2_pspline, "rand_tf2_pspline.csv")

#tf2 bin
`%>%` <- magrittr::`%>%`
if(!require(foreach)){install.packages("foreach");library(foreach)}
cv_tf2_Tmin_bin <- cv.pr.w(taxa,
                               modern_pollen$Tmin,
                               nPLS = 5,
                               TWAPLS.w2,
                               TWAPLS.predict.w,
                               pseudo_Tmin,
                               usefx = TRUE,
                               fx_method = "bin",
                               bin = 0.02,
                               cpus = 8,
                               test_mode = F)   
write.csv(cv_tf2_Tmin_bin, "cv_tf2_Tmin_bin.csv")

cv_tf2_gdd_bin <- cv.pr.w(taxa,
                              modern_pollen$gdd,
                              nPLS = 5,
                              TWAPLS.w2,
                              TWAPLS.predict.w,
                              pseudo_gdd,
                              usefx = TRUE,
                              fx_method = "bin",
                              bin = 20,
                              cpus = 8,
                              test_mode = F)   %>% fxTWAPLS::pb()
write.csv(cv_tf2_gdd_bin, "cv_tf2_gdd_bin.csv")

cv_tf2_alpha_bin <- cv.pr.w(taxa,
                                modern_pollen$alpha,
                                nPLS = 5,
                                TWAPLS.w2,
                                TWAPLS.predict.w,
                                pseudo_alpha,
                                usefx = TRUE,
                                fx_method = "bin",
                                bin = 0.002,
                                cpus = 8,
                                test_mode = F)   
write.csv(cv_tf2_alpha_bin, "cv_tf2_alpha_bin.csv")

cv_tf2_Tmax_bin <- cv.pr.w(taxa,
                               modern_pollen$Tmax,
                               nPLS = 5,
                               TWAPLS.w2,
                               TWAPLS.predict.w,
                               pseudo_Tmax,
                               usefx = TRUE,
                               fx_method = "bin",
                               bin = 0.02,
                               cpus = 8,
                               test_mode = F)   
write.csv(cv_tf2_Tmax_bin, "cv_tf2_Tmax_bin.csv")

rand_tf2_Tmin_bin <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmin_bin, n.perm = 999) 
rand_tf2_gdd_bin <- fxTWAPLS::rand.t.test.w(cv_tf2_gdd_bin, n.perm = 999) 
rand_tf2_alpha_bin <- fxTWAPLS::rand.t.test.w(cv_tf2_alpha_bin, n.perm = 999) 
rand_tf2_Tmax_bin <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmax_bin, n.perm = 999) 
rand_tf2_bin<-rbind.data.frame(rand_tf2_Tmin_bin,rand_tf2_Tmax_bin,rand_tf2_alpha_bin)
write.csv(rand_tf2_bin, "rand_tf2_bin.csv")

#tf1 pspline
if(!require(foreach)){install.packages("foreach");library(foreach)}
cv_tf1_Tmin_pspline <- cv.pr.w(taxa,
                               modern_pollen$Tmin,
                               nPLS = 5,
                               TWAPLS.w,
                               TWAPLS.predict.w,
                               pseudo_Tmin,
                               usefx = TRUE,
                               fx_method = "pspline",
                               bin = 0.02,
                               cpus = 8,
                               test_mode = F)   
write.csv(cv_tf1_Tmin_pspline, "cv_tf1_Tmin_pspline.csv")

cv_tf1_gdd_pspline <- cv.pr.w(taxa,
                              modern_pollen$gdd,
                              nPLS = 5,
                              TWAPLS.w,
                              TWAPLS.predict.w,
                              pseudo_gdd,
                              usefx = TRUE,
                              fx_method = "pspline",
                              bin = 20,
                              cpus = 8,
                              test_mode = F)   
write.csv(cv_tf1_gdd_pspline, "cv_tf1_gdd_pspline.csv")

cv_tf1_alpha_pspline <- cv.pr.w(taxa,
                                modern_pollen$alpha,
                                nPLS = 5,
                                TWAPLS.w,
                                TWAPLS.predict.w,
                                pseudo_alpha,
                                usefx = TRUE,
                                fx_method = "pspline",
                                bin = 0.002,
                                cpus = 8,
                                test_mode = F)   
write.csv(cv_tf1_alpha_pspline, "cv_tf1_alpha_pspline.csv")

cv_tf1_Tmax_pspline <- cv.pr.w(taxa,
                               modern_pollen$Tmax,
                               nPLS = 5,
                               TWAPLS.w,
                               TWAPLS.predict.w,
                               pseudo_Tmax,
                               usefx = TRUE,
                               fx_method = "pspline",
                               bin = 0.02,
                               cpus = 8,
                               test_mode = F)   
write.csv(cv_tf1_Tmax_pspline, "cv_tf1_Tmax_pspline.csv")

rand_tf1_Tmin_pspline <- fxTWAPLS::rand.t.test.w(cv_tf1_Tmin_pspline, n.perm = 999) 
rand_tf1_alpha_pspline <- fxTWAPLS::rand.t.test.w(cv_tf1_alpha_pspline, n.perm = 999) 
rand_tf1_Tmax_pspline <- fxTWAPLS::rand.t.test.w(cv_tf1_Tmax_pspline, n.perm = 999) 
rand_tf1_gdd_pspline <- fxTWAPLS::rand.t.test.w(cv_tf1_gdd_pspline, n.perm = 999) 

rand_tf1_pspline<-rbind.data.frame(rand_tf1_Tmin_pspline,rand_tf1_Tmax_pspline,rand_tf1_alpha_pspline,rand_tf1_gdd_pspline)
write.csv(rand_tf1_pspline, "rand_tf1_pspline.csv")

#tf1 bin
if(!require(foreach)){install.packages("foreach");library(foreach)}
cv_tf1_Tmin_bin <- cv.pr.w(taxa,
                           modern_pollen$Tmin,
                           nPLS = 5,
                           TWAPLS.w,
                           TWAPLS.predict.w,
                           pseudo_Tmin,
                           usefx = TRUE,
                           fx_method = "bin",
                           bin = 0.02,
                           cpus = 8,
                           test_mode = F)   
write.csv(cv_tf1_Tmin_bin, "cv_tf1_Tmin_bin.csv")

cv_tf1_gdd_bin <- cv.pr.w(taxa,
                          modern_pollen$gdd,
                          nPLS = 5,
                          TWAPLS.w,
                          TWAPLS.predict.w,
                          pseudo_gdd,
                          usefx = TRUE,
                          fx_method = "bin",
                          bin = 20,
                          cpus = 8,
                          test_mode = F)   
write.csv(cv_tf1_gdd_bin, "cv_tf1_gdd_bin.csv")

cv_tf1_alpha_bin <- cv.pr.w(taxa,
                            modern_pollen$alpha,
                            nPLS = 5,
                            TWAPLS.w,
                            TWAPLS.predict.w,
                            pseudo_alpha,
                            usefx = TRUE,
                            fx_method = "bin",
                            bin = 0.002,
                            cpus = 8,
                            test_mode = F)   
write.csv(cv_tf1_alpha_bin, "cv_tf1_alpha_bin.csv")

cv_tf1_Tmax_bin <- cv.pr.w(taxa,
                           modern_pollen$Tmax,
                           nPLS = 5,
                           TWAPLS.w,
                           TWAPLS.predict.w,
                           pseudo_Tmax,
                           usefx = TRUE,
                           fx_method = "bin",
                           bin = 0.02,
                           cpus = 8,
                           test_mode = F)   
write.csv(cv_tf1_Tmax_bin, "cv_tf1_Tmax_bin.csv")

rand_tf1_Tmin_bin <- fxTWAPLS::rand.t.test.w(cv_tf1_Tmin_bin, n.perm = 999) 
rand_tf1_gdd_bin <- fxTWAPLS::rand.t.test.w(cv_tf1_gdd_bin, n.perm = 999) 
rand_tf1_alpha_bin <- fxTWAPLS::rand.t.test.w(cv_tf1_alpha_bin, n.perm = 999) 
rand_tf1_Tmax_bin <- fxTWAPLS::rand.t.test.w(cv_tf1_Tmax_bin, n.perm = 999) 
rand_tf1_bin<-rbind.data.frame(rand_tf1_Tmin_bin,rand_tf1_Tmax_bin,rand_tf1_alpha_bin)
write.csv(rand_tf1_bin, "rand_tf1_bin.csv")
########################################################################################################
####################################  Reconstruction  ##################################################
########################################################################################################
Iberian_all <- read.csv(paste(wd,"/Master Project/Data/Input data/Iberia_Pollen_add_INTCAL20_0313.csv",sep=""))
Iberian<-Iberian_all[which(Iberian_all$`INTCAL2020.median`<=12000),]
#Age
colnames(Iberian);str(Iberian)
summary(Iberian[,"INTCAL2020.median"])
#Exclude samples with large age uncertainty
hist((Iberian[,"INTCAL2020_95"]-Iberian[,"INTCAL2020_5"])/(2*1.645))
which_age_keep<-which(abs((Iberian[,"INTCAL2020_95"]-Iberian[,"INTCAL2020_5"])/(2*1.645))<=100)
Iberian<-Iberian[which_age_keep,]
write.csv(Iberian,paste(wd,"/Master Project/Data/Input data/Iberia.csv",sep=""))

#Reconstruction
taxaColMin <- which(colnames(Iberian) == "Abies")
taxaColMax <- which(colnames(Iberian) == "Zygophyllaceae")
core0<-Iberian[,taxaColMin:taxaColMax]
core0[is.na(core0)]<-0;str(core0)
core<-core0

#get common taxa
colnames(core)[!(colnames(core)%in% colnames(taxa))]
colnames(taxa)[!(colnames(taxa)%in% colnames(core))]

colnames(core)[(colnames(core)%in% colnames(taxa))]

taxa_to_delete<-colnames(core)[!(colnames(core)%in% colnames(taxa))]
core[,taxa_to_delete]<-NULL
taxa_to_add<-colnames(taxa)[!(colnames(taxa)%in% colnames(core))]
core[,taxa_to_add]<-0
core<-core[,order(colnames(core))]

summary(rowSums(core))
core<-core/rowSums(core)
summary(rowSums(core))
str(core)

setwd(paste(wd,"/Master Project/Data/Output data/Core/Iberian",sep=""))
#MTCO
fossil_tf_Tmin<-fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin,core)
#GDD0
fossil_tf_gdd<-fxTWAPLS::TWAPLS.predict.w(fit_tf_gdd,core)
#alpha
fossil_tf_alpha<-fxTWAPLS::TWAPLS.predict.w(fit_tf_alpha,core)
#MTWA
fossil_tf_Tmax<-fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmax,core)


#Use the last significant number of components
core_sig<-cbind.data.frame(Iberian[,c("Site.Name","Latitude","Longitude","Elevation","INTCAL2020.median")],
                  fossil_tf_Tmin[["fit"]][,4],fossil_tf_gdd[["fit"]][,2],fossil_tf_alpha[["fit"]][,3], fossil_tf_Tmax[["fit"]][,4])
colnames(core_sig)<-c("site","lat","lon","elv","age","Tmin","gdd","alpha","Tmax")
write.csv(core_sig,paste(wd,"/Master Project/Data/Output data/Core/Iberian/Iberian_core_sig.csv",sep=""))

########################################################################################################
####################################### Sitelist Table S1 ##############################################
########################################################################################################
#Get information for sites
#sitelist<-Iberian[!duplicated(Iberian$Site.Name),]
sitename<-unique(Iberian$Site.Name)
sitelist<-as.data.frame(sitename)
colnames(Iberian)[1:18]
for(i in 1:nrow(sitelist)){
  #i=i+1
  siteinfo<-Iberian[which(Iberian$Site.Name==sitelist[i,"sitename"]),]
  sitelist[i,c("lon","lat","elv")]<-unique(siteinfo[,c("Longitude","Latitude","Elevation")])
  ref<-unique(siteinfo[,"Reference"])
  if(all(ref=="")){sitelist[i,"ref"]<-""}else{sitelist[i,"ref"]<-ref[which(ref!="")]}
  sitelist[i,"length_yr"]<-max(siteinfo[,"INTCAL2020.median"])-min(siteinfo[,"INTCAL2020.median"])
  sitelist[i,"n_sample"]<-sum(!is.na(siteinfo[,"INTCAL2020.median"]))
  
}
sitelist<-sitelist[order(sitelist$lon, sitelist$lat), ]
#elv_label
sitelist$elv_label<-NA
for(j in 1:nrow(sitelist)){
  elv<-as.numeric(sitelist[j,"elv"])
  if(elv<=1000){
    sitelist[j,"elv_label"]<-"low"
  }else if(elv>1000){
    sitelist[j,"elv_label"]<-"high"
  }
}
write.csv(sitelist,paste(wd,"/Master Project/Data/Output data/Core/Iberian/sitelist.csv",sep=""))

########################################################################################################
#######################################    CCA  Table 2 ###############################################
########################################################################################################
if(!require(vegan)){ install.packages("vegan");library(vegan)}
ord<-cca(formula=taxa~modern_pollen$Tmin+modern_pollen$Tmax+modern_pollen$alpha)
ord
vif.cca(ord)
round(coef(ord),digits = 3)
anova(ord)
anova(ord, by="term", permutations=999)
anova(ord, by="axis", permutations=999)
plot(ord)

ord<-cca(formula=core~core_sig$Tmin+core_sig$Tmax+core_sig$alpha,na.action=na.omit)
ord
vif.cca(ord)
round(coef(ord),digits = 3)
anova(ord)
anova(ord, by="term", permutations=999)
anova(ord, by="axis", permutations=999)
plot(ord)

########################################################################################################
###########################    Calculate insolation at each site  #######################################
########################################################################################################
plotdata<-core_sig

#exlcude outliers of the climate varaibles
plotdata<-plotdata[which(plotdata$alpha>0&plotdata$alpha<1.26),]

for(i in 1:nrow(plotdata)){
  Lat=plotdata[i,"lat"]
  Age<-plotdata[i,"age"]
  plotdata[i,c("insol_summer","insol_winter",
               "Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")]<-get_insolation(Lat,age_start =Age,age_end=Age, interval=0 )[,c("insol_summer","insol_winter",
                                                                                                                          "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                                                                                          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")]
}
write.csv(plotdata,paste(wd,"/Master Project/Data/Output data/Core/Iberian/plotdata.csv",sep=""))

########################################################################################################
###########################   Relationship between MTWA and alpha  #####################################
########################################################################################################
#Figure 4
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))
summary(plotdata$Tmax);summary(modern_pollen$Tmax)

p1<-ggplot(modern_pollen,aes(alpha,Tmax))+geom_point(size=0.8)+theme_bw()+
  labs(x= expression(alpha), y = "MTWA (°C)")+ylim(0,35)+
  theme(legend.position = "none",
        axis.text=element_text(size=14),axis.title=element_text(size=14))+
  annotate("text", y= 35, x =0,label="(a)",size=6)
p2<-ggplot(plotdata,aes(alpha,Tmax))+geom_point(size=0.8)+theme_bw()+
  labs(x= expression(alpha), y = "MTWA (°C)")+ylim(0,35)+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.text=element_text(size=14),axis.title=element_text(size=14))+
  annotate("text", y= 35, x =0,label="(b)",size=6)
p<-ggarrange(p1,p2,ncol=2)
ggsave(file="Relationship between MTWA and alpha.jpeg",p,width=12,height=6)

########################################################################################################
#####################################   Map of fossil sites and alpha  #################################
########################################################################################################
#Map at 0.5 ka
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))
background1 <- read.csv(paste(wd,"/Master Project/Data/Input data/add on region_1.csv",sep=""), row.names=1)
background2 <- read.csv(paste(wd,"/Master Project/Data/Input data/add on region_2.csv",sep=""))
background1<-background1[,c("lat","lon","elv","mtco","gdd","MI")]
background2<-background2[,c("latitude","longitude","elevation_m","MTCM","GDD0","MI")]
colnames(background2)=colnames(background1)=c("lat","lon","elv","Tmin","gdd","MI")
background<-rbind.data.frame(background1,background2)
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}
if(!require(sf)){ install.packages("sf");library(sf)}
world <- map_data("world") 
summary(Iberian$Longitude);summary(Iberian$Latitude)
region<-world[which(world$long>(-10)&world$long<5&world$lat>36&world$lat<44),]
Iberian_background<-background[which(background$lon>(-10)&background$lon<5&background$lat>36&background$lat<44),]

#tranfer MI to alpha
MI<-Iberian_background$MI
w<-1.5
fai<-1/MI
F<-1+fai-(1+fai^w)^(1/w)
alpha<-1.26*MI*F
hist(alpha)
plot(alpha~MI)
Iberian_background$alpha<-alpha

#get MTWA
for(i in 1:nrow(Iberian_background)){
  #i=i+1
  Iberian_background[i,"Tmax"]<-get_MTWA(MTCO=Iberian_background[i,"Tmin"],GDD0=Iberian_background[i,"gdd"])
}

write.csv(Iberian_background,paste(wd,"/Master Project/Data/Input data/Iberian background.csv",sep=""))

Iberian_background<-read.csv(paste(wd,"/Master Project/Data/Input data/Iberian background.csv",sep=""),row.names = 1)


#plot 
p1<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=Tmin),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour ="MTCO (°C)")+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 °W","5 °W","0 °E","5 °E"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 °N","38 °N","40 °N","42 °N","44 °N"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  annotate("text",x=-10,y=44,label="(a)")
p2<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=Tmax),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour ="MTWA (°C)")+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 °W","5 °W","0 °E","5 °E"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 °N","38 °N","40 °N","42 °N","44 °N"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  annotate("text",x=-10,y=44,label="(b)")
p3<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=alpha),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = rev(c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red")),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour =expression(alpha))+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 °W","5 °W","0 °E","5 °E"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 °N","38 °N","40 °N","42 °N","44 °N"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  annotate("text",x=-10,y=44,label="(c)")
p<-ggarrange(p1,p2,p3,ncol=1)
ggsave(file="Map of modern alpha.jpeg",p,width=6,height=10)

########################################################################################################
######################################   Get the gradient  #############################################
########################################################################################################
get_gradient<-function(plotdata,age,range){
  plotdata$site<-as.factor(plotdata$site)
  mean_env<-data.frame(matrix(NA,ncol=5+4+2,nrow=nlevels(plotdata$site)))
  for(m in 1:nlevels(plotdata$site)){
    plotdata_each<-plotdata[which(plotdata$site==levels(plotdata$site)[m]),]
    
    mean_Tmin<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"Tmin"]),na.rm = TRUE)
    mean_gdd<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"gdd"]),na.rm = TRUE)
    mean_alpha<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"alpha"]),na.rm = TRUE)
    mean_Tmax<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"Tmax"]),na.rm = TRUE)
    
    mean_insol_summer<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"insol_summer"]),na.rm = TRUE)
    mean_insol_winter<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"insol_winter"]),na.rm = TRUE)
    mean_env[m,]<-c(age,levels(plotdata$site)[m],unique(plotdata_each$lon),unique(plotdata_each$lat),unique(plotdata_each$elv),mean_Tmin,mean_gdd,mean_alpha,mean_Tmax,mean_insol_summer,mean_insol_winter)
  }
  colnames(mean_env)<-c("age","site","lon","lat","elv","mean_Tmin","mean_gdd","mean_alpha","mean_Tmax","mean_insol_summer","mean_insol_winter")
  return(mean_env)
}

# Get fossil gradient 
gradient_env<-rbind.data.frame(
  get_gradient(plotdata,500,500),get_gradient(plotdata,1500,500),
  get_gradient(plotdata,2500,500),get_gradient(plotdata,3500,500),
  get_gradient(plotdata,4500,500),get_gradient(plotdata,5500,500),
  get_gradient(plotdata,6500,500),get_gradient(plotdata,7500,500),
  get_gradient(plotdata,8500,500),get_gradient(plotdata,9500,500),
  get_gradient(plotdata,10500,500),get_gradient(plotdata,11500,500))
#get modern gradient
modern_gradient<-get_gradient(plotdata,1000,500)
modern_gradient$age<-0
for(i in 1:nrow(modern_gradient)){
  tryCatch({
    sitelon<-as.numeric(modern_gradient[i,"lon"])
    sitelat<-as.numeric(modern_gradient[i,"lat"])
    modern_gradient[i,"mean_Tmin"]<-mean(Iberian_background[which(Iberian_background$lon>=(sitelon-0.167)&
                                                                      Iberian_background$lon<=(sitelon+0.167)&
                                                                      Iberian_background$lat>=(sitelat-0.167)&
                                                                      Iberian_background$lat<=(sitelat+0.167)),"Tmin"])
    modern_gradient[i,"mean_gdd"]<-mean(Iberian_background[which(Iberian_background$lon>=(sitelon-0.167)&
                                                                     Iberian_background$lon<=(sitelon+0.167)&
                                                                     Iberian_background$lat>=(sitelat-0.167)&
                                                                     Iberian_background$lat<=(sitelat+0.167)),"gdd"])
    modern_gradient[i,"mean_alpha"]<-mean(Iberian_background[which(Iberian_background$lon>=(sitelon-0.167)&
                                                                       Iberian_background$lon<=(sitelon+0.167)&
                                                                       Iberian_background$lat>=(sitelat-0.167)&
                                                                       Iberian_background$lat<=(sitelat+0.167)),"alpha"])
    modern_gradient[i,"mean_Tmax"]<-mean(Iberian_background[which(Iberian_background$lon>=(sitelon-0.167)&
                                                                    Iberian_background$lon<=(sitelon+0.167)&
                                                                    Iberian_background$lat>=(sitelat-0.167)&
                                                                    Iberian_background$lat<=(sitelat+0.167)),"Tmax"])
    
    modern_gradient[i,c("mean_insol_summer","mean_insol_winter")]<-get_insolation(sitelat,age_start =0,age_end=0, interval=0 )[,c("insol_summer","insol_winter")]
    
  }, error=function(e){})  
}
#combine
str(modern_gradient)
colnames(modern_gradient)<-colnames(gradient_env)
for(j in c(1,3:ncol(modern_gradient))){
  modern_gradient[,j]<-as.numeric(modern_gradient[,j])
}

gradient_env<-rbind.data.frame(modern_gradient,gradient_env)
str(gradient_env)
for(j in c(1,3:ncol(gradient_env))){
  gradient_env[,j]<-as.numeric(gradient_env[,j])
}
gradient_env$agef<-factor(gradient_env$age,labels=c("0 ka","0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                    "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))


#get anomaly to 0 ka
gradient_env$Tmin_anomaly<-gradient_env$mean_Tmin-rep(modern_gradient[,"mean_Tmin"],nlevels(gradient_env$agef))
gradient_env$gdd_anomaly<-gradient_env$mean_gdd-rep(modern_gradient[,"mean_gdd"],nlevels(gradient_env$agef))
gradient_env$alpha_anomaly<-gradient_env$mean_alpha-rep(modern_gradient[,"mean_alpha"],nlevels(gradient_env$agef))
gradient_env$Tmax_anomaly<-gradient_env$mean_Tmax-rep(modern_gradient[,"mean_Tmax"],nlevels(gradient_env$agef))

gradient_env$insol_summer_anomaly<-gradient_env$mean_insol_summer-rep(modern_gradient[,"mean_insol_summer"],nlevels(gradient_env$agef))
gradient_env$insol_winter_anomaly<-gradient_env$mean_insol_winter-rep(modern_gradient[,"mean_insol_winter"],nlevels(gradient_env$agef))

gradient_env$Tmin_anomaly_0.5ka<-gradient_env$mean_Tmin-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_Tmin"],nlevels(gradient_env$agef)))
gradient_env$gdd_anomaly_0.5ka<-gradient_env$mean_gdd-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_gdd"],nlevels(gradient_env$agef)))
gradient_env$alpha_anomaly_0.5ka<-gradient_env$mean_alpha-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_alpha"],nlevels(gradient_env$agef)))
gradient_env$Tmax_anomaly_0.5ka<-gradient_env$mean_Tmax-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_Tmax"],nlevels(gradient_env$agef)))

gradient_env$insol_summer_anomaly_0.5ka<-gradient_env$mean_insol_summer-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_insol_summer"],nlevels(gradient_env$agef)))
gradient_env$insol_winter_anomaly_0.5ka<-gradient_env$mean_insol_winter-as.numeric(rep(get_gradient(plotdata,500,500)[,"mean_insol_winter"],nlevels(gradient_env$agef)))

gradient_env$elv_label<-NA
for(j in 1:nrow(gradient_env)){
  elv<-as.numeric(gradient_env[j,"elv"])
  if(elv<=1000){
    gradient_env[j,"elv_label"]<-"low"
  }else if(elv>1000){
    gradient_env[j,"elv_label"]<-"high"
  }
}
gradient_env$elv_label<-factor(gradient_env$elv_label,levels = c("low","high"))
write.csv(gradient_env,paste(wd,"/Master Project/Data/Output data/Core/Iberian/gradient_env.csv",sep=""))

########################################################################################################
#####################################   "Hovemoller" type plots  #######################################
########################################################################################################
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))
#"Hovemoller" type plots showing the anomalies in MTCO, MTWA and alpha relative to 0 ka; Figure 2
#anomaly to lon
summary(gradient_env$Tmin_anomaly)
summary(gradient_env$Tmax_anomaly)
summary(gradient_env$alpha_anomaly)

#Make sure 0="whit"
#min+(max-min)*(4/6)=0, so 6min+4max-4min=0, so 4max+2min=0, so min=-2max
p1<-ggplot(data=gradient_env[which(gradient_env$age!=0),],aes(x=fct_reorder(site, lon),y=age,color=Tmin_anomaly,shape=elv_label))+
  geom_point(size=2.5)+theme_dark()+  
  scale_y_continuous(trans = "reverse",breaks=1000*seq(0.5,11.5),labels=seq(0.5,11.5))+
  scale_colour_gradientn(colours = c("darkblue","blue","dodgerblue2","lightskyblue","white","orange","orangered"),
                         limits=c(-24,12),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="MTCO (°C)")+labs(y="Age (kyr BP)")+
  annotate("text", y= 0, x =3,label="(a)",size=5)+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')
#min+(max-min)*(3/6)=0
p2<-ggplot(data=gradient_env[which(gradient_env$age!=0),],aes(x=fct_reorder(site, lon),y=age,color=Tmax_anomaly,shape=elv_label))+
  geom_point(size=2.5)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(0.5,11.5),labels=seq(0.5,11.5))+
  scale_colour_gradientn(colours = c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red"),
                         limits=c(-9.9,9.9),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="MTWA (°C)")+labs(y="Age (kyr BP)")+
  annotate("text", y= 0, x =3,label="(b)",size=5)+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')
#min+(max-min)*(3/6)=0
p3<-ggplot(data=gradient_env[which(gradient_env$age!=0),],aes(x=fct_reorder(site, lon),y=age,color=alpha_anomaly,shape=elv_label))+
  geom_point(size=2.5)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(0.5,11.5),labels=seq(0.5,11.5))+
  scale_colour_gradientn(colours = rev(c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red")),
                         limits=c(-0.5,0.5),na.value = "transparent")+
  theme(axis.text.x = element_blank())+
  labs(x="site in the order of longitude")+
  labs(color=expression(alpha))+labs(y="Age (kyr BP)")+
  annotate("text", y= 0, x =3,label="(c)",size=5)+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')
p<-ggarrange(p1,p2,p3,ncol=1)

ggsave(file="anomaly to lon.jpeg",p,width=11,height=7)

########################################################################################################
#####################################   Time series  Figure 3 ##########################################
########################################################################################################
age_ka<-c(0,seq(0.5,11.5))
mean_site_env<-as.data.frame(age_ka)
mean_site_env$mean_alpha_anomaly<-aggregate(gradient_env$alpha_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$sd_alpha_anomaly<-aggregate(gradient_env$alpha_anomaly~gradient_env$age, FUN=sd)[,2]
mean_site_env$mean_Tmin_anomaly<-aggregate(gradient_env$Tmin_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$sd_Tmin_anomaly<-aggregate(gradient_env$Tmin_anomaly~gradient_env$age, FUN=sd)[,2]
mean_site_env$mean_Tmax_anomaly<-aggregate(gradient_env$Tmax_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$sd_Tmax_anomaly<-aggregate(gradient_env$Tmax_anomaly~gradient_env$age, FUN=sd)[,2]
mean_site_env$mean_insol_summer_anomaly<-aggregate(gradient_env$insol_summer_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$sd_insol_summer_anomaly<-aggregate(gradient_env$insol_summer_anomaly~gradient_env$age, FUN=sd)[,2]
mean_site_env$mean_insol_winter_anomaly<-aggregate(gradient_env$insol_winter_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$sd_insol_winter_anomaly<-aggregate(gradient_env$insol_winter_anomaly~gradient_env$age, FUN=sd)[,2]

setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))

#Time series
p1<-ggplot()+theme_bw()+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly,
                                         ymin=mean_Tmin_anomaly-sd_Tmin_anomaly,
                                         ymax=mean_Tmin_anomaly+sd_Tmin_anomaly))+
  labs(y="MTCO anomaly (°C)")+scale_x_continuous(breaks= seq(0.5,11.5))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  annotate("text", y= max(mean_site_env$mean_Tmin_anomaly+mean_site_env$sd_Tmin_anomaly,na.rm=T), x =0,label="(a)",size=5)

p2<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly,
                                         ymin=mean_Tmax_anomaly-sd_Tmax_anomaly,
                                         ymax=mean_Tmax_anomaly+sd_Tmax_anomaly))+
  labs(y="MTWA anomaly (°C)")+scale_x_continuous(breaks= seq(0.5,11.5))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  annotate("text", y= max(mean_site_env$mean_Tmax_anomaly+mean_site_env$sd_Tmax_anomaly,na.rm=T), x =0,label="(b)",size=5)

p3<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_alpha_anomaly,
                                         ymin=mean_alpha_anomaly-sd_alpha_anomaly,
                                         ymax=mean_alpha_anomaly+sd_alpha_anomaly))+
  labs(y=expression(alpha* " anomaly"),x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5))+
  annotate("text", y= max(mean_site_env$mean_alpha_anomaly+mean_site_env$sd_alpha_anomaly,na.rm=T), x =0,label="(c)",size=5)

p4<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly,
                                         ymin=mean_insol_summer_anomaly-sd_insol_summer_anomaly,
                                         ymax=mean_insol_summer_anomaly+sd_insol_summer_anomaly))+
  labs(y=bquote('Summer insolation anomaly (W'~m^-2~')'),x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5))+
  annotate("text", y= max(gradient_env$insol_summer_anomaly,na.rm=T), x =-1,label="(e)",size=5)

p5<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly,
                                         ymin=mean_insol_winter_anomaly-sd_insol_winter_anomaly,
                                         ymax=mean_insol_winter_anomaly+sd_insol_winter_anomaly))+
  labs(y=bquote('Winter insolation anmaly (W'~m^-2~')'))+
  scale_x_continuous(breaks= seq(0.5,11.5))+  
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  annotate("text", y= max(gradient_env$insol_winter_anomaly,na.rm=T), x =-1,label="(d)",size=5)

p<-ggarrange(p1,p5,p2,p4,p3,ncol=2)
ggsave(p,file="Time series.jpeg",width=10,height=10)

#######################################################################################################
##############################   Anomaly and longitude and elevation   ################################
#######################################################################################################
gradient_env$elv_km<-gradient_env$elv/1000
#if(!require(ggpubr)){ install.packages("ggpubr");library(ggpubr)}
#Tmin
df_model_Tmin_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,Tmin_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmin_lon<-rbind.data.frame(df_model_Tmin_lon,dat)
}
colnames(df_model_Tmin_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmin_lon[which(df_model_Tmin_lon$p.slope<0.01),"age"]

#Tmax
df_model_Tmax_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,Tmax_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmax_lon<-rbind.data.frame(df_model_Tmax_lon,dat)
}
colnames(df_model_Tmax_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmax_lon[which(df_model_Tmax_lon$p.slope<0.01),"age"]

#alpha
df_model_alpha_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,alpha_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=3),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_alpha_lon<-rbind.data.frame(df_model_alpha_lon,dat)
}
colnames(df_model_alpha_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_alpha_lon[which(df_model_alpha_lon$p.slope<0.01),"age"]

#Tmin
df_model_Tmin_elv<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,Tmin_anomaly~elv_km))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["elv_km","Estimate"],digits=2),
                        round(m[["coefficients"]]["elv_km","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["elv_km","Estimate"],digits=2))
  df_model_Tmin_elv<-rbind.data.frame(df_model_Tmin_elv,dat)
}
colnames(df_model_Tmin_elv)<-c("age","slope_elv","p.slope_elv","x0_elv")
df_model_Tmin_elv[which(df_model_Tmin_elv$p.slope<0.01),"age"]

#Tmax
df_model_Tmax_elv<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,Tmax_anomaly~elv_km))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["elv_km","Estimate"],digits=2),
                        round(m[["coefficients"]]["elv_km","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["elv_km","Estimate"],digits=2))
  df_model_Tmax_elv<-rbind.data.frame(df_model_Tmax_elv,dat)
}
colnames(df_model_Tmax_elv)<-c("age","slope_elv","p.slope_elv","x0_elv")
df_model_Tmax_elv[which(df_model_Tmax_elv$p.slope<0.01),"age"]

#alpha
df_model_alpha_elv<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env[which(gradient_env$age==age),]
  m<-summary(lm(data=sub,alpha_anomaly~elv_km))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["elv_km","Estimate"],digits=3),
                        round(m[["coefficients"]]["elv_km","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["elv_km","Estimate"],digits=2))
  df_model_alpha_elv<-rbind.data.frame(df_model_alpha_elv,dat)
}
colnames(df_model_alpha_elv)<-c("age","slope_elv","p.slope_elv","x0_elv")
df_model_alpha_elv[which(df_model_alpha_elv$p.slope<0.01),"age"]

#combine
df_model_lon<-rbind.data.frame(df_model_Tmin_lon,df_model_Tmax_lon,df_model_alpha_lon)
df_model_elv<-rbind.data.frame(df_model_Tmin_elv,df_model_Tmax_elv,df_model_alpha_elv)
df_model<-cbind.data.frame(df_model_lon,df_model_elv[,-1])
df_model$age<-df_model$age/1000
write.csv(df_model,paste(wd,"/Master Project/Data/Output data/Core/Iberian/df_model.csv",sep=""))

############ plot
plot_gradient<-gradient_env[which(gradient_env$age!=0),]
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="3.5 ka"|age=="4.5 ka"|age=="5.5 ka"|age=="6.5 ka"){
    plot_gradient[i,"label_alpha_lon"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_alpha_lon"]<-as.character(plot_gradient[i,"agef"])
  }
}
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="3.5 ka"|age=="4.5 ka"|age=="5.5 ka"|age=="6.5 ka"){
    plot_gradient[i,"label_Tmax_lon"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_Tmax_lon"]<-as.character(plot_gradient[i,"agef"])
  }
}
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="6.5 ka"|age=="7.5 ka"){
    plot_gradient[i,"label_alpha_elv"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_alpha_elv"]<-as.character(plot_gradient[i,"agef"])
  }
}
############ Figure 5 alpha to lon
dat_text <- data.frame(
  label = plot_gradient$label_alpha_lon,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,alpha_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Longitude (°E)",y=expression(alpha~"anomaly to 0 ka"))+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y = 0.33, label = label))
ggsave(file="alpha to lon.jpeg",p,width=8,height=8)

############ Figure 6 alpha to elv
dat_text <- data.frame(
  label = plot_gradient$label_alpha_elv,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,alpha_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Elevation (km)",y=expression(alpha~"anomaly to 0 ka"))+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = 0.33, label = label))
ggsave(file="alpha to elv.jpeg",p,width=8,height=8)

############ Figure S1 Tmax to lon
dat_text <- data.frame(
  label = plot_gradient$label_Tmax_lon,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,Tmax_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Longitude (°E)",y="MTWA (°C) anomaly to 0 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y = 7, label = label))
ggsave(file="Tmax to lon.jpeg",p,width=8,height=8)

############ Figure S2 Tmax to elv
dat_text <- data.frame(
  label = plot_gradient$agef,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,Tmax_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Elevation (km)",y="MTWA (°C) anomaly to 0 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = 7, label = label))
ggsave(file="Tmax to elv.jpeg",p,width=8,height=8)


############ Figure S3 Tmin to lon
dat_text <- data.frame(
  label = plot_gradient$agef,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,Tmin_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Longitude (°E)",y="MTCO (°C) anomaly to 0 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y =-20, label = label))
ggsave(file="Tmin to lon.jpeg",p,width=8,height=8)

############ Figure S4 Tmin to elv
dat_text <- data.frame(
  label = plot_gradient$agef,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,Tmin_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",se=F,col="red")+labs(x="Elevation (km)",y="MTCO (°C) anomaly to 0 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = -20, label = label))
ggsave(file="Tmin to elv.jpeg",p,width=8,height=8)


