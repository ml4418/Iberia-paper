rm(list = ls())
Sys.setenv(LANG = "en")

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
install.packages("remotes")
remotes::install_github("special-uor/fxTWAPLS@v0.1.0")
libary(fxTWAPLS)
source(paste(wd,"/Master Project/Script/Reconstruction paper/Iberian Peninsula reconstruction functions.R",sep=""))

modern_pollen<- read.csv(paste(wd,"/Master Project/Data/Input data/Modern_Pollen_gdd_alpha_Tmin.csv",sep=""), row.names=1)

taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

for(i in 1:nrow(modern_pollen)){
  #i=i+1
  modern_pollen[i,"Tmax"]<-get_MTWA(MTCO=modern_pollen[i,"Tmin"],GDD0=modern_pollen[i,"gdd"])
}
summary(modern_pollen$Tmin);
summary(modern_pollen$Tmax);
summary(modern_pollen$alpha);

# Training ---------------------------------------------

#compare training results
fit_tf1_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf2_Tmin <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf1_Tmin_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_tf2_Tmin_pspline <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

fit_tf1_Tmax <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf2_Tmax <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.02)
fit_tf1_Tmax_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_tf2_Tmax_pspline <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

fit_tf1_alpha <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.002)
fit_tf2_alpha <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="bin",bin=0.002)
fit_tf1_alpha_pspline <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)
fit_tf2_alpha_pspline <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)

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
    labs(y="fxTWA-PLS1 MTCO (¡ãC)",x="MTCO (¡ãC)")
  p_Tmin_tf2<-ggplot(data=train,aes(Tmin,Tmin_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmin,max_Tmin)+xlim(min_Tmin,max_Tmin)+
    annotate("text", y= max_Tmin, x = min_Tmin,label="(b)")+
    labs(y="fxTWA-PLS2 MTCO (¡ãC)",x="MTCO (¡ãC)")
  
  
  p_Tmax_tf1<-ggplot(data=train,aes(Tmax,Tmax_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmax,max_Tmax)+xlim(min_Tmax,max_Tmax)+
    annotate("text", y= max_Tmax, x = min_Tmax,label="(c)")+
    labs(y="fxTWA-PLS1 MTWA (¡ãC)",x="MTWA (¡ãC)")
  p_Tmax_tf2<-ggplot(data=train,aes(Tmax,Tmax_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+ylim(min_Tmax,max_Tmax)+xlim(min_Tmax,max_Tmax)+
    annotate("text", y= max_Tmax, x = min_Tmax,label="(d)")+
    labs(y="fxTWA-PLS2 MTWA (¡ãC)",x="MTWA (¡ãC)")
  
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
plot.resid.sig<-function(name,fit_tf1_Tmin,fit_tf2_Tmin,nsig_tf1_Tmin,nsig_tf2_Tmin,
                         fit_tf1_Tmax,fit_tf2_Tmax,nsig_tf1_Tmax,nsig_tf2_Tmax,
                         fit_tf1_alpha,fit_tf2_alpha,nsig_tf1_alpha,nsig_tf2_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #Get the residuals
  resid<-cbind.data.frame(fit_tf1_Tmin[["x"]],
                          fit_tf1_Tmin[["fit"]][,nsig_tf1_Tmin]-fit_tf1_Tmin[["x"]],
                          fit_tf2_Tmin[["fit"]][,nsig_tf2_Tmin]-fit_tf1_Tmin[["x"]],
                          fit_tf1_Tmax[["x"]],
                          fit_tf1_Tmax[["fit"]][,nsig_tf1_Tmax]-fit_tf1_Tmax[["x"]],
                          fit_tf2_Tmax[["fit"]][,nsig_tf2_Tmax]-fit_tf1_Tmax[["x"]],
                          fit_tf1_alpha[["x"]],
                          fit_tf1_alpha[["fit"]][,nsig_tf1_alpha]-fit_tf1_alpha[["x"]],
                          fit_tf2_alpha[["fit"]][,nsig_tf2_alpha]-fit_tf1_alpha[["x"]])
  
  colnames(resid)<-c("Tmin","Tmin_tf1","Tmin_tf2","Tmax","Tmax_tf1","Tmax_tf2","alpha","alpha_tf1","alpha_tf2")
  max_y_Tmin<-max(resid[,c("Tmin_tf1","Tmin_tf2")])
  min_y_Tmin<-min(resid[,c("Tmin_tf1","Tmin_tf2")])
  max_y_Tmax<-max(resid[,c("Tmax_tf1","Tmax_tf2")])
  min_y_Tmax<-min(resid[,c("Tmax_tf1","Tmax_tf2")])
  max_y_alpha<-max(resid[,c("alpha_tf1","alpha_tf2")])
  min_y_alpha<-min(resid[,c("alpha_tf1","alpha_tf2")])
  
  p_Tmin_tf1<-ggplot(data=resid,aes(Tmin,Tmin_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_Tmin,max_y_Tmin)+
    annotate("text", y= max_y_Tmin, x = min(resid[,"Tmin"]),label="(a)")+
    labs(y="fxTWA-PLS1 MTCO (¡ãC)",x="MTCO (¡ãC)")
  p_Tmin_tf2<-ggplot(data=resid,aes(Tmin,Tmin_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_Tmin,max_y_Tmin)+
    annotate("text", y= max_y_Tmin, x = min(resid[,"Tmin"]),label="(b)")+
    labs(y="fxTWA-PLS2 MTCO (¡ãC)",x="MTCO (¡ãC)")
  
  
  p_Tmax_tf1<-ggplot(data=resid,aes(Tmax,Tmax_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_Tmax,max_y_Tmax)+
    annotate("text", y= max_y_Tmax, x = min(resid[,"Tmax"]),label="(c)")+
    labs(y="fxTWA-PLS1 MTWA (¡ãC)",x="MTWA (¡ãC)")
  p_Tmax_tf2<-ggplot(data=resid,aes(Tmax,Tmax_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_Tmax,max_y_Tmax)+
    annotate("text", y= max_y_Tmax, x = min(resid[,"Tmax"]),label="(d)")+
    labs(y="fxTWA-PLS2 MTWA (¡ãC)",x="MTWA (¡ãC)")
  
  p_alpha_tf1<-ggplot(data=resid,aes(alpha,alpha_tf1))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_alpha,max_y_alpha)+
    annotate("text", y= max_y_alpha, x = min(resid[,"alpha"]),label="(e)")+
    labs(y= expression("fxTWA-PLS1   "*alpha), x = expression(alpha))
  
  p_alpha_tf2<-ggplot(data=resid,aes(alpha,alpha_tf2))+geom_point(size=0.8)+theme_bw()+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess', formula= y~x,color='red', se = FALSE)+ylim(min_y_alpha,max_y_alpha)+
    annotate("text", y= max_y_alpha, x = min(resid[,"alpha"]),label="(f)")+
    labs(y= expression("fxTWA-PLS2   "*alpha), x = expression(alpha))
  
  p<-ggarrange(p_Tmin_tf1,p_Tmin_tf2,p_Tmax_tf1,p_Tmax_tf2,p_alpha_tf1,p_alpha_tf2,  ncol = 2)
  ggsave(file=paste(name,"residual results.jpeg"),p,width=8,height=9)
  
}
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))

plot.train.sig("First and second version training results",
               fit_tf1_Tmin,fit_tf2_Tmin_pspline,4,4,
               fit_tf1_Tmax,fit_tf2_Tmax_pspline,2,4,
               fit_tf1_alpha,fit_tf2_alpha_pspline,3,3)

plot.resid.sig("First and second version residuals",
               fit_tf1_Tmin,fit_tf2_Tmin_pspline,4,4,
               fit_tf1_Tmax,fit_tf2_Tmax_pspline,2,4,
               fit_tf1_alpha,fit_tf2_alpha_pspline,3,3)
# In step 7 of training, climate variable is regressed to the components obtained,
# MTCO
fit_tf_Tmin <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

# GDD0
fit_tf_gdd <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$gdd, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=20)

# alpha
fit_tf_alpha <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)

#MTWA
fit_tf_Tmax <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmax, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)


########################################################################################################
####################################  Cross validation  Table 1 ########################################
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
`%>%` <- magrittr::`%>%`
if(!require(foreach)){install.packages("foreach");library(foreach)}
cv_tf2_Tmin_pspline <- fxTWAPLS::cv.pr.w(taxa,
                                         modern_pollen$Tmin,
                                         nPLS = 5,
                                         fxTWAPLS::TWAPLS.w2,
                                         fxTWAPLS::TWAPLS.predict.w,
                                         pseudo_Tmin,
                                         usefx = TRUE,
                                         fx_method = "pspline",
                                         bin = 0.02,
                                         cpus = 8,
                                         test_mode = F)  %>% fxTWAPLS::pb()  
write.csv(cv_tf2_Tmin_pspline, "cv_tf2_Tmin_pspline.csv")

cv_tf2_gdd_pspline <- fxTWAPLS::cv.pr.w(taxa,
                                        modern_pollen$gdd,
                                        nPLS = 5,
                                        fxTWAPLS::TWAPLS.w2,
                                        fxTWAPLS::TWAPLS.predict.w,
                                        pseudo_gdd,
                                        usefx = TRUE,
                                        fx_method = "pspline",
                                        bin = 20,
                                        cpus = 8,
                                        test_mode = F)   %>% fxTWAPLS::pb()  
write.csv(cv_tf2_gdd_pspline, "cv_tf2_gdd_pspline.csv")

cv_tf2_alpha_pspline <- fxTWAPLS::cv.pr.w(taxa,
                                          modern_pollen$alpha,
                                          nPLS = 5,
                                          fxTWAPLS::TWAPLS.w2,
                                          fxTWAPLS::TWAPLS.predict.w,
                                          pseudo_alpha,
                                          usefx = TRUE,
                                          fx_method = "pspline",
                                          bin = 0.002,
                                          cpus = 8,
                                          test_mode = F)   %>% fxTWAPLS::pb()  
write.csv(cv_tf2_alpha_pspline, "cv_tf2_alpha_pspline.csv")

cv_tf2_Tmax_pspline <- fxTWAPLS::cv.pr.w(taxa,
                                         modern_pollen$Tmax,
                                         nPLS = 5,
                                         fxTWAPLS::TWAPLS.w2,
                                         fxTWAPLS::TWAPLS.predict.w,
                                         pseudo_Tmax,
                                         usefx = TRUE,
                                         fx_method = "pspline",
                                         bin = 0.02,
                                         cpus = 8,
                                         test_mode = F)   %>% fxTWAPLS::pb()  
write.csv(cv_tf2_Tmax_pspline, "cv_tf2_Tmax_pspline.csv")

rand_tf2_Tmin_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmin_pspline, n.perm = 999) 
rand_tf2_alpha_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_alpha_pspline, n.perm = 999) 
rand_tf2_Tmax_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_Tmax_pspline, n.perm = 999) 
rand_tf2_gdd_pspline <- fxTWAPLS::rand.t.test.w(cv_tf2_gdd_pspline, n.perm = 999) 

rand_tf2_pspline<-rbind.data.frame(rand_tf2_Tmin_pspline,rand_tf2_Tmax_pspline,rand_tf2_alpha_pspline,rand_tf2_gdd_pspline)
write.csv(rand_tf2_pspline, "rand_tf2_pspline.csv")


########################################################################################################
####################################  Reconstruction  ##################################################
########################################################################################################
library(readxl)
Iberian_all <- read_excel(paste(wd,"/Master Project/Data/Input data/Iberia_pollen_records_v3_0307.xlsx",sep=""))

Iberian<-Iberian_all[which(Iberian_all$`INTCAL2020_median`<=12000),]
#Age
colnames(Iberian);str(Iberian)
summary(Iberian[,"INTCAL2020_median"])
#Exclude samples with large age uncertainty
summary(abs(Iberian[,"INTCAL2020_uncert_95"]-Iberian[,"INTCAL2020_uncert_5"])/(2*1.645))
which_age_keep<-which(abs((Iberian[,"INTCAL2020_uncert_95"]-Iberian[,"INTCAL2020_uncert_5"])/(2*1.645))<=100)
Iberian<-Iberian[which_age_keep,]
write.csv(Iberian,paste(wd,"/Master Project/Data/Input data/Iberia.csv",sep=""))

#get average resolution
resolution<-data.frame()
Iberian1<-Iberian
Iberian1$entity_name<-as.factor(Iberian1$entity_name)
for(k in 1:nlevels(Iberian1$entity_name)){
  eachsite<-Iberian1[which(Iberian1$entity_name==levels(Iberian1$entity_name)[k]),]
  resolution<-rbind.data.frame(resolution,eachsite[-1,"INTCAL2020_median"]-eachsite[-nrow(eachsite),"INTCAL2020_median"])
  print(c(k,summary(resolution$INTCAL2020_median)))
}
mean(resolution$INTCAL2020_median)
summary(resolution$INTCAL2020_median)

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
core_sig<-cbind.data.frame(Iberian[,c("site_name","latitude","longitude","elevation","INTCAL2020_median")],
                           fossil_tf_Tmin[["fit"]][,4],fossil_tf_gdd[["fit"]][,2],fossil_tf_alpha[["fit"]][,3], fossil_tf_Tmax[["fit"]][,4])
colnames(core_sig)<-c("site","lat","lon","elv","age","Tmin","gdd","alpha","Tmax")
write.csv(core_sig,paste(wd,"/Master Project/Data/Output data/Core/Iberian/Iberian_core_sig.csv",sep=""))

#Get the sample specific errors, use nboot=1000
`%>%` <- magrittr::`%>%`
#MTCO
sse_tf_Tmin<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmin,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w2,predictfun=fxTWAPLS::TWAPLS.predict.w,
                                  nboot=1000,nPLS=5,nsig=4,usefx=TRUE,fx_method = "pspline",bin=0.02)%>%fxTWAPLS::pb()
#GDD0
sse_tf_gdd<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$gdd,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w2,predictfun=fxTWAPLS::TWAPLS.predict.w,
                                 nboot=1000,nPLS=5,nsig=2,usefx=TRUE,fx_method = "pspline",bin=20)%>%fxTWAPLS::pb()
#alpha
sse_tf_alpha<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$alpha,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w2,predictfun=fxTWAPLS::TWAPLS.predict.w,
                                   nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx_method = "pspline",bin=0.002)%>%fxTWAPLS::pb()
#MTWA
sse_tf_Tmax<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmax,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w2,predictfun=fxTWAPLS::TWAPLS.predict.w,
                                  nboot=1000,nPLS=5,nsig=4,usefx=TRUE,fx_method = "pspline",bin=0.02)%>%fxTWAPLS::pb()

sse_core_sig<-cbind.data.frame(sse_tf_Tmin,sse_tf_gdd,sse_tf_alpha,sse_tf_Tmax)
colnames(sse_core_sig)<-c("sse_Tmin","sse_gdd","sse_alpha","sse_Tmax")
write.csv(sse_core_sig,paste(wd,"/Master Project/Data/Output data/Core/Iberian/Iberian_sse_core_sig.csv",sep=""))

########################################################################################################
####################################### Sitelist Table ##############################################
########################################################################################################
#Get information for sites
#sitelist<-Iberian[!duplicated(Iberian$site_name),]
sitename<-unique(Iberian$site_name)
sitelist<-as.data.frame(sitename)
for(i in 1:nrow(sitelist)){
  #i=i+1
  siteinfo<-Iberian[which(Iberian$site_name==sitelist[i,"sitename"]),]
  sitelist[i,c("source","lon","lat","elv")]<-unique(siteinfo[,c("source","longitude","latitude","elevation")])
  sitelist[i,c("entity")]<-paste(unique(siteinfo[,c("entity_name")]))
  
  sitelist[i,"start_age"]<-max(siteinfo[,"INTCAL2020_median"])
  sitelist[i,"end_age"]<-min(siteinfo[,"INTCAL2020_median"])
  
  sitelist[i,"length_yr"]<-max(siteinfo[,"INTCAL2020_median"])-min(siteinfo[,"INTCAL2020_median"])
  sitelist[i,"n_sample"]<-sum(!is.na(siteinfo[,"INTCAL2020_median"]))
  
  ref<-unique(siteinfo[,"reference"])
  if(all(ref=="")){
    sitelist[i,"ref"]<-""
  }else{
    sitelist[i,"ref"]<-paste(ref[which(ref!=""),])
  }
  
}
sitelist<-sitelist[order(sitelist$sitename), ]
sitelist$lon<-round(sitelist$lon,digits=2)
sitelist$lat<-round(sitelist$lat,digits=2)
sitelist$elv<-round(sitelist$elv,digits=0)
sitelist$start_age<-round(sitelist$start_age,digits=0)
sitelist$end_age<-round(sitelist$end_age,digits=0)
sitelist$length_yr<-round(sitelist$length_yr,digits=0)
rownames(sitelist)<-seq(1:nrow(sitelist))

#sitelist information from Yicheng Shen
date_control <- read_excel(paste(wd,"/Master Project/Data/Input data/Iberia_data_dates_v3.xlsx",sep=""))
date_control[which(date_control$`Site name`=="Las Pardillas Lake"),"Site name"]<-"Las Pardillas"
sitelist[which(!(sitelist$sitename%in%date_control$`Site name`)),"sitename"]
date_control[which(!(date_control$`Site name`%in%sitelist$sitename)),"Site name"]
date_control<-date_control[-which(!(date_control$`Site name`%in%sitelist$sitename)),]
setequal(sitelist$sitename,date_control$`Site name`)
date_control<-date_control[order(date_control$`Site name`), ]

#save them
if(!require(writexl)){ install.packages("writexl");library(writexl)}
sitelist$n_dating<-table(date_control$`Site name`)
#reorder
sitelist<-sitelist[,c("sitename","entity","lon","lat","elv", "start_age", "end_age","length_yr","n_sample","n_dating",  "source","ref")]
#check
sum(sitelist$source=="author");sum(sitelist$source=="EPD");sum(sitelist$source=="PANGAEA")
sum(sitelist$source=="author")+sum(sitelist$source=="EPD")+sum(sitelist$source=="PANGAEA")
sum(sitelist$n_sample)
write_xlsx(sitelist,paste(wd,"/Master Project/Data/Output data/Core/Iberian/sitelist.xlsx",sep=""))
write_xlsx(date_control,paste(wd,"/Master Project/Data/Output data/Core/Iberian/date_control.xlsx",sep=""))

########################################################################################################
#######################################    CCA  Table  ###############################################
########################################################################################################
if(!require(vegan)){ install.packages("vegan");library(vegan)}
ord<-cca(formula=taxa~modern_pollen$Tmin+modern_pollen$Tmax+modern_pollen$alpha)
ord
round(vif.cca(ord),digits=2)
round(ord[["CCA"]][["biplot"]],digits = 3)
anova(ord)
anova(ord, by="term", permutations=999)
anova(ord, by="axis", permutations=999)
plot(ord,display=c("species","bp"))
text(ord, "species", col="black", cex=0.8)



ord<-cca(formula=core~core_sig$Tmin+core_sig$Tmax+core_sig$alpha,na.action=na.omit)
ord
round(vif.cca(ord),digits=2)
round(ord[["CCA"]][["biplot"]],digits = 3)
anova(ord)
anova(ord, by="term", permutations=999)
anova(ord, by="axis", permutations=999)
plot(ord,display=c("bp"),xlim=c(-3.5,2.5),ylim=c(-3,4.5))
plot(ord,display=c("bp","sp"))
plot(ord,display=c("bp"))
points(ord,display=c("bp"),col="red")
text(ord, "species", col="black", cex=0.7)

########################################################################################################
##################################   Find out taxa contribution  #######################################
########################################################################################################
#the occurence of taxa in fossil pollen data
occur<-colSums(core != 0,na.rm=T)
summary(occur)
#u and t are both centered in fxTWAPLS
#Tmin
u_t2_Tmin<-as.data.frame(fit_tf_Tmin[["u"]]/fit_tf_Tmin[["t"]]^2)
para_Tmin<-fit_tf_Tmin[["alpha"]]
u_t2_Tmin$taxa<-colnames(taxa)
u_t2_Tmin$occur<-occur
u_t2_Tmin_keep<-u_t2_Tmin[which(u_t2_Tmin$occur>=median(occur)),]
u_t2_Tmin_keep<-u_t2_Tmin_keep[order(u_t2_Tmin_keep$V1),];
u_t2_Tmin_keep[1:10,"taxa"]
u_t2_Tmin_keep[(nrow(u_t2_Tmin_keep)-9):nrow(u_t2_Tmin_keep),"taxa"]

#Tmax
u_t2_Tmax<-as.data.frame(fit_tf_Tmax[["u"]]/fit_tf_Tmax[["t"]]^2)
para_Tmax<-fit_tf_Tmax[["alpha"]]
u_t2_Tmax$taxa<-colnames(taxa)
u_t2_Tmax$occur<-occur
u_t2_Tmax_keep<-u_t2_Tmax[which(u_t2_Tmax$occur>=median(occur)),]
u_t2_Tmax_keep<-u_t2_Tmax_keep[order(u_t2_Tmax_keep$V1),]
u_t2_Tmax_keep[1:10,"taxa"]
u_t2_Tmax_keep[(nrow(u_t2_Tmax_keep)-9):nrow(u_t2_Tmax_keep),"taxa"]

#alph
u_t2_alpha<-as.data.frame(fit_tf_alpha[["u"]]/fit_tf_alpha[["t"]]^2)
para_alpha<-fit_tf_alpha[["alpha"]]
u_t2_alpha$taxa<-colnames(taxa)
u_t2_alpha$occur<-occur
u_t2_alpha_keep<-u_t2_alpha[which(u_t2_alpha$occur>=median(occur)),]
u_t2_alpha_keep<-u_t2_alpha_keep[rev(order(u_t2_alpha_keep$V1)),]
u_t2_alpha_keep[1:10,"taxa"]
u_t2_alpha_keep[(nrow(u_t2_alpha_keep)-9):nrow(u_t2_alpha_keep),"taxa"]

contri_top10_up<-cbind.data.frame(u_t2_Tmin_keep[1:10,"taxa"],
                                  u_t2_Tmax_keep[1:10,"taxa"],
                                  u_t2_alpha_keep[1:10,"taxa"])
contri_top10_down<-cbind.data.frame(u_t2_Tmin_keep[(nrow(u_t2_Tmin_keep)-9):nrow(u_t2_Tmin_keep),"taxa"],
                                    u_t2_Tmax_keep[(nrow(u_t2_Tmax_keep)-9):nrow(u_t2_Tmax_keep),"taxa"],
                                    u_t2_alpha_keep[(nrow(u_t2_alpha_keep)-9):nrow(u_t2_alpha_keep),"taxa"])
write.csv(contri_top10_up,paste(wd,"/Master Project/Data/Output data/Core/Iberian/contri_top10_up.csv",sep=""))
write.csv(contri_top10_down,paste(wd,"/Master Project/Data/Output data/Core/Iberian/contri_top10_down.csv",sep=""))

########################################################################################################
###########################    Calculate insolation at each site  #######################################
########################################################################################################
plotdata<-core_sig

#exclude outliers of the climate variables
plotdata<-plotdata[which(plotdata$alpha>=0&plotdata$alpha<=1.26),]

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
unique(plotdata$site)
nrow(plotdata)
########################################################################################################
###########################   Relationship between MTWA and alpha  #####################################
########################################################################################################
#Figure 
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))
summary(plotdata$Tmax);summary(modern_pollen$Tmax)

p1<-ggplot(modern_pollen,aes(alpha,Tmax))+geom_point(size=0.8)+theme_bw()+
  labs(x= expression(alpha), y = "MTWA (¡ãC)")+ylim(0,35)+
  theme(legend.position = "none",
        axis.text=element_text(size=14),axis.title=element_text(size=14))+
  annotate("text", y= 35, x =0,label="(a)",size=6)
p2<-ggplot(plotdata,aes(alpha,Tmax))+geom_point(size=0.8)+theme_bw()+
  labs(x= expression(alpha), y = "MTWA (¡ãC)")+ylim(0,35)+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.text=element_text(size=14),axis.title=element_text(size=14))+
  annotate("text", y= 35, x =0,label="(b)",size=6)
p<-ggarrange(p1,p2,ncol=2)
ggsave(file="Relationship between MTWA and alpha.jpeg",p,width=12,height=6)

########################################################################################################
########################    Maps to show the training data  (Figure S1)  ################################
########################################################################################################
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))
#Modern sites
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}

world <- map_data("world") 
minLong<-min(modern_pollen$Long);maxLong<-max(modern_pollen$Long);
minLat<-min(modern_pollen$Lat);maxLat<-max(modern_pollen$Lat);

region<-world[which(world$long>minLong&world$long<maxLong&world$lat>minLat&world$lat<maxLat),]

xat <- pretty(region$long)
yat <- pretty(region$lat)

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}

p1<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray60',color='gray30') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= Tmin), size = 1)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  scalebar(data=region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  scale_y_continuous(breaks = yat, labels = paste0(yat,'¡ãN')) +
  scale_x_continuous(breaks = xat, labels = paste0(xat,'¡ãE')) +
  theme(axis.text=element_text(size=12),legend.position = "right")+
  labs(x=NULL,y="Latitude",colour ="(a) MTCO (¡ãC)")+
  theme(axis.text.x = element_blank())
p2<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray60',color='gray30') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= Tmax), size = 1)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  scalebar(data=region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  scale_y_continuous(breaks = yat, labels = paste0(yat,'¡ãN')) +
  scale_x_continuous(breaks = xat, labels = paste0(xat,'¡ãE')) +
  theme(axis.text=element_text(size=12),legend.position = "right")+
  labs(x=NULL,y="Latitude",colour ="(b) MTWA (¡ãC)")+
  theme(axis.text.x = element_blank())
p3<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray60',color='gray30') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= alpha), size = 1)+
  scale_colour_gradientn(colours = rev(c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red")),na.value = "grey")+
  scalebar(data=region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  scale_y_continuous(breaks = yat, labels = paste0(yat,'¡ãN')) +
  scale_x_continuous(breaks = xat, labels = paste0(xat,'¡ãE')) +
  theme(axis.text=element_text(size=12),legend.position = "right")+
  labs(x="Longitude",y="Latitude",colour =expression("(c)"~alpha))

p<-ggarrange(p1,p2,p3,ncol=1)
ggsave(file="Map of SMPDS climate.jpeg",p,width=10,height=12)

########################################################################################################
#################################    Maps of the climate space  ########################################
########################################################################################################
library(raster)
library(rgdal) # for spTransform
if(!require(ncdf4)){ install.packages("ncdf4");library(ncdf4)}

#get CRU background for climate space
#from Turner, M., Wei, D., Prentice, I., & Harrison, S. (2021). The impact of methodological decisions on climate reconstructions using WA-PLS. Quaternary Research, 99, 341-356. doi:10.1017/qua.2020.44
CRU_Tmin_gdd<-read.csv(paste(wd,"/Master Project/Data/Input data/CRU background/CRU Europe Sep18.csv",sep=""),row.names = 1)
colnames(CRU_Tmin_gdd)<-c("lat","lon","elv","Tmin","gdd0")
#ggplot(CRU_Tmin_gdd,aes(lon,lat,color=Tmin))+geom_point()
#ggplot(CRU_Tmin_gdd,aes(lon,lat,color=gdd0))+geom_point()

#from Roberto
CRU_MI<-brick(paste(wd,"/Master Project/Data/Input data/CRU background/cru_ts4.04-clim-1961-1990-mi.nc",sep=""))
CRU_MI<-rasterToPoints(CRU_MI, spatial=TRUE)
CRU_MI<-spTransform(CRU_MI, CRS(CRU_MI@proj4string@projargs))
CRU_MI <- data.frame(lon=coordinates(CRU_MI)[,1],lat=coordinates(CRU_MI)[,2],CRU_MI@data)  
CRU_MI<-CRU_MI[which(CRU_MI$lon>=(-21)&CRU_MI$lon<=150&CRU_MI$lat>=29&CRU_MI$lat<=82),]
colnames(CRU_MI)<-c("lon","lat","MI")
#ggplot(CRU_MI,aes(lon,lat,color=MI))+geom_point()

CRU_gdd<-brick(paste(wd,"/Master Project/Data/Input data/CRU background/cru_ts4.04-clim-1961-1990-daily.tmp-gdd0.nc",sep=""))
CRU_gdd<-rasterToPoints(CRU_gdd, spatial=TRUE)
CRU_gdd<-spTransform(CRU_gdd, CRS(CRU_gdd@proj4string@projargs))
CRU_gdd <- data.frame(lon=coordinates(CRU_gdd)[,1],lat=coordinates(CRU_gdd)[,2],CRU_gdd@data)  
CRU_gdd<-CRU_gdd[which(CRU_gdd$lon>=(-21)&CRU_gdd$lon<=150&CRU_gdd$lat>=29&CRU_gdd$lat<=82),]
colnames(CRU_gdd)<-c("lon","lat","gdd")
#ggplot(CRU_gdd,aes(lon,lat,color=gdd))+geom_point()

CRU<-merge(CRU_MI,CRU_gdd,by=c("lon","lat"))
CRU<-merge(CRU,CRU_Tmin_gdd,by=c("lon","lat"))
for(i in 1:nrow(CRU)){
  #i=i+1
  CRU[i,"Tmax"]<-get_MTWA(MTCO=CRU[i,"Tmin"],GDD0=CRU[i,"gdd"])
}
summary(CRU$Tmax);hist(CRU$Tmax)

MI<-CRU$MI
w<-3
fai<-1/MI
F<-1+fai-(1+fai^w)^(1/w)
alpha<-1.26*MI*F
hist(alpha)
plot(alpha~MI)
CRU$alpha<-alpha

summary(lm(CRU$gdd0~CRU$gdd)) #values from Mark and from Roberto are comparable

modern_Iberian_pollen<-modern_pollen[which(modern_pollen$Long>(-10)&modern_pollen$Long<5&modern_pollen$Lat>36&modern_pollen$Lat<44),]
#ggplot(modern_Iberian_pollen,aes(Long,Lat,color=Tmin))+geom_point()
#ggplot(modern_Iberian_pollen,aes(Long,Lat,color=Tmax))+geom_point()
#ggplot(modern_Iberian_pollen,aes(Long,Lat,color=alpha))+geom_point()

p1<-ggplot()+theme_bw()+geom_point(data=CRU,aes(Tmin,Tmax),color="grey",size=1)+
  geom_point(data=modern_pollen,aes(Tmin,Tmax),color="black",size=1)+
  geom_point(data=modern_Iberian_pollen,aes(Tmin,Tmax),color="red",size=1)+
  labs(y="MTWA (¡ãC)",x="MTCO (¡ãC)")+
  annotate("text", y= max(CRU$Tmax), x =min(CRU$Tmin),label="(a)",size=5)
p2<-ggplot()+theme_bw()+geom_point(data=CRU,aes(Tmax,alpha),color="grey",size=1)+
  geom_point(data=modern_pollen,aes(Tmax,alpha),color="black",size=1)+
  geom_point(data=modern_Iberian_pollen,aes(Tmax,alpha),color="red",size=1)+
  labs(y=expression(alpha),x="MTWA (¡ãC)")+
  annotate("text", y= max(CRU$alpha), x =min(CRU$Tmax),label="(b)",size=5)
p<-ggarrange(p1,p2,ncol=2)
ggsave(file="Climate space.jpeg",p,width=8,height=4)

########################################################################################################
#####################################   Map of fossil sites   #########################################
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
summary(Iberian$longitude);summary(Iberian$latitude)
region<-world[which(world$long>(-10)&world$long<5&world$lat>36&world$lat<44),]
Iberian_background<-background[which(background$lon>(-10)&background$lon<5&background$lat>36&background$lat<44),]

#tranfer MI to alpha
MI<-Iberian_background$MI
w<-3
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

sitelist$elv_label<-NA
for(j in 1:nrow(sitelist)){
  elv<-as.numeric(sitelist[j,"elv"])
  if(elv<=1000){
    sitelist[j,"elv_label"]<-"low"
  }else if(elv>1000){
    sitelist[j,"elv_label"]<-"high"
  }
}
#plot 
p1<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=Tmin),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour ="MTCO (¡ãC)")+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 ¡ãW","5 ¡ãW","0 ¡ãE","5 ¡ãE"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 ¡ãN","38 ¡ãN","40 ¡ãN","42 ¡ãN","44 ¡ãN"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  annotate("text",x=-10,y=44,label="(a)")
p2<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=Tmax),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red"),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour ="MTWA (¡ãC)")+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 ¡ãW","5 ¡ãW","0 ¡ãE","5 ¡ãE"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 ¡ãN","38 ¡ãN","40 ¡ãN","42 ¡ãN","44 ¡ãN"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  annotate("text",x=-10,y=44,label="(b)")
p3<-ggplot()+geom_point(data=Iberian_background,aes(lon,lat,color=alpha),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_point(data=sitelist,aes(lon,lat,shape=elv_label),size=2)+
  scale_colour_gradientn(colours = rev(c("blue","dodgerblue2","lightskyblue","white","gold","orangered","red")),na.value = "grey")+
  labs(x="Longitude",y="Latitude",colour =expression(alpha))+
  scale_x_continuous(breaks = c(-10,-5,0,5),labels=c("10 ¡ãW","5 ¡ãW","0 ¡ãE","5 ¡ãE"))+
  scale_y_continuous(breaks = c(36,38,40,42,44),labels=c("36 ¡ãN","38 ¡ãN","40 ¡ãN","42 ¡ãN","44 ¡ãN"))+ 
  scale_shape_manual(values=c("high"=17,"low"=15),guide = 'none')+
  annotate("text",x=-10,y=44,label="(c)")
p<-ggarrange(p1,p2,p3,ncol=1)
ggsave(file="Map of modern climate at Iberian Peninsula.jpeg",p,width=6,height=10)

########################################################################################################
######################################   Get the gradient  #############################################
########################################################################################################
plotdata<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Iberian/plotdata.csv",sep=""),row.names = 1)

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
modern_gradient<-get_gradient(plotdata,500,500)
#combine
str(modern_gradient)
colnames(modern_gradient)<-colnames(gradient_env)
for(j in c(1,3:ncol(modern_gradient))){
  modern_gradient[,j]<-as.numeric(modern_gradient[,j])
}

str(gradient_env)
for(j in c(1,3:ncol(gradient_env))){
  gradient_env[,j]<-as.numeric(gradient_env[,j])
}
gradient_env$agef<-factor(gradient_env$age,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                    "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))


#get anomaly to 0.5 ka
gradient_env$Tmin_anomaly<-gradient_env$mean_Tmin-rep(modern_gradient[,"mean_Tmin"],nlevels(gradient_env$agef))
gradient_env$gdd_anomaly<-gradient_env$mean_gdd-rep(modern_gradient[,"mean_gdd"],nlevels(gradient_env$agef))
gradient_env$alpha_anomaly<-gradient_env$mean_alpha-rep(modern_gradient[,"mean_alpha"],nlevels(gradient_env$agef))
gradient_env$Tmax_anomaly<-gradient_env$mean_Tmax-rep(modern_gradient[,"mean_Tmax"],nlevels(gradient_env$agef))

gradient_env$insol_summer_anomaly<-gradient_env$mean_insol_summer-rep(modern_gradient[,"mean_insol_summer"],nlevels(gradient_env$agef))
gradient_env$insol_winter_anomaly<-gradient_env$mean_insol_winter-rep(modern_gradient[,"mean_insol_winter"],nlevels(gradient_env$agef))

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
#"Hovemoller" type plots showing the anomalies in MTCO, MTWA and alpha relative to 0.5 ka; Figure 2
#anomaly to lon
summary(gradient_env$Tmin_anomaly)
summary(gradient_env$Tmax_anomaly)
summary(gradient_env$alpha_anomaly)

#High elv sites
#Make sure 0="white"
#min=-max
gradient_env_high<-gradient_env[which(gradient_env$elv_label=="high"),]
p1<-ggplot(data=gradient_env_high[which(!is.na(gradient_env_high$Tmin_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=Tmin_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+  
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red"),
                         limits=c(-18,18),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="(a) MTCO (¡ãC)")+labs(y="Age (kyr BP)")
#min=-max
p2<-ggplot(data=gradient_env_high[which(!is.na(gradient_env_high$Tmax_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=Tmax_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red"),
                         limits=c(-9,9),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="(b) MTWA (¡ãC)")+labs(y="Age (kyr BP)")
#min=-max
p3<-ggplot(data=gradient_env_high[which(!is.na(gradient_env_high$alpha_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=alpha_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = rev(c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red")),
                         limits=c(-0.5,0.5),na.value = "transparent")+
  theme(axis.text.x = element_blank())+
  labs(x="W                                 High elevation sites in the order of longitude                                 E")+
  labs(color=expression("(c)"~alpha))+labs(y="Age (kyr BP)")
#Low elv sites
#Make sure 0="white"
#min=-max
gradient_env_low<-gradient_env[which(gradient_env$elv_label=="low"),]
p4<-ggplot(data=gradient_env_low[which(!is.na(gradient_env_low$Tmin_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=Tmin_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+  
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red"),
                         limits=c(-18,18),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="(d) MTCO (¡ãC)")+labs(y="Age (kyr BP)")
#min=-max
p5<-ggplot(data=gradient_env_low[which(!is.na(gradient_env_low$Tmax_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=Tmax_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red"),
                         limits=c(-9,9),na.value = "transparent")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  labs(color="(e) MTWA (¡ãC)")+labs(y="Age (kyr BP)")
#min=-max
p6<-ggplot(data=gradient_env_low[which(!is.na(gradient_env_low$alpha_anomaly)),],aes(x=fct_reorder(site, lon),y=age,color=alpha_anomaly))+
  geom_point(size=2.5,shape=15)+theme_dark()+
  scale_y_continuous(trans = "reverse",breaks=1000*seq(1.5,11.5,by=2),labels=seq(1.5,11.5,by=2))+
  scale_colour_gradientn(colours = rev(c("darkblue","dodgerblue2","lightskyblue","white","orange","orangered","red")),
                         limits=c(-0.5,0.5),na.value = "transparent")+
  theme(axis.text.x = element_blank())+
  labs(x="W                                 Low elevation sites in the order of longitude                                 E")+
  labs(color=expression("(f)"~alpha))+labs(y="Age (kyr BP)")

p<-ggarrange(p1,p2,p3,p4,p5,p6,ncol=1)

ggsave(file="anomaly to lon.jpeg",p,width=8,height=9)

########################################################################################################
#####################################   Time series  Figure  ##########################################
########################################################################################################
if(!require(sf)){ install.packages("sf");library(sf)}

range<-function(data,nboot,variable){
  data$site<-as.factor(data$site)
  k_samples <- replicate(nboot, sample(1:nlevels(data$site),
                                       size = nlevels(data$site),
                                       replace = T))
  xboot<-data.frame()
  
  for(i in 1:nboot){
    k <- k_samples[, i]
    select<-levels(data$site)[k]
    datak<-data[levels(data$site) %in% select,]
    datak$age<-as.factor(datak$age)
    out<-aggregate(datak[,colnames(datak)==variable]~datak$age, FUN=mean)[,2]
    xboot<-rbind.data.frame(xboot,out)
  }
  return(xboot)
}
age_ka<-seq(0.5,11.5)
mean_site_env<-as.data.frame(age_ka)
mean_site_env$mean_alpha_anomaly<-aggregate(gradient_env$alpha_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$mean_Tmin_anomaly<-aggregate(gradient_env$Tmin_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$mean_Tmax_anomaly<-aggregate(gradient_env$Tmax_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$mean_insol_summer_anomaly<-aggregate(gradient_env$insol_summer_anomaly~gradient_env$age, FUN=mean)[,2]
mean_site_env$mean_insol_winter_anomaly<-aggregate(gradient_env$insol_winter_anomaly~gradient_env$age, FUN=mean)[,2]

mean_site_env$range_alpha_anomaly<-apply(range(gradient_env,nboot=1000,variable="alpha_anomaly"),2,sd,na.rm=T)
mean_site_env$range_Tmin_anomaly<-apply(range(gradient_env,nboot=1000,variable="Tmin_anomaly"),2,sd,na.rm=T)
mean_site_env$range_Tmax_anomaly<-apply(range(gradient_env,nboot=1000,variable="Tmax_anomaly"),2,sd,na.rm=T)
mean_site_env$range_insol_summer_anomaly<-apply(range(gradient_env,nboot=1000,variable="insol_summer_anomaly"),2,sd,na.rm=T)
mean_site_env$range_insol_winter_anomaly<-apply(range(gradient_env,nboot=1000,variable="insol_winter_anomaly"),2,sd,na.rm=T)

setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))

#Time series
p1<-ggplot()+theme_bw()+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly,
                                         ymin=mean_Tmin_anomaly-range_Tmin_anomaly,
                                         ymax=mean_Tmin_anomaly+range_Tmin_anomaly))+
  labs(y="MTCO anomaly (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5))+
  theme(axis.title.x = element_blank())+
  annotate("text", y= max(mean_site_env$mean_Tmin_anomaly+mean_site_env$range_Tmin_anomaly,na.rm=T), x =0,label="(a)",size=5)

p2<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly,
                                         ymin=mean_Tmax_anomaly-range_Tmax_anomaly,
                                         ymax=mean_Tmax_anomaly+range_Tmax_anomaly))+
  labs(y="MTWA anomaly (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5))+
  theme(axis.title.x = element_blank())+
  annotate("text", y= max(mean_site_env$mean_Tmax_anomaly+mean_site_env$range_Tmax_anomaly,na.rm=T), x =0,label="(b)",size=5)

p3<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_alpha_anomaly,
                                         ymin=mean_alpha_anomaly-range_alpha_anomaly,
                                         ymax=mean_alpha_anomaly+range_alpha_anomaly))+
  labs(y=expression(alpha* " anomaly"),x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5))+
  annotate("text", y= max(mean_site_env$mean_alpha_anomaly+mean_site_env$range_alpha_anomaly,na.rm=T), x =0,label="(c)",size=5)

p4<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_insol_summer_anomaly,
                                         ymin=mean_insol_summer_anomaly-range_insol_summer_anomaly,
                                         ymax=mean_insol_summer_anomaly+range_insol_summer_anomaly))+
  labs(y=bquote('Summer insolation anomaly (W'~m^-2~')'),x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5))+
  annotate("text", y= max(gradient_env$insol_summer_anomaly,na.rm=T), x =-0.5,label="(e)",size=5)

p5<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_insol_winter_anomaly,
                                         ymin=mean_insol_winter_anomaly-range_insol_winter_anomaly,
                                         ymax=mean_insol_winter_anomaly+range_insol_winter_anomaly))+
  labs(y=bquote('Winter insolation anomaly (W'~m^-2~')'))+
  scale_x_continuous(breaks= seq(0.5,11.5))+  
  theme(axis.title.x = element_blank())+
  annotate("text", y= max(gradient_env$insol_winter_anomaly,na.rm=T), x =-0.5,label="(d)",size=5)

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
df_model_Tmin_lon[which(df_model_Tmin_lon$p.slope<0.05),"age"]

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
df_model_Tmax_lon[which(df_model_Tmax_lon$p.slope<0.05),"age"]

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
df_model_alpha_lon[which(df_model_alpha_lon$p.slope<0.05),"age"]

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
df_model_Tmin_elv[which(df_model_Tmin_elv$p.slope<0.05),"age"]

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
df_model_Tmax_elv[which(df_model_Tmax_elv$p.slope<0.05),"age"]

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
df_model_alpha_elv[which(df_model_alpha_elv$p.slope<0.05),"age"]

#combine
df_model_lon<-rbind.data.frame(df_model_Tmin_lon,df_model_Tmax_lon,df_model_alpha_lon)
df_model_elv<-rbind.data.frame(df_model_Tmin_elv,df_model_Tmax_elv,df_model_alpha_elv)
df_model<-cbind.data.frame(df_model_lon,df_model_elv[,-1])
df_model$age<-df_model$age/1000
write.csv(df_model,paste(wd,"/Master Project/Data/Output data/Core/Iberian/df_model.csv",sep=""))

############ significance labels
plot_gradient<-gradient_env[which(gradient_env$age!=0),]
plot_gradient$agef<-factor(plot_gradient$agef,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                        "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="3.5 ka"|age=="4.5 ka"|age=="5.5 ka"|age=="6.5 ka"|age=="7.5 ka"|age=="9.5 ka"){
    plot_gradient[i,"label_Tmax_lon"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_Tmax_lon"]<-as.character(plot_gradient[i,"agef"])
  }
}
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="7.5 ka"|age=="8.5 ka"){
    plot_gradient[i,"label_Tmax_elv"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_Tmax_elv"]<-as.character(plot_gradient[i,"agef"])
  }
}
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="3.5 ka"|age=="4.5 ka"|age=="5.5 ka"|age=="6.5 ka"|age=="7.5 ka"|age=="8.5 ka"|age=="9.5 ka"){
    plot_gradient[i,"label_alpha_lon"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_alpha_lon"]<-as.character(plot_gradient[i,"agef"])
  }
}
for(i in 1:nrow(plot_gradient)){
  age<-plot_gradient[i,"agef"]
  if(age=="4.5 ka"|age=="5.5 ka"|age=="6.5 ka"|age=="7.5 ka"|age=="8.5 ka"|age=="9.5 ka"){
    plot_gradient[i,"label_alpha_elv"]<-as.character(paste(plot_gradient[i,"agef"]," significant"))
  }else{
    plot_gradient[i,"label_alpha_elv"]<-as.character(plot_gradient[i,"agef"])
  }
}
############ Figure alpha to lon
dat_text <- data.frame(
  label = plot_gradient$label_alpha_lon,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,alpha_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y=expression(alpha~"anomaly to 0.5 ka"))+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y = 0.33, label = label))
ggsave(file="alpha to lon.jpeg",p,width=8,height=8)

############ Figure alpha to elv
dat_text <- data.frame(
  label = plot_gradient$label_alpha_elv,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,alpha_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Elevation (km)",y=expression(alpha~"anomaly to 0.5 ka"))+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = 0.33, label = label))
ggsave(file="alpha to elv.jpeg",p,width=8,height=8)

############ Figure Tmax to lon
dat_text <- data.frame(
  label = plot_gradient$label_Tmax_lon,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,Tmax_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y="MTWA (¡ãC) anomaly to 0.5 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y = 7, label = label))
ggsave(file="Tmax to lon.jpeg",p,width=8,height=8)

############ Figure Tmax to elv
dat_text <- data.frame(
  label = plot_gradient$label_Tmax_elv,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,Tmax_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Elevation (km)",y="MTWA (¡ãC) anomaly to 0.5 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = 7, label = label))
ggsave(file="Tmax to elv.jpeg",p,width=8,height=8)


############ Figure Tmin to lon
dat_text <- data.frame(
  label = plot_gradient$agef,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(lon,Tmin_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y="MTCO (¡ãC) anomaly to 0.5 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = -2.5, y =-20, label = label))
ggsave(file="Tmin to lon.jpeg",p,width=8,height=8)

############ Figure Tmin to elv
dat_text <- data.frame(
  label = plot_gradient$agef,
  agef   = plot_gradient$agef
)
p<-ggplot(data=plot_gradient,aes(elv_km,Tmin_anomaly))+geom_point()+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Elevation (km)",y="MTCO (¡ãC) anomaly to 0.5 ka")+
  facet_wrap(vars(agef),nrow=3)+geom_abline(intercept=0,slope=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())+ 
  geom_text(data = dat_text, size=4,mapping = aes(x = 1.5, y = -20, label = label))
ggsave(file="Tmin to elv.jpeg",p,width=8,height=8)

#################################### west-east gradient change with elv
gradient_env$elv_km<-gradient_env$elv/1000
plot_gradient<-gradient_env[which(gradient_env$age!=0),]
plot_gradient$agef<-factor(plot_gradient$agef,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                       "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))

#Tmin
#low elv
df_model_Tmin_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_low[which(gradient_env_low$age==age),]
  m<-summary(lm(data=sub,Tmin_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmin_lon<-rbind.data.frame(df_model_Tmin_lon,dat)
}
colnames(df_model_Tmin_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmin_lon[which(df_model_Tmin_lon$p.slope<0.05),"age"]
#high elv
df_model_Tmin_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_high[which(gradient_env_high$age==age),]
  m<-summary(lm(data=sub,Tmin_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmin_lon<-rbind.data.frame(df_model_Tmin_lon,dat)
}
colnames(df_model_Tmin_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmin_lon[which(df_model_Tmin_lon$p.slope<0.05),"age"]

#Tmax
#low elv
df_model_Tmax_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_low[which(gradient_env_low$age==age),]
  m<-summary(lm(data=sub,Tmax_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmax_lon<-rbind.data.frame(df_model_Tmax_lon,dat)
}
colnames(df_model_Tmax_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmax_lon[which(df_model_Tmax_lon$p.slope<0.05),"age"]
#high elv
df_model_Tmax_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_high[which(gradient_env_high$age==age),]
  m<-summary(lm(data=sub,Tmax_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=2),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_Tmax_lon<-rbind.data.frame(df_model_Tmax_lon,dat)
}
colnames(df_model_Tmax_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_Tmax_lon[which(df_model_Tmax_lon$p.slope<0.05),"age"]

#alpha
#low elv
df_model_alpha_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_low[which(gradient_env_low$age==age),]
  m<-summary(lm(data=sub,alpha_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=3),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_alpha_lon<-rbind.data.frame(df_model_alpha_lon,dat)
}
colnames(df_model_alpha_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_alpha_lon[which(df_model_alpha_lon$p.slope<0.05),"age"]
#high elv
df_model_alpha_lon<-data.frame()
for(j in 1:12){
  age<-500+1000*(j-1)
  sub<-gradient_env_high[which(gradient_env_high$age==age),]
  m<-summary(lm(data=sub,alpha_anomaly~lon))
  dat<-cbind.data.frame(age,round(m[["coefficients"]]["lon","Estimate"],digits=3),
                        round(m[["coefficients"]]["lon","Pr(>|t|)"],digits=3),
                        round(-m[["coefficients"]]["(Intercept)","Estimate"]/m[["coefficients"]]["lon","Estimate"],digits=2))
  df_model_alpha_lon<-rbind.data.frame(df_model_alpha_lon,dat)
}
colnames(df_model_alpha_lon)<-c("age","slope_lon","p.slope_lon","x0_lon")
df_model_alpha_lon[which(df_model_alpha_lon$p.slope<0.05),"age"]

plot_gradient$elv_label<-factor(plot_gradient$elv_label,levels=c("high","low"))

p1<-ggplot(data=plot_gradient[which(plot_gradient$age!=500),],aes(lon,alpha_anomaly))+geom_point(size=1)+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y=expression(alpha~"anomaly to 0.5 ka"))+
  facet_grid(elv_label~agef)+geom_abline(intercept=0,slope=0)
ggsave(file="alpha west-east gradient to elv.jpeg",p1,width=8,height=6)

p2<-ggplot(data=plot_gradient[which(plot_gradient$age!=500),],aes(lon,Tmax_anomaly))+geom_point(size=1)+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y="MTWA (¡ãC) anomaly to 0.5 ka")+
  facet_grid(elv_label~agef)+geom_abline(intercept=0,slope=0)
ggsave(file="Tmax west-east gradient to elv.jpeg",p2,width=8,height=6)

p3<-ggplot(data=plot_gradient[which(plot_gradient$age!=500),],aes(lon,Tmin_anomaly))+geom_point(size=1)+theme_bw()+
  geom_smooth(method="lm",formula="y~x",col="red")+labs(x="Longitude (¡ãE)",y="MTCO (¡ãC) anomaly to 0.5 ka")+
  facet_grid(elv_label~agef)+geom_abline(intercept=0,slope=0)
ggsave(file="Tmin west-east gradient to elv.jpeg",p3,width=8,height=6)



#######################################################################################################
#############################   Reconstruction results to compare   ###################################
#######################################################################################################

########################### Mauri paper ###############################
if(!require(ncdf4)){ install.packages("ncdf4");library(ncdf4)}
library(raster)
library(rgdal) # for spTransform
range_grid<-function(data,nboot){
  k_samples <- replicate(nboot, sample(1:nrow(data),
                                       size = nrow(data),
                                       replace = T))
  xboot<-data.frame()
  
  for(i in 1:nboot){
    k <- k_samples[, i]
    datak<-data[k,]
    out<-mean(datak[,1])
    xboot<-rbind.data.frame(xboot,out)
  }
  return(xboot)
}
#winter temperature anomaly
d_winter<-data.frame()
for(i in c(100,c(1:12)*1000)){
  mauri_age<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Mauri paper/epoch2_tanom_djf_",i,".tif",sep=""))
  mauri_age <- rasterToPoints(mauri_age, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
  mauri_age <- data.frame(mauri_age@data, lon=coordinates(mauri_age)[,1],lat=coordinates(mauri_age)[,2])    # Assign coordinates to @data slot                     
  mauri_age_Iberian<-mauri_age[which(mauri_age$lon>=-10&mauri_age$lon<=5&mauri_age$lat>=36&mauri_age$lat<=44),]
  mauri_age_Iberian_mean<-colMeans(mauri_age_Iberian)[1]
  mauri_age_Iberian_sd<-apply(range_grid(mauri_age_Iberian,nboot=1000),2,sd,na.rm=T)
  d_age<-cbind.data.frame(i,mauri_age_Iberian_mean,mauri_age_Iberian_sd)
  colnames(d_age)<-c("age","mean","range")
  d_winter<-rbind.data.frame(d_winter,d_age)
}
d_winter$age_kyr<-d_winter$age/1000

#summer temperature anomaly
d_summer<-data.frame()
for(i in c(100,c(1:12)*1000)){
  mauri_age<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Mauri paper/epoch2_tanom_jja_",i,".tif",sep=""))
  mauri_age <- rasterToPoints(mauri_age, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
  mauri_age <- data.frame(mauri_age@data, lon=coordinates(mauri_age)[,1],lat=coordinates(mauri_age)[,2])    # Assign coordinates to @data slot                     
  mauri_age_Iberian<-mauri_age[which(mauri_age$lon>=-10&mauri_age$lon<=5&mauri_age$lat>=36&mauri_age$lat<=44),]
  mauri_age_Iberian_mean<-colMeans(mauri_age_Iberian)[1]
  mauri_age_Iberian_sd<-apply(range_grid(mauri_age_Iberian,nboot=1000),2,sd,na.rm=T)
  d_age<-cbind.data.frame(i,mauri_age_Iberian_mean,mauri_age_Iberian_sd)
  colnames(d_age)<-c("age","mean","range")
  d_summer<-rbind.data.frame(d_summer,d_age)
}
d_summer$age_kyr<-d_summer$age/1000

#alpha anomaly
d_alpha<-data.frame()
for(i in c(100,c(1:12)*1000)){
  mauri_age<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Mauri paper/epoch2_anom_alpha_",i,".tif",sep=""))
  mauri_age <- rasterToPoints(mauri_age, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
  mauri_age <- data.frame(mauri_age@data, lon=coordinates(mauri_age)[,1],lat=coordinates(mauri_age)[,2])    # Assign coordinates to @data slot                     
  mauri_age_Iberian<-mauri_age[which(mauri_age$lon>=-10&mauri_age$lon<=5&mauri_age$lat>=36&mauri_age$lat<=44),]
  #Mauri paper defines the ratio of actual to potential evapotranspiration as alpha, 
  #but our paper defines the ratio of actual to equilibrium evapotranspiration as alpha,
  #potential=1.26equilibrium, so alpha_this/alpha_mauri=potential/equilibrium=1.26,
  #so the alpha result of mauri paper should be multiplied by 1.26
  mauri_age_Iberian[,1]<-mauri_age_Iberian[,1]*1.26 
  mauri_age_Iberian_mean<-colMeans(mauri_age_Iberian)[1]
  mauri_age_Iberian_sd<-apply(range_grid(mauri_age_Iberian,nboot=1000),2,sd,na.rm=T)
  d_age<-cbind.data.frame(i,mauri_age_Iberian_mean,mauri_age_Iberian_sd)
  colnames(d_age)<-c("age","mean","range")
  d_alpha<-rbind.data.frame(d_alpha,d_age)
}
d_alpha$age_kyr<-d_alpha$age/1000

###################### Kaufman paper ############################
if(!require(devtools)){ install.packages("devtools");library(devtools)}
devtools::install_github("nickmckay/LiPD-Utilities", subdir = "R")
library(lipdR)

#Read data and find sites in Iberian Peninsula
D <- readLipd(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Kaufman paper/Temp12k_v1_0_0_LiPD",sep=""))

D_Iberian<-list();j=0
for(i in 1:length(D)){
  location<-D[[i]]$geo
  if(location$longitude>=-10&location$longitude<=5&location$latitude>=36&location$latitude<=44){
    j=j+1
    D_Iberian[[j]]<-D[[i]]
  }
}

#Plot where the sites are
location_list<-data.frame()
for(i in 1:length(D_Iberian)){
  location<-cbind.data.frame(D_Iberian[[i]]$geo$longitude,D_Iberian[[i]]$geo$latitude,D_Iberian[[i]]$geo$siteName)
  colnames(location)<-c("lon","lat","site")
  location_list<-rbind.data.frame(location_list,location)
}
ggplot(data=location_list,aes(lon,lat))+theme_bw()+geom_point()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') +
  geom_text(label=location_list$site,nudge_y = 0,nudge_x = 0)

#Get temperature seasonality type, age and values
D_Iberian_winter<-data.frame();D_Iberian_summer<-data.frame()
for(i in 1:length(D_Iberian)){
  L<-D_Iberian[[i]];print(i)
  for(j in 1:length((L[["paleoData"]][[1]][["measurementTable"]]))){
    Lj1<-L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature"]][["interpretation"]][[1]][["seasonalityGeneral"]]
    Lj2<-L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-1"]][["interpretation"]][[1]][["seasonalityGeneral"]]
    Lj3<-L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-2"]][["interpretation"]][[1]][["seasonalityGeneral"]]
    Lj4<-L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-3"]][["interpretation"]][[1]][["seasonalityGeneral"]]
    
    print(paste(L[["geo"]][["siteName"]],"temperature:",Lj1,"; temperature-1:",Lj2,"temperature-2:",Lj1,"; temperature-3:",Lj2))
    if(!is.null(Lj1)){
      if(Lj1=="winterOnly"|Lj1=="winter+"){
        winter_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature"]][["values"]])
        colnames(winter_sample)<-c("lon","lat","site","age","Tmin")
        D_Iberian_winter<-rbind.data.frame(D_Iberian_winter,winter_sample)
      }else if(Lj1=="summerOnly"|Lj1=="summer+"){
        summer_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature"]][["values"]])
        colnames(summer_sample)<-c("lon","lat","site","age","Tmax")
        D_Iberian_summer<-rbind.data.frame(D_Iberian_summer,summer_sample)
      }
    }
    if(!is.null(Lj2)){
      if(Lj2=="winterOnly"|Lj2=="winter+"){
        winter_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-1"]][["values"]])
        colnames(winter_sample)<-c("lon","lat","site","age","Tmin")
        D_Iberian_winter<-rbind.data.frame(D_Iberian_winter,winter_sample)
      }else if(Lj2=="summerOnly"|Lj2=="summer+"){
        summer_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-1"]][["values"]])
        colnames(summer_sample)<-c("lon","lat","site","age","Tmax")
        D_Iberian_summer<-rbind.data.frame(D_Iberian_summer,summer_sample)
      }
    }
    if(!is.null(Lj3)){
      if(Lj3=="winterOnly"|Lj3=="winter+"){
        winter_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-2"]][["values"]])
        colnames(winter_sample)<-c("lon","lat","site","age","Tmin")
        D_Iberian_winter<-rbind.data.frame(D_Iberian_winter,winter_sample)
      }else if(Lj3=="summerOnly"|Lj3=="summer+"){
        summer_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-2"]][["values"]])
        colnames(summer_sample)<-c("lon","lat","site","age","Tmax")
        D_Iberian_summer<-rbind.data.frame(D_Iberian_summer,summer_sample)
      }
    }
    if(!is.null(Lj4)){
      if(Lj4=="winterOnly"|Lj4=="winter+"){
        winter_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-3"]][["values"]])
        colnames(winter_sample)<-c("lon","lat","site","age","Tmin")
        D_Iberian_winter<-rbind.data.frame(D_Iberian_winter,winter_sample)
      }else if(Lj4=="summerOnly"|Lj4=="summer+"){
        summer_sample<-cbind.data.frame(L[["geo"]][["longitude"]],L[["geo"]][["latitude"]],L[["geo"]][["siteName"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["age"]][["values"]],
                                        L[["paleoData"]][[1]][["measurementTable"]][[j]][["temperature-3"]][["values"]])
        colnames(summer_sample)<-c("lon","lat","site","age","Tmax")
        D_Iberian_summer<-rbind.data.frame(D_Iberian_summer,summer_sample)
      }
    }
  }
}


#Apply +-500 year bins
D_Iberian_winter$site<-as.factor(D_Iberian_winter$site)
D_Iberian_summer$site<-as.factor(D_Iberian_summer$site)
D_Iberian_winter$age_kyr<-D_Iberian_winter$age/1000
D_Iberian_summer$age_kyr<-D_Iberian_summer$age/1000

ggplot(D_Iberian_winter,aes(age_kyr,Tmin,color=site))+geom_point()+geom_line()+theme(legend.position = "bottom")+
  labs(x="age (kyr)",y="MTCO (¡ãC)")
ggplot(D_Iberian_summer,aes(age_kyr,Tmax,color=site))+geom_point()+geom_line()+theme(legend.position = "bottom")+
  labs(x="age (kyr)",y="MTWA (¡ãC)")

get_gradient2<-function(plotdata,variable,age,range){
  plotdata$site<-as.factor(plotdata$site)
  mean_env<-data.frame(matrix(NA,ncol=5,nrow=nlevels(plotdata$site)))
  for(m in 1:nlevels(plotdata$site)){
    plotdata_each<-plotdata[which(plotdata$site==levels(plotdata$site)[m]),]
    
    mean<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),variable]),na.rm = TRUE)
    mean_env[m,]<-c(age,levels(plotdata$site)[m],unique(plotdata_each$lon),unique(plotdata_each$lat),mean)
    
  }
  colnames(mean_env)<-c("age","site","lon","lat","mean")
  return(mean_env)
}
#get gradient
gradient_winter_Kaufman<-rbind.data.frame(
  get_gradient2(D_Iberian_winter, "Tmin",500,500),get_gradient2(D_Iberian_winter, "Tmin",1500,500),
  get_gradient2(D_Iberian_winter, "Tmin",2500,500),get_gradient2(D_Iberian_winter, "Tmin",3500,500),
  get_gradient2(D_Iberian_winter, "Tmin",4500,500),get_gradient2(D_Iberian_winter, "Tmin",5500,500),
  get_gradient2(D_Iberian_winter, "Tmin",6500,500),get_gradient2(D_Iberian_winter, "Tmin",7500,500),
  get_gradient2(D_Iberian_winter, "Tmin",8500,500),get_gradient2(D_Iberian_winter, "Tmin",9500,500),
  get_gradient2(D_Iberian_winter, "Tmin",10500,500),get_gradient2(D_Iberian_winter, "Tmin",11500,500))
gradient_summer_Kaufman<-rbind.data.frame(
  get_gradient2(D_Iberian_summer, "Tmax",500,500),get_gradient2(D_Iberian_summer, "Tmax",1500,500),
  get_gradient2(D_Iberian_summer, "Tmax",2500,500),get_gradient2(D_Iberian_summer, "Tmax",3500,500),
  get_gradient2(D_Iberian_summer, "Tmax",4500,500),get_gradient2(D_Iberian_summer, "Tmax",5500,500),
  get_gradient2(D_Iberian_summer, "Tmax",6500,500),get_gradient2(D_Iberian_summer, "Tmax",7500,500),
  get_gradient2(D_Iberian_summer, "Tmax",8500,500),get_gradient2(D_Iberian_summer, "Tmax",9500,500),
  get_gradient2(D_Iberian_summer, "Tmax",10500,500),get_gradient2(D_Iberian_summer, "Tmax",11500,500))
#get modern gradient
modern_gradient_winter<-get_gradient2(D_Iberian_winter,"Tmin",500,500)
modern_gradient_summer<-get_gradient2(D_Iberian_summer,"Tmax",500,500)

#get anomalies
#winter
str(modern_gradient_winter)
for(j in c(1,3:ncol(modern_gradient_winter))){
  modern_gradient_winter[,j]<-as.numeric(modern_gradient_winter[,j])
}
str(gradient_winter_Kaufman)
for(j in c(1,3:ncol(gradient_winter_Kaufman))){
  gradient_winter_Kaufman[,j]<-as.numeric(gradient_winter_Kaufman[,j])
}
gradient_winter_Kaufman$agef<-factor(gradient_winter_Kaufman$age,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                                          "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))
gradient_winter_Kaufman$Tmin_anomaly<-gradient_winter_Kaufman$mean-rep(modern_gradient_winter[,"mean"],nlevels(gradient_winter_Kaufman$agef))
gradient_winter_Kaufman<-na.omit(gradient_winter_Kaufman)

#summer
str(modern_gradient_summer)
for(j in c(1,3:ncol(modern_gradient_summer))){
  modern_gradient_summer[,j]<-as.numeric(modern_gradient_summer[,j])
}
str(gradient_summer_Kaufman)
for(j in c(1,3:ncol(gradient_summer_Kaufman))){
  gradient_summer_Kaufman[,j]<-as.numeric(gradient_summer_Kaufman[,j])
}
gradient_summer_Kaufman$agef<-factor(gradient_summer_Kaufman$age,labels=c("0.5 ka","1.5 ka","2.5 ka","3.5 ka","4.5 ka","5.5 ka",
                                                                          "6.5 ka","7.5 ka","8.5 ka","9.5 ka","10.5 ka","11.5 ka"))

gradient_summer_Kaufman$Tmax_anomaly<-gradient_summer_Kaufman$mean-rep(modern_gradient_summer[,"mean"],nlevels(gradient_summer_Kaufman$agef))
gradient_summer_Kaufman<-na.omit(gradient_summer_Kaufman)

unique(gradient_summer_Kaufman$site)
unique(gradient_winter_Kaufman$site)

#remove marine sites
gradient_summer_Kaufman<-gradient_summer_Kaufman[-which(gradient_summer_Kaufman$site=="PP10-07"),]
gradient_winter_Kaufman<-gradient_winter_Kaufman[-which(gradient_winter_Kaufman$site=="PP10-07"),]

ggplot()+geom_point(data=gradient_summer_Kaufman,aes(lon,lat),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') 
ggplot()+geom_point(data=gradient_winter_Kaufman,aes(lon,lat),size=3,shape=15)+theme_bw()+
  geom_polygon(data = region,aes(x=long, y = lat, group = group),alpha=0,color='black') 

ggplot(gradient_winter_Kaufman,aes(age,Tmin_anomaly,color=site))+theme_bw()+geom_point()+geom_line()+theme(legend.position = "bottom")+
  labs(x="age (kyr)",y="MTCO (¡ãC)")
ggplot(gradient_summer_Kaufman,aes(age,Tmax_anomaly,color=site))+theme_bw()+geom_point()+geom_line()+theme(legend.position = "bottom")+
  labs(x="age (kyr)",y="MTWA (¡ãC)")

p1<-ggplot(D_Iberian_winter[which(D_Iberian_winter$site%in%unique(gradient_winter_Kaufman$site)),],aes(age_kyr,Tmin,color=site))+
  theme_bw()+geom_point()+geom_line()+
  labs(x=NULL,y="MTCO (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5),limits=c(0,12))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+
  annotate("text", y= 11, x =0,label="(a)",size=5)
p2<-ggplot(D_Iberian_summer[which(D_Iberian_summer$site%in%unique(gradient_summer_Kaufman$site)),],aes(age_kyr,Tmax,color=site))+
  theme_bw()+geom_point()+geom_line()+theme(legend.position = "bottom")+
  labs(x="Age (kyr)",y="MTWA (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5),limits=c(0,12))+
  annotate("text", y= 30, x =0,label="(b)",size=5)
p<-ggarrange(p1,p2,ncol=1)
ggsave(file="Kaufman individual sites.jpeg",p,width=8,height=8)

#Get mean values across sites
if(!require(sf)){ install.packages("sf");library(sf)}
rm(mean)
range<-function(data,nboot,variable){
  data$site<-as.factor(data$site)
  k_samples <- replicate(nboot, sample(1:nlevels(data$site),
                                       size = nlevels(data$site),
                                       replace = T))
  xboot<-data.frame()
  
  for(i in 1:nboot){
    k <- k_samples[, i]
    select<-levels(data$site)[k]
    datak<-data[levels(data$site) %in% select,]
    datak$age<-as.factor(datak$age)
    out<-aggregate(datak[,colnames(datak)==variable]~datak$age, FUN=mean)[,2]
    xboot<-rbind.data.frame(xboot,out)
  }
  return(xboot)
}
age_ka<-seq(0.5,10.5)
mean_site_winter_Kaufman<-as.data.frame(age_ka)
mean_site_winter_Kaufman$mean_Tmin_anomaly<-aggregate(gradient_winter_Kaufman$Tmin_anomaly~gradient_winter_Kaufman$age, FUN=mean)[,2]
mean_site_winter_Kaufman$range_Tmin_anomaly<-apply(range(gradient_winter_Kaufman,nboot=1000,variable="Tmin_anomaly"),2,sd,na.rm=T)

mean_site_summer_Kaufman<-as.data.frame(age_ka)
mean_site_summer_Kaufman$mean_Tmax_anomaly<-aggregate(gradient_summer_Kaufman$Tmax_anomaly~gradient_summer_Kaufman$age, FUN=mean)[,2]
mean_site_summer_Kaufman$range_Tmax_anomaly<-apply(range(gradient_summer_Kaufman,nboot=1000,variable="Tmax_anomaly"),2,sd,na.rm=T)

########################### Torroso paper ###############################
Torroso_Tmin_anomaly<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Torroso paper/climate_reconstruction.nc",sep=""),varname="tmin_jan_anomaly")
Torroso_Tmax_anomaly<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Torroso paper/climate_reconstruction.nc",sep=""),varname="tmax_jul_anomaly")
Torroso_pr_anomaly<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/Torroso paper/climate_reconstruction.nc",sep=""),varname="annprc_anomaly")

# Convert raster to dataframe
Torroso_Tmin_anomaly<-rasterToPoints(Torroso_Tmin_anomaly, spatial=TRUE)
Torroso_Tmin_anomaly<-spTransform(Torroso_Tmin_anomaly, CRS(Torroso_Tmin_anomaly@proj4string@projargs))
Torroso_Tmin_anomaly <- data.frame(Torroso_Tmin_anomaly@data)           

Torroso_Tmax_anomaly<-rasterToPoints(Torroso_Tmax_anomaly, spatial=TRUE)
Torroso_Tmax_anomaly<-spTransform(Torroso_Tmax_anomaly, CRS(Torroso_Tmax_anomaly@proj4string@projargs))
Torroso_Tmax_anomaly <- data.frame(Torroso_Tmax_anomaly@data)   

Torroso_pr_anomaly<-rasterToPoints(Torroso_pr_anomaly, spatial=TRUE)
Torroso_pr_anomaly<-spTransform(Torroso_pr_anomaly, CRS(Torroso_pr_anomaly@proj4string@projargs))
Torroso_pr_anomaly <- data.frame(Torroso_pr_anomaly@data)               

# get the average and sd
range_grid_allage<-function(data,nboot){
  k_samples <- replicate(nboot, sample(1:nrow(data),
                                       size = nrow(data),
                                       replace = T))
  xboot<-data.frame()
  
  for(i in 1:nboot){
    k <- k_samples[, i]
    datak<-data[k,]
    out<-colMeans(datak)
    xboot<-rbind.data.frame(xboot,out)
  }
  return(xboot)
}
age_ka<-c(3:15)
Torroso_Tmin_anomaly_mean<-colMeans(Torroso_Tmin_anomaly)
Torroso_Tmin_anomaly_sd<-apply(range_grid_allage(Torroso_Tmin_anomaly,nboot=1000),2,sd,na.rm=T)
Torroso_Tmax_anomaly_mean<-colMeans(Torroso_Tmax_anomaly)
Torroso_Tmax_anomaly_sd<-apply(range_grid_allage(Torroso_Tmax_anomaly,nboot=1000),2,sd,na.rm=T)
Torroso_pr_anomaly_mean<-colMeans(Torroso_pr_anomaly)
Torroso_pr_anomaly_sd<-apply(range_grid_allage(Torroso_pr_anomaly,nboot=1000),2,sd,na.rm=T)
Torroso_plot<-cbind.data.frame(age_ka,Torroso_Tmin_anomaly_mean,Torroso_Tmin_anomaly_sd,
                               Torroso_Tmax_anomaly_mean,Torroso_Tmax_anomaly_sd,
                               Torroso_pr_anomaly_mean,Torroso_pr_anomaly_sd)
colnames(Torroso_plot)<-c("age_ka","mean_Tmin_anomaly","range_Tmin_anomaly",
                          "mean_Tmax_anomaly","range_Tmax_anomaly",
                          "mean_pr_anomaly","range_pr_anomaly")
Torroso_plot<-Torroso_plot[which(Torroso_plot$age_ka<=12),]

############################### Compare composite curves 
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))

#Time series
min_Tmin<- -6;max_Tmin<-1.5
min_Tmax<--4;max_Tmax<-2
min_alpha<- -0.04;max_alpha<-0.1
p_this1<-ggplot()+theme_bw()+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmin_anomaly,
                                         ymin=mean_Tmin_anomaly-range_Tmin_anomaly,
                                         ymax=mean_Tmin_anomaly+range_Tmin_anomaly))+
  labs(y="MTCO anomaly (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmin,max_Tmin,by=2),limits=c(min_Tmin,max_Tmin))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  annotate("text", y= max_Tmin, x =0,label="(a)",size=5)

p_this2<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_Tmax_anomaly,
                                         ymin=mean_Tmax_anomaly-range_Tmax_anomaly,
                                         ymax=mean_Tmax_anomaly+range_Tmax_anomaly))+
  labs(y="MTWA anomaly (¡ãC)")+scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmax,max_Tmax),limits=c(min_Tmax,max_Tmax))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
  annotate("text", y= max_Tmax, x =0,label="(b)",size=5)

p_this3<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=1.5)+
  geom_line(data=mean_site_env,aes(age_ka,mean_alpha_anomaly),size=2)+
  geom_pointrange(data=mean_site_env,aes(age_ka,mean_alpha_anomaly,
                                         ymin=mean_alpha_anomaly-range_alpha_anomaly,
                                         ymax=mean_alpha_anomaly+range_alpha_anomaly))+
  labs(y=expression(alpha* " anomaly"),x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_alpha,max_alpha,by=0.04),limits=c(min_alpha,max_alpha))+
  annotate("text", y= max_alpha, x =0,label="(c)",size=5)

p_Mauri1<-ggplot()+theme_bw()+
  geom_point(data=d_winter,aes(age_kyr,mean),size=1.5)+geom_line(data=d_winter,aes(age_kyr,mean),size=2)+
  geom_pointrange(data=d_winter,aes(age_kyr,mean,ymin=mean-range,ymax=mean+range))+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmin,max_Tmin,by=2),limits=c(min_Tmin,max_Tmin))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmin, x =0,label="(d)",size=5)

p_Mauri2<-ggplot()+theme_bw()+
  geom_point(data=d_summer,aes(age_kyr,mean),size=1.5)+geom_line(data=d_summer,aes(age_kyr,mean),size=2)+
  geom_pointrange(data=d_summer,aes(age_kyr,mean,ymin=mean-range,ymax=mean+range))+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmax,max_Tmax),limits=c(min_Tmax,max_Tmax))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmax, x =0,label="(e)",size=5)

p_Mauri3<-ggplot()+theme_bw()+
  geom_point(data=d_alpha,aes(age_kyr,mean),size=1.5)+geom_line(data=d_alpha,aes(age_kyr,mean),size=2)+
  geom_pointrange(data=d_alpha,aes(age_kyr,mean,ymin=mean-range,ymax=mean+range))+
  labs(y=NULL,x="Age (kyr BP)")+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_alpha,max_alpha,by=0.04),limits=c(min_alpha,max_alpha))+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_alpha, x =0,label="(f)",size=5)

p_Kaufman1<-ggplot()+theme_bw()+
  geom_point(data=mean_site_winter_Kaufman,aes(age_ka,mean_Tmin_anomaly),size=1.5)+
  geom_line(data=mean_site_winter_Kaufman,aes(age_ka,mean_Tmin_anomaly),size=2)+
  geom_pointrange(data=mean_site_winter_Kaufman,aes(age_ka,mean_Tmin_anomaly,
                                                    ymin=mean_Tmin_anomaly-range_Tmin_anomaly,
                                                    ymax=mean_Tmin_anomaly+range_Tmin_anomaly))+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmin,max_Tmin,by=2),limits=c(min_Tmin,max_Tmin))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmin, x =0,label="(g)",size=5)

p_Kaufman2<-ggplot()+theme_bw()+
  theme(legend.position="none")+
  geom_point(data=mean_site_summer_Kaufman,aes(age_ka,mean_Tmax_anomaly),size=1.5)+
  geom_line(data=mean_site_summer_Kaufman,aes(age_ka,mean_Tmax_anomaly),size=2)+
  geom_pointrange(data=mean_site_summer_Kaufman,aes(age_ka,mean_Tmax_anomaly,
                                                    ymin=mean_Tmax_anomaly-range_Tmax_anomaly,
                                                    ymax=mean_Tmax_anomaly+range_Tmax_anomaly))+
  labs(x="Age (kyr BP)")+scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmax,max_Tmax),limits=c(min_Tmax,max_Tmax))+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmax, x =0,label="(h)",size=5)

p_Torroso1<-ggplot()+theme_bw()+
  geom_point(data=Torroso_plot,aes(age_ka,mean_Tmin_anomaly),size=1.5)+
  geom_line(data=Torroso_plot,aes(age_ka,mean_Tmin_anomaly),size=2)+
  geom_pointrange(data=Torroso_plot,aes(age_ka,mean_Tmin_anomaly,
                                        ymin=mean_Tmin_anomaly-range_Tmin_anomaly,
                                        ymax=mean_Tmin_anomaly+range_Tmin_anomaly))+
  scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmin,max_Tmin,by=2),limits=c(min_Tmin,max_Tmin))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmin, x =0,label="(i)",size=5)

p_Torroso2<-ggplot()+theme_bw()+
  geom_point(data=Torroso_plot,aes(age_ka,mean_Tmax_anomaly),size=1.5)+
  geom_line(data=Torroso_plot,aes(age_ka,mean_Tmax_anomaly),size=2)+
  geom_pointrange(data=Torroso_plot,aes(age_ka,mean_Tmax_anomaly,
                                        ymin=mean_Tmax_anomaly-range_Tmax_anomaly,
                                        ymax=mean_Tmax_anomaly+range_Tmax_anomaly))+
  labs(x="Age (kyr BP)")+scale_x_continuous(breaks= seq(0.5,11.5,by=2))+
  scale_y_continuous(breaks= seq(min_Tmax,max_Tmax),limits=c(min_Tmax,max_Tmax))+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  annotate("text", y= max_Tmax, x =0,label="(j)",size=5)

p<-ggarrange(p_this1,p_Mauri1,p_Kaufman1,p_Torroso1,p_this2,p_Mauri2,p_Kaufman2,p_Torroso2,p_this3,p_Mauri3,ncol=4)

ggsave(file="Compare composite curves.jpeg",p,width=11,height=10)


#######################################################################################################
##################################   Model results to compare   #######################################
#######################################################################################################
if(!require(ncdf4)){ install.packages("ncdf4");library(ncdf4)}
library(raster)
library(rgdal) # for spTransform
#PACMEDY
#pr = precipitation,tas = surface air temperature (2m)
#Import data
tas_MPI_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_MPI.nc",sep=""))
pr_MPI_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_MPI.nc",sep=""))

tas_AWI_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_AWI.nc",sep=""),varname="temp2")
pr_AWI_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_AWI.nc",sep=""))

tas_TR5AS_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR5AS.nc",sep=""))
pr_TR5AS_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR5AS.nc",sep=""))

tas_TR6AV_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR6AV.nc",sep=""))
pr_TR6AV_nc<-brick(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR6AV.nc",sep=""))

################################################# Extract Iberian Peninsula
#MPI
tas_MPI <- rasterToPoints(tas_MPI_nc)# Convert raster to SpatialPointsDataFrame
lat<-tas_MPI[,2];summary(lat)
lon<-tas_MPI[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
tas_MPI_Iberian<-tas_MPI[which(lon1>=-10&lon1<=5&lat>=36&lat<=44),]
write.csv(tas_MPI_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_MPI_Iberian.csv",sep=""))

pr_MPI <- rasterToPoints(pr_MPI_nc)# Convert raster to SpatialPointsDataFrame
lat<-pr_MPI[,2];summary(lat)
lon<-pr_MPI[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
pr_MPI_Iberian<-pr_MPI[which(lon1>=-10&lon1<=5&lat>=36&lat<=44),]
write.csv(pr_MPI_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_MPI_Iberian.csv",sep=""))

#AWI
tas_AWI <- rasterToPoints(tas_AWI_nc)# Convert raster to SpatialPointsDataFrame
lat<-tas_AWI[,2];summary(lat)
lon<-tas_AWI[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
tas_AWI_Iberian<-tas_AWI[which(lon1>=-10&lon1<=5&lat>=36&lat<=44),]
write.csv(tas_AWI_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_AWI_Iberian.csv",sep=""))

pr_AWI <- rasterToPoints(pr_AWI_nc)# Convert raster to SpatialPointsDataFrame
lat<-pr_AWI[,2];summary(lat)
lon<-pr_AWI[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
pr_AWI_Iberian<-pr_AWI[which(lon1>=-10&lon1<=5&lat>=36&lat<=44),]
write.csv(pr_AWI_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_AWI_Iberian.csv",sep=""))

#TR5AS
tas_TR5AS <- rasterToPoints(tas_TR5AS_nc)# Convert raster to SpatialPointsDataFrame
lat<-tas_TR5AS[,2];summary(lat)
lon<-tas_TR5AS[,1];summary(lon)
tas_TR5AS_Iberian<-tas_TR5AS[which(lon>=-10&lon<=5&lat>=36&lat<=44),]
write.csv(tas_TR5AS_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR5AS_Iberian.csv",sep=""))

pr_TR5AS <- rasterToPoints(pr_TR5AS_nc)# Convert raster to SpatialPointsDataFrame
lat<-pr_TR5AS[,2];summary(lat)
lon<-pr_TR5AS[,1];summary(lon)
pr_TR5AS_Iberian<-pr_TR5AS[which(lon>=-10&lon<=5&lat>=36&lat<=44),]
write.csv(pr_TR5AS_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR5AS_Iberian.csv",sep=""))

#TR6AV
tas_TR6AV <- rasterToPoints(tas_TR6AV_nc)# Convert raster to SpatialPointsDataFrame
lat<-tas_TR6AV[,2];summary(lat)
lon<-tas_TR6AV[,1];summary(lon)
tas_TR6AV_Iberian<-tas_TR6AV[which(lon>=-10&lon<=5&lat>=36&lat<=44),]
write.csv(tas_TR6AV_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR6AV_Iberian.csv",sep=""))

pr_TR6AV <- rasterToPoints(pr_TR6AV_nc)# Convert raster to SpatialPointsDataFrame
lat<-pr_TR6AV[,2];summary(lat)
lon<-pr_TR6AV[,1];summary(lon)
pr_TR6AV_Iberian<-pr_TR6AV[which(lon>=-10&lon<=5&lat>=36&lat<=44),]
write.csv(pr_TR6AV_Iberian,paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR6AV_Iberian.csv",sep=""))

############################################################ Get the model mean
#Get the mean across sites
if(!require(matrixStats)){ install.packages("matrixStats");library(matrixStats)}
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}

##################################### MPI ################################################
tas_MPI_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_MPI_Iberian.csv",sep=""),row.names=1)
pr_MPI_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_MPI_Iberian.csv",sep=""),row.names=1)

#Temperature
tas_Iberian<-tas_MPI_Iberian[,-c(1,2)]-273.15 #the unit is K, convert to degree celcius
age<-rev(rep(seq(101:7950),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),7850)
df1<-cbind.data.frame(age,age_100yrbin,month,colMeans(tas_Iberian),colSds(as.matrix(tas_Iberian)))
colnames(df1)<-c("age","age_100yrbin","month","T","sd_T")

df1_winter<-df1[which(df1$month==12|df1$month==1|df1$month==2),]
df1_winter<-aggregate(df1_winter,by=list(df1_winter$age_100yrbin),FUN=mean)
ggplot(df1_winter,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="MPI  MTCO (¡ãC)",x="Age (yr BP)")
df1_summer<-df1[which(df1$month==6|df1$month==7|df1$month==8),]
df1_summer<-aggregate(df1_summer,by=list(df1_summer$age_100yrbin),FUN=mean)
ggplot(df1_summer,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="MPI  MTWA (¡ãC)",x="Age (yr BP)")

#Precipitation
pr_Iberian<-pr_MPI_Iberian[,-c(1,2)]
age<-rev(rep(seq(101:7950),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),7850)
df1_pr<-cbind.data.frame(age,age_100yrbin,month,colMeans(pr_Iberian),colSds(as.matrix(pr_Iberian)))
colnames(df1_pr)<-c("age","age_100yrbin","month","prep","sd_prep")

df1_pr<-aggregate(df1_pr,by=list(df1_pr$age_100yrbin),FUN=mean)
ggplot(df1_pr,aes(age_100yrbin,prep))+geom_point()+geom_line()+theme_bw()+
  labs(y=bquote("MPI  mean monthly precipitation (kg"~s^-1~m^-2~")"),x="Age (yr BP)")

########################################### AWI ########################################
tas_AWI_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_AWI_Iberian.csv",sep=""),row.names=1)
pr_AWI_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_AWI_Iberian.csv",sep=""),row.names=1)

#Temperature
tas_Iberian<-tas_AWI_Iberian[,-c(1,2)]-273.15 #the unit is K, convert to degree celcius
age<-rev(rep(seq(0:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6001)
df2<-cbind.data.frame(age,age_100yrbin,month,colMeans(tas_Iberian),colSds(as.matrix(tas_Iberian)))
colnames(df2)<-c("age","age_100yrbin","month","T","sd_T")

df2_winter<-df2[which(df2$month==12|df2$month==1|df2$month==2),]
df2_winter<-aggregate(df2_winter,by=list(df2_winter$age_100yrbin),FUN=mean)
ggplot(df2_winter,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="AWI  MTCO (¡ãC)",x="Age (yr BP)")
df2_summer<-df2[which(df2$month==6|df2$month==7|df2$month==8),]
df2_summer<-aggregate(df2_summer,by=list(df2_summer$age_100yrbin),FUN=mean)
ggplot(df2_summer,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="AWI  MTWA (¡ãC)",x="Age (yr BP)")

#Precipitation
pr_Iberian<-pr_AWI_Iberian[,-c(1,2)]
age<-rev(rep(seq(0:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6001)
df2_pr<-cbind.data.frame(age,age_100yrbin,month,colMeans(pr_Iberian),colSds(as.matrix(pr_Iberian)))
colnames(df2_pr)<-c("age","age_100yrbin","month","prep","sd_prep")

df2_pr<-aggregate(df2_pr,by=list(df2_pr$age_100yrbin),FUN=mean)
ggplot(df2_pr,aes(age_100yrbin,prep))+geom_point()+geom_line()+theme_bw()+
  labs(y=bquote("AWI  mean monthly precipitation (kg"~s^-1~m^-2~")"),x="Age (yr BP)")

############################################ TR5AS ################################################
tas_TR5AS_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR5AS_Iberian.csv",sep=""),row.names=1)
pr_TR5AS_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR5AS_Iberian.csv",sep=""),row.names=1)

#Temperature
tas_Iberian<-tas_TR5AS_Iberian[,-c(1,2)]-273.15 #the unit is K, convert to degree celcius
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
df3<-cbind.data.frame(age,age_100yrbin,month,colMeans(tas_Iberian),colSds(as.matrix(tas_Iberian)))
colnames(df3)<-c("age","age_100yrbin","month","T","sd_T")

df3_winter<-df3[which(df3$month==12|df3$month==1|df3$month==2),]
df3_winter<-aggregate(df3_winter,by=list(df3_winter$age_100yrbin),FUN=mean)
ggplot(df3_winter,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="TR5AS  MTCO (¡ãC)",x="Age (yr BP)")
df3_summer<-df3[which(df3$month==6|df3$month==7|df3$month==8),]
df3_summer<-aggregate(df3_summer,by=list(df3_summer$age_100yrbin),FUN=mean)
ggplot(df3_summer,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="TR5AS  MTWA (¡ãC)",x="Age (yr BP)")

#Precipitation
pr_Iberian<-pr_TR5AS_Iberian[,-c(1,2)]
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
df3_pr<-cbind.data.frame(age,age_100yrbin,month,colMeans(pr_Iberian),colSds(as.matrix(pr_Iberian)))
colnames(df3_pr)<-c("age","age_100yrbin","month","prep","sd_prep")

df3_pr<-aggregate(df3_pr,by=list(df3_pr$age_100yrbin),FUN=mean)
ggplot(df3_pr,aes(age_100yrbin,prep))+geom_point()+geom_line()+theme_bw()+
  labs(y=bquote("TR5AS  mean monthly precipitation (kg"~s^-1~m^-2~")"),x="Age (yr BP)")

############################################ TR6AV #########################################
tas_TR6AV_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/tas_TR6AV_Iberian.csv",sep=""),row.names=1)
pr_TR6AV_Iberian<-read.csv(paste(wd,"/Master Project/Data/Output data/Core/Results to compare/PACMEDY/pr_TR6AV_Iberian.csv",sep=""),row.names=1)

#Temperature
tas_Iberian<-tas_TR6AV_Iberian[,-c(1,2)]-273.15 #the unit is K, convert to degree celcius
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
df4<-cbind.data.frame(age,age_100yrbin,month,colMeans(tas_Iberian),colSds(as.matrix(tas_Iberian)))
colnames(df4)<-c("age","age_100yrbin","month","T","sd_T")

df4_winter<-df4[which(df4$month==12|df4$month==1|df4$month==2),]
df4_winter<-aggregate(df4_winter,by=list(df4_winter$age_100yrbin),FUN=mean)
ggplot(df4_winter,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="TR6AV  MTCO (¡ãC)",x="Age (yr BP)")
df4_summer<-df4[which(df4$month==6|df4$month==7|df4$month==8),]
df4_summer<-aggregate(df4_summer,by=list(df4_summer$age_100yrbin),FUN=mean)
ggplot(df4_summer,aes(age_100yrbin,T))+geom_point()+geom_line()+theme_bw()+
  labs(y="TR6AV  MTWA (¡ãC)",x="Age (yr BP)")

#Precipitation
pr_Iberian<-pr_TR6AV_Iberian[,-c(1,2)]
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
df4_pr<-cbind.data.frame(age,age_100yrbin,month,colMeans(pr_Iberian),colSds(as.matrix(pr_Iberian)))
colnames(df4_pr)<-c("age","age_100yrbin","month","prep","sd_prep")

df4_pr<-aggregate(df4_pr,by=list(df4_pr$age_100yrbin),FUN=mean)
ggplot(df4_pr,aes(age_100yrbin,prep))+geom_point()+geom_line()+theme_bw()+
  labs(y=bquote("TR6AV  mean monthly precipitation (kg"~s^-1~m^-2~")"),x="Age (yr BP)")

################################# Plot them together
#1 MPI; 2 AWI; 3 TR5AS; 4 TR6AV
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots"))
p1<-ggplot()+theme_bw()+
  geom_point(data=df1_winter,aes(age_100yrbin,T))+geom_line(data=df1_winter,aes(age_100yrbin,T))+  
  geom_point(data=df2_winter,aes(age_100yrbin,T),col="red")+geom_line(data=df2_winter,aes(age_100yrbin,T),col="red")+
  geom_point(data=df3_winter,aes(age_100yrbin,T),col="blue")+geom_line(data=df3_winter,aes(age_100yrbin,T),col="blue")+
  geom_point(data=df4_winter,aes(age_100yrbin,T),col="gold2")+geom_line(data=df4_winter,aes(age_100yrbin,T),col="gold2")+
  labs(y="MTCO (¡ãC)",x=NULL)+annotate("text",x=0,y=10,label="(a)")+
  theme(axis.text.x=element_blank())

p2<-ggplot()+theme_bw()+
  geom_point(data=df1_summer,aes(age_100yrbin,T))+geom_line(data=df1_summer,aes(age_100yrbin,T))+  
  geom_point(data=df2_summer,aes(age_100yrbin,T),col="red")+geom_line(data=df2_summer,aes(age_100yrbin,T),col="red")+
  geom_point(data=df3_summer,aes(age_100yrbin,T),col="blue")+geom_line(data=df3_summer,aes(age_100yrbin,T),col="blue")+
  geom_point(data=df4_summer,aes(age_100yrbin,T),col="gold2")+geom_line(data=df4_summer,aes(age_100yrbin,T),col="gold2")+
  labs(y="MTWA (¡ãC)",x=NULL)+annotate("text",x=0,y=23,label="(b)")+
  theme(axis.text.x=element_blank())

p3<-ggplot()+theme_bw()+
  geom_point(data=df1_pr,aes(age_100yrbin,prep*10^5))+geom_line(data=df1_pr,aes(age_100yrbin,prep*10^5))+  
  geom_point(data=df2_pr,aes(age_100yrbin,prep*10^5),col="red")+geom_line(data=df2_pr,aes(age_100yrbin,prep*10^5),col="red")+
  geom_point(data=df3_pr,aes(age_100yrbin,prep*10^5),col="blue")+geom_line(data=df3_pr,aes(age_100yrbin,prep*10^5),col="blue")+
  geom_point(data=df4_pr,aes(age_100yrbin,prep*10^5),col="gold2")+geom_line(data=df4_pr,aes(age_100yrbin,prep*10^5),col="gold2")+
  labs(y=bquote("Mean daily precipitation ("~10^-5~kg~s^-1~m^-2~")"),x="Age (yr BP)")+
  annotate("text",x=0,y=2.6,label="(c)")

p<-ggarrange(p1,p2,p3,ncol=1)
ggsave(file="PACMEDY simulations.png",p,width=6,height=10)

######################################### Get the west-east gradient plot #######################################3
if(!require(reshape2)){ install.packages("reshape2");library(reshape2)}
############## MPI
pr_Iberian<-pr_MPI_Iberian
age<-rev(rep(seq(101:7950),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),7850)

lat<-pr_MPI_Iberian[,2];summary(lat)
lon<-pr_MPI_Iberian[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
pr_Iberian$x<-lon1
pr_Iberian$x<-as.factor(pr_Iberian$x);levels(pr_Iberian$x)

df1_pr<-cbind.data.frame(age,age_100yrbin,month)
colnames(df1_pr)<-c("age","age_100yrbin","month")

for(k in 1:nlevels(pr_Iberian$x)){
  pr_Iberian_k<-pr_Iberian[which(pr_Iberian$x==levels(pr_Iberian$x)[k]),]
  df1_pr_k<-as.data.frame(colMeans(pr_Iberian_k[,-c(1,2)]))
  colnames(df1_pr_k)<-levels(pr_Iberian$x)[k]
  df1_pr<-cbind.data.frame(df1_pr,df1_pr_k)
}

df1_pr<-aggregate(df1_pr,by=list(df1_pr$age_100yrbin),FUN=mean)
df1_pr$difference<-df1_pr$`-9.375`-df1_pr$`3.75`

################# AWI
pr_Iberian<-pr_AWI_Iberian
age<-rev(rep(seq(0:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6001)

lat<-pr_AWI_Iberian[,2];summary(lat)
lon<-pr_AWI_Iberian[,1];summary(lon)
lon1<-lon;for(i in 1:length(lon)){if(lon[i]>180){lon1[i]<-lon[i]-360}};summary(lon1)
pr_Iberian$x<-lon1
pr_Iberian$x<-as.factor(pr_Iberian$x);levels(pr_Iberian$x)

df2_pr<-cbind.data.frame(age,age_100yrbin,month)
colnames(df2_pr)<-c("age","age_100yrbin","month")

for(k in 1:nlevels(pr_Iberian$x)){
  pr_Iberian_k<-pr_Iberian[which(pr_Iberian$x==levels(pr_Iberian$x)[k]),]
  df2_pr_k<-as.data.frame(colMeans(pr_Iberian_k[,-c(1,2)]))
  colnames(df2_pr_k)<-levels(pr_Iberian$x)[k]
  df2_pr<-cbind.data.frame(df2_pr,df2_pr_k)
}

df2_pr<-aggregate(df2_pr,by=list(df2_pr$age_100yrbin),FUN=mean)
df2_pr$difference<-df2_pr$`-9.375`-df2_pr$`3.75`

################# TR5AS
pr_Iberian<-pr_TR5AS_Iberian
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
pr_Iberian$x<-as.factor(pr_Iberian$x);levels(pr_Iberian$x)

df3_pr<-cbind.data.frame(age,age_100yrbin,month)
colnames(df3_pr)<-c("age","age_100yrbin","month")

for(k in 1:nlevels(pr_Iberian$x)){
  pr_Iberian_k<-pr_Iberian[which(pr_Iberian$x==levels(pr_Iberian$x)[k]),]
  df3_pr_k<-as.data.frame(colMeans(pr_Iberian_k[,-c(1,2)]))
  colnames(df3_pr_k)<-levels(pr_Iberian$x)[k]
  df3_pr<-cbind.data.frame(df3_pr,df3_pr_k)
}

df3_pr<-aggregate(df3_pr,by=list(df3_pr$age_100yrbin),FUN=mean)
df3_pr$difference<-df3_pr$`-7.5`-df3_pr$`3.75`

####################### TR6AV
pr_Iberian<-pr_TR6AV_Iberian
age<-rev(rep(seq(1:6000),each=12))
age_100yrbin<-round(age,digits=-2)
month<-rep(seq(1:12),6000)
pr_Iberian$x<-as.factor(pr_Iberian$x);levels(pr_Iberian$x)

df4_pr<-cbind.data.frame(age,age_100yrbin,month)
colnames(df4_pr)<-c("age","age_100yrbin","month")

for(k in 1:nlevels(pr_Iberian$x)){
  pr_Iberian_k<-pr_Iberian[which(pr_Iberian$x==levels(pr_Iberian$x)[k]),]
  df4_pr_k<-as.data.frame(colMeans(pr_Iberian_k[,-c(1,2)]))
  colnames(df4_pr_k)<-levels(pr_Iberian$x)[k]
  df4_pr<-cbind.data.frame(df4_pr,df4_pr_k)
}

df4_pr<-aggregate(df4_pr,by=list(df4_pr$age_100yrbin),FUN=mean)
df4_pr$difference<-df4_pr$`-10`-df4_pr$`5`

p<-ggplot()+theme_bw()+
  geom_point(data=df1_pr,aes(age_100yrbin,difference/10^(-5)))+geom_line(data=df1_pr,aes(age_100yrbin,difference/10^(-5)))+  
  geom_point(data=df2_pr,aes(age_100yrbin,difference/10^(-5)),col="red")+geom_line(data=df2_pr,aes(age_100yrbin,difference/10^(-5)),col="red")+
  geom_point(data=df3_pr,aes(age_100yrbin,difference/10^(-5)),col="blue")+geom_line(data=df3_pr,aes(age_100yrbin,difference/10^(-5)),col="blue")+
  geom_point(data=df4_pr,aes(age_100yrbin,difference/10^(-5)),col="gold2")+geom_line(data=df4_pr,aes(age_100yrbin,difference/10^(-5)),col="gold2")+
  labs(y=bquote("The difference between the westmost and eastmost mean daily precipitation ("~10^-5~kg~s^-1~m^-2~")"),x="Age (yr BP)")
ggsave(p,file="PACMEDY simulated west-east prep difference.jpeg",width=8,height=7)

#######################################################################################################
##################################   Check Basa de la Mora   #######################################
#######################################################################################################
setwd(paste(wd,"/Master Project/Data/Output data/Core plots/Iberian plots",sep=""))

plotdata1<-cbind.data.frame(core_sig,sse_core_sig)
plotdata1<-plotdata1[which(plotdata1$alpha>=0&plotdata1$alpha<=1.26),]
sitename<-"Basa de la Mora"
xbreak<-2000*c(seq(0,6))
plotsite<-plotdata1[which(plotdata1$site==sitename),]

lon<-unique(plotsite$lon);  lat<-unique(plotsite$lat);  elv<-unique(plotsite$elv)

p_Tmax<-ggplot(data=plotsite)+theme_bw()+geom_point(aes(age,Tmax))+geom_line(aes(age,Tmax))+
  geom_ribbon(aes(x=age,ymin=Tmax-sse_Tmax*1.96,ymax=Tmax+sse_Tmax*1.96),alpha=0.36)+
  labs(y= "Reconstructed MTWA (¡ãC)", x = "Age (yr BP)")+
  scale_y_continuous(labels = function(x) sprintf("%g", x))+
  scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))

ggsave(file=paste("Lon",round(lon,digits=2),"Lat",round(lat,digits=2),sitename,"MTWA reconstruction results.jpeg"),p_Tmax,width=8,height=6)
