# Get insolation --------------------------------------------------------
get_insolation<-function(Lat,age_start,age_end,interval){
  if(!require(palinsol)){ install.packages("palinsol");library(palinsol)}
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  
  # get orbital parameters for past 
  time_BP <- seq(-age_start,-age_end,by=interval)  
  orbital_params <- data.frame(time_BP, t(sapply(time_BP, function(tt) astro(tt, ber78, degree=FALSE))))
  
  # define output variables
  length(time_BP)
  insol_month <- matrix(0, nrow=length(time_BP), ncol=12)
  
  # generate present-day mid-month day numbers
  mid_month_days_0ka <- seq(15.5, 345.5, by=30)
  mid_month_days_0ka
  
  # get true solar longitudes (TSL) for present-day mid-month day numbers
  tt_present = 0.0
  orbit_present <- astro(tt_present, ber78, degree = FALSE)
  mid_month_tsl_0ka <- day2l(orbit_present, mid_month_days_0ka) #/(pi/180)
  
  # mid-month insolation values
  for (month in seq(1:12)) {
    tsl <- mid_month_tsl_0ka[month] 
    insol_month[,month] <- as.numeric(Insol(orbital_params, long=tsl, lat=Lat*pi/180, S0=1365))
  }
  colnames(insol_month) <- c( "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  if(nrow(insol_month)==1){
    insol_summer<-mean(insol_month[,c("Jun","Jul","Aug")])
    insol_winter<-mean(insol_month[,c("Dec","Jan","Feb")])
  }else{
    insol_summer<-rowMeans(insol_month[,c("Jun","Jul","Aug")])
    insol_winter<-rowMeans(insol_month[,c("Dec","Jan","Feb")])
  }
  
  insol_out <- data.frame(-time_BP, insol_summer,insol_winter,insol_month)
  colnames(insol_out)<-c("age","insol_summer","insol_winter",
                         "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  return(insol_out)
  
}

#Plot the core results -------------------------
plot.core.sig.sse<-function(sitename, core_sig,sse_core_sig,zero_Tmin,zero_gdd,zero_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  xbreak<-2000*c(seq(0,6))
  plotsite<-plotdata[which(plotdata$site==sitename),]
  lon<-unique(plotsite$lon);  lat<-unique(plotsite$lat);  elv<-unique(plotsite$elv)
  min_insol_winter<-min(plotsite$insol_winter,na.rm=TRUE)
  min_insol_summer<-min(plotsite$insol_summer,na.rm=TRUE)
  
  coeff_Tmin<-(max(plotsite$insol_summer)-min(plotsite$insol_summer))/(max(plotsite$Tmin)-min(plotsite$Tmin))
  p_Tmin<-ggplot(data=plotsite)+theme_bw()+geom_point(aes(age,Tmin))+geom_line(aes(age,Tmin))+
    geom_ribbon(aes(x=age,ymin=Tmin-sse_Tmin*1.96,ymax=Tmin+sse_Tmin*1.96),alpha=0.36)+
    geom_point(aes(age,(insol_winter-min_insol_winter)/coeff_Tmin),col="red")+
    geom_line(aes(age,(insol_winter-min_insol_winter)/coeff_Tmin),col="red")+
    labs(y= "Reconstructed MTCO (°C)", x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(labels = function(x) sprintf("%g", x),
                       sec.axis = sec_axis(~(.*coeff_Tmin+min_insol_winter), name=bquote("Winter insolation (W m"^-2*")")))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    ggtitle(paste(sitename,"; Lon:",lon,"; Lat:",lat,"; Elv:",elv))
  
  coeff_gdd<-(max(plotsite$insol_summer)-min(plotsite$insol_summer))/(max(plotsite$gdd)-min(plotsite$gdd))
  p_gdd<-ggplot(data=plotsite)+theme_bw()+geom_point(aes(age,gdd))+geom_line(aes(age,gdd))+
    geom_ribbon(aes(x=age,ymin=gdd-sse_gdd*1.96,ymax=gdd+sse_gdd*1.96),alpha=0.36)+
    geom_point(aes(age,(insol_summer-min_insol_summer)/coeff_gdd),col="red")+
    geom_line(aes(age,(insol_summer-min_insol_summer)/coeff_gdd),col="red")+
    labs(y= bquote('Reconstructed '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(labels = function(x) sprintf("%g", x),
                       sec.axis = sec_axis(~(.*coeff_gdd+min_insol_summer), name=bquote("Summer insolation (W m"^-2*")")))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))
  
  coeff_alpha<-(max(plotsite$insol_summer)-min(plotsite$insol_summer))/(max(plotsite$alpha)-min(plotsite$alpha))
  p_alpha<-ggplot(data=plotsite)+theme_bw()+geom_point(aes(age,alpha))+geom_line(aes(age,alpha))+
    geom_ribbon(aes(x=age,ymin=alpha-sse_alpha*1.96,ymax=alpha+sse_alpha*1.96),alpha=0.36)+
    geom_point(aes(age,(insol_summer-min_insol_summer)/coeff_alpha),col="red")+
    geom_line(aes(age,(insol_summer-min_insol_summer)/coeff_alpha),col="red")+
    labs(y= expression("Reconstructed "*alpha),x = "Age (yr BP)")+
    scale_y_continuous(labels = function(x) sprintf("%0.2f", x),
                       sec.axis = sec_axis(~(.*coeff_alpha+min_insol_summer), name=bquote("Summer insolation (W m"^-2*")")))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))
  
  p<-ggarrange(p_Tmin ,p_gdd, p_alpha, ncol = 1,nrow=3)
  ggsave(file=paste("Lon",round(lon,digits=2),"Lat",round(lat,digits=2),sitename,"reconstruction results.jpeg"),p,width=10,height=9)
  gc()
  
}
# Get MTWA (mean temperature of the warmest month) from GDD0 and MTCO --------------------------------------------------------
get_MTWA<-function(MTCO,GDD0){
  if(!is.na(MTCO)&!is.na(GDD0)){
    if(MTCO<0){
      GDD0_rad<-GDD0*2*pi/365
      if(!require(rootSolve)){ install.packages("rootSolve");library(rootSolve)}
      #- GDD0 /Tmin = 2 {[(u/(1 - u)] cos-1 (-u) + ???[(1 + u)/(1 ¬- u)]}.
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

# Get the gradient --------------------------------------------------------
get_gradient<-function(plotdata,age,range){
  plotdata$site<-as.factor(plotdata$site)
  mean_env<-data.frame(matrix(NA,ncol=5+3+2,nrow=nlevels(plotdata$site)))
  for(m in 1:nlevels(plotdata$site)){
    plotdata_each<-plotdata[which(plotdata$site==levels(plotdata$site)[m]),]
    
    mean_Tmin<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"Tmin"]),na.rm = TRUE)
    mean_gdd<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"gdd"]),na.rm = TRUE)
    mean_alpha<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"alpha"]),na.rm = TRUE)
    
    mean_insol_summer<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"insol_summer"]),na.rm = TRUE)
    mean_insol_winter<-mean(as.matrix(plotdata_each[which(plotdata_each$age>age-range & plotdata_each$age<age+range),"insol_winter"]),na.rm = TRUE)
    mean_env[m,]<-c(age,levels(plotdata$site)[m],unique(plotdata_each$lon),unique(plotdata_each$lat),unique(plotdata_each$elv),mean_Tmin,mean_gdd,mean_alpha,mean_insol_summer,mean_insol_winter)
  }
  colnames(mean_env)<-c("age","site","lon","lat","elv","mean_Tmin","mean_gdd","mean_alpha","mean_insol_summer","mean_insol_winter")
  return(mean_env)
}
