#Changepoint analysis of Motus detections to determine the onset and end of diel activity
#For project 882 - Montreal Northern Cardinal Winter Movements
#Created by: Andrew Beauchamp
#R version 4.4.2
#Last updated 2025-04-28

Sys.setenv(TZ = "UTC")

getwd()
# rm(list=ls())
# gc()

#Installing packages and loading libraries####
# install.packages("remotes")
# library(remotes)
# remotes::update_packages()
# install.packages("Rtools")

# install.packages("motus",repos = c(birdscanada = 'https://birdscanada.r-universe.dev',CRAN = 'https://cloud.r-project.org'))
# install.packages("changepoint")

library(motus)
# packageVersion("motus")

library(dplyr)
library(lubridate)
library(ggplot2)
library(suncalc)
library(evaluate)
library(changepoint)
library(plotly)


#Downloading Motus data preparing time variables####
getwd()
# sql_motus <- tagme(projRecv = 882, new = TRUE, update=TRUE, dir = "./data/") #last run 2025-04-29
# Andrew_Beauchamp
# Motus_atb_WTS_1066-1944

motusLogout()
#sql_motus <- tagme(projRecv = 882, new = FALSE, update=TRUE, dir = "./data/") #last run 2025-04-29
sql_motus <- tagme(projRecv = 882, new = FALSE, update=FALSE, dir = "./data/")

tbl.alltags <- tbl(sql_motus, "alltags")

tbl.alltags %>% filter(motusTagID==68906) %>% dplyr::select(tagBI,speciesEN)


##Filtering to Montreal region, collecting as RDS####
df.proj882 <- tbl.alltags %>%
  filter(recvDeployLon>-74.2796,
         recvDeployLon<-71.9254,
         recvDeployLat>44.9655,
         recvDeployLat<46.132) %>% #filter to towers in Montreal region
  collect() %>%
  as.data.frame() %>%
  mutate(time_UTC = as_datetime(ts))

saveRDS(df.proj882, "./data/df_proj882_Mo.rds")


##Inspecting filtered data####
df.proj882<-readRDS("./data/df_proj882_Mo.rds")

#issue with missing metadata...needs at least tagBI for changepoint analysis
df.proj882 %>% filter(motusTagID==68906) %>% head()

#infilling metadata separately - done by batch now


glimpse(df.proj882) #71,612,279 rows
unique(df.proj882$speciesEN) #NA species name due to missing metadata


min(df.proj882$recvDeployLat)
max(df.proj882$recvDeployLat)

min(df.proj882$recvDeployLon)
max(df.proj882$recvDeployLon)

unique(df.proj882$recvDeployName)


###Count of unique bird-station-date - Not run####
#Computationally intense - commented out to avoid re-running
# glimpse(df.proj882)
# df.proj882a<-df.proj882 %>% mutate(days1970=floor(ts/60/60/24))
# recvdays<-df.proj882a %>% dplyr::select(motusTagID,recvDeployName,days1970)%>% 
#   unique()%>%
#   count(motusTagID,recvDeployName)
# 
# daysrecv<-df.proj882a %>% dplyr::select(motusTagID,recvDeployName,days1970)%>% 
#   unique()%>%
#   count(motusTagID,days1970)
# 
# daysrecv_n<-daysrecv %>% count(motusTagID, n,name="days")%>%rename(n.recv=n)
# 
# write.csv(recvdays,file="Proj882-recieverdays.csv")
# write.csv(daysrecv_n,file="Proj882-recieverdays_n.csv")
# 



#Getting sunrise and sunset for each day in the study####
glimpse(df.proj882)#Rows: 71,612,279

start<-date(min(df.proj882$time_UTC))
end<-date(max(df.proj882$time_UTC)) 

Dates<-seq(start,end, by = '1 day')

Montreal.sunriset<-data.frame(Dates)
Montreal.sunriset$sunrise.UTC <- getSunlightTimes(date=as.Date(Montreal.sunriset$Date),
                                             lon=-73.9385,
                                             lat=45.4307,
                                             tz="UTC",
                                             keep=c("sunrise"))$sunrise

Montreal.sunriset$sunset.UTC <- getSunlightTimes(date=as.Date(Montreal.sunriset$Date),
                                            lon=-73.9385,
                                            lat=45.4307,
                                            tz="UTC",
                                            keep=c("sunset"))$sunset
glimpse(Montreal.sunriset)

saveRDS(Montreal.sunriset,file="Montreal.sunriset.rds")

Montreal.sunriset<-readRDS("Montreal.sunriset.rds")

##Adding sunrise and sunset times to Motus data for filtering####
glimpse(Montreal.sunriset)
glimpse(df.proj882)#71,612,279 rows

df.proj882$Dates<-date(df.proj882$time_UTC)
df.proj882.MO<-left_join(df.proj882,Montreal.sunriset,by="Dates")

glimpse(df.proj882.MO)

nrow(df.proj882.MO)==nrow(df.proj882)
any(is.na(df.proj882.MO$sunrise.UTC)==TRUE)
head(df.proj882.MO)
rm(df.proj882)
gc()

saveRDS(df.proj882.MO,file="df.proj882.MO.rds")

rm(list=ls())
gc()

##

#Creating seperate RDS files for each batch####
#71,612,279 rows


Montreal.sunriset<-readRDS("Montreal.sunriset.rds")
df.proj882.MO<-readRDS("df.proj882.MO.rds")


length(unique(df.proj882.MO$motusTagID)) # Batches of 5 tags


taglist<-unique(df.proj882.MO$motusTagID) %>% sort()
taglist.a<-taglist[1:5]
taglist.b<-taglist[6:10]
taglist.c<-taglist[11:15]
taglist.d<-taglist[16:20]
taglist.e<-taglist[21:30]
taglist.f<-taglist[31:40]
taglist.g<-taglist[41:50]
taglist.h<-taglist[51:60]
taglist.i<-taglist[61:70]
taglist.j<-taglist[71:80]
taglist.k<-taglist[81:90]
taglist.l<-taglist[91:100]


df.proj882.MO.a<-df.proj882.MO %>% filter(motusTagID %in% taglist.a)
head(df.proj882.MO.a)
saveRDS(df.proj882.MO.a,file="df.proj882.MO.a.rds")
rm(df.proj882.MO.a)
gc()




df.proj882.MO.b<-df.proj882.MO %>% filter(motusTagID %in% taglist.b)
head(df.proj882.MO.b)
saveRDS(df.proj882.MO.b,file="df.proj882.MO.b.rds")
rm(df.proj882.MO.b)
gc()


df.proj882.MO.c<-df.proj882.MO %>% filter(motusTagID %in% taglist.c)
head(df.proj882.MO.c)
saveRDS(df.proj882.MO.c,file="df.proj882.MO.c.rds")
rm(df.proj882.MO.c)
gc()


df.proj882.MO.d<-df.proj882.MO %>% filter(motusTagID %in% taglist.d)
head(df.proj882.MO.d)
saveRDS(df.proj882.MO.d,file="df.proj882.MO.d.rds")
rm(df.proj882.MO.d)
gc()

df.proj882.MO.e<-df.proj882.MO %>% filter(motusTagID %in% taglist.e)
head(df.proj882.MO.e)
saveRDS(df.proj882.MO.e,file="df.proj882.MO.e.rds")
rm(df.proj882.MO.e)
gc()

df.proj882.MO.f<-df.proj882.MO %>% filter(motusTagID %in% taglist.f)
head(df.proj882.MO.f)
saveRDS(df.proj882.MO.f,file="df.proj882.MO.f.rds")
rm(df.proj882.MO.f)
gc()

df.proj882.MO.g<-df.proj882.MO %>% filter(motusTagID %in% taglist.g)
head(df.proj882.MO.g)
saveRDS(df.proj882.MO.g,file="df.proj882.MO.g.rds")
rm(df.proj882.MO.g)
gc()

df.proj882.MO.h<-df.proj882.MO %>% filter(motusTagID %in% taglist.h)
head(df.proj882.MO.h)
saveRDS(df.proj882.MO.h,file="df.proj882.MO.h.rds")
rm(df.proj882.MO.h)
gc()

df.proj882.MO.i<-df.proj882.MO %>% filter(motusTagID %in% taglist.i)
head(df.proj882.MO.i)
saveRDS(df.proj882.MO.i,file="df.proj882.MO.i.rds")
rm(df.proj882.MO.i)
gc()

df.proj882.MO.j<-df.proj882.MO %>% filter(motusTagID %in% taglist.j)
head(df.proj882.MO.j)
saveRDS(df.proj882.MO.j,file="df.proj882.MO.j.rds")
rm(df.proj882.MO.j)
gc()

df.proj882.MO.k<-df.proj882.MO %>% filter(motusTagID %in% taglist.k)
head(df.proj882.MO.k)
saveRDS(df.proj882.MO.k,file="df.proj882.MO.k.rds")
rm(df.proj882.MO.k)
gc()

df.proj882.MO.l<-df.proj882.MO %>% filter(motusTagID %in% taglist.l)
head(df.proj882.MO.l)
saveRDS(df.proj882.MO.l,file="df.proj882.MO.l.rds")
rm(df.proj882.MO.l)
gc()

##RDS FILELIST#####
rdsfilelist<-c("df.proj882.MO.a.rds",
               "df.proj882.MO.b.rds",
               "df.proj882.MO.c.rds",
               "df.proj882.MO.d.rds",
               "df.proj882.MO.e.rds",
               "df.proj882.MO.f.rds",
               "df.proj882.MO.g.rds",
               "df.proj882.MO.h.rds",
               "df.proj882.MO.i.rds",
               "df.proj882.MO.j.rds",
               "df.proj882.MO.k.rds",
               "df.proj882.MO.l.rds")


#NOT DONE - Visualizing activity for each bird across all days of stopover####
# #using a loop to call each separate file, generate plot, then remove data
# 
# rm(df.proj882)
# gc()
# 
##loop for each batch####
# starttime<-Sys.time()
# 
# for (k in 1:length(rdsfilelist)){
#     # k<-7  #for testing
#   filedata<-readRDS(rdsfilelist[k])
# # head(filedata)
# 
##Plotting signal ~ time for each bird-station combination####
# filedata<-filedata%>%arrange(motusTagID,ts)
# 
# filedata$motustagIDstation<-paste0(filedata$motusTagID,"-",filedata$recvDeployName)
# # unique(filedata$motustagIDstation)
# motustagIDlist<-unique(filedata$motustagIDstation)
# # length(motustagIDlist)
# # motustagIDlist<-c(unique(df.proj551.MBO3$motusTagID)[1:5],68698)     #FOR TESTING
# #	i<-1                                                                #FOR TESTING
# 
# 
# #list of plots
# signalplots=list()
# gc()
# 
###Loop to generate plots####
# for(i in 1:length(motustagIDlist)){
#   
#   targetbird<-motustagIDlist[i] #select target bird
#   
#   bird.loop.data<- filedata %>%
#     dplyr::filter(motustagIDstation == targetbird)
#   
#   sunriset.loop<-Montreal.sunriset %>%
#     filter(sunrise.UTC > (min(bird.loop.data$time_UTC)-86400) &  sunset.UTC < (max(bird.loop.data$time_UTC)+86400))
#   
#   signalplots[[i]]<-ggplot() +
#     geom_rect(data=sunriset.loop,
#               aes(xmin=sunrise.UTC,
#                   xmax=sunset.UTC,
#                   ymin=-Inf,
#                   ymax=Inf),
#               fill="white")+
#     ylab("Signal strength (dB)")+
#     xlab("Date (UTC)")+
#     xlim(c(min(sunriset.loop$sunrise.UTC),max(sunriset.loop$sunset.UTC)))+
#     ggtitle(paste0("MotusTagID - station: ", targetbird))+
#     geom_point(data=bird.loop.data,aes(x = time_UTC, y = sig, col = as.factor(port)))+
#     scale_color_discrete(name = "Antenna\nport") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background=element_rect(fill="grey80"),
#       panel.border = element_rect(colour = "black", fill=NA, size = 0.5),
#       plot.margin=unit(c(0.2,0.3,0.2,0.2), "cm"),
#       axis.text.x = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
#       axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
#       axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
#       axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain")
#     )
#   rm(bird.loop.data,sunriset.loop)
#   gc()
# }
# 
# 
# #saving - takes ~ 10 mins
# nplotfile<-8
# for(j in 1:ceiling((length(signalplots)/nplotfile))) {
#   sigplotj<-signalplots[(1+((j-1)*nplotfile)):(nplotfile+((j-1)*nplotfile))]
#   graphics.off()
#   pdf(paste0("./plots/Project882tag-stationplots","batch_",k,"_group",j,".pdf"),
#       width=24,
#       height=8,
#       onefile=TRUE)
#   for(l in sigplotj) {replay(l)
#     rm(sigplotj)
#   }
#   graphics.off()
# }
# 
# 
# rm(signalplots)
# rm(filedata)
# gc()
# 
# }
# 
# Sys.time()-starttime



#Change point analysis####
### Automated changepoint detection method ###
# Prepare data frame for analysis.
# For each port on each station, calculate lagged differences in sig between consecutive detections.
# Keep diffs only if consecutive detections occur with a short time lag.

##Onset of activity####
###Setting buffer around sunrise, and filtering#####
sunrisbuffer <- 3 # buffer (in hours) that the detection must be prior to or following sunrise/sunset

###Initiating loop for each batch#####
column.names <- c("batch","elapsed.t")
elapse.t.onset <- array(NA,dim=c(12,2))
dimnames(elapse.t.onset) <- list(NULL,column.names)
elapse.t.onset

metadata<-read.csv("./data/MotusTagsTemplate-882.csv")
metadata<-metadata %>% mutate(motusTagID=tagID,tagBI=period) %>% dplyr::select(motusTagID,tagBI)%>%unique()

start.t<-Sys.time()
for (k in 1:length(rdsfilelist)){
  # k<-1  #for testing
  filedata<-readRDS(rdsfilelist[k])


#infilling missing metadata
filedata.metadata<-filedata %>% filter(is.na(tagBI)==FALSE)
filedata.nometadata<-filedata %>% filter(is.na(tagBI)==TRUE) %>% dplyr::select(-c(tagBI))
filedata.nometadata<-left_join(filedata.nometadata,metadata,by="motusTagID")
filedata<-rbind(filedata.metadata,filedata.nometadata)


####Generate the variable to determine changepoints#####
#sqrt(abs(siglag-sig)) 
glimpse(filedata)
filedata$Dates

filedata2 <- filedata %>%
  group_by(motusTagID) %>%
  mutate(dos = as.numeric(Dates-min(Dates))+1) %>% #dos= Day of study, from 1 to max days for each bird
  ungroup() %>%
  group_by(motusTagID,recvDeployName,port) %>% 
  arrange(ts)%>% 							
  mutate(siglag=dplyr::lag(sig),
         tslag=dplyr::lag(ts),
         hour=hour(as.POSIXct(ts,tz="UTC"))+minute(as.POSIXct(ts,tz="UTC"))/60) %>%  	
  mutate(diff.abs=sqrt(abs(siglag-sig)),
         diff.ts=ts-tslag) %>% 		
  filter(diff.ts < (tagBI*2.2)) %>% #drops observations with a difference between consecutive detections greater than 2 times the tag burst interval + 10%
  as.data.frame()	


rm(filedata)
gc()

# glimpse(filedata2)	#FOR TESTING
# length(unique(filedata2$motusTagID)) #FOR TESTING

#validating
range(filedata2$diff.ts)# some with 0 seconds between bursts received on the same station and port...not possible

###Filter to remove false detections####
##just removing detections with a diff.ts less than the tagBI
filedata3<-filedata2 %>%
  filter(round(diff.ts,1)>=round(tagBI,1)) 

rm(filedata2)

###Dropping detections outside of buffer####
filedata3$sunrisediff<-as.numeric(difftime(filedata3$time_UTC,filedata3$sunrise.UTC, units="hours"))

batchonsetdat <- filedata3 %>%
  filter(sunrisediff> -sunrisbuffer &  #drops observations earlier than 3hrs before sunrise
           sunrisediff < sunrisbuffer) %>% #drops observations later than 3hrs after sunrise 
  as.data.frame()

# glimpse(batchonsetdat)	#FOR TESTING #Rows: 4,143,132
# length(unique(batchonsetdat$motusTagID)) #FOR TESTING #113 birds, 10 birds lacking detections around sunrise

#verifying
#  range(batchonsetdat$sunrisediff) #FOR TESTING - GOOD


###Loop to determine onset of activity#####
###list of birds with good signals

bird.list <- sort(unique(batchonsetdat$motusTagID))

batchonsetdat <- batchonsetdat %>% arrange(motusTagID, ts)   # Inter-digitate diff.abs values from all ports and stations for one sequence.

crit.nobs <- 800         # number of observations required in focal period to do the changepoint analysis.

dim.j <- max(batchonsetdat$dos)         # max number of days to look for changepoints for each bird
dim.i <- length(bird.list)              # number of unique bird ids

# initialize arrays

cpt.ts1.array <- array(NA,dim=c(dim.i,dim.j))	# hold changepoints
cpt.tsmax <- array(NA,dim=c(dim.i,dim.j))	  # hold maximum ts for that bird
cpt.tspre <- array(NA,dim=c(dim.i,dim.j))       # hold ts of changepoint - 1 observation
cpt.tspost <- array(NA,dim=c(dim.i,dim.j))      # hold ts of changepoint + 1 observation
cpt.meanpre <- array(NA,dim=c(dim.i,dim.j))     # hold mean of diff.abs before the changepoint
cpt.meanpost <- array(NA,dim=c(dim.i,dim.j))    # hold mean of diff.abs after the changepoint
cpt.varpre <- array(NA,dim=c(dim.i,dim.j))      # hold var of diff.abs before the changepoint
cpt.varpost <- array(NA,dim=c(dim.i,dim.j))     # hold var of diff.abs after the changepoint
cpt.motusTagID <- array(NA,dim=c(dim.i,dim.j))     	# hold id of bird

# assign names to arrays
dimnames(cpt.ts1.array) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.tsmax) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.tspre) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.tspost) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.meanpre) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.meanpost) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.varpre) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.varpost) <- list(bird.num=1:dim.i, dos=1:dim.j )
dimnames(cpt.motusTagID) <- list(bird.num=1:dim.i, dos=1:dim.j )

# i<-5 #for testing
# j<-12 #for testing
# rm(i,j) #for testing

for (i in 1:dim.i) {                           
  id <- bird.list[i]
  for (j in 1:dim.j) {
    print(paste("motusTagID=",bird.list[i],"i=",i,"dos=",j))
    temp <- batchonsetdat %>% filter(motusTagID==bird.list[i] & dos==j)													
    n <- length(temp$diff.abs)
    print(paste("n=",n))																							
    if (n > crit.nobs) {
      out1.pmv <- cpt.meanvar(temp$diff.abs,
                              method="AMOC",
                              penalty="Asymptotic",
                              pen.value=0.001,
                              minseglen=20,
                              param.estimates=T,
                              class=T)
      if (ncpts(out1.pmv) > 0) {         
        cpt1 <- cpts(out1.pmv)
        cpt.ts1 <- temp[cpt1,"ts"]
        cpt.ts1.array[i,j] <- cpt.ts1
        cpt.tsmax[i,j] <- max(temp$ts)
        cpt.tspre[i,j]  <- round(cpt.ts1-temp[cpt1-1,"ts"],2)
        cpt.tspost[i,j] <- round(temp[cpt1+1,"ts"]-cpt.ts1,2)
        cpt.meanpre[i,j] <- param.est(out1.pmv)$mean[1]
        cpt.meanpost[i,j] <- param.est(out1.pmv)$mean[2]
        cpt.varpre[i,j] <- param.est(out1.pmv)$variance[1]
        cpt.varpost[i,j] <- param.est(out1.pmv)$variance[2]
        cpt.motusTagID[i,j]<-bird.list[i]
        print(paste(cpt.ts1))
      } else {
        if (ncpts(out1.pmv) == 0) {
          cpt.ts1 <- max(temp$ts)
          cpt.ts1.array[i,j] <- NA
          print(paste(cpt.ts1,"***************************No change points***************************"))
        }}
    }
  }    # end j
}          # end i




# create data frames from the arrays
ts1 <- as.data.frame.table(cpt.ts1.array, responseName = "ts1") 
tsmax <- as.data.frame.table(cpt.tsmax, responseName = "tsmax") 
tspre <- as.data.frame.table(cpt.tspre, responseName = "tspre") 
tspost <- as.data.frame.table(cpt.tspost, responseName = "tspost") 
meanpre <- as.data.frame.table(cpt.meanpre, responseName = "meanpre")
meanpost <- as.data.frame.table(cpt.meanpost, responseName = "meanpost")
varpre <- as.data.frame.table(cpt.varpre, responseName = "varpre")
varpost <- as.data.frame.table(cpt.varpost, responseName = "varpost")
motusTagID<- as.data.frame.table(cpt.motusTagID, responseName = "motusTagID")


# merge the data frames
cpts <- merge(ts1,tspre,by=c("bird.num","dos"))
cpts <- merge(cpts,tsmax,by=c("bird.num","dos"))
cpts <- merge(cpts,tspost,by=c("bird.num","dos"))
cpts <- merge(cpts,meanpre,by=c("bird.num","dos"))
cpts <- merge(cpts,meanpost,by=c("bird.num","dos"))
cpts <- merge(cpts,varpre,by=c("bird.num","dos"))
cpts <- merge(cpts,varpost,by=c("bird.num","dos"))
cpts <- merge(cpts,motusTagID,by=c("bird.num","dos"))


cpts <- cpts %>% 
  mutate(cpt = as_datetime(ts1,tz="UTC")) %>%
  mutate(tsmax.dt = as_datetime(tsmax)) %>% 
  filter(ts1 > 0|is.na(ts1)==F) %>%   # keep change points, omit instances with no change points
  dplyr::select(bird.num,motusTagID,dos,tspre,tspost,tsmax.dt,meanpre,meanpost,varpre,varpost,cpt)


cpts<-cpts %>% arrange (bird.num,dos)
head(cpts,500)
names(cpts)
dim(cpts)                          
head(cpts)

cpts$batch<-k

###Saving changepoint batches####
saveRDS(cpts,file=paste0("./data/cpts_onset_raw_","batch_",k,".rds"))
write.csv(cpts,file=paste0("./data/cpts_onset_raw_","batch_",k,".csv"))

rm(filedata3)
gc()

cpt.ts1.array[i,j]
elapse.t.onset[k,"batch"]<-k
elapse.t.onset[k,"elapsed.t"]<- Sys.time()-start.t

}


###Inspecting onset changepoint metrics####
onset.cpts.all<-rbind(
                  readRDS("./data/cpts_onset_raw_batch_1.rds"),
                  readRDS("./data/cpts_onset_raw_batch_2.rds"),
                  readRDS("./data/cpts_onset_raw_batch_3.rds"),
                  readRDS("./data/cpts_onset_raw_batch_4.rds"),
                  readRDS("./data/cpts_onset_raw_batch_5.rds"),
                  readRDS("./data/cpts_onset_raw_batch_6.rds"),
                  readRDS("./data/cpts_onset_raw_batch_7.rds"),
                  readRDS("./data/cpts_onset_raw_batch_8.rds"),
                  readRDS("./data/cpts_onset_raw_batch_9.rds"),
                  readRDS("./data/cpts_onset_raw_batch_10.rds"),
                  readRDS("./data/cpts_onset_raw_batch_11.rds"),
                  readRDS("./data/cpts_onset_raw_batch_12.rds"))
glimpse(onset.cpts.all) #3,198 change points

onset.cpts.all%>%count(motusTagID,dos) %>%filter(n>1) # Check, should be no rows here

qplot(na.exclude(onset.cpts.all$meanpre),bins=50)  # View distribution of meanpre values  #good	- mean should be low (<1) prior to changepoint for onset	
qplot(na.exclude(onset.cpts.all$meanpost),bins=50) # View distribution of meanpost values #good - mean should be higher (1 to > 1) after changepoint when birds are more active
qplot(na.exclude(onset.cpts.all$meanpre-onset.cpts.all$meanpost),bins=50)  # View distribution of diff in mean (most < 0 for onset)	
qplot(na.exclude(onset.cpts.all$varpre),bins=50)  # View distribution of varpre values    		#good - should be generally low
qplot(na.exclude(onset.cpts.all$varpost),bins=50) # View distribution of varpost values 			#good - should be low, but generally higher than pre-changepoint
qplot(na.exclude(onset.cpts.all$varpre-onset.cpts.all$varpost),bins=50)  # View distribution of diff in var 	#good - should be low or negative

####Filter cpts based on criteria####
crit.mean1 <- 1   # critical meanpre value
crit.mean2 <- 0   # critical meanpost value 


#####View questionable change points for further investigation####
delete <- onset.cpts.all[which(onset.cpts.all$meanpre > crit.mean1 |
                       onset.cpts.all$meanpost < crit.mean2 |						
                       onset.cpts.all$meanpre > onset.cpts.all$meanpost | 
                       onset.cpts.all$varpre > onset.cpts.all$varpost),]
delete %>% arrange(bird.num,dos) %>% select(bird.num,cpt,tsmax.dt,meanpre,meanpost,varpre,varpost)

delete <- which(onset.cpts.all$meanpre > crit.mean1 |
                  onset.cpts.all$meanpost < crit.mean2 |
                  onset.cpts.all$meanpre > onset.cpts.all$meanpost | 
                  onset.cpts.all$varpre > onset.cpts.all$varpost)

nrow(onset.cpts.all)
length(delete)
#530 changepoints out of 3198

#Visualizing changepoint distribution after removals
onset.cpts.all.filtered<-onset.cpts.all[-delete,]
nrow(onset.cpts.all.filtered)

qplot(na.exclude(onset.cpts.all.filtered$meanpre),bins=50)  # View distribution of meanpre values  #good	- mean should be low (<1) prior to changepoint for onset	
qplot(na.exclude(onset.cpts.all.filtered$meanpost),bins=50) # View distribution of meanpost values #good - mean should be higher (1 to > 1) after changepoint when birds are more active
qplot(na.exclude(onset.cpts.all.filtered$meanpre-onset.cpts.all.filtered$meanpost),bins=50)  # View distribution of diff in mean (most < 0 for onset)	
qplot(na.exclude(onset.cpts.all.filtered$varpre),bins=50)  # View distribution of varpre values    		#good - should be generally low
qplot(na.exclude(onset.cpts.all.filtered$varpost),bins=50) # View distribution of varpost values 			#good - should be low, but generally higher than pre-changepoint
qplot(na.exclude(onset.cpts.all.filtered$varpre-onset.cpts.all.filtered$varpost),bins=50)  # View distribution of diff in var 	#good - should be low or negative



###Marking questionable changepoints for later removal#####
onset.cpts.all.2 <- onset.cpts.all
onset.cpts.all.2$postfilter<-"good"
onset.cpts.all.2[delete,"postfilter"]<-"bad"
glimpse(onset.cpts.all.2)          

onset.cpts.all.3 <- onset.cpts.all.2 %>%
  mutate(time = hour(cpt)+minute(cpt+(second(cpt)/60))/60) %>% arrange(time)


###Convert changepoints to minutes since sunset####
glimpse(onset.cpts.all.3)
onset.cpts.all.3$Dates<-date(onset.cpts.all.3$cpt)

onset.cpts.all.3 %>% dplyr::count(motusTagID,dos) %>%filter(n>1) # Check for duplicates, should be no rows here


onset.cpts.all.4<-left_join(onset.cpts.all.3,Montreal.sunriset,by="Dates") %>% dplyr::select(-sunset.UTC)
glimpse(onset.cpts.all.4)

onset.cpts.all.4$cpt.sr<-as.numeric(difftime(onset.cpts.all.4$cpt,onset.cpts.all.4$sunrise.UTC, units="mins")) #time1-time2, negative = prior to sunrise

###Visualizing distribution of onset times for good points

onset.cpts.all.4.good<-onset.cpts.all.4 %>% filter(postfilter=="good")
qplot(onset.cpts.all.4.good$cpt.sr,bins=30) #looks good


###Adding variables for plotting later####
onset.cpts.all.4$postfilter<-as.factor(onset.cpts.all.4$postfilter)

onset.cpts.all.4$postfiltercol<-"black"
onset.cpts.all.4[which(onset.cpts.all.4$postfilter=="bad"),"postfiltercol"]<-"red"

onset.cpts.all.4$postfilterlntyp<-1
onset.cpts.all.4[which(onset.cpts.all.4$postfilter=="bad"),"postfilterlntyp"]<-2
glimpse(onset.cpts.all.4)

###Saving onset CPTS####
saveRDS(onset.cpts.all.4,file="onset.cpts.882.rds")
write.csv(onset.cpts.all.4,file="onset.cpts.882.csv")



##End of activity####
###Setting buffer around sunrise, and filtering#####
sunsetbuffer <- 3 # buffer (in hours) that the detection must be prior to or following sunrise/sunset

###Initiating loop for each batch#####
column.names <- c("batch","elapsed.t")
elapse.t.end <- array(NA,dim=c(12,2))
dimnames(elapse.t.end) <- list(NULL,column.names)
elapse.t.onset


 
for (k in 1:length(rdsfilelist)){
  # k<-10  #for testing
  filedata<-readRDS(rdsfilelist[k])
  
  start.t<-Sys.time()
  
  #infilling missing metadata
  filedata.metadata<-filedata %>% filter(is.na(tagBI)==FALSE)
  filedata.nometadata<-filedata %>% filter(is.na(tagBI)==TRUE) %>% dplyr::select(-c(tagBI))
  filedata.nometadata<-left_join(filedata.nometadata,metadata,by="motusTagID")
  filedata<-rbind(filedata.metadata,filedata.nometadata)
  
  ####Generate the variable to determine changepoints#####
  #sqrt(abs(siglag-sig)) 
  glimpse(filedata)
  
  filedata2 <- filedata %>%
    group_by(motusTagID) %>%
    mutate(dos = as.numeric(Dates-min(Dates))+1) %>% #Day of study from 1 to max days for each bird
    ungroup() %>%
    group_by(motusTagID,recvDeployName,port) %>% 
    arrange(ts)%>% 							
    mutate(siglag=dplyr::lag(sig),
           tslag=dplyr::lag(ts),
           hour=hour(as.POSIXct(ts,tz="UTC"))+minute(as.POSIXct(ts,tz="UTC"))/60) %>%  	
    mutate(diff.abs=sqrt(abs(siglag-sig)),
           diff.ts=ts-tslag) %>% 		
    filter(diff.ts < (tagBI*2.2)) %>% #drops observations with a difference between consecutive detections greater than 2 times the tag burst interval + 10%
    as.data.frame()	
  
  rm(filedata)
  gc()
  
  #validating
  range(filedata2$diff.ts)# some with 0 seconds between bursts received on the same station and port...not possible
  
  ###Filter to remove false detections####
  ##just removing detections with a diff.ts less than the tagBI
  filedata3<-filedata2 %>%
    filter(round(diff.ts,1)>=round(tagBI,1)) 
  
  rm(filedata2)

  ###Dropping detections outside of buffer####
  ####Adding sunset time, accounting for previous night####
  glimpse(filedata3)
  
  filedata3$prev.night<-filedata3$Dates-1
  filedata3 %>% dplyr::select(prev.night,Dates)
  
  Montreal.sunriset2<-Montreal.sunriset %>% rename(prev.night=Dates,prev.night.sunset.UTC=sunset.UTC) %>% dplyr::select(prev.night,prev.night.sunset.UTC)
  
  filedata4<-left_join(filedata3,Montreal.sunriset2,by="prev.night")
  rm(filedata3)
  
  ####Getting difference between TS and sunset
  filedata4$sunsetdiff<-ifelse(filedata4$time_UTC < filedata4$sunrise.UTC,
                                     as.numeric(difftime(filedata4$time_UTC,filedata4$prev.night.sunset.UTC, units="hours")),
                                     as.numeric(difftime(filedata4$time_UTC,filedata4$sunset.UTC, units="hours")))
  
  batchenddat <- filedata4 %>%
    filter(sunsetdiff> -sunsetbuffer &  #drops observations earlier than 3hrs before sunrise
             sunsetdiff < sunsetbuffer) %>% #drops observations later than 3hrs after sunrise 
    as.data.frame()
  
  rm(filedata4)
  
  
  #verifying
  #  range(batchenddat$sunsetdiff) #FOR TESTING - GOOD - has +3 to -3 range
  
  ###Loop to determine end of activity#####
  ###list of birds with good signals
  
  bird.list <- sort(unique(batchenddat$motusTagID))
  
  batchenddat <- batchenddat %>% arrange(motusTagID, ts)   # Inter-digitate diff.abs values from all ports and stations for one sequence.
  
  crit.nobs <- 800         # number of observations required in focal period to do the changepoint analysis.
  
  dim.j <- max(batchenddat$dos)         # max number of days to look for changepoints for each bird
  dim.i <- length(bird.list)              # number of unique bird ids
  
  # initialize arrays
  
  cpt.ts1.array <- array(NA,dim=c(dim.i,dim.j))	# hold changepoints
  cpt.tsmax <- array(NA,dim=c(dim.i,dim.j))	  # hold maximum ts for that bird
  cpt.tspre <- array(NA,dim=c(dim.i,dim.j))       # hold ts of changepoint - 1 observation
  cpt.tspost <- array(NA,dim=c(dim.i,dim.j))      # hold ts of changepoint + 1 observation
  cpt.meanpre <- array(NA,dim=c(dim.i,dim.j))     # hold mean of diff.abs before the changepoint
  cpt.meanpost <- array(NA,dim=c(dim.i,dim.j))    # hold mean of diff.abs after the changepoint
  cpt.varpre <- array(NA,dim=c(dim.i,dim.j))      # hold var of diff.abs before the changepoint
  cpt.varpost <- array(NA,dim=c(dim.i,dim.j))     # hold var of diff.abs after the changepoint
  cpt.motusTagID <- array(NA,dim=c(dim.i,dim.j))     	# hold id of bird
  
  # assign names to arrays
  dimnames(cpt.ts1.array) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.tsmax) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.tspre) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.tspost) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.meanpre) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.meanpost) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.varpre) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.varpost) <- list(bird.num=1:dim.i, dos=1:dim.j )
  dimnames(cpt.motusTagID) <- list(bird.num=1:dim.i, dos=1:dim.j )
  
  # i<-5 #for testing
  # j<-12 #for testing
  # rm(i,j) #for testing
  
  for (i in 1:dim.i) {                           
    id <- bird.list[i]
    for (j in 1:dim.j) {
      print(paste("motusTagID=",bird.list[i],"i=",i,"dos=",j))
      temp <- batchenddat %>% filter(motusTagID==bird.list[i] & dos==j)													
      n <- length(temp$diff.abs)
      print(paste("n=",n))																							
      if (n > crit.nobs) {
        out1.pmv <- cpt.meanvar(temp$diff.abs,
                                method="AMOC",
                                penalty="Asymptotic",
                                pen.value=0.001,
                                minseglen=20,
                                param.estimates=T,
                                class=T)
        if (ncpts(out1.pmv) > 0) {         
          cpt1 <- cpts(out1.pmv)
          cpt.ts1 <- temp[cpt1,"ts"]
          cpt.ts1.array[i,j] <- cpt.ts1
          cpt.tsmax[i,j] <- max(temp$ts)
          cpt.tspre[i,j]  <- round(cpt.ts1-temp[cpt1-1,"ts"],2)
          cpt.tspost[i,j] <- round(temp[cpt1+1,"ts"]-cpt.ts1,2)
          cpt.meanpre[i,j] <- param.est(out1.pmv)$mean[1]
          cpt.meanpost[i,j] <- param.est(out1.pmv)$mean[2]
          cpt.varpre[i,j] <- param.est(out1.pmv)$variance[1]
          cpt.varpost[i,j] <- param.est(out1.pmv)$variance[2]
          cpt.motusTagID[i,j]<-bird.list[i]
          print(paste(cpt.ts1))
        } else {
          if (ncpts(out1.pmv) == 0) {
            cpt.ts1 <- max(temp$ts)
            cpt.ts1.array[i,j] <- NA
            print(paste(cpt.ts1,"***************************No change points***************************"))
          }}
      }
    }    # end j
  }          # end i
  
  
  
  
  # create data frames from the arrays
  ts1 <- as.data.frame.table(cpt.ts1.array, responseName = "ts1") 
  tsmax <- as.data.frame.table(cpt.tsmax, responseName = "tsmax") 
  tspre <- as.data.frame.table(cpt.tspre, responseName = "tspre") 
  tspost <- as.data.frame.table(cpt.tspost, responseName = "tspost") 
  meanpre <- as.data.frame.table(cpt.meanpre, responseName = "meanpre")
  meanpost <- as.data.frame.table(cpt.meanpost, responseName = "meanpost")
  varpre <- as.data.frame.table(cpt.varpre, responseName = "varpre")
  varpost <- as.data.frame.table(cpt.varpost, responseName = "varpost")
  motusTagID<- as.data.frame.table(cpt.motusTagID, responseName = "motusTagID")
  
  
  # merge the data frames
  cpts <- merge(ts1,tspre,by=c("bird.num","dos"))
  cpts <- merge(cpts,tsmax,by=c("bird.num","dos"))
  cpts <- merge(cpts,tspost,by=c("bird.num","dos"))
  cpts <- merge(cpts,meanpre,by=c("bird.num","dos"))
  cpts <- merge(cpts,meanpost,by=c("bird.num","dos"))
  cpts <- merge(cpts,varpre,by=c("bird.num","dos"))
  cpts <- merge(cpts,varpost,by=c("bird.num","dos"))
  cpts <- merge(cpts,motusTagID,by=c("bird.num","dos"))
  
  
  cpts <- cpts %>% 
    mutate(cpt = as_datetime(ts1,tz="UTC")) %>%
    mutate(tsmax.dt = as_datetime(tsmax)) %>% 
    filter(ts1 > 0|is.na(ts1)==F) %>%   # keep change points, omit instances with no change points
    dplyr::select(bird.num,motusTagID,dos,tspre,tspost,tsmax.dt,meanpre,meanpost,varpre,varpost,cpt)
  
  
  cpts<-cpts %>% arrange (bird.num,dos)
  head(cpts,500)
  names(cpts)
  dim(cpts)                          
  head(cpts)
  
  cpts$batch<-k
  
  ###Saving changepoint batches####
  saveRDS(cpts,file=paste0("./data/cpts_end_raw_","batch_",k,".rds"))
  write.csv(cpts,file=paste0("./data/cpts_end_raw_","batch_",k,".csv"))
  
  rm(batchenddat)
  gc()
  
  #cpt.ts1.array[i,j]
  elapse.t.end[k,"batch"]<-k
  elapse.t.end[k,"elapsed.t"]<- Sys.time()-start.t
  
}



###Inspecting end changepoint metrics####
end.cpts.all<-rbind(readRDS("./data/cpts_end_raw_batch_1.rds"),
                      readRDS("./data/cpts_end_raw_batch_2.rds"),
                      readRDS("./data/cpts_end_raw_batch_3.rds"),
                      readRDS("./data/cpts_end_raw_batch_4.rds"),
                      readRDS("./data/cpts_end_raw_batch_5.rds"),
                      readRDS("./data/cpts_end_raw_batch_6.rds"),
                      readRDS("./data/cpts_end_raw_batch_7.rds"),
                      readRDS("./data/cpts_end_raw_batch_8.rds"),
                      readRDS("./data/cpts_end_raw_batch_9.rds"),
                      readRDS("./data/cpts_end_raw_batch_10.rds"),
                      readRDS("./data/cpts_end_raw_batch_11.rds"),
                      readRDS("./data/cpts_end_raw_batch_12.rds"))

glimpse(end.cpts.all) #3,655 change points

qplot(na.exclude(end.cpts.all$meanpre),bins=50)  # View distribution of meanpre values - mean should be higher (>1) prior to changepoint for end	
qplot(na.exclude(end.cpts.all$meanpost),bins=50) # View distribution of meanpost values - mean should be lower  after changepoint when birds are mostly stationary
qplot(na.exclude(end.cpts.all$meanpre-end.cpts.all$meanpost),bins=50)  # View distribution of diff in mean (most >1 for end)	#some less than 1
qplot(na.exclude(end.cpts.all$varpre),bins=50)  # View distribution of varpre values - should be generally low, but generally higher than post-changepoint
qplot(na.exclude(end.cpts.all$varpost),bins=50) # View distribution of varpost values - should be low
qplot(na.exclude(end.cpts.all$varpre-end.cpts.all$varpost),bins=50)  # View distribution of diff in var - should be low but generally not negative...some potentially some poor changepoints here


####Filter cpts based on criteria####
crit.mean1 <- 0.25   # critical meanpre value
crit.mean2 <- 1   # critical meanpost value 


#####View questionable changepoints for further investigation####
delete <- end.cpts.all[which(end.cpts.all$meanpre < crit.mean1 |
                               end.cpts.all$meanpost > crit.mean2 |						
                               end.cpts.all$meanpre < end.cpts.all$meanpost | 
                               end.cpts.all$varpre < end.cpts.all$varpost),]
delete %>% arrange(bird.num,dos) %>% select(bird.num,cpt,tsmax.dt,meanpre,meanpost,varpre,varpost)

delete <- which(end.cpts.all$meanpre < crit.mean1 |
                  end.cpts.all$meanpost > crit.mean2 |						
                  end.cpts.all$meanpre < end.cpts.all$meanpost | 
                  end.cpts.all$varpre < end.cpts.all$varpost)

nrow(end.cpts.all)
length(delete)

#1199 bad change points out of a total 3655 determined

#Visualizing changepoint distribution after removals
end.cpts.all.filtered<-end.cpts.all[-delete,]
nrow(end.cpts.all.filtered)

qplot(na.exclude(end.cpts.all.filtered$meanpre),bins=50)  # View distribution of meanpre values - mean should be higher (>1) prior to changepoint for end	
qplot(na.exclude(end.cpts.all.filtered$meanpost),bins=50) # View distribution of meanpost values - mean should be lower  after changepoint when birds are more still
qplot(na.exclude(end.cpts.all.filtered$meanpre-end.cpts.all.filtered$meanpost),bins=50)  # View distribution of diff in mean (most >1 for end)
qplot(na.exclude(end.cpts.all.filtered$varpre),bins=50)  # View distribution of varpre values - should below, but generally higher than post-changepoint
qplot(na.exclude(end.cpts.all.filtered$varpost),bins=50) # View distribution of varpost values - should be low
qplot(na.exclude(end.cpts.all.filtered$varpre-end.cpts.all.filtered$varpost),bins=50)  # View distribution of diff in var - should be low but generally not negative...



###Marking questionable changepoints for later removal#####
end.cpts.all.2 <- end.cpts.all
end.cpts.all.2$postfilter<-"good"
end.cpts.all.2[delete,"postfilter"]<-"bad"
glimpse(end.cpts.all.2)          

end.cpts.all.3 <- end.cpts.all.2 %>%
  mutate(time = hour(cpt)+minute(cpt+(second(cpt)/60))/60) %>% arrange(time)


###Convert changepoints to minutes since sunset####
#need to account for cpts past midnight UTC

glimpse(end.cpts.all.3)
end.cpts.all.3$Dates<-date(end.cpts.all.3$cpt)
end.cpts.all.3$prev.night<-end.cpts.all.3$Dates-1
end.cpts.all.3 %>% dplyr::select(Dates, prev.night)

glimpse(Montreal.sunriset2)

end.cpts.all.4<-left_join(end.cpts.all.3,Montreal.sunriset2,by="prev.night") #adds sunset for previous night
end.cpts.all.5<-left_join(end.cpts.all.4,Montreal.sunriset,by="Dates") #adds sunrise and sunset for the changepoint date


end.cpts.all.5$cpt.ss<-ifelse(end.cpts.all.5$cpt < end.cpts.all.5$sunrise.UTC, #if the end of activity is less than sunrise on a given date (i.e. between 00:00 and sunrise)), then...
                                   as.numeric(difftime(end.cpts.all.5$cpt,end.cpts.all.5$prev.night.sunset.UTC, units="mins")), #calculate the difference between the changepoint and suneset the previous calendar day (changepoint will have occurred prior to sunrise on the next calendar day)
                                   as.numeric(difftime(end.cpts.all.5$cpt,end.cpts.all.5$sunset.UTC, units="mins"))) #calculate the difference between the changepoint and sunset in the current day (changepoint time will be greater than sunrise, but less than 24:00hrs)


###Visualizing distribution of onset times for good points
#good changepoints
end.cpts.all.5.good<-end.cpts.all.5 %>% filter(postfilter=="good")
qplot(end.cpts.all.5.good$cpt.ss,bins=30) #looks good

#compared with all changepoints
qplot(end.cpts.all.5$cpt.ss,bins=30) # many at extreme of 3hr buffer, likely not true changes to behaviour

cpts.end<-end.cpts.all.5
glimpse(cpts.end)


###Adding variables for plotting later####
end.cpts.all.5$postfilter<-as.factor(end.cpts.all.5$postfilter)

end.cpts.all.5$postfiltercol<-"black"
end.cpts.all.5[which(end.cpts.all.5$postfilter=="bad"),"postfiltercol"]<-"red"

end.cpts.all.5$postfilterlntyp<-1
end.cpts.all.5[which(end.cpts.all.5$postfilter=="bad"),"postfilterlntyp"]<-2
glimpse(end.cpts.all.5)


###Saving end CPTS####
saveRDS(end.cpts.all.5,file="end.cpts.882.rds")
write.csv(end.cpts.all.5,file="end.cpts.882.csv")






#Visualizing changepoints####
#Multiple stations per plot
#Need to round Motus data to reduce data load, then plot changepoints
#keep batch framework to prevent memory overload

column.names <- c("batch","elapsed.t")
elapse.t.roundplots <- array(NA,dim=c(12,2))
dimnames(elapse.t.roundplots) <- list(NULL,column.names)
elapse.t.roundplots

rdsfilelist

for (k in 1:length(rdsfilelist)){
  # k<-10  #for testing
  filedata<-readRDS(rdsfilelist[k])
  
  start.t<-Sys.time()
  
##Rounding Motus data to make it easier to visualize####
glimpse(filedata)
  
filedata2 <- filedata %>%
  dplyr::select(time_UTC,
                sig,
                port,
                motusTagID,
                recvDeployName)%>%
  mutate(time_UTC.r=ceiling_date(time_UTC,unit="minute")) %>%
  group_by(motusTagID,recvDeployName,port,time_UTC.r)%>%
  filter(sig==max(sig)) %>%
  ungroup() %>%
# distinct() %>% 
  arrange(motusTagID,recvDeployName,port,time_UTC) %>%
  as.data.frame()
  
# glimpse(filedata)
# glimpse(filedata2)
# 
# unique(filedata$recvDeployName) %>% sort()
# unique(filedata2$recvDeployName)%>% sort()
# 
# unique(filedata$motusTagID)%>% sort()
# unique(filedata2$motusTagID)%>% sort()


###Loop to generate plots####
motustagIDlist<-unique(filedata2$motusTagID)
motustagIDlist
# length(motustagIDlist)
# motustagIDlist<-c(unique(df.proj551.MBO3$motusTagID)[1:5],68698)     #FOR TESTING
#	i<-1                                                                #FOR TESTING

#list of plots
signalplots=list()
gc()

for(i in 1:length(motustagIDlist)){
  
  targetbird<-motustagIDlist[i] #select target bird
  
  bird.loop.data<- filedata2 %>%
    dplyr::filter(motusTagID == targetbird)
  
  loop.cpts.onset<-onset.cpts.all.4 %>%
    dplyr::filter(motusTagID == targetbird)

  loop.cpts.end<-end.cpts.all.5 %>%
    dplyr::filter(motusTagID == targetbird)
  
  sunriset.loop<-Montreal.sunriset %>%
    filter(sunrise.UTC > (min(bird.loop.data$time_UTC)-86400) &  sunset.UTC < (max(bird.loop.data$time_UTC)+86400))
  

  signalplots[[i]]<-ggplot() +
    geom_rect(data=sunriset.loop,
              aes(xmin=sunrise.UTC,
                  xmax=sunset.UTC,
                  ymin=-Inf,
                  ymax=Inf),
              fill="white")+
    ylab("Signal strength (dB)")+
    xlab("Date (UTC)")+
    #xlim(c(min(sunriset.loop$sunrise.UTC),max(sunriset.loop$sunset.UTC)))+
    
    scale_x_datetime(limits=c(min(sunriset.loop$sunrise.UTC),
                              max(sunriset.loop$sunset.UTC)),
                       date_breaks = "1 month", date_minor_breaks = "1 week",
                 date_labels = "%b-%y")+
    
    ggtitle(paste0("MotusTagID: ", targetbird))+
    

    #Plotting data colored by station
    geom_point(data=bird.loop.data,aes(x = time_UTC, y = sig, col = as.factor(recvDeployName)),alpha=0.3,size=0.7)+
    
    #Changepoints color coded by post-filtering status
    geom_vline(xintercept = loop.cpts.onset$cpt,col=loop.cpts.onset$postfiltercol,linewidth=0.5,alpha=1,linetype=loop.cpts.onset$postfilterlntyp)+
    geom_vline(xintercept = loop.cpts.end$cpt,col=loop.cpts.end$postfiltercol,linewidth=0.5,alpha=1,linetype=loop.cpts.end$postfilterlntyp)+
    
    scale_color_discrete(name = "Station") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill="grey80"),
      panel.border = element_rect(colour = "black", fill=NA, size = 0.5),
      plot.margin=unit(c(0.2,0.3,0.2,0.2), "cm"),
      axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain")
    )
  rm(bird.loop.data,sunriset.loop)
  gc()

}


####Saving#### 
nplotfile<-5
for(j in 1:ceiling((length(signalplots)/nplotfile))) {
  sigplotj<-signalplots[(1+((j-1)*nplotfile)):(nplotfile+((j-1)*nplotfile))]
  graphics.off()
  pdf(paste0("./plots/Project882tag_plots_rounded_","batch_",k,"_group",j,".pdf"),
      width=24,
      height=8,
      onefile=TRUE)
  for(l in sigplotj) {replay(l)
  }
  graphics.off()
}


rm(signalplots)
rm(filedata)
gc()


elapse.t.roundplots[k,"batch"]<-k
elapse.t.roundplots[k,"elapsed.t"]<- Sys.time()-start.t


}




