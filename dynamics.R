#This code can be used to assess laboratory test time-courses.
script.name <- "dynamics"

source("_vars.R")
source("_plot.R")

#################################################################################

#This function collects laboratory test time-course data and attributes it to
#a specific time-period.
#A desired lab test name have to be given (only one).
collectLabData <- function(labTest){
  #Load the desired laboratory data.
  if(length(labTest)!=1){
    print(paste0("The labTest length is not 1: ", length(labTest)))
    return(NULL)
  }
  lab <- getTEST(labTest)
  if(is.null(lab)){
    print(paste0("No such laboratory data:", labTest))
    return(NULL)
  }
  
  #The data will be attributed to a several specific time-periods:
  #From 0 to 1 day, from 2 to 4 day, ... , from 14 to 16 day
  #If a patients has several test measurement within the time period the nearest to its
  #middle will be chosen.
  
  #the lower period limits: -1 (actually 0 because no points measured on -1 day), 2, 5, 8, 11 and 14 day
  per.low <- seq(0,15,by=3)-1
  #the upper period limits: 1, 4, 7, 10, 13 and 16 day
  per.upp <- seq(0,15,by=3)+1
  #the period middle points: 0, 3, 6, 9, 12 and 15 day
  per.mid <- seq(0,15,by=3)
  #period titles: 0-1, 2-4, 5-7, 8-10, 11-13 and 14-16
  per.tit <- paste0(per.low, "-", per.upp)
  per.tit[[1]] <- "0-1"
  
  #here we attribute each point to a specific time-period.
  lab$DAY.P <- NA
  lab$MID <- NA
  for(i in seq_along(per.tit)){
    lab$DAY.P[lab$DAY>=per.low[[i]] & lab$DAY<per.upp[[i]]+1] <- per.tit[[i]]
    lab$MID[lab$DAY>=per.low[[i]] & lab$DAY<per.upp[[i]]+1] <- per.mid[[i]] + 0.5
  }
  #delete all points outside the periods
  lab <- lab[!is.na(lab$DAY.P),]
  
  #Leave the nearest to a period middle points
  #Calculate distances between the points and their middles
  lab$DIST <- abs(lab$DAY-lab$MID)
  #Add unique id to each period
  lab$ID <- paste0(lab$UNNUM, ":", lab$DAY.P)
  #Sort by name, then by mid and then by dist
  #Now the points with lowerest dists will appear upper in the data frame
  lab <- lab[order(lab$UNNUM, lab$MID, lab$DIST),]
  #Delete duplicate IDs
  lab <- lab[!duplicated(lab$ID),]
  
  #Rename the day.p column
  lab$DAY <- lab$DAY.P
  #Factorize day column
  lab$DAY <- factor(lab$DAY, levels = per.tit)
  
  rm("i","per.low","per.upp","per.mid","per.tit")
  return(lab[c("UNNUM","DAY","MID",labTest)])
}

#This function evaluates if a patient was ever in ICU or not and a patients outcome
collectStatData <- function(){
  stat <- getVAR("stat")[c("UNNUM","UNIT","OUT")]
  return(stat)
}

#This function prepares the final data.
#It deletes non-covid patients and the patients who received anti-covid treatment
#in another hospital prior to hospitalization or had no COVID19 infection.
#It also does some test-specific preparations
#A list of dataframes to merge have to be given.
prepareData <- function(labTest, lst){
  #Preliminary checks
  if(length(lst)==0){
    print("A list have to be given.")
    return(NULL)
  }
  if(length(lst)==1 & !is.data.frame(lst) & is.list(lst)){
    df <- lst[[1]]
  }
  if(!Reduce(`&`, sapply(lst, is.data.frame))){
    print("A list have to contain data frames only.")
    return(NULL)
  }
  #Merge the data frames
  df <- Reduce(function(x,y,...) merge(x,y,by="UNNUM",all=T), lst)
  
  #exclude patients who already received COVID19 treatment in another hospital
  stat <- getVAR("stat")
  stat <- stat[stat$PRIME==TRUE,]
  
  #Some tests are very sensitive to anticoagulant administration.
  #At the start of the pandemic anticoagulants were given on admission to hospital.
  #However, later it was decided to give anticoagulants immediately in an ambulance.
  #This dramatically changed the picture of TDX V at admission.
  #Thereby we include only the patients who were admitted without anticoagulant
  #administration (the patients from the first phase only)
  #Here it is important to exclude the first day points.
  if(labTest %in% c("V","VI","D","R","K","MA","ADEG")){
    df <- df[(df$UNNUM %in% stat$UNNUM[stat$PH < 200]) | df$DAY!="0-1",]
  }
  #TDX TSP which were detected after 30 minutes of the test start are due to prolonged measurements.
  #Prolonged measurements were done rarely and in one hospital only. The correct test duration is 30 min/
  if(labTest=="TSP"){
    df <- df[df$TSP<=30,]
  }
  #Hospitals < 104 used aPTT kits with another normal ranges. We excluded them.
  if(labTest %in% c("APTT","PT")){
    stat <- stat[stat$PH > 104,]
  }
  
  df <- df[df$UNNUM %in% stat$UNNUM,]
  df <- df[!is.na(df[[labTest]]),]
  
  rm("stat")
  return(df)
}

#################################################################################

test <- "V"
lab <- collectLabData(test)
stat <- collectStatData()
df <- prepareData(test, list(lab,stat))
df$TeT <- paste0(df$UNIT,":",df$OUT)
df$TeT[df$TeT=="TU:DEA"] <- "ICU:DEA"
df$TeT <- factor(df$TeT, levels=c("ICU:DEA","ICU:REC","TU:REC"))
rm("lab","stat")

params <- list("IND.VAR"="MID", 
               "LAB.TEST"=test, 
               "FILL"="TeT", 
               "IND.VAR.LABEL"="Stay day")
plotLabBoxData(df, params)
#Width - 6.85, height - 4.55 inches


un <- unique(df$UNNUM[df$DAY=="2-4" & df$V>=40])
df <- df[df$UNNUM %in% un,]
stat <- getVAR("stat")
df <- df[(df$UNNUM %in% stat$UNNUM[stat$PH<200]),]
rm("stat")

s1 <- tapply(df[[test]], df$MID, quantile, probs=0.25)
s1 <- melt(s1)
colnames(s1) <- c("MID","Q1")
s2 <- tapply(df[[test]], df$MID, median)
s2 <- melt(s2)
colnames(s2) <- c("MID","Q2")
s3 <- tapply(df[[test]], df$MID, quantile, probs=0.75)
s3 <- melt(s3)
colnames(s3) <- c("MID","Q3")
s <- merge(s2,s1,by="MID")
s <- merge(s,s3,by="MID")
rm("s1","s2","s3")


p <- ggplot(s, aes(x=MID)) +
  theme_prism(base_size = 18) +
  theme(aspect.ratio=6/8) +
  geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=29,ymin=20),fill="grey85",color="grey85") +
  geom_line(aes(y=Q2),linewidth=1.2) +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=0.5, linewidth=1.3) +
  geom_point(aes(y=Q2), size=4.5) +
  labs(x="Stay day", y="TDX-V, um/min") +
  expand_limits(x = 0, y = 0) +
  scale_colour_manual(values = c("black"))

print(p)
