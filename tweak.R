#This code evaluate trajectories of ICU patients with and without
#in-hospital thrombosis, it logic is similar to the code in dynamics.R
script.name <- "tweak"

source("_vars.R")
source("_plot.R")

collectLabData <- function(labTest){
  lab <- getTEST(labTest)
  per.low <- seq(0,18,by=3)-1
  per.upp <- seq(0,18,by=3)+1
  per.mid <- seq(0,18,by=3)
  per.tit <- paste0(per.low, "-", per.upp)
  per.tit[[1]] <- "0-1"
  
  lab$DAY.P <- NA
  lab$MID <- NA
  for(i in seq_along(per.mid)){
    lab$DAY.P[lab$DAY>=per.low[[i]] & lab$DAY<per.upp[[i]]+1] <- per.tit[[i]]
    lab$MID[lab$DAY>=per.low[[i]] & lab$DAY<per.upp[[i]]+1] <- per.mid[[i]]
  }
  lab <- lab[!is.na(lab$MID),]
  lab$DIST <- abs(lab$DAY-lab$MID)
  lab$ID <- paste0(lab$UNNUM, ":", lab$MID)
  lab <- lab[order(lab$UNNUM, lab$MID, lab$DIST),]
  lab <- lab[!duplicated(lab$ID),]
  
  lab$DAY <- lab$DAY.P
  lab$DAY <- factor(lab$DAY, levels = per.tit)
  
  rm("i","per.low","per.upp","per.mid")
  return(lab[c("UNNUM","DAY","MID",labTest)])
}

thrombPats <- function(){
  events <- getVAR("event")
  events <- events[events$TYPE %in% c("AR.THR","DV.THR","PE.THR","SV.THR"),]
  events <- events[order(events$UNNUM, events$DAY),]
  events <- events[!duplicated(events$UNNUM),]
  late <- events$UNNUM[events$DAY>=3]
  early <- events$UNNUM[events$DAY<3]
  return(list(unique(late),unique(early)))
}

filter <- function(labTest, df){
  stat <- getVAR("stat")
  stat <- stat[stat$PRIME,]
  stat <- stat[stat$UNIT=="ICU",]
  
  if(labTest %in% c("V","VI","D","R","K","MA","ADEG")){
    df <- df[(df$UNNUM %in% stat$UNNUM[stat$PH < 200]) | df$DAY!="0-1",]
  }
  if(labTest=="TSP"){
    df <- df[df$TSP<=30,]
  }
  if(labTest %in% c("APTT","PT")){
    stat <- stat[stat$PH > 104,]
  }
  df <- df[df$UNNUM %in% stat$UNNUM,]
  
  rm("stat")
  return(df)
}

test <- "AT3"

df <- filter(test, collectLabData(test))
tl <- thrombPats()
df$TAR <- FALSE
df <- df[!(df$UNNUM %in% tl[[2]]),]
df$TAR[df$UNNUM %in% tl[[1]]] <- TRUE

if(!(test %in% c("R","K","ADEG","MA","AT3")))
{
  df$SAVE <- T
  ndf <- data.frame(table(df$TAR, df$DAY))
  for(i in unique(df$DAY)){
    n <- ndf$Freq[ndf$Var1==T & ndf$Var2==i]
    nc <- ndf$Freq[ndf$Var1==F & ndf$Var2==i]
    if(nc>n){
      cg <- which(df$TAR==F & df$DAY==i)
      scg <- sample(cg,nc-n)
      df$SAVE[scg] <- F
      df <- df[df$SAVE,]
    } else if(n>nc){
      cg <- which(df$TAR==T & df$DAY==i)
      scg <- sample(cg,n-nc)
      df$SAVE[scg] <- F
      df <- df[df$SAVE,]
    } else {
      df <- df
    }
  }
  df$SAVE <- NULL
}

params <- list("IND.VAR"="MID", 
               "LAB.TEST"=test, 
               "FILL"="TAR",
               "IND.VAR.LABEL"="Stay day")
plotLabBoxData(df, params)
