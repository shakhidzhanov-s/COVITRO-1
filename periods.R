#This code evaluates the association of hypercoagulability with thromboembolism
source("_vars.R")
library("stringr")
library("ggplot2")
library("pROC")

#################################################################################

toFrame <-function(lst){
  #Special function
  #The function a frame from a list
  lst <- lst[!sapply(lst,is.null)]
  col_names <- colnames(lst[[1]])
  lst <- lapply(lst, function(x){x[col_names]})
  df <- do.call(rbind,lst)
  rownames(df) <- NULL
  return(df)
}
toPatList <- function(df){
  #Special function
  #The function makes a list of patients from a frame
  pnames <- unique(df$UNNUM)
  pnames <- pnames[!is.na(pnames)]
  ldf <- lapply(pnames, function(x) {df[df$UNNUM==x & !is.na(df$UNNUM),]})
  names(ldf) <- pnames
  return(ldf)
}
getThrombPeriods <- function(tillFirstEvent=T){
  #The function divides patient's time course into periods between ultrasound imagings.
  #The periods could end up with thrombosis or not.
  #Thromboses which could not be found using lower or upper limbs ultrasound imagings are considered
  # to start N days before thrombosis discovery.
  #After that lab values can be evaluated according to the period type.
  
  #Predefined constants
  PredefUncondPeriodLength <- 5 #days
  
  #Simple preparations.
  event <- getVAR("event")
  thrs <- event[event$TYPE %in% c("AC.ISH","AR.THR","DV.THR","NO.THR","PE.THR","SV.THR"),]
  thrs$TYPE <- str_remove(thrs$TYPE, "\\.THR|\\.ISH")
  rm("event")
  thrs <- thrs[order(thrs$UNNUM,thrs$DAY),]
  
  #exlude data which were acquired after an event (count first event only)
  if(tillFirstEvent){
    lthrs <- toPatList(thrs)
    lthrs <- lapply(lthrs, function(x){
      s <- nrow(x)
      if(any(x$TYPE!="NO")){
        s <- which(x$TYPE!="NO")[[1]]
      }
      return(x[1:s,])
    })
    thrs <- toFrame(lthrs)
    rm("lthrs")
  }

  #PE is counted together with lower extremity VT.
  thrs$LOCALIZ[thrs$TYPE=="PE"] <- "L.LIMBS"
  #DVT and sural vein thrombosis are considered together.
  thrs$TYPE[thrs$TYPE %in% c("DV","SV")] <- "VT"
  #Lower and upper limbs are considered separately.
  l.thr.l <- toPatList(thrs[thrs$LOCALIZ=="L.LIMBS",])
  u.thr.l <- toPatList(thrs[thrs$LOCALIZ=="U.LIMBS",])
  #The other thromboses are considered separately. These thrombosis types will not depend on
  # previously obtained instrumental result (if any exits) but will have predefined maximal length.
  #CORONAR events are excluded from the analysis, CHEST is mainly PE.
  e.thr.l <- toPatList(thrs[thrs$LOCALIZ %in% c("ABDOMEN", "BRAIN", "ELSE") & thrs$TYPE!="NO",])
  rm("thrs")
  
  periodsFromUS <- function(le){
    #The function creates ultrasound periods from two patient's sequential ultrasound
    # imaging results (if it is possible or it returns null)
    if(nrow(le)>1){
      le1 <- le[1:(nrow(le)-1),]
      le2 <- le[2:nrow(le),]
      unnum <- le1$UNNUM
      day1 <- as.numeric(le1$DAY)
      type1 <- le1$TYPE
      day2 <- as.numeric(le2$DAY)
      type2 <- le2$TYPE
      type <- paste0(type1,":",type2)
      durat <- day2 - day1
      loc <- le1$LOCALIZ
      return(data.frame("UNNUM"=unnum,"FDAY"=day1,"LDAY"=day2,"DURAT"=durat,"TYPE"=type,"LOCALIZ"=loc))
    }
    else{
      return(NULL)
    }
  }
  
  periodsFromPredefLength <- function(le){
    #The functions creates periods using predefined maximal period length; the periods end up
    # with an event but their beginning is dependent on PredefUncondPeriodLength constant only.
    unnum <- le$UNNUM
    day2 <- as.numeric(le$DAY)
    type2 <- le$TYPE
    day1 <- pmax(day2-PredefUncondPeriodLength,0)
    type <- paste0("NA:",type2)
    durat <- day2 - day1
    loc <- le$LOCALIZ
    return(data.frame("UNNUM"=unnum,"FDAY"=day1,"LDAY"=day2,"DURAT"=durat,"TYPE"=type,"LOCALIZ"=loc))
  }
  
  #apply functions to the lists of events to create the period frame.
  l.pers <- Filter(Negate(is.null), lapply(l.thr.l, periodsFromUS))
  u.pers <- Filter(Negate(is.null), lapply(u.thr.l, periodsFromUS))
  e.pers <- Filter(Negate(is.null), lapply(e.thr.l, periodsFromPredefLength))
  rm("l.thr.l","u.thr.l","e.thr.l")
  
  #create an empty data frame
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(df) <- c("UNNUM","FDAY","LDAY","DURAT","TYPE","LOCALIZ")
  #append l.pers and u.pers to df
  if(length(l.pers)>0){
    df <- rbind(df,toFrame(l.pers))
  }
  if(length(u.pers)>0){
    df <- rbind(df,toFrame(u.pers))
  }
  #exclude the No:No US periods which contain else_type of thrombosis
  if(length(e.pers)>0 & nrow(df)>0){
    df$SAVE <- T
    for(e in e.pers){
      for(i in 1:nrow(e)){
        df$SAVE[df$TYPE=="NO:NO" & df$UNNUM==e[i,"UNNUM"] & df$FDAY<=e[i,"LDAY"] & df$LDAY>=e[i,"LDAY"]] <- F
      }
    }
    df <- df[df$SAVE,]
    df$SAVE <- NULL
    rm("e","i")
  }
  #append e.pers to df
  if(length(e.pers)>0){
    df <- rbind(df,toFrame(e.pers))
  }
  
  rm("e.pers","l.pers","u.pers")
  
  #numerate periods
  ldf <- toPatList(df)
  ldf <- lapply(ldf, function(l){
    l$NP <- paste0(l$UNNUM,":",1:nrow(l))
    return(l)
  })
  df <- toFrame(ldf)
  rm("ldf")
  return(df)
}
filterPeriods <- function(periods, test, quality=0.01){
  #This function filter bad periods and calculates quality of measurements within each period
  periods <- periods[periods$TYPE %in% c("NO:NO","NO:VT","NO:PE","NA:VT","NA:AC","NA:AR","NO:AR"),]
  periods <- periods[periods$DURAT>=2 & periods$DURAT<=8,]
  periods <- periods[!(periods$TYPE %in% c("NA:VT","NA:AC","NA:AR")) | periods$FDAY>=3,]
  
  periods[paste0("Q_",test)] <- 0
  for(j in test){
    j.var <- getTEST(j)
    for(i in 1:nrow(periods)){
      un <- as.character(periods[i,"UNNUM"])
      fd <- as.numeric(periods[i,"FDAY"])
      ld <- as.numeric(periods[i,"LDAY"])
      dr <- as.numeric(periods[i,"DURAT"])
      n <- nrow(j.var[j.var$UNNUM==un & j.var$DAY>=fd & j.var$DAY<ld,])
      periods[i,paste0("Q_",j)] <- n/dr
    }
  }
  
  if(quality>0){
    for(j in test){
      periods <- periods[periods[paste0("Q_",j)]>=quality,]
    } 
  }
  
  return(periods)
}
preparePeriods <- function(periods, test){
  
  collectDataWithinPeriods <- function(periods, test){
    #collect lad data within periods
    if(nrow(periods)<1){
      print("Empty periods!")
      return(NULL)
    }
    
    #Simple preparations
    dat <- getTEST(test)
    
    #create empty data frame
    df <- data.frame(matrix(ncol = 5+length(test), nrow = 0))
    colnames(df) <- c("UNNUM","DAY",test,"TYPE","NP","DURAT")
    
    #collect lab data in the periods
    for(i in 1:nrow(periods)){
      un <- periods[i,"UNNUM"]
      fd <- periods[i,"FDAY"]
      ld <- periods[i,"LDAY"]
      tp <- periods[i,"TYPE"]
      np <- periods[i,"NP"]
      dur <- periods[i,"DURAT"]
      dt <- dat[dat$UNNUM==un & dat$DAY>=fd & dat$DAY<=ld,]
      if(nrow(dt)>0){
        dt$TYPE <- tp
        dt$NP <- np
        dt$DURAT <- dur
        df <- rbind(df,dt)
      }
    }
    rm("i","dt","un","fd","ld","tp","np")
    return(df[c("UNNUM","DAY","NP","TYPE","DURAT",test)])
  }
  
  dat <- collectDataWithinPeriods(periods, test)
  
  periods[c("MAX","MIN")] <- NA
  for(i in 1:nrow(periods)){
    s <- dat[dat$NP==as.character(periods[i,"NP"]),test]
    periods$MAX[i] <- max(s)
    periods$MIN[i] <- min(s)
  }
  rm("s","i","dat")
  
  periods <- periods[periods$TYPE %in% c("NA:AC","NA:AR","NO:AR","NO:NO","NO:PE","NO:VT"),]
  periods$TYPE <- str_remove(periods$TYPE, "NO:|NA:")
  periods$TYPE[periods$TYPE %in% c("AC","AR")] <- "AT"
  periods$MIDP <- (periods$FDAY + periods$LDAY)/2
  periods <- periods[c("UNNUM","MIDP","TYPE",paste0("Q_",test),"MAX","MIN","NP")]
  
  #TDX-V is a fast changing parameter, here we leave periods which has indicated
  #hypercoagulability of we sure that we made enough measurements to be sure that
  #we did not overlook it.
  if(test=="V"){
    periods <- periods[periods$MAX>=50 | periods$Q_V>=0.7,]
  }
  periods$THR <- F
  periods$THR[periods$TYPE!="NO"] <- T
  
  periods$TYPE <- factor(periods$TYPE, levels = c("NO","AT","PE","VT"))
  periods$THR <- factor(periods$THR, levels = c(T,F))
  
  return(periods)
}

test <- c("V")
pers <- getThrombPeriods()
pers <- filterPeriods(pers, test)
pers <- preparePeriods(pers, test)

p <- ggplot(data=pers, aes(x=TYPE, y=MAX)) + geom_boxplot()
print(p)

myroc <- roc(pers$THR, pers$MAX, plot=TRUE, print.auc=TRUE, show.thres=TRUE)

mycoords <- coords(myroc, "all")
plot(as.numeric(mycoords[["threshold"]]), as.numeric(mycoords[["specificity"]]), type="l", col="red")
lines(as.numeric(mycoords[["threshold"]]), as.numeric(mycoords[["sensitivity"]]), type="l", col="blue")
