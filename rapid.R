#This code shows more closely features of TDX V rapid decline after admission to hospital.
script.name <- "rapid"

library("ggplot2")
source("_vars.R")

#################################################################################

#This function collects TDX V data on 0-1 and 2-4 hospital stay day when the dramatic
# changes can be seen.
#Here the nearest to the 0th and 2nd stay day is chosen if a patient had more than
# 1 measurement in a period.
#The patients from the early pandemic stages are chosen because they administered
# their first anticoagulant doses on admission to hospital (they were not given in
# an ambulance or at home).
#The patients with known AC doses are included.
collectTDXVData <- function(){
  #COVID19 patients hospitalized from home from the first pandemic stages are chosen.
  #Load relevant data.
  stat <- getVAR("stat")
  stat <- stat[stat$PRIME,]
  stat <- stat[stat$PH<200,]
  #IDs of patients' subset.
  selpats <- stat$UNNUM
  rm("stat")
  
  #Load TDX and anticoagulation data.
  tdx <- getVAR("tdx")
  #Leave only relevant information.
  tdx <- tdx[tdx$UNNUM %in% selpats,c("UNNUM","BDAY","V","AC","THER","VALUE")]
  #Leave measurements performed up to the forth hospital stay day.
  tdx <- tdx[tdx$BDAY<5,]
  
  #Separate the tdx data frame into two - with the first measurement up to the first
  # hospital day and with the measurement performed after the second hospital stay day.
  #This date will will processed separately and merged together after.
  frst <- tdx[tdx$BDAY<1,]
  scnd <- tdx[tdx$BDAY>=2,]
  rm("tdx")
  
  #Leave only the patients who did not administer AC before the first TDV V measurement
  frst <- frst[frst$THER %in% c("NoTherapy","Unknown"), c("UNNUM","BDAY","V")]
  frst <- frst[order(frst$UNNUM, frst$BDAY),]
  frst <- frst[!duplicated(frst$UNNUM),]
  
  #Leave patients with known AC therapy
  scnd <- scnd[scnd$THER %in% c("Enoxaparin", "Nadroparin"), ]
  scnd <- scnd[!is.na(scnd$VALUE),]
  scnd <- scnd[order(scnd$UNNUM, scnd$BDAY),]
  scnd <- scnd[!duplicated(scnd$UNNUM),]

  #the first and the second data frames have to have the same set of patients
  set <- intersect(frst$UNNUM, scnd$UNNUM)
  frst <- frst[frst$UNNUM %in% set,]
  scnd <- scnd[scnd$UNNUM %in% set,]
  
  #Merge the data
  df <- merge(frst, scnd, by="UNNUM")
  df$DIFF <- df$V.y - df$V.x
  
  rm("frst","scnd")
  df$AC <- factor(df$AC, levels=c("Proph.","Inter.","Ther."))
  return(df)
}

#################################################################################

df <- collectTDXVData()

w <- ggplot(df, aes(x=V.x, y=DIFF)) + 
  geom_vline(xintercept = c(20,29), size=1, color="black") + 
  geom_hline(yintercept = c(0), size=1, color="black") +
  labs(x="TDX V on day 0, um/min", y="TDX V change (um/min)") + 
  geom_point(aes(fill=AC),size=1.5,color="black",pch=21) + 
  geom_smooth(method='lm')
print(w)

df$V.s <- (df$V.x - mean(df$V.x))/sd(df$V.x)
fit <- lm(df$DIFF~df$V.s)
summary(fit)
fit <- lm(df$DIFF~df$V.s+df$AC)
summary(fit)
