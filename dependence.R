script.name <- "dependence"
source("_vars.R")
source("_plot.R")
library("stringr")

#This code plots the dependence of TDX-V on anticoagulation dose
tdx <- getVAR("tdx")
tdx <- tdx[tdx$BDAY>=2,]

tdx <- tdx[tdx$THER %in% c("Enoxaparin","Nadroparin"),]
tdx <- tdx[!is.na(tdx$V),]

log <- tdx[tdx$THER=="Enoxaparin" & tdx$VALUE %in% c("0.4x1","0.4x2","0.6x2","0.8x2","1.0x2","1.2x2"),]
log$VALUE <- factor(log$VALUE,levels = c("0.4x1","0.4x2","0.6x2","0.8x2","1.0x2","1.2x2"))

s1 <- tapply(log$V, log$VALUE, quantile, probs=0.25)
s1 <- melt(s1)
colnames(s1) <- c("Dose","Q1")
s2 <- tapply(log$V, log$VALUE, median)
s2 <- melt(s2)
colnames(s2) <- c("Dose","Q2")
s3 <- tapply(log$V, log$VALUE, quantile, probs=0.75)
s3 <- melt(s3)
colnames(s3) <- c("Dose","Q3")
s <- merge(s2,s1,by="Dose")
s <- merge(s,s3,by="Dose")
rm("s1","s2","s3")

p <- ggplot(s, aes_string(x="Dose",y="Q2")) +
  theme_prism(base_size = 18) +
  theme(aspect.ratio=6/8) +
  geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=29,ymin=20),fill="grey85",color="grey85") +
  geom_line(linewidth=1.2,group = 1) +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=0.5, linewidth=1.3) +
  geom_point(size=4.5) +
  ylim(c(0,50)) +
  scale_x_discrete(guide = "prism_bracket")
p <- p + labs(x="Enoxaparin dosage (ml)", y="TDX V (um/min)")
print(p)



log <- tdx[tdx$THER=="Nadroparin" & tdx$VALUE %in% c("0.3x1","0.6x1","0.3x2","0.6x2","0.8x2","0.9x2","1.2x2"),]
log$VALUE <- factor(log$VALUE,levels = c("0.3x1","0.6x1","0.3x2","0.6x2","0.8x2","0.9x2","1.2x2"))

s1 <- tapply(log$V, log$VALUE, quantile, probs=0.25)
s1 <- melt(s1)
colnames(s1) <- c("Dose","Q1")
s2 <- tapply(log$V, log$VALUE, median)
s2 <- melt(s2)
colnames(s2) <- c("Dose","Q2")
s3 <- tapply(log$V, log$VALUE, quantile, probs=0.75)
s3 <- melt(s3)
colnames(s3) <- c("Dose","Q3")
s <- merge(s2,s1,by="Dose")
s <- merge(s,s3,by="Dose")
rm("s1","s2","s3")

p <- ggplot(s, aes_string(x="Dose",y="Q2")) +
  theme_prism(base_size = 18) +
  theme(aspect.ratio=6/8) +
  geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=29,ymin=20),fill="grey85",color="grey85") +
  geom_line(linewidth=1.2,group = 1) +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=0.5, linewidth=1.3) +
  geom_point(size=4.5) +
  ylim(c(0,50)) +
  scale_x_discrete(guide = "prism_bracket")
p <- p + labs(x="Nadroparin dosage (ml)", y="TDX V (um/min)")
print(p)


######################################################################################################################

#This code plots the dependence of TDX-V on anticoagulation dose per weight
stat <- getVAR("stat")
stat <- stat[!is.na(stat$W),]
stat <- stat[c("UNNUM","W","UNIT")]

tdx <- getVAR("tdx")
tdx <- tdx[tdx$THER %in% c("Enoxaparin","Nadroparin"),]
tdx <- merge(tdx,stat,by="UNNUM")

tdx$VAL <- as.numeric(str_extract(tdx$VALUE, "[\\d\\.]{3}"))
tdx$M <- as.numeric(str_remove(str_extract(tdx$VALUE, "x\\d"),"x"))
tdx <- tdx[tdx$M<=2,]
tdx$TEST <- tdx$VAL*tdx$M/2/tdx$W*10000
tdx <- tdx[!is.na(tdx$TEST),]
tdx$TEST[tdx$THER=="Nadroparin"] <- tdx$TEST[tdx$THER=="Nadroparin"]*9500/10000

tu <- tdx[tdx$UNIT=="TU",]
icu <- tdx[tdx$UNIT=="ICU",]
icu <- icu[sample(nrow(icu), nrow(tu)),]

a <- rbind(tu,icu)
a$THER <- factor(a$THER, levels = c("Enoxaparin","Nadroparin"))
a <- a[!is.na(a$V) & !is.na(a$TEST),]
a$V[a$V<=3] <- 3
model <- lm(1/V ~ TEST, data=a)
summary(model)
a <- a[sample(nrow(a), 1000),]
a$VT <- 1/(3.833e-02 + 4.404e-04*a$TEST)

q <- ggplot(a, aes(x=TEST*2, y=V)) + 
  theme_prism(base_size = 18) + 
  geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=29,ymin=20),fill="grey85",color="grey85") +
  scale_x_continuous(limits = c(10,300)) + 
  geom_point(aes(fill=THER),size=1.5,color="black",pch=21) +
  geom_line(aes(y=VT), color="black", size=1.2)
print(q)


######################################################################################################################

#This code plots the dependence of lab assays on anticoagulation dose

doses <- na.omit(getVAR("tdx")[c("UNNUM","TDAY","THER","VALUE")])
doses <- doses[doses$THER %in% c("Enoxaparin","Nadroparin"),]

test <- "FG"
lab <- getTEST(test)

doses[[test]] <- NA
for(i in seq_len(nrow(doses))){
  unnum <- doses[[i,"UNNUM"]]
  day <- doses[[i,"TDAY"]]
  x <- lab[[test]][lab$DAY>=day+(6/24) & lab$DAY<=day+(30/24) & lab$UNNUM==unnum]
  if(length(x)>0){
    doses[i, test] <- x[[1]]
  }
}
rm("unnum","day","x","i")
doses <- na.omit(doses)

log <- doses[doses$THER=="Enoxaparin" & doses$VALUE %in% c("0.4x1","0.4x2","0.6x2","0.8x2","1.0x2","1.2x2"),]
log$VALUE <- factor(log$VALUE,levels = c("0.4x1","0.4x2","0.6x2","0.8x2","1.0x2","1.2x2"))

#log <- doses[doses$THER=="Nadroparin" & doses$VALUE %in% c("0.3x1","0.6x1","0.3x2","0.6x2","0.8x2","0.9x2","1.2x2"),]
#log$VALUE <- factor(log$VALUE,levels = c("0.3x1","0.6x1","0.3x2","0.6x2","0.8x2","0.9x2","1.2x2"))

s1 <- tapply(log[[test]], log$VALUE, quantile, probs=0.25)
s1 <- melt(s1)
colnames(s1) <- c("Dose","Q1")
s2 <- tapply(log[[test]], log$VALUE, median)
s2 <- melt(s2)
colnames(s2) <- c("Dose","Q2")
s3 <- tapply(log[[test]], log$VALUE, quantile, probs=0.75)
s3 <- melt(s3)
colnames(s3) <- c("Dose","Q3")
s <- merge(s2,s1,by="Dose")
s <- merge(s,s3,by="Dose")
rm("s1","s2","s3")

p <- ggplot(s, aes_string(x="Dose",y="Q2")) +
  theme_prism(base_size = 18) +
  theme(aspect.ratio=6/8) +
  geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=2,ymin=4),fill="grey85",color="grey85") +
  geom_line(linewidth=1.2,group = 1) +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=0.5, linewidth=1.3) +
  geom_point(size=4.5) +
  ylim(c(0,8)) +
  scale_x_discrete(guide = "prism_bracket")
p <- p + labs(x="dosage (ml)", y=test)
if(test=="DD"){
  p <- p + scale_y_log10(limits = c(100,10000))
}
print(p)
