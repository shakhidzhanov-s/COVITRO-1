readVAR <- function(arg){
  return(readRDS(paste0("_vars/F.",arg,".RDS")))
}

varGenerator <- function(){
  stat <- readVAR("STAT")
  diag <- readVAR("DIAG")
  event <- readVAR("EVENT")
  tdx <- readVAR("TDX")
  coag <- readVAR("COAG")
  blood <- readVAR("BLOOD")
  inflam <- readVAR("INFLAM")
  kidney <- readVAR("KIDNEY")
  compl <- readVAR("COMPL")
  biochem <- readVAR("BIOCHEM")
  
  function(var){
    switch (var,
      "stat" = stat,
      "diag" = diag,
      "event" = event,
      "tdx" = tdx,
      "coag" = coag,
      "blood" = blood,
      "inflam" = inflam,
      "kidney" = kidney,
      "compl" = compl,
      "biochem" = biochem,
      NULL
    )
  }
}

if(!exists("getVAR")){
  getVAR <- varGenerator()
}

searchGenerator <- function(){
  coag.cols <- c("V","VI","D","TSP","APTT","PT","PTR","INR","DD","FG","TT","AT3","AXA","R","K","ADEG","MA")
  blood.cols <- c("RBC","PLT","WBC","NEUT","LYMP")
  inflam.cols <- c("CRP","FER")
  kidney.cols <- c("CRE","CKD.EPI")
  biochem.cols <- c("LDHh","LDHl")
  
  cols <- c(coag.cols, blood.cols, inflam.cols, biochem.cols, kidney.cols)
  vs <- c(rep("coag", length(coag.cols)),
          rep("blood", length(blood.cols)),
          rep("inflam", length(inflam.cols)),
          rep("biochem", length(biochem.cols)),
          rep("kidney", length(kidney.cols))
          )
  
  function(test){
    if(!Reduce(`&`, test %in% cols)){
      return(NULL)
    }
    VR <- data.frame()
    for(i in 1:length(test)){
      vr <- getVAR(vs[which(cols==test[i])])
      vr <- vr[!is.na(vr[[test[i]]]),c("UNNUM","DAY",test[i])]
      if(i==1){
        VR <- vr
      } else if(i>1) {
        VR <- merge(VR, vr, by=c("UNNUM","DAY"), all=TRUE)
      }
    }
    return(VR)
  }
}

if(!exists("getTEST")){
  getTEST <- searchGenerator()
}

rm("readVAR", "varGenerator", "searchGenerator")
