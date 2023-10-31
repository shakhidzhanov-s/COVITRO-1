library("ggplot2")
library("reshape2")
library("ggprism")

plot.nrange <- function(labTest){
  switch(labTest,
         "V" = c(20,29),
         "VI" = c(38,56),
         "D" = c(15000,32000),
         "PLT" = c(180,380),
         "FG" = c(2,4),
         "INR" = c(0.8,1.2),
         "PT" = c(9.4,12.5),
         "APTT" = c(24,36),
         "DD" = c(-Inf,500),
         "TT" = c(10.3,16.6),
         "R" = c(9,27),
         "K" = c(2,9),
         "ADEG" = c(22,58),
         "MA" = c(44,64),
         "CRP" = c(-Inf,5),
         "WBC" = c(4,9),
         "RBC" = c(4.1,5.1),
         "FER" = c(11,336),
         "LYMP" = c(0.8,3.3),
         "NEUT" = c(1.5,6.1),
         "PTR" = c(70,130),
         "AT3" = c(66,124),
         "CKD.EPI" = c(60,Inf),
         "LDHl" = c(-Inf,250),
         NULL
  )
}

plot.limits <- function(labTest){
  switch(labTest,
         "V" = c(0,100),
         "VI" = c(0,100),
         "TSP" = c(0,30),
         "FG" = c(0,15),
         "PLT" = c(0,800),
         "DD" = c(50,150000),
         "R" = c(0,40),
         "PT" = c(0,30),
         "PTR" = c(0,150),
         "APTT" = c(0,100),
         "INR" = c(0,3),
         "D" = c(0,50000),
         "K" = c(0,20),
         "ADEG" = c(0,80),
         "MA" = c(0,100),
         "TT" = c(0,60),
         "RBC" = c(0,8),
         "NEUT" = c(0,30),
         "CRP" = c(0,400),
         "FER" = c(0,3000),
         "LYMP" = c(0,10),
         "WBC" = c(0,40),
         NULL
         )
}

plot.labels <- function(labTest){
  switch(labTest,
         "V" = "TDX V, um/min",
         "VI" = "TDX VI, um/min",
         "D" = "TDX D, AU",
         "PLT" = "PLT, 10^9/l",
         "FG" = "Fibrinogen, g/l",
         "INR" = "INR",
         "PT" = "Prothrombin time, sec",
         "PTR" = "Prothrombin index, %",
         "APTT" = "aPTT, sec",
         "DD" = "D-dimer, ng/ml",
         "TSP" = "TDX TSP, min",
         "TT" = "Thrombin time, sec",
         "R" = "TEG R, min",
         "K" = "TEG K, min",
         "ADEG" = "TEG a, deg",
         "MA" = "TEG MA, mm",
         "RBC" = "RBC, 10^12/l",
         "NEUT" = "NEUT, 10^9/l",
         "CRP" = "CRP, mg/l",
         "WBC" = "WBC, 10^9/l",
         "LYMP" = "LYMPH, 10^9/l",
         "FER" = "Ferritin, ng/ml",
         "AT3" = "Antithrombin III, %",
         NULL
         )
}

plotLabBoxData <- function(dat, params){
  if(!exists("script.name")){
    stop("Variable 'script.name' does not exist!")
  }
  
  ind.var.dat <- params[["IND.VAR"]]
  lab.dat <- params[["LAB.TEST"]]
  fill.dat <- params[["FILL"]]
  ind.val.lbl <- params[["IND.VAR.LABEL"]]
  
  if(is.null(ind.var.dat) | is.null(lab.dat)){
    stop("No necessary data is given")
  }
  
  lowest.dat <- quantile(dat[[lab.dat]], probs = 0.01)[[1]]
  highest.dat <- quantile(dat[[lab.dat]], probs = 0.99)[[1]]
  
  if(!is.null(plot.limits(lab.dat))){
    plot.lims <- plot.limits(lab.dat)
  } else {
    plot.lims <- c(lowest.dat, highest.dat)
  }
  
  plot.lbl <- plot.labels(lab.dat)
  if(is.null(plot.lbl)){
    plot.lbl <- lab.dat
  }
  
  if(!is.null(fill.dat)){
    tab <- as.data.frame(table(dat[[fill.dat]],dat[[ind.var.dat]]))
    ns.x <- sapply(unique(dat[[ind.var.dat]]), function(ind.v){
      freqs <- tab$Freq[tab$Var2==ind.v]
      n.x <- paste0(ind.v,"\n","(n=",paste(freqs, collapse="\n"), ")")
      return(n.x)
    })
  }
  
  if(ind.var.dat=="MID"){
    s1 <- tapply(dat[[lab.dat]], list(dat[[ind.var.dat]], dat[[fill.dat]]), quantile, probs=0.25)
    s1 <- melt(s1)
    colnames(s1) <- c(ind.var.dat,fill.dat,"Q1")
    s2 <- tapply(dat[[lab.dat]], list(dat[[ind.var.dat]], dat[[fill.dat]]), median)
    s2 <- melt(s2)
    colnames(s2) <- c(ind.var.dat,fill.dat,"Q2")
    s3 <- tapply(dat[[lab.dat]], list(dat[[ind.var.dat]], dat[[fill.dat]]), quantile, probs=0.75)
    s3 <- melt(s3)
    colnames(s3) <- c(ind.var.dat,fill.dat,"Q3")
    s <- merge(s2,s1,by=c(ind.var.dat,fill.dat))
    s <- merge(s,s3,by=c(ind.var.dat,fill.dat))
    rm("s1","s2","s3")
    
    s[[ind.var.dat]] <- s[[ind.var.dat]] + (as.numeric(s[[fill.dat]]) - mean(as.numeric(s[[fill.dat]])))/4
    
    p <- ggplot(s, aes_string(x=ind.var.dat, color=fill.dat)) +
      theme_prism(base_size = 18) +
      theme(aspect.ratio=6/8) +
      geom_rect(aes(xmax=Inf,xmin=-Inf,ymax=plot.nrange(lab.dat)[2],ymin=plot.nrange(lab.dat)[1]),fill="grey85",color="grey85") +
      geom_hline(yintercept = plot.nrange(lab.dat), size=1, color="black") +
      geom_line(aes(y=Q2),linewidth=1.2) +
      geom_errorbar(aes(ymin=Q1, ymax=Q3), width=0.5, linewidth=1.3) +
      geom_point(aes(y=Q2), size=4.5) +
      labs(x=ind.val.lbl, y=plot.lbl) +
      expand_limits(x = 0, y = 0) +
      scale_colour_manual(values = c("black", "#bb0a21", "#3f88c5"))
    
  } else {
    p <- ggplot(dat, aes_string(x=ind.var.dat, y=lab.dat, fill=fill.dat)) +
      geom_hline(yintercept = plot.nrange(lab.dat), size=1, color="black") +
      stat_boxplot() +
      labs(x=ind.val.lbl, y=plot.lbl) +
      coord_cartesian(ylim = plot.lims)
  }
  
  if(test=="DD"){
    #p <- p + scale_y_log10()
  }
  print(p)
}