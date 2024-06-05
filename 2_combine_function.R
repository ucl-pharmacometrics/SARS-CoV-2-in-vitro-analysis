###################FORMAT THE DATA for the drug combinations#####################
extract.2combtherapy <- function(df, drugno_A, drugno_B, ConcUnit){
  # mergedcsv <- "Wuhan_merged.csv"
  # drugno_A <- 3
  # drugno_B <- 2
  platedf <- df
  # work out which of the 4 drugs the user has NOT chosen
  exclude <- c(1:4)
  exclude <- exclude[exclude != drugno_A]
  exclude <- exclude[exclude != drugno_B]
  # remove wells with other drugs
  platedf <- platedf[platedf[paste0("drug", exclude[1], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[2], "_conc")] == 0, ]
  # select nonzero drug wells
  platedf$conc1 <- platedf[paste0("drug", drugno_A, "_conc")]
  platedf$conc2 <- platedf[paste0("drug", drugno_B, "_conc")]
  analysisdf <- platedf[platedf$conc1 != 0, ]
  analysisdf <- analysisdf[analysisdf$conc2 != 0, ]
  # include the data with concentration of 0
  for(i in unique(analysisdf$conc1)){
    for(j in unique(analysisdf$conc2)){
      analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == i & platedf$conc2 == 0,])
      analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == 0 & platedf$conc2 == j,])
    }
  }
  #  # include the response of 0 when no drug is added
  df0 <- platedf[platedf$conc1 == 0 & platedf$conc2 ==0,]
  df0$viral_percent_inh <- platedf[1,]$virus_control_od
  analysisdf <- rbind(analysisdf, df0)
  
  # make a clean analysis dataframe
  df <- data.frame(matrix(nrow = nrow(analysisdf), ncol = 7))
  colnames(df) <- c("block_id", "drug1", "drug2", "conc1", "conc2", "response", "conc_unit") 
  
  #PairIndex
  df$block_id <- 1
  
  # drug name and concentration
  df$drug1 <- analysisdf[[paste0("drug", drugno_A, "_name")]]
  df$drug2 <- analysisdf[[paste0("drug", drugno_B, "_name")]]
  df$conc1 <- analysisdf[[paste0("drug", drugno_A, "_conc")]]
  df$conc2 <- analysisdf[[paste0("drug", drugno_B, "_conc")]]
  
  # drug response
  df$response<- analysisdf$viral_percent_inh
  
  # unit of concentration
  df$conc_unit <- ConcUnit
  
  df
}

############# plot the response
plot.response <- function(res,title){
  df <- data.frame(matrix(nrow = length(res$response_statistics$response_mean), ncol = 7))
  colnames(df) <- c("drug1","drug2","unit1","unit2","conc1","conc2","response")
  df$drug1 <- res$drug_pairs$drug1
  df$drug2 <- res$drug_pairs$drug2
  df$unit1 <- res$drug_pairs$conc_unit1
  df$unit2 <- res$drug_pairs$conc_unit2
  df$conc1 <- res$response_statistics$conc1
  df$conc2 <- res$response_statistics$conc2
  df$response <- round(res$response_statistics$response_mean)
  df$x <- 0
  df$y <- 0
  for(i in 1:nrow(df)){
    for(j in 1:length(unique(df$conc1))){
      if(df[i,]$conc1==unique(df$conc1)[j]){df[i,]$x <- j}
    }
  }
  for(i in 1:nrow(df)){
    for(j in 1:length(unique(df$conc2))){
      if(df[i,]$conc2==unique(df$conc2)[j]){df[i,]$y <- j}
    }
  }
  
  library(ggplot2)
  ggplot(df,aes(x=x, y=y, fill= response)) +
    # theme_classic()+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(colour = "black",size = 15),
          # legend.position = "none")+
          legend.key.width = unit(1.1, "lines"),
          legend.key.height = unit(1.3, "lines"),
          legend.key.size = unit(9,"pt"))+
    labs(x = paste0(df$drug1," (",df$unit1,")"), y = paste0(df$drug2," (",df$unit2,")"),
         title = paste(title,"Dose Response", sep = "\n")) + 
    geom_tile(colour = "black",size = .5) +
    # geom_text(aes(label = response), color = "black", size = 4) +
    scale_fill_gradient2(name = "% Inhibition",
                         low="#FBC85F",mid = "#FFFFFF", high="#757CBB",
                         breaks = seq(-100,100,50),limits = c(-150,150))+
    scale_x_continuous(breaks = c(1:length(unique(df$conc1))), labels = c(sort(round(unique(df$conc1),1))))+
    scale_y_continuous(breaks = c(1:length(unique(df$conc2))), labels = c(sort(round(unique(df$conc2),1))))
}

df.comb <- function(res,comb){
  df <- data.frame(matrix(nrow = 16, ncol = 7))
  colnames(df) <- c("comb","conc1","conc2","Bliss score", "Loewe score", "HSA score", "ZIP score")
  df$conc1 <- res$synergy_scores$conc1
  df$conc2 <- res$synergy_scores$conc2
  df$'Bliss score' <- res$synergy_scores$Bliss_synergy
  df$'HSA score' <- res$synergy_scores$HSA_synergy
  df$'Loewe score' <- res$synergy_scores$Loewe_synergy
  df$'ZIP score' <- res$synergy_scores$ZIP_synergy
  df$comb <- comb
  df
}
