####### function to analyse monotherapy data #########
#   For monotherapies we will only extract the data for each drug where it was given without
#   any of the other drugs
extract.monotherapy <- function(mergedcsv, drugno, virus, MW){
  # mergedcsv <- "Wuhan_merged.csv"
  # drugno <- 1
  platedf <- read.csv(mergedcsv)
  # work out which of the 4 drugs the user has NOT chosen
  exclude <- which(c(1:4) != drugno)
  # select nonzero drug wells
  platedf <- platedf[platedf[paste0("drug", drugno, "_conc")] != 0, ]
  # remove nonzero combination wells
  platedf <- platedf[platedf[paste0("drug", exclude[1], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[2], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[3], "_conc")] == 0, ]
  # make a clean analysis dataframe
  analysisdf <- data.frame(matrix(nrow = nrow(platedf), ncol = 8))
  colnames(analysisdf) <- c("virus", "drugno", "drug", "conc_mg", "conc","viral_percent_inh", 
                            "mean_percent_inh", "sem_percent_inh" ) 
  # name the concentration in mg/L to conc_mg
  analysisdf$conc_mg <-  platedf[[paste0("drug", drugno, "_conc")]]
  # calculate the concentration in uM
  analysisdf$conc <-  analysisdf$conc_mg/MW * 10^3
  analysisdf$viral_percent_inh <-  platedf$viral_percent_inh
  #define standard error of mean function
  std.error <- function(x) sd(x)/sqrt(length(x))
  # calculate the mean and sem for each concentration
  analysisdf$mean_percent_inh <- -999
  analysisdf$sem_percent_inh <- -999
  for(i in unique(analysisdf$conc)){
    analysisdf$mean_percent_inh[analysisdf$conc == i] <- mean(analysisdf$viral_percent_inh[analysisdf$conc == i])
    analysisdf$sem_percent_inh[analysisdf$conc == i] <- std.error(analysisdf$viral_percent_inh[analysisdf$conc == i]) 
  }
  # include other information
  analysisdf$drug <- platedf[[paste0("drug", drugno, "_name")]]
  analysisdf$drugno <- drugno
  analysisdf$virus <- virus
  analysisdf
}

####### estimate EC50 ########
# takes the merged data from 1 or more replicates, and calculates EC50 by:
#  1. estimating only EC50
#  2. estimating EC50 and gamma
#  3. estimating EC50 and Emax
#  4. estimating EC50, emax and gamma
est.hill <- function(dfmono, drugno, gam){
  # specify the 8 models:
  m1 <- y ~  100 * x / (EC50 + x)
  m2 <- y ~  100 * x ^ gam / (EC50 ^ gam + x ^ gam)
  m3 <- y ~  Emax * x / (EC50  + x)
  m4 <- y ~  E0 + 100 * x / (EC50  + x)
  m5 <- y ~  E0 + Emax * x / (EC50  + x)
  m6 <- y ~  Emax * x ^ gam / (EC50 ^ gam  + x ^gam)
  m7 <- y ~  E0 + 100 * x ^ gam / (EC50 ^ gam  + x ^ gam)
  m8 <- y ~  E0 + Emax * x ^ gam / (EC50 ^ gam  + x ^ gam)
  
  df <- dfmono[dfmono$drugno == drugno, ]
  # 1. take EC50 initial est as half of the conc range
  EC50_init1 <- median(df$conc)
  x <- df$conc
  y <- df$viral_percent_inh
  r1 <- try(nls(formula = m1, start = c(EC50 = EC50_init1)),
            silent = TRUE)
  # update the initial estimate for model 2
  EC50_init2 <- summary(r1)$coefficients[1]
  # model 2 estimate gamma
  r2 <- try(nls(formula = m2, start = c(EC50 = EC50_init2, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 3 estimate Emax
  r3 <- try(nls(formula = m3, start = c(EC50 = EC50_init2, Emax = 100), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 4 estimate E0
  r4 <- try(nls(formula = m4, start = c(EC50 = EC50_init2, E0 = 0), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 5 estimate Emax and E0
  r5 <- try(nls(formula = m5, start = c(EC50 = EC50_init2, E0 = 0, Emax = 100), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 6 estimate Emax and gamma
  r6 <- try(nls(formula = m6, start = c(EC50 = EC50_init2, Emax = 100, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 7 estimate E0 and gam
  r7 <- try(nls(formula = m7, start = c(EC50 = EC50_init2, E0 = 0, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  # model 8 estimate emax E0 and gam
  r8 <- try(nls(formula = m8, start = c(EC50 = EC50_init2, Emax = 100, E0 = 0, gam = gam)),
            silent = TRUE)
  # make a dataframe for the output
  out <- data.frame(matrix(nrow = 8, ncol = 17))
  colnames(out) <- c("drugno","drug","model", "AIC", "EC50","EC50_low", "EC50_high","EC50_se",
                     "gam", "gam_low", "gam_high","Emax","Emax_low", "Emax_high",
                     "E0", "E0_low", "E0_high")
  out$drugno <- drugno
  out$drug <- df[1,]$drug
  out$model <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8")
  out$AIC[1] <- ifelse(is.numeric(try(AIC(r1), silent = TRUE)) == TRUE, AIC(r1), "fail")
  out$AIC[2] <- ifelse(is.numeric(try(AIC(r2), silent = TRUE)) == TRUE, AIC(r2), "fail")
  out$AIC[3] <- ifelse(is.numeric(try(AIC(r3), silent = TRUE)) == TRUE, AIC(r3), "fail")
  out$AIC[4] <- ifelse(is.numeric(try(AIC(r4), silent = TRUE)) == TRUE, AIC(r4), "fail")
  out$AIC[5] <- ifelse(is.numeric(try(AIC(r5), silent = TRUE)) == TRUE, AIC(r5), "fail")
  out$AIC[6] <- ifelse(is.numeric(try(AIC(r6), silent = TRUE)) == TRUE, AIC(r6), "fail")
  out$AIC[7] <- ifelse(is.numeric(try(AIC(r7), silent = TRUE)) == TRUE, AIC(r7), "fail")
  out$AIC[8] <- ifelse(is.numeric(try(AIC(r8), silent = TRUE)) == TRUE, AIC(r8), "fail")
  # EC50
  out$EC50[1] <- ifelse(out$AIC[1] == "fail", "fail", summary(r1)$coefficients[1])
  out$EC50_se[1] <- ifelse(out$AIC[1] == "fail", "fail", summary(r1)$coefficients[1,2])
  out$EC50_low[1] <- ifelse(is.numeric(try(confint2(r1), silent = TRUE)) == TRUE, confint2(r1)[1], "fail")
  out$EC50_high[1] <- ifelse(is.numeric(try(confint2(r1), silent = TRUE)) == TRUE, confint2(r1)[2], "fail")
  out$EC50[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[1])
  out$EC50_se[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[1,2])
  out$EC50_low[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[1], "fail")
  out$EC50_high[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[3], "fail")
  out$EC50[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[1])
  out$EC50_se[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[1,2])
  out$EC50_low[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[1], "fail")
  out$EC50_high[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[3], "fail")
  out$EC50[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[1])
  out$EC50_se[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[1,2])
  out$EC50_low[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[1], "fail")
  out$EC50_high[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[3], "fail")
  out$EC50[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[1])
  out$EC50_se[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[1,2])
  out$EC50_low[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[1], "fail")
  out$EC50_high[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[4], "fail")
  out$EC50[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[1])
  out$EC50_se[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[1,2])
  out$EC50_low[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[1], "fail")
  out$EC50_high[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[4], "fail")
  out$EC50[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[1])
  out$EC50_se[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[1,2])
  out$EC50_low[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[1], "fail")
  out$EC50_high[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[4], "fail")
  out$EC50[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[1])
  out$EC50_se[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[1,2])
  out$EC50_low[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[1], "fail")
  out$EC50_high[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[5], "fail")
  # gam
  out$gam[1] <- ifelse(out$AIC[1] == "fail", "fail", 1)
  out$gam_low[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$gam_high[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$gam[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[2])
  out$gam_low[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[2], "fail")
  out$gam_high[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[4], "fail")
  out$gam[3] <- ifelse(out$AIC[3] == "fail", "fail", 1)
  out$gam_low[3] <- ifelse(out$AIC[3] == "fail", "fail", "fixed")
  out$gam_high[3] <- ifelse(out$AIC[3] == "fail", "fail", "fixed")
  out$gam[4] <- ifelse(out$AIC[4] == "fail", "fail", 1)
  out$gam_low[4] <- ifelse(out$AIC[4] == "fail", "fail", "fixed")
  out$gam_high[4] <- ifelse(out$AIC[4] == "fail", "fail", "fixed")
  out$gam[5] <- ifelse(out$AIC[5] == "fail", "fail", 1)
  out$gam_low[5] <- ifelse(out$AIC[5] == "fail", "fail", "fixed")
  out$gam_high[5] <- ifelse(out$AIC[5] == "fail", "fail", "fixed")
  out$gam[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[3])
  out$gam_low[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[3], "fail")
  out$gam_high[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[6], "fail")
  out$gam[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[3])
  out$gam_low[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[3], "fail")
  out$gam_high[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[6], "fail")
  out$gam[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[4])
  out$gam_low[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[4], "fail")
  out$gam_high[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[8], "fail")
  # Emax
  out$Emax[1] <- ifelse(out$AIC[1] == "fail", "fail", 100)
  out$Emax_low[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$Emax_high[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$Emax[2] <- ifelse(out$AIC[2] == "fail", "fail", 100)
  out$Emax_low[2] <- ifelse(out$AIC[2] == "fail", "fail", "fixed")
  out$Emax_high[2] <- ifelse(out$AIC[2] == "fail", "fail", "fixed")
  out$Emax[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[2])
  out$Emax_low[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[2], "fail")
  out$Emax_high[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[4], "fail")
  out$Emax[4] <- ifelse(out$AIC[4] == "fail", "fail", 100)
  out$Emax_low[4] <- ifelse(out$AIC[4] == "fail", "fail", "fixed")
  out$Emax_high[4] <- ifelse(out$AIC[4] == "fail", "fail", "fixed")
  out$Emax[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[3])
  out$Emax_low[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[3], "fail")
  out$Emax_high[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[6], "fail")
  out$Emax[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[2])
  out$Emax_low[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[2], "fail")
  out$Emax_high[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[5], "fail")
  out$Emax[7] <- ifelse(out$AIC[7] == "fail", "fail", 100)
  out$Emax_low[7] <- ifelse(out$AIC[7] == "fail", "fail", "fixed")
  out$Emax_high[7] <- ifelse(out$AIC[7] == "fail", "fail", "fixed")
  out$Emax[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[2])
  out$Emax_low[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[2], "fail")
  out$Emax_high[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[6], "fail")
  # E0
  out$E0[1] <- ifelse(out$AIC[1] == "fail", "fail", 0)
  out$E0_low[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$E0_high[1] <- ifelse(out$AIC[1] == "fail", "fail", "fixed")
  out$E0[2] <- ifelse(out$AIC[2] == "fail", "fail", 0)
  out$E0_low[2] <- ifelse(out$AIC[2] == "fail", "fail", "fixed")
  out$E0_high[2] <- ifelse(out$AIC[2] == "fail", "fail", "fixed")
  out$E0[3] <- ifelse(out$AIC[3] == "fail", "fail", 0)
  out$E0_low[3] <- ifelse(out$AIC[3] == "fail", "fail", "fixed")
  out$E0_high[3] <- ifelse(out$AIC[3] == "fail", "fail", "fixed")
  out$E0[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[2])
  out$E0_low[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[2], "fail")
  out$E0_high[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[4], "fail")
  out$E0[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[2])
  out$E0_low[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[2], "fail")
  out$E0_high[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[5], "fail")
  out$E0[6] <- ifelse(out$AIC[6] == "fail", "fail", 0)
  out$E0_low[6] <- ifelse(out$AIC[6] == "fail", "fail", "fixed")
  out$E0_high[6] <- ifelse(out$AIC[6] == "fail", "fail", "fixed")
  out$E0[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[2])
  out$E0_low[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[2], "fail")
  out$E0_high[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[5], "fail")
  out$E0[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[3])
  out$E0_low[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[3], "fail")
  out$E0_high[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[7], "fail")
  out
}

###### Function to plot the model #####
pl.mono <- function(drugno, table, MW){
  #drugno <- 2, table <- param
  # prepare the data for plot
  dfpred <- dfmono[dfmono$drugno == drugno, ]
  dftoxi <- dfcyto[dfcyto$drugno == drugno, ]
  dftoxi$drug <- "Cytotoxicity"
  dfplot <- rbind(dfpred, dftoxi)
  # use the drug name from the dataset as title for the plot
  drugname <- dfpred$drug[1]
  #  Generate a range of concentrations starting lower and ending higher than those observed
  cmin <- min(dfpred$conc) 
  cmax <- max(dfpred$conc) 
  cpreds <- seq(cmin/2, cmax*2, 0.1)
  crange <- cmax-cmin
  # select the model with lowest AIC
  tableparam <- table[table$drugno == drugno,]
  estparam <- tableparam[tableparam$AIC == min(tableparam$AIC),]
  EC50 <- as.numeric(estparam$EC50)
  E0 <- as.numeric(estparam$E0)
  Emax <- as.numeric(estparam$Emax)
  gam <- as.numeric(estparam$gam)
  # make a vector of predictions
  hill_pred <- E0 + Emax * cpreds ^ gam / (EC50 ^ gam + cpreds ^ gam)
  # prepare the annotation
  ec50 <- round(as.numeric(estparam$EC50),2)
  
  library(ggplot2)
  library(ggprism)
  pl.fit <- ggplot() +
    theme_prism() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
          text = element_text(size = 15, face = "bold"))+
    guides(colour = "none")+
    labs(y = "% viral inhibition",title = drugname) +
    geom_point(aes(x = conc, y = mean_percent_inh, colour = drug), data = dfplot, size=1.5)+
    geom_line(aes(x = conc, y = mean_percent_inh), data = dftoxi, colour = "#ED7170", linewidth=1)+
    geom_line(aes(x = cpreds, y = hill_pred), colour = "#757CBB", linewidth = 1)+
    geom_errorbar(aes(x = conc, y = mean_percent_inh, ymin = mean_percent_inh - sem_percent_inh, 
                      ymax = mean_percent_inh + sem_percent_inh, colour = drug), data = dfplot, width=.1) +
    scale_colour_manual(name = NULL, values = c("#ED7170", "#757CBB"))+
    scale_x_continuous(trans = 'log2',
                       breaks = 2^seq(log2(cmin),log2(cmax),1),
                       labels = ~paste(signif(.,2), round(.* MW / 10^3,2), sep = "\n"),
                       name = "Drug concentration")+
    scale_y_continuous(breaks = seq(0,100,50))+
    annotate("text", x = 0.03*crange, y = 115, size = 6,
             label = paste0("EC50 = ",round(as.numeric(estparam$EC50),1),"uM"))+
    annotate("text", x = cmax*2, y = -13, label = bquote(bold(paste("(",mu,"M)", sep = ""))), size = 5)+
    annotate("text", x = cmax*2, y = -21, label = "(mg/L)", size = 5, fontface = "bold")+
    coord_cartesian(ylim = c(0,120), clip = "off")+
    geom_hline(aes(yintercept = 50), colour="grey", linetype="dashed")
}
