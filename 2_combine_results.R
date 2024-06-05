#setwd() to relevant folder getwd() to check (CAN DO THIS VIA SESSION TAB)
getwd()
rm(list=ls())

####### prepare the function #########
# Prepare the functions for the analysis
# other functions used are from package "synergyfinder" doi: 10.1016/j.gpb.2022.01.004.
source("2_combine_function.R")

####### STEP 1 reshape the data  #########
# read in the data
df <- read.csv("Wuhan_merged.csv")
# use the function to format the 2 drug combination data:
# drugno_A, drugno_B - in accordance with the drug number(drugno) in formatted data,
# ConcUnit - concentration unit
df1 <- extract.2combtherapy(df, drugno_A = 3, drugno_B = 1, ConcUnit = "mg/L")
df2 <- extract.2combtherapy(df, drugno_A = 3, drugno_B = 4, ConcUnit = "mg/L")
df3 <- extract.2combtherapy(df, drugno_A = 4, drugno_B = 1, ConcUnit = "mg/L")

# reshape the data to use synergyfinder
library(synergyfinder)
res1 <- ReshapeData(data = df1, data_type = "inhibition", impute = F)
res2 <- ReshapeData(data = df2, data_type = "inhibition", impute = F)
res3 <- ReshapeData(data = df3, data_type = "inhibition", impute = F)

####### STEP 2 Calculate the synergy scores  #########
# Calculate the synergy score with all four models
res1 <- CalculateSynergy(data = res1, method =  c("ZIP", "HSA", "Bliss", "Loewe"))
res2 <- CalculateSynergy(data = res2, method =  c("ZIP", "HSA", "Bliss", "Loewe"))
res3 <- CalculateSynergy(data = res3, method =  c("ZIP", "HSA", "Bliss", "Loewe"))

# save the results for later use
saveRDS(res1, file="Wuhan_R_N")
saveRDS(res2, file="Wuhan_E_N")
saveRDS(res3, file="Wuhan_R_E")

####### STEP 3 plot the response #########
# read in the results from previous step
res1 <- readRDS("Wuhan_R_N")
res2 <- readRDS("Wuhan_E_N")
res3 <- readRDS("Wuhan_R_E")

# create a pdf to save the plots
pdf("F2_Wuhan_response.pdf",width = 4.8, height = 4.2)
pl.r1 <- plot.response(res = res1, title = "Nirmatrelvir - Remdesivir")
pl.r1
pl.r2 <- plot.response(res = res2, title = "Nirmatrelvir - EIDD-1931")
pl.r2
pl.r3 <- plot.response(res = res3, title = "EIDD-1931 - Remdesivir")
pl.r3
dev.off()

####### STEP 4 plot the synergy scores #########
pdf("F2_Wuhan_synergy.pdf",width =6.7, height = 6.2)
pl.s1 <- Plot2DrugSurface(data = res1, plot_block = 1, drugs = c(1, 2),
                          plot_title = "Nirmatrelvir - Remdesivir\nSynergy Score",
                          plot_value = "Bliss_synergy",dynamic = FALSE,
                          high_value_color = "#757CBB",
                          low_value_color = "#FBC85F",
                          color_range = c(-100,100),
                          # z_range= c(-100,100),
                          text_size_scale = 1.5,
                          summary_statistic = NULL)

pl.s2 <- Plot2DrugSurface(data = res2, plot_block = 1, drugs = c(1, 2),
                          plot_title = "Nirmatrelvir - EIDD-1931\nSynergy Score",
                          plot_value = "Bliss_synergy",dynamic = FALSE,
                          high_value_color = "#757CBB",
                          low_value_color = "#FBC85F",
                          color_range = c(-100,100),
                          # z_range= c(-100,100),
                          text_size_scale = 1.5,
                          summary_statistic = NULL)

pl.s3 <- Plot2DrugSurface(data = res3, plot_block = 1, drugs = c(1, 2),
                          plot_title = "EIDD-1931 - Remdesivir\nSynergy Score",
                          plot_value = "Bliss_synergy",dynamic = FALSE,
                          high_value_color = "#757CBB",
                          low_value_color = "#FBC85F",
                          color_range = c(-100,100),
                          # z_range= c(-100,100),
                          text_size_scale = 1.5,
                          summary_statistic = NULL)

dev.off()

####### STEP 5 save the synergy scores #########
# create a table saving all the synergy scores
# comb - drug combination name
df1 <- df.comb(res1, comb = "Nirmatrelvir - Remdesivir")
df2 <- df.comb(res2, comb = "Nirmatrelvir - EIDD-1931")
df3 <- df.comb(res3, comb = "EIDD-1931 - Remdesivir")
df <- rbind(df1,df2,df3)
# save the results
write.csv(df, "Synergy_score_Wuhan.csv")
