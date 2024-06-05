getwd()
rm(list=ls())
####### STEP 1 Prepare the formatted data #########
# The following code uses the output from STEP 1 data formatting (see e.g. Wuhan_merged.csv)
# The formatting of cytotoxicity keeps the same as viral results
# the unit of concentration in the file is mg/L

####### function to analyse monotherapy data #########
# Prepare the functions for the analysis
source("monotherapy_function.R")

####### STEP 2 format the data per drug######
# for each of 1-4 drugs makes a dataset of just the monotherapy data
# use the extract.monotherapy function to format the data:
# mergedcsv - the name of the formatted file, drugno - same as the drug number (drugno) in formatted data,
# virus - the virus strain used in the assay, MW - molecular weight to change the concentration into uM
df1 <- extract.monotherapy(mergedcsv = "Wuhan_merged.csv", drugno = 1, virus = "Wuhan", MW = 602.6)
df2 <- extract.monotherapy(mergedcsv = "Wuhan_merged.csv", drugno = 2, virus = "Wuhan", MW = 157.10)
df3 <- extract.monotherapy(mergedcsv = "Wuhan_merged.csv", drugno = 3, virus = "Wuhan", MW = 499.5)
df4 <- extract.monotherapy(mergedcsv = "Wuhan_merged.csv", drugno = 4, virus = "Wuhan", MW = 259.22)
# combine the 4 extracted datasets together
dfmono <- rbind(df1, df2, df3, df4)
# 
#for only 2 drugs:
#dfmono <- rbind(df1, df2)

# extract the monotherapy data for cytotoxicity
df1.c <- extract.monotherapy(mergedcsv = "Wuhan_cytotoxicity_merged.csv", drugno = 1, 
                           virus = "cytotoxicity", MW = 602.6)
df2.c <- extract.monotherapy(mergedcsv = "Wuhan_cytotoxicity_merged.csv", drugno = 2, 
                           virus = "cytotoxicity", MW = 157.10)
df3.c <- extract.monotherapy(mergedcsv = "Wuhan_cytotoxicity_merged.csv", drugno = 3, 
                           virus = "cytotoxicity", MW = 499.5)
df4.c <- extract.monotherapy(mergedcsv = "Wuhan_cytotoxicity_merged.csv", drugno = 4, 
                           virus = "cytotoxicity", MW = 259.22)
dfcyto <- rbind(df1.c, df2.c, df3.c, df4.c)

######## plot raw monotherapy data ########
# plots the data
library(ggplot2)
pl.drug <- ggplot() +
  theme_bw() +
  scale_x_log10() +
  labs(x = "Drug concentration (mg/L)", y = "% viral inhibition",
       title = "") +
  geom_point(aes(x = conc, y = viral_percent_inh), data = dfmono) +
  # geom_line(aes(x = conc, y = mean_percent_inh), data = dfmono) +
  # geom_errorbar(aes(x = conc, y = mean_percent_inh, ymin = mean_percent_inh - sem_percent_inh,
  # ymax = mean_percent_inh + sem_percent_inh), data = dfmono) +
  facet_wrap(~ drug)
pl.drug

####### STEP 3 estimate EC50 ########
#  Now we can apply the est.hill function for each drug and look at the output.
# dfmono - the monotherapy data we just define in STEP 2
# gam - initial estimate for the Hill coefficient (affect slope of the curve), 
#       usually start with 1

# Output: xx_high, xx_low - 95% confidence interval of the parameter,
#         EC50_se - standard error of the EC50
#  The best fitting model has the lowest AIC BUT we need to check if the parameters make sense.
#  Only fit the model to data that makes sense
table1 = est.hill(dfmono = dfmono, drugno = 1, gam = 5)
table1
table2 = est.hill(dfmono = dfmono, drugno = 2, gam = 1.5)
table2
table3 = est.hill(dfmono = dfmono, drugno = 3, gam = 5)
table3
table4 = est.hill(dfmono = dfmono, drugno = 4, gam = 3)
table4
# combine the estimated results together
param <- rbind(table1, table2, table3, table4)
# save the estimation for later use, no need to run this part each time
write.csv(param, "Wuhan_table.csv")

###### STEP 4 plot the results #####
# read in the results from STEP 3
param <- read.csv("Wuhan_table.csv")

# use the lowest AIC model to plot monotherapy curve
# table - the result name from previous step, MW - molecular weight
pl.1 <- pl.mono(drugno = 1, table = param, MW = 602.6)
pl.2 <- pl.mono(drugno = 2, table = param, MW = 157.10)
pl.3 <- pl.mono(drugno = 3, table = param, MW = 499.5)
pl.4 <- pl.mono(drugno = 4, table = param, MW = 259.22)

# create a pdf file to save the four plots
pdf("Wuhan_monotherapy.pdf", height = 5.3, width = 5.5)
pl.1
pl.2
pl.3
pl.4
dev.off()

