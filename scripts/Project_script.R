library(tidyverse)
library(naniar)

setwd()

## Project work
df <- read.csv("00-scripts-data/data/plasma.csv", stringsAsFactors = T)

#########
# 1. EDA
#########

#Ethnic 
factor(df$ethnic) # df contains only information about white people
length(df$ethnic)
sum(is.na(df$ethnic)) # number of NA data 10 times more that present

# gender
table(df$gender) # number of females is 10 times greater males
table(df$gender, df$Plate)

# Visualizing every glycans by plate
pdf(file="01-Plasma_block_rand/EDA/GP_vs_Plates.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$Plate,
          main=glycan,
          xlab = "plate",
          ylab = "glycan expression")
}
dev.off()

#Creating boxplot graphics for each glycans gender vs GP
pdf(file="01-Plasma_block_rand/EDA/GP_vs_age.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$age,
          data = df,
          main=glycan,
          xlab = "age",
          ylab = "glycan expression")
}
dev.off()

# removing males 
sub_df <- subset(df, gender == "F")

pdf(file="01-Plasma_block_rand/EDA/age_vs_Plates.pdf")
boxplot(df$age ~ df$Plate,
        xlab = "plate",
        ylab = "age")
dev.off()

pdf(file="01-Plasma_block_rand/EDA/age_vs_Plates(no_males).pdf")
boxplot(sub_df$age ~ sub_df$Plate,
        xlab = "plate",
        ylab = "age")
dev.off()

####################################
# 2. Blocking and randomization
####################################

#Blocking and randomization by gender and age
df_f <- df[df$gender == "F", ] # selected by Females
df_f <- df_f[order(df_f$age),] # sorted by age in Females dataframe
a = rep(seq(1:31), 42)[1:1278] # created sequence of 31 plates for 42 samples
nrow(df_f)
length(a)
df_f$new_plate = a # assign new plates

pdf(file="01-Plasma_block_rand/Block_rand/age_vs_Plates(F)_newplate.pdf")
boxplot(age~new_plate, 
        data = df_f)
dev.off()

df_m <- df[df$gender == "M", ] # selected by Males
df_m <- df_m[order(df_m$age),] # sorted by age in Males dataframe
a = rep(seq(1:31), 4)[1:114] # created sequence of 31 plates for 4 samples
nrow(df_m)
length(a)
df_m$new_plate = a # assign new plates

pdf(file="01-Plasma_block_rand/Block_rand/age_vs_Plates(M)_newplate.pdf")
boxplot(age ~ new_plate,
        data = df_m)
dev.off()

pdf(file="01-Plasma_block_rand/Block_rand/age_vs_Plates(Total)_newplate.pdf")
df_new <- rbind(df_f, df_m)
boxplot(age~new_plate,
        data=df_new)
dev.off()

##################################################
# 3. Making experiment design again: replication
##################################################

df <- read.csv("/media/rostyslav/Toshiba/Projects/RSSSO_glycans/data/plasma.csv")

0.1 * 1392 # 10% from original dataframe is 140

L<-list()
nrep <- 1000
for(ind in 1:nrep){
  #replication
  tmp <- sample (1:nrow(df), size = 140, replace = FALSE) # randomly pick up 140 samples
  dupl  <- df[tmp,] # duplicate picked samples
  dupl$gid <- paste0(dupl$gid, "_D")
  df_new <- rbind(df, dupl) # merge picked samples with origin dataframe
  
  nplates <- 16
  
  #blocking
  df_f <- df_new[df_new$gender == "F", ] # selected by Females from orign
  
  nfp  <- floor(nrow(df_f)/nplates) # cealing the number of Female samples by plates to be created
  df_f <- df_f[sample(1:nrow(df_f)), ] # shuffle dataframe of Female samples
  df_f$newplate <- rep(1:nplates, length=nrow(df_f)) # create new plates and put samples into them
  head(df_f)
  table(df_f$newplate)
  
  df_m <- df_new[df_new$gender == "M", ] # selected by Males from orign
  
  # calculate free space
  
  nfp  <- ceiling(nrow(df_m)/nplates) # cealing the number of Male samples by plates to be created
  df_m <- df_m[sample(1:nrow(df_m)), ] # shuffle dataframe of Male samples
  df_m$newplate <- rep(1:nplates, length=nrow(df_m)) # create new plates and put samples into them
  head(df_m)
  table(df_m$newplate)
  
  df_new_D <- rbind(df_f, df_m)
  
  L[[ind]] <- df_new_D
}

summary(L[1:10])
x <-as.data.frame(L[1])
summary(x)

lapply(L, quantile, probs = 1:3/4, n= TRUE)

x <- rep(0, nrep)
for (i in 1:length(L)){
  elem <- L[[i]]
  tmp <- elem %>%
    group_by(newplate) %>%
    summarise(m=mean(age)) %>%
    ungroup()
  
  x[i] <- var(tmp$m)
}

ind <- which.min(x)

final <- L[[ind]]

pdf(file="01-Plasma_block_rand/DoF_corrected/corrected_plasma.pdf")
boxplot(age~newplate,
        data=final)
dev.off()

write.csv(final, "01-Plasma_block_rand/DoF_corrected/corrected_plasma.csv")
