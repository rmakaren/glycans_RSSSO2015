# GP - glycoprotein
# This script performs some EDA on the GP data from blood plasma. The dataset consits of
# mesurements from 31 plates, with additional variables as ethnics, age and gender, batch
# In order to exclude observation bias from such variables considering their distribution
# among the plates, the data was randomised and blocked

#set the working directory
setwd(dir = "/home/rostyslav/Desktop/RSSSO2015/R project/data")

#upload the libraries
library(dplyr)

## upload the dataset
df <- read.csv("plasma.csv")

###########################
# EXPLOATORY DATA ANALYSIS#
###########################

#verify the type of the data
lapply(df, class)

# 1. Ethnic 
sum(is.na(df$"ethnic")) # number of NA in ethnic column is 10 times higher that actual data
df$ethnic <- NULL

# 2. gender
table(df$"gender") # the females:males ratio is biased 10 times towards females
tapply(df$gender, df$Plate, table)

# 3. plate
df_Plate <- table(df$Plate)
length(table(df$Plate)) #31 Plates

#4. age
hist(df$age)

#####################
# Data visualisation#
#####################

#5. Variance distribution for each GP vs plate
############################################

pdf(file="GP_vs_Plates.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$Plate,
          main=glycan)
}
dev.off()

#6. Variance distribution for each GP vs gender
############################################

pdf(file="GP_vs_age.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$age,
          main=glycan)
}
dev.off()

#7. Variance distribution for age within the plates
################################################

pdf(file="age_vs_Plates.pdf")
boxplot(df$age ~ df$Plate)
dev.off()

# Excluding the males from the analysis due to observation bias
sub_df <- subset(df, gender == "F")

pdf(file="age_vs_Plates(no_males).pdf")
boxplot(sub_df$age ~ sub_df$Plate)
dev.off()

#8. Blocking and randomisation by gender and age
################################################

L<-list()
nrep <- 1000
for(ind in 1:nrep){
  #replication
  tmp <- sample (1:nrow(df), size = 140, replace = FALSE) # randomly pick up 140 samples
  dupl  <- df[tmp,] # duplicate picked samples
  dupl$gid <- paste0(dupl$gid, "_D")
  df_new <- rbind(df, dupl) # merge picked samples with orign dataframe
  
  #nplates <- 16
  nplates <- 31
  
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
boxplot(final$age~final$newplate)