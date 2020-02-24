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

#Ethnic 
sum(is.na(df$"ethnic")) # number of NA in ethnic column is 10 times higher that actual data
df$ethnic <- NULL

# gender
table(df$"gender") # the females:males ratio is biased 10 times towards females
tapply(df$gender, df$Plate, table)

#plate
df_Plate <- table(df$Plate)
length(table(df$Plate)) #31 Plates

#age
hist(df$age)

#####################
# Data visualisation#
#####################

#Creating boxplot graphics for each glycans gender vs GP
pdf(file="GP_vs_Plates.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$Plate,
          main=glycan)
}
dev.off()

#Creating boxplot graphics for each glycans gender vs GP
pdf(file="GP_vs_age.pdf")
glycans <- grep("GP", names(df), value=TRUE)

for(glycan in glycans){
  boxplot(df[[glycan]] ~ df$age,
          main=glycan)
}
dev.off()

# Excluding the males from the analysis due to number bias
sub_df <- subset(df, gender == "F")

pdf(file="age_vs_Plates.pdf")
boxplot(df$age ~ df$Plate)
dev.off()

pdf(file="age_vs_Plates(no_males).pdf")
boxplot(sub_df$age ~ sub_df$Plate)
dev.off()

#Blocking and randomisation by gender and age
df_f <- df[df$gender == "F", ] # selected by Females
df_f <- df_f[order(df_f$age),] # sorted by age in Females dataframe
a = rep(seq(1:31), 42)[1:1278] # created sequence of 31 plates for 42 samples
nrow(df_f)
length(a)
df_f$new_plate = a # assign new plates
boxplot(df_f$age~df_f$new_plate)

df_m <- df[df$gender == "M", ] # selected by Males
df_m <- df_m[order(df_m$age),] # sorted by age in Males dataframe
a = rep(seq(1:31), 4)[1:114] # created sequence of 31 plates for 4 samples
nrow(df_m)
length(a)
df_m$new_plate = a # assign new plates
boxplot(df_m$age~df_m$new_plate)

df_new <- rbind(df_f, df_m)
boxplot(df_new$age~df_new$new_plate)

# Making experiment design again
# replication

df <- read.csv("plasma.csv")

0.1 * 1392 # 10% from orign dataframe is 140

L<-list()
nrep <- 1000
for(ind in 1:nrep){
  #replication
  tmp <- sample (1:nrow(df), size = 140, replace = FALSE) # randomly pick up 140 samples
  dupl  <- df[tmp,] # duplicate picked samples
  dupl$gid <- paste0(dupl$gid, "_D")
  df_new <- rbind(df, dupl) # merge picked samples with orign dataframe
  
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
?as.data.frame()

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
boxplot(final$age~final$newplate)