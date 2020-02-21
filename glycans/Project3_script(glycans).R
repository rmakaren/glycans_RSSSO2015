library(dplyr)
library(stringr)
library(tidyr) #spread
library(preprocessCore) # quantile
df <-read.csv(file="areas-20141202.csv", header=TRUE, sep=",", dec=".")
df$narea <- NULL

# ============ Normalisation ====================
## Normalization for E-samples

df_E <- df %>%
  filter(grepl("E", SampleName))
df_E <- df_E %>% spread(glycan, area) #spread
df_E$SampleName <- NULL

# Total Volume Normalization (TVN)

df_tvn_E <- apply (df_E, MARGIN = 1, FUN = function(x) {x/sum(x)})
df_tvn_E <- as.data.frame(df_tvn_E)

# quantile 

df_tmp_E <- data.matrix(df_E)
df_q_E <- normalize.quantiles(df_tmp_E)
df_q_E <- as.data.frame(df_q_E)
colnames(df_q_E) <- names(df_E)

# median

df_med_E <- apply(df_E, 1, function(x){(x-median(x))/mad(x)})
tmp1 <- t(df_med_E)
df_med_E <- as.data.frame(tmp1)


#log ratio for each GP for each sample! Now we need to create a matrix with GPs as glycans
#log TVN
all.pairs <- combn(as.data.frame(t(df_tvn_E)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

TVNlog_E <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log q

all.pairs <- combn(as.data.frame(t(df_q_E)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Qlog_E <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log med

all.pairs <- combn(as.data.frame(t(df_med_E)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Mlog_E <- all.pairs %>%
  dplyr::select(glycan, ratio)

# Plotting results

TVNlog_E$norm = "TVN_E" 
Mlog_E$norm = "Med_E"
Qlog_E$norm = "QN_E"

X <- rbind_list(TVNlog_E, Mlog_E, Qlog_E)

library(ggplot2)
p <- ggplot(X, aes(x=ratio))
print(
  p
  + geom_density(aes(group=norm, color=norm, y=..scaled..))
  + facet_wrap(~glycan)
)

## Normalization for D samples

df_D <- df %>%
  filter(grepl("D", SampleName))
df_D <- df_D %>% spread(glycan, area) #spread
df_D$SampleName <- NULL

# Normalization TVN

df_tvn_D <- apply (df_D, MARGIN = 2, FUN = function(x) {x/sum(x)})
df_tvn_D <- as.data.frame(df_tvn_D)

# quantile 

df_tmp_D <- data.matrix(df_D)
df_q_D <- normalize.quantiles(df_tmp_D)
df_q_D <- as.data.frame(df_q_D)
colnames(df_q_D) <- names(df_D)

# median

df_med_D <- apply(df_D, 1, function(x){(x-median(x))/mad(x)})
tmp1 <- t(df_med_D)
df_med_D <- as.data.frame(tmp1)

#log ratio for each GP for each sample! Now we need to create a matrix with GPs as glycans
#log TVN
all.pairs <- combn(as.data.frame(t(df_tvn_D)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

TVNlog_D <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log q

all.pairs <- combn(as.data.frame(t(df_q_D)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Qlog_D <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log med

all.pairs <- combn(as.data.frame(t(df_med_D)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Mlog_D <- all.pairs %>%
  dplyr::select(glycan, ratio)

TVNlog_D$norm = "TVN_D" 
Mlog_D$norm = "Med_D"
Qlog_D$norm = "QN_D"

X <- rbind_list(TVNlog_D, Mlog_D, Qlog_D)

library(ggplot2)
p <- ggplot(X, aes(x=ratio))
print(
  p
  + geom_density(aes(group=norm, color=norm, y=..scaled..))
  + facet_wrap(~glycan)
)

## Normalization for C samples

df_C <- df %>%
  filter(grepl("C", SampleName))
df_C <- df_C %>% spread(glycan, area) #spread
df_C$SampleName <- NULL

# Normalization TVN

df_tvn_C <- apply (df_C, MARGIN = 2, FUN = function(x) {x/sum(x)})
df_tvn_C <- as.data.frame(df_tvn_C)

# quantile 

df_tmp_C <- data.matrix(df_C)
df_q_C <- normalize.quantiles(df_tmp_C)
df_q_C <- as.data.frame(df_q_C)
colnames(df_q_C) <- names(df_C)

# median

df_med_C <- apply(df_C, 1, function(x){(x-median(x))/mad(x)})
tmp1 <- t(df_med_C)
df_med_C <- as.data.frame(tmp1)

#log ratio for each GP for each sample! Now we need to create a matrix with GPs as glycans
#log TVN_C
all.pairs <- combn(as.data.frame(t(df_tvn_C)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

TVNlog_C <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log Q

all.pairs <- combn(as.data.frame(t(df_q_C)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Qlog_C <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log med

all.pairs <- combn(as.data.frame(t(df_med_C)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Mlog_C <- all.pairs %>%
  dplyr::select(glycan, ratio)

TVNlog_C$norm = "TVN_C" 
Mlog_C$norm = "Med_C"
Qlog_C$norm = "QN_C"

X <- rbind_list(TVNlog_C, Mlog_C, Qlog_C)

library(ggplot2)
p <- ggplot(X, aes(x=ratio))
print(
  p
  + geom_density(aes(group=norm, color=norm, y=..scaled..))
  + facet_wrap(~glycan)
)

## Normalization for B samples

df_B <- df %>%
  filter(grepl("B", SampleName))
df_B <- df_B %>% spread(glycan, area) #spread
df_B$SampleName <- NULL

# Normalization TVN

df_tvn_B <- apply (df_B, MARGIN = 2, FUN = function(x) {x/sum(x)})
df_tvn_B <- as.data.frame(df_tvn_B)

# quantile 

df_tmp_B <- data.matrix(df_B)
df_q_B <- normalize.quantiles(df_tmp_B)
df_q_B <- as.data.frame(df_q_B)
colnames(df_q_B) <- names(df_B)

# median

df_med_B <- apply(df_B, 1, function(x){(x-median(x))/mad(x)})
tmp1 <- t(df_med_B)
df_med_B <- as.data.frame(tmp1)

#log TVN
all.pairs <- combn(as.data.frame(t(df_tvn_B)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

TVNlog_B <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log Q

all.pairs <- combn(as.data.frame(t(df_q_B)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Qlog_B <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log med

all.pairs <- combn(as.data.frame(t(df_med_B)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Mlog_B <- all.pairs %>%
  dplyr::select(glycan, ratio)

TVNlog_B$norm = "TVN_B" 
Mlog_B$norm = "Med_B"
Qlog_B$norm = "QN_B"

X <- rbind_list(TVNlog_B, Mlog_B, Qlog_B)

library(ggplot2)
p <- ggplot(X, aes(x=ratio))
print(
  p
  + geom_density(aes(group=norm, color=norm, y=..scaled..))
  + facet_wrap(~glycan)
)

## Normalization for A samples

df_A <- df %>%
  filter(grepl("A", SampleName))
df_A <- df_A %>% spread(glycan, area) #spread
df_A$SampleName <- NULL

# TVN

df_tvn_A <- apply (df_A, MARGIN = 2, FUN = function(x) {x/sum(x)})
df_tvn_A <- as.data.frame(df_tvn_A)

# quantile 

df_tmp_A <- data.matrix(df_A)
df_q_A <- normalize.quantiles(df_tmp_A)
df_q_A <- as.data.frame(df_q_A)
colnames(df_q_A) <- names(df_A)

# median

df_med_A <- apply(df_A, 1, function(x){(x-median(x))/mad(x)})
tmp1 <- t(df_med_A)
df_med_A <- as.data.frame(tmp1)

#log TVN
all.pairs <- combn(as.data.frame(t(df_tvn_A)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

TVNlog_A <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log Q

all.pairs <- combn(as.data.frame(t(df_q_A)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Qlog_A <- all.pairs %>%
  dplyr::select(glycan, ratio)

# log med

all.pairs <- combn(as.data.frame(t(df_med_A)), m=2, simplify=FALSE) #return the combination of each 

tmp <- lapply(all.pairs, 
              function(X){
                names(X) <- c("V1", "V2")
                X$glycan <- rownames(X)
                return(X)
              })

all.pairs <- rbind_all(tmp)

all.pairs$ratio <- log(abs(all.pairs$V1/all.pairs$V2))

Mlog_A <- all.pairs %>%
  dplyr::select(glycan, ratio)

TVNlog_A$norm = "TVN_A" 
Mlog_A$norm = "Med_A"
Qlog_A$norm = "QN_A"

X <- rbind_list(TVNlog_A, Mlog_A, Qlog_A)

library(ggplot2)
p <- ggplot(X, aes(x=ratio))
print(
  p
  + geom_density(aes(group=norm, color=norm, y=..scaled..))
  + facet_wrap(~glycan)
)
