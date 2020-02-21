library(ggplot2)

df <- read.csv("data2.csv")

# Total Volume Normalization (TVN):
df$Spot_alt <- NULL # removing spot_alt

df_new <- cbind.data.frame(df$X10m_Gel1_IS, df$X10m_Gel2_IS, df$X10m_Gel3_IS, df$X10m_Gel4_IS) #selecting only 10_IS)
colnames(df_new) <- c("X10m_Gel1_IS", "X10m_Gel2_IS", "X10m_Gel3_IS", "X10m_Gel4_IS") #assign the names to columns in df_new

df_tvn <- apply (df_new, MARGIN = 2, FUN = function(x) {x/sum(x)*median(x)}) # normalisation
df_tvn <- as.data.frame(df_tvn) # matrix into the dataframe

ggplot(df_tvn, aes(df_tvn$X10m_Gel1_IS, df_tvn$X10m_Gel2and4IS)) + 
  geom_point(aes(y = df_tvn$X10m_Gel2_IS, colour = "df_tvn$X10m_Gel2_IS")) +
  geom_point(aes(y = df_tvn$X10m_Gel3_IS, colour = "df_tvn$X10m_Gel3_IS")) +
  geom_point(aes(y = df_tvn$X10m_Gel4_IS, colour = "df_tvn$X10m_Gel4_IS")) 


# Median Normalisation by clusterSim package (Median)
# library(clusterSim)

# df_fin2 <- data.Normalization (df_med, type="n2",normalization="column")

#ggplot(df_med, aes(df_med$X10m_Gel1_IS, df_med$X10m_Gel2and4IS)) + 
  #geom_point(aes(y = df_med$X10m_Gel2_IS, colour = "df_med$X10m_Gel2_IS")) +
  #geom_point(aes(y = df_med$X10m_Gel3_IS, colour = "df_med$X10m_Gel3_IS")) +
  #geom_point(aes(y = df_med$X10m_Gel4_IS, colour = "df_med$X10m_Gel4_IS")) 

# code for Median

#(x-median)/mad

df_med <- apply(df_new, 1, function(x){x-median(x)/mad(x)})
df_med <- as.data.frame(t(df_med))

ggplot(df_med, aes(df_med$X10m_Gel1_IS, df_med$X10m_Gel2and4IS)) + 
  geom_point(aes(y = df_med$X10m_Gel2_IS, colour = "df_med$X10m_Gel2_IS")) +
  geom_point(aes(y = df_med$X10m_Gel3_IS, colour = "df_med$X10m_Gel3_IS")) +
  geom_point(aes(y = df_med$X10m_Gel4_IS, colour = "df_med$X10m_Gel4_IS")) 
  
# normalize.quantiles

df_rank <- apply(df_new,2,rank,ties.method="min")
df_sorted <- data.frame(apply(df_new, 2, sort))
df_mean <- apply(df_sorted, 1, mean)

index_to_mean <- function(my_index, my_mean){
  return(my_mean[my_index])
}

df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
rownames(df_final) <- toupper(letters[1:576])

quantile_normalisation <- function(df_new){
  df_rank <- apply(df_new,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df_new, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

df_q <- quantile_normalisation(df_new)
df_q <- data.frame(df_q)

ggplot(df_q, aes(df_q$X10m_Gel1_IS, df_q$X10m_Gel2and4IS)) + 
  geom_point(aes(y = df_q$X10m_Gel2_IS, colour = "df_q$X10m_Gel2_IS")) +
  geom_point(aes(y = df_q$X10m_Gel3_IS, colour = "df_q$X10m_Gel3_IS")) +
  geom_point(aes(y = df_q$X10m_Gel4_IS, colour = "df_q$X10m_Gel4_IS")) 

# VSN

# source("http://bioconductor.org/biocLite.R")
#biocLite("vsn")
library("vsn")

df_vsn <- vsn(df_new)

exprs(df_vsn[1:576,])
df_vsn <- as.data.frame(exprs(df_vsn[1:576,]))


ggplot(df_vsn, aes(df_vsn$X10m_Gel1_IS, df_vsn$X10m_Gel2and4IS)) + 
  geom_point(aes(y = df_vsn$X10m_Gel2_IS, colour = "df_vsn$X10m_Gel2_IS")) +
  geom_point(aes(y = df_vsn$X10m_Gel3_IS, colour = "df_vsn$X10m_Gel3_IS")) +
  geom_point(aes(y = df_vsn$X10m_Gel4_IS, colour = "df_vsn$X10m_Gel4_IS")) 

# logariphm for TVN df_tvn

r1r2 <- log(df_tvn$X10m_Gel1_IS/df_tvn$X10m_Gel2_IS)
r1r3 <- log(df_tvn$X10m_Gel1_IS/df_tvn$X10m_Gel3_IS)
r1r4 <- log(df_tvn$X10m_Gel1_IS/df_tvn$X10m_Gel4_IS)
r2r3 <- log(df_tvn$X10m_Gel2_IS/df_tvn$X10m_Gel3_IS)
r2r4 <- log(df_tvn$X10m_Gel2_IS/df_tvn$X10m_Gel4_IS)
r3r4 <- log(df_tvn$X10m_Gel3_IS/df_tvn$X10m_Gel4_IS)

log_tvn <- c(r1r2, r1r3, r1r4, r2r3, r2r4, r3r4)

# logariphm for Median df_tvn

r1r2 <- log(df_med$X10m_Gel1_IS/df_med$X10m_Gel2_IS)
r1r3 <- log(df_med$X10m_Gel1_IS/df_med$X10m_Gel3_IS)
r1r4 <- log(df_med$X10m_Gel1_IS/df_med$X10m_Gel4_IS)
r2r3 <- log(df_med$X10m_Gel2_IS/df_med$X10m_Gel3_IS)
r2r4 <- log(df_med$X10m_Gel2_IS/df_med$X10m_Gel4_IS)
r3r4 <- log(df_med$X10m_Gel3_IS/df_med$X10m_Gel4_IS)

log_med <- c(r1r2, r1r3, r1r4, r2r3, r2r4, r3r4)

# logariphm for Quantile df_q

r1r2 <- log(df_q$X10m_Gel1_IS/df_q$X10m_Gel2_IS)
r1r3 <- log(df_q$X10m_Gel1_IS/df_q$X10m_Gel3_IS)
r1r4 <- log(df_q$X10m_Gel1_IS/df_q$X10m_Gel4_IS)
r2r3 <- log(df_q$X10m_Gel2_IS/df_q$X10m_Gel3_IS)
r2r4 <- log(df_q$X10m_Gel2_IS/df_q$X10m_Gel4_IS)
r3r4 <- log(df_q$X10m_Gel3_IS/df_q$X10m_Gel4_IS)

log_q <- c(r1r2, r1r3, r1r4, r2r3, r2r4, r3r4)

# logarighm for VSN df_vsn

r1r2 <- log(df_vsn$X10m_Gel1_IS/df_vsn$X10m_Gel2_IS)
r1r3 <- log(df_vsn$X10m_Gel1_IS/df_vsn$X10m_Gel3_IS)
r1r4 <- log(df_vsn$X10m_Gel1_IS/df_vsn$X10m_Gel4_IS)
r2r3 <- log(df_vsn$X10m_Gel2_IS/df_vsn$X10m_Gel3_IS)
r2r4 <- log(df_vsn$X10m_Gel2_IS/df_vsn$X10m_Gel4_IS)
r3r4 <- log(df_vsn$X10m_Gel3_IS/df_vsn$X10m_Gel4_IS)

log_vsn <- c(r1r2, r1r3, r1r4, r2r3, r2r4, r3r4)

log_all <- c(log_tvn, log_med, log_q, log_vsn)

#  Plotting all of the log

plot(density(log_tvn), col = "blue")
lines(density(log_med), col = "red")
lines(density(log_q), col = "black")
lines(density(log_vsn), col = "orange")
legend(x = 2, y = 2, c("log_tvn", "log_med", "log_q", "log_vsn"), 
       col = par("blue", "red", "black", "orange"), cex = 8)

# VSN for data

df_carb <- cbind.data.frame(df$X10m_Gel1_Carbonyl..Cy5, df$X10m_Gel2_Carbonyl..Cy5, df$X10m_Gel3_Carbonyl..Cy5, df$X10m_Gel4_Carbonyl..Cy5)
library("vsn")
df_vsn_c10 <- justvsn(data.matrix(df_carb))

df_vsn_c10 <- as.data.frame(df_vsn_c10)
colnames(df_vsn_c10) <- c("X10m_Gel1_Carbonyl..Cy5", "X10m_Gel2_Carbonyl..Cy5", "X10m_Gel3_Carbonyl..Cy5", "X10m_Gel4_Carbonyl..Cy5")

ggplot(df_vsn_c10, aes(x=df_vsn_c10$X10m_Gel1_Carbonyl..Cy5)) + 
  geom_point(aes(y = df_vsn_c10$X10m_Gel2_Carbonyl..Cy5, colour = "df_vsn_c10$X10m_Gel2_Carbonyl..Cy5")) +
  geom_point(aes(y = df_vsn_c10$X10m_Gel3_Carbonyl..Cy5, colour = "df_vsn_c10$X10m_Gel3_Carbonyl..Cy5")) +
  geom_point(aes(y = df_vsn_c10$X10m_Gel4_Carbonyl..Cy5, colour = "df_vsn_c10$X10m_Gel4_Carbonyl..Cy5")) 
