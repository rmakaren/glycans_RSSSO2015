library(dplyr)
library(stringr)
library(tidyr) #spread
library(preprocessCore)
library(ggplot2)
library(sva)
library(lme4)
install.packages("lme4")
dat <- read.csv("ultimate_igg-areas-20150417.csv")
dat <- dat[, -(c(5:7))]

dat <- dat[complete.cases(dat),]

# Apply function for normalization

f <- function(dat){
#s  MedN <- medianNorm(dat)
  MedQM <- medianQuotientNorm(dat)
  QN <- quantileNorm(dat)
  TAN <- totalAreaNorm(dat)
  
  rbind_list( MedQM, QN, TAN)
}

####### grouping by A, B, C

tmp <- do.call(rbind, str_split(dat$gid, "_"))[,1]
dat$uid <- tmp

Norm <- dat %>%
  group_by(Plate) %>%
  do(f(.)) %>%
  ungroup()

dat$normalization <- "None"

Norm <- rbind(Norm, dat)


# plot the results

MedNr_A <- MedNr %>%
  filter(grepl("A", uid))

p <- ggplot(pNorm, aes(x=ratio))
print(
  p
  + geom_density(aes(group=normalization,
                     color=normalization, 
                     y=..scaled..))
  + facet_wrap(~glycan)
)

###### pair all logs

all.pairs <- function(X){ 
  glycans <- grepl("GP", names(X)) 
  tmpX <- X[, glycans] 
  tX <- as.data.frame(t(tmpX)) 
  X.pairs <- combn(tX, m=2, simplify=FALSE) # enforce same column names! 
  X.pairs <- lapply(X.pairs, function(Y){ 
    names(Y) <- c("V1", "V2") 
    Y$glycan <- rownames(Y) 
    Y 
    }) 
  X.pairs <- rbind_all(X.pairs) 
  X.pairs <- X.pairs %>% mutate(ratio=log(abs(V1/V2))) 
  tmp <- X[, !glycans] 
  X.pairs <- merge(X.pairs, tmp[1, c("uid", "normalization")]) 
  return(X.pairs) 
}

pNorm = Norm %>% 
  group_by(uid, normalization) %>% 
  do(all.pairs(.)) %>% 
  ungroup()


#########################
tmp <- subset(B_Norm, B_Norm$normalization == "QN")

# Plotting the results

boxplot(B_Norm$GP1~B_Norm$Plate)
boxplot(Norm$GP1~Norm$Plate)

head(B_Norm)

X <- B_Norm %>%
  gather(glycans, area, GP1:GP24)

p <- ggplot(X, aes(x=factor(Plate), y=area))
print(
  p
  + geom_violin()
  + facet_wrap(~glycans+normalization,
               scales="free",
               ncol=3)
)


MedNr_A <- MedNr %>%
  filter(grepl("A", uid))

p <- ggplot(pNorm, aes(x=ratio))
print(
  p
  + geom_density(aes(group=normalization,
                     color=normalization, 
                     y=..scaled..))
  + facet_wrap(~glycan)
)


###### searching for batch effect
#source("http://bioconductor.org/biocLite.R")
#biocLite("sva")

B_Norm <-  Norm %>%
  group_by(normalization) %>%
  do(empirical_bayes_BC(.)) %>%
  ungroup()


##### Linear mixed model
glycans <- paste("GP",1:24,sep="")
matrix()

td <- Norm %>%
  gather(glycan, area, GP1:GP24)
 
f <- function(X){
  m <- lmer(area ~ 1 + (1|labPeriod) + (1|labPerson) + (1|uid),
              data=X)
  data.frame(VarCorr(m))
}

final <- td %>%
  group_by(glycan, normalization) %>%
  do(f(.)) %>%
  ungroup()

final <- final %>%
  select(glycan, normalization, grp, sdcor)

###########

td <- normData %>%
  gather(glycan, area, GP1:GP24)

f <- function(X){
  m <-lmer(area ~1 + (1|uid) 
#            +(1|labPerson) 
#            + (1|labPeriod), 
            + (1|Plate),
           data=X)
  return(as.data.frame(VarCorr(m)))
}

tmp <- td %>%
  group_by(glycan, normalization) %>%
  do(f(.)) %>%
  ungroup()

tmp <- tmp %>%
  select(glycan, normalization, grp, vcov, sdcor)

tmp$grp <- factor(tmp$grp, 
                  # levels=c("uid", "labPerson", "labPeriod", "Residual"))
                levels=c("uid", "Plate", "Residual"))

tmp$normalization <- factor(tmp$normalization, 
                  levels=c("MQ", "QN", "TA", "None"))

p<-ggplot(tmp,aes(y=sdcor,x=factor(grp),fill=factor(grp)))
print(
  p
  + geom_bar(stat="identity")
  +facet_wrap(~glycan+normalization, scales="free", ncol=4)
)
