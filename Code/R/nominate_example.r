library(wnominate)
library(pscl)
library(ggplot2)
library(dplyr)

rm(list=ls())
s118<-read.csv("datasets/S118_votes.csv")

votes <-s118[,c("icpsr", "rollnumber", "cast_code")]

table(votes$cast_code)
votes$cast_code <- ifelse(votes$cast_code==1,
                          1, 
                          ifelse(votes$cast_code==6,
                                 0,
                                 NA))

voteMat <- reshape(votes, 
                   direction="wide",
                   idvar="icpsr",
                   timevar="rollnumber")


V <- as.matrix(voteMat[,-1])
rownames(V) <- voteMat[,1]
dim(V)
dropLeg<-rowSums(V,na.rm=TRUE)<20
voteMat[dropLeg,]$icpsr #Andy Kim (NJ) and Adam Schiff (CA)
V <- V[!dropLeg,]

which(rownames(V)=="41301") # Warren
which(rownames(V)=="42102") # Tuberville
V <- rbind(V[67,,drop=FALSE],
           V[-c(67,88),], 
           V[88,,drop=FALSE])
dim(V)

keepVote <- which(colMeans(V,na.rm=T) > 0.025 & colMeans(V,na.rm=T) <.975 )
V <- V[,keepVote]
dim(V)



rc118 <- rollcall(V,legis.names=rownames(V))
NOMout <- wnominate(rc118, polarity = nrow(V),dims=1)
summary(NOMout)


ideal <- data.frame(icpsr=rownames(NOMout$legislators),
                    ideal=NOMout$legislators$coord1D)

info <-read.csv("datasets/S118_info.csv")

ideal <- merge(ideal, info, by="icpsr")
ideal <- ideal %>% arrange(ideal) 
head(ideal)
tail(ideal)

ideal <- ideal %>% 
  mutate(Party=ifelse(party_code ==200, "GOP", "DEM"),
         last_name= gsub(pattern = ",.*$", replacement = "", x = bioname),
         Rank = cumsum((ideal - lag(ideal, default=-1))>0)+1)

ggplot(ideal)+
  geom_text(aes(x=ideal, y=Rank,label = last_name,color=Party), 
            position="jitter", alpha=.4)+
  scale_color_discrete(palette=c("blue","red"))
ggplot(ideal)+
  geom_histogram(aes(x=ideal, 
                     y=after_stat(density),
                     fill=Party),color="white", alpha=.7,
                 bins=12)+
  # scale_color_discrete(palette=c("blue","red"))+
  scale_fill_discrete(palette=c("blue","red"))


