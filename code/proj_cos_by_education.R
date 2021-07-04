# load data
obs <- read.csv("../data/adjmatrix-sample.csv",header=T)
obs[is.na(obs)] <- 0
users <- read.csv("../data/nodeattr-sample.csv",header=T)
obs.ext <- merge(obs,users,by="caseid")

# eduhsorless
obs.1 <- obs.ext[which(obs.ext$eduhsorless==1),2:12]  
names(obs.1) <- gsub("use","",names(obs.1))

# build 1 network
require(igraph)
g.1 <- graph.incidence(obs.1,mode="all")
summary(g.1)

# degree distributions
deg_users.1<-degree(g.1)[V(g.1)$type==FALSE]
png('../output/user_degree_edu1.png',width=500,height=300)
hist(deg_users.1,col="#1f78b4",xlab="Number of SNS per user with highschool degree or less",ylab="Frequency",main="")
dev.off()
pdf('../output/user_degree_edu1.pdf',width=5,height=3)
hist(deg_users.1,col="#1f78b4",xlab="Number of SNS per user with highschool degree or less",ylab="Frequency",main="")
dev.off()
deg_sns.1<-degree(g.1)[V(g.1)$type==TRUE]
png('../output/sns_degree_edu1.png',width=500,height=300)
hist(deg_sns.1,col="#1a9641",xlab="Number of users with highschool degree or less per SNS",ylab="Frequency",main="",breaks=14)
dev.off()
pdf('../output/sns_degree_edu1.pdf',width=5,height=3)
hist(deg_sns.1,col="#1a9641",xlab="Number of users with highschool degree or less per SNS",ylab="Frequency",main="",breaks=14)
dev.off()

# compute common neighbors
m <- as.matrix(obs.1)
m_cn <- t(m) %*% m
for (i in 1:nrow(m_cn)) {
  m_cn[i,i] <- 0
}

# compute cosine
degrees <- degree(g.1, v=V(g.1)[V(g.1)$type==TRUE])
degpairs <- degrees %*% t(degrees)
degpairs[degpairs==0] <- 1
m_cos <- m_cn/sqrt(degpairs)

# build projection graph
g_proj.1 <- graph.adjacency(m_cos, mode="undirected", weighted=TRUE)
summary(g_proj.1)

# map weights to colors
rbPal <- colorRampPalette(c('gray20','black'))
E(g_proj.1)$color <- rbPal(10)[as.numeric(cut(E(g_proj.1)$weight,breaks=10))]
E(g_proj.1)$thickness <- as.numeric(cut(E(g_proj.1)$weight,breaks=10))

# export projection graph
write.graph(g_proj.1,"../output/sns-associations-edu1.graphml",format="graphml")

top.1 <- as.data.frame(cbind(get.edgelist(g_proj.1),round(E(g_proj.1)$weight,2)))
names(top.1) <- c("SNS1","SNS2","cosine.association")
top.1 <- top.1[order(top.1$cosine.association,decreasing=T),]

# store data for legend
d_legend <- head(top.1,11)

# export top network
g2 <- delete.edges(g_proj.1,E(g_proj.1)[!weight%in%tail(sort(E(g_proj.1)$weight),11)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)

# export projection graph
write.graph(g2,"../output/sns-associations-edu1-top11.graphml",format="graphml")


top10 <- head(top.1,10)
library(gridExtra)
library(ggplot2)
rbPal <- colorRampPalette(c('gray80','gray0'))
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
top10$cosine.association <- as.numeric(as.character(top10$cosine.association))
top10$color <- rbPal(10)[as.numeric(cut(top10$cosine.association,breaks=10))]
top10 <- top10[order(-top10$cosine.association),]
p1 <- ggplot(top10,aes(SNS1,SNS2)) + geom_raster(aes(fill=as.factor(top10$cosine.association))) + 
  scale_fill_manual(values=rev(top10$color),breaks=top10$cosine.association,name="Cosine\nsimilarity")
legend <- g_legend(p1)
pdf('../output/legend-sns-associations-edu1-top10.pdf',width=1,height=3)
grid.arrange(legend,ncol=1,nrow=1,widths=c(1/6))
dev.off()


# edusc
obs.2 <- obs.ext[which(obs.ext$edusc==1),2:12]  
names(obs.2) <- gsub("use","",names(obs.2))

# build 2 network
require(igraph)
g.2 <- graph.incidence(obs.2,mode="all")
summary(g.2)

# degree distributions
deg_users.2<-degree(g.2)[V(g.2)$type==FALSE]
png('../output/user_degree_edu2.png',width=500,height=300)
hist(deg_users.2,col="#1f78b4",xlab="Number of SNS per user with some college",ylab="Frequency",main="")
dev.off()
pdf('../output/user_degree_edu2.pdf',width=5,height=3)
hist(deg_users.2,col="#1f78b4",xlab="Number of SNS per user with some college",ylab="Frequency",main="")
dev.off()
deg_sns.2<-degree(g.2)[V(g.2)$type==TRUE]
png('../output/sns_degree_edu2.png',width=500,height=300)
hist(deg_sns.2,col="#1a9641",xlab="Number of users with some college per SNS",ylab="Frequency",main="",breaks=14)
dev.off()
pdf('../output/sns_degree_edu2.pdf',width=5,height=3)
hist(deg_sns.2,col="#1a9641",xlab="Number of users with some college per SNS",ylab="Frequency",main="",breaks=14)
dev.off()

# compute common neighbors
m <- as.matrix(obs.2)
m_cn <- t(m) %*% m
for (i in 1:nrow(m_cn)) {
  m_cn[i,i] <- 0
}

# compute cosine
degrees <- degree(g.2, v=V(g.2)[V(g.2)$type==TRUE])
degpairs <- degrees %*% t(degrees)
degpairs[degpairs==0] <- 1
m_cos <- m_cn/sqrt(degpairs)

# build projection graph
g_proj.2 <- graph.adjacency(m_cos, mode="undirected", weighted=TRUE)
summary(g_proj.2)

# map weights to colors
rbPal <- colorRampPalette(c('gray20','black'))
E(g_proj.2)$color <- rbPal(10)[as.numeric(cut(E(g_proj.2)$weight,breaks=10))]
E(g_proj.2)$thickness <- as.numeric(cut(E(g_proj.2)$weight,breaks=10))

# export projection graph
write.graph(g_proj.2,"../output/sns-associations-edu2.graphml",format="graphml")

top.2 <- as.data.frame(cbind(get.edgelist(g_proj.2),round(E(g_proj.2)$weight,2)))
names(top.2) <- c("SNS1","SNS2","cosine.association")
top.2 <- top.2[order(top.2$cosine.association,decreasing=T),]

d_legend <- rbind(d_legend,head(top.2,10))

# export top network
g2 <- delete.edges(g_proj.2,E(g_proj.2)[!weight%in%tail(sort(E(g_proj.2)$weight),10)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)

# export projection graph
write.graph(g2,"../output/sns-associations-edu2-top10.graphml",format="graphml")

top10 <- head(top.2,10)
library(gridExtra)
library(ggplot2)
rbPal <- colorRampPalette(c('gray80','gray0'))
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
top10$cosine.association <- as.numeric(as.character(top10$cosine.association))
top10$color <- rbPal(10)[as.numeric(cut(top10$cosine.association,breaks=10))]
top10 <- top10[order(-top10$cosine.association),]
p1 <- ggplot(top10,aes(SNS1,SNS2)) + geom_raster(aes(fill=as.factor(top10$cosine.association))) + 
  scale_fill_manual(values=rev(top10$color),breaks=top10$cosine.association,name="Cosine\nsimilarity")
legend <- g_legend(p1)
pdf('../output/legend-sns-associations-edu2-top10.pdf',width=1,height=3)
grid.arrange(legend,ncol=1,nrow=1,widths=c(1/6))
dev.off()

# eduaormore
obs.3 <- obs.ext[which(obs.ext$edubaormore==1),2:12]  
names(obs.3) <- gsub("use","",names(obs.3))

# build 1 network
require(igraph)
g.3 <- graph.incidence(obs.3,mode="all")
summary(g.3)

# degree distributions
deg_users.3<-degree(g.3)[V(g.3)$type==FALSE]
png('../output/user_degree_edu3.png',width=500,height=300)
hist(deg_users.3,col="#1f78b4",xlab="Number of SNS per user with college or more",ylab="Frequency",main="")
dev.off()
pdf('../output/user_degree_edu3.pdf',width=5,height=3)
hist(deg_users.3,col="#1f78b4",xlab="Number of SNS per user with college or more",ylab="Frequency",main="")
dev.off()
deg_sns.3<-degree(g.3)[V(g.3)$type==TRUE]
png('../output/sns_degree_edu3.png',width=500,height=300)
hist(deg_sns.3,col="#1a9641",xlab="Number of users with college or more per SNS",ylab="Frequency",main="",breaks=14)
dev.off()
pdf('../output/sns_degree_edu3.pdf',width=5,height=3)
hist(deg_sns.3,col="#1a9641",xlab="Number of users with college or more per SNS",ylab="Frequency",main="",breaks=14)
dev.off()

# compute common neighbors
m <- as.matrix(obs.3)
m_cn <- t(m) %*% m
for (i in 1:nrow(m_cn)) {
  m_cn[i,i] <- 0
}

# compute cosine
degrees <- degree(g.3, v=V(g.3)[V(g.3)$type==TRUE])
degpairs <- degrees %*% t(degrees)
degpairs[degpairs==0] <- 1
m_cos <- m_cn/sqrt(degpairs)

# build projection graph
g_proj.3 <- graph.adjacency(m_cos, mode="undirected", weighted=TRUE)
summary(g_proj.3)

# map weights to colors
rbPal <- colorRampPalette(c('gray20','black'))
E(g_proj.3)$color <- rbPal(10)[as.numeric(cut(E(g_proj.3)$weight,breaks=10))]
E(g_proj.3)$thickness <- as.numeric(cut(E(g_proj.3)$weight,breaks=10))

# export projection graph
write.graph(g_proj.3,"../output/sns-associations-edu3.graphml",format="graphml")

top.3 <- as.data.frame(cbind(get.edgelist(g_proj.3),round(E(g_proj.3)$weight,2)))
names(top.3) <- c("SNS1","SNS2","cosine.association")
top.3 <- top.3[order(top.3$cosine.association,decreasing=T),]

d_legend <- rbind(d_legend,head(top.3,11))
write.csv(d_legend,"../output/links-education-for-legend.csv")

# export top network
g2 <- delete.edges(g_proj.3,E(g_proj.3)[!weight%in%tail(sort(E(g_proj.3)$weight),11)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)

# export projection graph
write.graph(g2,"../output/sns-associations-edu3-top11.graphml",format="graphml")

top10 <- head(top.3,10)
library(gridExtra)
library(ggplot2)
rbPal <- colorRampPalette(c('gray80','gray0'))
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
top10$cosine.association <- as.numeric(as.character(top10$cosine.association))
top10$color <- rbPal(10)[as.numeric(cut(top10$cosine.association,breaks=10))]
top10 <- top10[order(-top10$cosine.association),]
p1 <- ggplot(top10,aes(SNS1,SNS2)) + geom_raster(aes(fill=as.factor(top10$cosine.association))) + 
  scale_fill_manual(values=rev(top10$color),breaks=top10$cosine.association,name="Cosine\nsimilarity")
legend <- g_legend(p1)
pdf('../output/legend-sns-associations-edu3-top10.pdf',width=1,height=3)
grid.arrange(legend,ncol=1,nrow=1,widths=c(1/6))
dev.off()


t <- merge(top.1,top.2,by=c("SNS1","SNS2"))
t <- merge(t,top.3,by=c("SNS1","SNS2"))
t$cosine.association.x <- as.numeric(as.character(t$cosine.association.x))
t$cosine.association.y <- as.numeric(as.character(t$cosine.association.y))
t$cosine.association <- as.numeric(as.character(t$cosine.association))
library(ggplot2)
d <- t[which(abs(t$cosine.association.x-t$cosine.association.y)>0.1&
         abs(t$cosine.association.x-t$cosine.association)>0.1&
         abs(t$cosine.association-t$cosine.association.y)>0.1),]
tt <- rbind(top.1,top.2,top.3)
tt$type <- c(rep("Highschool or less",nrow(top.1)),
             rep("Some college",nrow(top.2)),
             rep("College or more", nrow(top.3)))
tt$pair <- paste(tt$SNS1,"-",tt$SNS2,sep="")
tt$cosine.association <- as.numeric(as.character(tt$cosine.association))
tt <- tt[which(tt$pair%in%c("facebook-linkedin","pinterest-instagram","pinterest-linkedin",
                            "myspace-googleplus","myspace-instagram",
                            "myspace-linkedin","myspace-twitter",
                            "tumblr-reddit")),]
tt$type <- factor(tt$type,levels=c("Highschool or less","Some college","College or more"))
p <- ggplot(tt,aes(x=tt$pair,y=tt$cosine.association)) + 
  geom_bar(stat="identity",aes(fill=type),position="dodge") +
  xlab("Pairs of SNS with the highest differences\namong the three education levels") + 
  ylab("Cosine similarity") + theme_bw() + 
  scale_fill_manual(values=c("red","grey","seagreen3"),name="Education level") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("../output/edudiff.pdf",width=8,height=5)

