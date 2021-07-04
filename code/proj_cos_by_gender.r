# load data
obs <- read.csv("../data/adjmatrix-sample.csv",header=T)
obs[is.na(obs)] <- 0
users <- read.csv("../data/nodeattr-sample.csv",header=T)
obs.ext <- merge(obs,users,by="caseid")

# MALE case
obs.male <- obs.ext[which(obs.ext$female==0),2:12]  
names(obs.male) <- gsub("use","",names(obs.male))

# build male network
require(igraph)
g.m <- graph.incidence(obs.male,mode="all")
summary(g.m)

# degree distributions
deg_users.m<-degree(g.m)[V(g.m)$type==FALSE]
png('../output/user_degree_male.png',width=500,height=300)
hist(deg_users.m,col="#1f78b4",xlab="Number of SNS per male user",ylab="Frequency",main="")
dev.off()
pdf('../output/user_degree_male.pdf',width=5,height=3)
hist(deg_users.m,col="#1f78b4",xlab="Number of SNS per male user",ylab="Frequency",main="")
dev.off()
deg_sns.m<-degree(g.m)[V(g.m)$type==TRUE]
png('../output/sns_degree_male.png',width=500,height=300)
hist(deg_sns.m,col="#1a9641",xlab="Number of male users per SNS",ylab="Frequency",main="",breaks=14)
dev.off()
pdf('../output/sns_degree_male.pdf',width=5,height=3)
hist(deg_sns.m,col="#1a9641",xlab="Number of male users per SNS",ylab="Frequency",main="",breaks=14)
dev.off()

# compute common neighbors
m <- as.matrix(obs.male)
m_cn <- t(m) %*% m
for (i in 1:nrow(m_cn)) {
  m_cn[i,i] <- 0
}

# compute cosine
degrees <- degree(g.m, v=V(g.m)[V(g.m)$type==TRUE])
degpairs <- degrees %*% t(degrees)
degpairs[degpairs==0] <- 1
m_cos <- m_cn/sqrt(degpairs)

# build projection graph
g_proj.m <- graph.adjacency(m_cos, mode="undirected", weighted=TRUE)
summary(g_proj.m)

# map weights to colors
rbPal <- colorRampPalette(c('gray20','black'))
E(g_proj.m)$color <- rbPal(10)[as.numeric(cut(E(g_proj.m)$weight,breaks=10))]
E(g_proj.m)$thickness <- as.numeric(cut(E(g_proj.m)$weight,breaks=10))

# export projection graph
write.graph(g_proj.m,"../output/sns-associations-male.graphml",format="graphml")

# export projections at different thresholds
hist(E(g_proj.m)$weight)
thr <- c(.3,.4,.5,.6)
for (t in thr) {
  thresholded_graph <- delete.edges(g_proj.m, E(g_proj.m)[E(g_proj.m)$weight<t])
  thresholded_graph <- delete.vertices(thresholded_graph, V(thresholded_graph)[degree(thresholded_graph)==0])
  summary(thresholded_graph)

  write.graph(thresholded_graph,paste("../output/sns-associations-male_thr",t,".graphml",sep=""),format="graphml")
}

#E(g_proj.m)[order(-E(g_proj.m)$weight)]
top.m <- as.data.frame(cbind(get.edgelist(g_proj.m),round(E(g_proj.m)$weight,2)))
names(top.m) <- c("SNS1","SNS2","cosine.association")
top.m <- top.m[order(top.m$cosine.association,decreasing=T),]
revtop.m <- top.m[order(top.m$cosine.association,decreasing=F),]
revtop.m[!(revtop.m$SNS1%in%c("myspace","flickr","googleplus")|revtop.m$SNS2%in%c("myspace","flickr","googleplus")),]

# export top network
g2 <- delete.edges(g_proj.m,E(g_proj.m)[!weight%in%tail(sort(E(g_proj.m)$weight),10)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)

# export projection graph
write.graph(g2,"../output/sns-associations-male-top10.graphml",format="graphml")

# store data for legend
d_legend <- head(top.m,10)

top10 <- head(top.m,10)
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
pdf('../output/legend-sns-associations-male-top10.pdf',width=1,height=3)
grid.arrange(legend,ncol=1,nrow=1,widths=c(1/6))
dev.off()


# FEMALE case
obs.female <- obs.ext[which(obs.ext$female==1),2:12]  
names(obs.female) <- gsub("use","",names(obs.female))

# build male network
require(igraph)
g.f <- graph.incidence(obs.female,mode="all")
summary(g.f)

# degree distributions
deg_users.f<-degree(g.f)[V(g.f)$type==FALSE]
png('../output/user_degree_female.png',width=500,height=300)
hist(deg_users.f,col="#1f78b4",xlab="Number of SNS per female user",ylab="Frequency",main="")
dev.off()
pdf('../output/user_degree_female.pdf',width=5,height=3)
hist(deg_users.f,col="#1f78b4",xlab="Number of SNS per female user",ylab="Frequency",main="")
dev.off()
deg_sns.f<-degree(g.f)[V(g.f)$type==TRUE]
png('../output/sns_degree_female.png',width=500,height=300)
hist(deg_sns.f,col="#1a9641",xlab="Number of female users per SNS",ylab="Frequency",main="",breaks=14)
dev.off()
pdf('../output/sns_degree_female.pdf',width=5,height=3)
hist(deg_sns.f,col="#1a9641",xlab="Number of female users per SNS",ylab="Frequency",main="",breaks=14)
dev.off()

# compute common neighbors
m <- as.matrix(obs.female)
m_cn <- t(m) %*% m
for (i in 1:nrow(m_cn)) {
  m_cn[i,i] <- 0
}

# compute cosine
degrees <- degree(g.f, v=V(g.f)[V(g.f)$type==TRUE])
degpairs <- degrees %*% t(degrees)
degpairs[degpairs==0] <- 1
m_cos <- m_cn/sqrt(degpairs)

# build projection graph
g_proj.f <- graph.adjacency(m_cos, mode="undirected", weighted=TRUE)
summary(g_proj.f)

# map weights to colors
rbPal <- colorRampPalette(c('gray20','black'))
E(g_proj.f)$color <- rbPal(10)[as.numeric(cut(E(g_proj.f)$weight,breaks=10))]
E(g_proj.f)$thickness <- as.numeric(cut(E(g_proj.f)$weight,breaks=10))

# export projection graph
write.graph(g_proj.f,"../output/sns-associations-female.graphml",format="graphml")

# export projections at different thresholds
hist(E(g_proj.f)$weight)
thr <- c(.3,.4,.5,.6)
for (t in thr) {
  thresholded_graph <- delete.edges(g_proj.f, E(g_proj.f)[E(g_proj.f)$weight<t])
  thresholded_graph <- delete.vertices(thresholded_graph, V(thresholded_graph)[degree(thresholded_graph)==0])
  summary(thresholded_graph)
  
  write.graph(thresholded_graph,paste("../output/sns-associations-female_thr",t,".graphml",sep=""),format="graphml")
}

#E(g_proj.f) <- E(g_proj.f)[order(-E(g_proj.f)$weight)]
top.f <- as.data.frame(cbind(get.edgelist(g_proj.f),round(E(g_proj.f)$weight,2)))
names(top.f) <- c("SNS1","SNS2","cosine.association")
top.f <- top.f[order(top.f$cosine.association,decreasing=T),]
revtop.f <- top.f[order(top.f$cosine.association,decreasing=F),]
revtop.f[!(revtop.f$SNS1%in%c("myspace","flickr","googleplus")|revtop.f$SNS2%in%c("myspace","flickr","googleplus")),]

#plot(g_proj,layout=layout.auto(g_proj),edge.width=E(g_proj)$weight)

g2 <- delete.edges(g_proj.f,E(g_proj.f)[!weight%in%tail(sort(E(g_proj.f)$weight),10)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)

# export projection graph
write.graph(g2,"../output/sns-associations-female-top10.graphml",format="graphml")

# store data for legend
d_legend <- rbind(d_legend,head(top.f,10))
write.csv(d_legend,"../output/links-gender-for-legend.csv")


top10 <- head(top.f,10)
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
pdf('../output/legend-sns-associations-female-top10.pdf',width=1,height=3)
grid.arrange(legend,ncol=1,nrow=1,widths=c(1/6))
dev.off()


t <- merge(top.f,top.m,by=c("SNS1","SNS2"))
t$cosine.association.x <- as.numeric(as.character(t$cosine.association.x))
t$cosine.association.y <- as.numeric(as.character(t$cosine.association.y))
library(ggplot2)
d <- t[which(abs(t$cosine.association.x-t$cosine.association.y)>0.13),]
p <- ggplot(t,aes(x=t$cosine.association.x,y=t$cosine.association.y)) + geom_point() +
  geom_abline(intercept=0,slope=1,linetype="dotted") + theme_bw() +
  xlab("Associations for females") + ylab("Associations for males") + xlim(0,1) + ylim(0,1) +
  annotate("text",x=d$cosine.association.x,y=d$cosine.association.y,
           label=paste(d$SNS1,"-",d$SNS2,sep=),size=3)
ggsave("../output/genderdiff.pdf",width=5,height=5)

