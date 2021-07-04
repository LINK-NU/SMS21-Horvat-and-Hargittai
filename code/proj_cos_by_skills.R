# load data
obs <- read.csv("../data/adjmatrix-sample.csv",header=T)
obs[is.na(obs)] <- 0
users <- read.csv("../data/nodeattr-sample.csv",header=T)
obs.ext <- merge(obs,users,by="caseid")
obs.ext$skills <- with(obs.ext, cut(skillsmean, 
                                breaks=quantile(skillsmean,probs=seq(0,1,by=0.24),na.rm=TRUE), 
                                include.lowest=TRUE))

# build network for lowest quartile
obs.1 <- obs.ext[which(obs.ext$skills=="[1,2.58]"),2:12]  
names(obs.1) <- gsub("use","",names(obs.1))
names(obs.1) <- c("Facebook","MySpace","Twitter","Tumblr","Pinterest","LinkedIn","GooglePlus",
                  "Instagram","Flickr","Reddit","Snapchat")

require(igraph)
g.1 <- graph.incidence(obs.1,mode="all")
summary(g.1)

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
write.graph(g_proj.1,"../output/sns-associations-lowskill.graphml",format="graphml")

top.1 <- as.data.frame(cbind(get.edgelist(g_proj.1),round(E(g_proj.1)$weight,2)))
names(top.1) <- c("SNS1","SNS2","cosine.association")
top.1 <- top.1[order(top.1$cosine.association,decreasing=T),]

# export top network
g2 <- delete.edges(g_proj.1,E(g_proj.1)[!weight%in%tail(sort(E(g_proj.1)$weight),10)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)
# export projection graph
write.graph(g2,"../output/sns-associations-lowskill-top10.graphml",format="graphml")

# store data for legend
d_legend <- head(top.1,10)


# build network for highest quartile
obs.2 <- obs.ext[which(obs.ext$skills=="(4.17,5]"),2:12]  
names(obs.2) <- gsub("use","",names(obs.2))
names(obs.2) <- c("Facebook","MySpace","Twitter","Tumblr","Pinterest","LinkedIn","GooglePlus",
                  "Instagram","Flickr","Reddit","Snapchat")
require(igraph)
g.2 <- graph.incidence(obs.2,mode="all")
summary(g.2)

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
write.graph(g_proj.2,"../output/sns-associations-highskill.graphml",format="graphml")

top.2 <- as.data.frame(cbind(get.edgelist(g_proj.2),round(E(g_proj.2)$weight,2)))
names(top.2) <- c("SNS1","SNS2","cosine.association")
top.2 <- top.2[order(top.2$cosine.association,decreasing=T),]

# export top network
g2 <- delete.edges(g_proj.2,E(g_proj.2)[!weight%in%tail(sort(E(g_proj.2)$weight),11)])
iso <- V(g2)[degree(g2)==0]
g2 <- delete.vertices(g2,iso)
summary(g2)
rbPal <- colorRampPalette(c('gray80','gray0'))
E(g2)$color <- rbPal(10)[as.numeric(cut(E(g2)$weight,breaks = 10))]
E(g2)$normalized <- round(1+(4/(max(E(g2)$weight)-min(E(g2)$weight))*(E(g2)$weight-min(E(g2)$weight))))
E(g2)$weight <- round(E(g2)$weight,2)
# export projection graph
write.graph(g2,"../output/sns-associations-highskill-top11.graphml",format="graphml")

# store data for legend
d_legend <- rbind(d_legend,head(top.2,11))
write.csv(d_legend,"../output/links-skills-for-legend.csv")
