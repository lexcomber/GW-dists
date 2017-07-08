######### Lex Comber
######### a.comber@leeds.ac.uk

library(GISTools)
library(foreign)
library(GWmodel)
library(car)
library(OpenStreetMap)
library(RgoogleMaps)
library(PBSmapping)
library(spgwr)
library(scales)
library(classInt)
library(tidyverse)
library(OSMscale)
library(reshape2)
library(repmis)

##### Part 1: Data prep
source_data("https://github.com/lexcomber/GW-dists/blob/master/datain.RData?raw=True")

# select variables (numeric)
index <- c(32,4,7,12,16,17,19,23,44,45)
summary(data.sp@data[,index])
data.sp <- data.sp[,index]
summary(data.sp@data)

# get rid of NAs
index <- as.vector(apply(data.sp@data, 1, function(x) sum(is.na(x))) == 0)
summary(index)
summary(data.sp@data[index,])
round(cor(data.sp@data[index,]), 2)
data.sp <- data.sp[index,]
as_tibble(data.sp@data)
# adjust distance matrices
dMat.ndp2 <- dMat.ndp2[index,index]
dim(dMat.ndp2)
dim(data.sp)

#### Figure 1. Map the data
png(filename = "F1.png", w = 5.5, h = 6, units = "in", res = 150)
par(mfrow = c(1,1))
plot(MyMap, removeMargin=T) 
plot(spTransform(data.sp, osm()),add = T, col="#25252580", pch = 19)
scaleBar(MyMap, x = 0.65, y = 0.035)
dev.off()

##### Example Distance metrics
dMat.ed <- gw.dist(coordinates(data.sp))
p.vals <- seq(0.25, 3, 0.25)
res.mat <- matrix(nrow = nrow(data.sp), ncol = 0) 
res.mat <- cbind(res.mat, dMat.ndp2[1,])
for (i in 1:length(p.vals)) {
  dMat.i <- gw.dist(coordinates(data.sp), p = p.vals[i])
  res.mat <- cbind(res.mat, dMat.i[,1])
  cat(i)
}  
colnames(res.mat) <- c("ND", as.character(c(p.vals)))
index <- order(res.mat[,1])
df <- as.data.frame(res.mat[index,])
df <- data.frame(x = 1:nrow(data.sp), df)
### plot functions
plot.func <- function(i){
  i <- i +2
  tit <- gsub("X", "",colnames(df)[i])
  df.tmp <- df[, c(1,2,i)]
  colnames(df.tmp)[3] <- "Y"
  p <- ggplot() +
    geom_point(data = df.tmp, aes(x,ND), col = "#25252580", 
               pch = 1, cex = 0.7) +
    geom_point(data = df.tmp, aes(x,Y), 
               col = add.alpha(brewer.pal(12, "Paired"), 0.7)[i-2], 
               pch = 19, cex = 0.7) +
    labs(subtitle = tit) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          legend.position = "none") 
  return(p)}
## multiplot function
multiplot2 <- function(plot.list, file, cols=3, layout=NULL) {
  library(grid)
  numPlots = length(plot.list)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plot.list[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
    }
  }
}
# now plot distances
for (i in 1:length(p.vals)) {
  p.i <- plot.func(i)
  tit <- paste("p", i, sep = "")
  assign(tit, p.i)
} 
png(filename = "F2.png", w = 8, h = 6, units = "in", res = 150)
multiplot2(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),cols = 4)
dev.off()

##### Part 2: Exploratory investigation: collinearity and regression
tab <- round(cor(data.sp@data[,-1]), 2)
write.csv(tab, file = "tab1.csv")

# set up a regression model with ALL variables
terms <- names(data.sp@data)
reg.mod <- paste(terms[1], "~")
for (i in 2:(length(terms))) {
  if (i != length(terms)) reg.mod <- paste(reg.mod, terms[i], "+")
  if (i == length(terms)) reg.mod <- paste(reg.mod, terms[i])
}
reg.mod <- as.formula(reg.mod)
# use for OLS regression
lm.1 <- lm(reg.mod, data = data.sp@data)
tab2 <- (summary(lm.1)[4])
tab2 <- round(as.matrix(tab2[[1]]), 3)
write.csv(tab2, file = "tab2.csv")

##### Table 1
## Global VIF and CN
tab <- summary(lm.1)[4]
tab <- round(tab[[1]], 3)
# VIFs > 10, VDP > 0.5 suggestive of collinearity
vif(lm.1) # all low
tab <- cbind(tab, append(NA, vif(lm.1)))
colnames(tab)[5] <- "VIF"
tab <- round(tab, 3)
tab
# setwd("~/Desktop/my_docs_mac/leeds_work/research/viet/Rplay/")
# write.csv(tab, "Table1.csv")
##### Condition number of the design matrix
X <- as.matrix(cbind(1,data.sp@data[2:10]))
BKWcn <- function(X) {
  p <- dim(X)[2]
  Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  Xsvd <- svd(Xscale)$d
  Xsvd[1] / Xsvd[p]
}
BKWcn(X)
# The BKW condition number (Belsley et al. 1980), is < 30
# indicative that collinearity is not a global problem for this data. 

##### Part 3: Exploratory GWR collinearity 
# not in the paper but helps to understand Part 4 and Part 5
# Euclidean: p = 2
# Manhattan: p = 1
# Minkowsski: p = 3
dMat.ed <- gw.dist(coordinates(data.sp))
# GWR euclidean and network distance model calibration
bw.ed <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV", dMat = dMat.ed) 
dMat.mh <- gw.dist(coordinates(data.sp), p = 1)
bw.mh <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV", dMat = dMat.mh) 
dMat.mk <- gw.dist(coordinates(data.sp), p = 3)
bw.mk <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV", dMat = dMat.mk) 
bw.nd <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
  adaptive = TRUE, approach = "CV", dMat = dMat.ndp2) 

bw.ed/nrow(data.sp)
bw.mh/nrow(data.sp)
bw.mk/nrow(data.sp)
bw.nd/nrow(data.sp)

# GW collinearity diagnostic
gw.col.ed <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.ed, 
  kernel = "bisquare", adaptive = TRUE, dMat = dMat.ed)
gw.col.mh <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.mh, 
  kernel = "bisquare", adaptive = TRUE, dMat = dMat.mh)
gw.col.mk <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.mk, 
  kernel = "bisquare", adaptive = TRUE, dMat = dMat.mk)
gw.col.nd <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.nd, 
  kernel = "bisquare", adaptive = TRUE, dMat = dMat.ndp2)

summary(gw.col.ed$SDF@data[,1:10])
summary(gw.col.mh$SDF@data[,1:10])
summary(gw.col.mk$SDF@data[,1:10])
summary(gw.col.nd$SDF@data[,1:10])

sum(gw.col.ed$SDF@data[,10] <30)/nrow(data.sp)
sum(gw.col.mh$SDF@data[,10] <30)/nrow(data.sp)
sum(gw.col.mk$SDF@data[,10] <30)/nrow(data.sp)
sum(gw.col.nd$SDF@data[,10] <30)/nrow(data.sp)

##### Part 4: Full GWR collinearity - varying p, theta = 0
# explore Mink distances
res.tab <- vector()
p.vals <- seq(0.05, 4, 0.05)
for (i in 1:length(p.vals)) {
  dMat.i <- gw.dist(coordinates(data.sp), p = p.vals[i])
  bw.i <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
    adaptive = TRUE, approach = "CV", dMat = dMat.i)
  gw.col.i <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.i, 
    kernel = "bisquare", adaptive = TRUE, dMat = dMat.i)
  res.i <- append(as.vector(summary(gw.col.i$SDF@data[,10])), 
    sum(gw.col.i$SDF@data[,10] < 30)/nrow(data.sp))
  res.i <- append(bw.i, res.i)
  res.tab <- append(res.tab, res.i)
}  
mat.pvals <- round(matrix(res.tab, ncol = 8, byrow = T), 2)  
rownames(mat.pvals) <- p.vals
colnames(mat.pvals) <- c("bw", names(summary(gw.col.i$SDF@data[,10])), "CN<30 rate")

df <- as.data.frame(mat.pvals)
df <- data.frame(df, x = as.numeric(as.character(rownames(mat.pvals))))
png(filename = "F3a.png", w = 8, h = 4, units = "in", res = 150)
ggplot(data = df, aes(x,Median)) +  
  geom_line(lwd =1.1)+
  xlab("Minkowski power value") +
  ylab("Median CN") +
  geom_ribbon(data=df,aes(ymin=Min.,ymax=Max.),fill = "indianred", alpha=0.4) +
  geom_ribbon(data=df,aes(ymin=X1st.Qu.,ymax=X3rd.Qu.),alpha=0.6) +
  ylim(c(0, max(df$Max.)))
dev.off()

png(filename = "F3b.png", w = 8, h = 4, units = "in", res = 150)
ggplot(data = df, aes(x,CN.30.rate)) +  
  geom_line(col = "indianred") +
  geom_point(data = df, aes(x,CN.30.rate), pch = 1, cex = 2) +
  xlab("Minkowski power value") +
  ylab("Proportion of CN < 30") 
dev.off()

##### Part 5: Full GWR collinearity - varying p, varying theta
# explore Mink distances
p.vals <- seq(0.05, 3, 0.05)
t.vals <- (seq(0,90,10) * (pi/180))
##### this will take a few minutes! #####
for (j in 1:length(t.vals)) {
  res.tab <- vector()
  for (i in 1:length(p.vals)) { 
    dMat.i <- gw.dist(coordinates(data.sp), p = p.vals[i], theta = t.vals[j])
    bw.i <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
                 adaptive = TRUE, approach = "CV", dMat = dMat.i)
    gw.col.i <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.i, 
                                kernel = "bisquare", adaptive = TRUE, dMat = dMat.i)
    res.i <- append(as.vector(summary(gw.col.i$SDF@data[,10])), 
                  sum(gw.col.i$SDF@data[,10] <30)/nrow(data.sp))
    res.i <- append(bw.i, res.i)
    res.tab <- append(res.tab, res.i)
  }
  tit <- paste("r", j, sep="")
  assign(tit, res.tab)
  cat("\n", j, "\n")
}  
#### Plot by Rotation degree - Figure 4
plot.func2 <- function(res.tab, t.val.i){
  mat <- round(matrix(res.tab, ncol = 8, byrow = T), 3)  
  rownames(mat) <- p.vals
  colnames(mat) <- c("bw", names(summary(gw.col.i$SDF@data[,10])), "CN<30 rate")
  df <- as.data.frame(mat)
  df <- data.frame(df, x = as.numeric(as.character(rownames(mat))))
  p.a<-ggplot(data = df, aes(x,Median)) +  
    geom_line(lwd =1.1)+
    #xlab("Minkowski power value") +
    labs(title=substitute(paste("",tit,degree),
      list(tit=as.character(t.val.i)))) +
    #ylab("Median CN") +
    geom_ribbon(data=df,aes(ymin=Min.,ymax=Max.),fill = "indianred", alpha=0.4) +
    geom_ribbon(data=df,aes(ymin=X1st.Qu.,ymax=X3rd.Qu.),alpha=0.6) +
    ylim(c(0, max(df$Max.))) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) 
  p.b<-ggplot(data = df, aes(x,CN.30.rate)) +  
    geom_line(col = "indianred") +
    geom_point(data = df, aes(x,CN.30.rate), pch = 1, cex = 1) +
    #xlab("Minkowski power value") +
    #ylab("Prop CN < 30") 
    ylim(c(0, 1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)
  return(list(p.a,p.b))
}
t.vals <- seq(0,90,10) 
for (i in 1:length(t.vals)){
  t.val.i <- t.vals[i]
  tit <- paste("r",i,sep = "")
  res.tab.i <- get(tit)
  tit <- paste("F4.",t.val.i,".png", sep = "")
  p.i <- plot.func2(res.tab.i, t.val.i)
  png(filename = tit, w = 5, h = 3, units = "in", res = 150)
  multiplot2(p.i, cols=1)
  dev.off()
}

#### Now plot by Power Figure 5
plot.func3 <- function(mat.i, p.val.k){
  df <- as.data.frame(mat.i)
  df <- data.frame(df, x = as.numeric(as.character(rownames(mat.i))))
  p.a<-ggplot(data = df, aes(x,Median)) +  
    geom_line(lwd =1.1)+
    xlab(substitute(paste("Rotation (",degree,")" ))) +
    labs(title = paste("p = ", p.val.k, sep = "")) +
    geom_ribbon(data=df,aes(ymin=Min.,ymax=Max.),fill = "indianred", alpha=0.4) +
    geom_ribbon(data=df,aes(ymin=X1st.Qu.,ymax=X3rd.Qu.),alpha=0.6) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())
  p.b<-ggplot(data = df, aes(x,CN.30.rate)) +  
    geom_line(col = "indianred") +
    geom_point(data = df, aes(x,CN.30.rate), pch = 1, cex = 1) +
    #xlab("Minkowski power value") +
    #ylab("Prop CN < 30") 
    ylim(c(0, 1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),)
  return(list(p.a,p.b))
}

p.vals.s <- seq(0.25, 3, 0.25)
p.vals.l <- seq(0.05, 3, 0.05)
for (k in 1: length(p.vals.s)) {
  p.val.k <- p.vals.s[k] 
  res.vec.i <- vector()
  for (i in 1:length(t.vals)){
    t.val.i <- t.vals[i]
    tit <- paste("r",i,sep = "")
    res.tab.i <- get(tit)
    mat.i <- round(matrix(res.tab.i, ncol = 8, byrow = T), 3)  
    rownames(mat.i) <- p.vals.l
    match.i <- match(p.val.k, rownames(mat.i))
    mat.i <- mat.i[match.i,]
    res.vec.i <- append(res.vec.i, mat.i)
  }
  mat.i <- round(matrix(res.vec.i, ncol = 8, byrow = T), 3)
  colnames(mat.i) <- c("bw", names(summary(gw.col.i$SDF@data[,10])), "CN<30 rate")
  rownames(mat.i) <- t.vals
  tit <- paste("F5.",p.val.k,".png", sep = "")
  p.i <- plot.func3(mat.i, p.val.k)
  png(filename = tit, w = 5, h = 3, units = "in", res = 150)
  multiplot2(p.i, cols=1)
  dev.off()
}

## Generate surfaces of CN and bandwidth (t.vals * p.vals) 
## Figures 6 and 7
t.vals <- seq(0,90,10) 
mat.CN <- matrix(ncol = 0, nrow = length(p.vals.l))
mat.bw <- matrix(ncol = 0, nrow = length(p.vals.l))
 for (i in 1:length(t.vals)){
  t.val.i <- t.vals[i]
  tit <- paste("r",i,sep = "")
  res.tab.i <- get(tit)
  mat.i <- round(matrix(res.tab.i, ncol = 8, byrow = T), 3)  
  rownames(mat.i) <- p.vals.l
  mat.CN <- cbind(mat.CN, mat.i[, 8])
  mat.bw <- cbind(mat.bw, mat.i[, 1])
 }  
rownames(mat.CN) <- p.vals.l
colnames(mat.CN) <- t.vals
rownames(mat.bw) <- p.vals.l
colnames(mat.bw) <- t.vals
image(mat.CN)
image(mat.bw)

tmp <- as.data.frame(mat.CN)
tmp$pv <- p.vals
aic_melt <- melt(tmp, id = "pv")
names(aic_melt) <- c("Power", "Rotation", "value")
aic_melt <- as.data.frame(apply(aic_melt, 2, as.double))
#### Figure 6
png(filename = "F6.png", w = 8, h = 4, units = "in", res = 150)
ggplot(data = na.omit(aic_melt)) +
  aes(x=Power, y=Rotation, fill=value, z = value) +
  geom_tile() +
  theme_bw() +
  #scale_fill_distiller(palette="Reds", direction = 1, na.value="white", name = "CN<30prop") +
  scale_fill_gradient(low = "white", high = "#A50F15", space = "Lab", name="") +  
  geom_contour(color = "black", alpha=0.5) +
 #            breaks = c(0.5, 0.1, 0.12, 0.15, 0.2, 0.3, 0.4), show.legend = F) +
  scale_x_continuous(breaks= seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(breaks = t.vals)
dev.off()

tmp <- as.data.frame(mat.bw)
tmp$bw <- p.vals
aic_melt <- melt(tmp, id = "bw")
names(aic_melt) <- c("Power", "Rotation", "value")
aic_melt <- as.data.frame(apply(aic_melt, 2, as.double))
#### Figure 7
png(filename = "F7.png", w = 8, h = 4, units = "in", res = 150)
ggplot(data = na.omit(aic_melt)) +
  aes(x=Power, y=Rotation, fill=value, z = value) +
  geom_tile() +
  theme_bw() +
  #scale_fill_distiller(palette="Reds", direction = 1, na.value="white", name = "CN<30prop") +
  scale_fill_gradient(low = "white", high = "#08519C", space = "Lab", name="") +  
  geom_contour(color = "black", alpha=0.5, show.legend = F) +
  scale_x_continuous(breaks= seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(breaks = t.vals)
dev.off()

### determine the MIN and MAX
# 1. MIN - most CN (prop of data points < 30)
# row index
i <- unlist(apply(mat.CN, 2, function(x) which(x == min(mat.CN))))
# col index
j <- unlist(apply(mat.CN, 1, function(x) which(x == min(mat.CN))))
# Min value
min(mat.CN) * 100
# Mink rotation
colnames(mat.CN)[j]
as.numeric(colnames(mat.CN)[j]) * (pi/180)
# Mink Power
rownames(mat.CN)[i]
# BW
mat.bw[i,j]
# 2. MAX - least CN (prop of data points < 30)
i <- unlist(apply(mat.CN, 2, function(x) which(x == max(mat.CN))))
# col index
j <- unlist(apply(mat.CN, 1, function(x) which(x == max(mat.CN))))
# Mink rotation
colnames(mat.CN)[j]
as.numeric(colnames(mat.CN)[j]) * (pi/180)
# Mink Power
rownames(mat.CN)[i]
# BW
mat.bw[i,j]

#### Figure 8
lm_eqn <- function(df){
  y <- df[,1]
  x <- df[,2]
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 3), 
                        b = format(coef(m)[2], digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
### Join 2 melts togther
tmp <- as.data.frame(mat.CN)
tmp$pv <- p.vals
aic_melt <- melt(tmp, id = "pv")
names(aic_melt) <- c("Power", "Rotation", "value")
tmp <- as.data.frame(mat.bw)
tmp$pv <- p.vals
aic_melt2 <- melt(tmp, id = "pv")
names(aic_melt2) <- c("Power", "Rotation", "value")
df <- data.frame(aic_melt, aic_melt2)
df$value.1 <- df$value.1 / nrow(data.sp)

head(df) 
df[which.min(df$value),] 
30 * (pi/180)

myColors <- add.alpha(rev(brewer.pal(10,"Spectral")), 0.5)
names(myColors) <- levels(df$Rotation)
colScale <- scale_colour_manual(name = "Rotation",values = myColors)
#df$Power <- rescale(df$Power, c(0.5, 3))
png(filename = "F8.png", w = 8, h = 5, units = "in", res = 300)
ggplot(data = df, aes(x = value.1, y = value, colour = Rotation, size = Power)) +
  geom_point(pch = 19) +
  colScale +  
  geom_text(y = 0.85, x = 0.15, label = lm_eqn(df[, c(3, 6)]), 
            parse = TRUE, size = 6, col = "black") +
  scale_x_continuous(name = "Adaptive Bandwidth (proportion of data points)") + 
  scale_y_continuous(name = "Proportion of local CNs < 30") +
  scale_fill_gradient(low = "white", high = "#08519C")   
dev.off()  

### Figure 9
# least to worst collinear
# best
d.mat.tmp <- gw.dist(coordinates(data.sp), p = 0.05, theta = (70 * pi / 180))
bw.i = 176
#bw.tmp <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
#	adaptive = TRUE, approach = "CV", dMat = d.mat.tmp)
gw.col.i <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.i, 
  kernel = "bisquare", adaptive = TRUE, dMat = d.mat.tmp)
summary(gw.col.i$SDF@data[,1:10])
(1- sum(gw.col.i$SDF@data[,10] < 30)/nrow(data.sp))*100

## 73.7%
i <- unlist(apply(mat.CN, 2, function(x) which(x == 0.737   )))
# col index
j <- unlist(apply(mat.CN, 1, function(x) which(x == 0.737   )))
# Mink rotation
colnames(mat.CN)[j]
as.numeric(colnames(mat.CN)[j]) * (pi/180)
# Mink Power
rownames(mat.CN)[i]
# BW
mat.bw[i,j]
d.mat.tmp <- gw.dist(coordinates(data.sp), p = 0.1, theta = (80 * pi / 180))
bw.ii = 101
#bw.tmp <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
#	adaptive = TRUE, approach = "CV", dMat = d.mat.tmp)
gw.col.ii <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.ii, 
  kernel = "bisquare", adaptive = TRUE, dMat = d.mat.tmp)
summary(gw.col.ii$SDF@data[,1:10])
(1- sum(gw.col.ii$SDF@data[,10] < 30)/nrow(data.sp))*100

## 12.5%
summary(as.vector(mat.CN))
# 1st Quartile ~12.5%
i <- unlist(apply(mat.CN, 2, function(x) which(x == summary(as.vector(mat.CN))[2]   )))
# col index
j <- unlist(apply(mat.CN, 1, function(x) which(x == summary(as.vector(mat.CN))[2]   )))
# Mink rotation
colnames(mat.CN)[j]
as.numeric(colnames(mat.CN)[j]) * (pi/180)
# Mink Power
rownames(mat.CN)[i]
# BW
mat.bw[i,j]
d.mat.tmp <- gw.dist(coordinates(data.sp), p = 1.75, theta = (30 * pi / 180))
bw.iii = 66
#bw.tmp <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
#	adaptive = TRUE, approach = "CV", dMat = d.mat.tmp)
gw.col.iii <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.iii, 
  kernel = "bisquare", adaptive = TRUE, dMat = d.mat.tmp)
summary(gw.col.iii$SDF@data[,1:10])
(1- sum(gw.col.iii$SDF@data[,10] < 30)/nrow(data.sp))*100

# worst
d.mat.tmp <- gw.dist(coordinates(data.sp), p = 0.7, theta = (30 * pi / 180))
bw.iv = 32  
#bw.tmp <- bw.gwr(reg.mod, data = data.sp, kernel = "bisquare", 
#	adaptive = TRUE, approach = "CV", dMat = d.mat.tmp)
gw.col.iv <- gwr.collin.diagno(reg.mod, data = data.sp, bw = bw.iv, 
  kernel = "bisquare", adaptive = TRUE, dMat = d.mat.tmp)
summary(gw.col.iv$SDF@data[,1:10])
(1- sum(gw.col.iv$SDF@data[,10] < 30)/nrow(data.sp))*100


map.cn <- function(gw.col = gw.col.iv) {
	index <- gw.col$SDF@data[,10] > 30
	summary(index)
	plot(MyMap, removeMargin=T) 
	plot(spTransform(data.sp[index,], osm()),add = T, col="#25252580", pch = 19)
}

png(filename = "F9a.png", w = 5.5, h = 6, units = "in", res = 150)
map.cn(gw.col.i)
dev.off()
png(filename = "F9b.png", w = 5.5, h = 6, units = "in", res = 150)
map.cn(gw.col.ii)
dev.off()
png(filename = "F9c.png", w = 5.5, h = 6, units = "in", res = 150)
map.cn(gw.col.iii)
dev.off()
png(filename = "F9d.png", w = 5.5, h = 6, units = "in", res = 150)
map.cn(gw.col.iv)
dev.off()


## END






