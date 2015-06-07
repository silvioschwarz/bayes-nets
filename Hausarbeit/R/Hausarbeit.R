# initilization ####

# Clear the memory 
rm(list = ls())

# Set working directory
path = getwd()
setwd(path)

#style
op = par()
settings = options()
par(pch = 19, ann = T, cex = 1)

# read in additional functions
pathnames <- list.files(pattern="[.]R$", path=".\\functions", full.names=TRUE);
sapply(pathnames, FUN=source);

packages = c("bnlearn", "Rgraphviz")

loadPackages(packages)

# load and visualize data ####

# read data
data <- read.csv(file = "data.dat") 

# matrix plot of all variables

cairo_pdf(file = "..\\Figures\\allvsall.pdf",
          height = 9,
          width = 16)

plot(data, 
     main = "Scatterplot Matrix of all Variables", 
     col  = rgb(0, 0, 100, 50, maxColorValue = 255))

Sys.sleep(600)

dev.off()
#histograms
#par(mfrow=c(3,2))
#hist(data$SD)
#box()
#hist(data$Q0)
#box()
#hist(data$kappa)
#box()
#hist(data$Vs30)
#box()
#hist(data$MAG)
#box()
#hist(data$DIST)
#box()


#par(op)
#options(settings)
#hist(data$PGA)

#hist(data$PGA, c(min(data$PGA),-5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, max(data$PGA)))

# divide into learn and test data ####

# make copy of data
data.disc = data

# Daten diskretisieren (geht sicher sch?ner...)
data.disc$SD = as.factor(findInterval(data$SD, c(0, 0.8792, 5.438, 14.92, 58, 500)))
data.disc$Q0 = as.factor(findInterval(data$Q0, c(0, 330, 5000)))
data.disc$kappa = as.factor(findInterval(data$kappa, c(0, 0.01053, 0.0345, 0.1)))
data.disc$Vs30 = as.factor(findInterval(data$Vs30, c(600, 1704.5, 2800)))
data.disc$MAG = as.factor(findInterval(data$MAG, c(5, 6.271, 7.5)))
data.disc$DIST = as.factor(findInterval(data$DIST, c(1, 4.38, 15.4885, 55.84, 200)))
data.disc$PGA = as.factor(findInterval(data$PGA, c(-Inf, -5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, Inf)))

# make folds
set.seed(5555)

k = 10 #number of folds
index = sample(1:nrow(data))
folds = split(index, 1:k)


# set up the nets ####

PGA_causal = empty.graph(names(data))
PGA_naive = empty.graph(names(data))

# causal network
arcs(PGA_causal) = matrix(c("SD","MAG","MAG","PGA", "DIST","PGA", "Vs30", "PGA", "kappa","Vs30", "Q0", "Vs30" ),ncol = 2, byrow = TRUE, 
                       dimnames = list(c(),c("from","to")))

graphviz.plot(PGA_causal,layout = "neato",main = "Causal Network")


# naive Bayes
PGA_naive = naive.bayes(data.disc, "PGA",names(data.disc)[-7])

graphviz.plot(PGA_naive,layout = "fdp",main = "Naive Bayes Network")



# Learning nets ####

# allocate lists
fit.causal = vector("list", length(folds))
fit.naive = vector("list", length(folds))

for(i in 1:length(folds)){
  fit.causal[[i]] = bn.fit(PGA_causal, data.disc[-folds[[i]],], method = "bayes")
  fit.naive[[i]] = bn.fit(PGA_naive, data.disc[-folds[[i]],], method = "bayes")
}

# Sturktur gelernt

#allocate lists
PGA_hc = vector("list", length(folds))
PGA_gs = vector("list", length(folds))
fit.hc = vector("list", length(folds))
fit.gs = vector("list", length(folds))

for(i in 1:length(folds)){
  # Score based hill-climber
  PGA_hc[[i]] = hc(data.disc[-folds[[i]],]) 
  
  # constraint based grow shrink 
  PGA_gs[[i]] = gs(data.disc[-folds[[i]],]) 
  
  # Parameter lernen 
  fit.hc[[i]] = bn.fit(PGA_hc[[i]], data.disc[-folds[[i]],])
  fit.gs[[i]] = bn.fit(PGA_gs[[i]], data.disc[-folds[[i]],])
}


# Validation#### 

pga_interval = c(min(data$PGA), -5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, max(data$PGA))

midpoints = 1:8
for(i in 1:8){
  midpoints[i] = mean(c(pga_interval[i], pga_interval[i+1]))
}


net.list = list(fit.causal, fit.naive, fit.hc, fit.gs)

numnets = length(net.list)
numdata = nrow(data.disc[folds[[1]],])
numpga = length(midpoints)
numfolds = length(fit.causal)

prob.array = array(NA, c(numnets,numfolds, numdata, numpga))

point.est = array(NA,  c(numnets,numfolds, numdata, 3))
error = array(NA, c(numnets,numfolds, numdata, 3))
prob = array(NA, c(numnets, numfolds,numdata))

# out of sample test and cv error estimates
total = numnets *numdata *numfolds
pb = winProgressBar(title = "Progressbar",
                    label = "Progressbar", 
                    min = 0, 
                    max = total, 
                    width = 300)

#net  = net, i = data fold = fold
for(fold in 1:numfolds){
  for(net in 1: numnets){
    for(i in 1:numdata){
      
      testdata = data.disc[folds[[fold]],]
      
      P = cpdist(net.list[[net]][[fold]], nodes="PGA", 
                 evidence = (SD == testdata$SD[i] & 
                               Q0 == testdata$Q0[i]  &
                               kappa == testdata$kappa[i]  &
                               Vs30 == testdata$Vs30[i]  &
                               MAG == testdata$MAG[i]  &
                               DIST == testdata$DIST[i]),
                 n = 5000*nparams(net.list[[net]][[fold]]))
      
      prob.array[net,fold, i,] = prop.table(table(P))
      
      weighted.mean = sum(midpoints *  prob.array[net,fold, i,] )
      modal = max(midpoints[which(max(prob.array[net, fold, i,] ) == prob.array[net,fold, i,])]) # if several modals pick PGA highest
      mediann = midpoints[min(which(cumsum( prob.array[net, fold,i,] ) >= 0.5))]
      
      
      point.est[[net, fold,i,1]] = weighted.mean
      point.est[[net,fold,i,2]] = modal
      point.est[[net,fold,i,3]] = mediann
      prob[[net,fold,i]] =  prob.array[net, fold, i, testdata$PGA[i]]
      
      
     
     percantage = (net-1)*numfolds*numdata+(fold-1)*numdata+i
     
    setWinProgressBar(pb,
                      percantage, 
                      title = paste(round((percentage/total*100, digits = 3), "% done"), 
                                    label ="Progressbar")
    }
    
  }
  testdata_cont = data[folds[[fold]],]$PGA
  
  testarray = array(rep(testdata_cont, numnets*3), c(numnets,numdata,3))
  
  error[[,fold,,]] = (point.est[[,fold,,]] - testarray)^2
}

load("cv.RData")
# Errors
# validation set
validationError = error[,1,,] #only use first fold
meanErrors = apply(validationError, c(1,3), mean) # mean across dimension 2 (data)

prob[which(prob ==0)] = NA
probability = rowMeans(log(prob), na.rm = T)

results = cbind(meanErrors, probability)
colnames(results) = c("mean","mode","median","probability")
results

# cv

meanErrorsCv = apply(error, c(1,4), mean) # mean across dimensions 2,3 (folds, data)

# exclude probability = 0 since log(0) = not defined
prob[which(prob == 0)] = NA 
probability =  rowMeans(log(prob), na.rm = T) 

resultsCV = cbind(meanErrorsCv, probability)
colnames(resultsCV) = c("mean","mode","median","probability")
resultsCV

# EXPERIMENTAL ####

# BUild avarage models across the folds
# This is really nasty but couldn't figure out a way to do this more elegantly
# avaraging across the parameters and combining them to a new BN object...
net.avg = list(1,2,3,4)

for(net in 1:numnets){
  
  pga = net.list[[net]][[1]]$PGA
  sd = net.list[[net]][[1]]$SD
  q0 =net.list[[net]][[1]]$Q0
  Kappa = net.list[[net]][[1]]$kappa
  vs30 = net.list[[net]][[1]]$Vs30
  dist = net.list[[net]][[1]]$DIST
  m = net.list[[net]][[1]]$MAG
  
  if(net == 3) m = net.list[[net]][[2]]$MAG # first net of hc has different structure

  
  
  for(i in 2: 10){
    pga$prob = pga$prob + net.list[[net]][[i]]$PGA$prob
    sd$prob = sd$prob +net.list[[net]][[i]]$SD$prob
    q0$prob = q0$prob + net.list[[net]][[i]]$Q0$prob
    Kappa$prob = Kappa$prob + net.list[[net]][[i]]$kappa$prob
    vs30$prob = vs30$prob + net.list[[net]][[i]]$Vs30$prob
    dist$prob = dist$prob + net.list[[net]][[i]]$DIST$prob
    m$prob = m$prob + net.list[[net]][[i]]$MAG$prob
    print(i)
  }
  
  net.avg[[net]] = list("PGA" = pga, "SD" = sd, "Q0" = q0, "kappa" = Kappa, "Vs30" = vs30, "DIST" = dist, "MAG" = m)
  
  
  net.avg[[net]]$PGA$prob = net.avg[[net]]$PGA$prob/length(folds)
  net.avg[[net]]$SD$prob = net.avg[[net]]$SD$prob/length(folds)
  net.avg[[net]]$Q0$prob = net.avg[[net]]$Q0$prob/length(folds)
  net.avg[[net]]$kappa$prob = net.avg[[net]]$kappa$prob/length(folds)
  net.avg[[net]]$Vs30$prob = net.avg[[net]]$Vs30$prob/length(folds)
  net.avg[[net]]$DIST$prob = net.avg[[net]]$DIST$prob/length(folds)
  net.avg[[net]]$MAG$prob = net.avg[[net]]$MAG$prob/length(folds)
  
  class(net.avg[[net]]) = class(net.list[[1]][[1]])
  
}

# make predictions
# append the average model to each BN class = 4*(10 folds + 1 average model)
# average model is the "11th fold"
nets.list = list(c(net.list[[1]], net.avg[1]),
                 c(net.list[[2]], net.avg[2]),
                 c(net.list[[3]], net.avg[3]),
                 c(net.list[[4]], net.avg[4]))


pga_interval = c(min(data$PGA), -5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, max(data$PGA))

midpoints = 1:8
for(i in 1:8){
  midpoints[i] = mean(c(pga_interval[i],pga_interval[i+1]))
}


#initialize arrays to hold results
numnets = length(nets.list)
numdata = nrow(data.disc)
numpga = length(midpoints)
numfolds = length(nets.list[[1]])

prob.array_bv = array(NA, c(numnets,numfolds, numdata, numpga))

point.est_bv = array(NA,  c(numnets,numfolds, numdata, 3))
error_bv = array(NA, c(numnets,numfolds, numdata, 3))
prob_bv = array(NA, c(numnets, numfolds,numdata))

# set loading bar...
total = numnets *numdata *numfolds
pb = winProgressBar(title = "Progressbar",
                    label = "Progressbar", 
                    min = 0, 
                    max = total, 
                    width = 300)

# i = data net = net type fold = fold
for(fold in 1: length(nets.list[[net]])){
  for(net in 1: length(nets.list)){
    for(i in 1:nrow(data.disc)){
      
      P = cpdist(nets.list[[net]][[fold]], nodes="PGA", 
                 evidence = (SD == data.disc$SD[i] & 
                               Q0 == data.disc$Q0[i]  &
                               kappa == data.disc$kappa[i]  &
                               Vs30 == data.disc$Vs30[i]  &
                               MAG == data.disc$MAG[i]  &
                               DIST == data.disc$DIST[i]),
                 n=5000*nparams(nets.list[[net]][[fold]]))
      
      prob.array_bv[net, fold,i,] = prop.table(table(P))
      weighted.mean = sum(midpoints * prob.array_bv[net, fold,i,])
      modal = max(midpoints[which(max(prob.array_bv[net, fold,i,])==
                                    prob.array_bv[net, fold,i,])]) # if several modals pick PGA highest
      mediann = midpoints[min(which(cumsum(prob.array_bv[net, fold,i,]) >= 0.5))]
      
      
      point.est_bv[[net, fold,i,1]] = weighted.mean
      point.est_bv[[net, fold,i,2]] = modal
      point.est_bv[[net, fold,i,3]] = mediann
      
      error_bv[[net, fold,i,1]] = (weighted.mean - data$PGA[i])^2
      error_bv[[net, fold,i,2]] = (modal - data$PGA[i])^2
      error_bv[[net, fold,i,3]] = (mediann - data$PGA[i])^2
      
      
      prob_bv[[net, fold, i]] = prob.array_bv[net, fold,i,data.disc$PGA[i]]
      
      percentage = (net - 1) * numfolds * numdata + (fold - 1) * numdata + i
      setWinProgressBar(pb, percentage,
                        title = paste(round((percentage/total*100, digits = 3), "% done"), 
                                      label = "Progressbar")
    }
  }
}


# bv
load("bv.RData")
# bias = mean error of average model
#variance = mean error between other models and avereage

observedPGA = array(rep(data$PGA, 12, each = 4), c(numnets, numdata, 3))
biasError = (point.est_bv[,11,,] - observedPGA)^2

hist(biasError[1,,1])

bias = apply(biasError, c(1,3), mean)
colnames(bias) = c("mean","mode", "median")
bias


varianceError = array(NA, c(numnets, 10, numdata, 3))
for(i in 1:10){
  varianceError[,i,,]= (point.est_bv[,i,,] - point.est_bv[,11,,])^2
}

variance = apply(varianceError, c(1,4), mean)
colnames(variance) = c("mean","mode", "median")
variance

# Graphics for Export####

#allvsall

cairo_pdf(file = "..\\Figures\\allvsall.pdf",
          height = 9,
          width = 16)

plot(data, 
     main = "Scatterplot Matrix of all Variables", 
     col  = rgb(0, 0, 100, 50, maxColorValue = 255))

Sys.sleep(60)

dev.off()

# causal network
arcs(PGA_causal) = matrix(c("SD","MAG","MAG","PGA", "DIST","PGA", "Vs30", "PGA", "kappa","Vs30", "Q0", "Vs30" ),ncol = 2, byrow = TRUE, 
                          dimnames = list(c(),c("from","to")))

cairo_pdf(file = "..\\Figures\\causal.pdf",
          height = 9,
          width = 16)
graphviz.plot(PGA_causal,layout = "neato",main = "Causal Network")
dev.off()

# naive Bayes
PGA_naive = naive.bayes(data.disc, "PGA",names(data.disc)[-7])

cairo_pdf(file = "..\\Figures\\naive.pdf",
          height = 9,
          width = 16)
graphviz.plot(PGA_naive,layout = "fdp",main = "Naive Bayes Network")
dev.off()


# learned structure networks

x1 = list()
x2 = list()

for(i in 1:length(folds)){
  x1[i] = graphviz.plot(PGA_hc[[i]], layout = "circo")
  x2[i] = graphviz.plot(PGA_gs[[i]])  
}

# hill-climber single
cairo_pdf(file = "..\\Figures\\hc.pdf",
          height = 9,
          width = 16)

par(mfrow = c(5,2))
par(mar = rep(1, 4))
for(i in 1:length(folds)){
  
  x = layoutGraph(x1[[i]], nodeAttrs = list(width = 25))
  nodeRenderInfo(x) = list(fontsize = 12)
  graph.par(list(graph = list(main = paste("Fold: ",as.character(i)), 
                              cex.main = 2)))
  renderGraph(x)
}
dev.off()

#grow-shrink single
cairo_pdf(file = "..\\Figures\\gs.pdf",
          height = 9,
          width = 16)

par(mfrow = c(5,2))
par(mar = rep(1, 4))

for(i in 1:length(folds)){  
  x = layoutGraph(x2[[i]], nodeAttrs = list(width = 20))
  nodeRenderInfo(x) = list(fontsize = 12)
  graph.par(list(graph = list(main = paste("Fold: ",as.character(i)), 
                              cex.main = 2)))
  renderGraph(x)
}
dev.off()

# 
cairo_pdf(file = "..\\Figures\\hc_one.pdf",
          height = 9,
          width = 16)

x = layoutGraph(x1[[1]], nodeAttrs = list(width = 25))
nodeRenderInfo(x) = list(fontsize = 12)
graph.par(list(graph = list(main = "Learned Structure", sub = "Method: Hill-Climber"), 
               cex.main = 2,  cex.sub = 2))
renderGraph(x)

dev.off()


cairo_pdf(file = "..\\Figures\\gs_one.pdf",
          height = 9,
          width = 16)


x = layoutGraph(x2[[1]], nodeAttrs = list(width = 20))
nodeRenderInfo(x) = list(fontsize = 12)
graph.par(list(graph = list(main = "Learned Structure", sub = "Method: Grow-Shrink"), 
               cex.main = 2, cex.sub = 2))
renderGraph(x)

dev.off()

# Animation####

for(i in 1:1000){
  values = test_cont.liste[[1]][i,]
  header = sprintf("SD = %3.2f, Q0 = %3.2f, Kappa = %3.2f, Vs30 = %3.2f, MAG = %3.2f, Dist = %3.2f, PGA = %3.2f\n", 
                    values$SD, values$Q0, values$kappa, values$Vs30, values$MAG, values$DIST, values$PGA)
  
mp = barplot(m[1,i,], 
             main = header,
             axes = F, 
             ylim = c(0,max(m[1,i,])))

axis(1, at = mp, labels = midpoints)
axis(2, seq(0,1,0.05), seq(0,1,0.05))

Sys.sleep(0.5)
#readline(prompt="")
}


