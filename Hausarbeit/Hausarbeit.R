# Loading data and initilization ####

# Set working directory

path = getwd()
setwd(path)

# Clear the memory 
rm(list = ls())

#style
op = par()
settings = options()
par(pch = 19, ann = T, cex = 1)


# load packages/install packages
#TODO
#install.packages("bnlearn")


library(bnlearn)
library(Rgraphviz)

# read in additional functions
source(".\\functions\\crossvalidation function.R")

# read data
data <- read.csv(file = "data.dat") 

# visualize data ####

plot(data, 
     main = "Scatterplot Matrix of all Variables", 
     col  = rgb(0, 0, 100, 50, maxColorValue = 255))

#options(ylim = c(0, 1))

par(mfrow=c(3,2))
hist(data$SD, freq = F)
hist(data$Q0, freq = F)
hist(data$kappa, freq = F)
hist(data$Vs30, freq = F)
hist(data$MAG, freq = F)
hist(data$DIST, freq = F)



par(mfrow=c(3,2))
hist(data$SD, breaks = c(0, 0.8792, 5.438, 14.92, 58, 500), freq = F)
hist(data$Q0, breaks = c(0, 330, 5000), freq = F)
hist(data$kappa, breaks = c(0, 0.01053, 0.0345, 0.1), freq = F)
hist(data$Vs30, breaks = c(600, 1704.5, 2800), freq = F)
hist(data$MAG, breaks = c(5, 6.271, 7.5), freq = F)
hist(data$DIST, breaks = c(1, 4.38, 15.4885, 55.84, 200), freq = F)
box()

par(op)
options(settings)
hist(data$PGA, c(min(data$PGA),-5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, max(data$PGA)))

# divide iinto learn and test data ####

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

liste = createFolds(data.disc, 10)

folds = liste$folds
learning.liste = liste$learning
test_disc.liste = liste$test
test_cont.liste = vector("list", length(folds))
for(i in 1:length(folds)){
  test_cont.liste[[i]] = data[folds[[i]],]
  
}
  
rm(liste)

#set up the nets ####

PGA_causal = empty.graph(names(data))
PGA_naive = empty.graph(names(data))

# causal network
arcs(PGA_causal) = matrix(c("SD","MAG","MAG","PGA", "DIST","PGA", "Vs30", "PGA", "kappa","Vs30", "Q0", "Vs30" ),ncol = 2, byrow = TRUE, 
                       dimnames = list(c(),c("from","to")))
graphviz.plot(PGA_causal,layout = "neato",main = "Causal Network")

# naive Bayes
PGA_naive = naive.bayes(data.disc, "PGA",names(data.disc)[-7])

graphviz.plot(PGA_naive,layout = "fdp",main = "Naive Bayes Network")


#Learning and testing ####

# allocate lists
fit.causal = vector("list", length(folds))
fit.naive = vector("list", length(folds))

for(i in 1:length(folds)){
  fit.causal[[i]] = bn.fit(PGA_causal, learning.liste[[i]], method = "bayes")
  fit.naive[[i]] = bn.fit(PGA_naive, learning.liste[[i]], method = "bayes")
}

# Sturktur gelernt

#allocate lists
PGA_hc = vector("list", length(folds))
PGA_gs = vector("list", length(folds))
fit.hc = vector("list", length(folds))
fit.gs = vector("list", length(folds))

for(i in 1:length(folds)){
  # Score based hill-climber
  PGA_hc[[i]] = hc(learning.liste[[i]]) 
  
  # constraint based grow shrink 
  PGA_gs[[i]] = gs(learning.liste[[i]]) 
  
  # Parameter lernen 
  fit.hc[[i]] = bn.fit(PGA_hc[[i]], learning.liste[[i]])
  fit.gs[[i]] = bn.fit(PGA_gs[[i]], learning.liste[[i]])
}

par(mfrow = c(5,2))
par(mar = rep(2, 4))
for(i in 1:length(folds)){
  
graphviz.plot(PGA_hc[[i]],layout = "dot", main = paste('Learned Structure, Fold: ',as.character(i)), sub = "Method: Hill-Climber")
#readline(prompt = "")
}

for(i in 1:length(folds)){  
graphviz.plot(PGA_gs[[i]],layout = "fdp",main = paste('Learned Structure, Fold: ',as.character(i)), sub = "Method: Growth-Shrink")
#readline(prompt = "")
}

par(op)

# Validation#### 

pga_interval = c(min(data$PGA), -5.135, -3.722, -2.627, -1.20742, 0.145, 1.657, 3.175, max(data$PGA))

midpoints = 1:8
for(i in 1:8){
  midpoints[i] = mean(c(pga_interval[i],pga_interval[i+1]))
}

# just one fold and model each
net.list = list(fit.causal, fit.naive, fit.hc, fit.gs)

numnets = length(net.list)
numdata = nrow(test_disc.liste[[1]])
numpga = length(midpoints)
numfolds = length(fit.causal)

prob.array = array(NA, c(numnets, numfolds, numdata, numpga))

point.est = array(NA,  c(numnets, numfolds, numdata, 3))
error = array(NA, c(numnets, numfolds, numdata, 3))
prob = array(NA, c(numnets, numfolds, numdata))

# out of sample test and cv error estimates

#net  = net, i = data fold = fold
for(net in 1: numnets){
  for(fold in 1: numfolds){
    for(i in 1:numdata){
      
      P = cpdist(net.list[[net]][[fold]], nodes="PGA", 
                 evidence = (SD == test_disc.liste[[fold]]$SD[i] & 
                               Q0 == test_disc.liste[[fold]]$Q0[i]  &
                               kappa == test_disc.liste[[fold]]$kappa[i]  &
                               Vs30 == test_disc.liste[[fold]]$Vs30[i]  &
                               MAG == test_disc.liste[[fold]]$MAG[i]  &
                               DIST == test_disc.liste[[fold]]$DIST[i]),
                 n=5000*nparams(net.list[[net]][[fold]]))
      
      prob.array[net, fold, i,]   = prop.table(table(P))
      
      weighted.mean = sum(midpoints *  prob.array[net,fold, i,] )
      modal = max(midpoints[which(max(prob.array[net, fold, i,] )== prob.array[net, fold, i,] )]) # if several modals pick PGA highest
      mediann = midpoints[min(which(cumsum( prob.array[net,fold,i,] ) >= 0.5))]
      
      
      point.est[[net, fold, i,1]] = weighted.mean
      point.est[[net,fold,i,2]] = modal
      point.est[[net,fold,i,3]] = mediann
      
      error[[net,fold,i,1]] = (weighted.mean - test_cont.liste[[fold]]$PGA[i])^2
      error[[net,fold,i,2]] = (modal - test_cont.liste[[fold]]$PGA[i])^2
      error[[net,fold,i,3]] = (mediann - test_cont.liste[[fold]]$PGA[i])^2
      
      prob[[net, fold,i]] =  prob.array[net,fold,i,test_disc.liste[[fold]]$PGA[i]]
      
    }
  }
}

# BUild avarage models across the folds
net.avg = list(1:length(net.list))

for(net in numnets){
  
  pga = net.list[[net]][[1]]$PGA
  sd = net.list[[net]][[1]]$SD
  q0 =net.list[[net]][[1]]$Q0
  Kappa = net.list[[net]][[1]]$kappa
  vs30 = net.list[[net]][[1]]$Vs30
  dist = net.list[[net]][[1]]$DIST
  mag = net.list[[net]][[1]]$MAG
  
  for(i in 2: 10){
    pga$prob = pga$prob + net.list[[net]][[i]]$PGA$prob
    sd$prob = sd$prob +net.list[[net]][[i]]$SD$prob
    q0$prob = q0$prob + net.list[[net]][[i]]$Q0$prob
    Kappa$prob = Kappa$prob + net.list[[net]][[i]]$kappa$prob
    vs30$prob = vs30$prob + net.list[[net]][[i]]$Vs30$prob
    dist$prob = dist$prob + net.list[[net]][[i]]$DIST$prob
    mag$prob = mag$prob + net.list[[net]][[i]]$MAG$prob
    
  }
  net.avg[[net]] = list ("PGA" = pga, "SD" = sd, "Q0" = q0, "kappa" = Kappa, "Vs30" = vs30, "DIST" = dist, "MAG" = mag)
  
  
  net.avg[[net]]$PGA$prob = net.avg[[net]]$PGA$prob/length(folds)
  net.avg[[net]]$SD$prob = net.avg[[net]]$SD$prob/length(folds)
  net.avg[[net]]$Q0$prob = net.avg[[net]]$Q0$prob/length(folds)
  net.avg[[net]]$kappa$prob = net.avg[[net]]$kappa$prob/length(folds)
  net.avg[[net]]$Vs30$prob = net.avg[[net]]$Vs30$prob/length(folds)
  net.avg[[net]]$DIST$prob = net.avg[[net]]$DIST$prob/length(folds)
  net.avg[[net]]$MAG$prob = net.avg[[net]]$MAG$prob/length(folds)
  
  class(net.avg[[net]]) = class(net.list[[net]][[1]])
  
}


# make predictions
nets.list = net.list
for(i in 1:length(nets.list)){
  nets.list[[i]] = c(nets.list[[i]],net.avg[i])
}




prob.array_bv = array(NA, c(numnets, numfolds, nrow(data), numpga))

point.est_bv = array(NA,  c(numnets, numfolds, nrow(data), 3))
error_bv = array(NA, c(numnets, numfolds, nrow(data), 3))
prob_bv = array(NA, c(numnets, numfolds, nrow(data)))

# i = data net = net type
for(net in 1: length(nets.list)){
  for(fold in 1: length(nets.list[[net]])){
    for(i in 1:nrow(data.disc)){
      
      P = cpdist(nets.list[[net]][[fold]], nodes="PGA", 
                 evidence = (SD == data.disc$SD[i] & 
                               Q0 == data.disc$Q0[i]  &
                               kappa == data.disc$kappa[i]  &
                               Vs30 == data.disc$Vs30[i]  &
                               MAG == data.disc$MAG[i]  &
                               DIST == data.disc$DIST[i]),
                 n=5000*nparams(nets.list[[net]][[fold]]))
      
      prob.array_bv[net, fold,i,]   = prop.table(table(P))
      weighted.mean = sum(midpoints *prob.array_bv[net, fold,i,])
      modal = max(midpoints[which(max(prob.array_bv[net, fold,i,])==prob.array_bv[net, fold,i,])]) # if several modals pick PGA highest
      mediann = midpoints[min(which(cumsum(prob.array_bv[net, fold,i,]) >= 0.5))]
      
      
      point.est_bv[[net, fold,i,1]] = weighted.mean
      point.est_bv[[net, fold,i,2]] = modal
      point.est_bv[[net, fold,i,3]] = mediann
      
      error_bv[[net, fold,i,1]] = (weighted.mean - data$PGA[i])^2
      error_bv[[net, fold,i,2]] = (modal - data$PGA[i])^2
      error_bv[[net, fold,i,3]] = (mediann - data$PGA[i])^2
      
      
      
      prob_bv[[net, fold, i]] = prob.array_bv[net, fold,i,data.disc$PGA[i]]
      
    }
    
  }
}

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


