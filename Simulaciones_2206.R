######################################################################
################### Simulations
######################################################################

# Authors: Santiago Ortiz (santiagoortiz00@usc.edu.co - sortiza2@eafit.edu.co - saortizar@unal.edu.co - santy_ortiz@hotmail.com)
# Jose Londoño (josclogo@gmail.com-jose.londono94076@u.icesi.edu.co)
# Date: 12/2023

set.seed(50)

## Required Packages
my_packages = c("expm","gridExtra","tidyverse","knitr","kableExtra",
                "IRdisplay","mrfDepth","e1071","pracma","MASS","mixtools",
                "doParallel","foreach","extraDistr","rrcov","Rfast",
                "dobin","FNN","parallel","depthTools","OutliersO3", "mvoutlier",
                "fds", "caret", "tictoc")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {
  library(my_packages[q], character.only = TRUE)
}

## Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Load outlier detection codes and auxiliary routines
source("OL_skew_self_V2_1805.R")

ini<-Sys.time()

iter = 1
alpha = c(0.1,0.2,0.3,0.4)
dimension = c(10,20,40)
cant_datos = 10*dimension
outl_dist = c(3,4)
outl_concent = c(0.5)
alpha = c(0.4)
dimension = c(10)
cant_datos = 10*dimension
outl_dist = c(4)
outl_concent = c(0.5)

#########################################################################################################
####### Multivariate normal distribution N(0,I) contaminated with (alpha/2)*N(0 + delta,lambda*I) #######
#########################################################################################################

registerDoParallel(detectCores())

#res4 = foreach (i2 = 1:length(dimension), .combine = rbind) %do% {
  #res3 = foreach (i1 = 1:length(alpha), .combine = rbind) %do% {
    #res2 = foreach (i3 = 1:length(outl_dist), .combine = rbind) %do% {
      #res1 = foreach (i4 = 1:length(outl_concent), .combine = rbind) %do% {
        #res0 = foreach (i5 = 1:iter, .combine = rbind) %dopar% {
resul = matrix(0, iter, 4)
resul.t = matrix(0, 0, 4)
for (i2 in 1:length(dimension)) {
  for (i1 in 1:length(alpha)) {
    for (i3 in 1:length(outl_dist)) {
      for (i4 in 1:length(outl_concent)) {
        for (i5 in 1:iter) {
          datos_limp = round(cant_datos[i2]*(1-alpha[i1]))
          datos_cont = cant_datos[i2]-datos_limp
          data.sim = GenAtip(cant_datos[i2],dimension[i2],c(alpha[i1],outl_dist[i3],outl_concent[i4]),1)
          label.X = as.factor(data.sim$lbl)
          X = data.sim$x
          
          ##
           pairs(X[,1:5], pch = 20)
          ##
          
          # RSP-SAO
          Sk = max_skew(X)
          Sk.proy = X %*% Sk$dv
          # hist(Sk.proy)
          Sk.dirs = Sk.Child(X, Sk.proy)
          Skproy.Child = X %*% Sk.dirs
          #Adj.SDO = apply(apply(cbind(Skproy.Child, Sk.proy), 2, adj.outly), 1, max)
          proy.t = cbind(Skproy.Child, Sk.proy)
          Adj.SDO = apply((proy.t - apply(proy.t, 2, median)) / apply(proy.t, 2, Qn), 1, max)
          # hist(Adj.SDO)
          out.points = Adj.Outlier(Adj.SDO, 1)
          W0.cm = confusionMatrix(as.factor(out.points), label.X)$table
          #W0.cm
          
          # SAO
          W1 = as.numeric(adjOutlyingness(X)$nonOut)
          W1[which(W1 == 1)] = 2
          W1[which(W1 == 0)] = 1
          W1[which(W1 == 2)] = 0
          W1.cm = confusionMatrix(as.factor(W1), label.X)$table
          
          #pairs(X[,1:5],
          pairs(X[,1:5],
                col=ifelse((out.points==1&W1==1),"green",ifelse(W1==1,"blue",ifelse(out.points==1,"red","black"))),
                pch=ifelse((out.points==1&W1==1),8,ifelse(W1==1,10,ifelse(out.points==1,12,14))),
                main="Gráfico de Pares con Datos Atípicos Resaltados")
    
          # RESULTS          
          resul[i5,] = c(W0.cm[2,2]/sum(W0.cm[,2]), W0.cm[2,1]/sum(W0.cm[,1]),
                    W1.cm[2,2]/sum(W1.cm[,2]), W1.cm[2,1]/sum(W1.cm[,1]))
        }
        resul1 = apply(resul, 2, mean)
      }
      resul.t = rbind(resul.t, resul1)
      cat(i2, i1, i3, i4)
    }
  }
}
#print("Experiment Symmetric TypeA has finished")
save(resul.t, file = "Symm.RData")
write.csv(resul.t, file = "Symm.csv")

stopImplicitCluster()

fin<-Sys.time();fin-ini

###########
### FIN ###
###########