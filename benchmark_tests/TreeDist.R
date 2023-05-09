library(treespace)
library(adegenet)
library(adegraphics)
library(rgl)
library(ape)
library(castor)
library(dplyr)
library(ggplot2)
library(formattable)
library(lubridate)
library(ggnewscale)
library(ggpubr)
library(TreeDist)

args = commandArgs(trailingOnly=TRUE)
print(args[1])
st1<-read.tree(args[1])
distances <- RobinsonFoulds(st1)
mapping <- cmdscale(distances, k = 3)

plot(mapping,
     asp = 1,
     ann = FALSE, axes = FALSE,
     col = c(1:args[2]), pch = 16
     )
print('Done!')
