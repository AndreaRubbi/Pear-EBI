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

args = commandArgs(trailingOnly=TRUE)
print(args[1])
st1<-read.tree(args[1])
res_seeds<-treespace(st1,"RF",nf=2)
df<-data.frame(A1=res_seeds$pco$li$A1,A2=res_seeds$pco$li$A2)
df["type"]<-replicate(args[2],"st1")
df["number"]<-seq(1,args[2],1)

a=1
fig<-ggplot(data=filter(df, type== "st1"),aes(x=A1,y=A2,col=number))+geom_point(alpha=1.0,size=a)+
  scale_color_gradient(low="yellow",high="red",name="Tree number")+
  ggtitle("st1")+
  theme(legend.title = element_text(color = "black", size = 17),
  legend.text = element_text(color = "black", size = 17),
  plot.title = element_text(color = "black", size = 17, face = "bold"),
  plot.subtitle = element_text(color = "black",size=17),
)
fig

print('Done!')
