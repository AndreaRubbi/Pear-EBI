library("ape", lib.loc="path/to/libloc")
setwd("random_trees")

args = commandArgs(trailingOnly=TRUE)

num_trees_list<-args[1]
num_tips_list<- args[2]

for (idx in 1:length(num_trees_list)){
for (jdx in 1:length(num_tips_list)){
num_trees = num_trees_list[idx]
num_tips = num_tips_list[jdx]
file_name = paste0("random_",toString(num_trees),"_trees_",toString(num_tips),"_tips")

start_time<-Sys.time()
tree<-rmtree(num_trees, num_tips)
write.tree(tree,file_name,append=FALSE)
end_time<-Sys.time()
print(c(num_tips,num_trees,as.double(end_time-start_time,units="secs")))
}}
