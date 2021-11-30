###CPS-merge analysis using simulation case1 as an example

###Load required packages and functions
library(Rcpp)
library(nnet)
library(OTclust)
library(leiden)
library(igraph)
library(mclust)
library(transport)
sourceCpp("huge.cpp")
sourceCpp("jaccard.cpp")
addnoise <- function(x,nrow,sd){
  x+rnorm(nrow,0,sd)
}

###Generate simulation data
set.seed(1230)
N=1000
#Generate component id 
component=sample(1:4,prob=c(0.25,0.25,0.25,0.25),size=N,replace=TRUE)
#View 1
mu1=matrix(0,ncol=10,nrow=4)
mu1[1,]=rep(5,10)
mu1[2,]=rep(5,10)
mu1[3,]=rep(10,10)
mu1[4,]=rep(10,10)
v1=apply(as.matrix(component,ncol=1),1,function(x){mu1[x,]+rnorm(10,0,1)})
v1=t(v1)
#View 2
mu2=matrix(0,ncol=10,nrow=4)
mu2[1,]=rep(10,10)
mu2[2,]=rep(5,10)
mu2[3,]=rep(10,10)
mu2[4,]=rep(5,10)
v2=apply(as.matrix(component,ncol=1),1,function(x){mu2[x,]+rnorm(10,0,1)})
v2=t(v2)

###Partition alignmnet under each view
#View 1
ref_1=kmeans(v1,4,iter.max=150,algorithm="MacQueen")$cluster
#Generate a collection of clustering results
save_1=matrix(0,ncol=101,nrow=length(ref_1))
save_1[,1]=ref_1
for(i in 2:101){
  inp=v1
  toy=apply(inp,2,addnoise,nrow=nrow(inp),sd=sqrt(0.1*mean(apply(inp,2,var))))
  cl=kmeans(toy,4,iter.max=150,algorithm="MacQueen")$cluster
  save_1[,i]=cl
}
#Align labels
cps1=align(save_1)
nbs=101-1
n=length(save_1)/(nbs+1)
K=cps1$numcls[-1]
new_1=matrix(0,nrow(save_1),ncol(save_1))
new_1[,1]=ref_1
for(i in 1:nbs){
  wtsub = cps1$weight[[i]]
  P_raw = matrix(0,n,K[i])
  for(j in 1:n){
    P_raw[j, save_1[j,i+1]] = 1
  }
  post=P_raw%*%(wtsub/ifelse(rowSums(wtsub)==0,1,rowSums(wtsub)))
  new_1[,i+1]=apply(post,1,which.is.max)
}
#View 2
ref_2=kmeans(v2,4,iter.max=150,algorithm="MacQueen")$cluster
#Generate a collection of clustering results
save_2=matrix(0,ncol=101,nrow=length(ref_2))
save_2[,1]=ref_2
for(i in 2:101){
  inp=v2
  toy=apply(inp,2,addnoise,nrow=nrow(inp),sd=sqrt(0.1*mean(apply(inp,2,var))))
  cl=kmeans(toy,4,iter.max=150,algorithm="MacQueen")$cluster
  save_2[,i]=cl
}
#Align labels
cps2=align(save_2)
nbs=101-1
n=length(save_1)/(nbs+1)
K=cps2$numcls[-1]
new_2=matrix(0,nrow(save_2),ncol(save_2))
new_2[,1]=ref_2
for(i in 1:nbs){
  wtsub = cps2$weight[[i]]
  P_raw = matrix(0,n,K[i])
  for(j in 1:n){
    P_raw[j, save_2[j,i+1]] = 1
  }
  post=P_raw%*%(wtsub/ifelse(rowSums(wtsub)==0,1,rowSums(wtsub)))
  new_2[,i+1]=apply(post,1,which.is.max)
}

###Bipartite clustering merging (it could be skipped in this example as the number of clusters is small)
#To skip this merging step, just run the following 2 lines before running "Tightness and Separability based Merging"
lab_combine=matrix(paste(as.vector(new_2),as.vector(new_1)),ncol=ncol(new_1))
lab_new=matrix(as.numeric(as.factor(as.vector(lab_combine))),ncol=ncol(lab_combine))
# nr=max(new_2)
# nc=max(new_1)
# w=matrix(0,ncol=nc,nrow=nr)
# for(i in 1:100){
#   wi=matrix(0,ncol=nc,nrow=nr)
#   can = sample(1:101, 2, replace=T)
#   rlab=new_1[,can[1]]
#   alab=new_2[,can[2]]
#   rs=as.numeric(names(table(rlab)))
#   as=as.numeric(names(table(alab)))
#   rlab=as.numeric(as.factor(rlab))
#   alab=as.numeric(as.factor(alab))
#   ref_i=cbind(rlab,alab)
#   icps=align(ref_i)
#   wj=icps$weight[[1]]
#   wj=(wj/ifelse(rowSums(wj)==0,1,rowSums(wj)))
#   wi[1:length(as),1:length(rs)]=wj
#   for(j in setdiff(seq(1,nr-1), as)){
#     l=sum(as>j)
#     wi[(j+1):(j+l),]=wi[j:(j+l-1),]
#     wi[j,]=0
#   }
#   for(j in setdiff(seq(1,nc-1), rs)){
#     l=sum(rs>j)
#     wi[,(j+1):(j+l)]=wi[,j:(j+l-1)]
#     wi[,j]=0
#   }
#   w=w+wi
# }
# w=w/100
# adjacency_matrix = w
# nc=ncol(adjacency_matrix)
# nr=nrow(adjacency_matrix)
# m0=matrix(0,nc,nc)
# m1=matrix(0,nr,nr)
# adjacency_matrix=rbind(cbind(m1,adjacency_matrix),cbind(t(adjacency_matrix),m0))
# bi_graph = graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
# bipartite_matrix <- igraph::as_adjacency_matrix(bi_graph, attr = "weight")
# partition <- leiden(bipartite_matrix,"ModularityVertexPartition")
# #Refine product clusters
# lab_combine=paste(as.vector(new_2),as.vector(new_1))
# comb=list()
# k=1
# for(i in 1:(max(partition)-1)){
#   for(j in (i+1):max(partition)){
#     c1=outer(which(partition[1:nr]==i),which(partition[(nr+1):(nr+nc)]==j),FUN = "paste")
#     c2=outer(which(partition[1:nr]==j),which(partition[(nr+1):(nr+nc)]==i),FUN = "paste")
#     comb[[k]]=c(c1)
#     comb[[k+1]]=c(c2)
#     k=k+2
#   }
# }
# for(i in 1:(k-1)){
#   lab_combine[lab_combine %in% comb[[i]]]=as.character(i)
# }
# lab_combine=matrix(lab_combine,ncol=ncol(new_1))
# lab_new=matrix(as.numeric(as.factor(as.vector(lab_combine))),ncol=ncol(lab_combine))

### Tightness and Separability based Merging
#Align again with combined reference partition
renumeric <- function(x){
  as.numeric(as.factor(x))
}
lab_new=apply(lab_new,2,renumeric)
nbs=101-1
n=length(lab_new)/(nbs+1)
lab_new_aligned=matrix(0,nrow(lab_new),ncol(lab_new))
lab_new_aligned[,1]=lab_new[,1]
k1 = max(lab_new[,1])
b=table(lab_new[,1])/n
for(i in 1:nbs){
  ki = max(lab_new[,i+1])
  d=Jac_dist(lab_new[,i+1]-1,lab_new[,1]-1,ki,k1)
  ot=transport(a=table(lab_new[,i+1])/n, b=b, method="networkflow", costm=d, fullreturn=TRUE)
  wtsub = ot$primal
  P_raw = matrix(0,n,ki)
  for(j in 1:n){
    P_raw[j, lab_new[j,i+1]] = 1
  }
  post=P_raw%*%(wtsub/ifelse(rowSums(wtsub)==0,1,rowSums(wtsub)))
  lab_new_aligned[,i+1]=apply(post,1,which.is.max)
}
#Merge rare clusters
lab_1=lab_new_aligned
lab_new=as.vector(lab_new_aligned)-1
fre=table(lab_new+1)
fre=fre[fre<100]
fre=as.numeric(names(fre))
for(i in 1:nrow(lab_1)){
  tt <- table(lab_1[i,])
  vote=names(tt[which.max(tt)])
  lab_1[i,lab_1[i,]%in% fre]=as.numeric(vote)
}

#CPS merging
lab_merge_t = as.numeric(as.factor(as.vector(lab_1+1)))-1
n_step_t = max(lab_merge_t)+1 - 4
ari_all_t = rep(0,n_step_t)
lab_all_t = matrix(0,nrow=nrow(lab_combine),ncol=n_step_t)
tit_all_t = matrix(0,nrow=2,ncol=n_step_t)
sep_all_t = rep(0,n_step_t)
for(i in 1:n_step_t){
  mcps_t=ACPS(lab_merge_t,101,1,1)
  min_tit_id=which.min(mcps_t$statistics[,4])
  dis_matrix_t=mcps_t$cap+diag(rep(2,ncol(mcps_t$cap)))
  merge_to_id=which.min(dis_matrix_t[min_tit_id,])
  min_tit_id = min_tit_id-1
  merge_to_id = merge_to_id-1
  lab_merge_t[lab_merge_t==min_tit_id]=merge_to_id
  lab_merge_t=as.numeric(as.factor(lab_merge_t))-1
  ari_all_t[i]=adjustedRandIndex(lab_merge_t[1:nrow(lab_combine)],component)
  lab_all_t[,i]=lab_merge_t[1:nrow(lab_combine)]
  sep_all_t[i] = dis_matrix_t[min_tit_id+1,merge_to_id+1]
  tit_all_t[1,i] = mcps_t$statistics[min_tit_id+1,4]
  tit_all_t[2,i] = mcps_t$statistics[merge_to_id+1,4]
  pen_t=mcps_t$match[,1]/apply(mcps_t$match,1,sum)
  tit_t=mcps_t$statistics[,4]*pen_t
}
#The ARI of final partition is 1
adjustedRandIndex(lab_all_t[,n_step_t],component)

### Cluster-wise contribution of each view
#Tightness-based contribution
v1_con=save_1-1
v1_con[,1]=lab_all_t[,n_step_t]
v1cps=ACPS(as.vector(v1_con),101,1,0)
pen_v1=v1cps$match[,1]/apply(v1cps$match,1,sum)
tit_v1=v1cps$statistics[,4]*pen_v1
v2_con=save_2-1
v2_con[,1]=lab_all_t[,n_step_t]
v2cps=ACPS(as.vector(v2_con),101,1,0)
pen_v2=v2cps$match[,1]/apply(v2cps$match,1,sum)
tit_v2=v2cps$statistics[,4]*pen_v2
con_matrix_tight = matrix(0,nrow=2,ncol=max(lab_all_t[,n_step_t]+1))
colnames(con_matrix_tight) = seq(1,max(lab_all_t[,n_step_t]+1),1)
rownames(con_matrix_tight) = c("View 1", "View 2")
con_matrix_tight[1,] = tit_v1/(tit_v1+tit_v2)
con_matrix_tight[2,] = tit_v2/(tit_v1+tit_v2)
con_matrix_tight[is.nan(con_matrix_tight)] = 0.5
con_matrix_tight

#Matching-weight-based contribution
v1align=align(v1_con+1)
w1=0
for(i in 1:100){
  wt = v1align$weight[[i]]
  wt=wt/ifelse(rowSums(wt)==0,1,rowSums(wt))
  w1=w1+apply(wt,2,mean)
}
w1=w1/100
v2align=align(v2_con+1)
w2=0
for(i in 1:100){
  wt = v2align$weight[[i]]
  wt=wt/ifelse(rowSums(wt)==0,1,rowSums(wt))
  w2=w2+apply(wt,2,mean)
}
w2=w2/100
con_matrix_weight = matrix(0,nrow=2,ncol=max(lab_all_t[,n_step_t]+1))
colnames(con_matrix_weight) = seq(1,max(lab_all_t[,n_step_t]+1),1)
rownames(con_matrix_weight) = c("View 1", "View 2")
con_matrix_weight[1,] = w1/(w1+w2)
con_matrix_weight[2,] = w2/(w1+w2)
con_matrix_weight[is.nan(con_matrix_weight)] = 0.5
con_matrix_weight