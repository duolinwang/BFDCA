##############################################################################
#Copyright (C) 2016 Duolin Wang

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
###############################################################################

BF_output_networks<- function(dataExp,class,classlabel,BFobt,mst2file="MST2.txt",bfthr=6,corthres=0.3,echo=FALSE)
{
    if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, with no default. It is a matrix or data frame containing expression data."); 
	}

   if(is.null(length(class))){stop("argument \"class\" is missing, with no default.")}

   if(is.null(length(classlabel))){stop("argument \"classlabel\" is missing, with no default.");}

   if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" must be a matrix or data frame containing expression data.");}   
   if(length(class)!=dim(dataExp)[1]){stop(" argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
   if(length(table(class))>2){stop(" by far package BFDCA can only accept binary conditions, but class contains more than 2 classes. ")}   
   if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}   
   if(!("BFobt" %in% class(BFobt))){stop("argument \"BFobt\" is not a BFobt object, please run BF_WGCNA to get BFobt first.");}

#library(Matrix);

if(is.null(mst2file)){stop("Specify a correct file name in argument 'mst2file' to output mst2 skeletons.");}

if (is.null(dim(BFobt[[1]]))) { 
		stop("modules in BFobt is missing, please check BFobt object.");
	}

modulesLabels<-BFobt[[1]]$labels;
modulesColors<-BFobt[[1]]$colors;
genelist<-BFobt[[1]]$gene;

if (is.null(dim(BFobt[[2]]))) { 
		stop("adjacency in BFobt is null.");
	}


if (is.null(dim(BFobt[[3]]))) { 
		stop("bf in BFobt is null.");
	}
gc();

ret<-NULL;
gsdeflist<-vector("list", max(modulesLabels)+1); 
for(i in c(1:dim(BFobt[[1]])[1]))
{
 gsdeflist[[modulesLabels[i]+1]]=append(gsdeflist[[modulesLabels[i]+1]],i);
}

if(length(gsdeflist)==1){stop("All genes are assigned to one group and they have no relationship with the class, change other parameters in function BF_WGCNA and try it again!");}

#library(igraph);
all<-igraph::graph.empty(n=0, directed=FALSE);

adj<-as.matrix(BFobt[[2]]);
gc();
colnames(adj)<-genelist;
rownames(adj)<-genelist;
for(i in c(2:length(gsdeflist)))
{ 
    distmat <- 1-adj[gsdeflist[[i]],gsdeflist[[i]]];
    gr <- igraph::graph.adjacency(as.matrix(distmat), weighted=TRUE, mode="undirected")
    first.mst <- igraph::minimum.spanning.tree(gr)
    mst1.matrix <- igraph::get.adjacency(first.mst, attr="weight", sparse=FALSE)
    distmat2 <- distmat - mst1.matrix
    gr2 <- igraph::graph.adjacency(as.matrix(distmat2), weighted=TRUE, mode="undirected")
    second.mst <- igraph::minimum.spanning.tree(gr2)
    all<-all+first.mst+second.mst;     
}
tryCatch(igraph::write_graph(all,mst2file, "ncol"),error=function(cond){stop(paste("Can't output mst2 skeletons into file ",mst2file,".",sep=""));})

rm(distmat);
rm(distmat2);
gc();
all<-NULL;
pr<-NULL;


for(i in c(2:length(gsdeflist)))
{
    object<-adj[gsdeflist[[i]],gsdeflist[[i]]];
    e <- eigen(as.matrix(object))
    p <- matrix(abs(e$vectors[,1]))
    p <- p * norm(p)
    pr<-data.frame(genelist[gsdeflist[[i]]],p);
    all<-rbind(all,pr);
}
rm(object);
gc();
rm(adj)
gc();
mean1<-NULL;
mean2<-NULL;

modules_weight<-cbind(BFobt[[1]],gene_weight=all[match(genelist,all[,1]),2])
ret$genegroups<-modules_weight[order(modules_weight[,3]),];
rownames(ret$genegroups)<-ret$genegroups$gene;
ret$genegroups[ret$genegroups$labels!=0,]->ret$genegroups;

dataExp<-as.matrix(dataExp);
 if(is.null(colnames(dataExp)))
   {
    colnames(dataExp)<-c(1:dim(dataExp)[2])
   }

datac1<-dataExp[class==classlabel[1],genelist];
datac2<-dataExp[class==classlabel[2],genelist];

medianc1<-NULL;
medianc2<-NULL;
for(i in c(1:dim(datac1)[2]))
{
medianc1[i]<-mean(datac1[,i]);
medianc2[i]<-mean(datac2[,i]);
}

medianc1[is.na(medianc1)]<-0;
medianc2[is.na(medianc2)]<-0;
if(echo==TRUE){echo=1}else{echo=0};
output_network (modulesLabels,as.numeric(modules_weight$gene_weight),as.matrix(BFobt[[3]]),as.matrix(BFobt[[4]]),genelist,datac1,datac2,medianc1,medianc2,gsdeflist[2:length(gsdeflist)],corthres,bfthr,echo)->ret$network;
rm(BFobt);
gc();
return (ret);
}
output_network = function(labels,weights,bf,type,genelist,datac1,datac2,medianc1,medianc2,gsdeflist,corthres,bfthr,echo){
.Call( 'output_network',labels=labels,weights=weights,bf=bf,type=type,genelist=genelist,datac1=datac1,datac2=datac2,medianc1=medianc1,medianc2=medianc2,gsdeflist=gsdeflist,corthres=corthres,bfthr=bfthr,echo=echo);
}

