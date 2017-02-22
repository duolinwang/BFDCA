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

sigDCpair_st1<- function(BFobt,mst2="mst2",bfthr=6,sparse=6,weight_cutoff=0.8)
{

if(!("BFobt" %in% class(BFobt))){stop("argument \"BFobt\" is not a BFobt object, please run BF_WGCNA to get BFobt first.");}
modules<-BFobt[[1]];
if (is.null(dim(modules))) { 
		stop("modules in BFobt is null.");
	}
as.matrix(read.table(mst2,header=FALSE,sep=" "))->mst2_data;

if(class(mst2_data[,1])!="character"){
mst2_data[,1]<-as.character(mst2_data[,1])
mst2_data[,2]<-as.character(mst2_data[,2])
}
modulesLabels<-modules$labels;
modulesColors<-modules$colors;
genelist<-modules$gene;
adjacency<-BFobt[[2]];
if (is.null(dim(adjacency))) { 
		stop("adjacency in BFobt is null.");
	}
bf<-BFobt[[3]];
rm(BFobt);
gc();
if (is.null(dim(bf))) { 
		stop("bf in BFobt is null.");
	}
ret<-NULL;
gsdeflist<-vector("list", max(modulesLabels)+1); 
for(i in c(1:dim(modules)[1]))
{
 gsdeflist[[modulesLabels[i]+1]]=append(gsdeflist[[modulesLabels[i]+1]],i);
}
if(length(gsdeflist)==1){stop("All genes are assigned to one group and they have no relationship with the class, change other parameters in function BF_WGCNA and try it again!");}
all<-NULL;
pr<-NULL;
for(i in c(2:length(gsdeflist)))
{
    object<-adjacency[gsdeflist[[i]],gsdeflist[[i]]]^sparse;
    e <- eigen(as.matrix(object))
    p <- matrix(abs(e$vectors[,1]))
    p <- p * norm(p)
    pr<-data.frame(modules$gene[gsdeflist[[i]]],p);
    all<-rbind(all,pr);
}
modules_weight<-cbind(modules,gene_weight=all[match(genelist,all[,1]),2])
bftest<-as.matrix(bf)
colnames(bftest)<-genelist;
rownames(bftest)<-genelist;
output_sigDC (mst2_data,modules_weight$labels,as.numeric(modules_weight$gene_weight),bftest,genelist,bfthr,weight_cutoff)->sig;
score=sig[,6]*sqrt(sig[,4]*sig[,5]);
ret<-data.frame(geneid1=sig[,1],geneid2=sig[,2],score=score);
ret<-ret[order(ret$score,decreasing = TRUE),];
return (ret);
}
output_sigDC = function(mst2_data,labels,weights,bf,genelist,bfthr,weight_cutoff){
.Call( 'output_sigDC',mst2_data=mst2_data,labels,weights=weights,bf=bf,genelist=genelist,bvthres=bfthr,cutoff=weight_cutoff);
}
