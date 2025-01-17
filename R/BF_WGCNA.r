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

BF_WGCNA <- function(dataExp, bfmatrix,bfthr=6, keepedges=dim(dataExp)[2],softPower=(ifelse(keepedges != 0, 1,6)),plotTree=TRUE,plotfile="Gene_dendrogram_and_module_colors.pdf",trueModule=NULL,cutHeight = NULL, minClusterSize = 20, method = "hybrid",deepSplit = (ifelse(method == "hybrid", TRUE,FALSE)), maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, minSplitHeight = NULL, minAbsSplitHeight = NULL, externalBranchSplitFnc = NULL, minExternalSplit = NULL, externalSplitOptions = list(), externalSplitFncNeedsDistance = NULL, assumeSimpleExternalSpecification = TRUE,pamStage = TRUE, pamRespectsDendro = TRUE, useMedoids = FALSE, maxDistToLabel = NULL, maxPamDist = cutHeight, respectSmallClusters = TRUE,verbose = 2, indent = 0){


 if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, with no default. It is be a matrix or data frame containing expression data."); 
	}
if (is.null(dim(bfmatrix))) { 
		stop("argument \"bfmatrix\" is missing, with no default.");
	}


if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" is a matrix or data frame containing expression data.");}
if(! ("bfmatrix" %in% class(bfmatrix)) ){stop("argument \"bfmatrix\" is not a bfmatrix object, run Compute_bf first to obtain object bfmatrix.");}

if(plotTree==TRUE)
{
if(is.null(plotfile)){stop("Specify a file to output dendrogram plot.");}
}

#library("Matrix");
options(stringsAsFactors = FALSE);

if(bfthr<0){warning("argument \"bfthr\" is not set properly. It should be set largher than 0. ");}

geneid1<-as.character(bfmatrix$geneid1);
geneid2<-as.character(bfmatrix$geneid2);
genelist<-unique(c(geneid1,geneid2));
index1<-match(geneid1,genelist);
index2<-match(geneid2,genelist);
bfmatrix$bf.value[which(bfmatrix$bf.value<0)]<-0;
max<-log(max(bfmatrix$bf.value),10);
min<-max(0,log(min(bfmatrix$bf.value),10));
if(keepedges!=0){

 if(keepedges>dim(bfmatrix)[1]){warning("argument \"keepedges\" is not set appropriately. It's larger than the remained edges. All the edges will be used."); keepedges=dim(bfmatrix)[1];}
 if(keepedges<1){stop("argument \"keepedges\" is not set appropriately.")} 
 bfthr=sort(bfmatrix[,3],decreasing = TRUE)[keepedges];
}
if(max<min){stop("max BF value cannot less than min BF value.");}
trans2matrix(geneid1,geneid2,bfmatrix$bf.value,bfmatrix$type,genelist,min,max,bfthr)->x;
rm(bfmatrix);
gc();
x[[1]]<-as(x[[1]],"sparseMatrix");
gc();
x[[2]]<-as(x[[2]],"sparseMatrix");
gc();
x[[3]]<-as(x[[3]],"sparseMatrix");
gc();


if(is.null(dim(x[[2]]))){stop("size of adjacencies is null.");}

#adjacency = WGCNA::adjacency.fromSimilarity(x[[1]], power = softPower);
adjacency = x[[2]]^softPower;
dissTOM = 1-WGCNA::TOMsimilarity(as.matrix(adjacency));
rm(adjacency);
gc();

# Call the hierarchical clustering function

geneTree = flashClust::flashClust(as.dist(dissTOM), method = "average");
dynamicMods = dynamicTreeCut::cutreeDynamic(dendro=geneTree, cutHeight = cutHeight, minClusterSize = minClusterSize, method = method, 
    distM = dissTOM, deepSplit = deepSplit, maxCoreScatter = maxCoreScatter, minGap = minGap, maxAbsCoreScatter = maxAbsCoreScatter, 
    minAbsGap = minAbsGap, minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight, 
    externalBranchSplitFnc = externalBranchSplitFnc, minExternalSplit = minExternalSplit, externalSplitOptions = externalSplitOptions, 
    externalSplitFncNeedsDistance = externalSplitFncNeedsDistance, assumeSimpleExternalSpecification = assumeSimpleExternalSpecification, 
    pamStage = pamStage, pamRespectsDendro = pamRespectsDendro, useMedoids = useMedoids, 
    maxDistToLabel = maxDistToLabel, maxPamDist = maxPamDist, respectSmallClusters = respectSmallClusters, verbose = verbose, indent = indent);
rm(dissTOM);
gc();
dynamicColors = WGCNA::labels2colors(dynamicMods)
modulesNumber=length(table(dynamicMods));
colorOrder = c("grey", WGCNA::standardColors(modulesNumber));
moduleLabels = match(dynamicColors, colorOrder)-1;
modules=data.frame(gene=genelist,colors=dynamicColors,labels=moduleLabels);

# Plot the dendrogram and colors underneath
if(plotTree==TRUE)
{
pdf(plotfile);
if(is.null(colnames(dataExp)))
{
   colnames(dataExp)<-c(1:dim(dataExp)[2])
}
if(is.null(trueModule))
{
WGCNA::plotDendroAndColors(geneTree, dynamicColors, "BFDCA",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors");
}else{ if(length(trueModule)!=dim(dataExp)[2]){stop("argument \"trueModule\" and \"dataExp\" are not matching. They must have same number of genes.")}
selected_truemodule=trueModule[match(genelist,colnames(dataExp))];
WGCNA::plotDendroAndColors(geneTree, colors=data.frame(truemodule=selected_truemodule,BFDCA=dynamicColors),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main=paste("Gene dendrogram and module colors with softpower=",softPower,sep=""));
}
dev.off();
}

BFobt<-list(modules=modules,adjacencies=x[[2]],bf=x[[1]],type=x[[3]],geneTree=geneTree);
class(BFobt)<-append(class(BFobt),"BFobt");
rm(modules);rm(x);rm(geneTree);gc();
return (BFobt);
}

trans2matrix <-function(gene1,gene2,bf,type,genelist,min,max,bfthr){
               .Call('trans2matrix',geneid1=gene1,geneid2=gene2,bf=bf,type=type,genelist=genelist,min=min,max=max,bfthr=bfthr);
}

trans2adjacency <-function(gene1,gene2,bf,genelist,min,max,bfthr){
               .Call('trans2adjacency',geneid1=gene1,geneid2=gene2,bf=bf,genelist=genelist,min=min,max=max,bfthr=bfthr);
}
