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

BFtrain<-function(dataExp,class,classlabel,edges,bfthr=6){
   if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, with no default. It must be a matrix or data frame containing expression data."); 
	}

   if(is.null(length(class))){stop("argument \"class\" is missing, with no default.")}

   if(is.null(length(classlabel))){stop("argument \"classlabel\" is missing, with no default.");}

   if(class(dataExp)=="data.frame"){as.matrix(dataExp)->dataExp;}
   if(is.null(dim(dataExp))){matrix(dataExp,nrow=1)->dataExp;} 
   if(length(class)!=dim(dataExp)[1]){stop(" argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
   if(length(table(class))>2){stop(" by far package BFDCA can only accept binary conditions, but class contains more than 2 classes. ")}   
   if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}  
   
   if(bfthr<0){warning("argument \"bfthr\" is not set properly. It should be set largher than 0. ");} 
   if(!("matrix" %in% class(edges)|"data.frame" %in% class(edges))){stop("argument \"edges\" must be a data frame or matrix containing at least two columns, representing ids (or indexes) for gene1 and ids (or indexes) for gene2. ");}

   if(dim(edges)[2]<2){stop("argument \"edges\" must be a data frame or matrix containing at least two columns, representing ids (or indexes) for gene1 and ids (or indexes) for gene2. ")}
   as.matrix(edges)->edges; # data.frame will be transferred into matrix. items will be characters.
   ret=NULL;
   rett=list();
   class_l<-array(NA,2);
   class_l[class==classlabel[1]]=1;
   class_l[class==classlabel[2]]=2;
   
   rett[[1]]=list(classlabel,edges);

if(is.null(colnames(dataExp)))
{
   colnames(dataExp)<-c(1:dim(dataExp)[2])
}   

   for(j in 1:dim(edges)[1])
   {
    tryCatch({traindata=dataExp[,c(edges[j,1],edges[j,2])]},error=function(cond){stop(paste("edge: ",edges[j,1]," ",edges[j,2]," is not in dataExp. Check colnames of the training data. ",sep=""))})
    ret=BF_train(traindata,class_l,bfthr);
    rett[[j+1]]=ret; 
   }
class(rett)<-append(class(rett),"DCmodel");
return (rett);
}

BF_train <- function(fixed.data,class,bfthr) {
	.Call('BF_train',dataExp=fixed.data,group=class,thres=bfthr)
}                        

