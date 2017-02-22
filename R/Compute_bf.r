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
Compute_bf <- function(dataExp, class, classlabel, bfthr=6,echo=FALSE){

   if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, with no default. It must be a matrix or data frame containing expression data."); 
	}

   if(is.null(length(class))){stop("argument \"class\" is missing, with no default.")}

   if(is.null(length(classlabel))){stop("argument \"classlabel\" is missing, with no default.");}

   if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" must be a matrix or data frame containing expression data.");}
   if(length(class)!=dim(dataExp)[1]){stop(" argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
   if(length(table(class))>2){stop(" by far package BFDCA can only accept binary conditions, but class contains more than 2 classes. ")}   
   if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}   
   if(is.null(bfthr)){bfthr=-Inf};
    
   class_l<-array(NA,2);
   class_l[class==classlabel[1]]=1;
   class_l[class==classlabel[2]]=2;
   
   dataExp<-as.matrix(dataExp);
   dataExp<-scale(dataExp);
   for(i in c(1:dim(dataExp)[2]))
   {
    if(var(dataExp[class_l==1,i])==0)stop(paste("the variance of gene with index ",i," with classlabel=\" ",classlabel[1],"\" is 0. Delete this gene and run again.",sep=""))
   }

   for(i in c(1:dim(dataExp)[2]))
   {
    if(var(dataExp[class_l==2,i])==0)stop(paste("the variance of gene with index ",i," with classlabel=\"",classlabel[2],"\" is 0. Delete this gene and run again.",sep=""))
   }


   if(!is.null(colnames(dataExp)))
   {
   geneid<-colnames(dataExp);
   }else{geneid<-c(1:dim(dataExp)[2])}

   ret=NULL;   
   if(echo==TRUE){echo=1;}else{echo=0};
   ret=BF1_pair(dataExp,class_l,bfthr,echo);   

rett<-data.frame(geneid1=geneid[ret[[1]]],geneid2=geneid[ret[[2]]],bf.value=ret[[3]],type=ret[[4]])
  
class(rett)<-append(class(rett),"bfmatrix");
return (rett);
}

colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
                     sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}


BF1_pair <- function(fixed.gs.data,class,bfthr,echo) {
	.Call('BF1_pair',data=as.matrix(fixed.gs.data),group=class,thres=bfthr,echo=echo)
}                        

