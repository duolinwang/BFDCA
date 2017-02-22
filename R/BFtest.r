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

BFtest <- function(testdata,model){
 #if (!("matrix" %in% class(testdata) | "data.frame" %in% class(testdata))) { 
#		stop("argument \"testdata\" must data frame or matrix containing expression of testing data. columns correspond to genes and rows to samples. the columns must be consistent with the training data. "); 
#	}

 if(is.null(dim(testdata))){matrix(testdata,nrow=1)->testdata;}
 if(!("DCmodel" %in% class(model))){stop("argument \"model\" is not a DCmodel object, please run BFtrain to get DCmodel first. ");}
 ret=NULL;
 class_l<-(model[[1]][1])[[1]];
 edges<-(model[[1]][2])[[1]];
 if(is.null(colnames(testdata)))
 {
   colnames(testdata)<-c(1:dim(testdata)[2])
 } 

for(s in c(1:dim(testdata)[1]))
{
p1<-0;
p2<-0;
for(i in 1:dim(edges)[1])
{
if(model[[i+1]][[1]]!=0)
{
e<-c(as.character(edges[i,1]),as.character(edges[i,2]));
d<-length(e);
has=1;
tryCatch({meanc1<-matrix(as.numeric(testdata[s,e]-model[[i+1]][[3]]),ncol=d)},error=function(cond){warning(paste("edge ",edges[i,1]," ",edges[i,2],"is not in data. This edge will be omitted in this procedure.",sep=""));has=0;});
tryCatch({meanc2<-matrix(as.numeric(testdata[s,e]-model[[i+1]][[4]]),ncol=d)},error=function(cond){warning(paste("edge ",edges[i,1]," ",edges[i,2],"is not in data. This edge will be omitted in this procedure.",sep=""));has=0;});

if(has==1)
{
if(model[[i+1]][[1]]==1)
{
sv1<-tryCatch({solve(model[[i+1]][[5]])},error=function(cond){return(NA)});
sv2<-tryCatch({solve(model[[i+1]][[6]])},error=function(cond){return(NA)});
if(!is.na(sv1[1]) && !is.na(sv2[1]))
{
p1=p1-((meanc1%*%sv1)%*%t(meanc1));
p2=p2-((meanc2%*%sv2)%*%t(meanc2));
}
}
if(model[[i+1]][[1]]==2)
{
sv1<-tryCatch({solve(model[[i+1]][[5]])},error=function(cond){return(NA)});
if(!is.na(sv1[1]))
{
 p1=p1-log(det(model[[i+1]][[5]]))-((meanc1%*%sv1)%*%t(meanc1)); 
 p2=p2-sum(log(model[[i+1]][[6]]))-sum(as.numeric(meanc2)^2/model[[i+1]][[6]]);
}
}
if(model[[i+1]][[1]]==3)
{
sv2<-tryCatch({solve(model[[i+1]][[6]])},error=function(cond) {return(NA)});
if(!is.na(sv2[1]))
{
p1=p1-sum(log(model[[i+1]][[5]]))-sum(as.numeric(meanc1)^2/model[[i+1]][[5]]);
p2=p2-log(det(model[[i+1]][[6]]))-(meanc2%*%sv2)%*%t(meanc2);
}
}
}
}
}
if(p1>p2){ret[s]=class_l[1];}else{ret[s]=class_l[2];}
}
return (ret);
}
