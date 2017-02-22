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

BFplot <- function(dataExp,class,classlabel,gene1,gene2){
   if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, it must be a matrix or data frame containing expression data."); 
	}
   if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" must be a matrix or data frame containing expression data.");} 

   if(length(class)!=dim(dataExp)[1]){stop("argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
   if(length(table(class))>2){stop("by far package BFDCA can only accept binary conditions, but argumenet \"class\" contains more than 2 classes. ")}   
   dataExp<-as.matrix(dataExp);
   if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}      
   if(is.null(colnames(dataExp)))
   {
   colnames(dataExp)<-c(1:dim(dataExp)[2])
   }  
   
   if(!(class(gene1)=="character" | class(gene1) == "numeric" | class(gene1) =="integer")){stop(paste("argument gene1 is not appropriate, it should be a numeric or character indicates gene indexes or gene ids for gene1 . ",sep=""))}  
   
   if(!(class(gene2)=="character" | class(gene2) == "numeric" | class(gene2) =="integer")){stop(paste("argument gene2 is not appropriate, it should be a numeric or character indicates gene indexes or gene ids for gene2 . ",sep=""))} 
 
   if(class(gene1)=="character"){if(!gene1%in%colnames(dataExp)){stop(paste("gene: ",gene1," is not in dataExp. ",sep=""))}}
   if(class(gene2)=="character"){if(!gene2%in%colnames(dataExp)){stop(paste("gene: ",gene2," is not in dataExp. ",sep=""))}}

   if(class(gene1)=="numeric"|class(gene1)=="integer"){if(!gene1%in%c(1:dim(dataExp)[2])){stop(paste("gene index: ",gene1," is not in dataExp. ",sep=""))}}
   if(class(gene2)=="numeric"|class(gene2)=="integer"){if(!gene2%in%c(1:dim(dataExp)[2])){stop(paste("gene index: ",gene2," is not in dataExp. ",sep=""))}}

   class_l<-array(NA,2);
   class_l[class==classlabel[1]]=1;
   class_l[class==classlabel[2]]=2;
   
   n=100;
   k=sqrt(6);

   #condition1   
   meanc1g1<-mean(dataExp[class_l==1,gene1]);
   meanc1g2<-mean(dataExp[class_l==1,gene2]);
   sigmac1<-cov(dataExp[class_l==1,c(gene1,gene2)]);
   evc1<-eigen(sigmac1);
   t=seq(0,2*pi,length.out=n);
   xy=rbind(cos(t),sin(t));
   w= (k*evc1$vectors%*%sqrt(matrix(c(evc1$values[1],0,0,evc1$values[2]),ncol=2)))%*%xy;
   xc1=meanc1g1+w[1,];
   yc1=meanc1g2+w[2,];
  
   #condition2
   meanc2g1<-mean(dataExp[class_l==2,gene1]);
   meanc2g2<-mean(dataExp[class_l==2,gene2]);
   sigmac2<-cov(dataExp[class_l==2,c(gene1,gene2)]);
   evc2<-eigen(sigmac2);
   w= (k*evc2$vectors%*%sqrt(matrix(c(evc2$values[1],0,0,evc2$values[2]),ncol=2)))%*%xy;
   xc2=meanc2g1+w[1,];
   yc2=meanc2g2+w[2,];

   xmin=min(c(xc1,xc2));
   xmax=max(c(xc1,xc2));
   ymin=min(c(yc1,yc2));
   ymax=max(c(yc1,yc2));
   
     
   plot(xc1,yc1,type="l",col="blue",lwd=2,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=gene1,ylab=gene2,main=paste(gene1,"~",gene2," ",classlabel[1],"=blue ",classlabel[2],"=red",sep=""));
   points(dataExp[class_l==1,gene1],dataExp[class_l==1,gene2],col="blue",pch=4);
   lines(xc2,yc2,type="l",col="red",lwd=2);

   points(dataExp[class_l==2,gene1],dataExp[class_l==2,gene2],col="red");
   points(meanc1g1,meanc1g2,pch=4,cex=2,lwd=3,col="blue");
   points(meanc2g1,meanc2g2,pch=4,cex=2,lwd=3,col="red"); 
   #legend("topleft",legend=classlabel,pch=c(4,1),col=c("blue","red"))
}                     
