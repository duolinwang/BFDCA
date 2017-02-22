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

sigDCpair_SFS<- function(dataExp,class,classlabel,sigDC,DC_acc="DC_pair_acc.txt",bfthr=6,by=50,LOOCV=TRUE,valid=NULL,valid_class=NULL)
{

    if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing, with no default. It must be a matrix or data frame containing expression data."); 
	}

   if(is.null(length(class))){stop("argument \"class\" is missing, with no default.")}

   if(is.null(length(classlabel))){stop("argument \"classlabel\" is missing, with no default.");}

   if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" must be a matrix or data frame containing expression of training data.");}
   if(length(class)!=dim(dataExp)[1]){stop(" argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
   if(length(table(class))>2){stop(" by far package BFDCA can only accept binary conditions, but class contains more than 2 classes. ")}   
   if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}    
 
selectpair_all<-sigDC[order(sigDC[,3],decreasing = TRUE),];
score=seq.int(1,dim(selectpair_all)[1],by=by);

if(LOOCV==FALSE){
    if(is.null(valid)){stop("argument \"valid\" is null. \"valid\" must be data frame or matrix containing expression of validation data with compatible columns (order of columns and dimensions of columns) to \"dataExp\".")}
    if(!("matrix" %in% class(valid)|"data.frame" %in% class(valid))){stop("argument \"valid\" has incompatible type. \"valid\" must be data frame or matrix containing expression of validation data with compatible columns (order of columns and dimensions of columns) to \"dataExp\"");}
    if(dim(dataExp)[2]!=dim(valid)[2]){stop(" argument \"dataExp\" and \"valid\" are not matching. They must have compatible columns (order of columns and dimensions of columns) to \"dataExp\". ");}
   if(length(valid_class)!=dim(valid)[1]){stop(" argument \"valid\" and \"valid_class\" are not matching. They must have same length. ");}
}

if(is.null(colnames(dataExp)))
{
   colnames(dataExp)<-c(1:dim(dataExp)[2])
}

out<-tryCatch({file(DC_acc,"w")},error=function(cond){return(NA)});
if(!is.na(out))
{
for(time in c(1:length(score)))
{
selectgene<-NULL;
edges=selectpair_all[1:score[time],];
e<-1;
ac=0;
 
 if(LOOCV==TRUE)
 {
   valid=dataExp;
   for(s in c(1:length(class)))
   {
    BFtrain(dataExp[-s,],class[-s],classlabel,edges,bfthr=bfthr)->model;
    BFtest(valid[s,],model)->tclass;
    if(class[s]==tclass){ac=ac+1;}
   }
   ac=ac/length(class);
 }else{
    BFtrain(dataExp,class,classlabel,edges,bfthr=bfthr)->model;
    BFtest(valid,model)->tclass;
    ac=sum(tclass==valid_class)/length(valid_class);
 }
writeLines(paste("top_",score[time],"_edges","\t",as.character(ac),sep=""),out);
}
close(out);
}else{stop(paste("can't create file ",DC_acc,sep=""));}

}
