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

BFplot2files <- function(dataExp,class,classlabel,edgefilename,plotfilename="Plotedges.pdf"){

    if (is.null(dim(dataExp))) { 
		stop("argument \"dataExp\" is missing,with no default. It is a matrix or data frame containing expression data."); 
	}
   if(!("matrix" %in% class(dataExp)|"data.frame" %in% class(dataExp))){stop("argument \"dataExp\" must be a matrix or data frame containing expression data.");} 

if(is.null(length(edgefilename))){stop("Specify a correct file name in argument \"edgefilename\" for edges to be ploted.");}
if(is.null(length(plotfilename))){stop("Specify a correct file name in argument \"plotfilename\" to output plots. For example: Plotedges.pdf. ");}
if(length(class)!=dim(dataExp)[1]){stop("argument \"dataExp\" and \"class\" are not matching. They must have same length. ");}
if(length(table(class))>2){stop("by far package BFDCA can only accept binary conditions, but argumenet \"class\" contains more than 2 classes. ")}   
dataExp<-as.matrix(dataExp);
if(length(classlabel)>2){stop(" by far package BFDCA can only accept binary conditions, but classlabel contains more than 2 characters. ")}

read.table(edgefilename,header=FALSE, sep="\t")->selectpair;

pdf(plotfilename);
for(i in c(1:dim(selectpair)[1]))
{
g1<- as.character(selectpair[i,1]);
g2<-as.character(selectpair[i,2]);
BFplot(dataExp,class,classlabel,g1,g2);
}
dev.off();

}
