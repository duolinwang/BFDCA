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

BF_similarity <- function(bfmatrix,softPower=1,bfthr=6,keepedges=0){


if(bfthr<0){warning("argument \"bfthr\" is not set properly. It should be set largher than 0. ");}
bfmatrix$bf.value[which(bfmatrix$bf.value<0)]<-0;
max<-log(max(bfmatrix$bf.value),10);
min<-max(0,log(min(bfmatrix$bf.value),10));
if(max<min){stop("max BF value cannot less than min BF value.");}
if(keepedges!=0){

 if(keepedges>dim(bfmatrix)[1]){stop("argument \"keepedges\" is not set appropriately. It's larger than the remained edges.")}
 bfthr=sort(bfmatrix[,3],decreasing = TRUE)[keepedges];
}

geneid1<-as.character(bfmatrix$geneid1);
geneid2<-as.character(bfmatrix$geneid2);
genelist<-unique(c(geneid1,geneid2));
trans2adjacency(geneid1,geneid2,bfmatrix$bf.value,genelist,min,max,bfthr)->x;
adjacencies<-x[[1]];
colnames(adjacencies)<-genelist;
rownames(adjacencies)<-genelist;
return(adjacencies);
}
