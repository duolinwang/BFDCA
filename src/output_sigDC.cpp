/*Copyright (C) 2016 Duolin Wang

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <gsl/gsl_statistics.h>
#include <exception>
#include "checkinter.h"

using namespace Rcpp;
using namespace std;

template<class InputIterator, class T>
int myfind2 (InputIterator first, InputIterator last, const T& val)
{
  int index=0;	
  while (first!=last) {
    if (*first==val) return index;
    {
      ++first;
      index++;
    }
  }
  return -1;
}

RcppExport SEXP output_sigDC(SEXP mst2_data,SEXP labels,SEXP weights,SEXP bf,SEXP genelist,SEXP bvthres,SEXP cutoff){

 try{
 Rcpp::NumericVector geneweight(weights); 
 Rcpp::NumericVector genegroup(labels); 
 Rcpp::NumericMatrix bvmatrix(bf); 
 Rcpp::CharacterMatrix mst2(mst2_data); 

 Rcpp::CharacterVector geneid(genelist); 

 double bv_thres=Rcpp::as<double>(bvthres);
 double weight_thres=Rcpp::as<double>(cutoff);
 Rcpp::CharacterVector geneid1,geneid2;
 Rcpp::NumericVector groupid;
 Rcpp::NumericVector gene1_weight,gene2_weight,bf_value; 
 int i,j;
 string g1,g2;
 double tempbv;
 int tmpg1,tmpg2;
 int edgenum=0;

  for(i=0;i<mst2.nrow()-1;i++)
  {
      edgenum++;
      if(edgenum%500==0)checkInterrupt();
      g1=mst2(i,0);
      g2=mst2(i,1);  
      tmpg1  = myfind2(geneid.begin(),geneid.end(),g1);
      tmpg2  = myfind2(geneid.begin(),geneid.end(),g2);
      tempbv=bvmatrix(tmpg1,tmpg2);
      if(tempbv>=bv_thres && geneweight(tmpg1)>=weight_thres && geneweight(tmpg2)>=weight_thres)
      {  
      bf_value.push_back(tempbv);
      geneid1.push_back(g1);
      geneid2.push_back(g2);
      groupid.push_back(genegroup(tmpg1));
      gene1_weight.push_back(geneweight(tmpg1));
      gene2_weight.push_back(geneweight(tmpg2));
      }   
    }//i

return Rcpp::DataFrame::create(Rcpp::Named("geneid1")=geneid1,Rcpp::Named("geneid2")=geneid2,Rcpp::Named("groupid")=groupid,Rcpp::Named("gene1_weight")=gene1_weight,Rcpp::Named("gene2_weight")=gene2_weight,Rcpp::Named("bf.value")=bf_value);
 }
  catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) output_sigDC: ");
   forward_exception_to_r(e);
  }
}
