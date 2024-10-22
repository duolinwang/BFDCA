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


RcppExport SEXP output_network(SEXP labels, SEXP weights,SEXP bf,SEXP type,SEXP genelist,SEXP datac1,SEXP datac2,SEXP medianc1,SEXP medianc2,SEXP gsdeflist,SEXP corthres,SEXP bvthres,SEXP echo){

 try{

 Rcpp::NumericVector genegroup(labels);
 Rcpp::NumericVector geneweight(weights); 
 Rcpp::NumericMatrix bvmatrix(bf); 
 Rcpp::NumericMatrix tymatrix(type); 

 Rcpp::CharacterVector geneid(genelist); 
 Rcpp::NumericMatrix data1(datac1); 
 Rcpp::NumericMatrix data2(datac2); 
 Rcpp::NumericVector median1(medianc1); 
 Rcpp::NumericVector median2(medianc2); 
 Rcpp::List gslist(gsdeflist);
 int ec=Rcpp::as<int>(echo);

 double cor_thres=Rcpp::as<double>(corthres);
 double bv_thres=Rcpp::as<double>(bvthres);


 Rcpp::CharacterVector geneid1,geneid2;
 Rcpp::NumericVector groupid;
 Rcpp::NumericVector gene1_weight,gene2_weight,bf_value,ty_value,pc1,pc2,median_g1_c1,median_g1_c2,median_g2_c1,median_g2_c2, median_foldchange1,median_foldchange2;
 Rcpp::CharacterVector correlation;
 

 int i,j,g1,g2;
 double tempbv;
 int tmpg1,tmpg2;
 string tmpgs1,tmpgs2;
 double cor1,cor2,fold_g1,fold_g2;
 vector<double> pdata1x(data1.nrow()),pdata1y(data1.nrow()), pdata2x(data2.nrow()),pdata2y(data2.nrow());
 int edgenum=0;
 for(i=0;i<gslist.size();i++)
 {
  Rcpp::NumericVector genes= Rcpp::as<Rcpp::NumericVector>(gslist(i));

  for(g1=0;g1<genes.size()-1;g1++)
  {
    for(g2=(g1+1);g2<genes.size();g2++)
    {
      edgenum++;      
      if(ec==1)if(edgenum%1000==0){Rcpp::Rcout<<edgenum<<" pairs ready..."<<std::endl;}
      if(edgenum%2000==0)checkInterrupt();

       
      tmpg1=genes(g1)-1;
      tmpg2=genes(g2)-1;
      tempbv=bvmatrix(tmpg1,tmpg2);
      if(tempbv>=bv_thres)
      {  
      bf_value.push_back(tempbv);
      ty_value.push_back(tymatrix(tmpg1,tmpg2));

      tmpgs1=geneid(tmpg1);
      tmpgs2=geneid(tmpg2);
      geneid1.push_back(tmpgs1);
      geneid2.push_back(tmpgs2);
      groupid.push_back(genegroup(tmpg1));
      gene1_weight.push_back(geneweight(tmpg1));
      gene2_weight.push_back(geneweight(tmpg2));
      
      median_g1_c1.push_back(median1(tmpg1));
      median_g1_c2.push_back(median2(tmpg1));
      median_g2_c1.push_back(median1(tmpg2));
      median_g2_c2.push_back(median2(tmpg2));
      
      fold_g1=median1(tmpg1)/median2(tmpg1);
      fold_g2=median1(tmpg2)/median2(tmpg2);
      median_foldchange1.push_back(fold_g1);
      median_foldchange2.push_back(fold_g2);
      
      for(j=0;j<data1.nrow();j++)
      {
       pdata1x[j]=data1(j,tmpg1);
       pdata1y[j]=data1(j,tmpg2);
      }
      cor1=  gsl_stats_covariance(pdata1x.data(), 1,pdata1y.data(), 1,data1.nrow())/(gsl_stats_sd(pdata1x.data(),1,data1.nrow())*gsl_stats_sd(pdata1y.data(),1,data1.nrow()));
      for(j=0;j<data2.nrow();j++)
      {
       pdata2x[j]=data2(j,tmpg1);
       pdata2y[j]=data2(j,tmpg2);
      }
      cor2=gsl_stats_covariance(pdata2x.data(), 1,pdata2y.data(), 1,data2.nrow())/(gsl_stats_sd(pdata2x.data(),1,data2.nrow())*gsl_stats_sd(pdata2y.data(),1,data2.nrow()));
      pc1.push_back(cor1);
      pc2.push_back(cor2);
      
      if(abs(cor1-cor2)>=cor_thres)
      {
      if(cor1>=cor2){correlation.push_back("down");}
      else{correlation.push_back("up");}
      }else{correlation.push_back("unchange");}


     }
     
         
     }//g1
    }//g2
  }//i

return Rcpp::DataFrame::create(Rcpp::Named("geneid1")=geneid1,Rcpp::Named("geneid2")=geneid2,Rcpp::Named("groupid")=groupid,Rcpp::Named("gene1_weight")=gene1_weight,Rcpp::Named("gene2_weight")=gene2_weight,Rcpp::Named("bf.value")=bf_value,Rcpp::Named("type")=ty_value,Rcpp::Named("pc1")=pc1,Rcpp::Named("pc2")=pc2,Rcpp::Named("cortypes")=correlation,Rcpp::Named("mean_g1_c1")=median_g1_c1,Rcpp::Named("mean_g1_c2")=median_g1_c2,Rcpp::Named("mean_foldchange_gene1")=median_foldchange1,Rcpp::Named("mean_g2_c1")=median_g2_c1,Rcpp::Named("mean_g2_c2")=median_g2_c2,Rcpp::Named("mean_foldchange_gene2")=median_foldchange2);
 }
  catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) output_network: ");
   forward_exception_to_r(e);
  }
}
