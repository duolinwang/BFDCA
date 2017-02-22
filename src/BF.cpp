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
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <gsl/gsl_cdf.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include<iostream> 
#include<string> 
#include<iostream> 
#include<string> 
#include <exception>
#include "checkinter.h"

using namespace Rcpp;
using namespace std;

double KAPPA0, MU0, NU0;

struct bfPairModel{
long double bf;
int bftype;
} tmo;


double logLikDisc( gsl_matrix *D){
	return 0;
}


void cov_calculate(gsl_matrix *r, gsl_matrix * m)
{
  gsl_vector_view a, b;
  size_t i, j;
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size1; j++) {
      double v;
      a = gsl_matrix_row (m, i);
      b = gsl_matrix_row (m, j);
      v = gsl_stats_covariance (a.vector.data, a.vector.stride,b.vector.data, b.vector.stride, a.vector.size);
      gsl_matrix_set(r, i, j, v);
    }
  }
}

double logLikCont( gsl_matrix *C){
	int i,j;
	int nObs = C->size2;
	int nVar = C->size1;
	double v0 = nVar+2;//nVar+2;//NU0;
	
	double kappa = 0.01;//KAPPA0;
	double mu ;//= MU0;
	double P;
	double temp;
        gsl_vector *v=gsl_vector_calloc(nObs);	
	if(nVar == 0) return 0;
	// one variate
	else if(nVar == 1) {
                gsl_matrix_get_row (v, C,0);
	        mu=gsl_stats_mean((*v).data,1,nObs);
	
		double sigma = pow(gsl_stats_variance((*v).data,1,nObs),0.5);
		//int sigma = 1;			
		double v0Var = v0*pow((double)sigma,2);
	
			v0Var += (nObs-1)*gsl_stats_variance((*v).data,1,nObs);	
			v0Var += kappa*nObs/(kappa+nObs) * pow((gsl_stats_mean((*v).data,1,nObs) - mu),2);
		// calculate log-likelihood with one col of data
		P = -1*log(M_PI) * (double)nObs/2 + log(kappa/(kappa+nObs))/2 + gsl_sf_lngamma((v0+nObs)/2);
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");
		}
		P += -1*gsl_sf_lngamma(v0/2) + log(v0*pow((double)sigma,2)) * v0/2 - log(v0Var) * (v0+nObs)/2;
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");			
		}
		
	} 
	// multiple variates
	else {
	
		gsl_vector *rv1=gsl_vector_calloc(nObs);
		gsl_vector *rv2=gsl_vector_calloc(nObs);

		gsl_matrix *lambda=gsl_matrix_calloc(nVar, nVar);
                cov_calculate(lambda,C);
              
		double *means = new double[nVar];
		for(i=0; i<nVar; i++){	
			gsl_matrix_get_row(rv1, C, i);
			if(nObs == 1) 
				means[i] = 0;//gsl_vector_get(rv1,0) - mu;
			else
			{means[i] = 0;			
			}			        
		}
		gsl_matrix *lambdaN=gsl_matrix_calloc(nVar, nVar);
		for(i=0; i<nVar; i++) {
			for(j=0; j<nVar; j++) {				
				
				gsl_matrix_set(lambdaN,i,j,((gsl_matrix_get(lambda,i,j) + kappa*nObs/(kappa+nObs) * means[i] * means[j])));
				
				if(nObs > 1) {
			                gsl_matrix_get_row(rv1, C, i);
			                gsl_matrix_get_row(rv2, C, j);					
					temp=gsl_matrix_get(lambdaN,i,j)+(nObs-1)*gsl_stats_covariance((*rv1).data,1,(*rv2).data,1,nObs);
				        gsl_matrix_set(lambdaN,i,j,temp);
				}
			}		
		}       
		
		// steps for calculating determinants of lambda and lambdaN
		int s1,s2;
		gsl_permutation *p1 = gsl_permutation_alloc(nVar);
		gsl_permutation *p2 = gsl_permutation_alloc(nVar);
		gsl_linalg_LU_decomp(lambda, p1, &s1);
		gsl_linalg_LU_decomp(lambdaN, p2, &s2);
		
		// calculate log-likelihood
		P = -log(M_PI) * (double)nObs*nVar/2 + log(kappa/(kappa+nObs)) * (double)nVar/2;
		if(isnan(P)){

                throw std::range_error("error in logLikCont: p is NA.");
		
	       	}
		P += log(gsl_linalg_LU_det(lambda, s1)) * v0/2;

		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");			
			}
		P -= log(gsl_linalg_LU_det(lambdaN, s2)) * (v0+nObs)/2;
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");	 
		}
		
		for(i=0; i<nVar; i++){
			P += gsl_sf_lngamma((v0+nObs)/2 - (double)i/2) - gsl_sf_lngamma(v0/2 - (double)i/2);
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");
		}
		}
		
		delete [] means;
		gsl_vector_free(v);		
		gsl_vector_free(rv1);
		gsl_vector_free(rv2);
		gsl_matrix_free(lambda);
                gsl_matrix_free(lambdaN); 
                gsl_permutation_free(p1);
                gsl_permutation_free(p2);
	}
	
	return P;
}

double logLikCont_vector( gsl_vector *v){
	
	int nObs = v->size;
	double v0 = 1+2;//1+2;1;//NU0;
	
	double kappa = 0.01;//KAPPA0;
	double mu;//MU0;
	double P;
                mu=gsl_stats_mean((*v).data,1,nObs);
		double sigma = pow(gsl_stats_variance((*v).data,1,nObs),0.5);
		double v0Var = v0*pow((double)sigma,2);
	
			v0Var += (nObs-1)*gsl_stats_variance((*v).data,1,nObs);	
			v0Var += kappa*nObs/(kappa+nObs) * pow((gsl_stats_mean((*v).data,1,nObs) - mu),2);
		// calculate log-likelihood with one col of data
		P = -1*log(M_PI) * (double)nObs/2 + log(kappa/(kappa+nObs))/2 + gsl_sf_lngamma((v0+nObs)/2);
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");			
			}
		P += -1*gsl_sf_lngamma(v0/2) + log(v0*pow((double)sigma,2)) * v0/2 - log(v0Var) * (v0+nObs)/2;
		if(isnan(P)){
                throw std::range_error("error in logLikCont: p is NA.");
		}
	return P;
}


struct bfPairModel * bf1_pair_compute( vector<int> &D,  gsl_matrix *C,int inx1,int inx2,double bfthres)//return bf1
{
try{	int nObs = C->size2;
	int nVarC = 2;//C->size1;
	int j,k;
	long double p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,p_up,p_down;
	long double p_c1_mvn=0,p_c2_mvn=0,p_c1_nmvn=0,p_c2_nmvn=0; // p4 p5 have no common part with others.
	vector<int> row_c1;
	vector<int> row_c2;
	gsl_vector *subC1g1,*subC1g2;
	gsl_vector *subC2g1,*subC2g2;
	gsl_vector *subg1,*subg2;
	gsl_matrix *subC,*subC1,*subC2;
	
	gsl_vector *index=gsl_vector_alloc(2);
        gsl_vector_set(index,0,inx1);
	gsl_vector_set(index,1,inx2);
	
	struct bfPairModel *tm;
	for(j=0; j<nObs; j++){
			if(D[j] == 1 ) row_c1.push_back(j);
			else  row_c2.push_back(j);
                }

	        subC1=gsl_matrix_alloc(nVarC, row_c1.size());
	        subC2=gsl_matrix_alloc(nVarC, row_c2.size());
	        subC1g1=gsl_vector_alloc(row_c1.size());
	        subC1g2=gsl_vector_alloc(row_c1.size());
	        subC2g1=gsl_vector_alloc(row_c2.size());
	        subC2g2=gsl_vector_alloc(row_c2.size());
		subg1=gsl_vector_alloc(nObs);
		subg2=gsl_vector_alloc(nObs);
	        subC=gsl_matrix_alloc(nVarC, nObs);
                
		for(k=0; k<nVarC; k++)
		{
		  //cout<<"inx="<<index[k]<<endl;		  
		  for(j=0; j<row_c1.size(); j++)
		  {
			  gsl_matrix_set(subC1,k,j,gsl_matrix_get(C,gsl_vector_get(index,k), row_c1[j])); 
			  if(k==0){
				  gsl_vector_set(subC1g1,j,gsl_matrix_get(C,gsl_vector_get(index,k),row_c1[j]));
			  }else{
				  gsl_vector_set(subC1g2,j,gsl_matrix_get(C,gsl_vector_get(index,k),row_c1[j]));
			  }
		   }
		  for(j=0; j<row_c2.size(); j++)
		  {
			  gsl_matrix_set(subC2,k,j,gsl_matrix_get(C,gsl_vector_get(index,k), row_c2[j]));
		          if(k==0){
				  gsl_vector_set(subC2g1,j,gsl_matrix_get(C,gsl_vector_get(index,k),row_c2[j]));
			  }else{
				  gsl_vector_set(subC2g2,j,gsl_matrix_get(C,gsl_vector_get(index,k),row_c2[j]));
			  }
		  }
		  //cout<<endl;
		}

	
		gsl_matrix_get_row (subg1,C,gsl_vector_get(index,0));
		gsl_matrix_get_row (subg2,C,gsl_vector_get(index,1));
		gsl_matrix_set_row(subC,0,subg1);
		gsl_matrix_set_row(subC,1,subg2);

		p_c1_mvn = logLikCont(subC1);
		p_c2_mvn = logLikCont(subC2);
                p_c1_nmvn = logLikCont_vector(subC1g1)+ logLikCont_vector(subC1g2);
		p_c2_nmvn = logLikCont_vector(subC2g1)+ logLikCont_vector(subC2g2);
                p1 = p_c1_mvn+p_c2_mvn;
                p2 = p_c1_mvn+p_c2_nmvn;
                p3 = p_c1_nmvn+p_c2_mvn;
                p4 = logLikCont(subC);
                p5 = logLikCont_vector(subg1)+logLikCont_vector(subg2);// conver vecotr into matrix 
                p6 = p_c1_nmvn+p_c2_nmvn;
 p_up=p1+log(1+(long double)exp(p2-p1)+(long double)exp(p3-p1));
 if(isnan(p_up))p_up=max(max(p1,p2),p3);		 
 p_down=p4+log(1+(long double)exp(p5-p4)+(long double)exp(p6-p4));
 if(isnan(p_down))p_down=max(max(p4,p5),p6);
 	      
 long double bf=(long double)2*(p_up-p_down);
 if(bf>bfthres)
 {
   if(p1>p2 && p1>p3){
         tmo.bftype=1;
	 tmo.bf=bf;
   }
   else if(p2>p1 && p2>p3){
         tmo.bftype=2;
	 tmo.bf=bf;
   }else if(p3>p1 && p3>p2){
         tmo.bftype=3;
	 tmo.bf=bf;
   }
 }else{
 tmo.bftype=0;
 tmo.bf=bf;
 }
 tm=&tmo;
  row_c1.clear();
		gsl_matrix_free(subC1);
		gsl_matrix_free(subC2);
		gsl_matrix_free(subC);
		gsl_vector_free(subC1g1);
		gsl_vector_free(subC1g2);
		gsl_vector_free(subC2g1);
		gsl_vector_free(subC2g2);
		gsl_vector_free(subg1);
		gsl_vector_free(subg2);
		row_c2.clear();
		gsl_vector_free(index);
vector<int>().swap( row_c1 );
vector<int>().swap( row_c2 );
		return tm;

}catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) bf1_pair_compute: ");
   forward_exception_to_r(e);
  }
}

RcppExport SEXP BF1_pair(SEXP data,SEXP group,SEXP thres,SEXP echo)
{
try{
	int i,j,t,nObs,varC;
        Rcpp::NumericMatrix dataC(data);
        Rcpp::NumericVector groupD(group);
	int ec=Rcpp::as<int>(echo);
        double bfthres  = Rcpp::as<double>(thres);
	vector<int> gene1id_vector;
	vector<int> gene2id_vector;
        vector<double> bf_vector;
	vector<int> bf_type;
	gsl_vector *wCC;
	nObs = dataC.nrow();
	varC = dataC.ncol();
	gsl_matrix *C=gsl_matrix_alloc(varC, nObs);

		for(i = 0; i < varC; i++)
		{
		  for(j=0;j<nObs;j++)
		  {
	           gsl_matrix_set(C,i,j,dataC(j,i));
		  }
		}
                  vector<int> D;
                  for(j=0;j<nObs;j++)
		  {
	           D.push_back(groupD(j));
		  }
	// define hyper-parameters
	KAPPA0 = 1;
	NU0 = varC + 1;
	double *means = new double[varC];
	gsl_vector *tmp=gsl_vector_alloc (nObs);
	for(i=0; i<varC; i++){
		gsl_matrix_get_row(tmp, C, i);
		double *r1 = new double[nObs];
               	for(j=0; j<nObs; j++)
		r1[j] = (double)gsl_vector_get(tmp,j);//row1(j);	
		means[i] =gsl_stats_mean(r1,1,nObs);
	        delete r1;
	}
	delete [] means;
	gsl_vector_free(tmp);	      

       vector<string> w;  
       string str;
       
       int s; 
       int passnum=0;
       double bf_large=0;
       int pair_num=1;
       double pvalue=0;
       struct bfPairModel * bf;
       //long double bf=0;
   
       for(i=0;i<varC-1;i++)    
       {
          for(j=i+1;j<varC;j++)
	  {
	   pair_num+=1;
	   if(pair_num%5000==0){
	   if(ec==1)Rcpp::Rcout<<pair_num<<" pairs ready ..."<<std::endl;
           checkInterrupt();
	   }
	   
	   bf=bf1_pair_compute(D,C,i,j,bfthres);
	   //bf=bf1_pair_compute(D,C,i,j);
	   
	   //if(bf>=bfthres)
	   if((*bf).bf>=bfthres)
	   {
                   gene1id_vector.push_back(i+1);
                   gene2id_vector.push_back(j+1);
                   bf_vector.push_back((*bf).bf); 
                   //bf_vector.push_back(6);
		   //bf_vector.push_back(bf);
		   bf_type.push_back((*bf).bftype);
		   //bf_type.push_back(1);
		   passnum+=1;
	           //free(bf);	   
           }//bftype
         }//j
       }//i

	      gsl_matrix_free(C);
	     
              SEXP Rgeneid1= PROTECT(Rf_allocVector(INTSXP, passnum));
              SEXP Rgeneid2= PROTECT(Rf_allocVector(INTSXP, passnum));
              SEXP Rbf= PROTECT(Rf_allocVector(REALSXP, passnum));
              SEXP Rbftype= PROTECT(Rf_allocVector(REALSXP, passnum));
	      
	      for(i=0;i<passnum;i++)
	      {
               INTEGER(Rgeneid1)[i]=gene1id_vector[i];
               INTEGER(Rgeneid2)[i]=gene2id_vector[i];
               REAL(Rbf)[i]=bf_vector[i];
	       REAL(Rbftype)[i]=bf_type[i];
              }
	      SEXP vec = PROTECT(Rf_allocVector(VECSXP, 4));
              SET_VECTOR_ELT(vec, 0, Rgeneid1);
              SET_VECTOR_ELT(vec, 1, Rgeneid2);
              SET_VECTOR_ELT(vec, 2, Rbf);
              SET_VECTOR_ELT(vec, 3, Rbftype);
              UNPROTECT(5);
	      return (vec);
  }
  catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) BF1_pair: ");
   forward_exception_to_r(e);
  }

}

struct trainmodule{
int Mselect;
long double Mbf;
gsl_vector * MmeanC1;
gsl_vector * MmeanC2;

union M1{
gsl_matrix * McovC1;
gsl_vector * MvarC1;
} Mc1;

union M2{
gsl_matrix * McovC2;
gsl_vector * MvarC2;
} Mc2;

};



struct trainmodule *bf1_module_compute_fortrain( vector<int> &D,  gsl_matrix *C,double bfthres)//return bf1
{
try{
	int nObs = C->size2;
	int nVarC = C->size1;//2;//C->size1;
	int i,j,k;
	long double p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,p_up,p_down;
	long double p_c1_mvn=0,p_c2_mvn=0,p_c1_nmvn=0,p_c2_nmvn=0; 
	vector<int> row_c1;
	vector<int> row_c2;
	gsl_vector *subC1g;//*subC1g2;
	gsl_vector *subC2g;//*subC2g2;
	gsl_vector *subg1;//*subg2;
	gsl_matrix *subC1,*subC2;

	for(j=0; j<nObs; j++){
			if(D[j] == 1 ) row_c1.push_back(j);
			else  row_c2.push_back(j);
                }

	        subC1=gsl_matrix_alloc(nVarC, row_c1.size());
	        subC2=gsl_matrix_alloc(nVarC, row_c2.size());
	       
	       	subC1g=gsl_vector_alloc(row_c1.size());
	        subC2g=gsl_vector_alloc(row_c2.size());
		subg1=gsl_vector_alloc(nObs);
		//subg2=gsl_vector_alloc(nObs);
		for(k=0; k<nVarC; k++)
		{
		  for(j=0; j<row_c1.size(); j++)
		  {
			 gsl_matrix_set(subC1,k,j,gsl_matrix_get(C,k, row_c1[j]));  
		         gsl_vector_set(subC1g,j,gsl_matrix_get(C,k,row_c1[j]));
		  }
	          p_c1_nmvn += logLikCont_vector(subC1g);
		  for(j=0; j<row_c2.size(); j++)
		  {
			  gsl_matrix_set(subC2,k,j,gsl_matrix_get(C,k, row_c2[j]));
	         	  gsl_vector_set(subC2g,j,gsl_matrix_get(C,k,row_c2[j]));
		  }
	          p_c2_nmvn += logLikCont_vector(subC2g);
		}
		for(k=0; k<nVarC; k++)
		{
		gsl_matrix_get_row(subg1,C,k);
		p5+= logLikCont_vector(subg1);
		}
		 
		p_c1_mvn = logLikCont(subC1);
		p_c2_mvn = logLikCont(subC2);
                p1 = p_c1_mvn+p_c2_mvn;
                p2 = p_c1_mvn+p_c2_nmvn;
                p3 = p_c1_nmvn+p_c2_mvn;
		p4 = logLikCont(C);
                p6 = p_c1_nmvn+p_c2_nmvn;
	
 p_up=p1+log(1+(long double)exp(p2-p1)+(long double)exp(p3-p1));
 if(isnan(p_up))p_up=max(max(p1,p2),p3);		 
 p_down=p4+log(1+(long double)exp(p5-p4)+(long double)exp(p6-p4));
 if(isnan(p_down))p_down=max(max(p4,p5),p6);
 
 struct trainmodule *tm;
 static struct trainmodule tmo;
 double bf=(long double)2*(p_up-p_down);
 gsl_matrix * cov1=gsl_matrix_alloc(nVarC, nVarC);
 gsl_matrix * cov2=gsl_matrix_alloc(nVarC, nVarC);
 double sigma1;
 double sigma2;
 double sigma;
 double v0;
 if(bf>bfthres)
 {
 tmo.Mbf=bf;
 tmo.MmeanC1=gsl_vector_alloc(nVarC);
 tmo.MmeanC2=gsl_vector_alloc(nVarC);
 if(p1>p2 && p1>p3){
       tmo.Mselect=1;
       v0=nVarC+2;
       tmo.Mc1.McovC1=gsl_matrix_alloc(nVarC, nVarC);
       tmo.Mc2.McovC2=gsl_matrix_alloc(nVarC, nVarC);

       for(k=0; k<nVarC; k++)
       {
        for(j=0; j<row_c1.size(); j++) gsl_vector_set(subC1g,j,gsl_matrix_get(C,k,row_c1[j]));	
        for(j=0; j<row_c2.size(); j++) gsl_vector_set(subC2g,j,gsl_matrix_get(C,k,row_c2[j]));	
        gsl_vector_set(tmo.MmeanC1,k,gsl_stats_mean((*subC1g).data,1,row_c1.size()));
        gsl_vector_set(tmo.MmeanC2,k,gsl_stats_mean((*subC2g).data,1,row_c2.size()));
       }
       cov_calculate(cov1,subC1);
       cov_calculate(cov2,subC2);

    for (i = 0; i < nVarC; i++) {
     for (j = 0; j < nVarC; j++){
       gsl_matrix_set(tmo.Mc1.McovC1, i, j,row_c1.size()/(v0+row_c1.size()+nVarC+2)*gsl_matrix_get(cov1,i,j));
       gsl_matrix_set(tmo.Mc2.McovC2, i, j, row_c2.size()/(v0+row_c2.size()+nVarC+2)*gsl_matrix_get(cov2,i,j));
    }}

 }
 else if(p2>p1 && p2>p3){
	 tmo.Mselect=2;
         tmo.Mc1.McovC1=gsl_matrix_alloc(nVarC, nVarC);
         tmo.Mc2.MvarC2=gsl_vector_alloc(nVarC);
	 v0=nVarC+2;
         for(k=0; k<nVarC; k++)
         {
            for(j=0; j<row_c1.size(); j++) gsl_vector_set(subC1g,j,gsl_matrix_get(C,k,row_c1[j]));	
            for(j=0; j<row_c2.size(); j++) gsl_vector_set(subC2g,j,gsl_matrix_get(C,k,row_c2[j]));	
            gsl_vector_set(tmo.MmeanC1,k,gsl_stats_mean((*subC1g).data,1,row_c1.size()));
            gsl_vector_set(tmo.MmeanC2,k,gsl_stats_mean((*subC2g).data,1,row_c2.size()));
            sigma=gsl_stats_variance((*subC2g).data,1,(*subC2g).size);
	    gsl_vector_set(tmo.Mc2.MvarC2,k,sigma);
         }
    cov_calculate(cov1,subC1);
    for (i = 0; i < nVarC; i++) {
     for (j = 0; j < nVarC; j++){
       gsl_matrix_set(tmo.Mc1.McovC1, i, j,row_c1.size()/(v0+row_c1.size()+nVarC+2)*gsl_matrix_get(cov1,i,j));
    }}
 }
 else if(p3>p1 && p3>p2){
       tmo.Mselect=3;
       tmo.Mc1.MvarC1=gsl_vector_alloc(nVarC);
       tmo.Mc2.McovC2=gsl_matrix_alloc(nVarC, nVarC);
       v0=nVarC+2;
       for(k=0; k<nVarC; k++)
       {
           for(j=0; j<row_c1.size(); j++) gsl_vector_set(subC1g,j,gsl_matrix_get(C,k,row_c1[j]));	
           for(j=0; j<row_c2.size(); j++) gsl_vector_set(subC2g,j,gsl_matrix_get(C,k,row_c2[j]));	
           gsl_vector_set(tmo.MmeanC1,k,gsl_stats_mean((*subC1g).data,1,row_c1.size()));
           gsl_vector_set(tmo.MmeanC2,k,gsl_stats_mean((*subC2g).data,1,row_c2.size()));
          sigma=gsl_stats_variance((*subC1g).data,1,(*subC1g).size);
	    gsl_vector_set(tmo.Mc1.MvarC1,k,sigma);
       }   

    cov_calculate(cov2,subC2);
    for (i = 0; i < nVarC; i++) {
     for (j = 0; j < nVarC; j++){
       gsl_matrix_set(tmo.Mc2.McovC2, i, j,row_c2.size()/(v0+row_c2.size()+nVarC+2)*gsl_matrix_get(cov2,i,j));
    }
   }
  }
}else{tmo.Mselect=0;}

                tm=&tmo;
	        row_c1.clear();
		gsl_matrix_free(subC1);
		gsl_matrix_free(subC2);
		gsl_vector_free(subC1g);
		gsl_vector_free(subC2g);
		gsl_vector_free(subg1);
		gsl_matrix_free(cov1);
		gsl_matrix_free(cov2);
		row_c2.clear();
return(tm);
}
catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) bf1_module_compute_fortrain: ");
   forward_exception_to_r(e);
  }
}

RcppExport SEXP BF_train(SEXP data,SEXP group,SEXP thres)
{
try{
	int i,j,t,nObs,varC;
        Rcpp::NumericMatrix dataC(data);
        Rcpp::NumericVector groupD(group);
        double bfthres  = Rcpp::as<double>(thres);
	vector<int> gene1id_vector;
	vector<int> gene2id_vector;
        vector<double> bf_vector;
	vector<double> pvalue_vector;
	gsl_vector *wCC;
	nObs = dataC.nrow();
	varC = dataC.ncol();
	gsl_matrix *C=gsl_matrix_alloc(varC, nObs);

		for(i = 0; i < varC; i++)
		{
		  for(j=0;j<nObs;j++)
		  {
	           gsl_matrix_set(C,i,j,dataC(j,i));
		  }
		}
                  vector<int> D;
                  for(j=0;j<nObs;j++)
		  {
	           D.push_back(groupD(j));
		  }

	      struct trainmodule * bf;
	      bf=bf1_module_compute_fortrain(D,C,bfthres);
              Rcpp:NumericVector  type(1);
	      type[0]=(*bf).Mselect;
	      if((*bf).Mselect==0)
	      {	
              SEXP vec = PROTECT(Rf_allocVector(VECSXP, 1));
              SET_VECTOR_ELT(vec, 0, type);
              UNPROTECT(1);
	      return (vec);
	      }else
	      {      
	        Rcpp::NumericVector meanc1(C->size1);
	        Rcpp::NumericVector meanc2(C->size1);
		Rcpp::NumericVector mbf(1);
		mbf[0]=(*bf).Mbf;
	      if((*bf).Mselect==1)
	      {
	        Rcpp::NumericMatrix mcovc1(C->size1,C->size1);
		Rcpp::NumericMatrix mcovc2(C->size1,C->size1);
	        for(i=0;i<C->size1;i++)
	        {
                meanc1(i)=gsl_vector_get((*bf).MmeanC1,i);
                meanc2(i)=gsl_vector_get((*bf).MmeanC2,i);
	        for(j=0;j<C->size1;j++)
	        {
                  mcovc1(i,j)=gsl_matrix_get((*bf).Mc1.McovC1,i,j);
                  mcovc2(i,j)=gsl_matrix_get((*bf).Mc2.McovC2,i,j);
                }
		}
              
	      SEXP vec = PROTECT(Rf_allocVector(VECSXP, 6));
              SET_VECTOR_ELT(vec, 0, type);
              SET_VECTOR_ELT(vec, 1, mbf);
              SET_VECTOR_ELT(vec, 2, meanc1);
              SET_VECTOR_ELT(vec, 3, meanc2);
              SET_VECTOR_ELT(vec, 4, mcovc1);
              SET_VECTOR_ELT(vec, 5, mcovc2);
              UNPROTECT(1);
	      return (vec);
	      }

              if((*bf).Mselect==2)
	      {
		Rcpp::NumericMatrix mcovc1(C->size1,C->size1);
		Rcpp::NumericVector mvarc2(C->size1);
	        for(i=0;i<C->size1;i++)
	        {
                meanc1(i)=gsl_vector_get((*bf).MmeanC1,i);
                meanc2(i)=gsl_vector_get((*bf).MmeanC2,i);
                mvarc2(i)=gsl_vector_get((*bf).Mc2.MvarC2,i);
	        for(j=0;j<C->size1;j++)
	        {
                  mcovc1(i,j)=gsl_matrix_get((*bf).Mc1.McovC1,i,j);
                }
		}
              SEXP vec = PROTECT(Rf_allocVector(VECSXP, 6));
              SET_VECTOR_ELT(vec, 0, type);
              SET_VECTOR_ELT(vec, 1, mbf);
              SET_VECTOR_ELT(vec, 2, meanc1);
              SET_VECTOR_ELT(vec, 3, meanc2);
              SET_VECTOR_ELT(vec, 4, mcovc1);
              SET_VECTOR_ELT(vec, 5, mvarc2);
              UNPROTECT(1);
	      return (vec);
	      }

              if((*bf).Mselect==3)
	      {
		Rcpp::NumericVector mvarc1(C->size1);
		Rcpp::NumericMatrix mcovc2(C->size1,C->size1);
	        for(i=0;i<C->size1;i++)
	        {
                meanc1(i)=gsl_vector_get((*bf).MmeanC1,i);
                meanc2(i)=gsl_vector_get((*bf).MmeanC2,i);
                mvarc1(i)=gsl_vector_get((*bf).Mc1.MvarC1,i);
	        for(j=0;j<C->size1;j++)
	        {
                  mcovc2(i,j)=gsl_matrix_get((*bf).Mc2.McovC2,i,j);
                }
	      }
              SEXP vec = PROTECT(Rf_allocVector(VECSXP, 6));
              SET_VECTOR_ELT(vec, 0, type);
              SET_VECTOR_ELT(vec, 1, mbf);
              SET_VECTOR_ELT(vec, 2, meanc1);
              SET_VECTOR_ELT(vec, 3, meanc2);
              SET_VECTOR_ELT(vec, 4, mvarc1);
              SET_VECTOR_ELT(vec, 5, mcovc2);
              UNPROTECT(1);
	      return (vec);
	      }	
	      }
}catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) BF_train: ");
   forward_exception_to_r(e);
  }
}
