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
#include <exception>
#include "checkinter.h"

using namespace Rcpp;
using namespace std;
template<class InputIterator, class T>

int myfind (InputIterator first, InputIterator last, const T& val)
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

RcppExport SEXP trans2adjacency(SEXP geneid1,SEXP geneid2,SEXP bf,SEXP genelist,SEXP min,SEXP max,SEXP bfthr)
{
 try
 {
 Rcpp::CharacterVector gene1(geneid1);
 Rcpp::CharacterVector gene2(geneid2); 
 Rcpp::NumericVector bv(bf);
 Rcpp::CharacterVector list(genelist); 
 double minn=Rcpp::as<double>(min);
 double maxx= Rcpp::as<double>(max);
 double thresbf= Rcpp::as<double>(bfthr);
 int i;
 string temp; 
 
 int findit1;
 int findit2;
 double bfactor;
 double ad;
 Rcpp::NumericMatrix adjacencies(list.size(), list.size());

 for (i=0;i<gene1.size();i++)
 {
  if(i%100==0){ checkInterrupt();}
  temp=gene1(i);
  findit1  = myfind(list.begin(),list.end(),temp);
  temp=gene2(i);
  findit2  = myfind(list.begin(),list.end(),temp);
  bfactor=bv(i);
  if(bfactor<=1||bfactor<=thresbf)
  {
  adjacencies(findit1,findit2)=0;
  adjacencies(findit2,findit1)=0;
  }else{
  ad=(log10(bfactor)-minn)/(maxx-minn);
  adjacencies(findit1,findit2)=ad;
  adjacencies(findit2,findit1)=ad;
  }
 }
  return Rcpp::List::create(Rcpp::Named("adjacency")=adjacencies);
 }
  catch (std::exception & e)
  {
   Rprintf("Error in (compiled code) trans2matrix: ");
   forward_exception_to_r(e);
  }
}
