#ifndef __checkinter_h__

#define __checkinter_h__

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <algorithm>
#include <exception>

using namespace Rcpp;
using namespace std;

static void chkIntFn(void *dummy);

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt();

#endif

