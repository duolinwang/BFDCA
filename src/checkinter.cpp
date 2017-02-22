#include "checkinter.h"

static void chkIntFn(void *dummy) {
	R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt() {
  if(R_ToplevelExec(chkIntFn, NULL) == FALSE) stop("user interuption");

}
