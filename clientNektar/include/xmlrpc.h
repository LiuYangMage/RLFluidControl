#include <stdlib.h>
#include <stdio.h>

#include <xmlrpc-c/base.h>
#include <xmlrpc-c/client.h>

//#include "config.h"  /* information about this build environment */
#include <xmlrpc-c/config.h>

#define NAME "Xmlrpc-c Test Client"
#define VERSION "1.0"

void dieIfFaultOccurred (xmlrpc_env * const envP);
