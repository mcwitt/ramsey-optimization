#include "ramsey.h"

#define HISTLEN  1000
#define HASHLEN  65536

int hist[HISTLEN][NED]; /* circular array for history queue */
int (*tabu)[HASHLEN];   /* hash table of solutions in the history queue */



int main()
{
    char filename[256];
}
