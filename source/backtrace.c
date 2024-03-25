#include <execinfo.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void getCallerC (int callerSize,char *caller) {
  void *array[10];
  char **strings;
  int size, i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);
  if (size >= 4 ) {
    i = 3;
    /* Skip over unnamed functions. */
    while ( i < size && strstr(strings[i],"()") ) {
	++i;
      }
      if ( i < size ) {
	strncpy(caller,strings[i],callerSize);
      } else {
	caller = "unknown";
      }
  } else {
    caller = "unknown";
  }
  free (strings);
}
