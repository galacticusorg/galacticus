/* Simple set of interfaces to glibc crypt function. */

#include <stdio.h>
#include <crypt.h>

/* Return the md5 of the input text */
void md5(int textLength, char *text, char *hash)
{
  char salt[] = "$1$........";

  /* Hash the text. */
  strcpy(hash,crypt(text, salt));

  return;
}
