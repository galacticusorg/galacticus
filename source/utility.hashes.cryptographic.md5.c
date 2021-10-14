/* Simple set of interfaces to glibc crypt function. */

#include <stdlib.h>
#include <stdio.h>
/* The crypt function is defined in crypt.h under Linux, but in unistd.h under MacOS */
#ifdef __linux__
#include <crypt.h>
#endif
#ifdef __APPLE__
#include <unistd.h>
#endif
#include <string.h>
#include <errno.h>

/* Return the md5 of the input text */
void md5(int textLength, char *text, char *hash)
{
  /* Hash the text. */
  char salt[] = "$1$........";
  /* Maximum length of input to the hash algorithm */
#ifdef CRYPT_MAX_PASSPHRASE_SIZE
  /* Use the maximum passphrase size specified for the crypt library */
  int maxLen   = CRYPT_MAX_PASSPHRASE_SIZE;
#else
  /* No maximum size specified, so just adopt something plausible */
  int maxLen   = 512;
#endif  
  int hashLen  = 22;            /* Length of the returned hash */
  int splitLen = maxLen-hashLen-1; /* Maximum size into which we can chunk the input text such that, when prepended with the previous hash, we are under the maximum size limit */
  /* Use a temporary hash which we evaluate for each chunk of text,
     and then append to the subsequent chunk. Initialize this to an
     string of whitespace. */
  char tmpHash[hashLen+1];
  int j;
  for(j=0;j<hashLen;++j) {
    tmpHash[j] = ' ';
  }
  tmpHash[hashLen] = '\0';
  /* Begin processing each chunk of the text */
  int i = 0;
  char *encrypted;
  char chunk[maxLen+1];
  while ( i < textLength ) {
    /* Determine how much to process in this chunk - either the
       maximum allowed, or the remaining amount of text */
    int j = splitLen;
    if ( i+j > textLength ) {
      j=textLength-i;
    }
    /* Join the previous chunk hash and the new chunk */
    char *concat = malloc(j+1+hashLen);
    strncpy(chunk,text+i,j);
    chunk[j] = '\0';
    strcpy(concat,tmpHash);
    strcat(concat,chunk);
    /* Encrypt it and check for errors */
    encrypted=crypt(concat, salt);
    if ( encrypted == NULL || strncmp(encrypted,"*",1) == 0 ) {
      int err = errno;
      if ( err == EINVAL ) {
        /* Salt has the wrong format */
        printf("md5(): salt has the wrong format\n");
        abort();
      } else if ( err == ERANGE ) {
        /* Text is too long */
        printf("md5(): text is too long\n");
        abort();      
      } else if ( err == ENOMEM ) {
        /* Failed to allocate memory */
        printf("md5(): failed to allocate scratch memory\n");
        abort();      
      } else if ( err == ENOSYS || err == EOPNOTSUPP ) {
        /* Crypt function is not implemented */
        printf("md5(): crypt function is not implemented\n");
        abort();      
      } else if ( err == EPERM ) {
        /* Weak encryption used */
        printf("md5(): weak encryption is not allowed\n");
        abort();      
      } else {
        printf("md5(): unknown error\n");
        abort();
      }
    } else {
      /* Extract the hash of this chunk for use in the next round */
      strncpy(tmpHash,encrypted+13,hashLen);
    }
    free(concat);
    i += splitLen;
  }
  /* Extract the full result */
  strcpy(hash,encrypted);  
  return;
}
