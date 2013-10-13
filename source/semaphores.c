#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include <stdlib.h>

sem_t * Semaphore_Open_C(const char *name, int initialValue) {
  sem_t *s;
  s = sem_open(name,O_CREAT,S_IRUSR | S_IWUSR,initialValue);
  if ( s == SEM_FAILED ) {
    if ( errno == EACCES       ) printf("sem_open() failed: access\n"         );
    if ( errno == EEXIST       ) printf("sem_open() failed: exists\n"         );
    if ( errno == EINVAL       ) printf("sem_open() failed: incorrect value\n");
    if ( errno == ENAMETOOLONG ) printf("sem_open() failed: name too long\n"  );
    if ( errno == ENOENT       ) printf("sem_open() failed: no such path\n"   );
    if ( errno == ENOSPC       ) printf("sem_open() failed: no space\n"       );
    abort();
  }
  return s;
}

void Semaphore_Close_C(sem_t *s) {
  int e;
  e = sem_close(s);
  if ( e == -1 ) {
    if ( errno == EINVAL ) printf("sem_close() failed: incorrect value\n");
    abort();
  }
}

void Semaphore_Wait_C(sem_t *s) {
  int e;
  e = -1;
  while ( e == -1 ) {
    e = sem_wait(s);
    if ( e == -1 ) {
      if ( errno == EINVAL ) {
	printf("sem_wait() failed: incorrect value\n");
	abort();
      }
      /* if ( errno == EINTR  ) printf("sem_wait() failed: interrupted\n"    ); */
    }
  }
}

void Semaphore_Post_C(sem_t *s) {
  int e;
  e = sem_post(s);
  if ( e == -1 ) {
    if ( errno == EINVAL ) printf("sem_post() failed: incorrect value\n");
    abort();
  }
}

void Semaphore_Unlink_C(const char *name) {
  int e;
  e = sem_unlink(name);
  if ( e == -1 ) {
    if ( errno == EACCES       ) printf("sem_unlink() failed: access\n"           );
    if ( errno == ENAMETOOLONG ) printf("sem_unlink() failed: name too long\n"    );
    if ( errno == ENOENT       ) printf("sem_unlink() failed: no such path\n"     );
    if ( errno == EFAULT       ) printf("sem_unlink() failed: incorrect address\n");
    abort();
  }
}
