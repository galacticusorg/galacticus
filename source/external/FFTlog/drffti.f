c-----------------------------------------------------------------------
c http://www.netlib.org/bihar/
c-----------------------------------------------------------------------
      subroutine drffti (n,wsave)
      integer n
c wsave is a work array which should be dimensioned at least 2*n+15
      real*8 wsave(*)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine drfti1 (n,wa,rfac)
      integer n, ifac(16)
      real*8 wa(n), rfac(8)
c
      real*8 TPI
      parameter (TPI = 6.2831853071 7958647692 5286766559 00577d0)
c
      integer i, ib, ido, ii, ip, ipm, is, j, k1, l1, l2, ld,
     &  nf, nfm1, nl, nq, nr, ntry
      integer ntryh(4)
      real*8 arg, argh, argld, fi
c
      data ntryh /4, 2, 3, 5/
c
      ifac=transfer(rfac,ifac)
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
      nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      rfac=transfer(ifac,rfac)
c
      argh = TPI/dble(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = dble(ld)*argh
            fi = 0.d0
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            is = is+ido
  109    continue
c
         l1 = l2
  110 continue
c
      rfac=transfer(ifac,rfac)
      return
      end
c
