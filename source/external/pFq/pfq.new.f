C     ****************************************************************  
C     *                                                              *  
C     *    SOLUTION TO THE GENERALIZED HYPERGEOMETRIC FUNCTION       *  
C     *                                                              *  
C     *                           by                                 *  
C     *                                                              *  
C     *                      W. F. PERGER,                           *  
C     *                                                              *  
C     *              MARK NARDIN  and ATUL BHALLA                    *  
C     *                                                              *  
C     *                                                              *  
C     *            Electrical Engineering Department                 *
C     *            Michigan Technological University                 *
C     *                  1400 Townsend Drive                         *
C     *                Houghton, MI  49931-1295   USA                *
C     *                     Copyright 1993                           *  
C     *                                                              *  
C     *               e-mail address: wfp@mtu.edu                    *
C     *                                                              *  
C     *  Description : A numerical evaluator for the generalized     *  
C     *    hypergeometric function for complex arguments with large  *  
C     *    magnitudes using a direct summation of the Gauss series.  *  
C     *    The method used allows an accuracy of up to thirteen      *  
C     *    decimal places through the use of large integer arrays    *  
C     *    and a single final division.                              *  
C     *    (original subroutines for the confluent hypergeometric    *  
C     *    written by Mark Nardin, 1989; modifications made to cal-  *  
C     *    culate the generalized hypergeometric function were       *  
C     *    written by W.F. Perger and A. Bhalla, June, 1990)         *  
C     *                                                              *
C     *  The evaluation of the pFq series is accomplished by a func- *
C     *  ion call to PFQ, which is a double precision complex func-  *
C     *  tion.  The required input is:                               *
C     *  1. Double precision complex arrays A and B.  These are the  *
C     *     arrays containing the parameters in the numerator and de-*
C     *     nominator, respectively.                                 *
C     *  2. Integers IP and IQ.  These integers indicate the number  *
C     *     of numerator and denominator terms, respectively (these  *
C     *     are p and q in the pFq function).                        *
C     *  3. Double precision complex argument Z.                     *
C     *  4. Integer LNPFQ.  This integer should be set to '1' if the *
C     *     result from PFQ is to be returned as the natural logaritm*
C     *     of the series, or '0' if not.  The user can generally set*
C     *     LNPFQ = '0' and change it if required.                   *
C     *  5. Integer IX.  This integer should be set to '0' if the    *
C     *     user desires the program PFQ to estimate the number of   *
C     *     array terms (in A and B) to be used, or an integer       *
C     *     greater than zero specifying the number of integer pos-  *
C     *     itions to be used.  This input parameter is escpecially  *
C     *     useful as a means to check the results of a given run.   *
C     *     Specificially, if the user obtains a result for a given  *
C     *     set of parameters, then changes IX and re-runs the eval- *
C     *     uator, and if the number of array positions was insuffi- *
C     *     cient, then the two results will likely differ.  The rec-*
C     *     commended would be to generally set IX = '0' and then set*
C     *     it to 100 or so for a second run.  Note that the LENGTH  *
C     *     parameter currently sets the upper limit on IX to 777,   *
C     *     but that can easily be changed (it is a single PARAMETER *
C     *     statement) and the program recompiled.                   *
C     *  6. Integer NSIGFIG.  This integer specifies the requested   *
C     *     number of significant figures in the final result.  If   *
C     *     the user attempts to request more than the number of bits*
C     *     in the mantissa allows, the program will abort with an   *
C     *     appropriate error message.  The recommended value is 10. *
C     *                                                              *  
C     *     Note: The variable NOUT is the file to which error mess- *
C     *           ages are written (default is 6).  This can be      *
C     *           changed in the FUNCTION PFQ to accomodate re-      *
C     *           of output to another file                          *
C     *                                                              * 
C     *  Subprograms called: HYPER.                                  *  
C     *                                                              *  
C     ****************************************************************  
C     
      
      FUNCTION PFQ (A,B,IP,IQ,Z,LNPFQ,IX,NSIGFIG)              
      
      INTEGER LNPFQ,IX,IP,IQ,NSIGFIG                               
      COMPLEX*16 HYPER,A,B,Z,PFQ
      COMPLEX*16 CGAMMA,A1(2),B1(1),GAM1,GAM2,GAM3,GAM4,HYPER1,
     :     HYPER2,GAM5,GAM6,GAM7,Z1
      DIMENSION A(IP),B(IQ)                                               
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      COMMON/IO/NOUT
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION DNUM,PRECIS,ARGI,ARGR,DIFF
      !$omp threadprivate(/CONSTS/,/IO/)
      
      IF (Z .EQ. DCMPLX(ZERO,ZERO)) THEN
         PFQ = DCMPLX(ONE,ZERO)
         RETURN
      END IF
      IF ((LNPFQ .NE. 0) .AND. (LNPFQ .NE. 1)) THEN 
         WRITE (NOUT,*) ' ERROR IN INPUT ARGUMENTS: LNPFQ /= 0 OR 1' 
         STOP 
      END IF
      IF ((IP .GT. IQ) .AND. (ABS(Z) .GT. ONE) .AND. IP .GT. 2) THEN
         WRITE (NOUT,300) IP,IQ,ABS(Z)
         STOP
      END IF
*
* For the 2F1 function with |z| > 1, use Abramowitz and Stegun, 15.3.7.
*
      IF (IP .EQ. 2 .AND. IQ .EQ. 1 .AND. ABS(Z) .GT. ONE) THEN
         Z1 = ONE/Z
         A1(1) = A(1)
         A1(2) = ONE-B(1)+A(1)
         B1(1) = ONE-A(2)+A(1)
         GAM1 = CGAMMA(B(1),LNPFQ)
         GAM2 = CGAMMA(A(2)-A(1),LNPFQ)
         GAM3 = CGAMMA(A(2),LNPFQ)
         GAM4 = CGAMMA(B(1)-A(1),LNPFQ)
         GAM5 = CGAMMA(A(1)-A(2),LNPFQ)
         GAM6 = CGAMMA(A(1),LNPFQ)
         GAM7 = CGAMMA(B(1)-A(2),LNPFQ)
         HYPER1=HYPER(A1,B1,IP,IQ,Z1,LNPFQ,IX,NSIGFIG)
         A1(1) = A(2)
         A1(2) = ONE-B(1)+A(2)
         B1(1) = ONE-A(1)+A(2)
         HYPER2=HYPER(A1,B1,IP,IQ,Z1,LNPFQ,IX,NSIGFIG)
         PFQ = GAM1*GAM2*((-Z)**(-A(1)))*HYPER1/(GAM3*GAM4) + 
     :         GAM1*GAM5*((-Z)**(-A(2)))*HYPER2/(GAM6*GAM7)
         RETURN
      END IF
      IF (IP .EQ. 2 .AND. IQ .EQ. 1 .AND. ABS(Z) .GT. 0.9) THEN
         IF (LNPFQ .EQ. 1) GO TO 30
*     
*     Check to see if the Gamma function arguments are o.k.; if not,
*     then the series will have to be used.
*     
*     PRECIS - MACHINE PRECISION
*     
         PRECIS = ONE
 10      PRECIS = PRECIS/TWO
         DNUM = PRECIS+ONE
         IF (DNUM .GT. ONE) GOTO 10
         PRECIS = TWO*PRECIS
         DO 20 I=1,6
            IF (I .EQ. 1) THEN
               ARGI=DIMAG(B(1))
               ARGR=DBLE(B(1))
            ELSE IF (I .EQ. 2) THEN
               ARGI=DIMAG(B(1)-A(1)-A(2))
               ARGR=DBLE(B(1)-A(1)-A(2))
            ELSE IF (I .EQ. 3) THEN
               ARGI=DIMAG(B(1)-A(1))
               ARGR=DBLE(B(1)-A(1))
            ELSE IF (I .EQ. 4) THEN
               ARGI=DIMAG(A(1)+A(2)-B(1))
               ARGR=DBLE(A(1)+A(2)-B(1))
            ELSE IF (I .EQ. 5) THEN
               ARGI=DIMAG(A(1))
               ARGR=DBLE(A(1))
            ELSE IF (I .EQ. 6) THEN
               ARGI=DIMAG(A(2))
               ARGR=DBLE(A(2))
            END IF
*     
*     CASES WHERE THE ARGUMENT IS REAL
*     
            IF (ARGI .EQ. ZERO) THEN
*     
*     CASES WHERE THE ARGUMENT IS REAL AND NEGATIVE
*     
               IF (ARGR .LE. ZERO) THEN
*     
*     USE THE SERIES EXPANSION IF THE ARGUMENT IS TOO NEAR A POLE
*     
                  DIFF = ABS (DBLE (NINT (ARGR))-ARGR)
                  IF (DIFF .LE. TWO*PRECIS) GO TO 30
               END IF
            END IF
 20      CONTINUE
         GAM1=CGAMMA(B(1),LNPFQ)
         GAM2=CGAMMA(B(1)-A(1)-A(2),LNPFQ)
         GAM3=CGAMMA(B(1)-A(1),LNPFQ)
         GAM4=CGAMMA(B(1)-A(2),LNPFQ)
         GAM5=CGAMMA(A(1)+A(2)-B(1),LNPFQ)
         GAM6=CGAMMA(A(1),LNPFQ)
         GAM7=CGAMMA(A(2),LNPFQ)
         A1(1)=A(1)
         A1(2)=A(2)
         B1(1)=A(1)+A(2)-B(1)+ONE
         Z1=ONE-Z
         HYPER1=HYPER(A1,B1,IP,IQ,Z1,LNPFQ,IX,NSIGFIG)
         A1(1)=B(1)-A(1)
         A1(2)=B(1)-A(2)
         B1(1)=B(1)-A(1)-A(2)+ONE
         HYPER2=HYPER(A1,B1,IP,IQ,Z1,LNPFQ,IX,NSIGFIG)
         PFQ=GAM1*GAM2*HYPER1/(GAM3*GAM4)+(ONE-Z)**(B(1)-A(1)-A(2))*
     :        GAM1*GAM5*HYPER2/(GAM6*GAM7)
         RETURN
      END IF
 30   CONTINUE 
      PFQ=HYPER(A,B,IP,IQ,Z,LNPFQ,IX,NSIGFIG)
      RETURN 
 300  FORMAT (/,1X,'IP=',1I2,3X,'IQ=',1I2,3X,'AND ABS(Z)=',
     :     1E12.5,2X,/,' WHICH IS GREATER THAN ONE--SERIES DOES',
     :     ' NOT CONVERGE')
      END                                                               
C     ****************************************************************  
C     *                                                              *  
C     *                   FUNCTION BITS                              *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Determines the number of significant figures  *  
C     *    of machine precision to arrive at the size of the array   *  
C     *    the numbers must be stored in to get the accuracy of the  *
C     *    solution.                                                 *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      DOUBLE PRECISION FUNCTION BITS() 
      
      DOUBLE PRECISION BIT,BIT2                                                   
      INTEGER COUNT                                                     
      
      BIT=1.0                                                           
      COUNT=0                                                           
 10   COUNT=COUNT+1                                                     
      BIT2=BIT*2.0                                                      
      BIT=BIT2+1.0                                                      
      IF ((BIT-BIT2) .NE. 0.0) GOTO 10                                  
      BITS=COUNT-3                                                      
      RETURN                                                            
      END                                                               
C     ****************************************************************  
C     *                                                              *  
C     *                   FUNCTION HYPER                             *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Function that sums the Gauss series.          *  
C     *                                                              *  
C     *  Subprograms called: ARMULT, ARYDIV, BITS, CMPADD, CMPMUL,   *  
C     *                      IPREMAX.                                *
C     *                                                              *  
C     ****************************************************************  
      
      FUNCTION HYPER (A,B,IP,IQ,Z,LNPFQ,IX,NSIGFIG) 
      
      PARAMETER (LENGTH=25000)
      INTEGER L,I,IBIT,LNPFQ,IP,IQ,IX,NSIGFIG,NMACH,ICOUNT,IXCNT
      INTEGER REXP,IR10,II10,LMAX,IPREMAX
      COMPLEX*16 A,B,Z,FINAL,TEMP,OLDTEMP,HYPER,FACTOR,CDUM1,CDUM2
      COMPLEX*16 TEMP1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SUMR, SUMI, 
     :     NUMR, NUMI, DENOMR,DENOMI, QR1, QR2, QI1, QI2, WK1,
     :     WK2, WK3, WK4, WK5, WK6
      DIMENSION A(IP),B(IQ),AR(10),AI(10),AR2(10),AI2(10),                
     :     CR(10),CI(10),CR2(10),CI2(10)                           
      DOUBLE PRECISION AR,AI,CR,CI,XR,XI,CNT,SIGFIG
      DOUBLE PRECISION RMAX,MX1,MX2,CREAL
      DOUBLE PRECISION AR2,AI2,CR2,CI2,XR2,XI2 
      DOUBLE PRECISION ACCY,BITS,EXPON,XL,X
      DOUBLE PRECISION DUM1,DUM2,RR10,RI10,LOG2
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      COMMON/IO/NOUT
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/,/IO/)

      allocate(SUMR(-1:LENGTH))
      allocate(SUMI(-1:LENGTH))
      allocate(NUMR(-1:LENGTH))
      allocate(NUMI(-1:LENGTH))
      allocate(DENOMR(-1:LENGTH))
      allocate(DENOMI(-1:LENGTH)) 
      allocate(QR1(-1:LENGTH))
      allocate(QR2(-1:LENGTH))
      allocate(QI1(-1:LENGTH))
      allocate(QI2(-1:LENGTH))         
      allocate(WK1(-1:LENGTH))
      allocate(WK2(-1:LENGTH))
      allocate(WK3(-1:LENGTH))
      allocate(WK4(-1:LENGTH))
      allocate(WK5(-1:LENGTH))
      allocate(WK6(-1:LENGTH))

      CREAL=0.0d0
      LOG2=LOG10(TWO)
      IBIT=INT(BITS())                                                       
      RMAX=TWO**(IBIT/2)                                               
      SIGFIG=TWO**(IBIT/4)
      TEMP1=DCMPLX(0.0D0,0.0D0)
*     
      DO 10 I1=1,IP                                                      
         AR2(I1)=DBLE(A(I1))*SIGFIG                                   
         AR(I1)=AINT(AR2(I1))                                          
         AR2(I1)=ANINT((AR2(I1)-AR(I1))*RMAX)                          
         AI2(I1)=DIMAG(A(I1))*SIGFIG                                   
         AI(I1)=AINT(AI2(I1))                                          
         AI2(I1)=ANINT((AI2(I1)-AI(I1))*RMAX)                          
 10   CONTINUE                                                          
      DO 20 I1=1,IQ                                                      
         CR2(I1)=DBLE(B(I1))*SIGFIG                                   
         CR(I1)=AINT(CR2(I1))                                          
         CR2(I1)=ANINT((CR2(I1)-CR(I1))*RMAX)                          
         CI2(I1)=DIMAG(B(I1))*SIGFIG                                   
         CI(I1)=AINT(CI2(I1))                                          
         CI2(I1)=ANINT((CI2(I1)-CI(I1))*RMAX)                          
 20   CONTINUE                                                          
      XR2=DBLE(Z)*SIGFIG                                               
      XR=AINT(XR2)                                                      
      XR2=ANINT((XR2-XR)*RMAX)                                          
      XI2=DIMAG(Z)*SIGFIG                                               
      XI=AINT(XI2)                                                      
      XI2=ANINT((XI2-XI)*RMAX)                                          
*     
*     WARN THE USER THAT THE INPUT VALUE WAS SO CLOSE TO ZERO THAT IT
*     WAS SET EQUAL TO ZERO.
*     
      DO 30 I1=1,IP
         IF ((DBLE(A(I1)) .NE. ZERO) .AND. (AR(I1) .EQ. ZERO)
     :        .AND. (AR2(I1) .EQ. ZERO))
     :        WRITE (NOUT,300) I1
         IF ((DIMAG(A(I1)) .NE. ZERO) .AND. (AI(I1) .EQ. ZERO)
     :        .AND. (AI2(I1) .EQ. ZERO))
     :        WRITE (NOUT,301) I1
 30   CONTINUE
      DO 40 I1=1,IQ
         IF ((DBLE(B(I1)) .NE. ZERO) .AND. (CR(I1) .EQ. ZERO)
     :        .AND. (CR2(I1) .EQ. ZERO))
     :        WRITE (NOUT,302) I1
         IF ((DIMAG(B(I1)) .NE. ZERO) .AND. (CI(I1) .EQ. ZERO)
     :        .AND. (CI2(I1) .EQ. ZERO))
     :        WRITE (NOUT,303) I1
 40   CONTINUE
      IF ((DBLE(Z).NE.ZERO) .AND. (XR.EQ.ZERO) .AND. (XR2.EQ.ZERO))
     :     THEN
         WRITE (NOUT,*) ' WARNING - REAL PART OF Z WAS SET TO ZERO'
         Z = DCMPLX(ZERO,DIMAG(Z))
      END IF
      IF ((DIMAG(Z).NE.ZERO) .AND. (XI.EQ.ZERO) .AND. (XI2.EQ.ZERO))
     :     THEN
         WRITE (NOUT,*) ' WARNING - IMAG PART OF Z WAS SET TO ZERO'
         Z = DCMPLX(DBLE(Z),ZERO)
      END IF
      
*     
*     SCREENING OF NUMERATOR ARGUMENTS FOR NEGATIVE INTEGERS OR ZERO.
*     ICOUNT WILL FORCE THE SERIES TO TERMINATE CORRECTLY.
*     
      NMACH=INT(LOG10(TWO**INT(BITS())))
      ICOUNT=-1
      DO 50 I1=1,IP
         IF ((AR2(I1) .EQ. ZERO) .AND. (AR(I1) .EQ. ZERO) .AND.
     :        (AI2(I1) .EQ. ZERO) .AND. (AI(I1) .EQ. ZERO)) THEN
            HYPER=DCMPLX(ONE,ZERO)
            RETURN
         END IF
         IF ((AI(I1) .EQ. ZERO) .AND. (AI2(I1) .EQ. ZERO) .AND.
     :        (DBLE(A(I1)) .LT. ZERO)) THEN
            IF (ABS(DBLE(A(I1))-DBLE(NINT(DBLE(A(I1)))))
     :           .LT. TEN**(-NMACH)) THEN
               IF (ICOUNT .NE. -1) THEN
                  ICOUNT=MIN(ICOUNT,-NINT(DBLE(A(I1))))
               ELSE
                  ICOUNT=-NINT(DBLE(A(I1)))
               END IF
            END IF
         END IF
 50   CONTINUE
*     
*     SCREENING OF DENOMINATOR ARGUMENTS FOR ZEROES OR NEGATIVE INTEGERS.
*     
      DO 60 I1=1,IQ
         IF ((CR(I1) .EQ. ZERO) .AND. (CR2(I1) .EQ. ZERO) .AND.
     :        (CI(I1) .EQ. ZERO) .AND. (CI2(I1) .EQ. ZERO)) THEN
            WRITE (NOUT,304) I1
            STOP
         END IF
         IF ((CI(I1) .EQ. ZERO) .AND. (CI2(I1) .EQ. ZERO) .AND.
     :        (DBLE(B(I1)) .LT. ZERO)) THEN
            IF ((ABS(DBLE(B(I1))-DBLE(NINT(DBLE(B(I1)))))
     :           .LT. TEN**(-NMACH)) .AND. 
     :           (ICOUNT .GE. -NINT(DBLE(B(I1))) .OR.
     :           ICOUNT .EQ. -1)) THEN
               WRITE (NOUT,305) I1
               STOP
            END IF
         END IF
 60   CONTINUE
*     
      NMACH=INT(LOG10(TWO**IBIT))
      NSIGFIG=MIN(NSIGFIG,INT(LOG10(TWO**IBIT)))
      ACCY = TEN**(-NSIGFIG)                                         
      L=IPREMAX(A,B,IP,IQ,Z)
      IF (L .EQ. 1) GO TO 110
*     
*     First, estimate the exponent of the maximum term in the pFq series.
*     
      EXPON=ZERO
      XL=DBLE(L)
      DO 70 I=1,IP
         EXPON=EXPON+DBLE(FACTOR(A(I)+XL-ONE))-DBLE(FACTOR(A(I)-ONE))
 70   CONTINUE
      DO 80 I=1,IQ
         EXPON=EXPON-DBLE(FACTOR(B(I)+XL-ONE))+DBLE(FACTOR(B(I)-ONE))
 80   CONTINUE
      EXPON = EXPON + XL*DBLE(LOG(Z)) - DBLE(FACTOR(DCMPLX(XL,ZERO)))
      LMAX=INT(LOG10(EXP(ONE))*EXPON)
      L=LMAX
*     
*     Now, estimate the exponent of where the pFq series will terminate.
*     
      TEMP1=DCMPLX(ONE,ZERO)
      CREAL=ONE
      DO 90 I1=1,IP
         TEMP1=TEMP1*DCMPLX(AR(I1),AI(I1))/SIGFIG
 90   CONTINUE
      DO 100 I1=1,IQ
         TEMP1=TEMP1/(DCMPLX(CR(I1),CI(I1))/SIGFIG)
         CREAL=CREAL*CR(I1)
 100  CONTINUE
      TEMP1=TEMP1*DCMPLX(XR,XI)
*     
*     Triple it to make sure.
*     
      L=3*L
*     
*     Divide the number of significant figures necessary by the number of
*     digits available per array position.
*     
      L=INT((2*L+NSIGFIG)/NMACH)+2
*     
*     Make sure there are at least 5 array positions used.
*     
 110  L=MAX(L,5)
      L=MAX(L,IX)
*     write (6,*) ' Estimated value of L=',L
      IF ((L .LT. 0) .OR. (L .GT. LENGTH)) THEN                          
         WRITE (NOUT,306) LENGTH
         STOP                                                          
      END IF                                                            
      IF (NSIGFIG .GT. NMACH) THEN
         WRITE (NOUT,307) NMACH
      END IF
      
      SUMR(-1)=ONE                                                    
      SUMI(-1)=ONE                                                    
      NUMR(-1)=ONE                                                    
      NUMI(-1)=ONE                                                    
      DENOMR(-1)=ONE                                                  
      DENOMI(-1)=ONE                                                  
      DO 120 I=0,L+1                                                    
         SUMR(I)=ZERO                                                   
         SUMI(I)=ZERO                                                   
         NUMR(I)=ZERO                                                   
         NUMI(I)=ZERO                                                   
         DENOMR(I)=ZERO                                                 
         DENOMI(I)=ZERO                                                 
 120  CONTINUE                                                          
      SUMR(1)=ONE                                                     
      NUMR(1)=ONE                                                     
      DENOMR(1)=ONE                                                   
      CNT=SIGFIG                                                        
      TEMP=DCMPLX(ZERO,ZERO)                                          
      OLDTEMP=TEMP                                                      
      IXCNT=0
      REXP=IBIT/2
      X=REXP*(SUMR(L+1)-2)                                                
      RR10=X*LOG2
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(SUMI(L+1)-2)                                                
      RI10=X*LOG2
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=SIGN(SUMR(1)*RMAX*RMAX+SUMR(2)*RMAX+SUMR(3),SUMR(-1))               
      DUM2=SIGN(SUMI(1)*RMAX*RMAX+SUMI(2)*RMAX+SUMI(3),SUMI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CDUM1=DCMPLX(DUM1,DUM2)
      X=REXP*(DENOMR(L+1)-2)                                                
      RR10=X*LOG2
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(DENOMI(L+1)-2)                                                
      RI10=X*LOG2
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=SIGN(DENOMR(1)*RMAX*RMAX+DENOMR(2)*RMAX+
     :     DENOMR(3),DENOMR(-1))               
      DUM2=SIGN(DENOMI(1)*RMAX*RMAX+DENOMI(2)*RMAX+
     :     DENOMI(3),DENOMI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CDUM2=DCMPLX(DUM1,DUM2)
      TEMP=CDUM1/CDUM2
      
*     130 IF (IP .GT. 0) THEN
 130  if (ip .lt. 0) then
         IF (SUMR(1) .LT. HALF) THEN
            MX1=SUMI(L+1)
         ELSE IF (SUMI(1) .LT. HALF) THEN                                   
            MX1=SUMR(L+1)                                                   
         ELSE                                                              
            MX1=DMAX1(SUMR(L+1),SUMI(L+1))                                  
         ENDIF                                                             
         IF (NUMR(1) .LT. HALF) THEN                                        
            MX2=NUMI(L+1)                                                   
         ELSE IF (NUMI(1) .LT. HALF) THEN                                   
            MX2=NUMR(L+1)                                                   
         ELSE                                                              
            MX2=DMAX1(NUMR(L+1),NUMI(L+1))                                  
         ENDIF                                                             
         IF (MX1-MX2 .GT.  2.0) THEN                                       
            IF (CREAL .GE. ZERO) THEN                                         
*     write (6,*) ' cdabs(temp1/cnt)=',cdabs(temp1/cnt)
               IF (CDABS(TEMP1/CNT) .LE. ONE) GOTO 240
            ENDIF                                                           
         ENDIF
      ELSE
         CALL ARYDIV(SUMR,SUMI,DENOMR,DENOMI,TEMP,L,LNPFQ,RMAX,IBIT)        
*     
*     First, estimate the exponent of the maximum term in the pFq series.
*     
         EXPON=ZERO
         XL=DBLE(ixcnt)
         DO 140 I=1,IP
            EXPON=EXPON+DBLE(FACTOR(A(I)+XL-ONE))-DBLE(FACTOR(A(I)-ONE))
 140     CONTINUE
         DO 150 I=1,IQ
            EXPON=EXPON-DBLE(FACTOR(B(I)+XL-ONE))+DBLE(FACTOR(B(I)-ONE))
 150     CONTINUE
         EXPON = EXPON + XL*DBLE(LOG(Z)) - DBLE(FACTOR(DCMPLX(XL,ZERO)))
         LMAX=INT(LOG10(EXP(ONE))*EXPON)
         IF (ABS(OLDTEMP-TEMP) .LT. ABS(TEMP*ACCY)) GO TO 240                
         OLDTEMP=TEMP                                                      
      END IF
      IF (IXCNT .EQ. ICOUNT) GO TO 240
      IXCNT=IXCNT+1
      DO 160 I1=1,IQ                                                      
*     
*     TAKE THE CURRENT SUM AND MULTIPLY BY THE DENOMINATOR OF THE NEXT
*     TERM, FOR BOTH THE MOST SIGNIFICANT HALF (CR,CI) AND THE LEAST
*     SIGNIFICANT HALF (CR2,CI2).
*     
         CALL CMPMUL(SUMR,SUMI,CR(I1),CI(I1),QR1,QI1,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)           
         CALL CMPMUL(SUMR,SUMI,CR2(I1),CI2(I1),QR2,QI2,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)         
         QR2(L+1)=QR2(L+1)-1                                           
         QI2(L+1)=QI2(L+1)-1                                           
*     
*     STORE THIS TEMPORARILY IN THE SUM ARRAYS.
*     
         CALL CMPADD(QR1,QI1,QR2,QI2,SUMR,SUMI,WK1,L,RMAX)                 
 160  CONTINUE                                                          
      
*     
*     MULTIPLY BY THE FACTORIAL TERM.
*     
      CALL ARMULT(SUMR,CNT,SUMR,WK6,L,RMAX)                                 
      CALL ARMULT(SUMI,CNT,SUMI,WK6,L,RMAX)                                 
*     
*     MULTIPLY BY THE SCALING FACTOR, SIGFIG, TO KEEP THE SCALE CORRECT.
*     
      DO 170 I1=1,IP-IQ                                                   
         CALL ARMULT(SUMR,SIGFIG,SUMR,WK6,L,RMAX)                          
         CALL ARMULT(SUMI,SIGFIG,SUMI,WK6,L,RMAX)                          
 170  CONTINUE                                                          
      DO 180 I1=1,IQ                                                     
*     
*     UPDATE THE DENOMINATOR.
*     
         CALL CMPMUL(DENOMR,DENOMI,CR(I1),CI(I1),QR1,QI1,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)       
         CALL CMPMUL(DENOMR,DENOMI,CR2(I1),CI2(I1),QR2,QI2,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)     
         QR2(L+1)=QR2(L+1)-1                                           
         QI2(L+1)=QI2(L+1)-1                                           
         CALL CMPADD(QR1,QI1,QR2,QI2,DENOMR,DENOMI,WK1,L,RMAX)             
 180  CONTINUE                                                          
      
*     
*     MULTIPLY BY THE FACTORIAL TERM.
*     
      CALL ARMULT(DENOMR,CNT,DENOMR,WK6,L,RMAX)                             
      CALL ARMULT(DENOMI,CNT,DENOMI,WK6,L,RMAX)                             
*     
*     MULTIPLY BY THE SCALING FACTOR, SIGFIG, TO KEEP THE SCALE CORRECT.
*     
      DO 190 I1=1,IP-IQ                                                  
         CALL ARMULT(DENOMR,SIGFIG,DENOMR,WK6,L,RMAX)                      
         CALL ARMULT(DENOMI,SIGFIG,DENOMI,WK6,L,RMAX)                      
 190  CONTINUE                                                          
*     
*     FORM THE NEXT NUMERATOR TERM BY MULTIPLYING THE CURRENT 
*     NUMERATOR TERM (AN ARRAY) WITH THE A ARGUMENT (A SCALAR).
*     
      DO 200 I1=1,IP                                                      
         CALL CMPMUL(NUMR,NUMI,AR(I1),AI(I1),QR1,QI1,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)           
         CALL CMPMUL(NUMR,NUMI,AR2(I1),AI2(I1),QR2,QI2,
     :        WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)         
         QR2(L+1)=QR2(L+1)-1                                           
         QI2(L+1)=QI2(L+1)-1                                           
         CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,WK1,L,RMAX)                 
 200  CONTINUE                                                          
*     
*     FINISH THE NEW NUMERATOR TERM BY MULTIPLYING BY THE Z ARGUMENT.
*     
      CALL CMPMUL(NUMR,NUMI,XR,XI,QR1,QI1,
     :     WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)                       
      CALL CMPMUL(NUMR,NUMI,XR2,XI2,QR2,QI2,
     :     WK1,WK2,WK3,WK4,WK5,WK6,L,RMAX)                     
      QR2(L+1)=QR2(L+1)-1                                               
      QI2(L+1)=QI2(L+1)-1                                               
      CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,WK1,L,RMAX)                     
*     
*     MULTIPLY BY THE SCALING FACTOR, SIGFIG, TO KEEP THE SCALE CORRECT.
*     
      DO 210 I1=1,IQ-IP                                                   
         CALL ARMULT(NUMR,SIGFIG,NUMR,WK6,L,RMAX)                          
         CALL ARMULT(NUMI,SIGFIG,NUMI,WK6,L,RMAX)                          
 210  CONTINUE                                                          
*     
*     FINALLY, ADD THE NEW NUMERATOR TERM WITH THE CURRENT RUNNING
*     SUM OF THE NUMERATOR AND STORE THE NEW RUNNING SUM IN SUMR, SUMI.
*     
      CALL CMPADD(SUMR,SUMI,NUMR,NUMI,SUMR,SUMI,WK1,L,RMAX)                 
*     
*     BECAUSE SIGFIG REPRESENTS "ONE" ON THE NEW SCALE, ADD SIGFIG
*     TO THE CURRENT COUNT AND, CONSEQUENTLY, TO THE IP ARGUMENTS
*     IN THE NUMERATOR AND THE IQ ARGUMENTS IN THE DENOMINATOR.
*     
      CNT=CNT+SIGFIG                                                    
      DO 220 I1=1,IP                                                     
         AR(I1)=AR(I1)+SIGFIG                                          
 220  CONTINUE                                                          
      DO 230 I1=1,IQ                                                     
         CR(I1)=CR(I1)+SIGFIG                                          
 230  CONTINUE                                                          
      GOTO 130                                                          
 240  CALL ARYDIV(SUMR,SUMI,DENOMR,DENOMI,FINAL,L,LNPFQ,RMAX,IBIT)       
*     write (6,*) 'Number of terms=',ixcnt
      HYPER=FINAL                                                       
      RETURN                                                            
 300  FORMAT (1X,'WARNING - REAL PART OF A(',1I2,') WAS SET TO ZERO')
 301  FORMAT (1X,'WARNING - IMAG PART OF A(',1I2,') WAS SET TO ZERO')
 302  FORMAT (1X,'WARNING - REAL PART OF B(',1I2,') WAS SET TO ZERO')
 303  FORMAT (1X,'WARNING - IMAG PART OF B(',1I2,') WAS SET TO ZERO')
 304  FORMAT (1X,'ERROR - ARGUMENT B(',1I2,') WAS EQUAL TO ZERO') 
 305  FORMAT (1X,'ERROR - ARGUMENT B(',1I2,') WAS A NEGATIVE',
     :     ' INTEGER')
 306  FORMAT (1X,'ERROR IN FN HYPER: L MUST BE < ',1I4)
 307  FORMAT (1X,' WARNING--THE NUMBER OF SIGNIFICANT FIGURES REQU',
     :     'ESTED',/,'IS GREATER THAN THE MACHINE PRECISION--',
     :     'FINAL ANSWER',/,'WILL BE ACCURATE TO ONLY',I3,' DIGITS')
      END                                                               
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ARADD                             *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Accepts two arrays of numbers and returns     *  
C     *    the sum of the array.  Each array is holding the value    *  
C     *    of one number in the series.  The parameter L is the      *  
C     *    size of the array representing the number and RMAX is     *  
C     *    the actual number of digits needed to give the numbers    *  
C     *    the desired accuracy.                                     *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ARADD(A,B,C,Z,L,RMAX) 
      
      INTEGER L                                                         
      DOUBLE PRECISION A,B,C,Z,RMAX
      INTEGER EDIFF,I,J                                                 
      DIMENSION A(-1:*),B(-1:*),C(-1:*),Z(-1:*)                 
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS

      DO 10 I=0,L+1                                                    
         Z(I)=ZERO                                                      
 10   CONTINUE                                                          
      EDIFF=INT(DNINT(A(L+1)-B(L+1)))
      IF (ABS(A(1)) .LT. HALF .OR. EDIFF .LE. -L) GOTO 20
      IF (ABS(B(1)) .LT. HALF .OR. EDIFF .GE. L) GOTO 40               
      GOTO 60                                                         
 20   DO 30 I=-1,L+1                                                   
         C(I)=B(I)                                                       
 30   CONTINUE                                                          
      GOTO 350                                                          
 40   DO 50 I=-1,L+1                                                   
         C(I)=A(I)                                                       
 50   CONTINUE                                                          
      GOTO 350                                                          
 60   Z(-1)=A(-1)                                                       
      IF (ABS(A(-1)-B(-1)) .LT. HALF) GOTO 80
      IF (EDIFF .GT. 0) THEN                                            
         Z(L+1)=A(L+1)                                                   
         GOTO 190                                                     
      ENDIF                                                             
      IF (EDIFF .LT. 0) THEN                                            
         Z(L+1)=B(L+1)                                                   
         Z(-1)=B(-1)                                                     
         GOTO 240                                                     
      ENDIF                                                             
      DO 70 I=1,L                                                      
         IF (A(I) .GT. B(I)) THEN                                        
            Z(L+1)=A(L+1)                                                 
            GOTO 190                                                      
         ENDIF                                                           
         IF (A(I) .LT. B(I)) THEN                                        
            Z(L+1)=B(L+1)                                                 
            Z(-1)=B(-1)                                                   
            GOTO 240                                                      
         ENDIF                                                           
 70   CONTINUE                                                          
      GOTO 300                                                          
      
 80   IF (EDIFF .GT. 0) GOTO 110                                        
      IF (EDIFF .LT. 0) GOTO 150
      Z(L+1)=A(L+1)                                                     
      DO 90 I=L,1,-1                                                   
         Z(I)=A(I)+B(I)+Z(I)                                             
         IF (Z(I) .GE. RMAX) THEN                                        
            Z(I)=Z(I)-RMAX                                                
            Z(I-1)=ONE                                                  
         ENDIF                                                           
 90   CONTINUE                                                          
      IF (Z(0) .GT. HALF) THEN                                           
         DO 100 I=L,1,-1                                                 
            Z(I)=Z(I-1)                                                   
 100     CONTINUE                                                        
         Z(L+1)=Z(L+1)+ONE                                             
         Z(0)=ZERO                                                      
      ENDIF                                                             
      GOTO 300                                                          
 110  Z(L+1)=A(L+1)                                                     
      DO 120 I=L,1+EDIFF,-1                                             
         Z(I)=A(I)+B(I-EDIFF)+Z(I)                                       
         IF (Z(I) .GE. RMAX) THEN                                        
            Z(I)=Z(I)-RMAX                                                
            Z(I-1)=ONE                                                  
         ENDIF                                                           
 120  CONTINUE                                                          
      DO 130 I=EDIFF,1,-1                                               
         Z(I)=A(I)+Z(I)                                                  
         IF (Z(I) .GE. RMAX) THEN                                        
            Z(I)=Z(I)-RMAX                                                
            Z(I-1)=ONE                                                  
         ENDIF                                                           
 130  CONTINUE                                                          
      IF (Z(0) .GT. HALF) THEN                                           
         DO 140 I=L,1,-1                                                 
            Z(I)=Z(I-1)                                                   
 140     CONTINUE                                                        
         Z(L+1)=Z(L+1)+1                                                 
         Z(0)=ZERO                                                      
      ENDIF                                                             
      GOTO 300                                                          
 150  Z(L+1)=B(L+1)                                                     
      DO 160 I=L,1-EDIFF,-1                                             
         Z(I)=A(I+EDIFF)+B(I)+Z(I)                                       
         IF (Z(I) .GE. RMAX) THEN                                        
            Z(I)=Z(I)-RMAX                                                
            Z(I-1)=ONE                                                  
         ENDIF                                                           
 160  CONTINUE                                                          
      DO 170 I=0-EDIFF,1,-1                                             
         Z(I)=B(I)+Z(I)                                                  
         IF (Z(I) .GE. RMAX) THEN                                        
            Z(I)=Z(I)-RMAX                                                
            Z(I-1)=ONE                                                  
         ENDIF                                                           
 170  CONTINUE                                                          
      IF (Z(0) .GT. HALF) THEN                                           
         DO 180 I=L,1,-1                                                 
            Z(I)=Z(I-1)                                                   
 180     CONTINUE                                                        
         Z(L+1)=Z(L+1)+ONE                                             
         Z(0)=ZERO                                                      
      ENDIF                                                             
      GOTO 300                                                          
      
 190  IF (EDIFF .GT. 0) GOTO 210                                        
      DO 200 I=L,1,-1                                                   
         Z(I)=A(I)-B(I)+Z(I)                                             
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 200  CONTINUE                                                          
      GOTO 290                                                          
 210  DO 220 I=L,1+EDIFF,-1                                             
         Z(I)=A(I)-B(I-EDIFF)+Z(I)                                       
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 220  CONTINUE                                                          
      DO 230 I=EDIFF,1,-1                                               
         Z(I)=A(I)+Z(I)                                                  
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 230  CONTINUE                                                          
      GOTO 290                                                          
      
 240  IF (EDIFF .LT. 0) GOTO 260                                        
      DO 250 I=L,1,-1                                                   
         Z(I)=B(I)-A(I)+Z(I)                                             
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 250  CONTINUE                                                          
      GOTO 290                                                          
 260  DO 270 I=L,1-EDIFF,-1                                             
         Z(I)=B(I)-A(I+EDIFF)+Z(I)                                       
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 270  CONTINUE                                                          
      DO 280 I=0-EDIFF,1,-1                                             
         Z(I)=B(I)+Z(I)                                                  
         IF (Z(I) .LT. ZERO) THEN                                       
            Z(I)=Z(I)+RMAX                                                
            Z(I-1)=-ONE                                                 
         ENDIF                                                           
 280  CONTINUE                                                          
      
 290  IF (Z(1) .GT. HALF) GOTO 330                      
      I=1                                                               
 300  I=I+1                                                           
      IF (Z(I) .LT. HALF .AND. I .LT. L+1) GOTO 300
      IF (I .EQ. L+1) THEN                                              
         Z(-1)=ONE                                                     
         Z(L+1)=ZERO                                                    
         GOTO 330                                      
      ENDIF                                                             
      DO 310 J=1,L+1-I                                                  
         Z(J)=Z(J+I-1)                                                   
 310  CONTINUE                                                          
      DO 320 J=L+2-I,L                                                  
         Z(J)=ZERO                                                      
 320  CONTINUE                                                          
      Z(L+1)=Z(L+1)-I+1                                                 
 330  DO 340 I=-1,L+1                                                   
         C(I)=Z(I)                                                       
 340  CONTINUE                                                          
 350  IF (C(1) .LT. HALF) THEN                                           
         C(-1)=ONE                                                     
         C(L+1)=ZERO                                                    
      ENDIF                                                             
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ARSUB                             *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Accepts two arrays and subtracts each element *  
C     *    in the second array from the element in the first array   *  
C     *    and returns the solution.  The parameters L and RMAX are  *  
C     *    the size of the array and the number of digits needed for *  
C     *    the accuracy, respectively.                               *  
C     *                                                              *  
C     *  Subprograms called: ARADD                                   *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ARSUB(A,B,C,WK1,WK2,L,RMAX)                                    
      
      INTEGER L,I                                                       
      DOUBLE PRECISION A,B,C,WK1,WK2,RMAX
      DIMENSION A(-1:*),B(-1:*),C(-1:*),WK1(-1:*),WK2(-1:*)
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      
      DO 10 I=-1,L+1                                                   
         WK2(I)=B(I)                                                      
 10   CONTINUE                                                          
      WK2(-1)=(-ONE)*WK2(-1)                                            
      CALL ARADD(A,WK2,C,WK1,L,RMAX)                                         
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ARMULT                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Accepts two arrays and returns the product.   *  
C     *    L and RMAX are the size of the arrays and the number of   *  
C     *    digits needed to represent the numbers with the required  *  
C     *    accuracy.                                                 *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ARMULT(A,B,C,Z,L,RMAX)                                   
      
      INTEGER L                                                         
      DOUBLE PRECISION A,B,C,Z,B2,CARRY,RMAX
      DIMENSION A(-1:*),C(-1:*),Z(-1:*)                           
      INTEGER I                                                         
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)

      Z(-1)=SIGN(ONE,B)*A(-1)                                        
      B2=ABS(B)                                                        
      Z(L+1)=A(L+1)                                                     
      DO 10 I=0,L                                                      
         Z(I)=ZERO                                                      
 10   CONTINUE                                                          
      IF (B2 .LE. EPS .OR. A(1) .LE. EPS) THEN                  
         Z(-1)=ONE                                                     
         Z(L+1)=ZERO                                                    
         GOTO 60                                                        
      ENDIF                                                             
      DO 20 I=L,1,-1                                                   
         Z(I)=A(I)*B2+Z(I)                                               
         IF (Z(I) .GE. RMAX) THEN                                        
            CARRY=AINT(Z(I)/RMAX)                                         
            Z(I)=Z(I)-CARRY*RMAX                                          
            Z(I-1)=CARRY                                                  
         ENDIF                                                           
 20   CONTINUE                                                          
      IF (Z(0) .LT. HALF) GOTO 50                                       
      DO 30 I=L,1,-1                                                   
         Z(I)=Z(I-1)                                                     
 30   CONTINUE                                                          
      Z(L+1)=Z(L+1)+ONE                                               
      IF (Z(1) .GE. RMAX) THEN
         DO 40 I=L,1,-1
            Z(I)=Z(I-1)
 40      CONTINUE
         CARRY=AINT(Z(1)/RMAX)
         Z(2)=Z(2)-CARRY*RMAX
         Z(1)=CARRY
         Z(L+1)=Z(L+1)+ONE
      END IF
      Z(0)=ZERO                                                        
 50   CONTINUE                                                          
 60   DO 70 I=-1,L+1                                                   
         C(I)=Z(I)                                                       
 70   CONTINUE                                                          
      IF (C(1) .LT. HALF) THEN                                           
         C(-1)=ONE                                                     
         C(L+1)=ZERO                                                    
      ENDIF                                                             
      RETURN                                                            
      END                                                               
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE CMPADD                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Takes two arrays representing one real and    *  
C     *    one imaginary part, and adds two arrays representing      *  
C     *    another complex number and returns two array holding the  *  
C     *    complex sum.                                              *  
C     *              (CR,CI) = (AR+BR, AI+BI)                        *  
C     *                                                              *  
C     *  Subprograms called: ARADD                                   *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE CMPADD(AR,AI,BR,BI,CR,CI,WK1,L,RMAX)                       
      
      INTEGER L                                                         
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI,RMAX,WK1
      DIMENSION AR(-1:*),AI(-1:*),BR(-1:*),BI(-1:*)             
      DIMENSION CR(-1:*),CI(-1:*),WK1(-1:*)
      
      CALL ARADD(AR,BR,CR,WK1,L,RMAX)                                       
      CALL ARADD(AI,BI,CI,WK1,L,RMAX)                                       
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE CMPSUB                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Takes two arrays representing one real and    *  
C     *    one imaginary part, and subtracts two arrays representing *  
C     *    another complex number and returns two array holding the  *  
C     *    complex sum.                                              *  
C     *              (CR,CI) = (AR+BR, AI+BI)                        *  
C     *                                                              *  
C     *  Subprograms called: ARADD                                   *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE CMPSUB(AR,AI,BR,BI,CR,CI,WK1,WK2,L,RMAX)                       
      
      INTEGER L                                                         
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI,RMAX,WK1,WK2
      DIMENSION AR(-1:*),AI(-1:*),BR(-1:*),BI(-1:*)             
      DIMENSION CR(-1:*),CI(-1:*),WK1(-1:*),WK2(-1:*)               
      
      CALL ARSUB(AR,BR,CR,WK1,WK2,L,RMAX)                                       
      CALL ARSUB(AI,BI,CI,WK1,WK2,L,RMAX)                                       
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE CMPMUL                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Takes two arrays representing one real and    *  
C     *    one imaginary part, and multiplies it with two arrays     *  
C     *    representing another complex number and returns the       *  
C     *    complex product.                                          *  
C     *                                                              *  
C     *  Subprograms called: ARMULT, ARSUB, ARADD                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE CMPMUL(AR,AI,BR,BI,CR,CI,WK1,WK2,CR2,D1,D2,WK6,L,RMAX)                       
      
      INTEGER L,I                                                       
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI,D1,D2,CR2,RMAX,WK1,WK2,WK6
      DIMENSION AR(-1:*),AI(-1:*),CR(-1:*),CI(-1:*),WK1(-1:*),WK2(-1:*)
      DIMENSION CR2(-1:*),D1(-1:*),D2(-1:*),WK6(-1:*)
      
      CALL ARMULT(AR,BR,D1,WK6,L,RMAX)                                      
      CALL ARMULT(AI,BI,D2,WK6,L,RMAX)                                      
      CALL ARSUB(D1,D2,CR2,WK1,WK2,L,RMAX)                                      
      CALL ARMULT(AR,BI,D1,WK6,L,RMAX)                                      
      CALL ARMULT(AI,BR,D2,WK6,L,RMAX)                                      
      CALL ARADD(D1,D2,CI,WK1,L,RMAX)                                       
      DO 10 I=-1,L+1                                                   
         CR(I)=CR2(I)                                                    
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ARYDIV                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Returns the double precision complex number   *  
C     *    resulting from the division of four arrays, representing  *  
C     *    two complex numbers.  The number returned will be in one  *  
C     *    of two different forms:  either standard scientific or as *  
C     *    the log (base 10) of the number.                          *  
C     *                                                              *  
C     *  Subprograms called: CONV21, CONV12, EADD, ECPDIV, EMULT.    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ARYDIV(AR,AI,BR,BI,C,L,LNPFQ,RMAX,IBIT)                 
      implicit none      
      INTEGER L,IBIT,REXP,IR10,II10,LNPFQ,ITNMAX
      COMPLEX*16 C,CDUM                                                 
      DOUBLE PRECISION AR,AI,BR,BI,PHI,N1,N2,N3,E1,E2,E3,RR10,RI10,X
      DOUBLE PRECISION AE,BE,X1,X2,DUM1,DUM2,CE,RMAX,TENMAX
      DOUBLE PRECISION DNUM
      DIMENSION AR(-1:*),AI(-1:*),BR(-1:*),BI(-1:*)             
      DIMENSION AE(2,2),BE(2,2),CE(2,2)                                 
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      
      REXP=IBIT/2                                                        
      X=REXP*(AR(L+1)-2)                                                
      RR10=X*LOG10(TWO)/LOG10(TEN)                               
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(AI(L+1)-2)                                                
      RI10=X*LOG10(TWO)/LOG10(TEN)                               
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=SIGN(AR(1)*RMAX*RMAX+AR(2)*RMAX+AR(3),AR(-1))               
      DUM2=SIGN(AI(1)*RMAX*RMAX+AI(2)*RMAX+AI(3),AI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CDUM=DCMPLX(DUM1,DUM2)
      CALL CONV12(CDUM,AE)                                 
      AE(1,2)=AE(1,2)+IR10                                              
      AE(2,2)=AE(2,2)+II10                                              
      X=REXP*(BR(L+1)-2)                                                
      RR10=X*LOG10(TWO)/LOG10(TEN)                               
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(BI(L+1)-2)                                                
      RI10=X*LOG10(TWO)/LOG10(TEN)                               
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=SIGN(BR(1)*RMAX*RMAX+BR(2)*RMAX+BR(3),BR(-1))               
      DUM2=SIGN(BI(1)*RMAX*RMAX+BI(2)*RMAX+BI(3),BI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CDUM=DCMPLX(DUM1,DUM2)
      CALL CONV12(CDUM,BE)                                 
      BE(1,2)=BE(1,2)+IR10                                              
      BE(2,2)=BE(2,2)+II10                                              
      CALL ECPDIV(AE,BE,CE)                                             
      IF (LNPFQ .EQ. 0) THEN                                            
         CALL CONV21(CE,C)                                               
      ELSE                                                              
         CALL EMULT(CE(1,1),CE(1,2),CE(1,1),CE(1,2),N1,E1)               
         CALL EMULT(CE(2,1),CE(2,2),CE(2,1),CE(2,2),N2,E2)               
         CALL EADD(N1,E1,N2,E2,N3,E3)                                    
         N1=CE(1,1)                                                      
         E1=CE(1,2)-CE(2,2)                                              
         X2=CE(2,1)                                                      
*     
*     TENMAX - MAXIMUM SIZE OF EXPONENT OF 10
*     THE FOLLOWING CODE CAN BE USED TO DETERMINE TENMAX, BUT IT
*     WILL LIKELY GENERATE AN IEEE FLOATING POINT UNDERFLOW ERROR
*     ON A SUN WORKSTATION.  REPLACE TENMAX WITH THE VALUE APPROPRIATE
*     FOR YOUR MACHINE.
*     
         TENMAX = 320
         ITNMAX = 1
         DNUM   = 0.1D0
 10      ITNMAX = ITNMAX+1
         DNUM = DNUM*0.1D0
         IF (DNUM .GT. ZERO) GOTO 10
         ITNMAX = ITNMAX-1
         TENMAX = DBLE(ITNMAX)
*     
         IF (E1 .GT. TENMAX) THEN
            X1=TENMAX
         ELSEIF (E1 .LT. -TENMAX) THEN
            X1=ZERO
         ELSE
            X1=N1*(TEN**E1)
         END IF
         IF (X2 .NE. ZERO) THEN
            PHI=ATAN2(X2,X1) 
         ELSE
            PHI = ZERO
         END IF
         C=DCMPLX(HALF*(LOG(N3)+E3*LOG(TEN)),PHI)                 
      ENDIF                                                             
      RETURN                                                            
      END                                                               
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE EMULT                             *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Takes one base and exponent and multiplies it *  
C     *    by another numbers base and exponent to give the product  *  
C     *    in the form of base and exponent.                         *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE EMULT(N1,E1,N2,E2,NF,EF)                               
      
      DOUBLE PRECISION N1,E1,N2,E2,NF,EF
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
     
      NF=N1*N2                                                          
      EF=E1+E2                                                          
      IF (ABS(NF) .GE. TEN) THEN                                    
         NF=NF/TEN                                                    
         EF=EF+ONE                                                     
      ENDIF                                                             
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE EDIV                              *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : returns the solution in the form of base and  *  
C     *    exponent of the division of two exponential numbers.      *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE EDIV(N1,E1,N2,E2,NF,EF)                                
      
      DOUBLE PRECISION N1,E1,N2,E2,NF,EF
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
 
      NF=N1/N2                                                          
      EF=E1-E2                                                          
      IF ((ABS(NF) .LT. ONE) .AND. (NF .NE. ZERO)) THEN             
         NF=NF*TEN                                                    
         EF=EF-ONE                                                     
      ENDIF                                                             
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE EADD                              *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Returns the sum of two numbers in the form    *  
C     *    of a base and an exponent.                                *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE EADD(N1,E1,N2,E2,NF,EF)                                
      
      DOUBLE PRECISION N1,E1,N2,E2,NF,EF,EDIFF
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      
      EDIFF=E1-E2                                                       
      IF (EDIFF .GT. 36.0D0) THEN                                       
         NF=N1                                                           
         EF=E1                                                           
      ELSE IF (EDIFF .LT. -36.0D0) THEN                                 
         NF=N2                                                           
         EF=E2                                                           
      ELSE                                                              
         NF=N1*(TEN**EDIFF)+N2                                        
         EF=E2                                                           
 10      IF (ABS(NF) .LT. TEN) GOTO 20                              
         NF=NF/TEN                                                  
         EF=EF+ONE                                                   
         GOTO 10                                                      
 20      IF ((ABS(NF) .GE. ONE) .OR. (NF .EQ. ZERO)) GOTO 30
         NF=NF*TEN                                                  
         EF=EF-ONE                                                   
         GOTO 20            
      ENDIF                                                             
 30   RETURN                                                            
      END                                                               
        
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ESUB                              *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Returns the solution to the subtraction of    *  
C     *    two numbers in the form of base and exponent.             *  
C     *                                                              *  
C     *  Subprograms called: EADD                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ESUB(N1,E1,N2,E2,NF,EF)                                
      
      DOUBLE PRECISION N1,E1,N2,E2,NF,EF
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      
      CALL EADD(N1,E1,N2*(-ONE),E2,NF,EF)                             
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE CONV12                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Converts a number from complex notation to a  *  
C     *    form of a 2x2 real array.                                 *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE CONV12(CN,CAE)                                         
      
      COMPLEX*16 CN                                                     
      DOUBLE PRECISION CAE
      DIMENSION CAE(2,2)                                                
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
     
      CAE(1,1)=DBLE(CN)                                                
      CAE(1,2)=ZERO                                                    
 10   IF (ABS(CAE(1,1)) .LT. TEN) GOTO 20
      CAE(1,1)=CAE(1,1)/TEN                                        
      CAE(1,2)=CAE(1,2)+ONE                                         
      GOTO 10                                                        
 20   IF ((ABS(CAE(1,1)) .GE. ONE) .OR. (CAE(1,1) .EQ. ZERO))       
     :     GOTO 30                                           
      CAE(1,1)=CAE(1,1)*TEN                                        
      CAE(1,2)=CAE(1,2)-ONE                                         
      GOTO 20                                                        
 30   CAE(2,1)=DIMAG(CN)                                                
      CAE(2,2)=ZERO                                                    
 40   IF (ABS(CAE(2,1)) .LT. TEN) GOTO 50
      CAE(2,1)=CAE(2,1)/TEN                                        
      CAE(2,2)=CAE(2,2)+ONE                                         
      GOTO 40                              
 50   IF ((ABS(CAE(2,1)) .GE. ONE) .OR. (CAE(2,1) .EQ. ZERO))       
     :     GOTO 60
      CAE(2,1)=CAE(2,1)*TEN                                        
      CAE(2,2)=CAE(2,2)-ONE                                         
      GOTO 50                                                    
 60   RETURN                                                            
      END                                                               
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE CONV21                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Converts a number represented in a 2x2 real   *  
C     *    array to the form of a complex number.                    *  
C     *                                                              *  
C     *  Subprograms called: none                                    *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE CONV21(CAE,CN)                                         
      
      DOUBLE PRECISION CAE,DNUM,TENMAX
      INTEGER ITNMAX
      COMPLEX*16 CN                                                     
      DIMENSION CAE(2,2)                                                
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      COMMON/IO/NOUT
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/,/IO/)
*     
*     TENMAX - MAXIMUM SIZE OF EXPONENT OF 10
*     
      ITNMAX = 1
      DNUM   = 0.1D0
 1    ITNMAX = ITNMAX+1
      DNUM = DNUM*0.1D0
      IF (DNUM .GT. ZERO) GOTO 1
      ITNMAX = ITNMAX-2
      TENMAX = DBLE(ITNMAX)
*     
      IF (CAE(1,2) .GT. TENMAX .OR. CAE(2,2) .GT. TENMAX) THEN                  
*     CN=DCMPLX(TENMAX,TENMAX)
         WRITE (NOUT,300) ITNMAX
         STOP
      ELSE IF (CAE(2,2) .LT. -TENMAX) THEN
         CN=DCMPLX(CAE(1,1)*(10**CAE(1,2)),ZERO)                          
      ELSE                                                              
         CN=DCMPLX(CAE(1,1)*(10**CAE(1,2)),CAE(2,1)*(10**CAE(2,2)))      
      ENDIF                                                             
      RETURN                                                            
 300  FORMAT (' ERROR - VALUE OF EXPONENT REQUIRED FOR SUMMATION',
     :     ' WAS LARGER',/,' THAN THE MAXIMUM MACHINE EXPONENT ',
     :     1I3,/,
     :     ' SUGGESTIONS:',/,' 1) RE-RUN USING LNPFQ=1.',/,
     :     ' 2) IF YOU ARE USING A VAX, TRY USING THE',
     :     ' FORTRAN/G_FLOATING OPTION')
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ECPMUL                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Multiplies two numbers which are each         *  
C     *    represented in the form of a two by two array and returns *  
C     *    the solution in the same form.                            *  
C     *                                                              *  
C     *  Subprograms called: EMULT, ESUB, EADD                       *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ECPMUL(A,B,C)                                          
      
      DOUBLE PRECISION A,B,C,N1,E1,N2,E2,C2
      DIMENSION A(2,2),B(2,2),C(2,2),C2(2,2)                            
      
      CALL EMULT(A(1,1),A(1,2),B(1,1),B(1,2),N1,E1)                     
      CALL EMULT(A(2,1),A(2,2),B(2,1),B(2,2),N2,E2)                     
      CALL ESUB(N1,E1,N2,E2,C2(1,1),C2(1,2))                            
      CALL EMULT(A(1,1),A(1,2),B(2,1),B(2,2),N1,E1)                     
      CALL EMULT(A(2,1),A(2,2),B(1,1),B(1,2),N2,E2)                     
      CALL EADD(N1,E1,N2,E2,C(2,1),C(2,2))                              
      C(1,1)=C2(1,1)                                                    
      C(1,2)=C2(1,2)                                                    
      RETURN                                                            
      END                                                               
      
      
C     ****************************************************************  
C     *                                                              *  
C     *                 SUBROUTINE ECPDIV                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Divides two numbers and returns the solution. *  
C     *    All numbers are represented by a 2x2 array.               *  
C     *                                                              *  
C     *  Subprograms called: EADD, ECPMUL, EDIV, EMULT               *  
C     *                                                              *  
C     ****************************************************************  
      
      SUBROUTINE ECPDIV(A,B,C)                                          
      
      DOUBLE PRECISION A,B,C,N1,E1,N2,E2,B2,N3,E3,C2
      DIMENSION A(2,2),B(2,2),C(2,2),B2(2,2),C2(2,2)                    
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
      
      B2(1,1)=B(1,1)                                                    
      B2(1,2)=B(1,2)                                                    
      B2(2,1)=-ONE*B(2,1)                                             
      B2(2,2)=B(2,2)                                                    
      CALL ECPMUL(A,B2,C2)                                              
      CALL EMULT(B(1,1),B(1,2),B(1,1),B(1,2),N1,E1)                     
      CALL EMULT(B(2,1),B(2,2),B(2,1),B(2,2),N2,E2)                     
      CALL EADD(N1,E1,N2,E2,N3,E3)                                      
      CALL EDIV(C2(1,1),C2(1,2),N3,E3,C(1,1),C(1,2))                    
      CALL EDIV(C2(2,1),C2(2,2),N3,E3,C(2,1),C(2,2))                    
      RETURN                                                            
      END                                                               
C     ****************************************************************  
C     *                                                              *  
C     *                   FUNCTION IPREMAX                           *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Predicts the maximum term in the pFq series   *  
C     *    via a simple scanning of arguments.                       *  
C     *                                                              *  
C     *  Subprograms called: none.                                   *  
C     *                                                              *  
C     ****************************************************************  
      
      FUNCTION IPREMAX(A,B,IP,IQ,Z)
      COMPLEX*16 A,B,Z,FACTOR
      INTEGER IP,IQ,J,IPREMAX
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      COMMON/IO/NOUT
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION EXPON,XL,XTERM,XMAX
      DIMENSION A(IP),B(IQ)
      !$omp threadprivate(/CONSTS/,/IO/)
*     
      XTERM=0
      DO 10 J=1,100000
*     
*     Estimate the exponent of the maximum term in the pFq series.
*     
         EXPON=ZERO
         XL=DBLE(J)
         DO 20 I=1,IP
            EXPON=EXPON+DBLE(FACTOR(A(I)+XL-ONE))-DBLE(FACTOR(A(I)-ONE))
 20      CONTINUE
         DO 30 I=1,IQ
            EXPON=EXPON-DBLE(FACTOR(B(I)+XL-ONE))+DBLE(FACTOR(B(I)-ONE))
 30      CONTINUE
         EXPON = EXPON + XL*DBLE(LOG(Z)) - DBLE(FACTOR(DCMPLX(XL,ZERO)))
         XMAX=LOG10(EXP(ONE))*EXPON
         IF ((XMAX .LT. XTERM) .AND. (J .GT. 2)) THEN
            IPREMAX=J
            RETURN
         END IF
         XTERM=MAX(XMAX,XTERM)
 10   CONTINUE
      WRITE (NOUT,*) ' ERROR IN IPREMAX--DID NOT FIND MAXIMUM EXPONENT'
      STOP
      END
C     ****************************************************************  
C     *                                                              *  
C     *                   FUNCTION FACTOR                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : This function is the log of the factorial.    *  
C     *                                                              *  
C     *  Subprograms called: none.                                   *  
C     *                                                              *  
C     ****************************************************************  
      
      FUNCTION FACTOR(Z)
      COMPLEX*16 Z,FACTOR
      DOUBLE PRECISION PI
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      !$omp threadprivate(/CONSTS/)
C     
      IF (((DBLE(Z) .EQ. ONE) .AND. (DIMAG(Z) .EQ. ZERO))
     :     .OR. (ABS(Z) .EQ. ZERO)) THEN
         FACTOR=DCMPLX(ZERO,ZERO)
         RETURN
      END IF
      PI=TWO*TWO*ATAN(ONE)
      FACTOR=HALF*LOG(TWO*PI)+(Z+HALF)*LOG(Z)-Z+(ONE/(12.0D0*Z))*
     :     (ONE-(ONE/(30.D0*Z*Z))*(ONE-(TWO/(7.0D0*Z*Z))))
      RETURN
      END
C     ****************************************************************
C     *                                                              *
C     *                   FUNCTION CGAMMA                            *
C     *                                                              *
C     *                                                              *
C     *  Description : Calculates the complex gamma function.  Based *
C     *     on a program written by F.A. Parpia published in Computer*
C     *     Physics Communications as the `GRASP2' program (public   *
C     *     domain).                                                 *
C     *                                                              *
C     *                                                              *
C     *  Subprograms called: none.                                   *
C     *                                                              *
C     ****************************************************************
      FUNCTION CGAMMA (ARG,LNPFQ)
*     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FIRST,NEGARG
      INTEGER LNPFQ
      COMPLEX*16 CGAMMA,ARG
*     
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION TENTH,PRECIS
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      COMMON/IO/NOUT
      !$omp threadprivate(/CONSTS/,/IO/)
      DIMENSION FN(7),FD(7)
*     
*----------------------------------------------------------------------*
*     *
*     THESE ARE THE BERNOULLI NUMBERS B02, B04, ..., B14, EXPRESSED AS   *
*     RATIONAL NUMBERS. FROM ABRAMOWITZ AND STEGUN, P. 810.              *
*     
      DATA FN/  1.0D 00,   -1.0D 00,    1.0D 00,
     :     -1.0D 00,    5.0D 00, -691.0D 00,
     :     7.0D 00/
      DATA FD/  6.0D 00,   30.0D 00,   42.0D 00,
     :     30.0D 00,   66.0D 00, 2730.0D 00,
     :     6.0D 00/
*     
*----------------------------------------------------------------------*
*     
      DATA HLNTPI/1.0D 00/
*     
      DATA FIRST/.TRUE./
*     
      DATA TENTH/0.1D 00/
*
      DATA PRECIS/0.0D0/
*
      DATA PI/0.0D0/
*
      DATA EXPMAX/0.0D0/
*     
      ARGR = DBLE(ARG)
      ARGI = DIMAG(ARG)
*     
*     ON THE FIRST ENTRY TO THIS ROUTINE, SET UP THE CONSTANTS REQUIRED
*     FOR THE REFLECTION FORMULA (CF. ABRAMOWITZ AND STEGUN 6.1.17) AND
*     STIRLING'S APPROXIMATION (CF. ABRAMOWITZ AND STEGUN 6.1.40).
*     
      IF (FIRST) THEN
         PI= 4.0D0*ATAN(ONE)
*     
*     SET THE MACHINE-DEPENDENT PARAMETERS:
*     
*     TENMAX - MAXIMUM SIZE OF EXPONENT OF 10
*     
         ITNMAX = 1
         DNUM   = TENTH
 10      ITNMAX = ITNMAX+1
         DNUM = DNUM*TENTH
         IF (DNUM .GT. ZERO) GOTO 10
         ITNMAX = ITNMAX-1
         TENMAX = DBLE (ITNMAX)
*     
*     EXPMAX - MAXIMUM SIZE OF EXPONENT OF E
*     
         DNUM = TENTH**ITNMAX
         EXPMAX = -LOG (DNUM)
*     
*     PRECIS - MACHINE PRECISION
*     
         PRECIS = ONE
 20      PRECIS = PRECIS/TWO
         DNUM = PRECIS+ONE
         IF (DNUM .GT. ONE) GOTO 20
         PRECIS = TWO*PRECIS
*     
         HLNTPI = HALF*LOG (TWO*PI)
*     
         DO 30 I = 1,7
            FN(I) = FN(I)/FD(I)
            TWOI = TWO*DBLE (I)
            FN(I) = FN(I)/(TWOI*(TWOI-ONE))
 30      CONTINUE
*     
         FIRST = .FALSE.
*     
      ENDIF
*     
*     CASES WHERE THE ARGUMENT IS REAL
*     
      IF (ARGI .EQ. ZERO) THEN
*     
*     CASES WHERE THE ARGUMENT IS REAL AND NEGATIVE
*     
         IF (ARGR .LE. ZERO) THEN
*     
*     STOP WITH AN ERROR MESSAGE IF THE ARGUMENT IS TOO NEAR A POLE
*     
            DIFF = ABS (DBLE (NINT (ARGR))-ARGR)
            IF (DIFF .LE. TWO*PRECIS) THEN
               WRITE (NOUT,300)
               WRITE (NOUT,301) ARGR,ARGI
               STOP '010801'
            ELSE
*     
*     OTHERWISE USE THE REFLECTION FORMULA (ABRAMOWITZ AND STEGUN 6.1.17)
*     TO ENSURE THAT THE ARGUMENT IS SUITABLE FOR STIRLING'S
*     FORMULA
*     
               ARGUM = PI/(-ARGR*SIN (PI*ARGR))
               IF (ARGUM .LT. ZERO) THEN
                  ARGUM = -ARGUM
                  CLNGI = PI
               ELSE
                  CLNGI = ZERO
               ENDIF
               FACNEG = LOG (ARGUM)
               ARGUR = -ARGR
               NEGARG = .TRUE.
*     
            ENDIF
*     
*     CASES WHERE THE ARGUMENT IS REAL AND POSITIVE
*     
         ELSE
*     
            CLNGI = ZERO
            ARGUR = ARGR
            NEGARG = .FALSE.
*     
         ENDIF
*     
*     USE ABRAMOWITZ AND STEGUN FORMULA 6.1.15 TO ENSURE THAT
*     THE ARGUMENT IN STIRLING'S FORMULA IS GREATER THAN 10
*     
         OVLFAC = ONE
 40      IF (ARGUR .LT. TEN) THEN
            OVLFAC = OVLFAC*ARGUR
            ARGUR = ARGUR+ONE
            GOTO 40
         ENDIF
*     
*     NOW USE STIRLING'S FORMULA TO COMPUTE LOG (GAMMA (ARGUM))
*     
         CLNGR = (ARGUR-HALF)*LOG (ARGUR)-ARGUR+HLNTPI
         FAC = ARGUR
         OBASQ = ONE/(ARGUR*ARGUR)
         DO 50 I = 1,7
            FAC = FAC*OBASQ
            CLNGR = CLNGR+FN(I)*FAC
 50      CONTINUE
*     
*     INCLUDE THE CONTRIBUTIONS FROM THE RECURRENCE AND REFLECTION
*     FORMULAE
*     
         CLNGR = CLNGR-LOG (OVLFAC)
         IF (NEGARG) CLNGR = FACNEG-CLNGR
*     
      ELSE
*     
*     CASES WHERE THE ARGUMENT IS COMPLEX
*     
         ARGUR = ARGR
         ARGUI = ARGI
         ARGUI2 = ARGUI*ARGUI
*     
*     USE THE RECURRENCE FORMULA (ABRAMOWITZ AND STEGUN 6.1.15)
*     TO ENSURE THAT THE MAGNITUDE OF THE ARGUMENT IN STIRLING'S
*     FORMULA IS GREATER THAN 10
*     
         OVLFR = ONE
         OVLFI = ZERO
 60      ARGUM = SQRT (ARGUR*ARGUR+ARGUI2)
         IF (ARGUM .LT. TEN) THEN
            TERMR = OVLFR*ARGUR-OVLFI*ARGUI
            TERMI = OVLFR*ARGUI+OVLFI*ARGUR
            OVLFR = TERMR
            OVLFI = TERMI
            ARGUR = ARGUR+ONE
            GOTO 60
         ENDIF
*     
*     NOW USE STIRLING'S FORMULA TO COMPUTE LOG (GAMMA (ARGUM))
*     
         ARGUR2 = ARGUR*ARGUR
         TERMR = HALF*LOG (ARGUR2+ARGUI2)
         TERMI = ATAN2 (ARGUI,ARGUR)
         CLNGR = (ARGUR-HALF)*TERMR-ARGUI*TERMI-ARGUR+HLNTPI
         CLNGI = (ARGUR-HALF)*TERMI+ARGUI*TERMR-ARGUI
         FAC = (ARGUR2+ARGUI2)**(-2)
         OBASQR = (ARGUR2-ARGUI2)*FAC
         OBASQI = -TWO*ARGUR*ARGUI*FAC
         ZFACR = ARGUR
         ZFACI = ARGUI
         DO 70 I = 1,7
            TERMR = ZFACR*OBASQR-ZFACI*OBASQI
            TERMI = ZFACR*OBASQI+ZFACI*OBASQR
            FAC = FN(I)
            CLNGR = CLNGR+TERMR*FAC
            CLNGI = CLNGI+TERMI*FAC
            ZFACR = TERMR
            ZFACI = TERMI
 70      CONTINUE
*     
*     ADD IN THE RELEVANT PIECES FROM THE RECURRENCE FORMULA
*     
         CLNGR = CLNGR-HALF*LOG (OVLFR*OVLFR+OVLFI*OVLFI)
         CLNGI = CLNGI-ATAN2 (OVLFI,OVLFR)
*     
      ENDIF
      IF (LNPFQ .EQ. 1) THEN
         CGAMMA = DCMPLX(CLNGR,CLNGI)
         RETURN
      END IF
*     
*     NOW EXPONENTIATE THE COMPLEX LOG GAMMA FUNCTION TO GET
*     THE COMPLEX GAMMA FUNCTION
*     
      IF ( (CLNGR .LE.  EXPMAX) .AND.
     :     (CLNGR .GE. -EXPMAX) ) THEN
         FAC = EXP (CLNGR)
      ELSE
         WRITE (NOUT,300)
         WRITE (NOUT,302) CLNGR
         STOP '010802'
      ENDIF
      RESR = FAC*COS (CLNGI)
      RESI = FAC*SIN (CLNGI)
      CGAMMA = DCMPLX(RESR,RESI)
*     
      RETURN
*     
 300  FORMAT (///' ***** ERROR IN SUBROUTINE CGAMMA *****')
 301  FORMAT (' ARGUMENT (',1P,1D14.7,',',1D14.7,') TOO CLOSE TO A',
     :     ' POLE.')
 302  FORMAT (' ARGUMENT TO EXPONENTIAL FUNCTION (',1P,1D14.7,
     :     ') OUT OF RANGE.')
*     
      END
      
C     ****************************************************************  
C     *                                                              *  
C     *                 BLOCK DATA BLDAT1                            *  
C     *                                                              *  
C     *                                                              *  
C     *  Description : Sets of frequently used numbers in a common   *  
C     *    block.  This makes it easier to convert the code to a     *  
C     *    single precision version.                                 *  
C     *                                                              *  
C     ****************************************************************  
      
      BLOCK DATA BLDAT1
C     
      COMMON/IO/NOUT
      COMMON/CONSTS/ZERO,HALF,ONE,TWO,TEN,EPS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,TEN,EPS
      DATA ZERO,HALF,ONE,TWO,TEN,EPS,NOUT/0.0D0,0.5D0,1.0D0,2.0D0,
     :     10.0D0,1.0D-10,6/
      !$omp threadprivate(/CONSTS/)
      !$omp threadprivate(/IO/)
      END
