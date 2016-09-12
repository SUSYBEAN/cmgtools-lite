**************************************************************
*
*     evolveNNPDF.f
*
*     Routine  for calculating
*     the value of all xpdfs at x and Q from replica KREP
*     in a (x,Q) point
* 
*     Analogous to the evolvePDF routine of the LHAPDF package
*
***************************************************************

      subroutine evolveNNPDF(X,Q,XPDF)
      IMPLICIT none
*
      INTEGER MXREP
      PARAMETER(MXREP=1e3)
      integer KREPWRAP,NREP
      common/CNREP/NREP
      common/CKREPWRAP/KREPWRAP
*
      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=100,NPQ2=50)
      parameter(NPL=5000)
      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,0:MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2
      
      integer I,IPDF,KREP
      double precision X,Q,Q2,XPDF(-6:6)
*     
*     For interpolation
      integer m,n,nmax,mmax,minq,maxq,midq,maxx,minx,midx
      parameter(m=4,n=2) ! NNPDF2.0 settings
      parameter(nmax=1e3,mmax=1e3)
      double precision dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      integer ix1a(m),ix2a(n),J

*     ---------------------------------------------------------


*     Set correct scale
      Q2=Q**2d0
      
*     Evolved PDF interpolation

      IF ( X.LT.XPDFMIN .OR. X.GT.XPDFMAX ) THEN
         WRITE(6,2000) 
 2000    FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE -- STOP')
         write(6,*) "X= ",X," XMAX, XMIN = ",XPDFMAX,XPDFMIN
         STOP
      ENDIF
*
      IF ( Q2.LT.Q2MIN .OR. Q2.GT.Q2MAX ) THEN
         WRITE(6,2001) 
 2001    FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE -- STOP')
         write(6,*) "Q2 ,Q2MIN, Q2MAX = ",Q2,Q2MIN,Q2MAX
         STOP
      ENDIF
*
*     FIND NEAREST POINTS IN THE GRID
*        
      MINX = 1
      MAXX = NPX+1
 10   CONTINUE
      MIDX = (MINX+MAXX)/2
      IF(X.LT.XG(MIDX)) THEN
         MAXX=MIDX
      ELSE
         MINX=MIDX
      END IF
      IF((MAXX-MINX).GT.1) GO TO 10
      IX = MINX

      MINQ = 1
      MAXQ = NPQ2+1
 20   CONTINUE
      MIDQ = (MINQ+MAXQ)/2
      IF(Q2.LT.Q2G(MIDQ)) THEN
         MAXQ=MIDQ
      ELSE
         MINQ=MIDQ
      END IF
      IF((MAXQ-MINQ).GT.1) GO TO 20
      IQ2 = MINQ

*
*     POLYNOMIAL INTERPOLATION
*        
      
*     Assign grid for interpolation. M, N -> order of polyN interpolation      
      do I=1,M
         if(IX.ge.M/2.and.IX.le.(NPX-M/2)) IX1A(I) = IX - M/2 + I
         if(IX.lt.M/2) IX1A(I) = I
         if(IX.gt.(NPX-M/2)) IX1A(I) = (NPX - M) + I
         
*     Check grids
         if(IX1A(I).le.0.or.IX1A(I).gt.NPX) then
            write(6,*) "Error in grids! "
            write(6,*) "I, IXIA(I) = ",I, IX1A(I)
            call exit(-10)
         endif
      enddo

      do J=1,N
         if(IQ2.ge.N/2.and.IQ2.le.(NPQ2-N/2)) IX2A(J) = IQ2 - N/2 + J
         if(IQ2.lt.N/2) IX2A(J) = J
         if(IQ2.gt.(NPQ2-N/2)) IX2A(J) = (NPQ2 - N) + J
*     Check grids
         if(IX2A(J).le.0.or.IX2A(J).gt.NPQ2) then
            write(6,*) "Error in grids! "
            write(6,*) "J, IXIA(J) = ",J,IX2A(J)
            call exit(-10)
         endif
      enddo
            
*     Define points where to evaluate interpolation
*     Choose between linear or logarithmic (x,Q2) interpolation

      IF(X.LT.XCH)THEN
         X1=dlog(X)          
      ELSE
         X1=X
      ENDIF
      X2=dlog(Q2)

*     Set replica
      KREP=KREPWRAP
      if(KREP.lt.0.or.KREP.gt.NREP) then
         write(6,*) "Invalid value for KREP = ",KREP
         call exit(-10)
      endif
      
      DO IPDF = -6,6,1
         
*     Choose between linear or logarithmic (x,Q2) interpolation
         
         DO I=1,M
            IF(X.LT.XCH)THEN
               X1A(I)= dlog(XG(IX1A(I)))
            ELSE
               X1A(I)= XG(IX1A(I))
            ENDIF
            DO J=1,N
               X2A(J) = dlog(Q2G(IX2A(J)))
               YA(I,J) = XPDFEV(IX1A(I),IX2A(J),IPDF,KREP)
            enddo
         enddo
         
*     2D polynomial interpolation
         call polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
         XPDF(IPDF) = y
         
      enddo                 
*
      RETURN
*
      END

*****************************************************
*     polin2.f
*
*     2D interpolation of arbitrary polinomial order
*     Uses polint
*     Given arrays x1a(1:m) and x2a(1:n) of independent variables,
*     and an m by n array of function values ya(1:m,1:n) tabulated
*     at the grid points defined by x1a,x2a; and given values x1,x2
*     of the independent variable, this routine returns
*     an interpolated function value y with error dy
*
*     Taken from NR fortran
*****************************************************

      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      implicit none
*
      integer m,n,nmax,mmax
      integer j,k
      parameter(nmax=1e3,mmax=1e3)

      DOUBLE PRECISION dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      DOUBLE PRECISION ymtmp(nmax),yntmp(nmax)

      do j=1,m
         do k=1,n
            yntmp(k)=ya(j,k)
         enddo
         call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      enddo
      call polint(x1a,ymtmp,m,x1,y,dy)
*
      return
      end

***********************************************
*     polint.f
*
*     Order N polynomial interpolation using Lagrange's formula
*     as descrived in Numerical Recipees:
*     Given arrays xa and ya each of length n, and given a value
*     x, this routine returns a value y and an error estimate dy.
*     If P(x) is the polynomial of degree N-1 such that 
*     P(xa_i)=ya_i,i=1,...,n, then the returned value is y=P(x)
*     The algorithm used is Neville's algorithm
*******************************************************

      subroutine POLINT(xa,ya,n,x,y,dy)
      implicit none
*
      integer n,NMAX
*     Largest anticipated value of n
      parameter(nmax=1e3)
      DOUBLE PRECISION dy,x,y,xa(nmax),ya(nmax)     
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
         dift=abs(x-xa(i))
         if(dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
 11   enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0) then
               write (*,*) 'failure in polint'
               read(5,*)
            endif
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if(2*ns.lt.(n-m)) then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo

      return
      end

***********************************************

      INTEGER FUNCTION JISEARCH(N,X,Y) 
!                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
!     Dynamical memory allocation                                       
      DOUBLE PRECISION X(*) 
!                                                                       
      MIN=1 
      MAX=N+1 
!                                                                       
   10 CONTINUE 
      MID=(MIN+MAX)/2 
      IF(Y.LT.X(MID)) THEN 
        MAX=MID 
      ELSE 
        MIN=MID 
      END IF 
      IF((MAX-MIN).GT.1) GOTO 10 
!                                                                       
      JISEARCH=MIN 
!                                                                       
      RETURN 
      END                                           



* ------------------------------------------------------


*     Analog routine to initPDF of LHAPDF
*     Selects the replica to be evaluated

      subroutine initNNPDF(IREP)
      implicit none

      integer MXREP,NREP, IREP
      parameter(MXREP=1e3)
      common/CNREP/NREP

      integer KREPWRAP
      common/CKREPWRAP/KREPWRAP

      if(IREP.lt.0.or.IREP.gt.NREP) then
         write(6,*) "Invalid value for IREP in initNNPDF = ",
     1        IREP, NREP
         call exit(-10)
      endif

      KREPWRAP = IREP


      return
      end


* ------------------------------------------------------


