***********************************************
*     NNPDFINTevolveLHA
*
*     Routine  for calculating
*     the value of all xpdfs at x and Q from replica KREP
*     in a (x,Q) point
* 
*     Similar to the evolvePDF routine of the LHAPDF package
*     with an extra input variable for the number of replica
*
************************************************

      subroutine NNPDFINTevolveLHA(X,Q,XPDF,KREP)
      IMPLICIT none
*
      INTEGER MXREP
      PARAMETER(MXREP=1e3)
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
      double precision X,Q,Q2,XDIST,QDIST,XPDF(-6:6)
*     
      integer m,n,nmax,mmax,minq,maxq,midq,maxx,minx,midx
      parameter(m=4,n=2)
      parameter(nmax=1e3,mmax=1e3)
      double precision dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      integer ix1a(m),ix2a(n),J

*     Set correct scale
      Q2=Q**2d0
      
*     Evolved PDF interpolation

      IF ( X.LT.XPDFMIN .OR. X.GT.XPDFMAX ) THEN
         WRITE(6,2000) 
 2000    FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE -- STOP')
         write(6,*) "X= ",X," XMAX, XMIN = ",XPDFMAX,XPDFMIN
         !STOP
      ENDIF
*
      IF ( Q2.LT.Q2MIN .OR. Q2.GT.Q2MAX ) THEN
         WRITE(6,2001) 
 2001    FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE -- STOP')
         write(6,*) "Q2 ,Q2MIN, Q2MAX = ",Q2,Q2MIN,Q2MAX
         !STOP
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

***********************************************
*
*
*     NNPDFwrapinit.f
*
*     This routine computes, from a grid of NNPDF
*     evolved pdfs, the value of all xpdfs at x and Q
*     (LHAPDF convention) from replica KREP
*     
*     To be used from an external wrapper
*
************************************************

      subroutine InitNNPDFwrap(wrapfile,NMEM)

*     Loads the array of evolved PDF grids for the
*     LHAPDF-like pdf interpolation

      implicit none
      
      integer MXREP
      parameter(MXREP=1e3)
      integer NREP
      common/CNREP/NREP

      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=100,NPQ2=50)
      parameter(NPL=5000)

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,0:MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2

      integer IPDF,KREP,NLINESDES,ILINES,NXTMP,NQ2TMP,NMEM
      parameter(NLINESDES=11)
      double precision  XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP

      character*53 info
      character*150 wrapfile

*     Read the NNPDF wrapper grid file

      write(*,*)"wrapfile=",wrapfile
      open(unit=10,status="old",file=wrapfile)

*     Read description of grid file
      do ILINES=1,NLINESDES
         read(10,*) info
         write(6,*) info
      enddo
      write(6,*) "  "
         
*     Read max and min extremes of grid
      read(10,*) XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP
 
*     Read and check grids
*     Read the grid in x
      read(10,*) NXTMP
      if(NXTMP.ne.NPX) then
         write(6,*) "Invalid value of NX!"
         call exit(-10)
      endif
      do IX=1,NPX
         read(10,*) XG(IX)
      enddo
*     Read the grid in Q2
      read(10,*) NQ2TMP
      if(NQ2TMP.ne.NPQ2) then
         write(6,*) "Invalid value of NQ2!"
         call exit(-10)
      endif
      do IQ2=1,NPQ2
         read(10,*) Q2G(IQ2)
      enddo

*     Read the number of replicas

      read(10,*) NREP
      NMEM=NREP

*     Read the evolved xpdf grid for each replica
      do KREP=0,NREP
         do IX=1,NPX
            do IQ2=1,NPQ2
               read(10,*) ( XPDFEV(IX,IQ2,IPDF,KREP), IPDF=-6,6,1 )
            enddo
         enddo
      enddo
      close(10)
*
      return
      end


* -----------------------------------------------
