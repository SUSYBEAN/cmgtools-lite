c       photon pdf
        
        double precision function retpdfph(x,iset)
C       implicit double precision(a-h,o-z)    
        implicit none
        double precision x,hh,rr,lambda
        integer L
        integer iset
        
        include 'pdfcontrol.f'
        
        hh=hhph
        
        lambda = (x-tau)*hh
        rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))
        
        L =int(rr)
        
        if ((L.gt.0).and.(L.lt.1501)) then
           retpdfph = CCPH(iset,L)+BBPH(iset,L)*x+AAPH(iset,L)*x**2
        else
           retpdfph = 0.0d0
        endif
        
        end
        

c       gluon pdf
	
	double precision function retpdfg(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset
	
        include 'pdfcontrol.f'
	
	hh=hhg
	
	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))
	
	L =int(rr) 
	
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfg = CCG(iset,L)+BBG(iset,L)*x+AAG(iset,L)*x**2
	else
	   retpdfg = 0.0d0
	endif
	
	end
	
	
c       up-valence pdf
	
	double precision function retpdfuv(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset
	
	include 'pdfcontrol.f'
	
	hh=hhuv
	
	lambda = (x-tau)*hh
	
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))
	
	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfuv = CCUV(iset,L)+BBUV(iset,L)*x+AAUV(iset,L)*x**2
	else
	   retpdfuv = 0.0d0
	endif
	
	end

	
	
c       down-valence pdf
	
	double precision function retpdfdv(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'
	
	hh=hhdv
	
	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))
	
	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfdv = CCDV(iset,L)+BBDV(iset,L)*x+AADV(iset,L)*x**2
	else
	   retpdfdv = 0.0d0
	endif

	end

	
	
c       up-sea pdf
	
	double precision function retpdfus(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'
	
	hh=hhus

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfus = CCUS(iset,L)+BBUS(iset,L)*x+AAUS(iset,L)*x**2
	else
	   retpdfus = 0.0d0
	endif

	end

	
	
c       down-sea pdf
	
	double precision function retpdfds(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'
	
	hh=hhds

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfds = CCDS(iset,L)+BBDS(iset,L)*x+AADS(iset,L)*x**2
	else
	   retpdfds = 0.0d0
	endif

	end

	
c       strange - astrange pdf ("valence")
	
	double precision function retpdfsv(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'

	hh=hhsv

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfsv = CCSV(iset,L)+BBSV(iset,L)*x+AASV(iset,L)*x**2
	else
	   retpdfsv = 0.0d0
	endif

	end


c       strange-sea pdf

        double precision function retpdfss(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'

	hh=hhss

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfss = CCSS(iset,L)+BBSS(iset,L)*x+AASS(iset,L)*x**2
	else
	   retpdfss = 0.0d0
	endif

	end

	
c       charm-sea pdf
	
	double precision function retpdfchm(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'

	hh=hhchm

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfchm = CCCHM(iset,L)+BBCHM(iset,L)*x+AACHM(iset,L)*x**2
	else
	   retpdfchm = 0.0d0
	endif

	end
	
	
	
c       bottom-sea pdf
	
	double precision function retpdfbtm(x,iset)
C	implicit double precision(a-h,o-z)    
	implicit none
	double precision x,hh,rr,lambda
	integer L
	integer iset

	include 'pdfcontrol.f'

	hh=hhbtm

	lambda = (x-tau)*hh
	rr = 0.5d0*(1.0d0+dsqrt(dabs(1+4d0*lambda)))

	L =int(rr) 
	if ((L.gt.0).and.(L.lt.1501)) then
	   retpdfbtm = CCBTM(iset,L)+BBBTM(iset,L)*x+AABTM(iset,L)*x**2
	else
	   retpdfbtm = 0.0d0
	endif

	end

CCCCC Make sure no accidental matches before we start

      SUBROUTINE flushCache()
      IMPLICIT NONE
      INCLUDE 'pdfcontrol.f'
      INTEGER i,j,k

      DO i=0,NPDFerror
         DO j=1,NPDFCACHE
            DO k=1,NPDFTYPE
               pdfStore(k,j,i) = 0d0
            ENDDO
            lastX(j,i) = -1d0
         ENDDO
         DO j=1,NLUMCACHE
            lumins(i,j) = 0d0
         ENDDO
      ENDDO
      DO i=1,NLUMCACHE
         lastLumx1(i) = -1d0
         lastLumx2(i) = -1d0
         lastFlag(i) = 0d0
         lastLum(i) = 0d0
      ENDDO
      RETURN
      END

CCCCC Bump old cache back

      SUBROUTINE bumpLumCache()
      IMPLICIT NONE
      INCLUDE 'pdfcontrol.f'
      INTEGER i,iset

      DO i=1,NLUMCACHE-1
         DO iset=0,NPDFERROR
            lumins(iset,NLUMCACHE-i+1) = lumins(iset,NLUMCACHE-i)
         ENDDO
         lastLumx1(NLUMCACHE-i+1) = lastLumx1(NLUMCACHE-i)
         lastLumx2(NLUMCACHE-i+1) = lastLumx2(NLUMCACHE-i)
         lastFlag(NLUMCACHE-i+1) = lastFlag(NLUMCACHE-i)
         lastLum(NLUMCACHE-i+1) = lastLum(NLUMCACHE-i)
      ENDDO

      RETURN
      END

CCCCC Look up luminosity cache, call lumin if necessary

      INTEGER FUNCTION findLumin(x1,x2,fpflag,lum)
      IMPLICIT NONE
      INCLUDE 'pdfcontrol.f'
      DOUBLE PRECISION x1,x2,fpflag,lum,lumin
      EXTERNAL lumin
      INTEGER i

      i = 1
      DO WHILE ((lastLumx1(i).ne.x1.or.
     &           lastLumx2(i).ne.x2.or.
     &           lastFlag(i).ne.fpflag.or.
     &           lastLum(i).ne.lum).and.
     &          i.le.NLUMCACHE)
         i = i + 1

      ENDDO

      IF(i.gt.NLUMCACHE) THEN ! didn't find, we have to compute luminosity
         CALL bumpLumCache()
         DO i=0,NPDFerror
            lumins(i,1) = lumin(x1,x2,fpflag,lum,i)
         ENDDO
         i = 1
         lastLumX1(1) = x1
         lastLumX2(1) = x2
         lastFlag(1) = fpflag
         lastLum(1) = lum            
      ENDIF

      findLumin = i
      END
            

CCCCC Subroutine to decide whether to use previously called PDF or new one

      SUBROUTINE getPDFs(x1,x2,fscale,iset)
      IMPLICIT NONE
      INCLUDE 'pdfcontrol.f'
      DOUBLE PRECISION x1,x2,fscale
      DOUBLE PRECISION pdfval(1:10)
      INTEGER i,iset,j
      double precision retpdfph,retpdfg,retpdfuv,retpdfdv,retpdfus,retpdfds,
     & retpdfsv,retpdfss,retpdfchm,retpdfbtm
      external retpdfph,retpdfg,retpdfuv,retpdfdv,retpdfus,retpdfds,
     &         retpdfsv,retpdfss,retpdfchm,retpdfbtm

      i = 1
      DO WHILE (i.le.NPDFCACHE.and.lastX(i,iset).ne.x1) ! x is not in this cache
         i = i + 1
      ENDDO

      IF(i.gt.NPDFCACHE) THEN ! we didn't find, bump cache down, look up new one
         DO i=1,NPDFCACHE-1
            DO j=1,NPDFTYPE
               pdfStore(j,NPDFCACHE-i+1,iset) = 
     &            pdfStore(j,NPDFCACHE-i,iset)
            ENDDO
            lastX(NPDFCACHE-i+1,iset) = lastX(NPDFCACHE-i,iset)
         ENDDO
         if (dabs(fitscale-fscale).le.0.5d0) then
         pdfStore(1,1,iset) = retpdfuv(x1,iset)
         pdfStore(2,1,iset) = retpdfdv(x1,iset)
         pdfStore(3,1,iset) = retpdfus(x1,iset)
         pdfStore(4,1,iset) = retpdfds(x1,iset)
         pdfStore(5,1,iset) = retpdfsv(x1,iset)
         pdfStore(6,1,iset) = retpdfss(x1,iset)
         pdfStore(7,1,iset) = retpdfchm(x1,iset)
         pdfStore(8,1,iset) = retpdfbtm(x1,iset)
         pdfStore(9,1,iset) = retpdfg(x1,iset)
C         pdfStore(10,1,iset) = retpdfph(x1,iset)
         else
         CALL ObtainPDFs(x1,fscale,iset,pdfval,0)
         pdfStore(1,1,iset) = pdfval(1) !uv
         pdfStore(2,1,iset) = pdfval(2) !dv
         pdfStore(3,1,iset) = pdfval(3) !us
         pdfStore(4,1,iset) = pdfval(4) !ds
         pdfStore(5,1,iset) = pdfval(5) !sv
         pdfStore(6,1,iset) = pdfval(6) !ss
         pdfStore(7,1,iset) = pdfval(7) !chm
         pdfStore(8,1,iset) = pdfval(8) !btm
         pdfStore(9,1,iset) = pdfval(9) !g
C         pdfStore(10,1,iset) = pdfval(10) !ph
         endif
         lastX(1,iset) = x1
         i = 1
      ENDIF

      uval1 = pdfStore(1,i,iset)
      dval1 = pdfStore(2,i,iset)
      usea1 = pdfStore(3,i,iset)
      dsea1 = pdfStore(4,i,iset)
      sval1 = pdfStore(5,i,iset)
      ssea1 = pdfStore(6,i,iset)
      chm1 = pdfStore(7,i,iset)
      btm1 = pdfStore(8,i,iset)
      glu1 = pdfStore(9,i,iset)
C      phot1 = pdfStore(10,i,iset)

      i = 1
      DO WHILE (i.le.NPDFCACHE.and.lastX(i,iset).ne.x2) ! x is not in this cache
         i = i + 1
      ENDDO

      IF(i.gt.NPDFCACHE) THEN ! we didn't find, bump cache down, look up new one
         DO i=1,NPDFCACHE-1
            DO j=1,NPDFTYPE
               pdfStore(j,NPDFCACHE-i+1,iset) = 
     &            pdfStore(j,NPDFCACHE-i,iset)
            ENDDO
            lastX(NPDFCACHE-i+1,iset) = lastX(NPDFCACHE-i,iset)
         ENDDO
         if (dabs(fitscale-fscale).le.0.5d0) then
         pdfStore(1,1,iset) = retpdfuv(x2,iset)
         pdfStore(2,1,iset) = retpdfdv(x2,iset)
         pdfStore(3,1,iset) = retpdfus(x2,iset)
         pdfStore(4,1,iset) = retpdfds(x2,iset)
         pdfStore(5,1,iset) = retpdfsv(x2,iset)
         pdfStore(6,1,iset) = retpdfss(x2,iset)
         pdfStore(7,1,iset) = retpdfchm(x2,iset)
         pdfStore(8,1,iset) = retpdfbtm(x2,iset)
         pdfStore(9,1,iset) = retpdfg(x2,iset)
C         pdfStore(10,1,iset) = retpdfg(x2,iset)
         else
         CALL ObtainPDFs(x2,fscale,iset,pdfval,0)
         pdfStore(1,1,iset) = pdfval(1) !uv
         pdfStore(2,1,iset) = pdfval(2) !dv
         pdfStore(3,1,iset) = pdfval(3) !us
         pdfStore(4,1,iset) = pdfval(4) !ds
         pdfStore(5,1,iset) = pdfval(5) !sv
         pdfStore(6,1,iset) = pdfval(6) !ss
         pdfStore(7,1,iset) = pdfval(7) !chm
         pdfStore(8,1,iset) = pdfval(8) !btm
         pdfStore(9,1,iset) = pdfval(9) !g
C         pdfStore(10,1,iset) = pdfval(10) !ph
         endif
         lastX(1,iset) = x2
         i = 1
      ENDIF

      uval2 = pdfStore(1,i,iset)
      dval2 = pdfStore(2,i,iset)
      usea2 = pdfStore(3,i,iset)
      dsea2 = pdfStore(4,i,iset)
      sval2 = pdfStore(5,i,iset)
      ssea2 = pdfStore(6,i,iset)
      chm2 = pdfStore(7,i,iset)
      btm2 = pdfStore(8,i,iset)
      glu2 = pdfStore(9,i,iset)
C      phot2 = pdfStore(10,i,iset)

      RETURN
      END
