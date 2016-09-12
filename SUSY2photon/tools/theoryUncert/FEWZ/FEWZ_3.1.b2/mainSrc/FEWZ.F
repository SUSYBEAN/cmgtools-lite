CCCCC Reads inputs, defines needed parameters
      subroutine initf(infile, len1, hfile, len2, pdir, len3, isectin)
      implicit none
      character*128 infile,hfile,pdir
      integer len1,len2,len3,isectin

      include 'pdfcontrol.f'
      include 'parameters.f'
      include 'histos.f'
      include 'declar.f'
      include 'cutcontrol.f'
      include 'sectorcontrol.f'

      double precision conv
      double precision pi
      parameter(pi = 3.141592654d0)

CCCCC Histogram parameters
      double precision niter,rescale,dd1,dd2

CCCCC PDF parameters
      integer is,Np
      external fitpdf

CCC   C strings are null terminated, Fortran strings are not
CCC   We have to tell Fortran how long the string is meant to be
CCC   and pad ther rest with blanks
CCC   otherwise we get problems when concatenating

      call Fix_C_Str(infile,len1,inputfile)
      call Fix_C_Str(hfile,len2,histofile)
      call Fix_C_Str(pdir,len3,ttmpdir)

CCC=================================================

CCC   fill Z/W specific parameters

      call inithelp(inputfile)

CCCC  ==============================================
CCC   Define the overall normalization

      conv = 0.3893792922d0*(10.0d0)**9
      norm = 4d0*pi/9d0*alphaqed**2*conv

CCC   Controls cutoff for arrangement of leptons

      ptbuf = 1d-2

CCC=====================================================
CCC   Initialize PDFs: set alpha_s value, pdf error eigenvector number

      call InitPDFs(ttmpdir)
      call CheckPDFRange(-1d0,mlower) ! minimum mass not below what the PDF set allows
      call CheckPDFRange(-1d0,fitscale) ! minimum mass not below what the PDF set allows

CCC   Initialize the histograms and its buffer

      call init_hist(histofile)
      call purgeFJBuf() ! first purge initialize the buffer

CCC   Current version of FEWZ only allows one Vegas component

      ncomp=1

CCC   Take sector number from command line now

      isect = isectin

CCC   Calls the PDF fitting routines

      Np = min(1000,NUMFIT)
      do is=0,NPDFerror
         call fitpdf(fitscale,is,Np)
      enddo

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Calls the integration routine

      subroutine forall(niter, done) !,
C     &     totsect,sectcnt)
      implicit none
      character*128 statefile,statefile2
      integer niter,done,fail
C      integer totsect,sectcnt
      integer numevals
      integer i,j,k,entexist,readtest
      character*128 tempsectnum,tempcontnum

      character*128 pepfx1,pepfx2,pepfx3,pepfx4,pepfx5,
     &     pepfx6,pepfx7
      character*128 peout

      double precision epsabs2,goalerr

      include 'declar.f'
      include 'parameters.f'
      include 'sectorcontrol.f'
      include 'histos.f'
c      include 'sectortiming.f'

      external vegasrout

CCC   Variable for rescaling error to each sector

      epsabs2 = epsabs

CCC   File for vegas saves
CCC   If the program is interupted, the user can restart 
CCC   from the last Vegas iteration: just rerun "fehip"
CCC   The data for the last iteration is stored in: 

      if (isect.lt.10) then
         write(tempsectnum,'(I1)') isect
      else if ((isect.gt.9).and.(isect.lt.100)) then
         write(tempsectnum,'(I2)') isect
      else
         write(tempsectnum,'(I3)') isect
      end if

      if (porder.eq.2) then
         statefile = "vegas_last_iter"//tempsectnum
         statefile2 = "entries"//tempsectnum
         epsabs2 = epsabs*goalerr(isect) ! rescale err target for this sector
      else
         statefile = "vegas_last_iter"

CCC   Auxiliary state file for histograms, pdfs

         statefile2 = "entries"
      end if

CCC   Vegas Integration 
CCC   Takes  minutes (LO,NLO), hours (NNLO inclusive cross-sections)
CCC   or days (NNLO with many cuts/large event rejection)
CCC  
      seed = seed + isect ! make sure sectors stat. indep.

CCC   Hard-code the number of dimensions and integrand components
CCC   The number of components will change as I understand the 
CCC   histogramming

      ndim = 11

CCC    Clear PDF/bin cumulators for the next run
      
      do i=0,NUMBINSTOT
      do j=0,NPDFS
         ENTcum(j,i)=0d0
         ENTsig(j,i)=0d0
      enddo
      enddo

      wgtsum=0d0
      nevit=0

CCC	Max evaluations set to 1 to force Vegas to quit after each 
C	iteration

      ntherm = 0
      imethod = 1

CCC    Read old PDF error file

      entexist = 1
      readtest = 1

      if (niter.eq.0) then ! First iteration, try to find old error file
         open(unit=11, file=statefile2, status='old', 
     &     form='unformatted', iostat=entexist)
         if (entexist.eq.0) then
            read(11, iostat=readtest), 
     &        ((oldcum(i,j),i=0,NPDFS),j=0,NUMBINSTOT)
            read(11, iostat=readtest), 
     &        ((oldsig(i,j),i=0,NPDFS),j=0,NUMBINSTOT)
            rewind(11)
            if (readtest.ne.0) then ! reading failed, start over
               close(11, status='delete')
               open(unit=11, file=statefile)
               close(11, status='delete') ! delete Vegas state file
               open(unit=11, file=statefile2, status='new',
     &           form='unformatted')
               entexist = 0
            end if
         else ! no PDF/histo file, kill any Vegas file and start new
            close(11, status='delete')
            open(unit=11, file=statefile)
            close(11, status='delete') ! delete Vegas state file
            open(unit=11, file=statefile2, status='new',
     &         form='unformatted')
         end if
      else
         open(unit=11, file=statefile2, status='old', 
     &     form='unformatted')
      end if

CCC   Main integration

      call Cuba(imethod,ndim,ncomp,
     &           vegasrout,
     $           0d0,0d0,seed,
     $           1,nstart,nincrease,statefile,
     $           integral,error,prob,neval,fail,maxeval)

CCC   Update averages and sigmas for histograms and PDF errors

      if ((niter.gt.0).or.(entexist.eq.0)) then
         do j=0,NPDFS ! use old method of combining iterations for total
            ENTsig(j,0)=(nevit*ENTsig(j,0)-ENTcum(j,0)**2)
     &      /(nevit-1)
            if (oldsig(j,0).gt.1d-20) then ! previous error is okay 		
            if (ENTsig(j,0).gt.0d0) then ! current error is okay
            oldcum(j,0)=(oldcum(j,0)/oldsig(j,0)
     &                    +ENTcum(j,0)/ENTsig(j,0))
     &                   /(1/oldsig(j,0)+1/ENTsig(j,0))
            oldsig(j,0)=1/(1/oldsig(j,0)+1/ENTsig(j,0))
            else ! current error is bad
               if (j.eq.0) then
                  write(*,*) "Warning: error for this iteration is bad:",
     &             ENTcum(0,0), ENTsig(0,0)
               endif 
            endif ! do nothing if current error is bad
            else ! old error not okay, consider this first iteration
            oldcum(j,0)=ENTcum(j,0)
            oldsig(j,0)=ENTsig(j,0)
            endif
         enddo
         do i=1,NUMBINSTOT ! rest of histogram bins uses old/new method
         do j=0,NPDFS
            ENTsig(j,i)=(nevit*ENTsig(j,i)-ENTcum(j,i)**2)
     &      /(nevit-1)
            if (oldsig(j,i).gt.1d-20.and.dabs(oldcum(j,i)).gt.0d0) then ! previous error is okay 		
            if (ENTsig(j,i).gt.0d0) then ! current error is okay
               if (histmethod.eq.0) then ! use new and more reliable method
               oldsig(j,i)=oldsig(j,i)/oldcum(j,i)**2 ! normalize old error
               oldcum(j,i)=oldcum(j,i)/oldsig(j,i)+ENTcum(j,i)**3/ENTsig(j,i) ! rescale old avg and add new term
               oldsig(j,i)=1/(1/oldsig(j,i)+ENTcum(j,i)**2/ENTsig(j,i)) ! new tot err
               oldcum(j,i)=oldcum(j,i)*oldsig(j,i) ! unnormalize new avg
               oldsig(j,i)=oldsig(j,i)*oldcum(j,i)**2 ! unnormalize new err
               else  ! use old method as for the total for consistency
               oldcum(j,i)=(oldcum(j,i)/oldsig(j,i)
     &                       +ENTcum(j,i)/ENTsig(j,i))
     &                      /(1/oldsig(j,i)+1/ENTsig(j,i))
               oldsig(j,i)=1/(1/oldsig(j,i)+1/ENTsig(j,i))
               endif
            endif ! do nothing if current error is bad
            else ! old error not okay, consider this first iteration
            oldcum(j,i)=ENTcum(j,i)
            oldsig(j,i)=ENTsig(j,i)
            endif
         enddo
         enddo

      else ! file didn't exist, this is first iteration
         do i=0,NUMBINSTOT
         do j=0,NPDFS
            ENTsig(j,i)=(nevit*ENTsig(j,i)-ENTcum(j,i)**2)
     &      /(nevit-1)
            if (ENTsig(j,i).le.0d0) then
               ENTsig(j,i) = 1d-20
            endif
            oldcum(j,i)=ENTcum(j,i)
            oldsig(j,i)=ENTsig(j,i)
         enddo
         enddo
      end if
      
CCC   Write new PDF error file
      write(11), ((oldcum(i,j),i=0,NPDFS),j=0,NUMBINSTOT)
      write(11), ((oldsig(i,j),i=0,NPDFS),j=0,NUMBINSTOT)
      close(11)

CCC	Check to see if error reached, or max evals

      numevals = neval
      
      done = 0

      if ((error(1).lt.epsabs2) .or. 
     &    (error(1).lt.epsrel/100d0*dabs(integral(1)))
     &    .or. (neval.gt.maxeval)) then
	done = 1
      end if

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Writes the output file

      subroutine ending(outfile, len4, outdir, len5, tdiff)
      implicit none
      character*128 outfile,outdir
      integer len4,len5
      character*128 tempsectnum2
      double precision tdiff
      integer noutf

      include 'pdfcontrol.f'
      include 'parameters.f'
      include 'histos.f'
      include 'declar.f'
      include 'cutcontrol.f'
      include 'sectorcontrol.f'

      character*133 fout,fpdf
      character*20 prefix
      character*20 pfx

      double precision SSC
      external SSC

      double precision ENTRIESp(0:NPDFerror,2*NHIST,BINPHIST),
     & histerrp(0:NPDFerror,2*NHIST,BINPHIST)
      integer bincount,counter1,counter2,counter3,i,c1t
      double precision rescale
      
      call Fix_C_Str(outfile,len4,outputfile)
      call Fix_C_Str(outdir,len5,outputdir)

CCC   Transfer the integration results to histogram format for output (histograms and cumulants alternate)

      bincount=1
      do counter1=1,2*NHIST
      do counter2=1,NUMBINS((counter1+1)/2)
      do counter3=0,NPDFerror
         ENTRIESp(counter3,counter1,counter2)
     &   =oldcum(counter3,bincount)
         histerrp(counter3,counter1,counter2)
     &   =oldsig(counter3,bincount)
      enddo
         bincount=bincount+1
      enddo
      enddo

CCC   Write output file with prefix indicating perturbative order

      if (isect.lt.10) then
         write(tempsectnum2,'(I1)') isect
      else if ((isect.gt.9).and.(isect.lt.100)) then
         write(tempsectnum2,'(I2)') isect
      else
         write(tempsectnum2,'(I3)') isect
      end if

      if (porder.ge.2) then
         pfx='NNLO.sect'//tempsectnum2
         prefix=pfx(1:len_trim(pfx))//'.'
      elseif ((porder.eq.1).or.(porder_ew.ge.1)) then
         prefix='NLO.'
      elseif ((porder.eq.0).and.(porder_ew.eq.0)) then
         prefix='LO.'
      endif
      fout=outputdir(1:len_trim(outputdir))//'/'//
     &     prefix(1:len_trim(prefix))//outputfile
      fpdf=outputdir(1:len_trim(outputdir))//'/'//
     &     prefix(1:len_trim(prefix))//'pdf.'//outputfile

CCC   Output all parameters and results

      if (NPDFerror.gt.0) then
         noutf=2
      else
         noutf=1
      endif

      do i=1,noutf

      if (i.eq.1) then
      call unlink(fout)
      open(21,FILE=fout,STATUS='NEW')      
      else
      call unlink(fpdf)
      open(21,FILE=fpdf,STATUS='NEW')
      endif

CC    File 21 should be open before calling this
CCC   Handle boson-specific output
      call endinghelp(i,tdiff)
  
      counter1=1

CCC   Write out histograms

      write(21,*) ''
      write(21,*) '===========   HISTOGRAMS    ================='
      write(21,*) '     bin           weight    numerical error'
C     &          //'      + scale var.     - scale var.'
C     &          //'       + pdf error      - pdf error'

      do while(counter1.le.NHIST)
         if (mod(ISHISTON(counter1),2).eq.1) then
         counter2=1
         c1t=2*counter1-1
c         binSize=(HIGHEDGE(counter1)-LOWEDGE(counter1))
c     &        /NUMBINS(counter1)
         write(21,*)
         write(21,102) HISTNAME(counter1)
         write(21,*)

         do while(counter2.le.NUMBINS(counter1))
            if (i.eq.1) then
               write(21,100) (binBounds(counter2,counter1)
     &                        +binBounds(counter2-1,counter1))/2,
     &                       ENTRIESp(0,c1t,counter2),
     &                       dsqrt(histerrp(0,c1t,counter2))
            else
               do counter3=1,NPDFerror
                  write(21,100) (binBounds(counter2,counter1)
     &                        +binBounds(counter2-1,counter1))/2,
     &                       ENTRIESp(counter3,c1t,counter2),
     &                       dsqrt(histerrp(counter3,c1t,counter2))
               enddo
            endif
            counter2=counter2+1
         enddo
         endif
         if (ISHISTON(counter1).ge.2) then
         counter2=1
         c1t=2*counter1
c         binSize=(HIGHEDGE(counter1)-LOWEDGE(counter1))
c     &        /NUMBINS(counter1)
         write(21,*)
         write(21,103) HISTNAME(counter1)
         write(21,*)

         do while(counter2.le.NUMBINS(counter1))

            if (i.eq.1) then
               write(21,100) (binBounds(counter2,counter1)
     &                        +binBounds(counter2-1,counter1))/2,
     &                       ENTRIESp(0,c1t,counter2),
     &                       dsqrt(histerrp(0,c1t,counter2))
            else
               do counter3=1,NPDFerror
                  write(21,100) (binBounds(counter2,counter1)
     &                        +binBounds(counter2-1,counter1))/2,
     &                       ENTRIESp(counter3,c1t,counter2),
     &                       dsqrt(histerrp(counter3,c1t,counter2))
               enddo
            endif
            counter2=counter2+1
         enddo
         endif
      counter1=counter1+1
      enddo

      close(21)

      enddo

 100  format(F9.2, G17.8, G17.8, G10.2)
 101  format(I3, F17.8, F17.8)
 102  format('  ----   ',A12,'   ----')
 103  format('  ----   ',A12,' (cuml)    ----')


      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Calls the cross section at the appropriate order

      subroutine vegasrout(nd,var,nc,ff,userdata,weight)
      implicit none
      INTEGER nd,nc,counter1,counter2,bincount,counter3
      double precision var(nd),ff(nc),userdata,weight

      external LOe0,NLOe0,NNLOe0

      include 'declar.f'
      include 'parameters.f'
      include 'histos.f'
      include 'pdfcontrol.f'

CCC   store nonleading digits of var(9) not used for later random use

      tmprand = dabs(2*var(9)-1)

CCC   Resetting histograms

      call clearhist()

CCC   Flush PDF/luminosity cache

      call flushCache()

CCC   Call cross section at appropriate order

      if (porder.eq.0) then 
         call LOe0(nd,var,nc,ff)
      elseif (porder.eq.1) then 
         call NLOe0(nd,var,nc,ff)
      elseif (porder.eq.2) then 
         call NNLOe0(nd,var,nc,ff)
      endif

CCC   Fill histograms that are still waiting

      CALL purgeFJBuf()

CCC   Transfer histogram entries to additional components of integrand
CCC   Allows an error estimate on each bin

      bincount=1

      if (dabs(ff(1)).lt.1d100) then ! some occasional floating point infinity bugs; prevent them from reaching Vegas (should be found and fixed later)

      if (ff(1).eq.ff(1)) then ! make sure ff(1) is a number

      ENTcum(0,0)=ENTcum(0,0)+weight*ff(1)
      ENTsig(0,0)=ENTsig(0,0)+weight**2*ff(1)**2

c      if (NPDFerror.gt.0) then
         do counter1=1,nPDFerror
            ENTcum(counter1,0)
     &      =ENTcum(counter1,0)+weight*pdfeigens(counter1)
            ENTsig(counter1,0)
     &      =ENTsig(counter1,0)+weight**2*pdfeigens(counter1)**2
         enddo
c      endif

      do counter1=1,NHIST
      do counter2=1,NUMBINS(counter1)
         if (ENTRIES(0,counter2,counter1).ne.0d0) then
            do counter3=0,nPDFerror
               ENTcum(counter3,bincount)
     &         =ENTcum(counter3,bincount)
     &         +weight*ENTRIES(counter3,counter2,counter1)
               ENTsig(counter3,bincount)
     &         =ENTsig(counter3,bincount)
     &         +weight**2*ENTRIES(counter3,counter2,counter1)**2
            enddo
         endif
         bincount=bincount+1
      enddo
      do counter2=1,NUMBINS(counter1) ! cumulants and histograms alternate in ENTcum/ENTsig
         if (CUMULANTS(0,counter2,counter1).ne.0d0) then
            do counter3=0,nPDFerror
               ENTcum(counter3,bincount)
     &         =ENTcum(counter3,bincount)
     &         +weight*CUMULANTS(counter3,counter2,counter1)
               ENTsig(counter3,bincount)
     &         =ENTsig(counter3,bincount)
     &         +weight**2*CUMULANTS(counter3,counter2,counter1)**2
            enddo
         endif
         bincount=bincount+1
      enddo
      enddo

      endif

      nevit=nevit+1

      else ! infinity, or at least very large from floating point error; zero out for Vegas and do not transfer histogram
         ff(1) = 0d0
      endif

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC
      SUBROUTINE Fix_C_Str (text, length, newtext)
      CHARACTER*(*) text
      INTEGER length
      character*(*) newtext

      INTEGER i


      i = 1
      DO WHILE ((i .LE. length))
         newtext(i:i)=text(i:i)
         i = i + 1
      END DO

C  Pad with blanks to the end of the string

      DO WHILE ((i .gt. length).and.(i.le.128))
        newtext(i:i) = ' '
        i=i+1
      END DO

      END
