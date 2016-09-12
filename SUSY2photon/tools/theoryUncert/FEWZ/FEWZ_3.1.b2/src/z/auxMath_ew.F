* auxMath.F
* auxiliary functions used by numerical evaluation

c comment the below function if has build-in isnan
      logical function isnan(testnum)
      double precision testnum
      if (testnum.ne.testnum) then
        isnan = .true.
      else
        isnan = .false.
      endif
      end

c comment the below fucntion if has build-in isinf
      logical function isinf(testnum)
      double precision testnum
      if ((1d0/testnum).eq.0d0) then
        isinf = .true.
      else
        isinf = .false.
      endif
      end

      logical function isbad(testnum)
      double precision testnum
!      logical isnan,isinf
!      if isnan(testnum).or.isinf(testnum) then
      if ((testnum*0d0).ne.0d0) then
        isbad = .true.
      else
        isbad = .false.
      endif
      end

************************************************************************

      logical function cisnan(testnum)
      double complex testnum
      logical isnan
      if (isnan(dble(testnum)).or.isnan(dimag(testnum))) then
        cisnan = .true.
      else
        cisnan = .false.
      endif
      end

      logical function cisinf(testnum)
      double complex testnum
      logical isinf
      if (isinf(dble(testnum)).or.isinf(dimag(testnum))) then
        cisinf = .true.
      else
        cisinf = .false.
      endif
      end

      logical function cisbad(testnum)
      double complex testnum
      logical isbad
      if (isbad(dble(testnum)).or.isbad(dimag(testnum))) then
        cisbad = .true.
      else
        cisbad = .false.
      endif
      end

************************************************************************

      double complex function ln(x, isig)
      implicit none

      double precision pi
      parameter(pi=3.14159265358979d0)

      double precision x, isig
      if( x .gt. 0 ) then
        ln = DCMPLX(log(x),0d0)
      else
       ln = dcmplx(log(-x)) + DCMPLX(0D0, sign(pi, isig))
      endif
      end

************************************************************************

      double complex function cln(z, isig)
      implicit none

      double precision pi
      parameter(pi=3.14159265358979d0)

      double complex z
      double precision isig

      if( DIMAG(z) .eq. 0d0 .and. DBLE(z) .le. 0d0 ) then
         cln = log(-z) + DCMPLX(0D0, sign(pi, isig))
      else
         cln = log(z)
      endif
      end

************************************************************************
      
      double complex function denspence(z, isig)
      implicit none
      double complex z
      double precision isig

      double precision pi,pisq,pisqo6
      parameter(pi=3.14159265358979d0,pisq=pi*pi,pisqo6=pisq/6d0)
      double complex cone
      parameter(cone=(1d0,0d0))

c      double precision acc
c      parameter(acc=1d-10)
c      parameter(acc=0d0)

      double complex z1
      double precision az1

      double complex li2series, cln
      external li2series, cln

      z1 = cone - z
      az1 = abs(z1)

c     let the compiler decides when on cut
c      if( isig .eq. 0d0 .and. 
c     &    DIMAG(z) .eq. 0d0 .and. abs(DBLE(z1)) .lt. acc )
c     &    print *, "denspence: argument on cut"

      if( az1 .lt. 1D-15 ) then
        denspence = dcmplx(pisqo6)
      else if( DBLE(z) .lt. .5D0 ) then
        if( abs(z) .lt. 1d0 ) then
          denspence = li2series(z, isig)
        else
          denspence = -dcmplx(pisqo6) -
     &        .5D0*cln(-z, -isig)**2 - li2series(1d0/z, -isig)
        endif
      else
        if( az1 .lt. 1d0 ) then
          denspence = dcmplx(pisqo6) -
     &        cln(z, isig)*cln(z1, -isig) - li2series(z1, -isig)
        else
          denspence = 2d0*dcmplx(pisqo6) +
     &        .5D0*cln(-z1, -isig)**2 - cln(z, isig)*cln(z1, -isig) +
     &        li2series(1d0/z1, isig)
        endif
      endif
      end

************************************************************************

      double complex function li2series(z, isig)
      implicit none

      double complex cone
      parameter(cone=(1d0,0d0))

      double complex z
      double precision isig

      double complex xm, x2, new
      integer j

      double complex cln
      external cln

* these are the even-n Bernoulli numbers, already divided by (n + 1)!
* as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 60, 2}]
      double precision b(30)
      data b /
     &    0.02777777777777777777777777777777777777777778774D0, 
     &    -0.000277777777777777777777777777777777777777777778D0, 
     &    4.72411186696900982615268329554043839758125472D-6, 
     &    -9.18577307466196355085243974132863021751910641D-8, 
     &    1.89788699889709990720091730192740293750394761D-9, 
     &    -4.06476164514422552680590938629196667454705711D-11, 
     &    8.92169102045645255521798731675274885151428361D-13, 
     &    -1.993929586072107568723644347793789705630694749D-14, 
     &    4.51898002961991819165047655285559322839681901D-16, 
     &    -1.035651761218124701448341154221865666596091238D-17, 
     &    2.39521862102618674574028374300098038167894899D-19, 
     &    -5.58178587432500933628307450562541990556705462D-21, 
     &    1.309150755418321285812307399186592301749849833D-22, 
     &    -3.087419802426740293242279764866462431595565203D-24, 
     &    7.31597565270220342035790560925214859103339899D-26, 
     &    -1.740845657234000740989055147759702545340841422D-27, 
     &    4.15763564461389971961789962077522667348825413D-29, 
     &    -9.96214848828462210319400670245583884985485196D-31, 
     &    2.394034424896165300521167987893749562934279156D-32, 
     &    -5.76834735536739008429179316187765424407233225D-34, 
     &    1.393179479647007977827886603911548331732410612D-35, 
     &    -3.372121965485089470468473635254930958979742891D-37, 
     &    8.17820877756210262176477721487283426787618937D-39, 
     &    -1.987010831152385925564820669234786567541858996D-40, 
     &    4.83577851804055089628705937311537820769430091D-42, 
c
     &    -1.1786937248718384326695767537213903193540705623D-43,
     &    2.87709640811725714500196673968866170968836865894D-45,
     &    -7.0320590981560280149649336675824257240019745278D-47,
     &    17208603145033146290899515161658919198419990913276265D-48,
     &    -4.216072390560445491680031819285991615693701229872D-50/

      xm = -cln(cone - z, -isig)
      x2 = xm**2
      li2series = xm - x2/4D0
      do j = 1, 30
        xm = xm*x2
        new = li2series + xm*b(j)
        if( new .eq. li2series ) return
        li2series = new
      enddo
      print *, "li2series: bad convergence"
      end

