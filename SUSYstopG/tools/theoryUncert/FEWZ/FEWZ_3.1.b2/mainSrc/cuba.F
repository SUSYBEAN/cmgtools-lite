* cuba.F
* Fortran chooser for the Cuba routines

CCC     Most of cuba.F is redundant, as FEWZ only supports Vegas currently

#define VEGAS 1

C	Passes back Vegas-specific neval, fail for use in updates when 
C	using Vegas

        subroutine Cuba(method, ndim, ncomp, integrand,
     &    epsrel,epsabs,seed,maxeval,nstart,nincrease,statefile,
     &    integral, error, prob, neval, fail, truemaxeval)
        implicit none
        integer method, ndim, ncomp
        external integrand
        double precision integral(*), error(*), prob(*)

        character*7 name(4)
        data name /"Vegas", "Suave", "Divonne", "Cuhre"/

        integer seed
        integer mineval
        integer maxeval, verbose, last
        double precision epsrel, epsabs
        parameter (verbose = 3)
        parameter (last = 0)
        parameter (mineval = 1)

        integer nstart, nincrease, nbatch
        parameter (nbatch=1000)

CCC     iprng=256: Level 1 randomness

        integer iprng
        parameter (iprng=256)

CCC     ismoothing= 0: additional smoothing to importance function applied, 
CCC     ismoothing=16: no additional smoothing
        integer ismoothing
        parameter(ismoothing=16)

CCC     savestate=32: Save Vegas state (we need this)
        integer savestate
        parameter(savestate=32)

        integer gridno
        parameter(gridno=1)

        integer nnew
        double precision flatness
        parameter (nnew = 1000)
        parameter (flatness = 25D0)

        integer nregions, neval, fail, truemaxeval
        integer c

        character*128 statefile

        if( method .eq. VEGAS ) then ! only Vegas currently supported

          call vegas(ndim, ncomp, integrand, 0, !userdata
     &      epsrel,epsabs,verbose+iprng+savestate,seed,mineval,maxeval,
     &      nstart, nincrease, nbatch, gridno, statefile,
     &      neval, fail, integral, error, prob)
          nregions = 1

        else

          print *, "invalid method ", method
          return

        endif

CCC	Print output if this is the last iteration
        if ((neval .gt. truemaxeval) .or. (fail .eq. 0)) then

        print *, "method   =", name(method)
        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(G20.12," +- ",G20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
        end if

        end
