
CCCCC File IO parameters
      character*128 outputfile
      character*128 inputfile
      character*128 histofile
      character*128 ttmpdir
      character*128 outputdir
      character*128 dummystring

CCCCC VEGAS parameters
      INTEGER ndim,nstart,nincrease,ncomp,maxeval,seed
      INTEGER ntherm,prodim,neval
      INTEGER imethod
      
      double precision epsrel,epsabs

      COMMON/needed/epsrel,epsabs,ndim,nstart,nincrease,seed,
     & ncomp,ntherm,imethod,maxeval,neval

      integer interdim
      parameter (interdim=1)

      double precision integral,error,prob
      dimension integral(interdim),error(interdim),prob(interdim)
      COMMON/results/integral,error,prob
