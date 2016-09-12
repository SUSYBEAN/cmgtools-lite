CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC Buffer for Fjet calls, saves 4-momenta, constraints, weight to be passed to histogram

      INTEGER NFJETBUF
      PARAMETER(NFJETBUF=44)
      INCLUDE 'particles.f' ! need to know how many particles
      INCLUDE 'pdfcontrol.f' ! need to know how many PDF eigenvectors

      DOUBLE PRECISION partBuf(NPART,4,NFJETBUF)
      LOGICAL constrBuf(NFJETBUF) ! stores previously calculated constraints
      INTEGER partflgBuf(NFJETBUF) ! stores the flag for whether to perform jet algorithm
                                   ! 0. no action; 1. combine two partons; 2. set first parton down beam;
                                   ! 3. set second parton down beam 4. set both partons down beam
      DOUBLE PRECISION weightBuf(0:NSMAX,NFJETBUF) ! 0 is central weight
      DOUBLE PRECISION newPdfVec(1:NSMAX) ! gets accumulated PDFs to send to histogram
      INTEGER bufPos ! how many particles we have stored
      COMMON /fjbuf/ partBuf,weightBuf,constrBuf,partflgBuf,bufPos
      SAVE /fjbuf/
