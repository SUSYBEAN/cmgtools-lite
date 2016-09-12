!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GJR08VFNS GRIDS (Phys. Lett. B664 (2008) 133 and arXiv:0801.3618)  !!
!!      (See also: Eur. Phys. J. C53 (2008) 355 and arXiv:0709.0614)  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This package contains the GJR08VFNS LO and NLO(MSbar) dynamical
!! parton distributions of the nucleon and their associated exact
!! iterative solutions for alpha_s. The heavy quark masses used as well
!! as the lambda_QCD parameters for the approximate 'asymptotic'
!! solution for alpha_s are given in the papers.
!!
!! The sets resulting from displacements in the parameter space with
!! respect to the central value of the NLO(MSbar) fit along the plus
!! (minus) directions of the (rescaled) eigenvectors of the hessian
!! matrix are included. The tolerance parameter for these displacements
!! was chosen to be T=4.654 for a total of 1739 fitted points. Since
!! alpha_s was included as a free parameter in the error estimation the
!! use of the provided alpha_s solution for each set is mandatory for
!! uncertainty studies (the difference on alpha_s for different
!! eigenvector sets can be as large as 10% at low scales).
!!
!! The grids are generated for the following range 10^-9 <= x <= 1 and
!! Qo^2 <= Q^2 <= 10^8 (GeV^2) where Qo^2 = 0.5 GeV^2 (0.3 GeV^2) for
!! the NLO (LO) fits. Outside these ranges the output is either set to
!! NaN or obtained through extrapolation and should NOT be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The package has been written in standard Fortran77 and tested using
!! gfortran(v4.2.1) and ifort(v9.1.043) and g77(v3.4.6). The routines
!! use a modification of the standard multidimensional linear
!! interpolation routine FINT (CERNLIB E104) distributed as the file
!! 'dfint.f'. A mimimum of about 24.4 MB of memory is needed to load the
!! grid. For questions, comments, problems etc please contact:
!! pjimenez@het.physik.uni-dortmund.de
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GJR08VFNSinit(set):
!!      Initialization routine of the package to be called before using
!!      any of the other subroutines. It reads the file
!!      './GJR08VFNS.grd', where ./ may be replaced by the path from the
!!      working directory to the grid file.
!!        set == Index specifying the set to be loaded:
!!             0 NLO(MSbar).
!!             1,2,...,13 set corresponding to a displacement +T with
!!                        respect to the set 0 in the direction of the
!!                        ith eigenvector.
!!             -1,-2,...,-13 the same for displacements -T.
!!             15 LO
!! GJR08VFNSx'parton'(x,Q2) with 'parton' = uv,dv,gl,ub,db,sb,cb,bb:
!!      Parton distribution 'parton' (times x) for the set specified in
!!      the last call to GJR08VFNSinit.
!!        x == Bjorken-x.
!!        Q2 == Q**2 (GeV**2) == Factorization scale
!!                            == Renormalization scale.
!! GJR08VFNSalphas(Q2):
!!      Value of alpha_s (no additional 2pi or 4pi factors.) for the set
!!      specified in the last call to GJR08VFNSinit.
!!        Q2 == Q**2 (GeV**2) == Renormalization scale.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GJR08VFNSinit(prefix,set)
       implicit none
       CHARACTER prefix*150
       integer ng(2),init,set,i,j,k,l
       double precision fgrid(-13:15,-5:3,118,99),grid(217),
     &                   xuv(118,99),xdv(118,99),xgl(118,99),
     &                   xub(118,99),xdb(118,99),xsb(118,99),
     &                   xcb(118,99),xbb(118,99),alphas(118,99)
       save init,fgrid
       common /gridc/ grid,ng
       common /xuvc/ xuv
       common /xdvc/ xdv
       common /xglc/ xgl
       common /xubc/ xub
       common /xdbc/ xdb
       common /xsbc/ xsb
       common /xcbc/ xcb
       common /xbbc/ xbb
       common /alphasc/ alphas
       if (init/=1) then
        init=1
        data ng /118,99/
        data grid
     &   /1d-9,1.25d-9,1.6d-9,2d-9,2.5d-9,3.16d-9,4d-9,5d-9,6.3d-9,8d-9,
     &    1d-8,1.25d-8,1.6d-8,2d-8,2.5d-8,3.16d-8,4d-8,5d-8,6.3d-8,8d-8,
     &    1d-7,1.25d-7,1.6d-7,2d-7,2.5d-7,3.16d-7,4d-7,5d-7,6.3d-7,8d-7,
     &    1d-6,1.25d-6,1.6d-6,2d-6,2.5d-6,3.16d-6,4d-6,5d-6,6.3d-6,8d-6,
     &    1d-5,1.25d-5,1.6d-5,2d-5,2.5d-5,3.16d-5,4d-5,5d-5,6.3d-5,8d-5,
     &    1d-4,1.25d-4,1.6d-4,2d-4,2.5d-4,3.16d-4,4d-4,5d-4,6.3d-4,8d-4,
     &    1d-3,1.25d-3,1.6d-3,2d-3,2.5d-3,3.16d-3,4d-3,5d-3,6.3d-3,8d-3,
     &    1d-2,1.25d-2,1.6d-2,2d-2,2.5d-2,3.16d-2,4d-2,5d-2,6.3d-2,8d-2,
     &    0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0,
     &    0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0,
     &    0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0,
     &    0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0,
     &    0.9d0,0.920d0,0.94d0,0.960d0,0.98d0,1d0,
     &    0.3d0,0.31d0,0.35d0,0.375d0,0.4d0,0.45d0,0.5d0,0.51d0,0.525d0,
     &    0.55d0,0.575d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0,
     &    1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0,
     &    1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1,
     &    1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2,
     &    1d3,1.25d3,1.6d3,2d3,2.5d3,3.16d3,4d3,5d3,6.3d3,8d3,
     &    1d4,1.25d4,1.6d4,2d4,2.5d4,3.16d4,4d4,5d4,6.3d4,8d4,
     &    1d5,1.25d5,1.6d5,2d5,2.5d5,3.16d5,4d5,5d5,6.3d5,8d5,
     &    1d6,1.25d6,1.6d6,2d6,2.5d6,3.16d6,4d6,5d6,6.3d6,8d6,
     &    1d7,1.25d7,1.6d7,2d7,2.5d7,3.16d7,4d7,5d7,6.3d7,8d7,1d8/
        open(10,file=prefix(1:len_trim(prefix))//'GJR08VFNS.grd')
        do 4 i=-13,13
         do 3 j=1,118
          do 2 k=1,99
           read(10,*) fgrid(i,-5,j,k),fgrid(i,-4,j,k),
     &                fgrid(i,-3,j,k),fgrid(i,-2,j,k),fgrid(i,-1,j,k),
     &                fgrid(i,0,j,k),fgrid(i,1,j,k),fgrid(i,2,j,k),
     &                fgrid(i,3,j,k)
           do 1 l=-3,3
            if (grid(118+k) < 0.5d0) then
             fgrid(i,l,j,k)=0d0
             fgrid(i,l,j,k)=0d0/fgrid(i,l,j,k)
            end if
    1      continue
    2     continue
    3    continue
    4   continue
        i=14
         do 7 j=1,118
          do 6 k=1,99
           do 5 l=-5,3
             fgrid(i,l,j,k)=0d0
             fgrid(i,l,j,k)=0d0/fgrid(i,l,j,k)
    5      continue
    6     continue
    7    continue
        i=15
         do 10 j=1,118
          do 9 k=1,99
           read(10,*) fgrid(i,-5,j,k),fgrid(i,-4,j,k),
     &                fgrid(i,-3,j,k),fgrid(i,-2,j,k),fgrid(i,-1,j,k),
     &                fgrid(i,0,j,k),fgrid(i,1,j,k),fgrid(i,2,j,k),
     &                fgrid(i,3,j,k)
    9     continue
   10    continue
        close(10)
       end if
       do 15 j=1,118
        do 14 k=1,99
         xuv(j,k) = fgrid(set,1,j,k)
         xdv(j,k) = fgrid(set,2,j,k)
         xgl(j,k) = fgrid(set,0,j,k)
         xub(j,k) = fgrid(set,-1,j,k)
         xdb(j,k) = fgrid(set,-2,j,k)
         xsb(j,k) = fgrid(set,-3,j,k)
         xcb(j,k) = fgrid(set,-4,j,k)
         xbb(j,k) = fgrid(set,-5,j,k)
         alphas(j,k) = fgrid(set,3,j,k)
   14   continue
   15  continue
       return
      end subroutine GJR08VFNSinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxuv(x,Q2)
       implicit none
       integer ng(2)
       double precision xuv(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xuvc/ xuv
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxuv = dfint(2,arg,ng,grid,xuv)
       return
      end function GJR08VFNSxuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxdv(x,Q2)
       implicit none
       integer ng(2)
       double precision xdv(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xdvc/ xdv
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxdv = dfint(2,arg,ng,grid,xdv)
       return
      end function GJR08VFNSxdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxgl(x,Q2)
       implicit none
       integer ng(2)
       double precision xgl(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xglc/ xgl
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxgl = dfint(2,arg,ng,grid,xgl)
       return
      end function GJR08VFNSxgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxub(x,Q2)
       implicit none
       integer ng(2)
       double precision xub(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xubc/ xub
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxub = dfint(2,arg,ng,grid,xub)
       return
      end function GJR08VFNSxub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxdb(x,Q2)
       implicit none
       integer ng(2)
       double precision xdb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xdbc/ xdb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxdb = dfint(2,arg,ng,grid,xdb)
       return
      end function GJR08VFNSxdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxsb(x,Q2)
       implicit none
       integer ng(2)
       double precision xsb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xsbc/ xsb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxsb = dfint(2,arg,ng,grid,xsb)
       return
      end function GJR08VFNSxsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxcb(x,Q2)
       implicit none
       integer ng(2)
       double precision xcb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xcbc/ xcb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxcb = dfint(2,arg,ng,grid,xcb)
       return
      end function GJR08VFNSxcb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxbb(x,Q2)
       implicit none
       integer ng(2)
       double precision xbb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xbbc/ xbb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxbb = dfint(2,arg,ng,grid,xbb)
       return
      end function GJR08VFNSxbb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSalphas(Q2)
       implicit none
       integer ng(2)
       double precision alphas(118,99),grid(217),arg(2),Q2,dfint
       common /gridc/ grid,ng
       common /alphasc/ alphas
       arg(1) = 1d-9
       arg(2) = Q2
       GJR08VFNSalphas = dfint(2,arg,ng,grid,alphas)
       return
      end function GJR08VFNSalphas

