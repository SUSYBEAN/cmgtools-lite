!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLO (To be published)
!!          (See also: Phys. Rev. D79 (2009) 074023 and arXiv:0810.4274)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This package contains the JR09 VFNS NNLO(MSbar) dynamical parton
!! distributions of the nucleon and their associated exact iterative
!! solutions for alpha_s.
!!
!! The sets resulting from displacements in the parameter space with
!! respect to the central value of the NNLO(MSbar) fit along the plus
!! (minus) directions of the (rescaled) eigenvectors of the hessian
!! matrix are included. The tolerance parameter for these displacements
!! was chosen to be T = 4.535 for a total of 1568 fitted points. Since
!! alpha_s was included as a free parameter in the error estimation the
!! use of the provided alpha_s solution for each set is mandatory for
!! uncertainty studies (the difference on alpha_s for different
!! eigenvector sets can be as large as 10% at low scales).
!!
!! The grids are generated for 10^-9 <= x <= 1 and Qo^2 <= Q^2 <= 10^8
!! (GeV^2) where Qo^2 = 0.55 GeV^2 for the NNLO distributions. Outside
!! these ranges the output is either set to NaN or obtained through
!! extrapolation and should NOT be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The routines use a modification of the standard multidimensional
!! linear interpolation routine FINT (CERNLIB E104) distributed as the
!! file 'dfint.f'.
!! The file './JR09VFNNLO.grd', where ./ means the path from the working
!! directory to the file, is read.
!! For questions, comments, problems etc please contact:
!! pjimenez@het.physik.uni-dortmund.de
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLOinit:
!!  Initialization routine of the package to be called (only once)
!!  before using any of the other subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLOx'parton'(x,Q2,set) with 'parton' = uv,dv,gl,ub,db,sb,cb,bb:
!!      Parton distribution 'parton' (times x).
!!        x == Bjorken-x.
!!        Q2 == Q**2 (GeV**2) == Factorization scale
!!                            == Renormalization scale.
!!        set == set to be used.
!! JR09VFNNLOalphas(Q2,set):
!!      Value of alpha_s (no additional 2pi or 4pi factors).
!!        Q2 == Q**2 (GeV**2) == Renormalization scale.
!!        set == set to be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set == Index specifying the set to be used:
!!        0 NNLO(MSbar).
!!        1,2,...,13 set corresponding to a displacement +T with respect
!!                 to the set 0 in the direction of the ith eigenvector.
!!       -1,-2,...,-13 the same for displacements -T.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      block data JR09VFNNLO
       implicit none
       integer shape(2)
       double precision grid(217)
       common /JR09VFNNLOgrid/ grid,shape
       data shape /118,99/
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
      end block data JR09VFNNLO

      subroutine JR09VFNNLOinit(prefix)
       implicit none
       CHARACTER prefix*150
       integer shape(2),i,j,k
       double precision NaN,grid(217),
     &                  alphas(99,-13:13),
     &                  xuv(118,99,-13:13),
     &                  xdv(118,99,-13:13),
     &                  xgl(118,99,-13:13),
     &                  xub(118,99,-13:13),
     &                  xdb(118,99,-13:13),
     &                  xsb(118,99,-13:13),
     &                  xcb(118,99,-13:13),
     &                  xbb(118,99,-13:13)
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOalphasc/ alphas
       common /JR09VFNNLOxuvc/ xuv
       common /JR09VFNNLOxdvc/ xdv
       common /JR09VFNNLOxglc/ xgl
       common /JR09VFNNLOxubc/ xub
       common /JR09VFNNLOxdbc/ xdb
       common /JR09VFNNLOxsbc/ xsb
       common /JR09VFNNLOxcbc/ xcb
       common /JR09VFNNLOxbbc/ xbb
       NaN=0d0
       NaN=0d0/NaN
       open(10,file=prefix(1:len_trim(prefix))//'JR09VFNNLO.grd')
       do 1 k=-13,13
          do 2 j=1,99
             read(10,*) alphas(j,k)
    2     continue
    1  continue
       do 10 k=-13,13
          do 11 j=1,99
             do 12 i=1,118
                read(10,*) xuv(i,j,k)
   12        continue
   11     continue
   10  continue
       do 20 k=-13,13
          do 21 j=1,99
             do 22 i=1,118
                read(10,*) xdv(i,j,k)
   22        continue
   21     continue
   20  continue
       do 30 k=-13,13
          do 31 j=1,99
             do 32 i=1,118
                read(10,*) xgl(i,j,k)
   32        continue
   31     continue
   30  continue
       do 40 k=-13,13
          do 41 j=1,99
             do 42 i=1,118
                read(10,*) xub(i,j,k)
   42        continue
   41     continue
   40  continue
       do 50 k=-13,13
          do 51 j=1,99
             do 52 i=1,118
                read(10,*) xdb(i,j,k)
   52        continue
   51     continue
   50  continue
       do 60 k=-13,13
          do 61 j=1,99
             do 62 i=1,118
                read(10,*) xsb(i,j,k)
   62        continue
   61     continue
   60  continue
       do 70 k=-13,13
          do 71 j=1,99
             do 72 i=1,118
                read(10,*) xcb(i,j,k)
   72        continue
   71     continue
   70  continue
       do 80 k=-13,13
          do 81 j=1,99
             do 82 i=1,118
                read(10,*) xbb(i,j,k)
   82        continue
   81     continue
   80  continue
       close(10)
       do 1000 k=-13,13
          do 1001 j=1,9
             do 1002 i=1,118
                xuv(i,j,k)=NaN
                xdv(i,j,k)=NaN
                xgl(i,j,k)=NaN
                xub(i,j,k)=NaN
                xdb(i,j,k)=NaN
                xsb(i,j,k)=NaN
                xcb(i,j,k)=NaN
                xbb(i,j,k)=NaN
 1002        continue
 1001     continue
 1000  continue
       return
      end subroutine JR09VFNNLOinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOalphas(Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),alphas(99,-13:13),arg,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOalphasc/ alphas
       arg = Q2
       JR09VFNNLOalphas = dfint(1,arg,shape(2),grid(119),alphas(1,set))
       return
      end function JR09VFNNLOalphas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxuv(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xuv(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxuvc/ xuv
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxuv = dfint(2,arg,shape,grid,xuv(1,1,set))
       return
      end function JR09VFNNLOxuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxdv(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xdv(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxdvc/ xdv
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxdv = dfint(2,arg,shape,grid,xdv(1,1,set))
       return
      end function JR09VFNNLOxdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxgl(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xgl(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxglc/ xgl
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxgl = dfint(2,arg,shape,grid,xgl(1,1,set))
       return
      end function JR09VFNNLOxgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxub(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xub(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxubc/ xub
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxub = dfint(2,arg,shape,grid,xub(1,1,set))
       return
      end function JR09VFNNLOxub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxdb(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xdb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxdbc/ xdb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxdb = dfint(2,arg,shape,grid,xdb(1,1,set))
       return
      end function JR09VFNNLOxdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxsb(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xsb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxsbc/ xsb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxsb = dfint(2,arg,shape,grid,xsb(1,1,set))
       return
      end function JR09VFNNLOxsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxcb(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xcb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxcbc/ xcb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxcb = dfint(2,arg,shape,grid,xcb(1,1,set))
       return
      end function JR09VFNNLOxcb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxbb(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(217),xbb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxbbc/ xbb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxbb = dfint(2,arg,shape,grid,xbb(1,1,set))
       return
      end function JR09VFNNLOxbb
