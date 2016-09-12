CCC part of parameters.f
      INTEGER porder,porder_ew,ptype
      INTEGER whichscheme,photflag,DISonQED,EWflag
      DOUBLE PRECISION mur,muf,mufflag,murflag
      double precision m2,Scm
      COMMON/params/ m2,Scm
      COMMON/perturbative/porder,porder_ew
      COMMON/scale/ mur,muf,mufflag,murflag
      COMMON/progmode/ ptype
      COMMON/schmparam/ whichscheme,photflag,DISonQED,EWflag

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE COMPLEX Ieps
      DOUBLE COMPLEX Mzsq
      DOUBLE COMPLEX Mwsq
      DOUBLE COMPLEX cw2
      DOUBLE COMPLEX sw2
      DOUBLE COMPLEX gvd
      DOUBLE COMPLEX gvl
      DOUBLE COMPLEX gvu
      DOUBLE COMPLEX gad
      DOUBLE COMPLEX gal
      DOUBLE COMPLEX gau
      DOUBLE PRECISION alf
      DOUBLE PRECISION ged
      DOUBLE PRECISION gel
      DOUBLE PRECISION geu
      DOUBLE PRECISION delc
      DOUBLE PRECISION dels
      DOUBLE PRECISION Ep
      DOUBLE PRECISION Ml
      DOUBLE PRECISION Mw
      DOUBLE PRECISION Mz
      DOUBLE PRECISION GAMw
      DOUBLE PRECISION GAMz
      DOUBLE PRECISION Mh
      DOUBLE PRECISION Mt
      DOUBLE PRECISION prej0
      DOUBLE PRECISION smax
      DOUBLE PRECISION smin
      DOUBLE PRECISION sllmax
      DOUBLE PRECISION sllmin
      DOUBLE PRECISION sqqmax
      DOUBLE PRECISION sqqmin
      COMMON/PARAM/ Ieps,Mzsq,Mwsq,cw2,sw2,gvu,gvl,gvd,
     &              gau,gal,gad,alf,geu,ged,gel,
     &              delc,dels,Ep,Ml,Mw,Mz,GAMw,GAMz,Mh,Mt,
     &              smax,smin,sllmax,sllmin,sqqmax,sqqmin,prej0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision gvusq,gvdsq,gvlsq,gausq,gadsq,galsq,gvugau,gvdgad,gvlgal
      common/auxcpl/ gvusq,gvdsq,gvlsq,gausq,gadsq,galsq,gvugau,gvdgad,gvlgal

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision GeuSq,GedSq,GelSq,GeuGel,GedGel
      common/gsetctrl/ GeuSq,GedSq,GelSq,GeuGel,GedGel

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE COMPLEX wC1
      DOUBLE COMPLEX wC2
      DOUBLE COMPLEX wC3
      DOUBLE COMPLEX wC4
      DOUBLE COMPLEX wC5
      DOUBLE COMPLEX wC6
      DOUBLE COMPLEX wC7
      DOUBLE COMPLEX wC8
      DOUBLE COMPLEX wC9
      DOUBLE COMPLEX wC10
      DOUBLE COMPLEX wC11
      DOUBLE COMPLEX wC12
      DOUBLE COMPLEX wC13
      DOUBLE COMPLEX wC14
      DOUBLE COMPLEX wC16
      DOUBLE COMPLEX wC18
      DOUBLE COMPLEX wC20
      DOUBLE COMPLEX wC22
      DOUBLE COMPLEX wC24
      DOUBLE COMPLEX wC26
      DOUBLE COMPLEX wC28
      COMMON/VARCONST/ wC1, wC2, wC3, wC4, wC5, wC6, wC7, wC8, wC9
     #, wC10, wC11, wC12, wC13, wC14, wC16, wC18, wC20, wC22, wC24, wC26
     #, wC28

