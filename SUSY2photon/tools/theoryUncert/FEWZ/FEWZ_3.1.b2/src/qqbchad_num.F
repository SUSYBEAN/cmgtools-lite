CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbcnum1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff1 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,T
     #coeff*tcoff1)

CCC Sum

      qqbcnum1e0 = tFjet1*tcoff1
      RETURN
      END

      double precision function qqbcnum1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff1 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,T
     #coeff*tcoff1)

CCC Sum

      qqbcnum1e1 = tFjet1*tcoff1
      RETURN
      END

      double precision function qqbcnum1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff1 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,T
     #coeff*tcoff1)

CCC Sum

      qqbcnum1e2 = tFjet1*tcoff1
      RETURN
      END

      double precision function qqbcnum1e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = 0

CCC Coefficients

      tcoff1 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,T
     #coeff*tcoff1)

CCC Sum

      qqbcnum1e3 = tFjet1*tcoff1
      RETURN
      END

      double precision function qqbcnum1e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux1 = 0

CCC Coefficients

      tcoff1 = taux1

CCC Jet functions

      tFjet1 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     #1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,T
     #coeff*tcoff1)

CCC Sum

      qqbcnum1e4 = tFjet1*tcoff1
      RETURN
      END

      double precision function qqbcnum2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux2 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff2 = taux2
      tcoff3 = taux2

CCC Jet functions

      tFjet2 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff3)

CCC Sum

      qqbcnum2e0 = tFjet2*tcoff2+tFjet3*tcoff3
      RETURN
      END

      double precision function qqbcnum2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux2 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff2 = taux2
      tcoff3 = taux2

CCC Jet functions

      tFjet2 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff3)

CCC Sum

      qqbcnum2e1 = tFjet2*tcoff2+tFjet3*tcoff3
      RETURN
      END

      double precision function qqbcnum2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux2 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff2 = taux2
      tcoff3 = taux2

CCC Jet functions

      tFjet2 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff3)

CCC Sum

      qqbcnum2e2 = tFjet2*tcoff2+tFjet3*tcoff3
      RETURN
      END

      double precision function qqbcnum2e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux2 = 0

CCC Coefficients

      tcoff2 = taux2
      tcoff3 = taux2

CCC Jet functions

      tFjet2 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff3)

CCC Sum

      qqbcnum2e3 = tFjet2*tcoff2+tFjet3*tcoff3
      RETURN
      END

      double precision function qqbcnum2e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux2 = 0

CCC Coefficients

      tcoff2 = taux2
      tcoff3 = taux2

CCC Jet functions

      tFjet2 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff2)
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff3)

CCC Sum

      qqbcnum2e4 = tFjet2*tcoff2+tFjet3*tcoff3
      RETURN
      END

      double precision function qqbcnum3e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t5 = -s1l1*s2l2-s1l2*s2l1
      t6 = dlog(z)
      t11 = z**2
      taux3 = -8*z*(1+z)*t5*t6+4.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = 32*z*(1+z)*t5*t6-16.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff4 = taux3
      tcoff5 = taux3
      tcoff6 = taux4
      tcoff7 = taux4

CCC Jet functions

      tFjet4 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff7)

CCC Sum

      qqbcnum3e0 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcof
     #f7
      RETURN
      END

      double precision function qqbcnum3e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t3 = dlog(z)
      t7 = z**2
      taux3 = -8*z*(1+z)*t3+4.D0/3.D0*(z-1)*(4*t7+7*z+4)
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = -96*z*(1+z)*t5*t6+16*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff4 = taux3
      tcoff5 = taux3
      tcoff6 = taux4
      tcoff7 = taux4

CCC Jet functions

      tFjet4 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff7)

CCC Sum

      qqbcnum3e1 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcof
     #f7
      RETURN
      END

      double precision function qqbcnum3e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux3 = 0
      t5 = -s1l2*s2l1+s1l1*s2l2
      t6 = dlog(z)
      t11 = z**2
      taux4 = 64*z*(1+z)*t5*t6-32.D0/3.D0*(z-1)*(4*t11+7*z+4)*t5

CCC Coefficients

      tcoff4 = taux3
      tcoff5 = taux3
      tcoff6 = taux4
      tcoff7 = taux4

CCC Jet functions

      tFjet4 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff7)

CCC Sum

      qqbcnum3e2 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcof
     #f7
      RETURN
      END

      double precision function qqbcnum3e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux3 = 0
      taux4 = 0

CCC Coefficients

      tcoff4 = taux3
      tcoff5 = taux3
      tcoff6 = taux4
      tcoff7 = taux4

CCC Jet functions

      tFjet4 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff7)

CCC Sum

      qqbcnum3e3 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcof
     #f7
      RETURN
      END

      double precision function qqbcnum3e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux3 = 0
      taux4 = 0

CCC Coefficients

      tcoff4 = taux3
      tcoff5 = taux3
      tcoff6 = taux4
      tcoff7 = taux4

CCC Jet functions

      tFjet4 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff4)
      tFjet5 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff5)
      tFjet6 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff*
     #tcoff6)
      tFjet7 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff*
     #tcoff7)

CCC Sum

      qqbcnum3e4 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcof
     #f7
      RETURN
      END

      double precision function qqbcnum4e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux6 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff8 = taux5
      tcoff9 = taux5
      tcoff10 = taux6
      tcoff11 = taux6

CCC Jet functions

      tFjet8 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff11)

CCC Sum

      qqbcnum4e0 = tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*t
     #coff11
      RETURN
      END

      double precision function qqbcnum4e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = -24
      taux6 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff8 = taux5
      tcoff9 = taux5
      tcoff10 = taux6
      tcoff11 = taux6

CCC Jet functions

      tFjet8 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff11)

CCC Sum

      qqbcnum4e1 = tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*t
     #coff11
      RETURN
      END

      double precision function qqbcnum4e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 0
      taux6 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff8 = taux5
      tcoff9 = taux5
      tcoff10 = taux6
      tcoff11 = taux6

CCC Jet functions

      tFjet8 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff11)

CCC Sum

      qqbcnum4e2 = tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*t
     #coff11
      RETURN
      END

      double precision function qqbcnum4e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 0
      taux6 = 0

CCC Coefficients

      tcoff8 = taux5
      tcoff9 = taux5
      tcoff10 = taux6
      tcoff11 = taux6

CCC Jet functions

      tFjet8 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff11)

CCC Sum

      qqbcnum4e3 = tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*t
     #coff11
      RETURN
      END

      double precision function qqbcnum4e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux5 = 0
      taux6 = 0

CCC Coefficients

      tcoff8 = taux5
      tcoff9 = taux5
      tcoff10 = taux6
      tcoff11 = taux6

CCC Jet functions

      tFjet8 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff*
     #tcoff8)
      tFjet9 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff*
     #tcoff9)
      tFjet10 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff10)
      tFjet11 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff11)

CCC Sum

      qqbcnum4e4 = tFjet8*tcoff8+tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*t
     #coff11
      RETURN
      END

      double precision function qqbcnum5e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux8 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff12 = taux7
      tcoff13 = taux7
      tcoff14 = taux8
      tcoff15 = taux8

CCC Jet functions

      tFjet12 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff15)

CCC Sum

      qqbcnum5e0 = tFjet12*tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet
     #15*tcoff15
      RETURN
      END

      double precision function qqbcnum5e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = -24
      taux8 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff12 = taux7
      tcoff13 = taux7
      tcoff14 = taux8
      tcoff15 = taux8

CCC Jet functions

      tFjet12 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff15)

CCC Sum

      qqbcnum5e1 = tFjet12*tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet
     #15*tcoff15
      RETURN
      END

      double precision function qqbcnum5e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 0
      taux8 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff12 = taux7
      tcoff13 = taux7
      tcoff14 = taux8
      tcoff15 = taux8

CCC Jet functions

      tFjet12 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff15)

CCC Sum

      qqbcnum5e2 = tFjet12*tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet
     #15*tcoff15
      RETURN
      END

      double precision function qqbcnum5e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 0
      taux8 = 0

CCC Coefficients

      tcoff12 = taux7
      tcoff13 = taux7
      tcoff14 = taux8
      tcoff15 = taux8

CCC Jet functions

      tFjet12 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff15)

CCC Sum

      qqbcnum5e3 = tFjet12*tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet
     #15*tcoff15
      RETURN
      END

      double precision function qqbcnum5e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux7 = 0
      taux8 = 0

CCC Coefficients

      tcoff12 = taux7
      tcoff13 = taux7
      tcoff14 = taux8
      tcoff15 = taux8

CCC Jet functions

      tFjet12 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff12)
      tFjet13 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff13)
      tFjet14 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff14)
      tFjet15 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff15)

CCC Sum

      qqbcnum5e4 = tFjet12*tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet
     #15*tcoff15
      RETURN
      END

      double precision function qqbcnum6e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux10 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff16 = taux9
      tcoff17 = taux9
      tcoff18 = taux10
      tcoff19 = taux10

CCC Jet functions

      tFjet16 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff19)

CCC Sum

      qqbcnum6e0 = tFjet16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet
     #19*tcoff19
      RETURN
      END

      double precision function qqbcnum6e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = -24
      taux10 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff16 = taux9
      tcoff17 = taux9
      tcoff18 = taux10
      tcoff19 = taux10

CCC Jet functions

      tFjet16 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff19)

CCC Sum

      qqbcnum6e1 = tFjet16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet
     #19*tcoff19
      RETURN
      END

      double precision function qqbcnum6e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 0
      taux10 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff16 = taux9
      tcoff17 = taux9
      tcoff18 = taux10
      tcoff19 = taux10

CCC Jet functions

      tFjet16 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff19)

CCC Sum

      qqbcnum6e2 = tFjet16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet
     #19*tcoff19
      RETURN
      END

      double precision function qqbcnum6e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 0
      taux10 = 0

CCC Coefficients

      tcoff16 = taux9
      tcoff17 = taux9
      tcoff18 = taux10
      tcoff19 = taux10

CCC Jet functions

      tFjet16 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff19)

CCC Sum

      qqbcnum6e3 = tFjet16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet
     #19*tcoff19
      RETURN
      END

      double precision function qqbcnum6e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux9 = 0
      taux10 = 0

CCC Coefficients

      tcoff16 = taux9
      tcoff17 = taux9
      tcoff18 = taux10
      tcoff19 = taux10

CCC Jet functions

      tFjet16 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff16)
      tFjet17 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff17)
      tFjet18 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff18)
      tFjet19 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff19)

CCC Sum

      qqbcnum6e4 = tFjet16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet
     #19*tcoff19
      RETURN
      END

      double precision function qqbcnum7e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux12 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff20 = taux11
      tcoff21 = taux12

CCC Jet functions

      tFjet20 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,Tco
     #eff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,Tco
     #eff*tcoff21)

CCC Sum

      qqbcnum7e0 = tFjet20*tcoff20+tFjet21*tcoff21
      RETURN
      END

      double precision function qqbcnum7e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = -24
      taux12 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff20 = taux11
      tcoff21 = taux12

CCC Jet functions

      tFjet20 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,Tco
     #eff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,Tco
     #eff*tcoff21)

CCC Sum

      qqbcnum7e1 = tFjet20*tcoff20+tFjet21*tcoff21
      RETURN
      END

      double precision function qqbcnum7e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 0
      taux12 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff20 = taux11
      tcoff21 = taux12

CCC Jet functions

      tFjet20 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,Tco
     #eff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,Tco
     #eff*tcoff21)

CCC Sum

      qqbcnum7e2 = tFjet20*tcoff20+tFjet21*tcoff21
      RETURN
      END

      double precision function qqbcnum7e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 0
      taux12 = 0

CCC Coefficients

      tcoff20 = taux11
      tcoff21 = taux12

CCC Jet functions

      tFjet20 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,Tco
     #eff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,Tco
     #eff*tcoff21)

CCC Sum

      qqbcnum7e3 = tFjet20*tcoff20+tFjet21*tcoff21
      RETURN
      END

      double precision function qqbcnum7e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux11 = 0
      taux12 = 0

CCC Coefficients

      tcoff20 = taux11
      tcoff21 = taux12

CCC Jet functions

      tFjet20 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,Tco
     #eff*tcoff20)
      tFjet21 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,Tco
     #eff*tcoff21)

CCC Sum

      qqbcnum7e4 = tFjet20*tcoff20+tFjet21*tcoff21
      RETURN
      END

      double precision function qqbcnum8e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux14 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff22 = taux13
      tcoff23 = taux13
      tcoff24 = taux14
      tcoff25 = taux14

CCC Jet functions

      tFjet22 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff25)

CCC Sum

      qqbcnum8e0 = tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff24+tFjet
     #25*tcoff25
      RETURN
      END

      double precision function qqbcnum8e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = -24
      taux14 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff22 = taux13
      tcoff23 = taux13
      tcoff24 = taux14
      tcoff25 = taux14

CCC Jet functions

      tFjet22 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff25)

CCC Sum

      qqbcnum8e1 = tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff24+tFjet
     #25*tcoff25
      RETURN
      END

      double precision function qqbcnum8e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 0
      taux14 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff22 = taux13
      tcoff23 = taux13
      tcoff24 = taux14
      tcoff25 = taux14

CCC Jet functions

      tFjet22 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff25)

CCC Sum

      qqbcnum8e2 = tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff24+tFjet
     #25*tcoff25
      RETURN
      END

      double precision function qqbcnum8e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 0
      taux14 = 0

CCC Coefficients

      tcoff22 = taux13
      tcoff23 = taux13
      tcoff24 = taux14
      tcoff25 = taux14

CCC Jet functions

      tFjet22 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff25)

CCC Sum

      qqbcnum8e3 = tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff24+tFjet
     #25*tcoff25
      RETURN
      END

      double precision function qqbcnum8e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux13 = 0
      taux14 = 0

CCC Coefficients

      tcoff22 = taux13
      tcoff23 = taux13
      tcoff24 = taux14
      tcoff25 = taux14

CCC Jet functions

      tFjet22 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff22)
      tFjet23 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff23)
      tFjet24 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff24)
      tFjet25 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff25)

CCC Sum

      qqbcnum8e4 = tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff24+tFjet
     #25*tcoff25
      RETURN
      END

      double precision function qqbcnum9e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux16 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff26 = taux15
      tcoff27 = taux15
      tcoff28 = taux16
      tcoff29 = taux16

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff26)
      t3 = (1-z)*l2+z
      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z/t3,t3,Tcoef
     #f*tcoff27)
      t3 = (1-z)*l2+z
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff28)
      t3 = (1-z)*l2+z
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z/t3,t3,Tcoef
     #f*tcoff29)

CCC Sum

      qqbcnum9e0 = tFjet26*tcoff26+tFjet27*tcoff27+tFjet28*tcoff28+tFjet
     #29*tcoff29
      RETURN
      END

      double precision function qqbcnum9e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = -24
      taux16 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff26 = taux15
      tcoff27 = taux15
      tcoff28 = taux16
      tcoff29 = taux16

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff26)
      t3 = (1-z)*l2+z
      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z/t3,t3,Tcoef
     #f*tcoff27)
      t3 = (1-z)*l2+z
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff28)
      t3 = (1-z)*l2+z
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z/t3,t3,Tcoef
     #f*tcoff29)

CCC Sum

      qqbcnum9e1 = tFjet26*tcoff26+tFjet27*tcoff27+tFjet28*tcoff28+tFjet
     #29*tcoff29
      RETURN
      END

      double precision function qqbcnum9e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 0
      taux16 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff26 = taux15
      tcoff27 = taux15
      tcoff28 = taux16
      tcoff29 = taux16

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff26)
      t3 = (1-z)*l2+z
      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z/t3,t3,Tcoef
     #f*tcoff27)
      t3 = (1-z)*l2+z
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff28)
      t3 = (1-z)*l2+z
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z/t3,t3,Tcoef
     #f*tcoff29)

CCC Sum

      qqbcnum9e2 = tFjet26*tcoff26+tFjet27*tcoff27+tFjet28*tcoff28+tFjet
     #29*tcoff29
      RETURN
      END

      double precision function qqbcnum9e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 0
      taux16 = 0

CCC Coefficients

      tcoff26 = taux15
      tcoff27 = taux15
      tcoff28 = taux16
      tcoff29 = taux16

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff26)
      t3 = (1-z)*l2+z
      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z/t3,t3,Tcoef
     #f*tcoff27)
      t3 = (1-z)*l2+z
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff28)
      t3 = (1-z)*l2+z
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z/t3,t3,Tcoef
     #f*tcoff29)

CCC Sum

      qqbcnum9e3 = tFjet26*tcoff26+tFjet27*tcoff27+tFjet28*tcoff28+tFjet
     #29*tcoff29
      RETURN
      END

      double precision function qqbcnum9e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux15 = 0
      taux16 = 0

CCC Coefficients

      tcoff26 = taux15
      tcoff27 = taux15
      tcoff28 = taux16
      tcoff29 = taux16

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet26 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff26)
      t3 = (1-z)*l2+z
      tFjet27 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z/t3,t3,Tcoef
     #f*tcoff27)
      t3 = (1-z)*l2+z
      tFjet28 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff28)
      t3 = (1-z)*l2+z
      tFjet29 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z/t3,t3,Tcoef
     #f*tcoff29)

CCC Sum

      qqbcnum9e4 = tFjet26*tcoff26+tFjet27*tcoff27+tFjet28*tcoff28+tFjet
     #29*tcoff29
      RETURN
      END

      double precision function qqbcnum10e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 24*s1l2*s2l1+24*s1l1*s2l2
      taux18 = -96*s1l2*s2l1+96*s1l1*s2l2

CCC Coefficients

      tcoff30 = taux17
      tcoff31 = taux18

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff30)
      t3 = (1-z)*l2+z
      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff31)

CCC Sum

      qqbcnum10e0 = tFjet30*tcoff30+tFjet31*tcoff31
      RETURN
      END

      double precision function qqbcnum10e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = -24
      taux18 = 288*s1l2*s2l1-288*s1l1*s2l2

CCC Coefficients

      tcoff30 = taux17
      tcoff31 = taux18

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff30)
      t3 = (1-z)*l2+z
      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff31)

CCC Sum

      qqbcnum10e1 = tFjet30*tcoff30+tFjet31*tcoff31
      RETURN
      END

      double precision function qqbcnum10e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 0
      taux18 = -192*s1l2*s2l1+192*s1l1*s2l2

CCC Coefficients

      tcoff30 = taux17
      tcoff31 = taux18

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff30)
      t3 = (1-z)*l2+z
      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff31)

CCC Sum

      qqbcnum10e2 = tFjet30*tcoff30+tFjet31*tcoff31
      RETURN
      END

      double precision function qqbcnum10e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 0
      taux18 = 0

CCC Coefficients

      tcoff30 = taux17
      tcoff31 = taux18

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff30)
      t3 = (1-z)*l2+z
      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff31)

CCC Sum

      qqbcnum10e3 = tFjet30*tcoff30+tFjet31*tcoff31
      RETURN
      END

      double precision function qqbcnum10e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux17 = 0
      taux18 = 0

CCC Coefficients

      tcoff30 = taux17
      tcoff31 = taux18

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet30 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,t3,z/t3,Tcoef
     #f*tcoff30)
      t3 = (1-z)*l2+z
      tFjet31 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,t3,z/t3,Tcoef
     #f*tcoff31)

CCC Sum

      qqbcnum10e4 = tFjet30*tcoff30+tFjet31*tcoff31
      RETURN
      END

      double precision function qqbcnum11e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = 128*s1l1*s2l1+64*s2l1+64*s1l1
      taux20 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff32 = taux19
      tcoff33 = taux19
      tcoff34 = taux20
      tcoff35 = taux20

CCC Jet functions

      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff35)

CCC Sum

      qqbcnum11e0 = tFjet32*tcoff32+tFjet33*tcoff33+tFjet34*tcoff34+tFje
     #t35*tcoff35
      RETURN
      END

      double precision function qqbcnum11e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = -96*s1l1-96*s2l1-64*s1l1*s2l1
      taux20 = 896*s2l1-896*s1l1

CCC Coefficients

      tcoff32 = taux19
      tcoff33 = taux19
      tcoff34 = taux20
      tcoff35 = taux20

CCC Jet functions

      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff35)

CCC Sum

      qqbcnum11e1 = tFjet32*tcoff32+tFjet33*tcoff33+tFjet34*tcoff34+tFje
     #t35*tcoff35
      RETURN
      END

      double precision function qqbcnum11e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = 128*s1l1*s2l1+96*s2l1+96*s1l1
      taux20 = -896*s2l1+896*s1l1

CCC Coefficients

      tcoff32 = taux19
      tcoff33 = taux19
      tcoff34 = taux20
      tcoff35 = taux20

CCC Jet functions

      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff35)

CCC Sum

      qqbcnum11e2 = tFjet32*tcoff32+tFjet33*tcoff33+tFjet34*tcoff34+tFje
     #t35*tcoff35
      RETURN
      END

      double precision function qqbcnum11e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = -64*s2l1-64*s1l1
      taux20 = 384*s2l1-384*s1l1

CCC Coefficients

      tcoff32 = taux19
      tcoff33 = taux19
      tcoff34 = taux20
      tcoff35 = taux20

CCC Jet functions

      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff35)

CCC Sum

      qqbcnum11e3 = tFjet32*tcoff32+tFjet33*tcoff33+tFjet34*tcoff34+tFje
     #t35*tcoff35
      RETURN
      END

      double precision function qqbcnum11e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux19 = 0
      taux20 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff32 = taux19
      tcoff33 = taux19
      tcoff34 = taux20
      tcoff35 = taux20

CCC Jet functions

      tFjet32 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,z,Tcoeff
     #*tcoff32)
      tFjet33 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,z,1.D0,Tcoeff
     #*tcoff33)
      tFjet34 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,z,Tcoeff
     #*tcoff34)
      tFjet35 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.
     #D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,z,1.D0,Tcoeff
     #*tcoff35)

CCC Sum

      qqbcnum11e4 = tFjet32*tcoff32+tFjet33*tcoff33+tFjet34*tcoff34+tFje
     #t35*tcoff35
      RETURN
      END

      double precision function qqbcnum12e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff36 = taux21

CCC Jet functions

      tFjet36 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff36)

CCC Sum

      qqbcnum12e0 = tFjet36*tcoff36
      RETURN
      END

      double precision function qqbcnum12e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = 896*s2l1-896*s1l1

CCC Coefficients

      tcoff36 = taux21

CCC Jet functions

      tFjet36 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff36)

CCC Sum

      qqbcnum12e1 = tFjet36*tcoff36
      RETURN
      END

      double precision function qqbcnum12e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = -896*s2l1+896*s1l1

CCC Coefficients

      tcoff36 = taux21

CCC Jet functions

      tFjet36 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff36)

CCC Sum

      qqbcnum12e2 = tFjet36*tcoff36
      RETURN
      END

      double precision function qqbcnum12e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = 384*s2l1-384*s1l1

CCC Coefficients

      tcoff36 = taux21

CCC Jet functions

      tFjet36 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff36)

CCC Sum

      qqbcnum12e3 = tFjet36*tcoff36
      RETURN
      END

      double precision function qqbcnum12e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux21 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff36 = taux21

CCC Jet functions

      tFjet36 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff36)

CCC Sum

      qqbcnum12e4 = tFjet36*tcoff36
      RETURN
      END

      double precision function qqbcnum13e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux22 = 128*s1l1*s2l1+64*s2l1+64*s1l1
      taux23 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff37 = taux22
      tcoff38 = taux23

CCC Jet functions

      tFjet37 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff37)
      tFjet38 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff38)

CCC Sum

      qqbcnum13e0 = tFjet37*tcoff37+tFjet38*tcoff38
      RETURN
      END

      double precision function qqbcnum13e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux22 = -96*s1l1-96*s2l1-64*s1l1*s2l1
      taux23 = 896*s2l1-896*s1l1

CCC Coefficients

      tcoff37 = taux22
      tcoff38 = taux23

CCC Jet functions

      tFjet37 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff37)
      tFjet38 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff38)

CCC Sum

      qqbcnum13e1 = tFjet37*tcoff37+tFjet38*tcoff38
      RETURN
      END

      double precision function qqbcnum13e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux22 = 128*s1l1*s2l1+96*s2l1+96*s1l1
      taux23 = -896*s2l1+896*s1l1

CCC Coefficients

      tcoff37 = taux22
      tcoff38 = taux23

CCC Jet functions

      tFjet37 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff37)
      tFjet38 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff38)

CCC Sum

      qqbcnum13e2 = tFjet37*tcoff37+tFjet38*tcoff38
      RETURN
      END

      double precision function qqbcnum13e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux22 = -64*s2l1-64*s1l1
      taux23 = 384*s2l1-384*s1l1

CCC Coefficients

      tcoff37 = taux22
      tcoff38 = taux23

CCC Jet functions

      tFjet37 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff37)
      tFjet38 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff38)

CCC Sum

      qqbcnum13e3 = tFjet37*tcoff37+tFjet38*tcoff38
      RETURN
      END

      double precision function qqbcnum13e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      taux22 = 0
      taux23 = -256*s2l1+256*s1l1

CCC Coefficients

      tcoff37 = taux22
      tcoff38 = taux23

CCC Jet functions

      tFjet37 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff37)
      tFjet38 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff38)

CCC Sum

      qqbcnum13e4 = tFjet37*tcoff37+tFjet38*tcoff38
      RETURN
      END

      double precision function qqbcnum14e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t1 = s1l1**2
      taux24 = 8.D0/3.D0*t1+4.D0/3.D0+8.D0/3.D0*s1l1
      taux25 = -32.D0/3.D0*s1l1-16.D0/3.D0

CCC Coefficients

      tcoff39 = taux24
      tcoff40 = taux25

CCC Jet functions

      tFjet39 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff39)
      tFjet40 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff40)

CCC Sum

      qqbcnum14e0 = tFjet39*tcoff39+tFjet40*tcoff40
      RETURN
      END

      double precision function qqbcnum14e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t1 = s1l1**2
      taux24 = t1*NF/6+21.D0/4.D0*t1+NF/12+31.D0/24.D0+s1l1*NF/6+21.D0/4
     #.D0*s1l1
      taux25 = -(2*s1l1+1)*(2*NF-33)/6

CCC Coefficients

      tcoff39 = taux24
      tcoff40 = taux25

CCC Jet functions

      tFjet39 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff39)
      tFjet40 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff40)

CCC Sum

      qqbcnum14e1 = tFjet39*tcoff39+tFjet40*tcoff40
      RETURN
      END

      double precision function qqbcnum14e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t1 = 0.3141592653589793D1**2
      t2 = s1l1**2
      taux24 = -103.D0/36.D0*t1*t2+27.D0/2.D0*t2+7.D0/9.D0*t2*NF-103.D0/
     #72.D0*t1+33.D0/8.D0+11.D0/36.D0*NF-103.D0/36.D0*t1*s1l1+27.D0/2.D0
     #*s1l1+7.D0/9.D0*s1l1*NF
      t4 = 0.3141592653589793D1**2
      taux25 = -(2*s1l1+1)*(10*NF-103*t4-81)/18

CCC Coefficients

      tcoff39 = taux24
      tcoff40 = taux25

CCC Jet functions

      tFjet39 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff39)
      tFjet40 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff40)

CCC Sum

      qqbcnum14e2 = tFjet39*tcoff39+tFjet40*tcoff40
      RETURN
      END

      double precision function qqbcnum14e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t6 = s1l1**2
      t9 = 0.3141592653589793D1**2
      t13 = t9*s1l1
      t15 = t6*t9
      taux24 = 269.D0/216.D0*NF-17.D0/36.D0*z3+1211.D0/72.D0*s1l1+353.D0
     #/108.D0*s1l1*NF+353.D0/108.D0*t6*NF-55.D0/144.D0*t9-17.D0/18.D0*s1
     #l1*z3-29.D0/8.D0*t13-11.D0/36.D0*t15*NF-11.D0/36.D0*t13*NF-17.D0/1
     #8.D0*t6*z3+1211.D0/72.D0*t6-29.D0/8.D0*t15-11.D0/72.D0*t9*NF+239.D
     #0/144.D0
      t5 = 0.3141592653589793D1**2
      taux25 = (2*s1l1+1)*(-202*NF+204*z3-1071*t5+2811+66*t5*NF)/108

CCC Coefficients

      tcoff39 = taux24
      tcoff40 = taux25

CCC Jet functions

      tFjet39 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff39)
      tFjet40 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff40)

CCC Sum

      qqbcnum14e3 = tFjet39*tcoff39+tFjet40*tcoff40
      RETURN
      END

      double precision function qqbcnum14e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)

CCC Auxiliary functions

      t2 = 0.3141592653589793D1**2
      t6 = t2**2
      t8 = s1l1*NF
      t10 = s1l1**2
      t11 = t10*NF
      t13 = t2*s1l1
      t15 = t2*t10
      t39 = -13531.D0/864.D0-77.D0/54.D0*t13*NF-121.D0/216.D0*t2*NF-13.D
     #0/9.D0*t11*z3-13.D0/9.D0*t8*z3+119.D0/135.D0*t6*t10+119.D0/135.D0*
     #t6*s1l1+41.D0/2.D0*s1l1*z3+41.D0/2.D0*t10*z3-6265.D0/432.D0*t10-13
     #.D0/18.D0*NF*z3
      taux24 = 5423.D0/1296.D0*NF-49.D0/144.D0*t2+193.D0/18.D0*z3-6265.D
     #0/432.D0*s1l1+119.D0/270.D0*t6+7541.D0/648.D0*t8+7541.D0/648.D0*t1
     #1-155.D0/36.D0*t13-155.D0/36.D0*t15-77.D0/54.D0*t15*NF+t39
      t7 = 0.3141592653589793D1**2
      t9 = t7**2
      taux25 = (2*s1l1+1)*(9360*NF*z3-17270*NF-151200*z3-45810*t7-5712*t
     #9+3300*t7*NF+242205)/3240

CCC Coefficients

      tcoff39 = taux24
      tcoff40 = taux25

CCC Jet functions

      tFjet39 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,3.D0,1.D0,1.D0,
     #Tcoeff*tcoff39)
      tFjet40 = Fjet(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
     #,1.D0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,6.D0,1.D0,1.D0,
     #Tcoeff*tcoff40)

CCC Sum

      qqbcnum14e4 = tFjet39*tcoff39+tFjet40*tcoff40
      RETURN
      END

      double precision function qqbcnum15e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp-1
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = -64*l1+64*t10+64-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux26 = 2.D0/3.D0*(1+t2)/yp*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1
     #l1+(t12*s2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = -zp+t6-l1
      t10 = -l1+t6+1-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = 1-zp
      taux27 = 2.D0/3.D0*(1+t2)/yp*((256*t7*s2l2+256*t10*t11)*s1l1+(-256
     #*t7*s2l1-256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff41 = taux26
      tcoff42 = taux26
      tcoff43 = taux27
      tcoff44 = -taux27

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f44)

CCC Sum

      qqbcnum15e0 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qqbcnum15e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t6 = zp*l1
      t8 = -64*l1+64*t6+64
      t25 = zp**2
      t26 = l1**2
      taux26 = 2.D0/3.D0*(1+t2)/yp*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(zp-1
     #)*(-1+l1)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64
     #*zp*(-t25*t26+2*t25*l1+2*zp-2*t25+2*zp*t26-1-t26-2*t6))
      t2 = (1-yp)**2
      t6 = 3*zp
      t7 = 3*l1
      t8 = zp*l1
      t9 = 3*t8
      t10 = t6+t7-t9+1
      t13 = t7-t9-3+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = -3-l1+t8+t6
      taux27 = 2.D0/3.D0*(1+t2)/yp*((256*t10*s2l2+256*t13*t14)*s1l1+(-25
     #6*t10*s2l1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff41 = taux26
      tcoff42 = taux26
      tcoff43 = taux27
      tcoff44 = -taux27

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f44)

CCC Sum

      qqbcnum15e1 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qqbcnum15e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = (1-yp)**2
      t8 = l1**2
      t11 = zp**2
      taux26 = -128.D0/3.D0*(1+t2)/yp*zp*(2*zp*l1-2*zp*t8+t11*t8+t8-2*l1
     #+1)
      t2 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = -3-l1+t6-t7
      t11 = -l1+t6+1-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = 1+3*l1-3*t6-t7
      taux27 = 2.D0/3.D0*(1+t2)/yp*((256*t8*s2l2+512*t11*t12)*s1l1+(-256
     #*t8*s2l1-512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff41 = taux26
      tcoff42 = taux26
      tcoff43 = taux27
      tcoff44 = -taux27

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f44)

CCC Sum

      qqbcnum15e2 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qqbcnum15e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux26 = 0
      t2 = (1-yp)**2
      t7 = 1-l1+zp*l1
      taux27 = 2.D0/3.D0*(1+t2)/yp*(512*t7*s2l2*s1l1-512*t7*s2l1*s1l2+51
     #2*t7*(-s1l2-s2l2-zp)*s2l1-512*t7*(-s1l1-s2l1-zp)*s2l2)

CCC Coefficients

      tcoff41 = taux26
      tcoff42 = taux26
      tcoff43 = taux27
      tcoff44 = -taux27

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f44)

CCC Sum

      qqbcnum15e3 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qqbcnum15e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux26 = 0
      taux27 = 0

CCC Coefficients

      tcoff41 = taux26
      tcoff42 = taux26
      tcoff43 = taux27
      tcoff44 = -taux27

CCC Jet functions

      tFjet41 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f41)
      tFjet42 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f42)
      tFjet43 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f43)
      tFjet44 = Fjet(var,ndim,z,ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,z
     #p,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f44)

CCC Sum

      qqbcnum15e4 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFje
     #t44*tcoff44
      RETURN
      END

      double precision function qqbcnum16e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z-1
      t6 = 64+64*z
      t8 = z*l1
      t10 = 64*l1+64*z-64*t8
      t11 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t25 = 64-64*l1+64*t8
      taux28 = (-128*t1*l1*s1l2+t6*s2l2+t10*t11)*s1l1+(t6*s2l1+t10*t16)*
     #s1l2+(128*t1*(-1+l1)*s2l2+t25*t11)*s2l1+t25*t16*s2l2
      t1 = 1+z
      t3 = z*l1
      t4 = t3-l1-z
      t5 = -s1l2-s2l2-zp
      t11 = -s1l1-s2l1-zp
      t16 = l1-t3-1
      taux29 = (256*t1*s2l2-256*t4*t5)*s1l1+(-256*t1*s2l1+256*t4*t11)*s1
     #l2+256*t16*t5*s2l1-256*t16*t11*s2l2

CCC Coefficients

      tcoff45 = taux28
      tcoff46 = taux29

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff45)
      tFjet46 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff46)

CCC Sum

      qqbcnum16e0 = tFjet45*tcoff45+tFjet46*tcoff46
      RETURN
      END

      double precision function qqbcnum16e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 64-64*z
      t3 = -s1l2-s2l2-zp
      t4 = t2*t3
      t6 = -s1l1-s2l1-zp
      t7 = t2*t6
      t13 = z**2
      t14 = l1**2
      taux28 = t4*s1l1+t7*s1l2+t4*s2l1+t7*s2l2+128*t6*t3+64*z*(-2*t13*t1
     #4+2*t13*l1-t13-2*t14+2*l1-4*z*l1+4*z*t14-1)
      t2 = -3-3*z
      t4 = 4*l1
      t6 = 4*z*l1
      t8 = t4-1-t6+4*z
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = t6-t4+3
      taux29 = (256*t2*s2l2-256*t8*t9)*s1l1+(-256*t2*s2l1+256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff45 = taux28
      tcoff46 = taux29

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff45)
      tFjet46 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff46)

CCC Sum

      qqbcnum16e1 = tFjet45*tcoff45+tFjet46*tcoff46
      RETURN
      END

      double precision function qqbcnum16e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      taux28 = 64*z*(t1-2*z+1)
      t1 = 1+z
      t5 = 4*l1
      t7 = 4*z*l1
      t8 = -5*z+3-t5+t7
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -z-1-t7+t5
      taux29 = (512*t1*s2l2-256*t8*t9)*s1l1+(-512*t1*s2l1+256*t8*t16)*s1
     #l2+256*t21*t9*s2l1-256*t21*t16*s2l2

CCC Coefficients

      tcoff45 = taux28
      tcoff46 = taux29

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff45)
      tFjet46 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff46)

CCC Sum

      qqbcnum16e2 = tFjet45*tcoff45+tFjet46*tcoff46
      RETURN
      END

      double precision function qqbcnum16e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux28 = 0
      t1 = z-1
      t3 = t1*(-s1l2-s2l2-zp)
      t6 = t1*(-s1l1-s2l1-zp)
      taux29 = -512*t3*s1l1+512*t6*s1l2+512*t3*s2l1-512*t6*s2l2

CCC Coefficients

      tcoff45 = taux28
      tcoff46 = taux29

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff45)
      tFjet46 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff46)

CCC Sum

      qqbcnum16e3 = tFjet45*tcoff45+tFjet46*tcoff46
      RETURN
      END

      double precision function qqbcnum16e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux28 = 0
      taux29 = 0

CCC Coefficients

      tcoff45 = taux28
      tcoff46 = taux29

CCC Jet functions

      tFjet45 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff45)
      tFjet46 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff46)

CCC Sum

      qqbcnum16e4 = tFjet45*tcoff45+tFjet46*tcoff46
      RETURN
      END

      double precision function qqbcnum17e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = 1+z
      t3 = z*l1
      t4 = t3-l1-z
      t5 = -s1l2-s2l2-zp
      t11 = -s1l1-s2l1-zp
      t16 = l1-t3-1
      taux30 = (256*t1*s2l2-256*t4*t5)*s1l1+(-256*t1*s2l1+256*t4*t11)*s1
     #l2+256*t16*t5*s2l1-256*t16*t11*s2l2

CCC Coefficients

      tcoff47 = taux30

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff47)

CCC Sum

      qqbcnum17e0 = tFjet47*tcoff47
      RETURN
      END

      double precision function qqbcnum17e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = -3-3*z
      t4 = 4*l1
      t6 = 4*z*l1
      t8 = t4-1-t6+4*z
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = t6-t4+3
      taux30 = (256*t2*s2l2-256*t8*t9)*s1l1+(-256*t2*s2l1+256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff47 = taux30

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff47)

CCC Sum

      qqbcnum17e1 = tFjet47*tcoff47
      RETURN
      END

      double precision function qqbcnum17e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = 1+z
      t5 = 4*l1
      t7 = 4*z*l1
      t8 = -5*z+3-t5+t7
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -z-1-t7+t5
      taux30 = (512*t1*s2l2-256*t8*t9)*s1l1+(-512*t1*s2l1+256*t8*t16)*s1
     #l2+256*t21*t9*s2l1-256*t21*t16*s2l2

CCC Coefficients

      tcoff47 = taux30

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff47)

CCC Sum

      qqbcnum17e2 = tFjet47*tcoff47
      RETURN
      END

      double precision function qqbcnum17e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z-1
      t3 = t1*(-s1l2-s2l2-zp)
      t6 = t1*(-s1l1-s2l1-zp)
      taux30 = -512*t3*s1l1+512*t6*s1l2+512*t3*s2l1-512*t6*s2l2

CCC Coefficients

      tcoff47 = taux30

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff47)

CCC Sum

      qqbcnum17e3 = tFjet47*tcoff47
      RETURN
      END

      double precision function qqbcnum17e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux30 = 0

CCC Coefficients

      tcoff47 = taux30

CCC Jet functions

      tFjet47 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff47)

CCC Sum

      qqbcnum17e4 = tFjet47*tcoff47
      RETURN
      END

      double precision function qqbcnum18e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z-1
      t6 = 64+64*z
      t8 = z*l1
      t10 = 64*l1+64*z-64*t8
      t11 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t25 = 64-64*l1+64*t8
      taux31 = (-128*t1*l1*s1l2+t6*s2l2+t10*t11)*s1l1+(t6*s2l1+t10*t16)*
     #s1l2+(128*t1*(-1+l1)*s2l2+t25*t11)*s2l1+t25*t16*s2l2
      t1 = 1+z
      t3 = z*l1
      t4 = t3-l1-z
      t5 = -s1l2-s2l2-zp
      t11 = -s1l1-s2l1-zp
      t16 = l1-t3-1
      taux32 = (256*t1*s2l2-256*t4*t5)*s1l1+(-256*t1*s2l1+256*t4*t11)*s1
     #l2+256*t16*t5*s2l1-256*t16*t11*s2l2

CCC Coefficients

      tcoff48 = taux31
      tcoff49 = taux31
      tcoff50 = taux32
      tcoff51 = taux32

CCC Jet functions

      tFjet48 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff48)
      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)

CCC Sum

      qqbcnum18e0 = tFjet48*tcoff48+tFjet49*tcoff49+tFjet50*tcoff50+tFje
     #t51*tcoff51
      RETURN
      END

      double precision function qqbcnum18e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 64-64*z
      t3 = -s1l2-s2l2-zp
      t4 = t2*t3
      t6 = -s1l1-s2l1-zp
      t7 = t2*t6
      t13 = z**2
      t14 = l1**2
      taux31 = t4*s1l1+t7*s1l2+t4*s2l1+t7*s2l2+128*t6*t3+64*z*(-2*t13*t1
     #4+2*t13*l1-t13-2*t14+2*l1-4*z*l1+4*z*t14-1)
      t2 = -3-3*z
      t4 = 4*l1
      t6 = 4*z*l1
      t8 = t4-1-t6+4*z
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = t6-t4+3
      taux32 = (256*t2*s2l2-256*t8*t9)*s1l1+(-256*t2*s2l1+256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff48 = taux31
      tcoff49 = taux31
      tcoff50 = taux32
      tcoff51 = taux32

CCC Jet functions

      tFjet48 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff48)
      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)

CCC Sum

      qqbcnum18e1 = tFjet48*tcoff48+tFjet49*tcoff49+tFjet50*tcoff50+tFje
     #t51*tcoff51
      RETURN
      END

      double precision function qqbcnum18e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      taux31 = 64*z*(t1-2*z+1)
      t1 = 1+z
      t5 = 4*l1
      t7 = 4*z*l1
      t8 = -5*z+3-t5+t7
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -z-1-t7+t5
      taux32 = (512*t1*s2l2-256*t8*t9)*s1l1+(-512*t1*s2l1+256*t8*t16)*s1
     #l2+256*t21*t9*s2l1-256*t21*t16*s2l2

CCC Coefficients

      tcoff48 = taux31
      tcoff49 = taux31
      tcoff50 = taux32
      tcoff51 = taux32

CCC Jet functions

      tFjet48 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff48)
      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)

CCC Sum

      qqbcnum18e2 = tFjet48*tcoff48+tFjet49*tcoff49+tFjet50*tcoff50+tFje
     #t51*tcoff51
      RETURN
      END

      double precision function qqbcnum18e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux31 = 0
      t1 = z-1
      t3 = t1*(-s1l2-s2l2-zp)
      t6 = t1*(-s1l1-s2l1-zp)
      taux32 = -512*t3*s1l1+512*t6*s1l2+512*t3*s2l1-512*t6*s2l2

CCC Coefficients

      tcoff48 = taux31
      tcoff49 = taux31
      tcoff50 = taux32
      tcoff51 = taux32

CCC Jet functions

      tFjet48 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff48)
      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)

CCC Sum

      qqbcnum18e3 = tFjet48*tcoff48+tFjet49*tcoff49+tFjet50*tcoff50+tFje
     #t51*tcoff51
      RETURN
      END

      double precision function qqbcnum18e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux31 = 0
      taux32 = 0

CCC Coefficients

      tcoff48 = taux31
      tcoff49 = taux31
      tcoff50 = taux32
      tcoff51 = taux32

CCC Jet functions

      tFjet48 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff48)
      tFjet49 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff49)
      tFjet50 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff50)
      tFjet51 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff51)

CCC Sum

      qqbcnum18e4 = tFjet48*tcoff48+tFjet49*tcoff49+tFjet50*tcoff50+tFje
     #t51*tcoff51
      RETURN
      END

      double precision function qqbcnum19e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t3 = zp-1
      t8 = 64+64*zp
      t10 = zp*l1
      t12 = 64*l1+64*zp-64*t10
      t13 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t27 = 64-64*l1+64*t10
      taux33 = (-2.D0/3.D0-2.D0/3.D0*yp)*((-128*t3*l1*s1l2+t8*s2l2+t12*t
     #13)*s1l1+(t8*s2l1+t12*t18)*s1l2+(128*t3*(-1+l1)*s2l2+t27*t13)*s2l1
     #+t27*t18*s2l2)
      t3 = 1+zp
      t5 = zp*l1
      t6 = t5-l1-zp
      t7 = -s1l2-s2l2-zp
      t13 = -s1l1-s2l1-zp
      t18 = l1-t5-1
      taux34 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t3*s2l2-256*t6*t7)*s1l1+(
     #-256*t3*s2l1+256*t6*t13)*s1l2+256*t18*t7*s2l1-256*t18*t13*s2l2)

CCC Coefficients

      tcoff52 = taux33
      tcoff53 = taux33
      tcoff54 = taux34
      tcoff55 = taux34

CCC Jet functions

      tFjet52 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f55)

CCC Sum

      qqbcnum19e0 = tFjet52*tcoff52+tFjet53*tcoff53+tFjet54*tcoff54+tFje
     #t55*tcoff55
      RETURN
      END

      double precision function qqbcnum19e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t4 = 64-64*zp
      t5 = -s1l2-s2l2-zp
      t6 = t4*t5
      t8 = -s1l1-s2l1-zp
      t9 = t4*t8
      t15 = zp**2
      t16 = l1**2
      taux33 = (-2.D0/3.D0-2.D0/3.D0*yp)*(t6*s1l1+t9*s1l2+t6*s2l1+t9*s2l
     #2+128*t8*t5+64*zp*(-2*t15*t16+2*t15*l1-t15-2*t16+2*l1-4*zp*l1+4*zp
     #*t16-1))
      t4 = -3-3*zp
      t6 = 4*l1
      t8 = 4*zp*l1
      t10 = t6-1-t8+4*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = t8-t6+3
      taux34 = (-2.D0/3.D0-2.D0/3.D0*yp)*((256*t4*s2l2-256*t10*t11)*s1l1
     #+(-256*t4*s2l1+256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2
     #)

CCC Coefficients

      tcoff52 = taux33
      tcoff53 = taux33
      tcoff54 = taux34
      tcoff55 = taux34

CCC Jet functions

      tFjet52 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f55)

CCC Sum

      qqbcnum19e1 = tFjet52*tcoff52+tFjet53*tcoff53+tFjet54*tcoff54+tFje
     #t55*tcoff55
      RETURN
      END

      double precision function qqbcnum19e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t4 = zp**2
      taux33 = 64*(-2.D0/3.D0-2.D0/3.D0*yp)*zp*(t4-2*zp+1)
      t3 = 1+zp
      t7 = 4*l1
      t9 = 4*zp*l1
      t10 = -5*zp+3-t7+t9
      t11 = -s1l2-s2l2-zp
      t18 = -s1l1-s2l1-zp
      t23 = -zp-1-t9+t7
      taux34 = (-2.D0/3.D0-2.D0/3.D0*yp)*((512*t3*s2l2-256*t10*t11)*s1l1
     #+(-512*t3*s2l1+256*t10*t18)*s1l2+256*t23*t11*s2l1-256*t23*t18*s2l2
     #)

CCC Coefficients

      tcoff52 = taux33
      tcoff53 = taux33
      tcoff54 = taux34
      tcoff55 = taux34

CCC Jet functions

      tFjet52 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f55)

CCC Sum

      qqbcnum19e2 = tFjet52*tcoff52+tFjet53*tcoff53+tFjet54*tcoff54+tFje
     #t55*tcoff55
      RETURN
      END

      double precision function qqbcnum19e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux33 = 0
      t3 = zp-1
      t5 = t3*(-s1l2-s2l2-zp)
      t8 = t3*(-s1l1-s2l1-zp)
      taux34 = (-2.D0/3.D0-2.D0/3.D0*yp)*(-512*t5*s1l1+512*t8*s1l2+512*t
     #5*s2l1-512*t8*s2l2)

CCC Coefficients

      tcoff52 = taux33
      tcoff53 = taux33
      tcoff54 = taux34
      tcoff55 = taux34

CCC Jet functions

      tFjet52 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f55)

CCC Sum

      qqbcnum19e3 = tFjet52*tcoff52+tFjet53*tcoff53+tFjet54*tcoff54+tFje
     #t55*tcoff55
      RETURN
      END

      double precision function qqbcnum19e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux33 = 0
      taux34 = 0

CCC Coefficients

      tcoff52 = taux33
      tcoff53 = taux33
      tcoff54 = taux34
      tcoff55 = taux34

CCC Jet functions

      tFjet52 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f52)
      tFjet53 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f53)
      tFjet54 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f54)
      tFjet55 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f55)

CCC Sum

      qqbcnum19e4 = tFjet52*tcoff52+tFjet53*tcoff53+tFjet54*tcoff54+tFje
     #t55*tcoff55
      RETURN
      END

      double precision function qqbcnum20e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z-1
      t6 = 64+64*z
      t8 = z*l1
      t10 = 64*l1+64*z-64*t8
      t11 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t25 = 64-64*l1+64*t8
      taux35 = (-128*t1*l1*s1l2+t6*s2l2+t10*t11)*s1l1+(t6*s2l1+t10*t16)*
     #s1l2+(128*t1*(-1+l1)*s2l2+t25*t11)*s2l1+t25*t16*s2l2
      t1 = 1+z
      t3 = z*l1
      t4 = t3-l1-z
      t5 = -s1l2-s2l2-zp
      t11 = -s1l1-s2l1-zp
      t16 = l1-t3-1
      taux36 = (256*t1*s2l2-256*t4*t5)*s1l1+(-256*t1*s2l1+256*t4*t11)*s1
     #l2+256*t16*t5*s2l1-256*t16*t11*s2l2

CCC Coefficients

      tcoff56 = taux35
      tcoff57 = taux35
      tcoff58 = taux36
      tcoff59 = taux36

CCC Jet functions

      tFjet56 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff59)

CCC Sum

      qqbcnum20e0 = tFjet56*tcoff56+tFjet57*tcoff57+tFjet58*tcoff58+tFje
     #t59*tcoff59
      RETURN
      END

      double precision function qqbcnum20e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 64-64*z
      t3 = -s1l2-s2l2-zp
      t4 = t2*t3
      t6 = -s1l1-s2l1-zp
      t7 = t2*t6
      t13 = z**2
      t14 = l1**2
      taux35 = t4*s1l1+t7*s1l2+t4*s2l1+t7*s2l2+128*t6*t3+64*z*(-2*t13*t1
     #4+2*t13*l1-t13-2*t14+2*l1-4*z*l1+4*z*t14-1)
      t2 = -3-3*z
      t4 = 4*l1
      t6 = 4*z*l1
      t8 = t4-1-t6+4*z
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = t6-t4+3
      taux36 = (256*t2*s2l2-256*t8*t9)*s1l1+(-256*t2*s2l1+256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff56 = taux35
      tcoff57 = taux35
      tcoff58 = taux36
      tcoff59 = taux36

CCC Jet functions

      tFjet56 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff59)

CCC Sum

      qqbcnum20e1 = tFjet56*tcoff56+tFjet57*tcoff57+tFjet58*tcoff58+tFje
     #t59*tcoff59
      RETURN
      END

      double precision function qqbcnum20e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = z**2
      taux35 = 64*z*(t1-2*z+1)
      t1 = 1+z
      t5 = 4*l1
      t7 = 4*z*l1
      t8 = -5*z+3-t5+t7
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -z-1-t7+t5
      taux36 = (512*t1*s2l2-256*t8*t9)*s1l1+(-512*t1*s2l1+256*t8*t16)*s1
     #l2+256*t21*t9*s2l1-256*t21*t16*s2l2

CCC Coefficients

      tcoff56 = taux35
      tcoff57 = taux35
      tcoff58 = taux36
      tcoff59 = taux36

CCC Jet functions

      tFjet56 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff59)

CCC Sum

      qqbcnum20e2 = tFjet56*tcoff56+tFjet57*tcoff57+tFjet58*tcoff58+tFje
     #t59*tcoff59
      RETURN
      END

      double precision function qqbcnum20e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux35 = 0
      t1 = z-1
      t3 = t1*(-s1l2-s2l2-zp)
      t6 = t1*(-s1l1-s2l1-zp)
      taux36 = -512*t3*s1l1+512*t6*s1l2+512*t3*s2l1-512*t6*s2l2

CCC Coefficients

      tcoff56 = taux35
      tcoff57 = taux35
      tcoff58 = taux36
      tcoff59 = taux36

CCC Jet functions

      tFjet56 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff59)

CCC Sum

      qqbcnum20e3 = tFjet56*tcoff56+tFjet57*tcoff57+tFjet58*tcoff58+tFje
     #t59*tcoff59
      RETURN
      END

      double precision function qqbcnum20e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux35 = 0
      taux36 = 0

CCC Coefficients

      tcoff56 = taux35
      tcoff57 = taux35
      tcoff58 = taux36
      tcoff59 = taux36

CCC Jet functions

      tFjet56 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff56)
      tFjet57 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,1.D0,Tcoeff*tco
     #ff57)
      tFjet58 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff58)
      tFjet59 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,1.D0,Tcoeff*tco
     #ff59)

CCC Sum

      qqbcnum20e4 = tFjet56*tcoff56+tFjet57*tcoff57+tFjet58*tcoff58+tFje
     #t59*tcoff59
      RETURN
      END

      double precision function qqbcnum21e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = zp-1
      t6 = 64+64*zp
      t8 = zp*l1
      t10 = 64*l1+64*zp-64*t8
      t11 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t25 = 64-64*l1+64*t8
      taux37 = (-128*t1*l1*s1l2+t6*s2l2+t10*t11)*s1l1+(t6*s2l1+t10*t16)*
     #s1l2+(128*t1*(-1+l1)*s2l2+t25*t11)*s2l1+t25*t16*s2l2
      t1 = 1+zp
      t3 = zp*l1
      t4 = t3-l1-zp
      t5 = -s1l2-s2l2-zp
      t11 = -s1l1-s2l1-zp
      t16 = l1-t3-1
      taux38 = (256*t1*s2l2-256*t4*t5)*s1l1+(-256*t1*s2l1+256*t4*t11)*s1
     #l2+256*t16*t5*s2l1-256*t16*t11*s2l2

CCC Coefficients

      tcoff60 = taux37
      tcoff61 = taux37
      tcoff62 = taux38
      tcoff63 = taux38

CCC Jet functions

      tFjet60 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f63)

CCC Sum

      qqbcnum21e0 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet62*tcoff62+tFje
     #t63*tcoff63
      RETURN
      END

      double precision function qqbcnum21e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t2 = 64-64*zp
      t3 = -s1l2-s2l2-zp
      t4 = t2*t3
      t6 = -s1l1-s2l1-zp
      t7 = t2*t6
      t13 = zp**2
      t14 = l1**2
      taux37 = t4*s1l1+t7*s1l2+t4*s2l1+t7*s2l2+128*t6*t3+64*zp*(-2*t13*t
     #14+2*t13*l1-t13-2*t14+2*l1-4*zp*l1+4*zp*t14-1)
      t2 = -3-3*zp
      t4 = 4*l1
      t6 = 4*zp*l1
      t8 = t4-1-t6+4*zp
      t9 = -s1l2-s2l2-zp
      t15 = -s1l1-s2l1-zp
      t20 = t6-t4+3
      taux38 = (256*t2*s2l2-256*t8*t9)*s1l1+(-256*t2*s2l1+256*t8*t15)*s1
     #l2+256*t20*t9*s2l1-256*t20*t15*s2l2

CCC Coefficients

      tcoff60 = taux37
      tcoff61 = taux37
      tcoff62 = taux38
      tcoff63 = taux38

CCC Jet functions

      tFjet60 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f63)

CCC Sum

      qqbcnum21e1 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet62*tcoff62+tFje
     #t63*tcoff63
      RETURN
      END

      double precision function qqbcnum21e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      t1 = zp**2
      taux37 = 64*zp*(t1-2*zp+1)
      t1 = 1+zp
      t5 = 4*l1
      t7 = 4*zp*l1
      t8 = -5*zp+3-t5+t7
      t9 = -s1l2-s2l2-zp
      t16 = -s1l1-s2l1-zp
      t21 = -zp-1-t7+t5
      taux38 = (512*t1*s2l2-256*t8*t9)*s1l1+(-512*t1*s2l1+256*t8*t16)*s1
     #l2+256*t21*t9*s2l1-256*t21*t16*s2l2

CCC Coefficients

      tcoff60 = taux37
      tcoff61 = taux37
      tcoff62 = taux38
      tcoff63 = taux38

CCC Jet functions

      tFjet60 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f63)

CCC Sum

      qqbcnum21e2 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet62*tcoff62+tFje
     #t63*tcoff63
      RETURN
      END

      double precision function qqbcnum21e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux37 = 0
      t1 = zp-1
      t3 = t1*(-s1l2-s2l2-zp)
      t6 = t1*(-s1l1-s2l1-zp)
      taux38 = -512*t3*s1l1+512*t6*s1l2+512*t3*s2l1-512*t6*s2l2

CCC Coefficients

      tcoff60 = taux37
      tcoff61 = taux37
      tcoff62 = taux38
      tcoff63 = taux38

CCC Jet functions

      tFjet60 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f63)

CCC Sum

      qqbcnum21e3 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet62*tcoff62+tFje
     #t63*tcoff63
      RETURN
      END

      double precision function qqbcnum21e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF,Li2
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2,qflag
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      s13 = -(1.0d0-l1)*(1.0d0-zp)
      s23 = -l1*(1.0d0-zp)
      s1v = s23
      s2v = s13
      call fillspin(zp,s13,0.0d0,s23,0.0d0,0.0d0,s1v,s2v,1.0d0,zp,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)

CCC Auxiliary functions

      taux37 = 0
      taux38 = 0

CCC Coefficients

      tcoff60 = taux37
      tcoff61 = taux37
      tcoff62 = taux38
      tcoff63 = taux38

CCC Jet functions

      tFjet60 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,yp,1.D0,Tcoeff*tcof
     #f60)
      tFjet61 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,3.D0,1.D0,yp,Tcoeff*tcof
     #f61)
      tFjet62 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,yp,1.D0,Tcoeff*tcof
     #f62)
      tFjet63 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,z
     #p,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,6.D0,1.D0,yp,Tcoeff*tcof
     #f63)

CCC Sum

      qqbcnum21e4 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet62*tcoff62+tFje
     #t63*tcoff63
      RETURN
      END

