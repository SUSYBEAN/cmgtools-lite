CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggcnum1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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

      t1 = l2**2
      t3 = t1*z
      t5 = z**2
      t6 = t1*t5
      t9 = l2*z
      t25 = -l2+t9-z
      t26 = t25**2
      taux1 = 3.D0/2.D0*z*(2*t1-4*t3+2*t6-2*l2+1+6*t9-2*z-4*l2*t5+2*t5)*
     #(t1-2*t3+t6+t5)*(s1l1*s2l2+s1l2*s2l1)*(-1+z)/t26/t25
      t1 = l2**2
      t3 = t1*z
      t5 = z**2
      t6 = t1*t5
      t9 = l2*z
      t25 = -l2+t9-z
      t26 = t25**2
      taux2 = 6*z*(2*t1-4*t3+2*t6-2*l2+1+6*t9-2*z-4*l2*t5+2*t5)*(t1-2*t3
     #+t6+t5)*(s1l1*s2l2-s1l2*s2l1)*(-1+z)/t26/t25

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,4.D0,t3,z/t3,Tcoeff
     #*tcoff1)
      t4 = (1-z)*l2+z
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,4.D0,t4,z/t4,Tcoe
     #ff*tcoff2)
      t3 = (1-z)*l2+z
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,8.D0,t3,z/t3,Tcoeff
     #*tcoff3)
      t4 = (1-z)*l2+z
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,8.D0,t4,z/t4,Tcoe
     #ff*tcoff4)

CCC Sum

      ggcnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function ggcnum1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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

      t1 = l2**2
      t3 = t1*z
      t5 = z**2
      t6 = t1*t5
      t9 = l2*z
      t21 = -l2+t9-z
      t22 = t21**2
      taux1 = -3.D0/2.D0*z*(2*t1-4*t3+2*t6-2*l2+1+6*t9-2*z-4*l2*t5+2*t5)
     #*(t1-2*t3+t6+t5)*(-1+z)/t22/t21
      t1 = l2**2
      t3 = t1*z
      t5 = z**2
      t6 = t1*t5
      t9 = l2*z
      t25 = -l2+t9-z
      t26 = t25**2
      taux2 = -18*z*(2*t1-4*t3+2*t6-2*l2+1+6*t9-2*z-4*l2*t5+2*t5)*(t1-2*
     #t3+t6+t5)*(s1l1*s2l2-s1l2*s2l1)*(-1+z)/t26/t25

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,4.D0,t3,z/t3,Tcoeff
     #*tcoff1)
      t4 = (1-z)*l2+z
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,4.D0,t4,z/t4,Tcoe
     #ff*tcoff2)
      t3 = (1-z)*l2+z
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,8.D0,t3,z/t3,Tcoeff
     #*tcoff3)
      t4 = (1-z)*l2+z
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,8.D0,t4,z/t4,Tcoe
     #ff*tcoff4)

CCC Sum

      ggcnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function ggcnum1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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
      t1 = l2**2
      t3 = t1*z
      t5 = z**2
      t6 = t1*t5
      t9 = l2*z
      t25 = -l2+t9-z
      t26 = t25**2
      taux2 = 12*z*(2*t1-4*t3+2*t6-2*l2+1+6*t9-2*z-4*l2*t5+2*t5)*(t1-2*t
     #3+t6+t5)*(s1l1*s2l2-s1l2*s2l1)*(-1+z)/t26/t25

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t3 = (1-z)*l2+z
      tFjet1 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,4.D0,t3,z/t3,Tcoeff
     #*tcoff1)
      t4 = (1-z)*l2+z
      tFjet2 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,4.D0,t4,z/t4,Tcoe
     #ff*tcoff2)
      t3 = (1-z)*l2+z
      tFjet3 = Fjet(var,ndim,z,ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1.D
     #0,1.D0,s1l1,s2l1,0.D0,0.D0,s1l2,s2l2,0.D0,0.D0,8.D0,t3,z/t3,Tcoeff
     #*tcoff3)
      t4 = (1-z)*l2+z
      tFjet4 = Fjet(var,ndim,z,1/ux,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,1
     #.D0,1.D0,s2l1,s1l1,0.D0,0.D0,s2l2,s1l2,0.D0,0.D0,8.D0,t4,z/t4,Tcoe
     #ff*tcoff4)

CCC Sum

      ggcnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function ggcnum2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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

      t1 = yp**2
      t3 = (1-yp)**2
      t6 = -1+zp
      t10 = zp*l1
      t12 = -64*l1-64*zp+64*t10
      t17 = 64-64*l1+64*t10-128*zp
      t18 = -s1l2-s2l2-zp
      t23 = -s1l1-s2l1-zp
      t28 = -64*t6
      taux3 = (t1/4+t3/4)*((128*t6*l1*s1l2+t12*s2l2+t17*t18)*s1l1+(t12*s
     #2l1+t17*t23)*s1l2+t28*t18*s2l1+t28*t23*s2l2+128*t23*t18)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp*l1
      t7 = -zp+t6-l1
      t10 = 1-l1+t6-2*zp
      t11 = -s1l2-s2l2-zp
      t17 = -s1l1-s2l1-zp
      t22 = 1-zp
      taux4 = (t1/4+t3/4)*((256*t7*s2l2+256*t10*t11)*s1l1+(-256*t7*s2l1-
     #256*t10*t17)*s1l2+256*t22*t11*s2l1-256*t22*t17*s2l2)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tcoff
     #5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tco
     #ff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tcoff
     #7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tco
     #ff8)

CCC Sum

      ggcnum2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff
     #8
      RETURN
      END

      double precision function ggcnum2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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

      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp*l1
      t8 = -64*l1+64*t6+64
      t25 = zp**2
      t29 = l1**2
      taux3 = (t1/4+t3/4)*(t8*s2l2*s1l1+t8*s2l1*s1l2+(128*(-1+zp)*(-1+l1
     #)*s2l2+t8*(-s1l2-s2l2-zp))*s2l1+t8*(-s1l1-s2l1-zp)*s2l2-64*zp*(-2*
     #t25+2*t25*l1-1-t25*t29+2*zp*t29+2*zp-2*t6-t29))
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp*l1
      t7 = 3*t6
      t8 = 3*l1
      t9 = 3*zp
      t10 = -t7+t8+t9+1
      t13 = -3+t8-t7+6*zp
      t14 = -s1l2-s2l2-zp
      t20 = -s1l1-s2l1-zp
      t25 = t9-l1+t6-3
      taux4 = (t1/4+t3/4)*((256*t10*s2l2+256*t13*t14)*s1l1+(-256*t10*s2l
     #1-256*t13*t20)*s1l2+256*t25*t14*s2l1-256*t25*t20*s2l2)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tcoff
     #5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tco
     #ff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tcoff
     #7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tco
     #ff8)

CCC Sum

      ggcnum2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff
     #8
      RETURN
      END

      double precision function ggcnum2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision s13,s23,s1v,s2v,zp,yp,qflag
      external fillspin
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

      t1 = yp**2
      t3 = (1-yp)**2
      t7 = l1**2
      t10 = zp**2
      taux3 = -64*(t1/4+t3/4)*zp*(-2*zp*t7+t10*t7-2*l1+2*zp*l1+1+t7)
      t1 = yp**2
      t3 = (1-yp)**2
      t6 = zp*l1
      t7 = 2*zp
      t8 = -3-l1+t6-t7
      t11 = 1-l1+t6-t7
      t12 = -s1l2-s2l2-zp
      t19 = -s1l1-s2l1-zp
      t26 = 1+3*l1-3*t6-t7
      taux4 = (t1/4+t3/4)*((256*t8*s2l2+512*t11*t12)*s1l1+(-256*t8*s2l1-
     #512*t11*t19)*s1l2+256*t26*t12*s2l1-256*t26*t19*s2l2)

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = taux4

CCC Jet functions

      tFjet5 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tcoff
     #5)
      tFjet6 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,4.D0,yp,1.D0,Tcoeff*tco
     #ff6)
      tFjet7 = Fjet(var,ndim,z,ux,s13,s23,0.D0,0.D0,0.D0,s23,s13,1.D0,zp
     #,s1l1,s2l1,s3l1,0.D0,s1l2,s2l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tcoff
     #7)
      tFjet8 = Fjet(var,ndim,z,1/ux,s23,s13,0.D0,0.D0,0.D0,s13,s23,1.D0,
     #zp,s2l1,s1l1,s3l1,0.D0,s2l2,s1l2,s3l2,0.D0,8.D0,yp,1.D0,Tcoeff*tco
     #ff8)

CCC Sum

      ggcnum2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff
     #8
      RETURN
      END

