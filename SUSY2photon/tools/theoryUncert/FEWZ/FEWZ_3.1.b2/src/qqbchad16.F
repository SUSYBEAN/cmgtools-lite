CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbc16e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum16e0
      external qqbcnum16e0
      double precision qqbcnum16e1
      external qqbcnum16e1
      double precision qqbcnum16e2
      external qqbcnum16e2
      double precision qqbcnum16e3
      external qqbcnum16e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = 2*NF-33
      t8 = -1/t1
      t11 = 0.3141592653589793D1**2
      t15 = 1/l1
      t16 = -1+l1
      t18 = t15/t16
      t26 = t5*t15
      t27 = dlog(l1)
      t28 = t27**2
      tncoff1 = -t3*t5*t8/384+t5*(t11-t11*l1)*t18/18432-z3*t5/4608+t11*t
     #5*t8/9216+t26*t28/3072+t2*t5*t16*t18*t8/384+t26*t8*t27/768
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = 2*NF-33
      t8 = -1/t1
      t11 = 0.3141592653589793D1**2
      t12 = t5*t11
      t13 = -1+l1
      t14 = 1/t13
      t19 = t5*t14
      t21 = dlog(-t13)
      t22 = t21**2
      tncoff2 = -t3*t5*t8/384+t12*t14/18432-z3*t5/4608-t19*t22/3072+t12*
     #t8/9216-t2*t5*t14*t8/384-t19*t8*t21/768
      t2 = 2*NF-33
      t3 = 0.3141592653589793D1**2
      t5 = 1/l1
      t6 = -1+l1
      t7 = 1/t6
      t8 = t5*t7
      t11 = t2*t5
      t13 = dlog(-t6)
      t15 = dlog(l1)
      t19 = t13**2
      t23 = t15**2
      t27 = 1-z
      t28 = dlog(t27)
      t31 = -1/t27
      t35 = t7*t31
      tncoff3 = -t2*t3*t8/18432+t11*t7*t13*t15/1536+t11*t7*t19/3072+t11*
     #t7*t23/3072+t28*t2*t8*t31/384+t11*t35*t15/768+t11*t35*t13/768
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = 2*NF-33
      t8 = -1/t1
      t12 = dlog(z)
      t13 = t12**2
      t16 = t2*t5
      t20 = 0.3141592653589793D1**2
      t24 = 1/l1
      t28 = 1-l1
      t30 = t28*t5*t24
      t33 = -1/t28*t8
      t34 = dlog(l1)
      tncoff4 = t3*t5*t8/384+t5*t8*t13/1536-t16*t8*t12/384-t20*t5*t8/921
     #6-t16*t24*t8/384+t30*t33*t34/768-t30*t33*t12/768
      t2 = 2*NF-33
      t3 = z-1
      t4 = 1/t3
      t6 = dlog(z)
      t7 = t6**2
      t11 = dlog(-t3)
      t12 = t11**2
      t16 = t11*t2
      t17 = t4*t6
      t20 = 0.3141592653589793D1**2
      t24 = -1+l1
      t25 = 1/t24
      t29 = t2*t25
      t31 = dlog(-t24)
      tncoff5 = t2*t4*t7/1536+t12*t2*t4/384-t16*t17/384-t20*t2*t4/9216+t
     #16*t25*t4/384+t29*t4*t31/768-t29*t17/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t6 = -1+l1
      t7 = 1/t6
      t8 = 1/l1
      t11 = -1/t1
      t15 = t4*t8
      t16 = t7*t11
      t17 = dlog(l1)
      t22 = dlog(-t6)
      t26 = dlog(z)
      tncoff6 = -t2*t4*t7*t8*t11/384-t15*t16*t17/768-t15*t16*t22/768+t15
     #*t16*t26/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t10 = 0.3141592653589793D1**2
      t13 = 1/l1
      t15 = dlog(l1)
      t18 = 1-l1
      tncoff7 = t2*t4*t7/384-t10*t4/18432-t4*t13*t15/1536-t4*t13*t7/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t10 = 0.3141592653589793D1**2
      t13 = -1+l1
      t15 = t4/t13
      t17 = dlog(-t13)
      tncoff8 = t2*t4*t7/384-t10*t4/18432+t15*t17/1536+t15*t7/768
      t4 = (2*NF-33)/l1
      t5 = -1+l1
      t6 = 1/t5
      t7 = dlog(l1)
      t12 = dlog(-t5)
      tncoff9 = -t4*t6*t7/1536-t4*t6*t12/1536-t4*t6/(z-1)/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t11 = dlog(z)
      tncoff10 = -t2*t4*t7/384+t4*t7*t11/768+t4/l1*t7/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t11 = dlog(z)
      tncoff11 = -t2*t4*t7/384+t4*t7*t11/768-t4/(-1+l1)*t7/768
      tncoff12 = (2*NF-33)/l1/(-1+l1)/(z-1)/768
      t2 = 2*NF-33
      t3 = -1+l1
      tncoff13 = t2/l1/1536-t2/(z-1)/768
      t2 = 2*NF-33
      tncoff14 = -t2/(-1+l1)/1536-t2/(z-1)/768
      tncoff15 = (2*NF-33)/l1/(-1+l1)/1536
      tncoff16 = (2*NF-33)/(z-1)/768
      tncoff17 = (2*NF-33)/(z-1)/768
      tncoff18 = NF/768-11.D0/512.D0
      tncoff19 = NF/768-11.D0/512.D0

CCC Numerator functions

      tcnum1 = qqbcnum16e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum16e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum16e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum16e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum16e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum16e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum16e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum16e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum16e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum16e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum16e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum16e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum16e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum16e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum16e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum16e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum16e2(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum16e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum16e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)

CCC Sum

      qqbc16e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc16e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum16e0
      external qqbcnum16e0
      double precision qqbcnum16e1
      external qqbcnum16e1
      double precision qqbcnum16e2
      external qqbcnum16e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t10 = 0.3141592653589793D1**2
      t13 = 1/l1
      t15 = dlog(l1)
      t18 = 1-l1
      tncoff1 = t2*t4*t7/384-t10*t4/18432-t4*t13*t15/1536-t4*t13*t7/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t10 = 0.3141592653589793D1**2
      t13 = -1+l1
      t15 = t4/t13
      t17 = dlog(-t13)
      tncoff2 = t2*t4*t7/384-t10*t4/18432+t15*t17/1536+t15*t7/768
      t4 = (2*NF-33)/l1
      t5 = -1+l1
      t6 = 1/t5
      t7 = dlog(l1)
      t12 = dlog(-t5)
      tncoff3 = -t4*t6*t7/1536-t4*t6*t12/1536-t4*t6/(z-1)/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t11 = dlog(z)
      tncoff4 = -t2*t4*t7/384+t4*t7*t11/768+t4/l1*t7/768
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2*NF-33
      t7 = -1/t1
      t11 = dlog(z)
      tncoff5 = -t2*t4*t7/384+t4*t7*t11/768-t4/(-1+l1)*t7/768
      tncoff6 = (2*NF-33)/l1/(-1+l1)/(z-1)/768
      t2 = 2*NF-33
      t7 = 1-l1
      tncoff7 = t2/l1/1536-t2/(z-1)/768
      t2 = 2*NF-33
      tncoff8 = -t2/(-1+l1)/1536-t2/(z-1)/768
      tncoff9 = (2*NF-33)/l1/(-1+l1)/1536
      tncoff10 = (2*NF-33)/(z-1)/768
      tncoff11 = (2*NF-33)/(z-1)/768
      tncoff12 = NF/768-11.D0/512.D0
      tncoff13 = NF/768-11.D0/512.D0

CCC Numerator functions

      tcnum1 = qqbcnum16e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum16e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum16e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum16e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum16e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum16e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum16e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum16e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum16e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum16e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum16e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum16e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum16e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)

CCC Sum

      qqbc16e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc16e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum16e0
      external qqbcnum16e0
      double precision qqbcnum16e1
      external qqbcnum16e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t2 = 2*NF-33
      t7 = 1-l1
      tncoff1 = -t2/(z-1)/768+t2/l1/1536
      t2 = 2*NF-33
      tncoff2 = -t2/(-1+l1)/1536-t2/(z-1)/768
      tncoff3 = (2*NF-33)/l1/(-1+l1)/1536
      tncoff4 = (2*NF-33)/(z-1)/768
      tncoff5 = (2*NF-33)/(z-1)/768
      tncoff6 = NF/768-11.D0/512.D0
      tncoff7 = NF/768-11.D0/512.D0

CCC Numerator functions

      tcnum1 = qqbcnum16e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum16e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum16e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum16e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum16e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum16e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum16e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qqbc16e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc16e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum16e0
      external qqbcnum16e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = NF/768-11.D0/512.D0
      tncoff2 = NF/768-11.D0/512.D0

CCC Numerator functions

      tcnum1 = qqbcnum16e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum16e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbc16e3int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc16e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      qqbc16e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

