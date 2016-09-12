CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbc18e0int(var,ndim,Tcoeff)

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
      double precision qqbcnum18e0
      external qqbcnum18e0
      double precision qqbcnum18e1
      external qqbcnum18e1
      double precision qqbcnum18e2
      external qqbcnum18e2
      double precision qqbcnum18e3
      external qqbcnum18e3
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

      t1 = 0.3141592653589793D1**2
      t4 = 1/l1
      t10 = z-1
      t11 = 1/t10
      t13 = dlog(-t10)
      t14 = t13**2
      t17 = t4*t11
      t20 = dlog(l1)
      t26 = t20**2
      tncoff1 = -(-t1+t1*l1)*t4/(-1+l1)/1536-t11*t14/32+t17*t13/32+t17*t
     #20/64-z3/384+t1*t11/768+t4*t26/256
      t1 = -1+l1
      t2 = 1/t1
      t4 = dlog(-t1)
      t5 = t4**2
      t8 = z-1
      t9 = 1/t8
      t11 = dlog(-t8)
      t12 = t11**2
      t16 = 0.3141592653589793D1**2
      t21 = t2*t9
      tncoff2 = -t2*t5/256-t9*t12/32-z3/384+t16*t9/768+t16*t2/1536-t21*t
     #4/64-t21*t11/32
      t1 = 1/l1
      t2 = -1+l1
      t3 = 1/t2
      t4 = t1*t3
      t5 = dlog(l1)
      t6 = t5**2
      t10 = dlog(-t2)
      t11 = t10**2
      t17 = 0.3141592653589793D1**2
      t21 = z-1
      t22 = 1/t21
      t30 = dlog(-t21)
      tncoff3 = t4*t6/256+t4*t11/256+t4*t10*t5/128-t17*t1*t3/1536+t4*t22
     #*t5/64+t4*t22*t10/64+t4*t22*t30/32
      t1 = 1-l1
      t3 = t1/l1
      t6 = z-1
      t7 = 1/t6
      t8 = -1/t1*t7
      t9 = dlog(l1)
      t14 = dlog(-t6)
      t15 = t14**2
      t18 = dlog(z)
      t19 = t18**2
      t31 = 0.3141592653589793D1**2
      tncoff4 = t3*t8*t9/64+t7*t15/32+t7*t19/128-t3*t8*t18/64+t3*t8*t14/
     #32-t7*t14*t18/32-t31*t7/768
      t1 = z-1
      t2 = 1/t1
      t3 = dlog(z)
      t4 = t3**2
      t8 = dlog(-t1)
      t9 = t8**2
      t12 = 0.3141592653589793D1**2
      t18 = -1+l1
      t20 = 1/t18*t2
      t24 = dlog(-t18)
      tncoff5 = t2*t4/128+t2*t9/32-t12*t2/768-t2*t8*t3/32-t20*t3/64+t20*
     #t24/64+t20*t8/32
      t2 = -1+l1
      t4 = 1/l1/t2
      t5 = z-1
      t6 = 1/t5
      t7 = dlog(z)
      t12 = dlog(-t2)
      t17 = dlog(-t5)
      t21 = dlog(l1)
      tncoff6 = t4*t6*t7/64-t4*t6*t12/64-t4*t6*t17/32-t4*t6*t21/64
      t1 = 1/l1
      t2 = dlog(l1)
      t5 = z-1
      t6 = 1/t5
      t8 = dlog(-t5)
      t11 = 0.3141592653589793D1**2
      t13 = 1-l1
      tncoff7 = -t1*t2/128+t6*t8/32-t11/1536-t1*t6/64
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(-t1)
      t7 = 0.3141592653589793D1**2
      t9 = -1+l1
      t10 = 1/t9
      t14 = dlog(-t9)
      tncoff8 = t2*t4/32-t7/1536+t10*t2/64+t10*t14/128
      t2 = -1+l1
      t4 = 1/l1/t2
      t5 = dlog(l1)
      t13 = dlog(-t2)
      tncoff9 = -t4*t5/128-t4/(z-1)/64-t4*t13/128
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(-t1)
      t10 = dlog(z)
      tncoff10 = -t2*t4/32+1/l1*t2/64+t10*t2/64
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(-t1)
      t7 = dlog(z)
      tncoff11 = -t2*t4/32+t2*t7/64-1/(-1+l1)*t2/64
      tncoff12 = 1/l1/(-1+l1)/(z-1)/64
      t1 = 1-l1
      tncoff13 = 1/l1/128-1/(z-1)/64
      tncoff14 = -1/(z-1)/64-1/(-1+l1)/128
      tncoff15 = 1/(-1+l1)/l1/128
      tncoff16 = 1/(z-1)/64
      tncoff17 = 1/(z-1)/64
      tncoff18 = 1.D0/128.D0
      tncoff19 = 1.D0/128.D0

CCC Numerator functions

      tcnum1 = qqbcnum18e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum18e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum18e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum18e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum18e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum18e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum18e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum18e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum18e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum18e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum18e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum18e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum18e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum18e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum18e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum18e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum18e2(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum18e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum18e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)

CCC Sum

      qqbc18e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc18e1int(var,ndim,Tcoeff)

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
      double precision qqbcnum18e0
      external qqbcnum18e0
      double precision qqbcnum18e1
      external qqbcnum18e1
      double precision qqbcnum18e2
      external qqbcnum18e2
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

      t1 = 0.3141592653589793D1**2
      t3 = z-1
      t4 = 1/t3
      t6 = dlog(-t3)
      t9 = 1/l1
      t10 = dlog(l1)
      tncoff1 = -t1/1536+t4*t6/32-t9*t10/128-t9*t4/64
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(-t1)
      t7 = 0.3141592653589793D1**2
      t9 = -1+l1
      t10 = 1/t9
      t14 = dlog(-t9)
      tncoff2 = t2*t4/32-t7/1536+t10*t2/64+t10*t14/128
      t2 = -1+l1
      t4 = 1/l1/t2
      t5 = dlog(l1)
      t13 = dlog(-t2)
      tncoff3 = -t4*t5/128-t4/(z-1)/64-t4*t13/128
      t1 = z-1
      t2 = 1/t1
      t3 = dlog(z)
      t7 = dlog(-t1)
      t10 = 1-l1
      tncoff4 = t2*t3/64-t2*t7/32+1/l1*t2/64
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(-t1)
      t7 = dlog(z)
      tncoff5 = -t2*t4/32+t2*t7/64-1/(-1+l1)*t2/64
      tncoff6 = 1/l1/(-1+l1)/(z-1)/64
      t4 = 1-l1
      tncoff7 = -1/(z-1)/64+1/l1/128
      tncoff8 = -1/(z-1)/64-1/(-1+l1)/128
      tncoff9 = 1/(-1+l1)/l1/128
      tncoff10 = 1/(z-1)/64
      tncoff11 = 1/(z-1)/64
      tncoff12 = 1.D0/128.D0
      tncoff13 = 1.D0/128.D0

CCC Numerator functions

      tcnum1 = qqbcnum18e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum18e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum18e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum18e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum18e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum18e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum18e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum18e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum18e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum18e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum18e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum18e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum18e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)

CCC Sum

      qqbc18e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc18e2int(var,ndim,Tcoeff)

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
      double precision qqbcnum18e0
      external qqbcnum18e0
      double precision qqbcnum18e1
      external qqbcnum18e1
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

      tncoff1 = -1/(z-1)/64+1/l1/128
      tncoff2 = -1/(z-1)/64-1/(-1+l1)/128
      tncoff3 = 1/(-1+l1)/l1/128
      tncoff4 = 1/(z-1)/64
      tncoff5 = 1/(z-1)/64
      tncoff6 = 1.D0/128.D0
      tncoff7 = 1.D0/128.D0

CCC Numerator functions

      tcnum1 = qqbcnum18e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum18e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum18e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum18e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum18e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum18e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum18e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qqbc18e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc18e3int(var,ndim,Tcoeff)

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
      double precision qqbcnum18e0
      external qqbcnum18e0
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

      tncoff1 = 1.D0/128.D0
      tncoff2 = 1.D0/128.D0

CCC Numerator functions

      tcnum1 = qqbcnum18e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum18e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbc18e3int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc18e4int(var,ndim,Tcoeff)

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

      qqbc18e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

