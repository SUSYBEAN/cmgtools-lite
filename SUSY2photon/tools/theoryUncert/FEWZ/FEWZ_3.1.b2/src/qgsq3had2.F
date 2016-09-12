CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgsq32e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgsq3num1e0
      external qgsq3num1e0
      double precision qgsq3num1e1
      external qgsq3num1e1
      double precision qgsq3num1e2
      external qgsq3num1e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1-l8
      t2 = 1/l2
      t5 = -1/t1
      t7 = 1/l8
      t8 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t9 = 1/t8
      t11 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t12 = 1/t11
      t13 = f1(1.D0,z,0.D0)
      t14 = t12*t13
      t18 = dlog(2.D0)
      t20 = 0.3141592653589793D1**2
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t27 = dlog(l1)
      t29 = t18**2
      t31 = dlog(z)
      t34 = dlog(1-z)
      t36 = f2(1.D0,z)
      t37 = dlog(t36)
      t39 = t27**2
      t41 = t34**2
      t43 = t22**2
      t45 = t25**2
      t57 = -3.D0/8192.D0*t18-t20/32768-3.D0/16384.D0*t22-3.D0/16384.D0*
     #t25-3.D0/16384.D0*t27+3.D0/8192.D0*t29+3.D0/8192.D0*t31-3.D0/4096.
     #D0*t34+3.D0/16384.D0*t37+3.D0/32768.D0*t39+3.D0/2048.D0*t41+3.D0/3
     #2768.D0*t43+3.D0/32768.D0*t45+3.D0/4096.D0*t22*t34-3.D0/8192.D0*t2
     #7*t31+3.D0/8192.D0*t22*t18-3.D0/16384.D0*t22*t37+3.D0/16384.D0*t25
     #*t27
      t90 = t37**2
      t92 = t31**2
      t94 = 3.D0/16384.D0+3.D0/16384.D0*t25*t22+3.D0/2048.D0*t34*t18-3.D
     #0/8192.D0*t22*t31+3.D0/16384.D0*t22*t27+3.D0/8192.D0*t25*t18-3.D0/
     #16384.D0*t25*t37+3.D0/4096.D0*t25*t34-3.D0/8192.D0*t25*t31-3.D0/16
     #384.D0*t27*t37+3.D0/4096.D0*t27*t34-3.D0/4096.D0*t18*t31-3.D0/4096
     #.D0*t37*t34-3.D0/8192.D0*t37*t18+3.D0/8192.D0*t37*t31-3.D0/2048.D0
     #*t34*t31+3.D0/8192.D0*t27*t18+3.D0/32768.D0*t90+3.D0/8192.D0*t92
      t98 = 1/l3
      t101 = t9*t12*t13
      t104 = dlog(l8)
      t105 = 2*t18
      t106 = 2*t31
      t107 = 4*t34
      t119 = t104+t105-t106+t22-t37+t107+l8-t25*l8-t22*l8-4*t34*l8+t37*l
     #8-2*t18*l8-t104*l8-t27*l8+2*t31*l8-1+t25+t27
      t124 = dlog(l3)
      t136 = dlog(l2)
      tncoff1 = -3.D0/8192.D0*t1*t2*t5*t7*t9*t14+(t57+t94)*t9*t14+3.D0/4
     #096.D0*t2*t98*t101+3.D0/16384.D0*t119*t5*t7*t101-3.D0/8192.D0*(t25
     #+t22+t105-t37+t107-1-t106+t27+2*t124)*t98*t101-3.D0/8192.D0*t1*t5*
     #t7*t98*t9*t14-3.D0/8192.D0*(t22+t25+t27-t37+t107+t105-1-t106+2*t13
     #6)*t2*t101
      t1 = 1/l2
      t2 = -1+l8
      t3 = 1/t2
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t7 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = 1/t7
      t10 = f1(1.D0,z,0.D0)
      t11 = t6*t8*t10
      t14 = 1/l3
      t18 = dlog(2.D0)
      t20 = 0.3141592653589793D1**2
      t23 = dlog(1-l4)
      t25 = dlog(l4)
      t29 = dlog(l1)
      t33 = f2(1.D0,z)
      t34 = dlog(t33)
      t37 = t18**2
      t41 = dlog(z)
      t45 = dlog(1-z)
      t51 = t29**2
      t55 = t45**2
      t57 = -3.D0/8192.D0*t18-t20/32768-3.D0/16384.D0*t23+3.D0/16384.D0*
     #t23*t25-3.D0/16384.D0*t25+3.D0/16384.D0*t23*t29-3.D0/16384.D0*t29-
     #3.D0/16384.D0*t25*t34+3.D0/8192.D0*t37+3.D0/8192.D0*t25*t18-3.D0/8
     #192.D0*t29*t41-3.D0/4096.D0*t45+3.D0/4096.D0*t25*t45+3.D0/8192.D0*
     #t41+3.D0/16384.D0*t34+3.D0/32768.D0*t51+3.D0/4096.D0*t29*t45+3.D0/
     #2048.D0*t55
      t62 = t23**2
      t66 = t25**2
      t90 = t34**2
      t92 = t41**2
      t94 = 3.D0/16384.D0-3.D0/16384.D0*t29*t34-3.D0/8192.D0*t23*t41+3.D
     #0/32768.D0*t62-3.D0/16384.D0*t23*t34+3.D0/32768.D0*t66+3.D0/4096.D
     #0*t23*t45+3.D0/16384.D0*t25*t29+3.D0/8192.D0*t23*t18-3.D0/8192.D0*
     #t25*t41+3.D0/2048.D0*t45*t18-3.D0/8192.D0*t34*t18+3.D0/8192.D0*t34
     #*t41-3.D0/4096.D0*t34*t45-3.D0/4096.D0*t18*t41+3.D0/8192.D0*t29*t1
     #8-3.D0/2048.D0*t45*t41+3.D0/32768.D0*t90+3.D0/8192.D0*t92
      t100 = dlog(-t2)
      t118 = 2*t18
      t119 = 4*t45
      t120 = 2*t41
      t121 = dlog(l3)
      t130 = dlog(l2)
      tncoff2 = -3.D0/8192.D0*t1*t3*t11+3.D0/4096.D0*t1*t14*t11+(t57+t94
     #)*t6*t8*t10+3.D0/16384.D0*(-l8+t100*l8+t29*l8-2*t41*l8+4*t45*l8-t3
     #4*l8+2*t18*l8+t25*l8+t23*l8)*t3/l8*t11-3.D0/8192.D0*(t23+t25+t118-
     #t34+t119-1-t120+t29+2*t121)*t14*t11-3.D0/8192.D0*t3*t14*t11-3.D0/8
     #192.D0*(t25+t23+t29-t34+t119+t118-1-t120+2*t130)*t1*t11
      t2 = -1+l8
      t3 = 1/t2
      t5 = 1/l8
      t7 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t12 = f1(1.D0,z,0.D0)
      t13 = 1/t7/t9*t12
      t16 = dlog(l4)
      t18 = dlog(1-l4)
      t19 = dlog(l8)
      t21 = dlog(-t2)
      t22 = dlog(l1)
      t24 = dlog(1-z)
      t26 = dlog(2.D0)
      t28 = dlog(z)
      t30 = f2(1.D0,z)
      t31 = dlog(t30)
      tncoff3 = 3.D0/8192.D0/l2*t3*t5*t13+3.D0/16384.D0*(1-t16-t18-t19-t
     #21-t22-4*t24-2*t26+2*t28+t31)*t3*t5*t13+3.D0/8192.D0*t3*t5/l3*t13
      t1 = -2+l3
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s23r2(z,0.D0,0.D0,t4,l4)
      t6 = 1/t5
      t8 = Dsquar321(z,l1,0.D0,t4,l4)
      t9 = 1/t8
      t12 = 1/l3
      t14 = f1(t4,z,0.D0)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t21 = f2(t4,z)
      t22 = dlog(t21)
      t23 = dlog(l1)
      t24 = dlog(z)
      t27 = dlog(1-z)
      t30 = dlog(-t1)
      t32 = dlog(l3)
      t41 = 1-l8
      tncoff4 = 3.D0/2048.D0*t2*t6*t9/l2*t12*t14-3.D0/4096.D0*(t19+t20-t
     #22+t23-1-2*t24+4*t27+2*t30+2*t32)*t2*t12*t9*t6*t14+3.D0/4096.D0*t2
     #*t12/l8*t9*t6*t14
      t1 = -2+l3
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s23r2(z,0.D0,0.D0,t4,l4)
      t6 = 1/t5
      t8 = Dsquar321(z,l1,0.D0,t4,l4)
      t9 = 1/t8
      t12 = 1/l3
      t14 = f1(t4,z,0.D0)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t21 = f2(t4,z)
      t22 = dlog(t21)
      t23 = dlog(l1)
      t24 = dlog(z)
      t27 = dlog(1-z)
      t30 = dlog(-t1)
      t32 = dlog(l3)
      t38 = t9*t6*t14
      tncoff5 = 3.D0/2048.D0*t2*t6*t9/l2*t12*t14-3.D0/4096.D0*(t19+t20-t
     #22+t23-1-2*t24+4*t27+2*t30+2*t32)*t2*t12*t38-3.D0/4096.D0*t2*t12/(
     #-1+l8)*t38
      t10 = 1-l3/2
      t11 = Dsquar321(z,l1,0.D0,t10,l4)
      t14 = s23r2(z,0.D0,0.D0,t10,l4)
      t16 = f1(t10,z,0.D0)
      tncoff6 = 3.D0/4096.D0/(-2+l3)/l3/l8/(-1+l8)/t11/t14*t16
      t1 = 1/l2
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t6 = l1*l2
      t7 = s23r2(z,t6,l2,1.D0,l4)
      t10 = f1(1.D0,z,t6)
      t11 = 1/t4/t7*t10
      t18 = dlog(l4)
      t20 = dlog(1-l4)
      t21 = f2(1.D0,z)
      t22 = dlog(t21)
      t24 = dlog(1-l2)
      t26 = dlog(1-z)
      t28 = f3(1.D0,z)
      t31 = dlog(1-t6*t28)
      t33 = dlog(1-t6)
      t34 = dlog(z)
      t36 = dlog(l1)
      t37 = dlog(2.D0)
      t39 = dlog(l2)
      t41 = t18+t20-t22+t24+4*t26+t31+t33-1-2*t34+t36+2*t37+2*t39
      tncoff7 = -3.D0/8192.D0*t1/l8*t11-3.D0/4096.D0*t1/l3*t11+3.D0/8192
     #.D0*t41*t1*t11
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = f2(1.D0,z)
      t5 = dlog(t4)
      t7 = dlog(1-l2)
      t9 = dlog(1-z)
      t11 = l1*l2
      t12 = f3(1.D0,z)
      t15 = dlog(1-t11*t12)
      t17 = dlog(1-t11)
      t18 = dlog(z)
      t20 = dlog(l1)
      t21 = dlog(2.D0)
      t23 = dlog(l2)
      t25 = t1+t3-t5+t7+4*t9+t15+t17-1-2*t18+t20+2*t21+2*t23
      t26 = 1/l2
      t28 = Dsquar321(z,l1,l2,1.D0,l4)
      t30 = s23r2(z,t11,l2,1.D0,l4)
      t33 = f1(1.D0,z,t11)
      t34 = 1/t28/t30*t33
      tncoff8 = 3.D0/8192.D0*t25*t26*t34+3.D0/8192.D0*t26/(-1+l8)*t34-3.
     #D0/4096.D0*t26/l3*t34
      t7 = Dsquar321(z,l1,l2,1.D0,l4)
      t9 = l1*l2
      t10 = s23r2(z,t9,l2,1.D0,l4)
      t13 = f1(1.D0,z,t9)
      tncoff9 = -3.D0/8192.D0/l2/(-1+l8)/l8/t7/t10*t13
      t2 = 1-l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t8 = s23r2(z,t3,l2,t2,l4)
      t11 = Dsquar321(z,l1,l2,t2,l4)
      tncoff10 = -3.D0/2048.D0*t4/(-2+l3)/t8/t11/l2/l3
      t2 = 1-l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t8 = s23r2(z,t3,l2,t2,l4)
      t11 = Dsquar321(z,l1,l2,t2,l4)
      tncoff11 = -3.D0/2048.D0*t4/(-2+l3)/t8/t11/l2/l3
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(l1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t24 = f1(1.D0,z,0.D0)
      t25 = 1/t22*t24
      tncoff12 = (3.D0/16384.D0-3.D0/16384.D0*t2-3.D0/16384.D0*t4-3.D0/1
     #6384.D0*t6+3.D0/16384.D0*t9+3.D0/8192.D0*t11-3.D0/8192.D0*t13-3.D0
     #/4096.D0*t16)*t20*t25+3.D0/16384.D0/l8*t20*t25+3.D0/8192.D0/l3*t20
     #*t25+3.D0/8192.D0/l2*t20*t25
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(l1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t24 = f1(1.D0,z,0.D0)
      t25 = 1/t22*t24
      tncoff13 = (3.D0/16384.D0-3.D0/16384.D0*t2-3.D0/16384.D0*t4-3.D0/1
     #6384.D0*t6+3.D0/16384.D0*t9+3.D0/8192.D0*t11-3.D0/8192.D0*t13-3.D0
     #/4096.D0*t16)*t20*t25-3.D0/16384.D0/(-1+l8)*t20*t25+3.D0/8192.D0/l
     #3*t20*t25+3.D0/8192.D0/l2*t20*t25
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = f1(1.D0,z,0.D0)
      tncoff14 = 3.D0/16384.D0/(-1+l8)/l8/t5/t7*t10
      t6 = 1-l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,0.D0,0.D0,t6,l4)
      t12 = f1(t6,z,0.D0)
      tncoff15 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t6 = 1-l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,0.D0,0.D0,t6,l4)
      t12 = f1(t6,z,0.D0)
      tncoff16 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,1.D0,l4)
      t8 = f1(1.D0,z,t5)
      tncoff17 = -3.D0/8192.D0/l2/t2/t6*t8
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,1.D0,l4)
      t8 = f1(1.D0,z,t5)
      tncoff18 = -3.D0/8192.D0/l2/t2/t6*t8
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff19 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff20 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgsq3num1e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgsq3num1e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgsq3num1e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgsq3num1e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgsq32e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq32e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgsq3num1e0
      external qgsq3num1e0
      double precision qgsq3num1e1
      external qgsq3num1e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(l1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t23 = 1/t22
      t24 = f1(1.D0,z,0.D0)
      t25 = t23*t24
      t27 = 1-l8
      tncoff1 = (3.D0/16384.D0-3.D0/16384.D0*t2-3.D0/16384.D0*t4-3.D0/16
     #384.D0*t6+3.D0/16384.D0*t9+3.D0/8192.D0*t11-3.D0/8192.D0*t13-3.D0/
     #4096.D0*t16)*t20*t25+3.D0/16384.D0/l8*t20*t23*t24+3.D0/8192.D0/l3*
     #t20*t25+3.D0/8192.D0/l2*t20*t25
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(l1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t24 = f1(1.D0,z,0.D0)
      t25 = 1/t22*t24
      tncoff2 = (3.D0/16384.D0-3.D0/16384.D0*t2-3.D0/16384.D0*t4-3.D0/16
     #384.D0*t6+3.D0/16384.D0*t9+3.D0/8192.D0*t11-3.D0/8192.D0*t13-3.D0/
     #4096.D0*t16)*t20*t25-3.D0/16384.D0/(-1+l8)*t20*t25+3.D0/8192.D0/l3
     #*t20*t25+3.D0/8192.D0/l2*t20*t25
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = f1(1.D0,z,0.D0)
      tncoff3 = 3.D0/16384.D0/(-1+l8)/l8/t5/t7*t10
      t6 = 1-l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,0.D0,0.D0,t6,l4)
      t12 = f1(t6,z,0.D0)
      tncoff4 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t6 = 1-l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,0.D0,0.D0,t6,l4)
      t12 = f1(t6,z,0.D0)
      tncoff5 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,1.D0,l4)
      t8 = f1(1.D0,z,t5)
      tncoff6 = -3.D0/8192.D0/l2/t2/t6*t8
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,1.D0,l4)
      t8 = f1(1.D0,z,t5)
      tncoff7 = -3.D0/8192.D0/l2/t2/t6*t8
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff8 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff9 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgsq32e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq32e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgsq3num1e0
      external qgsq3num1e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff1 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t6 = f1(1.D0,z,0.D0)
      tncoff2 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgsq32e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq32e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      qgsq32e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

