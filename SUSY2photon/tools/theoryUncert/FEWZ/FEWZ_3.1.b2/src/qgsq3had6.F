CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgsq36e0int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
      double precision qgsq3num2e1
      external qgsq3num2e1
      double precision qgsq3num2e2
      external qgsq3num2e2
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
      t2 = 1/l8
      t5 = -1/t1
      t7 = 1/l3
      t8 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t9 = 1/t8
      t11 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t12 = 1/t11
      t13 = f1(1.D0,z,1.D0)
      t14 = t12*t13
      t18 = dlog(2.D0)
      t19 = 2*t18
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t23 = dlog(1-z)
      t24 = 4*t23
      t25 = f3(1.D0,z)
      t27 = dlog(1-t25)
      t28 = dlog(z)
      t29 = 2*t28
      t30 = dlog(l1)
      t31 = dlog(l2)
      t33 = dlog(l4)
      t35 = dlog(1-l4)
      t37 = 1/l2
      t40 = t9*t12*t13
      t71 = -3.D0/8192.D0*t33*t28-3.D0/8192.D0*t30*t28+3.D0/8192.D0*t33*
     #t18-3.D0/16384.D0*t33*t21+3.D0/4096.D0*t33*t23-3.D0/2048.D0*t23*t2
     #8+3.D0/4096.D0*t30*t23+3.D0/16384.D0*t33*t27+3.D0/8192.D0*t27*t18+
     #3.D0/4096.D0*t35*t23+3.D0/16384.D0*t35*t27
      t80 = t30**2
      t82 = t23**2
      t84 = t35**2
      t86 = t33**2
      t88 = t28**2
      t90 = t21**2
      t92 = t27**2
      t95 = -3.D0/8192.D0*t35*t28+3.D0/16384.D0*t35*t33+3.D0/8192.D0*t30
     #*t18-3.D0/4096.D0*t28*t18+3.D0/32768.D0*t80+3.D0/2048.D0*t82+3.D0/
     #32768.D0*t84+3.D0/32768.D0*t86+3.D0/8192.D0*t88+3.D0/32768.D0*t90+
     #3.D0/32768.D0*t92-3.D0/8192.D0*t18
      t104 = 0.3141592653589793D1**2
      t106 = t18**2
      t110 = 3.D0/16384.D0+3.D0/16384.D0*t21-3.D0/4096.D0*t23-3.D0/16384
     #.D0*t27+3.D0/8192.D0*t28-3.D0/16384.D0*t30-3.D0/16384.D0*t33-3.D0/
     #16384.D0*t35-t104/32768+3.D0/8192.D0*t106+3.D0/2048.D0*t23*t18
      t135 = 3.D0/16384.D0*t35*t30+3.D0/8192.D0*t21*t28-3.D0/8192.D0*t27
     #*t28+3.D0/4096.D0*t23*t27-3.D0/4096.D0*t21*t23-3.D0/8192.D0*t21*t1
     #8-3.D0/16384.D0*t21*t27-3.D0/16384.D0*t35*t21+3.D0/8192.D0*t35*t18
     #+3.D0/16384.D0*t30*t27+3.D0/16384.D0*t33*t30-3.D0/16384.D0*t30*t21
      t150 = dlog(l8)
      t156 = t24+t35+t33+t30-t29+t27-t21+t19+l8-4*t23*l8-t27*l8+2*t28*l8
     #-t30*l8-t35*l8-t150*l8-2*t18*l8+t21*l8-t33*l8+t150-1
      t161 = dlog(l3)
      tncoff1 = -3.D0/8192.D0*t1*t2*t5*t7*t9*t14-3.D0/8192.D0*(t19-t21+t
     #24+t27-1-t29+t30+2*t31+t33+t35)*t37*t40-3.D0/8192.D0*t1*t37*t5*t2*
     #t9*t14+(t71+t95+t110+t135)*t9*t14+3.D0/4096.D0*t37*t7*t40+3.D0/163
     #84.D0*t156*t5*t2*t40-3.D0/8192.D0*(t30+t27+t24-t21+t19-1-t29+2*t16
     #1+t33+t35)*t7*t40
      t1 = dlog(l1)
      t2 = f3(1.D0,z)
      t4 = dlog(1-t2)
      t6 = dlog(1-z)
      t7 = 4*t6
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t10 = dlog(2.D0)
      t11 = 2*t10
      t12 = dlog(z)
      t13 = 2*t12
      t14 = dlog(l3)
      t16 = dlog(l4)
      t18 = dlog(1-l4)
      t20 = 1/l3
      t22 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t23 = 1/t22
      t24 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t25 = 1/t24
      t27 = f1(1.D0,z,1.D0)
      t28 = t23*t25*t27
      t31 = -1+l8
      t32 = 1/t31
      t36 = dlog(l2)
      t39 = 1/l2
      t71 = 3.D0/16384.D0*t18*t16-3.D0/8192.D0*t18*t12+3.D0/16384.D0*t18
     #*t4-3.D0/4096.D0*t12*t10+3.D0/8192.D0*t4*t10-3.D0/2048.D0*t6*t12+3
     #.D0/4096.D0*t16*t6+3.D0/16384.D0*t16*t4-3.D0/16384.D0*t16*t9+3.D0/
     #4096.D0*t1*t6+3.D0/4096.D0*t18*t6
      t84 = t1**2
      t94 = t6**2
      t96 = 3.D0/8192.D0*t16*t10+3.D0/4096.D0*t6*t4-3.D0/8192.D0*t4*t12+
     #3.D0/8192.D0*t9*t12+3.D0/16384.D0*t18*t1+3.D0/2048.D0*t6*t10+3.D0/
     #32768.D0*t84+3.D0/16384.D0*t16*t1-3.D0/16384.D0*t1*t9-3.D0/8192.D0
     #*t16*t12-3.D0/8192.D0*t1*t12+3.D0/2048.D0*t94
      t110 = t18**2
      t112 = t16**2
      t114 = t9**2
      t116 = t4**2
      t118 = 3.D0/16384.D0+3.D0/16384.D0*t1*t4+3.D0/8192.D0*t18*t10-3.D0
     #/16384.D0*t18*t9-3.D0/16384.D0*t9*t4-3.D0/8192.D0*t9*t10-3.D0/4096
     #.D0*t9*t6+3.D0/32768.D0*t110+3.D0/32768.D0*t112+3.D0/32768.D0*t114
     #+3.D0/32768.D0*t116
      t119 = t12**2
      t131 = 0.3141592653589793D1**2
      t133 = t10**2
      t135 = 3.D0/8192.D0*t119-3.D0/16384.D0*t1-3.D0/16384.D0*t4-3.D0/40
     #96.D0*t6+3.D0/16384.D0*t9+3.D0/8192.D0*t12-3.D0/8192.D0*t10-3.D0/1
     #6384.D0*t16-3.D0/16384.D0*t18+3.D0/8192.D0*t1*t10-t131/32768+3.D0/
     #8192.D0*t133
      t142 = dlog(-t31)
      tncoff2 = -3.D0/8192.D0*(t1+t4+t7-t9+t11-1-t13+2*t14+t16+t18)*t20*
     #t28-3.D0/8192.D0*t32*t20*t28-3.D0/8192.D0*(t11-t9+t7+t4-1-t13+t1+2
     #*t36+t16+t18)*t39*t28-3.D0/8192.D0*t39*t32*t28+3.D0/4096.D0*t39*t2
     #0*t28+(t71+t96+t118+t135)*t23*t25*t27+3.D0/16384.D0*(-l8+t142*l8+t
     #16*l8+t18*l8+2*t10*l8-t9*l8+4*t6*l8+t4*l8-2*t12*l8+t1*l8)*t32/l8*t
     #28
      t1 = 1/l8
      t2 = -1+l8
      t3 = 1/t2
      t7 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t12 = f1(1.D0,z,1.D0)
      t13 = 1/t7/t9*t12
      t21 = dlog(z)
      t23 = dlog(l1)
      t24 = dlog(2.D0)
      t27 = dlog(1-l4)
      t28 = dlog(l4)
      t29 = dlog(l8)
      t31 = dlog(-t2)
      t32 = f2(1.D0,z)
      t33 = dlog(t32)
      t35 = dlog(1-z)
      t37 = f3(1.D0,z)
      t39 = dlog(1-t37)
      t40 = 2*t21-t23-2*t24-t27-t28-t29-t31+t33-4*t35-t39+1
      tncoff3 = 3.D0/8192.D0*t1*t3/l3*t13+3.D0/8192.D0/l2*t3*t1*t13+3.D0
     #/16384.D0*t40*t3*t1*t13
      t1 = 1-l8
      t2 = -2+l3
      t3 = 1/t2
      t5 = 1/l3
      t12 = 1-l3/2
      t13 = Dsquar311(z,l1,0.D0,t12,l4)
      t14 = 1/t13
      t16 = s23r2(z,1.D0,1.D0,t12,l4)
      t17 = 1/t16
      t18 = f1(t12,z,1.D0)
      t31 = dlog(1-z)
      t33 = f2(t12,z)
      t34 = dlog(t33)
      t35 = dlog(l1)
      t36 = dlog(z)
      t39 = dlog(-t2)
      t41 = f3(t12,z)
      t43 = dlog(1-t41)
      t44 = dlog(l3)
      t46 = dlog(l4)
      t48 = dlog(1-l4)
      tncoff4 = 3.D0/4096.D0*t3*t5/l8*t14*t17*t18+3.D0/2048.D0*t3*t17*t1
     #4/l2*t5*t18-3.D0/4096.D0*(4*t31-t34+t35-2*t36-1+2*t39+t43+2*t44+t4
     #6+t48)*t3*t5*t14*t17*t18
      t2 = dlog(1-z)
      t5 = 1-l3/2
      t6 = f2(t5,z)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = dlog(z)
      t11 = 2-l3
      t12 = dlog(t11)
      t14 = f3(t5,z)
      t16 = dlog(1-t14)
      t17 = dlog(l3)
      t19 = dlog(l4)
      t21 = dlog(1-l4)
      t24 = -1/t11
      t26 = 1/l3
      t28 = Dsquar311(z,l1,0.D0,t5,l4)
      t29 = 1/t28
      t30 = s23r2(z,1.D0,1.D0,t5,l4)
      t31 = 1/t30
      t33 = f1(t5,z,1.D0)
      t34 = t29*t31*t33
      tncoff5 = -3.D0/4096.D0*(4*t2-t7+t8-2*t9-1+2*t12+t16+2*t17+t19+t21
     #)*t24*t26*t34-3.D0/4096.D0*t24*t26/(-1+l8)*t34+3.D0/2048.D0*t24*t3
     #1*t29/l2*t26*t33
      t10 = 1-l3/2
      t11 = Dsquar311(z,l1,0.D0,t10,l4)
      t14 = s23r2(z,1.D0,1.D0,t10,l4)
      t16 = f1(t10,z,1.D0)
      tncoff6 = 3.D0/4096.D0/(-2+l3)/l3/(-1+l8)/l8/t11/t14*t16
      t1 = dlog(l1)
      t2 = f2(1.D0,z)
      t3 = dlog(t2)
      t4 = 1-l2
      t5 = dlog(t4)
      t7 = dlog(1-z)
      t9 = l1*l2
      t10 = 1-t9
      t11 = dlog(t10)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t16 = f3(1.D0,z)
      t19 = dlog(1-t16+t16*t9)
      t20 = dlog(l2)
      t22 = dlog(l4)
      t24 = dlog(1-l4)
      t25 = t1-t3+t5+4*t7+t11+2*t12-1-2*t14+t19+2*t20+t22+t24
      t26 = 1/l2
      t28 = Dsquar311(z,l1,l2,1.D0,l4)
      t29 = 1/t28
      t30 = s23r2(z,t10,t4,1.D0,l4)
      t31 = 1/t30
      t33 = f1(1.D0,z,t10)
      t34 = t29*t31*t33
      t37 = 1-l8
      tncoff7 = 3.D0/8192.D0*t25*t26*t34-3.D0/8192.D0*t26/l8*t29*t31*t33
     #-3.D0/4096.D0*t26/l3*t34
      t1 = dlog(l1)
      t2 = f2(1.D0,z)
      t3 = dlog(t2)
      t4 = 1-l2
      t5 = dlog(t4)
      t7 = dlog(1-z)
      t9 = l1*l2
      t10 = 1-t9
      t11 = dlog(t10)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t16 = f3(1.D0,z)
      t19 = dlog(1-t16+t16*t9)
      t20 = dlog(l2)
      t22 = dlog(l4)
      t24 = dlog(1-l4)
      t25 = t1-t3+t5+4*t7+t11+2*t12-1-2*t14+t19+2*t20+t22+t24
      t26 = 1/l2
      t28 = Dsquar311(z,l1,l2,1.D0,l4)
      t30 = s23r2(z,t10,t4,1.D0,l4)
      t33 = f1(1.D0,z,t10)
      t34 = 1/t28/t30*t33
      tncoff8 = 3.D0/8192.D0*t25*t26*t34+3.D0/8192.D0*t26/(-1+l8)*t34-3.
     #D0/4096.D0*t26/l3*t34
      t7 = Dsquar311(z,l1,l2,1.D0,l4)
      t10 = 1-l1*l2
      t12 = s23r2(z,t10,1-l2,1.D0,l4)
      t15 = f1(1.D0,z,t10)
      tncoff9 = -3.D0/8192.D0/l2/(-1+l8)/l8/t7/t12*t15
      t2 = 1-l3/2
      t4 = 1-l1*l2
      t5 = f1(t2,z,t4)
      t10 = s23r2(z,t4,1-l2,t2,l4)
      t13 = Dsquar311(z,l1,l2,t2,l4)
      tncoff10 = -3.D0/2048.D0*t5/(-2+l3)/t10/t13/l2/l3
      t2 = 1-l3/2
      t4 = 1-l1*l2
      t5 = f1(t2,z,t4)
      t10 = s23r2(z,t4,1-l2,t2,l4)
      t13 = Dsquar311(z,l1,l2,t2,l4)
      tncoff11 = -3.D0/2048.D0*t5/(-2+l3)/t10/t13/l2/l3
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t7 = f1(1.D0,z,1.D0)
      t8 = 1/t5*t7
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t15 = dlog(1-z)
      t17 = f3(1.D0,z)
      t19 = dlog(1-t17)
      t21 = dlog(2.D0)
      t23 = dlog(z)
      t26 = dlog(1-l4)
      t28 = dlog(l1)
      t30 = dlog(l4)
      tncoff12 = 3.D0/8192.D0/l2*t3*t8+(3.D0/16384.D0+3.D0/16384.D0*t12-
     #3.D0/4096.D0*t15-3.D0/16384.D0*t19-3.D0/8192.D0*t21+3.D0/8192.D0*t
     #23-3.D0/16384.D0*t26-3.D0/16384.D0*t28-3.D0/16384.D0*t30)*t3*t8+3.
     #D0/16384.D0/l8*t3*t8+3.D0/8192.D0/l3*t3*t8
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t7 = f1(1.D0,z,1.D0)
      t8 = 1/t5*t7
      t15 = f2(1.D0,z)
      t16 = dlog(t15)
      t19 = dlog(1-z)
      t21 = f3(1.D0,z)
      t23 = dlog(1-t21)
      t25 = dlog(2.D0)
      t27 = dlog(z)
      t30 = dlog(1-l4)
      t32 = dlog(l1)
      t34 = dlog(l4)
      tncoff13 = 3.D0/8192.D0/l3*t3*t8+3.D0/8192.D0/l2*t3*t8+(3.D0/16384
     #.D0+3.D0/16384.D0*t16-3.D0/4096.D0*t19-3.D0/16384.D0*t23-3.D0/8192
     #.D0*t25+3.D0/8192.D0*t27-3.D0/16384.D0*t30-3.D0/16384.D0*t32-3.D0/
     #16384.D0*t34)*t3*t8-3.D0/16384.D0/(-1+l8)*t3*t8
      t5 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = f1(1.D0,z,1.D0)
      tncoff14 = 3.D0/16384.D0/(-1+l8)/l8/t5/t7*t10
      t6 = 1-l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,1.D0,1.D0,t6,l4)
      t12 = f1(t6,z,1.D0)
      tncoff15 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t6 = 1-l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,1.D0,1.D0,t6,l4)
      t12 = f1(t6,z,1.D0)
      tncoff16 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,1.D0,l4)
      t10 = f1(1.D0,z,t6)
      tncoff17 = -3.D0/8192.D0/l2/t2/t8*t10
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,1.D0,l4)
      t10 = f1(1.D0,z,t6)
      tncoff18 = -3.D0/8192.D0/l2/t2/t8*t10
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff19 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff20 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l
     #7,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l
     #7,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgsq3num2e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,
     #0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgsq3num2e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,
     #1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgsq3num2e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgsq3num2e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgsq36e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq36e1int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
      double precision qgsq3num2e1
      external qgsq3num2e1
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

      t1 = f2(1.D0,z)
      t2 = dlog(t1)
      t5 = dlog(1-z)
      t7 = f3(1.D0,z)
      t9 = dlog(1-t7)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t16 = dlog(1-l4)
      t18 = dlog(l1)
      t20 = dlog(l4)
      t23 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t27 = 1/t26
      t28 = f1(1.D0,z,1.D0)
      t29 = t27*t28
      t31 = 1-l8
      tncoff1 = (3.D0/16384.D0+3.D0/16384.D0*t2-3.D0/4096.D0*t5-3.D0/163
     #84.D0*t9-3.D0/8192.D0*t11+3.D0/8192.D0*t13-3.D0/16384.D0*t16-3.D0/
     #16384.D0*t18-3.D0/16384.D0*t20)*t24*t29+3.D0/16384.D0/l8*t24*t27*t
     #28+3.D0/8192.D0/l3*t24*t29+3.D0/8192.D0/l2*t24*t29
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t7 = f1(1.D0,z,1.D0)
      t8 = 1/t5*t7
      t15 = f2(1.D0,z)
      t16 = dlog(t15)
      t19 = dlog(1-z)
      t21 = f3(1.D0,z)
      t23 = dlog(1-t21)
      t25 = dlog(2.D0)
      t27 = dlog(z)
      t30 = dlog(1-l4)
      t32 = dlog(l1)
      t34 = dlog(l4)
      tncoff2 = 3.D0/8192.D0/l3*t3*t8+3.D0/8192.D0/l2*t3*t8+(3.D0/16384.
     #D0+3.D0/16384.D0*t16-3.D0/4096.D0*t19-3.D0/16384.D0*t23-3.D0/8192.
     #D0*t25+3.D0/8192.D0*t27-3.D0/16384.D0*t30-3.D0/16384.D0*t32-3.D0/1
     #6384.D0*t34)*t3*t8-3.D0/16384.D0/(-1+l8)*t3*t8
      t5 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = f1(1.D0,z,1.D0)
      tncoff3 = 3.D0/16384.D0/(-1+l8)/l8/t5/t7*t10
      t6 = 1-l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,1.D0,1.D0,t6,l4)
      t12 = f1(t6,z,1.D0)
      tncoff4 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t6 = 1-l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t9 = s23r2(z,1.D0,1.D0,t6,l4)
      t12 = f1(t6,z,1.D0)
      tncoff5 = 3.D0/4096.D0/(-2+l3)/l3/t7/t9*t12
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,1.D0,l4)
      t10 = f1(1.D0,z,t6)
      tncoff6 = -3.D0/8192.D0/l2/t2/t8*t10
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,1.D0,l4)
      t10 = f1(1.D0,z,t6)
      tncoff7 = -3.D0/8192.D0/l2/t2/t8*t10
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff8 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff9 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgsq36e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq36e2int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
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

      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff1 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t6 = f1(1.D0,z,1.D0)
      tncoff2 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgsq36e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq36e3int(var,ndim,Tcoeff)

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

      qgsq36e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

