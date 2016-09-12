CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgiq74e0int(var,ndim,Tcoeff)

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
      double precision qgiq7num1e0
      external qgiq7num1e0
      double precision qgiq7num1e1
      external qgiq7num1e1
      double precision qgiq7num1e2
      external qgiq7num1e2
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

      t1 = 1/l8
      t2 = 1/l2
      t4 = Diquad721(z,l1,0.D0,1.D0,l4)
      t5 = 1/t4
      t7 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t8 = 1/t7
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = 1/t9
      t12 = f1(1.D0,z,1.D0)
      t13 = t8*t10*t12
      t16 = 0.3141592653589793D1**2
      t18 = dlog(2.D0)
      t20 = dlog(l4)
      t22 = t18**2
      t24 = f2(1.D0,z)
      t25 = dlog(t24)
      t27 = dlog(l1)
      t30 = dlog(1-l4)
      t33 = dlog(1-z)
      t35 = dlog(z)
      t37 = t20**2
      t41 = -t16/65536-3.D0/32768.D0*t18-3.D0/32768.D0*t20+3.D0/65536.D0
     #*t22+3.D0/32768.D0*t25-3.D0/32768.D0*t27-3.D0/32768.D0*t30-3.D0/81
     #92.D0*t33+3.D0/16384.D0*t35+3.D0/65536.D0*t37-3.D0/32768.D0*t20*t2
     #5
      t42 = t30**2
      t44 = t35**2
      t46 = f3(1.D0,z)
      t48 = dlog(1-t46)
      t51 = t33**2
      t59 = t25**2
      t69 = 3.D0/65536.D0*t42+3.D0/16384.D0*t44+3.D0/8192.D0*t48*t33+3.D
     #0/4096.D0*t51+3.D0/32768.D0*t27*t18+3.D0/32768.D0*t30*t20+3.D0/327
     #68.D0*t20*t27+3.D0/65536.D0*t59+3.D0/32768.D0*t30*t27-3.D0/16384.D
     #0*t27*t35-3.D0/16384.D0*t35*t18-3.D0/16384.D0*t20*t35
      t91 = 3.D0/32768.D0+3.D0/8192.D0*t27*t33+3.D0/8192.D0*t20*t33+3.D0
     #/8192.D0*t33*t18+3.D0/32768.D0*t20*t18-3.D0/4096.D0*t33*t35-3.D0/1
     #6384.D0*t30*t35+3.D0/32768.D0*t30*t18+3.D0/8192.D0*t30*t33-3.D0/32
     #768.D0*t25*t18-3.D0/8192.D0*t25*t33
      t109 = t27**2
      t113 = t48**2
      t115 = -3.D0/16384.D0*t48*t35+3.D0/32768.D0*t30*t48-3.D0/32768.D0*
     #t30*t25+3.D0/32768.D0*t48*t18+3.D0/16384.D0*t25*t35+3.D0/32768.D0*
     #t27*t48-3.D0/32768.D0*t27*t25+3.D0/32768.D0*t20*t48-3.D0/32768.D0*
     #t48+3.D0/65536.D0*t109-3.D0/32768.D0*t48*t25+3.D0/65536.D0*t113
      t120 = 1/l3
      t133 = dlog(l8)
      t141 = 4*t33
      t142 = 2*t35
      t143 = t25*l8-t48*l8-4*t33*l8-t133*l8-t20*l8-t30*l8-t18*l8+2*t35*l
     #8+t133-t27*l8-1+t27+t20+t30+l8+t18-t25+t48+t141-t142
      t153 = dlog(l2)
      t160 = dlog(l3)
      tncoff1 = 3.D0/16384.D0*t1*t2*t5*t13+(t41+t69+t91+t115)*t5*t13+3.D
     #0/16384.D0*t1*t120*t5*t13+3.D0/8192.D0*t2*t120*t5*t13+3.D0/32768.D
     #0*t143*t1/(-1+l8)*t5*t8*t10*t12-3.D0/16384.D0*(2*t153+t27+t18+t141
     #-t142-t25+t48-1+t20+t30)*t2*t5*t13-3.D0/16384.D0*(t48+t27+t141-t14
     #2-t25+t18-1+2*t160+t30+t20)*t120*t5*t13
      t1 = f3(1.D0,z)
      t3 = dlog(1-t1)
      t4 = dlog(l1)
      t6 = dlog(1-z)
      t7 = 4*t6
      t8 = dlog(z)
      t9 = 2*t8
      t10 = f2(1.D0,z)
      t11 = dlog(t10)
      t12 = dlog(2.D0)
      t13 = dlog(l3)
      t16 = dlog(1-l4)
      t17 = dlog(l4)
      t19 = 1/l3
      t21 = Diquad721(z,l1,0.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t25 = 1/t24
      t26 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t27 = 1/t26
      t29 = f1(1.D0,z,1.D0)
      t30 = t25*t27*t29
      t33 = -1+l8
      t34 = 1/t33
      t35 = 1/l2
      t45 = dlog(-t33)
      t92 = 3.D0/32768.D0*t16*t17+3.D0/32768.D0*t16*t4-3.D0/16384.D0*t8*
     #t12+3.D0/32768.D0*t17*t4+3.D0/32768.D0*t17*t12-3.D0/4096.D0*t6*t8-
     #3.D0/16384.D0*t16*t8+3.D0/32768.D0*t16*t12+3.D0/8192.D0*t6*t12+3.D
     #0/8192.D0*t16*t6-3.D0/16384.D0*t17*t8
      t109 = -3.D0/16384.D0*t4*t8+3.D0/32768.D0*t4*t12+3.D0/8192.D0*t17*
     #t6-3.D0/32768.D0*t4-3.D0/32768.D0*t3-3.D0/32768.D0*t12+3.D0/32768.
     #D0*t11+3.D0/16384.D0*t8-3.D0/8192.D0*t6-3.D0/32768.D0*t17-3.D0/327
     #68.D0*t16-3.D0/32768.D0*t17*t11
      t111 = 0.3141592653589793D1**2
      t123 = t12**2
      t131 = 3.D0/32768.D0-t111/65536-3.D0/8192.D0*t11*t6-3.D0/32768.D0*
     #t11*t12-3.D0/32768.D0*t3*t11+3.D0/8192.D0*t3*t6-3.D0/32768.D0*t16*
     #t11+3.D0/65536.D0*t123+3.D0/32768.D0*t16*t3-3.D0/16384.D0*t3*t8+3.
     #D0/32768.D0*t3*t12
      t140 = t16**2
      t142 = t8**2
      t144 = t6**2
      t146 = t17**2
      t148 = t4**2
      t150 = t11**2
      t152 = t3**2
      t156 = 3.D0/32768.D0*t17*t3-3.D0/32768.D0*t4*t11+3.D0/32768.D0*t4*
     #t3+3.D0/16384.D0*t11*t8+3.D0/65536.D0*t140+3.D0/16384.D0*t142+3.D0
     #/4096.D0*t144+3.D0/65536.D0*t146+3.D0/65536.D0*t148+3.D0/65536.D0*
     #t150+3.D0/65536.D0*t152+3.D0/8192.D0*t4*t6
      t161 = dlog(l2)
      tncoff2 = -3.D0/16384.D0*(t3+t4+t7-t9-t11+t12-1+2*t13+t16+t17)*t19
     #*t22*t30-3.D0/16384.D0*t34*t35*t22*t30-3.D0/16384.D0*t34*t19*t22*t
     #30+3.D0/32768.D0*(t45*l8+t4*l8-t11*l8+t3*l8+4*t6*l8+t17*l8+t16*l8+
     #t12*l8-2*t8*l8-l8)/l8*t34*t22*t25*t27*t29+3.D0/8192.D0*t35*t19*t22
     #*t30+(t92+t109+t131+t156)*t22*t30-3.D0/16384.D0*(2*t161+t4+t12+t7-
     #t9-t11+t3-1+t17+t16)*t35*t22*t30
      t1 = 1/l8
      t2 = -1+l8
      t3 = 1/t2
      t4 = t1*t3
      t7 = Diquad721(z,l1,0.D0,1.D0,l4)
      t9 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t12 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t14 = f1(1.D0,z,1.D0)
      t16 = 1/t7/t9/t12*t14
      t23 = dlog(2.D0)
      t24 = f2(1.D0,z)
      t25 = dlog(t24)
      t26 = f3(1.D0,z)
      t28 = dlog(1-t26)
      t30 = dlog(1-z)
      t32 = dlog(z)
      t34 = dlog(l1)
      t35 = dlog(l4)
      t37 = dlog(1-l4)
      t39 = dlog(-t2)
      t40 = dlog(l8)
      t41 = -t23+t25-t28-4*t30+2*t32-t34-t35-t37-t39-t40+1
      tncoff3 = 3.D0/16384.D0*t4/l2*t16+3.D0/16384.D0*t4/l3*t16+3.D0/327
     #68.D0*t41*t1*t3*t16
      t1 = -2+l3
      t2 = 1/t1
      t3 = 1/l3
      t8 = 1-l3/2
      t9 = Diquad721(z,l1,0.D0,t8,l4)
      t11 = s34r2(z,1.D0,1.D0,t8,l4)
      t14 = s23r2(z,1.D0,1.D0,t8,l4)
      t16 = f1(t8,z,1.D0)
      t18 = 1/t9/t11/t14*t16
      t27 = dlog(1-z)
      t29 = dlog(l1)
      t30 = dlog(2.D0)
      t31 = dlog(z)
      t33 = f2(t8,z)
      t34 = dlog(t33)
      t35 = f3(t8,z)
      t37 = dlog(1-t35)
      t39 = dlog(-t1)
      t41 = dlog(l3)
      t44 = dlog(1-l4)
      t45 = dlog(l4)
      t46 = 4*t27+t29-t30-2*t31-t34+t37+2*t39-1+2*t41+t44+t45
      tncoff4 = 3.D0/8192.D0*t2*t3/l8*t18+3.D0/4096.D0*t3/l2*t2*t18-3.D0
     #/8192.D0*t46*t3*t2*t18
      t2 = dlog(1-z)
      t4 = dlog(l1)
      t5 = dlog(2.D0)
      t6 = dlog(z)
      t9 = 1-l3/2
      t10 = f2(t9,z)
      t11 = dlog(t10)
      t12 = f3(t9,z)
      t14 = dlog(1-t12)
      t15 = 2-l3
      t16 = dlog(t15)
      t18 = dlog(l3)
      t21 = dlog(1-l4)
      t22 = dlog(l4)
      t23 = 4*t2+t4-t5-2*t6-t11+t14+2*t16-1+2*t18+t21+t22
      t24 = 1/l3
      t27 = -1/t15
      t29 = Diquad721(z,l1,0.D0,t9,l4)
      t31 = s34r2(z,1.D0,1.D0,t9,l4)
      t34 = s23r2(z,1.D0,1.D0,t9,l4)
      t36 = f1(t9,z,1.D0)
      t38 = 1/t29/t31/t34*t36
      tncoff5 = -3.D0/8192.D0*t23*t24*t27*t38-3.D0/8192.D0*t27*t24/(-1+l
     #8)*t38+3.D0/4096.D0*t24/l2*t27*t38
      t11 = 1-l3/2
      t12 = Diquad721(z,l1,0.D0,t11,l4)
      t14 = s34r2(z,1.D0,1.D0,t11,l4)
      t17 = s23r2(z,1.D0,1.D0,t11,l4)
      t19 = f1(t11,z,1.D0)
      tncoff6 = 3.D0/8192.D0/(-2+l3)/l3/l8/(-1+l8)/t12/t14/t17*t19
      t1 = l1*l2
      t2 = -2+t1
      t3 = 1/t2
      t4 = 1/l2
      t8 = Diquad721(z,l1,l2,1.D0,l4)
      t11 = 1-t1/2
      t12 = 1-l2
      t13 = s34r2(z,t11,t12,1.D0,l4)
      t16 = s23r2(z,t11,t12,1.D0,l4)
      t18 = f1(1.D0,z,t11)
      t20 = 1/t8/t13/t16*t18
      t28 = dlog(l2)
      t30 = f3(1.D0,z)
      t35 = dlog(2-2*t30+t30*l1*l2)
      t36 = dlog(l1)
      t38 = dlog(1-z)
      t40 = dlog(z)
      t42 = dlog(2.D0)
      t44 = dlog(-t2)
      t45 = dlog(t12)
      t46 = f2(1.D0,z)
      t47 = dlog(t46)
      t48 = dlog(l4)
      t50 = dlog(1-l4)
      t51 = 2*t28+t35+t36+4*t38-2*t40-t42+t44+t45-t47-1+t48+t50
      tncoff7 = 3.D0/8192.D0*t3*t4/l8*t20+3.D0/4096.D0/l3*t3*t4*t20-3.D0
     #/8192.D0*t51*t3*t4*t20
      t1 = dlog(l2)
      t3 = f3(1.D0,z)
      t8 = dlog(2-2*t3+t3*l1*l2)
      t9 = dlog(l1)
      t11 = dlog(1-z)
      t13 = dlog(z)
      t15 = dlog(2.D0)
      t16 = l1*l2
      t17 = 2-t16
      t18 = dlog(t17)
      t19 = 1-l2
      t20 = dlog(t19)
      t21 = f2(1.D0,z)
      t22 = dlog(t21)
      t23 = dlog(l4)
      t25 = dlog(1-l4)
      t26 = 2*t1+t8+t9+4*t11-2*t13-t15+t18+t20-t22-1+t23+t25
      t28 = -1/t17
      t30 = 1/l2
      t32 = Diquad721(z,l1,l2,1.D0,l4)
      t35 = 1-t16/2
      t36 = s34r2(z,t35,t19,1.D0,l4)
      t39 = s23r2(z,t35,t19,1.D0,l4)
      t41 = f1(1.D0,z,t35)
      t43 = 1/t32/t36/t39*t41
      tncoff8 = -3.D0/8192.D0*t26*t28*t30*t43-3.D0/8192.D0*t28*t30/(-1+l
     #8)*t43+3.D0/4096.D0/l3*t28*t30*t43
      t1 = l1*l2
      t11 = Diquad721(z,l1,l2,1.D0,l4)
      t14 = 1-t1/2
      t15 = 1-l2
      t16 = s34r2(z,t14,t15,1.D0,l4)
      t19 = s23r2(z,t14,t15,1.D0,l4)
      t21 = f1(1.D0,z,t14)
      tncoff9 = 3.D0/8192.D0/(-2+t1)/l2/(-1+l8)/l8/t11/t16/t19*t21
      t2 = 1-l3/2
      t3 = l1*l2
      t5 = 1-t3/2
      t6 = f1(t2,z,t5)
      t7 = 1-l2
      t8 = s23r2(z,t5,t7,t2,l4)
      t11 = Diquad721(z,l1,l2,t2,l4)
      t19 = s34r2(z,t5,t7,t2,l4)
      tncoff10 = 3.D0/2048.D0*t6/t8/t11/(-2+t3)/(-2+l3)/t19/l2/l3
      t2 = 1-l3/2
      t3 = l1*l2
      t5 = 1-t3/2
      t6 = f1(t2,z,t5)
      t7 = 1-l2
      t8 = s23r2(z,t5,t7,t2,l4)
      t11 = Diquad721(z,l1,l2,t2,l4)
      t19 = s34r2(z,t5,t7,t2,l4)
      tncoff11 = 3.D0/2048.D0*t6/t8/t11/(-2+t3)/(-2+l3)/t19/l2/l3
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(1.D0,z)
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(2.D0)
      t13 = f3(1.D0,z)
      t15 = dlog(1-t13)
      t18 = dlog(1-z)
      t20 = dlog(l1)
      t23 = Diquad721(z,l1,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t28 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t31 = f1(1.D0,z,1.D0)
      t32 = 1/t26/t28*t31
      tncoff12 = (-3.D0/32768.D0*t2-3.D0/32768.D0*t4+3.D0/32768.D0*t7+3.
     #D0/16384.D0*t9-3.D0/32768.D0*t11-3.D0/32768.D0*t15-3.D0/8192.D0*t1
     #8-3.D0/32768.D0*t20+3.D0/32768.D0)*t24*t32+3.D0/32768.D0/l8*t24*t3
     #2+3.D0/16384.D0/l2*t24*t32+3.D0/16384.D0/l3*t24*t32
      t2 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = f1(1.D0,z,1.D0)
      t11 = 1/t5/t7*t10
      t20 = dlog(1-l4)
      t22 = dlog(l4)
      t24 = f2(1.D0,z)
      t25 = dlog(t24)
      t27 = dlog(z)
      t29 = dlog(2.D0)
      t31 = f3(1.D0,z)
      t33 = dlog(1-t31)
      t36 = dlog(1-z)
      t38 = dlog(l1)
      tncoff13 = 3.D0/16384.D0/l3*t3*t11-3.D0/32768.D0/(-1+l8)*t3*t11+(-
     #3.D0/32768.D0*t20-3.D0/32768.D0*t22+3.D0/32768.D0*t25+3.D0/16384.D
     #0*t27-3.D0/32768.D0*t29-3.D0/32768.D0*t33-3.D0/8192.D0*t36-3.D0/32
     #768.D0*t38+3.D0/32768.D0)*t3*t11+3.D0/16384.D0/l2*t3*t11
      t5 = Diquad721(z,l1,0.D0,1.D0,l4)
      t8 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff14 = 3.D0/32768.D0/l8/(-1+l8)/t5/t8/t10*t13
      t6 = 1-l3/2
      t7 = Diquad721(z,l1,0.D0,t6,l4)
      t10 = s34r2(z,1.D0,1.D0,t6,l4)
      t12 = s23r2(z,1.D0,1.D0,t6,l4)
      t15 = f1(t6,z,1.D0)
      tncoff15 = 3.D0/8192.D0/l3/(-2+l3)/t7/t10/t12*t15
      t6 = 1-l3/2
      t7 = Diquad721(z,l1,0.D0,t6,l4)
      t10 = s34r2(z,1.D0,1.D0,t6,l4)
      t12 = s23r2(z,1.D0,1.D0,t6,l4)
      t15 = f1(t6,z,1.D0)
      tncoff16 = 3.D0/8192.D0/l3/(-2+l3)/t7/t10/t12*t15
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,1.D0,l4)
      t14 = s23r2(z,t10,t11,1.D0,l4)
      t17 = f1(1.D0,z,t10)
      tncoff17 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,1.D0,l4)
      t14 = s23r2(z,t10,t11,1.D0,l4)
      t17 = f1(1.D0,z,t10)
      tncoff18 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff19 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff20 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5,l6
     #,l7,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5,l6
     #,l7,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l
     #7,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l
     #7,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgiq7num1e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgiq7num1e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgiq74e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq74e1int(var,ndim,Tcoeff)

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
      double precision qgiq7num1e0
      external qgiq7num1e0
      double precision qgiq7num1e1
      external qgiq7num1e1
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
      t6 = f2(1.D0,z)
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(2.D0)
      t13 = f3(1.D0,z)
      t15 = dlog(1-t13)
      t18 = dlog(1-z)
      t20 = dlog(l1)
      t23 = Diquad721(z,l1,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t28 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t31 = f1(1.D0,z,1.D0)
      t32 = 1/t26/t28*t31
      tncoff1 = (-3.D0/32768.D0*t2-3.D0/32768.D0*t4+3.D0/32768.D0*t7+3.D
     #0/16384.D0*t9-3.D0/32768.D0*t11-3.D0/32768.D0*t15-3.D0/8192.D0*t18
     #-3.D0/32768.D0*t20+3.D0/32768.D0)*t24*t32+3.D0/32768.D0/l8*t24*t32
     #+3.D0/16384.D0/l2*t24*t32+3.D0/16384.D0/l3*t24*t32
      t2 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = f1(1.D0,z,1.D0)
      t11 = 1/t5/t7*t10
      t20 = dlog(1-l4)
      t22 = dlog(l4)
      t24 = f2(1.D0,z)
      t25 = dlog(t24)
      t27 = dlog(z)
      t29 = dlog(2.D0)
      t31 = f3(1.D0,z)
      t33 = dlog(1-t31)
      t36 = dlog(1-z)
      t38 = dlog(l1)
      tncoff2 = 3.D0/16384.D0/l3*t3*t11-3.D0/32768.D0/(-1+l8)*t3*t11+(-3
     #.D0/32768.D0*t20-3.D0/32768.D0*t22+3.D0/32768.D0*t25+3.D0/16384.D0
     #*t27-3.D0/32768.D0*t29-3.D0/32768.D0*t33-3.D0/8192.D0*t36-3.D0/327
     #68.D0*t38+3.D0/32768.D0)*t3*t11+3.D0/16384.D0/l2*t3*t11
      t5 = Diquad721(z,l1,0.D0,1.D0,l4)
      t8 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff3 = 3.D0/32768.D0/l8/(-1+l8)/t5/t8/t10*t13
      t6 = 1-l3/2
      t7 = Diquad721(z,l1,0.D0,t6,l4)
      t10 = s34r2(z,1.D0,1.D0,t6,l4)
      t12 = s23r2(z,1.D0,1.D0,t6,l4)
      t15 = f1(t6,z,1.D0)
      tncoff4 = 3.D0/8192.D0/l3/(-2+l3)/t7/t10/t12*t15
      t6 = 1-l3/2
      t7 = Diquad721(z,l1,0.D0,t6,l4)
      t10 = s34r2(z,1.D0,1.D0,t6,l4)
      t12 = s23r2(z,1.D0,1.D0,t6,l4)
      t15 = f1(t6,z,1.D0)
      tncoff5 = 3.D0/8192.D0/l3/(-2+l3)/t7/t10/t12*t15
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,1.D0,l4)
      t14 = s23r2(z,t10,t11,1.D0,l4)
      t17 = f1(1.D0,z,t10)
      tncoff6 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,1.D0,l4)
      t14 = s23r2(z,t10,t11,1.D0,l4)
      t17 = f1(1.D0,z,t10)
      tncoff7 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff8 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff9 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgiq74e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq74e2int(var,ndim,Tcoeff)

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
      double precision qgiq7num1e0
      external qgiq7num1e0
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

      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff1 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,1.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,z,1.D0)
      tncoff2 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgiq74e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq74e3int(var,ndim,Tcoeff)

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

      qgiq74e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

