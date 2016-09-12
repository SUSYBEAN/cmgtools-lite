CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgiq73e0int(var,ndim,Tcoeff)

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

      t2 = dlog(1-l4)
      t4 = dlog(z)
      t7 = f3(0.D0,z)
      t9 = dlog(1-t7)
      t11 = dlog(2.D0)
      t14 = dlog(1-z)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t20 = dlog(l1)
      t22 = dlog(l8)
      t24 = dlog(l4)
      t26 = 4*t14
      t27 = 2*t4
      t28 = -t2*l8+2*t4*l8-t9*l8-t11*l8-4*t14*l8+t18*l8-t20*l8-t22*l8-t2
     #4*l8+t9+t11+t26+t2+t22+l8+t24+t20-t18-t27-1
      t29 = 1/l8
      t34 = Diquad721(z,l1,0.D0,0.D0,l4)
      t35 = 1/t34
      t36 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t37 = 1/t36
      t39 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t40 = 1/t39
      t41 = f1(0.D0,z,1.D0)
      t46 = dlog(l3)
      t49 = 1/l3
      t53 = t37*t40*t41
      t56 = t18**2
      t58 = t2**2
      t60 = t20**2
      t62 = t9**2
      t70 = 3.D0/32768.D0+3.D0/65536.D0*t56+3.D0/65536.D0*t58+3.D0/65536
     #.D0*t60+3.D0/65536.D0*t62+3.D0/16384.D0*t4-3.D0/32768.D0*t2-3.D0/3
     #2768.D0*t11-3.D0/32768.D0*t9-3.D0/8192.D0*t14+3.D0/32768.D0*t18
      t77 = t11**2
      t79 = 0.3141592653589793D1**2
      t93 = -3.D0/32768.D0*t20-3.D0/32768.D0*t24-3.D0/32768.D0*t9*t18+3.
     #D0/32768.D0*t24*t20+3.D0/65536.D0*t77-t79/65536+3.D0/32768.D0*t2*t
     #24+3.D0/32768.D0*t20*t11+3.D0/32768.D0*t20*t9+3.D0/8192.D0*t24*t14
     #-3.D0/8192.D0*t18*t14-3.D0/32768.D0*t18*t11
      t117 = 3.D0/16384.D0*t18*t4-3.D0/32768.D0*t20*t18+3.D0/8192.D0*t20
     #*t14-3.D0/16384.D0*t24*t4-3.D0/32768.D0*t24*t18+3.D0/8192.D0*t9*t1
     #4-3.D0/16384.D0*t4*t11-3.D0/16384.D0*t20*t4+3.D0/32768.D0*t2*t20+3
     #.D0/8192.D0*t2*t14+3.D0/32768.D0*t2*t11
      t136 = t24**2
      t138 = t14**2
      t140 = t4**2
      t142 = 3.D0/32768.D0*t2*t9-3.D0/16384.D0*t2*t4-3.D0/16384.D0*t9*t4
     #+3.D0/32768.D0*t9*t11-3.D0/4096.D0*t14*t4-3.D0/32768.D0*t2*t18+3.D
     #0/32768.D0*t24*t11+3.D0/32768.D0*t24*t9+3.D0/8192.D0*t14*t11+3.D0/
     #65536.D0*t136+3.D0/4096.D0*t138+3.D0/16384.D0*t140
      t147 = 1/l2
      t160 = dlog(l2)
      tncoff1 = 3.D0/32768.D0*t28*t29/(-1+l8)*t35*t37*t40*t41-3.D0/16384
     #.D0*(2*t46+t20-t27-1-t18+t9+t26+t11+t2+t24)*t49*t35*t53+(t70+t93+t
     #117+t142)*t35*t53+3.D0/16384.D0*t29*t147*t35*t53+3.D0/8192.D0*t147
     #*t49*t35*t53+3.D0/16384.D0*t29*t49*t35*t53-3.D0/16384.D0*(t2-1+t9+
     #t11+t26-t27-t18+t20+2*t160+t24)*t147*t35*t53
      t1 = dlog(z)
      t5 = dlog(1-z)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t12 = dlog(1-t10)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(l1)
      t19 = dlog(l4)
      t21 = 1-l8
      t22 = dlog(t21)
      t25 = dlog(1-l4)
      t31 = -1/t21
      t33 = Diquad721(z,l1,0.D0,0.D0,l4)
      t34 = 1/t33
      t35 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t36 = 1/t35
      t38 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t39 = 1/t38
      t40 = f1(0.D0,z,1.D0)
      t45 = 1/l2
      t49 = t36*t39*t40
      t52 = dlog(l3)
      t54 = 2*t1
      t55 = 4*t5
      t57 = 1/l3
      t92 = 3.D0/16384.D0*t15*t1-3.D0/32768.D0*t17*t15+3.D0/8192.D0*t17*
     #t5-3.D0/16384.D0*t19*t1-3.D0/32768.D0*t19*t15+3.D0/8192.D0*t12*t5-
     #3.D0/16384.D0*t1*t8-3.D0/16384.D0*t17*t1+3.D0/32768.D0*t25*t17+3.D
     #0/8192.D0*t25*t5+3.D0/32768.D0*t25*t8
      t115 = t19**2
      t117 = 3.D0/32768.D0*t25*t12-3.D0/16384.D0*t25*t1-3.D0/16384.D0*t1
     #2*t1+3.D0/8192.D0*t19*t5+3.D0/32768.D0*t17*t12+3.D0/32768.D0*t12*t
     #8-3.D0/4096.D0*t5*t1-3.D0/32768.D0*t25*t15+3.D0/32768.D0*t19*t8+3.
     #D0/32768.D0*t19*t12+3.D0/8192.D0*t5*t8+3.D0/65536.D0*t115
      t119 = t5**2
      t121 = t1**2
      t123 = t15**2
      t125 = t25**2
      t129 = t17**2
      t131 = t12**2
      t138 = 3.D0/32768.D0+3.D0/4096.D0*t119+3.D0/16384.D0*t121+3.D0/655
     #36.D0*t123+3.D0/65536.D0*t125+3.D0/32768.D0*t17*t8+3.D0/65536.D0*t
     #129+3.D0/65536.D0*t131+3.D0/32768.D0*t25*t19+3.D0/32768.D0*t19*t17
     #+3.D0/16384.D0*t1
      t139 = 0.3141592653589793D1**2
      t144 = t8**2
      t156 = -t139/65536-3.D0/8192.D0*t5-3.D0/32768.D0*t8-3.D0/32768.D0*
     #t12+3.D0/65536.D0*t144-3.D0/32768.D0*t17-3.D0/32768.D0*t19+3.D0/32
     #768.D0*t15-3.D0/32768.D0*t12*t15-3.D0/32768.D0*t25-3.D0/32768.D0*t
     #15*t8-3.D0/8192.D0*t15*t5
      t161 = dlog(l2)
      tncoff2 = 3.D0/32768.D0*(-2*t1*l8+4*t5*l8+t8*l8+t12*l8-t15*l8+t17*
     #l8+t19*l8+t22*l8+t25*l8-l8)/l8*t31*t34*t36*t39*t40-3.D0/16384.D0*t
     #31*t45*t34*t49-3.D0/16384.D0*(2*t52+t17-t54-1-t15+t12+t55+t8+t25+t
     #19)*t57*t34*t49+3.D0/8192.D0*t45*t57*t34*t49-3.D0/16384.D0*t31*t57
     #*t34*t49+(t92+t117+t138+t156)*t34*t49-3.D0/16384.D0*(t25-1+t12+t8+
     #t55-t54-t15+t17+2*t161+t19)*t45*t34*t49
      t1 = dlog(l8)
      t2 = 1-l8
      t3 = dlog(t2)
      t4 = dlog(l4)
      t5 = dlog(z)
      t7 = f3(0.D0,z)
      t9 = dlog(1-t7)
      t11 = dlog(1-z)
      t13 = dlog(2.D0)
      t14 = dlog(l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(1-l4)
      t19 = 1-t1-t3-t4+2*t5-t9-4*t11-t13-t14+t16-t18
      t20 = 1/l8
      t23 = -1/t2
      t25 = Diquad721(z,l1,0.D0,0.D0,l4)
      t27 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t30 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t32 = f1(0.D0,z,1.D0)
      t34 = 1/t25/t27/t30*t32
      t37 = t23*t20
      tncoff3 = 3.D0/32768.D0*t19*t20*t23*t34+3.D0/16384.D0*t37/l2*t34+3
     #.D0/16384.D0*t37/l3*t34
      t1 = dlog(l3)
      t3 = dlog(2.D0)
      t4 = l3/2
      t5 = f2(t4,z)
      t6 = dlog(t5)
      t7 = dlog(z)
      t9 = 2-l3
      t10 = dlog(t9)
      t12 = f3(t4,z)
      t14 = dlog(1-t12)
      t16 = dlog(1-z)
      t18 = dlog(l1)
      t20 = dlog(1-l4)
      t21 = dlog(l4)
      t22 = 2*t1-t3-t6-2*t7+2*t10+t14+4*t16-1+t18+t20+t21
      t24 = -1/t9
      t26 = 1/l3
      t28 = Diquad721(z,l1,0.D0,t4,l4)
      t30 = s34r2(z,1.D0,1.D0,t4,l4)
      t33 = s23r2(z,1.D0,1.D0,t4,l4)
      t35 = f1(t4,z,1.D0)
      t37 = 1/t28/t30/t33*t35
      t40 = t26*t24
      tncoff4 = -3.D0/8192.D0*t22*t24*t26*t37+3.D0/4096.D0*t40/l2*t37+3.
     #D0/8192.D0*t40/l8*t37
      t1 = 1/l3
      t2 = -2+l3
      t3 = 1/t2
      t4 = t1*t3
      t7 = l3/2
      t8 = Diquad721(z,l1,0.D0,t7,l4)
      t10 = s34r2(z,1.D0,1.D0,t7,l4)
      t13 = s23r2(z,1.D0,1.D0,t7,l4)
      t15 = f1(t7,z,1.D0)
      t17 = 1/t8/t10/t13*t15
      t20 = dlog(l3)
      t22 = dlog(2.D0)
      t23 = f2(t7,z)
      t24 = dlog(t23)
      t25 = dlog(z)
      t28 = dlog(-t2)
      t30 = f3(t7,z)
      t32 = dlog(1-t30)
      t34 = dlog(1-z)
      t36 = dlog(l1)
      t38 = dlog(1-l4)
      t39 = dlog(l4)
      t40 = 2*t20-t22-t24-2*t25+2*t28+t32+4*t34-1+t36+t38+t39
      tncoff5 = 3.D0/4096.D0*t4/l2*t17-3.D0/8192.D0*t40*t3*t1*t17-3.D0/8
     #192.D0*t4/(-1+l8)*t17
      t10 = l3/2
      t11 = Diquad721(z,l1,0.D0,t10,l4)
      t13 = s34r2(z,1.D0,1.D0,t10,l4)
      t16 = s23r2(z,1.D0,1.D0,t10,l4)
      t18 = f1(t10,z,1.D0)
      tncoff6 = 3.D0/8192.D0/(-2+l3)/l3/(-1+l8)/l8/t11/t13/t16*t18
      t1 = l1*l2
      t2 = -2+t1
      t3 = 1/t2
      t4 = 1/l2
      t8 = Diquad721(z,l1,l2,0.D0,l4)
      t11 = 1-t1/2
      t12 = 1-l2
      t13 = s34r2(z,t11,t12,0.D0,l4)
      t16 = s23r2(z,t11,t12,0.D0,l4)
      t18 = f1(0.D0,z,t11)
      t20 = 1/t8/t13/t16*t18
      t29 = dlog(1-l4)
      t30 = dlog(2.D0)
      t31 = dlog(l1)
      t32 = f2(0.D0,z)
      t33 = dlog(t32)
      t34 = dlog(z)
      t36 = dlog(t12)
      t38 = dlog(-t2)
      t39 = f3(0.D0,z)
      t44 = dlog(2-2*t39+t39*l1*l2)
      t46 = dlog(1-z)
      t48 = dlog(l2)
      t50 = dlog(l4)
      t51 = t29-t30+t31-t33-2*t34+t36+t38+t44-1+4*t46+2*t48+t50
      tncoff7 = 3.D0/8192.D0*t3*t4/l8*t20+3.D0/4096.D0/l3*t3*t4*t20-3.D0
     #/8192.D0*t51*t3*t4*t20
      t2 = dlog(1-l4)
      t3 = dlog(2.D0)
      t4 = dlog(l1)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t7 = dlog(z)
      t9 = 1-l2
      t10 = dlog(t9)
      t11 = l1*l2
      t12 = 2-t11
      t13 = dlog(t12)
      t14 = f3(0.D0,z)
      t19 = dlog(2-2*t14+t14*l1*l2)
      t21 = dlog(1-z)
      t23 = dlog(l2)
      t25 = dlog(l4)
      t26 = t2-t3+t4-t6-2*t7+t10+t13+t19-1+4*t21+2*t23+t25
      t28 = -1/t12
      t30 = 1/l2
      t32 = Diquad721(z,l1,l2,0.D0,l4)
      t35 = 1-t11/2
      t36 = s34r2(z,t35,t9,0.D0,l4)
      t39 = s23r2(z,t35,t9,0.D0,l4)
      t41 = f1(0.D0,z,t35)
      t43 = 1/t32/t36/t39*t41
      tncoff8 = -3.D0/8192.D0*t26*t28*t30*t43-3.D0/8192.D0*t28*t30/(-1+l
     #8)*t43+3.D0/4096.D0/l3*t28*t30*t43
      t1 = l1*l2
      t11 = Diquad721(z,l1,l2,0.D0,l4)
      t14 = 1-t1/2
      t15 = 1-l2
      t16 = s34r2(z,t14,t15,0.D0,l4)
      t19 = s23r2(z,t14,t15,0.D0,l4)
      t21 = f1(0.D0,z,t14)
      tncoff9 = 3.D0/8192.D0/(-2+t1)/l2/l8/(-1+l8)/t11/t16/t19*t21
      t1 = l3/2
      t2 = l1*l2
      t4 = 1-t2/2
      t5 = f1(t1,z,t4)
      t9 = 1-l2
      t10 = s34r2(z,t4,t9,t1,l4)
      t17 = Diquad721(z,l1,l2,t1,l4)
      t20 = s23r2(z,t4,t9,t1,l4)
      tncoff10 = 3.D0/2048.D0*t5/(-2+t2)/t10/l2/(-2+l3)/t17/t20/l3
      t1 = l3/2
      t2 = l1*l2
      t4 = 1-t2/2
      t5 = f1(t1,z,t4)
      t9 = 1-l2
      t10 = s34r2(z,t4,t9,t1,l4)
      t17 = Diquad721(z,l1,l2,t1,l4)
      t20 = s23r2(z,t4,t9,t1,l4)
      tncoff11 = 3.D0/2048.D0*t5/(-2+t2)/t10/l2/(-2+l3)/t17/t20/l3
      t2 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,z,1.D0)
      t11 = 1/t5/t7*t10
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t21 = dlog(z)
      t23 = f3(0.D0,z)
      t25 = dlog(1-t23)
      t28 = dlog(1-l4)
      t30 = dlog(l4)
      t32 = dlog(2.D0)
      t35 = dlog(1-z)
      t37 = dlog(l1)
      tncoff12 = 3.D0/32768.D0/l8*t3*t11+3.D0/16384.D0/l3*t3*t11+(3.D0/3
     #2768.D0+3.D0/32768.D0*t19+3.D0/16384.D0*t21-3.D0/32768.D0*t25-3.D0
     #/32768.D0*t28-3.D0/32768.D0*t30-3.D0/32768.D0*t32-3.D0/8192.D0*t35
     #-3.D0/32768.D0*t37)*t3*t11+3.D0/16384.D0/l2*t3*t11
      t3 = Diquad721(z,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,z,1.D0)
      t12 = 1/t6/t8*t11
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t22 = dlog(z)
      t24 = f3(0.D0,z)
      t26 = dlog(1-t24)
      t29 = dlog(1-l4)
      t31 = dlog(l4)
      t33 = dlog(2.D0)
      t36 = dlog(1-z)
      t38 = dlog(l1)
      tncoff13 = -3.D0/32768.D0/(-1+l8)*t4*t12+3.D0/16384.D0/l3*t4*t12+(
     #3.D0/32768.D0+3.D0/32768.D0*t20+3.D0/16384.D0*t22-3.D0/32768.D0*t2
     #6-3.D0/32768.D0*t29-3.D0/32768.D0*t31-3.D0/32768.D0*t33-3.D0/8192.
     #D0*t36-3.D0/32768.D0*t38)*t4*t12+3.D0/16384.D0/l2*t4*t12
      t5 = Diquad721(z,l1,0.D0,0.D0,l4)
      t8 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff14 = 3.D0/32768.D0/l8/(-1+l8)/t5/t8/t10*t13
      t5 = l3/2
      t6 = Diquad721(z,l1,0.D0,t5,l4)
      t9 = s34r2(z,1.D0,1.D0,t5,l4)
      t11 = s23r2(z,1.D0,1.D0,t5,l4)
      t14 = f1(t5,z,1.D0)
      tncoff15 = 3.D0/8192.D0/(-2+l3)/l3/t6/t9/t11*t14
      t5 = l3/2
      t6 = Diquad721(z,l1,0.D0,t5,l4)
      t9 = s34r2(z,1.D0,1.D0,t5,l4)
      t11 = s23r2(z,1.D0,1.D0,t5,l4)
      t14 = f1(t5,z,1.D0)
      tncoff16 = 3.D0/8192.D0/(-2+l3)/l3/t6/t9/t11*t14
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,0.D0,l4)
      t14 = s23r2(z,t10,t11,0.D0,l4)
      t17 = f1(0.D0,z,t10)
      tncoff17 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,0.D0,l4)
      t14 = s23r2(z,t10,t11,0.D0,l4)
      t17 = f1(0.D0,z,t10)
      tncoff18 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff19 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff20 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l6,l
     #7,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l6,l
     #7,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l
     #7,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l
     #7,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgiq7num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgiq7num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgiq73e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq73e1int(var,ndim,Tcoeff)

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

      t2 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,z,1.D0)
      t11 = 1/t5/t7*t10
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t21 = dlog(z)
      t23 = f3(0.D0,z)
      t25 = dlog(1-t23)
      t28 = dlog(1-l4)
      t30 = dlog(l4)
      t32 = dlog(2.D0)
      t35 = dlog(1-z)
      t37 = dlog(l1)
      tncoff1 = 3.D0/32768.D0/l8*t3*t11+3.D0/16384.D0/l3*t3*t11+(3.D0/32
     #768.D0+3.D0/32768.D0*t19+3.D0/16384.D0*t21-3.D0/32768.D0*t25-3.D0/
     #32768.D0*t28-3.D0/32768.D0*t30-3.D0/32768.D0*t32-3.D0/8192.D0*t35-
     #3.D0/32768.D0*t37)*t3*t11+3.D0/16384.D0/l2*t3*t11
      t3 = Diquad721(z,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,z,1.D0)
      t12 = 1/t6/t8*t11
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t22 = dlog(z)
      t24 = f3(0.D0,z)
      t26 = dlog(1-t24)
      t29 = dlog(1-l4)
      t31 = dlog(l4)
      t33 = dlog(2.D0)
      t36 = dlog(1-z)
      t38 = dlog(l1)
      tncoff2 = -3.D0/32768.D0/(-1+l8)*t4*t12+3.D0/16384.D0/l3*t4*t12+(3
     #.D0/32768.D0+3.D0/32768.D0*t20+3.D0/16384.D0*t22-3.D0/32768.D0*t26
     #-3.D0/32768.D0*t29-3.D0/32768.D0*t31-3.D0/32768.D0*t33-3.D0/8192.D
     #0*t36-3.D0/32768.D0*t38)*t4*t12+3.D0/16384.D0/l2*t4*t12
      t5 = Diquad721(z,l1,0.D0,0.D0,l4)
      t8 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff3 = 3.D0/32768.D0/l8/(-1+l8)/t5/t8/t10*t13
      t5 = l3/2
      t6 = Diquad721(z,l1,0.D0,t5,l4)
      t9 = s34r2(z,1.D0,1.D0,t5,l4)
      t11 = s23r2(z,1.D0,1.D0,t5,l4)
      t14 = f1(t5,z,1.D0)
      tncoff4 = 3.D0/8192.D0/(-2+l3)/l3/t6/t9/t11*t14
      t5 = l3/2
      t6 = Diquad721(z,l1,0.D0,t5,l4)
      t9 = s34r2(z,1.D0,1.D0,t5,l4)
      t11 = s23r2(z,1.D0,1.D0,t5,l4)
      t14 = f1(t5,z,1.D0)
      tncoff5 = 3.D0/8192.D0/(-2+l3)/l3/t6/t9/t11*t14
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,0.D0,l4)
      t14 = s23r2(z,t10,t11,0.D0,l4)
      t17 = f1(0.D0,z,t10)
      tncoff6 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = l1*l2
      t6 = Diquad721(z,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(z,t10,t11,0.D0,l4)
      t14 = s23r2(z,t10,t11,0.D0,l4)
      t17 = f1(0.D0,z,t10)
      tncoff7 = 3.D0/8192.D0/(-2+t1)/l2/t6/t12/t14*t17
      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff8 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff9 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgiq73e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq73e2int(var,ndim,Tcoeff)

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

      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff1 = 3.D0/32768.D0/t1/t3/t6*t8
      t1 = Diquad721(z,l1,0.D0,0.D0,l4)
      t3 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      tncoff2 = 3.D0/32768.D0/t1/t3/t6*t8

CCC Numerator functions

      tcnum1 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgiq73e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgiq73e3int(var,ndim,Tcoeff)

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

      qgiq73e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

