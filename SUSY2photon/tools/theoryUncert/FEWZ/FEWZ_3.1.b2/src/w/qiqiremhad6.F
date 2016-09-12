CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm6e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum3e0
      external qiqirmnum3e0
      double precision qiqirmnum3e1
      external qiqirmnum3e1
      double precision qiqirmnum3e2
      external qiqirmnum3e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t11 = dlog(1-l1)
      t14 = dlog(l4)
      t18 = dlog(-t1)
      t21 = dlog(l1)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t28 = dlog(z)
      t37 = dlog(2.D0)
      t40 = -6*t3*t7-6*t3*t11-6*t3*t14-24*t18*t11-6*t3*t21+24*t18*t25+12
     #*t11*t28+12*t3*t28-24*t14*t18+6*t3*t25-24*t3*t37
      t49 = t11**2
      t51 = 0.3141592653589793D1**2
      t52 = t7**2
      t56 = t18**2
      t58 = t37**2
      t60 = t25**2
      t62 = t28**2
      t64 = 48*t28*t37-6*t21*t7+6*t21*t25-24*t21*t37-3*t49+t51-3*t52-24*
     #t7*t37-48*t56-48*t58-3*t60-12*t62
      t66 = t21**2
      t68 = t3**2
      t84 = t14**2
      t88 = -3*t66-3*t68+6*t14*t25+12*t7*t28+12*t14*t28-24*t14*t37+24*t2
     #5*t37+12*t21*t28+6*t7*t25-3*t84-6*t7*t11
      t113 = -24*t11*t37-6*t14*t7-6*t14*t11+6*t11*t25-6*t14*t21-6*t21*t1
     #1-96*t18*t37-24*t18*t7+48*t18*t28-24*t21*t18-12*t25*t28-24*t3*t18
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = dlog(l8)
      t129 = 4*t37
      t130 = 4*t18
      t144 = 2*t28
      t145 = -t128+t25-t21-t129-t7-t130-t11-t14-t3+t11*l8+t128*l8+t14*l8
     #+t3*l8+4*t37*l8-2*t28*l8+4*t18*l8-t25*l8+t21*l8+t7*l8+t144
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = t1*(t40+t64+t88+t113)*t125/6144+t1*t145/(-1+l8)*t150*t11
     #9*t122*t124/1024-t1*(-t14-2*t156-t21-t130-t7-t11+t25-t129+t144-t3)
     #*t160*t125/512-t1*t160*t150*t125/512-t1*(-t21-t130-t7-t11+t25-t129
     #+t144-t168-t14-t3)*t171*t125/1024-t1*t150*t171*t125/1024-t1*t119*t
     #171*t160*t122*t124/512
      t1 = -1+z
      t2 = dlog(l4)
      t3 = t2**2
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t10 = dlog(-t1)
      t12 = dlog(1-l1)
      t16 = dlog(1-l4)
      t17 = dlog(l1)
      t22 = dlog(z)
      t31 = dlog(2.D0)
      t36 = 3*t3-6*t2*t6+24*t10*t12+6*t16*t17-24*t10*t6-12*t12*t22-12*t1
     #6*t22+24*t2*t10-6*t16*t6+24*t16*t31-48*t22*t31
      t37 = f3(0.D0,z)
      t40 = dlog(1-l1*t37)
      t65 = 6*t17*t40-6*t17*t6+24*t17*t31+24*t40*t31-12*t40*t22-12*t2*t2
     #2+24*t2*t31-24*t6*t31-12*t17*t22-6*t40*t6+6*t40*t12+24*t12*t31
      t87 = t22**2
      t89 = 6*t2*t40+6*t2*t12-6*t12*t6+6*t2*t17+6*t17*t12+96*t10*t31+12*
     #t6*t22+24*t17*t10-48*t10*t22+6*t16*t2+12*t87
      t92 = t12**2
      t94 = t40**2
      t98 = 0.3141592653589793D1**2
      t99 = t31**2
      t103 = t10**2
      t107 = t6**2
      t109 = t17**2
      t111 = t16**2
      t113 = 6*t16*t12+3*t92+3*t94+6*t16*t40-t98+48*t99+24*t16*t10+48*t1
     #03+24*t10*t40+3*t107+3*t109+3*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t133 = 1-l8
      t134 = dlog(t133)
      t146 = -1/t133
      t154 = 4*t10
      t155 = 4*t31
      t156 = 2*t22
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = -t1*(t36+t65+t89+t113)*t125/6144+t1*(-t16*l8-4*t10*l8-t1
     #7*l8-t2*l8-t134*l8-t12*l8+t6*l8-t40*l8-4*t31*l8+2*t22*l8)*t146/l8*
     #t119*t122*t124/1024+t1*(t154+t17+t12+t40-t6+t155-t156+2*t157+t2+t1
     #6)*t161*t125/512+t1*t161*t146*t125/512+t1*(t154+t17+t12+t40-t6+t15
     #5-t156+t169+t2+t16)*t172*t125/1024+t1*t172*t146*t125/1024-t1*t119*
     #t172*t161*t122*t124/512
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = 1-l8
      t10 = dlog(t9)
      t11 = dlog(l8)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t15 = dlog(1-l1)
      t16 = dlog(z)
      t18 = dlog(2.D0)
      t20 = dlog(l1)
      t22 = dlog(-t1)
      t24 = t2+t4+t8+t10+t11-t13+t15-2*t16+4*t18+t20+4*t22
      t27 = -1/t9
      t29 = 1/l8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t27*t39/1024-t1/l3*t27*t39/512-t1*t29/l2*t27*t32*
     #t38/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(l4)
      t10 = dlog(-t1)
      t13 = dlog(1-l1)
      t14 = dlog(l1)
      t15 = dlog(l3)
      t18 = dlog(-t2)
      t20 = f3(l3,z)
      t23 = dlog(1-l1*t20)
      t25 = dlog(1-l4)
      t26 = dlog(2.D0)
      t28 = -2*t4-t7+t8+4*t10+t13+t14+2*t15+2*t18+t23+t25+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(l3,z,l1)
      t40 = t30*t33*t37*t38
      t42 = 1-l8
      tncoff4 = t3*t28*t40/512-t3*t30/l8*t33*t37*t38/512-t3/l2*t40/512
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(l4)
      t10 = dlog(-t1)
      t13 = dlog(1-l1)
      t14 = dlog(l1)
      t15 = dlog(l3)
      t18 = dlog(-t2)
      t20 = f3(l3,z)
      t23 = dlog(1-l1*t20)
      t25 = dlog(1-l4)
      t26 = dlog(2.D0)
      t28 = -2*t4-t7+t8+4*t10+t13+t14+2*t15+2*t18+t23+t25+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = t3*t28*t40/512+t3*t30/(-1+l8)*t33*t39/512-t3/l2*t40/512
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/(-1+l8)/l8/t11/t15*t17/512
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = dlog(z)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(-t1)
      t13 = dlog(1-l1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(1-l2)
      t20 = dlog(l2)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = t2+4*t3-2*t5-t8+4*t10+t13+t17+t19+t20+t21+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t39 = 1-l8
      tncoff7 = -t1*t24*t26*t30*t33*t35/1024+t1/l8*t26*t30*t33*t35/1024+
     #t1*t33*t30*t26/l3*t35/512
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = dlog(z)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(-t1)
      t13 = dlog(1-l1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(1-l2)
      t20 = dlog(l2)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = t2+4*t3-2*t5-t8+4*t10+t13+t17+t19+t20+t21+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = -t1*t24*t26*t36/1024-t1*t26/(-1+l8)*t36/1024+t1*t33*t30*
     #t26/l3*t35/512
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = (-1+z)/l8/l2/(-1+l8)/t9/t13*t15/1024
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/512
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/512
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = dlog(l1)
      t8 = dlog(1-l1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(l4)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t20 = dlog(2.D0)
      t22 = dlog(z)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = t1*(24*t3+6*t5+6*t8+6*t13+6*t15-6*t18+24*t20-12*t22+6*t
     #25)*t37/6144-t1/l8*t31*t34*t36/1024-t1/l3*t37/512-t1/l2*t37/1024
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = dlog(l1)
      t8 = dlog(1-l1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(2.D0)
      t20 = dlog(z)
      t23 = dlog(1-l4)
      t25 = dlog(l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-24*t3-6*t5-6*t8-6*t13+6*t16-24*t18+12*t20-6*t23-6
     #*t25)*t37/6144+t1/(-1+l8)*t37/1024-t1/l3*t37/512-t1/l2*t37/1024
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff15 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff16 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff17 = (-1+z)/l2/t5/t8*t11/1024
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff18 = (-1+z)/l2/t5/t8*t11/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff19 = -(-1+z)/t3/t7*t9/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff20 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm6e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm6e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum3e0
      external qiqirmnum3e0
      double precision qiqirmnum3e1
      external qiqirmnum3e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = -1+z
      t3 = dlog(-t1)
      t6 = dlog(1-l1)
      t7 = dlog(l1)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-l4)
      t18 = dlog(z)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-4*t3-t6-t7-4*t8-t13+t15-t17+2*t18-t20)*t31/1024-t1
     #/l8*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(-t1)
      t7 = dlog(l4)
      t9 = dlog(1-l1)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t12 = dlog(l1)
      t13 = dlog(z)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t19 = dlog(2.D0)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t3+4*t5+t7+t9-t11+t12-2*t13+t18+4*t19)*t31/1024+t1/(
     #-1+l8)*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff4 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff5 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff6 = (-1+z)/l2/t5/t8*t11/1024
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff7 = (-1+z)/l2/t5/t8*t11/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -(-1+z)/t3/t7*t9/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm6e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm6e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum3e0
      external qiqirmnum3e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -(-1+z)/t3/t7*t9/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm6e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

