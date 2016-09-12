CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqjrm2e0int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum2e0
      external qiqjrmnum2e0
      double precision qiqjrmnum2e1
      external qiqjrmnum2e1
      double precision qiqjrmnum2e2
      external qiqjrmnum2e2
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
      t4 = dlog(z)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t14 = dlog(1-l1)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t22 = dlog(-t1)
      t25 = dlog(2.D0)
      t30 = dlog(l1)
      t37 = dlog(l4)
      t40 = -12*t3*t4+6*t3*t10+6*t3*t14-6*t3*t18-48*t22*t4+96*t22*t25-12
     #*t10*t4+6*t30*t14+24*t10*t25-24*t18*t22-6*t37*t18
      t65 = -12*t30*t4+24*t30*t25+24*t37*t25-6*t18*t10+24*t37*t22-12*t37
     #*t4+6*t37*t14+6*t37*t10+6*t30*t10+24*t30*t22+12*t18*t4+6*t37*t30
      t67 = t14**2
      t69 = t3**2
      t71 = t30**2
      t73 = t22**2
      t89 = 3*t67+3*t69+3*t71+48*t73+24*t3*t22+24*t14*t25+24*t10*t22+24*
     #t14*t22-12*t14*t4-6*t14*t18+6*t14*t10
      t100 = t10**2
      t102 = t4**2
      t104 = 0.3141592653589793D1**2
      t105 = t25**2
      t109 = t37**2
      t111 = t18**2
      t113 = -24*t18*t25-6*t30*t18-48*t4*t25+24*t3*t25+6*t3*t30+3*t100+1
     #2*t102-t104+48*t105+6*t3*t37+3*t109+3*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = dlog(l8)
      t137 = 2*t4
      t138 = 4*t25
      t139 = 4*t22
      t145 = -t128-t14+4*t22*l8+t14*l8+t30*l8+t37*l8-t18*l8+t10*l8+t128*
     #l8-t10-t37-t3+t137-t138-t139+t18-t30+t3*l8-2*t4*l8+4*t25*l8
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = -t1*(t40+t65+t89+t113)*t125/6144+t1*t145/(-1+l8)*t150*t1
     #19*t122*t124/1024+t1*(t37+t3+2*t156+t30+t10+t139-t137+t138-t18+t14
     #)*t160*t125/512-t1*t160*t150*t125/512+t1*(t3+t37+t168+t30+t10+t139
     #-t137+t138-t18+t14)*t171*t125/1024-t1*t171*t150*t125/1024-t160*t17
     #1*t1*t125/512
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = f3(0.D0,z)
      t9 = dlog(1-l1*t6)
      t13 = dlog(1-l4)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t19 = dlog(-t1)
      t23 = dlog(1-l1)
      t28 = dlog(z)
      t31 = dlog(l4)
      t40 = 24*t2*t3+6*t2*t9+6*t13*t2-24*t17*t19+24*t23*t3+24*t9*t19-48*
     #t19*t28-12*t31*t28+6*t13*t9-24*t17*t3+96*t19*t3
      t51 = 0.3141592653589793D1**2
      t52 = t3**2
      t58 = t17**2
      t60 = t31**2
      t62 = t13**2
      t64 = -12*t13*t28-6*t31*t17-12*t23*t28+6*t13*t31-12*t9*t28-t51+48*
     #t52+24*t9*t3-48*t28*t3+3*t58+3*t60+3*t62
      t70 = t23**2
      t76 = t9**2
      t78 = t19**2
      t80 = t2**2
      t82 = t28**2
      t88 = 6*t23*t9-6*t23*t17+3*t70-6*t2*t17+24*t2*t19+3*t76+48*t78+3*t
     #80+12*t82+24*t13*t3+24*t31*t19
      t113 = 12*t17*t28+6*t31*t23+24*t23*t19+6*t13*t23-12*t2*t28+6*t31*t
     #9+24*t13*t19+6*t31*t2-6*t17*t9+6*t2*t23-6*t13*t17+24*t31*t3
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t130 = 1-l8
      t131 = dlog(t130)
      t146 = -1/t130
      t154 = dlog(l3)
      t156 = 4*t19
      t157 = 2*t28
      t158 = 4*t3
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = -t1*(t40+t64+t88+t113)*t125/6144+t1*(-t13*l8-t2*l8-t131*
     #l8-t31*l8+t17*l8-t9*l8+2*t28*l8-4*t3*l8-4*t19*l8-t23*l8)*t146/l8*t
     #119*t122*t124/1024+t1*(t31+t13+2*t154+t2+t9+t156-t157+t158-t17+t23
     #)*t161*t125/512+t1*t161*t146*t125/512+t1*(t13+t31+t169+t2+t9+t156-
     #t157+t158-t17+t23)*t172*t125/1024+t1*t172*t146*t125/1024-t161*t172
     #*t1*t125/512
      t1 = -1+z
      t2 = dlog(2.D0)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = 1-l8
      t7 = dlog(t6)
      t8 = dlog(l4)
      t10 = dlog(1-l4)
      t11 = dlog(l8)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(l1)
      t18 = dlog(-t1)
      t20 = dlog(z)
      t23 = dlog(1-l1)
      t24 = 4*t2-t5+t7+t8+t10+t11+t15+t16+4*t18-2*t20+t23
      t27 = -1/t6
      t29 = 1/l8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t27*t39/1024-t1/l3*t27*t39/512-t1/l2*t29*t27*t32*
     #t38/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t6 = dlog(2.D0)
      t8 = f3(l3,z)
      t11 = dlog(1-l1*t8)
      t12 = dlog(l4)
      t14 = dlog(1-l4)
      t15 = dlog(l1)
      t16 = dlog(l3)
      t19 = dlog(-t2)
      t21 = f2(l3,z)
      t22 = dlog(t21)
      t24 = dlog(1-l1)
      t26 = dlog(-t1)
      t28 = -2*t4+4*t6+t11+t12+t14+t15+2*t16+2*t19-t22+t24+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff4 = t3*t28*t40/512-t3*t30/l8*t33*t39/512-t3/l2*t40/512
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t6 = dlog(2.D0)
      t8 = f3(l3,z)
      t11 = dlog(1-l1*t8)
      t12 = dlog(l4)
      t14 = dlog(1-l4)
      t15 = dlog(l1)
      t16 = dlog(l3)
      t19 = dlog(-t2)
      t21 = f2(l3,z)
      t22 = dlog(t21)
      t24 = dlog(1-l1)
      t26 = dlog(-t1)
      t28 = -2*t4+4*t6+t11+t12+t14+t15+2*t16+2*t19-t22+t24+4*t26
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
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(l2)
      t7 = dlog(1-l1)
      t8 = dlog(2.D0)
      t10 = dlog(l1)
      t11 = dlog(z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l2)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(-t1)
      t24 = t3+t4+t5+t7+4*t8+t10-2*t11-t14+t16+t20+4*t22
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t34 = 1/t29/t32
      t35 = f1(0.D0,z,l1)
      t36 = t34*t35
      t39 = t1*t26
      tncoff7 = -t1*t24*t26*t36/1024+t39/l8*t36/1024+t34/l3*t39*t35/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(l2)
      t7 = dlog(1-l1)
      t8 = dlog(2.D0)
      t10 = dlog(l1)
      t11 = dlog(z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l2)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(-t1)
      t24 = t3+t4+t5+t7+4*t8+t10-2*t11-t14+t16+t20+4*t22
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t34 = 1/t29/t32
      t35 = f1(0.D0,z,l1)
      t36 = t34*t35
      t39 = t1*t26
      tncoff8 = -t1*t24*t26*t36/1024-t39/(-1+l8)*t36/1024+t34/l3*t39*t35
     #/512
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = (-1+z)/l2/l8/(-1+l8)/t9/t13*t15/1024
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = (l3-1)*(-1+z)/l2/l3/t8/t12*t14/512
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = (l3-1)*(-1+z)/l2/l3/t8/t12*t14/512
      t1 = -1+z
      t2 = dlog(l1)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(1-l1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t16 = dlog(-t1)
      t18 = dlog(z)
      t20 = dlog(2.D0)
      t23 = dlog(1-l4)
      t25 = dlog(l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = -t1*(-6*t2+6*t5-6*t8-6*t13-24*t16+12*t18-24*t20-6*t23-6
     #*t25)*t37/6144-t1/l8*t31*t34*t36/1024-t1/l3*t37/512-t1/l2*t37/1024
      t1 = -1+z
      t2 = dlog(l1)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(1-l1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t16 = dlog(-t1)
      t18 = dlog(z)
      t20 = dlog(2.D0)
      t23 = dlog(1-l4)
      t25 = dlog(l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t2+6*t5-6*t8-6*t13-24*t16+12*t18-24*t20-6*t23-6
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

      tcnum1 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqjrmnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqjrmnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqjrmnum2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqjrmnum2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqjrm2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm2e1int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum2e0
      external qiqjrmnum2e0
      double precision qiqjrmnum2e1
      external qiqjrmnum2e1
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
      t5 = dlog(1-l1)
      t6 = dlog(l4)
      t7 = dlog(l1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t17 = dlog(-t1)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(t3+t5+t6+t7+t11-2*t12+4*t14+4*t17-t20)*t31/1024-t1/l
     #8*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(1-l1)
      t6 = dlog(l4)
      t7 = dlog(l1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t17 = dlog(-t1)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t3+t5+t6+t7+t11-2*t12+4*t14+4*t17-t20)*t31/1024+t1/(
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

      tcnum1 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqjrm2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm2e2int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum2e0
      external qiqjrmnum2e0
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

      tcnum1 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqjrm2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

