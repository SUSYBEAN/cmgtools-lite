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
      double precision qiqjrmnum1e0
      external qiqjrmnum1e0
      double precision qiqjrmnum1e1
      external qiqjrmnum1e1
      double precision qiqjrmnum1e2
      external qiqjrmnum1e2
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
      t2 = dlog(l1)
      t3 = t2**2
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t7 = t6**2
      t9 = dlog(2.D0)
      t11 = dlog(-t1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t20 = dlog(z)
      t28 = dlog(1-l4)
      t34 = dlog(1-l1)
      t39 = -3*t3-3*t7-72*t9*t11-24*t17*t11+12*t17*t20-6*t17*t2+6*t6*t2-
     #6*t28*t17-24*t11*t2-18*t9*t34-18*t9*t2
      t46 = dlog(l4)
      t65 = 24*t11*t6+12*t34*t20+18*t9*t6-6*t28*t46+12*t28*t20-24*t28*t1
     #1-6*t28*t34+12*t46*t20-6*t28*t2+6*t28*t6-18*t28*t9-24*t46*t11
      t89 = 6*t46*t6-6*t2*t34-24*t11*t34+48*t11*t20-6*t17*t34-18*t17*t9+
     #6*t17*t6-6*t46*t17+12*t2*t20+6*t6*t34-6*t46*t2
      t90 = 0.3141592653589793D1**2
      t95 = t9**2
      t97 = t17**2
      t99 = t46**2
      t101 = t28**2
      t105 = t34**2
      t109 = t20**2
      t111 = t11**2
      t113 = t90+36*t20*t9-12*t6*t20-27*t95-3*t97-3*t99-3*t101-6*t46*t34
     #-3*t105-18*t46*t9-12*t109-48*t111
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,1.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 4*t11
      t129 = 3*t9
      t130 = 2*t20
      t139 = dlog(l8)
      t145 = -t46-t128-t2+t6-t129+t130+4*t11*l8-t6*l8+t2*l8+t34*l8-2*t20
     #*l8+t28*l8+t139*l8+t46*l8+3*t9*l8+t17*l8-t28-t139-t34-t17
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = t1*(t39+t65+t89+t113)*t125/12288+t1*t145/(-1+l8)*t150*t1
     #19*t122*t124/2048+t1*(2*t156+t17+t129+t128-t6+t2+t34-t130+t46+t28)
     #*t160*t125/1024-t1*t160*t150*t125/1024+(t28+t46+t168+t128+t17+t129
     #+t2-t130+t34-t6)*t1*t171*t125/2048-t1*t150*t171*t125/2048-t1*t171*
     #t160*t125/1024
      t1 = -1+z
      t2 = dlog(l1)
      t4 = dlog(1-l1)
      t8 = dlog(-t1)
      t11 = dlog(z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t20 = dlog(2.D0)
      t23 = f2(0.D0,z)
      t24 = dlog(t23)
      t27 = dlog(l4)
      t38 = 6*t2*t4+24*t8*t4-48*t8*t11+6*t17*t4+18*t17*t20-6*t17*t24+6*t
     #27*t17-12*t2*t11-6*t24*t4+12*t24*t11-36*t11*t20
      t50 = dlog(1-l4)
      t65 = 72*t20*t8+24*t17*t8-12*t17*t11+6*t17*t2-6*t24*t2+6*t50*t17+2
     #4*t8*t2+18*t20*t4+18*t20*t2-24*t8*t24-12*t4*t11+24*t50*t8
      t67 = t11**2
      t69 = t8**2
      t71 = t2**2
      t75 = 0.3141592653589793D1**2
      t76 = t20**2
      t78 = t24**2
      t88 = 12*t67+48*t69+3*t71-18*t20*t24-t75+27*t76+3*t78-12*t50*t11+6
     #*t50*t4-12*t27*t11+6*t50*t2
      t103 = t17**2
      t105 = t27**2
      t109 = t4**2
      t111 = t50**2
      t113 = -6*t50*t24+18*t50*t20+24*t27*t8+6*t27*t4+6*t27*t2-6*t27*t24
     #+18*t27*t20+3*t103+3*t105+6*t50*t27+3*t109+3*t111
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,1.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t130 = 1-l8
      t131 = dlog(t130)
      t146 = -1/t130
      t154 = dlog(l3)
      t156 = 3*t20
      t157 = 4*t8
      t158 = 2*t11
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = -t1*(t38+t65+t88+t113)*t125/12288+t1*(-t50*l8-t17*l8-t13
     #1*l8-t27*l8-3*t20*l8-4*t8*l8+t24*l8-t2*l8-t4*l8+2*t11*l8)*t146/l8*
     #t119*t122*t124/2048+t1*(2*t154+t17+t156+t157-t24+t2+t4-t158+t27+t5
     #0)*t161*t125/1024+t1*t161*t146*t125/1024+(t50+t27+t169+t157+t17+t1
     #56+t2-t158+t4-t24)*t1*t172*t125/2048+t1*t146*t172*t125/2048-t1*t17
     #2*t161*t125/1024
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(z)
      t11 = dlog(1-l1)
      t12 = dlog(l1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(-t1)
      t18 = dlog(2.D0)
      t20 = dlog(l8)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+t7-2*t8+t11+t12-t14+4*t16+3*t18+t20+t22+t23
      t27 = -1/t2
      t29 = 1/l8
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      tncoff3 = t1*t24*t27*t29*t32*t38/2048-t1/l3*t29*t27*t32*t38/1024-t
     #1*t27*t29/l2*t32*t38/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(l3)
      t11 = dlog(l1)
      t12 = f3(l3,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t20 = dlog(l4)
      t22 = dlog(-t2)
      t25 = dlog(1-l4)
      t27 = dlog(1-l1)
      t28 = 4*t5-t8+2*t9+t11+t15+3*t16-2*t18+t20+2*t22+t25+t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(l3,z,l1)
      t40 = t30*t33*t37*t38
      t42 = 1-l8
      tncoff4 = t3*t28*t40/1024-t3*t30/l8*t33*t37*t38/1024-t3/l2*t40/102
     #4
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(l3)
      t11 = dlog(l1)
      t12 = f3(l3,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t20 = dlog(l4)
      t22 = dlog(-t2)
      t25 = dlog(1-l4)
      t27 = dlog(1-l1)
      t28 = 4*t5-t8+2*t9+t11+t15+3*t16-2*t18+t20+2*t22+t25+t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = t3*t28*t40/1024+t3*t30/(-1+l8)*t33*t39/1024-t3/l2*t40/10
     #24
      t10 = s23r2(z,l1,1.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,1.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/l8/(-1+l8)/t11/t15*t17/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(l2)
      t6 = dlog(2.D0)
      t9 = dlog(-t1)
      t11 = f3(0.D0,z)
      t14 = dlog(1-l1*t11)
      t16 = dlog(1-l1)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t19 = -l2+2
      t20 = dlog(t19)
      t21 = dlog(l1)
      t22 = dlog(z)
      t24 = t3+t4+t5+2*t6+4*t9+t14+t16-t18+t20+t21-2*t22
      t27 = -1/t19
      t29 = 1/l2
      t31 = 1-l2/2
      t32 = s23r2(z,l1,t31,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = s1vr2(z,l1,t31,0.D0,l4)
      t37 = t36**2
      t38 = 1/t37
      t39 = f1(0.D0,z,l1)
      t44 = 1-l8
      t52 = t34*t38
      tncoff7 = t1*t24*t27*t29*t34*t38*t39/1024-t1*t27*t29/l8*t52*t39/10
     #24-t1*t27*t29*t52/l3*t39/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(l2)
      t6 = dlog(2.D0)
      t9 = dlog(-t1)
      t11 = f3(0.D0,z)
      t14 = dlog(1-l1*t11)
      t16 = dlog(1-l1)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t19 = -l2+2
      t20 = dlog(t19)
      t21 = dlog(l1)
      t22 = dlog(z)
      t24 = t3+t4+t5+2*t6+4*t9+t14+t16-t18+t20+t21-2*t22
      t27 = -1/t19
      t29 = 1/l2
      t31 = 1-l2/2
      t32 = s23r2(z,l1,t31,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = s1vr2(z,l1,t31,0.D0,l4)
      t37 = t36**2
      t38 = 1/t37
      t39 = f1(0.D0,z,l1)
      t40 = t38*t39
      t45 = t1*t27*t29
      tncoff8 = t1*t24*t27*t29*t34*t40/1024+t45/(-1+l8)*t34*t40/1024-t45
     #*t38*t34/l3*t39/512
      t12 = 1-l2/2
      t13 = s23r2(z,l1,t12,0.D0,l4)
      t14 = t13**2
      t17 = s1vr2(z,l1,t12,0.D0,l4)
      t18 = t17**2
      t20 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/(l2-2)/l2/l8/(-1+l8)/t14/t18*t20/1024
      t11 = 1-l2/2
      t12 = s23r2(z,l1,t11,l3,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,l3,l4)
      t17 = t16**2
      t19 = f1(l3,z,l1)
      tncoff10 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t13/t17*t19/512
      t11 = 1-l2/2
      t12 = s23r2(z,l1,t11,l3,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,l3,l4)
      t17 = t16**2
      t19 = f1(l3,z,l1)
      tncoff11 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t13/t17*t19/512
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t5 = dlog(l4)
      t8 = dlog(1-l4)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(2.D0)
      t18 = dlog(-t1)
      t21 = dlog(1-l1)
      t23 = dlog(z)
      t25 = dlog(l1)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,1.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = t1*(-6*t3+6*t5+6*t8+6*t13+18*t15+24*t18+6*t21-12*t23+6*
     #t25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/204
     #8
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(l4)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t12 = dlog(2.D0)
      t15 = dlog(-t1)
      t18 = dlog(1-l1)
      t20 = dlog(z)
      t22 = dlog(l1)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,1.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t3-6*t5-6*t10-18*t12-24*t15-6*t18+12*t20-6*t22+
     #6*t25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff15 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff16 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff17 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff18 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff19 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff20 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqjrmnum1e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqjrmnum1e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqjrmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqjrmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

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
      double precision qiqjrmnum1e0
      external qiqjrmnum1e0
      double precision qiqjrmnum1e1
      external qiqjrmnum1e1
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
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t7 = dlog(-t1)
      t9 = dlog(l1)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t16 = dlog(l4)
      t18 = dlog(1-l1)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,1.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-t5-4*t7-t9-3*t10+2*t12+t15-t16-t18-t20)*t31/2048-t
     #1/l8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = dlog(l4)
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t8 = dlog(1-l4)
      t10 = dlog(-t1)
      t12 = dlog(l1)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,1.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t2+t6+t8+4*t10+t12+3*t13-2*t15-t18+t20)*t31/2048+t1/
     #(-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff4 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff5 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff6 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff7 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

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
      double precision qiqjrmnum1e0
      external qiqjrmnum1e0
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

      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqjrm2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

