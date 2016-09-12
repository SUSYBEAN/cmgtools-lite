CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqjrm1e0int(var,ndim,Tcoeff)

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
      t2 = dlog(l4)
      t3 = t2**2
      t6 = dlog(1-l4)
      t7 = t6**2
      t10 = dlog(-t1)
      t12 = dlog(1-l1)
      t15 = dlog(z)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t22 = dlog(2.D0)
      t27 = dlog(l1)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t38 = t12**2
      t40 = 3*t3+3*t7+24*t10*t12-48*t10*t15+18*t21*t22+24*t21*t10-12*t27
     #*t15-24*t10*t31-6*t31*t27+12*t31*t15+3*t38
      t65 = 24*t2*t10-6*t2*t31-12*t12*t15-6*t21*t31+72*t10*t22+18*t27*t2
     #2-12*t21*t15+6*t2*t21+6*t21*t27+6*t21*t12+24*t10*t27-12*t6*t15
      t89 = 6*t6*t2+6*t6*t27+6*t6*t12-6*t6*t31+18*t6*t22+6*t6*t21+24*t6*
     #t10+6*t27*t12-36*t15*t22+6*t2*t12-12*t2*t15
      t98 = t15**2
      t100 = t10**2
      t102 = t27**2
      t104 = t31**2
      t108 = 0.3141592653589793D1**2
      t109 = t22**2
      t111 = t21**2
      t113 = 18*t2*t22+6*t2*t27-18*t31*t22+18*t12*t22+12*t98+48*t100+3*t
     #102+3*t104-6*t31*t12-t108+27*t109+3*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t129 = 2*t15
      t130 = 4*t10
      t131 = 3*t22
      t132 = dlog(l8)
      t145 = -t27+t2*l8+t129-t130-t131-t6-t132+t6*l8-t31*l8+3*t22*l8+t13
     #2*l8+t27*l8+t31-t21-t12-t2+t12*l8+4*t10*l8+t21*l8-2*t15*l8
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = -t1*(t40+t65+t89+t113)*t125/12288+t1*t145/(-1+l8)*t150*t
     #119*t122*t124/2048+t1*(2*t156+t21-t31+t27+t12-t129+t130+t131+t6+t2
     #)*t160*t125/1024-t1*t160*t150*t125/1024+(t12+t27+t21+t130-t31-t129
     #+t131+t168+t2+t6)*t1*t171*t125/2048-t1*t150*t171*t125/2048-t1*t171
     #*t160*t125/1024
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(-t1)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t12 = dlog(1-l1)
      t13 = dlog(z)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t27 = dlog(1-l4)
      t28 = dlog(l1)
      t37 = dlog(2.D0)
      t40 = 24*t2*t4-6*t2*t8-12*t12*t13-6*t19*t8+24*t4*t12-48*t4*t13+6*t
     #27*t28+24*t19*t4+6*t2*t27-12*t27*t13+18*t19*t37
      t57 = t12**2
      t59 = t13**2
      t61 = t4**2
      t63 = t28**2
      t65 = 24*t4*t28+6*t19*t12+6*t19*t28+6*t2*t19-12*t19*t13+18*t28*t37
     #+72*t4*t37-12*t28*t13+3*t57+12*t59+48*t61+3*t63
      t67 = 0.3141592653589793D1**2
      t68 = t37**2
      t70 = t8**2
      t80 = t19**2
      t88 = -t67+27*t68+3*t70+12*t8*t13-6*t8*t28-24*t4*t8-6*t8*t12+3*t80
     #+18*t12*t37-18*t8*t37+6*t2*t28
      t97 = t2**2
      t109 = t27**2
      t113 = -12*t2*t13+18*t2*t37+6*t2*t12-36*t13*t37+3*t97+24*t27*t4+6*
     #t28*t12+6*t27*t19+18*t27*t37+6*t27*t12+3*t109-6*t27*t8
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t139 = 1-l8
      t140 = dlog(t139)
      t146 = -1/t139
      t154 = dlog(l3)
      t156 = 2*t13
      t157 = 4*t4
      t158 = 3*t37
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = -t1*(t40+t65+t88+t113)*t125/12288+t1*(t8*l8-3*t37*l8-4*t
     #4*l8-t19*l8-t28*l8-t12*l8-t27*l8+2*t13*l8-t140*l8-t2*l8)*t146/l8*t
     #119*t122*t124/2048+t1*(2*t154+t19-t8+t28+t12-t156+t157+t158+t27+t2
     #)*t161*t125/1024+t1*t161*t146*t125/1024+(t12+t28+t19+t157-t8-t156+
     #t158+t169+t2+t27)*t1*t172*t125/2048+t1*t146*t172*t125/2048-t1*t172
     #*t161*t125/1024
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t4 = dlog(2.D0)
      t6 = dlog(l1)
      t7 = 1-l8
      t8 = dlog(t7)
      t9 = dlog(l4)
      t10 = dlog(l8)
      t12 = dlog(1-l4)
      t14 = dlog(1-l1)
      t15 = dlog(z)
      t18 = dlog(-t1)
      t20 = f3(0.D0,z)
      t23 = dlog(1-l1*t20)
      t24 = -t3+3*t4+t6+t8+t9+t10+t12+t14-2*t15+4*t18+t23
      t27 = -1/t7
      t29 = 1/l8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      tncoff3 = t1*t24*t27*t29*t32*t38/2048-t1/l3*t29*t27*t32*t38/1024-t
     #1*t29*t27/l2*t32*t38/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l3)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = f3(l3,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l1)
      t16 = dlog(l1)
      t17 = dlog(z)
      t20 = dlog(1-l4)
      t22 = dlog(-t2)
      t24 = dlog(l4)
      t26 = dlog(-t1)
      t28 = 2*t4-t7+3*t8+t13+t15+t16-2*t17+t20+2*t22+t24+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff4 = t3*t28*t40/1024-t3*t30/l8*t33*t39/1024-t3/l2*t40/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l3)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = f3(l3,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l1)
      t16 = dlog(l1)
      t17 = dlog(z)
      t20 = dlog(1-l4)
      t22 = dlog(-t2)
      t24 = dlog(l4)
      t26 = dlog(-t1)
      t28 = 2*t4-t7+3*t8+t13+t15+t16-2*t17+t20+2*t22+t24+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = t3*t28*t40/1024+t3*t30/(-1+l8)*t33*t39/1024-t3/l2*t40/10
     #24
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/l8/(-1+l8)/t11/t15*t17/1024
      t1 = -1+z
      t2 = -l2+2
      t3 = dlog(t2)
      t4 = dlog(z)
      t7 = dlog(-t1)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t19 = dlog(1-l1)
      t20 = dlog(l2)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = t3-2*t4+4*t7+t12+2*t13-t16+t17+t19+t20+t21+t23
      t27 = -1/t2
      t29 = 1/l2
      t30 = l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,t30,0.D0,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(0.D0,z,l1)
      t39 = t37*t38
      tncoff7 = t1*t24*t27*t29*t33*t39/1024-t1*t29*t27/l8*t33*t39/1024-t
     #1*t37*t33*t27*t29/l3*t38/512
      t1 = -1+z
      t2 = -l2+2
      t3 = dlog(t2)
      t4 = dlog(z)
      t7 = dlog(-t1)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t19 = dlog(1-l1)
      t20 = dlog(l2)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = t3-2*t4+4*t7+t12+2*t13-t16+t17+t19+t20+t21+t23
      t27 = -1/t2
      t29 = 1/l2
      t30 = l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,t30,0.D0,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(0.D0,z,l1)
      t39 = t37*t38
      tncoff8 = t1*t24*t27*t29*t33*t39/1024+t1*t29/(-1+l8)*t27*t33*t39/1
     #024-t1*t37*t33*t27*t29/l3*t38/512
      t11 = l2/2
      t12 = s23r2(z,l1,t11,0.D0,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,0.D0,l4)
      t17 = t16**2
      t19 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/l2/(-1+l8)/(l2-2)/l8/t13/t17*t19/1024
      t10 = l2/2
      t11 = s23r2(z,l1,t10,l3,l4)
      t12 = t11**2
      t15 = s1vr2(z,l1,t10,l3,l4)
      t16 = t15**2
      t18 = f1(l3,z,l1)
      tncoff10 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t12/t16*t18/512
      t10 = l2/2
      t11 = s23r2(z,l1,t10,l3,l4)
      t12 = t11**2
      t15 = s1vr2(z,l1,t10,l3,l4)
      t16 = t15**2
      t18 = f1(l3,z,l1)
      tncoff11 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t12/t16*t18/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t11 = dlog(-t1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t21 = dlog(1-l1)
      t23 = dlog(z)
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
      tncoff12 = -t1*(-6*t3-6*t8-24*t11+6*t14-18*t16-6*t18-6*t21+12*t23-
     #6*t25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/2
     #048
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t11 = dlog(-t1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t21 = dlog(1-l1)
      t23 = dlog(z)
      t25 = dlog(l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t3-6*t8-24*t11+6*t14-18*t16-6*t18-6*t21+12*t23-
     #6*t25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff15 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff16 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff17 = -(-1+z)/(l2-2)/l2/t9/t12*t15/1024
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff18 = -(-1+z)/(l2-2)/l2/t9/t12*t15/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff19 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff20 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqjrmnum1e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqjrmnum1e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqjrmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqjrmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqjrm1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm1e1int(var,ndim,Tcoeff)

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
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(l1)
      t12 = dlog(z)
      t15 = dlog(-t1)
      t17 = dlog(2.D0)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(t2+t4+t8-t10+t11-2*t12+4*t15+3*t17+t20)*t31/2048-t1/
     #l8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(l1)
      t12 = dlog(z)
      t15 = dlog(-t1)
      t17 = dlog(2.D0)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t2+t4+t8-t10+t11-2*t12+4*t15+3*t17+t20)*t31/2048+t1/
     #(-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff4 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff5 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff6 = -(-1+z)/(l2-2)/l2/t9/t12*t15/1024
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff7 = -(-1+z)/(l2-2)/l2/t9/t12*t15/1024
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqjrm1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm1e2int(var,ndim,Tcoeff)

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

      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqjrm1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

