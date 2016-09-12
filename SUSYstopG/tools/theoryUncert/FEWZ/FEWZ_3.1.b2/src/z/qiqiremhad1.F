CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm1e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
      double precision qiqirmnum1e2
      external qiqirmnum1e2
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
      t5 = dlog(1-l1)
      t8 = dlog(l4)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t15 = dlog(2.D0)
      t20 = f2(0.D0,z)
      t21 = dlog(t20)
      t24 = dlog(z)
      t37 = 24*t3*t5+6*t8*t12+18*t8*t15+18*t5*t15-18*t21*t15-12*t8*t24-1
     #2*t5*t24-24*t3*t21-6*t8*t21+72*t3*t15+6*t8*t5
      t42 = 0.3141592653589793D1**2
      t43 = t15**2
      t45 = dlog(l1)
      t48 = t5**2
      t50 = t3**2
      t52 = t21**2
      t55 = dlog(1-l4)
      t56 = t55**2
      t58 = t45**2
      t60 = t8**2
      t62 = t12**2
      t64 = -48*t3*t24+24*t8*t3-t42+27*t43+24*t3*t45+3*t48+48*t50+3*t52+
     #3*t56+3*t58+3*t60+3*t62
      t66 = t24**2
      t88 = 12*t66-12*t45*t24+6*t8*t45+6*t45*t5-12*t55*t24+18*t55*t15+18
     #*t45*t15-6*t21*t5+18*t12*t15+12*t21*t24-6*t12*t21
      t113 = 24*t55*t3+6*t55*t12+6*t55*t45-6*t21*t45-6*t55*t21+6*t12*t45
     #+6*t55*t5-36*t24*t15+6*t55*t8+6*t12*t5+24*t12*t3-12*t12*t24
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 2*t24
      t129 = dlog(l8)
      t139 = 4*t3
      t140 = 3*t15
      t145 = t21+t128-t55+t129*l8+t45*l8+t5*l8-2*t24*l8+t55*l8+t8*l8+t12
     #*l8-t21*l8-t45-t139-t12-t140-t8-t5-t129+3*t15*l8+4*t3*l8
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      t175 = t1*t171
      tncoff1 = -t1*(t37+t64+t88+t113)*t125/12288+t1*t145/(-1+l8)*t150*t
     #119*t122*t124/2048+t1*(t55+t8+t45+t140-t21+t5+t12+t139-t128+2*t156
     #)*t160*t125/1024-t1*t160*t150*t125/1024+(t8+t55-t128+t45+t12+t5-t2
     #1+t139+t140+t168)*t1*t171*t125/2048-t175*t150*t125/2048-t175*t160*
     #t125/1024
      t1 = -1+z
      t2 = dlog(2.D0)
      t3 = t2**2
      t6 = dlog(-t1)
      t7 = t6**2
      t10 = dlog(1-l1)
      t11 = t10**2
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t15 = t14**2
      t17 = dlog(l1)
      t21 = dlog(1-l4)
      t22 = t21**2
      t24 = t17**2
      t28 = dlog(l4)
      t29 = f3(0.D0,z)
      t32 = dlog(1-l1*t29)
      t35 = t28**2
      t37 = t32**2
      t39 = -27*t3-48*t7-3*t11-3*t15-18*t17*t2-3*t22-3*t24-72*t6*t2-6*t2
     #8*t32-3*t35-3*t37
      t42 = dlog(z)
      t43 = t42**2
      t65 = 6*t28*t14-12*t43-6*t28*t10-24*t28*t6+6*t14*t17+6*t21*t14-6*t
     #21*t10-6*t32*t17+36*t42*t2-6*t21*t28+24*t6*t14-6*t32*t10
      t89 = 6*t32*t14-18*t21*t2-24*t6*t17-24*t6*t10+12*t10*t42+48*t6*t42
     #+6*t14*t10+12*t28*t42+18*t14*t2-18*t28*t2+12*t17*t42
      t102 = 0.3141592653589793D1**2
      t113 = -24*t32*t6-18*t32*t2-18*t10*t2+12*t32*t42-6*t21*t32-6*t21*t
     #17+t102+12*t21*t42-24*t21*t6-6*t17*t10-6*t28*t17-12*t14*t42
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t136 = 1-l8
      t137 = dlog(t136)
      t146 = -1/t136
      t154 = 4*t6
      t155 = 2*t42
      t156 = 3*t2
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t39+t65+t89+t113)*t125/12288+t1*(t14*l8-t17*l8-t32*l
     #8-4*t6*l8-t28*l8-3*t2*l8-t137*l8-t21*l8-t10*l8+2*t42*l8)*t146/l8*t
     #119*t122*t124/2048-t1*(-t28-t21-t10-t154-t32+t155-t17+t14-t156-2*t
     #157)*t161*t125/1024+t1*t161*t146*t125/1024-(-t28-t21-t156-t17-t10-
     #t32+t14+t155-t154-t169)*t1*t172*t125/2048+t1*t146*t172*t125/2048-t
     #1*t172*t161*t125/1024
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(2.D0)
      t9 = dlog(-t1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = dlog(l1)
      t15 = dlog(1-l1)
      t16 = dlog(z)
      t19 = dlog(1-l4)
      t20 = dlog(l8)
      t21 = dlog(l4)
      t22 = 1-l8
      t23 = dlog(t22)
      t24 = t5+3*t6+4*t9-t12+t13+t15-2*t16+t19+t20+t21+t23
      t27 = -1/t22
      t29 = 1/l8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t27*t39/2048-t1/l3*t29*t27*t32*t38/1024-t1*t27/l2
     #*t39/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(1-l4)
      t9 = f2(l3,z)
      t10 = dlog(t9)
      t11 = dlog(2.D0)
      t14 = dlog(-t2)
      t16 = dlog(l1)
      t17 = f3(l3,z)
      t20 = dlog(1-l1*t17)
      t21 = dlog(l4)
      t23 = dlog(1-l1)
      t24 = dlog(z)
      t26 = dlog(l3)
      t28 = 4*t5+t8-t10+3*t11+2*t14+t16+t20+t21+t23-2*t24+2*t26
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
      tncoff4 = t3*t28*t40/1024-t3*t30/l8*t33*t37*t38/1024-t3/l2*t40/102
     #4
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(1-l4)
      t9 = f2(l3,z)
      t10 = dlog(t9)
      t11 = dlog(2.D0)
      t14 = dlog(-t2)
      t16 = dlog(l1)
      t17 = f3(l3,z)
      t20 = dlog(1-l1*t17)
      t21 = dlog(l4)
      t23 = dlog(1-l1)
      t24 = dlog(z)
      t26 = dlog(l3)
      t28 = 4*t5+t8-t10+3*t11+2*t14+t16+t20+t21+t23-2*t24+2*t26
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
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = -l2+2
      t10 = dlog(t9)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(l1)
      t17 = dlog(-t1)
      t20 = dlog(1-l1)
      t21 = dlog(z)
      t23 = dlog(l2)
      t24 = t2+t4+t8+t10-t12+2*t13+t15+4*t17+t20-2*t21+t23
      t27 = -1/t9
      t29 = 1/l2
      t30 = l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,t30,0.D0,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(0.D0,z,l1)
      t43 = 1-l8
      tncoff7 = t1*t24*t27*t29*t33*t37*t38/1024-t1*t27*t29/l8*t33*t37*t3
     #8/1024-t37*t27*t29*t1*t33/l3*t38/512
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = -l2+2
      t10 = dlog(t9)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(l1)
      t17 = dlog(-t1)
      t20 = dlog(1-l1)
      t21 = dlog(z)
      t23 = dlog(l2)
      t24 = t2+t4+t8+t10-t12+2*t13+t15+4*t17+t20-2*t21+t23
      t27 = -1/t9
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
      tncoff8 = t1*t24*t27*t29*t33*t39/1024+t1*t27*t29/(-1+l8)*t33*t39/1
     #024-t37*t27*t29*t1*t33/l3*t38/512
      t11 = l2/2
      t12 = s23r2(z,l1,t11,0.D0,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,0.D0,l4)
      t17 = t16**2
      t19 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/(l2-2)/l2/l8/(-1+l8)/t13/t17*t19/1024
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
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t6 = dlog(1-l1)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(l1)
      t17 = dlog(z)
      t19 = dlog(l4)
      t22 = dlog(1-l4)
      t25 = dlog(-t1)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = -t1*(6*t3-6*t6-18*t8-6*t13-6*t15+12*t17-6*t19-6*t22-24*
     #t25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/204
     #8
      t1 = -1+z
      t2 = dlog(2.D0)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t9 = dlog(z)
      t12 = dlog(-t1)
      t15 = dlog(1-l4)
      t17 = dlog(l4)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t23 = dlog(1-l1)
      t25 = dlog(l1)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = t1*(18*t2+6*t7-12*t9+24*t12+6*t15+6*t17-6*t20+6*t23+6*t
     #25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
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

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum1e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum1e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm1e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
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
      t3 = dlog(z)
      t5 = dlog(2.D0)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t12 = dlog(1-l4)
      t14 = dlog(-t1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(t2-2*t3+3*t5+t10+t12+4*t14-t17+t18+t20)*t31/2048-t1/
     #l8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t7 = dlog(l1)
      t9 = dlog(1-l1)
      t10 = dlog(l4)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = -t1*(-4*t3+t6-t7-t9-t10+2*t11-3*t13-t18-t20)*t31/2048+t1
     #/(-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
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

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm1e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
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

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

