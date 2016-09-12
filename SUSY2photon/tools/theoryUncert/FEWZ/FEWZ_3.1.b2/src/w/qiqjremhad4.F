CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqjrm4e0int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum3e0
      external qiqjrmnum3e0
      double precision qiqjrmnum3e1
      external qiqjrmnum3e1
      double precision qiqjrmnum3e2
      external qiqjrmnum3e2
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
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t4 = dlog(z)
      t8 = dlog(-t1)
      t12 = dlog(1-l1)
      t15 = dlog(l1)
      t18 = dlog(2.D0)
      t23 = dlog(l4)
      t27 = dlog(1-l4)
      t36 = 12*t3*t4-24*t8*t3+24*t8*t12+24*t8*t15+96*t8*t18+24*t15*t18-1
     #2*t23*t4+6*t27*t23+24*t12*t18-12*t12*t4+6*t27*t15
      t39 = f3(0.D0,z)
      t42 = dlog(1-l1*t39)
      t65 = 6*t15*t12+6*t42*t15+24*t42*t18-12*t15*t4-6*t42*t3+24*t42*t8-
     #48*t8*t4-6*t3*t15-6*t3*t12-48*t4*t18-24*t3*t18+6*t23*t42
      t89 = 6*t23*t12+24*t23*t8-6*t23*t3+6*t23*t15+24*t23*t18+6*t27*t12+
     #6*t27*t42+6*t42*t12-6*t27*t3+24*t27*t8-12*t27*t4
      t94 = 0.3141592653589793D1**2
      t95 = t18**2
      t97 = t42**2
      t99 = t23**2
      t101 = t27**2
      t103 = t12**2
      t105 = t4**2
      t107 = t8**2
      t109 = t15**2
      t111 = t3**2
      t113 = 24*t27*t18-12*t42*t4-t94+48*t95+3*t97+3*t99+3*t101+3*t103+1
     #2*t105+48*t107+3*t109+3*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 2*t4
      t129 = dlog(l8)
      t141 = 4*t18
      t142 = 4*t8
      t145 = -t15+t128-t12-t42+t129*l8+4*t8*l8-t3*l8-2*t4*l8+4*t18*l8+t1
     #2*l8+t15*l8+t42*l8-t129-t23+t3-t27-t141-t142+t23*l8+t27*l8
      t147 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t164 = t1*t160
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = -t1*(t36+t65+t89+t113)*t125/12288+t1*t145*t147/(-1+l8)*t
     #119*t122*t124/2048+t1*(t23+t142-t3-t128+t141+t12+t15+t42+2*t156+t2
     #7)*t160*t125/1024-t164*t147*t125/1024+t1*(t27+t142-t3-t128+t141+t1
     #2+t15+t42+t168+t23)*t171*t125/2048-t1*t171*t147*t125/2048-t164*t17
     #1*t125/1024
      t1 = -1+z
      t2 = dlog(2.D0)
      t3 = t2**2
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t9 = dlog(z)
      t12 = t8**2
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-l1)
      t20 = dlog(l4)
      t21 = t20**2
      t24 = dlog(1-l4)
      t25 = t24**2
      t27 = t17**2
      t29 = t9**2
      t32 = dlog(-t1)
      t33 = t32**2
      t35 = dlog(l1)
      t36 = t35**2
      t40 = -48*t3+12*t8*t9-3*t12+6*t15*t17-3*t21-3*t25-3*t27-12*t29-48*
     #t33-3*t36-6*t20*t17
      t41 = t15**2
      t65 = -3*t41-24*t8*t2-6*t24*t20-6*t20*t8-24*t35*t2-12*t15*t9-6*t24
     #*t8-24*t20*t32-6*t35*t17+12*t24*t9+48*t9*t2-24*t32*t35
      t89 = -96*t32*t2+48*t32*t9-24*t8*t32+6*t15*t35+6*t8*t15-24*t17*t2-
     #24*t24*t2-6*t8*t35+12*t20*t9-24*t32*t17+24*t32*t15
      t112 = 0.3141592653589793D1**2
      t113 = 24*t15*t2+6*t20*t15-24*t20*t2-6*t20*t35+12*t35*t9-6*t24*t17
     #+6*t24*t15-24*t24*t32-6*t8*t17-6*t24*t35+12*t17*t9+t112
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
      t148 = -1/t130
      t154 = dlog(l3)
      t156 = 4*t32
      t157 = 4*t2
      t158 = 2*t9
      t161 = 1/l3
      t165 = t1*t161
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t40+t65+t89+t113)*t125/12288+t1*(-t8*l8-t20*l8-t131*
     #l8-t24*l8+t15*l8+2*t9*l8-4*t2*l8-4*t32*l8-t17*l8-t35*l8)/l8*t148*t
     #119*t122*t124/2048-t1*(-t20-t24-2*t154-t8-t156+t15-t157-t17-t35+t1
     #58)*t161*t125/1024+t165*t148*t125/1024-t1*(-t24-t20-t8-t156+t15-t1
     #57-t17-t35+t158-t169)*t172*t125/2048+t1*t172*t148*t125/2048-t165*t
     #172*t125/1024
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(l4)
      t8 = dlog(1-l4)
      t9 = dlog(l8)
      t10 = 1-l8
      t11 = dlog(t10)
      t13 = dlog(-t1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(z)
      t19 = dlog(2.D0)
      t22 = dlog(1-l1)
      t23 = dlog(l1)
      t24 = t5+t6+t8+t9+t11+4*t13-t16-2*t17+4*t19+t22+t23
      t26 = 1/l8
      t29 = -1/t10
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t26*t39/2048-t1/l3*t29*t26*t32*t38/1024-t1/l2*t26
     #*t39/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l4)
      t6 = dlog(-t1)
      t8 = dlog(z)
      t10 = dlog(l1)
      t12 = dlog(-t2)
      t15 = dlog(1-l1)
      t16 = dlog(2.D0)
      t18 = dlog(l3)
      t20 = f3(l3,z)
      t23 = dlog(1-l1*t20)
      t25 = dlog(1-l4)
      t26 = f2(l3,z)
      t27 = dlog(t26)
      t28 = t4+4*t6-2*t8+t10+2*t12+t15+4*t16+2*t18+t23+t25-t27
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
      t4 = f3(l3,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(z)
      t10 = f2(l3,z)
      t11 = dlog(t10)
      t12 = dlog(l1)
      t13 = dlog(l4)
      t14 = dlog(l3)
      t17 = dlog(1-l4)
      t19 = dlog(1-l1)
      t21 = dlog(-t1)
      t24 = dlog(-t2)
      t26 = dlog(2.D0)
      t28 = -t7+2*t8+t11-t12-t13-2*t14-t17-t19-4*t21-2*t24-4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = -t3*t28*t40/1024+t3*t30/(-1+l8)*t33*t39/1024-t3/l2*t40/1
     #024
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/(-1+l8)/l8/t11/t15*t17/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(1-l1)
      t8 = dlog(-t1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l2)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t19 = dlog(z)
      t21 = dlog(2.D0)
      t23 = dlog(l2)
      t24 = t3+t4+t6+4*t8+t13+t15-t17+t18-2*t19+4*t21+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t39 = 1-l8
      tncoff7 = -t1*t24*t26*t30*t33*t35/2048+t1*t26/l8*t30*t33*t35/2048+
     #t30/l3*t1*t26*t33*t35/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(1-l1)
      t8 = dlog(-t1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l2)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t19 = dlog(z)
      t21 = dlog(2.D0)
      t23 = dlog(l2)
      t24 = t3+t4+t6+4*t8+t13+t15-t17+t18-2*t19+4*t21+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = -t1*t24*t26*t36/2048-t1*t26/(-1+l8)*t36/2048+t30/l3*t1*t
     #26*t33*t35/1024
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = (-1+z)/l2/l8/(-1+l8)/t9/t13*t15/2048
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = (l3-1)*(-1+z)/l2/l3/t8/t12*t14/1024
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = (l3-1)*(-1+z)/l2/l3/t8/t12*t14/1024
      t1 = -1+z
      t2 = dlog(z)
      t4 = dlog(2.D0)
      t6 = dlog(l4)
      t9 = dlog(1-l4)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t15 = dlog(1-l1)
      t18 = dlog(-t1)
      t20 = f3(0.D0,z)
      t23 = dlog(1-l1*t20)
      t25 = dlog(l1)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = -t1*(12*t2-24*t4-6*t6-6*t9+6*t12-6*t15-24*t18-6*t23-6*t
     #25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/2048
      t1 = -1+z
      t2 = dlog(l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t9 = dlog(l1)
      t12 = dlog(1-l4)
      t15 = dlog(-t1)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t21 = dlog(1-l1)
      t23 = dlog(z)
      t25 = dlog(2.D0)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = t1*(6*t2+6*t7+6*t9+6*t12+24*t15-6*t18+6*t21-12*t23+24*t
     #25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/l8/(-1+l8)/t8/t11*t14/2048
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
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff17 = (-1+z)/l2/t5/t8*t11/2048
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff18 = (-1+z)/l2/t5/t8*t11/2048
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

      tcnum1 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqjrmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqjrmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqjrmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqjrmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqjrm4e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm4e1int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum3e0
      external qiqjrmnum3e0
      double precision qiqjrmnum3e1
      external qiqjrmnum3e1
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
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t6 = dlog(1-l4)
      t7 = dlog(2.D0)
      t10 = dlog(-t1)
      t12 = dlog(l1)
      t13 = dlog(z)
      t16 = dlog(1-l1)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(t2-t4+t6+4*t7+4*t10+t12-2*t13+t16+t20)*t31/2048-t1/l
     #8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = dlog(-t1)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t10 = dlog(z)
      t12 = dlog(2.D0)
      t15 = dlog(1-l1)
      t16 = dlog(l1)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = -t1*(-t2-t4-4*t6+t9+2*t10-4*t12-t15-t16-t20)*t31/2048+t1
     #/(-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/l8/(-1+l8)/t8/t11*t14/2048
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
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff6 = (-1+z)/l2/t5/t8*t11/2048
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff7 = (-1+z)/l2/t5/t8*t11/2048
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

      tcnum1 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqjrmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqjrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqjrm4e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrm4e2int(var,ndim,Tcoeff)

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
      double precision qiqjrmnum3e0
      external qiqjrmnum3e0
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

      tcnum1 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqjrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqjrm4e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

