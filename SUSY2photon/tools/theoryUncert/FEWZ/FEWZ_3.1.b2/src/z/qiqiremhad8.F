CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm8e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum4e0
      external qiqirmnum4e0
      double precision qiqirmnum4e1
      external qiqirmnum4e1
      double precision qiqirmnum4e2
      external qiqirmnum4e2
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
      t4 = t3**2
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = t7**2
      t11 = dlog(1-l4)
      t12 = t11**2
      t14 = dlog(l1)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t20 = f3(0.D0,z)
      t23 = dlog(1-l1*t20)
      t24 = t23**2
      t26 = dlog(z)
      t27 = t26**2
      t29 = dlog(2.D0)
      t38 = 48*t4+3*t8+3*t12+3*t15+3*t18+3*t24+12*t27+96*t3*t29+24*t14*t
     #29+24*t11*t3-48*t26*t29
      t46 = dlog(1-l1)
      t49 = t29**2
      t51 = t46**2
      t65 = -6*t7*t14-24*t7*t29+6*t17*t14+6*t14*t46+48*t49+3*t51-48*t3*t
     #26+6*t17*t46+6*t11*t23+24*t23*t3-6*t23*t7+6*t23*t46
      t89 = -12*t11*t26+12*t7*t26+24*t46*t29-12*t46*t26-6*t7*t46+24*t3*t
     #14+6*t23*t14+24*t23*t29-6*t17*t7-24*t3*t7+24*t3*t46
      t112 = 0.3141592653589793D1**2
      t113 = -12*t23*t26+6*t11*t46+6*t11*t14+24*t17*t29-12*t17*t26+6*t17
     #*t23+24*t17*t3-6*t11*t7+24*t11*t29-12*t14*t26+6*t11*t17-t112
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t130 = dlog(l8)
      t140 = 4*t29
      t141 = 2*t26
      t142 = 4*t3
      t145 = -2*t26*l8+t130*l8+4*t3*l8-t7*l8+t14*l8+t46*l8+4*t29*l8+t17*
     #l8-t11-t14-t17-t140+t141-t23-t142-t130+t7-t46+t23*l8+t11*l8
      t147 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = -t1*(t38+t65+t89+t113)*t125/12288+t1*t145*t147/(-1+l8)*t
     #119*t122*t124/2048+t1*(-t7+t140+t46+t23+t14-t141+t142+2*t156+t17+t
     #11)*t160*t125/1024-t1*t160*t147*t125/1024+t1*(t46+t140-t141+t14+t2
     #3+t142-t7+t168+t11+t17)*t171*t125/2048-t1*t147*t171*t125/2048-t171
     #*t1*t119*t160*t122*t124/1024
      t1 = -1+z
      t2 = 0.3141592653589793D1**2
      t3 = dlog(l4)
      t5 = dlog(-t1)
      t8 = dlog(2.D0)
      t9 = t8**2
      t12 = dlog(1-l1)
      t13 = t12**2
      t15 = t5**2
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t21 = dlog(z)
      t24 = dlog(l1)
      t27 = f2(0.D0,z)
      t28 = dlog(t27)
      t29 = t28**2
      t32 = dlog(1-l4)
      t33 = t32**2
      t35 = t24**2
      t37 = t3**2
      t39 = -t2+24*t3*t5+48*t9+3*t13+48*t15-12*t20*t21+24*t24*t8+3*t29+3
     #*t33+3*t35+3*t37
      t40 = t20**2
      t42 = t21**2
      t64 = 3*t40+12*t42-12*t3*t21+24*t32*t8+6*t32*t3-12*t24*t21-6*t32*t
     #28+6*t3*t20+6*t32*t24+24*t3*t8+6*t32*t12+24*t5*t12
      t88 = 24*t20*t8-6*t3*t28-24*t5*t28+6*t20*t24-6*t28*t12+24*t5*t24-1
     #2*t12*t21+12*t28*t21+24*t12*t8+6*t20*t12-12*t32*t21
      t113 = -6*t20*t28+6*t32*t20+24*t20*t5+6*t3*t12-48*t5*t21+96*t5*t8+
     #6*t3*t24+6*t24*t12-6*t28*t24-24*t28*t8+24*t32*t5-48*t21*t8
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
      t148 = -1/t139
      t154 = 4*t8
      t155 = 2*t21
      t156 = 4*t5
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      t176 = t1*t172
      tncoff2 = -t1*(t39+t64+t88+t113)*t125/12288+t1*(-t32*l8-t20*l8-4*t
     #5*l8+t28*l8-t24*l8-t12*l8-4*t8*l8+2*t21*l8-t140*l8-t3*l8)/l8*t148*
     #t119*t122*t124/2048+t1*(-t28+t154+t12+t20+t24-t155+t156+2*t157+t3+
     #t32)*t161*t125/1024+t1*t161*t148*t125/1024+t1*(t12+t154-t155+t24+t
     #20+t156-t28+t169+t32+t3)*t172*t125/2048+t176*t148*t125/2048-t176*t
     #119*t161*t122*t124/1024
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = 1-l8
      t6 = dlog(t5)
      t7 = dlog(l8)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(1-l1)
      t12 = dlog(l1)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(-t1)
      t24 = t2+t4+t6+t7-t9+t11+t12+4*t13-2*t15+t20+4*t22
      t26 = 1/l8
      t29 = -1/t5
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t46 = t26*t32*t38
      tncoff3 = t1*t24*t26*t29*t32*t38/2048-t1/l3*t29*t46/1024-t1/l2*t29
     #*t46/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = f2(l3,z)
      t5 = dlog(t4)
      t6 = dlog(2.D0)
      t8 = dlog(l1)
      t10 = dlog(1-l1)
      t11 = f3(l3,z)
      t14 = dlog(1-l1*t11)
      t15 = dlog(l3)
      t17 = dlog(z)
      t19 = dlog(l4)
      t21 = dlog(-t1)
      t24 = dlog(1-l4)
      t26 = dlog(-t2)
      t28 = -t5+4*t6+t8+t10+t14+2*t15-2*t17+t19+4*t21+t24+2*t26
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
      t4 = f2(l3,z)
      t5 = dlog(t4)
      t6 = dlog(2.D0)
      t8 = dlog(l1)
      t10 = dlog(1-l1)
      t11 = f3(l3,z)
      t14 = dlog(1-l1*t11)
      t15 = dlog(l3)
      t17 = dlog(z)
      t19 = dlog(l4)
      t21 = dlog(-t1)
      t24 = dlog(1-l4)
      t26 = dlog(-t2)
      t28 = -t5+4*t6+t8+t10+t14+2*t15-2*t17+t19+4*t21+t24+2*t26
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
      tncoff6 = -(-1+z)*(l3-1)/l3/(-1+l8)/l8/t11/t15*t17/1024
      t1 = -1+z
      t3 = dlog(1-l2)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t9 = dlog(-t1)
      t11 = dlog(z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l1)
      t17 = dlog(l1)
      t18 = dlog(2.D0)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+t7+4*t9-2*t11-t14+t16+t17+4*t18+t20+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff7 = -t1*t24*t26*t36/2048+t1*t26/l8*t36/2048+t26*t33*t1*t30/l
     #3*t35/1024
      t1 = -1+z
      t3 = dlog(1-l2)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t9 = dlog(-t1)
      t11 = dlog(z)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l1)
      t17 = dlog(l1)
      t18 = dlog(2.D0)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+t7+4*t9-2*t11-t14+t16+t17+4*t18+t20+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = -t1*t24*t26*t36/2048-t1*t26/(-1+l8)*t36/2048+t26*t33*t1*
     #t30/l3*t35/1024
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = (-1+z)/l2/(-1+l8)/l8/t9/t13*t15/2048
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/1024
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(z)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t13 = dlog(-t1)
      t15 = dlog(l4)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t21 = dlog(1-l1)
      t23 = dlog(l1)
      t25 = dlog(2.D0)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = -t1*(-6*t3+12*t5-6*t10-24*t13-6*t15+6*t18-6*t21-6*t23-2
     #4*t25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/2
     #048
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(z)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t13 = dlog(-t1)
      t15 = dlog(l4)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t21 = dlog(1-l1)
      t23 = dlog(l1)
      t25 = dlog(2.D0)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t3+12*t5-6*t10-24*t13-6*t15+6*t18-6*t21-6*t23-2
     #4*t25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
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

      tcnum1 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum4e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum4e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum4e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum4e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm8e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm8e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum4e0
      external qiqirmnum4e0
      double precision qiqirmnum4e1
      external qiqirmnum4e1
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
      t5 = dlog(1-l1)
      t7 = dlog(1-l4)
      t8 = dlog(l1)
      t9 = dlog(l4)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(-t1)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(-t3+t5+t7+t8+t9+4*t10-2*t12+t17+4*t19)*t31/2048-t1/l
     #8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t5 = dlog(1-l1)
      t7 = dlog(1-l4)
      t8 = dlog(l1)
      t9 = dlog(l4)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(-t1)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(-t3+t5+t7+t8+t9+4*t10-2*t12+t17+4*t19)*t31/2048+t1/(
     #-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
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

      tcnum1 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum4e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum4e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm8e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm8e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum4e0
      external qiqirmnum4e0
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

      tcnum1 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum4e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm8e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

