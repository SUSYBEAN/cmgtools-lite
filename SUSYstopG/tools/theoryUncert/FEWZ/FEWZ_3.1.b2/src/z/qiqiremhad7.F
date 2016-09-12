CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm7e0int(var,ndim,Tcoeff)

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
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(z)
      t10 = dlog(1-l4)
      t12 = dlog(1-l1)
      t15 = dlog(l1)
      t18 = dlog(l4)
      t19 = dlog(2.D0)
      t27 = dlog(-t1)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t36 = 0.3141592653589793D1**2
      t37 = t19**2
      t39 = -12*t5*t6+6*t10*t12+6*t10*t15+24*t18*t19-12*t18*t6+6*t18*t5+
     #24*t18*t27-6*t10*t31+24*t10*t19-t36+48*t37
      t40 = t12**2
      t46 = t27**2
      t48 = t31**2
      t50 = t10**2
      t52 = t15**2
      t54 = t18**2
      t56 = t5**2
      t58 = t6**2
      t64 = 3*t40-12*t15*t6+6*t10*t18+48*t46+3*t48+3*t50+3*t52+3*t54+3*t
     #56+12*t58+24*t15*t19+24*t10*t27
      t88 = -48*t6*t19-6*t31*t15-24*t31*t19+6*t15*t12+6*t18*t15+96*t27*t
     #19-48*t27*t6+6*t18*t12+6*t10*t5+24*t5*t27+6*t5*t12
      t113 = -6*t5*t31-12*t10*t6+12*t31*t6+24*t12*t19-12*t12*t6-6*t31*t1
     #2+24*t27*t15+6*t5*t15+24*t5*t19-6*t18*t31-24*t27*t31+24*t27*t12
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 4*t27
      t129 = 2*t6
      t130 = dlog(l8)
      t131 = 4*t19
      t145 = t15+t10+t128-t31-t129+t5+t130+t18+t12+t131-t18*l8-4*t19*l8+
     #2*t6*l8-t10*l8-t130*l8-t15*l8-t12*l8-t5*l8-4*t27*l8+t31*l8
      t147 = -1+l8
      t148 = 1/t147
      t153 = 1/l8*t119*t122*t124
      t156 = dlog(l3)
      t160 = 1/l3
      t165 = -t1*t147
      t170 = dlog(l2)
      t173 = 1/l2
      tncoff1 = -t1*(t39+t64+t88+t113)*t125/12288-t1*t145*t148*t153/2048
     #+t1*(-t31+t131+t12+t5+t15-t129+t128+2*t156+t18+t10)*t160*t125/1024
     #+t165*t160*t148*t153/1024+t1*(t12+t131-t129+t15+t5+t128-t31+t170+t
     #10+t18)*t173*t125/2048+t165*t173*t148*t153/2048-t173*t1*t119*t160*
     #t122*t124/1024
      t1 = -1+z
      t3 = dlog(-t1)
      t4 = dlog(l1)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t13 = dlog(2.D0)
      t16 = dlog(l4)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t24 = dlog(1-l1)
      t27 = dlog(z)
      t31 = dlog(1-l4)
      t40 = 24*t3*t4+6*t10*t4+24*t10*t13-6*t16*t18-24*t3*t18+24*t3*t24-1
     #2*t10*t27+6*t31*t24+6*t31*t4+24*t16*t13-12*t16*t27
      t65 = 6*t16*t10+24*t16*t3-6*t18*t31+24*t31*t13-48*t27*t13-6*t18*t4
     #-24*t18*t13+6*t4*t24+6*t16*t4+96*t3*t13-48*t3*t27+6*t16*t24
      t83 = 0.3141592653589793D1**2
      t84 = t13**2
      t86 = t24**2
      t88 = 6*t10*t31+24*t10*t3-6*t10*t18+6*t10*t24-12*t27*t31+12*t18*t2
     #7-6*t18*t24+24*t4*t13-t83+48*t84+3*t86
      t89 = t3**2
      t91 = t18**2
      t97 = t31**2
      t99 = t4**2
      t103 = t16**2
      t105 = t10**2
      t107 = t27**2
      t113 = 48*t89+3*t91-12*t4*t27+6*t31*t16+3*t97+3*t99-12*t24*t27+3*t
     #103+3*t105+12*t107+24*t24*t13+24*t31*t3
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
      t154 = 4*t13
      t155 = 2*t27
      t156 = 4*t3
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      t176 = t1*t172
      tncoff2 = -t1*(t40+t65+t88+t113)*t125/12288-t1*(4*t3*l8+t131*l8+t1
     #6*l8+t10*l8-t18*l8+t4*l8+t31*l8+t24*l8+4*t13*l8-2*t27*l8)*t146/l8*
     #t119*t122*t124/2048+t1*(-t18+t154+t24+t10+t4-t155+t156+2*t157+t16+
     #t31)*t161*t125/1024+t1*t161*t146*t125/1024+t1*(t24+t154-t155+t4+t1
     #0+t156-t18+t169+t31+t16)*t172*t125/2048+t176*t146*t125/2048-t176*t
     #119*t161*t122*t124/1024
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t4 = dlog(z)
      t6 = dlog(2.D0)
      t9 = dlog(1-l1)
      t10 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t18 = dlog(-t1)
      t21 = dlog(1-l4)
      t22 = dlog(l8)
      t23 = dlog(l4)
      t24 = -t3+2*t4-4*t6-t9-t10+t12-t16-4*t18-t21-t22-t23
      t27 = -1/t2
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t39 = 1/l8/t31/t35*t37
      tncoff3 = -t1*t24*t27*t39/2048-t1/l3*t27*t39/1024-t1/l2*t27*t39/20
     #48
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
      t37 = 1/t36
      t38 = f1(l3,z,l1)
      t40 = t30*t33*t37*t38
      t42 = 1-l8
      tncoff4 = t3*t28*t40/1024-t3*t30/l8*t33*t37*t38/1024-t3/l2*t40/102
     #4
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
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff12 = -t1*(-6*t3+12*t5-6*t10-24*t13-6*t15+6*t18-6*t21-6*t23-2
     #4*t25)*t37/12288-t1/l8*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
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

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum3e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm7e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm7e1int(var,ndim,Tcoeff)

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
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t6 = dlog(1-l4)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t11 = f3(0.D0,z)
      t14 = dlog(1-l1*t11)
      t16 = dlog(-t1)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = dlog(l1)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t25 = 1/t24
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,l1)
      t31 = t25*t28*t30
      t34 = 1-l8
      tncoff1 = t1*(t2+t4+t6+4*t7-2*t9+t14+4*t16-t19+t20)*t31/2048-t1/l8
     #*t25*t28*t30/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t6 = dlog(1-l4)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t11 = f3(0.D0,z)
      t14 = dlog(1-l1*t11)
      t16 = dlog(-t1)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = dlog(l1)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t2+t4+t6+4*t7-2*t9+t14+4*t16-t19+t20)*t31/2048+t1/(-
     #1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
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

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum3e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm7e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm7e2int(var,ndim,Tcoeff)

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
      tncoff1 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm7e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

