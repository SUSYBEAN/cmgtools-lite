CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqbirm1e0int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum1e0
      external qiqbirmnum1e0
      double precision qiqbirmnum1e1
      external qiqbirmnum1e1
      double precision qiqbirmnum1e2
      external qiqbirmnum1e2
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
      t3 = dlog(1-l1)
      t4 = dlog(2.D0)
      t7 = dlog(z)
      t11 = dlog(-t1)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t20 = f2(0.D0,z)
      t21 = dlog(t20)
      t25 = dlog(1-l4)
      t30 = dlog(l1)
      t37 = dlog(l4)
      t40 = -24*t3*t4+12*t3*t7-24*t11*t15-24*t11*t3+24*t11*t21-24*t25*t4
     #+6*t21*t3+6*t21*t30+24*t21*t4-12*t21*t7-6*t37*t15
      t65 = -24*t37*t11+12*t25*t7+6*t37*t21-6*t37*t30-6*t37*t3-24*t37*t4
     #-6*t25*t37-6*t25*t3-6*t25*t15+6*t25*t21-6*t25*t30-24*t11*t30
      t87 = t15**2
      t89 = -24*t25*t11-96*t11*t4+48*t11*t7-6*t15*t3+6*t15*t21+12*t15*t7
     #-6*t15*t30-6*t3*t30-24*t15*t4+12*t37*t7-3*t87
      t96 = t25**2
      t98 = 0.3141592653589793D1**2
      t99 = t4**2
      t101 = t7**2
      t103 = t37**2
      t105 = t21**2
      t107 = t3**2
      t109 = t30**2
      t111 = t11**2
      t113 = 12*t30*t7-24*t30*t4+48*t7*t4-3*t96+t98-48*t99-12*t101-3*t10
     #3-3*t105-3*t107-3*t109-48*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t139 = dlog(l8)
      t142 = 4*t11
      t143 = 4*t4
      t144 = 2*t7
      t145 = t21*l8-4*t11*l8-t15*l8-t30*l8+2*t7*l8-4*t4*l8-t37*l8-t25*l8
     #-t139*l8-t3*l8+t3+t139+t25+t142+t37+t15+t30-t21+t143-t144
      t147 = 1/l8
      t149 = -1+l8
      t150 = 1/t149
      t152 = t122*t124
      t153 = t150*t119*t152
      t156 = dlog(l3)
      t160 = 1/l3
      t165 = -t1*t149
      t170 = dlog(l2)
      t173 = 1/l2
      tncoff1 = t1*(t40+t65+t89+t113)*t125/6144-t1*t145*t147*t153/1024+t
     #1*(t37+t25+2*t156+t142+t15-t21+t3+t30-t144+t143)*t160*t125/512+t16
     #5*t160*t147*t153/512+t1*(t25+t142+t15-t21+t3+t30-t144+t143+t170+t3
     #7)*t173*t125/1024+t165*t150*t173*t147*t119*t152/1024-t1*t119*t173*
     #t160*t122*t124/512
      t1 = -1+z
      t2 = dlog(z)
      t3 = dlog(2.D0)
      t6 = dlog(l4)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t18 = dlog(-t1)
      t22 = dlog(1-l1)
      t29 = 0.3141592653589793D1**2
      t31 = dlog(1-l4)
      t34 = dlog(l1)
      t39 = 48*t2*t3+12*t6*t2+6*t12*t14+24*t18*t14+12*t22*t2-24*t18*t22-
     #24*t18*t12+t29+6*t31*t14-24*t18*t34-96*t18*t3
      t64 = -6*t6*t34-6*t6*t22+12*t31*t2-6*t22*t34-24*t34*t3-24*t12*t3-6
     #*t12*t34-24*t6*t18-6*t6*t12+6*t14*t22+24*t14*t3+6*t6*t14
      t76 = t3**2
      t78 = t2**2
      t80 = t6**2
      t82 = t14**2
      t84 = t22**2
      t86 = t34**2
      t88 = -12*t14*t2-24*t31*t3+6*t14*t34-24*t22*t3+12*t34*t2-48*t76-12
     #*t78-3*t80-3*t82-3*t84-3*t86
      t89 = t18**2
      t93 = t31**2
      t95 = t12**2
      t113 = -48*t89-24*t6*t3-3*t93-3*t95-6*t31*t6-6*t31*t22-6*t31*t12-6
     #*t31*t34-6*t12*t22-24*t31*t18+48*t18*t2+12*t12*t2
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 1-l8
      t129 = dlog(t128)
      t148 = -1/t128
      t154 = dlog(l3)
      t156 = 4*t18
      t157 = 2*t2
      t158 = 4*t3
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t39+t64+t88+t113)*t125/6144-t1*(t129*l8+t34*l8+t22*l
     #8+t6*l8+t31*l8+4*t18*l8+t12*l8-t14*l8-2*t2*l8+4*t3*l8)/l8*t148*t11
     #9*t122*t124/1024+t1*(t6+t31+2*t154+t156+t12-t14+t22+t34-t157+t158)
     #*t161*t125/512+t1*t161*t148*t125/512+t1*(t31+t156+t12-t14+t22+t34-
     #t157+t158+t169+t6)*t172*t125/1024+t1*t148*t172*t125/1024-t1*t119*t
     #172*t161*t122*t124/512
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t4 = dlog(l8)
      t5 = dlog(l4)
      t7 = dlog(-t1)
      t10 = dlog(1-l4)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t17 = dlog(z)
      t19 = dlog(2.D0)
      t22 = dlog(1-l1)
      t23 = dlog(l1)
      t24 = -t3-t4-t5-4*t7-t10+t12-t16+2*t17-4*t19-t22-t23
      t26 = 1/l8
      t29 = -1/t2
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = -t1*t24*t26*t39/1024-t1/l3*t26*t39/512-t1*t29/l2*t26*t32
     #*t38/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l1)
      t5 = dlog(l4)
      t6 = f3(l3,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l4)
      t12 = dlog(l3)
      t14 = f2(l3,z)
      t15 = dlog(t14)
      t16 = dlog(z)
      t18 = dlog(2.D0)
      t21 = dlog(-t2)
      t24 = dlog(-t1)
      t27 = dlog(1-l1)
      t28 = t4+t5+t9+t11+2*t12-t15-2*t16+4*t18+2*t21+4*t24+t27
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
      t4 = dlog(l1)
      t5 = dlog(l4)
      t6 = f3(l3,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l4)
      t12 = dlog(l3)
      t14 = f2(l3,z)
      t15 = dlog(t14)
      t16 = dlog(z)
      t18 = dlog(2.D0)
      t21 = dlog(-t2)
      t24 = dlog(-t1)
      t27 = dlog(1-l1)
      t28 = t4+t5+t9+t11+2*t12-t15-2*t16+4*t18+2*t21+4*t24+t27
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
      tncoff6 = -(-1+z)*(l3-1)/l3/l8/(-1+l8)/t11/t15*t17/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(z)
      t10 = dlog(2.D0)
      t12 = dlog(l1)
      t14 = dlog(-t1)
      t17 = dlog(1-l2)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t21 = dlog(1-l1)
      t22 = dlog(l2)
      t23 = dlog(l4)
      t24 = t3+t7-2*t8+4*t10+t12+4*t14+t17-t19+t21+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff7 = -t1*t24*t26*t36/1024+t1*t26/l8*t36/1024+t1*t33*t30*t26/l
     #3*t35/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(z)
      t10 = dlog(2.D0)
      t12 = dlog(l1)
      t14 = dlog(-t1)
      t17 = dlog(1-l2)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t21 = dlog(1-l1)
      t22 = dlog(l2)
      t23 = dlog(l4)
      t24 = t3+t7-2*t8+4*t10+t12+4*t14+t17-t19+t21+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = -t1*t24*t26*t36/1024-t1/(-1+l8)*t26*t36/1024+t1*t33*t30*
     #t26/l3*t35/512
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = (-1+z)/(-1+l8)/l2/l8/t9/t13*t15/1024
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
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t14 = dlog(1-l1)
      t16 = dlog(l1)
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
      tncoff12 = t1*(24*t3+6*t8-6*t11+6*t14+6*t16-12*t18+24*t20+6*t23+6*
     #t25)*t37/6144-t1/l8*t37/1024-t1/l3*t37/512-t1/l2*t37/1024
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t14 = dlog(1-l1)
      t16 = dlog(l1)
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
      tncoff13 = t1*(24*t3+6*t8-6*t11+6*t14+6*t16-12*t18+24*t20+6*t23+6*
     #t25)*t37/6144+t1/(-1+l8)*t37/1024-t1/l3*t37/512-t1/l2*t37/1024
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/l8/(-1+l8)/t8/t11*t14/1024
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

      tcnum1 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqbirmnum1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqbirmnum1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqbirmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqbirmnum1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqbirm1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum
     #4*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tnco
     #ff8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncof
     #f12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnc
     #off16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*t
     #ncoff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbirm1e1int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum1e0
      external qiqbirmnum1e0
      double precision qiqbirmnum1e1
      external qiqbirmnum1e1
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
      t6 = dlog(-t1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t14 = dlog(l1)
      t15 = dlog(z)
      t17 = dlog(2.D0)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t25 = 1/t24
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,l1)
      t31 = t25*t28*t30
      t34 = 1-l8
      tncoff1 = t1*(t2+t4+4*t6+t11-t13+t14-2*t15+4*t17+t20)*t31/1024-t1/
     #l8*t25*t28*t30/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t6 = dlog(-t1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t14 = dlog(l1)
      t15 = dlog(z)
      t17 = dlog(2.D0)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t2+t4+4*t6+t11-t13+t14-2*t15+4*t17+t20)*t31/1024+t1/
     #(-1+l8)*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/l8/(-1+l8)/t8/t11*t14/1024
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

      tcnum1 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqbirmnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqbirmnum1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqbirm1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum
     #4*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tnco
     #ff8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbirm1e2int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum1e0
      external qiqbirmnum1e0
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

      tcnum1 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqbirm1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

