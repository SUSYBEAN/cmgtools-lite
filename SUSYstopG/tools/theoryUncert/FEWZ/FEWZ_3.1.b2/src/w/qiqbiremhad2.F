CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqbirm2e0int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum2e0
      external qiqbirmnum2e0
      double precision qiqbirmnum2e1
      external qiqbirmnum2e1
      double precision qiqbirmnum2e2
      external qiqbirmnum2e2
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
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(l4)
      t10 = dlog(1-l1)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t23 = dlog(l1)
      t26 = dlog(z)
      t29 = dlog(2.D0)
      t32 = t26**2
      t34 = t8**2
      t37 = dlog(1-l4)
      t40 = 24*t3*t5-6*t8*t10-24*t8*t3-6*t8*t18+6*t8*t5-6*t8*t23+12*t8*t
     #26-24*t8*t29-12*t32-3*t34-24*t37*t3
      t65 = -24*t18*t29-6*t10*t23+12*t10*t26-24*t10*t29+6*t5*t23+6*t5*t1
     #0-12*t5*t26+12*t18*t26-24*t3*t10-24*t3*t23+48*t3*t26-6*t18*t23
      t81 = t5**2
      t83 = t10**2
      t85 = t23**2
      t87 = t3**2
      t89 = 6*t18*t5-6*t37*t8-6*t37*t18+6*t37*t5-6*t37*t10-6*t37*t23-24*
     #t37*t29-3*t81-3*t83-3*t85-48*t87
      t90 = t18**2
      t92 = t37**2
      t94 = 0.3141592653589793D1**2
      t95 = t29**2
      t113 = -3*t90-3*t92+t94-48*t95-96*t3*t29+12*t37*t26-6*t18*t10-24*t
     #23*t29+48*t26*t29+12*t23*t26+24*t5*t29-24*t3*t18
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 4*t29
      t129 = 4*t3
      t130 = 2*t26
      t131 = dlog(l8)
      t145 = -t37-t23-t18-t128-t129-t8-t10+t5+t130-t131-t5*l8-2*t26*l8+t
     #23*l8+t10*l8+t18*l8+t37*l8+t131*l8+t8*l8+4*t29*l8+4*t3*l8
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      tncoff1 = t1*(t40+t65+t89+t113)*t125/6144+t1*t145/(-1+l8)*t150*t11
     #9*t122*t124/1024-t1*(-t37-t129-t18+t5-t23-t10+t130-t128-2*t156-t8)
     #*t160*t125/512-t1*t150*t160*t125/512-t1*(-t168-t129-t18+t5-t23-t10
     #+t130-t128-t8-t37)*t171*t125/1024-t1*t171*t150*t125/1024-t1*t119*t
     #171*t160*t122*t124/512
      t1 = -1+z
      t2 = dlog(z)
      t3 = t2**2
      t6 = dlog(1-l1)
      t7 = dlog(l1)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t12 = dlog(2.D0)
      t16 = dlog(-t1)
      t17 = t16**2
      t19 = f3(0.D0,z)
      t22 = dlog(1-l1*t19)
      t23 = t22**2
      t26 = dlog(1-l4)
      t27 = t26**2
      t39 = -12*t3-6*t6*t7+24*t11*t12-48*t17-3*t23-3*t27+12*t6*t2-24*t6*
     #t12+6*t11*t7+6*t11*t6-12*t11*t2
      t52 = dlog(l4)
      t65 = 12*t22*t2-24*t16*t6-24*t16*t7+48*t16*t2-6*t22*t7+6*t22*t11-6
     #*t26*t52-6*t26*t22+6*t26*t11-6*t26*t6-6*t26*t7+12*t26*t2
      t87 = t52**2
      t89 = -24*t26*t12-24*t16*t22+24*t16*t11-6*t52*t6-24*t52*t16-6*t52*
     #t22+6*t52*t11-6*t52*t7+12*t52*t2-24*t52*t12-3*t87
      t90 = t11**2
      t92 = t6**2
      t94 = t7**2
      t110 = 0.3141592653589793D1**2
      t111 = t12**2
      t113 = -3*t90-3*t92-3*t94-24*t22*t12-24*t26*t16+12*t7*t2-6*t22*t6-
     #24*t7*t12+48*t2*t12-96*t16*t12+t110-48*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t129 = 1-l8
      t130 = dlog(t129)
      t146 = -1/t129
      t154 = 4*t16
      t155 = 2*t2
      t156 = 4*t12
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t39+t65+t89+t113)*t125/6144+t1*(-t26*l8-t130*l8-4*t1
     #6*l8-t22*l8+t11*l8-t6*l8-t7*l8+2*t2*l8-4*t12*l8-t52*l8)*t146/l8*t1
     #19*t122*t124/1024-t1*(-t26-t154-t22+t11-t7-t6+t155-t156-2*t157-t52
     #)*t161*t125/512+t1*t161*t146*t125/512-t1*(-t169-t154-t22+t11-t7-t6
     #+t155-t156-t52-t26)*t172*t125/1024+t1*t146*t172*t125/1024-t1*t119*
     #t172*t161*t122*t124/512
      t1 = -1+z
      t2 = dlog(2.D0)
      t4 = 1-l8
      t5 = dlog(t4)
      t6 = dlog(l4)
      t8 = dlog(1-l4)
      t9 = dlog(l8)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t13 = dlog(1-l1)
      t14 = dlog(l1)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t20 = dlog(-t1)
      t22 = dlog(z)
      t24 = 4*t2+t5+t6+t8+t9-t11+t13+t14+t18+4*t20-2*t22
      t27 = -1/t4
      t29 = 1/l8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t27*t39/1024-t1*t29/l3*t27*t32*t38/512-t1*t27/l2*
     #t39/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(1-l4)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t9 = dlog(-t1)
      t11 = dlog(l1)
      t12 = f3(l3,z)
      t15 = dlog(1-l1*t12)
      t17 = dlog(1-l1)
      t19 = dlog(-t2)
      t21 = dlog(l3)
      t23 = dlog(z)
      t25 = dlog(l4)
      t26 = dlog(2.D0)
      t28 = -t5+t7-4*t9-t11-t15-t17-2*t19-2*t21+2*t23-t25-4*t26
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff4 = -t3*t28*t40/512-t3/l8*t40/512-t3/l2*t40/512
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t7 = dlog(1-l4)
      t8 = dlog(2.D0)
      t11 = dlog(-t2)
      t14 = dlog(-t1)
      t16 = f3(l3,z)
      t19 = dlog(1-l1*t16)
      t21 = dlog(1-l1)
      t22 = dlog(l3)
      t24 = dlog(l1)
      t25 = dlog(l4)
      t26 = f2(l3,z)
      t27 = dlog(t26)
      t28 = -2*t4+t7+4*t8+2*t11+4*t14+t19+t21+2*t22+t24+t25-t27
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
      tncoff6 = -(-1+z)*(l3-1)/l8/l3/(-1+l8)/t11/t15*t17/512
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(1-l2)
      t5 = dlog(z)
      t7 = dlog(2.D0)
      t10 = dlog(1-l1)
      t12 = dlog(-t1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = dlog(l1)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = -t2-t4+2*t5-4*t7-t10-4*t12-t17+t19-t20-t21-t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff7 = t1*t24*t26*t36/1024+t1*t26/l8*t36/1024+t1*t33*t30*t26/l3
     #*t35/512
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(1-l2)
      t5 = dlog(z)
      t7 = dlog(2.D0)
      t10 = dlog(1-l1)
      t12 = dlog(-t1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = dlog(l1)
      t21 = dlog(l4)
      t23 = dlog(1-l4)
      t24 = -t2-t4+2*t5-4*t7-t10-4*t12-t17+t19-t20-t21-t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = t1*t24*t26*t36/1024-t1/(-1+l8)*t26*t36/1024+t1*t33*t30*t
     #26/l3*t35/512
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
      t22 = dlog(l4)
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
      tncoff12 = t1*(24*t3+6*t8-6*t11+6*t14+6*t16-12*t18+24*t20+6*t22+6*
     #t25)*t37/6144-t1/l8*t31*t34*t36/1024-t1/l3*t37/512-t1/l2*t37/1024
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
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = t1*(24*t3+6*t8-6*t11+6*t14+6*t16-12*t18+24*t20+6*t22+6*
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

      tcnum1 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqbirmnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqbirmnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqbirmnum2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqbirmnum2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqbirm2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum
     #4*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tnco
     #ff8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncof
     #f12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnc
     #off16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*t
     #ncoff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbirm2e1int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum2e0
      external qiqbirmnum2e0
      double precision qiqbirmnum2e1
      external qiqbirmnum2e1
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
      t6 = dlog(l1)
      t8 = dlog(-t1)
      t10 = dlog(z)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t15 = dlog(1-l1)
      t17 = dlog(1-l4)
      t18 = dlog(2.D0)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-t5-t6-4*t8+2*t10+t13-t15-t17-4*t18-t20)*t31/1024-t
     #1/l8*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(l1)
      t8 = dlog(-t1)
      t10 = dlog(z)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t15 = dlog(1-l1)
      t17 = dlog(1-l4)
      t18 = dlog(2.D0)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = -t1*(-t5-t6-4*t8+2*t10+t13-t15-t17-4*t18-t20)*t31/1024+t
     #1/(-1+l8)*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
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

      tcnum1 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqbirmnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqbirmnum2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqbirm2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum
     #4*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tnco
     #ff8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbirm2e2int(var,ndim,Tcoeff)

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
      double precision qiqbirmnum2e0
      external qiqbirmnum2e0
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

      tcnum1 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbirmnum2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqbirm2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

