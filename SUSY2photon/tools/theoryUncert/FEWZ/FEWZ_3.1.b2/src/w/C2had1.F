CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C21e0int(var,ndim,Tcoeff)

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
      double precision C2num1e0
      external C2num1e0
      double precision C2num1e1
      external C2num1e1
      double precision C2num1e2
      external C2num1e2
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
      t10 = dlog(1-l1)
      t11 = t10**2
      t13 = dlog(z)
      t14 = t13**2
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t20 = t19**2
      t22 = dlog(l4)
      t23 = t22**2
      t26 = dlog(-t1)
      t27 = t26**2
      t29 = 0.3141592653589793D1**2
      t31 = dlog(1-l4)
      t32 = t31**2
      t36 = dlog(2.D0)
      t39 = -3*t3-3*t7-3*t11-12*t14-3*t20-3*t23-48*t27+t29-3*t32-24*t26*
     #t10+48*t13*t36
      t40 = t36**2
      t64 = -48*t40+12*t22*t13-6*t31*t10+6*t31*t6-24*t31*t26-6*t31*t19+6
     #*t19*t6-6*t19*t2+12*t10*t13-24*t31*t36+12*t19*t13-24*t22*t26
      t88 = -24*t26*t2-6*t22*t19+12*t31*t13-6*t31*t22+24*t26*t6-24*t26*t
     #19-6*t19*t10-6*t10*t2-24*t10*t36-12*t6*t13+6*t10*t6
      t113 = 6*t6*t2+24*t6*t36-24*t2*t36+12*t2*t13-96*t26*t36+48*t26*t13
     #-6*t31*t2-6*t22*t10-24*t19*t36-24*t22*t36+6*t22*t6-6*t22*t2
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t134 = dlog(l8)
      t142 = 4*t26
      t143 = 4*t36
      t144 = 2*t13
      t145 = t2*l8+4*t36*l8-2*t13*l8+t22*l8+t134*l8+t31*l8+t19*l8+4*t26*
     #l8+t10*l8-t6*l8-t31-t134+t6-t19-t10-t142-t2-t143+t144-t22
      t147 = 1/l8
      t149 = -1+l8
      t150 = 1/t149
      t152 = t122*t124
      t153 = t150*t119*t152
      t156 = dlog(l3)
      t160 = 1/l3
      t165 = -t1*t149
      t172 = dlog(l2)
      t175 = 1/l2
      tncoff1 = t1*(t39+t64+t88+t113)*t125/6144+t1*t145*t147*t153/1024+t
     #1*(t22+t31+t19+t10-t6+t2+t143-t144+t142+2*t156)*t160*t125/512+t165
     #*t160*t150*t147*t119*t152/512+t1*(t31+t19+t10-t6+t2+t143-t144+t142
     #+t172+t22)*t175*t125/1024+t165*t175*t147*t153/1024-t1*t119*t175*t1
     #60*t122*t124/512
      t1 = -1+z
      t3 = dlog(1-l1)
      t4 = t3**2
      t6 = dlog(z)
      t7 = t6**2
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = t12**2
      t16 = dlog(1-l4)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t22 = dlog(-t1)
      t29 = dlog(l1)
      t34 = dlog(2.D0)
      t39 = -3*t4-12*t7-3*t13+6*t16*t18-24*t16*t22-6*t16*t12+6*t12*t18-6
     #*t12*t29+12*t3*t6-96*t22*t34+48*t22*t6
      t42 = dlog(l4)
      t43 = t42**2
      t57 = t22**2
      t65 = -6*t16*t29-3*t43-6*t42*t3-24*t12*t34-24*t42*t34+6*t42*t18-6*
     #t42*t29-6*t16*t3-48*t57+24*t18*t34-24*t29*t34+12*t29*t6
      t85 = t16**2
      t89 = -24*t22*t12-6*t12*t3-6*t3*t29-24*t3*t34-12*t18*t6+6*t3*t18+6
     #*t18*t29-24*t42*t22-24*t22*t29-3*t85-6*t42*t12
      t106 = t34**2
      t108 = 0.3141592653589793D1**2
      t109 = t29**2
      t111 = t18**2
      t113 = 12*t16*t6-6*t16*t42+24*t22*t18+12*t42*t6+48*t6*t34-24*t22*t
     #3-24*t16*t34+12*t12*t6-48*t106+t108-3*t109-3*t111
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
      t154 = 4*t34
      t155 = 2*t6
      t156 = 4*t22
      t157 = dlog(l3)
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t39+t65+t89+t113)*t125/6144+t1*(-t129*l8+2*t6*l8-t42
     #*l8-4*t22*l8-t16*l8-t12*l8-t3*l8+t18*l8-t29*l8-4*t34*l8)/l8*t148*t
     #119*t122*t124/1024+t1*(t42+t16+t12+t3-t18+t29+t154-t155+t156+2*t15
     #7)*t161*t125/512+t1*t161*t148*t125/512+t1*(t16+t12+t3-t18+t29+t154
     #-t155+t156+t169+t42)*t172*t125/1024+t1*t172*t148*t125/1024-t1*t119
     #*t172*t161*t122*t124/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(z)
      t7 = dlog(l8)
      t8 = 1-l8
      t9 = dlog(t8)
      t11 = dlog(-t1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t18 = dlog(1-l1)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t21 = dlog(l1)
      t22 = dlog(2.D0)
      t24 = t3+t4-2*t5+t7+t9+4*t11+t16+t18-t20+t21+4*t22
      t26 = 1/l8
      t29 = -1/t8
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t26*t39/1024-t1/l3*t29*t26*t32*t38/512-t1/l2*t26*
     #t39/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(1-l4)
      t6 = dlog(l4)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(2.D0)
      t12 = dlog(1-l1)
      t14 = dlog(-t2)
      t16 = f3(l3,z)
      t19 = dlog(1-l1*t16)
      t20 = dlog(z)
      t22 = dlog(l1)
      t23 = dlog(l3)
      t26 = dlog(-t1)
      t28 = t5+t6-t8+4*t9+t12+2*t14+t19-2*t20+t22+2*t23+4*t26
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
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = f2(l3,z)
      t9 = dlog(t8)
      t11 = dlog(-t2)
      t13 = dlog(l4)
      t14 = f3(l3,z)
      t17 = dlog(1-l1*t14)
      t18 = dlog(l1)
      t19 = dlog(l3)
      t22 = dlog(-t1)
      t25 = dlog(1-l4)
      t27 = dlog(1-l1)
      t28 = -4*t4+2*t6+t9-2*t11-t13-t17-t18-2*t19-4*t22-t25-t27
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = -t3*t28*t40/512+t3*t30/(-1+l8)*t33*t39/512-t3/l2*t40/512
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/(-1+l8)/l8/t11/t15*t17/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(1-l2)
      t6 = dlog(2.D0)
      t8 = dlog(z)
      t10 = dlog(l1)
      t12 = dlog(-t1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(1-l1)
      t20 = f2(0.D0,z)
      t21 = dlog(t20)
      t22 = dlog(l2)
      t23 = dlog(l4)
      t24 = t3+t5+4*t6-2*t8+t10+4*t12+t17+t19-t21+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff7 = -t1*t24*t26*t36/1024+t1*t26/l8*t36/1024+t33*t1*t30*t26/l
     #3*t35/512
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(1-l2)
      t5 = dlog(l1)
      t6 = dlog(2.D0)
      t8 = dlog(z)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t13 = dlog(-t1)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t20 = dlog(1-l1)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = -t2-t4-t5-4*t6+2*t8+t11-4*t13-t18-t20-t22-t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = t1*t24*t26*t36/1024-t1*t26/(-1+l8)*t36/1024+t33*t1*t30*t
     #26/l3*t35/512
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
      tncoff10 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/512
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = (-1+z)*(l3-1)/l2/l3/t8/t12*t14/512
      t1 = -1+z
      t3 = dlog(1-l1)
      t5 = dlog(l1)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t12 = dlog(1-l4)
      t14 = dlog(l4)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t20 = dlog(-t1)
      t22 = f3(0.D0,z)
      t25 = dlog(1-l1*t22)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = t1*(6*t3+6*t5+24*t7-12*t9+6*t12+6*t14-6*t17+24*t20+6*t2
     #5)*t37/6144-t1/l8*t31*t34*t36/1024-t1/l3*t37/512-t1/l2*t37/1024
      t1 = -1+z
      t3 = dlog(1-l1)
      t5 = dlog(l1)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t12 = dlog(1-l4)
      t14 = dlog(l4)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t20 = dlog(-t1)
      t22 = f3(0.D0,z)
      t25 = dlog(1-l1*t22)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = t1*(6*t3+6*t5+24*t7-12*t9+6*t12+6*t14-6*t17+24*t20+6*t2
     #5)*t37/6144+t1/(-1+l8)*t37/1024-t1/l3*t37/512-t1/l2*t37/1024
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

      tcnum1 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = C2num1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = C2num1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = C2num1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = C2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = C2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = C2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = C2num1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = C2num1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = C2num1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = C2num1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = C2num1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = C2num1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = C2num1e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = C2num1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = C2num1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = C2num1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = C2num1e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      C21e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C21e1int(var,ndim,Tcoeff)

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
      double precision C2num1e0
      external C2num1e0
      double precision C2num1e1
      external C2num1e1
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
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = f3(0.D0,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t15 = dlog(l1)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t25 = 1/t24
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = f1(0.D0,z,l1)
      t31 = t25*t28*t30
      t34 = 1-l8
      tncoff1 = t1*(t3-t5+t9+t11+4*t13+t15+4*t16-2*t18+t20)*t31/1024-t1/
     #l8*t25*t28*t30/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = f3(0.D0,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t15 = dlog(l1)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t3-t5+t9+t11+4*t13+t15+4*t16-2*t18+t20)*t31/1024+t1/
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

      tcnum1 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = C2num1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = C2num1e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = C2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = C2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = C2num1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = C2num1e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      C21e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C21e2int(var,ndim,Tcoeff)

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
      double precision C2num1e0
      external C2num1e0
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

      tcnum1 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num1e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      C21e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

