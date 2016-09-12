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
      t3 = dlog(1-l4)
      t5 = dlog(-t1)
      t8 = 0.3141592653589793D1**2
      t9 = dlog(2.D0)
      t10 = t9**2
      t12 = dlog(l4)
      t13 = t12**2
      t15 = t3**2
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t22 = dlog(1-l1)
      t23 = f3(0.D0,z)
      t26 = dlog(1-l1*t23)
      t31 = dlog(z)
      t34 = dlog(l1)
      t39 = -24*t3*t5+t8-48*t10-3*t13-3*t15+6*t3*t18-6*t22*t26-96*t5*t9+
     #48*t5*t31-24*t5*t34+24*t5*t18
      t48 = t22**2
      t50 = t18**2
      t52 = t5**2
      t54 = t34**2
      t64 = 12*t3*t31-24*t5*t26-6*t3*t22-6*t3*t26-3*t48-3*t50-48*t52-3*t
     #54+24*t18*t9+6*t34*t18-12*t18*t31-6*t34*t26
      t88 = 12*t26*t31-24*t26*t9+6*t26*t18-24*t3*t9-6*t12*t26-6*t3*t34-6
     #*t12*t22-24*t12*t5+48*t31*t9-24*t34*t9+12*t34*t31
      t91 = t26**2
      t95 = t31**2
      t113 = 12*t22*t31-3*t91+6*t12*t18-12*t95-6*t3*t12-24*t12*t9-24*t22
     #*t9-6*t12*t34-6*t22*t34+6*t22*t18-24*t22*t5+12*t12*t31
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t130 = 4*t5
      t133 = 2*t31
      t134 = dlog(l8)
      t144 = 4*t9
      t145 = 4*t5*l8-t130-t22+t22*l8+t12*l8-t26+t133+t18-t134+t134*l8+t3
     #4*l8-2*t31*l8+4*t9*l8-t18*l8+t26*l8+t3*l8-t34-t3-t144-t12
      t147 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      t175 = t1*t171
      tncoff1 = t1*(t39+t64+t88+t113)*t125/6144+t1*t145*t147/(-1+l8)*t11
     #9*t122*t124/1024-t1*(-t3-t12-2*t156-t34+t18-t26-t22+t133-t130-t144
     #)*t160*t125/512-t1*t160*t147*t125/512+t1*(t34-t18+t26+t22-t133+t13
     #0+t144+t168+t3+t12)*t171*t125/1024-t175*t147*t125/1024-t175*t160*t
     #125/512
      t1 = -1+z
      t2 = dlog(z)
      t3 = t2**2
      t5 = 0.3141592653589793D1**2
      t6 = dlog(2.D0)
      t7 = t6**2
      t10 = dlog(-t1)
      t11 = dlog(l1)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t19 = dlog(1-l4)
      t22 = f3(0.D0,z)
      t25 = dlog(1-l1*t22)
      t29 = dlog(1-l1)
      t36 = dlog(l4)
      t37 = t36**2
      t39 = -12*t3+t5-48*t7-24*t10*t11+24*t10*t15+12*t19*t2-24*t10*t25-6
     #*t19*t29-6*t19*t25+48*t10*t2-3*t37
      t40 = t19**2
      t42 = t29**2
      t48 = t15**2
      t50 = t10**2
      t52 = t11**2
      t54 = t25**2
      t64 = -3*t40-3*t42-24*t29*t10+6*t36*t15-3*t48-48*t50-3*t52-3*t54+2
     #4*t15*t6+6*t11*t15-12*t15*t2-6*t11*t25
      t88 = 12*t25*t2-24*t25*t6+6*t25*t15-24*t19*t6-6*t36*t25-6*t19*t11-
     #6*t36*t29-24*t36*t10+48*t2*t6-24*t11*t6+12*t11*t2
      t113 = -24*t19*t10+12*t36*t2-6*t19*t36+12*t29*t2+6*t19*t15-6*t29*t
     #25-96*t10*t6+6*t29*t15-24*t29*t6-24*t36*t6-6*t36*t11-6*t29*t11
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t137 = 1-l8
      t138 = dlog(t137)
      t148 = -1/t137
      t154 = dlog(l3)
      t156 = 2*t2
      t157 = 4*t10
      t158 = 4*t6
      t161 = 1/l3
      t165 = t1*t148
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t39+t64+t88+t113)*t125/6144+t1*(-t25*l8-4*t10*l8-t29
     #*l8-t11*l8-4*t6*l8+2*t2*l8-t138*l8-t36*l8+t15*l8-t19*l8)/l8*t148*t
     #119*t122*t124/1024-t1*(-t19-t36-2*t154-t11+t15-t25-t29+t156-t157-t
     #158)*t161*t125/512+t165*t161*t125/512+t1*(t11-t15+t25+t29-t156+t15
     #7+t158+t169+t19+t36)*t172*t125/1024+t165*t172*t125/1024-t1*t172*t1
     #61*t125/512
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t4 = dlog(2.D0)
      t7 = dlog(1-l4)
      t8 = dlog(l4)
      t9 = 1-l8
      t10 = dlog(t9)
      t11 = dlog(z)
      t14 = dlog(-t1)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t20 = dlog(l1)
      t21 = dlog(l8)
      t23 = dlog(1-l1)
      t24 = -t3+4*t4+t7+t8+t10-2*t11+4*t14+t19+t20+t21+t23
      t26 = 1/l8
      t29 = -1/t9
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t42 = t1*t29
      t46 = t26*t32*t38
      tncoff3 = t1*t24*t26*t29*t32*t38/1024-t42/l3*t46/512-t42/l2*t46/10
     #24
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(1-l4)
      t9 = dlog(l4)
      t10 = f2(l3,z)
      t11 = dlog(t10)
      t12 = dlog(l3)
      t14 = dlog(2.D0)
      t16 = f3(l3,z)
      t19 = dlog(1-l1*t16)
      t21 = dlog(-t2)
      t23 = dlog(z)
      t26 = dlog(1-l1)
      t27 = dlog(l1)
      t28 = 4*t5+t8+t9-t11+2*t12+4*t14+t19+2*t21-2*t23+t26+t27
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
      t5 = dlog(-t1)
      t8 = dlog(1-l4)
      t9 = dlog(l4)
      t10 = f2(l3,z)
      t11 = dlog(t10)
      t12 = dlog(l3)
      t14 = dlog(2.D0)
      t16 = f3(l3,z)
      t19 = dlog(1-l1*t16)
      t21 = dlog(-t2)
      t23 = dlog(z)
      t26 = dlog(1-l1)
      t27 = dlog(l1)
      t28 = 4*t5+t8+t9-t11+2*t12+4*t14+t19+2*t21-2*t23+t26+t27
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff5 = t3*t28*t40/512+t3/(-1+l8)*t40/512-t3/l2*t40/512
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/(-1+l8)/l3/l8/t11/t15*t17/512
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = dlog(2.D0)
      t7 = dlog(z)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(1-l1)
      t13 = dlog(l1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(1-l2)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = 4*t3+4*t5-2*t7-t10+t12+t13+t17+t19+t20+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t35 = f1(0.D0,z,l1)
      t36 = 1/t29/t32*t35
      t39 = t1*t26
      tncoff7 = -t1*t24*t26*t36/1024+t39/l8*t36/1024+t39/l3*t36/512
      t1 = -1+z
      t3 = dlog(-t1)
      t5 = dlog(2.D0)
      t7 = dlog(z)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(1-l1)
      t13 = dlog(l1)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = dlog(1-l2)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = 4*t3+4*t5-2*t7-t10+t12+t13+t17+t19+t20+t22+t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t35 = f1(0.D0,z,l1)
      t36 = 1/t29/t32*t35
      tncoff8 = -t1*t24*t26*t36/1024-t1/(-1+l8)*t26*t36/1024+t1*t26/l3*t
     #36/512
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
      t3 = dlog(1-l1)
      t6 = dlog(-t1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t13 = dlog(l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(2.D0)
      t20 = dlog(z)
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
      tncoff12 = t1*(6*t3+24*t6+6*t11+6*t13-6*t16+24*t18-12*t20+6*t22+6*
     #t25)*t37/6144-t1/l8*t31*t34*t36/1024-t1/l3*t37/512-t1/l2*t37/1024
      t1 = -1+z
      t3 = dlog(1-l1)
      t6 = dlog(-t1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t13 = dlog(l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(2.D0)
      t20 = dlog(z)
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = t1*(6*t3+24*t6+6*t11+6*t13-6*t16+24*t18-12*t20+6*t22+6*
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
      t3 = dlog(1-l1)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(-t1)
      t13 = dlog(z)
      t15 = dlog(l1)
      t17 = dlog(1-l4)
      t18 = dlog(2.D0)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-t3-t7+t9-4*t11+2*t13-t15-t17-4*t18-t20)*t31/1024-t
     #1/l8*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t3 = dlog(1-l1)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(-t1)
      t13 = dlog(z)
      t15 = dlog(l1)
      t17 = dlog(1-l4)
      t18 = dlog(2.D0)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = -t1*(-t3-t7+t9-4*t11+2*t13-t15-t17-4*t18-t20)*t31/1024+t
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

