CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C22e0int(var,ndim,Tcoeff)

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
      double precision C2num2e0
      external C2num2e0
      double precision C2num2e1
      external C2num2e1
      double precision C2num2e2
      external C2num2e2
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
      t3 = dlog(2.D0)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t18 = dlog(l1)
      t25 = dlog(z)
      t29 = dlog(1-l4)
      t36 = -24*t2*t3+6*t2*t7-24*t11*t13+6*t11*t7-6*t11*t18-6*t2*t18-24*
     #t11*t3+12*t2*t25-24*t29*t13+12*t18*t25-24*t18*t3
      t45 = f3(0.D0,z)
      t48 = dlog(1-l1*t45)
      t65 = 48*t25*t3-24*t2*t13-6*t2*t11-6*t29*t18-6*t2*t48-24*t29*t3+6*
     #t48*t7-24*t48*t3+12*t48*t25-6*t18*t48-12*t7*t25+6*t18*t7
      t89 = 24*t7*t3-6*t29*t48-6*t29*t11-24*t13*t48+12*t29*t25+24*t13*t7
     #-24*t13*t18+48*t13*t25-96*t13*t3-6*t11*t48+6*t29*t7
      t94 = 0.3141592653589793D1**2
      t95 = t3**2
      t97 = t2**2
      t99 = t29**2
      t101 = t11**2
      t103 = t7**2
      t105 = t13**2
      t107 = t18**2
      t109 = t48**2
      t111 = t25**2
      t113 = 12*t11*t25-6*t29*t2+t94-48*t95-3*t97-3*t99-3*t101-3*t103-48
     #*t105-3*t107-3*t109-12*t111
      t117 = s23r2(z,l1,0.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,0.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t130 = 4*t13
      t133 = 2*t25
      t134 = dlog(l8)
      t144 = 4*t3
      t145 = 4*t13*l8-t130-t11+t11*l8+t2*l8-t48+t133+t7-t134+t134*l8+t18
     #*l8-2*t25*l8+4*t3*l8-t7*l8+t48*l8+t29*l8-t18-t29-t144-t2
      t147 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      t175 = t1*t171
      tncoff1 = t1*(t36+t65+t89+t113)*t125/6144+t1*t145*t147/(-1+l8)*t11
     #9*t122*t124/1024-t1*(-t29-t2-2*t156-t18+t7-t48-t11+t133-t130-t144)
     #*t160*t125/512-t1*t160*t147*t125/512+t1*(t18-t7+t48+t11-t133+t130+
     #t144+t168+t29+t2)*t171*t125/1024-t175*t147*t125/1024-t175*t160*t12
     #5/512
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l1)
      t8 = dlog(1-l4)
      t9 = dlog(2.D0)
      t13 = dlog(-t1)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t22 = 0.3141592653589793D1**2
      t23 = t9**2
      t25 = t2**2
      t27 = t8**2
      t29 = dlog(z)
      t32 = t4**2
      t34 = f2(0.D0,z)
      t35 = dlog(t34)
      t36 = t35**2
      t38 = -6*t2*t4-24*t8*t9-24*t2*t13-24*t19*t9+t22-48*t23-3*t25-3*t27
     #+12*t19*t29-3*t32-3*t36
      t39 = dlog(l1)
      t42 = t13**2
      t44 = t39**2
      t50 = t19**2
      t62 = t29**2
      t64 = 6*t39*t35-48*t42-3*t44+24*t35*t9+12*t29*t2-3*t50+12*t8*t29+2
     #4*t13*t35+6*t4*t35-6*t8*t19-6*t2*t39-12*t62
      t88 = -24*t13*t19-6*t8*t39-24*t39*t9+6*t8*t35-6*t8*t4-24*t13*t39-6
     #*t39*t19-12*t35*t29-24*t4*t13+48*t29*t13+48*t29*t9
      t113 = 12*t4*t29+12*t39*t29-6*t2*t19-6*t4*t39-24*t4*t9-24*t2*t9-24
     #*t8*t13+6*t2*t35-6*t8*t2-96*t13*t9+6*t19*t35-6*t4*t19
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
      t156 = 2*t29
      t157 = 4*t13
      t158 = 4*t9
      t161 = 1/l3
      t165 = t1*t148
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = t1*(t38+t64+t88+t113)*t125/6144+t1*(-t19*l8-4*t13*l8-t4*
     #l8-t39*l8-4*t9*l8+2*t29*l8-t138*l8-t2*l8+t35*l8-t8*l8)/l8*t148*t11
     #9*t122*t124/1024-t1*(-t8-t2-2*t154-t39+t35-t19-t4+t156-t157-t158)*
     #t161*t125/512+t165*t161*t125/512+t1*(t39-t35+t19+t4-t156+t157+t158
     #+t169+t8+t2)*t172*t125/1024+t165*t172*t125/1024-t1*t172*t161*t125/
     #512
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
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      t39 = 1-l8
      tncoff7 = -t1*t24*t26*t36/1024+t1*t26/l8*t30*t33*t35/1024+t1*t26/l
     #3*t36/512
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

      tcnum1 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = C2num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = C2num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = C2num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = C2num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = C2num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = C2num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = C2num2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999.D0
     #,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = C2num2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999.D0
     #,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = C2num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = C2num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = C2num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = C2num2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = C2num2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = C2num2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = C2num2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = C2num2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = C2num2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      C22e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C22e1int(var,ndim,Tcoeff)

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
      double precision C2num2e0
      external C2num2e0
      double precision C2num2e1
      external C2num2e1
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
      t4 = dlog(l4)
      t6 = dlog(1-l1)
      t8 = dlog(-t1)
      t10 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-t3-t4-t6-4*t8-t10+t12-4*t13+2*t15-t20)*t31/1024-t1
     #/l8*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(1-l1)
      t8 = dlog(-t1)
      t10 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = -t1*(-t3-t4-t6-4*t8-t10+t12-4*t13+2*t15-t20)*t31/1024+t1
     #/(-1+l8)*t31/1024-t1/l3*t31/512-t1/l2*t31/1024
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

      tcnum1 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = C2num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = C2num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = C2num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = C2num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = C2num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = C2num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      C22e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C22e2int(var,ndim,Tcoeff)

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
      double precision C2num2e0
      external C2num2e0
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

      tcnum1 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = C2num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      C22e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

