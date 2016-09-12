CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function BC1e0int(var,ndim,Tcoeff)

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
      double precision BCnum1e0
      external BCnum1e0
      double precision BCnum1e1
      external BCnum1e1
      double precision BCnum1e2
      external BCnum1e2
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

      t1 = z-1
      t3 = dlog(-t1)
      t4 = t3**2
      t6 = dlog(z)
      t7 = t6**2
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = t12**2
      t16 = dlog(1-l1)
      t19 = dlog(l1)
      t22 = dlog(l4)
      t23 = t22**2
      t27 = t19**2
      t29 = f2(0.D0,z)
      t30 = dlog(t29)
      t33 = t30**2
      t37 = 48*t4+12*t7+3*t13-12*t16*t6-12*t19*t6+3*t23+24*t22*t3+3*t27-
     #6*t22*t30+3*t33+6*t19*t16
      t41 = dlog(1-l4)
      t42 = t41**2
      t44 = dlog(2.D0)
      t47 = t16**2
      t55 = 0.3141592653589793D1**2
      t56 = t44**2
      t64 = -6*t30*t19+3*t42+24*t22*t44+3*t47+6*t22*t19+6*t22*t12+6*t12*
     #t19-t55+48*t56-12*t22*t6-12*t41*t6-12*t12*t6
      t88 = 6*t12*t16+24*t12*t44+6*t22*t16+6*t41*t16+24*t41*t44+24*t3*t1
     #2-24*t3*t30+24*t3*t16+24*t3*t19+24*t19*t44+12*t30*t6
      t113 = -6*t30*t12-6*t30*t16-24*t30*t44+6*t41*t19+6*t41*t12-6*t41*t
     #30+24*t41*t3-48*t6*t44+6*t41*t22+96*t3*t44-48*t3*t6+24*t16*t44
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t126 = 4*t44
      t127 = 2*t6
      t128 = 4*t3
      t129 = 1-l2
      t130 = dlog(t129)
      t134 = -1/t129
      t138 = dlog(l8)
      t152 = -t22+t30-t12-t19-t16-t126+t127-t128-t41-t138+t22*l8+t41*l8+
     #4*t3*l8-t30*l8+t12*l8+t19*l8+t16*l8+4*t44*l8-2*t6*l8+t138*l8
      t154 = 1/l8
      t167 = dlog(l3)
      t171 = 1/l3
      tncoff1 = -t1*(t37+t64+t88+t113)*t123/6144-t1*(t22-t30+t12+t19+t16
     #+t126-t127+t128+t41+t130)*t134*t123/1024+t1*t152*t154/(-1+l8)*t118
     #*t120*t122/1024+t1*t134*t154*t123/1024+t1*(2*t167-t30+t12+t19+t16+
     #t126-t127+t128+t22+t41)*t171*t123/512+t1*t118*t134*t171*t120*t122/
     #512-t1*t171*t154*t123/512
      t1 = z-1
      t3 = dlog(1-l1)
      t4 = dlog(2.D0)
      t7 = dlog(z)
      t10 = dlog(l1)
      t13 = dlog(l4)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t20 = t7**2
      t23 = dlog(-t1)
      t24 = t23**2
      t28 = f3(0.D0,z)
      t31 = dlog(1-l1*t28)
      t32 = t31**2
      t34 = t13**2
      t36 = t17**2
      t38 = 24*t3*t4-12*t3*t7-12*t10*t7+6*t13*t10-6*t13*t17+12*t20+48*t2
     #4+6*t10*t3+3*t32+3*t34+3*t36
      t39 = t10**2
      t42 = dlog(1-l4)
      t43 = t42**2
      t45 = t3**2
      t51 = 0.3141592653589793D1**2
      t64 = 3*t39+3*t43+3*t45+24*t13*t4-6*t17*t10-t51+6*t31*t10-12*t13*t
     #7-12*t42*t7-12*t31*t7+24*t13*t23+6*t13*t31
      t66 = t4**2
      t88 = 48*t66+6*t31*t3+24*t31*t4+6*t13*t3+6*t42*t3+24*t42*t4+24*t23
     #*t31+24*t23*t3-24*t23*t17+24*t23*t10+24*t10*t4
      t113 = 12*t17*t7-6*t17*t31-6*t17*t3-24*t17*t4+6*t42*t10+6*t42*t31+
     #24*t42*t23-6*t42*t17-48*t7*t4+6*t42*t13+96*t23*t4-48*t23*t7
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t138 = 1-l8
      t139 = dlog(t138)
      t146 = -1/t138
      t152 = 4*t4
      t153 = 2*t7
      t154 = 4*t23
      t155 = 1-l2
      t156 = dlog(t155)
      t160 = -1/t155
      t168 = dlog(l3)
      t172 = 1/l3
      tncoff2 = -t1*(t38+t64+t88+t113)*t123/6144+t1*(-4*t23*l8+t17*l8-t3
     #1*l8-t10*l8-t3*l8-4*t4*l8+2*t7*l8-t13*l8-t42*l8-t139*l8)/l8*t146*t
     #118*t120*t122/1024-t1*(t13-t17+t31+t10+t3+t152-t153+t154+t42+t156)
     #*t160*t123/1024-t1*t160*t146*t123/1024+t1*(2*t168-t17+t31+t10+t3+t
     #152-t153+t154+t13+t42)*t172*t123/512+t1*t172*t146*t123/512+t1*t118
     #*t160*t172*t120*t122/512
      t1 = z-1
      t2 = 1-l8
      t3 = dlog(t2)
      t4 = dlog(l4)
      t6 = dlog(1-l4)
      t7 = dlog(l8)
      t8 = dlog(l1)
      t9 = dlog(2.D0)
      t12 = dlog(-t1)
      t14 = dlog(z)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t23 = dlog(1-l1)
      t24 = t3+t4+t6+t7+t8+4*t9+4*t12-2*t14-t17+t21+t23
      t26 = 1/l8
      t29 = -1/t2
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = 1/t30
      t33 = s1vr2(z,l1,1.D0,0.D0,l4)
      t35 = f1(0.D0,z,l1)
      t36 = 1/t33*t35
      t37 = t29*t31*t36
      tncoff3 = t1*t24*t26*t37/1024+t1/(l2-1)*t26*t37/1024-t1/l3*t29*t26
     #*t31*t36/512
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t7 = dlog(-t2)
      t9 = f3(l3,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(l3)
      t15 = dlog(2.D0)
      t18 = dlog(1-l1)
      t19 = dlog(l1)
      t20 = dlog(l4)
      t22 = dlog(-t1)
      t25 = dlog(1-l4)
      t26 = f2(l3,z)
      t27 = dlog(t26)
      t28 = -2*t4+2*t7+t12+2*t13+4*t15+t18+t19+t20+4*t22+t25-t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = f1(l3,z,l1)
      t37 = 1/t34*t36
      t38 = t30*t32*t37
      tncoff4 = t3*t28*t38/512+t3/(l2-1)*t38/512-t3*t30/l8*t32*t37/512
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(z)
      t7 = dlog(-t2)
      t9 = f3(l3,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(l3)
      t15 = dlog(2.D0)
      t18 = dlog(1-l1)
      t19 = dlog(l1)
      t20 = dlog(l4)
      t22 = dlog(-t1)
      t25 = dlog(1-l4)
      t26 = f2(l3,z)
      t27 = dlog(t26)
      t28 = -2*t4+2*t7+t12+2*t13+4*t15+t18+t19+t20+4*t22+t25-t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = f1(l3,z,l1)
      t37 = 1/t34*t36
      t38 = t30*t32*t37
      tncoff5 = t3*t28*t38/512+t3*t30/(-1+l8)*t32*t37/512+t3/(l2-1)*t38/
     #512
      t10 = s23r2(z,l1,1.D0,l3,l4)
      t13 = s1vr2(z,l1,1.D0,l3,l4)
      t15 = f1(l3,z,l1)
      tncoff6 = -(z-1)*(l3-1)/l3/(-1+l8)/l8/t10/t13*t15/512
      t1 = z-1
      t2 = dlog(l4)
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t7 = dlog(z)
      t9 = dlog(l2)
      t11 = dlog(-t1)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t19 = dlog(1-l1)
      t21 = dlog(1-l4)
      t22 = 1-l2
      t23 = dlog(t22)
      t24 = t2+t6-2*t7+t9+4*t11+4*t13-t16+t17+t19+t21+t23
      t27 = -1/t22
      t29 = s23r2(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t38 = 1-l8
      tncoff7 = t1*t24*t27*t30*t32*t34/1024-t1*t27/l8*t30*t32*t34/1024-t
     #1*t30*t27*t32/l3*t34/512
      t1 = z-1
      t2 = dlog(l4)
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t7 = dlog(z)
      t9 = dlog(l2)
      t11 = dlog(-t1)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t19 = dlog(1-l1)
      t21 = dlog(1-l4)
      t22 = 1-l2
      t23 = dlog(t22)
      t24 = t2+t6-2*t7+t9+4*t11+4*t13-t16+t17+t19+t21+t23
      t27 = -1/t22
      t29 = s23r2(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t35 = t30*t32*t34
      tncoff8 = t1*t24*t27*t35/1024+t1*t27/(-1+l8)*t35/1024-t1*t30*t27*t
     #32/l3*t34/512
      t9 = s23r2(z,l1,l2,0.D0,l4)
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t14 = f1(0.D0,z,l1)
      tncoff9 = -(z-1)/(l2-1)/l8/(-1+l8)/t9/t12*t14/1024
      t8 = s23r2(z,l1,l2,l3,l4)
      t11 = s1vr2(z,l1,l2,l3,l4)
      t13 = f1(l3,z,l1)
      tncoff10 = -(z-1)*(l3-1)/(l2-1)/l3/t8/t11*t13/512
      t8 = s23r2(z,l1,l2,l3,l4)
      t11 = s1vr2(z,l1,l2,l3,l4)
      t13 = f1(l3,z,l1)
      tncoff11 = -(z-1)*(l3-1)/(l2-1)/l3/t8/t11*t13/512
      t1 = z-1
      t2 = dlog(l4)
      t5 = dlog(-t1)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(l1)
      t12 = dlog(z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t20 = dlog(1-l1)
      t22 = dlog(2.D0)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t35 = t30*t32*t34
      t43 = 1-l8
      tncoff12 = -t1*(-6*t2-24*t5+6*t8-6*t10+12*t12-6*t17-6*t20-24*t22-6
     #*t25)*t35/6144+t1/(l2-1)*t35/1024-t1/l8*t30*t32*t34/1024-t1/l3*t35
     #/512
      t1 = z-1
      t2 = dlog(l4)
      t5 = dlog(-t1)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(l1)
      t12 = dlog(z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t20 = dlog(1-l1)
      t22 = dlog(2.D0)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t34 = f1(0.D0,z,l1)
      t35 = 1/t29/t31*t34
      tncoff13 = -t1*(-6*t2-24*t5+6*t8-6*t10+12*t12-6*t17-6*t20-24*t22-6
     #*t25)*t35/6144+t1/(-1+l8)*t35/1024+t1/(l2-1)*t35/1024-t1/l3*t35/51
     #2
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t9 = s1vr2(z,l1,1.D0,0.D0,l4)
      t12 = f1(0.D0,z,l1)
      tncoff14 = -(z-1)/l8/(-1+l8)/t7/t9*t12/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t8 = s1vr2(z,l1,1.D0,l3,l4)
      t11 = f1(l3,z,l1)
      tncoff15 = -(z-1)*(l3-1)/l3/t6/t8*t11/512
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t8 = s1vr2(z,l1,1.D0,l3,l4)
      t11 = f1(l3,z,l1)
      tncoff16 = -(z-1)*(l3-1)/l3/t6/t8*t11/512
      t5 = s23r2(z,l1,l2,0.D0,l4)
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff17 = -(z-1)/(l2-1)/t5/t7*t10/1024
      t5 = s23r2(z,l1,l2,0.D0,l4)
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff18 = -(z-1)/(l2-1)/t5/t7*t10/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff19 = -(z-1)/t2/t5*t7/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff20 = -(z-1)/t2/t5*t7/1024

CCC Numerator functions

      tcnum1 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = BCnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = BCnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = BCnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = BCnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = BCnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = BCnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = BCnum1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = BCnum1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = BCnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = BCnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = BCnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = BCnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = BCnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = BCnum1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = BCnum1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = BCnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = BCnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      BC1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function BC1e1int(var,ndim,Tcoeff)

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
      double precision BCnum1e0
      external BCnum1e0
      double precision BCnum1e1
      external BCnum1e1
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

      t1 = z-1
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t7 = dlog(1-l4)
      t9 = dlog(1-l1)
      t10 = dlog(l4)
      t11 = dlog(z)
      t14 = dlog(-t1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t19 = dlog(2.D0)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,l1)
      t29 = 1/t23/t25*t28
      tncoff1 = t1*(t5+t7+t9+t10-2*t11+4*t14-t17+t18+4*t19)*t29/1024+t1/
     #(l2-1)*t29/1024-t1/l8*t29/1024-t1/l3*t29/512
      t1 = z-1
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t7 = dlog(1-l4)
      t9 = dlog(1-l1)
      t10 = dlog(l4)
      t11 = dlog(z)
      t14 = dlog(-t1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t19 = dlog(2.D0)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,l1)
      t29 = 1/t23/t25*t28
      tncoff2 = t1*(t5+t7+t9+t10-2*t11+4*t14-t17+t18+4*t19)*t29/1024+t1/
     #(-1+l8)*t29/1024+t1/(l2-1)*t29/1024-t1/l3*t29/512
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t9 = s1vr2(z,l1,1.D0,0.D0,l4)
      t12 = f1(0.D0,z,l1)
      tncoff3 = -(z-1)/l8/(-1+l8)/t7/t9*t12/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t8 = s1vr2(z,l1,1.D0,l3,l4)
      t11 = f1(l3,z,l1)
      tncoff4 = -(z-1)*(l3-1)/l3/t6/t8*t11/512
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t8 = s1vr2(z,l1,1.D0,l3,l4)
      t11 = f1(l3,z,l1)
      tncoff5 = -(z-1)*(l3-1)/l3/t6/t8*t11/512
      t5 = s23r2(z,l1,l2,0.D0,l4)
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff6 = -(z-1)/(l2-1)/t5/t7*t10/1024
      t5 = s23r2(z,l1,l2,0.D0,l4)
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff7 = -(z-1)/(l2-1)/t5/t7*t10/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff8 = -(z-1)/t2/t5*t7/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff9 = -(z-1)/t2/t5*t7/1024

CCC Numerator functions

      tcnum1 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = BCnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = BCnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = BCnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = BCnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = BCnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = BCnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      BC1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function BC1e2int(var,ndim,Tcoeff)

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
      double precision BCnum1e0
      external BCnum1e0
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

      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff1 = -(z-1)/t2/t5*t7/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t5 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,l1)
      tncoff2 = -(z-1)/t2/t5*t7/1024

CCC Numerator functions

      tcnum1 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      BC1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

