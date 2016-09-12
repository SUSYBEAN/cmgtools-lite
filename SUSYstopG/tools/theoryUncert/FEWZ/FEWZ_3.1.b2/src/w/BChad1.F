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
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t4 = t3**2
      t6 = dlog(z)
      t7 = dlog(2.D0)
      t11 = dlog(-t1)
      t12 = dlog(l1)
      t23 = dlog(l4)
      t27 = dlog(1-l4)
      t32 = t27**2
      t34 = -3*t4+48*t6*t7-24*t11*t12-96*t11*t7-12*t3*t6+6*t3*t12+24*t3*
     #t7+6*t23*t3+12*t27*t6-6*t27*t23-3*t32
      t36 = dlog(1-l1)
      t43 = 0.3141592653589793D1**2
      t44 = t7**2
      t48 = f3(0.D0,z)
      t51 = dlog(1-l1*t48)
      t64 = -6*t36*t12+12*t12*t6+6*t36*t3+t43-48*t44+6*t27*t3+6*t51*t3-6
     #*t51*t12-24*t11*t36-24*t51*t7+12*t51*t6-6*t51*t36
      t86 = t51**2
      t88 = 12*t23*t6-24*t12*t7-6*t23*t51-6*t23*t36-24*t11*t51+12*t36*t6
     #+24*t11*t3-24*t23*t7-24*t23*t11-6*t27*t36-3*t86
      t89 = t12**2
      t91 = t36**2
      t93 = t23**2
      t95 = t11**2
      t97 = t6**2
      t113 = -3*t89-3*t91-3*t93-48*t95-12*t97-24*t27*t11-6*t27*t12-24*t2
     #7*t7-6*t27*t51-24*t36*t7-6*t23*t12+48*t11*t6
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t126 = 4*t7
      t127 = 2*t6
      t128 = 4*t11
      t129 = 1-l2
      t130 = dlog(t129)
      t134 = -1/t129
      t138 = dlog(l8)
      t152 = -t138+t23*l8-t27-t23-t51-t36+t3-t12-t126+t127-t128+t27*l8+t
     #138*l8+4*t11*l8+t51*l8+t36*l8-t3*l8+t12*l8+4*t7*l8-2*t6*l8
      t154 = 1/l8
      t163 = t1*t134
      t167 = dlog(l3)
      t171 = 1/l3
      tncoff1 = t1*(t34+t64+t88+t113)*t123/6144-t1*(t27+t23+t51+t36-t3+t
     #12+t126-t127+t128+t130)*t134*t123/1024+t1*t152*t154/(-1+l8)*t118*t
     #120*t122/1024+t163*t154*t123/1024+t1*(t27+t51+t36-t3+t12+t126-t127
     #+t128+2*t167+t23)*t171*t123/512+t163*t171*t123/512-t1*t171*t154*t1
     #23/512
      t1 = z-1
      t3 = dlog(-t1)
      t5 = dlog(1-l1)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t12 = dlog(2.D0)
      t15 = dlog(z)
      t20 = dlog(l4)
      t23 = dlog(l1)
      t34 = f2(0.D0,z)
      t35 = dlog(t34)
      t38 = -24*t3*t5-24*t11*t12+12*t11*t15-6*t11*t5+12*t20*t15-24*t23*t
     #12-6*t20*t11-6*t20*t5-24*t11*t3+12*t5*t15+24*t3*t35
      t44 = dlog(1-l4)
      t47 = t11**2
      t49 = t5**2
      t51 = t23**2
      t53 = t20**2
      t55 = t3**2
      t65 = -24*t20*t12-24*t20*t3-6*t44*t5-3*t47-3*t49-3*t51-3*t53-48*t5
     #5-6*t5*t23+12*t23*t15+6*t5*t35+6*t44*t35
      t67 = t15**2
      t69 = t35**2
      t89 = -12*t67-3*t69-6*t23*t11-24*t44*t3-6*t44*t23-24*t44*t12-6*t44
     #*t11-24*t5*t12-6*t20*t23+48*t3*t15+48*t15*t12
      t104 = 0.3141592653589793D1**2
      t105 = t12**2
      t111 = t44**2
      t113 = -24*t23*t3-96*t3*t12-12*t35*t15+6*t35*t23+24*t35*t12+6*t20*
     #t35-6*t44*t20+t104-48*t105+12*t44*t15+6*t11*t35-3*t111
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t138 = 1-l8
      t139 = dlog(t138)
      t146 = -1/t138
      t152 = 4*t12
      t153 = 2*t15
      t154 = 4*t3
      t155 = 1-l2
      t156 = dlog(t155)
      t160 = -1/t155
      t164 = t1*t160
      t168 = dlog(l3)
      t172 = 1/l3
      tncoff2 = t1*(t38+t65+t89+t113)*t123/6144+t1*(-t20*l8-t44*l8-4*t3*
     #l8-t11*l8-t5*l8+t35*l8-t23*l8-4*t12*l8+2*t15*l8-t139*l8)/l8*t146*t
     #118*t120*t122/1024-t1*(t44+t20+t11+t5-t35+t23+t152-t153+t154+t156)
     #*t160*t123/1024-t164*t146*t123/1024+t1*(t44+t11+t5-t35+t23+t152-t1
     #53+t154+2*t168+t20)*t172*t123/512+t1*t146*t172*t123/512+t164*t172*
     #t123/512
      t1 = z-1
      t2 = dlog(l4)
      t3 = 1-l8
      t4 = dlog(t3)
      t6 = dlog(1-l4)
      t7 = dlog(l8)
      t8 = f3(0.D0,z)
      t11 = dlog(1-l1*t8)
      t13 = dlog(1-l1)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t16 = dlog(l1)
      t17 = dlog(2.D0)
      t19 = dlog(z)
      t22 = dlog(-t1)
      t24 = t2+t4+t6+t7+t11+t13-t15+t16+4*t17-2*t19+4*t22
      t26 = 1/l8
      t29 = -1/t3
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = 1/t30
      t33 = s1vr2(z,l1,1.D0,0.D0,l4)
      t35 = f1(0.D0,z,l1)
      t36 = 1/t33*t35
      t45 = t26*t31*t36
      tncoff3 = t1*t24*t26*t29*t31*t36/1024+t1/(l2-1)*t29*t45/1024-t1*t2
     #9/l3*t45/512
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t12 = dlog(z)
      t15 = dlog(1-l4)
      t17 = dlog(-t2)
      t19 = dlog(l3)
      t21 = f3(l3,z)
      t24 = dlog(1-l1*t21)
      t25 = dlog(l4)
      t27 = dlog(1-l1)
      t28 = -4*t5+t8-4*t9-t11+2*t12-t15-2*t17-2*t19-t24-t25-t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = f1(l3,z,l1)
      t37 = 1/t34*t36
      t38 = t30*t32*t37
      tncoff4 = -t3*t28*t38/512+t3/(l2-1)*t38/512-t3*t30/l8*t32*t37/512
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t12 = dlog(z)
      t15 = dlog(1-l4)
      t17 = dlog(-t2)
      t19 = dlog(l3)
      t21 = f3(l3,z)
      t24 = dlog(1-l1*t21)
      t25 = dlog(l4)
      t27 = dlog(1-l1)
      t28 = -4*t5+t8-4*t9-t11+2*t12-t15-2*t17-2*t19-t24-t25-t27
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t34 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = f1(l3,z,l1)
      t38 = 1/l3/t31/t34*t36
      tncoff5 = -t3*t28*t38/512+t3/(-1+l8)*t38/512+t3/(l2-1)*t38/512
      t10 = s23r2(z,l1,1.D0,l3,l4)
      t13 = s1vr2(z,l1,1.D0,l3,l4)
      t15 = f1(l3,z,l1)
      tncoff6 = -(z-1)*(l3-1)/(-1+l8)/l3/l8/t10/t13*t15/512
      t1 = z-1
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(1-l1)
      t7 = dlog(z)
      t9 = dlog(l2)
      t11 = dlog(-t1)
      t13 = dlog(2.D0)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t21 = dlog(l1)
      t22 = 1-l2
      t23 = dlog(t22)
      t24 = t3+t4+t6-2*t7+t9+4*t11+4*t13+t18-t20+t21+t23
      t27 = -1/t22
      t29 = s23r2(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t35 = t30*t32*t34
      t38 = 1-l8
      tncoff7 = t1*t24*t27*t35/1024-t1*t27/l8*t30*t32*t34/1024-t1*t27/l3
     #*t35/512
      t1 = z-1
      t3 = dlog(1-l1)
      t4 = dlog(z)
      t6 = dlog(l2)
      t8 = dlog(-t1)
      t10 = dlog(2.D0)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(l1)
      t19 = 1-l2
      t20 = dlog(t19)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = -t3+2*t4-t6-4*t8-4*t10-t15+t17-t18-t20-t22-t23
      t27 = -1/t19
      t29 = s23r2(z,l1,l2,0.D0,l4)
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t34 = f1(0.D0,z,l1)
      t35 = 1/t29/t31*t34
      t38 = t1*t27
      tncoff8 = -t1*t24*t27*t35/1024+t38/(-1+l8)*t35/1024-t38/l3*t35/512
      t9 = s23r2(z,l1,l2,0.D0,l4)
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t14 = f1(0.D0,z,l1)
      tncoff9 = -(z-1)/(l2-1)/(-1+l8)/l8/t9/t12*t14/1024
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
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t14 = dlog(-t1)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t21 = dlog(2.D0)
      t23 = dlog(z)
      t25 = dlog(l1)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t35 = t30*t32*t34
      t43 = 1-l8
      tncoff12 = t1*(6*t2-6*t5+6*t8+6*t11+24*t14+6*t19+24*t21-12*t23+6*t
     #25)*t35/6144+t1/(l2-1)*t35/1024-t1/l8*t30*t32*t34/1024-t1/l3*t35/5
     #12
      t1 = z-1
      t2 = dlog(l4)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t14 = dlog(-t1)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t21 = dlog(2.D0)
      t23 = dlog(z)
      t25 = dlog(l1)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t34 = f1(0.D0,z,l1)
      t35 = 1/t29/t31*t34
      tncoff13 = t1*(6*t2-6*t5+6*t8+6*t11+24*t14+6*t19+24*t21-12*t23+6*t
     #25)*t35/6144+t1/(-1+l8)*t35/1024+t1/(l2-1)*t35/1024-t1/l3*t35/512
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
      t2 = dlog(2.D0)
      t4 = dlog(z)
      t7 = dlog(-t1)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(1-l4)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t17 = dlog(l4)
      t19 = dlog(1-l1)
      t20 = dlog(l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,l1)
      t29 = 1/t23/t25*t28
      tncoff1 = t1*(4*t2-2*t4+4*t7-t10+t12+t16+t17+t19+t20)*t29/1024+t1/
     #(l2-1)*t29/1024-t1/l8*t29/1024-t1/l3*t29/512
      t1 = z-1
      t2 = dlog(2.D0)
      t4 = dlog(z)
      t7 = dlog(-t1)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(1-l4)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t17 = dlog(l4)
      t19 = dlog(1-l1)
      t20 = dlog(l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,l1)
      t29 = 1/t23/t25*t28
      tncoff2 = t1*(4*t2-2*t4+4*t7-t10+t12+t16+t17+t19+t20)*t29/1024+t1/
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

