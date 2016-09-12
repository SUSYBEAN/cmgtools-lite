CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function BC2e0int(var,ndim,Tcoeff)

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
      double precision BCnum2e0
      external BCnum2e0
      double precision BCnum2e1
      external BCnum2e1
      double precision BCnum2e2
      external BCnum2e2
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
      t2 = 0.3141592653589793D1**2
      t4 = dlog(1-l4)
      t6 = dlog(-t1)
      t9 = dlog(2.D0)
      t10 = t9**2
      t14 = dlog(z)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t21 = dlog(l4)
      t28 = dlog(l1)
      t31 = f3(0.D0,z)
      t34 = dlog(1-l1*t31)
      t37 = -t2+24*t4*t6+48*t10+96*t6*t9-48*t6*t14-6*t4*t18+6*t4*t21+24*
     #t21*t6-6*t21*t18+6*t21*t28+6*t21*t34
      t53 = dlog(1-l1)
      t64 = 6*t4*t34-12*t4*t14-24*t18*t9+6*t4*t28+12*t18*t14+24*t28*t9-4
     #8*t14*t9+6*t4*t53+24*t4*t9-12*t34*t14+24*t53*t9-12*t53*t14
      t88 = 6*t28*t53+6*t21*t53-12*t28*t14-6*t18*t28-6*t18*t34-6*t18*t53
     #+24*t34*t9-24*t6*t18+24*t6*t28+24*t6*t53+24*t21*t9
      t97 = t14**2
      t99 = t6**2
      t101 = t34**2
      t103 = t21**2
      t105 = t18**2
      t107 = t28**2
      t109 = t4**2
      t111 = t53**2
      t113 = -12*t21*t14+24*t6*t34+6*t34*t53+6*t34*t28+12*t97+48*t99+3*t
     #101+3*t103+3*t105+3*t107+3*t109+3*t111
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t126 = 1-l2
      t127 = dlog(t126)
      t128 = 4*t6
      t129 = 4*t9
      t130 = 2*t14
      t134 = -1/t126
      t144 = dlog(l8)
      t152 = t28+t53-t28*l8-t53*l8-4*t9*l8+2*t14*l8-t144*l8-4*t6*l8+t18*
     #l8-t34*l8-t4*l8-t21*l8+t144+t21+t4-t18-t130+t129+t34+t128
      t154 = -1+l8
      t155 = 1/t154
      t157 = 1/l8
      t159 = t120*t122
      t160 = t157*t118*t159
      t164 = -t1*t154
      t169 = dlog(l3)
      t173 = 1/l3
      tncoff1 = -t1*(t37+t64+t88+t113)*t123/6144-t1*(t127+t128-t18+t34+t
     #28+t53+t129-t130+t4+t21)*t134*t123/1024-t1*t152*t155*t160/1024-t16
     #4*t134*t155*t160/1024+t1*(2*t169+t128-t18+t34+t28+t53+t129-t130+t2
     #1+t4)*t173*t123/512+t1*t118*t134*t173*t120*t122/512+t164*t157*t173
     #*t155*t118*t159/512
      t1 = z-1
      t2 = dlog(l4)
      t4 = dlog(-t1)
      t8 = dlog(1-l4)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t15 = dlog(z)
      t18 = dlog(2.D0)
      t23 = dlog(l1)
      t29 = dlog(1-l1)
      t36 = 24*t2*t4+6*t8*t2-6*t8*t12-48*t4*t15+96*t4*t18+24*t8*t4+6*t2*
     #t23-6*t2*t12+24*t4*t29+24*t4*t23-24*t4*t12
      t37 = f3(0.D0,z)
      t40 = dlog(1-l1*t37)
      t65 = 24*t40*t18-6*t12*t29-6*t12*t40-6*t12*t23-12*t23*t15+6*t2*t29
     #+6*t23*t29-12*t29*t15+24*t29*t18-12*t40*t15+24*t8*t18+6*t8*t29
      t89 = 24*t23*t18-48*t15*t18+12*t12*t15+6*t8*t23-24*t12*t18-12*t8*t
     #15+6*t8*t40+6*t2*t40+24*t4*t40+6*t40*t29+6*t40*t23
      t94 = 0.3141592653589793D1**2
      t95 = t18**2
      t97 = t15**2
      t99 = t4**2
      t101 = t40**2
      t103 = t2**2
      t105 = t12**2
      t107 = t23**2
      t109 = t8**2
      t111 = t29**2
      t113 = -12*t2*t15+24*t2*t18-t94+48*t95+12*t97+48*t99+3*t101+3*t103
     #+3*t105+3*t107+3*t109+3*t111
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = 1/t117
      t119 = s1vr2(z,l1,1.D0,0.D0,l4)
      t120 = 1/t119
      t122 = f1(0.D0,z,l1)
      t123 = t118*t120*t122
      t138 = 1-l8
      t139 = dlog(t138)
      t144 = -1/t138
      t152 = 1-l2
      t153 = dlog(t152)
      t154 = 4*t4
      t155 = 4*t18
      t156 = 2*t15
      t160 = -1/t152
      t168 = dlog(l3)
      t172 = 1/l3
      tncoff2 = -t1*(t36+t65+t89+t113)*t123/6144-t1*(t40*l8+t2*l8+4*t4*l
     #8-t12*l8+t8*l8+4*t18*l8+t23*l8+t29*l8-2*t15*l8+t139*l8)*t144/l8*t1
     #18*t120*t122/1024-t1*(t153+t154-t12+t40+t23+t29+t155-t156+t8+t2)*t
     #160*t123/1024-t1*t160*t144*t123/1024+t1*(2*t168+t154-t12+t40+t23+t
     #29+t155-t156+t2+t8)*t172*t123/512+t1*t172*t144*t123/512+t1*t118*t1
     #60*t172*t120*t122/512
      t1 = z-1
      t2 = dlog(l8)
      t3 = 1-l8
      t4 = dlog(t3)
      t5 = dlog(l1)
      t6 = dlog(l4)
      t8 = dlog(1-l4)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t16 = dlog(1-l1)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(-t1)
      t24 = -t2-t4-t5-t6-t8+t10+2*t11-4*t13-t16-t20-4*t22
      t27 = -1/t3
      t29 = 1/l8
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = 1/t30
      t33 = s1vr2(z,l1,1.D0,0.D0,l4)
      t35 = f1(0.D0,z,l1)
      t36 = 1/t33*t35
      t37 = t29*t31*t36
      tncoff3 = -t1*t24*t27*t37/1024+t1/(l2-1)*t27*t37/1024-t1*t29/l3*t2
     #7*t31*t36/512
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(-t2)
      t10 = f3(l3,z)
      t13 = dlog(1-l1*t10)
      t14 = dlog(z)
      t16 = dlog(l1)
      t17 = dlog(l3)
      t20 = dlog(1-l1)
      t21 = dlog(l4)
      t22 = f2(l3,z)
      t23 = dlog(t22)
      t25 = dlog(1-l4)
      t26 = dlog(2.D0)
      t28 = 4*t5+2*t8+t13-2*t14+t16+2*t17+t20+t21-t23+t25+4*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,l3,l4)
      t35 = 1/t34
      t36 = f1(l3,z,l1)
      t38 = t30*t32*t35*t36
      t44 = 1-l8
      tncoff4 = t3*t28*t38/512+t3/(l2-1)*t38/512-t3/l8*t30*t32*t35*t36/5
     #12
      t1 = z-1
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(-t2)
      t10 = f3(l3,z)
      t13 = dlog(1-l1*t10)
      t14 = dlog(z)
      t16 = dlog(l1)
      t17 = dlog(l3)
      t20 = dlog(1-l1)
      t21 = dlog(l4)
      t22 = f2(l3,z)
      t23 = dlog(t22)
      t25 = dlog(1-l4)
      t26 = dlog(2.D0)
      t28 = 4*t5+2*t8+t13-2*t14+t16+2*t17+t20+t21-t23+t25+4*t26
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
      tncoff6 = -(z-1)*(l3-1)/l8/l3/(-1+l8)/t10/t13*t15/512
      t1 = z-1
      t2 = 1-l2
      t3 = dlog(t2)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t10 = dlog(l2)
      t12 = dlog(-t1)
      t14 = dlog(z)
      t16 = dlog(l1)
      t18 = dlog(1-l1)
      t19 = dlog(2.D0)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+t7-t9+t10+4*t12-2*t14+t16+t18+4*t19+t22+t23
      t27 = -1/t2
      t29 = s23r2(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t38 = 1-l8
      tncoff7 = t1*t24*t27*t30*t32*t34/1024-t1*t27/l8*t30*t32*t34/1024-t
     #1*t30*t27*t32/l3*t34/512
      t1 = z-1
      t2 = 1-l2
      t3 = dlog(t2)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t10 = dlog(l2)
      t12 = dlog(-t1)
      t14 = dlog(z)
      t16 = dlog(l1)
      t18 = dlog(1-l1)
      t19 = dlog(2.D0)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+t7-t9+t10+4*t12-2*t14+t16+t18+4*t19+t22+t23
      t27 = -1/t2
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
      t2 = dlog(l1)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t10 = dlog(1-l1)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = dlog(-t1)
      t20 = dlog(1-l4)
      t22 = dlog(l4)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t34 = f1(0.D0,z,l1)
      t35 = 1/t29/t31*t34
      tncoff12 = -t1*(-6*t2-6*t7-6*t10-24*t12+12*t14-24*t17-6*t20-6*t22+
     #6*t25)*t35/6144+t1/(l2-1)*t35/1024-t1/l8*t35/1024-t1/l3*t35/512
      t1 = z-1
      t2 = dlog(l1)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t10 = dlog(1-l1)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = dlog(-t1)
      t20 = dlog(1-l4)
      t22 = dlog(l4)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t34 = f1(0.D0,z,l1)
      t35 = 1/t29/t31*t34
      tncoff13 = -t1*(-6*t2-6*t7-6*t10-24*t12+12*t14-24*t17-6*t20-6*t22+
     #6*t25)*t35/6144+t1/(-1+l8)*t35/1024+t1/(l2-1)*t35/1024-t1/l3*t35/5
     #12
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

      tcnum1 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = BCnum2e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = BCnum2e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = BCnum2e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = BCnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = BCnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = BCnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = BCnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = BCnum2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = BCnum2e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = BCnum2e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = BCnum2e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = BCnum2e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = BCnum2e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = BCnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = BCnum2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = BCnum2e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = BCnum2e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      BC2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function BC2e1int(var,ndim,Tcoeff)

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
      double precision BCnum2e0
      external BCnum2e0
      double precision BCnum2e1
      external BCnum2e1
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
      t4 = dlog(l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = dlog(1-l4)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t16 = dlog(-t1)
      t18 = dlog(l1)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = 1/t23
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t26 = 1/t25
      t28 = f1(0.D0,z,l1)
      t29 = t24*t26*t28
      t37 = 1-l8
      tncoff1 = t1*(-t3+t4+t8+t10+4*t11-2*t13+4*t16+t18+t20)*t29/1024+t1
     #/(l2-1)*t29/1024-t1/l8*t24*t26*t28/1024-t1/l3*t29/512
      t1 = z-1
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t4 = dlog(l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = dlog(1-l4)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t16 = dlog(-t1)
      t18 = dlog(l1)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t25 = s1vr2(z,l1,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,l1)
      t29 = 1/t23/t25*t28
      tncoff2 = t1*(-t3+t4+t8+t10+4*t11-2*t13+4*t16+t18+t20)*t29/1024+t1
     #/(-1+l8)*t29/1024+t1/(l2-1)*t29/1024-t1/l3*t29/512
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

      tcnum1 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = BCnum2e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = BCnum2e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = BCnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = BCnum2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = BCnum2e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = BCnum2e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      BC2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function BC2e2int(var,ndim,Tcoeff)

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
      double precision BCnum2e0
      external BCnum2e0
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

      tcnum1 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = BCnum2e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      BC2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

