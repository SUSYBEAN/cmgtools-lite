CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm2e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
      double precision qiqirmnum1e2
      external qiqirmnum1e2
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
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t9 = dlog(2.D0)
      t13 = dlog(1-l4)
      t18 = dlog(z)
      t21 = dlog(l1)
      t25 = dlog(1-l1)
      t36 = 6*t2*t6+18*t6*t9+18*t13*t9+6*t13*t6-12*t2*t18+6*t2*t21+6*t13
     #*t25+6*t13*t21-12*t13*t18+6*t6*t25+6*t6*t21
      t37 = t25**2
      t40 = dlog(-t1)
      t41 = t40**2
      t43 = f2(0.D0,z)
      t44 = dlog(t43)
      t45 = t44**2
      t47 = t13**2
      t49 = t21**2
      t65 = 3*t37+48*t41+3*t45+3*t47+3*t49-6*t6*t44+24*t6*t40+24*t13*t40
     #-6*t13*t44-12*t6*t18-6*t44*t21+6*t21*t25
      t89 = 18*t21*t9-12*t21*t18-6*t44*t25+72*t40*t9-12*t25*t18-24*t40*t
     #44+24*t40*t21-48*t40*t18+24*t40*t25+18*t2*t9-18*t44*t9
      t96 = 0.3141592653589793D1**2
      t97 = t2**2
      t99 = t6**2
      t101 = t9**2
      t103 = t18**2
      t113 = 12*t44*t18+18*t25*t9+6*t2*t25-t96+3*t97+3*t99+27*t101+12*t1
     #03+6*t13*t2-36*t18*t9+24*t2*t40-6*t2*t44
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,1.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t128 = 4*t40
      t129 = 3*t9
      t130 = 2*t18
      t131 = dlog(l8)
      t145 = -t128-t6-t129-t21-t25-t2+t130-t13+t44-t131+t13*l8+3*t9*l8+t
     #6*l8+t21*l8+t25*l8-2*t18*l8+t131*l8+t2*l8+4*t40*l8-t44*l8
      t150 = 1/l8
      t156 = dlog(l3)
      t160 = 1/l3
      t168 = dlog(l2)
      t171 = 1/l2
      t175 = t1*t171
      tncoff1 = -t1*(t36+t65+t89+t113)*t125/12288+t1*t145/(-1+l8)*t150*t
     #119*t122*t124/2048+t1*(2*t156+t128+t6+t129+t21+t25-t130-t44+t13+t2
     #)*t160*t125/1024-t1*t150*t160*t125/1024+(t168+t129+t25+t21+t6+t128
     #-t130-t44+t13+t2)*t1*t171*t125/2048-t175*t150*t125/2048-t175*t160*
     #t125/1024
      t1 = -1+z
      t2 = dlog(l1)
      t4 = dlog(1-l1)
      t7 = dlog(l4)
      t11 = dlog(1-l4)
      t12 = dlog(2.D0)
      t16 = dlog(-t1)
      t21 = f3(0.D0,z)
      t24 = dlog(1-l1*t21)
      t29 = dlog(z)
      t38 = 6*t2*t4+6*t7*t2+18*t11*t12+72*t16*t12+24*t16*t2+6*t24*t2+6*t
     #11*t4-12*t7*t29-12*t2*t29+18*t2*t12+6*t11*t24
      t41 = f2(0.D0,z)
      t42 = dlog(t41)
      t65 = 18*t24*t12-6*t42*t4-12*t24*t29-24*t16*t42+6*t7*t24+12*t42*t2
     #9+6*t11*t2-18*t42*t12+24*t24*t16-48*t16*t29+6*t24*t4+18*t7*t12
      t87 = 0.3141592653589793D1**2
      t88 = 18*t4*t12+24*t16*t4+24*t7*t16+6*t11*t7+6*t7*t4-36*t29*t12-6*
     #t7*t42-12*t11*t29-12*t4*t29+24*t11*t16-t87
      t89 = t12**2
      t91 = t4**2
      t93 = t16**2
      t95 = t42**2
      t97 = t11**2
      t99 = t2**2
      t103 = t7**2
      t105 = t24**2
      t107 = t29**2
      t113 = 27*t89+3*t91+48*t93+3*t95+3*t97+3*t99-6*t42*t2+3*t103+3*t10
     #5+12*t107-6*t11*t42-6*t24*t42
      t117 = s23r2(z,l1,1.D0,0.D0,l4)
      t118 = t117**2
      t119 = 1/t118
      t120 = s1vr2(z,l1,1.D0,0.D0,l4)
      t121 = t120**2
      t122 = 1/t121
      t124 = f1(0.D0,z,l1)
      t125 = t119*t122*t124
      t139 = 1-l8
      t140 = dlog(t139)
      t146 = -1/t139
      t154 = dlog(l3)
      t156 = 4*t16
      t157 = 3*t12
      t158 = 2*t29
      t161 = 1/l3
      t169 = dlog(l2)
      t172 = 1/l2
      tncoff2 = -t1*(t38+t65+t88+t113)*t125/12288+t1*(-t7*l8-t2*l8-t4*l8
     #+2*t29*l8-t24*l8-4*t16*l8+t42*l8-3*t12*l8-t140*l8-t11*l8)*t146/l8*
     #t119*t122*t124/2048+t1*(2*t154+t156+t24+t157+t2+t4-t158-t42+t11+t7
     #)*t161*t125/1024+t1*t161*t146*t125/1024+(t169+t157+t4+t2+t24+t156-
     #t158-t42+t11+t7)*t1*t172*t125/2048+t1*t146*t172*t125/2048-t1*t172*
     #t161*t125/1024
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(z)
      t6 = dlog(1-l1)
      t7 = dlog(l1)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(-t1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t17 = dlog(2.D0)
      t19 = dlog(l8)
      t21 = dlog(1-l4)
      t22 = 1-l8
      t23 = dlog(t22)
      t24 = t2-2*t3+t6+t7-t9+4*t11+t16+3*t17+t19+t21+t23
      t27 = -1/t22
      t29 = 1/l8
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,1.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      t39 = t29*t32*t38
      tncoff3 = t1*t24*t27*t39/2048-t1*t29/l3*t27*t32*t38/1024-t1*t27/l2
     #*t39/2048
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(2.D0)
      t6 = f3(l3,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t15 = dlog(l3)
      t17 = dlog(z)
      t19 = f2(l3,z)
      t20 = dlog(t19)
      t22 = dlog(1-l4)
      t23 = dlog(l1)
      t24 = dlog(l4)
      t26 = dlog(-t2)
      t28 = 3*t4+t9+t11+4*t13+2*t15-2*t17-t20+t22+t23+t24+2*t26
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff4 = t3*t28*t40/1024-t3/l8*t40/1024-t3/l2*t40/1024
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(2.D0)
      t6 = f3(l3,z)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t15 = dlog(l3)
      t17 = dlog(z)
      t19 = f2(l3,z)
      t20 = dlog(t19)
      t22 = dlog(1-l4)
      t23 = dlog(l1)
      t24 = dlog(l4)
      t26 = dlog(-t2)
      t28 = 3*t4+t9+t11+4*t13+2*t15-2*t17-t20+t22+t23+t24+2*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,1.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      t40 = t30*t33*t39
      tncoff5 = t3*t28*t40/1024+t3*t30/(-1+l8)*t33*t39/1024-t3/l2*t40/10
     #24
      t10 = s23r2(z,l1,1.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,1.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l8/l3/(-1+l8)/t11/t15*t17/1024
      t1 = -1+z
      t2 = dlog(l2)
      t3 = dlog(z)
      t5 = dlog(l1)
      t7 = dlog(-t1)
      t9 = dlog(2.D0)
      t12 = dlog(1-l1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t15 = -l2+2
      t16 = dlog(t15)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t2-2*t3+t5+4*t7+2*t9+t12-t14+t16+t20+t22+t23
      t27 = -1/t15
      t29 = 1/l2
      t31 = 1-l2/2
      t32 = s23r2(z,l1,t31,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t35 = t29*t34
      t36 = s1vr2(z,l1,t31,0.D0,l4)
      t37 = t36**2
      t38 = 1/t37
      t39 = f1(0.D0,z,l1)
      t44 = 1-l8
      tncoff7 = t1*t24*t27*t35*t38*t39/1024-t1*t27*t29/l8*t34*t38*t39/10
     #24-t1*t38*t27*t35/l3*t39/512
      t1 = -1+z
      t2 = dlog(l2)
      t3 = dlog(z)
      t5 = dlog(l1)
      t7 = dlog(-t1)
      t9 = dlog(2.D0)
      t12 = dlog(1-l1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t15 = -l2+2
      t16 = dlog(t15)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t2-2*t3+t5+4*t7+2*t9+t12-t14+t16+t20+t22+t23
      t27 = -1/t15
      t29 = 1/l2
      t31 = 1-l2/2
      t32 = s23r2(z,l1,t31,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t35 = t29*t34
      t36 = s1vr2(z,l1,t31,0.D0,l4)
      t37 = t36**2
      t38 = 1/t37
      t39 = f1(0.D0,z,l1)
      t40 = t38*t39
      tncoff8 = t1*t24*t27*t35*t40/1024+t1*t27*t29/(-1+l8)*t34*t40/1024-
     #t1*t38*t27*t35/l3*t39/512
      t12 = 1-l2/2
      t13 = s23r2(z,l1,t12,0.D0,l4)
      t14 = t13**2
      t17 = s1vr2(z,l1,t12,0.D0,l4)
      t18 = t17**2
      t20 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/(l2-2)/l2/l8/(-1+l8)/t14/t18*t20/1024
      t11 = 1-l2/2
      t12 = s23r2(z,l1,t11,l3,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,l3,l4)
      t17 = t16**2
      t19 = f1(l3,z,l1)
      tncoff10 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t13/t17*t19/512
      t11 = 1-l2/2
      t12 = s23r2(z,l1,t11,l3,l4)
      t13 = t12**2
      t16 = s1vr2(z,l1,t11,l3,l4)
      t17 = t16**2
      t19 = f1(l3,z,l1)
      tncoff11 = -(-1+z)*(l3-1)/(l2-2)/l2/l3/t13/t17*t19/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = dlog(2.D0)
      t13 = dlog(-t1)
      t15 = dlog(z)
      t17 = dlog(l4)
      t20 = dlog(1-l1)
      t22 = dlog(l1)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = s1vr2(z,l1,1.D0,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t37 = t31*t34*t36
      t40 = 1-l8
      tncoff12 = -t1*(-6*t3-6*t8-18*t10-24*t13+12*t15-6*t17-6*t20-6*t22+
     #6*t25)*t37/12288-t1/l8*t31*t34*t36/2048-t1/l3*t37/1024-t1/l2*t37/2
     #048
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = dlog(2.D0)
      t13 = dlog(-t1)
      t15 = dlog(z)
      t17 = dlog(l4)
      t20 = dlog(1-l1)
      t22 = dlog(l1)
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,1.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t3-6*t8-18*t10-24*t13+12*t15-6*t17-6*t20-6*t22+
     #6*t25)*t37/12288+t1/(-1+l8)*t37/2048-t1/l3*t37/1024-t1/l2*t37/2048
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff15 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff16 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff17 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff18 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff19 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff20 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum1e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum1e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm2e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
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
      t3 = dlog(z)
      t6 = dlog(1-l4)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(-t1)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,1.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = t1*(t2-2*t3+t6-t8+4*t10+t15+3*t16+t18+t20)*t31/2048-t1/l
     #8*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(z)
      t6 = dlog(1-l4)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(-t1)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t20 = dlog(1-l1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,1.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = t1*(t2-2*t3+t6-t8+4*t10+t15+3*t16+t18+t20)*t31/2048+t1/(
     #-1+l8)*t31/2048-t1/l3*t31/1024-t1/l2*t31/2048
      t7 = s23r2(z,l1,1.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t8/t11*t14/2048
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff4 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t6 = s23r2(z,l1,1.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,1.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff5 = -(-1+z)*(l3-1)/l3/t7/t10*t13/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff6 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t8 = 1-l2/2
      t9 = s23r2(z,l1,t8,0.D0,l4)
      t10 = t9**2
      t12 = s1vr2(z,l1,t8,0.D0,l4)
      t13 = t12**2
      t16 = f1(0.D0,z,l1)
      tncoff7 = -(-1+z)/(l2-2)/l2/t10/t13*t16/1024
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff9 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm2e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
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
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -(-1+z)/t3/t7*t9/2048
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -(-1+z)/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

