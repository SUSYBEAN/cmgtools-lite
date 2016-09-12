CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgrm9e0int(var,ndim,Tcoeff)

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
      double precision qgrmnum7e0
      external qgrmnum7e0
      double precision qgrmnum7e1
      external qgrmnum7e1
      double precision qgrmnum7e2
      external qgrmnum7e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l4)
      t6 = dlog(l1)
      t10 = dlog(1-l4)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t25 = dlog(z)
      t28 = dlog(2.D0)
      t32 = dlog(1-z)
      t36 = dlog(1-l1)
      t45 = 3.D0/16384.D0*t3+3.D0/8192.D0*t5*t6-3.D0/8192.D0*t2*t10+3.D0
     #/8192.D0*t10*t6+3.D0/8192.D0*t5*t18-3.D0/8192.D0*t5*t2-3.D0/8192.D
     #0*t2*t6+3.D0/4096.D0*t2*t25-9.D0/8192.D0*t2*t28-3.D0/2048.D0*t32*t
     #2+3.D0/8192.D0*t10*t36+3.D0/8192.D0*t10*t18-3.D0/4096.D0*t18*t25+9
     #.D0/8192.D0*t18*t28
      t74 = 3.D0/8192.D0*t36*t18-3.D0/8192.D0*t36*t2+3.D0/2048.D0*t36*t3
     #2+3.D0/8192.D0*t36*t6-3.D0/4096.D0*t36*t25-3.D0/4096.D0*t6*t25+3.D
     #0/2048.D0*t18*t32-3.D0/8192.D0*t18*t2+3.D0/8192.D0*t18*t6+9.D0/819
     #2.D0*t6*t28+3.D0/8192.D0*t5*t36+3.D0/2048.D0*t32*t6+9.D0/8192.D0*t
     #36*t28+9.D0/8192.D0*t5*t28
      t95 = t6**2
      t98 = t18**2
      t101 = 3.D0/2048.D0*t5*t32-3.D0/1024.D0*t32*t25-9.D0/4096.D0*t25*t
     #28+9.D0/2048.D0*t28*t32-3.D0/4096.D0*t10*t25+9.D0/8192.D0*t10*t28+
     #3.D0/2048.D0*t10*t32+3.D0/8192.D0*t10*t5-3.D0/4096.D0*t5*t25-3.D0/
     #8192.D0*t10+3.D0/16384.D0*t95-3.D0/8192.D0*t6+3.D0/16384.D0*t98-3.
     #D0/8192.D0*t18
      t102 = t36**2
      t105 = t5**2
      t107 = t10**2
      t109 = t25**2
      t111 = t32**2
      t114 = 0.3141592653589793D1**2
      t116 = t28**2
      t122 = 3.D0/8192.D0+3.D0/16384.D0*t102-3.D0/8192.D0*t36+3.D0/16384
     #.D0*t105+3.D0/16384.D0*t107+3.D0/4096.D0*t109+3.D0/1024.D0*t111-9.
     #D0/8192.D0*t28-t114/16384+27.D0/16384.D0*t116-3.D0/2048.D0*t32+3.D
     #0/8192.D0*t2+3.D0/4096.D0*t25-3.D0/8192.D0*t5
      t125 = s23r2(z,l1,1.D0,0.D0,l4)
      t126 = t125**2
      t127 = 1/t126
      t129 = s1vr2(z,l1,1.D0,0.D0,l4)
      t130 = t129**2
      t131 = 1/t130
      t132 = f1(0.D0,z,l1)
      t133 = t131*t132
      t135 = dlog(l8)
      t136 = 4*t32
      t137 = 2*t25
      t138 = 3*t28
      t139 = 1-l8-t135-t136+t2-t6+t137-t138-t36-t10-t5
      t153 = -t18+t10*l8+t5*l8+t18*l8+4*t32*l8-t2*l8+t6*l8-2*t25*l8+3*t2
     #8*l8+t36*l8+t135*l8
      t155 = -1+l8
      t156 = 1/t155
      t158 = 1/l8
      t161 = t127*t131*t132
      t164 = dlog(l3)
      t166 = -2*t164-t18-t136+t2-t6+1+t137-t138-t36-t10-t5
      t167 = 1/l3
      t171 = -t155
      t178 = dlog(l2)
      t179 = -t5+t137-t138+t2-t136-t6-t18-t36+1-t178-t10
      t180 = 1/l2
      tncoff1 = (t45+t74+t101+t122)*t127*t133-3.D0/8192.D0*(t139+t153)*t
     #156*t158*t161+3.D0/4096.D0*t166*t167*t161-3.D0/4096.D0*t171*t167*t
     #156*t158*t127*t133+3.D0/8192.D0*t179*t180*t161-3.D0/8192.D0*t171*t
     #180*t158*t156*t127*t133+3.D0/4096.D0*t180*t167*t161
      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l4)
      t6 = dlog(l1)
      t10 = dlog(1-l4)
      t15 = f3(0.D0,z)
      t18 = dlog(1-l1*t15)
      t25 = dlog(z)
      t28 = dlog(2.D0)
      t32 = dlog(1-z)
      t36 = dlog(1-l1)
      t43 = 3.D0/8192.D0+3.D0/16384.D0*t3+3.D0/8192.D0*t5*t6-3.D0/8192.D
     #0*t2*t10+3.D0/8192.D0*t10*t6+3.D0/8192.D0*t5*t18-3.D0/8192.D0*t5*t
     #2-3.D0/8192.D0*t2*t6+3.D0/4096.D0*t2*t25-9.D0/8192.D0*t2*t28-3.D0/
     #2048.D0*t32*t2+3.D0/8192.D0*t10*t36+3.D0/8192.D0*t10*t18-3.D0/4096
     #.D0*t18*t25
      t72 = 9.D0/8192.D0*t18*t28+3.D0/8192.D0*t36*t18-3.D0/8192.D0*t36*t
     #2+3.D0/2048.D0*t36*t32+3.D0/8192.D0*t36*t6-3.D0/4096.D0*t36*t25-3.
     #D0/4096.D0*t6*t25+3.D0/2048.D0*t18*t32-3.D0/8192.D0*t18*t2+3.D0/81
     #92.D0*t18*t6+9.D0/8192.D0*t6*t28+3.D0/8192.D0*t5*t36+3.D0/2048.D0*
     #t32*t6+9.D0/8192.D0*t36*t28
      t95 = t5**2
      t97 = t10**2
      t99 = t25**2
      t101 = 9.D0/8192.D0*t5*t28+3.D0/2048.D0*t5*t32-3.D0/1024.D0*t32*t2
     #5-9.D0/4096.D0*t25*t28+9.D0/2048.D0*t28*t32-3.D0/4096.D0*t10*t25+9
     #.D0/8192.D0*t10*t28+3.D0/2048.D0*t10*t32+3.D0/8192.D0*t10*t5-3.D0/
     #4096.D0*t5*t25-3.D0/8192.D0*t36+3.D0/16384.D0*t95+3.D0/16384.D0*t9
     #7+3.D0/4096.D0*t99
      t102 = t32**2
      t109 = t6**2
      t112 = t18**2
      t115 = t36**2
      t118 = 0.3141592653589793D1**2
      t120 = t28**2
      t122 = 3.D0/1024.D0*t102-3.D0/2048.D0*t32+3.D0/8192.D0*t2+3.D0/409
     #6.D0*t25-3.D0/8192.D0*t5-3.D0/8192.D0*t10+3.D0/16384.D0*t109-3.D0/
     #8192.D0*t6+3.D0/16384.D0*t112-3.D0/8192.D0*t18+3.D0/16384.D0*t115-
     #9.D0/8192.D0*t28-t118/16384+27.D0/16384.D0*t120
      t125 = s23r2(z,l1,1.D0,0.D0,l4)
      t126 = t125**2
      t127 = 1/t126
      t129 = s1vr2(z,l1,1.D0,0.D0,l4)
      t130 = t129**2
      t131 = 1/t130
      t132 = f1(0.D0,z,l1)
      t135 = 1-l8
      t136 = dlog(t135)
      t150 = -t136*l8-t18*l8-4*t32*l8+t2*l8-t6*l8+2*t25*l8-3*t28*l8-t36*
     #l8-t10*l8-t5*l8+l8
      t152 = -1/t135
      t157 = t127*t131*t132
      t160 = dlog(l3)
      t162 = 4*t32
      t163 = 2*t25
      t164 = 3*t28
      t165 = -2*t160-t18-t162+t2-t6+1+t163-t164-t36-t10-t5
      t166 = 1/l3
      t173 = dlog(l2)
      t174 = -t5+t163-t164+t2-t162-t6-t18-t36+1-t173-t10
      t175 = 1/l2
      tncoff2 = (t43+t72+t101+t122)*t127*t131*t132-3.D0/8192.D0*t150*t15
     #2/l8*t157+3.D0/4096.D0*t165*t166*t157-3.D0/4096.D0*t166*t152*t157+
     #3.D0/8192.D0*t174*t175*t157-3.D0/8192.D0*t175*t152*t157+3.D0/4096.
     #D0*t175*t166*t157
      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t4 = dlog(1-l1)
      t6 = dlog(1-z)
      t8 = 1-l8
      t9 = dlog(t8)
      t10 = dlog(l8)
      t12 = dlog(1-l4)
      t13 = dlog(l4)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t18 = dlog(l1)
      t19 = f3(0.D0,z)
      t22 = dlog(1-l1*t19)
      t23 = -1-t2+t4+4*t6+t9+t10+t12+t13+3*t14-2*t16+t18+t22
      t25 = -1/t8
      t27 = 1/l8
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,1.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff3 = -3.D0/8192.D0*t23*t25*t27*t37+3.D0/4096.D0/l3*t25*t27*t3
     #7+3.D0/8192.D0/l2*t27*t25*t37
      t1 = l3-1
      t2 = dlog(l3)
      t5 = dlog(1-l4)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(l4)
      t9 = dlog(l1)
      t11 = dlog(-t1)
      t13 = f3(l3,z)
      t16 = dlog(1-l1*t13)
      t18 = dlog(1-z)
      t21 = dlog(1-l1)
      t22 = dlog(z)
      t24 = dlog(2.D0)
      t26 = 2*t2-1+t5-t7+t8+t9+2*t11+t16+4*t18+t21-2*t22+3*t24
      t28 = 1/l3
      t30 = s23r2(z,l1,1.D0,l3,l4)
      t31 = t30**2
      t33 = s1vr2(z,l1,1.D0,l3,l4)
      t34 = t33**2
      t37 = f1(l3,z,l1)
      t38 = 1/t31/t34*t37
      tncoff4 = -3.D0/4096.D0*t1*t26*t28*t38+3.D0/4096.D0*t1*t28/l8*t38+
     #3.D0/4096.D0*t1/l2*t28*t38
      t1 = l3-1
      t2 = dlog(l3)
      t5 = dlog(1-l4)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t8 = dlog(l4)
      t9 = dlog(l1)
      t11 = dlog(-t1)
      t13 = f3(l3,z)
      t16 = dlog(1-l1*t13)
      t18 = dlog(1-z)
      t21 = dlog(1-l1)
      t22 = dlog(z)
      t24 = dlog(2.D0)
      t26 = 2*t2-1+t5-t7+t8+t9+2*t11+t16+4*t18+t21-2*t22+3*t24
      t28 = 1/l3
      t30 = s23r2(z,l1,1.D0,l3,l4)
      t31 = t30**2
      t33 = s1vr2(z,l1,1.D0,l3,l4)
      t34 = t33**2
      t37 = f1(l3,z,l1)
      t38 = 1/t31/t34*t37
      tncoff5 = -3.D0/4096.D0*t1*t26*t28*t38-3.D0/4096.D0*t1*t28/(-1+l8)
     #*t38+3.D0/4096.D0*t1/l2*t28*t38
      t8 = s23r2(z,l1,1.D0,l3,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,1.D0,l3,l4)
      t13 = t12**2
      t15 = f1(l3,z,l1)
      tncoff6 = 3.D0/4096.D0*(l3-1)/l3/(-1+l8)/l8/t9/t13*t15
      t1 = dlog(l4)
      t2 = dlog(2.D0)
      t5 = dlog(1-l1)
      t6 = f3(0.D0,z)
      t9 = dlog(1-l1*t6)
      t10 = dlog(z)
      t13 = dlog(1-z)
      t15 = dlog(l1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = 2-l2
      t19 = dlog(t18)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = t1+2*t2+t5+t9-2*t10+4*t13+t15-t17-1+t19+t20+t22
      t25 = -1/t18
      t27 = 1/l2
      t30 = 1-l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = s1vr2(z,l1,t30,0.D0,l4)
      t35 = t34**2
      t36 = 1/t35
      t38 = f1(0.D0,z,l1)
      t39 = t33*t36*t38
      tncoff7 = -3.D0/4096.D0*t23*t25*t27*t39+3.D0/4096.D0*t25*t27/l8*t3
     #9+3.D0/2048.D0*t25*t36*t33*t27/l3*t38
      t1 = dlog(l4)
      t2 = dlog(2.D0)
      t5 = dlog(1-l1)
      t6 = f3(0.D0,z)
      t9 = dlog(1-l1*t6)
      t10 = dlog(z)
      t13 = dlog(1-z)
      t15 = dlog(l1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = 2-l2
      t19 = dlog(t18)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = t1+2*t2+t5+t9-2*t10+4*t13+t15-t17-1+t19+t20+t22
      t25 = -1/t18
      t27 = 1/l2
      t30 = 1-l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = s1vr2(z,l1,t30,0.D0,l4)
      t35 = t34**2
      t36 = 1/t35
      t38 = f1(0.D0,z,l1)
      t39 = t33*t36*t38
      tncoff8 = -3.D0/4096.D0*t23*t25*t27*t39-3.D0/4096.D0*t25*t27/(-1+l
     #8)*t39+3.D0/2048.D0*t25*t36*t33*t27/l3*t38
      t10 = 1-l2/2
      t11 = s23r2(z,l1,t10,0.D0,l4)
      t12 = t11**2
      t15 = s1vr2(z,l1,t10,0.D0,l4)
      t16 = t15**2
      t18 = f1(0.D0,z,l1)
      tncoff9 = 3.D0/4096.D0/(l2-2)/l2/l8/(-1+l8)/t12/t16*t18
      t9 = 1-l2/2
      t10 = s23r2(z,l1,t9,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,t9,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff10 = 3.D0/2048.D0*(l3-1)/(l2-2)/l2/l3/t11/t15*t17
      t9 = 1-l2/2
      t10 = s23r2(z,l1,t9,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,t9,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff11 = 3.D0/2048.D0*(l3-1)/(l2-2)/l2/l3/t11/t15*t17
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = dlog(l4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t19 = dlog(1-z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(l1)
      t27 = s23r2(z,l1,1.D0,0.D0,l4)
      t28 = t27**2
      t29 = 1/t28
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = f1(0.D0,z,l1)
      t35 = t33*t34
      t37 = 1-l8
      tncoff12 = (3.D0/4096.D0*t1-9.D0/8192.D0*t3-3.D0/8192.D0*t5-3.D0/8
     #192.D0*t8-3.D0/8192.D0*t11-3.D0/8192.D0*t16-3.D0/2048.D0*t19+3.D0/
     #8192.D0*t22-3.D0/8192.D0*t24+3.D0/8192.D0)*t29*t35+3.D0/8192.D0/l8
     #*t29*t33*t34+3.D0/4096.D0/l3*t29*t35+3.D0/8192.D0/l2*t29*t35
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = dlog(l4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t19 = dlog(1-z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(l1)
      t27 = s23r2(z,l1,1.D0,0.D0,l4)
      t28 = t27**2
      t29 = 1/t28
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = t31**2
      t34 = f1(0.D0,z,l1)
      t35 = 1/t32*t34
      tncoff13 = (3.D0/4096.D0*t1-9.D0/8192.D0*t3-3.D0/8192.D0*t5-3.D0/8
     #192.D0*t8-3.D0/8192.D0*t11-3.D0/8192.D0*t16-3.D0/2048.D0*t19+3.D0/
     #8192.D0*t22-3.D0/8192.D0*t24+3.D0/8192.D0)*t29*t35-3.D0/8192.D0/(-
     #1+l8)*t29*t35+3.D0/4096.D0/l3*t29*t35+3.D0/8192.D0/l2*t29*t35
      t5 = s23r2(z,l1,1.D0,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(z,l1,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff14 = 3.D0/8192.D0/(-1+l8)/l8/t6/t9*t12
      t4 = s23r2(z,l1,1.D0,l3,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,l1)
      tncoff15 = 3.D0/4096.D0*(l3-1)/l3/t5/t8*t11
      t4 = s23r2(z,l1,1.D0,l3,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,l1)
      tncoff16 = 3.D0/4096.D0*(l3-1)/l3/t5/t8*t11
      t6 = 1-l2/2
      t7 = s23r2(z,l1,t6,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,t6,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff17 = 3.D0/4096.D0/(l2-2)/l2/t8/t11*t14
      t6 = 1-l2/2
      t7 = s23r2(z,l1,t6,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,t6,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff18 = 3.D0/4096.D0/(l2-2)/l2/t8/t11*t14
      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff19 = 3.D0/8192.D0/t2/t5*t8
      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff20 = 3.D0/8192.D0/t2/t5*t8

CCC Numerator functions

      tcnum1 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgrmnum7e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgrmnum7e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgrmnum7e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgrmnum7e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgrm9e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12
     #+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff
     #16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnco
     #ff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm9e1int(var,ndim,Tcoeff)

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
      double precision qgrmnum7e0
      external qgrmnum7e0
      double precision qgrmnum7e1
      external qgrmnum7e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = dlog(l4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t19 = dlog(1-z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(l1)
      t27 = s23r2(z,l1,1.D0,0.D0,l4)
      t28 = t27**2
      t29 = 1/t28
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = f1(0.D0,z,l1)
      t35 = t33*t34
      t37 = 1-l8
      tncoff1 = (3.D0/4096.D0*t1-9.D0/8192.D0*t3-3.D0/8192.D0*t5-3.D0/81
     #92.D0*t8-3.D0/8192.D0*t11-3.D0/8192.D0*t16-3.D0/2048.D0*t19+3.D0/8
     #192.D0*t22-3.D0/8192.D0*t24+3.D0/8192.D0)*t29*t35+3.D0/8192.D0/l8*
     #t29*t33*t34+3.D0/4096.D0/l3*t29*t35+3.D0/8192.D0/l2*t29*t35
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = dlog(l4)
      t8 = dlog(1-l4)
      t11 = dlog(1-l1)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t19 = dlog(1-z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(l1)
      t27 = s23r2(z,l1,1.D0,0.D0,l4)
      t28 = t27**2
      t29 = 1/t28
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = t31**2
      t34 = f1(0.D0,z,l1)
      t35 = 1/t32*t34
      tncoff2 = (3.D0/4096.D0*t1-9.D0/8192.D0*t3-3.D0/8192.D0*t5-3.D0/81
     #92.D0*t8-3.D0/8192.D0*t11-3.D0/8192.D0*t16-3.D0/2048.D0*t19+3.D0/8
     #192.D0*t22-3.D0/8192.D0*t24+3.D0/8192.D0)*t29*t35-3.D0/8192.D0/(-1
     #+l8)*t29*t35+3.D0/4096.D0/l3*t29*t35+3.D0/8192.D0/l2*t29*t35
      t5 = s23r2(z,l1,1.D0,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(z,l1,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff3 = 3.D0/8192.D0/(-1+l8)/l8/t6/t9*t12
      t4 = s23r2(z,l1,1.D0,l3,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,l1)
      tncoff4 = 3.D0/4096.D0*(l3-1)/l3/t5/t8*t11
      t4 = s23r2(z,l1,1.D0,l3,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,l1)
      tncoff5 = 3.D0/4096.D0*(l3-1)/l3/t5/t8*t11
      t6 = 1-l2/2
      t7 = s23r2(z,l1,t6,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,t6,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff6 = 3.D0/4096.D0/(l2-2)/l2/t8/t11*t14
      t6 = 1-l2/2
      t7 = s23r2(z,l1,t6,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,t6,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff7 = 3.D0/4096.D0/(l2-2)/l2/t8/t11*t14
      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff8 = 3.D0/8192.D0/t2/t5*t8
      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff9 = 3.D0/8192.D0/t2/t5*t8

CCC Numerator functions

      tcnum1 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgrmnum7e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgrmnum7e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgrm9e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm9e2int(var,ndim,Tcoeff)

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
      double precision qgrmnum7e0
      external qgrmnum7e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff1 = 3.D0/8192.D0/t2/t5*t8
      t1 = s23r2(z,l1,1.D0,0.D0,l4)
      t2 = t1**2
      t4 = s1vr2(z,l1,1.D0,0.D0,l4)
      t5 = t4**2
      t8 = f1(0.D0,z,l1)
      tncoff2 = 3.D0/8192.D0/t2/t5*t8

CCC Numerator functions

      tcnum1 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum7e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgrm9e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm9e3int(var,ndim,Tcoeff)

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
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      qgrm9e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

