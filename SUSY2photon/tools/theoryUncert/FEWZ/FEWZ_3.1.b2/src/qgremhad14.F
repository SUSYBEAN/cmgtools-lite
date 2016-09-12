CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgrm14e0int(var,ndim,Tcoeff)

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
      double precision qgrmnum11e0
      external qgrmnum11e0
      double precision qgrmnum11e1
      external qgrmnum11e1
      double precision qgrmnum11e2
      external qgrmnum11e2
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

      t1 = -1+z
      t3 = dlog(-t1)
      t4 = t3**2
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = t7**2
      t10 = dlog(z)
      t11 = t10**2
      t14 = dlog(1-l1)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t24 = f3(0.D0,z)
      t27 = dlog(1-l1*t24)
      t28 = t27**2
      t30 = dlog(l1)
      t31 = t30**2
      t43 = 6+48*t4+3*t8+12*t11+3*t15+3*t18+3*t22+3*t28+3*t31+6*t17*t30+
     #6*t21*t30+6*t21*t14+24*t14*t3+6*t14*t30
      t48 = dlog(2.D0)
      t73 = -12*t14*t10-12*t30*t10+24*t30*t48+6*t17*t14+24*t3*t30+24*t14
     #*t48+24*t17*t48+24*t17*t3-48*t3*t10-48*t10*t48+96*t48*t3-12*t21*t1
     #0+24*t21*t48+24*t21*t3
      t102 = 6*t21*t17-12*t17*t10+6*t14*t27-6*t14*t7+6*t17*t27-6*t17*t7-
     #6*t7*t30+12*t7*t10-24*t7*t48-24*t3*t7-6*t17+6*t21*t27-6*t21*t7+24*
     #t27*t3
      t118 = t48**2
      t120 = 0.3141592653589793D1**2
      t122 = -6*t27*t7+6*t27*t30-12*t27*t10+24*t27*t48-6*t14-6*t30+12*t1
     #0-24*t3+6*t7-6*t27-6*t21+48*t118-t120-24*t48
      t126 = s23r2(z,l1,0.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = s1vr2(z,l1,0.D0,0.D0,l4)
      t130 = t129**2
      t131 = 1/t130
      t133 = f1(0.D0,z,l1)
      t134 = t128*t131*t133
      t141 = dlog(l8)
      t151 = 1+t14*l8+t21*l8+t27*l8-t7*l8+t141*l8+4*t3*l8+t30*l8-2*t10*l
     #8+4*t48*l8+t17*l8
      t152 = 4*t3
      t153 = 4*t48
      t154 = 2*t10
      t155 = -t21+t7-t152-t27-t153-t30+t154-t14-t17-l8-t141
      t158 = 1/l8
      t167 = dlog(l3)
      t169 = t17+t21+t152+t27-t7+t30-t154-1+t153+t14+2*t167
      t171 = 1/l3
      t179 = dlog(l2)
      t180 = t17+t21+t179+t152+t27-t7+t30-t154-1+t153+t14
      t182 = 1/l2
      tncoff1 = -t1*(t43+t73+t102+t122)*t134/16384+3.D0/8192.D0*t1*(t151
     #+t155)*t158/(-1+l8)*t128*t131*t133+3.D0/4096.D0*t1*t169*t171*t134-
     #3.D0/4096.D0*t1*t158*t171*t134+3.D0/8192.D0*t1*t180*t182*t134-3.D0
     #/8192.D0*t1*t182*t158*t134-3.D0/4096.D0*t1*t128*t182*t171*t131*t13
     #3
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = t3**2
      t6 = dlog(z)
      t7 = t6**2
      t10 = dlog(-t1)
      t11 = t10**2
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t15 = t14**2
      t17 = dlog(l1)
      t18 = t17**2
      t20 = f3(0.D0,z)
      t23 = dlog(1-l1*t20)
      t24 = t23**2
      t27 = dlog(1-l1)
      t28 = t27**2
      t30 = dlog(l4)
      t31 = t30**2
      t38 = 6+3*t4+12*t7+48*t11+3*t15+3*t18+3*t24+3*t28+3*t31-6*t30-6*t2
     #7+12*t6-6*t17-6*t23
      t42 = dlog(2.D0)
      t43 = t42**2
      t45 = 0.3141592653589793D1**2
      t63 = 6*t14-24*t10-6*t3+48*t43-t45-24*t42-48*t10*t6-48*t6*t42+96*t
     #42*t10-12*t3*t6+24*t3*t42+24*t3*t10+6*t3*t30-12*t30*t6
      t93 = 24*t27*t10+6*t17*t27-12*t27*t6-12*t17*t6+24*t17*t42+6*t27*t3
     #0+24*t17*t10+24*t27*t42+24*t30*t42+24*t30*t10+6*t17*t30+6*t17*t3+6
     #*t27*t3-24*t14*t10
      t122 = 24*t23*t42+6*t27*t23-6*t14*t27+6*t30*t23-6*t14*t30-6*t17*t1
     #4+12*t14*t6-24*t14*t42+6*t3*t23-6*t14*t3+24*t23*t10-6*t23*t14+6*t2
     #3*t17-12*t23*t6
      t126 = s23r2(z,l1,0.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = s1vr2(z,l1,0.D0,0.D0,l4)
      t130 = t129**2
      t131 = 1/t130
      t133 = f1(0.D0,z,l1)
      t134 = t128*t131*t133
      t148 = 1-l8
      t149 = dlog(t148)
      t152 = l8-t30*l8-4*t10*l8-t23*l8+t14*l8-t17*l8+2*t6*l8-4*t42*l8-t2
     #7*l8-t149*l8-t3*l8
      t157 = -1/t148
      t163 = 4*t10
      t164 = 2*t6
      t165 = 4*t42
      t166 = dlog(l3)
      t168 = t30+t3+t163+t23-t14+t17-t164-1+t165+t27+2*t166
      t170 = 1/l3
      t178 = dlog(l2)
      t179 = t30+t3+t178+t163+t23-t14+t17-t164-1+t165+t27
      t181 = 1/l2
      tncoff2 = -t1*(t38+t63+t93+t122)*t134/16384+3.D0/8192.D0*t1*t152/l
     #8*t157*t128*t131*t133+3.D0/4096.D0*t1*t168*t170*t134+3.D0/4096.D0*
     #t1*t157*t170*t134+3.D0/8192.D0*t1*t179*t181*t134+3.D0/8192.D0*t1*t
     #181*t157*t134-3.D0/4096.D0*t1*t128*t181*t170*t131*t133
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t5 = dlog(1-l4)
      t6 = dlog(l4)
      t7 = dlog(l8)
      t9 = dlog(-t1)
      t11 = f3(0.D0,z)
      t14 = dlog(1-l1*t11)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t18 = dlog(z)
      t21 = dlog(1-l1)
      t22 = dlog(2.D0)
      t24 = -1+t3+t5+t6+t7+4*t9+t14-t16+t17-2*t18+t21+4*t22
      t26 = 1/l8
      t29 = -1/t2
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t37 = f1(0.D0,z,l1)
      t38 = 1/t35*t37
      tncoff3 = 3.D0/8192.D0*t1*t24*t26*t29*t32*t38-3.D0/4096.D0*t1*t29*
     #t26/l3*t32*t38-3.D0/8192.D0*t1/l2*t29*t26*t32*t38
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l4)
      t6 = dlog(1-l4)
      t7 = dlog(l3)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t12 = f3(l3,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(z)
      t18 = f2(l3,z)
      t19 = dlog(t18)
      t21 = dlog(-t1)
      t24 = dlog(-t2)
      t27 = dlog(1-l1)
      t28 = t4+t6+2*t7+4*t9+t11+t15-2*t16-t19+4*t21+2*t24+t27-1
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff4 = 3.D0/4096.D0*t3*t28*t40-3.D0/4096.D0*t3/l8*t40-3.D0/4096
     #.D0*t3/l2*t40
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l4)
      t6 = dlog(1-l4)
      t7 = dlog(l3)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t12 = f3(l3,z)
      t15 = dlog(1-l1*t12)
      t16 = dlog(z)
      t18 = f2(l3,z)
      t19 = dlog(t18)
      t21 = dlog(-t1)
      t24 = dlog(-t2)
      t27 = dlog(1-l1)
      t28 = t4+t6+2*t7+4*t9+t11+t15-2*t16-t19+4*t21+2*t24+t27-1
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff5 = 3.D0/4096.D0*t3*t28*t40+3.D0/4096.D0*t3/(-1+l8)*t40-3.D0
     #/4096.D0*t3/l2*t40
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t14 = s1vr2(z,l1,0.D0,l3,l4)
      t15 = t14**2
      t17 = f1(l3,z,l1)
      tncoff6 = -3.D0/4096.D0*(-1+z)*(l3-1)/(-1+l8)/l8/l3/t11/t15*t17
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = dlog(l2)
      t6 = dlog(2.D0)
      t8 = dlog(z)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t16 = dlog(1-l2)
      t17 = dlog(l1)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t24 = t2+t4+t5-1+4*t6-2*t8+t11+4*t13+t16+t17+t21-t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff7 = -3.D0/8192.D0*t1*t24*t26*t36+3.D0/8192.D0*t1*t26/l8*t36+
     #3.D0/4096.D0*t1*t33*t30*t26/l3*t35
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = dlog(l2)
      t6 = dlog(2.D0)
      t8 = dlog(z)
      t11 = dlog(1-l1)
      t13 = dlog(-t1)
      t16 = dlog(1-l2)
      t17 = dlog(l1)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t24 = t2+t4+t5-1+4*t6-2*t8+t11+4*t13+t16+t17+t21-t23
      t26 = 1/l2
      t28 = s23r2(z,l1,l2,0.D0,l4)
      t29 = t28**2
      t30 = 1/t29
      t31 = s1vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t36 = t30*t33*t35
      tncoff8 = -3.D0/8192.D0*t1*t24*t26*t36-3.D0/8192.D0*t1*t26/(-1+l8)
     #*t36+3.D0/4096.D0*t1*t33*t30*t26/l3*t35
      t8 = s23r2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = 3.D0/8192.D0*(-1+z)/l2/(-1+l8)/l8/t9/t13*t15
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = 3.D0/4096.D0*(-1+z)*(l3-1)/l2/l3/t8/t12*t14
      t7 = s23r2(z,l1,l2,l3,l4)
      t8 = t7**2
      t11 = s1vr2(z,l1,l2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = 3.D0/4096.D0*(-1+z)*(l3-1)/l2/l3/t8/t12*t14
      t1 = -1+z
      t3 = dlog(1-l1)
      t5 = dlog(l4)
      t8 = dlog(-t1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(l1)
      t20 = dlog(z)
      t22 = dlog(2.D0)
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
      tncoff12 = -t1*(-6*t3-6*t5-24*t8-6*t13+6*t16-6*t18+12*t20-24*t22-6
     #*t25+6)*t37/16384-3.D0/8192.D0*t1/l8*t31*t34*t36-3.D0/4096.D0*t1/l
     #3*t37-3.D0/8192.D0*t1/l2*t37
      t1 = -1+z
      t3 = dlog(1-l1)
      t5 = dlog(l4)
      t8 = dlog(-t1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t18 = dlog(l1)
      t20 = dlog(z)
      t22 = dlog(2.D0)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff13 = -t1*(-6*t3-6*t5-24*t8-6*t13+6*t16-6*t18+12*t20-24*t22-6
     #*t25+6)*t37/16384+3.D0/8192.D0*t1/(-1+l8)*t37-3.D0/4096.D0*t1/l3*t
     #37-3.D0/8192.D0*t1/l2*t37
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff14 = -3.D0/8192.D0*(-1+z)/l8/(-1+l8)/t8/t11*t14
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff15 = -3.D0/4096.D0*(-1+z)*(l3-1)/l3/t7/t10*t13
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff16 = -3.D0/4096.D0*(-1+z)*(l3-1)/l3/t7/t10*t13
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff17 = 3.D0/8192.D0*(-1+z)/l2/t5/t8*t11
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff18 = 3.D0/8192.D0*(-1+z)/l2/t5/t8*t11
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff19 = -3.D0/8192.D0*(-1+z)/t3/t7*t9
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff20 = -3.D0/8192.D0*(-1+z)/t3/t7*t9

CCC Numerator functions

      tcnum1 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qgrmnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qgrmnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qgrmnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,999.D0,Tcoeff*tncoff9)
      tcnum10 = qgrmnum11e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff10)
      tcnum11 = qgrmnum11e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,999.D0,Tcoeff*tncoff11)
      tcnum12 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff12)
      tcnum13 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff13)
      tcnum14 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff14)
      tcnum15 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff15)
      tcnum16 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff16)
      tcnum17 = qgrmnum11e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff17)
      tcnum18 = qgrmnum11e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,999.D0,Tcoeff*tncoff18)
      tcnum19 = qgrmnum11e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff19)
      tcnum20 = qgrmnum11e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,999.D0,Tcoeff*tncoff20)

CCC Sum

      qgrm14e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm14e1int(var,ndim,Tcoeff)

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
      double precision qgrmnum11e0
      external qgrmnum11e0
      double precision qgrmnum11e1
      external qgrmnum11e1
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

      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(l4)
      t9 = dlog(l1)
      t11 = dlog(-t1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l1)
      t17 = dlog(z)
      t19 = dlog(2.D0)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = 3.D0/8192.D0*t1*(t3+t7+t8+t9+4*t11-1-t14+t16-2*t17+4*t19
     #)*t31-3.D0/8192.D0*t1/l8*t31-3.D0/4096.D0*t1/l3*t31-3.D0/8192.D0*t
     #1/l2*t31
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = f3(0.D0,z)
      t7 = dlog(1-l1*t4)
      t8 = dlog(l4)
      t9 = dlog(l1)
      t11 = dlog(-t1)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(1-l1)
      t17 = dlog(z)
      t19 = dlog(2.D0)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff2 = 3.D0/8192.D0*t1*(t3+t7+t8+t9+4*t11-1-t14+t16-2*t17+4*t19
     #)*t31+3.D0/8192.D0*t1/(-1+l8)*t31-3.D0/4096.D0*t1/l3*t31-3.D0/8192
     #.D0*t1/l2*t31
      t7 = s23r2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,l1,0.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff3 = -3.D0/8192.D0*(-1+z)/l8/(-1+l8)/t8/t11*t14
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff4 = -3.D0/4096.D0*(-1+z)*(l3-1)/l3/t7/t10*t13
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff5 = -3.D0/4096.D0*(-1+z)*(l3-1)/l3/t7/t10*t13
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff6 = 3.D0/8192.D0*(-1+z)/l2/t5/t8*t11
      t4 = s23r2(z,l1,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(z,l1,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff7 = 3.D0/8192.D0*(-1+z)/l2/t5/t8*t11
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -3.D0/8192.D0*(-1+z)/t3/t7*t9
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff9 = -3.D0/8192.D0*(-1+z)/t3/t7*t9

CCC Numerator functions

      tcnum1 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff2)
      tcnum3 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff3)
      tcnum4 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff4)
      tcnum5 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff5)
      tcnum6 = qgrmnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff6)
      tcnum7 = qgrmnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,999.D0,Tcoeff*tncoff7)
      tcnum8 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff8)
      tcnum9 = qgrmnum11e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff9)

CCC Sum

      qgrm14e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm14e2int(var,ndim,Tcoeff)

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
      double precision qgrmnum11e0
      external qgrmnum11e0
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

      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -3.D0/8192.D0*(-1+z)/t3/t7*t9
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff2 = -3.D0/8192.D0*(-1+z)/t3/t7*t9

CCC Numerator functions

      tcnum1 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff1)
      tcnum2 = qgrmnum11e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,999.D0,Tcoeff*tncoff2)

CCC Sum

      qgrm14e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgrm14e3int(var,ndim,Tcoeff)

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

      qgrm14e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

