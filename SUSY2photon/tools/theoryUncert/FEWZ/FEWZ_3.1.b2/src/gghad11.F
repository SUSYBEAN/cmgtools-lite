CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm11e0int(var,ndim,Tcoeff)

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
      double precision ggrmnum3e0
      external ggrmnum3e0
      double precision ggrmnum3e1
      external ggrmnum3e1
      double precision ggrmnum3e2
      external ggrmnum3e2
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
      t2 = dlog(2.D0)
      t3 = dlog(l4)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = dlog(z)
      t14 = dlog(1-l1)
      t18 = dlog(1-l4)
      t21 = f3(0.D0,z)
      t24 = dlog(1-l1*t21)
      t27 = dlog(l1)
      t33 = dlog(-t1)
      t44 = 18+18*t2*t3-18*t2*t7-12*t10*t3+18*t14*t2-12*t10*t18+18*t2*t2
     #4-6*t7*t27+6*t24*t27+72*t2*t33+6*t27*t18+24*t33*t18-12*t10*t24-12*
     #t10*t14
      t73 = 6*t14*t18+24*t33*t14+24*t33*t3-24*t33*t7+6*t14*t27-6*t7*t14+
     #6*t24*t14+6*t24*t18+12*t10*t7+24*t33*t24+6*t27*t3-48*t10*t33+18*t2
     #*t27-12*t10*t27
      t83 = t7**2
      t85 = t24**2
      t87 = t14**2
      t89 = t27**2
      t91 = t18**2
      t93 = t3**2
      t95 = t10**2
      t97 = t33**2
      t100 = 0.3141592653589793D1**2
      t101 = -36*t2*t10-6*t24*t7-6*t7*t18+24*t33*t27+3*t83+3*t85+3*t87+3
     #*t89+3*t91+3*t93+12*t95+48*t97-36*t2-t100
      t102 = t2**2
      t122 = 27*t102-12*t18-12*t3+24*t10-12*t27-12*t14+18*t2*t18+6*t14*t
     #3+6*t24*t3-6*t7*t3+6*t3*t18-48*t33-12*t24+12*t7
      t126 = s23r2(z,l1,1.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = f1(0.D0,z,l1)
      t133 = 2*t10
      t134 = 4*t33
      t135 = 3*t2
      t137 = dlog(l8)
      t138 = t133-t134-t135-t24-t3+t7-t27-t18-t14-2*l8-t137
      t152 = 2+t14*l8+t27*l8-t7*l8+t3*l8+t18*l8+t24*l8+t137*l8+4*t33*l8-
     #2*t10*l8+3*t2*l8
      t155 = 1/l8
      t163 = dlog(l3)
      t165 = t18+t14-t7+t24+t135-2-t133+t27+t134+2*t163+t3
      t167 = 1/l3
      t169 = t167*t128*t129
      t174 = t155*t128*t129
      t177 = dlog(l2)
      t178 = -t133+t27+t14+t24-t7+t135+t134-2+t177+t18+t3
      t180 = 1/l2
      t185 = t1*t180
      tncoff1 = -3.D0/131072.D0*t1*(t44+t73+t101+t122)*t128*t129+9.D0/65
     #536.D0*t1*(t138+t152)*t155/(-1+l8)*t128*t129+9.D0/32768.D0*t1*t165
     #*t169-9.D0/32768.D0*t1*t167*t174+9.D0/65536.D0*t178*t1*t180*t128*t
     #129-9.D0/65536.D0*t185*t174-9.D0/32768.D0*t185*t169
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(l4)
      t10 = dlog(1-l1)
      t13 = dlog(2.D0)
      t15 = dlog(1-l4)
      t20 = dlog(-t1)
      t23 = dlog(l1)
      t26 = f2(0.D0,z)
      t27 = dlog(t26)
      t32 = dlog(z)
      t42 = 18+6*t5*t6+6*t10*t6+18*t13*t15-12*t10+72*t13*t20+6*t5*t23-6*
     #t27*t23+18*t13*t5-12*t32*t15+18*t10*t13-12*t32*t6-18*t13*t27-48*t2
     #0
      t71 = 18*t13*t6+6*t6*t15-6*t27*t6+24*t20*t10+6*t10*t15-12*t32*t5-1
     #2*t32*t10+24*t20*t15+6*t23*t15+24*t20*t23-6*t27*t15-6*t5*t27-36*t1
     #3*t32-12*t32*t23
      t97 = t5**2
      t99 = -12*t5+18*t13*t23-48*t32*t20+6*t23*t6+24*t20*t5+12*t32*t27+6
     #*t5*t15+6*t5*t10-6*t27*t10+6*t10*t23-24*t20*t27+24*t20*t6+12*t27+3
     #*t97
      t100 = t27**2
      t102 = t23**2
      t104 = t10**2
      t106 = t15**2
      t108 = t32**2
      t110 = t6**2
      t112 = t20**2
      t114 = t13**2
      t116 = 0.3141592653589793D1**2
      t122 = 3*t100+3*t102+3*t104+3*t106+12*t108+3*t110+48*t112+27*t114-
     #t116-36*t13-12*t15-12*t6+24*t32-12*t23
      t126 = s23r2(z,l1,1.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = f1(0.D0,z,l1)
      t134 = 1-l8
      t135 = dlog(t134)
      t149 = 2*l8-t135*l8-3*t13*l8+t27*l8-4*t20*l8-t15*l8-t23*l8-t5*l8-t
     #10*l8+2*t32*l8-t6*l8
      t154 = -1/t134
      t156 = t154*t128*t129
      t159 = 3*t13
      t160 = 2*t32
      t161 = 4*t20
      t162 = dlog(l3)
      t164 = t15+t10-t27+t5+t159-2-t160+t23+t161+2*t162+t6
      t166 = 1/l3
      t168 = t166*t128*t129
      t174 = dlog(l2)
      t175 = -t160+t23+t10+t5-t27+t159+t161-2+t174+t15+t6
      t177 = 1/l2
      t179 = t177*t128*t129
      tncoff2 = -3.D0/131072.D0*t1*(t42+t71+t99+t122)*t128*t129+9.D0/655
     #36.D0*t1*t149/l8*t156+9.D0/32768.D0*t1*t164*t168+9.D0/32768.D0*t1*
     #t166*t156+9.D0/65536.D0*t175*t1*t179+9.D0/65536.D0*t1*t154*t179-9.
     #D0/32768.D0*t1*t177*t168
      t1 = -1+z
      t2 = dlog(2.D0)
      t5 = dlog(1-l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t14 = dlog(-t1)
      t16 = dlog(z)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t21 = 1-l8
      t22 = dlog(t21)
      t23 = dlog(l8)
      t24 = 3*t2+t5-t7+t8+t12+4*t14-2*t16+t19+t20+t22+t23-2
      t26 = 1/l8
      t29 = -1/t21
      t30 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      t35 = t29*t32*t34
      tncoff3 = 9.D0/65536.D0*t1*t24*t26*t35-9.D0/32768.D0*t1/l3*t26*t35
     #-9.D0/65536.D0*t1*t29/l2*t26*t32*t34
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = f2(l3,z)
      t5 = dlog(t4)
      t7 = dlog(-t1)
      t10 = dlog(-t2)
      t12 = dlog(2.D0)
      t15 = dlog(1-l1)
      t16 = dlog(l1)
      t17 = dlog(l4)
      t18 = dlog(z)
      t20 = dlog(l3)
      t23 = dlog(1-l4)
      t24 = f3(l3,z)
      t27 = dlog(1-l1*t24)
      t28 = -t5+4*t7+2*t10+3*t12+t15+t16+t17-2*t18+2*t20-2+t23+t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(l3,z,l1)
      t36 = t30*t33*t35
      tncoff4 = 9.D0/32768.D0*t3*t28*t36-9.D0/32768.D0*t3*t30/l8*t33*t35
     #-9.D0/32768.D0*t3/l2*t36
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = f2(l3,z)
      t5 = dlog(t4)
      t7 = dlog(-t1)
      t10 = dlog(-t2)
      t12 = dlog(2.D0)
      t15 = dlog(1-l1)
      t16 = dlog(l1)
      t17 = dlog(l4)
      t18 = dlog(z)
      t20 = dlog(l3)
      t23 = dlog(1-l4)
      t24 = f3(l3,z)
      t27 = dlog(1-l1*t24)
      t28 = -t5+4*t7+2*t10+3*t12+t15+t16+t17-2*t18+2*t20-2+t23+t27
      t30 = 1/l3
      t31 = s23r2(z,l1,1.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(l3,z,l1)
      t36 = t30*t33*t35
      tncoff5 = 9.D0/32768.D0*t3*t28*t36+9.D0/32768.D0*t3*t30/(-1+l8)*t3
     #3*t35-9.D0/32768.D0*t3/l2*t36
      t10 = s23r2(z,l1,1.D0,l3,l4)
      t11 = t10**2
      t13 = f1(l3,z,l1)
      tncoff6 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/l8/(-1+l8)/t11*t13
      t1 = -1+z
      t2 = 2-l2
      t3 = dlog(t2)
      t5 = dlog(-t1)
      t7 = dlog(z)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t14 = dlog(1-l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t18 = dlog(2.D0)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+4*t5-2*t7+t12+t14-t16-2+t17+2*t18+t20+t22+t23
      t27 = -1/t2
      t29 = 1/l2
      t32 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t40 = 1-l8
      tncoff7 = 9.D0/32768.D0*t1*t24*t27*t29*t34*t36-9.D0/32768.D0*t1*t2
     #7*t29/l8*t34*t36-9.D0/16384.D0*t1*t27*t34*t29/l3*t36
      t1 = -1+z
      t2 = 2-l2
      t3 = dlog(t2)
      t5 = dlog(-t1)
      t7 = dlog(z)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t14 = dlog(1-l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = dlog(l1)
      t18 = dlog(2.D0)
      t20 = dlog(l2)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t24 = t3+4*t5-2*t7+t12+t14-t16-2+t17+2*t18+t20+t22+t23
      t27 = -1/t2
      t29 = 1/l2
      t32 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = f1(0.D0,z,l1)
      t40 = t1*t27
      tncoff8 = 9.D0/32768.D0*t1*t24*t27*t29*t34*t36+9.D0/32768.D0*t40*t
     #29/(-1+l8)*t34*t36-9.D0/16384.D0*t40*t34*t29/l3*t36
      t13 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t14 = t13**2
      t16 = f1(0.D0,z,l1)
      tncoff9 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/l8/(-1+l8)/t14*t16
      t12 = s23r2(z,l1,1-l2/2,l3,l4)
      t13 = t12**2
      t15 = f1(l3,z,l1)
      tncoff10 = -9.D0/16384.D0*(-1+z)*(l3-1)/(-2+l2)/l2/l3/t13*t15
      t12 = s23r2(z,l1,1-l2/2,l3,l4)
      t13 = t12**2
      t15 = f1(l3,z,l1)
      tncoff11 = -9.D0/16384.D0*(-1+z)*(l3-1)/(-2+l2)/l2/l3/t13*t15
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t5 = dlog(2.D0)
      t7 = dlog(z)
      t10 = dlog(1-l1)
      t12 = dlog(l1)
      t14 = dlog(l4)
      t17 = dlog(-t1)
      t19 = f3(0.D0,z)
      t22 = dlog(1-l1*t19)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = f1(0.D0,z,l1)
      t33 = t31*t32
      t36 = 1-l8
      tncoff12 = -3.D0/131072.D0*t1*(6*t3-18*t5+12*t7-6*t10-6*t12-6*t14-
     #24*t17-6*t22-6*t25+12)*t33-9.D0/65536.D0*t1/l8*t31*t32-9.D0/32768.
     #D0*t1/l3*t33-9.D0/65536.D0*t1/l2*t33
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t5 = dlog(2.D0)
      t7 = dlog(z)
      t10 = dlog(1-l1)
      t12 = dlog(l1)
      t14 = dlog(l4)
      t17 = dlog(-t1)
      t19 = f3(0.D0,z)
      t22 = dlog(1-l1*t19)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,1.D0,0.D0,l4)
      t30 = t29**2
      t32 = f1(0.D0,z,l1)
      t33 = 1/t30*t32
      tncoff13 = -3.D0/131072.D0*t1*(6*t3-18*t5+12*t7-6*t10-6*t12-6*t14-
     #24*t17-6*t22-6*t25+12)*t33+9.D0/65536.D0*t1/(-1+l8)*t33-9.D0/32768
     #.D0*t1/l3*t33-9.D0/65536.D0*t1/l2*t33
      t6 = s23r2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t10 = f1(0.D0,z,l1)
      tncoff14 = -9.D0/65536.D0*(-1+z)/l8/(-1+l8)/t7*t10
      t5 = s23r2(z,l1,1.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff15 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t5 = s23r2(z,l1,1.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff16 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t8 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff17 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t9*t12
      t8 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff18 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t9*t12
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff19 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff20 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ggrmnum3e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ggrmnum3e1(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ggrmnum3e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ggrmnum3e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      ggrm11e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm11e1int(var,ndim,Tcoeff)

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
      double precision ggrmnum3e0
      external ggrmnum3e0
      double precision ggrmnum3e1
      external ggrmnum3e1
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
      t5 = dlog(1-l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(l4)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t19 = dlog(-t1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = f1(0.D0,z,l1)
      t27 = 1/t24*t26
      tncoff1 = 9.D0/65536.D0*t1*(t3+t5-t7+t8+t12-2+t13+3*t14-2*t16+4*t1
     #9)*t27-9.D0/65536.D0*t1/l8*t27-9.D0/32768.D0*t1/l3*t27-9.D0/65536.
     #D0*t1/l2*t27
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(1-l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = f3(0.D0,z)
      t12 = dlog(1-l1*t9)
      t13 = dlog(l4)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t19 = dlog(-t1)
      t23 = s23r2(z,l1,1.D0,0.D0,l4)
      t24 = t23**2
      t26 = f1(0.D0,z,l1)
      t27 = 1/t24*t26
      tncoff2 = 9.D0/65536.D0*t1*(t3+t5-t7+t8+t12-2+t13+3*t14-2*t16+4*t1
     #9)*t27+9.D0/65536.D0*t1/(-1+l8)*t27-9.D0/32768.D0*t1/l3*t27-9.D0/6
     #5536.D0*t1/l2*t27
      t6 = s23r2(z,l1,1.D0,0.D0,l4)
      t7 = t6**2
      t10 = f1(0.D0,z,l1)
      tncoff3 = -9.D0/65536.D0*(-1+z)/l8/(-1+l8)/t7*t10
      t5 = s23r2(z,l1,1.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff4 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t5 = s23r2(z,l1,1.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff5 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t8 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff6 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t9*t12
      t8 = s23r2(z,l1,1-l2/2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,l1)
      tncoff7 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t9*t12
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff8 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff9 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum3e0(var,ndim,z,ux,l1,1-l2/2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum3e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      ggrm11e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm11e2int(var,ndim,Tcoeff)

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
      double precision ggrmnum3e0
      external ggrmnum3e0
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

      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff1 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,1.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff2 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm11e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

