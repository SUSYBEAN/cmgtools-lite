CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm10e0int(var,ndim,Tcoeff)

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
      t3 = dlog(-t1)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t8 = dlog(l1)
      t11 = dlog(2.D0)
      t12 = f3(0.D0,z)
      t15 = dlog(1-l1*t12)
      t21 = dlog(1-l4)
      t26 = dlog(z)
      t33 = dlog(l4)
      t37 = dlog(1-l1)
      t44 = 18-24*t3*t5+24*t3*t8+18*t11*t15+24*t3*t15+6*t8*t21-6*t5*t8-1
     #2*t26*t15-12*t26*t8+18*t11*t8+6*t8*t33-12*t26*t37+24*t3*t37-6*t5*t
     #37
      t73 = 6*t15*t8+6*t15*t21+18*t11*t33+6*t15*t33-6*t5*t21+6*t15*t37+6
     #*t37*t33-6*t15*t5-12*t26*t21+6*t33*t21-48*t3*t26+24*t3*t33-12*t26*
     #t33+24*t3*t21
      t94 = t11**2
      t96 = 0.3141592653589793D1**2
      t99 = -18*t11*t5+12*t26*t5+18*t11*t21+6*t37*t21-6*t5*t33+6*t37*t8+
     #72*t11*t3-36*t11*t26+18*t37*t11-48*t3+27*t94-t96-36*t11-12*t8
      t106 = t5**2
      t108 = t15**2
      t110 = t37**2
      t112 = t8**2
      t114 = t21**2
      t116 = t33**2
      t118 = t26**2
      t120 = t3**2
      t122 = -12*t33-12*t21+12*t5+24*t26-12*t37-12*t15+3*t106+3*t108+3*t
     #110+3*t112+3*t114+3*t116+12*t118+48*t120
      t126 = s23r2(z,l1,0.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = f1(0.D0,z,l1)
      t133 = dlog(l8)
      t147 = t133*l8+t8*l8+t37*l8-t133+t15*l8-t5*l8+3*t11*l8-2*t26*l8+4*
     #t3*l8-2*l8+t21*l8
      t149 = 2*t26
      t150 = 3*t11
      t151 = 4*t3
      t152 = 2+t33*l8+t149-t150-t151-t15-t8-t33-t21-t37+t5
      t155 = 1/l8
      t163 = dlog(l3)
      t165 = t21+t33+2*t163+t8+t15+t151-2+t150-t149-t5+t37
      t167 = 1/l3
      t169 = t167*t128*t129
      t177 = dlog(l2)
      t178 = t33+t21+t15+t8+t37-t5-2+t150-t149+t151+t177
      t180 = 1/l2
      t182 = t180*t128*t129
      tncoff1 = -3.D0/131072.D0*t1*(t44+t73+t99+t122)*t128*t129+9.D0/655
     #36.D0*t1*(t147+t152)*t155/(-1+l8)*t128*t129+9.D0/32768.D0*t1*t165*
     #t169-9.D0/32768.D0*t1*t167*t155*t128*t129+9.D0/65536.D0*t178*t1*t1
     #82-9.D0/65536.D0*t1*t155*t182-9.D0/32768.D0*t1*t180*t169
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(l4)
      t6 = dlog(2.D0)
      t8 = 0.3141592653589793D1**2
      t9 = t6**2
      t12 = dlog(-t1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t20 = dlog(1-l4)
      t22 = dlog(z)
      t25 = dlog(1-l1)
      t27 = f3(0.D0,z)
      t30 = dlog(1-l1*t27)
      t34 = 18+6*t2*t3-36*t6-t8+27*t9-48*t12-12*t2+12*t16-12*t3-12*t20+2
     #4*t22-12*t25-12*t30+6*t30*t2
      t47 = t30**2
      t49 = t16**2
      t51 = t25**2
      t55 = t2**2
      t57 = t20**2
      t63 = -6*t16*t25+24*t12*t25-12*t22*t30-48*t12*t22+6*t3*t20-12*t22*
     #t20+3*t47+3*t49+3*t51+18*t6*t20+3*t55+3*t57+18*t6*t2+6*t30*t3
      t67 = t3**2
      t71 = t22**2
      t77 = t12**2
      t93 = -6*t16*t2+3*t67+6*t25*t3+12*t71+72*t6*t12-36*t6*t22+48*t77+2
     #4*t12*t30+18*t6*t30-12*t22*t25+6*t30*t20+18*t25*t6+6*t30*t25+24*t1
     #2*t2
      t122 = -6*t16*t20+6*t25*t20-6*t16*t3+6*t25*t2+6*t2*t20-24*t12*t16+
     #18*t6*t3+12*t22*t16-18*t6*t16-6*t30*t16-12*t22*t2+24*t12*t20-12*t2
     #2*t3+24*t12*t3
      t126 = s23r2(z,l1,0.D0,0.D0,l4)
      t127 = t126**2
      t128 = 1/t127
      t129 = f1(0.D0,z,l1)
      t133 = 1-l8
      t134 = dlog(t133)
      t149 = -t134*l8+2*l8+2*t22*l8-t25*l8-t20*l8-t3*l8+t16*l8-t2*l8-t30
     #*l8-4*t12*l8-3*t6*l8
      t154 = -1/t133
      t156 = t154*t128*t129
      t159 = dlog(l3)
      t161 = 4*t12
      t162 = 3*t6
      t163 = 2*t22
      t164 = t20+t3+2*t159+t2+t30+t161-2+t162-t163-t16+t25
      t166 = 1/l3
      t168 = t166*t128*t129
      t174 = dlog(l2)
      t175 = t3+t20+t30+t2+t25-t16-2+t162-t163+t161+t174
      t177 = 1/l2
      t179 = t177*t128*t129
      tncoff2 = -3.D0/131072.D0*t1*(t34+t63+t93+t122)*t128*t129+9.D0/655
     #36.D0*t1*t149/l8*t156+9.D0/32768.D0*t1*t164*t168+9.D0/32768.D0*t1*
     #t166*t156+9.D0/65536.D0*t175*t1*t179+9.D0/65536.D0*t1*t154*t179-9.
     #D0/32768.D0*t1*t177*t168
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t5 = dlog(1-l4)
      t6 = dlog(l4)
      t7 = dlog(l8)
      t8 = dlog(l1)
      t10 = dlog(1-l1)
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(-t1)
      t24 = -2+t3+t5+t6+t7+t8+t10-2*t11+3*t13-t16+t20+4*t22
      t26 = 1/l8
      t29 = -1/t2
      t30 = s23r2(z,l1,0.D0,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = f1(0.D0,z,l1)
      tncoff3 = 9.D0/65536.D0*t1*t24*t26*t29*t32*t34-9.D0/32768.D0*t1/l3
     #*t29*t26*t32*t34-9.D0/65536.D0*t1*t26*t29/l2*t32*t34
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(1-l1)
      t9 = f2(l3,z)
      t10 = dlog(t9)
      t12 = dlog(-t2)
      t14 = f3(l3,z)
      t17 = dlog(1-l1*t14)
      t18 = dlog(l1)
      t20 = dlog(1-l4)
      t21 = dlog(l4)
      t22 = dlog(2.D0)
      t24 = dlog(l3)
      t26 = dlog(z)
      t28 = 4*t5+t8-t10+2*t12+t17-2+t18+t20+t21+3*t22+2*t24-2*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(l3,z,l1)
      t36 = t30*t33*t35
      t38 = 1-l8
      tncoff4 = 9.D0/32768.D0*t3*t28*t36-9.D0/32768.D0*t3*t30/l8*t33*t35
     #-9.D0/32768.D0*t3/l2*t36
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t5 = dlog(-t1)
      t8 = dlog(1-l1)
      t9 = f2(l3,z)
      t10 = dlog(t9)
      t12 = dlog(-t2)
      t14 = f3(l3,z)
      t17 = dlog(1-l1*t14)
      t18 = dlog(l1)
      t20 = dlog(1-l4)
      t21 = dlog(l4)
      t22 = dlog(2.D0)
      t24 = dlog(l3)
      t26 = dlog(z)
      t28 = 4*t5+t8-t10+2*t12+t17-2+t18+t20+t21+3*t22+2*t24-2*t26
      t30 = 1/l3
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(l3,z,l1)
      t36 = t30*t33*t35
      tncoff5 = 9.D0/32768.D0*t3*t28*t36+9.D0/32768.D0*t3*t30/(-1+l8)*t3
     #3*t35-9.D0/32768.D0*t3/l2*t36
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t11 = t10**2
      t13 = f1(l3,z,l1)
      tncoff6 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/(-1+l8)/l8/t11*t13
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = dlog(l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(-t1)
      t18 = dlog(1-l1)
      t19 = dlog(z)
      t21 = 2-l2
      t22 = dlog(t21)
      t23 = dlog(l2)
      t24 = t2+t4-2+t5-t7+2*t8+t13+4*t15+t18-2*t19+t22+t23
      t27 = -1/t21
      t29 = 1/l2
      t31 = s23r2(z,l1,l2/2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t39 = 1-l8
      tncoff7 = 9.D0/32768.D0*t1*t24*t27*t29*t33*t35-9.D0/32768.D0*t1*t2
     #7*t29/l8*t33*t35-9.D0/16384.D0*t1*t27*t33*t29/l3*t35
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = dlog(l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(-t1)
      t18 = dlog(1-l1)
      t19 = dlog(z)
      t21 = 2-l2
      t22 = dlog(t21)
      t23 = dlog(l2)
      t24 = t2+t4-2+t5-t7+2*t8+t13+4*t15+t18-2*t19+t22+t23
      t27 = -1/t21
      t29 = 1/l2
      t31 = s23r2(z,l1,l2/2,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = f1(0.D0,z,l1)
      t39 = t1*t27
      tncoff8 = 9.D0/32768.D0*t1*t24*t27*t29*t33*t35+9.D0/32768.D0*t39*t
     #29/(-1+l8)*t33*t35-9.D0/16384.D0*t39*t33*t29/l3*t35
      t12 = s23r2(z,l1,l2/2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff9 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/(-1+l8)/l8/t13*t15
      t11 = s23r2(z,l1,l2/2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff10 = -9.D0/16384.D0*(-1+z)*(l3-1)/(-2+l2)/l2/l3/t12*t14
      t11 = s23r2(z,l1,l2/2,l3,l4)
      t12 = t11**2
      t14 = f1(l3,z,l1)
      tncoff11 = -9.D0/16384.D0*(-1+z)*(l3-1)/(-2+l2)/l2/l3/t12*t14
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t7 = dlog(l4)
      t9 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(2.D0)
      t17 = dlog(-t1)
      t19 = dlog(z)
      t22 = dlog(1-l1)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = f1(0.D0,z,l1)
      t33 = t31*t32
      t36 = 1-l8
      tncoff12 = -3.D0/131072.D0*t1*(12-6*t5-6*t7-6*t9+6*t12-18*t14-24*t
     #17+12*t19-6*t22-6*t25)*t33-9.D0/65536.D0*t1/l8*t31*t32-9.D0/32768.
     #D0*t1/l3*t33-9.D0/65536.D0*t1/l2*t33
      t1 = -1+z
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t7 = dlog(l4)
      t9 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(2.D0)
      t17 = dlog(-t1)
      t19 = dlog(z)
      t22 = dlog(1-l1)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = f1(0.D0,z,l1)
      t33 = 1/t30*t32
      tncoff13 = -3.D0/131072.D0*t1*(12-6*t5-6*t7-6*t9+6*t12-18*t14-24*t
     #17+12*t19-6*t22-6*t25)*t33+9.D0/65536.D0*t1/(-1+l8)*t33-9.D0/32768
     #.D0*t1/l3*t33-9.D0/65536.D0*t1/l2*t33
      t6 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t10 = f1(0.D0,z,l1)
      tncoff14 = -9.D0/65536.D0*(-1+z)/l8/(-1+l8)/t7*t10
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff15 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff16 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t7 = s23r2(z,l1,l2/2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff17 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t8*t11
      t7 = s23r2(z,l1,l2/2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff18 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t8*t11
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff19 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff20 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ggrmnum3e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ggrmnum3e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ggrmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ggrmnum3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      ggrm10e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm10e1int(var,ndim,Tcoeff)

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
      t2 = dlog(l1)
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t8 = dlog(-t1)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-l1)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = f1(0.D0,z,l1)
      t27 = 1/t24*t26
      tncoff1 = 9.D0/65536.D0*t1*(t2+t6+4*t8-2+3*t10-2*t12-t15+t17+t19+t
     #20)*t27-9.D0/65536.D0*t1/l8*t27-9.D0/32768.D0*t1/l3*t27-9.D0/65536
     #.D0*t1/l2*t27
      t1 = -1+z
      t2 = dlog(l1)
      t3 = f3(0.D0,z)
      t6 = dlog(1-l1*t3)
      t8 = dlog(-t1)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-l1)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = f1(0.D0,z,l1)
      t27 = 1/t24*t26
      tncoff2 = 9.D0/65536.D0*t1*(t2+t6+4*t8-2+3*t10-2*t12-t15+t17+t19+t
     #20)*t27+9.D0/65536.D0*t1/(-1+l8)*t27-9.D0/32768.D0*t1/l3*t27-9.D0/
     #65536.D0*t1/l2*t27
      t6 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t10 = f1(0.D0,z,l1)
      tncoff3 = -9.D0/65536.D0*(-1+z)/l8/(-1+l8)/t7*t10
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff4 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t6 = t5**2
      t9 = f1(l3,z,l1)
      tncoff5 = -9.D0/32768.D0*(-1+z)*(l3-1)/l3/t6*t9
      t7 = s23r2(z,l1,l2/2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff6 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t8*t11
      t7 = s23r2(z,l1,l2/2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff7 = -9.D0/32768.D0*(-1+z)/(-2+l2)/l2/t8*t11
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff8 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff9 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      ggrm10e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm10e2int(var,ndim,Tcoeff)

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

      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff1 = -9.D0/65536.D0*(-1+z)/t3*t6
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = f1(0.D0,z,l1)
      tncoff2 = -9.D0/65536.D0*(-1+z)/t3*t6

CCC Numerator functions

      tcnum1 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm10e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

