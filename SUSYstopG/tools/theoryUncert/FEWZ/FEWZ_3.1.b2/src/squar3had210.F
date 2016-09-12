CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq310e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
      double precision sq3num5e2
      external sq3num5e2
      double precision sq3num5e3
      external sq3num5e3
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
      t2 = 1/t1
      t3 = 1/l2
      t5 = 1/l3
      t6 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t7 = 1/t6
      t8 = t5*t7
      t10 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = 1/t10
      t12 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t13 = 1/t12
      t15 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t16 = 1/t15
      t17 = f1(1.D0,1.D0,1.D0)
      t19 = t11*t13*t16*t17
      t22 = dlog(l1)
      t23 = dlog(2.D0)
      t24 = t22*t23
      t25 = f3(1.D0,1.D0)
      t27 = dlog(1-t25)
      t30 = f2(1.D0,1.D0)
      t31 = dlog(t30)
      t32 = t22*t31
      t35 = dlog(l4)
      t36 = t35*t23
      t41 = t35*t31
      t44 = t35*t22
      t48 = dlog(1-l4)
      t49 = t48*t35
      t52 = t27*t31
      t55 = t22*t27
      t59 = 0.3141592653589793D1**2
      t62 = t35**2
      t65 = t48**2
      t70 = t23**2
      t73 = t24*t27/4096-t32*t23/4096+t36*t22/4096+t36*t27/4096-t41*t23/
     #4096-t44*t31/8192+t49*t27/8192-t52*t23/4096-t55*t31/8192+z3/3072-t
     #59*t23/12288-t62*t31/16384+t65*t22/16384-t65*t31/16384+t48*t70/409
     #6
      t74 = t27**2
      t79 = t22**2
      t82 = t31**2
      t105 = t48*t23
      t108 = t48*t27
      t111 = t48*t74/16384+t62*t27/16384+t79*t27/16384+t48*t82/16384+t48
     #*t79/16384+t35*t70/4096+t35*t79/16384+t35*t82/16384+t62*t22/16384+
     #t22*t74/16384-t22*t59/24576+t22*t82/16384-t27*t59/24576+t70*t23/61
     #44+t105*t27/4096-t108*t31/8192
      t113 = t48*t22
      t120 = t35*t27
      t147 = -t113*t31/8192+t108*t22/8192+t49*t22/8192-t120*t31/8192+t79
     #*t22/49152+t79*t23/8192-t31*t70/4096+t27*t82/16384+t22*t70/4096+t3
     #1*t59/24576+t27*t70/4096-t74*t31/16384-t48*t59/24576+t48*t62/16384
     #-t79*t31/16384+t74*t23/8192
      t174 = t48*t31
      t181 = t62*t23/8192-t35*t59/24576+t65*t27/16384+t65*t23/8192+t65*t
     #35/16384+t35*t74/16384+t23*t82/8192+t49*t23/4096-t49*t31/8192+t65*
     #t48/49152+t62*t35/49152-t82*t31/49152+t74*t27/49152-t174*t23/4096+
     #t105*t22/4096+t120*t22/8192
      t189 = dlog(l3)
      t190 = 2*t189
      t191 = dlog(l2)
      t192 = 2*t191
      t193 = 2*t23
      t200 = dlog(-t1)
      t201 = 4*t200
      t204 = t2*t7
      t208 = 2*t59
      t209 = 3*t62
      t210 = 3*t65
      t211 = 12*t70
      t212 = 3*t74
      t213 = 3*t79
      t214 = 3*t82
      t215 = 6*t49
      t216 = 6*t113
      t217 = 6*t120
      t222 = 6*t52
      t223 = 12*t105
      t224 = -t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+24*t48*t
     #200+24*t35*t200-t222+t223
      t226 = 12*t23*t27
      t231 = 6*t174
      t232 = 6*t108
      t233 = 6*t44
      t234 = 12*t36
      t235 = 6*t41
      t236 = 6*t32
      t237 = 12*t24
      t238 = 6*t55
      t244 = 12*t23*t31
      t245 = t200**2
      t247 = t226+24*t200*t22+24*t200*t27-t231+t232+t233+t234-t235-t236+
     #t237+t238+48*t200*t23-24*t200*t31-t244+48*t245
      t265 = 12*t191*t22-t208+t209+t210+t211+t212+t213+t214+12*t191*t27+
     #12*t35*t191+12*t48*t191+24*t191*t23-12*t191*t31+t215
      t266 = t191**2
      t268 = t216+t217-t222+t223+t226-t231+t232+t233+t234-t235-t236+t237
     #+t238-t244+12*t266
      t286 = 12*t189*t27+12*t189*t22+12*t48*t189-12*t189*t31+12*t35*t189
     #+24*t189*t23-t208+t209+t210+t211+t212+t213+t214+t215
      t287 = t189**2
      t289 = t216+t217-t222+t223+t226-t231+t232+t233+t234-t235-t236+t237
     #+t238-t244+12*t287
      tncoff1 = t2*t3*t8*t19/512+(t73+t111+t147+t181)*t7*t11*t13*t16*t17
     #+(t35+t48+t190-t31+t192+t27+t193+t22)*t3*t8*t19/2048-(t201+t48+t19
     #3+t27+t22-t31+t190+t35)*t5*t204*t19/1024+(t224+t247)*t2*t7*t19/122
     #88-(t265+t268)*t3*t7*t19/24576-(t286+t289)*t5*t7*t19/24576-(t35+t4
     #8+t201+t22+t27-t31+t193+t192)*t3*t204*t19/1024
      t1 = -2+l3
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s1vr2(1.D0,1.D0,1.D0,t4,l4)
      t6 = 1/t5
      t8 = s2vr2(1.D0,1.D0,1.D0,t4,l4)
      t9 = 1/t8
      t10 = s23r2(1.D0,1.D0,1.D0,t4,l4)
      t11 = 1/t10
      t12 = t9*t11
      t14 = Dsquar311(1.D0,l1,0.D0,t4,l4)
      t15 = 1/t14
      t16 = 1/l2
      t18 = 1/l3
      t19 = -1+z
      t21 = t18/t19
      t22 = f1(t4,1.D0,1.D0)
      t28 = dlog(1-l4)
      t29 = dlog(l1)
      t32 = dlog(l3)
      t35 = dlog(l4)
      t36 = f2(t4,1.D0)
      t37 = dlog(t36)
      t41 = dlog(-t1)
      t46 = f3(t4,1.D0)
      t48 = dlog(1-t46)
      t67 = 6*t28*t29+12*t32*t29-6*t35*t37+12*t29*t41+6*t28*t35+6*t28*t4
     #8+6*t35*t29-6*t48*t37+12*t35*t32-6*t29*t37+24*t32*t41-6*t28*t37+12
     #*t32*t48+12*t28*t32
      t80 = 0.3141592653589793D1**2
      t82 = t41**2
      t86 = t32**2
      t88 = t48**2
      t90 = t37**2
      t92 = t28**2
      t94 = t35**2
      t96 = t29**2
      t98 = 12*t35*t41+12*t28*t41+12*t48*t41+6*t48*t35-12*t37*t41-12*t32
     #*t37-2*t80+12*t82+6*t29*t48+12*t86+3*t88+3*t90+3*t92+3*t94+3*t96
      t108 = dlog(-t19)
      t110 = 2*t41
      t111 = 2*t32
      t118 = t15*t9*t11*t6*t22
      t121 = dlog(l2)
      tncoff2 = t2*t6*t12*t15*t16*t21*t22/256-(t67+t98)*t2*t18*t15*t12*t
     #6*t22/12288-(4*t108+t28+t110+t48-t37+t29+t111+t35)*t2*t21*t118/512
     #+(t35+t28+t111+2*t121+t110-t37+t48+t29)*t2*t16*t18*t118/1024
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l2
      t6 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t7 = 1/t6
      t8 = 1/l3*t7
      t11 = 1-l1*l2
      t12 = 1-l2
      t13 = s2vr2(1.D0,t11,t12,1.D0,l4)
      t15 = s23r2(1.D0,t11,t12,1.D0,l4)
      t18 = s1vr2(1.D0,t11,t12,1.D0,l4)
      t20 = f1(1.D0,1.D0,t11)
      t22 = 1/t13/t15/t18*t20
      t25 = dlog(l2)
      t26 = dlog(2.D0)
      t29 = dlog(l4)
      t32 = 0.3141592653589793D1**2
      t34 = dlog(t12)
      t37 = t26**2
      t39 = dlog(t11)
      t40 = t39**2
      t45 = dlog(1-l4)
      t48 = t34**2
      t50 = f2(1.D0,1.D0)
      t51 = dlog(t50)
      t56 = 24*t25*t26+12*t29*t25-2*t32+6*t29*t34+12*t37+3*t40+6*t29*t39
     #+6*t45*t29+3*t48-6*t39*t51+12*t45*t25
      t59 = dlog(l1)
      t64 = f3(1.D0,1.D0)
      t68 = dlog(1-t64+t64*l1*l2)
      t73 = t25**2
      t87 = 12*t25*t39+12*t25*t59+12*t25*t34+12*t25*t68+6*t45*t59+12*t73
     #+6*t59*t68+12*t39*t26+6*t29*t68-6*t59*t51+6*t68*t39+6*t59*t39
      t111 = -12*t26*t51+6*t45*t39+12*t68*t26+6*t34*t39+12*t34*t26+12*t4
     #5*t26+12*t26*t59-6*t29*t51-12*t25*t51+6*t29*t59+12*t29*t26
      t118 = t45**2
      t120 = t29**2
      t122 = t59**2
      t128 = t68**2
      t134 = t51**2
      t136 = 6*t34*t68+6*t45*t68+6*t45*t34+3*t118+3*t120+3*t122-6*t34*t5
     #1+6*t59*t34+3*t128-6*t68*t51-6*t45*t51+3*t134
      t143 = dlog(l3)
      t145 = 2*t26
      t146 = 2*t25
      t153 = dlog(-t1)
      tncoff3 = -t2*t3*t8*t22/512+(t56+t87+t111+t136)*t3*t7*t22/24576-(t
     #29+t45+2*t143+t68+t59+t39+t34-t51+t145+t146)*t3*t8*t22/2048+(t29+t
     #45+4*t153+t59+t34+t68+t39-t51+t145+t146)*t3*t2*t7*t22/1024
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l3)
      t6 = dlog(l2)
      t8 = dlog(l1)
      t9 = 2-l3
      t10 = dlog(t9)
      t13 = 1-l1*l2
      t14 = dlog(t13)
      t15 = 1-l2
      t16 = dlog(t15)
      t18 = 1-l3/2
      t19 = f3(t18,1.D0)
      t23 = dlog(1-t19+t19*l1*l2)
      t24 = f2(t18,1.D0)
      t25 = dlog(t24)
      t28 = -1/t9
      t30 = 1/l2
      t31 = 1/l3
      t34 = Dsquar311(1.D0,l1,l2,t18,l4)
      t35 = 1/t34
      t36 = s2vr2(1.D0,t13,t15,t18,l4)
      t37 = 1/t36
      t39 = s23r2(1.D0,t13,t15,t18,l4)
      t40 = 1/t39
      t41 = s1vr2(1.D0,t13,t15,t18,l4)
      t42 = 1/t41
      t44 = f1(t18,1.D0,t13)
      tncoff4 = -(t1+t3+2*t4+2*t6+t8+2*t10+t14+t16+t23-t25)*t28*t30*t31*
     #t35*t37*t40*t42*t44/1024-t28*t30*t31*t42*t37*t40*t35/(-1+z)*t44/25
     #6
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = 1-z
      t5 = dlog(t4)
      t6 = 4*t5
      t7 = dlog(l1)
      t8 = f3(1.D0,z)
      t10 = dlog(1-t8)
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t13 = dlog(z)
      t14 = 2*t13
      t15 = dlog(l2)
      t17 = dlog(2.D0)
      t18 = 2*t17
      t20 = 1/l2
      t23 = -1/t4
      t24 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t25 = 1/t24
      t26 = t23*t25
      t28 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t30 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t33 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t35 = f1(1.D0,z,1.D0)
      t37 = 1/t28/t30/t33*t35
      t76 = 6*t1*t7+12*t1*t17+12*t17*t7+48*t5*t17-12*t3*t13+6*t1*t10-24*
     #t13*t17-6*t10*t12-12*t13*t7-12*t1*t13-6*t1*t12-6*t3*t12-48*t5*t13+
     #24*t5*t10+6*t3*t10-12*t10*t13-12*t12*t17-24*t5*t12
      t85 = t10**2
      t87 = t12**2
      t89 = t13**2
      t91 = t5**2
      t93 = t3**2
      t95 = t1**2
      t97 = t7**2
      t111 = 0.3141592653589793D1**2
      t113 = t17**2
      t115 = 6*t7*t10+12*t10*t17+12*t12*t13-6*t7*t12+3*t85+3*t87+12*t89+
     #48*t91+3*t93+3*t95+3*t97+24*t5*t7+6*t3*t1+6*t3*t7+24*t3*t5+24*t1*t
     #5+12*t3*t17-2*t111+12*t113
      t122 = 1/l3
      t127 = dlog(l3)
      tncoff5 = (t1+t3+t6+t7+t10-t12-t14+2*t15+t18)*t20*t26*t37/1024-(t7
     #6+t115)*t23*t25*t37/12288-t23*t20*t122*t25*t37/512+(t6+t3-t12+t7+2
     #*t127-t14+t18+t10+t1)*t122*t26*t37/1024
      t1 = -2+l3
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s1vr2(z,1.D0,1.D0,t4,l4)
      t6 = 1/t5
      t8 = s2vr2(z,1.D0,1.D0,t4,l4)
      t9 = 1/t8
      t10 = s23r2(z,1.D0,1.D0,t4,l4)
      t11 = 1/t10
      t14 = Dsquar311(z,l1,0.D0,t4,l4)
      t15 = 1/t14
      t19 = -1+z
      t21 = 1/l3/t19
      t22 = f1(t4,z,1.D0)
      t28 = dlog(-t19)
      t31 = dlog(1-l4)
      t32 = f2(t4,z)
      t33 = dlog(t32)
      t34 = f3(t4,z)
      t36 = dlog(1-t34)
      t37 = dlog(z)
      t40 = dlog(-t1)
      t42 = dlog(l1)
      t43 = dlog(l3)
      t45 = dlog(l4)
      tncoff6 = -t2*t6*t9*t11*t15/l2*t21*t22/256+(4*t28+t31-t33+t36-2*t3
     #7+2*t40+t42+2*t43+t45)*t2*t21*t15*t9*t11*t6*t22/512
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = 1-z
      t5 = dlog(t4)
      t7 = dlog(2.D0)
      t10 = 1-l1*l2
      t11 = dlog(t10)
      t12 = dlog(z)
      t14 = f3(1.D0,z)
      t18 = dlog(1-t14+t14*l1*l2)
      t19 = 1-l2
      t20 = dlog(t19)
      t21 = dlog(l1)
      t22 = f2(1.D0,z)
      t23 = dlog(t22)
      t24 = dlog(l2)
      t26 = t1+t3+4*t5+2*t7+t11-2*t12+t18+t20+t21-t23+2*t24
      t27 = 1/l2
      t30 = -1/t4
      t31 = Dsquar311(z,l1,l2,1.D0,l4)
      t32 = 1/t31
      t35 = s2vr2(z,t10,t19,1.D0,l4)
      t37 = s23r2(z,t10,t19,1.D0,l4)
      t40 = s1vr2(z,t10,t19,1.D0,l4)
      t42 = f1(1.D0,z,t10)
      t44 = 1/t35/t37/t40*t42
      tncoff7 = -t26*t27*t30*t32*t44/1024+t30*t27/l3*t32*t44/512
      t2 = 1-l3/2
      t4 = 1-l1*l2
      t5 = f1(t2,z,t4)
      t6 = 1-l2
      t7 = s1vr2(z,t4,t6,t2,l4)
      t10 = s2vr2(z,t4,t6,t2,l4)
      t16 = s23r2(z,t4,t6,t2,l4)
      t18 = Dsquar311(z,l1,l2,t2,l4)
      tncoff8 = t5/t7/t10/(-2+l3)/t16/t18/l2/l3/(-1+z)/256
      t1 = 1/l3
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t9 = 1/t8
      t10 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = 1/t10
      t13 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t14 = 1/t13
      t15 = f1(1.D0,1.D0,1.D0)
      t17 = t9*t11*t14*t15
      t20 = 1/l2
      t25 = f3(1.D0,1.D0)
      t27 = dlog(1-t25)
      t28 = 6*t27
      t29 = dlog(l4)
      t30 = 6*t29
      t32 = dlog(1-l4)
      t33 = 6*t32
      t34 = f2(1.D0,1.D0)
      t35 = dlog(t34)
      t36 = 6*t35
      t37 = dlog(2.D0)
      t38 = 12*t37
      t39 = dlog(l1)
      t40 = 6*t39
      t42 = dlog(-t2)
      t49 = dlog(l2)
      t56 = dlog(l3)
      t77 = 0.3141592653589793D1**2
      t79 = t35**2
      t81 = t27**2
      t83 = t37**2
      t85 = t32**2
      t87 = t29**2
      t89 = -t29*t39/8192+t29*t35/8192-t29*t27/8192+t39*t35/8192-t29*t37
     #/4096+t77/24576-t79/16384-t81/16384-t83/4096-t85/16384-t87/16384
      t90 = t39**2
      t112 = -t90/16384-t32*t37/4096-t32*t39/8192-t32*t27/8192-t32*t29/8
     #192+t32*t35/8192+t37*t35/4096-t37*t27/4096+t27*t35/8192-t27*t39/81
     #92-t37*t39/4096
      tncoff9 = t1*t3*t6*t17/1024-t20*t1*t6*t17/2048+(-t28-t30-t33+t36-t
     #38-t40-24*t42)*t3*t6*t17/12288-(-t28-t30-t33+t36-t38-t40-12*t49)*t
     #20*t6*t17/24576-(-t28-t30-t33+t36-t38-t40-12*t56)*t1*t6*t17/24576+
     #t20*t3*t6*t17/1024+(t89+t112)*t6*t9*t11*t14*t15
      t2 = 1-l3/2
      t3 = f3(t2,1.D0)
      t5 = dlog(1-t3)
      t7 = 2-l3
      t8 = dlog(t7)
      t10 = f2(t2,1.D0)
      t11 = dlog(t10)
      t13 = dlog(l1)
      t16 = dlog(1-l4)
      t18 = dlog(l4)
      t20 = dlog(l3)
      t24 = -1/t7
      t26 = 1/l3
      t27 = Dsquar311(1.D0,l1,0.D0,t2,l4)
      t28 = 1/t27
      t29 = t26*t28
      t31 = s2vr2(1.D0,1.D0,1.D0,t2,l4)
      t33 = s23r2(1.D0,1.D0,1.D0,t2,l4)
      t36 = s1vr2(1.D0,1.D0,1.D0,t2,l4)
      t38 = f1(t2,1.D0,1.D0)
      t40 = 1/t31/t33/t36*t38
      tncoff10 = -(-6*t5-12*t8+6*t11-6*t13-6*t16-6*t18-12*t20)*t24*t29*t
     #40/12288+t24*t26/(-1+z)*t28*t40/512-t24/l2*t29*t40/1024
      t1 = dlog(l2)
      t3 = dlog(l4)
      t6 = dlog(1-l4)
      t8 = dlog(2.D0)
      t11 = 1-l1*l2
      t12 = dlog(t11)
      t14 = f2(1.D0,1.D0)
      t15 = dlog(t14)
      t17 = f3(1.D0,1.D0)
      t21 = dlog(1-t17+t17*l1*l2)
      t23 = 1-l2
      t24 = dlog(t23)
      t26 = dlog(l1)
      t29 = 1/l2
      t31 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t32 = 1/t31
      t34 = s2vr2(1.D0,t11,t23,1.D0,l4)
      t36 = s23r2(1.D0,t11,t23,1.D0,l4)
      t39 = s1vr2(1.D0,t11,t23,1.D0,l4)
      t41 = f1(1.D0,1.D0,t11)
      t43 = 1/t34/t36/t39*t41
      tncoff11 = (-12*t1-6*t3-6*t6-12*t8-6*t12+6*t15-6*t21-6*t24-6*t26)*
     #t29*t32*t43/24576+t29/l3*t32*t43/2048-t29/(-1+z)*t32*t43/1024
      t7 = 1-l3/2
      t8 = Dsquar311(1.D0,l1,l2,t7,l4)
      t13 = 1-l1*l2
      t14 = 1-l2
      t15 = s2vr2(1.D0,t13,t14,t7,l4)
      t17 = s23r2(1.D0,t13,t14,t7,l4)
      t20 = s1vr2(1.D0,t13,t14,t7,l4)
      t22 = f1(t7,1.D0,t13)
      tncoff12 = 1/(-2+l3)/l2/l3/t8/t15/t17/t20*t22/1024
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t13 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t15 = f1(1.D0,z,1.D0)
      t17 = 1/t8/t10/t13*t15
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t23 = f3(1.D0,z)
      t25 = dlog(1-t23)
      t27 = dlog(l4)
      t30 = dlog(1-l4)
      t33 = dlog(-t2)
      t35 = dlog(l1)
      t37 = dlog(2.D0)
      t39 = dlog(z)
      tncoff13 = -1/l2*t3*t6*t17/1024-(6*t21-6*t25-6*t27-6*t30-24*t33-6*
     #t35-12*t37+12*t39)*t3*t6*t17/12288-1/l3*t3*t6*t17/1024
      t8 = 1-l3/2
      t9 = Dsquar311(z,l1,0.D0,t8,l4)
      t13 = s2vr2(z,1.D0,1.D0,t8,l4)
      t15 = s23r2(z,1.D0,1.D0,t8,l4)
      t18 = s1vr2(z,1.D0,1.D0,t8,l4)
      t20 = f1(t8,z,1.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t9/t13/t15/t18*t20/512
      t5 = Dsquar311(z,l1,l2,1.D0,l4)
      t9 = 1-l1*l2
      t10 = 1-l2
      t11 = s2vr2(z,t9,t10,1.D0,l4)
      t13 = s23r2(z,t9,t10,1.D0,l4)
      t16 = s1vr2(z,t9,t10,1.D0,l4)
      t18 = f1(1.D0,z,t9)
      tncoff15 = 1/l2/(-1+z)/t5/t11/t13/t16*t18/1024
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f3(1.D0,1.D0)
      t8 = dlog(1-t6)
      t10 = dlog(2.D0)
      t12 = f2(1.D0,1.D0)
      t13 = dlog(t12)
      t15 = dlog(l1)
      t18 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t26 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t29 = f1(1.D0,1.D0,1.D0)
      t30 = 1/t24/t26*t29
      tncoff16 = (t1/8192+t4/8192+t8/8192+t10/4096-t13/8192+t15/8192)*t1
     #9*t22*t30+1/(-1+z)*t19*t22*t30/2048-1/l3*t19*t22*t30/4096-1/l2*t19
     #*t22*t30/4096
      t6 = 1-l3/2
      t7 = Dsquar311(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = s1vr2(1.D0,1.D0,1.D0,t6,l4)
      t17 = f1(t6,1.D0,1.D0)
      tncoff17 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,1.D0,l4)
      t11 = s23r2(1.D0,t6,t7,1.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,1.D0,l4)
      t16 = f1(1.D0,1.D0,t6)
      tncoff18 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t14 = f1(1.D0,z,1.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      tncoff20 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,
     #l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num5e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = sq3num5e2(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = sq3num5e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = sq3num5e3(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      sq310e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12
     #+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff
     #16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnco
     #ff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq310e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
      double precision sq3num5e2
      external sq3num5e2
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

      t1 = dlog(l4)
      t2 = dlog(l1)
      t5 = f2(1.D0,1.D0)
      t6 = dlog(t5)
      t9 = f3(1.D0,1.D0)
      t11 = dlog(1-t9)
      t16 = dlog(2.D0)
      t19 = 0.3141592653589793D1**2
      t21 = t6**2
      t23 = t11**2
      t25 = t16**2
      t28 = dlog(1-l4)
      t29 = t28**2
      t31 = t1**2
      t33 = -t1*t2/8192+t1*t6/8192-t1*t11/8192+t2*t6/8192-t1*t16/4096+t1
     #9/24576-t21/16384-t23/16384-t25/4096-t29/16384-t31/16384
      t34 = t2**2
      t56 = -t34/16384-t28*t16/4096-t28*t2/8192-t28*t11/8192-t28*t1/8192
     #+t28*t6/8192+t16*t6/4096-t16*t11/4096+t11*t6/8192-t11*t2/8192-t16*
     #t2/4096
      t58 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t59 = 1/t58
      t61 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t62 = 1/t61
      t64 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t65 = 1/t64
      t66 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t67 = 1/t66
      t69 = f1(1.D0,1.D0,1.D0)
      t72 = 2*t16
      t73 = 1-z
      t74 = dlog(t73)
      t78 = -1/t73
      t83 = t62*t65*t67*t69
      t86 = dlog(l3)
      t89 = 1/l3
      t98 = dlog(l2)
      t101 = 1/l2
      tncoff1 = (t33+t56)*t59*t62*t65*t67*t69-(t1+t2+t11+t72-t6+t28+4*t7
     #4)*t78*t59*t83/2048+(t2+t11+t72-t6+2*t86+t1+t28)*t89*t59*t83/4096+
     #t89*t78*t59*t83/1024+(t1+t28+2*t98+t2+t11+t72-t6)*t101*t59*t83/409
     #6+t101*t78*t59*t83/1024-t101*t89*t59*t83/2048
      t1 = dlog(l1)
      t3 = 1-l3/2
      t4 = f2(t3,1.D0)
      t5 = dlog(t4)
      t6 = f3(t3,1.D0)
      t8 = dlog(1-t6)
      t9 = 2-l3
      t10 = dlog(t9)
      t12 = dlog(l3)
      t14 = dlog(l4)
      t16 = dlog(1-l4)
      t19 = -1/t9
      t21 = 1/l3
      t22 = Dsquar311(1.D0,l1,0.D0,t3,l4)
      t23 = 1/t22
      t24 = t21*t23
      t26 = s2vr2(1.D0,1.D0,1.D0,t3,l4)
      t28 = s23r2(1.D0,1.D0,1.D0,t3,l4)
      t31 = s1vr2(1.D0,1.D0,1.D0,t3,l4)
      t33 = f1(t3,1.D0,1.D0)
      t35 = 1/t26/t28/t31*t33
      tncoff2 = (t1-t5+t8+2*t10+2*t12+t14+t16)*t19*t24*t35/2048+t19*t21/
     #(-1+z)*t23*t35/512-t19/l2*t24*t35/1024
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t7 = 1-l2
      t8 = dlog(t7)
      t9 = f3(1.D0,1.D0)
      t13 = dlog(1-t9+t9*l1*l2)
      t15 = 1-l1*l2
      t16 = dlog(t15)
      t17 = f2(1.D0,1.D0)
      t18 = dlog(t17)
      t19 = dlog(2.D0)
      t22 = 1/l2
      t24 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(1.D0,t15,t7,1.D0,l4)
      t29 = s23r2(1.D0,t15,t7,1.D0,l4)
      t32 = s1vr2(1.D0,t15,t7,1.D0,l4)
      t34 = f1(1.D0,1.D0,t15)
      t36 = 1/t27/t29/t32*t34
      tncoff3 = -(t1+t3+2*t4+t6+t8+t13+t16-t18+2*t19)*t22*t25*t36/4096-t
     #22/(-1+z)*t25*t36/1024+t22/l3*t25*t36/2048
      t7 = 1-l3/2
      t8 = Dsquar311(1.D0,l1,l2,t7,l4)
      t13 = 1-l1*l2
      t14 = 1-l2
      t15 = s2vr2(1.D0,t13,t14,t7,l4)
      t17 = s23r2(1.D0,t13,t14,t7,l4)
      t20 = s1vr2(1.D0,t13,t14,t7,l4)
      t22 = f1(t7,1.D0,t13)
      tncoff4 = 1/(-2+l3)/l2/l3/t8/t15/t17/t20*t22/1024
      t1 = dlog(l4)
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = f2(1.D0,z)
      t6 = dlog(t5)
      t7 = dlog(z)
      t9 = f3(1.D0,z)
      t11 = dlog(1-t9)
      t13 = dlog(1-l4)
      t14 = 1-z
      t15 = dlog(t14)
      t19 = -1/t14
      t21 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t26 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t29 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t31 = f1(1.D0,z,1.D0)
      t33 = 1/t24/t26/t29*t31
      tncoff5 = (t1+t2+2*t3-t6-2*t7+t11+t13+4*t15)*t19*t22*t33/2048-1/l3
     #*t19*t22*t33/1024-1/l2*t19*t22*t33/1024
      t8 = 1-l3/2
      t9 = Dsquar311(z,l1,0.D0,t8,l4)
      t13 = s2vr2(z,1.D0,1.D0,t8,l4)
      t15 = s23r2(z,1.D0,1.D0,t8,l4)
      t18 = s1vr2(z,1.D0,1.D0,t8,l4)
      t20 = f1(t8,z,1.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t9/t13/t15/t18*t20/512
      t5 = Dsquar311(z,l1,l2,1.D0,l4)
      t9 = 1-l1*l2
      t10 = 1-l2
      t11 = s2vr2(z,t9,t10,1.D0,l4)
      t13 = s23r2(z,t9,t10,1.D0,l4)
      t16 = s1vr2(z,t9,t10,1.D0,l4)
      t18 = f1(1.D0,z,t9)
      tncoff7 = 1/l2/(-1+z)/t5/t11/t13/t16*t18/1024
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f3(1.D0,1.D0)
      t8 = dlog(1-t6)
      t10 = dlog(2.D0)
      t12 = f2(1.D0,1.D0)
      t13 = dlog(t12)
      t15 = dlog(l1)
      t18 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t26 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t29 = f1(1.D0,1.D0,1.D0)
      t30 = 1/t24/t26*t29
      tncoff8 = (t1/8192+t4/8192+t8/8192+t10/4096-t13/8192+t15/8192)*t19
     #*t22*t30+1/(-1+z)*t19*t22*t30/2048-1/l3*t19*t22*t30/4096-1/l2*t19*
     #t22*t30/4096
      t6 = 1-l3/2
      t7 = Dsquar311(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = s1vr2(1.D0,1.D0,1.D0,t6,l4)
      t17 = f1(t6,1.D0,1.D0)
      tncoff9 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,1.D0,l4)
      t11 = s23r2(1.D0,t6,t7,1.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,1.D0,l4)
      t16 = f1(1.D0,1.D0,t6)
      tncoff10 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t14 = f1(1.D0,z,1.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      tncoff12 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num5e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num5e1(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num5e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num5e2(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      sq310e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq310e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num5e0
      external sq3num5e0
      double precision sq3num5e1
      external sq3num5e1
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

      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = f3(1.D0,1.D0)
      t8 = dlog(1-t6)
      t10 = dlog(2.D0)
      t12 = f2(1.D0,1.D0)
      t13 = dlog(t12)
      t15 = dlog(l1)
      t18 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t22 = 1/t21
      t24 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t26 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t29 = f1(1.D0,1.D0,1.D0)
      t30 = 1/t24/t26*t29
      tncoff1 = (t1/8192+t4/8192+t8/8192+t10/4096-t13/8192+t15/8192)*t19
     #*t22*t30+1/(-1+z)*t19*t22*t30/2048-1/l3*t19*t22*t30/4096-1/l2*t19*
     #t22*t30/4096
      t6 = 1-l3/2
      t7 = Dsquar311(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = s1vr2(1.D0,1.D0,1.D0,t6,l4)
      t17 = f1(t6,1.D0,1.D0)
      tncoff2 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar311(1.D0,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(1.D0,t6,t7,1.D0,l4)
      t11 = s23r2(1.D0,t6,t7,1.D0,l4)
      t13 = s1vr2(1.D0,t6,t7,1.D0,l4)
      t16 = f1(1.D0,1.D0,t6)
      tncoff3 = 1/l2/t2/t8/t11/t13*t16/4096
      t3 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t14 = f1(1.D0,z,1.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      tncoff5 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num5e0(var,ndim,1.D0,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num5e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num5e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      sq310e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq310e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num5e0
      external sq3num5e0
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

      t1 = Dsquar311(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s1vr2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      tncoff1 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num5e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq310e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq310e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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

      sq310e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

