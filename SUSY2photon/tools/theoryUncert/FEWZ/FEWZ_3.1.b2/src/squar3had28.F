CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq38e0int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
      double precision sq3num4e1
      external sq3num4e1
      double precision sq3num4e2
      external sq3num4e2
      double precision sq3num4e3
      external sq3num4e3
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

      t2 = dlog(1-z)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = t4**2
      t9 = dlog(1-l4)
      t10 = dlog(l4)
      t11 = t9*t10
      t12 = dlog(2.D0)
      t16 = dlog(z)
      t17 = t4*t16
      t20 = t9**2
      t23 = t10**2
      t26 = t2**2
      t29 = t16**2
      t36 = 0.3141592653589793D1**2
      t41 = t12**2
      t48 = -t2*t5/1024-t11*t12/512-z3/768-t17*t12/256-t20*t9/12288-t23*
     #t10/12288-t26*t2/192+t29*t16/1536+t5*t4/12288+t26*t16/128+t10*t36/
     #6144-t10*t29/1024-t10*t41/256+t23*t4/4096+t9*t36/6144
      t81 = -t9*t41/256+t2*t36/1536-t2*t41/64-t9*t23/4096-t9*t26/256-t9*
     #t5/4096-t9*t29/1024+t20*t16/2048-t26*t12/64+t23*t16/2048+t16*t41/1
     #28+t4*t41/256-t29*t12/256-t36*t16/3072+t26*t4/256-t5*t12/1024
      t103 = t10*t16
      t106 = t10*t4
      t113 = t10*t2
      t116 = t9*t16
      t119 = t5*t16/2048-t23*t12/1024+t20*t4/4096-t20*t12/1024+t4*t29/10
     #24-t4*t36/6144-t23*t2/1024-t20*t10/4096-t20*t2/1024+t36*t12/1536+t
     #103*t12/256-t106*t16/1024+t106*t12/512-t10*t5/4096+t113*t4/512+t11
     #6*t12/256
      t120 = t2*t16
      t123 = t2*t4
      t126 = t9*t2
      t145 = t9*t4
      t156 = t120*t12/64+t123*t12/128-t126*t12/128+t126*t4/512-t123*t16/
     #256+t126*t16/256-t10*t26/256-t41*t12/192+t11*t16/1024+t11*t4/2048-
     #t11*t2/512-t145*t16/1024+t145*t12/512-t113*t12/128+t113*t16/256-t2
     #*t29/256
      t159 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t160 = 1/t159
      t162 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t163 = 1/t162
      t164 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t165 = 1/t164
      t167 = f1(1.D0,z,0.D0)
      t168 = t163*t165*t167
      t170 = 2*t36
      t172 = 24*t10*t12
      t173 = 1-l3
      t174 = dlog(t173)
      t175 = t174**2
      t177 = 48*t26
      t178 = 3*t5
      t179 = 12*t29
      t182 = 48*t120
      t183 = 24*t123
      t184 = 24*t126
      t185 = 6*t145
      t187 = 24*t9*t12
      t188 = 6*t11
      t190 = 24*t4*t12
      t191 = -t170+t172+12*t175+t177+t178+t179+12*t9*t174-t182-t183+t184
     #-t185+t187+t188-t190
      t193 = 48*t16*t12
      t195 = 96*t2*t12
      t198 = 12*t116
      t199 = 24*t113
      t200 = 6*t106
      t201 = 48*t41
      t202 = 3*t23
      t209 = 3*t20
      t212 = 12*t17
      t213 = 12*t103
      t214 = -t193+t195-24*t174*t16-t198+t199-t200+t201+t202+48*t174*t12
     #-12*t174*t4+12*t10*t174+t209+48*t174*t2+t212-t213
      t217 = -1/t173
      t222 = -t170+t172+t177+t178+t179-t182-t183+t184-t185+t187+t188-t19
     #0-t193+t195
      t223 = dlog(l2)
      t224 = t223**2
      t238 = -t198+t199-t200+t201+t202+t209+3*t224-6*t223*t4+24*t223*t12
     #-12*t223*t16+6*t10*t223+24*t223*t2+6*t9*t223+t212-t213
      t240 = 1/l2
      t241 = (t222+t238)*t240
      t242 = Dsquar322(z,0.D0,l2,1.D0,l4)
      t243 = 1/t242
      t249 = t165*t167
      t251 = 4*t2
      t252 = 4*t12
      t253 = 2*t16
      t254 = 2*t174
      t256 = (-t10-t9-t251+t4-t252+t253-t223-t254)*t240
      t259 = t217*t160
      t265 = dlog(l1)
      t278 = t265**2
      t280 = 12*t10*t265+12*t9*t265-24*t265*t16+48*t265*t12-12*t265*t4+4
     #8*t265*t2-t170+12*t278+t172+t177+t178+t179-t182-t183
      t281 = t184-t185+t187+t188-t190-t193+t195-t198+t199-t200+t201+t202
     #+t209+t212-t213
      t283 = 1/l1
      t288 = 2*t265
      t297 = (-t10-t9-t223-t288-t251+t4-t252+t253)*t283
      tncoff1 = (t48+t81+t119+t156)*t160*t168-(t191+t214)*t217*t160*t168
     #/6144+(-t241*t243/12288+t241*t160/12288)*t163*t249+(t256*t217*t243
     #/1024-t256*t259/1024)*t163*t249+(t280+t281)*t283*t160*t168/6144-(-
     #t9-t288-t251+t4-t252+t253-t254-t10)*t217*t283*t160*t163*t249/512+(
     #-t297*t240*t243/1024+t297*t240*t160/1024)*t163*t249+(t243*t283*t24
     #0*t217/512-t283*t240*t259/512)*t163*t249
      t1 = dlog(l4)
      t2 = t1**2
      t5 = dlog(1-l4)
      t6 = dlog(l3)
      t10 = dlog(1-z)
      t11 = f2(l3,z)
      t12 = dlog(t11)
      t17 = dlog(z)
      t24 = dlog(2.D0)
      t31 = 1-l3
      t32 = dlog(t31)
      t49 = 3*t2+12*t5*t6-24*t10*t12+6*t5*t1-24*t6*t17+48*t10*t6-6*t5*t1
     #2+96*t10*t24-12*t6*t12+24*t5*t10+48*t32*t10+12*t5*t32+24*t32*t6+48
     #*t6*t24+12*t1*t32-24*t12*t24+12*t1*t6-12*t32*t12
      t66 = 0.3141592653589793D1**2
      t68 = t24**2
      t70 = t6**2
      t72 = t32**2
      t74 = t12**2
      t82 = t17**2
      t84 = t10**2
      t86 = t5**2
      t88 = -48*t10*t17+24*t5*t24+48*t32*t24-12*t1*t17+24*t1*t24+24*t1*t
     #10+12*t12*t17-6*t1*t12-2*t66+48*t68+12*t70+12*t72+3*t74-12*t5*t17-
     #48*t17*t24-24*t32*t17+12*t82+48*t84+3*t86
      t91 = -1/t31
      t93 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t94 = 1/t93
      t96 = s2vr2(z,0.D0,0.D0,l3,l4)
      t97 = 1/t96
      t98 = s23r2(z,0.D0,0.D0,l3,l4)
      t99 = 1/t98
      t101 = f1(l3,z,0.D0)
      t105 = 4*t24
      t106 = 2*t6
      t107 = 4*t10
      t108 = 2*t17
      t109 = dlog(l2)
      t110 = 2*t32
      t112 = 1/l2
      t113 = (-t1-t5-t105-t106-t107+t12+t108-t109-t110)*t112
      t114 = Dsquar322(z,0.D0,l2,l3,l4)
      t115 = 1/t114
      t118 = t91*t94
      t123 = t99*t101
      t125 = dlog(l1)
      t129 = 1/l1
      tncoff2 = (t49+t88)*t91*t94*t97*t99*t101/6144+(-t113*t91*t115/1024
     #+t113*t118/1024)*t97*t123+(-t5+t12-t107-t106-t105-2*t125+t108-t110
     #-t1)*t91*t129*t94*t97*t123/512+(-t115*t129*t112*t91/512+t129*t112*
     #t118/512)*t97*t123
      t1 = dlog(l4)
      t2 = f2(1.D0,z)
      t3 = dlog(t2)
      t7 = dlog(1-z)
      t8 = dlog(2.D0)
      t14 = dlog(1-l4)
      t19 = dlog(z)
      t24 = f3(1.D0,z)
      t27 = dlog(1-t24*l1)
      t30 = dlog(l1)
      t37 = -6*t1*t3+96*t7*t8-24*t3*t8+6*t14*t1+24*t14*t7-12*t1*t19+24*t
     #1*t8+24*t27*t8-24*t30*t19+12*t3*t19-24*t7*t3
      t42 = 0.3141592653589793D1**2
      t44 = t8**2
      t47 = dlog(1-l1)
      t48 = t47**2
      t50 = t27**2
      t64 = 24*t14*t8-48*t7*t19-2*t42+48*t44+3*t48+3*t50+24*t27*t7-6*t27
     #*t3+48*t30*t7+48*t30*t8-12*t30*t3-6*t14*t3
      t74 = t3**2
      t80 = t19**2
      t86 = t14**2
      t88 = -48*t19*t8-12*t14*t19+6*t27*t47+6*t14*t47+3*t74+12*t30*t47+2
     #4*t1*t7+12*t80-12*t27*t19+12*t30*t27+3*t86
      t89 = t7**2
      t97 = t1**2
      t101 = t30**2
      t113 = 48*t89+6*t1*t47+6*t14*t27+6*t1*t27+3*t97+24*t47*t8+12*t101+
     #24*t7*t47+12*t1*t30+12*t14*t30-12*t47*t19-6*t47*t3
      t116 = 1/l1
      t118 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t119 = 1/t118
      t121 = s2vr2(z,l1,0.D0,1.D0,l4)
      t122 = 1/t121
      t123 = s23r2(z,l1,0.D0,1.D0,l4)
      t124 = 1/t123
      t126 = f1(1.D0,z,l1)
      t130 = 2*t30
      t131 = 2*t19
      t132 = 4*t7
      t133 = 4*t8
      t134 = 1-l3
      t135 = dlog(t134)
      t139 = -1/t134
      t143 = t124*t126
      t144 = t119*t122*t143
      t147 = dlog(l2)
      t150 = 1/l2
      tncoff3 = -(t37+t64+t88+t113)*t116*t119*t122*t124*t126/6144+(-t14-
     #t130+t3+t131-t132-t27-t47-t133-2*t135-t1)*t139*t116*t144/512-(-t1-
     #t14-t147-t130+t3+t131-t132-t27-t47-t133)*t116*t150*t144/1024+t116*
     #t122*t119*t150*t139*t143/512
      t2 = dlog(1-l4)
      t3 = dlog(l3)
      t6 = dlog(1-l1)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = f3(l3,z)
      t12 = dlog(1-l1*t9)
      t14 = dlog(1-z)
      t16 = dlog(z)
      t18 = dlog(l1)
      t20 = dlog(2.D0)
      t22 = 1-l3
      t23 = dlog(t22)
      t25 = dlog(l4)
      t26 = -t2-2*t3-t6+t8-t12-4*t14+2*t16-2*t18-4*t20-2*t23-t25
      t28 = -1/t22
      t30 = 1/l1
      t32 = Dsquar322(z,l1,0.D0,l3,l4)
      t33 = 1/t32
      t34 = s2vr2(z,l1,0.D0,l3,l4)
      t35 = 1/t34
      t37 = s23r2(z,l1,0.D0,l3,l4)
      t39 = f1(l3,z,l1)
      t40 = 1/t37*t39
      tncoff4 = -t26*t28*t30*t33*t35*t40/512-t30*t35*t33/l2*t28*t40/512
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l2)
      t5 = l1*l2
      t7 = dlog(1-t5)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(1-l1)
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t17 = dlog(1-z)
      t19 = dlog(l1)
      t21 = f3(1.D0,z)
      t24 = dlog(1-t21*l1)
      t25 = -t1-t3-t4-t7+t9-t11+2*t12-4*t14-4*t17-2*t19-t24
      t26 = 1/l1
      t28 = 1/l2
      t30 = Dsquar322(z,l1,l2,1.D0,l4)
      t31 = 1/t30
      t32 = s2vr2(z,l1,t5,1.D0,l4)
      t33 = 1/t32
      t35 = s23r2(z,l1,t5,1.D0,l4)
      t36 = 1/t35
      t37 = f1(1.D0,z,l1)
      tncoff5 = t25*t26*t28*t31*t33*t36*t37/1024-t33*t36*t31*t26*t28/(l3
     #-1)*t37/512
      t1 = l1*l2
      t2 = s2vr2(z,l1,t1,l3,l4)
      t4 = s23r2(z,l1,t1,l3,l4)
      t7 = Dsquar322(z,l1,l2,l3,l4)
      t15 = f1(l3,z,l1)
      tncoff6 = 1/t2/t4/t7/l1/l2/(l3-1)*t15/512
      t2 = dlog(1-l4)
      t4 = dlog(1-z)
      t7 = dlog(l4)
      t10 = dlog(z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = dlog(2.D0)
      t32 = t2*t4/512+t2*t7/2048-t2*t10/1024-t2*t14/2048+t2*t17/512-t4*t
     #14/512+t4*t17/128-t4*t10/256+t7*t4/512-t7*t10/1024-t14*t7/2048
      t41 = 0.3141592653589793D1**2
      t43 = t17**2
      t45 = t10**2
      t47 = t14**2
      t49 = t4**2
      t51 = t2**2
      t53 = t7**2
      t55 = t7*t17/512-t10*t17/256+t14*t10/1024-t14*t17/512-t41/6144+t43
     #/256+t45/1024+t47/4096+t49/256+t51/4096+t53/4096
      t57 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t61 = 1/t60
      t62 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t63 = 1/t62
      t65 = f1(1.D0,z,0.D0)
      t66 = t61*t63*t65
      t68 = 24*t17
      t69 = 24*t4
      t70 = 6*t14
      t71 = 12*t10
      t72 = 6*t7
      t73 = 6*t2
      t74 = 1-l3
      t75 = dlog(t74)
      t79 = -1/t74
      t84 = dlog(l2)
      t87 = 1/l2
      t88 = (-t68-t69+t70+t71-t72-t73-6*t84)*t87
      t89 = Dsquar322(z,0.D0,l2,1.D0,l4)
      t90 = 1/t89
      t96 = t63*t65
      t98 = t87*t79
      t105 = dlog(l1)
      t108 = 1/l1
      t117 = t108*t87
      tncoff7 = (t32+t55)*t58*t66-(-t68-t69+t70+t71-t72-t73-12*t75)*t79*
     #t58*t66/6144+(-t88*t90/12288+t88*t58/12288)*t61*t96+(t98*t90/1024-
     #t98*t58/1024)*t61*t96+(-t68-t69+t70+t71-t72-t73-12*t105)*t108*t58*
     #t66/6144-t79*t108*t58*t66/512+(-t117*t90/1024+t117*t58/1024)*t61*t
     #96
      t2 = dlog(1-l4)
      t4 = 1-l3
      t5 = dlog(t4)
      t8 = dlog(1-z)
      t10 = dlog(l4)
      t12 = f2(l3,z)
      t13 = dlog(t12)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = dlog(l3)
      t23 = -1/t4
      t25 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t26 = 1/t25
      t28 = s2vr2(z,0.D0,0.D0,l3,l4)
      t29 = 1/t28
      t30 = s23r2(z,0.D0,0.D0,l3,l4)
      t31 = 1/t30
      t33 = f1(l3,z,0.D0)
      t34 = t29*t31*t33
      t38 = 1/l2*t23
      t39 = Dsquar322(z,0.D0,l2,l3,l4)
      tncoff8 = (-6*t2-12*t5-24*t8-6*t10+6*t13-24*t15+12*t17-12*t19)*t23
     #*t26*t34/6144+(-t38/t39/1024+t38*t26/1024)*t29*t31*t33+t23/l1*t26*
     #t34/512
      t1 = dlog(z)
      t4 = dlog(1-l4)
      t6 = f3(1.D0,z)
      t9 = dlog(1-t6*l1)
      t11 = dlog(l4)
      t13 = dlog(l1)
      t16 = dlog(1-l1)
      t19 = dlog(1-z)
      t21 = f2(1.D0,z)
      t22 = dlog(t21)
      t24 = dlog(2.D0)
      t27 = 1/l1
      t29 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,l1,0.D0,1.D0,l4)
      t34 = s23r2(z,l1,0.D0,1.D0,l4)
      t37 = f1(1.D0,z,l1)
      t38 = 1/t32/t34*t37
      tncoff9 = -(12*t1-6*t4-6*t9-6*t11-12*t13-6*t16-24*t19+6*t22-24*t24
     #)*t27*t30*t38/6144+1/(l3-1)*t27*t30*t38/512-t27/l2*t30*t38/1024
      t5 = Dsquar322(z,l1,0.D0,l3,l4)
      t8 = s2vr2(z,l1,0.D0,l3,l4)
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t13 = f1(l3,z,l1)
      tncoff10 = -1/(l3-1)/l1/t5/t8/t10*t13/512
      t4 = Dsquar322(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,l1,t7,1.D0,l4)
      t10 = s23r2(z,l1,t7,1.D0,l4)
      t13 = f1(1.D0,z,l1)
      tncoff11 = 1/l1/l2/t4/t8/t10*t13/1024
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t6 = dlog(1-l4)
      t8 = dlog(l4)
      t11 = dlog(1-z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t21 = 1/t20
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t23 = 1/t22
      t25 = f1(1.D0,z,0.D0)
      t26 = t21*t23*t25
      t33 = 1/l2
      t34 = Dsquar322(z,0.D0,l2,1.D0,l4)
      tncoff12 = (-t1/512+t3/1024-t6/2048-t8/2048-t11/512+t14/2048)*t18*
     #t26-1/(l3-1)*t18*t26/1024+(-t33/t34/2048+t33*t18/2048)*t21*t23*t25
     #+1/l1*t18*t26/1024
      t3 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t6 = s2vr2(z,0.D0,0.D0,l3,l4)
      t8 = s23r2(z,0.D0,0.D0,l3,l4)
      t11 = f1(l3,z,0.D0)
      tncoff13 = 1/(l3-1)/t3/t6/t8*t11/1024
      t2 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,l1,0.D0,1.D0,l4)
      t10 = f1(1.D0,z,l1)
      tncoff14 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff15 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num4e0(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num4e0(var,ndim,z,ux,l1,l1*l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num4e0(var,ndim,z,ux,l1,l1*l2,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num4e1(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num4e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num4e1(var,ndim,z,ux,l1,l1*l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num4e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num4e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num4e2(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num4e3(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff15)

CCC Sum

      sq38e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq38e1int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
      double precision sq3num4e1
      external sq3num4e1
      double precision sq3num4e2
      external sq3num4e2
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

      t2 = dlog(1-l4)
      t4 = dlog(1-z)
      t7 = dlog(l4)
      t10 = dlog(z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = dlog(2.D0)
      t32 = t2*t4/512+t2*t7/2048-t2*t10/1024-t2*t14/2048+t2*t17/512-t4*t
     #14/512+t4*t17/128-t4*t10/256+t7*t4/512-t7*t10/1024-t7*t14/2048
      t41 = 0.3141592653589793D1**2
      t43 = t17**2
      t45 = t10**2
      t47 = t14**2
      t49 = t4**2
      t51 = t2**2
      t53 = t7**2
      t55 = t7*t17/512-t10*t17/256+t14*t10/1024-t14*t17/512-t41/6144+t43
     #/256+t45/1024+t47/4096+t49/256+t51/4096+t53/4096
      t57 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t61 = 1/t60
      t62 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t63 = 1/t62
      t65 = f1(1.D0,z,0.D0)
      t66 = t61*t63*t65
      t68 = 1-l3
      t69 = dlog(t68)
      t71 = 4*t4
      t72 = 4*t17
      t73 = 2*t10
      t76 = -1/t68
      t81 = dlog(l2)
      t83 = 1/l2
      t84 = (-t71+t14-t72+t73-t81-t7-t2)*t83
      t85 = Dsquar322(z,0.D0,l2,1.D0,l4)
      t86 = 1/t85
      t92 = t63*t65
      t94 = t83*t76
      t101 = dlog(l1)
      t104 = 1/l1
      t113 = t104*t83
      tncoff1 = (t32+t55)*t58*t66-(-2*t69-t2-t7-t71+t14-t72+t73)*t76*t58
     #*t66/1024+(-t84*t86/2048+t84*t58/2048)*t61*t92+(t94*t86/1024-t94*t
     #58/1024)*t61*t92+(-t71-t72+t73+t14-2*t101-t7-t2)*t104*t58*t66/1024
     #-t76*t104*t58*t66/512+(-t113*t86/1024+t113*t58/1024)*t61*t92
      t1 = 1-l3
      t2 = dlog(t1)
      t5 = dlog(1-l4)
      t6 = dlog(l4)
      t8 = dlog(1-z)
      t10 = dlog(l3)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t16 = f2(l3,z)
      t17 = dlog(t16)
      t20 = -1/t1
      t22 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t23 = 1/t22
      t25 = s2vr2(z,0.D0,0.D0,l3,l4)
      t26 = 1/t25
      t27 = s23r2(z,0.D0,0.D0,l3,l4)
      t28 = 1/t27
      t30 = f1(l3,z,0.D0)
      t31 = t26*t28*t30
      t35 = 1/l2*t20
      t36 = Dsquar322(z,0.D0,l2,l3,l4)
      tncoff2 = (-2*t2-t5-t6-4*t8-2*t10-4*t12+2*t14+t17)*t20*t23*t31/102
     #4+(-t35/t36/1024+t35*t23/1024)*t26*t28*t30+t20/l1*t23*t31/512
      t1 = f3(1.D0,z)
      t4 = dlog(1-t1*l1)
      t6 = dlog(1-z)
      t9 = dlog(1-l1)
      t10 = dlog(2.D0)
      t12 = dlog(z)
      t14 = f2(1.D0,z)
      t15 = dlog(t14)
      t16 = dlog(l1)
      t18 = dlog(l4)
      t20 = dlog(1-l4)
      t22 = 1/l1
      t24 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(z,l1,0.D0,1.D0,l4)
      t29 = s23r2(z,l1,0.D0,1.D0,l4)
      t32 = f1(1.D0,z,l1)
      t33 = 1/t27/t29*t32
      tncoff3 = -(-t4-4*t6-t9-4*t10+2*t12+t15-2*t16-t18-t20)*t22*t25*t33
     #/1024+1/(l3-1)*t22*t25*t33/512-t22/l2*t25*t33/1024
      t5 = Dsquar322(z,l1,0.D0,l3,l4)
      t8 = s2vr2(z,l1,0.D0,l3,l4)
      t10 = s23r2(z,l1,0.D0,l3,l4)
      t13 = f1(l3,z,l1)
      tncoff4 = -1/(l3-1)/l1/t5/t8/t10*t13/512
      t4 = Dsquar322(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,l1,t7,1.D0,l4)
      t10 = s23r2(z,l1,t7,1.D0,l4)
      t13 = f1(1.D0,z,l1)
      tncoff5 = 1/l1/l2/t4/t8/t10*t13/1024
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t6 = dlog(1-l4)
      t8 = dlog(l4)
      t11 = dlog(1-z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t21 = 1/t20
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t23 = 1/t22
      t25 = f1(1.D0,z,0.D0)
      t26 = t21*t23*t25
      t33 = 1/l2
      t34 = Dsquar322(z,0.D0,l2,1.D0,l4)
      tncoff6 = (-t1/512+t3/1024-t6/2048-t8/2048-t11/512+t14/2048)*t18*t
     #26-1/(l3-1)*t18*t26/1024+(-t33/t34/2048+t33*t18/2048)*t21*t23*t25+
     #1/l1*t18*t26/1024
      t3 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t6 = s2vr2(z,0.D0,0.D0,l3,l4)
      t8 = s23r2(z,0.D0,0.D0,l3,l4)
      t11 = f1(l3,z,0.D0)
      tncoff7 = 1/(l3-1)/t3/t6/t8*t11/1024
      t2 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,l1,0.D0,1.D0,l4)
      t10 = f1(1.D0,z,l1)
      tncoff8 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff9 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num4e0(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num4e0(var,ndim,z,ux,l1,l1*l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num4e1(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num4e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      sq38e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq38e2int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
      double precision sq3num4e1
      external sq3num4e1
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

      t1 = dlog(2.D0)
      t3 = dlog(z)
      t6 = dlog(1-l4)
      t8 = dlog(l4)
      t11 = dlog(1-z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t21 = 1/t20
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t23 = 1/t22
      t25 = f1(1.D0,z,0.D0)
      t26 = t21*t23*t25
      t33 = 1/l2
      t34 = Dsquar322(z,0.D0,l2,1.D0,l4)
      tncoff1 = (-t1/512+t3/1024-t6/2048-t8/2048-t11/512+t14/2048)*t18*t
     #26-1/(l3-1)*t18*t26/1024+(-t33/t34/2048+t33*t18/2048)*t21*t23*t25+
     #1/l1*t18*t26/1024
      t3 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t6 = s2vr2(z,0.D0,0.D0,l3,l4)
      t8 = s23r2(z,0.D0,0.D0,l3,l4)
      t11 = f1(l3,z,0.D0)
      tncoff2 = 1/(l3-1)/t3/t6/t8*t11/1024
      t2 = Dsquar322(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,l1,0.D0,1.D0,l4)
      t7 = s23r2(z,l1,0.D0,1.D0,l4)
      t10 = f1(1.D0,z,l1)
      tncoff3 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff4 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num4e0(var,ndim,z,ux,l1,0.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq38e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq38e3int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
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

      t1 = Dsquar322(z,0.D0,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff1 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq38e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq38e4int(var,ndim,Tcoeff)

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

      sq38e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

