CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function rm13e0int(var,ndim,Tcoeff)

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
      double precision rm1num2e0
      external rm1num2e0
      double precision rm1num2e1
      external rm1num2e1
      double precision rm1num2e2
      external rm1num2e2
      double precision rm1num2e3
      external rm1num2e3
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

      t2 = dlog(l4)
      t3 = t2**2
      t4 = dlog(2.D0)
      t8 = dlog(1-z)
      t11 = dlog(z)
      t15 = dlog(1-l4)
      t16 = t15**2
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t25 = t20**2
      t28 = t8**2
      t31 = t11**2
      t40 = 0.3141592653589793D1**2
      t45 = -z3/384-t3*t4/512-t3*t8/512+t3*t11/1024-t16*t8/512+t16*t20/2
     #048+t16*t11/1024-t15*t25/2048-t15*t28/128-t15*t31/512-t2*t25/2048-
     #t2*t28/128-t2*t31/512-t40*t11/1536-t31*t4/128
      t46 = t4**2
      t79 = t11*t46/64-t25*t8/512-t25*t4/512+t25*t11/1024-t28*t4/32+t28*
     #t11/64-t8*t31/128+t8*t40/768-t8*t46/32+t20*t28/128+t20*t31/512-t20
     #*t40/3072+t20*t46/128+t3*t20/2048-t15*t3/2048-t15*t46/128
      t105 = t2*t20
      t112 = t2*t4
      t115 = t15*t40/3072-t16*t2/2048-t16*t4/512-t2*t46/128+t2*t40/3072-
     #t28*t8/96+t31*t11/768-t3*t2/6144+t25*t20/6144-t16*t15/6144-t46*t4/
     #96+t40*t4/768+t105*t8/256+t105*t4/256-t105*t11/512-t112*t8/64
      t116 = t2*t11
      t121 = t8*t11
      t124 = t20*t8
      t129 = t20*t11
      t132 = t15*t20
      t139 = t15*t4
      t142 = t15*t2
      t145 = t15*t11
      t156 = t116*t8/128+t116*t4/128+t121*t4/32+t124*t4/64-t124*t11/128-
     #t129*t4/128+t132*t8/256+t132*t4/256-t132*t11/512-t139*t8/64-t142*t
     #4/256+t145*t8/128+t145*t4/128+t142*t20/1024-t142*t8/256+t142*t11/5
     #12
      t159 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t160 = 1/t159
      t162 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t163 = 1/t162
      t164 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t165 = t164**2
      t166 = 1/t165
      t168 = f1(0.D0,z,0.D0)
      t169 = t163*t166*t168
      t171 = 6*t142
      t172 = 24*t112
      t173 = 24*t139
      t174 = dlog(l3)
      t182 = 96*t8*t4
      t184 = 24*t15*t8
      t186 = 24*t2*t8
      t193 = 6*t132
      t194 = 12*t145
      t195 = t171+t172+t173+48*t174*t4+12*t15*t174+12*t2*t174+t182+t184+
     #t186-12*t174*t20+48*t174*t8-24*t174*t11-t193-t194
      t196 = 6*t105
      t197 = 12*t116
      t198 = 24*t124
      t199 = 48*t121
      t201 = 48*t11*t4
      t203 = 24*t20*t4
      t204 = 12*t129
      t205 = 3*t25
      t206 = 48*t28
      t207 = 3*t16
      t208 = t174**2
      t210 = 3*t3
      t211 = 48*t46
      t212 = 2*t40
      t213 = 12*t31
      t214 = -t196-t197-t198-t199-t201-t203+t204+t205+t206+t207+12*t208+
     #t210+t211-t212+t213
      t216 = 1/l3
      t221 = t171+t172+t173+t182+t184+t186-t193-t194-t196-t197-t198-t199
     #-t201-t203
      t222 = dlog(l2)
      t235 = t222**2
      t237 = t204+24*t222*t4+6*t2*t222+6*t15*t222+24*t8*t222-6*t20*t222-
     #12*t222*t11+3*t235+t205+t206+t207+t210+t211-t212+t213
      t239 = 1/l2
      t244 = 2*t11
      t245 = 4*t8
      t246 = 4*t4
      t247 = 2*t174
      t252 = t166*t168
      t253 = t160*t163*t252
      t256 = dlog(l1)
      t269 = t256**2
      t271 = 24*t256*t4+6*t2*t256+6*t15*t256-6*t256*t20+24*t256*t8-12*t2
     #56*t11+3*t269+t171+t172+t173+t182+t184+t186-t193
      t272 = -t194-t196-t197-t198-t199-t201-t203+t204+t205+t206+t207+t21
     #0+t211-t212+t213
      t274 = 1/l1
      tncoff1 = (t45+t79+t115+t156)*t160*t169+(t195+t214)*t216*t160*t169
     #/3072+(t221+t237)*t239*t160*t169/6144-(-t244+t222+t245-t20+t246+t2
     #+t15+t247)*t239*t216*t253/512+(t271+t272)*t274*t160*t169/6144-(t2+
     #t256-t20+t245+t246-t244+t247+t15)*t216*t274*t253/512-(t15+t2+t222+
     #t256-t20+t245+t246-t244)*t274*t239*t253/1024+t274*t160*t163*t239*t
     #216*t252/512
      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = dlog(1-z)
      t6 = t5**2
      t8 = dlog(z)
      t9 = t8**2
      t11 = 0.3141592653589793D1**2
      t14 = dlog(1-l4)
      t15 = dlog(l4)
      t22 = dlog(l3)
      t40 = dlog(-l3+1)
      t45 = f2(l3,z)
      t46 = dlog(t45)
      t49 = 48*t2+48*t6+12*t9-2*t11+6*t14*t15+24*t15*t1+24*t14*t1+48*t22
     #*t1+12*t14*t22+12*t15*t22+96*t5*t1+24*t14*t5+24*t15*t5+48*t22*t5-2
     #4*t22*t8+48*t40*t1-24*t8*t40-24*t46*t5
      t78 = t15**2
      t80 = t14**2
      t82 = t22**2
      t84 = t40**2
      t86 = t46**2
      t88 = 48*t5*t40+12*t8*t46-12*t46*t40-24*t46*t1-12*t22*t46+24*t22*t
     #40-6*t15*t46+12*t15*t40-6*t14*t46+12*t14*t40-12*t14*t8-12*t15*t8-4
     #8*t5*t8-48*t8*t1+3*t78+3*t80+12*t82+12*t84+3*t86
      t90 = 1/l3
      t92 = s34r2(z,0.D0,0.D0,l3,l4)
      t93 = 1/t92
      t95 = s23r2(z,0.D0,0.D0,l3,l4)
      t96 = 1/t95
      t97 = s1vr2(z,0.D0,0.D0,l3,l4)
      t98 = t97**2
      t99 = 1/t98
      t101 = f1(l3,z,0.D0)
      t105 = 2*t8
      t106 = dlog(l2)
      t107 = 2*t40
      t108 = 4*t1
      t109 = 4*t5
      t110 = 2*t22
      t112 = 1/l2
      t116 = t99*t101
      t117 = t93*t96*t116
      t120 = dlog(l1)
      t123 = 1/l1
      tncoff2 = -(t49+t88)*t90*t93*t96*t99*t101/3072+(-t105-t46+t106+t10
     #7+t108+t109+t15+t14+t110)*t112*t90*t117/512+(t15-t105+t109+t108-t4
     #6+t120+t107+t110+t14)*t90*t123*t117/512-t123*t93*t96*t112*t90*t116
     #/512
      t1 = dlog(z)
      t2 = t1**2
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = t5**2
      t9 = dlog(1-z)
      t10 = t9**2
      t12 = 0.3141592653589793D1**2
      t14 = dlog(2.D0)
      t15 = t14**2
      t17 = dlog(l4)
      t19 = dlog(1-l2)
      t22 = dlog(l2)
      t28 = dlog(1-l4)
      t49 = 12*t2+3*t6+48*t10-2*t12+48*t15+6*t17*t19+6*t17*t22+24*t17*t1
     #4+6*t28*t19+6*t28*t22+24*t28*t14+96*t9*t14+24*t28*t9+24*t17*t9+24*
     #t9*t19+24*t9*t22-6*t28*t5-12*t28*t1
      t80 = t19**2
      t82 = t22**2
      t84 = t17**2
      t86 = t28**2
      t88 = -6*t17*t5-12*t17*t1-24*t9*t5-48*t9*t1-6*t22*t5-12*t22*t1-48*
     #t14*t1-6*t5*t19-24*t14*t5+12*t5*t1-12*t19*t1+6*t28*t17+24*t14*t22+
     #6*t19*t22+24*t14*t19+3*t80+3*t82+3*t84+3*t86
      t90 = 1/l2
      t92 = s34r2(z,0.D0,l2,0.D0,l4)
      t93 = 1/t92
      t95 = s23r2(z,0.D0,l2,0.D0,l4)
      t96 = 1/t95
      t97 = s1vr2(z,0.D0,l2,0.D0,l4)
      t98 = t97**2
      t99 = 1/t98
      t101 = f1(0.D0,z,0.D0)
      t105 = 4*t14
      t106 = 2*t1
      t107 = 4*t9
      t108 = dlog(l3)
      t112 = 1/l3
      t116 = t93*t96*t99*t101
      t119 = dlog(l1)
      t121 = 1/l1
      tncoff3 = -(t49+t88)*t90*t93*t96*t99*t101/6144+(t22+t19+t105-t106-
     #t5+t107+t17+t28+2*t108)*t90*t112*t116/512+(t28+t17+t22+t119+t105+t
     #19-t5-t106+t107)*t121*t90*t116/1024-t93*t121*t96*t99*t90*t112*t101
     #/512
      t1 = f2(l3,z)
      t2 = dlog(t1)
      t4 = dlog(1-l2)
      t6 = dlog(-l3+1)
      t8 = dlog(2.D0)
      t10 = dlog(z)
      t12 = dlog(l2)
      t14 = dlog(1-z)
      t16 = dlog(l4)
      t18 = dlog(1-l4)
      t19 = dlog(l3)
      t22 = 1/l2
      t24 = 1/l3
      t26 = s34r2(z,0.D0,l2,l3,l4)
      t28 = s23r2(z,0.D0,l2,l3,l4)
      t30 = 1/t26/t28
      t31 = s1vr2(z,0.D0,l2,l3,l4)
      t32 = t31**2
      t33 = 1/t32
      t34 = f1(l3,z,0.D0)
      tncoff4 = -(-t2+t4+2*t6+4*t8-2*t10+t12+4*t14+t16+t18+2*t19)*t22*t2
     #4*t30*t33*t34/512+t30*t33/l1*t22*t24*t34/512
      t1 = dlog(l1)
      t3 = dlog(1-l1)
      t6 = dlog(l4)
      t10 = dlog(1-l4)
      t13 = dlog(2.D0)
      t29 = dlog(1-z)
      t32 = 6*t1*t3+6*t6*t3+6*t10*t3+24*t3*t13+6*t10*t6+24*t6*t13+24*t10
     #*t13+24*t1*t13+6*t6*t1+6*t10*t1+96*t29*t13
      t37 = f2(0.D0,z)
      t38 = dlog(t37)
      t43 = dlog(z)
      t50 = f3(0.D0,z)
      t53 = dlog(1-l1*t50)
      t64 = 24*t10*t29+24*t6*t29-6*t1*t38+24*t1*t29-12*t1*t43-12*t3*t43+
     #24*t29*t3+24*t53*t13-12*t53*t43-6*t53*t38+24*t53*t29+6*t53*t3
      t76 = t6**2
      t78 = t10**2
      t80 = t1**2
      t82 = t3**2
      t84 = t53**2
      t86 = t38**2
      t88 = -6*t38*t3+6*t1*t53+6*t6*t53+6*t10*t53-6*t10*t38+3*t76+3*t78+
     #3*t80+3*t82+3*t84+3*t86
      t89 = t43**2
      t91 = 0.3141592653589793D1**2
      t93 = t13**2
      t95 = t29**2
      t113 = 12*t89-2*t91+48*t93+48*t95-12*t10*t43-6*t6*t38-12*t6*t43-24
     #*t29*t38-48*t29*t43-48*t43*t13-24*t38*t13+12*t38*t43
      t116 = 1/l1
      t118 = s34r2(z,l1,0.D0,0.D0,l4)
      t119 = 1/t118
      t121 = s23r2(z,l1,0.D0,0.D0,l4)
      t122 = 1/t121
      t123 = s1vr2(z,l1,0.D0,0.D0,l4)
      t124 = t123**2
      t125 = 1/t124
      t127 = f1(0.D0,z,l1)
      t131 = 4*t29
      t132 = 4*t13
      t133 = 2*t43
      t134 = dlog(l3)
      t137 = 1/l3
      t141 = t125*t127
      t142 = t119*t122*t141
      t145 = dlog(l2)
      t148 = 1/l2
      tncoff5 = -(t32+t64+t88+t113)*t116*t119*t122*t125*t127/6144+(t6+t1
     #+t131+t53-t38+t3+t132-t133+2*t134+t10)*t137*t116*t142/512+(t10+t6+
     #t145+t1+t131+t53-t38+t3+t132-t133)*t116*t148*t142/1024-t116*t148*t
     #119*t122*t137*t141/512
      t1 = dlog(l4)
      t2 = dlog(l1)
      t3 = f3(l3,z)
      t6 = dlog(1-l1*t3)
      t8 = dlog(1-z)
      t11 = dlog(1-l1)
      t13 = dlog(-l3+1)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = f2(l3,z)
      t20 = dlog(t19)
      t21 = dlog(l3)
      t24 = dlog(1-l4)
      t25 = t1+t2+t6+4*t8+t11+2*t13+4*t15-2*t17-t20+2*t21+t24
      t26 = 1/l3
      t28 = 1/l1
      t30 = s34r2(z,l1,0.D0,l3,l4)
      t31 = 1/t30
      t32 = s23r2(z,l1,0.D0,l3,l4)
      t33 = 1/t32
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      tncoff6 = -t25*t26*t28*t31*t33*t39/512+t28*t31*t33/l2*t26*t39/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t5 = dlog(z)
      t8 = dlog(1-z)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l2)
      t16 = dlog(l1)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t20 = dlog(1-l1)
      t21 = dlog(2.D0)
      t23 = t2+t3+t4-2*t5+4*t8+t13+t15+t16-t18+t20+4*t21
      t24 = 1/l1
      t26 = 1/l2
      t28 = s34r2(z,l1,l2,0.D0,l4)
      t29 = 1/t28
      t30 = s23r2(z,l1,l2,0.D0,l4)
      t31 = 1/t30
      t33 = s1vr2(z,l1,l2,0.D0,l4)
      t34 = t33**2
      t35 = 1/t34
      t36 = f1(0.D0,z,l1)
      tncoff7 = -t23*t24*t26*t29*t31*t35*t36/1024+t24*t29*t31*t35*t26/l3
     #*t36/512
      t1 = s34r2(z,l1,l2,l3,l4)
      t3 = s23r2(z,l1,l2,l3,l4)
      t6 = s1vr2(z,l1,l2,l3,l4)
      t7 = t6**2
      t14 = f1(l3,z,l1)
      tncoff8 = -1/t1/t3/t7/l1/l2/l3*t14/512
      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t3 = t2**2
      t6 = dlog(1-z)
      t7 = t6**2
      t9 = dlog(z)
      t10 = t9**2
      t12 = 0.3141592653589793D1**2
      t14 = dlog(2.D0)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t32 = t3/2048+t7/128+t10/512-t12/3072+t15/128+t18/2048+t22/2048-t1
     #7*t2/1024+t17*t6/256+t17*t14/256-t17*t9/512
      t55 = -t6*t2/256-t2*t14/256+t2*t9/512+t6*t14/64-t6*t9/128-t9*t14/1
     #28-t21*t2/1024+t21*t6/256+t21*t14/256-t21*t9/512+t21*t17/1024
      t57 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t62 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t63 = t62**2
      t66 = f1(0.D0,z,0.D0)
      t67 = 1/t60/t63*t66
      t69 = 6*t2
      t70 = 24*t6
      t71 = 24*t14
      t72 = 12*t9
      t73 = 6*t17
      t74 = 6*t21
      t75 = dlog(l3)
      t78 = 1/l3
      t83 = dlog(l2)
      t86 = 1/l2
      t95 = dlog(l1)
      t98 = 1/l1
      tncoff9 = (t32+t55)*t58*t67+(t69-t70-t71+t72-t73-t74-12*t75)*t78*t
     #58*t67/3072+(t69-t70-t71+t72-t73-t74-6*t83)*t86*t58*t67/6144+t86*t
     #78*t58*t67/512+(t69-t70-t71+t72-t73-t74-6*t95)*t98*t58*t67/6144+t7
     #8*t98*t58*t67/512+t98*t86*t58*t67/1024
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = f2(l3,z)
      t9 = dlog(t8)
      t12 = dlog(-l3+1)
      t14 = dlog(l3)
      t16 = dlog(l4)
      t19 = dlog(1-l4)
      t22 = 1/l3
      t24 = s34r2(z,0.D0,0.D0,l3,l4)
      t25 = 1/t24
      t27 = s23r2(z,0.D0,0.D0,l3,l4)
      t29 = s1vr2(z,0.D0,0.D0,l3,l4)
      t30 = t29**2
      t33 = f1(l3,z,0.D0)
      t34 = 1/t27/t30*t33
      tncoff10 = -(-24*t2-24*t4+12*t6+6*t9-12*t12-12*t14-6*t16-6*t19)*t2
     #2*t25*t34/3072-1/l2*t22*t25*t34/512-t22/l1*t25*t34/512
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t6 = dlog(1-l2)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t12 = dlog(1-z)
      t14 = dlog(l2)
      t16 = dlog(l4)
      t19 = dlog(1-l4)
      t22 = 1/l2
      t24 = s34r2(z,0.D0,l2,0.D0,l4)
      t25 = 1/t24
      t27 = s23r2(z,0.D0,l2,0.D0,l4)
      t29 = s1vr2(z,0.D0,l2,0.D0,l4)
      t30 = t29**2
      t33 = f1(0.D0,z,0.D0)
      t34 = 1/t27/t30*t33
      tncoff11 = -(12*t1-24*t3-6*t6+6*t9-24*t12-6*t14-6*t16-6*t19)*t22*t
     #25*t34/6144-t22/l3*t25*t34/512-t22/l1*t25*t34/1024
      t4 = s34r2(z,0.D0,l2,l3,l4)
      t7 = s23r2(z,0.D0,l2,l3,l4)
      t9 = s1vr2(z,0.D0,l2,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,0.D0)
      tncoff12 = 1/l2/l3/t4/t7/t10*t13/512
      t1 = f3(0.D0,z)
      t4 = dlog(1-l1*t1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t13 = dlog(1-l1)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = dlog(l1)
      t21 = dlog(l4)
      t24 = dlog(1-l4)
      t27 = 1/l1
      t29 = s34r2(z,l1,0.D0,0.D0,l4)
      t30 = 1/t29
      t32 = s23r2(z,l1,0.D0,0.D0,l4)
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t35 = t34**2
      t38 = f1(0.D0,z,l1)
      t39 = 1/t32/t35*t38
      tncoff13 = -(-6*t4+6*t7-24*t10-6*t13-24*t15+12*t17-6*t19-6*t21-6*t
     #24)*t27*t30*t39/6144-1/l3*t27*t30*t39/512-t27/l2*t30*t39/1024
      t4 = s34r2(z,l1,0.D0,l3,l4)
      t7 = s23r2(z,l1,0.D0,l3,l4)
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff14 = 1/l3/l1/t4/t7/t10*t13/512
      t4 = s34r2(z,l1,l2,0.D0,l4)
      t7 = s23r2(z,l1,l2,0.D0,l4)
      t9 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,z,l1)
      tncoff15 = 1/l1/l2/t4/t7/t10*t13/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = t22**2
      t26 = f1(0.D0,z,0.D0)
      t27 = 1/t20/t23*t26
      tncoff16 = (-t2/1024-t4/1024+t7/1024-t10/256-t12/256+t14/512)*t18*
     #t27+1/l3*t18*t27/512+1/l2*t18*t27/1024+1/l1*t18*t27/1024
      t2 = s34r2(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,0.D0)
      tncoff17 = -1/l3/t2/t5/t8*t11/512
      t2 = s34r2(z,0.D0,l2,0.D0,l4)
      t5 = s23r2(z,0.D0,l2,0.D0,l4)
      t7 = s1vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,0.D0)
      tncoff18 = -1/l2/t2/t5/t8*t11/1024
      t2 = s34r2(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff19 = -1/l1/t2/t5/t8*t11/1024
      t1 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,0.D0)
      tncoff20 = 1/t1/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num2e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num2e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,-999.D0,-
     #999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num2e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num2e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num2e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = rm1num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = rm1num2e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = rm1num2e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = rm1num2e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = rm1num2e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = rm1num2e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = rm1num2e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = rm1num2e3(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      rm13e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff1
     #6+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncof
     #f20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm13e1int(var,ndim,Tcoeff)

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
      double precision rm1num2e0
      external rm1num2e0
      double precision rm1num2e1
      external rm1num2e1
      double precision rm1num2e2
      external rm1num2e2
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
      t6 = dlog(1-z)
      t7 = t6**2
      t9 = dlog(z)
      t10 = t9**2
      t12 = 0.3141592653589793D1**2
      t14 = dlog(2.D0)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t32 = t3/2048+t7/128+t10/512-t12/3072+t15/128+t18/2048+t22/2048-t1
     #7*t2/1024+t17*t6/256+t17*t14/256-t17*t9/512
      t55 = -t6*t2/256-t2*t14/256+t2*t9/512+t6*t14/64-t6*t9/128-t9*t14/1
     #28-t21*t2/1024+t21*t6/256+t21*t14/256-t21*t9/512+t21*t17/1024
      t57 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t62 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t63 = t62**2
      t66 = f1(0.D0,z,0.D0)
      t67 = 1/t60/t63*t66
      t69 = 4*t6
      t70 = 4*t14
      t71 = 2*t9
      t72 = dlog(l3)
      t75 = 1/l3
      t80 = dlog(l2)
      t82 = 1/l2
      t91 = dlog(l1)
      t93 = 1/l1
      tncoff1 = (t32+t55)*t58*t67-(t17+t21-t2+t69+t70-t71+2*t72)*t75*t58
     #*t67/512-(t80-t2+t69+t70-t71+t21+t17)*t82*t58*t67/1024+t82*t75*t58
     #*t67/512-(t21+t91-t2+t69+t70-t71+t17)*t93*t58*t67/1024+t75*t93*t58
     #*t67/512+t93*t82*t58*t67/1024
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = f2(l3,z)
      t5 = dlog(t4)
      t7 = dlog(1-z)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      t14 = dlog(-l3+1)
      t16 = dlog(l3)
      t19 = 1/l3
      t21 = s34r2(z,0.D0,0.D0,l3,l4)
      t22 = 1/t21
      t24 = s23r2(z,0.D0,0.D0,l3,l4)
      t26 = s1vr2(z,0.D0,0.D0,l3,l4)
      t27 = t26**2
      t30 = f1(l3,z,0.D0)
      t31 = 1/t24/t27*t30
      tncoff2 = (t1+t3-t5+4*t7+4*t9-2*t11+2*t14+2*t16)*t19*t22*t31/512-1
     #/l2*t19*t22*t31/512-t19/l1*t22*t31/512
      t1 = dlog(l2)
      t3 = dlog(1-l2)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(1-z)
      t14 = dlog(1-l4)
      t15 = dlog(l4)
      t17 = 1/l2
      t19 = s34r2(z,0.D0,l2,0.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,l2,0.D0,l4)
      t24 = s1vr2(z,0.D0,l2,0.D0,l4)
      t25 = t24**2
      t28 = f1(0.D0,z,0.D0)
      t29 = 1/t22/t25*t28
      tncoff3 = (t1+t3+4*t4-2*t6-t9+4*t11+t14+t15)*t17*t20*t29/1024-t17/
     #l3*t20*t29/512-1/l1*t17*t20*t29/1024
      t4 = s34r2(z,0.D0,l2,l3,l4)
      t7 = s23r2(z,0.D0,l2,l3,l4)
      t9 = s1vr2(z,0.D0,l2,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,0.D0)
      tncoff4 = 1/l2/l3/t4/t7/t10*t13/512
      t2 = dlog(1-l4)
      t3 = dlog(l1)
      t5 = dlog(1-z)
      t7 = f3(0.D0,z)
      t10 = dlog(1-l1*t7)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(1-l1)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = dlog(l4)
      t21 = 1/l1
      t23 = s34r2(z,l1,0.D0,0.D0,l4)
      t24 = 1/t23
      t26 = s23r2(z,l1,0.D0,0.D0,l4)
      t28 = s1vr2(z,l1,0.D0,0.D0,l4)
      t29 = t28**2
      t32 = f1(0.D0,z,l1)
      t33 = 1/t26/t29*t32
      tncoff5 = (t2+t3+4*t5+t10-t12+t14+4*t15-2*t17+t19)*t21*t24*t33/102
     #4-1/l3*t21*t24*t33/512-t21/l2*t24*t33/1024
      t4 = s34r2(z,l1,0.D0,l3,l4)
      t7 = s23r2(z,l1,0.D0,l3,l4)
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff6 = 1/l3/l1/t4/t7/t10*t13/512
      t4 = s34r2(z,l1,l2,0.D0,l4)
      t7 = s23r2(z,l1,l2,0.D0,l4)
      t9 = s1vr2(z,l1,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,z,l1)
      tncoff7 = 1/l1/l2/t4/t7/t10*t13/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = t22**2
      t26 = f1(0.D0,z,0.D0)
      t27 = 1/t20/t23*t26
      tncoff8 = (-t2/1024-t4/1024+t7/1024-t10/256-t12/256+t14/512)*t18*t
     #27+1/l3*t18*t27/512+1/l2*t18*t27/1024+1/l1*t18*t27/1024
      t2 = s34r2(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,0.D0)
      tncoff9 = -1/l3/t2/t5/t8*t11/512
      t2 = s34r2(z,0.D0,l2,0.D0,l4)
      t5 = s23r2(z,0.D0,l2,0.D0,l4)
      t7 = s1vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,0.D0)
      tncoff10 = -1/l2/t2/t5/t8*t11/1024
      t2 = s34r2(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff11 = -1/l1/t2/t5/t8*t11/1024
      t1 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,0.D0)
      tncoff12 = 1/t1/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num2e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num2e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num2e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num2e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num2e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num2e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num2e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      rm13e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm13e2int(var,ndim,Tcoeff)

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
      double precision rm1num2e0
      external rm1num2e0
      double precision rm1num2e1
      external rm1num2e1
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
      t4 = dlog(l4)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t12 = dlog(2.D0)
      t14 = dlog(z)
      t17 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = t22**2
      t26 = f1(0.D0,z,0.D0)
      t27 = 1/t20/t23*t26
      tncoff1 = (-t2/1024-t4/1024+t7/1024-t10/256-t12/256+t14/512)*t18*t
     #27+1/l3*t18*t27/512+1/l2*t18*t27/1024+1/l1*t18*t27/1024
      t2 = s34r2(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,0.D0)
      tncoff2 = -1/l3/t2/t5/t8*t11/512
      t2 = s34r2(z,0.D0,l2,0.D0,l4)
      t5 = s23r2(z,0.D0,l2,0.D0,l4)
      t7 = s1vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,0.D0)
      tncoff3 = -1/l2/t2/t5/t8*t11/1024
      t2 = s34r2(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,z,l1)
      tncoff4 = -1/l1/t2/t5/t8*t11/1024
      t1 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,0.D0)
      tncoff5 = 1/t1/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num2e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num2e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      rm13e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm13e3int(var,ndim,Tcoeff)

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
      double precision rm1num2e0
      external rm1num2e0
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

      t1 = s34r2(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,0.D0)
      tncoff1 = 1/t1/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      rm13e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm13e4int(var,ndim,Tcoeff)

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

      rm13e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

