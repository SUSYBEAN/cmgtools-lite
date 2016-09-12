CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function rm12e0int(var,ndim,Tcoeff)

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
      double precision rm1num1e0
      external rm1num1e0
      double precision rm1num1e1
      external rm1num1e1
      double precision rm1num1e2
      external rm1num1e2
      double precision rm1num1e3
      external rm1num1e3
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
      t3 = dlog(l4)
      t4 = t2*t3
      t6 = dlog(1-l2)
      t9 = dlog(l2)
      t12 = dlog(2.D0)
      t15 = t2*t6
      t18 = t2*t12
      t21 = t2*t9
      t24 = t3*t6
      t27 = t3*t12
      t30 = t3*t9
      t33 = t6*t9
      t36 = f2(1.D0,1.D0)
      t37 = dlog(t36)
      t38 = t2*t37
      t45 = t37*t6
      t50 = t4*t6/4096+t4*t9/4096+t4*t12/2048+t15*t9/4096+t6*t18/2048+t2
     #1*t12/2048+t24*t9/4096+t27*t6/2048+t30*t12/2048+t33*t12/2048-t38*t
     #6/4096-t38*t9/4096-t38*t12/2048-t45*t9/4096-t45*t12/2048
      t51 = t37*t9
      t54 = t3*t37
      t63 = 0.3141592653589793D1**2
      t66 = t9**2
      t69 = t3**2
      t74 = t12**2
      t79 = t6**2
      t82 = t2**2
      t91 = -t51*t12/2048-t54*t6/4096-t54*t9/4096-t54*t12/2048-t4*t37/40
     #96-t63*t12/6144+t12*t66/4096+t2*t69/8192+t2*t66/8192+t2*t74/2048-t
     #2*t63/12288+t2*t79/8192+t82*t3/8192+t82*t6/8192+t82*t9/8192+t82*t1
     #2/4096
      t123 = t37**2
      t126 = t3*t66/8192+t3*t74/2048-t3*t63/12288+t3*t79/8192+t9*t74/204
     #8-t9*t63/12288+t6*t66/8192+t6*t74/2048-t6*t63/12288+t79*t9/8192+t7
     #9*t12/4096+t69*t6/8192+t69*t9/8192+t69*t12/4096-t69*t37/8192+t3*t1
     #23/8192
      t158 = -t82*t37/8192+t2*t123/8192+t123*t6/8192+t123*t9/8192+t123*t
     #12/4096-t37*t79/8192-t37*t66/8192-t37*t74/2048+t37*t63/12288-t123*
     #t37/24576+z3/1536+t74*t12/3072+t82*t2/24576+t69*t3/24576+t79*t6/24
     #576+t66*t9/24576
      t161 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t162 = 1/t161
      t164 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t165 = t164**2
      t166 = 1/t165
      t167 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t168 = t167**2
      t169 = 1/t168
      t171 = f1(1.D0,1.D0,0.D0)
      t172 = t166*t169*t171
      t174 = -1+z
      t175 = 1/t174
      t176 = 1/l3
      t178 = 1/l1
      t182 = t162*t166*t169*t171
      t185 = 6*t4
      t187 = 12*t9*t12
      t188 = 6*t33
      t190 = 12*t6*t12
      t191 = 6*t24
      t192 = 6*t30
      t193 = 12*t27
      t194 = 6*t15
      t195 = 6*t21
      t196 = 12*t18
      t197 = dlog(l3)
      t206 = t185+t187+t188+t190+t191+t192+t193+t194+t195+t196+12*t197*t
     #6+12*t197*t9+24*t197*t12+12*t2*t197
      t209 = 6*t54
      t210 = 6*t38
      t213 = 6*t45
      t214 = 6*t51
      t216 = 12*t37*t12
      t217 = 3*t123
      t218 = t197**2
      t220 = 3*t79
      t221 = 3*t82
      t222 = 2*t63
      t223 = 3*t66
      t224 = 3*t69
      t225 = 12*t74
      t226 = 12*t3*t197-t209-t210-12*t197*t37-t213-t214-t216+t217+12*t21
     #8+t220+t221-t222+t223+t224+t225
      t232 = t185+t187+t188+t190+t191+t192+t193+t194+t195+t196-t209-t210
     #-t213-t214
      t233 = dlog(l1)
      t246 = t233**2
      t248 = -t216+6*t233*t6+6*t233*t9+12*t233*t12+6*t3*t233+6*t2*t233-6
     #*t233*t37+t217+3*t246+t220+t221-t222+t223+t224+t225
      t255 = dlog(-t174)
      t268 = t255**2
      t270 = -t216+48*t255*t12+24*t2*t255+24*t3*t255+24*t255*t6+24*t255*
     #t9-24*t255*t37+t217+48*t268+t220+t221-t222+t223+t224+t225
      t276 = 2*t197
      t277 = 2*t12
      t278 = 4*t255
      tncoff1 = (t50+t91+t126+t158)*t162*t172+t175*t176*t178*t182/512-(t
     #206+t226)*t176*t162*t172/12288-(t232+t248)*t178*t162*t172/24576+(t
     #232+t270)*t175*t162*t172/6144-(t276+t6-t37+t9+t277+t278+t2+t3)*t17
     #5*t176*t182/512-(t3+t233+t6-t37+t9+t277+t278+t2)*t178*t175*t182/10
     #24+(t276+t6+t277+t233+t9-t37+t2+t3)*t176*t178*t182/2048
      t1 = dlog(l3)
      t2 = 2*t1
      t4 = dlog(1-l2)
      t6 = 1-l3/2
      t7 = f2(t6,1.D0)
      t8 = dlog(t7)
      t9 = dlog(l2)
      t10 = -l3+2
      t11 = dlog(t10)
      t12 = 2*t11
      t13 = 1-z
      t14 = dlog(t13)
      t17 = dlog(1-l4)
      t18 = dlog(l4)
      t21 = -1/t10
      t23 = 1/l3
      t25 = -1/t13
      t28 = s34r2(1.D0,0.D0,l2,t6,l4)
      t29 = 1/t28
      t30 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t31 = t30**2
      t32 = 1/t31
      t34 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t35 = t34**2
      t36 = 1/t35
      t37 = f1(t6,1.D0,0.D0)
      t39 = t29*t32*t36*t37
      t70 = -6*t17*t8-6*t8*t4-6*t8*t9-12*t8*t11-12*t1*t8+12*t9*t11+12*t4
     #*t11+12*t17*t11+12*t18*t11+24*t1*t11+6*t17*t18+6*t4*t9+6*t18*t4+6*
     #t18*t9
      t85 = 0.3141592653589793D1**2
      t87 = t4**2
      t89 = t9**2
      t91 = t11**2
      t93 = t18**2
      t95 = t17**2
      t97 = t1**2
      t99 = t8**2
      t101 = 6*t17*t4+6*t17*t9+12*t1*t4+12*t1*t9+12*t17*t1+12*t18*t1-6*t
     #18*t8-2*t85+3*t87+3*t89+12*t91+3*t93+3*t95+12*t97+3*t99
      t107 = dlog(l1)
      t111 = 1/l1*t23
      tncoff2 = (-t2-t4+t8-t9-t12-4*t14-t17-t18)*t21*t23*t25*t39/256-(t7
     #0+t101)*t21*t23*t39/6144-(-t2+t8-t107-t4-t9-t12-t17-t18)*t21*t111*
     #t39/1024+t21*t29*t36*t32*t111*t25*t37/256
      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t3 = t2**2
      t6 = dlog(1-l2)
      t7 = t6**2
      t9 = dlog(l2)
      t10 = t9**2
      t12 = dlog(l4)
      t13 = t12**2
      t16 = dlog(1-l4)
      t17 = t16**2
      t19 = dlog(l1)
      t20 = t19**2
      t23 = dlog(1-l1)
      t24 = t23**2
      t26 = 0.3141592653589793D1**2
      t28 = dlog(2.D0)
      t29 = t28**2
      t31 = f3(1.D0,1.D0)
      t34 = dlog(1-l1*t31)
      t35 = t34**2
      t39 = 3*t3+3*t7+3*t10+3*t13+3*t17+3*t20+3*t24-2*t26+12*t29+3*t35+6
     #*t19*t23
      t64 = 6*t12*t23+6*t16*t23+6*t9*t23+12*t23*t28+6*t6*t23+6*t16*t12+1
     #2*t9*t28+6*t6*t9+12*t28*t6+6*t6*t12+6*t9*t12+12*t12*t28
      t88 = 6*t16*t6+6*t16*t9+12*t16*t28+6*t19*t6+6*t19*t9+12*t19*t28+6*
     #t12*t19+6*t16*t19-6*t2*t12-6*t19*t2-6*t16*t2
      t113 = -6*t2*t6-6*t2*t9-12*t2*t28+6*t34*t6-6*t34*t2+6*t34*t9+12*t3
     #4*t28+6*t34*t23-6*t2*t23+6*t19*t34+6*t12*t34+6*t16*t34
      t116 = 1/l1
      t118 = s34r2(1.D0,l1,l2,1.D0,l4)
      t119 = 1/t118
      t121 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t122 = t121**2
      t123 = 1/t122
      t124 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t125 = t124**2
      t126 = 1/t125
      t128 = f1(1.D0,1.D0,l1)
      t132 = 2*t28
      t133 = 1-z
      t134 = dlog(t133)
      t139 = -1/t133
      t143 = t119*t123*t126*t128
      t146 = dlog(l3)
      t149 = 1/l3
      tncoff3 = (t39+t64+t88+t113)*t116*t119*t123*t126*t128/24576+(t12+t
     #19-t2+t34+t6+t9+t132+t23+4*t134+t16)*t116*t139*t143/1024-(2*t146-t
     #2+t132+t19+t34+t6+t9+t23+t16+t12)*t149*t116*t143/2048-t139*t149*t1
     #16*t143/512
      t1 = dlog(l3)
      t3 = dlog(l1)
      t5 = 1-l3/2
      t6 = f2(t5,1.D0)
      t7 = dlog(t6)
      t8 = f3(t5,1.D0)
      t11 = dlog(1-l1*t8)
      t13 = dlog(1-l2)
      t15 = dlog(1-l1)
      t16 = dlog(l2)
      t17 = -l3+2
      t18 = dlog(t17)
      t21 = dlog(1-l4)
      t22 = dlog(l4)
      t25 = -1/t17
      t29 = 1/l1/l3
      t31 = s34r2(1.D0,l1,l2,t5,l4)
      t32 = 1/t31
      t33 = s2vr2(1.D0,l1,l2,t5,l4)
      t34 = t33**2
      t35 = 1/t34
      t37 = s1vr2(1.D0,l1,l2,t5,l4)
      t38 = t37**2
      t39 = 1/t38
      t40 = f1(t5,1.D0,l1)
      tncoff4 = -(2*t1+t3-t7+t11+t13+t15+t16+2*t18+t21+t22)*t25*t29*t32*
     #t35*t39*t40/1024-t32*t25*t39*t35*t29/(-1+z)*t40/256
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l3
      t5 = 1/l1
      t7 = s34r2(z,0.D0,l2,1.D0,l4)
      t8 = 1/t7
      t9 = s2vr2(z,0.D0,l2,1.D0,l4)
      t10 = t9**2
      t11 = 1/t10
      t13 = s1vr2(z,0.D0,l2,1.D0,l4)
      t14 = t13**2
      t15 = 1/t14
      t16 = f1(1.D0,z,0.D0)
      t18 = t8*t11*t15*t16
      t21 = dlog(2.D0)
      t22 = t21**2
      t25 = dlog(-t1)
      t26 = t25**2
      t28 = dlog(z)
      t29 = t28**2
      t31 = f2(1.D0,z)
      t32 = dlog(t31)
      t33 = t32**2
      t36 = dlog(1-l2)
      t37 = t36**2
      t39 = dlog(l2)
      t40 = t39**2
      t42 = dlog(l4)
      t43 = t42**2
      t46 = dlog(1-l4)
      t47 = t46**2
      t49 = 0.3141592653589793D1**2
      t69 = 12*t22+48*t26+12*t29+3*t33+3*t37+3*t40+3*t43+3*t47-2*t49-24*
     #t25*t32-6*t46*t32-6*t42*t32+12*t32*t28-6*t32*t36-6*t32*t39-12*t32*
     #t21+6*t46*t42+12*t39*t21
      t108 = 6*t36*t39+12*t21*t36+6*t42*t36+6*t42*t39+12*t42*t21+6*t46*t
     #36+6*t46*t39+12*t46*t21+48*t25*t21+24*t46*t25+24*t42*t25+24*t25*t3
     #6+24*t25*t39-12*t46*t28-12*t42*t28-48*t25*t28-12*t39*t28-24*t28*t2
     #1-12*t36*t28
      t116 = 2*t28
      t117 = 2*t21
      t118 = dlog(l3)
      t120 = 4*t25
      t126 = dlog(l1)
      tncoff5 = -t2*t3*t5*t18/512-(t69+t108)*t2*t8*t11*t15*t16/6144+(-t1
     #16+t117-t32+2*t118+t36+t39+t120+t46+t42)*t2*t3*t18/512+(t42+t126+t
     #39+t117-t116+t36-t32+t120+t46)*t5*t2*t18/1024
      t1 = l3-2
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s34r2(z,0.D0,l2,t4,l4)
      t6 = 1/t5
      t8 = s1vr2(z,0.D0,l2,t4,l4)
      t9 = t8**2
      t10 = 1/t9
      t11 = s2vr2(z,0.D0,l2,t4,l4)
      t12 = t11**2
      t13 = 1/t12
      t17 = 1/l3
      t19 = -1+z
      t20 = 1/t19
      t21 = f1(t4,z,0.D0)
      t26 = dlog(-t1)
      t29 = dlog(1-l2)
      t30 = dlog(l3)
      t32 = dlog(z)
      t34 = f2(t4,z)
      t35 = dlog(t34)
      t36 = dlog(l2)
      t38 = dlog(-t19)
      t41 = dlog(1-l4)
      t42 = dlog(l4)
      tncoff6 = -t2*t6*t10*t13/l1*t17*t20*t21/256-(-2*t26-t29-2*t30+2*t3
     #2+t35-t36-4*t38-t41-t42)*t2*t17*t20*t6*t13*t10*t21/256
      t1 = dlog(l4)
      t2 = dlog(l2)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = dlog(2.D0)
      t7 = dlog(z)
      t10 = dlog(1-l1)
      t11 = dlog(l1)
      t12 = f3(1.D0,z)
      t15 = dlog(1-l1*t12)
      t17 = dlog(1-l2)
      t18 = 1-z
      t19 = dlog(t18)
      t22 = dlog(1-l4)
      t23 = t1+t2-t4+2*t5-2*t7+t10+t11+t15+t17+4*t19+t22
      t24 = 1/l1
      t27 = -1/t18
      t29 = s34r2(z,l1,l2,1.D0,l4)
      t31 = s2vr2(z,l1,l2,1.D0,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,l2,1.D0,l4)
      t36 = t35**2
      t38 = f1(1.D0,z,l1)
      t40 = 1/t29/t32/t36*t38
      tncoff7 = -t23*t24*t27*t40/1024+t27/l3*t24*t40/512
      t2 = 1-l3/2
      t3 = f1(t2,z,l1)
      t4 = s34r2(z,l1,l2,t2,l4)
      t9 = s1vr2(z,l1,l2,t2,l4)
      t10 = t9**2
      t14 = s2vr2(z,l1,l2,t2,l4)
      t15 = t14**2
      tncoff8 = t3/t4/(l3-2)/t10/t15/l1/l3/(-1+z)/256
      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l2)
      t6 = t5**2
      t8 = dlog(2.D0)
      t9 = t8**2
      t11 = 0.3141592653589793D1**2
      t14 = dlog(1-l2)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t32 = -t3/8192-t6/8192-t9/2048+t11/12288-t15/8192-t18/8192-t22/819
     #2+t2*t14/4096+t2*t5/4096+t2*t8/2048-t14*t5/4096
      t55 = -t8*t14/2048-t5*t8/2048+t17*t2/4096-t17*t14/4096-t17*t5/4096
     #-t17*t8/2048-t21*t17/4096+t2*t21/4096-t21*t14/4096-t21*t5/4096-t21
     #*t8/2048
      t57 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t61 = t60**2
      t63 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t64 = t63**2
      t67 = f1(1.D0,1.D0,0.D0)
      t68 = 1/t61/t64*t67
      t70 = -1+z
      t71 = 1/t70
      t72 = 1/l3
      t77 = 1/l1
      t86 = dlog(l3)
      t88 = 6*t2
      t89 = 6*t14
      t90 = 6*t5
      t91 = 12*t8
      t92 = 6*t17
      t93 = 6*t21
      t99 = dlog(l1)
      t107 = dlog(-t70)
      tncoff9 = (t32+t55)*t58*t68+t71*t72*t58*t68/512+t77*t71*t58*t68/10
     #24-t72*t77*t58*t68/2048-(-12*t86+t88-t89-t90-t91-t92-t93)*t72*t58*
     #t68/12288-(t88-t89-t90-t91-t92-t93-6*t99)*t77*t58*t68/24576+(-24*t
     #107+t88-t89-t90-t91-t92-t93)*t71*t58*t68/6144
      t1 = dlog(l3)
      t4 = dlog(1-l4)
      t6 = dlog(l4)
      t9 = dlog(1-l2)
      t12 = 1-l3/2
      t13 = f2(t12,1.D0)
      t14 = dlog(t13)
      t16 = dlog(l2)
      t18 = -l3+2
      t19 = dlog(t18)
      t23 = -1/t18
      t25 = 1/l3
      t27 = s34r2(1.D0,0.D0,l2,t12,l4)
      t29 = s2vr2(1.D0,0.D0,l2,t12,l4)
      t30 = t29**2
      t33 = s1vr2(1.D0,0.D0,l2,t12,l4)
      t34 = t33**2
      t36 = f1(t12,1.D0,0.D0)
      t38 = 1/t27/t30/t34*t36
      tncoff10 = -(-12*t1-6*t4-6*t6-6*t9+6*t14-6*t16-12*t19)*t23*t25*t38
     #/6144-t23/l1*t25*t38/1024+t23*t25/(-1+z)*t38/256
      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t4 = dlog(l1)
      t6 = dlog(l4)
      t9 = dlog(1-l4)
      t11 = dlog(l2)
      t13 = dlog(2.D0)
      t16 = dlog(1-l1)
      t18 = f3(1.D0,1.D0)
      t21 = dlog(1-l1*t18)
      t24 = dlog(1-l2)
      t27 = 1/l1
      t29 = s34r2(1.D0,l1,l2,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t33 = t32**2
      t35 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t36 = t35**2
      t39 = f1(1.D0,1.D0,l1)
      t40 = 1/t33/t36*t39
      tncoff11 = (6*t2-6*t4-6*t6-6*t9-6*t11-12*t13-6*t16-6*t21-6*t24)*t2
     #7*t30*t40/24576-t27/(-1+z)*t30*t40/1024+1/l3*t27*t30*t40/2048
      t8 = 1-l3/2
      t9 = s34r2(1.D0,l1,l2,t8,l4)
      t11 = s2vr2(1.D0,l1,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(1.D0,l1,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,1.D0,l1)
      tncoff12 = 1/(l3-2)/l1/l3/t9/t12/t16*t18/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = 1-z
      t7 = dlog(t6)
      t9 = dlog(z)
      t12 = dlog(1-l2)
      t14 = dlog(l2)
      t16 = dlog(2.D0)
      t18 = f2(1.D0,z)
      t19 = dlog(t18)
      t23 = -1/t6
      t25 = s34r2(z,0.D0,l2,1.D0,l4)
      t26 = 1/t25
      t28 = s2vr2(z,0.D0,l2,1.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,0.D0,l2,1.D0,l4)
      t32 = t31**2
      t35 = f1(1.D0,z,0.D0)
      t36 = 1/t29/t32*t35
      tncoff13 = -(-6*t2-6*t4-24*t7+12*t9-6*t12-6*t14-12*t16+6*t19)*t23*
     #t26*t36/6144-t23/l3*t26*t36/512-1/l1*t23*t26*t36/1024
      t9 = 1-l3/2
      t10 = s34r2(z,0.D0,l2,t9,l4)
      t12 = s2vr2(z,0.D0,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(z,0.D0,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,z,0.D0)
      tncoff14 = -1/(l3-2)/l3/(-1+z)/t10/t13/t17*t19/256
      t5 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,l2,1.D0,l4)
      t12 = t11**2
      t15 = f1(1.D0,z,l1)
      tncoff15 = 1/l1/(-1+z)/t5/t9/t12*t15/1024
      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t5 = dlog(1-l2)
      t7 = dlog(l2)
      t9 = dlog(2.D0)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t21 = t20**2
      t23 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t24 = t23**2
      t27 = f1(1.D0,1.D0,0.D0)
      t28 = 1/t21/t24*t27
      tncoff16 = (-t2/4096+t5/4096+t7/4096+t9/2048+t11/4096+t14/4096)*t1
     #8*t28-1/l3*t18*t28/2048-1/l1*t18*t28/4096+1/(-1+z)*t18*t28/1024
      t6 = 1-l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t11 = t10**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t17 = f1(t6,1.D0,0.D0)
      tncoff17 = -1/(l3-2)/l3/t7/t11/t14*t17/1024
      t2 = s34r2(1.D0,l1,l2,1.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,1.D0,l1)
      tncoff18 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,z,0.D0)
      tncoff19 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = f1(1.D0,1.D0,0.D0)
      tncoff20 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num1e0(var,ndim,z,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = rm1num1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = rm1num1e1(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = rm1num1e1(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = rm1num1e2(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = rm1num1e2(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = rm1num1e3(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      rm12e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff1
     #6+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncof
     #f20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm12e1int(var,ndim,Tcoeff)

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
      double precision rm1num1e0
      external rm1num1e0
      double precision rm1num1e1
      external rm1num1e1
      double precision rm1num1e2
      external rm1num1e2
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

      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l2)
      t6 = t5**2
      t8 = dlog(2.D0)
      t9 = t8**2
      t11 = 0.3141592653589793D1**2
      t14 = dlog(1-l2)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t32 = -t3/8192-t6/8192-t9/2048+t11/12288-t15/8192-t18/8192-t22/819
     #2+t2*t14/4096+t2*t5/4096+t2*t8/2048-t14*t5/4096
      t55 = -t8*t14/2048-t5*t8/2048+t17*t2/4096-t17*t14/4096-t17*t5/4096
     #-t17*t8/2048-t21*t17/4096+t2*t21/4096-t21*t14/4096-t21*t5/4096-t21
     #*t8/2048
      t57 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t61 = t60**2
      t63 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t64 = t63**2
      t67 = f1(1.D0,1.D0,0.D0)
      t68 = 1/t61/t64*t67
      t70 = 1-z
      t71 = dlog(t70)
      t73 = 2*t8
      t76 = -1/t70
      t81 = dlog(l3)
      t84 = 1/l3
      t93 = dlog(l1)
      t95 = 1/l1
      tncoff1 = (t32+t55)*t58*t68-(t21+t17+4*t71+t14-t2+t5+t73)*t76*t58*
     #t68/1024+(t21+t14-t2+t5+t73+2*t81+t17)*t84*t58*t68/2048+t76*t84*t5
     #8*t68/512+(t17+t93+t14-t2+t5+t73+t21)*t95*t58*t68/4096+t95*t76*t58
     #*t68/1024-t84*t95*t58*t68/2048
      t2 = dlog(1-l4)
      t4 = dlog(1-l2)
      t6 = 1-l3/2
      t7 = f2(t6,1.D0)
      t8 = dlog(t7)
      t9 = dlog(l2)
      t10 = -l3+2
      t11 = dlog(t10)
      t13 = dlog(l3)
      t15 = dlog(l4)
      t18 = -1/t10
      t20 = 1/l3
      t22 = s34r2(1.D0,0.D0,l2,t6,l4)
      t24 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t25 = t24**2
      t28 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t29 = t28**2
      t31 = f1(t6,1.D0,0.D0)
      t33 = 1/t22/t25/t29*t31
      tncoff2 = -(-t2-t4+t8-t9-2*t11-2*t13-t15)*t18*t20*t33/1024+t18*t20
     #/(-1+z)*t33/256-t18/l1*t20*t33/1024
      t1 = dlog(l4)
      t2 = dlog(l1)
      t3 = f2(1.D0,1.D0)
      t4 = dlog(t3)
      t5 = f3(1.D0,1.D0)
      t8 = dlog(1-l1*t5)
      t10 = dlog(1-l2)
      t11 = dlog(l2)
      t12 = dlog(2.D0)
      t15 = dlog(1-l1)
      t17 = dlog(1-l4)
      t19 = 1/l1
      t21 = s34r2(1.D0,l1,l2,1.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t25 = t24**2
      t27 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t28 = t27**2
      t31 = f1(1.D0,1.D0,l1)
      t32 = 1/t25/t28*t31
      tncoff3 = -(t1+t2-t4+t8+t10+t11+2*t12+t15+t17)*t19*t22*t32/4096-t1
     #9/(-1+z)*t22*t32/1024+1/l3*t19*t22*t32/2048
      t8 = 1-l3/2
      t9 = s34r2(1.D0,l1,l2,t8,l4)
      t11 = s2vr2(1.D0,l1,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(1.D0,l1,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,1.D0,l1)
      tncoff4 = 1/(l3-2)/l1/l3/t9/t12/t16*t18/1024
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = 1-z
      t5 = dlog(t4)
      t8 = dlog(1-l2)
      t9 = f2(1.D0,z)
      t10 = dlog(t9)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t15 = dlog(l2)
      t18 = -1/t4
      t20 = s34r2(z,0.D0,l2,1.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,0.D0,l2,1.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,0.D0,l2,1.D0,l4)
      t27 = t26**2
      t30 = f1(1.D0,z,0.D0)
      t31 = 1/t24/t27*t30
      tncoff5 = (t2+t3+4*t5+t8-t10+2*t11-2*t13+t15)*t18*t21*t31/1024-t18
     #/l3*t21*t31/512-1/l1*t18*t21*t31/1024
      t9 = 1-l3/2
      t10 = s34r2(z,0.D0,l2,t9,l4)
      t12 = s2vr2(z,0.D0,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(z,0.D0,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,z,0.D0)
      tncoff6 = -1/(l3-2)/l3/(-1+z)/t10/t13/t17*t19/256
      t5 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,l2,1.D0,l4)
      t12 = t11**2
      t15 = f1(1.D0,z,l1)
      tncoff7 = 1/l1/(-1+z)/t5/t9/t12*t15/1024
      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t5 = dlog(1-l2)
      t7 = dlog(l2)
      t9 = dlog(2.D0)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t21 = t20**2
      t23 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t24 = t23**2
      t27 = f1(1.D0,1.D0,0.D0)
      t28 = 1/t21/t24*t27
      tncoff8 = (-t2/4096+t5/4096+t7/4096+t9/2048+t11/4096+t14/4096)*t18
     #*t28-1/l3*t18*t28/2048-1/l1*t18*t28/4096+1/(-1+z)*t18*t28/1024
      t6 = 1-l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t11 = t10**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t17 = f1(t6,1.D0,0.D0)
      tncoff9 = -1/(l3-2)/l3/t7/t11/t14*t17/1024
      t2 = s34r2(1.D0,l1,l2,1.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,1.D0,l1)
      tncoff10 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,z,0.D0)
      tncoff11 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = f1(1.D0,1.D0,0.D0)
      tncoff12 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      rm12e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm12e2int(var,ndim,Tcoeff)

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
      double precision rm1num1e0
      external rm1num1e0
      double precision rm1num1e1
      external rm1num1e1
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

      t1 = f2(1.D0,1.D0)
      t2 = dlog(t1)
      t5 = dlog(1-l2)
      t7 = dlog(l2)
      t9 = dlog(2.D0)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t21 = t20**2
      t23 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t24 = t23**2
      t27 = f1(1.D0,1.D0,0.D0)
      t28 = 1/t21/t24*t27
      tncoff1 = (-t2/4096+t5/4096+t7/4096+t9/2048+t11/4096+t14/4096)*t18
     #*t28-1/l3*t18*t28/2048-1/l1*t18*t28/4096+1/(-1+z)*t18*t28/1024
      t6 = 1-l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t11 = t10**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t17 = f1(t6,1.D0,0.D0)
      tncoff2 = -1/(l3-2)/l3/t7/t11/t14*t17/1024
      t2 = s34r2(1.D0,l1,l2,1.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,1.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,1.D0,l1)
      tncoff3 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,1.D0,l4)
      t10 = t9**2
      t13 = f1(1.D0,z,0.D0)
      tncoff4 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = f1(1.D0,1.D0,0.D0)
      tncoff5 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      rm12e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm12e3int(var,ndim,Tcoeff)

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
      double precision rm1num1e0
      external rm1num1e0
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

      t1 = s34r2(1.D0,0.D0,l2,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,1.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,1.D0,l4)
      t8 = t7**2
      t10 = f1(1.D0,1.D0,0.D0)
      tncoff1 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      rm12e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm12e4int(var,ndim,Tcoeff)

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

      rm12e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

