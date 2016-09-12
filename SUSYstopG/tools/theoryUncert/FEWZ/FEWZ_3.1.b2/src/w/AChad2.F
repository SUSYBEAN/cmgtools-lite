CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function AC2e0int(var,ndim,Tcoeff)

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
      double precision ACnum1e0
      external ACnum1e0
      double precision ACnum1e1
      external ACnum1e1
      double precision ACnum1e2
      external ACnum1e2
      double precision ACnum1e3
      external ACnum1e3
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

      t1 = dlog(2.D0)
      t2 = t1**2
      t6 = dlog(1-l4)
      t7 = t6**2
      t11 = dlog(l4)
      t12 = t11**2
      t15 = dlog(z)
      t16 = t15**2
      t20 = dlog(1-z)
      t21 = t20**2
      t24 = 0.3141592653589793D1**2
      t27 = f2(1.D0,z)
      t28 = dlog(t27)
      t29 = t28**2
      t46 = t2*t1/384+t7*t6/3072+z3/192+t12*t11/3072-t16*t15/384+t21*t20
     #/48-t11*t24/1536-t29*t15/512+t1*t29/512-t21*t15/32-t21*t28/64+t28*
     #t24/1536+t16*t1/128-t28*t16/256+t12*t20/256
      t51 = t11*t28
      t54 = t6*t11
      t57 = t6*t20
      t76 = t20*t15
      t79 = t11*t20
      t84 = t24*t15/768-t2*t28/256-t51*t1/256-t54*t15/256-t57*t28/128+t2
     #1*t1/32+t7*t1/512+t7*t20/256+t7*t11/1024-t12*t15/512-t12*t28/1024-
     #t20*t24/384+t11*t21/64-t76*t1/32-t79*t28/128+t79*t1/64
      t86 = t1*t28
      t91 = t11*t15
      t120 = t86*t15/128-t57*t15/64-t91*t1/128-t79*t15/64+t11*t16/256+t2
     #0*t16/64+t20*t29/256+t20*t2/64-t15*t2/128-t7*t15/512-t7*t28/1024+t
     #6*t2/256+t6*t21/64+t12*t1/512+t6*t12/1024-t6*t24/1536
      t135 = t20*t28
      t138 = t20*t1
      t141 = t6*t28
      t154 = t6*t15
      t157 = t6*t16/256+t54*t20/128+t6*t29/1024+t11*t29/1024+t11*t2/256-
     #t24*t1/768+t57*t1/64+t135*t15/64-t138*t28/64+t141*t15/256+t51*t15/
     #256-t29*t28/3072+t54*t1/256-t54*t28/512-t141*t1/256-t154*t1/128
      t160 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t161 = 1/t160
      t163 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t164 = 1/t163
      t166 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t167 = 1/t166
      t168 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t169 = t168**2
      t170 = 1/t169
      t172 = f1(1.D0,z,0.D0)
      t175 = 12*t154
      t176 = 6*t51
      t177 = 12*t86
      t178 = 1-l2
      t179 = dlog(t178)
      t190 = 12*t91
      t192 = 12*t28*t15
      t193 = 6*t141
      t196 = 6*t54
      t198 = 12*t11*t1
      t199 = -t175-t176-t177+24*t179*t20-12*t179*t15-6*t179*t28+6*t6*t17
     #9+6*t11*t179-t190+t192-t193+12*t179*t1+t196+t198
      t200 = 24*t79
      t201 = 24*t135
      t202 = 2*t24
      t203 = 3*t29
      t204 = 48*t21
      t205 = 12*t16
      t206 = 3*t12
      t207 = 12*t2
      t208 = 48*t76
      t209 = 3*t7
      t211 = 24*t15*t1
      t212 = t179**2
      t214 = 24*t57
      t216 = 12*t6*t1
      t217 = 48*t138
      t218 = t200-t201-t202+t203+t204+t205+t206+t207-t208+t209-t211+3*t2
     #12+t214+t216+t217
      t221 = -1/t178
      t226 = t164*t167*t170*t172
      t229 = -t175-t176-t177-t190+t192-t193+t196+t198+t200-t201-t202+t20
     #3+t204+t205
      t230 = dlog(l3)
      t231 = t230**2
      t245 = t206+12*t231+t207-t208+t209-12*t230*t28+12*t11*t230+24*t230
     #*t1-24*t230*t15-t211+48*t230*t20+12*t6*t230+t214+t216+t217
      t247 = 1/l3
      t252 = 2*t230
      t253 = 2*t1
      t254 = 2*t15
      t255 = 4*t20
      t258 = t221*t161
      t262 = dlog(l1)
      t275 = t262**2
      t277 = -t175+6*t11*t262+6*t6*t262-t176-t177+12*t262*t1+24*t262*t20
     #-12*t262*t15-6*t262*t28+3*t275-t190+t192-t193+t196
      t278 = t198+t200-t201-t202+t203+t204+t205+t206+t207-t208+t209-t211
     #+t214+t216+t217
      t280 = 1/l1
      tncoff1 = (t46+t84+t120+t157)*t161*t164*t167*t170*t172+(t199+t218)
     #*t221*t161*t226/3072-(t229+t245)*t247*t161*t226/1536-(t252+t6+t11+
     #t179+t253-t254+t255-t28)*t247*t258*t226/256-(t277+t278)*t280*t161*
     #t226/3072-(t179+t255+t262-t28-t254+t253+t11+t6)*t280*t258*t226/512
     #+(t262-t254+t253+t255-t28+t252+t6+t11)*t247*t280*t161*t226/256+t24
     #7*t280*t258*t226/256
      t2 = dlog(1-l4)
      t3 = t2**2
      t6 = dlog(1-z)
      t9 = dlog(z)
      t10 = t9**2
      t12 = dlog(l4)
      t13 = t12**2
      t15 = 0.3141592653589793D1**2
      t17 = t6**2
      t20 = 1-l3/2
      t21 = f2(t20,z)
      t22 = dlog(t21)
      t23 = t22**2
      t25 = dlog(l3)
      t26 = t25**2
      t28 = -l3+2
      t29 = dlog(t28)
      t30 = t29**2
      t42 = 3*t3+24*t6*t2+12*t10+3*t13-2*t15+48*t17+3*t23+12*t26+12*t30-
     #6*t12*t22+24*t25*t29-12*t25*t22-6*t2*t22-12*t12*t9
      t73 = 12*t22*t9+12*t2*t29-12*t22*t29+12*t12*t29-24*t9*t29-24*t6*t2
     #2+48*t6*t29+12*t12*t25+12*t2*t25+48*t25*t6-24*t25*t9+6*t2*t12+24*t
     #6*t12-48*t6*t9-12*t2*t9
      t76 = -1/t28
      t78 = 1/l3
      t79 = s34r2(z,0.D0,1.D0,t20,l4)
      t80 = 1/t79
      t83 = s2vr2(z,0.D0,1.D0,t20,l4)
      t84 = 1/t83
      t85 = s23r2(z,0.D0,1.D0,t20,l4)
      t86 = 1/t85
      t88 = s1vr2(z,0.D0,1.D0,t20,l4)
      t89 = t88**2
      t90 = 1/t89
      t91 = f1(t20,z,0.D0)
      t96 = 2*t25
      t97 = 4*t6
      t98 = 2*t9
      t99 = 2*t29
      t100 = 1-l2
      t101 = dlog(t100)
      t105 = -1/t100
      t108 = t80*t84
      t111 = t108*t86*t90*t91
      t114 = dlog(l1)
      t117 = 1/l1
      tncoff2 = -(t42+t73)*t76*t78*t80*t84*t86*t90*t91/768-(t96+t2+t12+t
     #97-t98-t22+t99+t101)*t76*t105*t78*t111/128+(-t98+t114+t99+t97-t22+
     #t96+t2+t12)*t76*t117*t78*t111/128+t76*t117*t108*t86*t105*t78*t90*t
     #91/128
      t1 = f2(1.D0,z)
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l4)
      t6 = dlog(l2)
      t10 = dlog(1-l4)
      t13 = dlog(z)
      t18 = 1-l2
      t19 = dlog(t18)
      t24 = dlog(2.D0)
      t28 = dlog(1-z)
      t39 = t24**2
      t41 = t19**2
      t43 = t6**2
      t49 = 3*t3+6*t5*t6+6*t10*t6-12*t6*t13-6*t6*t2+6*t19*t6+6*t10*t5+12
     #*t5*t24+24*t5*t28-12*t10*t13+48*t28*t24-24*t13*t24-48*t28*t13+12*t
     #39+3*t41+3*t43-12*t5*t13+24*t6*t28
      t74 = t13**2
      t76 = t5**2
      t78 = 0.3141592653589793D1**2
      t80 = t28**2
      t82 = t10**2
      t88 = 12*t6*t24-24*t28*t2-6*t5*t2-12*t2*t24+24*t19*t28-12*t19*t13-
     #6*t19*t2+6*t10*t19+6*t5*t19+12*t2*t13-6*t10*t2+12*t19*t24+12*t74+3
     #*t76-2*t78+48*t80+3*t82+12*t10*t24+24*t10*t28
      t91 = -1/t18
      t93 = s34r2(z,0.D0,l2,1.D0,l4)
      t94 = 1/t93
      t96 = s2vr2(z,0.D0,l2,1.D0,l4)
      t98 = s23r2(z,0.D0,l2,1.D0,l4)
      t101 = s1vr2(z,0.D0,l2,1.D0,l4)
      t102 = t101**2
      t104 = f1(1.D0,z,0.D0)
      t106 = 1/t96/t98/t102*t104
      t109 = dlog(l3)
      t111 = 2*t13
      t112 = 4*t28
      t113 = 2*t24
      t115 = 1/l3
      t117 = t91*t94
      t121 = dlog(l1)
      t123 = 1/l1
      tncoff3 = -(t49+t88)*t91*t94*t106/3072+(2*t109+t10+t5-t111+t19+t11
     #2-t2+t113+t6)*t115*t117*t106/256+(t19+t6+t112-t2+t121+t113-t111+t5
     #+t10)*t123*t117*t106/512-t115*t123*t117*t106/256
      t1 = dlog(l3)
      t4 = dlog(1-l4)
      t5 = dlog(l4)
      t6 = -l3+2
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t13 = 1-l3/2
      t14 = f2(t13,z)
      t15 = dlog(t14)
      t16 = 1-l2
      t17 = dlog(t16)
      t18 = dlog(z)
      t20 = dlog(l2)
      t23 = -1/t6
      t28 = -1/t16/l3
      t30 = s34r2(z,0.D0,l2,t13,l4)
      t31 = 1/t30
      t32 = s2vr2(z,0.D0,l2,t13,l4)
      t33 = 1/t32
      t35 = s23r2(z,0.D0,l2,t13,l4)
      t36 = 1/t35
      t37 = s1vr2(z,0.D0,l2,t13,l4)
      t38 = t37**2
      t39 = 1/t38
      t41 = f1(t13,z,0.D0)
      tncoff4 = (2*t1+t4+t5+2*t7+4*t10-t15+t17-2*t18+t20)*t23*t28*t31*t3
     #3*t36*t39*t41/128-t23*t33*t31*t39*t36/l1*t28*t41/128
      t1 = dlog(z)
      t2 = t1**2
      t4 = dlog(l4)
      t5 = t4**2
      t7 = 0.3141592653589793D1**2
      t10 = dlog(1-z)
      t11 = t10**2
      t16 = dlog(1-l4)
      t17 = t16**2
      t19 = f2(1.D0,z)
      t20 = dlog(t19)
      t21 = dlog(2.D0)
      t24 = t21**2
      t26 = f3(1.D0,z)
      t29 = dlog(1-l1*t26)
      t30 = t29**2
      t36 = 12*t2+3*t5-2*t7+48*t11-12*t4*t1+3*t17-12*t20*t21+12*t24+3*t3
     #0+12*t16*t21+12*t20*t1
      t37 = dlog(l1)
      t39 = dlog(1-l1)
      t46 = t37**2
      t50 = t39**2
      t64 = 6*t37*t39+6*t16*t37+6*t4*t39+3*t46+12*t37*t21+3*t50-12*t39*t
     #1-24*t10*t20+24*t4*t10+48*t10*t21-24*t1*t21-48*t10*t1
      t78 = t20**2
      t88 = -6*t4*t20-12*t16*t1+24*t16*t10-12*t37*t1-6*t16*t20+6*t16*t29
     #+3*t78-6*t29*t20+6*t29*t39+12*t4*t21+24*t10*t39
      t113 = 6*t16*t4+6*t16*t39-6*t39*t20-12*t29*t1+6*t4*t29+6*t37*t29+6
     #*t4*t37+24*t37*t10-6*t37*t20+12*t39*t21+12*t29*t21+24*t10*t29
      t116 = 1/l1
      t118 = s34r2(z,l1,1.D0,1.D0,l4)
      t119 = 1/t118
      t121 = s2vr2(z,l1,1.D0,1.D0,l4)
      t123 = s23r2(z,l1,1.D0,1.D0,l4)
      t126 = s1vr2(z,l1,1.D0,1.D0,l4)
      t127 = t126**2
      t129 = f1(1.D0,z,l1)
      t131 = 1/t121/t123/t127*t129
      t134 = 1-l2
      t135 = dlog(t134)
      t136 = 4*t10
      t137 = 2*t21
      t138 = 2*t1
      t143 = -1/t134*t119
      t147 = dlog(l3)
      t150 = 1/l3
      tncoff5 = (t36+t64+t88+t113)*t116*t119*t131/3072+(t135+t37+t136+t1
     #37-t20+t29+t39-t138+t4+t16)*t116*t143*t131/512-(t137-t138-t20+t136
     #+t29+t37+t39+2*t147+t16+t4)*t150*t116*t119*t131/256-t150*t116*t143
     #*t131/256
      t1 = dlog(z)
      t4 = 1-l3/2
      t5 = f2(t4,z)
      t6 = dlog(t5)
      t7 = dlog(l1)
      t8 = -l3+2
      t9 = dlog(t8)
      t12 = dlog(1-z)
      t15 = dlog(1-l1)
      t16 = f3(t4,z)
      t19 = dlog(1-l1*t16)
      t20 = dlog(l3)
      t23 = dlog(1-l4)
      t24 = dlog(l4)
      t27 = -1/t8
      t29 = 1/l1
      t30 = 1/l3
      t33 = s34r2(z,l1,1.D0,t4,l4)
      t35 = s2vr2(z,l1,1.D0,t4,l4)
      t37 = 1/t33/t35
      t38 = s23r2(z,l1,1.D0,t4,l4)
      t39 = 1/t38
      t40 = s1vr2(z,l1,1.D0,t4,l4)
      t41 = t40**2
      t42 = 1/t41
      t44 = f1(t4,z,l1)
      tncoff6 = -(-2*t1-t6+t7+2*t9+4*t12+t15+t19+2*t20+t23+t24)*t27*t29*
     #t30*t37*t39*t42*t44/128-t27*t29*t37*t39/(l2-1)*t30*t42*t44/128
      t1 = 1-l2
      t2 = dlog(t1)
      t3 = dlog(l2)
      t4 = f2(1.D0,z)
      t5 = dlog(t4)
      t6 = dlog(l1)
      t8 = dlog(1-z)
      t10 = f3(1.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l1)
      t16 = dlog(2.D0)
      t18 = dlog(z)
      t20 = dlog(l4)
      t22 = dlog(1-l4)
      t23 = t2+t3-t5+t6+4*t8+t13+t15+2*t16-2*t18+t20+t22
      t24 = 1/l1
      t28 = s34r2(z,l1,l2,1.D0,l4)
      t30 = -1/t1/t28
      t32 = s2vr2(z,l1,l2,1.D0,l4)
      t34 = s23r2(z,l1,l2,1.D0,l4)
      t37 = s1vr2(z,l1,l2,1.D0,l4)
      t38 = t37**2
      t40 = f1(1.D0,z,l1)
      t42 = 1/t32/t34/t38*t40
      tncoff7 = -t23*t24*t30*t42/512+1/l3*t24*t30*t42/256
      t2 = 1-l3/2
      t3 = s2vr2(z,l1,l2,t2,l4)
      t5 = s34r2(z,l1,l2,t2,l4)
      t8 = s1vr2(z,l1,l2,t2,l4)
      t9 = t8**2
      t11 = s23r2(z,l1,l2,t2,l4)
      t23 = f1(t2,z,l1)
      tncoff8 = 1/t3/t5/t9/t11/(l3-2)/l1/(l2-1)/l3*t23/128
      t1 = dlog(l4)
      t2 = t1**2
      t5 = dlog(1-l4)
      t6 = t5**2
      t8 = dlog(2.D0)
      t12 = dlog(1-z)
      t17 = dlog(z)
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t32 = -t2/1024-t6/1024-t5*t8/256-t5*t12/128-t5*t1/512+t1*t17/256+t
     #1*t21/512-t1*t8/256-t1*t12/128-t21*t17/256+t17*t8/128
      t45 = t12**2
      t47 = t21**2
      t49 = t8**2
      t51 = 0.3141592653589793D1**2
      t53 = t17**2
      t55 = t21*t8/256+t12*t17/64+t12*t21/128-t12*t8/64+t5*t17/256+t5*t2
     #1/512-t45/64-t47/1024-t49/256+t51/1536-t53/256
      t57 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t61 = 1/t60
      t63 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t64 = 1/t63
      t65 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t66 = t65**2
      t67 = 1/t66
      t69 = f1(1.D0,z,0.D0)
      t72 = 6*t5
      t73 = 12*t17
      t74 = 6*t21
      t75 = 12*t8
      t76 = 24*t12
      t77 = 6*t1
      t78 = 1-l2
      t79 = dlog(t78)
      t83 = -1/t78
      t88 = t61*t64*t67*t69
      t91 = dlog(l3)
      t94 = 1/l3
      t103 = dlog(l1)
      t106 = 1/l1
      tncoff9 = (t32+t55)*t58*t61*t64*t67*t69+(-t72+t73+t74-t75-t76-t77-
     #6*t79)*t83*t58*t88/3072-(-t72+t73+t74-t75-t76-t77-12*t91)*t94*t58*
     #t88/1536+t94*t83*t58*t88/256-(-t72+t73+t74-t75-t76-t77-6*t103)*t10
     #6*t58*t88/3072+t106*t83*t58*t88/512-t94*t106*t58*t88/256
      t1 = -l3+2
      t2 = dlog(t1)
      t4 = dlog(z)
      t7 = 1-l3/2
      t8 = f2(t7,z)
      t9 = dlog(t8)
      t12 = dlog(1-z)
      t14 = dlog(l3)
      t16 = dlog(l4)
      t19 = dlog(1-l4)
      t23 = -1/t1
      t26 = s34r2(z,0.D0,1.D0,t7,l4)
      t28 = 1/l3/t26
      t30 = s2vr2(z,0.D0,1.D0,t7,l4)
      t32 = s23r2(z,0.D0,1.D0,t7,l4)
      t35 = s1vr2(z,0.D0,1.D0,t7,l4)
      t36 = t35**2
      t38 = f1(t7,z,0.D0)
      t40 = 1/t30/t32/t36*t38
      tncoff10 = -(-12*t2+12*t4+6*t9-24*t12-12*t14-6*t16-6*t19)*t23*t28*
     #t40/768+t23/(l2-1)*t28*t40/128-t23/l1*t28*t40/128
      t1 = dlog(l4)
      t3 = 1-l2
      t4 = dlog(t3)
      t7 = dlog(1-l4)
      t10 = dlog(1-z)
      t12 = f2(1.D0,z)
      t13 = dlog(t12)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = dlog(l2)
      t23 = -1/t3
      t25 = s34r2(z,0.D0,l2,1.D0,l4)
      t26 = 1/t25
      t28 = s2vr2(z,0.D0,l2,1.D0,l4)
      t30 = s23r2(z,0.D0,l2,1.D0,l4)
      t33 = s1vr2(z,0.D0,l2,1.D0,l4)
      t34 = t33**2
      t36 = f1(1.D0,z,0.D0)
      t38 = 1/t28/t30/t34*t36
      tncoff11 = -(-6*t1-6*t4-6*t7-24*t10+6*t13-12*t15+12*t17-6*t19)*t23
     #*t26*t38/3072-1/l3*t23*t26*t38/256-1/l1*t23*t26*t38/512
      t8 = 1-l3/2
      t9 = s34r2(z,0.D0,l2,t8,l4)
      t13 = s2vr2(z,0.D0,l2,t8,l4)
      t15 = s23r2(z,0.D0,l2,t8,l4)
      t18 = s1vr2(z,0.D0,l2,t8,l4)
      t19 = t18**2
      t21 = f1(t8,z,0.D0)
      tncoff12 = -1/(l3-2)/(l2-1)/l3/t9/t13/t15/t19*t21/128
      t1 = dlog(l1)
      t4 = dlog(1-l4)
      t6 = dlog(l4)
      t8 = dlog(z)
      t11 = dlog(1-l1)
      t13 = dlog(2.D0)
      t15 = f3(1.D0,z)
      t18 = dlog(1-l1*t15)
      t21 = dlog(1-z)
      t23 = f2(1.D0,z)
      t24 = dlog(t23)
      t27 = 1/l1
      t29 = s34r2(z,l1,1.D0,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,l1,1.D0,1.D0,l4)
      t34 = s23r2(z,l1,1.D0,1.D0,l4)
      t37 = s1vr2(z,l1,1.D0,1.D0,l4)
      t38 = t37**2
      t40 = f1(1.D0,z,l1)
      t42 = 1/t32/t34/t38*t40
      tncoff13 = (-6*t1-6*t4-6*t6+12*t8-6*t11-12*t13-6*t18-24*t21+6*t24)
     #*t27*t30*t42/3072-t27/(l2-1)*t30*t42/512+1/l3*t27*t30*t42/256
      t7 = 1-l3/2
      t8 = s34r2(z,l1,1.D0,t7,l4)
      t12 = s2vr2(z,l1,1.D0,t7,l4)
      t14 = s23r2(z,l1,1.D0,t7,l4)
      t17 = s1vr2(z,l1,1.D0,t7,l4)
      t18 = t17**2
      t20 = f1(t7,z,l1)
      tncoff14 = 1/(l3-2)/l1/l3/t8/t12/t14/t18*t20/128
      t5 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t10 = s23r2(z,l1,l2,1.D0,l4)
      t13 = s1vr2(z,l1,l2,1.D0,l4)
      t14 = t13**2
      t16 = f1(1.D0,z,l1)
      tncoff15 = 1/l1/(l2-1)/t5/t8/t10/t14*t16/512
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(z)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(2.D0)
      t14 = dlog(1-z)
      t17 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t26 = t25**2
      t29 = f1(1.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff16 = (t2/512+t4/512-t6/256-t9/512+t11/256+t14/128)*t18*t21*t
     #30+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/
     #512
      t6 = 1-l3/2
      t7 = s34r2(z,0.D0,1.D0,t6,l4)
      t10 = s2vr2(z,0.D0,1.D0,t6,l4)
      t12 = s23r2(z,0.D0,1.D0,t6,l4)
      t15 = s1vr2(z,0.D0,1.D0,t6,l4)
      t16 = t15**2
      t18 = f1(t6,z,0.D0)
      tncoff17 = -1/(l3-2)/l3/t7/t10/t12/t16*t18/128
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t9 = s23r2(z,0.D0,l2,1.D0,l4)
      t11 = s1vr2(z,0.D0,l2,1.D0,l4)
      t12 = t11**2
      t15 = f1(1.D0,z,0.D0)
      tncoff18 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,1.D0,l4)
      t5 = s2vr2(z,l1,1.D0,1.D0,l4)
      t8 = s23r2(z,l1,1.D0,1.D0,l4)
      t10 = s1vr2(z,l1,1.D0,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,l1)
      tncoff19 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,z,0.D0)
      tncoff20 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e0(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ACnum1e0(var,ndim,z,ux,l1,1.D0,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ACnum1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ACnum1e0(var,ndim,z,ux,l1,l2,1-l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ACnum1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ACnum1e1(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ACnum1e1(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ACnum1e1(var,ndim,z,ux,l1,1.D0,1-l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ACnum1e1(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ACnum1e2(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ACnum1e2(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ACnum1e3(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      AC2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function AC2e1int(var,ndim,Tcoeff)

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
      double precision ACnum1e0
      external ACnum1e0
      double precision ACnum1e1
      external ACnum1e1
      double precision ACnum1e2
      external ACnum1e2
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

      t1 = dlog(l4)
      t2 = t1**2
      t5 = dlog(1-l4)
      t6 = t5**2
      t8 = dlog(2.D0)
      t12 = dlog(1-z)
      t17 = dlog(z)
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t32 = -t2/1024-t6/1024-t5*t8/256-t5*t12/128-t5*t1/512+t1*t17/256+t
     #1*t21/512-t1*t8/256-t1*t12/128-t21*t17/256+t17*t8/128
      t45 = t12**2
      t47 = t21**2
      t49 = t8**2
      t51 = 0.3141592653589793D1**2
      t53 = t17**2
      t55 = t21*t8/256+t12*t17/64+t12*t21/128-t12*t8/64+t5*t17/256+t5*t2
     #1/512-t45/64-t47/1024-t49/256+t51/1536-t53/256
      t57 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t61 = 1/t60
      t63 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t64 = 1/t63
      t65 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t66 = t65**2
      t67 = 1/t66
      t69 = f1(1.D0,z,0.D0)
      t72 = 2*t17
      t73 = 4*t12
      t74 = 2*t8
      t75 = 1-l2
      t76 = dlog(t75)
      t79 = -1/t75
      t84 = t61*t64*t67*t69
      t87 = dlog(l3)
      t90 = 1/l3
      t99 = dlog(l1)
      t101 = 1/l1
      tncoff1 = (t32+t55)*t58*t61*t64*t67*t69-(t1-t21-t72+t73+t74+t76+t5
     #)*t79*t58*t84/512+(t5+t1+2*t87-t21+t74+t73-t72)*t90*t58*t84/256+t9
     #0*t79*t58*t84/256+(t99+t74+t73-t21-t72+t1+t5)*t101*t58*t84/512+t10
     #1*t79*t58*t84/512-t90*t101*t58*t84/256
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l3)
      t7 = dlog(1-z)
      t10 = 1-l3/2
      t11 = f2(t10,z)
      t12 = dlog(t11)
      t13 = dlog(z)
      t15 = -l3+2
      t16 = dlog(t15)
      t20 = -1/t15
      t23 = s34r2(z,0.D0,1.D0,t10,l4)
      t25 = 1/l3/t23
      t27 = s2vr2(z,0.D0,1.D0,t10,l4)
      t29 = s23r2(z,0.D0,1.D0,t10,l4)
      t32 = s1vr2(z,0.D0,1.D0,t10,l4)
      t33 = t32**2
      t35 = f1(t10,z,0.D0)
      t37 = 1/t27/t29/t33*t35
      tncoff2 = (t2+t3+2*t4+4*t7-t12-2*t13+2*t16)*t20*t25*t37/128+t20/(l
     #2-1)*t25*t37/128-t20/l1*t25*t37/128
      t1 = dlog(l4)
      t2 = dlog(l2)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = dlog(z)
      t7 = dlog(2.D0)
      t10 = dlog(1-z)
      t12 = 1-l2
      t13 = dlog(t12)
      t15 = dlog(1-l4)
      t18 = -1/t12
      t20 = s34r2(z,0.D0,l2,1.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,0.D0,l2,1.D0,l4)
      t25 = s23r2(z,0.D0,l2,1.D0,l4)
      t28 = s1vr2(z,0.D0,l2,1.D0,l4)
      t29 = t28**2
      t31 = f1(1.D0,z,0.D0)
      t33 = 1/t23/t25/t29*t31
      tncoff3 = (t1+t2-t4-2*t5+2*t7+4*t10+t13+t15)*t18*t21*t33/512-1/l3*
     #t18*t21*t33/256-1/l1*t18*t21*t33/512
      t8 = 1-l3/2
      t9 = s34r2(z,0.D0,l2,t8,l4)
      t13 = s2vr2(z,0.D0,l2,t8,l4)
      t15 = s23r2(z,0.D0,l2,t8,l4)
      t18 = s1vr2(z,0.D0,l2,t8,l4)
      t19 = t18**2
      t21 = f1(t8,z,0.D0)
      tncoff4 = -1/(l3-2)/(l2-1)/l3/t9/t13/t15/t19*t21/128
      t1 = dlog(l1)
      t2 = f3(1.D0,z)
      t5 = dlog(1-l1*t2)
      t6 = dlog(2.D0)
      t9 = dlog(1-z)
      t12 = dlog(1-l1)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t15 = dlog(z)
      t17 = dlog(l4)
      t19 = dlog(1-l4)
      t21 = 1/l1
      t23 = s34r2(z,l1,1.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s2vr2(z,l1,1.D0,1.D0,l4)
      t28 = s23r2(z,l1,1.D0,1.D0,l4)
      t31 = s1vr2(z,l1,1.D0,1.D0,l4)
      t32 = t31**2
      t34 = f1(1.D0,z,l1)
      t36 = 1/t26/t28/t32*t34
      tncoff5 = -(t1+t5+2*t6+4*t9+t12-t14-2*t15+t17+t19)*t21*t24*t36/512
     #-t21/(l2-1)*t24*t36/512+1/l3*t21*t24*t36/256
      t7 = 1-l3/2
      t8 = s34r2(z,l1,1.D0,t7,l4)
      t12 = s2vr2(z,l1,1.D0,t7,l4)
      t14 = s23r2(z,l1,1.D0,t7,l4)
      t17 = s1vr2(z,l1,1.D0,t7,l4)
      t18 = t17**2
      t20 = f1(t7,z,l1)
      tncoff6 = 1/(l3-2)/l1/l3/t8/t12/t14/t18*t20/128
      t5 = s34r2(z,l1,l2,1.D0,l4)
      t8 = s2vr2(z,l1,l2,1.D0,l4)
      t10 = s23r2(z,l1,l2,1.D0,l4)
      t13 = s1vr2(z,l1,l2,1.D0,l4)
      t14 = t13**2
      t16 = f1(1.D0,z,l1)
      tncoff7 = 1/l1/(l2-1)/t5/t8/t10/t14*t16/512
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(z)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(2.D0)
      t14 = dlog(1-z)
      t17 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t26 = t25**2
      t29 = f1(1.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff8 = (t2/512+t4/512-t6/256-t9/512+t11/256+t14/128)*t18*t21*t3
     #0+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/5
     #12
      t6 = 1-l3/2
      t7 = s34r2(z,0.D0,1.D0,t6,l4)
      t10 = s2vr2(z,0.D0,1.D0,t6,l4)
      t12 = s23r2(z,0.D0,1.D0,t6,l4)
      t15 = s1vr2(z,0.D0,1.D0,t6,l4)
      t16 = t15**2
      t18 = f1(t6,z,0.D0)
      tncoff9 = -1/(l3-2)/l3/t7/t10/t12/t16*t18/128
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t9 = s23r2(z,0.D0,l2,1.D0,l4)
      t11 = s1vr2(z,0.D0,l2,1.D0,l4)
      t12 = t11**2
      t15 = f1(1.D0,z,0.D0)
      tncoff10 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,1.D0,l4)
      t5 = s2vr2(z,l1,1.D0,1.D0,l4)
      t8 = s23r2(z,l1,1.D0,1.D0,l4)
      t10 = s1vr2(z,l1,1.D0,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,l1)
      tncoff11 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,z,0.D0)
      tncoff12 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,0.D0,l2,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e0(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ACnum1e0(var,ndim,z,ux,l1,1.D0,1-l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ACnum1e0(var,ndim,z,ux,l1,l2,1.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ACnum1e1(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ACnum1e1(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      AC2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function AC2e2int(var,ndim,Tcoeff)

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
      double precision ACnum1e0
      external ACnum1e0
      double precision ACnum1e1
      external ACnum1e1
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

      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(z)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t11 = dlog(2.D0)
      t14 = dlog(1-z)
      t17 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t26 = t25**2
      t29 = f1(1.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff1 = (t2/512+t4/512-t6/256-t9/512+t11/256+t14/128)*t18*t21*t3
     #0+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/5
     #12
      t6 = 1-l3/2
      t7 = s34r2(z,0.D0,1.D0,t6,l4)
      t10 = s2vr2(z,0.D0,1.D0,t6,l4)
      t12 = s23r2(z,0.D0,1.D0,t6,l4)
      t15 = s1vr2(z,0.D0,1.D0,t6,l4)
      t16 = t15**2
      t18 = f1(t6,z,0.D0)
      tncoff2 = -1/(l3-2)/l3/t7/t10/t12/t16*t18/128
      t3 = s34r2(z,0.D0,l2,1.D0,l4)
      t6 = s2vr2(z,0.D0,l2,1.D0,l4)
      t9 = s23r2(z,0.D0,l2,1.D0,l4)
      t11 = s1vr2(z,0.D0,l2,1.D0,l4)
      t12 = t11**2
      t15 = f1(1.D0,z,0.D0)
      tncoff3 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,1.D0,l4)
      t5 = s2vr2(z,l1,1.D0,1.D0,l4)
      t8 = s23r2(z,l1,1.D0,1.D0,l4)
      t10 = s1vr2(z,l1,1.D0,1.D0,l4)
      t11 = t10**2
      t14 = f1(1.D0,z,l1)
      tncoff4 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,1.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,1.D0,l4)
      t9 = t8**2
      t12 = f1(1.D0,z,0.D0)
      tncoff5 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,l1,1.D0,1.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      AC2e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

