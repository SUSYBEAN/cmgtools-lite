CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function AC1e0int(var,ndim,Tcoeff)

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

      t1 = dlog(l4)
      t2 = t1**2
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t7 = dlog(2.D0)
      t10 = 0.3141592653589793D1**2
      t14 = dlog(1-l4)
      t18 = dlog(1-z)
      t19 = t4**2
      t22 = dlog(z)
      t27 = t7**2
      t32 = t22**2
      t43 = t18**2
      t46 = -t2*t4/1024+t2*t7/512-t1*t10/1536-t14*t10/1536+t18*t19/256-t
     #2*t22/512+t2*t18/256-t27*t4/256+t14*t19/1024+t14*t32/256+t19*t7/51
     #2-t19*t22/512-t18*t10/384+t18*t32/64-t43*t22/32
      t53 = t14**2
      t64 = t1*t18
      t79 = t14*t4
      t82 = t4*t10/1536-t43*t4/64+t43*t7/32+t53*t18/256-t4*t32/256+t10*t
     #22/768+t53*t1/1024+t18*t27/64-t64*t4/128-t64*t22/64+t53*t7/512+t53
     #*t14/3072+t2*t1/3072-t32*t22/384+t43*t18/48-t79*t7/256
      t84 = t14*t22
      t87 = t14*t18
      t90 = t14*t1
      t95 = t18*t4
      t98 = t18*t7
      t101 = t18*t22
      t108 = t1*t4
      t115 = t1*t22
      t122 = t4*t22
      t125 = -t84*t7/128+t87*t7/64-t90*t4/512+t79*t22/256+t95*t22/64-t98
     #*t4/64-t101*t7/32+t90*t18/128+t90*t7/256+t108*t22/256-t87*t22/64-t
     #108*t7/256-t115*t7/128-t87*t4/128-t90*t22/256+t122*t7/128
      t157 = t64*t7/64-t19*t4/3072+t27*t7/384+z3/192+t14*t43/64+t14*t27/
     #256-t10*t7/768-t22*t27/128+t1*t32/256+t1*t19/1024+t1*t43/64+t1*t27
     #/256+t32*t7/128-t53*t4/1024-t53*t22/512+t14*t2/1024
      t160 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t161 = 1/t160
      t163 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t164 = 1/t163
      t166 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t167 = 1/t166
      t168 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t169 = t168**2
      t170 = 1/t169
      t172 = f1(0.D0,z,0.D0)
      t175 = 3*t2
      t176 = 1-l2
      t177 = dlog(t176)
      t184 = 6*t79
      t186 = 12*t14*t7
      t191 = 48*t101
      t195 = 24*t22*t7
      t197 = 12*t4*t7
      t198 = 48*t98
      t199 = 12*t84
      t200 = t175+6*t14*t177+6*t1*t177+24*t177*t18-t184+t186-6*t177*t4+1
     #2*t177*t7-t191-12*t177*t22-t195-t197+t198-t199
      t201 = 24*t87
      t202 = 24*t95
      t203 = 12*t115
      t204 = 6*t90
      t205 = 6*t108
      t207 = 12*t1*t7
      t208 = 12*t122
      t209 = 24*t64
      t210 = 48*t43
      t211 = 12*t27
      t212 = 2*t10
      t213 = t177**2
      t215 = 12*t32
      t216 = 3*t19
      t217 = 3*t53
      t218 = t201-t202-t203+t204-t205+t207+t208+t209+t210+t211-t212+3*t2
     #13+t215+t216+t217
      t221 = -1/t176
      t226 = t164*t167*t170*t172
      t229 = dlog(l3)
      t230 = t229**2
      t232 = t175+12*t230-t184+t186-t191-t195-t197+t198-t199+t201-t202-t
     #203+t204-t205
      t245 = t207+t208+t209-24*t229*t22+24*t229*t7+12*t1*t229-12*t229*t4
     #+48*t229*t18+12*t14*t229+t210+t211-t212+t215+t216+t217
      t247 = 1/l3
      t252 = 2*t229
      t253 = 2*t22
      t254 = 4*t18
      t255 = 2*t7
      t258 = t221*t161
      t262 = dlog(l1)
      t263 = t262**2
      t271 = t175-t184+t186-t191-t195-t197+t198-t199+t201-t202+3*t263+6*
     #t1*t262+6*t14*t262-12*t262*t22
      t278 = 12*t262*t7+24*t262*t18-6*t262*t4-t203+t204-t205+t207+t208+t
     #209+t210+t211-t212+t215+t216+t217
      t280 = 1/l1
      t292 = t280*t161
      tncoff1 = (t46+t82+t125+t157)*t161*t164*t167*t170*t172+(t200+t218)
     #*t221*t161*t226/3072-(t232+t245)*t247*t161*t226/1536-(t252-t253+t2
     #54-t4+t255+t177+t1+t14)*t247*t258*t226/256-(t271+t278)*t280*t161*t
     #226/3072-(t14+t177+t262+t254+t255-t4-t253+t1)*t280*t258*t226/512+(
     #t252+t262+t255-t4-t253+t254+t14+t1)*t247*t292*t226/256+t247*t221*t
     #292*t226/256
      t1 = dlog(l3)
      t2 = l3/2
      t3 = f2(t2,z)
      t4 = dlog(t3)
      t7 = dlog(z)
      t10 = -l3+2
      t11 = dlog(t10)
      t15 = dlog(1-l4)
      t18 = dlog(l4)
      t25 = t7**2
      t29 = t18**2
      t31 = 0.3141592653589793D1**2
      t34 = dlog(1-z)
      t35 = t34**2
      t37 = t15**2
      t41 = 12*t1*t4+24*t1*t7+12*t4*t11-12*t15*t11-12*t18*t11+6*t18*t4+6
     #*t15*t4-12*t25-12*t1*t18-3*t29+2*t31-48*t35-3*t37-6*t15*t18
      t44 = t1**2
      t52 = t4**2
      t62 = t11**2
      t72 = 24*t7*t11-12*t44-48*t1*t34+12*t15*t7-24*t15*t34-3*t52-24*t18
     #*t34-12*t15*t1+12*t18*t7+24*t34*t4-12*t62-24*t1*t11-12*t4*t7-48*t3
     #4*t11+48*t34*t7
      t75 = -1/t10
      t77 = 1/l3
      t78 = s34r2(z,0.D0,1.D0,t2,l4)
      t79 = 1/t78
      t82 = s2vr2(z,0.D0,1.D0,t2,l4)
      t83 = 1/t82
      t84 = s23r2(z,0.D0,1.D0,t2,l4)
      t85 = 1/t84
      t87 = s1vr2(z,0.D0,1.D0,t2,l4)
      t88 = t87**2
      t89 = 1/t88
      t90 = f1(t2,z,0.D0)
      t95 = 2*t1
      t96 = 2*t11
      t97 = 4*t34
      t98 = 2*t7
      t99 = 1-l2
      t100 = dlog(t99)
      t104 = -1/t99
      t107 = t79*t83
      t110 = t107*t85*t89*t90
      t113 = dlog(l1)
      t116 = 1/l1
      tncoff2 = (t41+t72)*t75*t77*t79*t83*t85*t89*t90/768-(t95+t96+t97-t
     #98+t100-t4+t18+t15)*t75*t104*t77*t110/128-(-t95-t96-t113+t4+t98-t9
     #7-t15-t18)*t75*t116*t77*t110/128+t75*t116*t107*t85*t104*t77*t89*t9
     #0/128
      t2 = dlog(1-l4)
      t3 = t2**2
      t5 = dlog(2.D0)
      t6 = t5**2
      t8 = 1-l2
      t9 = dlog(t8)
      t10 = t9**2
      t12 = dlog(l2)
      t13 = t12**2
      t15 = dlog(l4)
      t18 = f2(0.D0,z)
      t19 = dlog(t18)
      t20 = dlog(z)
      t32 = dlog(1-z)
      t49 = 3*t3+12*t6+3*t10+3*t13+12*t15*t5+12*t19*t20+6*t2*t12+12*t2*t
     #5-6*t15*t19+12*t12*t5+24*t12*t32-12*t12*t20-6*t12*t19-6*t2*t19+6*t
     #2*t15+24*t9*t32-48*t32*t20-12*t9*t20
      t76 = t19**2
      t80 = t20**2
      t82 = t15**2
      t84 = 0.3141592653589793D1**2
      t86 = t32**2
      t88 = -12*t19*t5+48*t32*t5-12*t2*t20+24*t2*t32-24*t32*t19+12*t9*t5
     #-24*t20*t5-6*t9*t19+6*t15*t9+6*t2*t9+24*t15*t32+6*t15*t12+6*t9*t12
     #+3*t76-12*t15*t20+12*t80+3*t82-2*t84+48*t86
      t91 = -1/t8
      t93 = s34r2(z,0.D0,l2,0.D0,l4)
      t94 = 1/t93
      t96 = s2vr2(z,0.D0,l2,0.D0,l4)
      t98 = s23r2(z,0.D0,l2,0.D0,l4)
      t101 = s1vr2(z,0.D0,l2,0.D0,l4)
      t102 = t101**2
      t104 = f1(0.D0,z,0.D0)
      t106 = 1/t96/t98/t102*t104
      t109 = dlog(l3)
      t111 = 4*t32
      t112 = 2*t20
      t113 = 2*t5
      t115 = 1/l3
      t117 = t91*t94
      t121 = dlog(l1)
      t123 = 1/l1
      tncoff3 = -(t49+t88)*t91*t94*t106/3072+(2*t109+t9+t111-t112+t113+t
     #12-t19+t15+t2)*t115*t117*t106/256+(t2+t9+t12+t121+t111+t113-t19-t1
     #12+t15)*t123*t117*t106/512-t115*t91*t123*t94*t106/256
      t1 = dlog(l3)
      t3 = -l3+2
      t4 = dlog(t3)
      t7 = dlog(1-z)
      t9 = dlog(l2)
      t10 = dlog(z)
      t12 = l3/2
      t13 = f2(t12,z)
      t14 = dlog(t13)
      t15 = 1-l2
      t16 = dlog(t15)
      t17 = dlog(l4)
      t19 = dlog(1-l4)
      t22 = -1/t3
      t27 = -1/t15/l3
      t29 = s34r2(z,0.D0,l2,t12,l4)
      t30 = 1/t29
      t31 = s2vr2(z,0.D0,l2,t12,l4)
      t32 = 1/t31
      t34 = s23r2(z,0.D0,l2,t12,l4)
      t35 = 1/t34
      t36 = s1vr2(z,0.D0,l2,t12,l4)
      t37 = t36**2
      t38 = 1/t37
      t40 = f1(t12,z,0.D0)
      tncoff4 = (2*t1+2*t4+4*t7+t9-2*t10-t14+t16+t17+t19)*t22*t27*t30*t3
     #2*t35*t38*t40/128-t22*t32*t30*t38*t35/l1*t27*t40/128
      t1 = dlog(l1)
      t2 = f3(0.D0,z)
      t5 = dlog(1-l1*t2)
      t9 = dlog(1-l4)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t14 = dlog(2.D0)
      t17 = dlog(l4)
      t19 = dlog(1-l1)
      t32 = dlog(z)
      t36 = dlog(1-z)
      t39 = 6*t1*t5-6*t9*t11+12*t5*t14+6*t17*t19-6*t19*t11+6*t9*t19+6*t9
     #*t17-6*t17*t11+12*t17*t14+12*t11*t32+24*t17*t36
      t64 = 12*t9*t14-48*t36*t32-24*t14*t32-12*t11*t14+48*t36*t14-12*t9*
     #t32+24*t9*t36-24*t36*t11+6*t17*t1+6*t9*t1-12*t1*t32+12*t14*t1
      t70 = t11**2
      t72 = t32**2
      t74 = t17**2
      t76 = 0.3141592653589793D1**2
      t78 = t36**2
      t80 = t9**2
      t82 = t14**2
      t88 = 24*t1*t36-6*t1*t11+3*t70+12*t72+3*t74-2*t76+48*t78+3*t80+12*
     #t82-12*t17*t32+6*t19*t5
      t105 = t1**2
      t107 = t19**2
      t111 = t5**2
      t113 = -6*t5*t11+24*t36*t5+12*t14*t19+6*t1*t19-12*t32*t19-12*t32*t
     #5+24*t36*t19+6*t17*t5+3*t105+3*t107+6*t9*t5+3*t111
      t116 = 1/l1
      t118 = s34r2(z,l1,1.D0,0.D0,l4)
      t119 = 1/t118
      t121 = s2vr2(z,l1,1.D0,0.D0,l4)
      t123 = s23r2(z,l1,1.D0,0.D0,l4)
      t126 = s1vr2(z,l1,1.D0,0.D0,l4)
      t127 = t126**2
      t129 = f1(0.D0,z,l1)
      t131 = 1/t121/t123/t127*t129
      t134 = 1-l2
      t135 = dlog(t134)
      t136 = 4*t36
      t137 = 2*t14
      t138 = 2*t32
      t142 = -1/t134
      t147 = dlog(l3)
      t150 = 1/l3
      t152 = t116*t119
      tncoff5 = (t39+t64+t88+t113)*t116*t119*t131/3072+(t9+t135+t1+t136+
     #t5+t137-t11-t138+t19+t17)*t116*t142*t119*t131/512-(2*t147+t5+t19+t
     #136+t1+t137-t138-t11+t9+t17)*t150*t152*t131/256-t150*t142*t152*t13
     #1/256
      t1 = dlog(l3)
      t4 = dlog(1-l1)
      t5 = dlog(z)
      t7 = l3/2
      t8 = f2(t7,z)
      t9 = dlog(t8)
      t10 = -l3+2
      t11 = dlog(t10)
      t14 = dlog(1-z)
      t16 = f3(t7,z)
      t19 = dlog(1-l1*t16)
      t20 = dlog(l1)
      t22 = dlog(1-l4)
      t23 = dlog(l4)
      t26 = -1/t10
      t28 = 1/l1
      t29 = 1/l3
      t32 = s34r2(z,l1,1.D0,t7,l4)
      t33 = 1/t32
      t34 = s2vr2(z,l1,1.D0,t7,l4)
      t35 = 1/t34
      t37 = s23r2(z,l1,1.D0,t7,l4)
      t38 = 1/t37
      t39 = s1vr2(z,l1,1.D0,t7,l4)
      t40 = t39**2
      t41 = 1/t40
      t43 = f1(t7,z,l1)
      tncoff6 = (-2*t1-t4+2*t5+t9-2*t11-4*t14-t19-t20-t22-t23)*t26*t28*t
     #29*t33*t35*t38*t41*t43/128-t26*t28*t38*t33*t35/(l2-1)*t29*t41*t43/
     #128
      t2 = dlog(1-l4)
      t3 = 1-l2
      t4 = dlog(t3)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t8 = dlog(1-l1)
      t9 = dlog(l2)
      t10 = dlog(l1)
      t11 = dlog(2.D0)
      t14 = dlog(1-z)
      t16 = dlog(z)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t22 = dlog(l4)
      t23 = t2+t4-t6+t8+t9+t10+2*t11+4*t14-2*t16+t21+t22
      t24 = 1/l1
      t27 = -1/t3
      t28 = s34r2(z,l1,l2,0.D0,l4)
      t29 = 1/t28
      t32 = s2vr2(z,l1,l2,0.D0,l4)
      t34 = s23r2(z,l1,l2,0.D0,l4)
      t37 = s1vr2(z,l1,l2,0.D0,l4)
      t38 = t37**2
      t40 = f1(0.D0,z,l1)
      t42 = 1/t32/t34/t38*t40
      tncoff7 = -t23*t24*t27*t29*t42/512+1/l3*t27*t24*t29*t42/256
      t1 = l3/2
      t2 = s2vr2(z,l1,l2,t1,l4)
      t4 = s34r2(z,l1,l2,t1,l4)
      t7 = s1vr2(z,l1,l2,t1,l4)
      t8 = t7**2
      t14 = s23r2(z,l1,l2,t1,l4)
      t22 = f1(t1,z,l1)
      tncoff8 = 1/t2/t4/t8/(l3-2)/t14/l1/(l2-1)/l3*t22/128
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t7 = dlog(1-z)
      t10 = dlog(2.D0)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t17 = dlog(z)
      t26 = t7**2
      t28 = t10**2
      t30 = t14**2
      t32 = -t2*t3/512-t2*t7/128-t2*t10/256+t2*t14/512+t2*t17/256-t7*t10
     #/64+t7*t14/128+t7*t17/64-t26/64-t28/256-t30/1024
      t33 = t2**2
      t35 = t3**2
      t37 = 0.3141592653589793D1**2
      t39 = t17**2
      t55 = -t33/1024-t35/1024+t37/1536-t39/256+t14*t10/256+t17*t10/128-
     #t14*t17/256+t3*t17/256+t3*t14/512-t3*t10/256-t3*t7/128
      t57 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t61 = 1/t60
      t63 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t64 = 1/t63
      t65 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t66 = t65**2
      t67 = 1/t66
      t69 = f1(0.D0,z,0.D0)
      t72 = 12*t17
      t73 = 6*t3
      t74 = 6*t2
      t75 = 24*t7
      t76 = 12*t10
      t77 = 6*t14
      t78 = 1-l2
      t79 = dlog(t78)
      t83 = -1/t78
      t88 = t61*t64*t67*t69
      t91 = dlog(l3)
      t94 = 1/l3
      t103 = dlog(l1)
      t106 = 1/l1
      tncoff9 = (t32+t55)*t58*t61*t64*t67*t69+(t72-t73-t74-t75-t76+t77-6
     #*t79)*t83*t58*t88/3072-(t72-t73-t74-t75-t76+t77-12*t91)*t94*t58*t8
     #8/1536+t94*t83*t58*t88/256-(t72-t73-t74-t75-t76+t77-6*t103)*t106*t
     #58*t88/3072+t106*t83*t58*t88/512-t94*t106*t58*t88/256
      t2 = dlog(1-l4)
      t4 = dlog(z)
      t6 = l3/2
      t7 = f2(t6,z)
      t8 = dlog(t7)
      t10 = dlog(l3)
      t12 = dlog(l4)
      t15 = dlog(1-z)
      t17 = -l3+2
      t18 = dlog(t17)
      t22 = -1/t17
      t25 = s34r2(z,0.D0,1.D0,t6,l4)
      t27 = 1/l3/t25
      t29 = s2vr2(z,0.D0,1.D0,t6,l4)
      t31 = s23r2(z,0.D0,1.D0,t6,l4)
      t34 = s1vr2(z,0.D0,1.D0,t6,l4)
      t35 = t34**2
      t37 = f1(t6,z,0.D0)
      t39 = 1/t29/t31/t35*t37
      tncoff10 = (6*t2-12*t4-6*t8+12*t10+6*t12+24*t15+12*t18)*t22*t27*t3
     #9/768+t22/(l2-1)*t27*t39/128-t22/l1*t27*t39/128
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = 1-l2
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(l2)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t19 = dlog(1-z)
      t23 = -1/t6
      t25 = s34r2(z,0.D0,l2,0.D0,l4)
      t26 = 1/t25
      t28 = s2vr2(z,0.D0,l2,0.D0,l4)
      t30 = s23r2(z,0.D0,l2,0.D0,l4)
      t33 = s1vr2(z,0.D0,l2,0.D0,l4)
      t34 = t33**2
      t36 = f1(0.D0,z,0.D0)
      t38 = 1/t28/t30/t34*t36
      tncoff11 = -(-6*t1-6*t4-6*t7+12*t9-6*t11+6*t14-12*t16-24*t19)*t23*
     #t26*t38/3072-1/l3*t23*t26*t38/256-1/l1*t23*t26*t38/512
      t7 = l3/2
      t8 = s34r2(z,0.D0,l2,t7,l4)
      t12 = s2vr2(z,0.D0,l2,t7,l4)
      t14 = s23r2(z,0.D0,l2,t7,l4)
      t17 = s1vr2(z,0.D0,l2,t7,l4)
      t18 = t17**2
      t20 = f1(t7,z,0.D0)
      tncoff12 = -1/(l3-2)/(l2-1)/l3/t8/t12/t14/t18*t20/128
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = dlog(l1)
      t9 = dlog(1-l1)
      t12 = dlog(1-z)
      t14 = f3(0.D0,z)
      t17 = dlog(1-l1*t14)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t22 = dlog(2.D0)
      t24 = dlog(z)
      t27 = 1/l1
      t29 = s34r2(z,l1,1.D0,0.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,l1,1.D0,0.D0,l4)
      t34 = s23r2(z,l1,1.D0,0.D0,l4)
      t37 = s1vr2(z,l1,1.D0,0.D0,l4)
      t38 = t37**2
      t40 = f1(0.D0,z,l1)
      t42 = 1/t32/t34/t38*t40
      tncoff13 = (-6*t2-6*t4-6*t6-6*t9-24*t12-6*t17+6*t20-12*t22+12*t24)
     #*t27*t30*t42/3072-t27/(l2-1)*t30*t42/512+1/l3*t27*t30*t42/256
      t6 = l3/2
      t7 = s34r2(z,l1,1.D0,t6,l4)
      t11 = s2vr2(z,l1,1.D0,t6,l4)
      t13 = s23r2(z,l1,1.D0,t6,l4)
      t16 = s1vr2(z,l1,1.D0,t6,l4)
      t17 = t16**2
      t19 = f1(t6,z,l1)
      tncoff14 = 1/(l3-2)/l1/l3/t7/t11/t13/t17*t19/128
      t5 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t10 = s23r2(z,l1,l2,0.D0,l4)
      t13 = s1vr2(z,l1,l2,0.D0,l4)
      t14 = t13**2
      t16 = f1(0.D0,z,l1)
      tncoff15 = 1/l1/(l2-1)/t5/t8/t10/t14*t16/512
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t26 = t25**2
      t29 = f1(0.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff16 = (t2/128+t4/256-t7/512-t9/256+t11/512+t14/512)*t18*t21*t
     #30+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/
     #512
      t5 = l3/2
      t6 = s34r2(z,0.D0,1.D0,t5,l4)
      t9 = s2vr2(z,0.D0,1.D0,t5,l4)
      t11 = s23r2(z,0.D0,1.D0,t5,l4)
      t14 = s1vr2(z,0.D0,1.D0,t5,l4)
      t15 = t14**2
      t17 = f1(t5,z,0.D0)
      tncoff17 = -1/(l3-2)/l3/t6/t9/t11/t15*t17/128
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t9 = s23r2(z,0.D0,l2,0.D0,l4)
      t11 = s1vr2(z,0.D0,l2,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,0.D0)
      tncoff18 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,0.D0,l4)
      t5 = s2vr2(z,l1,1.D0,0.D0,l4)
      t8 = s23r2(z,l1,1.D0,0.D0,l4)
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff19 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,0.D0)
      tncoff20 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ACnum1e0(var,ndim,z,ux,l1,1.D0,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ACnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ACnum1e0(var,ndim,z,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ACnum1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ACnum1e1(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ACnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ACnum1e1(var,ndim,z,ux,l1,1.D0,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ACnum1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ACnum1e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ACnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ACnum1e3(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      AC1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function AC1e1int(var,ndim,Tcoeff)

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

      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t7 = dlog(1-z)
      t10 = dlog(2.D0)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t17 = dlog(z)
      t26 = t7**2
      t28 = t10**2
      t30 = t14**2
      t32 = -t2*t3/512-t2*t7/128-t2*t10/256+t2*t14/512+t2*t17/256-t7*t10
     #/64+t7*t14/128+t7*t17/64-t26/64-t28/256-t30/1024
      t33 = t2**2
      t35 = t3**2
      t37 = 0.3141592653589793D1**2
      t39 = t17**2
      t55 = -t33/1024-t35/1024+t37/1536-t39/256+t14*t10/256+t17*t10/128-
     #t14*t17/256+t3*t17/256+t3*t14/512-t3*t10/256-t3*t7/128
      t57 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t61 = 1/t60
      t63 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t64 = 1/t63
      t65 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t66 = t65**2
      t67 = 1/t66
      t69 = f1(0.D0,z,0.D0)
      t72 = 4*t7
      t73 = 2*t10
      t74 = 2*t17
      t75 = 1-l2
      t76 = dlog(t75)
      t79 = -1/t75
      t84 = t61*t64*t67*t69
      t87 = dlog(l3)
      t90 = 1/l3
      t99 = dlog(l1)
      t101 = 1/l1
      tncoff1 = (t32+t55)*t58*t61*t64*t67*t69-(t72+t73-t14-t74+t3+t76+t2
     #)*t79*t58*t84/512+(t3-t14+t73+t72-t74+2*t87+t2)*t90*t58*t84/256+t9
     #0*t79*t58*t84/256+(t3+t2+t99+t72+t73-t14-t74)*t101*t58*t84/512+t10
     #1*t79*t58*t84/512-t90*t101*t58*t84/256
      t1 = dlog(l4)
      t3 = dlog(1-z)
      t5 = -l3+2
      t6 = dlog(t5)
      t8 = l3/2
      t9 = f2(t8,z)
      t10 = dlog(t9)
      t11 = dlog(z)
      t13 = dlog(l3)
      t16 = dlog(1-l4)
      t19 = -1/t5
      t22 = s34r2(z,0.D0,1.D0,t8,l4)
      t24 = 1/l3/t22
      t26 = s2vr2(z,0.D0,1.D0,t8,l4)
      t28 = s23r2(z,0.D0,1.D0,t8,l4)
      t31 = s1vr2(z,0.D0,1.D0,t8,l4)
      t32 = t31**2
      t34 = f1(t8,z,0.D0)
      t36 = 1/t26/t28/t32*t34
      tncoff2 = -(-t1-4*t3-2*t6+t10+2*t11-2*t13-t16)*t19*t24*t36/128+t19
     #/(l2-1)*t24*t36/128-t19/l1*t24*t36/128
      t2 = dlog(1-z)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = dlog(z)
      t8 = dlog(l2)
      t9 = dlog(2.D0)
      t11 = dlog(l4)
      t12 = 1-l2
      t13 = dlog(t12)
      t15 = dlog(1-l4)
      t18 = -1/t12
      t20 = s34r2(z,0.D0,l2,0.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,0.D0,l2,0.D0,l4)
      t25 = s23r2(z,0.D0,l2,0.D0,l4)
      t28 = s1vr2(z,0.D0,l2,0.D0,l4)
      t29 = t28**2
      t31 = f1(0.D0,z,0.D0)
      t33 = 1/t23/t25/t29*t31
      tncoff3 = (4*t2-t5-2*t6+t8+2*t9+t11+t13+t15)*t18*t21*t33/512-1/l3*
     #t18*t21*t33/256-1/l1*t18*t21*t33/512
      t7 = l3/2
      t8 = s34r2(z,0.D0,l2,t7,l4)
      t12 = s2vr2(z,0.D0,l2,t7,l4)
      t14 = s23r2(z,0.D0,l2,t7,l4)
      t17 = s1vr2(z,0.D0,l2,t7,l4)
      t18 = t17**2
      t20 = f1(t7,z,0.D0)
      tncoff4 = -1/(l3-2)/(l2-1)/l3/t8/t12/t14/t18*t20/128
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l1)
      t5 = dlog(2.D0)
      t8 = dlog(1-z)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t15 = dlog(1-l1)
      t16 = f2(0.D0,z)
      t17 = dlog(t16)
      t18 = dlog(z)
      t21 = 1/l1
      t23 = s34r2(z,l1,1.D0,0.D0,l4)
      t24 = 1/t23
      t26 = s2vr2(z,l1,1.D0,0.D0,l4)
      t28 = s23r2(z,l1,1.D0,0.D0,l4)
      t31 = s1vr2(z,l1,1.D0,0.D0,l4)
      t32 = t31**2
      t34 = f1(0.D0,z,l1)
      t36 = 1/t26/t28/t32*t34
      tncoff5 = -(t1+t3+t4+2*t5+4*t8+t13+t15-t17-2*t18)*t21*t24*t36/512-
     #t21/(l2-1)*t24*t36/512+1/l3*t21*t24*t36/256
      t6 = l3/2
      t7 = s34r2(z,l1,1.D0,t6,l4)
      t11 = s2vr2(z,l1,1.D0,t6,l4)
      t13 = s23r2(z,l1,1.D0,t6,l4)
      t16 = s1vr2(z,l1,1.D0,t6,l4)
      t17 = t16**2
      t19 = f1(t6,z,l1)
      tncoff6 = 1/(l3-2)/l1/l3/t7/t11/t13/t17*t19/128
      t5 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t10 = s23r2(z,l1,l2,0.D0,l4)
      t13 = s1vr2(z,l1,l2,0.D0,l4)
      t14 = t13**2
      t16 = f1(0.D0,z,l1)
      tncoff7 = 1/l1/(l2-1)/t5/t8/t10/t14*t16/512
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t26 = t25**2
      t29 = f1(0.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff8 = (t2/128+t4/256-t7/512-t9/256+t11/512+t14/512)*t18*t21*t3
     #0+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/5
     #12
      t5 = l3/2
      t6 = s34r2(z,0.D0,1.D0,t5,l4)
      t9 = s2vr2(z,0.D0,1.D0,t5,l4)
      t11 = s23r2(z,0.D0,1.D0,t5,l4)
      t14 = s1vr2(z,0.D0,1.D0,t5,l4)
      t15 = t14**2
      t17 = f1(t5,z,0.D0)
      tncoff9 = -1/(l3-2)/l3/t6/t9/t11/t15*t17/128
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t9 = s23r2(z,0.D0,l2,0.D0,l4)
      t11 = s1vr2(z,0.D0,l2,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,0.D0)
      tncoff10 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,0.D0,l4)
      t5 = s2vr2(z,l1,1.D0,0.D0,l4)
      t8 = s23r2(z,l1,1.D0,0.D0,l4)
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff11 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,0.D0)
      tncoff12 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ACnum1e0(var,ndim,z,ux,l1,1.D0,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ACnum1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ACnum1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ACnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ACnum1e2(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      AC1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function AC1e2int(var,ndim,Tcoeff)

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

      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(l4)
      t14 = dlog(1-l4)
      t17 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t25 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t26 = t25**2
      t29 = f1(0.D0,z,0.D0)
      t30 = 1/t23/t26*t29
      tncoff1 = (t2/128+t4/256-t7/512-t9/256+t11/512+t14/512)*t18*t21*t3
     #0+1/(l2-1)*t18*t21*t30/512-1/l3*t18*t21*t30/256-1/l1*t18*t21*t30/5
     #12
      t5 = l3/2
      t6 = s34r2(z,0.D0,1.D0,t5,l4)
      t9 = s2vr2(z,0.D0,1.D0,t5,l4)
      t11 = s23r2(z,0.D0,1.D0,t5,l4)
      t14 = s1vr2(z,0.D0,1.D0,t5,l4)
      t15 = t14**2
      t17 = f1(t5,z,0.D0)
      tncoff2 = -1/(l3-2)/l3/t6/t9/t11/t15*t17/128
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t9 = s23r2(z,0.D0,l2,0.D0,l4)
      t11 = s1vr2(z,0.D0,l2,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,0.D0)
      tncoff3 = -1/(l2-1)/t3/t6/t9/t12*t15/512
      t2 = s34r2(z,l1,1.D0,0.D0,l4)
      t5 = s2vr2(z,l1,1.D0,0.D0,l4)
      t8 = s23r2(z,l1,1.D0,0.D0,l4)
      t10 = s1vr2(z,l1,1.D0,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,l1)
      tncoff4 = 1/l1/t2/t5/t8/t11*t14/512
      t1 = s34r2(z,0.D0,1.D0,0.D0,l4)
      t3 = s2vr2(z,0.D0,1.D0,0.D0,l4)
      t6 = s23r2(z,0.D0,1.D0,0.D0,l4)
      t8 = s1vr2(z,0.D0,1.D0,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,z,0.D0)
      tncoff5 = -1/t1/t3/t6/t9*t12/512

CCC Numerator functions

      tcnum1 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ACnum1e0(var,ndim,z,ux,0.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ACnum1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ACnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ACnum1e1(var,ndim,z,ux,0.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      AC1e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

