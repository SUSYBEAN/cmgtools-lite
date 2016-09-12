CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq36e0int(var,ndim,Tcoeff)

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
      double precision sq3num3e0
      external sq3num3e0
      double precision sq3num3e1
      external sq3num3e1
      double precision sq3num3e2
      external sq3num3e2
      double precision sq3num3e3
      external sq3num3e3
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
      t3 = dlog(z)
      t4 = t3**2
      t7 = dlog(l4)
      t8 = t7**2
      t12 = dlog(1-z)
      t13 = t2**2
      t16 = t12**2
      t17 = dlog(2.D0)
      t21 = dlog(1-l4)
      t30 = t17**2
      t33 = 0.3141592653589793D1**2
      t36 = t21**2
      t47 = t2*t4/1024+t8*t3/2048-t12*t13/1024-t16*t17/64-t21*t8/4096-t8
     #*t17/1024-t21*t4/1024-t8*t12/1024-t21*t30/256+t21*t33/6144-t36*t7/
     #4096-t21*t13/4096+t36*t3/2048-t21*t16/256-t36*t12/1024
      t74 = t21*t2
      t77 = t21*t12
      t80 = t21*t7
      t83 = -t36*t17/1024-t7*t4/1024-t7*t16/256+t7*t33/6144-t7*t13/4096-
     #t7*t30/256+t36*t2/4096-t13*t17/1024-t8*t7/12288-t16*t12/192-t36*t2
     #1/12288+t13*t2/12288+t4*t3/1536-t74*t3/1024+t77*t2/512-t80*t12/512
      t87 = t12*t2
      t92 = t2*t3
      t118 = -t77*t17/128+t87*t17/128-t87*t3/256-t92*t17/256-z3/768-t30*
     #t17/192+t33*t17/1536-t2*t33/6144+t13*t3/2048-t4*t17/256+t2*t30/256
     #-t33*t3/3072+t3*t30/128+t12*t33/1536-t12*t4/256+t16*t3/128
      t121 = t7*t12
      t128 = t7*t3
      t131 = t7*t2
      t150 = t21*t3
      t153 = t12*t3
      t156 = -t12*t30/64+t121*t2/512-t121*t17/128+t121*t3/256+t128*t17/2
     #56+t131*t17/512-t131*t3/1024-t80*t17/512+t77*t3/256+t8*t2/4096+t16
     #*t2/256+t80*t2/2048+t74*t17/512+t80*t3/1024+t150*t17/256+t153*t17/
     #64
      t159 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t160 = 1/t159
      t162 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t163 = 1/t162
      t164 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t165 = 1/t164
      t167 = f1(0.D0,z,0.D0)
      t168 = t163*t165*t167
      t170 = 12*t4
      t171 = 3*t13
      t172 = 48*t30
      t173 = 48*t16
      t174 = 3*t8
      t175 = dlog(l3)
      t186 = 3*t36
      t187 = t175**2
      t189 = 24*t121
      t192 = t170+t171+t172+t173+t174+48*t175*t17-12*t175*t2+12*t21*t175
     #-24*t175*t3+12*t7*t175+t186+12*t187+t189+48*t175*t12
      t194 = 96*t12*t17
      t195 = 48*t153
      t196 = 12*t150
      t197 = 24*t77
      t199 = 24*t21*t17
      t201 = 24*t2*t17
      t202 = 12*t128
      t203 = 6*t131
      t205 = 24*t7*t17
      t206 = 6*t74
      t207 = 6*t80
      t208 = 2*t33
      t209 = 12*t92
      t210 = 24*t87
      t212 = 48*t3*t17
      t213 = t194-t195-t196+t197+t199-t201-t202-t203+t205-t206+t207-t208
     #+t209-t210-t212
      t215 = 1/l3
      t220 = t170+t171+t172+t173+t174+t186+t189+t194-t195-t196+t197+t199
     #-t201-t202
      t221 = dlog(l2)
      t234 = t221**2
      t236 = -t203+t205-t206+t207-6*t221*t2+24*t221*t17+6*t7*t221+24*t22
     #1*t12+6*t21*t221-12*t221*t3+3*t234-t208+t209-t210-t212
      t238 = 1/l2
      t239 = (t220+t236)*t238
      t240 = Dsquar322(z,0.D0,l2,0.D0,l4)
      t241 = 1/t240
      t247 = t165*t167
      t249 = 4*t17
      t250 = 2*t3
      t251 = 4*t12
      t252 = 2*t175
      t254 = (t7+t221-t2+t249-t250+t251+t21+t252)*t238
      t257 = t215*t160
      t263 = dlog(l1)
      t264 = t263**2
      t266 = t170+t171+t172+t173+t174+t186+12*t264+t189+t194-t195-t196+t
     #197+t199-t201
      t279 = -t202-t203+t205-t206+t207-t208+t209-24*t263*t3+12*t7*t263+1
     #2*t21*t263-12*t263*t2+48*t263*t12+48*t263*t17-t210-t212
      t281 = 1/l1
      t286 = 2*t263
      t295 = (t286+t251-t2+t249-t250+t221+t21+t7)*t281
      tncoff1 = (t47+t83+t118+t156)*t160*t168+(t192+t213)*t215*t160*t168
     #/6144+(-t239*t241/12288+t239*t160/12288)*t163*t247+(t254*t215*t241
     #/1024-t254*t257/1024)*t163*t247+(t266+t279)*t281*t160*t168/6144-(t
     #252+t286+t251-t2+t249-t250+t21+t7)*t281*t215*t160*t163*t247/512+(t
     #295*t238*t241/1024-t295*t238*t160/1024)*t163*t247+(-t241*t281*t238
     #*t215/512+t281*t238*t257/512)*t163*t247
      t1 = dlog(l3)
      t2 = t1**2
      t5 = dlog(1-z)
      t9 = dlog(1-l3)
      t10 = t9**2
      t12 = f2(l3,z)
      t13 = dlog(t12)
      t14 = t13**2
      t16 = dlog(z)
      t17 = t16**2
      t19 = t5**2
      t21 = dlog(l4)
      t22 = t21**2
      t25 = dlog(1-l4)
      t26 = t25**2
      t36 = dlog(2.D0)
      t49 = 12*t2+48*t1*t5+12*t10+3*t14+12*t17+48*t19+3*t22+3*t26+12*t13
     #*t16+12*t21*t9-12*t1*t13-6*t25*t13+96*t5*t36-12*t9*t13-48*t5*t16-2
     #4*t1*t16+12*t25*t1+48*t1*t36
      t84 = 0.3141592653589793D1**2
      t86 = t36**2
      t88 = -48*t16*t36-12*t21*t16+24*t21*t5+24*t21*t36+6*t25*t21+24*t25
     #*t36+24*t25*t5-12*t25*t16-24*t13*t36+12*t21*t1+48*t9*t36-24*t5*t13
     #+12*t25*t9+48*t9*t5-24*t9*t16-6*t21*t13+24*t1*t9-2*t84+48*t86
      t90 = 1/l3
      t92 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t93 = 1/t92
      t95 = s23r2(z,0.D0,0.D0,l3,l4)
      t96 = 1/t95
      t97 = s1vr2(z,0.D0,0.D0,l3,l4)
      t98 = 1/t97
      t100 = f1(l3,z,0.D0)
      t104 = 4*t5
      t105 = 4*t36
      t106 = 2*t16
      t107 = 2*t9
      t108 = dlog(l2)
      t109 = 2*t1
      t111 = 1/l2
      t112 = (-t21-t104-t105+t106-t107-t108+t13-t25-t109)*t111
      t113 = Dsquar322(z,0.D0,l2,l3,l4)
      t114 = 1/t113
      t117 = t90*t93
      t122 = t98*t100
      t124 = dlog(l1)
      t127 = 1/l1
      tncoff2 = -(t49+t88)*t90*t93*t96*t98*t100/6144+(t112*t90*t114/1024
     #-t112*t117/1024)*t96*t122-(-t109-2*t124+t13-t104+t106-t105-t107-t2
     #5-t21)*t127*t90*t93*t96*t122/512+(t127*t114*t111*t90/512-t127*t111
     #*t117/512)*t96*t122
      t1 = dlog(l1)
      t2 = dlog(2.D0)
      t6 = dlog(1-l1)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t11 = dlog(l4)
      t15 = dlog(1-l4)
      t21 = dlog(1-z)
      t24 = dlog(z)
      t29 = 0.3141592653589793D1**2
      t31 = t2**2
      t33 = f3(0.D0,z)
      t36 = dlog(1-l1*t33)
      t37 = t36**2
      t39 = -48*t1*t2+6*t6*t8-12*t11*t1-12*t15*t1+12*t1*t8-48*t21*t1-12*
     #t8*t24+24*t1*t24+2*t29-48*t31-3*t37
      t64 = -96*t21*t2-6*t15*t6-24*t21*t6-6*t11*t6+12*t36*t24+24*t21*t8-
     #24*t21*t36-6*t11*t36-6*t36*t6-12*t1*t6+6*t36*t8+48*t24*t2
      t88 = 12*t11*t24-24*t15*t21+24*t8*t2-24*t6*t2-24*t36*t2-6*t15*t11+
     #6*t11*t8+12*t6*t24+6*t15*t8-24*t11*t2-24*t11*t21
      t99 = t24**2
      t101 = t6**2
      t103 = t1**2
      t105 = t11**2
      t107 = t15**2
      t109 = t21**2
      t111 = t8**2
      t113 = -24*t15*t2-12*t1*t36-6*t15*t36+48*t21*t24+12*t15*t24-12*t99
     #-3*t101-12*t103-3*t105-3*t107-48*t109-3*t111
      t116 = 1/l1
      t118 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t119 = 1/t118
      t121 = s23r2(z,l1,0.D0,0.D0,l4)
      t122 = 1/t121
      t123 = s1vr2(z,l1,0.D0,0.D0,l4)
      t124 = 1/t123
      t126 = f1(0.D0,z,l1)
      t130 = dlog(l3)
      t132 = 2*t1
      t133 = 4*t21
      t134 = 4*t2
      t135 = 2*t24
      t138 = 1/l3
      t141 = t124*t126
      t142 = t119*t122*t141
      t145 = dlog(l2)
      t148 = 1/l2
      tncoff3 = (t39+t64+t88+t113)*t116*t119*t122*t124*t126/6144-(-2*t13
     #0-t132-t6-t133-t134+t8+t135-t36-t15-t11)*t116*t138*t142/512-(-t132
     #-t6-t133-t134+t8+t135-t36-t145-t15-t11)*t116*t148*t142/1024-t116*t
     #122*t119*t148*t138*t141/512
      t1 = dlog(l3)
      t4 = dlog(1-l1)
      t5 = f3(l3,z)
      t8 = dlog(1-l1*t5)
      t9 = f2(l3,z)
      t10 = dlog(t9)
      t12 = dlog(1-z)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t18 = dlog(l1)
      t21 = dlog(1-l3)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = -2*t1-t4-t8+t10-4*t12-4*t14+2*t16-2*t18-2*t21-t24-t25
      t27 = 1/l1
      t29 = 1/l3
      t31 = Dsquar322(z,l1,0.D0,l3,l4)
      t32 = 1/t31
      t33 = s23r2(z,l1,0.D0,l3,l4)
      t34 = 1/t33
      t36 = s1vr2(z,l1,0.D0,l3,l4)
      t38 = f1(l3,z,l1)
      t39 = 1/t36*t38
      tncoff4 = t26*t27*t29*t32*t34*t39/512+t27*t32*t34/l2*t29*t39/512
      t2 = dlog(1-z)
      t4 = dlog(l1)
      t7 = dlog(1-l1)
      t8 = dlog(z)
      t10 = l1*l2
      t12 = dlog(1-t10)
      t13 = dlog(2.D0)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t21 = dlog(l2)
      t23 = dlog(1-l4)
      t24 = dlog(l4)
      t25 = 4*t2+2*t4+t7-2*t8+t12+4*t13-t16+t20+t21+t23+t24
      t26 = 1/l1
      t28 = 1/l2
      t30 = Dsquar322(z,l1,l2,0.D0,l4)
      t31 = 1/t30
      t32 = s23r2(z,l1,t10,0.D0,l4)
      t33 = 1/t32
      t35 = s1vr2(z,l1,t10,0.D0,l4)
      t36 = 1/t35
      t37 = f1(0.D0,z,l1)
      tncoff5 = -t25*t26*t28*t31*t33*t36*t37/1024+t36*t33*t31*t26*t28/l3
     #*t37/512
      t1 = Dsquar322(z,l1,l2,l3,l4)
      t5 = l1*l2
      t6 = s23r2(z,l1,t5,l3,l4)
      t9 = s1vr2(z,l1,t5,l3,l4)
      t14 = f1(l3,z,l1)
      tncoff6 = -1/t1/l1/t6/t9/l2/l3*t14/512
      t1 = dlog(l4)
      t2 = dlog(z)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t9 = dlog(2.D0)
      t13 = dlog(1-z)
      t29 = dlog(1-l4)
      t32 = -t1*t2/1024-t1*t6/2048+t1*t9/512+t1*t13/512-t6*t9/512+t6*t2/
     #1024-t2*t9/256+t9*t13/128-t13*t2/256-t13*t6/512+t29*t9/512
      t41 = t13**2
      t43 = t1**2
      t45 = t29**2
      t47 = t9**2
      t49 = 0.3141592653589793D1**2
      t51 = t2**2
      t53 = t6**2
      t55 = -t29*t2/1024-t29*t6/2048+t29*t13/512+t29*t1/2048+t41/256+t43
     #/4096+t45/4096+t47/256-t49/6144+t51/1024+t53/4096
      t57 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t61 = 1/t60
      t62 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t63 = 1/t62
      t65 = f1(0.D0,z,0.D0)
      t66 = t61*t63*t65
      t68 = 12*t2
      t69 = 24*t13
      t70 = dlog(l3)
      t72 = 6*t29
      t73 = 6*t1
      t74 = 6*t6
      t75 = 24*t9
      t77 = 1/l3
      t82 = dlog(l2)
      t85 = 1/l2
      t86 = (t68-t69-6*t82-t72-t73+t74-t75)*t85
      t87 = Dsquar322(z,0.D0,l2,0.D0,l4)
      t88 = 1/t87
      t94 = t63*t65
      t96 = t85*t77
      t103 = dlog(l1)
      t106 = 1/l1
      t115 = t106*t85
      tncoff7 = (t32+t55)*t58*t66+(t68-t69-12*t70-t72-t73+t74-t75)*t77*t
     #58*t66/6144+(-t86*t88/12288+t86*t58/12288)*t61*t94+(-t96*t88/1024+
     #t96*t58/1024)*t61*t94+(t68-t69-12*t103-t72-t73+t74-t75)*t106*t58*t
     #66/6144+t106*t77*t58*t66/512+(-t115*t88/1024+t115*t58/1024)*t61*t9
     #4
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t7 = dlog(1-l3)
      t10 = dlog(1-l4)
      t12 = dlog(l4)
      t14 = f2(l3,z)
      t15 = dlog(t14)
      t17 = dlog(z)
      t19 = dlog(l3)
      t22 = 1/l3
      t24 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t25 = 1/t24
      t27 = s23r2(z,0.D0,0.D0,l3,l4)
      t28 = 1/t27
      t29 = s1vr2(z,0.D0,0.D0,l3,l4)
      t30 = 1/t29
      t32 = f1(l3,z,0.D0)
      t33 = t28*t30*t32
      t37 = 1/l2*t22
      t38 = Dsquar322(z,0.D0,l2,l3,l4)
      tncoff8 = -(-24*t2-24*t4-12*t7-6*t10-6*t12+6*t15+12*t17-12*t19)*t2
     #2*t25*t33/6144+(t37/t38/1024-t37*t25/1024)*t28*t30*t32-1/l1*t22*t2
     #5*t33/512
      t1 = dlog(l4)
      t4 = dlog(1-l4)
      t6 = dlog(2.D0)
      t8 = dlog(l1)
      t11 = dlog(1-z)
      t13 = dlog(z)
      t16 = dlog(1-l1)
      t18 = f3(0.D0,z)
      t21 = dlog(1-l1*t18)
      t23 = f2(0.D0,z)
      t24 = dlog(t23)
      t27 = 1/l1
      t29 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t30 = 1/t29
      t32 = s23r2(z,l1,0.D0,0.D0,l4)
      t34 = s1vr2(z,l1,0.D0,0.D0,l4)
      t37 = f1(0.D0,z,l1)
      t38 = 1/t32/t34*t37
      tncoff9 = (6*t1+6*t4+24*t6+12*t8+24*t11-12*t13+6*t16+6*t21-6*t24)*
     #t27*t30*t38/6144-t27/l3*t30*t38/512-t27/l2*t30*t38/1024
      t4 = Dsquar322(z,l1,0.D0,l3,l4)
      t7 = s23r2(z,l1,0.D0,l3,l4)
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t12 = f1(l3,z,l1)
      tncoff10 = 1/l1/l3/t4/t7/t9*t12/512
      t4 = Dsquar322(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t10 = s1vr2(z,l1,t7,0.D0,l4)
      t13 = f1(0.D0,z,l1)
      tncoff11 = 1/l1/l2/t4/t8/t10*t13/1024
      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t7 = dlog(1-z)
      t9 = dlog(l4)
      t12 = dlog(1-l4)
      t14 = dlog(z)
      t17 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = 1/t22
      t25 = f1(0.D0,z,0.D0)
      t26 = t21*t23*t25
      t32 = 1/l2
      t33 = Dsquar322(z,0.D0,l2,0.D0,l4)
      tncoff12 = (t2/2048-t4/512-t7/512-t9/2048-t12/2048+t14/1024)*t18*t
     #26+1/l3*t18*t26/1024+(-t32/t33/2048+t32*t18/2048)*t21*t23*t25+1/l1
     #*t18*t26/1024
      t2 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t10 = f1(l3,z,0.D0)
      tncoff13 = -1/l3/t2/t5/t7*t10/1024
      t2 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff14 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff15 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num3e0(var,ndim,z,ux,l1,l1*l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num3e0(var,ndim,z,ux,l1,l1*l2,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num3e1(var,ndim,z,ux,l1,l1*l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num3e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num3e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num3e3(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff15)

CCC Sum

      sq36e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq36e1int(var,ndim,Tcoeff)

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
      double precision sq3num3e0
      external sq3num3e0
      double precision sq3num3e1
      external sq3num3e1
      double precision sq3num3e2
      external sq3num3e2
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
      t2 = dlog(z)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t9 = dlog(2.D0)
      t13 = dlog(1-z)
      t29 = dlog(1-l4)
      t32 = -t1*t2/1024-t1*t6/2048+t1*t9/512+t1*t13/512-t6*t9/512+t6*t2/
     #1024-t2*t9/256+t13*t9/128-t13*t2/256-t13*t6/512+t29*t9/512
      t41 = t13**2
      t43 = t1**2
      t45 = t29**2
      t47 = t9**2
      t49 = 0.3141592653589793D1**2
      t51 = t2**2
      t53 = t6**2
      t55 = -t29*t2/1024-t29*t6/2048+t29*t13/512+t29*t1/2048+t41/256+t43
     #/4096+t45/4096+t47/256-t49/6144+t51/1024+t53/4096
      t57 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t58 = 1/t57
      t60 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t61 = 1/t60
      t62 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t63 = 1/t62
      t65 = f1(0.D0,z,0.D0)
      t66 = t61*t63*t65
      t68 = 4*t13
      t69 = 4*t9
      t70 = 2*t2
      t71 = dlog(l3)
      t74 = 1/l3
      t79 = dlog(l2)
      t81 = 1/l2
      t82 = (t29+t79+t68-t6+t69-t70+t1)*t81
      t83 = Dsquar322(z,0.D0,l2,0.D0,l4)
      t84 = 1/t83
      t90 = t63*t65
      t92 = t81*t74
      t99 = dlog(l1)
      t102 = 1/l1
      t111 = t102*t81
      tncoff1 = (t32+t55)*t58*t66-(t68-t6+t69-t70+2*t71+t29+t1)*t74*t58*
     #t66/1024+(t82*t84/2048-t82*t58/2048)*t61*t90+(-t92*t84/1024+t92*t5
     #8/1024)*t61*t90-(2*t99+t68-t6+t69-t70+t29+t1)*t102*t58*t66/1024+t1
     #02*t74*t58*t66/512+(-t111*t84/1024+t111*t58/1024)*t61*t90
      t2 = dlog(1-l3)
      t5 = dlog(1-z)
      t7 = f2(l3,z)
      t8 = dlog(t7)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      t13 = dlog(l3)
      t16 = dlog(1-l4)
      t17 = dlog(l4)
      t19 = 1/l3
      t21 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t22 = 1/t21
      t24 = s23r2(z,0.D0,0.D0,l3,l4)
      t25 = 1/t24
      t26 = s1vr2(z,0.D0,0.D0,l3,l4)
      t27 = 1/t26
      t29 = f1(l3,z,0.D0)
      t30 = t25*t27*t29
      t34 = 1/l2*t19
      t35 = Dsquar322(z,0.D0,l2,l3,l4)
      tncoff2 = -(-2*t2-4*t5+t8-4*t9+2*t11-2*t13-t16-t17)*t19*t22*t30/10
     #24+(t34/t35/1024-t34*t22/1024)*t25*t27*t29-1/l1*t19*t22*t30/512
      t1 = dlog(l1)
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t5 = f3(0.D0,z)
      t8 = dlog(1-l1*t5)
      t10 = dlog(1-z)
      t13 = dlog(1-l1)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t19 = dlog(1-l4)
      t20 = dlog(l4)
      t22 = 1/l1
      t24 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t25 = 1/t24
      t27 = s23r2(z,l1,0.D0,0.D0,l4)
      t29 = s1vr2(z,l1,0.D0,0.D0,l4)
      t32 = f1(0.D0,z,l1)
      t33 = 1/t27/t29*t32
      tncoff3 = -(-2*t1+t4-t8-4*t10-t13-4*t14+2*t16-t19-t20)*t22*t25*t33
     #/1024-t22/l3*t25*t33/512-t22/l2*t25*t33/1024
      t4 = Dsquar322(z,l1,0.D0,l3,l4)
      t7 = s23r2(z,l1,0.D0,l3,l4)
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t12 = f1(l3,z,l1)
      tncoff4 = 1/l1/l3/t4/t7/t9*t12/512
      t4 = Dsquar322(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t10 = s1vr2(z,l1,t7,0.D0,l4)
      t13 = f1(0.D0,z,l1)
      tncoff5 = 1/l1/l2/t4/t8/t10*t13/1024
      t1 = f2(0.D0,z)
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t7 = dlog(1-z)
      t9 = dlog(l4)
      t12 = dlog(1-l4)
      t14 = dlog(z)
      t17 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = 1/t22
      t25 = f1(0.D0,z,0.D0)
      t26 = t21*t23*t25
      t32 = 1/l2
      t33 = Dsquar322(z,0.D0,l2,0.D0,l4)
      tncoff6 = (t2/2048-t4/512-t7/512-t9/2048-t12/2048+t14/1024)*t18*t2
     #6+1/l3*t18*t26/1024+(-t32/t33/2048+t32*t18/2048)*t21*t23*t25+1/l1*
     #t18*t26/1024
      t2 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t10 = f1(l3,z,0.D0)
      tncoff7 = -1/l3/t2/t5/t7*t10/1024
      t2 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff8 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff9 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num3e0(var,ndim,z,ux,l1,l1*l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num3e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      sq36e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq36e2int(var,ndim,Tcoeff)

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
      double precision sq3num3e0
      external sq3num3e0
      double precision sq3num3e1
      external sq3num3e1
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
      t4 = dlog(2.D0)
      t7 = dlog(1-z)
      t9 = dlog(l4)
      t12 = dlog(1-l4)
      t14 = dlog(z)
      t17 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = 1/t22
      t25 = f1(0.D0,z,0.D0)
      t26 = t21*t23*t25
      t32 = 1/l2
      t33 = Dsquar322(z,0.D0,l2,0.D0,l4)
      tncoff1 = (t2/2048-t4/512-t7/512-t9/2048-t12/2048+t14/1024)*t18*t2
     #6+1/l3*t18*t26/1024+(-t32/t33/2048+t32*t18/2048)*t21*t23*t25+1/l1*
     #t18*t26/1024
      t2 = Dsquar322(z,0.D0,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t10 = f1(l3,z,0.D0)
      tncoff2 = -1/l3/t2/t5/t7*t10/1024
      t2 = Dsquar322(z,l1,0.D0,0.D0,l4)
      t5 = s23r2(z,l1,0.D0,0.D0,l4)
      t7 = s1vr2(z,l1,0.D0,0.D0,l4)
      t10 = f1(0.D0,z,l1)
      tncoff3 = -1/l1/t2/t5/t7*t10/1024
      t1 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff4 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq36e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq36e3int(var,ndim,Tcoeff)

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
      double precision sq3num3e0
      external sq3num3e0
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

      t1 = Dsquar322(z,0.D0,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff1 = 1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq36e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq36e4int(var,ndim,Tcoeff)

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

      sq36e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

