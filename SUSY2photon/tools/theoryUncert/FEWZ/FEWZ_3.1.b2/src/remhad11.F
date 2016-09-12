CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function rm11e0int(var,ndim,Tcoeff)

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
      t15 = f2(0.D0,1.D0)
      t16 = dlog(t15)
      t17 = t2*t16
      t24 = t2*t6
      t27 = t2*t12
      t30 = t2*t9
      t33 = t3*t16
      t40 = t3*t6
      t43 = t3*t12
      t46 = t9*t3
      t49 = t4*t6/4096+t4*t9/4096+t4*t12/2048-t17*t6/4096-t17*t12/2048-t
     #17*t9/4096+t24*t9/4096+t27*t6/2048+t30*t12/2048-t33*t6/4096-t33*t9
     #/4096-t33*t12/2048+t40*t9/4096+t43*t6/2048+t46*t12/2048
      t50 = t6*t9
      t53 = t16*t6
      t56 = t16*t12
      t59 = t16*t9
      t64 = t12**2
      t67 = t2**2
      t70 = t3**2
      t73 = t6**2
      t76 = t9**2
      t80 = t16**2
      t83 = 0.3141592653589793D1**2
      t92 = t50*t12/2048-t53*t9/4096-t56*t6/2048-t59*t12/2048-t4*t16/409
     #6+t64*t12/3072+t67*t2/24576+t70*t3/24576+t73*t6/24576+t76*t9/24576
     #+z3/1536-t80*t16/24576-t83*t12/6144+t12*t76/4096+t2*t70/8192+t2*t8
     #0/8192
      t126 = t2*t76/8192+t2*t64/2048-t2*t83/12288+t2*t73/8192+t67*t3/819
     #2-t67*t16/8192+t67*t6/8192+t67*t9/8192+t67*t12/4096+t3*t80/8192+t3
     #*t76/8192+t64*t3/2048-t3*t83/12288+t3*t73/8192+t80*t6/8192+t80*t9/
     #8192
      t159 = t80*t12/4096+t9*t64/2048-t9*t83/12288+t6*t76/8192+t6*t64/20
     #48-t6*t83/12288-t16*t73/8192-t16*t76/8192-t16*t64/2048+t16*t83/122
     #88+t73*t9/8192+t73*t12/4096-t70*t16/8192+t70*t6/8192+t70*t9/8192+t
     #70*t12/4096
      t162 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t163 = 1/t162
      t165 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t166 = t165**2
      t167 = 1/t166
      t168 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t169 = t168**2
      t170 = 1/t169
      t172 = f1(0.D0,1.D0,0.D0)
      t173 = t167*t170*t172
      t175 = 2*t12
      t176 = dlog(l3)
      t177 = 2*t176
      t178 = 1-z
      t179 = dlog(t178)
      t180 = 4*t179
      t183 = -1/t178
      t185 = 1/l3
      t189 = t163*t167*t170*t172
      t192 = dlog(l1)
      t194 = 1/l1
      t208 = 2*t83
      t209 = 3*t76
      t210 = 3*t70
      t211 = 3*t80
      t212 = 12*t64
      t213 = 3*t73
      t214 = 3*t67
      t215 = 12*t56
      t216 = 6*t4
      t217 = 6*t53
      t219 = 12*t9*t12
      t220 = 6*t59
      t221 = 6*t50
      t223 = 12*t12*t6
      t224 = -t208+t209+t210+t211+t212+t213+t214-t215+t216-t217+t219-t22
     #0+t221+t223
      t225 = 6*t33
      t226 = 6*t40
      t227 = 6*t46
      t228 = 12*t43
      t229 = 6*t17
      t230 = 6*t24
      t231 = 6*t30
      t232 = 12*t27
      t245 = t176**2
      t247 = -t225+t226+t227+t228-t229+t230+t231+t232+12*t176*t6-12*t176
     #*t16+12*t176*t9+24*t176*t12+12*t2*t176+12*t3*t176+12*t245
      t265 = t192**2
      t267 = -t225+t226+t227+t228-t229+t230+t231+t232+6*t192*t6-6*t192*t
     #16+6*t192*t9+12*t192*t12+6*t3*t192+6*t2*t192+3*t265
      t285 = 48*t179*t12+24*t2*t179+24*t3*t179+24*t179*t6-24*t179*t16+24
     #*t179*t9-t208+t209+t210+t211+t212+t213+t214-t215
      t286 = t179**2
      t288 = t216-t217+t219-t220+t221+t223-t225+t226+t227+t228-t229+t230
     #+t231+t232+48*t286
      tncoff1 = (t49+t92+t126+t159)*t163*t173-(-t16+t6+t175+t9+t177+t180
     #+t2+t3)*t183*t185*t189/512-(t180+t192+t6-t16+t9+t175+t3+t2)*t194*t
     #183*t189/1024+(t177-t16+t175+t6+t9+t192+t3+t2)*t185*t194*t189/2048
     #+t183*t185*t194*t189/512-(t224+t247)*t185*t163*t173/12288-(t224+t2
     #67)*t194*t163*t173/24576+(t285+t288)*t183*t163*t173/6144
      t2 = dlog(1-l2)
      t3 = t2**2
      t5 = dlog(l2)
      t6 = t5**2
      t8 = -l3+2
      t9 = dlog(t8)
      t10 = t9**2
      t12 = l3/2
      t13 = f2(t12,1.D0)
      t14 = dlog(t13)
      t15 = t14**2
      t17 = dlog(l4)
      t18 = t17**2
      t21 = dlog(1-l4)
      t22 = t21**2
      t24 = dlog(l3)
      t25 = t24**2
      t27 = 0.3141592653589793D1**2
      t41 = 3*t3+3*t6+12*t10+3*t15+3*t18+3*t22+12*t25-2*t27-12*t14*t9+12
     #*t5*t9+12*t2*t9-6*t14*t2-6*t14*t5-6*t21*t14
      t72 = 12*t21*t9-6*t17*t14+12*t17*t9-12*t24*t14+24*t24*t9+6*t21*t17
     #+6*t2*t5+6*t17*t2+6*t17*t5+6*t21*t2+6*t21*t5+12*t24*t2+12*t24*t5+1
     #2*t21*t24+12*t17*t24
      t75 = -1/t8
      t77 = 1/l3
      t79 = s34r2(1.D0,0.D0,l2,t12,l4)
      t80 = 1/t79
      t81 = s2vr2(1.D0,0.D0,l2,t12,l4)
      t82 = t81**2
      t83 = 1/t82
      t85 = s1vr2(1.D0,0.D0,l2,t12,l4)
      t86 = t85**2
      t87 = 1/t86
      t88 = f1(t12,1.D0,0.D0)
      t90 = t80*t83*t87*t88
      t93 = 2*t24
      t94 = 2*t9
      t95 = 1-z
      t96 = dlog(t95)
      t101 = -1/t95
      t110 = 1/l1*t77
      t115 = dlog(l1)
      tncoff2 = -(t41+t72)*t75*t77*t90/6144-(t2-t14+t93+t5+t94+4*t96+t21
     #+t17)*t75*t77*t101*t90/256+t75*t80*t87*t83*t110*t101*t88/256+(t93+
     #t2-t14+t115+t5+t94+t17+t21)*t75*t110*t90/1024
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l1)
      t5 = f2(0.D0,1.D0)
      t6 = dlog(t5)
      t7 = f3(0.D0,1.D0)
      t10 = dlog(1-l1*t7)
      t12 = dlog(1-l2)
      t13 = dlog(l2)
      t15 = dlog(1-l1)
      t16 = dlog(2.D0)
      t17 = 2*t16
      t18 = dlog(l4)
      t20 = dlog(1-l4)
      t22 = 1/l1
      t25 = -1/t1
      t27 = s34r2(1.D0,l1,l2,0.D0,l4)
      t28 = 1/t27
      t29 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t33 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t34 = t33**2
      t35 = 1/t34
      t36 = f1(0.D0,1.D0,l1)
      t38 = t28*t31*t35*t36
      t41 = dlog(l3)
      t44 = 1/l3
      t71 = -6*t6*t13+6*t12*t13+12*t16*t12-6*t18*t6+6*t18*t12+6*t18*t13+
     #12*t18*t16-6*t20*t6+6*t20*t12+6*t20*t13+12*t16*t20
      t96 = 6*t4*t12-6*t4*t6+6*t4*t13+12*t4*t16+6*t18*t4+6*t20*t4+6*t4*t
     #10+6*t4*t15+6*t18*t10+6*t18*t15+6*t20*t10+6*t20*t15
      t120 = 6*t13*t15+6*t10*t12+12*t15*t16+6*t12*t15-6*t10*t6+6*t10*t13
     #+6*t10*t15+12*t10*t16-6*t6*t15-12*t16*t6+6*t20*t18
      t125 = 0.3141592653589793D1**2
      t127 = t16**2
      t129 = t6**2
      t131 = t12**2
      t133 = t13**2
      t135 = t18**2
      t137 = t20**2
      t139 = t4**2
      t141 = t10**2
      t143 = t15**2
      t145 = -6*t12*t6+12*t13*t16-2*t125+12*t127+3*t129+3*t131+3*t133+3*
     #t135+3*t137+3*t139+3*t141+3*t143
      tncoff3 = (4*t2+t4-t6+t10+t12+t13+t15+t17+t18+t20)*t22*t25*t38/102
     #4-(2*t41+t4+t15+t12+t10+t13+t17-t6+t18+t20)*t44*t22*t38/2048+(t71+
     #t96+t120+t145)*t22*t28*t31*t35*t36/24576-t25*t44*t22*t38/512
      t1 = dlog(l3)
      t3 = dlog(l2)
      t5 = dlog(1-l1)
      t6 = -l3+2
      t7 = dlog(t6)
      t10 = dlog(1-l2)
      t11 = l3/2
      t12 = f2(t11,1.D0)
      t13 = dlog(t12)
      t14 = f3(t11,1.D0)
      t17 = dlog(1-l1*t14)
      t18 = dlog(l1)
      t19 = dlog(l4)
      t21 = dlog(1-l4)
      t24 = -1/t6
      t28 = 1/l1/l3
      t30 = s34r2(1.D0,l1,l2,t11,l4)
      t31 = 1/t30
      t32 = s2vr2(1.D0,l1,l2,t11,l4)
      t33 = t32**2
      t34 = 1/t33
      t36 = s1vr2(1.D0,l1,l2,t11,l4)
      t37 = t36**2
      t38 = 1/t37
      t39 = f1(t11,1.D0,l1)
      tncoff4 = -(2*t1+t3+t5+2*t7+t10-t13+t17+t18+t19+t21)*t24*t28*t31*t
     #34*t38*t39/1024-t24*t31*t38*t34*t28/(-1+z)*t39/256
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l3
      t5 = 1/l1
      t7 = s34r2(z,0.D0,l2,0.D0,l4)
      t8 = 1/t7
      t9 = s2vr2(z,0.D0,l2,0.D0,l4)
      t10 = t9**2
      t11 = 1/t10
      t13 = s1vr2(z,0.D0,l2,0.D0,l4)
      t14 = t13**2
      t15 = 1/t14
      t16 = f1(0.D0,z,0.D0)
      t18 = t8*t11*t15*t16
      t21 = dlog(2.D0)
      t22 = 2*t21
      t23 = dlog(l2)
      t24 = dlog(z)
      t25 = 2*t24
      t26 = dlog(l3)
      t29 = dlog(1-l2)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t33 = dlog(-t1)
      t34 = 4*t33
      t36 = dlog(1-l4)
      t37 = dlog(l4)
      t43 = dlog(l1)
      t85 = 6*t36*t37+12*t23*t21+6*t29*t23+12*t21*t29+6*t37*t29+6*t37*t2
     #3+12*t37*t21+6*t36*t29+6*t36*t23+12*t36*t21+48*t33*t21+24*t36*t33+
     #24*t37*t33+24*t33*t29+24*t33*t23-6*t36*t31-12*t36*t24-6*t37*t31
      t106 = t29**2
      t108 = t23**2
      t110 = t37**2
      t112 = t36**2
      t114 = 0.3141592653589793D1**2
      t116 = t21**2
      t118 = t33**2
      t120 = t31**2
      t122 = t24**2
      t124 = -12*t37*t24-24*t33*t31-48*t33*t24-6*t31*t23-12*t23*t24-24*t
     #21*t24-6*t31*t29-12*t31*t21+12*t31*t24-12*t29*t24+3*t106+3*t108+3*
     #t110+3*t112-2*t114+12*t116+48*t118+3*t120+12*t122
      tncoff5 = -t2*t3*t5*t18/512+(t22+t23-t25+2*t26+t29-t31+t34+t36+t37
     #)*t2*t3*t18/512+(t34+t43-t25+t29+t23-t31+t22+t37+t36)*t5*t2*t18/10
     #24-(t85+t124)*t2*t8*t11*t15*t16/6144
      t1 = dlog(l2)
      t2 = dlog(l3)
      t4 = l3/2
      t5 = f2(t4,z)
      t6 = dlog(t5)
      t7 = -l3+2
      t8 = dlog(t7)
      t10 = dlog(z)
      t13 = dlog(1-l2)
      t14 = 1-z
      t15 = dlog(t14)
      t18 = dlog(1-l4)
      t19 = dlog(l4)
      t22 = -1/t7
      t24 = 1/l3
      t26 = -1/t14
      t29 = s34r2(z,0.D0,l2,t4,l4)
      t30 = 1/t29
      t31 = s2vr2(z,0.D0,l2,t4,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,0.D0,l2,t4,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(t4,z,0.D0)
      tncoff6 = (t1+2*t2-t6+2*t8-2*t10+t13+4*t15+t18+t19)*t22*t24*t26*t3
     #0*t33*t37*t38/256-t22*t30*t37*t33/l1*t24*t26*t38/256
      t1 = 1-z
      t2 = dlog(t1)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t7 = dlog(1-l2)
      t9 = dlog(1-l1)
      t10 = dlog(l2)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t15 = dlog(l1)
      t16 = f3(0.D0,z)
      t19 = dlog(1-l1*t16)
      t20 = dlog(l4)
      t22 = dlog(1-l4)
      t23 = 4*t2-t5+t7+t9+t10+2*t11-2*t13+t15+t19+t20+t22
      t24 = 1/l1
      t27 = -1/t1
      t29 = s34r2(z,l1,l2,0.D0,l4)
      t31 = s2vr2(z,l1,l2,0.D0,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,l2,0.D0,l4)
      t36 = t35**2
      t38 = f1(0.D0,z,l1)
      t40 = 1/t29/t32/t36*t38
      tncoff7 = -t23*t24*t27*t40/1024+t27/l3*t24*t40/512
      t1 = l3/2
      t2 = f1(t1,z,l1)
      t6 = s34r2(z,l1,l2,t1,l4)
      t8 = s1vr2(z,l1,l2,t1,l4)
      t9 = t8**2
      t13 = s2vr2(z,l1,l2,t1,l4)
      t14 = t13**2
      tncoff8 = t2/(l3-2)/t6/t9/t14/l1/l3/(-1+z)/256
      t2 = dlog(1-l4)
      t3 = t2**2
      t5 = f2(0.D0,1.D0)
      t6 = dlog(t5)
      t7 = t6**2
      t9 = dlog(l2)
      t10 = t9**2
      t12 = dlog(2.D0)
      t13 = t12**2
      t15 = 0.3141592653589793D1**2
      t18 = dlog(1-l2)
      t19 = t18**2
      t21 = dlog(l4)
      t22 = t21**2
      t32 = -t3/8192-t7/8192-t10/8192-t13/2048+t15/12288-t19/8192-t22/81
     #92-t21*t2/4096+t2*t6/4096-t2*t18/4096-t2*t9/4096
      t55 = -t2*t12/2048+t21*t6/4096-t21*t18/4096-t21*t9/4096-t21*t12/20
     #48+t18*t6/4096+t6*t9/4096+t12*t6/2048-t18*t9/4096-t18*t12/2048-t9*
     #t12/2048
      t57 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t61 = t60**2
      t63 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t64 = t63**2
      t67 = f1(0.D0,1.D0,0.D0)
      t68 = 1/t61/t64*t67
      t70 = dlog(l3)
      t72 = 6*t21
      t73 = 6*t2
      t74 = 6*t6
      t75 = 6*t18
      t76 = 6*t9
      t77 = 12*t12
      t79 = 1/l3
      t84 = dlog(l1)
      t87 = 1/l1
      t92 = 1-z
      t93 = dlog(t92)
      t97 = -1/t92
      tncoff9 = (t32+t55)*t58*t68-(-12*t70-t72-t73+t74-t75-t76-t77)*t79*
     #t58*t68/12288-(-t72-t73+t74-t75-t76-t77-6*t84)*t87*t58*t68/24576+(
     #-t72-t73+t74-t75-t76-t77-24*t93)*t97*t58*t68/6144+t97*t79*t58*t68/
     #512+t87*t97*t58*t68/1024-t79*t87*t58*t68/2048
      t1 = l3-2
      t2 = 1/t1
      t5 = 1/l3
      t7 = l3/2
      t8 = s34r2(1.D0,0.D0,l2,t7,l4)
      t10 = s2vr2(1.D0,0.D0,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,0.D0,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,0.D0)
      t19 = 1/t8/t11/t15*t17
      t29 = dlog(1-l2)
      t31 = f2(t7,1.D0)
      t32 = dlog(t31)
      t35 = dlog(-t1)
      t37 = dlog(l2)
      t40 = dlog(1-l4)
      t42 = dlog(l4)
      t44 = dlog(l3)
      tncoff10 = -t2/l1*t5*t19/1024+t2*t5/(-1+z)*t19/256-(-6*t29+6*t32-1
     #2*t35-6*t37-6*t40-6*t42-12*t44)*t2*t5*t19/6144
      t2 = dlog(1-l2)
      t4 = f2(0.D0,1.D0)
      t5 = dlog(t4)
      t7 = dlog(l2)
      t10 = dlog(1-l1)
      t12 = dlog(2.D0)
      t14 = f3(0.D0,1.D0)
      t17 = dlog(1-l1*t14)
      t19 = dlog(l1)
      t21 = dlog(l4)
      t24 = dlog(1-l4)
      t27 = 1/l1
      t29 = s34r2(1.D0,l1,l2,0.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t33 = t32**2
      t35 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t36 = t35**2
      t39 = f1(0.D0,1.D0,l1)
      t40 = 1/t33/t36*t39
      tncoff11 = (-6*t2+6*t5-6*t7-6*t10-12*t12-6*t17-6*t19-6*t21-6*t24)*
     #t27*t30*t40/24576-t27/(-1+z)*t30*t40/1024+1/l3*t27*t30*t40/2048
      t7 = l3/2
      t8 = s34r2(1.D0,l1,l2,t7,l4)
      t10 = s2vr2(1.D0,l1,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,l1,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,l1)
      tncoff12 = 1/(l3-2)/l1/l3/t8/t11/t15*t17/1024
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t6 = 1-z
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(l2)
      t16 = dlog(2.D0)
      t19 = dlog(1-l2)
      t23 = -1/t6
      t25 = s34r2(z,0.D0,l2,0.D0,l4)
      t26 = 1/t25
      t28 = s2vr2(z,0.D0,l2,0.D0,l4)
      t29 = t28**2
      t31 = s1vr2(z,0.D0,l2,0.D0,l4)
      t32 = t31**2
      t35 = f1(0.D0,z,0.D0)
      t36 = 1/t29/t32*t35
      tncoff13 = -(-6*t2-6*t4-24*t7+12*t9+6*t12-6*t14-12*t16-6*t19)*t23*
     #t26*t36/6144-t23/l3*t26*t36/512-1/l1*t23*t26*t36/1024
      t8 = l3/2
      t9 = s34r2(z,0.D0,l2,t8,l4)
      t11 = s2vr2(z,0.D0,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(z,0.D0,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,z,0.D0)
      tncoff14 = -1/(l3-2)/l3/(-1+z)/t9/t12/t16*t18/256
      t5 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,l2,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff15 = 1/l1/(-1+z)/t5/t9/t12*t15/1024
      t2 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,0.D0)
      t13 = 1/t6/t9*t12
      t26 = dlog(1-l4)
      t28 = dlog(l4)
      t30 = f2(0.D0,1.D0)
      t31 = dlog(t30)
      t34 = dlog(1-l2)
      t36 = dlog(l2)
      t38 = dlog(2.D0)
      tncoff16 = -1/l3*t3*t13/2048-1/l1*t3*t13/4096+1/(-1+z)*t3*t13/1024
     #+(t26/4096+t28/4096-t31/4096+t34/4096+t36/4096+t38/2048)*t3*t13
      t5 = l3/2
      t6 = s34r2(1.D0,0.D0,l2,t5,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t5,l4)
      t10 = t9**2
      t12 = s1vr2(1.D0,0.D0,l2,t5,l4)
      t13 = t12**2
      t16 = f1(t5,1.D0,0.D0)
      tncoff17 = -1/(l3-2)/l3/t6/t10/t13*t16/1024
      t2 = s34r2(1.D0,l1,l2,0.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,l1)
      tncoff18 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,z,0.D0)
      tncoff19 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = f1(0.D0,1.D0,0.D0)
      tncoff20 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num1e0(var,ndim,z,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = rm1num1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = rm1num1e1(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = rm1num1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = rm1num1e2(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = rm1num1e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = rm1num1e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      rm11e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff1
     #6+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncof
     #f20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm11e1int(var,ndim,Tcoeff)

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

      t2 = dlog(1-l4)
      t3 = t2**2
      t5 = f2(0.D0,1.D0)
      t6 = dlog(t5)
      t7 = t6**2
      t9 = dlog(l2)
      t10 = t9**2
      t12 = dlog(2.D0)
      t13 = t12**2
      t15 = 0.3141592653589793D1**2
      t18 = dlog(1-l2)
      t19 = t18**2
      t21 = dlog(l4)
      t22 = t21**2
      t32 = -t3/8192-t7/8192-t10/8192-t13/2048+t15/12288-t19/8192-t22/81
     #92-t2*t21/4096+t2*t6/4096-t2*t18/4096-t2*t9/4096
      t55 = -t2*t12/2048+t21*t6/4096-t21*t18/4096-t21*t9/4096-t21*t12/20
     #48+t6*t18/4096+t6*t9/4096+t6*t12/2048-t18*t9/4096-t12*t18/2048-t9*
     #t12/2048
      t57 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t58 = 1/t57
      t60 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t61 = t60**2
      t63 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t64 = t63**2
      t67 = f1(0.D0,1.D0,0.D0)
      t68 = 1/t61/t64*t67
      t70 = 1-z
      t71 = dlog(t70)
      t73 = 2*t12
      t76 = -1/t70
      t81 = dlog(l3)
      t84 = 1/l3
      t93 = dlog(l1)
      t95 = 1/l1
      tncoff1 = (t32+t55)*t58*t68-(4*t71+t18-t6+t9+t73+t21+t2)*t76*t58*t
     #68/1024+(t18-t6+t9+t73+2*t81+t2+t21)*t84*t58*t68/2048+t76*t84*t58*
     #t68/512+(t2+t93+t18-t6+t9+t73+t21)*t95*t58*t68/4096+t95*t76*t58*t6
     #8/1024-t84*t95*t58*t68/2048
      t2 = dlog(1-l2)
      t3 = l3/2
      t4 = f2(t3,1.D0)
      t5 = dlog(t4)
      t6 = dlog(l2)
      t7 = -l3+2
      t8 = dlog(t7)
      t10 = dlog(l3)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t17 = -1/t7
      t19 = 1/l3
      t21 = s34r2(1.D0,0.D0,l2,t3,l4)
      t23 = s2vr2(1.D0,0.D0,l2,t3,l4)
      t24 = t23**2
      t27 = s1vr2(1.D0,0.D0,l2,t3,l4)
      t28 = t27**2
      t30 = f1(t3,1.D0,0.D0)
      t32 = 1/t21/t24/t28*t30
      tncoff2 = (t2-t5+t6+2*t8+2*t10+t13+t14)*t17*t19*t32/1024+t17*t19/(
     #-1+z)*t32/256-t17/l1*t19*t32/1024
      t2 = dlog(1-l4)
      t3 = dlog(l1)
      t4 = f2(0.D0,1.D0)
      t5 = dlog(t4)
      t6 = f3(0.D0,1.D0)
      t9 = dlog(1-l1*t6)
      t11 = dlog(1-l2)
      t12 = dlog(l2)
      t14 = dlog(1-l1)
      t15 = dlog(2.D0)
      t17 = dlog(l4)
      t19 = 1/l1
      t21 = s34r2(1.D0,l1,l2,0.D0,l4)
      t22 = 1/t21
      t24 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t25 = t24**2
      t27 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t28 = t27**2
      t31 = f1(0.D0,1.D0,l1)
      t32 = 1/t25/t28*t31
      tncoff3 = -(t2+t3-t5+t9+t11+t12+t14+2*t15+t17)*t19*t22*t32/4096-t1
     #9/(-1+z)*t22*t32/1024+1/l3*t19*t22*t32/2048
      t7 = l3/2
      t8 = s34r2(1.D0,l1,l2,t7,l4)
      t10 = s2vr2(1.D0,l1,l2,t7,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,l1,l2,t7,l4)
      t15 = t14**2
      t17 = f1(t7,1.D0,l1)
      tncoff4 = 1/(l3-2)/l1/l3/t8/t11/t15*t17/1024
      t1 = 1-z
      t2 = dlog(t1)
      t5 = dlog(1-l2)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t10 = dlog(z)
      t12 = dlog(l2)
      t13 = dlog(l4)
      t15 = dlog(1-l4)
      t18 = -1/t1
      t20 = s34r2(z,0.D0,l2,0.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,0.D0,l2,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,0.D0,l2,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,0.D0)
      t31 = 1/t24/t27*t30
      tncoff5 = (4*t2+t5-t7+2*t8-2*t10+t12+t13+t15)*t18*t21*t31/1024-t18
     #/l3*t21*t31/512-1/l1*t18*t21*t31/1024
      t8 = l3/2
      t9 = s34r2(z,0.D0,l2,t8,l4)
      t11 = s2vr2(z,0.D0,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(z,0.D0,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,z,0.D0)
      tncoff6 = -1/(l3-2)/l3/(-1+z)/t9/t12/t16*t18/256
      t5 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,l2,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff7 = 1/l1/(-1+z)/t5/t9/t12*t15/1024
      t2 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,0.D0)
      t13 = 1/t6/t9*t12
      t26 = dlog(1-l4)
      t28 = dlog(l4)
      t30 = f2(0.D0,1.D0)
      t31 = dlog(t30)
      t34 = dlog(1-l2)
      t36 = dlog(l2)
      t38 = dlog(2.D0)
      tncoff8 = -1/l3*t3*t13/2048-1/l1*t3*t13/4096+1/(-1+z)*t3*t13/1024+
     #(t26/4096+t28/4096-t31/4096+t34/4096+t36/4096+t38/2048)*t3*t13
      t5 = l3/2
      t6 = s34r2(1.D0,0.D0,l2,t5,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t5,l4)
      t10 = t9**2
      t12 = s1vr2(1.D0,0.D0,l2,t5,l4)
      t13 = t12**2
      t16 = f1(t5,1.D0,0.D0)
      tncoff9 = -1/(l3-2)/l3/t6/t10/t13*t16/1024
      t2 = s34r2(1.D0,l1,l2,0.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,l1)
      tncoff10 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,z,0.D0)
      tncoff11 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = f1(0.D0,1.D0,0.D0)
      tncoff12 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = rm1num1e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = rm1num1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = rm1num1e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      rm11e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm11e2int(var,ndim,Tcoeff)

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

      t2 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,0.D0)
      t13 = 1/t6/t9*t12
      t26 = dlog(1-l4)
      t28 = dlog(l4)
      t30 = f2(0.D0,1.D0)
      t31 = dlog(t30)
      t34 = dlog(1-l2)
      t36 = dlog(l2)
      t38 = dlog(2.D0)
      tncoff1 = -1/l3*t3*t13/2048-1/l1*t3*t13/4096+1/(-1+z)*t3*t13/1024+
     #(t26/4096+t28/4096-t31/4096+t34/4096+t36/4096+t38/2048)*t3*t13
      t5 = l3/2
      t6 = s34r2(1.D0,0.D0,l2,t5,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t5,l4)
      t10 = t9**2
      t12 = s1vr2(1.D0,0.D0,l2,t5,l4)
      t13 = t12**2
      t16 = f1(t5,1.D0,0.D0)
      tncoff2 = -1/(l3-2)/l3/t6/t10/t13*t16/1024
      t2 = s34r2(1.D0,l1,l2,0.D0,l4)
      t5 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t6 = t5**2
      t8 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = f1(0.D0,1.D0,l1)
      tncoff3 = 1/l1/t2/t6/t9*t12/4096
      t3 = s34r2(z,0.D0,l2,0.D0,l4)
      t6 = s2vr2(z,0.D0,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(z,0.D0,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,z,0.D0)
      tncoff4 = -1/(-1+z)/t3/t7/t10*t13/1024
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = f1(0.D0,1.D0,0.D0)
      tncoff5 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      rm11e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm11e3int(var,ndim,Tcoeff)

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

      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t4 = t3**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = f1(0.D0,1.D0,0.D0)
      tncoff1 = -1/t1/t4/t8*t10/4096

CCC Numerator functions

      tcnum1 = rm1num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      rm11e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm11e4int(var,ndim,Tcoeff)

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

      rm11e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

