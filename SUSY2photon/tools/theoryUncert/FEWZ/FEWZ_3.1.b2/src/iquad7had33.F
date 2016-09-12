CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qbiq73e0int(var,ndim,Tcoeff)

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
      double precision qbiq7num1e0
      external qbiq7num1e0
      double precision qbiq7num1e1
      external qbiq7num1e1
      double precision qbiq7num1e2
      external qbiq7num1e2
      double precision qbiq7num1e3
      external qbiq7num1e3
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
      t2 = f2(0.D0,1.D0)
      t3 = dlog(t2)
      t4 = f3(0.D0,1.D0)
      t6 = dlog(1-t4)
      t7 = dlog(2.D0)
      t8 = dlog(l1)
      t9 = dlog(l2)
      t10 = 2*t9
      t12 = dlog(1-l4)
      t13 = 1-z
      t14 = dlog(t13)
      t15 = 4*t14
      t18 = -1/t13
      t20 = 1/l2
      t22 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t23 = 1/t22
      t24 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t25 = 1/t24
      t27 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t28 = 1/t27
      t29 = f1(0.D0,1.D0,1.D0)
      t31 = t23*t25*t28*t29
      t34 = dlog(l3)
      t35 = 2*t34
      t38 = 1/l3
      t42 = t6**2
      t43 = 3*t42
      t44 = t34**2
      t46 = t8**2
      t47 = 3*t46
      t48 = t3**2
      t49 = 3*t48
      t50 = t1**2
      t51 = 3*t50
      t52 = t12**2
      t53 = 3*t52
      t54 = t12*t1
      t55 = 6*t54
      t56 = t6*t8
      t57 = 6*t56
      t59 = 6*t7*t6
      t60 = t6*t3
      t61 = 6*t60
      t63 = 6*t3*t8
      t64 = t7*t3
      t65 = 6*t64
      t69 = 6*t1*t8
      t70 = t43+12*t44+t47+t49+t51+t53+t55+t57+t59-t61-t63-t65+12*t34*t8
     #+t69
      t71 = t1*t7
      t72 = 6*t71
      t73 = t1*t3
      t74 = 6*t73
      t75 = t12*t6
      t76 = 6*t75
      t78 = 6*t12*t8
      t79 = t12*t7
      t80 = 6*t79
      t81 = t12*t3
      t82 = 6*t81
      t83 = t1*t6
      t84 = 6*t83
      t86 = 6*t7*t8
      t97 = t7**2
      t98 = 3*t97
      t99 = 0.3141592653589793D1**2
      t100 = 2*t99
      t101 = t72-t74+t76+t78+t80-t82+t84+t86+12*t34*t7+12*t12*t34-12*t34
     #*t3+12*t1*t34+12*t34*t6+t98-t100
      t106 = t25*t28*t29
      t109 = t14**2
      t111 = 48*t109+t43+t47+t49+t51+t53+t55+t57+t59-t61-t63-t65+t69+t72
      t124 = -t74+t76+t78+t80-t82+t84+t86+24*t14*t6+24*t12*t14-24*t14*t3
     #+24*t14*t8+24*t14*t7+24*t1*t14+t98-t100
      t130 = t43+t47+t49+t51+t53+t55+t57+t59-t61-t63-t65+t69+t72-t74
      t143 = t9**2
      t145 = t76+t78+t80-t82+t84+t86-12*t9*t3+12*t9*t8+12*t1*t9+12*t7*t9
     #+12*t9*t6+12*t12*t9+12*t143+t98-t100
      t181 = t42*t7/32768-t42*t3/32768+t12*t48/32768-t3*t46/32768+t1*t48
     #/32768+t1*t42/32768+t7*t46/32768+t8*t97/32768-t8*t99/49152+t42*t8/
     #32768+t6*t46/32768-t99*t7/49152+t52*t1/32768-t52*t3/32768+t52*t6/3
     #2768
      t214 = t52*t8/32768+t52*t7/32768-t50*t3/32768+t1*t46/32768+t1*t97/
     #32768-t1*t99/49152+t12*t50/32768+t50*t6/32768+t50*t8/32768+t12*t42
     #/32768+t50*t7/32768+t12*t46/32768+t12*t97/32768-t12*t99/49152+t48*
     #t7/32768-t3*t97/32768
      t247 = t48*t8/32768+t6*t97/32768-t6*t99/49152+t6*t48/32768+t3*t99/
     #49152+t46*t8/98304+t42*t6/98304+t50*t1/98304+t52*t12/98304+z3/6144
     #+t97*t7/98304-t48*t3/98304-t60*t7/16384+t56*t7/16384+t54*t7/16384+
     #t54*t8/16384
      t264 = -t81*t8-t75*t3-t64*t8-t60*t8-t83*t3-t73*t8-t71*t3+t71*t8+t8
     #3*t8+t71*t6-t54*t3+t54*t6+t79*t8-t79*t3+t75*t8+t79*t6
      tncoff1 = -(t1-t3+t6+t7+t8+t10+t12+t15)*t18*t20*t31/2048+(t35-t3+t
     #10+t6+t8+t7+t12+t1)*t20*t38*t31/4096-(t70+t101)*t38*t23*t106/49152
     #+(t111+t124)*t18*t23*t106/24576-(t130+t145)*t20*t23*t106/49152+(t1
     #81+t214+t247+t264/16384)*t23*t106-(t35+t7+t6+t8-t3+t12+t1+t15)*t18
     #*t38*t31/2048+t18*t20*t38*t31/1024
      t1 = 0.3141592653589793D1**2
      t3 = dlog(2.D0)
      t4 = t3**2
      t7 = dlog(1-l4)
      t8 = l3/2
      t9 = f2(t8,1.D0)
      t10 = dlog(t9)
      t13 = f3(t8,1.D0)
      t15 = dlog(1-t13)
      t18 = dlog(l4)
      t21 = 2-l3
      t22 = dlog(t21)
      t29 = dlog(l3)
      t38 = dlog(l1)
      t51 = -2*t1+3*t4-6*t7*t10-6*t15*t3-6*t18*t10+12*t7*t22+6*t7*t15+6*
     #t7*t18-12*t29*t10+12*t18*t22+6*t18*t15-6*t10*t15+12*t38*t22-12*t3*
     #t22-12*t10*t22+12*t15*t22+6*t10*t3-6*t10*t38
      t76 = t22**2
      t78 = t10**2
      t80 = t29**2
      t82 = t38**2
      t84 = t15**2
      t86 = t18**2
      t88 = t7**2
      t90 = 6*t15*t38+24*t29*t22+12*t29*t15+12*t29*t38+6*t18*t38-6*t18*t
     #3+6*t7*t38-6*t7*t3-6*t3*t38-12*t29*t3+12*t7*t29+12*t18*t29+12*t76+
     #3*t78+12*t80+3*t82+3*t84+3*t86+3*t88
      t93 = -1/t21
      t95 = 1/l3
      t97 = Diquad721(1.D0,l1,0.D0,t8,l4)
      t99 = s34r2(1.D0,1.D0,1.D0,t8,l4)
      t102 = s23r2(1.D0,1.D0,1.D0,t8,l4)
      t104 = f1(t8,1.D0,1.D0)
      t106 = 1/t97/t99/t102*t104
      t109 = -1+z
      t111 = 1/t109*t95
      t113 = 1/l2*t93
      t117 = 2*t22
      t118 = 2*t29
      t120 = dlog(-t109)
      t127 = dlog(l2)
      tncoff2 = -(t51+t90)*t93*t95*t106/24576+t111*t113*t106/512-(-t10+t
     #117+t38+t118+t15-t3+t7+t18+4*t120)*t93*t111*t106/1024+(t15-t3+t38+
     #t117-t10+2*t127+t118+t7+t18)*t95*t113*t106/2048
      t1 = dlog(2.D0)
      t2 = t1**2
      t4 = l1*l2
      t5 = 2-t4
      t6 = dlog(t5)
      t7 = t6**2
      t9 = dlog(l1)
      t10 = t9**2
      t12 = f3(0.D0,1.D0)
      t17 = dlog(2-2*t12+t12*l1*l2)
      t18 = t17**2
      t20 = f2(0.D0,1.D0)
      t21 = dlog(t20)
      t22 = t21**2
      t24 = dlog(l4)
      t25 = t24**2
      t27 = 1-l2
      t28 = dlog(t27)
      t29 = t28**2
      t31 = 0.3141592653589793D1**2
      t36 = dlog(1-l4)
      t41 = 3*t2+3*t7+3*t10+3*t18+3*t22+3*t25+3*t29-2*t31-6*t28*t21+6*t3
     #6*t6+6*t36*t28
      t52 = dlog(l2)
      t67 = 6*t36*t17+6*t28*t6-6*t21*t6+6*t17*t6+6*t9*t6+12*t52*t17+6*t2
     #4*t17+12*t52*t6+12*t52*t28+6*t24*t6+6*t24*t28-12*t52*t21
      t91 = 12*t52*t9+12*t24*t52-12*t52*t1+6*t36*t24+12*t36*t52-6*t21*t9
     #+6*t1*t21+6*t24*t9-6*t24*t1-6*t24*t21+6*t36*t9
      t112 = t36**2
      t114 = t52**2
      t116 = -6*t36*t1-6*t36*t21-6*t1*t9-6*t28*t1-6*t17*t1+6*t9*t17+6*t2
     #8*t17-6*t21*t17-6*t1*t6+6*t28*t9+3*t112+12*t114
      t120 = -1/t5
      t122 = 1/l2
      t124 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t127 = 1-t4/2
      t128 = s34r2(1.D0,t127,t27,0.D0,l4)
      t131 = s23r2(1.D0,t127,t27,0.D0,l4)
      t133 = f1(0.D0,1.D0,t127)
      t135 = 1/t124/t128/t131*t133
      t138 = 2*t52
      t139 = dlog(l3)
      t142 = 1/l3
      t144 = t120*t122
      t148 = -1+z
      t149 = 1/t148
      t155 = dlog(-t148)
      tncoff3 = -(t41+t67+t91+t116)*t120*t122*t135/24576+(t28-t21+t17+t9
     #+t138-t1+t6+2*t139+t36+t24)*t142*t144*t135/2048+t149*t142*t144*t13
     #5/512-(t24-t1+t6+t28-t21+t17+t9+t138+t36+4*t155)*t120*t122*t149*t1
     #35/1024
      t1 = l3/2
      t2 = f3(t1,1.D0)
      t7 = dlog(2-2*t2+t2*l1*l2)
      t8 = 1-l2
      t9 = dlog(t8)
      t10 = dlog(l1)
      t11 = l1*l2
      t12 = 2-t11
      t13 = dlog(t12)
      t14 = dlog(2.D0)
      t16 = dlog(l2)
      t18 = 2-l3
      t19 = dlog(t18)
      t21 = f2(t1,1.D0)
      t22 = dlog(t21)
      t23 = dlog(l3)
      t26 = dlog(1-l4)
      t27 = dlog(l4)
      t28 = t7+t9+t10+t13-3*t14+2*t16+2*t19-t22+2*t23+t26+t27
      t30 = -1/t18
      t33 = -1/t12
      t34 = 1/l2
      t37 = 1/l3
      t38 = Diquad721(1.D0,l1,l2,t1,l4)
      t39 = 1/t38
      t42 = 1-t11/2
      t43 = s34r2(1.D0,t42,t8,t1,l4)
      t45 = s23r2(1.D0,t42,t8,t1,l4)
      t47 = 1/t43/t45
      t48 = f1(t1,1.D0,t42)
      tncoff4 = t28*t30*t33*t34*t37*t39*t47*t48/1024+t33*t30*t34*t37*t47
     #*t39/(-1+z)*t48/256
      t1 = dlog(l4)
      t2 = dlog(2.D0)
      t3 = dlog(l2)
      t5 = dlog(z)
      t6 = 2*t5
      t7 = f3(0.D0,z)
      t9 = dlog(1-t7)
      t10 = dlog(l1)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(1-l4)
      t15 = 1-z
      t16 = dlog(t15)
      t17 = 4*t16
      t20 = -1/t15
      t22 = 1/l2
      t24 = Diquad721(z,l1,0.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t27 = 1/t26
      t29 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t30 = 1/t29
      t31 = f1(0.D0,z,1.D0)
      t33 = t25*t27*t30*t31
      t36 = dlog(l3)
      t40 = 1/l3
      t44 = t9**2
      t46 = t16**2
      t48 = t10**2
      t50 = t12**2
      t52 = t5**2
      t54 = t1**2
      t56 = t14**2
      t58 = 0.3141592653589793D1**2
      t60 = t2**2
      t80 = 3*t44+48*t46+3*t48+3*t50+12*t52+3*t54+3*t56-2*t58+3*t60+24*t
     #14*t16+6*t1*t14+24*t16*t10+24*t16*t2+24*t1*t16+6*t1*t10+6*t1*t2+6*
     #t14*t10+6*t14*t2
      t119 = 6*t2*t10-12*t2*t5-6*t12*t14+6*t14*t9-12*t9*t5-12*t14*t5-6*t
     #12*t10-6*t12*t2+12*t12*t5-6*t12*t1-12*t5*t10-6*t12*t9+6*t9*t10+6*t
     #9*t2-12*t1*t5-24*t16*t12+24*t16*t9-48*t16*t5+6*t1*t9
      tncoff5 = (t1+t2+2*t3-t6+t9+t10-t12+t14+t17)*t20*t22*t33/2048+(2*t
     #36-t12-t6+t9+t10+t2+t14+t1+t17)*t20*t40*t33/2048-(t80+t119)*t20*t2
     #5*t27*t30*t31/24576-t20*t22*t40*t33/1024
      t1 = 2-l3
      t2 = dlog(t1)
      t4 = dlog(l3)
      t6 = l3/2
      t7 = f3(t6,z)
      t9 = dlog(1-t7)
      t10 = dlog(z)
      t12 = f2(t6,z)
      t13 = dlog(t12)
      t14 = dlog(l1)
      t15 = dlog(2.D0)
      t17 = dlog(1-l4)
      t18 = dlog(l4)
      t19 = 1-z
      t20 = dlog(t19)
      t24 = -1/t1
      t29 = -1/l3/t19
      t31 = Diquad721(z,l1,0.D0,t6,l4)
      t33 = s34r2(z,1.D0,1.D0,t6,l4)
      t36 = s23r2(z,1.D0,1.D0,t6,l4)
      t38 = f1(t6,z,1.D0)
      t40 = 1/t31/t33/t36*t38
      tncoff6 = (2*t2+2*t4+t9-2*t10-t13+t14-t15+t17+t18+4*t20)*t24*t29*t
     #40/1024-t29/l2*t24*t40/512
      t1 = -1+z
      t2 = 1/t1
      t5 = l1*l2
      t6 = -2+t5
      t7 = 1/t6
      t8 = 1/l2
      t11 = Diquad721(z,l1,l2,0.D0,l4)
      t14 = 1-t5/2
      t15 = 1-l2
      t16 = s34r2(z,t14,t15,0.D0,l4)
      t19 = s23r2(z,t14,t15,0.D0,l4)
      t21 = f1(0.D0,z,t14)
      t23 = 1/t11/t16/t19*t21
      t26 = dlog(l4)
      t27 = f2(0.D0,z)
      t28 = dlog(t27)
      t29 = dlog(2.D0)
      t30 = dlog(z)
      t32 = dlog(l2)
      t34 = f3(0.D0,z)
      t39 = dlog(2-2*t34+t34*l1*l2)
      t40 = dlog(t15)
      t41 = dlog(l1)
      t43 = dlog(-t6)
      t45 = dlog(1-l4)
      t47 = dlog(-t1)
      t49 = t26-t28-t29-2*t30+2*t32+t39+t40+t41+t43+t45+4*t47
      tncoff7 = -t2/l3*t7*t8*t23/512+t49*t7*t8*t2*t23/1024
      t1 = l3/2
      t2 = l1*l2
      t4 = 1-t2/2
      t5 = f1(t1,z,t4)
      t6 = 1-l2
      t7 = s23r2(z,t4,t6,t1,l4)
      t16 = s34r2(z,t4,t6,t1,l4)
      t18 = Diquad721(z,l1,l2,t1,l4)
      tncoff8 = -t5/t7/(-2+t2)/(-2+l3)/t16/t18/l2/l3/(-1+z)/256
      t1 = dlog(l2)
      t3 = dlog(l4)
      t4 = 6*t3
      t5 = f2(0.D0,1.D0)
      t6 = dlog(t5)
      t7 = 6*t6
      t8 = f3(0.D0,1.D0)
      t10 = dlog(1-t8)
      t11 = 6*t10
      t12 = dlog(l1)
      t13 = 6*t12
      t14 = dlog(2.D0)
      t15 = 6*t14
      t17 = dlog(1-l4)
      t18 = 6*t17
      t20 = 1/l2
      t22 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t23 = 1/t22
      t25 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t27 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t30 = f1(0.D0,1.D0,1.D0)
      t31 = 1/t25/t27*t30
      t34 = -1+z
      t35 = 1/t34
      t40 = 1/l3
      t49 = dlog(l3)
      t67 = -t17*t3+t17*t6-t17*t10-t17*t12-t17*t14+t6*t3-t10*t3-t12*t3-t
     #14*t3-t12*t14-t12*t10
      t71 = t12**2
      t73 = t17**2
      t75 = t3**2
      t79 = t14**2
      t81 = 0.3141592653589793D1**2
      t83 = t6**2
      t85 = t10**2
      t91 = t12*t6/16384-t71/32768-t73/32768-t75/32768+t14*t6/16384-t79/
     #32768+t81/49152-t83/32768-t85/32768-t14*t10/16384+t10*t6/16384
      t96 = dlog(-t34)
      tncoff9 = -(-12*t1-t4+t7-t11-t13-t15-t18)*t20*t23*t31/49152+t35*t2
     #0*t23*t31/2048-t20*t40*t23*t31/4096+t35*t40*t23*t31/2048-(-12*t49-
     #t4+t7-t11-t13-t15-t18)*t40*t23*t31/49152+(t67/16384+t91)*t23*t31+(
     #-t4-24*t96+t7-t11-t13-t15-t18)*t35*t23*t31/24576
      t1 = 2-l3
      t2 = dlog(t1)
      t5 = dlog(1-l4)
      t7 = dlog(l4)
      t9 = dlog(l3)
      t11 = dlog(2.D0)
      t13 = l3/2
      t14 = f2(t13,1.D0)
      t15 = dlog(t14)
      t17 = f3(t13,1.D0)
      t19 = dlog(1-t17)
      t21 = dlog(l1)
      t25 = -1/t1
      t27 = 1/l3
      t29 = Diquad721(1.D0,l1,0.D0,t13,l4)
      t31 = s34r2(1.D0,1.D0,1.D0,t13,l4)
      t34 = s23r2(1.D0,1.D0,1.D0,t13,l4)
      t36 = f1(t13,1.D0,1.D0)
      t38 = 1/t29/t31/t34*t36
      tncoff10 = -(-12*t2-6*t5-6*t7-12*t9+6*t11+6*t15-6*t19-6*t21)*t25*t
     #27*t38/24576+t25*t27/(-1+z)*t38/1024-t27/l2*t25*t38/2048
      t1 = dlog(l2)
      t4 = dlog(1-l4)
      t6 = dlog(l1)
      t8 = 1-l2
      t9 = dlog(t8)
      t11 = dlog(2.D0)
      t13 = l1*l2
      t14 = 2-t13
      t15 = dlog(t14)
      t17 = f2(0.D0,1.D0)
      t18 = dlog(t17)
      t20 = dlog(l4)
      t22 = f3(0.D0,1.D0)
      t27 = dlog(2-2*t22+t22*l1*l2)
      t31 = -1/t14
      t33 = 1/l2
      t35 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t38 = 1-t13/2
      t39 = s34r2(1.D0,t38,t8,0.D0,l4)
      t42 = s23r2(1.D0,t38,t8,0.D0,l4)
      t44 = f1(0.D0,1.D0,t38)
      t46 = 1/t35/t39/t42*t44
      tncoff11 = -(-12*t1-6*t4-6*t6-6*t9+6*t11-6*t15+6*t18-6*t20-6*t27)*
     #t31*t33*t46/24576-1/l3*t31*t33*t46/2048+t31*t33/(-1+z)*t46/1024
      t3 = l1*l2
      t11 = l3/2
      t12 = Diquad721(1.D0,l1,l2,t11,l4)
      t15 = 1-t3/2
      t16 = 1-l2
      t17 = s34r2(1.D0,t15,t16,t11,l4)
      t20 = s23r2(1.D0,t15,t16,t11,l4)
      t22 = f1(t11,1.D0,t15)
      tncoff12 = -1/(-2+l3)/(-2+t3)/l2/l3/t12/t17/t20*t22/1024
      t1 = -1+z
      t2 = 1/t1
      t5 = Diquad721(z,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      t14 = 1/t8/t10*t13
      t22 = f3(0.D0,z)
      t24 = dlog(1-t22)
      t26 = f2(0.D0,z)
      t27 = dlog(t26)
      t30 = dlog(-t1)
      t32 = dlog(l4)
      t35 = dlog(1-l4)
      t37 = dlog(l1)
      t39 = dlog(2.D0)
      t41 = dlog(z)
      tncoff13 = -t2/l3*t6*t14/2048-t2/l2*t6*t14/2048-(-6*t24+6*t27-24*t
     #30-6*t32-6*t35-6*t37-6*t39+12*t41)*t2*t6*t14/24576
      t8 = l3/2
      t9 = Diquad721(z,l1,0.D0,t8,l4)
      t11 = s34r2(z,1.D0,1.D0,t8,l4)
      t14 = s23r2(z,1.D0,1.D0,t8,l4)
      t16 = f1(t8,z,1.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t9/t11/t14*t16/1024
      t1 = l1*l2
      t9 = Diquad721(z,l1,l2,0.D0,l4)
      t12 = 1-t1/2
      t13 = 1-l2
      t14 = s34r2(z,t12,t13,0.D0,l4)
      t17 = s23r2(z,t12,t13,0.D0,l4)
      t19 = f1(0.D0,z,t12)
      tncoff15 = -1/(-2+t1)/l2/(-1+z)/t9/t14/t17*t19/1024
      t2 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,1.D0,1.D0)
      t11 = 1/t5/t7*t10
      t23 = dlog(l4)
      t25 = dlog(1-l4)
      t26 = f2(0.D0,1.D0)
      t27 = dlog(t26)
      t28 = dlog(2.D0)
      t29 = dlog(l1)
      t30 = f3(0.D0,1.D0)
      t32 = dlog(1-t30)
      tncoff16 = -1/l2*t3*t11/8192-1/l3*t3*t11/8192+1/(-1+z)*t3*t11/4096
     #+(t23/16384+t25/16384-t27/16384+t28/16384+t29/16384+t32/16384)*t3*
     #t11
      t5 = l3/2
      t6 = Diquad721(1.D0,l1,0.D0,t5,l4)
      t9 = s34r2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = f1(t5,1.D0,1.D0)
      tncoff17 = -1/(-2+l3)/l3/t6/t9/t11*t14/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,0.D0,l4)
      t14 = s23r2(1.D0,t10,t11,0.D0,l4)
      t17 = f1(0.D0,1.D0,t10)
      tncoff18 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,0.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,z,1.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,1.D0,1.D0)
      tncoff20 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qbiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qbiq7num1e2(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qbiq7num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qbiq7num1e3(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qbiq73e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq73e1int(var,ndim,Tcoeff)

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
      double precision qbiq7num1e0
      external qbiq7num1e0
      double precision qbiq7num1e1
      external qbiq7num1e1
      double precision qbiq7num1e2
      external qbiq7num1e2
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
      t5 = f2(0.D0,1.D0)
      t6 = dlog(t5)
      t8 = f3(0.D0,1.D0)
      t10 = dlog(1-t8)
      t12 = dlog(l1)
      t14 = dlog(2.D0)
      t22 = -t2*t3+t2*t6-t2*t10-t12*t2-t14*t2+t6*t3-t10*t3-t12*t3-t14*t3
     #-t12*t14-t12*t10
      t26 = t12**2
      t28 = t2**2
      t30 = t3**2
      t34 = t14**2
      t36 = 0.3141592653589793D1**2
      t38 = t6**2
      t40 = t10**2
      t46 = t12*t6/16384-t26/32768-t28/32768-t30/32768+t14*t6/16384-t34/
     #32768+t36/49152-t38/32768-t40/32768-t14*t10/16384+t10*t6/16384
      t48 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t49 = 1/t48
      t51 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t53 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t56 = f1(0.D0,1.D0,1.D0)
      t57 = 1/t51/t53*t56
      t59 = 1-z
      t60 = dlog(t59)
      t64 = -1/t59
      t69 = dlog(l3)
      t72 = 1/l3
      t81 = dlog(l2)
      t84 = 1/l2
      tncoff1 = (t22/16384+t46)*t49*t57-(4*t60+t3-t6+t10+t12+t14+t2)*t64
     #*t49*t57/4096+(2*t69-t6+t10+t12+t14+t2+t3)*t72*t49*t57/8192+t64*t7
     #2*t49*t57/2048+(t3+t14-t6+t12+t10+2*t81+t2)*t84*t49*t57/8192+t64*t
     #84*t49*t57/2048-t84*t72*t49*t57/4096
      t1 = dlog(l3)
      t3 = dlog(2.D0)
      t4 = 2-l3
      t5 = dlog(t4)
      t7 = l3/2
      t8 = f3(t7,1.D0)
      t10 = dlog(1-t8)
      t11 = dlog(l1)
      t12 = f2(t7,1.D0)
      t13 = dlog(t12)
      t15 = dlog(1-l4)
      t16 = dlog(l4)
      t19 = -1/t4
      t21 = 1/l3
      t23 = Diquad721(1.D0,l1,0.D0,t7,l4)
      t25 = s34r2(1.D0,1.D0,1.D0,t7,l4)
      t28 = s23r2(1.D0,1.D0,1.D0,t7,l4)
      t30 = f1(t7,1.D0,1.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff2 = (2*t1-t3+2*t5+t10+t11-t13+t15+t16)*t19*t21*t32/4096+t19*
     #t21/(-1+z)*t32/1024-t21/l2*t19*t32/2048
      t1 = dlog(l4)
      t2 = dlog(2.D0)
      t3 = l1*l2
      t4 = 2-t3
      t5 = dlog(t4)
      t6 = f3(0.D0,1.D0)
      t11 = dlog(2-2*t6+t6*l1*l2)
      t12 = dlog(l1)
      t13 = 1-l2
      t14 = dlog(t13)
      t15 = f2(0.D0,1.D0)
      t16 = dlog(t15)
      t17 = dlog(l2)
      t20 = dlog(1-l4)
      t23 = -1/t4
      t25 = 1/l2
      t27 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t30 = 1-t3/2
      t31 = s34r2(1.D0,t30,t13,0.D0,l4)
      t34 = s23r2(1.D0,t30,t13,0.D0,l4)
      t36 = f1(0.D0,1.D0,t30)
      t38 = 1/t27/t31/t34*t36
      tncoff3 = (t1-t2+t5+t11+t12+t14-t16+2*t17+t20)*t23*t25*t38/4096+t2
     #3*t25/(-1+z)*t38/1024-1/l3*t23*t25*t38/2048
      t3 = l1*l2
      t11 = l3/2
      t12 = Diquad721(1.D0,l1,l2,t11,l4)
      t15 = 1-t3/2
      t16 = 1-l2
      t17 = s34r2(1.D0,t15,t16,t11,l4)
      t20 = s23r2(1.D0,t15,t16,t11,l4)
      t22 = f1(t11,1.D0,t15)
      tncoff4 = -1/(-2+l3)/(-2+t3)/l2/l3/t12/t17/t20*t22/1024
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l4)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t7 = f3(0.D0,z)
      t9 = dlog(1-t7)
      t10 = dlog(2.D0)
      t11 = dlog(z)
      t13 = dlog(l1)
      t15 = dlog(1-l4)
      t18 = -1/t1
      t20 = Diquad721(z,l1,0.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t28 = f1(0.D0,z,1.D0)
      t29 = 1/t23/t25*t28
      tncoff5 = (4*t2+t4-t6+t9+t10-2*t11+t13+t15)*t18*t21*t29/4096-t18/l
     #3*t21*t29/2048-t18/l2*t21*t29/2048
      t8 = l3/2
      t9 = Diquad721(z,l1,0.D0,t8,l4)
      t11 = s34r2(z,1.D0,1.D0,t8,l4)
      t14 = s23r2(z,1.D0,1.D0,t8,l4)
      t16 = f1(t8,z,1.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t9/t11/t14*t16/1024
      t1 = l1*l2
      t9 = Diquad721(z,l1,l2,0.D0,l4)
      t12 = 1-t1/2
      t13 = 1-l2
      t14 = s34r2(z,t12,t13,0.D0,l4)
      t17 = s23r2(z,t12,t13,0.D0,l4)
      t19 = f1(0.D0,z,t12)
      tncoff7 = -1/(-2+t1)/l2/(-1+z)/t9/t14/t17*t19/1024
      t2 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,1.D0,1.D0)
      t11 = 1/t5/t7*t10
      t23 = dlog(l4)
      t25 = dlog(1-l4)
      t26 = f2(0.D0,1.D0)
      t27 = dlog(t26)
      t28 = dlog(2.D0)
      t29 = dlog(l1)
      t30 = f3(0.D0,1.D0)
      t32 = dlog(1-t30)
      tncoff8 = -1/l2*t3*t11/8192-1/l3*t3*t11/8192+1/(-1+z)*t3*t11/4096+
     #(t23/16384+t25/16384-t27/16384+t28/16384+t29/16384+t32/16384)*t3*t
     #11
      t5 = l3/2
      t6 = Diquad721(1.D0,l1,0.D0,t5,l4)
      t9 = s34r2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = f1(t5,1.D0,1.D0)
      tncoff9 = -1/(-2+l3)/l3/t6/t9/t11*t14/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,0.D0,l4)
      t14 = s23r2(1.D0,t10,t11,0.D0,l4)
      t17 = f1(0.D0,1.D0,t10)
      tncoff10 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,0.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,z,1.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,1.D0,1.D0)
      tncoff12 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,l3/2,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qbiq73e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq73e2int(var,ndim,Tcoeff)

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
      double precision qbiq7num1e0
      external qbiq7num1e0
      double precision qbiq7num1e1
      external qbiq7num1e1
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

      t2 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t7 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,1.D0,1.D0)
      t11 = 1/t5/t7*t10
      t23 = dlog(l4)
      t25 = dlog(1-l4)
      t26 = f2(0.D0,1.D0)
      t27 = dlog(t26)
      t28 = dlog(2.D0)
      t29 = dlog(l1)
      t30 = f3(0.D0,1.D0)
      t32 = dlog(1-t30)
      tncoff1 = -1/l2*t3*t11/8192-1/l3*t3*t11/8192+1/(-1+z)*t3*t11/4096+
     #(t23/16384+t25/16384-t27/16384+t28/16384+t29/16384+t32/16384)*t3*t
     #11
      t5 = l3/2
      t6 = Diquad721(1.D0,l1,0.D0,t5,l4)
      t9 = s34r2(1.D0,1.D0,1.D0,t5,l4)
      t11 = s23r2(1.D0,1.D0,1.D0,t5,l4)
      t14 = f1(t5,1.D0,1.D0)
      tncoff2 = -1/(-2+l3)/l3/t6/t9/t11*t14/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,0.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,0.D0,l4)
      t14 = s23r2(1.D0,t10,t11,0.D0,l4)
      t17 = f1(0.D0,1.D0,t10)
      tncoff3 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,0.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t11 = f1(0.D0,z,1.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,1.D0,1.D0)
      tncoff5 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,0.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qbiq73e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq73e3int(var,ndim,Tcoeff)

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
      double precision qbiq7num1e0
      external qbiq7num1e0
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

      t1 = Diquad721(1.D0,l1,0.D0,0.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,0.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,1.D0,1.D0)
      tncoff1 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qbiq73e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq73e4int(var,ndim,Tcoeff)

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

      qbiq73e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

