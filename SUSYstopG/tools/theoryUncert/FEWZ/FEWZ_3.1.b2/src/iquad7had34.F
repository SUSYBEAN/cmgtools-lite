CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qbiq74e0int(var,ndim,Tcoeff)

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

      t1 = 0.3141592653589793D1**2
      t2 = 2*t1
      t3 = dlog(2.D0)
      t4 = t3**2
      t5 = 3*t4
      t6 = f2(1.D0,1.D0)
      t7 = dlog(t6)
      t8 = t7**2
      t9 = 3*t8
      t11 = dlog(1-l4)
      t12 = 1-z
      t13 = dlog(t12)
      t16 = dlog(l4)
      t17 = t11*t16
      t18 = 6*t17
      t19 = dlog(l1)
      t27 = 6*t16*t19
      t28 = t16*t3
      t29 = 6*t28
      t31 = 6*t11*t19
      t32 = t11*t3
      t33 = 6*t32
      t35 = 6*t3*t19
      t38 = -t2+t5+t9+24*t11*t13+t18+24*t13*t19+24*t13*t3+24*t16*t13+t27
     #+t29+t31+t33+t35-24*t13*t7
      t39 = f3(1.D0,1.D0)
      t41 = dlog(1-t39)
      t45 = 6*t41*t19
      t46 = t3*t41
      t47 = 6*t46
      t48 = t16*t41
      t49 = 6*t48
      t50 = t16*t7
      t51 = 6*t50
      t52 = t11*t41
      t53 = 6*t52
      t54 = t11*t7
      t55 = 6*t54
      t57 = 6*t7*t19
      t58 = t7*t3
      t59 = 6*t58
      t60 = t7*t41
      t61 = 6*t60
      t62 = t13**2
      t64 = t19**2
      t65 = 3*t64
      t66 = t41**2
      t67 = 3*t66
      t68 = t16**2
      t69 = 3*t68
      t70 = t11**2
      t71 = 3*t70
      t72 = 24*t13*t41+t45+t47+t49-t51+t53-t55-t57-t59-t61+48*t62+t65+t6
     #7+t69+t71
      t75 = -1/t12
      t77 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t78 = 1/t77
      t80 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t81 = 1/t80
      t82 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t83 = 1/t82
      t85 = f1(1.D0,1.D0,1.D0)
      t86 = t81*t83*t85
      t89 = dlog(l3)
      t90 = t89**2
      t92 = -t2+t5+t9+12*t90+t18+t27+t29+t31+t33+t35+t45+t47+t49-t51
      t105 = t53-t55-t57-t59-t61+12*t89*t19+12*t89*t3+12*t11*t89+12*t16*
     #t89+12*t89*t41-12*t89*t7+t65+t67+t69+t71
      t107 = 1/l3
      t112 = dlog(l2)
      t113 = t112**2
      t115 = -t2+t5+t9+12*t113+t18+t27+t29+t31+t33+t35+t45+t47+t49-t51
      t128 = t53-t55-t57-t59-t61+12*t112*t19+12*t16*t112+12*t3*t112+12*t
     #11*t112-12*t112*t7+12*t112*t41+t65+t67+t69+t71
      t130 = 1/l2
      t139 = t78*t81*t83*t85
      t142 = 4*t13
      t143 = 2*t89
      t149 = 2*t112
      t190 = -t7*t8/98304+t64*t19/98304+t66*t41/98304+t68*t16/98304+t70*
     #t11/98304+t17*t3/16384+t17*t19/16384+t28*t19/16384+t32*t19/16384-t
     #54*t41/16384+t52*t19/16384+t32*t41/16384-t54*t19/16384-t32*t7/1638
     #4+t17*t41/16384
      t223 = -t17*t7/16384-t50*t41/16384+t48*t19/16384+t28*t41/16384-t50
     #*t19/16384-t28*t7/16384-t60*t19/16384+t46*t19/16384-t58*t19/16384-
     #t60*t3/16384+t4*t3/98304+t3*t64/32768+t19*t4/32768-t19*t1/49152-t1
     #*t3/49152+t70*t16/32768
      t257 = t70*t19/32768+t70*t3/32768+t16*t64/32768+t16*t4/32768-t1*t1
     #6/49152+t11*t68/32768+t68*t19/32768+t68*t3/32768+t11*t64/32768+t11
     #*t4/32768-t11*t1/49152-t66*t7/32768+t41*t4/32768+t8*t19/32768-t7*t
     #64/32768+t70*t41/32768
      t289 = -t70*t7/32768+t68*t41/32768-t68*t7/32768+t11*t66/32768+t11*
     #t8/32768+t16*t66/32768+t16*t8/32768+t66*t19/32768+t41*t64/32768-t4
     #1*t1/49152+t41*t8/32768+t3*t8/32768-t7*t4/32768+t1*t7/49152+t66*t3
     #/32768+z3/6144
      tncoff1 = (t38+t72)*t75*t78*t86/24576-(t92+t105)*t107*t78*t86/4915
     #2-(t115+t128)*t130*t78*t86/49152+t75*t130*t107*t139/1024-(t16+t142
     #+t11+t143+t41+t19+t3-t7)*t75*t107*t139/2048-(t3-t7+t19+t41+t149+t1
     #42+t11+t16)*t75*t130*t139/2048+(t11+t16+t143+t3+t41+t149-t7+t19)*t
     #130*t107*t139/4096+(t190+t223+t257+t289)*t78*t86
      t2 = dlog(1-l4)
      t3 = 2-l3
      t4 = dlog(t3)
      t7 = dlog(l4)
      t12 = dlog(l1)
      t15 = dlog(2.D0)
      t18 = dlog(l3)
      t40 = 1-l3/2
      t41 = f2(t40,1.D0)
      t42 = dlog(t41)
      t45 = f3(t40,1.D0)
      t47 = dlog(1-t45)
      t52 = 12*t2*t4+6*t2*t7+12*t7*t4+12*t12*t4-12*t15*t4+24*t18*t4+12*t
     #18*t12+6*t12*t7-6*t7*t15+6*t12*t2-6*t2*t15-6*t15*t12-12*t18*t15+12
     #*t2*t18+12*t7*t18-6*t42*t7+6*t7*t47+6*t2*t47
      t73 = 0.3141592653589793D1**2
      t75 = t15**2
      t77 = t4**2
      t79 = t18**2
      t81 = t12**2
      t83 = t42**2
      t85 = t47**2
      t87 = t7**2
      t89 = t2**2
      t91 = -6*t2*t42+12*t18*t47-12*t18*t42-6*t12*t42-6*t47*t42+6*t47*t1
     #2-6*t47*t15+6*t42*t15+12*t47*t4-12*t42*t4-2*t73+3*t75+12*t77+12*t7
     #9+3*t81+3*t83+3*t85+3*t87+3*t89
      t93 = 1/l3
      t96 = -1/t3
      t98 = Diquad721(1.D0,l1,0.D0,t40,l4)
      t100 = s34r2(1.D0,1.D0,1.D0,t40,l4)
      t103 = s23r2(1.D0,1.D0,1.D0,t40,l4)
      t105 = f1(t40,1.D0,1.D0)
      t107 = 1/t98/t100/t103*t105
      t110 = dlog(l2)
      t112 = 2*t4
      t113 = 2*t18
      t117 = 1/l2*t96
      t121 = -1+z
      t123 = 1/t121*t93
      t128 = dlog(-t121)
      tncoff2 = -(t52+t91)*t93*t96*t107/24576+(t2+t7-t42+t12+t47+2*t110-
     #t15+t112+t113)*t93*t117*t107/2048+t123*t117*t107/512-(t7+4*t128+t2
     #+t113-t15+t112-t42+t12+t47)*t96*t123*t107/1024
      t1 = 1-l2
      t2 = dlog(t1)
      t3 = t2**2
      t5 = 0.3141592653589793D1**2
      t7 = dlog(2.D0)
      t8 = t7**2
      t10 = f2(1.D0,1.D0)
      t11 = dlog(t10)
      t12 = t11**2
      t14 = l1*l2
      t15 = 2-t14
      t16 = dlog(t15)
      t17 = t16**2
      t19 = dlog(l1)
      t20 = t19**2
      t22 = f3(1.D0,1.D0)
      t27 = dlog(2-2*t22+t22*l1*l2)
      t28 = t27**2
      t30 = dlog(l4)
      t31 = t30**2
      t34 = dlog(1-l4)
      t35 = t34**2
      t37 = dlog(l2)
      t38 = t37**2
      t42 = 3*t3-2*t5+3*t8+3*t12+3*t17+3*t20+3*t28+3*t31+3*t35+12*t38-6*
     #t34*t11
      t67 = -6*t11*t19+6*t11*t7+6*t34*t16+6*t34*t2+6*t16*t2+6*t16*t19+12
     #*t37*t16+12*t37*t2+6*t30*t16+6*t30*t2+12*t37*t19+12*t30*t37
      t91 = -12*t7*t37+6*t34*t30+12*t34*t37+6*t30*t19-6*t30*t7+6*t34*t19
     #-6*t34*t7-6*t7*t19-12*t37*t11+12*t37*t27+6*t34*t27
      t104 = t30*t27-t2*t11+t27*t16-t11*t16-t11*t27+t2*t27-t7*t27+t19*t2
     #7-t2*t7-t16*t7+t2*t19-t30*t11
      t109 = -1/t15
      t111 = 1/l2
      t113 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t116 = 1-t14/2
      t117 = s34r2(1.D0,t116,t1,1.D0,l4)
      t120 = s23r2(1.D0,t116,t1,1.D0,l4)
      t122 = f1(1.D0,1.D0,t116)
      t124 = 1/t113/t117/t120*t122
      t127 = -1+z
      t128 = 1/t127
      t129 = 1/l3
      t131 = t109*t111
      t135 = 2*t37
      t137 = dlog(-t127)
      t145 = dlog(l3)
      tncoff3 = -(t42+t67+t91+6*t104)*t109*t111*t124/24576+t128*t129*t13
     #1*t124/512-(-t11+t16+t135+t19+t27-t7+t2+4*t137+t34+t30)*t109*t111*
     #t128*t124/1024+(t34+t30+t19+t2-t11+t135-t7+t16+t27+2*t145)*t129*t1
     #31*t124/2048
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l1)
      t5 = l1*l2
      t6 = 2-t5
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t11 = 1-l3/2
      t12 = f2(t11,1.D0)
      t13 = dlog(t12)
      t14 = 1-l2
      t15 = dlog(t14)
      t16 = dlog(l2)
      t18 = f3(t11,1.D0)
      t23 = dlog(2-2*t18+t18*l1*l2)
      t24 = 2-l3
      t25 = dlog(t24)
      t27 = dlog(l3)
      t29 = t2+t3+t4+t7-3*t8-t13+t15+2*t16+t23+2*t25+2*t27
      t31 = -1/t24
      t34 = -1/t6
      t35 = 1/l2
      t38 = 1/l3
      t39 = Diquad721(1.D0,l1,l2,t11,l4)
      t40 = 1/t39
      t43 = 1-t5/2
      t44 = s34r2(1.D0,t43,t14,t11,l4)
      t45 = 1/t44
      t46 = s23r2(1.D0,t43,t14,t11,l4)
      t47 = 1/t46
      t49 = f1(t11,1.D0,t43)
      tncoff4 = t29*t31*t34*t35*t38*t40*t45*t47*t49/1024+t45*t40*t31*t34
     #*t47*t35*t38/(-1+z)*t49/256
      t2 = dlog(1-l4)
      t3 = 1-z
      t4 = dlog(t3)
      t7 = dlog(l4)
      t10 = dlog(l1)
      t13 = dlog(2.D0)
      t18 = f3(1.D0,z)
      t20 = dlog(1-t18)
      t23 = f2(1.D0,z)
      t24 = dlog(t23)
      t35 = dlog(z)
      t50 = 24*t2*t4+6*t2*t7+24*t4*t10+24*t4*t13+24*t7*t4+6*t2*t20-6*t2*
     #t24+6*t7*t20-6*t7*t24+24*t4*t20-24*t4*t24-12*t20*t35-6*t20*t24-6*t
     #24*t13+12*t24*t35+6*t20*t13+6*t20*t10-6*t24*t10
      t71 = t24**2
      t73 = t4**2
      t75 = t10**2
      t77 = t20**2
      t79 = t35**2
      t81 = t7**2
      t83 = t2**2
      t85 = 0.3141592653589793D1**2
      t87 = t13**2
      t89 = 6*t7*t10+6*t13*t7+6*t2*t10+6*t2*t13+6*t13*t10-12*t13*t35-12*
     #t2*t35-12*t35*t10-12*t7*t35-48*t4*t35+3*t71+48*t73+3*t75+3*t77+12*
     #t79+3*t81+3*t83-2*t85+3*t87
      t92 = -1/t3
      t94 = Diquad721(z,l1,0.D0,1.D0,l4)
      t95 = 1/t94
      t97 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t98 = 1/t97
      t99 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t100 = 1/t99
      t102 = f1(1.D0,z,1.D0)
      t106 = 1/l2
      t108 = 1/l3
      t112 = t95*t98*t100*t102
      t115 = 4*t4
      t116 = dlog(l3)
      t118 = 2*t35
      t124 = dlog(l2)
      tncoff5 = -(t50+t89)*t92*t95*t98*t100*t102/24576-t92*t106*t108*t11
     #2/1024+(t7+t115+t2+2*t116+t20-t24+t13+t10-t118)*t92*t108*t112/2048
     #+(-t24+2*t124+t10+t20-t118+t13+t115+t2+t7)*t92*t106*t112/2048
      t1 = -1+z
      t4 = 1/t1/l3
      t6 = -2+l3
      t7 = 1/t6
      t11 = 1-l3/2
      t12 = Diquad721(z,l1,0.D0,t11,l4)
      t14 = s34r2(z,1.D0,1.D0,t11,l4)
      t17 = s23r2(z,1.D0,1.D0,t11,l4)
      t19 = f1(t11,z,1.D0)
      t21 = 1/t12/t14/t17*t19
      t24 = dlog(l4)
      t26 = dlog(-t1)
      t29 = dlog(1-l4)
      t30 = f3(t11,z)
      t32 = dlog(1-t30)
      t33 = f2(t11,z)
      t34 = dlog(t33)
      t35 = dlog(2.D0)
      t37 = dlog(-t6)
      t39 = dlog(z)
      t41 = dlog(l1)
      t42 = dlog(l3)
      tncoff6 = -t4/l2*t7*t21/512+(t24+4*t26+t29+t32-t34-t35+2*t37-2*t39
     #+t41+2*t42)*t7*t4*t21/1024
      t1 = -1+z
      t2 = 1/t1
      t5 = l1*l2
      t6 = -2+t5
      t7 = 1/t6
      t8 = 1/l2
      t11 = Diquad721(z,l1,l2,1.D0,l4)
      t14 = 1-t5/2
      t15 = 1-l2
      t16 = s34r2(z,t14,t15,1.D0,l4)
      t19 = s23r2(z,t14,t15,1.D0,l4)
      t21 = f1(1.D0,z,t14)
      t23 = 1/t11/t16/t19*t21
      t26 = dlog(l1)
      t27 = dlog(l2)
      t29 = f2(1.D0,z)
      t30 = dlog(t29)
      t31 = dlog(z)
      t33 = dlog(2.D0)
      t34 = f3(1.D0,z)
      t39 = dlog(2-2*t34+t34*l1*l2)
      t41 = dlog(-t6)
      t42 = dlog(t15)
      t44 = dlog(-t1)
      t47 = dlog(1-l4)
      t48 = dlog(l4)
      t49 = t26+2*t27-t30-2*t31-t33+t39+t41+t42+4*t44+t47+t48
      tncoff7 = -t2/l3*t7*t8*t23/512+t49*t7*t8*t2*t23/1024
      t2 = 1-l3/2
      t3 = l1*l2
      t5 = 1-t3/2
      t6 = f1(t2,z,t5)
      t7 = 1-l2
      t8 = s34r2(z,t5,t7,t2,l4)
      t11 = Diquad721(z,l1,l2,t2,l4)
      t19 = s23r2(z,t5,t7,t2,l4)
      tncoff8 = -t6/t8/t11/(-2+l3)/(-2+t3)/t19/l2/l3/(-1+z)/256
      t2 = dlog(1-l4)
      t3 = dlog(l1)
      t6 = dlog(2.D0)
      t9 = dlog(l4)
      t12 = f3(1.D0,1.D0)
      t14 = dlog(1-t12)
      t17 = f2(1.D0,1.D0)
      t18 = dlog(t17)
      t31 = t6**2
      t33 = -t2*t3/16384-t2*t6/16384-t2*t9/16384-t9*t14/16384+t9*t18/163
     #84-t9*t3/16384-t9*t6/16384-t14*t3/16384+t18*t14/16384+t6*t18/16384
     #-t31/32768
      t34 = 0.3141592653589793D1**2
      t36 = t18**2
      t38 = t14**2
      t40 = t3**2
      t42 = t2**2
      t44 = t9**2
      t56 = t34/49152-t36/32768-t38/32768-t40/32768-t42/32768-t44/32768-
     #t6*t14/16384-t6*t3/16384+t18*t3/16384-t2*t14/16384+t2*t18/16384
      t58 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t59 = 1/t58
      t61 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t63 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t66 = f1(1.D0,1.D0,1.D0)
      t67 = 1/t61/t63*t66
      t69 = 1-z
      t70 = dlog(t69)
      t72 = 6*t9
      t73 = 6*t2
      t74 = 6*t3
      t75 = 6*t6
      t76 = 6*t18
      t77 = 6*t14
      t80 = -1/t69
      t85 = dlog(l3)
      t88 = 1/l3
      t93 = dlog(l2)
      t96 = 1/l2
      tncoff9 = (t33+t56)*t59*t67+(-24*t70-t72-t73-t74-t75+t76-t77)*t80*
     #t59*t67/24576-(-t72-t73-t74-t75+t76-t77-12*t85)*t88*t59*t67/49152-
     #(-t72-t73-t74-t75+t76-t77-12*t93)*t96*t59*t67/49152+t80*t88*t59*t6
     #7/2048+t80*t96*t59*t67/2048-t96*t88*t59*t67/4096
      t1 = 1/l3
      t4 = -2+l3
      t5 = 1/t4
      t8 = 1-l3/2
      t9 = Diquad721(1.D0,l1,0.D0,t8,l4)
      t11 = s34r2(1.D0,1.D0,1.D0,t8,l4)
      t14 = s23r2(1.D0,1.D0,1.D0,t8,l4)
      t16 = f1(t8,1.D0,1.D0)
      t18 = 1/t9/t11/t14*t16
      t27 = dlog(2.D0)
      t29 = f2(t8,1.D0)
      t30 = dlog(t29)
      t32 = f3(t8,1.D0)
      t34 = dlog(1-t32)
      t37 = dlog(-t4)
      t39 = dlog(l1)
      t41 = dlog(l4)
      t44 = dlog(1-l4)
      t46 = dlog(l3)
      tncoff10 = -t1/l2*t5*t18/2048+t5*t1/(-1+z)*t18/1024-(6*t27+6*t30-6
     #*t34-12*t37-6*t39-6*t41-6*t44-12*t46)*t1*t5*t18/24576
      t2 = l1*l2
      t3 = -2+t2
      t4 = 1/t3
      t6 = 1/l2
      t8 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t11 = 1-t2/2
      t12 = 1-l2
      t13 = s34r2(1.D0,t11,t12,1.D0,l4)
      t16 = s23r2(1.D0,t11,t12,1.D0,l4)
      t18 = f1(1.D0,1.D0,t11)
      t20 = 1/t8/t13/t16*t18
      t29 = dlog(l2)
      t32 = dlog(1-l4)
      t34 = dlog(l4)
      t36 = dlog(l1)
      t38 = dlog(2.D0)
      t40 = dlog(t12)
      t42 = f2(1.D0,1.D0)
      t43 = dlog(t42)
      t45 = f3(1.D0,1.D0)
      t50 = dlog(2-2*t45+t45*l1*l2)
      t53 = dlog(-t3)
      tncoff11 = -1/l3*t4*t6*t20/2048+t4*t6/(-1+z)*t20/1024-(-12*t29-6*t
     #32-6*t34-6*t36+6*t38-6*t40+6*t43-6*t50-6*t53)*t4*t6*t20/24576
      t3 = l1*l2
      t12 = 1-l3/2
      t13 = Diquad721(1.D0,l1,l2,t12,l4)
      t16 = 1-t3/2
      t17 = 1-l2
      t18 = s34r2(1.D0,t16,t17,t12,l4)
      t21 = s23r2(1.D0,t16,t17,t12,l4)
      t23 = f1(t12,1.D0,t16)
      tncoff12 = -1/(-2+l3)/(-2+t3)/l2/l3/t13/t18/t21*t23/1024
      t1 = f3(1.D0,z)
      t3 = dlog(1-t1)
      t5 = f2(1.D0,z)
      t6 = dlog(t5)
      t8 = dlog(z)
      t10 = dlog(2.D0)
      t12 = dlog(l1)
      t15 = dlog(1-l4)
      t17 = dlog(l4)
      t19 = 1-z
      t20 = dlog(t19)
      t24 = -1/t19
      t26 = Diquad721(z,l1,0.D0,1.D0,l4)
      t27 = 1/t26
      t29 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t31 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t34 = f1(1.D0,z,1.D0)
      t35 = 1/t29/t31*t34
      tncoff13 = -(-6*t3+6*t6+12*t8-6*t10-6*t12-6*t15-6*t17-24*t20)*t24*
     #t27*t35/24576-t24/l3*t27*t35/2048-t24/l2*t27*t35/2048
      t9 = 1-l3/2
      t10 = Diquad721(z,l1,0.D0,t9,l4)
      t12 = s34r2(z,1.D0,1.D0,t9,l4)
      t15 = s23r2(z,1.D0,1.D0,t9,l4)
      t17 = f1(t9,z,1.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t10/t12/t15*t17/1024
      t1 = l1*l2
      t9 = Diquad721(z,l1,l2,1.D0,l4)
      t12 = 1-t1/2
      t13 = 1-l2
      t14 = s34r2(z,t12,t13,1.D0,l4)
      t17 = s23r2(z,t12,t13,1.D0,l4)
      t19 = f1(1.D0,z,t12)
      tncoff15 = -1/(-2+t1)/l2/(-1+z)/t9/t14/t17*t19/1024
      t3 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t4 = 1/t3
      t6 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      t12 = 1/t6/t8*t11
      t23 = dlog(2.D0)
      t24 = f2(1.D0,1.D0)
      t25 = dlog(t24)
      t26 = f3(1.D0,1.D0)
      t28 = dlog(1-t26)
      t30 = dlog(1-l4)
      t31 = dlog(l4)
      t32 = dlog(l1)
      tncoff16 = 1/(-1+z)*t4*t12/4096-1/l3*t4*t12/8192-1/l2*t4*t12/8192+
     #(t23/16384-t25/16384+t28/16384+t30/16384+t31/16384+t32/16384)*t4*t
     #12
      t6 = 1-l3/2
      t7 = Diquad721(1.D0,l1,0.D0,t6,l4)
      t10 = s34r2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = f1(t6,1.D0,1.D0)
      tncoff17 = -1/l3/(-2+l3)/t7/t10/t12*t15/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,1.D0,l4)
      t14 = s23r2(1.D0,t10,t11,1.D0,l4)
      t17 = f1(1.D0,1.D0,t10)
      tncoff18 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,1.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,z,1.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,1.D0,1.D0)
      tncoff20 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5,
     #l6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5,l6,
     #l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5
     #,l6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qbiq7num1e1(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l
     #7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qbiq7num1e2(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qbiq7num1e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qbiq7num1e3(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qbiq74e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq74e1int(var,ndim,Tcoeff)

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
      t3 = dlog(l1)
      t6 = dlog(2.D0)
      t9 = dlog(l4)
      t12 = f3(1.D0,1.D0)
      t14 = dlog(1-t12)
      t17 = f2(1.D0,1.D0)
      t18 = dlog(t17)
      t31 = t6**2
      t33 = -t2*t3/16384-t2*t6/16384-t2*t9/16384-t9*t14/16384+t9*t18/163
     #84-t9*t3/16384-t9*t6/16384-t14*t3/16384+t18*t14/16384+t6*t18/16384
     #-t31/32768
      t34 = 0.3141592653589793D1**2
      t36 = t18**2
      t38 = t14**2
      t40 = t3**2
      t42 = t2**2
      t44 = t9**2
      t56 = t34/49152-t36/32768-t38/32768-t40/32768-t42/32768-t44/32768-
     #t6*t14/16384-t6*t3/16384+t18*t3/16384-t2*t14/16384+t2*t18/16384
      t58 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t59 = 1/t58
      t61 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t63 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t66 = f1(1.D0,1.D0,1.D0)
      t67 = 1/t61/t63*t66
      t69 = 1-z
      t70 = dlog(t69)
      t74 = -1/t69
      t79 = dlog(l3)
      t82 = 1/l3
      t91 = dlog(l2)
      t94 = 1/l2
      tncoff1 = (t33+t56)*t59*t67-(t14+t3+t6-t18+4*t70+t9+t2)*t74*t59*t6
     #7/4096+(t2+t9+t14+t3+t6-t18+2*t79)*t82*t59*t67/8192+t74*t82*t59*t6
     #7/2048+(t2+t9+t6-t18+t3+t14+2*t91)*t94*t59*t67/8192+t74*t94*t59*t6
     #7/2048-t94*t82*t59*t67/4096
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(2.D0)
      t5 = 2-l3
      t6 = dlog(t5)
      t8 = dlog(l1)
      t10 = 1-l3/2
      t11 = f3(t10,1.D0)
      t13 = dlog(1-t11)
      t14 = f2(t10,1.D0)
      t15 = dlog(t14)
      t16 = dlog(l3)
      t19 = 1/l3
      t22 = -1/t5
      t24 = Diquad721(1.D0,l1,0.D0,t10,l4)
      t26 = s34r2(1.D0,1.D0,1.D0,t10,l4)
      t29 = s23r2(1.D0,1.D0,1.D0,t10,l4)
      t31 = f1(t10,1.D0,1.D0)
      t33 = 1/t24/t26/t29*t31
      tncoff2 = (t2+t3-t4+2*t6+t8+t13-t15+2*t16)*t19*t22*t33/4096+t22*t1
     #9/(-1+z)*t33/1024-t19/l2*t22*t33/2048
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = f2(1.D0,1.D0)
      t5 = dlog(t4)
      t6 = l1*l2
      t7 = 2-t6
      t8 = dlog(t7)
      t9 = dlog(l2)
      t11 = dlog(l1)
      t12 = f3(1.D0,1.D0)
      t17 = dlog(2-2*t12+t12*l1*l2)
      t18 = dlog(2.D0)
      t19 = 1-l2
      t20 = dlog(t19)
      t23 = -1/t7
      t25 = 1/l2
      t27 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t30 = 1-t6/2
      t31 = s34r2(1.D0,t30,t19,1.D0,l4)
      t34 = s23r2(1.D0,t30,t19,1.D0,l4)
      t36 = f1(1.D0,1.D0,t30)
      t38 = 1/t27/t31/t34*t36
      tncoff3 = (t2+t3-t5+t8+2*t9+t11+t17-t18+t20)*t23*t25*t38/4096+t23*
     #t25/(-1+z)*t38/1024-1/l3*t23*t25*t38/2048
      t3 = l1*l2
      t12 = 1-l3/2
      t13 = Diquad721(1.D0,l1,l2,t12,l4)
      t16 = 1-t3/2
      t17 = 1-l2
      t18 = s34r2(1.D0,t16,t17,t12,l4)
      t21 = s23r2(1.D0,t16,t17,t12,l4)
      t23 = f1(t12,1.D0,t16)
      tncoff4 = -1/(-2+l3)/(-2+t3)/l2/l3/t13/t18/t21*t23/1024
      t1 = f3(1.D0,z)
      t3 = dlog(1-t1)
      t4 = dlog(2.D0)
      t5 = dlog(z)
      t7 = f2(1.D0,z)
      t8 = dlog(t7)
      t9 = dlog(l1)
      t10 = 1-z
      t11 = dlog(t10)
      t13 = dlog(l4)
      t15 = dlog(1-l4)
      t18 = -1/t10
      t20 = Diquad721(z,l1,0.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = f1(1.D0,z,1.D0)
      t29 = 1/t23/t25*t28
      tncoff5 = (t3+t4-2*t5-t8+t9+4*t11+t13+t15)*t18*t21*t29/4096-t18/l3
     #*t21*t29/2048-t18/l2*t21*t29/2048
      t9 = 1-l3/2
      t10 = Diquad721(z,l1,0.D0,t9,l4)
      t12 = s34r2(z,1.D0,1.D0,t9,l4)
      t15 = s23r2(z,1.D0,1.D0,t9,l4)
      t17 = f1(t9,z,1.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t10/t12/t15*t17/1024
      t1 = l1*l2
      t9 = Diquad721(z,l1,l2,1.D0,l4)
      t12 = 1-t1/2
      t13 = 1-l2
      t14 = s34r2(z,t12,t13,1.D0,l4)
      t17 = s23r2(z,t12,t13,1.D0,l4)
      t19 = f1(1.D0,z,t12)
      tncoff7 = -1/(-2+t1)/l2/(-1+z)/t9/t14/t17*t19/1024
      t3 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t4 = 1/t3
      t6 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      t12 = 1/t6/t8*t11
      t23 = dlog(2.D0)
      t24 = f2(1.D0,1.D0)
      t25 = dlog(t24)
      t26 = f3(1.D0,1.D0)
      t28 = dlog(1-t26)
      t30 = dlog(1-l4)
      t31 = dlog(l4)
      t32 = dlog(l1)
      tncoff8 = 1/(-1+z)*t4*t12/4096-1/l3*t4*t12/8192-1/l2*t4*t12/8192+(
     #t23/16384-t25/16384+t28/16384+t30/16384+t31/16384+t32/16384)*t4*t1
     #2
      t6 = 1-l3/2
      t7 = Diquad721(1.D0,l1,0.D0,t6,l4)
      t10 = s34r2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = f1(t6,1.D0,1.D0)
      tncoff9 = -1/l3/(-2+l3)/t7/t10/t12*t15/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,1.D0,l4)
      t14 = s23r2(1.D0,t10,t11,1.D0,l4)
      t17 = f1(1.D0,1.D0,t10)
      tncoff10 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,1.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,z,1.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,1.D0,1.D0)
      tncoff12 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1-l3/2,l4,l5,
     #l6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qbiq7num1e0(var,ndim,z,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qbiq7num1e1(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l
     #6,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qbiq7num1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qbiq7num1e2(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qbiq74e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq74e2int(var,ndim,Tcoeff)

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

      t3 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t4 = 1/t3
      t6 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,1.D0)
      t12 = 1/t6/t8*t11
      t23 = dlog(2.D0)
      t24 = f2(1.D0,1.D0)
      t25 = dlog(t24)
      t26 = f3(1.D0,1.D0)
      t28 = dlog(1-t26)
      t30 = dlog(1-l4)
      t31 = dlog(l4)
      t32 = dlog(l1)
      tncoff1 = 1/(-1+z)*t4*t12/4096-1/l3*t4*t12/8192-1/l2*t4*t12/8192+(
     #t23/16384-t25/16384+t28/16384+t30/16384+t31/16384+t32/16384)*t4*t1
     #2
      t6 = 1-l3/2
      t7 = Diquad721(1.D0,l1,0.D0,t6,l4)
      t10 = s34r2(1.D0,1.D0,1.D0,t6,l4)
      t12 = s23r2(1.D0,1.D0,1.D0,t6,l4)
      t15 = f1(t6,1.D0,1.D0)
      tncoff2 = -1/l3/(-2+l3)/t7/t10/t12*t15/4096
      t1 = l1*l2
      t6 = Diquad721(1.D0,l1,l2,1.D0,l4)
      t10 = 1-t1/2
      t11 = 1-l2
      t12 = s34r2(1.D0,t10,t11,1.D0,l4)
      t14 = s23r2(1.D0,t10,t11,1.D0,l4)
      t17 = f1(1.D0,1.D0,t10)
      tncoff3 = -1/(-2+t1)/l2/t6/t12/t14*t17/4096
      t3 = Diquad721(z,l1,0.D0,1.D0,l4)
      t6 = s34r2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t11 = f1(1.D0,z,1.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t8*t11/4096
      t1 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,1.D0,1.D0)
      tncoff5 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7
     #,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qbiq7num1e0(var,ndim,1.D0,ux,1-l1*l2/2,1-l2,1.D0,l4,l5,l6
     #,l7,l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qbiq7num1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qbiq7num1e1(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qbiq74e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq74e3int(var,ndim,Tcoeff)

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

      t1 = Diquad721(1.D0,l1,0.D0,1.D0,l4)
      t3 = s34r2(1.D0,1.D0,1.D0,1.D0,l4)
      t6 = s23r2(1.D0,1.D0,1.D0,1.D0,l4)
      t8 = f1(1.D0,1.D0,1.D0)
      tncoff1 = -1/t1/t3/t6*t8/16384

CCC Numerator functions

      tcnum1 = qbiq7num1e0(var,ndim,1.D0,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qbiq74e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qbiq74e4int(var,ndim,Tcoeff)

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

      qbiq74e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

