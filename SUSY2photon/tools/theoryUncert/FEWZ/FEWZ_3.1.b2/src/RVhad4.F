CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV4e0int(var,ndim,Tcoeff)

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
      double precision RVnum3e0
      external RVnum3e0
      double precision RVnum3e1
      external RVnum3e1
      double precision RVnum3e2
      external RVnum3e2
      double precision RVnum3e3
      external RVnum3e3
      double precision RVnum3e4
      external RVnum3e4
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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

      t1 = dlog(l3)
      t2 = dlog(2.D0)
      t3 = t1*t2
      t4 = 24*t3
      t5 = dlog(l1)
      t6 = t5*t2
      t7 = t6*l1
      t8 = 24*t7
      t9 = 48*t6
      t10 = t5**2
      t11 = t10*l1
      t12 = 12*t11
      t13 = 0.3141592653589793D1**2
      t14 = t13*l1
      t15 = 4*t14
      t16 = t2**2
      t17 = t16*l1
      t18 = 12*t17
      t19 = t2*l1
      t20 = 24*t19
      t21 = t1*l1
      t23 = t5*l1
      t24 = 24*t23
      t25 = t1**2
      t28 = t1*t5
      t30 = -t4+t8-t9-t12+t15-t18-t20+12*t21+t24-3*t25*l1+24*t28
      t31 = 24*t16
      t32 = 8*t13
      t33 = 48*t2
      t34 = 24*t1
      t35 = 6*t25
      t38 = 24*l1
      t39 = 48*t5
      t40 = 24*t10
      t43 = 48+t31-t32+t33-t34+t35-12*t28*l1-t38-t39+t40+12*t3*l1
      t46 = 1/(l1-2)
      t48 = 1/l1
      t49 = 1/l3
      t50 = t48*t49
      t53 = 1-z
      t54 = dlog(t53)
      t55 = t54*l1
      t67 = -1/t53
      t71 = t13*t2
      t73 = 8*t71
      t80 = 16*z3
      t81 = t16*t2
      t82 = 8*t81
      t84 = t4+t73-4*t1*t13+12*t1*t16-6*t25*t2-t31+t32-t33-48+t34-t35+t8
     #0-t82+t25*t1
      t94 = t5*t13
      t97 = t5*t16
      t100 = t10*t2
      t106 = -48-t8-t24+t20+t18-t15+t12+t9+4*t94*l1-12*t97*l1+12*t100*l1
     #-4*t71*l1+t73-8*t94
      t110 = t10*t5
      t117 = 4*t81*l1-24*t100-4*t110*l1-8*z3*l1+24*t97+t38-t33-t40+t39+t
     #32-t31-t82+t80+8*t110
      t122 = t54*t2
      t123 = 48*t122
      t128 = t54**2
      t132 = 12*t16
      t133 = 4*t13
      t134 = 48*t54
      t135 = 24*t2
      t136 = 48*t128
      t141 = t123-8*t54*t13+24*t54*t16-48*t128*t2+4*t71-t132+t133+t134-t
     #135-t136-24+32*t128*t54+8*z3-4*t81
      t146 = t54*t5
      t159 = -24*t146*l1+24*t122*l1-t123+12*t7+24*t55-24*t6-6*t11+2*t14-
     #6*t17-12*t19+12*t23
      t166 = 24+48*t146-24*t128*l1+t132-t133-t134+t135+t136-12*l1-24*t5+
     #12*t10
      t172 = t16**2
      t174 = t13**2
      t181 = 24*t54*t1-12*t3-t123+t132-t133-t134+t135+t136+24-12*t1+3*t2
     #5
      tncoff1 = (t30+t43)*t46*t50/3072-(4*t55-2*t19+t21+2*t23-8*t54+4*t2
     #+4-2*t1-2*l1-4*t5)*t46*t50*t67/128+t71/384+t84*t49/6144-t16/128+t1
     #3/384-t2/64+t13*t16/768+t2*z3/192-(t106+t117)*t46*t48/1536-1.D0/64
     #.D0-t141*t67/768+z3/192-t81/384-(t159+t166)*t46*t48*t67/384-t172/1
     #536-t174/10240+t181*t49*t67/1536
      t1 = dlog(l1)
      t2 = dlog(2.D0)
      t3 = t1*t2
      t6 = dlog(l3)
      t7 = t6*t1
      t10 = 24*t6
      t12 = dlog(1-l3)
      t13 = t12*t2
      t16 = t1*t12
      t19 = t6*t12
      t22 = 24*t12
      t23 = 12*t19
      t24 = t1**2
      t27 = t2*l1
      t29 = t6*t2
      t30 = 24*t29
      t31 = t2**2
      t34 = 0.3141592653589793D1**2
      t38 = -48-24*t3*l1+12*t7*l1+t10+12*t13*l1-12*t16*l1-6*t19*l1-t22+t
     #23-24*t24-24*t7+24*t27+t30+12*t31*l1-4*t34*l1+24*t16
      t39 = t12*l1
      t41 = t12**2
      t44 = 24*t13
      t45 = t6**2
      t48 = t6*l1
      t52 = 6*t45
      t53 = 6*t41
      t57 = t1*l1
      t59 = 48*t2
      t60 = 8*t34
      t63 = 24*t31
      t64 = 12*t39+3*t41*l1-t44+3*t45*l1-12*t48+12*t24*l1-t52-t53-12*t29
     #*l1+48*t3-24*t57-t59+t60+48*t1+24*l1-t63
      t67 = 1/(l1-2)
      t70 = 1/l3
      t71 = 1/l1*t70
      t75 = 1-z
      t76 = dlog(t75)
      t93 = t76**2
      t95 = 12*t6-24*t76*t6+24*t76*t12+48*t76*t2+6*t19+12*t29-12*t13-24*
     #t2-12*t12-3*t45-3*t41-12*t31+4*t34-24+48*t76-48*t93
      t98 = -1/t75
      t111 = 2*t6+2*l1-4*t76*l1-t48-2*t57+t39+4*t1+2*t27-4*t2-2*t12-4+8*
     #t76
      t131 = 48+6*t45*t2+4*t6*t34-t10+t22+6*t41*t2-t23-12*t19*t2-4*t12*t
     #34-8*t34*t2-t45*t6-12*t6*t31
      t142 = t41*t12-3*t6*t41+8*t31*t2-t30-16*z3+12*t12*t31+t44+t52+t53+
     #3*t45*t12+t59-t60+t63
      tncoff2 = (t38+t64)*t67*t71/3072+t95*t70*t98/1536-t111*t67*t71*t98
     #/128+(t131+t142)*t70/6144
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l3)
      t6 = dlog(2.D0)
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = dlog(l1)
      t15 = -1/t8
      t17 = 1/l1
      t19 = t17/l3
      t21 = -1/t1
      t25 = 0.3141592653589793D1**2
      t26 = 8*t25
      t27 = t6**2
      t30 = 96*t27
      t31 = t11**2
      t32 = 24*t31
      t33 = t11*t6
      t34 = 96*t33
      t35 = t9**2
      t36 = 24*t35
      t37 = t11*t9
      t48 = 48-t26+64*t27*t6+t30+t32-t34+t36+96*t37*t6-16*t25*t6+48*t6*t
     #31+8*t11*t25-96*t11*t27
      t54 = 96*t6
      t55 = 48*t9
      t56 = 48*t11
      t57 = 48*t37
      t62 = t9*t6
      t63 = 96*t62
      t70 = -16*z3-8*t35*t9-8*t31*t11+t54-t55-t56+t57-24*t31*t9-24*t11*t
     #35-t63+8*t25*t9-96*t9*t27+48*t6*t35
      t78 = t4**2
      t84 = 24*t4-t54+t55+t56-t32+t26-t30-t57-24*t4*t9+t63+t34-t36-6*t78
     #+48*t4*t6-24*t4*t11-48
      t105 = t2**2
      t107 = 48*t2-48*t6+24*t9+24*t11-12*t31+4*t25-48*t27-24*t37+48*t62+
     #48*t33-12*t35+96*t6*t2-48*t2*t9-48*t2*t11-48*t105-24
      tncoff3 = -(8*t2+2*t4-8*t6+4*t9+4*t11-4)*t15*t19*t21/128-(t48+t70)
     #*t15*t17/1536+t84*t15*t19/3072-t107*t15*t17*t21/384
      t1 = dlog(l3)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(2.D0)
      t8 = dlog(l1)
      t14 = 1-l3
      t16 = dlog(4*t14)
      t19 = t1**2
      t21 = t8**2
      t25 = 48-24*t1-48*t4+48*t6-48*t8-24*t6*t1+24*t1*t4-12*t1*t16+6*t19
     #+24*t21-24*t8*t16
      t30 = t16**2
      t36 = t4**2
      t39 = 0.3141592653589793D1**2
      t41 = t6**2
      t47 = 48*t8*t4+24*t6*t16+6*t30-48*t4*t6-24*t4*t16+24*t36+24*t16-8*
     #t39+24*t41+24*t1*t8-48*t8*t6
      t50 = -1/t3
      t54 = 1/l3/l1
      t57 = 1-z
      t58 = dlog(t57)
      t63 = dlog(t14)
      tncoff4 = (t25+t47)*t50*t54/3072+(4-8*t58-2*t1-4*t4+8*t6+2*t63-4*t
     #8)*t50*t54/t57/128
      t1 = dlog(l3)
      t2 = dlog(z)
      t5 = dlog(2.D0)
      t8 = 1-z
      t9 = dlog(t8)
      t10 = t9*t5
      t11 = 48*t10
      t12 = t2*t9
      t13 = 48*t12
      t14 = t5*t2
      t15 = 24*t14
      t18 = t9**2
      t19 = 48*t18
      t20 = 48*t9
      t21 = 24*t5
      t22 = 24*t2
      t23 = t5**2
      t24 = 12*t23
      t25 = t2**2
      t26 = 12*t25
      t27 = 0.3141592653589793D1**2
      t28 = 4*t27
      t29 = t1**2
      t32 = 12*t1*t2+12*t1*t5+t11+t13-t15-24*t9*t1-t19+t20-t21-t22-t24-t
     #26+t28-3*t29+12*t1-24
      t33 = 1/l3
      t36 = -1/t8
      t42 = 24-t13+t15-t11+t21-t20+4*t23*t5+t19-8*z3+t24+t26-t28
      t67 = -4*t27*t5+8*t9*t27-24*t9*t23+48*t18*t5+12*t23*t2-24*t25*t9+4
     #8*t18*t2+12*t5*t25-4*t27*t2+t22-32*t18*t9+4*t25*t2-48*t10*t2
      t75 = dlog(l1)
      t79 = t75*l1
      t81 = t9*l1
      t83 = t5*l1
      t85 = t2*l1
      t87 = 8*t9-4*t5-4*t2+2*t1+4*t75+2*l1-t1*l1-2*t79-4*t81+2*t83+2*t85
     #-4
      t89 = 1/(l1-2)
      t91 = 1/l1
      t96 = t75*t5
      t99 = t75**2
      t105 = -24+t13-t15+t11-12*t96*l1-t21+t20-12*t99+12*t85+12*l1+24*t7
     #5-t19-t24-t26+t28-12*t79
      t108 = t9*t75
      t113 = t75*t2
      t133 = -24*t81+12*t83+24*t108*l1-24*t10*l1-12*t113*l1+12*t14*l1-24
     #*t12*l1+24*t113+24*t96+6*t23*l1+6*t99*l1-2*t27*l1-48*t108+24*t18*l
     #1+6*t25*l1-t22
      tncoff5 = t32*t33*t36/1536-(t42+t67)*t36/768-t87*t89*t91*t33*t36/1
     #28-(t105+t133)*t89*t91*t36/384
      t1 = dlog(l3)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t8 = dlog(2.D0)
      t12 = dlog(-z*(-1+l3))
      t15 = dlog(l1)
      t26 = -2*t1+2*t3-8*t6+4*t8+2*t12-t12*l1+2*t15*l1+t1*l1-2*t8*l1-t3*
     #l1+4*t6*l1-2*l1-4*t15+4
      t31 = 1/l3
      t34 = -1/t5
      t42 = t12**2
      t44 = 0.3141592653589793D1**2
      t46 = t8**2
      t48 = t3**2
      t50 = t6**2
      t60 = -24*t6*t3+24*t6*t1+3*t42-4*t44+12*t46+3*t48+48*t50-48*t8*t6-
     #24*t6*t12+6*t3*t12-6*t1*t12
      t65 = t1**2
      t76 = 24-12*t1*t8-12*t1+12*t3+3*t65-48*t6+24*t8+12*t12+12*t12*t8+1
     #2*t8*t3-6*t1*t3
      tncoff6 = -t26/(l1-2)/l1*t31*t34/128+(t60+t76)*t31*t34/1536
      t1 = z*l1
      t6 = dlog(-(t1+2-l1)*(t1-2*z-l1))
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = dlog(l1)
      t13 = dlog(l3)
      t15 = 1-z
      t16 = dlog(t15)
      t18 = dlog(2.D0)
      t20 = dlog(z)
      t24 = -1/t8
      t26 = 1/l1
      t30 = -1/t15
      t37 = t11**2
      t39 = t9**2
      t41 = t16**2
      t57 = 24+12*t6-24*t9-24*t11+12*t37+12*t39+48*t41-12*t11*t20+24*t11
     #*t9-24*t18*t9+48*t16*t9-24*t11*t18-48*t16*t18+12*t18*t6
      t77 = t20**2
      t79 = 0.3141592653589793D1**2
      t81 = t18**2
      t83 = t6**2
      t85 = -12*t20*t9-12*t9*t6+6*t20*t6-12*t11*t6-24*t6*t16+48*t16*t11-
     #24*t16*t20+12*t18*t20-48*t16+24*t18+12*t20+3*t77-4*t79+12*t81+3*t8
     #3
      tncoff7 = -(2*t6-4*t9-4*t11-2*t13+4-8*t16+4*t18+2*t20)*t24*t26/l3*
     #t30/128-(t57+t85)*t24*t26*t30/384
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l3)
      t6 = dlog(z)
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = dlog(l1)
      t13 = dlog(2.D0)
      t21 = z**2
      t24 = l1**2
      t29 = dlog(4*z-4*z*l3+2*l1-4*z*l1+2*t21*l1-t24+2*z*t24-t21*t24)
      tncoff8 = -(-4+8*t2+2*t4-2*t6+4*t9+4*t11-4*t13-2*t29)/t8/l1/l3/t1/
     #128
      t1 = 1-z
      t2 = dlog(t1)
      t3 = 24*t2
      t4 = dlog(l3)
      t6 = dlog(2.D0)
      t7 = 12*t6
      t9 = 1/l3
      t12 = -1/t1
      t17 = t6**2
      t19 = t2**2
      t21 = 0.3141592653589793D1**2
      t28 = 2-l1
      t30 = -1/t28
      t32 = 1/l1
      t33 = t32*t9
      t39 = 12*t4
      t40 = 12*t17
      t41 = 4*t21
      t44 = t4**2
      t46 = 24*t6
      t52 = 12*l1
      t53 = t6*l1
      t54 = 12*t53
      t55 = dlog(l1)
      t56 = t55*l1
      t57 = 12*t56
      t58 = 24*t55
      t59 = t55**2
      t63 = t55*t6
      t71 = -t40+t41+t52+t54-t57+t58+6*t59*l1-12*t59+24*t63+6*t17*l1-2*t
     #21*l1-12*t63*l1-24-t46
      tncoff9 = (-t3-6*t4+12+t7)*t9*t12/1536-(t3+24*t2*t6-6*t17-24*t19+2
     #*t21-12-t7)*t12/768-t17/256+t21/768-t28*t30*t33*t12/128-t17*t6/768
     #+(t39-t40+t41+12*t4*t6-3*t44-24-t46)*t9/6144+t21*t6/768-t71*t30*t3
     #2/1536+z3/384+(-t39-t52+6*t4*l1-t54+t57-t58+24+t46)*t30*t33/3072-(
     #-t3-6*l1-6*t53+6*t56-12*t55+12*t2*l1+12+t7)*t30*t32*t12/384-1.D0/1
     #28.D0-t6/128
      t1 = dlog(l3)
      t4 = dlog(1-l3)
      t6 = dlog(2.D0)
      t8 = 1-z
      t9 = dlog(t8)
      t12 = 1/l3
      t15 = -1/t8
      t24 = dlog(l1)
      t28 = 12*t1
      t29 = 12*t4
      t30 = 24*t6
      t37 = 1/l1*t12
      t42 = 0.3141592653589793D1**2
      t44 = t4**2
      t50 = t1**2
      t52 = t6**2
      t56 = -t28+t29+t30-4*t42+3*t44+24-6*t1*t4+12*t4*t6+3*t50+12*t52-12
     #*t1*t6
      tncoff10 = (6*t1-6*t4-12*t6+24*t9-12)*t12*t15/1536+(-6*t1*l1+6*t4*
     #l1+12*t6*l1+24*t24-12*t24*l1+t28-t29-t30-24+12*l1)/(l1-2)*t37/3072
     #-t37*t15/128+t56*t12/6144
      t1 = 2-l1
      t2 = dlog(t1)
      t3 = 24*t2
      t4 = dlog(l1)
      t5 = 24*t4
      t6 = dlog(2.D0)
      t7 = 48*t6
      t10 = 0.3141592653589793D1**2
      t12 = t6**2
      t14 = t4**2
      t16 = t2**2
      t22 = -t3-t5+t7-48*t2*t6-4*t10+48*t12+12*t14+12*t16-48*t4*t6+24*t4
     #*t2+24
      t24 = -1/t1
      t26 = 1/l1
      t30 = 1/l3
      t31 = -1+z
      t32 = 1/t31
      t36 = dlog(l3)
      t47 = dlog(-t31)
      tncoff11 = -t22*t24*t26/1536+t24*t26*t30*t32/64+(t3+12*t36+t5-t7-2
     #4)*t24*t26*t30/3072-(12*t2+12*t4-24*t6+24*t47-12)*t24*t26*t32/384
      t1 = l1-2
      t2 = 1/t1
      t3 = 1/l1
      t5 = 1/l3
      t11 = dlog(l3)
      t13 = dlog(l1)
      t15 = dlog(2.D0)
      t19 = dlog(4-4*l3)
      t22 = dlog(-t1)
      tncoff12 = -t2*t3*t5/(-1+z)/64+(-12*t11-24*t13+24*t15+12*t19-24*t2
     #2+24)*t2*t3*t5/3072
      t1 = dlog(z)
      t2 = 12*t1
      t3 = dlog(2.D0)
      t4 = 12*t3
      t5 = 1-z
      t6 = dlog(t5)
      t7 = 24*t6
      t8 = dlog(l3)
      t11 = 1/l3
      t14 = -1/t5
      t21 = t3**2
      t23 = t1**2
      t25 = 0.3141592653589793D1**2
      t29 = t6**2
      t31 = t2+t4-t7-24*t1*t6-24*t6*t3+6*t21+6*t23-2*t25+12*t3*t1+24*t29
     #+12
      t34 = 1/l1
      t38 = dlog(l1)
      tncoff13 = (-t2-t4+t7+6*t8-12)*t11*t14/1536-t31*t14/768-t34*t11*t1
     #4/128-(-6*t38*l1+6*t1*l1-t2-t4+t7+6*l1+12*t38-12*t6*l1-12+6*t3*l1)
     #/(l1-2)*t34*t14/384
      t3 = dlog(-z*(-1+l3))
      t5 = dlog(z)
      t7 = dlog(2.D0)
      t9 = dlog(l3)
      t11 = 1-z
      t12 = dlog(t11)
      t15 = 1/l3
      t18 = -1/t11
      t21 = 2-l1
      tncoff14 = (6*t3+6*t5+12*t7-6*t9-24*t12+12)*t15*t18/1536+1/l1*t15*
     #t18/128
      t1 = z*l1
      t6 = dlog(-(t1+2-l1)*(t1-2*z-l1))
      t8 = 1-z
      t9 = dlog(t8)
      t11 = dlog(z)
      t13 = 2-l1
      t14 = dlog(t13)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t22 = -1/t13
      t24 = 1/l1
      t26 = -1/t8
      tncoff15 = -(6*t6-24*t9+6*t11-12*t14+12*t16-12*t18+12)*t22*t24*t26
     #/384-t22*t24/l3*t26/64
      tncoff16 = 1/(l1-2)/l1/l3/(-1+z)/64
      t1 = 1/l3
      t2 = -1+z
      t3 = 1/t2
      t6 = dlog(2.D0)
      t9 = dlog(-t2)
      t15 = 0.3141592653589793D1**2
      t17 = 12*t6
      t18 = dlog(l3)
      t23 = 6*l1
      t26 = 1/(l1-2)
      t28 = 1/l1
      t38 = dlog(l1)
      t48 = t6**2
      tncoff17 = t3*t1/256-(-6*t6+12*t9-6)*t3/768-t6/256+t15/1536+(-t17+
     #6*t18-12)*t1/6144+(-t23+12)*t26*t28*t1/3072-(-3*l1+6)*t26*t28*t3/3
     #84-(-t17+t23+12*t38-6*t38*l1+6*t6*l1-12)*t26*t28/1536-t48/512-1.D0
     #/256.D0
      t1 = 1/l3
      t6 = dlog(l3)
      t8 = dlog(2.D0)
      t11 = dlog(1-l3)
      tncoff18 = -t1/(-1+z)/256+(12-6*t6+12*t8+6*t11)*t1/6144+(-12+6*l1)
     #/(l1-2)/l1*t1/3072
      t1 = l1-2
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t12 = dlog(l1)
      t15 = dlog(-t1)
      t17 = dlog(2.D0)
      tncoff19 = -t4/l3/256+t4/(-1+z)/64-(12-12*t12-12*t15+24*t17)*t2*t3
     #/1536
      tncoff20 = 1/(l1-2)/l1/l3/256
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff21 = -(6*t1+6*t3-12*t6+6)*t10/768-(-6+3*l1)/(l1-2)/l1*t10/38
     #4-1/l3*t10/256
      tncoff22 = 1/l3/(-1+z)/256
      tncoff23 = -1/(l1-2)/l1/(-1+z)/64
      t3 = dlog(2.D0)
      tncoff24 = -1/l3/1024-t3/512-(-6+3*l1)/(l1-2)/l1/1536-1.D0/512.D0+
     #1/(-1+z)/256
      tncoff25 = 1/l3/1024
      tncoff26 = -1/(l1-2)/l1/256
      tncoff27 = -1/(-1+z)/256
      tncoff28 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum3e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum3e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum3e1(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum3e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum3e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum3e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum3e2(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum3e2(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = RVnum3e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = RVnum3e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = RVnum3e2(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = RVnum3e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = RVnum3e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = RVnum3e3(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = RVnum3e3(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = RVnum3e4(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff28)

CCC Sum

      t15 = tcnum22*tncoff22+tcnum2*tncoff2+tcnum14*tncoff14+tcnum21*tnc
     #off21+tcnum27*tncoff27+tcnum24*tncoff24+tcnum28*tncoff28+tcnum11*t
     #ncoff11+tcnum26*tncoff26+tcnum5*tncoff5+tcnum19*tncoff19+tcnum25*t
     #ncoff25+tcnum10*tncoff10+tcnum18*tncoff18
      t30 = tcnum12*tncoff12+tcnum20*tncoff20+tcnum9*tncoff9+tcnum8*tnco
     #ff8+tcnum6*tncoff6+tcnum17*tncoff17+tcnum13*tncoff13+tcnum15*tncof
     #f15+tcnum7*tncoff7+tcnum23*tncoff23+tcnum1*tncoff1+tcnum4*tncoff4+
     #tcnum16*tncoff16+tcnum3*tncoff3
      RV4e0int = t15+t30

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV4e1int(var,ndim,Tcoeff)

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
      double precision RVnum3e0
      external RVnum3e0
      double precision RVnum3e1
      external RVnum3e1
      double precision RVnum3e2
      external RVnum3e2
      double precision RVnum3e3
      external RVnum3e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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
      t3 = dlog(2.D0)
      t6 = t3**2
      t12 = 1-z
      t13 = dlog(t12)
      t18 = t13**2
      t24 = -1/t12
      t27 = dlog(l3)
      t29 = 12*t6
      t30 = 4*t1
      t33 = t27**2
      t35 = 24*t3
      t37 = 1/l3
      t47 = t3*l1
      t49 = dlog(l1)
      t50 = t49*l1
      t53 = t49**2
      t57 = t49*t3
      t65 = -t29+t30+12*l1+12*t47-12*t50+24*t49+6*t53*l1-12*t53+24*t57+6
     #*t6*l1-2*t1*l1-12*t57*l1-24-t35
      t67 = 1/(l1-2)
      t69 = 1/l1
      t72 = 4*t3
      t74 = 4*t49
      t75 = 2*l1
      t76 = 2*t50
      t77 = 2*t47
      t89 = t69*t37
      tncoff1 = -1.D0/128.D0+t1/768+t1*t3/768-t6*t3/768-t6/256-t3/128+z3
     #/384-(24*t13+24*t13*t3-6*t6-24*t18+2*t1-12-12*t3)*t24/768+(12*t27-
     #t29+t30+12*t27*t3-3*t33-24-t35)*t37/6144+(-4*t13+2+2*t3-t27)*t37*t
     #24/256-t65*t67*t69/1536-(t72-8*t13-t74-t75+4+t76-t77+4*t13*l1)*t69
     #*t67*t24/128+(t72-2*t27-t74-t75+4+t76-t77+t27*l1)*t67*t89/512+t89*
     #t24/128
      t1 = dlog(l3)
      t4 = dlog(1-l3)
      t6 = dlog(2.D0)
      t8 = 0.3141592653589793D1**2
      t10 = t4**2
      t16 = t1**2
      t18 = t6**2
      t22 = -12*t1+12*t4+24*t6-4*t8+3*t10+24-6*t1*t4+12*t4*t6+3*t16+12*t
     #18-12*t1*t6
      t23 = 1/l3
      t26 = 1-z
      t27 = dlog(t26)
      t33 = -1/t26
      t37 = dlog(l1)
      t49 = l1-2
      t50 = 1/t49
      t53 = 1/l1*t23
      tncoff2 = t22*t23/6144+(4*t27-t4-2-2*t6+t1)*t23*t33/256+(-t1*l1-2*
     #t37*l1+t4*l1+2*t6*l1+2*t1-4*t6-2*t4+2*l1+4*t37-4)*t50*t53/512-t49*
     #t50*t53*t33/128
      t1 = 2-l1
      t2 = dlog(t1)
      t4 = dlog(l1)
      t6 = dlog(2.D0)
      t10 = 0.3141592653589793D1**2
      t12 = t6**2
      t14 = t4**2
      t16 = t2**2
      t22 = -24*t2-24*t4+48*t6-48*t2*t6-4*t10+48*t12+12*t14+12*t16-48*t4
     #*t6+24*t4*t2+24
      t24 = -1/t1
      t26 = 1/l1
      t29 = 4*t4
      t30 = 4*t2
      t31 = 8*t6
      t32 = 1-z
      t33 = dlog(t32)
      t38 = -1/t32
      t42 = dlog(l3)
      t46 = 1/l3
      tncoff3 = -t22*t24*t26/1536-(t29+t30-t31+8*t33-4)*t26*t24*t38/128+
     #(t29+t30-t31+2*t42-4)*t24*t26*t46/512+t24*t26*t46*t38/64
      t1 = 2-l1
      t2 = dlog(t1)
      t4 = dlog(l3)
      t8 = dlog(4-4*l3)
      t10 = dlog(2.D0)
      t12 = dlog(l1)
      t16 = -1/t1
      t18 = 1/l1
      t19 = 1/l3
      tncoff4 = (-4*t2-2*t4+2*t8+4*t10-4*t12+4)*t16*t18*t19/512-t16*t18*
     #t19/(-1+z)/64
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t12 = t3**2
      t14 = t1**2
      t16 = 0.3141592653589793D1**2
      t20 = t6**2
      t22 = 12*t1+12*t3-24*t6-24*t1*t6-24*t6*t3+6*t12+6*t14-2*t16+12*t1*
     #t3+24*t20+12
      t24 = -1/t5
      t29 = dlog(l3)
      t32 = 1/l3
      t36 = dlog(l1)
      t51 = 1/l1
      t53 = l1-2
      t54 = 1/t53
      tncoff5 = -t22*t24/768+(4*t6-2*t1-2+t29-2*t3)*t32*t24/256-(-2*t36*
     #l1-4*t6*l1+2*t3*l1+2*t1*l1-4-4*t1-4*t3+8*t6+4*t36+2*l1)*t51*t54*t2
     #4/128-t53*t54*t51*t32*t24/128
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t8 = dlog(-z*(-1+l3))
      t9 = dlog(z)
      t10 = dlog(l3)
      t12 = 1/l3
      t15 = -1/t1
      tncoff6 = (-4*t2+2+2*t4+t8+t9-t10)*t12*t15/256+1/l1*t12*t15/128
      t1 = dlog(l1)
      t3 = z*l1
      t8 = dlog(-(t3+2-l1)*(t3-2*z-l1))
      t10 = 1-z
      t11 = dlog(t10)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = 2-l1
      t18 = dlog(t17)
      t21 = 1/l1
      t24 = -1/t17
      t26 = -1/t10
      tncoff7 = -(-4*t1+2*t8-8*t11+4*t13+2*t15-4*t18+4)*t21*t24*t26/128-
     #t24*t21/l3*t26/64
      tncoff8 = 1/(l1-2)/l1/l3/(-1+z)/64
      t1 = 0.3141592653589793D1**2
      t3 = dlog(2.D0)
      t4 = t3**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t16 = 12*t3
      t17 = dlog(l3)
      t20 = 1/l3
      t26 = dlog(l1)
      t33 = l1-2
      t34 = 1/t33
      t36 = 1/l1
      t39 = -t33
      tncoff9 = -1.D0/256.D0+t1/1536-t4/512-t3/256-(-6*t3+12*t9-6)*t13/7
     #68+(-t16+6*t17-12)*t20/6144+t20*t13/256-(-t16+6*l1+12*t26-6*t26*l1
     #+6*t3*l1-12)*t34*t36/1536-t39*t36*t34*t13/128+t39*t34*t36*t20/512
      t1 = dlog(l3)
      t3 = dlog(2.D0)
      t6 = dlog(1-l3)
      t9 = 1/l3
      tncoff10 = (12-6*t1+12*t3+6*t6)*t9/6144-t9/(-1+z)/256+1/l1*t9/512
      t1 = l1-2
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t12 = dlog(l1)
      t15 = dlog(-t1)
      t17 = dlog(2.D0)
      tncoff11 = -t4/l3/256+t4/(-1+z)/64-(12-12*t12-12*t15+24*t17)*t2*t3
     #/1536
      tncoff12 = 1/(l1-2)/l1/l3/256
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff13 = -(6*t1+6*t3-12*t6+6)*t10/768-1/l3*t10/256-1/l1*t10/128
      tncoff14 = 1/l3/(-1+z)/256
      tncoff15 = -1/(l1-2)/l1/(-1+z)/64
      t3 = dlog(2.D0)
      tncoff16 = -1/l3/1024-t3/512-(-6+3*l1)/(l1-2)/l1/1536-1.D0/512.D0+
     #1/(-1+z)/256
      tncoff17 = 1/l3/1024
      tncoff18 = -1/(l1-2)/l1/256
      tncoff19 = -1/(-1+z)/256
      tncoff20 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum3e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum3e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum3e1(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum3e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum3e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum3e2(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum3e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum3e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      RV4e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV4e2int(var,ndim,Tcoeff)

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
      double precision RVnum3e0
      external RVnum3e0
      double precision RVnum3e1
      external RVnum3e1
      double precision RVnum3e2
      external RVnum3e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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
      t3 = dlog(2.D0)
      t4 = t3**2
      t7 = 2*t3
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t16 = dlog(l3)
      t18 = 1/l3
      t27 = dlog(l1)
      t32 = 1/l1
      t34 = l1-2
      t35 = 1/t34
      t38 = -t34
      tncoff1 = -1.D0/256.D0+t1/1536-t4/512-t3/256-(-2-t7+4*t9)*t13/256+
     #(-2-t7+t16)*t18/1024+t18*t13/256-(2*l1-4+2*t3*l1-4*t3-2*t27*l1+4*t
     #27)*t32*t35/512-t38*t32*t35*t13/128+t38*t35*t32*t18/512
      t1 = dlog(2.D0)
      t4 = dlog(1-l3)
      t5 = dlog(l3)
      t7 = 1/l3
      tncoff2 = (2+2*t1+t4-t5)*t7/1024-t7/(-1+z)/256+1/l1*t7/512
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t9 = 1/l1
      t12 = -1/t3
      t15 = t12*t9
      tncoff3 = -(8*t1-4*t4+4-4*t6)*t9*t12/512+t15/(-1+z)/64-t15/l3/256
      tncoff4 = 1/(l1-2)/l1/l3/256
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff5 = -(2+2*t1+2*t3-4*t6)*t10/256-1/l3*t10/256-1/l1*t10/128
      tncoff6 = 1/l3/(-1+z)/256
      tncoff7 = -1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff8 = -t1/512-1.D0/512.D0+1/(-1+z)/256-1/l3/1024-1/l1/512
      tncoff9 = 1/l3/1024
      tncoff10 = -1/(l1-2)/l1/256
      tncoff11 = -1/(-1+z)/256
      tncoff12 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV4e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV4e3int(var,ndim,Tcoeff)

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
      double precision RVnum3e0
      external RVnum3e0
      double precision RVnum3e1
      external RVnum3e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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
      tncoff1 = -t1/512-1.D0/512.D0+1/(-1+z)/256-1/l3/1024-1/l1/512
      tncoff2 = 1/l3/1024
      tncoff3 = -1/(l1-2)/l1/256
      tncoff4 = -1/(-1+z)/256
      tncoff5 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      RV4e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV4e4int(var,ndim,Tcoeff)

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
      double precision RVnum3e0
      external RVnum3e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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

      tncoff1 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV4e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

