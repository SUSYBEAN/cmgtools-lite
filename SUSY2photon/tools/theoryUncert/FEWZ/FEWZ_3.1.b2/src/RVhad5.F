CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV5e0int(var,ndim,Tcoeff)

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

      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(l3)
      t6 = dlog(2.D0)
      t7 = 24*t6
      t8 = 0.3141592653589793D1**2
      t9 = 4*t8
      t10 = t6**2
      t11 = 12*t10
      t12 = 48*t2
      t14 = t2**2
      t15 = 48*t14
      t16 = t2*t6
      t17 = 48*t16
      t18 = t3*t6
      t20 = t3**2
      t22 = -24*t2*t3-24-t7+t9-t11+t12+12*t3-t15+t17+12*t18-3*t20
      t23 = 1/l3
      t26 = -1/t1
      t29 = dlog(l1)
      t30 = t29*t6
      t31 = t30*l1
      t32 = 24*t31
      t33 = t3*t29
      t38 = 48*t6
      t39 = 8*t8
      t40 = 24*t10
      t41 = 24*t3
      t42 = 24*t18
      t43 = 6*t20
      t44 = t29**2
      t45 = 24*t44
      t46 = 48+t32-12*t33*l1+12*t18*l1+t38-t39+t40-t41-t42+t43+t45
      t47 = t6*l1
      t48 = 24*t47
      t49 = t29*l1
      t50 = 24*t49
      t51 = 48*t30
      t52 = t44*l1
      t53 = 12*t52
      t54 = t8*l1
      t55 = 4*t54
      t56 = t10*l1
      t57 = 12*t56
      t58 = 24*l1
      t59 = 48*t29
      t60 = t3*l1
      t65 = -t48+t50-t51-t53+t55-t57-t58-t59+12*t60+24*t33-3*t20*l1
      t68 = 1/(l1-2)
      t70 = 1/l1
      t77 = t2*t29
      t81 = 24+12*t31+24*t16*l1-24*t77*l1+t7-t9+t11-t12+t15-t17+12*t44
      t93 = t2*l1
      t95 = -12*t47+12*t49-24*t30-6*t52+2*t54-6*t56-24*t14*l1+48*t77-12*
     #l1-24*t29+24*t93
      t107 = 16*z3
      t108 = t10*t6
      t109 = 8*t108
      t114 = t6*t8
      t115 = 8*t114
      t116 = -48-t38+t39-t40+t41+t42+12*t3*t10-t43+t20*t3+t107-t109-4*t3
     #*t8-6*t20*t6+t115
      t124 = t29*t8
      t129 = t44*t6
      t132 = t44*t29
      t134 = -48+t53-t55+t57-4*t114*l1+4*t124*l1-8*z3*l1+12*t129*l1+t58+
     #t107+t39-t38-t40+8*t132
      t140 = t29*t10
      t145 = -24*t129-4*t132*l1+4*t108*l1+24*t140-12*t140*l1+t115-8*t124
     #-t109+t59-t32+t48-t50+t51-t45
      t164 = t8**2
      t166 = t10**2
      t179 = -24-t7+t9-t11+t12-t15+t17+8*z3-4*t108+4*t114+32*t14*t2-48*t
     #14*t6+24*t2*t10-8*t2*t8
      tncoff1 = -t22*t23*t26/1536+(t46+t65)*t68*t70*t23/3072-(t81+t95)*t
     #68*t70*t26/384-1.D0/64.D0-t6/64+t8/384-t10/128+t116*t23/6144+z3/19
     #2-t108/384+t114/384-(t134+t145)*t68*t70/1536-(4+4*t6-8*t2-2*t3-2*t
     #47+2*t49-2*l1-4*t29+4*t93+t60)*t70*t68*t23*t26/128-t164/10240-t166
     #/1536-t179*t26/768+t8*t10/768+t6*z3/192
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(1-l3)
      t7 = dlog(2.D0)
      t10 = dlog(l3)
      t14 = t2**2
      t16 = t4*t7
      t21 = t7**2
      t23 = 0.3141592653589793D1**2
      t25 = t10*t4
      t27 = t10*t7
      t29 = t4**2
      t31 = t10**2
      t33 = -24*t2*t4-48*t2*t7+24*t2*t10-48*t2+48*t14+12*t16+12*t4+24*t7
     #-12*t10+12*t21-4*t23-6*t25-12*t27+3*t29+3*t31+24
      t34 = 1/l3
      t37 = -1/t1
      t40 = dlog(l1)
      t41 = t40*t4
      t44 = t40*t7
      t51 = t10*l1
      t53 = t10*t40
      t59 = t40**2
      t65 = 24*t4
      t66 = 48*t7
      t67 = 24*t10
      t68 = -48-12*t41*l1-24*t44*l1+24*t41-4*t23*l1+48*t44-12*t51-24*t53
     #+3*t31*l1+12*t21*l1+12*t59*l1+12*t53*l1-24*t59-t65-t66+t67
      t69 = 12*t25
      t70 = 24*t27
      t71 = 24*t16
      t72 = t7*l1
      t74 = 24*t21
      t75 = 8*t23
      t76 = 6*t29
      t77 = t40*l1
      t79 = 6*t31
      t90 = t4*l1
      t92 = t69+t70-t71+24*t72-t74+t75-t76-24*t77-t79+3*t29*l1+24*l1-12*
     #t27*l1-6*t25*l1+12*t16*l1+48*t40+12*t90
      t95 = 1/(l1-2)
      t97 = 1/l1
      t121 = 48+12*t4*t21+6*t29*t7+t29*t4-8*t7*t23+4*t10*t23-12*t10*t21+
     #6*t31*t7+8*t21*t7-16*z3-4*t4*t23-3*t29*t10
      t127 = t65+t66-t67-t69-t70+t71+t74-t75+t76+t79-t31*t10-12*t25*t7+3
     #*t31*t4
      t141 = 8*t2-2*t4-4*t7+2*t10-4-2*t77+2*t72+t90-t51-4*t2*l1+4*t40+2*
     #l1
      tncoff2 = -t33*t34*t37/1536+(t68+t92)*t95*t97*t34/3072+(t121+t127)
     #*t34/6144-t141*t97*t95*t34*t37/128
      t1 = dlog(2.D0)
      t2 = 96*t1
      t3 = 0.3141592653589793D1**2
      t4 = 8*t3
      t5 = t1**2
      t6 = 96*t5
      t7 = dlog(l1)
      t8 = 2-l1
      t9 = dlog(t8)
      t10 = t7*t9
      t15 = 48*t7
      t16 = 48*t9
      t17 = t7**2
      t22 = t9**2
      t27 = 48+t2-t4+t6+96*t10*t1+64*t5*t1-t15-t16+48*t17*t1-96*t7*t5-24
     #*t7*t22-24*t17*t9
      t32 = t9*t1
      t33 = 96*t32
      t34 = 48*t10
      t41 = t7*t1
      t42 = 96*t41
      t43 = 24*t17
      t45 = 24*t22
      t50 = -96*t9*t5+48*t1*t22-t33+t34+8*t3*t9+8*t7*t3-16*t1*t3-t42+t43
     #-16*z3+t45-8*t22*t9-8*t17*t7
      t53 = -1/t8
      t55 = 1/l1
      t66 = 1-z
      t67 = dlog(t66)
      t75 = t67**2
      t79 = -24-12*t17-12*t22+4*t3-48*t5+48*t41-48*t1+24*t9+24*t7-48*t67
     #*t9+48*t32-24*t10-48*t67*t7+48*t67-48*t75+96*t67*t1
      t82 = -1/t66
      t86 = dlog(l3)
      t94 = t86**2
      t96 = -48-t43-t45+t4-t6+t42-t2+t16+t15+t33-t34+24*t86-24*t86*t7+48
     #*t86*t1-24*t86*t9-6*t94
      t98 = 1/l3
      tncoff3 = -(t27+t50)*t53*t55/1536-t79*t53*t55*t82/384+t96*t53*t55*
     #t98/3072-(-4-8*t1+4*t9+4*t7+8*t67+2*t86)*t55*t53*t98*t82/128
      t1 = 1-l3
      t2 = dlog(t1)
      t4 = dlog(l3)
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t13 = 1-z
      t14 = dlog(t13)
      t17 = 1/l1
      t20 = -1/t6
      t21 = 1/l3
      t32 = t4**2
      t38 = t7**2
      t41 = dlog(4*t1)
      t42 = t41**2
      t44 = t9**2
      t46 = 48-24*t4-48*t7+48*t9-48*t11+6*t32+24*t4*t11+24*t4*t7+24*t38+
     #6*t42+24*t44
      t47 = 0.3141592653589793D1**2
      t49 = t11**2
      t68 = -8*t47+24*t49+24*t41-12*t4*t41-24*t7*t41-48*t7*t9+48*t11*t7-
     #48*t11*t9+24*t9*t41-24*t4*t9-24*t11*t41
      tncoff4 = (4+2*t2-2*t4-4*t7+8*t9-4*t11-8*t14)*t17*t20*t21/t13/128+
     #(t46+t68)*t20*t17*t21/3072
      t1 = dlog(l3)
      t3 = dlog(2.D0)
      t5 = dlog(z)
      t7 = 1-z
      t8 = dlog(t7)
      t10 = dlog(l1)
      t13 = t8*l1
      t15 = t5*l1
      t17 = t3*l1
      t19 = t10*l1
      t22 = 2*t1-4*t3-4*t5+8*t8-4+4*t10+2*l1-4*t13+2*t15+2*t17-2*t19-t1*
     #l1
      t23 = 1/l1
      t26 = 1/(l1-2)
      t27 = 1/l3
      t30 = -1/t7
      t34 = t8**2
      t37 = t5**2
      t41 = t3**2
      t44 = 0.3141592653589793D1**2
      t47 = 24*t3
      t50 = 24*t5
      t51 = 48*t8
      t52 = 12*t37
      t55 = 24-32*t34*t8+4*t37*t5-8*z3-24*t8*t41+8*t8*t44+t47+4*t41*t3+t
     #50-t51+t52-4*t5*t44
      t56 = t3*t8
      t59 = 4*t44
      t60 = 12*t41
      t63 = t3*t5
      t64 = 24*t63
      t71 = t8*t5
      t72 = 48*t71
      t77 = 48*t56
      t78 = 48*t34
      t79 = -48*t56*t5-t59+t60+12*t5*t41+t64+12*t37*t3-24*t8*t37+48*t34*
     #t5-t72+48*t34*t3-4*t3*t44-t77+t78
      t88 = t1**2
      t92 = 24*t8*t1-12*t1+t47-t72-t77-12*t5*t1+t64+t50-t51+t78+t60-t59+
     #3*t88+t52-12*t1*t3+24
      t96 = t3*t10
      t98 = t10**2
      t112 = -24+24*t96-12*t98+6*t41*l1+6*t98*l1+24*t34*l1-24*t13+12*t15
     #+12*t17-12*t19-t47-t50+t51-t52+24*t10+12*l1
      t119 = t10*t5
      t122 = t8*t10
      t133 = -24*t56*l1-24*t71*l1+12*t63*l1-12*t119*l1+24*t122*l1-12*t96
     #*l1+t59-t60-t64+t72+t77+6*t37*l1+24*t119-2*t44*l1-48*t122-t78
      tncoff5 = -t22*t23*t26*t27*t30/128-(t55+t79)*t30/768-t92*t27*t30/1
     #536-(t112+t133)*t26*t23*t30/384
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t6 = dlog(l3)
      t9 = dlog(-z*(-1+l3))
      t14 = dlog(2.D0)
      t25 = t2**2
      t27 = 24*t2*t3+6*t6*t9+24*t2*t9+48*t14*t2-24*t2*t6+12*t14*t6-12*t9
     #-24*t14-12*t3+12*t6-48*t25
      t28 = t9**2
      t30 = t6**2
      t32 = t3**2
      t42 = 0.3141592653589793D1**2
      t44 = t14**2
      t47 = -24-3*t28-3*t30-3*t32+6*t3*t6-12*t14*t3-12*t9*t14-6*t3*t9+4*
     #t42-12*t44+48*t2
      t49 = 1/l3
      t52 = -1/t1
      t63 = dlog(l1)
      t72 = 2*t9+4*t14+2*t3-2*t6-8*t2+t6*l1-2*t14*l1-4*t63-2*l1-t3*l1+2*
     #t63*l1+4-t9*l1+4*t2*l1
      tncoff6 = -(t27+t47)*t49*t52/1536-t72/l1/(l1-2)*t49*t52/128
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2-l1
      t5 = dlog(t4)
      t7 = dlog(z)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t13 = z*l1
      t18 = dlog(-(t13+2-l1)*(t13-2*z-l1))
      t20 = dlog(l3)
      t23 = 1/l1
      t26 = -1/t4
      t30 = -1/t1
      t60 = 24+12*t9*t7-12*t11*t7-24*t2*t7+6*t7*t18+12*t9*t18-12*t11*t18
     #-12*t7*t5-24*t5*t9-12*t5*t18-24*t2*t18-48*t2*t9+24*t11*t5-24*t11*t
     #9
      t65 = t18**2
      t67 = t11**2
      t69 = t5**2
      t72 = t2**2
      t79 = t7**2
      t81 = 0.3141592653589793D1**2
      t83 = t9**2
      t85 = 48*t2*t11+48*t2*t5+3*t65+12*t67+12*t69-48*t2+48*t72-24*t5+12
     #*t7+24*t9-24*t11+12*t18+3*t79-4*t81+12*t83
      tncoff7 = -(-8*t2-4*t5+2*t7+4*t9-4*t11+2*t18-2*t20+4)*t23*t26/l3*t
     #30/128-(t60+t85)*t26*t23*t30/384
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = dlog(l3)
      t7 = dlog(l1)
      t9 = 2-l1
      t10 = dlog(t9)
      t18 = z**2
      t21 = l1**2
      t26 = dlog(4*z-4*z*l3+2*l1-4*z*l1+2*t18*l1-t21+2*z*t21-t18*t21)
      t28 = 1-z
      t29 = dlog(t28)
      tncoff8 = -(-2*t1-4*t3+2*t5+4*t7+4*t10-2*t26+8*t29-4)/l1/t9/l3/t28
     #/128
      t1 = dlog(l3)
      t2 = 12*t1
      t3 = t1**2
      t5 = dlog(2.D0)
      t8 = t5**2
      t9 = 12*t8
      t10 = 0.3141592653589793D1**2
      t11 = 4*t10
      t12 = 24*t5
      t14 = 1/l3
      t19 = dlog(l1)
      t20 = t19**2
      t22 = 24*t19
      t23 = t19*t5
      t27 = t19*l1
      t28 = 12*t27
      t33 = t5*l1
      t34 = 12*t33
      t35 = 12*l1
      t38 = -t9+t11-12*t20+t22+24*t23+6*t20*l1-t28-2*t10*l1+6*t8*l1+t34+
     #t35-12*t23*l1-24-t12
      t39 = l1-2
      t40 = 1/t39
      t42 = 1/l1
      t48 = -1+z
      t49 = 1/t48
      t59 = dlog(-t48)
      t60 = t59**2
      t64 = 12*t5
      t65 = 24*t59
      tncoff9 = (t2-3*t3+12*t1*t5-t9+t11-24-t12)*t14/6144-t8/256+t10/768
     #-t38*t40*t42/1536+t39*t42*t40*t14*t49/128-t8*t5/768+z3/384-(-6*t8+
     #2*t10-24*t60+24*t59*t5-12-t64+t65)*t49/768+t5*t10/768-(-12*t19+6*t
     #27-6*t33-6*l1+12*t59*l1+12+t64-t65)*t40*t42*t49/384+(-t2-t22+t28-t
     #34-t35+6*t1*l1+24+t12)*t40*t42*t14/3072-1.D0/128.D0-(6*t1-12-t64+t
     #65)*t14*t49/1536-t5/128
      t1 = dlog(2.D0)
      t4 = dlog(1-l3)
      t6 = dlog(l3)
      t8 = 1-z
      t9 = dlog(t8)
      t12 = 1/l3
      t15 = -1/t8
      t18 = 24*t1
      t21 = 12*t4
      t22 = 12*t6
      t23 = t6**2
      t25 = t1**2
      t27 = 0.3141592653589793D1**2
      t29 = t4**2
      t35 = t18-6*t6*t4+t21-t22+3*t23+12*t25-4*t27+3*t29+12*t4*t1-12*t6*
     #t1+24
      t40 = dlog(l1)
      t54 = 1/l1*t12
      tncoff10 = -(12*t1+6*t4-6*t6+12-24*t9)*t12*t15/1536+t35*t12/6144+(
     #-6*t6*l1+24*t40-12*t40*l1+12*t1*l1+6*t4*l1-t18-t21+t22-24+12*l1)/(
     #l1-2)*t54/3072-t54*t15/128
      t1 = dlog(l1)
      t2 = 24*t1
      t3 = t1**2
      t5 = 2-l1
      t6 = dlog(t5)
      t7 = 24*t6
      t8 = dlog(2.D0)
      t9 = 48*t8
      t10 = t6**2
      t12 = t8**2
      t20 = 0.3141592653589793D1**2
      t22 = -t2+12*t3-t7+t9+12*t10+48*t12-48*t1*t8+24*t6*t1-48*t6*t8-4*t
     #20+24
      t24 = -1/t5
      t26 = 1/l1
      t30 = 1/l3
      t31 = -1+z
      t32 = 1/t31
      t40 = dlog(-t31)
      t47 = dlog(l3)
      tncoff11 = -t22*t24*t26/1536+t26*t24*t30*t32/64-(12*t1+12*t6-24*t8
     #+24*t40-12)*t24*t26*t32/384+(t2+t7-t9-24+12*t47)*t24*t26*t30/3072
      t1 = 1/l1
      t2 = l1-2
      t3 = 1/t2
      t5 = 1/l3
      t11 = dlog(l1)
      t13 = dlog(2.D0)
      t16 = dlog(-t2)
      t20 = dlog(4-4*l3)
      t22 = dlog(l3)
      tncoff12 = -t1*t3*t5/(-1+z)/64+(-24*t11+24*t13-24*t16+12*t20-12*t2
     #2+24)*t3*t1*t5/3072
      t1 = 1-z
      t2 = dlog(t1)
      t3 = 24*t2
      t4 = dlog(l3)
      t6 = dlog(z)
      t7 = 12*t6
      t8 = dlog(2.D0)
      t9 = 12*t8
      t11 = 1/l3
      t14 = -1/t1
      t24 = dlog(l1)
      t32 = 1/l1
      t43 = t6**2
      t45 = t8**2
      t47 = 0.3141592653589793D1**2
      t49 = t2**2
      t53 = -24*t2*t6-t3+t7+t9+12+12*t6*t8+6*t43+6*t45-2*t47+24*t49-24*t
     #2*t8
      tncoff13 = -(-t3-6*t4+t7+t9+12)*t11*t14/1536-(t3-t7-t9-12*t2*l1+6*
     #l1-12+6*t6*l1+6*t8*l1-6*t24*l1+12*t24)/(l1-2)*t32*t14/384-t32*t11*
     #t14/128-t53*t14/768
      t1 = 1-z
      t2 = dlog(t1)
      t6 = dlog(-z*(-1+l3))
      t8 = dlog(2.D0)
      t10 = dlog(z)
      t12 = dlog(l3)
      t15 = 1/l3
      t18 = -1/t1
      t21 = 2-l1
      tncoff14 = -(24*t2-6*t6-12*t8-6*t10+6*t12-12)*t15*t18/1536+1/l1*t1
     #5*t18/128
      t1 = z*l1
      t6 = dlog(-(t1+2-l1)*(t1-2*z-l1))
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t15 = 1-z
      t16 = dlog(t15)
      t18 = dlog(l1)
      t22 = -1/t8
      t24 = 1/l1
      t26 = -1/t15
      tncoff15 = -(6*t6-12*t9+12*t11+6*t13-24*t16-12*t18+12)*t22*t24*t26
     #/384-t24*t22/l3*t26/64
      tncoff16 = 1/l1/(l1-2)/l3/(-1+z)/64
      t1 = 1/l3
      t2 = -1+z
      t3 = 1/t2
      t6 = dlog(2.D0)
      t7 = 12*t6
      t8 = dlog(l3)
      t17 = 1/(l1-2)
      t19 = 1/l1
      t23 = dlog(l1)
      t25 = 6*l1
      t39 = 0.3141592653589793D1**2
      t41 = t6**2
      t45 = dlog(-t2)
      tncoff17 = t1*t3/256-1.D0/256.D0+(-12-t7+6*t8)*t1/6144-t6/256-(6-3
     #*l1)*t17*t19*t3/384-(-12-t7+12*t23+t25-6*t23*l1+6*t6*l1)*t17*t19/1
     #536+(12-t25)*t17*t19*t1/3072+t39/1536-t41/512-(-6-6*t6+12*t45)*t3/
     #768
      t1 = 1/l3
      t6 = dlog(2.D0)
      t9 = dlog(1-l3)
      t11 = dlog(l3)
      tncoff18 = -t1/(-1+z)/256+(12*t6+6*t9-6*t11+12)*t1/6144+(-12+6*l1)
     #/(l1-2)/l1*t1/3072
      t1 = l1-2
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t13 = dlog(-t1)
      t15 = dlog(2.D0)
      t17 = dlog(l1)
      tncoff19 = t4/(-1+z)/64-t4/l3/256-(-12*t13+24*t15-12*t17+12)*t2*t3
     #/1536
      tncoff20 = 1/(l1-2)/l1/l3/256
      t2 = -1+z
      t3 = 1/t2
      t15 = dlog(z)
      t17 = dlog(2.D0)
      t20 = dlog(-t2)
      tncoff21 = -1/l3*t3/256-(-6+3*l1)/(l1-2)/l1*t3/384-(6+6*t15+6*t17-
     #12*t20)*t3/768
      tncoff22 = 1/l3/(-1+z)/256
      tncoff23 = -1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff24 = -1.D0/512.D0-t1/512+1/(-1+z)/256-1/l3/1024-(-6+3*l1)/(l
     #1-2)/l1/1536
      tncoff25 = 1/l3/1024
      tncoff26 = -1/(l1-2)/l1/256
      tncoff27 = -1/(-1+z)/256
      tncoff28 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum3e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum3e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum3e1(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum3e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum3e2(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum3e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum3e2(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum3e2(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = RVnum3e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = RVnum3e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = RVnum3e2(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = RVnum3e3(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = RVnum3e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = RVnum3e3(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = RVnum3e3(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = RVnum3e4(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff28)

CCC Sum

      t15 = tcnum8*tncoff8+tcnum16*tncoff16+tcnum20*tncoff20+tcnum22*tnc
     #off22+tcnum19*tncoff19+tcnum4*tncoff4+tcnum21*tncoff21+tcnum7*tnco
     #ff7+tcnum6*tncoff6+tcnum5*tncoff5+tcnum15*tncoff15+tcnum12*tncoff1
     #2+tcnum26*tncoff26+tcnum1*tncoff1
      t30 = tcnum2*tncoff2+tcnum3*tncoff3+tcnum24*tncoff24+tcnum25*tncof
     #f25+tcnum23*tncoff23+tcnum10*tncoff10+tcnum27*tncoff27+tcnum28*tnc
     #off28+tcnum18*tncoff18+tcnum17*tncoff17+tcnum11*tncoff11+tcnum13*t
     #ncoff13+tcnum14*tncoff14+tcnum9*tncoff9
      RV5e0int = t15+t30

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV5e1int(var,ndim,Tcoeff)

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
      t14 = 1-z
      t15 = dlog(t14)
      t16 = t15**2
      t24 = -1/t14
      t27 = dlog(l3)
      t29 = t27**2
      t33 = 12*t6
      t34 = 4*t1
      t35 = 24*t3
      t37 = 1/l3
      t46 = dlog(l1)
      t47 = t46**2
      t50 = t46*t3
      t54 = t46*l1
      t60 = t3*l1
      t65 = -t33+t34-12*t47+24*t46+24*t50+6*t47*l1-12*t54-2*t1*l1+6*t6*l
     #1+12*t60+12*l1-12*t50*l1-24-t35
      t66 = l1-2
      t67 = 1/t66
      t69 = 1/l1
      t72 = 2*t60
      t73 = 2*t54
      t74 = 4*t3
      t78 = 2*l1
      t79 = 4*t46
      tncoff1 = -1.D0/128.D0+t1/768+t3*t1/768-t6*t3/768-t6/256-t3/128+z3
     #/384-(-6*t6+2*t1-24*t16+24*t15*t3-12-12*t3+24*t15)*t24/768+(12*t27
     #-3*t29+12*t27*t3-t33+t34-24-t35)*t37/6144-(4*t15-2*t3-2+t27)*t37*t
     #24/256-t65*t67*t69/1536-(4-t72+t73+t74-8*t15+4*t15*l1-t78-t79)*t67
     #*t69*t24/128+(4-t72+t73+t74-2*t27+t27*l1-t78-t79)*t67*t69*t37/512+
     #t66*t69*t67*t37*t24/128
      t1 = dlog(2.D0)
      t3 = dlog(l3)
      t5 = dlog(1-l3)
      t10 = t3**2
      t12 = t1**2
      t14 = 0.3141592653589793D1**2
      t16 = t5**2
      t22 = 24*t1-6*t3*t5+12*t5-12*t3+3*t10+12*t12-4*t14+3*t16+12*t5*t1-
     #12*t3*t1+24
      t23 = 1/l3
      t26 = 1-z
      t27 = dlog(t26)
      t33 = -1/t26
      t39 = dlog(l1)
      t53 = 1/l1*t23
      tncoff2 = t22*t23/6144-(-4*t27+t5+2+2*t1-t3)*t23*t33/256+(-4-2*t5+
     #2*t3-4*t1-2*t39*l1+2*t1*l1+t5*l1-t3*l1+2*l1+4*t39)/(l1-2)*t53/512-
     #t53*t33/128
      t1 = dlog(l1)
      t3 = t1**2
      t5 = 2-l1
      t6 = dlog(t5)
      t8 = dlog(2.D0)
      t10 = t6**2
      t12 = t8**2
      t20 = 0.3141592653589793D1**2
      t22 = -24*t1+12*t3-24*t6+48*t8+12*t10+48*t12-48*t1*t8+24*t1*t6-48*
     #t6*t8-4*t20+24
      t24 = -1/t5
      t26 = 1/l1
      t29 = 1-z
      t30 = dlog(t29)
      t32 = 8*t8
      t33 = 4*t1
      t34 = 4*t6
      t38 = -1/t29
      t42 = dlog(l3)
      t46 = 1/l3
      tncoff3 = -t22*t24*t26/1536-(-4+8*t30-t32+t33+t34)*t24*t26*t38/128
     #+(-4+2*t42-t32+t33+t34)*t24*t26*t46/512+t26*t24*t46*t38/64
      t1 = dlog(l3)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(2.D0)
      t8 = dlog(l1)
      t12 = dlog(4-4*l3)
      t16 = -1/t3
      t18 = 1/l1
      t19 = 1/l3
      tncoff4 = (-2*t1-4*t4+4*t6-4*t8+2*t12+4)*t16*t18*t19/512-t18*t16*t
     #19/(-1+z)/64
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t8 = dlog(2.D0)
      t12 = t3**2
      t14 = t8**2
      t16 = 0.3141592653589793D1**2
      t18 = t2**2
      t22 = -24*t2*t3-24*t2+12*t3+12*t8+12+12*t8*t3+6*t12+6*t14-2*t16+24
     #*t18-24*t2*t8
      t24 = -1/t1
      t30 = dlog(l3)
      t32 = 1/l3
      t39 = dlog(l1)
      t54 = 1/l1
      tncoff5 = -t22*t24/768-(-4*t2+2*t8+2*t3-t30+2)*t32*t24/256-(-4*t8-
     #4*t3+8*t2+4*t39-4*t2*l1+2*t3*l1+2*t8*l1-2*t39*l1+2*l1-4)/(l1-2)*t5
     #4*t24/128-t54*t32*t24/128
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t8 = dlog(-z*(-1+l3))
      t9 = dlog(l3)
      t10 = dlog(z)
      t12 = 1/l3
      t15 = -1/t1
      t18 = 2-l1
      tncoff6 = -(4*t2-2*t4-t8-2+t9-t10)*t12*t15/256+1/l1*t12*t15/128
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 2-l1
      t5 = dlog(t4)
      t7 = dlog(z)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      t13 = z*l1
      t18 = dlog(-(t13+2-l1)*(t13-2*z-l1))
      t22 = -1/t4
      t24 = 1/l1
      t26 = -1/t1
      tncoff7 = -(4-8*t2-4*t5+2*t7+4*t9-4*t11+2*t18)*t22*t24*t26/128-t24
     #*t22/l3*t26/64
      tncoff8 = 1/l1/(l1-2)/l3/(-1+z)/64
      t1 = 0.3141592653589793D1**2
      t3 = dlog(2.D0)
      t4 = t3**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t16 = 12*t3
      t17 = dlog(l3)
      t20 = 1/l3
      t25 = dlog(l1)
      t33 = l1-2
      t34 = 1/t33
      t36 = 1/l1
      t40 = -t33*t34
      tncoff9 = -1.D0/256.D0+t1/1536-t4/512-t3/256-(-6-6*t3+12*t9)*t13/7
     #68+(-12-t16+6*t17)*t20/6144+t20*t13/256-(-12-t16+12*t25+6*l1-6*t25
     #*l1+6*t3*l1)*t34*t36/1536-t40*t36*t13/128+t40*t36*t20/512
      t1 = dlog(2.D0)
      t4 = dlog(1-l3)
      t6 = dlog(l3)
      t9 = 1/l3
      tncoff10 = (12*t1+6*t4-6*t6+12)*t9/6144-t9/(-1+z)/256+1/l1*t9/512
      t1 = l1-2
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t13 = dlog(-t1)
      t15 = dlog(2.D0)
      t17 = dlog(l1)
      tncoff11 = t4/(-1+z)/64-t4/l3/256-(-12*t13+24*t15-12*t17+12)*t2*t3
     #/1536
      tncoff12 = 1/(l1-2)/l1/l3/256
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff13 = -(6+6*t1+6*t3-12*t6)*t10/768-1/l3*t10/256-1/l1*t10/128
      tncoff14 = 1/l3/(-1+z)/256
      tncoff15 = -1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff16 = -1.D0/512.D0-t1/512+1/(-1+z)/256-1/l3/1024-(-6+3*l1)/(l
     #1-2)/l1/1536
      tncoff17 = 1/l3/1024
      tncoff18 = -1/(l1-2)/l1/256
      tncoff19 = -1/(-1+z)/256
      tncoff20 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum3e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum3e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum3e1(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum3e2(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum3e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum3e2(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum3e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum3e3(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      RV5e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV5e2int(var,ndim,Tcoeff)

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
      t32 = l1-2
      t33 = 1/t32
      t35 = 1/l1
      t39 = -t32*t33
      tncoff1 = -1.D0/256.D0+t1/1536-t4/512-t3/256-(-t7-2+4*t9)*t13/256+
     #(-t7-2+t16)*t18/1024+t18*t13/256-(2*l1-4+2*t3*l1-4*t3-2*t27*l1+4*t
     #27)*t33*t35/512-t39*t35*t13/128+t39*t35*t18/512
      t2 = dlog(1-l3)
      t3 = dlog(2.D0)
      t5 = dlog(l3)
      t7 = 1/l3
      tncoff2 = (t2+2+2*t3-t5)*t7/1024-t7/(-1+z)/256+1/l1*t7/512
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t10 = -1/t3
      t12 = 1/l1
      t15 = t10*t12
      tncoff3 = -(8*t1-4*t4+4-4*t6)*t10*t12/512+t15/(-1+z)/64-t15/l3/256
      tncoff4 = 1/(l1-2)/l1/l3/256
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff5 = -(2*t1+2+2*t3-4*t6)*t10/256-1/l3*t10/256-1/l1*t10/128
      tncoff6 = 1/l3/(-1+z)/256
      tncoff7 = -1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff8 = -1.D0/512.D0-t1/512+1/(-1+z)/256-1/l3/1024-1/l1/512
      tncoff9 = 1/l3/1024
      tncoff10 = -1/(l1-2)/l1/256
      tncoff11 = -1/(-1+z)/256
      tncoff12 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum3e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum3e1(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum3e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum3e2(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV5e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV5e3int(var,ndim,Tcoeff)

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
      tncoff1 = -1.D0/512.D0-t1/512+1/(-1+z)/256-1/l3/1024-1/l1/512
      tncoff2 = 1/l3/1024
      tncoff3 = -1/(l1-2)/l1/256
      tncoff4 = -1/(-1+z)/256
      tncoff5 = -1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum3e0(var,ndim,1.D0,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum3e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum3e1(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      RV5e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV5e4int(var,ndim,Tcoeff)

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

      tcnum1 = RVnum3e0(var,ndim,1.D0,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV5e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

