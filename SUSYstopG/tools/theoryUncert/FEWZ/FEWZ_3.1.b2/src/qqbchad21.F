CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbc21e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum21e0
      external qqbcnum21e0
      double precision qqbcnum21e1
      external qqbcnum21e1
      double precision qqbcnum21e2
      external qqbcnum21e2
      double precision qqbcnum21e3
      external qqbcnum21e3
      double precision qqbcnum21e4
      external qqbcnum21e4
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t4 = -1+l2
      t5 = 1/t4
      t7 = 1/l2
      t10 = 1-z
      t11 = dlog(t10)
      t12 = 0.3141592653589793D1**2
      t15 = -1/t10
      t21 = 1/l1
      t23 = -1+l1
      t24 = 1/t23
      t27 = -t23
      t31 = dlog(l1)
      t35 = t21*t15
      t36 = t11**2
      t40 = -1-l1*l2+l1+l2
      t43 = t5*t7
      t48 = -t4
      t55 = t31**2
      t59 = -t40
      t61 = t59*t7*t21
      t62 = t24*t5
      t63 = dlog(l2)
      t68 = t12*l2
      t69 = -t12+t68
      t78 = t7*t15
      t81 = (-4*z3+4*z3*l2)*t5*t7/1152+t11*t12*t15/576-(4*z3-4*z3*l1)*t2
     #1*t24/2304-t11*t27*t21*t24*t15*t31/48+t35*t36/48+t40*t24*t21*t43*t
     #15*t31/48-t48*t7*t5*t15*t36/24-t27*t24*t35*t55/192-t61*t62*t63*t31
     #/48-t69*t5*t7*t63/576-z3*t15/288-t69*t5*t78/576
      t82 = t12*l1
      t83 = -t12+t82
      t89 = t59*t24*t21
      t94 = t15*t11
      t101 = t48*t5
      t102 = t63**2
      t127 = t12**2
      t136 = -t83*t21*t24*t31/2304-t89*t43*t15*t63/24-t89*t43*t94/24-t61
     #*t62*t55/192-t101*t78*t102/24+t40*t7*t21*t62*t102/48-t101*t7*t94*t
     #63/12-t83*t24*t35/1152-t36*t11*t15/72+t21*t55*t31/1152+t7*t102*t63
     #/72-t127/276480+(-t68+t82*l2+t12-t82)*t7*t21*t24*t5/1152
      tncoff1 = t81+t136
      t1 = 0.3141592653589793D1**2
      t6 = 1/l2
      t8 = 1/l1
      t10 = 1/(-1+l1)
      t13 = 1/(-1+l2)
      t18 = l1*l2-l2
      t21 = t13*t6
      t22 = z-1
      t23 = 1/t22
      t24 = dlog(l1)
      t31 = dlog(-t22)
      t32 = t31**2
      t35 = -t18
      t39 = t24**2
      tncoff2 = (-t1*l1*l2+t1*l2)*t6*t8*t10*t13/1152+t18*t10*t8*t21*t23*
     #t24/48-t13*t23*t32/24-t35*t6*t8*t10*t13*t39/192-z3*t13/288+t1*t13*
     #t23/576-t35*t10*t8*t21*t23*t31/24
      t1 = 0.3141592653589793D1**2
      t3 = 1/(-1+l2)
      t4 = t1*t3
      t5 = 1/l2
      t6 = dlog(l2)
      t10 = 1/l1
      t11 = t5*t10
      t12 = dlog(l1)
      t13 = t12**2
      t17 = z-1
      t18 = 1/t17
      t19 = t5*t18
      t29 = t3*t5
      t30 = t6**2
      t34 = 1-l1
      t36 = -1/t34
      t44 = dlog(-t17)
      t55 = t10*t3
      t72 = t44**2
      tncoff3 = -t4*t5*t6/576-t11*t3*t13/192-t4*t19/576+z3*t3*t5/288-t11
     #*t3*t6*t12/48+t29*t18*t30/24+t34*t36*t10*t29*t18*t12/48+t29*t18*t4
     #4*t6/12+t34*t5*t10*t36*t3*t30/48-t55*t19*t44/24-t55*t19*t6/24+(-t1
     #+t1*l1)*t5*t10*t36*t3/1152+t29*t30*t6/72+t29*t18*t72/24
      t1 = -l2+1
      t3 = -1/t1
      t4 = t1*t3
      t5 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t10 = dlog(-t7)
      t11 = t8*t10
      t12 = dlog(l2)
      t17 = -l1+l1*l2
      t19 = 1/l1
      t20 = t17*t5*t19
      t21 = -1+l1
      t22 = 1/t21
      t23 = t22*t3
      t24 = t12**2
      t29 = dlog(-t21)
      t34 = t10**2
      t40 = t29**2
      t44 = 0.3141592653589793D1**2
      t47 = (-t44*l2+t44)*t3
      t53 = -t17*t22*t19
      t54 = t3*t5
      t59 = t44**2
      t61 = t5*t8
      t67 = t22*t8
      t70 = -t4*t5*t11*t12/12+t20*t23*t24/48+t20*t23*t29*t12/48-t34*t10*
     #t8/72-z3*t8/288+t20*t23*t40/192+t47*t5*t12/576-t53*t54*t8*t29/48-t
     #59/276480+t47*t61/576-t53*t54*t11/24-t67*t34/48
      t71 = t44*t22
      t110 = t44*l1
      t118 = t71*t8/1152-t67*t10*t29/48+(-4*z3+4*z3*l2)*t3*t5/1152+t5*t2
     #4*t12/72-t22*t40*t29/1152-t1*t5*t3*t8*t34/24-z3*t22/576-t53*t54*t8
     #*t12/24+t71*t29/2304-t67*t40/192+t10*t44*t8/576-t4*t61*t24/24+(t11
     #0-t110*l2)*t5*t19*t22*t3/1152
      tncoff4 = t70+t118
      t1 = 0.3141592653589793D1**2
      t3 = 1/(-1+l2)
      t5 = z-1
      t6 = 1/t5
      t11 = dlog(-t5)
      t12 = t11**2
      t15 = -1+l1
      t16 = 1/t15
      t17 = t16*t3
      t19 = dlog(-t15)
      t20 = t19**2
      tncoff5 = t1*t3*t6/576-t3*t6*t12/24-t17*t20/192-t17*t6*t11/24+t1*t
     #16*t3/1152-z3*t3/288-t17*t6*t19/48
      t2 = 1/(-1+l2)
      t4 = 1/l2
      t7 = -1+l1
      t8 = 1/t7
      t9 = t8*t2
      t10 = z-1
      t11 = 1/t10
      t12 = t4*t11
      t14 = dlog(-t7)
      t18 = t2*t4
      t19 = dlog(l2)
      t20 = t19**2
      t24 = t4*t8
      t30 = dlog(-t10)
      t35 = t30**2
      t39 = 0.3141592653589793D1**2
      t40 = t39*t2
      t46 = t14**2
      tncoff6 = z3*t2*t4/288+t9*t12*t14/48+t18*t11*t20/24+t24*t2*t14*t19
     #/48+t18*t11*t30*t19/12+t18*t11*t35/24-t40*t12/576-t40*t4*t19/576+t
     #24*t2*t46/192+t24*t2*t20/48+t18*t20*t19/72+t9*t12*t30/24+t9*t12*t1
     #9/24-t39*t4*t9/1152
      t1 = -1+l1
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t5 = z-1
      t6 = 1/t5
      t8 = dlog(-t1)
      t13 = dlog(-t5)
      t14 = t13*t3
      t15 = t2*t6
      t19 = dlog(l1)
      t21 = t19**2
      t28 = 0.3141592653589793D1**2
      t29 = t28*t3
      t33 = 1/l2
      t34 = t33*t3
      t38 = t33*t6
      t42 = t8**2
      t49 = dlog(l2)
      t53 = -l2+1
      t55 = t53*t33*t3
      t57 = -1/t53
      t58 = t2*t57
      t67 = t13**2
      t71 = (t4*t6*t8/96+t14*t15/48)*t19+t4*t6*t21/192+z3*t3*t2/576-t29*
     #t2*t8/2304-t34*t2*t21/192-t4*t38*t8/48+t4*t6*t42/192+t14*t15*t8/48
     #+(-t34*t2*t49/48+t55*t58*t8/96)*t19+t4*t21*t19/1152+t4*t6*t67/48
      t93 = t49**2
      t118 = t55*t58*t8*t49/48-t28*t2*t3*t6/1152+t53*t2*t3*t57*t33*t6*t1
     #9/48+t4*t8*t21/384+t55*t58*t42/192+t55*t58*t93/48-t4*t38*t13/24+(t
     #4*t42/384-t29*t2/2304)*t19-t4*t38*t49/24+t4*t42*t8/1152+(-t28+t28*
     #l2)*t33*t4*t57/1152
      tncoff7 = t71+t118
      t1 = 1/l1
      t2 = -1+l1
      t3 = 1/t2
      t4 = t1*t3
      t6 = 1/(-1+l2)
      t7 = dlog(l1)
      t8 = t7**2
      t12 = z-1
      t14 = t6/t12
      t16 = dlog(-t2)
      t27 = t16**2
      t32 = dlog(-t12)
      t36 = 0.3141592653589793D1**2
      tncoff8 = t4*t6*t8/192+t4*t14*t16/48+t4*t6*t16*t7/96+t4*t14*t7/48+
     #t4*t6*t27/192+t4*t14*t32/24-t36*t1*t3*t6/1152
      t1 = -1+l1
      t2 = 1/t1
      t3 = 1/l1
      t6 = 1/(-1+l2)
      t7 = t2*t3*t6
      t8 = 1/l2
      t9 = z-1
      t11 = t8/t9
      t13 = dlog(-t9)
      t17 = dlog(l2)
      t21 = 0.3141592653589793D1**2
      t25 = t8*t3
      t26 = t2*t6
      t27 = dlog(l1)
      t28 = t27**2
      t33 = dlog(-t1)
      t50 = t33**2
      t57 = t17**2
      tncoff9 = -t7*t11*t13/24-t7*t11*t17/24+t21*t8*t7/1152-t25*t26*t28/
     #192-t7*t11*t33/48-t25*t2*t6*t33*t17/48+(-t25*t26*t17/48-t25*t26*t3
     #3/96)*t27-t25*t26*t50/192-t7*t11*t27/48-t25*t26*t57/48
      t1 = l2*z
      t2 = z*l1
      t4 = -z+t1+t2-t2*l2
      t6 = -1+l1
      t7 = 1/t6
      t9 = 1/l1
      t12 = 1/(-1+l2)
      t13 = 1/l2
      t14 = t12*t13
      t15 = z-1
      t16 = 1/t15
      t17 = dlog(z)
      t18 = t16*t17
      t22 = -t2+z
      t24 = t9*t16
      t26 = dlog(-t15)
      t27 = t26**2
      t31 = -z+t1
      t34 = t13*t16
      t35 = t17**2
      t39 = 0.3141592653589793D1**2
      t40 = z*t39
      t53 = t4*t7*t9
      t58 = t31*t12
      t59 = dlog(l2)
      t60 = t59**2
      t75 = dlog(l1)
      t80 = z*t16
      t84 = t4*t7*t9*t14*t18/24+t22*t7*t24*t27/48-t31*t12*t34*t35/24-(-t
     #40*l1+t40)*t7*t24/1152+t26*t31*t12*t34*t17/12-t53*t14*t16*t26/24-t
     #58*t34*t60/24-t53*t14*t16*t59/24+t27*t26*z*t16/72-t4*t7*t9*t14*t16
     #*t75/48-t80*t35*t17/72
      t92 = t26*z
      t116 = -t22
      t125 = t75**2
      t131 = -z*t6*t7
      t144 = (t80*t13*t17/12-t58*t34*t26/12)*t59+t92*t16*t35/24+(-t27*z*
     #t16/24+t40*t16/576)*t17-t92*t39*t16/576-t31*t13*t12*t16*t27/24+z3*
     #z*t16/288+(z*t9*t18/48-t26*t116*t9*t7*t16/48)*t75-t116*t7*t24*t125
     #/192+t131*t24*t35/48-t131*t24*t26*t17/24+(t40*l2-t40)*t12*t34/576
      tncoff10 = t84+t144
      t2 = l1*l2-l2
      t4 = 1/(-1+l1)
      t6 = 1/l1
      t7 = t2*t4*t6
      t9 = 1/(-1+l2)
      t11 = t9/l2
      t12 = z-1
      t13 = 1/t12
      t15 = dlog(-t12)
      t23 = dlog(l1)
      t28 = t9*t13
      t29 = t15**2
      t32 = 0.3141592653589793D1**2
      t36 = dlog(z)
      t37 = t13*t36
      t41 = t36**2
      tncoff11 = -t7*t11*t13*t15/24-t2*t4*t6*t11*t13*t23/48+t28*t29/24-t
     #32*t9*t13/576+t7*t11*t37/48+t28*t41/96-t15*t9*t37/24
      t1 = l2*z
      t2 = l2-t1+z
      t3 = 1/l2
      t6 = 1/(-1+l2)
      t7 = z-1
      t8 = 1/t7
      t11 = dlog(-t7)
      t12 = t11**2
      t16 = -t2
      t17 = 1-l1
      t20 = -1/t17
      t21 = 1/l1
      t23 = t16*t17*t20*t21
      t24 = t6*t3
      t27 = dlog(-z/t16)
      t32 = t16*t6
      t34 = dlog(t2)
      t35 = t8*t34
      t39 = t3*t8
      t40 = t34**2
      t44 = 0.3141592653589793D1**2
      t45 = z*t44
      t52 = z*l1
      t55 = t52*l2-t52-l1*l2-t1+z+l2
      t57 = t55*t20*t21
      t62 = t39*t11
      t65 = t32*t39*t34
      t72 = t2*t6
      t77 = dlog(l2)
      t79 = t27**2
      t87 = t77**2
      t94 = dlog(l1)
      tncoff12 = -t2*t3*t6*t8*t12/24-t23*t24*t8*t27/48-t32*t3*t35*t11/12
     #+t32*t39*t40/24+(t45-t45*l2+t44*l2)*t6*t39/576-t57*t24*t8*t11/24+(
     #-t32*t62/24+t65/24)*t27+(-t32*t39*t27/24-t72*t62/12-t65/12)*t77+t3
     #2*t39*t79/96-t57*t24*t8*t77/24-t72*t39*t87/24-t55*t20*t21*t24*t8*t
     #94/48-t23*t24*t35/24
      t1 = z*l1
      t3 = t1*l2-t1
      t4 = -1+l1
      t5 = 1/t4
      t7 = 1/l1
      t8 = t3*t5*t7
      t10 = 1/(-1+l2)
      t11 = 1/l2
      t12 = t10*t11
      t13 = z-1
      t14 = 1/t13
      t16 = dlog(-t13)
      t17 = t14*t16
      t21 = dlog(l2)
      t26 = z*t5
      t27 = dlog(z)
      t28 = t27**2
      t29 = t14*t28
      t33 = dlog(-t4)
      t34 = t14*t33
      t44 = t16*z
      t47 = z*t14
      t59 = -z+l2*z
      t62 = t11*t14
      t71 = t59*t10
      t77 = -t8*t12*t17/24-t8*t12*t14*t21/24+t26*t29/48-t8*t12*t34/48+(-
     #t26*t34/48-t26*t17/24)*t27+t44*t29/24-t47*t28*t27/72+t3*t5*t7*t12*
     #t14*t27/24-t59*t10*t62*t28/24+t26*t17*t33/48+(t47*t11*t27/12-t71*t
     #62*t16/12)*t21
      t84 = t16**2
      t88 = 0.3141592653589793D1**2
      t89 = z*t88
      t98 = t21**2
      t102 = t33**2
      t128 = t16*t59*t10*t62*t27/12+(-t84*z*t14/24+t89*t14/576)*t27+t84*
     #t16*z*t14/72-t71*t62*t98/24+t26*t14*t102/192+t26*t14*t84/48-t59*t1
     #1*t10*t14*t84/24+z3*z*t14/288-t44*t88*t14/576+(t89*l2-t89)*t10*t62
     #/576-t89*t5*t14/1152
      tncoff13 = t77+t128
      t1 = -1+l1
      t4 = 1/(-1+l2)
      t5 = 1/t1*t4
      t6 = z-1
      t7 = 1/t6
      t9 = dlog(-t6)
      t14 = dlog(-t1)
      t18 = t4*t7
      t19 = dlog(z)
      t20 = t19**2
      t23 = t7*t19
      t29 = 0.3141592653589793D1**2
      t33 = t9**2
      tncoff14 = t5*t7*t9/24+t5*t7*t14/48+t18*t20/96-t5*t23/48-t9*t4*t23
     #/24-t29*t4*t7/576+t18*t33/24
      t2 = l2-l2*z+z
      t4 = 1/(-1+l2)
      t5 = t2*t4
      t6 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t9 = t6*t8
      t10 = dlog(l2)
      t11 = t10**2
      t15 = -t2
      t16 = t15*t4
      t19 = dlog(-z/t15)
      t20 = t9*t19
      t24 = dlog(-t7)
      t25 = t9*t24
      t28 = dlog(t2)
      t29 = t9*t28
      t30 = t16*t29
      t34 = t19**2
      t47 = 0.3141592653589793D1**2
      t48 = z*t47
      t57 = t24**2
      t61 = t28**2
      t65 = z*l1
      t69 = -1+l1
      t70 = 1/t69
      t73 = (-t65*l2+l1*l2+t65)*t70/l1
      t74 = t4*t6
      t84 = t15*t70*t4
      t88 = dlog(-t69)
      tncoff15 = -t5*t9*t11/24+(-t16*t20/24-t5*t25/12-t30/12)*t10+t16*t9
     #*t34/96+(-t16*t25/24+t30/24)*t19-t16*t6*t8*t28*t24/12+(t48-t48*l2+
     #t47*l2)*t4*t9/576-t2*t6*t4*t8*t57/24+t16*t9*t61/24-t73*t74*t8*t24/
     #24-t73*t74*t8*t10/24-t84*t29/24-t73*t74*t8*t88/48-t84*t20/48
      t1 = -1+l1
      t2 = 1/t1
      t3 = z*t2
      t4 = 1/l1
      t5 = z-1
      t6 = 1/t5
      t7 = t4*t6
      t8 = dlog(z)
      t9 = t8**2
      t14 = -l2*z+z
      t16 = t14*t2*t4
      t20 = 1/(-1+l2)/l2
      t22 = dlog(-t1)
      t27 = dlog(l1)
      t28 = t27**2
      t33 = t3*t7*t22
      t36 = dlog(-t5)
      t38 = t3*t7*t36
      t42 = t22**2
      t47 = t6*t36
      t58 = 0.3141592653589793D1**2
      t64 = t36**2
      t78 = dlog(l2)
      tncoff16 = -t3*t7*t9/48-t16*t20*t6*t22/48-t3*t7*t28/192+(t33/48+t3
     #8/24)*t8-t3*t7*t42/192-t3*t4*t47*t22/48+t14*t2*t4*t20*t6*t8/24+z*t
     #58*t2*t4*t6/1152-t3*t7*t64/48-t14*t2*t4*t20*t6*t27/48-t16*t20*t47/
     #24-t16*t20*t6*t78/24+(t3*t7*t8/48-t33/96-t38/48)*t27
      t1 = -1+l1
      t4 = 1/t1/l1
      t7 = z-1
      t9 = 1/(-1+l2)/t7
      t11 = dlog(-t1)
      t15 = dlog(z)
      t19 = dlog(l1)
      t24 = dlog(-t7)
      tncoff17 = -t4*t9*t11/48+t4*t9*t15/48-t4*t9*t19/48-t4*t9*t24/24
      t2 = l2-l2*z+z
      t3 = -1+l1
      t4 = 1/t3
      t6 = 1/l1
      t11 = 1/(-1+l2)/l2
      t12 = z-1
      t13 = 1/t12
      t14 = dlog(l1)
      t19 = -t2
      t21 = t19*t4*t6
      t23 = dlog(-t12)
      t28 = dlog(l2)
      t33 = dlog(t2)
      t39 = dlog(-t3)
      t46 = dlog(-z/t19)
      tncoff18 = t2*t4*t6*t11*t13*t14/48-t21*t11*t13*t23/24-t21*t11*t13*
     #t28/24+t21*t11*t13*t33/24-t21*t11*t13*t39/48+t21*t11*t13*t46/48
      t2 = -l1-l2+l1*l2+1
      t3 = 1/l2
      t5 = 1/l1
      t6 = t2*t3*t5
      t7 = -1+l1
      t8 = 1/t7
      t9 = -1+l2
      t10 = 1/t9
      t11 = t8*t10
      t12 = dlog(l1)
      t16 = 1-z
      t17 = dlog(t16)
      t18 = -t7
      t22 = -1/t16
      t27 = -t9*t10
      t28 = t3*t22
      t29 = dlog(l2)
      t36 = t12**2
      t39 = t17**2
      t49 = t29**2
      t52 = 0.3141592653589793D1**2
      tncoff19 = t6*t11*t12/96+t17*t18*t5*t8*t22/48+t27*t28*t29/24+t6*t1
     #1*t29/48-t5*t36/384+t39*t22/48+z3/576+t2*t8*t5*t10*t3*t22/48-t3*t4
     #9/48-(t52-t52*l1)*t5*t8/2304-t52*t22/1152+(-t52+t52*l2)*t10*t3/115
     #2+t18*t8*t5*t22*t12/96+t27*t28*t17/24
      t2 = l2-l1*l2
      t4 = 1/(-1+l1)
      t6 = 1/l1
      t9 = 1/(-1+l2)
      t10 = 1/l2
      t12 = z-1
      t13 = 1/t12
      t17 = 0.3141592653589793D1**2
      t21 = dlog(-t12)
      t28 = dlog(l1)
      tncoff20 = t2*t4*t6*t9*t10*t13/48-t17*t9/1152+t21*t9*t13/24+t2*t10
     #*t6*t4*t9*t28/96
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t5 = z-1
      t6 = 1/t5
      t8 = dlog(-t5)
      t12 = 0.3141592653589793D1**2
      t16 = 1/l1
      t17 = t3*t16
      t18 = dlog(l2)
      t22 = dlog(l1)
      t29 = t18**2
      tncoff21 = -t4*t6*t8/24+t12*t2*t3/1152+t17*t2*t18/48+t17*t2*t22/96
     #-t4*t6*t18/24-t4*t29/48+t16*t2*t3*t6/48
      t1 = -1+l1
      t2 = 1/t1
      t3 = z-1
      t4 = 1/t3
      t5 = t2*t4
      t7 = dlog(-t3)
      t10 = 1/l2
      t11 = t10*t4
      t14 = 0.3141592653589793D1**2
      t18 = -l1*l2+l1
      t20 = 1/l1
      t23 = 1/(-1+l2)
      t28 = dlog(l2)
      t29 = t28**2
      t33 = dlog(-t1)
      t34 = t33**2
      t39 = t7**2
      t49 = t18*t10*t20
      t50 = t2*t23
      tncoff22 = t5*t7/48-t11*t7/24-t14*t4/1152+t18*t2*t20*t23*t10*t4/48
     #-t10*t29/48+t2*t34/384-t14*t2/2304+t39*t4/48+(-t14+t14*l2)*t23*t10
     #/1152+z3/576+t49*t50*t28/48+t5*t33/96-t11*t28/24+t49*t50*t33/96
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 1/(-1+l2)
      t7 = -1/t1
      t10 = -1+l1
      t12 = 1/t10*t4
      t14 = dlog(-t10)
      t19 = 0.3141592653589793D1**2
      tncoff23 = t2*t4*t7/24+t12*t14/96+t12*t7/48-t19*t4/1152
      t1 = 0.3141592653589793D1**2
      t3 = 1/(-1+l2)
      t5 = 1/l2
      t8 = t3*t5
      t9 = z-1
      t10 = 1/t9
      t11 = dlog(l2)
      t16 = dlog(-t9)
      t20 = -1+l1
      t21 = 1/t20
      t22 = t5*t21
      t26 = t11**2
      t30 = dlog(-t20)
      tncoff24 = t1*t3*t5/1152-t8*t10*t11/24-t8*t10*t16/24-t22*t3*t11/48
     #-t8*t26/48-t22*t3*t30/96-t21*t3*t5*t10/48
      t1 = -1+l1
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t5 = 1/l2
      t6 = z-1
      t7 = 1/t6
      t11 = 0.3141592653589793D1**2
      t15 = dlog(l1)
      t20 = dlog(-t1)
      t21 = t20**2
      t24 = t3*t5
      t25 = dlog(l2)
      t36 = dlog(-t6)
      t41 = t15**2
      tncoff25 = t4*t5*t7/48+t11*t3*t2/2304-t4*t7*t15/96-t4*t21/384+t24*
     #t2*t25/48-t4*t7*t20/96+t24*t2*t15/96-t36*t3*t2*t7/48-t4*t41/384+t2
     #4*t2*t20/96-t4*t20*t15/192
      t2 = -1+l1
      t4 = 1/l1/t2
      t6 = 1/(-1+l2)
      t8 = dlog(-t2)
      t17 = dlog(l1)
      tncoff26 = -t4*t6*t8/96-t4*t6/(z-1)/48-t4*t6*t17/96
      t1 = 1/l2
      t2 = 1/l1
      t3 = t1*t2
      t4 = -1+l1
      t5 = 1/t4
      t7 = 1/(-1+l2)
      t8 = t5*t7
      t9 = dlog(l2)
      t13 = dlog(l1)
      t18 = dlog(-t4)
      tncoff27 = t3*t8*t9/48+t3*t8*t13/96+t3*t8*t18/96+t5*t2*t7*t1/(z-1)
     #/48
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t6 = -1/t1
      t9 = z*l1
      t10 = t9-z
      t12 = 1/l1
      t14 = 1/(-1+l1)
      t19 = l2*z
      t20 = -z+t19
      t22 = 1/(-1+l2)
      t23 = t20*t22
      t24 = 1/l2
      t25 = t24*t6
      t26 = dlog(l2)
      t31 = dlog(z)
      t32 = t31**2
      t36 = t6*t31
      t47 = 0.3141592653589793D1**2
      t53 = dlog(l1)
      tncoff28 = -t3*z*t6/48+t2*t10*t12*t14*t6/48+t23*t25*t26/24-z*t6*t3
     #2/48+t2*z*t36/24+(-z+t19+t9-t9*l2)*t14*t12*t22*t24*t6/48+z*t47*t6/
     #1152+t10*t14*t12*t6*t53/96+t23*t25*t2/24-z*t12*t36/48-t20*t22*t25*
     #t31/24
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 1/(-1+l2)
      t7 = -1/t1
      t11 = dlog(z)
      tncoff29 = -t2*t4*t7/24+t4*t7*t11/48+(l1*l2-l2)/(-1+l1)/l1*t4/l2*t
     #7/48
      t1 = l2*z
      t2 = l2-t1+z
      t4 = 1/(-1+l2)
      t5 = t2*t4
      t6 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t9 = t6*t8
      t11 = dlog(-t7)
      t15 = -t2
      t16 = t15*t4
      t17 = dlog(t2)
      t23 = dlog(-z/t15)
      t27 = dlog(l2)
      t31 = z*l1
      tncoff30 = t5*t9*t11/24+t16*t9*t17/24+t16*t9*t23/48+t5*t9*t27/24+(
     #t31*l2-t31-l1*l2-t1+z+l2)/(-1+l1)/l1*t4*t6*t8/48
      t2 = -l2*z+z
      t4 = 1/(-1+l2)
      t5 = t2*t4
      t6 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t9 = t6*t8
      t11 = dlog(-t7)
      t15 = -1+l1
      t16 = 1/t15
      t17 = z*t16
      t21 = 0.3141592653589793D1**2
      t25 = z*l1
      t35 = t11**2
      t39 = dlog(z)
      t40 = t8*t39
      t44 = t39**2
      t48 = dlog(-t15)
      t52 = dlog(l2)
      tncoff31 = -t5*t9*t11/24-t17*t8*t11/48+z*t21*t8/1152+(t25*l2-t25)*
     #t16/l1*t4*t6*t8/48-t35*z*t8/48+t17*t40/48-z*t8*t44/48-t17*t8*t48/9
     #6-t5*t9*t52/24+t2*t4*t9*t39/24+t11*z*t40/24
      t4 = 1/(-1+l2)
      t6 = z-1
      t7 = 1/t6
      t11 = dlog(z)
      t15 = dlog(-t6)
      tncoff32 = -1/(-1+l1)*t4*t7/48+t4*t7*t11/48-t4*t7*t15/24
      t2 = -l2+l2*z-z
      t4 = 1/(-1+l2)
      t5 = t2*t4
      t6 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t9 = t6*t8
      t10 = dlog(l2)
      t16 = dlog(-z/t2)
      t21 = dlog(-t7)
      t25 = z*l1
      t39 = dlog(-t2)
      tncoff33 = -t5*t9*t10/24+t5*t9*t16/48-t5*t9*t21/24+(-t25*l2+l1*l2+
     #t25)/(-1+l1)/l1*t4*t6*t8/48+t5*t9*t39/24
      t3 = -1+l1
      t4 = 1/t3
      t6 = 1/l1
      t12 = z-1
      t13 = 1/t12
      t17 = z*t4
      t18 = t6*t13
      t19 = dlog(l1)
      t23 = dlog(z)
      t28 = dlog(-t3)
      t33 = dlog(-t12)
      tncoff34 = (-l2*z+z)*t4*t6/(-1+l2)/l2*t13/48+t17*t18*t19/96-t17*t1
     #8*t23/48+t17*t18*t28/96+t17*t18*t33/48
      tncoff35 = 1/(-1+l1)/l1/(-1+l2)/(z-1)/48
      tncoff36 = (-l2+l2*z-z)/(-1+l1)/l1/(-1+l2)/l2/(z-1)/48
      t1 = 1-l1
      t2 = 12*t1
      t4 = -1/t1
      t6 = 1/l1
      t7 = z-1
      t8 = 1/t7
      t12 = 0.3141592653589793D1**2
      t17 = 1/l2
      t20 = -1+l2
      t21 = 1/t20
      t26 = dlog(-t7)
      t30 = -12*t20
      t32 = dlog(l2)
      t42 = dlog(l1)
      tncoff37 = -t2*t4*t6*t8/1152+t12/2304+(12*l2-12*l1*l2-12+12*l1)*t1
     #7*t6*t4*t21/1152-t8*t26/48-t30*t21*t17*t32/576-t30*t21*t17*t8/576-
     #t2*t6*t4*t42/2304
      t11 = 1/(-1+l2)
      tncoff38 = (-12*l2+12*l1*l2)/l2/l1/(-1+l1)*t11/1152-t11/(z-1)/48
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t9 = 1-l1
      t19 = dlog(l2)
      tncoff39 = t4/(z-1)/48-t3/l1*t2/96+t4*t19/48
      t1 = 1-z
      t2 = dlog(t1)
      t4 = -1/t1
      t7 = -1+l2
      t9 = 1/t7
      t10 = 12*t7*t9
      t11 = 1/l2
      t15 = -1+l1
      t16 = 1/t15
      t19 = dlog(l2)
      t33 = dlog(-t15)
      t36 = 0.3141592653589793D1**2
      tncoff40 = -t2*t4/48+t10*t11*t4/576-t16*t4/96+t10*t11*t19/576+(12*
     #l1*l2-12*l1)*t11/l1*t16*t9/1152-t16*t33/192+t36/2304
      t2 = 1/(-1+l2)
      tncoff41 = -t2/(z-1)/48-1/(-1+l1)*t2/96
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t9 = dlog(l2)
      tncoff42 = t4/(z-1)/48+t4*t9/48+t3/(-1+l1)*t2/96
      t2 = -1+l1
      t4 = 1/l1/t2
      t6 = dlog(-t2)
      t9 = dlog(l1)
      t16 = -l2+1
      tncoff43 = t4*t6/192+t4*t9/192+t4/(z-1)/96-1/l2*t4/96
      tncoff44 = 1/l1/(-1+l1)/(-1+l2)/96
      tncoff45 = -1/l2/l1/(-1+l1)/(-1+l2)/96
      t8 = z-1
      t9 = 1/t8
      t14 = dlog(z)
      t18 = dlog(-t8)
      tncoff46 = -(12*z*l1-12*z)/(-1+l1)/l1*t9/1152-z*t9*t14/48+t18*z*t9
     #/48+(12*z-12*l2*z)/(-1+l2)/l2*t9/576
      tncoff47 = 1/(-1+l2)/(z-1)/48
      tncoff48 = (12*l2*z-12*l2-12*z)/(-1+l2)/l2/(z-1)/576
      t8 = z-1
      t9 = 1/t8
      t19 = dlog(z)
      t23 = dlog(-t8)
      tncoff49 = (12*z-12*l2*z)/(-1+l2)/l2*t9/576+z/(-1+l1)*t9/96-z*t9*t
     #19/48+t23*z*t9/48
      tncoff50 = 1/(-1+l2)/(z-1)/48
      tncoff51 = (12*l2*z-12*l2-12*z)/(-1+l2)/l2/(z-1)/576
      tncoff52 = -z/(-1+l1)/l1/(z-1)/96
      t1 = -1+l1
      t8 = -l2+1
      tncoff53 = -1/l1/192-1/l2/96+1/(z-1)/96
      tncoff54 = 1/(-1+l2)/96
      tncoff55 = -1/(-1+l2)/l2/96
      t4 = -l2+1
      tncoff56 = 1/(z-1)/96-1/l2/96+1/(-1+l1)/192
      tncoff57 = 1/(-1+l2)/96
      tncoff58 = -1/(-1+l2)/l2/96
      tncoff59 = -1/(-1+l1)/l1/192
      tncoff60 = -z/(z-1)/96
      tncoff61 = -z/(z-1)/96
      tncoff62 = -1.D0/192.D0
      tncoff63 = -1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum21e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum21e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum21e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum21e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum21e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum21e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum21e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum21e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum21e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum21e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum21e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = qqbcnum21e1(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = qqbcnum21e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = qqbcnum21e1(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = qqbcnum21e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff28)
      tcnum29 = qqbcnum21e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff29)
      tcnum30 = qqbcnum21e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff30)
      tcnum31 = qqbcnum21e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff31)
      tcnum32 = qqbcnum21e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff32)
      tcnum33 = qqbcnum21e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff33)
      tcnum34 = qqbcnum21e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff34)
      tcnum35 = qqbcnum21e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff35)
      tcnum36 = qqbcnum21e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff36)
      tcnum37 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff37)
      tcnum38 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff38)
      tcnum39 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff39)
      tcnum40 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff40)
      tcnum41 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff41)
      tcnum42 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff42)
      tcnum43 = qqbcnum21e2(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff43)
      tcnum44 = qqbcnum21e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff44)
      tcnum45 = qqbcnum21e2(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff45)
      tcnum46 = qqbcnum21e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff46)
      tcnum47 = qqbcnum21e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff47)
      tcnum48 = qqbcnum21e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff48)
      tcnum49 = qqbcnum21e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff49)
      tcnum50 = qqbcnum21e2(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff50)
      tcnum51 = qqbcnum21e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff51)
      tcnum52 = qqbcnum21e2(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff52)
      tcnum53 = qqbcnum21e3(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff53)
      tcnum54 = qqbcnum21e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff54)
      tcnum55 = qqbcnum21e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff55)
      tcnum56 = qqbcnum21e3(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff56)
      tcnum57 = qqbcnum21e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff57)
      tcnum58 = qqbcnum21e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff58)
      tcnum59 = qqbcnum21e3(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff59)
      tcnum60 = qqbcnum21e3(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff60)
      tcnum61 = qqbcnum21e3(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff61)
      tcnum62 = qqbcnum21e4(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff62)
      tcnum63 = qqbcnum21e4(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff63)

CCC Sum

      t16 = tcnum48*tncoff48+tcnum50*tncoff50+tcnum32*tncoff32+tcnum54*t
     #ncoff54+tcnum24*tncoff24+tcnum1*tncoff1+tcnum63*tncoff63+tcnum60*t
     #ncoff60+tcnum43*tncoff43+tcnum25*tncoff25+tcnum34*tncoff34+tcnum61
     #*tncoff61+tcnum40*tncoff40+tcnum42*tncoff42+tcnum15*tncoff15
      t33 = tcnum44*tncoff44+tcnum37*tncoff37+tcnum49*tncoff49+tcnum52*t
     #ncoff52+tcnum33*tncoff33+tcnum41*tncoff41+tcnum2*tncoff2+tcnum59*t
     #ncoff59+tcnum26*tncoff26+tcnum53*tncoff53+tcnum36*tncoff36+tcnum13
     #*tncoff13+tcnum35*tncoff35+tcnum17*tncoff17+tcnum38*tncoff38+tcnum
     #51*tncoff51
      t51 = tcnum12*tncoff12+tcnum30*tncoff30+tcnum47*tncoff47+tcnum21*t
     #ncoff21+tcnum7*tncoff7+tcnum20*tncoff20+tcnum22*tncoff22+tcnum8*tn
     #coff8+tcnum5*tncoff5+tcnum29*tncoff29+tcnum58*tncoff58+tcnum55*tnc
     #off55+tcnum4*tncoff4+tcnum39*tncoff39+tcnum3*tncoff3+tcnum28*tncof
     #f28
      t68 = tcnum6*tncoff6+tcnum27*tncoff27+tcnum19*tncoff19+tcnum14*tnc
     #off14+tcnum9*tncoff9+tcnum18*tncoff18+tcnum11*tncoff11+tcnum56*tnc
     #off56+tcnum31*tncoff31+tcnum46*tncoff46+tcnum62*tncoff62+tcnum10*t
     #ncoff10+tcnum45*tncoff45+tcnum23*tncoff23+tcnum16*tncoff16+tcnum57
     #*tncoff57
      qqbc21e0int = t16+t33+t51+t68

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc21e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum21e0
      external qqbcnum21e0
      double precision qqbcnum21e1
      external qqbcnum21e1
      double precision qqbcnum21e2
      external qqbcnum21e2
      double precision qqbcnum21e3
      external qqbcnum21e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1/l2
      t2 = dlog(l2)
      t3 = t2**2
      t6 = -l2+1
      t8 = -1/t6
      t9 = t6*t8
      t10 = z-1
      t11 = 1/t10
      t12 = t1*t11
      t17 = dlog(-t10)
      t21 = t17**2
      t24 = 0.3141592653589793D1**2
      t27 = 1/l1
      t28 = t27*t11
      t32 = -l1-l2+l1*l2+1
      t35 = -1+l1
      t36 = 1/t35
      t38 = dlog(l1)
      t52 = t38**2
      t58 = t8*t1
      tncoff1 = -t1*t3/48+t9*t12*t2/24+t9*t12*t17/24+t21*t11/48-t24*t11/
     #1152-t28*t17/48+t32*t1*t27*t36*t8*t38/96-(t24-t24*l1)*t27*t36/2304
     #-t35*t36*t28*t38/96-t27*t52/384+t32*t36*t27*t58*t2/48+t32*t36*t27*
     #t58*t11/48+(-t24+t24*l2)*t8*t1/1152+z3/576
      t2 = l2-l1*l2
      t4 = 1/(-1+l1)
      t6 = 1/l1
      t9 = 1/(-1+l2)
      t10 = 1/l2
      t12 = z-1
      t13 = 1/t12
      t17 = 0.3141592653589793D1**2
      t21 = dlog(-t12)
      t28 = dlog(l1)
      tncoff2 = t2*t4*t6*t9*t10*t13/48-t17*t9/1152+t21*t9*t13/24+t2*t10*
     #t6*t4*t9*t28/96
      t1 = 1/l2
      t2 = 1/l1
      t5 = 1/(-1+l2)
      t6 = dlog(l1)
      t10 = t5*t1
      t11 = dlog(l2)
      t12 = t11**2
      t15 = 1-l1
      t24 = z-1
      t25 = 1/t24
      t29 = 0.3141592653589793D1**2
      t37 = dlog(-t24)
      tncoff3 = t1*t2*t5*t6/96-t10*t12/48+t2*t10*t11/48+t2*t5*t1*t25/48+
     #t29*t5*t1/1152-t10*t25*t11/24-t10*t25*t37/24
      t2 = -l1+l1*l2
      t3 = -1+l1
      t4 = 1/t3
      t6 = 1/l1
      t7 = t2*t4*t6
      t8 = -1+l2
      t9 = 1/t8
      t10 = 1/l2
      t11 = t9*t10
      t12 = dlog(l2)
      t17 = dlog(-t3)
      t24 = z-1
      t25 = 1/t24
      t29 = 0.3141592653589793D1**2
      t38 = t4*t25
      t40 = dlog(-t24)
      t44 = -t8*t9
      t45 = t10*t25
      t49 = t40**2
      t59 = t12**2
      t62 = t17**2
      tncoff4 = -t7*t11*t12/48-t7*t11*t17/96-t2*t4*t6*t11*t25/48-(-t29*l
     #2+t29)*t9*t10/1152+z3/576-t29*t25/1152+t38*t40/48+t44*t45*t12/24+t
     #49*t25/48+t44*t45*t40/24-t29*t4/2304+t38*t17/96-t10*t59/48+t4*t62/
     #384
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 1/(-1+l2)
      t7 = -1/t1
      t10 = -1+l1
      t12 = 1/t10*t4
      t14 = dlog(-t10)
      t19 = 0.3141592653589793D1**2
      tncoff5 = t2*t4*t7/24+t12*t14/96+t12*t7/48-t19*t4/1152
      t1 = 0.3141592653589793D1**2
      t3 = 1/(-1+l2)
      t5 = 1/l2
      t8 = t3*t5
      t9 = z-1
      t10 = 1/t9
      t11 = dlog(l2)
      t16 = dlog(-t9)
      t20 = -1+l1
      t21 = 1/t20
      t22 = t5*t21
      t26 = t11**2
      t30 = dlog(-t20)
      tncoff6 = t1*t3*t5/1152-t8*t10*t11/24-t8*t10*t16/24-t22*t3*t11/48-
     #t8*t26/48-t22*t3*t30/96-t21*t3*t5*t10/48
      t1 = -l2+1
      t2 = -1+l1
      t3 = 1/t2
      t5 = 1/l1
      t8 = -1/t1
      t9 = 1/l2
      t11 = dlog(l2)
      t19 = dlog(-t2)
      t23 = t3*t5
      t24 = z-1
      t25 = 1/t24
      t29 = t19**2
      t32 = dlog(l1)
      t33 = t32**2
      t37 = dlog(-t24)
      t46 = 0.3141592653589793D1**2
      tncoff7 = -t1*t3*t5*t8*t9*t11/48-t1*t9*t5*t3*t8*t19/96+t23*t9*t25/
     #48-t23*t29/384-t23*t33/384-t37*t5*t3*t25/48+t9*t5*t3*t32/96+t46*t5
     #*t3/2304-t23*t19*t32/192-t23*t25*t32/96-t23*t25*t19/96
      t2 = -1+l1
      t4 = 1/l1/t2
      t6 = 1/(-1+l2)
      t8 = dlog(-t2)
      t17 = dlog(l1)
      tncoff8 = -t4*t6*t8/96-t4*t6/(z-1)/48-t4*t6*t17/96
      t1 = 1/l2
      t2 = 1/l1
      t3 = t1*t2
      t4 = -1+l1
      t5 = 1/t4
      t7 = 1/(-1+l2)
      t8 = t5*t7
      t9 = dlog(l2)
      t13 = dlog(l1)
      t18 = dlog(-t4)
      tncoff9 = t3*t8*t9/48+t3*t8*t13/96+t3*t8*t18/96+t5*t2*t7*t1/(z-1)/
     #48
      t1 = 0.3141592653589793D1**2
      t3 = z-1
      t4 = 1/t3
      t8 = dlog(z)
      t9 = t8**2
      t12 = z*l1
      t13 = t12-z
      t14 = -1+l1
      t15 = 1/t14
      t17 = 1/l1
      t18 = t17*t4
      t19 = dlog(l1)
      t30 = dlog(-t3)
      t35 = l2*z
      t41 = 1/(-1+l2)
      t42 = 1/l2
      t47 = -z+t35
      t48 = t47*t41
      t49 = t42*t4
      t50 = dlog(l2)
      t59 = t30**2
      tncoff10 = z*t1*t4/1152-z*t4*t9/48+t13*t15*t18*t19/96-z*t14*t15*t1
     #8*t8/48+t30*z*t4*t8/24+(-z+t35+t12-t12*l2)*t15*t17*t41*t42*t4/48+t
     #48*t49*t50/24-t47*t41*t49*t8/24-t59*z*t4/48+t48*t49*t30/24+t13*t15
     #*t18*t30/48
      t1 = 1-z
      t2 = dlog(t1)
      t4 = 1/(-1+l2)
      t7 = -1/t1
      t11 = dlog(z)
      tncoff11 = -t2*t4*t7/24+t4*t7*t11/48+(l1*l2-l2)/(-1+l1)/l1*t4/l2*t
     #7/48
      t1 = l2*z
      t2 = l2-t1+z
      t4 = 1/(-1+l2)
      t5 = t2*t4
      t6 = 1/l2
      t7 = z-1
      t8 = 1/t7
      t9 = t6*t8
      t11 = dlog(-t7)
      t15 = -t2
      t16 = t15*t4
      t17 = dlog(t2)
      t23 = dlog(-z/t15)
      t27 = dlog(l2)
      t31 = z*l1
      tncoff12 = t5*t9*t11/24+t16*t9*t17/24+t16*t9*t23/48+t5*t9*t27/24+(
     #t31*l2-t31-l1*l2-t1+z+l2)/(-1+l1)/l1*t4*t6*t8/48
      t1 = z*l1
      t4 = -1+l1
      t5 = 1/t4
      t10 = 1/(-1+l2)
      t11 = 1/l2
      t13 = z-1
      t14 = 1/t13
      t19 = -z+l2*z
      t20 = t19*t10
      t21 = t11*t14
      t22 = dlog(l2)
      t27 = dlog(-t13)
      t28 = t27**2
      t32 = z*t5
      t38 = dlog(z)
      t42 = 0.3141592653589793D1**2
      t47 = t14*t38
      t51 = t38**2
      t60 = dlog(-t4)
      tncoff13 = (t1*l2-t1)*t5/l1*t10*t11*t14/48+t20*t21*t22/24-t28*z*t1
     #4/48-t32*t14*t27/48-t19*t10*t21*t38/24+z*t42*t14/1152+t27*z*t47/24
     #-z*t14*t51/48+t20*t21*t27/24+t32*t47/48-t32*t14*t60/96
      t4 = 1/(-1+l2)
      t6 = z-1
      t7 = 1/t6
      t11 = dlog(z)
      t15 = dlog(-t6)
      tncoff14 = -1/(-1+l1)*t4*t7/48+t4*t7*t11/48-t15*t4*t7/24
      t1 = z*l1
      t11 = 1/(-1+l2)
      t12 = 1/l2
      t14 = z-1
      t15 = 1/t14
      t20 = l2-l2*z+z
      t21 = t20*t11
      t22 = t12*t15
      t23 = dlog(l2)
      t27 = -t20
      t28 = t27*t11
      t31 = dlog(-z/t27)
      t36 = dlog(-t14)
      t40 = dlog(t20)
      tncoff15 = (-t1*l2+l1*l2+t1)/(-1+l1)/l1*t11*t12*t15/48+t21*t22*t23
     #/24+t28*t22*t31/48+t21*t22*t36/24+t28*t22*t40/24
      t3 = -1+l1
      t4 = 1/t3
      t6 = 1/l1
      t12 = z-1
      t13 = 1/t12
      t17 = z*t4
      t18 = t6*t13
      t19 = dlog(l1)
      t23 = dlog(z)
      t28 = dlog(-t3)
      t33 = dlog(-t12)
      tncoff16 = (-l2*z+z)*t4*t6/(-1+l2)/l2*t13/48+t17*t18*t19/96-t17*t1
     #8*t23/48+t17*t18*t28/96+t17*t18*t33/48
      tncoff17 = 1/(-1+l1)/l1/(-1+l2)/(z-1)/48
      tncoff18 = (-l2+l2*z-z)/(-1+l1)/l1/(-1+l2)/l2/(z-1)/48
      t3 = -1+l1
      t4 = 1/t3
      t6 = 1/l1
      t7 = -1+l2
      t8 = 1/t7
      t10 = 1/l2
      t14 = dlog(l2)
      t17 = 0.3141592653589793D1**2
      t21 = z-1
      t22 = 1/t21
      t27 = dlog(-t21)
      t35 = dlog(l1)
      tncoff19 = -(-l1-l2+l1*l2+1)*t4*t6*t8*t10/96+t10*t14/48+t17/2304+t
     #3*t4*t6*t22/96-t27*t22/48+t7*t8*t10*t22/48+t6*t35/192
      t8 = 1/(-1+l2)
      tncoff20 = -(l2-l1*l2)/(-1+l1)/l1*t8/l2/96-t8/(z-1)/48
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t13 = dlog(l2)
      tncoff21 = t4/(z-1)/48-1/l1*t2*t3/96+t4*t13/48
      t3 = -1+l1
      t4 = 1/t3
      t10 = 1/l2
      t14 = dlog(l2)
      t17 = z-1
      t18 = 1/t17
      t24 = dlog(-t3)
      t28 = dlog(-t17)
      t31 = 0.3141592653589793D1**2
      tncoff22 = -(-l1*l2+l1)*t4/l1/(-1+l2)*t10/96+t10*t14/48+t10*t18/48
     #-t4*t18/96-t4*t24/192-t28*t18/48+t31/2304
      t2 = 1/(-1+l2)
      tncoff23 = -t2/(z-1)/48-1/(-1+l1)*t2/96
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t9 = dlog(l2)
      tncoff24 = t4/(z-1)/48+t4*t9/48+t3/(-1+l1)*t2/96
      t1 = -1+l1
      t4 = 1/t1/l1
      t5 = dlog(l1)
      t16 = dlog(-t1)
      tncoff25 = t4*t5/192-t4/l2/96+t4/(z-1)/96+t4*t16/192
      tncoff26 = 1/l1/(-1+l1)/(-1+l2)/96
      tncoff27 = -1/l2/l1/(-1+l1)/(-1+l2)/96
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(z)
      t17 = dlog(-t1)
      tncoff28 = -z*t2*t4/48-(z*l1-z)/(-1+l1)/l1*t2/96+t17*z*t2/48-(-z+l
     #2*z)/(-1+l2)/l2*t2/48
      tncoff29 = 1/(-1+l2)/(z-1)/48
      tncoff30 = -(l2-l2*z+z)/(-1+l2)/l2/(z-1)/48
      t7 = z-1
      t8 = 1/t7
      t18 = dlog(z)
      t22 = dlog(-t7)
      tncoff31 = (-l2*z+z)/(-1+l2)/l2*t8/48+z/(-1+l1)*t8/96-z*t8*t18/48+
     #t22*z*t8/48
      tncoff32 = 1/(-1+l2)/(z-1)/48
      tncoff33 = (-l2+l2*z-z)/(-1+l2)/l2/(z-1)/48
      tncoff34 = -z/(-1+l1)/l1/(z-1)/96
      t1 = -1+l1
      t8 = -l2+1
      tncoff35 = -1/l1/192-1/l2/96+1/(z-1)/96
      tncoff36 = 1/(-1+l2)/96
      tncoff37 = -1/(-1+l2)/l2/96
      t1 = -1+l2
      tncoff38 = -1/l2/96+1/(z-1)/96+1/(-1+l1)/192
      tncoff39 = 1/(-1+l2)/96
      tncoff40 = -1/(-1+l2)/l2/96
      tncoff41 = -1/(-1+l1)/l1/192
      tncoff42 = -z/(z-1)/96
      tncoff43 = -z/(z-1)/96
      tncoff44 = -1.D0/192.D0
      tncoff45 = -1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum21e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum21e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum21e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum21e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum21e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum21e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum21e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum21e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum21e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum21e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum21e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = qqbcnum21e1(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = qqbcnum21e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = qqbcnum21e1(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = qqbcnum21e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff28)
      tcnum29 = qqbcnum21e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff29)
      tcnum30 = qqbcnum21e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff30)
      tcnum31 = qqbcnum21e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff31)
      tcnum32 = qqbcnum21e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff32)
      tcnum33 = qqbcnum21e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff33)
      tcnum34 = qqbcnum21e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff34)
      tcnum35 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff35)
      tcnum36 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff36)
      tcnum37 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff37)
      tcnum38 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff38)
      tcnum39 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff39)
      tcnum40 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff40)
      tcnum41 = qqbcnum21e2(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff41)
      tcnum42 = qqbcnum21e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff42)
      tcnum43 = qqbcnum21e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff43)
      tcnum44 = qqbcnum21e3(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff44)
      tcnum45 = qqbcnum21e3(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff45)

CCC Sum

      t12 = tcnum36*tncoff36+tcnum45*tncoff45+tcnum30*tncoff30+tcnum35*t
     #ncoff35+tcnum31*tncoff31+tcnum7*tncoff7+tcnum32*tncoff32+tcnum4*tn
     #coff4+tcnum33*tncoff33+tcnum3*tncoff3+tcnum8*tncoff8
      t24 = tcnum6*tncoff6+tcnum5*tncoff5+tcnum1*tncoff1+tcnum38*tncoff3
     #8+tcnum44*tncoff44+tcnum12*tncoff12+tcnum19*tncoff19+tcnum17*tncof
     #f17+tcnum13*tncoff13+tcnum28*tncoff28+tcnum14*tncoff14
      t37 = tcnum40*tncoff40+tcnum15*tncoff15+tcnum16*tncoff16+tcnum39*t
     #ncoff39+tcnum41*tncoff41+tcnum11*tncoff11+tcnum10*tncoff10+tcnum9*
     #tncoff9+tcnum18*tncoff18+tcnum2*tncoff2+tcnum42*tncoff42
      t50 = tcnum20*tncoff20+tcnum21*tncoff21+tcnum26*tncoff26+tcnum22*t
     #ncoff22+tcnum27*tncoff27+tcnum43*tncoff43+tcnum23*tncoff23+tcnum25
     #*tncoff25+tcnum24*tncoff24+tcnum29*tncoff29+tcnum37*tncoff37+tcnum
     #34*tncoff34
      qqbc21e1int = t12+t24+t37+t50

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc21e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum21e0
      external qqbcnum21e0
      double precision qqbcnum21e1
      external qqbcnum21e1
      double precision qqbcnum21e2
      external qqbcnum21e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t6 = -1/t3
      t9 = 1/l2
      t10 = dlog(l2)
      t13 = -l2+1
      t15 = -1/t13
      t20 = 1/l1
      t21 = dlog(l1)
      tncoff1 = t1/2304-t4*t6/48+t9*t10/48-t13*t15*t9*t6/48+t20*t21/192+
     #t20*t6/96+(-1-l1*l2+l1+l2)/(-1+l1)*t20*t15*t9/96
      t2 = 1/(-1+l2)
      tncoff2 = -t2/(z-1)/48+(l1*l2-l2)/(-1+l1)/l1*t2/l2/96
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t5 = dlog(l2)
      t12 = 1-l1
      tncoff3 = t4*t5/48+t4/(z-1)/48-1/l1*t2*t3/96
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t6 = -1/t3
      t9 = -1+l1
      t10 = 1/t9
      t12 = dlog(-t9)
      t17 = 1/l2
      t18 = dlog(l2)
      t21 = -l2+1
      t23 = -1/t21
      tncoff4 = t1/2304-t4*t6/48-t10*t12/192-t10*t6/96+t17*t18/48-t21*t2
     #3*t17*t6/48+(-l1+l1*l2)*t10/l1*t23*t17/96
      t2 = 1/(-1+l2)
      tncoff5 = -t2/(z-1)/48-1/(-1+l1)*t2/96
      t2 = 1/(-1+l2)
      t3 = 1/l2
      t4 = t2*t3
      t9 = dlog(l2)
      tncoff6 = t4/(z-1)/48+t4*t9/48+t3/(-1+l1)*t2/96
      t1 = -1+l1
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t5 = dlog(l1)
      t9 = dlog(-t1)
      t16 = -l2+1
      tncoff7 = t4*t5/192+t4*t9/192+t4/(z-1)/96-t2*t3/l2/96
      tncoff8 = 1/l1/(-1+l1)/(-1+l2)/96
      tncoff9 = -1/l2/l1/(-1+l1)/(-1+l2)/96
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(z)
      t8 = dlog(-t1)
      tncoff10 = -z*t2*t4/48+t8*z*t2/48-(-z+l2*z)/(-1+l2)/l2*t2/48+(-z*l
     #1+z)/(-1+l1)/l1*t2/96
      tncoff11 = 1/(-1+l2)/(z-1)/48
      tncoff12 = -(l2-l2*z+z)/(-1+l2)/l2/(z-1)/48
      t1 = z-1
      t2 = 1/t1
      t4 = dlog(z)
      t8 = dlog(-t1)
      tncoff13 = -z*t2*t4/48+t8*z*t2/48+z/(-1+l1)*t2/96-(-z+l2*z)/(-1+l2
     #)/l2*t2/48
      tncoff14 = 1/(-1+l2)/(z-1)/48
      tncoff15 = -(l2-l2*z+z)/(-1+l2)/l2/(z-1)/48
      tncoff16 = -z/(-1+l1)/l1/(z-1)/96
      tncoff17 = 1/(z-1)/96-1/l2/96-1/l1/192
      tncoff18 = 1/(-1+l2)/96
      tncoff19 = -1/(-1+l2)/l2/96
      tncoff20 = 1/(z-1)/96+1/(-1+l1)/192-1/l2/96
      tncoff21 = 1/(-1+l2)/96
      tncoff22 = -1/(-1+l2)/l2/96
      tncoff23 = -1/(-1+l1)/l1/192
      tncoff24 = -z/(z-1)/96
      tncoff25 = -z/(z-1)/96
      tncoff26 = -1.D0/192.D0
      tncoff27 = -1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum21e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum21e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum21e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum21e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum21e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum21e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum21e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum21e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum21e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = qqbcnum21e1(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = qqbcnum21e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = qqbcnum21e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = qqbcnum21e2(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = qqbcnum21e2(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff27)

CCC Sum

      t14 = tcnum20*tncoff20+tcnum22*tncoff22+tcnum16*tncoff16+tcnum13*t
     #ncoff13+tcnum14*tncoff14+tcnum18*tncoff18+tcnum23*tncoff23+tcnum7*
     #tncoff7+tcnum21*tncoff21+tcnum5*tncoff5+tcnum27*tncoff27+tcnum25*t
     #ncoff25+tcnum10*tncoff10
      t29 = tcnum1*tncoff1+tcnum11*tncoff11+tcnum8*tncoff8+tcnum9*tncoff
     #9+tcnum4*tncoff4+tcnum19*tncoff19+tcnum12*tncoff12+tcnum17*tncoff1
     #7+tcnum15*tncoff15+tcnum26*tncoff26+tcnum2*tncoff2+tcnum6*tncoff6+
     #tcnum3*tncoff3+tcnum24*tncoff24
      qqbc21e2int = t14+t29

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc21e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum21e0
      external qqbcnum21e0
      double precision qqbcnum21e1
      external qqbcnum21e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t6 = 1-l1
      tncoff1 = 1/(z-1)/96-1/l2/96-1/l1/192
      tncoff2 = 1/(-1+l2)/96
      tncoff3 = -1/(-1+l2)/l2/96
      tncoff4 = 1/(z-1)/96+1/(-1+l1)/192-1/l2/96
      tncoff5 = 1/(-1+l2)/96
      tncoff6 = -1/(-1+l2)/l2/96
      tncoff7 = -1/(-1+l1)/l1/192
      tncoff8 = -z/(z-1)/96
      tncoff9 = -z/(z-1)/96
      tncoff10 = -1.D0/192.D0
      tncoff11 = -1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum21e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum21e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum21e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum21e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff11)

CCC Sum

      qqbc21e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc21e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum21e0
      external qqbcnum21e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -1.D0/192.D0
      tncoff2 = -1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum21e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum21e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbc21e4int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

