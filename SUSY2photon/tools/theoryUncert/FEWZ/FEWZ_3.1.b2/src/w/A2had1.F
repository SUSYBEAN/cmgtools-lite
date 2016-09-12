CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function A21e0int(var,ndim,Tcoeff)

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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
      double precision A2num1e2
      external A2num1e2
      double precision A2num1e3
      external A2num1e3
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

      t1 = dlog(l3)
      t2 = 2*t1
      t4 = dlog(1-l2)
      t5 = dlog(2.D0)
      t6 = 2*t5
      t7 = dlog(l2)
      t8 = f2(0.D0,1.D0)
      t9 = dlog(t8)
      t10 = 1-z
      t11 = dlog(t10)
      t12 = 4*t11
      t13 = dlog(l4)
      t15 = dlog(1-l4)
      t18 = 1/l3
      t20 = -1/t10
      t23 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t24 = 1/t23
      t25 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t26 = t25**2
      t27 = 1/t26
      t29 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t30 = t29**2
      t31 = 1/t30
      t32 = f1(0.D0,1.D0,0.D0)
      t34 = t24*t27*t31*t32
      t37 = dlog(l1)
      t40 = 1/l1
      t46 = 6*t4*t9
      t47 = t5**2
      t48 = 12*t47
      t49 = t4*t5
      t50 = 12*t49
      t51 = 0.3141592653589793D1**2
      t52 = 2*t51
      t54 = 6*t7*t9
      t56 = 12*t9*t5
      t57 = t15*t4
      t58 = 6*t57
      t59 = t7*t5
      t60 = 12*t59
      t61 = t15*t5
      t62 = 12*t61
      t63 = t15*t13
      t64 = 6*t63
      t65 = t7*t4
      t66 = 6*t65
      t68 = 6*t15*t9
      t69 = t15*t7
      t70 = 6*t69
      t73 = -t46+t48+t50-t52-t54-t56+t58+t60+t62+t64+t66-t68+t70+12*t37*
     #t5
      t75 = 6*t13*t9
      t76 = t13*t4
      t77 = 6*t76
      t78 = t13*t5
      t79 = 12*t78
      t82 = t13*t7
      t83 = 6*t82
      t92 = t4**2
      t93 = 3*t92
      t94 = t9**2
      t95 = 3*t94
      t96 = t13**2
      t97 = 3*t96
      t98 = t15**2
      t99 = 3*t98
      t100 = t7**2
      t101 = 3*t100
      t102 = t37**2
      t104 = -t75+t77+t79+6*t13*t37+t83+6*t37*t4+6*t15*t37+6*t37*t7-6*t3
     #7*t9+t93+t95+t97+t99+t101+3*t102
      t111 = t18*t40
      t115 = -t46+t48+t50-t52-t54-t56+t58+t60+t62+t64+t66-t68+t70-t75
      t126 = t1**2
      t130 = t77+t79+t83+t93+t95+t97+t99+t101+12*t1*t4-12*t1*t9+24*t1*t5
     #+12*t1*t7+12*t15*t1+12*t126+12*t13*t1
      t141 = t11**2
      t143 = -t46+48*t141+t48+t50-t52-t54-t56+t58+t60+t62+t64+t66-t68+t7
     #0
      t156 = -t75+t77+t79+t83+24*t11*t4+24*t11*t7+24*t15*t11+24*t13*t11+
     #t93+t95+t97+t99+t101-24*t11*t9+48*t11*t5
      t191 = -6*t82*t9+12*t78*t4+12*t82*t5+t96*t13+6*t5*t92-2*t4*t51+2*t
     #9*t51+3*t98*t13-2*t7*t51+8*t47*t5-3*t92*t9+3*t98*t7+6*t98*t5+3*t98
     #*t4-3*t98*t9
      t220 = 12*t7*t47+3*t100*t4+6*t100*t5+3*t13*t100+12*t13*t47+3*t13*t
     #94-3*t96*t9+12*t69*t5+t92*t4+t100*t7-t94*t9+t98*t15-6*t65*t9-12*t7
     #8*t9-6*t76*t9+6*t69*t4
      t253 = -6*t69*t9+3*t13*t92-2*t13*t51+3*t7*t94+3*t7*t92-12*t59*t9+6
     #*t63*t4-6*t63*t9+6*t5*t94-4*t51*t5+3*t96*t7+6*t96*t5+16*z3+3*t96*t
     #4-3*t100*t9+3*t15*t96
      t286 = 12*t47*t4-12*t47*t9+3*t4*t94+3*t15*t100+12*t15*t47+3*t15*t9
     #4+3*t15*t92-2*t15*t51+12*t61*t4-12*t61*t9-6*t57*t9+6*t63*t7+12*t63
     #*t5-12*t49*t9+6*t82*t4+12*t59*t4
      tncoff1 = -(t2+t4+t6+t7-t9+t12+t13+t15)*NF*t18*t20*t34/256-NF*(t12
     #+t15+t13+t37+t7+t4-t9+t6)*t40*t20*t34/512-NF*(t73+t104)*t40*t34/12
     #288+NF*t20*t111*t34/256-(t115+t130)*NF*t18*t34/6144+(t15+t13-t9+t4
     #+t7+t37+t6+t2)*NF*t111*t34/1024+NF*(t143+t156)*t20*t34/3072+NF*(t1
     #91+t220+t253+t286)*t24*t27*t31*t32/12288
      t1 = l3-2
      t2 = 1/t1
      t4 = 1/l1
      t5 = 1/l3
      t8 = l3/2
      t9 = s2vr2(1.D0,0.D0,l2,t8,l4)
      t10 = t9**2
      t11 = 1/t10
      t12 = s34r2(1.D0,0.D0,l2,t8,l4)
      t13 = 1/t12
      t14 = t13*t11
      t15 = s1vr2(1.D0,0.D0,l2,t8,l4)
      t16 = t15**2
      t17 = 1/t16
      t18 = -1+z
      t21 = f1(t8,1.D0,0.D0)
      t23 = t14*t17/t18*t21
      t26 = f2(t8,1.D0)
      t27 = dlog(t26)
      t28 = dlog(l3)
      t29 = 2*t28
      t31 = dlog(-t1)
      t32 = 2*t31
      t33 = dlog(l2)
      t35 = dlog(1-l2)
      t37 = dlog(-t18)
      t39 = dlog(l4)
      t41 = dlog(1-l4)
      t44 = t2*t5
      t48 = dlog(l1)
      t59 = t31**2
      t61 = t27**2
      t63 = t28**2
      t65 = t35**2
      t69 = t33**2
      t87 = 12*t59+3*t61+12*t63+3*t65+12*t35*t31+3*t69+6*t33*t35+6*t41*t
     #39+6*t41*t35+6*t41*t33+6*t39*t35+6*t39*t33+12*t39*t31-6*t39*t27
      t112 = t41**2
      t114 = t39**2
      t116 = 0.3141592653589793D1**2
      t118 = 24*t28*t31-6*t41*t27+12*t41*t31+12*t28*t35+12*t28*t33+12*t4
     #1*t28+12*t39*t28-6*t35*t27-6*t33*t27-12*t27*t31+12*t33*t31-12*t28*
     #t27+3*t112+3*t114-2*t116
      tncoff2 = NF*t2*t4*t5*t23/128-NF*(-t27+t29+t32+t33+t35+4*t37+t39+t
     #41)*t44*t23/128+NF*(t41+t39+t48+t35-t27+t33+t32+t29)*t2*t4*t5*t13*
     #t11*t17*t21/512-NF*(t87+t118)*t44*t14*t17*t21/3072
      t1 = -1+z
      t2 = 1/t1
      t5 = 1/l1
      t6 = 1/l3*t5
      t8 = s34r2(1.D0,l1,l2,0.D0,l4)
      t10 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t11 = t10**2
      t14 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t15 = t14**2
      t17 = f1(0.D0,1.D0,l1)
      t19 = 1/t8/t11/t15*t17
      t23 = dlog(1-l4)
      t24 = dlog(2.D0)
      t27 = dlog(l2)
      t30 = dlog(l4)
      t31 = f2(0.D0,1.D0)
      t32 = dlog(t31)
      t36 = dlog(1-l2)
      t43 = dlog(l1)
      t54 = 12*t23*t24+6*t23*t27-6*t30*t32+6*t30*t36+12*t30*t24+6*t30*t2
     #7+12*t43*t24+6*t30*t43+6*t43*t36+6*t23*t43+6*t43*t27
      t58 = dlog(1-l1)
      t61 = f3(0.D0,1.D0)
      t64 = dlog(1-l1*t61)
      t85 = -6*t43*t32+6*t43*t58+6*t43*t64+6*t30*t58+6*t30*t64+6*t23*t58
     #+6*t23*t64-6*t32*t58+12*t64*t24+12*t58*t24+6*t27*t58+6*t36*t64
      t101 = 0.3141592653589793D1**2
      t103 = t24**2
      t109 = 6*t36*t58+6*t64*t58+6*t27*t64-6*t32*t64-6*t23*t32-12*t32*t2
     #4+6*t23*t36-2*t101+12*t103+6*t27*t36-6*t36*t32
      t114 = t58**2
      t116 = t64**2
      t118 = t32**2
      t120 = t36**2
      t122 = t27**2
      t124 = t23**2
      t126 = t30**2
      t128 = t43**2
      t134 = 12*t27*t24-6*t27*t32+3*t114+3*t116+3*t118+3*t120+3*t122+3*t
     #124+3*t126+3*t128+6*t23*t30+12*t36*t24
      t142 = dlog(-t1)
      t144 = 2*t24
      t151 = dlog(l3)
      tncoff3 = -NF*t2*t6*t19/256+NF*(t54+t85+t109+t134)*t5*t19/12288+NF
     #*(4*t142+t23+t30+t43-t32+t36+t144+t58+t27+t64)*t5*t2*t19/512-(t23+
     #t30+t27-t32+t144+t43+t64+t58+t36+2*t151)*NF*t6*t19/1024
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t5 = dlog(1-l1)
      t6 = l3/2
      t7 = f3(t6,1.D0)
      t10 = dlog(1-l1*t7)
      t11 = dlog(l1)
      t12 = -l3+2
      t13 = dlog(t12)
      t15 = f2(t6,1.D0)
      t16 = dlog(t15)
      t18 = dlog(1-l2)
      t19 = dlog(l2)
      t20 = dlog(l3)
      t25 = -1/t12
      t26 = 1/l1
      t29 = 1/l3
      t30 = s34r2(1.D0,l1,l2,t6,l4)
      t31 = 1/t30
      t33 = s2vr2(1.D0,l1,l2,t6,l4)
      t34 = t33**2
      t35 = 1/t34
      t36 = s1vr2(1.D0,l1,l2,t6,l4)
      t37 = t36**2
      t38 = 1/t37
      t40 = f1(t6,1.D0,l1)
      tncoff4 = -NF*(t2+t3+t5+t10+t11+2*t13-t16+t18+t19+2*t20)*t25*t26*t
     #29*t31*t35*t38*t40/512-NF*t25*t26*t29*t35*t31*t38/(-1+z)*t40/128
      t2 = dlog(1-l4)
      t3 = dlog(2.D0)
      t6 = dlog(z)
      t7 = t6**2
      t10 = dlog(1-l2)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t19 = dlog(l2)
      t26 = dlog(l4)
      t31 = 1-z
      t32 = dlog(t31)
      t49 = 12*t2*t3+12*t7-6*t10*t12-24*t6*t3-12*t3*t12-6*t19*t12-6*t2*t
     #12-12*t2*t6-6*t26*t12-12*t26*t6-24*t32*t12-48*t32*t6-12*t10*t6-12*
     #t19*t6+12*t12*t6+48*t32*t3+24*t32*t10+24*t32*t19
      t56 = t10**2
      t58 = t32**2
      t60 = t19**2
      t62 = t12**2
      t64 = 0.3141592653589793D1**2
      t66 = t3**2
      t68 = t2**2
      t70 = t26**2
      t88 = 24*t2*t32+24*t26*t32+6*t2*t10+3*t56+48*t58+3*t60+3*t62-2*t64
     #+12*t66+3*t68+3*t70+6*t26*t19+12*t26*t3+6*t26*t10+6*t2*t19+6*t2*t2
     #6+12*t19*t3+6*t19*t10+12*t10*t3
      t92 = -1/t31
      t94 = s34r2(z,0.D0,l2,0.D0,l4)
      t96 = s2vr2(z,0.D0,l2,0.D0,l4)
      t97 = t96**2
      t100 = s1vr2(z,0.D0,l2,0.D0,l4)
      t101 = t100**2
      t103 = f1(0.D0,z,0.D0)
      t105 = 1/t94/t97/t101*t103
      t109 = 1/l3
      t110 = 1/l1
      t115 = 4*t32
      t116 = 2*t3
      t117 = dlog(l1)
      t118 = 2*t6
      t125 = dlog(l3)
      tncoff5 = -NF*(t49+t88)*t92*t105/3072-NF*t92*t109*t110*t105/256-NF
     #*(-t115-t2-t26-t19-t116-t10-t117+t118+t12)*t110*t92*t105/512-(-t10
     #+t12+t118-2*t125-t116-t19-t115-t26-t2)*NF*t109*t92*t105/256
      t1 = l3-2
      t2 = 1/t1
      t4 = l3/2
      t5 = s2vr2(z,0.D0,l2,t4,l4)
      t6 = t5**2
      t7 = 1/t6
      t8 = s34r2(z,0.D0,l2,t4,l4)
      t9 = 1/t8
      t12 = s1vr2(z,0.D0,l2,t4,l4)
      t13 = t12**2
      t14 = 1/t13
      t17 = 1/l3
      t18 = -1+z
      t19 = 1/t18
      t21 = f1(t4,z,0.D0)
      t25 = dlog(l2)
      t26 = dlog(l3)
      t28 = dlog(z)
      t31 = dlog(1-l2)
      t32 = f2(t4,z)
      t33 = dlog(t32)
      t35 = dlog(-t1)
      t38 = dlog(-t18)
      t40 = dlog(l4)
      t42 = dlog(1-l4)
      tncoff6 = -NF*t2*t7*t9*t14/l1*t17*t19*t21/128-NF*(-t25-2*t26+2*t28
     #-t31+t33-2*t35-4*t38-t40-t42)*t2*t17*t19*t9*t7*t14*t21/128
      t1 = 1-z
      t2 = dlog(t1)
      t5 = dlog(1-l4)
      t6 = dlog(l4)
      t8 = dlog(1-l2)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(z)
      t13 = dlog(l2)
      t14 = dlog(l1)
      t15 = dlog(2.D0)
      t18 = dlog(1-l1)
      t19 = f3(0.D0,z)
      t22 = dlog(1-l1*t19)
      t23 = -4*t2-t5-t6-t8+t10+2*t11-t13-t14-2*t15-t18-t22
      t25 = 1/l1
      t27 = -1/t1
      t30 = s34r2(z,l1,l2,0.D0,l4)
      t32 = s2vr2(z,l1,l2,0.D0,l4)
      t33 = t32**2
      t36 = s1vr2(z,l1,l2,0.D0,l4)
      t37 = t36**2
      t39 = f1(0.D0,z,l1)
      t41 = 1/t30/t33/t37*t39
      tncoff7 = NF*t23*t25*t27*t41/512+NF*t27/l3*t25*t41/256
      t1 = l3/2
      t2 = f1(t1,z,l1)
      t4 = s2vr2(z,l1,l2,t1,l4)
      t5 = t4**2
      t11 = s34r2(z,l1,l2,t1,l4)
      t13 = s1vr2(z,l1,l2,t1,l4)
      t14 = t13**2
      tncoff8 = NF*t2/t5/(l3-2)/t11/t14/l1/l3/(-1+z)/128
      t2 = dlog(1-l4)
      t3 = 6*t2
      t4 = 1-z
      t5 = dlog(t4)
      t7 = dlog(l4)
      t8 = 6*t7
      t9 = dlog(l2)
      t10 = 6*t9
      t11 = dlog(2.D0)
      t12 = 12*t11
      t14 = dlog(1-l2)
      t15 = 6*t14
      t16 = f2(0.D0,1.D0)
      t17 = dlog(t16)
      t18 = 6*t17
      t22 = -1/t4
      t24 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t25 = 1/t24
      t26 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t27 = t26**2
      t28 = 1/t27
      t30 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t31 = t30**2
      t32 = 1/t31
      t33 = f1(0.D0,1.D0,0.D0)
      t35 = t25*t28*t32*t33
      t38 = 1/l1
      t43 = dlog(l1)
      t50 = 1/l3
      t51 = NF*t50
      t55 = dlog(l3)
      t65 = t9**2
      t67 = t11**2
      t69 = t17**2
      t71 = t14**2
      t73 = 0.3141592653589793D1**2
      t75 = t7**2
      t87 = -3*t65-12*t67-3*t69-3*t71+2*t73-3*t75-6*t2*t9-12*t2*t11-12*t
     #9*t11-6*t9*t14+6*t9*t17
      t108 = t2**2
      t110 = -12*t14*t11+12*t17*t11+6*t14*t17-12*t7*t11-6*t7*t14+6*t7*t1
     #7-6*t2*t14+6*t2*t17-6*t2*t7-6*t7*t9-3*t108
      tncoff9 = NF*(-t3-24*t5-t8-t10-t12-t15+t18)*t22*t35/3072+NF*t38*t2
     #2*t35/512-NF*(-t3-t8-t10-t12-t15+t18-6*t43)*t38*t35/12288+t51*t22*
     #t35/256-(-t3-t8-t10-t12-t15+t18-12*t55)*NF*t50*t35/6144-t51*t38*t3
     #5/1024+NF*(t87+t110)*t25*t28*t32*t33/12288
      t1 = l3-2
      t2 = 1/t1
      t3 = NF*t2
      t5 = 1/l3
      t8 = l3/2
      t9 = s34r2(1.D0,0.D0,l2,t8,l4)
      t11 = s2vr2(1.D0,0.D0,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(1.D0,0.D0,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,1.D0,0.D0)
      t20 = 1/t9/t12/t16*t18
      t24 = dlog(-t1)
      t26 = f2(t8,1.D0)
      t27 = dlog(t26)
      t30 = dlog(1-l2)
      t32 = dlog(l2)
      t34 = dlog(l3)
      t37 = dlog(1-l4)
      t39 = dlog(l4)
      tncoff10 = -t3/l1*t5*t20/512-NF*(-12*t24+6*t27-6*t30-6*t32-12*t34-
     #6*t37-6*t39)*t2*t5*t20/3072+t3*t5/(-1+z)*t20/128
      t1 = 1/l1
      t6 = s34r2(1.D0,l1,l2,0.D0,l4)
      t8 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,1.D0,l1)
      t17 = 1/t6/t9/t13*t15
      t25 = dlog(l1)
      t27 = dlog(l4)
      t30 = dlog(1-l4)
      t32 = dlog(l2)
      t35 = dlog(1-l1)
      t37 = f3(0.D0,1.D0)
      t40 = dlog(1-l1*t37)
      t42 = f2(0.D0,1.D0)
      t43 = dlog(t42)
      t45 = dlog(2.D0)
      t48 = dlog(1-l2)
      tncoff11 = -NF*t1/(-1+z)*t17/512+NF/l3*t1*t17/1024+NF*(-6*t25-6*t2
     #7-6*t30-6*t32-6*t35-6*t40+6*t43-12*t45-6*t48)*t1*t17/12288
      t8 = l3/2
      t9 = s34r2(1.D0,l1,l2,t8,l4)
      t11 = s2vr2(1.D0,l1,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(1.D0,l1,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,1.D0,l1)
      tncoff12 = NF/(l3-2)/l1/l3/t9/t12/t16*t18/512
      t1 = dlog(l4)
      t3 = 1-z
      t4 = dlog(t3)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(l2)
      t14 = dlog(1-l2)
      t16 = dlog(z)
      t19 = dlog(1-l4)
      t24 = -1/t3
      t26 = s34r2(z,0.D0,l2,0.D0,l4)
      t28 = s2vr2(z,0.D0,l2,0.D0,l4)
      t29 = t28**2
      t32 = s1vr2(z,0.D0,l2,0.D0,l4)
      t33 = t32**2
      t35 = f1(0.D0,z,0.D0)
      t37 = 1/t26/t29/t33*t35
      tncoff13 = -NF*(-6*t1-24*t4+6*t7-12*t9-6*t11-6*t14+12*t16-6*t19)*t
     #24*t37/3072-NF/l1*t24*t37/512-NF/l3*t24*t37/256
      t9 = l3/2
      t10 = s34r2(z,0.D0,l2,t9,l4)
      t12 = s2vr2(z,0.D0,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(z,0.D0,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,z,0.D0)
      tncoff14 = -NF/(l3-2)/l3/(-1+z)/t10/t13/t17*t19/128
      t6 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff15 = NF/l1/(-1+z)/t6/t9/t13*t15/512
      t3 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t4 = 1/t3
      t6 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,1.D0,0.D0)
      t14 = 1/t7/t10*t13
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t27 = dlog(l2)
      t29 = dlog(2.D0)
      t32 = dlog(1-l2)
      t34 = f2(0.D0,1.D0)
      t35 = dlog(t34)
      tncoff16 = -NF/l1*t4*t14/2048-NF/l3*t4*t14/1024+NF*(6*t22+6*t25+6*
     #t27+12*t29+6*t32-6*t35)*t4*t14/12288+NF/(-1+z)*t4*t14/512
      t6 = l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t10 = t9**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t16 = f1(t6,1.D0,0.D0)
      tncoff17 = -NF/(l3-2)/l3/t7/t10/t14*t16/512
      t3 = s34r2(1.D0,l1,l2,0.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,1.D0,l1)
      tncoff18 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,0.D0,l4)
      t7 = s2vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,0.D0)
      tncoff19 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff20 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = A2num1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = A2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = A2num1e0(var,ndim,z,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = A2num1e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = A2num1e1(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = A2num1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = A2num1e1(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = A2num1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = A2num1e2(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = A2num1e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = A2num1e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      A21e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function A21e1int(var,ndim,Tcoeff)

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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
      double precision A2num1e2
      external A2num1e2
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
      t3 = t2**2
      t6 = dlog(1-l2)
      t7 = f2(0.D0,1.D0)
      t8 = dlog(t7)
      t11 = t8**2
      t13 = dlog(l2)
      t14 = t13**2
      t16 = dlog(l4)
      t17 = t16**2
      t19 = 0.3141592653589793D1**2
      t21 = t6**2
      t23 = dlog(2.D0)
      t32 = 3*t3-6*t6*t8+3*t11+3*t14+3*t17-2*t19+3*t21-12*t8*t23+12*t6*t
     #23-6*t13*t8+6*t13*t6
      t53 = t23**2
      t55 = 12*t13*t23+12*t2*t23+6*t2*t16-6*t2*t8+12*t16*t23+6*t16*t13+6
     #*t2*t6+6*t2*t13-6*t16*t8+6*t16*t6+12*t53
      t58 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t59 = 1/t58
      t61 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t62 = t61**2
      t63 = 1/t62
      t64 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t65 = t64**2
      t66 = 1/t65
      t68 = f1(0.D0,1.D0,0.D0)
      t72 = 2*t23
      t73 = 1-z
      t74 = dlog(t73)
      t79 = -1/t73
      t83 = t59*t63*t66*t68
      t86 = dlog(l3)
      t90 = 1/l3
      t94 = NF*t90
      t98 = dlog(l1)
      t101 = 1/l1
      tncoff1 = -NF*(t32+t55)*t59*t63*t66*t68/12288-NF*(t13+t6-t8+t72+4*
     #t74+t2+t16)*t79*t83/512+(t2+t16+2*t86+t6+t72+t13-t8)*NF*t90*t83/10
     #24+t94*t79*t83/256+NF*(t2+t16+t98-t8+t13+t6+t72)*t101*t83/2048+NF*
     #t101*t79*t83/512-t94*t101*t83/1024
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l3)
      t6 = dlog(l2)
      t8 = dlog(1-l2)
      t9 = l3/2
      t10 = f2(t9,1.D0)
      t11 = dlog(t10)
      t12 = -l3+2
      t13 = dlog(t12)
      t18 = -1/t12
      t19 = 1/l3
      t22 = s34r2(1.D0,0.D0,l2,t9,l4)
      t24 = s2vr2(1.D0,0.D0,l2,t9,l4)
      t25 = t24**2
      t28 = s1vr2(1.D0,0.D0,l2,t9,l4)
      t29 = t28**2
      t31 = f1(t9,1.D0,0.D0)
      t33 = 1/t22/t25/t29*t31
      t36 = NF*t18
      tncoff2 = NF*(t2+t3+2*t4+t6+t8-t11+2*t13)*t18*t19*t33/512+t36*t19/
     #(-1+z)*t33/128-t36/l1*t19*t33/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = f2(0.D0,1.D0)
      t5 = dlog(t4)
      t6 = f3(0.D0,1.D0)
      t9 = dlog(1-l1*t6)
      t10 = dlog(2.D0)
      t13 = dlog(1-l1)
      t15 = dlog(1-l2)
      t16 = dlog(l2)
      t17 = dlog(l1)
      t20 = 1/l1
      t22 = s34r2(1.D0,l1,l2,0.D0,l4)
      t24 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t25 = t24**2
      t28 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t29 = t28**2
      t31 = f1(0.D0,1.D0,l1)
      t33 = 1/t22/t25/t29*t31
      tncoff3 = -NF*(t2+t3-t5+t9+2*t10+t13+t15+t16+t17)*t20*t33/2048-NF*
     #t20/(-1+z)*t33/512+NF/l3*t20*t33/1024
      t8 = l3/2
      t9 = s34r2(1.D0,l1,l2,t8,l4)
      t11 = s2vr2(1.D0,l1,l2,t8,l4)
      t12 = t11**2
      t15 = s1vr2(1.D0,l1,l2,t8,l4)
      t16 = t15**2
      t18 = f1(t8,1.D0,l1)
      tncoff4 = NF/(l3-2)/l1/l3/t9/t12/t16*t18/512
      t1 = dlog(2.D0)
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t5 = dlog(l2)
      t7 = dlog(1-l2)
      t8 = dlog(z)
      t10 = 1-z
      t11 = dlog(t10)
      t14 = dlog(1-l4)
      t15 = dlog(l4)
      t19 = -1/t10
      t21 = s34r2(z,0.D0,l2,0.D0,l4)
      t23 = s2vr2(z,0.D0,l2,0.D0,l4)
      t24 = t23**2
      t27 = s1vr2(z,0.D0,l2,0.D0,l4)
      t28 = t27**2
      t30 = f1(0.D0,z,0.D0)
      t32 = 1/t21/t24/t28*t30
      tncoff5 = -NF*(-2*t1+t4-t5-t7+2*t8-4*t11-t14-t15)*t19*t32/512-NF/l
     #3*t19*t32/256-NF/l1*t19*t32/512
      t9 = l3/2
      t10 = s34r2(z,0.D0,l2,t9,l4)
      t12 = s2vr2(z,0.D0,l2,t9,l4)
      t13 = t12**2
      t16 = s1vr2(z,0.D0,l2,t9,l4)
      t17 = t16**2
      t19 = f1(t9,z,0.D0)
      tncoff6 = -NF/(l3-2)/l3/(-1+z)/t10/t13/t17*t19/128
      t6 = s34r2(z,l1,l2,0.D0,l4)
      t8 = s2vr2(z,l1,l2,0.D0,l4)
      t9 = t8**2
      t12 = s1vr2(z,l1,l2,0.D0,l4)
      t13 = t12**2
      t15 = f1(0.D0,z,l1)
      tncoff7 = NF/l1/(-1+z)/t6/t9/t13*t15/512
      t1 = dlog(l2)
      t3 = dlog(l4)
      t6 = dlog(1-l4)
      t8 = f2(0.D0,1.D0)
      t9 = dlog(t8)
      t12 = dlog(1-l2)
      t14 = dlog(2.D0)
      t18 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t22 = t21**2
      t24 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t25 = t24**2
      t28 = f1(0.D0,1.D0,0.D0)
      t29 = 1/t22/t25*t28
      tncoff8 = -NF*(-6*t1-6*t3-6*t6+6*t9-6*t12-12*t14)*t19*t29/12288+NF
     #/(-1+z)*t19*t29/512-NF/l3*t19*t29/1024-NF/l1*t19*t29/2048
      t6 = l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t10 = t9**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t16 = f1(t6,1.D0,0.D0)
      tncoff9 = -NF/(l3-2)/l3/t7/t10/t14*t16/512
      t3 = s34r2(1.D0,l1,l2,0.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,1.D0,l1)
      tncoff10 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,0.D0,l4)
      t7 = s2vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,0.D0)
      tncoff11 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff12 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,1.D0,ux,l1,l2,l3/2,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = A2num1e0(var,ndim,z,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = A2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.D0,
     #-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = A2num1e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = A2num1e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = A2num1e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      A21e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function A21e2int(var,ndim,Tcoeff)

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
      double precision A2num1e0
      external A2num1e0
      double precision A2num1e1
      external A2num1e1
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
      t3 = dlog(1-l2)
      t4 = f2(0.D0,1.D0)
      t5 = dlog(t4)
      t6 = dlog(l2)
      t7 = dlog(2.D0)
      t10 = dlog(1-l4)
      t13 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t14 = 1/t13
      t16 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t17 = t16**2
      t19 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t20 = t19**2
      t23 = f1(0.D0,1.D0,0.D0)
      t24 = 1/t17/t20*t23
      tncoff1 = NF*(t1+t3-t5+t6+2*t7+t10)*t14*t24/2048+NF/(-1+z)*t14*t24
     #/512-NF/l3*t14*t24/1024-NF/l1*t14*t24/2048
      t6 = l3/2
      t7 = s34r2(1.D0,0.D0,l2,t6,l4)
      t9 = s2vr2(1.D0,0.D0,l2,t6,l4)
      t10 = t9**2
      t13 = s1vr2(1.D0,0.D0,l2,t6,l4)
      t14 = t13**2
      t16 = f1(t6,1.D0,0.D0)
      tncoff2 = -NF/(l3-2)/l3/t7/t10/t14*t16/512
      t3 = s34r2(1.D0,l1,l2,0.D0,l4)
      t6 = s2vr2(1.D0,l1,l2,0.D0,l4)
      t7 = t6**2
      t9 = s1vr2(1.D0,l1,l2,0.D0,l4)
      t10 = t9**2
      t13 = f1(0.D0,1.D0,l1)
      tncoff3 = NF/l1/t3/t7/t10*t13/2048
      t4 = s34r2(z,0.D0,l2,0.D0,l4)
      t7 = s2vr2(z,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t10 = s1vr2(z,0.D0,l2,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,0.D0)
      tncoff4 = -NF/(-1+z)/t4/t8/t11*t14/512
      t1 = s34r2(1.D0,0.D0,l2,0.D0,l4)
      t4 = s2vr2(1.D0,0.D0,l2,0.D0,l4)
      t5 = t4**2
      t7 = s1vr2(1.D0,0.D0,l2,0.D0,l4)
      t8 = t7**2
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff5 = -NF/t1/t5/t8*t11/2048

CCC Numerator functions

      tcnum1 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = A2num1e0(var,ndim,1.D0,ux,0.D0,l2,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = A2num1e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = A2num1e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = A2num1e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      A21e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

