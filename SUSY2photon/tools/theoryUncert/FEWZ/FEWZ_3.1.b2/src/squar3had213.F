CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq313e0int(var,ndim,Tcoeff)

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
      double precision sq3num6e0
      external sq3num6e0
      double precision sq3num6e1
      external sq3num6e1
      double precision sq3num6e2
      external sq3num6e2
      double precision sq3num6e3
      external sq3num6e3
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
      t2 = 1-z
      t3 = dlog(t2)
      t4 = 4*t3
      t6 = dlog(1-l4)
      t7 = dlog(l1)
      t8 = dlog(2.D0)
      t9 = 2*t8
      t10 = f2(0.D0,1.D0)
      t11 = dlog(t10)
      t12 = dlog(l3)
      t13 = 2*t12
      t15 = 1/l3
      t18 = -1/t2
      t19 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t20 = 1/t19
      t21 = t18*t20
      t23 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t24 = 1/t23
      t25 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t26 = 1/t25
      t28 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t29 = 1/t28
      t30 = f1(0.D0,1.D0,0.D0)
      t32 = t24*t26*t29*t30
      t35 = t6**2
      t38 = t1**2
      t43 = 0.3141592653589793D1**2
      t46 = t8**2
      t49 = t7**2
      t54 = t11**2
      t81 = t35*t6/49152+t38*t1/49152+t38*t7/16384+t11*t43/24576-t11*t46
     #/4096+t1*t49/16384+t1*t46/4096+t1*t54/16384-t1*t43/24576+t35*t7/16
     #384+t54*t8/8192+t35*t8/8192+t35*t1/16384-t35*t11/16384+t6*t38/1638
     #4+t38*t8/8192-t38*t11/16384+t49*t8/8192-t49*t11/16384+t7*t46/4096
      t102 = t6*t1
      t107 = t1*t7
      t110 = t6*t7
      t117 = t7*t8
      t124 = t6*t8
      t127 = t1*t8
      t130 = t102*t8/4096-t102*t11/8192+t107*t8/4096+t110*t8/4096+t49*t7
     #/49152+t102*t7/8192-t117*t11/4096-t107*t11/8192-t110*t11/8192-t124
     #*t11/4096-t127*t11/4096
      t138 = dlog(l2)
      t139 = 2*t138
      t141 = 1/l2
      t147 = 6*t7*t11
      t148 = 12*t117
      t153 = 6*t110
      t155 = 12*t8*t11
      t162 = 12*t124
      t164 = 6*t6*t11
      t165 = -t147+t148+48*t3*t8-24*t3*t11+t153-t155+24*t6*t3+24*t1*t3+2
     #4*t3*t7+t162-t164
      t166 = 6*t102
      t167 = 6*t107
      t168 = 12*t127
      t170 = 6*t1*t11
      t171 = 3*t49
      t172 = t3**2
      t174 = 3*t54
      t175 = 12*t46
      t176 = 2*t43
      t177 = 3*t38
      t178 = 3*t35
      t179 = t166+t167+t168-t170+t171+48*t172+t174+t175-t176+t177+t178
      t187 = t15*t20
      t193 = -t147+t148+t153-t155+t162-t164+t166+t167+t168-t170+12*t138*
     #t7
      t202 = t138**2
      t204 = 12*t6*t138+12*t1*t138+24*t138*t8-12*t138*t11+t171+t174+t175
     #-t176+t177+t178+12*t202
      t216 = -t147+t148+t153-t155+t162-t164+t166+t167+t168-t170+12*t6*t1
     #2
      t225 = t12**2
      t227 = 24*t12*t8+12*t12*t7-12*t12*t11+12*t1*t12+t171+t174+t175-t17
     #6+t177+t178+12*t225
      tncoff1 = -(t1+t4+t6+t7+t9-t11+t13)*t15*t21*t32/1024+(t81+t7*t54/1
     #6384-t7*t43/24576-t54*t11/49152+t46*t8/6144-t43*t8/12288+t6*t46/40
     #96+t6*t54/16384+t6*t49/16384-t6*t43/24576+z3/3072+t130)*t20*t24*t2
     #6*t29*t30-(t4+t139+t7+t9-t11+t6+t1)*t141*t21*t32/1024+(t165+t179)*
     #t18*t20*t32/12288+(t9-t11+t7+t139+t13+t6+t1)*t141*t187*t32/2048-(t
     #193+t204)*t141*t20*t32/24576+t18*t141*t187*t32/512-(t216+t227)*t15
     #*t20*t32/24576
      t1 = dlog(l4)
      t2 = 1-z
      t3 = dlog(t2)
      t6 = dlog(1-l4)
      t7 = dlog(l3)
      t8 = 2*t7
      t9 = l3/2
      t10 = f2(t9,1.D0)
      t11 = dlog(t10)
      t12 = 2-l3
      t13 = dlog(t12)
      t14 = 2*t13
      t15 = dlog(l1)
      t18 = -1/t12
      t20 = 1/l3
      t23 = -t20/t2
      t25 = Dsquar321(1.D0,l1,0.D0,t9,l4)
      t26 = 1/t25
      t27 = s2vr2(1.D0,0.D0,0.D0,t9,l4)
      t28 = 1/t27
      t30 = s23r2(1.D0,0.D0,0.D0,t9,l4)
      t31 = 1/t30
      t32 = s1vr2(1.D0,0.D0,0.D0,t9,l4)
      t33 = 1/t32
      t35 = f1(t9,1.D0,0.D0)
      t37 = t26*t28*t31*t33*t35
      t40 = dlog(l2)
      t44 = 1/l2
      t50 = t28*t31
      t59 = t15**2
      t61 = t13**2
      t63 = 0.3141592653589793D1**2
      t65 = t11**2
      t67 = t6**2
      t69 = t1**2
      t79 = -6*t1*t11+3*t59+12*t61-2*t63+3*t65+3*t67+3*t69+12*t1*t13+6*t
     #1*t15-12*t11*t13-6*t15*t11
      t100 = t7**2
      t102 = 12*t15*t13-6*t6*t11+12*t6*t13+6*t6*t15+6*t6*t1+12*t1*t7-12*
     #t7*t11+24*t7*t13+12*t7*t15+12*t6*t7+12*t100
      tncoff2 = -(t1+4*t3+t6+t8-t11+t14+t15)*t18*t23*t37/512+(t15-t11+2*
     #t40+t14+t8+t6+t1)*t18*t44*t20*t37/1024+t18*t33*t50*t26*t44*t23*t35
     #/256-(t79+t102)*t18*t20*t26*t50*t33*t35/12288
      t1 = dlog(l1)
      t3 = dlog(1-l2)
      t6 = t3**2
      t8 = dlog(l4)
      t9 = f3(0.D0,1.D0)
      t13 = dlog(1-t9*l1*l2)
      t16 = l1*l2
      t18 = dlog(1-t16)
      t19 = f2(0.D0,1.D0)
      t20 = dlog(t19)
      t24 = dlog(1-l4)
      t27 = dlog(2.D0)
      t38 = dlog(l2)
      t39 = t38**2
      t41 = 6*t1*t3+3*t6+6*t8*t13-6*t18*t20+6*t24*t1+12*t24*t27+6*t24*t3
     #+6*t1*t13+12*t3*t27+6*t24*t18+12*t39
      t66 = 6*t3*t18+6*t24*t13+6*t13*t3+12*t13*t27+12*t38*t3+12*t38*t13+
     #6*t8*t3+24*t38*t27+12*t24*t38-6*t8*t20-6*t13*t20+6*t1*t18
      t90 = 6*t8*t1+6*t13*t18+12*t18*t27+12*t8*t38+12*t8*t27-12*t38*t20-
     #12*t27*t20+12*t38*t1-6*t24*t20-6*t1*t20+6*t8*t18
      t95 = t24**2
      t97 = t8**2
      t99 = t1**2
      t101 = t13**2
      t103 = t20**2
      t109 = 0.3141592653589793D1**2
      t111 = t27**2
      t113 = t18**2
      t115 = 12*t1*t27+6*t24*t8+3*t95+3*t97+3*t99+3*t101+3*t103-6*t3*t20
     #+12*t38*t18-2*t109+12*t111+3*t113
      t118 = 1/l2
      t120 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t121 = 1/t120
      t123 = s2vr2(1.D0,t16,l2,0.D0,l4)
      t125 = s23r2(1.D0,t16,l2,0.D0,l4)
      t128 = s1vr2(1.D0,t16,l2,0.D0,l4)
      t130 = f1(0.D0,1.D0,t16)
      t132 = 1/t123/t125/t128*t130
      t135 = 1-z
      t136 = dlog(t135)
      t138 = 2*t38
      t139 = 2*t27
      t143 = -1/t135
      t148 = dlog(l3)
      t153 = 1/l3*t121
      tncoff3 = (t41+t66+t90+t115)*t118*t121*t132/24576+(4*t136+t138+t1+
     #t139+t13+t3+t18-t20+t24+t8)*t118*t143*t121*t132/1024-(t1+t18+t3+t1
     #3+t139+t138-t20+2*t148+t24+t8)*t118*t153*t132/2048-t143*t118*t153*
     #t132/512
      t1 = dlog(l1)
      t2 = l3/2
      t3 = f2(t2,1.D0)
      t4 = dlog(t3)
      t5 = dlog(l2)
      t8 = dlog(1-l2)
      t9 = 2-l3
      t10 = dlog(t9)
      t12 = f3(t2,1.D0)
      t16 = dlog(1-t12*l1*l2)
      t17 = l1*l2
      t19 = dlog(1-t17)
      t20 = dlog(l3)
      t23 = dlog(1-l4)
      t24 = dlog(l4)
      t27 = -1/t9
      t29 = 1/l2
      t30 = 1/l3
      t33 = Dsquar321(1.D0,l1,l2,t2,l4)
      t34 = 1/t33
      t35 = s2vr2(1.D0,t17,l2,t2,l4)
      t36 = 1/t35
      t38 = s23r2(1.D0,t17,l2,t2,l4)
      t39 = 1/t38
      t40 = s1vr2(1.D0,t17,l2,t2,l4)
      t41 = 1/t40
      t43 = f1(t2,1.D0,t17)
      tncoff4 = -(t1-t4+2*t5+t8+2*t10+t16+t19+2*t20+t23+t24)*t27*t29*t30
     #*t34*t36*t39*t41*t43/1024-t27*t29*t30*t41*t36*t39*t34/(-1+z)*t43/2
     #56
      t1 = dlog(l4)
      t2 = 1-z
      t3 = dlog(t2)
      t4 = 4*t3
      t6 = dlog(1-l4)
      t7 = dlog(z)
      t8 = 2*t7
      t9 = dlog(2.D0)
      t10 = 2*t9
      t11 = dlog(l1)
      t12 = dlog(l3)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = 1/l3
      t20 = -1/t2
      t21 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t22 = 1/t21
      t23 = t20*t22
      t25 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t27 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t30 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t32 = f1(0.D0,z,0.D0)
      t34 = 1/t25/t27/t30*t32
      t37 = dlog(l2)
      t40 = 1/l2
      t47 = 0.3141592653589793D1**2
      t53 = t9**2
      t57 = t3**2
      t63 = t7**2
      t65 = t15**2
      t71 = t6**2
      t73 = 24*t3*t11-2*t47+24*t6*t3+6*t6*t11+12*t53+24*t1*t3+48*t57+12*
     #t11*t9+48*t3*t9+12*t63+3*t65-6*t1*t15-48*t3*t7+3*t71
      t80 = t1**2
      t82 = t11**2
      t104 = -12*t1*t7-12*t6*t7-24*t3*t15+3*t80+3*t82-12*t7*t11-6*t6*t15
     #+12*t1*t9+12*t6*t9+6*t6*t1-12*t15*t9-6*t15*t11-24*t7*t9+12*t7*t15+
     #6*t1*t11
      tncoff5 = (t1+t4+t6-t8+t10+t11+2*t12-t15)*t17*t23*t34/1024+(t4+t10
     #-t8+2*t37+t11-t15+t6+t1)*t40*t23*t34/1024-(t73+t104)*t20*t22*t34/1
     #2288-t20*t40*t17*t22*t34/512
      t1 = dlog(l4)
      t2 = 1-z
      t3 = dlog(t2)
      t6 = dlog(1-l4)
      t7 = dlog(z)
      t9 = l3/2
      t10 = f2(t9,z)
      t11 = dlog(t10)
      t12 = dlog(l3)
      t14 = 2-l3
      t15 = dlog(t14)
      t17 = dlog(l1)
      t20 = -1/t14
      t25 = -1/l3/t2
      t27 = Dsquar321(z,l1,0.D0,t9,l4)
      t28 = 1/t27
      t29 = s2vr2(z,0.D0,0.D0,t9,l4)
      t30 = 1/t29
      t32 = s23r2(z,0.D0,0.D0,t9,l4)
      t33 = 1/t32
      t34 = s1vr2(z,0.D0,0.D0,t9,l4)
      t35 = 1/t34
      t37 = f1(t9,z,0.D0)
      tncoff6 = (t1+4*t3+t6-2*t7-t11+2*t12+2*t15+t17)*t20*t25*t28*t30*t3
     #3*t35*t37/512-t20*t35*t30*t33*t28/l2*t25*t37/256
      t1 = 1-z
      t2 = dlog(t1)
      t4 = f3(0.D0,z)
      t8 = dlog(1-t4*l1*l2)
      t10 = dlog(1-l2)
      t11 = dlog(2.D0)
      t13 = f2(0.D0,z)
      t14 = dlog(t13)
      t15 = l1*l2
      t17 = dlog(1-t15)
      t18 = dlog(l1)
      t19 = dlog(l2)
      t21 = dlog(z)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = 4*t2+t8+t10+2*t11-t14+t17+t18+2*t19-2*t21+t24+t25
      t27 = 1/l2
      t30 = -1/t1
      t31 = Dsquar321(z,l1,l2,0.D0,l4)
      t32 = 1/t31
      t35 = s2vr2(z,t15,l2,0.D0,l4)
      t37 = s23r2(z,t15,l2,0.D0,l4)
      t40 = s1vr2(z,t15,l2,0.D0,l4)
      t42 = f1(0.D0,z,t15)
      t44 = 1/t35/t37/t40*t42
      tncoff7 = -t26*t27*t30*t32*t44/1024+t30*t27/l3*t32*t44/512
      t1 = l3/2
      t2 = l1*l2
      t3 = f1(t1,z,t2)
      t4 = s1vr2(z,t2,l2,t1,l4)
      t7 = s2vr2(z,t2,l2,t1,l4)
      t13 = s23r2(z,t2,l2,t1,l4)
      t15 = Dsquar321(z,l1,l2,t1,l4)
      tncoff8 = t3/t4/t7/(-2+l3)/t13/t15/l2/l3/(-1+z)/256
      t1 = 1/l3
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t9 = 1/t8
      t10 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t11 = 1/t10
      t13 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t14 = 1/t13
      t15 = f1(0.D0,1.D0,0.D0)
      t17 = t9*t11*t14*t15
      t21 = dlog(1-l4)
      t22 = dlog(l4)
      t25 = t21**2
      t27 = t22**2
      t29 = dlog(l1)
      t30 = t29**2
      t32 = f2(0.D0,1.D0)
      t33 = dlog(t32)
      t34 = t33**2
      t36 = dlog(2.D0)
      t37 = t36**2
      t39 = 0.3141592653589793D1**2
      t59 = -t21*t22/8192-t25/16384-t27/16384-t30/16384-t34/16384-t37/40
     #96+t39/24576-t21*t29/8192+t21*t33/8192-t21*t36/4096+t36*t33/4096+t
     #22*t33/8192-t22*t36/4096+t29*t33/8192-t29*t36/4096-t22*t29/8192
      t65 = 1/l2
      t70 = 6*t29
      t71 = 12*t36
      t73 = dlog(-t2)
      t75 = 6*t33
      t76 = 6*t21
      t77 = 6*t22
      t87 = dlog(l2)
      t94 = dlog(l3)
      tncoff9 = t1*t3*t6*t17/1024+t59*t6*t9*t11*t14*t15+t65*t3*t6*t17/10
     #24+(-t70-t71-24*t73+t75-t76-t77)*t3*t6*t17/12288-t65*t1*t6*t17/204
     #8-(-t70-t71+t75-t76-t77-12*t87)*t65*t6*t17/24576-(-t70-t71+t75-t76
     #-t77-12*t94)*t1*t6*t17/24576
      t1 = dlog(l1)
      t3 = 2-l3
      t4 = dlog(t3)
      t6 = l3/2
      t7 = f2(t6,1.D0)
      t8 = dlog(t7)
      t10 = dlog(l3)
      t12 = dlog(l4)
      t15 = dlog(1-l4)
      t19 = -1/t3
      t21 = 1/l3
      t22 = Dsquar321(1.D0,l1,0.D0,t6,l4)
      t23 = 1/t22
      t24 = t21*t23
      t26 = s2vr2(1.D0,0.D0,0.D0,t6,l4)
      t28 = s23r2(1.D0,0.D0,0.D0,t6,l4)
      t31 = s1vr2(1.D0,0.D0,0.D0,t6,l4)
      t33 = f1(t6,1.D0,0.D0)
      t35 = 1/t26/t28/t31*t33
      tncoff10 = -(-6*t1-12*t4+6*t8-12*t10-6*t12-6*t15)*t19*t24*t35/1228
     #8+t19*t21/(-1+z)*t23*t35/512-t19/l2*t24*t35/1024
      t2 = dlog(1-l4)
      t4 = f3(0.D0,1.D0)
      t8 = dlog(1-t4*l1*l2)
      t10 = dlog(2.D0)
      t12 = dlog(l1)
      t15 = dlog(1-l2)
      t17 = f2(0.D0,1.D0)
      t18 = dlog(t17)
      t20 = l1*l2
      t22 = dlog(1-t20)
      t24 = dlog(l4)
      t26 = dlog(l2)
      t29 = 1/l2
      t31 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t32 = 1/t31
      t34 = s2vr2(1.D0,t20,l2,0.D0,l4)
      t36 = s23r2(1.D0,t20,l2,0.D0,l4)
      t39 = s1vr2(1.D0,t20,l2,0.D0,l4)
      t41 = f1(0.D0,1.D0,t20)
      t43 = 1/t34/t36/t39*t41
      tncoff11 = (-6*t2-6*t8-12*t10-6*t12-6*t15+6*t18-6*t22-6*t24-12*t26
     #)*t29*t32*t43/24576-t29/(-1+z)*t32*t43/1024+t29/l3*t32*t43/2048
      t6 = l3/2
      t7 = Dsquar321(1.D0,l1,l2,t6,l4)
      t11 = l1*l2
      t12 = s2vr2(1.D0,t11,l2,t6,l4)
      t14 = s23r2(1.D0,t11,l2,t6,l4)
      t17 = s1vr2(1.D0,t11,l2,t6,l4)
      t19 = f1(t6,1.D0,t11)
      tncoff12 = 1/(-2+l3)/l2/l3/t7/t12/t14/t17*t19/1024
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t10 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t13 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t15 = f1(0.D0,z,0.D0)
      t17 = 1/t8/t10/t13*t15
      t21 = dlog(-t2)
      t24 = dlog(1-l4)
      t26 = dlog(l4)
      t28 = dlog(z)
      t30 = dlog(l1)
      t32 = f2(0.D0,z)
      t33 = dlog(t32)
      t35 = dlog(2.D0)
      tncoff13 = -1/l2*t3*t6*t17/1024-(-24*t21-6*t24-6*t26+12*t28-6*t30+
     #6*t33-12*t35)*t3*t6*t17/12288-1/l3*t3*t6*t17/1024
      t7 = l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t8/t12/t14/t17*t19/512
      t5 = Dsquar321(z,l1,l2,0.D0,l4)
      t8 = l1*l2
      t9 = s2vr2(z,t8,l2,0.D0,l4)
      t11 = s23r2(z,t8,l2,0.D0,l4)
      t14 = s1vr2(z,t8,l2,0.D0,l4)
      t16 = f1(0.D0,z,t8)
      tncoff15 = 1/l2/(-1+z)/t5/t9/t11/t14*t16/1024
      t2 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = dlog(2.D0)
      t24 = f2(0.D0,1.D0)
      t25 = dlog(t24)
      t27 = dlog(l4)
      t30 = dlog(1-l4)
      t32 = dlog(l1)
      tncoff16 = -1/l2*t3*t6*t14/4096-1/l3*t3*t6*t14/4096+(t22/4096-t25/
     #8192+t27/8192+t30/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/2048
      t5 = l3/2
      t6 = Dsquar321(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,0.D0,0.D0,t5,l4)
      t11 = s23r2(1.D0,0.D0,0.D0,t5,l4)
      t14 = s1vr2(1.D0,0.D0,0.D0,t5,l4)
      t16 = f1(t5,1.D0,0.D0)
      tncoff17 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,0.D0,l4)
      t9 = s23r2(1.D0,t5,l2,0.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,0.D0,l4)
      t14 = f1(0.D0,1.D0,t5)
      tncoff18 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t14 = f1(0.D0,z,0.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff20 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num6e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = sq3num6e2(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = sq3num6e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = sq3num6e3(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      sq313e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12
     #+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff
     #16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnco
     #ff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq313e1int(var,ndim,Tcoeff)

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
      double precision sq3num6e0
      external sq3num6e0
      double precision sq3num6e1
      external sq3num6e1
      double precision sq3num6e2
      external sq3num6e2
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
      t6 = t2**2
      t8 = t3**2
      t10 = dlog(l1)
      t11 = t10**2
      t13 = f2(0.D0,1.D0)
      t14 = dlog(t13)
      t15 = t14**2
      t17 = dlog(2.D0)
      t18 = t17**2
      t20 = 0.3141592653589793D1**2
      t40 = -t2*t3/8192-t6/16384-t8/16384-t11/16384-t15/16384-t18/4096+t
     #20/24576-t2*t10/8192+t2*t14/8192-t2*t17/4096+t17*t14/4096+t3*t14/8
     #192-t3*t17/4096+t10*t14/8192-t10*t17/4096-t3*t10/8192
      t41 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t42 = 1/t41
      t44 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t45 = 1/t44
      t47 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t48 = 1/t47
      t49 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t50 = 1/t49
      t52 = f1(0.D0,1.D0,0.D0)
      t55 = 1-z
      t56 = dlog(t55)
      t58 = 2*t17
      t61 = -1/t55
      t66 = t45*t48*t50*t52
      t69 = dlog(l3)
      t72 = 1/l3
      t81 = dlog(l2)
      t84 = 1/l2
      tncoff1 = t40*t42*t45*t48*t50*t52-(4*t56+t3+t10+t58-t14+t2)*t61*t4
     #2*t66/2048+(t3+2*t69+t10+t58-t14+t2)*t72*t42*t66/4096+t72*t61*t42*
     #t66/1024+(t2+2*t81+t10+t58-t14+t3)*t84*t42*t66/4096+t84*t61*t42*t6
     #6/1024-t84*t72*t42*t66/2048
      t1 = dlog(l4)
      t2 = dlog(l3)
      t4 = l3/2
      t5 = f2(t4,1.D0)
      t6 = dlog(t5)
      t7 = dlog(l1)
      t8 = 2-l3
      t9 = dlog(t8)
      t12 = dlog(1-l4)
      t15 = -1/t8
      t17 = 1/l3
      t18 = Dsquar321(1.D0,l1,0.D0,t4,l4)
      t19 = 1/t18
      t20 = t17*t19
      t22 = s2vr2(1.D0,0.D0,0.D0,t4,l4)
      t24 = s23r2(1.D0,0.D0,0.D0,t4,l4)
      t27 = s1vr2(1.D0,0.D0,0.D0,t4,l4)
      t29 = f1(t4,1.D0,0.D0)
      t31 = 1/t22/t24/t27*t29
      tncoff2 = (t1+2*t2-t6+t7+2*t9+t12)*t15*t20*t31/2048+t15*t17/(-1+z)
     #*t19*t31/512-t15/l2*t20*t31/1024
      t2 = dlog(1-l4)
      t3 = dlog(l2)
      t6 = dlog(1-l2)
      t7 = dlog(l1)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,1.D0)
      t14 = dlog(1-t10*l1*l2)
      t15 = l1*l2
      t17 = dlog(1-t15)
      t18 = f2(0.D0,1.D0)
      t19 = dlog(t18)
      t20 = dlog(l4)
      t22 = 1/l2
      t24 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(1.D0,t15,l2,0.D0,l4)
      t29 = s23r2(1.D0,t15,l2,0.D0,l4)
      t32 = s1vr2(1.D0,t15,l2,0.D0,l4)
      t34 = f1(0.D0,1.D0,t15)
      t36 = 1/t27/t29/t32*t34
      tncoff3 = -(t2+2*t3+t6+t7+2*t8+t14+t17-t19+t20)*t22*t25*t36/4096-t
     #22/(-1+z)*t25*t36/1024+t22/l3*t25*t36/2048
      t6 = l3/2
      t7 = Dsquar321(1.D0,l1,l2,t6,l4)
      t11 = l1*l2
      t12 = s2vr2(1.D0,t11,l2,t6,l4)
      t14 = s23r2(1.D0,t11,l2,t6,l4)
      t17 = s1vr2(1.D0,t11,l2,t6,l4)
      t19 = f1(t6,1.D0,t11)
      tncoff4 = 1/(-2+l3)/l2/l3/t7/t12/t14/t17*t19/1024
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l4)
      t5 = dlog(l1)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(z)
      t10 = dlog(2.D0)
      t13 = dlog(1-l4)
      t16 = -1/t1
      t18 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t23 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t26 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t28 = f1(0.D0,z,0.D0)
      t30 = 1/t21/t23/t26*t28
      tncoff5 = (4*t2+t4+t5-t7-2*t8+2*t10+t13)*t16*t19*t30/2048-1/l3*t16
     #*t19*t30/1024-1/l2*t16*t19*t30/1024
      t7 = l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t8/t12/t14/t17*t19/512
      t5 = Dsquar321(z,l1,l2,0.D0,l4)
      t8 = l1*l2
      t9 = s2vr2(z,t8,l2,0.D0,l4)
      t11 = s23r2(z,t8,l2,0.D0,l4)
      t14 = s1vr2(z,t8,l2,0.D0,l4)
      t16 = f1(0.D0,z,t8)
      tncoff7 = 1/l2/(-1+z)/t5/t9/t11/t14*t16/1024
      t2 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = dlog(2.D0)
      t24 = f2(0.D0,1.D0)
      t25 = dlog(t24)
      t27 = dlog(l4)
      t30 = dlog(1-l4)
      t32 = dlog(l1)
      tncoff8 = -1/l2*t3*t6*t14/4096-1/l3*t3*t6*t14/4096+(t22/4096-t25/8
     #192+t27/8192+t30/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/2048
      t5 = l3/2
      t6 = Dsquar321(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,0.D0,0.D0,t5,l4)
      t11 = s23r2(1.D0,0.D0,0.D0,t5,l4)
      t14 = s1vr2(1.D0,0.D0,0.D0,t5,l4)
      t16 = f1(t5,1.D0,0.D0)
      tncoff9 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,0.D0,l4)
      t9 = s23r2(1.D0,t5,l2,0.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,0.D0,l4)
      t14 = f1(0.D0,1.D0,t5)
      tncoff10 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t14 = f1(0.D0,z,0.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff12 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      sq313e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq313e2int(var,ndim,Tcoeff)

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
      double precision sq3num6e0
      external sq3num6e0
      double precision sq3num6e1
      external sq3num6e1
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

      t2 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = dlog(2.D0)
      t24 = f2(0.D0,1.D0)
      t25 = dlog(t24)
      t27 = dlog(l4)
      t30 = dlog(1-l4)
      t32 = dlog(l1)
      tncoff1 = -1/l2*t3*t6*t14/4096-1/l3*t3*t6*t14/4096+(t22/4096-t25/8
     #192+t27/8192+t30/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/2048
      t5 = l3/2
      t6 = Dsquar321(1.D0,l1,0.D0,t5,l4)
      t9 = s2vr2(1.D0,0.D0,0.D0,t5,l4)
      t11 = s23r2(1.D0,0.D0,0.D0,t5,l4)
      t14 = s1vr2(1.D0,0.D0,0.D0,t5,l4)
      t16 = f1(t5,1.D0,0.D0)
      tncoff2 = -1/(-2+l3)/l3/t6/t9/t11/t14*t16/2048
      t2 = Dsquar321(1.D0,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,0.D0,l4)
      t9 = s23r2(1.D0,t5,l2,0.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,0.D0,l4)
      t14 = f1(0.D0,1.D0,t5)
      tncoff3 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t14 = f1(0.D0,z,0.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff5 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      sq313e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq313e3int(var,ndim,Tcoeff)

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
      double precision sq3num6e0
      external sq3num6e0
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

      t1 = Dsquar321(1.D0,l1,0.D0,0.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,0.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,0.D0,l4)
      t11 = f1(0.D0,1.D0,0.D0)
      tncoff1 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq313e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq313e4int(var,ndim,Tcoeff)

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

      sq313e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

