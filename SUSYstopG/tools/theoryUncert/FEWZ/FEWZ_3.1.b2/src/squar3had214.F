CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq314e0int(var,ndim,Tcoeff)

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

      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t5 = 2*t4
      t6 = f2(1.D0,1.D0)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = dlog(2.D0)
      t10 = 2*t9
      t11 = dlog(l3)
      t12 = 2*t11
      t14 = 1/l2
      t16 = 1/l3
      t17 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t18 = 1/t17
      t19 = t16*t18
      t21 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t22 = 1/t21
      t23 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t27 = 1/t26
      t28 = f1(1.D0,1.D0,0.D0)
      t30 = t22*t24*t27*t28
      t33 = -1+z
      t34 = 1/t33
      t39 = 0.3141592653589793D1**2
      t40 = 2*t39
      t41 = t2**2
      t42 = 3*t41
      t43 = t11**2
      t45 = t8**2
      t46 = 3*t45
      t47 = t3**2
      t48 = 3*t47
      t49 = t9**2
      t50 = 12*t49
      t51 = t7**2
      t52 = 3*t51
      t54 = 12*t2*t9
      t56 = 12*t8*t9
      t57 = t8*t7
      t58 = 6*t57
      t60 = 12*t9*t7
      t61 = -t40+t42+12*t43+t46+t48+t50+t52+t54+t56-t58-t60
      t62 = t2*t8
      t63 = 6*t62
      t70 = t2*t7
      t71 = 6*t70
      t76 = t2*t3
      t77 = 6*t76
      t78 = t3*t8
      t79 = 6*t78
      t81 = 12*t3*t9
      t82 = t3*t7
      t83 = 6*t82
      t84 = t63+24*t11*t9-12*t11*t7+12*t11*t8-t71+12*t2*t11+12*t3*t11+t7
     #7+t79+t81-t83
      t129 = z3/3072+t76*t9/4096+t9*t51/8192+t45*t8/49152+t41*t2/49152-t
     #62*t7/8192-t82*t9/4096+t47*t3/49152-t51*t7/49152-t2*t39/24576-t3*t
     #39/24576+t41*t3/16384+t2*t45/16384+t41*t8/16384-t7*t49/4096-t57*t9
     #/4096+t78*t9/4096-t76*t7/8192+t76*t8/8192+t62*t9/4096
      t173 = t8*t49/4096+t3*t49/4096-t78*t7/8192-t70*t9/4096+t49*t9/6144
     #+t2*t49/4096+t3*t45/16384+t3*t51/16384+t2*t51/16384+t2*t47/16384-t
     #39*t9/12288
      t181 = t4**2
      t183 = -t40+t42+12*t181+t46+t48+t50+t52+t54+t56-t58-t60
      t194 = t63-t71+12*t4*t8+12*t2*t4+12*t3*t4+24*t4*t9-12*t4*t7+t77+t7
     #9+t81-t83
      t201 = dlog(-t33)
      t202 = t201**2
      t204 = -t40+t42+t46+t48+t50+t52+48*t202+t54+t56-t58-t60
      t215 = t63-t71+t77+t79+t81-t83+24*t3*t201+24*t201*t8-24*t201*t7+24
     #*t2*t201+48*t201*t9
      t221 = 4*t201
      t224 = t34*t18
      tncoff1 = (t2+t3+t5-t7+t8+t10+t12)*t14*t19*t30/2048+t34*t14*t19*t3
     #0/512-(t61+t84)*t16*t18*t30/24576+(t129+t7*t39/24576+t41*t9/8192-t
     #41*t7/16384+t47*t8/16384+t45*t9/8192-t47*t7/16384-t8*t39/24576+t8*
     #t51/16384+t47*t9/8192-t45*t7/16384+t173)*t18*t22*t24*t27*t28-(t183
     #+t194)*t14*t18*t30/24576+(t204+t215)*t34*t18*t30/12288-(t3+t12+t8+
     #t10-t7+t2+t221)*t16*t224*t30/1024-(t10+t8-t7+t5+t2+t221+t3)*t14*t2
     #24*t30/1024
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t7 = 1-l3/2
      t8 = f2(t7,1.D0)
      t9 = dlog(t8)
      t10 = 2-l3
      t11 = dlog(t10)
      t12 = 2*t11
      t13 = dlog(l1)
      t14 = dlog(l3)
      t15 = 2*t14
      t18 = -1/t10
      t20 = 1/l2
      t21 = 1/l3
      t24 = Dsquar321(1.D0,l1,0.D0,t7,l4)
      t25 = 1/t24
      t26 = s2vr2(1.D0,0.D0,0.D0,t7,l4)
      t27 = 1/t26
      t29 = s23r2(1.D0,0.D0,0.D0,t7,l4)
      t30 = 1/t29
      t31 = s1vr2(1.D0,0.D0,0.D0,t7,l4)
      t32 = 1/t31
      t34 = f1(t7,1.D0,0.D0)
      t36 = t25*t27*t30*t32*t34
      t39 = t2**2
      t61 = 3*t39+12*t2*t14-6*t2*t9+12*t2*t11-12*t11*t9+6*t2*t3+12*t13*t
     #11-6*t13*t9+12*t3*t14+12*t3*t11+6*t2*t13
      t72 = t14**2
      t74 = t9**2
      t76 = t11**2
      t78 = 0.3141592653589793D1**2
      t80 = t13**2
      t82 = t3**2
      t84 = 12*t14*t13-6*t3*t9+6*t3*t13+24*t14*t11-12*t14*t9+12*t72+3*t7
     #4+12*t76-2*t78+3*t80+3*t82
      t89 = t27*t30
      t97 = -1+z
      t99 = t21/t97
      t105 = dlog(-t97)
      tncoff2 = (t2+t3+2*t4-t9+t12+t13+t15)*t18*t20*t21*t36/1024-(t61+t8
     #4)*t18*t21*t25*t89*t32*t34/12288+t18*t32*t89*t25*t20*t99*t34/256-(
     #t3+t15+t12-t9+t13+t2+4*t105)*t18*t99*t36/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = f3(1.D0,1.D0)
      t8 = dlog(1-t4*l1*l2)
      t9 = l1*l2
      t11 = dlog(1-t9)
      t13 = dlog(1-l2)
      t14 = f2(1.D0,1.D0)
      t15 = dlog(t14)
      t16 = dlog(l1)
      t17 = dlog(l2)
      t18 = 2*t17
      t19 = dlog(2.D0)
      t20 = 2*t19
      t21 = dlog(l3)
      t24 = 1/l2
      t27 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t28 = 1/t27
      t29 = 1/l3*t28
      t31 = s2vr2(1.D0,t9,l2,1.D0,l4)
      t33 = s23r2(1.D0,t9,l2,1.D0,l4)
      t36 = s1vr2(1.D0,t9,l2,1.D0,l4)
      t38 = f1(1.D0,1.D0,t9)
      t40 = 1/t31/t33/t36*t38
      t43 = -1+z
      t44 = 1/t43
      t49 = t8**2
      t51 = t2**2
      t53 = t3**2
      t65 = t15**2
      t67 = 0.3141592653589793D1**2
      t71 = 3*t49+3*t51+3*t53+6*t3*t16+12*t3*t19-6*t2*t15+12*t2*t19+6*t2
     #*t16+3*t65-2*t67+6*t2*t8
      t72 = t11**2
      t74 = t17**2
      t78 = t16**2
      t96 = 3*t72+12*t74+6*t2*t3+3*t78+6*t3*t11+12*t16*t19-6*t3*t15-6*t1
     #6*t15+12*t3*t17+24*t17*t19-12*t17*t15+12*t17*t16
      t120 = 12*t2*t17+12*t17*t11+12*t8*t19+6*t11*t8-12*t19*t15+12*t17*t
     #13-6*t8*t15+6*t13*t11+6*t16*t8+6*t16*t13-6*t13*t15
      t141 = t19**2
      t143 = t13**2
      t145 = 6*t16*t11-6*t11*t15+12*t17*t8+12*t13*t19+6*t8*t13+12*t11*t1
     #9+6*t2*t11+6*t3*t13+6*t3*t8+6*t2*t13+12*t141+3*t143
      t153 = dlog(-t43)
      tncoff3 = -(t2+t3+t8+t11+t13-t15+t16+t18+t20+2*t21)*t24*t29*t40/20
     #48-t44*t24*t29*t40/512+(t71+t96+t120+t145)*t24*t28*t40/24576+(t20+
     #t18-t15+t16+t13+t11+t8+t2+4*t153+t3)*t24*t44*t28*t40/1024
      t1 = -2+l3
      t2 = 1/t1
      t3 = 1/l2
      t5 = 1/l3
      t7 = 1-l3/2
      t8 = Dsquar321(1.D0,l1,l2,t7,l4)
      t9 = 1/t8
      t12 = l1*l2
      t13 = s23r2(1.D0,t12,l2,t7,l4)
      t14 = 1/t13
      t15 = s2vr2(1.D0,t12,l2,t7,l4)
      t16 = 1/t15
      t18 = s1vr2(1.D0,t12,l2,t7,l4)
      t19 = 1/t18
      t23 = f1(t7,1.D0,t12)
      t29 = dlog(1-l4)
      t30 = dlog(l4)
      t32 = dlog(-t1)
      t34 = dlog(l1)
      t35 = f3(t7,1.D0)
      t39 = dlog(1-t35*l1*l2)
      t40 = dlog(l2)
      t42 = f2(t7,1.D0)
      t43 = dlog(t42)
      t45 = dlog(1-t12)
      t47 = dlog(1-l2)
      t48 = dlog(l3)
      tncoff4 = -t2*t3*t5*t9*t14*t16*t19/(-1+z)*t23/256-(t29+t30+2*t32+t
     #34+t39+2*t40-t43+t45+t47+2*t48)*t2*t3*t5*t9*t16*t14*t19*t23/1024
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l2
      t5 = 1/l3
      t6 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t7 = 1/t6
      t10 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = f1(1.D0,z,0.D0)
      t19 = 1/t10/t12/t15*t17
      t22 = dlog(l4)
      t23 = dlog(l3)
      t25 = dlog(z)
      t26 = 2*t25
      t27 = dlog(l1)
      t28 = f2(1.D0,z)
      t29 = dlog(t28)
      t30 = dlog(2.D0)
      t31 = 2*t30
      t33 = dlog(1-l4)
      t35 = dlog(-t1)
      t36 = 4*t35
      t39 = t2*t7
      t65 = 0.3141592653589793D1**2
      t69 = t30**2
      t71 = -24*t25*t30+6*t33*t27-6*t22*t29+24*t33*t35+12*t22*t30+6*t22*
     #t27+12*t33*t30+6*t33*t22+12*t27*t30+24*t35*t27+48*t35*t30-2*t65+24
     #*t22*t35+12*t69
      t74 = t35**2
      t76 = t29**2
      t78 = t25**2
      t82 = t33**2
      t84 = t22**2
      t86 = t27**2
      t102 = -48*t35*t25+48*t74+3*t76+12*t78-12*t27*t25+3*t82+3*t84+3*t8
     #6-6*t33*t29-12*t22*t25-24*t35*t29+12*t29*t25-12*t33*t25-6*t27*t29-
     #12*t29*t30
      t108 = dlog(l2)
      tncoff5 = -t2*t3*t5*t7*t19/512+(t22+2*t23-t26+t27-t29+t31+t33+t36)
     #*t5*t39*t19/1024-(t71+t102)*t2*t7*t19/12288+(2*t108-t29-t26+t27+t3
     #1+t33+t36+t22)*t3*t39*t19/1024
      t1 = -2+l3
      t2 = 1/t1
      t4 = 1-l3/2
      t5 = s1vr2(z,0.D0,0.D0,t4,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,t4,l4)
      t9 = 1/t8
      t10 = s23r2(z,0.D0,0.D0,t4,l4)
      t11 = 1/t10
      t14 = Dsquar321(z,l1,0.D0,t4,l4)
      t15 = 1/t14
      t19 = -1+z
      t21 = 1/l3/t19
      t22 = f1(t4,z,0.D0)
      t27 = dlog(l4)
      t28 = dlog(l1)
      t30 = dlog(-t1)
      t32 = f2(t4,z)
      t33 = dlog(t32)
      t34 = dlog(l3)
      t36 = dlog(z)
      t39 = dlog(1-l4)
      t41 = dlog(-t19)
      tncoff6 = -t2*t6*t9*t11*t15/l2*t21*t22/256+(t27+t28+2*t30-t33+2*t3
     #4-2*t36+t39+4*t41)*t2*t21*t15*t9*t11*t6*t22/512
      t1 = -1+z
      t2 = 1/t1
      t3 = 1/l2
      t6 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = 1/t6
      t10 = l1*l2
      t11 = s2vr2(z,t10,l2,1.D0,l4)
      t13 = s23r2(z,t10,l2,1.D0,l4)
      t16 = s1vr2(z,t10,l2,1.D0,l4)
      t18 = f1(1.D0,z,t10)
      t20 = 1/t11/t13/t16*t18
      t23 = dlog(l2)
      t25 = f3(1.D0,z)
      t29 = dlog(1-t25*l1*l2)
      t31 = dlog(1-l2)
      t33 = dlog(1-t10)
      t34 = dlog(z)
      t36 = f2(1.D0,z)
      t37 = dlog(t36)
      t38 = dlog(l1)
      t39 = dlog(2.D0)
      t42 = dlog(1-l4)
      t44 = dlog(-t1)
      t46 = dlog(l4)
      t47 = 2*t23+t29+t31+t33-2*t34-t37+t38+2*t39+t42+4*t44+t46
      tncoff7 = t2*t3/l3*t7*t20/512-t47*t3*t2*t7*t20/1024
      t2 = 1-l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t5 = Dsquar321(z,l1,l2,t2,l4)
      t8 = s23r2(z,t3,l2,t2,l4)
      t14 = s2vr2(z,t3,l2,t2,l4)
      t16 = s1vr2(z,t3,l2,t2,l4)
      tncoff8 = t4/t5/t8/(-2+l3)/t14/t16/l2/l3/(-1+z)/256
      t2 = dlog(1-l4)
      t3 = 6*t2
      t4 = dlog(l4)
      t5 = 6*t4
      t6 = dlog(l3)
      t8 = dlog(l1)
      t9 = 6*t8
      t10 = dlog(2.D0)
      t11 = 12*t10
      t12 = f2(1.D0,1.D0)
      t13 = dlog(t12)
      t14 = 6*t13
      t16 = 1/l3
      t18 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t22 = 1/t21
      t23 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t27 = 1/t26
      t28 = f1(1.D0,1.D0,0.D0)
      t30 = t22*t24*t27*t28
      t37 = t10**2
      t39 = 0.3141592653589793D1**2
      t41 = t4**2
      t43 = t8**2
      t45 = t2**2
      t51 = t13**2
      t65 = -t2*t10/4096-t2*t8/8192-t37/4096+t39/24576-t41/16384-t43/163
     #84-t45/16384+t8*t13/8192+t10*t13/4096-t51/16384+t2*t13/8192-t2*t4/
     #8192+t4*t13/8192-t4*t10/4096-t8*t10/4096-t4*t8/8192
      t71 = dlog(l2)
      t74 = 1/l2
      t79 = 1-z
      t80 = dlog(t79)
      t84 = -1/t79
      tncoff9 = -(-t3-t5-12*t6-t9-t11+t14)*t16*t19*t30/24576+t65*t19*t22
     #*t24*t27*t28-(-t3-t5-t9-t11+t14-12*t71)*t74*t19*t30/24576+(-t3-t5-
     #t9-t11+t14-24*t80)*t84*t19*t30/12288+t16*t84*t19*t30/1024+t74*t84*
     #t19*t30/1024-t74*t16*t19*t30/2048
      t1 = -2+l3
      t2 = 1/t1
      t5 = 1/l3
      t7 = 1-l3/2
      t8 = Dsquar321(1.D0,l1,0.D0,t7,l4)
      t9 = 1/t8
      t10 = t5*t9
      t12 = s2vr2(1.D0,0.D0,0.D0,t7,l4)
      t14 = s23r2(1.D0,0.D0,0.D0,t7,l4)
      t17 = s1vr2(1.D0,0.D0,0.D0,t7,l4)
      t19 = f1(t7,1.D0,0.D0)
      t21 = 1/t12/t14/t17*t19
      t24 = dlog(l1)
      t27 = dlog(-t1)
      t29 = f2(t7,1.D0)
      t30 = dlog(t29)
      t32 = dlog(l3)
      t35 = dlog(1-l4)
      t37 = dlog(l4)
      tncoff10 = -t2/l2*t10*t21/1024-(-6*t24-12*t27+6*t30-12*t32-6*t35-6
     #*t37)*t2*t10*t21/12288+t2*t5/(-1+z)*t9*t21/512
      t1 = 1/l2
      t4 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t5 = 1/t4
      t7 = l1*l2
      t8 = s2vr2(1.D0,t7,l2,1.D0,l4)
      t10 = s23r2(1.D0,t7,l2,1.D0,l4)
      t13 = s1vr2(1.D0,t7,l2,1.D0,l4)
      t15 = f1(1.D0,1.D0,t7)
      t17 = 1/t8/t10/t13*t15
      t21 = dlog(1-l2)
      t23 = dlog(2.D0)
      t25 = dlog(l1)
      t27 = f2(1.D0,1.D0)
      t28 = dlog(t27)
      t30 = f3(1.D0,1.D0)
      t34 = dlog(1-t30*l1*l2)
      t37 = dlog(1-t7)
      t40 = dlog(1-l4)
      t42 = dlog(l4)
      t44 = dlog(l2)
      tncoff11 = t1/l3*t5*t17/2048+(-6*t21-12*t23-6*t25+6*t28-6*t34-6*t3
     #7-6*t40-6*t42-12*t44)*t1*t5*t17/24576-t1/(-1+z)*t5*t17/1024
      t7 = 1-l3/2
      t8 = Dsquar321(1.D0,l1,l2,t7,l4)
      t12 = l1*l2
      t13 = s2vr2(1.D0,t12,l2,t7,l4)
      t15 = s23r2(1.D0,t12,l2,t7,l4)
      t18 = s1vr2(1.D0,t12,l2,t7,l4)
      t20 = f1(t7,1.D0,t12)
      tncoff12 = 1/(-2+l3)/l2/l3/t8/t13/t15/t18*t20/1024
      t2 = -1+z
      t3 = 1/t2
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t10 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t13 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t15 = f1(1.D0,z,0.D0)
      t17 = 1/t8/t10/t13*t15
      t20 = dlog(2.D0)
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t28 = dlog(-t2)
      t30 = dlog(l1)
      t32 = dlog(z)
      t34 = f2(1.D0,z)
      t35 = dlog(t34)
      tncoff13 = -1/l3*t3*t6*t17/1024-(-12*t20-6*t22-6*t25-24*t28-6*t30+
     #12*t32+6*t35)*t3*t6*t17/12288-1/l2*t3*t6*t17/1024
      t8 = 1-l3/2
      t9 = Dsquar321(z,l1,0.D0,t8,l4)
      t13 = s2vr2(z,0.D0,0.D0,t8,l4)
      t15 = s23r2(z,0.D0,0.D0,t8,l4)
      t18 = s1vr2(z,0.D0,0.D0,t8,l4)
      t20 = f1(t8,z,0.D0)
      tncoff14 = -1/(-2+l3)/l3/(-1+z)/t9/t13/t15/t18*t20/512
      t5 = Dsquar321(z,l1,l2,1.D0,l4)
      t8 = l1*l2
      t9 = s2vr2(z,t8,l2,1.D0,l4)
      t11 = s23r2(z,t8,l2,1.D0,l4)
      t14 = s1vr2(z,t8,l2,1.D0,l4)
      t16 = f1(1.D0,z,t8)
      tncoff15 = 1/l2/(-1+z)/t5/t9/t11/t14*t16/1024
      t2 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = f2(1.D0,1.D0)
      t23 = dlog(t22)
      t25 = dlog(2.D0)
      t27 = dlog(l1)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff16 = -1/l3*t3*t6*t14/4096-1/l2*t3*t6*t14/4096+(-t23/8192+t25
     #/4096+t27/8192+t29/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/204
     #8
      t6 = 1-l3/2
      t7 = Dsquar321(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,0.D0,0.D0,t6,l4)
      t12 = s23r2(1.D0,0.D0,0.D0,t6,l4)
      t15 = s1vr2(1.D0,0.D0,0.D0,t6,l4)
      t17 = f1(t6,1.D0,0.D0)
      tncoff17 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,1.D0,l4)
      t9 = s23r2(1.D0,t5,l2,1.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,1.D0,l4)
      t14 = f1(1.D0,1.D0,t5)
      tncoff18 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t14 = f1(1.D0,z,0.D0)
      tncoff19 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff20 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = sq3num6e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #l8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = sq3num6e2(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = sq3num6e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = sq3num6e3(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      sq314e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12
     #+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff
     #16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnco
     #ff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq314e1int(var,ndim,Tcoeff)

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
      t3 = dlog(2.D0)
      t6 = dlog(l1)
      t9 = t3**2
      t11 = 0.3141592653589793D1**2
      t13 = dlog(l4)
      t14 = t13**2
      t16 = t6**2
      t18 = t2**2
      t20 = f2(1.D0,1.D0)
      t21 = dlog(t20)
      t26 = t21**2
      t40 = -t2*t3/4096-t2*t6/8192-t9/4096+t11/24576-t14/16384-t16/16384
     #-t18/16384+t6*t21/8192+t3*t21/4096-t26/16384+t2*t21/8192-t2*t13/81
     #92+t13*t21/8192-t13*t3/4096-t6*t3/4096-t13*t6/8192
      t41 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t42 = 1/t41
      t44 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t45 = 1/t44
      t47 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t48 = 1/t47
      t49 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t50 = 1/t49
      t52 = f1(1.D0,1.D0,0.D0)
      t55 = 2*t3
      t56 = 1-z
      t57 = dlog(t56)
      t61 = -1/t56
      t66 = t45*t48*t50*t52
      t69 = dlog(l3)
      t72 = 1/l3
      t81 = dlog(l2)
      t84 = 1/l2
      tncoff1 = t40*t42*t45*t48*t50*t52-(t13+t6+t55-t21+t2+4*t57)*t61*t4
     #2*t66/2048+(t2+t13+t6+t55-t21+2*t69)*t72*t42*t66/4096+t72*t61*t42*
     #t66/1024+(t2+t13+2*t81+t6+t55-t21)*t84*t42*t66/4096+t84*t61*t42*t6
     #6/1024-t84*t72*t42*t66/2048
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t5 = 1-l3/2
      t6 = f2(t5,1.D0)
      t7 = dlog(t6)
      t8 = dlog(l1)
      t9 = 2-l3
      t10 = dlog(t9)
      t12 = dlog(l3)
      t16 = -1/t9
      t18 = 1/l3
      t19 = Dsquar321(1.D0,l1,0.D0,t5,l4)
      t20 = 1/t19
      t21 = t18*t20
      t23 = s2vr2(1.D0,0.D0,0.D0,t5,l4)
      t25 = s23r2(1.D0,0.D0,0.D0,t5,l4)
      t28 = s1vr2(1.D0,0.D0,0.D0,t5,l4)
      t30 = f1(t5,1.D0,0.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff2 = (t2+t3-t7+t8+2*t10+2*t12)*t16*t21*t32/2048+t16*t18/(-1+z
     #)*t20*t32/512-t16/l2*t21*t32/1024
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t7 = dlog(1-l2)
      t8 = dlog(l1)
      t9 = f3(1.D0,1.D0)
      t13 = dlog(1-t9*l1*l2)
      t14 = dlog(2.D0)
      t16 = l1*l2
      t18 = dlog(1-t16)
      t19 = f2(1.D0,1.D0)
      t20 = dlog(t19)
      t22 = 1/l2
      t24 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t25 = 1/t24
      t27 = s2vr2(1.D0,t16,l2,1.D0,l4)
      t29 = s23r2(1.D0,t16,l2,1.D0,l4)
      t32 = s1vr2(1.D0,t16,l2,1.D0,l4)
      t34 = f1(1.D0,1.D0,t16)
      t36 = 1/t27/t29/t32*t34
      tncoff3 = -(t2+t3+2*t4+t7+t8+t13+2*t14+t18-t20)*t22*t25*t36/4096-t
     #22/(-1+z)*t25*t36/1024+t22/l3*t25*t36/2048
      t7 = 1-l3/2
      t8 = Dsquar321(1.D0,l1,l2,t7,l4)
      t12 = l1*l2
      t13 = s2vr2(1.D0,t12,l2,t7,l4)
      t15 = s23r2(1.D0,t12,l2,t7,l4)
      t18 = s1vr2(1.D0,t12,l2,t7,l4)
      t20 = f1(t7,1.D0,t12)
      tncoff4 = 1/(-2+l3)/l2/l3/t8/t13/t15/t18*t20/1024
      t1 = dlog(l4)
      t2 = dlog(l1)
      t3 = f2(1.D0,z)
      t4 = dlog(t3)
      t5 = dlog(z)
      t7 = dlog(2.D0)
      t10 = dlog(1-l4)
      t11 = 1-z
      t12 = dlog(t11)
      t16 = -1/t11
      t18 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t19 = 1/t18
      t21 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t23 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t26 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t28 = f1(1.D0,z,0.D0)
      t30 = 1/t21/t23/t26*t28
      tncoff5 = (t1+t2-t4-2*t5+2*t7+t10+4*t12)*t16*t19*t30/2048-1/l3*t16
     #*t19*t30/1024-1/l2*t16*t19*t30/1024
      t8 = 1-l3/2
      t9 = Dsquar321(z,l1,0.D0,t8,l4)
      t13 = s2vr2(z,0.D0,0.D0,t8,l4)
      t15 = s23r2(z,0.D0,0.D0,t8,l4)
      t18 = s1vr2(z,0.D0,0.D0,t8,l4)
      t20 = f1(t8,z,0.D0)
      tncoff6 = -1/(-2+l3)/l3/(-1+z)/t9/t13/t15/t18*t20/512
      t5 = Dsquar321(z,l1,l2,1.D0,l4)
      t8 = l1*l2
      t9 = s2vr2(z,t8,l2,1.D0,l4)
      t11 = s23r2(z,t8,l2,1.D0,l4)
      t14 = s1vr2(z,t8,l2,1.D0,l4)
      t16 = f1(1.D0,z,t8)
      tncoff7 = 1/l2/(-1+z)/t5/t9/t11/t14*t16/1024
      t2 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = f2(1.D0,1.D0)
      t23 = dlog(t22)
      t25 = dlog(2.D0)
      t27 = dlog(l1)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff8 = -1/l3*t3*t6*t14/4096-1/l2*t3*t6*t14/4096+(-t23/8192+t25/
     #4096+t27/8192+t29/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/2048
      t6 = 1-l3/2
      t7 = Dsquar321(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,0.D0,0.D0,t6,l4)
      t12 = s23r2(1.D0,0.D0,0.D0,t6,l4)
      t15 = s1vr2(1.D0,0.D0,0.D0,t6,l4)
      t17 = f1(t6,1.D0,0.D0)
      tncoff9 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,1.D0,l4)
      t9 = s23r2(1.D0,t5,l2,1.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,1.D0,l4)
      t14 = f1(1.D0,1.D0,t5)
      tncoff10 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t14 = f1(1.D0,z,0.D0)
      tncoff11 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff12 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num6e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = sq3num6e1(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = sq3num6e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = sq3num6e2(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      sq314e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq314e2int(var,ndim,Tcoeff)

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

      t2 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = 1/t2
      t5 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,1.D0,0.D0)
      t14 = 1/t8/t10*t13
      t22 = f2(1.D0,1.D0)
      t23 = dlog(t22)
      t25 = dlog(2.D0)
      t27 = dlog(l1)
      t29 = dlog(l4)
      t32 = dlog(1-l4)
      tncoff1 = -1/l3*t3*t6*t14/4096-1/l2*t3*t6*t14/4096+(-t23/8192+t25/
     #4096+t27/8192+t29/8192+t32/8192)*t3*t6*t14+1/(-1+z)*t3*t6*t14/2048
      t6 = 1-l3/2
      t7 = Dsquar321(1.D0,l1,0.D0,t6,l4)
      t10 = s2vr2(1.D0,0.D0,0.D0,t6,l4)
      t12 = s23r2(1.D0,0.D0,0.D0,t6,l4)
      t15 = s1vr2(1.D0,0.D0,0.D0,t6,l4)
      t17 = f1(t6,1.D0,0.D0)
      tncoff2 = -1/(-2+l3)/l3/t7/t10/t12/t15*t17/2048
      t2 = Dsquar321(1.D0,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(1.D0,t5,l2,1.D0,l4)
      t9 = s23r2(1.D0,t5,l2,1.D0,l4)
      t11 = s1vr2(1.D0,t5,l2,1.D0,l4)
      t14 = f1(1.D0,1.D0,t5)
      tncoff3 = 1/l2/t2/t6/t9/t11*t14/4096
      t3 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t11 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t14 = f1(1.D0,z,0.D0)
      tncoff4 = -1/(-1+z)/t3/t6/t9/t11*t14/2048
      t1 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff5 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num6e0(var,ndim,1.D0,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num6e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num6e1(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      sq314e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq314e3int(var,ndim,Tcoeff)

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

      t1 = Dsquar321(1.D0,l1,0.D0,1.D0,l4)
      t3 = s2vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(1.D0,0.D0,0.D0,1.D0,l4)
      t8 = s1vr2(1.D0,0.D0,0.D0,1.D0,l4)
      t11 = f1(1.D0,1.D0,0.D0)
      tncoff1 = -1/t1/t3/t6/t8*t11/8192

CCC Numerator functions

      tcnum1 = sq3num6e0(var,ndim,1.D0,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq314e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq314e4int(var,ndim,Tcoeff)

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

      sq314e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

