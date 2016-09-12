CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function CD1e0int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
      double precision CDnum1e1
      external CDnum1e1
      double precision CDnum1e2
      external CDnum1e2
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

      t1 = -1+z
      t2 = 1/l8
      t5 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t7 = 1/l3*t6
      t9 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t10 = 1/t9
      t11 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t12 = 1/t11
      t14 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t15 = 1/t14
      t16 = f1(0.D0,z,1.D0)
      t18 = t10*t12*t15*t16
      t21 = dlog(l2)
      t23 = dlog(l1)
      t25 = dlog(-t1)
      t26 = 4*t25
      t27 = f3(0.D0,z)
      t29 = dlog(1-t27)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t32 = dlog(2.D0)
      t33 = 2*t32
      t34 = dlog(z)
      t35 = 2*t34
      t36 = dlog(l4)
      t38 = dlog(1-l4)
      t41 = 1/l2
      t46 = dlog(l8)
      t59 = -t23-t26-t29+t31-t33+t35-t36-t38-t46-2*t34*l8+t29*l8-t31*l8+
     #t36*l8+4*t25*l8+t46*l8+t23*l8+2*t32*l8+t38*l8
      t71 = t1*t41
      t107 = t29**2
      t109 = t23**2
      t111 = t25**2
      t113 = t38**2
      t115 = -6*t36*t31+6*t36*t29-12*t38*t34+24*t38*t25-12*t29*t34+12*t3
     #8*t32-12*t31*t32+6*t38*t29+6*t38*t23+12*t23*t32+6*t23*t29+6*t36*t2
     #3-24*t34*t32-6*t23*t31+3*t107+3*t109+48*t111+3*t113
      t116 = t36**2
      t118 = t32**2
      t122 = 0.3141592653589793D1**2
      t123 = t31**2
      t125 = t34**2
      t153 = 3*t116+12*t118-24*t25*t31-t122+3*t123+12*t125+24*t25*t29-6*
     #t29*t31-6*t38*t31+12*t31*t34-48*t25*t34+48*t25*t32+24*t36*t25-12*t
     #36*t34-12*t23*t34+24*t23*t25+12*t29*t32+6*t38*t36+12*t36*t32
      t159 = dlog(l3)
      tncoff1 = -t1*t2*t7*t18/1024+t1*(2*t21+t23+t26+t29-t31+t33-t35+t36
     #+t38)*t41*t6*t18/1024+t1*t59/(-1+l8)*t2*t6*t10*t12*t15*t16/2048-t7
     #1*t2*t6*t18/1024-t71*t7*t18/512-t1*(t115+t153)*t6*t18/12288+(t38+t
     #36+2*t159+t23+t33+t26+t29-t31-t35)*t1*t7*t18/1024
      t1 = -1+z
      t2 = 1/l3
      t4 = -1+l8
      t5 = 1/t4
      t6 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t7 = 1/t6
      t10 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t11 = 1/t10
      t12 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t13 = 1/t12
      t15 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t16 = 1/t15
      t17 = f1(0.D0,z,1.D0)
      t19 = t11*t13*t16*t17
      t22 = dlog(l2)
      t24 = dlog(l1)
      t26 = dlog(-t1)
      t27 = 4*t26
      t28 = f3(0.D0,z)
      t30 = dlog(1-t28)
      t31 = f2(0.D0,z)
      t32 = dlog(t31)
      t33 = dlog(2.D0)
      t34 = 2*t33
      t35 = dlog(z)
      t36 = 2*t35
      t37 = dlog(l4)
      t39 = dlog(1-l4)
      t42 = 1/l2
      t43 = t42*t7
      t51 = dlog(-t4)
      t77 = t2*t7
      t81 = dlog(l3)
      t92 = 0.3141592653589793D1**2
      t93 = t33**2
      t95 = t32**2
      t97 = t35**2
      t123 = -12*t24*t35-12*t37*t35-t92+12*t93+3*t95+12*t97+24*t37*t26+1
     #2*t30*t33+6*t39*t37+12*t37*t33+6*t37*t24-6*t37*t32+6*t37*t30-12*t3
     #9*t35+24*t39*t26-12*t30*t35+12*t39*t33-12*t32*t33
      t136 = t30**2
      t138 = t24**2
      t140 = t37**2
      t142 = t26**2
      t144 = t39**2
      t162 = 6*t39*t30+6*t39*t24+12*t24*t33+6*t24*t30-24*t35*t33-6*t24*t
     #32+3*t136+3*t138+3*t140+48*t142+3*t144-24*t26*t32+24*t26*t30-6*t30
     #*t32-6*t39*t32+12*t32*t35-48*t26*t35+48*t26*t33+24*t24*t26
      tncoff2 = t1*t2*t5*t7*t19/1024+t1*(2*t22+t24+t27+t30-t32+t34-t36+t
     #37+t39)*t43*t19/1024+t1*(2*t35*l8-t37*l8-t51*l8-t39*l8-2*t33*l8-t3
     #0*l8+t32*l8-t24*l8-4*t26*l8)*t5/l8*t7*t11*t13*t16*t17/2048+t1*t5*t
     #43*t19/1024-t1*t42*t77*t19/512+(t39+t37+2*t81+t24+t34+t27+t30-t32-
     #t36)*t1*t77*t19/1024-t1*(t123+t162)*t7*t19/12288
      t1 = -1+z
      t2 = 1/l8
      t5 = -1+l8
      t6 = 1/t5
      t9 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t11 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t14 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t16 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t19 = f1(0.D0,z,1.D0)
      t21 = 1/t9/t11/t14/t16*t19
      t24 = dlog(z)
      t26 = dlog(2.D0)
      t28 = f2(0.D0,z)
      t29 = dlog(t28)
      t31 = dlog(-t5)
      t32 = f3(0.D0,z)
      t34 = dlog(1-t32)
      t36 = dlog(-t1)
      t38 = dlog(l1)
      t40 = dlog(1-l4)
      t41 = dlog(l8)
      t42 = dlog(l4)
      tncoff3 = -t1*t2/l3*t6*t21/1024+t1*(-2*t24+2*t26-t29+t31+t34+4*t36
     #+t38+t40+t41+t42)*t6*t2*t21/2048-t1*t6/l2*t2*t21/1024
      t1 = -1+z
      t2 = -2+l3
      t3 = 1/t2
      t4 = t1*t3
      t5 = 1/l3
      t9 = l3/2
      t10 = Dsquar311(z,l1,0.D0,t9,l4)
      t12 = s2vr2(z,1.D0,1.D0,t9,l4)
      t15 = s23r2(z,1.D0,1.D0,t9,l4)
      t17 = s1vr2(z,1.D0,1.D0,t9,l4)
      t20 = f1(t9,z,1.D0)
      t22 = 1/t10/t12/t15/t17*t20
      t31 = dlog(1-l4)
      t32 = dlog(l4)
      t33 = dlog(l3)
      t36 = dlog(-t1)
      t38 = f3(t9,z)
      t40 = dlog(1-t38)
      t41 = dlog(l1)
      t43 = dlog(-t2)
      t45 = dlog(z)
      t47 = f2(t9,z)
      t48 = dlog(t47)
      tncoff4 = -t4*t5/l8*t22/512-t4/l2*t5*t22/256+t1*(t31+t32+2*t33+4*t
     #36+t40+t41+2*t43-2*t45-t48)*t3*t5*t22/512
      t1 = -1+z
      t5 = -2+l3
      t6 = 1/t5
      t7 = 1/l3
      t8 = t6*t7
      t10 = l3/2
      t11 = Dsquar311(z,l1,0.D0,t10,l4)
      t13 = s2vr2(z,1.D0,1.D0,t10,l4)
      t16 = s23r2(z,1.D0,1.D0,t10,l4)
      t18 = s1vr2(z,1.D0,1.D0,t10,l4)
      t21 = f1(t10,z,1.D0)
      t23 = 1/t11/t13/t16/t18*t21
      t33 = dlog(1-l4)
      t34 = dlog(l4)
      t35 = dlog(l3)
      t38 = dlog(-t1)
      t40 = f3(t10,z)
      t42 = dlog(1-t40)
      t43 = dlog(l1)
      t45 = dlog(-t5)
      t47 = dlog(z)
      t49 = f2(t10,z)
      t50 = dlog(t49)
      tncoff5 = t1/(-1+l8)*t8*t23/512-t1*t6/l2*t7*t23/256+t1*(t33+t34+2*
     #t35+4*t38+t42+t43+2*t45-2*t47-t50)*t8*t23/512
      t12 = l3/2
      t13 = Dsquar311(z,l1,0.D0,t12,l4)
      t15 = s2vr2(z,1.D0,1.D0,t12,l4)
      t18 = s23r2(z,1.D0,1.D0,t12,l4)
      t20 = s1vr2(z,1.D0,1.D0,t12,l4)
      t23 = f1(t12,z,1.D0)
      tncoff6 = -(-1+z)/(-1+l8)/(-2+l3)/l3/l8/t13/t15/t18/t20*t23/512
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(l1)
      t6 = 1-l1*l2
      t7 = dlog(t6)
      t8 = f3(0.D0,z)
      t12 = dlog(1-t8+t8*l1*l2)
      t13 = 1-l2
      t14 = dlog(t13)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t20 = dlog(-t1)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t24 = dlog(l4)
      t26 = dlog(1-l4)
      t27 = 2*t2+t4+t7+t12+t14+2*t15-2*t17+4*t20-t23+t24+t26
      t29 = 1/l2
      t30 = Dsquar311(z,l1,l2,0.D0,l4)
      t31 = 1/t30
      t34 = s2vr2(z,t6,t13,0.D0,l4)
      t36 = s23r2(z,t6,t13,0.D0,l4)
      t39 = s1vr2(z,t6,t13,0.D0,l4)
      t41 = f1(0.D0,z,t6)
      t43 = 1/t34/t36/t39*t41
      t46 = t1*t29
      tncoff7 = -t1*t27*t29*t31*t43/1024+t46/l8*t31*t43/1024+t46/l3*t31*
     #t43/512
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(l1)
      t6 = 1-l1*l2
      t7 = dlog(t6)
      t8 = f3(0.D0,z)
      t12 = dlog(1-t8+t8*l1*l2)
      t13 = 1-l2
      t14 = dlog(t13)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t20 = dlog(-t1)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t24 = dlog(l4)
      t26 = dlog(1-l4)
      t27 = 2*t2+t4+t7+t12+t14+2*t15-2*t17+4*t20-t23+t24+t26
      t29 = 1/l2
      t30 = Dsquar311(z,l1,l2,0.D0,l4)
      t31 = 1/t30
      t32 = t29*t31
      t34 = s2vr2(z,t6,t13,0.D0,l4)
      t36 = s23r2(z,t6,t13,0.D0,l4)
      t39 = s1vr2(z,t6,t13,0.D0,l4)
      t41 = f1(0.D0,z,t6)
      t43 = 1/t34/t36/t39*t41
      tncoff8 = -t1*t27*t32*t43/1024-t1/(-1+l8)*t32*t43/1024+t1*t29/l3*t
     #31*t43/512
      t9 = Dsquar311(z,l1,l2,0.D0,l4)
      t12 = 1-l1*l2
      t13 = 1-l2
      t14 = s2vr2(z,t12,t13,0.D0,l4)
      t17 = s23r2(z,t12,t13,0.D0,l4)
      t19 = s1vr2(z,t12,t13,0.D0,l4)
      t22 = f1(0.D0,z,t12)
      tncoff9 = (-1+z)/(-1+l8)/l2/l8/t9/t14/t17/t19*t22/1024
      t9 = l3/2
      t10 = Dsquar311(z,l1,l2,t9,l4)
      t13 = 1-l1*l2
      t14 = 1-l2
      t15 = s2vr2(z,t13,t14,t9,l4)
      t18 = s23r2(z,t13,t14,t9,l4)
      t20 = s1vr2(z,t13,t14,t9,l4)
      t23 = f1(t9,z,t13)
      tncoff10 = (-1+z)/(-2+l3)/l2/l3/t10/t15/t18/t20*t23/256
      t9 = l3/2
      t10 = Dsquar311(z,l1,l2,t9,l4)
      t13 = 1-l1*l2
      t14 = 1-l2
      t15 = s2vr2(z,t13,t14,t9,l4)
      t18 = s23r2(z,t13,t14,t9,l4)
      t20 = s1vr2(z,t13,t14,t9,l4)
      t23 = f1(t9,z,t13)
      tncoff11 = (-1+z)/(-2+l3)/l2/l3/t10/t15/t18/t20*t23/256
      t1 = -1+z
      t4 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = 1/t7
      t9 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = 1/t9
      t12 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = 1/t12
      t14 = f1(0.D0,z,1.D0)
      t16 = t8*t10*t13*t14
      t19 = 1-l8
      t33 = dlog(-t1)
      t35 = f2(0.D0,z)
      t36 = dlog(t35)
      t38 = dlog(2.D0)
      t40 = dlog(z)
      t42 = f3(0.D0,z)
      t44 = dlog(1-t42)
      t46 = dlog(l4)
      t49 = dlog(1-l4)
      t51 = dlog(l1)
      tncoff12 = -t1/l2*t5*t16/1024-t1/l8*t5*t8*t10*t13*t14/2048-t1*(-24
     #*t33+6*t36-12*t38+12*t40-6*t44-6*t46-6*t49-6*t51)*t5*t16/12288-t1/
     #l3*t5*t16/1024
      t1 = -1+z
      t4 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t12 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t14 = f1(0.D0,z,1.D0)
      t16 = 1/t7/t9/t12*t14
      t31 = dlog(-t1)
      t33 = f2(0.D0,z)
      t34 = dlog(t33)
      t36 = dlog(2.D0)
      t38 = dlog(z)
      t40 = f3(0.D0,z)
      t42 = dlog(1-t40)
      t44 = dlog(l4)
      t47 = dlog(1-l4)
      t49 = dlog(l1)
      tncoff13 = -t1/l2*t5*t16/1024+t1/(-1+l8)*t5*t16/2048-t1/l3*t5*t16/
     #1024-t1*(-24*t31+6*t34-12*t36+12*t38-6*t42-6*t44-6*t47-6*t49)*t5*t
     #16/12288
      t6 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t12 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t15 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t17 = f1(0.D0,z,1.D0)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t6 = l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,1.D0,1.D0,t6,l4)
      t13 = s23r2(z,1.D0,1.D0,t6,l4)
      t16 = s1vr2(z,1.D0,1.D0,t6,l4)
      t18 = f1(t6,z,1.D0)
      tncoff15 = -(-1+z)/(-2+l3)/l3/t7/t11/t13/t16*t18/512
      t6 = l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,1.D0,1.D0,t6,l4)
      t13 = s23r2(z,1.D0,1.D0,t6,l4)
      t16 = s1vr2(z,1.D0,1.D0,t6,l4)
      t18 = f1(t6,z,1.D0)
      tncoff16 = -(-1+z)/(-2+l3)/l3/t7/t11/t13/t16*t18/512
      t4 = Dsquar311(z,l1,l2,0.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,0.D0,l4)
      t12 = s23r2(z,t8,t9,0.D0,l4)
      t15 = s1vr2(z,t8,t9,0.D0,l4)
      t17 = f1(0.D0,z,t8)
      tncoff17 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t4 = Dsquar311(z,l1,l2,0.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,0.D0,l4)
      t12 = s23r2(z,t8,t9,0.D0,l4)
      t15 = s1vr2(z,t8,t9,0.D0,l4)
      t17 = f1(0.D0,z,t8)
      tncoff18 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = CDnum1e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = CDnum1e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = CDnum1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = CDnum1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      CD1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function CD1e1int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
      double precision CDnum1e1
      external CDnum1e1
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

      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(-t1)
      t7 = dlog(l4)
      t8 = dlog(l1)
      t9 = dlog(2.D0)
      t11 = f3(0.D0,z)
      t13 = dlog(1-t11)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t16 = dlog(z)
      t20 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t30 = f1(0.D0,z,1.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff1 = t1*(t3+4*t5+t7+t8+2*t9+t13-t15-2*t16)*t21*t32/2048-t1/l8
     #*t21*t32/2048-t1/l3*t21*t32/1024-t1/l2*t21*t32/1024
      t1 = -1+z
      t3 = dlog(1-l4)
      t5 = dlog(-t1)
      t7 = dlog(l4)
      t8 = dlog(l1)
      t9 = dlog(2.D0)
      t11 = f3(0.D0,z)
      t13 = dlog(1-t11)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t16 = dlog(z)
      t20 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t30 = f1(0.D0,z,1.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff2 = t1*(t3+4*t5+t7+t8+2*t9+t13-t15-2*t16)*t21*t32/2048+t1/(-
     #1+l8)*t21*t32/2048-t1/l3*t21*t32/1024-t1/l2*t21*t32/1024
      t6 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t12 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t15 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t17 = f1(0.D0,z,1.D0)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t6 = l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,1.D0,1.D0,t6,l4)
      t13 = s23r2(z,1.D0,1.D0,t6,l4)
      t16 = s1vr2(z,1.D0,1.D0,t6,l4)
      t18 = f1(t6,z,1.D0)
      tncoff4 = -(-1+z)/(-2+l3)/l3/t7/t11/t13/t16*t18/512
      t6 = l3/2
      t7 = Dsquar311(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,1.D0,1.D0,t6,l4)
      t13 = s23r2(z,1.D0,1.D0,t6,l4)
      t16 = s1vr2(z,1.D0,1.D0,t6,l4)
      t18 = f1(t6,z,1.D0)
      tncoff5 = -(-1+z)/(-2+l3)/l3/t7/t11/t13/t16*t18/512
      t4 = Dsquar311(z,l1,l2,0.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,0.D0,l4)
      t12 = s23r2(z,t8,t9,0.D0,l4)
      t15 = s1vr2(z,t8,t9,0.D0,l4)
      t17 = f1(0.D0,z,t8)
      tncoff6 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t4 = Dsquar311(z,l1,l2,0.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,0.D0,l4)
      t12 = s23r2(z,t8,t9,0.D0,l4)
      t15 = s1vr2(z,t8,t9,0.D0,l4)
      t17 = f1(0.D0,z,t8)
      tncoff7 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      CD1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function CD1e2int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
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

      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,0.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t13 = f1(0.D0,z,1.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      CD1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

