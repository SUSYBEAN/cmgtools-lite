CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm1e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
      double precision qiqirmnum1e2
      external qiqirmnum1e2
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
      t2 = dlog(l1)
      t4 = dlog(z)
      t5 = 2*t4
      t6 = dlog(l4)
      t7 = f2(0.D0,z)
      t8 = dlog(t7)
      t10 = dlog(1-l4)
      t11 = dlog(l8)
      t12 = dlog(2.D0)
      t13 = 2*t12
      t15 = dlog(-t1)
      t16 = 4*t15
      t27 = t2*l8-t2+t5-t6+t8-t10-t11-t13-t16+t6*l8+t10*l8+t11*l8-2*t4*l
     #8+4*t15*l8-t8*l8+2*t12*l8
      t29 = -1+l8
      t30 = 1/t29
      t31 = 1/l8
      t32 = t30*t31
      t34 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t35 = 1/t34
      t36 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t37 = 1/t36
      t39 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t40 = 1/t39
      t41 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t42 = 1/t41
      t44 = f1(0.D0,z,0.D0)
      t46 = t35*t37*t40*t42*t44
      t50 = -t1*t29
      t51 = 1/l2
      t58 = 1/l3
      t59 = t58*t35
      t63 = t37*t40*t42*t44
      t66 = dlog(l3)
      t101 = -6*t2*t8-12*t8*t12+12*t6*t12-48*t15*t4-24*t4*t12+6*t10*t6+1
     #2*t8*t4-12*t10*t4+6*t10*t2-6*t6*t8-6*t10*t8+6*t6*t2+48*t15*t12+12*
     #t10*t12
      t114 = t2**2
      t116 = t10**2
      t118 = t6**2
      t122 = 0.3141592653589793D1**2
      t123 = t12**2
      t125 = t15**2
      t127 = t8**2
      t129 = t4**2
      t131 = 24*t10*t15-24*t15*t8-12*t6*t4+24*t6*t15+12*t2*t12+24*t2*t15
     #+3*t114+3*t116+3*t118-12*t2*t4-t122+12*t123+48*t125+3*t127+12*t129
      t141 = dlog(l2)
      tncoff1 = t1*t27*t32*t46/2048+t50*t31*t51*t30*t46/1024-t1*t51*t59*
     #t63/512-(-t10-t6-t13-t16+t8+t5-t2-2*t66)*t1*t59*t63/1024-t1*(t101+
     #t131)*t35*t63/12288+t50*t32*t58*t46/1024-t1*(-t6-t10-t2-t13-t16+t8
     #+t5-2*t141)*t51*t35*t63/1024
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t6 = dlog(-t1)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(z)
      t16 = dlog(1-l4)
      t18 = dlog(l4)
      t20 = dlog(l1)
      t22 = dlog(2.D0)
      t28 = -1/t2
      t32 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t33 = 1/t32
      t34 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t35 = 1/t34
      t37 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t38 = 1/t37
      t39 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t40 = 1/t39
      t42 = f1(0.D0,z,0.D0)
      t47 = 1/l2
      t48 = t1*t47
      t53 = t35*t38*t40*t42
      t57 = 1/l3*t33
      t61 = 2*t22
      t62 = 4*t6
      t63 = 2*t12
      t64 = dlog(l3)
      t81 = 0.3141592653589793D1**2
      t84 = t22**2
      t86 = t6**2
      t88 = t10**2
      t90 = t12**2
      t92 = t20**2
      t94 = t16**2
      t96 = t18**2
      t102 = 24*t16*t6+48*t6*t22+12*t16*t22-t81+6*t18*t20+12*t84+48*t86+
     #3*t88+12*t90+3*t92+3*t94+3*t96-6*t18*t10-6*t16*t10
      t133 = 6*t16*t20+12*t10*t12-12*t16*t12+6*t16*t18-48*t6*t12-24*t12*
     #t22+12*t18*t22-12*t10*t22-6*t20*t10+24*t20*t6-12*t20*t12+12*t20*t2
     #2+24*t18*t6-12*t18*t12-24*t10*t6
      t139 = dlog(l2)
      tncoff2 = t1*(-t3*l8-4*t6*l8+t10*l8+2*t12*l8-t16*l8-t18*l8-t20*l8-
     #2*t22*l8)*t28/l8*t33*t35*t38*t40*t42/2048+t48*t28*t33*t53/1024-t48
     #*t57*t53/512-(-t16-t18-t61-t62+t10+t63-t20-2*t64)*t1*t57*t53/1024+
     #t1*t28*t57*t53/1024-t1*(t102+t133)*t33*t53/12288-t1*(-t18-t16-t20-
     #t61-t62+t10+t63-2*t139)*t47*t33*t53/1024
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = dlog(-t1)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t10 = dlog(z)
      t12 = dlog(l8)
      t14 = dlog(1-l4)
      t15 = dlog(l4)
      t16 = 1-l8
      t17 = dlog(t16)
      t21 = -1/t16
      t22 = 1/l8
      t25 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t27 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t30 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t32 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t35 = f1(0.D0,z,0.D0)
      t37 = 1/t25/t27/t30/t32*t35
      t40 = t1*t22
      tncoff3 = t1*(t2+2*t3+4*t6-t9-2*t10+t12+t14+t15+t17)*t21*t22*t37/2
     #048-t40/l2*t21*t37/1024-t40*t21/l3*t37/1024
      t1 = -1+z
      t2 = l3-2
      t3 = 1/t2
      t4 = t1*t3
      t5 = l3/2
      t6 = s1vr2(z,0.D0,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t10 = 1/t6/t8
      t12 = s2vr2(z,0.D0,0.D0,t5,l4)
      t14 = Dsquar321(z,l1,0.D0,t5,l4)
      t16 = 1/t12/t14
      t18 = 1/l3
      t20 = f1(t5,z,0.D0)
      t26 = dlog(1-l4)
      t27 = dlog(l4)
      t28 = dlog(l1)
      t30 = dlog(-t2)
      t32 = f2(t5,z)
      t33 = dlog(t32)
      t34 = dlog(z)
      t37 = dlog(-t1)
      t39 = dlog(l3)
      t46 = t16*t10*t20
      tncoff4 = -t4*t10*t16/l2*t18*t20/256+t1*(t26+t27+t28+2*t30-t33-2*t
     #34+4*t37+2*t39)*t3*t18*t46/512-t4*t18/l8*t46/512
      t1 = -1+z
      t2 = l3-2
      t3 = 1/t2
      t4 = t1*t3
      t5 = l3/2
      t6 = s1vr2(z,0.D0,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t10 = 1/t6/t8
      t12 = s2vr2(z,0.D0,0.D0,t5,l4)
      t14 = Dsquar321(z,l1,0.D0,t5,l4)
      t16 = 1/t12/t14
      t18 = 1/l3
      t20 = f1(t5,z,0.D0)
      t26 = dlog(1-l4)
      t27 = dlog(l4)
      t28 = dlog(l1)
      t30 = dlog(-t2)
      t32 = f2(t5,z)
      t33 = dlog(t32)
      t34 = dlog(z)
      t37 = dlog(-t1)
      t39 = dlog(l3)
      t46 = t16*t10*t20
      tncoff5 = -t4*t10*t16/l2*t18*t20/256+t1*(t26+t27+t28+2*t30-t33-2*t
     #34+4*t37+2*t39)*t3*t18*t46/512+t4*t18/(-1+l8)*t46/512
      t12 = l3/2
      t13 = Dsquar321(z,l1,0.D0,t12,l4)
      t15 = s2vr2(z,0.D0,0.D0,t12,l4)
      t18 = s23r2(z,0.D0,0.D0,t12,l4)
      t20 = s1vr2(z,0.D0,0.D0,t12,l4)
      t23 = f1(t12,z,0.D0)
      tncoff6 = -(-1+z)/(l3-2)/l3/(-1+l8)/l8/t13/t15/t18/t20*t23/512
      t1 = -1+z
      t2 = 1-l8
      t5 = 1/l2
      t11 = Dsquar321(z,l1,l2,0.D0,l4)
      t12 = 1/t11
      t13 = l1*l2
      t14 = s2vr2(z,t13,l2,0.D0,l4)
      t15 = 1/t14
      t17 = s23r2(z,t13,l2,0.D0,l4)
      t18 = 1/t17
      t19 = s1vr2(z,t13,l2,0.D0,l4)
      t20 = 1/t19
      t22 = f1(0.D0,z,t13)
      t33 = t15*t18*t20*t22
      t36 = dlog(l4)
      t38 = dlog(1-l4)
      t39 = dlog(l1)
      t41 = dlog(-t1)
      t43 = f3(0.D0,z)
      t46 = dlog(1-t13*t43)
      t47 = dlog(2.D0)
      t49 = f2(0.D0,z)
      t50 = dlog(t49)
      t52 = dlog(1-l2)
      t53 = dlog(z)
      t56 = dlog(1-t13)
      t57 = dlog(l2)
      t59 = t36+t38+t39+4*t41+t46+2*t47-t50+t52-2*t53+t56+2*t57
      tncoff7 = t1/l8*t5*t12*t15*t18*t20*t22/1024+t1*t5/l3*t12*t33/512-t
     #1*t59*t5*t12*t33/1024
      t1 = -1+z
      t2 = 1/l2
      t3 = t1*t2
      t6 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = 1/t6
      t10 = l1*l2
      t11 = s2vr2(z,t10,l2,0.D0,l4)
      t13 = s23r2(z,t10,l2,0.D0,l4)
      t16 = s1vr2(z,t10,l2,0.D0,l4)
      t18 = f1(0.D0,z,t10)
      t20 = 1/t11/t13/t16*t18
      t28 = dlog(l4)
      t30 = dlog(1-l4)
      t31 = dlog(l1)
      t33 = dlog(-t1)
      t35 = f3(0.D0,z)
      t38 = dlog(1-t10*t35)
      t39 = dlog(2.D0)
      t41 = f2(0.D0,z)
      t42 = dlog(t41)
      t44 = dlog(1-l2)
      t45 = dlog(z)
      t48 = dlog(1-t10)
      t49 = dlog(l2)
      t51 = t28+t30+t31+4*t33+t38+2*t39-t42+t44-2*t45+t48+2*t49
      tncoff8 = -t3/(-1+l8)*t7*t20/1024+t3/l3*t7*t20/512-t1*t51*t2*t7*t2
     #0/1024
      t9 = Dsquar321(z,l1,l2,0.D0,l4)
      t11 = l1*l2
      t12 = s2vr2(z,t11,l2,0.D0,l4)
      t15 = s23r2(z,t11,l2,0.D0,l4)
      t17 = s1vr2(z,t11,l2,0.D0,l4)
      t20 = f1(0.D0,z,t11)
      tncoff9 = (-1+z)/l8/l2/(-1+l8)/t9/t12/t15/t17*t20/1024
      t2 = l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t6 = s1vr2(z,t3,l2,t2,l4)
      t8 = s23r2(z,t3,l2,t2,l4)
      t12 = s2vr2(z,t3,l2,t2,l4)
      t17 = Dsquar321(z,l1,l2,t2,l4)
      tncoff10 = (-1+z)*t4/t6/t8/t12/(l3-2)/t17/l2/l3/256
      t2 = l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t6 = s1vr2(z,t3,l2,t2,l4)
      t8 = s23r2(z,t3,l2,t2,l4)
      t12 = s2vr2(z,t3,l2,t2,l4)
      t17 = Dsquar321(z,l1,l2,t2,l4)
      tncoff11 = (-1+z)*t4/t6/t8/t12/(l3-2)/t17/l2/l3/256
      t1 = -1+z
      t2 = 1-l8
      t9 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t10 = 1/t9
      t11 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t12 = 1/t11
      t14 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t15 = 1/t14
      t16 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t17 = 1/t16
      t19 = f1(0.D0,z,0.D0)
      t29 = t12*t15*t17*t19
      t33 = dlog(-t1)
      t35 = dlog(2.D0)
      t37 = f2(0.D0,z)
      t38 = dlog(t37)
      t40 = dlog(l1)
      t42 = dlog(l4)
      t45 = dlog(1-l4)
      t47 = dlog(z)
      tncoff12 = -t1/l8*t10*t12*t15*t17*t19/2048-t1/l3*t10*t29/1024-t1*(
     #-24*t33-12*t35+6*t38-6*t40-6*t42-6*t45+12*t47)*t10*t29/12288-t1/l2
     #*t10*t29/1024
      t1 = -1+z
      t5 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t10 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t13 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t15 = f1(0.D0,z,0.D0)
      t17 = 1/t8/t10/t13*t15
      t26 = dlog(-t1)
      t28 = dlog(2.D0)
      t30 = f2(0.D0,z)
      t31 = dlog(t30)
      t33 = dlog(l1)
      t35 = dlog(l4)
      t38 = dlog(1-l4)
      t40 = dlog(z)
      tncoff13 = t1/(-1+l8)*t6*t17/2048-t1/l3*t6*t17/1024-t1*(-24*t26-12
     #*t28+6*t31-6*t33-6*t35-6*t38+12*t40)*t6*t17/12288-t1/l2*t6*t17/102
     #4
      t6 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t17 = f1(0.D0,z,0.D0)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff15 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/512
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff16 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/512
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff17 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff18 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum1e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum1e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum1e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum1e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm1e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
      double precision qiqirmnum1e1
      external qiqirmnum1e1
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
      t2 = dlog(2.D0)
      t5 = dlog(-t1)
      t7 = dlog(l1)
      t8 = dlog(z)
      t10 = dlog(l4)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(1-l4)
      t17 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t23 = 1/t22
      t25 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t26 = 1/t25
      t27 = f1(0.D0,z,0.D0)
      t29 = t21*t23*t26*t27
      t32 = 1-l8
      tncoff1 = -t1*(-2*t2-4*t5-t7+2*t8-t10+t12-t14)*t18*t29/2048-t1/l8*
     #t18*t21*t23*t26*t27/2048-t1/l3*t18*t29/1024-t1/l2*t18*t29/1024
      t1 = -1+z
      t2 = dlog(2.D0)
      t5 = dlog(-t1)
      t7 = dlog(l1)
      t8 = dlog(z)
      t10 = dlog(l4)
      t11 = f2(0.D0,z)
      t12 = dlog(t11)
      t14 = dlog(1-l4)
      t17 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t27 = f1(0.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff2 = -t1*(-2*t2-4*t5-t7+2*t8-t10+t12-t14)*t18*t29/2048+t1/(-1
     #+l8)*t18*t29/2048-t1/l3*t18*t29/1024-t1/l2*t18*t29/1024
      t6 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t17 = f1(0.D0,z,0.D0)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff4 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/512
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff5 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/512
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff6 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff7 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm1e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum1e0
      external qiqirmnum1e0
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

      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

