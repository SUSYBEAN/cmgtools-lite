CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm3e0int(var,ndim,Tcoeff)

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
      double precision qiqirmnum2e0
      external qiqirmnum2e0
      double precision qiqirmnum2e1
      external qiqirmnum2e1
      double precision qiqirmnum2e2
      external qiqirmnum2e2
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
      t2 = dlog(l3)
      t4 = 1-z
      t5 = dlog(t4)
      t6 = 4*t5
      t7 = dlog(2.D0)
      t8 = 2*t7
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(z)
      t12 = 2*t11
      t13 = dlog(l1)
      t15 = dlog(1-l4)
      t17 = -t4
      t19 = 1/l3
      t20 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = t19*t21
      t24 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t27 = 1/t26
      t29 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t30 = 1/t29
      t31 = f1(0.D0,z,0.D0)
      t33 = t25*t27*t30*t31
      t36 = 1-l8
      t37 = t36*t17
      t38 = 1/l8
      t40 = -1/t36
      t41 = t38*t40
      t47 = t21*t25*t27*t30*t31
      t50 = dlog(l2)
      t54 = 1/l2
      t71 = 0.3141592653589793D1**2
      t72 = t7**2
      t74 = t5**2
      t76 = t10**2
      t78 = t15**2
      t80 = t13**2
      t82 = t1**2
      t86 = -48*t5*t11+12*t1*t7-6*t1*t10+48*t5*t7-12*t1*t11+6*t1*t13-t71
     #+12*t72+48*t74+3*t76+3*t78+3*t80+3*t82+12*t10*t11
      t99 = t11**2
      t117 = 12*t15*t7-6*t15*t10-12*t15*t11-12*t7*t10+6*t15*t13-24*t7*t1
     #1+12*t99+24*t15*t5+6*t15*t1-12*t13*t11+24*t1*t5-6*t13*t10+12*t13*t
     #7+24*t13*t5-24*t5*t10
      t132 = dlog(l8)
      t144 = -t13-t15+t12-t8+t10-t132-t6-t1-2*t11*l8+4*t5*l8+t13*l8+t1*l
     #8+t15*l8-t10*l8+t132*l8+2*t7*l8
      tncoff1 = (t1+2*t2+t6+t8-t10-t12+t13+t15)*t17*t22*t33/2048+t37*t41
     #*t19*t47/2048-t17*(-2*t50-t13-t8-t6+t10+t12-t15-t1)*t54*t21*t33/20
     #48-t17*(t86+t117)*t21*t33/24576+t37*t54*t40*t38*t47/2048-t17*t54*t
     #22*t33/1024+t17*t144*t41*t47/4096
      t1 = dlog(l4)
      t2 = dlog(l3)
      t4 = 1-z
      t5 = dlog(t4)
      t6 = 4*t5
      t7 = dlog(2.D0)
      t8 = 2*t7
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(z)
      t12 = 2*t11
      t13 = dlog(l1)
      t15 = dlog(1-l4)
      t17 = -t4
      t20 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t21 = 1/t20
      t22 = 1/l3*t21
      t24 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t25 = 1/t24
      t26 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t27 = 1/t26
      t29 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t30 = 1/t29
      t31 = f1(0.D0,z,0.D0)
      t33 = t25*t27*t30*t31
      t36 = -1+l8
      t37 = 1/t36
      t42 = dlog(l2)
      t46 = 1/l2
      t59 = 0.3141592653589793D1**2
      t62 = t7**2
      t64 = t5**2
      t66 = t10**2
      t68 = t15**2
      t74 = t11**2
      t76 = t13**2
      t78 = 12*t10*t11+12*t15*t7-6*t15*t10+24*t15*t5-t59-24*t7*t11+12*t6
     #2+48*t64+3*t66+3*t68+6*t15*t13-12*t7*t10+12*t74+3*t76
      t79 = t1**2
      t109 = 3*t79-12*t15*t11+6*t15*t1-12*t13*t11+24*t1*t5-6*t13*t10+12*
     #t13*t7+24*t13*t5-24*t5*t10-48*t5*t11+12*t1*t7-6*t1*t10+48*t5*t7-12
     #*t1*t11+6*t1*t13
      t115 = t17*t46
      t128 = dlog(-t36)
      tncoff2 = (t1+2*t2+t6+t8-t10-t12+t13+t15)*t17*t22*t33/2048+t17*t37
     #*t22*t33/2048-t17*(-2*t42-t13-t8-t6+t10+t12-t15-t1)*t46*t21*t33/20
     #48-t17*(t78+t109)*t21*t33/24576+t115*t37*t21*t33/2048-t115*t22*t33
     #/1024+t17*(t10*l8+2*t11*l8-t1*l8-t128*l8-t15*l8-2*t7*l8-4*t5*l8-t1
     #3*l8)/l8*t37*t21*t25*t27*t30*t31/4096
      t1 = -1+z
      t2 = 1/l8
      t4 = -1+l8
      t5 = 1/t4
      t9 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t11 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t14 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t16 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t19 = f1(0.D0,z,0.D0)
      t21 = 1/t9/t11/t14/t16*t19
      t26 = t5*t2
      t31 = dlog(-t4)
      t32 = dlog(l8)
      t34 = dlog(1-l4)
      t36 = dlog(-t1)
      t38 = f2(0.D0,z)
      t39 = dlog(t38)
      t40 = dlog(z)
      t42 = dlog(l4)
      t43 = dlog(l1)
      t44 = dlog(2.D0)
      tncoff3 = -t1*t2*t5/l3*t21/2048-t1/l2*t26*t21/2048+t1*(t31+t32+t34
     #+4*t36-t39-2*t40+t42+t43+2*t44)*t26*t21/4096
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(l3)
      t5 = dlog(l1)
      t6 = 2-l3
      t7 = dlog(t6)
      t9 = dlog(z)
      t12 = dlog(-t1)
      t14 = l3/2
      t15 = f2(t14,z)
      t16 = dlog(t15)
      t18 = dlog(1-l4)
      t22 = -1/t6
      t23 = 1/l3
      t24 = t22*t23
      t26 = Dsquar321(z,l1,0.D0,t14,l4)
      t27 = 1/t26
      t28 = s2vr2(z,0.D0,0.D0,t14,l4)
      t29 = 1/t28
      t31 = s23r2(z,0.D0,0.D0,t14,l4)
      t32 = 1/t31
      t33 = s1vr2(z,0.D0,0.D0,t14,l4)
      t34 = 1/t33
      t36 = f1(t14,z,0.D0)
      t37 = t32*t34*t36
      t41 = 1-l8
      tncoff4 = t1*(t2+2*t3+t5+2*t7-2*t9+4*t12-t16+t18)*t24*t27*t29*t37/
     #1024-t1*t24/l8*t27*t29*t37/1024-t1*t22*t29*t34*t32*t27/l2*t23*t36/
     #512
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(l3)
      t5 = dlog(l1)
      t6 = 2-l3
      t7 = dlog(t6)
      t9 = dlog(z)
      t12 = dlog(-t1)
      t14 = l3/2
      t15 = f2(t14,z)
      t16 = dlog(t15)
      t18 = dlog(1-l4)
      t22 = -1/t6
      t23 = 1/l3
      t26 = Dsquar321(z,l1,0.D0,t14,l4)
      t27 = 1/t26
      t28 = s2vr2(z,0.D0,0.D0,t14,l4)
      t29 = 1/t28
      t31 = s23r2(z,0.D0,0.D0,t14,l4)
      t32 = 1/t31
      t33 = s1vr2(z,0.D0,0.D0,t14,l4)
      t34 = 1/t33
      t36 = f1(t14,z,0.D0)
      t38 = t27*t29*t32*t34*t36
      t41 = t1*t22
      tncoff5 = t1*(t2+2*t3+t5+2*t7-2*t9+4*t12-t16+t18)*t22*t23*t38/1024
     #+t41*t23/(-1+l8)*t38/1024-t41*t29*t34*t32*t27/l2*t23*t36/512
      t12 = l3/2
      t13 = Dsquar321(z,l1,0.D0,t12,l4)
      t15 = s2vr2(z,0.D0,0.D0,t12,l4)
      t18 = s23r2(z,0.D0,0.D0,t12,l4)
      t20 = s1vr2(z,0.D0,0.D0,t12,l4)
      t23 = f1(t12,z,0.D0)
      tncoff6 = -(-1+z)/(l3-2)/l3/(-1+l8)/l8/t13/t15/t18/t20*t23/1024
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(l1)
      t5 = l1*l2
      t6 = f3(0.D0,z)
      t9 = dlog(1-t5*t6)
      t11 = dlog(-t1)
      t14 = dlog(1-l2)
      t15 = dlog(2.D0)
      t18 = dlog(1-t5)
      t19 = dlog(z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = -2*t2-t4-t9-4*t11-t14-2*t15-t18+2*t19+t22-t24-t25
      t28 = 1/l2
      t29 = Dsquar321(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t33 = s2vr2(z,t5,l2,0.D0,l4)
      t35 = s23r2(z,t5,l2,0.D0,l4)
      t38 = s1vr2(z,t5,l2,0.D0,l4)
      t40 = f1(0.D0,z,t5)
      t42 = 1/t33/t35/t38*t40
      t45 = t1*t28
      tncoff7 = t1*t26*t28*t30*t42/2048+t45/l8*t30*t42/2048+t45/l3*t30*t
     #42/1024
      t1 = -1+z
      t2 = dlog(l2)
      t4 = dlog(l1)
      t5 = l1*l2
      t6 = f3(0.D0,z)
      t9 = dlog(1-t5*t6)
      t11 = dlog(-t1)
      t14 = dlog(1-l2)
      t15 = dlog(2.D0)
      t18 = dlog(1-t5)
      t19 = dlog(z)
      t21 = f2(0.D0,z)
      t22 = dlog(t21)
      t24 = dlog(1-l4)
      t25 = dlog(l4)
      t26 = -2*t2-t4-t9-4*t11-t14-2*t15-t18+2*t19+t22-t24-t25
      t28 = 1/l2
      t29 = Dsquar321(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t33 = s2vr2(z,t5,l2,0.D0,l4)
      t35 = s23r2(z,t5,l2,0.D0,l4)
      t38 = s1vr2(z,t5,l2,0.D0,l4)
      t40 = f1(0.D0,z,t5)
      t42 = 1/t33/t35/t38*t40
      t45 = t1*t28
      tncoff8 = t1*t26*t28*t30*t42/2048-t45/(-1+l8)*t30*t42/2048+t45/l3*
     #t30*t42/1024
      t9 = Dsquar321(z,l1,l2,0.D0,l4)
      t11 = l1*l2
      t12 = s2vr2(z,t11,l2,0.D0,l4)
      t15 = s23r2(z,t11,l2,0.D0,l4)
      t17 = s1vr2(z,t11,l2,0.D0,l4)
      t20 = f1(0.D0,z,t11)
      tncoff9 = (-1+z)/l2/(-1+l8)/l8/t9/t12/t15/t17*t20/2048
      t2 = l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t6 = s2vr2(z,t3,l2,t2,l4)
      t8 = s1vr2(z,t3,l2,t2,l4)
      t12 = s23r2(z,t3,l2,t2,l4)
      t17 = Dsquar321(z,l1,l2,t2,l4)
      tncoff10 = (-1+z)*t4/t6/t8/t12/(l3-2)/t17/l2/l3/512
      t2 = l3/2
      t3 = l1*l2
      t4 = f1(t2,z,t3)
      t6 = s2vr2(z,t3,l2,t2,l4)
      t8 = s1vr2(z,t3,l2,t2,l4)
      t12 = s23r2(z,t3,l2,t2,l4)
      t17 = Dsquar321(z,l1,l2,t2,l4)
      tncoff11 = (-1+z)*t4/t6/t8/t12/(l3-2)/t17/l2/l3/512
      t1 = -1+z
      t4 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = 1/t7
      t9 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = 1/t9
      t12 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = 1/t12
      t14 = f1(0.D0,z,0.D0)
      t16 = t8*t10*t13*t14
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t27 = dlog(2.D0)
      t29 = dlog(l1)
      t32 = dlog(1-l4)
      t35 = dlog(-t1)
      t37 = dlog(z)
      t39 = dlog(l4)
      t46 = 1-l8
      tncoff12 = -t1/l3*t5*t16/2048-t1/l2*t5*t16/2048-t1*(6*t25-12*t27-6
     #*t29-6*t32-24*t35+12*t37-6*t39)*t5*t16/24576-t1/l8*t5*t8*t10*t13*t
     #14/4096
      t1 = -1+z
      t4 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t9 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t12 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t14 = f1(0.D0,z,0.D0)
      t16 = 1/t7/t9/t12*t14
      t24 = f2(0.D0,z)
      t25 = dlog(t24)
      t27 = dlog(2.D0)
      t29 = dlog(l1)
      t32 = dlog(1-l4)
      t35 = dlog(-t1)
      t37 = dlog(z)
      t39 = dlog(l4)
      tncoff13 = -t1/l3*t5*t16/2048-t1/l2*t5*t16/2048-t1*(6*t25-12*t27-6
     #*t29-6*t32-24*t35+12*t37-6*t39)*t5*t16/24576+t1/(-1+l8)*t5*t16/409
     #6
      t6 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t17 = f1(0.D0,z,0.D0)
      tncoff14 = -(-1+z)/l8/(-1+l8)/t6/t10/t12/t15*t17/4096
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff15 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/1024
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff16 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/1024
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff17 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff18 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum2e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum2e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum2e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum2e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm3e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm3e1int(var,ndim,Tcoeff)

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
      double precision qiqirmnum2e0
      external qiqirmnum2e0
      double precision qiqirmnum2e1
      external qiqirmnum2e1
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
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = dlog(l1)
      t7 = dlog(z)
      t10 = dlog(-t1)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
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
      tncoff1 = -t1*(-2*t2+t5-t6+2*t7-4*t10-t13-t14)*t18*t29/4096-t1/l8*
     #t18*t21*t23*t26*t27/4096-t1/l3*t18*t29/2048-t1/l2*t18*t29/2048
      t1 = -1+z
      t2 = dlog(2.D0)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t6 = dlog(l1)
      t7 = dlog(z)
      t10 = dlog(-t1)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t17 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t27 = f1(0.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff2 = -t1*(-2*t2+t5-t6+2*t7-4*t10-t13-t14)*t18*t29/4096+t1/(-1
     #+l8)*t18*t29/4096-t1/l3*t18*t29/2048-t1/l2*t18*t29/2048
      t6 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t17 = f1(0.D0,z,0.D0)
      tncoff3 = -(-1+z)/l8/(-1+l8)/t6/t10/t12/t15*t17/4096
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff4 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/1024
      t6 = l3/2
      t7 = Dsquar321(z,l1,0.D0,t6,l4)
      t11 = s2vr2(z,0.D0,0.D0,t6,l4)
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t16 = s1vr2(z,0.D0,0.D0,t6,l4)
      t18 = f1(t6,z,0.D0)
      tncoff5 = -(-1+z)/(l3-2)/l3/t7/t11/t13/t16*t18/1024
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff6 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t4 = Dsquar321(z,l1,l2,0.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,0.D0,l4)
      t10 = s23r2(z,t7,l2,0.D0,l4)
      t13 = s1vr2(z,t7,l2,0.D0,l4)
      t15 = f1(0.D0,z,t7)
      tncoff7 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm3e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm3e2int(var,ndim,Tcoeff)

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
      double precision qiqirmnum2e0
      external qiqirmnum2e0
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
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t13 = f1(0.D0,z,0.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm3e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

