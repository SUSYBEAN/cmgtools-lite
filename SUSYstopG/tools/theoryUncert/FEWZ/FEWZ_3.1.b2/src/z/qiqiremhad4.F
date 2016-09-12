CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm4e0int(var,ndim,Tcoeff)

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

      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t4 = 2*t3
      t6 = dlog(-t1)
      t7 = 4*t6
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t10 = dlog(z)
      t11 = 2*t10
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t16 = dlog(l8)
      t27 = -t2-t4-t7+t9+t11-t13+t14*l8+t16*l8+t2*l8+4*t6*l8-t9*l8+2*t3*
     #l8-2*t10*l8-t14-t16+t13*l8
      t31 = 1/l8
      t34 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t35 = 1/t34
      t36 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t37 = 1/t36
      t39 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t40 = 1/t39
      t41 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t42 = 1/t41
      t44 = f1(1.D0,z,0.D0)
      t49 = dlog(l3)
      t54 = 1/l3*t35
      t58 = t37*t40*t42*t44
      t65 = dlog(l2)
      t69 = 1/l2
      t78 = t14**2
      t80 = t10**2
      t84 = t9**2
      t86 = t13**2
      t88 = t2**2
      t102 = 12*t2*t3-6*t13*t9+3*t78+12*t80+12*t9*t10+3*t84+3*t86+3*t88-
     #48*t6*t10+12*t14*t3+48*t6*t3-12*t14*t10+6*t14*t2-6*t2*t9
      t103 = 0.3141592653589793D1**2
      t104 = t3**2
      t130 = t6**2
      t132 = -t103+12*t104+24*t2*t6-6*t14*t9-12*t9*t3+12*t13*t3-12*t13*t
     #10+6*t13*t2-24*t3*t10+24*t13*t6+6*t13*t14-12*t2*t10+24*t14*t6-24*t
     #6*t9+48*t130
      t138 = t1*t69
      tncoff1 = t1*t27/(-1+l8)*t31*t35*t37*t40*t42*t44/4096+(t14+t13+2*t
     #49+t2+t4+t7-t9-t11)*t1*t54*t58/2048-t1*t31*t54*t58/2048+t1*(t2+t7+
     #t4-t9-t11+2*t65+t13+t14)*t69*t35*t58/2048-t1*(t102+t132)*t35*t58/2
     #4576-t138*t31*t35*t58/2048-t138*t54*t58/1024
      t1 = -1+z
      t2 = 1/l2
      t3 = t1*t2
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t7 = 1/l3*t6
      t9 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t10 = 1/t9
      t11 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t12 = 1/t11
      t14 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t15 = 1/t14
      t16 = f1(1.D0,z,0.D0)
      t18 = t10*t12*t15*t16
      t22 = dlog(1-l4)
      t24 = 1-l8
      t25 = dlog(t24)
      t27 = dlog(l1)
      t29 = dlog(l4)
      t32 = dlog(-t1)
      t35 = f2(1.D0,z)
      t36 = dlog(t35)
      t38 = dlog(2.D0)
      t41 = dlog(z)
      t47 = -1/t24
      t57 = dlog(l3)
      t59 = 2*t38
      t60 = 4*t32
      t61 = 2*t41
      t71 = dlog(l2)
      t111 = -12*t36*t38-6*t29*t36-24*t32*t36-6*t22*t36-6*t27*t36+12*t36
     #*t41+6*t22*t27-12*t22*t41+12*t22*t38+6*t29*t27-12*t29*t41+48*t32*t
     #38+12*t29*t38-48*t32*t41
      t126 = 0.3141592653589793D1**2
      t127 = t38**2
      t129 = t32**2
      t131 = t22**2
      t133 = t27**2
      t135 = t29**2
      t137 = t41**2
      t139 = t36**2
      t141 = 24*t27*t32+12*t27*t38+24*t29*t32-12*t27*t41+6*t22*t29+24*t2
     #2*t32-24*t38*t41-t126+12*t127+48*t129+3*t131+3*t133+3*t135+12*t137
     #+3*t139
      tncoff2 = -t3*t7*t18/1024+t1*(-t22*l8-t25*l8-t27*l8-t29*l8-4*t32*l
     #8+t36*l8-2*t38*l8+2*t41*l8)*t47/l8*t6*t10*t12*t15*t16/4096+(t29+t2
     #2+2*t57+t27+t59+t60-t36-t61)*t1*t7*t18/2048+t1*t47*t7*t18/2048+t1*
     #(t27+t60+t59-t36-t61+2*t71+t22+t29)*t2*t6*t18/2048+t3*t47*t6*t18/2
     #048-t1*(t111+t141)*t6*t18/24576
      t1 = -1+z
      t2 = f2(1.D0,z)
      t3 = dlog(t2)
      t4 = dlog(l1)
      t6 = dlog(-t1)
      t8 = dlog(2.D0)
      t10 = dlog(z)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t15 = dlog(l8)
      t16 = 1-l8
      t17 = dlog(t16)
      t21 = -1/t16
      t22 = 1/l8
      t23 = t21*t22
      t25 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t27 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t30 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t32 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t35 = f1(1.D0,z,0.D0)
      t37 = 1/t25/t27/t30/t32*t35
      tncoff3 = t1*(-t3+t4+4*t6+2*t8-2*t10+t13+t14+t15+t17)*t23*t37/4096
     #-t1*t21*t22/l3*t37/2048-t1/l2*t23*t37/2048
      t1 = -1+z
      t2 = l3-2
      t3 = 1/t2
      t4 = t1*t3
      t6 = 1-l3/2
      t7 = s2vr2(z,0.D0,0.D0,t6,l4)
      t8 = 1/t7
      t9 = s1vr2(z,0.D0,0.D0,t6,l4)
      t10 = 1/t9
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t14 = 1/t13
      t15 = Dsquar321(z,l1,0.D0,t6,l4)
      t16 = 1/t15
      t19 = 1/l3
      t21 = f1(t6,z,0.D0)
      t26 = dlog(l4)
      t28 = dlog(1-l4)
      t29 = dlog(l3)
      t31 = dlog(l1)
      t33 = dlog(-t1)
      t35 = dlog(z)
      t38 = dlog(-t2)
      t40 = f2(t6,z)
      t41 = dlog(t40)
      t49 = t16*t8*t14*t10*t21
      tncoff4 = -t4*t8*t10*t14*t16/l2*t19*t21/512+t1*(t26+t28+2*t29+t31+
     #4*t33-2*t35+2*t38-t41)*t3*t19*t49/1024-t4*t19/l8*t49/1024
      t1 = -1+z
      t2 = l3-2
      t3 = 1/t2
      t4 = t1*t3
      t6 = 1-l3/2
      t7 = s2vr2(z,0.D0,0.D0,t6,l4)
      t8 = 1/t7
      t9 = s1vr2(z,0.D0,0.D0,t6,l4)
      t10 = 1/t9
      t13 = s23r2(z,0.D0,0.D0,t6,l4)
      t14 = 1/t13
      t15 = Dsquar321(z,l1,0.D0,t6,l4)
      t16 = 1/t15
      t19 = 1/l3
      t21 = f1(t6,z,0.D0)
      t26 = dlog(l4)
      t28 = dlog(1-l4)
      t29 = dlog(l3)
      t31 = dlog(l1)
      t33 = dlog(-t1)
      t35 = dlog(z)
      t38 = dlog(-t2)
      t40 = f2(t6,z)
      t41 = dlog(t40)
      t49 = t16*t8*t14*t10*t21
      tncoff5 = -t4*t8*t10*t14*t16/l2*t19*t21/512+t1*(t26+t28+2*t29+t31+
     #4*t33-2*t35+2*t38-t41)*t3*t19*t49/1024+t4*t19/(-1+l8)*t49/1024
      t13 = 1-l3/2
      t14 = Dsquar321(z,l1,0.D0,t13,l4)
      t16 = s2vr2(z,0.D0,0.D0,t13,l4)
      t19 = s23r2(z,0.D0,0.D0,t13,l4)
      t21 = s1vr2(z,0.D0,0.D0,t13,l4)
      t24 = f1(t13,z,0.D0)
      tncoff6 = -(-1+z)/(l3-2)/l3/l8/(-1+l8)/t14/t16/t19/t21*t24/1024
      t1 = -1+z
      t2 = 1/l2
      t5 = Dsquar321(z,l1,l2,1.D0,l4)
      t6 = 1/t5
      t9 = l1*l2
      t10 = s2vr2(z,t9,l2,1.D0,l4)
      t11 = 1/t10
      t12 = s23r2(z,t9,l2,1.D0,l4)
      t13 = 1/t12
      t15 = s1vr2(z,t9,l2,1.D0,l4)
      t16 = 1/t15
      t17 = f1(1.D0,z,t9)
      t19 = t11*t13*t16*t17
      t22 = dlog(l1)
      t23 = f3(1.D0,z)
      t26 = dlog(1-t9*t23)
      t27 = dlog(2.D0)
      t29 = f2(1.D0,z)
      t30 = dlog(t29)
      t32 = dlog(-t1)
      t35 = dlog(1-t9)
      t36 = dlog(z)
      t39 = dlog(1-l2)
      t40 = dlog(l2)
      t43 = dlog(1-l4)
      t44 = dlog(l4)
      t45 = t22+t26+2*t27-t30+4*t32+t35-2*t36+t39+2*t40+t43+t44
      t51 = 1-l8
      tncoff7 = t1*t2/l3*t6*t19/1024-t1*t45*t2*t6*t19/2048+t1*t2/l8*t6*t
     #11*t13*t16*t17/2048
      t1 = -1+z
      t2 = 1/l2
      t3 = t1*t2
      t5 = Dsquar321(z,l1,l2,1.D0,l4)
      t6 = 1/t5
      t9 = l1*l2
      t10 = s2vr2(z,t9,l2,1.D0,l4)
      t12 = s23r2(z,t9,l2,1.D0,l4)
      t15 = s1vr2(z,t9,l2,1.D0,l4)
      t17 = f1(1.D0,z,t9)
      t19 = 1/t10/t12/t15*t17
      t22 = dlog(l1)
      t23 = f3(1.D0,z)
      t26 = dlog(1-t9*t23)
      t27 = dlog(2.D0)
      t29 = f2(1.D0,z)
      t30 = dlog(t29)
      t32 = dlog(-t1)
      t35 = dlog(1-t9)
      t36 = dlog(z)
      t39 = dlog(1-l2)
      t40 = dlog(l2)
      t43 = dlog(1-l4)
      t44 = dlog(l4)
      t45 = t22+t26+2*t27-t30+4*t32+t35-2*t36+t39+2*t40+t43+t44
      tncoff8 = t3/l3*t6*t19/1024-t1*t45*t2*t6*t19/2048-t3/(-1+l8)*t6*t1
     #9/2048
      t9 = Dsquar321(z,l1,l2,1.D0,l4)
      t11 = l1*l2
      t12 = s2vr2(z,t11,l2,1.D0,l4)
      t15 = s23r2(z,t11,l2,1.D0,l4)
      t17 = s1vr2(z,t11,l2,1.D0,l4)
      t20 = f1(1.D0,z,t11)
      tncoff9 = (-1+z)/l2/(-1+l8)/l8/t9/t12/t15/t17*t20/2048
      t3 = 1-l3/2
      t4 = l1*l2
      t5 = f1(t3,z,t4)
      t7 = s2vr2(z,t4,l2,t3,l4)
      t13 = s1vr2(z,t4,l2,t3,l4)
      t15 = s23r2(z,t4,l2,t3,l4)
      t18 = Dsquar321(z,l1,l2,t3,l4)
      tncoff10 = (-1+z)*t5/t7/(l3-2)/t13/t15/t18/l2/l3/512
      t3 = 1-l3/2
      t4 = l1*l2
      t5 = f1(t3,z,t4)
      t7 = s2vr2(z,t4,l2,t3,l4)
      t13 = s1vr2(z,t4,l2,t3,l4)
      t15 = s23r2(z,t4,l2,t3,l4)
      t18 = Dsquar321(z,l1,l2,t3,l4)
      tncoff11 = (-1+z)*t5/t7/(l3-2)/t13/t15/t18/l2/l3/512
      t1 = -1+z
      t2 = 1-l8
      t9 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t10 = 1/t9
      t11 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = 1/t11
      t14 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = 1/t14
      t16 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = 1/t16
      t19 = f1(1.D0,z,0.D0)
      t29 = t12*t15*t17*t19
      t37 = dlog(l1)
      t39 = dlog(l4)
      t42 = dlog(1-l4)
      t44 = f2(1.D0,z)
      t45 = dlog(t44)
      t47 = dlog(z)
      t49 = dlog(2.D0)
      t52 = dlog(-t1)
      tncoff12 = -t1/l8*t10*t12*t15*t17*t19/4096-t1/l3*t10*t29/2048-t1/l
     #2*t10*t29/2048-t1*(-6*t37-6*t39-6*t42+6*t45+12*t47-12*t49-24*t52)*
     #t10*t29/24576
      t1 = -1+z
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t10 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t13 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t15 = f1(1.D0,z,0.D0)
      t17 = 1/t8/t10/t13*t15
      t30 = dlog(l1)
      t32 = dlog(l4)
      t35 = dlog(1-l4)
      t37 = f2(1.D0,z)
      t38 = dlog(t37)
      t40 = dlog(z)
      t42 = dlog(2.D0)
      t45 = dlog(-t1)
      tncoff13 = t1/(-1+l8)*t6*t17/4096-t1/l3*t6*t17/2048-t1/l2*t6*t17/2
     #048-t1*(-6*t30-6*t32-6*t35+6*t38+12*t40-12*t42-24*t45)*t6*t17/2457
     #6
      t6 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = f1(1.D0,z,0.D0)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/4096
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff15 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/1024
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff16 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/1024
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff17 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff18 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum2e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum2e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum2e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum2e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm4e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm4e1int(var,ndim,Tcoeff)

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
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = dlog(-t1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t10 = dlog(z)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t17 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff1 = t1*(t2+2*t3+4*t6-t9-2*t10+t13+t14)*t18*t29/4096-t1/l8*t1
     #8*t29/4096-t1/l3*t18*t29/2048-t1/l2*t18*t29/2048
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = dlog(-t1)
      t8 = f2(1.D0,z)
      t9 = dlog(t8)
      t10 = dlog(z)
      t13 = dlog(1-l4)
      t14 = dlog(l4)
      t17 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff2 = t1*(t2+2*t3+4*t6-t9-2*t10+t13+t14)*t18*t29/4096+t1/(-1+l
     #8)*t18*t29/4096-t1/l3*t18*t29/2048-t1/l2*t18*t29/2048
      t6 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = f1(1.D0,z,0.D0)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/4096
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff4 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/1024
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff5 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/1024
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff6 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff7 = (-1+z)/l2/t4/t8/t10/t13*t15/2048
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum2e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum2e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm4e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm4e2int(var,ndim,Tcoeff)

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

      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/4096
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/4096

CCC Numerator functions

      tcnum1 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum2e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm4e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

