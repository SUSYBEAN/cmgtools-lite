CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqirm2e0int(var,ndim,Tcoeff)

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
      t2 = dlog(z)
      t6 = dlog(-t1)
      t9 = dlog(2.D0)
      t12 = dlog(l1)
      t14 = dlog(l4)
      t16 = dlog(l8)
      t19 = dlog(1-l4)
      t21 = 4*t6
      t22 = 2*t2
      t23 = f2(1.D0,z)
      t24 = dlog(t23)
      t26 = 2*t9
      t27 = -2*t2*l8+4*t6*l8+2*t9*l8+t12*l8+t14*l8+t16*l8+t19*l8-t19-t21
     #-t12+t22+t24-t24*l8-t14-t26-t16
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
      t67 = 0.3141592653589793D1**2
      t68 = t9**2
      t70 = t24**2
      t72 = t2**2
      t74 = t6**2
      t92 = 48*t9*t6-t67+12*t68+3*t70+12*t72+48*t74-24*t9*t2+12*t19*t9+2
     #4*t19*t6-12*t19*t2+24*t12*t6-6*t12*t24-12*t12*t2-24*t6*t24
      t105 = t19**2
      t107 = t12**2
      t109 = t14**2
      t123 = -48*t6*t2+6*t19*t14-6*t14*t24+12*t12*t9-12*t14*t2+12*t24*t2
     #+3*t105+3*t107+3*t109+24*t14*t6-6*t19*t24-12*t9*t24+12*t14*t9+6*t1
     #4*t12+6*t19*t12
      t129 = dlog(l2)
      t133 = 1/l2
      t138 = t1*t133
      tncoff1 = t1*t27/(-1+l8)*t31*t35*t37*t40*t42*t44/2048+(t14+2*t49+t
     #26+t21-t24-t22+t12+t19)*t1*t54*t58/1024-t1*t31*t54*t58/1024-t1*(t9
     #2+t123)*t35*t58/12288+t1*(t14+t19+2*t129+t26+t21-t24-t22+t12)*t133
     #*t35*t58/1024-t138*t31*t35*t58/1024-t138*t54*t58/512
      t1 = -1+z
      t2 = 1-l8
      t3 = dlog(t2)
      t5 = dlog(z)
      t8 = dlog(l4)
      t10 = dlog(l1)
      t12 = dlog(2.D0)
      t16 = dlog(-t1)
      t19 = f2(1.D0,z)
      t20 = dlog(t19)
      t23 = dlog(1-l4)
      t28 = -1/t2
      t32 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t33 = 1/t32
      t34 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t35 = 1/t34
      t37 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t38 = 1/t37
      t39 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t40 = 1/t39
      t42 = f1(1.D0,z,0.D0)
      t47 = 1/l2
      t48 = t1*t47
      t50 = 1/l3*t33
      t54 = t35*t38*t40*t42
      t57 = dlog(l3)
      t59 = 2*t12
      t60 = 4*t16
      t61 = 2*t5
      t71 = dlog(l2)
      t95 = 0.3141592653589793D1**2
      t100 = t12**2
      t102 = t20**2
      t104 = t5**2
      t106 = -6*t23*t20+48*t12*t16+6*t23*t10+6*t8*t10+12*t8*t12-12*t10*t
     #5-6*t10*t20+24*t10*t16-t95-12*t23*t5+24*t23*t16+12*t100+3*t102+12*
     #t104
      t109 = t16**2
      t131 = t10**2
      t133 = t23**2
      t135 = t8**2
      t137 = 12*t23*t12+48*t109-24*t12*t5-12*t12*t20-24*t16*t20-48*t16*t
     #5+6*t23*t8-6*t8*t20+12*t10*t12-12*t8*t5+12*t20*t5+24*t8*t16+3*t131
     #+3*t133+3*t135
      tncoff2 = t1*(-t3*l8+2*t5*l8-t8*l8-t10*l8-2*t12*l8-4*t16*l8+t20*l8
     #-t23*l8)*t28/l8*t33*t35*t38*t40*t42/2048-t48*t50*t54/512+(t8+2*t57
     #+t59+t60-t20-t61+t10+t23)*t1*t50*t54/1024+t1*t28*t50*t54/1024+t1*(
     #t8+t23+2*t71+t59+t60-t20-t61+t10)*t47*t33*t54/1024-t1*(t106+t137)*
     #t33*t54/12288+t48*t28*t33*t54/1024
      t1 = -1+z
      t2 = f2(1.D0,z)
      t3 = dlog(t2)
      t4 = dlog(l1)
      t5 = dlog(l4)
      t6 = 1-l8
      t7 = dlog(t6)
      t9 = dlog(-t1)
      t11 = dlog(2.D0)
      t13 = dlog(z)
      t16 = dlog(1-l4)
      t17 = dlog(l8)
      t21 = -1/t6
      t22 = 1/l8
      t23 = t21*t22
      t25 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t27 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t30 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t32 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t35 = f1(1.D0,z,0.D0)
      t37 = 1/t25/t27/t30/t32*t35
      tncoff3 = t1*(-t3+t4+t5+t7+4*t9+2*t11-2*t13+t16+t17)*t23*t37/2048-
     #t1*t21*t22/l3*t37/1024-t1/l2*t23*t37/1024
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(l3)
      t5 = dlog(l1)
      t6 = 2-l3
      t7 = dlog(t6)
      t10 = 1-l3/2
      t11 = f2(t10,z)
      t12 = dlog(t11)
      t13 = dlog(z)
      t16 = dlog(-t1)
      t19 = dlog(1-l4)
      t23 = -1/t6
      t24 = 1/l3
      t27 = Dsquar321(z,l1,0.D0,t10,l4)
      t29 = s2vr2(z,0.D0,0.D0,t10,l4)
      t31 = 1/t27/t29
      t32 = s23r2(z,0.D0,0.D0,t10,l4)
      t34 = s1vr2(z,0.D0,0.D0,t10,l4)
      t36 = 1/t32/t34
      t37 = f1(t10,z,0.D0)
      t39 = t31*t36*t37
      t42 = t1*t23
      tncoff4 = t1*(t2+2*t3+t5+2*t7-t12-2*t13+4*t16+t19)*t23*t24*t39/512
     #-t42*t24/l8*t39/512-t42*t36*t31/l2*t24*t37/256
      t1 = -1+z
      t2 = l3-2
      t3 = 1/t2
      t4 = t1*t3
      t6 = 1-l3/2
      t7 = s1vr2(z,0.D0,0.D0,t6,l4)
      t9 = s23r2(z,0.D0,0.D0,t6,l4)
      t11 = 1/t7/t9
      t13 = s2vr2(z,0.D0,0.D0,t6,l4)
      t15 = Dsquar321(z,l1,0.D0,t6,l4)
      t17 = 1/t13/t15
      t19 = 1/l3
      t21 = f1(t6,z,0.D0)
      t26 = dlog(l4)
      t27 = dlog(l3)
      t29 = dlog(l1)
      t31 = dlog(-t2)
      t33 = f2(t6,z)
      t34 = dlog(t33)
      t35 = dlog(z)
      t38 = dlog(-t1)
      t41 = dlog(1-l4)
      t47 = t17*t11*t21
      tncoff5 = -t4*t11*t17/l2*t19*t21/256+t1*(t26+2*t27+t29+2*t31-t34-2
     #*t35+4*t38+t41)*t3*t19*t47/512+t4*t19/(-1+l8)*t47/512
      t13 = 1-l3/2
      t14 = Dsquar321(z,l1,0.D0,t13,l4)
      t16 = s2vr2(z,0.D0,0.D0,t13,l4)
      t19 = s23r2(z,0.D0,0.D0,t13,l4)
      t21 = s1vr2(z,0.D0,0.D0,t13,l4)
      t24 = f1(t13,z,0.D0)
      tncoff6 = -(-1+z)/(l3-2)/l3/l8/(-1+l8)/t14/t16/t19/t21*t24/512
      t1 = -1+z
      t2 = dlog(l4)
      t4 = dlog(1-l4)
      t5 = dlog(l2)
      t7 = l1*l2
      t8 = f3(1.D0,z)
      t11 = dlog(1-t7*t8)
      t13 = dlog(1-t7)
      t14 = f2(1.D0,z)
      t15 = dlog(t14)
      t17 = dlog(-t1)
      t19 = dlog(z)
      t21 = dlog(l1)
      t22 = dlog(2.D0)
      t25 = dlog(1-l2)
      t26 = t2+t4+2*t5+t11+t13-t15+4*t17-2*t19+t21+2*t22+t25
      t28 = 1/l2
      t29 = Dsquar321(z,l1,l2,1.D0,l4)
      t30 = 1/t29
      t33 = s2vr2(z,t7,l2,1.D0,l4)
      t35 = s23r2(z,t7,l2,1.D0,l4)
      t38 = s1vr2(z,t7,l2,1.D0,l4)
      t40 = f1(1.D0,z,t7)
      t42 = 1/t33/t35/t38*t40
      t45 = t1*t28
      tncoff7 = -t1*t26*t28*t30*t42/1024+t45/l8*t30*t42/1024+t45/l3*t30*
     #t42/512
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
      t22 = dlog(l4)
      t24 = dlog(1-l4)
      t25 = dlog(l2)
      t27 = f3(1.D0,z)
      t30 = dlog(1-t9*t27)
      t32 = dlog(1-t9)
      t33 = f2(1.D0,z)
      t34 = dlog(t33)
      t36 = dlog(-t1)
      t38 = dlog(z)
      t40 = dlog(l1)
      t41 = dlog(2.D0)
      t44 = dlog(1-l2)
      t45 = t22+t24+2*t25+t30+t32-t34+4*t36-2*t38+t40+2*t41+t44
      tncoff8 = t3/l3*t6*t19/512-t1*t45*t2*t6*t19/1024-t3/(-1+l8)*t6*t19
     #/1024
      t9 = Dsquar321(z,l1,l2,1.D0,l4)
      t11 = l1*l2
      t12 = s2vr2(z,t11,l2,1.D0,l4)
      t15 = s23r2(z,t11,l2,1.D0,l4)
      t17 = s1vr2(z,t11,l2,1.D0,l4)
      t20 = f1(1.D0,z,t11)
      tncoff9 = (-1+z)/l2/(-1+l8)/l8/t9/t12/t15/t17*t20/1024
      t3 = 1-l3/2
      t4 = l1*l2
      t5 = f1(t3,z,t4)
      t7 = s1vr2(z,t4,l2,t3,l4)
      t9 = s23r2(z,t4,l2,t3,l4)
      t13 = s2vr2(z,t4,l2,t3,l4)
      t18 = Dsquar321(z,l1,l2,t3,l4)
      tncoff10 = (-1+z)*t5/t7/t9/t13/(l3-2)/t18/l2/l3/256
      t3 = 1-l3/2
      t4 = l1*l2
      t5 = f1(t3,z,t4)
      t7 = s1vr2(z,t4,l2,t3,l4)
      t9 = s23r2(z,t4,l2,t3,l4)
      t13 = s2vr2(z,t4,l2,t3,l4)
      t18 = Dsquar321(z,l1,l2,t3,l4)
      tncoff11 = (-1+z)*t5/t7/t9/t13/(l3-2)/t18/l2/l3/256
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
      t32 = dlog(2.D0)
      t34 = f2(1.D0,z)
      t35 = dlog(t34)
      t37 = dlog(z)
      t40 = dlog(-t1)
      t42 = dlog(l4)
      t44 = dlog(l1)
      t47 = dlog(1-l4)
      tncoff12 = -t1/l8*t10*t12*t15*t17*t19/2048-t1/l3*t10*t29/1024-t1*(
     #-12*t32+6*t35+12*t37-24*t40-6*t42-6*t44-6*t47)*t10*t29/12288-t1/l2
     #*t10*t29/1024
      t1 = -1+z
      t5 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t8 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t10 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t13 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t15 = f1(1.D0,z,0.D0)
      t17 = 1/t8/t10/t13*t15
      t30 = dlog(2.D0)
      t32 = f2(1.D0,z)
      t33 = dlog(t32)
      t35 = dlog(z)
      t38 = dlog(-t1)
      t40 = dlog(l4)
      t42 = dlog(l1)
      t45 = dlog(1-l4)
      tncoff13 = t1/(-1+l8)*t6*t17/2048-t1/l3*t6*t17/1024-t1/l2*t6*t17/1
     #024-t1*(-12*t30+6*t33+12*t35-24*t38-6*t40-6*t42-6*t45)*t6*t17/1228
     #8
      t6 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = f1(1.D0,z,0.D0)
      tncoff14 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff15 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/512
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff16 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/512
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff17 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff18 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,
     #1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qiqirmnum1e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qiqirmnum1e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qiqirmnum1e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qiqirmnum1e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.
     #D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qiqirm2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tn
     #coff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm2e1int(var,ndim,Tcoeff)

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
      t2 = dlog(l4)
      t3 = dlog(2.D0)
      t6 = dlog(1-l4)
      t8 = dlog(-t1)
      t10 = dlog(l1)
      t11 = dlog(z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff1 = t1*(t2+2*t3+t6+4*t8+t10-2*t11-t14)*t18*t29/2048-t1/l8*t1
     #8*t29/2048-t1/l3*t18*t29/1024-t1/l2*t18*t29/1024
      t1 = -1+z
      t2 = dlog(l4)
      t3 = dlog(2.D0)
      t6 = dlog(1-l4)
      t8 = dlog(-t1)
      t10 = dlog(l1)
      t11 = dlog(z)
      t13 = f2(1.D0,z)
      t14 = dlog(t13)
      t17 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t18 = 1/t17
      t20 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t22 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t25 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t29 = 1/t20/t22/t25*t27
      tncoff2 = t1*(t2+2*t3+t6+4*t8+t10-2*t11-t14)*t18*t29/2048+t1/(-1+l
     #8)*t18*t29/2048-t1/l3*t18*t29/1024-t1/l2*t18*t29/1024
      t6 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t12 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t15 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t17 = f1(1.D0,z,0.D0)
      tncoff3 = -(-1+z)/(-1+l8)/l8/t6/t10/t12/t15*t17/2048
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff4 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/512
      t7 = 1-l3/2
      t8 = Dsquar321(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,0.D0,0.D0,t7,l4)
      t14 = s23r2(z,0.D0,0.D0,t7,l4)
      t17 = s1vr2(z,0.D0,0.D0,t7,l4)
      t19 = f1(t7,z,0.D0)
      tncoff5 = -(-1+z)/(l3-2)/l3/t8/t12/t14/t17*t19/512
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff6 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t4 = Dsquar321(z,l1,l2,1.D0,l4)
      t7 = l1*l2
      t8 = s2vr2(z,t7,l2,1.D0,l4)
      t10 = s23r2(z,t7,l2,1.D0,l4)
      t13 = s1vr2(z,t7,l2,1.D0,l4)
      t15 = f1(1.D0,z,t7)
      tncoff7 = (-1+z)/l2/t4/t8/t10/t13*t15/1024
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qiqirmnum1e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qiqirmnum1e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qiqirm2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqirm2e2int(var,ndim,Tcoeff)

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

      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t8 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t10 = s1vr2(z,0.D0,0.D0,1.D0,l4)
      t13 = f1(1.D0,z,0.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqirmnum1e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qiqirm2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

