CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgsq35e0int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
      double precision qgsq3num2e1
      external qgsq3num2e1
      double precision qgsq3num2e2
      external qgsq3num2e2
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

      t1 = dlog(z)
      t2 = 2*t1
      t3 = f3(0.D0,z)
      t5 = dlog(1-t3)
      t7 = dlog(1-z)
      t8 = 4*t7
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t11 = dlog(2.D0)
      t12 = 2*t11
      t13 = dlog(l8)
      t14 = dlog(l1)
      t16 = dlog(1-l4)
      t17 = dlog(l4)
      t30 = -t2+t5+t8-t10+t12+t13+t14+t16+t17+l8-1-2*t11*l8+2*t1*l8-4*t7
     #*l8-t5*l8-t16*l8-t13*l8-t14*l8+t10*l8-t17*l8
      t31 = 1/l8
      t33 = -1+l8
      t34 = 1/t33
      t36 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t37 = 1/t36
      t38 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t39 = 1/t38
      t41 = f1(0.D0,z,1.D0)
      t42 = t37*t39*t41
      t45 = dlog(l3)
      t48 = 1/l3
      t52 = -t33
      t56 = t39*t41
      t60 = dlog(l2)
      t63 = 1/l2
      t87 = t10**2
      t93 = t1**2
      t95 = -3.D0/8192.D0*t16*t1+3.D0/8192.D0*t10*t1+3.D0/16384.D0*t16*t
     #17+3.D0/2048.D0*t7*t11-3.D0/16384.D0*t10*t5+3.D0/4096.D0*t7*t5+3.D
     #0/8192.D0*t5*t11+3.D0/32768.D0*t87-3.D0/2048.D0*t7*t1-3.D0/8192.D0
     #*t10*t11+3.D0/8192.D0*t93
      t104 = 0.3141592653589793D1**2
      t106 = t11**2
      t112 = 3.D0/8192.D0*t1-3.D0/16384.D0*t5-3.D0/4096.D0*t7+3.D0/16384
     #.D0*t10-3.D0/8192.D0*t11-3.D0/16384.D0*t14-3.D0/16384.D0*t16-3.D0/
     #16384.D0*t17-t104/32768+3.D0/8192.D0*t106-3.D0/8192.D0*t5*t1+3.D0/
     #16384.D0*t16*t5
      t114 = t14**2
      t116 = t7**2
      t122 = t16**2
      t126 = t17**2
      t130 = t5**2
      t134 = 3.D0/16384.D0+3.D0/32768.D0*t114+3.D0/2048.D0*t116-3.D0/163
     #84.D0*t14*t10+3.D0/8192.D0*t14*t11+3.D0/32768.D0*t122-3.D0/8192.D0
     #*t14*t1+3.D0/32768.D0*t126+3.D0/8192.D0*t16*t11+3.D0/32768.D0*t130
     #+3.D0/16384.D0*t17*t5
      t159 = 3.D0/4096.D0*t16*t7-3.D0/4096.D0*t10*t7-3.D0/16384.D0*t17*t
     #10+3.D0/8192.D0*t17*t11+3.D0/4096.D0*t17*t7+3.D0/16384.D0*t14*t5-3
     #.D0/8192.D0*t17*t1+3.D0/16384.D0*t17*t14+3.D0/4096.D0*t14*t7-3.D0/
     #16384.D0*t16*t10+3.D0/16384.D0*t16*t14-3.D0/4096.D0*t1*t11
      tncoff1 = 3.D0/16384.D0*t30*t31*t34*t42-3.D0/8192.D0*(t17+t16+2*t4
     #5+t12+t14-t2-1+t8+t5-t10)*t48*t42-3.D0/8192.D0*t52*t34*t31*t48*t37
     #*t56-3.D0/8192.D0*(t17+t16+t14+t12-t10+t8+t5-1-t2+2*t60)*t63*t42-3
     #.D0/8192.D0*t52*t63*t31*t34*t37*t56+(t95+t112+t134+t159)*t37*t56+3
     #.D0/4096.D0*t63*t48*t42
      t2 = dlog(1-z)
      t5 = f2(0.D0,z)
      t6 = dlog(t5)
      t8 = dlog(l1)
      t10 = 1-l8
      t11 = dlog(t10)
      t13 = dlog(l4)
      t16 = dlog(1-l4)
      t18 = dlog(z)
      t21 = f3(0.D0,z)
      t23 = dlog(1-t21)
      t25 = dlog(2.D0)
      t32 = -1/t10
      t34 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t35 = 1/t34
      t36 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t37 = 1/t36
      t39 = f1(0.D0,z,1.D0)
      t40 = t35*t37*t39
      t43 = dlog(l3)
      t45 = 2*t25
      t46 = 2*t18
      t47 = 4*t2
      t49 = 1/l3
      t56 = dlog(l2)
      t59 = 1/l2
      t69 = t6**2
      t71 = t18**2
      t75 = 0.3141592653589793D1**2
      t77 = t25**2
      t79 = t2**2
      t81 = t8**2
      t85 = t16**2
      t89 = 3.D0/32768.D0*t69+3.D0/8192.D0*t71-3.D0/16384.D0*t23-3.D0/81
     #92.D0*t25-t75/32768+3.D0/8192.D0*t77+3.D0/2048.D0*t79+3.D0/32768.D
     #0*t81-3.D0/8192.D0*t8*t18+3.D0/32768.D0*t85+3.D0/8192.D0*t8*t25
      t90 = t13**2
      t92 = t23**2
      t114 = 3.D0/32768.D0*t90+3.D0/32768.D0*t92-3.D0/16384.D0*t8*t6-3.D
     #0/8192.D0*t23*t18+3.D0/16384.D0*t16*t23-3.D0/4096.D0*t6*t2-3.D0/16
     #384.D0*t13*t6+3.D0/4096.D0*t16*t2+3.D0/16384.D0*t13*t23+3.D0/8192.
     #D0*t16*t25-3.D0/8192.D0*t13*t18+3.D0/4096.D0*t13*t2
      t134 = 3.D0/16384.D0+3.D0/8192.D0*t13*t25+3.D0/16384.D0*t13*t8+3.D
     #0/4096.D0*t8*t2-3.D0/16384.D0*t16*t6+3.D0/16384.D0*t8*t23+3.D0/163
     #84.D0*t16*t8-3.D0/4096.D0*t2+3.D0/8192.D0*t18-3.D0/4096.D0*t18*t25
     #-3.D0/8192.D0*t16*t18
      t155 = 3.D0/8192.D0*t6*t18+3.D0/16384.D0*t16*t13+3.D0/2048.D0*t2*t
     #25-3.D0/16384.D0*t6*t23+3.D0/4096.D0*t2*t23+3.D0/8192.D0*t23*t25-3
     #.D0/2048.D0*t2*t18-3.D0/8192.D0*t6*t25-3.D0/16384.D0*t16+3.D0/1638
     #4.D0*t6-3.D0/16384.D0*t8-3.D0/16384.D0*t13
      tncoff2 = 3.D0/16384.D0*(-l8+4*t2*l8-t6*l8+t8*l8+t11*l8+t13*l8+t16
     #*l8-2*t18*l8+t23*l8+2*t25*l8)/l8*t32*t40-3.D0/8192.D0*(t13+t16+2*t
     #43+t45+t8-t46-1+t47+t23-t6)*t49*t40-3.D0/8192.D0*t32*t49*t40-3.D0/
     #8192.D0*(t13+t16+t8+t45-t6+t47+t23-1-t46+2*t56)*t59*t40-3.D0/8192.
     #D0*t59*t32*t40+3.D0/4096.D0*t59*t49*t40+(t89+t114+t134+t155)*t35*t
     #37*t39
      t1 = dlog(2.D0)
      t3 = f3(0.D0,z)
      t5 = dlog(1-t3)
      t6 = dlog(z)
      t8 = dlog(l4)
      t10 = dlog(1-l4)
      t12 = dlog(1-z)
      t14 = dlog(l1)
      t15 = f2(0.D0,z)
      t16 = dlog(t15)
      t17 = 1-l8
      t18 = dlog(t17)
      t19 = dlog(l8)
      t20 = -2*t1-t5+2*t6-t8-t10-4*t12-t14+t16-t18-t19+1
      t21 = 1/l8
      t24 = -1/t17
      t26 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t28 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t31 = f1(0.D0,z,1.D0)
      t32 = 1/t26/t28*t31
      tncoff3 = 3.D0/16384.D0*t20*t21*t24*t32+3.D0/8192.D0*t24*t21/l3*t3
     #2+3.D0/8192.D0/l2*t21*t24*t32
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l3)
      t6 = 2-l3
      t7 = dlog(t6)
      t10 = dlog(1-z)
      t12 = l3/2
      t13 = f3(t12,z)
      t15 = dlog(1-t13)
      t16 = dlog(l1)
      t17 = dlog(z)
      t19 = f2(t12,z)
      t20 = dlog(t19)
      t23 = -1/t6
      t25 = 1/l3
      t27 = Dsquar311(z,l1,0.D0,t12,l4)
      t28 = 1/t27
      t29 = s23r2(z,1.D0,1.D0,t12,l4)
      t30 = 1/t29
      t32 = f1(t12,z,1.D0)
      t33 = t28*t30*t32
      tncoff4 = -3.D0/4096.D0*(t1+t3+2*t4+2*t7+4*t10-1+t15+t16-2*t17-t20
     #)*t23*t25*t33+3.D0/4096.D0*t23*t25/l8*t33+3.D0/2048.D0*t23*t30*t28
     #/l2*t25*t32
      t2 = dlog(1-l4)
      t3 = dlog(l3)
      t6 = dlog(1-z)
      t8 = dlog(l1)
      t9 = dlog(z)
      t11 = l3/2
      t12 = f3(t11,z)
      t14 = dlog(1-t12)
      t15 = f2(t11,z)
      t16 = dlog(t15)
      t17 = 2-l3
      t18 = dlog(t17)
      t20 = dlog(l4)
      t23 = -1/t17
      t25 = 1/l3
      t27 = Dsquar311(z,l1,0.D0,t11,l4)
      t28 = 1/t27
      t29 = s23r2(z,1.D0,1.D0,t11,l4)
      t30 = 1/t29
      t32 = f1(t11,z,1.D0)
      t33 = t28*t30*t32
      tncoff5 = 3.D0/4096.D0*(-t2-2*t3-4*t6-t8+2*t9-t14+t16+1-2*t18-t20)
     #*t23*t25*t33-3.D0/4096.D0*t23*t25/(-1+l8)*t33+3.D0/2048.D0*t23*t30
     #*t28/l2*t25*t32
      t9 = l3/2
      t10 = Dsquar311(z,l1,0.D0,t9,l4)
      t13 = s23r2(z,1.D0,1.D0,t9,l4)
      t15 = f1(t9,z,1.D0)
      tncoff6 = 3.D0/4096.D0/(-2+l3)/l3/l8/(-1+l8)/t10/t13*t15
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l1)
      t5 = dlog(2.D0)
      t7 = l1*l2
      t8 = 1-t7
      t9 = dlog(t8)
      t10 = 1-l2
      t11 = dlog(t10)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t14 = f3(0.D0,z)
      t17 = dlog(1-t14+t7*t14)
      t18 = dlog(z)
      t21 = dlog(1-z)
      t23 = dlog(l2)
      t25 = t1+t3+t4+2*t5+t9+t11-t13+t17-1-2*t18+4*t21+2*t23
      t26 = 1/l2
      t28 = Dsquar311(z,l1,l2,0.D0,l4)
      t30 = s23r2(z,t8,t10,0.D0,l4)
      t33 = f1(0.D0,z,t8)
      t34 = 1/t28/t30*t33
      tncoff7 = 3.D0/8192.D0*t25*t26*t34-3.D0/8192.D0*t26/l8*t34-3.D0/40
     #96.D0*t26/l3*t34
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l1)
      t5 = dlog(2.D0)
      t7 = l1*l2
      t8 = 1-t7
      t9 = dlog(t8)
      t10 = 1-l2
      t11 = dlog(t10)
      t12 = f2(0.D0,z)
      t13 = dlog(t12)
      t14 = f3(0.D0,z)
      t17 = dlog(1-t14+t7*t14)
      t18 = dlog(z)
      t21 = dlog(1-z)
      t23 = dlog(l2)
      t25 = t1+t3+t4+2*t5+t9+t11-t13+t17-1-2*t18+4*t21+2*t23
      t26 = 1/l2
      t28 = Dsquar311(z,l1,l2,0.D0,l4)
      t30 = s23r2(z,t8,t10,0.D0,l4)
      t33 = f1(0.D0,z,t8)
      t34 = 1/t28/t30*t33
      tncoff8 = 3.D0/8192.D0*t25*t26*t34+3.D0/8192.D0*t26/(-1+l8)*t34-3.
     #D0/4096.D0*t26/l3*t34
      t7 = Dsquar311(z,l1,l2,0.D0,l4)
      t10 = 1-l1*l2
      t12 = s23r2(z,t10,1-l2,0.D0,l4)
      t15 = f1(0.D0,z,t10)
      tncoff9 = -3.D0/8192.D0/l2/l8/(-1+l8)/t7/t12*t15
      t1 = l3/2
      t3 = 1-l1*l2
      t4 = f1(t1,z,t3)
      t6 = s23r2(z,t3,1-l2,t1,l4)
      t12 = Dsquar311(z,l1,l2,t1,l4)
      tncoff10 = -3.D0/2048.D0*t4/t6/(-2+l3)/t12/l2/l3
      t1 = l3/2
      t3 = 1-l1*l2
      t4 = f1(t1,z,t3)
      t6 = s23r2(z,t3,1-l2,t1,l4)
      t12 = Dsquar311(z,l1,l2,t1,l4)
      tncoff11 = -3.D0/2048.D0*t4/t6/(-2+l3)/t12/l2/l3
      t2 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t7 = f1(0.D0,z,1.D0)
      t8 = 1/t5*t7
      t19 = dlog(l1)
      t21 = dlog(l4)
      t24 = dlog(1-l4)
      t26 = dlog(2.D0)
      t28 = f3(0.D0,z)
      t30 = dlog(1-t28)
      t33 = dlog(1-z)
      t35 = f2(0.D0,z)
      t36 = dlog(t35)
      t38 = dlog(z)
      tncoff12 = 3.D0/16384.D0/l8*t3*t8+3.D0/8192.D0/l3*t3*t8+3.D0/8192.
     #D0/l2*t3*t8+(-3.D0/16384.D0*t19-3.D0/16384.D0*t21-3.D0/16384.D0*t2
     #4-3.D0/8192.D0*t26-3.D0/16384.D0*t30-3.D0/4096.D0*t33+3.D0/16384.D
     #0*t36+3.D0/8192.D0*t38+3.D0/16384.D0)*t3*t8
      t3 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      t9 = 1/t6*t8
      t20 = dlog(l1)
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t27 = dlog(2.D0)
      t29 = f3(0.D0,z)
      t31 = dlog(1-t29)
      t34 = dlog(1-z)
      t36 = f2(0.D0,z)
      t37 = dlog(t36)
      t39 = dlog(z)
      tncoff13 = -3.D0/16384.D0/(-1+l8)*t4*t9+3.D0/8192.D0/l3*t4*t9+3.D0
     #/8192.D0/l2*t4*t9+(-3.D0/16384.D0*t20-3.D0/16384.D0*t22-3.D0/16384
     #.D0*t25-3.D0/8192.D0*t27-3.D0/16384.D0*t31-3.D0/4096.D0*t34+3.D0/1
     #6384.D0*t37+3.D0/8192.D0*t39+3.D0/16384.D0)*t4*t9
      t5 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,z,1.D0)
      tncoff14 = 3.D0/16384.D0/l8/(-1+l8)/t5/t7*t10
      t5 = l3/2
      t6 = Dsquar311(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,1.D0,1.D0,t5,l4)
      t11 = f1(t5,z,1.D0)
      tncoff15 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t5 = l3/2
      t6 = Dsquar311(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,1.D0,1.D0,t5,l4)
      t11 = f1(t5,z,1.D0)
      tncoff16 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,0.D0,l4)
      t10 = f1(0.D0,z,t6)
      tncoff17 = -3.D0/8192.D0/l2/t2/t8*t10
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,0.D0,l4)
      t10 = f1(0.D0,z,t6)
      tncoff18 = -3.D0/8192.D0/l2/t2/t8*t10
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff19 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff20 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l
     #8,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,
     #0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,l3/2,l4,l5,l6,l7,
     #1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgsq3num2e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,
     #0.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgsq3num2e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,
     #1.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgsq3num2e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgsq3num2e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgsq35e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq35e1int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
      double precision qgsq3num2e1
      external qgsq3num2e1
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

      t1 = dlog(l1)
      t3 = dlog(l4)
      t6 = dlog(1-l4)
      t8 = dlog(2.D0)
      t10 = f3(0.D0,z)
      t12 = dlog(1-t10)
      t15 = dlog(1-z)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t20 = dlog(z)
      t23 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t24 = 1/t23
      t26 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t27 = 1/t26
      t28 = f1(0.D0,z,1.D0)
      t29 = t27*t28
      t31 = 1-l8
      tncoff1 = (-3.D0/16384.D0*t1-3.D0/16384.D0*t3-3.D0/16384.D0*t6-3.D
     #0/8192.D0*t8-3.D0/16384.D0*t12-3.D0/4096.D0*t15+3.D0/16384.D0*t18+
     #3.D0/8192.D0*t20+3.D0/16384.D0)*t24*t29+3.D0/16384.D0/l8*t24*t27*t
     #28+3.D0/8192.D0/l3*t24*t29+3.D0/8192.D0/l2*t24*t29
      t3 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t4 = 1/t3
      t6 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t8 = f1(0.D0,z,1.D0)
      t9 = 1/t6*t8
      t20 = dlog(l1)
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t27 = dlog(2.D0)
      t29 = f3(0.D0,z)
      t31 = dlog(1-t29)
      t34 = dlog(1-z)
      t36 = f2(0.D0,z)
      t37 = dlog(t36)
      t39 = dlog(z)
      tncoff2 = -3.D0/16384.D0/(-1+l8)*t4*t9+3.D0/8192.D0/l3*t4*t9+3.D0/
     #8192.D0/l2*t4*t9+(-3.D0/16384.D0*t20-3.D0/16384.D0*t22-3.D0/16384.
     #D0*t25-3.D0/8192.D0*t27-3.D0/16384.D0*t31-3.D0/4096.D0*t34+3.D0/16
     #384.D0*t37+3.D0/8192.D0*t39+3.D0/16384.D0)*t4*t9
      t5 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t7 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t10 = f1(0.D0,z,1.D0)
      tncoff3 = 3.D0/16384.D0/l8/(-1+l8)/t5/t7*t10
      t5 = l3/2
      t6 = Dsquar311(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,1.D0,1.D0,t5,l4)
      t11 = f1(t5,z,1.D0)
      tncoff4 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t5 = l3/2
      t6 = Dsquar311(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,1.D0,1.D0,t5,l4)
      t11 = f1(t5,z,1.D0)
      tncoff5 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,0.D0,l4)
      t10 = f1(0.D0,z,t6)
      tncoff6 = -3.D0/8192.D0/l2/t2/t8*t10
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t8 = s23r2(z,t6,1-l2,0.D0,l4)
      t10 = f1(0.D0,z,t6)
      tncoff7 = -3.D0/8192.D0/l2/t2/t8*t10
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff8 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff9 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num2e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgsq35e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq35e2int(var,ndim,Tcoeff)

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
      double precision qgsq3num2e0
      external qgsq3num2e0
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

      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff1 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = f1(0.D0,z,1.D0)
      tncoff2 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num2e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgsq35e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq35e3int(var,ndim,Tcoeff)

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

      qgsq35e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

