CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgsq31e0int(var,ndim,Tcoeff)

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
      double precision qgsq3num1e0
      external qgsq3num1e0
      double precision qgsq3num1e1
      external qgsq3num1e1
      double precision qgsq3num1e2
      external qgsq3num1e2
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

      t2 = dlog(1-z)
      t5 = dlog(z)
      t8 = f2(0.D0,z)
      t9 = dlog(t8)
      t11 = dlog(l8)
      t13 = dlog(2.D0)
      t17 = dlog(1-l4)
      t19 = dlog(l4)
      t20 = dlog(l1)
      t21 = 2*t13
      t24 = 4*t2
      t25 = 2*t5
      t26 = -4*t2*l8+2*t5*l8+t9*l8-t11*l8-2*t13*l8-t17*l8+t19+t17+t11+t2
     #0+t21-t9-1-t20*l8-t19*l8+t24-t25+l8
      t27 = 1/l8
      t32 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t33 = 1/t32
      t34 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t35 = 1/t34
      t37 = f1(0.D0,z,0.D0)
      t38 = t33*t35*t37
      t41 = dlog(l3)
      t44 = 1/l3
      t51 = dlog(l2)
      t54 = 1/l2
      t58 = 0.3141592653589793D1**2
      t60 = t13**2
      t62 = t2**2
      t64 = t20**2
      t66 = t19**2
      t68 = t17**2
      t70 = t9**2
      t72 = t5**2
      t94 = -t58/32768+3.D0/8192.D0*t60+3.D0/2048.D0*t62+3.D0/32768.D0*t
     #64+3.D0/32768.D0*t66+3.D0/32768.D0*t68+3.D0/32768.D0*t70+3.D0/8192
     #.D0*t72+3.D0/4096.D0*t17*t2+3.D0/8192.D0*t17*t13+3.D0/2048.D0*t2*t
     #13+3.D0/4096.D0*t20*t2+3.D0/8192.D0*t9*t5+3.D0/16384.D0*t17*t19+3.
     #D0/16384.D0*t17*t20-3.D0/16384.D0*t17*t9-3.D0/16384.D0*t20*t9-3.D0
     #/4096.D0*t5*t13
      t124 = 3.D0/16384.D0-3.D0/4096.D0*t9*t2-3.D0/8192.D0*t20*t5+3.D0/8
     #192.D0*t20*t13-3.D0/16384.D0*t19*t9+3.D0/16384.D0*t19*t20+3.D0/409
     #6.D0*t19*t2-3.D0/8192.D0*t19*t5+3.D0/8192.D0*t19*t13-3.D0/2048.D0*
     #t2*t5-3.D0/8192.D0*t17*t5-3.D0/8192.D0*t9*t13-3.D0/4096.D0*t2+3.D0
     #/8192.D0*t5+3.D0/16384.D0*t9-3.D0/8192.D0*t13-3.D0/16384.D0*t20-3.
     #D0/16384.D0*t19-3.D0/16384.D0*t17
      tncoff1 = 3.D0/16384.D0*t26*t27/(-1+l8)*t38-3.D0/8192.D0*(2*t41+t2
     #0+t21-t9-1+t24-t25+t19+t17)*t44*t38+3.D0/8192.D0*t44*t27*t38-3.D0/
     #8192.D0*(t17+t19+2*t51+t20-t9+t24-1+t21-t25)*t54*t38+(t94+t124)*t3
     #3*t35*t37+3.D0/8192.D0*t54*t27*t38+3.D0/4096.D0*t54*t44*t38
      t2 = dlog(1-z)
      t3 = dlog(z)
      t6 = f2(0.D0,z)
      t7 = dlog(t6)
      t8 = dlog(2.D0)
      t15 = dlog(l1)
      t18 = dlog(l4)
      t22 = dlog(1-l4)
      t36 = 0.3141592653589793D1**2
      t41 = t8**2
      t43 = -3.D0/2048.D0*t2*t3-3.D0/8192.D0*t7*t8-3.D0/4096.D0*t3*t8-3.
     #D0/4096.D0*t7*t2+3.D0/8192.D0*t15*t8+3.D0/16384.D0*t18*t15+3.D0/16
     #384.D0*t22*t18+3.D0/2048.D0*t2*t8-3.D0/16384.D0*t22*t7+3.D0/8192.D
     #0*t22*t8+3.D0/8192.D0*t7*t3+3.D0/4096.D0*t22*t2-3.D0/8192.D0*t8-t3
     #6/32768-3.D0/16384.D0*t18-3.D0/16384.D0*t22-3.D0/16384.D0*t15+3.D0
     #/8192.D0*t41
      t46 = t15**2
      t49 = t22**2
      t51 = t2**2
      t53 = t18**2
      t73 = t3**2
      t75 = t7**2
      t77 = 3.D0/16384.D0+3.D0/16384.D0*t7-3.D0/4096.D0*t2+3.D0/32768.D0
     #*t46+3.D0/8192.D0*t3+3.D0/32768.D0*t49+3.D0/2048.D0*t51+3.D0/32768
     #.D0*t53-3.D0/8192.D0*t18*t3+3.D0/8192.D0*t18*t8-3.D0/8192.D0*t22*t
     #3-3.D0/16384.D0*t15*t7-3.D0/8192.D0*t15*t3+3.D0/4096.D0*t15*t2+3.D
     #0/16384.D0*t22*t15+3.D0/4096.D0*t18*t2-3.D0/16384.D0*t18*t7+3.D0/8
     #192.D0*t73+3.D0/32768.D0*t75
      t79 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t80 = 1/t79
      t82 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t83 = 1/t82
      t84 = f1(0.D0,z,0.D0)
      t96 = 1-l8
      t97 = dlog(t96)
      t104 = -1/t96
      t107 = t80*t83*t84
      t110 = dlog(l3)
      t112 = 2*t8
      t113 = 4*t2
      t114 = 2*t3
      t116 = 1/l3
      t123 = dlog(l2)
      t126 = 1/l2
      tncoff2 = (t43+t77)*t80*t83*t84+3.D0/16384.D0*(-l8-2*t3*l8+2*t8*l8
     #+4*t2*l8-t7*l8+t15*l8+t22*l8+t97*l8+t18*l8)/l8*t104*t107-3.D0/8192
     #.D0*(2*t110+t15+t112-t7-1+t113-t114+t18+t22)*t116*t107-3.D0/8192.D
     #0*t104*t116*t107-3.D0/8192.D0*(t22+t18+2*t123+t15-t7+t113-1+t112-t
     #114)*t126*t107-3.D0/8192.D0*t126*t104*t107+3.D0/4096.D0*t126*t116*
     #t107
      t1 = dlog(l1)
      t2 = dlog(2.D0)
      t5 = dlog(1-z)
      t7 = dlog(z)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t12 = dlog(1-l4)
      t13 = dlog(l8)
      t14 = 1-l8
      t15 = dlog(t14)
      t16 = dlog(l4)
      t18 = 1/l8
      t21 = -1/t14
      t23 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t25 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t28 = f1(0.D0,z,0.D0)
      t29 = 1/t23/t25*t28
      tncoff3 = 3.D0/16384.D0*(-t1-2*t2-4*t5+2*t7+t10-t12-t13-t15-t16+1)
     #*t18*t21*t29+3.D0/8192.D0*t21/l3*t18*t29+3.D0/8192.D0/l2*t18*t21*t
     #29
      t1 = dlog(l3)
      t3 = 2-l3
      t4 = dlog(t3)
      t6 = l3/2
      t7 = f2(t6,z)
      t8 = dlog(t7)
      t9 = dlog(l1)
      t11 = dlog(1-z)
      t13 = dlog(z)
      t15 = dlog(l4)
      t17 = dlog(1-l4)
      t20 = -1/t3
      t22 = 1/l3
      t24 = Dsquar321(z,l1,0.D0,t6,l4)
      t25 = 1/t24
      t26 = s23r2(z,0.D0,0.D0,t6,l4)
      t27 = 1/t26
      t29 = f1(t6,z,0.D0)
      t30 = t25*t27*t29
      tncoff4 = -3.D0/4096.D0*(2*t1+2*t4-t8+t9+4*t11-2*t13-1+t15+t17)*t2
     #0*t22*t30+3.D0/4096.D0*t20*t22/l8*t30+3.D0/2048.D0*t20*t27*t25/l2*
     #t22*t29
      t1 = dlog(l3)
      t3 = 2-l3
      t4 = dlog(t3)
      t6 = l3/2
      t7 = f2(t6,z)
      t8 = dlog(t7)
      t9 = dlog(l1)
      t11 = dlog(1-z)
      t13 = dlog(z)
      t15 = dlog(l4)
      t17 = dlog(1-l4)
      t20 = -1/t3
      t22 = 1/l3
      t24 = Dsquar321(z,l1,0.D0,t6,l4)
      t25 = 1/t24
      t26 = s23r2(z,0.D0,0.D0,t6,l4)
      t27 = 1/t26
      t29 = f1(t6,z,0.D0)
      t30 = t25*t27*t29
      tncoff5 = -3.D0/4096.D0*(2*t1+2*t4-t8+t9+4*t11-2*t13-1+t15+t17)*t2
     #0*t22*t30-3.D0/4096.D0*t20*t22/(-1+l8)*t30+3.D0/2048.D0*t20*t27*t2
     #5/l2*t22*t29
      t9 = l3/2
      t10 = Dsquar321(z,l1,0.D0,t9,l4)
      t13 = s23r2(z,0.D0,0.D0,t9,l4)
      t15 = f1(t9,z,0.D0)
      tncoff6 = 3.D0/4096.D0/(-2+l3)/l3/l8/(-1+l8)/t10/t13*t15
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t7 = dlog(2.D0)
      t9 = l1*l2
      t11 = dlog(1-t9)
      t13 = dlog(1-l2)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-z)
      t19 = dlog(z)
      t21 = f3(0.D0,z)
      t24 = dlog(1-t9*t21)
      t25 = t2+t3+2*t4+t6+2*t7+t11+t13-t15+4*t17-1-2*t19+t24
      t26 = 1/l2
      t28 = Dsquar321(z,l1,l2,0.D0,l4)
      t30 = s23r2(z,t9,l2,0.D0,l4)
      t33 = f1(0.D0,z,t9)
      t34 = 1/t28/t30*t33
      tncoff7 = 3.D0/8192.D0*t25*t26*t34-3.D0/8192.D0*t26/l8*t34-3.D0/40
     #96.D0*t26/l3*t34
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t7 = dlog(2.D0)
      t9 = l1*l2
      t11 = dlog(1-t9)
      t13 = dlog(1-l2)
      t14 = f2(0.D0,z)
      t15 = dlog(t14)
      t17 = dlog(1-z)
      t19 = dlog(z)
      t21 = f3(0.D0,z)
      t24 = dlog(1-t9*t21)
      t25 = t2+t3+2*t4+t6+2*t7+t11+t13-t15+4*t17-1-2*t19+t24
      t26 = 1/l2
      t28 = Dsquar321(z,l1,l2,0.D0,l4)
      t30 = s23r2(z,t9,l2,0.D0,l4)
      t33 = f1(0.D0,z,t9)
      t34 = 1/t28/t30*t33
      tncoff8 = 3.D0/8192.D0*t25*t26*t34+3.D0/8192.D0*t26/(-1+l8)*t34-3.
     #D0/4096.D0*t26/l3*t34
      t7 = Dsquar321(z,l1,l2,0.D0,l4)
      t9 = l1*l2
      t10 = s23r2(z,t9,l2,0.D0,l4)
      t13 = f1(0.D0,z,t9)
      tncoff9 = -3.D0/8192.D0/l2/l8/(-1+l8)/t7/t10*t13
      t1 = l3/2
      t2 = l1*l2
      t3 = f1(t1,z,t2)
      t4 = s23r2(z,t2,l2,t1,l4)
      t10 = Dsquar321(z,l1,l2,t1,l4)
      tncoff10 = -3.D0/2048.D0*t3/t4/(-2+l3)/t10/l2/l3
      t1 = l3/2
      t2 = l1*l2
      t3 = f1(t1,z,t2)
      t4 = s23r2(z,t2,l2,t1,l4)
      t10 = Dsquar321(z,l1,l2,t1,l4)
      tncoff11 = -3.D0/2048.D0*t3/t4/(-2+l3)/t10/l2/l3
      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t7 = f1(0.D0,z,0.D0)
      t8 = 1/t5*t7
      t20 = dlog(1-z)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t25 = dlog(2.D0)
      t27 = dlog(z)
      t29 = dlog(l1)
      t32 = dlog(1-l4)
      t34 = dlog(l4)
      tncoff12 = 3.D0/16384.D0/l8*t3*t8+3.D0/8192.D0/l3*t3*t8+3.D0/8192.
     #D0/l2*t3*t8+(-3.D0/4096.D0*t20+3.D0/16384.D0*t23-3.D0/8192.D0*t25+
     #3.D0/8192.D0*t27-3.D0/16384.D0*t29-3.D0/16384.D0*t32-3.D0/16384.D0
     #*t34+3.D0/16384.D0)*t3*t8
      t2 = dlog(1-z)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t11 = dlog(l1)
      t14 = dlog(1-l4)
      t16 = dlog(l4)
      t19 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t24 = f1(0.D0,z,0.D0)
      t25 = 1/t22*t24
      tncoff13 = (-3.D0/4096.D0*t2+3.D0/16384.D0*t5-3.D0/8192.D0*t7+3.D0
     #/8192.D0*t9-3.D0/16384.D0*t11-3.D0/16384.D0*t14-3.D0/16384.D0*t16+
     #3.D0/16384.D0)*t20*t25-3.D0/16384.D0/(-1+l8)*t20*t25+3.D0/8192.D0/
     #l3*t20*t25+3.D0/8192.D0/l2*t20*t25
      t5 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t7 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = f1(0.D0,z,0.D0)
      tncoff14 = 3.D0/16384.D0/l8/(-1+l8)/t5/t7*t10
      t5 = l3/2
      t6 = Dsquar321(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t11 = f1(t5,z,0.D0)
      tncoff15 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t5 = l3/2
      t6 = Dsquar321(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t11 = f1(t5,z,0.D0)
      tncoff16 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = f1(0.D0,z,t5)
      tncoff17 = -3.D0/8192.D0/l2/t2/t6*t8
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = f1(0.D0,z,t5)
      tncoff18 = -3.D0/8192.D0/l2/t2/t6*t8
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff19 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff20 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qgsq3num1e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qgsq3num1e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qgsq3num1e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qgsq3num1e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      qgsq31e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq31e1int(var,ndim,Tcoeff)

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
      double precision qgsq3num1e0
      external qgsq3num1e0
      double precision qgsq3num1e1
      external qgsq3num1e1
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

      t2 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = 1/t2
      t5 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t7 = f1(0.D0,z,0.D0)
      t8 = 1/t5*t7
      t20 = dlog(1-z)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t25 = dlog(2.D0)
      t27 = dlog(z)
      t29 = dlog(l1)
      t32 = dlog(1-l4)
      t34 = dlog(l4)
      tncoff1 = 3.D0/16384.D0/l8*t3*t8+3.D0/8192.D0/l3*t3*t8+3.D0/8192.D
     #0/l2*t3*t8+(-3.D0/4096.D0*t20+3.D0/16384.D0*t23-3.D0/8192.D0*t25+3
     #.D0/8192.D0*t27-3.D0/16384.D0*t29-3.D0/16384.D0*t32-3.D0/16384.D0*
     #t34+3.D0/16384.D0)*t3*t8
      t2 = dlog(1-z)
      t4 = f2(0.D0,z)
      t5 = dlog(t4)
      t7 = dlog(2.D0)
      t9 = dlog(z)
      t11 = dlog(l1)
      t14 = dlog(1-l4)
      t16 = dlog(l4)
      t19 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t24 = f1(0.D0,z,0.D0)
      t25 = 1/t22*t24
      tncoff2 = (-3.D0/4096.D0*t2+3.D0/16384.D0*t5-3.D0/8192.D0*t7+3.D0/
     #8192.D0*t9-3.D0/16384.D0*t11-3.D0/16384.D0*t14-3.D0/16384.D0*t16+3
     #.D0/16384.D0)*t20*t25-3.D0/16384.D0/(-1+l8)*t20*t25+3.D0/8192.D0/l
     #3*t20*t25+3.D0/8192.D0/l2*t20*t25
      t5 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t7 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t10 = f1(0.D0,z,0.D0)
      tncoff3 = 3.D0/16384.D0/l8/(-1+l8)/t5/t7*t10
      t5 = l3/2
      t6 = Dsquar321(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t11 = f1(t5,z,0.D0)
      tncoff4 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t5 = l3/2
      t6 = Dsquar321(z,l1,0.D0,t5,l4)
      t8 = s23r2(z,0.D0,0.D0,t5,l4)
      t11 = f1(t5,z,0.D0)
      tncoff5 = 3.D0/4096.D0/(-2+l3)/l3/t6/t8*t11
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = f1(0.D0,z,t5)
      tncoff6 = -3.D0/8192.D0/l2/t2/t6*t8
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = f1(0.D0,z,t5)
      tncoff7 = -3.D0/8192.D0/l2/t2/t6*t8
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff8 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff9 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgsq3num1e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgsq3num1e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      qgsq31e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq31e2int(var,ndim,Tcoeff)

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
      double precision qgsq3num1e0
      external qgsq3num1e0
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

      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff1 = 3.D0/16384.D0/t1/t3*t6
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = f1(0.D0,z,0.D0)
      tncoff2 = 3.D0/16384.D0/t1/t3*t6

CCC Numerator functions

      tcnum1 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgsq3num1e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qgsq31e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsq31e3int(var,ndim,Tcoeff)

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

      qgsq31e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

