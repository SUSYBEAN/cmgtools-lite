CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV6e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum4e0
      external qgRVnum4e0
      double precision qgRVnum4e1
      external qgRVnum4e1
      double precision qgRVnum4e2
      external qgRVnum4e2
      double precision qgRVnum4e3
      external qgRVnum4e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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

      t3 = dlog(1-z)
      t4 = t3**2
      t7 = dlog(z)
      t9 = dlog(2.D0)
      t11 = t9*t7
      t14 = 0.3141592653589793D1**2
      t16 = t9**2
      t20 = dlog(l1)
      t21 = t20*t7
      t24 = t3*t9
      t27 = t20*t3
      t30 = 36*t7
      t31 = 36*t9
      t32 = 24*t11
      t35 = 72*t3
      t36 = 4*t14
      t37 = 12*t16
      t38 = t20**2
      t44 = t3*l1
      t46 = -42-12*t21*l1-24*t24*l1+24*t27*l1-t30-t31-t32+6*t16*l1+t35+t
     #36-t37+6*t38*l1-12*t38+21*l1+36*t20-36*t44
      t47 = t9*l1
      t49 = t20*l1
      t51 = t7*l1
      t53 = 48*t4
      t54 = t7**2
      t55 = 12*t54
      t56 = t3*t7
      t57 = 48*t56
      t58 = 48*t24
      t63 = t20*t9
      t75 = 18*t47-18*t49+18*t51-t53-t55+t57+t58+6*t54*l1+24*t21-48*t27+
     #24*t63-12*t63*l1+12*t11*l1-24*t56*l1-2*t14*l1+24*t4*l1
      t78 = 1/(-2+l1)
      t80 = 1/l1
      t88 = dlog(l3)
      t97 = 4*t7+4*t9-8*t3+6-2*t88-3*l1-4*t20+4*t44-2*t47+2*t49-2*t51+t8
     #8*l1
      t99 = 1/l3
      t105 = -z3/768-t4*t3/192+7.D0/1024.D0*t7+7.D0/1024.D0*t9+3.D0/512.
     #D0*t11-7.D0/512.D0*t3-t14/1024+3.D0/1024.D0*t16+t16*t9/1536+(t46+t
     #75)*t78*t80/3072-t24*t7/128-t97*t78*t80*t99/1024+3.D0/256.D0*t4+3.
     #D0/1024.D0*t54
      t106 = t88**2
      t115 = 3*t106+t30+t31-t35+t37+42+t53+t55-t36-12*t88*t7+24*t88*t3-1
     #2*t88*t9-18*t88+t32-t58-t57
      t140 = 15.D0/2048.D0+t115*t99/12288+t14*t3/768-t14*t9/1536+t4*t9/1
     #28-t54*t3/256+t54*t9/512-t14*t7/1536+t4*t7/128-t3*t16/256+t16*t7/5
     #12-3.D0/256.D0*t56-3.D0/256.D0*t24+t54*t7/1536
      tncoff1 = t105+t140
      t2 = dlog(1-z)
      t6 = dlog(-z*(-1+l3))
      t8 = dlog(2.D0)
      t10 = dlog(z)
      t12 = dlog(l3)
      t14 = dlog(l1)
      t26 = -6+8*t2-2*t6-4*t8-2*t10+2*t12+4*t14+3*l1+t10*l1-4*t2*l1+t6*l
     #1-t12*l1+2*t8*l1-2*t14*l1
      t31 = 1/l3
      t35 = t10**2
      t37 = 0.3141592653589793D1**2
      t39 = t8**2
      t41 = t6**2
      t43 = t12**2
      t45 = t2**2
      t55 = -42-3*t35+4*t37-12*t39-3*t41-3*t43-48*t45-12*t8*t6-6*t10*t6+
     #24*t2*t6-12*t8*t10
      t73 = -24*t12*t2+12*t12*t8+6*t12*t6+48*t2*t8+24*t2*t10+6*t12*t10+7
     #2*t2-18*t6-36*t8-18*t10+18*t12
      tncoff2 = -t26/(-2+l1)/l1*t31/1024+(t55+t73)*t31/12288
      t1 = z*l1
      t6 = dlog(-(t1+2-l1)*(-2*z+t1-l1))
      t8 = dlog(l1)
      t10 = 2-l1
      t11 = dlog(t10)
      t13 = dlog(2.D0)
      t16 = dlog(1-z)
      t18 = dlog(l3)
      t20 = dlog(z)
      t24 = -1/t10
      t26 = 1/l1
      t31 = t8**2
      t49 = t11**2
      t51 = t16**2
      t55 = 42+12*t31+24*t8*t11+48*t16*t11-12*t8*t6-24*t16*t6-12*t6*t11+
     #6*t20*t6-12*t20*t11+12*t13*t6+12*t49+48*t51+18*t6-36*t8
      t74 = 0.3141592653589793D1**2
      t76 = t13**2
      t78 = t6**2
      t80 = t20**2
      t82 = -36*t11+36*t13-72*t16-24*t16*t20-48*t16*t13+12*t13*t20-12*t8
     #*t20+48*t8*t16-24*t8*t13-24*t13*t11+18*t20-4*t74+12*t76+3*t78+3*t8
     #0
      tncoff3 = -(-6-2*t6+4*t8+4*t11-4*t13+8*t16+2*t18-2*t20)*t24*t26/l3
     #/1024+(t55+t82)*t24*t26/3072
      t7 = z**2
      t10 = l1**2
      t15 = dlog(4*z-4*z*l3+2*l1-4*z*l1+2*t7*l1-t10+2*z*t10-t7*t10)
      t17 = dlog(l1)
      t19 = dlog(l3)
      t21 = dlog(z)
      t23 = 2-l1
      t24 = dlog(t23)
      t26 = dlog(2.D0)
      t29 = dlog(1-z)
      tncoff4 = (6+2*t15-4*t17-2*t19+2*t21-4*t24+4*t26-8*t29)/t23/l1/l3/
     #1024
      t1 = 2-l1
      t3 = -1/t1
      t5 = 1/l1
      t6 = 1/l3
      t10 = dlog(l3)
      t13 = dlog(1-z)
      t14 = 24*t13
      t15 = dlog(z)
      t16 = 12*t15
      t17 = dlog(2.D0)
      t18 = 12*t17
      t22 = t17**2
      t24 = t13**2
      t26 = t15**2
      t39 = dlog(l1)
      t52 = 0.3141592653589793D1**2
      tncoff5 = -t1*t3*t5*t6/1024+7.D0/2048.D0+(18-6*t10-t14+t16+t18)*t6
     #/12288+t22/1024+t24/256+t26/1024-t13*t15/256+t17*t15/512+(-18+6*t1
     #5*l1+t14-t16-t18-12*t13*l1+9*l1+6*t17*l1-6*t39*l1+12*t39)*t3*t5/30
     #72-t13*t17/256-3.D0/512.D0*t13+3.D0/1024.D0*t15+3.D0/1024.D0*t17-t
     #52/3072
      t2 = 1/l3
      t6 = dlog(1-z)
      t8 = dlog(z)
      t10 = dlog(l3)
      t12 = dlog(2.D0)
      t16 = dlog(-z*(-1+l3))
      tncoff6 = -1/l1*t2/1024+(-18+24*t6-6*t8+6*t10-12*t12-6*t16)*t2/122
     #88
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(z)
      t9 = dlog(1-z)
      t11 = z*l1
      t16 = dlog(-(t11+2-l1)*(-2*z+t11-l1))
      t18 = dlog(l1)
      t22 = -1/t3
      t24 = 1/l1
      tncoff7 = (12*t1-12*t4+6*t6-24*t9+6*t16-12*t18+18)*t22*t24/3072+t2
     #2*t24/l3/512
      tncoff8 = -1/(-2+l1)/l1/l3/512
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t8 = dlog(z)
      tncoff9 = 3.D0/2048.D0-t2/512+t4/1024+1/l3/2048+t8/1024+(-6+3*l1)/
     #(-2+l1)/l1/3072
      tncoff10 = -1/l3/2048
      tncoff11 = 1/(-2+l1)/l1/512
      tncoff12 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum4e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum4e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum4e1(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum4e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgRVnum4e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgRVnum4e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgRVnum4e2(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgRVnum4e3(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgRV6e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV6e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum4e0
      external qgRVnum4e0
      double precision qgRVnum4e1
      external qgRVnum4e1
      double precision qgRVnum4e2
      external qgRVnum4e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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

      t1 = 0.3141592653589793D1**2
      t3 = dlog(z)
      t4 = t3**2
      t7 = dlog(1-z)
      t8 = t7**2
      t10 = dlog(2.D0)
      t19 = t10**2
      t22 = dlog(l3)
      t27 = 1/l3
      t40 = dlog(l1)
      t48 = 1/l1
      tncoff1 = -t1/3072+7.D0/2048.D0+t4/1024+t8/256-t7*t10/256-t7*t3/25
     #6-3.D0/512.D0*t7+t10*t3/512+3.D0/1024.D0*t3+t19/1024+3.D0/1024.D0*
     #t10-(t22-2*t10-3-2*t3+4*t7)*t27/2048+(-6+3*l1-4*t3+2*t3*l1+8*t7-4*
     #t7*l1-4*t10+2*t10*l1+4*t40-2*t40*l1)/(-2+l1)*t48/1024+t48*t27/1024
      t1 = dlog(l3)
      t2 = dlog(z)
      t3 = dlog(2.D0)
      t7 = dlog(-z*(-1+l3))
      t9 = dlog(1-z)
      t12 = 1/l3
      t15 = 2-l1
      tncoff2 = -(-t1+3+t2+2*t3+t7-4*t9)*t12/2048-1/l1*t12/1024
      t1 = dlog(z)
      t4 = dlog(1-z)
      t6 = dlog(2.D0)
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = z*l1
      t16 = dlog(-(t11+2-l1)*(-2*z+t11-l1))
      t18 = dlog(l1)
      t22 = -1/t8
      t24 = 1/l1
      tncoff3 = (2*t1-8*t4+4*t6-4*t9+2*t16+6-4*t18)*t22*t24/1024+t22*t24
     #/l3/512
      tncoff4 = -1/(-2+l1)/l1/l3/512
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      tncoff5 = -t2/512+t4/1024+3.D0/2048.D0+t6/1024+1/l3/2048+1/l1/1024
      tncoff6 = -1/l3/2048
      tncoff7 = 1/(-2+l1)/l1/512
      tncoff8 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum4e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum4e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum4e1(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum4e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgRV6e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV6e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum4e0
      external qgRVnum4e0
      double precision qgRVnum4e1
      external qgRVnum4e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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
      t4 = dlog(2.D0)
      t6 = dlog(z)
      tncoff1 = -t2/512+t4/1024+3.D0/2048.D0+t6/1024+1/l3/2048+1/l1/1024
      tncoff2 = -1/l3/2048
      tncoff3 = 1/(-2+l1)/l1/512
      tncoff4 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,1-l1/2,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgRV6e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV6e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum4e0
      external qgRVnum4e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
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

      tncoff1 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV6e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

