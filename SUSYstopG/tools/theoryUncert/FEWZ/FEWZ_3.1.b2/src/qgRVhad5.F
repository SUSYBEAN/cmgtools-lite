CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV5e0int(var,ndim,Tcoeff)

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

      t1 = dlog(l3)
      t3 = dlog(2.D0)
      t6 = dlog(1-z)
      t8 = dlog(z)
      t11 = dlog(l1)
      t13 = t8*l1
      t15 = t6*l1
      t17 = t3*l1
      t19 = t11*l1
      t22 = -2*t1+4*t3-8*t6+4*t8-3*l1-4*t11-2*t13+4*t15-2*t17+2*t19+t1*l
     #1+6
      t24 = 1/(-2+l1)
      t26 = 1/l1
      t27 = 1/l3
      t34 = 0.3141592653589793D1**2
      t37 = t8**2
      t42 = t11*t3
      t45 = t6*t3
      t49 = 36*t8
      t50 = 36*t3
      t52 = 72*t6
      t53 = t11*t6
      t56 = t11**2
      t64 = -42-12*t42*l1-24*t45*l1+18*t13-t49-t50+21*l1+t52-48*t53+24*t
     #42+6*t56*l1-12*t56+36*t11-36*t15+18*t17-18*t19
      t65 = t6*t8
      t66 = 48*t65
      t67 = 48*t45
      t68 = t3*t8
      t69 = 24*t68
      t70 = 4*t34
      t71 = t3**2
      t72 = 12*t71
      t73 = 12*t37
      t74 = t6**2
      t75 = 48*t74
      t78 = t11*t8
      t94 = t66+t67-t69+t70-t72-t73-t75+12*t68*l1-12*t78*l1+24*t53*l1-24
     #*t65*l1+6*t37*l1+24*t74*l1+6*t71*l1-2*t34*l1+24*t78
      t106 = -t22*t24*t26*t27/1024+7.D0/1024.D0*t8+7.D0/1024.D0*t3-7.D0/
     #512.D0*t6+t34*t6/768+t37*t3/512-t34*t8/1536+(t64+t94)*t24*t26/3072
     #-3.D0/256.D0*t65-3.D0/256.D0*t45+3.D0/512.D0*t68-t34/1024+t71*t8/5
     #12+3.D0/1024.D0*t71
      t120 = t1**2
      t122 = -18*t1+t72-t70+t73+t50-t52+t49+t75-12*t1*t8+24*t1*t6-12*t1*
     #t3+t69-t66-t67+42+3*t120
      t140 = 15.D0/2048.D0-t6*t71/256+t71*t3/1536+3.D0/1024.D0*t37+3.D0/
     #256.D0*t74+t122*t27/12288-t37*t6/256-t34*t3/1536+t74*t3/128+t74*t8
     #/128-t45*t8/128+t37*t8/1536-z3/768-t74*t6/192
      tncoff1 = t106+t140
      t1 = 0.3141592653589793D1**2
      t3 = dlog(z)
      t4 = t3**2
      t8 = dlog(-z*(-1+l3))
      t9 = t8**2
      t11 = dlog(l3)
      t12 = t11**2
      t15 = dlog(1-z)
      t16 = t15**2
      t18 = dlog(2.D0)
      t19 = t18**2
      t26 = 4*t1-3*t4-3*t9-3*t12-48*t16-12*t19+18*t11+72*t15-18*t8-36*t1
     #8-18*t3
      t47 = -42+24*t15*t3+48*t15*t18-12*t18*t3+12*t11*t18+6*t11*t3+6*t11
     #*t8-24*t11*t15-12*t18*t8-6*t3*t8+24*t15*t8
      t49 = 1/l3
      t58 = dlog(l1)
      t69 = 2*t11+8*t15-2*t8-4*t18-2*t3+3*l1+4*t58-4*t15*l1+t8*l1-t11*l1
     #+2*t18*l1+t3*l1-2*t58*l1-6
      tncoff2 = (t26+t47)*t49/12288-t69/(-2+l1)/l1*t49/1024
      t1 = dlog(l1)
      t3 = z*l1
      t8 = dlog(-(t3+2-l1)*(-2*z+t3-l1))
      t11 = dlog(1-z)
      t13 = dlog(l3)
      t15 = dlog(z)
      t17 = 2-l1
      t18 = dlog(t17)
      t20 = dlog(2.D0)
      t24 = -1/t17
      t26 = 1/l1
      t33 = t15**2
      t38 = t1**2
      t42 = t18**2
      t48 = t11**2
      t51 = 42+12*t20*t15+3*t33-36*t1-48*t11*t20+12*t38-24*t11*t15+12*t4
     #2+18*t8-72*t11+18*t15-36*t18+48*t48+36*t20
      t52 = 0.3141592653589793D1**2
      t54 = t20**2
      t56 = t8**2
      t82 = -4*t52+12*t54+3*t56-12*t1*t15+48*t1*t11-24*t1*t20-24*t11*t8-
     #12*t15*t18+6*t15*t8+48*t11*t18-24*t20*t18-12*t18*t8+12*t20*t8+24*t
     #1*t18-12*t1*t8
      tncoff3 = -(4*t1-2*t8+8*t11+2*t13-2*t15+4*t18-4*t20-6)*t24*t26/l3/
     #1024+(t51+t82)*t24*t26/3072
      t1 = 2-l1
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t7 = dlog(1-z)
      t9 = dlog(l1)
      t17 = z**2
      t20 = l1**2
      t25 = dlog(4*z-4*z*l3+2*l1-4*z*l1+2*t17*l1-t20+2*z*t20-t17*t20)
      t27 = dlog(z)
      t29 = dlog(l3)
      tncoff4 = (-4*t2+4*t4-8*t7-4*t9+2*t25+2*t27-2*t29+6)/t1/l1/l3/1024
      t1 = dlog(l1)
      t4 = dlog(z)
      t8 = dlog(1-z)
      t11 = dlog(2.D0)
      t14 = 12*t11
      t15 = 24*t8
      t16 = 12*t4
      t20 = -2+l1
      t21 = 1/t20
      t23 = 1/l1
      t34 = 1/l3
      t38 = t4**2
      t40 = 0.3141592653589793D1**2
      t44 = t8**2
      t46 = t11**2
      t48 = dlog(l3)
      tncoff5 = (-6*t1*l1+6*t4*l1-12*t8*l1+6*t11*l1-t14+t15-t16+9*l1+12*
     #t1-18)*t21*t23/3072-t8*t4/256+t11*t4/512-t11*t8/256+t20*t21*t23*t3
     #4/1024+t38/1024-t40/3072+3.D0/1024.D0*t11-3.D0/512.D0*t8+t44/256+t
     #46/1024+(t14-t15-6*t48+t16+18)*t34/12288+3.D0/1024.D0*t4+7.D0/2048
     #.D0
      t2 = 1/l3
      t5 = dlog(l3)
      t7 = dlog(2.D0)
      t10 = dlog(1-z)
      t12 = dlog(z)
      t16 = dlog(-z*(-1+l3))
      tncoff6 = -1/l1*t2/1024+(-18+6*t5-12*t7+24*t10-6*t12-6*t16)*t2/122
     #88
      t1 = -2+l1
      t2 = 1/t1
      t3 = 1/l1
      t8 = dlog(2.D0)
      t11 = dlog(-t1)
      t13 = dlog(z)
      t15 = dlog(l1)
      t17 = z*l1
      t22 = dlog(-(t17+2-l1)*(-2*z+t17-l1))
      t25 = dlog(1-z)
      tncoff7 = t2*t3/l3/512+(18+12*t8-12*t11+6*t13-12*t15+6*t22-24*t25)
     #*t2*t3/3072
      tncoff8 = -1/(-2+l1)/l1/l3/512
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t8 = dlog(1-z)
      tncoff9 = t1/1024+t3/1024+1/l3/2048-t8/512+3.D0/2048.D0+(-6+3*l1)/
     #(-2+l1)/l1/3072
      tncoff10 = -1/l3/2048
      tncoff11 = 1/(-2+l1)/l1/512
      tncoff12 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum4e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum4e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum4e1(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum4e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgRVnum4e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgRVnum4e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgRVnum4e2(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgRVnum4e3(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgRV5e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV5e1int(var,ndim,Tcoeff)

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
      t25 = dlog(l3)
      t27 = 1/l3
      t30 = dlog(l1)
      t48 = 1/l1
      tncoff1 = -t1/3072+7.D0/2048.D0+t4/1024+t8/256-t7*t10/256-t7*t3/25
     #6-3.D0/512.D0*t7+t10*t3/512+3.D0/1024.D0*t3+t19/1024+3.D0/1024.D0*
     #t10-(-3-2*t3+4*t7-2*t10+t25)*t27/2048+(-2*t30*l1+4*t30-6+3*l1-4*t3
     #+2*t3*l1+8*t7-4*t7*l1-4*t10+2*t10*l1)/(-2+l1)*t48/1024+t48*t27/102
     #4
      t1 = dlog(z)
      t2 = dlog(2.D0)
      t6 = dlog(-z*(-1+l3))
      t8 = dlog(1-z)
      t10 = dlog(l3)
      t12 = 1/l3
      t15 = 2-l1
      tncoff2 = -(3+t1+2*t2+t6-4*t8-t10)*t12/2048-1/l1*t12/1024
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
      tncoff3 = (4*t1-4*t4+2*t6-8*t9+2*t16+6-4*t18)*t22*t24/1024+t22*t24
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

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum4e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum4e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum4e1(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum4e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgRV5e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV5e2int(var,ndim,Tcoeff)

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

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum4e0(var,ndim,z,ux,l1/2,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum4e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgRV5e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV5e3int(var,ndim,Tcoeff)

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

      tcnum1 = qgRVnum4e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV5e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

