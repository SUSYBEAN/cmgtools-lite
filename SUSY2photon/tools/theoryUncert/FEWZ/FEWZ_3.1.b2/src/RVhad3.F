CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV3e0int(var,ndim,Tcoeff)

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
      double precision RVnum2e0
      external RVnum2e0
      double precision RVnum2e1
      external RVnum2e1
      double precision RVnum2e2
      external RVnum2e2
      double precision RVnum2e3
      external RVnum2e3
      double precision RVnum2e4
      external RVnum2e4
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

      t1 = dlog(2.D0)
      t2 = 0.3141592653589793D1**2
      t3 = t1*t2
      t5 = t1**2
      t6 = t5*t1
      t18 = 1-z
      t19 = dlog(t18)
      t20 = t19*t1
      t23 = t19**2
      t33 = 48-6*t3+8*t6-16*z3+24*t5+48*t1-6*t2-96*t20-96*t19+96*t23+96*
     #t23*t1-64*t23*t19+12*t19*t2-48*t19*t5
      t35 = -1/t18
      t38 = t2**2
      t40 = t5**2
      t42 = dlog(l1)
      t44 = t5*l1
      t46 = t1*l1
      t48 = t42*l1
      t50 = t42*t1
      t52 = t42**2
      t53 = t52*l1
      t55 = t2*l1
      t59 = -16+16*t42+4*t44+8*t46-8*t48+16*t50+4*t53-t55-8*t52+8*l1-8*t
     #5
      t63 = t19*t42
      t71 = t50*l1
      t77 = -16*t1+2*t2+32*t20-32*t63+16*t23*l1-16*t19*l1+32*t19-32*t23-
     #8*t71-16*t20*l1+16*t63*l1
      t80 = 1/(l1-2)
      t82 = 1/l1
      t86 = t42*t2
      t89 = t42*t5
      t94 = t52*t1
      t98 = t52*t42
      t108 = -96+6*t86*l1-24*t89*l1-12*t86+48*t89-48*t94+12*t3-96*t1-8*t
     #98*l1-16*z3*l1+8*t6*l1+16*t98+12*t2-16*t6
      t125 = -48*t5+24*t44+48*l1+96*t42+48*t46-48*t48+96*t50+24*t53-48*t
     #52-6*t55-48*t71+32*z3-6*t3*l1+24*t94*l1
      tncoff1 = 1.D0/64.D0-t3/512+t6/384-z3/192+t5/128+t1/64-t2/512-t33*
     #t35/1536+t38/122880+t40/1536-(t59+t77)*t80*t82*t35/256+(t108+t125)
     #*t80*t82/3072-t1*z3/192-t5*t2/1024
      t1 = 2-l1
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(l1)
      t6 = t5**2
      t8 = dlog(2.D0)
      t11 = t5*t2
      t13 = t5*t8
      t15 = t8*t2
      t17 = t8**2
      t19 = 0.3141592653589793D1**2
      t22 = 1-z
      t23 = dlog(t22)
      t30 = t23**2
      t33 = 16+8*t3+8*t6+32*t8-16*t2+16*t11-32*t13-32*t15+32*t17-2*t19-1
     #6*t5-64*t23*t8+32*t23*t5+32*t23*t2+32*t30-32*t23
      t35 = -1/t1
      t37 = 1/l1
      t64 = 96-16*t6*t5-16*t3*t2-12*t19-192*t17*t2+96*t8*t3+12*t19*t2-48
     #*t6*t2-48*t5*t3+12*t19*t5+192*t13*t2-24*t8*t19
      t81 = -192*t15+48*t3-192*t13+128*t17*t8-32*z3+96*t6*t8-96*t2+48*t6
     #-96*t5+96*t11+192*t8+192*t17-192*t5*t17
      tncoff2 = t33*t35*t37/t22/256+(t64+t81)*t35*t37/3072
      t1 = dlog(2.D0)
      t2 = dlog(z)
      t3 = t1*t2
      t5 = 1-z
      t6 = dlog(t5)
      t7 = t6**2
      t10 = t2**2
      t21 = 0.3141592653589793D1**2
      t23 = t1**2
      t25 = -48-48*t3+64*t7*t6-8*t10*t2+48*t6*t10+96*t6-48*t1-24*t10-48*
     #t2-96*t7*t2+6*t21-24*t23
      t31 = t6*t2
      t42 = t6*t1
      t50 = -24*t1*t10-8*t23*t1-96*t7+96*t31+6*t21*t2+16*z3-24*t23*t2+6*
     #t1*t21+48*t6*t23+96*t42-96*t7*t1-12*t6*t21+96*t42*t2
      t53 = -1/t5
      t62 = dlog(l1)
      t73 = t62*t2
      t77 = t62*t1
      t79 = t62**2
      t84 = 16+t21*l1-16*t7*l1+16*t3+16*t6*l1+8*t62*l1-8*l1-16*t62-8*t1*
     #l1-4*t10*l1-4*t23*l1-16*t73-8*t2*l1-16*t77-4*t79*l1+16*t42*l1
      t85 = t6*t62
      t107 = -16*t85*l1+8*t77*l1+32*t85-32*t6+16*t1+8*t10+16*t2-2*t21+8*
     #t23+32*t7-32*t31-8*t3*l1+8*t73*l1+16*t31*l1-32*t42+8*t79
      tncoff3 = -(t25+t50)*t53/1536-(t84+t107)/(l1-2)/l1*t53/256
      t1 = dlog(2.D0)
      t2 = t1**2
      t4 = dlog(z)
      t6 = 2-l1
      t7 = dlog(t6)
      t8 = t7**2
      t10 = 1-z
      t11 = dlog(t10)
      t12 = dlog(l1)
      t19 = 0.3141592653589793D1**2
      t27 = -16-32*t2-16*t4-8*t8-32*t11*t12+64*t11*t1-32*t11*t7+2*t19+32
     #*t11*t4+16*t12*t4+32*t12*t1
      t38 = t11**2
      t42 = t4**2
      t44 = t12**2
      t46 = -16*t12*t7-32*t1*t4+32*t7*t1+16*t4*t7+16*t7-32*t1-32*t38+32*
     #t11+16*t12-8*t42-8*t44
      tncoff4 = -(t27+t46)/t6/l1/t10/256
      t1 = 0.3141592653589793D1**2
      t3 = dlog(2.D0)
      t4 = t3**2
      t7 = dlog(l1)
      t8 = t7*l1
      t11 = t3*l1
      t14 = 1-z
      t15 = dlog(t14)
      t22 = 1/(l1-2)
      t24 = 1/l1
      t26 = -1/t14
      t37 = t15**2
      t49 = t7*t3
      t51 = t7**2
      t63 = 12*t4*l1-24*t8+24*l1-48+24*t11+48*t49+12*t51*l1-3*t1*l1-24*t
     #51+6*t1-24*t4-48*t3-24*t49*l1+48*t7
      tncoff5 = 1.D0/128.D0-t1/1024+t4/256+t3/128-(-4*t8+4*l1-8+4*t11-8*
     #t3-8*t15*l1+16*t15+8*t7)*t22*t24*t26/256-z3/384+t4*t3/768-(24-3*t1
     #+12*t4+24*t3-48*t15+48*t37-48*t15*t3)*t26/1536+t63*t22*t24/3072-t3
     #*t1/1024
      t1 = dlog(2.D0)
      t2 = 2-l1
      t3 = dlog(t2)
      t7 = dlog(l1)
      t14 = t3**2
      t16 = t1**2
      t18 = 0.3141592653589793D1**2
      t20 = t7**2
      t22 = -96*t3*t1-48*t3-48*t7+96*t1-96*t7*t1+48*t7*t3+24*t14+96*t16-
     #6*t18+24*t20+48
      t24 = -1/t2
      t26 = 1/l1
      t32 = 1-z
      t33 = dlog(t32)
      tncoff6 = t22*t24*t26/3072+(-8*t3-8*t7+16*t1-16*t33+8)*t24*t26/t32
     #/256
      t1 = 0.3141592653589793D1**2
      t3 = dlog(z)
      t5 = dlog(2.D0)
      t7 = 1-z
      t8 = dlog(t7)
      t10 = t8**2
      t12 = t3**2
      t20 = t5**2
      t22 = -24+3*t1-24*t3-24*t5+48*t8-48*t10-12*t12+48*t8*t3+48*t5*t8-2
     #4*t5*t3-12*t20
      t24 = -1/t7
      t31 = dlog(l1)
      tncoff7 = -t22*t24/1536-(8-4*t5*l1-4*t3*l1+4*t31*l1+8*t8*l1-4*l1+8
     #*t3+8*t5-16*t8-8*t31)/(l1-2)/l1*t24/256
      t1 = dlog(z)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(2.D0)
      t8 = 1-z
      t9 = dlog(t8)
      t11 = dlog(l1)
      tncoff8 = -(-8-8*t1+8*t4-16*t6+16*t9+8*t11)/t3/l1/t8/256
      t1 = dlog(2.D0)
      t2 = t1**2
      t4 = 0.3141592653589793D1**2
      t6 = 1-z
      t7 = dlog(t6)
      t12 = -1/t6
      t15 = dlog(l1)
      t25 = 1/(l1-2)
      t27 = 1/l1
      tncoff9 = t2/512-t4/2048-(-24*t7+12+12*t1)*t12/1536+1.D0/256.D0+(-
     #24-12*t15*l1+12*t1*l1+12*l1-24*t1+24*t15)*t25*t27/3072+t1/256-(-4+
     #2*l1)*t25*t27*t12/256
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t10 = -1/t3
      t12 = 1/l1
      tncoff10 = (48*t1-24*t4-24*t6+24)*t10*t12/3072-t10*t12/(-1+z)/64
      t7 = -1+z
      t8 = 1/t7
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t17 = dlog(-t7)
      tncoff11 = -(-2*l1+4)/(l1-2)/l1*t8/256-(-12*t12-12*t14+24*t17-12)*
     #t8/1536
      tncoff12 = 1/(l1-2)/l1/(-1+z)/64
      t9 = dlog(2.D0)
      tncoff13 = (6*l1-12)/(l1-2)/l1/3072+1.D0/512.D0+t9/512-1/(-1+z)/25
     #6
      tncoff14 = 1/(l1-2)/l1/256
      tncoff15 = 1/(-1+z)/256
      tncoff16 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum2e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum2e2(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum2e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum2e2(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum2e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum2e3(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum2e3(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum2e4(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)

CCC Sum

      RV3e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV3e1int(var,ndim,Tcoeff)

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
      double precision RVnum2e0
      external RVnum2e0
      double precision RVnum2e1
      external RVnum2e1
      double precision RVnum2e2
      external RVnum2e2
      double precision RVnum2e3
      external RVnum2e3
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

      t1 = dlog(2.D0)
      t2 = t1**2
      t6 = 0.3141592653589793D1**2
      t12 = 1-z
      t13 = dlog(t12)
      t14 = t13**2
      t23 = -1/t12
      t27 = dlog(l1)
      t28 = t27**2
      t32 = t27*t1
      t39 = t1*l1
      t43 = t27*l1
      t47 = 8*l1-8*t28+16*t27+2*t6-8*t32*l1-16-8*t2-t6*l1+4*t2*l1+8*t39+
     #4*t28*l1-8*t43+16*t32-16*t1
      t48 = 1/l1
      t51 = 1/(l1-2)
      tncoff1 = 1.D0/128.D0+t2*t1/768+t2/256-t1*t6/1024+t1/128-t6/1024-z
     #3/384-(16*t14-16*t13-t6+8+4*t2-16*t13*t1+8*t1)*t23/512+t47*t48*t51
     #/1024-(8*t13+2*l1+4*t27-4+2*t39-4*t13*l1-2*t43-4*t1)*t51*t48*t23/1
     #28
      t1 = dlog(l1)
      t3 = 0.3141592653589793D1**2
      t5 = dlog(2.D0)
      t6 = t5**2
      t8 = t1**2
      t10 = 2-l1
      t11 = dlog(t10)
      t12 = t11**2
      t22 = 16-16*t1-2*t3+32*t6+8*t8+8*t12-16*t11+32*t5-32*t1*t5-32*t5*t
     #11+16*t1*t11
      t23 = 1/l1
      t26 = -1/t10
      t29 = 1-z
      t30 = dlog(t29)
      tncoff2 = t22*t23*t26/1024+(4-8*t30-4*t1-4*t11+8*t5)*t26*t23/t29/1
     #28
      t1 = 0.3141592653589793D1**2
      t2 = dlog(z)
      t3 = t2**2
      t5 = 1-z
      t6 = dlog(t5)
      t8 = t6**2
      t11 = dlog(2.D0)
      t13 = t11**2
      t21 = -8+t1-4*t3+16*t6-16*t8-8*t2-8*t11-4*t13-8*t11*t2+16*t6*t2+16
     #*t6*t11
      t23 = -1/t5
      t28 = dlog(l1)
      tncoff3 = -t21*t23/512-(4-8*t6-2*l1-4*t28+4*t2+4*t11-2*t11*l1+2*t2
     #8*l1+4*t6*l1-2*t2*l1)/(l1-2)/l1*t23/128
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(z)
      t8 = dlog(l1)
      t10 = 1-z
      t11 = dlog(t10)
      tncoff4 = -(-4-8*t1+4*t4-4*t6+4*t8+8*t11)/t3/l1/t10/128
      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t18 = dlog(l1)
      t25 = 1/l1
      tncoff5 = 1.D0/256.D0+t2/512+t1/256-t5/2048-(4+4*t1-8*t9)*t13/512+
     #(-8-8*t1+4*l1+8*t18-4*t18*l1+4*t1*l1)*t25/(l1-2)/1024-t25*t13/128
      t1 = dlog(l1)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(2.D0)
      t9 = 1/l1
      t12 = -1/t3
      tncoff6 = (8-8*t1-8*t4+16*t6)*t9*t12/1024-t12*t9/(-1+z)/64
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      t13 = 2-l1
      tncoff7 = -(-4-4*t1-4*t3+8*t6)*t10/512+1/l1*t10/128
      tncoff8 = 1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff9 = 1.D0/512.D0+t1/512-1/(-1+z)/256+(-4+2*l1)/l1/(l1-2)/1024
      tncoff10 = 1/(l1-2)/l1/256
      tncoff11 = 1/(-1+z)/256
      tncoff12 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum2e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum2e2(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum2e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum2e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV3e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV3e2int(var,ndim,Tcoeff)

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
      double precision RVnum2e0
      external RVnum2e0
      double precision RVnum2e1
      external RVnum2e1
      double precision RVnum2e2
      external RVnum2e2
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

      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t20 = dlog(l1)
      t28 = 1/l1
      tncoff1 = 1.D0/256.D0+t2/512+t1/256-t5/2048-(2+2*t1-4*t9)*t13/256+
     #(2*t1*l1-4*t1+2*l1-4-2*t20*l1+4*t20)/(l1-2)*t28/512-t28*t13/128
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t10 = -1/t3
      t12 = 1/l1
      tncoff2 = (8*t1-4*t4+4-4*t6)*t10*t12/512-t10*t12/(-1+z)/64
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      t13 = 2-l1
      tncoff3 = -(-2*t1-2*t3-2+4*t6)*t10/256+1/l1*t10/128
      tncoff4 = 1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff5 = 1.D0/512.D0+t1/512-1/(-1+z)/256+1/l1/512
      tncoff6 = 1/(l1-2)/l1/256
      tncoff7 = 1/(-1+z)/256
      tncoff8 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV3e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV3e3int(var,ndim,Tcoeff)

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
      double precision RVnum2e0
      external RVnum2e0
      double precision RVnum2e1
      external RVnum2e1
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

      t1 = dlog(2.D0)
      tncoff1 = 1.D0/512.D0+t1/512-1/(-1+z)/256+1/l1/512
      tncoff2 = 1/(l1-2)/l1/256
      tncoff3 = 1/(-1+z)/256
      tncoff4 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV3e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV3e4int(var,ndim,Tcoeff)

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
      double precision RVnum2e0
      external RVnum2e0
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

      tncoff1 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV3e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

