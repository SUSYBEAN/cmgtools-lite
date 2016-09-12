CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV2e0int(var,ndim,Tcoeff)

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
      t4 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t10 = dlog(l1)
      t11 = t9*t10
      t14 = t10*t1
      t15 = t14*l1
      t17 = t9*t1
      t21 = t10**2
      t28 = 16-16*t11*l1+8*t15+16*t17*l1-16*t10+8*t21-32*t9-2*t5+16*t1+8
     #*t4-8*l1
      t29 = t21*l1
      t31 = t10*l1
      t34 = t5*l1
      t35 = t9**2
      t43 = t1*l1
      t45 = t4*l1
      t47 = -4*t29+8*t31-32*t17+t34-16*t35*l1+16*t9*l1+32*t11+32*t35-16*
     #t14-8*t43-4*t45
      t50 = 1/(l1-2)
      t52 = 1/l1
      t54 = -1/t8
      t62 = t4*t1
      t66 = t1*t5
      t72 = t10*t5
      t75 = t21*t10
      t80 = t10*t4
      t82 = -96+32*z3-96*t1+12*t66+24*t29-48*t4+12*t5+48*l1-12*t72+96*t1
     #0+16*t75-48*t15-8*t75*l1+48*t80
      t83 = t21*t1
      t104 = -48*t83+8*t62*l1-48*t21+6*t72*l1+24*t83*l1-24*t80*l1-6*t66*
     #l1-48*t31-6*t34+96*t14+48*t43+24*t45-16*z3*l1-16*t62
      t110 = t4**2
      t112 = t5**2
      t131 = 48*t9*t4-12*t9*t5+64*t35*t9-48+96*t9+6*t5-48*t1-24*t4+16*z3
     #-8*t62+96*t17+6*t66-96*t35-96*t35*t1
      tncoff1 = -t1*z3/192-t4*t5/1024+(t28+t47)*t50*t52*t54/256+1.D0/64.
     #D0-t5/512+t1/64+t4/128-z3/192+t62/384+(t82+t104)*t50*t52/3072-t66/
     #512+t110/1536+t112/122880+t131*t54/1536
      t1 = dlog(2.D0)
      t3 = dlog(l1)
      t4 = t3*t1
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = t1**2
      t11 = t7**2
      t13 = 0.3141592653589793D1**2
      t16 = 1-z
      t17 = dlog(t16)
      t18 = t17**2
      t21 = t3**2
      t23 = t3*t7
      t25 = t1*t7
      t33 = -32*t1+32*t4-16+16*t7-32*t9-8*t11+2*t13+16*t3-32*t18+32*t17-
     #8*t21-16*t23+32*t25-32*t17*t3-32*t17*t7+64*t17*t1
      t35 = -1/t6
      t37 = 1/l1
      t61 = 96+128*t9*t1+192*t4*t7+12*t13*t7+192*t1-96*t3-96*t7-192*t4-4
     #8*t21*t7-192*t9*t7+96*t1*t11-16*t11*t7
      t81 = -32*z3+192*t9+48*t11-12*t13-16*t21*t3+48*t21+96*t23-192*t25-
     #24*t1*t13+12*t3*t13-192*t3*t9+96*t21*t1-48*t3*t11
      tncoff2 = -t33*t35*t37/t16/256+(t61+t81)*t35*t37/3072
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t8 = dlog(z)
      t9 = t8**2
      t13 = dlog(2.D0)
      t15 = t2*t13
      t19 = t13**2
      t22 = 0.3141592653589793D1**2
      t27 = 48+96*t3-16*z3-64*t3*t2+8*t9*t8-96*t2+48*t13-96*t15*t8-96*t1
     #5-48*t2*t19+12*t2*t22-6*t13*t22
      t30 = t2*t8
      t32 = t13*t8
      t50 = 96*t3*t13-96*t30+48*t32+24*t19*t8+24*t13*t9-6*t22*t8+96*t3*t
     #8-48*t2*t9+24*t9+48*t8-6*t22+24*t19+8*t19*t13
      t53 = -1/t1
      t57 = dlog(l1)
      t58 = t2*t57
      t63 = t57*t13
      t66 = t57*t8
      t75 = t57**2
      t86 = -16-32*t3+16*t58*l1-16*t15*l1-8*t63*l1-8*t66*l1-16*t30*l1+8*
     #t32*l1+32*t2-16*t13+4*t75*l1-8*t57*l1-t22*l1+16*t3*l1-16*t2*l1-32*
     #t58
      t107 = 16*t63+8*t13*l1+4*t19*l1+16*t66+4*t9*l1+8*t8*l1+32*t15+32*t
     #30-16*t32-8*t9-16*t8+2*t22+8*l1-8*t75+16*t57-8*t19
      tncoff3 = (t27+t50)*t53/1536+(t86+t107)/(l1-2)/l1*t53/256
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 2-l1
      t6 = dlog(t5)
      t8 = dlog(l1)
      t10 = 1-z
      t11 = dlog(t10)
      t13 = t11**2
      t25 = 32*t1+16*t3-16*t6-16*t8-32*t11+32*t13-16*t3*t6+32*t3*t1-32*t
     #1*t6-32*t8*t1-16*t3*t8
      t36 = t3**2
      t38 = t6**2
      t40 = t8**2
      t42 = t1**2
      t44 = 0.3141592653589793D1**2
      t46 = 16+16*t8*t6-32*t3*t11+32*t11*t6+32*t11*t8-64*t1*t11+8*t36+8*
     #t38+8*t40+32*t42-2*t44
      tncoff4 = (t25+t46)/t5/l1/t10/256
      t1 = dlog(2.D0)
      t2 = 0.3141592653589793D1**2
      t5 = 1-z
      t6 = dlog(t5)
      t10 = dlog(l1)
      t12 = t10*l1
      t14 = t1*l1
      t20 = 1/(l1-2)
      t22 = 1/l1
      t24 = -1/t5
      t30 = t10*t1
      t41 = t1**2
      t45 = t10**2
      t50 = -24*t30*l1+24*l1-48*t1+48*t10+6*t2-48-24*t12+24*t14-3*t2*l1+
     #12*t41*l1-24*t41-24*t45+48*t30+12*t45*l1
      t58 = t6**2
      tncoff5 = -t1*t2/1024+(-16*t6-4*l1+8*t1-8*t10+8+4*t12-4*t14+8*t6*l
     #1)*t20*t22*t24/256+t1/128-t2/1024+1.D0/128.D0+t50*t20*t22/3072+t41
     #/256-z3/384+t41*t1/768+(-48*t58+48*t1*t6+48*t6-24*t1+3*t2-24-12*t4
     #1)*t24/1536
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t8 = t4**2
      t10 = t1**2
      t18 = 0.3141592653589793D1**2
      t20 = t6**2
      t22 = 96*t1-48*t4-48*t6+24*t8+96*t10-96*t1*t4-96*t6*t1+48*t6*t4-6*
     #t18+24*t20+48
      t24 = -1/t3
      t26 = 1/l1
      t32 = 1-z
      t33 = dlog(t32)
      tncoff6 = t22*t24*t26/3072-(-16*t1+8*t4+8*t6+16*t33-8)*t24*t26/t32
     #/256
      t1 = dlog(2.D0)
      t2 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t9 = t6**2
      t11 = t2**2
      t13 = 0.3141592653589793D1**2
      t15 = t1**2
      t22 = 24*t1*t2-48*t6*t2+48*t9+12*t11-3*t13+12*t15-48*t6*t1-48*t6+2
     #4+24*t2+24*t1
      t24 = -1/t5
      t27 = dlog(l1)
      tncoff7 = t22*t24/1536+(8*t27-4*t27*l1-8*t6*l1+16*t6-8-8*t2-8*t1+4
     #*t1*l1+4*t2*l1+4*l1)/(l1-2)/l1*t24/256
      t1 = dlog(l1)
      t3 = dlog(z)
      t5 = dlog(2.D0)
      t7 = 2-l1
      t8 = dlog(t7)
      t10 = 1-z
      t11 = dlog(t10)
      tncoff8 = (-8*t1+8*t3+16*t5-8*t8-16*t11+8)/t7/l1/t10/256
      t1 = dlog(2.D0)
      t7 = dlog(l1)
      t13 = 1/(l1-2)
      t15 = 1/l1
      t18 = t1**2
      t20 = 0.3141592653589793D1**2
      t23 = 1-z
      t24 = dlog(t23)
      t28 = -1/t23
      tncoff9 = t1/256+(12*l1-24*t1+12*t1*l1-12*t7*l1+24*t7-24)*t13*t15/
     #3072+t18/512-t20/2048+(-12*t1+24*t24-12)*t28/1536+(-2*l1+4)*t13*t1
     #5*t28/256+1.D0/256.D0
      t1 = 2-l1
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t6 = dlog(l1)
      t10 = -1/t1
      t12 = 1/l1
      tncoff10 = (-24*t2+48*t4-24*t6+24)*t10*t12/3072-t10*t12/(-1+z)/64
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff11 = (12+12*t1+12*t3-24*t6)*t10/1536+(2*l1-4)/(l1-2)/l1*t10/
     #256
      tncoff12 = 1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff13 = t1/512-1/(-1+z)/256+(6*l1-12)/(l1-2)/l1/3072+1.D0/512.D
     #0
      tncoff14 = 1/(l1-2)/l1/256
      tncoff15 = 1/(-1+z)/256
      tncoff16 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum2e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum2e2(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum2e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum2e2(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum2e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum2e3(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum2e3(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum2e4(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)

CCC Sum

      RV2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV2e1int(var,ndim,Tcoeff)

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
      t19 = t13**2
      t23 = -1/t12
      t26 = dlog(l1)
      t27 = t26*t1
      t33 = t1*l1
      t36 = t26**2
      t40 = t26*l1
      t47 = -8*t27*l1+4*t2*l1+16*t27+8*t33-16+16*t26-8*t36+4*t36*l1-8*t4
     #0-t6*l1+8*l1-8*t2+2*t6-16*t1
      t49 = 1/(l1-2)
      t51 = 1/l1
      tncoff1 = 1.D0/128.D0+t2*t1/768+t2/256-t1*t6/1024+t1/128-t6/1024-z
     #3/384+(16*t13*t1-8-4*t2+t6-8*t1+16*t13-16*t19)*t23/512+t47*t49*t51
     #/1024+(-2*t33+4-4*t26+2*t40+4*t13*l1-2*l1+4*t1-8*t13)*t49*t51*t23/
     #128
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = t1**2
      t8 = t4**2
      t10 = 0.3141592653589793D1**2
      t12 = dlog(l1)
      t14 = t12**2
      t22 = 16+32*t1-16*t4+32*t6+8*t8-2*t10-16*t12+8*t14-32*t12*t1+16*t1
     #2*t4-32*t1*t4
      t24 = -1/t3
      t26 = 1/l1
      t32 = 1-z
      t33 = dlog(t32)
      tncoff2 = t22*t24*t26/1024-(-4-8*t1+4*t4+4*t12+8*t33)*t24*t26/t32/
     #128
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t6 = dlog(2.D0)
      t10 = 0.3141592653589793D1**2
      t11 = t6**2
      t14 = t2**2
      t16 = t3**2
      t21 = -16*t2*t3-16*t2*t6+8*t3-t10+4*t11-16*t2+16*t14+4*t16+8*t6+8*
     #t6*t3+8
      t23 = -1/t1
      t30 = dlog(l1)
      tncoff3 = t21*t23/512+(-4*t3+8*t2-4*t6+2*l1+4*t30+2*t6*l1+2*t3*l1-
     #2*t30*l1-4*t2*l1-4)/(l1-2)/l1*t23/128
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 2-l1
      t6 = dlog(t5)
      t8 = dlog(l1)
      t10 = 1-z
      t11 = dlog(t10)
      tncoff4 = (8*t1+4*t3-4*t6-4*t8-8*t11+4)/t5/l1/t10/128
      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t17 = dlog(l1)
      t25 = l1-2
      t26 = 1/t25
      t28 = 1/l1
      tncoff5 = 1.D0/256.D0+t2/512+t1/256-t5/2048+(-4-4*t1+8*t9)*t13/512
     #+(-8-8*t1-4*t17*l1+4*t1*l1+4*l1+8*t17)*t26*t28/1024-t25*t26*t28*t1
     #3/128
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t10 = -1/t3
      t12 = 1/l1
      tncoff6 = (8+16*t1-8*t4-8*t6)*t10*t12/1024-t10*t12/(-1+z)/64
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(z)
      t6 = dlog(2.D0)
      t10 = -1/t1
      tncoff7 = (-8*t2+4*t4+4*t6+4)*t10/512+1/l1*t10/128
      tncoff8 = 1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff9 = 1.D0/512.D0+t1/512-1/(-1+z)/256+(2*l1-4)/(l1-2)/l1/1024
      tncoff10 = 1/(l1-2)/l1/256
      tncoff11 = 1/(-1+z)/256
      tncoff12 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum2e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum2e2(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum2e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum2e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV2e2int(var,ndim,Tcoeff)

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
      t16 = dlog(l1)
      t25 = l1-2
      t26 = 1/t25
      t28 = 1/l1
      tncoff1 = 1.D0/256.D0+t2/512+t1/256-t5/2048+(-2-2*t1+4*t9)*t13/256
     #+(-2*t16*l1+4*t16+2*t1*l1-4*t1+2*l1-4)*t26*t28/512-t25*t26*t28*t13
     #/128
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t10 = -1/t5
      t12 = 1/l1
      tncoff2 = (-4*t1+8*t3-4*t6+4)*t10*t12/512-t10*t12/(-1+z)/64
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      tncoff3 = (2+2*t1+2*t3-4*t6)*t10/256+1/l1*t10/128
      tncoff4 = 1/(l1-2)/l1/(-1+z)/64
      t1 = dlog(2.D0)
      tncoff5 = 1.D0/512.D0+t1/512-1/(-1+z)/256+1/l1/512
      tncoff6 = 1/(l1-2)/l1/256
      tncoff7 = 1/(-1+z)/256
      tncoff8 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum2e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum2e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum2e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum2e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV2e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV2e3int(var,ndim,Tcoeff)

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

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum2e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum2e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum2e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV2e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV2e4int(var,ndim,Tcoeff)

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

      tcnum1 = RVnum2e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV2e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

