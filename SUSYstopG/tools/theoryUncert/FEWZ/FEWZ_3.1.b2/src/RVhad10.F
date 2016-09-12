CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV10e0int(var,ndim,Tcoeff)

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
      double precision RVnum7e0
      external RVnum7e0
      double precision RVnum7e1
      external RVnum7e1
      double precision RVnum7e2
      external RVnum7e2
      double precision RVnum7e3
      external RVnum7e3
      double precision RVnum7e4
      external RVnum7e4
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

      t2 = 0.3141592653589793D1**2
      t3 = t2**2
      t6 = 1-z
      t7 = dlog(t6)
      t10 = 24*t7
      t11 = 4*t2
      t13 = t7**2
      t16 = 12*t13
      t19 = -1/t6
      t22 = 8*t2
      t23 = 16*z3
      t24 = dlog(l3)
      t25 = t24**2
      t32 = 1/l3
      t35 = 3*t25
      t36 = 12*t24
      t43 = dlog(l1)
      t44 = t43**2
      t52 = 1/l1
      t55 = 24*t43
      t56 = 12*t44
      tncoff1 = -1.D0/96.D0+z3/288-t3/15360+t2/576+(4*t7*t2-t10-t11-8*z3
     #-4*t13*t7+t16+24)*t19/2304+(t22+t23-6*t25+24*t24+t25*t24-48-4*t24*
     #t2)*t32/9216+(-t10-t11+t16+t35-t36+24+12*t7*t24)*t32*t19/4608-(-t2
     #2-t23+48-8*t44*t43-48*t43+24*t44+8*t43*t2)*t52/4608-(t10+t11-t16-2
     #4+t55-t56-24*t7*t43)*t52*t19/2304-(-t11+t35-t36+24-t55+t56+12*t24*
     #t43)*t32*t52/4608-(2*t7+t24-2+2*t43)*t32*t52*t19/384
      t1 = dlog(l3)
      t3 = 0.3141592653589793D1**2
      t5 = t1**2
      t12 = 1/l3
      t15 = 12*t1
      t16 = 4*t3
      t17 = 3*t5
      t18 = 1-z
      t19 = dlog(t18)
      t21 = t19**2
      t28 = -1/t18
      t31 = dlog(l1)
      t32 = t31**2
      t39 = 1/l1
      tncoff2 = (48-24*t1-8*t3+6*t5-16*z3-t5*t1+4*t1*t3)*t12/9216+(-24+t
     #15+t16-t17+24*t19-12*t21-12*t19*t1)*t12*t28/4608-(-24+t15+t16-t17-
     #12*t32+24*t31-12*t1*t31)*t12*t39/4608-(2-t1-2*t19-2*t31)*t12*t39*t
     #28/384
      t2 = dlog(1-l1)
      t3 = 0.3141592653589793D1**2
      t6 = dlog(l1)
      t7 = t6**2
      t12 = t2**2
      t15 = t6*t2
      t26 = -48-4*t2*t3+12*t7*t2-8*t6*t3+6*t6*t12-24*t15+16*z3+8*t7*t6-2
     #4*t7-6*t12+48*t6+24*t2+8*t3+t12*t2
      t27 = 1/l1
      t30 = 1-z
      t31 = dlog(t30)
      t36 = 12*t15
      t37 = 12*t7
      t38 = 3*t12
      t39 = 24*t6
      t40 = 12*t2
      t41 = 4*t3
      t42 = t31**2
      t45 = 24+12*t31*t2+24*t31*t6+t36+t37+t38-t39-t40-t41+12*t42-24*t31
      t48 = -1/t30
      t51 = dlog(l3)
      t56 = t51**2
      t59 = -24-t36-12*t51*t6-6*t51*t2-t37-t38+t39+t40+t41-3*t56+12*t51
      t60 = 1/l3
      tncoff3 = -t26*t27/4608-t45*t27*t48/2304-t59*t60*t27/4608-(2-2*t6-
     #t2-2*t31-t51)*t60*t27*t48/384
      t1 = dlog(l3)
      t3 = dlog(1-l1)
      t6 = dlog(l1)
      t7 = t6**2
      t13 = t1**2
      t18 = 0.3141592653589793D1**2
      t20 = t3**2
      t22 = 24+6*t1*t3+12*t7+12*t6*t3+12*t1*t6+3*t13-12*t3-12*t1-24*t6-4
     #*t18+3*t20
      t23 = 1/l3
      t25 = 1/l1
      t29 = 1-z
      t30 = dlog(t29)
      tncoff4 = -t22*t23*t25/4608+(-2+t3+t1+2*t6+2*t30)*t23*t25/t29/384
      t1 = dlog(z)
      t2 = 24*t1
      t3 = 1-z
      t4 = dlog(t3)
      t5 = t4**2
      t6 = 12*t5
      t7 = t1**2
      t10 = 24*t4
      t14 = 24*t4*t1
      t17 = 0.3141592653589793D1**2
      t18 = 4*t17
      t19 = 12*t7
      t27 = -24-t2-t6-4*t7*t1+t10+12*t4*t7+t14-12*t5*t1+t18-t19+8*z3+4*t
     #17*t1-4*t4*t17+4*t5*t4
      t29 = -1/t3
      t32 = dlog(l3)
      t33 = t32**2
      t40 = -3*t33+12*t32*t1-24-t2-t6+t10+t14+t18-t19-12*t4*t32+12*t32
      t41 = 1/l3
      t45 = dlog(l1)
      t46 = t45**2
      t53 = 12*t46-24*t45*t1+24*t4*t45-24*t45+24+t2+t6-t10-t14-t18+t19
      t54 = 1/l1
      tncoff5 = t27*t29/2304+t40*t41*t29/4608-t53*t54*t29/2304-(-2*t45+2
     #+2*t1-2*t4-t32)*t41*t54*t29/384
      t1 = dlog(z)
      t3 = dlog(l3)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = 0.3141592653589793D1**2
      t12 = t3**2
      t18 = t1**2
      t20 = t6**2
      t22 = 24+24*t1-12*t3-24*t6-24*t1*t6-4*t10+3*t12+12*t6*t3-12*t1*t3+
     #12*t18+12*t20
      t23 = 1/l3
      t26 = -1/t5
      t31 = dlog(l1)
      tncoff6 = t22*t23*t26/4608-(-2-2*t1+t3+2*t6+2*t31)*t23/l1*t26/384
      t1 = dlog(l1)
      t2 = dlog(z)
      t6 = dlog(1-l1)
      t11 = 0.3141592653589793D1**2
      t13 = t2**2
      t15 = t6**2
      t18 = 1-z
      t19 = dlog(t18)
      t22 = t19**2
      t30 = t1**2
      t33 = -24+24*t1*t2-12*t1*t6+12*t6*t2+4*t11-12*t13-3*t15+24*t1+24*t
     #19-24*t2-12*t22-12*t19*t6+24*t19*t2-24*t19*t1-12*t30+12*t6
      t34 = 1/l1
      t37 = -1/t18
      t43 = dlog(l3)
      tncoff7 = -t33*t34*t37/2304-(-2+2*t1+2*t19-2*t2+t6+t43)/l3*t34*t37
     #/384
      t1 = dlog(l3)
      t2 = 1-z
      t3 = dlog(t2)
      t6 = dlog(1-l1)
      t7 = dlog(z)
      t9 = dlog(l1)
      tncoff8 = (2-t1-2*t3-t6+2*t7-2*t9)/l3/l1/t2/384
      t1 = 0.3141592653589793D1**2
      t5 = 1-z
      t6 = dlog(t5)
      t7 = t6**2
      t9 = 12*t6
      t12 = -1/t5
      t15 = 4*t1
      t16 = dlog(l3)
      t18 = t16**2
      t21 = 1/l3
      t24 = 6*t16
      t29 = dlog(l1)
      t30 = t29**2
      t34 = 1/l1
      tncoff9 = -1.D0/192.D0+t1/1152+z3/576+(12-2*t1+6*t7-t9)*t12/2304+(
     #-24+t15+12*t16-3*t18)*t21/9216+(12-t9-t24)*t21*t12/4608-(12*t30-24
     #*t29+24-t15)*t34/4608-(12*t29-12+12*t6)*t34*t12/2304-(-12*t29+12-t
     #24)*t21*t34/4608+t21*t34*t12/384
      t1 = dlog(l3)
      t3 = 0.3141592653589793D1**2
      t5 = t1**2
      t8 = 1/l3
      t11 = 6*t1
      t12 = 1-z
      t13 = dlog(t12)
      t18 = -1/t12
      t21 = dlog(l1)
      t25 = 1/l1
      tncoff10 = (24-12*t1-4*t3+3*t5)*t8/9216+(-12+t11+12*t13)*t8*t18/46
     #08-(-12+t11+12*t21)*t8*t25/4608-t8*t25*t18/384
      t1 = dlog(l1)
      t2 = t1**2
      t5 = dlog(1-l1)
      t8 = t5**2
      t12 = 0.3141592653589793D1**2
      t15 = 1/l1
      t18 = 6*t5
      t19 = 12*t1
      t20 = 1-z
      t21 = dlog(t20)
      t26 = -1/t20
      t29 = dlog(l3)
      t32 = 1/l3
      tncoff11 = -(-24-12*t2+12*t5+24*t1-3*t8-12*t5*t1+4*t12)*t15/4608-(
     #12-t18-t19-12*t21)*t15*t26/2304-(-12+t18+t19+6*t29)*t32*t15/4608-t
     #32*t15*t26/384
      t2 = dlog(1-l1)
      t4 = dlog(l3)
      t6 = dlog(l1)
      t9 = 1/l3
      t11 = 1/l1
      tncoff12 = -(12-6*t2-6*t4-12*t6)*t9*t11/4608+t9*t11/(-1+z)/384
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = dlog(z)
      t8 = 0.3141592653589793D1**2
      t10 = t5**2
      t12 = 12*t5
      t13 = 12*t2
      t16 = -1/t1
      t19 = dlog(l3)
      t22 = 1/l3
      t26 = dlog(l1)
      t29 = 1/l1
      tncoff13 = (-12-6*t3+12*t2*t5+2*t8-6*t10-t12+t13)*t16/2304+(-12-t1
     #2+t13+6*t19)*t22*t16/4608-(12+12*t5-12*t2-12*t26)*t29*t16/2304-t22
     #*t29*t16/384
      t1 = dlog(z)
      t3 = 1-z
      t4 = dlog(t3)
      t6 = dlog(l3)
      t9 = 1/l3
      t12 = -1/t3
      tncoff14 = (12+12*t1-12*t4-6*t6)*t9*t12/4608+t9/l1*t12/384
      t2 = dlog(1-l1)
      t4 = dlog(z)
      t6 = dlog(l1)
      t8 = 1-z
      t9 = dlog(t8)
      t12 = 1/l1
      t15 = -1/t8
      tncoff15 = -(-12+6*t2-12*t4+12*t6+12*t9)*t12*t15/2304+1/l3*t12*t15
     #/384
      tncoff16 = -1/l3/l1/(-1+z)/384
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l3)
      t14 = 1/l3
      t19 = dlog(l1)
      t22 = 1/l1
      tncoff17 = t1/2304-1.D0/384.D0+(-6*t4+6)*t8/2304+(6*t11-12)*t14/92
     #16+t14*t8/768-(-12*t19+12)*t22/4608+t22*t8/384-t14*t22/768
      t1 = dlog(l3)
      t4 = 1/l3
      tncoff18 = (12-6*t1)*t4/9216-t4/(-1+z)/768+t4/l1/768
      t2 = dlog(1-l1)
      t4 = dlog(l1)
      t7 = 1/l1
      tncoff19 = -(-12+6*t2+12*t4)*t7/4608-t7/(-1+z)/384+1/l3*t7/768
      tncoff20 = -1/l3/l1/768
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff21 = (-6-6*t1+6*t3)*t7/2304-1/l3*t7/768-1/l1*t7/384
      tncoff22 = 1/l3/(-1+z)/768
      tncoff23 = 1/l1/(-1+z)/384
      tncoff24 = -1.D0/768.D0+1/(-1+z)/768-1/l3/1536-1/l1/768
      tncoff25 = 1/l3/1536
      tncoff26 = 1/l1/768
      tncoff27 = -1/(-1+z)/768
      tncoff28 = -1.D0/1536.D0

CCC Numerator functions

      tcnum1 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum7e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum7e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum7e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,999
     #.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum7e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum7e1(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum7e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum7e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum7e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum7e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,99
     #9.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum7e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum7e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum7e2(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum7e2(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = RVnum7e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = RVnum7e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = RVnum7e2(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = RVnum7e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = RVnum7e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = RVnum7e3(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = RVnum7e3(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = RVnum7e4(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff28)

CCC Sum

      t15 = tcnum19*tncoff19+tcnum15*tncoff15+tcnum20*tncoff20+tcnum6*tn
     #coff6+tcnum5*tncoff5+tcnum3*tncoff3+tcnum12*tncoff12+tcnum8*tncoff
     #8+tcnum14*tncoff14+tcnum17*tncoff17+tcnum11*tncoff11+tcnum23*tncof
     #f23+tcnum2*tncoff2+tcnum7*tncoff7
      t30 = tcnum25*tncoff25+tcnum24*tncoff24+tcnum28*tncoff28+tcnum18*t
     #ncoff18+tcnum26*tncoff26+tcnum21*tncoff21+tcnum22*tncoff22+tcnum13
     #*tncoff13+tcnum27*tncoff27+tcnum4*tncoff4+tcnum16*tncoff16+tcnum9*
     #tncoff9+tcnum1*tncoff1+tcnum10*tncoff10
      RV10e0int = t15+t30

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV10e1int(var,ndim,Tcoeff)

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
      double precision RVnum7e0
      external RVnum7e0
      double precision RVnum7e1
      external RVnum7e1
      double precision RVnum7e2
      external RVnum7e2
      double precision RVnum7e3
      external RVnum7e3
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
      t4 = 1-z
      t5 = dlog(t4)
      t6 = t5**2
      t12 = -1/t4
      t15 = dlog(l3)
      t17 = 4*t1
      t18 = t15**2
      t21 = 1/l3
      t29 = dlog(l1)
      t31 = t29**2
      t34 = 1/l1
      tncoff1 = -1.D0/192.D0+t1/1152+z3/576-(-6*t6+2*t1-12+12*t5)*t12/23
     #04-(24-12*t15-t17+3*t18)*t21/9216-(t15-2+2*t5)*t21*t12/768+(24*t29
     #-12*t31-24+t17)*t34/4608-(-2+2*t29+2*t5)*t34*t12/384+(2*t29-2+t15)
     #*t34*t21/768+t21*t34*t12/384
      t1 = 0.3141592653589793D1**2
      t3 = dlog(l3)
      t5 = t3**2
      t8 = 1/l3
      t11 = 1-z
      t12 = dlog(t11)
      t17 = -1/t11
      t20 = dlog(l1)
      t23 = 1/l1
      tncoff2 = -(-24+4*t1+12*t3-3*t5)*t8/9216-(-t3+2-2*t12)*t8*t17/768+
     #(-2*t20+2-t3)*t23*t8/768-t8*t23*t17/384
      t1 = dlog(l1)
      t3 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = dlog(1-l1)
      t12 = t8**2
      t15 = 1/l1
      t18 = 2*t1
      t19 = 1-z
      t20 = dlog(t19)
      t25 = -1/t19
      t28 = dlog(l3)
      t31 = 1/l3
      tncoff3 = (-24*t1+24+12*t3-4*t5-12*t8+12*t1*t8+3*t12)*t15/4608-(-t
     #8+2-t18-2*t20)*t15*t25/384+(-t18-t8+2-t28)*t15*t31/768-t31*t15*t25
     #/384
      t1 = dlog(l1)
      t4 = dlog(1-l1)
      t5 = dlog(l3)
      t7 = 1/l1
      t9 = 1/l3
      tncoff4 = (2*t1-2+t4+t5)*t7*t9/768+t9*t7/(-1+z)/384
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = 0.3141592653589793D1**2
      t7 = dlog(z)
      t8 = t7**2
      t16 = -1/t1
      t19 = dlog(l3)
      t23 = 1/l3
      t27 = dlog(l1)
      t30 = 1/l1
      tncoff5 = -(6*t3-2*t5+12+6*t8+12*t7-12*t2-12*t2*t7)*t16/2304-(-t19
     #+2*t7+2-2*t2)*t23*t16/768-(2-2*t27+2*t7-2*t2)*t30*t16/384-t23*t30*
     #t16/384
      t1 = dlog(z)
      t3 = dlog(l3)
      t4 = 1-z
      t5 = dlog(t4)
      t8 = 1/l3
      t11 = -1/t4
      tncoff6 = -(-2*t1-2+t3+2*t5)*t8*t11/768+t8/l1*t11/384
      t2 = dlog(1-l1)
      t3 = dlog(z)
      t5 = dlog(l1)
      t7 = 1-z
      t8 = dlog(t7)
      t11 = 1/l1
      t14 = -1/t7
      tncoff7 = -(t2-2*t3+2*t5-2+2*t8)*t11*t14/384+1/l3*t11*t14/384
      tncoff8 = -1/l3/l1/(-1+z)/384
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l3)
      t14 = 1/l3
      t19 = dlog(l1)
      t22 = 1/l1
      tncoff9 = t1/2304-1.D0/384.D0-(-6+6*t4)*t8/2304-(12-6*t11)*t14/921
     #6+t14*t8/768+(-12+12*t19)*t22/4608+t22*t8/384-t14*t22/768
      t1 = dlog(l3)
      t4 = 1/l3
      tncoff10 = -(6*t1-12)*t4/9216-t4/(-1+z)/768+t4/l1/768
      t2 = dlog(1-l1)
      t4 = dlog(l1)
      t7 = 1/l1
      tncoff11 = (-6*t2-12*t4+12)*t7/4608-t7/(-1+z)/384+1/l3*t7/768
      tncoff12 = -1/l3/l1/768
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff13 = -(6+6*t1-6*t3)*t7/2304-1/l3*t7/768-1/l1*t7/384
      tncoff14 = 1/l3/(-1+z)/768
      tncoff15 = 1/l1/(-1+z)/384
      tncoff16 = -1.D0/768.D0+1/(-1+z)/768-1/l3/1536-1/l1/768
      tncoff17 = 1/l3/1536
      tncoff18 = 1/l1/768
      tncoff19 = -1/(-1+z)/768
      tncoff20 = -1.D0/1536.D0

CCC Numerator functions

      tcnum1 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum7e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum7e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum7e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,999
     #.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum7e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum7e1(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum7e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum7e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum7e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum7e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = RVnum7e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = RVnum7e2(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = RVnum7e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = RVnum7e3(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      RV10e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+
     #tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff1
     #6+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncof
     #f20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV10e2int(var,ndim,Tcoeff)

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
      double precision RVnum7e0
      external RVnum7e0
      double precision RVnum7e1
      external RVnum7e1
      double precision RVnum7e2
      external RVnum7e2
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
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l3)
      t13 = 1/l3
      t18 = dlog(l1)
      t21 = 1/l1
      tncoff1 = t1/2304-1.D0/384.D0-(-2+2*t4)*t8/768+(-2+t11)*t13/1536+t
     #13*t8/768+(-2+2*t18)*t21/768+t21*t8/384-t13*t21/768
      t1 = dlog(l3)
      t3 = 1/l3
      tncoff2 = (2-t1)*t3/1536-t3/(-1+z)/768+t3/l1/768
      t2 = dlog(1-l1)
      t3 = dlog(l1)
      t6 = 1/l1
      tncoff3 = (-t2+2-2*t3)*t6/768-t6/(-1+z)/384+1/l3*t6/768
      tncoff4 = -1/l3/l1/768
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff5 = -(2+2*t1-2*t3)*t7/768-1/l3*t7/768-1/l1*t7/384
      tncoff6 = 1/l3/(-1+z)/768
      tncoff7 = 1/l1/(-1+z)/384
      tncoff8 = -1.D0/768.D0+1/(-1+z)/768-1/l3/1536-1/l1/768
      tncoff9 = 1/l3/1536
      tncoff10 = 1/l1/768
      tncoff11 = -1/(-1+z)/768
      tncoff12 = -1.D0/1536.D0

CCC Numerator functions

      tcnum1 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum7e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum7e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum7e1(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum7e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum7e2(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV10e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV10e3int(var,ndim,Tcoeff)

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
      double precision RVnum7e0
      external RVnum7e0
      double precision RVnum7e1
      external RVnum7e1
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

      tncoff1 = -1.D0/768.D0+1/(-1+z)/768-1/l3/1536-1/l1/768
      tncoff2 = 1/l3/1536
      tncoff3 = 1/l1/768
      tncoff4 = -1/(-1+z)/768
      tncoff5 = -1.D0/1536.D0

CCC Numerator functions

      tcnum1 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum7e0(var,ndim,1.D0,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum7e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum7e1(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      RV10e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV10e4int(var,ndim,Tcoeff)

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
      double precision RVnum7e0
      external RVnum7e0
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

      tncoff1 = -1.D0/1536.D0

CCC Numerator functions

      tcnum1 = RVnum7e0(var,ndim,1.D0,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV10e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

