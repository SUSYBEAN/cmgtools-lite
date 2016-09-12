CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbc20e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum20e0
      external qqbcnum20e0
      double precision qqbcnum20e1
      external qqbcnum20e1
      double precision qqbcnum20e2
      external qqbcnum20e2
      double precision qqbcnum20e3
      external qqbcnum20e3
      double precision qqbcnum20e4
      external qqbcnum20e4
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t4 = 1/l1
      t6 = -1+l1
      t7 = 1/t6
      t10 = 0.3141592653589793D1**2
      t14 = dlog(l1)
      t18 = z-1
      t19 = 1/t18
      t21 = dlog(-t18)
      t22 = t21**2
      t29 = -t6
      t36 = t10**2
      t38 = t14**2
      tncoff1 = (4*z3-4*z3*l1)*t4*t7/2304-(t10-t10*l1)*t4*t7*t14/2304-t1
     #9*t22*t21/24+t21*t10*t19/576-t29*t7*t4*t19*t21*t14/48+t36/276480-t
     #4*t38*t14/1152-t22*t29*t4*t7*t19/24
      t1 = -1+l1
      t2 = 1/t1
      t4 = dlog(-t1)
      t5 = t4**2
      t9 = 0.3141592653589793D1**2
      t15 = 1-z
      t16 = dlog(t15)
      t17 = t16**2
      t20 = -1/t15
      t27 = t9**2
      tncoff2 = t2*t5*t4/1152-t9*t2*t4/2304+z3*t2/576-t17*t2*t20/24-t2*t
     #20*t16*t4/48+t27/276480-t20*t17*t16/24+t16*t9*t20/576
      t1 = 1/l1
      t2 = -1+l1
      t3 = 1/t2
      t4 = t1*t3
      t6 = dlog(-t2)
      t7 = dlog(l1)
      t8 = t7**2
      t12 = t6**2
      t18 = 0.3141592653589793D1**2
      t19 = t18*t1
      t30 = z-1
      t31 = 1/t30
      t33 = dlog(-t30)
      t34 = t31*t33
      t38 = t33**2
      tncoff3 = -t4*t6*t8/384-t4*t12*t6/1152+(-t4*t12/384+t19*t3/2304)*t
     #7+t19*t3*t6/2304-z3*t1*t3/576+t4*t34*t7/48+t38*t1*t3*t31/24+t4*t34
     #*t6/48-t4*t8*t7/1152
      t1 = 1-z
      t2 = dlog(t1)
      t4 = -1/t1
      t5 = t2*t4
      t6 = dlog(z)
      t7 = t6**2
      t10 = 1-l1
      t14 = 1/l1
      t19 = t2**2
      t27 = dlog(l1)
      t34 = 0.3141592653589793D1**2
      tncoff4 = t5*t7/96+t14*t5*t6/48-t19*t14*t4/24-t19*t4*t6/24-t14*t4*
     #t2*t27/48+t4*t19*t2/24-t2*t34*t4/576
      t1 = -1+l1
      t2 = 1/t1
      t3 = z-1
      t4 = 1/t3
      t5 = t2*t4
      t7 = dlog(-t3)
      t8 = dlog(z)
      t13 = t8**2
      t16 = t7**2
      t20 = 0.3141592653589793D1**2
      t28 = dlog(-t1)
      tncoff5 = -t5*t7*t8/48+t7*t4*t13/96+t4*t16*t7/24-t7*t20*t4/576+t16
     #*t2*t4/24+t5*t7*t28/48-t16*t4*t8/24
      t1 = -1+l1
      t2 = 1/t1
      t3 = 1/l1
      t4 = t2*t3
      t5 = z-1
      t6 = 1/t5
      t8 = dlog(-t5)
      t9 = t6*t8
      t10 = dlog(l1)
      t14 = dlog(z)
      t18 = t8**2
      t24 = dlog(-t1)
      tncoff6 = -t4*t9*t10/48+t4*t9*t14/48-t18*t3*t2*t6/24-t4*t9*t24/48
      t1 = 1/l1
      t2 = dlog(l1)
      t3 = t2**2
      t6 = 0.3141592653589793D1**2
      t15 = z-1
      t16 = 1/t15
      t19 = dlog(-t15)
      t22 = t19**2
      tncoff7 = t1*t3/384+(t6-t6*l1)*t1/(-1+l1)/2304-z3/576-t1*t16*t19/4
     #8+t22*t16/24
      t1 = -1+l1
      t2 = 1/t1
      t4 = dlog(-t1)
      t5 = t4**2
      t8 = 0.3141592653589793D1**2
      t11 = z-1
      t12 = 1/t11
      t15 = dlog(-t11)
      t19 = t15**2
      tncoff8 = -t2*t5/384+t8*t2/2304+t2*t12*t15/48-z3/576+t19*t12/24
      t1 = 1/l1
      t2 = -1+l1
      t3 = 1/t2
      t4 = t1*t3
      t6 = dlog(-t2)
      t7 = dlog(l1)
      t11 = t6**2
      t14 = t7**2
      t17 = 0.3141592653589793D1**2
      t21 = z-1
      t24 = dlog(-t21)
      tncoff9 = t4*t6*t7/192+t4*t11/384+t4*t14/384-t17*t1*t3/2304-t4/t21
     #*t24/48
      t1 = 1-z
      t2 = dlog(t1)
      t4 = -1/t1
      t6 = dlog(z)
      t9 = t2**2
      t12 = 1-l1
      tncoff10 = t2*t4*t6/48-t9*t4/24+1/l1*t4*t2/48
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = -1/t1
      t14 = dlog(z)
      tncoff11 = -t3*t5/24-1/(-1+l1)*t5*t2/48+t2*t5*t14/48
      t5 = z-1
      t8 = dlog(-t5)
      tncoff12 = 1/(-1+l1)/l1/t5*t8/48
      t1 = -1+l1
      t6 = dlog(l1)
      t10 = 0.3141592653589793D1**2
      t12 = 1-z
      t13 = dlog(t12)
      tncoff13 = -1/l1*t6/192-t10/2304+t13/t12/48
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t9 = -1+l1
      t12 = dlog(-t9)
      tncoff14 = -t1/2304+t4/t3/48+1/t9*t12/192
      t1 = -1+l1
      t4 = 1/t1/l1
      t5 = dlog(l1)
      t8 = dlog(-t1)
      tncoff15 = -t4*t5/192-t4*t8/192
      t1 = 1-z
      t2 = dlog(t1)
      tncoff16 = -t2/t1/48
      t1 = 1-z
      t2 = dlog(t1)
      tncoff17 = -t2/t1/48
      t1 = -1+l1
      tncoff18 = 1/l1/192
      tncoff19 = -1/(-1+l1)/192
      tncoff20 = 1/(-1+l1)/l1/192
      tncoff21 = 1.D0/192.D0
      tncoff22 = 1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum20e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum20e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum20e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum20e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum20e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum20e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum20e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum20e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum20e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum20e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum20e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum20e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum20e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum20e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum20e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum20e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum20e2(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum20e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum20e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbcnum20e3(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbcnum20e4(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbcnum20e4(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff22)

CCC Sum

      t12 = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tncoff4+
     #tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+tcnum9
     #*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11
      t24 = tcnum12*tncoff12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*t
     #ncoff15+tcnum16*tncoff16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19
     #*tncoff19+tcnum20*tncoff20+tcnum21*tncoff21+tcnum22*tncoff22
      qqbc20e0int = t12+t24

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc20e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum20e0
      external qqbcnum20e0
      double precision qqbcnum20e1
      external qqbcnum20e1
      double precision qqbcnum20e2
      external qqbcnum20e2
      double precision qqbcnum20e3
      external qqbcnum20e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t2 = 1-z
      t3 = dlog(t2)
      t4 = t3**2
      t6 = -1/t2
      t9 = 1/l1
      t10 = dlog(l1)
      t11 = t10**2
      t14 = 0.3141592653589793D1**2
      t18 = -1+l1
      t19 = 1/t18
      tncoff1 = -z3/576+t4*t6/24+t9*t11/384+(t14-t14*l1)*t9*t19/2304-t18
     #*t19*t9*t6*t3/48
      t1 = -1+l1
      t2 = 1/t1
      t4 = dlog(-t1)
      t5 = t4**2
      t8 = 0.3141592653589793D1**2
      t11 = z-1
      t12 = 1/t11
      t15 = dlog(-t11)
      t19 = t15**2
      tncoff2 = -t2*t5/384+t8*t2/2304+t2*t12*t15/48-z3/576+t19*t12/24
      t1 = 1/l1
      t2 = -1+l1
      t3 = 1/t2
      t4 = t1*t3
      t6 = dlog(-t2)
      t7 = dlog(l1)
      t11 = t6**2
      t14 = t7**2
      t17 = 0.3141592653589793D1**2
      t21 = z-1
      t24 = dlog(-t21)
      tncoff3 = t4*t6*t7/192+t4*t11/384+t4*t14/384-t17*t1*t3/2304-t4/t21
     #*t24/48
      t1 = 1-z
      t2 = dlog(t1)
      t4 = -1/t1
      t6 = dlog(z)
      t9 = t2**2
      tncoff4 = t2*t4*t6/48-t9*t4/24+1/l1*t4*t2/48
      t1 = 1-z
      t2 = dlog(t1)
      t3 = t2**2
      t5 = -1/t1
      t14 = dlog(z)
      tncoff5 = -t3*t5/24-1/(-1+l1)*t5*t2/48+t2*t5*t14/48
      t5 = z-1
      t8 = dlog(-t5)
      tncoff6 = 1/(-1+l1)/l1/t5*t8/48
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t10 = dlog(l1)
      tncoff7 = -t1/2304+t4/t3/48-1/l1*t10/192
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t9 = -1+l1
      t12 = dlog(-t9)
      tncoff8 = -t1/2304+t4/t3/48+1/t9*t12/192
      t1 = -1+l1
      t4 = 1/t1/l1
      t5 = dlog(l1)
      t8 = dlog(-t1)
      tncoff9 = -t4*t5/192-t4*t8/192
      t1 = 1-z
      t2 = dlog(t1)
      tncoff10 = -t2/t1/48
      t1 = 1-z
      t2 = dlog(t1)
      tncoff11 = -t2/t1/48
      t1 = -1+l1
      tncoff12 = 1/l1/192
      tncoff13 = -1/(-1+l1)/192
      tncoff14 = 1/(-1+l1)/l1/192
      tncoff15 = 1.D0/192.D0
      tncoff16 = 1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum20e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum20e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum20e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum20e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum20e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum20e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum20e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum20e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum20e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum20e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum20e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum20e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum20e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum20e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum20e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum20e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff16)

CCC Sum

      qqbc20e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc20e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum20e0
      external qqbcnum20e0
      double precision qqbcnum20e1
      external qqbcnum20e1
      double precision qqbcnum20e2
      external qqbcnum20e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t10 = dlog(l1)
      tncoff1 = -t1/2304+t4/t3/48-1/l1*t10/192
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t9 = -1+l1
      t12 = dlog(-t9)
      tncoff2 = -t1/2304+t4/t3/48+1/t9*t12/192
      t1 = -1+l1
      t4 = 1/t1/l1
      t5 = dlog(l1)
      t8 = dlog(-t1)
      tncoff3 = -t4*t5/192-t4*t8/192
      t1 = 1-z
      t2 = dlog(t1)
      tncoff4 = -t2/t1/48
      t1 = 1-z
      t2 = dlog(t1)
      tncoff5 = -t2/t1/48
      t1 = 1-l1
      tncoff6 = 1/l1/192
      tncoff7 = -1/(-1+l1)/192
      tncoff8 = 1/(-1+l1)/l1/192
      tncoff9 = 1.D0/192.D0
      tncoff10 = 1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum20e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum20e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum20e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum20e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum20e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum20e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum20e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum20e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum20e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum20e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff10)

CCC Sum

      qqbc20e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc20e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum20e0
      external qqbcnum20e0
      double precision qqbcnum20e1
      external qqbcnum20e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1-l1
      tncoff1 = 1/l1/192
      tncoff2 = -1/(-1+l1)/192
      tncoff3 = 1/(-1+l1)/l1/192
      tncoff4 = 1.D0/192.D0
      tncoff5 = 1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum20e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum20e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum20e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum20e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum20e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qqbc20e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc20e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qqbcnum20e0
      external qqbcnum20e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = 1.D0/192.D0
      tncoff2 = 1.D0/192.D0

CCC Numerator functions

      tcnum1 = qqbcnum20e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum20e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbc20e4int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

