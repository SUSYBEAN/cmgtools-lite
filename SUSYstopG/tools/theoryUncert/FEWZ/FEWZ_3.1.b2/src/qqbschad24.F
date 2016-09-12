CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbsc24e0int(var,ndim,Tcoeff)

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
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qqbscnum24e0
      external qqbscnum24e0
      double precision qqbscnum24e1
      external qqbscnum24e1
      double precision qqbscnum24e2
      external qqbscnum24e2
      double precision qqbscnum24e3
      external qqbscnum24e3
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
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t3 = 1-z
      t4 = dlog(t3)
      t5 = t4**2
      t8 = -1/t3
      t11 = 1/l1
      t12 = Lf*t11
      t13 = dlog(l1)
      t14 = t13**2
      t17 = 0.3141592653589793D1**2
      tncoff1 = Lf*z3/576-t5*Lf*t8/24-t12*t14/384+Lf*(-t17+t17*l1)*t11/(
     #-1+l1)/2304+t12*t8*t4/48
      t3 = 1-z
      t4 = dlog(t3)
      t5 = t4**2
      t8 = -1/t3
      t11 = -1+l1
      t12 = 1/t11
      t13 = Lf*t12
      t15 = dlog(-t11)
      t16 = t15**2
      t19 = 0.3141592653589793D1**2
      tncoff2 = Lf*z3/576-t5*Lf*t8/24+t13*t16/384-Lf*t19*t12/2304-t13*t8
     #*t4/48
      t1 = 1/l1
      t2 = Lf*t1
      t3 = -1+l1
      t4 = 1/t3
      t5 = dlog(l1)
      t6 = t5**2
      t11 = dlog(-t3)
      t16 = t11**2
      t20 = 0.3141592653589793D1**2
      t26 = z-1
      t30 = dlog(-t26)
      tncoff3 = -t2*t4*t6/384-t2*t4*t11*t5/192-t2*t4*t16/384+Lf*t20*t1*t
     #4/2304+Lf*t4*t1/t26*t30/48
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t6 = dlog(z)
      t10 = t2**2
      t14 = 1-l1
      tncoff4 = -t2*Lf*t5*t6/48+t10*Lf*t5/24-Lf/l1*t5*t2/48
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t6 = dlog(z)
      t10 = t2**2
      tncoff5 = -t2*Lf*t5*t6/48+t10*Lf*t5/24+Lf/(-1+l1)*t5*t2/48
      t5 = z-1
      t9 = dlog(-t5)
      tncoff6 = -Lf/(-1+l1)/l1/t5*t9/48
      t1 = 0.3141592653589793D1**2
      t4 = 1-z
      t5 = dlog(t4)
      t13 = dlog(l1)
      tncoff7 = Lf*t1/2304-t5*Lf/t4/48+Lf/l1*t13/192
      t1 = 0.3141592653589793D1**2
      t4 = 1-z
      t5 = dlog(t4)
      t11 = -1+l1
      t15 = dlog(-t11)
      tncoff8 = Lf*t1/2304-t5*Lf/t4/48-Lf/t11*t15/192
      t2 = Lf/l1
      t3 = -1+l1
      t4 = 1/t3
      t5 = dlog(l1)
      t9 = dlog(-t3)
      tncoff9 = t2*t4*t5/192+t2*t4*t9/192
      t1 = 1-z
      t2 = dlog(t1)
      tncoff10 = t2*Lf/t1/48
      t1 = 1-z
      t2 = dlog(t1)
      tncoff11 = t2*Lf/t1/48
      t1 = 1-l1
      tncoff12 = -Lf/l1/192
      tncoff13 = Lf/(-1+l1)/192
      tncoff14 = -Lf/l1/(-1+l1)/192
      tncoff15 = -Lf/192
      tncoff16 = -Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum24e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum24e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum24e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum24e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum24e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbscnum24e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbscnum24e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbscnum24e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbscnum24e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbscnum24e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbscnum24e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbscnum24e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbscnum24e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbscnum24e2(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbscnum24e3(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbscnum24e3(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff16)

CCC Sum

      qqbsc24e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff
     #12+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tnco
     #ff16

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc24e1int(var,ndim,Tcoeff)

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
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qqbscnum24e0
      external qqbscnum24e0
      double precision qqbscnum24e1
      external qqbscnum24e1
      double precision qqbscnum24e2
      external qqbscnum24e2
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
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 0.3141592653589793D1**2
      t4 = 1-z
      t5 = dlog(t4)
      t13 = dlog(l1)
      tncoff1 = Lf*t1/2304-t5*Lf/t4/48+Lf/l1*t13/192
      t1 = 0.3141592653589793D1**2
      t4 = 1-z
      t5 = dlog(t4)
      t11 = -1+l1
      t15 = dlog(-t11)
      tncoff2 = Lf*t1/2304-t5*Lf/t4/48-Lf/t11*t15/192
      t2 = Lf/l1
      t3 = -1+l1
      t4 = 1/t3
      t5 = dlog(l1)
      t9 = dlog(-t3)
      tncoff3 = t2*t4*t5/192+t2*t4*t9/192
      t1 = 1-z
      t2 = dlog(t1)
      tncoff4 = t2*Lf/t1/48
      t1 = 1-z
      t2 = dlog(t1)
      tncoff5 = t2*Lf/t1/48
      tncoff6 = -Lf/l1/192
      tncoff7 = Lf/(-1+l1)/192
      tncoff8 = -Lf/l1/(-1+l1)/192
      tncoff9 = -Lf/192
      tncoff10 = -Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum24e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum24e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum24e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum24e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum24e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbscnum24e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbscnum24e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbscnum24e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbscnum24e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbscnum24e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff10)

CCC Sum

      qqbsc24e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc24e2int(var,ndim,Tcoeff)

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
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qqbscnum24e0
      external qqbscnum24e0
      double precision qqbscnum24e1
      external qqbscnum24e1
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
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -Lf/l1/192
      tncoff2 = Lf/(-1+l1)/192
      tncoff3 = -Lf/l1/(-1+l1)/192
      tncoff4 = -Lf/192
      tncoff5 = -Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum24e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum24e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum24e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum24e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum24e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qqbsc24e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc24e3int(var,ndim,Tcoeff)

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
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qqbscnum24e0
      external qqbscnum24e0
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
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -Lf/192
      tncoff2 = -Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum24e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum24e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbsc24e3int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc24e4int(var,ndim,Tcoeff)

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
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
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
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      qqbsc24e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

