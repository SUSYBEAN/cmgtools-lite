CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbsc25e0int(var,ndim,Tcoeff)

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
      double precision qqbscnum25e0
      external qqbscnum25e0
      double precision qqbscnum25e1
      external qqbscnum25e1
      double precision qqbscnum25e2
      external qqbscnum25e2
      double precision qqbscnum25e3
      external qqbscnum25e3
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

      t3 = 0.3141592653589793D1**2
      t5 = z-1
      t6 = 1/t5
      t11 = dlog(-t5)
      t12 = t11**2
      t18 = 1/l2
      t20 = 1/(-1+l2)
      t28 = 1/(-1+l1)
      t29 = 1/l1
      t33 = Lf*t18
      t34 = dlog(l2)
      t35 = t34**2
      t41 = t6*t11
      t44 = Lf*t29
      t45 = dlog(l1)
      t46 = t45**2
      t57 = Lf*(-l2-l1+1+l1*l2)*t29
      t58 = t28*t20
      tncoff1 = -Lf*z3/576+Lf*t3*t6/1152-Lf*t6*t12/48-Lf*(-t3+t3*l2)*t18
     #*t20/1152-Lf*(-t3+t3*l1)*t28*t29/2304+t33*t35/48+t33*t6*t34/24+t33
     #*t41/24+t44*t46/384+t44*t6*t45/96+t44*t41/48-t57*t58*t18*t6/48-t57
     #*t58*t18*t34/48-t57*t58*t18*t45/96
      t2 = 1/(-1+l2)
      t4 = z-1
      t5 = 1/t4
      t7 = dlog(-t4)
      t15 = Lf*(l2-l1*l2)/l1
      t18 = 1/(-1+l1)*t2
      t19 = 1/l2
      t24 = dlog(l1)
      t29 = 0.3141592653589793D1**2
      tncoff2 = -Lf*t2*t5*t7/24-t15*t18*t19*t5/48-t15*t18*t19*t24/96+Lf*
     #t29*t2/1152
      t2 = Lf/l1
      t4 = 1/(-1+l2)
      t5 = 1/l2
      t6 = t4*t5
      t7 = z-1
      t8 = 1/t7
      t12 = dlog(l2)
      t16 = dlog(l1)
      t20 = 0.3141592653589793D1**2
      t24 = Lf*t5
      t25 = t12**2
      t29 = t4*t8
      t34 = dlog(-t7)
      tncoff3 = -t2*t6*t8/48-t2*t6*t12/48-t2*t6*t16/96-Lf*t20*t6/1152+t2
     #4*t4*t25/48+t24*t29*t12/24+t24*t29*t34/24
      t1 = 1/l2
      t2 = Lf*t1
      t3 = z-1
      t4 = 1/t3
      t6 = dlog(-t3)
      t7 = t4*t6
      t10 = dlog(l2)
      t18 = Lf*(l1-l1*l2)/l1
      t19 = -1+l1
      t20 = 1/t19
      t22 = 1/(-1+l2)
      t23 = t20*t22
      t29 = dlog(-t19)
      t35 = t6**2
      t38 = 0.3141592653589793D1**2
      t39 = Lf*t38
      t52 = Lf*t20
      t53 = t29**2
      t59 = t10**2
      tncoff4 = t2*t7/24+t2*t4*t10/24-t18*t23*t1*t10/48-t18*t23*t1*t29/9
     #6-Lf*t4*t35/48+t39*t20/2304-Lf*z3/576+t39*t4/1152+Lf*(t38-t38*l2)*
     #t1*t22/1152-t52*t53/384-t52*t4*t29/96+t2*t59/48-t52*t7/48-t18*t23*
     #t1*t4/48
      t2 = 1/(-1+l2)
      t4 = z-1
      t5 = 1/t4
      t7 = dlog(-t4)
      t11 = -1+l1
      t13 = Lf/t11
      t15 = dlog(-t11)
      t19 = 0.3141592653589793D1**2
      tncoff5 = -Lf*t2*t5*t7/24-t13*t2*t15/96+Lf*t19*t2/1152-t13*t2*t5/4
     #8
      t1 = -1+l1
      t3 = Lf/t1
      t5 = 1/(-1+l2)
      t6 = 1/l2
      t7 = t5*t6
      t9 = dlog(-t1)
      t13 = 0.3141592653589793D1**2
      t17 = Lf*t6
      t18 = dlog(l2)
      t19 = t18**2
      t23 = z-1
      t24 = 1/t23
      t28 = t5*t24
      t30 = dlog(-t23)
      tncoff6 = t3*t7*t9/96-Lf*t13*t7/1152+t17*t5*t19/48+t3*t7*t24/48+t1
     #7*t28*t30/24+t17*t28*t18/24+t3*t7*t18/48
      t1 = 1/l1
      t2 = Lf*t1
      t3 = -1+l1
      t4 = 1/t3
      t6 = t4/l2
      t7 = dlog(l2)
      t11 = dlog(l1)
      t16 = dlog(-t3)
      t20 = 0.3141592653589793D1**2
      t29 = Lf*t4
      t30 = t16**2
      t34 = t11**2
      t38 = z-1
      t39 = 1/t38
      t40 = t1*t39
      t48 = dlog(-t38)
      tncoff7 = -t2*t6*t7/48-t2*t6*t11/96-t2*t6*t16/96-Lf*t20*t4*t1/2304
     #+t2*t4*t16*t11/192+t29*t1*t30/384+t29*t1*t34/384+t29*t40*t16/96+t2
     #9*t40*t11/96+t29*t40*t48/48-t2*t6*t39/48
      t2 = Lf/l1
      t3 = -1+l1
      t7 = 1/t3/(-1+l2)
      t9 = dlog(-t3)
      t18 = dlog(l1)
      tncoff8 = t2*t7*t9/96+t2*t7/(z-1)/48+t2*t7*t18/96
      t3 = -1+l1
      t5 = Lf/l1/t3
      t9 = 1/(-1+l2)/l2
      t15 = dlog(l2)
      t19 = dlog(l1)
      t24 = dlog(-t3)
      tncoff9 = -t5*t9/(z-1)/48-t5*t9*t15/48-t5*t9*t19/96-t5*t9*t24/96
      t1 = z*Lf
      t2 = 0.3141592653589793D1**2
      t3 = z-1
      t4 = 1/t3
      t9 = dlog(-t3)
      t12 = dlog(z)
      t16 = t12**2
      t20 = t9**2
      t24 = l2*z
      t25 = -t24+z
      t27 = 1/l2
      t28 = Lf*t25*t27
      t30 = 1/(-1+l2)
      t31 = t30*t4
      t32 = dlog(l2)
      t46 = l1*z
      t49 = -1+l1
      t50 = 1/t49
      t51 = Lf*(z-t46)*t50
      t52 = 1/l1
      t53 = t52*t4
      t54 = dlog(l1)
      tncoff10 = -t1*t2*t4/1152-t9*Lf*z*t4*t12/24+t1*t4*t16/48+t1*t4*t20
     #/48+t28*t31*t32/24-Lf*t25*t27*t31*t12/24+t28*t31*t9/24+t51*t53*t54
     #/96+t1*t49*t50*t52*t4*t12/48+t51*t53*t9/48-Lf*(-t46*l2-z+t24+t46)*
     #t52*t50*t30*t27*t4/48
      t2 = 1/(-1+l2)
      t3 = Lf*t2
      t4 = z-1
      t5 = 1/t4
      t6 = dlog(z)
      t11 = dlog(-t4)
      tncoff11 = -t3*t5*t6/48+t3*t5*t11/24-Lf*(-l2+l1*l2)/l1/(-1+l1)*t2/
     #l2*t5/48
      t2 = l1*z
      t3 = l2*z
      t12 = 1/(-1+l2)
      t14 = 1/l2
      t15 = z-1
      t16 = 1/t15
      t21 = -l2+t3-z
      t23 = Lf*t21*t14
      t24 = t12*t16
      t25 = dlog(l2)
      t31 = dlog(-z/t21)
      t36 = dlog(-t15)
      t41 = dlog(-t21)
      tncoff12 = -Lf*(-l1*l2-t2-t3+z+l2+t2*l2)/l1/(-1+l1)*t12*t14*t16/48
     #+t23*t24*t25/24-t23*t24*t31/48+t23*t24*t36/24-t23*t24*t41/24
      t1 = z*Lf
      t2 = z-1
      t3 = 1/t2
      t5 = dlog(-t2)
      t6 = t5**2
      t10 = 0.3141592653589793D1**2
      t16 = dlog(z)
      t20 = t16**2
      t24 = -1+l1
      t25 = 1/t24
      t26 = t25*t3
      t28 = dlog(-t24)
      t38 = l1*z
      t45 = 1/(-1+l2)
      t47 = 1/l2
      t53 = -l2*z+z
      t55 = Lf*t53*t47
      t56 = t45*t3
      t60 = dlog(l2)
      tncoff13 = t1*t3*t6/48-t1*t10*t3/1152-t5*Lf*z*t3*t16/24+t1*t3*t20/
     #48+t1*t26*t28/96-t1*t26*t16/48+t1*t26*t5/48-Lf*(-t38+t38*l2)/l1*t2
     #5*t45*t47*t3/48+t55*t56*t5/24+t55*t56*t60/24-Lf*t53*t47*t56*t16/24
      t5 = 1/(-1+l2)
      t6 = z-1
      t7 = 1/t6
      t11 = Lf*t5
      t13 = dlog(-t6)
      t17 = dlog(z)
      tncoff14 = Lf/(-1+l1)*t5*t7/48+t11*t7*t13/24-t11*t7*t17/48
      t2 = -l2+l2*z-z
      t4 = 1/l2
      t5 = Lf*t2*t4
      t7 = 1/(-1+l2)
      t8 = z-1
      t9 = 1/t8
      t10 = t7*t9
      t11 = dlog(l2)
      t17 = dlog(-z/t2)
      t22 = dlog(-t2)
      t26 = l1*z
      t41 = dlog(-t8)
      tncoff15 = t5*t10*t11/24-t5*t10*t17/48-t5*t10*t22/24-Lf*(-t26*l2+l
     #1*l2+t26)/l1/(-1+l1)*t7*t4*t9/48+t5*t10*t41/24
      t2 = -1+l1
      t3 = 1/t2
      t4 = Lf*z*t3
      t5 = 1/l1
      t6 = z-1
      t7 = 1/t6
      t8 = t5*t7
      t10 = dlog(-t2)
      t14 = dlog(l1)
      t18 = dlog(z)
      t23 = dlog(-t6)
      tncoff16 = -t4*t8*t10/96-t4*t8*t14/96+t4*t8*t18/48-t4*t8*t23/48-Lf
     #*(-l2*z+z)*t5*t3/(-1+l2)/l2*t7/48
      tncoff17 = -Lf/l1/(-1+l1)/(-1+l2)/(z-1)/48
      tncoff18 = -Lf*(-l2+l2*z-z)/l1/(-1+l1)/(-1+l2)/l2/(z-1)/48
      t4 = 1/l1
      t11 = 1/l2
      t15 = 0.3141592653589793D1**2
      t18 = 1-z
      t19 = dlog(t18)
      t22 = -1/t18
      t25 = Lf*t11
      t26 = dlog(l2)
      t31 = Lf*t4
      t32 = dlog(l1)
      tncoff19 = Lf*(-l2-l1+1+l1*l2)*t4/(-1+l1)/(-1+l2)*t11/96-Lf*t15/23
     #04+t19*Lf*t22/48-t25*t26/48-t25*t22/48-t31*t32/192-t31*t22/96
      t9 = 1/(-1+l2)
      tncoff20 = Lf*(l2-l1*l2)/l1/(-1+l1)*t9/l2/96+Lf*t9/(z-1)/48
      t4 = 1/(-1+l2)
      t5 = 1/l2
      t9 = Lf*t5
      t10 = dlog(l2)
      tncoff21 = Lf/l1*t4*t5/96-t9*t4*t10/48-t9*t4/(z-1)/48
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t13 = -1+l1
      t14 = 1/t13
      t18 = 1/l2
      t22 = 0.3141592653589793D1**2
      t25 = Lf*t14
      t27 = dlog(-t13)
      t32 = Lf*t18
      t33 = dlog(l2)
      tncoff22 = t2*Lf*t5/48+Lf*(l1-l1*l2)/l1*t14/(-1+l2)*t18/96-Lf*t22/
     #2304+t25*t27/192+t25*t5/96-t32*t33/48-t32*t5/48
      t2 = 1/(-1+l2)
      tncoff23 = Lf*t2/(z-1)/48+Lf/(-1+l1)*t2/96
      t1 = 1/l2
      t2 = Lf*t1
      t4 = 1/(-1+l2)
      t5 = dlog(l2)
      tncoff24 = -t2*t4*t5/48-t2*t4/(z-1)/48-Lf/(-1+l1)*t4*t1/96
      t1 = 1/l1
      t2 = Lf*t1
      t3 = -1+l1
      t4 = 1/t3
      t9 = dlog(l1)
      t14 = dlog(-t3)
      tncoff25 = t2*t4/l2/96-t2*t4*t9/192-t2*t4*t14/192-Lf*t4*t1/(z-1)/9
     #6
      tncoff26 = -Lf/l1/(-1+l1)/(-1+l2)/96
      tncoff27 = Lf/l1/(-1+l1)/(-1+l2)/l2/96
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t10 = dlog(z)
      tncoff28 = -t2*Lf*z*t5/48+Lf*z*t5*t10/48-Lf*(-l2*z+z)/l2/(-1+l2)*t
     #5/48-Lf*(z-l1*z)/(-1+l1)/l1*t5/96
      tncoff29 = -Lf/(-1+l2)/(z-1)/48
      tncoff30 = -Lf*(-l2+l2*z-z)/l2/(-1+l2)/(z-1)/48
      t1 = Lf*z
      t4 = z-1
      t5 = 1/t4
      t20 = dlog(-t4)
      t25 = dlog(z)
      tncoff31 = -t1/(-1+l1)*t5/96-Lf*(-l2*z+z)/l2/(-1+l2)*t5/48-t20*Lf*
     #z*t5/48+t1*t5*t25/48
      tncoff32 = -Lf/(-1+l2)/(z-1)/48
      tncoff33 = -Lf*(-l2+l2*z-z)/l2/(-1+l2)/(z-1)/48
      tncoff34 = Lf*z/(-1+l1)/l1/(z-1)/96
      t5 = 1-l2
      t14 = 1-l1
      tncoff35 = -Lf/(z-1)/96+Lf/l2/96+Lf/l1/192
      tncoff36 = -Lf/(-1+l2)/96
      tncoff37 = Lf/l2/(-1+l2)/96
      t9 = -1+l2
      tncoff38 = -Lf/(z-1)/96-Lf/(-1+l1)/192+Lf/l2/96
      tncoff39 = -Lf/(-1+l2)/96
      tncoff40 = Lf/l2/(-1+l2)/96
      tncoff41 = Lf/l1/(-1+l1)/192
      tncoff42 = Lf*z/(z-1)/96
      tncoff43 = Lf*z/(z-1)/96
      tncoff44 = Lf/192
      tncoff45 = Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbscnum25e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbscnum25e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbscnum25e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbscnum25e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbscnum25e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbscnum25e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbscnum25e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbscnum25e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbscnum25e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbscnum25e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbscnum25e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbscnum25e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = qqbscnum25e1(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = qqbscnum25e1(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = qqbscnum25e1(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff27)
      tcnum28 = qqbscnum25e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff28)
      tcnum29 = qqbscnum25e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff29)
      tcnum30 = qqbscnum25e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff30)
      tcnum31 = qqbscnum25e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff31)
      tcnum32 = qqbscnum25e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff32)
      tcnum33 = qqbscnum25e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff33)
      tcnum34 = qqbscnum25e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff34)
      tcnum35 = qqbscnum25e2(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff35)
      tcnum36 = qqbscnum25e2(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff36)
      tcnum37 = qqbscnum25e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff37)
      tcnum38 = qqbscnum25e2(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff38)
      tcnum39 = qqbscnum25e2(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff39)
      tcnum40 = qqbscnum25e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff40)
      tcnum41 = qqbscnum25e2(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff41)
      tcnum42 = qqbscnum25e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff42)
      tcnum43 = qqbscnum25e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff43)
      tcnum44 = qqbscnum25e3(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff44)
      tcnum45 = qqbscnum25e3(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff45)

CCC Sum

      t12 = tcnum21*tncoff21+tcnum44*tncoff44+tcnum24*tncoff24+tcnum22*t
     #ncoff22+tcnum42*tncoff42+tcnum13*tncoff13+tcnum43*tncoff43+tcnum14
     #*tncoff14+tcnum38*tncoff38+tcnum11*tncoff11+tcnum12*tncoff12
      t24 = tcnum26*tncoff26+tcnum7*tncoff7+tcnum39*tncoff39+tcnum15*tnc
     #off15+tcnum8*tncoff8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum16*tncof
     #f16+tcnum25*tncoff25+tcnum29*tncoff29+tcnum28*tncoff28
      t37 = tcnum36*tncoff36+tcnum27*tncoff27+tcnum31*tncoff31+tcnum30*t
     #ncoff30+tcnum37*tncoff37+tcnum34*tncoff34+tcnum1*tncoff1+tcnum35*t
     #ncoff35+tcnum6*tncoff6+tcnum32*tncoff32+tcnum19*tncoff19
      t50 = tcnum33*tncoff33+tcnum20*tncoff20+tcnum45*tncoff45+tcnum17*t
     #ncoff17+tcnum23*tncoff23+tcnum2*tncoff2+tcnum3*tncoff3+tcnum18*tnc
     #off18+tcnum4*tncoff4+tcnum5*tncoff5+tcnum41*tncoff41+tcnum40*tncof
     #f40
      qqbsc25e0int = t12+t24+t37+t50

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc25e1int(var,ndim,Tcoeff)

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
      double precision qqbscnum25e0
      external qqbscnum25e0
      double precision qqbscnum25e1
      external qqbscnum25e1
      double precision qqbscnum25e2
      external qqbscnum25e2
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

      t4 = 1/l1
      t11 = 1/l2
      t15 = 0.3141592653589793D1**2
      t18 = 1-z
      t19 = dlog(t18)
      t22 = -1/t18
      t25 = Lf*t11
      t26 = dlog(l2)
      t31 = Lf*t4
      t32 = dlog(l1)
      tncoff1 = Lf*(-l2-l1+1+l1*l2)*t4/(-1+l1)/(-1+l2)*t11/96-Lf*t15/230
     #4+t19*Lf*t22/48-t25*t26/48-t25*t22/48-t31*t32/192-t31*t22/96
      t9 = 1/(-1+l2)
      tncoff2 = Lf*(l2-l1*l2)/l1/(-1+l1)*t9/l2/96+Lf*t9/(z-1)/48
      t4 = 1/(-1+l2)
      t5 = 1/l2
      t9 = Lf*t5
      t10 = dlog(l2)
      tncoff3 = Lf/l1*t4*t5/96-t9*t4*t10/48-t9*t4/(z-1)/48
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t13 = -1+l1
      t14 = 1/t13
      t18 = 1/l2
      t22 = 0.3141592653589793D1**2
      t25 = Lf*t14
      t27 = dlog(-t13)
      t32 = Lf*t18
      t33 = dlog(l2)
      tncoff4 = t2*Lf*t5/48+Lf*(l1-l1*l2)/l1*t14/(-1+l2)*t18/96-Lf*t22/2
     #304+t25*t27/192+t25*t5/96-t32*t33/48-t32*t5/48
      t2 = 1/(-1+l2)
      tncoff5 = Lf*t2/(z-1)/48+Lf/(-1+l1)*t2/96
      t1 = 1/l2
      t2 = Lf*t1
      t4 = 1/(-1+l2)
      t5 = dlog(l2)
      tncoff6 = -t2*t4*t5/48-t2*t4/(z-1)/48-Lf/(-1+l1)*t4*t1/96
      t1 = 1/l1
      t2 = Lf*t1
      t3 = -1+l1
      t4 = 1/t3
      t9 = dlog(l1)
      t14 = dlog(-t3)
      tncoff7 = t2*t4/l2/96-t2*t4*t9/192-t2*t4*t14/192-Lf*t4*t1/(z-1)/96
      tncoff8 = -Lf/l1/(-1+l1)/(-1+l2)/96
      tncoff9 = Lf/l1/(-1+l1)/(-1+l2)/l2/96
      t1 = 1-z
      t2 = dlog(t1)
      t5 = -1/t1
      t10 = dlog(z)
      tncoff10 = -t2*Lf*z*t5/48+Lf*z*t5*t10/48-Lf*(-l2*z+z)/l2/(-1+l2)*t
     #5/48-Lf*(z-l1*z)/(-1+l1)/l1*t5/96
      tncoff11 = -Lf/(-1+l2)/(z-1)/48
      tncoff12 = -Lf*(-l2+l2*z-z)/l2/(-1+l2)/(z-1)/48
      t1 = Lf*z
      t4 = z-1
      t5 = 1/t4
      t20 = dlog(-t4)
      t25 = dlog(z)
      tncoff13 = -t1/(-1+l1)*t5/96-Lf*(-l2*z+z)/l2/(-1+l2)*t5/48-t20*Lf*
     #z*t5/48+t1*t5*t25/48
      tncoff14 = -Lf/(-1+l2)/(z-1)/48
      tncoff15 = -Lf*(-l2+l2*z-z)/l2/(-1+l2)/(z-1)/48
      tncoff16 = Lf*z/(-1+l1)/l1/(z-1)/96
      tncoff17 = -Lf/(z-1)/96+Lf/l2/96+Lf/l1/192
      tncoff18 = -Lf/(-1+l2)/96
      tncoff19 = Lf/l2/(-1+l2)/96
      tncoff20 = -Lf/(z-1)/96-Lf/(-1+l1)/192+Lf/l2/96
      tncoff21 = -Lf/(-1+l2)/96
      tncoff22 = Lf/l2/(-1+l2)/96
      tncoff23 = Lf/l1/(-1+l1)/192
      tncoff24 = Lf*z/(z-1)/96
      tncoff25 = Lf*z/(z-1)/96
      tncoff26 = Lf/192
      tncoff27 = Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbscnum25e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbscnum25e0(var,ndim,1.D0,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbscnum25e0(var,ndim,1.D0,ux,l1,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbscnum25e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbscnum25e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbscnum25e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbscnum25e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbscnum25e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbscnum25e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbscnum25e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff20)
      tcnum21 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff21)
      tcnum22 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff22)
      tcnum23 = qqbscnum25e1(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff23)
      tcnum24 = qqbscnum25e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff24)
      tcnum25 = qqbscnum25e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff25)
      tcnum26 = qqbscnum25e2(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff26)
      tcnum27 = qqbscnum25e2(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff27)

CCC Sum

      t14 = tcnum10*tncoff10+tcnum13*tncoff13+tcnum20*tncoff20+tcnum21*t
     #ncoff21+tcnum15*tncoff15+tcnum22*tncoff22+tcnum26*tncoff26+tcnum27
     #*tncoff27+tcnum11*tncoff11+tcnum4*tncoff4+tcnum6*tncoff6+tcnum5*tn
     #coff5+tcnum1*tncoff1
      t29 = tcnum2*tncoff2+tcnum3*tncoff3+tcnum7*tncoff7+tcnum8*tncoff8+
     #tcnum9*tncoff9+tcnum19*tncoff19+tcnum23*tncoff23+tcnum25*tncoff25+
     #tcnum24*tncoff24+tcnum17*tncoff17+tcnum18*tncoff18+tcnum12*tncoff1
     #2+tcnum16*tncoff16+tcnum14*tncoff14
      qqbsc25e1int = t14+t29

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc25e2int(var,ndim,Tcoeff)

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
      double precision qqbscnum25e0
      external qqbscnum25e0
      double precision qqbscnum25e1
      external qqbscnum25e1
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

      tncoff1 = -Lf/(z-1)/96+Lf/l2/96+Lf/l1/192
      tncoff2 = -Lf/(-1+l2)/96
      tncoff3 = Lf/l2/(-1+l2)/96
      t9 = 1-l2
      tncoff4 = -Lf/(z-1)/96-Lf/(-1+l1)/192+Lf/l2/96
      tncoff5 = -Lf/(-1+l2)/96
      tncoff6 = Lf/l2/(-1+l2)/96
      tncoff7 = Lf/l1/(-1+l1)/192
      tncoff8 = Lf*z/(z-1)/96
      tncoff9 = Lf*z/(z-1)/96
      tncoff10 = Lf/192
      tncoff11 = Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbscnum25e0(var,ndim,1.D0,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbscnum25e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbscnum25e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbscnum25e1(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbscnum25e1(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l
     #8,999.D0,999.D0,-999.D0,Tcoeff*tncoff11)

CCC Sum

      qqbsc25e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4
     #*tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncof
     #f8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc25e3int(var,ndim,Tcoeff)

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
      double precision qqbscnum25e0
      external qqbscnum25e0
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

      tncoff1 = Lf/192
      tncoff2 = Lf/192

CCC Numerator functions

      tcnum1 = qqbscnum25e0(var,ndim,1.D0,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbscnum25e0(var,ndim,1.D0,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbsc25e3int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbsc25e4int(var,ndim,Tcoeff)

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

      qqbsc25e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

