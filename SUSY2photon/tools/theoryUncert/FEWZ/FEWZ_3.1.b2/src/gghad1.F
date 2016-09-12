CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm1e0int(var,ndim,Tcoeff)

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
      double precision ggrmnum1e0
      external ggrmnum1e0
      double precision ggrmnum1e1
      external ggrmnum1e1
      double precision ggrmnum1e2
      external ggrmnum1e2
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t5 = dlog(1-l4)
      t7 = dlog(1-l1)
      t10 = dlog(l1)
      t13 = dlog(l4)
      t16 = dlog(2.D0)
      t19 = dlog(z)
      t22 = dlog(t1)
      t26 = dlog(-t2)
      t35 = 18+12*t5*t7+24*t10*t7+6*t5*t13+48*t10*t16-24*t7*t19-24*t22*t
     #19+48*t10*t26+48*t22*t16+24*t10*t22+24*t5*t26
      t43 = t16**2
      t45 = 0.3141592653589793D1**2
      t53 = -12*t13*t19-12*t5*t19+48*t7*t26-48*t16+48*t43-t45-12*t5-12*t
     #13+24*t19+12*t13*t22-24*t10-24*t22
      t56 = t10**2
      t58 = t7**2
      t61 = t5**2
      t63 = t13**2
      t65 = t22**2
      t67 = t19**2
      t69 = t26**2
      t75 = -24*t7+12*t56+12*t58-48*t26+3*t61+3*t63+12*t65+12*t67+48*t69
     #+24*t5*t16+24*t13*t26
      t100 = 12*t13*t10+24*t7*t22+48*t22*t26+12*t5*t22-48*t26*t19-24*t10
     #*t19+48*t7*t16-48*t19*t16+96*t26*t16+12*t13*t7+24*t13*t16+12*t5*t1
     #0
      t111 = dlog(l8)
      t128 = 2-t13+2*t19-t5-2*t10-2*t7-2*t22-4*t26-4*t16-t111-2*l8+4*t26
     #*l8+t13*l8+4*t16*l8-2*t19*l8+2*t22*l8+2*t7*l8+2*t10*l8+t111*l8+t5*
     #l8
      t130 = 1/(-1+l8)
      t132 = 1/l8
      t136 = t1*t13
      t138 = 2*t1*t10
      t140 = 2*t1*t7
      t142 = 2*t1*t22
      t144 = 4*t1*t26
      t145 = 2*t1
      t147 = 4*t1*t16
      t149 = 2*t1*t19
      t150 = dlog(l3)
      t152 = t1*t5
      t155 = 1/l3
      t160 = t2*(-t1+t1*l8)
      t165 = dlog(l2)
      t169 = 1/l2
      tncoff1 = -3.D0/65536.D0*t3*(t35+t53+t75+t100)+9.D0/32768.D0*t3*t1
     #28*t130*t132+9.D0/32768.D0*t2*(t136+t138+t140+t142+t144-t145+t147-
     #t149+t150*t1+t152)*t155-9.D0/32768.D0*t160*t155*t132*t130+9.D0/327
     #68.D0*t2*(t136+t152+t165*t1+t140+t142+t144-t145+t147-t149+t138)*t1
     #69-9.D0/32768.D0*t160*t169*t132*t130-9.D0/32768.D0*t3*t169*t155
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t4 = dlog(2.D0)
      t6 = 0.3141592653589793D1**2
      t7 = t4**2
      t9 = dlog(l4)
      t12 = dlog(1-l4)
      t14 = dlog(z)
      t16 = dlog(l1)
      t19 = dlog(1-l1)
      t21 = dlog(t1)
      t23 = t16**2
      t25 = 18-48*t4-t6+48*t7-12*t9-12*t12+24*t14-24*t16-24*t19-24*t21+1
     #2*t23
      t26 = t19**2
      t29 = dlog(-t2)
      t31 = t12**2
      t33 = t9**2
      t35 = t14**2
      t37 = t21**2
      t39 = t29**2
      t51 = 12*t26-48*t29+3*t31+3*t33+12*t35+12*t37+48*t39+24*t19*t21+96
     #*t29*t4-48*t29*t14+12*t12*t21-12*t9*t14
      t75 = -24*t21*t14-12*t12*t14-48*t14*t4+12*t9*t16-24*t16*t14+48*t19
     #*t4+24*t12*t29+48*t19*t29+48*t16*t29+24*t9*t4+12*t12*t16
      t100 = 12*t12*t19+12*t9*t21+24*t12*t4+24*t9*t29+48*t21*t29+6*t12*t
     #9+24*t19*t16+48*t16*t4-24*t19*t14+24*t16*t21+48*t21*t4+12*t9*t19
      t117 = 1-l8
      t118 = dlog(t117)
      t125 = -1/t117
      t131 = t1*t9
      t133 = 2*t1*t16
      t135 = 2*t1*t19
      t137 = 2*t1*t21
      t139 = 4*t1*t29
      t140 = 2*t1
      t142 = 4*t1*t4
      t144 = 2*t1*t14
      t145 = dlog(l3)
      t147 = t1*t12
      t150 = 1/l3
      t156 = dlog(l2)
      t160 = 1/l2
      tncoff2 = -3.D0/65536.D0*t3*(t25+t51+t75+t100)+9.D0/32768.D0*t3*(-
     #4*t29*l8+2*t14*l8-4*t4*l8-2*t21*l8-2*t19*l8-2*t16*l8-t118*l8-t9*l8
     #+2*l8-t12*l8)*t125/l8+9.D0/32768.D0*t2*(t131+t133+t135+t137+t139-t
     #140+t142-t144+t145*t1+t147)*t150+9.D0/32768.D0*t3*t150*t125+9.D0/3
     #2768.D0*t2*(t131+t147+t156*t1+t135+t137+t139-t140+t142-t144+t133)*
     #t160+9.D0/32768.D0*t3*t160*t125-9.D0/32768.D0*t3*t160*t150
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t4 = dlog(t1)
      t6 = dlog(z)
      t8 = dlog(l1)
      t11 = dlog(1-l1)
      t13 = 1-l8
      t14 = dlog(t13)
      t15 = dlog(l4)
      t16 = dlog(l8)
      t18 = dlog(1-l4)
      t19 = dlog(2.D0)
      t22 = dlog(-t2)
      t24 = -2+2*t4-2*t6+2*t8+2*t11+t14+t15+t16+t18+4*t19+4*t22
      t26 = -1/t13
      t28 = 1/l8
      tncoff3 = 9.D0/32768.D0*t3*t24*t26*t28-9.D0/32768.D0*t3/l3*t28*t26
     #-9.D0/32768.D0*t3/l2*t28*t26
      t1 = -1+z
      t2 = dlog(l4)
      t3 = k1(z,l1,1.D0,l3,l4)
      t6 = dlog(1-l1)
      t10 = dlog(-t1)
      t14 = dlog(1-l3)
      t17 = dlog(t3)
      t20 = dlog(z)
      t23 = dlog(l1)
      t26 = dlog(2.D0)
      t29 = dlog(l3)
      t32 = dlog(1-l4)
      t34 = -t2*t3-2*t3*t6-4*t3*t10-t3*t14+2*t3-2*t3*t17+2*t3*t20-2*t3*t
     #23-4*t3*t26-t29*t3-t32*t3
      t36 = 1/l3
      tncoff4 = 9.D0/32768.D0*t1*t34*t36-9.D0/32768.D0*t1*(t3-t3*l8)*t36
     #/l8/(-1+l8)+9.D0/32768.D0*t1*t3/l2*t36
      t1 = -1+z
      t2 = dlog(l4)
      t3 = k1(z,l1,1.D0,l3,l4)
      t6 = dlog(1-l1)
      t10 = dlog(-t1)
      t14 = dlog(1-l3)
      t17 = dlog(t3)
      t20 = dlog(z)
      t23 = dlog(l1)
      t26 = dlog(2.D0)
      t29 = dlog(l3)
      t32 = dlog(1-l4)
      t34 = -t2*t3-2*t3*t6-4*t3*t10-t3*t14+2*t3-2*t3*t17+2*t3*t20-2*t3*t
     #23-4*t3*t26-t29*t3-t32*t3
      t36 = 1/l3
      t38 = t1*t3
      tncoff5 = 9.D0/32768.D0*t1*t34*t36-9.D0/32768.D0*t38*t36/(-1+l8)+9
     #.D0/32768.D0*t38/l2*t36
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff6 = 9.D0/32768.D0*(-1+z)*t2/l3/l8/(-1+l8)
      t1 = -1+z
      t2 = dlog(l4)
      t3 = 1-l2
      t4 = k1(z,l1,t3,0.D0,l4)
      t7 = dlog(1-l4)
      t9 = dlog(l2)
      t12 = dlog(1-l1)
      t15 = dlog(t4)
      t19 = dlog(-t1)
      t23 = dlog(t3)
      t25 = dlog(z)
      t28 = dlog(l1)
      t31 = dlog(2.D0)
      t34 = -t2*t4-t7*t4-t9*t4-2*t12*t4-2*t4*t15-4*t19*t4+2*t4-t4*t23+2*
     #t4*t25-2*t4*t28-4*t4*t31
      t36 = 1/l2
      tncoff7 = 9.D0/32768.D0*t1*t34*t36-9.D0/32768.D0*t1*(t4-t4*l8)*t36
     #/l8/(-1+l8)+9.D0/32768.D0*t1*t4*t36/l3
      t1 = -1+z
      t2 = dlog(l4)
      t3 = 1-l2
      t4 = k1(z,l1,t3,0.D0,l4)
      t7 = dlog(1-l4)
      t9 = dlog(l2)
      t12 = dlog(1-l1)
      t15 = dlog(t4)
      t19 = dlog(-t1)
      t23 = dlog(t3)
      t25 = dlog(z)
      t28 = dlog(l1)
      t31 = dlog(2.D0)
      t34 = -t2*t4-t7*t4-t9*t4-2*t12*t4-2*t4*t15-4*t19*t4+2*t4-t4*t23+2*
     #t4*t25-2*t4*t28-4*t4*t31
      t36 = 1/l2
      t38 = t1*t4
      tncoff8 = 9.D0/32768.D0*t1*t34*t36-9.D0/32768.D0*t38*t36/(-1+l8)+9
     #.D0/32768.D0*t38*t36/l3
      t3 = k1(z,l1,1-l2,0.D0,l4)
      tncoff9 = 9.D0/32768.D0*(-1+z)*t3/l2/l8/(-1+l8)
      t3 = k1(z,l1,1-l2,l3,l4)
      tncoff10 = -9.D0/32768.D0*(-1+z)*t3/l2/l3
      t3 = k1(z,l1,1-l2,l3,l4)
      tncoff11 = -9.D0/32768.D0*(-1+z)*t3/l2/l3
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = dlog(l1)
      t11 = dlog(1-l1)
      t13 = dlog(l4)
      t16 = dlog(1-l4)
      t18 = dlog(t1)
      t21 = dlog(-t2)
      t26 = 1-l8
      tncoff12 = -3.D0/65536.D0*t3*(12-24*t4+12*t6-12*t8-12*t11-6*t13-6*
     #t16-12*t18-24*t21)-9.D0/32768.D0*t3/l8-9.D0/32768.D0*t3/l3-9.D0/32
     #768.D0*t3/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = dlog(l1)
      t11 = dlog(1-l1)
      t13 = dlog(l4)
      t16 = dlog(1-l4)
      t18 = dlog(t1)
      t21 = dlog(-t2)
      tncoff13 = -3.D0/65536.D0*t3*(12-24*t4+12*t6-12*t8-12*t11-6*t13-6*
     #t16-12*t18-24*t21)+9.D0/32768.D0*t3/(-1+l8)-9.D0/32768.D0*t3/l3-9.
     #D0/32768.D0*t3/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff14 = -9.D0/32768.D0*t1*(-1+z)/(-1+l8)/l8
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff15 = 9.D0/32768.D0*(-1+z)*t2/l3
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff16 = 9.D0/32768.D0*(-1+z)*t2/l3
      t3 = k1(z,l1,1-l2,0.D0,l4)
      tncoff17 = 9.D0/32768.D0*(-1+z)*t3/l2
      t3 = k1(z,l1,1-l2,0.D0,l4)
      tncoff18 = 9.D0/32768.D0*(-1+z)*t3/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff19 = -9.D0/32768.D0*t1*(-1+z)
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff20 = -9.D0/32768.D0*t1*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,l3,l4,l5,l6,l7,0.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,l3,l4,l5,l6,l7,1.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ggrmnum1e1(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ggrmnum1e1(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ggrmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ggrmnum1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      ggrm1e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12
     #+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff
     #16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnco
     #ff20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm1e1int(var,ndim,Tcoeff)

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
      double precision ggrmnum1e0
      external ggrmnum1e0
      double precision ggrmnum1e1
      external ggrmnum1e1
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t5 = dlog(1-l4)
      t6 = dlog(2.D0)
      t8 = dlog(l4)
      t9 = dlog(z)
      t11 = dlog(l1)
      t14 = dlog(1-l1)
      t16 = dlog(t1)
      t19 = dlog(-t2)
      tncoff1 = 9.D0/32768.D0*t3*(t5+4*t6+t8-2*t9+2*t11+2*t14+2*t16+4*t1
     #9-2)-9.D0/32768.D0*t3/l8-9.D0/32768.D0*t3/l3-9.D0/32768.D0*t3/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t5 = dlog(1-l4)
      t6 = dlog(2.D0)
      t8 = dlog(l4)
      t9 = dlog(z)
      t11 = dlog(l1)
      t14 = dlog(1-l1)
      t16 = dlog(t1)
      t19 = dlog(-t2)
      tncoff2 = 9.D0/32768.D0*t3*(t5+4*t6+t8-2*t9+2*t11+2*t14+2*t16+4*t1
     #9-2)+9.D0/32768.D0*t3/(-1+l8)-9.D0/32768.D0*t3/l3-9.D0/32768.D0*t3
     #/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff3 = -9.D0/32768.D0*t1*(-1+z)/(-1+l8)/l8
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff4 = 9.D0/32768.D0*(-1+z)*t2/l3
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff5 = 9.D0/32768.D0*(-1+z)*t2/l3
      t3 = k1(z,l1,1-l2,0.D0,l4)
      tncoff6 = 9.D0/32768.D0*(-1+z)*t3/l2
      t3 = k1(z,l1,1-l2,0.D0,l4)
      tncoff7 = 9.D0/32768.D0*(-1+z)*t3/l2
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff8 = -9.D0/32768.D0*t1*(-1+z)
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff9 = -9.D0/32768.D0*t1*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum1e0(var,ndim,z,ux,l1,1-l2,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      ggrm1e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm1e2int(var,ndim,Tcoeff)

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
      double precision ggrmnum1e0
      external ggrmnum1e0
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff1 = -9.D0/32768.D0*t1*(-1+z)
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff2 = -9.D0/32768.D0*t1*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm1e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

