CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm14e0int(var,ndim,Tcoeff)

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
      double precision ggrmnum6e0
      external ggrmnum6e0
      double precision ggrmnum6e1
      external ggrmnum6e1
      double precision ggrmnum6e2
      external ggrmnum6e2
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

      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t8 = dlog(1-l3)
      t13 = dlog(-t3)
      t18 = dlog(l3)
      t25 = dlog(z)
      t28 = 18-36*t5-12*t8+18*t5*t8+24*t8*t13+72*t5*t13+18*t5*t18-48*t13
     #-12*t18+24*t18*t13-12*t8*t25
      t29 = dlog(t1)
      t38 = 0.3141592653589793D1**2
      t40 = dlog(1-l4)
      t42 = dlog(l4)
      t52 = 12*t29*t8-48*t13*t25-24*t29*t25-24*t29+24*t25-t38-12*t40-12*
     #t42+6*t40*t42+18*t40*t5-36*t5*t25+12*t18*t29
      t60 = t5**2
      t76 = 48*t29*t13-12*t42*t25+24*t40*t13+27*t60+18*t42*t5+6*t42*t18+
     #6*t42*t8+12*t42*t29+6*t40*t8-12*t40*t25+6*t40*t18
      t83 = t8**2
      t85 = t29**2
      t89 = t40**2
      t93 = t18**2
      t95 = t42**2
      t97 = t25**2
      t99 = t13**2
      t101 = 24*t42*t13+36*t5*t29+6*t18*t8+3*t83+12*t85-12*t18*t25+3*t89
     #+12*t40*t29+3*t93+3*t95+12*t97+48*t99
      t106 = dlog(l8)
      t125 = -2-t106*l8-t42*l8-t18*l8-t8*l8-3*t5*l8-t40*l8-4*t13*l8+2*t2
     #5*l8-2*t29*l8+t40+t42-2*t25+4*t13+t8+2*t29+t18+t106+2*l8+3*t5
      t127 = 1/(-1+l8)
      t129 = 1/l8
      t133 = t2*t25
      t136 = t2*t13
      t139 = t2*t8
      t141 = t2*t18
      t143 = t2*t29
      t146 = t2*l2
      t148 = t2*t5
      t151 = dlog(l2)
      t152 = t151*t2
      t154 = t42*t2
      t156 = t40*t2
      t168 = -2*t133*l2+4*t136*l2+t139*l2+t141*l2+2*t143*l2-2*t146+3*t14
     #8*l2+t152*l2+t154*l2+t156*l2-2*t154-2*t156-8*t136-6*t148-2*t141-4*
     #t143+4*t133-2*t139-2*t152+4*t2
      t171 = 1/(-2+l2)
      t172 = 1/l2
      t176 = 2*t2
      t177 = t2*l8
      t191 = dlog(l1)
      t196 = 1/l1
      tncoff1 = -3.D0/131072.D0*t4*(t28+t52+t76+t101)-9.D0/65536.D0*t4*t
     #125*t127*t129+9.D0/65536.D0*t3*t168*t171*t172-9.D0/65536.D0*t3*(t1
     #76-2*t177-t146+t146*l8)*t171*t172*t129*t127+9.D0/32768.D0*t3*(t154
     #+t156+4*t136+3*t148+t141+2*t143-2*t133+t139-t176+2*t191*t2)*t196+9
     #.D0/32768.D0*t3*(t2-t177)*t196*t129*t127-9.D0/32768.D0*t3*(-t176+t
     #146)*t171*t196*t172
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t8 = dlog(1-l3)
      t11 = dlog(-t3)
      t13 = dlog(l3)
      t15 = 0.3141592653589793D1**2
      t16 = dlog(z)
      t19 = dlog(t1)
      t27 = dlog(1-l4)
      t28 = dlog(l4)
      t31 = 18-36*t5-12*t8-48*t11-12*t13-t15-12*t8*t16+12*t19*t8-48*t11*
     #t16-24*t19*t16+6*t27*t28
      t52 = 18*t27*t5-36*t5*t16+12*t13*t19+48*t19*t11+24*t16+6*t28*t13-2
     #4*t19+18*t28*t5-12*t27-12*t28+24*t27*t11-12*t28*t16
      t56 = t5**2
      t74 = t8**2
      t76 = 6*t28*t8+27*t56+12*t28*t19+6*t27*t8-12*t27*t16+6*t27*t13+36*
     #t5*t19+6*t8*t13+24*t28*t11+18*t5*t13+3*t74
      t77 = t19**2
      t81 = t13**2
      t83 = t27**2
      t87 = t28**2
      t89 = t16**2
      t91 = t11**2
      t101 = 12*t77-12*t13*t16+3*t81+3*t83+12*t27*t19+3*t87+12*t89+48*t9
     #1+18*t5*t8+24*t13*t11+72*t5*t11+24*t8*t11
      t109 = 1-l8
      t110 = dlog(t109)
      t124 = -1/t109
      t126 = 1/l8
      t130 = t2*t16
      t133 = t2*t11
      t136 = t2*t8
      t138 = t2*t13
      t140 = t2*t19
      t143 = t2*l2
      t145 = t2*t5
      t148 = dlog(l2)
      t149 = t148*t2
      t151 = t28*t2
      t153 = t27*t2
      t165 = -2*t130*l2+4*t133*l2+t136*l2+t138*l2+2*t140*l2-2*t143+3*t14
     #5*l2+t149*l2+t151*l2+t153*l2-2*t151-2*t153-8*t133-6*t145-2*t138-4*
     #t140+4*t130-2*t136-2*t149+4*t2
      t168 = 1/(-2+l2)
      t169 = 1/l2
      t187 = 2*t2
      t188 = dlog(l1)
      t193 = 1/l1
      tncoff2 = -3.D0/131072.D0*t4*(t31+t52+t76+t101)-9.D0/65536.D0*t4*(
     #t8*l8+t13*l8+t27*l8+t110*l8+3*t5*l8+2*t19*l8+4*t11*l8-2*t16*l8+t28
     #*l8-2*l8)*t124*t126+9.D0/65536.D0*t3*t165*t168*t169-9.D0/65536.D0*
     #t3*(2*t2*l8-t143*l8)*t168*t169*t126*t124+9.D0/32768.D0*t3*(t151+t1
     #53+4*t133+3*t145+t138+2*t140-2*t130+t136-t187+2*t188*t2)*t193+9.D0
     #/32768.D0*t4*t193*t124-9.D0/32768.D0*t3*(-t187+t143)*t168*t193*t16
     #9
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(l8)
      t7 = dlog(-t3)
      t9 = dlog(t1)
      t11 = dlog(l3)
      t12 = dlog(2.D0)
      t15 = dlog(1-l3)
      t16 = dlog(z)
      t19 = dlog(1-l4)
      t20 = 1-l8
      t21 = dlog(t20)
      t22 = dlog(l4)
      t23 = 2-t5-4*t7-2*t9-t11-3*t12-t15+2*t16-t19-t21-t22
      t25 = -1/t20
      t27 = 1/l8
      tncoff3 = -9.D0/65536.D0*t4*t23*t25*t27-9.D0/65536.D0*t3*(-2*t2+t2
     #*l2)/(-2+l2)/l2*t27*t25-9.D0/32768.D0*t4/l1*t27*t25
      t1 = -1+z
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      t6 = dlog(-t1)
      t9 = dlog(z)
      t13 = dlog(1-l3)
      t16 = 2-l2
      t17 = dlog(t16)
      t20 = dlog(t3)
      t23 = dlog(2.D0)
      t26 = dlog(l3)
      t30 = dlog(1-l4)
      t33 = dlog(l2)
      t36 = dlog(l4)
      t40 = 8*t4*t6-4*t4*t9+2*t4*t13+2*t4*t17+4*t4*t20+4*t4*t23+2*t4*t26
     #+2*t30*t4+2*t33*t4+2*t36*t4-4*t4
      t43 = -1/t16
      t44 = 1/l2
      tncoff4 = 9.D0/65536.D0*t1*t40*t43*t44-9.D0/65536.D0*t1*(-2*t4+2*t
     #4*l8)*t43*t44/l8/(-1+l8)-9.D0/16384.D0*t1*t4*t43/l1*t44
      t1 = -1+z
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      t6 = dlog(-t1)
      t9 = dlog(z)
      t13 = dlog(1-l3)
      t16 = 2-l2
      t17 = dlog(t16)
      t20 = dlog(t3)
      t23 = dlog(2.D0)
      t26 = dlog(l3)
      t30 = dlog(1-l4)
      t33 = dlog(l2)
      t36 = dlog(l4)
      t40 = 8*t4*t6-4*t4*t9+2*t4*t13+2*t4*t17+4*t4*t20+4*t4*t23+2*t4*t26
     #+2*t30*t4+2*t33*t4+2*t36*t4-4*t4
      t43 = -1/t16
      t44 = 1/l2
      t45 = t43*t44
      t48 = t1*t4
      tncoff5 = 9.D0/65536.D0*t1*t40*t45+9.D0/32768.D0*t48*t45/(-1+l8)-9
     #.D0/16384.D0*t48*t43/l1*t44
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      tncoff6 = -9.D0/32768.D0*(-1+z)*t4/(-2+l2)/l2/l8/(-1+l8)
      t1 = -1+z
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      t4 = t3*l1
      t6 = dlog(1-l3)
      t9 = dlog(-t1)
      t12 = 2*t3
      t13 = dlog(2.D0)
      t16 = dlog(l3)
      t18 = dlog(t2)
      t23 = dlog(1-l1)
      t29 = dlog(1-l4)
      t30 = t29*t3
      t31 = dlog(l4)
      t32 = t31*t3
      t33 = t4*t6+4*t4*t9+t12-3*t3*t13-t3*t16-2*t3*t18-t3*t6-2*t3*t23-4*
     #t3*t9-t30-t32
      t34 = dlog(l1)
      t35 = t34*t3
      t37 = dlog(z)
      t53 = 2*t4
      t54 = -2*t35+2*t3*t37+t32*l1+2*t35*l1-2*t4*t37+3*t4*t13+2*t4*t23+t
     #4*t16+2*t4*t18+t30*l1-t53
      t57 = 1/l1
      tncoff7 = 9.D0/32768.D0*t1*(t33+t54)*t57+9.D0/32768.D0*t1*(t4-t4*l
     #8-t3+t3*l8)*t57/l8/(-1+l8)-9.D0/32768.D0*t1*(-t53+t4*l2+t12-t3*l2)
     #/(-2+l2)*t57/l2
      t1 = -1+z
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      t4 = t3*l1
      t6 = dlog(1-l3)
      t9 = dlog(-t1)
      t12 = 2*t3
      t13 = dlog(2.D0)
      t16 = dlog(l3)
      t18 = dlog(t2)
      t23 = dlog(1-l1)
      t29 = dlog(1-l4)
      t30 = t29*t3
      t31 = dlog(l4)
      t32 = t31*t3
      t33 = t4*t6+4*t4*t9+t12-3*t3*t13-t3*t16-2*t3*t18-t3*t6-2*t3*t23-4*
     #t3*t9-t30-t32
      t34 = dlog(l1)
      t35 = t34*t3
      t37 = dlog(z)
      t53 = 2*t4
      t54 = -2*t35+2*t3*t37+t32*l1+2*t35*l1-2*t4*t37+3*t4*t13+2*t4*t23+t
     #4*t16+2*t4*t18+t30*l1-t53
      t57 = 1/l1
      tncoff8 = 9.D0/32768.D0*t1*(t33+t54)*t57+9.D0/32768.D0*t1*(t4*l8-t
     #3*l8)*t57/l8/(-1+l8)-9.D0/32768.D0*t1*(-t53+t4*l2+t12-t3*l2)/(-2+l
     #2)*t57/l2
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      tncoff9 = 9.D0/32768.D0*(-1+z)*(-t3*l1+t3)/l1/l8/(-1+l8)
      t3 = k1(z,l1,l2/2,l3,l4)
      t4 = t3**2
      tncoff10 = -9.D0/32768.D0*(-1+z)*(2*t4*l1-2*t4)/(-2+l2)/l1/l2
      t3 = k1(z,l1,l2/2,l3,l4)
      t4 = t3**2
      tncoff11 = -9.D0/32768.D0*(-1+z)*(2*t4*l1-2*t4)/(-2+l2)/l1/l2
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(l3)
      t7 = dlog(t1)
      t10 = dlog(-t3)
      t13 = dlog(1-l3)
      t15 = dlog(z)
      t18 = dlog(1-l4)
      t20 = dlog(l4)
      t22 = dlog(2.D0)
      tncoff12 = -3.D0/131072.D0*t4*(-6*t5-12*t7-24*t10-6*t13+12*t15-6*t
     #18-6*t20-18*t22+12)-9.D0/65536.D0*t4/l8+9.D0/65536.D0*t3*(2*t2-t2*
     #l2)/(-2+l2)/l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(l3)
      t7 = dlog(t1)
      t10 = dlog(-t3)
      t13 = dlog(1-l3)
      t15 = dlog(z)
      t18 = dlog(1-l4)
      t20 = dlog(l4)
      t22 = dlog(2.D0)
      tncoff13 = -3.D0/131072.D0*t4*(-6*t5-12*t7-24*t10-6*t13+12*t15-6*t
     #18-6*t20-18*t22+12)+9.D0/65536.D0*t4/(-1+l8)+9.D0/65536.D0*t3*(2*t
     #2-t2*l2)/(-2+l2)/l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff14 = -9.D0/65536.D0*t2*(-1+z)/(-1+l8)/l8
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      tncoff15 = -9.D0/32768.D0*(-1+z)*t4/(-2+l2)/l2
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      tncoff16 = -9.D0/32768.D0*(-1+z)*t4/(-2+l2)/l2
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      tncoff17 = 9.D0/32768.D0*(-1+z)*(-t3*l1+t3)/l1
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      tncoff18 = 9.D0/32768.D0*(-1+z)*(-t3*l1+t3)/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff19 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff20 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum6e0(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum6e0(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum6e0(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum6e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum6e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum6e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ggrmnum6e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,0.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ggrmnum6e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,1.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ggrmnum6e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ggrmnum6e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ggrmnum6e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ggrmnum6e1(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ggrmnum6e1(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ggrmnum6e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ggrmnum6e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ggrmnum6e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ggrmnum6e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      ggrm14e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm14e1int(var,ndim,Tcoeff)

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
      double precision ggrmnum6e0
      external ggrmnum6e0
      double precision ggrmnum6e1
      external ggrmnum6e1
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

      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t7 = dlog(l3)
      t8 = dlog(t1)
      t11 = dlog(1-l3)
      t13 = dlog(-t3)
      t15 = dlog(z)
      t17 = dlog(l4)
      t19 = dlog(1-l4)
      t23 = 1-l8
      tncoff1 = 9.D0/65536.D0*t4*(3*t5+t7+2*t8+t11+4*t13-2-2*t15+t17+t19
     #)-9.D0/65536.D0*t4/l8-9.D0/65536.D0*t3*(-2*t2+t2*l2)/(-2+l2)/l2-9.
     #D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t7 = dlog(l3)
      t8 = dlog(t1)
      t11 = dlog(1-l3)
      t13 = dlog(-t3)
      t15 = dlog(z)
      t17 = dlog(l4)
      t19 = dlog(1-l4)
      tncoff2 = 9.D0/65536.D0*t4*(3*t5+t7+2*t8+t11+4*t13-2-2*t15+t17+t19
     #)+9.D0/65536.D0*t4/(-1+l8)-9.D0/65536.D0*t3*(-2*t2+t2*l2)/(-2+l2)/
     #l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff3 = -9.D0/65536.D0*t2*(-1+z)/(-1+l8)/l8
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      tncoff4 = -9.D0/32768.D0*(-1+z)*t4/(-2+l2)/l2
      t3 = k1(z,0.D0,l2/2,l3,l4)
      t4 = t3**2
      tncoff5 = -9.D0/32768.D0*(-1+z)*t4/(-2+l2)/l2
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      tncoff6 = -9.D0/32768.D0*(-1+z)*(t3*l1-t3)/l1
      t2 = k1(z,l1,0.D0,l3,l4)
      t3 = t2**2
      tncoff7 = -9.D0/32768.D0*(-1+z)*(t3*l1-t3)/l1
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff8 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff9 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum6e0(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum6e0(var,ndim,z,ux,0.D0,l2/2,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum6e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum6e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum6e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum6e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      ggrm14e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm14e2int(var,ndim,Tcoeff)

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
      double precision ggrmnum6e0
      external ggrmnum6e0
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

      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff1 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,0.D0,l3,l4)
      t2 = t1**2
      tncoff2 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm14e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

