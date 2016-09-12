CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm15e0int(var,ndim,Tcoeff)

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

      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t7 = dlog(1-l3)
      t11 = dlog(-t3)
      t17 = dlog(l3)
      t25 = 0.3141592653589793D1**2
      t26 = 18+18*t5*t7-48*t11+24*t7*t11+72*t5*t11+18*t5*t17+24*t17*t11-
     #12*t7-36*t5-12*t17-t25
      t27 = dlog(z)
      t33 = dlog(1-l4)
      t34 = dlog(l4)
      t52 = -12*t7*t27-48*t11*t27+6*t33*t34+18*t33*t5-36*t5*t27-12*t34*t
     #27+24*t33*t11+24*t27-12*t33-12*t34+18*t34*t5+6*t34*t17
      t64 = t5**2
      t66 = t7**2
      t70 = dlog(t1)
      t77 = 6*t34*t7+6*t33*t7-12*t33*t27+6*t33*t17+6*t17*t7+27*t64+3*t66
     #+24*t34*t11+36*t5*t70-12*t17*t27-24*t70*t27
      t82 = t70**2
      t85 = t17**2
      t93 = t33**2
      t95 = t34**2
      t97 = t27**2
      t99 = t11**2
      t101 = 12*t33*t70+12*t17*t70+12*t82-24*t70+3*t85+12*t70*t7+48*t70*
     #t11+12*t34*t70+3*t93+3*t95+12*t97+48*t99
      t118 = dlog(l8)
      t125 = t33-2*t27+t7-2+2*t70-t17*l8-2*t70*l8-4*t11*l8-t34*l8-t33*l8
     #-t7*l8-3*t5*l8-t118*l8+2*t27*l8+4*t11+t34+t17+3*t5+t118+2*l8
      t127 = 1/(-1+l8)
      t129 = 1/l8
      t133 = t2*t17
      t135 = t2*t5
      t138 = t2*t70
      t141 = t2*t27
      t144 = t2*t7
      t146 = t2*t11
      t149 = t33*t2
      t151 = t2*l2
      t153 = dlog(l2)
      t154 = t153*t2
      t156 = t34*t2
      t168 = t133*l2+3*t135*l2+2*t138*l2-2*t141*l2+t144*l2+4*t146*l2+t14
     #9*l2-2*t151+t154*l2+t156*l2-6*t135-4*t138-8*t146-2*t144-2*t133-2*t
     #154-2*t149-2*t156+4*t141+4*t2
      t171 = 1/(-2+l2)
      t172 = 1/l2
      t173 = t171*t172
      t176 = 2*t2
      t177 = t2*l8
      t191 = dlog(l1)
      t196 = 1/l1
      tncoff1 = -3.D0/131072.D0*t4*(t26+t52+t77+t101)-9.D0/65536.D0*t4*t
     #125*t127*t129+9.D0/65536.D0*t3*t168*t173-9.D0/65536.D0*t3*(t176-2*
     #t177-t151+t151*l8)*t171*t172*t127*t129+9.D0/32768.D0*t3*(3*t135+2*
     #t138+4*t146+t144+t133+t149+t156-2*t141+2*t191*t2-t176)*t196-9.D0/3
     #2768.D0*t3*(-t2+t177)*t127*t196*t129-9.D0/32768.D0*(-t176+t151)*t3
     #*t173*t196
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(l3)
      t7 = dlog(1-l3)
      t10 = dlog(l4)
      t12 = dlog(-t3)
      t15 = dlog(2.D0)
      t16 = dlog(t1)
      t19 = dlog(z)
      t25 = dlog(1-l4)
      t28 = t7**2
      t30 = t16**2
      t33 = t5**2
      t35 = 18+6*t5*t7+24*t10*t12+36*t15*t16-12*t5*t19-24*t16*t19+12*t25
     #*t16+3*t28+12*t30-24*t16+3*t33
      t45 = t25**2
      t47 = t10**2
      t49 = t19**2
      t51 = t12**2
      t58 = 12*t5*t16+12*t16*t7+48*t16*t12+12*t10*t16-36*t15+3*t45+3*t47
     #+12*t49+48*t51-12*t7+18*t15*t7+24*t7*t12
      t72 = 0.3141592653589793D1**2
      t78 = 72*t15*t12+18*t15*t5+24*t5*t12-12*t7*t19-48*t12*t19-48*t12-1
     #2*t5-t72+24*t19+6*t25*t10+18*t25*t15
      t99 = t15**2
      t101 = -36*t15*t19-12*t10*t19+24*t25*t12+18*t10*t15+6*t10*t5+6*t10
     #*t7+6*t25*t7-12*t25*t19+6*t25*t5-12*t25-12*t10+27*t99
      t116 = 1-l8
      t117 = dlog(t116)
      t124 = -1/t116
      t126 = 1/l8
      t130 = t2*t5
      t132 = t2*t15
      t135 = t2*t16
      t138 = t2*t19
      t141 = t2*t7
      t143 = t2*t12
      t146 = t25*t2
      t148 = t2*l2
      t150 = dlog(l2)
      t151 = t150*t2
      t153 = t10*t2
      t165 = t130*l2+3*t132*l2+2*t135*l2-2*t138*l2+t141*l2+4*t143*l2+t14
     #6*l2-2*t148+t151*l2+t153*l2-6*t132-4*t135-8*t143-2*t141-2*t130-2*t
     #151-2*t146-2*t153+4*t138+4*t2
      t168 = 1/(-2+l2)
      t169 = 1/l2
      t170 = t168*t169
      t187 = dlog(l1)
      t190 = 2*t2
      t193 = 1/l1
      tncoff2 = -3.D0/131072.D0*t4*(t35+t58+t78+t101)-9.D0/65536.D0*t4*(
     #3*t15*l8+4*t12*l8+t7*l8-2*t19*l8+t25*l8+t5*l8+t10*l8+t117*l8+2*t16
     #*l8-2*l8)*t124*t126+9.D0/65536.D0*t3*t165*t170-9.D0/65536.D0*t3*(2
     #*t2*l8-t148*l8)*t168*t169*t124*t126+9.D0/32768.D0*t3*(3*t132+2*t13
     #5+4*t143+t141+t130+t146+t153-2*t138+2*t187*t2-t190)*t193+9.D0/3276
     #8.D0*t4*t124*t193-9.D0/32768.D0*(-t190+t148)*t3*t170*t193
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(z)
      t8 = dlog(1-l3)
      t10 = dlog(-t3)
      t12 = dlog(t1)
      t14 = dlog(l3)
      t15 = dlog(2.D0)
      t17 = dlog(l8)
      t18 = dlog(l4)
      t19 = 1-l8
      t20 = dlog(t19)
      t22 = dlog(1-l4)
      t23 = 2+2*t5-t8-4*t10-2*t12-t14-3*t15-t17-t18-t20-t22
      t25 = -1/t19
      t27 = 1/l8
      tncoff3 = -9.D0/65536.D0*t4*t23*t25*t27-9.D0/65536.D0*t3*(-2*t2+t2
     #*l2)/(-2+l2)/l2*t25*t27-9.D0/32768.D0*t4*t25/l1*t27
      t1 = -1+z
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      t7 = 2-l2
      t8 = dlog(t7)
      t11 = dlog(l3)
      t14 = dlog(t4)
      t18 = dlog(1-l3)
      t22 = dlog(1-l4)
      t25 = dlog(l2)
      t29 = dlog(-t1)
      t32 = dlog(z)
      t35 = dlog(2.D0)
      t38 = dlog(l4)
      t41 = -4*t5+2*t5*t8+2*t5*t11+4*t5*t14+2*t5*t18+2*t22*t5+2*t25*t5+8
     #*t5*t29-4*t5*t32+4*t5*t35+2*t38*t5
      t44 = -1/t7
      t45 = 1/l2
      t46 = t44*t45
      tncoff4 = 9.D0/65536.D0*t1*t41*t46-9.D0/65536.D0*t1*(-2*t5+2*t5*l8
     #)*t44*t45/(-1+l8)/l8-9.D0/16384.D0*t5*t1*t46/l1
      t1 = -1+z
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      t7 = 2-l2
      t8 = dlog(t7)
      t11 = dlog(l3)
      t14 = dlog(t4)
      t18 = dlog(1-l3)
      t22 = dlog(1-l4)
      t25 = dlog(l2)
      t29 = dlog(-t1)
      t32 = dlog(z)
      t35 = dlog(2.D0)
      t38 = dlog(l4)
      t41 = -4*t5+2*t5*t8+2*t5*t11+4*t5*t14+2*t5*t18+2*t22*t5+2*t25*t5+8
     #*t5*t29-4*t5*t32+4*t5*t35+2*t38*t5
      t46 = -1/t7/l2
      t49 = t1*t5
      tncoff5 = 9.D0/65536.D0*t1*t41*t46+9.D0/32768.D0*t49*t46/(-1+l8)-9
     #.D0/16384.D0*t49*t46/l1
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      tncoff6 = -9.D0/32768.D0*(-1+z)*t5/(-2+l2)/l2/(-1+l8)/l8
      t1 = -1+z
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      t4 = t3*l1
      t5 = dlog(l3)
      t8 = dlog(1-l1)
      t11 = dlog(2.D0)
      t14 = dlog(t2)
      t17 = 2*t3
      t19 = dlog(1-l3)
      t21 = dlog(l4)
      t22 = t21*t3
      t26 = dlog(-t1)
      t29 = dlog(z)
      t33 = t4*t5+2*t4*t8+3*t4*t11+2*t4*t14+t17-t3*t19-t22-2*t3*t8-4*t3*
     #t26+2*t3*t29-t3*t5
      t39 = dlog(1-l4)
      t40 = t39*t3
      t41 = 2*t4
      t42 = dlog(l1)
      t43 = t42*t3
      t54 = -2*t3*t14-3*t3*t11-t40-t41-2*t43-2*t4*t29+t4*t19+4*t4*t26+t2
     #2*l1+t40*l1+2*t43*l1
      t57 = 1/l1
      tncoff7 = 9.D0/32768.D0*t1*(t33+t54)*t57-9.D0/32768.D0*t1*(-t4+t4*
     #l8+t3-t3*l8)/(-1+l8)*t57/l8-9.D0/32768.D0*(-t41+t4*l2+t17-t3*l2)*t
     #1/(-2+l2)/l2*t57
      t1 = -1+z
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      t4 = t3*l1
      t5 = dlog(l3)
      t8 = dlog(1-l1)
      t11 = dlog(2.D0)
      t14 = dlog(t2)
      t17 = 2*t3
      t19 = dlog(1-l3)
      t21 = dlog(l4)
      t22 = t21*t3
      t26 = dlog(-t1)
      t29 = dlog(z)
      t33 = t4*t5+2*t4*t8+3*t4*t11+2*t4*t14+t17-t3*t19-t22-2*t3*t8-4*t3*
     #t26+2*t3*t29-t3*t5
      t39 = dlog(1-l4)
      t40 = t39*t3
      t41 = 2*t4
      t42 = dlog(l1)
      t43 = t42*t3
      t54 = -2*t3*t14-3*t3*t11-t40-t41-2*t43-2*t4*t29+t4*t19+4*t4*t26+t2
     #2*l1+t40*l1+2*t43*l1
      t57 = 1/l1
      tncoff8 = 9.D0/32768.D0*t1*(t33+t54)*t57-9.D0/32768.D0*t1*(-t4*l8+
     #t3*l8)/(-1+l8)*t57/l8-9.D0/32768.D0*(-t41+t4*l2+t17-t3*l2)*t1/(-2+
     #l2)/l2*t57
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      tncoff9 = -9.D0/32768.D0*(-1+z)*(t3*l1-t3)/(-1+l8)/l1/l8
      t3 = k1(z,l1,1-l2/2,l3,l4)
      t4 = t3**2
      tncoff10 = -9.D0/32768.D0*(2*t4*l1-2*t4)*(-1+z)/(-2+l2)/l2/l1
      t3 = k1(z,l1,1-l2/2,l3,l4)
      t4 = t3**2
      tncoff11 = -9.D0/32768.D0*(2*t4*l1-2*t4)*(-1+z)/(-2+l2)/l2/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(t1)
      t7 = dlog(l4)
      t10 = dlog(1-l4)
      t12 = dlog(z)
      t15 = dlog(-t3)
      t18 = dlog(1-l3)
      t20 = dlog(l3)
      t22 = dlog(2.D0)
      tncoff12 = -3.D0/131072.D0*t4*(12-12*t5-6*t7-6*t10+12*t12-24*t15-6
     #*t18-6*t20-18*t22)-9.D0/65536.D0*t4/l8+9.D0/65536.D0*t3*(2*t2-t2*l
     #2)/(-2+l2)/l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(t1)
      t7 = dlog(l4)
      t10 = dlog(1-l4)
      t12 = dlog(z)
      t15 = dlog(-t3)
      t18 = dlog(1-l3)
      t20 = dlog(l3)
      t22 = dlog(2.D0)
      tncoff13 = -3.D0/131072.D0*t4*(12-12*t5-6*t7-6*t10+12*t12-24*t15-6
     #*t18-6*t20-18*t22)+9.D0/65536.D0*t4/(-1+l8)+9.D0/65536.D0*t3*(2*t2
     #-t2*l2)/(-2+l2)/l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff14 = -9.D0/65536.D0*t2*(-1+z)/(-1+l8)/l8
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      tncoff15 = -9.D0/32768.D0*(-1+z)*t5/(-2+l2)/l2
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      tncoff16 = -9.D0/32768.D0*(-1+z)*t5/(-2+l2)/l2
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      tncoff17 = 9.D0/32768.D0*(-1+z)*(t3-t3*l1)/l1
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      tncoff18 = 9.D0/32768.D0*(-1+z)*(t3-t3*l1)/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff19 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff20 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum6e0(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,0.D0,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum6e0(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,1.D0,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum6e0(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum6e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum6e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum6e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = ggrmnum6e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = ggrmnum6e0(var,ndim,z,ux,l1,1-l2/2,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = ggrmnum6e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = ggrmnum6e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = ggrmnum6e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = ggrmnum6e1(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,0.D0
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = ggrmnum6e1(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,1.D0
     #,999.D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = ggrmnum6e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = ggrmnum6e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = ggrmnum6e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = ggrmnum6e2(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      ggrm15e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tnc
     #off20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm15e1int(var,ndim,Tcoeff)

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

      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t7 = dlog(l3)
      t8 = dlog(t1)
      t11 = dlog(-t3)
      t14 = dlog(1-l3)
      t15 = dlog(z)
      t18 = dlog(1-l4)
      t19 = dlog(l4)
      tncoff1 = 9.D0/65536.D0*t4*(3*t5+t7+2*t8+4*t11+t14-2-2*t15+t18+t19
     #)-9.D0/65536.D0*t4/l8-9.D0/65536.D0*t3*(-2*t2+t2*l2)/(-2+l2)/l2-9.
     #D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      t3 = -1+z
      t4 = t2*t3
      t5 = dlog(2.D0)
      t7 = dlog(l3)
      t8 = dlog(t1)
      t11 = dlog(-t3)
      t14 = dlog(1-l3)
      t15 = dlog(z)
      t18 = dlog(1-l4)
      t19 = dlog(l4)
      tncoff2 = 9.D0/65536.D0*t4*(3*t5+t7+2*t8+4*t11+t14-2-2*t15+t18+t19
     #)+9.D0/65536.D0*t4/(-1+l8)-9.D0/65536.D0*t3*(-2*t2+t2*l2)/(-2+l2)/
     #l2-9.D0/32768.D0*t4/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff3 = -9.D0/65536.D0*t2*(-1+z)/(-1+l8)/l8
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      tncoff4 = -9.D0/32768.D0*(-1+z)*t5/(-2+l2)/l2
      t4 = k1(z,0.D0,1-l2/2,l3,l4)
      t5 = t4**2
      tncoff5 = -9.D0/32768.D0*(-1+z)*t5/(-2+l2)/l2
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      tncoff6 = -9.D0/32768.D0*(-1+z)*(t3*l1-t3)/l1
      t2 = k1(z,l1,1.D0,l3,l4)
      t3 = t2**2
      tncoff7 = -9.D0/32768.D0*(-1+z)*(t3*l1-t3)/l1
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff8 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff9 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum6e0(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,0.D0,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum6e0(var,ndim,z,ux,0.D0,1-l2/2,l3,l4,l5,l6,l7,1.D0,
     #999.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum6e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,0.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum6e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,1.D0,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = ggrmnum6e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = ggrmnum6e1(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      ggrm15e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm15e2int(var,ndim,Tcoeff)

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

      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff1 = -9.D0/65536.D0*t2*(-1+z)
      t1 = k1(z,0.D0,1.D0,l3,l4)
      t2 = t1**2
      tncoff2 = -9.D0/65536.D0*t2*(-1+z)

CCC Numerator functions

      tcnum1 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,0.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum6e0(var,ndim,z,ux,0.D0,1.D0,l3,l4,l5,l6,l7,1.D0,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm15e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

