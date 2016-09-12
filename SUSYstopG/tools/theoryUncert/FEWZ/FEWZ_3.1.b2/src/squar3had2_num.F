CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq3num1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = -1+z
      t3 = t1*t2
      t27 = t2**2
      t29 = z*l3
      t30 = s1vr**2
      t46 = t1**2
      t47 = t27*t2*t46
      t52 = -t1*t2
      t53 = s4l2*s3l1
      t56 = s3l2*s4l1
      t66 = t2*l3
      t85 = l3*t27
      t89 = -512.D0/3.D0*t85*t1*s34r-512.D0/3.D0+512.D0/3.D0*z
      t91 = t1*t27
      t102 = 2*z-1
      t136 = l3**2
      t145 = z*t27
      t157 = 1/s234r
      t158 = t157*s3l2
      t159 = t3*t158
      t160 = t157*s4l2
      t161 = t3*t160
      t163 = -64.D0/3.D0*t159-64.D0/3.D0*t161
      t165 = t27*t46
      t166 = t157*s34r
      t168 = t165*t166*s2l2
      t173 = 64.D0/3.D0*t165*t158+64.D0/3.D0*t165*t160
      t179 = t157*s3l1
      t181 = t157*s4l1
      t184 = -64.D0/3.D0*t3*t179-64.D0/3.D0*t3*t181
      t192 = 64.D0/3.D0*t165*t179+64.D0/3.D0*t165*t181
      t198 = t85*t1
      t199 = z*t157
      t244 = t1*t157*s34r
      t246 = t2*t157
      t248 = -64.D0/3.D0*t85*t244+64.D0/3.D0*t246
      t250 = t160*s34r
      t278 = t181*s34r
      t294 = s23r**2
      t307 = -128.D0/3.D0*t158-64.D0/3.D0*t160
      t309 = t181*s3l2
      t310 = 64.D0/3.D0*t309
      t330 = t3*t157
      s2 = ((-512.D0/3.D0*t3*s4l2*s1vr+(-512.D0/3.D0*t3*s2l2-512.D0/3.D0
     #*t3*s3l2-512.D0/3.D0*t3*s4l2)*s23r)*s1l1+(-512.D0/3.D0*t3*s4l1*s1v
     #r+(-512.D0/3.D0*t3*s2l1-512.D0/3.D0*t3*s3l1-512.D0/3.D0*t3*s4l1)*s
     #23r)*s1l2-512.D0/3.D0*t1*t27*t29*t30+(1024.D0/3.D0*t1*t27*t29*s23r
     #-512.D0/3.D0*t3*s4l2*s2l1-512.D0/3.D0*t3*s4l1*s2l2+512.D0/3.D0*t47
     #*t29*s34r+512.D0/3.D0*t52*t53+512.D0/3.D0*t52*t56-1024.D0/3.D0*t3*
     #s4l2*s4l1)*s1vr)*s13r+((1024.D0/3.D0*t66*t30-1024.D0/3.D0*t66*s23r
     #*s1vr)*s1l2+(512.D0/3.D0*t2*s2l2*l3+1024.D0/3.D0*t66*s3l2+1024.D0/
     #3.D0*t66*s4l2)*t30+(-512.D0/3.D0*t66*s3l2*s23r+t89*s2l2+(-512.D0/3
     #.D0*t91*l3*s3l2-512.D0/3.D0*t91*l3*s4l2)*s34r+512.D0/3.D0*t2*s3l2+
     #512.D0/3.D0*t102*s4l2)*s1vr)*s1l1
      s1 = s2+((512.D0/3.D0*t66*s2l1+1024.D0/3.D0*t66*s3l1+1024.D0/3.D0*
     #t66*s4l1)*t30+(-512.D0/3.D0*t66*s3l1*s23r+t89*s2l1+(-512.D0/3.D0*t
     #91*l3*s3l1-512.D0/3.D0*t91*l3*s4l1)*s34r+512.D0/3.D0*t2*s3l1+512.D
     #0/3.D0*t102*s4l1)*s1vr)*s1l2+(-512.D0/3.D0*t27*t136*z*s23r+512.D0/
     #3.D0*t66*t53+512.D0/3.D0*t66*t56)*t30+(-512.D0/3.D0*t145*l3*s23r+5
     #12.D0/3.D0*t29*t91*s34r+512.D0/3.D0*z*s4l2*s3l1+512.D0/3.D0*z*s4l1
     #*s3l2)*s1vr
      s2 = s1+((t163*t30+(64.D0/3.D0*t168+t173*s34r)*s1vr)*s1l1+(t184*t3
     #0+(64.D0/3.D0*t165*t166*s2l1+t192*s34r)*s1vr)*s1l2+(64.D0/3.D0*t19
     #8*t199*s23r-t163*s2l1-t184*s2l2)*t30+(-64.D0/3.D0*t47*l3*t199*s34r
     #*s23r+(-128.D0/3.D0*t168-t173*s34r)*s2l1-t192*s34r*s2l2)*s1vr)*s13
     #r+(128.D0/3.D0*t66*t157*s23r*t30*s1l2+((64.D0/3.D0*t66*t157*s2l2+6
     #4.D0/3.D0*t66*t158+128.D0/3.D0*t66*t160)*s23r-128.D0/3.D0*t160)*t3
     #0+((t248*s2l2-64.D0/3.D0*t198*t250+64.D0/3.D0*t199*s3l2+64.D0/3.D0
     #*t199*s4l2)*s23r+64.D0/3.D0*t3*t250)*s1vr)*s1l1
      taux1 = s2+(((64.D0/3.D0*t66*t157*s2l1+64.D0/3.D0*t66*t179+128.D0/
     #3.D0*t66*t181)*s23r-128.D0/3.D0*t181)*t30+((t248*s2l1-64.D0/3.D0*t
     #198*t278+64.D0/3.D0*t199*s3l1+64.D0/3.D0*t199*s4l1)*s23r+64.D0/3.D
     #0*t3*t278)*s1vr)*s1l2+(-64.D0/3.D0*t145*t136*t157*t294+(-64.D0/3.D
     #0*t66*t158*s2l1-64.D0/3.D0*t66*t179*s2l2)*s23r+t307*s3l1-t310)*t30
     #+(-64.D0/3.D0*t145*l3*t157*t294+(-64.D0/3.D0*t246*s3l2*s2l1-64.D0/
     #3.D0*t246*s3l1*s2l2+64.D0/3.D0*t29*t27*t244-t307*s3l1+t310)*s23r+6
     #4.D0/3.D0*t330*s3l2*s34r*s2l1+64.D0/3.D0*t330*s3l1*s34r*s2l2+((128
     #.D0/3.D0*t159+64.D0/3.D0*t161)*s3l1+64.D0/3.D0*t3*t309)*s34r)*s1vr
      t1 = 1-z
      t2 = l3-1
      t3 = t1*t2
      t7 = -t1
      t8 = t2*t7
      t28 = s3l2*s2l1
      t30 = s3l1*s2l2
      t32 = s4l2*s3l1
      t34 = s3l2*s4l1
      t44 = t7*l3
      t50 = s1vr**2
      t55 = 2048.D0/3.D0*t7
      t57 = t7**2
      t58 = t2*t57
      t63 = t1*s3l2
      t66 = -2*z+1
      t90 = t1*s3l1
      t98 = t1*l3
      t115 = 1/s234r
      t116 = t115*s3l2
      t118 = t115*s4l2
      t119 = t3*t118
      t121 = -256.D0/3.D0*t8*t116+256.D0/3.D0*t119
      t123 = t2**2
      t124 = t57*t123
      t125 = t115*s34r
      t132 = 256.D0/3.D0*t124*t116+256.D0/3.D0*t124*t118
      t138 = t115*s3l1
      t140 = t115*s4l1
      t141 = t3*t140
      t143 = 256.D0/3.D0*t8*t138-256.D0/3.D0*t141
      t151 = -256.D0/3.D0*t124*t138-256.D0/3.D0*t124*t140
      t186 = t98*t118
      t196 = s23r**2
      t198 = t1*t115
      t201 = l3*t115
      t205 = -t2*t57
      t211 = z*t115
      t230 = t98*t140
      t263 = t118*s3l1
      t264 = t140*s3l2
      t268 = t3*t115
      s2 = ((2048.D0/3.D0*t3*s4l2*s1vr+(-2048.D0/3.D0*t8*s2l2-2048.D0/3.
     #D0*t8*s3l2-2048.D0/3.D0*t8*s4l2)*s23r)*s1l1+(-2048.D0/3.D0*t3*s4l1
     #*s1vr+(2048.D0/3.D0*t8*s2l1+2048.D0/3.D0*t8*s3l1+2048.D0/3.D0*t8*s
     #4l1)*s23r)*s1l2+(-2048.D0/3.D0*t3*t28+2048.D0/3.D0*t3*t30+2048.D0/
     #3.D0*t3*t32-2048.D0/3.D0*t3*t34)*s1vr)*s13r+((2048.D0/3.D0*t7*s2l2
     #*l3+2048.D0/3.D0*t44*s3l2+4096.D0/3.D0*t44*s4l2)*t50+(-2048.D0/3.D
     #0*t44*s4l2*s23r+t55*s2l2-2048.D0/3.D0*t58*l3*s4l2*s34r-2048.D0/3.D
     #0*t63-2048.D0/3.D0*t66*s4l2)*s1vr)*s1l1
      s1 = s2+((-2048.D0/3.D0*t44*s2l1-2048.D0/3.D0*t44*s3l1-4096.D0/3.D
     #0*t44*s4l1)*t50+(2048.D0/3.D0*t44*s4l1*s23r-t55*s2l1+2048.D0/3.D0*
     #t58*l3*s4l1*s34r+2048.D0/3.D0*t90+2048.D0/3.D0*t66*s4l1)*s1vr)*s1l
     #2+(2048.D0/3.D0*t98*t28-2048.D0/3.D0*t98*t30+2048.D0/3.D0*t44*t32-
     #2048.D0/3.D0*t44*t34)*t50+(2048.D0/3.D0*t63*s2l1-2048.D0/3.D0*t90*
     #s2l2+2048.D0/3.D0*t7*s4l2*s3l1-2048.D0/3.D0*t7*s4l1*s3l2)*s1vr
      s2 = s1+((t121*t50+(256.D0/3.D0*t124*t125*s2l2+t132*s34r)*s1vr)*s1
     #l1+(t143*t50+(-256.D0/3.D0*t124*t125*s2l1+t151*s34r)*s1vr)*s1l2+((
     #-256.D0/3.D0*t3*t116-256.D0/3.D0*t119)*s2l1+(256.D0/3.D0*t3*t138+2
     #56.D0/3.D0*t141)*s2l2)*t50+((t121*s2l1+t143*s2l2)*s23r-t132*s34r*s
     #2l1-t151*s34r*s2l2)*s1vr)*s13r+(((-256.D0/3.D0*t98*t115*s2l2+512.D
     #0/3.D0*t44*t116-512.D0/3.D0*t186)*s23r-512.D0/3.D0*t118)*t50+((256
     #.D0/3.D0*t98*t116+256.D0/3.D0*t186)*t196+(-256.D0/3.D0*t198*s2l2+(
     #-256.D0/3.D0*t58*t201*s3l2+256.D0/3.D0*t205*t201*s4l2)*s34r+256.D0
     #/3.D0*t211*s3l2+256.D0/3.D0*t211*s4l2)*s23r-256.D0/3.D0*t3*t118*s3
     #4r)*s1vr)*s1l1
      taux2 = s2+(((256.D0/3.D0*t98*t115*s2l1-512.D0/3.D0*t44*t138+512.D
     #0/3.D0*t230)*s23r+512.D0/3.D0*t140)*t50+((-256.D0/3.D0*t98*t138-25
     #6.D0/3.D0*t230)*t196+(256.D0/3.D0*t198*s2l1+(256.D0/3.D0*t58*t201*
     #s3l1-256.D0/3.D0*t205*t201*s4l1)*s34r-256.D0/3.D0*t211*s3l1-256.D0
     #/3.D0*t211*s4l1)*s23r+256.D0/3.D0*t3*t140*s34r)*s1vr)*s1l2+(-256.D
     #0/3.D0*t263+256.D0/3.D0*t264)*t50+(256.D0/3.D0*t268*s3l2*s34r*s2l1
     #-256.D0/3.D0*t268*s3l1*s34r*s2l2+(256.D0/3.D0*t8*t263-256.D0/3.D0*
     #t8*t264)*s34r)*s1vr

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux2
      tcoff6 = taux2
      tcoff7 = -taux2
      tcoff8 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet2 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet3 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet4 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet6 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet7 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet8 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      sq3num1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function sq3num1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = z*t2
      t4 = l3-1
      t5 = t4**2
      t7 = s13r**2
      t12 = -t4*t1
      t14 = t12*s3l2
      t15 = t4*t1
      t16 = t15*s4l2
      t20 = -t1
      t22 = 2*t20*t4
      t32 = t12*s3l1
      t33 = t15*s4l1
      t45 = 3*t4
      t47 = z*l3
      t48 = s1vr**2
      t63 = t2*t1
      t64 = t63*t5
      t65 = t47*s34r
      t74 = 2*z
      t75 = -t74+1
      t81 = z**2
      t82 = t81*t1
      t88 = l3*t1
      t98 = t88*s3l2
      t100 = t88*s4l2
      t104 = t20*l3
      t109 = l3*t2
      t114 = 1024.D0/3.D0*t109*t4*s34r+512.D0/3.D0-512.D0/3.D0*z
      t116 = t4*t2
      t127 = -t75
      t136 = t88*s3l1
      t163 = l3**2
      t164 = t2*t163
      t165 = t48*s1vr
      t172 = s3l2*s2l1
      t175 = s3l1*s2l2
      t188 = s3l2*s4l1
      t199 = t20*s3l2
      t211 = -1-z
      t222 = 1/s234r
      t223 = t222*s4l2
      t227 = t222*s2l2
      t230 = t2*t5
      t231 = t223*s34r
      t238 = t222*s4l1
      t242 = t222*s2l1
      t245 = t238*s34r
      t252 = t116*l3
      t253 = z*t222
      t257 = t253*s23r
      t262 = 64.D0/3.D0*t15*t223*s2l1
      t265 = 64.D0/3.D0*t15*t238*s2l2
      t266 = t64*l3
      t270 = t222*s3l2
      t274 = 64.D0/3.D0*t15*t223
      t277 = t238*s3l2
      t279 = 64.D0/3.D0*t15*t277
      t287 = s23r**2
      t291 = t223*s3l1
      t299 = t230*t222
      t301 = s3l2*s34r*s2l1
      t305 = s3l1*s34r*s2l2
      t308 = t2**2
      t312 = l3*t222
      t313 = s34r**2
      t350 = t88*t270
      t352 = t88*t223
      t366 = t1*t222
      t372 = t4*t222
      t382 = t222*t313
      t397 = t222*s3l1
      t398 = t88*t397
      t400 = t88*t238
      t441 = t270*s3l1
      t495 = t20*t222
      t507 = t116*t222
      s3 = -512.D0/3.D0*t3*t5*s1vr*t7
      s4 = (((512.D0/3.D0*t12*s2l2+512.D0/3.D0*t14+512.D0/3.D0*t16)*s1vr
     #+(-512.D0/3.D0*t22*s2l2-512.D0/3.D0*t22*s3l2-512.D0/3.D0*t22*s4l2)
     #*s23r)*s1l1+((512.D0/3.D0*t12*s2l1+512.D0/3.D0*t32+512.D0/3.D0*t33
     #)*s1vr+(-512.D0/3.D0*t22*s2l1-512.D0/3.D0*t22*s3l1-512.D0/3.D0*t22
     #*s4l1)*s23r)*s1l2+512.D0/3.D0*t45*t2*t47*t48+(-512.D0/3.D0*t45*t2*
     #t47*s23r+(512.D0/3.D0*t14+512.D0/3.D0*t16)*s2l1+(512.D0/3.D0*t32+5
     #12.D0/3.D0*t33)*s2l2-512*t64*t65+1024.D0/3.D0*t12*s3l2*s3l1+1024.D
     #0/3.D0*t15*s4l2*s4l1-512.D0/3.D0*t75*z*t15)*s1vr+512.D0/3.D0*t82*t
     #4*s23r)*s13r+((-1024.D0/3.D0*t88*t48+1024.D0/3.D0*t88*s23r*s1vr)*s
     #1l2+(-512.D0/3.D0*t1*s2l2*l3-1024.D0/3.D0*t98-1024.D0/3.D0*t100)*t
     #48+((512.D0/3.D0*t98+512.D0/3.D0*t104*s4l2)*s23r+t114*s2l2+(1024.D
     #0/3.D0*t116*l3*s3l2+512.D0/3.D0*t116*l3*s4l2)*s34r-512.D0/3.D0*t1*
     #s3l2-512.D0/3.D0*t127*s4l2)*s1vr)*s1l1
      s2 = s3+s4
      s1 = s2+((-512.D0/3.D0*t88*s2l1-1024.D0/3.D0*t136-1024.D0/3.D0*t88
     #*s4l1)*t48+((512.D0/3.D0*t136+512.D0/3.D0*t104*s4l1)*s23r+t114*s2l
     #1+(1024.D0/3.D0*t116*l3*s3l1+512.D0/3.D0*t116*l3*s4l1)*s34r-512.D0
     #/3.D0*t1*s3l1-512.D0/3.D0*t127*s4l1)*s1vr)*s1l2-1024.D0/3.D0*t164*
     #z*t165+(512*t164*z*s23r-512.D0/3.D0*t104*t172-512.D0/3.D0*t104*t17
     #5+1024.D0/3.D0*z*t4*t63*t163*s34r+(-1024.D0/3.D0*t104*s3l2-512.D0/
     #3.D0*t100)*s3l1-512.D0/3.D0*t88*t188+1024.D0/3.D0*t75*t1*t47)*t48
      s3 = s1+(512*t3*l3*s23r-512.D0/3.D0*t199*s2l1-512.D0/3.D0*t20*s3l1
     #*s2l2+512.D0/3.D0*(-3+t74)*t4*t2*t65+(-1024.D0/3.D0*t199+512.D0/3.
     #D0*t211*s4l2)*s3l1+512.D0/3.D0*t211*s4l1*s3l2-1024.D0/3.D0*t82)*s1
     #vr
      s4 = s3
      s6 = ((64.D0/3.D0*t15*t223*t48+(64.D0/3.D0*t15*t227*s23r-64.D0/3.D
     #0*t230*t231)*s1vr)*s1l1+(64.D0/3.D0*t15*t238*t48+(64.D0/3.D0*t15*t
     #242*s23r-64.D0/3.D0*t230*t245)*s1vr)*s1l2+128.D0/3.D0*t252*t253*t1
     #65+(128.D0/3.D0*t252*t257-t262-t265-256.D0/3.D0*t266*t253*s34r+(-1
     #28.D0/3.D0*t15*t270-t274)*s3l1-t279-64.D0/3.D0*(3-z)*z*t15*t222)*t
     #48+(-128.D0/3.D0*t252*t253*t287+(t262+t265+64.D0/3.D0*t15*t291+t27
     #9+128.D0/3.D0*t15*t223*s4l1)*s23r+64.D0/3.D0*t299*t301+64.D0/3.D0*
     #t299*t305+128.D0/3.D0*z*t308*t5*t4*t312*t313+((128.D0/3.D0*t230*t2
     #70+64.D0/3.D0*t230*t223)*s3l1+64.D0/3.D0*t230*t277+64.D0/3.D0*t20*
     #z*t299)*s34r)*s1vr)*s13r
      s7 = ((-128.D0/3.D0*t88*t222*s23r*t48+(128.D0/3.D0*t88*t222*t287+1
     #28.D0/3.D0*t252*t222*s34r*s23r)*s1vr)*s1l2+((-64.D0/3.D0*t88*t227-
     #64.D0/3.D0*t350-128.D0/3.D0*t352)*s23r+64.D0/3.D0*t252*t270*s34r+1
     #28.D0/3.D0*t223)*t48+((128.D0/3.D0*t350+64.D0/3.D0*t352)*t287+(-64
     #.D0/3.D0*t366*s2l2+(64.D0/3.D0*t116*t312*s3l2+128.D0/3.D0*t109*t37
     #2*s4l2)*s34r-64.D0/3.D0*t253*s4l2)*s23r-64.D0/3.D0*t266*t382*s2l2-
     #64.D0/3.D0*t266*t270*t313-64.D0/3.D0*t15*t231)*s1vr)*s1l1
      s5 = s6+s7
      s2 = s4+s5
      s3 = s2+(((-64.D0/3.D0*t88*t242-64.D0/3.D0*t398-128.D0/3.D0*t400)*
     #s23r+64.D0/3.D0*t252*t397*s34r+128.D0/3.D0*t238)*t48+((128.D0/3.D0
     #*t398+64.D0/3.D0*t400)*t287+(-64.D0/3.D0*t366*s2l1+(64.D0/3.D0*t11
     #6*t312*s3l1+128.D0/3.D0*t109*t372*s4l1)*s34r-64.D0/3.D0*t253*s4l1)
     #*s23r-64.D0/3.D0*t266*t382*s2l1-64.D0/3.D0*t266*t397*t313-64.D0/3.
     #D0*t15*t245)*s1vr)*s1l2
      taux1 = s3+(-128.D0/3.D0*t164*t257+128.D0/3.D0*t88*t441+128.D0/3.D
     #0*t47*t366)*t165+((-64.D0/3.D0*t88*t270*s2l1-64.D0/3.D0*t88*t397*s
     #2l2-64.D0/3.D0*t88*t291-64.D0/3.D0*t88*t277+64.D0/3.D0*(1-3*z)*z*t
     #88*t222)*s23r-64.D0/3.D0*t252*t270*s34r*s2l1-64.D0/3.D0*t252*t397*
     #s34r*s2l2+(-128.D0/3.D0*t252*t441-128.D0/3.D0*t116*t47*t222)*s34r+
     #(128.D0/3.D0*t366*s3l2+64.D0/3.D0*t223)*s3l1+64.D0/3.D0*t277+128.D
     #0/3.D0*t1*z*t222)*t48+((-64.D0/3.D0*t366*t172-64.D0/3.D0*t366*t175
     #+128.D0/3.D0*t47*t2*t372*s34r+(128.D0/3.D0*t270+64.D0/3.D0*t495*s4
     #l2)*s3l1+64.D0/3.D0*t495*t188+64.D0/3.D0*t211*t1*t253)*s23r-64.D0/
     #3.D0*t507*t301-64.D0/3.D0*t507*t305+128.D0/3.D0*t266*t253*t313+((-
     #128.D0/3.D0*t116*t270+t274)*s3l1+t279-64.D0/3.D0*t3*t372)*s34r)*s1
     #vr
      t1 = -l3+1
      t2 = -1+z
      t3 = t1*t2
      t6 = 4*t2
      t7 = -t1
      t8 = t6*t7
      t13 = -t2
      t15 = 5*t13*t7
      t38 = s3l2*s2l1
      t40 = s3l1*s2l2
      t42 = s4l2*s3l1
      t44 = s3l2*s4l1
      t54 = t2*l3
      t60 = s1vr**2
      t65 = 8192.D0/3.D0*t13
      t67 = t2**2
      t68 = t7*t67
      t73 = t6*s3l2
      t76 = 8*z-4
      t100 = t6*s3l1
      t108 = t6*l3
      t121 = -t6
      t130 = 3*t13*t7
      t131 = 1/s234r
      t132 = t131*s3l2
      t134 = t131*s4l2
      t135 = t8*t134
      t139 = t2*t7
      t140 = t131*s2l2
      t144 = t7**2
      t145 = t67*t144
      t146 = t131*s34r
      t150 = t145*t132
      t152 = t145*t134
      t160 = t131*s3l1
      t162 = t131*s4l1
      t163 = t8*t162
      t167 = t131*s2l1
      t174 = t145*t160
      t176 = t145*t162
      t192 = t134*s3l1
      t195 = t162*s3l2
      t228 = t54*t132
      t234 = t68*l3
      t242 = 2*t2*l3
      t246 = s23r**2
      t254 = 512.D0/3.D0*t68*l3*t131*s34r-256.D0/3.D0*t6*t131
      t257 = 2*t1*l3
      t258 = t67*t131
      t262 = 3*t7*l3
      t268 = z*t131
      t277 = t67*t2*t144*l3
      t278 = s34r**2
      t279 = t131*t278
      t295 = t54*t160
      t345 = t13*l3
      t353 = t132*s34r*s2l1
      t357 = t160*s34r*s2l2
      t364 = t2*t131
      t377 = (2-z)*t2*t7
      s2 = (((2048.D0/3.D0*t3*s2l2+2048.D0/3.D0*t3*s3l2+2048.D0/3.D0*t8*
     #s4l2)*s1vr+(-2048.D0/3.D0*t15*s2l2-2048.D0/3.D0*t15*s3l2-2048.D0/3
     #.D0*t15*s4l2)*s23r)*s1l1+((-2048.D0/3.D0*t3*s2l1-2048.D0/3.D0*t3*s
     #3l1-2048.D0/3.D0*t8*s4l1)*s1vr+(2048.D0/3.D0*t15*s2l1+2048.D0/3.D0
     #*t15*s3l1+2048.D0/3.D0*t15*s4l1)*s23r)*s1l2+(-2048.D0/3.D0*t8*t38+
     #2048.D0/3.D0*t8*t40+2048.D0/3.D0*t8*t42-2048.D0/3.D0*t8*t44)*s1vr)
     #*s13r+((-8192.D0/3.D0*t2*s2l2*l3-8192.D0/3.D0*t54*s3l2-16384.D0/3.
     #D0*t54*s4l2)*t60+(8192.D0/3.D0*t54*s4l2*s23r+t65*s2l2+8192.D0/3.D0
     #*t68*l3*s4l2*s34r-2048.D0/3.D0*t73-2048.D0/3.D0*t76*s4l2)*s1vr)*s1
     #l1
      s1 = s2+((8192.D0/3.D0*t54*s2l1+8192.D0/3.D0*t54*s3l1+16384.D0/3.D
     #0*t54*s4l1)*t60+(-8192.D0/3.D0*t54*s4l1*s23r-t65*s2l1-8192.D0/3.D0
     #*t68*l3*s4l1*s34r+2048.D0/3.D0*t100+2048.D0/3.D0*t76*s4l1)*s1vr)*s
     #1l2+(2048.D0/3.D0*t108*t38-2048.D0/3.D0*t108*t40-8192.D0/3.D0*t54*
     #t42+8192.D0/3.D0*t54*t44)*t60+(2048.D0/3.D0*t73*s2l1-2048.D0/3.D0*
     #t100*s2l2+2048.D0/3.D0*t121*s4l2*s3l1-2048.D0/3.D0*t121*s4l1*s3l2)
     #*s1vr
      s3 = s1
      s5 = (((-256.D0/3.D0*t130*t132+256.D0/3.D0*t135)*t60+(256.D0/3.D0*
     #t139*t140*s23r-256*t145*t146*s2l2+(-256*t150-1024.D0/3.D0*t152)*s3
     #4r)*s1vr)*s1l1+((256.D0/3.D0*t130*t160-256.D0/3.D0*t163)*t60+(-256
     #.D0/3.D0*t139*t167*s23r+256*t145*t146*s2l1+(256*t174+1024.D0/3.D0*
     #t176)*s34r)*s1vr)*s1l2+((-256.D0/3.D0*t139*t132-256.D0/3.D0*t135)*
     #s2l1+(256.D0/3.D0*t139*t160+256.D0/3.D0*t163)*s2l2-256.D0/3.D0*t3*
     #t192+256.D0/3.D0*t3*t195)*t60+((256.D0/3.D0*t139*t134*s2l1-256.D0/
     #3.D0*t139*t162*s2l2+256.D0/3.D0*t139*t192-256.D0/3.D0*t139*t195)*s
     #23r+(512.D0/3.D0*t150+512.D0/3.D0*t152)*s34r*s2l1+(-512.D0/3.D0*t1
     #74-512.D0/3.D0*t176)*s34r*s2l2+(-256.D0/3.D0*t145*t192+256.D0/3.D0
     #*t145*t195)*s34r)*s1vr)*s13r
      s6 = (((-256.D0/3.D0*t108*t140-1024.D0/3.D0*t228-512.D0/3.D0*t108*
     #t134)*s23r-256.D0/3.D0*t234*t132*s34r+2048.D0/3.D0*t134)*t60+((256
     #.D0/3.D0*t228+256.D0/3.D0*t242*t134)*t246+(t254*s2l2+(-256.D0/3.D0
     #*t257*t258*s3l2+256.D0/3.D0*t262*t258*s4l2)*s34r-256*t268*s3l2-102
     #4.D0/3.D0*t268*s4l2)*s23r+256.D0/3.D0*t277*t279*s2l2+256.D0/3.D0*t
     #277*t132*t278-256.D0/3.D0*t8*t134*s34r)*s1vr)*s1l1
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+(((256.D0/3.D0*t108*t167+1024.D0/3.D0*t295+512.D0/3.D0*
     #t108*t162)*s23r+256.D0/3.D0*t234*t160*s34r-2048.D0/3.D0*t162)*t60+
     #((-256.D0/3.D0*t295-256.D0/3.D0*t242*t162)*t246+(-t254*s2l1+(256.D
     #0/3.D0*t257*t258*s3l1-256.D0/3.D0*t262*t258*s4l1)*s34r+256*t268*s3
     #l1+1024.D0/3.D0*t268*s4l1)*s23r-256.D0/3.D0*t277*t279*s2l1-256.D0/
     #3.D0*t277*t160*t278+256.D0/3.D0*t8*t162*s34r)*s1vr)*s1l2+((256.D0/
     #3.D0*t54*t132*s2l1-256.D0/3.D0*t54*t160*s2l2+512.D0/3.D0*t345*t192
     #-512.D0/3.D0*t345*t195)*s23r-256.D0/3.D0*t234*t353+256.D0/3.D0*t23
     #4*t357+1024.D0/3.D0*t192-1024.D0/3.D0*t195)*t60+((256.D0/3.D0*t364
     #*t38-256.D0/3.D0*t364*t40-512.D0/3.D0*t268*t42+512.D0/3.D0*t268*t4
     #4)*s23r+256.D0/3.D0*t377*t353-256.D0/3.D0*t377*t357+(256.D0/3.D0*t
     #130*t192-256.D0/3.D0*t130*t195)*s34r)*s1vr

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux2
      tcoff6 = taux2
      tcoff7 = -taux2
      tcoff8 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet2 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet3 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet4 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet6 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet7 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet8 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      sq3num1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function sq3num1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = z*t2
      t4 = l3-1
      t5 = t4**2
      t7 = s13r**2
      t11 = t4*t1
      t12 = t11*s2l2
      t13 = t11*s3l2
      t17 = t11*s4l2
      t23 = t11*s2l1
      t24 = t11*s3l1
      t37 = z*l3
      t38 = s1vr**2
      t42 = s3l2*s2l1
      t45 = s3l1*s2l2
      t48 = t2*t1
      t49 = t48*t5
      t50 = t37*s34r
      t57 = s3l2*s4l1
      t66 = z**2
      t67 = t66*t1
      t73 = l3*t1
      t76 = l3*t2
      t77 = t4*s34r
      t80 = t4*t2
      t99 = l3**2
      t101 = t38*s1vr
      t112 = z*t4
      t117 = s3l2*s3l1
      t121 = 2*z-1
      t130 = t1*s3l2
      t149 = 1/s234r
      t154 = t5*t4
      t156 = t149*s34r
      t162 = t37*t2
      t163 = t4*t149
      t167 = t76*t4
      t168 = z*t149
      t172 = t149*s3l2
      t175 = t149*s4l2
      t177 = 64.D0/3.D0*t11*t175
      t178 = -128.D0/3.D0*t11*t172-t177
      t180 = t149*s4l1
      t181 = t180*s3l2
      t183 = 64.D0/3.D0*t11*t181
      t186 = t11*t149
      t191 = s23r**2
      t195 = t172*s2l1
      t198 = t149*s3l1
      t199 = t198*s2l2
      t202 = t49*l3
      t208 = 1+z
      t214 = t2*t5
      t215 = t214*t149
      t217 = s3l2*s34r*s2l1
      t221 = s3l1*s34r*s2l2
      t224 = t2**2
      t227 = l3*t149
      t228 = s34r**2
      t250 = s34r*t38
      t269 = t149*t228
      t308 = t99*t149
      t311 = t172*s3l1
      t312 = t73*t311
      t346 = t80*t149
      t351 = t1*t149
      s2 = 1024.D0/3.D0*t3*t5*s1vr*t7+(((512.D0/3.D0*t12+512.D0/3.D0*t13
     #)*s1vr+(-512.D0/3.D0*t12-512.D0/3.D0*t13-512.D0/3.D0*t17)*s23r)*s1
     #l1+((512.D0/3.D0*t23+512.D0/3.D0*t24)*s1vr+(-512.D0/3.D0*t23-512.D
     #0/3.D0*t24-512.D0/3.D0*t11*s4l1)*s23r)*s1l2-512*t4*t2*t37*t38+(512
     #.D0/3.D0*t11*t42+512.D0/3.D0*t11*t45+512.D0/3.D0*t49*t50+(1024.D0/
     #3.D0*t13+512.D0/3.D0*t17)*s3l1+512.D0/3.D0*t11*t57-1024.D0/3.D0*t1
     #*z*t11)*s1vr-1024.D0/3.D0*t67*t4*s23r)*s13r+(512.D0/3.D0*t73*s4l2*
     #s23r-512.D0/3.D0*t76*t77*s2l2-512.D0/3.D0*t80*l3*s3l2*s34r)*s1vr*s
     #1l1
      s1 = s2+(512.D0/3.D0*t73*s4l1*s23r-512.D0/3.D0*t76*t77*s2l1-512.D0
     #/3.D0*t80*l3*s3l1*s34r)*s1vr*s1l2+1024.D0/3.D0*t2*t99*z*t101+(-512
     #.D0/3.D0*t3*t99*s23r-512.D0/3.D0*t73*t42-512.D0/3.D0*t73*t45-1024.
     #D0/3.D0*t112*t48*t99*s34r-1024.D0/3.D0*t73*t117+1024.D0/3.D0*t121*
     #t1*t37)*t38+(-512.D0/3.D0*t3*l3*s23r-512.D0/3.D0*t130*s2l1-512.D0/
     #3.D0*t1*s3l1*s2l2-512.D0/3.D0*t121*t4*t2*t50+(-1024.D0/3.D0*t130+5
     #12.D0/3.D0*s4l2)*s3l1+512.D0/3.D0*t57+1024.D0/3.D0*t67)*s1vr
      s2 = s1+(64.D0/3.D0*t3*t5*t149*t38-64.D0/3.D0*z*t48*t154*t156*s1vr
     #)*t7+(-128.D0/3.D0*t162*t163*t101+(-128.D0/3.D0*t167*t168*s23r+t17
     #8*s3l1-t183-64.D0/3.D0*(-3+z)*z*t186)*t38+(128.D0/3.D0*t162*t163*t
     #191+(64.D0/3.D0*t11*t195+64.D0/3.D0*t11*t199+320.D0/3.D0*t202*t168
     #*s34r-t178*s3l1+t183-64.D0/3.D0*t208*z*t186)*s23r+64.D0/3.D0*t215*
     #t217+64.D0/3.D0*t215*t221+128.D0/3.D0*z*t224*t154*t227*t228+((128.
     #D0/3.D0*t214*t172+64.D0/3.D0*t214*t175)*s3l1+64.D0/3.D0*t214*t181+
     #64.D0/3.D0*(-2+z)*z*t215)*s34r)*s1vr)*s13r+(64.D0/3.D0*t167*t172*t
     #250+(64.D0/3.D0*t73*t175*t191+(-64.D0/3.D0*t167*t156*s2l2+(-64.D0/
     #3.D0*t76*t163*s3l2+64.D0/3.D0*t76*t163*s4l2)*s34r)*s23r-64.D0/3.D0
     #*t202*t269*s2l2-64.D0/3.D0*t202*t172*t228)*s1vr)*s1l1
      s3 = s2+(64.D0/3.D0*t167*t198*t250+(64.D0/3.D0*t73*t180*t191+(-64.
     #D0/3.D0*t167*t156*s2l1+(-64.D0/3.D0*t76*t163*s3l1+64.D0/3.D0*t76*t
     #163*s4l1)*s34r)*s23r-64.D0/3.D0*t202*t269*s2l1-64.D0/3.D0*t202*t19
     #8*t228)*s1vr)*s1l2
      taux1 = s3+(128.D0/3.D0*t3*t308*s23r+128.D0/3.D0*t312-128.D0/3.D0*
     #t73*t168)*t101+(-64*t3*t308*t191+(-64.D0/3.D0*t73*t195-64.D0/3.D0*
     #t73*t199-256.D0/3.D0*t112*t48*t308*s34r-128.D0/3.D0*t312+64.D0/3.D
     #0*(-1+3*z)*z*t73*t149)*s23r-64.D0/3.D0*t167*t172*s34r*s2l1-64.D0/3
     #.D0*t167*t198*s34r*s2l2+(-128.D0/3.D0*t167*t311+128.D0/3.D0*t37*t3
     #46)*s34r+128.D0/3.D0*t351*t117-128.D0/3.D0*t1*z*t149)*t38+(-64*t3*
     #t227*t191+(-64.D0/3.D0*t351*t42-64.D0/3.D0*t351*t45-64.D0/3.D0*(-5
     #+4*z)*z*l3*t80*t156+(-128.D0/3.D0*t351*s3l2+64.D0/3.D0*t175)*s3l1+
     #64.D0/3.D0*t181+64.D0/3.D0*t208*t1*t168)*s23r-64.D0/3.D0*t346*t217
     #-64.D0/3.D0*t346*t221+128.D0/3.D0*t202*t168*t228+((-128.D0/3.D0*t8
     #0*t172+t177)*s3l1+t183+64.D0/3.D0*t3*t163)*s34r)*s1vr
      t1 = l3-1
      t3 = -1+z
      t4 = 4*t3*t1
      t7 = -t3
      t8 = 5*t7
      t9 = t8*t1
      t15 = 9*t3*t1
      t38 = 4*t7
      t39 = t38*t1
      t40 = s3l2*s2l1
      t42 = s3l1*s2l2
      t44 = s4l2*s3l1
      t46 = s3l2*s4l1
      t56 = t3*l3
      t62 = s1vr**2
      t67 = t3**2
      t68 = l3*t67
      t73 = -2048.D0/3.D0*t68*t1*s34r-10240.D0/3.D0+10240.D0/3.D0*z
      t75 = t1*t67
      t87 = -10*z+5
      t123 = t38*l3
      t139 = -4+5*z
      t147 = t3*t1
      t148 = 1/s234r
      t149 = t148*s3l2
      t151 = 512.D0/3.D0*t147*t149
      t152 = t148*s4l2
      t154 = 256.D0/3.D0*t9*t152
      t157 = -t1
      t159 = 3*t157*t3
      t160 = t148*s2l2
      t164 = t1**2
      t165 = t67*t164
      t166 = t148*s34r
      t171 = 512.D0/3.D0*t165*t149
      t172 = t165*t152
      t180 = t148*s3l1
      t182 = 512.D0/3.D0*t147*t180
      t183 = t148*s4l1
      t185 = 256.D0/3.D0*t9*t183
      t188 = t148*s2l1
      t196 = 512.D0/3.D0*t165*t180
      t197 = t165*t183
      t223 = t152*s3l1
      t226 = t183*s3l2
      t243 = t8*l3
      t254 = 2*t7
      t255 = t254*l3
      t259 = s23r**2
      t267 = -256*t68*t1*t148*s34r-256.D0/3.D0*t8*t148
      t269 = t157*t67
      t270 = l3*t148
      t274 = 5*t157*l3
      t275 = t67*t148
      t281 = z*t148
      t347 = t38*t148
      t352 = (-1+4*z)*t148
      t360 = -t254*t1*t148
      s2 = (((2048.D0/3.D0*t4*s2l2+2048.D0/3.D0*t4*s3l2+2048.D0/3.D0*t9*
     #s4l2)*s1vr+(-2048.D0/3.D0*t15*s2l2-2048.D0/3.D0*t15*s3l2-2048.D0/3
     #.D0*t15*s4l2)*s23r)*s1l1+((-2048.D0/3.D0*t4*s2l1-2048.D0/3.D0*t4*s
     #3l1-2048.D0/3.D0*t9*s4l1)*s1vr+(2048.D0/3.D0*t15*s2l1+2048.D0/3.D0
     #*t15*s3l1+2048.D0/3.D0*t15*s4l1)*s23r)*s1l2+(-2048.D0/3.D0*t39*t40
     #+2048.D0/3.D0*t39*t42+2048.D0/3.D0*t39*t44-2048.D0/3.D0*t39*t46)*s
     #1vr)*s13r+((10240.D0/3.D0*t3*s2l2*l3+10240.D0/3.D0*t56*s3l2+20480.
     #D0/3.D0*t56*s4l2)*t62+(-8192.D0/3.D0*t56*s4l2*s23r+t73*s2l2+(-2048
     #.D0/3.D0*t75*l3*s3l2-10240.D0/3.D0*t75*l3*s4l2)*s34r-2048.D0/3.D0*
     #t8*s3l2-2048.D0/3.D0*t87*s4l2)*s1vr)*s1l1
      s1 = s2+((-10240.D0/3.D0*t56*s2l1-10240.D0/3.D0*t56*s3l1-20480.D0/
     #3.D0*t56*s4l1)*t62+(8192.D0/3.D0*t56*s4l1*s23r-t73*s2l1+(2048.D0/3
     #.D0*t75*l3*s3l1+10240.D0/3.D0*t75*l3*s4l1)*s34r+2048.D0/3.D0*t8*s3
     #l1+2048.D0/3.D0*t87*s4l1)*s1vr)*s1l2+(2048.D0/3.D0*t123*t40-2048.D
     #0/3.D0*t123*t42+10240.D0/3.D0*t56*t44-10240.D0/3.D0*t56*t46)*t62+(
     #2048.D0/3.D0*t38*s3l2*s2l1-2048.D0/3.D0*t38*s3l1*s2l2+2048.D0/3.D0
     #*t139*s4l2*s3l1-2048.D0/3.D0*t139*s4l1*s3l2)*s1vr
      s3 = s1
      s5 = (((-t151+t154)*t62+(256.D0/3.D0*t159*t160*s23r+512.D0/3.D0*t1
     #65*t166*s2l2+(t171+1280.D0/3.D0*t172)*s34r)*s1vr)*s1l1+((t182-t185
     #)*t62+(-256.D0/3.D0*t159*t188*s23r-512.D0/3.D0*t165*t166*s2l1+(-t1
     #96-1280.D0/3.D0*t197)*s34r)*s1vr)*s1l2+((-t151-t154)*s2l1+(t182+t1
     #85)*s2l2)*t62+(((-256.D0/3.D0*t9*t149+512.D0/3.D0*t147*t152)*s2l1+
     #(256.D0/3.D0*t9*t180-512.D0/3.D0*t147*t183)*s2l2+256.D0/3.D0*t159*
     #t223-256.D0/3.D0*t159*t226)*s23r+(t171-256.D0/3.D0*t172)*s34r*s2l1
     #+(-t196+256.D0/3.D0*t197)*s34r*s2l2)*s1vr)*s13r
      s6 = (((-256.D0/3.D0*t243*t160-512.D0/3.D0*t243*t152)*s23r-2560.D0
     #/3.D0*t152)*t62+((256.D0/3.D0*t56*t149+256.D0/3.D0*t255*t152)*t259
     #+(t267*s2l2+(-256.D0/3.D0*t269*t270*s3l2+256.D0/3.D0*t274*t275*s4l
     #2)*s34r+512.D0/3.D0*t281*s3l2+1280.D0/3.D0*t281*s4l2)*s23r-256.D0/
     #3.D0*t9*t152*s34r)*s1vr)*s1l1
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+(((256.D0/3.D0*t243*t188+512.D0/3.D0*t243*t183)*s23r+25
     #60.D0/3.D0*t183)*t62+((-256.D0/3.D0*t56*t180-256.D0/3.D0*t255*t183
     #)*t259+(-t267*s2l1+(256.D0/3.D0*t269*t270*s3l1-256.D0/3.D0*t274*t2
     #75*s4l1)*s34r-512.D0/3.D0*t281*s3l1-1280.D0/3.D0*t281*s4l1)*s23r+2
     #56.D0/3.D0*t9*t183*s34r)*s1vr)*s1l2+((256.D0/3.D0*t123*t149*s2l1-2
     #56.D0/3.D0*t123*t180*s2l2+1024.D0/3.D0*t56*t223-1024.D0/3.D0*t56*t
     #226)*s23r-1280.D0/3.D0*t223+1280.D0/3.D0*t226)*t62+((256.D0/3.D0*t
     #347*t40-256.D0/3.D0*t347*t42+256.D0/3.D0*t352*t44-256.D0/3.D0*t352
     #*t46)*s23r+256.D0/3.D0*t360*s3l2*s34r*s2l1-256.D0/3.D0*t360*s3l1*s
     #34r*s2l2+(256.D0/3.D0*t147*t223-256.D0/3.D0*t147*t226)*s34r)*s1vr

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux2
      tcoff6 = taux2
      tcoff7 = -taux2
      tcoff8 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet2 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet3 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet4 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet6 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet7 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet8 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      sq3num1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function sq3num1e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = z*t2
      t4 = l3-1
      t5 = t4**2
      t7 = s13r**2
      t11 = t4*t2
      t12 = l3*z
      t13 = s1vr**2
      t19 = t2*t1
      t28 = z**2
      t35 = l3**2
      t47 = l3*t2
      t49 = 1/s234r
      t53 = t12*t19
      t54 = t5*t49
      t60 = s34r*s23r
      t63 = t2**2
      t67 = l3*t49
      t68 = s34r**2
      t69 = t67*t68
      t76 = t35*t49
      t77 = s23r**2
      taux1 = -512.D0/3.D0*t3*t5*s1vr*t7+(512.D0/3.D0*t11*t12*t13+(512.D
     #0/3.D0*t11*t12*s23r+512.D0/3.D0*t19*t5*t12*s34r-512.D0/3.D0*z*t1*t
     #4)*s1vr+512.D0/3.D0*t28*t1*t4*s23r)*s13r-512.D0/3.D0*t3*t35*s23r*t
     #13+(-512.D0/3.D0*t3*l3*s23r+512.D0/3.D0*t12*t11*s34r)*s1vr+((64.D0
     #/3.D0*t47*t4*z*t49*s23r-64.D0/3.D0*t53*t54*s34r)*t13+(64.D0/3.D0*t
     #53*t54*t60+64.D0/3.D0*z*t63*t5*t4*t69)*s1vr)*s13r+(-64.D0/3.D0*t3*
     #t76*t77-64.D0/3.D0*z*t4*t19*t76*t60)*t13+(-64.D0/3.D0*t3*t67*t77-6
     #4.D0/3.D0*(-2+z)*z*t47*t4*t49*t60+64.D0/3.D0*z*t19*t5*t69)*s1vr
      t1 = -l3+1
      t2 = 5*t1
      t3 = -1+z
      t4 = t2*t3
      t9 = -t1
      t10 = t9*t3
      t15 = -t3
      t17 = 7*t15*t9
      t42 = t15*t9
      t43 = s3l2*s2l1
      t45 = s3l1*s2l2
      t47 = s4l2*s3l1
      t49 = s3l2*s4l1
      t59 = l3*t3
      t65 = s1vr**2
      t70 = t3**2
      t76 = 2048*l3*t70*t9*s34r+4096.D0/3.D0-4096.D0/3.D0*z
      t78 = t70*t9
      t89 = 2*z
      t90 = t89-1
      t126 = t15*l3
      t141 = -1-t89
      t149 = 1/s234r
      t150 = t149*s4l2
      t154 = t149*s2l2
      t157 = t9**2
      t158 = t70*t157
      t159 = t150*s34r
      t166 = t149*s4l1
      t170 = t149*s2l1
      t173 = t166*s34r
      t187 = 3*t9*t3
      t188 = t150*s3l1
      t191 = t166*s3l2
      t196 = t149*s3l2
      t197 = t196*s2l1
      t199 = t149*s3l1
      t200 = t199*s2l2
      t202 = t10*t188
      t203 = t10*t191
      t236 = t78*l3
      t245 = 3*t3
      t246 = t245*l3
      t250 = s23r**2
      t252 = l3*t149
      t258 = -256.D0/3.D0*t78*t252*s34r-512.D0/3.D0*t3*t149
      t263 = -l3*t2
      t264 = t70*t149
      t270 = z*t149
      t277 = t70*t3*t157*l3
      t278 = s34r**2
      t279 = t149*t278
      t341 = t196*s34r*s2l1
      t345 = t199*s34r*s2l2
      t352 = t245*t149
      t359 = t245*t3*t9
      s2 = (((2048.D0/3.D0*t4*s2l2+2048.D0/3.D0*t4*s3l2+4096.D0/3.D0*t10
     #*s4l2)*s1vr+(-2048.D0/3.D0*t17*s2l2-2048.D0/3.D0*t17*s3l2-2048.D0/
     #3.D0*t17*s4l2)*s23r)*s1l1+((-2048.D0/3.D0*t4*s2l1-2048.D0/3.D0*t4*
     #s3l1-4096.D0/3.D0*t10*s4l1)*s1vr+(2048.D0/3.D0*t17*s2l1+2048.D0/3.
     #D0*t17*s3l1+2048.D0/3.D0*t17*s4l1)*s23r)*s1l2+(-2048.D0/3.D0*t42*t
     #43+2048.D0/3.D0*t42*t45+2048.D0/3.D0*t42*t47-2048.D0/3.D0*t42*t49)
     #*s1vr)*s13r+((-4096.D0/3.D0*t3*s2l2*l3-4096.D0/3.D0*t59*s3l2-8192.
     #D0/3.D0*t59*s4l2)*t65+(-2048.D0/3.D0*t59*s4l2*s23r+t76*s2l2+(2048*
     #t78*l3*s3l2+4096.D0/3.D0*t78*l3*s4l2)*s34r-4096.D0/3.D0*t3*s3l2-40
     #96.D0/3.D0*t90*s4l2)*s1vr)*s1l1
      s1 = s2+((4096.D0/3.D0*t59*s2l1+4096.D0/3.D0*t59*s3l1+8192.D0/3.D0
     #*t59*s4l1)*t65+(2048.D0/3.D0*t59*s4l1*s23r-t76*s2l1+(-2048*t78*l3*
     #s3l1-4096.D0/3.D0*t78*l3*s4l1)*s34r+4096.D0/3.D0*t3*s3l1+4096.D0/3
     #.D0*t90*s4l1)*s1vr)*s1l2+(2048.D0/3.D0*t126*t43-2048.D0/3.D0*t126*
     #t45-4096.D0/3.D0*t59*t47+4096.D0/3.D0*t59*t49)*t65+(2048.D0/3.D0*t
     #15*s3l2*s2l1-2048.D0/3.D0*t15*s3l1*s2l2+2048.D0/3.D0*t141*s4l2*s3l
     #1-2048.D0/3.D0*t141*s4l1*s3l2)*s1vr
      s3 = s1
      s5 = ((512.D0/3.D0*t10*t150*t65+(512.D0/3.D0*t10*t154*s23r-512.D0/
     #3.D0*t158*t159)*s1vr)*s1l1+(-512.D0/3.D0*t10*t166*t65+(-512.D0/3.D
     #0*t10*t170*s23r+512.D0/3.D0*t158*t173)*s1vr)*s1l2+(-512.D0/3.D0*t1
     #0*t150*s2l1+512.D0/3.D0*t10*t166*s2l2-256.D0/3.D0*t187*t188+256.D0
     #/3.D0*t187*t191)*t65+((-256.D0/3.D0*t10*t197+256.D0/3.D0*t10*t200+
     #256.D0/3.D0*t202-256.D0/3.D0*t203)*s23r+(-256*t158*t196+512.D0/3.D
     #0*t158*t150)*s34r*s2l1+(256*t158*t199-512.D0/3.D0*t158*t166)*s34r*
     #s2l2+(256*t158*t188-256*t158*t191)*s34r)*s1vr)*s13r
      s6 = (((-512.D0/3.D0*t59*t154-1024.D0/3.D0*t59*t150)*s23r+256*t236
     #*t196*s34r+1024.D0/3.D0*t150)*t65+((512.D0/3.D0*t59*t196+256.D0/3.
     #D0*t246*t150)*t250+(t258*s2l2+(-256.D0/3.D0*t78*t252*s3l2+256.D0/3
     #.D0*t263*t264*s4l2)*s34r-512.D0/3.D0*t270*s4l2)*s23r-256*t277*t279
     #*s2l2-256*t277*t196*t278-512.D0/3.D0*t10*t159)*s1vr)*s1l1
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+(((512.D0/3.D0*t59*t170+1024.D0/3.D0*t59*t166)*s23r-256
     #*t236*t199*s34r-1024.D0/3.D0*t166)*t65+((-512.D0/3.D0*t59*t199-256
     #.D0/3.D0*t246*t166)*t250+(-t258*s2l1+(256.D0/3.D0*t78*t252*s3l1-25
     #6.D0/3.D0*t263*t264*s4l1)*s34r+512.D0/3.D0*t270*s4l1)*s23r+256*t27
     #7*t279*s2l1+256*t277*t199*t278+512.D0/3.D0*t10*t173)*s1vr)*s1l2+((
     #256.D0/3.D0*t246*t197-256.D0/3.D0*t246*t200)*s23r+256*t236*t341-25
     #6*t236*t345+512.D0/3.D0*t188-512.D0/3.D0*t191)*t65+((256.D0/3.D0*t
     #352*t43-256.D0/3.D0*t352*t45+256.D0/3.D0*t188-256.D0/3.D0*t191)*s2
     #3r+256.D0/3.D0*t359*t341-256.D0/3.D0*t359*t345+(256.D0/3.D0*t202-2
     #56.D0/3.D0*t203)*s34r)*s1vr

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux2
      tcoff6 = taux2
      tcoff7 = -taux2
      tcoff8 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet2 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet3 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet4 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet6 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet7 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet8 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      sq3num1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function sq3num1e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux1 = 0
      t1 = l3-1
      t2 = -1+z
      t3 = t1*t2
      t4 = t3*s2l2
      t5 = t3*s3l2
      t15 = t3*s2l1
      t16 = t3*s3l1
      t26 = s3l2*s2l1
      t28 = s3l1*s2l2
      t30 = s4l2*s3l1
      t32 = s3l2*s4l1
      t39 = l3*t2
      t42 = t2**2
      t43 = l3*t42
      t44 = t1*s34r
      t47 = t42*t1
      t70 = s1vr**2
      t79 = 1/s234r
      t80 = t79*s4l2
      t81 = t80*s3l1
      t82 = t3*t81
      t83 = t79*s4l1
      t84 = t83*s3l2
      t85 = t3*t84
      t87 = -512.D0/3.D0*t82+512.D0/3.D0*t85
      t89 = t79*s3l2
      t90 = t89*s2l1
      t92 = t79*s3l1
      t93 = t92*s2l2
      t98 = t1**2
      t99 = t42*t98
      t100 = t99*t79
      t102 = s3l2*s34r*s2l1
      t106 = s3l1*s34r*s2l2
      t118 = t43*t1
      t119 = s34r*t70
      t123 = s23r**2
      t127 = t79*s34r
      t131 = l3*t1
      t132 = t42*t79
      t144 = t42*t2*t98*l3
      t145 = s34r**2
      t146 = t79*t145
      t200 = t2*t79
      t206 = t47*t79
      s1 = (((4096.D0/3.D0*t4+4096.D0/3.D0*t5)*s1vr+(-4096.D0/3.D0*t4-40
     #96.D0/3.D0*t5-4096.D0/3.D0*t3*s4l2)*s23r)*s1l1+((-4096.D0/3.D0*t15
     #-4096.D0/3.D0*t16)*s1vr+(4096.D0/3.D0*t15+4096.D0/3.D0*t16+4096.D0
     #/3.D0*t3*s4l1)*s23r)*s1l2+(-4096.D0/3.D0*t3*t26+4096.D0/3.D0*t3*t2
     #8+4096.D0/3.D0*t3*t30-4096.D0/3.D0*t3*t32)*s1vr)*s13r+(4096.D0/3.D
     #0*t39*s4l2*s23r-4096.D0/3.D0*t43*t44*s2l2-4096.D0/3.D0*t47*l3*s3l2
     #*s34r)*s1vr*s1l1+(-4096.D0/3.D0*t39*s4l1*s23r+4096.D0/3.D0*t43*t44
     #*s2l1+4096.D0/3.D0*t47*l3*s3l1*s34r)*s1vr*s1l2+(4096.D0/3.D0*t39*t
     #26-4096.D0/3.D0*t39*t28)*t70+(4096.D0/3.D0*t2*s3l2*s2l1-4096.D0/3.
     #D0*t2*s3l1*s2l2+4096.D0/3.D0*t30-4096.D0/3.D0*t32)*s1vr
      s2 = s1+(t87*t70+((-512.D0/3.D0*t3*t90+512.D0/3.D0*t3*t93+512.D0/3
     #.D0*t82-512.D0/3.D0*t85)*s23r-512.D0/3.D0*t100*t102+512.D0/3.D0*t1
     #00*t106+(512.D0/3.D0*t99*t81-512.D0/3.D0*t99*t84)*s34r)*s1vr)*s13r
     #+(512.D0/3.D0*t118*t89*t119+(512.D0/3.D0*t39*t80*t123+(-512.D0/3.D
     #0*t118*t127*s2l2+(-512.D0/3.D0*t131*t132*s3l2+512.D0/3.D0*t131*t13
     #2*s4l2)*s34r)*s23r-512.D0/3.D0*t144*t146*s2l2-512.D0/3.D0*t144*t89
     #*t145)*s1vr)*s1l1
      taux2 = s2+(-512.D0/3.D0*t118*t92*t119+(-512.D0/3.D0*t39*t83*t123+
     #(512.D0/3.D0*t118*t127*s2l1+(512.D0/3.D0*t131*t132*s3l1-512.D0/3.D
     #0*t131*t132*s4l1)*s34r)*s23r+512.D0/3.D0*t144*t146*s2l1+512.D0/3.D
     #0*t144*t92*t145)*s1vr)*s1l2+((512.D0/3.D0*t39*t90-512.D0/3.D0*t39*
     #t93)*s23r+512.D0/3.D0*t118*t89*s34r*s2l1-512.D0/3.D0*t118*t92*s34r
     #*s2l2)*t70+((512.D0/3.D0*t200*t26-512.D0/3.D0*t200*t28+512.D0/3.D0
     #*t81-512.D0/3.D0*t84)*s23r+512.D0/3.D0*t206*t102-512.D0/3.D0*t206*
     #t106-t87*s34r)*s1vr

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux2
      tcoff6 = taux2
      tcoff7 = -taux2
      tcoff8 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet2 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet3 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet4 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet6 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet7 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet8 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      sq3num1e4 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function sq3num2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3**2
      t2 = l3-1
      t3 = t1*t2
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t7 = t6*s3l2
      t9 = l1**2
      t10 = r**2
      t12 = 1+t10
      t13 = t12**2
      t14 = r+t
      t15 = t14**2
      t19 = 1+r*t
      t20 = t19**2
      t22 = t13/t15/t20
      t23 = t9*t10*t22
      t26 = t1*t5
      t28 = t26*s3l2
      t35 = 2*l3
      t36 = t35-1
      t37 = t36*t5
      t38 = l3*s3l2
      t39 = t37*t38
      t47 = t12/t14/t19
      t49 = t4*s2l2
      t52 = t26*s3l2*s34r
      t53 = t4*l3
      t54 = t53*s3l2
      t58 = l3*t5
      t62 = s34r**2
      t69 = 512.D0/3.D0*t54
      t72 = 1024.D0/3.D0*t54
      t77 = t53*s4l2
      t78 = 512.D0/3.D0*t77
      t79 = 1024.D0/3.D0*t39-t78
      t81 = s3l2*s4l1
      t83 = 512.D0/3.D0*t53*t81
      t84 = z*l3
      t92 = s4l1*s2l2
      t93 = s4l2*s2l1
      t94 = t38*s2l1
      t95 = l3*s3l1
      t96 = t95*s2l2
      t97 = l3*s2l2
      t98 = t97*s2l1
      t99 = 2*t98
      t100 = s1l1*l3
      t101 = t100*s2l2
      t102 = l3*s2l1
      t104 = t102*s2l2*z
      t105 = 2*t104
      t106 = z*s3l2
      t107 = t102*t106
      t108 = z*s3l1
      t109 = t97*t108
      t110 = s1l1*s2l2
      t111 = t110*t84
      t112 = s1l2*z
      t113 = t112*t102
      t114 = s2l2*s2l1
      t115 = t114*z
      t116 = 2*t115
      t117 = s3l2*s2l1
      t118 = t117*z
      t119 = s3l1*s2l2
      t120 = t119*z
      t121 = t100*s3l2
      t122 = t100*t106
      t123 = s1l2*s2l1
      t124 = t123*l3
      t125 = s1l2*l3
      t126 = t125*s3l1
      t127 = t125*t108
      t128 = 2*t114
      t129 = -t92-t93+t94+t96+t99+t101-t105-t107-t109-t111-t113+t116+t11
     #8+t120+t121-t122+t124+t126-t127-t128
      t130 = 1/s134r
      t134 = t3*t7*s3l1
      t135 = 1024.D0/3.D0*t134
      t136 = z*t1
      t137 = t6*t2
      t139 = 512.D0/3.D0*t136*t137
      t142 = s1l1*s1l2
      t144 = 2*t142*z
      t145 = t142*t84
      t146 = 4*t145
      t147 = 2*t142
      t148 = 2*t101
      t149 = 2*t111
      t150 = 2*t113
      t151 = 2*t124
      t152 = s3l1*s1l2
      t153 = t144-t146-t147+t128+t105-t99-t148+t149+t150-t116-t121+t122+
     #t127-t151-t126+t152
      t154 = s3l2*s1l1
      t155 = z**2
      t156 = l3*t155
      t157 = s1l2*s4l1
      t158 = s1l1*s4l2
      t159 = t110*z
      t160 = t112*s2l1
      t161 = t152*z
      t162 = t154*z
      t163 = t142*l3
      t164 = 4*t163
      t165 = z+t154-t117-t155+t156-t84-t119+t123+t157+t110+t158-t159-t16
     #0-t161-t162+t164
      t171 = l3*t1
      t174 = (-t147+t110+t123)*t130
      t177 = -t92-t93-t117+t154+t152+t156-t155-t84-t119+t123+t157+t110+t
     #158-t147+z
      t179 = t5*t130
      t180 = t179*s34r
      t182 = -t92-t93-t117+t154+t152+t156-t155-t84-t119+t157+t158+z
      t184 = t4*t130
      t189 = t2*t171
      t190 = t5**2
      t191 = t190*z
      t199 = t1**2
      t201 = t190*t4
      t203 = z*t130
      t208 = t190*t171
      t213 = t9*l1
      t215 = t10*r
      t223 = t13*t12/t15/t14/t20/t19
      t226 = t26*s4l2
      t228 = -1024.D0/3.D0*t28-512.D0/3.D0*t226
      t230 = t26*s3l1
      t232 = t26*s4l1
      t241 = 512.D0/3.D0*t26*t81
      t242 = 1-l3+t84-z
      t243 = 2*t242
      t245 = t5*z
      t254 = t6*s4l2
      t260 = t6*s4l1
      t263 = t53*s3l1
      s1 = (512.D0/3.D0*t3*t7*t23+((-512.D0/3.D0*t26*s2l2-512.D0/3.D0*t2
     #8)*s23r-1024.D0/3.D0*t3*t7*s34r-512.D0/3.D0*t39)*l1*r*t47+(512.D0/
     #3.D0*t49*l3+512.D0/3.D0*t52+512.D0/3.D0*t54)*s23r-512.D0/3.D0*t58*
     #s34r*s2l2+512.D0/3.D0*t3*t7*t62+512.D0/3.D0*t37*t38*s34r+t69)*s1l1
     #+(t72-512.D0/3.D0*s4l2)*s3l1+(t79*s3l1-t83-512.D0/3.D0*t84*t2*t5)*
     #s34r-512.D0/3.D0*t81-64.D0/3.D0*t129*t130
      t276 = s1+(t135-t139)*t62+64.D0/3.D0*t26*(t153+t165)*t130*t62+(64.
     #D0/3.D0*t6*t171*t174*t62+64.D0/3.D0*t177*t1*t180+64.D0/3.D0*t182*l
     #3*t184)*s23r+(-512.D0/3.D0*t189*t191*s23r+t135)*t9*t10*t22+(-64.D0
     #/3.D0*t2*t199*t201*t203*s23r+64.D0/3.D0*(-t110-t154-t123-t152-t128
     #-t117-t119)*t2*t208*t130)*t213*t215*t223+((t228*s2l1+(-1024.D0/3.D
     #0*t230-512.D0/3.D0*t232)*s2l2+512.D0/3.D0*t189*t191*s34r+t228*s3l1
     #-t241+512.D0/3.D0*t243*t1*t245)*s23r+(1024.D0/3.D0*t3*t6*s34r*s2l2
     #+512.D0/3.D0*t3*t254*s34r+t69)*s2l1+(512.D0/3.D0*t3*t260*s34r+512.
     #D0/3.D0*t263)*s2l2+(-2048.D0/3.D0*t134+t139)*s34r-t79*s3l1+t83)*l1
     #*r*t47
      t291 = t26*s3l1*s34r
      t294 = t53*s4l1
      t310 = -512.D0/3.D0*t242*z*t53
      t316 = t58*t2*s34r
      t330 = t2*z
      t335 = 3*t156
      t336 = 3*t155
      t337 = 3*t84
      t338 = 3*z
      t339 = -t92-t93-t117+t154+t152+t335-t336-t337-t119+t157+t158+t338
      t346 = 2*t154
      t347 = 2*t152
      t348 = 4*t114
      t349 = 2*t123
      t350 = 2*t110
      t354 = t171*t130
      t355 = t354*s34r
      t358 = 6*t114
      t359 = 6*t104
      t360 = 3*t107
      t361 = 3*t94
      t362 = 3*t96
      t363 = 6*t98
      t364 = 3*t101
      t365 = 3*t109
      t366 = 3*t111
      t367 = 3*t113
      t368 = 6*t115
      t369 = 3*t118
      t370 = 3*t120
      t371 = -t92-t93-t358-t359-t360+t361+t362+t363+t364-t365-t366-t367+
     #t368+t369+t370
      t372 = 3*t121
      t373 = 3*t122
      t374 = 3*t127
      t375 = 3*t124
      t376 = 3*t126
      t377 = 2*t117
      t378 = 2*t119
      t379 = 2*t159
      t380 = 2*t160
      t383 = t372-t373-t374+t375+t376-t347-t346-t377-t378-t349-t350+t379
     #+t380+2*t161+2*t162
      t392 = t2**2
      t394 = t201*t171
      t406 = t392*t130
      t407 = t406*s34r
      t410 = 2*t155
      t411 = 2*z
      t412 = -t92-t93-t377+t335-t410-t337-t378-t128-t123+t157-t110+t158+
     #t411
      t425 = -t155+z-t84+t156
      t436 = t3*t6
      t437 = 2*t92
      t438 = 2*t93
      t439 = 3*t117
      t441 = 4*t84
      t442 = 3*t119
      t443 = 2*t157
      t444 = 2*t158
      t445 = -t437-t438-t439+4*t156-t336-t441-t442-t348-t123+t443-t110+t
     #444+t411
      t451 = 2*s4l1*s4l2
      t452 = t1*t155
      t454 = t155*z
      t455 = t1*t454
      t457 = 4*t104
      t461 = -t451-6*t452+3*t455+t454-t81-t437-t438-t348-t457-4*t107+4*t
     #94+4*t96
      t462 = 4*t98
      t464 = 4*t115
      t465 = s4l2*s3l1
      t466 = t157*z
      t467 = t158*z
      t468 = t462+t148-4*t109-t149-t150+t464+t369+t370+t151-t465-t466-t4
     #67
      t470 = t158*l3
      t471 = 2*t470
      t474 = l3*s4l1*s2l2
      t475 = 2*t474
      t477 = l3*s4l2*s2l1
      t478 = 2*t477
      t479 = t157*l3
      t480 = 2*t479
      t483 = -t471+z-4*t117+t475+t478-t480-t410+8*t156-t441-4*t119-t123+
     #t157
      t485 = t102*s4l2*z
      t486 = 2*t485
      t488 = t97*s4l1*z
      t489 = 2*t488
      t490 = t157*t84
      t491 = 2*t490
      t492 = l3*t454
      t495 = t92*z
      t496 = t93*z
      t497 = t158*t84
      t498 = 2*t497
      t499 = -t110+t158+t159+t160-t486-t489+t491-4*t492+3*t136+t495+t496
     #+t498
      t525 = 2*t118
      t526 = 2*t120
      t527 = t361+t362+t462+t101-t365-t111-t113+t464+t525+t526+t124
      t530 = 2*t84
      t531 = -t465-t466-t467-t471-t439+t475+t478-t480-t155+5*t156-t530
      t534 = -t442+t157+t158-t486-t489+t491-3*t492+2*t136+t495+t496+t498
      t542 = 2*t107
      t543 = 2*t94
      t544 = 2*t96
      t545 = 2*t109
      t546 = -t451-2*t452+t455-t81-t92-t93-t128-t105-t542+t543+t544+t99+
     #t101-t545-t111-t113+t116
      t547 = 2*t156
      t548 = t118+t120+t124-t465-t470-t377+t474+t477-t479+t547-t84-t378-
     #t485-t488+t490-t492+t136+t497
      t556 = -2*t145+t92+t93+t348+t457+t107-t94-t96-t462-t148+t109+t149+
     #t150-t464-t118-t120-2*t121
      t561 = 2*t122+2*t127-t151-2*t126+t152+z+t154-t117-t155+t156-t84-t1
     #19+t157+t158-t161-t162+2*t163
      t567 = -t92-t93-t117+t154+t152+t547-t410-t530-t119+t123+t157+t110+
     #t158-t147+t411
      t570 = -t437-t438-t377+t346+t347+t335-t336-t337-t378+t443+t444+t33
     #8
      t576 = t208*t2
      t591 = t144-t146-t147+t92+t93+8*t114+8*t104+t542-t543-t544-8*t98-4
     #*t101+t545+4*t111+4*t113-8*t115-t525-t526-4*t121
      t600 = 4*t122+4*t127-4*t124-4*t126+3*t152+t411+3*t154-t410+t547-t5
     #30+t349+t157+t350+t158-t379-t380-3*t161-3*t162+t164
      t605 = -t437-t438-t358-t359-t360+t361+t362+t363+t364-t365-t366-t36
     #7+t368+t369+t370
      t606 = t372-t373-t374+t375+t376-t152-t154-t117-t119-t123-t110+t159
     #+t160+t161+t162
      t624 = t191*t1
      t627 = l3*t2
      t628 = z-t155-t530+t547
      t645 = t4*t2
      t651 = s13r**2
      t653 = t6*s3l1
      t687 = t58*t2*s2l2
      t689 = l1*r*t47
      t696 = 512.D0/3.D0*t53+512.D0/3.D0*t316
      t710 = t190*t392
      t718 = t2*t5
      t720 = 512.D0/3.D0*t718*t38
      t723 = 1024.D0/3.D0*t58*t2*s4l2
      t731 = -512.D0/3.D0*t718*t95-1024.D0/3.D0*t58*t2*s4l1
      t736 = z+t35-t530-1
      t763 = -t4
      t764 = t763*l3
      t768 = -t763*t36
      s2 = (-512.D0/3.D0*t3*t260*t62+(-512.D0/3.D0*t53*t108-512.D0/3.D0*
     #t232)*s34r-512.D0/3.D0*s3l1)*s2l2+((512.D0/3.D0*t52+t72+t78)*s2l1+
     #(512.D0/3.D0*t291+1024.D0/3.D0*t263+512.D0/3.D0*t294)*s2l2+(-t228*
     #s3l1+t241-512.D0/3.D0*t245*t1*t242)*s34r+(t72+t78)*s3l1+t83+t310)*
     #s23r
      s3 = s2+((-1024.D0/3.D0*t3*t6*t62-1024.D0/3.D0*t316)*s2l2-512.D0/3
     #.D0*t3*t254*t62+(-512.D0/3.D0*t53*t106-512.D0/3.D0*t226)*s34r-512.
     #D0/3.D0*s3l2)*s2l1
      s4 = s3
      s6 = ((64.D0/3.D0*t330*t201*t199*t130*s34r+64.D0/3.D0*t339*t171*t6
     #*t130)*s23r+64.D0/3.D0*(t117+t346+t347+t119+t348+t349+t350-t147+z)
     #*t2*t190*t355-64.D0/3.D0*(t371+t383)*t1*t179)*t9*t10*t22
      s7 = (64.D0/3.D0*t392*z*t394*t130*t213*t215*t223+(64.D0/3.D0*t330*
     #t190*t354*s23r-128.D0/3.D0*t394*z*t407-64.D0/3.D0*t412*t1*t137*t13
     #0)*t9*t10*t22+((-64.D0/3.D0*t208*z*t2*t130*s34r-128.D0/3.D0*t425*t
     #1*t179)*s23r+64.D0/3.D0*t394*t392*t203*t62+64.D0/3.D0*t436*t445*t1
     #30*s34r+64.D0/3.D0*(t461+t468+t483+t499)*l3*t184)*l1*r*t47+(64.D0/
     #3.D0*t26*t425*t130*s34r+64.D0/3.D0*t425*l3*t184)*s23r+64.D0/3.D0*t
     #436*(t92+t93+t117-t156+t155+t84+t119+t128-t157-t158)*t130*t62-64.D
     #0/3.D0*t53*(-t451-4*t452+2*t455+t454-t81-t437-t438-t348-t457-t360+
     #t527+t531+t534)*t130*s34r-64.D0/3.D0*(t546+t548)*t130)*s13r
      s5 = s6+s7
      s1 = s4+s5
      s2 = s1+64.D0/3.D0*(t556+t561)*l3*t184*s34r+((-64.D0/3.D0*t567*t6*
     #t355-64.D0/3.D0*t570*t1*t179)*s23r-64.D0/3.D0*t576*(t154+t152+t128
     #+t349+t350-4*t142+z)*t130*t62-64.D0/3.D0*(t591+t600)*t1*t180+64.D0
     #/3.D0*(t605+t606)*l3*t184)*l1*r*t47
      s3 = s2+64.D0/3.D0*t576*t174*t62*s34r
      s4 = s3+(-64.D0/3.D0*t191*t1*t392*t130*t23+(128.D0/3.D0*t624*t407+
     #64.D0/3.D0*t627*t5*t628*t130)*l1*r*t47-64.D0/3.D0*t624*t406*t62-64
     #.D0/3.D0*t627*t5*t628*t130*s34r-64.D0/3.D0*t645*(-t84+t156)*t130)*
     #t651
      s5 = s4
      s7 = (512.D0/3.D0*t3*t653*t23+((-512.D0/3.D0*t26*s2l1-512.D0/3.D0*
     #t230)*s23r-1024.D0/3.D0*t3*t653*s34r-512.D0/3.D0*t37*t95)*l1*r*t47
     #+(512.D0/3.D0*t53*s2l1+512.D0/3.D0*t291+512.D0/3.D0*t263)*s23r-512
     #.D0/3.D0*t58*s34r*s2l1+512.D0/3.D0*t3*t653*t62+512.D0/3.D0*t37*t95
     #*s34r+512.D0/3.D0*t263)*s1l2
      s8 = ((-512.D0/3.D0*t687*t689+512.D0/3.D0*t49*l3*s23r+t696*s2l2)*s
     #1l1+(-512.D0/3.D0*t58*t2*s2l1*t689+512.D0/3.D0*t53*s2l1*s23r+t696*
     #s2l1)*s1l2-512.D0/3.D0*t136*t710*t23+(1024.D0/3.D0*t136*t137*s23r+
     #(-1024.D0/3.D0*t687-t720-t723)*s2l1+t731*s2l2+1024.D0/3.D0*t136*t7
     #10*s34r-512.D0/3.D0*t84*t718*t736)*l1*r*t47+((512.D0/3.D0*t54+512.
     #D0/3.D0*t77)*s2l1+(512.D0/3.D0*t263+512.D0/3.D0*t294)*s2l2-1024.D0
     #/3.D0*t136*t137*s34r-512.D0/3.D0*t84*t4*t243)*s23r+((1024.D0/3.D0*
     #t316+1024.D0/3.D0*t645)*s2l2+(t720+t723)*s34r-512.D0/3.D0*t764*s3l
     #2+512.D0/3.D0*t768*s4l2)*s2l1+(-t731*s34r-512.D0/3.D0*t764*s3l1+51
     #2.D0/3.D0*t768*s4l1)*s2l2-512.D0/3.D0*t136*t710*t62+512.D0/3.D0*t8
     #4*t5*t2*t736*s34r+t310)*s13r
      s6 = s7+s8
      t791 = s5+s6
      taux3 = t276+t791
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t8 = r**2
      t9 = 1+t8
      t10 = r+t
      t14 = 1+r*t
      t16 = t9/t10/t14
      t17 = l1*r*t16
      t24 = l3*t1
      t28 = -2048.D0/3.D0*t24-2048.D0/3.D0*t3*t4*s34r
      t43 = t4*t2
      t44 = l3*s3l2
      t47 = t4*s4l2
      t50 = -4096.D0/3.D0*t43*t44-2048.D0/3.D0*t3*t47
      t52 = l3*s3l1
      t55 = t4*s4l1
      t58 = 4096.D0/3.D0*t43*t52+2048.D0/3.D0*t3*t55
      t65 = t24*s4l2
      t70 = t24*s4l1
      t78 = 2*l3
      t80 = t1*(t78-1)
      t95 = l3**2
      t96 = t95*t2
      t98 = t96*s4l2
      t102 = t4*t95
      t103 = t2*t1
      t104 = t103*s34r
      t112 = -t1
      t113 = t112*l3
      t121 = s34r**2
      t126 = 2048.D0/3.D0*t102*t103*t121+2048.D0/3.D0*t96*s34r
      t131 = t96*s4l1
      t153 = t95*t103
      t154 = t4*s3l2
      t157 = t4*s3l1
      t162 = l1**2
      t165 = t9**2
      t166 = t10**2
      t169 = t14**2
      t171 = t165/t166/t169
      t183 = s4l2*s3l1
      t184 = t96*t183
      t186 = s3l2*s4l1
      t187 = t96*t186
      t196 = 2*z
      t197 = z*l3
      t200 = (1-t196-t78+2*t197)*t1
      t221 = 2*l3*t1
      t245 = t196-t197-1+l3
      t251 = -l3+t197-z
      t269 = s3l1*s2l2
      t270 = s4l2*s2l1
      t271 = s3l2*s2l1
      t272 = s4l1*s2l2
      t273 = -t269+t270+t271-t272
      t275 = 1/s134r
      t280 = t4*t275
      t287 = s13r**2
      t289 = s3l2*s1l1
      t290 = s3l1*s1l2
      t291 = s1l2*s2l1
      t292 = s1l2*s4l1
      t293 = s1l1*s2l2
      t294 = s1l1*s4l2
      t303 = -t290-t292+t289+t294
      t304 = t303*t95
      t305 = t2*t275
      t308 = t153*t4
      t310 = 2*t293
      t314 = 2*t291
      t319 = t293*z
      t320 = s1l1*l3
      t321 = t320*s2l2
      t322 = 2*t321
      t323 = t271*z
      t324 = t269*z
      t325 = t320*s3l2
      t326 = 2*t325
      t327 = t291*l3
      t328 = 2*t327
      t329 = s1l2*z
      t330 = t329*s2l1
      t331 = s1l2*l3
      t332 = t331*s3l1
      t333 = 2*t332
      t334 = t290*z
      t335 = t289*z
      t336 = t292*z
      t337 = t292*l3
      t338 = 2*t337
      t339 = z*s3l2
      t340 = t320*t339
      t341 = 2*t340
      t342 = z*s3l1
      t343 = t331*t342
      t344 = 2*t343
      t345 = t319+t322+t323-t324+t326-t328-t330-t333-t334+t335-t336-t338
     #-t341+t344
      t346 = t294*l3
      t347 = 2*t346
      t348 = t292*t197
      t349 = 2*t348
      t350 = t294*t197
      t351 = 2*t350
      t352 = t294*z
      t353 = t272*z
      t354 = t270*z
      t355 = t293*t197
      t356 = 2*t355
      t357 = l3*s2l1
      t358 = t329*t357
      t359 = 2*t358
      t360 = t347+t349-t351+t352-t353+t354-t356+t359-t289+t290+t291+t292
     #-t293-t294
      t363 = t1*t275
      t370 = t305*s34r
      t384 = l3*s2l2
      t390 = t44*s2l1-t52*s2l2-t319-t322+t384*t342-t326+t328+t330+t333+t
     #334-t335+t336+t338-l3*s4l1*s2l2+t341-t344-t347-t357*s4l2*z
      t396 = t384*s4l1*z-t270+t269+l3*s4l2*s2l1-t271+t272-t349+t351-t352
     #+t356-t359+t289-t290-t291-t292+t293+t294-t357*t339
      t402 = t321-t355+t358+t323-t324+t325-t340-t327-t332+t343+t346-t337
     #-t353+t354+t348-t350
      t407 = -t291+t293
      t408 = t407*t4
      t409 = t95*l3
      t410 = t2**2
      t429 = t103*t275
      t434 = 3*t407*t4*t409
      t435 = t410*t275
      t438 = 2*t319
      t439 = 3*t321
      t440 = 3*t355
      t441 = 3*t358
      t442 = 3*t327
      t443 = 2*t330
      t444 = -t272+t270-t271+t438+t439-t440+t441-t442-t443+t269-t186+t18
     #3+t314-t310
      t452 = t272-t270-t271+t269+t291-t293
      t465 = 2*t271
      t472 = 2*t269
      t475 = -t272+t270-t465+4*t319+6*t321-6*t355+6*t358-6*t327-4*t330+t
     #472-t186+t183+4*t291+t292-4*t293-t294
      t483 = -2*t272+2*t270-t465+t319+t439-t440+t441-t442-t330+t472-2*t1
     #86+2*t183+t291-t293
      t500 = l3*t275
      t510 = t443-t269+t310+t271-t314-t292-t439+t294-t441+t442+t440-t438
      t515 = -t272+t270-t465+t319+t439-t440+t441-t442-t330+t472-t186+t18
     #3+t291+t292-t293-t294
      s3 = -((2048.D0/3.D0*t3*t4*s2l2*t17-2048.D0/3.D0*t1*s2l2*l3*s23r+t
     #28*s2l2)*s1l1+(-2048.D0/3.D0*t3*t4*s2l1*t17+2048.D0/3.D0*t24*s2l1*
     #s23r-t28*s2l1)*s1l2+(t50*s2l1+t58*s2l2)*l1*r*t16+((2048.D0/3.D0*t2
     #4*s3l2+2048.D0/3.D0*t65)*s2l1+(-2048.D0/3.D0*t24*s3l1-2048.D0/3.D0
     #*t70)*s2l2)*s23r+(-t50*s34r+2048.D0/3.D0*t80*s3l2+2048.D0/3.D0*t65
     #)*s2l1+(-t58*s34r-2048.D0/3.D0*t80*s3l1-2048.D0/3.D0*t70)*s2l2)*s1
     #3r-(((2048.D0/3.D0*t96*s2l2+2048.D0/3.D0*t98)*s23r-2048.D0/3.D0*t1
     #02*t104*s2l2)*l1*r*t16+(2048.D0/3.D0*t113*s2l2-2048.D0/3.D0*t96*s4
     #l2*s34r+2048.D0/3.D0*t113*s4l2)*s23r+t126*s2l2)*s1l1
      s2 = s3-(((-2048.D0/3.D0*t96*s2l1-2048.D0/3.D0*t131)*s23r+2048.D0/
     #3.D0*t102*t104*s2l1)*l1*r*t16+(-2048.D0/3.D0*t113*s2l1+2048.D0/3.D
     #0*t96*s4l1*s34r-2048.D0/3.D0*t113*s4l1)*s23r-t126*s2l1)*s1l2-(2048
     #.D0/3.D0*t153*t154*s2l1-2048.D0/3.D0*t153*t157*s2l2)*t162*t8*t171
      s3 = s2-(((-4096.D0/3.D0*t96*s3l2-2048.D0/3.D0*t98)*s2l1+(4096.D0/
     #3.D0*t96*s3l1+2048.D0/3.D0*t131)*s2l2+2048.D0/3.D0*t184-2048.D0/3.
     #D0*t187)*s23r+((-2048.D0/3.D0*t153*t154+2048.D0/3.D0*t153*t47)*s34
     #r-2048.D0/3.D0*t200*t44)*s2l1+((2048.D0/3.D0*t153*t157-2048.D0/3.D
     #0*t153*t55)*s34r+2048.D0/3.D0*t200*t52)*s2l2)*l1*r*t16-((2048.D0/3
     #.D0*t96*s3l2*s34r+2048.D0/3.D0*t221*s3l2+2048.D0/3.D0*t65)*s2l1+(-
     #2048.D0/3.D0*t96*s3l1*s34r-2048.D0/3.D0*t221*s3l1-2048.D0/3.D0*t70
     #)*s2l2+(-2048.D0/3.D0*t184+2048.D0/3.D0*t187)*s34r+2048.D0/3.D0*t1
     #13*t183-2048.D0/3.D0*t113*t186)*s23r
      s1 = s3-(-2048.D0/3.D0*t102*t103*s4l2*t121+(-2048.D0/3.D0*t24*t245
     #*s3l2-2048.D0/3.D0*t98)*s34r+2048.D0/3.D0*t251*s3l2)*s2l1-(2048.D0
     #/3.D0*t102*t103*s4l1*t121+(2048.D0/3.D0*t24*t245*s3l1+2048.D0/3.D0
     #*t131)*s34r-2048.D0/3.D0*t251*s3l1)*s2l2-(256.D0/3.D0*t3*t273*t4*t
     #275*t17-256.D0/3.D0*t3*t273*t280*s34r+256.D0/3.D0*t112*t273*t280)*
     #t287
      s2 = s1-(-256.D0/3.D0*(t289-t290-t291-t292+t293+t294)*t4*t153*t275
     #*t162*t8*t171+(256.D0/3.D0*t304*t305*s23r+256.D0/3.D0*t308*(-t271+
     #t269+t272-t270-2*t290+t310+2*t294+2*t289-2*t292-t314)*t275*s34r-25
     #6.D0/3.D0*(t345+t360)*l3*t363)*l1*r*t16+(-256.D0/3.D0*t304*t370-25
     #6.D0/3.D0*t303*l3*t363)*s23r-256.D0/3.D0*t308*(-t271+t269+t272-t27
     #0-t290-t292+t294+t289-t291+t293)*t275*t121-256.D0/3.D0*t24*(t390+t
     #396)*t275*s34r+256.D0/3.D0*t402*t275)*s13r-256.D0/3.D0*t408*t409*t
     #410*t275*t162*l1*t8*r*t165*t9/t166/t10/t169/t14-(-256.D0/3.D0*t273
     #*t409*t429*s23r-256.D0/3.D0*t434*t435*s34r+256.D0/3.D0*t444*t95*t3
     #05)*t162*t8*t171-((-256.D0/3.D0*t452*t409*t429*s34r+512.D0/3.D0*t2
     #73*t95*t305)*s23r+256.D0/3.D0*t434*t435*t121-256.D0/3.D0*t475*t95*
     #t370-256.D0/3.D0*t483*l3*t363)*l1*r*t16
      taux4 = s2-(-256.D0/3.D0*t409*t103*t407*t275*t121+256.D0/3.D0*t452
     #*t95*t370-256.D0/3.D0*t273*t1*t500)*s23r+256.D0/3.D0*t408*t409*t43
     #5*t121*s34r+256.D0/3.D0*t96*t510*t275*t121-256.D0/3.D0*t515*t1*t50
     #0*s34r-256.D0/3.D0*(-t272+t270-t271+t321-t355+t358-t327+t269-t186+
     #t183)*t275

CCC Coefficients

      tcoff9 = taux3
      tcoff10 = taux3
      tcoff11 = taux3
      tcoff12 = taux3
      tcoff13 = taux4
      tcoff14 = taux4
      tcoff15 = -taux4
      tcoff16 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet10 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet11 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet12 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet13 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet14 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet15 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet16 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff16)

CCC Sum

      sq3num2e0 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*
     #tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tc
     #off16
      RETURN
      END

      double precision function sq3num2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = l3**2
      t3 = t2*l3
      t4 = t1*t3
      t5 = -1+z
      t6 = t5**2
      t7 = t6**2
      t8 = t7*z
      t12 = t7*t5
      t13 = t3*t12
      t14 = t1**2
      t15 = z*t14
      t19 = t6*t5
      t20 = t2*t19
      t22 = t20*t1*s3l2
      t24 = t1*s4l2
      t25 = t20*t24
      t28 = (-1024.D0/3.D0*t22+512.D0/3.D0*t25)*s3l1
      t29 = t1*s4l1
      t31 = t20*t29*s3l2
      t32 = 512.D0/3.D0*t31
      t34 = l1**2
      t36 = r**2
      t38 = 1+t36
      t39 = t38**2
      t40 = r+t
      t41 = t40**2
      t45 = 1+r*t
      t46 = t45**2
      t48 = t39/t41/t46
      t50 = t1*z
      t52 = t2**2
      t53 = 1/s134r
      t58 = 3*z
      t59 = s4l1*s4l2
      t60 = 2*t59
      t61 = s3l1*s3l2
      t62 = 2*t61
      t63 = s1l1*s2l2
      t64 = s1l2*s2l1
      t65 = s3l1*s2l2
      t66 = s3l2*s2l1
      t67 = z**2
      t68 = l3*t67
      t69 = 6*t68
      t70 = z*l3
      t71 = 6*t70
      t72 = 5*t67
      t80 = s3l2*s4l1
      t81 = 2*t80
      t82 = s4l2*s3l1
      t83 = 2*t82
      t84 = s3l1*s1l2
      t85 = s1l2*s4l1
      t86 = 2*t85
      t87 = s3l2*s1l1
      t88 = s1l1*s4l2
      t89 = 2*t88
      t90 = 4*z
      t91 = 4*t61
      t92 = s4l1*s2l2
      t93 = 2*t92
      t94 = s4l2*s2l1
      t95 = 2*t94
      t96 = -t81-t83-t84-t86-t87-t89-t90+t91+t65+t93+t66+t95
      t99 = t3*t53
      t100 = t99*s34r
      t103 = l3*s2l1
      t104 = s4l2*z
      t105 = t103*t104
      t106 = 3*t105
      t107 = l3*s2l2
      t108 = s4l1*z
      t109 = t107*t108
      t110 = 3*t109
      t111 = t85*t70
      t112 = 3*t111
      t113 = t88*t70
      t114 = 3*t113
      t115 = l3*s3l2
      t116 = t115*t108
      t117 = 3*t116
      t118 = l3*s3l1
      t119 = t118*t104
      t120 = 3*t119
      t121 = t92*z
      t122 = 3*t121
      t123 = t94*z
      t124 = 3*t123
      t125 = t82*z
      t127 = l3*s4l1
      t128 = t127*s3l2
      t129 = 3*t128
      t130 = t118*s4l2
      t131 = 3*t130
      t132 = s3l2*z
      t133 = t118*t132
      t134 = 6*t133
      t135 = 4*t67
      t136 = 3*t94
      t137 = t80*z
      t139 = t61*z
      t141 = -t106-t110+t112+t114+t117+t120+t122+t124-2*t125-t129-t131+t
     #90-t134-t135-t136-2*t137+4*t139
      t142 = t118*s3l2
      t143 = 6*t142
      t144 = t85*z
      t146 = t88*z
      t148 = t88*l3
      t149 = 3*t148
      t150 = t127*s2l2
      t151 = 3*t150
      t153 = l3*s4l2*s2l1
      t154 = 3*t153
      t155 = t85*l3
      t156 = 3*t155
      t157 = 3*t92
      t158 = t143-2*t144-2*t146-t149+t151+t154-t156-t60-t91-t65-t66+t69-
     #t71+t80+t82+t86+t89-t157
      t161 = t6*t53
      t168 = t3*t19
      t169 = s23r**2
      t173 = t2*t6
      t174 = t173*s3l2
      t176 = t173*s4l2
      t180 = t173*s3l1
      t182 = t173*s4l1
      t190 = 1024.D0/3.D0*t173*t61
      t192 = 1024.D0/3.D0*t173*t59
      t193 = 1-l3+t70-z
      t194 = 2*t193
      t196 = t6*z
      t199 = 2*l3
      t200 = 2*t70
      t201 = t58-1+t199-t200
      t207 = t1*t2
      t212 = t19*s4l2
      t216 = l3*t5
      t217 = l3*t6
      t218 = t216-t217
      t223 = t19*s4l1
      t230 = s34r**2
      t239 = 1+t200-t199
      t242 = 512.D0/3.D0*t20*t50*t239
      t244 = 512.D0/3.D0*t20*t50
      t247 = t199-1
      t248 = t247*t6
      t249 = t248*t115
      t252 = (-3+z+t199-t200)*l3
      t256 = 1024.D0/3.D0*t249+512.D0/3.D0*t252*t5*s4l2
      t261 = 512.D0/3.D0*t252*t5*s4l1*s3l2
      t268 = t38/t40/t45
      t273 = z*s3l1
      t275 = -t217-t173
      t280 = 2-z
      t287 = (-512.D0/3.D0*t174+512.D0/3.D0*t176)*s34r
      t288 = t216*s3l2
      t290 = t216*s4l2
      t296 = (-512.D0/3.D0*t180+512.D0/3.D0*t182)*s34r
      t297 = t216*s3l1
      t299 = t216*s4l1
      t320 = -l3*t5
      t331 = t34*l1
      t332 = t36*r
      t340 = t39*t38/t41/t40/t46/t45
      t349 = t14*t53
      t350 = t349*s34r
      t354 = s2l2*s2l1
      t355 = 2*t354
      t357 = 12*t70
      t358 = 7*t67
      t359 = -t84+t85-t87+t88+9*z-t62-t157-t355-t136+12*t68-t357-t358
      t369 = 4*t68
      t370 = 4*t70
      t371 = t80+t82+t84+t87-t90+t62+t63+t64-t369+t370+t135
      t377 = z*t53
      t378 = t377*t230
      t381 = t207*t19
      t382 = 2*t84
      t383 = 2*t87
      t388 = -t382+t86-t383+t89+12*z-t91-t157-t355-t136+16*t68-16*t70-10
     #*t67
      t395 = 2*t111
      t396 = 2*t113
      t399 = 5*z
      t400 = 4*t133
      t402 = t107*s2l1
      t405 = t103*s2l2*z
      t407 = -6*t105-6*t109+t395+t396+6*t121+6*t123+t399-t400-8*t67+4*t4
     #02-4*t405
      t408 = t67*z
      t409 = 3*t408
      t411 = 2*t139
      t412 = 4*t142
      t413 = 2*t148
      t416 = 2*t155
      t418 = t409-6*t94+t411+t412-t144-t146-t413+6*t150+6*t153-t416-4*t3
     #54
      t420 = t354*z
      t422 = t2*t408
      t424 = t2*t67
      t426 = z*t2
      t428 = l3*t408
      t432 = 4*t420+12*t422-24*t424+12*t426-14*t428-t60-t62+32*t68-18*t7
     #0-t84+t85
      t433 = t66*z
      t434 = t65*z
      t435 = s1l1*l3
      t436 = t435*s3l2
      t437 = 2*t436
      t438 = s1l2*l3
      t439 = t438*s3l1
      t440 = 2*t439
      t441 = t84*z
      t442 = t87*z
      t444 = t435*t132
      t445 = 2*t444
      t446 = t438*t273
      t447 = 2*t446
      t448 = -t87+t88+t433+t434+t437+t440+t441+t442-6*t92-t445-t447
      t452 = t5*t53
      t459 = t4*t7
      t467 = 2*z
      t468 = 2*t68
      t469 = 2*t67
      t470 = t80+t82+t84+t87-t467+t62+t63+t64-t468+t200+t469
      t476 = 3*t67
      t483 = 2*t402
      t484 = 2*t405
      t485 = -t106-t110+t395+t396+t122+t124+t58-t400-6*t67+t483-t484
      t486 = t409-t136+t411+t412-t144-t146-t413+t151+t154-t416-t355
      t488 = 2*t420
      t494 = t488+8*t422-16*t424+8*t426-10*t428-t60-t62+22*t68-t357-t84+
     #t85
      t495 = -t87+t88+t433+t434+t437+t440+t441+t442-t157-t445-t447
      t502 = 2*t133
      t503 = 2*t142
      t504 = -t106-t110+t111+t113+t122+t124-t502+t483-t484-t136+t503-t14
     #8+t151+t154-t155
      t509 = 8*t68
      t511 = -t355+t488+4*t422-8*t424+4*t426-3*t428-t60+t509-5*t70+t433+
     #t434+t436+t439-t157-t444-t446
      t520 = t217*t1*s34r
      t536 = t230*s34r
      t541 = t217*t1*s2l2
      t543 = l1*r*t268
      t551 = -1024.D0/3.D0*t216-1024.D0/3.D0*t520
      t569 = t1*t6
      t571 = 512.D0/3.D0*t569*t115
      t573 = 512*t217*t24
      t580 = 512.D0/3.D0*t569*t118+512*t217*t29
      t584 = -t217*t50*t239
      t588 = t217*t50*(z+t199-t200-1)
      t603 = -l3+t70+1
      t612 = t5*t1
      t618 = 1-z-t216
      t623 = t320-t5*t247
      t650 = t34*t36*t48
      t657 = t19*s3l2
      t668 = 512.D0/3.D0*t173*s34r-512.D0/3.D0*t216
      t670 = 512.D0/3.D0*t288
      t687 = -2+l3-t70
      s2 = (1024.D0/3.D0*t4*t8*s23r+512.D0/3.D0*t13*t15*s34r+t28+t32)*t3
     #4*t36*t48+((256.D0/3.D0*t50*t12*t52*t53*s34r+64.D0/3.D0*(t58-t60+t
     #62+t63+t64+t65+t66+t69-t71-t72)*t3*t19*t53)*s23r+64.D0/3.D0*t96*t1
     #*t7*t100-64.D0/3.D0*(t141+t158)*t2*t161)*t34*t36*t48+(-512.D0/3.D0
     #*t168*z*t169+((512*t174+512.D0/3.D0*t176)*s2l1+(512*t180+512.D0/3.
     #D0*t182)*s2l2-2048.D0/3.D0*t4*t8*s34r+t190-t192-512.D0/3.D0*t194*t
     #2*t196+512.D0/3.D0*t201*t2*t196)*s23r+(-1024.D0/3.D0*t207*t19*s34r
     #*s2l2-512.D0/3.D0*t207*t212*s34r-512.D0/3.D0*t218*s3l2)*s2l1+(-512
     #.D0/3.D0*t207*t223*s34r-512.D0/3.D0*t218*s3l1)*s2l2-1024.D0/3.D0*t
     #13*t15*t230+((2048.D0/3.D0*t22-1024.D0/3.D0*t25)*s3l1-1024.D0/3.D0
     #*t31-t242-t244)*s34r+t256*s3l1+t261)*l1*r*t268
      s3 = s2+(512.D0/3.D0*t207*t223*t230+(512.D0/3.D0*t216*t273-512.D0/
     #3.D0*t275*s4l1)*s34r+512.D0/3.D0*t280*s3l1)*s2l2
      s4 = s3
      s6 = ((t287-512*t288-512.D0/3.D0*t290)*s2l1+(t296-512*t297-512.D0/
     #3.D0*t299)*s2l2+1024.D0/3.D0*t4*t8*t230+(-t190+t192-512.D0/3.D0*t1
     #73*z*(-2+3*l3-3*t70+t58)+512.D0/3.D0*t173*z*t193)*s34r-1024.D0/3.D
     #0*t216*t61-1024.D0/3.D0*t320*t59-512.D0/3.D0*t201*z*t216)*s23r
      s7 = (256.D0/3.D0*t15*t13*t53*t331*t332*t340+(-128.D0/3.D0*t50*t7*
     #t99*s23r-512.D0/3.D0*t13*z*t350-64.D0/3.D0*t359*t2*t1*t19*t53)*t34
     #*t36*t48+(-64.D0/3.D0*t371*t2*t161*s23r+256.D0/3.D0*t13*t14*t378+6
     #4.D0/3.D0*t381*t388*t53*s34r+64.D0/3.D0*(t407+t418+t432+t448)*l3*t
     #452)*l1*r*t268+(128.D0/3.D0*t459*t378+64.D0/3.D0*t173*(t62+t82+t80
     #+t64+t84+t63+t87)*t53*s34r+64.D0/3.D0*t470*l3*t452)*s23r+64.D0/3.D
     #0*t381*(t84-t85+t87-t88-t58+t62-t369+t370+t476)*t53*t230-64.D0/3.D
     #0*t216*(t485+t486+t494+t495)*t53*s34r-64.D0/3.D0*(t504+t511)*t53)*
     #s13r
      s5 = s6+s7
      s1 = s4+s5
      s2 = s1+((1024.D0/3.D0*t207*t19*t230+1024.D0/3.D0*t520)*s2l2+512.D
     #0/3.D0*t207*t212*t230+(512.D0/3.D0*t216*t132-512.D0/3.D0*t275*s4l2
     #)*s34r+512.D0/3.D0*t280*s3l2)*s2l1+512.D0/3.D0*t13*t15*t536
      s3 = s2+((1024.D0/3.D0*t541*t543-1024.D0/3.D0*t5*s2l2*l3*s23r+t551
     #*s2l2)*s1l1+(1024.D0/3.D0*t217*t1*s2l1*t543-1024.D0/3.D0*t216*s2l1
     #*s23r+t551*s2l1)*s1l2+(-512*t20*t50*s23r+(1024.D0/3.D0*t541+t571+t
     #573)*s2l1+t580*s2l2-512.D0/3.D0*t584+512.D0/3.D0*t588)*l1*r*t268+(
     #(-1024.D0/3.D0*t288-1024.D0/3.D0*t290)*s2l1+(-1024.D0/3.D0*t297-10
     #24.D0/3.D0*t299)*s2l2+512*t20*t50*s34r-512.D0/3.D0*t70*t5*t603+102
     #4.D0/3.D0*t70*t5*t194)*s23r+((-1024.D0/3.D0*t520-1024.D0/3.D0*t612
     #)*s2l2+(-t571-t573)*s34r+512.D0/3.D0*t618*s3l2+512.D0/3.D0*t623*s4
     #l2)*s2l1+(-t580*s34r+512.D0/3.D0*t618*s3l1+512.D0/3.D0*t623*s4l1)*
     #s2l2-512.D0/3.D0*(-t584+t588)*s34r-512.D0/3.D0*t70*t5*t603+512.D0/
     #3.D0*z*t193*t216)*s13r
      t704 = s3+(-512.D0/3.D0*t22*t650+((512.D0/3.D0*t173*s2l2+512.D0/3.
     #D0*t174-512.D0/3.D0*t176)*s23r+1024.D0/3.D0*t207*t657*s34r+512.D0/
     #3.D0*t249)*l1*r*t268+(t668*s2l2+t287-t670-512.D0/3.D0*t320*s4l2)*s
     #23r+1024.D0/3.D0*t217*s34r*s2l2-512.D0/3.D0*t207*t657*t230-512.D0/
     #3.D0*t248*t115*s34r-t670)*s1l1+(-1024.D0/3.D0*t288-512.D0/3.D0*t68
     #7*s4l2)*s3l1+(-t256*s3l1-t261+512.D0/3.D0*t196*l3*(t426-t2+l3)+512
     #.D0/3.D0*t196*l3*t1)*s34r
      t711 = -t105-t109+t111+t113+t116+t119+t121+t123-t128-t130-t502-t94
      t712 = t503-t148+t150+t153-t155-t60-t65-t66+t468-t200-t80-t82-t92
      t727 = 2*t105+2*t109-t395-t396-2*t116-2*t119-2*t121-2*t123+t125+2*
     #t128+2*t130
      t728 = t103*t132
      t729 = t107*t273
      t730 = -z+t400+t67+t728+t729+t95+t137-t411-t412+t144+t146
      t736 = t413-2*t150-2*t153+t416+t60+t62+3*t65+3*t66-t369+t370+t81
      t737 = t115*s2l1
      t738 = t118*s2l2
      t739 = t83-t85-t88-t737-t738-t433-t434+t436+t439+t93-t444-t446
      t767 = -t84+t85-t87+t88-z-t60+t62+t63+t64+t65+t66+t369-t370-t476
      t769 = t161*s34r
      t777 = s1l2*z
      t779 = t105+t109-t111-t113-t116-t119-t121-t123+t125+t128+t130-z+t7
     #77*s2l1+t502
      t783 = t67+t728+t729+t63*z+t435*s2l2+t64*l3+t94+t137-t411-t503+t14
     #4+t146+t148-t150
      t787 = 2*t65
      t788 = 2*t66
      t789 = -t153+t155-t63*t70-t777*t103+t62-t63-t64+t787+t788-t468+t20
     #0+t80+t82-t84
      t790 = -t85-t87-t88-t737-t738-t433-t434+t437+t440+t441+t442+t92-t4
     #45-t447
      t799 = t1*t53
      t803 = 8*t70
      t804 = -t84+t85-t87+t88+t58-t60+t62+t63+t64+t65+t66+t509-t803-t358
      t808 = 4*t59
      t817 = -t80-t82-t382-t85-t383-t88-t467+t62-t63-t64+t65+t92+t66+t94
      t833 = 4*t105+4*t109-4*t111-4*t113-4*t116-4*t119-4*t121-4*t123+3*t
     #125+4*t128+4*t130
      t844 = -t399+8*t133+t72+2*t728+2*t729+4*t94+3*t137-6*t139-8*t142+3
     #*t144+3*t146+4*t148
      t853 = -4*t150-4*t153+4*t155+t60+6*t61+4*t65+4*t66-t509+t803-t84-3
     #*t85-t87
      t860 = -3*t88-2*t737-2*t738-2*t433-2*t434+t437+t440+t441+t442+4*t9
     #2-t445-t447
      t866 = -t106-t110+t112+t114+t117+t120+t122+t124-t125-t129-t131+t46
     #7-t134-t469-t136-t137+t411
      t867 = t143-t144-t146-t149+t151+t154-t156-t808-t62-t787-t788+t69-t
     #71-t80-t82+t85+t88-t157
      t886 = t20*z
      t890 = t8*t2
      t893 = l3*t1
      t894 = t467-t469-t355-t66-t95-t65-t93-t370+t369
      t923 = s13r**2
      t925 = t19*s3l1
      t943 = 512.D0/3.D0*t297
      s1 = (512.D0/3.D0*t168*z*s34r+512.D0/3.D0*t173*z)*t169-64.D0/3.D0*
     #(t711+t712)*t53-512.D0/3.D0*t687*s4l1*s3l2+64.D0/3.D0*(t727+t730+t
     #736+t739)*l3*t452*s34r+(t28+t32+t242+t244)*t230+(-128.D0/3.D0*t1*t
     #52*t12*t377*s23r+64.D0/3.D0*(t88+t85-t62+t82-t94-t92+t80+t467)*t1*
     #t7*t3*t53)*t331*t332*t340
      s2 = s1+(64.D0/3.D0*t168*(-t87+t88-t84+t85-t200+t468-t469)*t53*t23
     #0+64.D0/3.D0*t767*t2*t769+64.D0/3.D0*(-z-t60+t62+t63+t64+t65+t66+t
     #468-t200-t67)*l3*t452)*s23r+64.D0/3.D0*t173*(t779+t783+t789+t790)*
     #t53*t230
      s3 = s2+((-128.D0/3.D0*t52*t12*z*t799*t230-64.D0/3.D0*t804*t19*t10
     #0-64.D0/3.D0*(-t808+t91+2*t63+2*t64+t787+t788+t69-t71-t135)*t2*t16
     #1)*s23r-64.D0/3.D0*t459*t817*t53*t230-64.D0/3.D0*(t833+t844+t853+t
     #860)*t2*t769+64.D0/3.D0*(t866+t867)*l3*t452)*l1*r*t268
      t959 = s3+64.D0/3.D0*t459*(-t64-t84-t63-t87)*t53*t536+(-128.D0/3.D
     #0*t8*t2*t14*t53*t650+(128.D0/3.D0*t886*t799*s23r+256.D0/3.D0*t890*
     #t350+64.D0/3.D0*t893*t6*t894*t53)*l1*r*t268+(-128.D0/3.D0*t886*t79
     #9*s34r-128.D0/3.D0*t217*t50*t53)*s23r-128.D0/3.D0*t890*t349*t230-6
     #4.D0/3.D0*t893*t6*t894*t53*s34r-64.D0/3.D0*t612*(-t200+t468-t66-t3
     #55-t95-t93-t65)*t53)*t923+(-512.D0/3.D0*t207*t925*t650+((512.D0/3.
     #D0*t173*s2l1+512.D0/3.D0*t180-512.D0/3.D0*t182)*s23r+1024.D0/3.D0*
     #t207*t925*s34r+512.D0/3.D0*t248*t118)*l1*r*t268+(t668*s2l1+t296-t9
     #43-512.D0/3.D0*t320*s4l1)*s23r+1024.D0/3.D0*t217*s34r*s2l1-512.D0/
     #3.D0*t207*t925*t230-512.D0/3.D0*t248*t118*s34r-t943)*s1l2
      taux3 = t704+t959
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t8 = r**2
      t9 = 1+t8
      t10 = r+t
      t11 = 1/t10
      t14 = 1+r*t
      t15 = 1/t14
      t16 = t9*t11*t15
      t17 = l1*r*t16
      t24 = l3*t1
      t28 = 10240.D0/3.D0*t24+10240.D0/3.D0*t3*t4*s34r
      t43 = t4*t2
      t44 = l3*s3l2
      t47 = t4*s4l2
      t50 = 6144*t43*t44+10240.D0/3.D0*t3*t47
      t52 = l3*s3l1
      t55 = t4*s4l1
      t58 = -6144*t43*t52-10240.D0/3.D0*t3*t55
      t65 = t24*s4l2
      t70 = t24*s4l1
      t81 = -t1
      t83 = 2*l3
      t86 = t1*(1-3*l3)+3*t81*(t83-1)
      t101 = l3**2
      t102 = t101*t2
      t104 = t102*s4l2
      t108 = t4*t101
      t109 = t2*t1
      t110 = t109*s34r
      t118 = t102*s34r
      t120 = 4*l3*t1
      t122 = -2048.D0/3.D0*t118+2048.D0/3.D0*t120
      t131 = s34r**2
      t135 = -10240.D0/3.D0*t108*t109*t131-10240.D0/3.D0*t118
      t140 = t102*s4l1
      t164 = t101*t109
      t165 = t4*s3l2
      t168 = t4*s3l1
      t173 = l1**2
      t176 = t9**2
      t177 = t10**2
      t178 = 1/t177
      t180 = t14**2
      t181 = 1/t180
      t182 = t176*t178*t181
      t184 = t102*s3l2
      t189 = t102*s3l1
      t194 = s4l2*s3l1
      t195 = t102*t194
      t197 = s3l2*s4l1
      t198 = t102*t197
      t207 = 2*z
      t208 = z*l3
      t209 = 2*t208
      t214 = 3*z
      t220 = -3*(1-t207-t83+t209)*t1*l3+(-2+t214+4*l3-4*t208)*t1*l3
      t243 = 9*t81*l3
      t276 = t24*(-t214+t209+2-t83)-3*t24*(t207-t208-1+l3)
      t285 = 5*l3-5*t208+4*z
      t303 = s4l1*s2l2
      t304 = s4l2*s2l1
      t305 = -t303+t304
      t307 = s3l1*s2l2
      t308 = s3l2*s2l1
      t309 = -t307+t304+t308-t303
      t313 = (-t3*t305+2*t3*t309)*t4
      t314 = 1/s134r
      t331 = s13r**2
      t333 = s3l2*s1l1
      t334 = s3l1*s1l2
      t335 = s1l1*s2l2
      t336 = s1l1*s4l2
      t337 = s1l2*s2l1
      t338 = s1l2*s4l1
      t343 = 4*t307
      t344 = 2*t303
      t345 = 4*t308
      t346 = 2*t304
      t347 = t197-t194-t334-t338+t333+t336-t335+t337-t343-t344+t345+t346
      t358 = -t333+t334-t337+t194-t197+t335
      t359 = t358*t101
      t360 = t2*t314
      t363 = 2*t336
      t364 = 2*t333
      t365 = 2*t194
      t366 = 2*t334
      t367 = 2*t338
      t368 = 7*t307
      t369 = 7*t308
      t370 = 2*t197
      t374 = 2*t335
      t375 = 2*t337
      t383 = t197*z
      t384 = t335*z
      t385 = 3*t384
      t386 = l3*s2l1
      t387 = z*s3l2
      t388 = t386*t387
      t390 = l3*s2l2
      t391 = z*s3l1
      t392 = t390*t391
      t394 = t44*s2l1
      t396 = t52*s2l2
      t398 = s1l1*l3
      t399 = t398*s2l2
      t401 = t308*z
      t403 = t307*z
      t405 = t398*s3l2
      t406 = 2*t405
      t407 = t383-t385+t333-t334-8*t388+8*t392+8*t394-8*t396-6*t399+4*t4
     #01-4*t403-t406
      t408 = t337*l3
      t410 = t335*t208
      t412 = s1l2*z
      t413 = t412*t386
      t415 = t398*t387
      t416 = 2*t415
      t417 = 3*t335
      t418 = 3*t337
      t419 = 3*t197
      t420 = 3*t194
      t421 = t412*s2l1
      t422 = 3*t421
      t423 = s1l2*l3
      t424 = t423*s3l1
      t425 = 2*t424
      t426 = 6*t408+6*t410-6*t413+t336+t416+t417-t418-t338-t419+t420+t42
     #2+t425
      t428 = t334*z
      t429 = t333*z
      t430 = t338*z
      t431 = t336*z
      t432 = t336*l3
      t433 = 2*t432
      t434 = t423*t391
      t435 = 2*t434
      t436 = z*s4l2
      t437 = t386*t436
      t439 = z*s4l1
      t440 = t390*t439
      t442 = l3*s4l1
      t443 = t442*s2l2
      t446 = l3*s4l2*s2l1
      t448 = t338*l3
      t449 = 2*t448
      t450 = t303*z
      t451 = t428-t429+t430-t431-t433-t435-4*t437+4*t440-4*t443+4*t446+t
     #449-t450
      t452 = t304*z
      t453 = t338*t208
      t454 = 2*t453
      t455 = t336*t208
      t456 = 2*t455
      t457 = t44*t439
      t458 = 2*t457
      t459 = t52*t436
      t460 = 2*t459
      t461 = t194*z
      t462 = t442*s3l2
      t463 = 2*t462
      t464 = t52*s4l2
      t465 = 2*t464
      t470 = t452-t454+t456-t458+t460-t461+t463-t465+4*t303-4*t304+8*t30
     #7-8*t308
      t474 = t1*t314
      t481 = t360*s34r
      t488 = 3*t307
      t489 = 3*t308
      t507 = -t383-t384+t333-t334+7*t388-7*t392-7*t394+7*t396-5*t401+5*t
     #403-t406
      t508 = t336+t416+t335-t337-t338+t419-t420+t421+t425+t428-t429
      t510 = 2*t437
      t511 = 2*t440
      t512 = 2*t443
      t513 = 2*t446
      t514 = t430-t431-t433-t435+t510-t511+t512-t513+t449+t450-t452
      t515 = -t454+t456+t458-t460+t461-t463+t465-t344+t346-t368+t369
      t523 = -t384+t333-t334-t388+t392+t394-t396-2*t399-t406+2*t408+2*t4
     #10-2*t413+t336+t416+t335-t337-t338+t421
      t524 = t425+t428-t429+t430-t431-t433-t435-t437+t440-t443+t446+t449
     #-t454+t456+t303-t304+t307-t308
      t536 = 3*t399
      t537 = 3*t408
      t538 = 3*t410
      t539 = 3*t413
      t540 = -4*t388+4*t392+4*t394-4*t396-t536-t405+t537+t538-t539+t415-
     #t370+t365+t424-t432-t434-t510
      t541 = t511-t512+t513+t448+t450-t452-t453+t455-t457+t459+t462-t464
     #+t344-t346+t343-t345
      t547 = -t337+t335
      t549 = t101*l3
      t550 = t2**2
      t551 = t549*t550
      t553 = 2*t547*t4*t551
      t554 = 3*t334
      t555 = 3*t333
      t556 = 2*t308
      t557 = 2*t307
      t580 = t370-t365+t554+t338-t555-t336-t335+t337-t557-t303+t556+t304
      t589 = 6*t547*t4*t551
      t601 = 6*t388
      t602 = 6*t392
      t603 = 6*t394
      t604 = 6*t396
      t605 = 9*t399
      t606 = 6*t401
      t607 = 6*t403
      t608 = 9*t405
      t609 = 9*t408
      t610 = 9*t410
      t611 = 9*t413
      t612 = -6*t384-6*t333+6*t334-t601+t602+t603-t604-t605+t606-t607+t6
     #08+t609+t610-t611
      t613 = 9*t415
      t617 = 9*t424
      t620 = 9*t434
      t621 = 3*t303
      t622 = 3*t304
      t623 = -t613+6*t335-6*t337+t370-t365+6*t421-t617-6*t428+6*t429+t62
     #0+t621-t622+t557-t556
      t636 = 5*t334
      t637 = 5*t333
      t645 = t621-t622-t345+t343-t370+t365-t554-t338+t555+t336+t335-t337
      t651 = t4*t549
      t669 = -8*t384-11*t333+11*t334-t601+t602+t603-t604-14*t399+t606-t6
     #07+16*t405+14*t408+14*t410-14*t413+t336
      t678 = -16*t415+8*t335-8*t337-t338+t419-t420+8*t421-16*t424-11*t42
     #8+11*t429+16*t434+t621-t622-t557+t556
      t683 = -t385-t555+t554-t601+t602+t603-t604-t605+t606-t607+t608+t60
     #9+t610-t611
      t686 = 3*t428
      t687 = 3*t429
      t690 = -t613+t417-t418+4*t197-4*t194+t422-t617-t686+t687+t620+6*t3
     #03-6*t304-t557+t556
      t699 = t549*t109
      t712 = l3*t314
      t728 = -t197-t308-t337+t335+t637+t194+7*t424-t401-t394-t392+7*t415
     #-t636
      t733 = t307+t403-5*t429+t396-t338+t388+t336+5*t428-7*t405-7*t434+t
     #421-t384
      t738 = t336+t537+t308-t375+t538-t539-t307+2*t421-t536+t374-2*t384-
     #t338
      t756 = -t384-t555+t554-3*t388+3*t392+3*t394-3*t396-7*t399+3*t401-3
     #*t403+8*t405+7*t408+7*t410-7*t413+t336
      t762 = -8*t415+t335-t337-t338+t419-t420+t421-8*t424-t686+t687+8*t4
     #34+t621-t622-5*t307+5*t308
      t778 = t621-t622+t556+2*t394-2*t396-t536-2*t388+2*t392+t538-t539+2
     #*t401-2*t403+3*t405-3*t415+t537-3*t424+3*t434-t557+t370-t365
      s3 = -((-10240.D0/3.D0*t3*t4*s2l2*t17+10240.D0/3.D0*t1*s2l2*l3*s23
     #r+t28*s2l2)*s1l1+(10240.D0/3.D0*t3*t4*s2l1*t17-10240.D0/3.D0*t24*s
     #2l1*s23r-t28*s2l1)*s1l2+(t50*s2l1+t58*s2l2)*l1*r*t16+((-10240.D0/3
     #.D0*t24*s3l2-10240.D0/3.D0*t65)*s2l1+(10240.D0/3.D0*t24*s3l1+10240
     #.D0/3.D0*t70)*s2l2)*s23r+(-t50*s34r+2048.D0/3.D0*t86*s3l2-10240.D0
     #/3.D0*t65)*s2l1+(-t58*s34r-2048.D0/3.D0*t86*s3l1+10240.D0/3.D0*t70
     #)*s2l2)*s13r-(((-8192.D0/3.D0*t102*s2l2-8192.D0/3.D0*t104)*s23r+10
     #240.D0/3.D0*t108*t110*s2l2)*l1*r*t16+(t122*s2l2+8192.D0/3.D0*t102*
     #s4l2*s34r+2048.D0/3.D0*t120*s4l2)*s23r+t135*s2l2)*s1l1
      s2 = s3-(((8192.D0/3.D0*t102*s2l1+8192.D0/3.D0*t140)*s23r-10240.D0
     #/3.D0*t108*t110*s2l1)*l1*r*t16+(-t122*s2l1-8192.D0/3.D0*t102*s4l1*
     #s34r-2048.D0/3.D0*t120*s4l1)*s23r-t135*s2l1)*s1l2-(-10240.D0/3.D0*
     #t164*t165*s2l1+10240.D0/3.D0*t164*t168*s2l2)*t173*t8*t182
      s3 = s2-(((6144*t184+8192.D0/3.D0*t104)*s2l1+(-6144*t189-8192.D0/3
     #.D0*t140)*s2l2-8192.D0/3.D0*t195+8192.D0/3.D0*t198)*s23r+((10240.D
     #0/3.D0*t164*t165-10240.D0/3.D0*t164*t47)*s34r-2048.D0/3.D0*t220*s3
     #l2)*s2l1+((-10240.D0/3.D0*t164*t168+10240.D0/3.D0*t164*t55)*s34r+2
     #048.D0/3.D0*t220*s3l1)*s2l2)*l1*r*t16-(((-8192.D0/3.D0*t184+2048.D
     #0/3.D0*t104)*s34r+2048.D0/3.D0*t243*s3l2-8192.D0/3.D0*t65)*s2l1+((
     #8192.D0/3.D0*t189-2048.D0/3.D0*t140)*s34r-2048.D0/3.D0*t243*s3l1+8
     #192.D0/3.D0*t70)*s2l2+(8192.D0/3.D0*t195-8192.D0/3.D0*t198)*s34r+2
     #048.D0/3.D0*t120*t194-2048.D0/3.D0*t120*t197)*s23r
      s1 = s3-(10240.D0/3.D0*t108*t109*s4l2*t131+(-2048.D0/3.D0*t276*s3l
     #2+10240.D0/3.D0*t104)*s34r+2048.D0/3.D0*t285*s3l2)*s2l1-(-10240.D0
     #/3.D0*t108*t109*s4l1*t131+(2048.D0/3.D0*t276*s3l1-10240.D0/3.D0*t1
     #40)*s34r-2048.D0/3.D0*t285*s3l1)*s2l2-(-256.D0/3.D0*t313*t314*l1*r
     #*t9*t11*t15+256.D0/3.D0*t313*t314*s34r+256.D0/3.D0*(t81*t305+2*t1*
     #t309)*t4*t314)*t331
      s3 = s1-(256.D0/3.D0*(2*(t333-t334+t335+t336-t337-t338)*t4*t164-t3
     #47*t4*t164)*t314*t173*t8*t176*t178*t181+(256.D0/3.D0*t359*t360*s23
     #r-256.D0/3.D0*(-t164*t4*(t363+t346+t364-t365-t344-t366-t367-t368+t
     #369+t370)+2*t164*t4*(-t308+t307+t303-t304-t366+t374+t363+t364-t367
     #-t375))*t314*s34r-256.D0/3.D0*(t407+t426+t451+t470)*l3*t474)*l1*r*
     #t16+(-256.D0/3.D0*t359*t481-256.D0/3.D0*t358*l3*t474)*s23r+256.D0/
     #3.D0*(-t164*t4*(t336-t334+t333-t194-t488+t335-t338-t337+t489+t197)
     #+2*t164*t4*(-t308+t307+t303-t304-t334-t338+t336+t333-t337+t335))*t
     #314*t131+256.D0/3.D0*(-t24*(t507+t508+t514+t515)+2*t24*(t523+t524)
     #)*t314*s34r+256.D0/3.D0*(t540+t541)*t314)*s13r
      s2 = s3+256.D0/3.D0*(t553-(t337-t554-t335+t555+t556-t557)*t4*t551)
     #*t314*t173*l1*t8*r*t176*t9/t177/t10/t180/t14-(-256.D0/3.D0*(-2*t30
     #9*t549*t109-t580*t549*t109)*t314*s23r+256.D0/3.D0*(t589-(-t335+8*t
     #333+t337-8*t334-t488+t489)*t4*t551)*t314*s34r+256.D0/3.D0*(t612+t6
     #23)*t101*t360)*t173*t8*t182-((256.D0/3.D0*(2*(t303-t304-t308+t337+
     #t307-t335)*t549*t109-(t370-t365+t636-t637-t557-t303+t556+t304)*t54
     #9*t109)*t314*s34r+512.D0/3.D0*t645*t101*t360)*s23r-256.D0/3.D0*(-t
     #651*t550*(t308+7*t333-t307-7*t334)+t589)*t314*t131-256.D0/3.D0*(t6
     #69+t678)*t101*t481-256.D0/3.D0*(t683+t690)*l3*t474)*l1*r*t16
      taux4 = s2-(-256.D0/3.D0*(-t699*(-t337+t366-t364-t338+t335+t336)-2
     #*t699*t547)*t314*t131-256.D0/3.D0*(t621-t622-t345+t375+t343-t374-t
     #370+t365-t636+t637)*t101*t481-256.D0/3.D0*t645*t1*t712)*s23r-256.D
     #0/3.D0*(-t651*t550*(2*t333-2*t334)+t553)*t314*t131*s34r-256.D0/3.D
     #0*(-t102*(t728+t733)+2*t102*t738)*t314*t131-256.D0/3.D0*(t756+t762
     #)*t1*t712*s34r-256.D0/3.D0*t778*t314

CCC Coefficients

      tcoff9 = taux3
      tcoff10 = taux3
      tcoff11 = taux3
      tcoff12 = taux3
      tcoff13 = taux4
      tcoff14 = taux4
      tcoff15 = -taux4
      tcoff16 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet10 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet11 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet12 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet13 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet14 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet15 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet16 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff16)

CCC Sum

      sq3num2e1 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*
     #tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tc
     #off16
      RETURN
      END

      double precision function sq3num2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = t1**2
      t3 = l3**2
      t4 = t3**2
      t5 = t2*t4
      t6 = -1+z
      t7 = t6**2
      t8 = t7**2
      t9 = t8*t7
      t11 = 1/s134r
      t14 = l1**2
      t15 = t14**2
      t16 = r**2
      t17 = t16**2
      t19 = 1+t16
      t20 = t19**2
      t21 = t20**2
      t22 = r+t
      t23 = t22**2
      t24 = t23**2
      t28 = 1+r*t
      t29 = t28**2
      t30 = t29**2
      t36 = z*l3
      t37 = -l3+t36+1
      t41 = t3*l3
      t42 = t41*t8
      t43 = t42*t1
      t44 = s1l2*s2l1
      t45 = 2*t36
      t46 = z**2
      t47 = l3*t46
      t48 = 2*t47
      t49 = s1l1*s2l2
      t50 = 2*t46
      t53 = s34r**2
      t54 = t53*s34r
      t58 = t7*t6
      t60 = t1*s4l2
      t62 = t3*t58*t60*s3l1
      t63 = t3*t1
      t66 = t63*t58*s4l1*s3l2
      t67 = z*t3
      t68 = t58*t1
      t70 = -2*z+2-2*l3+2*t36
      t72 = t67*t68*t70
      t73 = 2*l3
      t74 = 1+t45-t73
      t76 = t67*t68*t74
      t83 = t41*t58
      t84 = s23r**2
      t88 = t3*t7
      t89 = s3l2*s2l1
      t92 = s3l1*s2l2
      t95 = z*t1
      t99 = s4l2*s3l1
      t101 = 512.D0/3.D0*t88*t99
      t102 = s3l2*s4l1
      t104 = 512.D0/3.D0*t88*t102
      t105 = s4l1*s4l2
      t107 = 1024.D0/3.D0*t88*t105
      t108 = 3*z
      t111 = t7*z
      t114 = t67*t58
      t118 = t7*l3
      t130 = (-z+2-t73+t45)*l3
      t133 = t130*t6*s4l2*s3l1
      t137 = t130*t6*s4l1*s3l2
      t145 = t19/t22/t28
      t147 = t4*t8
      t148 = z*t11
      t152 = t8*t6
      t157 = s1l1*s4l2
      t158 = s1l2*s4l1
      t159 = t157+t158
      t166 = t5*t9
      t167 = t148*t53
      t171 = 2*t92
      t172 = 2*t89
      t173 = 18*t36
      t179 = t41*t11
      t180 = t179*s34r
      t184 = t46*z
      t185 = t3*t184
      t187 = t3*t46
      t189 = l3*t184
      t192 = 3*t92
      t194 = l3*s2l1
      t196 = t194*s3l2*z
      t197 = 3*t196
      t198 = l3*s2l2
      t200 = t198*s3l1*z
      t201 = 3*t200
      t202 = l3*s3l1
      t203 = t202*s2l2
      t204 = 3*t203
      t205 = l3*s3l2
      t206 = t205*s2l1
      t207 = 3*t206
      t208 = t89*z
      t209 = 3*t208
      t210 = t92*z
      t211 = 3*t210
      t212 = 3*t89
      t214 = 6*z
      t215 = -t99-12*t46+12*t185-24*t187-18*t189+12*t67-t192+6*t184-t197
     #-t201-t102+t204+t207+t209+t211-t212-t173+36*t47+t214
      t217 = t7*t11
      t226 = t20/t23/t29
      t228 = z*t58
      t237 = l1*r*t145
      t240 = s34r*s2l2
      t242 = s4l2*s34r
      t243 = t88*t242
      t244 = l3*t6
      t245 = t244*s4l2
      t254 = t6*s3l2
      t258 = s4l1*s34r
      t260 = t6*s3l1
      s1 = -128.D0/3.D0*t5*t9*z*t11*t15*t17*t21/t24/t30-512.D0/3.D0*t37*
     #s4l1*s3l2+64.D0/3.D0*t43*(-t44+t45-t48-t49+z+t50)*t11*t54+(-512.D0
     #/3.D0*t62-512.D0/3.D0*t66+512.D0/3.D0*t72-512.D0/3.D0*t76)*t53-512
     #.D0/3.D0*t37*s4l2*s3l1+(1024.D0/3.D0*t83*z*t84+(-512.D0/3.D0*t88*t
     #89-512.D0/3.D0*t88*t92+512*t42*t95*s34r+t101+t104+t107-512.D0/3.D0
     #*(t108-t45+t73-1)*t3*t111+512.D0/3.D0*t114)*s23r-512.D0/3.D0*t118*
     #t89-512.D0/3.D0*t118*t92+(1024.D0/3.D0*t62+1024.D0/3.D0*t66-512.D0
     #/3.D0*t72+512.D0/3.D0*t76)*s34r+512.D0/3.D0*t133+512.D0/3.D0*t137)
     #*l1*r*t145
      t273 = s1+(64.D0/3.D0*t147*t148*t84+(64.D0/3.D0*t95*t152*t4*t11*s3
     #4r+64.D0/3.D0*t159*t41*t58*t11)*s23r-128*t166*t167+64.D0/3.D0*(14*
     #t46-t49-t44+t171+t172+t173-18*t47-12*z)*t1*t8*t180-64.D0/3.D0*t215
     #*t3*t217)*t14*t16*t226+(64.D0/3.D0*t228*t180+64.D0/3.D0*t88*t148)*
     #t84+(512.D0/3.D0*t88*s4l2*s23r*t237+(-512.D0/3.D0*t88*t240-512.D0/
     #3.D0*t243-512.D0/3.D0*t245)*s23r-512.D0/3.D0*t118*t240)*s1l1+(-512
     #.D0/3.D0*t118*t242+512.D0/3.D0*t254)*s2l1+(-512.D0/3.D0*t118*t258+
     #512.D0/3.D0*t260)*s2l2+(-512.D0/3.D0*t133-512.D0/3.D0*t137+512.D0/
     #3.D0*t111*l3*(t67-1+t73-t3-t45)-512.D0/3.D0*t111*l3*(t67-t3+l3))*s
     #34r
      t274 = t148*s34r
      t277 = 6*t46
      t278 = 8*t36
      t286 = t14*l1
      t288 = t16*r
      t296 = t20*t19/t23/t22/t29/t28
      t298 = t1*t41
      t299 = t8*z
      t307 = s1l2*z
      t308 = t307*s2l1
      t309 = 2*t308
      t310 = t44*l3
      t311 = 2*t310
      t313 = s1l1*l3*s2l2
      t314 = 2*t313
      t315 = t49*z
      t316 = 2*t315
      t317 = 2*t44
      t318 = 5*t36
      t319 = 2*t49
      t320 = -t210+t309+t311+t314+t316-t208-t206-t203-t317+t318+t200+t19
     #6-t319
      t322 = t307*t194
      t323 = 2*t322
      t324 = t49*t36
      t325 = 2*t324
      t326 = 6*t67
      t327 = 10*t189
      t328 = 12*t187
      t329 = 6*t185
      t331 = 3*t46
      t332 = -15*t47-t323-t325+t99+t89+t92+t102-t326+t327+t328-t329+z-4*
     #t184+t331
      t345 = t152*t41
      t356 = t345*z
      t357 = t2*t11
      t358 = t357*s34r
      t361 = 4*t46
      t362 = 6*t36
      t363 = 6*t47
      t364 = 4*z
      t375 = t1*t11
      t379 = 2*t105
      t380 = 5*t46
      t381 = 4*t36
      t382 = 4*t47
      t389 = t63*t58
      t390 = 2*t99
      t391 = 2*t102
      t392 = s4l2*s2l1
      t393 = s4l1*s2l2
      t394 = -t390-t380+t49+t44-t192-t391+t392+t393-t212-t362+t363+t364
      t399 = 4*t210
      t400 = 4*t208
      t401 = 4*t206
      t402 = 4*t203
      t403 = 4*t200
      t404 = 4*t196
      t406 = t102*z
      t407 = t202*s4l2
      t408 = 2*t407
      t409 = l3*s4l1
      t410 = t409*s3l2
      t411 = 2*t410
      t412 = t99*z
      t413 = t399+t400+t401+t402-t278-t403-t404+16*t47+t406+t408+t411+t4
     #12
      t414 = z*s4l2
      t415 = t202*t414
      t416 = 2*t415
      t417 = z*s4l1
      t418 = t205*t417
      t419 = 2*t418
      t420 = 4*t89
      t421 = 4*t92
      t422 = 8*t189
      t423 = 2*z
      t424 = 2*t184
      t425 = -t416-t419-t390-t420-t421-t391+t326-t422-t328+t329+t423+t42
     #4-t361
      t428 = t6*t11
      t458 = t211-t308-t310-t313-t315+t209+t207+t204+t44-t381-t201-t197+
     #t49+9*t47+t322+t324+t406+t408+t411+t412
      t471 = -t192-t391+3*t67-5*t189-6*t187+3*t185+z+t424+t393-t331+t392
      t478 = 2*t185
      t479 = 4*t187
      t480 = 2*t189
      t481 = 2*t67
      t482 = 2*t196
      t483 = 2*t200
      t484 = 2*t203
      t485 = 2*t206
      t486 = 2*t208
      t487 = 2*t210
      t488 = -t418-t415-t99+t478-t479-t480+t481-t171-t482-t483-t102+t484
     #+t485+t486+t487-t172-t45+t382+t407+t410
      t504 = 512.D0/3.D0*t244+512.D0/3.D0*t118*t1*s34r
      t518 = t8*t2
      t521 = t14*t16*t226
      t527 = t1*s4l1
      t534 = l3*t1
      t536 = t534*t7*t46
      t540 = -t534*t111*t74
      t546 = t244*s3l2
      t550 = t244*s3l1
      t551 = t244*s4l1
      t585 = (-512.D0/3.D0*t118*t1*s2l2*t237+512.D0/3.D0*t6*s2l2*l3*s23r
     #+t504*s2l2)*s1l1+(-512.D0/3.D0*t118*t1*s2l1*t237+512.D0/3.D0*t244*
     #s2l1*s23r+t504*s2l1)*s1l2+512.D0/3.D0*t67*t518*t521+(-512.D0/3.D0*
     #t118*t60*s2l1-512.D0/3.D0*t118*t527*s2l2-1024.D0/3.D0*t67*t518*s34
     #r-512.D0/3.D0*t536+512.D0/3.D0*t540)*l1*r*t145+((512.D0/3.D0*t546+
     #512.D0/3.D0*t245)*s2l1+(512.D0/3.D0*t550+512.D0/3.D0*t551)*s2l2+10
     #24.D0/3.D0*t36*t6*t37-512.D0/3.D0*t36*t6*t70)*s23r+(512.D0/3.D0*t1
     #18*t60*s34r+512.D0/3.D0*t254+512.D0/3.D0*t245)*s2l1+(512.D0/3.D0*t
     #118*t527*s34r+512.D0/3.D0*t260+512.D0/3.D0*t551)*s2l2+512.D0/3.D0*
     #t67*t518*t53-512.D0/3.D0*(-t536+t540)*s34r+512.D0/3.D0*t47*t6+512.
     #D0/3.D0*t36*t6*t37
      t587 = -t99+t478-t479-t480+t481-t92-t196-t200-t102+t203+t206+t208+
     #t210-t89-t45+t382
      t593 = t88*t258
      t628 = 2*t157
      t629 = 2*t158
      t632 = t217*s34r
      t639 = -t487+t308+t310+t313+t315-t486-t485-t484-t44+t362+t483+t482
     #-t49
      t641 = -14*t47-t322-t324+t390+t172+t171+t391-t326+t422+t328-t329-t
     #424+t50
      t679 = -t399+t309+t311+t314+t316-t400-t401-t402-t317+24*t36+t403+t
     #404-t319
      t687 = -52*t47-t323-t325+t390+t420+t421+t391-18*t67+28*t189+36*t18
     #7-18*t185-t214-10*t184+16*t46
      t697 = -t390-t361+8*t185-16*t187-t327+8*t67-t192+t424-t197-t201-t3
     #91+t204+t207+t209+t211-t212-10*t36+20*t47+t423
      t709 = s34r*s2l1
      t726 = t299*t3
      t729 = z-t46-t89-t92-t45+t48
      t758 = s13r**2
      s2 = (128*t166*t274+64.D0/3.D0*(-t89-t92+t214-t277-t278+8*t47)*t1*
     #t42*t11)*t286*t288*t296+(-512.D0/3.D0*t298*t299*s23r-512.D0/3.D0*t
     #62-512.D0/3.D0*t66)*t14*t16*t226+64.D0/3.D0*t88*(t320+t332)*t11*t5
     #3
      s3 = s2+(-1024.D0/3.D0*t83*z*s34r-1024.D0/3.D0*t88*z)*t84
      s4 = s3
      s7 = 128.D0/3.D0*t2*z*t345*t11*t286*t288*t296+(128.D0/3.D0*t95*t8*
     #t179*s23r-64*t356*t358-64.D0/3.D0*(-t99-t361-t171-t102-t172-t362+t
     #363+t364)*t3*t68*t11)*t14*t16*t226+((-320.D0/3.D0*t42*z*t375*s34r-
     #64.D0/3.D0*(t379+t99+t157+t158-t380+t102-t381+t382+t108)*t3*t217)*
     #s23r+64.D0/3.D0*t389*t394*t11*s34r+64.D0/3.D0*(t413+t425)*l3*t428)
     #*l1*r*t145+(64*t43*t167+64.D0/3.D0*t88*(t99+t102+t379-t277+t364-t3
     #18+5*t47+t158+t157)*t11*s34r+64.D0/3.D0*(t379+t99+t157+t158-t331+t
     #102-t45+t48+z)*l3*t428)*s23r+64.D0/3.D0*t356*t357*t54+64.D0/3.D0*t
     #389*(t99+t46-t49-t44+t92+t102-t392-t393+t89)*t11*t53-64.D0/3.D0*t2
     #44*(t458-t416-t419+t194*t414+t198*t417-t392*z-t393*z-l3*s4l2*s2l1-
     #t409*s2l2-t390-t212+t471)*t11*s34r-64.D0/3.D0*t488*t11
      s8 = s13r
      s6 = s7*s8
      s7 = t585*s13r
      s5 = s6+s7
      s1 = s4+s5
      s2 = s1-64.D0/3.D0*t587*t11+((-512.D0/3.D0*t243+512.D0/3.D0*t546)*
     #s2l1+(-512.D0/3.D0*t593+512.D0/3.D0*t550)*s2l2-1024.D0/3.D0*t298*t
     #299*t53+(-t101-t104-t107-512.D0/3.D0*t88*z*(t423-1)+512.D0/3.D0*t8
     #8*z*(-2+3*l3-3*t36+t108))*s34r-512.D0/3.D0*t244*t99-512.D0/3.D0*t2
     #44*t102-1024.D0/3.D0*t244*t105+512.D0/3.D0*(-t36+l3+z+1)*z*t244)*s
     #23r+(64.D0/3.D0*t83*(t157+t158-t36+t47+t423-t46)*t11*t53+64.D0/3.D
     #0*(t628+t629-t46-t36+t47+t423)*t3*t632+64.D0/3.D0*t159*l3*t428)*s2
     #3r
      t760 = s2+64.D0/3.D0*(t639+t641)*l3*t428*s34r+((-64.D0/3.D0*t147*t
     #274-128.D0/3.D0*t228*t179)*t84+(-64.D0/3.D0*t4*t152*z*t375*t53-64.
     #D0/3.D0*(t628+t629-t50-t45+t48+t108)*t58*t180-128.D0/3.D0*t159*t3*
     #t217)*s23r+128.D0/3.D0*t166*t148*t54-64.D0/3.D0*t43*(10*t46-t319-t
     #317+t92+t89+12*t36-12*t47-5*z)*t11*t53-64.D0/3.D0*(t679+t687)*t3*t
     #632+64.D0/3.D0*t697*l3*t428)*l1*r*t145+(512.D0/3.D0*t88*s4l1*s23r*
     #t237+(-512.D0/3.D0*t88*t709-512.D0/3.D0*t593-512.D0/3.D0*t551)*s23
     #r-512.D0/3.D0*t118*t709)*s1l2+(-64.D0/3.D0*t299*t3*t2*t11*t521+(-1
     #28.D0/3.D0*t114*t375*s23r+128.D0/3.D0*t726*t358+64.D0/3.D0*t534*t7
     #*t729*t11)*l1*r*t145+(128.D0/3.D0*t389*t274+128.D0/3.D0*t118*t95*t
     #11)*s23r-64.D0/3.D0*t726*t357*t53-64.D0/3.D0*t534*t7*t729*t11*s34r
     #-64.D0/3.D0*t6*t1*(-t36+t47-t89-t92)*t11)*t758
      taux3 = t273+t760
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t8 = r**2
      t9 = 1+t8
      t10 = r+t
      t11 = 1/t10
      t14 = 1+r*t
      t15 = 1/t14
      t16 = t9*t11*t15
      t17 = l1*r*t16
      t24 = l3*t1
      t28 = -6144*t24-6144*t3*t4*s34r
      t43 = t4*s3l2
      t46 = t4*s4l2
      t49 = -26624.D0/3.D0*t3*t43-6144*t3*t46
      t51 = t4*s3l1
      t54 = t4*s4l1
      t57 = 26624.D0/3.D0*t3*t51+6144*t3*t54
      t64 = t24*s4l2
      t69 = t24*s4l1
      t77 = -t1
      t82 = 2*t1
      t83 = 2*l3
      t86 = -1+z+3*t77*(1-3*l3)+t82*(t83-1)
      t101 = l3**2
      t102 = t101*t2
      t105 = t102*s4l2
      t106 = 8192.D0/3.D0*t105
      t109 = t101*t4
      t110 = t2*t1
      t111 = t110*s34r
      t124 = 8192.D0/3.D0*t102*s34r+10240.D0/3.D0*t77*l3
      t130 = 4*t77*l3
      t135 = t101*t110
      t136 = s34r**2
      t141 = t3+8*t102
      t144 = 16384.D0/3.D0*t135*t4*t136+2048.D0/3.D0*t141*s34r
      t150 = t102*s4l1
      t151 = 8192.D0/3.D0*t150
      t181 = t135*t46*s3l1
      t185 = t109*t110*s4l1*s3l2
      t188 = l1**2
      t191 = t9**2
      t192 = t10**2
      t193 = 1/t192
      t195 = t14**2
      t196 = 1/t195
      t197 = t191*t193*t196
      t199 = t102*s3l2
      t204 = t102*s3l1
      t209 = s4l2*s3l1
      t210 = t102*t209
      t212 = s3l2*s4l1
      t213 = t102*t212
      t222 = 2*z
      t223 = z*l3
      t224 = 2*t223
      t229 = 3*z
      t236 = 2*(1-t222-t83+t224)*t1*l3-3*(-2+t229+4*l3-4*t223)*t1*l3+t3
      t250 = -t181+t185
      t254 = (-z+2-t83+t224)*l3
      t257 = t254*t1*s4l2*s3l1
      t261 = t254*t1*s4l1*s3l2
      t271 = 14*l3*t1
      t303 = -3*t24*(-t229+t224+2-t83)+2*t24*(t222-t223-1+l3)
      t312 = -8*l3+8*t223-4*z-1
      t335 = -l3+t223+1
      s2 = -((6144*t3*t4*s2l2*t17-6144*t1*s2l2*l3*s23r+t28*s2l2)*s1l1+(-
     #6144*t3*t4*s2l1*t17+6144*t24*s2l1*s23r-t28*s2l1)*s1l2+(t49*s2l1+t5
     #7*s2l2)*l1*r*t16+((6144*t24*s3l2+6144*t64)*s2l1+(-6144*t24*s3l1-61
     #44*t69)*s2l2)*s23r+(-t49*s34r+2048.D0/3.D0*t86*s3l2+6144*t64)*s2l1
     #+(-t57*s34r-2048.D0/3.D0*t86*s3l1-6144*t69)*s2l2)*s13r-(((10240.D0
     #/3.D0*t102*s2l2+t106)*s23r-16384.D0/3.D0*t109*t111*s2l2)*l1*r*t16+
     #(t124*s2l2-8192.D0/3.D0*t102*s4l2*s34r+2048.D0/3.D0*t130*s4l2)*s23
     #r+t144*s2l2)*s1l1
      s1 = s2-(((-10240.D0/3.D0*t102*s2l1-t151)*s23r+16384.D0/3.D0*t109*
     #t111*s2l1)*l1*r*t16+(-t124*s2l1+8192.D0/3.D0*t102*s4l1*s34r-2048.D
     #0/3.D0*t130*s4l1)*s23r-t144*s2l1)*s1l2-(16384.D0/3.D0*t135*t43*s2l
     #1-16384.D0/3.D0*t135*t51*s2l2+2048.D0/3.D0*t181-2048.D0/3.D0*t185)
     #*t188*t8*t197-(((-28672.D0/3.D0*t199-10240.D0/3.D0*t105)*s2l1+(286
     #72.D0/3.D0*t204+10240.D0/3.D0*t150)*s2l2+8192.D0/3.D0*t210-8192.D0
     #/3.D0*t213)*s23r+((-16384.D0/3.D0*t135*t43+16384.D0/3.D0*t135*t46)
     #*s34r-2048.D0/3.D0*t236*s3l2)*s2l1+((16384.D0/3.D0*t135*t51-16384.
     #D0/3.D0*t135*t54)*s34r+2048.D0/3.D0*t236*s3l1)*s2l2+4096.D0/3.D0*t
     #250*s34r-2048.D0/3.D0*t257+2048.D0/3.D0*t261)*l1*r*t16
      t339 = s1-(((10240.D0/3.D0*t199-t106)*s34r+2048.D0/3.D0*t271*s3l2+
     #10240.D0/3.D0*t64)*s2l1+((-10240.D0/3.D0*t204+t151)*s34r-2048.D0/3
     #.D0*t271*s3l1-10240.D0/3.D0*t69)*s2l2+(-8192.D0/3.D0*t210+8192.D0/
     #3.D0*t213)*s34r+2048.D0/3.D0*t130*t209-2048.D0/3.D0*t130*t212)*s23
     #r-(-16384.D0/3.D0*t135*t46*t136+(-2048.D0/3.D0*t303*s3l2-2048.D0/3
     #.D0*t141*s4l2)*s34r+2048.D0/3.D0*t312*s3l2)*s2l1-(16384.D0/3.D0*t1
     #35*t54*t136+(2048.D0/3.D0*t303*s3l1+2048.D0/3.D0*t141*s4l1)*s34r-2
     #048.D0/3.D0*t312*s3l1)*s2l2+2048.D0/3.D0*t250*t136-(2048.D0/3.D0*t
     #257-2048.D0/3.D0*t261)*s34r-2048.D0/3.D0*t335*s4l2*s3l1
      t343 = s4l2*s2l1
      t344 = s4l1*s2l2
      t345 = s3l2*s2l1
      t346 = 2*t345
      t347 = s3l1*s2l2
      t348 = 2*t347
      t349 = t343-t344+t346-t348
      t351 = -t344+t343
      t355 = (-t3*t349+2*t3*t351)*t4
      t356 = 1/s134r
      t372 = s13r**2
      t374 = 2*t343
      t375 = s1l1*s2l2
      t376 = 4*t345
      t377 = s1l2*s2l1
      t378 = s1l1*s4l2
      t379 = 4*t347
      t380 = 2*t344
      t381 = s3l2*s1l1
      t382 = s3l1*s1l2
      t383 = s1l2*s4l1
      t384 = t212+t374-t375+t376-t209+t377+t378-t379-t380+t381-t382-t383
      t399 = 3*t377
      t400 = 2*t382
      t402 = 3*t212
      t403 = 3*t375
      t404 = 3*t209
      t405 = 2*t381
      t407 = t399+t400+5*t383+t402-t403-t404-t405-5*t378
      t408 = t407*t101
      t409 = t2*t356
      t412 = 2*t383
      t413 = 3*t347
      t414 = 2*t378
      t415 = 2*t212
      t416 = 2*t209
      t417 = 3*t345
      t421 = 7*t347
      t422 = 7*t345
      t430 = t375*t223
      t432 = s1l2*z
      t433 = l3*s2l1
      t434 = t432*t433
      t436 = t378*t223
      t437 = 2*t436
      t438 = t432*s2l1
      t439 = 2*t438
      t440 = l3*s4l1
      t441 = t440*s2l2
      t442 = 4*t441
      t444 = l3*s4l2*s2l1
      t445 = 4*t444
      t446 = t344*z
      t447 = 2*t446
      t448 = t383*t223
      t449 = 2*t448
      t450 = z*s4l2
      t451 = t433*t450
      t452 = 4*t451
      t453 = t405-t400-t383-4*t430+4*t434+t437-t439+t442-t445+t447-t449+
     #t452
      t454 = l3*s2l2
      t455 = z*s4l1
      t456 = t454*t455
      t457 = 4*t456
      t458 = l3*s3l2
      t459 = t458*t455
      t460 = 2*t459
      t461 = t343*z
      t462 = 2*t461
      t463 = t209*z
      t464 = l3*s3l1
      t465 = t464*t450
      t466 = 2*t465
      t467 = t440*s3l2
      t468 = 2*t467
      t469 = t464*s4l2
      t470 = 2*t469
      t475 = t381*z
      t477 = -t457+t460-t462+t463-t466-t468+t470+4*t343-4*t344+12*t345-1
     #2*t347-2*t475
      t480 = s1l1*l3
      t481 = t480*s3l2
      t483 = t382*z
      t485 = s1l2*l3
      t486 = t485*s3l1
      t488 = z*s3l1
      t489 = t485*t488
      t491 = z*s3l2
      t492 = t480*t491
      t494 = 2*t375
      t495 = t383*z
      t496 = t383*l3
      t497 = 2*t496
      t498 = t378*z
      t499 = t378*l3
      t500 = 2*t499
      t502 = 4*t212-4*t481+2*t483+4*t486-4*t489+4*t492-t494+t495+t497-t4
     #98-t500-4*t209
      t503 = 2*t377
      t504 = t212*z
      t505 = t375*z
      t506 = 2*t505
      t507 = t458*s2l1
      t509 = t464*s2l2
      t511 = t433*t491
      t513 = t454*t488
      t515 = t480*s2l2
      t517 = t345*z
      t519 = t347*z
      t521 = t377*l3
      t523 = t503+t378-t504+t506-12*t507+12*t509+12*t511-12*t513+4*t515-
     #7*t517+7*t519-4*t521
      t527 = t1*t356
      t534 = t409*s34r
      t552 = -t383+t430-t434+t437+t438+t441-t444+t446-t449+t451-t456+t46
     #0-t461+t463-t466-t468+t470+t343-t344+t417
      t553 = 3*t507
      t554 = 3*t509
      t555 = 3*t511
      t556 = 3*t513
      t559 = -t413+t415+t375+t495+t497-t498-t500-t416-t377+t378-t504-t50
     #5-t553+t554+t555-t556-t515-4*t517+4*t519+t521
      t562 = 2*t441
      t563 = 2*t444
      t564 = 2*t451
      t565 = 2*t456
      t566 = t381-t382-t383+t437+t438+t562-t563+t446-t449+t564-t565
      t567 = t460-t461+t463-t466-t468+t470+t374-t380+t422-t421-t475
      t569 = 2*t481
      t570 = 2*t486
      t571 = 2*t489
      t572 = 2*t492
      t573 = t402-t569+t483+t570-t571+t572+t375+t495+t497-t498-t500
      t580 = -t404-t377+t378-t504-t505-7*t507+7*t509+7*t511-7*t513-5*t51
     #7+5*t519
      t589 = 2*t430
      t590 = 2*t434
      t593 = -t589+t590+t436+t562-t563-t448+t564-t565+t459-t465-t467+t46
     #9+t374-t380+6*t345-6*t347
      t594 = 6*t507
      t595 = 6*t509
      t596 = 6*t511
      t597 = 6*t513
      t598 = 2*t515
      t599 = 2*t521
      t600 = t402-t569+t570-t571+t572+t496-t499-t404-t594+t595+t596-t597
     #+t598-t517+t519-t599
      t606 = 3*t382
      t607 = 3*t381
      t610 = t101*l3
      t611 = t2**2
      t612 = t610*t611
      t633 = t415-t416+t606-t607-t378+t377+t383-t375+t346-t344+t343-t348
      t649 = t375-t405+t414-t377+t400-t412-t413-t380+t416-t415+t417+t374
      t657 = 6*t430
      t658 = 6*t434
      t659 = 3*t436
      t661 = 3*t441
      t662 = 3*t444
      t663 = 3*t446
      t664 = 3*t448
      t665 = 14*t381-14*t382+t412-t657+t658-t659-4*t438-t661+t662-t663+t
     #664
      t666 = 3*t451
      t667 = 3*t456
      t668 = 3*t459
      t669 = 3*t461
      t671 = 3*t465
      t672 = 3*t467
      t673 = 3*t469
      t674 = -t666+t667+t668+t669+2*t463-t671-t672+t673-t343+t344+t376-t
     #379
      t677 = 21*t481
      t679 = 21*t486
      t680 = 21*t489
      t681 = 21*t492
      t684 = 3*t496
      t686 = 3*t499
      t687 = -14*t475-t677+14*t483+t679-t680+t681-4*t375-2*t495-t684+2*t
     #498+t686
      t691 = 9*t507
      t692 = 9*t509
      t693 = 9*t511
      t694 = 9*t513
      t695 = 6*t515
      t696 = 9*t517
      t697 = 9*t519
      t698 = 6*t521
      t699 = 4*t377-t414-2*t504+4*t505-t691+t692+t693-t694+t695-t696+t69
     #7-t698
      t709 = 5*t381
      t710 = 5*t382
      t715 = 3*t378
      t716 = 3*t383
      t724 = 7*t382
      t725 = 7*t381
      t726 = 5*t345
      t727 = 5*t347
      t728 = t402-t404+t724-t725-t715+t399+t716-t403+t726-t380+t374-t727
      t734 = t4*t610
      t735 = -t381+t343+t382-t344+t346-t399+t209+t403-t212-t383-t348+t37
     #8
      t751 = 25*t381-25*t382+t383-t657+t658-4*t436-t439-t442+t445-4*t446
     #+4*t448-t452
      t759 = t457+4*t459+4*t461+3*t463-4*t465-4*t467+4*t469-t374+t380-t3
     #76+t379-25*t475
      t770 = -t402-36*t481+25*t483+36*t486-36*t489+36*t492-t494-3*t495-4
     #*t496+3*t498+4*t499+t404
      t774 = t503-t378-3*t504+t506-t594+t595+t596-t597+t695-6*t517+6*t51
     #9-t698
      t780 = t725-t724+t383-t657+t658-t659-t439-t661+t662-t663+t664-t666
      t781 = 7*t475
      t782 = t667+t668+t669+t463-t671-t672+t673+t343-t344-t345+t347-t781
      t784 = 7*t483
      t785 = -t402-t677+t784+t679-t680+t681-t494-t495-t684+t498+t686+t40
     #4
      t786 = t503-t378-t504+t506-t691+t692+t693-t694+t695-t696+t697-t698
      t809 = -11*t381+t402-t404+11*t382+t726-t727-t380+t374-t715+t716-t3
     #75+t377
      t813 = l3*t356
      t832 = 3*t430
      t833 = 3*t434
      t835 = -t381+t382-t383+t832-t833+t436+3*t438+t441-t444+t446-t448+t
     #451
      t836 = -t456-t459-t461-t463+t465+t467-t469+t343-t344+t417-t413+t47
     #5
      t838 = t212+t481-t483-t486+t489-t492+t403+t495+t496-t498-t499-t209
      t844 = 3*t515
      t847 = 3*t521
      t848 = -t399+t378+t504-3*t505-2*t507+2*t509+2*t511-2*t513-t844-2*t
     #517+2*t519+t847
      t854 = 7*t481
      t855 = 5*t483+t709+t519-5*t475-t507+t509-t212+t438-t854+t347-t710-
     #t383
      t856 = 7*t492
      t857 = 7*t486
      t858 = 7*t489
      t859 = t856+t378-t345+t375-t517-t505+t209+t857-t377+t511-t513-t858
      t867 = t725-t724-t383-t832+t833-t437+t438-t562+t563-t447+t449
      t869 = -t564+t565+t460+t462+t463-t466-t468+t470-t422+t421-t781-5*t
     #212
      t876 = -18*t481+t784+18*t486-18*t489+18*t492+t375-t495-t497+t498+t
     #500+5*t209
      t877 = 3*t517
      t878 = 3*t519
      t879 = -t377+t378-t504-t505-t553+t554+t555-t556+t844-t877+t878-t84
     #7
      t886 = -t589+t590-t436-t441+t444-t446+t448-t451+t456+t459+t461-t46
     #5-t467+t469+t343-t344-t346
      t887 = t348-t415-t854+t857-t858+t856-t496+t499+t416-t553+t554+t555
     #-t556+t598-t877+t878-t599
      s2 = 2048.D0/3.D0*t335*s4l1*s3l2-(-256.D0/3.D0*t355*t356*l1*r*t9*t
     #11*t15+256.D0/3.D0*t355*t356*s34r+256.D0/3.D0*(t77*t349+t82*t351)*
     #t4*t356)*t372
      s3 = s2-(256.D0/3.D0*(2*t384*t4*t135-(-t209+t378-t383-t348+t212+t3
     #74-t380+t346)*t4*t135)*t356*t188*t8*t191*t193*t196+(256.D0/3.D0*t4
     #08*t409*s23r-256.D0/3.D0*(-t135*t4*(-t412-t413+t414+t343-t377-t344
     #+t415+t375-t416+t417)+2*t135*t4*(t414+t374+t405-t416-t380-t400-t41
     #2-t421+t422+t415))*t356*s34r-256.D0/3.D0*(t453+t477+t502+t523)*l3*
     #t527)*l1*r*t16+(-256.D0/3.D0*t408*t534-256.D0/3.D0*t407*l3*t527)*s
     #23r+256.D0/3.D0*(-t135*t4*(t378-t377-t347+t212+t375-t343-t209+t345
     #+t344-t383)+2*t135*t4*(t378-t382+t381-t209-t413+t375-t383-t377+t41
     #7+t212))*t356*t136+256.D0/3.D0*(-t24*(t552+t559)+2*t24*(t566+t567+
     #t573+t580))*t356*s34r+256.D0/3.D0*(t593+t600)*t356)*s13r
      s1 = s3+256.D0/3.D0*(2*(t377-t606+t607-t375-t348+t346)*t4*t612-(t3
     #82-t383-t381+t378+t345+t343-t347-t344+t209-t212)*t4*t612)*t356*t18
     #8*l1*t8*r*t191*t9/t192/t10/t195/t14-(-256.D0/3.D0*(2*t633*t610*t11
     #0-(-t209+t378-t383+t375-t377-t382+t381+t347+t212-t345)*t610*t110)*
     #t356*s23r+256.D0/3.D0*(2*(-8*t382+8*t381-t413+t377-t375+t417)*t4*t
     #612-t649*t4*t612)*t356*s34r+256.D0/3.D0*(t665+t674+t687+t699)*t101
     #*t409)*t188*t8*t197
      s2 = s1-((256.D0/3.D0*(2*(-t709+t415-t416+t710+t346-t348-t344+t343
     #)*t610*t110-(-t209+t715-t716+t375-t377-t382+t381+t347+t212-t345)*t
     #610*t110)*t356*s34r+512.D0/3.D0*t728*t101*t409)*s23r-256.D0/3.D0*(
     #-t734*t611*t735+2*t734*t611*(t345+t725-t347-t724))*t356*t136-256.D
     #0/3.D0*(t751+t759+t770+t774)*t101*t534-256.D0/3.D0*(t780+t782+t785
     #+t786)*l3*t527)*l1*r*t16-(-256.D0/3.D0*(-(2*t378-2*t383)*t610*t110
     #+2*t610*t110*(-t377+t400-t405-t383+t375+t378))*t356*t136-256.D0/3.
     #D0*t809*t101*t534-256.D0/3.D0*t728*t1*t813)*s23r
      t891 = s2-256.D0/3.D0*(-t734*t611*(2*t375-2*t377)+2*t734*t611*(2*t
     #381-2*t382))*t356*t136*s34r-256.D0/3.D0*(-t102*(t835+t836+t838+t84
     #8)+2*t102*(t855+t859))*t356*t136-256.D0/3.D0*(t867+t869+t876+t879)
     #*t1*t813*s34r-256.D0/3.D0*(t886+t887)*t356
      taux4 = t339+t891

CCC Coefficients

      tcoff9 = taux3
      tcoff10 = taux3
      tcoff11 = taux3
      tcoff12 = taux3
      tcoff13 = taux4
      tcoff14 = taux4
      tcoff15 = -taux4
      tcoff16 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet10 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet11 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet12 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet13 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet14 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet15 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet16 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff16)

CCC Sum

      sq3num2e2 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*
     #tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tc
     #off16
      RETURN
      END

      double precision function sq3num2e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3**2
      t2 = t1**2
      t3 = -1+z
      t4 = t3**2
      t5 = t4**2
      t6 = t5*t3
      t9 = l3-1
      t10 = 1/s134r
      t11 = t9*t10
      t12 = s34r**2
      t15 = z**2
      t16 = l3*t15
      t17 = z*l3
      t18 = t16-t17-t15+z
      t20 = t4*t3
      t22 = t1*l3
      t23 = t22*t10
      t26 = 3*t18
      t28 = t4*t10
      t33 = t22*t5
      t34 = t33*t9
      t35 = z*t10
      t36 = t35*t12
      t38 = 128.D0/3.D0*t34*t36
      t39 = 2*t18
      t40 = t39*t1
      t42 = t40*t28*s34r
      t47 = r**2
      t48 = 1+t47
      t49 = r+t
      t53 = 1+r*t
      t55 = t48/t49/t53
      t57 = t12*s34r
      t58 = t35*t57
      t61 = t1*t4
      t66 = t5*z
      t67 = t9**2
      t71 = l1**2
      t73 = t48**2
      t74 = t49**2
      t77 = t53**2
      t79 = t73/t74/t77
      t80 = t71*t47*t79
      t83 = t66*t1
      t84 = t67*t10
      t85 = t84*s34r
      t88 = l3*t9
      t89 = 2*t17
      t90 = 2*t16
      t91 = z-t15-t89+t90
      t114 = s13r**2
      t117 = t6*z
      t120 = t71*l1
      t122 = t47*r
      t130 = t73*t48/t74/t49/t77/t53
      t134 = t67*z
      t135 = t6*t22
      t142 = t9*z
      t148 = 3*t16
      t150 = 2*t15
      t151 = 2*z
      t154 = t9*t20
      t164 = t33*z*t11*s34r
      t170 = t135*t67
      t174 = t1*t9*t20
      t180 = t15*z
      t183 = t1*t180
      t185 = t1*t15
      t187 = z*t1
      t189 = l3*t180
      t193 = t3*t10
      t204 = t18*l3
      t205 = t204*t193
      t217 = 2*t185
      t233 = t22*t20
      t234 = s23r**2
      t245 = t187*t20
      t248 = t9*(2*t17-2*z+2-2*l3)
      t280 = t4*z
      t287 = t4*t15
      s2 = ((-64.D0/3.D0*t2*t6*z*t11*t12-256.D0/3.D0*t18*t20*t23*s34r-64
     #.D0/3.D0*t26*t1*t28)*s23r-t38-64.D0/3.D0*t42)*l1*r*t55+64.D0/3.D0*
     #t34*t58+64.D0/3.D0*t61*t39*t10*t12+(-64.D0/3.D0*t66*t1*t67*t10*t80
     #+(128.D0/3.D0*t83*t85+64.D0/3.D0*t88*t4*t91*t10)*l1*r*t55-64.D0/3.
     #D0*t83*t84*t12-64.D0/3.D0*t88*t4*t91*t10*s34r-64.D0/3.D0*t3*t9*(-t
     #17+t16)*t10)*t114
      s1 = s2-64.D0/3.D0*t9*t2*t117*t10*s23r*t120*t122*t130+(64.D0/3.D0*
     #t134*t135*t10*t120*t122*t130+(64.D0/3.D0*t142*t5*t23*s23r-64.D0/3.
     #D0*t135*z*t85-64.D0/3.D0*(t148-3*t17-t150+t151)*t1*t154*t10)*t71*t
     #47*t79+((-64*t164-64.D0/3.D0*t40*t28)*s23r-64.D0/3.D0*t170*t36+64.
     #D0/3.D0*t174*(t90-t89-t150+z)*t10*s34r+64.D0/3.D0*(t180+8*t16-4*t1
     #7-t150+z+3*t183-6*t185+3*t187-4*t189)*l3*t193)*l1*r*t55+(t38+64.D0
     #/3.D0*t61*t26*t10*s34r+64.D0/3.D0*t205)*s23r+64.D0/3.D0*t170*t58+6
     #4.D0/3.D0*t174*(t16-t17+z)*t10*t12-64.D0/3.D0*l3*t3*(t180+t148-t17
     #-t15+t183-t217+t187-2*t189)*t10*s34r-64.D0/3.D0*(t90-t17+t183-t217
     #+t187-t189)*t10)*s13r-512.D0/3.D0*t135*t134*s34r*t80+(-512.D0/3.D0
     #*t233*z*t234-512.D0/3.D0*t20*t1*z*s23r+1024.D0/3.D0*t135*t134*t12+
     #512.D0/3.D0*t245*t248*s34r)*l1*r*t55
      taux3 = s1-512.D0/3.D0*t22*t67*t117*t57+((128.D0/3.D0*t142*t6*t2*t
     #10*s34r+64.D0/3.D0*t26*t22*t20*t10)*s23r+64.D0/3.D0*t164)*t71*t47*
     #t79+(512.D0/3.D0*t61*z*(t151-1)*s34r+512.D0/3.D0*t280*l3)*s23r+((5
     #12.D0/3.D0*t187*t154*s23r+512.D0/3.D0*t88*t287)*l1*r*t55+(-512.D0/
     #3.D0*t187*t154*s34r-512.D0/3.D0*t17*t3*(-l3+t17+1))*s23r-512.D0/3.
     #D0*t88*t287*s34r-512.D0/3.D0*t16*t3)*s13r-512.D0/3.D0*t17*t4*(t187
     #-t89-1+2*l3-t1)*s34r+(512.D0/3.D0*t233*z*s34r+512.D0/3.D0*t280*t1)
     #*t234+64.D0/3.D0*t204*t193*s34r-512.D0/3.D0*t245*t248*t12+(64.D0/3
     #.D0*t233*t18*t10*t12+64.D0/3.D0*t42+64.D0/3.D0*t205)*s23r
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t8 = r**2
      t9 = 1+t8
      t10 = r+t
      t11 = 1/t10
      t14 = 1+r*t
      t15 = 1/t14
      t16 = t9*t11*t15
      t17 = l1*r*t16
      t24 = l3*t1
      t25 = t4*s34r
      t28 = 14336.D0/3.D0*t24+14336.D0/3.D0*t3*t25
      t43 = t4*s3l2
      t49 = 4096*t3*t43+14336.D0/3.D0*t3*t4*s4l2
      t51 = t4*s3l1
      t57 = -4096*t3*t51-14336.D0/3.D0*t3*t4*s4l1
      t64 = t24*s4l2
      t69 = t24*s4l1
      t77 = 3*z
      t78 = 2*t1
      t82 = 3-t77+t78*(1-3*l3)
      t97 = l3**2
      t98 = t97*t2
      t101 = t98*s4l2
      t105 = t2*t1
      t106 = t97*t105
      t117 = -10240.D0/3.D0*t98*s34r+4096.D0/3.D0*t24
      t122 = -t1
      t123 = t122*l3
      t128 = t97*t4
      t129 = s34r**2
      t133 = 3*t3
      t135 = -t133-4*t98
      t138 = -8192.D0/3.D0*t128*t105*t129+2048.D0/3.D0*t135*s34r
      t144 = t98*s4l1
      t174 = t105*s4l2
      t176 = t128*t174*s3l1
      t178 = t105*s4l1
      t180 = t128*t178*s3l2
      t183 = l1**2
      t186 = t9**2
      t187 = t10**2
      t188 = 1/t187
      t190 = t14**2
      t191 = 1/t190
      t192 = t186*t188*t191
      t194 = t98*s3l2
      t199 = t98*s3l1
      t204 = s4l2*s3l1
      t205 = t98*t204
      t207 = s3l2*s4l1
      t208 = t98*t207
      t217 = 4*l3
      t218 = z*l3
      t219 = 4*t218
      t224 = 2*(-2+t77+t217-t219)*t1*l3-t133
      t238 = t176-t180
      t241 = 2*l3
      t242 = 2*t218
      t244 = (-z+2-t241+t242)*t1
      t245 = l3*s3l1
      t246 = t245*s4l2
      t247 = t244*t246
      t249 = l3*s4l1
      t250 = t249*s3l2
      t251 = t244*t250
      t262 = 9*t122*l3
      t289 = -t77+t242+2-t241
      t297 = t217-t219-z+3
      t323 = 3*l3-3*t218-3
      s2 = -((-14336.D0/3.D0*t3*t4*s2l2*t17+14336.D0/3.D0*t1*s2l2*l3*s23
     #r+t28*s2l2)*s1l1+(14336.D0/3.D0*t3*t4*s2l1*t17-14336.D0/3.D0*t24*s
     #2l1*s23r-t28*s2l1)*s1l2+(t49*s2l1+t57*s2l2)*l1*r*t16+((-14336.D0/3
     #.D0*t24*s3l2-14336.D0/3.D0*t64)*s2l1+(14336.D0/3.D0*t24*s3l1+14336
     #.D0/3.D0*t69)*s2l2)*s23r+(-t49*s34r+2048.D0/3.D0*t82*s3l2-14336.D0
     #/3.D0*t64)*s2l1+(-t57*s34r-2048.D0/3.D0*t82*s3l1+14336.D0/3.D0*t69
     #)*s2l2)*s13r-(((-4096.D0/3.D0*t98*s2l2+2048.D0/3.D0*t101)*s23r+819
     #2.D0/3.D0*t106*t25*s2l2)*l1*r*t16+(t117*s2l2-2048.D0/3.D0*t98*s4l2
     #*s34r+2048.D0/3.D0*t123*s4l2)*s23r+t138*s2l2)*s1l1
      s1 = s2-(((4096.D0/3.D0*t98*s2l1-2048.D0/3.D0*t144)*s23r-8192.D0/3
     #.D0*t106*t25*s2l1)*l1*r*t16+(-t117*s2l1+2048.D0/3.D0*t98*s4l1*s34r
     #-2048.D0/3.D0*t123*s4l1)*s23r-t138*s2l1)*s1l2-(-8192.D0/3.D0*t106*
     #t43*s2l1+8192.D0/3.D0*t106*t51*s2l2-2048*t176+2048*t180)*t183*t8*t
     #192-(((6144*t194+4096.D0/3.D0*t101)*s2l1+(-6144*t199-4096.D0/3.D0*
     #t144)*s2l2+2048.D0/3.D0*t205-2048.D0/3.D0*t208)*s23r+((8192.D0/3.D
     #0*t106*t43-8192.D0/3.D0*t128*t174)*s34r-2048.D0/3.D0*t224*s3l2)*s2
     #l1+((-8192.D0/3.D0*t106*t51+8192.D0/3.D0*t128*t178)*s34r+2048.D0/3
     #.D0*t224*s3l1)*s2l2+4096*t238*s34r+2048*t247-2048*t251)*l1*r*t16
      t327 = s1-(((-4096.D0/3.D0*t194+10240.D0/3.D0*t101)*s34r+2048.D0/3
     #.D0*t262*s3l2-4096.D0/3.D0*t64)*s2l1+((4096.D0/3.D0*t199-10240.D0/
     #3.D0*t144)*s34r-2048.D0/3.D0*t262*s3l1+4096.D0/3.D0*t69)*s2l2+(-20
     #48.D0/3.D0*t205+2048.D0/3.D0*t208)*s34r+2048.D0/3.D0*t123*t204-204
     #8.D0/3.D0*t123*t207)*s23r-(8192.D0/3.D0*t128*t174*t129+(-4096.D0/3
     #.D0*t24*t289*s3l2-2048.D0/3.D0*t135*s4l2)*s34r+2048.D0/3.D0*t297*s
     #3l2)*s2l1-(-8192.D0/3.D0*t128*t178*t129+(4096.D0/3.D0*t24*t289*s3l
     #1+2048.D0/3.D0*t135*s4l1)*s34r-2048.D0/3.D0*t297*s3l1)*s2l2+2048*t
     #238*t129-(-2048*t247+2048*t251)*s34r-2048.D0/3.D0*t323*s4l2*s3l1
      t331 = s3l2*s2l1
      t332 = s3l1*s2l2
      t333 = t331-t332
      t335 = s4l2*s2l1
      t336 = s4l1*s2l2
      t337 = 2*t331
      t338 = 2*t332
      t339 = t335-t336+t337-t338
      t343 = (-t3*t333+2*t3*t339)*t4
      t344 = 1/s134r
      t360 = s13r**2
      t362 = s1l2*s4l1
      t363 = 2*t336
      t364 = s1l1*s4l2
      t365 = 2*t335
      t381 = s1l2*s2l1
      t382 = 2*t381
      t383 = s1l1*s2l2
      t384 = 2*t383
      t385 = -t362-t382+t204+t384-t207+t364
      t386 = t385*t97
      t387 = t2*t344
      t390 = 3*t331
      t391 = 2*t207
      t392 = 2*t204
      t393 = 3*t332
      t397 = 2*t362
      t398 = 2*t364
      t406 = t362*z
      t407 = 2*t406
      t408 = t362*l3
      t409 = 4*t408
      t410 = t364*z
      t411 = 2*t410
      t412 = t364*l3
      t413 = 4*t412
      t414 = l3*s3l2
      t415 = t414*s2l1
      t416 = 4*t415
      t417 = t245*s2l2
      t418 = 4*t417
      t419 = t331*z
      t421 = t332*z
      t423 = t249*s2l2
      t424 = 8*t423
      t426 = l3*s4l2*s2l1
      t427 = 8*t426
      t428 = t336*z
      t429 = 8*t428
      t430 = t335*z
      t431 = 8*t430
      t432 = t407+t409-t392+t391-t397+t398-t411-t413-t416+t418-6*t419+6*
     #t421+t424-t427+t429-t431
      t433 = l3*s2l1
      t434 = s3l2*z
      t435 = t433*t434
      t436 = 4*t435
      t437 = l3*s2l2
      t438 = s3l1*z
      t439 = t437*t438
      t440 = 4*t439
      t441 = t364*t218
      t442 = 4*t441
      t443 = t362*t218
      t444 = 4*t443
      t445 = s4l2*z
      t446 = t433*t445
      t447 = 8*t446
      t448 = s4l1*z
      t449 = t437*t448
      t450 = 8*t449
      t451 = t204*z
      t452 = t207*z
      t453 = 2*t246
      t454 = 2*t250
      t455 = t414*t448
      t456 = 2*t455
      t457 = t245*t445
      t458 = 2*t457
      t461 = 8*t335
      t462 = 8*t336
      t463 = t436-t440+t442-t444+t447-t450+t451-t452+t453-t454+t456-t458
     #-4*t332+4*t331+t461-t462
      t466 = t1*t344
      t473 = t387*s34r
      t491 = t383*z
      t492 = 3*t415
      t493 = 3*t417
      t494 = 3*t419
      t495 = 3*t421
      t496 = s1l1*l3
      t497 = t496*s2l2
      t498 = t381*l3
      t499 = s1l2*z
      t500 = t499*s2l1
      t501 = 3*t435
      t502 = 3*t439
      t503 = -t491-t392+t391-t492+t493-t494+t495-t423+t426-t428+t430-t49
     #7+t498+t500+t501-t502
      t504 = t383*t218
      t505 = t499*t433
      t506 = -t446+t449+t504-t505+t451-t452+t453-t454+t456-t458-t381-t39
     #3+t390+t383-t335+t336
      t509 = 2*t408
      t510 = 2*t412
      t511 = 4*t419
      t512 = 4*t421
      t513 = -t491+t406+t509-t392+t391-t362+t364-t410-t510-t492+t493-t51
     #1+t512+t423-t426+t428-t430-t497+t498+t500
      t514 = 2*t441
      t515 = 2*t443
      t516 = t501-t502+t514-t515+t446-t449+t504-t505+t451-t452+t453-t454
     #+t456-t458-t381-t393+t390+t383+t335-t336
      t524 = 2*t415
      t525 = 2*t417
      t526 = 4*t423
      t527 = 4*t426
      t528 = 4*t428
      t529 = 4*t430
      t530 = 2*t435
      t531 = t509-t204+t207-t510-t524+t525-t511+t512+t526-t527+t528-t529
     #+t530
      t532 = 2*t439
      t533 = 4*t446
      t534 = 4*t449
      t535 = 4*t335
      t536 = 4*t336
      t537 = -t532+t514-t515+t533-t534+t246-t250+t455-t457-t338+t337+t53
     #5-t536
      t543 = s3l2*s1l1
      t544 = s3l1*s1l2
      t547 = t97*l3
      t548 = t2**2
      t549 = t547*t548
      t573 = -t362+t364
      t575 = t573*t547*t105
      t579 = 2*t544
      t580 = 2*t543
      t581 = t579-t397+t383-t393-t580-t363+t398-t391-t381+t390+t392+t365
      t591 = t543*z
      t593 = t496*s3l2
      t594 = 6*t593
      t595 = t544*z
      t597 = s1l2*l3
      t598 = t597*s3l1
      t599 = 6*t598
      t601 = 6*t408
      t602 = 4*t362
      t603 = 4*t364
      t604 = t496*t434
      t605 = 6*t604
      t607 = 6*t412
      t608 = 6*t423
      t609 = 6*t426
      t610 = 6*t428
      t611 = 4*t591+t594-4*t595-t599+4*t406+t601+t204-t207-t602+t603-t60
     #5-4*t410-t607-t492+t493-t494+t495+t608-t609+t610
      t612 = 6*t430
      t613 = 6*t441
      t614 = 6*t443
      t615 = 6*t446
      t616 = 6*t449
      t619 = 6*t246
      t620 = 6*t250
      t621 = 6*t455
      t622 = 6*t457
      t627 = t597*t438
      t628 = 6*t627
      t629 = 6*t335
      t630 = 6*t336
      t631 = -t612+t501-t502+t613-t614+t615-t616-4*t451+4*t452-t619+t620
     #-t621+t622-5*t332+5*t331-4*t543+4*t544+t628+t629-t630
      t648 = 2*t573*t547*t105
      t653 = -t382-t392+t384+t391-t579+t580+t338-t362+t364-t337
      t659 = t4*t547
      t663 = 3*t381
      t664 = 3*t383
      t665 = -t543+t335+t544-t336+t337-t663+t204+t664-t207-t362-t338+t36
     #4
      t673 = 2*t491
      t684 = -t673+6*t591+8*t593-6*t595-8*t598+6*t406+8*t408-6*t362+6*t3
     #64-8*t604-6*t410
      t690 = 2*t497
      t691 = 2*t498
      t692 = 2*t500
      t693 = -8*t412-8*t415+8*t417-8*t419+8*t421+t424-t427+t429-t431-t69
     #0+t691+t692
      t699 = 2*t504
      t700 = 2*t505
      t704 = 8*t435-8*t439+8*t441-8*t443+t447-t450+t699-t700-6*t451+6*t4
     #52-8*t246
      t713 = 8*t250-8*t455+8*t457-t382-12*t332+12*t331+t384-6*t543+6*t54
     #4+8*t627+t461-t462
      t719 = 2*t591
      t720 = 2*t595
      t721 = 4*t204
      t722 = 4*t207
      t723 = t719+t594-t720-t599+t407+t601-t721+t722-t397+t398-t605-t411
     #-t607-t492+t493-t494+t495+t608-t609+t610
      t724 = 2*t451
      t725 = 2*t452
      t728 = -t612+t501-t502+t613-t614+t615-t616-t724+t725-t619+t620-t62
     #1+t622-7*t332+7*t331-t580+t579+t628+t629-t630
      t745 = l3*t344
      t750 = -t381+t383
      t762 = t417-t419+t207+t691+t699-t382-t690-t415-t700-t673+t421+t384
     #+t435+t331-t204+t692-t439-t332
      t765 = -3*t491+t591+t593-t595-t598+t406+t408-t204+t207-t362+t364-t
     #604
      t770 = -t410-t412-t524+t525-2*t419+2*t421+t423-t426+t428-t430-3*t4
     #97+3*t498
      t775 = 3*t500+t530-t532+t441-t443+t446-t449+3*t504-3*t505-t451+t45
     #2-t246
      t776 = t250-t455+t457-t663-t393+t390+t664-t543+t544+t627+t335-t336
      t788 = -t491+t719+4*t593-t720-4*t598+t407+t409-t721+t722-t397+t398
     #-4*t604
      t789 = -t411-t413-t416+t418-t511+t512+t526-t527+t528-t529-t497+t49
     #8
      t792 = t500+t436-t440+t442-t444+t533-t534+t504-t505-t724+t725-4*t2
     #46
      t799 = 4*t250-4*t455+4*t457-t381-8*t332+8*t331+t383-t580+t579+4*t6
     #27+t535-t536
      t815 = 2*t593-2*t598+t509-3*t204+3*t207-2*t604-t510-t415+t417-t419
     #+t421+2*t423-2*t426+2*t428-2*t430
      t819 = t435-t439+t514-t515+2*t446-2*t449-t453+t454-t456+t458-t393+
     #t390+2*t627+t365-t363
      s2 = 2048.D0/3.D0*t323*s4l1*s3l2-(-256.D0/3.D0*t343*t344*l1*r*t9*t
     #11*t15+256.D0/3.D0*t343*t344*s34r+256.D0/3.D0*(t122*t333+t78*t339)
     #*t4*t344)*t360
      s3 = s2-(256.D0/3.D0*(2*(-t204+t207-t362-t363+t364-t338+t337+t365)
     #*t4*t106-(t207-t204+t337-t338)*t4*t106)*t344*t183*t8*t186*t188*t19
     #1+(256.D0/3.D0*t386*t387*s23r-256.D0/3.D0*(-t128*t105*(t383-t381+t
     #390+t336+t391-t392-t393-t335)+2*t128*t105*(-t397-t393+t398+t335-t3
     #81-t336+t391+t383-t392+t390))*t344*s34r-256.D0/3.D0*(t432+t463)*l3
     #*t466)*l1*r*t16+(-256.D0/3.D0*t386*t473-256.D0/3.D0*t385*l3*t466)*
     #s23r+256.D0/3.D0*(-t128*t105*(t383+t207+t331+t336-t381-t204-t332-t
     #335)+2*t128*t105*(t364-t381-t332+t207+t383-t335-t204+t331+t336-t36
     #2))*t344*t129+256.D0/3.D0*(-t24*(t503+t506)+2*t24*(t513+t516))*t34
     #4*s34r+256.D0/3.D0*(t531+t537)*t344)*s13r
      s1 = s3+256.D0/3.D0*(2*(-t362+t364+t331-t543+t335-t332+t544-t336+t
     #204-t207)*t4*t549-t333*t4*t549)*t344*t183*l1*t8*r*t186*t9/t187/t10
     #/t190/t14-(-256.D0/3.D0*(2*(-t381-t204+t383+t207-t544+t543+t332-t3
     #62+t364-t331)*t547*t105-t575)*t344*s23r+256.D0/3.D0*(2*t581*t4*t54
     #9-(t383-t381-t338+t337)*t4*t549)*t344*s34r+256.D0/3.D0*(t611+t631)
     #*t97*t387)*t183*t8*t192
      t823 = s1-((256.D0/3.D0*(2*(-t331+t383-t381-t544+t332+t207-t204+3*
     #t364-3*t362+t543)*t547*t105-t648)*t344*s34r+512.D0/3.D0*t653*t97*t
     #387)*s23r-256.D0/3.D0*(-t659*t548*(-t382-t332+t384+t331)+2*t659*t5
     #48*t665)*t344*t129-256.D0/3.D0*(t684+t693+t704+t713)*t97*t473-256.
     #D0/3.D0*(t723+t728)*l3*t466)*l1*r*t16-(-256.D0/3.D0*(-t575+2*t648)
     #*t344*t129-256.D0/3.D0*(-t337+t384-t382-t579+t338+t391-t392+t603-t
     #602+t580)*t97*t473-256.D0/3.D0*t653*t1*t745)*s23r-256.D0/3.D0*(-t7
     #50*t4*t549+4*t659*t548*t750)*t344*t129*s34r-256.D0/3.D0*(-t98*t762
     #+2*t98*(t765+t770+t775+t776))*t344*t129-256.D0/3.D0*(t788+t789+t79
     #2+t799)*t1*t745*s34r-256.D0/3.D0*(t815+t819)*t344
      taux4 = t327+t823

CCC Coefficients

      tcoff9 = taux3
      tcoff10 = taux3
      tcoff11 = taux3
      tcoff12 = taux3
      tcoff13 = taux4
      tcoff14 = taux4
      tcoff15 = -taux4
      tcoff16 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet10 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet11 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet12 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet13 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet14 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet15 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = -l3+1
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet16 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff16)

CCC Sum

      sq3num2e3 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*
     #tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tc
     #off16
      RETURN
      END

      double precision function sq3num2e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux3 = 0
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t8 = r**2
      t9 = 1+t8
      t10 = r+t
      t14 = 1+r*t
      t16 = t9/t10/t14
      t17 = l1*r*t16
      t20 = l3*t1
      t27 = -4096.D0/3.D0*t20-4096.D0/3.D0*t3*t4*s34r
      t42 = t4*s4l2
      t45 = t4*s4l1
      t53 = t20*s3l2
      t54 = t20*s4l2
      t58 = t20*s3l1
      t59 = t20*s4l1
      t67 = t1*s3l2
      t73 = t1*s3l1
      t79 = l3**2
      t80 = t79*t2
      t85 = s34r*s2l2
      t87 = s4l2*s34r
      t88 = t80*t87
      t100 = s34r*s2l1
      t102 = s4l1*s34r
      t103 = t80*t102
      t111 = t79*t4
      t112 = t2*t1
      t119 = t111*t112*s4l2*s3l1-t111*t112*s4l1*s3l2
      t120 = 4096.D0/3.D0*t119
      t121 = l1**2
      t124 = t9**2
      t125 = t10**2
      t128 = t14**2
      t130 = t124/t125/t128
      t132 = s3l2*s2l1
      t134 = s3l1*s2l2
      t136 = s4l2*s3l1
      t137 = t80*t136
      t138 = s3l2*s4l1
      t139 = t80*t138
      t151 = z*l3
      t154 = (-z+2-2*l3+2*t151)*t1
      t155 = l3*s3l1
      t156 = t155*s4l2
      t157 = t154*t156
      t159 = l3*s4l1
      t160 = t159*s3l2
      t161 = t154*t160
      t190 = s34r**2
      t195 = -l3+t151+1
      s2 = -((4096.D0/3.D0*t3*t4*s2l2*t17-4096.D0/3.D0*t20*s2l2*s23r+t27
     #*s2l2)*s1l1+(-4096.D0/3.D0*t3*t4*s2l1*t17+4096.D0/3.D0*t20*s2l1*s2
     #3r-t27*s2l1)*s1l2+(-4096.D0/3.D0*t3*t42*s2l1+4096.D0/3.D0*t3*t45*s
     #2l2)*l1*r*t16+((4096.D0/3.D0*t53+4096.D0/3.D0*t54)*s2l1+(-4096.D0/
     #3.D0*t58-4096.D0/3.D0*t59)*s2l2)*s23r+(4096.D0/3.D0*t3*t42*s34r+40
     #96.D0/3.D0*t67+4096.D0/3.D0*t54)*s2l1+(-4096.D0/3.D0*t3*t45*s34r-4
     #096.D0/3.D0*t73-4096.D0/3.D0*t59)*s2l2)*s13r-(-4096.D0/3.D0*t80*s4
     #l2*s23r*t17+(4096.D0/3.D0*t80*t85+4096.D0/3.D0*t88+4096.D0/3.D0*t5
     #4)*s23r+4096.D0/3.D0*t3*t85)*s1l1
      s1 = s2-(4096.D0/3.D0*t80*s4l1*s23r*t17+(-4096.D0/3.D0*t80*t100-40
     #96.D0/3.D0*t103-4096.D0/3.D0*t59)*s23r-4096.D0/3.D0*t3*t100)*s1l2-
     #t120*t121*t8*t130-((-4096.D0/3.D0*t80*t132+4096.D0/3.D0*t80*t134-4
     #096.D0/3.D0*t137+4096.D0/3.D0*t139)*s23r-4096.D0/3.D0*t3*t132+4096
     #.D0/3.D0*t3*t134-8192.D0/3.D0*t119*s34r-4096.D0/3.D0*t157+4096.D0/
     #3.D0*t161)*l1*r*t16
      t199 = s1-((-4096.D0/3.D0*t88+4096.D0/3.D0*t53)*s2l1+(4096.D0/3.D0
     #*t103-4096.D0/3.D0*t58)*s2l2+(4096.D0/3.D0*t137-4096.D0/3.D0*t139)
     #*s34r+4096.D0/3.D0*t20*t136-4096.D0/3.D0*t20*t138)*s23r-(-4096.D0/
     #3.D0*t3*t87+4096.D0/3.D0*t67)*s2l1-(4096.D0/3.D0*t3*t102-4096.D0/3
     #.D0*t73)*s2l2-t120*t190-(4096.D0/3.D0*t157-4096.D0/3.D0*t161)*s34r
     #-4096.D0/3.D0*t195*s4l2*s3l1
      t203 = t132-t134
      t204 = t203*t4
      t205 = 1/s134r
      t210 = t4*t205
      t217 = s13r**2
      t219 = 2*t132
      t220 = 2*t134
      t230 = s1l2*s4l1
      t231 = s1l1*s4l2
      t232 = -t230+t136-t138+t231
      t233 = t232*t79
      t234 = t2*t205
      t237 = t111*t112
      t238 = 2*t138
      t239 = 2*t136
      t240 = 3*t132
      t241 = 3*t134
      t242 = s1l1*s2l2
      t243 = s4l1*s2l2
      t244 = s1l2*s2l1
      t245 = s4l2*s2l1
      t250 = 4*t132
      t251 = 4*t134
      t252 = t136*z
      t253 = t132*z
      t254 = 4*t253
      t255 = t138*z
      t256 = t134*z
      t257 = 4*t256
      t258 = l3*s3l2
      t259 = t258*s2l1
      t260 = 4*t259
      t261 = 2*t156
      t262 = t155*s2l2
      t263 = 4*t262
      t264 = 2*t160
      t265 = l3*s2l2
      t267 = t265*s3l1*z
      t268 = 4*t267
      t269 = z*s4l2
      t270 = t155*t269
      t271 = 2*t270
      t272 = z*s4l1
      t273 = t258*t272
      t274 = 2*t273
      t275 = l3*s2l1
      t277 = t275*s3l2*z
      t278 = 4*t277
      t279 = t238-t239+t250-t251+t252-t254-t255+t257-t260+t261+t263-t264
     #-t268-t271+t274+t278
      t281 = t1*t205
      t288 = t234*s34r
      t301 = s1l1*l3*s2l2
      t302 = t242*z
      t304 = s1l2*z
      t305 = t304*s2l1
      t306 = t244*l3
      t307 = t304*t275
      t308 = t242*t151
      t313 = -t239+t238-t301-t302-t159*s2l2+t305+t306-t307+t308-t275*t26
     #9+t242+t243-t244-t245+t265*t272+l3*s4l2*s2l1
      t316 = 3*t253
      t317 = 3*t256
      t318 = 3*t259
      t319 = 3*t262
      t320 = 3*t267
      t321 = 3*t277
      t322 = t245*z-t243*z+t240-t241+t252-t316-t255+t317-t318+t261+t319-
     #t264-t320-t271+t274+t321
      t328 = 2*t253
      t329 = 2*t256
      t330 = 2*t259
      t331 = 2*t262
      t332 = 2*t267
      t333 = 2*t277
      t334 = t138-t136+t219-t220-t328+t329-t330+t156+t331-t160-t332-t270
     #+t273+t333
      t339 = t79*l3
      t340 = t2**2
      t357 = -t230+t231
      t359 = t112*t205
      t376 = 2*t357
      t387 = t4*t339*t340
      t388 = 2*t244
      t389 = 2*t242
      t395 = 2*t305
      t396 = 2*t306
      t397 = 2*t307
      t398 = 2*t308
      t399 = 2*t301
      t400 = 2*t302
      t401 = t238-t239+t250-t251+t389-t254+t257-t260+t263-t268+t278-t388
     #+t395+t396-t397+t398-t399-t400
      t420 = l3*t205
      t431 = t138-t136+t132-t134+t389-t253+t256-t259+t262-t267+t277-t388
     #+t395+t396-t397+t398-t399-t400
      t436 = t238-t239+t219-t220+t242-t328+t329-t330+t331-t332+t333-t244
     #+t305+t306-t307+t308-t301-t302
      s1 = 4096.D0/3.D0*t195*s4l1*s3l2-(-512.D0/3.D0*t3*t204*t205*t17+51
     #2.D0/3.D0*t3*t203*t210*s34r+512.D0/3.D0*t1*t203*t210)*t217-(512.D0
     #/3.D0*(t138-t136+t219-t220)*t4*t79*t112*t205*t121*t8*t130+(512.D0/
     #3.D0*t233*t234*s23r-512.D0/3.D0*t237*(t238-t239+t240-t241+t242+t24
     #3-t244-t245)*t205*s34r-512.D0/3.D0*t279*l3*t281)*l1*r*t16+(-512.D0
     #/3.D0*t233*t288-512.D0/3.D0*t232*l3*t281)*s23r+512.D0/3.D0*t237*(t
     #138-t136+t132-t134+t242+t243-t244-t245)*t205*t190+512.D0/3.D0*t20*
     #(t313+t322)*t205*s34r+512.D0/3.D0*t334*t205)*s13r+512.D0/3.D0*t204
     #*t339*t340*t205*t121*l1*t8*r*t124*t9/t125/t10/t128/t14-(-512.D0/3.
     #D0*t357*t339*t359*s23r+512.D0/3.D0*(t242-t244-t220+t219)*t4*t339*t
     #340*t205*s34r+512.D0/3.D0*(t138-t136+t240-t241-t316+t317-t318+t319
     #-t320+t321)*t79*t234)*t121*t8*t130
      t444 = s1-((512.D0/3.D0*t376*t339*t359*s34r+1024.D0/3.D0*t357*t79*
     #t234)*s23r-512.D0/3.D0*t387*(t132-t134-t388+t389)*t205*t190-512.D0
     #/3.D0*t401*t79*t288-512.D0/3.D0*(t238-t239+t240-t241-t316+t317-t31
     #8+t319-t320+t321)*l3*t281)*l1*r*t16-(-512.D0/3.D0*t339*t112*t357*t
     #205*t190-512.D0/3.D0*t376*t79*t288-512.D0/3.D0*t357*t1*t420)*s23r-
     #512.D0/3.D0*t387*(t242-t244)*t205*t190*s34r-512.D0/3.D0*t80*t431*t
     #205*t190-512.D0/3.D0*t436*t1*t420*s34r-512.D0/3.D0*(t138-t136+t132
     #-t134-t253+t256-t259+t262-t267+t277)*t205
      taux4 = t199+t444

CCC Coefficients

      tcoff9 = taux3
      tcoff10 = taux3
      tcoff11 = taux3
      tcoff12 = taux3
      tcoff13 = taux4
      tcoff14 = taux4
      tcoff15 = -taux4
      tcoff16 = -taux4

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+z
     #-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet10 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet11 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet12 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet13 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet14 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet15 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet16 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff16)

CCC Sum

      sq3num2e4 = tFjet9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*
     #tcoff12+tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tc
     #off16
      RETURN
      END

      double precision function sq3num3e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = 1/s234r
      t5 = t4*s3l2
      t6 = t3*t5
      t7 = t4*s4l2
      t8 = t3*t7
      t10 = -64.D0/3.D0*t6-64.D0/3.D0*t8
      t12 = t1**2
      t13 = t2**2
      t14 = t12*t13
      t15 = t4*s34r
      t17 = t14*t15*s2l2
      t22 = 64.D0/3.D0*t14*t5+64.D0/3.D0*t14*t7
      t26 = t4*s3l1
      t28 = t4*s4l1
      t31 = -64.D0/3.D0*t3*t26-64.D0/3.D0*t3*t28
      t39 = 64.D0/3.D0*t14*t26+64.D0/3.D0*t14*t28
      t43 = l3*t12
      t44 = t43*t2
      t45 = z*t4
      t55 = z*l3
      t73 = l3*t1
      t92 = t2*t4*s34r
      t94 = t1*t4
      t96 = -64.D0/3.D0*t43*t92+64.D0/3.D0*t94
      t98 = t7*s34r
      t124 = t28*s34r
      t137 = z*t12
      t138 = l3**2
      t140 = s23r**2
      t153 = -128.D0/3.D0*t5-64.D0/3.D0*t7
      t155 = t28*s3l2
      t156 = 64.D0/3.D0*t155
      t176 = t3*t4
      s2 = ((t10*s1vr+64.D0/3.D0*t17+t22*s34r)*s1l1+(t31*s1vr+64.D0/3.D0
     #*t14*t15*s2l1+t39*s34r)*s1l2+(64.D0/3.D0*t44*t45*s23r-t10*s2l1-t31
     #*s2l2)*s1vr-64.D0/3.D0*t55*t12*t1*t13*t4*s34r*s23r+(-128.D0/3.D0*t
     #17-t22*s34r)*s2l1-t39*s34r*s2l2)*s13r+(128.D0/3.D0*t73*t4*s23r*s1v
     #r*s1l2+((64.D0/3.D0*t73*t4*s2l2+64.D0/3.D0*t73*t5+128.D0/3.D0*t73*
     #t7)*s23r-128.D0/3.D0*t7)*s1vr+(t96*s2l2-64.D0/3.D0*t44*t98+64.D0/3
     #.D0*t45*s3l2+64.D0/3.D0*t45*s4l2)*s23r+64.D0/3.D0*t3*t98)*s1l1
      s1 = s2+(((64.D0/3.D0*t73*t4*s2l1+64.D0/3.D0*t73*t26+128.D0/3.D0*t
     #73*t28)*s23r-128.D0/3.D0*t28)*s1vr+(t96*s2l1-64.D0/3.D0*t44*t124+6
     #4.D0/3.D0*t45*s3l1+64.D0/3.D0*t45*s4l1)*s23r+64.D0/3.D0*t3*t124)*s
     #1l2+(-64.D0/3.D0*t137*t138*t4*t140+(-64.D0/3.D0*t73*t5*s2l1-64.D0/
     #3.D0*t73*t26*s2l2)*s23r+t153*s3l1-t156)*s1vr
      taux5 = s1-64.D0/3.D0*t137*l3*t4*t140+(-64.D0/3.D0*t94*s3l2*s2l1-6
     #4.D0/3.D0*t94*s3l1*s2l2+64.D0/3.D0*t55*t12*t92-t153*s3l1+t156)*s23
     #r+64.D0/3.D0*t176*s3l2*s34r*s2l1+64.D0/3.D0*t176*s3l1*s34r*s2l2+((
     #128.D0/3.D0*t6+64.D0/3.D0*t8)*s3l1+64.D0/3.D0*t3*t155)*s34r
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = 1/s234r
      t5 = t4*s3l2
      t7 = -t1
      t8 = t7*t2
      t9 = t4*s4l2
      t10 = t8*t9
      t12 = -256.D0/3.D0*t3*t5+256.D0/3.D0*t10
      t14 = t1**2
      t15 = t2**2
      t16 = t14*t15
      t17 = t4*s34r
      t24 = 256.D0/3.D0*t16*t5+256.D0/3.D0*t16*t9
      t28 = t4*s3l1
      t30 = t4*s4l1
      t31 = t8*t30
      t33 = 256.D0/3.D0*t3*t28-256.D0/3.D0*t31
      t41 = -256.D0/3.D0*t16*t28-256.D0/3.D0*t16*t30
      t67 = t7*l3
      t71 = t1*l3
      t74 = t67*t9
      t84 = s23r**2
      t86 = t7*t4
      t89 = l3*t14
      t90 = t2*t4
      t94 = -t2*l3
      t95 = t14*t4
      t101 = z*t4
      t118 = t67*t30
      t149 = t9*s3l1
      t150 = t30*s3l2
      t154 = t8*t4
      s2 = ((t12*s1vr+256.D0/3.D0*t16*t17*s2l2+t24*s34r)*s1l1+(t33*s1vr-
     #256.D0/3.D0*t16*t17*s2l1+t41*s34r)*s1l2+((-256.D0/3.D0*t8*t5-256.D
     #0/3.D0*t10)*s2l1+(256.D0/3.D0*t8*t28+256.D0/3.D0*t31)*s2l2)*s1vr+(
     #t12*s2l1+t33*s2l2)*s23r-t24*s34r*s2l1-t41*s34r*s2l2)*s13r
      s3 = (((-256.D0/3.D0*t67*t4*s2l2+512.D0/3.D0*t71*t5-512.D0/3.D0*t7
     #4)*s23r-512.D0/3.D0*t9)*s1vr+(256.D0/3.D0*t67*t5+256.D0/3.D0*t74)*
     #t84+(-256.D0/3.D0*t86*s2l2+(-256.D0/3.D0*t89*t90*s3l2+256.D0/3.D0*
     #t94*t95*s4l2)*s34r+256.D0/3.D0*t101*s3l2+256.D0/3.D0*t101*s4l2)*s2
     #3r-256.D0/3.D0*t8*t9*s34r)*s1l1+(((256.D0/3.D0*t67*t4*s2l1-512.D0/
     #3.D0*t71*t28+512.D0/3.D0*t118)*s23r+512.D0/3.D0*t30)*s1vr+(-256.D0
     #/3.D0*t67*t28-256.D0/3.D0*t118)*t84+(256.D0/3.D0*t86*s2l1+(256.D0/
     #3.D0*t89*t90*s3l1-256.D0/3.D0*t94*t95*s4l1)*s34r-256.D0/3.D0*t101*
     #s3l1-256.D0/3.D0*t101*s4l1)*s23r+256.D0/3.D0*t8*t30*s34r)*s1l2
      s1 = s2+s3
      taux6 = s1+(-256.D0/3.D0*t149+256.D0/3.D0*t150)*s1vr+256.D0/3.D0*t
     #154*s3l2*s34r*s2l1-256.D0/3.D0*t154*s3l1*s34r*s2l2+(256.D0/3.D0*t3
     #*t149-256.D0/3.D0*t3*t150)*s34r

CCC Coefficients

      tcoff17 = taux5
      tcoff18 = taux5
      tcoff19 = taux5
      tcoff20 = taux5
      tcoff21 = taux6
      tcoff22 = taux6
      tcoff23 = -taux6
      tcoff24 = -taux6

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet17 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet19 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet20 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet23 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet24 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      sq3num3e0 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet2
     #0*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*
     #tcoff24
      RETURN
      END

      double precision function sq3num3e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = 1/s234r
      t5 = t4*s4l2
      t8 = t4*s2l2
      t11 = t1**2
      t12 = t2**2
      t13 = t11*t12
      t14 = t5*s34r
      t19 = t4*s4l1
      t22 = t4*s2l1
      t25 = t19*s34r
      t30 = z*l3
      t31 = t30*t11
      t32 = t2*t4
      t33 = s1vr**2
      t37 = l3*t11
      t38 = t37*t2
      t39 = z*t4
      t45 = 64.D0/3.D0*t3*t5*s2l1
      t48 = 64.D0/3.D0*t3*t19*s2l2
      t49 = t11*t1
      t50 = t30*t49
      t51 = t12*t4
      t55 = t4*s3l2
      t59 = 64.D0/3.D0*t3*t5
      t62 = t19*s3l2
      t64 = 64.D0/3.D0*t3*t62
      t72 = s23r**2
      t76 = t5*s3l1
      t84 = t13*t4
      t86 = s3l2*s34r*s2l1
      t90 = s3l1*s34r*s2l2
      t93 = t11**2
      t98 = s34r**2
      t110 = -t1
      t118 = t1*l3
      t132 = t118*t55
      t134 = t118*t5
      t148 = t1*t4
      t164 = t49*t12*l3
      t165 = t4*t98
      t178 = t4*s3l1
      t179 = t118*t178
      t181 = t118*t19
      t219 = z*t11
      t220 = l3**2
      t224 = t55*s3l1
      t253 = t11*t2
      t254 = t253*t4
      t280 = t110*t4
      s3 = ((64.D0/3.D0*t3*t5*s1vr+64.D0/3.D0*t3*t8*s23r-64.D0/3.D0*t13*
     #t14)*s1l1+(64.D0/3.D0*t3*t19*s1vr+64.D0/3.D0*t3*t22*s23r-64.D0/3.D
     #0*t13*t25)*s1l2+128.D0/3.D0*t31*t32*t33+(128.D0/3.D0*t38*t39*s23r-
     #t45-t48-256.D0/3.D0*t50*t51*s34r+(-128.D0/3.D0*t3*t55-t59)*s3l1-t6
     #4-64.D0/3.D0*(3-z)*z*t3*t4)*s1vr-128.D0/3.D0*t31*t32*t72+(t45+t48+
     #64.D0/3.D0*t3*t76+t64+128.D0/3.D0*t3*t5*s4l1)*s23r+64.D0/3.D0*t84*
     #t86+64.D0/3.D0*t84*t90+128.D0/3.D0*z*t93*t12*t2*l3*t4*t98+((128.D0
     #/3.D0*t13*t55+64.D0/3.D0*t13*t5)*s3l1+64.D0/3.D0*t13*t62+64.D0/3.D
     #0*t110*z*t84)*s34r)*s13r
      s4 = ((-128.D0/3.D0*t118*t4*s23r*s1vr+128.D0/3.D0*t118*t4*t72+128.
     #D0/3.D0*t38*t4*s34r*s23r)*s1l2+((-64.D0/3.D0*t118*t8-64.D0/3.D0*t1
     #32-128.D0/3.D0*t134)*s23r+64.D0/3.D0*t38*t55*s34r+128.D0/3.D0*t5)*
     #s1vr+(128.D0/3.D0*t132+64.D0/3.D0*t134)*t72+(-64.D0/3.D0*t148*s2l2
     #+(64.D0/3.D0*t37*t32*s3l2+128.D0/3.D0*t37*t32*s4l2)*s34r-64.D0/3.D
     #0*t39*s4l2)*s23r-64.D0/3.D0*t164*t165*s2l2-64.D0/3.D0*t164*t55*t98
     #-64.D0/3.D0*t3*t14)*s1l1
      s2 = s3+s4
      s1 = s2+(((-64.D0/3.D0*t118*t22-64.D0/3.D0*t179-128.D0/3.D0*t181)*
     #s23r+64.D0/3.D0*t38*t178*s34r+128.D0/3.D0*t19)*s1vr+(128.D0/3.D0*t
     #179+64.D0/3.D0*t181)*t72+(-64.D0/3.D0*t148*s2l1+(64.D0/3.D0*t37*t3
     #2*s3l1+128.D0/3.D0*t37*t32*s4l1)*s34r-64.D0/3.D0*t39*s4l1)*s23r-64
     #.D0/3.D0*t164*t165*s2l1-64.D0/3.D0*t164*t178*t98-64.D0/3.D0*t3*t25
     #)*s1l2+(-128.D0/3.D0*t219*t220*t4*s23r+128.D0/3.D0*t118*t224+128.D
     #0/3.D0*t30*t148)*t33+((-64.D0/3.D0*t118*t55*s2l1-64.D0/3.D0*t118*t
     #178*s2l2-64.D0/3.D0*t118*t76-64.D0/3.D0*t118*t62+64.D0/3.D0*(1-3*z
     #)*z*t118*t4)*s23r-64.D0/3.D0*t38*t55*s34r*s2l1-64.D0/3.D0*t38*t178
     #*s34r*s2l2+(-128.D0/3.D0*t38*t224-128.D0/3.D0*t30*t254)*s34r+(128.
     #D0/3.D0*t148*s3l2+64.D0/3.D0*t5)*s3l1+64.D0/3.D0*t62+128.D0/3.D0*t
     #1*z*t4)*s1vr
      taux5 = s1+(-64.D0/3.D0*t148*s3l2*s2l1-64.D0/3.D0*t148*s3l1*s2l2+1
     #28.D0/3.D0*t31*t32*s34r+(128.D0/3.D0*t55+64.D0/3.D0*t280*s4l2)*s3l
     #1+64.D0/3.D0*t280*s4l1*s3l2+64.D0/3.D0*(-1-z)*t1*t39)*s23r-64.D0/3
     #.D0*t254*t86-64.D0/3.D0*t254*t90+128.D0/3.D0*t50*t51*t98+((-128.D0
     #/3.D0*t253*t55+t59)*s3l1+t64-64.D0/3.D0*t219*t32)*s34r
      t1 = 1-z
      t3 = l3-1
      t4 = 3*t1*t3
      t5 = 1/s234r
      t6 = t5*s3l2
      t8 = -t1
      t9 = 4*t8
      t10 = t9*t3
      t11 = t5*s4l2
      t12 = t10*t11
      t16 = t8*t3
      t17 = t5*s2l2
      t21 = t8**2
      t22 = t3**2
      t23 = t21*t22
      t24 = t5*s34r
      t28 = t23*t6
      t30 = t23*t11
      t36 = t5*s3l1
      t38 = t5*s4l1
      t39 = t10*t38
      t43 = t5*s2l1
      t50 = t23*t36
      t52 = t23*t38
      t66 = -t3
      t67 = t66*t8
      t68 = t11*s3l1
      t71 = t38*s3l2
      t100 = t9*l3
      t103 = t8*l3
      t104 = t103*t6
      t110 = l3*t21
      t111 = t110*t3
      t119 = 2*t8*l3
      t123 = s23r**2
      t131 = 512.D0/3.D0*t110*t3*t5*s34r-256.D0/3.D0*t9*t5
      t134 = 2*t66*l3
      t135 = t21*t5
      t139 = 3*t3*l3
      t145 = z*t5
      t154 = t21*t8*t22*l3
      t155 = s34r**2
      t156 = t5*t155
      t170 = t103*t36
      t218 = t1*l3
      t226 = t6*s34r*s2l1
      t230 = t36*s34r*s2l2
      t237 = t8*t5
      t254 = (2-z)*t8*t3
      s3 = (((-256.D0/3.D0*t4*t6+256.D0/3.D0*t12)*s1vr+256.D0/3.D0*t16*t
     #17*s23r-256*t23*t24*s2l2+(-256*t28-1024.D0/3.D0*t30)*s34r)*s1l1+((
     #256.D0/3.D0*t4*t36-256.D0/3.D0*t39)*s1vr-256.D0/3.D0*t16*t43*s23r+
     #256*t23*t24*s2l1+(256*t50+1024.D0/3.D0*t52)*s34r)*s1l2+((-256.D0/3
     #.D0*t16*t6-256.D0/3.D0*t12)*s2l1+(256.D0/3.D0*t16*t36+256.D0/3.D0*
     #t39)*s2l2-256.D0/3.D0*t67*t68+256.D0/3.D0*t67*t71)*s1vr+(256.D0/3.
     #D0*t16*t11*s2l1-256.D0/3.D0*t16*t38*s2l2+256.D0/3.D0*t16*t68-256.D
     #0/3.D0*t16*t71)*s23r+(512.D0/3.D0*t28+512.D0/3.D0*t30)*s34r*s2l1+(
     #-512.D0/3.D0*t50-512.D0/3.D0*t52)*s34r*s2l2+(-256.D0/3.D0*t23*t68+
     #256.D0/3.D0*t23*t71)*s34r)*s13r
      s4 = (((-256.D0/3.D0*t100*t17-1024.D0/3.D0*t104-512.D0/3.D0*t100*t
     #11)*s23r-256.D0/3.D0*t111*t6*s34r+2048.D0/3.D0*t11)*s1vr+(256.D0/3
     #.D0*t104+256.D0/3.D0*t119*t11)*t123+(t131*s2l2+(-256.D0/3.D0*t134*
     #t135*s3l2+256.D0/3.D0*t139*t135*s4l2)*s34r-256*t145*s3l2-1024.D0/3
     #.D0*t145*s4l2)*s23r+256.D0/3.D0*t154*t156*s2l2+256.D0/3.D0*t154*t6
     #*t155-256.D0/3.D0*t10*t11*s34r)*s1l1
      s2 = s3+s4
      s1 = s2+(((256.D0/3.D0*t100*t43+1024.D0/3.D0*t170+512.D0/3.D0*t100
     #*t38)*s23r+256.D0/3.D0*t111*t36*s34r-2048.D0/3.D0*t38)*s1vr+(-256.
     #D0/3.D0*t170-256.D0/3.D0*t119*t38)*t123+(-t131*s2l1+(256.D0/3.D0*t
     #134*t135*s3l1-256.D0/3.D0*t139*t135*s4l1)*s34r+256*t145*s3l1+1024.
     #D0/3.D0*t145*s4l1)*s23r-256.D0/3.D0*t154*t156*s2l1-256.D0/3.D0*t15
     #4*t36*t155+256.D0/3.D0*t10*t38*s34r)*s1l2+((256.D0/3.D0*t103*t6*s2
     #l1-256.D0/3.D0*t103*t36*s2l2+512.D0/3.D0*t218*t68-512.D0/3.D0*t218
     #*t71)*s23r-256.D0/3.D0*t111*t226+256.D0/3.D0*t111*t230+1024.D0/3.D
     #0*t68-1024.D0/3.D0*t71)*s1vr
      taux6 = s1+(256.D0/3.D0*t237*s3l2*s2l1-256.D0/3.D0*t237*s3l1*s2l2-
     #512.D0/3.D0*t145*s4l2*s3l1+512.D0/3.D0*t145*s4l1*s3l2)*s23r+256.D0
     #/3.D0*t254*t226-256.D0/3.D0*t254*t230+(256.D0/3.D0*t4*t68-256.D0/3
     #.D0*t4*t71)*s34r

CCC Coefficients

      tcoff17 = taux5
      tcoff18 = taux5
      tcoff19 = taux5
      tcoff20 = taux5
      tcoff21 = taux6
      tcoff22 = taux6
      tcoff23 = -taux6
      tcoff24 = -taux6

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet17 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet19 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet20 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet23 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet24 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      sq3num3e1 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet2
     #0*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*
     #tcoff24
      RETURN
      END

      double precision function sq3num3e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = z*t2
      t4 = l3-1
      t5 = t4**2
      t6 = 1/s234r
      t7 = t5*t6
      t10 = t2*t1
      t12 = t5*t4
      t18 = s13r**2
      t20 = z*l3
      t21 = t20*t2
      t22 = t4*t6
      t23 = s1vr**2
      t27 = l3*t2
      t28 = t27*t4
      t29 = z*t6
      t33 = t1*t4
      t34 = t6*s3l2
      t37 = t6*s4l2
      t39 = 64.D0/3.D0*t33*t37
      t40 = -128.D0/3.D0*t33*t34-t39
      t42 = t6*s4l1
      t43 = t42*s3l2
      t45 = 64.D0/3.D0*t33*t43
      t48 = t33*t6
      t53 = s23r**2
      t57 = t34*s2l1
      t60 = t6*s3l1
      t61 = t60*s2l2
      t64 = t20*t10
      t70 = 1+z
      t76 = t2*t5
      t77 = t76*t6
      t79 = s3l2*s34r*s2l1
      t83 = s3l1*s34r*s2l2
      t86 = t2**2
      t89 = l3*t6
      t90 = s34r**2
      t110 = s34r*s1vr
      t114 = l3*t1
      t118 = t6*s34r
      t132 = t10*t5*l3
      t133 = t6*t90
      t168 = l3**2
      t169 = t168*t6
      t172 = t34*s3l1
      t173 = t114*t172
      t208 = t2*t4
      t209 = t208*t6
      t214 = t1*t6
      s2 = (64.D0/3.D0*t3*t7*s1vr-64.D0/3.D0*z*t10*t12*t6*s34r)*t18+(-12
     #8.D0/3.D0*t21*t22*t23+(-128.D0/3.D0*t28*t29*s23r+t40*s3l1-t45-64.D
     #0/3.D0*(-3+z)*z*t48)*s1vr+128.D0/3.D0*t21*t22*t53+(64.D0/3.D0*t33*
     #t57+64.D0/3.D0*t33*t61+320.D0/3.D0*t64*t7*s34r-t40*s3l1+t45-64.D0/
     #3.D0*t70*z*t48)*s23r+64.D0/3.D0*t77*t79+64.D0/3.D0*t77*t83+128.D0/
     #3.D0*z*t86*t12*t89*t90+((128.D0/3.D0*t76*t34+64.D0/3.D0*t76*t37)*s
     #3l1+64.D0/3.D0*t76*t43+64.D0/3.D0*(-2+z)*z*t77)*s34r)*s13r+(64.D0/
     #3.D0*t28*t34*t110+64.D0/3.D0*t114*t37*t53+(-64.D0/3.D0*t28*t118*s2
     #l2+(-64.D0/3.D0*t27*t22*s3l2+64.D0/3.D0*t27*t22*s4l2)*s34r)*s23r-6
     #4.D0/3.D0*t132*t133*s2l2-64.D0/3.D0*t132*t34*t90)*s1l1
      s1 = s2+(64.D0/3.D0*t28*t60*t110+64.D0/3.D0*t114*t42*t53+(-64.D0/3
     #.D0*t28*t118*s2l1+(-64.D0/3.D0*t27*t22*s3l1+64.D0/3.D0*t27*t22*s4l
     #1)*s34r)*s23r-64.D0/3.D0*t132*t133*s2l1-64.D0/3.D0*t132*t60*t90)*s
     #1l2+(128.D0/3.D0*t3*t169*s23r+128.D0/3.D0*t173-128.D0/3.D0*t114*t2
     #9)*t23+(-64*t3*t169*t53+(-64.D0/3.D0*t114*t57-64.D0/3.D0*t114*t61-
     #256.D0/3.D0*z*t4*t10*t169*s34r-128.D0/3.D0*t173+64.D0/3.D0*(-1+3*z
     #)*z*t114*t6)*s23r-64.D0/3.D0*t28*t34*s34r*s2l1-64.D0/3.D0*t28*t60*
     #s34r*s2l2+(-128.D0/3.D0*t28*t172+128.D0/3.D0*t20*t209)*s34r+128.D0
     #/3.D0*t214*s3l2*s3l1-128.D0/3.D0*t1*z*t6)*s1vr
      taux5 = s1-64*t3*t89*t53+(-64.D0/3.D0*t214*s3l2*s2l1-64.D0/3.D0*t2
     #14*s3l1*s2l2-64.D0/3.D0*(-5+4*z)*z*l3*t208*t118+(-128.D0/3.D0*t214
     #*s3l2+64.D0/3.D0*t37)*s3l1+64.D0/3.D0*t43+64.D0/3.D0*t70*t1*t29)*s
     #23r-64.D0/3.D0*t209*t79-64.D0/3.D0*t209*t83+128.D0/3.D0*t64*t7*t90
     #+((-128.D0/3.D0*t208*t34+t39)*s3l1+t45+64.D0/3.D0*t3*t22)*s34r
      t1 = -1+z
      t2 = l3-1
      t3 = t2*t1
      t4 = 1/s234r
      t5 = t4*s3l2
      t7 = 512.D0/3.D0*t3*t5
      t8 = -t1
      t9 = 5*t8
      t10 = t9*t2
      t11 = t4*s4l2
      t13 = 256.D0/3.D0*t10*t11
      t16 = -t2
      t18 = 3*t16*t1
      t19 = t4*s2l2
      t23 = t1**2
      t24 = t2**2
      t25 = t23*t24
      t26 = t4*s34r
      t31 = 512.D0/3.D0*t25*t5
      t32 = t25*t11
      t38 = t4*s3l1
      t40 = 512.D0/3.D0*t3*t38
      t41 = t4*s4l1
      t43 = 256.D0/3.D0*t10*t41
      t46 = t4*s2l1
      t54 = 512.D0/3.D0*t25*t38
      t55 = t25*t41
      t79 = t11*s3l1
      t82 = t41*s3l2
      t97 = t9*l3
      t107 = l3*t1
      t109 = 2*t8
      t110 = t109*l3
      t114 = s23r**2
      t123 = -256*l3*t23*t2*t4*s34r-256.D0/3.D0*t9*t4
      t125 = t16*l3
      t126 = t23*t4
      t130 = 5*t16*l3
      t136 = z*t4
      t182 = 4*t8
      t183 = t182*l3
      t200 = t182*t4
      t207 = (-1+4*z)*t4
      t217 = -t109*t2*t4
      s2 = (((-t7+t13)*s1vr+256.D0/3.D0*t18*t19*s23r+512.D0/3.D0*t25*t26
     #*s2l2+(t31+1280.D0/3.D0*t32)*s34r)*s1l1+((t40-t43)*s1vr-256.D0/3.D
     #0*t18*t46*s23r-512.D0/3.D0*t25*t26*s2l1+(-t54-1280.D0/3.D0*t55)*s3
     #4r)*s1l2+((-t7-t13)*s2l1+(t40+t43)*s2l2)*s1vr+((-256.D0/3.D0*t10*t
     #5+512.D0/3.D0*t3*t11)*s2l1+(256.D0/3.D0*t10*t38-512.D0/3.D0*t3*t41
     #)*s2l2+256.D0/3.D0*t18*t79-256.D0/3.D0*t18*t82)*s23r+(t31-256.D0/3
     #.D0*t32)*s34r*s2l1+(-t54+256.D0/3.D0*t55)*s34r*s2l2)*s13r+(((-256.
     #D0/3.D0*t97*t19-512.D0/3.D0*t97*t11)*s23r-2560.D0/3.D0*t11)*s1vr+(
     #256.D0/3.D0*t107*t5+256.D0/3.D0*t110*t11)*t114+(t123*s2l2+(-256.D0
     #/3.D0*t125*t126*s3l2+256.D0/3.D0*t130*t126*s4l2)*s34r+512.D0/3.D0*
     #t136*s3l2+1280.D0/3.D0*t136*s4l2)*s23r-256.D0/3.D0*t10*t11*s34r)*s
     #1l1
      s1 = s2+(((256.D0/3.D0*t97*t46+512.D0/3.D0*t97*t41)*s23r+2560.D0/3
     #.D0*t41)*s1vr+(-256.D0/3.D0*t107*t38-256.D0/3.D0*t110*t41)*t114+(-
     #t123*s2l1+(256.D0/3.D0*t125*t126*s3l1-256.D0/3.D0*t130*t126*s4l1)*
     #s34r-512.D0/3.D0*t136*s3l1-1280.D0/3.D0*t136*s4l1)*s23r+256.D0/3.D
     #0*t10*t41*s34r)*s1l2+((256.D0/3.D0*t183*t5*s2l1-256.D0/3.D0*t183*t
     #38*s2l2+1024.D0/3.D0*t107*t79-1024.D0/3.D0*t107*t82)*s23r-1280.D0/
     #3.D0*t79+1280.D0/3.D0*t82)*s1vr
      taux6 = s1+(256.D0/3.D0*t200*s3l2*s2l1-256.D0/3.D0*t200*s3l1*s2l2+
     #256.D0/3.D0*t207*s4l2*s3l1-256.D0/3.D0*t207*s4l1*s3l2)*s23r+256.D0
     #/3.D0*t217*s3l2*s34r*s2l1-256.D0/3.D0*t217*s3l1*s34r*s2l2+(256.D0/
     #3.D0*t3*t79-256.D0/3.D0*t3*t82)*s34r

CCC Coefficients

      tcoff17 = taux5
      tcoff18 = taux5
      tcoff19 = taux5
      tcoff20 = taux5
      tcoff21 = taux6
      tcoff22 = taux6
      tcoff23 = -taux6
      tcoff24 = -taux6

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet17 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet19 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet20 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet23 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet24 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      sq3num3e2 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet2
     #0*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*
     #tcoff24
      RETURN
      END

      double precision function sq3num3e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = 1/s234r
      t11 = t2*t1
      t12 = z*l3*t11
      t13 = t4**2
      t14 = t13*t6
      t20 = s34r*s23r
      t24 = t2**2
      t28 = l3*t6
      t29 = s34r**2
      t30 = t28*t29
      t35 = z*t2
      t36 = l3**2
      t37 = t36*t6
      t38 = s23r**2
      taux5 = ((64.D0/3.D0*t3*t4*z*t6*s23r-64.D0/3.D0*t12*t14*s34r)*s1vr
     #+64.D0/3.D0*t12*t14*t20+64.D0/3.D0*z*t24*t13*t4*t30)*s13r+(-64.D0/
     #3.D0*t35*t37*t38-64.D0/3.D0*z*t4*t11*t37*t20)*s1vr-64.D0/3.D0*t35*
     #t28*t38-64.D0/3.D0*(-2+z)*z*t3*t4*t6*t20+64.D0/3.D0*z*t11*t13*t30
      t1 = -1+z
      t2 = l3-1
      t3 = t2*t1
      t4 = 1/s234r
      t5 = t4*s4l2
      t8 = t4*s2l2
      t11 = t1**2
      t12 = t2**2
      t13 = t11*t12
      t14 = t5*s34r
      t19 = t4*s4l1
      t22 = t4*s2l1
      t25 = t19*s34r
      t37 = 3*t2*t1
      t38 = t5*s3l1
      t41 = t19*s3l2
      t46 = t4*s3l2
      t47 = t46*s2l1
      t49 = t4*s3l1
      t50 = t49*s2l2
      t52 = t3*t38
      t53 = t3*t41
      t78 = l3*t1
      t85 = l3*t11
      t86 = t85*t2
      t95 = 3*t1
      t96 = t95*l3
      t100 = s23r**2
      t102 = t2*t4
      t108 = -256.D0/3.D0*t85*t102*s34r-512.D0/3.D0*t1*t4
      t113 = 5*l3*t2
      t114 = t11*t4
      t120 = z*t4
      t127 = t11*t1*t12*l3
      t128 = s34r**2
      t129 = t4*t128
      t187 = t46*s34r*s2l1
      t191 = t49*s34r*s2l2
      t198 = t95*t4
      t207 = t95*t1*t2
      s2 = ((512.D0/3.D0*t3*t5*s1vr+512.D0/3.D0*t3*t8*s23r-512.D0/3.D0*t
     #13*t14)*s1l1+(-512.D0/3.D0*t3*t19*s1vr-512.D0/3.D0*t3*t22*s23r+512
     #.D0/3.D0*t13*t25)*s1l2+(-512.D0/3.D0*t3*t5*s2l1+512.D0/3.D0*t3*t19
     #*s2l2-256.D0/3.D0*t37*t38+256.D0/3.D0*t37*t41)*s1vr+(-256.D0/3.D0*
     #t3*t47+256.D0/3.D0*t3*t50+256.D0/3.D0*t52-256.D0/3.D0*t53)*s23r+(-
     #256*t13*t46+512.D0/3.D0*t13*t5)*s34r*s2l1+(256*t13*t49-512.D0/3.D0
     #*t13*t19)*s34r*s2l2+(256*t13*t38-256*t13*t41)*s34r)*s13r+(((-512.D
     #0/3.D0*t78*t8-1024.D0/3.D0*t78*t5)*s23r+256*t86*t46*s34r+1024.D0/3
     #.D0*t5)*s1vr+(512.D0/3.D0*t78*t46+256.D0/3.D0*t96*t5)*t100+(t108*s
     #2l2+(-256.D0/3.D0*t85*t102*s3l2+256.D0/3.D0*t113*t114*s4l2)*s34r-5
     #12.D0/3.D0*t120*s4l2)*s23r-256*t127*t129*s2l2-256*t127*t46*t128-51
     #2.D0/3.D0*t3*t14)*s1l1
      s1 = s2+(((512.D0/3.D0*t78*t22+1024.D0/3.D0*t78*t19)*s23r-256*t86*
     #t49*s34r-1024.D0/3.D0*t19)*s1vr+(-512.D0/3.D0*t78*t49-256.D0/3.D0*
     #t96*t19)*t100+(-t108*s2l1+(256.D0/3.D0*t85*t102*s3l1-256.D0/3.D0*t
     #113*t114*s4l1)*s34r+512.D0/3.D0*t120*s4l1)*s23r+256*t127*t129*s2l1
     #+256*t127*t49*t128+512.D0/3.D0*t3*t25)*s1l2+((256.D0/3.D0*t96*t47-
     #256.D0/3.D0*t96*t50)*s23r+256*t86*t187-256*t86*t191+512.D0/3.D0*t3
     #8-512.D0/3.D0*t41)*s1vr
      taux6 = s1+(256.D0/3.D0*t198*s3l2*s2l1-256.D0/3.D0*t198*s3l1*s2l2+
     #256.D0/3.D0*t38-256.D0/3.D0*t41)*s23r+256.D0/3.D0*t207*t187-256.D0
     #/3.D0*t207*t191+(256.D0/3.D0*t52-256.D0/3.D0*t53)*s34r

CCC Coefficients

      tcoff17 = taux5
      tcoff18 = taux5
      tcoff19 = taux5
      tcoff20 = taux5
      tcoff21 = taux6
      tcoff22 = taux6
      tcoff23 = -taux6
      tcoff24 = -taux6

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet17 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet19 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet20 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet23 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet24 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      sq3num3e3 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet2
     #0*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*
     #tcoff24
      RETURN
      END

      double precision function sq3num3e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux5 = 0
      t1 = l3-1
      t2 = -1+z
      t3 = t1*t2
      t4 = 1/s234r
      t5 = t4*s4l2
      t6 = t5*s3l1
      t7 = t3*t6
      t8 = t4*s4l1
      t9 = t8*s3l2
      t10 = t3*t9
      t12 = -512.D0/3.D0*t7+512.D0/3.D0*t10
      t14 = t4*s3l2
      t15 = t14*s2l1
      t17 = t4*s3l1
      t18 = t17*s2l2
      t23 = t2**2
      t24 = t1**2
      t25 = t23*t24
      t26 = t25*t4
      t28 = s3l2*s34r*s2l1
      t32 = s3l1*s34r*s2l2
      t43 = l3*t23*t1
      t44 = s34r*s1vr
      t48 = l3*t2
      t49 = s23r**2
      t53 = t4*s34r
      t57 = l3*t1
      t58 = t23*t4
      t70 = t23*t2*t24*l3
      t71 = s34r**2
      t72 = t4*t71
      t122 = t2*t4
      t131 = t23*t1*t4
      s1 = (t12*s1vr+(-512.D0/3.D0*t3*t15+512.D0/3.D0*t3*t18+512.D0/3.D0
     #*t7-512.D0/3.D0*t10)*s23r-512.D0/3.D0*t26*t28+512.D0/3.D0*t26*t32+
     #(512.D0/3.D0*t25*t6-512.D0/3.D0*t25*t9)*s34r)*s13r+(512.D0/3.D0*t4
     #3*t14*t44+512.D0/3.D0*t48*t5*t49+(-512.D0/3.D0*t43*t53*s2l2+(-512.
     #D0/3.D0*t57*t58*s3l2+512.D0/3.D0*t57*t58*s4l2)*s34r)*s23r-512.D0/3
     #.D0*t70*t72*s2l2-512.D0/3.D0*t70*t14*t71)*s1l1+(-512.D0/3.D0*t43*t
     #17*t44-512.D0/3.D0*t48*t8*t49+(512.D0/3.D0*t43*t53*s2l1+(512.D0/3.
     #D0*t57*t58*s3l1-512.D0/3.D0*t57*t58*s4l1)*s34r)*s23r+512.D0/3.D0*t
     #70*t72*s2l1+512.D0/3.D0*t70*t17*t71)*s1l2+((512.D0/3.D0*t48*t15-51
     #2.D0/3.D0*t48*t18)*s23r+512.D0/3.D0*t43*t14*s34r*s2l1-512.D0/3.D0*
     #t43*t17*s34r*s2l2)*s1vr
      taux6 = s1+(512.D0/3.D0*t122*s3l2*s2l1-512.D0/3.D0*t122*s3l1*s2l2+
     #512.D0/3.D0*t6-512.D0/3.D0*t9)*s23r+512.D0/3.D0*t131*t28-512.D0/3.
     #D0*t131*t32-t12*s34r

CCC Coefficients

      tcoff17 = taux5
      tcoff18 = taux5
      tcoff19 = taux5
      tcoff20 = taux5
      tcoff21 = taux6
      tcoff22 = taux6
      tcoff23 = -taux6
      tcoff24 = -taux6

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet17 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet19 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet20 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet23 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet24 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      sq3num3e4 = tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet2
     #0*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*
     #tcoff24
      RETURN
      END

      double precision function sq3num4e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = l3*t3
      t5 = l3-1
      t6 = t5**2
      t8 = 1/s134r
      t12 = r**2
      t13 = 1+t12
      t15 = 1/(r+t)
      t19 = 1/(1+r*t)
      t20 = t13*t15*t19
      t23 = l3*t2
      t28 = z*t8
      t33 = s13r**2
      t35 = t5*t8
      t37 = t23*t35*s2l2
      t38 = l3*t5
      t39 = t2*t8
      t41 = t38*t39*s3l2
      t47 = l3*t1
      t51 = t38*t2
      t52 = t8*s3l2
      t53 = t52*s34r
      t56 = t47*t52
      t63 = t38*t39*s3l1
      t72 = t8*s3l1
      t73 = t72*s34r
      t76 = t47*t72
      t80 = l3**2
      t89 = t23*t35*s4l2
      t95 = t23*t35*s4l1
      t104 = z*l3
      t110 = t35*s34r
      t121 = (l3-t104-1)*t8
      t127 = (-z+2-2*l3+2*t104)*t8
      t145 = t8*s4l2
      t148 = (128.D0/3.D0*t52+64.D0/3.D0*t145)*s3l1
      t149 = t8*s4l1
      t150 = t149*s3l2
      t151 = 64.D0/3.D0*t150
      t154 = t80*t2
      t163 = -64.D0/3.D0*t154*t52-64.D0/3.D0*t154*t145
      t168 = r*t13*t15*t19
      t170 = t8*s34r
      t173 = t47*t145
      t184 = -64.D0/3.D0*t154*t72-64.D0/3.D0*t154*t149
      t190 = t47*t149
      s2 = (-64.D0/3.D0*t4*t6*z*t8*l1*r*t20+64.D0/3.D0*t23*z*t5*t8+64.D0
     #/3.D0*t4*t6*t28*s34r)*t33+(((-64.D0/3.D0*t37+64.D0/3.D0*t41)*l1*r*
     #t20+64.D0/3.D0*t47*t8*s2l2-64.D0/3.D0*t51*t53-64.D0/3.D0*t56)*s1l1
     #+((-64.D0/3.D0*t23*t35*s2l1+64.D0/3.D0*t63)*l1*r*t20+64.D0/3.D0*t4
     #7*t8*s2l1-64.D0/3.D0*t51*t73-64.D0/3.D0*t76)*s1l2+(64.D0/3.D0*t80*
     #t5*t3*t28*s23r+(-128.D0/3.D0*t37-64.D0/3.D0*t41-128.D0/3.D0*t89)*s
     #2l1+(-64.D0/3.D0*t63-128.D0/3.D0*t95)*s2l2)*l1*r*t20-64.D0/3.D0*t4
     #7*z*(-l3+t104+1-z)*t8*s23r+((128.D0/3.D0*t23*t110+128.D0/3.D0*t5*t
     #1*t8)*s2l2+(64.D0/3.D0*t41+64.D0/3.D0*t89)*s34r-64.D0/3.D0*t121*s3
     #l2+64.D0/3.D0*t127*s4l2)*s2l1+((64.D0/3.D0*t63+64.D0/3.D0*t95)*s34
     #r-64.D0/3.D0*t121*s3l1+64.D0/3.D0*t127*s4l1)*s2l2+64.D0/3.D0*t23*z
     #*t110+t148+t151)*s13r
      s1 = s2+(-128.D0/3.D0*t154*t8*s34r*s23r*s1l2+t163*s23r*l1*t168+(64
     #.D0/3.D0*t154*t170*s2l2+64.D0/3.D0*t56+64.D0/3.D0*t173)*s23r-64.D0
     #/3.D0*t47*t53)*s1l1+(t184*s23r*l1*t168+(64.D0/3.D0*t154*t170*s2l1+
     #64.D0/3.D0*t76+64.D0/3.D0*t190)*s23r-64.D0/3.D0*t47*t73)*s1l2+((-t
     #163*s2l1-t184*s2l2)*s23r-128.D0/3.D0*t47*t145*s2l1-128.D0/3.D0*t47
     #*t149*s2l2+(-128.D0/3.D0*t56-64.D0/3.D0*t173)*s3l1-64.D0/3.D0*t47*
     #t150)*l1*r*t20
      taux7 = s1+((-64.D0/3.D0*t56-64.D0/3.D0*t173)*s2l1+(-64.D0/3.D0*t7
     #6-64.D0/3.D0*t190)*s2l2)*s23r+(64.D0/3.D0*t47*t145*s34r+128.D0/3.D
     #0*t145)*s2l1+(64.D0/3.D0*t47*t149*s34r+128.D0/3.D0*t149)*s2l2+t148
     #+t151
      t1 = l3-1
      t2 = -1+z
      t3 = t1*t2
      t4 = 1/s134r
      t5 = t4*s3l2
      t7 = t4*s4l2
      t12 = t4*s3l1
      t14 = t4*s4l1
      t20 = s13r**2
      t22 = t2**2
      t23 = l3*t22
      t24 = t1*t4
      t28 = r**2
      t29 = 1+t28
      t31 = 1/(r+t)
      t35 = 1/(1+r*t)
      t36 = t29*t31*t35
      t37 = l1*r*t36
      t40 = l3*t2
      t41 = t40*t5
      t49 = -l3*t2
      t52 = -256.D0/3.D0*t23*t24*s34r+256.D0/3.D0*t49*t4
      t60 = t40*t12
      t69 = l3*t1
      t70 = t22*t4
      t75 = -t69*t70*s3l2-t23*t24*s4l2
      t82 = t69*t70*s3l1+t23*t24*s4l1
      t95 = -z+2-2*l3+2*z*l3
      t96 = t95*t4
      t100 = -t95*t4
      t116 = l3**2
      t117 = t116*t22
      t121 = -256.D0/3.D0*t117*t5-256.D0/3.D0*t117*t7
      t126 = r*t29*t31*t35
      t128 = t4*s34r
      t131 = t49*t7
      t143 = 256.D0/3.D0*t117*t12+256.D0/3.D0*t117*t14
      t149 = t49*t14
      t170 = t7*s3l1
      t173 = t14*s3l2
      s2 = -((256.D0/3.D0*t3*t5+256.D0/3.D0*t3*t7)*s2l1+(-256.D0/3.D0*t3
     #*t12-256.D0/3.D0*t3*t14)*s2l2)*t20-((256.D0/3.D0*t23*t24*s2l2*t37+
     #(256.D0/3.D0*t41+256.D0/3.D0*t40*t7)*s23r+t52*s2l2)*s1l1+(-256.D0/
     #3.D0*t23*t24*s2l1*t37+(-256.D0/3.D0*t60-256.D0/3.D0*t40*t14)*s23r-
     #t52*s2l1)*s1l2+(512.D0/3.D0*t75*s2l1+512.D0/3.D0*t82*s2l2)*l1*r*t3
     #6+(-256.D0/3.D0*t75*s34r+256.D0/3.D0*t96*s3l2-256.D0/3.D0*t100*s4l
     #2)*s2l1+(-256.D0/3.D0*t82*s34r-256.D0/3.D0*t96*s3l1+256.D0/3.D0*t1
     #00*s4l1)*s2l2)*s13r
      s1 = s2-(t121*s23r*l1*t126+(-256.D0/3.D0*t117*t128*s2l2+256.D0/3.D
     #0*t41-256.D0/3.D0*t131)*s23r+256.D0/3.D0*t40*t5*s34r)*s1l1-(t143*s
     #23r*l1*t126+(256.D0/3.D0*t117*t128*s2l1-256.D0/3.D0*t60+256.D0/3.D
     #0*t149)*s23r-256.D0/3.D0*t40*t12*s34r)*s1l2-((-t121*s2l1-t143*s2l2
     #)*s23r+512.D0/3.D0*t49*t7*s2l1-512.D0/3.D0*t49*t14*s2l2+256.D0/3.D
     #0*t49*t170-256.D0/3.D0*t49*t173)*l1*r*t36
      taux8 = s1-((-256.D0/3.D0*t41+256.D0/3.D0*t131)*s2l1+(256.D0/3.D0*
     #t60-256.D0/3.D0*t149)*s2l2)*s23r-(-256.D0/3.D0*t49*t7*s34r+512.D0/
     #3.D0*t7)*s2l1-(256.D0/3.D0*t49*t14*s34r-512.D0/3.D0*t14)*s2l2-256.
     #D0/3.D0*t170+256.D0/3.D0*t173

CCC Coefficients

      tcoff25 = taux7
      tcoff26 = taux7
      tcoff27 = taux7
      tcoff28 = taux7
      tcoff29 = taux8
      tcoff30 = taux8
      tcoff31 = -taux8
      tcoff32 = -taux8

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet26 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet27 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet28 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet29 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet30 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet31 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet32 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      sq3num4e0 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet2
     #8*tcoff28+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet32*
     #tcoff32
      RETURN
      END

      double precision function sq3num4e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = t2*l3
      t4 = t3*z
      t5 = l3-1
      t6 = 1/s134r
      t7 = t5*t6
      t11 = t5*t1
      t12 = t6*s2l2
      t14 = 128.D0/3.D0*t11*t12
      t15 = t6*s3l2
      t18 = t6*s4l2
      t23 = t6*s3l1
      t26 = t6*s4l1
      t32 = s13r**2
      t35 = t3*t7*s2l2
      t36 = l3*t5
      t37 = t2*t6
      t39 = t36*t37*s3l2
      t44 = r**2
      t45 = 1+t44
      t46 = r+t
      t50 = 1+r*t
      t52 = t45/t46/t50
      t54 = t1*l3
      t56 = t54*t18
      t64 = -64.D0/3.D0*t3*t7*s34r-64.D0/3.D0*t54*t6
      t66 = t36*t2
      t67 = t15*s34r
      t70 = t54*t15
      t77 = t36*t37*s3l1
      t85 = t54*t26
      t90 = t23*s34r
      t97 = l3**2
      t98 = t2**2
      t99 = t97*t98
      t100 = t5**2
      t104 = l1**2
      t106 = t45**2
      t107 = t46**2
      t110 = t50**2
      t112 = t106/t107/t110
      t116 = t97*t5
      t117 = t2*t1
      t118 = t116*t117
      t119 = z*t6
      t120 = t119*s23r
      t136 = t99*z
      t137 = t100*t6
      t141 = t18*s3l1
      t142 = t66*t141
      t144 = t26*s3l2
      t145 = t66*t144
      t147 = 4*l3
      t148 = z*l3
      t149 = 4*t148
      t159 = t119*s34r
      t164 = t54*t144
      t165 = 64.D0/3.D0*t164
      t169 = -l3+t148+1-z
      t177 = t169*t6
      t180 = 2*l3
      t181 = 2*t148
      t183 = (-z+2-t180+t181)*t6
      t193 = s34r**2
      t205 = l3-t148
      t206 = t205*t6
      s1 = ((64.D0/3.D0*t35+64.D0/3.D0*t39)*l1*r*t52+(-64.D0/3.D0*t54*t1
     #2-64.D0/3.D0*t56)*s23r+t64*s2l2-64.D0/3.D0*t66*t67-64.D0/3.D0*t70)
     #*s1l1+((64.D0/3.D0*t3*t7*s2l1+64.D0/3.D0*t77)*l1*r*t52+(-64.D0/3.D
     #0*t54*t6*s2l1-64.D0/3.D0*t85)*s23r+t64*s2l1-64.D0/3.D0*t66*t90-64.
     #D0/3.D0*t54*t23)*s1l2-128.D0/3.D0*t99*z*t100*t6*t104*t44*t112+(128
     #.D0/3.D0*t118*t120+(128.D0/3.D0*t35+64.D0/3.D0*t39+128.D0/3.D0*t3*
     #t7*s4l2)*s2l1+(64.D0/3.D0*t77+128.D0/3.D0*t3*t7*s4l1)*s2l2+256.D0/
     #3.D0*t136*t137*s34r+64.D0/3.D0*t142+64.D0/3.D0*t145+64.D0/3.D0*t4*
     #t5*(-t147+t149+3-z)*t6)*l1*r*t52+(-128.D0/3.D0*t118*t159-64.D0/3.D
     #0*t54*t141-t165-128.D0/3.D0*t54*t18*s4l1-128.D0/3.D0*t54*z*t169*t6
     #)*s23r+(-t14-64.D0/3.D0*t177*s3l2-64.D0/3.D0*t183*s4l2)*s2l1
      t219 = s1+(-64.D0/3.D0*t177*s3l1-64.D0/3.D0*t183*s4l1)*s2l2-128.D0
     #/3.D0*t136*t137*t193+(-64.D0/3.D0*t142-64.D0/3.D0*t145-64.D0/3.D0*
     #t66*z*(-t147+t149-z+1)*t6)*s34r+(128.D0/3.D0*t15+64.D0/3.D0*t206*s
     #4l2)*s3l1+64.D0/3.D0*t206*s4l1*s3l2-64.D0/3.D0*t148*t1*(-t180+t181
     #+3-z)*t6
      t221 = t97*t2
      t224 = t118*t67
      t230 = t221*t15
      t231 = t221*t18
      t235 = 64.D0/3.D0*t56
      t238 = t6*t193
      t245 = -t205
      t246 = t54*t245
      t253 = t118*t90
      t282 = t118*t15*s3l1
      t285 = t116*t117*z*t6
      t291 = t18*s2l1
      t294 = t26*s2l2
      t325 = t1*t6
      t334 = t2*(t180-1)*t148*t6
      t357 = -t39-t56
      s2 = (-128.D0/3.D0*t4*t7*s23r+(-t14-128.D0/3.D0*t11*t15-64.D0/3.D0
     #*t11*t18)*s2l1+(-128.D0/3.D0*t11*t23-64.D0/3.D0*t11*t26)*s2l2)*t32
     #+t219*s13r+((64.D0/3.D0*t221*t18*s23r-64.D0/3.D0*t224)*l1*r*t52+((
     #64.D0/3.D0*t230-64.D0/3.D0*t231)*s34r-t235)*s23r+64.D0/3.D0*t118*t
     #238*s2l2+64.D0/3.D0*t118*t15*t193+64.D0/3.D0*t246*t67)*s1l1
      s1 = s2+((64.D0/3.D0*t221*t26*s23r-64.D0/3.D0*t253)*l1*r*t52+((64.
     #D0/3.D0*t221*t23-64.D0/3.D0*t221*t26)*s34r-64.D0/3.D0*t85)*s23r+64
     #.D0/3.D0*t118*t238*s2l1+64.D0/3.D0*t118*t23*t193+64.D0/3.D0*t246*t
     #90)*s1l2+(128.D0/3.D0*t5*t97*l3*t98*t120-128.D0/3.D0*t282-128.D0/3
     #.D0*t285)*t104*t44*t112+((-64.D0/3.D0*t221*t291-64.D0/3.D0*t221*t2
     #94+(128.D0/3.D0*t230+64.D0/3.D0*t231)*s3l1+64.D0/3.D0*t221*t144+64
     #.D0/3.D0*(-1+3*z+t147-t149)*t2*t97*z*t6)*s23r+(64.D0/3.D0*t224+128
     #.D0/3.D0*t56)*s2l1+(64.D0/3.D0*t253+128.D0/3.D0*t85)*s2l2+(128.D0/
     #3.D0*t282+128.D0/3.D0*t285)*s34r+(128.D0/3.D0*(1-z-t180+t181)*l3*t
     #325*s3l2+t235)*s3l1+t165+128.D0/3.D0*t334)*l1*r*t52+(64.D0/3.D0*t5
     #4*t291+64.D0/3.D0*t54*t294-128.D0/3.D0*t221*t159+(-128.D0/3.D0*t70
     #-t235)*s3l1-t165-64.D0/3.D0*(-t181+t180+1+z)*l3*z*t1*t6)*s23r
      taux7 = s1+(64.D0/3.D0*t357*s34r-128.D0/3.D0*t18)*s2l1+((-64.D0/3.
     #D0*t77-64.D0/3.D0*t85)*s34r-128.D0/3.D0*t26)*s2l2-128.D0/3.D0*t118
     #*t119*t193+(128.D0/3.D0*t357*s3l1-128.D0/3.D0*t164-64.D0/3.D0*t334
     #)*s34r+(-128.D0/3.D0*t245*t6*s3l2-64.D0/3.D0*t18)*s3l1-64.D0/3.D0*
     #t144-128.D0/3.D0*t148*t325
      t1 = l3-1
      t2 = -1+z
      t3 = t2*t1
      t4 = 1/s134r
      t5 = t4*s3l2
      t8 = t4*s4l2
      t13 = t4*s3l1
      t16 = t4*s4l1
      t22 = s13r**2
      t24 = t2**2
      t25 = l3*t24
      t26 = t1*t4
      t30 = t1*l3
      t31 = t24*t4
      t33 = t30*t31*s3l2
      t38 = r**2
      t46 = (1+t38)/(r+t)/(1+r*t)
      t48 = t2*l3
      t59 = 4*t2*l3
      t62 = 512.D0/3.D0*t25*t26*s34r+256.D0/3.D0*t59*t4
      t64 = t30*t24
      t65 = t5*s34r
      t69 = 256.D0/3.D0*t48*t5
      t76 = t30*t31*s3l1
      t90 = t13*s34r
      t94 = 256.D0/3.D0*t48*t13
      t99 = t25*t26*s4l2
      t105 = t25*t26*s4l1
      t109 = t8*s3l1
      t110 = t64*t109
      t112 = t16*s3l2
      t113 = t64*t112
      t122 = -256.D0/3.D0*t48*t109+256.D0/3.D0*t48*t112
      t129 = z*l3
      t132 = (z+4*l3-4*t129-4)*t4
      t139 = (-4*z+8-8*l3+8*t129)*t4
      t159 = (2*l3-2*t129-2)*t4
      t168 = l3**2
      t169 = t168*t24
      t170 = t169*t5
      t171 = 256.D0/3.D0*t170
      t172 = t169*t8
      t173 = 1024.D0/3.D0*t172
      t178 = t168*t1*t24*t2
      t180 = 256.D0/3.D0*t178*t65
      t185 = t4*s34r
      t192 = t59*t8
      t196 = s34r**2
      t197 = t4*t196
      t207 = (-t48*(-l3+t129)+t48)*t4
      t213 = t169*t13
      t214 = 256.D0/3.D0*t213
      t215 = t169*t16
      t216 = 1024.D0/3.D0*t215
      t220 = 256.D0/3.D0*t178*t90
      t231 = t59*t16
      t272 = -t2
      t274 = 3*t272*l3
      t283 = t272*l3
      s2 = -((-256.D0/3.D0*t3*t5-512.D0/3.D0*t3*t8)*s2l1+(256.D0/3.D0*t3
     #*t13+512.D0/3.D0*t3*t16)*s2l2)*t22-(((-1024.D0/3.D0*t25*t26*s2l2-2
     #56.D0/3.D0*t33)*l1*r*t46+(256.D0/3.D0*t48*t4*s2l2-256.D0/3.D0*t48*
     #t8)*s23r+t62*s2l2+256.D0/3.D0*t64*t65+t69)*s1l1+((1024.D0/3.D0*t25
     #*t26*s2l1+256.D0/3.D0*t76)*l1*r*t46+(-256.D0/3.D0*t48*t4*s2l1+256.
     #D0/3.D0*t48*t16)*s23r-t62*s2l1-256.D0/3.D0*t64*t90-t94)*s1l2+((102
     #4.D0/3.D0*t33+2048.D0/3.D0*t99)*s2l1+(-1024.D0/3.D0*t76-2048.D0/3.
     #D0*t105)*s2l2+512.D0/3.D0*t110-512.D0/3.D0*t113)*l1*r*t46+t122*s23
     #r+((-512.D0/3.D0*t33-1280.D0/3.D0*t99)*s34r+256.D0/3.D0*t132*s3l2-
     #256.D0/3.D0*t139*s4l2)*s2l1+((512.D0/3.D0*t76+1280.D0/3.D0*t105)*s
     #34r-256.D0/3.D0*t132*s3l1+256.D0/3.D0*t139*s4l1)*s2l2+(-512.D0/3.D
     #0*t110+512.D0/3.D0*t113)*s34r+256.D0/3.D0*t159*s4l2*s3l1-256.D0/3.
     #D0*t159*s4l1*s3l2)*s13r
      s1 = s2-(((t171+t173)*s23r-t180)*l1*r*t46+(256*t169*t185*s2l2+(t17
     #1-512.D0/3.D0*t172)*s34r-t69-256.D0/3.D0*t192)*s23r+256.D0/3.D0*t1
     #78*t197*s2l2+256.D0/3.D0*t178*t5*t196-256.D0/3.D0*t207*s3l2*s34r)*
     #s1l1-(((-t214-t216)*s23r+t220)*l1*r*t46+(-256*t169*t185*s2l1+(-t21
     #4+512.D0/3.D0*t215)*s34r+t94+256.D0/3.D0*t231)*s23r-256.D0/3.D0*t1
     #78*t197*s2l1-256.D0/3.D0*t178*t13*t196+256.D0/3.D0*t207*s3l1*s34r)
     #*s1l2-(((-256*t170-t173)*s2l1+(256*t213+t216)*s2l2-256.D0/3.D0*t16
     #9*t109+256.D0/3.D0*t169*t112)*s23r+(-t180+512.D0/3.D0*t192)*s2l1+(
     #t220-512.D0/3.D0*t231)*s2l2+256.D0/3.D0*t59*t109-256.D0/3.D0*t59*t
     #112)*l1*r*t46
      taux8 = s1-((-256.D0/3.D0*t274*t5+256.D0/3.D0*t192)*s2l1+(256.D0/3
     #.D0*t274*t13-256.D0/3.D0*t231)*s2l2-256.D0/3.D0*t283*t109+256.D0/3
     #.D0*t283*t112)*s23r-((256.D0/3.D0*t33-256.D0/3.D0*t192)*s34r-2048.
     #D0/3.D0*t8)*s2l1-((-256.D0/3.D0*t76+256.D0/3.D0*t231)*s34r+2048.D0
     #/3.D0*t16)*s2l2-t122*s34r+1024.D0/3.D0*t109-1024.D0/3.D0*t112

CCC Coefficients

      tcoff25 = taux7
      tcoff26 = taux7
      tcoff27 = taux7
      tcoff28 = taux7
      tcoff29 = taux8
      tcoff30 = taux8
      tcoff31 = -taux8
      tcoff32 = -taux8

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = t2*l3
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet26 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet27 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet28 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet29 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet30 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet31 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet32 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      sq3num4e1 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet2
     #8*tcoff28+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet32*
     #tcoff32
      RETURN
      END

      double precision function sq3num4e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = l3*t3
      t5 = l3-1
      t6 = t5**2
      t8 = 1/s134r
      t9 = t6*z*t8
      t12 = r**2
      t13 = 1+t12
      t14 = r+t
      t18 = 1+r*t
      t20 = t13/t14/t18
      t21 = l1*r*t20
      t24 = l3*t2
      t26 = t5*t8
      t30 = t5*t1
      t31 = t8*s4l2
      t35 = t8*s4l1
      t39 = t5*z
      t44 = z*t8
      t45 = t44*s34r
      t49 = s13r**2
      t51 = l3*t5
      t52 = t2*t8
      t54 = t51*t52*s3l2
      t56 = l3*t1
      t57 = s3l2*t8
      t60 = t24*t5
      t61 = t8*s34r
      t64 = t57*s34r
      t66 = t56*t57
      t71 = t51*t52*s3l1
      t73 = s3l1*t8
      t78 = t73*s34r
      t84 = l3**2
      t85 = t2**2
      t88 = l1**2
      t90 = t13**2
      t91 = t14**2
      t94 = t18**2
      t96 = t90/t91/t94
      t100 = t84*t5
      t101 = t100*t3
      t102 = t44*s23r
      t105 = t57*s3l1
      t107 = 128.D0/3.D0*t60*t105
      t109 = 4*l3
      t110 = z*l3
      t111 = 4*t110
      t125 = 64.D0/3.D0*t56*t31
      t127 = (-128.D0/3.D0*t66-t125)*s3l1
      t128 = t35*s3l2
      t130 = 64.D0/3.D0*t56*t128
      t131 = 2*l3
      t132 = 2*t110
      t133 = t131-t132+z-3
      t154 = s34r**2
      t155 = t44*t154
      t165 = -l3+t110
      t166 = t165*t8
      t178 = (64.D0/3.D0*t54*t21-64.D0/3.D0*t56*t57*s23r-64.D0/3.D0*t60*
     #t61*s2l2-64.D0/3.D0*t60*t64-64.D0/3.D0*t66)*s1l1+(64.D0/3.D0*t71*t
     #21-64.D0/3.D0*t56*t73*s23r-64.D0/3.D0*t60*t61*s2l1-64.D0/3.D0*t60*
     #t78-64.D0/3.D0*t56*t73)*s1l2+128.D0/3.D0*t84*t85*t9*t88*t12*t96+(-
     #128.D0/3.D0*t101*t102+t107-64.D0/3.D0*t39*l3*t2*(-t109+t111+3-z)*t
     #8)*l1*r*t20+(-64*t101*t45+t127-t130-64.D0/3.D0*t56*z*t133*t8)*s23r
     #+(-64.D0/3.D0*t54+64.D0/3.D0*t24*t26*s4l2)*s34r*s2l1+(-64.D0/3.D0*
     #t71+64.D0/3.D0*t24*t26*s4l1)*s34r*s2l2-128.D0/3.D0*t84*t6*t85*t155
     #+(-t107-64.D0/3.D0*t60*z*(z-4)*t8)*s34r+(-128.D0/3.D0*t166*s3l2+64
     #.D0/3.D0*t31)*s3l1+64.D0/3.D0*t128-64.D0/3.D0*t110*t1*t133*t8
      t180 = t3*t8
      t185 = s34r*l1*r*t20
      t187 = t84*t2
      t188 = t187*t8
      t192 = t8*t154
      t197 = t56*t165
      t216 = t84*l3
      t218 = t5*t216*t85
      t220 = t101*t105
      t221 = t3*z
      t223 = t100*t221*t8
      t230 = s23r**2
      t254 = t57*s34r*s2l1
      t258 = t73*s34r*s2l2
      t272 = t2*(t131-1)*t110*t8
      t290 = z*t1
      s2 = (-64*t4*t9*t21+128.D0/3.D0*t24*z*t26*s23r-64.D0/3.D0*t30*t31*
     #s2l1-64.D0/3.D0*t30*t35*s2l2+64*t24*t39*t8+64*t4*t6*t45)*t49+t178*
     #s13r+(-64.D0/3.D0*t100*t180*s3l2*t185+64.D0/3.D0*t188*s3l2*s34r*s2
     #3r+64.D0/3.D0*t101*t192*s2l2+64.D0/3.D0*t101*t57*t154+64.D0/3.D0*t
     #197*t64)*s1l1
      s1 = s2+(-64.D0/3.D0*t100*t180*s3l1*t185+64.D0/3.D0*t188*s3l1*s34r
     #*s23r+64.D0/3.D0*t101*t192*s2l1+64.D0/3.D0*t101*t73*t154+64.D0/3.D
     #0*t197*t78)*s1l2+(-128.D0/3.D0*t218*t102-128.D0/3.D0*t220+128.D0/3
     #.D0*t223)*t88*t12*t96+(64.D0/3.D0*t221*t216*t8*t230+(256.D0/3.D0*t
     #218*t45+(128.D0/3.D0*t187*t57+64.D0/3.D0*t187*t31)*s3l1+64.D0/3.D0
     #*t187*t128+64.D0/3.D0*(1-3*z-t109+t111)*t2*t84*z*t8)*s23r+64.D0/3.
     #D0*t101*t254+64.D0/3.D0*t101*t258+(128.D0/3.D0*t220-128.D0/3.D0*t2
     #23)*s34r+128.D0/3.D0*(1-z-t131+t132)*l3*t1*t105-128.D0/3.D0*t272)*
     #l1*r*t20-64.D0/3.D0*t187*t44*t230
      taux7 = s1+(-64.D0/3.D0*t187*z*(-t109+t111-4*z+3)*t8*s34r+t127-t13
     #0-64.D0/3.D0*(t132-t131-1-z)*l3*t290*t8)*s23r-64.D0/3.D0*t60*t254-
     #64.D0/3.D0*t60*t258-128.D0/3.D0*t101*t155+((-128.D0/3.D0*t54-t125)
     #*s3l1-t130+64.D0/3.D0*t272)*s34r-128.D0/3.D0*t166*s3l2*s3l1+128.D0
     #/3.D0*t290*l3*t8
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = 1/s134r
      t5 = t4*s3l2
      t8 = t4*s4l2
      t13 = t4*s3l1
      t16 = t4*s4l1
      t22 = s13r**2
      t24 = t1**2
      t25 = t24*l3
      t26 = t2*t4
      t30 = l3*t2
      t31 = t24*t4
      t33 = t30*t31*s3l2
      t38 = r**2
      t39 = 1+t38
      t41 = 1/(r+t)
      t45 = 1/(1+r*t)
      t46 = t39*t41*t45
      t48 = l3*t1
      t52 = t48*t5
      t63 = -5*l3*t1
      t66 = -512.D0/3.D0*t25*t26*s34r+256.D0/3.D0*t63*t4
      t68 = t30*t24
      t69 = t5*s34r
      t79 = t30*t31*s3l1
      t88 = t48*t13
      t96 = t13*s34r
      t102 = t8*s2l1
      t105 = t16*s2l2
      t108 = t8*s3l1
      t109 = t68*t108
      t111 = t16*s3l2
      t112 = t68*t111
      t120 = t48*t108-t48*t111
      t129 = z*t4
      t134 = z*l3
      t137 = (5*z-10+10*l3-10*t134)*t4
      t160 = (3-4*l3+4*t134)*t4
      t169 = l3**2
      t170 = t169*t24
      t172 = 512.D0/3.D0*t170*t5
      t174 = 1280.D0/3.D0*t170*t8
      t180 = r*t39*t41*t45
      t182 = t4*s34r
      t186 = t8*s34r
      t189 = 512.D0/3.D0*t52
      t191 = 256.D0/3.D0*t63*t8
      t199 = 512.D0/3.D0*t170*t13
      t201 = 1280.D0/3.D0*t170*t16
      t209 = t16*s34r
      t212 = 512.D0/3.D0*t88
      t214 = 256.D0/3.D0*t63*t16
      s3 = -((-256.D0/3.D0*t3*t5+512.D0/3.D0*t3*t8)*s2l1+(256.D0/3.D0*t3
     #*t13-512.D0/3.D0*t3*t16)*s2l2)*t22
      s4 = -(((1280.D0/3.D0*t25*t26*s2l2+1024.D0/3.D0*t33)*l1*r*t46+(-25
     #6*t48*t4*s2l2-1280.D0/3.D0*t52-512.D0/3.D0*t48*t8)*s23r+t66*s2l2-1
     #024.D0/3.D0*t68*t69-1024.D0/3.D0*t52)*s1l1+((-1280.D0/3.D0*t25*t26
     #*s2l1-1024.D0/3.D0*t79)*l1*r*t46+(256*t48*t4*s2l1+1280.D0/3.D0*t88
     #+512.D0/3.D0*t48*t16)*s23r-t66*s2l1+1024.D0/3.D0*t68*t96+1024.D0/3
     #.D0*t88)*s1l2+(-2560.D0/3.D0*t68*t102+2560.D0/3.D0*t68*t105-1024.D
     #0/3.D0*t109+1024.D0/3.D0*t112)*l1*r*t46+256*t120*s23r+((256.D0/3.D
     #0*t33+1280.D0/3.D0*t25*t26*s4l2)*s34r+512.D0/3.D0*t129*s3l2-256.D0
     #/3.D0*t137*s4l2)*s2l1+((-256.D0/3.D0*t79-1280.D0/3.D0*t25*t26*s4l1
     #)*s34r-512.D0/3.D0*t129*s3l1+256.D0/3.D0*t137*s4l1)*s2l2+(1024.D0/
     #3.D0*t109-1024.D0/3.D0*t112)*s34r+256.D0/3.D0*t160*s4l2*s3l1-256.D
     #0/3.D0*t160*s4l1*s3l2)*s13r
      s2 = s3+s4
      s1 = s2-((t172-t174)*s23r*l1*t180+(-512.D0/3.D0*t170*t182*s2l2+102
     #4.D0/3.D0*t170*t186-t189-t191)*s23r-512.D0/3.D0*t48*t69)*s1l1-((-t
     #199+t201)*s23r*l1*t180+(512.D0/3.D0*t170*t182*s2l1-1024.D0/3.D0*t1
     #70*t209+t212+t214)*s23r+512.D0/3.D0*t48*t96)*s1l2-(((t172+t174)*s2
     #l1+(-t199-t201)*s2l2)*s23r+512.D0/3.D0*t63*t102-512.D0/3.D0*t63*t1
     #05+256.D0/3.D0*t63*t108-256.D0/3.D0*t63*t111)*l1*r*t46
      taux8 = s1-((-t189+t191)*s2l1+(t212-t214)*s2l2)*s23r-(-256.D0/3.D0
     #*t63*t186+2560.D0/3.D0*t8)*s2l1-(256.D0/3.D0*t63*t209-2560.D0/3.D0
     #*t16)*s2l2-1024.D0/3.D0*t120*s34r-1280.D0/3.D0*t108+1280.D0/3.D0*t
     #111

CCC Coefficients

      tcoff25 = taux7
      tcoff26 = taux7
      tcoff27 = taux7
      tcoff28 = taux7
      tcoff29 = taux8
      tcoff30 = taux8
      tcoff31 = -taux8
      tcoff32 = -taux8

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet26 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet27 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet28 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet29 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet30 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet31 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet32 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      sq3num4e2 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet2
     #8*tcoff28+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet32*
     #tcoff32
      RETURN
      END

      double precision function sq3num4e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = -1+z
      t2 = t1**2
      t3 = t1*t2
      t4 = l3*t3
      t5 = l3-1
      t6 = t5**2
      t8 = 1/s134r
      t12 = r**2
      t20 = (1+t12)/(r+t)/(1+r*t)
      t23 = l3*t2
      t28 = z*t8
      t29 = t28*s34r
      t33 = s13r**2
      t35 = l3**2
      t37 = t35*t5*t3
      t41 = t2**2
      t42 = t35*t6*t41
      t52 = z*l3
      t53 = -l3+t52+1-z
      t60 = s34r**2
      t61 = t28*t60
      taux7 = (-64.D0/3.D0*t4*t6*z*t8*l1*r*t20+64.D0/3.D0*t23*t5*z*t8+64
     #.D0/3.D0*t4*t6*t29)*t33+((64.D0/3.D0*t37*t28*s23r+64.D0/3.D0*t42*t
     #29)*l1*r*t20+(-128.D0/3.D0*t37*t29-64.D0/3.D0*l3*t1*z*t53*t8)*s23r
     #-64.D0/3.D0*t42*t61-64.D0/3.D0*t23*t5*z*(-l3+t52-1)*t8*s34r)*s13r+
     #64.D0/3.D0*t35*l3*t41*t5*t29*s23r*l1*r*t20-64.D0/3.D0*t2*t35*z*t53
     #*t8*s34r*s23r-64.D0/3.D0*t37*t61
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = 1/s134r
      t5 = t4*s3l2
      t8 = t4*s4l2
      t13 = t4*s3l1
      t16 = t4*s4l1
      t22 = s13r**2
      t24 = t1**2
      t25 = t24*l3
      t26 = t2*t4
      t30 = l3*t2
      t31 = t24*t4
      t33 = t30*t31*s3l2
      t34 = 256*t33
      t38 = r**2
      t46 = (1+t38)/(r+t)/(1+r*t)
      t48 = l3*t1
      t52 = t48*t5
      t61 = 256*t25*t26*s34r+512.D0/3.D0*t48*t4
      t63 = t30*t24
      t64 = t5*s34r
      t74 = t30*t31*s3l1
      t75 = 256*t74
      t83 = t48*t13
      t89 = t13*s34r
      t95 = t8*s2l1
      t97 = t16*s2l2
      t104 = t8*s3l1
      t105 = t48*t104
      t106 = t16*s3l2
      t107 = t48*t106
      t108 = -t105+t107
      t117 = z*l3
      t120 = (-z+2-2*l3+2*t117)*t4
      t138 = l3**2
      t139 = t24*t138
      t145 = t138*t2*t24*t1
      t147 = 256*t145*t64
      t154 = t48*t8
      t155 = 512.D0/3.D0*t154
      t158 = s34r**2
      t159 = t4*t158
      t167 = t48*(-l3+t117)
      t176 = 256*t145*t89
      t183 = t48*t16
      t184 = 512.D0/3.D0*t183
      t224 = 3*l3*t1
      s2 = -((-512.D0/3.D0*t3*t5-256*t3*t8)*s2l1+(512.D0/3.D0*t3*t13+256
     #*t3*t16)*s2l2)*t22-(((-512.D0/3.D0*t25*t26*s2l2-t34)*l1*r*t46+(512
     #.D0/3.D0*t48*t4*s2l2+256.D0/3.D0*t52)*s23r+t61*s2l2+256*t63*t64+25
     #6*t52)*s1l1+((512.D0/3.D0*t25*t26*s2l1+t75)*l1*r*t46+(-512.D0/3.D0
     #*t48*t4*s2l1-256.D0/3.D0*t83)*s23r-t61*s2l1-256*t63*t89-256*t83)*s
     #1l2+(1024.D0/3.D0*t63*t95-1024.D0/3.D0*t63*t97)*l1*r*t46+256.D0/3.
     #D0*t108*s23r+((-256.D0/3.D0*t33+256.D0/3.D0*t25*t26*s4l2)*s34r-512
     #.D0/3.D0*t120*s4l2)*s2l1+((256.D0/3.D0*t74-256.D0/3.D0*t25*t26*s4l
     #1)*s34r+512.D0/3.D0*t120*s4l1)*s2l2+256.D0/3.D0*t104-256.D0/3.D0*t
     #106)*s13r
      s1 = s2-((512.D0/3.D0*t139*t8*s23r+t147)*l1*r*t46+(-256*t139*t64-t
     #155)*s23r-256*t145*t159*s2l2-256*t145*t5*t158-256*t167*t64)*s1l1-(
     #(-512.D0/3.D0*t139*t16*s23r-t176)*l1*r*t46+(256*t139*t89+t184)*s23
     #r+256*t145*t159*s2l1+256*t145*t13*t158+256*t167*t89)*s1l2-((-512.D
     #0/3.D0*t139*t95+512.D0/3.D0*t139*t97+256*t139*t104-256*t139*t106)*
     #s23r+(t147+1024.D0/3.D0*t154)*s2l1+(-t176-1024.D0/3.D0*t183)*s2l2+
     #512.D0/3.D0*t105-512.D0/3.D0*t107)*l1*r*t46
      taux8 = s1-(512.D0/3.D0*t48*t95-512.D0/3.D0*t48*t97-256.D0/3.D0*t2
     #24*t104+256.D0/3.D0*t224*t106)*s23r-((-t34-t155)*s34r-1024.D0/3.D0
     #*t8)*s2l1-((t75+t184)*s34r+1024.D0/3.D0*t16)*s2l2-256*t108*s34r+51
     #2.D0/3.D0*t104-512.D0/3.D0*t106

CCC Coefficients

      tcoff25 = taux7
      tcoff26 = taux7
      tcoff27 = taux7
      tcoff28 = taux7
      tcoff29 = taux8
      tcoff30 = taux8
      tcoff31 = -taux8
      tcoff32 = -taux8

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet26 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet27 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet28 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet29 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet30 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet31 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet32 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      sq3num4e3 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet2
     #8*tcoff28+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet32*
     #tcoff32
      RETURN
      END

      double precision function sq3num4e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux7 = 0
      t1 = l3-1
      t2 = -1+z
      t3 = t1*t2
      t4 = 1/s134r
      t5 = t4*s4l2
      t8 = t4*s4l1
      t13 = s13r**2
      t15 = l3*t1
      t16 = t2**2
      t17 = t16*t4
      t19 = t15*t17*s3l2
      t21 = r**2
      t29 = (1+t21)/(r+t)/(1+r*t)
      t30 = l1*r*t29
      t32 = l3*t2
      t33 = s3l2*t4
      t36 = t16*l3
      t37 = t36*t1
      t38 = t4*s34r
      t41 = t33*s34r
      t48 = t15*t17*s3l1
      t50 = s3l1*t4
      t55 = t50*s34r
      t61 = t5*s3l1
      t63 = t8*s3l2
      t66 = -512.D0/3.D0*t32*t61+512.D0/3.D0*t32*t63
      t67 = t66*s23r
      t68 = t1*t4
      t85 = l3**2
      t86 = t85*t1
      t87 = t16*t2
      t88 = t87*t4
      t93 = s34r*l1*r*t29
      t95 = t85*t16
      t96 = t95*t4
      t100 = t86*t87
      t101 = s34r**2
      t102 = t4*t101
      t109 = t32*(-l3+z*l3)
      t134 = t33*s34r*s2l1
      t138 = t50*s34r*s2l2
      s1 = -(-512.D0/3.D0*t3*t5*s2l1+512.D0/3.D0*t3*t8*s2l2)*t13-((-512.
     #D0/3.D0*t19*t30+512.D0/3.D0*t32*t33*s23r+512.D0/3.D0*t37*t38*s2l2+
     #512.D0/3.D0*t37*t41+512.D0/3.D0*t32*t33)*s1l1+(512.D0/3.D0*t48*t30
     #-512.D0/3.D0*t32*t50*s23r-512.D0/3.D0*t37*t38*s2l1-512.D0/3.D0*t37
     #*t55-512.D0/3.D0*t32*t50)*s1l2+t67+(-512.D0/3.D0*t19+512.D0/3.D0*t
     #36*t68*s4l2)*s34r*s2l1+(512.D0/3.D0*t48-512.D0/3.D0*t36*t68*s4l1)*
     #s34r*s2l2+512.D0/3.D0*t61-512.D0/3.D0*t63)*s13r-(512.D0/3.D0*t86*t
     #88*s3l2*t93-512.D0/3.D0*t96*s3l2*s34r*s23r-512.D0/3.D0*t100*t102*s
     #2l2-512.D0/3.D0*t100*t33*t101-512.D0/3.D0*t109*t41)*s1l1-(-512.D0/
     #3.D0*t86*t88*s3l1*t93+512.D0/3.D0*t96*s3l1*s34r*s23r+512.D0/3.D0*t
     #100*t102*s2l1+512.D0/3.D0*t100*t50*t101+512.D0/3.D0*t109*t55)*s1l2
      taux8 = s1-((512.D0/3.D0*t95*t61-512.D0/3.D0*t95*t63)*s23r+512.D0/
     #3.D0*t100*t134-512.D0/3.D0*t100*t138)*l1*r*t29-t67+512.D0/3.D0*t37
     #*t134-512.D0/3.D0*t37*t138-t66*s34r

CCC Coefficients

      tcoff25 = taux7
      tcoff26 = taux7
      tcoff27 = taux7
      tcoff28 = taux7
      tcoff29 = taux8
      tcoff30 = taux8
      tcoff31 = -taux8
      tcoff32 = -taux8

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet26 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet27 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet28 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet29 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet30 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet31 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet32 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      sq3num4e4 = tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tFjet2
     #8*tcoff28+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet32*
     #tcoff32
      RETURN
      END

      double precision function sq3num5e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t7 = t1**2
      t8 = t3*t2
      t9 = t7*t8
      t10 = l3*s34r
      t15 = l3*t3
      t16 = t1*s2l2
      t17 = t15*t16
      t18 = t1*s3l2
      t19 = t15*t18
      t25 = t1*s4l2
      t26 = t15*t25
      t27 = 512.D0/3.D0*t26
      t30 = 1+z
      t31 = t2*t30
      t34 = l3*s3l2
      t38 = t1*t2
      t45 = t1*s2l1
      t46 = t15*t45
      t47 = t1*s3l1
      t48 = t15*t47
      t54 = t1*s4l1
      t55 = t15*t54
      t71 = l3**2
      t72 = z*t71
      t73 = t3**2
      t74 = t73*t7
      t77 = 1024.D0/3.D0*t72*t74*s34r
      t82 = t15*t54*s3l2
      t83 = 512.D0/3.D0*t82
      t84 = l3*z
      t85 = -t2
      t92 = t19+t26
      t109 = l3*s4l2
      t113 = t2*z
      t116 = z*s4l2
      t121 = l3*s4l1
      t125 = z*s4l1
      t133 = s34r**2
      t156 = t71*t1*t8*s34r
      t157 = t84*t2
      t161 = t71*t8
      t162 = t1*s34r
      t166 = t71*t7
      t168 = t166*t73*t133
      t171 = t84*t4*s34r
      t175 = t71*t3
      t176 = t175*s2l2
      t177 = t175*s3l2
      t180 = s1vr**2
      t182 = t175*s4l2
      t186 = l3*t2
      t188 = 1024.D0/3.D0*t156+1024.D0/3.D0*t186
      t190 = 2*z
      t191 = t190-1
      t199 = -512.D0/3.D0*t156-1024.D0/3.D0*t157
      t204 = -t191
      t215 = -512.D0/3.D0*t168+512.D0/3.D0*t171+512.D0/3.D0*z*t30
      t221 = l3*t191
      t227 = z*s3l2
      t229 = 512.D0/3.D0*t116
      t232 = t175*s2l1
      t233 = t175*s3l1
      t269 = z*s3l1
      t276 = s3l2*s3l1
      t293 = s3l2*s4l1
      t308 = 512.D0/3.D0*t186*s4l2
      t313 = t186*s4l1
      t320 = t8*t1
      t353 = z**2
      s5 = ((-1024.D0/3.D0*t4*l3*s1vr+1024.D0/3.D0*t9*t10)*s1l2+(-512.D0
     #/3.D0*t17-512.D0/3.D0*t19)*s1vr+(1024.D0/3.D0*t17+512.D0/3.D0*t19+
     #t27)*s23r-512.D0/3.D0*t31*t16+512.D0/3.D0*t9*t34*s34r-512.D0/3.D0*
     #t38*s3l2-512.D0/3.D0*t38*s4l2)*s1l1+((-512.D0/3.D0*t46-512.D0/3.D0
     #*t48)*s1vr+(1024.D0/3.D0*t46+512.D0/3.D0*t48+512.D0/3.D0*t55)*s23r
     #-512.D0/3.D0*t31*t45+512.D0/3.D0*t9*l3*s3l1*s34r-512.D0/3.D0*t38*s
     #3l1-512.D0/3.D0*t38*s4l1)*s1l2+(-t77+(-1024.D0/3.D0*t19-t27)*s3l1-
     #t83+512.D0/3.D0*t84*t3*t85*t1)*s1vr+(512.D0/3.D0*t92*s2l1+(512.D0/
     #3.D0*t48+512.D0/3.D0*t55)*s2l2+t77+1024.D0/3.D0*t92*s3l1+1024.D0/3
     #.D0*t82+1024.D0/3.D0*t15*t25*s4l1)*s23r
      s4 = s5+(-1024.D0/3.D0*t9*t10*s2l2-512.D0/3.D0*t9*t109*s34r-512.D0
     #/3.D0*t113*t18-512.D0/3.D0*t38*t116)*s2l1+(-512.D0/3.D0*t9*t121*s3
     #4r-512.D0/3.D0*t113*t47-512.D0/3.D0*t38*t125)*s2l2+1024.D0/3.D0*t7
     #2*t73*t2*t7*t1*t133+((1024.D0/3.D0*t9*t34+512.D0/3.D0*t9*t109)*s3l
     #1+512.D0/3.D0*t9*t121*s3l2+512.D0/3.D0*t84*t9*t30)*s34r
      s5 = s13r
      s3 = s4*s5
      s4 = (((1024.D0/3.D0*t156+1024.D0/3.D0*t157)*s1vr-1024.D0/3.D0*t16
     #1*t162*s23r-1024.D0/3.D0*t168-1024.D0/3.D0*t171)*s1l2+(-512.D0/3.D
     #0*t176+512.D0/3.D0*t177)*t180+((-512.D0/3.D0*t176-512.D0/3.D0*t177
     #-512.D0/3.D0*t182)*s23r+t188*s2l2+512.D0/3.D0*t186*t191*s3l2)*s1vr
     #+(t199*s2l2-512.D0/3.D0*t161*t18*s34r+512.D0/3.D0*t186*t204*s3l2-5
     #12.D0/3.D0*t186*t191*s4l2)*s23r+t215*s2l2-512.D0/3.D0*t166*t73*s3l
     #2*t133+(-512.D0/3.D0*t4*t221*s3l2+512.D0/3.D0*t26)*s34r+512.D0/3.D
     #0*t227+t229)*s1l1
      s2 = s3+s4
      s3 = s2+((-512.D0/3.D0*t232+512.D0/3.D0*t233)*t180+((-512.D0/3.D0*
     #t232-512.D0/3.D0*t233-512.D0/3.D0*t175*s4l1)*s23r+t188*s2l1+512.D0
     #/3.D0*t186*t191*s3l1)*s1vr+(t199*s2l1-512.D0/3.D0*t161*t47*s34r+51
     #2.D0/3.D0*t186*t204*s3l1-512.D0/3.D0*t186*t191*s4l1)*s23r+t215*s2l
     #1-512.D0/3.D0*t166*t73*s3l1*t133+(-512.D0/3.D0*t4*t221*s3l1+512.D0
     #/3.D0*t55)*s34r+512.D0/3.D0*t269+512.D0/3.D0*t125)*s1l2
      s1 = s3+(-1024.D0/3.D0*t175*s2l2*s2l1+1024.D0/3.D0*t175*t276+1024.
     #D0/3.D0*t72*t3)*t180+(((1024.D0/3.D0*t176-512.D0/3.D0*t177)*s2l1-5
     #12.D0/3.D0*t175*s3l1*s2l2+(-1024.D0/3.D0*t177-512.D0/3.D0*t182)*s3
     #l1-512.D0/3.D0*t175*t293+512.D0/3.D0*t72*t8)*s23r+((1024.D0/3.D0*t
     #156-1024.D0/3.D0*t221*t2)*s2l2+512.D0/3.D0*t186*t30*s3l2+t308)*s2l
     #1+(512.D0/3.D0*t186*t30*s3l1+512.D0/3.D0*t313)*s2l2+(-1024.D0/3.D0
     #*t161*t18*s3l1-2048.D0/3.D0*t72*t320)*s34r+1024.D0/3.D0*t15*t276-1
     #024.D0/3.D0*t84*t2*t85)*s1vr
      taux9 = s1+((1024.D0/3.D0*t15*s2l2-512.D0/3.D0*t186*t227-t308)*s2l
     #1+(-512.D0/3.D0*t186*t269-512.D0/3.D0*t313)*s2l2+1024.D0/3.D0*t72*
     #t320*s34r+(-1024.D0/3.D0*t15*s3l2-512.D0/3.D0*t15*s4l2)*s3l1-512.D
     #0/3.D0*t15*t293+512.D0/3.D0*t84*t2*(t353-t190+1))*s23r+((1024.D0/3
     #.D0*t8*l3*t162-1024.D0/3.D0*t113)*s2l2-512.D0/3.D0*t15*t25*s34r+51
     #2.D0/3.D0*t353*s3l2+t229)*s2l1+(-512.D0/3.D0*t15*t54*s34r+512.D0/3
     #.D0*t353*s3l1+512.D0/3.D0*t125)*s2l2+1024.D0/3.D0*t72*t74*t133+((-
     #1024.D0/3.D0*t320*t34+t27)*s3l1+t83+512.D0/3.D0*t84*t3*(-3*z+1)*t1
     #)*s34r
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t7 = t1**2
      t8 = t3*t2
      t9 = t7*t8
      t10 = l3*s34r
      t15 = l3*t3
      t16 = t1*s2l2
      t17 = t15*t16
      t18 = t1*s3l2
      t19 = t15*t18
      t25 = t1*s4l2
      t26 = t15*t25
      t27 = 512.D0/3.D0*t26
      t30 = 1+z
      t31 = t2*t30
      t34 = l3*s3l2
      t38 = t1*t2
      t45 = t1*s2l1
      t46 = t15*t45
      t47 = t1*s3l1
      t48 = t15*t47
      t54 = t1*s4l1
      t55 = t15*t54
      t71 = l3**2
      t72 = z*t71
      t73 = t3**2
      t74 = t73*t7
      t77 = 1024.D0/3.D0*t72*t74*s34r
      t82 = t15*t54*s3l2
      t83 = 512.D0/3.D0*t82
      t84 = l3*z
      t85 = -t2
      t92 = t19+t26
      t109 = l3*s4l2
      t113 = t2*z
      t116 = z*s4l2
      t121 = l3*s4l1
      t125 = z*s4l1
      t133 = s34r**2
      t156 = t71*t1*t8*s34r
      t157 = t84*t2
      t161 = t71*t8
      t162 = t1*s34r
      t166 = t71*t7
      t168 = t166*t73*t133
      t171 = t84*t4*s34r
      t175 = t71*t3
      t176 = t175*s2l2
      t177 = t175*s3l2
      t180 = s1vr**2
      t182 = t175*s4l2
      t186 = l3*t2
      t188 = 1024.D0/3.D0*t156+1024.D0/3.D0*t186
      t190 = 2*z
      t191 = t190-1
      t199 = -512.D0/3.D0*t156-1024.D0/3.D0*t157
      t204 = -t191
      t215 = -512.D0/3.D0*t168+512.D0/3.D0*t171+512.D0/3.D0*z*t30
      t221 = l3*t191
      t227 = z*s3l2
      t229 = 512.D0/3.D0*t116
      t232 = t175*s2l1
      t233 = t175*s3l1
      t269 = z*s3l1
      t276 = s3l2*s3l1
      t293 = s3l2*s4l1
      t308 = 512.D0/3.D0*t186*s4l2
      t313 = t186*s4l1
      t320 = t8*t1
      t353 = z**2
      s5 = ((-1024.D0/3.D0*t4*l3*s1vr+1024.D0/3.D0*t9*t10)*s1l2+(-512.D0
     #/3.D0*t17-512.D0/3.D0*t19)*s1vr+(1024.D0/3.D0*t17+512.D0/3.D0*t19+
     #t27)*s23r-512.D0/3.D0*t31*t16+512.D0/3.D0*t9*t34*s34r-512.D0/3.D0*
     #t38*s3l2-512.D0/3.D0*t38*s4l2)*s1l1+((-512.D0/3.D0*t46-512.D0/3.D0
     #*t48)*s1vr+(1024.D0/3.D0*t46+512.D0/3.D0*t48+512.D0/3.D0*t55)*s23r
     #-512.D0/3.D0*t31*t45+512.D0/3.D0*t9*l3*s3l1*s34r-512.D0/3.D0*t38*s
     #3l1-512.D0/3.D0*t38*s4l1)*s1l2+(-t77+(-1024.D0/3.D0*t19-t27)*s3l1-
     #t83+512.D0/3.D0*t84*t3*t85*t1)*s1vr+(512.D0/3.D0*t92*s2l1+(512.D0/
     #3.D0*t48+512.D0/3.D0*t55)*s2l2+t77+1024.D0/3.D0*t92*s3l1+1024.D0/3
     #.D0*t82+1024.D0/3.D0*t15*t25*s4l1)*s23r
      s4 = s5+(-1024.D0/3.D0*t9*t10*s2l2-512.D0/3.D0*t9*t109*s34r-512.D0
     #/3.D0*t113*t18-512.D0/3.D0*t38*t116)*s2l1+(-512.D0/3.D0*t9*t121*s3
     #4r-512.D0/3.D0*t113*t47-512.D0/3.D0*t38*t125)*s2l2+1024.D0/3.D0*t7
     #2*t73*t2*t7*t1*t133+((1024.D0/3.D0*t9*t34+512.D0/3.D0*t9*t109)*s3l
     #1+512.D0/3.D0*t9*t121*s3l2+512.D0/3.D0*t84*t9*t30)*s34r
      s5 = s13r
      s3 = s4*s5
      s4 = (((1024.D0/3.D0*t156+1024.D0/3.D0*t157)*s1vr-1024.D0/3.D0*t16
     #1*t162*s23r-1024.D0/3.D0*t168-1024.D0/3.D0*t171)*s1l2+(-512.D0/3.D
     #0*t176+512.D0/3.D0*t177)*t180+((-512.D0/3.D0*t176-512.D0/3.D0*t177
     #-512.D0/3.D0*t182)*s23r+t188*s2l2+512.D0/3.D0*t186*t191*s3l2)*s1vr
     #+(t199*s2l2-512.D0/3.D0*t161*t18*s34r+512.D0/3.D0*t186*t204*s3l2-5
     #12.D0/3.D0*t186*t191*s4l2)*s23r+t215*s2l2-512.D0/3.D0*t166*t73*s3l
     #2*t133+(-512.D0/3.D0*t4*t221*s3l2+512.D0/3.D0*t26)*s34r+512.D0/3.D
     #0*t227+t229)*s1l1
      s2 = s3+s4
      s3 = s2+((-512.D0/3.D0*t232+512.D0/3.D0*t233)*t180+((-512.D0/3.D0*
     #t232-512.D0/3.D0*t233-512.D0/3.D0*t175*s4l1)*s23r+t188*s2l1+512.D0
     #/3.D0*t186*t191*s3l1)*s1vr+(t199*s2l1-512.D0/3.D0*t161*t47*s34r+51
     #2.D0/3.D0*t186*t204*s3l1-512.D0/3.D0*t186*t191*s4l1)*s23r+t215*s2l
     #1-512.D0/3.D0*t166*t73*s3l1*t133+(-512.D0/3.D0*t4*t221*s3l1+512.D0
     #/3.D0*t55)*s34r+512.D0/3.D0*t269+512.D0/3.D0*t125)*s1l2
      s1 = s3+(-1024.D0/3.D0*t175*s2l2*s2l1+1024.D0/3.D0*t175*t276+1024.
     #D0/3.D0*t72*t3)*t180+(((1024.D0/3.D0*t176-512.D0/3.D0*t177)*s2l1-5
     #12.D0/3.D0*t175*s3l1*s2l2+(-1024.D0/3.D0*t177-512.D0/3.D0*t182)*s3
     #l1-512.D0/3.D0*t175*t293+512.D0/3.D0*t72*t8)*s23r+((1024.D0/3.D0*t
     #156-1024.D0/3.D0*t221*t2)*s2l2+512.D0/3.D0*t186*t30*s3l2+t308)*s2l
     #1+(512.D0/3.D0*t186*t30*s3l1+512.D0/3.D0*t313)*s2l2+(-1024.D0/3.D0
     #*t161*t18*s3l1-2048.D0/3.D0*t72*t320)*s34r+1024.D0/3.D0*t15*t276-1
     #024.D0/3.D0*t84*t2*t85)*s1vr
      taux10 = s1+((1024.D0/3.D0*t15*s2l2-512.D0/3.D0*t186*t227-t308)*s2
     #l1+(-512.D0/3.D0*t186*t269-512.D0/3.D0*t313)*s2l2+1024.D0/3.D0*t72
     #*t320*s34r+(-1024.D0/3.D0*t15*s3l2-512.D0/3.D0*t15*s4l2)*s3l1-512.
     #D0/3.D0*t15*t293+512.D0/3.D0*t84*t2*(t353-t190+1))*s23r+((1024.D0/
     #3.D0*t8*l3*t162-1024.D0/3.D0*t113)*s2l2-512.D0/3.D0*t15*t25*s34r+5
     #12.D0/3.D0*t353*s3l2+t229)*s2l1+(-512.D0/3.D0*t15*t54*s34r+512.D0/
     #3.D0*t353*s3l1+512.D0/3.D0*t125)*s2l2+1024.D0/3.D0*t72*t74*t133+((
     #-1024.D0/3.D0*t320*t34+t27)*s3l1+t83+512.D0/3.D0*t84*t3*(-3*z+1)*t
     #1)*s34r
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t7 = t4*s3l2
      t8 = t3*t7
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t25 = 1+z
      t29 = 4096.D0/3.D0*t21*l3*s34r-2048.D0/3.D0*t1*t25*t4
      t35 = t4*t1
      t36 = t35*s3l2
      t38 = t35*s4l2
      t43 = t3*t4*s2l1
      t44 = t4*s3l1
      t45 = t3*t44
      t51 = t4*s4l1
      t52 = t3*t51
      t62 = t35*s3l1
      t64 = t35*s4l1
      t73 = 2048.D0/3.D0*t3*t14*s3l1-2048.D0/3.D0*t3*t51*s3l2
      t83 = l3*s4l2
      t89 = l3*s4l1
      t104 = l3**2
      t105 = t104*t2
      t106 = t105*s2l2
      t107 = t105*s3l2
      t110 = s1vr**2
      t118 = t104*t4*t20*s34r
      t119 = z*l3
      t120 = t119*t1
      t122 = -4096.D0/3.D0*t118+4096.D0/3.D0*t120
      t124 = t104*t20
      t126 = t124*t7*s34r
      t128 = l3*t1
      t129 = z*s3l2
      t130 = t128*t129
      t131 = 2048.D0/3.D0*t130
      t136 = 2048.D0/3.D0*t118-4096.D0/3.D0*t120
      t139 = t1*s4l2
      t144 = t104*t19
      t145 = t2**2
      t146 = s34r**2
      t150 = 2*z-1
      t156 = 2048.D0/3.D0*t144*t145*t146-2048.D0/3.D0*t3*t150*t4*s34r+20
     #48.D0/3.D0*z*t25
      t162 = z*t4
      t169 = z*s4l2
      t173 = t105*s2l1
      t174 = t105*s3l1
      t185 = t124*t44*s34r
      t187 = z*s3l1
      t188 = t128*t187
      t189 = 2048.D0/3.D0*t188
      t195 = t1*s4l1
      t212 = z*s4l1
      t220 = s4l2*s3l1
      t222 = s3l2*s4l1
      t229 = t128*s4l2
      t235 = t128*s4l1
      s3 = (((-2048.D0/3.D0*t6-2048.D0/3.D0*t8)*s1vr+(4096.D0/3.D0*t6+20
     #48.D0/3.D0*t8+2048.D0/3.D0*t15)*s23r+t29*s2l2+2048.D0/3.D0*t21*l3*
     #s3l2*s34r-2048.D0/3.D0*t36-2048.D0/3.D0*t38)*s1l1+((2048.D0/3.D0*t
     #43+2048.D0/3.D0*t45)*s1vr+(-4096.D0/3.D0*t43-2048.D0/3.D0*t45-2048
     #.D0/3.D0*t52)*s23r-t29*s2l1-2048.D0/3.D0*t21*l3*s3l1*s34r+2048.D0/
     #3.D0*t62+2048.D0/3.D0*t64)*s1l2+t73*s1vr+((-2048.D0/3.D0*t8-2048.D
     #0/3.D0*t15)*s2l1+(2048.D0/3.D0*t45+2048.D0/3.D0*t52)*s2l2)*s23r+(-
     #2048.D0/3.D0*t21*t83*s34r+2048.D0/3.D0*t36+2048.D0/3.D0*t38)*s2l1+
     #(2048.D0/3.D0*t21*t89*s34r-2048.D0/3.D0*t62-2048.D0/3.D0*t64)*s2l2
     #+(-2048.D0/3.D0*t21*t83*s3l1+2048.D0/3.D0*t21*t89*s3l2)*s34r)*s13r
      s4 = ((2048.D0/3.D0*t106+2048.D0/3.D0*t107)*t110+((-2048.D0/3.D0*t
     #106-2048.D0/3.D0*t107-2048.D0/3.D0*t105*s4l2)*s23r+t122*s2l2-4096.
     #D0/3.D0*t126+t131)*s1vr+(t136*s2l2+2048.D0/3.D0*t126-t131-2048.D0/
     #3.D0*t119*t139)*s23r+t156*s2l2+2048.D0/3.D0*t144*t145*s3l2*t146+(-
     #2048.D0/3.D0*t3*t162*s3l2+2048.D0/3.D0*t15)*s34r+2048.D0/3.D0*t129
     #+2048.D0/3.D0*t169)*s1l1
      s2 = s3+s4
      s1 = s2+((-2048.D0/3.D0*t173-2048.D0/3.D0*t174)*t110+((2048.D0/3.D
     #0*t173+2048.D0/3.D0*t174+2048.D0/3.D0*t105*s4l1)*s23r-t122*s2l1+40
     #96.D0/3.D0*t185-t189)*s1vr+(-t136*s2l1-2048.D0/3.D0*t185+t189+2048
     #.D0/3.D0*t119*t195)*s23r-t156*s2l1-2048.D0/3.D0*t144*t145*s3l1*t14
     #6+(2048.D0/3.D0*t3*t162*s3l1-2048.D0/3.D0*t52)*s34r-2048.D0/3.D0*t
     #187-2048.D0/3.D0*t212)*s1l2+((2048.D0/3.D0*t105*s3l2*s2l1-2048.D0/
     #3.D0*t105*s3l1*s2l2-2048.D0/3.D0*t105*t220+2048.D0/3.D0*t105*t222)
     #*s23r+(-4096.D0/3.D0*t128*s3l2-2048.D0/3.D0*t229)*s2l1+(4096.D0/3.
     #D0*t128*s3l1+2048.D0/3.D0*t235)*s2l2)*s1vr+((2048.D0/3.D0*t130+204
     #8.D0/3.D0*t229)*s2l1+(-2048.D0/3.D0*t188-2048.D0/3.D0*t235)*s2l2-2
     #048.D0/3.D0*t3*t220+2048.D0/3.D0*t3*t222)*s23r
      taux11 = s1+(2048.D0/3.D0*t3*t14*s34r-2048.D0/3.D0*t150*s3l2-2048.
     #D0/3.D0*t169)*s2l1+(-2048.D0/3.D0*t3*t51*s34r+2048.D0/3.D0*t150*s3
     #l1+2048.D0/3.D0*t212)*s2l2+t73*s34r+2048.D0/3.D0*t139*s3l1-2048.D0
     #/3.D0*t195*s3l2
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t7 = t4*s3l2
      t8 = t3*t7
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t25 = 1+z
      t29 = -4096.D0/3.D0*t21*l3*s34r+2048.D0/3.D0*t1*t25*t4
      t35 = t4*t1
      t36 = t35*s3l2
      t38 = t35*s4l2
      t43 = t3*t4*s2l1
      t44 = t4*s3l1
      t45 = t3*t44
      t51 = t4*s4l1
      t52 = t3*t51
      t62 = t35*s3l1
      t64 = t35*s4l1
      t73 = -2048.D0/3.D0*t3*t14*s3l1+2048.D0/3.D0*t3*t51*s3l2
      t83 = l3*s4l2
      t89 = l3*s4l1
      t104 = l3**2
      t105 = t104*t2
      t106 = t105*s2l2
      t107 = t105*s3l2
      t110 = s1vr**2
      t118 = t104*t4*t20*s34r
      t119 = z*l3
      t120 = t119*t1
      t122 = 4096.D0/3.D0*t118-4096.D0/3.D0*t120
      t124 = t104*t20
      t126 = t124*t7*s34r
      t128 = l3*t1
      t129 = z*s3l2
      t130 = t128*t129
      t131 = 2048.D0/3.D0*t130
      t136 = -2048.D0/3.D0*t118+4096.D0/3.D0*t120
      t139 = t1*s4l2
      t144 = t104*t19
      t145 = t2**2
      t146 = s34r**2
      t150 = 2*z-1
      t156 = -2048.D0/3.D0*t144*t145*t146+2048.D0/3.D0*t3*t150*t4*s34r-2
     #048.D0/3.D0*z*t25
      t162 = z*t4
      t169 = z*s4l2
      t173 = t105*s2l1
      t174 = t105*s3l1
      t185 = t124*t44*s34r
      t187 = z*s3l1
      t188 = t128*t187
      t189 = 2048.D0/3.D0*t188
      t195 = t1*s4l1
      t212 = z*s4l1
      t220 = s4l2*s3l1
      t222 = s3l2*s4l1
      t229 = t128*s4l2
      t235 = t128*s4l1
      s3 = (((2048.D0/3.D0*t6+2048.D0/3.D0*t8)*s1vr+(-4096.D0/3.D0*t6-20
     #48.D0/3.D0*t8-2048.D0/3.D0*t15)*s23r+t29*s2l2-2048.D0/3.D0*t21*l3*
     #s3l2*s34r+2048.D0/3.D0*t36+2048.D0/3.D0*t38)*s1l1+((-2048.D0/3.D0*
     #t43-2048.D0/3.D0*t45)*s1vr+(4096.D0/3.D0*t43+2048.D0/3.D0*t45+2048
     #.D0/3.D0*t52)*s23r-t29*s2l1+2048.D0/3.D0*t21*l3*s3l1*s34r-2048.D0/
     #3.D0*t62-2048.D0/3.D0*t64)*s1l2+t73*s1vr+((2048.D0/3.D0*t8+2048.D0
     #/3.D0*t15)*s2l1+(-2048.D0/3.D0*t45-2048.D0/3.D0*t52)*s2l2)*s23r+(2
     #048.D0/3.D0*t21*t83*s34r-2048.D0/3.D0*t36-2048.D0/3.D0*t38)*s2l1+(
     #-2048.D0/3.D0*t21*t89*s34r+2048.D0/3.D0*t62+2048.D0/3.D0*t64)*s2l2
     #+(2048.D0/3.D0*t21*t83*s3l1-2048.D0/3.D0*t21*t89*s3l2)*s34r)*s13r
      s4 = ((-2048.D0/3.D0*t106-2048.D0/3.D0*t107)*t110+((2048.D0/3.D0*t
     #106+2048.D0/3.D0*t107+2048.D0/3.D0*t105*s4l2)*s23r+t122*s2l2+4096.
     #D0/3.D0*t126-t131)*s1vr+(t136*s2l2-2048.D0/3.D0*t126+t131+2048.D0/
     #3.D0*t119*t139)*s23r+t156*s2l2-2048.D0/3.D0*t144*t145*s3l2*t146+(2
     #048.D0/3.D0*t3*t162*s3l2-2048.D0/3.D0*t15)*s34r-2048.D0/3.D0*t129-
     #2048.D0/3.D0*t169)*s1l1
      s2 = s3+s4
      s1 = s2+((2048.D0/3.D0*t173+2048.D0/3.D0*t174)*t110+((-2048.D0/3.D
     #0*t173-2048.D0/3.D0*t174-2048.D0/3.D0*t105*s4l1)*s23r-t122*s2l1-40
     #96.D0/3.D0*t185+t189)*s1vr+(-t136*s2l1+2048.D0/3.D0*t185-t189-2048
     #.D0/3.D0*t119*t195)*s23r-t156*s2l1+2048.D0/3.D0*t144*t145*s3l1*t14
     #6+(-2048.D0/3.D0*t3*t162*s3l1+2048.D0/3.D0*t52)*s34r+2048.D0/3.D0*
     #t187+2048.D0/3.D0*t212)*s1l2+((-2048.D0/3.D0*t105*s3l2*s2l1+2048.D
     #0/3.D0*t105*s3l1*s2l2+2048.D0/3.D0*t105*t220-2048.D0/3.D0*t105*t22
     #2)*s23r+(4096.D0/3.D0*t128*s3l2+2048.D0/3.D0*t229)*s2l1+(-4096.D0/
     #3.D0*t128*s3l1-2048.D0/3.D0*t235)*s2l2)*s1vr+((-2048.D0/3.D0*t130-
     #2048.D0/3.D0*t229)*s2l1+(2048.D0/3.D0*t188+2048.D0/3.D0*t235)*s2l2
     #+2048.D0/3.D0*t3*t220-2048.D0/3.D0*t3*t222)*s23r
      taux12 = s1+(-2048.D0/3.D0*t3*t14*s34r+2048.D0/3.D0*t150*s3l2+2048
     #.D0/3.D0*t169)*s2l1+(2048.D0/3.D0*t3*t51*s34r-2048.D0/3.D0*t150*s3
     #l1-2048.D0/3.D0*t212)*s2l2+t73*s34r-2048.D0/3.D0*t139*s3l1+2048.D0
     #/3.D0*t195*s3l2

CCC Coefficients

      tcoff33 = taux9
      tcoff34 = taux9
      tcoff35 = taux10
      tcoff36 = taux10
      tcoff37 = taux11
      tcoff38 = taux11
      tcoff39 = taux12
      tcoff40 = taux12

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet33 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet34 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff34)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet35 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet36 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet38 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet39 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet40 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff40)

CCC Sum

      sq3num5e0 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet3
     #6*tcoff36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet40*
     #tcoff40
      RETURN
      END

      double precision function sq3num5e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t11 = t6*t3
      t25 = t3**2
      t32 = s13r**2
      t34 = t3*l3
      t35 = t5*s2l2
      t37 = t5*s3l2
      t38 = t34*t37
      t39 = t5*t3
      t40 = l3*s4l2
      t46 = t7*l3*s3l2
      t47 = t7*t40
      t55 = t5*s3l1
      t56 = t34*t55
      t57 = l3*s4l1
      t63 = t7*l3*s3l1
      t64 = t7*t57
      t70 = l3**2
      t71 = z*t70
      t72 = t4*t5
      t73 = s1vr**2
      t82 = t71*t25*t6*s34r
      t84 = 3*z
      t93 = 1+z
      t108 = t2*z
      t111 = t5*t2
      t131 = s3l2*s3l1
      t134 = s4l1*s4l2
      t137 = z**2
      t141 = ((512.D0/3.D0*t34*t35-512.D0/3.D0*t38+512.D0/3.D0*t39*t40)*
     #s1vr+(512.D0/3.D0*t46-512.D0/3.D0*t47)*s34r)*s1l1+((512.D0/3.D0*t3
     #4*t5*s2l1-512.D0/3.D0*t56+512.D0/3.D0*t39*t57)*s1vr+(512.D0/3.D0*t
     #63-512.D0/3.D0*t64)*s34r)*s1l2-512*t71*t72*t73+(-1024.D0/3.D0*t34*
     #t35*s2l1+512*t82+512.D0/3.D0*t1*t3*(-t84+5)*t5)*s1vr+(-512.D0/3.D0
     #*t82-512*t1*t93*t3*t5)*s23r+((-1024.D0/3.D0*t7*l3*s34r-1024.D0/3.D
     #0*t39)*s2l2+(-512.D0/3.D0*t46-512.D0/3.D0*t47)*s34r-1024.D0/3.D0*t
     #108*t37+1024.D0/3.D0*t111*s4l2)*s2l1+((-512.D0/3.D0*t63-512.D0/3.D
     #0*t64)*s34r-1024.D0/3.D0*t108*t55+1024.D0/3.D0*t111*s4l1)*s2l2-512
     #.D0/3.D0*t1*t4*t6*(-t84+4)*s34r-1024.D0/3.D0*t111*t131+1024.D0/3.D
     #0*t111*t134+1024.D0/3.D0*t137*t2*t5
      t143 = t70*t3
      t150 = t70*t5*t4*s34r
      t151 = 1024.D0/3.D0*t150
      t154 = s23r**2
      t157 = t70*t4
      t158 = t5*s34r
      t165 = 512.D0/3.D0*t143*s2l2
      t166 = t143*s3l2
      t168 = t143*s4l2
      t177 = -t151+512.D0/3.D0*t34
      t179 = t157*t37
      t181 = t70*s4l2
      t182 = t72*t181
      t186 = l3*t2
      t187 = z-2
      t200 = 512.D0/3.D0*t150-512.D0/3.D0*t34
      t207 = z*s4l2
      t212 = t70*t6
      t213 = s34r**2
      t216 = t4*l3
      t219 = 512.D0/3.D0*t212*t25*t213-512.D0/3.D0*t216*t158
      t228 = l3*t187
      t239 = 512.D0/3.D0*t143*s2l1
      t240 = t143*s3l1
      t242 = t143*s4l1
      t251 = t157*t55
      t253 = t70*s4l1
      t254 = t72*t253
      t297 = t70*l3
      t298 = t297*t4
      t335 = 2*z
      t336 = t335-1
      t353 = 1024.D0/3.D0*t186*t131
      t355 = 1024.D0/3.D0*t186*t134
      s2 = (512*t1*t7*s1vr+(2048.D0/3.D0*t11*s2l2+1024.D0/3.D0*t11*s3l2+
     #1024.D0/3.D0*t11*s4l2)*s2l1+(1024.D0/3.D0*t11*s3l1+1024.D0/3.D0*t1
     #1*s4l1)*s2l2-512*t1*t25*t6*t5*s34r)*t32+t141*s13r+((1024.D0/3.D0*t
     #143*t73+(-1024*t143*s23r-t151)*s1vr+2048.D0/3.D0*t143*t154+2048.D0
     #/3.D0*t157*t158*s23r)*s1l2+(t165+1024.D0/3.D0*t166+512.D0/3.D0*t16
     #8)*t73+((-t165-2048.D0/3.D0*t166-2048.D0/3.D0*t168)*s23r+t177*s2l2
     #+(-512*t179-1024.D0/3.D0*t182)*s34r+512.D0/3.D0*t186*t187*s3l2+512
     #.D0/3.D0*t186*t93*s4l2)*s1vr+(1024.D0/3.D0*t166+1024.D0/3.D0*t168)
     #*t154+(t200*s2l2+(512*t179+512*t182)*s34r+1024.D0/3.D0*t186*s3l2-1
     #024.D0/3.D0*t186*t207)*s23r+t219*s2l2+(512.D0/3.D0*t212*t25*s3l2+5
     #12.D0/3.D0*t212*t25*s4l2)*t213+(-512.D0/3.D0*t39*t228*s3l2-512.D0/
     #3.D0*t39*t1*s4l2)*s34r)*s1l1
      s3 = s2+((t239+1024.D0/3.D0*t240+512.D0/3.D0*t242)*t73+((-t239-204
     #8.D0/3.D0*t240-2048.D0/3.D0*t242)*s23r+t177*s2l1+(-512*t251-1024.D
     #0/3.D0*t254)*s34r+512.D0/3.D0*t186*t187*s3l1+512.D0/3.D0*t186*t93*
     #s4l1)*s1vr+(1024.D0/3.D0*t240+1024.D0/3.D0*t242)*t154+(t200*s2l1+(
     #512*t251+512*t254)*s34r+1024.D0/3.D0*t186*s3l1-1024.D0/3.D0*t186*z
     #*s4l1)*s23r+t219*s2l1+(512.D0/3.D0*t212*t25*s3l1+512.D0/3.D0*t212*
     #t25*s4l1)*t213+(-512.D0/3.D0*t39*t228*s3l1-512.D0/3.D0*t39*t1*s4l1
     #)*s34r)*s1l2
      s1 = s3+(512*t298*z*s23r-512.D0/3.D0*t143*s4l2*s2l1-512.D0/3.D0*t1
     #43*s4l1*s2l2-1024.D0/3.D0*t71*t3)*t73+(-512*t298*z*t154+((-512.D0/
     #3.D0*t166+512.D0/3.D0*t168)*s2l1+(-512.D0/3.D0*t240+512.D0/3.D0*t2
     #42)*s2l2-1024.D0/3.D0*t297*t5*t25*z*s34r+512.D0/3.D0*t71*t3*(6*z-2
     #))*s23r+(512.D0/3.D0*t72*t181*s34r-512.D0/3.D0*t186*t336*s4l2)*s2l
     #1+(512.D0/3.D0*t72*t253*s34r-512.D0/3.D0*t186*t336*s4l1)*s2l2+512.
     #D0/3.D0*t71*t72*s34r-t353+t355-1024.D0/3.D0*t1*t2*t336)*s1vr-512*t
     #71*t4*t154
      taux9 = s1+((-512.D0/3.D0*t34*s3l2+512.D0/3.D0*t34*s4l2)*s2l1+(-51
     #2.D0/3.D0*t34*s3l1+512.D0/3.D0*t34*s4l1)*s2l2-512.D0/3.D0*t71*t4*(
     #t335-3)*t5*s34r+t353-t355+512.D0/3.D0*t1*t2*(3*t137-1))*s23r+((512
     #.D0/3.D0*t38+512.D0/3.D0*t216*t5*s4l2)*s34r-512.D0/3.D0*t108*s4l2)
     #*s2l1+((512.D0/3.D0*t56+512.D0/3.D0*t216*t5*s4l1)*s34r-512.D0/3.D0
     #*t108*s4l1)*s2l2+(1024.D0/3.D0*t34*t37*s3l1-1024.D0/3.D0*t39*t40*s
     #4l1+512.D0/3.D0*t1*t72)*s34r+1024.D0/3.D0*t131*z+1024.D0/3.D0*t207
     #*s4l1-1024.D0/3.D0*t137*z
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t11 = t6*t3
      t25 = t3**2
      t32 = s13r**2
      t34 = t3*l3
      t35 = t5*s2l2
      t37 = t5*s3l2
      t38 = t34*t37
      t39 = t5*t3
      t40 = l3*s4l2
      t46 = t7*l3*s3l2
      t47 = t7*t40
      t55 = t5*s3l1
      t56 = t34*t55
      t57 = l3*s4l1
      t63 = t7*l3*s3l1
      t64 = t7*t57
      t70 = l3**2
      t71 = z*t70
      t72 = t4*t5
      t73 = s1vr**2
      t82 = t71*t25*t6*s34r
      t84 = 3*z
      t93 = 1+z
      t108 = t2*z
      t111 = t5*t2
      t131 = s3l2*s3l1
      t134 = s4l1*s4l2
      t137 = z**2
      t141 = ((512.D0/3.D0*t34*t35-512.D0/3.D0*t38+512.D0/3.D0*t39*t40)*
     #s1vr+(512.D0/3.D0*t46-512.D0/3.D0*t47)*s34r)*s1l1+((512.D0/3.D0*t3
     #4*t5*s2l1-512.D0/3.D0*t56+512.D0/3.D0*t39*t57)*s1vr+(512.D0/3.D0*t
     #63-512.D0/3.D0*t64)*s34r)*s1l2-512*t71*t72*t73+(-1024.D0/3.D0*t34*
     #t35*s2l1+512*t82+512.D0/3.D0*t1*t3*(-t84+5)*t5)*s1vr+(-512.D0/3.D0
     #*t82-512*t1*t93*t3*t5)*s23r+((-1024.D0/3.D0*t7*l3*s34r-1024.D0/3.D
     #0*t39)*s2l2+(-512.D0/3.D0*t46-512.D0/3.D0*t47)*s34r-1024.D0/3.D0*t
     #108*t37+1024.D0/3.D0*t111*s4l2)*s2l1+((-512.D0/3.D0*t63-512.D0/3.D
     #0*t64)*s34r-1024.D0/3.D0*t108*t55+1024.D0/3.D0*t111*s4l1)*s2l2-512
     #.D0/3.D0*t1*t4*t6*(-t84+4)*s34r-1024.D0/3.D0*t111*t131+1024.D0/3.D
     #0*t111*t134+1024.D0/3.D0*t137*t2*t5
      t143 = t70*t3
      t150 = t70*t5*t4*s34r
      t151 = 1024.D0/3.D0*t150
      t154 = s23r**2
      t157 = t70*t4
      t158 = t5*s34r
      t165 = 512.D0/3.D0*t143*s2l2
      t166 = t143*s3l2
      t168 = t143*s4l2
      t177 = -t151+512.D0/3.D0*t34
      t179 = t157*t37
      t181 = t70*s4l2
      t182 = t72*t181
      t186 = t2*l3
      t187 = z-2
      t200 = 512.D0/3.D0*t150-512.D0/3.D0*t34
      t207 = z*s4l2
      t212 = t70*t6
      t213 = s34r**2
      t216 = t4*l3
      t219 = 512.D0/3.D0*t212*t25*t213-512.D0/3.D0*t216*t158
      t228 = l3*t187
      t239 = 512.D0/3.D0*t143*s2l1
      t240 = t143*s3l1
      t242 = t143*s4l1
      t251 = t157*t55
      t253 = t70*s4l1
      t254 = t72*t253
      t297 = t70*l3
      t298 = t297*t4
      t335 = 2*z
      t336 = t335-1
      t353 = 1024.D0/3.D0*t186*t131
      t355 = 1024.D0/3.D0*t186*t134
      s2 = (512*t1*t7*s1vr+(2048.D0/3.D0*t11*s2l2+1024.D0/3.D0*t11*s3l2+
     #1024.D0/3.D0*t11*s4l2)*s2l1+(1024.D0/3.D0*t11*s3l1+1024.D0/3.D0*t1
     #1*s4l1)*s2l2-512*t1*t25*t6*t5*s34r)*t32+t141*s13r+((1024.D0/3.D0*t
     #143*t73+(-1024*t143*s23r-t151)*s1vr+2048.D0/3.D0*t143*t154+2048.D0
     #/3.D0*t157*t158*s23r)*s1l2+(t165+1024.D0/3.D0*t166+512.D0/3.D0*t16
     #8)*t73+((-t165-2048.D0/3.D0*t166-2048.D0/3.D0*t168)*s23r+t177*s2l2
     #+(-512*t179-1024.D0/3.D0*t182)*s34r+512.D0/3.D0*t186*t187*s3l2+512
     #.D0/3.D0*t186*t93*s4l2)*s1vr+(1024.D0/3.D0*t166+1024.D0/3.D0*t168)
     #*t154+(t200*s2l2+(512*t179+512*t182)*s34r+1024.D0/3.D0*t186*s3l2-1
     #024.D0/3.D0*t186*t207)*s23r+t219*s2l2+(512.D0/3.D0*t212*t25*s3l2+5
     #12.D0/3.D0*t212*t25*s4l2)*t213+(-512.D0/3.D0*t39*t228*s3l2-512.D0/
     #3.D0*t39*t1*s4l2)*s34r)*s1l1
      s3 = s2+((t239+1024.D0/3.D0*t240+512.D0/3.D0*t242)*t73+((-t239-204
     #8.D0/3.D0*t240-2048.D0/3.D0*t242)*s23r+t177*s2l1+(-512*t251-1024.D
     #0/3.D0*t254)*s34r+512.D0/3.D0*t186*t187*s3l1+512.D0/3.D0*t186*t93*
     #s4l1)*s1vr+(1024.D0/3.D0*t240+1024.D0/3.D0*t242)*t154+(t200*s2l1+(
     #512*t251+512*t254)*s34r+1024.D0/3.D0*t186*s3l1-1024.D0/3.D0*t186*z
     #*s4l1)*s23r+t219*s2l1+(512.D0/3.D0*t212*t25*s3l1+512.D0/3.D0*t212*
     #t25*s4l1)*t213+(-512.D0/3.D0*t39*t228*s3l1-512.D0/3.D0*t39*t1*s4l1
     #)*s34r)*s1l2
      s1 = s3+(512*t298*z*s23r-512.D0/3.D0*t143*s4l2*s2l1-512.D0/3.D0*t1
     #43*s4l1*s2l2-1024.D0/3.D0*t71*t3)*t73+(-512*t298*z*t154+((-512.D0/
     #3.D0*t166+512.D0/3.D0*t168)*s2l1+(-512.D0/3.D0*t240+512.D0/3.D0*t2
     #42)*s2l2-1024.D0/3.D0*t297*t5*t25*z*s34r+512.D0/3.D0*t71*t3*(6*z-2
     #))*s23r+(512.D0/3.D0*t72*t181*s34r-512.D0/3.D0*t186*t336*s4l2)*s2l
     #1+(512.D0/3.D0*t72*t253*s34r-512.D0/3.D0*t186*t336*s4l1)*s2l2+512.
     #D0/3.D0*t71*t72*s34r-t353+t355-1024.D0/3.D0*t1*t2*t336)*s1vr-512*t
     #71*t4*t154
      taux10 = s1+((-512.D0/3.D0*t34*s3l2+512.D0/3.D0*t34*s4l2)*s2l1+(-5
     #12.D0/3.D0*t34*s3l1+512.D0/3.D0*t34*s4l1)*s2l2-512.D0/3.D0*t71*t4*
     #(t335-3)*t5*s34r+t353-t355+512.D0/3.D0*t1*t2*(3*t137-1))*s23r+((51
     #2.D0/3.D0*t38+512.D0/3.D0*t216*t5*s4l2)*s34r-512.D0/3.D0*t108*s4l2
     #)*s2l1+((512.D0/3.D0*t56+512.D0/3.D0*t216*t5*s4l1)*s34r-512.D0/3.D
     #0*t108*s4l1)*s2l2+(1024.D0/3.D0*t34*t37*s3l1-1024.D0/3.D0*t39*t40*
     #s4l1+512.D0/3.D0*t1*t72)*s34r+1024.D0/3.D0*t131*z+1024.D0/3.D0*t20
     #7*s4l1-1024.D0/3.D0*t137*z
      t1 = -1+z
      t2 = t1**2
      t3 = t2*l3
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t11 = t4*t2
      t12 = l3*s4l2
      t13 = t11*t12
      t22 = t4**2
      t23 = t2*t1
      t24 = t22*t23
      t28 = 1+z
      t32 = -8192.D0/3.D0*t24*l3*s34r+2048*t1*t28*t4
      t35 = t24*l3*s3l2
      t37 = t24*t12
      t41 = -t1
      t43 = 3*t41*t4
      t51 = t3*t4*s2l1
      t53 = t4*s3l1
      t54 = t3*t53
      t56 = l3*s4l1
      t57 = t11*t56
      t69 = t24*l3*s3l1
      t71 = t24*t56
      t81 = -t9-t13
      t84 = t54+t57
      t87 = t12*s3l1
      t88 = t11*t87
      t90 = t56*s3l2
      t91 = t11*t90
      t108 = 2*z
      t110 = (t41*z+2-t108)*t4
      t134 = l3**2
      t135 = t134*t2
      t137 = 2048*t135*s2l2
      t138 = t135*s3l2
      t141 = s1vr**2
      t144 = t135*s4l2
      t150 = t134*t4*t23*s34r
      t152 = z*l3
      t153 = t152*t1
      t154 = 4096*t153
      t155 = 8192.D0/3.D0*t150-t154
      t157 = t134*t23
      t159 = t157*t8*s34r
      t160 = 4096.D0/3.D0*t159
      t161 = l3*t1
      t162 = 2*t153
      t163 = -t161-t162
      t169 = -2048.D0/3.D0*t150+t154
      t175 = z*s4l2
      t180 = t134*t22
      t181 = t2**2
      t182 = s34r**2
      t187 = t108-1
      t196 = -2048.D0/3.D0*t180*t181*t182-2048.D0/3.D0*(l3*t23-3*t3*t187
     #)*t4*s34r-2048*z*t28
      t202 = z*t4
      t210 = 2048*t175
      t214 = 2048*t135*s2l1
      t215 = t135*s3l1
      t220 = t135*s4l1
      t227 = t157*t53*s34r
      t228 = 4096.D0/3.D0*t227
      t239 = z*s4l1
      t257 = 2048*t239
      t275 = s4l2*s3l1
      t278 = s3l2*s4l1
      t287 = -t161*(-1+3*z)-4*t161
      t304 = -l3*t187*t1-t162
      t307 = -t153-2*t161
      t323 = t187*t4
      t330 = z**2
      t332 = -t330-4*z+2
      s3 = (((2048*t6+4096.D0/3.D0*t9-2048.D0/3.D0*t13)*s1vr+(-4096*t6-2
     #048*t9-2048*t13)*s23r+t32*s2l2+(-4096.D0/3.D0*t35+2048.D0/3.D0*t37
     #)*s34r-2048.D0/3.D0*t43*s3l2-2048.D0/3.D0*t43*s4l2)*s1l1+((-2048*t
     #51-4096.D0/3.D0*t54+2048.D0/3.D0*t57)*s1vr+(4096*t51+2048*t54+2048
     #*t57)*s23r-t32*s2l1+(4096.D0/3.D0*t69-2048.D0/3.D0*t71)*s34r+2048.
     #D0/3.D0*t43*s3l1+2048.D0/3.D0*t43*s4l1)*s1l2+(2048.D0/3.D0*t81*s2l
     #1+2048.D0/3.D0*t84*s2l2-4096.D0/3.D0*t88+4096.D0/3.D0*t91)*s1vr+(-
     #2048*t81*s2l1-2048*t84*s2l2)*s23r+((2048.D0/3.D0*t35+4096.D0/3.D0*
     #t37)*s34r+2048.D0/3.D0*t110*s3l2+2048.D0/3.D0*t110*s4l2)*s2l1+((-2
     #048.D0/3.D0*t69-4096.D0/3.D0*t71)*s34r-2048.D0/3.D0*t110*s3l1-2048
     #.D0/3.D0*t110*s4l1)*s2l2+(4096.D0/3.D0*t24*t87-4096.D0/3.D0*t24*t9
     #0)*s34r)*s13r
      s4 = ((-t137-2048.D0/3.D0*t138)*t141+((t137+4096.D0/3.D0*t138+4096
     #.D0/3.D0*t144)*s23r+t155*s2l2+t160+2048.D0/3.D0*t163*s3l2)*s1vr+(t
     #169*s2l2-2048.D0/3.D0*t159+2048*t152*t1*s3l2+2048*t161*t175)*s23r+
     #t196*s2l2-2048.D0/3.D0*t180*t181*s3l2*t182+(4096.D0/3.D0*t3*t202*s
     #3l2-4096.D0/3.D0*t13)*s34r-2048*z*s3l2-t210)*s1l1
      s2 = s3+s4
      s1 = s2+((t214+2048.D0/3.D0*t215)*t141+((-t214-4096.D0/3.D0*t215-4
     #096.D0/3.D0*t220)*s23r-t155*s2l1-t228-2048.D0/3.D0*t163*s3l1)*s1vr
     #+(-t169*s2l1+2048.D0/3.D0*t227-2048*t152*t1*s3l1-2048*t161*t239)*s
     #23r-t196*s2l1+2048.D0/3.D0*t180*t181*s3l1*t182+(-4096.D0/3.D0*t3*t
     #202*s3l1+4096.D0/3.D0*t57)*s34r+2048*z*s3l1+t257)*s1l2+(4096.D0/3.
     #D0*t135*s3l2*s2l1-4096.D0/3.D0*t135*s3l1*s2l2)*t141+(((-8192.D0/3.
     #D0*t138-2048.D0/3.D0*t144)*s2l1+(8192.D0/3.D0*t215+2048.D0/3.D0*t2
     #20)*s2l2+4096.D0/3.D0*t135*t275-4096.D0/3.D0*t135*t278)*s23r+(-t16
     #0-2048.D0/3.D0*t287*s3l2+2048*t161*s4l2)*s2l1+(t228+2048.D0/3.D0*t
     #287*s3l1-2048*t161*s4l1)*s2l2)*s1vr
      taux11 = s1+((2048.D0/3.D0*t304*s3l2+2048.D0/3.D0*t307*s4l2)*s2l1+
     #(-2048.D0/3.D0*t304*s3l1-2048.D0/3.D0*t307*s4l1)*s2l2+4096.D0/3.D0
     #*t3*t275-4096.D0/3.D0*t3*t278)*s23r+((-2048.D0/3.D0*t3*t323*s3l2-4
     #096.D0/3.D0*t13)*s34r-2048.D0/3.D0*t332*s3l2+t210)*s2l1+((2048.D0/
     #3.D0*t3*t323*s3l1+4096.D0/3.D0*t57)*s34r+2048.D0/3.D0*t332*s3l1-t2
     #57)*s2l2+(-4096.D0/3.D0*t88+4096.D0/3.D0*t91)*s34r-4096.D0/3.D0*t1
     #*s4l2*s3l1+4096.D0/3.D0*t1*s4l1*s3l2
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t11 = t4*t2
      t12 = l3*s4l2
      t13 = t11*t12
      t22 = t4**2
      t23 = t2*t1
      t24 = t22*t23
      t28 = 1+z
      t32 = 8192.D0/3.D0*t24*l3*s34r-2048*t1*t28*t4
      t35 = t24*l3*s3l2
      t37 = t24*t12
      t41 = -t1
      t43 = 3*t41*t4
      t51 = t3*t4*s2l1
      t53 = t4*s3l1
      t54 = t3*t53
      t56 = l3*s4l1
      t57 = t11*t56
      t69 = t24*l3*s3l1
      t71 = t24*t56
      t81 = t9+t13
      t84 = -t54-t57
      t87 = t12*s3l1
      t88 = t11*t87
      t90 = t56*s3l2
      t91 = t11*t90
      t108 = 2*z
      t110 = (t41*z+2-t108)*t4
      t134 = l3**2
      t135 = t134*t2
      t137 = 2048*t135*s2l2
      t138 = t135*s3l2
      t141 = s1vr**2
      t144 = t135*s4l2
      t150 = t134*t4*t23*s34r
      t152 = z*l3
      t153 = t152*t1
      t154 = 4096*t153
      t155 = -8192.D0/3.D0*t150+t154
      t157 = t134*t23
      t159 = t157*t8*s34r
      t160 = 4096.D0/3.D0*t159
      t161 = l3*t1
      t162 = 2*t153
      t163 = -t161-t162
      t169 = 2048.D0/3.D0*t150-t154
      t175 = z*s4l2
      t180 = t134*t22
      t181 = t2**2
      t182 = s34r**2
      t187 = t108-1
      t196 = 2048.D0/3.D0*t180*t181*t182+2048.D0/3.D0*(l3*t23-3*t3*t187)
     #*t4*s34r+2048*z*t28
      t202 = z*t4
      t210 = 2048*t175
      t214 = 2048*t135*s2l1
      t215 = t135*s3l1
      t220 = t135*s4l1
      t227 = t157*t53*s34r
      t228 = 4096.D0/3.D0*t227
      t239 = z*s4l1
      t257 = 2048*t239
      t275 = s4l2*s3l1
      t278 = s3l2*s4l1
      t287 = -t161*(-1+3*z)-4*t161
      t304 = -l3*t187*t1-t162
      t307 = -t153-2*t161
      t323 = t187*t4
      t330 = z**2
      t332 = -t330-4*z+2
      s3 = (((-2048*t6-4096.D0/3.D0*t9+2048.D0/3.D0*t13)*s1vr+(4096*t6+2
     #048*t9+2048*t13)*s23r+t32*s2l2+(4096.D0/3.D0*t35-2048.D0/3.D0*t37)
     #*s34r+2048.D0/3.D0*t43*s3l2+2048.D0/3.D0*t43*s4l2)*s1l1+((2048*t51
     #+4096.D0/3.D0*t54-2048.D0/3.D0*t57)*s1vr+(-4096*t51-2048*t54-2048*
     #t57)*s23r-t32*s2l1+(-4096.D0/3.D0*t69+2048.D0/3.D0*t71)*s34r-2048.
     #D0/3.D0*t43*s3l1-2048.D0/3.D0*t43*s4l1)*s1l2+(2048.D0/3.D0*t81*s2l
     #1+2048.D0/3.D0*t84*s2l2+4096.D0/3.D0*t88-4096.D0/3.D0*t91)*s1vr+(-
     #2048*t81*s2l1-2048*t84*s2l2)*s23r+((-2048.D0/3.D0*t35-4096.D0/3.D0
     #*t37)*s34r-2048.D0/3.D0*t110*s3l2-2048.D0/3.D0*t110*s4l2)*s2l1+((2
     #048.D0/3.D0*t69+4096.D0/3.D0*t71)*s34r+2048.D0/3.D0*t110*s3l1+2048
     #.D0/3.D0*t110*s4l1)*s2l2+(-4096.D0/3.D0*t24*t87+4096.D0/3.D0*t24*t
     #90)*s34r)*s13r
      s4 = ((t137+2048.D0/3.D0*t138)*t141+((-t137-4096.D0/3.D0*t138-4096
     #.D0/3.D0*t144)*s23r+t155*s2l2-t160-2048.D0/3.D0*t163*s3l2)*s1vr+(t
     #169*s2l2+2048.D0/3.D0*t159-2048*t152*t1*s3l2-2048*t161*t175)*s23r+
     #t196*s2l2+2048.D0/3.D0*t180*t181*s3l2*t182+(-4096.D0/3.D0*t3*t202*
     #s3l2+4096.D0/3.D0*t13)*s34r+2048*z*s3l2+t210)*s1l1
      s2 = s3+s4
      s1 = s2+((-t214-2048.D0/3.D0*t215)*t141+((t214+4096.D0/3.D0*t215+4
     #096.D0/3.D0*t220)*s23r-t155*s2l1+t228+2048.D0/3.D0*t163*s3l1)*s1vr
     #+(-t169*s2l1-2048.D0/3.D0*t227+2048*t152*t1*s3l1+2048*t161*t239)*s
     #23r-t196*s2l1-2048.D0/3.D0*t180*t181*s3l1*t182+(4096.D0/3.D0*t3*t2
     #02*s3l1-4096.D0/3.D0*t57)*s34r-2048*z*s3l1-t257)*s1l2+(-4096.D0/3.
     #D0*t135*s3l2*s2l1+4096.D0/3.D0*t135*s3l1*s2l2)*t141+(((8192.D0/3.D
     #0*t138+2048.D0/3.D0*t144)*s2l1+(-8192.D0/3.D0*t215-2048.D0/3.D0*t2
     #20)*s2l2-4096.D0/3.D0*t135*t275+4096.D0/3.D0*t135*t278)*s23r+(t160
     #+2048.D0/3.D0*t287*s3l2-2048*t161*s4l2)*s2l1+(-t228-2048.D0/3.D0*t
     #287*s3l1+2048*t161*s4l1)*s2l2)*s1vr
      taux12 = s1+((-2048.D0/3.D0*t304*s3l2-2048.D0/3.D0*t307*s4l2)*s2l1
     #+(2048.D0/3.D0*t304*s3l1+2048.D0/3.D0*t307*s4l1)*s2l2-4096.D0/3.D0
     #*t3*t275+4096.D0/3.D0*t3*t278)*s23r+((2048.D0/3.D0*t3*t323*s3l2+40
     #96.D0/3.D0*t13)*s34r+2048.D0/3.D0*t332*s3l2-t210)*s2l1+((-2048.D0/
     #3.D0*t3*t323*s3l1-4096.D0/3.D0*t57)*s34r-2048.D0/3.D0*t332*s3l1+t2
     #57)*s2l2+(4096.D0/3.D0*t88-4096.D0/3.D0*t91)*s34r+4096.D0/3.D0*t1*
     #s4l2*s3l1-4096.D0/3.D0*t1*s4l1*s3l2

CCC Coefficients

      tcoff33 = taux9
      tcoff34 = taux9
      tcoff35 = taux10
      tcoff36 = taux10
      tcoff37 = taux11
      tcoff38 = taux11
      tcoff39 = taux12
      tcoff40 = taux12

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet33 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet34 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff34)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet35 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet36 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet38 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet39 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet40 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff40)

CCC Sum

      sq3num5e1 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet3
     #6*tcoff36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet40*
     #tcoff40
      RETURN
      END

      double precision function sq3num5e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3*t3
      t22 = 1024.D0/3.D0*t19*t5*s2l2
      t23 = t5*s3l2
      t24 = t19*t23
      t25 = 512.D0/3.D0*t24
      t28 = t5*t3
      t29 = l3*s4l2
      t30 = t28*t29
      t31 = 512.D0/3.D0*t30
      t34 = l3*s34r
      t38 = l3*s3l2
      t46 = 1024.D0/3.D0*t19*t5*s2l1
      t47 = t5*s3l1
      t48 = t19*t47
      t49 = 512.D0/3.D0*t48
      t52 = l3*s4l1
      t53 = t28*t52
      t66 = -t24-t30
      t67 = 512.D0/3.D0*t66
      t70 = -512.D0/3.D0*t48-512.D0/3.D0*t53
      t72 = l3**2
      t73 = z*t72
      t76 = t73*t10*t6*s34r
      t81 = t52*s3l2
      t82 = t28*t81
      t83 = 512.D0/3.D0*t82
      t108 = t7*l3
      t110 = s3l2*s34r*s2l1
      t114 = s3l1*s34r*s2l2
      t117 = t10*t2
      t119 = s34r**2
      t137 = t72*t5
      t138 = t4*s34r
      t141 = t72*t4
      t143 = t141*t23*s34r
      t145 = 512.D0/3.D0*t137*t138*s2l2+512.D0/3.D0*t143
      t149 = t72*t6
      t150 = t10*t119
      t163 = t141*t47*s34r
      t165 = 512.D0/3.D0*t137*t138*s2l1+512.D0/3.D0*t163
      t178 = t72*l3
      t180 = s1vr**2
      t185 = z*t178
      t189 = t72*t3
      t190 = s3l2*s2l1
      t192 = 512.D0/3.D0*t189*t190
      t193 = s3l1*s2l2
      t195 = 512.D0/3.D0*t189*t193
      t199 = t178*t5*t10*z*s34r
      t201 = s3l2*s3l1
      t204 = 2*z
      t211 = s23r**2
      t212 = t4*t211
      t222 = s3l2*s4l1
      t226 = 6*z-4
      t232 = t19*s3l2
      t257 = z**2
      t258 = t257-z
      t283 = 2*t258
      t289 = l3*t4
      t290 = t289*t5
      s2 = (1024.D0/3.D0*t1*t7*s1vr-1024.D0/3.D0*t1*t10*t11*s34r)*t17+((
     #(-t22-t25)*s1vr+(t22+t25+t31)*s23r+512.D0/3.D0*t7*t34*s2l2+512.D0/
     #3.D0*t7*t38*s34r)*s1l1+((-t46-t49)*s1vr+(t46+t49+512.D0/3.D0*t53)*
     #s23r+512.D0/3.D0*t7*t34*s2l1+512.D0/3.D0*t7*l3*s3l1*s34r)*s1l2+(t6
     #7*s2l1+t70*s2l2-512*t76+(-1024.D0/3.D0*t24-t31)*s3l1-t83+1024.D0/3
     #.D0*t1*t28)*s1vr+(-t67*s2l1-t70*s2l2+2048.D0/3.D0*t76-1024.D0/3.D0
     #*t66*s3l1+1024.D0/3.D0*t82+1024.D0/3.D0*t19*t5*s4l2*s4l1-1024.D0/3
     #.D0*t1*(1+z)*t3*t5)*s23r+512.D0/3.D0*t108*t110+512.D0/3.D0*t108*t1
     #14+512*t73*t117*t11*t119+((1024.D0/3.D0*t7*t38+512.D0/3.D0*t7*t29)
     #*s3l1+512.D0/3.D0*t7*t81-1024.D0/3.D0*t1*t7)*s34r)*s13r+(t145*s1vr
     #-t145*s23r-512.D0/3.D0*t149*t150*s2l2-512.D0/3.D0*t149*t10*s3l2*t1
     #19)*s1l1
      s1 = s2+(t165*s1vr-t165*s23r-512.D0/3.D0*t149*t150*s2l1-512.D0/3.D
     #0*t149*t10*s3l1*t119)*s1l2-1024.D0/3.D0*t178*t4*z*t180*s1vr+(2048.
     #D0/3.D0*t185*t4*s23r+t192+t195+2048.D0/3.D0*t199+1024.D0/3.D0*t189
     #*t201-1024.D0/3.D0*t73*t3*(t204-1))*t180
      taux9 = s1+(-1024.D0/3.D0*t185*t212+(-t192-t195-2560.D0/3.D0*t199+
     #(-1024.D0/3.D0*t189*s3l2-512.D0/3.D0*t189*s4l2)*s3l1-512.D0/3.D0*t
     #189*t222+512.D0/3.D0*t73*t3*t226)*s23r+(-512.D0/3.D0*t143+512.D0/3
     #.D0*t232)*s2l1+(-512.D0/3.D0*t163+512.D0/3.D0*t19*s3l1)*s2l2-1024.
     #D0/3.D0*t178*t6*t117*z*t119+(-1024.D0/3.D0*t141*t23*s3l1+512.D0/3.
     #D0*t73*t4*t226*t5)*s34r+1024.D0/3.D0*t19*t201-1024.D0/3.D0*t1*t2*t
     #258)*s1vr-1024.D0/3.D0*t73*t212+(-512.D0/3.D0*t19*t190-512.D0/3.D0
     #*t19*t193-2560.D0/3.D0*t73*t4*t2*t5*s34r+(-1024.D0/3.D0*t232-512.D
     #0/3.D0*t19*s4l2)*s3l1-512.D0/3.D0*t19*t222+512.D0/3.D0*t1*t2*t283)
     #*s23r-512.D0/3.D0*t290*t110-512.D0/3.D0*t290*t114+512.D0/3.D0*t73*
     #t10*t6*(-t204+3)*t119+((-1024.D0/3.D0*t289*t23+t31)*s3l1+t83+512.D
     #0/3.D0*t1*t3*t283*t5)*s34r
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3*t3
      t22 = 1024.D0/3.D0*t19*t5*s2l2
      t23 = t5*s3l2
      t24 = t19*t23
      t25 = 512.D0/3.D0*t24
      t28 = t5*t3
      t29 = l3*s4l2
      t30 = t28*t29
      t31 = 512.D0/3.D0*t30
      t34 = l3*s34r
      t38 = l3*s3l2
      t46 = 1024.D0/3.D0*t19*t5*s2l1
      t47 = t5*s3l1
      t48 = t19*t47
      t49 = 512.D0/3.D0*t48
      t52 = l3*s4l1
      t53 = t28*t52
      t66 = -t24-t30
      t67 = 512.D0/3.D0*t66
      t70 = -512.D0/3.D0*t48-512.D0/3.D0*t53
      t72 = l3**2
      t73 = z*t72
      t76 = t73*t10*t6*s34r
      t81 = t52*s3l2
      t82 = t28*t81
      t83 = 512.D0/3.D0*t82
      t108 = t7*l3
      t110 = s3l2*s34r*s2l1
      t114 = s3l1*s34r*s2l2
      t117 = t10*t2
      t119 = s34r**2
      t137 = t72*t5
      t138 = t4*s34r
      t141 = t72*t4
      t143 = t141*t23*s34r
      t145 = 512.D0/3.D0*t137*t138*s2l2+512.D0/3.D0*t143
      t149 = t72*t6
      t150 = t10*t119
      t163 = t141*t47*s34r
      t165 = 512.D0/3.D0*t137*t138*s2l1+512.D0/3.D0*t163
      t178 = t72*l3
      t180 = s1vr**2
      t185 = z*t178
      t189 = t72*t3
      t190 = s3l2*s2l1
      t192 = 512.D0/3.D0*t189*t190
      t193 = s3l1*s2l2
      t195 = 512.D0/3.D0*t189*t193
      t199 = t178*t5*t10*z*s34r
      t201 = s3l2*s3l1
      t204 = 2*z
      t211 = s23r**2
      t212 = t4*t211
      t222 = s3l2*s4l1
      t226 = 6*z-4
      t232 = t19*s3l2
      t257 = z**2
      t258 = t257-z
      t283 = 2*t258
      t289 = l3*t4
      t290 = t289*t5
      s2 = (1024.D0/3.D0*t1*t7*s1vr-1024.D0/3.D0*t1*t10*t11*s34r)*t17+((
     #(-t22-t25)*s1vr+(t22+t25+t31)*s23r+512.D0/3.D0*t7*t34*s2l2+512.D0/
     #3.D0*t7*t38*s34r)*s1l1+((-t46-t49)*s1vr+(t46+t49+512.D0/3.D0*t53)*
     #s23r+512.D0/3.D0*t7*t34*s2l1+512.D0/3.D0*t7*l3*s3l1*s34r)*s1l2+(t6
     #7*s2l1+t70*s2l2-512*t76+(-1024.D0/3.D0*t24-t31)*s3l1-t83+1024.D0/3
     #.D0*t1*t28)*s1vr+(-t67*s2l1-t70*s2l2+2048.D0/3.D0*t76-1024.D0/3.D0
     #*t66*s3l1+1024.D0/3.D0*t82+1024.D0/3.D0*t19*t5*s4l2*s4l1-1024.D0/3
     #.D0*t1*(1+z)*t3*t5)*s23r+512.D0/3.D0*t108*t110+512.D0/3.D0*t108*t1
     #14+512*t73*t117*t11*t119+((1024.D0/3.D0*t7*t38+512.D0/3.D0*t7*t29)
     #*s3l1+512.D0/3.D0*t7*t81-1024.D0/3.D0*t1*t7)*s34r)*s13r+(t145*s1vr
     #-t145*s23r-512.D0/3.D0*t149*t150*s2l2-512.D0/3.D0*t149*t10*s3l2*t1
     #19)*s1l1
      s1 = s2+(t165*s1vr-t165*s23r-512.D0/3.D0*t149*t150*s2l1-512.D0/3.D
     #0*t149*t10*s3l1*t119)*s1l2-1024.D0/3.D0*t178*t4*z*t180*s1vr+(2048.
     #D0/3.D0*t185*t4*s23r+t192+t195+2048.D0/3.D0*t199+1024.D0/3.D0*t189
     #*t201-1024.D0/3.D0*t73*t3*(t204-1))*t180
      taux10 = s1+(-1024.D0/3.D0*t185*t212+(-t192-t195-2560.D0/3.D0*t199
     #+(-1024.D0/3.D0*t189*s3l2-512.D0/3.D0*t189*s4l2)*s3l1-512.D0/3.D0*
     #t189*t222+512.D0/3.D0*t73*t3*t226)*s23r+(-512.D0/3.D0*t143+512.D0/
     #3.D0*t232)*s2l1+(-512.D0/3.D0*t163+512.D0/3.D0*t19*s3l1)*s2l2-1024
     #.D0/3.D0*t178*t6*t117*z*t119+(-1024.D0/3.D0*t141*t23*s3l1+512.D0/3
     #.D0*t73*t4*t226*t5)*s34r+1024.D0/3.D0*t19*t201-1024.D0/3.D0*t1*t2*
     #t258)*s1vr-1024.D0/3.D0*t73*t212+(-512.D0/3.D0*t19*t190-512.D0/3.D
     #0*t19*t193-2560.D0/3.D0*t73*t4*t2*t5*s34r+(-1024.D0/3.D0*t232-512.
     #D0/3.D0*t19*s4l2)*s3l1-512.D0/3.D0*t19*t222+512.D0/3.D0*t1*t2*t283
     #)*s23r-512.D0/3.D0*t290*t110-512.D0/3.D0*t290*t114+512.D0/3.D0*t73
     #*t10*t6*(-t204+3)*t119+((-1024.D0/3.D0*t289*t23+t31)*s3l1+t83+512.
     #D0/3.D0*t1*t3*t283*t5)*s34r
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 2048.D0/3.D0*t9
      t11 = t4*t2
      t12 = l3*s4l2
      t13 = t11*t12
      t21 = t4**2
      t22 = t2*t1
      t23 = t21*t22
      t27 = 1+z
      t31 = -2048.D0/3.D0*t23*l3*s34r-4096.D0/3.D0*t1*t27*t4
      t40 = t4*t1
      t48 = t3*t4*s2l1
      t50 = t4*s3l1
      t51 = t3*t50
      t52 = 2048.D0/3.D0*t51
      t53 = l3*s4l1
      t54 = t11*t53
      t67 = l3*t22
      t68 = t21*s4l1
      t79 = t9+t13
      t82 = -t51-t54
      t85 = t12*s3l1
      t88 = t53*s3l2
      t104 = t1*z
      t107 = t4*s4l2
      t117 = t4*s4l1
      t129 = l3**2
      t130 = t129*t2
      t131 = t130*s2l2
      t133 = t130*s3l2
      t135 = t130*s4l2
      t136 = 2048.D0/3.D0*t135
      t138 = s1vr**2
      t145 = t129*t4*t22*s34r
      t147 = 2048.D0/3.D0*t3
      t148 = z*l3
      t149 = t148*t1
      t150 = 8192.D0/3.D0*t149
      t151 = 2048.D0/3.D0*t145-t147+t150
      t153 = t129*t22
      t154 = t153*t8
      t156 = t153*t107
      t157 = 4096.D0/3.D0*t156
      t160 = l3*t1
      t161 = 2*t160
      t162 = -t149+t161
      t170 = -4096.D0/3.D0*t145+t147-t150
      t176 = -t3+2*t149
      t183 = t129*t21
      t184 = t2**2
      t185 = s34r**2
      t190 = 2*z
      t191 = t190-1
      t193 = 2*t3*t191
      t200 = -2048.D0/3.D0*t183*t184*t185-2048.D0/3.D0*(-3*t67+t193)*t4*
     #s34r+4096.D0/3.D0*z*t27
      t204 = t129*t184
      t221 = t130*s2l1
      t223 = t130*s3l1
      t225 = t130*s4l1
      t226 = 2048.D0/3.D0*t225
      t234 = t153*t50
      t236 = t153*t117
      t237 = 4096.D0/3.D0*t236
      t276 = 2048*t133
      t279 = 2048*t223
      t290 = s4l2*s3l1
      t293 = s3l2*s4l1
      t306 = -t3+2*t160*(-1+3*z)
      t310 = l3*t191*t1
      t311 = -t310-t161
      t329 = -t3+2*t310
      t353 = (-2*t3*t1+t193)*t4
      t360 = z**2
      t363 = -t104-t190
      s3 = (((-2048.D0/3.D0*t6+t10+4096.D0/3.D0*t13)*s1vr+(4096.D0/3.D0*
     #t6+t10+2048.D0/3.D0*t13)*s23r+t31*s2l2+(-2048.D0/3.D0*t23*l3*s3l2-
     #4096.D0/3.D0*t23*t12)*s34r-4096.D0/3.D0*t40*s3l2-4096.D0/3.D0*t40*
     #s4l2)*s1l1+((2048.D0/3.D0*t48-t52-4096.D0/3.D0*t54)*s1vr+(-4096.D0
     #/3.D0*t48-t52-2048.D0/3.D0*t54)*s23r-t31*s2l1+(2048.D0/3.D0*t23*l3
     #*s3l1+4096.D0/3.D0*t67*t68)*s34r+4096.D0/3.D0*t40*s3l1+4096.D0/3.D
     #0*t40*s4l1)*s1l2+(4096.D0/3.D0*t79*s2l1+4096.D0/3.D0*t82*s2l2+2048
     #.D0/3.D0*t11*t85-2048.D0/3.D0*t11*t88)*s1vr+(-2048.D0/3.D0*t79*s2l
     #1-2048.D0/3.D0*t82*s2l2)*s23r+(-2048.D0/3.D0*t23*t12*s34r+4096.D0/
     #3.D0*t104*t8+4096.D0/3.D0*t104*t107)*s2l1+(2048.D0/3.D0*t67*t68*s3
     #4r-4096.D0/3.D0*t104*t50-4096.D0/3.D0*t104*t117)*s2l2+(-2048.D0/3.
     #D0*t23*t85+2048.D0/3.D0*t23*t88)*s34r)*s13r
      s4 = ((2048.D0/3.D0*t131-4096.D0/3.D0*t133-t136)*t138+((-2048.D0/3
     #.D0*t131+2048.D0/3.D0*t133+2048.D0/3.D0*t135)*s23r+t151*s2l2+(2048
     #*t154+t157)*s34r+2048.D0/3.D0*t162*s3l2-2048.D0/3.D0*t3*s4l2)*s1vr
     #+(t170*s2l2+(-2048.D0/3.D0*t154-t157)*s34r-2048.D0/3.D0*t176*s3l2-
     #2048.D0/3.D0*t176*s4l2)*s23r+t200*s2l2+(-2048.D0/3.D0*t183*t184*s3
     #l2-2048.D0/3.D0*t204*t21*s4l2)*t185+(2048.D0/3.D0*t67*t8+2048.D0/3
     #.D0*t67*t107)*s34r+4096.D0/3.D0*z*s3l2+4096.D0/3.D0*z*s4l2)*s1l1
      s2 = s3+s4
      s1 = s2+((-2048.D0/3.D0*t221+4096.D0/3.D0*t223+t226)*t138+((2048.D
     #0/3.D0*t221-2048.D0/3.D0*t223-2048.D0/3.D0*t225)*s23r-t151*s2l1+(-
     #2048*t234-t237)*s34r-2048.D0/3.D0*t162*s3l1+2048.D0/3.D0*t3*s4l1)*
     #s1vr+(-t170*s2l1+(2048.D0/3.D0*t234+t237)*s34r+2048.D0/3.D0*t176*s
     #3l1+2048.D0/3.D0*t176*s4l1)*s23r-t200*s2l1+(2048.D0/3.D0*t183*t184
     #*s3l1+2048.D0/3.D0*t204*t68)*t185+(-2048.D0/3.D0*t67*t50-2048.D0/3
     #.D0*t67*t117)*s34r-4096.D0/3.D0*z*s3l1-4096.D0/3.D0*z*s4l1)*s1l2+(
     #(-t276-t136)*s2l1+(t279+t226)*s2l2)*t138
      taux11 = s1+(((t276+4096.D0/3.D0*t135)*s2l1+(-t279-4096.D0/3.D0*t2
     #25)*s2l2-2048.D0/3.D0*t130*t290+2048.D0/3.D0*t130*t293)*s23r+((409
     #6.D0/3.D0*t154+2048.D0/3.D0*t156)*s34r-2048.D0/3.D0*t306*s3l2+2048
     #.D0/3.D0*t311*s4l2)*s2l1+((-4096.D0/3.D0*t234-2048.D0/3.D0*t236)*s
     #34r+2048.D0/3.D0*t306*s3l1-2048.D0/3.D0*t311*s4l1)*s2l2)*s1vr+((20
     #48.D0/3.D0*t329*s3l2+4096.D0/3.D0*t148*t1*s4l2)*s2l1+(-2048.D0/3.D
     #0*t329*s3l1-4096.D0/3.D0*t148*t1*s4l1)*s2l2-2048.D0/3.D0*t3*t290+2
     #048.D0/3.D0*t3*t293)*s23r+((2048.D0/3.D0*t353*s3l2+2048.D0/3.D0*t1
     #48*t11*s4l2)*s34r-4096.D0/3.D0*t360*s3l2+2048.D0/3.D0*t363*s4l2)*s
     #2l1+((-2048.D0/3.D0*t353*s3l1-2048.D0/3.D0*t148*t11*s4l1)*s34r+409
     #6.D0/3.D0*t360*s3l1-2048.D0/3.D0*t363*s4l1)*s2l2
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 2048.D0/3.D0*t9
      t11 = t4*t2
      t12 = l3*s4l2
      t13 = t11*t12
      t21 = t4**2
      t22 = t2*t1
      t23 = t21*t22
      t27 = 1+z
      t31 = 2048.D0/3.D0*t23*l3*s34r+4096.D0/3.D0*t1*t27*t4
      t40 = t4*t1
      t48 = t3*t4*s2l1
      t50 = t4*s3l1
      t51 = t3*t50
      t52 = 2048.D0/3.D0*t51
      t53 = l3*s4l1
      t54 = t11*t53
      t67 = l3*t22
      t68 = t21*s4l1
      t79 = -t9-t13
      t82 = t51+t54
      t85 = t12*s3l1
      t88 = t53*s3l2
      t104 = t1*z
      t107 = t4*s4l2
      t117 = t4*s4l1
      t129 = l3**2
      t130 = t129*t2
      t131 = t130*s2l2
      t133 = t130*s3l2
      t135 = t130*s4l2
      t136 = 2048.D0/3.D0*t135
      t138 = s1vr**2
      t145 = t129*t4*t22*s34r
      t147 = 2048.D0/3.D0*t3
      t148 = z*l3
      t149 = t148*t1
      t150 = 8192.D0/3.D0*t149
      t151 = -2048.D0/3.D0*t145+t147-t150
      t153 = t129*t22
      t154 = t153*t8
      t156 = t153*t107
      t157 = 4096.D0/3.D0*t156
      t160 = l3*t1
      t161 = 2*t160
      t162 = -t149+t161
      t170 = 4096.D0/3.D0*t145-t147+t150
      t176 = -t3+2*t149
      t183 = t129*t21
      t184 = t2**2
      t185 = s34r**2
      t190 = 2*z
      t191 = t190-1
      t193 = 2*t3*t191
      t200 = 2048.D0/3.D0*t183*t184*t185+2048.D0/3.D0*(-3*t67+t193)*t4*s
     #34r-4096.D0/3.D0*z*t27
      t204 = t129*t184
      t221 = t130*s2l1
      t223 = t130*s3l1
      t225 = t130*s4l1
      t226 = 2048.D0/3.D0*t225
      t234 = t153*t50
      t236 = t153*t117
      t237 = 4096.D0/3.D0*t236
      t276 = 2048*t133
      t279 = 2048*t223
      t290 = s4l2*s3l1
      t293 = s3l2*s4l1
      t306 = -t3+2*t160*(-1+3*z)
      t310 = l3*t191*t1
      t311 = -t310-t161
      t329 = -t3+2*t310
      t353 = (-2*t3*t1+t193)*t4
      t360 = z**2
      t363 = -t104-t190
      s3 = (((2048.D0/3.D0*t6-t10-4096.D0/3.D0*t13)*s1vr+(-4096.D0/3.D0*
     #t6-t10-2048.D0/3.D0*t13)*s23r+t31*s2l2+(2048.D0/3.D0*t23*l3*s3l2+4
     #096.D0/3.D0*t23*t12)*s34r+4096.D0/3.D0*t40*s3l2+4096.D0/3.D0*t40*s
     #4l2)*s1l1+((-2048.D0/3.D0*t48+t52+4096.D0/3.D0*t54)*s1vr+(4096.D0/
     #3.D0*t48+t52+2048.D0/3.D0*t54)*s23r-t31*s2l1+(-2048.D0/3.D0*t23*l3
     #*s3l1-4096.D0/3.D0*t67*t68)*s34r-4096.D0/3.D0*t40*s3l1-4096.D0/3.D
     #0*t40*s4l1)*s1l2+(4096.D0/3.D0*t79*s2l1+4096.D0/3.D0*t82*s2l2-2048
     #.D0/3.D0*t11*t85+2048.D0/3.D0*t11*t88)*s1vr+(-2048.D0/3.D0*t79*s2l
     #1-2048.D0/3.D0*t82*s2l2)*s23r+(2048.D0/3.D0*t23*t12*s34r-4096.D0/3
     #.D0*t104*t8-4096.D0/3.D0*t104*t107)*s2l1+(-2048.D0/3.D0*t67*t68*s3
     #4r+4096.D0/3.D0*t104*t50+4096.D0/3.D0*t104*t117)*s2l2+(2048.D0/3.D
     #0*t23*t85-2048.D0/3.D0*t23*t88)*s34r)*s13r
      s4 = ((-2048.D0/3.D0*t131+4096.D0/3.D0*t133+t136)*t138+((2048.D0/3
     #.D0*t131-2048.D0/3.D0*t133-2048.D0/3.D0*t135)*s23r+t151*s2l2+(-204
     #8*t154-t157)*s34r-2048.D0/3.D0*t162*s3l2+2048.D0/3.D0*t3*s4l2)*s1v
     #r+(t170*s2l2+(2048.D0/3.D0*t154+t157)*s34r+2048.D0/3.D0*t176*s3l2+
     #2048.D0/3.D0*t176*s4l2)*s23r+t200*s2l2+(2048.D0/3.D0*t183*t184*s3l
     #2+2048.D0/3.D0*t204*t21*s4l2)*t185+(-2048.D0/3.D0*t67*t8-2048.D0/3
     #.D0*t67*t107)*s34r-4096.D0/3.D0*z*s3l2-4096.D0/3.D0*z*s4l2)*s1l1
      s2 = s3+s4
      s1 = s2+((2048.D0/3.D0*t221-4096.D0/3.D0*t223-t226)*t138+((-2048.D
     #0/3.D0*t221+2048.D0/3.D0*t223+2048.D0/3.D0*t225)*s23r-t151*s2l1+(2
     #048*t234+t237)*s34r+2048.D0/3.D0*t162*s3l1-2048.D0/3.D0*t3*s4l1)*s
     #1vr+(-t170*s2l1+(-2048.D0/3.D0*t234-t237)*s34r-2048.D0/3.D0*t176*s
     #3l1-2048.D0/3.D0*t176*s4l1)*s23r-t200*s2l1+(-2048.D0/3.D0*t183*t18
     #4*s3l1-2048.D0/3.D0*t204*t68)*t185+(2048.D0/3.D0*t67*t50+2048.D0/3
     #.D0*t67*t117)*s34r+4096.D0/3.D0*z*s3l1+4096.D0/3.D0*z*s4l1)*s1l2+(
     #(t276+t136)*s2l1+(-t279-t226)*s2l2)*t138
      taux12 = s1+(((-t276-4096.D0/3.D0*t135)*s2l1+(t279+4096.D0/3.D0*t2
     #25)*s2l2+2048.D0/3.D0*t130*t290-2048.D0/3.D0*t130*t293)*s23r+((-40
     #96.D0/3.D0*t154-2048.D0/3.D0*t156)*s34r+2048.D0/3.D0*t306*s3l2-204
     #8.D0/3.D0*t311*s4l2)*s2l1+((4096.D0/3.D0*t234+2048.D0/3.D0*t236)*s
     #34r-2048.D0/3.D0*t306*s3l1+2048.D0/3.D0*t311*s4l1)*s2l2)*s1vr+((-2
     #048.D0/3.D0*t329*s3l2-4096.D0/3.D0*t148*t1*s4l2)*s2l1+(2048.D0/3.D
     #0*t329*s3l1+4096.D0/3.D0*t148*t1*s4l1)*s2l2+2048.D0/3.D0*t3*t290-2
     #048.D0/3.D0*t3*t293)*s23r+((-2048.D0/3.D0*t353*s3l2-2048.D0/3.D0*t
     #148*t11*s4l2)*s34r+4096.D0/3.D0*t360*s3l2-2048.D0/3.D0*t363*s4l2)*
     #s2l1+((2048.D0/3.D0*t353*s3l1+2048.D0/3.D0*t148*t11*s4l1)*s34r-409
     #6.D0/3.D0*t360*s3l1+2048.D0/3.D0*t363*s4l1)*s2l2

CCC Coefficients

      tcoff33 = taux9
      tcoff34 = taux9
      tcoff35 = taux10
      tcoff36 = taux10
      tcoff37 = taux11
      tcoff38 = taux11
      tcoff39 = taux12
      tcoff40 = taux12

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet33 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet34 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff34)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet35 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet36 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet38 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet39 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet40 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff40)

CCC Sum

      sq3num5e2 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet3
     #6*tcoff36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet40*
     #tcoff40
      RETURN
      END

      double precision function sq3num5e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3**2
      t20 = t19*z
      t21 = t4*t5
      t22 = s1vr**2
      t30 = t10*t6
      t42 = s34r**2
      t52 = z*t19*l3
      t57 = s23r**2
      t58 = t4*t57
      t64 = t20*t4
      taux9 = (512.D0/3.D0*t1*t7*s1vr-512.D0/3.D0*t1*t10*t11*s34r)*t17+(
     #-512.D0/3.D0*t20*t21*t22+1024.D0/3.D0*t1*t3*t5*s1vr+(512.D0/3.D0*t
     #20*t30*s34r-512.D0/3.D0*t1*(1+z)*t3*t5)*s23r+512.D0/3.D0*t20*t10*t
     #2*t11*t42-512.D0/3.D0*t1*t7*s34r)*s13r+512.D0/3.D0*t52*t4*s23r*t22
     #+(-512.D0/3.D0*t52*t58+(-512.D0/3.D0*t52*t10*t5*s34r+512.D0/3.D0*t
     #64)*s23r-512.D0/3.D0*t20*t21*s34r)*s1vr-512.D0/3.D0*t20*t58-512.D0
     #/3.D0*t64*(z-2)*t5*s34r*s23r+512.D0/3.D0*t20*t30*t42
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3**2
      t20 = t19*z
      t21 = t4*t5
      t22 = s1vr**2
      t30 = t10*t6
      t42 = s34r**2
      t52 = z*t19*l3
      t57 = s23r**2
      t58 = t4*t57
      t64 = t20*t4
      taux10 = (512.D0/3.D0*t1*t7*s1vr-512.D0/3.D0*t1*t10*t11*s34r)*t17+
     #(-512.D0/3.D0*t20*t21*t22+1024.D0/3.D0*t1*t3*t5*s1vr+(512.D0/3.D0*
     #t20*t30*s34r-512.D0/3.D0*t1*(1+z)*t3*t5)*s23r+512.D0/3.D0*t20*t10*
     #t2*t11*t42-512.D0/3.D0*t1*t7*s34r)*s13r+512.D0/3.D0*t52*t4*s23r*t2
     #2+(-512.D0/3.D0*t52*t58+(-512.D0/3.D0*t52*t10*t5*s34r+512.D0/3.D0*
     #t64)*s23r-512.D0/3.D0*t20*t21*s34r)*s1vr-512.D0/3.D0*t20*t58-512.D
     #0/3.D0*t64*(z-2)*t5*s34r*s23r+512.D0/3.D0*t20*t30*t42
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t5 = t4*s3l2
      t12 = t3*t5
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t33 = t4*s3l1
      t40 = t3*t33
      t42 = t4*s4l1
      t43 = t3*t42
      t50 = t19*s3l1
      t58 = (-2048.D0/3.D0*t12-2048.D0/3.D0*t15)*s2l1
      t61 = (2048.D0/3.D0*t40+2048.D0/3.D0*t43)*s2l2
      t74 = t19*s4l2
      t82 = t19*s4l1
      t97 = l3**2
      t98 = t97*t2
      t99 = t98*s2l2
      t100 = t98*s4l2
      t103 = s1vr**2
      t105 = t98*s3l2
      t111 = t97*t4*t20*s34r
      t113 = 4096.D0/3.D0*t3
      t114 = -2048*t111+t113
      t116 = t97*t20
      t117 = t116*t5
      t119 = t116*t14
      t121 = 2048.D0/3.D0*t117-8192.D0/3.D0*t119
      t123 = t1*l3
      t128 = 4096.D0/3.D0*t3*s4l2
      t132 = 2048.D0/3.D0*t111-t113
      t136 = t3*s3l2
      t141 = t2**2
      t142 = s34r**2
      t149 = -2048.D0/3.D0*t97*t19*t141*t142-4096.D0/3.D0*t26*t4*s34r
      t151 = t97*t141
      t165 = t98*s2l1
      t166 = t98*s4l1
      t170 = t98*s3l1
      t176 = t116*t33
      t178 = t116*t42
      t180 = -2048.D0/3.D0*t176+8192.D0/3.D0*t178
      t186 = 4096.D0/3.D0*t3*s4l1
      t193 = t3*s3l1
      t222 = s3l2*s2l1
      t225 = s3l1*s2l2
      t228 = s4l2*s3l1
      t231 = s3l2*s4l1
      t242 = 2*z-1
      t274 = (-t26+4*t3*t1)*t4
      t277 = z*t4
      t283 = z*t1
      s3 = ((-2048.D0/3.D0*t3*t5*s1vr+(4096.D0/3.D0*t3*t4*s2l2+2048.D0/3
     #.D0*t12+2048.D0/3.D0*t15)*s23r+2048.D0/3.D0*t21*t22*s2l2+2048.D0/3
     #.D0*t26*t27*s34r)*s1l1+(2048.D0/3.D0*t3*t33*s1vr+(-4096.D0/3.D0*t3
     #*t4*s2l1-2048.D0/3.D0*t40-2048.D0/3.D0*t43)*s23r-2048.D0/3.D0*t21*
     #t22*s2l1-2048.D0/3.D0*t26*t50*s34r)*s1l2+(t58+t61-4096.D0/3.D0*t3*
     #t14*s3l1+4096.D0/3.D0*t3*t42*s3l2)*s1vr+(t58+t61)*s23r+(-2048*t26*
     #t27+4096.D0/3.D0*t26*t74)*s34r*s2l1+(2048*t26*t50-4096.D0/3.D0*t26
     #*t82)*s34r*s2l2+(4096.D0/3.D0*t26*t74*s3l1-4096.D0/3.D0*t26*t82*s3
     #l2)*s34r)*s13r
      s4 = ((4096.D0/3.D0*t99+4096.D0/3.D0*t100)*t103+((-4096.D0/3.D0*t9
     #9-4096.D0/3.D0*t105-4096.D0/3.D0*t100)*s23r+t114*s2l2+t121*s34r+40
     #96.D0/3.D0*t123*z*s3l2+t128)*s1vr+(t132*s2l2-t121*s34r-4096.D0/3.D
     #0*t136-t128)*s23r+t149*s2l2+(-2048.D0/3.D0*t151*t27+4096.D0/3.D0*t
     #151*t74)*t142+(-4096.D0/3.D0*t26*t5-4096.D0/3.D0*t26*t14)*s34r)*s1
     #l1
      s2 = s3+s4
      s1 = s2+((-4096.D0/3.D0*t165-4096.D0/3.D0*t166)*t103+((4096.D0/3.D
     #0*t165+4096.D0/3.D0*t170+4096.D0/3.D0*t166)*s23r-t114*s2l1+t180*s3
     #4r-4096.D0/3.D0*t123*z*s3l1-t186)*s1vr+(-t132*s2l1-t180*s34r+4096.
     #D0/3.D0*t193+t186)*s23r-t149*s2l1+(2048.D0/3.D0*t151*t50-4096.D0/3
     #.D0*t151*t82)*t142+(4096.D0/3.D0*t26*t33+4096.D0/3.D0*t26*t42)*s34
     #r)*s1l2+((-2048.D0/3.D0*t105+4096.D0/3.D0*t100)*s2l1+(2048.D0/3.D0
     #*t170-4096.D0/3.D0*t166)*s2l2)*t103
      taux11 = s1+((2048.D0/3.D0*t98*t222-2048.D0/3.D0*t98*t225+4096.D0/
     #3.D0*t98*t228-4096.D0/3.D0*t98*t231)*s23r+((2048*t117-4096.D0/3.D0
     #*t119)*s34r-2048.D0/3.D0*t136+4096.D0/3.D0*t123*t242*s4l2)*s2l1+((
     #-2048*t176+4096.D0/3.D0*t178)*s34r+2048.D0/3.D0*t193-4096.D0/3.D0*
     #t123*t242*s4l1)*s2l2)*s1vr+(2048.D0/3.D0*t3*t222-2048.D0/3.D0*t3*t
     #225+4096.D0/3.D0*t3*t228-4096.D0/3.D0*t3*t231)*s23r+((2048.D0/3.D0
     #*t274*s3l2-4096.D0/3.D0*t3*t277*s4l2)*s34r+4096.D0/3.D0*t283*s4l2)
     #*s2l1+((-2048.D0/3.D0*t274*s3l1+4096.D0/3.D0*t3*t277*s4l1)*s34r-40
     #96.D0/3.D0*t283*s4l1)*s2l2
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t5 = t4*s3l2
      t12 = t3*t5
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t33 = t4*s3l1
      t40 = t3*t33
      t42 = t4*s4l1
      t43 = t3*t42
      t50 = t19*s3l1
      t58 = (2048.D0/3.D0*t12+2048.D0/3.D0*t15)*s2l1
      t61 = (-2048.D0/3.D0*t40-2048.D0/3.D0*t43)*s2l2
      t74 = t19*s4l2
      t82 = t19*s4l1
      t97 = l3**2
      t98 = t97*t2
      t99 = t98*s2l2
      t100 = t98*s4l2
      t103 = s1vr**2
      t105 = t98*s3l2
      t111 = t97*t4*t20*s34r
      t113 = 4096.D0/3.D0*t3
      t114 = 2048*t111-t113
      t116 = t97*t20
      t117 = t116*t5
      t119 = t116*t14
      t121 = -2048.D0/3.D0*t117+8192.D0/3.D0*t119
      t123 = t1*l3
      t128 = 4096.D0/3.D0*t3*s4l2
      t132 = -2048.D0/3.D0*t111+t113
      t136 = t3*s3l2
      t141 = t2**2
      t142 = s34r**2
      t149 = 2048.D0/3.D0*t97*t19*t141*t142+4096.D0/3.D0*t26*t4*s34r
      t151 = t97*t141
      t165 = t98*s2l1
      t166 = t98*s4l1
      t170 = t98*s3l1
      t176 = t116*t33
      t178 = t116*t42
      t180 = 2048.D0/3.D0*t176-8192.D0/3.D0*t178
      t186 = 4096.D0/3.D0*t3*s4l1
      t193 = t3*s3l1
      t222 = s3l2*s2l1
      t225 = s3l1*s2l2
      t228 = s4l2*s3l1
      t231 = s3l2*s4l1
      t242 = 2*z-1
      t274 = (-t26+4*t3*t1)*t4
      t277 = z*t4
      t283 = z*t1
      s3 = ((2048.D0/3.D0*t3*t5*s1vr+(-4096.D0/3.D0*t3*t4*s2l2-2048.D0/3
     #.D0*t12-2048.D0/3.D0*t15)*s23r-2048.D0/3.D0*t21*t22*s2l2-2048.D0/3
     #.D0*t26*t27*s34r)*s1l1+(-2048.D0/3.D0*t3*t33*s1vr+(4096.D0/3.D0*t3
     #*t4*s2l1+2048.D0/3.D0*t40+2048.D0/3.D0*t43)*s23r+2048.D0/3.D0*t21*
     #t22*s2l1+2048.D0/3.D0*t26*t50*s34r)*s1l2+(t58+t61+4096.D0/3.D0*t3*
     #t14*s3l1-4096.D0/3.D0*t3*t42*s3l2)*s1vr+(t58+t61)*s23r+(2048*t26*t
     #27-4096.D0/3.D0*t26*t74)*s34r*s2l1+(-2048*t26*t50+4096.D0/3.D0*t26
     #*t82)*s34r*s2l2+(-4096.D0/3.D0*t26*t74*s3l1+4096.D0/3.D0*t26*t82*s
     #3l2)*s34r)*s13r
      s4 = ((-4096.D0/3.D0*t99-4096.D0/3.D0*t100)*t103+((4096.D0/3.D0*t9
     #9+4096.D0/3.D0*t105+4096.D0/3.D0*t100)*s23r+t114*s2l2+t121*s34r-40
     #96.D0/3.D0*t123*z*s3l2-t128)*s1vr+(t132*s2l2-t121*s34r+4096.D0/3.D
     #0*t136+t128)*s23r+t149*s2l2+(2048.D0/3.D0*t151*t27-4096.D0/3.D0*t1
     #51*t74)*t142+(4096.D0/3.D0*t26*t5+4096.D0/3.D0*t26*t14)*s34r)*s1l1
      s2 = s3+s4
      s1 = s2+((4096.D0/3.D0*t165+4096.D0/3.D0*t166)*t103+((-4096.D0/3.D
     #0*t165-4096.D0/3.D0*t170-4096.D0/3.D0*t166)*s23r-t114*s2l1+t180*s3
     #4r+4096.D0/3.D0*t123*z*s3l1+t186)*s1vr+(-t132*s2l1-t180*s34r-4096.
     #D0/3.D0*t193-t186)*s23r-t149*s2l1+(-2048.D0/3.D0*t151*t50+4096.D0/
     #3.D0*t151*t82)*t142+(-4096.D0/3.D0*t26*t33-4096.D0/3.D0*t26*t42)*s
     #34r)*s1l2+((2048.D0/3.D0*t105-4096.D0/3.D0*t100)*s2l1+(-2048.D0/3.
     #D0*t170+4096.D0/3.D0*t166)*s2l2)*t103
      taux12 = s1+((-2048.D0/3.D0*t98*t222+2048.D0/3.D0*t98*t225-4096.D0
     #/3.D0*t98*t228+4096.D0/3.D0*t98*t231)*s23r+((-2048*t117+4096.D0/3.
     #D0*t119)*s34r+2048.D0/3.D0*t136-4096.D0/3.D0*t123*t242*s4l2)*s2l1+
     #((2048*t176-4096.D0/3.D0*t178)*s34r-2048.D0/3.D0*t193+4096.D0/3.D0
     #*t123*t242*s4l1)*s2l2)*s1vr+(-2048.D0/3.D0*t3*t222+2048.D0/3.D0*t3
     #*t225-4096.D0/3.D0*t3*t228+4096.D0/3.D0*t3*t231)*s23r+((-2048.D0/3
     #.D0*t274*s3l2+4096.D0/3.D0*t3*t277*s4l2)*s34r-4096.D0/3.D0*t283*s4
     #l2)*s2l1+((2048.D0/3.D0*t274*s3l1-4096.D0/3.D0*t3*t277*s4l1)*s34r+
     #4096.D0/3.D0*t283*s4l1)*s2l2

CCC Coefficients

      tcoff33 = taux9
      tcoff34 = taux9
      tcoff35 = taux10
      tcoff36 = taux10
      tcoff37 = taux11
      tcoff38 = taux11
      tcoff39 = taux12
      tcoff40 = taux12

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet33 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet34 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff34)
      t1 = 1-z
      t2 = t1*l3
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet35 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet36 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet38 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet39 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet40 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff40)

CCC Sum

      sq3num5e3 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet3
     #6*tcoff36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet40*
     #tcoff40
      RETURN
      END

      double precision function sq3num5e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux9 = 0
      taux10 = 0
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t7 = 8192.D0/3.D0*t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 4096.D0/3.D0*t9
      t14 = t3*t4*s4l2
      t18 = t4**2
      t19 = t2*t1
      t20 = t18*t19
      t21 = l3*s34r
      t25 = l3*t19
      t26 = t18*s3l2
      t34 = 8192.D0/3.D0*t3*t4*s2l1
      t35 = t4*s3l1
      t36 = t3*t35
      t37 = 4096.D0/3.D0*t36
      t41 = t3*t4*s4l1
      t48 = t18*s3l1
      t55 = 4096.D0/3.D0*t9+4096.D0/3.D0*t14
      t58 = -4096.D0/3.D0*t36-4096.D0/3.D0*t41
      t68 = t25*t18
      t70 = s3l2*s34r*s2l1
      t74 = s3l1*s34r*s2l2
      t79 = l3**2
      t80 = t79*t4
      t81 = t19*s34r
      t84 = t79*t19
      t86 = t84*t8*s34r
      t88 = 4096.D0/3.D0*t80*t81*s2l2+4096.D0/3.D0*t86
      t92 = t79*t18
      t93 = t2**2
      t94 = s34r**2
      t95 = t93*t94
      t99 = t79*t93
      t108 = t84*t35*s34r
      t110 = -4096.D0/3.D0*t80*t81*s2l1-4096.D0/3.D0*t108
      t122 = t79*t2
      t123 = s3l2*s2l1
      t125 = s3l1*s2l2
      t128 = -4096.D0/3.D0*t122*t123+4096.D0/3.D0*t122*t125
      t129 = s1vr**2
      t148 = t25*t4
      taux11 = (((-t7-t10)*s1vr+(t7+t10+4096.D0/3.D0*t14)*s23r+4096.D0/3
     #.D0*t20*t21*s2l2+4096.D0/3.D0*t25*t26*s34r)*s1l1+((t34+t37)*s1vr+(
     #-t34-t37-4096.D0/3.D0*t41)*s23r-4096.D0/3.D0*t20*t21*s2l1-4096.D0/
     #3.D0*t25*t48*s34r)*s1l2+(t55*s2l1+t58*s2l2)*s1vr+(-t55*s2l1-t58*s2
     #l2)*s23r-4096.D0/3.D0*t68*t70+4096.D0/3.D0*t68*t74)*s13r+(t88*s1vr
     #-t88*s23r-4096.D0/3.D0*t92*t95*s2l2-4096.D0/3.D0*t99*t26*t94)*s1l1
     #+(t110*s1vr-t110*s23r+4096.D0/3.D0*t92*t95*s2l1+4096.D0/3.D0*t99*t
     #48*t94)*s1l2+t128*t129+(-t128*s23r+(4096.D0/3.D0*t86-4096.D0/3.D0*
     #t3*s3l2)*s2l1+(-4096.D0/3.D0*t108+4096.D0/3.D0*t3*s3l1)*s2l2)*s1vr
     #+(4096.D0/3.D0*t3*t123-4096.D0/3.D0*t3*t125)*s23r+4096.D0/3.D0*t14
     #8*t70-4096.D0/3.D0*t148*t74
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t7 = 8192.D0/3.D0*t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 4096.D0/3.D0*t9
      t14 = t3*t4*s4l2
      t18 = t4**2
      t19 = t2*t1
      t20 = t18*t19
      t21 = l3*s34r
      t25 = l3*t19
      t26 = t18*s3l2
      t34 = 8192.D0/3.D0*t3*t4*s2l1
      t35 = t4*s3l1
      t36 = t3*t35
      t37 = 4096.D0/3.D0*t36
      t41 = t3*t4*s4l1
      t48 = t18*s3l1
      t55 = -4096.D0/3.D0*t9-4096.D0/3.D0*t14
      t58 = 4096.D0/3.D0*t36+4096.D0/3.D0*t41
      t68 = t25*t18
      t70 = s3l2*s34r*s2l1
      t74 = s3l1*s34r*s2l2
      t79 = l3**2
      t80 = t79*t4
      t81 = t19*s34r
      t84 = t79*t19
      t86 = t84*t8*s34r
      t88 = -4096.D0/3.D0*t80*t81*s2l2-4096.D0/3.D0*t86
      t92 = t79*t18
      t93 = t2**2
      t94 = s34r**2
      t95 = t93*t94
      t99 = t79*t93
      t108 = t84*t35*s34r
      t110 = 4096.D0/3.D0*t80*t81*s2l1+4096.D0/3.D0*t108
      t122 = t79*t2
      t123 = s3l2*s2l1
      t125 = s3l1*s2l2
      t128 = 4096.D0/3.D0*t122*t123-4096.D0/3.D0*t122*t125
      t129 = s1vr**2
      t148 = t25*t4
      taux12 = (((t7+t10)*s1vr+(-t7-t10-4096.D0/3.D0*t14)*s23r-4096.D0/3
     #.D0*t20*t21*s2l2-4096.D0/3.D0*t25*t26*s34r)*s1l1+((-t34-t37)*s1vr+
     #(t34+t37+4096.D0/3.D0*t41)*s23r+4096.D0/3.D0*t20*t21*s2l1+4096.D0/
     #3.D0*t25*t48*s34r)*s1l2+(t55*s2l1+t58*s2l2)*s1vr+(-t55*s2l1-t58*s2
     #l2)*s23r+4096.D0/3.D0*t68*t70-4096.D0/3.D0*t68*t74)*s13r+(t88*s1vr
     #-t88*s23r+4096.D0/3.D0*t92*t95*s2l2+4096.D0/3.D0*t99*t26*t94)*s1l1
     #+(t110*s1vr-t110*s23r-4096.D0/3.D0*t92*t95*s2l1-4096.D0/3.D0*t99*t
     #48*t94)*s1l2+t128*t129+(-t128*s23r+(-4096.D0/3.D0*t86+4096.D0/3.D0
     #*t3*s3l2)*s2l1+(4096.D0/3.D0*t108-4096.D0/3.D0*t3*s3l1)*s2l2)*s1vr
     #+(-4096.D0/3.D0*t3*t123+4096.D0/3.D0*t3*t125)*s23r-4096.D0/3.D0*t1
     #48*t70+4096.D0/3.D0*t148*t74

CCC Coefficients

      tcoff33 = taux9
      tcoff34 = taux9
      tcoff35 = taux10
      tcoff36 = taux10
      tcoff37 = taux11
      tcoff38 = taux11
      tcoff39 = taux12
      tcoff40 = taux12

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet33 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet34 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff34)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet35 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet36 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet38 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet39 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet40 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff40)

CCC Sum

      sq3num5e4 = tFjet33*tcoff33+tFjet34*tcoff34+tFjet35*tcoff35+tFjet3
     #6*tcoff36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet39*tcoff39+tFjet40*
     #tcoff40
      RETURN
      END

      double precision function sq3num6e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t7 = t1**2
      t8 = t3*t2
      t9 = t7*t8
      t11 = t9*l3*s34r
      t15 = 1+z
      t19 = -64.D0/3.D0*t11-64.D0/3.D0*t2*t15*t1
      t21 = t1*t2
      t23 = 64.D0/3.D0*t21*s3l2
      t25 = 64.D0/3.D0*t21*s4l2
      t29 = t21*s3l1
      t31 = t21*s4l1
      t35 = l3*t3
      t37 = s2l2*s2l1
      t48 = z*l3
      t54 = l3**2
      t55 = t54*t3
      t57 = l3*t2
      t62 = t1*s34r
      t66 = t35*t62
      t75 = 64.D0/3.D0*t54*t1*t8*s34r+64.D0/3.D0*t57*t15
      t78 = 64.D0/3.D0*t57*s3l2
      t80 = 64.D0/3.D0*t57*s4l2
      t84 = -64.D0/3.D0*t66-64.D0/3.D0-64.D0/3.D0*z
      t86 = 64.D0/3.D0*s4l2
      t87 = 64.D0/3.D0*s3l2
      t91 = t57*s3l1
      t93 = t57*s4l1
      taux13 = (((-128.D0/3.D0*t4*l3*s1vr+128.D0/3.D0*t11)*s1l2+t19*s2l2
     #-t23-t25)*s1l1+(t19*s2l1-64.D0/3.D0*t29-64.D0/3.D0*t31)*s1l2+128.D
     #0/3.D0*t35*t1*t37*s1vr+(128.D0/3.D0*t4*s2l2-t23-t25)*s2l1+(-64.D0/
     #3.D0*t29-64.D0/3.D0*t31)*s2l2+128.D0/3.D0*t48*t9*s34r)*s13r+(((128
     #.D0/3.D0*t55*s23r-128.D0/3.D0*t57)*s1vr-128.D0/3.D0*t54*t8*t62*s23
     #r+128.D0/3.D0*t66)*s1l2+(t75*s2l2+t78+t80)*s23r+t84*s2l2-t86-t87)*
     #s1l1+((t75*s2l1+64.D0/3.D0*t91+64.D0/3.D0*t93)*s23r+t84*s2l1-64.D0
     #/3.D0*s3l1-64.D0/3.D0*s4l1)*s1l2+(-128.D0/3.D0*t55*t37*s23r+128.D0
     #/3.D0*t2*s2l2*l3*s2l1)*s1vr+((-128.D0/3.D0*t35*s2l2+t78+t80)*s2l1+
     #(64.D0/3.D0*t91+64.D0/3.D0*t93)*s2l2-128.D0/3.D0*t54*z*t8*t1*s34r)
     #*s23r+(128.D0/3.D0*t2*s2l2-t87-t86)*s2l1+(-64.D0/3.D0*s3l1-64.D0/3
     #.D0*s4l1)*s2l2+128.D0/3.D0*t48*t4*s34r
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t7 = t1**2
      t8 = t3*t2
      t9 = t7*t8
      t11 = t9*l3*s34r
      t15 = 1+z
      t19 = -64.D0/3.D0*t11-64.D0/3.D0*t2*t15*t1
      t21 = t1*t2
      t23 = 64.D0/3.D0*t21*s3l2
      t25 = 64.D0/3.D0*t21*s4l2
      t29 = t21*s3l1
      t31 = t21*s4l1
      t35 = l3*t3
      t37 = s2l2*s2l1
      t48 = z*l3
      t54 = l3**2
      t55 = t54*t3
      t57 = l3*t2
      t62 = t1*s34r
      t66 = t35*t62
      t75 = 64.D0/3.D0*t54*t1*t8*s34r+64.D0/3.D0*t57*t15
      t78 = 64.D0/3.D0*t57*s3l2
      t80 = 64.D0/3.D0*t57*s4l2
      t84 = -64.D0/3.D0*t66-64.D0/3.D0-64.D0/3.D0*z
      t86 = 64.D0/3.D0*s4l2
      t87 = 64.D0/3.D0*s3l2
      t91 = t57*s3l1
      t93 = t57*s4l1
      taux14 = (((-128.D0/3.D0*t4*l3*s1vr+128.D0/3.D0*t11)*s1l2+t19*s2l2
     #-t23-t25)*s1l1+(t19*s2l1-64.D0/3.D0*t29-64.D0/3.D0*t31)*s1l2+128.D
     #0/3.D0*t35*t1*t37*s1vr+(128.D0/3.D0*t4*s2l2-t23-t25)*s2l1+(-64.D0/
     #3.D0*t29-64.D0/3.D0*t31)*s2l2+128.D0/3.D0*t48*t9*s34r)*s13r+(((128
     #.D0/3.D0*t55*s23r-128.D0/3.D0*t57)*s1vr-128.D0/3.D0*t54*t8*t62*s23
     #r+128.D0/3.D0*t66)*s1l2+(t75*s2l2+t78+t80)*s23r+t84*s2l2-t86-t87)*
     #s1l1+((t75*s2l1+64.D0/3.D0*t91+64.D0/3.D0*t93)*s23r+t84*s2l1-64.D0
     #/3.D0*s3l1-64.D0/3.D0*s4l1)*s1l2+(-128.D0/3.D0*t55*t37*s23r+128.D0
     #/3.D0*t2*s2l2*l3*s2l1)*s1vr+((-128.D0/3.D0*t35*s2l2+t78+t80)*s2l1+
     #(64.D0/3.D0*t91+64.D0/3.D0*t93)*s2l2-128.D0/3.D0*t54*z*t8*t1*s34r)
     #*s23r+(128.D0/3.D0*t2*s2l2-t87-t86)*s2l1+(-64.D0/3.D0*s3l1-64.D0/3
     #.D0*s4l1)*s2l2+128.D0/3.D0*t48*t4*s34r
      t1 = l3-1
      t2 = t1**2
      t3 = -1+z
      t4 = t3**2
      t5 = t4*t3
      t9 = 1+z
      t13 = 256.D0/3.D0*t2*t5*l3*s34r-256.D0/3.D0*t3*t9*t1
      t15 = t1*t3
      t16 = t15*s3l2
      t18 = t15*s4l2
      t24 = t15*s3l1
      t26 = t15*s4l1
      t38 = l3**2
      t42 = l3*t3
      t46 = -256.D0/3.D0*t38*t1*t5*s34r+256.D0/3.D0*t42*t9
      t48 = t42*s3l2
      t50 = t42*s4l2
      t58 = 256.D0/3.D0*l3*t4*t1*s34r-256.D0/3.D0-256.D0/3.D0*z
      t66 = t42*s3l1
      t68 = t42*s4l1
      taux15 = ((t13*s2l2-256.D0/3.D0*t16-256.D0/3.D0*t18)*s1l1+(-t13*s2
     #l1+256.D0/3.D0*t24+256.D0/3.D0*t26)*s1l2+(256.D0/3.D0*t16+256.D0/3
     #.D0*t18)*s2l1+(-256.D0/3.D0*t24-256.D0/3.D0*t26)*s2l2)*s13r+((t46*
     #s2l2+256.D0/3.D0*t48+256.D0/3.D0*t50)*s23r+t58*s2l2-256.D0/3.D0*s3
     #l2-256.D0/3.D0*s4l2)*s1l1+((-t46*s2l1-256.D0/3.D0*t66-256.D0/3.D0*
     #t68)*s23r-t58*s2l1+256.D0/3.D0*s3l1+256.D0/3.D0*s4l1)*s1l2+((-256.
     #D0/3.D0*t48-256.D0/3.D0*t50)*s2l1+(256.D0/3.D0*t66+256.D0/3.D0*t68
     #)*s2l2)*s23r+(256.D0/3.D0*s3l2+256.D0/3.D0*s4l2)*s2l1+(-256.D0/3.D
     #0*s3l1-256.D0/3.D0*s4l1)*s2l2
      t1 = l3-1
      t2 = t1**2
      t3 = -1+z
      t4 = t3**2
      t5 = t4*t3
      t9 = 1+z
      t13 = -256.D0/3.D0*t2*t5*l3*s34r+256.D0/3.D0*t3*t9*t1
      t15 = t1*t3
      t16 = t15*s3l2
      t18 = t15*s4l2
      t24 = t15*s3l1
      t26 = t15*s4l1
      t38 = l3**2
      t42 = l3*t3
      t46 = 256.D0/3.D0*t38*t1*t5*s34r-256.D0/3.D0*t42*t9
      t48 = t42*s3l2
      t50 = t42*s4l2
      t58 = -256.D0/3.D0*l3*t4*t1*s34r+256.D0/3.D0+256.D0/3.D0*z
      t66 = t42*s3l1
      t68 = t42*s4l1
      taux16 = ((t13*s2l2+256.D0/3.D0*t16+256.D0/3.D0*t18)*s1l1+(-t13*s2
     #l1-256.D0/3.D0*t24-256.D0/3.D0*t26)*s1l2+(-256.D0/3.D0*t16-256.D0/
     #3.D0*t18)*s2l1+(256.D0/3.D0*t24+256.D0/3.D0*t26)*s2l2)*s13r+((t46*
     #s2l2-256.D0/3.D0*t48-256.D0/3.D0*t50)*s23r+t58*s2l2+256.D0/3.D0*s3
     #l2+256.D0/3.D0*s4l2)*s1l1+((-t46*s2l1+256.D0/3.D0*t66+256.D0/3.D0*
     #t68)*s23r-t58*s2l1-256.D0/3.D0*s3l1-256.D0/3.D0*s4l1)*s1l2+((256.D
     #0/3.D0*t48+256.D0/3.D0*t50)*s2l1+(-256.D0/3.D0*t66-256.D0/3.D0*t68
     #)*s2l2)*s23r+(-256.D0/3.D0*s3l2-256.D0/3.D0*s4l2)*s2l1+(256.D0/3.D
     #0*s3l1+256.D0/3.D0*s4l1)*s2l2

CCC Coefficients

      tcoff41 = taux13
      tcoff42 = taux13
      tcoff43 = taux14
      tcoff44 = taux14
      tcoff45 = taux15
      tcoff46 = taux15
      tcoff47 = taux16
      tcoff48 = taux16

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet42 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet43 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet44 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet45 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet46 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet47 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet48 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff48)

CCC Sum

      sq3num6e0 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet4
     #4*tcoff44+tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet48*
     #tcoff48
      RETURN
      END

      double precision function sq3num6e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t11 = t6*t3
      t25 = t3**2
      t26 = t6*t5
      t32 = s13r**2
      t34 = t3*t5
      t36 = t34*l3*s2l2
      t38 = l3*t3
      t39 = t5*s3l2
      t40 = t38*t39
      t41 = 128.D0/3.D0*t40
      t42 = l3*s4l2
      t43 = t34*t42
      t44 = 64.D0/3.D0*t43
      t49 = -128.D0/3.D0*t36-64.D0/3.D0*t40-t44
      t51 = l3*s34r
      t55 = l3*t6
      t58 = 128.D0/3.D0*t55*t4*s3l2
      t60 = 64.D0/3.D0*t7*t42
      t66 = t34*l3*s2l1
      t68 = t5*s3l1
      t69 = t38*t68
      t71 = l3*s4l1
      t72 = t34*t71
      t73 = 64.D0/3.D0*t72
      t95 = 64.D0/3.D0*t69+64.D0/3.D0*t72
      t97 = l3**2
      t98 = z*t97
      t99 = t25*t6
      t105 = t71*s3l2
      t106 = t34*t105
      t112 = -t40-t43
      t119 = 128.D0/3.D0*t106
      t123 = 1+z
      t137 = t5*t2
      t138 = z*s3l2
      t157 = s34r**2
      t169 = s3l2*s3l1
      t172 = s4l1*s4l2
      t175 = z**2
      s1 = ((64.D0/3.D0*t36+t41-t44)*s1vr+t49*s23r-64.D0/3.D0*t7*t51*s2l
     #2+(-t58+t60)*s34r)*s1l1+((64.D0/3.D0*t66+128.D0/3.D0*t69-t73)*s1vr
     #+(-128.D0/3.D0*t66-64.D0/3.D0*t69-t73)*s23r-64.D0/3.D0*t7*t51*s2l1
     #+(-128.D0/3.D0*t55*t4*s3l1+64.D0/3.D0*t7*t71)*s34r)*s1l2+(-t49*s2l
     #1+t95*s2l2+64.D0/3.D0*t98*t99*s34r+(t41+t44)*s3l1+64.D0/3.D0*t106-
     #128.D0/3.D0*t1*t34)*s1vr+(64.D0/3.D0*t112*s2l1-t95*s2l2+128.D0/3.D
     #0*t112*s3l1-t119-128.D0/3.D0*t34*t42*s4l1+128.D0/3.D0*t1*t123*t3*t
     #5)*s23r+((128.D0/3.D0*t7*t51+128.D0/3.D0*t34)*s2l2+64.D0/3.D0*t7*t
     #42*s34r+128.D0/3.D0*t137*t138-128.D0/3.D0*t137*s4l2)*s2l1
      t179 = s1+(64.D0/3.D0*t7*t71*s34r+128.D0/3.D0*t137*z*s3l1-128.D0/3
     #.D0*t137*s4l1)*s2l2-64.D0/3.D0*t98*t25*t2*t26*t157+((-t58-t60)*s3l
     #1-64.D0/3.D0*t7*t105+128.D0/3.D0*t1*t7)*s34r+128.D0/3.D0*t137*t169
     #-128.D0/3.D0*t137*t172+128.D0/3.D0*t175*t2*t5
      t181 = t97*t3
      t185 = s23r**2
      t188 = t97*t4
      t196 = t181*s3l2
      t197 = t181*s4l2
      t202 = t188*t39*s34r
      t204 = t2*s4l2
      t205 = t204*l3
      t216 = t5*t4
      t222 = l3*t2
      t231 = t97*t6
      t232 = t25*t157
      t236 = t97*t25
      t247 = t181*s3l1
      t248 = t181*s4l1
      t253 = t188*t68*s34r
      t256 = t2*s4l1*l3
      t292 = s1vr**2
      t296 = t97*l3
      t298 = z*t185
      t301 = 128.D0/3.D0*t196
      t302 = 64.D0/3.D0*t197
      t311 = t25*t5*s34r
      t316 = s3l2*s4l1
      t339 = 2*t2
      t342 = t38*s4l2
      t354 = z-2
      t355 = t354*s3l2
      t372 = t38*t5
      t376 = 64.D0/3.D0*t2
      s3 = (-128.D0/3.D0*t1*t7*s1vr+(-256.D0/3.D0*t11*s2l2-128.D0/3.D0*t
     #11*s3l2-128.D0/3.D0*t11*s4l2)*s2l1+(-128.D0/3.D0*t11*s3l1-128.D0/3
     #.D0*t11*s4l1)*s2l2+128.D0/3.D0*t1*t25*t26*s34r)*t32+t179*s13r
      s2 = s3+((128.D0/3.D0*t181*s23r*s1vr-256.D0/3.D0*t181*t185-256.D0/
     #3.D0*t188*t5*s34r*s23r)*s1l2+((64.D0/3.D0*t181*s2l2+64.D0/3.D0*t19
     #6+64.D0/3.D0*t197)*s23r-64.D0/3.D0*t202-64.D0/3.D0*t205)*s1vr+(-12
     #8.D0/3.D0*t196-128.D0/3.D0*t197)*t185+(64.D0/3.D0*t38*s2l2+(-64.D0
     #/3.D0*t188*t39-128.D0/3.D0*t216*t97*s4l2)*s34r-64.D0/3.D0*t222*t12
     #3*s3l2+64.D0/3.D0*t222*t123*s4l2)*s23r+64.D0/3.D0*t231*t232*s2l2+6
     #4.D0/3.D0*t236*t6*s3l2*t157+(-64.D0/3.D0*t40+64.D0/3.D0*t43)*s34r)
     #*s1l1+(((64.D0/3.D0*t181*s2l1+64.D0/3.D0*t247+64.D0/3.D0*t248)*s23
     #r-64.D0/3.D0*t253-64.D0/3.D0*t256)*s1vr+(-128.D0/3.D0*t247-128.D0/
     #3.D0*t248)*t185+(64.D0/3.D0*t38*s2l1+(-64.D0/3.D0*t188*t68-128.D0/
     #3.D0*t216*t97*s4l1)*s34r-64.D0/3.D0*t222*t123*s3l1+64.D0/3.D0*t222
     #*t123*s4l1)*s23r+64.D0/3.D0*t231*t232*s2l1+64.D0/3.D0*t236*t6*s3l1
     #*t157+(-64.D0/3.D0*t69+64.D0/3.D0*t72)*s34r)*s1l2
      s1 = s2-128.D0/3.D0*t181*t169*t292+(128.D0/3.D0*t296*t4*t298+((t30
     #1-t302)*s2l1+(128.D0/3.D0*t247-64.D0/3.D0*t248)*s2l2+64.D0/3.D0*z*
     #t296*t311+(t301+t302)*s3l1+64.D0/3.D0*t181*t316-128.D0/3.D0*t181*z
     #)*s23r+(64.D0/3.D0*t202+64.D0/3.D0*t205)*s2l1+(64.D0/3.D0*t253+64.
     #D0/3.D0*t256)*s2l2+128.D0/3.D0*t188*t5*t169*s34r-128.D0/3.D0*t38*t
     #169)*s1vr+128.D0/3.D0*t188*t298+((64.D0/3.D0*t222*t339*s3l2-64.D0/
     #3.D0*t342)*s2l1+(64.D0/3.D0*t222*t339*s3l1-64.D0/3.D0*t38*s4l1)*s2
     #l2+64.D0/3.D0*t98*t311+(128.D0/3.D0*t222*t355+64.D0/3.D0*t342)*s3l
     #1+64.D0/3.D0*t38*t316+128.D0/3.D0*t222*t172-64.D0/3.D0*t1*t2*(4*z-
     #2))*s23r
      taux13 = s1+(64.D0/3.D0*t372*t355*s34r+t376*s4l2)*s2l1+(64.D0/3.D0
     #*t372*t354*s3l1*s34r+t376*s4l1)*s2l2-64.D0/3.D0*t98*t99*t157+((128
     #.D0/3.D0*t38*t5*t354*s3l2-128.D0/3.D0*t43)*s3l1-t119)*s34r+(-128.D
     #0/3.D0*t138+64.D0/3.D0*t204)*s3l1+64.D0/3.D0*t2*s3l2*s4l1-128.D0/3
     #.D0*t172+128.D0/3.D0*t175
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t11 = t6*t3
      t25 = t3**2
      t26 = t6*t5
      t32 = s13r**2
      t34 = t3*t5
      t36 = t34*l3*s2l2
      t38 = l3*t3
      t39 = t5*s3l2
      t40 = t38*t39
      t41 = 128.D0/3.D0*t40
      t42 = l3*s4l2
      t43 = t34*t42
      t44 = 64.D0/3.D0*t43
      t49 = -128.D0/3.D0*t36-64.D0/3.D0*t40-t44
      t51 = l3*s34r
      t55 = l3*t6
      t58 = 128.D0/3.D0*t55*t4*s3l2
      t60 = 64.D0/3.D0*t7*t42
      t66 = t34*l3*s2l1
      t68 = t5*s3l1
      t69 = t38*t68
      t71 = l3*s4l1
      t72 = t34*t71
      t73 = 64.D0/3.D0*t72
      t95 = 64.D0/3.D0*t69+64.D0/3.D0*t72
      t97 = l3**2
      t98 = z*t97
      t99 = t25*t6
      t105 = t71*s3l2
      t106 = t34*t105
      t112 = -t40-t43
      t119 = 128.D0/3.D0*t106
      t123 = 1+z
      t137 = t5*t2
      t138 = z*s3l2
      t157 = s34r**2
      t169 = s3l2*s3l1
      t172 = s4l1*s4l2
      t175 = z**2
      s1 = ((64.D0/3.D0*t36+t41-t44)*s1vr+t49*s23r-64.D0/3.D0*t7*t51*s2l
     #2+(-t58+t60)*s34r)*s1l1+((64.D0/3.D0*t66+128.D0/3.D0*t69-t73)*s1vr
     #+(-128.D0/3.D0*t66-64.D0/3.D0*t69-t73)*s23r-64.D0/3.D0*t7*t51*s2l1
     #+(-128.D0/3.D0*t55*t4*s3l1+64.D0/3.D0*t7*t71)*s34r)*s1l2+(-t49*s2l
     #1+t95*s2l2+64.D0/3.D0*t98*t99*s34r+(t41+t44)*s3l1+64.D0/3.D0*t106-
     #128.D0/3.D0*t1*t34)*s1vr+(64.D0/3.D0*t112*s2l1-t95*s2l2+128.D0/3.D
     #0*t112*s3l1-t119-128.D0/3.D0*t34*t42*s4l1+128.D0/3.D0*t1*t123*t3*t
     #5)*s23r+((128.D0/3.D0*t7*t51+128.D0/3.D0*t34)*s2l2+64.D0/3.D0*t7*t
     #42*s34r+128.D0/3.D0*t137*t138-128.D0/3.D0*t137*s4l2)*s2l1
      t179 = s1+(64.D0/3.D0*t7*t71*s34r+128.D0/3.D0*t137*z*s3l1-128.D0/3
     #.D0*t137*s4l1)*s2l2-64.D0/3.D0*t98*t25*t2*t26*t157+((-t58-t60)*s3l
     #1-64.D0/3.D0*t7*t105+128.D0/3.D0*t1*t7)*s34r+128.D0/3.D0*t137*t169
     #-128.D0/3.D0*t137*t172+128.D0/3.D0*t175*t2*t5
      t181 = t97*t3
      t185 = s23r**2
      t188 = t97*t4
      t196 = t181*s3l2
      t197 = t181*s4l2
      t202 = t188*t39*s34r
      t204 = t2*s4l2
      t205 = t204*l3
      t216 = t5*t4
      t222 = l3*t2
      t231 = t97*t6
      t232 = t25*t157
      t236 = t97*t25
      t247 = t181*s3l1
      t248 = t181*s4l1
      t253 = t188*t68*s34r
      t256 = t2*s4l1*l3
      t292 = s1vr**2
      t296 = t97*l3
      t298 = z*t185
      t301 = 128.D0/3.D0*t196
      t302 = 64.D0/3.D0*t197
      t311 = t25*t5*s34r
      t316 = s3l2*s4l1
      t339 = 2*t2
      t342 = t38*s4l2
      t354 = z-2
      t355 = t354*s3l2
      t372 = t38*t5
      t376 = 64.D0/3.D0*t2
      s3 = (-128.D0/3.D0*t1*t7*s1vr+(-256.D0/3.D0*t11*s2l2-128.D0/3.D0*t
     #11*s3l2-128.D0/3.D0*t11*s4l2)*s2l1+(-128.D0/3.D0*t11*s3l1-128.D0/3
     #.D0*t11*s4l1)*s2l2+128.D0/3.D0*t1*t25*t26*s34r)*t32+t179*s13r
      s2 = s3+((128.D0/3.D0*t181*s23r*s1vr-256.D0/3.D0*t181*t185-256.D0/
     #3.D0*t188*t5*s34r*s23r)*s1l2+((64.D0/3.D0*t181*s2l2+64.D0/3.D0*t19
     #6+64.D0/3.D0*t197)*s23r-64.D0/3.D0*t202-64.D0/3.D0*t205)*s1vr+(-12
     #8.D0/3.D0*t196-128.D0/3.D0*t197)*t185+(64.D0/3.D0*t38*s2l2+(-64.D0
     #/3.D0*t188*t39-128.D0/3.D0*t216*t97*s4l2)*s34r-64.D0/3.D0*t222*t12
     #3*s3l2+64.D0/3.D0*t222*t123*s4l2)*s23r+64.D0/3.D0*t231*t232*s2l2+6
     #4.D0/3.D0*t236*t6*s3l2*t157+(-64.D0/3.D0*t40+64.D0/3.D0*t43)*s34r)
     #*s1l1+(((64.D0/3.D0*t181*s2l1+64.D0/3.D0*t247+64.D0/3.D0*t248)*s23
     #r-64.D0/3.D0*t253-64.D0/3.D0*t256)*s1vr+(-128.D0/3.D0*t247-128.D0/
     #3.D0*t248)*t185+(64.D0/3.D0*t38*s2l1+(-64.D0/3.D0*t188*t68-128.D0/
     #3.D0*t216*t97*s4l1)*s34r-64.D0/3.D0*t222*t123*s3l1+64.D0/3.D0*t222
     #*t123*s4l1)*s23r+64.D0/3.D0*t231*t232*s2l1+64.D0/3.D0*t236*t6*s3l1
     #*t157+(-64.D0/3.D0*t69+64.D0/3.D0*t72)*s34r)*s1l2
      s1 = s2-128.D0/3.D0*t181*t169*t292+(128.D0/3.D0*t296*t4*t298+((t30
     #1-t302)*s2l1+(128.D0/3.D0*t247-64.D0/3.D0*t248)*s2l2+64.D0/3.D0*z*
     #t296*t311+(t301+t302)*s3l1+64.D0/3.D0*t181*t316-128.D0/3.D0*t181*z
     #)*s23r+(64.D0/3.D0*t202+64.D0/3.D0*t205)*s2l1+(64.D0/3.D0*t253+64.
     #D0/3.D0*t256)*s2l2+128.D0/3.D0*t188*t5*t169*s34r-128.D0/3.D0*t38*t
     #169)*s1vr+128.D0/3.D0*t188*t298+((64.D0/3.D0*t222*t339*s3l2-64.D0/
     #3.D0*t342)*s2l1+(64.D0/3.D0*t222*t339*s3l1-64.D0/3.D0*t38*s4l1)*s2
     #l2+64.D0/3.D0*t98*t311+(128.D0/3.D0*t222*t355+64.D0/3.D0*t342)*s3l
     #1+64.D0/3.D0*t38*t316+128.D0/3.D0*t222*t172-64.D0/3.D0*t1*t2*(4*z-
     #2))*s23r
      taux14 = s1+(64.D0/3.D0*t372*t355*s34r+t376*s4l2)*s2l1+(64.D0/3.D0
     #*t372*t354*s3l1*s34r+t376*s4l1)*s2l2-64.D0/3.D0*t98*t99*t157+((128
     #.D0/3.D0*t38*t5*t354*s3l2-128.D0/3.D0*t43)*s3l1-t119)*s34r+(-128.D
     #0/3.D0*t138+64.D0/3.D0*t204)*s3l1+64.D0/3.D0*t2*s3l2*s4l1-128.D0/3
     #.D0*t172+128.D0/3.D0*t175
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t5 = l3*s4l2
      t12 = l3*t3
      t13 = t1*s3l2
      t14 = t12*t13
      t16 = t5*t4
      t20 = t1**2
      t21 = t3*t2
      t22 = t20*t21
      t26 = 1+z
      t30 = -512.D0/3.D0*t22*l3*s34r+256*t2*t26*t1
      t33 = t22*t5*s34r
      t35 = -t2
      t37 = 3*t35*t1
      t44 = l3*s4l1
      t51 = t1*s3l1
      t52 = t12*t51
      t54 = t4*t44
      t61 = t22*t44*s34r
      t70 = -256.D0/3.D0*t14-256.D0/3.D0*t16
      t73 = 256.D0/3.D0*t52+256.D0/3.D0*t54
      t75 = t5*s3l1
      t76 = t4*t75
      t78 = t44*s3l2
      t79 = t4*t78
      t90 = 2*z
      t92 = (t35*z-t90+2)*t1
      t110 = l3**2
      t111 = t110*t3
      t117 = t110*t21
      t119 = t117*t13*s34r
      t120 = 256.D0/3.D0*t119
      t121 = t12*s3l2
      t129 = l3*t2
      t135 = 512.D0/3.D0*t110*t1*t21*s34r-256*t129*t26
      t138 = 256*t129*s3l2
      t139 = t2*s4l2
      t145 = t3**2
      t146 = s34r**2
      t157 = -256.D0/3.D0*t110*t20*t145*t146+256.D0/3.D0*(-t12*z-2*t12)*
     #t1*s34r+256+256*z
      t159 = t110*t145
      t164 = z*l3
      t168 = -256.D0/3.D0*t164*t4*s3l2+256.D0/3.D0*t16
      t171 = 256*s4l2
      t180 = t117*t51*s34r
      t181 = 256.D0/3.D0*t180
      t182 = t12*s3l1
      t189 = 256*t129*s3l1
      t190 = t2*s4l1
      t204 = 256.D0/3.D0*t164*t4*s3l1-256.D0/3.D0*t54
      t207 = 256*s4l1
      t214 = s4l2*s3l1
      t216 = s3l2*s4l1
      t237 = -t164*t2-2*t129
      t254 = z**2
      t255 = t90-t254+2
      s2 = ((-256.D0/3.D0*t4*t5*s1vr+(-512.D0/3.D0*t4*l3*s2l2-256.D0/3.D
     #0*t14-256.D0/3.D0*t16)*s23r+t30*s2l2+256.D0/3.D0*t33-256.D0/3.D0*t
     #37*s3l2-256.D0/3.D0*t37*s4l2)*s1l1+(256.D0/3.D0*t4*t44*s1vr+(512.D
     #0/3.D0*t4*l3*s2l1+256.D0/3.D0*t52+256.D0/3.D0*t54)*s23r-t30*s2l1-2
     #56.D0/3.D0*t61+256.D0/3.D0*t37*s3l1+256.D0/3.D0*t37*s4l1)*s1l2+(t7
     #0*s2l1+t73*s2l2-256.D0/3.D0*t76+256.D0/3.D0*t79)*s1vr+(-t70*s2l1-t
     #73*s2l2)*s23r+(256.D0/3.D0*t33+256.D0/3.D0*t92*s3l2+256.D0/3.D0*t9
     #2*s4l2)*s2l1+(-256.D0/3.D0*t61-256.D0/3.D0*t92*s3l1-256.D0/3.D0*t9
     #2*s4l1)*s2l2+(256.D0/3.D0*t22*t75-256.D0/3.D0*t22*t78)*s34r)*s13r+
     #(((256.D0/3.D0*t111*s3l2+256.D0/3.D0*t111*s4l2)*s23r+t120+256.D0/3
     #.D0*t121)*s1vr+(t135*s2l2-t120-t138-256*t139*l3)*s23r+t157*s2l2-25
     #6.D0/3.D0*t159*t20*s3l2*t146+t168*s34r+256*s3l2+t171)*s1l1
      s1 = s2+(((-256.D0/3.D0*t111*s3l1-256.D0/3.D0*t111*s4l1)*s23r-t181
     #-256.D0/3.D0*t182)*s1vr+(-t135*s2l1+t181+t189+256*t190*l3)*s23r-t1
     #57*s2l1+256.D0/3.D0*t159*t20*s3l1*t146+t204*s34r-256*s3l1-t207)*s1
     #l2+((256.D0/3.D0*t111*s4l2*s2l1-256.D0/3.D0*t111*s4l1*s2l2+256.D0/
     #3.D0*t111*t214-256.D0/3.D0*t111*t216)*s23r+(256.D0/3.D0*t119+256.D
     #0/3.D0*t121)*s2l1+(-256.D0/3.D0*t180-256.D0/3.D0*t182)*s2l2+512.D0
     #/3.D0*t139*l3*s3l1-512.D0/3.D0*t190*l3*s3l2)*s1vr+((t138-256.D0/3.
     #D0*t237*s4l2)*s2l1+(-t189+256.D0/3.D0*t237*s4l1)*s2l2+256.D0/3.D0*
     #t12*t214-256.D0/3.D0*t12*t216)*s23r
      taux15 = s1+(-t168*s34r-256.D0/3.D0*t255*s3l2-t171)*s2l1+(-t204*s3
     #4r+256.D0/3.D0*t255*s3l1+t207)*s2l2+(-256.D0/3.D0*t76+256.D0/3.D0*
     #t79)*s34r-256.D0/3.D0*t35*s4l2*s3l1+256.D0/3.D0*t35*s4l1*s3l2
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t5 = l3*s4l2
      t12 = l3*t3
      t13 = t1*s3l2
      t14 = t12*t13
      t16 = t5*t4
      t20 = t1**2
      t21 = t3*t2
      t22 = t20*t21
      t26 = 1+z
      t30 = 512.D0/3.D0*t22*l3*s34r-256*t2*t26*t1
      t33 = t22*t5*s34r
      t35 = -t2
      t37 = 3*t35*t1
      t44 = l3*s4l1
      t51 = t1*s3l1
      t52 = t12*t51
      t54 = t4*t44
      t61 = t22*t44*s34r
      t70 = 256.D0/3.D0*t14+256.D0/3.D0*t16
      t73 = -256.D0/3.D0*t52-256.D0/3.D0*t54
      t75 = t5*s3l1
      t76 = t4*t75
      t78 = t44*s3l2
      t79 = t4*t78
      t90 = 2*z
      t92 = (t35*z-t90+2)*t1
      t110 = l3**2
      t111 = t110*t3
      t117 = t110*t21
      t119 = t117*t13*s34r
      t120 = 256.D0/3.D0*t119
      t121 = t12*s3l2
      t129 = l3*t2
      t135 = -512.D0/3.D0*t110*t1*t21*s34r+256*t129*t26
      t138 = 256*t129*s3l2
      t139 = t2*s4l2
      t145 = t3**2
      t146 = s34r**2
      t157 = 256.D0/3.D0*t110*t20*t145*t146-256.D0/3.D0*(-t12*z-2*t12)*t
     #1*s34r-256-256*z
      t159 = t110*t145
      t164 = z*l3
      t168 = 256.D0/3.D0*t164*t4*s3l2-256.D0/3.D0*t16
      t171 = 256*s4l2
      t180 = t117*t51*s34r
      t181 = 256.D0/3.D0*t180
      t182 = t12*s3l1
      t189 = 256*t129*s3l1
      t190 = t2*s4l1
      t204 = -256.D0/3.D0*t164*t4*s3l1+256.D0/3.D0*t54
      t207 = 256*s4l1
      t214 = s4l2*s3l1
      t216 = s3l2*s4l1
      t237 = -t164*t2-2*t129
      t254 = z**2
      t255 = t90-t254+2
      s2 = ((256.D0/3.D0*t4*t5*s1vr+(512.D0/3.D0*t4*l3*s2l2+256.D0/3.D0*
     #t14+256.D0/3.D0*t16)*s23r+t30*s2l2-256.D0/3.D0*t33+256.D0/3.D0*t37
     #*s3l2+256.D0/3.D0*t37*s4l2)*s1l1+(-256.D0/3.D0*t4*t44*s1vr+(-512.D
     #0/3.D0*t4*l3*s2l1-256.D0/3.D0*t52-256.D0/3.D0*t54)*s23r-t30*s2l1+2
     #56.D0/3.D0*t61-256.D0/3.D0*t37*s3l1-256.D0/3.D0*t37*s4l1)*s1l2+(t7
     #0*s2l1+t73*s2l2+256.D0/3.D0*t76-256.D0/3.D0*t79)*s1vr+(-t70*s2l1-t
     #73*s2l2)*s23r+(-256.D0/3.D0*t33-256.D0/3.D0*t92*s3l2-256.D0/3.D0*t
     #92*s4l2)*s2l1+(256.D0/3.D0*t61+256.D0/3.D0*t92*s3l1+256.D0/3.D0*t9
     #2*s4l1)*s2l2+(-256.D0/3.D0*t22*t75+256.D0/3.D0*t22*t78)*s34r)*s13r
     #+(((-256.D0/3.D0*t111*s3l2-256.D0/3.D0*t111*s4l2)*s23r-t120-256.D0
     #/3.D0*t121)*s1vr+(t135*s2l2+t120+t138+256*t139*l3)*s23r+t157*s2l2+
     #256.D0/3.D0*t159*t20*s3l2*t146+t168*s34r-256*s3l2-t171)*s1l1
      s1 = s2+(((256.D0/3.D0*t111*s3l1+256.D0/3.D0*t111*s4l1)*s23r+t181+
     #256.D0/3.D0*t182)*s1vr+(-t135*s2l1-t181-t189-256*t190*l3)*s23r-t15
     #7*s2l1-256.D0/3.D0*t159*t20*s3l1*t146+t204*s34r+256*s3l1+t207)*s1l
     #2+((-256.D0/3.D0*t111*s4l2*s2l1+256.D0/3.D0*t111*s4l1*s2l2-256.D0/
     #3.D0*t111*t214+256.D0/3.D0*t111*t216)*s23r+(-256.D0/3.D0*t119-256.
     #D0/3.D0*t121)*s2l1+(256.D0/3.D0*t180+256.D0/3.D0*t182)*s2l2-512.D0
     #/3.D0*t139*l3*s3l1+512.D0/3.D0*t190*l3*s3l2)*s1vr+((-t138+256.D0/3
     #.D0*t237*s4l2)*s2l1+(t189-256.D0/3.D0*t237*s4l1)*s2l2-256.D0/3.D0*
     #t12*t214+256.D0/3.D0*t12*t216)*s23r
      taux16 = s1+(-t168*s34r+256.D0/3.D0*t255*s3l2+t171)*s2l1+(-t204*s3
     #4r-256.D0/3.D0*t255*s3l1-t207)*s2l2+(256.D0/3.D0*t76-256.D0/3.D0*t
     #79)*s34r+256.D0/3.D0*t35*s4l2*s3l1-256.D0/3.D0*t35*s4l1*s3l2

CCC Coefficients

      tcoff41 = taux13
      tcoff42 = taux13
      tcoff43 = taux14
      tcoff44 = taux14
      tcoff45 = taux15
      tcoff46 = taux15
      tcoff47 = taux16
      tcoff48 = taux16

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet42 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet43 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet44 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet45 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet46 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet47 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet48 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff48)

CCC Sum

      sq3num6e1 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet4
     #4*tcoff44+tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet48*
     #tcoff48
      RETURN
      END

      double precision function sq3num6e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3*t3
      t20 = t5*s3l2
      t24 = t5*t3
      t28 = t19*t20
      t31 = t24*l3*s4l2
      t32 = 64.D0/3.D0*t31
      t35 = l3*s34r
      t39 = l3*t6
      t40 = t4*s3l2
      t46 = t5*s3l1
      t53 = t19*t46
      t55 = l3*s4l1
      t56 = t24*t55
      t69 = l3**2
      t70 = t69*t4
      t71 = z*t5
      t72 = s1vr**2
      t76 = z*t69
      t77 = t10*t6
      t79 = t76*t77*s34r
      t84 = t55*s3l2
      t85 = t24*t84
      t86 = 64.D0/3.D0*t85
      t93 = -t28-t31
      t114 = t39*t4
      t115 = s3l2*s34r
      t116 = t115*s2l1
      t119 = s3l1*s34r
      t120 = t119*s2l2
      t125 = s34r**2
      t146 = t70*t5
      t150 = t69*t5
      t151 = t4*s34r
      t159 = t69*t6
      t160 = t10*t125
      t164 = t69*t10
      t190 = t69*l3
      t194 = t69*t3
      t195 = s3l2*s3l1
      t201 = z*t190
      t202 = s23r**2
      t203 = t4*t202
      t206 = s3l2*s2l1
      t209 = s3l1*s2l2
      t222 = s3l2*s4l1
      t269 = z**2
      t276 = l3*t4
      t277 = t276*t5
      s3 = (-64*t1*t7*s1vr+64*t1*t10*t11*s34r)*t17
      s5 = ((64.D0/3.D0*t19*t20*s1vr+(-128.D0/3.D0*t24*l3*s2l2-64.D0/3.D
     #0*t28-t32)*s23r-64.D0/3.D0*t7*t35*s2l2-64.D0/3.D0*t39*t40*s34r)*s1
     #l1+(64.D0/3.D0*t19*t46*s1vr+(-128.D0/3.D0*t24*l3*s2l1-64.D0/3.D0*t
     #53-64.D0/3.D0*t56)*s23r-64.D0/3.D0*t7*t35*s2l1-64.D0/3.D0*t39*t4*s
     #3l1*s34r)*s1l2+128.D0/3.D0*t70*t71*t72+(64.D0/3.D0*t79+(128.D0/3.D
     #0*t28+t32)*s3l1+t86+64.D0/3.D0*t1*t24*(-3+z))*s1vr+(64.D0/3.D0*t93
     #*s2l1+(-64.D0/3.D0*t53-64.D0/3.D0*t56)*s2l2-320.D0/3.D0*t79+128.D0
     #/3.D0*t93*s3l1-128.D0/3.D0*t85-128.D0/3.D0*t19*t5*s4l2*s4l1+64*t1*
     #(1+z)*t3*t5)*s23r-64.D0/3.D0*t114*t116-64.D0/3.D0*t114*t120-64*t76
     #*t10*t2*t11*t125+((-128.D0/3.D0*t39*t40-64.D0/3.D0*t39*t4*s4l2)*s3
     #l1-64.D0/3.D0*t7*t84+64.D0/3.D0*t1*t7*(-z+4))*s34r)*s13r
      s6 = (-64.D0/3.D0*t146*t115*s1vr+(64.D0/3.D0*t150*t151*s2l2+64.D0/
     #3.D0*t70*t20*s34r)*s23r+64.D0/3.D0*t159*t160*s2l2+64.D0/3.D0*t164*
     #t6*s3l2*t125)*s1l1
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-64.D0/3.D0*t146*t119*s1vr+(64.D0/3.D0*t150*t151*s2l1+64.
     #D0/3.D0*t70*t46*s34r)*s23r+64.D0/3.D0*t159*t160*s2l1+64.D0/3.D0*t1
     #64*t6*s3l1*t125)*s1l2+(-128.D0/3.D0*t190*t4*z*s23r-128.D0/3.D0*t19
     #4*t195+128.D0/3.D0*t194*z)*t72+(64*t201*t203+(64.D0/3.D0*t194*t206
     #+64.D0/3.D0*t194*t209+320.D0/3.D0*t201*t10*t5*s34r+(128.D0/3.D0*t1
     #94*s3l2+64.D0/3.D0*t194*s4l2)*s3l1+64.D0/3.D0*t194*t222-64.D0/3.D0
     #*t76*t3*(-1+3*z))*s23r+64.D0/3.D0*t146*t116+64.D0/3.D0*t146*t120+(
     #128.D0/3.D0*t70*t20*s3l1-128.D0/3.D0*t70*t71)*s34r-128.D0/3.D0*t19
     #*t195+128.D0/3.D0*t19*z)*s1vr
      taux13 = s1+64*t76*t203+(64.D0/3.D0*t19*t206+64.D0/3.D0*t19*t209+6
     #4.D0/3.D0*t76*t4*(5*z-6)*t5*s34r+(128.D0/3.D0*t19*s3l2+64.D0/3.D0*
     #t19*s4l2)*s3l1+64.D0/3.D0*t19*t222-64.D0/3.D0*t1*t2*(t269-1))*s23r
     #+64.D0/3.D0*t277*t116+64.D0/3.D0*t277*t120-64*t76*t77*t125+((128.D
     #0/3.D0*t276*t20-t32)*s3l1-t86-64.D0/3.D0*t1*t4*t5)*s34r
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3*t3
      t20 = t5*s3l2
      t24 = t5*t3
      t28 = t19*t20
      t31 = t24*l3*s4l2
      t32 = 64.D0/3.D0*t31
      t35 = l3*s34r
      t39 = l3*t6
      t40 = t4*s3l2
      t46 = t5*s3l1
      t53 = t19*t46
      t55 = l3*s4l1
      t56 = t24*t55
      t69 = l3**2
      t70 = t69*t4
      t71 = z*t5
      t72 = s1vr**2
      t76 = z*t69
      t77 = t10*t6
      t79 = t76*t77*s34r
      t84 = t55*s3l2
      t85 = t24*t84
      t86 = 64.D0/3.D0*t85
      t93 = -t28-t31
      t114 = t39*t4
      t115 = s3l2*s34r
      t116 = t115*s2l1
      t119 = s3l1*s34r
      t120 = t119*s2l2
      t125 = s34r**2
      t146 = t70*t5
      t150 = t69*t5
      t151 = t4*s34r
      t159 = t69*t6
      t160 = t10*t125
      t164 = t69*t10
      t190 = t69*l3
      t194 = t69*t3
      t195 = s3l2*s3l1
      t201 = z*t190
      t202 = s23r**2
      t203 = t4*t202
      t206 = s3l2*s2l1
      t209 = s3l1*s2l2
      t222 = s3l2*s4l1
      t269 = z**2
      t276 = l3*t4
      t277 = t276*t5
      s3 = (-64*t1*t7*s1vr+64*t1*t10*t11*s34r)*t17
      s5 = ((64.D0/3.D0*t19*t20*s1vr+(-128.D0/3.D0*t24*l3*s2l2-64.D0/3.D
     #0*t28-t32)*s23r-64.D0/3.D0*t7*t35*s2l2-64.D0/3.D0*t39*t40*s34r)*s1
     #l1+(64.D0/3.D0*t19*t46*s1vr+(-128.D0/3.D0*t24*l3*s2l1-64.D0/3.D0*t
     #53-64.D0/3.D0*t56)*s23r-64.D0/3.D0*t7*t35*s2l1-64.D0/3.D0*t39*t4*s
     #3l1*s34r)*s1l2+128.D0/3.D0*t70*t71*t72+(64.D0/3.D0*t79+(128.D0/3.D
     #0*t28+t32)*s3l1+t86+64.D0/3.D0*t1*t24*(-3+z))*s1vr+(64.D0/3.D0*t93
     #*s2l1+(-64.D0/3.D0*t53-64.D0/3.D0*t56)*s2l2-320.D0/3.D0*t79+128.D0
     #/3.D0*t93*s3l1-128.D0/3.D0*t85-128.D0/3.D0*t19*t5*s4l2*s4l1+64*t1*
     #(1+z)*t3*t5)*s23r-64.D0/3.D0*t114*t116-64.D0/3.D0*t114*t120-64*t76
     #*t10*t2*t11*t125+((-128.D0/3.D0*t39*t40-64.D0/3.D0*t39*t4*s4l2)*s3
     #l1-64.D0/3.D0*t7*t84+64.D0/3.D0*t1*t7*(-z+4))*s34r)*s13r
      s6 = (-64.D0/3.D0*t146*t115*s1vr+(64.D0/3.D0*t150*t151*s2l2+64.D0/
     #3.D0*t70*t20*s34r)*s23r+64.D0/3.D0*t159*t160*s2l2+64.D0/3.D0*t164*
     #t6*s3l2*t125)*s1l1
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-64.D0/3.D0*t146*t119*s1vr+(64.D0/3.D0*t150*t151*s2l1+64.
     #D0/3.D0*t70*t46*s34r)*s23r+64.D0/3.D0*t159*t160*s2l1+64.D0/3.D0*t1
     #64*t6*s3l1*t125)*s1l2+(-128.D0/3.D0*t190*t4*z*s23r-128.D0/3.D0*t19
     #4*t195+128.D0/3.D0*t194*z)*t72+(64*t201*t203+(64.D0/3.D0*t194*t206
     #+64.D0/3.D0*t194*t209+320.D0/3.D0*t201*t10*t5*s34r+(128.D0/3.D0*t1
     #94*s3l2+64.D0/3.D0*t194*s4l2)*s3l1+64.D0/3.D0*t194*t222-64.D0/3.D0
     #*t76*t3*(-1+3*z))*s23r+64.D0/3.D0*t146*t116+64.D0/3.D0*t146*t120+(
     #128.D0/3.D0*t70*t20*s3l1-128.D0/3.D0*t70*t71)*s34r-128.D0/3.D0*t19
     #*t195+128.D0/3.D0*t19*z)*s1vr
      taux14 = s1+64*t76*t203+(64.D0/3.D0*t19*t206+64.D0/3.D0*t19*t209+6
     #4.D0/3.D0*t76*t4*(5*z-6)*t5*s34r+(128.D0/3.D0*t19*s3l2+64.D0/3.D0*
     #t19*s4l2)*s3l1+64.D0/3.D0*t19*t222-64.D0/3.D0*t1*t2*(t269-1))*s23r
     #+64.D0/3.D0*t277*t116+64.D0/3.D0*t277*t120-64*t76*t77*t125+((128.D
     #0/3.D0*t276*t20-t32)*s3l1-t86-64.D0/3.D0*t1*t4*t5)*s34r
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t6 = t4*l3*s2l2
      t8 = l3*t3
      t9 = t1*s3l2
      t10 = t8*t9
      t12 = l3*s4l2
      t13 = t4*t12
      t22 = t1**2
      t23 = t3*t2
      t24 = t22*t23
      t27 = 1+z
      t31 = 512.D0/3.D0*t24*l3*s34r-512.D0/3.D0*t2*t27*t1
      t33 = l3*t22
      t35 = t33*t23*s3l2
      t36 = t23*s4l2
      t37 = t33*t36
      t41 = t1*t2
      t49 = t4*l3*s2l1
      t51 = t1*s3l1
      t52 = t8*t51
      t54 = l3*s4l1
      t55 = t4*t54
      t67 = t33*t23*s3l1
      t70 = l3*t23*t22*s4l1
      t80 = t10+t13
      t83 = -t52-t55
      t87 = t4*t12*s3l1
      t89 = t54*s3l2
      t90 = t4*t89
      t106 = t2*z
      t109 = t1*s4l2
      t120 = t1*s4l1
      t133 = l3**2
      t134 = t133*t3
      t137 = t134*s3l2
      t139 = t134*s4l2
      t143 = t133*t23
      t146 = 256.D0/3.D0*t143*t9*s34r
      t148 = 512.D0/3.D0*t8*s3l2
      t149 = l3*t2
      t150 = t149*s4l2
      t151 = 256.D0/3.D0*t150
      t161 = -256.D0/3.D0*t133*t1*t23*s34r-256.D0/3.D0*t8+512.D0/3.D0*t1
     #49*t27
      t174 = t3**2
      t175 = s34r**2
      t179 = z*l3
      t184 = -256.D0/3.D0*t133*t22*t174*t175+512.D0/3.D0*t179*t4*s34r-51
     #2.D0/3.D0-512.D0/3.D0*z
      t186 = t133*t174
      t202 = t134*s3l1
      t204 = t134*s4l1
      t210 = 256.D0/3.D0*t143*t51*s34r
      t212 = 512.D0/3.D0*t8*s3l1
      t213 = t149*s4l1
      t214 = 256.D0/3.D0*t213
      t250 = s4l2*s3l1
      t253 = s3l2*s4l1
      t271 = -2*t149*t2+2*t149
      t297 = (-3*t8*t2+2*t8*z)*t1
      t303 = z**2
      t304 = -2*z+t303
      t307 = -t27
      s3 = (((-256.D0/3.D0*t6-512.D0/3.D0*t10+512.D0/3.D0*t13)*s1vr+(204
     #8.D0/3.D0*t6+1024.D0/3.D0*t10+1024.D0/3.D0*t13)*s23r+t31*s2l2+(512
     #.D0/3.D0*t35-512.D0/3.D0*t37)*s34r-512.D0/3.D0*t41*s3l2-512.D0/3.D
     #0*t41*s4l2)*s1l1+((256.D0/3.D0*t49+512.D0/3.D0*t52-512.D0/3.D0*t55
     #)*s1vr+(-2048.D0/3.D0*t49-1024.D0/3.D0*t52-1024.D0/3.D0*t55)*s23r-
     #t31*s2l1+(-512.D0/3.D0*t67+512.D0/3.D0*t70)*s34r+512.D0/3.D0*t41*s
     #3l1+512.D0/3.D0*t41*s4l1)*s1l2+(512.D0/3.D0*t80*s2l1+512.D0/3.D0*t
     #83*s2l2+256.D0/3.D0*t87-256.D0/3.D0*t90)*s1vr+(-1024.D0/3.D0*t80*s
     #2l1-1024.D0/3.D0*t83*s2l2)*s23r+((-256*t35-256.D0/3.D0*t37)*s34r+5
     #12.D0/3.D0*t106*t9+512.D0/3.D0*t106*t109)*s2l1+((256*t67+256.D0/3.
     #D0*t70)*s34r-512.D0/3.D0*t106*t51-512.D0/3.D0*t106*t120)*s2l2+(-25
     #6.D0/3.D0*t33*t36*s3l1+256.D0/3.D0*t24*t89)*s34r)*s13r
      s4 = (((-256.D0/3.D0*t134*s2l2-512.D0/3.D0*t137-512.D0/3.D0*t139)*
     #s23r+t146-t148+t151)*s1vr+(t161*s2l2+(-256.D0/3.D0*t143*t9+256.D0/
     #3.D0*t143*t109)*s34r+512.D0/3.D0*t149*s3l2+512.D0/3.D0*t150)*s23r+
     #t184*s2l2-256.D0/3.D0*t186*t22*s3l2*t175+(512.D0/3.D0*t179*t4*s3l2
     #-512.D0/3.D0*t13)*s34r-512.D0/3.D0*s3l2-512.D0/3.D0*s4l2)*s1l1
      s2 = s3+s4
      s1 = s2+(((256.D0/3.D0*t134*s2l1+512.D0/3.D0*t202+512.D0/3.D0*t204
     #)*s23r-t210+t212-t214)*s1vr+(-t161*s2l1+(256.D0/3.D0*t143*t51-256.
     #D0/3.D0*t143*t120)*s34r-512.D0/3.D0*t149*s3l1-512.D0/3.D0*t213)*s2
     #3r-t184*s2l1+256.D0/3.D0*t186*t22*s3l1*t175+(-512.D0/3.D0*t179*t4*
     #s3l1+512.D0/3.D0*t55)*s34r+512.D0/3.D0*s3l1+512.D0/3.D0*s4l1)*s1l2
     #+(((512.D0/3.D0*t137-512.D0/3.D0*t139)*s2l1+(-512.D0/3.D0*t202+512
     #.D0/3.D0*t204)*s2l2-256.D0/3.D0*t134*t250+256.D0/3.D0*t134*t253)*s
     #23r+(t146-t148+t151)*s2l1+(-t210+t212-t214)*s2l2-512.D0/3.D0*t149*
     #t250+512.D0/3.D0*t149*t253)*s1vr+((-256.D0/3.D0*t271*s3l2-512.D0/3
     #.D0*t149*z*s4l2)*s2l1+(256.D0/3.D0*t271*s3l1+512.D0/3.D0*t149*z*s4
     #l1)*s2l2-256.D0/3.D0*t8*t250+256.D0/3.D0*t8*t253)*s23r
      taux15 = s1+((-256.D0/3.D0*t297*s3l2+256.D0/3.D0*t13)*s34r-512.D0/
     #3.D0*t304*s3l2-256.D0/3.D0*t307*s4l2)*s2l1+((256.D0/3.D0*t297*s3l1
     #-256.D0/3.D0*t55)*s34r+512.D0/3.D0*t304*s3l1+256.D0/3.D0*t307*s4l1
     #)*s2l2+(256.D0/3.D0*t87-256.D0/3.D0*t90)*s34r-256.D0/3.D0*t2*s4l2*
     #s3l1+256.D0/3.D0*t2*s4l1*s3l2
      t1 = l3-1
      t2 = -1+z
      t3 = t2**2
      t4 = t1*t3
      t6 = t4*l3*s2l2
      t8 = l3*t3
      t9 = t1*s3l2
      t10 = t8*t9
      t12 = l3*s4l2
      t13 = t4*t12
      t22 = t1**2
      t23 = t3*t2
      t24 = t22*t23
      t27 = 1+z
      t31 = -512.D0/3.D0*t24*l3*s34r+512.D0/3.D0*t2*t27*t1
      t33 = l3*t22
      t35 = t33*t23*s3l2
      t36 = t23*s4l2
      t37 = t33*t36
      t41 = t1*t2
      t49 = t4*l3*s2l1
      t51 = t1*s3l1
      t52 = t8*t51
      t54 = l3*s4l1
      t55 = t4*t54
      t67 = t33*t23*s3l1
      t70 = l3*t23*t22*s4l1
      t80 = -t10-t13
      t83 = t52+t55
      t87 = t4*t12*s3l1
      t89 = t54*s3l2
      t90 = t4*t89
      t106 = t2*z
      t109 = t1*s4l2
      t120 = t1*s4l1
      t133 = l3**2
      t134 = t133*t3
      t137 = t134*s3l2
      t139 = t134*s4l2
      t143 = t133*t23
      t146 = 256.D0/3.D0*t143*t9*s34r
      t148 = 512.D0/3.D0*t8*s3l2
      t149 = l3*t2
      t150 = t149*s4l2
      t151 = 256.D0/3.D0*t150
      t161 = 256.D0/3.D0*t133*t1*t23*s34r+256.D0/3.D0*t8-512.D0/3.D0*t14
     #9*t27
      t174 = t3**2
      t175 = s34r**2
      t179 = z*l3
      t184 = 256.D0/3.D0*t133*t22*t174*t175-512.D0/3.D0*t179*t4*s34r+512
     #.D0/3.D0+512.D0/3.D0*z
      t186 = t133*t174
      t202 = t134*s3l1
      t204 = t134*s4l1
      t210 = 256.D0/3.D0*t143*t51*s34r
      t212 = 512.D0/3.D0*t8*s3l1
      t213 = t149*s4l1
      t214 = 256.D0/3.D0*t213
      t250 = s4l2*s3l1
      t253 = s3l2*s4l1
      t271 = -2*t149*t2+2*t149
      t297 = (-3*t8*t2+2*t8*z)*t1
      t303 = z**2
      t304 = -2*z+t303
      t307 = -t27
      s3 = (((256.D0/3.D0*t6+512.D0/3.D0*t10-512.D0/3.D0*t13)*s1vr+(-204
     #8.D0/3.D0*t6-1024.D0/3.D0*t10-1024.D0/3.D0*t13)*s23r+t31*s2l2+(-51
     #2.D0/3.D0*t35+512.D0/3.D0*t37)*s34r+512.D0/3.D0*t41*s3l2+512.D0/3.
     #D0*t41*s4l2)*s1l1+((-256.D0/3.D0*t49-512.D0/3.D0*t52+512.D0/3.D0*t
     #55)*s1vr+(2048.D0/3.D0*t49+1024.D0/3.D0*t52+1024.D0/3.D0*t55)*s23r
     #-t31*s2l1+(512.D0/3.D0*t67-512.D0/3.D0*t70)*s34r-512.D0/3.D0*t41*s
     #3l1-512.D0/3.D0*t41*s4l1)*s1l2+(512.D0/3.D0*t80*s2l1+512.D0/3.D0*t
     #83*s2l2-256.D0/3.D0*t87+256.D0/3.D0*t90)*s1vr+(-1024.D0/3.D0*t80*s
     #2l1-1024.D0/3.D0*t83*s2l2)*s23r+((256*t35+256.D0/3.D0*t37)*s34r-51
     #2.D0/3.D0*t106*t9-512.D0/3.D0*t106*t109)*s2l1+((-256*t67-256.D0/3.
     #D0*t70)*s34r+512.D0/3.D0*t106*t51+512.D0/3.D0*t106*t120)*s2l2+(256
     #.D0/3.D0*t33*t36*s3l1-256.D0/3.D0*t24*t89)*s34r)*s13r
      s4 = (((256.D0/3.D0*t134*s2l2+512.D0/3.D0*t137+512.D0/3.D0*t139)*s
     #23r-t146+t148-t151)*s1vr+(t161*s2l2+(256.D0/3.D0*t143*t9-256.D0/3.
     #D0*t143*t109)*s34r-512.D0/3.D0*t149*s3l2-512.D0/3.D0*t150)*s23r+t1
     #84*s2l2+256.D0/3.D0*t186*t22*s3l2*t175+(-512.D0/3.D0*t179*t4*s3l2+
     #512.D0/3.D0*t13)*s34r+512.D0/3.D0*s3l2+512.D0/3.D0*s4l2)*s1l1
      s2 = s3+s4
      s1 = s2+(((-256.D0/3.D0*t134*s2l1-512.D0/3.D0*t202-512.D0/3.D0*t20
     #4)*s23r+t210-t212+t214)*s1vr+(-t161*s2l1+(-256.D0/3.D0*t143*t51+25
     #6.D0/3.D0*t143*t120)*s34r+512.D0/3.D0*t149*s3l1+512.D0/3.D0*t213)*
     #s23r-t184*s2l1-256.D0/3.D0*t186*t22*s3l1*t175+(512.D0/3.D0*t179*t4
     #*s3l1-512.D0/3.D0*t55)*s34r-512.D0/3.D0*s3l1-512.D0/3.D0*s4l1)*s1l
     #2+(((-512.D0/3.D0*t137+512.D0/3.D0*t139)*s2l1+(512.D0/3.D0*t202-51
     #2.D0/3.D0*t204)*s2l2+256.D0/3.D0*t134*t250-256.D0/3.D0*t134*t253)*
     #s23r+(-t146+t148-t151)*s2l1+(t210-t212+t214)*s2l2+512.D0/3.D0*t149
     #*t250-512.D0/3.D0*t149*t253)*s1vr+((256.D0/3.D0*t271*s3l2+512.D0/3
     #.D0*t149*z*s4l2)*s2l1+(-256.D0/3.D0*t271*s3l1-512.D0/3.D0*t149*z*s
     #4l1)*s2l2+256.D0/3.D0*t8*t250-256.D0/3.D0*t8*t253)*s23r
      taux16 = s1+((256.D0/3.D0*t297*s3l2-256.D0/3.D0*t13)*s34r+512.D0/3
     #.D0*t304*s3l2+256.D0/3.D0*t307*s4l2)*s2l1+((-256.D0/3.D0*t297*s3l1
     #+256.D0/3.D0*t55)*s34r-512.D0/3.D0*t304*s3l1-256.D0/3.D0*t307*s4l1
     #)*s2l2+(-256.D0/3.D0*t87+256.D0/3.D0*t90)*s34r+256.D0/3.D0*t2*s4l2
     #*s3l1-256.D0/3.D0*t2*s4l1*s3l2

CCC Coefficients

      tcoff41 = taux13
      tcoff42 = taux13
      tcoff43 = taux14
      tcoff44 = taux14
      tcoff45 = taux15
      tcoff46 = taux15
      tcoff47 = taux16
      tcoff48 = taux16

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet42 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet43 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet44 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet45 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet46 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet47 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet48 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff48)

CCC Sum

      sq3num6e2 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet4
     #4*tcoff44+tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet48*
     #tcoff48
      RETURN
      END

      double precision function sq3num6e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3**2
      t20 = t19*z
      t26 = t10*t6
      t38 = s34r**2
      t48 = z*t19*l3
      t49 = s23r**2
      t50 = t4*t49
      taux13 = (-64.D0/3.D0*t1*t7*s1vr+64.D0/3.D0*t1*t10*t11*s34r)*t17+(
     #64.D0/3.D0*t20*t10*t6*s34r*s1vr+(-64.D0/3.D0*t20*t26*s34r+64.D0/3.
     #D0*t1*(1+z)*t3*t5)*s23r-64.D0/3.D0*t20*t10*t2*t11*t38+64.D0/3.D0*t
     #1*t7*s34r)*s13r+(64.D0/3.D0*t48*t50+64.D0/3.D0*t48*t10*t5*s34r*s23
     #r)*s1vr+64.D0/3.D0*t20*t50+64.D0/3.D0*t20*t4*(z-2)*t5*s34r*s23r-64
     #.D0/3.D0*t20*t26*t38
      t1 = z*l3
      t2 = -1+z
      t3 = t2**2
      t4 = t3*t2
      t5 = l3-1
      t6 = t5**2
      t7 = t4*t6
      t10 = t3**2
      t11 = t6*t5
      t17 = s13r**2
      t19 = l3**2
      t20 = t19*z
      t26 = t10*t6
      t38 = s34r**2
      t48 = z*t19*l3
      t49 = s23r**2
      t50 = t4*t49
      taux14 = (-64.D0/3.D0*t1*t7*s1vr+64.D0/3.D0*t1*t10*t11*s34r)*t17+(
     #64.D0/3.D0*t20*t10*t6*s34r*s1vr+(-64.D0/3.D0*t20*t26*s34r+64.D0/3.
     #D0*t1*(1+z)*t3*t5)*s23r-64.D0/3.D0*t20*t10*t2*t11*t38+64.D0/3.D0*t
     #1*t7*s34r)*s13r+(64.D0/3.D0*t48*t50+64.D0/3.D0*t48*t10*t5*s34r*s23
     #r)*s1vr+64.D0/3.D0*t20*t50+64.D0/3.D0*t20*t4*(z-2)*t5*s34r*s23r-64
     #.D0/3.D0*t20*t26*t38
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 256*t9
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t34 = t3*t4*s2l1
      t36 = t4*s3l1
      t37 = t3*t36
      t38 = 256*t37
      t42 = t4*s4l1
      t43 = t3*t42
      t50 = t19*s3l1
      t61 = 512.D0/3.D0*t3*t14*s3l1-512.D0/3.D0*t3*t42*s3l2
      t73 = t19*s4l2
      t81 = t19*s4l1
      t96 = l3**2
      t97 = t96*t2
      t101 = t96*t20
      t104 = 1280.D0/3.D0*t101*t8*s34r
      t105 = l3*t1
      t107 = 512.D0/3.D0*t105*s4l2
      t115 = 256.D0/3.D0*t96*t4*t20*s34r+512.D0/3.D0*t3
      t125 = t96*t19
      t126 = t2**2
      t127 = s34r**2
      t128 = t126*t127
      t132 = t96*t126
      t143 = 1280.D0/3.D0*t101*t36*s34r
      t145 = 512.D0/3.D0*t105*s4l1
      t172 = s4l2*s3l1
      t175 = s3l2*s4l1
      t193 = -t3+4*t105*t1
      t210 = (-t26+6*t3*t1)*t4
      t216 = t1*s4l2
      t225 = t1*s4l1
      s2 = (((512.D0/3.D0*t6+t10)*s1vr+(-512*t6-t10-256*t15)*s23r-256*t2
     #1*t22*s2l2-256*t26*t27*s34r)*s1l1+((-512.D0/3.D0*t34-t38)*s1vr+(51
     #2*t34+t38+256*t43)*s23r+256*t21*t22*s2l1+256*t26*t50*s34r)*s1l2+t6
     #1*s1vr+((256*t9+256*t15)*s2l1+(-256*t37-256*t43)*s2l2)*s23r+(1280.
     #D0/3.D0*t26*t27-512.D0/3.D0*t26*t73)*s34r*s2l1+(-1280.D0/3.D0*t26*
     #t50+512.D0/3.D0*t26*t81)*s34r*s2l2+(-512.D0/3.D0*t26*t73*s3l1+512.
     #D0/3.D0*t26*t81*s3l2)*s34r)*s13r+((512.D0/3.D0*t97*s2l2*s23r-t104-
     #t107)*s1vr+(t115*s2l2+(1280.D0/3.D0*t101*t8-512.D0/3.D0*t101*t14)*
     #s34r)*s23r+1280.D0/3.D0*t125*t128*s2l2+1280.D0/3.D0*t132*t27*t127)
     #*s1l1
      s1 = s2+((-512.D0/3.D0*t97*s2l1*s23r+t143+t145)*s1vr+(-t115*s2l1+(
     #-1280.D0/3.D0*t101*t36+512.D0/3.D0*t101*t42)*s34r)*s23r-1280.D0/3.
     #D0*t125*t128*s2l1-1280.D0/3.D0*t132*t50*t127)*s1l2+((-256*t97*s3l2
     #*s2l1+256*t97*s3l1*s2l2-512.D0/3.D0*t97*t172+512.D0/3.D0*t97*t175)
     #*s23r+(-t104-t107)*s2l1+(t143+t145)*s2l2-1024.D0/3.D0*t105*t172+10
     #24.D0/3.D0*t105*t175)*s1vr+(-256.D0/3.D0*t193*s3l2*s2l1+256.D0/3.D
     #0*t193*s3l1*s2l2-512.D0/3.D0*t3*t172+512.D0/3.D0*t3*t175)*s23r
      taux15 = s1+((-256.D0/3.D0*t210*s3l2+512.D0/3.D0*t15)*s34r-512.D0/
     #3.D0*t216)*s2l1+((256.D0/3.D0*t210*s3l1-512.D0/3.D0*t43)*s34r+512.
     #D0/3.D0*t225)*s2l2+t61*s34r-512.D0/3.D0*t216*s3l1+512.D0/3.D0*t225
     #*s3l2
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t6 = t3*t4*s2l2
      t8 = t4*s3l2
      t9 = t3*t8
      t10 = 256*t9
      t14 = t4*s4l2
      t15 = t3*t14
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t34 = t3*t4*s2l1
      t36 = t4*s3l1
      t37 = t3*t36
      t38 = 256*t37
      t42 = t4*s4l1
      t43 = t3*t42
      t50 = t19*s3l1
      t61 = -512.D0/3.D0*t3*t14*s3l1+512.D0/3.D0*t3*t42*s3l2
      t73 = t19*s4l2
      t81 = t19*s4l1
      t96 = l3**2
      t97 = t96*t2
      t101 = t96*t20
      t104 = 1280.D0/3.D0*t101*t8*s34r
      t105 = l3*t1
      t107 = 512.D0/3.D0*t105*s4l2
      t115 = -256.D0/3.D0*t96*t4*t20*s34r-512.D0/3.D0*t3
      t125 = t96*t19
      t126 = t2**2
      t127 = s34r**2
      t128 = t126*t127
      t132 = t96*t126
      t143 = 1280.D0/3.D0*t101*t36*s34r
      t145 = 512.D0/3.D0*t105*s4l1
      t172 = s4l2*s3l1
      t175 = s3l2*s4l1
      t193 = -t3+4*t105*t1
      t210 = (-t26+6*t3*t1)*t4
      t216 = t1*s4l2
      t225 = t1*s4l1
      s2 = (((-512.D0/3.D0*t6-t10)*s1vr+(512*t6+t10+256*t15)*s23r+256*t2
     #1*t22*s2l2+256*t26*t27*s34r)*s1l1+((512.D0/3.D0*t34+t38)*s1vr+(-51
     #2*t34-t38-256*t43)*s23r-256*t21*t22*s2l1-256*t26*t50*s34r)*s1l2+t6
     #1*s1vr+((-256*t9-256*t15)*s2l1+(256*t37+256*t43)*s2l2)*s23r+(-1280
     #.D0/3.D0*t26*t27+512.D0/3.D0*t26*t73)*s34r*s2l1+(1280.D0/3.D0*t26*
     #t50-512.D0/3.D0*t26*t81)*s34r*s2l2+(512.D0/3.D0*t26*t73*s3l1-512.D
     #0/3.D0*t26*t81*s3l2)*s34r)*s13r+((-512.D0/3.D0*t97*s2l2*s23r+t104+
     #t107)*s1vr+(t115*s2l2+(-1280.D0/3.D0*t101*t8+512.D0/3.D0*t101*t14)
     #*s34r)*s23r-1280.D0/3.D0*t125*t128*s2l2-1280.D0/3.D0*t132*t27*t127
     #)*s1l1
      s1 = s2+((512.D0/3.D0*t97*s2l1*s23r-t143-t145)*s1vr+(-t115*s2l1+(1
     #280.D0/3.D0*t101*t36-512.D0/3.D0*t101*t42)*s34r)*s23r+1280.D0/3.D0
     #*t125*t128*s2l1+1280.D0/3.D0*t132*t50*t127)*s1l2+((256*t97*s3l2*s2
     #l1-256*t97*s3l1*s2l2+512.D0/3.D0*t97*t172-512.D0/3.D0*t97*t175)*s2
     #3r+(t104+t107)*s2l1+(-t143-t145)*s2l2+1024.D0/3.D0*t105*t172-1024.
     #D0/3.D0*t105*t175)*s1vr+(256.D0/3.D0*t193*s3l2*s2l1-256.D0/3.D0*t1
     #93*s3l1*s2l2+512.D0/3.D0*t3*t172-512.D0/3.D0*t3*t175)*s23r
      taux16 = s1+((256.D0/3.D0*t210*s3l2-512.D0/3.D0*t15)*s34r+512.D0/3
     #.D0*t216)*s2l1+((-256.D0/3.D0*t210*s3l1+512.D0/3.D0*t43)*s34r-512.
     #D0/3.D0*t225)*s2l2+t61*s34r+512.D0/3.D0*t216*s3l1-512.D0/3.D0*t225
     #*s3l2

CCC Coefficients

      tcoff41 = taux13
      tcoff42 = taux13
      tcoff43 = taux14
      tcoff44 = taux14
      tcoff45 = taux15
      tcoff46 = taux15
      tcoff47 = taux16
      tcoff48 = taux16

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet42 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet43 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet44 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet45 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet46 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet47 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet48 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff48)

CCC Sum

      sq3num6e3 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet4
     #4*tcoff44+tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet48*
     #tcoff48
      RETURN
      END

      double precision function sq3num6e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux13 = 0
      taux14 = 0
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t5 = t4*s3l2
      t12 = t3*t5
      t15 = t3*t4*s4l2
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t33 = t4*s3l1
      t40 = t3*t33
      t43 = t3*t4*s4l1
      t50 = t19*s3l1
      t64 = t26*t19
      t65 = s3l2*s34r
      t66 = t65*s2l1
      t69 = s3l1*s34r
      t70 = t69*s2l2
      t75 = l3**2
      t76 = t75*t20
      t77 = t76*t4
      t81 = t75*t4
      t82 = t20*s34r
      t90 = t75*t19
      t91 = t2**2
      t92 = s34r**2
      t93 = t91*t92
      t97 = t75*t91
      t121 = t75*t2
      t122 = s3l2*s2l1
      t124 = s3l1*s2l2
      t140 = t26*t4
      s1 = ((512.D0/3.D0*t3*t5*s1vr+(-1024.D0/3.D0*t3*t4*s2l2-512.D0/3.D
     #0*t12-512.D0/3.D0*t15)*s23r-512.D0/3.D0*t21*t22*s2l2-512.D0/3.D0*t
     #26*t27*s34r)*s1l1+(-512.D0/3.D0*t3*t33*s1vr+(1024.D0/3.D0*t3*t4*s2
     #l1+512.D0/3.D0*t40+512.D0/3.D0*t43)*s23r+512.D0/3.D0*t21*t22*s2l1+
     #512.D0/3.D0*t26*t50*s34r)*s1l2+((512.D0/3.D0*t12+512.D0/3.D0*t15)*
     #s2l1+(-512.D0/3.D0*t40-512.D0/3.D0*t43)*s2l2)*s23r+512.D0/3.D0*t64
     #*t66-512.D0/3.D0*t64*t70)*s13r+(-512.D0/3.D0*t77*t65*s1vr+(512.D0/
     #3.D0*t81*t82*s2l2+512.D0/3.D0*t76*t5*s34r)*s23r+512.D0/3.D0*t90*t9
     #3*s2l2+512.D0/3.D0*t97*t27*t92)*s1l1+(512.D0/3.D0*t77*t69*s1vr+(-5
     #12.D0/3.D0*t81*t82*s2l1-512.D0/3.D0*t76*t33*s34r)*s23r-512.D0/3.D0
     #*t90*t93*s2l1-512.D0/3.D0*t97*t50*t92)*s1l2
      taux15 = s1+((-512.D0/3.D0*t121*t122+512.D0/3.D0*t121*t124)*s23r-5
     #12.D0/3.D0*t77*t66+512.D0/3.D0*t77*t70)*s1vr+(-512.D0/3.D0*t3*t122
     #+512.D0/3.D0*t3*t124)*s23r-512.D0/3.D0*t140*t66+512.D0/3.D0*t140*t
     #70
      t1 = -1+z
      t2 = t1**2
      t3 = l3*t2
      t4 = l3-1
      t5 = t4*s3l2
      t12 = t3*t5
      t15 = t3*t4*s4l2
      t19 = t4**2
      t20 = t2*t1
      t21 = t19*t20
      t22 = l3*s34r
      t26 = l3*t20
      t27 = t19*s3l2
      t33 = t4*s3l1
      t40 = t3*t33
      t43 = t3*t4*s4l1
      t50 = t19*s3l1
      t64 = t26*t19
      t65 = s3l2*s34r
      t66 = t65*s2l1
      t69 = s3l1*s34r
      t70 = t69*s2l2
      t75 = l3**2
      t76 = t75*t20
      t77 = t76*t4
      t81 = t75*t4
      t82 = t20*s34r
      t90 = t75*t19
      t91 = t2**2
      t92 = s34r**2
      t93 = t91*t92
      t97 = t75*t91
      t121 = t75*t2
      t122 = s3l2*s2l1
      t124 = s3l1*s2l2
      t140 = t26*t4
      s1 = ((-512.D0/3.D0*t3*t5*s1vr+(1024.D0/3.D0*t3*t4*s2l2+512.D0/3.D
     #0*t12+512.D0/3.D0*t15)*s23r+512.D0/3.D0*t21*t22*s2l2+512.D0/3.D0*t
     #26*t27*s34r)*s1l1+(512.D0/3.D0*t3*t33*s1vr+(-1024.D0/3.D0*t3*t4*s2
     #l1-512.D0/3.D0*t40-512.D0/3.D0*t43)*s23r-512.D0/3.D0*t21*t22*s2l1-
     #512.D0/3.D0*t26*t50*s34r)*s1l2+((-512.D0/3.D0*t12-512.D0/3.D0*t15)
     #*s2l1+(512.D0/3.D0*t40+512.D0/3.D0*t43)*s2l2)*s23r-512.D0/3.D0*t64
     #*t66+512.D0/3.D0*t64*t70)*s13r+(512.D0/3.D0*t77*t65*s1vr+(-512.D0/
     #3.D0*t81*t82*s2l2-512.D0/3.D0*t76*t5*s34r)*s23r-512.D0/3.D0*t90*t9
     #3*s2l2-512.D0/3.D0*t97*t27*t92)*s1l1+(-512.D0/3.D0*t77*t69*s1vr+(5
     #12.D0/3.D0*t81*t82*s2l1+512.D0/3.D0*t76*t33*s34r)*s23r+512.D0/3.D0
     #*t90*t93*s2l1+512.D0/3.D0*t97*t50*t92)*s1l2
      taux16 = s1+((512.D0/3.D0*t121*t122-512.D0/3.D0*t121*t124)*s23r+51
     #2.D0/3.D0*t77*t66-512.D0/3.D0*t77*t70)*s1vr+(512.D0/3.D0*t3*t122-5
     #12.D0/3.D0*t3*t124)*s23r+512.D0/3.D0*t140*t66-512.D0/3.D0*t140*t70

CCC Coefficients

      tcoff41 = taux13
      tcoff42 = taux13
      tcoff43 = taux14
      tcoff44 = taux14
      tcoff45 = taux15
      tcoff46 = taux15
      tcoff47 = taux16
      tcoff48 = taux16

CCC Jet functions

      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet42 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet43 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,3.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet44 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,3.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s23r
      t7 = t5*s1vr
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet45 = Fjet(var,ndim,z,ux,t4,t6,-1-t4-t7+z,t7-t6-t12,t12,t7,-1+
     #z-t7+t12,-t4-t6-t12+z,1+t4+t6,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet46 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = 1-l3
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet47 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t
     #7+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s
     #4l2,6.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet48 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,6.D0,1.D0,1.D0,Tcoeff*tcoff48)

CCC Sum

      sq3num6e4 = tFjet41*tcoff41+tFjet42*tcoff42+tFjet43*tcoff43+tFjet4
     #4*tcoff44+tFjet45*tcoff45+tFjet46*tcoff46+tFjet47*tcoff47+tFjet48*
     #tcoff48
      RETURN
      END

