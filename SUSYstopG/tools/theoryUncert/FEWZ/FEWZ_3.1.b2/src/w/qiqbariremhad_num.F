CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbirmnum1e0(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t5 = -1+z
      t6 = t5**2
      t7 = t6*t5
      t8 = t3*z*t7
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t19 = t6*t2
      t23 = t12*s23r
      t27 = t10*s34r
      t28 = t27*t11
      t32 = s1l1**2
      t34 = t10*s2l1
      t37 = t6*z
      t38 = t2*t10
      t44 = t2*z
      t52 = l3*t7
      t53 = t52*z
      t54 = t11*s23r
      t58 = t6*t10
      t59 = t58*s2l1
      t62 = t58*s3l1
      t65 = 2*z
      t66 = t65+1
      t75 = z**2
      t76 = t2*t75
      t80 = t75*t7
      t86 = t1*t6
      t87 = t86*l3
      t88 = s1vr**2
      t94 = l3*t2
      t95 = t7*t10
      t97 = t94*t95*s34r
      t99 = z+1
      t101 = t5*t10
      t104 = 768*t97-256*t99*t1*t101
      t108 = s23r**2
      t109 = t12*t108
      t113 = t1*t5
      t119 = l3*t3
      t120 = t6**2
      t122 = s34r**2
      t131 = l3*t10
      t132 = t131*s2l1
      t133 = t86*t132
      t135 = t131*s3l1
      t136 = t86*t135
      t137 = 256*t136
      t140 = t86*l3*z*t10
      t159 = t10*s3l1
      t162 = 2+z
      t164 = z*t5
      t165 = t164*t10
      t172 = z*t1
      t173 = t172*t6
      t174 = t11*t108
      t175 = t131*t174
      t178 = t38*s34r
      t179 = t53*t178
      t187 = z*t10
      t196 = l3**2
      t197 = t196*t7
      t198 = t197*z
      t199 = t1*t10
      t200 = t88*s1vr
      t208 = s2l1**2
      t209 = t131*t208
      t215 = t196*t120
      t216 = t215*z
      t218 = 128*t216*t178
      t219 = t75*t1
      t222 = t219*t6*l3*t10
      t223 = 256*t222
      t234 = t87*t187*s3l1
      t249 = z*t196
      t250 = t120*t5
      t258 = t75*l3*t7*t178
      t260 = s3l1**2
      t267 = 2*t75
      t268 = 1+t267+t65
      t276 = t219*t6
      t285 = t75*t120
      t287 = t131*t122
      t296 = t6*t196
      t297 = t296*t10
      t301 = t199*s34r
      t302 = t197*t301
      t303 = l3*t5
      t304 = t303*t187
      t306 = -256*t302+256*t304
      t309 = t131*s34r
      t310 = t86*t309
      t312 = 256*t310-256*t187
      t317 = t12*t108*s23r
      t320 = t303*t10
      t325 = t38*t122
      t343 = t249*t58
      t384 = -t5
      t398 = t198*t301
      t403 = t216*t325
      t405 = t187*s34r
      t421 = t196*l3
      t422 = z*t7*t421
      t439 = 128*t1*t120*t421*t405
      t440 = t10*t260
      t459 = t75*t5
      t460 = t459*t131
      t478 = -t99*z
      t496 = t384*z
      t533 = t196*t10
      s3 = -64*t8*t12*s1vr*t13*s13r
      s5 = ((-256*t19*t12*s1vr+128*t19*t23+256*t3*t7*t28)*t32+((-256*t19
     #*t34-384*t37*t38)*t11*s1vr+128*t44*t6*t23+256*t8*t28)*s1l1+(64*t53
     #*t38*t54+(-128*t44*t59-128*t44*t62-64*t66*t2*t37*t10)*t11)*s1vr+64
     #*t76*t6*t23+128*t80*t3*t28)*t13
      s8 = (-256*t87*t12*t88+(256*t87*t23+t104*t11)*s1vr-256*t87*t109+(-
     #512*t97+256*t113*t10)*t11*s23r+(-512*t119*t120*t10*t122+512*t19*t2
     #7)*t11)*t32+((-512*t133+t137-384*t140)*t11*t88+((256*t133-t137+512
     #*t140)*t11*s23r+(t104*s2l1+(-256*t94*t95*s3l1+768*t52*t44*t10)*s34
     #r-256*t86*t159-256*t162*t1*t165)*t11)*s1vr-256*t173*t175+(-512*t17
     #9+256*t172*t101)*t11*s23r+(-512*t119*t120*t187*t122+512*t37*t178)*
     #t11)*s1l1-64*t198*t199*t11*t200+(128*t198*t199*t54+(-256*t86*t209+
     #(256*t136-256*t140)*s2l1+t218-t223)*t11)*t88
      s7 = s8+(-64*t198*t199*t174+(128*t173*t132-t218+128*t234+t223)*t11
     #*s23r+(-128*t172*t101*t208+(256*t179-128*t162*t1*t165)*s2l1-128*t2
     #49*t250*t3*t10*t122+384*t258-128*t172*t101*t260-128*t219*t101*s3l1
     #-64*t268*t1*t165)*t11)*s1vr-128*t276*t175+(-256*t258+128*t219*t101
     #)*t11*s23r+(-256*t285*t3*t287+256*t76*t58*s34r)*t11
      s8 = s13r
      s6 = s7*s8
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(128*t297*t54*t88+(t306*t11*s23r+t312*t11)*s1vr+128*t296*t
     #317+(256*t302-256*t320)*t11*t108+(256*t215*t325-512*t310+128*t10)*
     #t11*s23r+(-512*t94*t95*t122+256*t113*t27)*t11)*t32+(((256*t296*t34
     #+256*t343)*t11*s23r+(256*t303*t159-128*t304)*t11)*t88+((256*t296*t
     #159-128*t343)*t11*t108+(t306*s2l1+(256*t197*t199*s3l1-256*t197*t17
     #2*t10)*s34r+256*(z-2)*t5*t135+256*t99*z*t320)*t11*s23r+(t312*s2l1+
     #(-768*t136+256*t140)*s34r+256*t384*t10*s3l1-128*t66*z*t10)*t11)*s1
     #vr+128*t37*t196*t317+(256*t398-256*t304)*t11*t108+(256*t403-512*t8
     #7*t405+128*t187)*t11*s23r+(-512*t53*t325+256*t113*t405)*t11)*s1l1+
     #(64*t422*t23-64*t249*t58*t11)*t200
      taux1 = s1+(-128*t422*t109+(128*t296*t10*t208+128*t249*t59-t439+12
     #8*t296*t440+128*t249*t62+64*(2+3*z)*z*t297)*t11*s23r+(256*t303*t15
     #9*s2l1+128*t398-256*t303*t440-128*t460)*t11)*t88+(64*t422*t317+(t4
     #39-64*t66*z*t297)*t11*t108+(128*t164*t209-128*t478*t5*t132+128*z*t
     #250*t2*t421*t10*t122+128*t478*t196*t7*t1*t27+128*t164*t131*t260-12
     #8*t496*t5*t135+64*t268*z*t320)*t11*s23r+(-128*t187*t208+128*t478*t
     #34-128*t403+(-256*t234+128*t222)*s34r-128*t187*t260+128*t496*t159-
     #64*(t267+1)*z*t10)*t11)*s1vr+64*t75*t6*t196*t317+(128*t219*t7*t533
     #*s34r-128*t460)*t11*t108+(128*t285*t2*t533*t122-256*t276*t309+64*t
     #75*t10)*t11*s23r+(-256*t80*t2*t287+128*t459*t301)*t11
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t5 = -1+z
      t6 = t5**2
      t7 = t6*t5
      t8 = t3*z*t7
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t19 = t6*z
      t20 = t19*t2
      t24 = t12*s23r
      t27 = t10*s34r
      t28 = t27*t11
      t34 = l3*t7*z
      t35 = t2*t10
      t36 = t11*s23r
      t40 = t2*z
      t41 = t6*t10
      t42 = t41*s2l1
      t45 = t41*s3l1
      t48 = 2*z
      t58 = z**2
      t59 = t2*t58
      t63 = t58*t7
      t70 = t1*t6*l3
      t71 = z*t10
      t72 = s1vr**2
      t80 = t35*s34r
      t81 = t34*t80
      t83 = t58*t1
      t84 = t5*t10
      t85 = t83*t84
      t91 = -t1
      t92 = t91*z
      t93 = t92*t6
      t94 = l3*t10
      t95 = s23r**2
      t96 = t11*t95
      t97 = t94*t96
      t107 = t6**2
      t109 = s34r**2
      t118 = l3**2
      t120 = t118*t7*z
      t121 = t1*t10
      t122 = t72*s1vr
      t132 = t94*s2l1
      t133 = z*t1*t6*t132
      t136 = t118*t107*z
      t138 = 512*t136*t80
      t141 = t83*t6*l3*t10
      t142 = 1024*t141
      t151 = t94*s3l1
      t158 = 512*t85
      t161 = z*t118
      t162 = t107*t5
      t170 = t58*l3*t7*t80
      t177 = 2*t58
      t178 = -t48+1-t177
      t196 = t58*t107
      t198 = t94*t109
      t207 = t161*t6
      t210 = l3*t5
      t216 = t12*t95
      t219 = t121*s34r
      t220 = t120*t219
      t221 = t58*t5
      t222 = t221*t94
      t223 = t220-t222
      t227 = t71*s34r
      t228 = t70*t227
      t239 = t12*t95*s23r
      t247 = t35*t109
      t248 = t136*t247
      t265 = t118*l3
      t266 = z*t7*t265
      t280 = 512*t1*t107*t265*t227
      t287 = t6*t118*t10
      t307 = t58-z
      t308 = t307*t5
      t364 = t118*t10
      s2 = -256*t8*t12*s1vr*t13*s13r+((512*t20*t12*s1vr-512*t20*t24-1024
     #*t8*t28)*s1l1+(256*t34*t35*t36+(512*t40*t42-512*t40*t45-256*(1-t48
     #)*t2*t19*t10)*t11)*s1vr-256*t59*t6*t24-512*t63*t3*t28)*t13
      s3 = s2+((1536*t70*t71*t11*t72+(-1024*t70*t71*t36+(-3072*t81+1024*
     #t85)*t11)*s1vr-1024*t93*t97+(2048*t81+1024*t92*t84)*t11*s23r+(2048
     #*l3*t3*t107*t71*t109-2048*t19*t80)*t11)*s1l1+256*t120*t121*t11*t12
     #2+(-512*t120*t121*t36+(1024*t133-t138+t142)*t11)*t72+(256*t120*t12
     #1*t96+(-512*t133+t138-512*t93*t151-t142)*t11*s23r+((-1024*t81+t158
     #)*s2l1+512*t161*t162*t3*t10*t109-1536*t170+512*(t58-t48)*t1*t84*s3
     #l1-256*t178*t1*z*t5*t10)*t11)*s1vr-512*t91*t58*t6*t97+(1024*t170-t
     #158)*t11*s23r+(1024*t196*t3*t198-1024*t59*t41*s34r)*t11)*s13r
      s1 = s3+((-1024*t207*t24+512*t210*t71*t11)*t72+(512*t207*t216+1024
     #*t223*t11*s23r+(-1024*t228-256*(t48-4*t58)*t10)*t11)*s1vr-512*t207
     #*t239+(-1024*t220+1024*t210*t71)*t11*t95+(-1024*t248+2048*t228-512
     #*t71)*t11*s23r+(2048*t34*t247-1024*t1*t5*t227)*t11)*s1l1+(-256*t26
     #6*t24+256*t161*t41*t11)*t122
      taux2 = s1+(512*t266*t216+(-512*t161*t42+t280-512*t161*t45+256*(-2
     #-3*z)*z*t287)*t11*s23r-512*t223*t11)*t72+(-256*t266*t239+(-t280+25
     #6*(t48+1)*z*t287)*t11*t95+(-512*t308*t132-512*z*t162*t2*t265*t10*t
     #109+512*(z+1)*z*t7*t1*t118*t27-512*t308*t151+256*t178*z*t210*t10)*
     #t11*s23r+(512*z*t5*t10*s2l1+512*t248+(1024*z*t1*t6*t151-512*t141)*
     #s34r+512*t307*t10*s3l1-256*(-t177+1)*z*t10)*t11)*s1vr-256*t58*t6*t
     #118*t239+(-512*t83*t7*t364*s34r+512*t222)*t11*t95+(-512*t196*t2*t3
     #64*t109+1024*t83*t6*t94*s34r-256*t58*t10)*t11*s23r+(1024*t63*t2*t1
     #98-512*t221*t219)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = -taux2
      tcoff4 = -taux2
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux1
      tcoff10 = taux1
      tcoff11 = taux2
      tcoff12 = taux2
      tcoff13 = taux1
      tcoff14 = taux1
      tcoff15 = taux2
      tcoff16 = taux2
      tcoff17 = taux1
      tcoff18 = taux1
      tcoff19 = -taux2
      tcoff20 = -taux2
      tcoff21 = taux1
      tcoff22 = taux1
      tcoff23 = taux2
      tcoff24 = taux2
      tcoff25 = taux1
      tcoff26 = taux1
      tcoff27 = taux2
      tcoff28 = taux2
      tcoff29 = taux1
      tcoff30 = taux1
      tcoff31 = -taux2
      tcoff32 = -taux2

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet1 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet2 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet3 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,12.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet4 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet6 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,14.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet8 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,17.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet10 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet11 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet12 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet13 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet14 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet15 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet16 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet17 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet18 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet19 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet20 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet22 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet23 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet24 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff24)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet26 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet27 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet28 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet29 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet30 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet31 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet32 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      t17 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet29*tcoff29+tFjet28*tcof
     #f28+tFjet2*tcoff2+tFjet32*tcoff32+tFjet31*tcoff31+tFjet17*tcoff17+
     #tFjet18*tcoff18+tFjet5*tcoff5+tFjet4*tcoff4+tFjet3*tcoff3+tFjet19*
     #tcoff19+tFjet6*tcoff6+tFjet20*tcoff20+tFjet21*tcoff21
      t34 = tFjet7*tcoff7+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24*tcoff2
     #4+tFjet8*tcoff8+tFjet9*tcoff9+tFjet1*tcoff1+tFjet25*tcoff25+tFjet1
     #1*tcoff11+tFjet12*tcoff12+tFjet13*tcoff13+tFjet26*tcoff26+tFjet27*
     #tcoff27+tFjet14*tcoff14+tFjet30*tcoff30+tFjet10*tcoff10
      qiqbirmnum1e0 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum1e1(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t5 = -1+z
      t6 = t5**2
      t7 = t6*t5
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t21 = t11*s23r
      t22 = s1l1**2
      t26 = t2*z
      t34 = t2*t10
      t38 = t6*t10
      t42 = t6**2
      t45 = z*t10
      t46 = t45*s34r
      t52 = t6*z
      t59 = z**2
      t60 = t2*t59
      t62 = t12*s23r
      t65 = t59*t7
      t73 = t1*t6
      t74 = t73*l3
      t75 = s23r*s1vr
      t78 = s23r**2
      t79 = t12*t78
      t86 = l3*t10
      t91 = l3*z
      t98 = z*t1
      t100 = t11*t78
      t101 = t86*t100
      t109 = l3**2
      t111 = t109*t7*z
      t112 = t1*t10
      t113 = s1vr**2
      t114 = t113*s1vr
      t122 = t45*s2l1
      t123 = t74*t122
      t125 = t109*t42*z
      t126 = t34*s34r
      t127 = t125*t126
      t128 = t45*s3l1
      t129 = t74*t128
      t130 = 2*z
      t135 = t6*l3*t10
      t156 = t5*t10
      t157 = s2l1**2
      t161 = t5*z
      t163 = t161*t112*s3l1
      t168 = z*t109
      t169 = t42*t5
      t172 = s34r**2
      t178 = t59*l3*t7*t126
      t180 = s3l1**2
      t196 = t59*t1*t6
      t208 = t59*t42
      t210 = t86*t172
      t219 = t6*t109
      t220 = t219*t10
      t225 = t5*l3
      t231 = t12*t78*s23r
      t239 = t10*s2l1
      t240 = t219*t239
      t242 = t10*s3l1
      t243 = t219*t242
      t245 = t168*t38
      t257 = t225*t45
      t275 = t109*l3
      t276 = z*t7*t275
      t294 = 128*t1*t42*t275*t46
      t298 = t59*t6
      t299 = t109*t10
      t306 = t112*s34r
      t309 = t91*t156*s3l1
      t342 = 128*t59*t109*t7*t306
      t367 = t86*s34r
      t385 = t59*t5
      s3 = 128*t3*z*t7*t12*s1vr*t13*s13r
      s5 = (-128*t2*t6*t10*t21*t22-128*t26*t6*t12*s23r*s1l1+(-256*l3*t7*
     #z*t34*t21+(128*t26*t38*s2l1-128*t3*t42*l3*t46+128*t26*t38*s3l1+256
     #*t52*t34)*t11)*s1vr-128*t60*t6*t62-128*t65*t3*t10*s34r*t11)*t13
      s6 = ((256*t74*t12*t75+256*t74*t79-256*t1*t5*t62)*t22+((256*t73*t8
     #6*s2l1+256*t73*t86*s3l1+256*t73*t91*t10)*t11*t75+256*t98*t6*t101-2
     #56*t98*t5*t62)*s1l1+128*t111*t112*t11*t114+(-128*t111*t112*t21+(-1
     #28*t123-128*t127-128*t129+128*(t130-1)*z*t1*t135)*t11)*t113+(256*t
     #111*t112*t100+(-128*t123+256*t127-128*t129-128*(2+z)*z*t1*t135)*t1
     #1*s23r+(128*t98*t156*t157+(256*t163+256*t98*t156)*s2l1+128*t168*t1
     #69*t3*t10*t172-384*t178+128*t98*t156*t180+256*t163-64*(-2*t59-4+t1
     #30)*t1*t161*t10)*t11)*s1vr+256*t196*t101+(256*t178-256*t59*t1*t156
     #)*t11*s23r+(256*t208*t3*t210-256*t60*t38*s34r)*t11)*s13r
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-128*t220*t21*t113+(-256*t219*t79+256*t225*t62)*s1vr-128*
     #t219*t231+256*t225*t79-128*t62)*t22+((-256*t240-256*t243-128*t245)
     #*t11*s23r*t113+((-256*t240-256*t243-256*t245)*t11*t78+(256*t225*t2
     #39+256*t225*t242+256*t257)*t11*s23r)*s1vr-128*t52*t109*t231+256*t2
     #25*z*t79-128*t45*t21)*s1l1+(-128*t276*t62+128*t168*t38*t11)*t114
      taux1 = s1+(128*t276*t79+(-128*t219*t10*t157-256*t219*t242*s2l1+t2
     #94-128*t219*t10*t180-256*t298*t299)*t11*s23r+(-128*t225*t122-128*t
     #111*t306-128*t309+128*t52*t86)*t11)*t113+(-128*t276*t231+(-t294+64
     #*(2+2*z)*z*t220)*t11*t78+(-128*t161*t86*t157+(-256*t309-128*t257)*
     #s2l1-128*z*t169*t2*t275*t10*t172+t342-128*t161*t86*t180-128*t309-6
     #4*(2*t59+2)*z*t225*t10)*t11*s23r+(128*t45*t157+(256*t128+128*t45)*
     #s2l1+128*t125*t34*t172-128*t7*t1*z*t367+128*t45*t180+128*t128-64*(
     #-2+2*z-2*t59)*z*t10)*t11)*s1vr-128*t298*t109*t231+(-t342+256*t385*
     #t86)*t11*t78+(-128*t208*t2*t299*t172+256*t196*t367-128*t59*t10)*t1
     #1*s23r+(256*t65*t2*t210-128*t385*t306)*t11
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t16 = t3*z
      t17 = t16*t6
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t35 = t5**2
      t36 = t35*t3
      t37 = l3*t9
      t38 = t37*s34r
      t40 = 1024*t36*t38
      t41 = 6*z
      t42 = -1+t41
      t44 = t5*t9
      t51 = t2*z
      t56 = t9*s34r*t11
      t64 = l3*t2
      t69 = l3*t6
      t70 = t51*t9
      t76 = t69*z
      t77 = t2*t9
      t82 = 3*z
      t89 = t36*l3
      t90 = z*t9
      t91 = t90*s34r
      t94 = t44*s3l1
      t97 = -7+t41
      t99 = z*t5
      t107 = z**2
      t108 = t2*t107
      t112 = t107*t6
      t118 = t1*t6
      t119 = l3**2
      t120 = t118*t119
      t121 = t28*s1vr
      t125 = -t1
      t126 = t125*t6
      t130 = t2*t35
      t131 = t119*t9
      t132 = t131*s34r
      t133 = t130*t132
      t135 = t1*t5
      t144 = s23r**2
      t145 = t18*t144
      t149 = 8*z
      t152 = t135*t9
      t158 = t35*t4
      t159 = t3*t158
      t160 = s34r**2
      t161 = t131*t160
      t162 = t159*t161
      t164 = t77*s34r
      t165 = t76*t164
      t167 = 12*t107
      t170 = t4*t9
      t177 = 4*t1
      t178 = t177*z
      t179 = t178*t5
      t180 = t11*t144
      t181 = t37*t180
      t198 = t131*s2l1
      t201 = 3*t125
      t203 = t131*s3l1
      t206 = t119*t6
      t209 = t206*z*t1*t9
      t214 = t1*t9
      t229 = t130*t203
      t230 = t119*t35
      t231 = t230*t70
      t239 = t107*t1
      t241 = t5*l3*t9
      t248 = t206*z
      t259 = t230*z
      t262 = t37*s3l1
      t265 = 7*z
      t278 = 3*t107
      t307 = 6*t107
      t331 = t107*t35
      t333 = t37*t160
      t342 = t119*l3
      t343 = t6*t342
      t345 = t5*t119
      t352 = t125*t35
      t353 = t342*t9
      t354 = t353*s34r
      t355 = t352*t354
      t356 = 1024*t355
      t357 = z*t119
      t358 = t357*t44
      t363 = t126*t132
      t365 = 1-t82
      t374 = t18*t144*s23r
      t377 = 10*z
      t386 = t353*t160
      t389 = t214*s34r
      t391 = 3072*t248*t389
      t400 = t77*t160
      t401 = t230*t400
      t417 = l3*t4
      t425 = t201*z
      t442 = t9*s2l1
      t444 = 256*t343*t442
      t445 = t9*s3l1
      t447 = 256*t343*t445
      t448 = z*t6
      t449 = t448*t353
      t479 = 5*z
      t486 = t345*t9
      t493 = z+1
      t510 = t417*t9
      t520 = -z+2
      t521 = t520*t5
      t614 = t201*t107
      t618 = t107*t4
      s2 = (-256*t3*t6*t9*t11*s1vr*s1l1+768*t17*t18*s1vr)*t23*s13r+((-51
     #2*t27*t18*t28+(768*t27*t32+(t40-256*t42*t2*t44)*t11)*s1vr+2048*t51
     #*t5*t32+3072*t17*t56)*s1l1+(768*t26*t37*s2l1+256*t64*t6*t9*s3l1-25
     #6*t69*t70)*t11*t28+(-512*t76*t77*t11*s23r+((-t40-256*(-2+t82)*t2*t
     #44)*s2l1+512*t89*t91+1792*t51*t94-256*t97*t2*t99*t9)*t11)*s1vr+102
     #4*t108*t5*t32+1536*t112*t3*t56)*t23
      s4 = s2
      s7 = (-256*t120*t18*t121+(-1024*t126*t119*t32+(1024*t133-4608*t135
     #*l3*z*t9)*t11)*t28+(-768*t120*t145+(-2048*t133-512*(-t149+1)*l3*t1
     #52)*t11*s23r+(-1024*t162+9216*t165-256*(t167-1)*t1*t170)*t11)*s1vr
     #-1024*t179*t181+(-6144*t165+1024*t178*t170)*t11*s23r+(-6144*t89*t9
     #0*t160+6144*t99*t164)*t11)*s1l1+(-256*t118*t198-256*t201*t6*t203-1
     #024*t209)*t11*t121+((-512*t206*t214*s3l1+2304*t209)*t11*s23r+(256*
     #(2-11*z)*l3*t5*t214*s2l1+(-1536*t229+1536*t231)*s34r+768*t135*l3*t
     #90*s3l1-3072*t239*t241)*t11)*t28
      s6 = s7+(-1280*t248*t214*t180+((512*t133+1024*l3*t4*t152)*s2l1-204
     #8*t259*t164-512*t179*t262-512*(2-t265)*t1*z*t241)*t11*s23r+((512*t
     #162-512*(4-t265)*t6*t2*t38+512*(-t278+1+2*z)*t1*t170)*s2l1+(512*t1
     #59*t203-1536*t119*t158*t16*t9)*t160+(-512*t76*t77*s3l1+4608*t112*t
     #64*t9)*s34r+512*(t149-1-t278)*t1*t170*s3l1-256*(t307+t41-7)*t1*z*t
     #4*t9)*t11)*s1vr-512*t177*t107*t5*t181+(-3072*t107*l3*t6*t164+2048*
     #t239*t170)*t11*s23r+(-3072*t331*t3*t333+3072*t108*t44*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((256*t343*t32-256*t345*t18)*t121+(-512*t343*t145+(t356+35
     #84*t358)*t11*s23r+(-1024*t363+512*t365*t4*t37)*t11)*t28+(256*t343*
     #t374+(-t356+256*(1-t377)*t5*t131)*t11*t144+(1024*t2*t158*t386-t391
     #-256*(-t167+1-4*z)*t4*t37)*t11*s23r+(-1024*t401-1024*t365*l3*t5*t3
     #89-256*(t167+1-t41)*t9)*t11)*s1vr+2048*t99*t119*t374+(t391-4096*t4
     #17*t90)*t11*t144+(3072*t259*t400+2048*t425*t5*t38+2048*t90)*t11*s2
     #3r+(-6144*t76*t400+3072*t1*t4*t91)*t11)*s1l1+((t444-t447+1024*t449
     #)*t11*s23r+(-256*t345*t442+256*t345*t445-1024*t358)*t11)*t121
      s3 = s1
      s5 = ((-t444+t447-2048*t449)*t11*t144+((512*t355+1792*t358)*s2l1+(
     #-512*t352*t353*s3l1+512*t425*t35*t342*t9)*s34r+256*(t479-2)*t5*t20
     #3+256*(t377+7)*z*t486)*t11*s23r+((-512*t363+256*t493*t4*t37)*s2l1+
     #(1024*t126*t203+1536*t209)*s34r+256*(t82+1)*t4*t262-256*t42*z*t510
     #)*t11)*t28
      s6 = (1024*t448*t342*t374+(256*t521*t198-512*t425*t35*t354+256*t35
     #7*t94+256*(-t149-3)*z*t486)*t11*t144+((512*t520*t1*t6*t132-512*(-t
     #278+t82+1)*t4*t37)*s2l1+1536*z*t158*t2*t386+(512*t493*t1*t6*t203+5
     #12*(-t82-2)*z*t6*t1*t119*t9)*s34r-512*(-t278+t479-1)*t4*t262-256*(
     #-t149+4-t307)*z*t510)*t11*s23r+((512*t401-512*t521*t1*t38-256*t97*
     #z*t9)*s2l1+(512*t229-1536*t231)*t160+(3584*t125*z*t5*t262-512*t365
     #*t1*z*t241)*s34r+256*(9*z-t307-2)*t9*s3l1-256*(t307-3)*z*t9)*t11)*
     #s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+1024*t107*t5*t119*t374+(-512*t614*t6*t132-2048*t618*t37
     #)*t11*t144+(1536*t331*t2*t161+1024*t614*t5*t38+1024*t107*t9)*t11*s
     #23r+(-3072*t112*t2*t333+1536*t618*t389)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = -taux2
      tcoff4 = -taux2
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux1
      tcoff10 = taux1
      tcoff11 = taux2
      tcoff12 = taux2
      tcoff13 = taux1
      tcoff14 = taux1
      tcoff15 = taux2
      tcoff16 = taux2
      tcoff17 = taux1
      tcoff18 = taux1
      tcoff19 = -taux2
      tcoff20 = -taux2
      tcoff21 = taux1
      tcoff22 = taux1
      tcoff23 = taux2
      tcoff24 = taux2
      tcoff25 = taux1
      tcoff26 = taux1
      tcoff27 = taux2
      tcoff28 = taux2
      tcoff29 = taux1
      tcoff30 = taux1
      tcoff31 = -taux2
      tcoff32 = -taux2

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet1 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet2 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet3 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,12.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet4 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet6 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,14.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet8 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,17.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet10 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet11 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet12 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet13 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet14 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet15 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet16 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet17 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet18 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet19 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet20 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet22 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet23 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet24 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff24)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet26 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet27 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet28 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet29 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet30 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet31 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet32 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      t17 = tFjet21*tcoff21+tFjet20*tcoff20+tFjet17*tcoff17+tFjet18*tcof
     #f18+tFjet1*tcoff1+tFjet23*tcoff23+tFjet24*tcoff24+tFjet19*tcoff19+
     #tFjet32*tcoff32+tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+tF
     #jet22*tcoff22+tFjet2*tcoff2+tFjet3*tcoff3+tFjet28*tcoff28
      t34 = tFjet4*tcoff4+tFjet5*tcoff5+tFjet29*tcoff29+tFjet30*tcoff30+
     #tFjet31*tcoff31+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*t
     #coff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet14*tcoff14+tFjet12*tcof
     #f12+tFjet13*tcoff13+tFjet15*tcoff15+tFjet16*tcoff16
      qiqbirmnum1e1 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum1e2(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = l3-1
      t2 = t1**2
      t5 = -1+z
      t6 = t5**2
      t7 = t6*t5
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t20 = t2*t7*z
      t21 = l3*t10
      t22 = s1vr**2
      t27 = t11*s23r
      t28 = t21*t27
      t30 = t6*z
      t37 = z**2
      t40 = t12*s23r
      t45 = l3**2
      t47 = z*t45*t7
      t48 = t1*t10
      t49 = t22*s1vr
      t56 = t1*t6
      t59 = z*t10*t11
      t64 = s23r**2
      t65 = t11*t64
      t81 = t37*t1
      t93 = z*t7*t45*l3
      t95 = t45*t10
      t101 = t12*t64
      t116 = t12*t64*s23r
      t120 = (2*z-3)*z
      taux1 = -64*t2*t1*z*t7*t12*s1vr*t13*s13r+(128*t20*t21*t11*t22+(192
     #*t20*t28-192*t30*t2*t10*t11)*s1vr+64*t2*t37*t6*t40)*t13+(-64*t47*t
     #48*t11*t49+(-256*t47*t48*t27+256*t56*l3*t59)*t22+(-192*t47*t48*t65
     #-128*(z-3)*z*t56*t28-192*t1*z*t5*t10*t11)*s1vr-128*t81*t6*t21*t65+
     #128*t81*t5*t40)*s13r+(64*t93*t40-64*t30*t95*t11)*t49+(128*t93*t101
     #+64*(-4+z)*z*t6*t95*t27+128*l3*t5*t59)*t22+(64*t93*t116+64*t120*t6
     #*t95*t65-64*t120*t5*t28-64*t59)*s1vr+64*t37*t6*t45*t116-128*t37*t5
     #*l3*t101+64*t37*t10*t27
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t16 = t3*z
      t17 = t16*t6
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = t18*s23r
      t31 = t5**2
      t33 = l3*t9
      t34 = t33*s34r
      t36 = 2048*t31*t3*t34
      t37 = -1-z
      t40 = t5*t9
      t47 = t5*z
      t52 = t9*s34r*t11
      t60 = l3*t2
      t67 = s1vr**2
      t69 = t2*z
      t81 = z*t31
      t85 = t40*s3l1
      t88 = 4*z
      t98 = z**2
      t99 = t2*t98
      t103 = t98*t6
      t109 = -t1
      t110 = 4*t109
      t111 = t110*t6
      t112 = l3**2
      t113 = t111*t112
      t114 = t67*s1vr
      t118 = t2*t31
      t119 = t112*t9
      t120 = t119*s34r
      t121 = t118*t120
      t123 = 6*z
      t126 = t5*t1
      t127 = t126*t9
      t133 = s23r**2
      t138 = t123+4
      t145 = t31*t4
      t146 = t3*t145
      t147 = s34r**2
      t148 = t119*t147
      t149 = t146*t148
      t151 = l3*t6
      t152 = t151*z
      t153 = t2*t9
      t154 = t153*s34r
      t155 = t152*t154
      t157 = 8*t98
      t160 = t4*t9
      t167 = 5*t109
      t168 = t167*z
      t169 = t168*t5
      t170 = t11*t133
      t171 = t33*t170
      t182 = z*t9
      t191 = t119*s2l1
      t196 = t119*s3l1
      t199 = t112*t6
      t200 = z*t1
      t202 = t199*t200*t9
      t207 = t1*t6
      t219 = t1*t9
      t223 = t118*t196
      t225 = t112*t31
      t227 = t225*t69*t9
      t231 = -t138
      t234 = t33*s3l1
      t237 = 2*z
      t242 = t5*l3*t9
      t249 = z*t112
      t250 = t249*t6
      t266 = 5*z
      t275 = 3*z
      t293 = t182*s3l1
      t306 = 4*t98
      t321 = t103*t2
      t324 = t98*t1
      t330 = t98*t31
      t332 = t33*t147
      t341 = t112*l3
      t342 = t6*t341
      t344 = t5*t112
      t349 = t1*t31
      t350 = t341*t9
      t351 = t350*s34r
      t352 = t349*t351
      t353 = 2048*t352
      t354 = 7*z
      t362 = t219*s34r
      t363 = t199*t362
      t374 = t18*t133*s23r
      t377 = 8*z
      t386 = t350*t147
      t390 = 2048*t250*t362
      t398 = t153*t147
      t399 = t225*t398
      t413 = l3*t4
      t422 = t126*l3
      t423 = t182*s34r
      t439 = t9*s2l1
      t442 = t342*t182
      t448 = 1024*t344*t439
      t484 = t344*t9
      t508 = t413*t9
      t531 = z-2
      t536 = 2*t98
      t609 = t98*t4
      s2 = (1024*t3*t6*t9*t11*s1vr*s1l1-256*t17*t18*s1vr)*t23*s13r+(((-3
     #072*t27*t28+(-t36-1024*t37*t2*t40)*t11)*s1vr-2560*t47*t2*t28-2048*
     #t17*t52)*s1l1+(-2048*t26*t33*s2l1-1024*t60*t6*t9*s3l1)*t11*t67+(-7
     #68*t69*t6*t33*t11*s23r+((t36-512*(2+z)*t2*t40)*s2l1-1024*t81*t3*t3
     #4-1536*t69*t85-256*(7-t88)*t2*t47*t9)*t11)*s1vr-1280*t99*t5*t28-10
     #24*t103*t3*t52)*t23
      s4 = s2
      s7 = (-256*t113*t18*t114+(-2048*t121+512*(-4+t123)*l3*t127)*t11*t6
     #7+(-768*t113*t18*t133+(4096*t121-512*t138*l3*t127)*t11*s23r+(2048*
     #t149-6144*t155-256*(-4-t157)*t1*t160)*t11)*s1vr-1024*t169*t171+(40
     #96*t155+1024*t168*t160)*t11*s23r+(4096*l3*t3*t31*t182*t147-4096*t4
     #7*t154)*t11)*s1l1+(-256*t111*t191+256*t110*t6*t196+1280*t202)*t11*
     #t114+((1024*t207*t191-512*t111*t196-1536*t202)*t11*s23r+(256*(-8+t
     #123)*l3*t5*t219*s2l1+(3072*t223-1024*t227)*s34r+256*t231*t1*t5*t23
     #4+1024*(t237-1)*t1*z*t242)*t11)*t67
      s6 = s7+(2304*t250*t219*t170+((-1024*t121-512*(-z-4)*l3*t127)*s2l1
     #+2048*t81*t2*t120-512*t169*t234-512*(t266-1)*t1*z*t242)*t11*s23r+(
     #(-1024*t149-1024*(-4+t275)*t6*t2*t34+1024*(t98-1-t237)*t1*t160)*s2
     #l1+(-1024*t146*t196+1024*t112*t145*t16*t9)*t147+(1024*t27*t293-307
     #2*t103*t60*t9)*s34r+1024*(t98-t88+1)*t1*t160*s3l1-256*(7-t306-t88)
     #*t1*z*t4*t9)*t11)*s1vr-512*t167*t98*t5*t171+(2048*t321*t34-2560*t3
     #24*t160)*t11*s23r+(2048*t330*t3*t332-2048*t99*t40*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((-1024*t342*t28+1024*t344*t18)*t114+((t353+512*(4-t354)*t
     #5*t119)*t11*s23r+(-2048*t363+512*(-4+t237)*t4*t33)*t11)*t67+(-1024
     #*t342*t374+(-t353+256*(t377+4)*t5*t119)*t11*t133+(-2048*t2*t145*t3
     #86+t390-256*(t88+4+t157)*t4*t33)*t11*s23r+(2048*t399-2048*t151*t36
     #2-256*(-4-t157+t88)*t9)*t11)*s1vr-2560*t47*t112*t374+(-t390+5120*t
     #413*t182)*t11*t133+(-2048*t225*z*t398+4096*t422*t423-2560*t182)*t1
     #1*s23r+(4096*t152*t398-2048*t1*t4*t423)*t11)*s1l1+((-1024*t342*t43
     #9-1280*t442)*t11*s23r+(t448+1280*t249*t40)*t11)*t114
      s3 = s1
      s5 = ((-1024*t342*t9*s3l1+1536*t442)*t11*t133+((1024*t352+256*(4-t
     #377)*t5*t119)*s2l1+(-1024*t349*t350*s3l1+1024*t200*t31*t341*t9)*s3
     #4r+256*(8-t88)*t5*t196+256*(-2-9*z)*z*t484)*t11*s23r+((-1024*t363+
     #256*(-4-t237)*t4*t33)*s2l1+(2048*t199*t219*s3l1-1024*t202)*s34r+25
     #6*t231*t4*t234+1024*z*t4*t508)*t11)*t67
      s6 = (-1280*z*t6*t341*t374+(-t448-1024*t200*t31*t351-1024*t249*t85
     #+256*(t123+5)*z*t484)*t11*t133+((1024*t531*t1*t6*t120-512*(-t275+t
     #536-2)*t4*t33)*s2l1-1024*z*t145*t2*t386+(1024*t37*t1*t6*t196+1024*
     #t103*t1*t112*t9)*s34r-512*(-t354+t536+2)*t4*t234-256*(t123-1+t306)
     #*z*t508)*t11*s23r+((-1024*t399-1024*t531*t5*t1*t34+512*(t237-3)*z*
     #t9)*s2l1+(-1024*t223+1024*t227)*t147+(3072*t422*t293-1024*t207*z*l
     #3*t9)*s34r+512*(2-t266+t536)*t9*s3l1-256*(-t306+1)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2-1280*t98*t5*t112*t374+(-1024*t324*t6*t120+2560*t609*t33
     #)*t11*t133+(-1024*t330*t2*t148+2048*t324*t5*t34-1280*t98*t9)*t11*s
     #23r+(2048*t321*t332-1024*t609*t362)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = -taux2
      tcoff4 = -taux2
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux1
      tcoff10 = taux1
      tcoff11 = taux2
      tcoff12 = taux2
      tcoff13 = taux1
      tcoff14 = taux1
      tcoff15 = taux2
      tcoff16 = taux2
      tcoff17 = taux1
      tcoff18 = taux1
      tcoff19 = -taux2
      tcoff20 = -taux2
      tcoff21 = taux1
      tcoff22 = taux1
      tcoff23 = taux2
      tcoff24 = taux2
      tcoff25 = taux1
      tcoff26 = taux1
      tcoff27 = taux2
      tcoff28 = taux2
      tcoff29 = taux1
      tcoff30 = taux1
      tcoff31 = -taux2
      tcoff32 = -taux2

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet1 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet2 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet3 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,12.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet4 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet6 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,14.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet8 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,17.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet10 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet11 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet12 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet13 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet14 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet15 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet16 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet17 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet18 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet19 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet20 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet22 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet23 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet24 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff24)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet26 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet27 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet28 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet29 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet30 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet31 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet32 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      t17 = tFjet27*tcoff27+tFjet20*tcoff20+tFjet26*tcoff26+tFjet3*tcoff
     #3+tFjet2*tcoff2+tFjet19*tcoff19+tFjet1*tcoff1+tFjet18*tcoff18+tFje
     #t32*tcoff32+tFjet31*tcoff31+tFjet17*tcoff17+tFjet16*tcoff16+tFjet8
     #*tcoff8+tFjet6*tcoff6+tFjet30*tcoff30+tFjet7*tcoff7
      t34 = tFjet29*tcoff29+tFjet23*tcoff23+tFjet24*tcoff24+tFjet9*tcoff
     #9+tFjet25*tcoff25+tFjet11*tcoff11+tFjet10*tcoff10+tFjet12*tcoff12+
     #tFjet15*tcoff15+tFjet22*tcoff22+tFjet14*tcoff14+tFjet13*tcoff13+tF
     #jet5*tcoff5+tFjet28*tcoff28+tFjet21*tcoff21+tFjet4*tcoff4
      qiqbirmnum1e2 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum1e3(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux1 = 0
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t39 = t2*z
      t45 = l3*t9
      t46 = t45*s2l1
      t60 = t11*s23r
      t63 = t5*z
      t70 = z**2
      t77 = t1*t6
      t78 = l3**2
      t79 = t77*t78
      t80 = t28*s1vr
      t85 = l3*t5
      t86 = t1*t9
      t92 = s23r**2
      t93 = t18*t92
      t96 = 2*z
      t108 = t1*z
      t110 = t11*t92
      t111 = t45*t110
      t119 = t78*t9
      t120 = t119*s2l1
      t121 = t77*t120
      t123 = t78*t6
      t125 = t123*t86*s3l1
      t128 = t123*t108*t9
      t137 = t1*t5
      t139 = t45*s3l1
      t149 = z*t78
      t154 = t137*l3
      t155 = z*t9
      t156 = t155*s2l1
      t175 = t1*t70
      t185 = t6*t78*l3
      t187 = t5*t78
      t194 = -4+z
      t199 = t4*l3
      t205 = t18*t92*s23r
      t209 = -6+4*z
      t231 = t9*s2l1
      t232 = t185*t231
      t234 = t9*s3l1
      t235 = t185*t234
      t237 = t185*t155
      t246 = t5*t9
      t258 = (t96-4)*t5
      t263 = t187*t9
      t286 = (-3+t96)*z
      s2 = (-1024*t3*t6*t9*t11*s1vr*s1l1-512*t3*z*t6*t18*s1vr)*t23*s13r+
     #((2048*t27*t18*t28+(3072*t27*t32-3072*t2*t5*t18)*s1vr+1024*t39*t5*
     #t32)*s1l1+(1024*t26*t46+1024*l3*t2*t6*t9*s3l1+1024*l3*t6*t39*t9)*t
     #11*t28+(1536*t39*t6*t45*t60-1536*t63*t2*t9*t11)*s1vr+512*t2*t70*t5
     #*t32)*t23
      s3 = s2+((-1024*t79*t18*t80+(-4096*t79*t32+4096*t85*t86*t11)*t28+(
     #-3072*t79*t93-1024*(-6+t96)*l3*t5*t86*t60-3072*t1*t4*t18)*s1vr-204
     #8*t108*t5*t111+2048*t108*t4*t32)*s1l1+(-1024*t121-1024*t125-512*t1
     #28)*t11*t80+((-2048*t121-2048*t125-2048*t128)*t11*s23r+(2048*t137*
     #t46+2048*t137*t139+2048*t137*l3*z*t9)*t11)*t28+(-1536*t149*t6*t86*
     #t110+(-1024*t154*t156-1024*t154*t155*s3l1-1024*(-3+z)*t1*z*t85*t9)
     #*t11*s23r-1536*t108*t4*t9*t11)*s1vr-1024*t175*t5*t111+1024*t175*t4
     #*t32)*s13r
      s1 = s3+((1024*t185*t32-1024*t187*t18)*t80+(2048*t185*t93+1024*t19
     #4*t5*t78*t32+2048*t199*t18)*t28+(1024*t185*t205+512*t209*t5*t78*t9
     #3-512*t209*t4*l3*t32-1024*t18)*s1vr+1024*t63*t78*t205-2048*t199*z*
     #t93+1024*t155*t60)*s1l1+((1024*t232+1024*t235+512*t237)*t11*s23r+(
     #-1024*t187*t231-1024*t187*t234-512*t149*t246)*t11)*t80
      taux2 = s1+((1024*t232+1024*t235+1024*t237)*t11*t92+(512*t258*t120
     #+512*t258*t119*s3l1+512*t194*z*t263)*t11*s23r+(1024*t199*t231+1024
     #*t199*t234+1024*t199*t155)*t11)*t28+(512*t185*z*t205+(1024*t63*t12
     #0+1024*t149*t246*s3l1+512*t286*t263)*t11*t92+(-1024*t199*t156-1024
     #*z*t4*t139-512*t286*t199*t9)*t11*s23r-512*t155*t11)*s1vr+512*t70*t
     #5*t78*t205-1024*t70*t4*l3*t93+512*t70*t9*t60

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = -taux2
      tcoff4 = -taux2
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux1
      tcoff10 = taux1
      tcoff11 = taux2
      tcoff12 = taux2
      tcoff13 = taux1
      tcoff14 = taux1
      tcoff15 = taux2
      tcoff16 = taux2
      tcoff17 = taux1
      tcoff18 = taux1
      tcoff19 = -taux2
      tcoff20 = -taux2
      tcoff21 = taux1
      tcoff22 = taux1
      tcoff23 = taux2
      tcoff24 = taux2
      tcoff25 = taux1
      tcoff26 = taux1
      tcoff27 = taux2
      tcoff28 = taux2
      tcoff29 = taux1
      tcoff30 = taux1
      tcoff31 = -taux2
      tcoff32 = -taux2

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet1 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet2 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet3 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,12.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet4 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet6 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,14.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet8 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,17.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet10 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet11 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet12 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet13 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet14 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet15 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet16 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet17 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet18 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet19 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet20 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet22 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet23 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet24 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff24)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet26 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet27 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet28 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet29 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet30 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet31 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet32 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      t17 = tFjet15*tcoff15+tFjet16*tcoff16+tFjet17*tcoff17+tFjet20*tcof
     #f20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet12*tcoff12+tFjet13*tcoff1
     #3+tFjet14*tcoff14+tFjet32*tcoff32+tFjet4*tcoff4+tFjet5*tcoff5+tFje
     #t6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet9*tcoff9
      t34 = tFjet18*tcoff18+tFjet19*tcoff19+tFjet10*tcoff10+tFjet11*tcof
     #f11+tFjet29*tcoff29+tFjet30*tcoff30+tFjet31*tcoff31+tFjet23*tcoff2
     #3+tFjet24*tcoff24+tFjet25*tcoff25+tFjet26*tcoff26+tFjet27*tcoff27+
     #tFjet28*tcoff28+tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3
      qiqbirmnum1e3 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum1e4(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux1 = 0
      taux2 = 0

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = -taux2
      tcoff4 = -taux2
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux1
      tcoff10 = taux1
      tcoff11 = taux2
      tcoff12 = taux2
      tcoff13 = taux1
      tcoff14 = taux1
      tcoff15 = taux2
      tcoff16 = taux2
      tcoff17 = taux1
      tcoff18 = taux1
      tcoff19 = -taux2
      tcoff20 = -taux2
      tcoff21 = taux1
      tcoff22 = taux1
      tcoff23 = taux2
      tcoff24 = taux2
      tcoff25 = taux1
      tcoff26 = taux1
      tcoff27 = taux2
      tcoff28 = taux2
      tcoff29 = taux1
      tcoff30 = taux1
      tcoff31 = -taux2
      tcoff32 = -taux2

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet1 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet2 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet3 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t3
     #+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s3
     #l2,12.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet4 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-1
     #+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet6 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,14.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet8 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet9 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,17.D0,1.D0,1.D0,Tcoeff*tcoff9)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet10 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet11 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet12 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff12)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet13 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff13)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet14 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet15 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet16 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet17 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet18 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet19 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet20 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet21 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet22 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet23 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet24 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff24)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet25 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff25)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet26 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff26)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet27 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff27)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet28 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff28)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet29 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff29)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet30 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff30)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet31 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff31)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet32 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff32)

CCC Sum

      t17 = tFjet10*tcoff10+tFjet11*tcoff11+tFjet18*tcoff18+tFjet12*tcof
     #f12+tFjet19*tcoff19+tFjet1*tcoff1+tFjet20*tcoff20+tFjet2*tcoff2+tF
     #jet21*tcoff21+tFjet22*tcoff22+tFjet3*tcoff3+tFjet23*tcoff23+tFjet2
     #4*tcoff24+tFjet4*tcoff4+tFjet25*tcoff25+tFjet31*tcoff31
      t34 = tFjet32*tcoff32+tFjet9*tcoff9+tFjet26*tcoff26+tFjet13*tcoff1
     #3+tFjet5*tcoff5+tFjet27*tcoff27+tFjet6*tcoff6+tFjet7*tcoff7+tFjet2
     #8*tcoff28+tFjet14*tcoff14+tFjet15*tcoff15+tFjet16*tcoff16+tFjet29*
     #tcoff29+tFjet8*tcoff8+tFjet17*tcoff17+tFjet30*tcoff30
      qiqbirmnum1e4 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum2e0(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s1l1**2
      t3 = z**2
      t8 = r*t
      t9 = -1+t8
      t10 = t9**2
      t12 = r**2
      t13 = 1/t12
      t15 = t12*r
      t16 = t15*t
      t17 = l1**2
      t18 = dsqrt(l2)
      t19 = t17*t18
      t21 = t**2
      t22 = t12*t21
      t26 = l1*t18
      t27 = t12**2
      t33 = t18*t15
      t35 = t18*t12
      t37 = t18*r
      t44 = (1+t12)**2
      t47 = t12*l1
      t48 = t47*t
      t50 = l1*r
      t52 = l1*t
      t56 = dsqrt(-l1*(l1-1)*(1-l2)*t44*r*(t12*t+t48+r+r*t21-t50-t50*t21
     #+t+t52))
      t60 = t16*t19-t19*t22-t19*t12+t8*t19+t26*t27+2*t26*t22+4*t26*t12+t
     #26-t33*t-t35*t21-t35-t37*t+2*t56-4*t56*l4
      t61 = t60**2
      t62 = r+t
      t63 = t62**2
      t66 = 1/t44
      t78 = l1*l2
      t82 = l2*t12
      t90 = l2*t17
      t95 = t56*t18
      t99 = -l1-2*t47-l1*t27+t78+4*t78*t12+t78*t27-t82-l2*t15*t-l2*r*t-t
     #82*t21+2*t78*t22-t90*t12+t90*t16+t90*t8-t90*t22-4*t95*l4+2*t95
      t105 = dsqrt((r-t)*s1l1*s2l1*t99*t9/r*t66)
      t109 = t18*s1l1
      t116 = s1l1*l1
      t121 = -t33*s1l1-2*t15*t105*c3l1+t109*t48+t37*s1l1*t21-2*t105*c3l1
     #*r-t37*t116-t37*t116*t21+t109*t52
      t130 = s1vr**2
      t131 = t121**2
      t137 = s4vr**2
      taux3 = ((-256*t1-128*t3-256*z*s1l1)*t10*t13*t61/t63*t66+256*s1vr*
     #(2*s1l1+z)*t121*t66*t13*t9*t60/t62-256*t130*t131*t66*t13)/t137
      t3 = z**2
      t5 = r*t
      t6 = -1+t5
      t7 = t6**2
      t9 = r**2
      t10 = 1/t9
      t12 = t9*r
      t13 = t12*t
      t14 = l1**2
      t15 = dsqrt(l2)
      t16 = t14*t15
      t18 = t**2
      t19 = t9*t18
      t23 = l1*t15
      t24 = t9**2
      t30 = t15*t12
      t32 = t15*t9
      t34 = t15*r
      t41 = (1+t9)**2
      t44 = t9*l1
      t45 = t44*t
      t47 = l1*r
      t49 = l1*t
      t53 = dsqrt(-l1*(l1-1)*(1-l2)*t41*r*(t9*t+t45+r+t18*r-t47-t47*t18+
     #t+t49))
      t57 = t13*t16-t16*t19-t16*t9+t5*t16+t23*t24+2*t23*t19+4*t23*t9+t23
     #-t30*t-t32*t18-t32-t34*t+2*t53-4*t53*l4
      t58 = t57**2
      t59 = r+t
      t60 = t59**2
      t63 = 1/t41
      t78 = l1*l2
      t82 = l2*t9
      t90 = l2*t14
      t95 = t53*t15
      t99 = -l1-2*t44-l1*t24+t78+4*t78*t9+t78*t24-t82-l2*t12*t-l2*r*t-t8
     #2*t18+2*t78*t19-t90*t9+t90*t13+t90*t5-t90*t19-4*t95*l4+2*t95
      t105 = dsqrt((r-t)*s1l1*s2l1*t99*t6/r*t63)
      t109 = t15*s1l1
      t116 = s1l1*l1
      t127 = s4vr**2
      taux4 = (-512*(-2*z*s1l1-t3)*t7*t10*t58/t60*t63-1024*s1vr*t6*t10*t
     #57/t59*t63*z*(-t30*s1l1-2*t12*t105*c3l1+t109*t45+t34*s1l1*t18-2*t1
     #05*c3l1*r-t34*t116-t34*t116*t18+t109*t49))/t127

CCC Coefficients

      tcoff33 = taux3
      tcoff34 = taux3
      tcoff35 = -taux4
      tcoff36 = -taux4
      tcoff37 = taux3
      tcoff38 = taux3
      tcoff39 = taux4
      tcoff40 = taux4
      tcoff41 = taux3
      tcoff42 = taux3
      tcoff43 = taux4
      tcoff44 = taux4
      tcoff45 = taux3
      tcoff46 = taux3
      tcoff47 = taux4
      tcoff48 = taux4
      tcoff49 = taux3
      tcoff50 = taux3
      tcoff51 = -taux4
      tcoff52 = -taux4
      tcoff53 = taux3
      tcoff54 = taux3
      tcoff55 = taux4
      tcoff56 = taux4
      tcoff57 = taux3
      tcoff58 = taux3
      tcoff59 = taux4
      tcoff60 = taux4
      tcoff61 = taux3
      tcoff62 = taux3
      tcoff63 = -taux4
      tcoff64 = -taux4

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet33 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet34 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff34)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet36 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet38 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet39 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet40 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff40)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,17.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet42 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet43 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet44 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet45 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet46 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet47 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet48 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff48)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet49 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff49)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet50 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff50)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet51 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff51)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet52 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff52)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet53 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff53)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet54 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff54)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet55 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff55)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet56 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff56)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet57 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff57)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet58 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff58)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet59 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff59)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet60 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff60)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet61 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff61)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet62 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff62)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet63 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff63)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet64 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff64)

CCC Sum

      t17 = tFjet49*tcoff49+tFjet34*tcoff34+tFjet35*tcoff35+tFjet36*tcof
     #f36+tFjet37*tcoff37+tFjet38*tcoff38+tFjet50*tcoff50+tFjet39*tcoff3
     #9+tFjet40*tcoff40+tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+
     #tFjet48*tcoff48+tFjet59*tcoff59+tFjet60*tcoff60+tFjet61*tcoff61
      t34 = tFjet62*tcoff62+tFjet63*tcoff63+tFjet45*tcoff45+tFjet64*tcof
     #f64+tFjet46*tcoff46+tFjet47*tcoff47+tFjet33*tcoff33+tFjet41*tcoff4
     #1+tFjet54*tcoff54+tFjet42*tcoff42+tFjet55*tcoff55+tFjet43*tcoff43+
     #tFjet56*tcoff56+tFjet44*tcoff44+tFjet57*tcoff57+tFjet58*tcoff58
      qiqbirmnum2e0 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum2e1(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = z**2
      t2 = r*t
      t4 = (-1+t2)**2
      t6 = r**2
      t9 = t6*r
      t11 = l1**2
      t12 = dsqrt(l2)
      t13 = t11*t12
      t15 = t**2
      t16 = t6*t15
      t20 = l1*t12
      t21 = t6**2
      t29 = t12*t6
      t38 = (1+t6)**2
      t44 = l1*r
      t50 = dsqrt(-l1*(l1-1)*(1-l2)*t38*r*(t6*t+t6*l1*t+r+r*t15-t44-t44*
     #t15+t+l1*t))
      t54 = t9*t*t13-t13*t16-t13*t6+t2*t13+t20*t21+2*t20*t16+4*t20*t6+t2
     #0-t12*t9*t-t29*t15-t29-t12*r*t+2*t50-4*t50*l4
      t55 = t54**2
      t57 = (r+t)**2
      t61 = s4vr**2
      taux3 = 128*t1*t4/t6*t55/t57/t38/t61
      t3 = l3*s1vr
      t4 = s2l1*s13
      t11 = z**2
      t14 = r*t
      t15 = -1+t14
      t16 = t15**2
      t18 = r**2
      t19 = 1/t18
      t21 = t18*r
      t22 = t21*t
      t23 = l1**2
      t24 = dsqrt(l2)
      t25 = t23*t24
      t27 = t**2
      t28 = t18*t27
      t32 = l1*t24
      t33 = t18**2
      t39 = t24*t21
      t41 = t24*t18
      t43 = t24*r
      t50 = (1+t18)**2
      t53 = t18*l1
      t54 = t53*t
      t56 = l1*r
      t58 = l1*t
      t62 = dsqrt(-l1*(l1-1)*(1-l2)*t50*r*(t18*t+t54+r+r*t27-t56-t56*t27
     #+t+t58))
      t66 = t22*t25-t25*t28-t25*t18+t14*t25+t32*t33+2*t32*t28+4*t32*t18+
     #t32-t39*t-t41*t27-t41-t43*t+2*t62-4*t62*l4
      t67 = t66**2
      t68 = r+t
      t69 = t68**2
      t72 = 1/t50
      t76 = s1vr*t15*t19
      t78 = t66/t68
      t79 = s13*l3
      t101 = l1*l2
      t105 = l2*t18
      t113 = l2*t23
      t118 = t62*t24
      t122 = -l1-2*t53-l1*t33+t101+4*t101*t18+t101*t33-t105-l2*t21*t-l2*
     #r*t-t105*t27+2*t101*t28-t113*t18+t113*t22+t113*t14-t113*t28-4*t118
     #*l4+2*t118
      t128 = dsqrt((r-t)*s1l1*s2l1*t122*t15/r*t72)
      t132 = t24*s1l1
      t139 = s1l1*l1
      t148 = s4vr**2
      taux4 = (-512*(6*z*s1l1+t3*t4+t3*s2l1-t3*s2l1*z-t3*t4*z+3*t11)*t16
     #*t19*t67/t69*t72-512*(-t76*t78*t72*(-t79*s1vr-s13+s34*s13-t3+t79*s
     #1vr*z+s34+2*z-1+t3*z)-4*t76*t78*t72*z)*(-t39*s1l1-2*t21*t128*c3l1+
     #t132*t54+t43*s1l1*t27-2*t128*c3l1*r-t43*t139-t43*t139*t27+t132*t58
     #))/t148

CCC Coefficients

      tcoff33 = taux3
      tcoff34 = taux3
      tcoff35 = -taux4
      tcoff36 = -taux4
      tcoff37 = taux3
      tcoff38 = taux3
      tcoff39 = taux4
      tcoff40 = taux4
      tcoff41 = taux3
      tcoff42 = taux3
      tcoff43 = taux4
      tcoff44 = taux4
      tcoff45 = taux3
      tcoff46 = taux3
      tcoff47 = taux4
      tcoff48 = taux4
      tcoff49 = taux3
      tcoff50 = taux3
      tcoff51 = -taux4
      tcoff52 = -taux4
      tcoff53 = taux3
      tcoff54 = taux3
      tcoff55 = taux4
      tcoff56 = taux4
      tcoff57 = taux3
      tcoff58 = taux3
      tcoff59 = taux4
      tcoff60 = taux4
      tcoff61 = taux3
      tcoff62 = taux3
      tcoff63 = -taux4
      tcoff64 = -taux4

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet33 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet34 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff34)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet36 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet38 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet39 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet40 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff40)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,17.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet42 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet43 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet44 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet45 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet46 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet47 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet48 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff48)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet49 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff49)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet50 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff50)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet51 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff51)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet52 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff52)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet53 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff53)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet54 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff54)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet55 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff55)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet56 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff56)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet57 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff57)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet58 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff58)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet59 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff59)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet60 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff60)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet61 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff61)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet62 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff62)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet63 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff63)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet64 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff64)

CCC Sum

      t17 = tFjet34*tcoff34+tFjet33*tcoff33+tFjet49*tcoff49+tFjet48*tcof
     #f48+tFjet47*tcoff47+tFjet42*tcoff42+tFjet41*tcoff41+tFjet53*tcoff5
     #3+tFjet52*tcoff52+tFjet51*tcoff51+tFjet50*tcoff50+tFjet44*tcoff44+
     #tFjet43*tcoff43+tFjet40*tcoff40+tFjet39*tcoff39+tFjet38*tcoff38
      t34 = tFjet37*tcoff37+tFjet36*tcoff36+tFjet35*tcoff35+tFjet64*tcof
     #f64+tFjet63*tcoff63+tFjet62*tcoff62+tFjet61*tcoff61+tFjet46*tcoff4
     #6+tFjet45*tcoff45+tFjet57*tcoff57+tFjet56*tcoff56+tFjet55*tcoff55+
     #tFjet54*tcoff54+tFjet60*tcoff60+tFjet59*tcoff59+tFjet58*tcoff58
      qiqbirmnum2e1 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum2e2(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux3 = 0
      t3 = l3*s1vr
      t4 = s2l1*s13
      t8 = z**2
      t13 = r*t
      t14 = -1+t13
      t15 = t14**2
      t17 = r**2
      t18 = 1/t17
      t20 = t17*r
      t21 = t20*t
      t22 = l1**2
      t23 = dsqrt(l2)
      t24 = t22*t23
      t26 = t**2
      t27 = t17*t26
      t31 = l1*t23
      t32 = t17**2
      t38 = t23*t20
      t40 = t23*t17
      t42 = t23*r
      t49 = (1+t17)**2
      t52 = t17*l1
      t53 = t52*t
      t55 = l1*r
      t57 = l1*t
      t61 = dsqrt(-l1*(l1-1)*(1-l2)*t49*r*(t17*t+t53+r+r*t26-t55-t55*t26
     #+t+t57))
      t65 = t21*t24-t24*t27-t24*t17+t13*t24+t31*t32+2*t31*t27+4*t31*t17+
     #t31-t38*t-t40*t26-t40-t42*t+2*t61-4*t61*l4
      t66 = t65**2
      t67 = r+t
      t68 = t67**2
      t71 = 1/t49
      t79 = s13*l3
      t93 = l1*l2
      t97 = l2*t17
      t105 = l2*t22
      t110 = t61*t23
      t114 = -l1-2*t52-l1*t32+t93+4*t93*t17+t93*t32-t97-l2*t20*t-l2*r*t-
     #t97*t26+2*t93*t27-t105*t17+t105*t21+t105*t13-t105*t27-4*t110*l4+2*
     #t110
      t120 = dsqrt((r-t)*s1l1*s2l1*t114*t14/r*t71)
      t124 = t23*s1l1
      t131 = s1l1*l1
      t142 = s4vr**2
      taux4 = (-1024*(-2*z*s1l1-t3*t4+t3*t4*z-t8-t3*s2l1+t3*s2l1*z)*t15*
     #t18*t66/t68*t71-1024*s1vr*t14*t18*t65/t67*t71*(-t79*s1vr+t79*s1vr*
     #z-s13+s34*s13-t3+t3*z+s34+2*z-1)*(-t38*s1l1-2*t20*t120*c3l1+t124*t
     #53+t42*s1l1*t26-2*t120*c3l1*r-t42*t131-t42*t131*t26+t124*t57))/t14
     #2

CCC Coefficients

      tcoff33 = taux3
      tcoff34 = taux3
      tcoff35 = -taux4
      tcoff36 = -taux4
      tcoff37 = taux3
      tcoff38 = taux3
      tcoff39 = taux4
      tcoff40 = taux4
      tcoff41 = taux3
      tcoff42 = taux3
      tcoff43 = taux4
      tcoff44 = taux4
      tcoff45 = taux3
      tcoff46 = taux3
      tcoff47 = taux4
      tcoff48 = taux4
      tcoff49 = taux3
      tcoff50 = taux3
      tcoff51 = -taux4
      tcoff52 = -taux4
      tcoff53 = taux3
      tcoff54 = taux3
      tcoff55 = taux4
      tcoff56 = taux4
      tcoff57 = taux3
      tcoff58 = taux3
      tcoff59 = taux4
      tcoff60 = taux4
      tcoff61 = taux3
      tcoff62 = taux3
      tcoff63 = -taux4
      tcoff64 = -taux4

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet33 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet34 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff34)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet36 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet38 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet39 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet40 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff40)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,17.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet42 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet43 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet44 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet45 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet46 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet47 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet48 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff48)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet49 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff49)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet50 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff50)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet51 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff51)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet52 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff52)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet53 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff53)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet54 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff54)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet55 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff55)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet56 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff56)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet57 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff57)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet58 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff58)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet59 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff59)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet60 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff60)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet61 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff61)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet62 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff62)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet63 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff63)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet64 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff64)

CCC Sum

      t17 = tFjet56*tcoff56+tFjet35*tcoff35+tFjet55*tcoff55+tFjet34*tcof
     #f34+tFjet50*tcoff50+tFjet36*tcoff36+tFjet37*tcoff37+tFjet46*tcoff4
     #6+tFjet51*tcoff51+tFjet42*tcoff42+tFjet43*tcoff43+tFjet38*tcoff38+
     #tFjet57*tcoff57+tFjet47*tcoff47+tFjet58*tcoff58+tFjet59*tcoff59
      t34 = tFjet60*tcoff60+tFjet61*tcoff61+tFjet52*tcoff52+tFjet53*tcof
     #f53+tFjet54*tcoff54+tFjet62*tcoff62+tFjet63*tcoff63+tFjet48*tcoff4
     #8+tFjet39*tcoff39+tFjet40*tcoff40+tFjet64*tcoff64+tFjet44*tcoff44+
     #tFjet33*tcoff33+tFjet45*tcoff45+tFjet41*tcoff41+tFjet49*tcoff49
      qiqbirmnum2e2 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum2e3(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux3 = 0
      taux4 = 0

CCC Coefficients

      tcoff33 = taux3
      tcoff34 = taux3
      tcoff35 = -taux4
      tcoff36 = -taux4
      tcoff37 = taux3
      tcoff38 = taux3
      tcoff39 = taux4
      tcoff40 = taux4
      tcoff41 = taux3
      tcoff42 = taux3
      tcoff43 = taux4
      tcoff44 = taux4
      tcoff45 = taux3
      tcoff46 = taux3
      tcoff47 = taux4
      tcoff48 = taux4
      tcoff49 = taux3
      tcoff50 = taux3
      tcoff51 = -taux4
      tcoff52 = -taux4
      tcoff53 = taux3
      tcoff54 = taux3
      tcoff55 = taux4
      tcoff56 = taux4
      tcoff57 = taux3
      tcoff58 = taux3
      tcoff59 = taux4
      tcoff60 = taux4
      tcoff61 = taux3
      tcoff62 = taux3
      tcoff63 = -taux4
      tcoff64 = -taux4

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet33 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet34 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff34)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet36 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet38 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet39 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet40 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff40)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,17.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet42 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet43 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet44 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet45 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet46 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet47 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet48 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff48)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet49 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff49)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet50 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff50)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet51 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff51)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet52 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff52)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet53 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff53)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet54 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff54)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet55 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff55)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet56 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff56)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet57 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff57)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet58 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff58)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet59 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff59)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet60 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff60)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet61 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff61)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet62 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff62)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet63 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff63)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet64 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff64)

CCC Sum

      t17 = tFjet54*tcoff54+tFjet42*tcoff42+tFjet57*tcoff57+tFjet43*tcof
     #f43+tFjet55*tcoff55+tFjet44*tcoff44+tFjet56*tcoff56+tFjet45*tcoff4
     #5+tFjet58*tcoff58+tFjet59*tcoff59+tFjet46*tcoff46+tFjet60*tcoff60+
     #tFjet61*tcoff61+tFjet47*tcoff47+tFjet62*tcoff62+tFjet48*tcoff48
      t34 = tFjet35*tcoff35+tFjet34*tcoff34+tFjet33*tcoff33+tFjet37*tcof
     #f37+tFjet36*tcoff36+tFjet49*tcoff49+tFjet38*tcoff38+tFjet63*tcoff6
     #3+tFjet50*tcoff50+tFjet39*tcoff39+tFjet40*tcoff40+tFjet51*tcoff51+
     #tFjet64*tcoff64+tFjet52*tcoff52+tFjet41*tcoff41+tFjet53*tcoff53
      qiqbirmnum2e3 = t17+t34
      RETURN
      END

      double precision function qiqbirmnum2e4(var,ndim,z,ux,l1,
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
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,
     & c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux3 = 0
      taux4 = 0

CCC Coefficients

      tcoff33 = taux3
      tcoff34 = taux3
      tcoff35 = -taux4
      tcoff36 = -taux4
      tcoff37 = taux3
      tcoff38 = taux3
      tcoff39 = taux4
      tcoff40 = taux4
      tcoff41 = taux3
      tcoff42 = taux3
      tcoff43 = taux4
      tcoff44 = taux4
      tcoff45 = taux3
      tcoff46 = taux3
      tcoff47 = taux4
      tcoff48 = taux4
      tcoff49 = taux3
      tcoff50 = taux3
      tcoff51 = -taux4
      tcoff52 = -taux4
      tcoff53 = taux3
      tcoff54 = taux3
      tcoff55 = taux4
      tcoff56 = taux4
      tcoff57 = taux3
      tcoff58 = taux3
      tcoff59 = taux4
      tcoff60 = taux4
      tcoff61 = taux3
      tcoff62 = taux3
      tcoff63 = -taux4
      tcoff64 = -taux4

CCC Jet functions

      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet33 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff33)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet34 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff34)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff35)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet36 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff36)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet37 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff37)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet38 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff38)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet39 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff39)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet40 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff40)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet41 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,17.D0,1.D0,1.D0,Tcoeff*tcoff41)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet42 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,17.D0,1.D0,1.D0,Tcoeff*tcoff42)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet43 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,18.D0,1.D0,1.D0,Tcoeff*tcoff43)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet44 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,18.D0,1.D0,1.D0,Tcoeff*tcoff44)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet45 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,15.D0,1.D0,1.D0,Tcoeff*tcoff45)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet46 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,15.D0,1.D0,1.D0,Tcoeff*tcoff46)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet47 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,16.D0,1.D0,1.D0,Tcoeff*tcoff47)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet48 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,16.D0,1.D0,1.D0,Tcoeff*tcoff48)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet49 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,21.D0,1.D0,1.D0,Tcoeff*tcoff49)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet50 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,21.D0,1.D0,1.D0,Tcoeff*tcoff50)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet51 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,22.D0,1.D0,1.D0,Tcoeff*tcoff51)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet52 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,22.D0,1.D0,1.D0,Tcoeff*tcoff52)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet53 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,23.D0,1.D0,1.D0,Tcoeff*tcoff53)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet54 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,23.D0,1.D0,1.D0,Tcoeff*tcoff54)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet55 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,24.D0,1.D0,1.D0,Tcoeff*tcoff55)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet56 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,24.D0,1.D0,1.D0,Tcoeff*tcoff56)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet57 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,19.D0,1.D0,1.D0,Tcoeff*tcoff57)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet58 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,19.D0,1.D0,1.D0,Tcoeff*tcoff58)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet59 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,20.D0,1.D0,1.D0,Tcoeff*tcoff59)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet60 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,20.D0,1.D0,1.D0,Tcoeff*tcoff60)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet61 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,25.D0,1.D0,1.D0,Tcoeff*tcoff61)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet62 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,25.D0,1.D0,1.D0,Tcoeff*tcoff62)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s1vr
      t4 = t2*s23r
      t5 = 1-l3
      t7 = t1**2
      t9 = l3*t5*t7*s34r
      t12 = t5*t1*s13r
      tFjet63 = Fjet(var,ndim,z,ux,t3-t4-t9,-1-t12-t3+z,t4,t12,t9,-1+z-t
     #3+t9,t3,1+t12+t4,-t12-t4-t9+z,s2l1,s1l1,s4l1,s3l1,s2l2,s1l2,s4l2,s
     #3l2,26.D0,1.D0,1.D0,Tcoeff*tcoff63)
      t2 = 1-l3
      t3 = 1-z
      t5 = t2*t3*s13r
      t6 = l3*t3
      t7 = t6*s23r
      t8 = t6*s1vr
      t11 = t3**2
      t13 = l3*t2*t11*s34r
      tFjet64 = Fjet(var,ndim,z,1/ux,t5,t7,-1-t5-t8+z,t8-t7-t13,t13,t8,-
     #1+z-t8+t13,-t5-t7-t13+z,1+t5+t7,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2
     #,s4l2,26.D0,1.D0,1.D0,Tcoeff*tcoff64)

CCC Sum

      t17 = tFjet33*tcoff33+tFjet60*tcoff60+tFjet40*tcoff40+tFjet34*tcof
     #f34+tFjet51*tcoff51+tFjet52*tcoff52+tFjet53*tcoff53+tFjet49*tcoff4
     #9+tFjet39*tcoff39+tFjet50*tcoff50+tFjet59*tcoff59+tFjet41*tcoff41+
     #tFjet61*tcoff61+tFjet42*tcoff42+tFjet62*tcoff62+tFjet63*tcoff63
      t34 = tFjet35*tcoff35+tFjet43*tcoff43+tFjet64*tcoff64+tFjet44*tcof
     #f44+tFjet45*tcoff45+tFjet54*tcoff54+tFjet55*tcoff55+tFjet56*tcoff5
     #6+tFjet36*tcoff36+tFjet46*tcoff46+tFjet47*tcoff47+tFjet37*tcoff37+
     #tFjet57*tcoff57+tFjet48*tcoff48+tFjet38*tcoff38+tFjet58*tcoff58
      qiqbirmnum2e4 = t17+t34
      RETURN
      END

