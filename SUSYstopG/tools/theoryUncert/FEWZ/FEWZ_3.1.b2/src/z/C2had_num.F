CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C2num1e0(var,ndim,z,ux,l1,
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
      t19 = t2*t6
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
      t65 = 2*z
      t66 = 1+t65
      t75 = z**2
      t76 = t2*t75
      t86 = t1*t6
      t87 = t86*l3
      t88 = s1vr**2
      t94 = t2*t7
      t95 = l3*t10
      t96 = t95*s34r
      t97 = t94*t96
      t99 = z+1
      t101 = t5*t10
      t104 = 768*t97-256*t99*t1*t101
      t108 = s23r**2
      t109 = t12*t108
      t113 = t1*t5
      t119 = l3*t3
      t120 = t6**2
      t122 = s34r**2
      t123 = t120*t10*t122
      t131 = t95*s2l1
      t132 = t86*t131
      t134 = t95*s3l1
      t135 = t86*t134
      t136 = 256*t135
      t137 = l3*z
      t139 = t86*t137*t10
      t157 = t10*s3l1
      t160 = z+2
      t162 = z*t5
      t163 = t162*t10
      t170 = t137*t6
      t171 = t1*t10
      t173 = t171*t11*t108
      t176 = t38*s34r
      t177 = t53*t176
      t179 = t1*z
      t186 = z*t10
      t187 = t186*t122
      t195 = l3**2
      t196 = t195*t7
      t197 = t196*z
      t198 = t88*s1vr
      t206 = s2l1**2
      t207 = t95*t206
      t216 = 128*t195*t120*z*t176
      t217 = t75*l3
      t219 = t217*t86*t10
      t220 = 256*t219
      t230 = t186*s3l1
      t231 = t87*t230
      t247 = z*t120*t5
      t249 = t195*t10
      t253 = t76*t7
      t254 = t253*t96
      t256 = s3l1**2
      t260 = t75*t1
      t264 = 2*t75
      t265 = 1+t264+t65
      t273 = t217*t6
      t284 = t95*t122
      t286 = t75*t6
      t293 = t6*t195
      t294 = t293*t10
      t298 = t171*s34r
      t299 = t196*t298
      t300 = l3*t5
      t301 = t300*t186
      t303 = -256*t299+256*t301
      t306 = t86*t96
      t308 = 256*t306-256*t186
      t313 = t12*t108*s23r
      t316 = t300*t10
      t341 = t293*t186
      t382 = -t5
      t396 = t197*t298
      t402 = t2*t120*t195
      t403 = t402*t187
      t405 = t186*s34r
      t422 = t195*l3
      t423 = z*t7*t422
      t441 = 128*t1*t120*t422*t405
      t442 = t10*t256
      t461 = t75*t5
      t462 = t461*t95
      t480 = -t99*z
      t497 = t382*z
      t539 = t75*t10
      s3 = -64*t8*t12*s1vr*t13*s13r
      s5 = ((-256*t19*t12*s1vr+128*t19*t23+256*t3*t7*t28)*t32+((-256*t19
     #*t34-384*t37*t38)*t11*s1vr+128*t44*t6*t23+256*t8*t28)*s1l1+(64*t53
     #*t38*t54+(-128*t44*t58*s2l1-128*t44*t58*s3l1-64*t66*t2*t37*t10)*t1
     #1)*s1vr+64*t76*t6*t23+128*t75*t7*t3*t28)*t13
      s8 = (-256*t87*t12*t88+(256*t87*t23+t104*t11)*s1vr-256*t87*t109+(-
     #512*t97+256*t113*t10)*t11*s23r+(-512*t119*t123+512*t19*t27)*t11)*t
     #32+((-512*t132+t136-384*t139)*t11*t88+((256*t132-t136+512*t139)*t1
     #1*s23r+(t104*s2l1+(-256*t94*t134+768*t52*t44*t10)*s34r-256*t86*t15
     #7-256*t160*t1*t163)*t11)*s1vr-256*t170*t173+(-512*t177+256*t179*t1
     #01)*t11*s23r+(-512*t119*t120*t187+512*t37*t176)*t11)*s1l1-64*t197*
     #t171*t11*t198+(128*t197*t171*t54+(-256*t86*t207+(256*t135-256*t139
     #)*s2l1+t216-t220)*t11)*t88
      s7 = s8+(-64*t197*t173+(128*t170*t171*s2l1-t216+128*t231+t220)*t11
     #*s23r+(-128*t179*t101*t206+(256*t177-128*t160*t1*t163)*s2l1-128*t2
     #47*t3*t249*t122+384*t254-128*t179*t101*t256-128*t260*t101*s3l1-64*
     #t265*t1*t163)*t11)*s1vr-128*t273*t173+(-256*t254+128*t260*t101)*t1
     #1*s23r+(-256*t75*t120*t3*t284+256*t286*t176)*t11
      s8 = s13r
      s6 = s7*s8
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(128*t294*t54*t88+(t303*t11*s23r+t308*t11)*s1vr+128*t293*t
     #313+(256*t299-256*t316)*t11*t108+(256*t195*t2*t123-512*t306+128*t1
     #0)*t11*s23r+(-512*l3*t2*t7*t10*t122+256*t113*t27)*t11)*t32+(((256*
     #t293*t34+256*t341)*t11*s23r+(256*t300*t157-128*t301)*t11)*t88+((25
     #6*t293*t157-128*t341)*t11*t108+(t303*s2l1+(256*t196*t171*s3l1-256*
     #t196*t179*t10)*s34r+256*(z-2)*t5*t134+256*t99*z*t316)*t11*s23r+(t3
     #08*s2l1+(-768*t135+256*t139)*s34r+256*t382*t10*s3l1-128*t66*z*t10)
     #*t11)*s1vr+128*t37*t195*t313+(256*t396-256*t301)*t11*t108+(256*t40
     #3-512*t87*t405+128*t186)*t11*s23r+(-512*t53*t38*t122+256*t113*t405
     #)*t11)*s1l1+(64*t423*t23-64*t293*t186*t11)*t198
      taux1 = s1+(-128*t423*t109+(128*t293*t10*t206+128*t293*t186*s2l1-t
     #441+128*t293*t442+128*t293*t230+64*(2+3*z)*z*t294)*t11*s23r+(256*t
     #300*t157*s2l1+128*t396-256*t300*t442-128*t462)*t11)*t88+(64*t423*t
     #313+(t441-64*t66*z*t294)*t11*t108+(128*t162*t207-128*t480*t5*t131+
     #128*t247*t2*t422*t10*t122+128*t480*t195*t7*t1*t27+128*t162*t95*t25
     #6-128*t497*t5*t134+64*t265*z*t316)*t11*s23r+(-128*t186*t206+128*t4
     #80*t34-128*t403+(-256*t231+128*t219)*s34r-128*t186*t256+128*t497*t
     #157-64*(1+t264)*z*t10)*t11)*s1vr+64*t286*t195*t313+(128*t260*t7*t2
     #49*s34r-128*t462)*t11*t108+(128*t402*t539*t122-256*t273*t298+64*t5
     #39)*t11*s23r+(-256*t253*t284+128*t461*t298)*t11
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
      t48 = 2*z
      t59 = z**2
      t60 = t2*t59
      t70 = t1*t6
      t71 = t70*l3
      t72 = z*t10
      t73 = s1vr**2
      t81 = t35*s34r
      t82 = t34*t81
      t84 = t1*t59
      t85 = t5*t10
      t86 = t84*t85
      t93 = z*l3*t6
      t94 = t1*t10
      t95 = s23r**2
      t97 = t94*t11*t95
      t101 = -t1
      t109 = t6**2
      t111 = s34r**2
      t112 = t72*t111
      t120 = l3**2
      t122 = t120*t7*z
      t123 = t73*s1vr
      t132 = t93*t94*s2l1
      t137 = 512*t120*t109*z*t81
      t138 = t59*l3
      t140 = t138*t70*t10
      t141 = 1024*t140
      t149 = t72*s3l1
      t150 = t71*t149
      t160 = z*t109*t5
      t162 = t120*t10
      t166 = t60*t7
      t167 = l3*t10
      t169 = t166*t167*s34r
      t176 = 2*t59
      t177 = 1-t48-t176
      t187 = t138*t6
      t199 = t167*t111
      t201 = t59*t6
      t208 = t120*t6
      t209 = t208*z
      t212 = l3*t5
      t213 = t72*t11
      t218 = t12*t95
      t221 = t94*s34r
      t222 = t122*t221
      t223 = t59*t5
      t224 = t223*t167
      t225 = t222-t224
      t229 = t72*s34r
      t230 = t71*t229
      t241 = t12*t95*s23r
      t250 = t2*t109*t120
      t251 = t250*t112
      t269 = t120*l3
      t270 = z*t7*t269
      t284 = 512*t1*t109*t269*t229
      t290 = t208*t10
      t310 = -z+t59
      t311 = t310*t5
      t370 = t59*t10
      s2 = -256*t8*t12*s1vr*t13*s13r+((512*t20*t12*s1vr-512*t20*t24-1024
     #*t8*t28)*s1l1+(256*t34*t35*t36+(512*t40*t41*s2l1-512*t40*t41*s3l1-
     #256*(-t48+1)*z*t6*t2*t10)*t11)*s1vr-256*t60*t6*t24-512*t59*t7*t3*t
     #28)*t13
      s3 = s2+((1536*t71*t72*t11*t73+(-1024*t71*t72*t36+(-3072*t82+1024*
     #t86)*t11)*s1vr+1024*t93*t97+(2048*t82+1024*t101*z*t85)*t11*s23r+(2
     #048*l3*t3*t109*t112-2048*t19*t81)*t11)*s1l1+256*t122*t94*t11*t123+
     #(-512*t122*t94*t36+(1024*t132-t137+t141)*t11)*t73+(256*t122*t97+(-
     #512*t132+t137+512*t150-t141)*t11*s23r+((-1024*t82+512*t86)*s2l1+51
     #2*t160*t3*t162*t111-1536*t169+512*(-t48+t59)*t5*t94*s3l1-256*t177*
     #t1*z*t5*t10)*t11)*s1vr+512*t187*t97+(1024*t169+512*t101*t59*t85)*t
     #11*s23r+(1024*t59*t109*t3*t199-1024*t201*t81)*t11)*s13r
      s1 = s3+((-1024*t209*t24+512*t212*t213)*t73+(512*t209*t218+1024*t2
     #25*t11*s23r+(-1024*t230-256*(t48-4*t59)*t10)*t11)*s1vr-512*t209*t2
     #41+(-1024*t222+1024*t212*t72)*t11*t95+(-1024*t251+2048*t230-512*t7
     #2)*t11*s23r+(2048*t34*t35*t111-1024*t1*t5*t229)*t11)*s1l1+(-256*t2
     #70*t24+256*t208*t213)*t123
      taux2 = s1+(512*t270*t218+(-512*t208*t72*s2l1+t284-512*t208*t149+2
     #56*(-2-3*z)*z*t290)*t11*s23r-512*t225*t11)*t73+(-256*t270*t241+(-t
     #284+256*(1+t48)*z*t290)*t11*t95+(-512*t311*t167*s2l1-512*t160*t2*t
     #269*t10*t111+512*(z+1)*z*t120*t7*t1*t27-512*t311*t167*s3l1+256*t17
     #7*z*t212*t10)*t11*s23r+(512*z*t5*t10*s2l1+512*t251+(1024*t150-512*
     #t140)*s34r+512*t310*t10*s3l1-256*(1-t176)*z*t10)*t11)*s1vr-256*t20
     #1*t120*t241+(-512*t84*t7*t162*s34r+512*t224)*t11*t95+(-512*t250*t3
     #70*t111+1024*t187*t221-256*t370)*t11*s23r+(1024*t166*t199-512*t223
     #*t221)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet2 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet3 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet4 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      C2num1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function C2num1e1(var,ndim,z,ux,l1,
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
      t33 = l3*t7*z
      t34 = t2*t10
      t38 = t6*t10
      t42 = t6**2
      t45 = z*t10
      t46 = t45*s34r
      t52 = t6*z
      t59 = z**2
      t60 = t2*t59
      t62 = t12*s23r
      t74 = l3*t6*t1
      t75 = s23r*s1vr
      t78 = s23r**2
      t79 = t12*t78
      t81 = t1*t5
      t86 = t1*t6
      t87 = l3*t10
      t92 = l3*z
      t99 = t92*t6
      t100 = t1*t10
      t102 = t100*t11*t78
      t105 = t1*z
      t111 = l3**2
      t113 = t111*t7*z
      t114 = s1vr**2
      t115 = t114*s1vr
      t124 = t99*t100*s2l1
      t127 = t34*s34r
      t128 = t111*t42*z*t127
      t129 = t45*s3l1
      t130 = t74*t129
      t131 = 2*z
      t135 = t86*t10
      t155 = t5*t10
      t156 = s2l1**2
      t160 = t81*t129
      t166 = z*t42*t5
      t168 = t111*t10
      t169 = s34r**2
      t173 = t60*t7
      t175 = t173*t87*s34r
      t177 = s3l1**2
      t185 = z*t5
      t194 = t59*l3*t6
      t208 = t87*t169
      t210 = t59*t6
      t217 = t6*t111
      t218 = t217*t10
      t223 = l3*t5
      t229 = t12*t78*s23r
      t237 = t10*s2l1
      t238 = t217*t237
      t240 = t10*s3l1
      t241 = t217*t240
      t243 = t217*t45
      t255 = t223*t45
      t273 = t111*l3
      t274 = z*t7*t273
      t292 = 128*t1*t42*t273*t46
      t303 = t100*s34r
      t305 = t223*t129
      t338 = 128*t1*t59*t7*t168*s34r
      t359 = t2*t42*t111
      t380 = t59*t5
      t386 = t59*t10
      s3 = 128*t3*z*t7*t12*s1vr*t13*s13r
      s5 = (-128*t2*t6*t10*t21*t22-128*t26*t6*t12*s23r*s1l1+(-256*t33*t3
     #4*t21+(128*t26*t38*s2l1-128*t3*t42*l3*t46+128*t26*t38*s3l1+256*t52
     #*t34)*t11)*s1vr-128*t60*t6*t62-128*t59*t7*t3*t10*s34r*t11)*t13
      s6 = ((256*t74*t12*t75+256*t74*t79-256*t81*t62)*t22+((256*t86*t87*
     #s2l1+256*t86*t87*s3l1+256*t86*t92*t10)*t11*t75+256*t99*t102-256*t1
     #05*t5*t62)*s1l1+128*t113*t100*t11*t115+(-128*t113*t100*t21+(-128*t
     #124-128*t128-128*t130+128*(-1+t131)*z*l3*t135)*t11)*t114+(256*t113
     #*t102+(-128*t124+256*t128-128*t130-128*(z+2)*z*l3*t135)*t11*s23r+(
     #128*t105*t155*t156+(256*t160+256*t105*t155)*s2l1+128*t166*t3*t168*
     #t169-384*t175+128*t105*t155*t177+256*t160-64*(t131-2*t59-4)*t1*t18
     #5*t10)*t11)*s1vr+256*t194*t102+(256*t175-256*t1*t59*t155)*t11*s23r
     #+(256*t59*t42*t3*t208-256*t210*t127)*t11)*s13r
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-128*t218*t21*t114+(-256*t217*t79+256*t223*t62)*s1vr-128*
     #t217*t229+256*t223*t79-128*t62)*t22+((-256*t238-256*t241-128*t243)
     #*t11*s23r*t114+((-256*t238-256*t241-256*t243)*t11*t78+(256*t223*t2
     #37+256*t223*t240+256*t255)*t11*s23r)*s1vr-128*t52*t111*t229+256*t2
     #23*z*t79-128*t45*t21)*s1l1+(-128*t274*t62+128*t217*t45*t11)*t115
      taux1 = s1+(128*t274*t79+(-128*t217*t10*t156-256*t217*t240*s2l1+t2
     #92-128*t217*t10*t177-256*t210*t168)*t11*s23r+(-128*t223*t45*s2l1-1
     #28*t113*t303-128*t305+128*t52*t87)*t11)*t114+(-128*t274*t229+(-t29
     #2+64*(2+2*z)*z*t218)*t11*t78+(-128*t185*t87*t156+(-256*t305-128*t2
     #55)*s2l1-128*t166*t2*t273*t10*t169+t338-128*t185*t87*t177-128*t305
     #-64*(2+2*t59)*z*t223*t10)*t11*s23r+(128*t45*t156+(256*t129+128*t45
     #)*s2l1+128*t359*t45*t169-128*t33*t303+128*t45*t177+128*t129-64*(-2
     #*t59+2*z-2)*z*t10)*t11)*s1vr-128*t210*t111*t229+(-t338+256*t380*t8
     #7)*t11*t78+(-128*t359*t386*t169+256*t194*t303-128*t386)*t11*s23r+(
     #256*t173*t208-128*t380*t303)*t11
      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = l3-1
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t17 = t6*z*t3
      t18 = t9*t11
      t23 = s13r**2
      t26 = t5*t3
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t35 = t2**2
      t40 = 1024*l3*t35*t6*t9*s34r
      t41 = 6*z
      t42 = t41-1
      t44 = t5*t9
      t51 = t2*z
      t56 = t9*s34r*t11
      t61 = l3*t9
      t65 = t61*s3l1
      t68 = l3*t3
      t69 = z*t5
      t76 = t68*z
      t81 = 3*z
      t89 = t6*t35*l3
      t90 = z*t9
      t91 = t90*s34r
      t98 = t41-7
      t108 = z**2
      t109 = t5*t108
      t119 = t4*t3
      t120 = l3**2
      t121 = t119*t120
      t122 = t28*s1vr
      t126 = -t4
      t127 = t126*t3
      t131 = t5*t35
      t132 = t120*t9
      t133 = t132*s34r
      t134 = t131*t133
      t136 = t4*t2
      t137 = l3*z
      t145 = s23r**2
      t146 = t18*t145
      t150 = 8*z
      t153 = t136*t9
      t159 = t35*t1
      t160 = t6*t159
      t161 = s34r**2
      t163 = t160*t132*t161
      t165 = t44*s34r
      t166 = t76*t165
      t168 = 12*t108
      t171 = t1*t9
      t179 = t4*t9
      t181 = t179*t11*t145
      t185 = 4*t4
      t192 = t90*t161
      t200 = t132*s2l1
      t203 = 3*t126
      t205 = t132*s3l1
      t208 = t120*t3
      t210 = z*t4*t9
      t211 = t208*t210
      t216 = t179*s3l1
      t230 = t131*t205
      t233 = t131*t120*z*t9
      t237 = t136*l3
      t238 = t90*s3l1
      t239 = t237*t238
      t241 = t108*l3
      t248 = t208*z
      t263 = 7*z
      t275 = t61*s34r
      t278 = 3*t108
      t287 = z*t159
      t306 = 6*t108
      t317 = t241*t2
      t320 = t109*t3
      t331 = t61*t161
      t333 = t108*t2
      t340 = t120*l3
      t341 = t3*t340
      t343 = t2*t120
      t350 = t126*t35
      t351 = t340*t9
      t353 = t350*t351*s34r
      t354 = 1024*t353
      t355 = t343*t90
      t360 = t127*t133
      t362 = 1-t81
      t371 = t18*t145*s23r
      t374 = 10*z
      t383 = t351*t161
      t386 = t179*s34r
      t388 = 3072*t248*t386
      t400 = t120*t5*t35*t9*t161
      t416 = l3*t1
      t422 = t131*t120
      t441 = t9*s2l1
      t443 = 256*t341*t441
      t444 = t9*s3l1
      t446 = 256*t341*t444
      t447 = t341*t90
      t477 = 5*z
      t484 = t343*t9
      t491 = z+1
      t508 = t416*t9
      t519 = 2-z
      t536 = t519*t4
      t612 = t108*t1
      t618 = t108*t9
      s2 = (-256*t3*t6*t9*t11*s1vr*s1l1+768*t17*t18*s1vr)*t23*s13r+((-51
     #2*t27*t18*t28+(768*t27*t32+(t40-256*t42*t2*t44)*t11)*s1vr+2048*t51
     #*t5*t32+3072*t17*t56)*s1l1+(768*t26*t61*s2l1+256*t26*t65-256*t68*t
     #69*t9)*t11*t28+(-512*t76*t44*t11*s23r+((-t40-256*(t81-2)*t2*t44)*s
     #2l1+512*t89*t91+1792*t69*t2*t9*s3l1-256*t98*z*t2*t5*t9)*t11)*s1vr+
     #1024*t109*t2*t32+1536*t108*t3*t6*t56)*t23
      s4 = s2
      s7 = (-256*t121*t18*t122+(-1024*t127*t120*t32+(1024*t134-4608*t136
     #*t137*t9)*t11)*t28+(-768*t121*t146+(-2048*t134-512*(1-t150)*l3*t15
     #3)*t11*s23r+(-1024*t163+9216*t166-256*(t168-1)*t4*t171)*t11)*s1vr-
     #4096*t137*t2*t181+(-6144*t166+1024*t185*z*t171)*t11*s23r+(-6144*t8
     #9*t192+6144*t51*t165)*t11)*s1l1+(-256*t119*t200-256*t203*t3*t205-1
     #024*t211)*t11*t122+((-512*t208*t216+2304*t211)*t11*s23r+(256*(-11*
     #z+2)*l3*t2*t179*s2l1+(-1536*t230+1536*t233)*s34r+768*t239-3072*t24
     #1*t153)*t11)*t28
      s6 = s7+(-1280*t248*t181+((512*t134+1024*l3*t1*t153)*s2l1-2048*t12
     #0*t35*z*t165-2048*t239-512*(-t263+2)*z*l3*t153)*t11*s23r+((512*t16
     #3-512*(4-t263)*t5*t3*t275+512*(2*z-t278+1)*t1*t179)*s2l1+(512*t160
     #*t205-1536*t287*t6*t120*t9)*t161+(-512*t76*t44*s3l1+4608*t109*t68*
     #t9)*s34r+512*(-1+t150-t278)*t1*t216-256*(-7+t41+t306)*t4*z*t1*t9)*
     #t11)*s1vr-2048*t317*t181+(-3072*t320*t275+512*t185*t108*t171)*t11*
     #s23r+(-3072*t108*t35*t6*t331+3072*t333*t165)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((256*t341*t32-256*t343*t18)*t122+(-512*t341*t146+(t354+35
     #84*t355)*t11*s23r+(-1024*t360+512*t362*t1*t61)*t11)*t28+(256*t341*
     #t371+(-t354+256*(-t374+1)*t2*t132)*t11*t145+(1024*t5*t159*t383-t38
     #8-256*(-t168-4*z+1)*t1*t61)*t11*s23r+(-1024*t400-1024*t362*l3*t2*t
     #386-256*(-t41+1+t168)*t9)*t11)*s1vr+2048*t51*t120*t371+(t388-4096*
     #t416*t90)*t11*t145+(3072*t422*t192-6144*t237*t91+2048*t90)*t11*s23
     #r+(-6144*t76*t44*t161+3072*t4*t1*t91)*t11)*s1l1+((t443-t446+1024*t
     #447)*t11*s23r+(-256*t343*t441+256*t343*t444-1024*t355)*t11)*t122
      s3 = s1
      s5 = ((-t443+t446-2048*t447)*t11*t145+((512*t353+1792*t355)*s2l1+(
     #-512*t350*t351*s3l1-1536*t340*t35*t210)*s34r+256*(-2+t477)*t2*t205
     #+256*(7+t374)*z*t484)*t11*s23r+((-512*t360+256*t491*t1*t61)*s2l1+(
     #1024*t127*t205+1536*t211)*s34r+256*(1+t81)*t1*t65-256*t42*z*t508)*
     #t11)*t28
      s6 = (1024*z*t3*t340*t371+(256*t519*t2*t200+1536*t4*t35*t340*t91+2
     #56*t343*t238+256*(-3-t150)*z*t484)*t11*t145+((512*t536*t3*t133-512
     #*(-t278+1+t81)*t1*t61)*s2l1+1536*t287*t5*t383+(512*t491*t4*t3*t205
     #+512*(-t81-2)*z*t120*t119*t9)*s34r-512*(-1-t278+t477)*t1*t65-256*(
     #4-t306-t150)*z*t508)*t11*s23r+((512*t400-512*t536*t2*t275-256*t98*
     #z*t9)*s2l1+(512*t230-1536*t233)*t161+(-3584*t239-512*t362*z*l3*t15
     #3)*s34r+256*(-t306+9*z-2)*t9*s3l1-256*(-3+t306)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+1024*t333*t120*t371+(-512*t203*t108*t3*t133-2048*t612*t
     #61)*t11*t145+(1536*t422*t618*t161-3072*t317*t386+1024*t618)*t11*s2
     #3r+(-3072*t320*t331+1536*t612*t386)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet2 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet3 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet4 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      C2num1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function C2num1e2(var,ndim,z,ux,l1,
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
      t54 = t48*t27
      t59 = z*t10*t11
      t64 = s23r**2
      t65 = t11*t64
      t66 = t48*t65
      t94 = z*t7*t45*l3
      t96 = t45*t10
      t102 = t12*t64
      t117 = t12*t64*s23r
      t121 = (-3+2*z)*z
      taux1 = -64*t2*t1*z*t7*t12*s1vr*t13*s13r+(128*t20*t21*t11*t22+(192
     #*t20*t28-192*t30*t2*t10*t11)*s1vr+64*t2*t37*t6*t40)*t13+(-64*t47*t
     #48*t11*t49+(-256*t47*t54+256*t1*t6*l3*t59)*t22+(-192*t47*t66-128*(
     #-3+z)*z*l3*t6*t54-192*t1*z*t5*t10*t11)*s1vr-128*t37*l3*t6*t66+128*
     #t1*t37*t5*t40)*s13r+(64*t94*t40-64*t30*t96*t11)*t49+(128*t94*t102+
     #64*(-4+z)*z*t6*t96*t27+128*l3*t5*t59)*t22+(64*t94*t117+64*t121*t6*
     #t96*t65-64*t121*t5*t28-64*t59)*s1vr+64*t37*t6*t45*t117-128*t37*t5*
     #l3*t102+64*t37*t10*t27
      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = l3-1
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t17 = t6*z*t3
      t18 = t9*t11
      t23 = s13r**2
      t26 = t5*t3
      t28 = t18*s23r
      t31 = t2**2
      t33 = t6*t9
      t34 = t33*s34r
      t36 = 2048*l3*t31*t34
      t37 = -1-z
      t40 = t5*t9
      t47 = t2*z
      t52 = t9*s34r*t11
      t57 = l3*t9
      t61 = t57*s3l1
      t66 = s1vr**2
      t69 = t5*z*t3
      t87 = 4*z
      t98 = z**2
      t99 = t5*t98
      t109 = -t4
      t110 = 4*t109
      t111 = t110*t3
      t112 = l3**2
      t113 = t111*t112
      t114 = t66*s1vr
      t118 = t5*t31
      t119 = t112*t9
      t120 = t119*s34r
      t121 = t118*t120
      t123 = 6*z
      t126 = t2*t4
      t127 = t126*t9
      t133 = s23r**2
      t138 = t123+4
      t145 = t31*t1
      t147 = s34r**2
      t149 = t6*t145*t119*t147
      t151 = l3*t3
      t152 = t151*z
      t153 = t40*s34r
      t154 = t152*t153
      t156 = 8*t98
      t159 = t1*t9
      t166 = t126*l3
      t167 = z*t9
      t168 = t11*t133
      t173 = 5*t109
      t182 = t167*t147
      t190 = t119*s2l1
      t195 = t119*s3l1
      t198 = t112*t3
      t200 = z*t4*t9
      t201 = t198*t200
      t206 = t4*t3
      t218 = t4*t9
      t222 = t118*t195
      t224 = t112*z
      t226 = t118*t224*t9
      t230 = -t138
      t235 = 2*z
      t245 = t224*t3
      t246 = t218*t168
      t259 = t167*s3l1
      t260 = t166*t259
      t262 = 5*z
      t271 = 3*z
      t275 = t57*s34r
      t302 = t218*s3l1
      t305 = 4*t98
      t317 = t98*l3*t2
      t320 = t99*t3
      t331 = t57*t147
      t333 = t98*t2
      t340 = t112*l3
      t341 = t3*t340
      t343 = t2*t112
      t348 = t4*t31
      t349 = t340*t9
      t351 = t348*t349*s34r
      t352 = 2048*t351
      t353 = 7*z
      t361 = t218*s34r
      t362 = t198*t361
      t373 = t18*t133*s23r
      t376 = 8*z
      t384 = t5*t145
      t389 = 2048*t245*t361
      t400 = t112*t5*t31*t9*t147
      t414 = l3*t1
      t420 = t118*t112
      t423 = t167*s34r
      t440 = t9*s2l1
      t443 = t341*t167
      t449 = 1024*t343*t440
      t472 = t340*t31
      t484 = t343*t9
      t507 = t414*t9
      t514 = z*t3
      t531 = (z-2)*t4
      t535 = 2*t98
      t605 = t98*t1
      t611 = t98*t9
      s2 = (1024*t3*t6*t9*t11*s1vr*s1l1-256*t17*t18*s1vr)*t23*s13r+(((-3
     #072*t26*l3*t28+(-t36-1024*t37*t2*t40)*t11)*s1vr-2560*t47*t5*t28-20
     #48*t17*t52)*s1l1+(-2048*t26*t57*s2l1-1024*t26*t61)*t11*t66+(-768*t
     #69*t57*t11*s23r+((t36-512*(z+2)*t2*t40)*s2l1-1024*z*l3*t31*t34-153
     #6*t47*t40*s3l1-256*(7-t87)*z*t2*t5*t9)*t11)*s1vr-1280*t99*t2*t28-1
     #024*t98*t3*t6*t52)*t23
      s4 = s2
      s7 = (-256*t113*t18*t114+(-2048*t121+512*(-4+t123)*l3*t127)*t11*t6
     #6+(-768*t113*t18*t133+(4096*t121-512*t138*l3*t127)*t11*s23r+(2048*
     #t149-6144*t154-256*(-4-t156)*t4*t159)*t11)*s1vr+5120*t166*t167*t16
     #8+(4096*t154+1024*t173*z*t159)*t11*s23r+(4096*l3*t6*t31*t182-4096*
     #t47*t153)*t11)*s1l1+(-256*t111*t190+256*t110*t3*t195+1280*t201)*t1
     #1*t114+((1024*t206*t190-512*t111*t195-1536*t201)*t11*s23r+(256*(-8
     #+t123)*l3*t2*t218*s2l1+(3072*t222-1024*t226)*s34r+256*t230*t4*t2*t
     #61+1024*(-1+t235)*z*l3*t127)*t11)*t66
      s6 = s7+(2304*t245*t246+((-1024*t121-512*(-z-4)*l3*t127)*s2l1+2048
     #*t224*t31*t153+2560*t260-512*(-1+t262)*z*l3*t127)*t11*s23r+((-1024
     #*t149-1024*(-4+t271)*t5*t3*t275+1024*(-t235+t98-1)*t1*t218)*s2l1+(
     #-1024*t112*t145*t33*s3l1+1024*z*t145*t6*t112*t9)*t147+(1024*t69*t6
     #1-3072*t99*t151*t9)*s34r+1024*(t98+1-t87)*t1*t302-256*(7-t87-t305)
     #*t4*z*t1*t9)*t11)*s1vr+2560*t317*t246+(2048*t320*t275+512*t173*t98
     #*t159)*t11*s23r+(2048*t98*t31*t6*t331-2048*t333*t153)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((-1024*t341*t28+1024*t343*t18)*t114+((t352+512*(4-t353)*t
     #2*t119)*t11*s23r+(-2048*t362+512*(-4+t235)*t1*t57)*t11)*t66+(-1024
     #*t341*t373+(-t352+256*(t376+4)*t2*t119)*t11*t133+(-2048*t384*t349*
     #t147+t389-256*(t87+4+t156)*t1*t57)*t11*s23r+(2048*t400-2048*t151*t
     #361-256*(-4+t87-t156)*t9)*t11)*s1vr-2560*t47*t112*t373+(-t389+5120
     #*t414*t167)*t11*t133+(-2048*t420*t182+4096*t166*t423-2560*t167)*t1
     #1*s23r+(4096*t152*t40*t147-2048*t4*t1*t423)*t11)*s1l1+((-1024*t341
     #*t440-1280*t443)*t11*s23r+(t449+1280*t343*t167)*t11)*t114
      s3 = s1
      s5 = ((-1024*t341*t9*s3l1+1536*t443)*t11*t133+((1024*t351+256*(4-t
     #376)*t2*t119)*s2l1+(-1024*t348*t349*s3l1+1024*t472*t200)*s34r+256*
     #(8-t87)*t2*t195+256*(-2-9*z)*z*t484)*t11*s23r+((-1024*t362+256*(-4
     #-t235)*t1*t57)*s2l1+(2048*t198*t302-1024*t201)*s34r+256*t230*t1*t6
     #1+1024*z*t1*t507)*t11)*t66
      s6 = (-1280*t514*t340*t373+(-t449-1024*t472*t4*t423-1024*t343*t259
     #+256*(5+t123)*z*t484)*t11*t133+((1024*t531*t3*t120-512*(-t271+t535
     #-2)*t1*t57)*s2l1-1024*t384*t340*t182+(1024*t37*t4*t3*t195+1024*t98
     #*t112*t206*t9)*s34r-512*(-t353+2+t535)*t1*t61-256*(-1+t123+t305)*z
     #*t507)*t11*s23r+((-1024*t400-1024*t531*t2*t275+512*(-3+t235)*z*t9)
     #*s2l1+(-1024*t222+1024*t226)*t147+(3072*t260-1024*t514*l3*t4*t9)*s
     #34r+512*(2+t535-t262)*t9*s3l1-256*(-t305+1)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2-1280*t333*t112*t373+(-1024*t4*t98*t3*t120+2560*t605*t57
     #)*t11*t133+(-1024*t420*t611*t147+2048*t317*t361-1280*t611)*t11*s23
     #r+(2048*t320*t331-1024*t605*t361)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet2 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet3 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet4 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      C2num1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function C2num1e3(var,ndim,z,ux,l1,
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
      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = l3-1
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t18 = t9*t11
      t23 = s13r**2
      t26 = t5*t3
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t39 = z*t2
      t45 = l3*t9
      t46 = t45*s2l1
      t48 = t45*s3l1
      t58 = z*l3
      t60 = t5*t9
      t61 = t11*s23r
      t69 = z**2
      t70 = t69*t2
      t76 = t4*t3
      t77 = l3**2
      t78 = t76*t77
      t79 = t28*s1vr
      t85 = t4*t9
      t91 = s23r**2
      t92 = t18*t91
      t95 = 2*z
      t109 = t85*t11*t91
      t112 = t4*z
      t118 = t77*t9
      t119 = t118*s2l1
      t120 = t76*t119
      t122 = t77*t3
      t124 = t122*t85*s3l1
      t127 = t122*t112*t9
      t136 = t4*t2
      t150 = t136*l3
      t151 = z*t9
      t152 = t151*s2l1
      t154 = t151*s3l1
      t182 = t3*t77*l3
      t184 = t2*t77
      t191 = -4+z
      t196 = l3*t1
      t202 = t18*t91*s23r
      t206 = -6+4*z
      t228 = t9*s2l1
      t229 = t182*t228
      t231 = t9*s3l1
      t232 = t182*t231
      t234 = t182*t151
      t254 = (-4+t95)*t2
      t259 = t184*t9
      t281 = (-3+t95)*z
      s2 = (-1024*t3*t6*t9*t11*s1vr*s1l1-512*t6*z*t3*t18*s1vr)*t23*s13r+
     #((2048*t27*t18*t28+(3072*t27*t32-3072*t2*t5*t18)*s1vr+1024*t39*t5*
     #t32)*s1l1+(1024*t26*t46+1024*t26*t48+1024*l3*t3*z*t5*t9)*t11*t28+(
     #1536*t58*t3*t60*t61-1536*t39*t60*t11)*s1vr+512*t70*t5*t32)*t23
      s3 = s2+((-1024*t78*t18*t79+(-4096*t78*t32+4096*l3*t2*t85*t11)*t28
     #+(-3072*t78*t92-1024*(-6+t95)*l3*t2*t85*t61-3072*t4*t1*t18)*s1vr-2
     #048*t58*t2*t109+2048*t112*t1*t32)*s1l1+(-1024*t120-1024*t124-512*t
     #127)*t11*t79+((-2048*t120-2048*t124-2048*t127)*t11*s23r+(2048*t136
     #*t46+2048*t136*t48+2048*t136*t58*t9)*t11)*t28+(-1536*z*t77*t3*t109
     #+(-1024*t150*t152-1024*t150*t154-1024*(-3+z)*z*l3*t136*t9)*t11*s23
     #r-1536*t112*t1*t9*t11)*s1vr-1024*t69*l3*t2*t109+1024*t4*t69*t1*t32
     #)*s13r
      s1 = s3+((1024*t182*t32-1024*t184*t18)*t79+(2048*t182*t92+1024*t19
     #1*t2*t77*t32+2048*t196*t18)*t28+(1024*t182*t202+512*t206*t2*t77*t9
     #2-512*t206*t1*l3*t32-1024*t18)*s1vr+1024*t39*t77*t202-2048*t196*z*
     #t92+1024*t151*t61)*s1l1+((1024*t229+1024*t232+512*t234)*t11*s23r+(
     #-1024*t184*t228-1024*t184*t231-512*t184*t151)*t11)*t79
      taux2 = s1+((1024*t229+1024*t232+1024*t234)*t11*t91+(512*t254*t119
     #+512*t254*t118*s3l1+512*t191*z*t259)*t11*s23r+(1024*t196*t228+1024
     #*t196*t231+1024*t196*t151)*t11)*t28+(512*t182*z*t202+(1024*t39*t11
     #9+1024*t184*t154+512*t281*t259)*t11*t91+(-1024*t196*t152-1024*z*t1
     #*t48-512*t281*t196*t9)*t11*s23r-512*t151*t11)*s1vr+512*t70*t77*t20
     #2-1024*t69*t1*l3*t92+512*t69*t9*t61

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet2 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet3 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet4 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      C2num1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function C2num1e4(var,ndim,z,ux,l1,
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
      tcoff3 = taux2
      tcoff4 = -taux2

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet1 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet2 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet3 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet4 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      C2num1e4 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff4
      RETURN
      END

      double precision function C2num2e0(var,ndim,z,ux,l1,
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

      t3 = z**2
      t5 = s1l1**2
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
      taux3 = ((-256*s1l1*z-128*t3-256*t5)*t10*t13*t61/t63*t66+256*s1vr*
     #(z+2*s1l1)*t121*t66*t13*t9*t60/t62-256*t130*t131*t66*t13)/t137
      t1 = r*t
      t2 = -1+t1
      t4 = r**2
      t5 = 1/t4
      t6 = t4*r
      t7 = t6*t
      t8 = l1**2
      t9 = dsqrt(l2)
      t10 = t8*t9
      t12 = t**2
      t13 = t4*t12
      t17 = l1*t9
      t18 = t4**2
      t24 = t9*t6
      t26 = t9*t4
      t28 = t9*r
      t35 = (1+t4)**2
      t38 = t4*l1
      t39 = t38*t
      t41 = l1*r
      t43 = l1*t
      t47 = dsqrt(-l1*(l1-1)*(1-l2)*t35*r*(t4*t+t39+r+t12*r-t41-t41*t12+
     #t+t43))
      t51 = t7*t10-t10*t13-t10*t4+t1*t10+t17*t18+2*t17*t13+4*t17*t4+t17-
     #t24*t-t26*t12-t26-t28*t+2*t47-4*t47*l4
      t54 = r+t
      t56 = 1/t35
      t64 = l1*l2
      t68 = l2*t4
      t76 = l2*t8
      t81 = t47*t9
      t85 = -l1-2*t38-l1*t18+t64+4*t64*t4+t64*t18-t68-l2*t6*t-l2*r*t-t68
     #*t12+2*t64*t13-t76*t4+t76*t7+t76*t1-t76*t13-4*t81*l4+2*t81
      t91 = dsqrt((r-t)*s1l1*s2l1*t85*t2/r*t56)
      t95 = t9*s1l1
      t102 = s1l1*l1
      t114 = z**2
      t116 = t2**2
      t119 = t51**2
      t120 = t54**2
      t127 = s4vr**2
      taux4 = (-1024*s1vr*t2*t5*t51/t54*t56*z*(-t24*s1l1-2*t6*t91*c3l1+t
     #95*t39+t28*s1l1*t12-2*t91*c3l1*r-t28*t102-t28*t102*t12+t95*t43)-51
     #2*(-2*s1l1*z-t114)*t116*t5*t119/t120*t56)/t127

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet6 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet8 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      C2num2e0 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function C2num2e1(var,ndim,z,ux,l1,
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
      t20 = t12*l1
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
      t1 = r*t
      t2 = -1+t1
      t4 = r**2
      t5 = 1/t4
      t6 = s1vr*t2*t5
      t7 = t4*r
      t8 = t7*t
      t9 = l1**2
      t10 = dsqrt(l2)
      t11 = t9*t10
      t13 = t**2
      t14 = t4*t13
      t18 = l1*t10
      t19 = t4**2
      t25 = t7*t10
      t27 = t10*t4
      t29 = t10*r
      t36 = (1+t4)**2
      t39 = t4*l1
      t40 = t39*t
      t42 = l1*r
      t44 = l1*t
      t48 = dsqrt(-l1*(l1-1)*(1-l2)*t36*r*(t4*t+t40+r+r*t13-t42-t42*t13+
     #t+t44))
      t52 = t8*t11-t11*t14-t11*t4+t1*t11+t18*t19+2*t18*t14+4*t18*t4+t18-
     #t25*t-t27*t13-t27-t29*t+2*t48-4*t48*l4
      t53 = r+t
      t55 = t52/t53
      t56 = 1/t36
      t57 = l3*s1vr
      t63 = t57*z
      t79 = l1*l2
      t83 = l2*t4
      t91 = l2*t9
      t96 = t48*t10
      t100 = -l1-2*t39-l1*t19+t79+4*t79*t4+t79*t19-t83-l2*t7*t-l2*r*t-t8
     #3*t13+2*t79*t14-t91*t4+t91*t8+t91*t1-t91*t14-4*t96*l4+2*t96
      t106 = dsqrt((r-t)*s1l1*s2l1*t100*t2/r*t56)
      t110 = t10*s1l1
      t117 = s1l1*l1
      t127 = s13*s2l1
      t129 = z**2
      t135 = t2**2
      t138 = t52**2
      t139 = t53**2
      t146 = s4vr**2
      taux4 = (-512*(-t6*t55*t56*(-t57*s13+s34*s13-s13-t57+t57*s13*z+s34
     #-1+2*z+t63)-4*t6*t55*t56*z)*(-t25*s1l1-2*t7*t106*c3l1+t110*t40+t29
     #*s1l1*t13-2*t106*c3l1*r-t29*t117-t29*t117*t13+t110*t44)-512*(t57*s
     #2l1-t57*s2l1*z-t127*t63+3*t129+6*s1l1*z+t127*t57)*t135*t5*t138/t13
     #9*t56)/t146

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet6 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet8 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      C2num2e1 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function C2num2e2(var,ndim,z,ux,l1,
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
      t1 = r*t
      t2 = -1+t1
      t4 = r**2
      t5 = 1/t4
      t6 = t4*r
      t7 = t6*t
      t8 = l1**2
      t9 = dsqrt(l2)
      t10 = t8*t9
      t12 = t**2
      t13 = t4*t12
      t17 = l1*t9
      t18 = t4**2
      t24 = t9*t6
      t26 = t9*t4
      t28 = t9*r
      t35 = (1+t4)**2
      t38 = t4*l1
      t39 = t38*t
      t41 = l1*r
      t43 = l1*t
      t47 = dsqrt(-l1*(l1-1)*(1-l2)*t35*r*(t4*t+t39+r+r*t12-t41-t41*t12+
     #t+t43))
      t51 = t7*t10-t10*t13-t10*t4+t1*t10+t17*t18+2*t17*t13+4*t17*t4+t17-
     #t24*t-t26*t12-t26-t28*t+2*t47-4*t47*l4
      t54 = r+t
      t56 = 1/t35
      t58 = l3*s1vr
      t64 = t58*z
      t72 = l1*l2
      t76 = l2*t4
      t84 = l2*t8
      t89 = t47*t9
      t93 = -l1-2*t38-l1*t18+t72+4*t72*t4+t72*t18-t76-l2*t6*t-l2*r*t-t76
     #*t12+2*t72*t13-t84*t4+t84*t7+t84*t1-t84*t13-4*t89*l4+2*t89
      t99 = dsqrt((r-t)*s1l1*s2l1*t93*t2/r*t56)
      t103 = t9*s1l1
      t110 = s1l1*l1
      t124 = s13*s2l1
      t127 = z**2
      t129 = t2**2
      t132 = t51**2
      t133 = t54**2
      t140 = s4vr**2
      taux4 = (-1024*s1vr*t2*t5*t51/t54*t56*(-t58*s13+s34*s13-s13-t58+t5
     #8*s13*z+s34-1+2*z+t64)*(-t24*s1l1-2*t6*t99*c3l1+t103*t39+t28*s1l1*
     #t12-2*t99*c3l1*r-t28*t110-t28*t110*t12+t103*t43)-1024*(-t58*s2l1+t
     #58*s2l1*z-2*s1l1*z-t124*t58+t124*t64-t127)*t129*t5*t132/t133*t56)/
     #t140

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet6 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet8 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      C2num2e2 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function C2num2e3(var,ndim,z,ux,l1,
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

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet6 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet8 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      C2num2e3 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

      double precision function C2num2e4(var,ndim,z,ux,l1,
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

      tcoff5 = taux3
      tcoff6 = taux3
      tcoff7 = taux4
      tcoff8 = -taux4

CCC Jet functions

      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet5 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet6 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet7 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = 1-z
      t2 = l3*t1
      t3 = t2*s23r
      t4 = -l3+1
      t6 = t4*t1*s13r
      t7 = t2*s1vr
      t9 = t1**2
      t11 = l3*t4*t9*s34r
      tFjet8 = Fjet(var,ndim,z,ux,t3,t6,t7-t3-t11,-1-t6-t7+z,t11,-1+z-t7
     #+t11,t7,-t6-t3-t11+z,1+t6+t3,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff8)

CCC Sum

      C2num2e4 = tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8
      RETURN
      END

