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
      t59 = t58*s2l1
      t62 = t58*s3l1
      t65 = 2*z
      t66 = t65+1
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
      t101 = t1*t10
      t104 = 768*t97-256*t99*t5*t101
      t108 = s23r**2
      t109 = t12*t108
      t113 = t5*t1
      t114 = t113*t10
      t119 = l3*t3
      t120 = t6**2
      t122 = s34r**2
      t131 = t95*s2l1
      t132 = t86*t131
      t134 = t95*s3l1
      t135 = t86*t134
      t136 = 256*t135
      t137 = l3*z
      t139 = t86*t137*t10
      t157 = t10*s3l1
      t160 = z+2
      t168 = t137*t6
      t169 = t11*t108
      t170 = t101*t169
      t173 = t38*s34r
      t174 = t53*t173
      t176 = z*t5
      t183 = z*t10
      t192 = l3**2
      t193 = t192*t7
      t194 = t193*t1
      t195 = t88*s1vr
      t203 = s2l1**2
      t204 = t95*t203
      t210 = t192*t120
      t211 = t210*z
      t213 = 128*t211*t173
      t214 = t75*l3
      t216 = t214*t86*t10
      t217 = 256*t216
      t229 = t87*t183*s3l1
      t245 = z*t120*t5
      t247 = t192*t10
      t251 = t76*t7
      t252 = t251*t96
      t254 = s3l1**2
      t258 = t75*t5
      t262 = 2*t75
      t263 = 1+t262+t65
      t271 = t214*t6
      t275 = t1*t75
      t284 = t95*t122
      t293 = t6*t192
      t294 = t293*t10
      t298 = t1*t7
      t299 = t247*s34r
      t300 = t298*t299
      t301 = l3*t5
      t302 = t301*t183
      t304 = -256*t300+256*t302
      t307 = t86*t96
      t309 = 256*t307-256*t183
      t314 = t12*t108*s23r
      t317 = t301*t10
      t322 = t38*t122
      t339 = z*t192
      t340 = t339*t58
      t382 = -t5
      t396 = t183*s34r
      t397 = t194*t396
      t402 = t211*t322
      t419 = t192*l3
      t420 = z*t7*t419
      t437 = 128*t1*t120*t419*t396
      t438 = t10*t254
      t457 = t258*t95
      t475 = -t99*z
      t491 = t382*z
      t536 = t101*s34r
      s3 = -64*t8*t12*s1vr*t13*s13r
      s5 = ((-256*t19*t12*s1vr+128*t19*t23+256*t3*t7*t28)*t32+((-256*t19
     #*t34-384*t37*t38)*t11*s1vr+128*t44*t6*t23+256*t8*t28)*s1l1+(64*t53
     #*t38*t54+(-128*t44*t59-128*t44*t62-64*t66*t2*t37*t10)*t11)*s1vr+64
     #*t76*t6*t23+128*t3*t75*t7*t28)*t13
      s8 = (-256*t87*t12*t88+(256*t87*t23+t104*t11)*s1vr-256*t87*t109+(-
     #512*t97+256*t114)*t11*s23r+(-512*t119*t120*t10*t122+512*t19*t27)*t
     #11)*t32+((-512*t132+t136-384*t139)*t11*t88+((256*t132-t136+512*t13
     #9)*t11*s23r+(t104*s2l1+(-256*t94*t134+768*t52*t44*t10)*s34r-256*t8
     #6*t157-256*t160*z*t114)*t11)*s1vr-256*t168*t170+(-512*t174+256*t17
     #6*t101)*t11*s23r+(-512*t119*t120*t183*t122+512*t37*t173)*t11)*s1l1
     #-64*t194*t183*t11*t195+(128*t194*t183*t54+(-256*t86*t204+(256*t135
     #-256*t139)*s2l1+t213-t217)*t11)*t88
      s7 = s8+(-64*t194*t183*t169+(128*t168*t101*s2l1-t213+128*t229+t217
     #)*t11*s23r+(-128*t176*t101*t203+(256*t174-128*t160*z*t114)*s2l1-12
     #8*t245*t3*t247*t122+384*t252-128*t176*t101*t254-128*t258*t101*s3l1
     #-64*t263*z*t114)*t11)*s1vr-128*t271*t170+(-256*t252+128*t275*t5*t1
     #0)*t11*s23r+(-256*t75*t120*t3*t284+256*t76*t58*s34r)*t11
      s8 = s13r
      s6 = s7*s8
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(128*t294*t54*t88+(t304*t11*s23r+t309*t11)*s1vr+128*t293*t
     #314+(256*t300-256*t317)*t11*t108+(256*t210*t322-512*t307+128*t10)*
     #t11*s23r+(-512*t94*t284+256*t113*t27)*t11)*t32+(((256*t293*t34+256
     #*t340)*t11*s23r+(256*t301*t157-128*t302)*t11)*t88+((256*t293*t157-
     #128*t340)*t11*t108+(t304*s2l1+(256*t298*t247*s3l1-256*t193*t1*z*t1
     #0)*s34r+256*(z-2)*t5*t134+256*t99*z*t317)*t11*s23r+(t309*s2l1+(-76
     #8*t135+256*t139)*s34r+256*t382*t10*s3l1-128*t66*z*t10)*t11)*s1vr+1
     #28*t37*t192*t314+(256*t397-256*t302)*t11*t108+(256*t402-512*t87*t3
     #96+128*t183)*t11*s23r+(-512*t53*t322+256*t113*t396)*t11)*s1l1+(64*
     #t420*t23-64*t339*t58*t11)*t195
      taux1 = s1+(-128*t420*t109+(128*t293*t10*t203+128*t339*t59-t437+12
     #8*t293*t438+128*t339*t62+64*(3*z+2)*z*t294)*t11*s23r+(256*t301*t15
     #7*s2l1+128*t397-256*t301*t438-128*t457)*t11)*t88+(64*t420*t314+(t4
     #37-64*t66*z*t294)*t11*t108+(128*t176*t204-128*t475*t5*t131+128*t24
     #5*t2*t419*t10*t122+128*t475*t192*t298*t27+128*t176*t95*t254-128*t4
     #91*t5*t134+64*t263*z*t317)*t11*s23r+(-128*t183*t203+128*t475*t34-1
     #28*t402+(-256*t229+128*t216)*s34r-128*t183*t254+128*t491*t157-64*(
     #t262+1)*z*t10)*t11)*s1vr+64*t75*t6*t192*t314+(128*t275*t7*t299-128
     #*t457)*t11*t108+(128*t210*t75*t322-256*t271*t536+64*t75*t10)*t11*s
     #23r+(-256*t251*t284+128*t258*t536)*t11
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
      t69 = t1*t6
      t70 = t69*l3
      t71 = z*t10
      t72 = s1vr**2
      t77 = t71*t36
      t80 = t35*s34r
      t81 = t34*t80
      t83 = t1*t58
      t84 = t5*t10
      t85 = t83*t84
      t92 = z*l3*t6
      t93 = t1*t10
      t94 = s23r**2
      t95 = t11*t94
      t96 = t93*t95
      t100 = -t1
      t108 = t6**2
      t110 = s34r**2
      t119 = l3**2
      t121 = t119*t7*t1
      t122 = t72*s1vr
      t130 = t92*t93*s2l1
      t132 = t119*t108
      t133 = t132*z
      t135 = 512*t133*t80
      t136 = t58*l3
      t138 = t136*t69*t10
      t139 = 1024*t138
      t149 = t70*t71*s3l1
      t159 = z*t108*t5
      t161 = t119*t10
      t165 = t59*t7
      t166 = l3*t10
      t168 = t165*t166*s34r
      t175 = 2*t58
      t176 = 1-t175-t48
      t186 = t136*t6
      t198 = t166*t110
      t207 = z*t119
      t208 = t207*t6
      t211 = l3*t5
      t217 = t12*t94
      t220 = t71*s34r
      t221 = t121*t220
      t222 = t58*t5
      t223 = t222*t166
      t224 = t221-t223
      t228 = t70*t220
      t239 = t12*t94*s23r
      t247 = t35*t110
      t248 = t133*t247
      t265 = t119*l3
      t266 = z*t7*t265
      t280 = 512*t1*t108*t265*t220
      t287 = t6*t119*t10
      t307 = t58-z
      t308 = t307*t5
      t371 = t93*s34r
      s2 = -256*t8*t12*s1vr*t13*s13r+((512*t20*t12*s1vr-512*t20*t24-1024
     #*t8*t28)*s1l1+(256*t34*t35*t36+(512*t40*t42-512*t40*t45-256*(-t48+
     #1)*t2*t19*t10)*t11)*s1vr-256*t59*t6*t24-512*t3*t58*t7*t28)*t13
      s3 = s2+((1536*t70*t71*t11*t72+(-1024*t70*t77+(-3072*t81+1024*t85)
     #*t11)*s1vr+1024*t92*t96+(2048*t81+1024*t100*z*t84)*t11*s23r+(2048*
     #l3*t3*t108*t71*t110-2048*t19*t80)*t11)*s1l1+256*t121*t71*t11*t122+
     #(-512*t121*t77+(1024*t130-t135+t139)*t11)*t72+(256*t121*t71*t95+(-
     #512*t130+t135+512*t149-t139)*t11*s23r+((-1024*t81+512*t85)*s2l1+51
     #2*t159*t3*t161*t110-1536*t168+512*(t58-t48)*t1*t84*s3l1-256*t176*t
     #1*z*t5*t10)*t11)*s1vr+512*t186*t96+(1024*t168+512*t100*t58*t84)*t1
     #1*s23r+(1024*t58*t108*t3*t198-1024*t59*t41*s34r)*t11)*s13r
      s1 = s3+((-1024*t208*t24+512*t211*t71*t11)*t72+(512*t208*t217+1024
     #*t224*t11*s23r+(-1024*t228-256*(t48-4*t58)*t10)*t11)*s1vr-512*t208
     #*t239+(-1024*t221+1024*t211*t71)*t11*t94+(-1024*t248+2048*t228-512
     #*t71)*t11*s23r+(2048*t34*t247-1024*t1*t5*t220)*t11)*s1l1+(-256*t26
     #6*t24+256*t207*t41*t11)*t122
      taux2 = s1+(512*t266*t217+(-512*t207*t42+t280-512*t207*t45+256*(-3
     #*z-2)*z*t287)*t11*s23r-512*t224*t11)*t72+(-256*t266*t239+(-t280+25
     #6*(t48+1)*z*t287)*t11*t94+(-512*t308*t166*s2l1-512*t159*t2*t265*t1
     #0*t110+512*(z+1)*z*t119*t7*t1*t27-512*t308*t166*s3l1+256*t176*z*t2
     #11*t10)*t11*s23r+(512*z*t5*t10*s2l1+512*t248+(1024*t149-512*t138)*
     #s34r+512*t307*t10*s3l1-256*(1-t175)*z*t10)*t11)*s1vr-256*t58*t6*t1
     #19*t239+(-512*t83*t7*t161*s34r+512*t223)*t11*t94+(-512*t132*t58*t2
     #47+1024*t186*t371-256*t58*t10)*t11*s23r+(1024*t165*t198-512*t222*t
     #371)*t11

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
      t42 = l3*z
      t43 = t6**2
      t52 = t6*z
      t59 = z**2
      t60 = t2*t59
      t62 = t12*s23r
      t74 = l3*t6*t1
      t75 = s23r*s1vr
      t78 = s23r**2
      t79 = t12*t78
      t81 = t5*t1
      t86 = t1*t6
      t87 = l3*t10
      t98 = t42*t6
      t99 = t1*t10
      t100 = t11*t78
      t101 = t99*t100
      t104 = z*t5
      t110 = l3**2
      t112 = t110*t7*t1
      t113 = z*t10
      t114 = s1vr**2
      t115 = t114*s1vr
      t120 = t113*t21
      t124 = t98*t99*s2l1
      t125 = t110*t43
      t126 = t125*z
      t128 = t126*t34*s34r
      t129 = t113*s3l1
      t130 = t74*t129
      t131 = 2*z
      t135 = t86*t10
      t156 = s2l1**2
      t161 = t104*t99*s3l1
      t167 = z*t43*t5
      t169 = t110*t10
      t170 = s34r**2
      t174 = t60*t7
      t176 = t174*t87*s34r
      t178 = s3l1**2
      t194 = t59*l3*t6
      t209 = t87*t170
      t218 = t6*t110
      t219 = t218*t10
      t224 = t5*l3
      t230 = t12*t78*s23r
      t238 = t10*s2l1
      t239 = t218*t238
      t241 = t10*s3l1
      t242 = t218*t241
      t244 = z*t110
      t245 = t244*t38
      t257 = t224*t113
      t274 = t110*l3
      t275 = z*t7*t274
      t292 = t113*s34r
      t294 = 128*t1*t43*t274*t292
      t298 = t59*t6
      t307 = t224*t129
      t340 = 128*t1*t59*t7*t169*s34r
      t360 = t34*t170
      t363 = t99*s34r
      t381 = t59*t5
      s3 = 128*t3*z*t7*t12*s1vr*t13*s13r
      s5 = (-128*t2*t6*t10*t21*t22-128*t26*t6*t12*s23r*s1l1+(-256*t33*t3
     #4*t21+(128*t26*t38*s2l1-128*t42*t43*t3*t10*s34r+128*t26*t38*s3l1+2
     #56*t52*t34)*t11)*s1vr-128*t60*t6*t62-128*t3*t59*t7*t10*s34r*t11)*t
     #13
      s6 = ((256*t74*t12*t75+256*t74*t79-256*t81*t62)*t22+((256*t86*t87*
     #s2l1+256*t86*t87*s3l1+256*t86*t42*t10)*t11*t75+256*t98*t101-256*t1
     #04*t1*t62)*s1l1+128*t112*t113*t11*t115+(-128*t112*t120+(-128*t124-
     #128*t128-128*t130+128*(-1+t131)*z*l3*t135)*t11)*t114+(256*t112*t11
     #3*t100+(-128*t124+256*t128-128*t130-128*(z+2)*z*l3*t135)*t11*s23r+
     #(128*t104*t99*t156+(256*t161+256*t104*t99)*s2l1+128*t167*t3*t169*t
     #170-384*t176+128*t104*t99*t178+256*t161-64*(-2*t59+t131-4)*z*t81*t
     #10)*t11)*s1vr+256*t194*t101+(256*t176-256*t1*t59*t5*t10)*t11*s23r+
     #(256*t59*t43*t3*t209-256*t60*t38*s34r)*t11)*s13r
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-128*t219*t21*t114+(-256*t218*t79+256*t224*t62)*s1vr-128*
     #t218*t230+256*t224*t79-128*t62)*t22+((-256*t239-256*t242-128*t245)
     #*t11*s23r*t114+((-256*t239-256*t242-256*t245)*t11*t78+(256*t224*t2
     #38+256*t224*t241+256*t257)*t11*s23r)*s1vr-128*t52*t110*t230+256*t2
     #24*z*t79-128*t120)*s1l1+(-128*t275*t62+128*t244*t38*t11)*t115
      taux1 = s1+(128*t275*t79+(-128*t218*t10*t156-256*t218*t241*s2l1+t2
     #94-128*t218*t10*t178-256*t298*t169)*t11*s23r+(-128*t224*t113*s2l1-
     #128*t112*t292-128*t307+128*t42*t38)*t11)*t114+(-128*t275*t230+(-t2
     #94+64*(2*z+2)*z*t219)*t11*t78+(-128*t104*t87*t156+(-256*t307-128*t
     #257)*s2l1-128*t167*t2*t274*t10*t170+t340-128*t104*t87*t178-128*t30
     #7-64*(2*t59+2)*z*t224*t10)*t11*s23r+(128*t113*t156+(256*t129+128*t
     #113)*s2l1+128*t126*t360-128*t33*t363+128*t113*t178+128*t129-64*(2*
     #z-2*t59-2)*z*t10)*t11)*s1vr-128*t298*t110*t230+(-t340+256*t381*t87
     #)*t11*t78+(-128*t125*t59*t360+256*t194*t363-128*t59*t10)*t11*s23r+
     #(256*t174*t209-128*t381*t363)*t11
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t17 = t3*z*t6
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t35 = t5**2
      t37 = l3*t9
      t38 = t37*s34r
      t40 = 1024*t3*t35*t38
      t41 = 6*z
      t42 = t41-1
      t44 = t5*t9
      t51 = z*t5
      t56 = t9*s34r*t11
      t64 = t37*s3l1
      t67 = l3*t6
      t68 = t2*z
      t69 = t68*t9
      t75 = t67*z
      t76 = t2*t9
      t81 = 3*z
      t88 = l3*z
      t94 = t44*s3l1
      t97 = t41-7
      t106 = z**2
      t107 = t2*t106
      t117 = t1*t6
      t118 = l3**2
      t119 = t117*t118
      t120 = t28*s1vr
      t124 = -t1
      t125 = t124*t6
      t129 = t2*t35
      t130 = t118*t9
      t131 = t130*s34r
      t132 = t129*t131
      t134 = t5*t1
      t142 = s23r**2
      t143 = t18*t142
      t147 = 8*z
      t150 = t134*t9
      t156 = t35*t4
      t157 = t3*t156
      t158 = s34r**2
      t160 = t157*t130*t158
      t162 = t76*s34r
      t163 = t75*t162
      t165 = 12*t106
      t168 = t4*t9
      t176 = t1*t9
      t177 = t11*t142
      t178 = t176*t177
      t182 = 4*t1
      t191 = z*t9
      t200 = t130*s2l1
      t203 = 3*t124
      t205 = t130*s3l1
      t211 = t118*t6*t1*z*t9
      t229 = t129*t205
      t230 = t118*t35
      t231 = t230*t69
      t235 = t134*l3
      t237 = t235*t191*s3l1
      t239 = t106*l3
      t256 = t230*z
      t260 = 7*z
      t273 = 3*t106
      t283 = z*t156
      t303 = 6*t106
      t314 = t239*t5
      t317 = t107*t6
      t328 = t37*t158
      t337 = t118*l3
      t338 = t6*t337
      t340 = t5*t118
      t347 = t124*t35
      t348 = t337*t9
      t350 = t347*t348*s34r
      t351 = 1024*t350
      t352 = z*t118
      t353 = t352*t44
      t358 = t125*t131
      t360 = 1-t81
      t369 = t18*t142*s23r
      t372 = 10*z
      t381 = t348*t158
      t384 = t191*s34r
      t386 = 3072*t119*t384
      t395 = t76*t158
      t396 = t230*t395
      t400 = t176*s34r
      t413 = l3*t4
      t436 = t9*s2l1
      t438 = 256*t338*t436
      t439 = t9*s3l1
      t441 = 256*t338*t439
      t442 = t338*t191
      t467 = t1*t35
      t474 = 5*z
      t481 = t340*t9
      t488 = z+1
      t505 = t413*t9
      t516 = 2-z
      t532 = t516*t1
      t609 = t106*t4
      s2 = (-256*t3*t6*t9*t11*s1vr*s1l1+768*t17*t18*s1vr)*t23*s13r+((-51
     #2*t27*t18*t28+(768*t27*t32+(t40-256*t42*t2*t44)*t11)*s1vr+2048*t51
     #*t2*t32+3072*t17*t56)*s1l1+(768*t26*t37*s2l1+256*t26*t64-256*t67*t
     #69)*t11*t28+(-512*t75*t76*t11*s23r+((-t40-256*(t81-2)*t2*t44)*s2l1
     #+512*t88*t35*t3*t9*s34r+1792*t68*t94-256*t97*t2*t51*t9)*t11)*s1vr+
     #1024*t107*t5*t32+1536*t3*t106*t6*t56)*t23
      s4 = s2
      s7 = (-256*t119*t18*t120+(-1024*t125*t118*t32+(1024*t132-4608*t134
     #*t88*t9)*t11)*t28+(-768*t119*t143+(-2048*t132-512*(1-t147)*l3*t150
     #)*t11*s23r+(-1024*t160+9216*t163-256*(-1+t165)*t1*t168)*t11)*s1vr-
     #4096*t88*t5*t178+(-6144*t163+1024*t182*z*t168)*t11*s23r+(-6144*l3*
     #t3*t35*t191*t158+6144*t51*t162)*t11)*s1l1+(-256*t117*t200-256*t203
     #*t6*t205-1024*t211)*t11*t120+((-512*t117*t205+2304*t211)*t11*s23r+
     #(256*(-11*z+2)*l3*t5*t176*s2l1+(-1536*t229+1536*t231)*s34r+768*t23
     #7-3072*t239*t150)*t11)*t28
      s6 = s7+(-1280*t119*t191*t177+((512*t132+1024*l3*t4*t150)*s2l1-204
     #8*t256*t162-2048*t237-512*(-t260+2)*z*l3*t150)*t11*s23r+((512*t160
     #-512*(-t260+4)*t2*t6*t38+512*(1-t273+2*z)*t1*t168)*s2l1+(512*t157*
     #t205-1536*t283*t3*t118*t9)*t158+(-512*t75*t76*s3l1+4608*t107*t67*t
     #9)*s34r+512*(-1-t273+t147)*t1*t168*s3l1-256*(t41+t303-7)*t1*z*t4*t
     #9)*t11)*s1vr-2048*t314*t178+(-3072*t317*t38+512*t182*t106*t168)*t1
     #1*s23r+(-3072*t106*t35*t3*t328+3072*t107*t44*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((256*t338*t32-256*t340*t18)*t120+(-512*t338*t143+(t351+35
     #84*t353)*t11*s23r+(-1024*t358+512*t360*t4*t37)*t11)*t28+(256*t338*
     #t369+(-t351+256*(-t372+1)*t5*t130)*t11*t142+(1024*t2*t156*t381-t38
     #6-256*(1-t165-4*z)*t4*t37)*t11*s23r+(-1024*t396-1024*t360*l3*t5*t4
     #00-256*(1-t41+t165)*t9)*t11)*s1vr+2048*t51*t118*t369+(t386-4096*t4
     #13*t191)*t11*t142+(3072*t256*t395-6144*t235*t384+2048*t191)*t11*s2
     #3r+(-6144*t75*t395+3072*t1*t4*t384)*t11)*s1l1+((t438-t441+1024*t44
     #2)*t11*s23r+(-256*t340*t436+256*t340*t439-1024*t353)*t11)*t120
      s3 = s1
      s5 = ((-t438+t441-2048*t442)*t11*t142+((512*t350+1792*t353)*s2l1+(
     #-512*t347*t348*s3l1-1536*t467*t337*z*t9)*s34r+256*(-2+t474)*t5*t20
     #5+256*(t372+7)*z*t481)*t11*s23r+((-512*t358+256*t488*t4*t37)*s2l1+
     #(1024*t125*t205+1536*t211)*s34r+256*(1+t81)*t4*t64-256*t42*z*t505)
     #*t11)*t28
      s6 = (1024*z*t6*t337*t369+(256*t516*t5*t200+1536*t467*t337*t384+25
     #6*t352*t94+256*(-t147-3)*z*t481)*t11*t142+((512*t532*t6*t131-512*(
     #t81+1-t273)*t4*t37)*s2l1+1536*t283*t2*t381+(512*t488*t1*t6*t205+51
     #2*(-2-t81)*z*t118*t117*t9)*s34r-512*(-1+t474-t273)*t4*t64-256*(-t3
     #03-t147+4)*z*t505)*t11*s23r+((512*t396-512*t532*t5*t38-256*t97*z*t
     #9)*s2l1+(512*t229-1536*t231)*t158+(-3584*t237-512*t360*z*l3*t150)*
     #s34r+256*(-2+9*z-t303)*t9*s3l1-256*(-3+t303)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+1024*t106*t5*t118*t369+(-512*t203*t106*t6*t131-2048*t60
     #9*t37)*t11*t142+(1536*t230*t106*t395-3072*t314*t400+1024*t106*t9)*
     #t11*s23r+(-3072*t317*t328+1536*t609*t400)*t11

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
      t21 = l3*t10
      t22 = s1vr**2
      t29 = t2*t10
      t30 = t11*s23r
      t33 = t6*z
      t39 = z**2
      t42 = t12*s23r
      t47 = l3**2
      t49 = z*t47*t7
      t50 = t1*t10
      t51 = t22*s1vr
      t56 = t50*t30
      t60 = z*t10
      t61 = t60*t11
      t68 = s23r**2
      t69 = t11*t68
      t98 = z*t7*t47*l3
      t100 = t47*t10
      t106 = t12*t68
      t121 = t12*t68*s23r
      t125 = (2*z-3)*z
      taux1 = -64*t2*t1*z*t7*t12*s1vr*t13*s13r+(128*t2*t7*z*t21*t11*t22+
     #(192*l3*t7*z*t29*t30-192*t33*t29*t11)*s1vr+64*t2*t39*t6*t42)*t13+(
     #-64*t49*t50*t11*t51+(-256*t49*t56+256*t1*t6*l3*t61)*t22+(-192*t47*
     #t7*t1*t60*t69-128*(z-3)*z*l3*t6*t56-192*z*t5*t50*t11)*s1vr-128*t39
     #*l3*t6*t50*t69+128*t1*t39*t5*t42)*s13r+(64*t98*t42-64*t33*t100*t11
     #)*t51+(128*t98*t106+64*(-4+z)*z*t6*t100*t30+128*l3*t5*t61)*t22+(64
     #*t98*t121+64*t125*t6*t100*t69-64*t125*t5*t21*t30-64*t61)*s1vr+64*t
     #39*t6*t47*t121-128*t39*t5*l3*t106+64*t39*t10*t30
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t17 = t3*z*t6
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t28 = t18*s23r
      t31 = t5**2
      t33 = l3*t9
      t34 = t33*s34r
      t36 = 2048*t3*t31*t34
      t37 = -1-z
      t40 = t5*t9
      t47 = t5*z
      t52 = t9*s34r*t11
      t60 = t33*s3l1
      t65 = s1vr**2
      t67 = l3*t6
      t68 = t67*z
      t69 = t2*t9
      t80 = l3*z
      t86 = t2*z
      t87 = t40*s3l1
      t90 = 4*z
      t100 = z**2
      t101 = t2*t100
      t111 = -t1
      t112 = 4*t111
      t113 = t112*t6
      t114 = l3**2
      t115 = t113*t114
      t116 = t65*s1vr
      t120 = t2*t31
      t121 = t114*t9
      t122 = t121*s34r
      t123 = t120*t122
      t125 = 6*z
      t128 = t5*t1
      t129 = t128*t9
      t135 = s23r**2
      t140 = t125+4
      t147 = t31*t4
      t148 = t3*t147
      t149 = s34r**2
      t151 = t148*t121*t149
      t153 = t69*s34r
      t154 = t68*t153
      t156 = 8*t100
      t159 = t4*t9
      t167 = t1*t9
      t168 = t11*t135
      t169 = t167*t168
      t173 = 5*t111
      t182 = z*t9
      t191 = t121*s2l1
      t196 = t121*s3l1
      t199 = t114*t6
      t202 = t199*t1*z*t9
      t207 = t1*t6
      t222 = t120*t196
      t224 = t114*t31
      t226 = t224*t86*t9
      t230 = -t140
      t235 = 2*z
      t245 = t199*t1
      t256 = t224*z
      t259 = t128*l3
      t261 = t259*t182*s3l1
      t263 = 5*z
      t272 = 3*z
      t284 = z*t147
      t304 = 4*t100
      t316 = t100*l3*t5
      t319 = t101*t6
      t330 = t33*t149
      t339 = t114*l3
      t340 = t6*t339
      t342 = t5*t114
      t347 = t1*t31
      t348 = t339*t9
      t350 = t347*t348*s34r
      t351 = 2048*t350
      t352 = 7*z
      t360 = t207*t122
      t371 = t18*t135*s23r
      t374 = 8*z
      t383 = t348*t149
      t386 = t182*s34r
      t388 = 2048*t245*t386
      t396 = t69*t149
      t397 = t224*t396
      t399 = t167*s34r
      t412 = l3*t4
      t435 = t9*s2l1
      t438 = t340*t182
      t444 = 1024*t342*t435
      t445 = z*t114
      t481 = t342*t9
      t504 = t412*t9
      t511 = t6*z
      t528 = (z-2)*t1
      t532 = 2*t100
      t603 = t100*t4
      s2 = (1024*t3*t6*t9*t11*s1vr*s1l1-256*t17*t18*s1vr)*t23*s13r+(((-3
     #072*t26*l3*t28+(-t36-1024*t37*t2*t40)*t11)*s1vr-2560*t47*t2*t28-20
     #48*t17*t52)*s1l1+(-2048*t26*t33*s2l1-1024*t26*t60)*t11*t65+(-768*t
     #68*t69*t11*s23r+((t36-512*(z+2)*t2*t40)*s2l1-1024*t80*t31*t3*t9*s3
     #4r-1536*t86*t87-256*(7-t90)*t2*t47*t9)*t11)*s1vr-1280*t101*t5*t28-
     #1024*t3*t100*t6*t52)*t23
      s4 = s2
      s7 = (-256*t115*t18*t116+(-2048*t123+512*(-4+t125)*l3*t129)*t11*t6
     #5+(-768*t115*t18*t135+(4096*t123-512*t140*l3*t129)*t11*s23r+(2048*
     #t151-6144*t154-256*(-4-t156)*t1*t159)*t11)*s1vr+5120*t80*t5*t169+(
     #4096*t154+1024*t173*z*t159)*t11*s23r+(4096*l3*t3*t31*t182*t149-409
     #6*t47*t153)*t11)*s1l1+(-256*t113*t191+256*t112*t6*t196+1280*t202)*
     #t11*t116+((1024*t207*t191-512*t113*t196-1536*t202)*t11*s23r+(256*(
     #-8+t125)*l3*t5*t167*s2l1+(3072*t222-1024*t226)*s34r+256*t230*t1*t5
     #*t60+1024*(-1+t235)*z*l3*t129)*t11)*t65
      s6 = s7+(2304*t245*t182*t168+((-1024*t123-512*(-z-4)*l3*t129)*s2l1
     #+2048*t256*t153+2560*t261-512*(t263-1)*z*l3*t129)*t11*s23r+((-1024
     #*t151-1024*(t272-4)*t2*t6*t34+1024*(-1+t100-t235)*t1*t159)*s2l1+(-
     #1024*t148*t196+1024*t284*t3*t114*t9)*t149+(1024*t68*t69*s3l1-3072*
     #t101*t67*t9)*s34r+1024*(-t90+1+t100)*t1*t159*s3l1-256*(7-t90-t304)
     #*t1*z*t4*t9)*t11)*s1vr+2560*t316*t169+(2048*t319*t34+512*t173*t100
     #*t159)*t11*s23r+(2048*t100*t31*t3*t330-2048*t101*t40*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((-1024*t340*t28+1024*t342*t18)*t116+((t351+512*(-t352+4)*
     #t5*t121)*t11*s23r+(-2048*t360+512*(-4+t235)*t4*t33)*t11)*t65+(-102
     #4*t340*t371+(-t351+256*(4+t374)*t5*t121)*t11*t135+(-2048*t2*t147*t
     #383+t388-256*(4+t90+t156)*t4*t33)*t11*s23r+(2048*t397-2048*t67*t39
     #9-256*(-4+t90-t156)*t9)*t11)*s1vr-2560*t47*t114*t371+(-t388+5120*t
     #412*t182)*t11*t135+(-2048*t256*t396+4096*t259*t386-2560*t182)*t11*
     #s23r+(4096*t68*t396-2048*t1*t4*t386)*t11)*s1l1+((-1024*t340*t435-1
     #280*t438)*t11*s23r+(t444+1280*t445*t40)*t11)*t116
      s3 = s1
      s5 = ((-1024*t340*t9*s3l1+1536*t438)*t11*t135+((1024*t350+256*(4-t
     #374)*t5*t121)*s2l1+(-1024*t347*t348*s3l1+1024*t347*t339*z*t9)*s34r
     #+256*(8-t90)*t5*t196+256*(-2-9*z)*z*t481)*t11*s23r+((-1024*t360+25
     #6*(-4-t235)*t4*t33)*s2l1+(2048*t207*t196-1024*t202)*s34r+256*t230*
     #t4*t60+1024*z*t4*t504)*t11)*t65
      s6 = (-1280*t511*t339*t371+(-t444-1024*t347*t339*t386-1024*t445*t8
     #7+256*(t125+5)*z*t481)*t11*t135+((1024*t528*t6*t122-512*(-t272-2+t
     #532)*t4*t33)*s2l1-1024*t284*t2*t383+(1024*t37*t1*t6*t196+1024*t100
     #*t114*t207*t9)*s34r-512*(-t352+2+t532)*t4*t60-256*(t125-1+t304)*z*
     #t504)*t11*s23r+((-1024*t397-1024*t528*t5*t34+512*(t235-3)*z*t9)*s2
     #l1+(-1024*t222+1024*t226)*t149+(3072*t261-1024*t511*l3*t1*t9)*s34r
     #+512*(t532+2-t263)*t9*s3l1-256*(1-t304)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2-1280*t100*t5*t114*t371+(-1024*t1*t100*t6*t122+2560*t603
     #*t33)*t11*t135+(-1024*t224*t100*t396+2048*t316*t399-1280*t100*t9)*
     #t11*s23r+(2048*t319*t330-1024*t603*t399)*t11

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
      t39 = t5*z
      t45 = l3*t9
      t46 = t45*s2l1
      t48 = t45*s3l1
      t50 = l3*t6
      t59 = t2*t9
      t60 = t11*s23r
      t68 = z**2
      t75 = t1*t6
      t76 = l3**2
      t77 = t75*t76
      t78 = t28*s1vr
      t84 = t1*t9
      t90 = s23r**2
      t91 = t18*t90
      t94 = 2*z
      t106 = z*l3
      t108 = t11*t90
      t109 = t84*t108
      t112 = t1*z
      t118 = t76*t9
      t119 = t118*s2l1
      t120 = t75*t119
      t122 = t118*s3l1
      t123 = t75*t122
      t127 = t76*t6*t112*t9
      t136 = t1*t5
      t146 = z*t9
      t150 = t136*l3
      t151 = t146*s2l1
      t180 = t76*l3
      t181 = t6*t180
      t183 = t5*t76
      t190 = -4+z
      t195 = t4*l3
      t201 = t18*t90*s23r
      t205 = 4*z-6
      t227 = t9*s2l1
      t228 = t181*t227
      t230 = t9*s3l1
      t231 = t181*t230
      t233 = t181*t146
      t242 = z*t76
      t243 = t5*t9
      t255 = (-4+t94)*t5
      t259 = t183*t9
      t283 = (t94-3)*z
      s2 = (-1024*t3*t6*t9*t11*s1vr*s1l1-512*t3*z*t6*t18*s1vr)*t23*s13r+
     #((2048*t27*t18*t28+(3072*t27*t32-3072*t2*t5*t18)*s1vr+1024*t39*t2*
     #t32)*s1l1+(1024*t26*t46+1024*t26*t48+1024*t50*z*t2*t9)*t11*t28+(15
     #36*t50*z*t59*t60-1536*t39*t59*t11)*s1vr+512*t2*t68*t5*t32)*t23
      s3 = s2+((-1024*t77*t18*t78+(-4096*t77*t32+4096*l3*t5*t84*t11)*t28
     #+(-3072*t77*t91-1024*(-6+t94)*l3*t5*t84*t60-3072*t1*t4*t18)*s1vr-2
     #048*t106*t5*t109+2048*t112*t4*t32)*s1l1+(-1024*t120-1024*t123-512*
     #t127)*t11*t78+((-2048*t120-2048*t123-2048*t127)*t11*s23r+(2048*t13
     #6*t46+2048*t136*t48+2048*t136*t106*t9)*t11)*t28+(-1536*t77*t146*t1
     #08+(-1024*t150*t151-1024*t150*t146*s3l1-1024*(-3+z)*z*l3*t136*t9)*
     #t11*s23r-1536*t112*t4*t9*t11)*s1vr-1024*t68*l3*t5*t109+1024*t1*t68
     #*t4*t32)*s13r
      s1 = s3+((1024*t181*t32-1024*t183*t18)*t78+(2048*t181*t91+1024*t19
     #0*t5*t76*t32+2048*t195*t18)*t28+(1024*t181*t201+512*t205*t5*t76*t9
     #1-512*t205*t4*l3*t32-1024*t18)*s1vr+1024*t39*t76*t201-2048*t195*z*
     #t91+1024*t146*t60)*s1l1+((1024*t228+1024*t231+512*t233)*t11*s23r+(
     #-1024*t183*t227-1024*t183*t230-512*t242*t243)*t11)*t78
      taux2 = s1+((1024*t228+1024*t231+1024*t233)*t11*t90+(512*t255*t119
     #+512*t255*t122+512*t190*z*t259)*t11*s23r+(1024*t195*t227+1024*t195
     #*t230+1024*t195*t146)*t11)*t28+(512*z*t6*t180*t201+(1024*t39*t119+
     #1024*t242*t243*s3l1+512*t283*t259)*t11*t90+(-1024*t195*t151-1024*z
     #*t4*t48-512*t283*t195*t9)*t11*s23r-512*t146*t11)*s1vr+512*t68*t5*t
     #76*t201-1024*t68*t4*l3*t91+512*t68*t9*t60

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

      t1 = z**2
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
      taux3 = ((-128*t1-256*z*s1l1-256*t5)*t10*t13*t61/t63*t66+256*s1vr*
     #(2*s1l1+z)*t121*t66*t13*t9*t60/t62-256*t130*t131*t66*t13)/t137
      t1 = z**2
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
      taux4 = (-512*(-t1-2*z*s1l1)*t7*t10*t58/t60*t63-1024*s1vr*t6*t10*t
     #57/t59*t63*z*(-t30*s1l1-2*t12*t105*c3l1+t109*t45+t34*s1l1*t18-2*t1
     #05*c3l1*r-t34*t116-t34*t116*t18+t109*t49))/t127

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
      t50 = dsqrt(-l1*(l1-1)*(1-l2)*t38*r*(t6*t+t6*l1*t+r+t15*r-t44-t44*
     #t15+t+l1*t))
      t54 = t9*t*t13-t13*t16-t13*t6+t2*t13+t20*t21+2*t20*t16+4*t20*t6+t2
     #0-t12*t9*t-t29*t15-t29-t12*r*t+2*t50-4*t50*l4
      t55 = t54**2
      t57 = (r+t)**2
      t61 = s4vr**2
      taux3 = 128*t1*t4/t6*t55/t57/t38/t61
      t3 = s13*l3
      t4 = s1vr*s2l1
      t6 = l3*s1vr
      t12 = z**2
      t15 = r*t
      t16 = -1+t15
      t17 = t16**2
      t19 = r**2
      t20 = 1/t19
      t22 = t19*r
      t23 = t22*t
      t24 = l1**2
      t25 = dsqrt(l2)
      t26 = t24*t25
      t28 = t**2
      t29 = t19*t28
      t33 = l1*t25
      t34 = t19**2
      t40 = t25*t22
      t42 = t25*t19
      t44 = t25*r
      t51 = (1+t19)**2
      t54 = t19*l1
      t55 = t54*t
      t57 = l1*r
      t59 = l1*t
      t63 = dsqrt(-l1*(l1-1)*(1-l2)*t51*r*(t19*t+t55+r+r*t28-t57-t57*t28
     #+t+t59))
      t67 = t23*t26-t26*t29-t26*t19+t15*t26+t33*t34+2*t33*t29+4*t33*t19+
     #t33-t40*t-t42*t28-t42-t44*t+2*t63-4*t63*l4
      t68 = t67**2
      t69 = r+t
      t70 = t69**2
      t73 = 1/t51
      t77 = s1vr*t16*t20
      t79 = t67/t69
      t101 = l1*l2
      t105 = l2*t19
      t113 = l2*t24
      t118 = t63*t25
      t122 = -l1-2*t54-l1*t34+t101+4*t101*t19+t101*t34-t105-l2*t22*t-l2*
     #r*t-t105*t28+2*t101*t29-t113*t19+t113*t23+t113*t15-t113*t29-4*t118
     #*l4+2*t118
      t128 = dsqrt((r-t)*s1l1*s2l1*t122*t16/r*t73)
      t132 = t25*s1l1
      t139 = s1l1*l1
      t148 = s4vr**2
      taux4 = (-512*(6*z*s1l1+t3*t4+t6*s2l1-t6*s2l1*z-t3*t4*z+3*t12)*t17
     #*t20*t68/t70*t73-512*(-t77*t79*t73*(-t3*s1vr+s34*s13-s13-t6+t3*s1v
     #r*z+2*z+s34-1+t6*z)-4*t77*t79*t73*z)*(-t40*s1l1-2*t22*t128*c3l1+t1
     #32*t55+t44*s1l1*t28-2*t128*c3l1*r-t44*t139-t44*t139*t28+t132*t59))
     #/t148

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
      t3 = s13*l3
      t4 = s1vr*s2l1
      t8 = z**2
      t9 = l3*s1vr
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
      t66 = t25*t22-t25*t28-t25*t18+t14*t25+t32*t33+2*t32*t28+4*t32*t18+
     #t32-t39*t-t41*t27-t41-t43*t+2*t62-4*t62*l4
      t67 = t66**2
      t68 = r+t
      t69 = t68**2
      t72 = 1/t50
      t93 = l1*l2
      t97 = l2*t18
      t105 = l2*t23
      t110 = t62*t24
      t114 = -l1-2*t53-l1*t33+t93+4*t93*t18+t93*t33-t97-l2*t21*t-l2*r*t-
     #t97*t27+2*t93*t28-t105*t18+t105*t22+t105*t14-t105*t28-4*t110*l4+2*
     #t110
      t120 = dsqrt((r-t)*s1l1*s2l1*t114*t15/r*t72)
      t124 = t24*s1l1
      t131 = s1l1*l1
      t142 = s4vr**2
      taux4 = (-1024*(-2*z*s1l1-t3*t4+t3*t4*z-t8-t9*s2l1+t9*s2l1*z)*t16*
     #t19*t67/t69*t72-1024*s1vr*t15*t19*t66/t68*t72*(-t3*s1vr+s34*s13-s1
     #3-t9+t3*s1vr*z+2*z+s34-1+t9*z)*(-t39*s1l1-2*t21*t120*c3l1+t124*t54
     #+t43*s1l1*t27-2*t120*c3l1*r-t43*t131-t43*t131*t27+t124*t58))/t142

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

