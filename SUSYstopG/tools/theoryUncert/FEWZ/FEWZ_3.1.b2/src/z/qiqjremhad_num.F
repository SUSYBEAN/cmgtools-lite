CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqjrmnum1e0(var,ndim,z,ux,l1,
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

      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = z*t3
      t5 = l3-1
      t6 = t5**2
      t7 = t6*t5
      t8 = t4*t7
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
      t37 = z*t10
      t43 = z*t6
      t51 = t6*t3
      t52 = t51*l3
      t53 = t11*s23r
      t54 = t37*t53
      t57 = t2*t10
      t64 = 2*z
      t65 = t64+1
      t66 = t65*z
      t74 = z**2
      t75 = t6*t74
      t85 = t5*t2
      t86 = t85*l3
      t87 = s1vr**2
      t93 = l3*t3
      t94 = t6*t10
      t95 = t94*s34r
      t96 = t93*t95
      t98 = 1+z
      t100 = t1*t10
      t103 = 768*t96-256*t98*t5*t100
      t107 = s23r**2
      t108 = t12*t107
      t112 = t5*t1
      t113 = t112*t10
      t119 = t2**2
      t121 = s34r**2
      t122 = t119*t10*t121
      t130 = l3*t10
      t131 = t130*s2l1
      t132 = t85*t131
      t134 = t130*s3l1
      t135 = t85*t134
      t136 = 256*t135
      t137 = l3*z
      t138 = t137*t10
      t139 = t85*t138
      t157 = t10*s3l1
      t160 = 2+z
      t168 = z*t5
      t169 = t168*t2
      t170 = t11*t107
      t171 = t130*t170
      t174 = t37*s34r
      t175 = t52*t174
      t190 = t5*t3
      t191 = l3**2
      t192 = t190*t191
      t193 = t87*s1vr
      t200 = s2l1**2
      t201 = t130*t200
      t208 = t6*t119*t191
      t210 = 128*t208*t174
      t211 = t74*t5
      t214 = t211*t2*l3*t10
      t215 = 256*t214
      t225 = t37*s3l1
      t226 = t86*t225
      t242 = z*t119*t1
      t249 = t74*l3*t3
      t250 = t249*t95
      t252 = s3l1**2
      t259 = 2*t74
      t260 = 1+t259+t64
      t268 = t211*t2
      t288 = t2*t191
      t289 = t288*t10
      t293 = t191*t3
      t294 = t5*t10
      t295 = t294*s34r
      t296 = t293*t295
      t297 = l3*t1
      t298 = t297*t37
      t300 = -256*t296+256*t298
      t303 = t130*s34r
      t304 = t85*t303
      t306 = 256*t304-256*t37
      t311 = t12*t107*s23r
      t314 = t297*t10
      t319 = t191*t119
      t320 = t94*t121
      t337 = t288*t37
      t379 = -t1
      t393 = t192*t174
      t398 = t37*t121
      t399 = t208*t398
      t415 = t191*l3
      t416 = t4*t415
      t434 = 128*t5*t119*t415*t174
      t435 = t10*t252
      t455 = t74*t1*t130
      t470 = z*t1
      t474 = -t98*z
      t490 = t379*z
      t538 = t74*t10
      s3 = -64*t8*t12*s1vr*t13*s13r
      s5 = ((-256*t19*t12*s1vr+128*t19*t23+256*t3*t7*t28)*t32+((-256*t19
     #*t34-384*t19*t37)*t11*s1vr+128*t43*t2*t23+256*t8*t28)*s1l1+(64*t52
     #*t54+(-128*t43*t57*s2l1-128*t43*t57*s3l1-64*t66*t19*t10)*t11)*s1vr
     #+64*t75*t2*t23+128*t74*t3*t7*t28)*t13
      s8 = (-256*t86*t12*t87+(256*t86*t23+t103*t11)*s1vr-256*t86*t108+(-
     #512*t96+256*t113)*t11*s23r+(-512*l3*t7*t122+512*t19*t27)*t11)*t32+
     #((-512*t132+t136-384*t139)*t11*t87+((256*t132-t136+512*t139)*t11*s
     #23r+(t103*s2l1+(-256*t93*t94*s3l1+768*t51*t138)*s34r-256*t85*t157-
     #256*t160*z*t113)*t11)*s1vr-256*t169*t171+(-512*t175+256*t168*t100)
     #*t11*s23r+(-512*t137*t7*t122+512*t19*t174)*t11)*s1l1-64*t192*t37*t
     #11*t193+(128*t192*t54+(-256*t85*t201+(256*t135-256*t139)*s2l1+t210
     #-t215)*t11)*t87
      s7 = s8+(-64*t192*t37*t170+(128*t169*t131-t210+128*t226+t215)*t11*
     #s23r+(-128*t168*t100*t200+(256*t175-128*t160*z*t113)*s2l1-128*t242
     #*t7*t191*t10*t121+384*t250-128*t168*t100*t252-128*t211*t100*s3l1-6
     #4*t260*z*t113)*t11)*s1vr-128*t268*t171+(-256*t250+128*t211*t100)*t
     #11*s23r+(-256*t74*t119*t7*t130*t121+256*t75*t57*s34r)*t11
      s8 = s13r
      s6 = s7*s8
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(128*t289*t53*t87+(t300*t11*s23r+t306*t11)*s1vr+128*t288*t
     #311+(256*t296-256*t314)*t11*t107+(256*t319*t320-512*t304+128*t10)*
     #t11*s23r+(-512*t93*t320+256*t112*t27)*t11)*t32+(((256*t288*t34+256
     #*t337)*t11*s23r+(256*t297*t157-128*t298)*t11)*t87+((256*t288*t157-
     #128*t337)*t11*t107+(t300*s2l1+(256*t293*t294*s3l1-256*t190*t191*z*
     #t10)*s34r+256*(z-2)*t1*t134+256*t98*z*t314)*t11*s23r+(t306*s2l1+(-
     #768*t135+256*t139)*s34r+256*t379*t10*s3l1-128*t66*t10)*t11)*s1vr+1
     #28*z*t2*t191*t311+(256*t393-256*t298)*t11*t107+(256*t399-512*t86*t
     #174+128*t37)*t11*s23r+(-512*t52*t398+256*t112*t174)*t11)*s1l1+(64*
     #t416*t23-64*t288*t37*t11)*t193
      taux1 = s1+(-128*t416*t108+(128*t288*t10*t200+128*t288*t37*s2l1-t4
     #34+128*t288*t435+128*t288*t225+64*(2+3*z)*z*t289)*t11*s23r+(256*t2
     #97*t157*s2l1+128*t393-256*t297*t435-128*t455)*t11)*t87+(64*t416*t3
     #11+(t434-64*t65*z*t289)*t11*t107+(128*t470*t201-128*t474*t1*t131+1
     #28*t242*t6*t415*t10*t121+128*t474*t191*t190*t27+128*t470*t130*t252
     #-128*t490*t1*t134+64*t260*z*t314)*t11*s23r+(-128*t37*t200+128*t474
     #*t34-128*t399+(-256*t226+128*t214)*s34r-128*t37*t252+128*t490*t157
     #-64*(1+t259)*z*t10)*t11)*s1vr+64*t74*t2*t191*t311+(128*t74*t191*t3
     #*t295-128*t455)*t11*t107+(128*t319*t74*t320-256*t268*t303+64*t538)
     #*t11*s23r+(-256*t249*t320+128*t112*t538*s34r)*t11
      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = z*t3
      t5 = l3-1
      t6 = t5**2
      t7 = t6*t5
      t8 = t4*t7
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t19 = t2*t6
      t20 = t19*z
      t24 = t12*s23r
      t27 = t10*s34r
      t28 = t27*t11
      t34 = t6*t3*l3
      t35 = z*t10
      t36 = t11*s23r
      t37 = t35*t36
      t40 = z*t6
      t41 = t2*t10
      t48 = 2*z
      t58 = z**2
      t59 = t6*t58
      t70 = t5*t2*l3
      t71 = s1vr**2
      t78 = t35*s34r
      t79 = t34*t78
      t81 = t58*t5
      t82 = t1*t10
      t83 = t81*t82
      t90 = z*t5*t2
      t91 = l3*t10
      t92 = s23r**2
      t93 = t11*t92
      t94 = t91*t93
      t99 = -z*t5
      t107 = t2**2
      t109 = s34r**2
      t119 = l3**2
      t120 = t5*t3*t119
      t121 = t71*s1vr
      t129 = t119*t10
      t133 = t91*s2l1
      t134 = t90*t133
      t137 = t6*t107*t119
      t139 = 512*t137*t78
      t142 = t81*t2*l3*t10
      t143 = 1024*t142
      t153 = t35*s3l1
      t154 = t70*t153
      t160 = 512*t83
      t164 = z*t107*t1
      t170 = t58*l3*t3
      t171 = t6*t10
      t173 = t170*t171*s34r
      t177 = t5*t10
      t181 = 2*t58
      t182 = -t48-t181+1
      t192 = t81*t2
      t210 = t119*t2
      t211 = t210*z
      t214 = l3*t1
      t215 = t35*t11
      t220 = t12*t92
      t223 = t120*t78
      t225 = t58*t1*t91
      t226 = t223-t225
      t230 = t70*t78
      t241 = t12*t92*s23r
      t249 = t35*t109
      t250 = t137*t249
      t259 = t5*t1
      t266 = t119*l3
      t267 = t4*t266
      t281 = 512*t5*t107*t266*t78
      t287 = t210*t10
      t307 = -z+t58
      t308 = t307*t1
      t370 = t171*t109
      t376 = t58*t10
      s2 = -256*t8*t12*s1vr*t13*s13r+((512*t20*t12*s1vr-512*t20*t24-1024
     #*t8*t28)*s1l1+(256*t34*t37+(512*t40*t41*s2l1-512*t40*t41*s3l1-256*
     #(-t48+1)*z*t19*t10)*t11)*s1vr-256*t59*t2*t24-512*t58*t3*t7*t28)*t1
     #3
      s3 = s2+((1536*t70*t35*t11*t71+(-1024*t70*t37+(-3072*t79+1024*t83)
     #*t11)*s1vr+1024*t90*t94+(2048*t79+1024*t99*t82)*t11*s23r+(2048*z*l
     #3*t7*t107*t10*t109-2048*t19*t78)*t11)*s1l1+256*t120*t35*t11*t121+(
     #-512*z*t5*t3*t129*t36+(1024*t134-t139+t143)*t11)*t71+(-256*t99*t3*
     #t129*t93+(-512*t134+t139+512*t154-t143)*t11*s23r+((-1024*t79+t160)
     #*s2l1+512*t164*t7*t129*t109-1536*t173+512*(t58-t48)*t1*t177*s3l1-2
     #56*t182*t5*z*t1*t10)*t11)*s1vr+512*t192*t94+(1024*t173-t160)*t11*s
     #23r+(1024*t58*t107*t7*t91*t109-1024*t59*t41*s34r)*t11)*s13r
      s1 = s3+((-1024*t211*t24+512*t214*t215)*t71+(512*t211*t220+1024*t2
     #26*t11*s23r+(-1024*t230-256*(t48-4*t58)*t10)*t11)*s1vr-512*t211*t2
     #41+(-1024*t223+1024*t214*t35)*t11*t92+(-1024*t250+2048*t230-512*t3
     #5)*t11*s23r+(2048*t34*t249-1024*t259*t78)*t11)*s1l1+(-256*t267*t24
     #+256*t210*t215)*t121
      taux2 = s1+(512*t267*t220+(-512*t210*t35*s2l1+t281-512*t210*t153+2
     #56*(-2-3*z)*z*t287)*t11*s23r-512*t226*t11)*t71+(-256*t267*t241+(-t
     #281+256*(t48+1)*z*t287)*t11*t92+(-512*t308*t133-512*t164*t6*t266*t
     #10*t109+512*(1+z)*t5*z*t3*t119*t27-512*t308*t91*s3l1+256*t182*z*t2
     #14*t10)*t11*s23r+(512*z*t1*t10*s2l1+512*t250+(1024*t154-512*t142)*
     #s34r+512*t307*t10*s3l1-256*(1-t181)*z*t10)*t11)*s1vr-256*t58*t2*t1
     #19*t241+(-512*t58*t119*t3*t177*s34r+512*t225)*t11*t92+(-512*t119*t
     #107*t58*t370+1024*t192*t91*s34r-256*t376)*t11*s23r+(1024*t170*t370
     #-512*t259*t376*s34r)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux2
      tcoff10 = taux2
      tcoff11 = taux2
      tcoff12 = taux2

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet1 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff1)
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
     #l2,1.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet3 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet4 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet5 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
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
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet7 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet8 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,10.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet9 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff9)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet11 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet12 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff12)

CCC Sum

      qiqjrmnum1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet
     #9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*tcoff12
      RETURN
      END

      double precision function qiqjrmnum1e1(var,ndim,z,ux,l1,
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

      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = z*t3
      t5 = l3-1
      t6 = t5**2
      t7 = t6*t5
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t19 = t6*t2
      t21 = t11*s23r
      t22 = s1l1**2
      t26 = z*t6
      t34 = z*t10
      t35 = t34*t21
      t38 = t2*t10
      t42 = t2**2
      t45 = l3*t10
      t46 = t45*s34r
      t58 = z**2
      t59 = t6*t58
      t61 = t12*s23r
      t72 = t5*t2
      t73 = t72*l3
      t74 = s23r*s1vr
      t77 = s23r**2
      t78 = t12*t77
      t80 = t5*t1
      t85 = t45*s2l1
      t89 = l3*z
      t96 = z*t5
      t97 = t96*t2
      t98 = t11*t77
      t99 = t45*t98
      t108 = l3**2
      t109 = t5*t3*t108
      t110 = s1vr**2
      t111 = t110*s1vr
      t118 = t97*t85
      t120 = t6*t42*t108
      t121 = t34*s34r
      t122 = t120*t121
      t123 = t34*s3l1
      t124 = t73*t123
      t125 = 2*z
      t130 = t2*l3*t10
      t151 = t1*t10
      t152 = s2l1**2
      t156 = t80*t123
      t162 = z*t42*t1
      t164 = t108*t10
      t165 = s34r**2
      t170 = t58*l3*t3
      t171 = t6*t10
      t173 = t170*t171*s34r
      t175 = s3l1**2
      t191 = t58*t5*t2
      t214 = t2*t108
      t215 = t214*t10
      t220 = t1*l3
      t226 = t12*t77*s23r
      t234 = t10*s2l1
      t235 = t214*t234
      t237 = t10*s3l1
      t238 = t214*t237
      t240 = t214*t34
      t252 = t220*t34
      t259 = z*t2
      t269 = t108*l3
      t270 = t4*t269
      t288 = 128*t5*t42*t269*t121
      t292 = t58*t2
      t302 = t89*t151*s3l1
      t319 = z*t1
      t337 = 128*t58*t108*t3*t5*t10*s34r
      t386 = t171*t165
      t391 = t58*t10
      s3 = 128*t4*t7*t12*s1vr*t13*s13r
      s5 = (-128*t19*t10*t21*t22-128*t26*t2*t12*s23r*s1l1+(-256*t6*t3*l3
     #*t35+(128*t26*t38*s2l1-128*t7*t42*z*t46+128*t26*t38*s3l1+256*t19*t
     #34)*t11)*s1vr-128*t59*t2*t61-128*t58*t3*t7*t10*s34r*t11)*t13
      s6 = ((256*t73*t12*t74+256*t73*t78-256*t80*t61)*t22+((256*t72*t85+
     #256*t72*t45*s3l1+256*t72*t89*t10)*t11*t74+256*t97*t99-256*t96*t1*t
     #61)*s1l1+128*t109*t34*t11*t111+(-128*t109*t35+(-128*t118-128*t122-
     #128*t124+128*(t125-1)*z*t5*t130)*t11)*t110+(256*t109*t34*t98+(-128
     #*t118+256*t122-128*t124-128*(2+z)*z*t5*t130)*t11*s23r+(128*t96*t15
     #1*t152+(256*t156+256*t96*t151)*s2l1+128*t162*t7*t164*t165-384*t173
     #+128*t96*t151*t175+256*t156-64*(t125-2*t58-4)*z*t80*t10)*t11)*s1vr
     #+256*t191*t99+(256*t173-256*t58*t5*t151)*t11*s23r+(256*t58*t42*t7*
     #t45*t165-256*t59*t38*s34r)*t11)*s13r
      s4 = s5+s6
      s2 = s3+s4
      s1 = s2+(-128*t215*t21*t110+(-256*t214*t78+256*t220*t61)*s1vr-128*
     #t214*t226+256*t220*t78-128*t61)*t22+((-256*t235-256*t238-128*t240)
     #*t11*s23r*t110+((-256*t235-256*t238-256*t240)*t11*t77+(256*t220*t2
     #34+256*t220*t237+256*t252)*t11*s23r)*s1vr-128*t259*t108*t226+256*t
     #220*z*t78-128*t35)*s1l1+(-128*t270*t61+128*t214*t34*t11)*t111
      taux1 = s1+(128*t270*t78+(-128*t214*t10*t152-256*t214*t237*s2l1+t2
     #88-128*t214*t10*t175-256*t292*t164)*t11*s23r+(-128*t220*t34*s2l1-1
     #28*t109*t121-128*t302+128*t259*t45)*t11)*t110+(-128*t270*t226+(-t2
     #88+64*(2*z+2)*z*t215)*t11*t77+(-128*t319*t45*t152+(-256*t302-128*t
     #252)*s2l1-128*t162*t6*t269*t10*t165+t337-128*t319*t45*t175-128*t30
     #2-64*(2+2*t58)*z*t220*t10)*t11*s23r+(128*t34*t152+(256*t123+128*t3
     #4)*s2l1+128*t120*t34*t165-128*t4*t5*t46+128*t34*t175+128*t123-64*(
     #-2*t58-2+2*z)*z*t10)*t11)*s1vr-128*t292*t108*t226+(-t337+256*t58*t
     #1*t45)*t11*t77+(-128*t108*t42*t58*t386+256*t191*t46-128*t391)*t11*
     #s23r+(256*t170*t386-128*t80*t391*s34r)*t11
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t16 = z*t6
      t17 = t16*t3
      t18 = t9*t11
      t23 = s13r**2
      t26 = t6*t2
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t35 = t5**2
      t36 = t3*t35
      t37 = l3*t9
      t38 = t37*s34r
      t40 = 1024*t36*t38
      t41 = 6*z
      t42 = t41-1
      t44 = t5*t9
      t51 = t5*t2
      t56 = t9*s34r*t11
      t61 = t37*s2l1
      t65 = t2*t9
      t69 = l3*z
      t70 = t69*t9
      t76 = z*t9
      t81 = 3*z
      t95 = -7+t41
      t104 = z**2
      t105 = t2*t104
      t109 = t104*t6
      t115 = t1*t6
      t116 = l3**2
      t117 = t115*t116
      t118 = t28*s1vr
      t122 = -t1
      t123 = t122*t6
      t127 = t2*t35
      t128 = t116*t9
      t129 = t128*s34r
      t130 = t127*t129
      t132 = t5*t1
      t139 = s23r**2
      t140 = t18*t139
      t144 = 8*z
      t148 = t1*l3*t9
      t154 = t35*t4
      t155 = t3*t154
      t156 = s34r**2
      t158 = t155*t128*t156
      t160 = t76*s34r
      t161 = t27*t160
      t163 = 12*t104
      t166 = t1*t9
      t175 = t11*t139
      t176 = t37*t175
      t182 = t4*t9
      t198 = t128*s2l1
      t201 = 3*t122
      t203 = t128*s3l1
      t207 = t116*z*t9
      t208 = t115*t207
      t213 = t116*t6
      t214 = t166*s3l1
      t219 = t213*t9
      t231 = t127*t203
      t232 = t127*t207
      t236 = t132*l3
      t237 = t76*s3l1
      t238 = t236*t237
      t240 = t104*t5
      t260 = t127*t116
      t264 = 7*z
      t278 = 3*t104
      t303 = 6*t104
      t315 = t104*t1*t5
      t319 = t104*l3*t6
      t323 = t104*t4
      t340 = t116*l3
      t341 = t6*t340
      t343 = t5*t116
      t350 = t122*t35
      t351 = t340*t9
      t353 = t350*t351*s34r
      t354 = 1024*t353
      t355 = t343*t76
      t360 = t123*t129
      t362 = 1-t81
      t371 = t18*t139*s23r
      t374 = 10*z
      t383 = t351*t156
      t386 = t201*z
      t389 = 1024*t386*t6*t129
      t398 = t116*t35
      t399 = t65*t156
      t400 = t398*t399
      t417 = l3*t4
      t423 = t76*t156
      t441 = t9*s2l1
      t443 = 256*t341*t441
      t444 = t9*s3l1
      t446 = 256*t341*t444
      t447 = t16*t351
      t479 = 5*z
      t486 = t343*t9
      t493 = 1+z
      t506 = t37*s3l1
      t511 = t417*t9
      t521 = -z+2
      t538 = t521*t1
      t625 = t104*t9
      s2 = (-256*t3*t6*t9*t11*s1vr*s1l1+768*t17*t18*s1vr)*t23*s13r+((-51
     #2*t27*t18*t28+(768*t27*t32+(t40-256*t42*t2*t44)*t11)*s1vr+2048*t51
     #*z*t32+3072*t17*t56)*s1l1+(768*t26*t61+256*l3*t6*t65*s3l1-256*t26*
     #t70)*t11*t28+(-512*t27*t76*t11*s23r+((-t40-256*(-2+t81)*t2*t44)*s2
     #l1+512*t36*z*t38+1792*z*t2*t44*s3l1-256*t95*z*t51*t9)*t11)*s1vr+10
     #24*t105*t5*t32+1536*t109*t3*t56)*t23
      s4 = s2
      s7 = (-256*t117*t18*t118+(-1024*t123*t116*t32+(1024*t130-4608*t132
     #*t70)*t11)*t28+(-768*t117*t140+(-2048*t130-512*(1-t144)*t5*t148)*t
     #11*s23r+(-1024*t158+9216*t161-256*(-1+t163)*t4*t166)*t11)*s1vr-409
     #6*z*t1*t5*t176+(-6144*t161+4096*z*t1*t182)*t11*s23r+(-6144*t69*t3*
     #t35*t9*t156+6144*t51*t160)*t11)*s1l1+(-256*t115*t198-256*t201*t6*t
     #203-1024*t208)*t11*t118+((-512*t213*t214-2304*t122*z*t219)*t11*s23
     #r+(256*(2-11*z)*t5*t1*t61+(-1536*t231+1536*t232)*s34r+768*t238-307
     #2*t240*t148)*t11)*t28
      s6 = s7+(-1280*z*t1*t6*t128*t175+((512*t130+1024*t5*t4*t148)*s2l1-
     #2048*t260*t160-2048*t238-512*(2-t264)*z*t5*t148)*t11*s23r+((512*t1
     #58-512*(-t264+4)*t6*t2*t38+512*(2*z-t278+1)*t4*t166)*s2l1+(512*t15
     #5*t203-1536*t155*t207)*t156+(-512*t27*t237+4608*t109*t2*l3*t9)*s34
     #r+512*(t144-t278-1)*t4*t214-256*(t303-7+t41)*t1*z*t4*t9)*t11)*s1vr
     #-2048*t315*t176+(-3072*t319*t65*s34r+2048*t323*t166)*t11*s23r+(-30
     #72*t104*t35*t3*t37*t156+3072*t105*t44*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((256*t341*t32-256*t343*t18)*t118+(-512*t341*t140+(t354+35
     #84*t355)*t11*s23r+(-1024*t360+512*t362*t4*t37)*t11)*t28+(256*t341*
     #t371+(-t354+256*(-t374+1)*t5*t128)*t11*t139+(1024*t2*t154*t383+t38
     #9-256*(-t163+1-4*z)*t4*t37)*t11*s23r+(-1024*t400-1024*t362*t5*t1*t
     #38-256*(1+t163-t41)*t9)*t11)*s1vr+2048*z*t5*t116*t371+(-t389-4096*
     #t417*t76)*t11*t139+(3072*t260*t423-6144*t236*t160+2048*t76)*t11*s2
     #3r+(-6144*t27*t423-1024*t386*t182*s34r)*t11)*s1l1+((t443-t446+1024
     #*t447)*t11*s23r+(-256*t343*t441+256*t343*t444-1024*t355)*t11)*t118
      s3 = s1
      s5 = ((-t443+t446-2048*t447)*t11*t139+((512*t353+1792*t355)*s2l1+(
     #-512*t350*t351*s3l1-1536*t340*t1*t35*z*t9)*s34r+256*(t479-2)*t5*t2
     #03+256*(t374+7)*z*t486)*t11*s23r+((-512*t360+256*t493*t4*t37)*s2l1
     #+(1024*t123*t203+1536*t208)*s34r+256*(t81+1)*t4*t506-256*t42*z*t51
     #1)*t11)*t28
      s6 = (1024*t16*t340*t371+(256*t521*t5*t198+1536*t1*t35*t340*t160+2
     #56*t343*t237+256*(-3-t144)*z*t486)*t11*t139+((512*t538*t6*t129-512
     #*(-t278+t81+1)*t4*t37)*s2l1+1536*z*t154*t2*t383+(512*t493*t1*t6*t2
     #03+512*(-t81-2)*t1*z*t219)*s34r-512*(t479-t278-1)*t4*t506-256*(-t1
     #44-t303+4)*z*t511)*t11*s23r+((512*t400-512*t538*t5*t38-256*t95*z*t
     #9)*s2l1+(512*t231-1536*t232)*t156+(-3584*t238-512*t362*z*t5*t148)*
     #s34r+256*(-t303-2+9*z)*t9*s3l1-256*(-3+t303)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2+1024*t240*t116*t371+(1536*t104*t116*t6*t166*s34r-2048*t
     #323*t37)*t11*t139+(1536*t398*t104*t399-3072*t315*t38+1024*t625)*t1
     #1*s23r+(-3072*t319*t399+1536*t1*t4*t625*s34r)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux2
      tcoff10 = taux2
      tcoff11 = taux2
      tcoff12 = taux2

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet1 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff1)
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
     #l2,1.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet3 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet4 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet5 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
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
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet7 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet8 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,10.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet9 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff9)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet11 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet12 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff12)

CCC Sum

      qiqjrmnum1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet
     #9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*tcoff12
      RETURN
      END

      double precision function qiqjrmnum1e2(var,ndim,z,ux,l1,
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

      t1 = -1+z
      t2 = t1**2
      t3 = t2*t1
      t4 = z*t3
      t5 = l3-1
      t6 = t5**2
      t9 = s4vr**2
      t10 = 1/t9
      t11 = s23r2(z,l1,l2,l3,l4)
      t12 = t10*t11
      t13 = s13r**2
      t22 = s1vr**2
      t29 = z*t10
      t30 = t11*s23r
      t34 = t29*t11
      t39 = z**2
      t42 = t12*s23r
      t47 = l3**2
      t49 = z*t47*t3
      t50 = t5*t10
      t51 = t22*s1vr
      t58 = t5*t2
      t64 = s23r**2
      t65 = t11*t64
      t72 = l3*t10
      t73 = t72*t30
      t88 = t39*t1
      t95 = t4*t47*l3
      t98 = t47*t10
      t104 = t12*t64
      t119 = t12*t64*s23r
      t123 = (2*z-3)*z
      taux1 = -64*t4*t6*t5*t12*s1vr*t13*s13r+(128*z*l3*t3*t6*t10*t11*t22
     #+(192*t6*t3*l3*t29*t30-192*t6*t2*t34)*s1vr+64*t6*t39*t2*t42)*t13+(
     #-64*t49*t50*t11*t51+(-256*t49*t50*t30+256*t58*l3*t34)*t22+(-192*t4
     #9*t50*t65-128*(-3+z)*z*t58*t73-192*z*t5*t1*t10*t11)*s1vr-128*t39*t
     #5*t2*t72*t65+128*t88*t5*t42)*s13r+(64*t95*t42-64*z*t2*t98*t11)*t51
     #+(128*t95*t104+64*(-4+z)*z*t2*t98*t30+128*l3*t1*t34)*t22+(64*t95*t
     #119+64*t123*t2*t98*t65-64*t123*t1*t73-64*t34)*s1vr+64*t39*t2*t47*t
     #119-128*t88*l3*t104+64*t39*t10*t30
      t1 = l3-1
      t2 = t1**2
      t3 = t2*t1
      t4 = -1+z
      t5 = t4**2
      t6 = t5*t4
      t8 = s4vr**2
      t9 = 1/t8
      t11 = s23r2(z,l1,l2,l3,l4)
      t16 = z*t6
      t17 = t16*t3
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = t18*s23r
      t31 = t5**2
      t33 = l3*t9
      t34 = t33*s34r
      t36 = 2048*t3*t31*t34
      t37 = -z-1
      t40 = t5*t9
      t47 = t2*t5
      t52 = t9*s34r*t11
      t57 = t33*s2l1
      t61 = t2*t9
      t62 = t61*s3l1
      t67 = s1vr**2
      t69 = z*t9
      t80 = l3*z
      t90 = 4*z
      t100 = z**2
      t101 = t2*t100
      t105 = t100*t6
      t111 = -t1
      t112 = 4*t111
      t113 = t112*t6
      t114 = l3**2
      t115 = t113*t114
      t116 = t67*s1vr
      t120 = t2*t31
      t121 = t114*t9
      t122 = t121*s34r
      t123 = t120*t122
      t125 = 6*z
      t129 = t1*l3*t9
      t135 = s23r**2
      t140 = 4+t125
      t147 = t31*t4
      t148 = t3*t147
      t149 = s34r**2
      t150 = t121*t149
      t151 = t148*t150
      t153 = t69*s34r
      t154 = t27*t153
      t156 = 8*t100
      t159 = t1*t9
      t167 = t1*t5*l3
      t168 = t11*t135
      t191 = t121*s2l1
      t196 = t121*s3l1
      t199 = t1*t6
      t201 = t114*z*t9
      t202 = t199*t201
      t213 = t6*t114
      t214 = t213*t9
      t225 = t120*t196
      t227 = t120*t201
      t231 = -t140
      t234 = t33*s3l1
      t237 = 2*z
      t260 = z*t31
      t264 = t69*s3l1
      t265 = t167*t264
      t267 = 5*z
      t276 = 3*z
      t303 = t159*s3l1
      t306 = 4*t100
      t317 = t100*t1
      t318 = t317*t5
      t322 = t105*t2
      t325 = t100*t4
      t331 = t100*t31
      t333 = t33*t149
      t342 = t114*l3
      t343 = t6*t342
      t345 = t5*t114
      t350 = t1*t31
      t351 = t342*t9
      t353 = t350*t351*s34r
      t354 = 2048*t353
      t355 = 7*z
      t364 = t213*t159*s34r
      t375 = t18*t135*s23r
      t378 = 8*z
      t392 = 2048*t199*t114*t153
      t402 = t114*t31*t61*t149
      t417 = l3*t4
      t424 = t69*t149
      t435 = t1*t4
      t442 = t9*s2l1
      t445 = t16*t351
      t451 = 1024*t345*t442
      t474 = t342*t1
      t487 = t345*t9
      t510 = t417*t9
      t533 = (z-2)*t1
      t537 = 2*t100
      t616 = t100*t9
      s2 = (1024*t3*t6*t9*t11*s1vr*s1l1-256*t17*t18*s1vr)*t23*s13r+(((-3
     #072*t27*t28+(-t36-1024*t37*t2*t40)*t11)*s1vr-2560*t47*z*t28-2048*t
     #17*t52)*s1l1+(-2048*t26*t57-1024*l3*t6*t62)*t11*t67+(-768*t27*t69*
     #t11*s23r+((t36-512*(2+z)*t2*t40)*s2l1-1024*t80*t31*t3*t9*s34r-1536
     #*z*t2*t40*s3l1-256*(7-t90)*z*t47*t9)*t11)*s1vr-1280*t101*t5*t28-10
     #24*t105*t3*t52)*t23
      s4 = s2
      s7 = (-256*t115*t18*t116+(-2048*t123+512*(-4+t125)*t5*t129)*t11*t6
     #7+(-768*t115*t18*t135+(4096*t123-512*t140*t5*t129)*t11*s23r+(2048*
     #t151-6144*t154-256*(-4-t156)*t4*t159)*t11)*s1vr+5120*t167*t69*t168
     #+(4096*t154+5120*t111*z*t4*t9)*t11*s23r+(4096*t80*t3*t31*t9*t149-4
     #096*t47*t153)*t11)*s1l1+(-256*t113*t191+256*t112*t6*t196+1280*t202
     #)*t11*t116+((1024*t199*t191-512*t113*t196-1536*t1*z*t214)*t11*s23r
     #+(256*(-8+t125)*t5*t1*t57+(3072*t225-1024*t227)*s34r+256*t231*t1*t
     #5*t234+1024*(-1+t237)*z*t5*t129)*t11)*t67
      s6 = s7+(-2304*t111*z*t6*t121*t168+((-1024*t123-512*(-z-4)*t5*t129
     #)*s2l1+2048*t260*t2*t122+2560*t265-512*(-1+t267)*z*t5*t129)*t11*s2
     #3r+((-1024*t151-1024*(t276-4)*t6*t2*t34+1024*(-t237+t100-1)*t4*t15
     #9)*s2l1+(-1024*t148*t196+1024*t148*t201)*t149+(1024*t80*t6*t62-307
     #2*t105*t2*l3*t9)*s34r+1024*(1-t90+t100)*t4*t303-256*(7-t306-t90)*t
     #1*z*t4*t9)*t11)*s1vr+2560*t318*t33*t168+(2048*t322*t34-2560*t325*t
     #159)*t11*s23r+(2048*t331*t3*t333-2048*t101*t40*s34r)*t11
      s7 = s13r
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+((-1024*t343*t28+1024*t345*t18)*t116+((t354+512*(4-t355)*t
     #5*t121)*t11*s23r+(-2048*t364+512*(-4+t237)*t4*t33)*t11)*t67+(-1024
     #*t343*t375+(-t354+256*(4+t378)*t5*t121)*t11*t135+(-2048*t2*t147*t3
     #51*t149+t392-256*(4+t90+t156)*t4*t33)*t11*s23r+(2048*t402-2048*t19
     #9*t34-256*(-4-t156+t90)*t9)*t11)*s1vr-2560*z*t5*t114*t375+(-t392+5
     #120*t417*t69)*t11*t135+(-2048*t120*t114*t424+4096*t167*t153-2560*t
     #69)*t11*s23r+(4096*t27*t424-2048*t435*t153)*t11)*s1l1+((-1024*t343
     #*t442-1280*t445)*t11*s23r+(t451+1280*t345*t69)*t11)*t116
      s3 = s1
      s5 = ((-1024*t343*t9*s3l1+1536*t445)*t11*t135+((1024*t353+256*(4-t
     #378)*t5*t121)*s2l1+(-1024*t350*t351*s3l1+1024*t474*t260*t9)*s34r+2
     #56*(8-t90)*t5*t196+256*(-2-9*z)*z*t487)*t11*s23r+((-1024*t364+256*
     #(-4-t237)*t4*t33)*s2l1+(2048*t213*t303-1024*t202)*s34r+256*t231*t4
     #*t234+1024*z*t4*t510)*t11)*t67
      s6 = (-1280*t16*t342*t375+(-t451-1024*t474*t31*t153-1024*t345*t264
     #+256*(t125+5)*z*t487)*t11*t135+((1024*t533*t6*t122-512*(-t276+t537
     #-2)*t4*t33)*s2l1-1024*z*t342*t2*t147*t9*t149+(1024*t37*t1*t6*t196+
     #1024*t317*t214)*s34r-512*(-t355+t537+2)*t4*t234-256*(t125-1+t306)*
     #z*t510)*t11*s23r+((-1024*t402-1024*t533*t5*t34+512*(t237-3)*z*t9)*
     #s2l1+(-1024*t225+1024*t227)*t149+(3072*t265-1024*t16*t129)*s34r+51
     #2*(-t267+t537+2)*t9*s3l1-256*(1-t306)*z*t9)*t11)*s1vr
      s4 = s5+s6
      s2 = s3+s4
      taux2 = s2-1280*t100*t5*t114*t375+(-1024*t105*t1*t122+2560*t325*t3
     #3)*t11*t135+(-1024*t331*t2*t150+2048*t318*t34-1280*t616)*t11*s23r+
     #(2048*t322*t333-1024*t435*t616*s34r)*t11

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux2
      tcoff10 = taux2
      tcoff11 = taux2
      tcoff12 = taux2

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet1 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff1)
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
     #l2,1.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet3 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet4 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet5 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
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
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet7 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet8 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,10.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet9 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff9)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet11 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet12 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff12)

CCC Sum

      qiqjrmnum1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet
     #9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*tcoff12
      RETURN
      END

      double precision function qiqjrmnum1e3(var,ndim,z,ux,l1,
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
      t16 = z*t6
      t18 = t9*t11
      t23 = s13r**2
      t26 = t2*t6
      t27 = t26*l3
      t28 = s1vr**2
      t32 = t18*s23r
      t34 = t2*t5
      t45 = l3*t9
      t46 = t45*s2l1
      t53 = l3*z*t9
      t60 = t11*s23r
      t61 = t45*t60
      t63 = z*t9
      t64 = t63*t11
      t69 = z**2
      t76 = t1*t6
      t77 = l3**2
      t78 = t76*t77
      t79 = t28*s1vr
      t84 = t5*t1
      t90 = s23r**2
      t91 = t18*t90
      t94 = 2*z
      t105 = z*t5
      t107 = t11*t90
      t108 = t45*t107
      t111 = t1*z
      t117 = t77*t9
      t118 = t117*s2l1
      t119 = t76*t118
      t124 = t77*t6*t1*t9*s3l1
      t128 = t76*t77*z*t9
      t138 = t45*s3l1
      t150 = t84*l3
      t151 = t63*s2l1
      t153 = t63*s3l1
      t171 = t69*t5
      t175 = t69*t4
      t181 = t77*l3
      t182 = t6*t181
      t184 = t5*t77
      t191 = z-4
      t196 = t4*l3
      t202 = t18*t90*s23r
      t206 = 4*z-6
      t228 = t9*s2l1
      t229 = t182*t228
      t231 = t9*s3l1
      t232 = t182*t231
      t235 = t16*t181*t9
      t255 = (t94-4)*t5
      t260 = t184*t9
      t282 = (t94-3)*z
      s2 = (-1024*t3*t6*t9*t11*s1vr*s1l1-512*t16*t3*t18*s1vr)*t23*s13r+(
     #(2048*t27*t18*t28+(3072*t27*t32-3072*t34*t18)*s1vr+1024*z*t2*t5*t3
     #2)*s1l1+(1024*t26*t46+1024*l3*t6*t2*t9*s3l1+1024*t26*t53)*t11*t28+
     #(1536*t16*t2*t61-1536*t34*t64)*s1vr+512*t2*t69*t5*t32)*t23
      s3 = s2+((-1024*t78*t18*t79+(-4096*t78*t32+4096*t84*t45*t11)*t28+(
     #-3072*t78*t91-1024*(t94-6)*t5*t1*t61-3072*t4*t1*t18)*s1vr-2048*t10
     #5*t1*t108+2048*t111*t4*t32)*s1l1+(-1024*t119-1024*t124-512*t128)*t
     #11*t79+((-2048*t119-2048*t124-2048*t128)*t11*s23r+(2048*t84*t46+20
     #48*t84*t138+2048*t84*t53)*t11)*t28+(-1536*t111*t6*t117*t107+(-1024
     #*t150*t151-1024*t150*t153-1024*(z-3)*z*t5*t1*l3*t9)*t11*s23r-1536*
     #t111*t4*t9*t11)*s1vr-1024*t171*t1*t108+1024*t175*t1*t32)*s13r
      s1 = s3+((1024*t182*t32-1024*t184*t18)*t79+(2048*t182*t91+1024*t19
     #1*t5*t77*t32+2048*t196*t18)*t28+(1024*t182*t202+512*t206*t5*t77*t9
     #1-512*t206*t4*l3*t32-1024*t18)*s1vr+1024*t105*t77*t202-2048*t196*z
     #*t91+1024*t63*t60)*s1l1+((1024*t229+1024*t232+512*t235)*t11*s23r+(
     #-1024*t184*t228-1024*t184*t231-512*t184*t63)*t11)*t79
      taux2 = s1+((1024*t229+1024*t232+1024*t235)*t11*t90+(512*t255*t118
     #+512*t255*t117*s3l1+512*t191*z*t260)*t11*s23r+(1024*t196*t228+1024
     #*t196*t231+1024*t196*t63)*t11)*t28+(512*t16*t181*t202+(1024*t105*t
     #118+1024*t184*t153+512*t282*t260)*t11*t90+(-1024*t196*t151-1024*z*
     #t4*t138-512*t282*t196*t9)*t11*s23r-512*t64)*s1vr+512*t171*t77*t202
     #-1024*t175*l3*t91+512*t69*t9*t60

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux2
      tcoff10 = taux2
      tcoff11 = taux2
      tcoff12 = taux2

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet1 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff1)
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
     #l2,1.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet3 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet4 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet5 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
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
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet7 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet8 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,10.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet9 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff9)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet11 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet12 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff12)

CCC Sum

      qiqjrmnum1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet
     #9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*tcoff12
      RETURN
      END

      double precision function qiqjrmnum1e4(var,ndim,z,ux,l1,
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
      tcoff3 = taux1
      tcoff4 = taux1
      tcoff5 = taux1
      tcoff6 = taux1
      tcoff7 = taux2
      tcoff8 = taux2
      tcoff9 = taux2
      tcoff10 = taux2
      tcoff11 = taux2
      tcoff12 = taux2

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet1 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff1)
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
     #l2,1.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet3 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet4 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,11.D0,1.D0,1.D0,Tcoeff*tcoff4)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet5 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff5)
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
     #l2,13.D0,1.D0,1.D0,Tcoeff*tcoff6)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet7 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff7)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet8 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+z
     #-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,10.D0,1.D0,1.D0,Tcoeff*tcoff8)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet9 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z-
     #t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2,
     #s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff9)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff10)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet11 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff11)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet12 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff12)

CCC Sum

      qiqjrmnum1e4 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tc
     #off4+tFjet5*tcoff5+tFjet6*tcoff6+tFjet7*tcoff7+tFjet8*tcoff8+tFjet
     #9*tcoff9+tFjet10*tcoff10+tFjet11*tcoff11+tFjet12*tcoff12
      RETURN
      END

      double precision function qiqjrmnum2e0(var,ndim,z,ux,l1,
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
      t5 = z**2
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
      taux3 = ((-256*t1-256*z*s1l1-128*t5)*t10*t13*t61/t63*t66+256*s1vr*
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

      tcoff13 = taux3
      tcoff14 = taux3
      tcoff15 = taux3
      tcoff16 = taux3
      tcoff17 = taux3
      tcoff18 = taux3
      tcoff19 = taux4
      tcoff20 = taux4
      tcoff21 = taux4
      tcoff22 = taux4
      tcoff23 = taux4
      tcoff24 = taux4

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet13 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff13)
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
     #3l2,1.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet15 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet16 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet17 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff17)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet18 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet19 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet20 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,10.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet21 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff21)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet22 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = -l3+1
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet23 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = -l3+1
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet24 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      qiqjrmnum2e0 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24
     #*tcoff24
      RETURN
      END

      double precision function qiqjrmnum2e1(var,ndim,z,ux,l1,
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
      t5 = s13*l3
      t6 = s1vr*s2l1
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
      taux4 = (-512*(6*z*s1l1+t3*s2l1+t5*t6-t5*t6*z-t3*s2l1*z+3*t12)*t17
     #*t20*t68/t70*t73-512*(-t77*t79*t73*(-t3*s13+s34*s13-s13-t3+t3*s13*
     #z-1+s34+2*z+t3*z)-4*t77*t79*t73*z)*(-t40*s1l1-2*t22*t128*c3l1+t132
     #*t55+t44*s1l1*t28-2*t128*c3l1*r-t44*t139-t44*t139*t28+t132*t59))/t
     #148

CCC Coefficients

      tcoff13 = taux3
      tcoff14 = taux3
      tcoff15 = taux3
      tcoff16 = taux3
      tcoff17 = taux3
      tcoff18 = taux3
      tcoff19 = taux4
      tcoff20 = taux4
      tcoff21 = taux4
      tcoff22 = taux4
      tcoff23 = taux4
      tcoff24 = taux4

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet13 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff13)
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
     #3l2,1.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet15 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet16 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet17 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff17)
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
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet19 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet20 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,10.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet21 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff21)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet23 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet24 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      qiqjrmnum2e1 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24
     #*tcoff24
      RETURN
      END

      double precision function qiqjrmnum2e2(var,ndim,z,ux,l1,
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
      t7 = z**2
      t8 = s13*l3
      t9 = s1vr*s2l1
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
      taux4 = (-1024*(-2*z*s1l1-t3*s2l1+t3*s2l1*z-t7-t8*t9+t8*t9*z)*t16*
     #t19*t67/t69*t72-1024*s1vr*t15*t19*t66/t68*t72*(-t3*s13+t3*s13*z+s3
     #4*s13-s13-t3+t3*z-1+s34+2*z)*(-t39*s1l1-2*t21*t120*c3l1+t124*t54+t
     #43*s1l1*t27-2*t120*c3l1*r-t43*t131-t43*t131*t27+t124*t58))/t142

CCC Coefficients

      tcoff13 = taux3
      tcoff14 = taux3
      tcoff15 = taux3
      tcoff16 = taux3
      tcoff17 = taux3
      tcoff18 = taux3
      tcoff19 = taux4
      tcoff20 = taux4
      tcoff21 = taux4
      tcoff22 = taux4
      tcoff23 = taux4
      tcoff24 = taux4

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet13 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff13)
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
     #3l2,1.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet15 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet16 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet17 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff17)
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
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet19 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet20 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,10.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet21 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff21)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet23 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet24 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      qiqjrmnum2e2 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24
     #*tcoff24
      RETURN
      END

      double precision function qiqjrmnum2e3(var,ndim,z,ux,l1,
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

      tcoff13 = taux3
      tcoff14 = taux3
      tcoff15 = taux3
      tcoff16 = taux3
      tcoff17 = taux3
      tcoff18 = taux3
      tcoff19 = taux4
      tcoff20 = taux4
      tcoff21 = taux4
      tcoff22 = taux4
      tcoff23 = taux4
      tcoff24 = taux4

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet13 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff13)
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
     #3l2,1.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet15 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet16 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet17 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff17)
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
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet19 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet20 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,10.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet21 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff21)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet23 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet24 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      qiqjrmnum2e3 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24
     #*tcoff24
      RETURN
      END

      double precision function qiqjrmnum2e4(var,ndim,z,ux,l1,
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

      tcoff13 = taux3
      tcoff14 = taux3
      tcoff15 = taux3
      tcoff16 = taux3
      tcoff17 = taux3
      tcoff18 = taux3
      tcoff19 = taux4
      tcoff20 = taux4
      tcoff21 = taux4
      tcoff22 = taux4
      tcoff23 = taux4
      tcoff24 = taux4

CCC Jet functions

      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet13 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,1.D0,1.D0,1.D0,Tcoeff*tcoff13)
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
     #3l2,1.D0,1.D0,1.D0,Tcoeff*tcoff14)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet15 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,11.D0,1.D0,1.D0,Tcoeff*tcoff15)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet16 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,11.D0,1.D0,1.D0,Tcoeff*tcoff16)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet17 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,13.D0,1.D0,1.D0,Tcoeff*tcoff17)
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
     #3l2,13.D0,1.D0,1.D0,Tcoeff*tcoff18)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet19 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,10.D0,1.D0,1.D0,Tcoeff*tcoff19)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet20 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,10.D0,1.D0,1.D0,Tcoeff*tcoff20)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet21 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,12.D0,1.D0,1.D0,Tcoeff*tcoff21)
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
     #3l2,12.D0,1.D0,1.D0,Tcoeff*tcoff22)
      t2 = 1-z
      t3 = l3*t2
      t4 = t3*s23r
      t5 = 1-l3
      t7 = t5*t2*s13r
      t8 = t3*s1vr
      t10 = t2**2
      t12 = l3*t5*t10*s34r
      tFjet23 = Fjet(var,ndim,z,1/ux,t4,t7,t8-t4-t12,-1-t7-t8+z,t12,-1+z
     #-t8+t12,t8,-t7-t4-t12+z,1+t7+t4,s2l1,s1l1,s3l1,s4l1,s2l2,s1l2,s3l2
     #,s4l2,14.D0,1.D0,1.D0,Tcoeff*tcoff23)
      t1 = 1-l3
      t2 = 1-z
      t4 = t1*t2*s13r
      t5 = l3*t2
      t6 = t5*s1vr
      t8 = t5*s23r
      t10 = t2**2
      t12 = l3*t1*t10*s34r
      tFjet24 = Fjet(var,ndim,z,ux,-1-t4-t6+z,t6-t8-t12,t4,t8,t12,t6,-1+
     #z-t6+t12,1+t4+t8,-t4-t8-t12+z,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s
     #3l2,14.D0,1.D0,1.D0,Tcoeff*tcoff24)

CCC Sum

      qiqjrmnum2e4 = tFjet13*tcoff13+tFjet14*tcoff14+tFjet15*tcoff15+tFj
     #et16*tcoff16+tFjet17*tcoff17+tFjet18*tcoff18+tFjet19*tcoff19+tFjet
     #20*tcoff20+tFjet21*tcoff21+tFjet22*tcoff22+tFjet23*tcoff23+tFjet24
     #*tcoff24
      RETURN
      END

