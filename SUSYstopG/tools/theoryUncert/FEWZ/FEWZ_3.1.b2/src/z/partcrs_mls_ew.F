
      SUBROUTINE CrsecSeqz_mls()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      t1 = CsTh**2
      t2 = 1+t1
      A1e0c1p1t1z = 12*t2
      t3 = s**2
      t4 = t3*Mzsq
      t5 = abs(DCMPLX(Mzsq))
      t6 = t5**2
      t8 = Mzsq**2
      t12 = 1/(-t4+s*t6+s*t8-t6*Mzsq)
      A1e0c2p1t1z = -12*t4*t2*t12
      t18 = 1/(-s+Mzsq)
      A1e0c3p1t1z = -24*s*t2*t18
      A1e0c1p1t2z = A1e0c1p1t1z
      A1e0c2p1t2z = A1e0c2p1t1z
      A1e0c3p1t2z = A1e0c3p1t1z
      A2e0c1p1t1z = -96*t4*CsTh*t12
      A2e0c2p1t1z = -48*s*CsTh*t18
      A2e0c1p1t2z = A2e0c1p1t1z
      A2e0c2p1t2z = A2e0c2p1t1z

      RETURN
      END

      SUBROUTINE CrsecSeqs_mls()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      DOUBLE COMPLEX spf7
      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX DENSPENCE
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t2 = log(DCMPLX(1/s))
      t3 = log(DCMPLX(dels))
      t4 = t2*t3
      t5 = 12*t4
      t6 = CsTh**2
      t9 = PI**2
      t10 = t9*t6
      t11 = t3**2
      t14 = t2**2
      t19 = 12*t11
      t20 = 3*t14
      t22 = 1/t9
      t23 = (t5+12*t4*t6+t10-12*t11*t6-3*t14*t6-24*t3+12*t2+t9-t19-t20)*
     #t22
      A1e0c1p1t1s = -t23/2
      t25 = 1-CsTh
      spf7 = DENSPENCE(DCMPLX(t25/2),spcZip)
      t27 = log(DCMPLX(t25))
      t28 = t27*t3
      t32 = t27*t2
      t36 = 1+CsTh
      t37 = log(DCMPLX(t36))
      t38 = t37*t2
      t42 = t37*t3
      t46 = log(DCMPLX(2.D0))
      t47 = t37*t46
      t51 = t10+t9+12*t28+12*t28*t6-6*t32-6*t32*t6+6*t38+6*t38*t6-12*t42
     #-12*t42*t6+12*t47+12*t47*t6
      t55 = t46**2
      t57 = t27**2
      t59 = t37**2
      t74 = 12*t27-12*t37-12*spf7-6*t55+3*t57-3*t59-6*t55*t6+3*t57*t6-3*
     #t59*t6-12*spf7*t6-6*t27*t37-6*t6*t27*t37
      t75 = t51+t74
      A1e0c2p1t1s = t75*t22
      t76 = s**2
      t77 = t76*Mzsq
      t78 = abs(DCMPLX(Mzsq))
      t79 = t78**2
      t81 = Mzsq**2
      t85 = 1/(-t77+s*t79+s*t81-t79*Mzsq)
      A1e0c3p1t1s = t77*t23*t85/2
      A1e0c4p1t1s = -t77*t75*t22*t85
      t91 = s*t6
      t94 = conjg(Mzsq)
      t95 = t94*t6
      t98 = t6*Mzsq
      t101 = t11*t94
      t102 = 12*t101
      t103 = t3*Mzsq
      t105 = t3*t94
      t107 = t11*Mzsq
      t108 = 12*t107
      t109 = t3*s
      t111 = t11*s
      t112 = 24*t111
      t113 = t2*Mzsq
      t115 = t2*t94
      t117 = t14*Mzsq
      t118 = 3*t117
      t119 = t14*t94
      t120 = 3*t119
      t121 = t2*s
      t123 = t14*s
      t124 = 6*t123
      t125 = -24*t4*t91+12*t4*t95+12*t4*t98-t102-24*t103-24*t105-t108+48
     #*t109+t112+12*t113+12*t115-t118-t120-24*t121+t124
      t126 = t9*t94
      t127 = t9*Mzsq
      t128 = s*t9
      t129 = 2*t128
      t143 = 24*t4*s
      t145 = 12*t4*t94
      t147 = 12*t4*Mzsq
      t152 = t126+t127-t129+24*t111*t6-12*t101*t6-12*t107*t6-3*t117*t6-3
     #*t119*t6+6*t123*t6-t143+t145+t147+t126*t6+t127*t6-2*t128*t6
      t159 = 1/(t76-s*t94-s*Mzsq+Mzsq*t94)
      t160 = t22*t159
      A1e0c5p1t1s = s*(t125+t152)*t160/2
      spf8 = DENSPENCE(DCMPLX(t36/2),spcZip)
      t163 = t27*t46
      t168 = t27*s
      t169 = t46*t6
      t188 = t37*s
      t200 = -2*t163*t98-2*t163*t95+4*t168*t169+4*t28*t95+4*t28*t98-8*t2
     #8*t91-2*t32*t95-2*t32*t98+4*t32*t91+2*t47*t95+2*t47*t98-4*t188*t16
     #9-4*t42*t98-4*t42*t95+8*t188*t3*t6+2*t38*t98
      t206 = spf7*s
      t207 = 4*t206
      t208 = spf7*Mzsq
      t209 = 2*t208
      t210 = spf7*t94
      t211 = 2*t210
      t212 = t57*s
      t213 = 2*t212
      t215 = t27*Mzsq
      t217 = t27*t94
      t219 = t57*t94
      t220 = t57*Mzsq
      t221 = t59*Mzsq
      t222 = t59*t94
      t223 = t37*Mzsq
      t225 = t37*t94
      t227 = s*t59
      t228 = 2*t227
      t230 = 2*t38*t95-4*t188*t2*t6+t207-t209-t211-t213-8*t168+4*t215+4*
     #t217+t219+t220-t221-t222-4*t223-4*t225+t228+8*t188
      t232 = spf8*t94
      t233 = 2*t232
      t234 = spf8*s
      t235 = 4*t234
      t236 = spf8*Mzsq
      t237 = 2*t236
      t239 = 4*t42*Mzsq
      t243 = 2*t163*t94
      t245 = 2*t38*Mzsq
      t251 = 4*t28*t94
      t253 = 4*t28*Mzsq
      t255 = 8*t28*s
      t257 = 2*t163*Mzsq
      t259 = 4*t168*t46
      t261 = 2*t32*t94
      t262 = t233-t235+t237-t239+4*t206*t6-t243+t245-2*t212*t6+t219*t6+t
     #220*t6+t251+t253-t255-t257+t259-t261
      t264 = 2*t32*Mzsq
      t266 = 4*t32*s
      t272 = 4*t42*t94
      t274 = 8*t188*t3
      t276 = 2*t47*t94
      t278 = 2*t47*Mzsq
      t280 = 4*t188*t46
      t282 = 2*t38*t94
      t284 = 4*t188*t2
      t295 = -t264+t266+2*t227*t6-t221*t6-t222*t6-t272+t274+t276+t278-t2
     #80+t282-t284-2*t208*t6-2*t210*t6-4*t234*t6+2*t236*t6+2*t232*t6
      A1e0c6p1t1s = -3*s*(t200+t230+t262+t295)*t160
      A1e0c1p1t2s = A1e0c1p1t1s
      A1e0c2p1t2s = A1e0c2p1t1s
      A1e0c3p1t2s = A1e0c3p1t1s
      A1e0c4p1t2s = A1e0c4p1t1s
      A1e0c5p1t2s = A1e0c5p1t1s
      A1e0c6p1t2s = A1e0c6p1t1s
      t301 = t77*CsTh
      A2e0c1p1t1s = 4*t301*(-12+t5+t9-t19-t20-36*t3+18*t2)*t22*t85
      t318 = -2*spf7-6*t37+2*spf8+6*t27-t59+t57+2*t47-2*t163+2*t38-2*t32
     #-4*t42+4*t28
      A2e0c2p1t1s = -24*t301*t318*t85*t22
      t323 = s*CsTh
      t334 = -t102-t108+t126-t118-t120+t127+t112+t124-t129-12*t94-12*Mzs
     #q
      A2e0c3p1t1s = t323*(-36*t105+24*s-36*t103+72*t109+18*t113+18*t115-
     #36*t121-t143+t145+t147+t334)*t22*t159
      t344 = t207-t209-t211-t213-12*t168+6*t215+6*t217+t219+t220-t221-t2
     #22-6*t223-6*t225+t228+12*t188+t233-t235+t237
      t345 = -t239-t243+t245+t251+t253-t255-t257+t259-t261-t264+t266-t27
     #2+t274+t276+t278-t280+t282-t284
      A2e0c4p1t1s = -6*t323*(t344+t345)*t22*t159
      A2e0c1p1t2s = A2e0c1p1t1s
      A2e0c2p1t2s = A2e0c2p1t1s
      A2e0c3p1t2s = A2e0c3p1t1s
      A2e0c4p1t2s = A2e0c4p1t1s

      RETURN
      END

      SUBROUTINE CrsecSeqc_mls()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      DOUBLE COMPLEX spf10
      DOUBLE COMPLEX spf9
      DOUBLE COMPLEX DENSPENCE
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t1 = 1-x1
      spf9 = DENSPENCE(DCMPLX(t1),spcZip)
      t2 = 1-x2
      spf10 = DENSPENCE(DCMPLX(t2),spcZip)
      t3 = CsTh**2
      t4 = t3*j4
      t5 = t4*x1
      t7 = 1/(-1+cz1)
      t9 = log(DCMPLX(-cz1*t7))
      t14 = cz1**2
      t19 = 1/cz1
      if (DISonQED .eq. 1) then
      t22 = (4*t9-3-4*cz1*t9+10*cz1+2*t14*t9-4*t14)*t19/2
      else
      t22 = 0
      end if
      t23 = t22*cz1
      t24 = t23*cz2
      t26 = t22*t14
      t27 = t26*cz2
      t29 = cz2**2
      t30 = t23*t29
      t32 = t26*t29
      t34 = t4*x2
      t36 = 1/(-1+cz2)
      t38 = log(DCMPLX(-cz2*t36))
      t47 = 1/cz2
      if (DISonQED .eq. 1) then
      t50 = (4*t38-3-4*cz2*t38+10*cz2+2*t29*t38-4*t29)*t47/2
      else
      t50 = 0
      end if
      t51 = t50*cz2
      t52 = t51*cz1
      t54 = t50*t29
      t55 = t54*cz1
      t57 = t51*t14
      t59 = t54*t14
      t61 = cz1*cz2
      t63 = t14*cz2
      t65 = cz1*t29
      t67 = t14*t29
      t70 = log(DCMPLX(1/s))
      t71 = log(DCMPLX(dels))
      t72 = t70*t71
      t73 = t3*cz1
      t77 = t3*t14
      t87 = -t5*t24+t5*t27+t5*t30-t5*t32-t34*t52+t34*t55+t34*t57-t34*t59
     #-6*t61+6*t63+6*t65-6*t67+4*t72*t73*cz2-4*t72*t77*cz2-4*t72*t73*t29
     #+4*t72*t77*t29
      t88 = muf**2
      t89 = 1/t88
      t91 = log(DCMPLX(s*t89))
      t92 = t3*t91
      t93 = t71*cz1
      t94 = t93*cz2
      t97 = t71*t14
      t98 = t97*cz2
      t101 = t93*t29
      t104 = t97*t29
      t107 = j4*x1
      t108 = t107*t24
      t109 = t107*t27
      t112 = t107*t30
      t113 = t107*t32
      t116 = j4*x2
      t117 = t116*t52
      t118 = t116*t55
      t121 = 4*t92*t94-4*t92*t98-4*t92*t101+4*t92*t104-t108+t109+t4*t24-
     #t4*t27+t112-t113-t4*t30+t4*t32-t117+t118+t4*t52-t4*t55
      t123 = t116*t57
      t124 = t116*t59
      t127 = delc*s
      t131 = log(DCMPLX(-t127*cz1*t89*t7))
      t132 = x1*t131
      t142 = log(DCMPLX(-t127*cz2*t89*t36))
      t143 = x2*t142
      t150 = -t2
      t155 = log(DCMPLX(-t127*t150*t89/x2))
      t156 = t155**2
      t157 = t156*t29
      t158 = t157*cz1
      t159 = -t1
      t164 = log(DCMPLX(-t127*t159*t89/x1))
      t165 = t164**2
      t166 = t165*cz1
      t167 = t166*t29
      t168 = j4*t131
      t170 = 2*t168*t29
      t171 = PI**2
      t172 = t171/3
      t173 = x1**2
      t175 = log(DCMPLX(t1))
      t180 = log(DCMPLX(x1))
      t185 = t175**2
      if (DISonQED .eq. 1) then
      t187 = t172+t173+7.D0/2.D0*x1+3*t175-t175*x1-t175*t173/2+t180*x1+t
     #180*t173/2-2*spf9-t185
      else
      t187 = 0
      end if
      t188 = t187*t14
      t189 = t188*cz2
      t190 = t156*cz2
      t191 = t190*cz1
      t192 = t166*cz2
      t194 = t157*t14
      t195 = t123-t124-t4*t57+t4*t59+2*t4*t132*cz2-2*t4*t132*t29+2*t4*t1
     #43*cz1-2*t4*t143*t14-t158-t167+t170-t189+t191+t192-8*t94+t194
      t199 = t70*cz1
      t201 = 3*t199*cz2
      t202 = t70*t14
      t204 = 3*t202*cz2
      t206 = 3*t199*t29
      t208 = 3*t202*t29
      t209 = t165*t14
      t210 = t209*cz2
      t211 = x2**2
      t213 = log(DCMPLX(t2))
      t218 = log(DCMPLX(x2))
      t223 = t213**2
      if (DISonQED .eq. 1) then
      t225 = t172+t211+7.D0/2.D0*x2+3*t213-t213*x2-t213*t211/2+t218*x2+t
     #218*t211/2-2*spf10-t223
      else
      t225 = 0
      end if
      t226 = t225*cz2
      t227 = t226*t14
      t228 = t209*t29
      t229 = t190*t14
      t231 = 2*t168*cz2
      t232 = t91*cz1
      t234 = 3*t232*cz2
      t235 = t91*t14
      t237 = 3*t235*cz2
      t239 = 3*t232*t29
      t241 = 3*t235*t29
      t242 = 8*t98+8*t101-8*t104+t201-t204-t206+t208-t210-t227+t228-t229
     #-t231+t234-t237-t239+t241
      t247 = log(DCMPLX(t127*t89*dels))
      t248 = t247**2
      t249 = t248*cz1
      t251 = 2*t249*cz2
      t252 = t248*t14
      t254 = 2*t252*cz2
      t256 = 2*t249*t29
      t258 = 2*t252*t29
      t259 = j4*t142
      t261 = 2*t259*cz1
      t263 = 2*t259*t14
      t264 = t188*t29
      t265 = t187*cz1
      t266 = t265*cz2
      t267 = t265*t29
      t268 = t225*t29
      t269 = t268*cz1
      t270 = t268*t14
      t271 = t226*cz1
      t272 = t3*t70
      t281 = -t251+t254+t256-t258-t261+t263+t264+t266-t267-t269+t270+t27
     #1+3*t272*t61-3*t272*t63-3*t272*t65+3*t272*t67
      t283 = 4*t72*t61
      t285 = 4*t72*t63
      t287 = 4*t72*t65
      t289 = 4*t72*t67
      t298 = t91*t71
      t300 = 4*t298*t61
      t302 = 4*t298*t63
      t304 = 4*t298*t65
      t306 = 4*t298*t67
      t307 = t3*t248
      t316 = t283-t285-t287+t289+3*t92*t61-3*t92*t63-3*t92*t65+3*t92*t67
     #+t300-t302-t304+t306-2*t307*t61+2*t307*t63+2*t307*t65-2*t307*t67
      t318 = j4*t22
      t319 = t318*t61
      t320 = t318*t63
      t321 = t318*t65
      t322 = t318*t67
      t323 = j4*t50
      t324 = t323*t61
      t325 = t323*t65
      t326 = t323*t63
      t327 = t323*t67
      t328 = t3*t165
      t333 = t3*t156
      t338 = t319-t320-t321+t322+t324-t325-t326+t327+t328*t61-t328*t63-t
     #328*t65+t328*t67+t333*t61-t333*t65-t333*t63+t333*t67
      t339 = t131*cz2
      t341 = 2*t107*t339
      t344 = t131*t29
      t346 = 2*t107*t344
      t349 = t142*cz1
      t351 = 2*t116*t349
      t354 = t142*t14
      t356 = 2*t116*t354
      t359 = t3*t187
      t364 = t3*t225
      t369 = t341-2*t4*t339-t346+2*t4*t344+t351-2*t4*t349-t356+2*t4*t354
     #-t359*t63+t359*t61+t359*t67-t359*t65+t364*t61-t364*t63+t364*t67-t3
     #64*t65
      t374 = (t87+t121+t195+t242+t281+t316+t338+t369)*t19*t7
      t375 = 1/t171
      A1e0c1p1t1c = 3.D0/2.D0*t374*t375*t47*t36
      t379 = 2*t171
      t380 = t71**2
      t381 = 6*t380
      t387 = log(DCMPLX(delc))
      t389 = 12*t387*t71
      t390 = 9*t387
      t391 = 12*t72
      t392 = 9*t70
      t393 = t3*t387
      t401 = -3+t379+t381-21*t3+2*t171*t3+6*t380*t3+t389+t390-t391-t392+
     #12*t393*t71+9*t393-12*t72*t3-9*t272+24*t71
      A1e0c2p1t1c = -t401*t375/2
      A1e0c3p1t1c = 0
      t404 = s**2
      t405 = t404*Mzsq
      t408 = 1/(-s+Mzsq)
      t409 = t375*t408
      t411 = abs(DCMPLX(Mzsq))
      t412 = t411**2
      t414 = 1/(s*Mzsq-t412)
      t416 = t414*t47*t36
      A1e0c4p1t1c = -3.D0/2.D0*t405*t374*t409*t416
      t421 = t409*t414
      A1e0c5p1t1c = t405*t401*t421/2
      t423 = conjg(Mzsq)
      t426 = s*(Mzsq+t423-2*s)
      t429 = 1/(-s+t423)
      t430 = t375*t429
      t432 = t408*t47*t36
      A1e0c6p1t1c = -3.D0/2.D0*t426*t374*t430*t432
      t437 = t430*t408
      A1e0c7p1t1c = t426*t401*t437/2
      A1e0c1p1t2c = A1e0c1p1t1c
      A1e0c2p1t2c = A1e0c2p1t1c
      A1e0c3p1t2c = 0
      A1e0c4p1t2c = A1e0c4p1t1c
      A1e0c5p1t2c = A1e0c5p1t1c
      A1e0c6p1t2c = A1e0c6p1t1c
      A1e0c7p1t2c = A1e0c7p1t1c
      t444 = 2*t142-2*t142*cz2+t142*t29+t29-t51
      t445 = j4*t150*t444
      t446 = 1+t3
      t447 = t446*t47
      A1e0c1p2t1c = 3.D0/2.D0*t445*t447*t36*t375
      A1e0c2p2t1c = 0
      A1e0c3p2t1c = 0
      t452 = t447*t36
      A1e0c4p2t1c = -3.D0/2.D0*t405*t445*t452*t421
      A1e0c5p2t1c = 0
      A1e0c6p2t1c = -3.D0/2.D0*t426*t445*t452*t437
      A1e0c7p2t1c = 0
      A1e0c1p2t2c = A1e0c1p2t1c
      A1e0c2p2t2c = 0
      A1e0c3p2t2c = 0
      A1e0c4p2t2c = A1e0c4p2t1c
      A1e0c5p2t2c = 0
      A1e0c6p2t2c = A1e0c6p2t1c
      A1e0c7p2t2c = 0
      t465 = 2*t131-2*t131*cz1+t131*t14+t14-t23
      t466 = j4*t159*t465
      t467 = t446*t19
      t468 = t7*t375
      A1e0c1p3t1c = 3.D0/2.D0*t466*t467*t468
      A1e0c2p3t1c = 0
      A1e0c3p3t1c = 0
      t472 = t467*t7
      A1e0c4p3t1c = -3.D0/2.D0*t405*t466*t472*t421
      A1e0c5p3t1c = 0
      A1e0c6p3t1c = -3.D0/2.D0*t426*t466*t472*t437
      A1e0c7p3t1c = 0
      A1e0c1p3t2c = A1e0c1p3t1c
      A1e0c2p3t2c = 0
      A1e0c3p3t2c = 0
      A1e0c4p3t2c = A1e0c4p3t1c
      A1e0c5p3t2c = 0
      A1e0c6p3t2c = A1e0c6p3t1c
      A1e0c7p3t2c = 0
      t484 = -9*t61+9*t63+9*t65-9*t67-t108+t109+t112-t113-t117+t118+t123
     #-t124-t158-t167+t170-t189+t191
      t489 = t192-12*t94+t194+12*t98+12*t101-12*t104+t201-t204-t206+t208
     #-t210-t227+t228-t229-t231+t234-t237
      t491 = -t239+t241-t251+t254+t256-t258-t261+t263+t264+t266-t267-t26
     #9+t270+t271+t283-t285-t287
      t492 = t289+t300-t302-t304+t306+t319-t320-t321+t322+t324-t325-t326
     #+t327+t341-t346+t351-t356
      t496 = CsTh*(t484+t489+t491+t492)*t19
      A2e0c1p1t1c = -12*t405*t496*t468*t408*t416
      t502 = Mzsq*CsTh
      t506 = (6+t379+t381+t389+t390-t391-t392+36*t71)*t375
      A2e0c2p1t1c = 4*t502*t404*t506*t408*t414
      A2e0c3p1t1c = 0
      A2e0c4p1t1c = -3*t426*t496*t468*t429*t432
      A2e0c5p1t1c = t426*CsTh*t506*t429*t408
      A2e0c1p1t2c = A2e0c1p1t1c
      A2e0c2p1t2c = A2e0c2p1t1c
      A2e0c3p1t2c = 0
      A2e0c4p1t2c = A2e0c4p1t1c
      A2e0c5p1t2c = A2e0c5p1t1c
      t518 = t404*j4
      t522 = t444*t47*t36
      A2e0c1p2t1c = -12*t502*t518*t150*t522*t421
      A2e0c2p2t1c = 0
      A2e0c3p2t1c = 0
      t526 = CsTh*j4
      A2e0c4p2t1c = -3*t426*t526*t150*t522*t437
      A2e0c5p2t1c = 0
      A2e0c1p2t2c = A2e0c1p2t1c
      A2e0c2p2t2c = 0
      A2e0c3p2t2c = 0
      A2e0c4p2t2c = A2e0c4p2t1c
      A2e0c5p2t2c = 0
      t535 = t465*t19*t7
      A2e0c1p3t1c = -12*t502*t518*t159*t535*t421
      A2e0c2p3t1c = 0
      A2e0c3p3t1c = 0
      A2e0c4p3t1c = -3*t426*t526*t159*t535*t437
      A2e0c5p3t1c = 0
      A2e0c1p3t2c = A2e0c1p3t1c
      A2e0c2p3t2c = 0
      A2e0c3p3t2c = 0
      A2e0c4p3t2c = A2e0c4p3t1c
      A2e0c5p3t2c = 0

      RETURN
      END

      SUBROUTINE CrsecSeqv_mls()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      DOUBLE COMPLEX spf10
      DOUBLE COMPLEX spf11
      DOUBLE COMPLEX spf12
      DOUBLE COMPLEX spf13
      DOUBLE COMPLEX spf16
      DOUBLE COMPLEX spf17
      DOUBLE COMPLEX spf18
      DOUBLE COMPLEX spf19
      DOUBLE COMPLEX spf20
      DOUBLE COMPLEX spf21
      DOUBLE COMPLEX spf22
      DOUBLE COMPLEX spf23
      DOUBLE COMPLEX spf24
      DOUBLE COMPLEX spf25
      DOUBLE COMPLEX spf26
      DOUBLE COMPLEX spf27
      DOUBLE COMPLEX spf28
      DOUBLE COMPLEX spf29
      DOUBLE COMPLEX spf30
      DOUBLE COMPLEX spf31
      DOUBLE COMPLEX spf32
      DOUBLE COMPLEX spf33
      DOUBLE COMPLEX spf34
      DOUBLE COMPLEX spf35
      DOUBLE COMPLEX spf36
      DOUBLE COMPLEX spf37
      DOUBLE COMPLEX spf4
      DOUBLE COMPLEX spf5
      DOUBLE COMPLEX spf6
      DOUBLE COMPLEX spf7
      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX spf9
      DOUBLE COMPLEX vI100pc1
      DOUBLE COMPLEX vI101pc1
      DOUBLE COMPLEX vI102pc1
      DOUBLE COMPLEX vI103pc1
      DOUBLE COMPLEX vI104pc1
      DOUBLE COMPLEX vI105pc1
      DOUBLE COMPLEX vI106pc1
      DOUBLE COMPLEX vI107pc1
      DOUBLE COMPLEX vI108pc1
      DOUBLE COMPLEX vI109pc1
      DOUBLE COMPLEX vI11pc1
      DOUBLE COMPLEX vI25pc1
      DOUBLE COMPLEX vI38pc1
      DOUBLE COMPLEX vI39pc1
      DOUBLE COMPLEX vI40pc1
      DOUBLE COMPLEX vI41pc1
      DOUBLE COMPLEX vI42pc1
      DOUBLE COMPLEX vI43pc1
      DOUBLE COMPLEX vI44pc1
      DOUBLE COMPLEX vI45pc1
      DOUBLE COMPLEX vI46pc1
      DOUBLE COMPLEX vI47pc1
      DOUBLE COMPLEX vI48pc1
      DOUBLE COMPLEX vI49pc1
      DOUBLE COMPLEX vI50pc1
      DOUBLE COMPLEX vI51pc1
      DOUBLE COMPLEX vI52pc1
      DOUBLE COMPLEX vI53pc1
      DOUBLE COMPLEX vI54pc1
      DOUBLE COMPLEX vI55pc1
      DOUBLE COMPLEX vI56pc1
      DOUBLE COMPLEX vI57pc1
      DOUBLE COMPLEX vI58pc1
      DOUBLE COMPLEX vI59pc1
      DOUBLE COMPLEX vI60pc1
      DOUBLE COMPLEX vI61pc1
      DOUBLE COMPLEX vI62pc1
      DOUBLE COMPLEX vI63pc1
      DOUBLE COMPLEX vI64pc1
      DOUBLE COMPLEX vI65pc1
      DOUBLE COMPLEX vI66pc1
      DOUBLE COMPLEX vI67pc1
      DOUBLE COMPLEX vI68pc1
      DOUBLE COMPLEX vI69pc1
      DOUBLE COMPLEX vI70pc1
      DOUBLE COMPLEX vI71pc1
      DOUBLE COMPLEX vI72pc1
      DOUBLE COMPLEX vI73pc1
      DOUBLE COMPLEX vI74pc1
      DOUBLE COMPLEX vI75pc1
      DOUBLE COMPLEX vI76pc1
      DOUBLE COMPLEX vI77pc1
      DOUBLE COMPLEX vI78pc1
      DOUBLE COMPLEX vI79pc1
      DOUBLE COMPLEX vI80pc1
      DOUBLE COMPLEX vI81pc1
      DOUBLE COMPLEX vI82pc1
      DOUBLE COMPLEX vI83pc1
      DOUBLE COMPLEX vI84pc1
      DOUBLE COMPLEX vI85pc1
      DOUBLE COMPLEX vI86pc1
      DOUBLE COMPLEX vI87pc1
      DOUBLE COMPLEX vI88pc1
      DOUBLE COMPLEX vI89pc1
      DOUBLE COMPLEX vI90pc1
      DOUBLE COMPLEX vI91pc1
      DOUBLE COMPLEX vI92pc1
      DOUBLE COMPLEX vI93pc1
      DOUBLE COMPLEX vI94pc1
      DOUBLE COMPLEX vI95pc1
      DOUBLE COMPLEX vI96pc1
      DOUBLE COMPLEX vI97pc1
      DOUBLE COMPLEX vI98pc1
      DOUBLE COMPLEX vI99pc1
      DOUBLE COMPLEX vV1
      DOUBLE COMPLEX vV2
      DOUBLE COMPLEX vV3
      DOUBLE COMPLEX vV4
      DOUBLE COMPLEX vV5
      DOUBLE COMPLEX vV6
      DOUBLE COMPLEX vV7
      DOUBLE COMPLEX vV8
      DOUBLE COMPLEX vV9
      DOUBLE COMPLEX DENSPENCE
      DOUBLE COMPLEX IM
      PARAMETER (IM=(0.D0,1.D0))
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t1 = -Mzsq+Ieps
      t2 = 1/s
      t3 = t1*t2
      spf7 = log(DCMPLX(-t3))
      t4 = -Mzsq+s+Ieps
      spf10 = DENSPENCE(DCMPLX(s/t4),spcZip)
      spf13 = log(DCMPLX(-t4*t2))
      t8 = 1/t1
      spf19 = log(DCMPLX(t4*t8))
      t18 = spf7**2
      t21 = spf13**2
      t26 = spf19*s
      t36 = log(DCMPLX(1/Mzsq))
      t37 = t36*s
      t39 = t36*spf19
      t44 = PI**2
      t46 = t36**2
      t49 = 24*s+6*Mzsq*spf10+12*spf19*Mzsq+6*Mzsq*spf7*spf19+3*Mzsq*t18
     #-3*Mzsq*t21-6*s*spf10-12*t26-6*s*spf7*spf19-3*s*t18+3*s*t21+12*t37
     #+6*t39*Mzsq-6*t39*s+s*t44+3*t46*s
      vV1 = t49*t2/144
      spf4 = log(DCMPLX(t1))
      t51 = 1+CsTh
      t52 = 1/t51
      spf21 = log(DCMPLX(2*t3*t52))
      t55 = s*t51
      t57 = t55*t8/2
      spf26 = log(DCMPLX(-t57))
      spf27 = log(DCMPLX(t57))
      t58 = 2*Mzsq
      t59 = 2*Ieps
      t60 = s*CsTh
      t61 = -t58+t59+s+t60
      spf33 = log(DCMPLX(2*t1/t61))
      t66 = t61*t8/2
      spf35 = log(DCMPLX(t66))
      spf36 = DENSPENCE(DCMPLX(t66),spcZip)
      t78 = log(DCMPLX(t55/2))
      vV2 = -(6*spf33*spf4+t44-6*spf36-6*spf35*spf26+6*IM*spf35*PI+6*spf
     #35*spf27-6*spf33*t78-6*spf33*spf21)*t2*t52/18
      t87 = -1+CsTh
      spf20 = log(DCMPLX(-2*t3/t87))
      t93 = s*t87*t8/2
      spf24 = log(DCMPLX(t93))
      spf25 = log(DCMPLX(-t93))
      t94 = t58-t59-s+t60
      spf32 = log(DCMPLX(-2*t1/t94))
      t99 = t94*t8/2
      spf34 = log(DCMPLX(-t99))
      spf37 = DENSPENCE(DCMPLX(-t99),spcZip)
      t110 = -t87
      t113 = log(DCMPLX(s*t110/2))
      vV3 = -(6*spf32*spf4+t44-6*spf37-6*spf34*spf24+6*IM*spf34*PI+6*spf
     #34*spf25-6*spf32*t113-6*spf32*spf20)*t2/t110/18
      vV4 = vV1
      vV5 = vV2
      vV6 = vV3
      vV7 = vV4
      vV8 = vV7
      vV9 = vV5
      spf6 = log(DCMPLX(-t8))
      t123 = s**2
      t124 = Mzsq*t123
      t125 = t124*CsTh
      t127 = 1/(-s+Mzsq)
      t128 = s*Mzsq
      t129 = abs(DCMPLX(Mzsq))
      t130 = t129**2
      t132 = 1/(t128-t130)
      t133 = t127*t132
      t134 = 1/t44
      t137 = t125*t133*t134*spf6
      t138 = t133*t134
      t139 = t125*t138
      vI11pc1 = -t137/6-t139/6
      t141 = t134*t127
      t142 = t141*spf6
      t143 = t60*t142
      t144 = t60*t141
      vI25pc1 = -t143/6-t144/6
      t146 = s+Ieps
      spf5 = log(DCMPLX(-1/t146))
      t148 = t134/2
      t149 = t134*spf5
      t151 = spf5**2
      t152 = t134*t151
      t154 = CsTh**2
      t155 = t134*t154
      t157 = t149*t154
      t159 = t155*t151
      vI38pc1 = t148+t149/4+t152/16+t155/2+t157/4+t159/16
      vI39pc1 = t148+5.D0/12.D0*t149+t157/12
      t163 = s+t60-t59
      t164 = 1/t163
      spf12 = log(DCMPLX(t164))
      spf22 = log(DCMPLX(-t146*t164))
      t166 = log(DCMPLX(2.D0))
      t167 = t134*t166
      t168 = t167*spf12
      t170 = spf12**2
      t171 = t134*t170
      t173 = t152/6
      t176 = spf22**2
      t179 = CsTh*t134
      t180 = t166*spf12
      t181 = t179*t180
      t183 = t179*t170
      t185 = t179*t151
      t186 = t185/6
      t187 = t166*spf22
      t192 = CsTh/6
      t197 = t159/6
      t202 = t154/6
      vI40pc1 = -t168/3-t171/6-t173+t167*spf22/3+t134*t176/6+1.D0/6.D0-t
     #181/3-t183/6-t186+t179*t187/3+t179*t176/6+t192-t155*t180/3-t155*t1
     #70/6-t197+t155*t187/3+t155*t176/6+t202
      t203 = t134*spf12
      vI41pc1 = t167+t203+t149
      t204 = -s+t60+t59
      t205 = 1/t204
      spf11 = log(DCMPLX(-t205))
      spf23 = log(DCMPLX(t146*t205))
      t207 = t167*spf11
      t209 = spf11**2
      t210 = t134*t209
      t214 = spf23**2
      t217 = t166*spf11
      t218 = t179*t217
      t220 = t179*t209
      t222 = t166*spf23
      vI42pc1 = -t207/3-t210/6-t173+t167*spf23/3+t134*t214/6+1.D0/6.D0+t
     #218/3+t220/6+t186-t179*t222/3-t179*t214/6-t192-t155*t217/3-t155*t2
     #09/6-t197+t155*t222/3+t155*t214/6+t202
      t235 = t134*spf11
      vI43pc1 = t167+t235+t149
      t236 = t179*spf5
      vI44pc1 = -t179-t236/2-t185/8
      t239 = t179*t166
      t240 = t239/2
      t241 = t179*spf12
      t243 = t166**2
      t245 = t179*t243/8
      vI45pc1 = -t179-t240-t241/2-t245-t181/4-t183/8
      t248 = t179*spf11
      vI46pc1 = -t179-t240-t248/2-t245-t218/4-t220/8
      vI47pc1 = 3.D0/16.D0*t239+3.D0/32.D0*t241+5.D0/16.D0*t236+t179+3.D
     #0/32.D0*t248
      t256 = t124*t138
      t257 = t256/2
      t258 = t124*t127
      t259 = t132*t134
      t260 = t259*spf5
      t261 = t258*t260
      t263 = t259*t151
      t266 = t259*t154
      t271 = t258*t259*spf5*t154
      t273 = t154*t151
      vI48pc1 = -t257-t261/4-t258*t263/16-t258*t266/2-t271/4-t258*t259*t
     #273/16
      vI49pc1 = -t257-5.D0/12.D0*t261-t271/12
      spf9 = log(DCMPLX(t163))
      t279 = -Mzsq+t59+s
      spf16 = DENSPENCE(DCMPLX(t146/t279),spcZip)
      spf17 = log(DCMPLX(t279*t8))
      spf28 = log(DCMPLX(-t163*t8))
      spf31 = DENSPENCE(DCMPLX((s+t60-t58)*t164),spcZip)
      t286 = Mzsq**2
      t287 = t286*t127
      t288 = t287*t132
      t289 = t134*s
      t290 = spf28**2
      t292 = t288*t289*t290
      t294 = t286*Mzsq
      t295 = t294*t127
      t296 = t295*t132
      t298 = t296*t167*spf28
      t300 = t132*s
      t301 = t287*t300
      t302 = t301/432
      t303 = t127*Mzsq
      t304 = t132*t123
      t306 = t303*t304/432
      t307 = t287*t259
      t308 = spf17*t166
      t310 = t307*t60*t308
      t311 = t310/36
      t312 = t303*t259
      t313 = t123*CsTh
      t315 = t312*t313*t308
      t316 = t315/36
      t317 = t123*t154
      t320 = t312*t317*t308/36
      t321 = spf17*spf9
      t323 = t307*t60*t321
      t326 = t312*t313*t321
      t331 = t166*spf28
      t333 = t307*t60*t331
      t336 = t312*t313*t331
      t341 = t292/144+t298/72+t302-t306+t311-t316-t320-t323/36+t326/36+t
     #312*t317*t321/36-t333/72+t336/72+t312*t317*t331/72
      t342 = t166*spf9
      t344 = t307*t60*t342
      t347 = t312*t313*t342
      t352 = t296/432
      t353 = t303*t132
      t354 = t134*t123
      t355 = CsTh*spf16
      t357 = t353*t354*t355
      t358 = t357/18
      t362 = t353*t354*t154*spf16/18
      t364 = t288*t289*t355
      t365 = t364/18
      t366 = CsTh*t290
      t368 = t353*t354*t366
      t375 = t288*t289*t366
      t381 = t288*t289*t331
      t383 = CsTh*spf31
      t385 = t353*t354*t383
      t391 = t344/72-t347/72-t312*t317*t342/72-t352-t358-t362+t365-t368/
     #144-t353*t354*t154*t290/144+t375/144+t353*t354*t331/72-t381/72-t38
     #5/72-t353*t354*t154*spf31/72
      t394 = t288*t289*t383
      t396 = spf9**2
      t397 = CsTh*t396
      t399 = t353*t354*t397
      t406 = t288*t289*t397
      t412 = t288*t289*t342
      t418 = t288*t289*t308
      t419 = t418/36
      t422 = t353*t354*t308/36
      t424 = t288*t289*t321
      t427 = t295*t259*t290
      t430 = t295*t259*spf16
      t431 = t430/18
      t433 = t295*t259*spf31
      t435 = t394/72+t399/144+t353*t354*t154*t396/144-t406/144-t353*t354
     #*t342/72+t412/72+t353*t354*t321/36+t419-t422-t424/36-t427/144-t431
     #-t433/72
      t437 = t295*t259*t396
      t440 = t303*t304*CsTh
      t441 = t440/432
      t444 = t303*t304*t154/432
      t446 = t287*t300*CsTh
      t447 = t446/432
      t448 = t134*spf17
      t450 = t296*t448*spf9
      t456 = t288*t289*t396
      t459 = t296*t167*spf9
      t462 = t296*t448*t166
      t463 = t462/36
      t468 = t288*t289*spf31
      t471 = t288*t289*spf16
      t472 = t471/18
      t475 = t353*t354*spf16/18
      t479 = t437/144-t441-t444+t447+t450/36+t353*t354*t396/144-t456/144
     #-t459/72-t463-t353*t354*spf31/72+t468/72+t472-t475-t353*t354*t290/
     #144
      vI50pc1 = t341+t391+t435+t479
      t481 = t279**2
      t482 = t1**2
      spf18 = log(DCMPLX(t481/t482))
      t485 = t259*t166
      t486 = t295*t485
      t487 = t486/12
      t490 = t295*t259*spf9/12
      t493 = t295*t259*spf18/12
      t495 = t288*t289*t166
      t496 = t495/18
      t498 = t288*t289*spf9
      t501 = t288*t289*spf18
      t502 = t501/18
      t504 = t353*t354*t166
      t505 = t504/12
      t507 = t353*t354*spf9
      t510 = t353*t354*spf18
      t511 = t510/12
      vI51pc1 = t487-t490-t493-t496+t498/18+t502+t505-t507/12-t511
      spf8 = log(DCMPLX(-t204))
      spf29 = log(DCMPLX(t204*t8))
      spf30 = DENSPENCE(DCMPLX((-s+t60+t58)*t205),spcZip)
      t516 = spf17*spf8
      t520 = t166*spf8
      t522 = t307*t60*t520
      t525 = t312*t313*t520
      t531 = t307*t60*t516
      t534 = t312*t313*t516
      t536 = t166*spf29
      t538 = t307*t60*t536
      t541 = t312*t313*t536
      t543 = t302-t306-t311+t316-t320+t312*t317*t516/36-t522/72+t525/72-
     #t312*t317*t520/72+t531/36-t534/36+t538/72-t541/72
      t547 = spf29**2
      t548 = CsTh*t547
      t550 = t353*t354*t548
      t557 = t288*t289*t548
      t560 = t288*t289*t536
      t565 = CsTh*spf30
      t567 = t353*t354*t565
      t574 = t288*t289*t565
      t576 = spf8**2
      t577 = CsTh*t576
      t579 = t353*t354*t577
      t586 = t288*t289*t577
      t589 = t288*t289*t520
      t594 = t312*t317*t536/72+t550/144-t353*t354*t154*t547/144-t557/144
     #-t560/72+t353*t354*t536/72+t567/72-t353*t354*t154*spf30/72-t574/72
     #-t579/144+t353*t354*t154*t576/144+t586/144+t589/72-t353*t354*t520/
     #72
      t600 = t288*t289*t516
      t603 = t295*t259*t547
      t605 = t353*t354*t516/36-t600/36-t352+t358-t362-t365+t419-t422-t43
     #1+t441-t444-t447-t603/144
      t607 = t295*t259*spf30
      t610 = t295*t259*t576
      t613 = t288*t289*spf30
      t619 = t296*t448*spf8
      t622 = t288*t289*t576
      t628 = t296*t167*spf8
      t631 = t288*t289*t547
      t637 = t296*t167*spf29
      t639 = -t607/72+t610/144-t463+t472-t475+t613/72-t353*t354*spf30/72
     #+t619/36-t622/144+t353*t354*t576/144-t628/72+t631/144-t353*t354*t5
     #47/144+t637/72
      vI52pc1 = t543+t594+t605+t639
      t643 = t295*t259*spf8/12
      t645 = t288*t289*spf8
      t648 = t353*t354*spf8
      vI53pc1 = t487-t643-t493-t496+t645/18+t502+t505-t648/12-t511
      t651 = spf6**2
      t657 = t139/6
      vI54pc1 = -t440/36-t125*t133*t134*t651/12-t137/6-t657
      t658 = t286*s
      t660 = t134*spf19
      t661 = t133*t660
      t666 = t134*t36
      t671 = t286*CsTh*t127
      t677 = t259*spf19
      vI55pc1 = -t657-t658*CsTh*t661/10+t125*t661/12-t125*t133*t666/12+t
     #671*t259*s/30+t294*CsTh*t127*t677/60+t671*t259*t37/60
      vI56pc1 = t125*vV1*t127*t259
      t685 = t286*t132
      t686 = t685*t134
      t688 = t685*t167/2
      t689 = t685*t203
      t693 = t685*t134*t243/8
      t698 = CsTh*t132
      t699 = t698*t134
      t700 = t128*t699
      t701 = t128*CsTh
      t702 = t701*t485
      t703 = t702/2
      t704 = t259*spf12
      t705 = t701*t704
      t708 = Mzsq*t132*t134
      t709 = t60*t243
      t711 = t708*t709/8
      t713 = t708*t60*t180
      t715 = t60*t170
      t716 = t708*t715
      vI57pc1 = t686+t688+t689/2+t693+t685*t168/4+t685*t171/8-t700-t703-
     #t705/2-t711-t713/4-t716/8
      t718 = t685*t235
      t724 = t259*spf11
      t725 = t701*t724
      t728 = t708*t60*t217
      t730 = t60*t209
      t731 = t708*t730
      vI58pc1 = -t686-t688-t718/2-t693-t685*t207/4-t685*t210/8-t700-t703
     #-t725/2-t711-t728/4-t731/8
      t736 = t124*vV2
      vI59pc1 = -t658*vV2*t132*t134/2+t736*t699/2-t658*vV2*t699/2+t736*t
     #266/2
      t745 = t124*vV3
      vI60pc1 = t658*vV3*t132*t134/2+t745*t699/2-t658*vV3*t699/2-t745*t2
     #66/2
      t753 = t128*t724
      t756 = t128*t704
      vI61pc1 = 7.D0/8.D0*t718-7.D0/8.D0*t689+t700-t753/8+t703+t725/4+t7
     #56/8+t705/4
      t759 = conjg(Mzsq)
      t761 = 1/(-s+t759)
      t762 = t134*t761
      t763 = t762*t127
      t764 = t313*t763
      t765 = t313*t127
      t766 = t762*spf5
      t767 = t765*t766
      t769 = t762*t151
      t772 = t60*t134
      t773 = t761*t127
      t774 = t773*Mzsq
      t775 = t772*t774
      t777 = t60*t127
      t780 = t777*t762*spf5*Mzsq
      t786 = t773*t759
      t787 = t772*t786
      t791 = t777*t762*spf5*t759
      vI62pc1 = -t764-t767/2-t765*t769/8+t775/2+t780/4+t777*t762*Mzsq*t1
     #51/16+t787/2+t791/4+t777*t762*t759*t151/16
      vI63pc1 = -t780/2-5.D0/14.D0*t775-t791/2+5.D0/7.D0*t764-5.D0/14.D0
     #*t787+t767
      t802 = t289*t786
      t804 = t141*s
      t805 = t761*spf5
      t807 = t804*t805*t759
      t809 = t761*t759
      t813 = t354*t773
      t815 = t123*t127
      t816 = t815*t766
      t820 = t289*t774
      t823 = t804*t805*Mzsq
      t825 = t761*Mzsq
      t829 = t289*t761
      t836 = t804*t805*t759*t154
      t846 = t804*t805*Mzsq*t154
      t854 = t815*t134
      t856 = t854*t805*t154
      vI64pc1 = -t802/4-t807/8-t804*t809*t151/32+t813/2+t816/4+t815*t769
     #/16-t820/4-t823/8-t804*t825*t151/32-t829*t127*t759*t154/4-t836/8-t
     #804*t809*t273/32-t829*t303*t154/4-t846/8-t804*t825*t273/32+t354*t7
     #73*t154/2+t856/4+t854*t761*t154*t151/16
      vI65pc1 = -t802/6+t813/3-t820/6+5.D0/18.D0*t816+t856/18-5.D0/36.D0
     #*t807-5.D0/36.D0*t823-t836/36-t846/36
      t871 = Mzsq*t134
      t874 = t871*t127*spf6/6
      t875 = t871*t127
      t876 = t875/6
      vI66pc1 = t874+t876+t141*s*spf6/24+t804/24
      t882 = t701*t263
      t883 = t882/2
      t885 = t701*t259*t187
      t886 = t259*t176
      t887 = t701*t886
      t889 = t128*t698
      t890 = t889/2
      vI67pc1 = t713+t716/2+t883-t885-t887/2-t890
      t891 = t128*t485
      t892 = t891/3
      t894 = t128*t260
      t895 = t894/3
      t896 = 2.D0/3.D0*t702
      t898 = t701*t260
      t899 = 2.D0/3.D0*t898
      vI68pc1 = t892+t756/3+t895-t896-2.D0/3.D0*t705-t899
      t901 = t141*Mzsq*t290
      t904 = t141*Mzsq*spf16
      t905 = t904/12
      t907 = t141*Mzsq*t396
      t910 = t141*Mzsq*spf31
      t912 = t303*CsTh
      t913 = t912/288
      t914 = t777/288
      t916 = t127*t2*t286
      t917 = t916/288
      t918 = Mzsq*t166
      t920 = t141*t918*spf28
      t922 = t2*t286
      t924 = t141*t922*t290
      t927 = t141*t60*t290
      t929 = Mzsq*CsTh
      t931 = t141*t929*t290
      t934 = t141*t922*spf16
      t935 = t934/12
      t937 = t141*t60*spf16
      t938 = t937/12
      t940 = t141*t929*spf16
      t941 = t940/12
      t943 = t141*t922*spf31
      t946 = t141*t60*spf31
      t949 = t141*t929*spf31
      t951 = Mzsq*spf17
      t953 = t141*t951*t166
      t954 = t953/24
      t955 = t901/96+t905-t907/96+t910/48+t913+t914-t917-t920/48-t924/96
     #+t927/96+t931/96-t935+t938+t941-t943/48+t946/48+t949/48+t954
      t957 = t141*t951*spf9
      t960 = t141*t922*t396
      t963 = t141*t60*t396
      t966 = t141*t929*t396
      t969 = t141*t918*spf9
      t971 = CsTh*t166
      t972 = t971*spf28
      t973 = t875*t972
      t975 = t141*t2
      t976 = t286*t166
      t978 = t975*t976*spf28
      t980 = t804*t972
      t982 = CsTh*spf17
      t983 = t982*t166
      t984 = t875*t983
      t985 = t984/24
      t986 = t286*spf17
      t988 = t975*t986*t166
      t989 = t988/24
      t990 = t804*t983
      t991 = t990/24
      t992 = t982*spf9
      t993 = t875*t992
      t996 = t975*t986*spf9
      t998 = t804*t992
      t1000 = t971*spf9
      t1001 = t875*t1000
      t1004 = t975*t976*spf9
      t1006 = t804*t1000
      t1008 = t303/288
      t1009 = -t957/24+t960/96-t963/96-t966/96+t969/48-t973/48+t978/48-t
     #980/48+t985-t989+t991-t993/24+t996/24-t998/24+t1001/48-t1004/48+t1
     #006/48+t1008
      vI69pc1 = t955+t1009
      t1011 = t141*t922*t166
      t1012 = t1011/10
      t1014 = t141*t922*spf9
      t1017 = t141*t922*spf18
      t1018 = t1017/10
      t1019 = t141*t918
      t1020 = t1019/15
      t1022 = t141*Mzsq*spf9
      t1025 = t141*Mzsq*spf18
      t1026 = t1025/15
      t1029 = t141*t929*t166/15
      t1035 = t141*t929*spf18/15
      t1037 = t141*s*t166
      t1038 = t1037/30
      t1040 = t141*s*spf9
      t1043 = t141*s*spf18
      t1044 = t1043/30
      t1045 = t60*t166
      t1047 = 2.D0/15.D0*t141*t1045
      t1053 = 2.D0/15.D0*t141*t60*spf18
      vI70pc1 = -t1012+t1014/10+t1018+t1020-t1022/15-t1026+t1029-t141*t9
     #29*spf9/15-t1035+t1038-t1040/30-t1044+t1047-2.D0/15.D0*t141*t60*sp
     #f9-t1053
      t1056 = t701*t259*t222
      t1057 = t259*t214
      t1058 = t701*t1057
      vI71pc1 = -t728-t731/2-t883+t1056+t1058/2+t890
      vI72pc1 = t892+t753/3+t895+t896+2.D0/3.D0*t725+t899
      t1063 = t141*Mzsq*spf30
      t1066 = t141*Mzsq*t547
      t1069 = t141*Mzsq*t576
      t1072 = t141*t922*t547
      t1075 = t141*t929*t547
      t1078 = t141*t918*spf29
      t1081 = t141*t60*spf30
      t1084 = t141*t922*spf30
      t1087 = t141*t929*spf30
      t1090 = t141*t951*spf8
      t1092 = t905+t1063/48+t1066/96-t1069/96-t913-t914-t917-t935-t938-t
     #941+t954-t1072/96-t1075/96-t1078/48-t1081/48-t1084/48-t1087/48-t10
     #90/24
      t1094 = t141*t60*t576
      t1097 = t141*t922*t576
      t1100 = t141*t929*t576
      t1103 = t141*t918*spf8
      t1106 = t141*t60*t547
      t1108 = t982*spf8
      t1109 = t804*t1108
      t1111 = t875*t1108
      t1114 = t975*t986*spf8
      t1116 = t971*spf8
      t1117 = t875*t1116
      t1120 = t975*t976*spf8
      t1122 = t804*t1116
      t1124 = t971*spf29
      t1125 = t804*t1124
      t1127 = t875*t1124
      t1130 = t975*t976*spf29
      t1132 = t1094/96+t1097/96+t1100/96+t1103/48-t1106/96-t985-t989-t99
     #1+t1109/24+t1111/24+t1114/24-t1117/48-t1120/48-t1122/48+t1125/48+t
     #1127/48+t1130/48+t1008
      vI73pc1 = t1092+t1132
      t1134 = t141*t922*spf8
      t1137 = t141*Mzsq*spf8
      t1143 = t141*s*spf8
      vI74pc1 = -t1012+t1134/10+t1018+t1020-t1137/15-t1026-t1029+t141*t9
     #29*spf8/15+t1035+t1038-t1143/30-t1044-t1047+2.D0/15.D0*t141*t60*sp
     #f8+t1053
      vI75pc1 = t303/36+t871*t127*t651/12+t874+t876
      vI76pc1 = t876+t287*t134*t2*spf19/12-t303*t660/18+t303*t666/12+t80
     #4/18-t141*t26/36+t141*t37/36
      vI77pc1 = Mzsq*vV4*t141
      t1167 = t762*Mzsq/2
      t1169 = t762*t918/4
      t1171 = t762*Mzsq*spf12
      t1175 = t762*Mzsq*t243/16
      t1185 = t762*t2*Mzsq*t759/2
      t1186 = t762*t2
      t1187 = Mzsq*t759
      t1188 = t1187*t166
      t1190 = t1186*t1188/4
      t1191 = t1187*spf12
      t1192 = t1186*t1191
      t1196 = t1186*t1187*t243/16
      t1203 = t1167+t1169+t1171/4+t1175+t762*t918*spf12/8+t762*Mzsq*t170
     #/16-t1185-t1190-t1192/4-t1196-t1186*t1187*t180/8-t1186*t1187*t170/
     #16
      t1204 = t762*t60
      t1205 = t762*t1045
      t1206 = t1205/2
      t1208 = t762*t60*spf12
      t1211 = t762*t709/8
      t1212 = t971*spf12
      t1217 = t759*CsTh
      t1218 = t762*t1217
      t1219 = t1218/2
      t1221 = t762*t1217*t166
      t1222 = t1221/4
      t1224 = t762*t1217*spf12
      t1228 = t762*t1217*t243/16
      t1229 = t762*t759
      t1235 = -t1204-t1206-t1208/2-t1211-t829*t1212/4-t762*t715/8+t1219+
     #t1222+t1224/4+t1228+t1229*t1212/8+t762*t1217*t170/16
      vI78pc1 = t1203+t1235
      t1237 = t762*Mzsq*spf11
      t1245 = t1187*spf11
      t1246 = t1186*t1245
      t1254 = t1167+t1169+t1237/4+t1175+t762*t918*spf11/8+t762*Mzsq*t209
     #/16-t1185-t1190-t1246/4-t1196-t1186*t1187*t217/8-t1186*t1187*t209/
     #16
      t1256 = t762*t60*spf11
      t1258 = t971*spf11
      t1264 = t762*t1217*spf11
      t1271 = t1204+t1206+t1256/2+t1211+t829*t1258/4+t762*t730/8-t1219-t
     #1222-t1264/4-t1228-t1229*t1258/8-t762*t1217*t209/16
      vI79pc1 = t1254+t1271
      t1272 = vV5*t134
      t1276 = s*t154
      vI80pc1 = t1272*Mzsq/2-t1272*t60/2+t1272*t929/2-t1272*t1276/2
      t1279 = vV6*t134
      vI81pc1 = t1279*Mzsq/2+t1279*t60/2-t1279*t929/2-t1279*t1276/2
      t1298 = t761*CsTh
      t1313 = t761*t286
      t1318 = t759*spf11
      t1322 = CsTh*spf11
      t1326 = 4.D0/15.D0*t1186*t287*t759+13.D0/30.D0*t804*t825*t166-t763
     #*t1188/3-t763*t1245/6+13.D0/60.D0*t804*t825*spf11+17.D0/60.D0*t854
     #*t1298*spf11+13.D0/60.D0*t804*t825*spf12-t763*t1191/6-17.D0/60.D0*
     #t854*t1298*spf12-4.D0/15.D0*t813-4.D0/15.D0*t762*t287+t975*t1313*t
     #759*t166/3+t975*t1313*t1318/6+t763*t1187*t1322/10
      t1333 = t759*spf12
      t1337 = CsTh*spf12
      t1350 = t123*t761
      t1357 = t815*spf11
      t1369 = -17.D0/60.D0*t804*t825*t1322-t804*t809*t1322/10+t975*t1313
     #*t1333/6-t763*t1187*t1337/10+17.D0/60.D0*t804*t825*t1337+t804*t809
     #*t1337/10-t141*t1313*t166/3-t1350*t141*t166/10-t141*t1313*spf11/6-
     #t762*t1357/20-t141*t1313*spf12/6-t1350*t141*spf12/20+4.D0/15.D0*t8
     #02-t816/30+t823/30
      vI82pc1 = t1326+t1369
      t1370 = t128*t134
      t1371 = t132*t166
      t1379 = t128*t263/6
      t1386 = t128*t132/6
      t1389 = t882/6
      t1392 = t889/6
      t1393 = t132*t154
      t1402 = t1370*t1393*t151/6
      t1410 = t128*t1393/6
      vI83pc1 = -t1370*t1371*spf12/3-t128*t259*t170/6-t1379+t1370*t1371*
     #spf22/3+t128*t886/6+t1386-t713/3-t716/6-t1389+t885/3+t887/6+t1392-
     #t1370*t1393*t180/3-t1370*t1393*t170/6-t1402+t1370*t1393*t187/3+t13
     #70*t1393*t176/6+t1410
      vI84pc1 = t891/2+t756/2+t894/2
      t1418 = t154*t166
      t1422 = t154*spf17
      t1425 = t804*t1422*t166/36
      t1433 = t904/18
      t1438 = t289*t127*spf16/18
      t1442 = t912/432
      t1443 = -t289*t127*t290/144-t289*t127*spf31/72+t804*t1418*spf28/72
     #-t1425+t804*t1422*spf9/36-t804*t1418*spf9/72+t901/144+t1433-t907/1
     #44+t910/72-t1438+t289*t127*t396/144+t1442
      t1444 = t127*t154
      t1448 = t127*t166
      t1457 = t289*t1444*spf16/18
      t1458 = t127*spf17
      t1461 = t289*t1458*t166/36
      t1471 = t777/432
      t1472 = s*t127
      t1474 = t1472*t154/432
      t1475 = t916/432
      t1479 = -t289*t1444*t290/144+t289*t1448*spf28/72-t289*t1444*spf31/
     #72-t1457-t1461+t289*t1458*spf9/36+t289*t1444*t396/144-t289*t1448*s
     #pf9/72-t1471-t1474-t1475-t920/72-t924/144-t927/144
      t1482 = t934/18
      t1483 = t937/18
      t1484 = t940/18
      t1488 = t953/36
      t1494 = t931/144-t1482-t1483+t1484-t943/72-t946/72+t949/72+t1488-t
     #957/36+t960/144+t963/144-t966/144+t969/72
      t1495 = t1472/432
      t1499 = t984/36
      t1500 = t988/36
      t1501 = t990/36
      t1508 = t303/432
      t1509 = -t1495-t973/72+t978/72+t980/72+t1499-t1500-t1501-t993/36+t
     #996/36+t998/36+t1001/72-t1004/72-t1006/72+t1508
      vI85pc1 = t1443+t1479+t1494+t1509
      t1511 = t1011/12
      t1513 = t1017/12
      t1514 = t1019/18
      t1516 = t1025/18
      t1517 = t1037/12
      t1519 = t1043/12
      vI86pc1 = t1511-t1014/12-t1513-t1514+t1022/18+t1516+t1517-t1040/12
     #-t1519
      vI87pc1 = -t1370*t1371*spf11/3-t128*t259*t209/6-t1379+t1370*t1371*
     #spf23/3+t128*t1057/6+t1386+t728/3+t731/6+t1389-t1056/3-t1058/6-t13
     #92-t1370*t1393*t217/3-t1370*t1393*t209/6-t1402+t1370*t1393*t222/3+
     #t1370*t1393*t214/6+t1410
      vI88pc1 = t891/2+t753/2+t894/2
      t1569 = t804*t1422*spf8/36-t804*t1418*spf8/72+t804*t1418*spf29/72-
     #t289*t127*spf30/72-t289*t127*t547/144+t289*t127*t576/144-t1425+t14
     #33+t1063/72+t1066/144-t1069/144-t1438-t1442
      t1588 = -t1457-t1461+t289*t1458*spf8/36-t289*t1444*spf30/72-t289*t
     #1444*t547/144+t289*t1448*spf29/72+t289*t1444*t576/144-t289*t1448*s
     #pf8/72+t1471-t1474-t1475-t1482+t1483-t1484
      t1602 = t1488-t1072/144-t1075/144-t1078/72+t1081/72-t1084/72-t1087
     #/72-t1090/36-t1094/144+t1097/144+t1100/144+t1103/72+t1106/144
      t1612 = -t1495-t1499-t1500+t1501-t1109/36+t1111/36+t1114/36-t1117/
     #72-t1120/72+t1122/72-t1125/72+t1127/72+t1130/72+t1508
      vI89pc1 = t1569+t1588+t1602+t1612
      vI90pc1 = t1511-t1134/12-t1513-t1514+t1137/18+t1516+t1517-t1143/12
     #-t1519
      t1618 = t141*t651
      t1622 = t144/6
      vI91pc1 = -t777/36-t60*t1618/12-t143/6-t1622
      vI92pc1 = -t700-t898/2-t882/8
      t1625 = CsTh*t127
      t1629 = t660*s
      t1632 = t666*s
      t1639 = t127*spf19
      vI93pc1 = -t1622-t1625*t660*Mzsq/10+t1625*t1629/12-t1625*t1632/12+
     #t179*t303/30+t2*CsTh*t134*t1639*t286/60+t179*t303*t36/60
      vI94pc1 = t60*vV7*t134*t127
      vI95pc1 = -vI78pc1
      vI96pc1 = -t1205/4-t762*s*spf12/40+t762*s*spf11/40-7.D0/40.D0*t123
     #7+7.D0/40.D0*t1171+t1224/20+t1264/20+t1218/5-t1204+7.D0/40.D0*t124
     #6-7.D0/40.D0*t1192-t762*t1318/40+t762*t1333/40-t1208/8-t1256/8+t12
     #21/10-t762*s*spf5*CsTh/4
      vI97pc1 = t139+t125*t133*t149/2+t125*t259*t127*t151/8
      t1683 = 1/(-t124+s*t130+t658-t130*Mzsq)
      vI98pc1 = -5.D0/7.D0*t313*t1683*Mzsq*t134-t313*t1683*t871*spf5
      t1698 = t658*t127*t259*spf6/6
      t1699 = t658*t138
      t1700 = t1699/6
      vI99pc1 = t124*t132*t142/42+t256/42+t1698+t1700
      t1703 = t301/288
      t1704 = t310/24
      t1705 = t315/24
      t1712 = t296/288
      t1713 = t357/12
      t1714 = t364/12
      t1719 = t292/96+t298/48+t1703+t1704+t1705-t323/24-t326/24-t333/48-
     #t336/48+t344/48+t347/48-t1712+t1713+t1714+t368/96+t375/96-t381/48+
     #t385/48
      t1724 = t418/24
      t1727 = t430/12
      t1730 = t440/288
      t1731 = t446/288
      t1735 = t462/24
      t1737 = t471/12
      t1738 = t394/48-t399/96-t406/96+t412/48+t1724-t424/24-t427/96-t172
     #7-t433/48+t437/96+t1730+t1731+t450/24-t456/96-t459/48-t1735+t468/4
     #8+t1737
      vI100pc1 = t1719+t1738
      t1739 = t495/15
      t1741 = t501/15
      t1744 = t288*t289*t971/10
      t1745 = CsTh*spf9
      t1749 = CsTh*spf18
      t1752 = t288*t289*t1749/10
      t1753 = t504/60
      t1755 = t510/60
      t1758 = t353*t354*t971/10
      t1764 = t353*t354*t1749/10
      vI101pc1 = t487-t490-t493-t1739+t498/15+t1741-t1744+t288*t289*t174
     #5/10+t1752-t1753+t507/60+t1755-t1758+t353*t354*t1745/10+t1764
      t1780 = t1703-t1704-t1705-t522/48-t525/48+t531/24+t534/24+t538/48+
     #t541/48-t550/96-t557/96-t560/48-t567/48-t574/48+t579/96+t586/96+t5
     #89/48-t600/24
      t1790 = -t1712-t1713-t1714+t1724-t1727-t1730-t1731-t603/96-t607/48
     #+t610/96-t1735+t1737+t613/48+t619/24-t622/96-t628/48+t631/96+t637/
     #48
      vI102pc1 = t1780+t1790
      t1792 = CsTh*spf8
      vI103pc1 = t487-t643-t493-t1739+t645/15+t1741+t1744-t288*t289*t179
     #2/10-t1752-t1753+t648/60+t1755+t1758-t353*t354*t1792/10-t1764
      vI104pc1 = t301/36+t658*t132*t1618/12+t1698+t1700
      vI105pc1 = t1700+t295*t677/12-t288*t1629/16+t288*t1632/12+t256/24-
     #t708*t1639*t123/48+t708*t127*t36*t123/48
      vI106pc1 = t658*vV8*t138
      vI107pc1 = -vI57pc1
      t1822 = t124*vV9
      vI108pc1 = t658*vV9*t134*t132/2-t1822*t699/2+t658*vV9*t699/2-t1822
     #*t266/2
      vI109pc1 = t294*t132*t141/2+9.D0/16.D0*t486+9.D0/32.D0*t295*t724+9
     #.D0/32.D0*t295*t704-t1699/8+t288*t289*t1322/4-t288*t289*t1337/4-5.
     #D0/16.D0*t288*t289*spf12-5.D0/16.D0*t288*t289*spf11-5.D0/8.D0*t495
     #+t708*t1357/32+3.D0/8.D0*t256-t353*t354*t1322/4+t708*t815*spf12/32
     #+t504/16+t353*t354*t1337/4
      A1e0c1p1t1v = -24*vI38pc1+18*vI39pc1
      A1e0c2p1t1v = -18*vI40pc1-9*vI41pc1+18*vI42pc1+9*vI43pc1+48*vI44pc
     #1+24*vI45pc1+24*vI46pc1+96*vI47pc1
      A1e0c3p1t1v = -24*vI48pc1+18*vI49pc1
      A1e0c4p1t1v = -144*vI11pc1-432*vI50pc1+108*vI51pc1+432*vI52pc1-108
     #*vI53pc1+144*vI54pc1+720*vI55pc1+576*vI56pc1+24*vI57pc1+24*vI58pc1
     #+36*vI59pc1+36*vI60pc1+24*vI61pc1
      A1e0c5p1t1v = -48*vI64pc1+54*vI65pc1
      t1904 = 48*vI79pc1
      t1905 = 36*vI80pc1
      t1906 = 36*vI81pc1
      A1e0c6p1t1v = -144*vI25pc1-18*vI83pc1-18*vI84pc1-432*vI85pc1+108*v
     #I86pc1+18*vI87pc1+18*vI88pc1+432*vI89pc1-108*vI90pc1+144*vI91pc1+4
     #8*vI92pc1+720*vI93pc1+576*vI94pc1+48*vI95pc1+t1904-t1905+t1906+120
     #*vI96pc1
      A1e0c1p1t2v = A1e0c1p1t1v
      A1e0c2p1t2v = A1e0c2p1t1v
      A1e0c3p1t2v = A1e0c3p1t1v
      A1e0c4p1t2v = A1e0c4p1t1v
      A1e0c5p1t2v = A1e0c5p1t1v
      A1e0c6p1t2v = A1e0c6p1t1v
      A2e0c1p1t1v = 96*vI97pc1+84*vI98pc1
      A2e0c2p1t1v = 2016*vI99pc1+1152*vI100pc1-720*vI101pc1+1152*vI102pc
     #1-720*vI103pc1-576*vI104pc1-4608*vI105pc1+2304*vI106pc1+96*vI107pc
     #1+96*vI58pc1+144*vI108pc1+144*vI60pc1+384*vI109pc1
      A2e0c3p1t1v = 48*vI62pc1+42*vI63pc1
      A2e0c4p1t1v = 576*vI66pc1+6*vI67pc1+18*vI68pc1+288*vI69pc1+180*vI7
     #0pc1+6*vI71pc1+18*vI72pc1+288*vI73pc1+180*vI74pc1-144*vI75pc1-1296
     #*vI76pc1+576*vI77pc1+48*vI78pc1+t1904+t1905+t1906+180*vI82pc1
      A2e0c1p1t2v = A2e0c1p1t1v
      A2e0c2p1t2v = A2e0c2p1t1v
      A2e0c3p1t2v = A2e0c3p1t1v
      A2e0c4p1t2v = A2e0c4p1t1v

      RETURN
      END

      SUBROUTINE CrsecSeqn_mls()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      t1 = Eq**2
      t2 = t1*Eq
      t3 = t2*Er
      t4 = El**2
      t5 = t3*t4
      t6 = SnThl**2
      t7 = SnThr**2
      t8 = t6*t7
      t9 = cos(DBLE(phil))
      t10 = t9**2
      t11 = t8*t10
      t13 = 32*t5*t11
      t14 = t1**2
      t15 = t14*Er
      t17 = SnThr*t9
      t20 = 64*t15*SnThl*t17*El
      t21 = t2*El
      t22 = t21*SnThl
      t23 = t17*sqq
      t24 = t22*t23
      t26 = t4*t6
      t27 = t26*t7
      t28 = t10*sqq
      t30 = t27*t28*t1
      t32 = sqq**2
      t33 = t32*El
      t34 = t33*SnThl
      t36 = t34*t17*Eq
      t37 = 4*t36
      t38 = CsThl**2
      t39 = t14*t38
      t40 = CsThr**2
      t41 = t40*sqq
      t42 = t39*t41
      t44 = t32*CsThl
      t45 = t1*CsThr
      t46 = t44*t45
      t47 = 4*t46
      t48 = t40*t2
      t49 = Er**2
      t50 = t49*Er
      t51 = t50*t38
      t53 = 8*t48*t51
      t54 = t1*t49
      t55 = t54*t41
      t57 = t14*t49
      t58 = CsThl*CsThr
      t59 = t57*t58
      t60 = 32*t59
      t61 = t38*t40
      t62 = t57*t61
      t63 = 32*t62
      t64 = t38*sqq
      t65 = t54*t64
      t66 = 2*t65
      t67 = t14*Eq
      t68 = t67*Er
      t69 = t68*t58
      t70 = 64*t69
      t71 = t68*t61
      t72 = 32*t71
      t73 = t14*sqq
      t74 = t58*t73
      t76 = t32*t1
      t77 = 4*t76
      t78 = 8*t73
      t79 = t13+t20-16*t24-8*t30+t37-8*t42-t47+t53-2*t55+t60-t63-t66-t70
     #+t72+16*t74+t77-t78
      t80 = t2*t50
      t81 = 8*t80
      t82 = 32*t68
      t83 = t1*Er
      t84 = t83*t58
      t85 = sqq*El
      t86 = SnThl*SnThr
      t87 = t86*t9
      t88 = t85*t87
      t89 = t84*t88
      t91 = t2*t49
      t92 = t91*t58
      t93 = t9*El
      t94 = t86*t93
      t96 = 32*t92*t94
      t97 = t58*t21
      t98 = t9*sqq
      t99 = t86*t98
      t100 = t97*t99
      t102 = t15*t58
      t104 = 64*t102*t94
      t105 = t40*t1
      t106 = t49*t38
      t107 = t106*sqq
      t109 = 2*t105*t107
      t110 = t58*sqq
      t111 = t54*t110
      t112 = 8*t111
      t113 = t3*t110
      t115 = t61*sqq
      t116 = t3*t115
      t118 = Eq*Er
      t119 = t58*t32
      t121 = 2*t118*t119
      t122 = t32*sqq
      t123 = t48*t50
      t124 = 8*t123
      t125 = t118*t32
      t126 = 6*t125
      t127 = t80*t38
      t128 = 8*t127
      t129 = t54*sqq
      t131 = t3*sqq
      t132 = 16*t131
      t134 = t83*sqq*t94
      t135 = 16*t134
      t136 = t81+t82-8*t89+t96+16*t100-t104-t109-t112+8*t113+8*t116+t121
     #-t122+t124+t126+t128-10*t129-t132-t135
      t137 = t79+t136
      t138 = 1/t49
      t140 = 1/t1
      t141 = sll**2
      t142 = 1/t141
      t145 = 1/(-1+t40)
      A1e0c1p1t1n = 24*t137*t138*t140*t142*t145
      t148 = sll*El
      t149 = t148*t87
      t152 = t2*sll
      t153 = t152*t58
      t156 = CsThr*Eq
      t157 = t156*Er
      t158 = sqq*CsThl
      t159 = t158*sll
      t160 = t157*t159
      t162 = t152*El
      t165 = sll*t1
      t169 = 4*t65
      t170 = t38*t1
      t171 = t49*sll
      t172 = t170*t171
      t174 = t40*sll
      t175 = t54*t174
      t177 = sll*t14
      t182 = t13+t20+8*t84*t149-16*t153*t94-2*t160+16*t162*t87+8*t165*t4
     #*t11+t53+t60-t63+t169-t70+t72-2*t172+2*t175+8*t177*t61-16*t177*t58
     #-t81+t82
      t183 = 8*t177
      t184 = sll*t32
      t185 = 16*t113
      t186 = t49*t40
      t190 = t40*Er
      t191 = t38*t2
      t192 = t191*sll
      t195 = CsThr*t2
      t196 = Er*CsThl
      t197 = t196*sll
      t198 = t195*t197
      t200 = t165*t110
      t202 = sqq*sll
      t203 = t202*t1
      t204 = 4*t203
      t205 = t54*sll
      t207 = 2*t125
      t208 = 4*t129
      t210 = t202*Eq*t94
      t212 = t183+t184+t96-t104-t112+t185+2*t170*t186*sll-8*t190*t192+8*
     #t198+4*t200-t204-2*t205+t124+t207-t128+t208-t132-t135-4*t210
      t213 = t182+t212
      t215 = 1/t32
      t218 = 1/(-1+t38)
      A1e0c2p1t1n = -24*t213*t138*t140*t215*t218
      t229 = t9*t40
      t233 = t14*t1
      t234 = 32*t233
      t235 = t38*CsThl
      t236 = t49*t235
      t240 = t40*CsThr
      t241 = t158*t240
      t246 = CsThl*t240
      t247 = t246*sll
      t250 = sll*CsThl
      t251 = t250*CsThr
      t254 = t235*t240
      t255 = t254*sqq
      t258 = t41*sll
      t260 = t61*t32
      t267 = t14*t4
      t268 = t267*t11
      t270 = t67*El
      t273 = t4*El
      t274 = t2*t273
      t275 = t6*SnThl
      t276 = t7*SnThr
      t277 = t275*t276
      t278 = t10*t9
      t279 = t277*t278
      t282 = 48*t24+24*t30-8*t36-8*t91*t61*t94+4*t83*t85*t86*t229-t234-2
     #*t45*t236*sqq+2*t54*t241-4*t3*t241-4*t54*t247+4*t54*t251-2*t54*t25
     #5+t118*t258+2*t118*t260+4*t3*t247+8*t3*t255-96*t268-96*t270*t87-32
     #*t274*t279
      t283 = t91*El
      t284 = t283*t87
      t286 = 40*t42
      t287 = 12*t46
      t293 = 56*t74
      t294 = t190*t152
      t295 = 4*t294
      t296 = t80*t58
      t300 = t14*CsThr
      t307 = t118*t202
      t308 = t61*t76
      t309 = 4*t308
      t310 = t14*t235
      t313 = 8*t310*t240*sqq
      t314 = 8*t284+t286+t287+6*t55-8*t59-8*t62+t66-32*t69+64*t71-t293-t
     #295+16*t296-16*t80*t246+8*t300*t236+8*t57*t254-32*t68*t254-t307-t3
     #09-t313
      t319 = 4*t3*t41
      t320 = t32*t40
      t322 = 8*t76
      t323 = 24*t73
      t331 = t57*t40
      t333 = t61*t233
      t334 = 96*t333
      t335 = CsThr*t122
      t336 = t335*CsThl
      t337 = t57*t38
      t339 = 24*t57*t246+t319-t118*t320-t322+t323+8*t57+16*t89-64*t100-t
     #104+t109-6*t111+12*t113-24*t116-4*t198+t122-24*t331-t334-t336-8*t3
     #37
      t340 = t58*t3
      t341 = t7*t10
      t342 = t26*t341
      t351 = Eq*t32
      t352 = El*SnThl
      t357 = t1*sqq
      t362 = sqq*t2
      t369 = t58*t233
      t370 = 96*t369
      t371 = t3*sll
      t372 = 4*t371
      t375 = 32*t233*t235*t240
      t377 = 4*t131
      t385 = t14*CsThl
      t386 = t385*CsThr
      t392 = t9*t38
      t393 = t86*t392
      t396 = t83*sll
      t399 = -32*t340*t342+64*t15*t61*t94-4*t83*t174*t94+4*t351*t352*t17
     #*t58-8*t357*t26*t341*t58+16*t61*t362*t94-8*t83*t61*t88+t370+t372-t
     #125+t375-2*t129+t377+4*t134-24*t91*t40*t94+192*t58*t67*t94+96*t386
     #*t342-96*t61*t67*t94-8*t283*t393+4*t396*t94
      t402 = 1/sll
      t404 = 1/sqq
      t406 = t140*t138
      t407 = t218*t145
      A1e0c3p1t1n = -48*(t282+t314+t339+t399)*t402*t404*t406*t407
      t414 = 1/(-sll+Mzsq)
      t416 = sll*Mzsq
      t417 = abs(DCMPLX(Mzsq))
      t418 = t417**2
      t420 = 1/(t416-t418)
      t422 = t140*t414*t420*t145
      A1e0c4p1t1n = -24*Mzsq*t137*t138*t422
      t427 = sqq*Mzsq
      t429 = 1/(t427-t418)
      t432 = 1/(-sqq+Mzsq)
      t434 = t140*t429*t432*t218
      A1e0c5p1t1n = 24*Mzsq*t213*t138*t434
      t436 = t58*Mzsq
      t437 = conjg(Mzsq)
      t438 = t32*t437
      t439 = t438*t1
      t441 = 24*t436*t439
      t442 = t39*t40
      t443 = t437*sqq
      t444 = t443*sll
      t445 = t442*t444
      t447 = Mzsq*t437
      t448 = t447*sqq
      t449 = t442*t448
      t451 = t3*CsThl
      t452 = CsThr*sqq
      t453 = t452*t141
      t456 = CsThr*sll
      t457 = t456*t32
      t458 = t451*t457
      t460 = t67*CsThl
      t461 = t460*Er
      t462 = t456*sqq
      t463 = t461*t462
      t465 = CsThr*Mzsq
      t466 = t465*t141
      t469 = t68*t235
      t470 = t240*Mzsq
      t471 = t470*sqq
      t474 = t240*sll
      t475 = t474*Mzsq
      t478 = t474*sqq
      t481 = t67*t38
      t482 = t481*Er
      t483 = t482*t258
      t485 = Mzsq*t67
      t486 = t485*Er
      t487 = sll*t38
      t488 = t487*t40
      t489 = t486*t488
      t491 = t486*t251
      t493 = t441-40*t445+80*t449-8*t451*t453+24*t458-64*t463+4*t451*t46
     #6+32*t469*t471+32*t469*t475-64*t469*t478+128*t483-64*t489+32*t491
      t494 = t118*t40
      t495 = sqq*t141
      t496 = t495*Mzsq
      t497 = t494*t496
      t498 = t118*t38
      t499 = t40*t122
      t500 = t499*Mzsq
      t502 = 2*t498*t500
      t503 = t174*t122
      t505 = 4*t498*t503
      t506 = t61*Mzsq
      t507 = t165*t32
      t509 = 4*t506*t507
      t510 = t254*t14
      t511 = t427*sll
      t514 = t57*sll
      t515 = t514*t110
      t517 = t57*Mzsq
      t518 = t517*t488
      t520 = t54*Mzsq
      t521 = t520*t119
      t523 = t54*CsThl
      t524 = t523*t457
      t526 = t54*t235
      t527 = CsThr*t32
      t528 = t527*Mzsq
      t533 = t57*t235
      t534 = t465*sqq
      t537 = t465*sll
      t541 = 12*t436*t507
      t542 = -t497-t502+t505+t509+8*t510*t511-16*t515+8*t518+6*t521-12*t
     #524+2*t526*t528-4*t526*t457-8*t533*t534-8*t533*t537-t541
      t544 = t386*t511
      t546 = t442*t511
      t548 = t385*t49
      t549 = t548*t537
      t551 = Mzsq*t9
      t553 = t86*t202
      t554 = t21*t551*t553
      t557 = t91*t38*El
      t558 = t9*t437
      t559 = t558*sqq
      t560 = t86*t559
      t563 = t551*sqq
      t564 = t86*t563
      t567 = t86*sll
      t568 = t551*t567
      t571 = t558*t567
      t574 = t91*t202
      t576 = t352*t17*t38
      t579 = El*Mzsq
      t580 = t83*t579
      t583 = t58*El
      t584 = t83*t583
      t585 = t9*SnThl
      t586 = SnThr*sll
      t588 = t585*t586*t32
      t591 = t15*t583
      t593 = 64*t591*t560
      t594 = t58*t4
      t595 = t3*t594
      t596 = t10*t437
      t597 = t596*sqq
      t598 = t8*t597
      t601 = 56*t544-40*t546+8*t549-48*t554+8*t557*t560+8*t557*t564+8*t5
     #57*t568+8*t557*t571-16*t574*t576+8*t580*t571+32*t584*t588+t593+32*
     #t595*t598
      t602 = t10*Mzsq
      t603 = t602*sqq
      t607 = t83*t119
      t608 = t17*t437
      t609 = t352*t608
      t612 = t447*t11
      t615 = t447*t87
      t617 = 128*t591*t615
      t618 = t15*t251
      t620 = 64*t618*t609
      t622 = 128*t618*t88
      t624 = t21*t558*t553
      t626 = El*t437
      t628 = t586*t351
      t630 = 8*t626*t585*t628
      t631 = t4*Mzsq
      t634 = t631*t596*t8*t357
      t639 = 16*t579*t558*t86*t351
      t640 = t21*t437
      t642 = t586*sqq
      t643 = t585*t642
      t646 = t21*Mzsq
      t647 = t61*t646
      t652 = t40*El
      t653 = t91*t652
      t656 = 32*t595*t8*t603-16*t607*t609-64*t595*t612-t617+t620-t622-48
     #*t624+t630+48*t634-t639-16*t61*t640*t643-16*t647*t643+32*t647*t560
     #+24*t653*t568
      t659 = t267*t437
      t660 = t10*t6
      t661 = t7*sll
      t662 = t660*t661
      t665 = Mzsq*t32
      t666 = t665*El
      t669 = t666*t86*t9*t2
      t671 = t83*t626
      t672 = SnThr*t141
      t673 = t585*t672
      t678 = t58*t270
      t681 = t58*t267
      t695 = t86*t32
      t696 = t558*t695
      t697 = t97*t696
      t699 = t58*t427
      t704 = 96*t659*t662-48*t669-4*t671*t673-4*t580*t673+384*t678*t615+
     #192*t681*t612-192*t58*t485*t88+384*t678*t643+192*t681*t660*t661*sq
     #q-192*t678*t560+64*t697-96*t699*t268-96*t681*t598
      t705 = t8*sll
      t706 = t596*t705
      t709 = t58*t665
      t711 = t352*t17*t2
      t712 = t709*t711
      t714 = t97*t588
      t723 = t58*t646
      t724 = t723*t643
      t726 = t723*t560
      t729 = t58*t640*t643
      t734 = t61*t270
      t743 = -96*t681*t706+64*t712-128*t714-192*t678*t568-192*t678*t571-
     #96*t681*t602*t705+64*t724-128*t726+64*t729+96*t61*t485*t88+96*t734
     #*t568+96*t734*t571-192*t734*t643-192*t734*t615
      t747 = t233*Mzsq
      t748 = t747*sll
      t751 = t233*sll
      t752 = t751*sqq
      t755 = t747*sqq
      t762 = t665*t14
      t763 = t58*t762
      t765 = Mzsq*t122
      t766 = t765*t1
      t768 = 12*t58*t766
      t769 = sll*t122
      t770 = t769*t1
      t772 = 24*t58*t770
      t773 = t177*t32
      t774 = t58*t773
      t778 = t61*t773
      t780 = t61*t762
      t782 = 96*t734*t560-32*t254*t748+64*t254*t752-32*t254*t755-96*t58*
     #t755-96*t58*t748+56*t763-t768+t772-112*t774+192*t58*t752+80*t778-4
     #0*t780
      t789 = t416*t336
      t790 = t3*t511
      t791 = 8*t790
      t792 = t118*t500
      t794 = 2*t118*t503
      t796 = 8*t61*t770
      t797 = Mzsq*t40
      t798 = t797*t32
      t800 = 4*t3*t798
      t801 = t416*t32
      t803 = 2*t118*t801
      t804 = t118*t496
      t805 = t41*t141
      t807 = 8*t3*t805
      t808 = t797*t141
      t810 = 4*t3*t808
      t811 = t174*t32
      t813 = 8*t3*t811
      t814 = 96*t61*t748+96*t61*t755-192*t61*t752+t789-t791+t792-t794-t7
     #96-t800+t803+t804-t807+t810+t813
      t818 = t320*t141
      t820 = 2*t118*t818
      t821 = t54*t511
      t824 = t54*t184*t38
      t825 = 4*t824
      t826 = t797*sqq
      t827 = t57*t826
      t829 = t797*sll
      t830 = t57*t829
      t832 = t57*t258
      t834 = t38*Mzsq
      t835 = t834*sll
      t836 = t57*t835
      t838 = t834*sqq
      t839 = t57*t838
      t842 = 6*t54*t798
      t844 = 12*t54*t811
      t845 = t38*t32
      t846 = t845*Mzsq
      t847 = t54*t846
      t848 = 2*t847
      t849 = t487*sqq
      t850 = t57*t849
      t852 = t233*t437
      t853 = t852*sll
      t856 = t820+2*t821+t825+24*t827+24*t830-48*t832+8*t836+8*t839-t842
     #+t844-t848-16*t850-32*t254*t853
      t857 = t747*t437
      t860 = t852*sqq
      t863 = t437*t122
      t864 = t863*t1
      t866 = 12*t58*t864
      t867 = t14*t437
      t868 = t867*t32
      t869 = t58*t868
      t883 = t61*t868
      t886 = 2*t447*t336
      t887 = t437*sll
      t888 = t887*t336
      t889 = t3*t444
      t890 = 8*t889
      t891 = 64*t254*t857-32*t254*t860-t866+56*t869+192*t58*t857-96*t58*
     #t860-96*t58*t853-192*t61*t857+96*t61*t853+96*t61*t860-40*t883-t886
     #+t888-t890
      t893 = t40*t437
      t894 = t893*t122
      t895 = t118*t894
      t897 = 4*t61*t864
      t898 = t437*t141
      t899 = t898*sqq
      t900 = t118*t899
      t901 = t887*t32
      t903 = 2*t118*t901
      t904 = t447*t32
      t905 = t118*t904
      t906 = 2*t905
      t907 = t447*sll
      t909 = 8*t3*t907
      t910 = t3*t448
      t911 = 8*t910
      t912 = t320*t437
      t914 = 4*t3*t912
      t915 = t893*t141
      t917 = 4*t3*t915
      t918 = t797*t437
      t919 = t57*t918
      t921 = t893*sll
      t922 = t57*t921
      t924 = t893*sqq
      t925 = t57*t924
      t927 = t38*t437
      t928 = t927*sqq
      t929 = t57*t928
      t931 = t895+t897+t900+t903-t906+t909+t911-t914+t917-48*t919+24*t92
     #2+24*t925+8*t929
      t932 = t927*sll
      t933 = t57*t932
      t935 = t834*t437
      t936 = t57*t935
      t938 = t54*t444
      t940 = t54*t448
      t943 = 6*t54*t912
      t944 = t845*t437
      t945 = t54*t944
      t946 = 2*t945
      t950 = 4*t61*t766
      t957 = t341*sqq
      t961 = t83*t665
      t962 = t961*t94
      t965 = t26*t341*Mzsq
      t968 = 8*t933-16*t936+2*t938-4*t940-t943-t946+8*t254*t868+t950+8*t
     #254*t762-16*t254*t773-16*t557*t615-64*t198*t26*t957-4*t962+32*t198
     #*t965
      t971 = t341*t437
      t972 = t26*t971
      t975 = t83*t436
      t976 = t33*t87
      t979 = t447*t9
      t980 = t86*sqq
      t981 = t979*t980
      t984 = t558*SnThl
      t985 = t984*t642
      t989 = t202*El*t87
      t992 = t15*t506
      t995 = t61*El
      t996 = t15*t995
      t999 = t15*t488
      t1008 = t83*t652
      t1014 = 32*t198*t972-16*t975*t976+32*t584*t981-16*t584*t985-16*t97
     #5*t989-64*t992*t149+128*t996*t615-64*t999*t609+128*t999*t88-64*t99
     #6*t560-64*t992*t88+8*t1008*t588-4*t83*t320*t609
      t1019 = t86*t141
      t1020 = t558*t1019
      t1023 = t551*t1019
      t1027 = t83*t652*Mzsq
      t1032 = t83*t995
      t1035 = t83*t260
      t1038 = t83*t506
      t1043 = t631*t660
      t1044 = t661*t357
      t1045 = t1043*t1044
      t1049 = 8*t579*t585*t628
      t1050 = t4*t437
      t1051 = t1050*t660
      t1052 = t1051*t1044
      t1054 = t54*t40
      t1055 = t1054*t511
      t1057 = t54*t38
      t1058 = t1057*t511
      t1060 = -8*t1008*t585*t672*sqq+4*t1008*t1020+4*t1008*t1023-8*t1027
     #*t571+8*t1027*t560-16*t1032*t588+8*t1035*t609+8*t1038*t976-16*t103
     #2*t981-24*t1045+t1049-24*t1052-6*t1055-2*t1058
      t1062 = t235*CsThr
      t1066 = t3*t38
      t1071 = t195*t50
      t1074 = t80*CsThl
      t1079 = t80*t240
      t1080 = CsThl*Mzsq
      t1081 = t1080*sqq
      t1084 = t1080*sll
      t1089 = t548*t534
      t1091 = t57*t240
      t1100 = 16*t514*t1062*sqq-48*t1066*t811+24*t1066*t798+32*t1071*t15
     #9-16*t1074*t537-16*t1074*t534+16*t1079*t1081+16*t1079*t1084-32*t10
     #79*t159+8*t1089-24*t1091*t1081-24*t1091*t1084+48*t1091*t159-8*t533
     #*t475
      t1105 = t337*t826
      t1107 = t337*t258
      t1109 = t1057*t811
      t1111 = t1057*t798
      t1117 = t240*t32
      t1118 = t1117*Mzsq
      t1121 = t240*t141
      t1122 = t1121*sqq
      t1125 = t470*t141
      t1128 = t1117*sll
      t1135 = -8*t533*t471+16*t533*t478+8*t1105-16*t1107+4*t1109-2*t1111
     #+8*t523*t453-4*t523*t466-2*t523*t1118-8*t523*t1122+4*t523*t1125+4*
     #t523*t1128-4*t526*t1128+2*t526*t1118
      t1140 = t486*t115
      t1142 = Mzsq*t2
      t1143 = t1142*Er
      t1144 = t1143*t119
      t1146 = t246*t32
      t1149 = t486*t110
      t1157 = t3*t235
      t1162 = t1054*t444
      t1164 = t1057*t444
      t1167 = 2*t494*t904
      t1170 = -64*t1140-12*t1144+4*t1143*t1146+32*t1149-4*t451*t1125+8*t
     #451*t1122-8*t451*t1128+16*t1157*t1128-8*t1157*t1118-6*t1162-2*t116
     #4-t1167+24*t1066*t912
      t1171 = CsThr*t437
      t1172 = t1171*sqq
      t1175 = t1171*sll
      t1178 = t465*t437
      t1181 = t1080*t437
      t1184 = CsThl*t437
      t1185 = t1184*sll
      t1188 = t1184*sqq
      t1191 = t548*t1172
      t1193 = t548*t1175
      t1195 = t548*t1178
      t1203 = t240*t437
      t1204 = t1203*sqq
      t1207 = t1203*sll
      t1210 = -16*t1074*t1172-16*t1074*t1175+32*t1074*t1178-32*t1079*t11
     #81+16*t1079*t1185+16*t1079*t1188+8*t1191+8*t1193-16*t1195+48*t1091
     #*t1181-24*t1091*t1185-24*t1091*t1188-8*t533*t1204-8*t533*t1207
      t1212 = t470*t437
      t1215 = t337*t918
      t1217 = t337*t921
      t1219 = t337*t924
      t1221 = t1057*t912
      t1224 = 12*t1054*t448
      t1225 = t1057*t448
      t1226 = 4*t1225
      t1227 = t898*CsThr
      t1230 = t1203*t141
      t1233 = t1117*t437
      t1238 = t3*t40
      t1240 = 8*t1238*t448
      t1242 = 8*t1238*t907
      t1243 = 16*t533*t1212-16*t1215+8*t1217+8*t1219-2*t1221+t1224+t1226
     #-4*t523*t1227+4*t523*t1230-2*t523*t1233+2*t526*t1233+t1240-t1242
      t1244 = t118*Mzsq
      t1246 = 2*t1244*t444
      t1253 = t461*t1172
      t1255 = t527*t437
      t1256 = t451*t1255
      t1260 = t461*t1175
      t1262 = t461*t1178
      t1270 = t482*t924
      t1272 = t482*t918
      t1274 = -t1246-4*t451*t1230+4*t451*t1233-8*t1157*t1233+32*t1253-12
     #*t1256+4*t451*t1227+32*t1260-64*t1262+32*t469*t1204+32*t469*t1207-
     #64*t469*t1212-64*t1270+128*t1272
      t1277 = t482*t921
      t1279 = t494*t899
      t1281 = 2*t498*t894
      t1282 = t61*t437
      t1284 = 4*t1282*t507
      t1286 = 8*t506*t439
      t1291 = t523*t1255
      t1301 = t58*t437
      t1303 = 12*t1301*t507
      t1304 = -64*t1277-t1279-t1281+t1284-t1286-16*t510*t448+8*t510*t444
     #+6*t1291+2*t526*t1255+16*t533*t1178-8*t533*t1175-8*t533*t1172-t130
     #3
      t1305 = t386*t444
      t1307 = t386*t448
      t1309 = t665*t54
      t1312 = 8*t416*t76
      t1313 = Mzsq*t14
      t1314 = t1313*t202
      t1315 = 24*t1314
      t1317 = 8*t887*t76
      t1318 = t867*t202
      t1319 = 24*t1318
      t1320 = t1313*t443
      t1321 = 48*t1320
      t1323 = 16*t447*t76
      t1324 = t32**2
      t1325 = sll*t1324
      t1327 = 2*t1325*t58
      t1328 = Mzsq*t1324
      t1329 = t1328*t58
      t1330 = t437*t1324
      t1331 = t1330*t58
      t1332 = t57*t443
      t1334 = t57*t887
      t1336 = 56*t1305-112*t1307+2*t1309+t1312-t1315+t1317-t1319+t1321-t
     #1323-t1327+t1329+t1331-8*t1332-8*t1334
      t1338 = t57*t447
      t1340 = t54*t438
      t1342 = t3*t438
      t1343 = 4*t1342
      t1345 = 4*t3*t898
      t1346 = t118*t863
      t1347 = t57*t202
      t1349 = t57*t416
      t1351 = t57*t427
      t1353 = t54*t184
      t1355 = t118*t765
      t1356 = t32*t141
      t1358 = 2*t118*t1356
      t1359 = t118*t769
      t1360 = 2*t1359
      t1362 = 8*t3*t495
      t1363 = 16*t1338+2*t1340-t1343-t1345+t1346+16*t1347-8*t1349-8*t135
     #1-4*t1353+t1355-t1358-t1360+t1362
      t1364 = t3*t184
      t1365 = 8*t1364
      t1366 = t3*t665
      t1367 = 4*t1366
      t1368 = t141*Mzsq
      t1370 = 4*t3*t1368
      t1371 = 2*t1325
      t1372 = t91*t995
      t1385 = t83*t93
      t1386 = t86*t184
      t1387 = t1385*t1386
      t1389 = t83*t33
      t1390 = t86*t558
      t1391 = t1389*t1390
      t1393 = t1365-t1367-t1370+t1371-t1330-t1328+8*t1372*t568-16*t1372*
     #t615+24*t653*t564-48*t653*t643+24*t653*t560+24*t653*t571+8*t1387-4
     #*t1391
      t1400 = 64*t857
      t1401 = 16*t770
      t1402 = 32*t748
      t1403 = 32*t755
      t1404 = t416*t122
      t1405 = 24*t868
      t1406 = 32*t853
      t1408 = t21*t447*t99
      t1410 = t671*t643
      t1412 = t83*t416
      t1413 = t1412*t88
      t1415 = Eq*El
      t1416 = t1415*t558
      t1417 = t86*t336
      t1419 = 4*t1416*t1417
      t1420 = t1415*t551
      t1422 = 4*t1420*t1417
      t1423 = 8*t1385*t86*t495-t1400-t1401+t1402+t1403-t1404-t1405+t1406
     #+96*t1408-8*t1410-8*t1413-t1419-t1422
      t1427 = 8*t1415*t585*t586*t336
      t1428 = t58*t507
      t1431 = t7*t32
      t1432 = t1*CsThl
      t1433 = t1432*CsThr
      t1434 = t1431*t1433
      t1440 = t586*t119
      t1442 = 4*t1415*t984*t1440
      t1446 = 4*t1415*t551*SnThl*t1440
      t1450 = 8*t1415*t979*t86*t119
      t1457 = t7*sqq
      t1461 = t61*t21
      t1466 = t61*t665
      t1471 = t1427-16*t342*t1428+8*t1051*t1434+8*t1043*t1434-t1442-t144
     #6+t1450+8*t972*t200+8*t965*t200-16*t631*t596*t6*t1457*t1433+32*t14
     #61*t588-16*t1461*t696-16*t1466*t711-48*t653*t615
      t1473 = t447*t122
      t1474 = 2*t1473
      t1475 = 8*t864
      t1476 = 8*t766
      t1477 = 64*t752
      t1478 = 24*t762
      t1479 = 32*t860
      t1480 = 48*t773
      t1481 = t887*t122
      t1487 = t352*t17*t40
      t1492 = 64*t15*t436*t88
      t1493 = t580*t560
      t1495 = t1474+t1475+t1476-t1477-t1478+t1479+t1480-t1481+8*t1032*t9
     #85+8*t1038*t989-4*t961*t1487+t1492+8*t1493
      t1496 = t15*t1084
      t1497 = CsThr*El
      t1498 = t1497*t87
      t1500 = 64*t1496*t1498
      t1503 = t91*t352
      t1504 = t17*t427
      t1507 = t91*t579
      t1508 = t585*t586
      t1516 = t17*t443
      t1527 = t270*Mzsq
      t1530 = t267*Mzsq
      t1531 = t596*t8
      t1537 = t1500+16*t574*t94-8*t1503*t1504-8*t1507*t1508-8*t91*t626*t
     #1508+16*t1507*t1390-8*t1503*t1516+8*t1372*t571+8*t1372*t564+8*t137
     #2*t560-16*t1372*t643-192*t1527*t1390-192*t1530*t1531+96*t485*sqq*t
     #94
      t1547 = t270*t437
      t1550 = SnThr*t32
      t1551 = t585*t1550
      t1552 = t640*t1551
      t1555 = t14*t6
      t1556 = t1555*t341
      t1559 = t660*t1457
      t1563 = t21*t9*t1386
      t1571 = t274*t437
      t1572 = t278*t275
      t1574 = t1572*t276*sll
      t1577 = t274*Mzsq
      t1583 = t4*t10*t6*t661*t76
      t1585 = -192*t270*t9*t553-192*t267*t10*t8*t202+96*t1547*t99-48*t15
     #52+96*t427*t4*t1556+96*t659*t1559+96*t1563+96*t1527*t1508+96*t1547
     #*t1508+96*t1530*t662+32*t1571*t1574+32*t1577*t1574+48*t1583
      t1587 = t1572*t276*sqq
      t1595 = t122*Eq
      t1598 = 16*t93*SnThl*t586*t1595
      t1600 = t86*t1595
      t1602 = 8*t579*t9*t1600
      t1606 = t797*t202
      t1609 = t470*t202
      t1614 = t1057*t1606
      t1616 = t250*t452
      t1617 = t520*t1616
      t1619 = t1143*t1616
      t1625 = 2*t498*t797*t184
      t1627 = t8*t76
      t1628 = t631*t10*t1627
      t1630 = 32*t1577*t1587-64*t1577*t437*t278*t277-t1598+t1602+2*t526*
     #t465*t202+24*t1066*t1606+2*t523*t1609+2*t526*t1609-2*t1614+2*t1617
     #-8*t1619-8*t1157*t1609-t1625-24*t1628
      t1632 = t626*t9
      t1634 = 8*t1632*t1600
      t1636 = t1050*t10*t1627
      t1644 = t1171*t202
      t1647 = t465*t443
      t1650 = t523*t1647
      t1652 = t893*t202
      t1655 = t797*t443
      t1658 = t447*t456
      t1661 = t470*t443
      t1664 = t470*t887
      t1667 = t1203*t202
      t1672 = t1634-24*t1636-64*t274*t278*t277*t202+32*t1571*t1587+2*t52
     #6*t1644-4*t526*t1647-12*t1650+24*t1066*t1652-48*t1066*t1655+8*t523
     #*t1658+4*t523*t1661-8*t523*t1664+2*t523*t1667+2*t526*t1667
      t1675 = t1057*t1652
      t1677 = t1057*t1655
      t1679 = t523*t1644
      t1689 = t447*t202
      t1691 = 2*t494*t1689
      t1694 = 2*t498*t893*t184
      t1695 = t797*t438
      t1697 = 4*t498*t1695
      t1700 = t451*t1647
      t1702 = t451*t1644
      t1704 = -4*t526*t1661-2*t1675+4*t1677+2*t1679-8*t451*t1661+8*t451*
     #t1664-8*t1157*t1667+16*t1157*t1661+t1691-t1694+t1697-8*t451*t1658+
     #24*t1700-8*t1702
      t1712 = 1/(-sqq+t437)
      t1714 = 1/(-sll+t437)
      t1715 = t1712*t1714
      t1720 = t432*t414*t138*t145*t218
      A1e0c6p1t1n = -24*(t1495+t1537+t931+t1274+t1243+t1704+t1210+t1672+
     #t1170+t1630+t1135+t1100+t542+t1060+t493+t1014+t1471+t968+t891+t656
     #+t1423+t601+t1393+t1363+t856+t1585+t704+t1336+t814+t782+t743+t1304
     #)*t140*t1715*t1720
      t1726 = 16*t1066*t258
      t1728 = 4*t1057*t258
      t1729 = t205*t110
      t1731 = t371*t110
      t1732 = 16*t1731
      t1735 = t523*t1172
      t1736 = 8*t1735
      t1739 = t118*CsThl
      t1741 = 2*t1739*t1255
      t1742 = t451*t1172
      t1744 = t520*t110
      t1745 = 8*t1744
      t1746 = t1143*t110
      t1750 = -2*t1057*t924-t1726+t1728+16*t1729-t1732-2*t1057*t826-t173
     #6+8*t1066*t924+t1741+8*t1742-t1745+8*t1746+8*t1143*t115
      t1752 = 2*t1244*t119
      t1754 = 4*t1739*t457
      t1755 = Mzsq*t1
      t1756 = Er*sqq
      t1758 = t1755*t1756*t94
      t1759 = 16*t1758
      t1761 = t83*t202*t94
      t1763 = t8*t602
      t1765 = 32*t5*t1763
      t1767 = 32*t5*t1531
      t1770 = 64*t15*t437*t94
      t1773 = t15*sll
      t1778 = 64*t1313*Er*t94
      t1781 = 8*t1415*t9*t1386
      t1783 = t17*Mzsq*Eq
      t1785 = 4*t34*t1783
      t1788 = t1752-t1754-t1759+32*t1761+t1765+t1767+t1770-64*t371*t342-
     #128*t1773*t94+t1778-t1781+t1785-16*t22*t1516
      t1799 = t17*t437*Eq
      t1801 = 4*t34*t1799
      t1803 = 32*t162*t99
      t1804 = sll*t4
      t1805 = t1804*t6
      t1808 = 16*t1805*t341*t357
      t1814 = 64*t1313*t196*t1498
      t1818 = 32*t92*t609
      t1824 = t49*CsThl
      t1827 = 32*t1142*t1824*t1498
      t1829 = 32*t153*t88
      t1830 = -16*t22*t1504-8*t27*t596*t357-8*t27*t602*t357+t1801+t1803+
     #t1808+16*t1080*t195*t88-t1814-8*t975*t88+t1818-64*t91*t456*t352*t1
     #7*CsThl+t1827-t1829
      t1836 = t83*t251
      t1838 = 16*t1836*t88
      t1842 = 64*t102*t609
      t1843 = t171*t357
      t1845 = t54*t258
      t1846 = 4*t1845
      t1847 = t54*t849
      t1848 = 4*t1847
      t1856 = 64*t57*t251
      t1857 = t1313*t115
      t1859 = 16*t97*t560+128*t15*t250*t1498+t1838-8*t584*t560-t1842+20*
     #t1843+t1846+t1848-64*t68*t488+128*t68*t251+64*t57*t488-t1856-8*t18
     #57
      t1862 = t665*t1433
      t1863 = 4*t1862
      t1864 = t186*sqq
      t1869 = 32*t1313*t1824*CsThr
      t1872 = 32*t1313*t106*t40
      t1873 = t1755*t107
      t1875 = t196*CsThr
      t1877 = 64*t485*t1875
      t1878 = Er*t38
      t1879 = t1878*t40
      t1881 = 32*t485*t1879
      t1882 = t300*sqq
      t1883 = t1080*t1882
      t1885 = 8*t1428
      t1886 = t177*sqq
      t1888 = 32*t58*t1886
      t1890 = 16*t61*t1886
      t1891 = t867*sqq
      t1892 = t61*t1891
      t1894 = t58*t439
      t1895 = 4*t1894
      t1896 = -t1863-2*t1755*t1864+t1869-t1872-2*t1873-t1877+t1881+16*t1
     #883+t1885-t1888+t1890-8*t1892-t1895
      t1897 = t58*t1891
      t1900 = 32*t68*t1282
      t1902 = 64*t68*t1301
      t1903 = t54*t928
      t1908 = 32*t57*t1282
      t1910 = 32*t57*t1301
      t1911 = t80*t437
      t1913 = 8*t61*t1911
      t1914 = t80*sll
      t1917 = t1142*t50
      t1919 = 8*t61*t1917
      t1920 = t80*t927
      t1921 = 8*t1920
      t1922 = t665*t118
      t1923 = 6*t1922
      t1924 = Mzsq*t49
      t1925 = t1924*t357
      t1927 = 16*t1897+t1900-t1902-2*t1903-2*t54*t924-t1908+t1910+t1913-
     #16*t61*t1914+t1919+t1921+t1923-10*t1925
      t1929 = t437*t49
      t1930 = t1929*t357
      t1932 = t1142*t51
      t1933 = 8*t1932
      t1934 = t1142*t1756
      t1935 = 16*t1934
      t1936 = t50*t437
      t1937 = t48*t1936
      t1938 = 8*t1937
      t1939 = t3*t202
      t1940 = 32*t1939
      t1942 = 12*t118*t184
      t1943 = t797*t80
      t1944 = 8*t1943
      t1945 = t3*t443
      t1946 = 16*t1945
      t1947 = t118*t438
      t1948 = 6*t1947
      t1949 = t80*t487
      t1951 = t50*sll
      t1954 = 2*t769
      t1955 = -10*t1930+t1933-t1935+t1938+t1940-t1942+t1944-t1946+t1948-
     #16*t1949-16*t48*t1951-t765+t1954
      t1956 = 8*t507
      t1957 = t1313*sqq
      t1958 = 8*t1957
      t1959 = 4*t439
      t1960 = 16*t1886
      t1961 = 8*t1891
      t1962 = t671*t99
      t1963 = 16*t1962
      t1964 = t665*t1
      t1965 = 4*t1964
      t1966 = 8*t1917
      t1968 = 32*t68*t437
      t1969 = 32*t486
      t1971 = t68*sll
      t1972 = 64*t1971
      t1973 = 8*t1911
      t1974 = -t863-t1956-t1958+t1959+t1960-t1961-t1963+t1965+t1966+t196
     #8+t1969-16*t1914-t1972+t1973
      t1980 = t414*t1714
      t1982 = t1980*t402*t145
      A1e0c7p1t1n = -24*(t1750+t1788+t1830+t1859+t1896+t1927+t1955+t1974
     #)*t138*t140*t1982
      t1987 = t191*t40
      t1988 = Er*sll
      t1989 = t1988*Mzsq
      t1992 = t2*CsThl
      t1993 = t1992*CsThr
      t1996 = t165*CsThl
      t1999 = t61*t1
      t2000 = t1924*sll
      t2007 = t1726-t1728-t1732-t1736+16*t1742-t1745+16*t1746+t1754-8*t1
     #987*t1989+8*t1993*t1989+4*t1996*t534+2*t1999*t2000-16*t177*t1301-2
     #*t170*t2000
      t2008 = t1929*sll
      t2011 = Eq*CsThl
      t2012 = t2011*CsThr
      t2018 = t86*t551
      t2024 = t437*t1
      t2028 = 2*t1999*t2008-t1759+t1765+t1767+t1770+t1778+t1781-t1803-t1
     #808-2*t2012*t1756*t416+16*t162*t1390+16*t162*t2018+8*t1805*t341*t1
     #755+8*t1805*t341*t2024
      t2037 = 4*t1873
      t2038 = 32*t83*t32*t94-128*t15*El*t99-t1814+t1818+t1827+t1829-t183
     #8-t1842+4*t1843-t1846+t1848+t1869-t1872+t2037
      t2039 = 4*t1903
      t2042 = t105*t2000
      t2044 = t190*sqq
      t2047 = -t1877+t1881-t1885+t1888-t1890+t1900-t1902+t2039-t1908+t19
     #10+t1913+t1919-2*t170*t2008+2*t2042-64*t481*t2044
      t2050 = CsThr*Er
      t2051 = t2050*t32
      t2052 = t1992*t2051
      t2054 = t2050*sqq
      t2061 = t50*sqq
      t2067 = t49*t32
      t2068 = t2067*CsThl
      t2071 = t49*sqq
      t2072 = t2071*CsThl
      t2075 = t105*t2008
      t2087 = t579*t87
      t2092 = -32*t2052+128*t460*t2054+8*t177*t506-16*t177*t436-16*t48*t
     #2061*t38+64*t39*t1864+16*t45*t2068-64*t300*t2072+2*t2075+8*t177*t1
     #282-4*t1420*t553-64*t195*t49*El*t86*t98*CsThl-16*t153*t2087-16*t15
     #3*t609
      t2099 = 2*t1922
      t2102 = 2*t1947
      t2103 = t416*t54
      t2105 = 128*t102*t88+8*t584*t571+8*t584*t568-t1921+t2099+4*t1925+4
     #*t1930-t1933-t1935+t1938+t1944-t1946+t2102-2*t2103
      t2107 = t887*t54
      t2110 = 4*t2024*t202
      t2112 = 4*t1755*t202
      t2113 = t48*t2061
      t2121 = Er*t437
      t2127 = t3*t32
      t2128 = 32*t2127
      t2130 = 8*t177*Mzsq
      t2131 = -2*t2107-t2110-t2112-16*t2113+16*t80*t64-8*t54*t845-64*t5*
     #t1559-2*t2012*t2121*t202-t1954-4*t1416*t553+t1956+t2128+t901+t2130
      t2132 = t68*sqq
      t2133 = 64*t2132
      t2134 = t80*sqq
      t2137 = 8*t177*t437
      t2138 = t2067*t1
      t2140 = t118*t122
      t2141 = 4*t2140
      t2142 = t2121*sll
      t2149 = -t1960-t1963-t2133+16*t2134+t2137-t1966+t1968+t1969+t801-8
     #*t2138-t1973-t2141-8*t1987*t2142+8*t1993*t2142+4*t1996*t1172
      t2155 = t432*t1712
      t2156 = t404*t218
      A1e0c8p1t1n = 24*(t2007+t2028+t2038+t2047+t2092+t2105+t2131+t2149)
     #*t138*t140*t2155*t2156
      t2159 = t57*t487
      t2160 = t2159*t1655
      t2162 = t76*t26
      t2163 = Mzsq**2
      t2164 = sll*t2163
      t2165 = t341*t2164
      t2166 = t2162*t2165
      t2169 = t58*t267*t6
      t2176 = t898*Mzsq
      t2180 = t300*t4
      t2181 = t250*t2180
      t2182 = t602*t32
      t2186 = t10*t2163
      t2187 = t2186*sqq
      t2191 = t2186*t437
      t2192 = t8*t2191
      t2195 = t437**2
      t2196 = t10*t2195
      t2197 = t2196*sqq
      t2201 = t2196*Mzsq
      t2202 = t27*t2201
      t2205 = t352*SnThr
      t2207 = t2205*t551*t32
      t2210 = t9*t2163
      t2211 = t2210*sqq
      t2212 = t2205*t2211
      t2218 = t2210*t437
      t2219 = t2205*t2218
      t2223 = t2205*t558*t32
      t2226 = 32*t2160-48*t2166+288*t2169*t341*t496+288*t2169*t341*t899-
     #192*t2169*t341*t2176+288*t2181*t8*t2182-192*t2181*t8*t2187+96*t218
     #1*t2192+64*t198*t27*t2197-32*t198*t2202+192*t999*t2207-128*t999*t2
     #212-256*t1773*t995*t981+64*t999*t2219+192*t999*t2223
      t2227 = t9*t2195
      t2229 = t2205*t2227*sqq
      t2232 = t2227*Mzsq
      t2233 = t2205*t2232
      t2237 = t2205*t551*t122
      t2241 = t2205*t2210*t32
      t2244 = t396*t583
      t2245 = t87*t904
      t2248 = t141*CsThr
      t2250 = t158*t447
      t2251 = t54*t2248*t2250
      t2253 = t233*sqq
      t2255 = 96*t2253*t1368
      t2256 = t67*t32
      t2257 = t2256*t995
      t2264 = t67*sqq
      t2265 = t2264*t995
      t2269 = t86*t2227*sll
      t2272 = t86*t2232
      t2275 = t2264*t583
      t2278 = t87*t907
      t2283 = -128*t999*t2229+64*t999*t2233+48*t1836*t2237-32*t1836*t224
     #1-64*t2244*t2245-4*t2251-t2255-288*t2257*t568-288*t2257*t571+192*t
     #2257*t615-288*t2265*t1020+192*t2265*t2269-96*t2265*t2272+576*t2275
     #*t1020-768*t2275*t2278-384*t2275*t2269
      t2290 = t86*t2210*sll
      t2293 = t86*t2218
      t2296 = t141*El
      t2298 = t91*t2296*SnThl
      t2305 = t2*t32
      t2306 = t2305*t352
      t2307 = t17*t907
      t2308 = t2306*t2307
      t2310 = t14*t141
      t2311 = t2310*t1879
      t2312 = t2205*t979
      t2319 = t320*El
      t2320 = t83*t2319
      t2323 = t15*t110
      t2325 = 64*t2323*t2219
      t2327 = 64*t2323*t2233
      t2328 = t1389*t2272
      t2331 = t83*t765*t609
      t2333 = t665*CsThl
      t2335 = t1497*SnThl
      t2338 = 128*t15*t2333*t2335*t608
      t2339 = 192*t2275*t2272+576*t2275*t1023-384*t2275*t2290+192*t2275*
     #t2293-24*t2298*t17*t928-24*t2298*t17*t838-192*t2308-128*t2311*t231
     #2-8*t83*t499*El*t615+4*t2320*t2293-t2325-t2327+4*t2328-8*t2331+t23
     #38
      t2340 = t1*t141
      t2342 = t86*t122
      t2343 = t93*t2342
      t2346 = t190*t122
      t2350 = Er*Mzsq
      t2351 = t2350*t122
      t2355 = Er*El
      t2356 = t2340*t2355
      t2357 = t2227*t980
      t2360 = Eq*t141
      t2361 = Er*t2163
      t2362 = t2361*t32
      t2364 = 2*t2360*t2362
      t2365 = t141*sll
      t2366 = Eq*t2365
      t2367 = Er*t32
      t2368 = t2367*Mzsq
      t2370 = 3*t2366*t2368
      t2372 = 4*t2366*t2346
      t2373 = t1328*t437
      t2374 = t118*t2373
      t2375 = 2*t2374
      t2376 = t122*t2195
      t2377 = t2376*Mzsq
      t2378 = t118*t2377
      t2379 = t122*t2163
      t2380 = t2379*t437
      t2381 = t118*t2380
      t2382 = t3*t1473
      t2383 = 8*t2382
      t2384 = t32*t2195
      t2385 = t2384*Mzsq
      t2386 = t3*t2385
      t2387 = 4*t2386
      t2388 = t32*t2163
      t2389 = t2388*t437
      t2390 = t3*t2389
      t2391 = 4*t2390
      t2392 = t2367*t437
      t2394 = 3*t2366*t2392
      t2395 = Er*t2195
      t2396 = t2395*t32
      t2398 = 2*t2360*t2396
      t2399 = t54*t2385
      t2401 = -16*t2340*t190*t2343+12*t165*t2346*t609+12*t165*t2351*t148
     #7-8*t2356*t2357+t2364-t2370-t2372+t2375-t2378-t2381-t2383+t2387+t2
     #391-t2394+t2398-2*t2399
      t2404 = sqq*t2195
      t2405 = t2404*Mzsq
      t2406 = t57*t2405
      t2408 = sqq*t2163
      t2409 = t2408*t437
      t2410 = t57*t2409
      t2412 = t54*t2389
      t2414 = t54*t1473
      t2416 = t57*t904
      t2418 = t32*t14
      t2419 = t2418*t907
      t2420 = 96*t2419
      t2421 = t2195*t1324
      t2424 = 2*t456*t2421*CsThl
      t2425 = Er*t122
      t2426 = t2425*t437
      t2427 = t152*t2426
      t2428 = 12*t2427
      t2429 = t152*t2396
      t2430 = 8*t2429
      t2432 = 6*t2360*t2426
      t2433 = t2*t141
      t2434 = t2433*t2392
      t2435 = 24*t2434
      t2436 = sll*Eq
      t2437 = Er*t1324
      t2439 = t2436*t2437*t437
      t2440 = 3*t2439
      t2442 = t2436*t2425*t2195
      t2443 = 2*t2442
      t2444 = t2395*sqq
      t2446 = 8*t2433*t2444
      t2447 = t2395*Mzsq
      t2449 = 4*t2433*t2447
      t2450 = 8*t2406+8*t2410-2*t2412+4*t2414-16*t2416-t2420+t2424+t2428
     #-t2430-t2432+t2435-t2440+t2443-t2446+t2449
      t2451 = t2*t2365
      t2452 = t2121*sqq
      t2454 = 12*t2451*t2452
      t2455 = t2121*Mzsq
      t2457 = 8*t2451*t2455
      t2458 = t2361*t437
      t2460 = 4*t2433*t2458
      t2461 = t2310*t448
      t2462 = 48*t2461
      t2464 = 16*t2340*t904
      t2465 = t177*t2405
      t2466 = 24*t2465
      t2468 = 8*t165*t2385
      t2469 = t177*t2409
      t2470 = 24*t2469
      t2472 = 8*t165*t2389
      t2474 = t54*t898*t32
      t2476 = t57*t901
      t2478 = sll*t2195
      t2479 = t2478*Mzsq
      t2480 = t57*t2479
      t2482 = t2164*t437
      t2483 = t57*t2482
      t2486 = t54*t184*t2195
      t2488 = t2478*sqq
      t2489 = t57*t2488
      t2491 = t54*t1481
      t2493 = t2454-t2457+t2460-t2462+t2464+t2466-t2468+t2470-t2472-6*t2
     #474+24*t2476+8*t2480+8*t2483+4*t2486-16*t2489-6*t2491
      t2495 = t57*t899
      t2497 = t57*t2176
      t2499 = t233*t32
      t2500 = t58*t141
      t2503 = t14*t122
      t2504 = t2503*t2500
      t2506 = t61*t141
      t2507 = t2503*t2506
      t2509 = t254*t141
      t2512 = t1*t1324
      t2514 = 16*t2512*t2506
      t2517 = t2512*t141
      t2519 = 48*t58*t2517
      t2520 = t2499*t141
      t2523 = t1328*CsThl
      t2525 = 3*t2248*t2523
      t2526 = t1324*sqq
      t2527 = Mzsq*t2526
      t2528 = t2527*CsThl
      t2530 = 3*t456*t2528
      t2531 = t2163*t1324
      t2534 = 2*t456*t2531*CsThl
      t2535 = t2433*t2368
      t2536 = 24*t2535
      t2537 = t152*t2351
      t2538 = 12*t2537
      t2540 = 6*t2360*t2351
      t2541 = 24*t2495-16*t2497-384*t2499*t2500+224*t2504-160*t2507-128*
     #t2499*t2509+t2514+32*t2503*t2509-t2519+384*t61*t2520-t2525-t2530+t
     #2534+t2536+t2538-t2540
      t2543 = t2436*t2425*t2163
      t2544 = 2*t2543
      t2545 = t152*t2362
      t2546 = 8*t2545
      t2549 = 4*t2360*t190*t1324
      t2551 = t2436*t2350*t1324
      t2552 = 3*t2551
      t2553 = t2361*sqq
      t2555 = 8*t2433*t2553
      t2556 = t1756*Mzsq
      t2558 = 12*t2451*t2556
      t2560 = 16*t2433*t2346
      t2561 = t190*t32
      t2563 = 16*t2451*t2561
      t2564 = t2164*t32
      t2565 = t54*t2564
      t2567 = t57*t496
      t2569 = t2164*sqq
      t2570 = t57*t2569
      t2573 = t54*t1368*t32
      t2575 = t141*t38
      t2576 = t2575*t32
      t2577 = t57*t2576
      t2579 = t54*t1404
      t2581 = t57*t801
      t2583 = t57*t818
      t2585 = t2544-t2546+t2549-t2552-t2555+t2558-t2560+t2563+4*t2565+24
     #*t2567-16*t2570-6*t2573+32*t2577-6*t2579+24*t2581+96*t2583
      t2589 = t141*t122
      t2591 = t54*t2589*t38
      t2592 = 8*t2591
      t2593 = t141*t40
      t2596 = 24*t54*t2593*t122
      t2597 = t122*t1
      t2599 = 32*t2597*t907
      t2601 = 128*t2253*t907
      t2603 = 2*t1171*t2528
      t2604 = CsThr*t2195
      t2605 = t2604*t2523
      t2606 = CsThr*t2163
      t2607 = t1330*CsThl
      t2608 = t2606*t2607
      t2610 = 3*t2248*t2607
      t2614 = 3*t456*t437*t2526*CsThl
      t2618 = t57*t141
      t2619 = t2618*t826
      t2621 = t2618*t260
      t2624 = 18*t205*t500
      t2625 = t54*t141
      t2626 = t2625*t336
      t2628 = t514*t798
      t2630 = -t2592-t2596+t2599+t2601+t2603-t2605-t2608-t2610-t2614+12*
     #t54*t1121*t2250-72*t2619+32*t2621+t2624+24*t2626-72*t2628
      t2634 = t2618*t838
      t2636 = t765*t38
      t2637 = t205*t2636
      t2638 = 6*t2637
      t2639 = t2388*t38
      t2640 = t205*t2639
      t2641 = 4*t2640
      t2642 = t320*t2163
      t2644 = 12*t205*t2642
      t2645 = t2625*t798
      t2650 = t61*t122
      t2651 = t2625*t2650
      t2653 = t40*t2163
      t2654 = t2653*sqq
      t2655 = t514*t2654
      t2657 = t514*t846
      t2659 = t38*t2163
      t2660 = t2659*sqq
      t2661 = t514*t2660
      t2663 = t54*t2365
      t2666 = t2625*t846
      t2670 = t254*t122
      t2673 = t246*t122
      t2676 = -32*t2618*t1062*t32-24*t2634+t2638-t2641-t2644+18*t2645+8*
     #t2625*t1062*t122-8*t2651+48*t2655-24*t2657+16*t2661-16*t2663*t119+
     #6*t2666+16*t2663*t1146+8*t2625*t2670-8*t2625*t2673
      t2678 = t254*t32
      t2681 = t80*t141
      t2690 = 128*t177*t2361*CsThl*t2335*t23
      t2693 = t1595*t2205
      t2696 = 4*t2693*t2227*t436
      t2697 = t2340*t1879
      t2700 = t254*t438
      t2703 = t2163*t437
      t2704 = t1062*t2703
      t2707 = t1062*t438
      t2710 = t61*t438
      t2711 = t2625*t2710
      t2713 = t58*t2703
      t2714 = t514*t2713
      t2716 = t68*t32
      t2717 = t834*t893
      t2718 = t2716*t2717
      t2720 = t3*t122
      t2721 = t1080*t1171
      t2722 = t2720*t2721
      t2724 = t2716*t2721
      t2726 = -32*t2618*t2678-64*t2681*t119+64*t2681*t1146+t2690+4*t2356
     #*t2272+t2696+32*t2697*t2343+24*t514*t2700+8*t514*t2704-6*t2625*t27
     #07+6*t2711-8*t2714-128*t2718-24*t2722+64*t2724
      t2727 = t2195*Mzsq
      t2728 = t254*t2727
      t2731 = t2127*t2713
      t2733 = t246*t447
      t2736 = t246*t2703
      t2739 = t254*t447
      t2742 = t58*t2727
      t2743 = t2127*t2742
      t2745 = t246*t2727
      t2753 = 8*t2340*t38*t1695
      t2754 = t177*CsThl
      t2755 = t2604*t427
      t2756 = t2754*t2755
      t2758 = t177*t38
      t2759 = t40*t2195
      t2760 = t2759*t427
      t2761 = t2758*t2760
      t2763 = t2606*t443
      t2764 = t2754*t2763
      t2766 = t2653*t443
      t2767 = t2758*t2766
      t2769 = t177*t235
      t2770 = t240*t2195
      t2771 = t2770*t427
      t2774 = t165*t38
      t2775 = t320*t2727
      t2777 = 4*t2774*t2775
      t2778 = t240*t2163
      t2779 = t2778*t443
      t2782 = -32*t2132*t2728+12*t2731+8*t2720*t2733-4*t2127*t2736+64*t2
     #716*t2739+12*t2743-4*t2127*t2745+16*t2310*t235*t1661+t2753-56*t275
     #6+40*t2761-56*t2764+40*t2767-8*t2769*t2771-t2777-8*t2769*t2779
      t2785 = t320*t2703
      t2787 = 4*t2774*t2785
      t2788 = t2310*CsThl
      t2789 = t2788*t1647
      t2791 = t2310*t38
      t2792 = t2791*t1655
      t2794 = t1432*t1171
      t2796 = 48*t1404*t2794
      t2798 = 12*t2564*t2794
      t2799 = t1356*t447
      t2801 = 24*t1433*t2799
      t2805 = t2132*t2713
      t2807 = t61*t2703
      t2808 = t2132*t2807
      t2810 = t254*t2703
      t2819 = t2132*t2742
      t2821 = t61*t2727
      t2822 = t2132*t2821
      t2824 = -t2787+112*t2789-80*t2792-t2796+t2798-t2801+384*t488*t747*
     #t443-32*t2805+64*t2808-32*t2132*t2810-16*t2720*t2739+8*t2127*t2728
     #+8*t2127*t2810-32*t2819+64*t2822
      t2825 = t118*t1324
      t2827 = 4*t2825*t2717
      t2829 = 2*t2140*t2821
      t2831 = 2*t2140*t2807
      t2840 = t80*t32
      t2847 = t246*t438
      t2850 = t246*t2404
      t2853 = t1971*t2807
      t2855 = t246*t863
      t2858 = t371*t1695
      t2862 = -t2827+t2829+t2831+16*t2134*t2713-16*t2134*t2736+16*t2134*
     #t2742-16*t2134*t2745-32*t2840*t2721+32*t2840*t2733+24*t514*t2736+7
     #2*t514*t2847-48*t514*t2850+64*t2853-12*t371*t2855-8*t2858-32*t1971
     #*t2810
      t2866 = t61*t863
      t2869 = t61*t2384
      t2872 = t2433*Er
      t2875 = t2436*Er
      t2878 = 6*t2875*t61*t1330
      t2881 = 4*t2875*t61*t2376
      t2882 = t1971*t2821
      t2884 = t2360*Er
      t2886 = 6*t2884*t2866
      t2887 = t254*t863
      t2890 = t254*t2384
      t2897 = t54*t122
      t2900 = t57*t32
      t2903 = t57*sqq
      t2906 = t1062*t2727
      t2909 = -96*t1971*t2700-72*t371*t2866+48*t371*t2869+24*t2872*t2700
     #+t2878-t2881+64*t2882+t2886+24*t371*t2887-16*t371*t2890-2*t2138*t2
     #810+2*t2138*t2736-4*t2897*t2733-48*t2900*t2733+8*t2903*t2728+8*t29
     #03*t2906
      t2910 = t1062*t447
      t2915 = t2138*t2742
      t2917 = t2897*t2717
      t2919 = t2138*t2821
      t2921 = t2903*t2742
      t2923 = t2903*t2821
      t2931 = t2903*t2713
      t2933 = t2138*t2713
      t2935 = t2138*t2807
      t2937 = t2903*t2807
      t2939 = t2897*t2721
      t2943 = 4*t2897*t2910-2*t2138*t2906-6*t2915-4*t2917+2*t2919-8*t292
     #1-8*t2923+8*t2903*t2810+8*t2903*t2704-2*t2138*t2704-8*t2931-6*t293
     #3+2*t2935-8*t2937+12*t2939-16*t2900*t2910
      t2958 = t58*t863
      t2959 = t205*t2958
      t2961 = t58*t2384
      t2962 = t205*t2961
      t2964 = t205*t1695
      t2966 = t246*t443
      t2973 = t246*t2384
      t2976 = t58*t438
      t2977 = t2872*t2976
      t2979 = t67*t141
      t2980 = t2979*Er
      t2981 = t61*t443
      t2982 = t2980*t2981
      t2984 = t371*t2958
      t2986 = -16*t2900*t2739+4*t2897*t2739-2*t2138*t2728+2*t2138*t2745-
     #6*t2625*t2847-18*t2959+12*t2962-24*t2964+72*t2618*t2966-48*t2618*t
     #2733+24*t514*t2745+8*t371*t2973+24*t2977+192*t2982+36*t2984
      t2987 = t371*t2961
      t2989 = t1971*t2710
      t2991 = t61*t2404
      t2992 = t1971*t2991
      t2994 = t254*t2404
      t3000 = 2*t2884*t1695
      t3001 = t2980*t2717
      t3009 = t2451*Er
      t3014 = t254*t443
      t3019 = t58*t443
      t3024 = -24*t2987+192*t2989-128*t2992+64*t1971*t2994-32*t1971*t272
     #8-t3000-128*t3001-72*t2872*t2710+4*t2872*t2736-4*t2872*t2713+12*t3
     #009*t2966-8*t3009*t2733-96*t2980*t3014+64*t2980*t2739-12*t3009*t30
     #19+8*t3009*t2721
      t3026 = t2872*t1655
      t3032 = t58*t2404
      t3037 = t58*t2408
      t3038 = t514*t3037
      t3040 = t514*t709
      t3042 = t514*t1466
      t3044 = t254*t427
      t3051 = t246*t427
      t3054 = t254*t665
      t3057 = t246*t665
      t3060 = t246*t2408
      t3063 = 8*t3026-8*t2872*t2850+4*t2872*t2745+8*t2872*t3032-4*t2872*
     #t2742+16*t3038-24*t3040-24*t3042+24*t2618*t3044+12*t2663*t699-8*t2
     #625*t3037-12*t2663*t3051-6*t2625*t3054-6*t2625*t3057+8*t2625*t3060
      t3065 = t54*t141*t235
      t3068 = t152*t196
      t3069 = t3068*t2763
      t3071 = Er*t240
      t3078 = t58*t765
      t3079 = t205*t3078
      t3085 = t371*t3078
      t3087 = t1971*t1466
      t3097 = t2659*t41
      t3098 = t1971*t3097
      t3100 = t2388*t58
      t3101 = t371*t3100
      t3103 = t1971*t709
      t3105 = 4*t3065*t1661+8*t3069-8*t2433*t3071*t2250+64*t80*t474*t225
     #0-18*t3079-48*t514*t3060+72*t2618*t3051+36*t3085+192*t3087+8*t2872
     #*t3037-8*t2872*t3060+12*t3009*t3051-12*t3009*t699-128*t3098-24*t31
     #01-96*t3103
      t3108 = t2980*t3019
      t3110 = t1971*t2976
      t3112 = t1971*t3032
      t3138 = -96*t3108-96*t3110+64*t3112+48*t1914*t2976-32*t1914*t3032+
     #16*t1914*t2713-16*t1914*t2736+16*t1914*t2742-16*t1914*t2745-48*t19
     #14*t2847+32*t1914*t2850+48*t2681*t3019-32*t2681*t2721-48*t2681*t29
     #66+32*t2681*t2733
      t3139 = t499*t447
      t3141 = 2*t2875*t3139
      t3142 = t2418*t235
      t3145 = t2597*t38
      t3146 = t893*t416
      t3148 = 16*t3145*t3146
      t3149 = t2253*t235
      t3152 = t2253*CsThl
      t3155 = t2418*CsThl
      t3156 = t3155*t1658
      t3158 = t845*t14
      t3159 = t3158*t3146
      t3161 = t1971*t2713
      t3163 = t1971*t2742
      t3169 = t2900*t2721
      t3171 = t2900*t2717
      t3173 = t2618*t3019
      t3178 = t2618*t2981
      t3180 = t3141+32*t3142*t1664+t3148-128*t3149*t1664-384*t3152*t1658
     #+224*t3156-160*t3159-32*t3161-32*t3163+48*t2720*t2717-24*t2127*t28
     #07+16*t3169+16*t3171-24*t3173+24*t2618*t1062*t443-24*t3178
      t3182 = t2618*t2721
      t3185 = t205*t438*t834
      t3187 = t205*t2760
      t3189 = t38*t2195
      t3191 = t205*t3189*t427
      t3199 = t205*t2766
      t3202 = t205*t2659*t443
      t3204 = t205*t2866
      t3206 = t205*t2869
      t3208 = t2618*t2717
      t3220 = 16*t3182-8*t3185+6*t3187+2*t3191-6*t205*t1062*t863+4*t205*
     #t1062*t2384+6*t3199+2*t3202+6*t3204-4*t3206+16*t3208-6*t205*t2887+
     #4*t205*t2890+6*t205*t2855-4*t205*t2973-16*t514*t2994
      t3223 = t2625*t1655
      t3230 = t834*t443
      t3231 = t514*t3230
      t3233 = t2625*t3230
      t3234 = 4*t3233
      t3235 = t514*t2742
      t3237 = t514*t2821
      t3244 = 192*t618*t2223
      t3246 = 128*t618*t2229
      t3247 = t2310*t1875
      t3249 = 256*t3247*t976
      t3250 = t2205*t559
      t3252 = 192*t3247*t3250
      t3254 = t165*t2368*t609
      t3259 = 192*t1496*t2335*t17*t32
      t3260 = 8*t514*t2728-12*t3223-16*t514*t1062*t2404+8*t514*t2906+32*
     #t3231-t3234-8*t3235-8*t3237+24*t514*t2707+8*t514*t2810-t3244+t3246
     #+t3249-t3252-24*t3254-t3259
      t3266 = 256*t1773*t436*t3250
      t3268 = t585*t1550*Mzsq
      t3269 = t2356*t3268
      t3271 = Eq*t1324
      t3272 = t3271*t2205
      t3274 = t9*sll*t436
      t3276 = 12*t3272*t3274
      t3277 = t558*t251
      t3279 = 12*t3272*t3277
      t3280 = t558*t436
      t3282 = 8*t3272*t3280
      t3285 = 8*t2693*t2227*t251
      t3286 = t2418*t26
      t3287 = t341*t887
      t3290 = sll*t235
      t3291 = t54*t3290
      t3297 = t54*t250
      t3298 = t3297*t2755
      t3300 = t172*t1695
      t3302 = t172*t2760
      t3306 = t3297*t2763
      t3308 = t341*t898
      t3309 = t2162*t3308
      t3311 = t3266+24*t3269+t3276+t3279-t3282-t3285-288*t3286*t3287+8*t
     #3291*t1171*t665-2*t3291*t2755-2*t3298-8*t3300+2*t3302-2*t3291*t276
     #3-2*t3306+72*t3309
      t3312 = t3271*t352
      t3313 = t17*t887
      t3315 = 24*t3312*t3313
      t3316 = t17*t447
      t3318 = 16*t3312*t3316
      t3319 = t1595*t352
      t3320 = t17*t2478
      t3322 = 16*t3319*t3320
      t3323 = t17*t2727
      t3325 = 8*t3319*t3323
      t3326 = t273*t275
      t3327 = t362*t3326
      t3328 = t276*t278
      t3332 = t2264*t352
      t3337 = t73*t26
      t3338 = t341*t2478
      t3344 = t15*t115
      t3349 = t2256*t583
      t3352 = t1992*t1951
      t3355 = t91*t488
      t3358 = t148*SnThl
      t3359 = t91*t3358
      t3366 = -t3315+t3318+t3322-t3325+64*t3327*t3328*t2478+192*t3332*t3
     #320-96*t3332*t3323+192*t3337*t3338-128*t15*t260*t2312+64*t3344*t22
     #19+4*t2320*t2272-384*t3349*t615-64*t3352*t1647-8*t3355*t2233-24*t3
     #359*t17*t944-24*t3359*t17*t846
      t3371 = t341*t2727
      t3374 = t17*t1368
      t3375 = t2306*t3374
      t3377 = t17*t2164
      t3378 = t2306*t3377
      t3381 = t165*t2355*SnThl
      t3382 = t17*t2405
      t3383 = t3381*t3382
      t3385 = t1*t2365
      t3386 = t3385*t2355
      t3394 = 128*t3247*t2312
      t3396 = t416*t1993*t2343
      t3398 = t2164*t1993
      t3399 = t3398*t976
      t3401 = t3398*t3250
      t3403 = t76*t27
      t3409 = 16*t3312*t17*t2500
      t3413 = 16*t3359*t17*t2660-96*t3337*t3371+144*t3375-96*t3378+8*t33
     #83+12*t3386*t564+12*t3386*t560-8*t3386*t615+t3394-192*t3396+128*t3
     #399-64*t3401-24*t3403*t602*t2500-t3409+16*t3403*t2186*t251
      t3414 = t602*t437
      t3421 = t165*t1879
      t3422 = t1632*t2342
      t3429 = El*t2195*t9*t695
      t3436 = t396*t995
      t3439 = t2340*Er
      t3443 = t2232*t980
      t3448 = t2218*t980
      t3451 = t2163*CsThl
      t3452 = t195*El
      t3454 = t3451*t3452*t696
      t3459 = t158*t2180
      t3460 = t8*t2201
      t3464 = t17*t448
      t3465 = t2433*t352*t3464
      t3467 = 32*t3403*t3414*t251-8*t3403*t2186*t1301-24*t3421*t3422-24*
     #t2697*t2223+16*t3421*t3429-24*t3421*t2237-24*t2697*t2207+32*t3436*
     #t2245+16*t3439*t995*t981-8*t3436*t3443+16*t3421*t2241-8*t3436*t344
     #8-64*t3454-192*t2169*t341*t904+96*t3459*t3460-96*t3465
      t3472 = t2340*t26*t341*t448
      t3474 = t165*t26
      t3476 = t3474*t341*t2409
      t3478 = t165*t27
      t3482 = t2436*t33
      t3484 = 8*t3482*t2293
      t3485 = t152*t3326
      t3486 = t3328*t2703
      t3489 = sll*t67
      t3490 = t3489*t352
      t3491 = t17*t2703
      t3494 = t177*t26
      t3495 = t341*t2703
      t3515 = t17*t898
      t3518 = -48*t3472+24*t3476-8*t3478*t2197*t436-t3484-32*t3485*t3486
     #-96*t3490*t3491-96*t3494*t3495-96*t3327*t3328*t898-96*t3327*t3328*
     #t1368+64*t3327*t3328*t2164-32*t3327*t3486-288*t3332*t3374+192*t333
     #2*t3377-96*t3332*t3491-288*t3332*t3515
      t3521 = t341*t1368
      t3530 = t352*t17*t141
      t3546 = t2306*t3491
      t3548 = t2306*t3320
      t3550 = t2306*t3323
      t3552 = t2256*t352
      t3553 = t17*t416
      t3556 = t341*t907
      t3557 = t2162*t3556
      t3559 = -288*t3337*t3308-288*t3337*t3521+192*t3337*t2165-96*t3337*
     #t3495-768*t2256*t58*t3530+576*t3349*t571-288*t2265*t1023+192*t2265
     #*t2290+384*t2265*t2278-96*t2265*t2293+384*t2256*t61*t3530+48*t3546
     #-96*t3548+48*t3550-288*t3552*t3553-96*t3557
      t3561 = t2900*t935
      t3563 = t2900*t918
      t3566 = 12*t2897*t918
      t3567 = t2872*t448
      t3568 = 24*t3567
      t3569 = t2759*sqq
      t3571 = 8*t2872*t3569
      t3572 = t2759*Mzsq
      t3574 = 4*t2872*t3572
      t3576 = 12*t3009*t924
      t3578 = 8*t3009*t918
      t3579 = t2653*t437
      t3581 = 4*t2872*t3579
      t3582 = t2875*t1473
      t3583 = 6*t3582
      t3585 = 2*t2875*t2385
      t3587 = 2*t2875*t2389
      t3590 = 2*t2875*t499*t2195
      t3593 = 3*t2875*t893*t1324
      t3594 = t371*t904
      t3595 = 24*t3594
      t3596 = t2884*t904
      t3597 = 6*t3596
      t3598 = 16*t3561+48*t3563-t3566-t3568+t3571-t3574-t3576+t3578-t358
     #1+t3583-t3585-t3587+t3590-t3593-t3595+t3597
      t3599 = t141*t233
      t3606 = t2366*Er
      t3608 = 2*t3606*t448
      t3609 = t2884*t2405
      t3610 = t2884*t2409
      t3611 = t2759*t32
      t3613 = 2*t2884*t3611
      t3615 = 3*t3606*t912
      t3616 = t2512*t38
      t3618 = 12*t3616*t921
      t3620 = 8*t3616*t918
      t3621 = t2503*t235
      t3626 = t2770*sll
      t3629 = t2770*Mzsq
      t3632 = t2759*sll
      t3634 = 8*t3145*t3632
      t3636 = 4*t3145*t3572
      t3637 = t2778*t437
      t3640 = -64*t3599*t235*t1212-192*t3599*CsThl*t1178+t3608-t3609-t36
     #10-t3613+t3615-t3618+t3620-24*t3621*t1207+16*t3621*t1212+16*t3142*
     #t3626-8*t3142*t3629+t3634-t3636-8*t3142*t3637
      t3649 = 12*t3145*t915
      t3651 = 4*t3145*t3579
      t3652 = t751*t235
      t3655 = t751*CsThl
      t3656 = t2606*t437
      t3661 = t2604*Mzsq
      t3664 = t2024*CsThr
      t3666 = 24*t2523*t3664
      t3669 = 12*t2379*CsThl*t3664
      t3670 = t233*t2163
      t3671 = t3670*t437
      t3674 = t233*t2195
      t3675 = t3674*Mzsq
      t3681 = t2195*t1*Mzsq
      t3683 = 12*t336*t3681
      t3685 = 4*t2127*t3572
      t3686 = t2625*t448
      t3687 = 4*t3686
      t3688 = -24*t3142*t1230-t3649-t3651+32*t3652*t3637+96*t3655*t3656+
     #32*t3652*t3629+96*t3655*t3661-t3666+t3669-96*t115*t3671-96*t115*t3
     #675+192*t333*t904+t3683+t3685+t3687
      t3689 = t2980*t119
      t3692 = 3*t3606*t798
      t3694 = 2*t2884*t2642
      t3696 = 12*t3616*t829
      t3701 = t2778*sll
      t3705 = 12*t3145*t808
      t3706 = t2653*sll
      t3708 = 8*t3145*t3706
      t3714 = t2305*t995
      t3725 = 128*t3689+t3692-t3694-t3696-24*t3621*t475-24*t3142*t1125+1
     #6*t3142*t3701-t3705+t3708-96*t2618*t1146-8*t165*t2362*t1487+48*t37
     #14*t1023-32*t3714*t2290-64*t3714*t2278+16*t3714*t2293+48*t3714*t10
     #20
      t3731 = t2597*t26
      t3737 = 12*t2693*t551*t2500
      t3740 = 8*t2693*t2210*t251
      t3743 = 4*t2693*t2210*t1301
      t3746 = 12*t2693*t558*t2500
      t3747 = t2*t122
      t3751 = t3747*t995
      t3768 = -32*t3714*t2269+16*t3714*t2272+32*t3731*t341*t2500+t3737-t
     #3740+t3743+t3746-64*t3747*t61*t3530+48*t3751*t568+48*t3751*t571-32
     #*t3751*t615+16*t83*t2650*t2312-8*t1035*t2233-8*t1035*t2219-16*t91*
     #t33*t615
      t3770 = t91*t174*El
      t3771 = t86*t2211
      t3774 = t91*t2296
      t3782 = t91*t2593*El
      t3796 = 64*t2499*t447
      t3797 = t2418*t898
      t3798 = 72*t3797
      t3800 = 24*t2597*t898
      t3801 = t2418*t2703
      t3802 = 24*t3801
      t3804 = 8*t2597*t2703
      t3806 = 32*t751*t2703
      t3807 = 48*t3770*t3771-16*t3774*t615-32*t3359*t3464+48*t3770*t2357
     #-72*t3782*t564+48*t3782*t615-72*t3782*t560+96*t3770*t981+24*t3774*
     #t564+24*t3774*t560+t3796+t3798-t3800+t3802-t3804-t3806
      t3811 = 32*t751*t2727
      t3813 = 64*t3599*t447
      t3815 = 4*t1325*t447
      t3816 = t769*t2703
      t3817 = t769*t2727
      t3819 = 2*t2589*t447
      t3822 = 4*t2248*t2526*CsThl
      t3827 = 4*t2366*t2425
      t3828 = t54*t2589
      t3830 = t57*t1356
      t3833 = 64*t2253*t2164
      t3835 = 24*t2597*t1368
      t3837 = 16*t2597*t2164
      t3838 = t2418*t1368
      t3839 = 72*t3838
      t3840 = -t3811+t3813-t3815+t3816+t3817-t3819+t3822-96*t57*t474*t22
     #50+t3827+8*t3828-32*t3830+t3833-t3835+t3837+t3839
      t3841 = t2418*t2164
      t3842 = 48*t3841
      t3844 = 24*t2512*t416
      t3845 = t2503*t416
      t3846 = 72*t3845
      t3848 = 96*t2499*t416
      t3849 = t2433*t2425
      t3850 = 16*t3849
      t3851 = t2360*t2437
      t3852 = 4*t3851
      t3854 = 16*t2451*t2367
      t3856 = 32*t2253*t2703
      t3858 = 96*t2253*t898
      t3860 = 64*t2253*t2478
      t3862 = 32*t2253*t2727
      t3863 = t2418*t2478
      t3864 = 48*t3863
      t3865 = t2418*t2727
      t3866 = 24*t3865
      t3868 = 16*t2597*t2478
      t3870 = 8*t2597*t2727
      t3871 = t2503*t887
      t3872 = 72*t3871
      t3873 = -t3842-t3844+t3846-t3848-t3850+t3852-t3854-t3856-t3858+t38
     #60-t3862-t3864+t3866+t3868-t3870+t3872
      t3875 = t2503*t447
      t3876 = 48*t3875
      t3878 = 24*t2512*t887
      t3880 = 16*t2512*t447
      t3882 = 96*t2499*t887
      t3884 = t10*t141
      t3885 = t8*t3884
      t3886 = t2597*t4*t3885
      t3890 = 32*t3271*El*t673
      t3895 = t3747*El*t673
      t3917 = -t3876-t3878+t3880-t3882-96*t3886+t3890+384*t2256*El*t673-
     #192*t3895+384*t2418*t4*t3885+128*t2305*t273*t277*t278*t141+48*t191
     #4*t709-32*t1914*t3037+32*t1914*t3060-48*t1914*t3057+48*t2681*t699-
     #48*t2681*t3051
      t3918 = t2872*t709
      t3920 = t61*t427
      t3921 = t2980*t3920
      t3923 = t1971*t3037
      t3925 = t2980*t699
      t3927 = t514*t3097
      t3929 = t2618*t699
      t3934 = t205*t3100
      t3942 = t61*t765
      t3943 = t205*t3942
      t3945 = t61*t2388
      t3946 = t205*t3945
      t3948 = t2618*t3920
      t3950 = t254*t765
      t3953 = t254*t2388
      t3956 = t2388*t246
      t3959 = 24*t3918+192*t3921+64*t3923-96*t3925+16*t3927-24*t3929+24*
     #t2618*t1062*t427+12*t3934-6*t205*t1062*t765+4*t205*t1062*t2388+6*t
     #3943-4*t3946-24*t3948-6*t205*t3950+4*t205*t3953-4*t205*t3956
      t3963 = t765*t246
      t3966 = t2625*t709
      t3968 = t1062*t665
      t3971 = t2625*t1466
      t3980 = t254*t2408
      t3997 = 6*t2875*t61*t1328
      t3998 = 6*t205*t3963-6*t3966-6*t2625*t3968+6*t3971+24*t514*t3968-1
     #6*t514*t1062*t2408+24*t514*t3054-16*t514*t3980+72*t514*t3057-12*t3
     #71*t3963+8*t371*t3956+64*t1971*t3980+48*t371*t3945+24*t2872*t3054+
     #t3997
      t4001 = 4*t2875*t61*t2379
      t4011 = 6*t2884*t3942
      t4018 = 12*t251*t2384*t1755
      t4021 = t2980*t2721
      t4023 = t2884*t2766
      t4024 = t2884*t2760
      t4026 = 2*t3606*t1655
      t4031 = -t4001-72*t2872*t1466-96*t2980*t3044-72*t371*t3942-96*t197
     #1*t3054+t4011+24*t371*t3950-16*t371*t3953+t4018-24*t2127*t2821+64*
     #t4021+t4023+t4024-t4026+24*t2903*t2745+24*t2903*t2736
      t4033 = t514*t2807
      t4035 = t514*t2976
      t4037 = t2625*t2976
      t4039 = t514*t1655
      t4041 = t514*t3032
      t4043 = t514*t2710
      t4045 = t514*t2991
      t4063 = -8*t4033-24*t4035-6*t4037+96*t4039+16*t4041-24*t4043+16*t4
     #045+24*t2618*t3014-16*t2618*t2739-16*t2618*t2910-8*t2625*t3032+4*t
     #2625*t2742+12*t2663*t3019-8*t2663*t2721+4*t2625*t2713
      t4077 = 16*t3319*t3377
      t4079 = 8*t3319*t3491
      t4081 = 24*t3319*t3515
      t4084 = t58*t22
      t4086 = t4084*t17*t1473
      t4089 = t4084*t17*t2385
      t4101 = t152*t995
      t4104 = -12*t2663*t2966+8*t2663*t2733+8*t2625*t2850-4*t2625*t2745-
     #6*t2625*t2700-4*t2625*t2736+t4077-t4079-t4081+96*t3459*t2192+128*t
     #4086-64*t4089+16*t2340*t27*t603*t1301-8*t3478*t2187*t1301-32*t2433
     #*t995*t981+16*t4101*t3448
      t4109 = t2436*t34
      t4113 = 4*t4109*t17*t2163*t1301
      t4114 = t17*t2195
      t4117 = 4*t4109*t4114*t436
      t4119 = 4*t141*t2526
      t4122 = t2305*t3326
      t4123 = t3328*t447
      t4128 = 3*sll*t2526*t437
      t4129 = t152*t1878
      t4138 = Er*t235
      t4142 = t3068*t2755
      t4148 = 4*t2360*t1878*t1695
      t4150 = t2433*t196*t1647
      t4152 = 16*t4101*t3443+t4113+t4117-t4119+576*t3349*t568+64*t4122*t
     #4123+t4128+96*t4129*t1695-24*t4129*t2766+8*t152*t3071*t1080*t438+1
     #28*t3489*t4138*t1661+8*t4142-24*t4129*t2760-t4148-8*t4150
      t4154 = t54*t2593*t3230
      t4156 = t2103*t2976
      t4158 = t57*t3290
      t4163 = t1203*t665
      t4174 = t385*t171
      t4175 = t4174*t1647
      t4177 = t152*t4138
      t4182 = t68*t416
      t4183 = t4182*t2981
      t4186 = t152*t2350*t2976
      t4188 = t4182*t3019
      t4192 = -4*t4154+16*t4156-32*t4158*t1647-32*t4158*t1661+8*t3291*t4
     #163-2*t3291*t2771-2*t3297*t2771-2*t3291*t2779-2*t3297*t2779+32*t41
     #75+8*t4177*t2771+8*t4177*t2779-256*t4183-40*t4186+128*t4188+64*t33
     #44*t2233
      t4202 = t165*t2425*t609
      t4206 = t165*t2367*t352*t4114
      t4208 = t1412*t2343
      t4210 = t2356*t696
      t4212 = t190*El
      t4213 = t165*t4212
      t4216 = t2340*t4212
      t4222 = 64*t618*t2219
      t4224 = 64*t618*t2233
      t4227 = t2205*t563
      t4229 = 192*t3247*t4227
      t4230 = t2433*t1875
      t4238 = -32*t83*t336*t2312+16*t607*t2219+16*t607*t2233+12*t4202-8*
     #t4206+12*t4208+24*t4210-8*t4213*t2245+8*t4216*t2357-4*t4216*t2272-
     #t4222-t4224-8*t2356*t3771-t4229+128*t4230*t26*t341*t32-96*t4230*t2
     #7*t603
      t4242 = t27*t3414
      t4251 = t2618*t119
      t4256 = 128*t2520
      t4258 = 2*t1325*t2195
      t4259 = t3328*t2727
      t4262 = t2306*t3515
      t4268 = t3747*t352
      t4269 = t4268*t3553
      t4271 = t4268*t3313
      t4273 = t4268*t3316
      t4275 = -96*t4230*t27*t597+64*t4230*t4242-256*t2311*t976+192*t2311
     #*t4227+192*t2311*t3250+32*t4251-24*t3403*t596*t2500+t4256-t4258-32
     #*t3485*t4259+144*t4262-288*t3552*t3313+192*t3552*t3316+144*t4269+1
     #44*t4271-96*t4273
      t4280 = t341*t416
      t4283 = t341*t447
      t4286 = t2162*t3371
      t4294 = t3731*t4280
      t4296 = t3731*t3287
      t4298 = t3731*t4283
      t4301 = 32*t3319*t2307
      t4309 = t2597*t27
      t4320 = -288*t3286*t4280+192*t3286*t4283+24*t4286-96*t4122*t3328*t
     #416-96*t4122*t3328*t887+72*t4294+72*t4296-48*t4298+t4301+128*t3327
     #*t3328*t907+384*t3332*t2307+384*t3337*t3556-24*t4309*t10*sll*t436-
     #24*t4309*t596*t251+16*t4309*t596*t436
      t4323 = 16*t2693*t979*t251
      t4324 = t250*t3452
      t4327 = t4324*t86*t558*t122
      t4329 = t4324*t2245
      t4332 = t86*t2227*t32
      t4333 = t4324*t4332
      t4335 = t4324*t3443
      t4338 = t58*t2433*t2343
      t4341 = t58*t2433*El
      t4342 = t4341*t3268
      t4344 = t4341*t696
      t4346 = t4341*t981
      t4353 = t86*t392*t437
      t4357 = t91*t2408*El
      t4361 = t86*t229*t437
      t4367 = t91*t115
      t4372 = -t4323-192*t4327+256*t4329+128*t4333-64*t4335+256*t4338-19
     #2*t4342-192*t4344+128*t4346+48*t91*t2319*t615+16*t91*t666*t4353-8*
     #t4357*t4353-24*t4357*t4361+16*t91*t260*t2312-8*t4367*t2219-8*t4367
     #*t2233
      t4380 = t91*t2164*El
      t4385 = t2531*t437
      t4386 = t141*t1324
      t4388 = 3*t4386*t437
      t4390 = 2*t2164*t1324
      t4391 = t2503*t141
      t4392 = 96*t4391
      t4393 = t2421*Mzsq
      t4395 = 2*t2527*t437
      t4397 = 3*t4386*Mzsq
      t4399 = 3*t416*t2526
      t4400 = 32*t2517
      t4410 = 32*t91*t416*El*t87*t928-8*t4380*t4353-24*t4380*t4361+t4385
     #+t4388-t4390-t4392+t4393-t4395+t4397+t4399+t4400+24*t91*t416*t976+
     #96*t91*t2593*t976-72*t91*t811*t609
      t4414 = t91*t148
      t4437 = t585*SnThr*t38*Mzsq
      t4441 = t91*t2478*El
      t4449 = t86*t229*Mzsq
      t4459 = t91*t2404*El
      t4464 = -72*t91*t801*t1487-16*t4414*t2357+32*t91*t1356*t576-32*t37
     #74*t1551+24*t91*t184*t609+8*t91*t2164*t609-16*t4414*t3771+8*t91*t2
     #478*t2087+16*t91*t2296*t437*t4437+16*t4441*t86*t392*sqq-8*t4441*t4
     #437-24*t4441*t4449+8*t91*t2408*t609+8*t91*t2404*t2087-8*t4459*t443
     #7-24*t4459*t4449
      t4467 = t1389*t2293
      t4471 = t27*t2191
      t4476 = t91*t2506
      t4489 = t91*sll
      t4501 = 4*t4467+64*t2052*t4242-32*t113*t4471-32*t113*t2202+32*t447
     #6*t976-24*t4476*t4227-24*t4476*t3250+16*t4476*t2312-24*t3355*t2207
     #+16*t3355*t2212+32*t4489*t995*t981-8*t3355*t2219-24*t3355*t2223+16
     #*t3355*t2229+16*t2244*t3448
      t4502 = t2340*t1875
      t4521 = t3385*t4212
      t4536 = 8*t2360*t34*t17*Mzsq*t1301
      t4540 = 36*t2500*t766
      t4541 = -64*t4502*t2343+48*t4502*t2207+48*t4502*t2223-32*t3439*t58
     #3*t981+48*t1836*t3422-32*t1836*t3429+16*t2244*t3443+16*t3385*t190*
     #t976-12*t4521*t564-12*t4521*t560+8*t4521*t615+8*t4216*t3771+8*t421
     #6*t981-t4536-4*t4216*t2293+t4540
      t4545 = t1328*t1
      t4547 = 36*t251*t4545
      t4548 = t2379*t1
      t4550 = 24*t251*t4548
      t4563 = t2606*sll
      t4566 = t2499*CsThl
      t4569 = t2503*CsThl
      t4570 = t4569*t537
      t4572 = t2503*t38
      t4573 = t4572*t829
      t4575 = t2499*t235
      t4578 = t3155*t466
      t4580 = t3155*t4563
      t4582 = t3158*t808
      t4584 = -288*t333*t496+t4547-t4550-288*t488*t747*t32+192*t488*t367
     #0*sqq+96*t3149*t1125-64*t3149*t3701+288*t3152*t466-192*t3152*t4563
     #+288*t4566*t537-168*t4570+120*t4573+96*t4575*t475-168*t4578+112*t4
     #580+120*t4582
      t4585 = t3158*t3706
      t4587 = t2872*t336
      t4589 = t2980*t260
      t4592 = 12*t3009*t826
      t4602 = 8*t2872*t2654
      t4605 = 3*t2875*t1328*t40
      t4608 = 2*t2875*t2379*t40
      t4611 = 8*t2884*t61*t1324
      t4617 = 12*t371*t500
      t4619 = 8*t371*t2642
      t4620 = -80*t4585-48*t4587-256*t4589-t4592+96*t2872*t2650-16*t3009
     #*t1146+16*t2872*t2673+16*t3009*t119+t4602-t4605+t4608-t4611+128*t2
     #980*t2678-32*t2872*t2670+t4617-t4619
      t4630 = 2*t2248*t122*t1181
      t4631 = t371*t2405
      t4632 = 8*t4631
      t4634 = 12*t371*t894
      t4636 = 8*t371*t3611
      t4637 = t371*t2409
      t4638 = 8*t4637
      t4640 = 8*t2720*t918
      t4642 = 4*t2127*t3579
      t4644 = 4*t1175*t2523
      t4645 = t456*t2195
      t4646 = t765*CsThl
      t4647 = t4645*t4646
      t4648 = t863*CsThl
      t4649 = t4563*t4648
      t4650 = t205*t904
      t4653 = 12*t205*t3611
      t4654 = t514*t912
      t4656 = -288*t2506*t860+192*t2506*t857+t4630+t4632+t4634-t4636+t46
     #38-t4640+t4642+t4644-t4647-t4649+8*t4650-t4653-72*t4654
      t4657 = t2618*t918
      t4659 = t2618*t924
      t4661 = t514*t3569
      t4663 = t2625*t912
      t4665 = t514*t448
      t4668 = 18*t205*t894
      t4669 = t2618*t928
      t4671 = t2618*t935
      t4673 = t863*t38
      t4674 = t205*t4673
      t4675 = 6*t4674
      t4676 = t205*t2405
      t4678 = t3189*t32
      t4679 = t205*t4678
      t4680 = 4*t4679
      t4681 = t205*t2409
      t4683 = t3189*sqq
      t4684 = t514*t4683
      t4686 = t3189*Mzsq
      t4687 = t514*t4686
      t4689 = t514*t3572
      t4691 = t2659*t437
      t4692 = t514*t4691
      t4694 = 48*t4657-72*t4659+48*t4661+18*t4663-32*t4665+t4668-24*t466
     #9+16*t4671+t4675-2*t4676-t4680-2*t4681+16*t4684-8*t4687-24*t4689-8
     #*t4692
      t4696 = t514*t3579
      t4698 = t514*t944
      t4700 = t2625*t944
      t4702 = t2589*t437
      t4704 = 36*t1433*t4702
      t4712 = 36*t251*t1330*t1
      t4715 = 24*t251*t2376*t1
      t4733 = -24*t4696-24*t4698+6*t4700+t4704-96*t488*t3671-288*t488*t8
     #52*t32+t4712-t4715+192*t488*t3674*sqq-96*t488*t3675+32*t3149*t3637
     #+96*t3149*t1230+288*t3152*t1227+96*t3152*t3656-64*t3149*t3626+32*t
     #3149*t3629
      t4742 = t4569*t1175
      t4744 = t4569*t1178
      t4746 = t4572*t921
      t4748 = t4572*t918
      t4754 = t3155*t4645
      t4756 = t3155*t3661
      t4758 = t3158*t3632
      t4760 = t3158*t3572
      t4762 = t3155*t1227
      t4764 = t3155*t3656
      t4766 = -192*t3152*t4645+96*t3152*t3661+288*t4566*t1175-192*t4566*
     #t1178-168*t4742+112*t4744+120*t4746-80*t4748+96*t4575*t1207-64*t45
     #75*t1212+112*t4754-56*t4756-80*t4758+40*t4760-168*t4762-56*t4764
      t4769 = t3158*t3579
      t4771 = t3158*t915
      t4773 = t2140*t3572
      t4775 = 2*t2825*t918
      t4776 = t2140*t3579
      t4777 = t2897*t935
      t4778 = 4*t4777
      t4779 = t2903*t4686
      t4782 = 6*t2138*t3572
      t4783 = t2138*t4686
      t4784 = 2*t4783
      t4785 = t2903*t3572
      t4787 = t2903*t4691
      t4790 = 6*t2138*t3579
      t4791 = t2138*t4691
      t4792 = 2*t4791
      t4793 = t2903*t3579
      t4796 = 24*t3312*t3553
      t4797 = 40*t4769+120*t4771-t4773+t4775-t4776-t4778-8*t4779+t4782+t
     #4784-24*t4785-8*t4787+t4790+t4792-24*t4793-t4796
      t4800 = t2162*t3521
      t4804 = t2162*t3338
      t4807 = t585*SnThr*t122
      t4808 = t2356*t4807
      t4813 = t172*t2766
      t4820 = t2436*t1878
      t4822 = 8*t4820*t3139
      t4824 = 2*t4820*t2775
      t4826 = 2*t4820*t2785
      t4830 = t3489*t995
      t4838 = -32*t3327*t4259+72*t4800-32*t4177*t4163-48*t4804-16*t4808+
     #48*t2433*t1878*t1655+2*t4813+4*t3065*t1647-16*t2433*t4138*t1661-t4
     #822+t4824+t4826+192*t2979*t995*t615-96*t4830*t2293-96*t4830*t2272-
     #384*t2979*t583*t615
      t4840 = t3489*t583
      t4846 = 8*t3482*t2272
      t4851 = t152*t352
      t4852 = t17*t2409
      t4853 = t4851*t4852
      t4855 = t4851*t3382
      t4858 = t3474*t341*t2405
      t4862 = 16*t2360*t33*t615
      t4873 = t165*t2361*t976
      t4879 = t2340*t2355*Mzsq*t560
      t4883 = 192*t4840*t2272+192*t4840*t2293-t4846-96*t3490*t3323-96*t3
     #494*t3371+48*t4853+48*t4855+24*t4858+t4862+64*t2433*t3326*t4123+19
     #2*t2979*t352*t3316+192*t2310*t26*t4283-8*t4873-16*t3386*t1551-24*t
     #4879+4*t2356*t2293
      t4884 = t3381*t4852
      t4893 = t11*t448
      t4898 = t596*t32
      t4911 = 24*t3319*t3374
      t4912 = t2162*t3495
      t4928 = 8*t4884-96*t198*t27*t2182+64*t198*t27*t2187+128*t371*t594*
     #t4893-32*t198*t4471-96*t198*t27*t4898-8*t4213*t4332+16*t3403*t2196
     #*t251-8*t3403*t2196*t436-t4911+24*t4912-384*t2181*t4893+288*t2181*
     #t8*t4898-192*t2181*t8*t2197+96*t2181*t3460-384*t681*t8*t3884*t32
      t4936 = t402*t1712
      t4937 = t4936*t432
      t4940 = t138*t1714*t414
      t4942 = t4940*t2156*t145
      A1e0c9p1t1n = 24*(t4797+t3873+t2862+t2676+t3640+t4766+t3840+t4541+
     #t4733+t4501+t3807+t4694+t2630+t4464+t4410+t3768+t2943+t2585+t3725+
     #t2450+t4372+t3366+t4320+t2401+t2824+t3598+t4063+t4031+t4238+t3998+
     #t3959+t4928+t3917+t4883+t4152+t4104+t3180+t3138+t2782+t3105+t3063+
     #t2726+t4192+t2909+t3518+t4656+t3467+t3559+t4620+t3413+t2339+t3311+
     #t2541+t4584+t3220+t2283+t4838+t3260+t3024+t2226+t2493+t4275+t3688+
     #t2986)*t140*t4937*t4942
      A1e0c1p1t2n = A1e0c1p1t1n
      A1e0c2p1t2n = A1e0c2p1t1n
      A1e0c3p1t2n = A1e0c3p1t1n
      A1e0c4p1t2n = A1e0c4p1t1n
      A1e0c5p1t2n = A1e0c5p1t1n
      A1e0c6p1t2n = A1e0c6p1t1n
      A1e0c7p1t2n = A1e0c7p1t1n
      A1e0c8p1t2n = A1e0c8p1t1n
      A1e0c9p1t2n = A1e0c9p1t1n
      t4947 = 32*t57
      t4949 = 24*t131
      t4950 = 8*t129
      t4952 = -t121+16*t111-24*t113-32*t296+t47+t60+t122-t4947+24*t134+t
     #4949+t4950-32*t284-t126-t37-t77
      A2e0c1p1t1n = 96*Mzsq*t4952*t138*t422
      t4965 = t83*El
      t4969 = Eq*sqq*Er
      t4970 = 2*t4969
      t4971 = 4*t165
      t4972 = 8*t3
      A2e0c2p1t1n = 96*t427*(8*t340+4*t165*t58-2*t156*t197-4*t1415*SnThl
     #*t17*sll-8*t4965*t87+t4970-t4971+t202-t4972)*t138*t434
      t4987 = -t441+16*t445-32*t449+32*t458-128*t463+64*t483-32*t489+64*
     #t491+t497+t502-t505-t509+64*t515+32*t518+8*t521
      t4992 = t3*t26
      t4996 = t2575*sqq
      t5002 = t54*t915
      t5007 = -16*t524+t541-32*t544+16*t546-32*t549+64*t4992*t4283+32*t5
     #54+4*t1054*t4996+t593-t617+t620-t622+2*t54*t1368*t38-2*t5002+2*t54
     #*t898*t38
      t5009 = t80*t921
      t5011 = t80*t924
      t5013 = t80*t918
      t5024 = 64*t68*t447
      t5026 = 32*t68*t887
      t5028 = 32*t68*t443
      t5029 = t80*t443
      t5031 = t80*t887
      t5033 = -8*t5009-8*t5011+16*t5013-16*t80*t935+8*t80*t932+8*t80*t92
     #8+32*t624-t630-32*t634+t639+t5024-t5026-t5028+8*t5029+8*t5031
      t5034 = t80*t447
      t5036 = t91*t583
      t5050 = -16*t5034-32*t5036*t560+32*t669-32*t697-32*t712+64*t714-32
     #*t724+64*t726-32*t729-32*t763+t768-t772+64*t774-32*t778+16*t780
      t5053 = 16*t790
      t5054 = 4*t821
      t5056 = -t789+t5053-t792+t794+t796-t803-t804-t820-t5054+t825+t842-
     #t844-t848+t866-32*t869
      t5058 = 16*t889
      t5059 = 32*t910
      t5060 = 4*t938
      t5062 = 16*t883+t886-t888+t5058-t895-t897-t900-t903+t906-t5059-t50
     #60+12*t940+t943-t946-t950
      t5070 = t91*t2205
      t5087 = t205*t918
      t5089 = t1368*t54
      t5091 = t898*t54
      t5093 = 32*t4182
      t5094 = 4*t205*t2717-32*t5036*t564+64*t5036*t643-32*t5070*t3274-8*
     #t127*t829-8*t123*t838-4*t205*t935-2*t2625*t1282-8*t123*t928-8*t127
     #*t921+16*t127*t918+4*t5087+2*t5089+2*t5091-t5093
      t5096 = 32*t68*t427
      t5098 = 64*t68*t202
      t5099 = t54*t495
      t5101 = t80*t202
      t5103 = t80*t427
      t5105 = t80*t416
      t5110 = 4*t1055
      t5116 = -t5096+t5098-4*t5099-16*t5101+8*t5103+8*t5105+16*t962+16*t
     #1045-t1049+16*t1052+t5110-32*t1089+32*t1105-64*t1107+12*t1109-6*t1
     #111
      t5123 = 4*t1162
      t5131 = -32*t1140-16*t1144+64*t1149+t5123+t1167-32*t1191-32*t1193+
     #64*t1195-64*t1215+32*t1217+32*t1219-6*t1221-t1224+t1226+t1246
      t5142 = 64*t1253-16*t1256+64*t1260-128*t1262-32*t1270+64*t1272-32*
     #t1277+t1279+t1281-t1284+t1286+8*t1291+t1303-32*t1305+64*t1307
      t5145 = 16*t1314
      t5146 = 16*t1318
      t5147 = 32*t1320
      t5149 = 16*t1342
      t5151 = -6*t1309-t1312+t5145-t1317+t5146-t5147+t1323+t1327-t1329-t
     #1331-6*t1340+t5149-t1346+12*t1353-t1355
      t5152 = 32*t1364
      t5153 = 16*t1366
      t5162 = t15*t352
      t5169 = t1358+t1360-t5152+t5153-t1371+t1330+t1328-32*t1387+16*t139
     #1+16*t123*t849-32*t5070*t3277+64*t5070*t3280-64*t5162*t3313-32*t49
     #92*t4280-32*t4992*t3287
      t5185 = t3*t1804
      t5190 = 16*t868
      t5194 = -64*t5162*t1516-64*t5162*t1504+128*t15*t148*t99-32*t4992*t
     #341*t443-32*t4992*t341*t427+64*t5185*t1559-2*t2625*t506+t1401+t140
     #4+t5190-64*t1408+16*t1410+16*t1413+t1419+t1422
      t5195 = 16*t762
      t5196 = 32*t773
      t5199 = t80*t826
      t5201 = t80*t829
      t5203 = t80*t258
      t5207 = -t1427+t1442+t1446-t1450-t1474-t1475-t1476+t5195-t5196+t14
     #81-16*t80*t849-8*t5199-8*t5201+16*t5203+8*t80*t835
      t5213 = t447*t205
      t5215 = t54*t808
      t5217 = t54*t805
      t5227 = 8*t80*t838-4*t54*t4996-4*t5213-2*t5215+4*t5217+t1492-32*t1
     #493+t1500-64*t5162*t3553+128*t5162*t3316+32*t1552-64*t1563-32*t158
     #3+t1598-t1602
      t5239 = -8*t1614+8*t1617-16*t1619+t1625+16*t1628-t1634+16*t1636-16
     #*t1650-8*t1675+12*t1677+8*t1679-t1691+t1694-t1697+32*t1700-16*t170
     #2
      A2e0c3p1t1n = 96*(t4987+t5007+t5033+t5050+t5056+t5062+t5094+t5116+
     #t5131+t5142+t5151+t5169+t5194+t5207+t5227+t5239)*t140*t1715*t1720
      t5256 = -32*t1729+48*t1731+16*t1735-t1741-24*t1742+16*t1744-24*t17
     #46-t1752+t1754+24*t1758+64*t4489*t94
      t5264 = 64*t514
      t5265 = 16*t1843
      t5266 = -32*t91*Mzsq*t94-32*t91*t437*t94-48*t1761+t5264+t1781-t178
     #5-t1801-t5265-t1856+t1863+t1869
      t5268 = t57*t437
      t5269 = 32*t5268
      t5270 = 8*t1925
      t5271 = 8*t1930
      t5272 = 24*t1934
      t5273 = 48*t1939
      t5274 = 24*t1945
      t5275 = -t1885+t1895+t1910-t5269-t1923+t5270+t5271+t5272-t5273+t19
     #42+t5274
      t5277 = t50*CsThr
      t5287 = 32*t517
      t5288 = -t1948+t765-t1954+t863+t1956-t1959+24*t1962-t1965-32*t1992
     #*t5277*Mzsq+64*t1992*t5277*sll-32*t1992*t5277*t437-t5287
      A2e0c4p1t1n = 24*(t5256+t5266+t5275+t5288)*t138*t140*t1982
      t5295 = 2*t118*t427
      t5297 = 2*t118*t443
      t5312 = 2*t184
      t5314 = 8*t203
      t5315 = t5295+t5297-4*t3358*t1799-4*t3358*t1783-8*t4965*t2018-8*t4
     #965*t1390-2*t2012*t1989-2*t1875*t2436*t437+4*t160-t5312-t185-8*t20
     #0+t5314
      t5316 = 8*t1143
      t5317 = t1755*sll
      t5318 = 4*t5317
      t5319 = t2024*sll
      t5320 = 4*t5319
      t5322 = 8*t3*t437
      t5323 = 4*t125
      t5324 = t465*t1
      t5334 = t511+t444-t5316-t5318-t5320-t5322-t5323+t132+4*t250*t5324+
     #8*t1142*t1875+4*t250*t3664+8*t3*t1301+t135+8*t210
      A2e0c5p1t1n = 24*(t5315+t5334)*t138*t140*t2155*t218
      t5346 = 32*t2382
      t5347 = 16*t2386
      t5348 = 128*t2160+32*t2166+16*t2251+128*t2308-t2325-t2327-16*t2328
     #+32*t2331+t2338-t2364+t2370+t2372-t2375+t2378+t2381+t5346-t5347
      t5349 = 16*t2390
      t5353 = 64*t2419
      t5354 = 48*t2427
      t5355 = 32*t2429
      t5356 = 48*t2434
      t5357 = 32*t2461
      t5358 = 16*t2465
      t5359 = -t5349+t2394-t2398+6*t2399+6*t2412-12*t2414+t5353-t2424-t5
     #354+t5355+t2432-t5356+t2440-t2443+t5357-t2464-t5358
      t5361 = 16*t2469
      t5362 = 12*t2474
      t5367 = 48*t2535
      t5368 = 48*t2537
      t5369 = t2468-t5361+t2472+t5362-12*t2486+18*t2491-128*t2504+64*t25
     #07-t2514+t2519+t2525+t2530-t2534-t5367-t5368+t2540-t2544
      t5370 = 32*t2545
      t5372 = 12*t2573
      t5382 = t5370-t2549+t2552-12*t2565+t5372+18*t2579-t2592+t2596-t259
     #9-t2603+t2605+t2608+t2610+t2614+8*t2134*t2807-16*t1914*t2991+8*t19
     #14*t2821+8*t2134*t2821
      t5391 = 12*t2645
      t5393 = t80*t1356
      t5399 = -16*t2840*t2717-16*t1914*t3097+128*t2621-t2624+32*t2626+t2
     #638-t2641+t2644-t5391-24*t2651+32*t5393+t2690-t2696+24*t2711+32*t2
     #714-64*t2718-32*t2722
      t5412 = 128*t2724+16*t2731+16*t2743-t2753+32*t2756-16*t2761+32*t27
     #64-16*t2767+t2777+t2787-64*t2789+32*t2792+t2796-t2798+t2801-64*t28
     #05+32*t2808-64*t2819
      t5425 = 32*t2822+t2827-t2829-t2831+32*t2853-t2878+t2881+32*t2882-t
     #2886-8*t2915-12*t2917+6*t2919+32*t2921-32*t2923+32*t2931-8*t2933+6
     #*t2935
      t5443 = -32*t2937+16*t2939-24*t2959+16*t2962+20*t2964+48*t2977+96*
     #t2982+48*t2984-32*t2987+96*t2989-64*t2992+t3000-64*t3001-64*t3038+
     #96*t3040-96*t3042+16*t3069-24*t3079
      t5462 = 48*t3085+96*t3087-64*t3098-32*t3101-192*t3103-192*t3108-19
     #2*t3110+128*t3112-t3141-t3148-128*t3156+64*t3159-64*t3161-64*t3163
     #-64*t3169+64*t3171+96*t3173
      t5466 = 4*t3187
      t5467 = 4*t3199
      t5474 = Er*t4
      t5475 = t362*t5474
      t5480 = t2418*t2355
      t5493 = -96*t3178-64*t3182-4*t3185-t5466-t5467+18*t3204-12*t3206+6
     #4*t3208+4*t3223+t3234+32*t3235-32*t3237+96*t5475*t8*t602*t141-128*
     #t5480*t615+96*t1364*t965+32*t5185*t3460+64*t177*t2395*t2087+64*t17
     #7*t2361*t609
      t5500 = t2305*t5474
      t5503 = t91*t251
      t5523 = t73*t2395
      t5528 = t73*t2361
      t5538 = 32*t5185*t2192+192*t2418*t1988*t2087+96*t5500*t706+96*t550
     #3*t2207-64*t91*t2248*El*t984*SnThr*CsThl*Mzsq+32*t5503*t2233+32*t5
     #503*t2219-64*t5503*t2212-64*t5475*t8*t2196*sll-128*t5523*t149+64*t
     #5523*t2087+64*t5528*t609+192*t73*t2350*t3530-128*t5528*t149+32*t54
     #75*t2192-t3244+t3246
      t5553 = t91*t2500
      t5556 = t3249-t3252+64*t3254-t3259+t3266-48*t3269-t3276-t3279+t328
     #2+t3285-8*t3298-256*t73*t2121*El*t568-64*t5500*t612+192*t5480*t571
     #-28*t3300+8*t3302-8*t3306+96*t5553*t3250
      t5579 = t2681*t826
      t5581 = t2663*t826
      t5595 = t1914*t1655
      t5597 = -128*t4489*t583*t981-64*t5503*t2229+96*t5553*t4227+96*t550
     #3*t2223-128*t5553*t976-64*t91*t119*t2312-128*t362*t5474*t6*t3556-2
     #4*t2681*t838+24*t5579+6*t5581-6*t2663*t838+24*t2681*t2981-16*t2681
     #*t2717-4*t2625*t2991+2*t2625*t2821+24*t1914*t2710-32*t5595
      t5615 = Er*t141
      t5617 = 128*t2256*t5615
      t5621 = 32*t1914*t3230+8*t1914*t2807+6*t2663*t2981-4*t2663*t2717+2
     #4*t2681*t3920-4*t2625*t3097+6*t2663*t3920+24*t1914*t1466-48*t3309+
     #t3315-t3318-t3322+t3325-t5617-96*t3375+64*t3378-16*t3383+t3394
      t5634 = 32*t3567
      t5635 = 96*t3396-64*t3399+32*t3401+t3409+32*t3454+64*t3465+32*t347
     #2-16*t3476+t3484-32*t3546+64*t3548-32*t3550+64*t3557+t3566+t5634-t
     #3583+t3585
      t5636 = 64*t3594
      t5637 = t3587-t3590+t3593+t5636-t3597-t3608+t3609+t3610+t3613-t361
     #5+t3618-t3620-t3634+t3636+t3649+t3651+t3666-t3669
      t5648 = t1914*t798
      t5654 = t2625*t2654
      t5656 = t1914*t2654
      t5662 = -t3683-t3687+256*t3689-t3692+t3694+t3696+t3705-t3708-128*t
     #2305*t5615*t342-24*t1914*t846+24*t5648+4*t2625*t2660-8*t2663*t260-
     #4*t5654-16*t5656+16*t1914*t2660+4*t2663*t935
      t5663 = t2663*t924
      t5672 = 128*t2132*t907
      t5677 = t2681*t924
      t5679 = t2681*t918
      t5681 = t2625*t3569
      t5683 = t2663*t918
      t5685 = t2625*t3572
      t5691 = t1914*t912
      t5693 = t2625*t3579
      t5695 = t2134*t3579
      t5697 = 6*t5663+4*t2625*t4683-2*t2625*t4686-6*t2663*t928-t5672+16*
     #t2681*t935-24*t2681*t928+24*t5677-16*t5679-4*t5681-4*t5683+2*t5685
     #-24*t1914*t944-2*t2625*t4691+24*t5691+2*t5693+8*t5695
      t5703 = t2134*t3572
      t5705 = t1914*t3579
      t5709 = t2840*t918
      t5717 = 48*t3797
      t5718 = 16*t3801
      t5719 = -8*t2134*t4691-8*t2134*t4686+8*t5703+8*t5705+16*t2840*t935
     #-16*t5709-8*t1914*t4691-8*t1914*t4686+32*t5475*t3460-t3737+t3740-t
     #3743-t3746-t5717+t3800-t5718+t3804
      t5720 = t1914*t3569
      t5722 = t1914*t3572
      t5724 = t1914*t448
      t5744 = t91*t110
      t5755 = t3815-t3816-t3817+t3819-t3822-16*t5720+8*t5722+32*t5724+16
     #*t1914*t4683-32*t2681*t260-64*t5475*t8*t2186*sll+96*t5475*t8*t596*
     #t141+192*t73*t2121*t3530-256*t2418*t5615*t94+32*t5744*t2233+32*t57
     #44*t2219-128*t2310*t2121*t2087-64*t2433*t5474*t612
      t5759 = 48*t3838
      t5760 = 32*t3841
      t5761 = 48*t3845
      t5762 = 64*t3849
      t5763 = 32*t3863
      t5764 = 16*t3865
      t5765 = 48*t3871
      t5766 = 32*t3875
      t5767 = -t3827-24*t3828+t3835-t3837-t5759+t5760+t3844-t5761+t5762-
     #t3852+t5763-t5764-t3868+t3870-t5765+t5766+t3878
      t5782 = -t3880+64*t3886-t3890+128*t3895+48*t3918+96*t3921+128*t392
     #3-192*t3925+64*t3927+96*t3929+16*t3934+18*t3943-12*t3946-96*t3948-
     #24*t3966+24*t3971-t3997+t4001
      t5793 = -t4011-t4018+128*t4021-t4023-t4024+t4026-32*t4033+96*t4035
     #-24*t4037-64*t4041-96*t4043+64*t4045-t4077+t4079+t4081-64*t4086+32
     #*t4089
      t5808 = -t4113-t4117-32*t1949*t1655+t4119-t4128+16*t4142+t4148-32*
     #t4150-20*t4154+32*t4156-128*t4175-128*t4183-64*t4186+256*t4188-48*
     #t4202+32*t4206-48*t4208-48*t4210
      t5814 = t80*t2405
      t5816 = t80*t2409
      t5819 = t54*t898*t2163
      t5821 = t141*t2195
      t5823 = t54*t5821*Mzsq
      t5825 = t2365*t437
      t5827 = t54*t5825*Mzsq
      t5829 = t2365*sqq
      t5831 = t54*t5829*t437
      t5834 = t54*t5821*sqq
      t5837 = 32*t2264*t2458
      t5839 = 64*t2979*t2455
      t5840 = t80*t2176
      t5842 = t80*t904
      t5845 = 32*t3489*t2447
      t5846 = t80*t899
      t5848 = t80*t496
      t5850 = t80*t801
      t5852 = t80*t2569
      t5854 = 2*t2625*t2807-8*t5814-8*t5816-2*t5819-2*t5823+4*t5827-6*t5
     #831+4*t5834+t5837-t5839+16*t5840+16*t5842+t5845-24*t5846-24*t5848-
     #24*t5850+16*t5852
      t5857 = 96*t2264*t5615*Mzsq
      t5860 = 64*t2264*t1988*t2163
      t5871 = t5857-t5860-t4222-t4224-t4229-128*t4251+t4258-96*t4262-96*
     #t4269-96*t4271+64*t4273-16*t4286-48*t4294-48*t4296+32*t4298-t4301+
     #t4323+96*t4327
      t5876 = t54*t2365*t32
      t5883 = 64*t4391
      t5884 = -128*t4329-64*t4333+8*t5876+32*t4335-128*t4338+96*t4342+96
     #*t4344-64*t4346-t4385-t4388+t4390+t5883-t4393+t4395-t4397-t4399-t4
     #400
      t5891 = t54*t495*t2163
      t5894 = t54*t5829*Mzsq
      t5898 = t54*t2365*t40*t32
      t5901 = 96*t2256*t1989
      t5904 = 64*t2264*t2395*sll
      t5906 = 32*t3489*t2458
      t5908 = 96*t2256*t2142
      t5910 = 64*t2256*t2455
      t5911 = t80*t2482
      t5913 = t80*t901
      t5915 = t80*t2488
      t5917 = t80*t2479
      t5921 = 96*t2264*t2121*t141
      t5923 = 32*t2264*t2447
      t5926 = t80*t818
      t5928 = -16*t4467+8*t54*t2365*t38*t32+4*t5891-6*t5894-8*t5898+t590
     #1-t5904+t5906+t5908-t5910-8*t5911-24*t5913+16*t5915-8*t5917+t5921+
     #t5923+32*t80*t2576-32*t5926
      t5939 = 16*t4631
      t5940 = t4536-t4540-t4547+t4550+96*t4570-48*t4573+96*t4578-64*t458
     #0-48*t4582+32*t4585-64*t4587-128*t4589+t4605-t4608+t4611-t4630-t59
     #39
      t5941 = 16*t4637
      t5943 = 12*t4663
      t5944 = 4*t4676
      t5945 = 4*t4681
      t5949 = -t5941-t4644+t4647+t4649-20*t4650+t4653-t5943-t4668+t4675+
     #t5944-t4680+t5945-t4704-t4712+t4715+96*t4742-64*t4744-48*t4746
      t5960 = 32*t4748-64*t4754+32*t4756+32*t4758-16*t4760+96*t4762+32*t
     #4764-16*t4769-48*t4771+t4773-t4775+t4776-t4778-t4782+t4784-t4790+t
     #4792
      t5972 = t4796-48*t4800+32*t4804+64*t4808+8*t4813+t4822-t4824-t4826
     #+t4846-32*t4853-32*t4855-16*t4858-t4862+32*t4873+32*t4879-16*t4884
     #+t4911-16*t4912
      A2e0c6p1t1n = -24*(t5755+t5493+t5621+t5949+t5940+t5928+t5884+t5871
     #+t5854+t5808+t5793+t5782+t5767+t5597+t5556+t5538+t5462+t5443+t5972
     #+t5960+t5425+t5412+t5399+t5382+t5369+t5719+t5697+t5662+t5637+t5359
     #+t5348+t5635)*t140*t4937*t4942
      A2e0c1p1t2n = A2e0c1p1t1n
      A2e0c2p1t2n = A2e0c2p1t1n
      A2e0c3p1t2n = A2e0c3p1t1n
      A2e0c4p1t2n = A2e0c4p1t1n
      A2e0c5p1t2n = A2e0c5p1t1n
      A2e0c6p1t2n = A2e0c6p1t1n
      t5983 = (-sqq+sllmin)*t404
      t5986 = log(DCMPLX(-t5983/dels))
      t5987 = t5986*Eq
      t5988 = t5987*Er
      t5990 = -j2+4*t5988
      t5994 = t58*t357
      t5997 = 8*t14
      t5999 = t86*t98*t1
      t6001 = t14*SnThl
      t6005 = 4*t357
      t6006 = 8*t442-16*t386*t87+4*t5994-16*t386+t5997+t32-4*t5999+16*t6
     #001*t17+8*t1556-t6005
      t6008 = t5990*t6006*t140
      t6010 = 1/j2
      t6011 = t6010*t138
      A1e0c1p1t1n1 = -24*t6008*t145*t404*t6011
      A1e0c2p1t1n1 = -24*t6008*t2156*t6011
      t6018 = t233*SnThl
      t6021 = t233*t6
      t6033 = t28*t58
      t6036 = t1*SnThl
      t6038 = t9*t32
      t6039 = t6038*t58
      t6040 = t6036*SnThr*t6039
      t6042 = t386*t99
      t6043 = 64*t6042
      t6044 = -96*t6018*t17-96*t6021*t341-96*t333*t87+192*t369*t87+96*t3
     #69*t11-t234+16*t442*t99-8*t1555*t7*t6033+4*t6040-t6043+t286+t287-t
     #293
      t6048 = t6001*t23
      t6050 = t1555*t957
      t6052 = t6038*t1
      t6053 = t86*t6052
      t6054 = 8*t6053
      t6055 = -t309-t313-32*t233*t275*t3328-t322+t323+48*t6048+24*t6050-
     #t6054+t122-t334-t336+t370+t375
      t6058 = t215*t138
      t6060 = t6010*t140
      A1e0c3p1t1n1 = -48*t5990*(t6044+t6055)*t6058*t6060*t407
      t6067 = t5990*sqq*Mzsq*t6006*t138
      t6068 = t429*t432
      t6070 = t6060*t6068*t145
      A1e0c4p1t1n1 = 24*t6067*t6070
      t6073 = t6060*t6068*t218
      A1e0c5p1t1n1 = 24*t6067*t6073
      t6075 = t14*j2
      t6076 = t61*t6075
      t6081 = t5986*t2
      t6084 = t5986*t67
      t6087 = j2*CsThl
      t6088 = t233*CsThr
      t6091 = j2*t38
      t6092 = t233*t40
      t6095 = j2*CsThr
      t6096 = t122*CsThl
      t6098 = j2*t233
      t6102 = CsThr*SnThl
      t6111 = t14*t2
      t6112 = t5986*t6111
      t6113 = t6112*t196
      t6115 = CsThr*t6*t341
      t6118 = t58*t6075
      t6124 = -16*t6076*t99+4*t5987*t2425-32*t6081*t2367+96*t6084*t1756-
     #96*t6087*t6088+96*t6091*t6092+t6095*t6096-32*t6098*t254-256*t6084*
     #t196*t6102*t23-32*t6084*Er*t6*t341*t110+384*t6113*t6115+64*t6118*t
     #99+8*t6075*t8*t6033
      t6125 = t6084*Er
      t6128 = t6081*Er
      t6139 = t6112*Er
      t6150 = j2*SnThl*SnThr
      t6158 = 160*t6125*t115+48*t6128*t119-4*t5988*t336-16*t6128*t260-32
     #*t6125*t255-224*t6125*t110-384*t6139*t87-384*t6139*t11-128*t6139*t
     #279-48*t6075*SnThl*t23+8*t6150*t6052-24*t6075*t6*t957+192*t6125*t9
     #9
      t6167 = t6098*CsThl
      t6168 = t6102*t17
      t6179 = j2*t1
      t6194 = -32*t6128*t1551+96*t6125*t1559+96*t6098*t40*t393-192*t6167
     #*t6168-96*t6167*t6115+64*t6084*t1878*t40*SnThl*t23-j2*t122-4*t6179
     #*t86*t6039-384*t6112*t190*t393+768*t6113*t6168+96*t6098*t11+32*t60
     #98*t279+96*t6098*t87
      t6204 = t32*j2
      t6225 = 384*t6112*t1875-384*t6112*t1879+128*t6112*t4138*t240-40*t6
     #075*t115-12*t6204*t1433+4*t6091*t320*t1+8*t6075*t255+56*t6087*t188
     #2+16*t6081*Er*SnThl*t17*t119+32*t6098-128*t6139+8*t6204*t1-24*t607
     #5*sqq
      t6229 = t138*t432
      t6233 = t140*t218*t145
      t6234 = t429*t6010*t6233
      A1e0c6p1t1n1 = 48*Mzsq*(t6124+t6158+t6194+t6225)*t6229*t6234
      t6236 = t58*t1313
      t6238 = t61*t1313
      t6240 = t58*t867
      t6242 = t61*t867
      t6244 = t558*t1
      t6245 = t980*t6244
      t6248 = t1755*SnThl*t23
      t6251 = t386*t1390
      t6253 = t385*Mzsq
      t6254 = t6253*t6168
      t6256 = 8*t867
      t6257 = 8*t1313
      t6261 = -16*t6236+8*t6238-16*t6240+8*t6242-4*t6245-4*t6248+32*t604
     #2-16*t6251-16*t6254+t6256+t6257+t665-16*t42-8*t46+32*t74
      t6262 = 16*t73
      t6265 = t427*t1
      t6266 = t58*t6265
      t6268 = t558*t14
      t6269 = t86*t6268
      t6271 = t1555*t971
      t6273 = t1313*t87
      t6275 = t1313*t11
      t6277 = t2024*sqq
      t6278 = t58*t6277
      t6280 = 2*t122
      t6281 = 4*t6265
      t6282 = 4*t6277
      t6283 = t322-t6262+t438-32*t6048-16*t6050+t6054+4*t6266+16*t6269+8
     #*t6271+16*t6273+8*t6275+4*t6278-t6280-t6281-t6282
      t6285 = (t6261+t6283)*t5990
      t6286 = t140*t1712
      A1e0c7p1t1n1 = 24*t6285*t6286*t6229*t6010*t145
      A1e0c8p1t1n1 = 24*t6285*t2155*t406*t6010*t218
      t6305 = t38*t122
      t6335 = -80*t61*t2418+112*t58*t2418-24*t58*t2597+16*t14*t240*t235*
     #t32-t2958+8*t105*t6305+192*t61*t2253+32*t254*t747+96*t58*t747-96*t
     #61*t747-64*t254*t2253-192*t58*t2253-t3078+96*t58*t852-96*t61*t852+
     #32*t254*t852-8*t665*SnThl*t17*t1+48*t1313*SnThl*t23+48*t6268*t980
      t6373 = 24*t867*t10*t8*sqq+24*t1957*t11-64*t6240*t99-8*t6244*t695+
     #40*t1857+12*t1862-56*t1883+40*t1892+12*t1894-56*t1897+128*t386*t15
     #51+16*t386*t660*t1431-8*t1433*t4807-32*t442*t1551-192*t369*t1559-9
     #6*t333*t1390+192*t369*t1390+96*t369*t1531+192*t333*t99-384*t369*t9
     #9
      t6381 = 48*t2418
      t6382 = 16*t2597
      t6383 = 32*t747
      t6384 = 32*t852
      t6385 = 64*t2253
      t6386 = t58*t1755
      t6397 = t14*t10
      t6401 = t14*t9
      t6411 = 192*t369*t2018+96*t369*t1763-96*t333*t2018-t6381+t6382-t63
     #83-t6384+t6385+4*t6386*t1551+16*t6238*t99+t765+t863-4*t170*t912-8*
     #t310*t1204-8*t254*t1957-48*t6397*t8*t32-96*t6401*t695+16*t1*t9*t23
     #42+64*t233*t278*t277*sqq
      t6430 = 24*t1957
      t6431 = 8*t439
      t6432 = 24*t1891
      t6435 = 2*t1324
      t6436 = 8*t1964
      t6447 = 192*t6021*t957-96*t852*t87-32*t852*t279-96*t852*t11+192*t6
     #018*t23-96*t747*t87-32*t747*t279-96*t747*t11-4*t61*t1964+t6430-t64
     #31+t6432+2*t58*t1324-t6435-t6436-64*t6236*t99-8*t6236*t1559+4*t279
     #4*t1551-8*t6240*t1559+16*t6242*t99
      A1e0c9p1t1n1 = 48*t5990*(t6335+t6373+t6411+t6447)*t2155*t404*t6011
     #*t6233
      A1e0c1p1t2n1 = A1e0c1p1t1n1
      A1e0c2p1t2n1 = A1e0c2p1t1n1
      A1e0c3p1t2n1 = A1e0c3p1t1n1
      A1e0c4p1t2n1 = A1e0c4p1t1n1
      A1e0c5p1t2n1 = A1e0c5p1t1n1
      A1e0c6p1t2n1 = A1e0c6p1t1n1
      A1e0c7p1t2n1 = A1e0c7p1t1n1
      A1e0c8p1t2n1 = A1e0c8p1t1n1
      A1e0c9p1t2n1 = A1e0c9p1t1n1
      t6456 = 4*t1
      t6463 = t5990*Mzsq*t32*(sqq-t6456+4*t1433-4*t6036*t17)*t138
      A2e0c1p1t1n1 = 96*t6463*t6070
      A2e0c2p1t1n1 = 96*t6463*t6073
      t6483 = t5986*t32
      t6504 = 16*t6081*t86*t98*t1875-128*t6084*t58*t86*t9*Er-4*t99*t58*t
     #6179-16*t6076+32*t6118+t58*t6204-32*t6081*t1756+4*t6483*t118-16*t6
     #081*t38*t2044+48*t6081*sqq*t1875-4*t5986*CsThr*t32*t1739+64*t6084*
     #t6*t341*Er+128*t6084*SnThl*t17*Er
      t6506 = sqq*j2*t1
      t6533 = 8*t87*t6506-32*t6081*SnThl*t17*t1756+32*t386*t585*SnThr*j2
     #-16*t6075-t6204-12*t58*t6506+4*t61*t6506-128*t6084*t1875+64*t6084*
     #t1879-32*t6401*t6150-16*t6397*t8*j2+8*t6506+64*t6125
      A2e0c3p1t1n1 = 192*t427*(t6504+t6533)*t138*t432*t6234
      t6540 = 2*t32
      t6546 = 4*t2024
      t6547 = 4*t1755
      t6550 = 8*t357
      t6551 = -8*t5994+t443+t427-t6540+8*t5999-4*t1755*t87-4*t2024*t87-t
     #6546-t6547+4*t6386+4*t2794+t6550
      t6552 = sqq*t6551
      A2e0c4p1t1n1 = 24*t6552*t5990*t140*t2155*t6011*t145
      t6560 = t432*t140
      A2e0c5p1t1n1 = 24*t6552*t5990*t1712*t6560*t6011*t218
      t6579 = 16*t867
      t6580 = 16*t1313
      t6582 = -32*t6236+16*t6238-32*t6240+16*t6242-t709-t2976+4*t2794*t9
     #9-8*t6245-8*t6248+4*t6386*t99-4*t170*t924-8*t6040+t6043-32*t6251-3
     #2*t6254+t6579+t6580+t665-32*t42
      t6586 = 16*t76
      t6587 = 32*t73
      t6600 = 8*t6265
      t6601 = 8*t6277
      t6602 = -24*t46+64*t74+8*t308+t6586-t6587+t438-64*t6048-32*t6050+1
     #6*t6053+12*t6266+32*t6269+16*t6271+32*t6273+16*t6275+12*t6278-4*t6
     #1*t6265-t6280+2*t336-t6600-t6601
      A2e0c6p1t1n1 = 48*t5990*(t6582+t6602)*t6060*t6229*t1712*t218*t145
      A2e0c1p1t2n1 = A2e0c1p1t1n1
      A2e0c2p1t2n1 = A2e0c2p1t1n1
      A2e0c3p1t2n1 = A2e0c3p1t1n1
      A2e0c4p1t2n1 = A2e0c4p1t1n1
      A2e0c5p1t2n1 = A2e0c5p1t1n1
      A2e0c6p1t2n1 = A2e0c6p1t1n1
      t6610 = t158*t1
      t6611 = 4*t6610
      t6612 = 16*t385
      t6613 = 8*t39
      t6615 = j2*j4
      t6616 = log(DCMPLX(-t5983))
      t6617 = t6616**2
      t6619 = 2*t118*t6617
      t6620 = t6617*CsThr
      t6622 = 2*t118*t6620
      t6623 = log(DCMPLX(dels))
      t6624 = t6623**2
      t6626 = 2*t118*t6624
      t6627 = t6624*CsThr
      t6629 = 2*t118*t6627
      t6630 = log(DCMPLX(delc))
      t6631 = t5986*t6630
      t6633 = 4*t118*t6631
      t6634 = t6631*CsThr
      t6636 = 4*t118*t6634
      t6637 = Er*j4
      t6639 = 4*t5987*t6637
      t6644 = log(DCMPLX(Er/delc/Eq))
      t6645 = t6644*j2
      t6646 = t6645*CsThr
      t6647 = t6615+t6619-t6622-t6626+t6629-t6633+t6636-t6639-t6645+t664
     #6
      t6648 = (-t6005+t5997+t32+t6611-t6612+t6613)*t6647
      t6649 = t140*t404
      t6652 = 1/(-1+CsThr)
      t6653 = t138*t6652
      t6654 = 1/j4
      t6655 = t6010*t6654
      A1e0c1p1t1n2 = -12*t6648*t6649*t6653*t6655
      A1e0c2p1t1n2 = 0
      t6659 = t233*CsThl
      t6660 = 64*t6659
      t6661 = t76*CsThl
      t6662 = 4*t6661
      t6663 = t39*sqq
      t6664 = 8*t6663
      t6665 = t233*t38
      t6666 = 32*t6665
      t6667 = t385*sqq
      t6668 = 32*t6667
      t6672 = 1/(1+CsThl)
      t6675 = t6652*t6010
      t6676 = t6675*t6654
      A1e0c3p1t1n2 = -24*t6647*(-t323+t234-t6660+t322-t6662-t6664+t6666-
     #t122+t6668)*t6672*t215*t406*t6676
      t6682 = t6068*t138
      t6683 = t6682*t6676
      A1e0c4p1t1n2 = 12*t427*t6648*t140*t6683
      A1e0c5p1t1n2 = 0
      t6686 = 24*t6615*t73
      t6688 = 64*t6615*t6659
      t6690 = 8*t6615*t76
      t6692 = 32*t6615*t6665
      t6694 = 24*t6645*t73
      t6695 = t6645*t6659
      t6696 = 64*t6695
      t6698 = 8*t6645*t76
      t6700 = 32*t6645*t6665
      t6701 = t6645*t6088
      t6702 = 32*t6701
      t6703 = t6645*t335
      t6706 = 48*t68*t6624*sqq
      t6707 = -t6686-t6688+t6690+t6692+t6694+t6696-t6698-t6700+t6702-t67
     #03+t6706
      t6708 = t6111*Er
      t6709 = t6624*CsThl
      t6710 = t6708*t6709
      t6711 = 128*t6710
      t6712 = t6624*t32
      t6714 = 16*t3*t6712
      t6715 = t6624*t38
      t6717 = 64*t6708*t6715
      t6720 = 2*t118*t6624*t122
      t6721 = t6708*t6627
      t6722 = 64*t6721
      t6725 = 48*t68*t6617*sqq
      t6726 = t6617*CsThl
      t6727 = t6708*t6726
      t6728 = 128*t6727
      t6729 = t6617*t32
      t6731 = 16*t3*t6729
      t6732 = t6617*t38
      t6734 = 64*t6708*t6732
      t6737 = 2*t118*t6617*t122
      t6738 = t6708*t6620
      t6739 = 64*t6738
      t6740 = t6711-t6714-t6717+t6720+t6722-t6725-t6728+t6731+t6734-t673
     #7-t6739
      t6743 = 128*t6112*t6637
      t6745 = 128*t6708*t6631
      t6748 = t6128*t6630*t32*CsThl
      t6749 = 16*t6748
      t6753 = 32*t6125*t6630*t38*sqq
      t6754 = t6630*CsThl
      t6756 = t6125*t6754*sqq
      t6757 = 128*t6756
      t6758 = t6630*CsThr
      t6760 = t6125*t6758*sqq
      t6761 = 96*t6760
      t6764 = 256*t6139*t6758*CsThl
      t6766 = t6128*t6758*t32
      t6767 = 32*t6766
      t6770 = 128*t6139*t6758*t38
      t6773 = 4*t5988*t6758*t122
      t6777 = 16*t6128*j4*t32*CsThl
      t6778 = -t6743-t6745+t6749+t6753-t6757-t6761-t6764+t6767+t6770-t67
     #73+t6777
      t6782 = 32*t6125*j4*t38*sqq
      t6786 = 128*t6125*j4*CsThl*sqq
      t6787 = t3*t6624
      t6789 = 8*t6787*t119
      t6790 = t68*t6624
      t6792 = CsThr*t38*sqq
      t6794 = 16*t6790*t6792
      t6796 = 64*t6790*t110
      t6797 = t68*t6617
      t6799 = 64*t6797*t110
      t6800 = t3*t6617
      t6802 = 8*t6800*t119
      t6804 = 16*t6797*t6792
      t6806 = 4*t6646*t6661
      t6808 = 8*t6646*t6663
      t6810 = 32*t6646*t6667
      t6812 = t3*t6729*CsThl
      t6813 = 8*t6812
      t6814 = t6782-t6786-t6789-t6794+t6796-t6799+t6802+t6804-t6806-t680
     #8+t6810-t6813
      t6819 = 16*t68*t6732*sqq
      t6821 = t68*t6726*sqq
      t6822 = 64*t6821
      t6824 = t68*t6620*sqq
      t6825 = 48*t6824
      t6828 = 128*t6708*t6620*CsThl
      t6830 = t3*t6620*t32
      t6831 = 16*t6830
      t6834 = 64*t6708*t6620*t38
      t6837 = 2*t118*t6620*t122
      t6838 = t6645*t6661
      t6839 = 4*t6838
      t6841 = 8*t6645*t6663
      t6842 = t6645*t6667
      t6843 = 32*t6842
      t6844 = t6645*t1882
      t6845 = 24*t6844
      t6846 = -t6819+t6822+t6825+t6828-t6831-t6834+t6837+t6839+t6841-t68
     #43-t6845
      t6848 = 64*t6645*t369
      t6849 = t76*CsThr
      t6850 = t6645*t6849
      t6851 = 8*t6850
      t6854 = 32*t6645*t6088*t38
      t6856 = 4*t6615*t6661
      t6858 = 8*t6615*t6663
      t6860 = 32*t6615*t6667
      t6863 = 4*t5987*t6637*t122
      t6866 = 96*t6084*t6637*sqq
      t6869 = 256*t6112*t6637*CsThl
      t6872 = 32*t6081*t6637*t32
      t6875 = 128*t6112*t6637*t38
      t6876 = -t6848+t6851+t6854-t6856-t6858+t6860+t6863+t6866+t6869-t68
     #72-t6875
      t6880 = 96*t68*t6631*sqq
      t6881 = t6631*CsThl
      t6882 = t6708*t6881
      t6883 = 256*t6882
      t6886 = 32*t3*t6631*t32
      t6889 = 128*t6708*t6631*t38
      t6892 = 4*t118*t6631*t122
      t6893 = t6708*t6634
      t6894 = 128*t6893
      t6896 = t3*t6712*CsThl
      t6897 = 8*t6896
      t6900 = 16*t68*t6715*sqq
      t6902 = t68*t6709*sqq
      t6903 = 64*t6902
      t6905 = t68*t6627*sqq
      t6906 = 48*t6905
      t6909 = 128*t6708*t6627*CsThl
      t6910 = t6880+t6883-t6886-t6889+t6892+t6894+t6897+t6900-t6903-t690
     #6-t6909
      t6912 = t3*t6627*t32
      t6913 = 16*t6912
      t6916 = 64*t6708*t6627*t38
      t6919 = 2*t118*t6627*t122
      t6922 = 32*t6125*t6758*t64
      t6925 = 16*t6128*t6758*t44
      t6928 = 128*t6125*t6758*t158
      t6930 = 32*t6615*t233
      t6931 = t6615*t122
      t6933 = 64*t6708*t6624
      t6935 = 64*t6708*t6617
      t6937 = 32*t6645*t233
      t6938 = t6645*t122
      t6939 = t6913+t6916-t6919-t6922-t6925+t6928+t6930-t6931-t6933+t693
     #5-t6937+t6938
      t6944 = t6672*t6010
      t6948 = t429*t6654*t140
      A1e0c6p1t1n2 = 24*Mzsq*(t6707+t6740+t6778+t6814+t6846+t6876+t6910+
     #t6939)*t6944*t6652*t6229*t6948
      t6951 = 8*t6661
      t6952 = 16*t6663
      t6954 = 8*t39*Mzsq
      t6955 = 16*t6253
      t6957 = 4*t427*t1432
      t6959 = 4*t2024*t158
      t6961 = 8*t39*t437
      t6963 = 16*t385*t437
      t6964 = t6256-t6951-t6952+t6668-t6262+t6257-t6282-t6280-t6281+t438
     #+t6954-t6955+t665+t6957+t6959+t6961-t6963+t322
      t6968 = t6229*t1712
      t6969 = t6675*t6968
      A1e0c7p1t1n2 = 12*t6964*t6647*t140*t6654*t6969
      A1e0c8p1t1n2 = 0
      t6971 = t1313*t158
      t6972 = 32*t6971
      t6974 = 4*t1755*t44
      t6976 = 8*t1313*t64
      t6978 = 8*t39*t443
      t6980 = 4*t1432*t438
      t6981 = t385*t443
      t6982 = 32*t6981
      t6983 = t6972-t6974-t6976-t6978-t6980+t6982+t6381-t6382+t6383+t638
     #4-t6385-t765-t863
      t6985 = 8*t2597*CsThl
      t6986 = 64*t3155
      t6987 = t6665*sqq
      t6988 = 64*t6987
      t6989 = 128*t3152
      t6991 = 32*t6665*Mzsq
      t6993 = 64*t747*CsThl
      t6994 = 16*t3158
      t6996 = 32*t6665*t437
      t6998 = 64*t852*CsThl
      t6999 = -t6430+t6431-t6432+t6985-t6986-t6988+t6989+t6991-t6993+t69
     #94+t6996-t6998+t6435+t6436
      A1e0c9p1t1n2 = 24*(t6983+t6999)*t6647*t404*t6672*t1712*t6560*t138*
     #t6676
      A1e0c1p1t2n2 = A1e0c1p1t1n2
      A1e0c2p1t2n2 = 0
      A1e0c3p1t2n2 = A1e0c3p1t1n2
      A1e0c4p1t2n2 = A1e0c4p1t1n2
      A1e0c5p1t2n2 = 0
      A1e0c6p1t2n2 = A1e0c6p1t1n2
      A1e0c7p1t2n2 = A1e0c7p1t1n2
      A1e0c8p1t2n2 = 0
      A1e0c9p1t2n2 = A1e0c9p1t1n2
      t7008 = 4*t1432
      A2e0c1p1t1n2 = 48*t665*(sqq-t6456+t7008)*t6647*t140*t6683
      A2e0c2p1t1n2 = 0
      t7014 = t158*t6800
      t7015 = 8*t7014
      t7016 = Er*t6617
      t7017 = t7016*CsThr
      t7018 = t362*t7017
      t7019 = 16*t7018
      t7021 = 2*t351*t7017
      t7023 = 16*t385*t6646
      t7026 = t158*t1*t6644*j2
      t7027 = 4*t7026
      t7028 = t357*t6646
      t7029 = 8*t7028
      t7032 = 4*t158*t6179*j4
      t7033 = t5986*Er
      t7034 = t7033*j4
      t7036 = 64*t460*t7034
      t7038 = 32*t362*t7034
      t7041 = 4*t6483*t118*j4
      t7042 = t7033*t6630
      t7044 = 64*t460*t7042
      t7046 = 32*t362*t7042
      t7047 = -t7015-t7019+t7021+t7023+t7027+t7029-t7032-t7036-t7038+t70
     #41-t7044-t7046
      t7049 = 64*t68*t6634
      t7051 = 4*t351*t7042
      t7052 = Er*t6624
      t7053 = t7052*CsThr
      t7055 = 32*t460*t7053
      t7056 = t158*t6787
      t7057 = 8*t7056
      t7058 = t362*t7053
      t7059 = 16*t7058
      t7061 = 2*t351*t7053
      t7063 = 32*t460*t7016
      t7065 = 16*t362*t7016
      t7067 = 32*t68*t6620
      t7069 = 2*t351*t7016
      t7071 = 32*t460*t7052
      t7073 = 16*t362*t7052
      t7075 = 32*t68*t6627
      t7076 = -t7049+t7051+t7055+t7057+t7059-t7061+t7063+t7065+t7067-t70
     #69-t7071-t7073-t7075
      t7079 = 2*t351*t7052
      t7081 = 16*t385*t6645
      t7083 = 8*t357*t6645
      t7084 = t14*t6644
      t7086 = 16*t7084*t6095
      t7087 = t32*t6644
      t7088 = t7087*t6095
      t7090 = 16*t385*t6615
      t7092 = 8*t357*t6615
      t7094 = 64*t68*t6631
      t7096 = 64*t6084*t6637
      t7097 = t158*t2
      t7099 = 16*t7097*t7034
      t7101 = 64*t461*t6634
      t7102 = t7097*t7042
      t7103 = 16*t7102
      t7104 = t7079-t7081-t7083-t7086-t7088+t7090+t7092+t7094+t7096+t709
     #9+t7101+t7103
      t7105 = t131*t6634
      t7106 = 32*t7105
      t7108 = 4*t125*t6634
      t7110 = 8*t7097*t7053
      t7112 = 8*t7097*t7017
      t7114 = 4*t6610*t6646
      t7117 = 16*t7097*t7033*t6758
      t7119 = 16*t6075*j4
      t7120 = t6204*j4
      t7121 = 32*t6790
      t7122 = t7087*j2
      t7124 = 16*t7084*j2
      t7126 = 32*t460*t7017
      t7127 = 32*t6797
      t7128 = t7106-t7108-t7110+t7112-t7114-t7117-t7119-t7120+t7121+t712
     #2+t7124-t7126-t7127
      A2e0c3p1t1n2 = 96*t427*(t7047+t7076+t7104+t7128)*t6672*t6010*t6653
     #*t432*t6948
      t7137 = 8*t6610
      t7138 = t1080*t1
      t7139 = 4*t7138
      t7140 = t1432*t437
      t7141 = 4*t7140
      A2e0c4p1t1n2 = 12*sqq*(-t7137-t6546-t6540-t6547+t443+t427+t7139+t7
     #141+t6550)*t6647*t140*t6654*t6969
      A2e0c5p1t1n2 = 0
      t7148 = -t6668+t6955-t6959-t6957+t6951+t6963-t6580-t6579+t6601-t66
     #5+t6587+t6600+t6280-t6586-t438
      t7150 = t6672*t140
      A2e0c6p1t1n2 = 24*t7148*t6647*t7150*t6654*t6969
      A2e0c1p1t2n2 = A2e0c1p1t1n2
      A2e0c2p1t2n2 = 0
      A2e0c3p1t2n2 = A2e0c3p1t1n2
      A2e0c4p1t2n2 = A2e0c4p1t1n2
      A2e0c5p1t2n2 = 0
      A2e0c6p1t2n2 = A2e0c6p1t1n2
      A1e0c1p1t1n3 = 0
      t7154 = t14*t40
      t7155 = 8*t7154
      t7156 = 16*t300
      t7157 = t452*t1
      t7158 = 4*t7157
      t7160 = j2*j3
      t7162 = 2*t118*t6726
      t7164 = 2*t118*t6709
      t7166 = 4*t118*t6881
      t7167 = Er*j3
      t7169 = 4*t5987*t7167
      t7170 = t6645*CsThl
      t7171 = t7160+t6619-t7162-t6626+t7164-t6633+t7166-t7169-t6645+t717
     #0
      t7172 = (-t6005+t7155-t7156+t5997+t7158+t32)*t7171
      t7175 = 1/(-1+CsThl)
      t7176 = t138*t7175
      t7177 = 1/j3
      t7178 = t6010*t7177
      A1e0c2p1t1n3 = -12*t7172*t6649*t7176*t7178
      t7182 = 32*t6092
      t7183 = 64*t6088
      t7184 = t7154*sqq
      t7185 = 8*t7184
      t7186 = 4*t6849
      t7187 = 32*t1882
      t7191 = 1/(1+CsThr)
      t7192 = t7191*t215
      t7194 = t7175*t6010
      t7195 = t7194*t7177
      A1e0c3p1t1n3 = -24*t7171*(-t323+t7182-t7183+t234-t7185+t322-t7186+
     #t7187-t122)*t7192*t406*t7195
      A1e0c4p1t1n3 = 0
      t7201 = t6682*t7195
      A1e0c5p1t1n3 = 12*t427*t7172*t140*t7201
      t7204 = 8*t6645*t7184
      t7207 = 32*t6645*t6659*t40
      t7208 = t6617*t40
      t7211 = 16*t68*t7208*sqq
      t7214 = 64*t6708*t6726*t40
      t7217 = 2*t118*t6726*t122
      t7218 = 32*t6695
      t7219 = 64*t6701
      t7220 = 64*t6710
      t7221 = t7204+t7207-t7211-t7214+t7217+t6694+t7218-t6698+t7219+t670
     #6+t7220
      t7222 = 128*t6721
      t7223 = 64*t6727
      t7224 = 128*t6738
      t7226 = 64*t6708*t7208
      t7227 = t6624*t40
      t7229 = 64*t6708*t7227
      t7230 = -t6714+t6720+t7222-t6725-t7223+t6731-t6737-t7224-t6745+t72
     #26-t7229
      t7233 = 32*t6645*t6092
      t7234 = t6645*t6096
      t7236 = 24*t7160*t73
      t7238 = 32*t7160*t6092
      t7240 = 64*t7160*t6088
      t7242 = 8*t7160*t76
      t7244 = 128*t6112*t7167
      t7245 = 32*t6748
      t7246 = 96*t6756
      t7247 = 128*t6760
      t7248 = -t7233-t7234-t7236+t7238-t7240+t7242-t7244+t7245-t7246-t72
     #47-t6764
      t7249 = 16*t6766
      t7253 = 32*t6125*t6630*t40*sqq
      t7256 = 128*t6139*t6754*t40
      t7259 = 4*t5988*t6754*t122
      t7263 = 32*t6125*j3*t40*sqq
      t7267 = 16*t6128*j3*t32*CsThr
      t7268 = t7249-t6789+t6796-t6799+t6802-t6806+t6810+t7253+t7256-t725
     #9+t7263+t7267
      t7274 = 128*t6125*j3*CsThr*sqq
      t7276 = CsThl*t40*sqq
      t7278 = 16*t6790*t7276
      t7280 = 16*t6797*t7276
      t7282 = 8*t7170*t7184
      t7283 = 16*t6812
      t7284 = 48*t6821
      t7285 = 64*t6824
      t7286 = 8*t6830
      t7287 = 8*t6838
      t7288 = 24*t6842
      t7289 = -t7274-t7278+t7280-t7282-t7283+t7284+t7285+t6828-t7286+t72
     #87-t7288
      t7290 = 32*t6844
      t7291 = 4*t6850
      t7292 = 128*t6882
      t7293 = 256*t6893
      t7294 = 16*t6896
      t7295 = 48*t6902
      t7296 = 64*t6905
      t7297 = -t7290-t6848+t7291+t6880+t7292-t6886+t6892+t7293+t7294-t72
     #95-t7296
      t7299 = 8*t6912
      t7302 = 32*t6125*t6754*t41
      t7303 = t7160*t122
      t7305 = 32*t7160*t233
      t7306 = -t6909+t7299-t6925+t6928-t7302-t7303+t7305-t6933+t6935-t69
     #37+t6938
      t7309 = 128*t6708*t6631*t40
      t7312 = 96*t6084*t7167*sqq
      t7315 = 128*t6112*t7167*t40
      t7318 = 256*t6112*t7167*CsThr
      t7321 = 32*t6081*t7167*t32
      t7324 = 4*t5987*t7167*t122
      t7327 = 16*t68*t7227*sqq
      t7330 = 64*t6708*t6709*t40
      t7333 = 2*t118*t6709*t122
      t7335 = 8*t7160*t7184
      t7337 = 4*t7160*t6849
      t7339 = 32*t7160*t1882
      t7340 = -t7309+t7312-t7315+t7318-t7321+t7324+t7327+t7330-t7333-t73
     #35-t7337+t7339
      t7345 = t7191*t7175
      t7348 = t7178*t140
      A1e0c6p1t1n3 = 24*Mzsq*(t7221+t7230+t7248+t7268+t7289+t7297+t7306+
     #t7340)*t7345*t138*t6068*t7348
      A1e0c7p1t1n3 = 0
      t7352 = 8*t7154*t437
      t7354 = 8*t7154*Mzsq
      t7356 = 16*t300*t437
      t7357 = 8*t6849
      t7359 = 16*t300*Mzsq
      t7360 = 16*t7184
      t7362 = 4*t2024*t452
      t7364 = 4*t427*t45
      t7365 = t6256+t7352-t6262+t7354-t7356+t6257-t7357-t6282-t7359-t628
     #0-t6281-t7360+t7187+t438+t7362+t665+t7364+t322
      t7369 = t7194*t6968
      A1e0c8p1t1n3 = 12*t7365*t7171*t140*t7177*t7369
      t7372 = 32*t1313*t452
      t7374 = 4*t1755*t527
      t7376 = 8*t1313*t41
      t7378 = 8*t7154*t443
      t7380 = 32*t300*t443
      t7382 = 4*t45*t438
      t7383 = t7372-t7374-t7376-t7378+t7380-t7382+t6381-t6382+t6383+t638
     #4-t6385-t765-t863
      t7385 = 8*t2597*CsThr
      t7386 = t2418*CsThr
      t7387 = 64*t7386
      t7389 = 128*t2253*CsThr
      t7390 = t6092*sqq
      t7391 = 64*t7390
      t7393 = 64*t747*CsThr
      t7395 = 32*t6092*Mzsq
      t7397 = 64*t852*CsThr
      t7399 = 32*t6092*t437
      t7401 = 16*t320*t14
      t7402 = t7385-t7387+t7389-t7391-t7393+t7395-t7397+t7399+t7401-t643
     #0+t6431-t6432+t6435+t6436
      t7406 = t7191*t1712*t432
      t7408 = t6649*t138
      A1e0c9p1t1n3 = 24*(t7383+t7402)*t7171*t7406*t7408*t7195
      A1e0c1p1t2n3 = 0
      A1e0c2p1t2n3 = A1e0c2p1t1n3
      A1e0c3p1t2n3 = A1e0c3p1t1n3
      A1e0c4p1t2n3 = 0
      A1e0c5p1t2n3 = A1e0c5p1t1n3
      A1e0c6p1t2n3 = A1e0c6p1t1n3
      A1e0c7p1t2n3 = 0
      A1e0c8p1t2n3 = A1e0c8p1t1n3
      A1e0c9p1t2n3 = A1e0c9p1t1n3
      A2e0c1p1t1n3 = 0
      t7411 = 4*t45
      A2e0c2p1t1n3 = 48*t665*(t7411-t6456+sqq)*t7171*t140*t7201
      t7417 = 16*t7014
      t7418 = 8*t7018
      t7419 = 8*t7026
      t7420 = 4*t7028
      t7421 = 16*t7056
      t7422 = 8*t7058
      t7423 = -t7417-t7418+t7023+t7419+t7420-t7044-t7046-t7049+t7051+t70
     #55+t7421+t7422
      t7424 = t7087*t6087
      t7425 = t7063+t7065+t7067-t7069-t7071-t7073-t7075+t7079-t7081-t708
     #3-t7086+t7094-t7424
      t7428 = 16*t300*t7160
      t7430 = 8*t357*t7160
      t7432 = 64*t6084*t7167
      t7433 = 32*t7102
      t7434 = 16*t7105
      t7436 = 4*t125*t6881
      t7438 = t7033*j3
      t7440 = 16*t452*t2*t7438
      t7441 = t7428+t7430+t7432+t7101+t7433+t7434-t7110+t7112-t7114-t743
     #6+t7440-t7117
      t7443 = 16*t6075*j3
      t7444 = t6204*j3
      t7447 = 2*t351*t7016*CsThl
      t7450 = 4*t452*t6179*j3
      t7451 = t67*CsThr
      t7453 = 64*t7451*t7438
      t7455 = 32*t362*t7438
      t7458 = 4*t6483*t118*j3
      t7461 = 2*t351*t7052*CsThl
      t7462 = -t7443-t7444+t7121+t7122+t7124+t7447-t7450-t7453-t7455+t74
     #58-t7461-t7126-t7127
      A2e0c3p1t1n3 = 96*t427*(t7423+t7425+t7441+t7462)*t7191*t7175*t6682
     #*t7348
      A2e0c4p1t1n3 = 0
      t7470 = 8*t7157
      t7471 = 4*t3664
      t7472 = 4*t5324
      A2e0c5p1t1n3 = 12*sqq*(-t7470-t6546-t6540-t6547+t443+t7471+t427+t7
     #472+t6550)*t7171*t140*t7177*t7369
      t7479 = -t7187+t7359-t7362-t7364+t7357+t7356-t6580-t6579+t6601-t66
     #5+t6587+t6600+t6280-t6586-t438
      t7481 = t7191*t140
      A2e0c6p1t1n3 = 24*t7479*t7171*t7481*t7177*t7369
      A2e0c1p1t2n3 = 0
      A2e0c2p1t2n3 = A2e0c2p1t1n3
      A2e0c3p1t2n3 = A2e0c3p1t1n3
      A2e0c4p1t2n3 = 0
      A2e0c5p1t2n3 = A2e0c5p1t1n3
      A2e0c6p1t2n3 = A2e0c6p1t1n3
      t7486 = -t6615-t6619-t6622+t6626+t6629+t6633+t6636+t6639+t6645+t66
     #46
      t7487 = (-t6005+t5997+t32-t6611+t6612+t6613)*t7486
      t7489 = t7191*t138
      A1e0c1p1t1n4 = -12*t7487*t6649*t7489*t6655
      A1e0c2p1t1n4 = 0
      t7495 = t7175*t140
      t7497 = t6011*t6654
      A1e0c3p1t1n4 = 24*t7486*(-t122+t234+t6660-t6668+t6662-t6664+t6666+
     #t322-t323)*t7495*t7192*t7497
      t7502 = t6068*t7191
      t7503 = t7502*t7497
      A1e0c4p1t1n4 = 12*t427*t7487*t140*t7503
      A1e0c5p1t1n4 = 0
      t7505 = -t6686+t6688+t6690+t6692+t6694-t6696-t6698-t6700-t6702+t67
     #03+t6706
      t7506 = -t6711-t6714-t6717+t6720-t6722-t6725+t6728+t6731+t6734-t67
     #37+t6739
      t7508 = -t6743-t6745-t6749+t6753+t6757+t6761-t6764-t6767-t6770+t67
     #73-t6777
      t7509 = t6782+t6786-t6789+t6794+t6796-t6799+t6802-t6804-t6806+t680
     #8+t6810+t6813
      t7512 = -t6819-t6822-t6825+t6828+t6831+t6834-t6837-t6839+t6841+t68
     #43+t6845
      t7513 = -t6848-t6851-t6854+t6856-t6858-t6860+t6863+t6866-t6869-t68
     #72-t6875
      t7515 = t6880-t6883-t6886-t6889+t6892-t6894-t6897+t6900+t6903+t690
     #6-t6909
      t7516 = -t6913-t6916+t6919+t6922-t6925+t6928+t6930-t6931-t6933+t69
     #35-t6937+t6938
      t7524 = t406*t6010
      A1e0c6p1t1n4 = 24*Mzsq*(t7505+t7506+t7508+t7509+t7512+t7513+t7515+
     #t7516)*t6654*t7175*t7191*t6068*t7524
      t7527 = t6256+t6951-t6952-t6668-t6262+t6257-t6282-t6280-t6281+t438
     #+t6954+t6955+t665-t6957-t6959+t6961+t6963+t322
      t7528 = -t7486
      t7532 = t432*t7191*t7497
      A1e0c7p1t1n4 = -12*t7527*t7528*t6286*t7532
      A1e0c8p1t1n4 = 0
      t7535 = -t6972+t6974-t6976-t6978+t6980-t6982+t6381-t6382+t6383+t63
     #84-t6385-t765-t863
      t7536 = -t6430+t6431-t6432-t6985+t6986-t6988-t6989+t6991+t6993+t69
     #94+t6996+t6998+t6435+t6436
      A1e0c9p1t1n4 = 24*(t7535+t7536)*t7528*t7175*t1712*t432*t6649*t7191
     #*t7497
      A1e0c1p1t2n4 = A1e0c1p1t1n4
      A1e0c2p1t2n4 = 0
      A1e0c3p1t2n4 = A1e0c3p1t1n4
      A1e0c4p1t2n4 = A1e0c4p1t1n4
      A1e0c5p1t2n4 = 0
      A1e0c6p1t2n4 = A1e0c6p1t1n4
      A1e0c7p1t2n4 = A1e0c7p1t1n4
      A1e0c8p1t2n4 = 0
      A1e0c9p1t2n4 = A1e0c9p1t1n4
      A2e0c1p1t1n4 = 48*t665*(-sqq+t6456+t7008)*t7528*t140*t7503
      A2e0c2p1t1n4 = 0
      t7550 = -t7015-t7019+t7021-t7023+t7027+t7029-t7032-t7036+t7038-t70
     #41-t7044+t7046
      t7551 = -t7049-t7051-t7055+t7057+t7059-t7061+t7063-t7065+t7067+t70
     #69-t7071+t7073-t7075
      t7553 = -t7079-t7081+t7083-t7086-t7088+t7090-t7092-t7094-t7096+t70
     #99-t7101+t7103
      t7554 = t7106-t7108+t7110-t7112+t7114+t7117+t7119+t7120-t7121-t712
     #2-t7124+t7126+t7127
      A2e0c3p1t1n4 = -96*t427*(t7550+t7551+t7553+t7554)*t6654*t7175*t750
     #2*t7524
      A2e0c4p1t1n4 = 12*sqq*(-t7137+t6546+t6540+t6547-t443-t427+t7139+t7
     #141-t6550)*t7528*t140*t1712*t7532
      A2e0c5p1t1n4 = 0
      t7569 = -t6957+t6963-t6959+t6951+t6955-t6668+t6586-t6587+t6580+t66
     #5-t6600-t6601+t438-t6280+t6579
      A2e0c6p1t1n4 = -24*t7569*t7528*t7497*t7345*t2155*t140
      A2e0c1p1t2n4 = A2e0c1p1t1n4
      A2e0c2p1t2n4 = 0
      A2e0c3p1t2n4 = A2e0c3p1t1n4
      A2e0c4p1t2n4 = A2e0c4p1t1n4
      A2e0c5p1t2n4 = 0
      A2e0c6p1t2n4 = A2e0c6p1t1n4
      A1e0c1p1t1n5 = 0
      t7577 = t7160+t6619+t7162-t6626-t7164-t6633-t7166-t7169-t6645-t717
     #0
      t7578 = (-t6005+t7155+t7156+t5997-t7158+t32)*t7577
      t7580 = t138*t6672
      A1e0c2p1t1n5 = 12*t7578*t6649*t7580*t7178
      t7585 = t6652*t140
      t7587 = t6944*t7177
      A1e0c3p1t1n5 = -24*t7577*(t7186+t322-t122-t7187-t7185+t7183+t7182+
     #t234-t323)*t7585*t6058*t7587
      A1e0c4p1t1n5 = 0
      t7593 = t6682*t7587
      A1e0c5p1t1n5 = -12*t427*t7578*t140*t7593
      t7596 = t7204-t7207-t7211+t7214-t7217+t6694-t7218-t6698-t7219+t670
     #6-t7220
      t7597 = -t6714+t6720-t7222-t6725+t7223+t6731-t6737+t7224-t6745+t72
     #26-t7229
      t7599 = -t7233+t7234-t7236+t7238+t7240+t7242-t7244-t7245+t7246+t72
     #47-t6764
      t7600 = -t7249-t6789+t6796-t6799+t6802-t6806+t6810+t7253-t7256+t72
     #59+t7263-t7267
      t7603 = t7274+t7278-t7280+t7282+t7283-t7284-t7285+t6828+t7286-t728
     #7+t7288
      t7604 = t7290-t6848-t7291+t6880-t7292-t6886+t6892-t7293-t7294+t729
     #5+t7296
      t7606 = -t6909-t7299-t6925+t6928+t7302-t7303+t7305-t6933+t6935-t69
     #37+t6938
      t7607 = -t7309+t7312-t7315-t7318-t7321+t7324+t7327-t7330+t7333-t73
     #35+t7337-t7339
      A1e0c6p1t1n5 = 24*Mzsq*(t7596+t7597+t7599+t7600+t7603+t7604+t7606+
     #t7607)*t7580*t432*t429*t6652*t7348
      A1e0c7p1t1n5 = 0
      t7617 = t6256+t7352-t6262+t7354+t7356+t6257+t7357-t6282+t7359-t628
     #0-t6281-t7360-t7187+t438-t7362+t665-t7364+t322
      t7620 = t6229*t7587
      A1e0c8p1t1n5 = -12*t7617*t7577*t6286*t7620
      t7623 = -t7372+t7374-t7376-t7378-t7380+t7382+t6381-t6382+t6383+t63
     #84-t6385-t765-t863
      t7624 = -t7385+t7387-t7389-t7391+t7393+t7395+t7397+t7399+t7401-t64
     #30+t6431-t6432+t6435+t6436
      t7627 = t6652*t1712
      A1e0c9p1t1n5 = 24*(t7623+t7624)*t7577*t7627*t432*t7408*t7587
      A1e0c1p1t2n5 = 0
      A1e0c2p1t2n5 = A1e0c2p1t1n5
      A1e0c3p1t2n5 = A1e0c3p1t1n5
      A1e0c4p1t2n5 = 0
      A1e0c5p1t2n5 = A1e0c5p1t1n5
      A1e0c6p1t2n5 = A1e0c6p1t1n5
      A1e0c7p1t2n5 = 0
      A1e0c8p1t2n5 = A1e0c8p1t1n5
      A1e0c9p1t2n5 = A1e0c9p1t1n5
      A2e0c1p1t1n5 = 0
      A2e0c2p1t1n5 = 48*t665*(t7411+t6456-sqq)*t7577*t140*t7593
      t7637 = -t7417-t7418-t7023+t7419+t7420-t7044+t7046-t7049-t7051-t70
     #55+t7421+t7422
      t7638 = t7063-t7065+t7067+t7069-t7071+t7073-t7075-t7079-t7081+t708
     #3-t7086-t7094-t7424
      t7640 = t7428-t7430-t7432-t7101+t7433+t7434+t7110-t7112+t7114-t743
     #6+t7440+t7117
      t7641 = t7443+t7444-t7121-t7122-t7124+t7447-t7450-t7453+t7455-t745
     #8-t7461+t7126+t7127
      A2e0c3p1t1n5 = -96*t427*(t7637+t7638+t7640+t7641)*t6672*t138*t6068
     #*t6652*t7348
      A2e0c4p1t1n5 = 0
      A2e0c5p1t1n5 = 12*sqq*(-t7470+t6546+t6540+t6547-t443+t7471-t427+t7
     #472-t6550)*t7577*t140*t1712*t7620
      t7657 = -t7364+t7356-t7362+t7357+t7359-t7187+t6586-t6587+t6580+t66
     #5-t6600-t6601+t438-t6280+t6579
      t7660 = t7627*t140
      A2e0c6p1t1n5 = -24*t7657*t7577*t7587*t6229*t7660
      A2e0c1p1t2n5 = 0
      A2e0c2p1t2n5 = A2e0c2p1t1n5
      A2e0c3p1t2n5 = A2e0c3p1t1n5
      A2e0c4p1t2n5 = 0
      A2e0c5p1t2n5 = A2e0c5p1t1n5
      A2e0c6p1t2n5 = A2e0c6p1t1n5
      t7664 = t6644*CsThr
      t7668 = log(DCMPLX(Eq/(Eq-Er)))
      t7669 = t7668*CsThr
      t7670 = j4-t6644+t7664-t7668+t7669
      t7671 = 16*t6667
      t7672 = 32*t337
      t7673 = 16*t80
      t7674 = t1992*t1756
      t7675 = 8*t7674
      t7678 = 2*t44*t118
      t7680 = 8*t362*t1878
      t7681 = 32*t548
      t7682 = t158*t54
      t7683 = 8*t7682
      t7684 = 16*t127
      t7685 = 12*t129
      t7686 = 64*t461
      t7687 = 32*t482
      t7688 = t7678+t7680-t169+t7681-t7683+t7684-t7685-t132+t82-t7686+t7
     #687
      t7689 = -t6662-t122-t6664+t7671+t77-t78-t7672+t126+t7673+t7675+t76
     #88
      t7692 = t6652*t6654
      A1e0c1p1t1n6 = -12*t7670*t7689*t142*t406*t7692
      A1e0c2p1t1n6 = 0
      t7696 = 16*t57
      t7697 = 16*t337
      t7698 = 16*t7674
      t7699 = 8*t131
      t7700 = 32*t461
      t7701 = -t6662+t7696+t234-t122-t6664+t6668+t322-t323+t7697+t207-t7
     #698+t7678+t7680-t169-t208-t7699-t6660+t6666+t7700-t7687
      t7705 = t138*t404
      t7706 = t7585*t6654
      A1e0c3p1t1n6 = 24*t7670*t7701*t6672*t402*t7705*t7706
      t7713 = t418*Mzsq
      t7717 = t140/(-t1368+t2164+t418*sll-t7713)*t7692
      A1e0c4p1t1n6 = 12*t7670*Mzsq*t7689*t138*t7717
      A1e0c5p1t1n6 = 0
      t7720 = 2*t1185*t125
      t7722 = 32*t460*t2556
      t7723 = t202*Er
      t7725 = 64*t460*t7723
      t7727 = 8*t39*t444
      t7728 = t385*t448
      t7729 = 64*t7728
      t7731 = 16*t39*t448
      t7733 = 8*t1080*t439
      t7735 = 64*t460*t2455
      t7737 = 32*t460*t2142
      t7739 = 32*t460*t2452
      t7741 = 8*t250*t2597
      t7743 = 128*t6659*t202
      t7745 = 64*t6659*t416
      t7746 = t385*t184
      t7747 = 64*t7746
      t7749 = 64*t747*t158
      t7750 = -t7720-t7722+t7725+t7727+t7729-t7731-t7733+t7735-t7737-t77
     #39-t7741-t7743+t7745+t7747+t7749
      t7752 = 64*t6665*t202
      t7754 = 32*t6665*t416
      t7756 = 16*t39*t184
      t7758 = 32*t747*t64
      t7760 = 4*t765*t1432
      t7762 = 8*t1313*t845
      t7763 = t1313*t44
      t7764 = 32*t7763
      t7766 = 64*t6659*t887
      t7768 = 64*t6659*t443
      t7770 = 128*t6659*t447
      t7772 = 32*t6665*t887
      t7774 = 32*t852*t64
      t7776 = 64*t6665*t447
      t7778 = 8*t845*t867
      t7780 = 4*t863*t1432
      t7781 = t7752-t7754-t7756-t7758+t7760+t7762-t7764+t7766+t7768-t777
     #0-t7772-t7774+t7776+t7778+t7780
      t7783 = t44*t867
      t7784 = 32*t7783
      t7785 = 8*t824
      t7786 = 16*t836
      t7787 = 16*t839
      t7788 = 4*t847
      t7789 = 32*t850
      t7790 = 4*t905
      t7791 = 16*t910
      t7792 = 16*t929
      t7793 = 16*t933
      t7794 = -t7784+t791-t803+t5054-t7785-t7786-t7787+t7788+t7789+t890-
     #t903+t7790-t7791-t7792-t7793
      t7803 = 32*t936+t5060-8*t940+4*t945+4*t1058+4*t1164-8*t1225+4*t130
     #9-t1312+t1315-t1317+t1319-t1321+t1323-16*t1332
      t7810 = 2*t1346
      t7815 = 2*t1355
      t7816 = 4*t1359
      t7819 = -16*t1334+32*t1338+4*t1340+8*t1342-t7810+32*t1347-16*t1349
     #-16*t1351-8*t1353-t7815+t7816-16*t1364+8*t1366-t1371+t1330
      t7820 = t1328+t1400+t1401-t1402-t1403+t1404+t1405-t1406-t1474-t147
     #5-t1476+t1477+t1478-t1479-t1480
      t7821 = t7819+t7820
      t7823 = 4*t1181*t125
      t7826 = 32*t1992*Mzsq*t2452
      t7829 = 16*t1992*t437*t7723
      t7832 = 16*t191*Mzsq*t2452
      t7835 = 8*t191*t437*t7723
      t7837 = 16*t7097*t1989
      t7838 = t44*Er
      t7841 = 2*t7838*t2436*Mzsq
      t7844 = 8*t362*t38*t1989
      t7846 = 2*t863*t1739
      t7848 = 4*t1184*t507
      t7850 = 32*t460*t1989
      t7852 = 32*t481*t1989
      t7854 = 4*t76*t1084
      t7855 = t385*t511
      t7856 = 32*t7855
      t7857 = t1481+t7823-t7826+t7829+t7832-t7835+t7837-t7841-t7844-t784
     #6+t7848-t7850+t7852+t7854-t7856
      t7859 = 8*t39*t511
      t7860 = t385*t444
      t7861 = 32*t7860
      t7863 = 8*t191*t2392
      t7865 = 16*t1992*t2392
      t7867 = 32*t481*t2452
      t7869 = 32*t481*t2142
      t7871 = 64*t481*t2455
      t7872 = t184*Er
      t7874 = 16*t191*t7872
      t7876 = 8*t191*t2368
      t7878 = 16*t1992*t2368
      t7880 = 4*t250*t2140
      t7882 = 2*t765*t1739
      t7884 = 32*t1992*t7872
      t7886 = 64*t481*t7723
      t7888 = 32*t481*t2556
      t7889 = t7859-t7861-t7863+t7865+t7867+t7869-t7871+t7874-t7876+t787
     #8+t7880-t7882-t7884-t7886+t7888
      t7894 = t7150*t1712
      t7897 = t1714*t432*t414
      t7899 = t7897*t6653*t6654
      A1e0c6p1t1n6 = 12*(t7750+t7781+t7794+t7803+t7821+t7857+t7889)*t767
     #0*t7894*t7899
      t7901 = 16*t6971
      t7902 = 16*t6981
      t7903 = t427*t523
      t7904 = 8*t7903
      t7905 = t443*t523
      t7906 = 8*t7905
      t7907 = t1992*t2452
      t7908 = 8*t7907
      t7910 = 2*t1184*t125
      t7912 = 8*t191*t2452
      t7913 = t1992*t7723
      t7914 = 16*t7913
      t7917 = 16*t362*t1878*sll
      t7919 = 4*t44*t2875
      t7920 = t158*t205
      t7921 = 16*t7920
      t7922 = t7901-t6974-t6976-t6978-t6980+t7902-t7904-t7906+t7908+t791
     #0+t7912-t7914-t7917-t7919+t7921
      t7923 = t1992*t2556
      t7924 = 8*t7923
      t7926 = 8*t191*t2556
      t7928 = 2*t1080*t125
      t7930 = 32*t1184*t57
      t7932 = 32*t39*t1929
      t7934 = 32*t481*t2121
      t7936 = 64*t460*t2121
      t7938 = 32*t1080*t57
      t7939 = 64*t4174
      t7941 = 8*t1432*t184
      t7943 = 32*t385*t202
      t7945 = 16*t39*t202
      t7947 = 64*t460*t2350
      t7949 = 32*t481*t2350
      t7951 = 32*t39*t1924
      t7952 = 64*t2159
      t7953 = t7924+t7926+t7928+t7930-t7932+t7934-t7936+t7938-t7939+t794
     #1-t7943+t7945-t7947+t7949-t7951+t7952
      t7956 = 64*t481*t1988
      t7958 = 128*t460*t1988
      t7959 = 24*t1843
      t7960 = 8*t1847
      t7961 = 16*t1920
      t7962 = 12*t1925
      t7963 = 12*t1930
      t7964 = 16*t1932
      t7965 = -t7956+t7958+t7959+t7960-t2037-t2039+t7961+t1923-t7962-t79
     #63+t7964-t1935+t1940-t1942-t1946+t1948
      t7967 = 16*t1917
      t7969 = 16*t1911
      t7970 = -32*t1949-t765+t1954-t863-t1956-t1958+t1959+t1960-t1961+t1
     #965+t7967+t1968+t1969-32*t1914-t1972+t7969
      t7978 = t402/(t141-t416-t887+t447)*t7692
      A1e0c7p1t1n6 = 12*t7670*(t7922+t7953+t7965+t7970)*t406*t7978
      A1e0c8p1t1n6 = 0
      t7980 = 4*t2374
      t7981 = 2*t2378
      t7982 = 2*t2381
      t7994 = t2255-t7980+t7981+t7982+16*t2382-8*t2386-8*t2390-4*t2399+1
     #6*t2406+16*t2410-4*t2412+8*t2414-32*t2416+t2420-24*t2427+16*t2429+
     #t2432
      t7995 = 6*t2439
      t7996 = 4*t2442
      t8005 = -t2435+t7995-t7996+t2462-t2464-t2466+t2468-t2470+t2472-t53
     #62+48*t2476+16*t2480+16*t2483+8*t2486-32*t2489-12*t2491+48*t2495-3
     #2*t2497
      t8006 = t7994+t8005
      t8008 = 4*t2543
      t8010 = 6*t2551
      t8019 = -t2536-24*t2537+t2540-t8008+16*t2545+t8010+8*t2565+48*t256
     #7-32*t2570-t5372-64*t2577-12*t2579+48*t2581+16*t2591-t2599-t2601+4
     #8*t2634
      t8020 = 12*t2637
      t8021 = 8*t2640
      t8022 = 48*t2657
      t8023 = 32*t2661
      t8024 = 12*t2666
      t8026 = 12*t1325*t7138
      t8029 = 8*t769*t3451*t1
      t8030 = 16*t3185
      t8031 = 4*t3191
      t8032 = 4*t3202
      t8033 = 64*t3231
      t8034 = 8*t3233
      t8035 = t2310*t2333
      t8036 = 96*t8035
      t8037 = t6665*t32
      t8039 = 96*t416*t8037
      t8041 = 64*t2164*t6987
      t8042 = CsThl*t141
      t8044 = 12*t8042*t766
      t8045 = CsThl*t2195
      t8046 = t8045*Mzsq
      t8048 = 64*t2253*t8046
      t8050 = 16*t2503*t935
      t8051 = -t8020+t8021+t8022-t8023-t8024-t8026+t8029+t8030-t8031-t80
     #32-t8033+t8034+t8036+t8039-t8041-t8044-t8048+t8050
      t8055 = 128*t2499*t1181
      t8056 = t2418*t8046
      t8057 = 32*t8056
      t8059 = 8*t2418*t4686
      t8061 = 8*t2418*t4691
      t8064 = 128*t751*t8045*sqq
      t8066 = 64*t751*t8046
      t8068 = 192*t3599*t1188
      t8070 = 128*t3599*t1181
      t8072 = 64*t751*t4683
      t8074 = 32*t751*t4686
      t8076 = 96*t887*t8037
      t8078 = 32*t751*t4691
      t8081 = 12*t1184*t2589*t1
      t8084 = 8*t769*t8045*t1
      t8086 = 96*t6665*t899
      t8088 = 64*t6665*t2176
      t8090 = 12*t1325*t7140
      t8091 = t8055+t8057-t8059-t8061+t8064-t8066-t8068+t8070-t8072+t807
     #4+t8076+t8078-t8081+t8084+t8086-t8088-t8090
      t8093 = 24*t2310*t944
      t8094 = t1184*t32
      t8095 = t2310*t8094
      t8096 = 96*t8095
      t8098 = 24*t177*t4673
      t8099 = t3451*t437
      t8101 = 64*t751*t8099
      t8102 = t447*t2367
      t8104 = 32*t192*t8102
      t8105 = t2703*t1756
      t8107 = 8*t192*t8105
      t8108 = t3489*CsThl
      t8109 = t447*t1756
      t8111 = 128*t8108*t8109
      t8112 = t152*CsThl
      t8114 = 64*t8112*t8102
      t8116 = 16*t8112*t8105
      t8117 = t3489*t38
      t8119 = 128*t8117*t8109
      t8120 = t2404*t2350
      t8122 = 8*t192*t8120
      t8124 = 16*t8112*t8120
      t8125 = t2433*t38
      t8127 = 16*t8125*t8109
      t8128 = t2433*CsThl
      t8130 = 32*t8128*t8109
      t8131 = t447*t769
      t8133 = 8*t1739*t8131
      t8135 = 4*t1739*t2799
      t8136 = t2727*t184
      t8138 = 2*t1739*t8136
      t8139 = t2703*t184
      t8141 = 2*t1739*t8139
      t8142 = -t8093+t8096-t8098-t8101-t8104+t8107-t8111+t8114-t8116+t81
     #19+t8122-t8124-t8127+t8130-t8133-t8135+t8138+t8141
      t8145 = 16*t8042*t2512
      t8147 = 128*t6665*t1356
      t8149 = 32*t2310*t6305
      t8151 = 256*t3599*t44
      t8152 = t2310*t6096
      t8153 = 128*t8152
      t8155 = 4*t3451*t864
      t8157 = 24*t2310*t846
      t8159 = 192*t3599*t1081
      t8160 = 32*t3561
      t8161 = 16*t3567
      t8162 = 8*t3582
      t8163 = 32*t3594
      t8164 = 4*t3596
      t8165 = 8*t3686
      t8167 = 24*t192*t2351
      t8168 = t8145-t8147+t8149+t8151-t8153-t8155-t8157-t8159-t8160+t816
     #1-t8162+t3585+t3587+t8163-t8164+t8165+t8167
      t8169 = t416*t1324
      t8171 = 6*t1739*t8169
      t8172 = t1368*t122
      t8174 = 6*t1739*t8172
      t8175 = t2164*t122
      t8177 = 4*t1739*t8175
      t8178 = t2979*t38
      t8180 = 96*t8178*t2556
      t8182 = 64*t8108*t2444
      t8184 = 32*t8108*t2447
      t8187 = 8*t1184*t141*t1964
      t8190 = 4*t184*CsThl*t3681
      t8192 = 128*t853*t838
      t8194 = 16*t1481*t7138
      t8195 = t2264*t38
      t8197 = 32*t8195*t2447
      t8199 = 96*t8178*t2452
      t8201 = 64*t8178*t2455
      t8203 = 32*t773*t935
      t8205 = 256*t3655*t448
      t8207 = 16*t192*t2396
      t8209 = 32*t8112*t2396
      t8211 = 64*t8117*t2444
      t8212 = t8171+t8174-t8177-t8180-t8182+t8184+t8187-t8190-t8192+t819
     #4-t8197-t8199+t8201+t8203+t8205-t8207+t8209+t8211
      t8217 = 32*t8117*t2447
      t8218 = t2754*t2405
      t8219 = 32*t8218
      t8221 = 8*t2758*t2405
      t8223 = 8*t2758*t2409
      t8225 = 24*t192*t2426
      t8227 = 32*t8108*t2458
      t8229 = 96*t8108*t2392
      t8230 = t2754*t904
      t8231 = 128*t8230
      t8232 = t2754*t2409
      t8233 = 32*t8232
      t8235 = 48*t8112*t2426
      t8237 = 32*t8117*t2458
      t8239 = 96*t8117*t2392
      t8241 = 24*t8125*t2392
      t8243 = 16*t2791*t448
      t8244 = t2979*CsThl
      t8246 = 96*t8244*t2452
      t8248 = 64*t8244*t2455
      t8249 = t2788*t448
      t8250 = 64*t8249
      t8251 = -t8217+t8219-t8221-t8223+t8225+t8227+t8229-t8231+t8233-t82
     #35-t8237-t8239+t8241+t8243+t8246-t8248-t8250
      t8253 = 48*t8128*t2392
      t8256 = 16*t3747*t38*t2455
      t8257 = t2305*t38
      t8259 = 8*t8257*t2458
      t8262 = 64*t2256*CsThl*t2455
      t8265 = 32*t3747*CsThl*t2455
      t8266 = t2305*CsThl
      t8268 = 16*t8266*t2458
      t8271 = 64*t2256*t38*t2455
      t8273 = 8*t8257*t2447
      t8275 = 16*t8266*t2447
      t8278 = 4*t250*t2163*t439
      t8279 = t2264*CsThl
      t8281 = 32*t8279*t2447
      t8282 = -t8253-t8256+t8259-t8262+t8265-t8268+t8271+t8273-t8275-t82
     #78+t8281-t3796-t3798+t3800-t3802+t3804+t3806+t3811
      t8287 = 8*t3851
      t8288 = -t3813+t3815-t3816-t3817+t3819+16*t3828-64*t3830-t3833+t38
     #35-t3837-t3839+t3842+t3844-t3846+t3848+32*t3849-t8287
      t8290 = 192*t751*t8094
      t8291 = t177*t4648
      t8292 = 96*t8291
      t8294 = t177*t8045*t32
      t8295 = 64*t8294
      t8297 = 16*t177*t4678
      t8299 = 32*t2253*t4686
      t8300 = t3856+t3858-t3860+t3862+t3864-t3866-t3868+t3870-t3872+t387
     #6+t3878-t3880+t3882-t8290+t8292-t8295+t8297+t8299
      t8304 = 32*t8279*t2458
      t8306 = 32*t8195*t2458
      t8308 = 16*t192*t2362
      t8310 = 48*t8112*t2351
      t8312 = 32*t8112*t2362
      t8314 = 96*t8108*t2368
      t8316 = 64*t8108*t2553
      t8318 = 96*t8117*t2368
      t8320 = 64*t8117*t2553
      t8322 = 48*t8128*t2368
      t8324 = 96*t8244*t2556
      t8326 = 24*t8125*t2368
      t8327 = t887*t1324
      t8329 = 6*t1739*t8327
      t8330 = t4119-t4128-t4256+t4258+t8304-t8306-t8308-t8310+t8312+t831
     #4-t8316-t8318+t8320-t8322+t8324+t8326+t8329
      t8332 = 6*t1739*t4702
      t8333 = t2478*t122
      t8335 = 4*t1739*t8333
      t8337 = 4*t1739*t2373
      t8339 = 2*t1739*t2377
      t8341 = 2*t1739*t2380
      t8342 = 16*t4650
      t8343 = 64*t4665
      t8344 = t8332-t8335-t8337+t8339+t8341-t4385-t4388+t4390+t4392-t439
     #3+t4395-t4397-t4399-t4400-t4632-t4638+t8342-t8343
      t8346 = 48*t4669
      t8347 = 32*t4671
      t8348 = 12*t4674
      t8349 = 8*t4679
      t8350 = 32*t4684
      t8351 = 16*t4687
      t8352 = 16*t4692
      t8353 = 48*t4698
      t8354 = 12*t4700
      t8355 = 8*t4777
      t8356 = 16*t4779
      t8357 = 4*t4783
      t8358 = 16*t4787
      t8359 = 4*t4791
      t8361 = 128*t2979*t7838
      t8362 = t8346-t8347-t8348-t5944+t8349-t5945-t8350+t8351+t8352+t835
     #3-t8354+t8355+t8356-t8357+t8358-t8359-t8361
      t8365 = 32*t2433*t6305*Er
      t8368 = 128*t2979*t845*Er
      t8370 = 24*t177*t2636
      t8372 = 16*t177*t2639
      t8374 = 192*t751*t2333
      t8377 = 128*t751*t2408*CsThl
      t8378 = t177*t4646
      t8379 = 96*t8378
      t8381 = t177*t2388*CsThl
      t8382 = 64*t8381
      t8383 = t2503*t1181
      t8384 = 64*t8383
      t8385 = t2418*t8099
      t8386 = 32*t8385
      t8388 = 32*t2253*t4691
      t8390 = 4*t6096*t3681
      t8392 = 64*t852*t846
      t8394 = 8*t1330*t7138
      t8396 = 64*t2253*t8099
      t8399 = 64*t2433*t6096*Er
      t8400 = t5615*t1324
      t8402 = 8*t2011*t8400
      t8404 = 96*t6665*t496
      t8405 = -t8365+t8368-t8370+t8372-t8374+t8377+t8379-t8382-t8384+t83
     #86+t8388-t8390-t8392+t8394-t8396+t8399-t8402+t8404
      t8414 = t7706*t7705*t1980
      A1e0c9p1t1n6 = -12*(t8006+t8019+t8051+t8091+t8142+t8168+t8212+t825
     #1+t8282+t8288+t8300+t8330+t8344+t8362+t8405)*t7670*t6672*t4937*t84
     #14
      A1e0c1p1t2n6 = A1e0c1p1t1n6
      A1e0c2p1t2n6 = 0
      A1e0c3p1t2n6 = A1e0c3p1t1n6
      A1e0c4p1t2n6 = A1e0c4p1t1n6
      A1e0c5p1t2n6 = 0
      A1e0c6p1t2n6 = A1e0c6p1t1n6
      A1e0c7p1t2n6 = A1e0c7p1t1n6
      A1e0c8p1t2n6 = 0
      A1e0c9p1t2n6 = A1e0c9p1t1n6
      t8417 = 16*t7682
      t8418 = 24*t7674
      t8419 = 32*t1074
      t8420 = t8417+t4949-t8418-t8419+t4950+t7681-t4947-t126+t6662-t77-t
     #7678+t122
      A2e0c1p1t1n6 = -48*t7670*t8420*Mzsq*t138*t7717
      A2e0c2p1t1n6 = 0
      t8426 = 32*t7728
      t8427 = 32*t7746
      t8428 = 16*t7763
      t8429 = 16*t7783
      t8430 = t7720-t5024+t5026+t5028-t7722+t7725-t8426+t7733+t7735-t773
     #7-t7739+t7741-t8427-t7760+t8428-t7780+t8429-t5053
      t8431 = t803-t5058+t903-t7790+t5059+t5093+t5096-t5098+t1312-t5145+
     #t1317-t5146+t5147-t1323-t5149+t7810+t7815-t7816
      t8434 = 16*t1181*t129
      t8436 = 8*t1185*t129
      t8438 = 8*t2072*t5317
      t8439 = t5152-t5153+t1371-t1330-t1328-t1401-t1404-t5190+t8434-t843
     #6-t8438+t1474+t1475+t1476-t5195+t5196-t1481-t7823
      t8442 = t202*t49
      t8457 = t447*t49
      t8460 = t2024*t49
      t8463 = t7841+t7846-t7848-t7850-t7854+16*t7855+16*t7860-t7880+t788
     #2-64*t385*t8442+32*t427*t548+16*t250*t2138-8*t665*t523+32*t385*t20
     #00+32*t158*t5268+32*t385*t2008-64*t385*t8457-8*t44*t8460
      A2e0c3p1t1n6 = -48*(t8430+t8431+t8439+t8463)*t7670*t7894*t7899
      t8482 = t6974+t6980+16*t7903+16*t7905-24*t7907-t7910+48*t7913+t791
     #9-32*t7920-24*t7923-t7928+64*t3352-32*t1992*t50*Mzsq-32*t1992*t193
     #6+t7930+t7938-t7939-t7941
      t8483 = t5264-t5265-t5269-t1923+t5270+t5271+t5272-t5273+t1942+t527
     #4-t1948+t765-t1954+t863+t1956-t1959-t1965-t5287
      A2e0c4p1t1n6 = -12*t7670*(t8482+t8483)*t406*t7978
      A2e0c5p1t1n6 = 0
      t8490 = 64*t2788*t8457
      t8491 = t1929*t32
      t8493 = 96*t2754*t8491
      t8495 = t122*t49*CsThl
      t8497 = 24*t5319*t8495
      t8498 = t2163*t49
      t8499 = t8498*t437
      t8501 = 32*t6667*t8499
      t8502 = t2195*t49
      t8503 = t8502*Mzsq
      t8505 = 32*t6667*t8503
      t8506 = t1924*CsThl
      t8508 = 16*t864*t8506
      t8511 = 8*t76*t2195*t8506
      t8513 = t1929*CsThl
      t8515 = 8*t76*t2163*t8513
      t8517 = 32*t2754*t8499
      t8519 = -t7981-t7982-t5346+t5347+t5349-t5353+t5354-t5355-t2432+t53
     #56-t7995
      t8521 = t7996-t5357+t2464+t5358-t2468+t5361-t2472+t5367+t5368-t254
     #0
      t8522 = t165*t2163
      t8524 = 16*t8522*t2068
      t8526 = 24*t5317*t8495
      t8527 = t1924*t32
      t8529 = 96*t2754*t8527
      t8530 = t8498*sqq
      t8532 = 64*t2754*t8530
      t8533 = t1929*sqq
      t8535 = 96*t2788*t8533
      t8536 = t8502*sqq
      t8538 = 64*t2754*t8536
      t8540 = 32*t2754*t8503
      t8541 = t8008-t5370-t8010+t2599-t8524+t8526-t8529+t8532-t8535+t853
     #8-t8540
      t8546 = t165*t2195
      t8549 = t2340*t32
      t8554 = t2071*Mzsq
      t8567 = t8081-t8084+t8090-48*t8095-t8111+t8133+t8135-t8138-t8141-3
     #2*t5319*t2067*t1080+8*t8546*t2071*t1080
      t8568 = 64*t3155*t8457-16*t8546*t2068+24*t8549*t8513+24*t8549*t850
     #6-96*t2788*t8554+t8026-t8029-48*t8035+t8044-16*t8056+t8567
      t8571 = 8*t8522*t2071*t1184
      t8572 = t1924*t443
      t8574 = 128*t2754*t8572
      t8575 = 64*t8152
      t8579 = 16*t2340*sqq*t1929*t1080
      t8581 = -t3587-t5636+t8164-t8171-t8174+t8177-t8182+t8184-t8187+t81
     #90-t8194
      t8589 = -16*t8218+t8227+t8229+64*t8230-16*t8232+t8246-t8248+32*t82
     #49-t8262+t8278
      t8590 = t8281+t5672+t5717-t3800+t5718-t3804-t3815+t3816+t3817-t381
     #9-t3835
      t8592 = t3837+t5759-t5760-t3844+t5761-t5762+t8287-t5763+t5764+t386
     #8
      t8593 = 48*t8291
      t8594 = 32*t8294
      t8595 = -t3870+t5765-t5766-t3878+t3880-t8593+t8594-t4119+t4128-t58
     #37+t5839
      t8599 = t8335+t8337-t8339-t8341+t4385+t4388-t4390-t5883+t4393-t439
     #5+t4397
      t8601 = t4399+t4400-t5901+t5904-t5906-t5908+t5910-t5921-t5923+t593
     #9
      t8602 = 48*t8378
      t8603 = 32*t8381
      t8604 = 32*t8383
      t8605 = 16*t8385
      t8607 = 128*t2310*t2068
      t8609 = 32*t2340*t8495
      t8610 = t5941-t8361-t8602+t8603+t8604-t8605+t8390-t8394+t8402+t860
     #7-t8609
      A2e0c6p1t1n6 = 12*(-t4258-t8332-t8329+t8324+t5860-t5857-t5845-t831
     #6+t8314+t8304+t5617-t3585-t5634+t8610+t8601+t8595+t8599-t8501+t859
     #2+t8590+t8574+t8575-t8579+t8581+t8589+t8497+t8568+t8571-t8493+t816
     #2+t8155-t8145+t8490+t8541+t8521+t8519+t8515-t8517-t8505-t8508+t851
     #1+t7980)*t7670*t6672*t4937*t8414
      A2e0c1p1t2n6 = A2e0c1p1t1n6
      A2e0c2p1t2n6 = 0
      A2e0c3p1t2n6 = A2e0c3p1t1n6
      A2e0c4p1t2n6 = A2e0c4p1t1n6
      A2e0c5p1t2n6 = 0
      A2e0c6p1t2n6 = A2e0c6p1t1n6
      A1e0c1p1t1n7 = 0
      t8619 = t6644*CsThl
      t8620 = j3-t6644+t8619
      t8621 = 32*t331
      t8622 = 8*t294
      t8623 = 16*t123
      t8624 = 4*t205
      t8626 = 32*t300*t49
      t8627 = 4*t175
      t8629 = 8*t452*t54
      t8630 = -t8621-t8622+t8623-t8624+t8626+t183-t204+t8627+t207-t7673-
     #t8629
      t8631 = t177*t40
      t8632 = 8*t8631
      t8633 = t456*t4969
      t8634 = 2*t8633
      t8636 = 8*t195*t1988
      t8637 = t300*sll
      t8638 = 16*t8637
      t8639 = t67*t40
      t8641 = 32*t8639*Er
      t8642 = t7451*Er
      t8643 = 64*t8642
      t8645 = 16*t195*t1756
      t8646 = t456*t357
      t8647 = 4*t8646
      t8648 = t4950-t132+t82+t184+t8632-t8634+t8636-t8638+t8641-t8643+t8
     #645+t8647
      t8649 = t8630+t8648
      t8652 = t7175*t7177
      t8653 = t406*t8652
      A1e0c2p1t1n7 = 12*t8620*t8649*t215*t8653
      t8655 = t527*t118
      t8657 = 4*t45*t171
      t8658 = t8629-t8645+t8655-t8657+t234-t7183-t8627+t295+t307+t319+t7
     #182+t322-t323
      t8659 = 32*t8642
      t8660 = 16*t1071
      t8661 = -t7185-t7186+t7187-t8641+t8659-t8660-t122+t8621-t8623-t372
     #+t125-t377+t8633
      A1e0c3p1t1n7 = 24*t8620*(t8658+t8661)*t7191*t404*t138*t402*t7495*t
     #7177
      A1e0c4p1t1n7 = 0
      t8670 = t8620*Mzsq
      t8675 = 1/(-t665+t2408+t418*sqq-t7713)
      A1e0c5p1t1n7 = -12*t8670*t8649*t138*t140*t8675*t8652
      t8681 = 8*t1238*t511
      t8683 = 2*t457*t1244
      t8684 = t195*sll
      t8686 = 16*t8684*t2556
      t8688 = 4*t1175*t129
      t8690 = 16*t1178*t129
      t8693 = 16*t195*t437*t7723
      t8696 = 32*t195*Mzsq*t2452
      t8699 = 8*t48*t437*t7723
      t8700 = t54*CsThr
      t8702 = 8*t907*t8700
      t8704 = 2*t1175*t125
      t8706 = 2*t1178*t125
      t8707 = t157*t899
      t8709 = 4*t537*t129
      t8710 = t157*t496
      t8711 = 4*t5002
      t8712 = 16*t5009
      t8713 = 16*t5011
      t8714 = 32*t5013
      t8715 = -t8681-t8683+t8686-t8688+t8690+t8693-t8696-t8699-t8702-t87
     #04+t8706-t8707-t8709-t8710+t8711+t8712+t8713-t8714
      t8717 = 2*t157*t1689
      t8719 = 8*t7154*t665
      t8720 = t300*t665
      t8721 = 32*t8720
      t8723 = 4*t335*t1755
      t8725 = 64*t747*t456
      t8727 = 8*t456*t2597
      t8728 = t300*t184
      t8729 = 64*t8728
      t8731 = 32*t747*t41
      t8733 = 32*t747*t174
      t8735 = 128*t6088*t202
      t8737 = 16*t7154*t184
      t8739 = 64*t6092*t202
      t8741 = 64*t747*t452
      t8743 = 4*t1171*t2597
      t8744 = t300*t438
      t8745 = 32*t8744
      t8747 = 64*t6092*t447
      t8749 = 32*t6092*t443
      t8751 = 64*t6088*t887
      t8752 = t8717+t8719-t8721+t8723+t8725-t8727+t8729-t8731-t8733-t873
     #5-t8737+t8739+t8741+t8743-t8745+t8747-t8749+t8751
      t8755 = 128*t6088*t447
      t8757 = 8*t7154*t438
      t8759 = 32*t6092*t887
      t8761 = 64*t6088*t443
      t8762 = t1171*t2140
      t8764 = 8*t177*t826
      t8765 = t300*t511
      t8766 = 32*t8765
      t8768 = 8*t527*t520
      t8770 = 32*t8639*t2556
      t8772 = 32*t7451*t2556
      t8774 = 16*t195*t2368
      t8776 = 4*t456*t1964
      t8778 = 64*t7451*t7723
      t8780 = 32*t7451*t1989
      t8782 = 2*t456*t2140
      t8784 = 32*t195*t7872
      t8785 = t335*t1244
      t8787 = 16*t80*t537
      t8788 = -t8755+t8757-t8759+t8761-t8762+t8764-t8766-t8768+t8770-t87
     #72+t8774+t8776+t8778-t8780+t8782-t8784-t8785+t8787
      t8790 = 16*t80*t534
      t8792 = 32*t80*t462
      t8794 = 4*t1368*t8700
      t8796 = 8*t495*t8700
      t8797 = 32*t827
      t8798 = 32*t830
      t8799 = 64*t832
      t8800 = t8790-t8792+t8794-t8796+t791-t800-t803-t804+t807-t810+t813
     #-t8797-t8798+t8799+t890-t900-t903+t906-t909
      t8807 = -t911-t914-t917+64*t919-32*t922-32*t925-8*t5087+t5110+t512
     #3+t1240+t1242+t1246-t1312+t1315-t1317+t1319-t1321+t1323
      t8808 = t1343+t1345-t1346-t1355+t1358+t1360-t1362-t1365+t1367+t137
     #0-t1371+t1330+t1328+t1400+t1401-t1402-t1403+t1404
      t8809 = t8807+t8808
      t8810 = 16*t5199
      t8811 = 16*t5201
      t8812 = 32*t5203
      t8813 = 4*t5215
      t8814 = 8*t5217
      t8816 = 64*t8639*t7723
      t8818 = 32*t8639*t1989
      t8820 = 16*t527*t205
      t8821 = t1405-t1406-t1474-t1475-t1476+t1477+t1478-t1479-t1480+t148
     #1+t8810+t8811-t8812+t8813-t8814-t8816+t8818+t8820
      t8824 = 2*t156*t5615*t32
      t8826 = 16*t195*t2392
      t8828 = 32*t7451*t2452
      t8830 = 32*t7451*t2142
      t8832 = 64*t8639*t2455
      t8834 = 32*t8639*t2142
      t8835 = t300*t444
      t8836 = 32*t8835
      t8837 = t300*t448
      t8838 = 64*t8837
      t8840 = 8*t465*t439
      t8842 = 8*t7154*t444
      t8844 = 16*t7154*t448
      t8846 = 4*t1171*t507
      t8848 = 8*t527*t8460
      t8850 = 32*t8639*t2452
      t8852 = 64*t7451*t2455
      t8854 = 16*t80*t1172
      t8856 = 16*t80*t1175
      t8858 = 32*t80*t1178
      t8860 = 4*t898*t8700
      t8861 = t8824+t8826-t8828-t8830-t8832+t8834-t8836+t8838-t8840+t884
     #2-t8844+t8846-t8848+t8850+t8852+t8854+t8856-t8858+t8860
      t8866 = t7481*t1712
      t8869 = t7897*t7176*t7177
      A1e0c6p1t1n7 = 12*(t8715+t8752+t8788+t8800+t8809+t8821+t8861)*t862
     #0*t8866*t8869
      A1e0c7p1t1n7 = 0
      t8872 = 2*t537*t4969
      t8874 = 2*t157*t444
      t8876 = 16*t300*t887
      t8878 = 8*t7154*t887
      t8880 = 8*t45*t184
      t8882 = 16*t45*t2067
      t8884 = 32*t195*t2367
      t8886 = 32*t300*t202
      t8888 = 64*t300*t2071
      t8890 = 128*t7451*t1756
      t8892 = 16*t7154*t202
      t8894 = 64*t7154*t2071
      t8896 = 64*t8639*t1756
      t8898 = 32*t1313*t186
      t8900 = 32*t465*t57
      t8902 = 32*t485*t190
      t8904 = 8*t1313*t174
      t8905 = -t8872-t8874-t8876+t8878-t8880+t8882-t8884+t8886-t8888+t88
     #90-t8892+t8894-t8896-t8898+t8900+t8902+t8904
      t8907 = 64*t485*t2050
      t8909 = 16*t1313*t456
      t8911 = 32*t1171*t57
      t8913 = 32*t8639*t2121
      t8915 = 32*t7154*t1929
      t8917 = 64*t7451*t2121
      t8918 = 8*t1843
      t8920 = 8*t465*t129
      t8922 = 4*t45*t444
      t8924 = 4*t156*t7872
      t8926 = 16*t195*t7723
      t8928 = 16*t48*t7723
      t8931 = 8*t1142*t174*Er
      t8933 = 16*t1142*t2054
      t8936 = 8*t1142*t456*Er
      t8938 = 4*t465*t203
      t8940 = 8*t1171*t129
      t8941 = -t8907-t8909+t8911+t8913-t8915-t8917+t8918-t8920+t8922+t89
     #24-t8926+t8928-t8931+t8933+t8936+t8938-t8940
      t8944 = 16*t195*t2452
      t8946 = 8*t48*t2142
      t8948 = 8*t195*t2142
      t8949 = 8*t1845
      t8950 = 4*t2042
      t8951 = 4*t2075
      t8952 = 16*t1937
      t8953 = 16*t1943
      t8954 = 4*t2103
      t8955 = 4*t2107
      t8956 = t8944-t8946+t8948-t8949+t8950+t8951+t2099+t5270+t5271-t193
     #5+t8952+t8953-t1946+t2102-t8954-t8955-t2110
      t8960 = -t2112-32*t2113-t1954+t1956+t2128+t901+t2130-t1960-t2133+3
     #2*t2134+t2137-t7967+t1968+t1969+t801-16*t2138-t7969-t2141
      t8965 = 1/(t32-t427-t443+t447)
      A1e0c8p1t1n7 = -12*t8620*(t8905+t8941+t8956+t8960)*t138*t8965*t664
     #9*t8652
      t8972 = 3*t157*t8327
      t8974 = 6*t157*t4702
      t8977 = 3*t157*t5825*t32
      t8979 = 2*t157*t8333
      t8982 = 2*t157*t5821*t32
      t8984 = 2*t157*t2373
      t8985 = t157*t2377
      t8986 = t157*t2380
      t8987 = t2604*t32
      t8989 = 8*t8987*t520
      t8990 = t3489*t40
      t8992 = 32*t8990*t2458
      t8994 = 96*t8990*t2392
      t8995 = t8637*t904
      t8996 = 128*t8995
      t8997 = t8637*t2409
      t8998 = 32*t8997
      t8999 = t2264*CsThr
      t9001 = 32*t8999*t2458
      t9002 = t2264*t40
      t9004 = 32*t9002*t2458
      t9005 = t2305*CsThr
      t9007 = 16*t9005*t2458
      t9009 = 3*t157*t8169
      t9011 = 6*t157*t8172
      t9015 = 3*t157*Mzsq*t2365*t32
      t9017 = 2*t157*t8175
      t9018 = t8972+t8974+t8977-t8979-t8982-t8984+t8985+t8986+t8989-t899
     #2-t8994-t8996+t8998+t9001-t9004-t9007+t9009+t9011+t9015-t9017
      t9022 = 2*t157*t2163*t141*t32
      t9023 = t2979*t40
      t9025 = 96*t9023*t2452
      t9027 = -t2391+t2394-t2398+t2420-t2428+t2430+t2432-t2435+t2440-t24
     #43+t2446
      t9032 = -t2536-t2538+t2540-t2544+t2546+t2552+t2555-t2558-t2560-t25
     #63-128*t2583
      t9033 = -t2449-t2454+t2457-t2460+t2462-t2464-t2466+t2468-t2470+t24
     #72+t9032
      t9034 = 96*t2619
      t9035 = 96*t2628
      t9036 = 64*t2655
      t9037 = 24*t2858
      t9038 = 8*t2964
      t9039 = 24*t3026
      t9041 = 12*t1325*t3664
      t9045 = 8*t769*t2604*t1
      t9046 = t2503*t1178
      t9047 = 64*t9046
      t9048 = t2418*t3656
      t9049 = 32*t9048
      t9051 = 16*t2503*t918
      t9053 = 8*t2418*t3579
      t9054 = t2418*t3661
      t9055 = 32*t9054
      t9057 = 128*t2499*t1178
      t9059 = 32*t751*t3579
      t9060 = t2310*t528
      t9061 = 96*t9060
      t9063 = 12*t2248*t766
      t9064 = t3489*CsThr
      t9066 = 128*t9064*t8109
      t9067 = t9045-t9047+t9049+t9051-t9053+t9055+t9057+t9059+t9061-t906
     #3-t9066
      t9072 = 128*t8990*t8109
      t9073 = t152*t40
      t9075 = 8*t9073*t8120
      t9077 = 16*t8684*t8105
      t9079 = 8*t9073*t8105
      t9081 = 2*t157*t8139
      t9083 = t157*t2703*t495
      t9085 = 6*t157*t8131
      t9087 = 6*t157*t2799
      t9090 = 2*t157*t447*t5829
      t9092 = 2*t157*t8136
      t9094 = t157*t2727*t495
      t9095 = 16*t3223
      t9098 = 8*t45*t141*t904
      t9100 = 32*t9002*t2447
      t9102 = 16*t2134*t3661
      t9103 = t2979*CsThr
      t9105 = 96*t9103*t2556
      t9107 = 96*t8990*t2368
      t9109 = 96*t9064*t2368
      t9110 = t9072+t9075-t9077+t9079+t9081+t9083-t9085-t9087-t9090+t909
     #2+t9094-t5466-t5467+t9095+t9098-t9100-t9102+t9105-t9107+t9109
      t9112 = 64*t9064*t2553
      t9113 = t2433*CsThr
      t9115 = 48*t9113*t2368
      t9116 = t2433*t40
      t9118 = 24*t9116*t2368
      t9119 = t2604*sqq
      t9121 = 32*t1914*t9119
      t9122 = t2606*t32
      t9124 = 8*t9122*t8460
      t9126 = 8*t8631*t2405
      t9128 = 64*t751*t3569
      t9130 = 32*t751*t3572
      t9132 = 4*t335*t3681
      t9136 = 8*t1755*t1324*CsThr*t437
      t9139 = t335*t437
      t9141 = 4*t1*t2163*t9139
      t9143 = 16*t5317*t9139
      t9145 = 4*t8522*t1255
      t9148 = 32*t3747*CsThr*t2455
      t9150 = 16*t9005*t2447
      t9153 = 64*t2256*CsThr*t2455
      t9156 = 64*t2256*t40*t2455
      t9157 = 48*t5579
      t9158 = 12*t5581
      t9160 = 128*t6092*t1356
      t9162 = 16*t2248*t2512
      t9163 = -t9141+t9143-t9145+t9148-t9150-t9153+t9156-t9157-t9158-t91
     #60+t9162
      t9166 = t2310*t335
      t9167 = 128*t9166
      t9169 = 256*t3599*t527
      t9171 = 32*t2310*t499
      t9172 = 64*t5595
      t9174 = 32*t9113*t8109
      t9176 = 64*t8684*t8102
      t9178 = 16*t8684*t8120
      t9180 = 64*t1914*t1647
      t9181 = t1175*t1309
      t9182 = 24*t9181
      t9184 = 4*t4645*t1925
      t9187 = 4*t4563*t1930
      t9189 = 16*t1914*t3661
      t9190 = t2340*t2195
      t9191 = t49*CsThr
      t9194 = 4*t9190*t9191*Mzsq
      t9196 = 96*t9103*t2452
      t9199 = 128*t751*t40*t448
      t9201 = 64*t9103*t2455
      t9203 = 16*t9139*t520
      t9205 = 16*t2134*t3656
      t9207 = 32*t8999*t2447
      t9208 = 64*t3563
      t9209 = t9187-t9189-t9194+t9196-t9199-t9201-t9203-t9205+t9207-t920
     #8+t3568
      t9219 = -t3608+t3609+t3610+t3685-48*t5648+8*t5654+32*t5656-12*t566
     #3-48*t5677+32*t5679+8*t5681
      t9220 = -t3571+t3574+t3576-t3578+t3581-t3583+t3585+t3587+t3595-t35
     #97+t9219
      t9232 = 8*t5683-4*t5685-48*t5691-4*t5693-16*t5695-16*t5703-16*t570
     #5+32*t5709-t3796-t3798+t3800-t3802+t3804+t3806+t3811-t3813+t3815-t
     #3816-t3817+t3819
      t9236 = t3848+t3850-t3852+t3854+t3856+t3858-t3860+t3862+t3864-t386
     #6-t3868
      t9238 = t9232+32*t5720-16*t5722-t3827-t3833+t3835-t3837-t3839+t384
     #2+t3844-t3846+t9236
      t9239 = 128*t4039
      t9240 = t9191*t32
      t9242 = 16*t3385*t9240
      t9245 = 64*t2433*t5277*t32
      t9247 = 96*t6092*t496
      t9250 = 96*t751*t912
      t9252 = 96*t3599*t924
      t9254 = 64*t3599*t918
      t9256 = 64*t6092*t904
      t9258 = 12*t45*t4702
      t9260 = 64*t2253*t3656
      t9262 = 64*t2253*t3661
      t9264 = 32*t2253*t3579
      t9266 = 32*t2253*t3572
      t9268 = 192*t3599*t1172
      t9270 = 128*t3599*t1178
      t9271 = t9250+t9252-t9254-t9256-t9258-t9260-t9262+t9264+t9266-t926
     #8+t9270
      t9274 = 24*t2310*t912
      t9275 = t2310*t1255
      t9276 = 96*t9275
      t9277 = t177*t9139
      t9278 = 96*t9277
      t9280 = 192*t751*t1255
      t9282 = 128*t751*t9119
      t9284 = 64*t751*t3656
      t9285 = t177*t8987
      t9286 = 64*t9285
      t9288 = 12*t456*t4545
      t9290 = 8*t456*t4548
      t9293 = 96*t416*t6092*t32
      t9296 = 64*t2164*t7390
      t9298 = 4*t156*t8400
      t9302 = 4*t156*Er*t2365*t122
      t9304 = 24*t177*t894
      t9306 = 16*t177*t3611
      t9308 = 8*t2418*t3572
      t9310 = 64*t751*t3661
      t9312 = 32*t335*t2625
      t9314 = 24*t177*t500
      t9316 = 16*t177*t2642
      t9318 = t177*t465*t122
      t9319 = 96*t9318
      t9320 = -t9296-t9298-t9302-t9304+t9306-t9308-t9310-t9312-t9314+t93
     #16+t9319
      t9324 = t177*t9122
      t9325 = 64*t9324
      t9327 = 192*t751*t528
      t9328 = t2606*sqq
      t9330 = 128*t751*t9328
      t9332 = 32*t8684*t2396
      t9334 = 32*t9064*t2458
      t9336 = 32*t8631*t904
      t9338 = 8*t8631*t2409
      t9340 = 96*t9064*t2392
      t9342 = -t4256+t4258-t4385-t4388+t4390+t4392-t4393+t4395-t4397-t43
     #99-t4400
      t9344 = t8637*t2405
      t9345 = 32*t9344
      t9348 = 256*t751*CsThr*t448
      t9349 = t2310*CsThr
      t9350 = t9349*t448
      t9351 = 64*t9350
      t9353 = 64*t9064*t2444
      t9355 = 32*t9064*t2447
      t9357 = 48*t8684*t2426
      t9359 = 192*t3599*t534
      t9361 = 24*t2310*t798
      t9363 = 128*t2979*t2051
      t9365 = 128*t2979*t2561
      t9369 = 64*t2433*t2050*t122
      t9370 = 16*t5898
      t9371 = 64*t5926
      t9372 = t9369+t9370+t9371+t4592-t4602+t4617-t4619-t4632+t4634-t463
     #6-t4638
      t9382 = -t4640+t4642+96*t4654-64*t4657+96*t4659-64*t4661-t5943+32*
     #t4689+32*t4696+32*t4785
      t9383 = 32*t4793
      t9385 = 64*t9023*t2455
      t9388 = 16*t2310*t40*t448
      t9390 = 48*t9113*t2392
      t9392 = 24*t9116*t2392
      t9393 = t2340*t437
      t9395 = 12*t9393*t9240
      t9397 = 16*t1914*t3656
      t9399 = 64*t8990*t2444
      t9401 = 32*t8990*t2447
      t9403 = 96*t9023*t2556
      t9405 = 48*t8684*t2351
      t9406 = t9383+t9385+t9388-t9390+t9392+t9395-t9397+t9399-t9401-t940
     #3-t9405
      t9411 = 48*t1914*t528
      t9413 = 32*t1914*t9328
      t9415 = 12*t2663*t534
      t9417 = 48*t2681*t534
      t9420 = 12*t2340*Mzsq*t9240
      t9421 = t2340*t2163
      t9422 = t9191*sqq
      t9424 = 8*t9421*t9422
      t9426 = 32*t2681*t1178
      t9429 = 4*t9421*t9191*t437
      t9433 = 32*t8684*t2362-t9411+t9413-t9415-t9417+t9420+t9424+t9426-t
     #9429+24*t456*t122*t520
      t9436 = 16*t457*t54*t2163
      t9438 = 48*t1914*t1255
      t9440 = 12*t2663*t1172
      t9442 = 8*t2663*t1178
      t9444 = 48*t2681*t1172
      t9446 = 8*t9190*t9422
      t9448 = 64*t8990*t2553
      t9450 = 4*t457*t3681
      t9452 = 24*t1175*t2897
      t9454 = 16*t4645*t2138
      t9456 = 32*t2840*t1178
      t9457 = -t9436-t9438-t9440+t9442-t9444+t9446+t9448-t9450+t9452-t94
     #54+t9456
      t9466 = t1980*t404
      t9467 = t9466*t8653
      A1e0c9p1t1n7 = -12*(t9027-t3872+t9018-t9022-t9025-t9115+t9118+t912
     #1+t9124-t9126-t9128+t9130-t9132+t9136-t9112-t4128-t2599-t2601-t238
     #7+t2383+t2381+t2378-t2375+t3870+t2370+t9247+t4119+t3876+t3878-t388
     #0+t3882-t5391+t9242+t9245+t9372+t9382-t9357-t9359-t9361-t9363+t936
     #5-t9351-t9353+t9355+t9340+t9342+t9345+t9348+t9334+t9336-t9338+t923
     #8-t9239-t9327+t9330+t9332+t9320-t9325+t9290+t9293+t9220-t9284-t928
     #6-t9288+t9278-t9280+t9282+t9276+t9209+t9271-t9274+t9171+t9172+t917
     #4+t9176-t9178+t9180-t9182+t9184-t2364-t9167+t9169+t2255+t9163+t911
     #0+t9067+t9038-t9039-t9041+t9035-t9036-t9037+t9457+t9033+t9034+t943
     #3+t9406)*t8620*t402*t7406*t9467
      A1e0c1p1t2n7 = 0
      A1e0c2p1t2n7 = A1e0c2p1t1n7
      A1e0c3p1t2n7 = A1e0c3p1t1n7
      A1e0c4p1t2n7 = 0
      A1e0c5p1t2n7 = A1e0c5p1t1n7
      A1e0c6p1t2n7 = A1e0c6p1t1n7
      A1e0c7p1t2n7 = 0
      A1e0c8p1t2n7 = A1e0c8p1t1n7
      A1e0c9p1t2n7 = A1e0c9p1t1n7
      A2e0c1p1t1n7 = 0
      t9471 = 2*t456*t118
      t9473 = 4*t45*sll
      t9475 = 8*t195*Er
      A2e0c2p1t1n7 = -48*t8670*sqq*(t202+t4970-t4972-t9471-t4971+t9473+t
     #9475)*t406*t8675*t7175*t7177
      t9484 = t8683-t8688+t8702+t8704-t8706+t8707-t8709+t8710-t8717-t502
     #4+t5026+t5028
      t9485 = 16*t5029
      t9486 = 16*t5031
      t9487 = 32*t5034
      t9488 = 16*t8720
      t9489 = 32*t8728
      t9490 = 16*t8744
      t9491 = 16*t8765
      t9492 = -t9485-t9486+t9487+t9488-t8723+t8727-t9489-t8743+t9490+t87
     #62+t9491-t8772-t8776
      t9494 = t8778-t8780-t8782+t8785-t8787-t8790+t8792-t8794+t8796-t505
     #3+t803+t804+t5054
      t9499 = -t5058+t900+t903-t906+t5059+t5060-16*t940-4*t5089-4*t5091+
     #t5093+t5096-t5098+8*t5099
      t9507 = 32*t5101-16*t5103-16*t5105-t1246+8*t1309+t1312-t5145+t1317
     #-t5146+t5147-t1323+8*t1340
      t9508 = 16*t1353
      t9510 = 64*t300*t8442
      t9512 = 32*t465*t2903
      t9514 = 32*t465*t514
      t9516 = 32*t300*t8533
      t9518 = 32*t300*t2008
      t9519 = -t5149+t1346-t9508+t1355-t1358-t1360+t5152-t5153-t9510+t95
     #12+t9514+t9516+t9518
      t9522 = 64*t465*t5268
      t9523 = -t9522+t1371-t1330-t1328-t1401-t1404-t5190+t1474+t1475+t14
     #76-t5195+t5196-t1481
      t9524 = 8*t5213
      t9525 = 16*t8835
      t9526 = 32*t8837
      t9527 = t9524-t8824-t8828-t8830+t9525-t9526+t8840-t8846+t8852-t885
     #4-t8856+t8858-t8860
      A2e0c3p1t1n7 = -48*(t9484+t9492+t9494+t9499+t9507+t9519+t9523+t952
     #7)*t8620*t8866*t8869
      A2e0c4p1t1n7 = 0
      t9536 = 2*t465*t2875
      t9538 = 4*t1755*t456
      t9540 = 8*t465*t3
      t9543 = 4*t45*t887
      t9545 = 8*t1171*t3
      t9546 = 4*t8633
      t9547 = 8*t8646
      t9549 = 2*t1171*t2875
      t9550 = t9543+t444-t5312+t9545-t5323-t5322-t5320+t9546-t8645-t9547
     #-t9549
      A2e0c5p1t1n7 = -12*t8620*(t511+t5295-t5316-t9536-t5318+t9538+t9540
     #+t132+t5314+t5297+t9550)*t8965*t8653
      t9558 = -t8972-t8974-t8977+t8979+t8982+t8984-t8985-t8986+64*t8995-
     #16*t8997+t9001-t9009-t9011-t9015
      t9559 = 8*t2399
      t9560 = 8*t2412
      t9561 = 16*t2414
      t9562 = t9017+t9022+t2364-t2370+t2375-t2378-t2381-t5346+t5347+t534
     #9-t2394+t2398-t9559-t9560+t9561
      t9564 = -t5353+t5354-t5355-t2432+t5356-t2440+t2443-t5357+t2464+t53
     #58-t2468+t5361-t2472-t5362
      t9565 = 16*t2486
      t9566 = 24*t2491
      t9567 = 16*t2565
      t9568 = 24*t2579
      t9569 = 64*t5393
      t9570 = t9565-t9566+t5367+t5368-t2540+t2544-t5370-t2552+t9567-t537
     #2-t9568+t2599-t9569+t9041-t9045
      t9577 = 32*t9046-16*t9048-16*t9054-48*t9060+t9063-t9066-t9081-t908
     #3+t9085+t9087+t9090-t9092-t9094-t9098
      t9580 = t9102+t9105+t9109-t9112-t9121+t9132-t9136+t9141-t9143+t914
     #5-t9153-t9162+64*t9166-t9180-8*t9181
      t9581 = t9577+t9580
      t9582 = t9184+t9187+t5617+t9189+t9194+t9196-t9201+t9205+t9207-t563
     #4+t3583-t3585-t3587-t5636
      t9583 = 64*t5724
      t9585 = 96*t9349*t8554
      t9586 = t3597+t3608-t3609-t3610+t5672+t5717-t3800+t5718-t3804-t381
     #5+t3816+t3817-t3819-t9583-t9585
      t9591 = 32*t1882*t8499
      t9593 = 32*t1882*t8503
      t9595 = 64*t9349*t8457
      t9597 = 32*t8637*t8499
      t9599 = 96*t8637*t8491
      t9600 = 32*t3828
      t9601 = -t9591-t9593+t9595-t9597-t9599+t3827+t9600-t3835+t3837+t57
     #59-t5760-t3844+t5761-t5762
      t9602 = 48*t9275
      t9603 = 48*t9277
      t9604 = 32*t9285
      t9605 = t3852-t5763+t5764+t3868-t3870+t5765-t5766-t3878+t3880-t924
     #2-t9245+t9258-t9602-t9603+t9604
      t9607 = 48*t9318
      t9608 = 32*t9324
      t9610 = 128*t2310*t9240
      t9611 = 16*t5814
      t9612 = 16*t5816
      t9613 = 4*t5819
      t9614 = 4*t5823
      t9615 = 8*t5827
      t9616 = t9288-t9290+t9298+t9302-t9607+t9608+t9610-t4119+t4128+t961
     #1+t9612+t9613+t9614-t9615
      t9617 = 12*t5831
      t9618 = 8*t5834
      t9619 = 32*t5840
      t9620 = 32*t5842
      t9621 = 48*t5846
      t9622 = 48*t5848
      t9623 = 48*t5850
      t9624 = 32*t5852
      t9625 = t9617-t9618-t5837+t5839-t9619-t9620-t5845+t9621+t9622+t962
     #3-t9624-t5857+t5860+t9334+t9340
      t9630 = 16*t9393*t9191*t427
      t9632 = 128*t8637*t8572
      t9633 = 16*t5876
      t9634 = 16*t9344
      t9635 = -t4258-t9630+t9632-t9633+t4385+t4388-t4390-t5883+t4393-t43
     #95+t4397+t4399+t4400-t9634
      t9636 = 32*t9350
      t9637 = 8*t5891
      t9638 = 12*t5894
      t9639 = 16*t5911
      t9640 = 48*t5913
      t9641 = 32*t5915
      t9642 = 16*t5917
      t9643 = t9636-t9353+t9355-t9363-t9637+t9638-t5901+t5904-t5906-t590
     #8+t5910+t9639+t9640-t9641+t9642
      t9645 = 24*t4650
      t9646 = -t5921-t5923+t5939+t5941+t9645-t5944-t5945+t9395+t9397+t94
     #11-t9413+t9415+t9417+t9420-t9424
      t9659 = -t9426+t9429+t9438+t9440-t9442+t9444-t9446-96*t9349*t8533+
     #64*t8637*t8536-32*t8637*t8503+64*t7386*t8457-96*t8637*t8527+64*t86
     #37*t8530+t9450-t9456
      A2e0c6p1t1n7 = 12*(t9558+t9562+t9564+t9570+t9581+t9582+t9586+t9601
     #+t9605+t9616+t9625+t9635+t9643+t9646+t9659)*t8620*t402*t7406*t9467
      A2e0c1p1t2n7 = 0
      A2e0c2p1t2n7 = A2e0c2p1t1n7
      A2e0c3p1t2n7 = A2e0c3p1t1n7
      A2e0c4p1t2n7 = 0
      A2e0c5p1t2n7 = A2e0c5p1t1n7
      A2e0c6p1t2n7 = A2e0c6p1t1n7
      t9668 = -j4+t6644+t7664+t7668+t7669
      t9670 = -t7678+t7680-t169-t7681+t7683+t7684-t7685-t132+t82+t7686+t
     #7687
      t9671 = t6662-t122-t6664-t7671+t77-t78-t7672+t126+t7673-t7675+t967
     #0
      t9674 = t7191*t6654
      A1e0c1p1t1n8 = -12*t9668*t9671*t142*t406*t9674
      A1e0c2p1t1n8 = 0
      t9678 = t6662+t7696+t234-t122-t6664-t6668+t322-t323+t7697+t207+t76
     #98-t7678+t7680-t169-t208-t7699+t6660+t6666-t7700-t7687
      t9683 = t7489*t6654
      t9684 = t402*t140*t9683
      A1e0c3p1t1n8 = -24*t9668*t9678*t7175*t404*t9684
      t9688 = t9668*t414
      t9692 = t420*t6654*t7489*t140
      A1e0c4p1t1n8 = 12*t9671*Mzsq*t9688*t9692
      A1e0c5p1t1n8 = 0
      t9694 = t7720+t7722-t7725+t7727-t7729-t7731+t7733-t7735+t7737+t773
     #9+t7741+t7743-t7745-t7747-t7749
      t9695 = t7752-t7754-t7756-t7758-t7760+t7762+t7764-t7766-t7768+t777
     #0-t7772-t7774+t7776+t7778-t7780
      t9697 = t7784+t791-t803+t5054-t7785-t7786-t7787+t7788+t7789+t890-t
     #903+t7790-t7791-t7792-t7793
      t9700 = t1481-t7823+t7826-t7829+t7832-t7835-t7837+t7841-t7844+t784
     #6-t7848+t7850+t7852-t7854+t7856
      t9701 = t7859+t7861-t7863-t7865+t7867+t7869-t7871+t7874-t7876-t787
     #8-t7880+t7882+t7884-t7886+t7888
      t9706 = t7495*t1712
      t9708 = t7897*t9683
      A1e0c6p1t1n8 = -12*(t9694+t9695+t9697+t7803+t7821+t9700+t9701)*t96
     #68*t9706*t9708
      t9711 = -t7901+t6974-t6976-t6978+t6980-t7902+t7904+t7906-t7908-t79
     #10+t7912+t7914-t7917+t7919-t7921
      t9712 = -t7924+t7926-t7928-t7930-t7932+t7934+t7936-t7938+t7939-t79
     #41+t7943+t7945+t7947+t7949-t7951+t7952
      t9714 = -t7956-t7958+t7959+t7960-t2037-t2039+t7961+t1923-t7962-t79
     #63+t7964-t1935+t1940-t1942-t1946+t1948
      A1e0c7p1t1n8 = 12*(t9711+t9712+t9714+t7970)*t9668*t1980*t9684
      A1e0c8p1t1n8 = 0
      t9720 = -t8020+t8021+t8022-t8023-t8024+t8026-t8029+t8030-t8031-t80
     #32-t8033+t8034-t8036+t8039-t8041+t8044+t8048+t8050
      t9723 = -t8055-t8057-t8059-t8061-t8064+t8066+t8068-t8070-t8072+t80
     #74+t8076+t8078+t8081-t8084+t8086-t8088+t8090
      t9724 = -t8093-t8096-t8098+t8101-t8104+t8107+t8111-t8114+t8116+t81
     #19+t8122+t8124-t8127-t8130+t8133+t8135-t8138-t8141
      t9726 = -t8145-t8147+t8149-t8151+t8153+t8155-t8157+t8159-t8160+t81
     #61-t8162+t3585+t3587+t8163-t8164+t8165+t8167
      t9727 = -t8171-t8174+t8177-t8180+t8182-t8184-t8187+t8190-t8192-t81
     #94-t8197-t8199+t8201+t8203-t8205-t8207-t8209+t8211
      t9731 = -t8217-t8219-t8221-t8223+t8225-t8227-t8229+t8231-t8233+t82
     #35-t8237-t8239+t8241+t8243-t8246+t8248+t8250
      t9732 = t8253-t8256+t8259+t8262-t8265+t8268+t8271+t8273+t8275+t827
     #8-t8281-t3796-t3798+t3800-t3802+t3804+t3806+t3811
      t9734 = t3856+t3858-t3860+t3862+t3864-t3866-t3868+t3870-t3872+t387
     #6+t3878-t3880+t3882+t8290-t8292+t8295+t8297+t8299
      t9737 = t4119-t4128-t4256+t4258-t8304-t8306-t8308+t8310-t8312-t831
     #4+t8316-t8318+t8320+t8322-t8324+t8326-t8329
      t9738 = -t8332+t8335+t8337-t8339-t8341-t4385-t4388+t4390+t4392-t43
     #93+t4395-t4397-t4399-t4400-t4632-t4638+t8342-t8343
      t9740 = t8346-t8347-t8348-t5944+t8349-t5945-t8350+t8351+t8352+t835
     #3-t8354+t8355+t8356-t8357+t8358-t8359+t8361
      t9741 = -t8365+t8368-t8370+t8372+t8374-t8377-t8379+t8382+t8384-t83
     #86+t8388+t8390-t8392-t8394+t8396-t8399+t8402+t8404
      t9748 = t4936*t7175
      t9751 = t4940*t6649*t9674
      A1e0c9p1t1n8 = 12*(t8006+t8019+t9720+t9723+t9724+t9726+t9727+t9731
     #+t9732+t8288+t9734+t9737+t9738+t9740+t9741)*t9668*t432*t9748*t9751
      A1e0c1p1t2n8 = A1e0c1p1t1n8
      A1e0c2p1t2n8 = 0
      A1e0c3p1t2n8 = A1e0c3p1t1n8
      A1e0c4p1t2n8 = A1e0c4p1t1n8
      A1e0c5p1t2n8 = 0
      A1e0c6p1t2n8 = A1e0c6p1t1n8
      A1e0c7p1t2n8 = A1e0c7p1t1n8
      A1e0c8p1t2n8 = 0
      A1e0c9p1t2n8 = A1e0c9p1t1n8
      t9753 = -t8418+t6662-t4950-t7678+t8417-t4949+t7681-t8419-t122+t494
     #7+t126+t77
      A2e0c1p1t1n8 = 48*t9753*Mzsq*t9688*t9692
      A2e0c2p1t1n8 = 0
      t9757 = t7720+t5024-t5026-t5028-t7722+t7725-t8426+t7733+t7735-t773
     #7-t7739+t7741-t8427-t7760+t8428-t7780+t8429+t5053
      t9760 = -t5152+t5153-t1371+t1330+t1328+t1401+t1404+t5190+t8434-t84
     #36-t8438-t1474-t1475-t1476+t5195-t5196+t1481-t7823
      A2e0c3p1t1n8 = -48*(t9757-t8431+t9760+t8463)*t9668*t9706*t9708
      A2e0c4p1t1n8 = 12*(t8482-t8483)*t9668*t1980*t9684
      A2e0c5p1t1n8 = 0
      t9776 = -t8008+t5370+t8010-t2599-t8524+t8526-t8529+t8532-t8535+t85
     #38-t8540
      t9780 = t3587+t5636-t8164-t8171-t8174+t8177-t8182+t8184-t8187+t819
     #0-t8194
      t9784 = t8281-t5672-t5717+t3800-t5718+t3804+t3815-t3816-t3817+t381
     #9+t3835
      t9787 = t3870-t5765+t5766+t3878-t3880-t8593+t8594+t4119-t4128+t583
     #7-t5839
      t9791 = t8335+t8337-t8339-t8341-t4385-t4388+t4390+t5883-t4393+t439
     #5-t4397
      t9794 = -t5941-t8361-t8602+t8603+t8604-t8605+t8390-t8394+t8402+t86
     #07-t8609
      A2e0c6p1t1n8 = 12*(t4258-t8332-t8329+t8324-t5860+t5857+t5845-t8316
     #+t8314+t8304-t5617+t3585+t5634+t9776+t9780+t9784+t9787+t9791+t9794
     #-t8601-t8501-t8592+t8574+t8575-t8579+t8589+t8497+t8568+t8571-t8493
     #-t8162+t8155-t8145+t8490-t8521-t8519+t8515-t8517-t8505-t8508+t8511
     #-t7980)*t9668*t432*t9748*t9751
      A2e0c1p1t2n8 = A2e0c1p1t1n8
      A2e0c2p1t2n8 = 0
      A2e0c3p1t2n8 = A2e0c3p1t1n8
      A2e0c4p1t2n8 = A2e0c4p1t1n8
      A2e0c5p1t2n8 = 0
      A2e0c6p1t2n8 = A2e0c6p1t1n8
      A1e0c1p1t1n9 = 0
      t9803 = -j3+t6644+t8619
      t9804 = -t8621-t8622+t8623-t8624-t8626+t183-t204+t8627+t207-t7673+
     #t8629
      t9805 = t4950-t132+t82+t184+t8632+t8634-t8636+t8638+t8641+t8643-t8
     #645-t8647
      t9806 = t9804+t9805
      t9810 = t406*t6672*t7177
      A1e0c2p1t1n9 = 12*t9803*t9806*t215*t9810
      t9812 = -t8629+t8645-t8655+t8657+t234+t7183-t8627+t295+t307+t319+t
     #7182+t322-t323
      t9813 = -t7185+t7186-t7187-t8641-t8659+t8660-t122+t8621-t8623-t372
     #+t125-t377-t8633
      t9819 = t7705*t7150*t7177
      A1e0c3p1t1n9 = -24*t9803*(t9812+t9813)*t6652*t402*t9819
      A1e0c4p1t1n9 = 0
      t9827 = t429*t7177*t7150*t138
      A1e0c5p1t1n9 = -12*t9806*Mzsq*t9803*t432*t9827
      t9830 = -t8681+t8683-t8686+t8688-t8690-t8693+t8696-t8699+t8702+t87
     #04-t8706+t8707+t8709+t8710+t8711+t8712+t8713-t8714
      t9831 = -t8717+t8719+t8721-t8723-t8725+t8727-t8729-t8731-t8733+t87
     #35-t8737+t8739-t8741-t8743+t8745+t8747-t8749-t8751
      t9833 = t8755+t8757-t8759-t8761+t8762+t8764+t8766+t8768+t8770+t877
     #2-t8774-t8776-t8778+t8780-t8782+t8784+t8785-t8787
      t9834 = -t8790+t8792-t8794+t8796+t791-t800-t803-t804+t807-t810+t81
     #3-t8797-t8798+t8799+t890-t900-t903+t906-t909
      t9837 = t1405-t1406-t1474-t1475-t1476+t1477+t1478-t1479-t1480+t148
     #1+t8810+t8811-t8812+t8813-t8814-t8816+t8818-t8820
      t9838 = -t8824-t8826+t8828+t8830-t8832+t8834+t8836-t8838+t8840+t88
     #42-t8844-t8846+t8848+t8850-t8852-t8854-t8856+t8858-t8860
      t9845 = t7897*t7580*t7177
      A1e0c6p1t1n9 = -12*(t9830+t9831+t9833+t9834+t8809+t9837+t9838)*t98
     #03*t7660*t9845
      A1e0c7p1t1n9 = 0
      t9848 = t8872+t8874+t8876+t8878+t8880-t8882+t8884-t8886+t8888-t889
     #0-t8892+t8894-t8896-t8898-t8900+t8902+t8904
      t9849 = t8907+t8909-t8911+t8913-t8915+t8917+t8918+t8920-t8922-t892
     #4+t8926+t8928-t8931-t8933-t8936-t8938+t8940
      t9851 = -t8944-t8946-t8948-t8949+t8950+t8951+t2099+t5270+t5271-t19
     #35+t8952+t8953-t1946+t2102-t8954-t8955-t2110
      A1e0c8p1t1n9 = -12*(t9848+t9849+t9851+t8960)*t9803*t2155*t9819
      t9858 = -t8972-t8974-t8977+t8979+t8982+t8984-t8985-t8986-t8989-t89
     #92-t8994+t8996-t8998-t9001-t9004+t9007-t9009-t9011-t9015+t9017
      t9863 = -t9045+t9047-t9049+t9051-t9053-t9055-t9057+t9059-t9061+t90
     #63+t9066
      t9867 = t9072+t9075+t9077+t9079-t9081-t9083+t9085+t9087+t9090-t909
     #2-t9094-t5466-t5467+t9095-t9098-t9100+t9102-t9105-t9107-t9109
      t9869 = t9141-t9143+t9145-t9148+t9150+t9153+t9156-t9157-t9158-t916
     #0-t9162
      t9873 = -t9187+t9189+t9194-t9196-t9199+t9201+t9203+t9205-t9207-t92
     #08+t3568
      t9879 = t9250+t9252-t9254-t9256+t9258+t9260+t9262+t9264+t9266+t926
     #8-t9270
      t9882 = -t9296+t9298+t9302-t9304+t9306-t9308+t9310+t9312-t9314+t93
     #16-t9319
      t9889 = -t9369+t9370+t9371+t4592-t4602+t4617-t4619-t4632+t4634-t46
     #36-t4638
      t9892 = t9383+t9385+t9388+t9390+t9392-t9395+t9397+t9399-t9401-t940
     #3+t9405
      t9895 = t9436+t9438+t9440-t9442+t9444-t9446+t9448+t9450-t9452+t945
     #4-t9456
      t9903 = t4936*t6652
      t9905 = t9466*t9810
      A1e0c9p1t1n9 = 12*(t9027-t3872+t9022-t9025+t9115+t9118-t9121-t9124
     #-t9126-t9128+t9130+t9132-t9136+t9112-t4128-t2599-t2601-t2387+t2383
     #+t2381+t2378-t2375+t3870+t2370+t9247+t4119+t3876+t3878-t3880+t3882
     #-t5391-t9242-t9245+t9382+t9357+t9359-t9361+t9363+t9365+t9351+t9353
     #-t9355-t9340+t9342-t9345-t9348-t9334+t9336-t9338+t9238-t9239+t9882
     #+t9873+t9879+t9858+t9863+t9867+t9869+t9327-t9330-t9332+t9325-t9290
     #+t9293+t9220+t9284+t9286+t9288-t9278+t9280-t9282-t9276-t9274+t9171
     #+t9172-t9174-t9176+t9178-t9180+t9182-t9184-t2364+t9167-t9169+t2255
     #+t9889+t9892+t9895+t9038-t9039+t9041+t9035-t9036-t9037+t9033+t9034
     #-t9433)*t9803*t432*t9903*t9905
      A1e0c1p1t2n9 = 0
      A1e0c2p1t2n9 = A1e0c2p1t1n9
      A1e0c3p1t2n9 = A1e0c3p1t1n9
      A1e0c4p1t2n9 = 0
      A1e0c5p1t2n9 = A1e0c5p1t1n9
      A1e0c6p1t2n9 = A1e0c6p1t1n9
      A1e0c7p1t2n9 = 0
      A1e0c8p1t2n9 = A1e0c8p1t1n9
      A1e0c9p1t2n9 = A1e0c9p1t1n9
      A2e0c1p1t1n9 = 0
      A2e0c2p1t1n9 = 48*(-t4970-t202+t4972+t4971+t9473+t9475-t9471)*sqq*
     #Mzsq*t9803*t432*t9827
      t9913 = t8683-t8688+t8702+t8704-t8706+t8707-t8709+t8710-t8717+t502
     #4-t5026-t5028
      t9914 = t9485+t9486-t9487+t9488-t8723+t8727-t9489-t8743+t9490+t876
     #2+t9491-t8772-t8776
      t9916 = t8778-t8780-t8782+t8785-t8787-t8790+t8792-t8794+t8796+t505
     #3-t803-t804-t5054
      t9921 = t5149-t1346+t9508-t1355+t1358+t1360-t5152+t5153-t9510+t951
     #2+t9514+t9516+t9518
      t9923 = -t9522-t1371+t1330+t1328+t1401+t1404+t5190-t1474-t1475-t14
     #76+t5195-t5196+t1481
      t9924 = -t9524-t8824-t8828-t8830+t9525-t9526+t8840-t8846+t8852-t88
     #54-t8856+t8858-t8860
      A2e0c3p1t1n9 = -48*(t9913+t9914+t9916-t9499-t9507+t9921+t9923+t992
     #4)*t9803*t7660*t9845
      A2e0c4p1t1n9 = 0
      t9933 = t5320+t9543+t9545-t9549-t132+t5323-t5314+t9546+t5312-t8645
     #-t9547
      A2e0c5p1t1n9 = 12*(-t5295-t511+t5316+t5318+t9538+t9540-t9536-t5297
     #-t444+t5322+t9933)*t9803*t2155*t9810
      t9938 = t9017+t9022-t2364+t2370-t2375+t2378+t2381+t5346-t5347-t534
     #9+t2394-t2398+t9559+t9560-t9561
      t9941 = -t9565+t9566-t5367-t5368+t2540-t2544+t5370+t2552-t9567+t53
     #72+t9568-t2599+t9569+t9041-t9045
      t9944 = t9184+t9187-t5617+t9189+t9194+t9196-t9201+t9205+t9207+t563
     #4-t3583+t3585+t3587+t5636
      t9945 = -t3597-t3608+t3609+t3610-t5672-t5717+t3800-t5718+t3804+t38
     #15-t3816-t3817+t3819+t9583-t9585
      t9949 = -t9591-t9593+t9595-t9597-t9599-t3827-t9600+t3835-t3837-t57
     #59+t5760+t3844-t5761+t5762
      t9950 = -t3852+t5763-t5764-t3868+t3870-t5765+t5766+t3878-t3880-t92
     #42-t9245+t9258-t9602-t9603+t9604
      t9952 = t9288-t9290+t9298+t9302-t9607+t9608+t9610+t4119-t4128-t961
     #1-t9612-t9613-t9614+t9615
      t9953 = -t9617+t9618+t5837-t5839+t9619+t9620+t5845-t9621-t9622-t96
     #23+t9624+t5857-t5860+t9334+t9340
      t9956 = t4258-t9630+t9632+t9633-t4385-t4388+t4390+t5883-t4393+t439
     #5-t4397-t4399-t4400-t9634
      t9957 = t9636-t9353+t9355-t9363+t9637-t9638+t5901-t5904+t5906+t590
     #8-t5910-t9639-t9640+t9641-t9642
      t9959 = t5921+t5923-t5939-t5941-t9645+t5944+t5945+t9395+t9397+t941
     #1-t9413+t9415+t9417+t9420-t9424
      A2e0c6p1t1n9 = 12*(t9558+t9938-t9564+t9941+t9581+t9944+t9945+t9949
     #+t9950+t9952+t9953+t9956+t9957+t9959+t9659)*t9803*t432*t9903*t9905
      A2e0c1p1t2n9 = 0
      A2e0c2p1t2n9 = A2e0c2p1t1n9
      A2e0c3p1t2n9 = A2e0c3p1t1n9
      A2e0c4p1t2n9 = 0
      A2e0c5p1t2n9 = A2e0c5p1t1n9
      A2e0c6p1t2n9 = A2e0c6p1t1n9

      RETURN
      END
