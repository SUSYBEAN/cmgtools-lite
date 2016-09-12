
      SUBROUTINE CrsecSeqz_msv()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      t1 = Ml**2
      t2 = 4*t1
      t3 = BCsTh**2
      t4 = s*t3
      t5 = s+t2+t4
      A1e0c1p1t1z = 12*t5/s
      t8 = Mzsq*s
      t9 = s**2
      t11 = abs(DCMPLX(Mzsq))
      t12 = t11**2
      t14 = Mzsq**2
      t18 = 1/(-t9*Mzsq+s*t12+s*t14-t12*Mzsq)
      A1e0c2p1t1z = -12*t8*t5*t18
      A1e0c3p1t1z = -12*t8*(-t2+s+t4)*t18
      t27 = 1/(-s+Mzsq)
      A1e0c4p1t1z = -24*t5*t27
      A1e0c1p1t2z = A1e0c1p1t1z
      A1e0c2p1t2z = A1e0c2p1t1z
      A1e0c3p1t2z = A1e0c3p1t1z
      A1e0c4p1t2z = A1e0c4p1t1z
      A2e0c1p1t1z = -96*t18*Mzsq*t9*BCsTh
      A2e0c2p1t1z = -48*s*t27*BCsTh
      A2e0c1p1t2z = A2e0c1p1t1z
      A2e0c2p1t2z = A2e0c2p1t1z

      RETURN
      END

      SUBROUTINE CrsecSeqs_msv()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      LOGICAL CIsBad
      LOGICAL IsBad
      DOUBLE COMPLEX spf10
      DOUBLE COMPLEX spf11
      DOUBLE COMPLEX spf12
      DOUBLE COMPLEX spf13
      DOUBLE COMPLEX spf4
      DOUBLE COMPLEX spf5
      DOUBLE COMPLEX spf6
      DOUBLE COMPLEX spf7
      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX spf9
      DOUBLE COMPLEX DENSPENCE
      DOUBLE PRECISION epsMl
      PARAMETER (epsMl=1d-12)
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t1 = 1/s
      t2 = log(DCMPLX(t1))
      t3 = log(DCMPLX(dels))
      t4 = t2*t3
      t5 = Ml**2
      t7 = 48*t4*t5
      t8 = BCsTh**2
      t9 = s*t8
      t11 = 12*t4*t9
      t12 = t3*s
      t13 = 24*t12
      t14 = t2*s
      t15 = 12*t14
      t16 = PI**2
      t17 = t16*s
      t18 = t3**2
      t19 = t18*s
      t20 = 12*t19
      t21 = t2**2
      t22 = t21*s
      t23 = 3*t22
      t24 = t18*t5
      t25 = 48*t24
      t26 = t21*t5
      t27 = 12*t26
      t28 = t16*t5
      t29 = 4*t28
      t30 = t17*t8
      t32 = 12*t19*t8
      t34 = 3*t22*t8
      t35 = t4*s
      t36 = 12*t35
      t37 = t7+t11-t13+t15+t17-t20-t23-t25-t27+t29+t30-t32-t34+t36
      t39 = 1/t16
      A1e0c1p1t1s = -t1*t37*t39/2
      t42 = 1+BtQ
      t43 = 1/t42
      spf4 = DENSPENCE(DCMPLX(2*BtQ*t43),spcZip)
      t46 = 1/t5
      t49 = t42*t46*s/2
      t50 = -1+BtQ
      t51 = 1/t50
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf6 = t49
      else
      spf6 = -t42*t51
      end if
      if ((ABS(DBLE(spf6)).eq.0d0).or.(CIsBad(spf6))) spf6 = t49
      spf6 = log(spf6)
      t53 = s**2
      t55 = t3*BtQ
      t56 = t5*s
      t58 = 16*t55*t56
      t59 = t53*t8
      t61 = 4*t55*t59
      t62 = t2*BtQ
      t64 = 8*t62*t56
      t66 = 2*t62*t59
      t67 = spf4*t5
      t69 = 8*t67*t9
      t70 = spf6**2
      t71 = t70*t5
      t73 = 2*t71*t9
      t74 = t3*spf6
      t75 = t74*t56
      t78 = 4*t74*t59
      t79 = spf6*t2
      t80 = t79*t56
      t83 = 2*t79*t59
      t84 = t55*t53
      t85 = 4*t84
      t86 = t62*t53
      t87 = 2*t86
      t88 = spf4*t53
      t90 = 4*t88*t8
      t91 = t67*s
      t93 = t58+t61-t64-t66-t69-t73-8*t75-t78+4*t80+t83+t85-t87+t90+8*t9
     #1
      t94 = spf6*t53
      t96 = 2*t94*t8
      t97 = t70*t53
      t98 = t97*t8
      t99 = t71*s
      t101 = t74*t53
      t102 = 4*t101
      t103 = t5**2
      t105 = 32*t74*t103
      t106 = t79*t53
      t107 = 2*t106
      t109 = 16*t79*t103
      t110 = t56*t8
      t112 = 8*t74*t110
      t114 = 4*t79*t110
      t115 = t53*BtQ
      t116 = 4*t115
      t117 = spf4*t103
      t118 = 32*t117
      t119 = 4*t88
      t120 = t70*t103
      t121 = 8*t120
      t122 = 6*t94
      t123 = -t96+t98+2*t99-t102+t105+t107-t109+t112-t114+t116-t118+t119
     #-t121-t122+t97
      t124 = t93+t123
      t126 = 1/BtQ
      t127 = t39*t126
      A1e0c2p1t1s = -3.D0/2.D0/t53*t124*t127
      t130 = CsTh*s
      t132 = 2*CsTh*t5
      t133 = s-t130+t132
      t134 = 1/t133
      t136 = 2*t5*t134
      t137 = -1+BCsTh
      t138 = 1/t137
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf5 = t136
      else
      spf5 = t50*t138
      end if
      if ((ABS(DBLE(spf5)).eq.0d0).or.(CIsBad(spf5))) spf5 = t136
      spf5 = log(spf5)
      t140 = t53**2
      t141 = t140*s
      t142 = CsTh*t141
      t143 = t5*t140
      t144 = 2*t143
      t146 = 2*t143*CsTh
      t147 = t53*s
      t148 = t103*t147
      t149 = 2*t148
      t151 = 2*t148*CsTh
      t153 = t103*t5*t53
      t154 = 4*t153
      t156 = 4*t153*CsTh
      t157 = t103**2
      t158 = t157*s
      t159 = 10*t158
      t161 = 10*t158*CsTh
      t162 = t157*t5
      t163 = 28*t162
      t165 = 28*t162*CsTh
      t166 = -t141-t142+t144+t146+t149+t151+t154+t156+t159+t161+t163+t16
     #5
      t167 = 1/t140
      t168 = t166*t167
      t170 = t168*t46/2
      t171 = BtQ+BCsTh
      if (ABS(DBLE(Ml**2/s)**5).lt.epsMl) then
      spf7 = t170
      else
      spf7 = t171*t51
      end if
      if (CIsBad(spf7)) spf7 = t170
      spf7 = DENSPENCE(spf7,spcZip)
      t173 = -s-t130+t132
      t174 = 1/t173
      t176 = 2*t5*t174
      t177 = 1+BCsTh
      t178 = 1/t177
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf8 = -t176
      else
      spf8 = -t50*t178
      end if
      if ((ABS(DBLE(spf8)).eq.0d0).or.(CIsBad(spf8))) spf8 = -t176
      spf8 = log(spf8)
      t180 = t141-t142-t144+t146-t149+t151-t154+t156-t159+t161-t163+t165
      t181 = t180*t167
      t182 = t181*t174
      t183 = -BtQ+BCsTh
      if (ABS(DBLE(Ml**2/s)**5).lt.epsMl) then
      spf9 = t182
      else
      spf9 = t183*t178
      end if
      if (CIsBad(spf9)) spf9 = t182
      spf9 = DENSPENCE(spf9,spcZip)
      t185 = t168*t134
      if (ABS(DBLE(Ml**2/s)**5).lt.epsMl) then
      spf10 = t185
      else
      spf10 = t171*t138
      end if
      if (CIsBad(spf10)) spf10 = t185
      spf10 = DENSPENCE(spf10,spcZip)
      t188 = t181*t46/2
      if (ABS(DBLE(Ml**2/s)**5).lt.epsMl) then
      spf11 = -t188
      else
      spf11 = -t183*t51
      end if
      if (CIsBad(spf11)) spf11 = -t188
      spf11 = DENSPENCE(spf11,spcZip)
      t190 = t133**2
      t192 = t46*t43
      t194 = t190*t1*t192/2
      t195 = t137**2
      t196 = BtQ**2
      t198 = 1/(-1+t196)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf12 = t194
      else
      spf12 = -t195*t198
      end if
      if ((ABS(DBLE(spf12)).eq.0d0).or.(CIsBad(spf12))) spf12 = t194
      spf12 = log(spf12)
      t200 = t173**2
      t203 = t200*t1*t192/2
      t204 = t177**2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf13 = t203
      else
      spf13 = -t204*t198
      end if
      if ((ABS(DBLE(spf13)).eq.0d0).or.(CIsBad(spf13))) spf13 = t203
      spf13 = log(spf13)
      t206 = t3*spf12
      t208 = 2*t206*t9
      t209 = spf12*t2
      t210 = t209*t9
      t211 = t3*spf13
      t213 = 2*t211*t9
      t214 = spf13*t2
      t215 = t214*t9
      t217 = 4*t209*t5
      t218 = spf9*s
      t220 = 2*t218*t8
      t221 = spf7*s
      t223 = 2*t221*t8
      t224 = t206*s
      t225 = 2*t224
      t227 = 8*t206*t5
      t228 = t209*s
      t229 = t211*s
      t230 = 2*t229
      t232 = 8*t211*t5
      t233 = t214*s
      t235 = 4*t214*t5
      t236 = spf8**2
      t237 = t236*s
      t238 = t237*t8
      t239 = spf5**2
      t240 = t239*s
      t241 = t240*t8
      t242 = -t208+t210+t213-t215+t217-t220+t223-t225-t227+t228+t230+t23
     #2-t233-t235+t238-t241
      t243 = spf10*s
      t245 = 2*t243*t8
      t246 = spf11*s
      t248 = 2*t246*t8
      t249 = t236*t5
      t250 = 4*t249
      t251 = spf11*t5
      t252 = 8*t251
      t253 = 2*t246
      t254 = spf13*s
      t255 = 2*t254
      t256 = spf10*t5
      t257 = 8*t256
      t258 = 2*t243
      t259 = spf12*s
      t260 = 2*t259
      t261 = t239*t5
      t262 = 4*t261
      t263 = spf9*t5
      t264 = 8*t263
      t265 = 2*t218
      t266 = spf7*t5
      t267 = 8*t266
      t268 = 2*t221
      t269 = t245-t248+t250+t237-t252-t253+t255+t257+t258-t260-t262-t240
     #-t264-t265+t267+t268
      t271 = (t242+t269)*t39
      A1e0c3p1t1s = -3*t271*t1
      t274 = Mzsq*s
      t276 = t53*Mzsq
      t277 = abs(DCMPLX(Mzsq))
      t278 = t277**2
      t280 = Mzsq**2
      t284 = 1/(-t276+s*t278+s*t280-t278*Mzsq)
      A1e0c4p1t1s = t274*t37*t39*t284/2
      t288 = t127*t284
      A1e0c5p1t1s = 3.D0/2.D0*Mzsq*t124*t288
      A1e0c6p1t1s = 3*t274*t271*t284
      t296 = t25+t27-t29+t30-t32-t34+t36-t7+t11-t13+t15+t17-t20-t23+96*t
     #3*t5-48*t2*t5
      A1e0c7p1t1s = t274*t296*t39*t284/2
      t305 = spf6*t5
      t306 = t305*s
      t308 = -t58+t61+t64-t66-t69-t73+24*t75-t78-12*t80+t83-16*t56*BtQ+t
     #85-t87+t90-24*t91+32*t306
      t312 = -t96+t98-6*t99-t102-t105+t107+t109+t112-t114+t116+t118+t119
     #+t121-t122+t97-32*spf6*t103
      A1e0c8p1t1s = 3.D0/2.D0*Mzsq*(t308+t312)*t288
      t316 = -t208+t210+t213-t215-t217-t220+t223-t225+t227+t228+t230-t23
     #2-t233+t235+t238-t241+t245
      t321 = -t248-t250+t237+t252-t253+t255-t257+t258-t260+t262-t240+t26
     #4-t265-t267+t268-8*spf13*t5+8*spf12*t5
      A1e0c9p1t1s = 3*t274*(t316+t321)*t39*t284
      A1e0c10p1t1s = A1e0c8p1t1s
      t326 = t16*Mzsq
      t328 = conjg(Mzsq)
      t329 = t16*t328
      t331 = t18*t328
      t334 = t18*Mzsq
      t337 = t21*Mzsq
      t340 = t21*t328
      t347 = t5*Mzsq
      t350 = t5*t328
      t356 = s*t328
      t363 = t3*Mzsq
      t366 = t21*t53
      t369 = t16*t53
      t378 = 12*t4*t274+12*t4*t356-8*t28*s-3*t340*s-24*t363*s+6*t366*t8-
     #2*t369*t8+4*t28*Mzsq+t326*s+t329*s+4*t28*t328
      t382 = t18*t53
      t399 = t2*Mzsq
      t411 = t274*t8
      t414 = t356*t8
      t424 = 12*t14*t328-3*t337*s-12*t26*t328-24*t4*t53+12*t4*t411+12*t4
     #*t414+48*t3*t53+24*t382-24*t2*t53+6*t366-2*t369
      t430 = 1/(t53-t356-t274+Mzsq*t328)
      A1e0c11p1t1s = (t326*t9+t329*t9-12*t331*t9-12*t334*t9-3*t337*t9-3*
     #t340*t9-96*t4*t56-24*t4*t59+48*t4*t347+48*t4*t350+t378+96*t24*s+24
     #*t382*t8-48*t24*Mzsq-24*t12*t328-12*t331*s-12*t334*s-48*t24*t328+2
     #4*t26*s-12*t26*Mzsq+12*t399*s+t424)*t39*t430/2
      t432 = t5*t53
      t437 = t147*t8
      t440 = t328*t53
      t458 = t328*t8
      t461 = Mzsq*t8
      t468 = spf6*t328
      t471 = spf6*Mzsq
      t478 = t103*s
      t481 = 8*t67*t356+4*t88*t458+4*t88*t461+4*t71*t59+2*t71*t274-2*t46
     #8*t59-2*t471*t59+2*t71*t356+t97*t461+t97*t458-64*t74*t478
      t487 = t103*Mzsq
      t490 = t103*t328
      t507 = 16*t74*t432+8*t74*t437+32*t74*t487+32*t74*t490-4*t74*t276-4
     #*t71*t53+4*t276*BtQ-16*t67*t53+16*t120*s+64*t117*s-8*t55*t147
      t512 = spf4*t147
      t525 = t70*t147
      t532 = 4*t62*t147-32*t117*Mzsq-8*t512*t8-8*t120*Mzsq+4*t115*t328-3
     #2*t117*t328+4*t88*t328+4*t88*Mzsq-2*t525*t8-6*t468*t53-6*t471*t53
      t535 = spf6*t147
      t546 = t56*Mzsq
      t549 = t56*t328
      t552 = t276*t8
      t555 = t440*t8
      t560 = 4*t535*t8-8*t120*t328+t97*Mzsq+t97*t328+8*t74*t147-4*t79*t1
     #47+16*t55*t546+16*t55*t549+4*t55*t552+4*t55*t555-8*t62*t546
      t576 = t432*t8
      t585 = -8*t62*t549-2*t62*t555-2*t62*t552-8*t67*t414-8*t67*t411-2*t
     #70*t8*t546-2*t71*t414-16*t74*t576-8*t74*t549-8*t74*t546-4*t74*t552
      t609 = -4*t74*t555+8*t79*t576+4*t79*t546+4*t79*t549+2*t79*t552+2*t
     #79*t555-4*t74*t440+32*t79*t478-8*t79*t432-4*t79*t437-16*t79*t487
      t616 = t74*t5
      t621 = t79*t5
      t631 = -16*t79*t490+2*t79*t440+2*t79*t276+8*t616*t411+8*t616*t414-
     #4*t621*t411-4*t621*t414-8*t147*BtQ-8*t512+12*t535-2*t525
      t636 = t127*t430
      A1e0c12p1t1s = 3.D0/2.D0*t1*(-32*t55*t432+4*t55*t276-8*t55*t437+4*
     #t55*t440+16*t62*t432-2*t62*t440-2*t62*t276+4*t62*t437+16*t67*t59+8
     #*t67*t274+t481+t507+t532+t560+t585+t609+t631)*t636
      t658 = spf9*t53
      t661 = -2*t206*t411-2*t206*t414+t209*t411+t209*t414+8*t261*s+2*t25
     #4*Mzsq+4*t249*Mzsq+16*t263*s-8*t263*Mzsq-8*t263*t328-2*t218*Mzsq+4
     #*t658*t8
      t662 = spf9*t328
      t678 = t236*t328
      t680 = t236*t53
      t687 = -2*t662*s-2*t259*Mzsq-2*t259*t328+4*t206*t53-2*t209*t53-8*t
     #249*s+4*t249*t328+t237*Mzsq+t678*s-2*t680*t8-16*t256*s+8*t256*Mzsq
      t693 = spf10*t328
      t696 = spf10*t53
      t710 = t239*t53
      t713 = t239*t328
      t715 = 8*t256*t328+2*t243*Mzsq+2*t693*s-4*t696*t8+2*t254*t328-4*t2
     #11*t53+2*t214*t53-4*t261*Mzsq-4*t261*t328-t240*Mzsq+2*t710*t8-t713
     #*s
      t724 = spf11*t328
      t727 = spf11*t53
      t738 = spf7*t53
      t741 = spf7*t328
      t744 = 16*t251*s-8*t251*Mzsq-8*t251*t328-2*t246*Mzsq-2*t724*s+4*t7
     #27*t8-16*t266*s+8*t266*Mzsq+8*t266*t328+2*t221*Mzsq-4*t738*t8+2*t7
     #41*s
      t767 = 2*t211*t411+2*t211*t414-t214*t411-t214*t414+t678*t9+t237*t4
     #61+2*t243*t461+2*t693*t9-2*t662*t9-2*t218*t461+16*t206*t56-8*t206*
     #t347
      t790 = -8*t206*t350-2*t206*t274+4*t206*t59-2*t206*t356-8*t209*t56+
     #4*t209*t347+4*t209*t350+t209*t274-2*t209*t59+t209*t356-16*t211*t56
     #+8*t211*t347
      t814 = 8*t211*t350+2*t211*t274-4*t211*t59+2*t211*t356+8*t214*t56-4
     #*t214*t347-4*t214*t350-t214*t274+2*t214*t59-t214*t356-2*t246*t461-
     #2*t724*t9
      t831 = -t240*t461-t713*t9+2*t221*t461+2*t741*t9-2*t680+4*t727-4*sp
     #f13*t53-4*t696+4*spf12*t53+2*t710+4*t658-4*t738
      A1e0c13p1t1s = 3*(t661+t687+t715+t744+t767+t790+t814+t831)*t39*t43
     #0
      A1e0c1p1t2s = A1e0c1p1t1s
      A1e0c2p1t2s = A1e0c2p1t1s
      A1e0c3p1t2s = A1e0c3p1t1s
      A1e0c4p1t2s = A1e0c4p1t1s
      A1e0c5p1t2s = A1e0c5p1t1s
      A1e0c6p1t2s = A1e0c6p1t1s
      A1e0c7p1t2s = A1e0c7p1t1s
      A1e0c8p1t2s = A1e0c10p1t1s
      A1e0c9p1t2s = A1e0c9p1t1s
      A1e0c10p1t2s = A1e0c8p1t2s
      A1e0c11p1t2s = A1e0c11p1t1s
      A1e0c12p1t2s = A1e0c12p1t1s
      A1e0c13p1t2s = A1e0c13p1t1s
      t837 = t276*BCsTh
      A2e0c1p1t1s = 4*t837*(12*t4-12+t16-12*t18-3*t21+18*t2-36*t3)*t39*t
     #284
      t854 = s*BtQ
      t858 = t70*s
      t859 = spf4*s
      t868 = -8*spf6*s+4*t55*s-2*t62*s+6*t854+8*t616-4*t621+t858+4*t859-
     #2*t71-8*t67-4*t74*s+2*t79*s+12*t305
      A2e0c2p1t1s = 12*t274*BCsTh*t868*t39*t126*t284
      t881 = 3*spf13+2*spf7-2*spf9+2*spf10-2*spf11-3*spf12+t236-t239-t21
     #4+t209+2*t211-2*t206
      A2e0c3p1t1s = 24*t837*t881*t39*t284
      t885 = s*BCsTh
      t910 = t326+t329-12*t331-3*t337-12*t334-3*t340-12*t328-36*t363+18*
     #t399+18*t2*t328-36*t3*t328
      A2e0c4p1t1s = t885*(-12*Mzsq-24*t35+24*s+6*t22+24*t19-2*t17+72*t12
     #-36*t14+12*t4*Mzsq+12*t4*t328+t910)*t39*t430
      t945 = -16*t75+8*t80+6*t274*BtQ+6*t854*t328+4*t859*t328-8*t67*t328
     #-8*t67*Mzsq+4*t859*Mzsq-2*t71*Mzsq-8*t468*s-8*t471*s+t858*Mzsq+12*
     #t305*Mzsq+12*t305*t328+t858*t328-2*t71*t328-8*t84+4*t86+16*t91
      t978 = -24*t306+4*t99+8*t101-4*t106+4*t55*t274+4*t55*t356-2*t62*t3
     #56-2*t62*t274+8*t74*t347-4*t74*t274-4*t74*t356+8*t74*t350-4*t79*t3
     #47+2*t79*t274+2*t79*t356-4*t79*t350-12*t115-8*t88+16*t94-2*t97
      A2e0c5p1t1s = 3*BCsTh*(t945+t978)*t636
      t1004 = -2*t206*Mzsq-2*t206*t328+t209*Mzsq+t209*t328+2*t211*Mzsq+2
     #*t211*t328-t214*Mzsq-t214*t328+4*t224-2*t228-4*t229+2*t233-2*t237+
     #4*t246-6*t254-4*t243+6*t259+2*t240
      t1029 = 4*t218-4*t221-2*spf11*Mzsq-2*t724+3*spf13*Mzsq+3*spf13*t32
     #8+2*spf10*Mzsq+2*t693-3*spf12*Mzsq-3*spf12*t328+t236*Mzsq+t678-t23
     #9*Mzsq-t713-2*spf9*Mzsq-2*t662+2*spf7*Mzsq+2*t741
      A2e0c6p1t1s = 6*t885*(t1004+t1029)*t39*t430
      A2e0c1p1t2s = A2e0c1p1t1s
      A2e0c2p1t2s = A2e0c2p1t1s
      A2e0c3p1t2s = A2e0c3p1t1s
      A2e0c4p1t2s = A2e0c4p1t1s
      A2e0c5p1t2s = A2e0c5p1t1s
      A2e0c6p1t2s = A2e0c6p1t1s

      RETURN
      END

      SUBROUTINE CrsecSeqc_msv()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX spf9
      DOUBLE COMPLEX DENSPENCE
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t1 = 1-x1
      spf8 = DENSPENCE(DCMPLX(t1),spcZip)
      t2 = 1-x2
      spf9 = DENSPENCE(DCMPLX(t2),spcZip)
      t3 = Ml**2
      t4 = PI**2
      t5 = t4/3
      t6 = x2**2
      t8 = log(DCMPLX(t2))
      t13 = log(DCMPLX(x2))
      t18 = t8**2
      if (DISonQED .eq. 1) then
      t20 = t5+t6+7.D0/2.D0*x2+3*t8-t8*x2-t8*t6/2+t13*x2+t13*t6/2-2*spf9
     #-t18
      else
      t20 = 0
      end if
      t21 = t3*t20
      t22 = cz1*cz2
      t24 = 4*t21*t22
      t25 = s*t20
      t26 = t25*t22
      t27 = cz1**2
      t28 = cz2**2
      t29 = t27*t28
      t31 = 4*t21*t29
      t32 = t25*t29
      t33 = t27*cz2
      t35 = 4*t21*t33
      t36 = t25*t33
      t37 = x1**2
      t39 = log(DCMPLX(t1))
      t44 = log(DCMPLX(x1))
      t49 = t39**2
      if (DISonQED .eq. 1) then
      t51 = t5+t37+7.D0/2.D0*x1+3*t39-t39*x1-t39*t37/2+t44*x1+t44*t37/2-
     #2*spf8-t49
      else
      t51 = 0
      end if
      t52 = t3*t51
      t54 = 4*t52*t33
      t55 = s*t51
      t56 = t55*t33
      t58 = 4*t52*t22
      t59 = t55*t22
      t61 = 4*t52*t29
      t62 = t24+t26+t31+t32-t35-t36-t54-t56+t58+t59+t61
      t63 = t55*t29
      t64 = cz1*t28
      t66 = 4*t52*t64
      t67 = t55*t64
      t69 = 4*t21*t64
      t70 = t25*t64
      t71 = 1/s
      t72 = log(DCMPLX(t71))
      t73 = t72*s
      t75 = 3*t73*t33
      t76 = t72*t3
      t78 = 12*t76*t33
      t80 = 3*t73*t22
      t82 = 12*t76*t22
      t84 = 3*t73*t29
      t86 = 12*t76*t29
      t87 = log(DCMPLX(dels))
      t88 = t87*s
      t90 = 8*t88*t33
      t91 = t63-t66-t67-t69-t70-t75-t78+t80+t82+t84+t86+t90
      t94 = 8*t88*t22
      t96 = 8*t88*t29
      t98 = 8*t88*t64
      t99 = t72*t87
      t100 = t99*s
      t101 = BCsTh**2
      t102 = t101*t27
      t105 = 4*t100*t102*cz2
      t106 = t101*cz1
      t109 = 4*t100*t106*cz2
      t110 = delc*s
      t111 = muf**2
      t112 = 1/t111
      t115 = log(DCMPLX(t110*t112*dels))
      t116 = t115**2
      t117 = s*t116
      t119 = 2*t117*t22
      t120 = t3*t116
      t122 = 8*t120*t29
      t124 = 2*t117*t29
      t126 = 8*t120*t64
      t128 = 2*t117*t64
      t130 = 8*t120*t33
      t132 = 2*t117*t33
      t133 = -t94-t96+t98-t105+t109-t119-t122-t124+t126+t128+t130+t132
      t135 = 8*t120*t22
      t137 = 3*t73*t64
      t139 = 12*t76*t64
      t141 = log(DCMPLX(s*t112))
      t142 = t3*t141
      t144 = 12*t142*t64
      t145 = s*t141
      t147 = 3*t145*t64
      t149 = 12*t142*t33
      t151 = 3*t145*t33
      t153 = 12*t142*t22
      t155 = 3*t145*t22
      t157 = 12*t142*t29
      t159 = 3*t145*t29
      t160 = s*j4
      t163 = 1/(-1+cz2)
      t166 = log(DCMPLX(-t110*cz2*t112*t163))
      t167 = t166*cz1
      t169 = 2*t160*t167
      t170 = -t135-t137-t139-t144-t147-t149-t151+t153+t155+t157+t159-t16
     #9
      t173 = t3*j4
      t174 = t166*t27
      t176 = 8*t173*t174
      t178 = 2*t160*t174
      t181 = 1/(-1+cz1)
      t184 = log(DCMPLX(-t110*cz1*t112*t181))
      t185 = t184*cz2
      t187 = 8*t173*t185
      t189 = 2*t160*t185
      t190 = t184*t28
      t192 = 8*t173*t190
      t194 = 2*t160*t190
      t195 = -t1
      t200 = log(DCMPLX(-t110*t195*t112/x1))
      t201 = t200**2
      t202 = t3*t201
      t204 = 4*t202*t33
      t205 = s*t201
      t206 = t205*t33
      t208 = 4*t202*t22
      t209 = t205*t22
      t211 = 4*t202*t29
      t212 = t176+t178-t187-t189+t192+t194-t204-t206+t208+t209+t211
      t213 = t205*t29
      t215 = 4*t202*t64
      t216 = t205*t64
      t217 = -t2
      t222 = log(DCMPLX(-t110*t217*t112/x2))
      t223 = t222**2
      t224 = t3*t223
      t226 = 4*t224*t64
      t227 = s*t223
      t228 = t227*t64
      t230 = 4*t224*t22
      t231 = t227*t22
      t233 = 4*t224*t29
      t234 = t227*t29
      t236 = 4*t224*t33
      t237 = t227*t33
      t239 = 8*t173*t167
      t240 = t213-t215-t216-t226-t228+t230+t231+t233+t234-t236-t237-t239
      t242 = s*cz1
      t243 = t242*t28
      t244 = 6*t243
      t245 = s*t27
      t246 = t245*t28
      t247 = 6*t246
      t248 = s*t101
      t249 = t248*j4
      t251 = log(DCMPLX(-cz1*t181))
      t260 = 1/cz1
      if (DISonQED .eq. 1) then
      t263 = (4*t251-3-4*cz1*t251+10*cz1+2*t27*t251-4*t27)*t260/2
      else
      t263 = 0
      end if
      t264 = t263*cz1
      t265 = t264*t28
      t266 = t249*t265
      t267 = t263*t27
      t268 = t267*t28
      t269 = t249*t268
      t270 = t173*x1
      t272 = 4*t270*t265
      t274 = 4*t270*t268
      t275 = t160*x1
      t276 = t275*t265
      t277 = t275*t268
      t278 = t3*cz1
      t279 = t278*t28
      t281 = 16*t99*t279
      t283 = 4*t99*t243
      t284 = t278*cz2
      t286 = 16*t99*t284
      t287 = t242*cz2
      t289 = 4*t99*t287
      t290 = t244-t247-t266+t269+t272-t274+t276-t277-t281-t283+t286+t289
      t291 = t245*cz2
      t293 = 4*t99*t291
      t294 = t3*t27
      t295 = t294*cz2
      t297 = 16*t99*t295
      t298 = x2*t166
      t299 = t298*cz1
      t301 = 8*t173*t299
      t303 = 2*t160*t299
      t304 = t298*t27
      t306 = 8*t173*t304
      t308 = 2*t160*t304
      t309 = j4*t166
      t310 = t309*t27
      t312 = 2*t248*t310
      t313 = t309*cz1
      t315 = 2*t248*t313
      t316 = 6*t291
      t318 = log(DCMPLX(-cz2*t163))
      t327 = 1/cz2
      if (DISonQED .eq. 1) then
      t330 = (4*t318-3-4*cz2*t318+10*cz2+2*t28*t318-4*t28)*t327/2
      else
      t330 = 0
      end if
      t331 = t330*cz2
      t332 = t331*cz1
      t334 = 4*t173*t332
      t335 = t330*t28
      t336 = t335*cz1
      t338 = 4*t173*t336
      t339 = t160*t332
      t340 = -t293-t297+t301+t303-t306-t308+t312-t315+t316+t334-t338+t33
     #9
      t344 = t160*t336
      t345 = t331*t27
      t347 = 4*t173*t345
      t348 = t335*t27
      t350 = 4*t173*t348
      t351 = t160*t345
      t352 = t160*t348
      t353 = t223*t28
      t354 = t353*cz1
      t355 = t248*t354
      t356 = t223*cz2
      t357 = t356*cz1
      t358 = t248*t357
      t359 = t353*t27
      t360 = t248*t359
      t361 = t356*t27
      t362 = t248*t361
      t363 = t264*cz2
      t364 = t160*t363
      t365 = t267*cz2
      t366 = t160*t365
      t367 = -t344-t347+t350-t351+t352-t355+t358+t360-t362+t364-t366
      t369 = 4*t173*t265
      t371 = 4*t173*t268
      t372 = t160*x2
      t373 = t372*t345
      t374 = t372*t348
      t375 = x1*t184
      t376 = t375*cz2
      t378 = 2*t249*t376
      t379 = t375*t28
      t381 = 2*t249*t379
      t383 = 2*t249*t299
      t385 = 2*t249*t304
      t387 = 4*t270*t363
      t389 = 4*t270*t365
      t390 = t275*t363
      t391 = t275*t365
      t392 = -t369+t371+t373-t374+t378-t381+t383-t385-t387+t389-t390+t39
     #1
      t396 = 4*t100*t102*t28
      t399 = 4*t100*t106*t28
      t400 = t249*t363
      t401 = t249*t365
      t402 = t248*t141
      t403 = t87*t27
      t404 = t403*cz2
      t406 = 4*t402*t404
      t407 = t87*cz1
      t408 = t407*cz2
      t410 = 4*t402*t408
      t411 = t403*t28
      t413 = 4*t402*t411
      t414 = t407*t28
      t416 = 4*t402*t414
      t417 = t249*t332
      t418 = t249*t336
      t419 = t173*x2
      t421 = 4*t419*t332
      t423 = 4*t419*t336
      t424 = t396-t399+t400-t401-t406+t410+t413-t416+t417-t418-t421+t423
      t425 = t372*t332
      t426 = t372*t336
      t427 = t249*t345
      t428 = t249*t348
      t430 = 4*t419*t345
      t432 = 4*t419*t348
      t433 = 6*t287
      t434 = x1*t263
      t436 = t249*t434*t64
      t438 = t249*t434*t22
      t440 = t249*t434*t33
      t442 = t249*t434*t29
      t443 = x2*t330
      t445 = t249*t443*t22
      t446 = -t425+t426-t427+t428+t430-t432-t433+t436-t438+t440-t442-t44
     #5
      t450 = t249*t443*t64
      t452 = t249*t443*t33
      t454 = t249*t443*t29
      t455 = t201*cz1
      t456 = t455*t28
      t457 = t248*t456
      t458 = t201*t27
      t459 = t458*cz2
      t460 = t248*t459
      t461 = t455*cz2
      t462 = t248*t461
      t463 = t458*t28
      t464 = t248*t463
      t466 = 16*t142*t411
      t467 = t160*t265
      t468 = t160*t268
      t469 = t20*t28
      t470 = t469*t27
      t471 = t248*t470
      t472 = t450+t452-t454-t457-t460+t462+t464+t466-t467+t468+t471
      t473 = t20*cz2
      t474 = t473*t27
      t475 = t248*t474
      t476 = t469*cz1
      t477 = t248*t476
      t478 = t473*cz1
      t479 = t248*t478
      t480 = t141*t27
      t481 = t480*cz2
      t483 = 3*t248*t481
      t484 = t141*cz1
      t485 = t484*cz2
      t487 = 3*t248*t485
      t488 = t480*t28
      t490 = 3*t248*t488
      t492 = 8*t173*t376
      t494 = 2*t160*t376
      t495 = j4*t184
      t496 = t495*t28
      t498 = 2*t248*t496
      t500 = 8*t173*t379
      t502 = 2*t160*t379
      t503 = t72*cz1
      t504 = t503*t28
      t506 = 3*t248*t504
      t507 = -t475-t477+t479-t483+t487+t490+t492+t494+t498-t500-t502-t50
     #6
      t509 = t72*t27
      t510 = t509*cz2
      t512 = 3*t248*t510
      t513 = t503*cz2
      t515 = 3*t248*t513
      t516 = t509*t28
      t518 = 3*t248*t516
      t520 = 4*t173*t363
      t522 = 4*t173*t365
      t524 = 4*t145*t411
      t526 = 16*t142*t404
      t528 = 4*t145*t404
      t530 = 16*t142*t408
      t532 = 4*t145*t408
      t534 = 16*t142*t414
      t536 = 4*t145*t414
      t537 = -t512+t515+t518+t520-t522+t524-t526-t528+t530+t532-t534-t53
     #6
      t538 = t484*t28
      t540 = 3*t248*t538
      t541 = t294*t28
      t543 = 16*t99*t541
      t545 = 4*t99*t246
      t546 = t51*cz1
      t547 = t546*t28
      t548 = t248*t547
      t549 = t51*t27
      t550 = t549*cz2
      t551 = t248*t550
      t552 = t546*cz2
      t553 = t248*t552
      t554 = t549*t28
      t555 = t248*t554
      t556 = t495*cz2
      t558 = 2*t248*t556
      t559 = t116*cz1
      t560 = t559*t28
      t562 = 2*t248*t560
      t563 = t116*t27
      t564 = t563*cz2
      t566 = 2*t248*t564
      t567 = t559*cz2
      t569 = 2*t248*t567
      t570 = t563*t28
      t572 = 2*t248*t570
      t573 = -t540+t543+t545-t548-t551+t553+t555-t558+t562+t566-t569-t57
     #2
      t577 = t62+t91+t133+t170+t212+t240+t290+t340+t367+t392+t424+t446+t
     #472+t507+t537+t573
      t579 = t577*t260*t181
      t580 = 1/t4
      t581 = t580*t71
      t583 = t581*t327*t163
      A1e0c1p1t1c = 3.D0/2.D0*t579*t583
      A1e0c2p1t1c = 0
      t585 = Mzsq*s
      t588 = 1/(-s+Mzsq)
      t589 = t580*t588
      t590 = abs(DCMPLX(Mzsq))
      t591 = t590**2
      t593 = 1/(t585-t591)
      t595 = t593*t327*t163
      t596 = t589*t595
      A1e0c3p1t1c = -3.D0/2.D0*t585*t579*t596
      t599 = -t24+t26-t31+t32+t35-t36+t54-t56-t58+t59-t61+t63
      t600 = t66-t67+t69-t70-t75+t78+t80-t82+t84-t86+t90-t94
      t602 = t87*t3
      t609 = -t96+t98-32*t602*t33+32*t602*t22+32*t602*t29-t105+t109-t119
     #+t122-t124-t126+t128
      t610 = -t130+t132+t135-t137+t139+t144-t147+t149-t151-t153+t155-t15
     #7+t159
      t613 = -t169-t176+t178+t187-t189-t192+t194+t204-t206-t208+t209-t21
     #1
      t614 = t213+t215-t216+t226-t228-t230+t231-t233+t234+t236-t237+t239
      t616 = t244-t247-t266+t269-t272+t274+t276-t277+t281-t283-t286+t289
      t620 = -t293+t297-t301+t303+t306-t308+t312-t315+24*t284-24*t295+t3
     #16+24*t541-t334
      t624 = t338+t339-t344+t347-t350-t351+t352-t355+t358+t360-t362+t364
      t625 = -t366+t369-t371+t373-t374+t378-t381+t383-t385+t387-t389-t39
     #0
      t627 = t391+t396-t399+t400-t401-t406+t410+t413-t416+t417-t418+t421
      t629 = -t423-t425+t426-t427+t428-t430+t432-24*t279-t433+t436-t438+
     #t440-t442
      t632 = -t445+t450+t452-t454-t457-t460+t462+t464-t466-t467+t468+t47
     #1
      t633 = -t475-t477+t479-t483+t487+t490-t492+t494+t498+t500-t502-t50
     #6
      t635 = -t512+t515+t518-t520+t522+t524+t526-t528-t530+t532+t534-t53
     #6
      t638 = -t540-t543+t545-t548-t551+t553+t555-t558+t562+t566-t569-t57
     #2-32*t602*t64
      A1e0c4p1t1c = -3.D0/2.D0*t585*(t599+t600+t609+t610+t613+t614+t616+
     #t620+t624+t625+t627+t629+t632+t633+t635+t638)*t260*t181*t596
      t648 = conjg(Mzsq)
      t650 = Mzsq+t648-2*s
      t652 = t260*t181
      t655 = 1/(-s+t648)
      t656 = t580*t655
      t658 = t588*t327*t163
      t659 = t656*t658
      A1e0c5p1t1c = -3.D0/2.D0*t650*t577*t652*t659
      A1e0c1p1t2c = A1e0c1p1t1c
      A1e0c2p1t2c = 0
      A1e0c3p1t2c = A1e0c3p1t1c
      A1e0c4p1t2c = A1e0c4p1t1c
      A1e0c5p1t2c = A1e0c5p1t1c
      t662 = j4*t217
      t667 = 2*t166-2*t166*cz2+t166*t28+t28-t331
      t668 = 4*t3
      t669 = s+t668+t248
      A1e0c1p2t1c = 3.D0/2.D0*t662*t667*t669*t583
      A1e0c2p2t1c = 0
      t674 = t585*t662*t667
      t677 = t589*t593
      A1e0c3p2t1c = -3.D0/2.D0*t674*t669*t327*t163*t677
      t681 = -t668+s+t248
      A1e0c4p2t1c = -3.D0/2.D0*t674*t681*t327*t163*t677
      t687 = t650*j4
      A1e0c5p2t1c = -3.D0/2.D0*t687*t217*t667*t669*t659
      A1e0c1p2t2c = A1e0c1p2t1c
      A1e0c2p2t2c = 0
      A1e0c3p2t2c = A1e0c3p2t1c
      A1e0c4p2t2c = A1e0c4p2t1c
      A1e0c5p2t2c = A1e0c5p2t1c
      t693 = j4*t195
      t698 = 2*t184-2*t184*cz1+t184*t27+t27-t264
      A1e0c1p3t1c = 3.D0/2.D0*t693*t698*t669*t652*t581
      A1e0c2p3t1c = 0
      t704 = t585*t693*t698
      A1e0c3p3t1c = -3.D0/2.D0*t704*t669*t260*t181*t677
      A1e0c4p3t1c = -3.D0/2.D0*t704*t681*t260*t181*t677
      t718 = t656*t588
      A1e0c5p3t1c = -3.D0/2.D0*t687*t195*t698*t669*t652*t718
      A1e0c1p3t2c = A1e0c1p3t1c
      A1e0c2p3t2c = 0
      A1e0c3p3t2c = A1e0c3p3t1c
      A1e0c4p3t2c = A1e0c4p3t1c
      A1e0c5p3t2c = A1e0c5p3t1c
      t722 = Mzsq*BCsTh
      t723 = s**2
      t737 = -2*t313+2*t560-2*t567+2*t564+3*t513-3*t510+3*t516-3*t504+3*
     #t485-3*t481-12*t408+12*t404-12*t411-t476+t478-t474-t547
      t745 = -t550+t552-2*t556+2*t496+t461-t459+t463-t456+3*t488-t354-2*
     #t570+t357+2*t310+t554-t361+12*t414-3*t538
      t747 = j4*x1
      t752 = j4*x2
      t757 = j4*t330
      t767 = t359+t470-t747*t363+t747*t365-t747*t268+t747*t265-t752*t348
     #+t752*t345-t752*t332+t752*t336+t757*t29-t757*t33+t757*t22-t757*t64
     #+4*t99*t29-4*t99*t64+9*t64
      t771 = j4*t263
      t780 = t141*t87
      t797 = 9*t33-9*t22-9*t29+t771*t22-t771*t64+t771*t29-t771*t33+4*t99
     #*t22-4*t99*t33+4*t780*t29-4*t780*t33+4*t780*t22-4*t780*t64-2*t752*
     #t174+2*t752*t167-2*t747*t190+2*t747*t185
      t799 = t737+t745+t767+t797
      t803 = t181*t580
      A2e0c1p1t1c = -12*t722*t723*t799*t260*t803*t588*t595
      A2e0c2p1t1c = 0
      t808 = t650*s
      A2e0c3p1t1c = -3*t808*BCsTh*t799*t260*t803*t655*t658
      A2e0c1p1t2c = A2e0c1p1t1c
      A2e0c2p1t2c = 0
      A2e0c3p1t2c = A2e0c3p1t1c
      t816 = t723*j4
      t820 = t667*t327*t163
      A2e0c1p2t1c = -12*t722*t816*t217*t820*t677
      A2e0c2p2t1c = 0
      t824 = BCsTh*j4
      A2e0c3p2t1c = -3*t808*t824*t217*t820*t718
      A2e0c1p2t2c = A2e0c1p2t1c
      A2e0c2p2t2c = 0
      A2e0c3p2t2c = A2e0c3p2t1c
      t833 = t698*t260*t181
      A2e0c1p3t1c = -12*t722*t816*t195*t833*t677
      A2e0c2p3t1c = 0
      A2e0c3p3t1c = -3*t808*t824*t195*t833*t718
      A2e0c1p3t2c = A2e0c1p3t1c
      A2e0c2p3t2c = 0
      A2e0c3p3t2c = A2e0c3p3t1c

      RETURN
      END

      SUBROUTINE CrsecSeqv_msv()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      LOGICAL CIsBad
      LOGICAL IsBad
      DOUBLE COMPLEX spf10
      DOUBLE COMPLEX spf100
      DOUBLE COMPLEX spf101
      DOUBLE COMPLEX spf102
      DOUBLE COMPLEX spf103
      DOUBLE COMPLEX spf104
      DOUBLE COMPLEX spf105
      DOUBLE COMPLEX spf106
      DOUBLE COMPLEX spf107
      DOUBLE COMPLEX spf108
      DOUBLE COMPLEX spf109
      DOUBLE COMPLEX spf11
      DOUBLE COMPLEX spf110
      DOUBLE COMPLEX spf111
      DOUBLE COMPLEX spf112
      DOUBLE COMPLEX spf113
      DOUBLE COMPLEX spf114
      DOUBLE COMPLEX spf115
      DOUBLE COMPLEX spf116
      DOUBLE COMPLEX spf117
      DOUBLE COMPLEX spf118
      DOUBLE COMPLEX spf119
      DOUBLE COMPLEX spf12
      DOUBLE COMPLEX spf120
      DOUBLE COMPLEX spf121
      DOUBLE COMPLEX spf122
      DOUBLE COMPLEX spf123
      DOUBLE COMPLEX spf124
      DOUBLE COMPLEX spf125
      DOUBLE COMPLEX spf126
      DOUBLE COMPLEX spf127
      DOUBLE COMPLEX spf128
      DOUBLE COMPLEX spf129
      DOUBLE COMPLEX spf13
      DOUBLE COMPLEX spf130
      DOUBLE COMPLEX spf131
      DOUBLE COMPLEX spf132
      DOUBLE COMPLEX spf133
      DOUBLE COMPLEX spf134
      DOUBLE COMPLEX spf135
      DOUBLE COMPLEX spf136
      DOUBLE COMPLEX spf137
      DOUBLE COMPLEX spf138
      DOUBLE COMPLEX spf139
      DOUBLE COMPLEX spf14
      DOUBLE COMPLEX spf140
      DOUBLE COMPLEX spf141
      DOUBLE COMPLEX spf142
      DOUBLE COMPLEX spf143
      DOUBLE COMPLEX spf144
      DOUBLE COMPLEX spf145
      DOUBLE COMPLEX spf146
      DOUBLE COMPLEX spf147
      DOUBLE COMPLEX spf148
      DOUBLE COMPLEX spf149
      DOUBLE COMPLEX spf15
      DOUBLE COMPLEX spf150
      DOUBLE COMPLEX spf151
      DOUBLE COMPLEX spf152
      DOUBLE COMPLEX spf153
      DOUBLE COMPLEX spf154
      DOUBLE COMPLEX spf155
      DOUBLE COMPLEX spf156
      DOUBLE COMPLEX spf157
      DOUBLE COMPLEX spf158
      DOUBLE COMPLEX spf159
      DOUBLE COMPLEX spf16
      DOUBLE COMPLEX spf160
      DOUBLE COMPLEX spf161
      DOUBLE COMPLEX spf162
      DOUBLE COMPLEX spf163
      DOUBLE COMPLEX spf164
      DOUBLE COMPLEX spf165
      DOUBLE COMPLEX spf166
      DOUBLE COMPLEX spf167
      DOUBLE COMPLEX spf168
      DOUBLE COMPLEX spf169
      DOUBLE COMPLEX spf17
      DOUBLE COMPLEX spf170
      DOUBLE COMPLEX spf171
      DOUBLE COMPLEX spf172
      DOUBLE COMPLEX spf173
      DOUBLE COMPLEX spf174
      DOUBLE COMPLEX spf175
      DOUBLE COMPLEX spf176
      DOUBLE COMPLEX spf177
      DOUBLE COMPLEX spf178
      DOUBLE COMPLEX spf179
      DOUBLE COMPLEX spf18
      DOUBLE COMPLEX spf180
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
      DOUBLE COMPLEX spf38
      DOUBLE COMPLEX spf39
      DOUBLE COMPLEX spf4
      DOUBLE COMPLEX spf40
      DOUBLE COMPLEX spf41
      DOUBLE COMPLEX spf42
      DOUBLE COMPLEX spf43
      DOUBLE COMPLEX spf44
      DOUBLE COMPLEX spf45
      DOUBLE COMPLEX spf46
      DOUBLE COMPLEX spf47
      DOUBLE COMPLEX spf48
      DOUBLE COMPLEX spf49
      DOUBLE COMPLEX spf5
      DOUBLE COMPLEX spf50
      DOUBLE COMPLEX spf51
      DOUBLE COMPLEX spf52
      DOUBLE COMPLEX spf53
      DOUBLE COMPLEX spf54
      DOUBLE COMPLEX spf55
      DOUBLE COMPLEX spf56
      DOUBLE COMPLEX spf57
      DOUBLE COMPLEX spf58
      DOUBLE COMPLEX spf59
      DOUBLE COMPLEX spf60
      DOUBLE COMPLEX spf61
      DOUBLE COMPLEX spf62
      DOUBLE COMPLEX spf63
      DOUBLE COMPLEX spf64
      DOUBLE COMPLEX spf65
      DOUBLE COMPLEX spf66
      DOUBLE COMPLEX spf67
      DOUBLE COMPLEX spf68
      DOUBLE COMPLEX spf69
      DOUBLE COMPLEX spf7
      DOUBLE COMPLEX spf70
      DOUBLE COMPLEX spf71
      DOUBLE COMPLEX spf72
      DOUBLE COMPLEX spf73
      DOUBLE COMPLEX spf74
      DOUBLE COMPLEX spf75
      DOUBLE COMPLEX spf76
      DOUBLE COMPLEX spf77
      DOUBLE COMPLEX spf78
      DOUBLE COMPLEX spf79
      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX spf80
      DOUBLE COMPLEX spf81
      DOUBLE COMPLEX spf82
      DOUBLE COMPLEX spf83
      DOUBLE COMPLEX spf84
      DOUBLE COMPLEX spf85
      DOUBLE COMPLEX spf86
      DOUBLE COMPLEX spf87
      DOUBLE COMPLEX spf88
      DOUBLE COMPLEX spf89
      DOUBLE COMPLEX spf9
      DOUBLE COMPLEX spf90
      DOUBLE COMPLEX spf91
      DOUBLE COMPLEX spf92
      DOUBLE COMPLEX spf93
      DOUBLE COMPLEX spf94
      DOUBLE COMPLEX spf95
      DOUBLE COMPLEX spf96
      DOUBLE COMPLEX spf97
      DOUBLE COMPLEX spf98
      DOUBLE COMPLEX spf99
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
      DOUBLE COMPLEX vI110pc1
      DOUBLE COMPLEX vI111pc1
      DOUBLE COMPLEX vI112pc1
      DOUBLE COMPLEX vI113pc1
      DOUBLE COMPLEX vI114pc1
      DOUBLE COMPLEX vI115pc1
      DOUBLE COMPLEX vI116pc1
      DOUBLE COMPLEX vI117pc1
      DOUBLE COMPLEX vI118pc1
      DOUBLE COMPLEX vI119pc1
      DOUBLE COMPLEX vI120pc1
      DOUBLE COMPLEX vI121pc1
      DOUBLE COMPLEX vI122pc1
      DOUBLE COMPLEX vI123pc1
      DOUBLE COMPLEX vI124pc1
      DOUBLE COMPLEX vI125pc1
      DOUBLE COMPLEX vI126pc1
      DOUBLE COMPLEX vI127pc1
      DOUBLE COMPLEX vI128pc1
      DOUBLE COMPLEX vI129pc1
      DOUBLE COMPLEX vI130pc1
      DOUBLE COMPLEX vI131pc1
      DOUBLE COMPLEX vI132pc1
      DOUBLE COMPLEX vI133pc1
      DOUBLE COMPLEX vI134pc1
      DOUBLE COMPLEX vI135pc1
      DOUBLE COMPLEX vI136pc1
      DOUBLE COMPLEX vI137pc1
      DOUBLE COMPLEX vI138pc1
      DOUBLE COMPLEX vI139pc1
      DOUBLE COMPLEX vI140pc1
      DOUBLE COMPLEX vI141pc1
      DOUBLE COMPLEX vI142pc1
      DOUBLE COMPLEX vI143pc1
      DOUBLE COMPLEX vI144pc1
      DOUBLE COMPLEX vI145pc1
      DOUBLE COMPLEX vI146pc1
      DOUBLE COMPLEX vI147pc1
      DOUBLE COMPLEX vI148pc1
      DOUBLE COMPLEX vI149pc1
      DOUBLE COMPLEX vI150pc1
      DOUBLE COMPLEX vI151pc1
      DOUBLE COMPLEX vI152pc1
      DOUBLE COMPLEX vI153pc1
      DOUBLE COMPLEX vI154pc1
      DOUBLE COMPLEX vI155pc1
      DOUBLE COMPLEX vI156pc1
      DOUBLE COMPLEX vI157pc1
      DOUBLE COMPLEX vI158pc1
      DOUBLE COMPLEX vI159pc1
      DOUBLE COMPLEX vI160pc1
      DOUBLE COMPLEX vI161pc1
      DOUBLE COMPLEX vI162pc1
      DOUBLE COMPLEX vI163pc1
      DOUBLE COMPLEX vI164pc1
      DOUBLE COMPLEX vI165pc1
      DOUBLE COMPLEX vI166pc1
      DOUBLE COMPLEX vI167pc1
      DOUBLE COMPLEX vI168pc1
      DOUBLE COMPLEX vI169pc1
      DOUBLE COMPLEX vI170pc1
      DOUBLE COMPLEX vI171pc1
      DOUBLE COMPLEX vI172pc1
      DOUBLE COMPLEX vI173pc1
      DOUBLE COMPLEX vI174pc1
      DOUBLE COMPLEX vI175pc1
      DOUBLE COMPLEX vI176pc1
      DOUBLE COMPLEX vI177pc1
      DOUBLE COMPLEX vI178pc1
      DOUBLE COMPLEX vI179pc1
      DOUBLE COMPLEX vI180pc1
      DOUBLE COMPLEX vI181pc1
      DOUBLE COMPLEX vI182pc1
      DOUBLE COMPLEX vI183pc1
      DOUBLE COMPLEX vI184pc1
      DOUBLE COMPLEX vI185pc1
      DOUBLE COMPLEX vI186pc1
      DOUBLE COMPLEX vI187pc1
      DOUBLE COMPLEX vI188pc1
      DOUBLE COMPLEX vI189pc1
      DOUBLE COMPLEX vI190pc1
      DOUBLE COMPLEX vI191pc1
      DOUBLE COMPLEX vI192pc1
      DOUBLE COMPLEX vI193pc1
      DOUBLE COMPLEX vI194pc1
      DOUBLE COMPLEX vI195pc1
      DOUBLE COMPLEX vI196pc1
      DOUBLE COMPLEX vI197pc1
      DOUBLE COMPLEX vI198pc1
      DOUBLE COMPLEX vI199pc1
      DOUBLE COMPLEX vI200pc1
      DOUBLE COMPLEX vI201pc1
      DOUBLE COMPLEX vI202pc1
      DOUBLE COMPLEX vI203pc1
      DOUBLE COMPLEX vI204pc1
      DOUBLE COMPLEX vI205pc1
      DOUBLE COMPLEX vI206pc1
      DOUBLE COMPLEX vI207pc1
      DOUBLE COMPLEX vI208pc1
      DOUBLE COMPLEX vI209pc1
      DOUBLE COMPLEX vI210pc1
      DOUBLE COMPLEX vI211pc1
      DOUBLE COMPLEX vI212pc1
      DOUBLE COMPLEX vI213pc1
      DOUBLE COMPLEX vI214pc1
      DOUBLE COMPLEX vI215pc1
      DOUBLE COMPLEX vI216pc1
      DOUBLE COMPLEX vI217pc1
      DOUBLE COMPLEX vI218pc1
      DOUBLE COMPLEX vI219pc1
      DOUBLE COMPLEX vI220pc1
      DOUBLE COMPLEX vI221pc1
      DOUBLE COMPLEX vI222pc1
      DOUBLE COMPLEX vI223pc1
      DOUBLE COMPLEX vI224pc1
      DOUBLE COMPLEX vI225pc1
      DOUBLE COMPLEX vI226pc1
      DOUBLE COMPLEX vI227pc1
      DOUBLE COMPLEX vI228pc1
      DOUBLE COMPLEX vI229pc1
      DOUBLE COMPLEX vI230pc1
      DOUBLE COMPLEX vI231pc1
      DOUBLE COMPLEX vI232pc1
      DOUBLE COMPLEX vI233pc1
      DOUBLE COMPLEX vI234pc1
      DOUBLE COMPLEX vI235pc1
      DOUBLE COMPLEX vI236pc1
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
      DOUBLE COMPLEX vV100
      DOUBLE COMPLEX vV101
      DOUBLE COMPLEX vV102
      DOUBLE COMPLEX vV103
      DOUBLE COMPLEX vV104
      DOUBLE COMPLEX vV105
      DOUBLE COMPLEX vV106
      DOUBLE COMPLEX vV107
      DOUBLE COMPLEX vV108
      DOUBLE COMPLEX vV109
      DOUBLE COMPLEX vV110
      DOUBLE COMPLEX vV111
      DOUBLE COMPLEX vV112
      DOUBLE COMPLEX vV113
      DOUBLE COMPLEX vV114
      DOUBLE COMPLEX vV115
      DOUBLE COMPLEX vV116
      DOUBLE COMPLEX vV117
      DOUBLE COMPLEX vV118
      DOUBLE COMPLEX vV119
      DOUBLE COMPLEX vV120
      DOUBLE COMPLEX vV121
      DOUBLE COMPLEX vV122
      DOUBLE COMPLEX vV123
      DOUBLE COMPLEX vV124
      DOUBLE COMPLEX vV125
      DOUBLE COMPLEX vV126
      DOUBLE COMPLEX vV127
      DOUBLE COMPLEX vV128
      DOUBLE COMPLEX vV129
      DOUBLE COMPLEX vV130
      DOUBLE COMPLEX vV131
      DOUBLE COMPLEX vV132
      DOUBLE COMPLEX vV133
      DOUBLE COMPLEX vV134
      DOUBLE COMPLEX vV135
      DOUBLE COMPLEX vV136
      DOUBLE COMPLEX vV137
      DOUBLE COMPLEX vV138
      DOUBLE COMPLEX vV139
      DOUBLE COMPLEX vV140
      DOUBLE COMPLEX vV141
      DOUBLE COMPLEX vV142
      DOUBLE COMPLEX vV143
      DOUBLE COMPLEX vV144
      DOUBLE COMPLEX vV145
      DOUBLE COMPLEX vV146
      DOUBLE COMPLEX vV147
      DOUBLE COMPLEX vV148
      DOUBLE COMPLEX vV149
      DOUBLE COMPLEX vV150
      DOUBLE COMPLEX vV151
      DOUBLE COMPLEX vV152
      DOUBLE COMPLEX vV153
      DOUBLE COMPLEX vV154
      DOUBLE COMPLEX vV155
      DOUBLE COMPLEX vV156
      DOUBLE COMPLEX vV157
      DOUBLE COMPLEX vV158
      DOUBLE COMPLEX vV159
      DOUBLE COMPLEX vV160
      DOUBLE COMPLEX vV161
      DOUBLE COMPLEX vV162
      DOUBLE COMPLEX vV163
      DOUBLE COMPLEX vV164
      DOUBLE COMPLEX vV165
      DOUBLE COMPLEX vV166
      DOUBLE COMPLEX vV167
      DOUBLE COMPLEX vV168
      DOUBLE COMPLEX vV169
      DOUBLE COMPLEX vV170
      DOUBLE COMPLEX vV171
      DOUBLE COMPLEX vV172
      DOUBLE COMPLEX vV173
      DOUBLE COMPLEX vV174
      DOUBLE COMPLEX vV175
      DOUBLE COMPLEX vV176
      DOUBLE COMPLEX vV177
      DOUBLE COMPLEX vV178
      DOUBLE COMPLEX vV179
      DOUBLE COMPLEX vV180
      DOUBLE COMPLEX vV181
      DOUBLE COMPLEX vV182
      DOUBLE COMPLEX vV183
      DOUBLE COMPLEX vV184
      DOUBLE COMPLEX vV185
      DOUBLE COMPLEX vV186
      DOUBLE COMPLEX vV187
      DOUBLE COMPLEX vV188
      DOUBLE COMPLEX vV189
      DOUBLE COMPLEX vV190
      DOUBLE COMPLEX vV191
      DOUBLE COMPLEX vV192
      DOUBLE COMPLEX vV193
      DOUBLE COMPLEX vV194
      DOUBLE COMPLEX vV55
      DOUBLE COMPLEX vV56
      DOUBLE COMPLEX vV57
      DOUBLE COMPLEX vV58
      DOUBLE COMPLEX vV59
      DOUBLE COMPLEX vV60
      DOUBLE COMPLEX vV61
      DOUBLE COMPLEX vV62
      DOUBLE COMPLEX vV63
      DOUBLE COMPLEX vV64
      DOUBLE COMPLEX vV65
      DOUBLE COMPLEX vV66
      DOUBLE COMPLEX vV67
      DOUBLE COMPLEX vV68
      DOUBLE COMPLEX vV69
      DOUBLE COMPLEX vV70
      DOUBLE COMPLEX vV71
      DOUBLE COMPLEX vV72
      DOUBLE COMPLEX vV73
      DOUBLE COMPLEX vV74
      DOUBLE COMPLEX vV75
      DOUBLE COMPLEX vV76
      DOUBLE COMPLEX vV77
      DOUBLE COMPLEX vV78
      DOUBLE COMPLEX vV79
      DOUBLE COMPLEX vV80
      DOUBLE COMPLEX vV81
      DOUBLE COMPLEX vV82
      DOUBLE COMPLEX vV83
      DOUBLE COMPLEX vV84
      DOUBLE COMPLEX vV85
      DOUBLE COMPLEX vV86
      DOUBLE COMPLEX vV87
      DOUBLE COMPLEX vV88
      DOUBLE COMPLEX vV89
      DOUBLE COMPLEX vV90
      DOUBLE COMPLEX vV91
      DOUBLE COMPLEX vV92
      DOUBLE COMPLEX vV93
      DOUBLE COMPLEX vV94
      DOUBLE COMPLEX vV95
      DOUBLE COMPLEX vV96
      DOUBLE COMPLEX vV97
      DOUBLE COMPLEX vV98
      DOUBLE COMPLEX vV99
      DOUBLE COMPLEX DENSPENCE
      DOUBLE COMPLEX IM
      PARAMETER (IM=(0.D0,1.D0))
      DOUBLE PRECISION epsMl
      PARAMETER (epsMl=1d-12)
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      t1 = -s+Ieps
      t2 = 1/t1
      spf10 = log(DCMPLX(-t2))
      t3 = sqrt(s)
      t4 = Ml**2
      t5 = 4*t4
      t6 = 4*Ieps
      t7 = s-t5+t6
      t8 = sqrt(t7)
      t9 = t3*t8
      t10 = s+t9
      t11 = 1/s
      t12 = t10*t11
      spf27 = log(DCMPLX(-t12))
      t13 = t1*t4
      t14 = s**2
      t15 = 1/t14
      t17 = 2*t13*t15
      t18 = -s+t9
      t19 = t18*t11
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf28 = t17
      else
      spf28 = t19
      end if
      if ((ABS(DBLE(spf28)).eq.0d0).or.(CIsBad(spf28))) spf28 = t17
      spf28 = log(spf28)
      t20 = t1*t11
      t21 = t4*s
      t22 = 4*t21
      t23 = t4*Ieps
      t24 = 4*t23
      t26 = sqrt(t14-t22+t24)
      t27 = s+t26
      t31 = 2*t20*t4/t27
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf67 = t31
      else
      spf67 = t18/t10
      end if
      if ((ABS(DBLE(spf67)).eq.0d0).or.(CIsBad(spf67))) spf67 = t31
      spf67 = log(spf67)
      t35 = 1/t4
      t38 = t27*t2*s*t35/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf68 = t38
      else
      spf68 = t10/t18
      end if
      if ((ABS(DBLE(spf68)).eq.0d0).or.(CIsBad(spf68))) spf68 = t38
      spf68 = log(spf68)
      t42 = log(DCMPLX(2.D0))
      t44 = spf27/4
      t46 = sqrt(t7*t11)
      t50 = spf28/4
      vV55 = 1.D0/2.D0+spf68/8+t42/2-t44+spf67*t46/8+spf10/4-t50+spf67/8
     #-spf68*t46/8
      t54 = -t4+Ieps
      spf19 = log(DCMPLX(-t54*t11))
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf29 = -t17
      else
      spf29 = -t19
      end if
      if ((ABS(DBLE(spf29)).eq.0d0).or.(CIsBad(spf29))) spf29 = -t17
      spf29 = log(spf29)
      spf30 = log(DCMPLX(t12))
      t58 = t20*t4/t26
      t59 = 2*t58
      t60 = 1/t3
      t62 = 1/t8
      t63 = t18*t60*t62
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf59 = t59
      else
      spf59 = t63
      end if
      if ((ABS(DBLE(spf59)).eq.0d0).or.(CIsBad(spf59))) spf59 = t59
      spf59 = log(spf59)
      t65 = t10*t60*t62
      spf60 = log(DCMPLX(t65))
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf65 = t58
      else
      spf65 = t63/2
      end if
      if (CIsBad(spf65)) spf65 = t58
      spf65 = DENSPENCE(spf65,spcZip)
      spf66 = DENSPENCE(DCMPLX(t65/2),spcZip)
      t68 = spf10*spf67
      t71 = t46*spf29
      t74 = t46*t42
      t77 = spf28**2
      t87 = 1+t68*t46/8+t71*spf60/8+t74*spf59/4-t44-t46*t77/16-spf28*t42
     #/4+t68/8+spf10*t42/2+spf10/2+t71/4
      t88 = spf29**2
      t93 = t46*spf30
      t95 = spf30**2
      t108 = spf10*spf68
      t110 = t46*spf27
      t112 = -t46*t88/16+spf29*t42/4-t93/4+t46*t95/16+spf30*t42/4+spf65*
     #t46/4-spf66*t46/4-spf29*spf60/8+spf28*spf60/8+t108/8-t110/4
      t114 = spf27**2
      t127 = t46*spf28
      t131 = spf10**2
      t133 = PI**2
      t135 = t46*t114/16-spf27*t42/4-spf10*spf27/4-spf30*spf59/8+spf27*s
     #pf59/8-spf10*spf28/4+t127/4+t114/16+t77/16+t131/8+t133/24
      t136 = spf19**2
      t152 = t136/8-t95/16-t88/16+t42-spf29/4-spf30/4-t50-t74*spf60/4+t1
     #27*spf60/8-t108*t46/8-t110*spf59/8-t93*spf59/8
      vV56 = t87+t112+t135+t152
      t154 = 1/t54
      spf14 = log(DCMPLX(-t154))
      vV57 = t4*(spf14+1)/6
      t157 = spf14**2
      vV58 = t4*(t133+3*t157+6*spf14+6)/36
      spf4 = log(DCMPLX(t35))
      t162 = CsTh*t14
      t163 = t4*CsTh
      t165 = 2*t163*s
      t166 = 2*t23
      t167 = -t14-t162+t165+t166
      t170 = t167*t35*t2/2
      t171 = s*BCsTh
      t172 = 2*Ieps
      t173 = s+t171-t172
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf40 = t170
      else
      spf40 = -t173*t154/2
      end if
      if ((ABS(DBLE(spf40)).eq.0d0).or.(CIsBad(spf40))) spf40 = t170
      spf40 = log(spf40)
      t176 = 1+BCsTh
      t177 = s*t176
      t178 = 2*t4
      t179 = spf40*t4
      t181 = -t177/2+t4
      t182 = spf40*t181
      t183 = spf4*t181
      t185 = 1/t181
      vV59 = (-t177+t178+t179-t182+t183)*t185/12
      t187 = 2*t21
      t188 = t4**2
      t189 = CsTh*t188
      t190 = 2*t189
      t191 = -t14-t162+t165+t187+t190
      t192 = 1/t167
      t193 = t191*t192
      t194 = s+t171-t178
      t195 = 1/t173
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf50 = t193
      else
      spf50 = t194*t195
      end if
      if (CIsBad(spf50)) spf50 = t193
      spf50 = DENSPENCE(spf50,spcZip)
      t197 = 1/t191
      t199 = 2*t13*t197
      t200 = 1/t194
      t202 = 2*t54*t200
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf54 = -t199
      else
      spf54 = t202
      end if
      if ((ABS(DBLE(spf54)).eq.0d0).or.(CIsBad(spf54))) spf54 = -t199
      spf54 = log(spf54)
      t203 = t167*t197
      t204 = t173*t200
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf62 = -t203
      else
      spf62 = -t204
      end if
      if ((ABS(DBLE(spf62)).eq.0d0).or.(CIsBad(spf62))) spf62 = -t203
      spf62 = log(spf62)
      t206 = 24*t4
      t213 = spf54**2
      t216 = spf62**2
      t230 = spf4*spf40
      t236 = spf4**2
      t239 = -12*t177+t206+6*t4*spf50+12*t179+6*t4*spf54*spf40+3*t4*t213
     #-3*t4*t216-6*t181*spf50-12*t182-6*t181*spf54*spf40-3*t181*t213+3*t
     #181*t216+12*t183+6*t230*t4-6*t230*t181+t133*t181+3*t236*t181
      vV60 = t239*t185/144
      t241 = t14-t162+t165-t166
      t244 = t241*t35*t2/2
      t245 = -s+t171+t172
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf39 = -t244
      else
      spf39 = t245*t154/2
      end if
      if ((ABS(DBLE(spf39)).eq.0d0).or.(CIsBad(spf39))) spf39 = -t244
      spf39 = log(spf39)
      t248 = 1-BCsTh
      t249 = s*t248
      t250 = spf39*t4
      t252 = -t249/2+t4
      t253 = spf39*t252
      t254 = spf4*t252
      t256 = 1/t252
      vV61 = (-t249+t178+t250-t253+t254)*t256/12
      t258 = t14-t162+t165-t187+t190
      t259 = 1/t241
      t260 = t258*t259
      t261 = -s+t171+t178
      t262 = 1/t245
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf49 = t260
      else
      spf49 = t261*t262
      end if
      if (CIsBad(spf49)) spf49 = t260
      spf49 = DENSPENCE(spf49,spcZip)
      t264 = 1/t258
      t266 = 2*t13*t264
      t267 = 1/t261
      t269 = 2*t54*t267
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf53 = t266
      else
      spf53 = -t269
      end if
      if ((ABS(DBLE(spf53)).eq.0d0).or.(CIsBad(spf53))) spf53 = t266
      spf53 = log(spf53)
      t270 = t241*t264
      t271 = t245*t267
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf61 = -t270
      else
      spf61 = -t271
      end if
      if ((ABS(DBLE(spf61)).eq.0d0).or.(CIsBad(spf61))) spf61 = -t270
      spf61 = log(spf61)
      t279 = spf53**2
      t282 = spf61**2
      t296 = spf4*spf39
      t304 = -12*t249+t206+6*t4*spf49+12*t250+6*t4*spf53*spf39+3*t4*t279
     #-3*t4*t282-6*t252*spf49-12*t253-6*t252*spf53*spf39-3*t252*t279+3*t
     #252*t282+12*t254+6*t296*t4-6*t296*t252+t133*t252+3*t236*t252
      vV62 = t304*t256/144
      t312 = sqrt((t14+5*t23-t22)/(t14+t23))
      t317 = 2*t1*t15*t4/(1+t312)
      t320 = s+Ieps
      t321 = 1/t320
      t323 = sqrt((s+5*Ieps-t5)*t321)
      t327 = (-1+t323)/(1+t323)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf69 = t317
      else
      spf69 = t327
      end if
      if ((ABS(DBLE(spf69)).eq.0d0).or.(CIsBad(spf69))) spf69 = t317
      spf69 = log(spf69)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf70 = t317
      else
      spf70 = t327
      end if
      if (CIsBad(spf70)) spf70 = t317
      spf70 = DENSPENCE(spf70,spcZip)
      t330 = spf69**2
      vV63 = (4*t133+12*spf70+3*t330)/t323*t11/72
      spf22 = log(DCMPLX(t320*t154))
      t340 = (-t14-t162+t24+t165)*t35*t2/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf38 = t340
      else
      spf38 = -(-t6+s+t171)*t154/2
      end if
      if ((ABS(DBLE(spf38)).eq.0d0).or.(CIsBad(spf38))) spf38 = t340
      spf38 = log(spf38)
      t349 = 2*spf14*spf22
      t350 = 2*t157
      t353 = 1/t176
      vV64 = (-4*spf38*spf22+t133+4*spf14*spf38+t349-t350)*t15*t353/4
      t356 = 2*spf14
      vV65 = (2*spf38+spf22-t356)*t15*t353/2
      t363 = (t14-t162-t24+t165)*t35*t2/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf37 = -t363
      else
      spf37 = (t6-s+t171)*t154/2
      end if
      if ((ABS(DBLE(spf37)).eq.0d0).or.(CIsBad(spf37))) spf37 = -t363
      spf37 = log(spf37)
      t373 = 1/t248
      vV66 = -(4*spf37*spf22-t133-4*spf14*spf37-t349+t350)*t15*t373/4
      vV67 = (2*spf37+spf22-t356)*t15*t373/2
      vV68 = vV57
      vV69 = vV58
      spf11 = log(DCMPLX(-t321))
      t381 = spf11**2
      vV70 = 1+spf11/2+t381/8
      vV71 = vV70
      vV72 = vV55
      vV73 = vV56
      vV74 = vV68
      vV75 = vV69
      vV76 = vV71
      vV77 = vV72
      vV78 = vV73
      vV79 = vV74
      vV80 = vV75
      vV81 = vV76
      t383 = -Mzsq+Ieps
      t384 = 1/t383
      spf13 = log(DCMPLX(-t384))
      vV82 = Mzsq*(spf13+1)/6
      vV83 = vV59
      vV84 = vV60
      vV85 = vV61
      vV86 = vV62
      t387 = Mzsq*s
      t388 = -t387+t23
      t389 = t388*t11
      t390 = CsTh*s
      t391 = 2*t163
      t392 = -s-t390+t391
      t393 = 1/t392
      t395 = 2*t389*t393
      t396 = t383*t11
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf26 = -t395
      else
      spf26 = 2*t396*t353
      end if
      if ((ABS(DBLE(spf26)).eq.0d0).or.(CIsBad(spf26))) spf26 = -t395
      spf26 = log(spf26)
      t399 = s*t392
      t400 = 1/t388
      t402 = t399*t400/2
      t404 = t177*t384/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf33 = t402
      else
      spf33 = -t404
      end if
      if ((ABS(DBLE(spf33)).eq.0d0).or.(CIsBad(spf33))) spf33 = t402
      spf33 = log(spf33)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf34 = -t402
      else
      spf34 = t404
      end if
      if ((ABS(DBLE(spf34)).eq.0d0).or.(CIsBad(spf34))) spf34 = -t402
      spf34 = log(spf34)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf42 = t199
      else
      spf42 = -t202
      end if
      if ((ABS(DBLE(spf42)).eq.0d0).or.(CIsBad(spf42))) spf42 = t199
      spf42 = log(spf42)
      t405 = 2*Mzsq
      t406 = -t405+t172+s+t171
      spf46 = log(DCMPLX(2*t383/t406))
      t411 = t406*t384/2
      spf48 = log(DCMPLX(t411))
      spf51 = DENSPENCE(DCMPLX(t411),spcZip)
      t412 = Mzsq**2
      t413 = t4*Mzsq
      t414 = 4*t413
      t416 = sqrt(t412-t414+t24)
      t417 = Mzsq-t416
      if (ABS(DBLE(Ml**2/DBLE(Mzsq))).lt.epsMl) then
      spf55 = 2
      else
      spf55 = t417*t35
      end if
      if ((ABS(DBLE(spf55)).eq.0d0).or.(CIsBad(spf55))) spf55 = 2
      spf55 = log(spf55)
      t419 = Mzsq+t416
      spf56 = log(DCMPLX(t419*t35))
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf58 = t203
      else
      spf58 = t204
      end if
      if ((ABS(DBLE(spf58)).eq.0d0).or.(CIsBad(spf58))) spf58 = t203
      spf58 = log(spf58)
      t421 = -Mzsq-t416+t178
      spf63 = log(DCMPLX(-t421*t35))
      t423 = 1/Mzsq
      t428 = t423*(Ieps-DBLE(Mzsq))*t4/DBLE(Mzsq)
      t429 = 2*t428
      t430 = -Mzsq+t416+t178
      if (ABS(DBLE(Ml**2/DBLE(Mzsq))**2).lt.epsMl) then
      spf64 = -t429
      else
      spf64 = -t430*t35
      end if
      if ((ABS(DBLE(spf64)).eq.0d0).or.(CIsBad(spf64))) spf64 = -t429
      spf64 = log(spf64)
      t433 = t14*Mzsq
      t434 = 2*t433
      t436 = 2*t433*CsTh
      t437 = t14*s
      t438 = CsTh**2
      t439 = t437*t438
      t441 = 2*CsTh*t437
      t442 = t4*t14
      t444 = 4*t442*CsTh
      t445 = t413*s
      t446 = 4*t445
      t447 = Mzsq*CsTh
      t449 = 4*t21*t447
      t451 = 4*t442*t438
      t453 = 4*t189*s
      t454 = t188*Mzsq
      t456 = 4*t454*CsTh
      t457 = t188*Ieps
      t458 = 4*t457
      t459 = t434+t436-t439-t437-t441+t444-t446-t449+t451+t453-t456+t458
      t460 = 1/t459
      t461 = t11*t460
      t463 = 2*t388*t191*t461
      t465 = 2*t387
      t466 = t387*BCsTh
      t467 = 2*t466
      t468 = t14*BCsTh
      t469 = 2*t468
      t470 = BCsTh**2
      t471 = t14*t470
      t472 = -t465-t467+t414+t14+t469+t471-t24
      t473 = 1/t472
      t475 = 2*t383*t194*t473
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf74 = t463
      else
      spf74 = t475
      end if
      if (CIsBad(spf74)) spf74 = t463
      spf74 = DENSPENCE(spf74,spcZip)
      t476 = -t465+t166+t14+t468
      t478 = t191*t476*t461
      t480 = t194*t406*t473
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf77 = t478
      else
      spf77 = t480
      end if
      if (CIsBad(spf77)) spf77 = t478
      spf77 = DENSPENCE(spf77,spcZip)
      t483 = 2*t13*t392*t460
      t484 = t54*s
      t485 = t176*t473
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf86 = -t483
      else
      spf86 = -2*t484*t485
      end if
      if ((ABS(DBLE(spf86)).eq.0d0).or.(CIsBad(spf86))) spf86 = -t483
      spf86 = log(spf86)
      t488 = -t405+s+t390
      t489 = 1/t488
      t490 = t389*t489
      t491 = 2*t490
      t492 = t4*t383
      t493 = t416*s
      t494 = t493*BCsTh
      t495 = -t414+t24+t387+t466-t493-t494
      t496 = 1/t495
      t497 = t492*t496
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf87 = t491
      else
      spf87 = 4*t497
      end if
      if (CIsBad(spf87)) spf87 = t491
      spf87 = DENSPENCE(spf87,spcZip)
      t499 = t4*t388
      t500 = t14*t412
      t501 = t500*CsTh
      t502 = t412*t4
      t504 = 2*t502*s
      t505 = t412*CsTh
      t507 = 2*t21*t505
      t508 = t442*Mzsq
      t509 = t442*t447
      t511 = 2*t457*Mzsq
      t512 = t188*t412
      t514 = 2*t512*CsTh
      t515 = t188*s
      t516 = t515*Ieps
      t517 = t188*t14
      t518 = Ieps*s
      t519 = t189*t518
      t520 = t189*t14
      t522 = 2*t454*t390
      t523 = t500+t501-t504-t507-t508-t509+t511-t514+t516-t517+t519-t520
     #+t522
      t524 = 1/t523
      t525 = Mzsq*t524
      t526 = t499*t525
      t527 = 2*t526
      t528 = -t414+t24+t387+t466+t493+t494
      t529 = 1/t528
      t530 = t492*t529
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf88 = t527
      else
      spf88 = 4*t530
      end if
      if (CIsBad(spf88)) spf88 = t527
      spf88 = DENSPENCE(spf88,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf94 = t463
      else
      spf94 = t475
      end if
      if ((ABS(DBLE(spf94)).eq.0d0).or.(CIsBad(spf94))) spf94 = t463
      spf94 = log(spf94)
      t533 = t459*t393*t197
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf97 = -t533
      else
      spf97 = t472*t11*t353*t200
      end if
      if ((ABS(DBLE(spf97)).eq.0d0).or.(CIsBad(spf97))) spf97 = -t533
      spf97 = log(spf97)
      t538 = t167*t392*t460
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf99 = -t538
      else
      spf99 = t173*s*t485
      end if
      if ((ABS(DBLE(spf99)).eq.0d0).or.(CIsBad(spf99))) spf99 = -t538
      spf99 = log(spf99)
      t542 = t392*t191*t460
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf101 = -t542
      else
      spf101 = t177*t194*t473
      end if
      if ((ABS(DBLE(spf101)).eq.0d0).or.(CIsBad(spf101))) spf101 = -t542
      spf101 = log(spf101)
      t546 = t476*t11*t489
      t547 = t4*t406
      t548 = t547*t496
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf102 = t546
      else
      spf102 = 2*t548
      end if
      if (CIsBad(spf102)) spf102 = t546
      spf102 = DENSPENCE(spf102,spcZip)
      t551 = t4*t476*t525
      t552 = t547*t529
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf103 = t551
      else
      spf103 = 2*t552
      end if
      if (CIsBad(spf103)) spf103 = t551
      spf103 = DENSPENCE(spf103,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf107 = t478
      else
      spf107 = t480
      end if
      if ((ABS(DBLE(spf107)).eq.0d0).or.(CIsBad(spf107))) spf107 = t478
      spf107 = log(spf107)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf114 = t490
      else
      spf114 = 2*t497
      end if
      if ((ABS(DBLE(spf114)).eq.0d0).or.(CIsBad(spf114))) spf114 = t490
      spf114 = log(spf114)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf115 = t526
      else
      spf115 = 2*t530
      end if
      if ((ABS(DBLE(spf115)).eq.0d0).or.(CIsBad(spf115))) spf115 = t526
      spf115 = log(spf115)
      t557 = 2*t488*t393
      t559 = t353*t35
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf118 = -t557
      else
      spf118 = t495*t11*t559
      end if
      if ((ABS(DBLE(spf118)).eq.0d0).or.(CIsBad(spf118))) spf118 = -t557
      spf118 = log(spf118)
      t565 = 2*t423*t523*t11*t393*t35
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf119 = -t565
      else
      spf119 = t528*t11*t559
      end if
      if ((ABS(DBLE(spf119)).eq.0d0).or.(CIsBad(spf119))) spf119 = -t565
      spf119 = log(spf119)
      t568 = t392*t489
      t569 = t568/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf122 = -t569
      else
      spf122 = t177*t4*t496
      end if
      if ((ABS(DBLE(spf122)).eq.0d0).or.(CIsBad(spf122))) spf122 = -t569
      spf122 = log(spf122)
      t574 = t399*t413*t524/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf123 = -t574
      else
      spf123 = t177*t4*t529
      end if
      if ((ABS(DBLE(spf123)).eq.0d0).or.(CIsBad(spf123))) spf123 = -t574
      spf123 = log(spf123)
      t577 = t546/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf128 = t577
      else
      spf128 = t548
      end if
      if ((ABS(DBLE(spf128)).eq.0d0).or.(CIsBad(spf128))) spf128 = t577
      spf128 = log(spf128)
      t578 = t551/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf129 = t578
      else
      spf129 = t552
      end if
      if ((ABS(DBLE(spf129)).eq.0d0).or.(CIsBad(spf129))) spf129 = t578
      spf129 = log(spf129)
      t579 = t417*s
      t580 = t176*t496
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf142 = -t568
      else
      spf142 = t579*t580
      end if
      if ((ABS(DBLE(spf142)).eq.0d0).or.(CIsBad(spf142))) spf142 = -t568
      spf142 = log(spf142)
      t582 = s*t412
      t584 = sqrt(t582-t446+t458)
      t585 = t3*t584
      t586 = t387+t585
      t589 = t586*t392*t525/2
      t590 = t419*s
      t591 = t176*t529
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf144 = -t589
      else
      spf144 = t590*t591
      end if
      if ((ABS(DBLE(spf144)).eq.0d0).or.(CIsBad(spf144))) spf144 = -t589
      spf144 = log(spf144)
      t593 = -t387-t585+t187
      t596 = t593*t392*t525/2
      t597 = t421*s
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf148 = t596
      else
      spf148 = -t597*t591
      end if
      if ((ABS(DBLE(spf148)).eq.0d0).or.(CIsBad(spf148))) spf148 = t596
      spf148 = log(spf148)
      t600 = t423*t1*t4
      t603 = t600*t11*t392*t489
      t604 = t430*s
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf150 = t603
      else
      spf150 = -t604*t580
      end if
      if ((ABS(DBLE(spf150)).eq.0d0).or.(CIsBad(spf150))) spf150 = t603
      spf150 = log(spf150)
      t615 = -spf103-spf102+spf46*spf63+spf129*spf63+spf129*spf123+spf48
     #*spf33-spf48*spf34+spf94*spf58-spf148*t42-spf148*spf129-spf46*spf1
     #19
      t623 = t42*spf64
      t628 = -spf46*spf42+spf46*spf26+spf46*spf97+spf142*spf114+spf142*t
     #42+spf46*spf64+t623+spf144*spf115-spf46*spf118+spf144*t42-spf115*s
     #pf123
      t631 = t42*spf56
      t632 = t42*spf55
      t642 = -t631-t632-spf107*spf42-spf107*spf101+spf107*spf86-spf114*s
     #pf122+spf128*spf64+spf128*spf122-spf150*spf128-spf150*t42-spf115*s
     #pf56
      t650 = t42*spf63
      t651 = 6*t650
      t663 = -6*spf114*spf55-6*spf94*spf99+6*spf94*spf101+t651+6*spf77+6
     #*spf87+6*spf88-6*spf74+6*spf51-t133-6*IM*spf46*PI-6*IM*spf48*PI
      vV87 = (6*t615+6*t628+6*t642+t663)*t11*t353/18
      t668 = s-t390+t391
      t669 = 1/t668
      t671 = 2*t389*t669
      t672 = -t248
      t673 = 1/t672
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf25 = t671
      else
      spf25 = -2*t396*t673
      end if
      if ((ABS(DBLE(spf25)).eq.0d0).or.(CIsBad(spf25))) spf25 = t671
      spf25 = log(spf25)
      t676 = s*t668
      t678 = t676*t400/2
      t679 = s*t672
      t681 = t679*t384/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf31 = -t678
      else
      spf31 = t681
      end if
      if ((ABS(DBLE(spf31)).eq.0d0).or.(CIsBad(spf31))) spf31 = -t678
      spf31 = log(spf31)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf32 = t678
      else
      spf32 = -t681
      end if
      if ((ABS(DBLE(spf32)).eq.0d0).or.(CIsBad(spf32))) spf32 = t678
      spf32 = log(spf32)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf41 = -t266
      else
      spf41 = t269
      end if
      if ((ABS(DBLE(spf41)).eq.0d0).or.(CIsBad(spf41))) spf41 = -t266
      spf41 = log(spf41)
      t682 = t405-t172-s+t171
      spf45 = log(DCMPLX(-2*t383/t682))
      t687 = t682*t384/2
      spf47 = log(DCMPLX(-t687))
      spf52 = DENSPENCE(DCMPLX(-t687),spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf57 = t270
      else
      spf57 = t271
      end if
      if ((ABS(DBLE(spf57)).eq.0d0).or.(CIsBad(spf57))) spf57 = t270
      spf57 = log(spf57)
      t689 = t434-t436-t439-t437+t441-t444-t446+t449+t451-t453+t456+t458
      t690 = 1/t689
      t691 = t11*t690
      t693 = 2*t388*t258*t691
      t695 = -t465+t467+t414+t14-t469+t471-t24
      t696 = 1/t695
      t698 = 2*t383*t261*t696
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf73 = -t693
      else
      spf73 = -t698
      end if
      if (CIsBad(spf73)) spf73 = -t693
      spf73 = DENSPENCE(spf73,spcZip)
      t699 = t465-t166-t14+t468
      t701 = t258*t699*t691
      t703 = t261*t682*t696
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf76 = t701
      else
      spf76 = t703
      end if
      if (CIsBad(spf76)) spf76 = t701
      spf76 = DENSPENCE(spf76,spcZip)
      t706 = 2*t13*t668*t690
      t707 = t672*t696
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf85 = t706
      else
      spf85 = 2*t484*t707
      end if
      if ((ABS(DBLE(spf85)).eq.0d0).or.(CIsBad(spf85))) spf85 = t706
      spf85 = log(spf85)
      t710 = -t500+t501+t504-t507+t508-t509-t511-t514-t516+t517+t519-t52
     #0+t522
      t711 = 1/t710
      t712 = Mzsq*t711
      t713 = t499*t712
      t714 = 2*t713
      t715 = t414-t24-t387+t466-t493+t494
      t716 = 1/t715
      t717 = t492*t716
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf89 = -t714
      else
      spf89 = -4*t717
      end if
      if (CIsBad(spf89)) spf89 = -t714
      spf89 = DENSPENCE(spf89,spcZip)
      t719 = t405-s+t390
      t720 = 1/t719
      t721 = t389*t720
      t722 = 2*t721
      t723 = t414-t24-t387+t466+t493-t494
      t724 = 1/t723
      t725 = t492*t724
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf90 = -t722
      else
      spf90 = -4*t725
      end if
      if (CIsBad(spf90)) spf90 = -t722
      spf90 = DENSPENCE(spf90,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf93 = -t693
      else
      spf93 = -t698
      end if
      if ((ABS(DBLE(spf93)).eq.0d0).or.(CIsBad(spf93))) spf93 = -t693
      spf93 = log(spf93)
      t728 = t689*t669*t264
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf96 = -t728
      else
      spf96 = t695*t11*t673*t267
      end if
      if ((ABS(DBLE(spf96)).eq.0d0).or.(CIsBad(spf96))) spf96 = -t728
      spf96 = log(spf96)
      t733 = t241*t668*t690
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf98 = -t733
      else
      spf98 = t245*s*t707
      end if
      if ((ABS(DBLE(spf98)).eq.0d0).or.(CIsBad(spf98))) spf98 = -t733
      spf98 = log(spf98)
      t737 = t668*t258*t690
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf100 = -t737
      else
      spf100 = t679*t261*t696
      end if
      if ((ABS(DBLE(spf100)).eq.0d0).or.(CIsBad(spf100))) spf100 = -t737
      spf100 = log(spf100)
      t741 = t4*t699*t712
      t742 = t4*t682
      t743 = t742*t716
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf104 = t741
      else
      spf104 = 2*t743
      end if
      if (CIsBad(spf104)) spf104 = t741
      spf104 = DENSPENCE(spf104,spcZip)
      t746 = t699*t11*t720
      t747 = t742*t724
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf105 = t746
      else
      spf105 = 2*t747
      end if
      if (CIsBad(spf105)) spf105 = t746
      spf105 = DENSPENCE(spf105,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf106 = t701
      else
      spf106 = t703
      end if
      if ((ABS(DBLE(spf106)).eq.0d0).or.(CIsBad(spf106))) spf106 = t701
      spf106 = log(spf106)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf112 = -t721
      else
      spf112 = -2*t725
      end if
      if ((ABS(DBLE(spf112)).eq.0d0).or.(CIsBad(spf112))) spf112 = -t721
      spf112 = log(spf112)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf113 = -t713
      else
      spf113 = -2*t717
      end if
      if ((ABS(DBLE(spf113)).eq.0d0).or.(CIsBad(spf113))) spf113 = -t713
      spf113 = log(spf113)
      t752 = 2*t719*t669
      t754 = t673*t35
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf116 = -t752
      else
      spf116 = t723*t11*t754
      end if
      if ((ABS(DBLE(spf116)).eq.0d0).or.(CIsBad(spf116))) spf116 = -t752
      spf116 = log(spf116)
      t760 = 2*t423*t710*t11*t669*t35
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf117 = -t760
      else
      spf117 = t715*t11*t754
      end if
      if ((ABS(DBLE(spf117)).eq.0d0).or.(CIsBad(spf117))) spf117 = -t760
      spf117 = log(spf117)
      t763 = t668*t720
      t764 = t763/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf120 = -t764
      else
      spf120 = t679*t4*t724
      end if
      if ((ABS(DBLE(spf120)).eq.0d0).or.(CIsBad(spf120))) spf120 = -t764
      spf120 = log(spf120)
      t769 = t676*t413*t711/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf121 = -t769
      else
      spf121 = t679*t4*t716
      end if
      if ((ABS(DBLE(spf121)).eq.0d0).or.(CIsBad(spf121))) spf121 = -t769
      spf121 = log(spf121)
      t772 = t746/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf126 = t772
      else
      spf126 = t747
      end if
      if ((ABS(DBLE(spf126)).eq.0d0).or.(CIsBad(spf126))) spf126 = t772
      spf126 = log(spf126)
      t773 = t741/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf127 = t773
      else
      spf127 = t743
      end if
      if ((ABS(DBLE(spf127)).eq.0d0).or.(CIsBad(spf127))) spf127 = t773
      spf127 = log(spf127)
      t774 = t672*t724
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf141 = -t763
      else
      spf141 = t579*t774
      end if
      if ((ABS(DBLE(spf141)).eq.0d0).or.(CIsBad(spf141))) spf141 = -t763
      spf141 = log(spf141)
      t778 = t586*t668*t712/2
      t779 = t672*t716
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf143 = -t778
      else
      spf143 = t590*t779
      end if
      if ((ABS(DBLE(spf143)).eq.0d0).or.(CIsBad(spf143))) spf143 = -t778
      spf143 = log(spf143)
      t783 = t593*t668*t712/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf147 = t783
      else
      spf147 = -t597*t779
      end if
      if ((ABS(DBLE(spf147)).eq.0d0).or.(CIsBad(spf147))) spf147 = t783
      spf147 = log(spf147)
      t787 = t600*t11*t668*t720
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf149 = t787
      else
      spf149 = -t604*t774
      end if
      if ((ABS(DBLE(spf149)).eq.0d0).or.(CIsBad(spf149))) spf149 = t787
      spf149 = log(spf149)
      t800 = spf147*t42+spf147*spf127-spf141*spf112-spf141*t42-spf47*spf
     #31-spf45*spf25+spf45*spf41+spf45*spf117-spf127*spf63-spf127*spf121
     #-spf93*spf57
      t813 = -spf143*spf113-spf143*t42+spf93*spf98+spf112*spf55-spf45*sp
     #f63+spf106*spf100-spf106*spf85+spf106*spf41+spf45*spf116-spf93*spf
     #100+spf112*spf120
      t825 = spf113*spf121-spf45*spf96-spf45*spf64+spf47*spf32+spf149*t4
     #2+spf149*spf126+spf113*spf56-spf126*spf120-spf126*spf64-t623+t631
      t841 = 6*t632-t651+6*IM*spf47*PI+6*spf105+6*spf104-6*spf76-6*spf89
     #-6*spf90+6*spf73-6*spf52+t133+6*IM*spf45*PI
      vV88 = -(6*t800+6*t813+6*t825+t841)*t11*t373/18
      spf7 = log(DCMPLX(t383))
      spf9 = log(DCMPLX(t54))
      spf12 = log(DCMPLX(t396))
      spf15 = log(DCMPLX(t54*t35))
      t848 = -Mzsq+s+Ieps
      t849 = t848*t11
      spf16 = log(DCMPLX(t849))
      t850 = Ieps*t14
      t851 = t21*Ieps
      t852 = 2*t851
      t853 = Ieps**2
      t855 = 2*t445
      t859 = t442-t850+t852+t502+t4*t853-t855+t387*Ieps-2*t413*Ieps
      t861 = -s+t5
      t862 = 1/t861
      spf75 = log(DCMPLX(t859*t11*t862))
      t864 = -t861
      t865 = sqrt(t864)
      t866 = Mzsq*t865
      t867 = t866*t3
      t869 = Ieps*t865*t3
      t870 = t3*s
      t871 = t865*t870
      t872 = -t387-t867-t518+t869+t871+t14
      t873 = 1/t872
      spf78 = log(DCMPLX(t873*t865*t870))
      t876 = -s+Mzsq+Ieps
      t877 = 1/t876
      t878 = t877*t35
      t880 = t871*t878/2
      t881 = t871+t387-t867+t518+t869-t14
      t882 = 1/t881
      t883 = t871*t882
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf79 = t880
      else
      spf79 = t883
      end if
      if ((ABS(DBLE(spf79)).eq.0d0).or.(CIsBad(spf79))) spf79 = t880
      spf79 = log(spf79)
      t884 = 1/t865
      t886 = 1/t870
      spf80 = log(DCMPLX(t872*t884*t886))
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf81 = -t880
      else
      spf81 = -t883
      end if
      if ((ABS(DBLE(spf81)).eq.0d0).or.(CIsBad(spf81))) spf81 = -t880
      spf81 = log(spf81)
      t888 = t876*t4
      t891 = 2*t888*t884*t886
      t892 = t881*t884
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf82 = -t891
      else
      spf82 = -t892*t886
      end if
      if ((ABS(DBLE(spf82)).eq.0d0).or.(CIsBad(spf82))) spf82 = -t891
      spf82 = log(spf82)
      t895 = t3*t873
      spf83 = log(DCMPLX(t383*t865*t895))
      spf84 = log(DCMPLX(t848*t865*t895))
      t899 = t865*t3
      t900 = t899*t877
      t901 = -s+t899
      t904 = t901*t865*t3*t882
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf91 = -t900
      else
      spf91 = t904
      end if
      if ((ABS(DBLE(spf91)).eq.0d0).or.(CIsBad(spf91))) spf91 = -t900
      spf91 = log(spf91)
      t906 = t876*t884*t60
      t908 = t60/t901
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf92 = -t906
      else
      spf92 = t892*t908
      end if
      if ((ABS(DBLE(spf92)).eq.0d0).or.(CIsBad(spf92))) spf92 = -t906
      spf92 = log(spf92)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf95 = t900
      else
      spf95 = -t904
      end if
      if ((ABS(DBLE(spf95)).eq.0d0).or.(CIsBad(spf95))) spf95 = t900
      spf95 = log(spf95)
      t910 = 2*s
      t911 = -t910+Mzsq+Ieps
      t912 = t911*t877
      t913 = t387-t22-t867+t518+t869
      t914 = t913*t882
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf108 = t912
      else
      spf108 = t914
      end if
      if ((ABS(DBLE(spf108)).eq.0d0).or.(CIsBad(spf108))) spf108 = t912
      spf108 = log(spf108)
      t915 = Mzsq*t870
      t916 = t866*s
      t917 = t23*t3
      t918 = t23*t865
      t919 = t865*t14
      t920 = t3*t14
      t921 = t915-t916+t917+t918-t919-t920
      t922 = t921*t60
      t923 = 1/t911
      t926 = t922*t923*t35/2
      t927 = t387-t867+t518+t869-t871-t14
      t928 = 1/t913
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf109 = t926
      else
      spf109 = t927*t928
      end if
      if ((ABS(DBLE(spf109)).eq.0d0).or.(CIsBad(spf109))) spf109 = t926
      spf109 = log(spf109)
      t930 = t876*t923
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf110 = t930
      else
      spf110 = t881*t928
      end if
      if ((ABS(DBLE(spf110)).eq.0d0).or.(CIsBad(spf110))) spf110 = t930
      spf110 = log(spf110)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf111 = t912
      else
      spf111 = t914
      end if
      if (CIsBad(spf111)) spf111 = t912
      spf111 = DENSPENCE(spf111,spcZip)
      t933 = t922*t878/2
      t934 = t927*t882
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf130 = t933
      else
      spf130 = t934
      end if
      if ((ABS(DBLE(spf130)).eq.0d0).or.(CIsBad(spf130))) spf130 = t933
      spf130 = log(spf130)
      t938 = 2*t888*t3/t921
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf131 = t938
      else
      spf131 = t881/t927
      end if
      if ((ABS(DBLE(spf131)).eq.0d0).or.(CIsBad(spf131))) spf131 = t938
      spf131 = log(spf131)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf132 = t933
      else
      spf132 = t934
      end if
      if (CIsBad(spf132)) spf132 = t933
      spf132 = DENSPENCE(spf132,spcZip)
      t941 = t927*t873
      spf133 = log(DCMPLX(-t941))
      t946 = 2*t888*t3/(-t915-t916-t917+t918+t919+t920)
      t947 = t881*t873
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf134 = -t946
      else
      spf134 = -t947
      end if
      if ((ABS(DBLE(spf134)).eq.0d0).or.(CIsBad(spf134))) spf134 = -t946
      spf134 = log(spf134)
      spf135 = DENSPENCE(DCMPLX(-t941),spcZip)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf136 = -t946
      else
      spf136 = -t947
      end if
      if (CIsBad(spf136)) spf136 = -t946
      spf136 = DENSPENCE(spf136,spcZip)
      t948 = -s+Mzsq
      t949 = 1/t948
      t950 = t899*t949
      t952 = t413*t899
      t954 = t23*t899
      t955 = t871*t4
      t957 = t445-t952-3*t851+t954+t955-t442-t871*Ieps+t850
      t958 = 1/t957
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf145 = -t950
      else
      spf145 = t899*t901*t4*t958
      end if
      if ((ABS(DBLE(spf145)).eq.0d0).or.(CIsBad(spf145))) spf145 = -t950
      spf145 = log(spf145)
      t962 = t948*t884*t60
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf146 = -t962
      else
      spf146 = t957*t884*t908*t35
      end if
      if ((ABS(DBLE(spf146)).eq.0d0).or.(CIsBad(spf146))) spf146 = -t962
      spf146 = log(spf146)
      t966 = t14**2
      t971 = Mzsq*t437
      t973 = Ieps*Mzsq
      t981 = 1/t865/t864
      t984 = t4*(t966-t437*Ieps+2*t23*t14+t500+t188*t853-2*t971+t973*t14
     #-2*t413*t518)/t920*t981/2
      t987 = 2*t4*t865*t3
      t990 = t981*t60
      t991 = t901**2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf151 = -t984
      else
      spf151 = t859*(t987-t22-t871+t14)*t990/t991
      end if
      if ((ABS(DBLE(spf151)).eq.0d0).or.(CIsBad(spf151))) spf151 = -t984
      spf151 = log(spf151)
      t997 = s+t899
      t998 = t997**2
      spf152 = log(DCMPLX(t859*(t987+t22-t871-t14)*t990/t998))
      t1004 = t1*t865*t60*t949
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf153 = -t1004
      else
      spf153 = t54*t865*t3*t901*t958
      end if
      if ((ABS(DBLE(spf153)).eq.0d0).or.(CIsBad(spf153))) spf153 = -t100
     #4
      spf153 = log(spf153)
      t1009 = t911*t949
      t1010 = t913*t4
      t1011 = t1010*t958
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf154 = t1009
      else
      spf154 = t1011
      end if
      if ((ABS(DBLE(spf154)).eq.0d0).or.(CIsBad(spf154))) spf154 = t1009
      spf154 = log(spf154)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf155 = t1009
      else
      spf155 = t1011
      end if
      if (CIsBad(spf155)) spf155 = t1009
      spf155 = DENSPENCE(spf155,spcZip)
      t1012 = t876*t949
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf156 = t1012
      else
      spf156 = t881*t4*t958
      end if
      if (CIsBad(spf156)) spf156 = t1012
      spf156 = DENSPENCE(spf156,spcZip)
      t1016 = t35*Mzsq
      t1017 = t412+t973-t518+t14-t465
      t1018 = 1/t1017
      t1020 = t899*t997*t1016*t1018
      t1021 = t1020/4
      t1022 = t997*t4
      t1023 = 2*t952
      t1024 = 2*t954
      t1025 = t871*Mzsq
      t1026 = t416*t14
      t1029 = 4*t4*t416*s
      t1030 = t871*t416
      t1031 = 2*t442
      t1032 = 2*t955
      t1033 = -t855-t1023+t852+t1024+t433+t1025-t1026+t1029-t1030-t1031-
     #t1032
      t1034 = 1/t1033
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf157 = t1021
      else
      spf157 = t899*t1022*t1034
      end if
      if ((ABS(DBLE(spf157)).eq.0d0).or.(CIsBad(spf157))) spf157 = t1021
      spf157 = log(spf157)
      t1037 = -t855-t1023+t852+t1024+t433+t1025+t1026-t1029+t1030-t1031-
     #t1032
      t1038 = 1/t1037
      spf158 = log(DCMPLX(t899*t1022*t1038))
      t1044 = 1/t997
      t1047 = 4*t423*t1017*t4*t884*t60*t1044
      t1050 = t60*t1044*t35
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf159 = t1047
      else
      spf159 = t1033*t884*t1050
      end if
      if ((ABS(DBLE(spf159)).eq.0d0).or.(CIsBad(spf159))) spf159 = t1047
      spf159 = log(spf159)
      spf160 = log(DCMPLX(t1037*t884*t1050))
      t1055 = t911*Mzsq*t1018
      t1056 = t1055/2
      t1057 = t1010*t1034
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf161 = t1056
      else
      spf161 = t1057
      end if
      if ((ABS(DBLE(spf161)).eq.0d0).or.(CIsBad(spf161))) spf161 = t1056
      spf161 = log(spf161)
      t1078 = -2*t413*t870-2*t413*t865*s+2*t457*t3+2*t457*t865+t920*Mzsq
     #+t919*Mzsq+t14*t584-4*t21*t584+t871*t584-2*t4*t920-2*t919*t4
      t1081 = t911*t188*t3/t1078
      t1082 = 2*t1081
      t1083 = t1010*t1038
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf162 = t1082
      else
      spf162 = t1083
      end if
      if ((ABS(DBLE(spf162)).eq.0d0).or.(CIsBad(spf162))) spf162 = t1082
      spf162 = log(spf162)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf163 = t1055
      else
      spf163 = 2*t1057
      end if
      if (CIsBad(spf163)) spf163 = t1055
      spf163 = DENSPENCE(spf163,spcZip)
      t1085 = 4*t1081
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf164 = t1085
      else
      spf164 = 2*t1083
      end if
      if (CIsBad(spf164)) spf164 = t1085
      spf164 = DENSPENCE(spf164,spcZip)
      t1087 = t1020/2
      t1088 = t997*t865
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf165 = t1087
      else
      spf165 = t1088*t3*t417*t1034
      end if
      if ((ABS(DBLE(spf165)).eq.0d0).or.(CIsBad(spf165))) spf165 = t1087
      spf165 = log(spf165)
      spf166 = log(DCMPLX(t1088*t3*t419*t1038))
      t1098 = t921*Mzsq*t60*t35*t1018
      t1099 = t1098/4
      t1100 = t4*t927
      t1101 = t1100*t1034
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf167 = t1099
      else
      spf167 = t1101
      end if
      if ((ABS(DBLE(spf167)).eq.0d0).or.(CIsBad(spf167))) spf167 = t1099
      spf167 = log(spf167)
      t1102 = t1100*t1038
      spf168 = log(DCMPLX(t1102))
      t1103 = t1098/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf169 = t1103
      else
      spf169 = 2*t1101
      end if
      if (CIsBad(spf169)) spf169 = t1103
      spf169 = DENSPENCE(spf169,spcZip)
      spf170 = DENSPENCE(DCMPLX(2*t1102),spcZip)
      spf171 = log(DCMPLX(-t421*t997*t899*t1038))
      t1112 = t1*t60*t1088*t1018/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf172 = -t1112
      else
      spf172 = -t430*t997*t899*t1034
      end if
      if ((ABS(DBLE(spf172)).eq.0d0).or.(CIsBad(spf172))) spf172 = -t111
     #2
      spf172 = log(spf172)
      t1133 = spf109*spf159-spf161*spf64+spf110*spf151-spf161*spf157+spf
     #168*spf158-spf166*spf168+spf168*spf56-spf165*spf167+spf167*spf55+s
     #pf167*spf157-spf162*spf63-spf162*spf158-spf110*spf146+spf131*spf80
     #-2*spf131*t42+spf131*spf7
      t1151 = -spf131*spf12+spf131*spf82-spf131*spf75+spf108*spf91-spf10
     #8*spf95+spf110*t42+spf109*spf160+spf110*spf15-spf110*spf92-spf110*
     #spf9+spf154*spf145-spf154*spf153+spf154*spf15-spf130*spf81+spf130*
     #spf79-spf109*t42-spf109*spf64
      t1167 = spf109*spf9-spf109*spf63+spf134*spf16+spf134*spf78-spf134*
     #spf84+spf171*t42+spf171*spf162-spf133*spf78+spf133*spf83-spf133*sp
     #f12-spf109*spf152-spf165*t42+spf172*t42+spf172*spf161-t623+t631
      t1177 = t632-t650-spf166*t42-spf170-spf169+IM*spf110*PI-IM*spf130*
     #PI+IM*spf108*PI-IM*spf131*PI+spf164+spf163+spf135-spf136-spf155+sp
     #f156+spf132-spf111
      t1182 = sqrt(-t861*t11)
      vV89 = -(t1133+t1151+t1167+t1177)*t11/t1182/4
      t1186 = spf13**2
      vV90 = Mzsq*(t133+3*t1186+6*spf13+6)/36
      t1191 = 1/t848
      spf21 = log(DCMPLX(-Mzsq*t1191))
      spf23 = DENSPENCE(DCMPLX(t320*t1191),spcZip)
      t1194 = t167*t11
      t1196 = t1194*t423/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf36 = -t1196
      else
      spf36 = t173*t423/2
      end if
      if ((ABS(DBLE(spf36)).eq.0d0).or.(CIsBad(spf36))) spf36 = -t1196
      spf36 = log(spf36)
      t1199 = Ml*s
      t1201 = 2*t1199*t192
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf44 = -t1201
      else
      spf44 = 2*Ml*t195
      end if
      if ((ABS(DBLE(spf44)).eq.0d0).or.(CIsBad(spf44))) spf44 = -t1201
      spf44 = log(spf44)
      t1204 = Mzsq*DBLE(Mzsq)
      t1205 = sqrt(Mzsq)
      t1206 = t1205*Ml
      t1208 = 2*t1206*DBLE(Mzsq)
      t1213 = sqrt((t23-t1204-t1208)/(t23-t1204+t1208))
      t1214 = t1213+1
      t1216 = 2/t1214
      t1217 = 1/Ml
      t1218 = t1205*t1217
      t1219 = 2*t1206
      t1224 = sqrt((Ieps-Mzsq-t1219)/(Ieps-Mzsq+t1219))
      t1225 = t1224-1
      t1226 = t1224+1
      t1227 = 1/t1226
      if (ABS(DBLE(Ml**2/DBLE(Mzsq))).lt.epsMl) then
      spf124 = t1216
      else
      spf124 = t1218*t1225*t1227
      end if
      if ((ABS(DBLE(spf124)).eq.0d0).or.(CIsBad(spf124))) spf124 = t1216
      spf124 = log(spf124)
      t1232 = Mzsq*t1214*t35/2
      t1233 = 1/t1225
      if (ABS(DBLE(Ml**2/DBLE(Mzsq))).lt.epsMl) then
      spf125 = t1232
      else
      spf125 = t1218*t1226*t1233
      end if
      if ((ABS(DBLE(spf125)).eq.0d0).or.(CIsBad(spf125))) spf125 = t1232
      spf125 = log(spf125)
      t1237 = 2*t1206*s
      t1242 = sqrt((t23-t387-t1237)/(t23-t387+t1237))
      t1243 = t1242+1
      t1244 = 1/t1243
      t1245 = t423*t1244
      t1246 = t1194*t1245
      t1248 = 1/t1205
      t1249 = t1248*t1217
      t1250 = t1249*t1227
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf139 = -t1246
      else
      spf139 = t173*t1225*t1250/2
      end if
      if ((ABS(DBLE(spf139)).eq.0d0).or.(CIsBad(spf139))) spf139 = -t124
     #6
      spf139 = log(spf139)
      t1253 = t1243*t35
      t1255 = t1194*t1253/4
      t1257 = t1249*t1233
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf140 = -t1255
      else
      spf140 = t173*t1226*t1257/2
      end if
      if ((ABS(DBLE(spf140)).eq.0d0).or.(CIsBad(spf140))) spf140 = -t125
     #5
      spf140 = log(spf140)
      t1261 = t488*t423*t1244
      t1263 = 2*t1206*t1224
      t1264 = s*t1224
      t1265 = t171*t1224
      t1267 = 2*Ieps*t1224
      t1270 = t1217*t1227
      t1272 = (-t1263-t1219+t1264-s+t1265-t171-t1267+t172)*t1248*t1270/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf173 = -t1261
      else
      spf173 = -t1272
      end if
      if (CIsBad(spf173)) spf173 = -t1261
      spf173 = DENSPENCE(spf173,spcZip)
      t1273 = t1205*t412
      t1274 = 1/t1273
      t1275 = t4*Ml
      t1276 = t1275*Mzsq
      t1278 = 2*t1276*t162
      t1279 = t1276*t518
      t1280 = t14*Ml
      t1281 = t1280*t412
      t1284 = 2*t1273*t188*CsTh
      t1285 = t1205*Mzsq
      t1286 = t1285*t188
      t1288 = 4*t1286*s
      t1290 = 4*t1286*Ieps
      t1291 = t1205*t188
      t1293 = 2*t1291*t518
      t1294 = t4*t1285
      t1295 = t1294*t14
      t1296 = t4*t1273
      t1298 = 2*t1296*Ieps
      t1300 = 2*t1296*s
      t1301 = t1275*t412
      t1303 = 2*t1301*s
      t1305 = 2*t1301*Ieps
      t1306 = t1278+t1279+t1281-t1284-t1288-t1290+t1293+t1295-t1298-t130
     #0-t1303-t1305
      t1308 = 3*t1291*t162
      t1309 = t162*t1273
      t1310 = t1280*t505
      t1311 = t390*Ieps
      t1312 = t1276*t1311
      t1313 = t1294*t162
      t1315 = 2*t1296*t390
      t1317 = 2*t1301*t390
      t1319 = 2*t1291*t1311
      t1321 = 2*t1286*t390
      t1323 = 2*t1276*t14
      t1325 = 3*t1291*t14
      t1326 = t1273*t14
      t1327 = t1308+t1309+t1310+t1312+t1313-t1315-t1317+t1319-t1321+t132
     #3+t1325+t1326
      t1330 = t35*t11
      t1332 = t1274*(t1306+t1327)*t1330/2
      t1335 = t1217*t1233
      t1337 = (-t1263+t1219+t1264+s+t1265+t171-t1267-t172)*t1248*t1335/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf174 = -t1332
      else
      spf174 = -t1337
      end if
      if (CIsBad(spf174)) spf174 = -t1332
      spf174 = DENSPENCE(spf174,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf178 = -t1332
      else
      spf178 = -t1337
      end if
      if ((ABS(DBLE(spf178)).eq.0d0).or.(CIsBad(spf178))) spf178 = -t133
     #2
      spf178 = log(spf178)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf180 = -t1261
      else
      spf180 = -t1272
      end if
      if ((ABS(DBLE(spf180)).eq.0d0).or.(CIsBad(spf180))) spf180 = -t126
     #1
      spf180 = log(spf180)
      t1340 = spf44**2
      t1342 = 24*spf23
      t1357 = t133-24*spf44*spf21-12*t1340+t1342-12*spf173-12*spf180*spf
     #139+12*spf180*spf36+12*spf180*spf124-12*spf174-12*spf178*spf140+12
     #*spf178*spf36+12*spf178*spf125
      t1359 = t11*t353
      vV91 = -t1357*t949*t1359/144
      vV92 = -(spf44+spf21)*t949*t1359/6
      t1366 = t241*t11
      t1368 = t1366*t423/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf35 = t1368
      else
      spf35 = -t245*t423/2
      end if
      if ((ABS(DBLE(spf35)).eq.0d0).or.(CIsBad(spf35))) spf35 = t1368
      spf35 = log(spf35)
      t1372 = 2*t1199*t259
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf43 = t1372
      else
      spf43 = -2*Ml*t262
      end if
      if ((ABS(DBLE(spf43)).eq.0d0).or.(CIsBad(spf43))) spf43 = t1372
      spf43 = log(spf43)
      t1375 = t1366*t1245
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf137 = t1375
      else
      spf137 = -t245*t1225*t1250/2
      end if
      if ((ABS(DBLE(spf137)).eq.0d0).or.(CIsBad(spf137))) spf137 = t1375
      spf137 = log(spf137)
      t1380 = t1366*t1253/4
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf138 = t1380
      else
      spf138 = -t245*t1226*t1257/2
      end if
      if ((ABS(DBLE(spf138)).eq.0d0).or.(CIsBad(spf138))) spf138 = t1380
      spf138 = log(spf138)
      t1384 = -t1323-t1325-t1326+t1313-t1315-t1317+t1319-t1321+t1310+t13
     #12-t1279+t1278
      t1385 = -t1281+t1309-t1295+t1298+t1300+t1303+t1305+t1308-t1284+t12
     #90+t1288-t1293
      t1389 = t1274*(t1384+t1385)*t1330/2
      t1393 = (t1263-t1219-t1264-s+t1265+t171+t1267+t172)*t1248*t1335/2
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf175 = t1389
      else
      spf175 = t1393
      end if
      if (CIsBad(spf175)) spf175 = t1389
      spf175 = DENSPENCE(spf175,spcZip)
      t1395 = t719*t423*t1244
      t1399 = (t1263+t1219-t1264+s+t1265-t171+t1267-t172)*t1248*t1270/2
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf176 = t1395
      else
      spf176 = t1399
      end if
      if (CIsBad(spf176)) spf176 = t1395
      spf176 = DENSPENCE(spf176,spcZip)
      if (ABS(DBLE(Ml**2/s)**2).lt.epsMl) then
      spf177 = t1389
      else
      spf177 = t1393
      end if
      if ((ABS(DBLE(spf177)).eq.0d0).or.(CIsBad(spf177))) spf177 = t1389
      spf177 = log(spf177)
      if (ABS(DBLE(Ml**2/s)).lt.epsMl) then
      spf179 = t1395
      else
      spf179 = t1399
      end if
      if ((ABS(DBLE(spf179)).eq.0d0).or.(CIsBad(spf179))) spf179 = t1395
      spf179 = log(spf179)
      t1402 = spf43**2
      t1418 = -t133+24*spf43*spf21+12*t1402-t1342+12*spf176+12*spf179*sp
     #f137-12*spf179*spf35-12*spf179*spf124+12*spf175+12*spf177*spf138-1
     #2*spf177*spf35-12*spf177*spf125
      t1420 = t11*t373
      vV93 = t1418*t949*t1420/144
      vV94 = (spf43+spf21)*t949*t1420/6
      vV95 = vV79
      vV96 = vV80
      spf71 = log(DCMPLX(-t421/t419))
      if (ABS(DBLE(Ml**2/DBLE(Mzsq))**2).lt.epsMl) then
      spf72 = -t428
      else
      spf72 = -t430/t417
      end if
      if ((ABS(DBLE(spf72)).eq.0d0).or.(CIsBad(spf72))) spf72 = -t428
      spf72 = log(spf72)
      vV97 = (t5+spf72*Mzsq+2*t4*spf14+spf71*Mzsq+spf71*t416-spf72*t416-
     #2*t4*spf64+4*t4*t42-2*t4*spf63)*t35/48
      spf24 = log(DCMPLX(t848*t384))
      t1444 = spf24*Mzsq
      t1445 = spf24*s
      t1446 = log(DCMPLX(t423))
      t1447 = t1446*s
      vV98 = (t910+t1444-t1445+t1447)*t11/12
      spf17 = log(DCMPLX(-t396))
      spf18 = DENSPENCE(DCMPLX(s*t1191),spcZip)
      spf20 = log(DCMPLX(-t849))
      t1458 = spf17**2
      t1461 = spf20**2
      t1475 = t1446*spf24
      t1481 = t1446**2
      t1484 = 24*s+6*Mzsq*spf18+12*t1444+6*Mzsq*spf17*spf24+3*Mzsq*t1458
     #-3*Mzsq*t1461-6*s*spf18-12*t1445-6*s*spf17*spf24-3*s*t1458+3*s*t14
     #61+12*t1447+6*t1475*Mzsq-6*t1475*s+t133*s+3*t1481*s
      vV99 = t1484*t11/144
      vV100 = vV83
      vV101 = vV84
      vV102 = vV85
      vV103 = vV86
      vV104 = vV87
      vV105 = vV88
      vV106 = vV89
      vV107 = vV91
      vV108 = vV92
      vV109 = vV93
      vV110 = vV94
      vV111 = vV95
      vV112 = vV96
      vV113 = vV97
      vV114 = vV98
      vV115 = vV77
      vV116 = vV78
      vV117 = vV111
      vV118 = vV112
      vV119 = vV115
      vV120 = vV116
      vV121 = vV117
      vV122 = vV118
      vV123 = vV82
      vV124 = vV100
      vV125 = vV101
      vV126 = vV102
      vV127 = vV103
      vV128 = vV104
      vV129 = vV105
      vV130 = vV63
      vV131 = vV106
      vV132 = vV90
      vV133 = vV64
      vV134 = vV65
      vV135 = vV107
      vV136 = vV108
      vV137 = vV66
      vV138 = vV67
      vV139 = vV109
      vV140 = vV110
      vV141 = vV121
      vV142 = vV122
      vV143 = vV113
      vV144 = vV114
      vV145 = vV99
      vV146 = vV123
      vV147 = vV124
      vV148 = vV125
      vV149 = vV126
      vV150 = vV127
      vV151 = vV128
      vV152 = vV129
      vV153 = vV130
      vV154 = vV131
      vV155 = vV132
      vV156 = vV133
      vV157 = vV134
      vV158 = vV135
      vV159 = vV136
      vV160 = vV137
      vV161 = vV138
      vV162 = vV139
      vV163 = vV140
      vV164 = vV141
      vV165 = vV142
      vV166 = vV143
      vV167 = vV81
      vV168 = vV144
      vV169 = vV145
      vV170 = vV167
      vV171 = vV170
      vV172 = vV119
      vV173 = vV120
      vV174 = vV164
      vV175 = vV165
      vV176 = vV171
      vV177 = vV146
      vV178 = vV147
      vV179 = vV148
      vV180 = vV149
      vV181 = vV150
      vV182 = vV151
      vV183 = vV152
      vV184 = vV154
      vV185 = vV155
      vV186 = vV158
      vV187 = vV159
      vV188 = vV162
      vV189 = vV163
      vV190 = vV174
      vV191 = vV175
      vV192 = vV166
      vV193 = vV168
      vV194 = vV169
      t1486 = 1/t133
      t1487 = vV55*t1486
      t1488 = t862*t4
      t1491 = t862*s
      vI79pc1 = t1487*t1488/7+5.D0/28.D0*t1487*t1491+t1487*t1491*t470/28
      t1497 = vV56*t862
      t1498 = t1486*t4
      t1501 = s*vV56
      t1502 = t862*t1486
      t1505 = t1502*t470
      t1509 = t1502*t188
      vI80pc1 = t1497*t1498/8+t1501*t1502/16+t1501*t1505/16-t11*vV56*t15
     #09/2-t1497*t1498*t470/8
      t1516 = t1486*t470
      t1520 = t1488*t1486
      t1523 = s*vV57
      t1525 = t35*t1486
      t1526 = t1525*t470
      t1529 = t862*t35
      t1530 = t1529*t1486
      vI81pc1 = vV57*t862*t1516/10+2.D0/5.D0*t11*vV57*t1520-t1523*t862*t
     #1526/40+t1523*t1530/40
      t1533 = vV58*t862
      t1536 = s*vV58
      vI82pc1 = t1533*t1486/8+t1536*t1530/16+t1536*t862*t1526/16-t11*vV5
     #8*t1520/2-t1533*t1516/8
      spf5 = log(DCMPLX(t4))
      t1558 = t1486*t188
      t1562 = t1486*spf5
      t1568 = t4*spf5
      vI83pc1 = -t1488*t1516/16+t1520/16+t1491*t1516/96-t188*t11*t1502/4
     #+t1491*t1486/96+t11*t862*t1558*spf5/4-t1491*t1562*t470/64-t1491*t1
     #562/64+t1502*t1568*t470/16
      t1572 = vV59*t353
      t1573 = t1486*t15
      t1574 = t1573*t188
      t1577 = t1486*t11
      t1578 = t1577*t4
      t1582 = t11*t4
      t1583 = t1582*BCsTh
      t1586 = t1486*BCsTh
      vI84pc1 = 4.D0/3.D0*t1572*t1574-7.D0/6.D0*t1572*t1578-4.D0/3.D0*t1
     #572*t1486*t1583+t1572*t1586/4+t1572*t1516/4
      t1591 = vV60*t353
      vI85pc1 = -t1591*t1578+t1591*t1586/2-2*t1591*t1486*t1583+t1591*t15
     #16/2+2*t1591*t1574
      t1602 = vV61*t673
      vI86pc1 = -4.D0/3.D0*t1602*t1574+7.D0/6.D0*t1602*t1578-4.D0/3.D0*t
     #1602*t1486*t1583-t1602*t1516/4+t1602*t1586/4
      t1614 = vV62*t673
      vI87pc1 = t1614*t1578+t1614*t1586/2-2*t1614*t1486*t1583-t1614*t151
     #6/2-2*t1614*t1574
      t1625 = BCsTh*vV63
      t1626 = t1502*t14
      t1630 = t1498*s
      vI88pc1 = t1625*t1626/4-2*t1625*t862*t1630+2*t1625*t1509
      t1635 = s*vV64
      t1638 = t14*vV64
      t1645 = t1586*t4
      t1648 = t470*BCsTh
      t1649 = t1486*t1648
      vI89pc1 = -t1635*t1498/16-t1638*t1486/32-t1638*t1586/16-t1638*t151
     #6/16-t1635*t1645/16-t1638*t1649/32
      t1652 = t14*vV65
      vI90pc1 = -t1652*t1486/6-t1652*t1586/6
      t1656 = s*vV66
      t1659 = t14*vV66
      vI91pc1 = -t1656*t1498/16-t1659*t1486/32+t1659*t1586/16-t1659*t151
     #6/16+t1656*t1645/16+t1659*t1649/32
      t1670 = t14*vV67
      vI92pc1 = -t1670*t1486/6+t1670*t1586/6
      t1674 = BCsTh*vV68
      t1677 = 1/(-1+t470)
      t1678 = t862*t1677
      t1686 = t1502*t1677
      t1691 = t1525*t1678
      vI93pc1 = -8*t1674*t15*t1558*t1678+4*t1674*t11*t1498*t1678-t1674*t
     #1686/2+t1648*vV68*s*t1691/4-t1674*s*t1691/4
      t1697 = BCsTh*vV69
      t1698 = t1525*t1677
      t1708 = t1498*t1677
      vI94pc1 = -t1697*t1698/3+4.D0/3.D0*t1697*t1577*t1677+t1648*vV69*t1
     #698/3+8.D0/3.D0*t1697*t15*t1708
      vI95pc1 = -BCsTh*vV70*t1486
      spf8 = spf11
      t1713 = t171*t1486
      t1714 = t1678*spf8
      t1717 = t11*BCsTh
      t1722 = t1678*t4
      t1726 = t1677*spf8
      t1727 = t1726*t4
      t1734 = s*t1648
      vI96pc1 = t1713*t1714/2+4*t1717*t1486*t1678*t188+10.D0/3.D0*t1649*
     #t1722+4.D0/3.D0*t1649*t862*t1727-4.D0/3.D0*t1586*t862*t1727+t171*t
     #1686-t1734*t1686-13.D0/3.D0*t1586*t1722-t1734*t1486*t1714/2
      t1744 = vV71*t1486
      vI97pc1 = t11*vV71*t1498/3+t1744/12+t1744*t470/12
      t1748 = t1486*spf8
      t1751 = t4*spf8
      vI98pc1 = 5.D0/28.D0*t1748+3.D0/14.D0*t1486+t1577*t1751/7+t1748*t4
     #70/28
      t1757 = abs(DCMPLX(Mzsq))
      t1758 = t1757**2
      t1760 = 1/(t387-t1758)
      t1761 = t1760*t949
      t1762 = t1761*t1498
      t1765 = t433*vV72
      t1766 = t1761*t1486
      t1769 = t1761*t1516
      vI99pc1 = 3.D0/2.D0*t387*vV72*t1762-5.D0/8.D0*t1765*t1766-t1765*t1
     #769/8
      t1773 = t1486*t949
      t1774 = t1760*t862
      t1775 = t1773*t1774
      t1778 = vV73*t1486
      t1779 = t433*t1778
      t1780 = t1761*t1488
      t1784 = t862*t470
      t1785 = t1761*t1784
      t1789 = t862*t188
      t1790 = t1761*t1789
      t1794 = t1761*t1488*t470
      vI100pc1 = -t971*vV73*t1775/4+3.D0/2.D0*t1779*t1780-t971*t1778*t17
     #85/4-2*t387*t1778*t1790+t1779*t1794/2
      t1797 = vV74*t1486
      t1798 = t971*t1797
      t1799 = t1761*t1529
      t1806 = t1761*t1529*t470
      vI101pc1 = -t1798*t1799/4+3.D0/2.D0*t433*vV74*t1775+t1798*t1806/4-
     #2*t387*t1797*t1780-t433*t1797*t1785/2
      t1815 = vV75*t1486
      t1816 = t971*t1815
      vI102pc1 = -t1816*t1799/4+3.D0/2.D0*t433*vV75*t1775-t1816*t1806/4-
     #2*t387*t1815*t1780+t433*t1815*t1785/2
      t1833 = 1/(-t433+s*t1758+t582-t1758*Mzsq)
      t1834 = t442*t1833
      t1835 = Mzsq*t862
      t1837 = t1834*t1835*t1516
      t1839 = t437*t1833
      t1840 = t1839*Mzsq
      t1843 = t1840*t1502*spf5/64
      t1845 = t1835*t1486
      t1846 = t515*t1833*t1845
      t1848 = t1834*t1845
      t1851 = t1840*t1505/96
      t1857 = s*t1833*Mzsq*t1502*t188*spf5/4
      t1859 = t1839*t1845/96
      t1860 = spf5*t470
      t1863 = t1840*t1502*t1860/64
      t1864 = t14*t1833
      t1866 = t1502*t1568
      t1872 = t1864*t1835*t1498*t1860/16
      vI103pc1 = t1837/12+t1843-t1846/6+t1848/12-t1851+t1857-t1859+t1863
     #-t1864*Mzsq*t1866/8-t1872
      t1873 = t433*vV76
      vI104pc1 = t1873*t1766/4-t387*vV76*t1762+t1873*t1769/4
      t1880 = t433*t1760
      t1882 = t1880*t1773*spf8
      t1884 = t433*t1766
      t1886 = t387*t1760
      t1890 = t1886*t1773*t1751
      t1892 = spf8*t470
      t1894 = t1880*t1773*t1892
      vI105pc1 = 5.D0/24.D0*t1882+t1884/4-t1886*t1773*t4-5.D0/6.D0*t1890
     #+t1894/24
      t1896 = vV77*t1486
      vI106pc1 = t433*t1896*t1780/19+5.D0/76.D0*t971*vV77*t1775+t971*t18
     #96*t1785/76
      t1906 = vV78*t1486
      t1907 = t433*t1906
      vI107pc1 = -t1907*t1780/8-t971*vV78*t1775/16-t971*t1906*t1785/16+t
     #387*t1906*t1790/2+t1907*t1794/8
      t1921 = vV79*t1486
      t1928 = t971*t1921
      vI108pc1 = -t433*t1921*t1785/10-2.D0/5.D0*t387*t1921*t1780+t1928*t
     #1806/40-t1928*t1799/40
      t1936 = vV80*t1486
      t1937 = t971*t1936
      vI109pc1 = -t433*vV80*t1775/8-t1937*t1799/16-t1937*t1806/16+t387*t
     #1936*t1780/2+t433*t1936*t1785/8
      vI110pc1 = t1837/16-t1848/16-t1851+t1846/4-t1859-t1857+t1863+t1843
     #-t1872
      t1954 = t433*vV81
      vI111pc1 = -t387*vV81*t1762/3-t1954*t1766/12-t1954*t1769/12
      vI112pc1 = -3.D0/14.D0*t1884-5.D0/28.D0*t1882-t1890/7-t1894/28
      vI113pc1 = t468*vV82*t1766/2
      t1966 = Mzsq*vV83*t1760
      t1968 = t1577*t353*t188
      t1971 = t1486*t353
      t1972 = t1971*t4
      t1975 = t4*BCsTh
      t1976 = t1971*t1975
      t1979 = t412*vV83
      t1980 = t1760*t1486
      t1981 = t1980*t353
      t1984 = t1979*t1760
      t1985 = t1971*BCsTh
      t1988 = t353*t4
      t1989 = t1577*t1988
      t1992 = t1486*s
      t1993 = t1992*t353
      t1996 = t353*t470
      t1997 = t1992*t1996
      t2000 = t353*BCsTh
      t2001 = t1992*t2000
      vI114pc1 = 16.D0/9.D0*t1966*t1968-5.D0/9.D0*t1966*t1972-8.D0/9.D0*
     #t1966*t1976+7.D0/18.D0*t1979*t1981+7.D0/18.D0*t1984*t1985-7.D0/9.D
     #0*t1984*t1989-t1966*t1993/18-t1966*t1997/9-t1966*t2001/6
      t2005 = Mzsq*vV84*t1760
      t2009 = t412*vV84
      t2014 = t2009*t1760
      vI115pc1 = -t2005*t1972-t2005*t2001/2+t2009*t1981/2-t2005*t1997/2+
     #t2014*t1985/2+2*t2005*t1968-t2014*t1989
      t2021 = Mzsq*vV85*t1760
      t2023 = t1577*t673*t188
      t2026 = t1486*t673
      t2027 = t2026*t4
      t2030 = t2026*t1975
      t2033 = t1992*t673
      t2036 = t673*BCsTh
      t2037 = t1992*t2036
      t2040 = t673*t470
      t2041 = t1992*t2040
      t2044 = t412*vV85
      t2045 = t1980*t673
      t2048 = t2044*t1760
      t2049 = t2026*BCsTh
      t2052 = t673*t4
      t2053 = t1577*t2052
      vI116pc1 = -16.D0/9.D0*t2021*t2023+5.D0/9.D0*t2021*t2027-8.D0/9.D0
     #*t2021*t2030+t2021*t2033/18-t2021*t2037/6+t2021*t2041/9-7.D0/18.D0
     #*t2044*t2045+7.D0/18.D0*t2048*t2049+7.D0/9.D0*t2048*t2053
      t2057 = Mzsq*vV86*t1760
      t2061 = t412*vV86
      t2066 = t2061*t1760
      vI117pc1 = t2057*t2027-t2057*t2037/2-t2061*t2045/2+t2057*t2041/2+t
     #2066*t2049/2-2*t2057*t2023+t2066*t2053
      t2072 = Mzsq*vV87
      t2073 = t2072*t1980
      t2074 = s*t470
      t2075 = t1988*t2074
      t2077 = t1988*t171
      t2081 = t412*vV87*t1760
      t2086 = t2072*t1760
      t2087 = t1971*t21
      t2092 = t1971*t468
      t2095 = t1971*t471
      t2098 = t14*t1648
      t2099 = t1971*t2098
      vI118pc1 = t2073*t2075+2.D0/3.D0*t2073*t2077+t2081*t2001/3+t2081*t
     #1997/6-t2086*t2087/3+t2081*t1993/6-t2086*t2092/6-t2086*t2095/3-t20
     #86*t2099/6-4.D0/3.D0*t2081*t1976
      t2104 = Mzsq*vV88
      t2105 = t2104*t1980
      t2106 = t2052*t2074
      t2108 = t2052*t171
      t2112 = t412*vV88*t1760
      t2117 = t2104*t1760
      t2118 = t2026*t21
      t2123 = t2026*t468
      t2126 = t2026*t471
      t2129 = t2026*t2098
      vI119pc1 = -t2105*t2106+2.D0/3.D0*t2105*t2108+t2112*t2037/3-t2112*
     #t2041/6+t2117*t2118/3-t2112*t2033/6-t2117*t2123/6+t2117*t2126/3-t2
     #117*t2129/6-4.D0/3.D0*t2112*t2030
      t2134 = t433*BCsTh
      t2135 = vV89*t1486
      vI120pc1 = -t2134*t2135*t1760+2*t466*t2135*t1760*t4
      vI121pc1 = -t468*vV90*t1766/2
      t2145 = t433*vV91
      t2146 = t1980*t4
      t2149 = t582*vV91
      t2153 = vV91*t1486*t1760
      t2156 = t971*vV91
      t2157 = t1980*BCsTh
      t2160 = t1980*t470
      t2165 = t500*vV91
      t2168 = t412*Mzsq
      t2169 = t2168*s
      t2172 = t1980*t1975
      t2177 = t1980*t1648
      vI122pc1 = -t2145*t2146/7-t2149*t2146/7+t971*t2153/14+t2156*t2157/
     #7+t2156*t2160/7-t500*t2153/14-t2165*t2157/7+t2169*t2153/14-t2145*t
     #2172/7-t2149*t2172/7+t2156*t2177/14-t2165*t2160/14+t2169*vV91*t215
     #7/14
      t2186 = vV92*t1486*t1760
      t2193 = t582*vV92
      t2196 = t433*vV92
      vI123pc1 = 3.D0/28.D0*t2169*t2186-t500*t2186/14+3.D0/28.D0*t971*t2
     #186-2.D0/7.D0*t2193*t2146-2.D0/7.D0*t2196*t2146+3.D0/28.D0*t2169*v
     #V92*t2157-t500*vV92*t2157/14+3.D0/28.D0*t971*vV92*t2157-2.D0/7.D0*
     #t2193*t2172-2.D0/7.D0*t2196*t2172
      t2212 = t582*vV93
      t2215 = t433*vV93
      t2219 = vV93*t1486*t1760
      t2222 = t500*vV93
      t2227 = t971*vV93
      vI124pc1 = -t2212*t2146/7-t2215*t2146/7-t500*t2219/14+t2222*t2157/
     #7+t971*t2219/14-t2227*t2157/7+t2227*t2160/7+t2169*t2219/14+t2212*t
     #2172/7+t2215*t2172/7-t2222*t2160/14-t2227*t2177/14-t2169*vV93*t215
     #7/14
      t2246 = vV94*t1486*t1760
      t2253 = t582*vV94
      t2256 = t433*vV94
      vI125pc1 = 3.D0/28.D0*t2169*t2246-t500*t2246/14+3.D0/28.D0*t971*t2
     #246-2.D0/7.D0*t2253*t2146-2.D0/7.D0*t2256*t2146-3.D0/28.D0*t2169*v
     #V94*t2157+t500*vV94*t2157/14-3.D0/28.D0*t971*vV94*t2157+2.D0/7.D0*
     #t2253*t2172+2.D0/7.D0*t2256*t2172
      t2272 = Mzsq*BCsTh
      t2273 = vV95*s
      t2275 = t35*t1677
      t2276 = t1980*t2275
      t2280 = t1980*t1677
      t2283 = Mzsq*t1648
      t2289 = t4*t1677
      t2290 = t1980*t2289
      t2293 = t412*BCsTh
      t2296 = t1486*t1677
      t2297 = t11*t1760*t2296
      vI126pc1 = t2272*t2273*t2276/8+t2272*vV95*t2280/4-t2283*t2273*t227
     #6/8+3*t2272*vV95*t11*t2290-5.D0/4.D0*t2293*vV95*t2297
      t2300 = vV96*s
      vI127pc1 = -t2272*t2300*t2276+4*t2272*vV96*t2280+t2283*t2300*t2276
     #-8*t2272*vV96*t11*t2290+4*t2293*vV96*t2297
      t2315 = t413*BCsTh
      vI128pc1 = 2*t2315*vV97*t1760*t2296
      t2320 = vV98*t1760*t1773
      vI129pc1 = -t2134*t2320/2+t582*BCsTh*t2320/6
      vI130pc1 = t2134*vV99*t1760*t1773/2
      t2331 = BCsTh*t1760*t2296
      vI131pc1 = 16.D0/7.D0*t454*t11*t2331+6.D0/7.D0*t2315*t2280-6.D0/7.
     #D0*t502*t11*t2331
      t2340 = Mzsq*vV100*t1760
      t2347 = t412*vV100
      t2350 = t2347*t1760
      vI132pc1 = 8.D0/7.D0*t2340*t1968-11.D0/7.D0*t2340*t1972-8.D0/7.D0*
     #t2340*t1976-t2347*t1981/2-t2350*t1985/2+t2350*t1989+t2340*t1993/14
     #+3.D0/14.D0*t2340*t2001+t2340*t1997/7
      t2361 = Mzsq*vV101*t1760
      t2365 = t412*vV101
      t2372 = t2365*t1760
      vI133pc1 = -t2361*t1972+t2361*t2001/2-t2365*t1981/2-2*t2361*t1976+
     #t2361*t1997/2-t2372*t1985/2+2*t2361*t1968+t2372*t1989
      t2379 = Mzsq*vV102*t1760
      t2392 = t412*vV102
      t2395 = t2392*t1760
      vI134pc1 = -8.D0/7.D0*t2379*t2023+11.D0/7.D0*t2379*t2027-8.D0/7.D0
     #*t2379*t2030-t2379*t2033/14-t2379*t2041/7+3.D0/14.D0*t2379*t2037+t
     #2392*t2045/2-t2395*t2049/2-t2395*t2053
      t2400 = Mzsq*vV103*t1760
      t2404 = t412*vV103
      t2411 = t2404*t1760
      vI135pc1 = t2400*t2027+t2400*t2037/2+t2404*t2045/2-2*t2400*t2030-t
     #2400*t2041/2-t2411*t2049/2-2*t2400*t2023-t2411*t2053
      t2417 = Mzsq*vV104
      t2418 = t2417*t1980
      t2424 = t412*vV104*t1760
      t2429 = t2417*t1760
      vI136pc1 = t2418*t2075/5+2.D0/5.D0*t2418*t2077+t2424*t2001/5+t2424
     #*t1997/10+t2429*t2087/5+4.D0/5.D0*t2424*t1972+t2424*t1993/10-t2429
     #*t2092/10-t2429*t2095/5-t2429*t2099/10
      t2442 = Mzsq*vV105
      t2443 = t2442*t1980
      t2449 = t412*vV105*t1760
      t2454 = t2442*t1760
      vI137pc1 = -t2443*t2106/5+2.D0/5.D0*t2443*t2108+t2449*t2037/5-t244
     #9*t2041/10-t2454*t2118/5-4.D0/5.D0*t2449*t2027-t2449*t2033/10-t245
     #4*t2123/10+t2454*t2126/5-t2454*t2129/10
      t2467 = t971*BCsTh
      t2472 = BCsTh*vV106
      t2474 = t1774*t1498
      t2478 = t1774*t1558
      vI138pc1 = t2467*vV106*t862*t1980/4-2*t433*t2472*t2474+2*t387*t247
     #2*t2478+t582*t2472*t2474/2
      t2484 = t433*vV107
      t2488 = vV107*t1486*t1760
      t2491 = t971*vV107
      t2498 = t500*vV107
      t2501 = t582*vV107
      vI139pc1 = -t2484*t2146/16-t971*t2488/32-t2491*t2157/16-t2491*t216
     #0/16+t500*t2488/32+t2498*t2157/16-t2501*t2146/16-t2169*t2488/32-t2
     #484*t2172/16-t2491*t2177/32+t2498*t2160/32-t2501*t2172/16-t2169*vV
     #107*t2157/32
      t2518 = vV108*t1486*t1760
      vI140pc1 = -t2169*t2518/16+t500*t2518/24-t971*t2518/16-t2169*vV108
     #*t2157/16+t500*vV108*t2157/24-t971*vV108*t2157/16
      t2534 = t433*vV109
      t2538 = vV109*t1486*t1760
      t2541 = t500*vV109
      t2544 = t582*vV109
      t2549 = t971*vV109
      vI141pc1 = -t2534*t2146/16+t500*t2538/32-t2541*t2157/16-t2544*t214
     #6/16-t971*t2538/32+t2549*t2157/16-t2549*t2160/16-t2169*t2538/32+t2
     #534*t2172/16+t2541*t2160/32+t2544*t2172/16+t2549*t2177/32+t2169*vV
     #109*t2157/32
      t2568 = vV110*t1486*t1760
      vI142pc1 = -t2169*t2568/16+t500*t2568/24-t971*t2568/16+t2169*vV110
     #*t2157/16-t500*vV110*t2157/24+t971*vV110*t2157/16
      t2584 = vV111*t11
      t2586 = t1558*t1677
      t2596 = t1774*t2296
      t2606 = vV111*t14
      t2608 = t1774*t1698
      vI143pc1 = -32.D0/5.D0*t2272*t2584*t1774*t2586+16.D0/5.D0*t2272*vV
     #111*t862*t2290-2.D0/5.D0*t2272*vV111*s*t2596+2*t2293*vV111*t2596-8
     #*t2293*t2584*t1774*t1708-t2272*t2606*t2608/5+t2283*t2606*t2608/5
      t2614 = vV112*s
      vI144pc1 = -t2272*t2614*t2276/3+4.D0/3.D0*t2272*vV112*t2280+t2283*
     #t2614*t2276/3+8.D0/3.D0*t2272*vV112*t11*t2290+4.D0/3.D0*t2293*vV11
     #2*t2297
      t2631 = BCsTh*vV113
      t2634 = t1980*t1678*s
      vI145pc1 = -2*t413*t2631*t2634+16.D0/3.D0*t454*t2631*t2596+2.D0/3.
     #D0*t413*t1648*vV113*t2634
      t2644 = BCsTh*vV114
      t2647 = t1502*t1761*t4
      vI146pc1 = -2*t433*t2644*t2647+3.D0/4.D0*t971*t2644*t1775-t500*t26
     #44*t1775/4
      vI147pc1 = 48.D0/11.D0*t454*BCsTh*t2596+12.D0/11.D0*t502*BCsTh*t25
     #96-20.D0/11.D0*t413*t171*t2596+8.D0/11.D0*t413*t1734*t2596-48.D0/1
     #1.D0*t512*t1717*t2596
      t2671 = vV115*t949
      t2672 = t468*t2671
      t2673 = conjg(Mzsq)
      t2675 = 1/(-s+t2673)
      t2676 = t1486*t2675
      t2677 = t2676*t1835
      t2680 = t437*BCsTh
      t2682 = t2675*t862
      t2683 = t1773*t2682
      t2686 = t171*t2671
      t2687 = t1488*Mzsq
      t2688 = t2676*t2687
      t2691 = t2676*t1488
      t2695 = t2676*t1488*t2673
      t2698 = t862*t2673
      t2699 = t2676*t2698
      vI148pc1 = 7.D0/32.D0*t2672*t2677-7.D0/16.D0*t2680*vV115*t2683-t26
     #86*t2688/4+3.D0/4.D0*t2672*t2691-t2686*t2695/2+7.D0/32.D0*t2672*t2
     #699
      t2705 = vV116*t949
      t2706 = t468*t2705
      t2712 = t171*t2705
      vI149pc1 = t2680*vV116*t2683/2-t2706*t2677/4-t2706*t2699/4-t2706*t
     #2691+t2712*t2688/2+t2712*t2695/2
      t2717 = vV117*t949
      t2718 = t468*t2717
      t2720 = t1525*t2682*Mzsq
      t2724 = t1525*t2682
      t2734 = t1525*t2682*t2673
      vI150pc1 = -t2718*t2720/6+t2680*t2717*t2724/3-t468*vV117*t2683/3+t
     #171*t2717*t2699/3-t2718*t2734/6
      t2737 = vV118*t949
      t2741 = t468*t2737
      t2748 = t171*t2737
      vI151pc1 = -t2680*t2737*t2724/2+t2741*t2720/4+t2741*t2734/4+t468*v
     #V118*t2683-t2748*t2677/2-t2748*t2699/2
      t2753 = t862*t949
      t2754 = t171*t2753
      t2758 = t468*t862
      t2759 = t2675*Mzsq
      t2760 = t1773*t2759
      t2763 = t2680*t862
      t2764 = t2676*t949
      t2767 = t2675*t2673
      t2768 = t1773*t2767
      t2771 = t4*t2673
      t2775 = t171*t2682
      t2784 = t468*t2682
      t2788 = spf5*t2673
      vI152pc1 = t2754*t2676*t413/6-7.D0/96.D0*t2758*t2760+t2763*t2764/2
     #4+t2758*t2768/32-t2754*t2676*t2771/6-t2775*t1773*t1568*t2673/4-t27
     #63*t2676*t949*spf5/16+t2784*t1773*t1568/4+t2784*t1773*t2788/16
      t2793 = t437*vV119*t949
      t2797 = t14*vV119
      t2798 = t2797*t949
      t2806 = s*vV119*t1773
      t2816 = t2797*t1773
      t2817 = Mzsq*t470
      t2821 = t470*t2673
      vI153pc1 = -5.D0/76.D0*t2793*t2676*t862+5.D0/152.D0*t2798*t2677+5.
     #D0/152.D0*t2798*t2699-t2798*t2691/19+t2806*t2682*t413/38+t2806*t26
     #82*t2771/38-t2793*t2676*t1784/76+t2816*t2682*t2817/152+t2816*t2682
     #*t2821/152
      t2825 = vV120*t862
      t2826 = t2825*t2675
      t2830 = t2825*t2676
      t2831 = t949*s
      t2835 = t949*t4
      t2836 = t2673*s
      t2843 = t1773*t433
      t2846 = t2673*t14
      t2847 = t1773*t2846
      t2850 = t437*t470
      t2851 = t1773*t2850
      t2854 = t949*Mzsq
      t2858 = t949*t14
      t2868 = t188*t2673
      t2875 = t2074*Mzsq
      vI154pc1 = t2826*t1773*t442/28-t2830*t2831*t413/56-t2830*t2835*t28
     #36/56+t2826*t1773*t437/56-t2826*t2843/112-t2826*t2847/112+t2826*t2
     #851/56-t2830*t2854*t471/112-t2830*t2858*t2821/112-t2826*t1773*t515
     #/7+t2826*t1773*t454/14+t2826*t1773*t2868/14-t2830*t2835*t471/28+t2
     #830*t2835*t2875/56+t2830*t2835*t2074*t2673/56
      t2883 = vV121*t2675
      t2885 = t1773*t471
      t2888 = t2883*t2753
      t2889 = t1992*t2817
      t2892 = t1992*t2821
      t2896 = vV121*t4*t2675
      t2897 = t2753*t1992
      t2900 = t1486*Mzsq
      t2901 = t2753*t2900
      t2905 = t2753*t1486*t2673
      t2908 = vV121*t35
      t2911 = t2753*t1486*t437
      t2914 = t2908*t2682
      t2922 = t1773*t433*t470
      t2926 = t1773*t471*t2673
      vI155pc1 = -t2883*t862*t2885/24+t2888*t2889/48+t2888*t2892/48-t289
     #6*t2897/6+t2896*t2901/12+t2896*t2905/12-t2908*t2675*t2911/96+t2914
     #*t2843/192+t2914*t2847/192+t2914*t2851/96-t2914*t2922/192-t2914*t2
     #926/192
      t2929 = vV122*t2675
      t2930 = t2929*t862
      t2931 = t1773*t14
      t2940 = vV122*t35
      t2944 = t2940*t2682
      t2956 = vV122*t4*t2675
      t2965 = t2929*t2753
      vI156pc1 = t2930*t2931/28-t2930*t1773*t387/56-t2930*t1773*t2836/56
     #+t2940*t2675*t2911/56-t2944*t2843/112-t2944*t2847/112+t2944*t2851/
     #56-t2944*t2922/112-t2944*t2926/112-t2956*t2897/7+t2956*t2901/14+t2
     #956*t2905/14-t2930*t2885/28+t2965*t2889/56+t2965*t2892/56
      t2996 = t862*spf5
      t3006 = 3.D0/256.D0*t2764*t1491*t413+t2764*t1835*t471/512-t2764*t1
     #835*t14/256-t2764*t1835*t188/64-t2683*t21*t2817/256+t2764*t862*t43
     #7*t470/256-3.D0/128.D0*t2683*t2788*t21*t470+3.D0/128.D0*t2764*t148
     #8*t14-3.D0/512.D0*t2764*t2996*t2850+3.D0/512.D0*t2683*t2788*t471+7
     #.D0/64.D0*t2764*t1789*t2673
      t3041 = 3.D0/512.D0*t2764*t2996*t2846+t2676*t2753*t437/256-3.D0/51
     #2.D0*t2764*t2996*t437-3.D0/32.D0*t2764*t1789*t2788-9.D0/256.D0*t27
     #64*t1488*t2836+7.D0/256.D0*t2683*t21*t2821-3.D0/128.D0*t2764*t1488
     #*t471-3.D0/32.D0*t2764*t1789*s+3.D0/128.D0*t2683*t1568*t471+3.D0/3
     #2.D0*t2764*t1789*s*spf5-3.D0/512.D0*t2764*t2698*t471
      vI157pc1 = t3006+t3041
      vI158pc1 = -vV123*t1486*t949+vV123*t423*t1773*s/9-2.D0/9.D0*vV123/
     #t412*t2931
      t3052 = vV124*t2675
      t3053 = t3052*t1486
      t3054 = t353*t423
      t3058 = t3052*s
      t3059 = t423*t2673
      t3063 = t3052*t11
      t3067 = t3052*t14
      t3076 = Mzsq*t2673
      t3099 = 2.D0/27.D0*t3053*t3054*t2771-t3058*t1971*t3059/27-t3063*t1
     #971*t2771/3+7.D0/54.D0*t3067*t1971*t423*BCsTh+2.D0/27.D0*t3067*t19
     #71*t423*t470-11.D0/54.D0*t3063*t1971*t3076-t3058*t1971*t423*t4/9-1
     #1.D0/27.D0*t3063*t1971*t413+t3052*t1971*t2673/18+11.D0/54.D0*t3052
     #*t1971*Mzsq-t3052*t1993/9+13.D0/27.D0*t3052*t1972-25.D0/54.D0*t305
     #8*t1985
      t3100 = t353*Mzsq
      t3101 = t3100*BCsTh
      t3104 = t1971*t470
      t3107 = t353*t2673
      t3114 = t2000*t4
      t3124 = t3052*t1992
      t3131 = t3052*t1577
      t3132 = t2673*BCsTh
      t3146 = 11.D0/54.D0*t3053*t3101-19.D0/54.D0*t3058*t3104+11.D0/54.D
     #0*t3053*t3107*BCsTh+4.D0/27.D0*t3053*t3107*t470+20.D0/27.D0*t3053*
     #t3114+t3067*t1971*t423/18+2.D0/27.D0*t3052*t1971*t3059*t1975-t3124
     #*t3054*t2821/27-4.D0/27.D0*t3124*t3054*t1975-11.D0/54.D0*t3131*t31
     #00*t3132-2.D0/27.D0*t3124*t3054*t3132+11.D0/27.D0*t3052*t1573*t310
     #0*t2771-10.D0/27.D0*t3131*t1988*t3132
      vI159pc1 = t3099+t3146
      t3147 = vV125*t353
      t3148 = t3147*t1486
      t3149 = s*t2675
      t3150 = t3149*BCsTh
      t3153 = t3149*t470
      t3157 = t2675*BCsTh*t4
      t3159 = t2676*Mzsq
      t3162 = t2759*BCsTh
      t3165 = t11*t2675
      t3166 = t3165*t413
      t3169 = t3165*t3076
      t3172 = t3147*t1577
      t3173 = t2759*t3132
      t3177 = t2675*t4
      t3178 = t3177*t3076
      t3181 = t2767*BCsTh
      t3184 = t2767*t470
      t3187 = t2767*t1975
      vI160pc1 = -t3148*t3150/2-t3148*t3153/2+t3148*t3157+t3147*t3159/4+
     #t3148*t3162/4-t3148*t3166/2-t3148*t3169/4-t3172*t3173/4+t3147*t157
     #3*t3178/2+t3148*t3181/4+t3148*t3184/4-t3172*t3187/2
      t3190 = vV126*t2675
      t3191 = t3190*t14
      t3192 = t1486*t423
      t3196 = t3190*s
      t3197 = t673*t2673
      t3204 = t3190*t11
      t3205 = t2900*t3197
      t3208 = t2026*t2771
      t3211 = t3190*t1486
      t3212 = t423*t673
      t3216 = t2900*t2052
      t3222 = t3190*t1577
      t3223 = t2052*t3132
      t3226 = t2026*t2673
      t3231 = Mzsq*t673
      t3232 = t3231*BCsTh
      t3235 = t3197*BCsTh
      t3238 = -2.D0/27.D0*t3191*t3192*t2040+t3196*t3192*t3197/27+t3196*t
     #3192*t2052/9+11.D0/54.D0*t3204*t3205+t3204*t3208/3-2.D0/27.D0*t321
     #1*t3212*t2771+11.D0/27.D0*t3204*t3216+7.D0/54.D0*t3191*t3192*t2036
     #-10.D0/27.D0*t3222*t3223-t3190*t3226/18-13.D0/27.D0*t3190*t2027+11
     #.D0/54.D0*t3211*t3232+11.D0/54.D0*t3211*t3235
      t3241 = t2026*t470
      t3247 = t2036*t4
      t3250 = t3190*t1992
      t3258 = t3231*t2771
      t3261 = t2900*t673
      t3264 = t3197*t470
      t3270 = t3231*t3132
      t3278 = -25.D0/54.D0*t3196*t2049+19.D0/54.D0*t3196*t3241-t3191*t31
     #92*t673/18+20.D0/27.D0*t3211*t3247+t3250*t3212*t2821/27-4.D0/27.D0
     #*t3250*t3212*t1975-11.D0/27.D0*t3190*t1573*t3258-11.D0/54.D0*t3190
     #*t3261-4.D0/27.D0*t3211*t3264+2.D0/27.D0*t3190*t3192*t3223-11.D0/5
     #4.D0*t3222*t3270-2.D0/27.D0*t3250*t3212*t3132+t3190*t2033/9
      vI161pc1 = t3238+t3278
      t3279 = vV127*t2675
      t3280 = t3279*s
      t3285 = t3279*t1486
      t3291 = t3279*t11
      t3296 = t3279*t1577
      vI162pc1 = -t3280*t2049/2+t3280*t3241/2+t3285*t3247-t3279*t3261/4+
     #t3285*t3232/4+t3291*t3216/2+t3291*t3205/4-t3296*t3270/4-t3279*t157
     #3*t3258/2+t3285*t3235/4-t3285*t3264/4-t3296*t3223/2
      t3308 = s*vV128
      t3313 = t1971*t1648
      t3316 = vV128*t1486
      t3321 = t1988*Mzsq
      t3328 = t3100*t470
      vI163pc1 = -t3308*t1985/26-t3308*t3104/13-t3308*t3313/26-4.D0/13.D
     #0*t3316*t1988+4.D0/13.D0*t11*vV128*t1486*t3321+t3316*t3100/26-4.D0
     #/13.D0*t3316*t3114+t3316*t3328/26+t3316*t3101/13
      t3333 = s*vV129
      t3338 = t2026*t1648
      t3341 = vV129*t1486
      t3346 = t2052*Mzsq
      t3353 = t3231*t470
      vI164pc1 = -t3333*t2049/26+t3333*t3241/13-t3333*t3338/26+4.D0/13.D
     #0*t3341*t2052-4.D0/13.D0*t11*vV129*t1486*t3346-t3341*t3231/26-4.D0
     #/13.D0*t3341*t3247-t3341*t3353/26+t3341*t3232/13
      vI165pc1 = t445*vV130*t1760*t1486/2
      vI166pc1 = t2900*vV131
      vI167pc1 = -vV132*t1486*t949
      t3364 = vV133*t1760*t1486
      t3366 = t971*t470
      vI168pc1 = t2467*t3364/6+t3366*t3364/6
      t3369 = t966*vV134
      t3372 = t437*vV134
      t3373 = t1980*Mzsq
      t3378 = t3372*t1760
      t3379 = t2900*BCsTh
      t3384 = t2900*t470
      vI169pc1 = 3.D0/68.D0*t3369*t2157-t3372*t3373/34+t3369*t1980/68+t3
     #378*t3379/34+t3369*t2160/34+t3378*t3384/17
      t3387 = vV135*t1486
      vI170pc1 = t3387*t412/8-t3387*t468/8-t3387*t387/8-t3387*t466/4+t33
     #87*t2293/8-t3387*t471/8-t3387*t2875/8
      t3402 = vV136*t1486
      t3403 = t423*t437
      t3404 = t3403*BCsTh
      t3419 = t3403*t470
      vI171pc1 = t3402*t3404/26+t3402*t387/13-3.D0/104.D0*t3402*t14+t340
     #2*t3403/52+t3402*t468/104-7.D0/104.D0*t3402*t412+7.D0/52.D0*t3402*
     #t466+t3402*t3419/52+t3402*t471/26-7.D0/104.D0*t3402*t2293+3.D0/52.
     #D0*t3402*t2875
      t3429 = vV137*t1760*t1486
      vI172pc1 = -t2467*t3429/6+t3366*t3429/6
      t3433 = t966*vV138
      t3436 = t437*vV138
      t3441 = t3436*t1760
      vI173pc1 = -3.D0/68.D0*t3433*t2157-t3436*t3373/34+t3433*t1980/68-t
     #3441*t3379/34+t3433*t2160/34+t3441*t3384/17
      t3448 = vV139*t1486
      vI174pc1 = -t3448*t387/8+t3448*t466/4+t3448*t468/8+t3448*t412/8-t3
     #448*t2875/8-t3448*t471/8-t3448*t2293/8
      t3463 = vV140*t1486
      vI175pc1 = -t3463*t3404/26+t3463*t387/13-3.D0/104.D0*t3463*t14+t34
     #63*t3403/52-t3463*t468/104-7.D0/104.D0*t3463*t412-7.D0/52.D0*t3463
     #*t466+t3463*t3419/52+t3463*t471/26+7.D0/104.D0*t3463*t2293+3.D0/52
     #.D0*t3463*t2875
      t3486 = t14*vV141
      t3491 = t11*vV141
      t3492 = t3491*t2676
      t3493 = t1677*t2673
      t3494 = t1016*t3493
      t3497 = t3491*t2675
      t3498 = t2296*t2821
      t3501 = vV141*t2675
      t3502 = t3501*t1486
      t3503 = t2275*t2673
      t3506 = t2900*t1677
      t3509 = t1677*t470
      t3510 = t1016*t3509
      t3515 = t423*t1677
      t3519 = s*vV141
      t3520 = t3519*t2676
      t3521 = t423*t35
      t3522 = t3493*t470
      t3529 = t3519*t2675
      t3533 = 3.D0/130.D0*t3486*t2675*t3192*t2275-9.D0/130.D0*t3492*t349
     #4+16.D0/65.D0*t3497*t3498+t3502*t3503/130-18.D0/65.D0*t3497*t3506-
     #9.D0/130.D0*t3502*t3510+18.D0/65.D0*t3501*t2296+4.D0/65.D0*t3502*t
     #3515*t2673+t3520*t3521*t3522/65-t3520*t3521*t3493/65+8.D0/65.D0*t3
     #529*t3192*t3509
      t3542 = t2900*t3493
      t3545 = t1016*t3522
      t3551 = t2296*t470
      t3554 = t1525*t3509
      t3564 = t1016*t1677
      t3569 = -2.D0/13.D0*t3497*t2296*t2673-4.D0/65.D0*t3502*t3515*t2821
     #+18.D0/65.D0*t15*vV141*t2675*t3542+9.D0/130.D0*t3492*t3545-t3502*t
     #2275*t2821/130-32.D0/65.D0*t3501*t3551+3.D0/130.D0*t3529*t3554-6.D
     #0/65.D0*t3529*t3192*t1677-3.D0/130.D0*t3486*t2676*t3521*t3509+9.D0
     #/130.D0*t3502*t3564-3.D0/130.D0*t3529*t1698
      vI176pc1 = t3533+t3569
      t3570 = vV142*t2675
      t3571 = t3570*t1486
      t3574 = t11*vV142
      t3575 = t3574*t2676
      t3578 = t3574*t2675
      vI177pc1 = t3571*t3564/8-t3575*t3494/8+t3578*t3498/2+t15*vV142*t26
     #75*t3542/2+t3575*t3545/8-t3578*t3506/2-t3571*t3510/8-t3570*t3551
      vI178pc1 = -2*t1582*vV143*t1486*t1677
      t3596 = vV144*t1486
      t3600 = t423*t949
      vI179pc1 = -t3596*t2854-t3596*t2831/9-2.D0/9.D0*t3596*t3600*t14
      vI180pc1 = -Mzsq*vV145*t1773
      t3606 = t1677*Mzsq
      t3610 = s*t1677
      t3617 = t4*t470
      t3618 = t3610*t3617
      t3621 = t3149*t1486
      t3622 = t949*t1677
      t3623 = t3622*t2821
      t3626 = t14*t2675
      t3627 = t3626*t1486
      t3634 = t437*t2675
      t3635 = t3634*t1486
      t3642 = t11*t1677
      t3653 = t3622*t1892
      t3656 = t2675*t949
      t3658 = t1677*t412
      t3662 = t3656*t1577
      t3672 = 19.D0/920.D0*t2764*t3606*t2821-3.D0/115.D0*t2764*t3610*t28
     #17+2.D0/115.D0*t2764*t2289*t2821-4.D0/115.D0*t2764*t3618-11.D0/184
     #0.D0*t3621*t3623-t3627*t3600*t3493/184-3.D0/460.D0*t3627*t3600*t22
     #89-13.D0/1840.D0*t3635*t3600*t3509-t3635*t3600*t1726/1840+7.D0/460
     #.D0*t2764*t3642*t502+21.D0/1840.D0*t2764*t3642*t412*t2673+3.D0/115
     #.D0*t2764*t3606*t3617-t3627*t3653/1840-7.D0/460.D0*t3656*t1573*t36
     #58*t2771-21.D0/1840.D0*t3662*t3658*t2821+13.D0/1840.D0*t3635*t3600
     #*t1677+11.D0/1840.D0*t3621*t3622*t2673
      t3676 = t2289*t2673
      t3679 = t3606*t2673
      t3682 = t3606*t4
      t3685 = t3622*t4
      t3688 = t3622*t470
      t3697 = t3626*t3192
      t3700 = t2771*t470
      t3704 = t3149*t3192
      t3708 = t3622*t3617
      t3714 = t1773*t1677
      t3717 = t2296*t412
      t3726 = t3627*t3622*spf8/1840-7.D0/460.D0*t2764*t3676-19.D0/920.D0
     #*t2764*t3679-4.D0/115.D0*t2764*t3682+3.D0/115.D0*t3621*t3685+3.D0/
     #230.D0*t3627*t3688+21.D0/1840.D0*t2764*t3658*t470+3.D0/115.D0*t276
     #4*t3610*Mzsq+t3697*t3623/184-3.D0/230.D0*t3662*t3606*t3700+t3704*t
     #3622*t2771/230+t3697*t3708/115-t3704*t3622*t3700/230-3.D0/230.D0*t
     #3626*t3714-21.D0/1840.D0*t3656*t3717+t3634*t3192*t3653/1840+3.D0/1
     #15.D0*t3662*t3606*t2771
      vI181pc1 = t3672+t3726
      t3728 = t3192*t949
      vI182pc1 = t171*vV146*t3728
      t3729 = vV147*t353
      t3730 = t3729*t1486
      t3731 = t3165*t188
      t3739 = t15*t2675*t2868
      t3742 = t2676*t4
      t3758 = t3165*t2771
      t3764 = t3729*t1577
      vI183pc1 = -16.D0/35.D0*t3730*t3731-3.D0/35.D0*t3730*t3150-t3730*t
     #3153/14+8.D0/35.D0*t3730*t3739+18.D0/35.D0*t3729*t3742-t3730*t3166
     #/5+t3729*t3159/10-t3729*t3621/70+16.D0/35.D0*t3730*t3157+t3729*t15
     #73*t3178/5+t3730*t3162/10-11.D0/35.D0*t3730*t3758+t3729*t2676*t267
     #3/70-t3764*t3173/10-8.D0/35.D0*t3764*t3187-t3730*t3169/10+3.D0/70.
     #D0*t3730*t3181+t3730*t3184/35
      t3775 = vV148*t353
      t3776 = t3775*t1486
      t3793 = t3775*t1577
      vI184pc1 = t3776*t3150/4+t3776*t3153/4-t3776*t3157-t3775*t3159/8-t
     #3776*t3162/8+t3776*t3166/4-t3775*t3742/2+t3776*t3731+t3776*t3169/8
     #+t3793*t3173/8-t3775*t1573*t3178/4+t3776*t3758/4+t3793*t3187/2-t37
     #76*t3739/2-t3776*t3181/8-t3776*t3184/8
      t3809 = vV149*t2675
      t3810 = t3809*t11
      t3811 = t2026*t188
      t3814 = t3809*s
      t3820 = t2026*t2868
      t3831 = t3809*t1486
      t3843 = t3809*t1577
      vI185pc1 = 16.D0/35.D0*t3810*t3811-3.D0/35.D0*t3814*t2049+t3814*t3
     #241/14-8.D0/35.D0*t3809*t15*t3820-18.D0/35.D0*t3809*t2027+t3810*t3
     #216/5-t3809*t3261/10+t3809*t2033/70+16.D0/35.D0*t3831*t3247-t3809*
     #t1573*t3258/5+t3831*t3232/10+11.D0/35.D0*t3810*t3208-t3809*t3226/7
     #0-t3843*t3270/10-8.D0/35.D0*t3843*t3223+t3810*t3205/10+3.D0/70.D0*
     #t3831*t3235-t3831*t3264/35
      t3854 = vV150*t2675
      t3855 = t3854*s
      t3860 = t3854*t1486
      t3866 = t3854*t11
      t3874 = t3854*t1577
      vI186pc1 = -t3855*t2049/4+t3855*t3241/4+t3860*t3247-t3854*t3261/8+
     #t3860*t3232/8+t3866*t3216/4-t3854*t2027/2+t3866*t3811+t3866*t3205/
     #8-t3874*t3270/8-t3854*t1573*t3258/4+t3866*t3208/4-t3874*t3223/2-t3
     #854*t15*t3820/2+t3860*t3235/8-t3860*t3264/8
      t3891 = s*vV151
      t3898 = vV151*t1486
      vI187pc1 = -t3891*t1985/10-t3891*t3104/5-t3891*t3313/10+t3898*t198
     #8/5+4.D0/5.D0*t11*vV151*t1486*t3321+t3898*t3100/10+2.D0/5.D0*t3898
     #*t3114+t3898*t1996*t4/5+t3898*t3328/10+t3898*t3101/5
      t3916 = s*vV152
      t3923 = vV152*t1486
      vI188pc1 = -t3916*t2049/10+t3916*t3241/5-t3916*t3338/10-t3923*t205
     #2/5-4.D0/5.D0*t11*vV152*t1486*t3346-t3923*t3231/10+2.D0/5.D0*t3923
     #*t3247-t3923*t2040*t4/5-t3923*t3353/10+t3923*t3232/5
      t3945 = BCsTh*vV153
      vI189pc1 = t2467*vV153*t1486*t1774/4-2*t433*t3945*t2474+2*t387*t39
     #45*t2478
      t3952 = BCsTh*vV154
      t3955 = t3952*t1486
      vI190pc1 = t3952*t1626/4-2*t3955*t1488*s+2*t3952*t1509+t3955*t2687
     #/2
      vI191pc1 = -t171*vV155*t3728
      t3965 = t433*vV156
      t3972 = t971*vV156
      vI192pc1 = -t3965*t2146/16-t971*vV156*t1760*t1486/32-t3972*t2157/1
     #6-t3972*t2160/16-t3965*t2172/16-t3972*t2177/32
      vI193pc1 = -t971*vV157*t1760*t1486/16-t971*vV157*t2157/16
      t3987 = vV158*t1486
      t4004 = t171*t4
      vI194pc1 = -t3987*t21/16-t3987*t412/32-t3987*t14/32-t3987*t468/16-
     #t3987*t471/16+t3987*t387/32+t3987*t466/16-t3987*t413/16-t3987*t400
     #4/16-t3987*t2293/32-t3987*t2098/32+t3987*t2875/32-t3987*t2315/16
      t4015 = vV159*t1486
      vI195pc1 = -t4015*t412/16+t4015*t387/24-t4015*t14/16-t4015*t2293/1
     #6+t4015*t466/24-t4015*t468/16
      t4028 = t433*vV160
      t4035 = t971*vV160
      vI196pc1 = -t4028*t2146/16-t971*vV160*t1760*t1486/32+t4035*t2157/1
     #6-t4035*t2160/16+t4028*t2172/16+t4035*t2177/32
      vI197pc1 = -t971*vV161*t1760*t1486/16+t971*vV161*t2157/16
      t4050 = vV162*t1486
      vI198pc1 = -t4050*t21/16-t4050*t14/32+t4050*t468/16-t4050*t471/16+
     #t4050*t387/32-t4050*t466/16-t4050*t413/16-t4050*t412/32+t4050*t400
     #4/16+t4050*t2098/32+t4050*t2875/32+t4050*t2315/16+t4050*t2293/32
      t4077 = vV163*t1486
      vI199pc1 = -t4077*t412/16+t4077*t387/24-t4077*t14/16+t4077*t2293/1
     #6-t4077*t466/24+t4077*t468/16
      t4090 = BCsTh*vV164
      t4092 = t188*t1677
      t4101 = t4090*t1486
      t4105 = t4090*t1577
      t4115 = t2682*t3503
      t4118 = t14*t1486
      t4120 = t2682*t2275
      t4123 = t1648*vV164
      vI200pc1 = -2*t4090*t1573*t2682*t4092*t2673+t4090*s*t2676*t1678/4-
     #2*t4101*t2682*t2289+4*t4105*t2682*t4092-5.D0/2.D0*t4090*t1573*t267
     #5*t1488*t3679-t4090*t1992*t4115/16+t4090*t4118*t4120/8-t4123*t4118
     #*t4120/8-5.D0/8.D0*t4101*t2682*t3606+5.D0/2.D0*t4105*t2682*t3682+t
     #4105*t2682*t3676-t4101*t2682*t3493/8+t4123*t1992*t4115/16+5.D0/8.D
     #0*t4105*t2682*t3679
      t4144 = BCsTh*vV165
      t4145 = t4144*t11
      t4146 = t2676*t2289
      t4149 = t4144*t1573
      t4153 = t1648*vV165
      t4156 = t2675*t35*t3493
      t4159 = t2676*t1677
      t4166 = t2676*t2275
      vI201pc1 = -4.D0/3.D0*t4145*t4146+2.D0/3.D0*t4149*t3177*t3493+t415
     #3*t1486*t4156/12-2.D0/3.D0*t4144*t4159-t4145*t2676*t3606/3+t4144*s
     #*t4166/6-t4153*s*t4166/6-t4144*t1486*t4156/12+t4149*t2675*t1677*t3
     #076/3+t4145*t2676*t3493/3
      vI202pc1 = -2*t1975*vV166*t1686+16.D0/3.D0*t188*BCsTh*t11*vV166*t1
     #486*t1678+2.D0/3.D0*t4*t1648*vV166*t1686
      vI203pc1 = -t466*vV167*t1760*t1486
      t4201 = t171*vV168
      vI204pc1 = t468*vV168*t1502*t949-8.D0/3.D0*t4201*t1502*t2835-t4201
     #*t2901/3
      vI205pc1 = t171*vV169*t949*t1486
      t4209 = t2098*t862
      t4212 = t1717*t2682
      t4218 = BCsTh*t862*t2675
      t4219 = t2296*t2771
      t4235 = t2296*t1751
      t4248 = t2676*t1726
      vI206pc1 = t4209*t4159/4+t4212*t2296*t413*t2673/4-5.D0/12.D0*t4218
     #*t4219+t4212*t2296*t2868+t1648*t862*t2675*t4219/6-t2758*t4159/4-t4
     #218*t2296*t413/4+t4212*t2296*t454+t2775*t4235/3-t1734*t2682*t4235/
     #3+3.D0/2.D0*t171*t862*t4146-2*t4218*t2586-t1734*t862*t4146+t4209*t
     #4248/8-t2758*t4248/8-t15*BCsTh*t2682*t2296*t454*t2673
      t4261 = t171*vV170
      vI207pc1 = -t468*vV170*t2764/2+t4261*t2760/4+t4261*t2768/4
      t4266 = spf8*Mzsq
      t4270 = t1773*Mzsq
      t4273 = spf8*t2673
      vI208pc1 = -7.D0/48.D0*t1713*t3656*t4266-5.D0/48.D0*t3150*t4270-7.
     #D0/48.D0*t1713*t3656*t4273+7.D0/24.D0*t468*t1486*t3656*spf8-5.D0/4
     #8.D0*t3150*t1773*t2673+5.D0/24.D0*t468*t2764
      t4286 = vV171*t1486
      t4287 = t4286*t2675
      vI209pc1 = -t4287*t2835*s/7+t4287*t2835*Mzsq/14+t4287*t2835*t2673/
     #14-t4286*t3656*t14/28+t4287*t2854*s/56+t4287*t949*t2673*s/56-t4287
     #*t2858*t470/28+t4287*t2854*t2074/56+t4287*t2831*t2821/56
      vI210pc1 = 3.D0/76.D0*t1773*t2759*s-3.D0/38.D0*t1773*t3626+3.D0/76
     #.D0*t1773*t3149*t2673+t2764*t4273*t4/38+t2764*t4273*t2074/152+t276
     #4*t4266*t2074/152+t2764*t1751*Mzsq/38-t2764*t1751*s/19+5.D0/152.D0
     #*t2764*t4273*s+5.D0/152.D0*t2764*t4266*s-5.D0/76.D0*t1773*t2675*sp
     #f8*t14-t2764*spf8*t14*t470/76
      t4353 = BCsTh*vV172
      vI211pc1 = 7.D0/16.D0*t971*t4353*t1775-3.D0/4.D0*t433*t4353*t2647
      t4360 = BCsTh*vV173
      vI212pc1 = -t971*t4360*t1775/2+t433*t4360*t2647
      t4366 = BCsTh*vV174
      t4370 = t1773*t1760*t35*t862
      vI213pc1 = -t971*t4366*t4370/3+t433*t4366*t1775/3
      t4375 = BCsTh*vV175
      vI214pc1 = t971*t4375*t4370/2-t433*t4375*t1775
      t4381 = t1833*Mzsq
      t4385 = t2680*t1833
      vI215pc1 = -t468*t4381*t1866/4-t4385*t1845/24+t4385*t1835*t1562/16
      vI216pc1 = t2134*vV176*t1760*t1773/2
      vI217pc1 = -7.D0/24.D0*t468*t1833*t2900*spf8-5.D0/24.D0*t468*t4381
     #*t1486
      vI218pc1 = -3.D0/58.D0*t437*vV177*t423*t1766-s*vV177*Mzsq*t1766/2
      t4409 = vV178*t1760
      t4410 = t4409*t1486
      t4411 = s*t353
      t4441 = t14*t353
      vI219pc1 = 9.D0/11.D0*t4410*t4411*t2817+t4410*t4411*t2272-37.D0/44
     #.D0*t4409*t1971*t412+3.D0/22.D0*t4410*t4411*t4+37.D0/22.D0*t4410*t
     #1359*t502-12.D0/11.D0*t4410*t3321+2.D0/11.D0*t4410*t4411*Mzsq-37.D
     #0/44.D0*t4410*t353*t412*BCsTh-3.D0/44.D0*t4409*t4118*t353-18.D0/11
     #.D0*t4410*t1988*t2272-7.D0/44.D0*t4410*t4441*BCsTh-t4410*t4441*t47
     #0/11+2.D0/11.D0*t4410*t2077
      t4450 = t412*vV179
      t4453 = Mzsq*vV179*t1760
      t4455 = t4450*t1760
      vI220pc1 = -t4450*t1981+t4453*t2001-t4455*t1985+t4453*t1997+2*t445
     #5*t1989-2*t4453*t1976
      t4462 = vV180*t1760
      t4463 = t4462*t1486
      t4464 = s*t673
      t4495 = t14*t673
      vI221pc1 = -9.D0/11.D0*t4463*t4464*t2817+t4463*t4464*t2272+37.D0/4
     #4.D0*t4462*t2026*t412-3.D0/22.D0*t4463*t4464*t4-37.D0/22.D0*t4463*
     #t11*t673*t502+12.D0/11.D0*t4463*t3346-2.D0/11.D0*t4463*t4464*Mzsq-
     #37.D0/44.D0*t4463*t673*t412*BCsTh+3.D0/44.D0*t4462*t4118*t673-18.D
     #0/11.D0*t4463*t2052*t2272-7.D0/44.D0*t4463*t4495*BCsTh+t4463*t4495
     #*t470/11+2.D0/11.D0*t4463*t2108
      t4505 = Mzsq*vV181*t1760
      t4507 = t412*vV181
      t4510 = t4507*t1760
      vI222pc1 = t4505*t2037+t4507*t2045-t4505*t2041-t4510*t2049-2*t4505
     #*t2030-2*t4510*t2053
      t4517 = t412*vV182*t1760
      t4522 = Mzsq*vV182
      t4523 = t4522*t1760
      vI223pc1 = t4517*t1997/14+t4517*t2001/7-2.D0/7.D0*t4523*t2087+2.D0
     #/7.D0*t4517*t1972+t4517*t1993/14-t4523*t2099/14-2.D0/7.D0*t4517*t1
     #976-t4523*t2092/14-t4523*t2095/7-2.D0/7.D0*t4522*t1980*t2077
      t4542 = t412*vV183*t1760
      t4547 = Mzsq*vV183
      t4548 = t4547*t1760
      vI224pc1 = -t4542*t2041/14+t4542*t2037/7+2.D0/7.D0*t4548*t2118-2.D
     #0/7.D0*t4542*t2027-t4542*t2033/14-t4548*t2129/14-2.D0/7.D0*t4542*t
     #2030-t4548*t2123/14+t4548*t2126/7-2.D0/7.D0*t4547*t1980*t2108
      t4566 = t412*vV184
      vI225pc1 = t4566*t2146/2-t4566*t1980*s/2+Mzsq*vV184*t1760*t1630/2
      vI226pc1 = -t387*vV185*t1766/2
      t4578 = vV186*t1760*t1486
      t4581 = t971*vV186
      t4586 = t500*vV186
      vI227pc1 = t2169*t4578/8-t4581*t2157/8-t500*t4578/8-t4586*t2157/4+
     #t2169*vV186*t2157/8-t4581*t2160/8-t4586*t2160/8
      t4596 = t437*vV187
      t4597 = t4596*t1760
      t4600 = t14*vV187
      t4601 = t1980*t412
      t4606 = t4600*t1760
      t4607 = t1486*t412
      t4608 = t4607*BCsTh
      t4611 = t966*vV187
      t4616 = s*vV187
      t4617 = t1980*t2168
      t4622 = t4607*t470
      t4628 = t1586*t2168
      vI228pc1 = -t4597*t3379/14-3.D0/26.D0*t4600*t4601+3.D0/182.D0*t459
     #6*t3373-7.D0/26.D0*t4606*t4608-3.D0/182.D0*t4611*t1980-t4611*t2157
     #/26+3.D0/26.D0*t4616*t4617-8.D0/91.D0*t4597*t3384-2.D0/13.D0*t4606
     #*t4622-2.D0/91.D0*t4611*t2160+3.D0/26.D0*t4616*t1760*t4628
      t4632 = vV188*t1760*t1486
      t4635 = t500*vV188
      t4638 = t971*vV188
      vI229pc1 = -t500*t4632/8+t4635*t2157/4+t4638*t2157/8+t2169*t4632/8
     #-t4635*t2160/8-t4638*t2160/8-t2169*vV188*t2157/8
      t4650 = t437*vV189
      t4651 = t4650*t1760
      t4654 = t14*vV189
      t4659 = t4654*t1760
      t4662 = t966*vV189
      t4667 = s*vV189
      vI230pc1 = t4651*t3379/14-3.D0/26.D0*t4654*t4601+3.D0/182.D0*t4650
     #*t3373+7.D0/26.D0*t4659*t4608-3.D0/182.D0*t4662*t1980+t4662*t2157/
     #26+3.D0/26.D0*t4667*t4617-8.D0/91.D0*t4651*t3384-2.D0/13.D0*t4659*
     #t4622-2.D0/91.D0*t4662*t2160-3.D0/26.D0*t4667*t1760*t4628
      t4679 = s*vV190
      t4682 = vV190*t1760
      t4683 = t4682*t1486
      t4693 = t14*vV190*t1760
      vI231pc1 = 3.D0/22.D0*t4679*t2280+29.D0/88.D0*t4683*t2275*t412*t47
     #0+29.D0/22.D0*t11*vV190*t1760*t3717-3.D0/88.D0*t4693*t1698+3.D0/88
     #.D0*t4693*t3554-6.D0/11.D0*t4682*t3506-2.D0/11.D0*t4679*t1760*t355
     #1+16.D0/11.D0*t4683*t3606*t470-29.D0/88.D0*t4683*t2275*t412
      t4714 = t412*vV191*t1760
      vI232pc1 = 2*Mzsq*vV191*t1760*t3551-t4714*t1698/2+2*t412*t11*vV191
     #*t2280+t4714*t3554/2
      t4723 = t413*vV192
      vI233pc1 = t4723*t2280+t4723*t1980*t3509/3
      vI234pc1 = -s*vV193*t1760*t1773*t412/2-4.D0/29.D0*t14*vV193*t1760*
     #t4270-3.D0/58.D0*t437*vV193*t1766
      vI235pc1 = -t582*vV194*t1766/2
      t4765 = t14*t1677
      vI236pc1 = -3.D0/44.D0*t2168*t1760*t3714+t4617*t949*t11*t2289/11+3
     #.D0/44.D0*t4617*t3688+5.D0/44.D0*t4601*t2831*t1677-5.D0/44.D0*t460
     #1*t2831*t3509+8.D0/77.D0*t4601*t3708-t4601*t3685/7+5.D0/77.D0*t176
     #6*t3610*t413-3.D0/308.D0*t1766*t4765*Mzsq+3.D0/308.D0*t1766*t4765*
     #t2817-4.D0/33.D0*Mzsq*t1760*t1773*t3618-5.D0/308.D0*t1766*t437*t16
     #77*t470+5.D0/308.D0*t1980*t949*t437*t1677+4.D0/231.D0*t1766*t4765*
     #t3617-t1766*t4765*t4/77
      A1e0c1p1t1v = -144*vI97pc1+42*vI98pc1
      A1e0c2p1t1v = -168*vI79pc1+192*vI80pc1+720*vI81pc1-288*vI82pc1+288
     #*vI83pc1
      A1e0c3p1t1v = 432*vI84pc1-288*vI85pc1-432*vI86pc1+288*vI87pc1+288*
     #vI88pc1+192*vI89pc1-54*vI90pc1-192*vI91pc1+54*vI92pc1+144*vI93pc1+
     #108*vI94pc1+24*vI95pc1+36*vI96pc1
      A1e0c4p1t1v = -144*vI111pc1+42*vI112pc1
      A1e0c5p1t1v = 456*vI106pc1+192*vI107pc1+720*vI108pc1-288*vI109pc1+
     #288*vI110pc1
      t4817 = 144*vI113pc1
      t4825 = 144*vI121pc1
      t4834 = 576*vI130pc1
      A1e0c6p1t1v = t4817+504*vI132pc1-288*vI133pc1-504*vI134pc1+288*vI1
     #35pc1-180*vI136pc1+180*vI137pc1+192*vI138pc1+t4825+2304*vI139pc1+8
     #64*vI140pc1-2304*vI141pc1+864*vI142pc1+180*vI143pc1+108*vI144pc1+8
     #64*vI145pc1+576*vI146pc1+t4834+33*vI147pc1
      A1e0c7p1t1v = 48*vI104pc1-36*vI105pc1
      A1e0c8p1t1v = 48*vI99pc1+48*vI100pc1+72*vI101pc1-72*vI102pc1+288*v
     #I103pc1
      A1e0c9p1t1v = t4817-648*vI114pc1+288*vI115pc1+648*vI116pc1-288*vI1
     #17pc1-108*vI118pc1+108*vI119pc1+48*vI120pc1+t4825-1008*vI122pc1-50
     #4*vI123pc1+1008*vI124pc1-504*vI125pc1+288*vI126pc1+36*vI127pc1+576
     #*vI128pc1+864*vI129pc1+t4834+42*vI131pc1
      A1e0c10p1t1v = 672*vI209pc1-228*vI210pc1
      A1e0c11p1t1v = 912*vI153pc1+1344*vI154pc1-3456*vI155pc1-2016*vI156
     #pc1+1536*vI157pc1
      t4878 = 144*vI206pc1+288*vI205pc1+432*vI204pc1-2304*vI198pc1+864*v
     #I199pc1+576*vI200pc1+432*vI201pc1+864*vI202pc1+24*vI203pc1+144*vI1
     #97pc1-192*vI196pc1+864*vI195pc1
      t4892 = 2304*vI194pc1+192*vI192pc1-144*vI193pc1+72*vI191pc1+192*vI
     #190pc1+180*vI188pc1+288*vI189pc1-180*vI187pc1+1152*vI186pc1+1152*v
     #I184pc1-2520*vI185pc1+2520*vI183pc1+72*vI182pc1
      A1e0c12p1t1v = t4878+t4892
      A1e0c1p1t2v = A1e0c1p1t1v
      A1e0c2p1t2v = A1e0c2p1t1v
      A1e0c3p1t2v = A1e0c3p1t1v
      A1e0c4p1t2v = A1e0c4p1t1v
      A1e0c5p1t2v = A1e0c5p1t1v
      A1e0c6p1t2v = A1e0c6p1t1v
      A1e0c7p1t2v = A1e0c7p1t1v
      A1e0c8p1t2v = A1e0c8p1t1v
      A1e0c9p1t2v = A1e0c9p1t1v
      A1e0c10p1t2v = A1e0c10p1t1v
      A1e0c11p1t2v = A1e0c11p1t1v
      A1e0c12p1t2v = A1e0c12p1t1v
      A2e0c1p1t1v = 192*vI216pc1+288*vI217pc1
      A2e0c2p1t1v = 768*vI211pc1+192*vI212pc1+864*vI213pc1+288*vI214pc1+
     #576*vI215pc1
      A2e0c3p1t1v = -2088*vI218pc1-1584*vI219pc1+576*vI220pc1-1584*vI221
     #pc1+576*vI222pc1+1008*vI223pc1+1008*vI224pc1+384*vI225pc1+576*vI22
     #6pc1+2304*vI227pc1+3276*vI228pc1+2304*vI229pc1-3276*vI230pc1-1584*
     #vI231pc1+288*vI232pc1+3456*vI233pc1+4176*vI234pc1-2304*vI235pc1-27
     #72*vI236pc1
      A2e0c4p1t1v = 96*vI207pc1+144*vI208pc1
      A2e0c5p1t1v = 384*vI148pc1+96*vI149pc1+432*vI150pc1+144*vI151pc1+2
     #88*vI152pc1
      t4938 = -324*vI158pc1-1944*vI159pc1+576*vI160pc1-1944*vI161pc1+576
     #*vI162pc1+468*vI163pc1+468*vI164pc1+288*vI165pc1-48*vI166pc1+72*vI
     #167pc1-36*vI168pc1+204*vI169pc1
      t4951 = 576*vI170pc1-1872*vI171pc1-36*vI172pc1+204*vI173pc1+576*vI
     #174pc1+1872*vI175pc1-2340*vI176pc1+288*vI177pc1-576*vI178pc1+648*v
     #I179pc1-288*vI180pc1+5520*vI181pc1
      A2e0c6p1t1v = t4938+t4951
      A2e0c1p1t2v = A2e0c1p1t1v
      A2e0c2p1t2v = A2e0c2p1t1v
      A2e0c3p1t2v = A2e0c3p1t1v
      A2e0c4p1t2v = A2e0c4p1t1v
      A2e0c5p1t2v = A2e0c5p1t1v
      A2e0c6p1t2v = A2e0c6p1t1v

      RETURN
      END

      SUBROUTINE CrsecSeqn_msv()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      t1 = Eq**2
      t2 = t1**2
      t3 = t2*Er
      t4 = BCsThl*CsThr
      t5 = t3*t4
      t6 = BSnThl*SnThr
      t7 = cos(DBLE(phil))
      t8 = t7*El
      t9 = t6*t8
      t11 = 64*t5*t9
      t12 = t1*Er
      t13 = t12*t4
      t14 = sqq*El
      t15 = t6*t7
      t16 = t14*t15
      t17 = t13*t16
      t19 = t1*Eq
      t20 = t19*El
      t21 = t4*t20
      t22 = t7*sqq
      t23 = t6*t22
      t24 = t21*t23
      t26 = Er**2
      t27 = t19*t26
      t28 = t27*t4
      t30 = 32*t28*t9
      t31 = sqq**2
      t32 = t31*sqq
      t33 = t31*t1
      t34 = 4*t33
      t35 = t2*sqq
      t36 = 8*t35
      t37 = Ml**2
      t38 = t37*t31
      t39 = 2*t38
      t40 = t26*Er
      t41 = t19*t40
      t42 = 8*t41
      t43 = t2*Eq
      t44 = t43*Er
      t45 = 32*t44
      t46 = CsThr**2
      t47 = BCsThl**2
      t48 = t46*t47
      t49 = t44*t48
      t50 = 32*t49
      t51 = t44*t4
      t52 = 64*t51
      t53 = Eq*Er
      t54 = t37*sqq
      t55 = t53*t54
      t56 = 8*t55
      t57 = t1*t26
      t58 = t47*sqq
      t59 = t57*t58
      t60 = 2*t59
      t61 = t4*t35
      t63 = t46*t19
      t65 = t63*t40*t47
      t67 = t2*t47
      t68 = t46*sqq
      t69 = t67*t68
      t71 = t46*t1
      t72 = t26*t37
      t73 = t71*t72
      t75 = t31*BCsThl
      t76 = CsThr*t1
      t77 = t75*t76
      t78 = 4*t77
      t79 = -t11-8*t17+16*t24+t30-t32+t34-t36-t39+t42+t45+t50-t52+t56-t6
     #0+16*t61+8*t65-8*t69-8*t73-t78
      t80 = t2*t26
      t81 = t80*t48
      t82 = 32*t81
      t83 = t80*t4
      t84 = 32*t83
      t85 = t57*t68
      t87 = t12*sqq
      t88 = t87*t9
      t89 = 16*t88
      t90 = t63*t40
      t91 = 8*t90
      t92 = t20*BSnThl
      t93 = SnThr*t7
      t94 = t93*sqq
      t95 = t92*t94
      t97 = t19*Er
      t98 = El**2
      t99 = t97*t98
      t100 = BSnThl**2
      t101 = SnThr**2
      t102 = t100*t101
      t103 = t7**2
      t104 = t102*t103
      t106 = 32*t99*t104
      t107 = t3*BSnThl
      t108 = t93*El
      t110 = 64*t107*t108
      t111 = t31*El
      t112 = t111*BSnThl
      t114 = t112*t93*Eq
      t115 = 4*t114
      t116 = t98*t100
      t117 = t116*t101
      t118 = t103*sqq
      t120 = t117*t118*t1
      t122 = t4*sqq
      t123 = t97*t122
      t125 = t4*t31
      t127 = 2*t53*t125
      t128 = t48*sqq
      t129 = t97*t128
      t131 = t57*t122
      t132 = 8*t131
      t133 = t26*t47
      t136 = 2*t71*t133*sqq
      t137 = t41*t47
      t139 = t57*t37
      t141 = t53*t31
      t142 = 6*t141
      t143 = t97*sqq
      t144 = 16*t143
      t145 = t57*sqq
      t147 = -t82+t84-2*t85-t89+t91-16*t95+t106+t110+t115-8*t120+8*t123+
     #t127+8*t129-t132-t136+8*t137-8*t139+t142-t144-10*t145
      t148 = t79+t147
      t149 = 1/t26
      t151 = 1/t1
      t152 = sll**2
      t153 = 1/t152
      t156 = 1/(-1+t46)
      A1e0c1p1t1n = 24*t148*t149*t151*t153*t156
      t159 = t47*BCsThl
      t160 = t19*t159
      t161 = CsThr*sll
      t164 = 16*t160*t161*t9
      t165 = t12*t14
      t167 = t6*t7*t47
      t169 = 16*t165*t167
      t170 = t1*t159
      t171 = CsThr*Er
      t173 = sll*El
      t174 = t173*t15
      t176 = 8*t170*t171*t174
      t177 = t7*t37
      t178 = t6*t177
      t180 = 64*t21*t178
      t181 = t19*BCsThl
      t184 = 16*t181*t161*t9
      t186 = 8*t13*t174
      t187 = t2*t159
      t188 = t187*t171
      t190 = 64*t188*t9
      t191 = sqq*sll
      t192 = Eq*El
      t194 = t191*t192*t167
      t195 = 4*t194
      t196 = CsThr*t159
      t197 = t27*t196
      t199 = 32*t197*t9
      t200 = 4*t38
      t201 = t37*t2
      t202 = 32*t201
      t203 = t37**2
      t205 = 8*t203*sqq
      t206 = sll*t2
      t207 = 8*t206
      t208 = sll*t31
      t209 = t47**2
      t210 = t43*t209
      t211 = t46*Er
      t213 = 32*t210*t211
      t214 = t11-t30-t164-t169+t176-t180+t184-t186-t190-t195+t199+t200+t
     #42-t45+t202+t205-t207-t208+t213
      t216 = 32*t201*t48
      t217 = t43*t159
      t219 = 64*t217*t171
      t220 = t2*BCsThl
      t222 = 16*t220*t161
      t224 = 16*t187*t161
      t225 = CsThr*t37
      t227 = 64*t220*t225
      t228 = t46*t209
      t230 = 32*t80*t228
      t231 = t40*t209
      t233 = 8*t63*t231
      t234 = t46*sll
      t235 = t57*t234
      t236 = 2*t235
      t237 = t209*sqq
      t239 = 4*t57*t237
      t242 = 2*t57*t209*sll
      t243 = t37*sll
      t245 = 2*t243*t58
      t246 = t31*t47
      t247 = t53*t246
      t248 = 2*t247
      t249 = t1*t47
      t250 = t191*t249
      t251 = 4*t250
      t252 = t80*t196
      t253 = 32*t252
      t255 = 8*t48*t206
      t257 = 8*t228*t206
      t258 = t97*t58
      t259 = 16*t258
      t260 = t216-t219+t222-t224-t227-t230+t233-t236+t239-t242+t245+t248
     #-t251+t253-t255+t257-t50+t52-t259+t82
      t262 = sll*t1
      t263 = t262*t98
      t265 = t102*t103*t47
      t267 = 8*t263*t265
      t269 = 32*t99*t265
      t271 = t54*Eq*t9
      t272 = 16*t271
      t273 = t191*Eq
      t274 = t273*t9
      t275 = 4*t274
      t276 = El*t47
      t279 = 64*t107*t93*t276
      t280 = sll*t19
      t281 = t280*El
      t283 = 16*t281*t167
      t284 = t196*sqq
      t285 = t57*t284
      t286 = 8*t285
      t287 = t171*sqq
      t289 = 16*t160*t287
      t290 = t171*sll
      t292 = 8*t160*t290
      t293 = t1*BCsThl
      t294 = t225*sqq
      t295 = t293*t294
      t296 = 16*t295
      t297 = CsThr*sqq
      t298 = t297*sll
      t299 = t293*t298
      t300 = 4*t299
      t301 = t170*t298
      t302 = 4*t301
      t305 = 2*t57*t228*sll
      t306 = t19*t47
      t307 = t211*sll
      t309 = 8*t306*t307
      t310 = t19*t209
      t312 = 8*t310*t307
      t313 = t181*t290
      t314 = 8*t313
      t315 = Eq*BCsThl
      t316 = t315*CsThr
      t317 = Er*sqq
      t318 = t317*sll
      t319 = t316*t318
      t320 = 2*t319
      t321 = -t84+t267+t269-t272+t275+t279+t283+t89-t91-t286+t289+t292+t
     #296-t300+t302+t305+t309-t312-t314+t320
      t322 = Eq*t159
      t323 = t322*CsThr
      t324 = t323*t318
      t325 = 2*t324
      t326 = t206*t47
      t327 = 8*t326
      t328 = t37*t19
      t331 = 64*t328*El*t15
      t333 = 8*t41*t209
      t334 = t54*t1
      t335 = 16*t334
      t336 = t57*sll
      t337 = 2*t336
      t338 = t208*t47
      t340 = 16*t281*t15
      t342 = 8*t263*t104
      t343 = t37*t1
      t346 = 32*t343*t98*t104
      t347 = t44*t47
      t348 = 32*t347
      t349 = 16*t123
      t350 = t191*t1
      t351 = 4*t350
      t352 = t243*sqq
      t353 = 2*t352
      t354 = 2*t141
      t355 = 4*t145
      t356 = -t325+t327+t331-t333-t335+t337+t338-t106-t110-t340-t342+t34
     #6+t348-t349+t132+t351-t353-t354+t144-t355
      t358 = t214+t260+t321+t356
      t360 = -1+t47
      t361 = t360**2
      t362 = 1/t361
      t363 = 1/t31
      A1e0c2p1t1n = -24*t358*t151*t362*t363*t149
      t370 = t2*t1
      t371 = 32*t370
      t375 = t46*El
      t376 = t12*t375
      t385 = t4*t97
      t386 = t101*t103
      t387 = t116*t386
      t391 = t1*t98
      t396 = sqq*t19
      t403 = Eq*t31
      t404 = El*BSnThl
      t417 = t46*CsThr
      t418 = t159*t417
      t421 = t97*t234
      t423 = t46*t37
      t428 = -8*t391*t102*t118*t4+16*t48*t396*t9-8*t12*t48*t16+4*t403*t4
     #04*t93*t4+4*t165*t6*t7*t46-8*t27*t48*t9+8*t80-32*t44*t418-4*t421+1
     #6*t97*t423-4*t48*t33
      t430 = t417*sqq
      t435 = 8*t33
      t436 = 24*t35
      t437 = t41*t4
      t439 = t417*BCsThl
      t445 = t31*t46
      t446 = t53*t445
      t447 = t97*t68
      t449 = -8*t187*t430-2*t225*t75-t435+t436+t39+16*t437-16*t41*t439+8
     #*t80*t418-t53*t191-t446+4*t447
      t461 = 24*t80*t439+8*t252+64*t49-32*t51+t60-56*t61+40*t69+12*t77-8
     #*t81-8*t83+6*t85
      t464 = t12*sll
      t473 = t4*t2
      t476 = t27*El
      t481 = t439*sqq
      t487 = t418*sqq
      t490 = t68*sll
      t492 = 4*t464*t9-24*t27*t46*t9+192*t4*t43*t9+96*t473*t387-8*t476*t
     #167-8*t271+4*t88-4*t97*t481-16*t97*t439*t37+8*t97*t487+t53*t490
      t493 = t48*t31
      t496 = t2*t98
      t504 = BCsThl*sll
      t505 = t504*CsThr
      t508 = t439*sll
      t515 = t98*El
      t516 = t19*t515
      t517 = t100*BSnThl
      t518 = t101*SnThr
      t519 = t517*t518
      t520 = t103*t7
      t524 = t43*El
      t528 = 2*t53*t493-96*t496*t104-4*t53*t54*t46+2*t57*t481+4*t57*t505
     #-4*t57*t508-2*t57*t487+4*t97*t508-32*t516*t519*t520-96*t524*t15-2*
     #t285
      t532 = Er*t37
      t533 = t532*sqq
      t536 = t476*t15
      t538 = t48*t370
      t540 = 8*t334
      t544 = CsThr*t32
      t545 = t544*BCsThl
      t546 = t4*t370
      t548 = 8*t295-4*t313+4*t316*t533+8*t536-96*t538-t540+48*t95-8*t114
     #+24*t120-t545+96*t546
      t550 = t370*t159*t417
      t555 = t80*t46
      t557 = t80*t47
      t560 = t97*sll
      t563 = 32*t550+12*t123-24*t129-6*t131+t136-24*t555-8*t557-t141+4*t
     #143+4*t560-2*t145
      t567 = 1/sqq
      t569 = 1/sll
      t571 = t151*t149
      t572 = 1/t360
      A1e0c3p1t1n = -48*(-t11+16*t17-64*t24-t371+t32-4*t12*t234*t9+16*t3
     #76*t178-96*t48*t43*t9+64*t3*t48*t9-32*t385*t387+t428+t449+t461+t49
     #2+t528+t548+t563)*t567*t569*t571*t572*t156
      t580 = 1/(-sll+Mzsq)
      t582 = sll*Mzsq
      t583 = abs(DCMPLX(Mzsq))
      t584 = t583**2
      t586 = 1/(t582-t584)
      t588 = t151*t580*t586*t156
      A1e0c4p1t1n = -24*Mzsq*t148*t149*t588
      t594 = 1/(-sqq+Mzsq)
      t595 = t149*t594
      t596 = sqq*Mzsq
      t598 = 1/(t596-t584)
      t600 = t595*t598*t151
      A1e0c5p1t1n = 24*Mzsq*t358*t362*t600
      t603 = conjg(Mzsq)
      t605 = t6*t22*t603
      t608 = t47*Mzsq
      t609 = t27*t608
      t612 = t48*Mzsq
      t613 = t27*t612
      t614 = t7*t603
      t615 = t614*SnThr
      t616 = t404*t615
      t621 = t48*El
      t622 = t27*t621
      t623 = t7*sll
      t625 = t6*t623*sqq
      t628 = t27*t596
      t633 = t47*BSnThl
      t635 = t603*SnThr
      t636 = t635*sll
      t637 = t8*t636
      t640 = t46*Mzsq
      t641 = t27*t640
      t648 = t27*t375
      t651 = t48*BSnThl
      t652 = t27*t651
      t655 = t27*t404
      t656 = t603*sqq
      t657 = t93*t656
      t662 = 8*t27*t276*t605+8*t609*t16-16*t613*t616+8*t609*t174-16*t622
     #*t625-8*t628*t9+8*t613*t16+8*t27*t633*t637-48*t641*t616-16*t609*t6
     #16+8*t622*t605-48*t648*t625+8*t652*t637-8*t655*t657+24*t641*t174
      t663 = t37*t603
      t664 = t663*t32
      t665 = 2*t664
      t666 = t27*t191
      t668 = t404*t93*t47
      t671 = SnThr*sll
      t672 = t614*t671
      t675 = t2*Mzsq
      t676 = t675*t100
      t677 = t98*t103
      t678 = t603*t101
      t679 = t677*t678
      t682 = Mzsq*BSnThl
      t683 = t682*El
      t684 = t32*Eq
      t685 = t93*t684
      t688 = t404*t7
      t689 = t31*t19
      t696 = t19*t517
      t697 = t696*t515
      t698 = t520*t518
      t702 = t520*t603
      t703 = t518*sqq
      t704 = t702*t703
      t707 = t518*sll
      t711 = t19*Mzsq
      t712 = t711*t517
      t713 = t515*t520
      t714 = t713*t703
      t717 = t27*t682
      t718 = t8*t671
      t721 = t713*t707
      t724 = t603*t518
      t725 = t713*t724
      t728 = t671*t684
      t731 = -t665-16*t666*t668-8*t655*t672-192*t676*t679+8*t683*t685+96
     #*t688*t671*t689-48*t688*t635*t689-64*t697*t698*t191+32*t697*t704+3
     #2*t697*t702*t707+32*t712*t714-8*t717*t718+32*t712*t721-64*t712*t72
     #5-16*t688*t728
      t733 = t635*t684
      t736 = t711*BSnThl
      t737 = SnThr*t31
      t738 = t8*t737
      t741 = t116*t103
      t742 = t101*sll
      t749 = t2*t100
      t750 = t749*t98
      t751 = t103*t603
      t755 = t43*BSnThl
      t756 = t755*El
      t757 = t93*t191
      t764 = t43*Mzsq
      t765 = t764*BSnThl
      t766 = SnThr*sqq
      t767 = t8*t766
      t770 = t677*t742
      t779 = Mzsq*t100
      t780 = t779*t98
      t781 = t386*t33
      t786 = t8*t635
      t789 = 8*t688*t733-48*t736*t738-192*t741*t742*t35+96*t741*t678*t35
     #+96*t750*t751*t742-192*t756*t757+96*t756*t657+96*t756*t672+96*t765
     #*t767+96*t676*t770+48*t741*t742*t33-24*t741*t678*t33-24*t780*t781+
     #96*t765*t718-192*t765*t786
      t790 = t57*BCsThl
      t791 = Mzsq*t603
      t792 = t791*t161
      t795 = t57*t417
      t796 = BCsThl*t603
      t800 = BCsThl*Mzsq
      t801 = t800*t656
      t804 = t603*sll
      t805 = t800*t804
      t808 = t97*BCsThl
      t811 = CsThr*Mzsq
      t812 = t811*t656
      t815 = t53*BCsThl
      t816 = t31*t603
      t817 = t225*t816
      t820 = t97*t46
      t821 = t608*t656
      t824 = t804*t58
      t827 = t53*t46
      t828 = t243*t656
      t831 = Mzsq*t37
      t832 = t831*t656
      t837 = t57*t159
      t838 = CsThr*t603
      t839 = t838*t191
      t844 = t791*t191
      t847 = t101*sqq
      t848 = t677*t847
      t851 = 8*t790*t792+2*t795*t796*t191+4*t795*t801-8*t795*t805-8*t808
     #*t792+24*t808*t812-4*t815*t817-48*t820*t821+24*t820*t824+4*t827*t8
     #28-8*t827*t832-12*t790*t812+2*t837*t839-4*t837*t812+2*t827*t844+96
     #*t676*t848
      t856 = El*t37
      t857 = t682*t856
      t858 = t93*t273
      t861 = t93*t403
      t864 = t779*t677
      t865 = t1*sqq
      t866 = t742*t865
      t869 = t80*t656
      t871 = t678*t865
      t874 = t682*t8
      t875 = t671*t403
      t878 = t80*t804
      t880 = t635*t403
      t884 = t614*t766*Eq
      t887 = t711*t404
      t890 = t404*t177
      t896 = t116*t751
      t899 = t12*t404
      t900 = SnThr*t152
      t904 = 16*t717*t786+8*t857*t858+8*t857*t861-24*t864*t866-8*t869+48
     #*t864*t871+8*t874*t875-8*t878-16*t874*t880-16*t857*t884-48*t887*t7
     #57+8*t890*t880+8*t890*t635*t273-24*t896*t866-4*t899*t614*t900
      t905 = t31**2
      t906 = sll*t905
      t909 = t46*BSnThl
      t910 = t12*t909
      t911 = t900*sqq
      t912 = t8*t911
      t916 = t93*t31
      t917 = t404*t916
      t922 = Mzsq*t905
      t925 = t19*BSnThl
      t926 = t925*t8
      t931 = t404*t614
      t934 = t101*t31
      t935 = t76*BCsThl
      t936 = t934*t935
      t939 = t152*sqq
      t940 = t93*t939
      t943 = Mzsq*t32
      t944 = t943*t1
      t945 = 8*t944
      t946 = 2*t906
      t947 = t635*t152
      t948 = t8*t947
      t951 = t93*t816
      t954 = t12*t640
      t957 = -2*t906*t4-8*t910*t912-16*t48*t711*t917-16*t890*t875-t922+9
     #6*t887*t657-48*t926*t635*t191+t922*t4+8*t931*t875+8*t896*t936+8*t8
     #99*t940+t945+t946+4*t910*t948-4*t899*t951-4*t954*t917
      t959 = t603*t905
      t960 = t671*t31
      t961 = t8*t960
      t964 = Mzsq*t31
      t965 = t12*t964
      t968 = t12*t682
      t969 = t8*t900
      t972 = t640*BSnThl
      t973 = t12*t972
      t974 = t635*sqq
      t975 = t8*t974
      t978 = t48*sll
      t979 = t3*t978
      t982 = t3*t621
      t987 = t909*El
      t988 = t12*t987
      t989 = t177*t636
      t994 = t423*El
      t995 = t12*t994
      t998 = t93*t152
      t999 = t404*t998
      t1002 = t7*t31
      t1004 = t6*t1002*t603
      t1007 = t53*t4
      t1012 = -t959+8*t910*t961-4*t965*t9-4*t968*t969+8*t973*t975-64*t97
     #9*t616-64*t982*t605+8*t968*t975-16*t988*t989+16*t666*t9+32*t995*t6
     #25+4*t954*t999-4*t376*t1004-4*t1007*t828+8*t1007*t832
      t1015 = t93*t208
      t1018 = t3*t612
      t1031 = t4*t98
      t1032 = t97*t1031
      t1034 = t102*t118*t603
      t1041 = t4*El
      t1042 = t12*t1041
      t1045 = t386*t603
      t1046 = t116*t1045
      t1049 = t4*BSnThl
      t1050 = t12*t1049
      t1054 = t116*t386*Mzsq
      t1057 = t3*t1041
      t1060 = t4*Mzsq
      t1061 = t12*t1060
      t1064 = 8*t968*t637+8*t899*t1015+128*t1018*t616-16*t995*t605-64*t1
     #018*t174-8*t973*t637+128*t979*t16-64*t1018*t16+32*t1032*t1034-16*t
     #12*t640*t37*t16-16*t1042*t1004+32*t313*t1046+32*t1050*t961+32*t313
     #*t1054+64*t1057*t605-16*t1061*t917
      t1068 = t3*t1060
      t1071 = t3*t505
      t1074 = t386*sqq
      t1075 = t116*t1074
      t1078 = t688*t974
      t1082 = t191*El*t15
      t1085 = t97*t1060
      t1088 = t8*t603
      t1089 = t671*sqq
      t1090 = t1088*t1089
      t1093 = t32*t603
      t1094 = t53*t1093
      t1095 = t603*t152
      t1098 = t97*t816
      t1101 = t80*t791
      t1103 = t57*t816
      t1106 = 24*t675*t191
      t1108 = 8*t831*t33
      t1109 = 64*t1068*t16-128*t1071*t16-64*t313*t1075+32*t1061*t1078-16
     #*t1061*t1082+32*t1085*t1075-16*t1050*t1090+t1094-4*t97*t1095-4*t10
     #98+t959*t4+16*t1101+2*t1103-t1106+t1108
      t1110 = t4*t37
      t1111 = t804*t865
      t1114 = t53*t47
      t1115 = t46*t603
      t1119 = t640*t816
      t1122 = t97*t417
      t1123 = BCsThl*t37
      t1129 = t1123*t804
      t1132 = t1123*t791
      t1137 = t159*t603
      t1138 = t1137*t191
      t1141 = t159*Mzsq
      t1142 = t1141*t656
      t1147 = t38*t262
      t1148 = 16*t1147
      t1149 = t582*t33
      t1150 = 8*t1149
      t1151 = t791*t33
      t1152 = 16*t1151
      t1155 = t57*t46
      t1158 = -8*t1110*t1111-2*t1114*t1115*t208+4*t1114*t1119+16*t1122*t
     #1123*t656-8*t1122*t801+16*t1122*t1129-32*t1122*t1132+8*t1122*t805-
     #8*t1122*t1138+16*t1122*t1142-8*t808*t839-t1148+t1150-t1152+2*t790*
     #t839-2*t1155*t824
      t1166 = t663*t865
      t1170 = 48*t675*t656
      t1171 = t663*t33
      t1172 = 8*t1171
      t1173 = t811*t191
      t1176 = t243*t596
      t1179 = t191*t608
      t1182 = t225*t208
      t1185 = t811*t38
      t1188 = t831*t816
      t1189 = 4*t1188
      t1190 = t159*sll
      t1191 = t1190*t596
      t1194 = t97*Mzsq
      t1195 = t4*t191
      t1198 = t1123*t582
      t1201 = 4*t1155*t821+2*t795*t1138-4*t795*t1142+16*t1060*t1166+t117
     #0+t1172+2*t837*t1173+4*t827*t1176+24*t820*t1179+8*t815*t1182-4*t81
     #5*t1185+t1189-8*t1122*t1191-8*t1194*t1195+16*t1122*t1198
      t1208 = t640*t208
      t1211 = BCsThl*t152
      t1212 = t1211*sqq
      t1215 = t75*sll
      t1218 = t75*Mzsq
      t1221 = t246*Mzsq
      t1226 = t80*t191
      t1228 = t80*t582
      t1230 = t939*CsThr
      t1233 = Mzsq*t152
      t1234 = t1233*CsThr
      t1237 = t800*t152
      t1240 = t57*t208
      t1242 = t57*t964
      t1244 = t80*t596
      t1246 = -32*t1122*t1123*t191+16*t1122*t800*t54-2*t1114*t1208-8*t79
     #5*t1212+4*t795*t1215-2*t795*t1218-2*t1155*t1221+4*t1155*t338+16*t1
     #226-8*t1228+8*t790*t1230-4*t790*t1234+4*t795*t1237-4*t1240+2*t1242
     #-8*t1244
      t1249 = t53*t943
      t1250 = t243*t865
      t1257 = t57*Mzsq
      t1260 = t800*t191
      t1263 = t80*t417
      t1264 = t1141*sll
      t1267 = t1190*sqq
      t1270 = t504*sqq
      t1273 = t608*sqq
      t1276 = t608*sll
      t1279 = t191*t47
      t1286 = t44*t417
      t1289 = t1249-8*t1060*t1250+2*t795*t1191-2*t1155*t1179+2*t1257*t11
     #95+2*t795*t1260-8*t1263*t1264+16*t1263*t1267+48*t1263*t1270+8*t555
     #*t1273+8*t555*t1276-16*t555*t1279+8*t97*t939+4*t808*t1234-64*t1286
     #*t1267
      t1292 = t1141*sqq
      t1295 = t44*t46
      t1298 = t1233*sqq
      t1300 = t640*t32
      t1303 = t97*t208
      t1305 = t234*t32
      t1308 = t220*t26
      t1311 = t811*sll
      t1314 = t800*sqq
      t1317 = t800*sll
      t1320 = t811*sqq
      t1325 = t97*t1233
      t1327 = t48*t2
      t1328 = t191*Mzsq
      t1331 = t208*t1
      t1334 = 32*t1286*t1264+32*t1286*t1292+128*t1295*t1279-t827*t1298-2
     #*t1114*t1300+8*t1303+4*t1114*t1305-16*t1308*t298+8*t1308*t1311-24*
     #t1263*t1314-24*t1263*t1317+8*t1308*t1320-8*t1263*t1292-4*t1325-40*
     #t1327*t1328+16*t1110*t1331
      t1336 = t44*Mzsq
      t1343 = t43*BCsThl
      t1344 = t1343*Er
      t1347 = t161*t31
      t1351 = 2*t663*t208
      t1352 = t41*BCsThl
      t1355 = t243*t31
      t1358 = t53*Mzsq
      t1359 = t38*t46
      t1364 = t38*t1
      t1369 = t804*t33
      t1370 = 8*t1369
      t1371 = t159*t31
      t1372 = t1371*sll
      t1375 = t1371*Mzsq
      t1378 = -64*t1336*t978+32*t1336*t505-8*t808*t1230-64*t1344*t298+24
     #*t808*t1347-t1351-16*t1352*t1311+2*t1355*t1060+4*t1358*t1359-12*t1
     #060*t1331-8*t1060*t1364+56*t473*t1328+t1370-4*t795*t1372+2*t795*t1
     #375
      t1379 = t41*t417
      t1390 = t2*t603
      t1392 = 24*t1390*t191
      t1395 = t80*CsThr
      t1405 = t38*t582
      t1406 = 2*t1405
      t1413 = t12*t621
      t1416 = -32*t1379*t1270+16*t1379*t1314+16*t1379*t1317+32*t1352*t29
     #8-16*t1352*t1320-t1392-6*t1155*t1328-8*t1395*t1292-8*t1395*t1264-1
     #2*t790*t1347+16*t80*t159*t298-t1406+24*t820*t1221-48*t820*t338-8*t
     #827*t1355+8*t1413*t1004
      t1421 = t12*t612
      t1424 = t12*t651
      t1433 = sll*t32
      t1434 = t53*t1433
      t1436 = t152*t31
      t1445 = Mzsq*Eq
      t1446 = t1445*t688
      t1447 = t671*t125
      t1450 = t97*t964
      t1452 = t48*t92
      t1453 = t614*t1089
      t1456 = t635*t125
      t1459 = Eq*BSnThl
      t1460 = t1459*t8
      t1464 = 8*t1421*t917-16*t1424*t961+64*t1071*t616-128*t1068*t616-64
     #*t1085*t1046-2*t1434-2*t53*t1436+8*t1421*t1082+8*t1424*t1090-16*t1
     #421*t1078-4*t1446*t1447-4*t1450-16*t1452*t1453+8*t1446*t1456-4*t14
     #60*t635*t545
      t1470 = t48*t755
      t1475 = t8*t1089
      t1478 = t48*t764
      t1490 = t4*t711
      t1493 = t93*sll
      t1494 = t856*t1493
      t1497 = t4*t116
      t1502 = t4*t755
      t1509 = -8*t12*t582*t16+6*t1257*t125+96*t1470*t637+96*t1470*t975-1
     #92*t1470*t1475+96*t1478*t174-192*t1478*t616+96*t1478*t16-96*t4*t74
     #9*t677*t678*sll+64*t1490*t917-16*t973*t1494-96*t1497*t751*t847*t2+
     #384*t1502*t1475-192*t1502*t975+64*t1068*t174
      t1511 = t206*sqq
      t1515 = t4*t764
      t1520 = t4*t675
      t1521 = t386*sll
      t1522 = t116*t1521
      t1533 = t48*t736
      t1536 = t1459*t1088
      t1541 = t4*t1331
      t1544 = t48*t404
      t1546 = t93*t208*t19
      t1553 = 192*t1497*t386*t1511-192*t1515*t16-192*t1515*t174-96*t1520
     #*t1522-192*t1502*t637+384*t1515*t616-96*t1520*t1075+8*t864*t936+32
     #*t1533*t975-4*t1536*t1447+8*t1054*t299-16*t387*t1541+32*t1544*t154
     #6-16*t1533*t1475+8*t1046*t299
      t1557 = t677*t603
      t1559 = t847*t935
      t1563 = t614*t737*t19
      t1567 = t93*t545
      t1572 = t234*t31
      t1575 = t640*t152
      t1578 = t46*t152
      t1579 = t1578*sqq
      t1582 = t640*t31
      t1585 = t1578*t31
      t1589 = t582*t31
      t1591 = 2*t53*t1589
      t1592 = t370*t603
      t1593 = t1592*sqq
      t1596 = t1592*sll
      t1599 = t1115*t32
      t1603 = 8*t1460*t671*t545-16*t779*t1557*t1559-16*t1544*t1563-4*t14
     #45*t404*t1567-8*t899*t1453+8*t97*t1572+4*t97*t1575-8*t97*t1579-4*t
     #97*t1582+2*t53*t1585+t53*t1298+t1591+96*t48*t1593+96*t48*t1596+t53
     #*t1599-2*t791*t545
      t1609 = t191*t603
      t1611 = 8*t97*t1609
      t1612 = t445*t603
      t1615 = t1115*t152
      t1618 = t57*t1609
      t1620 = t791*sqq
      t1621 = t57*t1620
      t1623 = t1115*sll
      t1626 = t640*t603
      t1629 = t1115*sqq
      t1632 = t608*t603
      t1633 = t80*t1632
      t1635 = t47*t603
      t1636 = t1635*sll
      t1637 = t80*t1636
      t1639 = t1635*sqq
      t1640 = t80*t1639
      t1642 = t370*Mzsq
      t1643 = t1642*t603
      t1646 = t816*t2
      t1649 = t804*t545+2*t663*t545-t1611-4*t97*t1612+4*t97*t1615+2*t161
     #8-4*t1621+24*t80*t1623-48*t80*t1626+24*t80*t1629-16*t1633+8*t1637+
     #8*t1640+192*t4*t1643+56*t4*t1646
      t1652 = t1093*t1
      t1667 = t804*t31
      t1669 = 2*t53*t1667
      t1670 = t1095*sqq
      t1672 = t791*t31
      t1673 = t53*t1672
      t1675 = t97*t1620
      t1677 = t804*Mzsq
      t1686 = -96*t4*t1596-12*t4*t1652-96*t4*t1593-32*t418*t1596+64*t418
     #*t1643-32*t418*t1593-192*t48*t1643-40*t48*t1646+t1669+t53*t1670-2*
     #t1673+8*t1675+8*t97*t1677+8*t418*t1646+4*t48*t1652-6*t57*t1612
      t1688 = t246*t603
      t1689 = t57*t1688
      t1691 = t796*t152
      t1694 = t75*t603
      t1697 = t57*t47
      t1698 = t1697*t1620
      t1705 = 8*t97*t1328
      t1707 = t37*t32
      t1712 = t582*t545
      t1713 = t370*sll
      t1714 = t1713*sqq
      t1717 = t1642*sll
      t1720 = t1433*t1
      t1723 = t208*t2
      t1726 = -2*t1689+4*t795*t1691-2*t795*t1694+4*t1698+12*t1155*t1620-
     #2*t1155*t1688-t1705+t53*t1300+2*t1707*t1060-4*t1707*t505+t1712+64*
     #t418*t1714-32*t418*t1717+24*t4*t1720-112*t4*t1723
      t1729 = t675*t31
      t1732 = t1642*sqq
      t1737 = t54*t603
      t1754 = t816*t1
      t1757 = t48*t603
      t1760 = t418*t2
      t1763 = t1371*t603
      t1766 = 192*t4*t1714+56*t4*t1729-96*t4*t1732-96*t4*t1717-16*t820*t
     #1737-40*t48*t1729-12*t4*t944-192*t48*t1714+96*t48*t1717+80*t48*t17
     #23+96*t48*t1732-2*t53*t1305-8*t612*t1754+4*t1757*t1331+8*t1760*t16
     #09-8*t1122*t1763
      t1774 = t663*Mzsq
      t1781 = t663*t31
      t1784 = t1697*t1609
      t1788 = CsThr*t31
      t1789 = t1788*t603
      t1794 = 64*t1714
      t1798 = 8*t831*t350
      t1799 = t1*t603
      t1800 = t1799*sqq
      t1802 = 16*t831*t1800
      t1804 = 8*t663*t350
      t1805 = 4*t1122*t1694-4*t1122*t1691+32*t820*t1774+24*t820*t1688-2*
     #t827*t1672+4*t827*t1781-2*t1784-6*t1155*t1609+2*t837*t1789-16*t176
     #0*t1620-t1794-16*t418*t1723+t1798-t1802+t1804
      t1806 = t57*t338
      t1808 = t57*t1328
      t1814 = t57*t1221
      t1818 = t640*sqq
      t1821 = t640*sll
      t1824 = t80*t1273
      t1826 = t80*t1276
      t1828 = t80*t1279
      t1836 = t4*t603
      t1841 = 4*t1806+2*t1808+12*t57*t1572-6*t57*t1582-2*t1814-48*t80*t4
     #90+24*t80*t1818+24*t80*t1821+8*t1824+8*t1826-16*t1828-32*t418*t173
     #2+24*t1060*t1754+56*t473*t1609-12*t1836*t1331-8*t1110*t1754
      t1847 = t663*sll
      t1856 = t838*sll
      t1859 = t838*sqq
      t1862 = t800*t603
      t1865 = t811*t603
      t1874 = t1137*sqq
      t1877 = t1137*sll
      t1880 = -8*t820*t1677+8*t820*t1620-16*t820*t1847+4*t48*t944+8*t418
     #*t1729-8*t48*t1720-16*t1352*t1856-16*t1352*t1859-32*t1379*t1862+32
     #*t1352*t1865-4*t1774*t125+2*t1847*t125-112*t473*t1620-8*t1395*t187
     #4-8*t1395*t1877
      t1881 = t1141*t603
      t1888 = t796*sqq
      t1891 = t796*sll
      t1908 = t838*t152
      t1917 = 16*t1395*t1881+6*t790*t1789+2*t795*t1763+16*t1379*t1888+16
     #*t1379*t1891-16*t1308*t1865-8*t1263*t1877+16*t1263*t1881-24*t1263*
     #t1888-24*t1263*t1891-16*t555*t1632+8*t555*t1639+4*t808*t1908-64*t1
     #295*t1639+128*t1295*t1632-64*t1295*t1636
      t1949 = 32*t1286*t1877+32*t1286*t1874-64*t1286*t1881-2*t1358*t1609
     #-t827*t1670-2*t1114*t1599+8*t1308*t1856+8*t1308*t1859+48*t1263*t18
     #62-8*t1263*t1874-64*t1336*t128+80*t1327*t1620-40*t1327*t1609-64*t1
     #344*t1865+32*t1344*t1856
      t1958 = t54*Mzsq
      t1961 = t243*Mzsq
      t1966 = t439*t31
      t1971 = t1697*t1328
      t1975 = t1788*Mzsq
      t1986 = -12*t808*t1789+32*t1344*t1859+8*t1122*t1212+32*t820*t352-1
     #6*t820*t1958-16*t820*t1961-12*t1194*t125+4*t1194*t1966+32*t1336*t1
     #22-2*t1971-4*t837*t1347+2*t837*t1975+8*t1760*t1328+4*t612*t1331-8*
     #t1122*t1375+16*t1122*t1372
      t1994 = t582*t32
      t1995 = 16*t1720
      t1996 = 32*t1596
      t1997 = 32*t1593
      t1998 = 32*t1717
      t1999 = 8*t1652
      t2000 = t1707*Mzsq
      t2001 = 2*t2000
      t2002 = t804*t32
      t2003 = t791*t32
      t2004 = 2*t2003
      t2005 = 24*t1729
      t2006 = 24*t1646
      t2007 = 64*t1643
      t2008 = -8*t1122*t1215-4*t1122*t1237-4*t1007*t1176-t1994-t1995+t19
     #96+t1997+t1998+t1999-t2001-t2002+t2004-t2005-t2006-t2007
      t2009 = 32*t1732
      t2010 = 48*t1723
      t2011 = t1707*sll
      t2012 = 4*t2011
      t2013 = t4*t736
      t2018 = t4*t404
      t2025 = t4*t92
      t2037 = t856*t615
      t2044 = t2009+t2010+t2012+64*t2013*t1475-128*t2013*t975+64*t2018*t
     #1563+8*t613*t174-128*t2018*t1546+64*t2025*t1453+24*t648*t605+192*t
     #1520*t1046+24*t641*t16+24*t27*t909*t637+32*t973*t2037+8*t555*t1636
     #-4*t790*t1908
      t2052 = 1/(-sqq+t603)
      t2054 = 1/(-sll+t603)
      t2057 = t594*t580
      t2059 = t149*t156*t572
      t2060 = t2057*t2059
      A1e0c6p1t1n = -24*(t1109+t1603+t1553+t1509+t1464+t957+t1064+t1805+
     #t1766+t1726+t1686+t1158+t1649+t1246+t2008+t1986+t1949+t1917+t1880+
     #t1841+t1416+t1378+t2044+t1334+t1289+t662+t731+t904+t851+t1201+t789
     #+t1012)*t151*t2052*t2054*t2060
      t2063 = Mzsq**2
      t2064 = t31*t2063
      t2066 = 6*t53*t2064
      t2067 = t32*t2063
      t2068 = t97*t2063
      t2074 = t404*SnThr
      t2075 = Eq*t2063
      t2079 = t1*t2063
      t2083 = sqq*t2063
      t2087 = t35*t2063
      t2088 = 8*t2087
      t2089 = t2063*El
      t2093 = t2089*t15
      t2097 = t6*t22*t2063
      t2100 = t4*t2063
      t2101 = t12*t2100
      t2108 = t2066-t2067+32*t2068*t387+64*t3*t2063*t9+4*t2074*t1002*t20
     #75-8*t117*t118*t2079-16*t92*t93*t2083-t2088-16*t12*t2089*t23-64*t5
     #*t2093+16*t21*t2097-8*t2101*t16+32*t28*t2093-16*t1155*t1961
      t2113 = t48*t1
      t2115 = t26*t2063*sqq
      t2118 = t297*t2063
      t2123 = t1788*t2063
      t2128 = t53*t37
      t2129 = t2128*t1328
      t2130 = 16*t2129
      t2131 = 4*t1405
      t2132 = t41*t243
      t2134 = t41*t831
      t2136 = t37*t2063
      t2137 = t57*t2136
      t2139 = t57*t2083
      t2141 = t37*t152
      t2142 = t57*t2141
      t2147 = 16*t1155*t1958-8*t1155*t352-2*t2113*t2115+8*t808*t2118+8*t
     #2068*t128+2*t815*t2123-8*t790*t2118+t2130-t2131+32*t2132-64*t2134+
     #8*t2137-10*t2139+8*t2142+8*t63*t40*t2063
      t2149 = t58*t2063
      t2150 = t57*t2149
      t2152 = t54*t2063
      t2154 = 8*t53*t2152
      t2157 = t2063*t46
      t2158 = t2157*t47
      t2161 = t2141*sqq
      t2163 = 8*t53*t2161
      t2164 = t4*t2087
      t2166 = t33*t2063
      t2169 = t53*t1355
      t2170 = 16*t2169
      t2171 = t831*t31
      t2172 = t53*t2171
      t2173 = 32*t2172
      t2180 = t40*t37
      t2187 = t41*t2063
      t2190 = t68*t2063
      t2193 = -2*t2150-t2154-64*t44*t2100+32*t44*t2158-t2163+16*t2164-4*
     #t4*t2166+t2170-t2173+8*t71*t72*t152+8*t71*t72*t2063+32*t63*t2180*s
     #ll-64*t63*t2180*Mzsq+8*t48*t2187-8*t67*t2190
      t2194 = t57*t1961
      t2196 = t57*t1958
      t2198 = t57*t352
      t2206 = t47*t2063
      t2207 = t41*t2206
      t2209 = t38*t2063
      t2210 = 2*t2209
      t2211 = t38*t152
      t2212 = 2*t2211
      t2213 = t44*t2063
      t2214 = 32*t2213
      t2215 = 4*t2000
      t2216 = 4*t2166
      t2217 = 2*t2011
      t2219 = t97*t2083
      t2220 = 16*t2219
      t2221 = -16*t2194+80*t2196-40*t2198-2*t57*t2190+32*t80*t2100-32*t8
     #0*t2158+8*t2207+t2210+t2212+t2214+t2215+t2216-t2217+8*t2187-t2220
      t2226 = t580*t586
      t2227 = 1/Mzsq
      t2228 = t2227*t156
      A1e0c7p1t1n = -24*(t2108+t2147+t2193+t2221)*t149*t151*t2226*t2228
      t2232 = t37*t47
      t2235 = t11-t30-t164-t169+t176-t180+t184-t186-t190-t195+t199+16*t5
     #7*t2232+t200+t42-t45+t202-t205-t207-t208+t213
      t2238 = -t245+t248-t251+t253-t255+t257-t50+t52+16*t73-t259+t82
      t2245 = -t286+t289+t292+t296-t300+t302+t305+t309-t312-t314+t320
      t2248 = 16*t139
      t2249 = -t342+t346+t348-t349+t132+t351+t353-t2248-t354+t144-t355
      A1e0c8p1t1n = 24*Mzsq*(-t110-t106-t91+t89-t84+t279+t283+t275+t269-
     #t272+t267-t236+t239-t242+t233-t227-t230-t219+t222-t224+t216-16*t57
     #*t423*t47+t2245+t2249+t2238+t2235+t337+t338-t340-t335-t333+t331+t3
     #27-t325)*t362*t600
      t2256 = t663*t939
      t2264 = t608*t816
      t2267 = t2206*t656
      t2270 = t225*t1095
      t2273 = CsThr*t2063
      t2274 = t2273*t656
      t2277 = t811*t816
      t2280 = t1697*t844
      t2286 = t1137*t582
      t2293 = t2136*t656
      t2300 = t2063*t603
      t2306 = -4*t1007*t2256-16*t820*t832+12*t815*t225*t1093+24*t820*t22
     #64-48*t820*t2267+16*t808*t2270-12*t790*t2274+6*t790*t2277-2*t2280-
     #6*t1155*t844-8*t1395*t1142-8*t1395*t2286+8*t1760*t844+4*t1757*t114
     #9-8*t827*t2293+4*t827*t1188-t827*t1095*t596+2*t827*t2300*t191-12*t
     #808*t2277
      t2311 = t838*t1233
      t2316 = t2273*t804
      t2319 = t640*t1093
      t2322 = t2157*t816
      t2325 = t1141*t816
      t2330 = t57*CsThr
      t2331 = t159*t2063
      t2332 = t2331*t656
      t2337 = t663*t582
      t2346 = t1635*t582
      t2349 = t2273*t663
      t2352 = t796*t1233
      t2357 = 24*t808*t2274+32*t1344*t812+4*t808*t2311+32*t1344*t792-8*t
     #808*t2316-2*t1114*t2319+4*t1114*t2322-8*t1122*t2325-48*t808*t817-4
     #*t2330*t2332+2*t2330*t2325-16*t820*t2337+32*t1286*t1142+32*t1286*t
     #2286-64*t1295*t821-64*t1295*t2346+64*t808*t2349-4*t1122*t2352-8*t1
     #263*t2286
      t2365 = BCsThl*t2063
      t2366 = t2365*t804
      t2373 = t2365*t663
      t2378 = t800*t816
      t2381 = t2365*t656
      t2388 = t225*t804
      t2391 = t225*t791
      t2405 = -16*t790*t817-2*t1155*t2264+4*t1155*t2267-8*t795*t2366-4*t
     #790*t2311+16*t1122*t2332-32*t1122*t2373+8*t1122*t2366+4*t1122*t237
     #8-8*t1122*t2381-8*t1110*t1151-24*t1263*t805+64*t1308*t2388-128*t13
     #08*t2391+64*t1308*t225*t656+8*t1308*t792-8*t1263*t1142+8*t1308*t81
     #2-24*t1263*t801
      t2414 = t2128*t844
      t2443 = t1095*t865
      t2450 = 8*t555*t821+32*t795*t2373+4*t795*t2381-2*t795*t2378-16*t24
     #14-4*t795*t2332+2*t795*t2325+8*t790*t2316-32*t790*t2349+12*t1007*t
     #663*t1328-2*t53*t48*t804*t964+8*t555*t2346+16*t1379*t805-16*t1352*
     #t792-16*t1352*t812+16*t1379*t801-16*t2100*t1166-8*t1110*t2443-40*t
     #1327*t844-12*t1836*t1149
      t2455 = t943*t4
      t2463 = t225*t1233
      t2469 = sll*t2063
      t2470 = t58*t2469
      t2481 = t297*t2469
      t2488 = t54*t2469
      t2495 = t811*t208
      t2502 = 56*t473*t844+12*t2128*t2455+4*t815*t225*t2064+96*t808*t118
     #2+16*t808*t2463-32*t808*t225*t939+24*t820*t2470-48*t808*t1185+8*t8
     #15*t225*t1436-24*t2128*t1433*t4-8*t808*t2481+32*t820*t1176-8*t827*
     #t1405+4*t827*t2488+128*t1295*t1179-64*t1344*t1173+24*t808*t2495+16
     #*t1395*t1191-12*t790*t2495
      t2503 = t2331*t191
      t2506 = t1141*t208
      t2509 = t97*t47
      t2516 = t800*t208
      t2519 = BCsThl*sqq
      t2520 = t2519*t1233
      t2523 = t2365*t243
      t2526 = t2365*t54
      t2537 = t2273*t243
      t2540 = t2273*t54
      t2543 = t297*t1233
      t2549 = t417*Mzsq
      t2560 = 2*t2330*t2503-4*t2330*t2506-48*t2509*t1208+16*t1263*t1191-
     #4*t795*t2506-8*t1122*t2516+8*t1122*t2520+16*t1122*t2523+16*t1122*t
     #2526-8*t1122*t2503+4*t1358*t1433*t48-2*t1114*t445*t2469-32*t808*t2
     #537-32*t808*t2540-8*t808*t2543-64*t1286*t1191+16*t97*t159*t2549*t2
     #08+2*t795*t2365*t191-16*t795*t2526+4*t795*t2516
      t2570 = t608*t208
      t2604 = -8*t795*t2520-16*t795*t2523+2*t790*t2481-2*t1155*t2470+4*t
     #1155*t2570+32*t790*t1182-16*t790*t1185-128*t1308*t225*t191+64*t130
     #8*t225*t582+8*t790*t2543+16*t790*t2540+16*t790*t2537+2*t795*t2503-
     #32*t1379*t1260+32*t1352*t1173-16*t555*t1179+48*t1263*t1260-16*t130
     #8*t1173+64*t1308*t811*t54
      t2609 = Mzsq*t1
      t2613 = t19*t2063
      t2614 = t2613*t517
      t2617 = t2063*BSnThl
      t2618 = t2617*El
      t2621 = t2063*t100
      t2622 = t2621*t98
      t2625 = t2*t2063
      t2626 = t2625*t100
      t2631 = t2613*BSnThl
      t2640 = t43*t2063
      t2641 = t2640*BSnThl
      t2648 = t97*t439
      t2662 = 8*t2100*t1250-16*t1110*t1149-8*t1110*t939*t2609+32*t2614*t
     #714+8*t2618*t685-24*t2622*t781+96*t2626*t848+96*t2626*t770-48*t263
     #1*t738-192*t2626*t679-64*t2614*t725+32*t2614*t721+96*t2641*t718-19
     #2*t2641*t786+96*t2641*t767+16*t2648*t2337+16*t2648*t832-64*t385*t2
     #337+2*t57*t196*t844+8*t53*t663*t208*t4
      t2673 = t177*t603
      t2675 = sqq*Eq
      t2676 = t2675*Mzsq
      t2680 = t2617*t856
      t2683 = t403*Mzsq
      t2684 = t671*t2683
      t2689 = t755*t8
      t2705 = t2617*t8
      t2714 = t635*t1328
      t2720 = t2621*t677
      t2723 = -8*t1007*t2293+24*t97*t48*t844-32*t385*t828-24*t404*t2673*
     #t671*t2676-8*t2680*t861+16*t890*t2684+8*t931*t2684+96*t2689*t635*t
     #582+96*t2689*t635*t596+48*t387*t1149+16*t2680*t884-8*t2680*t858+8*
     #t890*t900*t2676+8*t2705*t875-16*t899*t177*t947+4*t12*t2157*t999-48
     #*t926*t2714-192*t2689*t671*t596+8*t2720*t936
      t2729 = t900*Mzsq
      t2735 = t2075*t688
      t2738 = t4*t2640
      t2741 = t57*t4
      t2744 = t57*t439
      t2751 = t804*t943
      t2752 = t27*t2617
      t2755 = t749*t677
      t2761 = t2300*sqq
      t2766 = t2331*t603
      t2769 = t2365*t603
      t2774 = t97*t603
      t2776 = 8*t2774*t1328
      t2777 = -4*t2075*t404*t1567+64*t385*t832-16*t899*t177*t2729-16*t27
     #05*t880-4*t2735*t1447-192*t2738*t16-48*t2741*t1176+32*t2744*t1176-
     #32*t2648*t1176+32*t1344*t2118-t2751-8*t2752*t718+96*t2755*t678*t59
     #6-96*t546*t1677-112*t473*t2761+80*t1327*t2761+16*t1395*t2766-32*t1
     #379*t2769-40*t1327*t1672-t2776
      t2780 = 16*t139*t1609
      t2781 = t2136*t603
      t2785 = t663*t152
      t2805 = t2300*t31
      t2810 = t2300*sll
      t2819 = t2780+4*t2781*t125+t2002*t1060+2*t2785*t125+2*t664*t1060-3
     #2*t550*t1620-32*t550*t1677+8*t1110*t1652+24*t2100*t1754-12*t1836*t
     #944+56*t473*t1672-96*t546*t1620+t827*t2003-2*t827*t2805-4*t820*t16
     #72-8*t820*t2810+8*t820*t2761+96*t538*t1677+96*t538*t1620
      t2830 = t1095*Mzsq
      t2834 = 8*t2128*t1667
      t2835 = t53*t2063
      t2836 = t2835*t1737
      t2838 = t1358*t1781
      t2842 = t2273*t603
      t2851 = t2206*t603
      t2854 = t53*t603
      t2856 = t328*Er
      t2858 = 32*t2856*t1609
      t2859 = t1194*t1737
      t2865 = 4*t1757*t944+8*t1760*t1672-16*t1760*t2761-8*t2158*t1754+32
     #*t820*t2781+4*t820*t2830-t2834+16*t2836+16*t2838-2*t2835*t1609+32*
     #t1352*t2842+16*t1263*t2766-16*t1308*t2842-2*t827*t1994+128*t1295*t
     #2851+t2854*t1298+t2858-64*t2859-2*t1355*t2100+4*t2011*t1060
      t2877 = t57*t603*t1328
      t2882 = 4*t2128*t1670
      t2883 = t2856*t1677
      t2888 = 2*t2854*t1589
      t2893 = t2469*sqq
      t2909 = 2*t2211*t1060+48*t1263*t2769+24*t555*t1620+24*t555*t1677+2
     #*t2877+32*t385*t1176+t2882+64*t2883-16*t555*t2851+t2888-64*t1344*t
     #2842-64*t1286*t2766+56*t473*t2893+192*t546*t1328+16*t1110*t1436*t1
     #-112*t473*t1589-16*t1110*t1720+24*t505*t944+8*t2100*t1364
      t2922 = t1371*t2063
      t2928 = 32*t1257*t1737
      t2931 = t1697*t2761
      t2933 = t1697*t1672
      t2935 = t557*t1620
      t2937 = t557*t1677
      t2939 = t159*sqq
      t2940 = t2939*t2063
      t2943 = t1697*t2893
      t2947 = t2331*sll
      t2951 = 4*t2128*t1298
      t2952 = t2128*t2893
      t2954 = -12*t2100*t1331+8*t1110*t944+80*t1327*t1589-40*t1327*t2893
     #-192*t538*t1328+6*t790*t2123+2*t2330*t2922-6*t1155*t1672-t2928+12*
     #t1155*t2761+4*t2931-2*t2933+8*t2935+8*t2937-8*t1395*t2940-2*t2943-
     #6*t1155*t2893-8*t1395*t2947+t2951-8*t2952
      t2956 = t75*t2063
      t2959 = t2365*t152
      t2970 = t2609*sqq
      t2972 = 24*t1847*t2970
      t2979 = t2136*sll
      t2984 = t54*t1233
      t2989 = t161*t2063
      t2992 = t2273*t152
      t3001 = 4*t1122*t2956-4*t1122*t2959-16*t1760*t1589+8*t1760*t2893+4
     #*t2158*t1331-8*t978*t944-t2972-8*t820*t1298+8*t820*t1589-16*t820*t
     #2152-16*t820*t2979+4*t1007*t2488-4*t1007*t2984+4*t827*t2209+32*t13
     #44*t2989+4*t808*t2992+8*t2128*t1589-32*t1194*t352-64*t347*t2190
      t3006 = t46*t32
      t3007 = t3006*t2063
      t3014 = t47*sll
      t3015 = t3014*t2063
      t3018 = t2063*t152
      t3019 = t3018*sqq
      t3021 = t1436*Mzsq
      t3031 = SnThr*t32
      t3032 = t3031*t1445
      t3035 = t2136*t816
      t3036 = 4*t3035
      t3037 = t663*t943
      t3038 = 2*t3037
      t3040 = 48*t2625*t656
      t3041 = t32*t1
      t3043 = 8*t663*t3041
      t3045 = 8*t1707*t2609
      t3046 = t2469*t33
      t3047 = 8*t3046
      t3049 = 8*t2136*t33
      t3050 = t1433*t2609
      t3051 = 16*t3050
      t3052 = -8*t1122*t2922-12*t808*t2123-2*t1114*t3007+32*t1286*t2947+
     #32*t1286*t2940-64*t1295*t3015-t827*t3019+2*t827*t3021-24*t2720*t86
     #6-4*t899*t93*t2064-8*t385*t844+8*t931*t3032-t3036-t3038+t3040-t304
     #3-t3045+t3047-t3049-t3051
      t3058 = 2*t38*t1233
      t3060 = 16*t1707*t262
      t3062 = 2*t38*t2469
      t3063 = t663*t1436
      t3064 = 2*t3063
      t3066 = 24*t1390*t964
      t3067 = t2300*t33
      t3068 = 16*t3067
      t3070 = 4*t1707*t582
      t3072 = 48*t206*t964
      t3073 = t152*t1
      t3075 = 16*t38*t3073
      t3077 = 64*t1713*t596
      t3079 = 24*t2625*t191
      t3080 = t1592*t582
      t3081 = 32*t3080
      t3083 = 32*t1592*t596
      t3084 = t1093*t2609
      t3085 = 8*t3084
      t3086 = t80*t2469
      t3088 = t57*t2064
      t3090 = t97*t2064
      t3092 = t97*t3018
      t3094 = t2063*t905
      t3096 = -t3058+t3060+t3062-t3064-t3066-t3068-t3070+t3072-t3075-t30
     #77-t3079+t3081+t3083+t3085-8*t3086+2*t3088-4*t3090-4*t3092+t3094*t
     #4
      t3097 = t80*t2300
      t3100 = 16*t139*t1328
      t3101 = t557*t1328
      t3105 = t1697*t1589
      t3107 = t246*t2063
      t3120 = t2365*sll
      t3123 = t2519*t2063
      t3138 = 16*t3097+t3100-16*t3101+12*t1155*t1589+4*t3105-2*t1155*t31
     #07-4*t790*t2992+4*t795*t2959-2*t795*t2956+2*t795*t2922+8*t555*t214
     #9-24*t1263*t3120-24*t1263*t3123-48*t555*t1328+8*t1308*t2118-8*t126
     #3*t2947-8*t1263*t2940+8*t555*t3015+8*t1308*t2989
      t3140 = t12*t2617
      t3150 = t2613*t404
      t3153 = t1002*t2063
      t3154 = t6*t3153
      t3157 = t635*t31
      t3158 = t8*t3157
      t3172 = t856*t94
      t3185 = t27*t2206
      t3188 = t614*t2729
      t3193 = 8*t3140*t975+8*t3140*t637-16*t48*t2089*t6*t1002*t19-48*t31
     #50*t757-4*t376*t3154-4*t968*t3158+8*t968*t961+96*t926*t671*t964+16
     #*t1379*t3123-16*t1352*t2118-16*t1352*t2989+32*t3140*t3172+32*t3140
     #*t1494+32*t899*t177*t911+16*t1379*t3120-48*t926*t635*t964+8*t3185*
     #t174-4*t899*t3188-64*t3140*t2037
      t3195 = CsThr*BSnThl
      t3196 = t3195*El
      t3198 = t615*t1328
      t3201 = t37*t905
      t3206 = t93*t1589
      t3210 = t614*t737*Mzsq
      t3213 = t48*t756
      t3214 = t93*t1328
      t3217 = t671*Mzsq
      t3218 = t614*t3217
      t3221 = t4*t756
      t3226 = t766*Mzsq
      t3227 = t614*t3226
      t3236 = t27*t972
      t3239 = t4*t2631
      t3244 = t4*t750
      t3246 = t751*t847*Mzsq
      t3254 = t742*Mzsq
      t3260 = -16*t12*BCsThl*t3196*t3198-2*t3201*t1060+64*t2025*t3198-12
     #8*t2025*t3206+64*t2025*t3210-192*t3213*t3214+96*t3213*t3218-192*t3
     #221*t3218+384*t3221*t3214+96*t3213*t3227-192*t3221*t3227+8*t12*t47
     #*t987*t3198-48*t3236*t1475-128*t3239*t975-8*t3140*t1475-96*t3244*t
     #3246+192*t3244*t386*t1328+64*t3239*t1475-96*t3244*t751*t3254+8*t61
     #3*t1078
      t3266 = t1088*t3217
      t3269 = t27*t2158
      t3286 = t177*t1089
      t3291 = t177*t974
      t3297 = t93*t1298
      t3306 = t635*Mzsq
      t3313 = 24*t3236*t975+64*t3*t1049*t3266+8*t3269*t174+32*t2101*t107
     #8+16*t717*t1475-16*t613*t1082+8*t652*t3266-8*t717*t975-8*t655*t321
     #8+96*t3150*t657+128*t655*t3286-64*t655*t989-64*t655*t3291-192*t275
     #5*t742*t596+8*t899*t3297+48*t968*t856*t916-96*t899*t177*t960+128*t
     #655*t177*t3306-64*t655*t177*t3217
      t3322 = t696*t713
      t3332 = t27*t608*BSnThl
      t3337 = t688*t3157
      t3348 = t27*t2157
      t3355 = t370*t2063
      t3356 = t3355*t603
      t3360 = 128*t80*t1774
      t3361 = t57*t2761
      t3363 = t57*t1672
      t3365 = t80*t2851
      t3367 = 48*t899*t177*t3157+96*t2755*t678*t582-16*t3185*t616+32*t33
     #22*t724*t582+8*t3185*t16-64*t655*t177*t3226-16*t3332*t1475-16*t210
     #1*t1082-16*t1061*t3337-16*t3269*t616+8*t3332*t975+8*t27*t633*El*t3
     #218-48*t3348*t616+24*t3348*t16+24*t3348*t174-192*t48*t3356+t3360-4
     #*t3361+2*t3363-16*t3365
      t3369 = t2157*t603
      t3372 = t80*t1677
      t3375 = 64*t80*t1737
      t3376 = t3355*sqq
      t3379 = t3355*sll
      t3382 = t2625*t31
      t3387 = t2067*t1
      t3399 = 64*t80*t1958
      t3400 = t57*t2893
      t3403 = 64*t80*t1961
      t3404 = t80*t1328
      t3406 = t57*t1589
      t3408 = t80*t2149
      t3411 = 128*t80*t352
      t3412 = -48*t80*t3369-8*t3372-t3375+96*t48*t3376+96*t48*t3379-40*t
     #48*t3382-32*t418*t3379-12*t4*t3387-96*t4*t3376-96*t4*t3379+56*t4*t
     #3382-32*t418*t3376-t3399+2*t3400-t3403+16*t3404-4*t3406+8*t3408+t3
     #411
      t3413 = t57*t3107
      t3415 = t445*t2063
      t3418 = t80*t3015
      t3422 = t2157*sll
      t3425 = t97*t2761
      t3431 = t3*t2100
      t3439 = t12*t683
      t3444 = t404*t37
      t3445 = t12*t3444
      t3450 = t97*t2100
      t3453 = t12*t2158
      t3460 = t688*t960
      t3463 = -2*t3413-6*t57*t3415+8*t3418+24*t80*t2190+24*t80*t3422+8*t
     #3425+32*t1032*t102*t118*t2063+64*t3431*t174-128*t3431*t616-t3094+3
     #2*t3322*t724*t596-64*t3439*t3291-32*t3439*t3286+64*t3445*t3218+4*t
     #795*t2352+32*t3450*t1522+8*t3453*t1082-16*t3453*t1078+8*t3269*t16+
     #32*t1061*t3460
      t3473 = t678*sqq
      t3483 = t742*sqq
      t3499 = t80*t2083
      t3506 = 48*t97*t1781
      t3507 = t97*t1672
      t3511 = t905*CsThr
      t3512 = t3511*BCsThl
      t3515 = 48*t2720*t871-16*t1421*t3460+8*t1421*t3337+32*t1085*t741*t
     #3473+48*t2741*t832+64*t1068*t1078-20*t1007*t1188-64*t1085*t741*t34
     #83+24*t27*t987*t3218+128*t1018*t1082+32*t97*t4*t100*t1557*t3254-12
     #8*t1068*t1082-8*t3499+8*t97*t2810-4*t97*t2830+t3506-4*t3507-4*t214
     #1*t545+4*t243*t3512
      t3518 = t2469*t545
      t3523 = 8*t97*t2893
      t3524 = t97*t2781
      t3526 = t53*t2003
      t3528 = 16*t97*t2785
      t3530 = 12*t53*t664
      t3531 = t53*t2805
      t3534 = 16*t57*t1781
      t3541 = 64*t80*t1847
      t3542 = t80*t1620
      t3548 = t964*t1
      t3550 = 8*t804*t3548
      t3552 = 6*t663*t1589
      t3553 = -2*t1707*t2100+t3518-2*t906*t1060+t53*t3007-t3523-64*t3524
     #+t3526-t3528-t3530-2*t3531+t3534-2*t663*t3512-2*t2300*t545+t959*t1
     #060-t3541-8*t3542+192*t4*t3356+64*t418*t3356+t3550+t3552
      t3556 = 16*t2136*t1800
      t3558 = 8*t663*t3548
      t3560 = 24*t1390*t1328
      t3563 = 8*t663*t939*t1
      t3565 = 8*t2136*t350
      t3567 = 32*t57*t1355
      t3569 = 16*t57*t2171
      t3571 = 24*t53*t2011
      t3573 = 8*t53*t2211
      t3575 = 12*t53*t2000
      t3578 = t2157*t152
      t3586 = 8*t2141*t2970
      t3588 = 16*t243*t3548
      t3594 = t3556+t3558-t3560+t3563-t3565-t3567+t3569+t3571-t3573-t357
     #5+8*t97*t1298+4*t97*t3578-2*t53*t3021-4*t97*t3415+t3586+t3588+4*t4
     #8*t3387+8*t418*t3382+t53*t3019
      t3596 = 96*t97*t1355
      t3597 = t97*t2979
      t3599 = t97*t2152
      t3603 = 16*t97*t2141*Mzsq
      t3605 = 32*t97*t2161
      t3606 = t97*t1589
      t3609 = 48*t97*t2171
      t3610 = t208*t2063
      t3612 = 2*t53*t3610
      t3613 = t53*t1994
      t3615 = t53*t2209
      t3624 = t243*t905
      t3625 = 4*t3624
      t3628 = t3*t2158
      t3631 = t964*t4
      t3640 = -t3596+32*t3597+32*t3599-t3603+t3605+8*t3606+t3609+t3612-2
     #*t3613-8*t3615+32*t3445*t1453-64*t1018*t1078-64*t3*t651*t3266-t362
     #5+64*t550*t1328-64*t3628*t174-6*t1847*t3631+64*t1057*t2097-16*t104
     #2*t3154-64*t3450*t1046
      t3646 = t12*t2157*BSnThl
      t3661 = t4*t2625
      t3672 = t856*t7
      t3676 = t48*t2631
      t3689 = 128*t3628*t616+32*t3646*t2037+4*t988*t3188-16*t3646*t1494+
     #8*t3646*t975-8*t3646*t637-64*t982*t2097+8*t973*t961+192*t3661*t104
     #6+8*t1413*t3154-16*t3646*t3172-4*t973*t3158-8*t988*t3297-16*t973*t
     #3672*t974+32*t3676*t975-16*t2621*t1557*t1559-16*t1452*t3210+32*t14
     #52*t3206+32*t973*t3672*t1089
      t3693 = t4*t3548
      t3699 = t2621*t677*t101
      t3710 = t505*t2970
      t3718 = t1459*t108
      t3728 = t48*t2640
      t3733 = 32*t3379
      t3736 = t53*t2067
      t3741 = -16*t988*t2673*t3217+8*t1046*t3693-16*t1522*t3693+8*t3699*
     #t299-8*t97*t418*t844+8*t2735*t1456-4*t1536*t3031*t1060+8*t1046*t37
     #10-16*t3676*t1475-4*t1536*t960*t1060+8*t3718*t1712-8*t12*t688*t271
     #4-16*t1452*t3198+384*t2738*t616+96*t3728*t174+16*t2752*t786+t3733+
     #16*t890*t728+t3736+8*t890*t635*t939*Eq
      t3764 = 32*t3376
      t3765 = t1707*t2063
      t3766 = 2*t3765
      t3767 = 64*t3356
      t3768 = t906*Mzsq
      t3769 = 2*t3768
      t3770 = t1433*t1445
      t3780 = t2300*t32
      t3781 = 2*t3780
      t3782 = t3201*Mzsq
      t3783 = 2*t3782
      t3784 = t663*t905
      t3785 = 2*t3784
      t3786 = -24*t896*t742*t2970+8*t890*t635*t2683-64*t3322*t707*t596-1
     #6*t890*t900*t403-8*t890*t733-24*t896*t934*t2609+24*t820*t3107-8*t8
     #90*t3032+t3764+t3766-t3767+t3769-16*t9*t3770-96*t3661*t1075-96*t36
     #61*t1522+64*t4*t2613*t917+t3781+t3783+t3785
      t3787 = 24*t3382
      t3788 = t2141*t32
      t3789 = 4*t3788
      t3790 = 8*t3387
      t3808 = t2469*t32
      t3820 = t959*Mzsq
      t3827 = t582*t865
      t3830 = -t3787+t3789+t3790-16*t2744*t832-8*t27*t2089*t23-4*t3140*t
     #969-192*t3728*t616+96*t3728*t16-192*t2738*t174+2*t2741*t844+2*t274
     #4*t844-t3808-16*t2741*t828-2*t57*t48*t844+2*t57*t418*t844+16*t57*t
     #1123*t792-t3820-16*t2744*t2337+4*t53*t423*t844+24*t4*t663*t3827
      t3838 = t2227*t2052*t2054
      A1e0c9p1t1n = -24*(t2604+t2560+t2405+t2502+t3553+t3515+t2450+t3463
     #+t3412+t3367+t3741+t3689+t3640+t3594+t2723+t2777+t2819+t2865+t2909
     #+t3830+t3786+t3052+t3001+t2954+t2357+t3313+t2306+t3260+t3193+t3138
     #+t3096+t2662)*t151*t3838*t2060
      t3842 = 2*t1781
      t3846 = 64*t5*t616
      t3847 = t4*t925
      t3851 = 32*t3847*t1475
      t3855 = t404*t93*Mzsq
      t3857 = 64*t5*t3855
      t3861 = 16*t1050*t1475
      t3867 = 32*t28*t616
      t3868 = t12*t596
      t3869 = t3868*t9
      t3870 = 16*t3869
      t3872 = 32*t28*t3855
      t3873 = -t3842+128*t5*t174-t3846+16*t3847*t975-t3851+16*t1490*t16-
     #t3857-64*t28*t174+t3861-8*t1061*t16-8*t1050*t975+t3867-t3870+t3872
      t3875 = 32*t1194*t387
      t3876 = t3*sll
      t3879 = t102*t751
      t3881 = 32*t99*t3879
      t3883 = 32*t281*t23
      t3884 = sll*t98
      t3885 = t3884*t100
      t3886 = t386*t865
      t3888 = 16*t3885*t3886
      t3891 = t6*t7*Eq
      t3893 = 4*t964*El*t3891
      t3894 = t711*sqq
      t3902 = t603*Eq
      t3905 = 4*t2074*t1002*t3902
      t3907 = 8*t688*t875
      t3912 = 64*t3*t603*t9
      t3915 = t6*t7*Mzsq
      t3917 = 64*t3*El*t3915
      t3918 = t103*sll
      t3922 = t3875-128*t3876*t9+t3881+t3883+t3888+t3893-16*t3894*t9-8*t
     #117*t118*t1799-16*t92*t657+t3905-t3907-8*t780*t3886+t3912+t3917-64
     #*t99*t102*t3918
      t3924 = t35*t603
      t3925 = 8*t3924
      t3926 = t899*t657
      t3927 = 16*t3926
      t3928 = t899*t757
      t3930 = t40*t603
      t3932 = 8*t63*t3930
      t3933 = t41*t3014
      t3935 = t57*t191
      t3937 = t97*t191
      t3938 = 32*t3937
      t3940 = 12*t53*t208
      t3941 = t57*t663
      t3943 = t41*t1635
      t3945 = t53*t816
      t3946 = 6*t3945
      t3947 = t97*t656
      t3948 = 16*t3947
      t3951 = -t3925-t3927-t943+32*t3928+t3932-16*t3933+20*t3935+t3938-t
     #3940-8*t3941+8*t3943+t3946-t3948+8*t820*t1639
      t3952 = t57*t831
      t3954 = t57*t243
      t3956 = t41*t608
      t3958 = t40*sll
      t3959 = t63*t3958
      t3961 = t40*Mzsq
      t3963 = 8*t63*t3961
      t3964 = t57*t656
      t3966 = t97*t596
      t3967 = 16*t3966
      t3968 = t53*t964
      t3969 = 6*t3968
      t3970 = t57*t596
      t3972 = t26*Mzsq
      t3973 = t3972*sqq
      t3976 = t790*t1320
      t3977 = 8*t3976
      t3979 = 2*t815*t1975
      t3981 = 4*t815*t1347
      t3982 = t808*t1320
      t3986 = -8*t3952+16*t3954+8*t3956-16*t3959+t3963-10*t3964-t3967+t3
     #969-10*t3970-2*t2113*t3973-t3977+t3979-t3981+8*t3982+8*t820*t1273
      t3990 = 2*t815*t1789
      t3991 = t790*t1859
      t3992 = 8*t3991
      t3993 = t808*t1859
      t3995 = t790*t298
      t4000 = 16*t820*t1279
      t4001 = t808*t298
      t4002 = 16*t4001
      t4005 = 2*t1433
      t4008 = t57*t1639
      t4011 = 32*t80*t1757
      t4013 = 32*t80*t1836
      t4016 = t3990-t3992+8*t3993+16*t3995+4*t1155*t1279-t4000-t4002-2*t
     #1155*t1639+t4005-2*t57*t1629-2*t4008-t4011+t4013-8*t48*t3924
      t4019 = t4*t1754
      t4020 = 4*t4019
      t4022 = 8*t53*t1737
      t4024 = 64*t44*t1836
      t4026 = 32*t44*t1757
      t4028 = 32*t44*t612
      t4030 = 64*t44*t1060
      t4032 = 8*t53*t1958
      t4033 = t57*t1273
      t4036 = 32*t80*t1060
      t4040 = 32*t80*t612
      t4041 = t675*sqq
      t4047 = t41*t603
      t4050 = 16*t4*t3924-t4020+t4022-t4024+t4026+t4028-t4030+t4032-2*t4
     #033+t4036-2*t57*t1818-t4040-8*t48*t4041-8*t71*t72*t603+8*t48*t4047
      t4053 = 64*t80*t505
      t4057 = 16*t48*t1511
      t4059 = 32*t4*t1511
      t4060 = 8*t1541
      t4061 = t57*t1279
      t4064 = 4*t57*t490
      t4066 = 16*t53*t352
      t4071 = 4*t3693
      t4072 = t4*t4041
      t4080 = t41*sll
      t4083 = -t4053+64*t80*t978+t4057-t4059+t4060+4*t4061+t4064-t4066+1
     #28*t44*t505-64*t44*t978-t4071+16*t4072+16*t71*t72*sll-8*t71*t72*Mz
     #sq-16*t48*t4080
      t4084 = t41*Mzsq
      t4087 = 2*t2171
      t4089 = t44*sll
      t4090 = 64*t4089
      t4091 = 8*t4047
      t4092 = 8*t4084
      t4094 = 32*t44*t603
      t4095 = 8*t4041
      t4096 = 8*t1331
      t4097 = 4*t3548
      t4098 = 16*t1511
      t4099 = 4*t1754
      t4100 = 32*t1336
      t4101 = 4*t1355
      t4102 = 8*t48*t4084-t1093-t4087-16*t4080-t4090+t4091+t4092+t4094-t
     #4095-t4096+t4097+t4098+t4099+t4100+t4101
      t4108 = t580*t2054
      t4110 = t4108*t569*t156
      A1e0c10p1t1n = -24*(t3873+t3922+t3951+t3986+t4016+t4050+t4083+t410
     #2)*t149*t151*t4110
      t4114 = Er*sll
      t4115 = t4114*t596
      t4120 = t386*t1799
      t4123 = t386*t2609
      t4126 = t102*t118
      t4129 = t12*El
      t4130 = t6*t1002
      t4133 = t3*sqq
      t4136 = t831*El
      t4142 = t37*t98*t100
      t4145 = 4*t1781+t3846-t3851+t3857+t3861-2*t323*t4115+2*t316*t4115-
     #8*t3885*t4120-8*t3885*t4123+64*t99*t4126-32*t4129*t4130+128*t4133*
     #t9+64*t4136*t6*t7*t19+32*t4142*t4123
      t4154 = Er*t603
      t4155 = t4154*t191
      t4166 = t6*t614
      t4171 = 32*t4142*t4120-64*t4142*t3886+64*t3444*t93*t603*t19-t3867+
     #t3870-t3872+2*t316*t4155-2*t323*t4155-128*t3444*t93*t396+32*t38*El
     #*t3891-16*t281*t4166-16*t281*t3915-t3875-t3881+t3883
      t4183 = t97*t125
      t4185 = t206*Mzsq
      t4188 = t831*t2
      t4193 = t3888-t3907-t3912-t3917+t3927-t1667-32*t4*t1364-64*t80*t28
     #4+64*t44*t128+2*t2232*t1609-128*t44*t122+32*t4183-8*t48*t4185+32*t
     #48*t4188-64*t80*t128
      t4212 = t671*t2675
      t4215 = t856*t6
      t4222 = t280*t404
      t4223 = t93*t608
      t4228 = 16*t203*t31
      t4229 = t97*t116
      t4233 = 64*t80*t122-16*t57*t125+16*t57*t196*t31-64*t220*t811*t37+1
     #6*t220*t1311-16*t187*t1311-16*t228*t1511+8*t228*t4185+4*t931*t4212
     #-16*t4215*t614*t2675-16*t857*t93*t2675+16*t4222*t4223-8*t1707-t422
     #8-64*t4229*t386*t58
      t4239 = t201*sqq
      t4248 = t4114*t31
      t4253 = t191*t209
      t4256 = t209*Mzsq
      t4257 = t4256*sll
      t4260 = t3884*t102
      t4261 = t103*Mzsq
      t4265 = t47*Eq
      t4269 = t171*El
      t4275 = t160*CsThr
      t4276 = Er*Mzsq
      t4277 = t4276*sll
      t4280 = t310*t46
      t4283 = 16*t4222*t93*t1635-64*t48*t4239+16*t1110*t1800+8*t170*t171
     #*BSnThl*t637+4*t323*t4248-8*t2330*t1292-4*t1155*t4253+2*t1155*t425
     #7+8*t4260*t4261*t249-4*t4265*t682*t1475-16*t170*t4269*t625-8*t1061
     #*t174+8*t4275*t4277-8*t4280*t4277
      t4284 = t306*t46
      t4313 = t171*Mzsq
      t4317 = t4276*sqq
      t4322 = t181*CsThr
      t4325 = 8*t4284*t4277+16*t4280*t318-32*t926*t671*t58-4*t931*t671*t
     #2675*t47-8*t1050*t637+8*t4260*t751*t249-16*t387*t250-64*t188*t616-
     #64*t188*t3855-16*t160*t3195*t637+32*t4229*t386*t1635+8*t170*t4313*
     #t174+16*t4275*t4317-16*t4275*t318-8*t4322*t4277
      t4338 = t12*t111
      t4341 = t4154*sll
      t4351 = t209*t603
      t4352 = t4351*sll
      t4357 = t4154*sqq
      t4364 = t32*t47
      t4365 = t53*t4364
      t4367 = 128*t188*t16-16*t160*t811*t174-128*t5*t16+32*t160*CsThr*El
     #*t625+32*t4338*t167+8*t4284*t4341+4*t170*CsThr*t1609-4*t935*t1609-
     #8*t2330*t1874+2*t1155*t4352-8*t808*t1856+16*t4275*t4357+8*t4275*t4
     #341-8*t4280*t4341-4*t4365
      t4368 = t831*t865
      t4369 = 16*t4368
      t4395 = 16*t1166
      t4397 = -t4369-64*t44*t58+32*t97*t246+8*t206*t608+8*t208*t249+8*t2
     #06*t1635+t804*t246+2*t582*t57+2*t804*t57+16*t41*t237-8*t41*t4256-8
     #*t57*t209*t31+16*t63*t40*sqq-t4395+4*t3827
      t4407 = 32*t44*t1635
      t4413 = 16*t67*t191
      t4415 = 32*t44*t608
      t4419 = t2570-4*t2232*t208+4*t1111-2*t1176-2*t828+t4407-8*t41*t435
     #1-t3932-4*t3935-2*t3945+t3948-t4413+t4415+4*t196*sll*t2970
      t4431 = 16*t1060*t334-4*t3710-t3963-4*t3964+t3967-2*t3968-4*t3970+
     #t3977-t3981-16*t3982+t3992-16*t3993-t4000+t4002-64*t1490*t856*t15
      t4433 = t603*El
      t4445 = Mzsq*t98
      t4469 = 64*t3*t4433*t167+16*t3847*t637+16*t1490*t174-16*t3868*t668
     #+64*t28*t16+32*t97*t4445*t265-64*t181*t225*t616+128*t21*t6*t177*sq
     #q-16*t899*t93*t1639-64*t197*t16+8*t9*t208*t4265-128*t3*t14*t167+t4
     #005-t1589+32*t197*t616
      t4472 = t306*t4317
      t4476 = 32*t197*t3855-16*t4472+t4011-t4013+t4024-t4026-t4028+t4030
     #-t4036+t4040+t4057-t4059+t4060+t4064+2*t608*t352
      t4490 = t663*t2
      t4495 = t206*t603
      t4503 = t97*t1639
      t4508 = t211*Mzsq
      t4516 = t2*t46
      t4517 = t26*t209
      t4521 = -16*t63*t231*sqq+8*t63*t231*Mzsq+128*t217*t287-32*t160*t17
     #1*t31+32*t48*t4490-16*t187*t1856+8*t228*t4495-8*t48*t4495+32*t80*t
     #196*t603-16*t4503-64*t210*t211*sqq+32*t210*t4508-64*t217*t4313+4*t
     #57*t4351*sqq-32*t4516*t4517*t603
      t4524 = t26*t603
      t4531 = t225*t603
      t4553 = t2*CsThr
      t4563 = 2*t53*t1688-2*t71*t4524*sll-4*t804*t865*t47-64*t220*t4531+
     #16*t220*t1856+8*t63*t231*t603-64*t217*t171*t603+32*t210*t211*t603-
     #2*t57*t4352-4*t3014*t2970+4*t57*t4253-2*t57*t4257+32*t4553*t26*t15
     #9*Mzsq-8*t170*t1347+4*t57*t4256*sqq
      t4573 = t26*sll
      t4574 = t4573*Mzsq
      t4589 = t57*t31
      t4591 = t53*t32
      t4593 = t97*t31
      t4595 = -32*t4516*t4517*Mzsq+64*t4516*t4517*sqq+2*t53*t1221-2*t71*
     #t4574+32*t187*t298+128*t220*t294+4*t874*t4212+64*t3*t404*t4223+4*t
     #2171+t4091-8*t4495-8*t4185+8*t4589+4*t4591-32*t4593
      t4596 = 32*t1364
      t4597 = 64*t4239
      t4598 = 32*t4188
      t4599 = 32*t4490
      t4602 = 8*t203*Mzsq*sqq
      t4605 = 8*t203*t603*sqq
      t4606 = t41*sqq
      t4608 = t44*sqq
      t4610 = t1433*t47
      t4612 = t4596-t4597+t4598+t4599+t4602+t4605-16*t4606+64*t4608-2*t4
     #610+t4092-t4094-t4096+t4098-t4100+t4101
      t4619 = t567*t151
      t4620 = t594*t2052
      A1e0c11p1t1n = 24*(t4145+t4171+t4193+t4233+t4283+t4325+t4367+t4397
     #+t4419+t4431+t4469+t4476+t4521+t4563+t4595+t4612)*t362*t149*t4619*
     #t4620
      t4623 = t27*t2149
      t4626 = t1578*El
      t4627 = t12*t4626
      t4628 = t603**2
      t4629 = t7*t4628
      t4630 = t4629*Mzsq
      t4631 = t6*t4630
      t4635 = t6*t614*t2063
      t4638 = sqq*BSnThl
      t4640 = t27*t4638*El
      t4645 = t27*t31*BSnThl
      t4646 = t8*t3306
      t4648 = 16*t4645*t4646
      t4653 = t27*t68*El
      t4660 = t31*t4628
      t4664 = t152*El
      t4671 = t27*t2190
      t4674 = t614*t152
      t4675 = t6*t4674
      t4678 = t27*t493
      t4679 = t614*Mzsq
      t4680 = t2074*t4679
      t4683 = t15*t1677
      t4688 = t152*sll
      t4693 = -8*t4623*t616-4*t4627*t4631-4*t4627*t4635+32*t4640*t615*t1
     #276-t4648-24*t27*t1273*t999-24*t4653*t4631-24*t27*t234*El*t4631-96
     #*t4222*t93*t4660+4*t12*t4664*t4631-72*t27*t1818*t999+48*t4671*t174
     #-72*t4653*t4675+16*t4678*t4680+96*t4653*t4683-24*t4671*t616-8*t12*
     #t4688*BSnThl*t4646
      t4694 = t4688*t46
      t4697 = t6*t4679
      t4700 = t152*BSnThl
      t4710 = t623*Mzsq
      t4711 = t2074*t4710
      t4716 = t614*sll
      t4717 = t2074*t4716
      t4721 = t27*t445*El
      t4722 = t6*t4716
      t4730 = t6*t4710
      t4735 = t206*t1093
      t4736 = 72*t4735
      t4740 = t905*sqq
      t4741 = t4740*t603
      t4743 = 2*t4741*Mzsq
      t4750 = t2157*t663
      t4753 = t603*t4688
      t4754 = t439*t4753
      t4757 = 8*t12*t4694*El*t4697-16*t27*t4700*t4646-8*t27*t3014*El*t46
     #31+32*t4678*t999-24*t4678*t4711+32*t27*sqq*t612*t4717-72*t4721*t47
     #22+48*t4721*t4697+96*t27*t445*t999-72*t4721*t4730-24*t4678*t4717+t
     #4736+144*t4222*t93*t1093-t4743-24*t27*t3422*t616-8*t27*t3015*t616+
     #16*t143*t4750+12*t143*t4754
      t4759 = t4628*Mzsq
      t4760 = t48*t4759
      t4763 = t206*t417
      t4768 = t1211*CsThr
      t4769 = t3*t4768
      t4774 = t8*t635*t2063
      t4777 = t27*t128
      t4778 = t4629*sll
      t4779 = t2074*t4778
      t4782 = t2618*t615
      t4790 = t683*t998
      t4793 = t2618*t1493
      t4796 = t2074*t4630
      t4799 = t2074*t4674
      t4802 = t1578*t47
      t4806 = t1071*t4782
      t4808 = t2622*t1045
      t4811 = t2063*t98
      t4812 = t206*t4811
      t4815 = Er*t417
      t4816 = t689*t4815
      t4819 = sll*t43
      t4820 = t4819*t2100
      t4823 = 64*t4608*t4760-24*t4763*t159*t32*t603+128*t4769*t4680+4*t1
     #2*t4700*t4774+16*t4777*t4779-8*t4777*t4782+48*t27*t4626*t4697+64*t
     #979*t4782-24*t4777*t4790+16*t4777*t4793-8*t4777*t4796-24*t4777*t47
     #99-128*t3*t4802*t4680-64*t4806-32*t313*t4808+192*t4812*t4126-48*t4
     #816*t1198-384*t4820*t16
      t4825 = t751*Mzsq
      t4826 = t117*t4825
      t4829 = t103*t4628
      t4830 = t4829*Mzsq
      t4831 = t117*t4830
      t4834 = CsThr*t19
      t4837 = sqq*t4628
      t4838 = t4837*Mzsq
      t4839 = t15*t4838
      t4842 = t206*t4445
      t4847 = t4819*t2089
      t4850 = t27*t2083
      t4853 = t152*t43
      t4854 = t4853*t1041
      t4859 = t4628*t1
      t4860 = t4859*t964
      t4865 = t423*t4759
      t4868 = t33*t4269
      t4869 = t4628*BCsThl
      t4870 = t4869*Mzsq
      t4874 = t423*t1095
      t4879 = t804*t2166
      t4880 = 8*t4879
      t4883 = t418*t4759
      t4886 = 64*t97*t4768*t4826-32*t313*t4831-64*t504*t4834*El*t4839+38
     #4*t4842*t1034-96*t4812*t3879+192*t4847*t23-16*t4850*t174+576*t4854
     #*t605-64*t1071*t4796+8*t1110*t4860-16*t1110*t3084+16*t560*t4865+16
     #*t4868*t15*t4870+48*t143*t4874+192*t4820*t616-t4880+64*t4089*t4760
     #-32*t4089*t4883
      t4892 = t418*t2469
      t4895 = t80*t31
      t4896 = t4*t582
      t4899 = t4819*t1060
      t4905 = t865*t211*t2063
      t4909 = t97*t152
      t4910 = t4*t2300
      t4913 = t4*t4759
      t4916 = t640*t663
      t4919 = t3073*t47
      t4925 = t44*t152
      t4926 = t418*t791
      t4929 = t97*t4688
      t4930 = t4*t791
      t4933 = t439*t791
      t4936 = t48*t791
      t4939 = t80*sqq
      t4940 = t4628*sll
      t4941 = t48*t4940
      t4944 = t4*t804
      t4947 = -16*t1110*t4859*t208+4*t4589*t4892-24*t4895*t4896+576*t489
     #9*t917-768*t4899*t1078-32*t4905*t3444*t1493-4*t4909*t4910-4*t4909*
     #t4913-32*t4909*t4916+8*t4919*t1119+24*t1110*t1799*t1436+64*t4925*t
     #4926+8*t4929*t4930-8*t4929*t4933-128*t4925*t4936+16*t4939*t4941-24
     #*t4895*t4944
      t4948 = t2939*t4759
      t4951 = t262*t47
      t4952 = t445*t4759
      t4955 = t4951*t2322
      t4959 = t418*t1095
      t4962 = t57*t32
      t4963 = t418*t804
      t4968 = t418*t2300
      t4971 = t418*t1233
      t4974 = t418*t582
      t4979 = t439*t4759
      t4982 = t439*t4940
      t4987 = t439*t804
      t4994 = t152*t2
      t4995 = t4994*t417
      t4998 = -8*t4763*t4948-4*t4951*t4952-4*t4955+32*t4763*t2325-6*t458
     #9*t4959-6*t4962*t4963+4*t4962*t4926-2*t4589*t4968-6*t4589*t4971-6*
     #t4962*t4974+2*t4589*t4760+2*t4589*t4979-4*t4589*t4982-4*t4962*t493
     #3+6*t4962*t4987-8*t4939*t4760+16*t4951*t2319+16*t4995*t1142
      t5002 = t196*t4940
      t5005 = t196*t4759
      t5008 = t4*t4940
      t5013 = t4962*t4930
      t5019 = t418*t4940
      t5022 = t48*t804
      t5025 = t439*t2300
      t5030 = t48*t2300
      t5032 = 64*t4089*t5030
      t5037 = t41*t31
      t5038 = t5037*t4930
      t5042 = t4*t816
      t5043 = t4080*t5042
      t5045 = -8*t4763*t2332-16*t4939*t5002+8*t4939*t5005+12*t4589*t5008
     #-6*t4589*t4913+12*t5013-18*t4962*t4944-2*t4589*t4883+4*t4589*t5019
     #+6*t4962*t5022+4*t4909*t5025-32*t4089*t4968+t5032-16*t4606*t4979+1
     #6*t4606*t4913-32*t5038+32*t5037*t4933+48*t5043
      t5046 = t4*t4837
      t5061 = t41*t152
      t5062 = t4*t656
      t5068 = t4089*t4910
      t5076 = t161*t37
      t5080 = t505*t3067
      t5085 = t104*t1677
      t5088 = t386*t152
      t5095 = -32*t4080*t5046+32*t4080*t439*t4837-48*t4080*t439*t816-16*
     #t4606*t5025+16*t4080*t4910-16*t4080*t4979+48*t5061*t5062-48*t5061*
     #t439*t656-32*t5068-32*t4089*t4913+32*t5061*t4933-16*t4080*t5025-2*
     #t5076*t4660*t800+12*t5080+12*t505*t4860+128*t143*t1031*t5085-96*t1
     #23*t780*t5088-24*t125*t1799*t1233
      t5105 = t53*t152
      t5106 = t68*t2300
      t5107 = t5105*t5106
      t5108 = t68*t4759
      t5112 = t43*t31
      t5113 = t5112*Er
      t5116 = t206*t46
      t5117 = t58*t4759
      t5120 = t206*BCsThl
      t5121 = t5120*t2274
      t5123 = t4994*BCsThl
      t5126 = t4*t1095
      t5129 = t4608*t4910
      t5131 = t1713*BCsThl
      t5138 = t297*t4759
      t5141 = t1713*t417
      t5146 = -2*t5076*t2064*t796+4*t4531*t1436*t800+t5107+t5105*t5108+6
     #4*t4925*t4930+64*t5113*t4930+40*t5116*t5117-56*t5121+112*t5123*t81
     #2-96*t4608*t5126-32*t5129-384*t5131*t812-96*t5113*t4944+224*t5120*
     #t2277-56*t5120*t5138-128*t5141*t1142-32*t4608*t4913
      t5159 = t4*t2083
      t5168 = t4080*t3631
      t5172 = t5116*t2267
      t5178 = t152*t19
      t5180 = t7*t32
      t5181 = t6*t5180
      t5188 = t520*t31
      t5189 = t519*t5188
      t5192 = t152*Eq
      t5197 = t103*t31
      t5198 = t102*t5197
      t5201 = t4*t2469
      t5202 = t4608*t5201
      t5204 = 24*t1110*t3041*t804-48*t545*t1799*t582+32*t4080*t439*t2083
     #-48*t5061*t439*t596-32*t4080*t5159+48*t5061*t4*t596-48*t4080*t439*
     #t964+48*t5168+16*t4606*t4910+40*t5172-160*t5116*t2264-96*t5113*t48
     #96-192*t5178*El*t5181+384*t4853*El*t4130+128*t5178*t515*t5189+32*t
     #5192*t905*t9+384*t4994*t98*t5198+64*t5202
      t5206 = t4*t1233
      t5209 = t2157*t243
      t5212 = t640*t2141
      t5216 = t103*t32
      t5217 = t102*t5216
      t5220 = t4*t3018
      t5223 = t4688*Mzsq
      t5224 = t4*t5223
      t5227 = t48*t2469
      t5229 = 128*t4608*t5227
      t5230 = t48*t1233
      t5233 = t439*t3018
      t5236 = t439*t5223
      t5239 = t423*t582
      t5240 = t4591*t5239
      t5249 = t141*t5209
      t5251 = t141*t5212
      t5256 = t4593*t5201
      t5258 = -96*t4608*t5206-32*t143*t5209+48*t143*t5212-96*t3073*t98*t
     #5217+8*t143*t5220-12*t143*t5224-t5229+192*t4608*t5230-8*t143*t5233
     #+12*t143*t5236-12*t5240+64*t4608*t4892-96*t4608*t4971-16*t4591*t4*
     #t2141+8*t5249-12*t5251+64*t4593*t439*t2141-24*t5256
      t5265 = t48*t582
      t5268 = t19*t32
      t5269 = t5268*Er
      t5272 = Er*t47
      t5273 = t5272*BSnThl
      t5274 = t33*t5273
      t5275 = t375*Mzsq
      t5276 = t4629*SnThr
      t5280 = t33*Er
      t5296 = t4591*t5227
      t5300 = Eq*t905
      t5301 = t5300*Er
      t5304 = t439*t582
      t5307 = t439*t2469
      t5310 = 24*t4593*t5206+48*t4593*t5227-72*t4593*t5230+192*t5113*t52
     #65+36*t5269*t4896-8*t5274*t5275*t5276+32*t5280*t651*t3266+24*t1110
     #*t3050-72*t5269*t5265-96*t5113*t4974+24*t4593*t4971+24*t5269*t4974
     #-16*t4593*t4892-4*t5296+6*t4591*t5230+6*t5301*t5265-12*t5269*t5304
     #+8*t4593*t5307
      t5315 = t375*t7
      t5316 = t4628*SnThr
      t5321 = t3041*t5273
      t5325 = Er*BCsThl
      t5326 = t684*t5325
      t5329 = t196*t582
      t5337 = t4589*t5201
      t5339 = t196*t2469
      t5342 = t196*t1233
      t5348 = t439*t1233
      t5363 = 48*t4593*t5239+16*t5274*t5315*t5316*sll-24*t5321*t5275*t14
     #93+12*t5326*t2388+24*t4895*t5329+24*t1110*t1436*t2609-16*t2100*t11
     #47+12*t5337-16*t4939*t5339+24*t4939*t5342-24*t5274*t5275*t998-6*t4
     #589*t5348-4*t4589*t5307-12*t145*t5236+8*t145*t5233+12*t145*t5224-8
     #*t145*t5220-24*t4895*t5265
      t5378 = t4589*t5206
      t5382 = t4939*t5201
      t5384 = t4939*t5206
      t5386 = t4589*t5227
      t5389 = 64*t1713*t4837
      t5390 = t375*t2063
      t5394 = t31*t2
      t5395 = t5325*CsThr
      t5396 = t5394*t5395
      t5399 = Er*BSnThl
      t5400 = t33*t5399
      t5405 = Er*t2063
      t5411 = 72*t4895*t5304+24*t4895*t4974-6*t4962*t5329+4*t4589*t5339-
     #6*t4589*t5342-16*t4939*t4892+24*t4939*t4971-6*t5378+6*t4589*t5230+
     #16*t5382-24*t5384-4*t5386+t5389+16*t5274*t5390*t1493-192*t5396*t47
     #11-16*t5400*t8*SnThr*t4688-8*t865*t5405*t999-18*t4962*t4896
      t5421 = t4939*t5227
      t5423 = t4939*t5230
      t5441 = t27*t58*El
      t5442 = t6*t4778
      t5445 = Mzsq*El
      t5446 = t4819*t5445
      t5449 = sll*t417
      t5453 = t4640*t3218
      t5455 = t4645*t969
      t5462 = 6*t4962*t5304+6*t4962*t5265-48*t4939*t5307+72*t4939*t5348+
     #16*t5421-24*t5423-16*t3041*t5399*t969+12*t865*t5399*t8*t635*t4688-
     #24*t4995*t1375+16*t4080*t4913+32*t27*t246*t999+16*t5441*t5442-288*
     #t5446*t4130+64*t41*t5449*t801-32*t5453-32*t5455+8*t57*t1371*t792+1
     #6*t560*t4750
      t5465 = t544*t603
      t5470 = t4*t831
      t5473 = t27*t14
      t5478 = t53*t4688
      t5479 = t68*t791
      t5494 = t75*t53
      t5498 = t152*t417
      t5506 = t4553*t98
      t5507 = t2519*t5506
      t5508 = t4829*sll
      t5517 = t396*t4815
      t5520 = -2*t5105*t1119-168*t5120*t5465-32*t5061*t4930-32*t5470*t13
     #69+8*t5473*t4631+128*t5396*t4680-2*t5478*t5479-32*t4868*t15*t3120+
     #48*t4868*t15*t1691-64*t5280*t3196*t615*t1317-8*t4338*t5442-16*t549
     #4*t225*t1677+4*t57*t5498*t1142-32*t123*t4808+24*t5400*t948-192*t55
     #07*t102*t5508+16*t1707*t404*t93*t3902*Mzsq-16*t5517*t2373
      t5531 = t47*t43
      t5536 = t1115*t582
      t5546 = t211*t47
      t5547 = t35*t5546
      t5554 = t404*t93*t4688
      t5571 = t4660*Mzsq
      t5572 = t93*t5571
      t5583 = 192*t524*t6*t1002*t791+192*t3213*t3210-96*t68*t5531*El*t46
     #31+2*t4591*t5536+16*t143*t4865+32*t5141*t2766+64*t123*t2622*t1521-
     #128*t5547*t4779+4*t4338*t4631+16*t33*t211*t5554+16*t4868*t15*t2769
     #-8*t2680*t93*t403*t603-8*t4215*t4629*t2683-96*t68*t5531*t2063*t616
     #-64*t2025*t5572+16*t116*t386*t32*t76*t1862+256*t4133*t1049*t3266
      t5586 = t35*t5395
      t5587 = t5586*t4793
      t5589 = t5400*t4774
      t5591 = t1115*t1233
      t5594 = t4628*t152
      t5595 = t4*t5594
      t5611 = t211*El
      t5612 = t33*t5611
      t5616 = 96*t1713*t816
      t5617 = CsThr*Eq
      t5618 = t2067*El
      t5622 = t5617*t5618*t6*t614*BCsThl
      t5628 = CsThr*t43
      t5633 = t4819*t621
      t5635 = t6*t4629*sqq
      t5642 = 64*t4608*t5008+128*t5587+4*t5589+8*t143*t5591+8*t143*t5595
     #+112*t5120*t1788*t4628-192*t5396*t4717-24*t4763*t1141*t32+8*t4905*
     #t999+48*t4868*t15*t1237+4*t5612*t4635-t5616+4*t5622+4*t5617*t2074*
     #t5180*t4870+192*t2519*t5628*t2063*t616+192*t5633*t5635-96*t5633*t4
     #631-384*t4854*t4697
      t5644 = t2613*t31
      t5654 = t4819*t1041
      t5657 = sll*Eq
      t5661 = t93*t4838
      t5664 = t4819*t2158
      t5668 = t3444*t615
      t5673 = t1578*t37
      t5677 = t47*t4628
      t5678 = t5677*t582
      t5681 = t27*t978
      t5684 = t15*t1620
      t5687 = t4962*t4768
      t5690 = t926*t671*t2064
      t5696 = t33*t5405
      t5697 = t5696*t174
      t5699 = t93*t2003
      t5702 = 16*t48*t5644*t616-768*t4853*t4*t917+384*t4853*t48*t917-384
     #*t5654*t5635+16*t5657*t2136*t917+48*t4222*t5661-96*t5664*t616+16*t
     #12*t3422*t5668+192*t5664*t16-32*t12*t5673*t4680+2*t85*t5678-8*t568
     #1*t4796+384*t5633*t5684+24*t5687-96*t5690-8*t5612*t6*t623*t2063-8*
     #t5697-32*t1452*t5699
      t5709 = t6*t7*CsThr*BCsThl
      t5710 = t5657*t5618*t5709
      t5712 = t37*t4628
      t5720 = t4853*t621
      t5722 = t6*t22*Mzsq
      t5727 = t262*t116
      t5732 = t57*t75*t792
      t5740 = t3014*t2300
      t5741 = t85*t5740
      t5743 = t80*t2519
      t5744 = t417*t603
      t5745 = t5744*t582
      t5748 = t57*t1994
      t5758 = t865*t5611
      t5762 = -48*t5517*t1123*t1095-8*t5710+4*t5494*t5712*t811-8*t3699*t
     #4019-8*t5326*t2391-288*t5720*t5722+192*t5720*t4697-48*t5727*t386*t
     #4660+16*t5732+16*t12*t234*t37*t4796-288*t5720*t605+2*t5741-96*t574
     #3*t5745-6*t5748-288*t4819*t612*t917+16*t27*t4802*t4680-8*t5681*t47
     #82+8*t5758*t15*t2830
      t5773 = t181*t3958
      t5776 = t1436*t47
      t5777 = t80*t5776
      t5780 = t80*t68*t2346
      t5786 = 32*t5743*t792
      t5791 = t6*t4629*t152
      t5799 = t386*t31
      t5800 = t116*t5799
      t5804 = t5399*El
      t5819 = -12*t5758*t6*t614*t4688+64*t979*t4796-288*t5633*t1004-64*t
     #5773*t812+32*t5777+32*t5780-8*t57*t445*t2346+t5786-16*t3041*t211*t
     #999+8*t5758*t5791-8*t404*t93*t905*t5617*t1862-8*t5800*t4859*t1060-
     #24*t33*t5804*t3218+24*t965*t999+192*t5654*t4631-8*t5657*t3444*t566
     #1-96*t4183*t117*t3918*Mzsq
      t5820 = t5394*t5546
      t5825 = t116*t5088
      t5837 = 16*t262*t2209
      t5839 = 24*t262*t2000
      t5841 = 24*t3073*t2171
      t5844 = t4869*sll
      t5848 = t53*t234
      t5849 = t5712*t596
      t5858 = t4940*Mzsq
      t5859 = t80*t5858
      t5861 = t80*t2810
      t5863 = t5712*t31
      t5865 = 16*t262*t5863
      t5866 = -128*t5820*t4680+48*t5269*t4936+128*t4183*t5825+192*t5820*
     #t4717+64*t4183*t4826-96*t5113*t4963+576*t5654*t1004+t5837-t5839-t5
     #841+192*t5820*t4711-32*t4868*t15*t5844-4*t5848*t5849+4*t97*t5594*M
     #zsq-8*t97*t4753*Mzsq+8*t5859+8*t5861+t5865
      t5868 = t48*t32
      t5871 = t152*t370
      t5877 = t48*t905
      t5880 = t3041*t4269
      t5884 = t4114*t2063
      t5885 = t684*t5884
      t5887 = t80*t2830
      t5889 = t211*t152
      t5892 = t211*t4688
      t5895 = t5300*t4277
      t5897 = t689*t5884
      t5899 = Er*t152
      t5900 = t5899*Mzsq
      t5901 = t689*t5900
      t5902 = 24*t5901
      t5903 = t5268*t4277
      t5905 = t5300*t4341
      t5908 = 6*t684*t5900
      t5909 = t5405*t152
      t5912 = Er*t4628
      t5913 = t5912*Mzsq
      t5914 = t689*t5913
      t5916 = t5912*t152
      t5919 = -160*t4994*t5868-128*t5871*t418*t31-384*t5871*t125+16*t307
     #3*t5877+48*t5880*t15*t1891+2*t5885-16*t5887-16*t5268*t5889+16*t689
     #*t5892-3*t5895-8*t5897+t5902+12*t5903-3*t5905-t5908-8*t396*t5909+4
     #*t5914+2*t403*t5916
      t5920 = t1390*t2893
      t5921 = 24*t5920
      t5924 = t4276*t4688
      t5929 = t2469*t3512
      t5931 = sll*t4740
      t5936 = t3073*t905
      t5939 = t5871*t31
      t5942 = t1095*t2063
      t5947 = t57*t3610
      t5949 = t80*t1589
      t5951 = t80*t2893
      t5953 = t57*t3021
      t5955 = t80*t1298
      t5957 = t32*t152
      t5959 = t57*t5957*t47
      t5962 = t57*t3006*t152
      t5964 = t5921-4*t684*t5892-3*t403*t5924+8*t2141*t3512+2*t5929-3*t5
     #931*t1060-3*t3511*t1237-48*t4*t5936+384*t48*t5939+4*t97*t5942+96*t
     #80*t1585+4*t5947+24*t5949-16*t5951-6*t5953+24*t5955-8*t5959-24*t59
     #62
      t5979 = t684*t5913
      t5980 = t4276*t603
      t5981 = t5300*t5980
      t5983 = t4154*t2063
      t5984 = t689*t5983
      t5986 = t5268*t5980
      t5988 = t5912*sll
      t5989 = t689*t5988
      t5991 = t4154*t152
      t5992 = t689*t5991
      t5993 = 24*t5992
      t5995 = 6*t684*t5991
      t6001 = t865*t211*t37
      t6006 = t48*t1095
      t6009 = t370*t4628
      t6010 = t6009*Mzsq
      t6013 = 4*t3511*Mzsq*t1891+48*t5880*t15*t1317+48*t5612*t15*t1961+6
     #4*t5113*t4926-t5979+2*t5981+4*t5984-8*t5986-8*t5989+t5993-t5995+12
     #*t396*t5924+4*t5300*t5889+16*t6001*t4796-32*t6001*t4779-72*t4593*t
     #6006+32*t487*t6010
      t6014 = t4740*CsThr
      t6019 = t3511*t2769
      t6026 = t4154*t4688
      t6029 = t684*t5983
      t6030 = t684*t5988
      t6033 = t80*t4837*sll
      t6035 = t80*t4838
      t6037 = t5268*t4341
      t6040 = 4*t663*t3021
      t6042 = 2*t243*t5571
      t6044 = 2*t2469*t1781
      t6048 = 8*t2136*t1754
      t6050 = -3*t6014*t1891+2*t6014*t1862-t6019-3*t3511*t1691+2*t3511*t
     #5844-8*t396*t5916+12*t396*t6026-t6029+2*t6030-16*t6033+8*t6035+12*
     #t6037-t6040+t6042+t6044-3*t403*t6026-t6048-t3511*t4870
      t6053 = 8*t663*t1994
      t6055 = 8*t5712*t3548
      t6056 = t80*t1672
      t6057 = 16*t6056
      t6058 = t57*t5571
      t6060 = t80*t2761
      t6062 = t80*t1670
      t6064 = t57*t2003
      t6066 = t80*t1667
      t6068 = t1436*t603
      t6069 = t57*t6068
      t6081 = t206*t4838
      t6082 = 24*t6081
      t6084 = 24*t3073*t1781
      t6085 = t4994*t1620
      t6086 = 48*t6085
      t6088 = 16*t3073*t1672
      t6089 = t57*t2805
      t6091 = t4660*sll
      t6092 = t57*t6091
      t6094 = -t6053-t6055-t6057-2*t6058+8*t6060+24*t6062+4*t6064+24*t60
     #66-6*t6069-64*t5612*t6*t177*t152+48*t5612*t15*t1847-32*t5612*t15*t
     #1774+t6082-t6084-t6086+t6088-2*t6089+4*t6092
      t6095 = t57*t2002
      t6097 = t206*t1672
      t6098 = 96*t6097
      t6099 = t1799*Mzsq
      t6101 = 16*t1707*t6099
      t6106 = 8*t262*t5571
      t6108 = 32*t262*t2003
      t6110 = 128*t1713*t1620
      t6116 = t2100*t1720
      t6121 = t5957*t1
      t6132 = t27*t4638*t948
      t6136 = -6*t6095-t6098+t6101-64*t87*t994*t4683-t6106+t6108+t6110+4
     #8*t6001*t4799-32*t5880*t15*t1862-24*t6116-24*t27*t1221*t174-32*t11
     #10*t6121-t5105*t4838-288*t538*t1589-16*t5473*t5442+8*t4850*t616+24
     #*t6132+224*t4994*t545
      t6145 = t4591*t3369
      t6146 = t46*t4628
      t6147 = t6146*t152
      t6150 = t1115*t4688
      t6153 = t5105*t2761
      t6159 = t905*t1
      t6160 = t6159*Mzsq
      t6170 = t27*t112
      t6175 = t97*t1190
      t6183 = t152*t47
      t6184 = t53*t6183
      t6189 = 12*t5269*t1623-8*t5269*t1626+4*t4593*t3369-t6145+8*t143*t6
     #147-12*t143*t6150-t6153+2*t5478*t1620-288*t128*t1642*t152+36*t505*
     #t6160+16*t4905*t5668+24*t27*t964*t174+24*t4645*t637-24*t6170*t614*
     #t671*t47+8*t6175*t2549*t4837+8*t6175*t417*t2063*t656-4*t6184*t1119
     #+192*t538*t1672
      t6190 = t1799*sll
      t6195 = t1799*t152
      t6200 = t6146*sll
      t6203 = t6146*Mzsq
      t6209 = t143*t5858
      t6210 = 8*t6209
      t6211 = t143*t2830
      t6213 = t4591*t1677
      t6215 = t53*t3014
      t6219 = 24*t262*t664
      t6220 = t865*t4508
      t6224 = t3041*t5611
      t6236 = 36*t3512*t6190-24*t3512*t6099+36*t545*t6195-64*t550*t1672+
     #2*t4591*t6200-t4591*t6203-3*t5301*t1623+2*t5301*t1626+t6210-24*t62
     #11+6*t6213+2*t6215*t4952-t6219+48*t6220*t3444*t998+12*t6224*t4730-
     #64*t5880*t6*t7*t152*BCsThl-8*t6224*t4697+4*t5612*t4631
      t6248 = t4859*sll
      t6251 = t4859*Mzsq
      t6254 = t1390*Mzsq
      t6259 = t4853*t404
      t6260 = t93*t791
      t6263 = t2063*t515
      t6264 = t280*t6263
      t6269 = t4994*t116
      t6270 = t386*t596
      t6275 = t4593*t1677
      t6285 = -8*t5612*t5442-96*t978*t3356-192*t546*t1672-56*t473*t5571+
     #40*t1327*t5571-24*t545*t6248+12*t545*t6251-80*t5868*t6254-288*t538
     #*t1667+192*t6259*t6260+64*t6264*t519*t520*sqq-288*t6269*t6270+12*t
     #6224*t4722-24*t6275-8*t4593*t6200+4*t4593*t6203+3*t141*t6150-2*t14
     #1*t6147
      t6297 = t515*t517
      t6298 = t5178*t6297
      t6304 = t280*t6297
      t6311 = t1123*t4759
      t6314 = t3201*CsThr
      t6317 = t225*t4628
      t6318 = t1433*BCsThl
      t6323 = t473*t2805
      t6325 = t2079*t603
      t6326 = t545*t6325
      t6328 = t1327*t2805
      t6336 = 96*t122*t6010+192*t128*t6009*sll-96*t128*t6010-96*t280*Mzs
     #q*t515*t5189-96*t6298*t698*t596+32*t3770*t616+128*t6304*t698*t1620
     #+64*t6298*t698*t791-16*t5517*t6311+4*t6314*t1862+4*t6317*t6318+112
     #*t545*t6254-56*t6323+12*t6326+40*t6328-8*t1760*t5571-8*t1760*t2805
     #-96*t6298*t704
      t6342 = t32*El
      t6343 = t5192*t6342
      t6346 = t386*t791
      t6349 = t5657*t404
      t6357 = t418*t152
      t6364 = t943*BCsThl
      t6367 = t418*t32
      t6370 = t2113*t3780
      t6372 = t80*t152
      t6377 = t386*t656
      t6380 = t206*t116
      t6381 = t386*t816
      t6386 = t93*t596
      t6389 = t386*t4759
      t6392 = -24*t6343*t3915+192*t6269*t6346-24*t6349*t93*t959+8*t504*t
     #225*t4859*t596-32*t5269*t6357+16*t4591*t5673-8*t5301*t4802-2*t6317
     #*t6364+16*t6367*t6254-4*t6370+48*t6372*t1626-6*t6314*t1891-288*t62
     #69*t6377-288*t6380*t6381-288*t6259*t657-288*t6259*t6386-96*t6380*t
     #6389
      t6393 = t4819*t404
      t6397 = t143*t2810
      t6398 = 8*t6397
      t6399 = t4*t4688
      t6408 = t640*t4688
      t6413 = t4591*t3422
      t6428 = t262*t117
      t6434 = 32*t1713*t4759
      t6438 = 192*t6393*t93*t4837+t6398+16*t4593*t6399-8*t4593*t3422-3*t
     #5301*t1821+12*t5269*t1821-12*t143*t6408+8*t143*t3578+2*t6413-12*t6
     #220*t5554-96*t4847*t4166+384*t5446*t605-192*t5586*t4790-8*t6349*t5
     #572+64*t6304*t698*t4837-8*t6428*t118*t4628*t1060-t6434-192*t2025*t
     #93*t6068
      t6441 = t698*t4759
      t6448 = t15*t1836
      t6449 = t5657*t2064*El*t6448
      t6451 = t5657*t2074
      t6452 = t4628*CsThr
      t6453 = t6452*BCsThl
      t6457 = t6342*BSnThl
      t6458 = t5192*t6457
      t6459 = t93*t1836
      t6462 = t3073*t117
      t6470 = 64*t1713*t2083
      t6471 = t206*t943
      t6472 = 72*t6471
      t6473 = t57*t5957
      t6475 = t80*t1436
      t6478 = 6*t663*t5957
      t6480 = 2*t2136*t1093
      t6481 = t370*t31
      t6483 = 64*t6481*t791
      t6484 = t5394*t4759
      t6485 = 24*t6484
      t6487 = 96*t1713*t964
      t6489 = 24*t262*t922
      t6491 = 32*t3073*t1707
      t6492 = -32*t396*t6297*t6441-384*t3221*t3210+4*t6449-8*t6451*t5180
     #*t6453+12*t6458*t6459+16*t6462*t118*CsThr*t1862-32*t6304*t6441+t64
     #70+t6472+8*t6473-32*t6475+t6478+t6480+t6483+t6485-t6487-t6489+t649
     #1
      t6493 = t4994*t964
      t6494 = 72*t6493
      t6496 = 24*t3073*t943
      t6498 = 96*t5871*t596
      t6499 = t1390*t2064
      t6500 = 24*t6499
      t6501 = t804*t2067
      t6502 = t1093*t2079
      t6503 = 8*t6502
      t6504 = t1433*t2079
      t6505 = 16*t6504
      t6506 = t208*t2625
      t6507 = 48*t6506
      t6509 = 96*t5871*t656
      t6511 = 24*t3073*t1093
      t6512 = t4994*t816
      t6513 = 72*t6512
      t6514 = t32*t4628
      t6516 = 8*t6514*t2609
      t6517 = t32*t2
      t6518 = t6517*t791
      t6519 = 48*t6518
      t6521 = 16*t6159*t791
      t6523 = 4*t5712*t1433
      t6525 = 4*t3201*t791
      t6527 = 6*t3201*t804
      t6529 = 2*t943*t1095
      t6530 = t6494-t6496-t6498+t6500+t6501-t6503+t6505-t6507-t6509-t651
     #1+t6513-t6516-t6519+t6521-t6523-t6525+t6527-t6529
      t6533 = t370*sqq
      t6535 = 32*t6533*t2300
      t6537 = 2*t5712*t943
      t6539 = 4*t922*t804
      t6540 = t5268*t5899
      t6542 = t5300*t5899
      t6544 = Er*t4688
      t6547 = t206*t4660
      t6548 = 48*t6547
      t6560 = t4*t5644
      t6562 = 128*t6560*t174
      t6563 = t865*t5804
      t6566 = t6563*t614*t671*t2063
      t6571 = t102*t4830
      t6575 = t519*t702
      t6582 = -t6535+t6537-t6539-16*t6540+4*t6542-16*t689*t6544-t6548+4*
     #t684*t6544-24*t6462*t5197*t1836-192*t2025*t93*t1994-96*t35*t4811*t
     #3879+t6562+8*t6566+16*t9*t5300*t791+96*t5507*t6571-32*t396*t6263*t
     #6575+64*t516*t519*t5188*t791
      t6583 = t20*t6
      t6587 = t386*t1672
      t6602 = 64*t6560*t616
      t6612 = t6215*t2322
      t6620 = sqq*t43
      t6625 = t4*t711*t31
      t6633 = 48*t6583*t1002*t4759-192*t3244*t6587+192*t496*t102*t5197*t
     #791-192*t5586*t4799+256*t5396*t999+288*t3244*t386*t1589-t6602+192*
     #t2519*t5628*El*t4631+8*t6563*t93*t5858-24*t6563*t3188+2*t6612+96*t
     #2164*t1046+288*t4072*t5825-192*t2164*t1522-96*t6620*t2089*t4166-19
     #2*t6625*t999+128*t2025*t5699+288*t3244*t386*t1667
      t6647 = t5586*t4782
      t6657 = t896*t742*t865*t2063
      t6660 = 32*t6533*t4759
      t6663 = t280*t943
      t6666 = t2519*t97
      t6673 = t926*t635*t2064
      t6678 = t97*t5449
      t6681 = 256*t6625*t4717-192*t2025*t93*t2002-8*t9*t6514*t1445+24*t3
     #87*t4860-48*t387*t3084-64*t6647+4*t4909*t4979+8*t2100*t1171-32*t46
     #08*t4883+24*t6657-t6660+576*t4854*t5722+144*t6663*t9-8*t6666*t2311
     #-8*t6666*t5744*t1233+48*t6673+32*t5517*t1123*t4940-16*t6678*t2373
      t6682 = t152*t4740
      t6683 = 4*t6682
      t6691 = t439*t5594
      t6694 = t141*t4750
      t6696 = t262*t37
      t6698 = 8*t6696*t2761
      t6702 = t6452*t582
      t6705 = t57*t1190
      t6710 = t152*t159
      t6714 = t53*t504
      t6730 = -t6683+12*t5494*t2270+4*t5494*t2349+72*t262*t4445*t5217-8*
     #t143*t6691-4*t6694-t6698-8*t5494*t5712*t161+8*t6666*t6702-2*t6705*
     #t5138-2*t6705*t2274+4*t57*t6710*t812+4*t6714*t225*t2761-24*t421*t2
     #267-16*t6678*t6311-16*t97*t6710*t2549*t656-64*t447*t2337-48*t4816*
     #t1129
      t6754 = t53*t1578
      t6755 = t6754*t832
      t6763 = t57*t5449
      t6774 = t4939*t4910
      t6776 = 8*t689*t5325*t5745+48*t97*t1578*t821-8*t53*t1211*t225*t162
     #0+32*t97*t5498*t1132+32*t4994*t6367+4*t6714*t5712*t1320-24*t421*t5
     #117+8*t6755-16*t4939*t5019+8*t5201*t1166-16*t5470*t2443-2*t6763*t2
     #332-2*t6763*t4948+8*t4962*t6357-4*t5848*t2293+24*t4939*t5025-8*t67
     #74
      t6792 = t5192*t856
      t6795 = t6183*t791
      t6798 = t3073*t116
      t6804 = t5657*t856
      t6805 = t4629*t31
      t6806 = t6*t6805
      t6817 = t57*t430
      t6824 = 6*t3201*t582
      t6825 = 96*t4939*t5536+24*t4939*t4979-48*t4939*t4982-8*t4939*t4913
     #+8*t4939*t4883-6*t3788*t1060-192*t3244*t386*t2830+32*t6792*t5181-4
     #*t85*t6795-48*t6798*t3246+72*t5727*t386*t1093+16*t6804*t6806+96*t5
     #04*t5506*t6571-48*t5517*t800*t2141+32*t5517*t2523-2*t6817*t4869*t5
     #82+36*t1060*t6121+t6824
      t6831 = t439*t1095
      t6834 = t4939*t5030
      t6840 = t4939*t6006
      t6842 = t280*t964
      t6852 = t386*t4838
      t6858 = t6*t614*t32
      t6863 = t2519*t57
      t6868 = t403*t5909
      t6871 = 24*t262*t959
      t6872 = 16*t4939*t5008+24*t4939*t4959+72*t4939*t6831-8*t6834+8*t49
     #39*t4968-4*t4962*t4936-24*t6840-192*t6842*t616-8*t5657*t2136*El*t6
     #05+32*t5657*t4136*t1004+24*t5727*t6852+72*t6798*t6381-24*t6804*t68
     #58-24*t6792*t1004-2*t6863*t6702-2*t6817*t2366+2*t6868-t6871
      t6876 = 3*t959*t152
      t6879 = t4589*t2851
      t6890 = t4962*t1623
      t6892 = t4589*t1677
      t6898 = t5677*Mzsq
      t6899 = t4939*t6898
      t6901 = t5677*sll
      t6902 = t4939*t6901
      t6908 = t53*sll
      t6910 = 2*t6908*t5571
      t6913 = 12*t6817*t2352+t6876+18*t4589*t1615+2*t6879-72*t4939*t1615
     #+96*t122*t3356-12*t4962*t1626-288*t128*t1592*t152+18*t6890+8*t6892
     #-24*t4939*t6203+48*t4939*t6200-8*t6899+16*t6902+6*t4589*t6203-12*t
     #4589*t6200-t6910+8*t4929*t1626
      t6921 = t751*t152
      t6928 = t1433*t4759
      t6933 = t5105*t1672
      t6938 = t1093*BCsThl
      t6944 = t4364*t603
      t6963 = -384*t4*t496*t102*t5197*t152+288*t5507*t102*t6921-96*t6583
     #*t5180*t791+t6928-8*t5612*t4683-4*t4909*t3369+6*t6933+32*t487*t335
     #6-2*t225*t2063*t6938-6*t4531*t5957*BCsThl+120*t5116*t6944-192*t513
     #1*t6452*sqq+96*t5131*t6452*Mzsq+6*t5301*t5022-128*t5113*t4936-24*t
     #5269*t4930+36*t5269*t4944+192*t5113*t5022
      t6968 = t159*t4628
      t6975 = t544*Mzsq
      t6980 = t161*t32
      t6981 = t6980*t2769
      t6983 = t4589*t6901
      t6985 = t4939*t1677
      t6987 = t4962*t1632
      t6989 = t4589*t6898
      t6991 = t4895*t1632
      t7001 = t196*t1095
      t7004 = -16*t4593*t5019+24*t4593*t4959-64*t5141*t6968*sqq+32*t5141
     #*t6968*Mzsq+2*t6975*t1691-96*t978*t6010-t6981-t6980*t4870-4*t6983-
     #32*t6985-4*t6987+2*t6989+16*t6991+128*t5113*t4768+24*t5269*t4963+8
     #*t4593*t4968+8*t4593*t4883+24*t4939*t7001
      t7008 = t196*t791
      t7014 = 8*t6696*t4838
      t7015 = t80*sll
      t7018 = t196*t2300
      t7023 = t196*t804
      t7036 = t4589*t4910
      t7040 = t57*t4688
      t7043 = t336*t1119
      t7045 = t336*t2264
      t7047 = 4*t4589*t5002-16*t4895*t7008+4*t4962*t7008-t7014+8*t7015*t
     #5005+8*t7015*t7018-16*t6372*t7008+24*t4895*t7023+384*t128*t3080-16
     #*t5269*t4926-48*t6372*t4933+16*t6372*t4930-6*t4962*t7023-6*t7036+8
     #*t4939*t7018+8*t7040*t4933-24*t7043-8*t7045
      t7058 = t57*t152
      t7069 = t439*t4688
      t7074 = t4589*t3422
      t7076 = t4589*t3015
      t7078 = t608*t152
      t7079 = t4589*t7078
      t7085 = t439*t152
      t7088 = -288*t4842*t5198+16*t6372*t4936-16*t6372*t4926+8*t7015*t49
     #68+24*t7015*t5025+4*t7058*t4910-8*t7040*t4930-4*t7058*t5025-128*t4
     #608*t4941-8*t4962*t4802+16*t4589*t7069-16*t4589*t6399-12*t7074-4*t
     #7076+6*t7079+48*t4939*t3422-72*t4939*t1575-8*t4962*t7085
      t7093 = t4962*t1821
      t7107 = 64*t504*t4834*t2063*t1078
      t7108 = t93*t4759
      t7113 = t7015*t4910
      t7119 = t7015*t5030
      t7127 = t608*t32
      t7132 = 18*t7093+18*t4589*t1575-64*t5061*t125+64*t5061*t1966+12*t6
     #84*t532*t4896-t7107-96*t6393*t7108-4*t7058*t4979-8*t7113-8*t7015*t
     #4760-6*t4589*t7001-8*t7119+4*t7058*t4913-2*t4589*t7018-6*t3624*t10
     #60+120*t5116*t7127+192*t128*t3379
      t7133 = t4939*t7078
      t7139 = t4962*t1276
      t7141 = t4895*t4802
      t7146 = 32*t4895*t4768
      t7147 = t4895*t1276
      t7149 = t5871*BCsThl
      t7152 = t4994*t46
      t7155 = t5871*t417
      t7170 = t196*t152
      t7173 = -24*t7133-32*t4895*t6357-96*t4895*t7085+6*t7139+32*t7141-7
     #2*t4895*t1821+t7146-24*t7147+288*t7149*t1320+120*t7152*t1221+96*t7
     #155*t1292+288*t5131*t1975-8*t5494*t2537+12*t5494*t2463-192*t5131*t
     #2118+96*t5141*t1375-64*t5141*t2940+8*t4962*t7170
      t7175 = t5116*t3107
      t7181 = t4951*t3007
      t7190 = t796*t4688*CsThr
      t7193 = t145*t5740
      t7196 = 16*t4895*t4936
      t7199 = t46*t905
      t7200 = t7199*Mzsq
      t7205 = t5120*t2123
      t7207 = t386*t545
      t7219 = -80*t7175-32*t4895*t7170-12*t4919*t1300+8*t7181+24*t7015*t
     #4979-8*t7015*t4913+8*t7015*t4883+12*t145*t7190+2*t7193+t7196-24*t4
     #895*t5022-12*t4951*t7200-168*t5120*t6975+112*t7205+32*t6798*t7207+
     #128*t4*t3894*t4799+96*t504*t4553*t2063*t1046-4*t6863*t2311
      t7220 = t4939*t3015
      t7230 = t4895*t4930
      t7234 = t234*t2300
      t7235 = t145*t7234
      t7239 = t145*t5678
      t7241 = t145*t6795
      t7251 = t145*t5591
      t7253 = t5178*t404
      t7254 = t93*t964
      t7259 = 16*t7220+72*t4895*t4987-16*t4895*t4926+24*t4895*t4963-48*t
     #4895*t4933+16*t7230-2*t4589*t5005+6*t7235+4*t2011*t2100+2*t7239-4*
     #t7241-12*t145*t4754-8*t145*t5595+6*t4589*t6006+2*t4589*t5025-12*t7
     #251+144*t7253*t7254+144*t7253*t951
      t7270 = t6*t1002*Mzsq
      t7273 = t386*t964
      t7277 = t93*t32
      t7278 = t404*t7277
      t7287 = t6146*t582
      t7292 = t4589*t5126
      t7296 = t4939*t2346
      t7298 = t4939*t5126
      t7300 = t4589*t5030
      t7302 = t6863*t2316
      t7308 = -96*t7253*t3227-96*t5727*t6587+16*t5192*t3444*t3227-24*t67
     #92*t7270+72*t6798*t7273-24*t5657*t831*t7278+256*t21*t6*t5180*t152-
     #4*t4589*t4941+6*t145*t7287+8*t145*t6691-6*t7292-6*t4589*t6831+32*t
     #7296-24*t7298+2*t7300-2*t7302-384*t5507*t5085-4*t4909*t6203
      t7309 = t6514*Mzsq
      t7314 = t1635*t152
      t7315 = t4589*t7314
      t7319 = t4939*t2851
      t7323 = t4589*t3369
      t7333 = t262*t2622
      t7338 = t5657*t943*El
      t7347 = t423*t804
      t7352 = t4591*t5030
      t7354 = -4*t2113*t7309+8*t5877*t6099+6*t7315-24*t4939*t3369-8*t731
     #9-96*t128*t3356+6*t7323+16*t6428*t5197*t6453+12*t3868*t5554-24*t64
     #62*t5197*t1060-8*t7333*t1074*t1836-16*t7338*t6448-80*t7152*t821+48
     #*t4593*t4941-32*t4593*t4916+48*t4593*t7347-24*t4593*t5030+2*t7352
      t7364 = t4593*t4910
      t7380 = t246*t4628
      t7383 = t6908*t2805
      t7384 = 2*t7383
      t7391 = t4591*t4916
      t7395 = 2*t4591*t4760-4*t5301*t4936-4*t4591*t4941+6*t4591*t6006+12
     #*t7364-4*t4593*t4979+8*t4593*t4982+12*t4593*t4913-24*t4593*t5008+4
     #*t6451*t6805*t1060-24*t6428*t5216*t1836-80*t5116*t7380-t7384+192*t
     #4608*t6006+64*t4608*t5019-4*t141*t4865+8*t7391-12*t4591*t7347
      t7415 = 32*t6696*t1672
      t7421 = t4895*t1636
      t7423 = t145*t5858
      t7425 = t4962*t1636
      t7431 = 16*t3073*t37*t1620
      t7434 = -32*t4608*t4968-96*t4608*t4959+24*t4593*t5126+8*t5269*t493
     #3-12*t5269*t4987-8*t4593*t5536-4*t4593*t5025-72*t5269*t5022-24*t45
     #93*t4760+t7415+96*t5131*t2842+16*t4763*t1371*t4628-24*t7421-2*t742
     #3+6*t7425+288*t7149*t1859+t7431+96*t7155*t1874
      t7438 = t4939*t7314
      t7446 = t145*t2810
      t7448 = t145*t2830
      t7453 = t423*t4940
      t7457 = 64*t4608*t5030
      t7472 = t7015*t6898
      t7476 = -24*t7438+48*t4895*t1626-72*t4895*t1623+192*t538*t2830-2*t
     #7446+4*t7448+8*t4531*t1433*t800+8*t141*t7453+t7457-12*t141*t4874-2
     #56*t4133*t612*t4717-64*t5586*t4796+64*t123*t117*t5508-24*t7015*t62
     #03-12*t4919*t1599-8*t7472-168*t5123*t1789
      t7485 = t280*t612
      t7507 = t6666*t2316
      t7511 = t7015*t2851
      t7513 = t3006*t4628
      t7518 = 120*t7152*t1688-64*t7155*t1881-192*t7149*t1865+288*t5131*t
     #1789+48*t7485*t7278+12*t5657*t922*El*t5709-96*t123*t117*t6921+128*
     #t5586*t4779+192*t5547*t4790+16*t446*t2337-32*t143*t7453-12*t143*t7
     #190+96*t5141*t1763+8*t7507-24*t7015*t3369-8*t7511+8*t4951*t7513-24
     #*t4995*t1763
      t7520 = t7199*t603
      t7525 = t93*CsThr*t1862
      t7528 = t262*t780
      t7540 = t93*t1060
      t7570 = -12*t4951*t7520-8*t5192*t112*t7525+32*t7528*t5799*t1836-16
     #8*t5123*t1975+128*t6620*t4815*t2286+16*t7333*t386*t125+12*t6458*t7
     #540-32*t123*t4831+64*t5547*t4782-24*t7528*t7207-8*t165*t5791-8*t52
     #74*t5390*t615+32*t5321*t375*t998-24*t5321*t5315*t636-256*t5820*t99
     #9-96*t4183*t117*t751*sll-128*t5547*t4793+64*t5547*t4796
      t7572 = 64*t5871*t791
      t7583 = t141*t3578
      t7603 = t3041*t4276
      t7615 = t7572+192*t5547*t4799-24*t5274*t5315*t947+16*t5321*t5275*t
     #615+16*t4763*t2922-2*t7583-16*t4593*t7069-256*t5113*t4802-48*t5269
     #*t4768+128*t5113*t6357-64*t4593*t5673+16*t5269*t7085+96*t5269*t480
     #2-16*t5192*t905*El*t5709+12*t7603*t174-32*t80*t430*t2286-24*t5657*
     #t922*t9+16*t6349*t93*t6514
      t7627 = t926*t635*t2893
      t7629 = t280*t2158
      t7638 = t280*t621
      t7641 = Er*El
      t7647 = t5178*t621
      t7656 = t9*t1433*t2075
      t7658 = t387*t3046
      t7661 = t896*t934*t2079
      t7668 = 192*t6380*t386*t4837+3*t141*t6408-256*t6620*t211*t2346+48*
     #t7627-32*t7629*t917+12*t6451*t7*t905*t1836-32*t6264*t6575+48*t7638
     #*t6858+12*t3041*t7641*t4722+16*t7629*t1078+48*t7647*t7270-64*t7485
     #*t3337-64*t5178*t48*t7278+16*t7656-48*t7658+24*t7661+96*t689*t211*
     #t2346-32*t7638*t6806
      t7679 = 32*t1713*t2300
      t7690 = Er*t159
      t7694 = 32*t5936
      t7695 = 128*t5939
      t7696 = t4994*t32
      t7697 = 96*t7696
      t7698 = t905*t4628
      t7700 = 2*t7698*sll
      t7701 = t7698*Mzsq
      t7703 = 3*t4741*sll
      t7704 = t5931*Mzsq
      t7705 = 3*t7704
      t7707 = 8*t2141*t905
      t7708 = 48*t7647*t1004-288*t6393*t951+16*t7638*t4839-96*t6304*t698
     #*t816-t7679+64*t5517*t1123*t1677-8*t7603*t616+16*t5192*t111*t4697-
     #24*t6343*t4166-32*t689*t7690*t5745+t7694+t7695-t7697-t7700+t7701+t
     #7703+t7705-t7707
      t7711 = 16*t262*t6514
      t7712 = t906*t2063
      t7713 = 2*t7712
      t7714 = t959*t2063
      t7716 = t931*t3031*t2075
      t7720 = t931*t671*t403*t2063
      t7722 = t1233*t905
      t7723 = 3*t7722
      t7731 = t6372*t1632
      t7736 = 6*t1707*t1233
      t7741 = t464*t651
      t7742 = t5445*t7
      t7743 = t5316*sqq
      t7752 = t97*t504
      t7758 = t7711-t7713+t7714-8*t7716-8*t7720+t7723+16*t1452*t5572+32*
     #t4816*t1132-32*t80*t2939*t792+16*t7731-32*t7647*t5684+t7736+8*t57*
     #t31*t417*t2286-8*t7741*t7742*t7743-96*t35*t116*t6389-8*t4365*t5536
     #-40*t7752*t2277-96*t6620*t404*t7108
      t7766 = 4*t1707*t2469
      t7767 = t464*t3196
      t7784 = t12*t152
      t7786 = t7742*t974
      t7803 = t614*SnThr*t47*Mzsq
      t7806 = t628*t999
      t7816 = 4*t6014*t1211-24*t4640*t614*t900*t47-t7766+16*t7767*t5276*
     #t1314+128*t6620*t5325*t792+128*t2025*t93*t6091+16*t7767*t94*t2769-
     #8*t7741*t2089*t7*t974+16*t7784*t651*t7786-32*t7784*t3196*t94*t1862
     #+8*t27*t2469*t616+8*t27*t173*t4631+16*t27*t4700*El*t7803+24*t7806+
     #16*t4623*t174+48*t4653*t5442+16*t6170*t7803-8*t5441*t4631
      A1e0c12p1t1n = 24*(t6633+t6582+t6285+t5642+t5363+t6236+t7004+t6963
     #+t6913+t5310+t6872+t5045+t5258+t6825+t5204+t4998+t7047+t7308+t7259
     #+t7219+t5964+t7173+t7132+t7088+t5919+t4947+t5583+t5462+t5520+t6681
     #+t5411+t5866+t7615+t7570+t7518+t6730+t7476+t7434+t7395+t7354+t4757
     #+t4693+t6136+t6094+t6530+t6050+t6013+t7816+t7758+t7708+t7668+t6492
     #+t6438+t4823+t4886+t5819+t5146+t5762+t6776+t6189+t6392+t5702+t5095
     #+t6336)*t2054*t580*t567*t151*t149*t569*t2052*t594*t156*t572
      A1e0c1p1t2n = A1e0c1p1t1n
      A1e0c2p1t2n = A1e0c2p1t1n
      A1e0c3p1t2n = A1e0c3p1t1n
      A1e0c4p1t2n = A1e0c4p1t1n
      A1e0c5p1t2n = A1e0c5p1t1n
      A1e0c6p1t2n = A1e0c6p1t1n
      A1e0c7p1t2n = A1e0c7p1t1n
      A1e0c8p1t2n = A1e0c8p1t1n
      A1e0c9p1t2n = A1e0c9p1t1n
      A1e0c10p1t2n = A1e0c10p1t1n
      A1e0c11p1t2n = A1e0c11p1t1n
      A1e0c12p1t2n = A1e0c12p1t1n
      t7837 = 24*t143
      t7838 = 8*t145
      t7840 = 32*t80
      t7841 = t78-24*t123+16*t131-t127+t84-32*t437-t142-t34-32*t536+t783
     #7-t115+t32+t7838+24*t88-t7840
      A2e0c1p1t1n = 96*Mzsq*t7841*t149*t588
      t7863 = -2*t322*t290+8*t160*t171+4*t170*t161-4*t5657*El*t167+2*t53
     #*t58+t1279-8*t4129*t167-4*t4951-8*t2509+16*t293*t225-8*t385
      t7868 = 4*t54
      t7871 = t2675*Er
      t7877 = t192*BSnThl
      t7882 = 16*t343
      t7883 = -4*t293*t161+2*t315*t290+t7868+4*t262+8*t97-2*t7871-16*Eq*
     #t37*El*t15+4*t7877*t1493-t191+8*t4129*t15-t7882
      t7885 = 1+BCsThl
      t7886 = t7885**2
      t7887 = 1/t7886
      t7890 = -1+BCsThl
      t7891 = t7890**2
      t7892 = 1/t7891
      t7894 = t598*t594
      t7895 = t7894*t149
      A2e0c2p1t1n = 96*t596*(t7863+t7883)*t7887*t151*t7892*t7895
      t7898 = t2063*Mzsq
      t7899 = t7898*t31
      t7900 = t7899*t1
      t7903 = t7898*Eq
      t7904 = t7903*Er
      t7907 = t43*t40
      t7908 = t7907*t46
      t7917 = t47*t7898
      t7918 = t7917*sqq
      t7925 = t41*t37
      t7928 = t2*t19
      t7929 = t7928*t40
      t7930 = t7929*BCsThl
      t7933 = t7929*t46
      t7940 = t7907*BCsThl
      t7943 = t246*t2
      t7944 = t26**2
      t7945 = t7944*t46
      t7949 = t67*sqq
      t7954 = 8*t7900*t978-4*t7904*t5868+32*t7908*t1221+8*t1352*t6975-4*
     #t137*t1589-16*t1352*t6980+16*t90*t7918+4*t137*t1300-8*t137*t1305-8
     #0*t7925*t1328+128*t7930*t1311-64*t7933*t1276+128*t7933*t1279-256*t
     #7930*t298+128*t7940*t1347+32*t7943*t7945*sll+16*t7949*t7944*sll*Mz
     #sq
      t7979 = t7917*sll
      t7982 = t7898*sll
      t7983 = t7982*sqq
      t7986 = CsThr*t7898
      t7993 = t2141*t2063
      t7996 = t7917*t31
      t7999 = -64*t7933*t1273-32*t90*t1355-4*t90*t1589-32*t90*t2152-4*t9
     #0*t1298+16*t90*t2171-32*t90*t2979+32*t137*t2893+4*t137*t1298-64*t7
     #908*t338-32*t90*t2893+8*t90*t5776+16*t90*t7979-8*t1155*t7983-16*t7
     #90*t7986*t31-16*t7943*t7945*Mzsq-8*t1155*t7993+12*t1155*t7996
      t8001 = t370*t26
      t8002 = t8001*t46
      t8009 = t4*t7898
      t8012 = t7898*t2
      t8013 = t8012*sqq
      t8028 = t54*t7898
      t8033 = t8001*t152
      t8037 = t4688*sqq*Mzsq
      t8045 = 32*t8002*t6183*sqq+64*t1308*t7986*sqq-24*t8009*t1331+64*t5
     #05*t8013-32*t978*t8013-8*t557*t1300+16*t557*t1305-4*t1155*t2000+24
     #*t1155*t2209+8*t1155*t2011-16*t1155*t8028-8*t1155*t2211-16*t8033*t
     #612-t1155*t8037+32*t8033*t1060+128*t7930*t1320-64*t7940*t1975
      t8054 = 16*t4879
      t8067 = t4648-128*t4806-8*t1155*t7917*t804+16*t27*BCsThl*t3196*t31
     #98+t8054+8*t4955+4*t5013+t5032-16*t5038+8*t5043-128*t5068-24*t5080
     #-2*t5107+64*t5121-128*t5129+8*t5168-32*t5172+256*t5202
      t8084 = t7898*t905
      t8088 = 16*t8001*t1095
      t8089 = -t5229+8*t5240-16*t5249-8*t5251-64*t5256+8*t5296+32*t5337+
     #2*t5378-128*t5382-16*t5384-24*t5386+128*t5421+8*t5423-16*t5453-16*
     #t5455+2*t8084*t4-t8088
      t8091 = 8*t80*t1093
      t8093 = 16*t8001*t816
      t8094 = t7898*BSnThl
      t8095 = t8094*t8
      t8098 = t7907*t116
      t8102 = 128*t7929*t791
      t8104 = 64*t7929*t804
      t8106 = 64*t7929*t656
      t8107 = t7898*t43
      t8109 = 128*t8107*t4154
      t8111 = 2*t57*t906
      t8114 = t7898*t152
      t8118 = 32*t8001*t939
      t8122 = 32*t8001*t208
      t8126 = t7898*t100
      t8127 = t8126*t677
      t8130 = -t8091+t8093-32*t8095*t880+128*t8098*t6346+t8102-t8104-t81
     #06-t8109-t8111+256*t5587-32*t5589-4*t57*t8114+t8118-64*t7940*t1173
     #-t8122+12*t57*t7899-8*t5622-32*t8127*t866
      t8140 = 16*t80*t1433
      t8144 = 16*t8001*t964
      t8145 = t27*t6457
      t8148 = t7899*t7877
      t8152 = 64*t8107*t4114
      t8155 = t7898*t603
      t8156 = t8155*sll
      t8160 = t7898*t37*t603
      t8163 = t8155*sqq
      t8166 = t32*t7898
      t8172 = 32*t8001*t1620
      t8175 = 64*t8127*t871+80*t90*t1176-4*t90*t58*t1233+t8140-8*t41*t14
     #36+t8144+8*t8145*t6459+16*t8148*t6459+t8152-4*t57*t5942+8*t57*t815
     #6-32*t57*t8160-24*t57*t8163-4*t8166*t1836+t57*t1635*t905-t8172+4*t
     #57*t664
      t8177 = 32*t8001*t1677
      t8182 = t7917*t603
      t8186 = 64*t7907*t1620
      t8192 = t7917*t152
      t8215 = t8177-16*t41*t2761-8*t41*t6944+32*t41*t8182-t8186-2*t57*t3
     #014*t905+4*t57*t7996-4*t57*t8192-24*t57*t2209-8*t57*t2011+8*t57*t2
     #211+4*t57*t2000+32*t41*t2893+16*t41*t4610+32*t41*t2152-8*t41*t5776
     #-8*t41*t7127-16*t41*t7979
      t8223 = t46*t7898
      t8224 = t8223*sqq
      t8229 = t8223*sll
      t8235 = 4*t7903*t4248
      t8240 = 64*t8107*t317
      t8241 = t7898*t19
      t8242 = t8241*Er
      t8247 = Er*t32
      t8248 = t7903*t8247
      t8260 = -16*t41*t7918+4*t41*t1298+16*t41*t8224+8*t41*t1585+16*t41*
     #t8229+8*t8166*t2113+t8235-16*t7949*t7945*t582+t8240-64*t8242*t5062
     #-32*t90*t2337+2*t8248+8*t90*t844+16*t90*t828-8*t137*t844-64*t7908*
     #t821+16*t4080*t5030
      t8270 = t532*t31
      t8273 = t532*t32
      t8276 = t2640*t318
      t8277 = 128*t8276
      t8282 = t8223*t31
      t8296 = t243*t7898
      t8301 = 32*t7908*t824+32*t41*t8155-4*t5687+128*t5690+2*t7903*t5899
     #*sqq-8*t7903*t8270+4*t2075*t8273-t8277+t57*t8037+2*t57*t4694*t31-1
     #2*t57*t8282+4*t57*t8223*t152+t57*t47*t905*Mzsq+64*t5697+8*t57*t301
     #9+16*t57*t8028+16*t57*t8296+8*t57*t7983
      t8305 = t7899*t2
      t8308 = t8166*t1
      t8327 = 8*t41*t1433
      t8332 = 32*t7907*t964
      t8334 = 128*t7929*t191
      t8337 = 64*t7929*t582
      t8338 = Er*t100
      t8339 = t8241*t8338
      t8342 = 8*t57*t7993+64*t4*t8305-24*t4*t8308-32*t48*t8305-4*t90*t16
     #35*t939+16*t90*t2267-64*t7940*t839+128*t7940*t812-8*t137*t1119-64*
     #t90*t832+t8327-8*t1155*t6183*t2083+t8332+t8334+16*t5710-t8337+64*t
     #8339*t770
      t8347 = 32*t8241*t318
      t8348 = t8166*t505
      t8357 = 64*t7907*t208
      t8358 = t8094*El
      t8364 = t8241*BSnThl
      t8371 = 64*t7929*t596
      t8372 = t48*t816
      t8375 = t423*t656
      t8378 = t1358*t2256
      t8380 = t7917*t656
      t8383 = 2*t7903*t211*t32-t8347+2*t8348-4*t5732+16*t5741+32*t41*t13
     #55-16*t41*t2171-t8357-t5748+16*t8358*t685-32*t8126*t98*t781-64*t83
     #64*t738+4*t7904*t1115*t191-t8371+8*t7904*t8372-16*t7904*t8375-4*t8
     #378-24*t1155*t8380
      t8394 = 4*t41*t943
      t8396 = 32*t7907*t1328
      t8397 = t2*t7944
      t8411 = t7986*t656
      t8420 = t8001*BCsThl
      t8427 = t8001*sll
      t8430 = -8*t1155*t663*t1233-8*t139*t844+4*t1155*t2256-t8394+t8396+
     #16*t246*t8397*Mzsq-32*t246*t8397*sll-16*t41*t2810-64*t41*t2781+16*
     #t246*t8397*t603+32*t790*t8411+32*t7949*t7945*t791-16*t7949*t7945*t
     #804+64*t8420*t812-32*t8002*t821-64*t1155*t2293+32*t8427*t4936
      t8433 = t43*t26
      t8434 = t8433*t1060
      t8444 = t1635*t3018
      t8451 = t7898*sqq
      t8456 = t8223*t603
      t8460 = 32*t7907*t1609
      t8464 = 64*t8241*t4357
      t8465 = t4154*t31
      t8466 = t7903*t8465
      t8468 = t2640*t4357
      t8469 = 64*t8468
      t8470 = -16*t41*t7982+32*t8434*t999+64*t8339*t848-16*t5780+t5786+3
     #2*t1155*t1188-64*t8427*t4930+4*t1155*t8444+32*t7908*t1179-32*t90*t
     #2470-16*t41*t8451-16*t41*t1781-32*t41*t8456+t8460+4*t41*t1670+t846
     #4-4*t8466+t8469
      t8472 = t2640*t4341
      t8473 = 64*t8472
      t8474 = t4753*sqq
      t8477 = 8*t41*t1672
      t8481 = t4688*t31
      t8486 = t8433*t1041
      t8491 = t8433*t4664
      t8502 = t27*t6342
      t8505 = t8433*t1049
      t8508 = t80*t116
      t8511 = t8433*t404
      t8517 = t8473+t57*t8474+t8477-8*t8148*t93*t505-2*t57*t8481-32*t843
     #4*t917-32*t8486*t1004-16*t6170*t7525-32*t8491*t3915-32*t8491*t4166
     #+64*t8434*t1078+64*t27*t2469*El*t6448-8*t8502*t4166+64*t8505*t961+
     #16*t8508*t6381+32*t8511*t7254-64*t7907*t582*t387
      t8521 = 4*t5885
      t8522 = 64*t5897
      t8523 = 32*t5920
      t8526 = 8*t5955
      t8529 = 32*t5984
      t8531 = 2*t6029
      t8532 = 8*t6062
      t8533 = 2*t6064
      t8535 = -64*t7907*t3884*t3879-t8521+t8522-t8523-4*t5929-24*t5947-t
     #5953+t8526+2*t5959-2*t5962-t8529+2*t6019+t8531+t6057+t8532-t8533-t
     #6069+12*t6089
      t8545 = t80*t152*t98
      t8550 = t370*t40
      t8556 = t8012*t5399
      t8563 = t7898*t1
      t8564 = t8563*t5399
      t8570 = t8166*t1459
      t8580 = -t6095-128*t8339*t679+32*t8486*t4675-16*t8545*t3879+48*t61
     #16+8*t6132-128*t8550*t173*t4166+2*t6145+2*t6153+128*t8556*t718+32*
     #t41*t2979+128*t8556*t767-32*t8564*t738-256*t2625*t5399*t1475-8*t85
     #70*t8*SnThr*CsThr*BCsThl+128*t2625*t4114*t616+64*t6323
      t8586 = 32*t6397
      t8590 = 16*t6475
      t8597 = 32*t6499
      t8598 = 2*t6501
      t8599 = 16*t6502
      t8600 = 32*t6504
      t8601 = 64*t6506
      t8602 = -24*t6326-32*t6328+8*t6370-256*t8556*t786-t8586-4*t6413-8*
     #t6449+4*t6473-t8590+64*t8564*t975-64*t8434*t4717+8*t8145*t7540-t85
     #97-t8598+t8599-t8600+t8601-t6562
      t8606 = t8550*t404
      t8618 = t2*t40
      t8619 = t8618*t404
      t8632 = t27*t1049
      t8642 = -32*t6566+t6602-4*t6612+256*t8606*t6260+32*t8511*t951-32*t
     #8508*t3246-64*t8511*t3227-128*t8550*t582*t9-64*t8619*t3227+32*t861
     #8*t682*t1475+32*t8619*t1453+64*t8511*t940+8*t73*t844-128*t6647-8*t
     #8632*t1088*t911-32*t8508*t386*t208+32*t8508*t386*t939
      t8662 = t8550*t1041
      t8665 = t57*t922
      t8668 = t8618*t1060
      t8687 = 32*t80*t3884*t102*t4825+16*t8632*t8*t900*t31-8*t8632*t108*
     #t1298+16*t655*t93*t1433-16*t8632*t8*t671*t32-32*t6657-256*t8662*t4
     #697+t8665-64*t8511*t1015-32*t8668*t917+128*t8662*t5722+128*t8662*t
     #605-8*t8502*t3915+128*t8662*t4722+48*t1155*t2488-16*t790*t7986*t19
     #1-64*t6673+32*t8242*t1195
      t8724 = 2*t8084
      t8727 = t2835*t1781
      t8729 = 16*t1155*t7917*t191-48*t1155*t1405+20*t1155*t2984+4*t2835*
     #t423*t939+8*t7904*t423*t191-4*t7904*t48*t208+16*t6694-8*t1155*t301
     #9-256*t7930*t1865-4*t7904*t1612+32*t7925*t1677-16*t1155*t8296-128*
     #t2613*t8338*t677*t3483+64*t2613*t4114*t1046-t8724+8*t1697*t8156-16
     #*t8727
      t8730 = t1358*t664
      t8732 = t8107*Er
      t8757 = t7944*t603
      t8768 = 4*t8730+256*t8732*t1836-128*t8732*t1757-4*t7904*t1609+16*t
     #7904*t1737+4*t6755+64*t6774+4*t1155*t5942-8*t1155*t8156+32*t1155*t
     #8160-8*t1697*t8163+24*t1155*t8163-t1155*t8474-16*t7943*t7945*t603+
     #16*t7949*t8757*sll-32*t139*t1672-2*t2330*t959*BCsThl+8*t139*t2830
      t8776 = t57*t2063
      t8806 = 32*t8033*t1836-16*t8033*t1757-4*t139*t1670+64*t8776*t1737-
     #4*t1155*t664-32*t8420*t1789+16*t8002*t1688-8*t557*t1599+16*t1308*t
     #5465-128*t1308*t7986*t603-128*t473*t8163+48*t8009*t1754+64*t1327*t
     #8163-16*t7925*t1609+16*t90*t2761+16*t90*t1781+4*t137*t1670
      t8824 = t7929*sll
      t8843 = 16*t137*t1672-16*t137*t2810-4*t1697*t5942+128*t555*t8182-6
     #4*t6834+8*t6840-16*t7900*t1757-32*t7949*t8757*Mzsq+128*t7933*t1632
     #+128*t8824*t1836-64*t8824*t1757+64*t90*t2781+16*t90*t2810-8*t4047*
     #t1328+48*t139*t1589-2*t1257*t3512+4*t336*t3512-4*t90*t1667
      t8882 = -4*t90*t1670-32*t90*t8182-16*t137*t2761-64*t7929*sqq*t1757
     #+64*t7925*t1620+8*t1352*t5465-64*t7940*t1789+32*t7908*t1688+4*t137
     #*t1599+128*t7930*t1859-32*t1308*t6980+64*t8420*t1347-32*t8420*t197
     #5-64*t555*t7918-64*t555*t7979+16*t8002*t1221+16*t1308*t6975
      t8884 = t423*t32
      t8903 = 4*t6879
      t8906 = 16*t6985
      t8907 = t8012*t5395
      t8914 = -4*t6868-4*t2835*t8884+8*t7904*t1359+8*t1697*t3019+4*t1155
     #*t8192-32*t8002*t338-20*t139*t1298-48*t139*t2893+64*t8732*t978-2*t
     #7904*t1579+t8903+t6890+4*t6892+2*t6981-t8906-128*t8907*t174+256*t8
     #907*t616-32*t8564*t1475
      t8916 = t27*t8009
      t8924 = t4*t8364
      t8927 = t27*t2100
      t8938 = 16*t8001*t1233
      t8944 = Er*t98
      t8956 = 64*t8916*t16+64*t8916*t174-2*t6987-128*t8916*t616+64*t8924
     #*t1475+64*t8927*t1078-128*t8927*t1082+64*t4*t8241*t917-128*t8924*t
     #975-t8938-64*t8098*t6377+128*t2625*t7641*t605+64*t2613*t8944*t1034
     #+256*t8606*t757-64*t8098*t6270-128*t8606*t657-128*t8606*t6386
      t8960 = t8241*t404
      t8963 = Er*t31
      t8965 = 32*t8241*t8963
      t8974 = 8*t7076
      t8977 = t8618*t1049
      t8984 = -16*t7036-2*t7043+2*t7045-64*t8960*t757-t8965+32*t8619*t95
     #1+128*t8098*t386*t191-64*t8619*t1015+24*t7074-t8974-t7079+t7093+t7
     #107+64*t7113-64*t7119-32*t8977*t1090+32*t8619*t7254+64*t8668*t1078
      t9011 = -32*t8668*t1082+64*t8977*t961-256*t8550*t1049*t1475-128*t8
     #907*t16+128*t8662*t4730-32*t8618*t1041*t1004-16*t7141+t7146+64*t71
     #75-16*t7181+t7196-128*t7205-32*t7230-8*t7235+2*t7251-8*t7904*t352+
     #2*t7292
      t9025 = t2835*t1355
      t9027 = t1358*t2011
      t9029 = t1358*t2211
      t9039 = 4*t7383
      t9040 = -16*t7298+12*t7300-16*t7302-4*t2835*t2161-128*t8732*t122+3
     #2*t8242*t125-128*t8732*t505+64*t8732*t128+16*t9025-8*t9027+8*t9029
     #-t7315-12*t7323+64*t1308*t7986*sll-64*t8420*t1230-4*t7352+32*t7364
     #+t9039
      t9049 = 32*t7907*t816
      t9057 = t7982*t32
      t9058 = 2*t9057
      t9059 = t8155*t32
      t9060 = 4*t9059
      t9061 = 32*t8305
      t9062 = 16*t8308
      t9063 = -4*t7391+8*t7446-2*t7448+t7457+32*t7507+128*t8960*t657+t90
     #49+8*t65*t844+4*t7583-64*t7627-32*t7656+64*t7658-32*t7661-t9058+t9
     #060-t9061+t9062
      t9064 = 4*t7712
      t9065 = 2*t7714
      t9069 = 8*t80*t943
      t9071 = 16*t7982*t33
      t9073 = 32*t7982*t35
      t9088 = t57*t959
      t9090 = 4*t41*t1093
      t9092 = 32*t8155*t33
      t9095 = 64*t8012*t656
      t9096 = t9064-t9065+16*t7716+16*t7720-t9069+t9071-t9073+16*t8508*t
     #7273-16*t8545*t102*t4261-4*t137*t1667-64*t8505*t912+64*t8433*t582*
     #t616+16*t8095*t875+t9088-t9090-t9092+8*t7806+t9095
      t9104 = t2054*t594*t580
      t9106 = 1/t2063
      A2e0c3p1t1n = -48*(t9063+t9040+t9011+t8843+t8882+t8914+t8956+t8984
     #+t8806+t8301+t8517+t8470+t9096+t8430+t8602+t8642+t7999+t8687+t8045
     #+t8729+t8067+t8768+t7954+t8383+t8580+t8535+t8260+t8215+t8175+t8130
     #+t8342+t8089)*t2052*t9104*t9106*t151*t2059
      t9122 = t40*CsThr
      t9132 = 24*t3869-t3893-t3905+t3907+24*t3926-32*t27*Mzsq*t9-32*t27*
     #t603*t9+64*t476*t6*t623-32*t181*t9122*t603-32*t181*t9122*Mzsq+64*t
     #181*t9122*sll
      t9134 = 32*t80*Mzsq
      t9136 = 16*t3935
      t9137 = 48*t3937
      t9138 = 24*t3947
      t9139 = 8*t3964
      t9140 = 24*t3966
      t9141 = -t9134+t943-48*t3928-t9136-t9137+t3940-t3946+t9138+t9139+t
     #9140-t3969
      t9143 = 8*t3970
      t9150 = t9143+16*t3976-t3979+t3981-24*t3982-t3990+16*t3991-24*t399
     #3-32*t3995+48*t4001-t4005
      t9151 = 64*t7015
      t9153 = 32*t80*t603
      t9154 = t4013+t4020+t4036-t4053-t4060+t4071+t1093+t4096+t9151-t409
     #7-t4099-t9153
      A2e0c4p1t1n = 24*(t9132+t9141+t9150+t9154)*t149*t151*t4110
      t9164 = t5617*Er
      t9176 = t671*t4265
      t9186 = t93*t3902
      t9192 = 2*t316*t4277-2*t323*t4277-2*t9164*t1877+2*t9164*t1891-8*t1
     #2*t5445*t167+4*t6190+t1179-8*t12*t4433*t167-4*t931*t9176+t169-2*t5
     #3*t596-4*t874*t9176+8*t194-16*t4136*t3891-16*t3444*t9186+4*t2074*t
     #614*t5657
      t9212 = t831*t1
      t9215 = t2609*sll
      t9222 = t4154*t47
      t9232 = 4*t683*t93*t5657+8*t12*t603*t9+8*t12*Mzsq*t9+16*t293*t4531
     #+8*t4834*t7690*t603-8*t97*t1836+8*t160*t4313-8*t1085+16*t4*t9212-4
     #*t4*t9215+4*t196*t9215+2*t53*t1273+2*t2675*t9222+4*t170*t1856-4*t2
     #93*t1856+8*t2774-8*t97*t1635
      t9238 = 8*t38
      t9248 = 4*t1737
      t9249 = 4*t1958
      t9250 = -4*t804*t249+t824-2*t2675*t4154-t9238+2*t208-4*t247+8*t250
     #+t259-8*t306*t4276+32*t271-8*t274-t89-4*t2609*t3014+t9248-t1609+t9
     #249
      t9253 = 16*t343*t603
      t9254 = 16*t9212
      t9260 = 32*t334
      t9265 = 8*t1194-t9253-t9254-t1328-t289-32*t295+8*t299-8*t301-4*t31
     #9+4*t324+t9260-2*t338+4*t9215+t349-8*t350+4*t141-t144
      A2e0c5p1t1n = 24*(t9192+t9232+t9250+t9265)*t7887*t7892*t4620*t571
      t9291 = t8166*t46
      t9295 = t3094*t46
      t9304 = t57*t504
      t9305 = t2273*t4837
      t9308 = t57*t1211
      t9311 = Er*t7898
      t9321 = t246*t2300
      t9328 = 3*t6908*t9295+16*t7040*t1359-4*t5105*t7200-8*t7040*t1582-8
     #*t9304*t9305+16*t9308*t2274+32*t4819*t9311*t603+96*t280*t4276*t100
     #*t677*t678*t31-28*t235*t9321+8*t235*t8380+16*t7752*t9305
      t9330 = t4276*BSnThl
      t9347 = t8338*t98
      t9356 = t811*t1093
      t9362 = t48*t964
      t9374 = t640*t38
      t9377 = t57*t1578
      t9382 = t2273*t816
      t9386 = t220*t4573
      t9391 = t80*t234
      t9394 = 4*t5192*t532*t5108-4*t57*t6183*t5108+64*t3959*t832-32*t408
     #0*t9374+24*t9377*t2264-20*t9377*t2267+32*t9304*t9382-t1433*t8155+9
     #6*t9386*t2277-64*t9386*t5138+64*t9391*t5117
      t9405 = t48*t8155
      t9408 = t4*t8155
      t9411 = t2063*t4628
      t9412 = t4*t9411
      t9417 = t423*t4660
      t9423 = t608*t4660
      t9428 = t608*t1093
      t9443 = t5178*t5325
      t9451 = -4*t336*t9423-40*t336*t5849+6*t336*t9428-4*t336*t9321-8*t9
     #304*t8411-192*t4853*t5325*t812+96*t4853*t211*t821-48*t1436*t2625+4
     #8*t9443*t2277+8*t3968*t4874+96*t80*t1211*t812
      t9468 = t68*t8155
      t9496 = t2157*t54
      t9500 = t92*t615
      t9503 = -24*t4606*t6795+8*t235*t2206*t4837-96*t80*t1578*t821+40*t2
     #35*t5849-32*t9443*t2274-6*t6184*t2319+4*t6754*t2293-96*t9391*t2264
     #+128*t9391*t2267+32*t4080*t9496+96*t1436*t1060*t9500
      t9507 = t4819*t211
      t9510 = t423*t596
      t9517 = t101*t1
      t9519 = t677*t9517*t603
      t9523 = t404*t5276
      t9529 = t925*t108
      t9544 = t92*t5276
      t9553 = t4688*Eq
      t9560 = t8247*Mzsq
      t9563 = t1233*t32
      t9568 = t5405*t31
      t9571 = t5405*t32
      t9574 = 64*t939*t2613*t616+32*t191*t2100*t9544-2*t5657*t9311*t32-8
     #*t5192*t532*t905+8*t9553*t8273-4*t5192*Er*t905*Mzsq-4*t9553*t9560-
     #24*t57*t9563-12*t57*t3780+96*t4819*t9568-48*t280*t9571
      t9579 = t9311*t31
      t9605 = t5405*BSnThl
      t9620 = t262*t5804
      t9636 = 48*t57*t3788-256*t206*t9605*t975-128*t280*t5405*t100*t677*
     #t3473-64*t280*t9347*t6852+192*t4994*t5804*t3227+32*t9620*t5572+2*t
     #906*t7898-64*t2067*t4322*t616-4*t2067*t316*t9523+8*t3094*t316*t616
     #+64*t6842*t9523
      t9655 = t677*t9517
      t9685 = t152*t905
      t9688 = t8*SnThr*Eq
      t9693 = t116*t386*t4628
      t9701 = t8*t635*Eq
      t9712 = -4*t7058*t8224-48*t1433*t4445*t102*t751*t1-32*t191*t8241*t
     #616-16*t3770*t9523-32*t9685*t682*t9688-16*t191*t2079*t9693-96*t143
     #6*t711*t616-16*t1436*t2617*t9701+64*t262*t9605*t3158-8*t9553*t532*
     #t5479+4*t6184*t2322
      t9720 = t811*t4660
      t9729 = t1343*t4114
      t9734 = t640*t6514
      t9743 = t2157*t4837
      t9750 = t640*t959
      t9766 = -2*t6215*t2157*t4660-20*t5848*t3035-4*t336*t9743+256*t9729
     #*t2274+96*t9507*t2264-6*t6215*t9750+8*t6215*t2157*t1093-192*t9729*
     #t2277+8*t57*t8481*Mzsq-16*t191*t8563*t1046+80*t336*t1188
      t9771 = t5712*t964
      t9774 = t423*t1093
      t9788 = t2157*t5677
      t9791 = t48*t8451
      t9809 = t4*t2067
      t9820 = t48*t2083
      t9825 = -32*t5037*t5230+24*t5037*t5536+2*t7058*t9405-32*t4606*t723
     #4-24*t336*t9809+32*t4606*t5740+8*t4606*t9405-24*t5037*t2346-8*t704
     #0*t9362+6*t7040*t9820+24*t5037*t5227
      t9827 = t4*t7899
      t9835 = t2074*t3153
      t9841 = t48*t2067
      t9844 = t48*t2064
      t9847 = t4*t8451
      t9850 = t4*t2064
      t9858 = t2157*t38
      t9861 = t640*t1707
      t9881 = 72*t7058*t9374-56*t336*t9858+24*t336*t9861-20*t5105*t9861-
     #6*t5105*t9841+8*t5105*t922*t48+96*t4089*t9844+48*t560*t9809-192*t4
     #089*t9850-128*t6372*t3631+96*t7015*t9850
      t9892 = t683*t916
      t9895 = t2618*t94
      t9899 = t2074*t4629*t2063
      t9916 = t423*t4837
      t9921 = t48*t4859
      t9936 = 16*t7040*t4916+6*t7040*t5479-4*t7058*t5108+8*t7058*t9916-8
     #*t7058*t4865-8*t1994*t9921-8*t7904*t445*t243+256*t3876*t2100*t1078
     #+128*t4769*t4782-4*t5105*t5849-64*t4089*t9408
      t9961 = t423*t816
      t9978 = t2*t4628
      t9985 = -32*t4080*t9961-64*t4089*t9412+32*t4089*t9788+8*t5478*t832
     #+32*t4089*t9405-16*t3094*t1799+32*t2067*t1390-16*t7899*t1390-16*t2
     #064*t9978+4*t906*t2300+2*t5957*t2300
      t10003 = t8358*t615
      t10024 = t4*t1390
      t10034 = t1436*t2063
      t10035 = t4*t1799
      t10038 = t48*t1799
      t10041 = -64*t4080*t832-32*t4080*t4750-32*t4080*t4865+8*t4080*t978
     #8-192*t4994*Er*t1049*t7786+32*t7983*t10024-24*t7058*t48*t943+24*t7
     #058*t9844-16*t7040*t9510+24*t10034*t10035-8*t10034*t10038
      t10047 = t48*t1390
      t10064 = t4*t4859
      t10068 = t208*t7898
      t10083 = t4*t9978
      t10088 = t48*t9978
      t10093 = -12*t10068*t10035+64*t3610*t10047-128*t9386*t2274+2*t7058
     #*t9788-8*t1358*t3063-36*t3768*t10035+24*t1994*t10064-64*t1589*t100
     #83+32*t2893*t10083+32*t1589*t10088+8*t7904*t9961
      t10145 = -2*t7904*t3006*t1635-16*t2835*t9774-16*t7983*t10047+48*t3
     #808*t10035-8*t4589*t9408+32*t4939*t9408-32*t4939*t9405-8*t4589*t94
     #12+64*t4895*t5030+3*t6908*t9750-12*t4962*t5030
      t10190 = -24*t6908*t3037+16*t4962*t4910-3*t5478*t1119+2*t5478*t510
     #6+4*t10068*t10038+4*t3610*t9921+2*t5105*t2322-4*t5105*t9774-128*t3
     #610*t10024+96*t4925*t9820-4*t6908*t640*t3201
      t10207 = t4837*t2063
      t10216 = t4688*t2063
      t10247 = -16*t5061*t5030+4*t6908*t48*t8166-6*t6908*t3094*t48-64*t4
     #909*t2455+48*t4909*t9850-128*t4925*t9362-64*t6372*t4910+64*t6372*t
     #5030+48*t7058*t832-12*t7058*t1119+40*t7058*t9961
      t10266 = t8963*Mzsq
      t10270 = t5405*sqq
      t10291 = t102*t103*t1
      t10297 = 96*t4853*t10270-8*t41*t8156-8*t41*t4940*t2063+16*t41*t594
     #2+t3094*t6453-24*t41*t3019+16*t41*t7983+32*t41*t3021-96*t6663*t616
     #-48*t1433*t4811*t10291+8*t208*t7903*t616
      t10315 = t2063*t4740
      t10321 = t2625*t603
      t10338 = t2064*t4628
      t10354 = -16*t191*t8012*t603-48*t1433*t6254+32*t208*t675*t4628-16*
     #t191*t2625*t4628-3*t5931*t2063-16*t97*t10338+32*t97*t3780-64*t44*t
     #2805+4*t57*t10216*t603-2*t57*t8114*t603-2*t57*t5594*t2063
      t10361 = t7899*t603
      t10397 = t7877*t615
      t10406 = 32*t939*t10321-2*t53*t7714+32*t44*t8163+32*t44*t10207-48*
     #t1436*t779*t9519+64*t208*t2621*t9519+8*t8348*t3718-18*t336*t3007+1
     #6*t3518*t10397-32*t1433*t2617*t9701+24*t5957*t2617*t9688
      t10420 = t27*t4768
      t10433 = t41*t596
      t10445 = t206*t7641
      t10454 = t22*t7898
      t10455 = t2074*t10454
      t10462 = 24*t10433*t6006+64*t2067*t925*t786-16*t3094*t404*t9186+64
     #*t3*t2083*t9523+192*t10445*t3154+18*t235*t9428-12*t235*t9423-12*t5
     #7*t10068+128*t1071*t10455+64*t5061*t1359-32*t4080*t2152
      t10480 = t93*t316
      t10513 = 8*t5105*t423*t905+4*t5478*t1300-3*t5478*t3415-20*t6908*t3
     #765+4*t6908*t3782-64*t3090*t1046-36*t7712*t935+64*t9563*t1327-64*t
     #10068*t473+24*t9057*t935-16*t7722*t2113
      t10522 = t5712*sqq
      t10525 = t5712*Mzsq
      t10549 = t5677*t2063
      t10560 = t262*t7641
      t10565 = -72*t7058*t2171-16*t145*t8296-4*t7058*t2761-40*t7058*t178
     #1-2*t7058*t8182-2*t7058*t10549-6*t7040*t1620+8*t7058*t2781-24*t336
     #*t2000+56*t336*t2209+32*t10560*t6*t1002*t7898
      t10569 = t1436*t2100
      t10572 = t208*t8009
      t10578 = t208*t2100
      t10584 = t6146*t2063
      t10603 = t615*t316
      t10608 = t5957*t2063
      t10619 = 8*t7904*t1355-128*t9563*t473+96*t3808*t473-48*t3808*t1327
     #-12*t906*t683*t10603-48*t10034*t1327-36*t10608*t935+96*t10034*t473
     #-8*t9057*t2113-48*t560*t2003+4*t336*t8163
      t10632 = t2067*BCsThl
      t10637 = t76*t603
      t10640 = t4553*t603
      t10676 = -20*t336*t2805-12*t336*t5571+18*t336*t2003-4*t7040*t3369+
     #2*t10315*t603-3*t9685*t2063+32*t7899*t4322*t616+32*t2067*t116*t412
     #0+48*t7722*t935-32*t5037*t1575+24*t4606*t3578
      t10710 = t32*BCsThl
      t10720 = t3196*t5276
      t10723 = t8451*BCsThl
      t10725 = t3196*t615
      t10738 = 12*t6908*t9771-32*t2068*t10710*t838-2*t2835*t3006*t804+32
     #*t1149*t9693+32*t27*t3123*t10720+32*t27*t10723*t10725-16*t1433*t78
     #98*El*t3891+t7904*t1599-2*t2835*t7520+t2835*t7513+16*t2835*t664
      t10786 = t27*sll
      t10787 = t10786*t1060
      t10793 = 2*t4589*t10549-32*t4589*t2781+2*t4589*t8182+16*t145*t8160
     #+16*t145*t2136*t4628+8*t10569*t10397-64*t27*t2956*t10725-t3094*t46
     #28-16*t2064*t391*t102*t4829-64*t10787*t688*t7743-64*t8732*t5062
      t10805 = t67*t603
      t10836 = t27*t505
      t10841 = -48*t4909*t1672+96*t4925*t1620+6*t5105*t2003-6*t5105*t280
     #5-128*t4089*t2761-64*t4089*t4838+64*t560*t2805-16*t2893*t10088+24*
     #t5848*t3037+32*t10836*t9899-64*t10420*t4782
      t10896 = 3*t5478*t1672-2*t5478*t2761-2*t6908*t7309+2*t6908*t10338+
     #4*t6908*t3784-2*t5105*t5571+t5105*t10207-4*t10578*t7877*t5276+64*t
     #208*t8241*t9-64*t4853*t5983+128*t208*t2613*t616
      t10941 = t5105*t8163-8*t4606*t8182-16*t5061*t3369+32*t4080*t10525+
     #32*t4080*t2781+8*t4080*t10584-8*t4080*t8182+8*t4080*t8456-8*t4080*
     #t10549-128*t10420*t9892-64*t10836*t10455
      t10944 = t280*t8944
      t11007 = 32*t280*t9311*t1046+32*t10944*t102*t4829*t2063+96*t10787*
     #t3337-64*t10944*t102*t118*t7898+96*t3092*t1075+192*t4994*t5405*t16
     #-4*t7712*t1836-2*t10608*t1836+64*t5061*t1737+32*t4080*t1781+4*t300
     #7*t249*t4628
      t11009 = t249*t603
      t11053 = 3*t7722*t1836+8*t7058*t2152+12*t336*t8282+t9057*t1836+t38
     #08*t6453+16*t5037*t2851-16*t5037*t3369+16*t4606*t5858-24*t4606*t28
     #30-32*t41*t2083*t5022-16*t10433*t4941
      t11062 = t4620*t149
      s5 = -384*t7752*t8411-3072*t9729*t5138+24*t8084*t603+96*t7040*t503
     #0+144*t4589*t8456+1152*t1994*t10047-768*t560*t5571+384*t560*t8163+
     #288*t5848*t9771+1536*t9507*t5117-96*t6215*t9734+1536*t7752*t9382+1
     #152*t1436*t6254+72*t9685*t791-192*t8166*t1799+576*t9304*t9356-576*
     #t4606*t5591-192*t4606*t9788+384*t560*t10207-96*t5105*t5863-768*t30
     #07*t10805+384*t8282*t10805+1536*t10632*t10640-288*t7058*t1672-96*t
     #7040*t2851-96*t336*t10207+576*t336*t664-48*t906*t4759+72*t5931*t79
     #1+384*t1433*t8563-768*t208*t8012-768*t4606*t2810+576*t5037*t1677+4
     #8*t3768*t6453-72*t7704*t1836-192*t4606*t8456+192*t4606*t10549-96*t
     #9291*t11009+192*t9295*t11009-288*t10608*t2113+192*t5478*t2171+192*
     #t5478*t8884-48*t5105*t8282+3072*t10578*t9500
      s4 = -48*t7058*t10584-384*t7040*t1958+576*t5037*t3015-576*t5037*t3
     #422-576*t336*t5863-384*t5061*t2851-768*t4080*t2171-432*t336*t9841+
     #2304*t7015*t9844+1536*t7015*t9847+4608*t1071*t9835-384*t336*t9827+
     #144*t7040*t821-192*t4606*t10584-48*t5105*t4952-72*t9685*t2100-72*t
     #5931*t2100+384*t97*t10361+48*t906*t8009+96*t6682*t1060-2304*t10420
     #*t9895+1152*t9377*t832+576*t9308*t2277-3072*t4089*t9847-2304*t1083
     #6*t9835-768*t2219*t9693-384*t9304*t9720-72*t6908*t3820+144*t6908*t
     #3780-48*t6908*t10361+768*t4080*t10522+1536*t5061*t1774-576*t906*t6
     #099+384*t1433*t6251-384*t8242*t5042+384*t1436*t6325+192*t41*t8163+
     #192*t41*t10207+1536*t4895*t4910-384*t41*t2805+1536*t10572*t9529-28
     #8*t3768*t10038+288*t3610*t10064-144*t57*t10338+s5
      s5 = -192*t8570*t786-3072*t2213*t5042+1536*t2213*t8372-144*t57*t10
     #361+1536*t2213*t5046+768*t1433*t6325-768*t7015*t9412-768*t7015*t94
     #08+768*t7015*t9405-96*t7058*t5106-576*t5957*t2079-1344*t235*t2293+
     #1920*t235*t1188-192*t2067*t4859+1536*t1071*t10003+1536*t5061*t8375
     #+384*t4606*t8229-288*t9563*t10038+1536*t3019*t10024+1152*t3021*t10
     #047-2304*t3021*t10024+864*t9563*t10035-192*t2835*t9417-768*t10068*
     #t1327-288*t7712*t2113-48*t7058*t8456-576*t7058*t1300+288*t7058*t34
     #15-144*t7040*t2190+1344*t336*t2293+576*t336*t9417+192*t7058*t4750+
     #96*t10572*t10397-768*t6560*t9523-480*t5105*t2000+96*t5105*t2209-38
     #4*t7040*t1737+384*t7040*t1774-96*t7058*t4838+192*t7058*t10522-192*
     #t7058*t10525-1536*t5061*t1958-768*t5037*t7078-384*t4606*t7979
      s6 = -2304*t10569*t9529-24*t10793-24*t10619-24*t10676-24*t10738-24
     #*t10565-24*t10513-24*t10462-24*t10406-24*t10354-24*t10297+144*t57*
     #t10216*sqq-96*t57*t8114*sqq+384*t9620*t93*t10207-192*t208*t2075*t9
     #523-768*t939*t2621*t9519-192*t208*t8563*t603+1536*t4819*t9311*sqq-
     #72*t5657*t5405*t905-768*t4819*t5912*t2063+96*t8166*t316*t616-3072*
     #t5957*t711*t9+768*t191*t2613*t9523
      s3 = s6-1536*t3073*t4276*t7278+3072*t4994*t7641*t4635-144*t336*t43
     #64*t2063+48*t6215*t8223*t816-2304*t5178*t9347*t3246+1536*t939*t210
     #0*t9500+2304*t1436*t2613*t9-768*t191*t8009*t9500+3072*t5957*t1060*
     #t9529-4608*t206*t9330*t3158+96*t6908*t423*t959+384*t145*t8223*t663
     #+384*t145*t2157*t5712-384*t2068*t75*t6452-768*t8732*t68*t1635-768*
     #t2213*t68*t5677+48*t2835*t3006*t5677-24*t53*t2067*t4628+2304*t1433
     #*t2613*t9-576*t906*t2089*t3891-576*t5957*t682*t9701+s5+s4
      s6 = 1152*t1436*t2621*t9655-1536*t5957*t779*t9655+1152*t3073*t5405
     #*t917+1536*t208*t1060*t9544+1152*t3073*t5804*t3210-96*t6908*t423*t
     #6514-384*t145*t8223*t243+288*t336*t48*t7899+4608*t3876*t1060*t3337
     #+384*t4606*t8223*t3014-576*t4606*t2157*t6183-96*t2835*t7199*t1635-
     #3072*t3876*t1041*t4839+3072*t10445*t6*t10454+6144*t4994*t4276*t917
     #+3072*t10786*t2100*t1078-480*t6908*t2157*t1707-192*t208*t2079*t462
     #8+576*t4606*t2206*t152-576*t41*t964*t5022+384*t3415*t67*t4628-768*
     #t2956*t4553*t4628
      s5 = s6-384*t9685*t683*t10480+288*t5957*t2100*t3718+288*t906*t2618
     #*t10480-576*t3094*BCsThl*t10637-1536*t206*t5912*t2093-1536*t206*t9
     #311*t616+768*t7899*t925*t786+288*t8166*BCsThl*t10637+288*t10632*t7
     #6*t4628-3072*t3*t2956*t10725+1536*t3*t3123*t10720+1536*t3*t10723*t
     #10725+288*t5957*t683*t10603-192*t7338*t15*t6453-768*t12*t2067*t616
     #-768*t7899*BCsThl*t10640-1536*t3*t8451*t616+96*t6908*t5712*t32+384
     #*t12*t7899*t616-768*t97*t8451*t1046+3072*t3*t2064*t616-24*t11053
      s4 = -24*t11007-24*t10941-24*t10896-24*t10841-24*t10247-24*t10190-
     #24*t10145-24*t10093-24*t10041-24*t9985-576*t906*t5445*t6*t614*Eq-2
     #4*t9936+1536*t5178*t8338*t677*t678*t2063+1152*t262*t9330*t8*t635*t
     #32-24*t9766+768*t7752*t9720-288*t4962*t3369+144*t4589*t10584+96*t4
     #962*t2851+192*t2835*t5863+192*t7904*t1781-768*t2068*t1670-2304*t40
     #89*t1672+24*t1433*t9411+384*t7040*t8375+1152*t1433*t2625+768*t9685
     #*t2609-1536*t5957*t675-192*t4080*t9405-768*t4080*t9916+96*t5105*t2
     #293-480*t6908*t3035+4608*t4769*t9895+1536*t1071*t9899-96*t5105*t98
     #58-192*t5478*t9374+2304*t6372*t9820-1536*t7015*t9791-6144*t4769*t9
     #892-768*t7058*t2455-2304*t6372*t5159+576*t7058*t9850+192*t7058*t94
     #96+96*t336*t7996+s5
      s6 = 384*t262*t9311*BSnThl*t975-96*t7058*t7918-192*t7040*t1221+144
     #*t7040*t2149-768*t3019*t10047-2304*t1994*t10024-96*t5105*t664-768*
     #t10836*t10003+384*t5696*t9523+96*t7058*t9791-576*t336*t9774+192*t5
     #478*t1781-480*t336*t2322-96*t7058*t5117+384*t5037*t5030-288*t336*t
     #4952+432*t336*t2319-576*t906*t2079-384*t4606*t5678-768*t208*t7898*
     #t98*t10291+1152*t10560*t6*t5180*t2063+768*t2064*t19*t4628*t9
      s5 = s6+384*t7899*t1*t100*t679-192*t2067*t192*t6*t4629-2304*t27*t1
     #52*t1041*t5684-2304*t10944*t102*t5197*t2063-768*t12*t2083*BSnThl*t
     #948+384*t57*t4688*t37*t31+3072*t206*t5912*El*t5722-144*t57*t4688*t
     #47*t5479+384*t4606*t7287+3072*t4853*t10266+576*t41*t3610-144*t5192
     #*t9571+1152*t5178*t9568-1536*t5178*t9560+1536*t4925*t5030-96*t145*
     #t8444+768*t7015*t9788+24*t5105*t9468-3072*t4925*t4910-1536*t5061*t
     #4916-288*t57*t10034+768*t560*t9827+1536*t4089*t9791
      s2 = 4608*t4925*t5159-6144*t4925*t3631+24*t5105*t9743-192*t5478*t9
     #961+384*t3808*t10038+48*t6908*t9734-768*t4589*t4750-768*t4939*t941
     #2+96*t5105*t9417-144*t4589*t9405+768*t4939*t9788-144*t4589*t9788+1
     #92*t7058*t7127+48*t6908*t9291+1536*t41*t2211-3072*t6372*t9362-432*
     #t57*t3808+48*t5192*t9579-768*t280*t9579-72*t9553*t9568-96*t6682*Mz
     #sq-24*t9825-1152*t7752*t9356-24*t8084*t1836-2304*t1712*t9500-576*t
     #5957*t6099-2304*t3518*t9529-1536*t208*t10321+48*t10315*t1836-24*t5
     #3*t9059+1152*t7058*t8884+3072*t9507*t2267+1536*t5061*t9510+3072*t1
     #325*t5800+96*t336*t9468-24*t9881+s5-24*t9712-24*t9636-24*t9574+s4-
     #24*t9503-24*t9451-24*t9394+s3-24*t9328
      s3 = t569*t151
      s1 = s2*s3
      s2 = t11062*t4108*t567*t2228*t572
      A2e0c6p1t1n = s1*s2
      A2e0c1p1t2n = A2e0c1p1t1n
      A2e0c2p1t2n = A2e0c2p1t1n
      A2e0c3p1t2n = A2e0c3p1t1n
      A2e0c4p1t2n = A2e0c4p1t1n
      A2e0c5p1t2n = A2e0c5p1t1n
      A2e0c6p1t2n = A2e0c6p1t1n
      t11070 = (-sqq+sllmin)*t567
      t11073 = log(DCMPLX(-t11070/dels))
      t11074 = t11073*Eq
      t11075 = t11074*Er
      t11077 = -j2+4*t11075
      t11078 = t11077*BtQlq
      t11079 = BCsThlq**2
      t11080 = t2*t11079
      t11081 = t11080*t46
      t11082 = 8*t11081
      t11083 = BCsThlq*CsThr
      t11084 = t11083*t865
      t11085 = 4*t11084
      t11086 = t2*BCsThlq
      t11087 = t11086*CsThr
      t11088 = BSnThlq*SnThr
      t11089 = t11088*t7
      t11091 = 16*t11087*t11089
      t11092 = 16*t11087
      t11093 = 8*t2
      t11094 = 4*t865
      t11095 = t22*t1
      t11096 = t11088*t11095
      t11097 = 4*t11096
      t11098 = 2*t54
      t11099 = BSnThlq**2
      t11100 = t2*t11099
      t11102 = 8*t11100*t386
      t11103 = t2*BSnThlq
      t11105 = 16*t11103*t93
      t11106 = t11082+t11085-t11091-t11092+t31+t11093-t11094-t11097+t110
     #98+t11102+t11105
      t11110 = 1/BtQl
      t11111 = 1/j2
      t11112 = t11110*t11111
      t11113 = t11112*t149
      A1e0c1p1t1n1 = -24*t11078*t11106*t151*t156*t567*t11113
      t11117 = t11083*t334
      t11118 = 16*t11117
      t11119 = t11086*t297
      t11120 = 16*t11119
      t11121 = t11088*t22
      t11122 = t11087*t11121
      t11123 = 16*t11122
      t11124 = t37*BSnThlq
      t11125 = t11124*t93
      t11127 = 64*t11087*t11125
      t11128 = t11079*BCsThlq
      t11129 = t2*t11128
      t11130 = t11129*CsThr
      t11132 = 16*t11130*t11121
      t11133 = t11128*CsThr
      t11134 = t11133*t33
      t11135 = 4*t11134
      t11137 = 16*t11129*t297
      t11138 = t11079**2
      t11141 = 8*t2*t11138*t68
      t11142 = t11083*t33
      t11143 = 4*t11142
      t11144 = t11080*t68
      t11145 = 8*t11144
      t11146 = t11079*t46
      t11148 = 32*t201*t11146
      t11149 = t11086*t225
      t11150 = 64*t11149
      t11151 = t11118+t11120-t32+t11123-t11127-t11132+t11135-t11137+t111
     #41-t11143-t11145+t11148-t11150+t34-t36+t39
      t11152 = t11103*t94
      t11153 = 16*t11152
      t11154 = t1002*t1
      t11155 = t11088*t11154
      t11156 = 4*t11155
      t11157 = t11099*t101
      t11158 = t11157*t103
      t11160 = 32*t201*t11158
      t11161 = t11100*t1074
      t11162 = 8*t11161
      t11164 = 64*t201*t11089
      t11165 = t32*t11079
      t11167 = t11124*SnThr*t11095
      t11168 = 16*t11167
      t11169 = t11100*t101
      t11172 = 8*t11169*t118*t11079
      t11176 = 16*t11103*SnThr*t22*t11079
      t11177 = t33*t11079
      t11178 = t11089*t11177
      t11179 = 4*t11178
      t11180 = 4*t11177
      t11181 = t35*t11079
      t11182 = 8*t11181
      t11183 = t11079*t37
      t11185 = 2*t11183*t31
      t11186 = t202+t205-t11153+t11156+t11160-t11162+t11164+t11165-t1116
     #8+t11172+t11176-t11179-t11180+t11182+t11185-t335
      t11187 = t11151+t11186
      t11190 = -1+t11079
      t11191 = t11190**2
      A1e0c2p1t1n1 = -24*t11078*t11187*t363/t11191*t151*t11113
      t11201 = t1*t11079
      t11202 = t11201*t445
      t11204 = t38*t11083
      t11206 = t370*t11099
      t11209 = t11099*BSnThlq
      t11210 = t370*t11209
      t11213 = t370*BSnThlq
      t11216 = 64*t11122
      t11217 = t118*t11083
      t11222 = t1*BSnThlq
      t11224 = t1002*t11083
      t11225 = t11222*SnThr*t11224
      t11228 = 8*t11117-56*t11119-8*t11129*t430-4*t11202-2*t11204-96*t11
     #206*t386-32*t11210*t698-96*t11213*t93-t371+t32-t11216-8*t11169*t11
     #217+16*t11081*t11121+4*t11225+12*t11142
      t11231 = 8*t11155
      t11234 = t370*BCsThlq
      t11235 = t11234*CsThr
      t11240 = t370*t11079
      t11241 = t11240*t46
      t11249 = t544*BCsThlq
      t11250 = 40*t11144-t435+t436+t39+48*t11152-t11231+24*t11161-8*t111
     #67+96*t11235*t11158+192*t11235*t11089-96*t11241*t11089-96*t11241+3
     #2*t370*t11128*t417+96*t11235-t11249-t540
      t11255 = 1/t11190
      t11257 = t151*t11255*t156
      A1e0c3p1t1n1 = -48*t11078*(t11228+t11250)*t363*t149*t11112*t11257
      t11261 = t11078*sqq
      t11267 = t11112*t151*t7894*t156
      A1e0c4p1t1n1 = 24*t11261*Mzsq*t11106*t149*t11267
      t11270 = Mzsq*t11077*BtQlq
      t11272 = -1+BCsThlq
      t11273 = t11272**2
      t11274 = 1/t11273
      t11277 = 1+BCsThlq
      t11278 = t11277**2
      t11279 = 1/t11278
      t11280 = t11112*t11279
      t11281 = t7895*t11280
      A1e0c5p1t1n1 = 24*t11270*t11187*t151*t11274*t11281
      t11283 = Mzsq*BtQlq
      t11284 = t11073*t43
      t11286 = 96*t11284*t317
      t11287 = t11073*t19
      t11289 = 32*t11287*t8963
      t11291 = 4*t11074*t8247
      t11292 = j2*t37
      t11294 = 8*t11292*t865
      t11295 = j2*t370
      t11297 = 96*t11295*t11083
      t11298 = t11128*t417
      t11300 = 32*t11295*t11298
      t11302 = 96*t11295*t11146
      t11303 = j2*t32
      t11304 = t11303*t11083
      t11307 = t1788*BCsThlq
      t11310 = 16*t11287*Er*BSnThlq*t93*t11307
      t11311 = t11073*t7928
      t11312 = Er*BCsThlq
      t11313 = t11311*t11312
      t11314 = CsThr*BSnThlq
      t11315 = t11314*t93
      t11317 = 768*t11313*t11315
      t11318 = t11311*Er
      t11319 = 128*t11318
      t11320 = j2*t2
      t11321 = t11320*t11083
      t11323 = 64*t11321*t11121
      t11324 = t11320*t11146
      t11326 = 16*t11324*t11121
      t11327 = 32*t11295
      t11329 = 96*t11295*t11089
      t11330 = t11286-t11289+t11291+t11294-t11297-t11300+t11302+t11304+t
     #11310+t11317-t11319+t11323-t11326+t11327+t11329
      t11334 = 8*t11292*BSnThlq*t93*t865
      t11335 = t11295*BCsThlq
      t11337 = 192*t11335*t11315
      t11339 = t46*BSnThlq
      t11340 = t11339*t93
      t11342 = 96*t11295*t11079*t11340
      t11343 = CsThr*t11099
      t11344 = t11343*t386
      t11346 = 96*t11335*t11344
      t11347 = t11284*Er
      t11348 = t11157*t118
      t11350 = 96*t11347*t11348
      t11352 = 192*t11347*t11121
      t11353 = t11287*Er
      t11354 = t11088*t1002
      t11356 = 32*t11353*t11354
      t11358 = 8*t11075*t11204
      t11361 = 32*t11353*t11083*t54
      t11364 = 8*t11320*t11157*t11217
      t11365 = j2*t1
      t11368 = 4*t11365*t11088*t11224
      t11370 = 96*t11295*t11158
      t11372 = t11209*t518*t520
      t11374 = 32*t11295*t11372
      t11375 = t11298*sqq
      t11377 = 8*t11320*t11375
      t11378 = t11284*t11312
      t11379 = t11314*t94
      t11381 = 256*t11378*t11379
      t11383 = 384*t11313*t11344
      t11384 = t11334-t11337+t11342-t11346+t11350+t11352-t11356-t11358+t
     #11361+t11364-t11368+t11370+t11374+t11377-t11381+t11383
      t11386 = Er*t11079
      t11390 = 64*t11284*t11386*t11339*t94
      t11393 = t297*BCsThlq
      t11396 = 32*t11284*Er*t11099*t386*t11393
      t11399 = 32*t11287*t532*t11121
      t11402 = 384*t11311*t11386*t11340
      t11404 = 32*t11287*t533
      t11406 = 8*t11074*t8270
      t11407 = t11386*t46
      t11409 = 384*t11311*t11407
      t11410 = t11146*t31
      t11412 = 4*t11365*t11410
      t11414 = 2*t11292*t11307
      t11415 = j2*BCsThlq
      t11416 = t1788*t1
      t11418 = 12*t11415*t11416
      t11419 = t11146*sqq
      t11421 = 40*t11320*t11419
      t11423 = 56*t11320*t11393
      t11427 = 128*t11311*Er*t11128*t417
      t11428 = t11312*CsThr
      t11430 = 384*t11311*t11428
      t11431 = t11390-t11396-t11399-t11402-t11404+t11406-t11409+t11412+t
     #11414-t11418-t11421+t11423+t11427+t11430-t11303
      t11434 = 8*t11415*CsThr*t334
      t11437 = 24*t11320*t11099*t1074
      t11439 = j2*BSnThlq*SnThr
      t11441 = 8*t11439*t11154
      t11444 = 48*t11320*BSnThlq*t94
      t11446 = 160*t11347*t11419
      t11448 = 224*t11347*t11393
      t11450 = 128*t11318*t11372
      t11452 = 384*t11318*t11158
      t11454 = 384*t11318*t11089
      t11456 = 4*t11075*t11249
      t11458 = 32*t11347*t11375
      t11460 = 16*t11353*t11410
      t11462 = 48*t11353*t11307
      t11463 = j2*t31
      t11465 = 8*t11463*t1
      t11467 = 2*t11292*t31
      t11469 = 24*t11320*sqq
      t11470 = -t11434-t11437+t11441-t11444+t11446-t11448-t11450-t11452-
     #t11454-t11456-t11458-t11460+t11462+t11465-t11467-t11469
      t11479 = t7894*t151*t11111*t11255*t156
      A1e0c6p1t1n1 = 48*t11283*(t11330+t11384+t11431+t11470)*t11110*t149
     #*t11479
      t11481 = t11085-t11094+t11105+t31+t11093+t11102-t11098-t11091+t110
     #82-t11092-t11097
      A1e0c7p1t1n1 = 24*t11261*Mzsq*t11481*t149*t11267
      t11487 = t11118+t11120-t32+t11123-t11127-t11132+t11135-t11137+t111
     #41-t11143-t11145+t11148-t11150+t34-t36+6*t38
      t11488 = t202-t205-t11153+t11156+t11160-t11162+t11164+t11165-t1116
     #8+t11172+t11176-t11179-t11180+t11182-t11185-t335
      A1e0c8p1t1n1 = 24*t11270*(t11487+t11488)*t151*t11274*t11281
      t11494 = t11286-t11289+t11291-t11294-t11297-t11300+t11302+t11304+t
     #11310+t11317-t11319+t11323-t11326+t11327+t11329
      t11495 = -t11334-t11337+t11342-t11346+t11350+t11352-t11356+t11358-
     #t11361+t11364-t11368+t11370+t11374+t11377-t11381+t11383
      t11497 = t11390-t11396+t11399-t11402+t11404-t11406-t11409+t11412-t
     #11414-t11418-t11421+t11423+t11427+t11430-t11303
      t11498 = t11434-t11437+t11441-t11444+t11446-t11448-t11450-t11452-t
     #11454-t11456-t11458-t11460+t11462+t11465+t11467-t11469
      A1e0c9p1t1n1 = 48*t11283*(t11494+t11495+t11497+t11498)*t11110*t149
     #*t11479
      t11505 = 8*t675
      t11506 = 8*t1390
      t11507 = t11086*t811
      t11510 = t4516*t11079*Mzsq
      t11512 = t11079*t603
      t11513 = t4516*t11512
      t11515 = t11086*t838
      t11518 = 2*t32
      t11519 = Mzsq*BSnThlq
      t11521 = t11087*t11519*t93
      t11523 = t11088*t614
      t11524 = t11087*t11523
      t11527 = 8*t11142
      t11529 = t11505+t11506-16*t11507+8*t11510+8*t11513-16*t11515+32*t1
     #1119+t816-t11518-16*t11521-16*t11524+32*t11122-t11527-16*t11144+t9
     #64+t435
      t11530 = 16*t35
      t11531 = sqq*BCsThlq
      t11532 = t11531*t10637
      t11533 = 4*t11532
      t11534 = t11100*t1045
      t11536 = t675*t11089
      t11538 = t675*t11158
      t11540 = t11103*t615
      t11542 = t2609*t11393
      t11543 = 4*t11542
      t11546 = 2*t1737
      t11547 = t2609*BSnThlq
      t11548 = t11547*t94
      t11549 = 4*t11548
      t11551 = sqq*BSnThlq*t7
      t11553 = t11551*t635*t1
      t11554 = 4*t11553
      t11555 = 4*t2970
      t11556 = 4*t1800
      t11557 = 2*t1958
      t11558 = -t11530-t200+t11533+8*t11534+16*t11536+8*t11538+16*t11540
     #+t11543-32*t11152+t11231-16*t11161+t11546-t11549-t11554-t11555-t11
     #556+t11557
      t11561 = BtQlq*t151
      t11562 = t11561*t2052
      t11564 = t11112*t156
      A1e0c10p1t1n1 = 24*(t11529+t11558)*t11077*t11562*t595*t11564
      t11571 = 2*t905
      t11574 = t31*t11079
      t11577 = sqq*t11079
      t11579 = 8*t1390*t11577
      t11580 = t4516*t11574
      t11582 = t11086*t1788
      t11587 = 8*t675*t11577
      t11588 = t3842-64*t11149*t11523+128*t11149*t11121+t11571-t3925+2*t
     #38*t11512-4*t1799*t11574+t11579+16*t11580-32*t11582-4*t964*t11201+
     #t11587
      t11589 = t1*t11128
      t11600 = t11086*t11314*t657
      t11602 = 4*t1707
      t11607 = t11086*t1859
      t11609 = BCsThlq*t1
      t11610 = t11609*t1789
      t11612 = -8*t11589*t544-16*t4516*t11138*t31+32*t11129*t1788+2*t831
     #*t11574+16*t11600-t943-t11602-t4228-t4369+8*t675*t11138*t46*sqq+16
     #*t11607-4*t11610
      t11614 = t11512*sqq
      t11615 = t4516*t11614
      t11617 = Mzsq*BCsThlq
      t11618 = t11617*t11416
      t11620 = t675*t11393
      t11622 = t675*t11419
      t11631 = t11100*t5799
      t11637 = -8*t11615-4*t11618+16*t11620-8*t11622-64*t11617*t4553*t37
     #-t4395-2*t905*t11079+16*t5394-8*t3041+16*t11631+128*t11086*t294-64
     #*t201*t11419
      t11654 = t11133*sqq
      t11657 = t11083*t1364
      t11662 = t11507*t11121
      t11664 = t11222*t7277
      t11666 = t11103*t916
      t11668 = 32*t11080*t423*t603-64*t11086*t4531+8*t4516*t11138*t603*s
     #qq+4*t11589*t1789-16*t11129*t1859+4*Mzsq*t11128*t11416-16*t675*t11
     #654-32*t11657-16*t675*t11133*t11121+16*t11662-8*t11664+32*t11666
      t11671 = t11103*t7
      t11675 = t11222*t7
      t11679 = t11100*t103
      t11683 = t201*t11099
      t11688 = t201*BSnThlq
      t11693 = t11124*t7
      t11695 = t11693*t737*t1
      t11697 = t11671*t974
      t11699 = t675*t11099
      t11700 = t11699*t1074
      t11702 = t11679*t3473
      t11704 = t11675*t3157
      t11706 = -32*t11671*t737*t11079+8*t11675*t3031*t11079-16*t11679*t9
     #34*t11079+32*t11683*t1045-64*t11683*t1074+64*t11688*t615-128*t1168
     #8*t94+32*t11695-16*t11697-8*t11700-8*t11702+4*t11704
      t11707 = t11519*SnThr
      t11708 = t11707*t11154
      t11711 = t11617*CsThr*t334
      t11714 = t11083*t37*t1800
      t11719 = t675*BSnThlq
      t11720 = t11719*t94
      t11726 = t11222*t93
      t11727 = t1958*t11726
      t11730 = t11693*t635*t865
      t11732 = t11087*t11354
      t11737 = 4*t11708+16*t11711+16*t11714-t1093+32*t675*t11183*t46-16*
     #t11720+64*t4188*t11089+32*t4188*t11158-16*t11727-16*t11730-32*t117
     #32+8*t11679*t678*t11577
      t11747 = t11609*t544
      t11751 = 16*t5394*t11079
      t11757 = -4*t11675*t635*t11574+32*t11130*t11354-4*t11707*t1002*t11
     #201+8*t11747+t4087+t943*t11079-t11751+8*t3041*t11079-4*t1707*t1107
     #9+t1093*t11079+t4596-t4597
      t11764 = t93*t11577
      t11772 = t4598+t4599+t4602+t4605+16*t11671*t635*t11577+8*t11699*t3
     #86*t11577+16*t11719*t11764-t4095+t4097-64*t11507*t11125-16*t11129*
     #t11314*t657+t4099
      t11780 = t4619*t149
      A1e0c11p1t1n1 = 24*t11078*(t11588+t11612+t11637+t11668+t11706+t117
     #37+t11757+t11772)*t11274*t2052*t594*t11780*t11280
      t11783 = 24*t3924
      t11792 = t370*t417
      t11793 = t11128*sqq
      t11800 = t3842-t11571+t11783-t11083*t1093+8*t71*t11165+4*t11083*t1
     #707-t11617*t544+192*t11240*t68-64*t11792*t11793-192*t11234*t297+32
     #*t1642*t11298
      t11810 = t11128*t603
      t11813 = t2*t417
      t11820 = t11609*t811
      t11823 = 96*t1642*t11083-96*t1642*t11146-96*t370*t46*t11512+96*t11
     #234*t838+32*t11792*t11810+16*t11813*t11128*t31-80*t11580+112*t1158
     #2-64*t11600+t943+4*t11820*t11354-t11602
      t11831 = t1788*t37
      t11834 = t816*t11079
      t11850 = -4*t2609*t11410-2*t11083*t1781-8*t675*t11375-2*t11617*t11
     #831-4*t71*t11834-8*t11813*t11810*sqq-96*t1642*t11089-96*t1642*t111
     #58-32*t1642*t11372+192*t11206*t1074+192*t11213*t94
      t11861 = 8*t4368
      t11868 = 8*t1166
      t11869 = 64*t11210*t698*sqq-96*t11213*t615-96*t11206*t1045-32*t112
     #10*t702*t518-t11861-56*t11607+12*t11610+40*t11615+12*t11618-56*t11
     #620+40*t11622-t11868
      t11872 = 48*t5394
      t11873 = 16*t3041
      t11875 = 64*t6533
      t11876 = 32*t1642
      t11877 = 32*t1592
      t11883 = -t11872+t11873-48*t11631+t11875-t11876-t11877-16*t11657-6
     #4*t11662+16*t11664-96*t11666+16*t11695
      t11896 = t11099*t103
      t11900 = 48*t11697+24*t11700+24*t11702-8*t11704-8*t11708+8*t11711+
     #8*t11714+t1093+48*t11720+192*t11235*t11523-192*t11235*t11348+96*t1
     #1235*t11896*t678
      t11903 = t11079*BSnThlq
      t11907 = t1642*BCsThlq
      t11916 = t11609*CsThr
      t11931 = -96*t1642*t46*t11903*t93+192*t11907*t11315+96*t11907*t113
     #44-96*t11241*t11523+192*t11241*t11121-8*t11916*BSnThlq*t7*t3031-32
     #*t11081*t11354+16*t11087*t11896*t934-384*t11235*t11121-8*t11727-8*
     #t11730+128*t11732
      t11939 = t11609*t11314
      t11948 = 16*t1364
      t11949 = 24*t4041
      t11950 = 8*t3548
      t11951 = 8*t1754
      t11952 = -24*t11747-8*t11086*t11343*t6377+t4087+16*t4516*t11903*t6
     #57+4*t11939*t951-8*t11507*t11348+16*t11510*t11121+2*t11083*t905+t1
     #1948+t11949-t11950-t11951
      A1e0c12p1t1n1 = 48*t11078*(t11800+t11823+t11850+t11869+t11883+t119
     #00+t11931+t11952)*t4620*t567*t11113*t11257
      A1e0c1p1t2n1 = A1e0c1p1t1n1
      A1e0c2p1t2n1 = A1e0c2p1t1n1
      A1e0c3p1t2n1 = A1e0c3p1t1n1
      A1e0c4p1t2n1 = A1e0c4p1t1n1
      A1e0c5p1t2n1 = A1e0c5p1t1n1
      A1e0c6p1t2n1 = A1e0c6p1t1n1
      A1e0c7p1t2n1 = A1e0c7p1t1n1
      A1e0c8p1t2n1 = A1e0c8p1t1n1
      A1e0c9p1t2n1 = A1e0c9p1t1n1
      A1e0c10p1t2n1 = A1e0c10p1t1n1
      A1e0c11p1t2n1 = A1e0c11p1t1n1
      A1e0c12p1t2n1 = A1e0c12p1t1n1
      t11961 = 4*t1
      A2e0c1p1t1n1 = 96*t11270*t31*(-t11961-4*t11726+sqq+4*t11916)*t149*
     #t11267
      t11969 = t596*BtQlq
      t11974 = t865*t11079
      t11980 = t7868+t11094-16*t343*t11089+4*t11133*t865-t7882+t11574-4*
     #t11974+16*t11609*t225+t11097-t11085-4*t11089*t11974-t31
      t11986 = t149*t11111
      A2e0c2p1t1n1 = 96*t11969*t11077*t11980*t151*t11274*t598*t594*t1198
     #6*t11110*t11279
      t12013 = t11415*t1788+32*t11321-16*t11324-32*t11287*t317+4*t11074*
     #t8963-16*t11320-t11463+8*t11365*sqq+128*t11347*t11089+16*t11287*t1
     #1312*t11379-4*t11365*t11083*t11121-128*t11378*t11315+64*t11284*t11
     #407
      t12041 = -16*t11320*t11158-32*t11320*t11089+4*t11365*t11419-12*t11
     #415*t297*t1-128*t11284*t11428-16*t11353*t11419+48*t11353*t11393-4*
     #t11075*t11307+64*t11347*t11158+8*t11439*t11095-32*t11353*t11121+32
     #*t11320*BCsThlq*t11315+64*t11347
      A2e0c3p1t1n1 = 192*t11969*(t12013+t12041)*t11110*t149*t11479
      t12048 = 8*t865
      t12049 = 4*t2609
      t12050 = 4*t1799
      t12054 = 2*t31
      t12060 = 4*t11820+t656+t596+t12048-t12049-t12050+8*t11096-4*t11222
     #*t615-t12054+4*t11609*t838-8*t11084-4*t2609*t11089
      A2e0c4p1t1n1 = 24*sqq*t12060*t11078*t151*t11062*t11564
      t12082 = -32*t11117-4*t2609*t11577-4*t656*t11201-t816+t11518-8*t11
     #134+t11527+4*t11793*t10637+4*t2609*t11654+16*t11609*t4531-t964+16*
     #t11617*t76*t37-t435-t9238-t11533-t11543-t11231-2*t11165
      t12097 = -16*t9212*t11089-16*t343*BSnThlq*t615+t9248-4*t11547*t117
     #64-4*t11551*t635*t11201+t11549+t11554+32*t11167+8*t11178+t11834+t9
     #64*t11079+t11555+t11556+8*t11177+t9249-t9253-t9254+t9260
      A2e0c5p1t1n1 = 24*t11078*(t12082+t12097)*t11274*t2052*t594*t151*t1
     #49*t11280
      t12106 = 16*t675
      t12107 = 16*t1390
      t12109 = t31*BCsThlq
      t12124 = t12106+t12107-t11617*t1788-t12109*t838-32*t11507+16*t1151
     #0+16*t11513-32*t11515+64*t11119+8*t11202+t816-t11518-32*t11521-32*
     #t11524+t11216-8*t11225-24*t11142-32*t11144-4*t71*t11614
      t12129 = 16*t33
      t12130 = 32*t35
      t12144 = 8*t2970
      t12145 = 8*t1800
      t12147 = -4*t2609*t11419+t964+4*t11820*t11121+t12129-t12130+12*t11
     #532+16*t11534+32*t11536+16*t11538+32*t11540+12*t11542-64*t11152+16
     #*t11155-32*t11161-8*t11548-8*t11553+4*t11939*t657-t12144-t12145+2*
     #t11249
      A2e0c6p1t1n1 = 48*t11078*(t12124+t12147)*t151*t11110*t11986*t594*t
     #2052*t11255*t156
      A2e0c1p1t2n1 = A2e0c1p1t1n1
      A2e0c2p1t2n1 = A2e0c2p1t1n1
      A2e0c3p1t2n1 = A2e0c3p1t1n1
      A2e0c4p1t2n1 = A2e0c4p1t1n1
      A2e0c5p1t2n1 = A2e0c5p1t1n1
      A2e0c6p1t2n1 = A2e0c6p1t1n1
      t12157 = t11609*sqq
      t12158 = 4*t12157
      t12159 = 16*t11086
      t12160 = 8*t11080
      t12161 = -t11094+t11098+t12158+t11093+t31-t12159+t12160
      t12162 = j2*j4
      t12163 = log(DCMPLX(-t11070))
      t12164 = t12163**2
      t12165 = t53*t12164
      t12166 = 2*t12165
      t12167 = t12164*CsThr
      t12169 = 2*t53*t12167
      t12170 = log(DCMPLX(dels))
      t12171 = t12170**2
      t12172 = t53*t12171
      t12173 = 2*t12172
      t12174 = t12171*CsThr
      t12176 = 2*t53*t12174
      t12177 = log(DCMPLX(delc))
      t12178 = t11073*t12177
      t12180 = 4*t53*t12178
      t12181 = t12178*CsThr
      t12183 = 4*t53*t12181
      t12184 = Er*j4
      t12186 = 4*t11074*t12184
      t12191 = log(DCMPLX(Er/delc/Eq))
      t12192 = t12191*j2
      t12193 = t12192*CsThr
      t12194 = t12162+t12166-t12169-t12173+t12176-t12180+t12183-t12186-t
     #12192+t12193
      t12196 = t11561*t567
      t12199 = 1/(-1+CsThr)
      t12200 = t149*t12199
      t12201 = 1/j4
      t12202 = t11111*t12201
      t12203 = t12202*t11110
      A1e0c1p1t1n2 = -12*t12161*t12194*t12196*t12200*t12203
      A1e0c2p1t1n2 = 0
      t12207 = t12194*BtQlq
      t12208 = 64*t11234
      t12209 = t12109*t1
      t12210 = 4*t12209
      t12211 = 32*t11240
      t12212 = t11086*sqq
      t12213 = 32*t12212
      t12214 = -t436+t540-t11182+t371-t12208+t435-t12210+t12211-t32+t122
     #13-t39
      t12215 = 1/t11277
      A1e0c3p1t1n2 = -24*t12207*t12214*t12215*t363*t571*t12199*t12203
      t12224 = t12207*t151
      t12227 = t12201*t11110
      t12228 = t12199*t11111*t12227
      t12229 = t7895*t12228
      A1e0c4p1t1n2 = 12*t596*t12161*t12224*t12229
      A1e0c5p1t1n2 = 0
      t12231 = t12177*CsThr
      t12234 = 128*t11347*t12231*t11531
      t12236 = 8*t12192*t33
      t12238 = 32*t12192*t11240
      t12240 = 2*t12192*t38
      t12243 = 32*t12192*CsThr*t370
      t12244 = t12192*t544
      t12246 = 24*t12192*t35
      t12249 = 2*t53*t12164*t32
      t12250 = t7928*Er
      t12252 = 64*t12250*t12167
      t12254 = 128*t11311*t12184
      t12256 = 128*t12250*t12178
      t12259 = 2*t53*t12171*t32
      t12262 = 64*t12250*t12171*t11079
      t12263 = t12234-t12236-t12238+t12240+t12243-t12244+t12246-t12249-t
     #12252-t12254-t12256+t12259-t12262
      t12266 = 16*t97*t12171*t31
      t12267 = t12171*BCsThlq
      t12269 = 128*t12250*t12267
      t12270 = t12171*sqq
      t12272 = 48*t44*t12270
      t12274 = 64*t12250*t12174
      t12276 = 64*t12192*t11234
      t12278 = 2*t12162*t38
      t12280 = 32*t12162*t11240
      t12282 = 8*t12162*t33
      t12284 = 64*t12162*t11234
      t12286 = 24*t12162*t35
      t12287 = t12164*sqq
      t12289 = 48*t44*t12287
      t12290 = t12164*BCsThlq
      t12292 = 128*t12250*t12290
      t12295 = 16*t97*t12164*t31
      t12298 = 64*t12250*t12164*t11079
      t12299 = -t12266+t12269+t12272+t12274+t12276-t12278+t12280+t12282-
     #t12284-t12286-t12289-t12292+t12295+t12298
      t12303 = 8*t97*t12290*t31
      t12306 = 64*t44*t12267*sqq
      t12307 = t12171*t37
      t12310 = 16*t97*t12307*sqq
      t12312 = 4*t12162*t12209
      t12314 = 8*t12162*t11181
      t12316 = 8*t12162*t334
      t12318 = 32*t12162*t12212
      t12320 = 8*t12192*t11416
      t12322 = 8*t12192*t334
      t12325 = 4*t53*t12307*t31
      t12328 = 16*t97*t12174*t31
      t12331 = 64*t12250*t12174*t11079
      t12334 = 2*t53*t12174*t32
      t12337 = 16*t44*t12270*t11079
      t12338 = -t12303-t12306-t12310-t12312-t12314+t12316+t12318+t12320-
     #t12322+t12325+t12328+t12331-t12334+t12337
      t12341 = 8*t97*t12267*t31
      t12344 = 48*t44*t12174*sqq
      t12347 = 128*t12250*t12174*BCsThlq
      t12350 = 32*t11287*t12184*t31
      t12353 = 128*t11311*t12184*t11079
      t12356 = 4*t11074*t12184*t32
      t12360 = 32*t12192*CsThr*t11079*t370
      t12362 = 2*t12192*t11831
      t12364 = 32*t12192*t12212
      t12366 = 4*t12192*t12209
      t12368 = 8*t12192*t11181
      t12371 = 24*t12192*t4553*sqq
      t12373 = 64*t12192*t11235
      t12376 = 4*t53*t12178*t32
      t12377 = t12341-t12344-t12347-t12350-t12353+t12356+t12360-t12362-t
     #12364+t12366+t12368-t12371-t12373+t12376
      t12381 = 128*t12250*t12181
      t12384 = 128*t12250*t12178*t11079
      t12387 = 96*t44*t12178*sqq
      t12390 = 256*t12250*t12178*BCsThlq
      t12393 = 32*t97*t12178*t31
      t12396 = 96*t11284*t12184*sqq
      t12399 = 256*t11311*t12184*BCsThlq
      t12400 = t12164*t37
      t12403 = 16*t97*t12400*sqq
      t12406 = 16*t44*t12287*t11079
      t12409 = 128*t12250*t12167*BCsThlq
      t12412 = 64*t44*t12290*sqq
      t12415 = 4*t53*t12400*t31
      t12418 = 48*t44*t12167*sqq
      t12419 = t12381-t12384+t12387+t12390-t12393+t12396+t12399+t12403-t
     #12406+t12409+t12412-t12415+t12418
      t12422 = 16*t97*t12167*t31
      t12425 = 64*t12250*t12167*t11079
      t12428 = 2*t53*t12167*t32
      t12429 = t97*t12164
      t12431 = 16*t12429*t294
      t12432 = t44*t12164
      t12433 = t297*t11079
      t12435 = 16*t12432*t12433
      t12437 = 8*t12429*t11307
      t12439 = 64*t12432*t11393
      t12440 = t44*t12171
      t12442 = 64*t12440*t11393
      t12444 = 4*t12193*t12209
      t12446 = 32*t12193*t12212
      t12448 = 8*t12193*t334
      t12450 = 8*t12193*t11181
      t12452 = 4*t12165*t11831
      t12455 = 8*t11075*t12231*t38
      t12456 = -t12422-t12425+t12428-t12431+t12435+t12437-t12439+t12442-
     #t12444+t12446+t12448-t12450+t12452-t12455
      t12460 = 16*t11353*t12231*t12109
      t12463 = 32*t11347*t12231*t11577
      t12466 = 32*t11353*t12231*t54
      t12469 = 32*t11353*t12231*t31
      t12472 = 128*t11318*t12231*t11079
      t12475 = 4*t11075*t12231*t32
      t12476 = t12177*BCsThlq
      t12479 = 128*t11347*t12476*sqq
      t12483 = 32*t11347*t12177*sqq*t11079
      t12484 = t12177*t37
      t12487 = 32*t11353*t12484*sqq
      t12490 = 8*t11075*t12484*t31
      t12493 = 16*t11353*t12476*t31
      t12495 = 4*t12172*t11831
      t12497 = 16*t12440*t12433
      t12498 = t97*t12171
      t12500 = 8*t12498*t11307
      t12501 = -t12460-t12463+t12466+t12469+t12472-t12475-t12479+t12483-
     #t12487+t12490+t12493-t12495-t12497-t12500
      t12503 = 16*t12498*t294
      t12504 = j4*t37
      t12507 = 32*t11353*t12504*sqq
      t12510 = 8*t11075*t12504*t31
      t12511 = j4*BCsThlq
      t12514 = 128*t11347*t12511*sqq
      t12518 = 32*t11347*j4*sqq*t11079
      t12521 = 16*t11353*t12511*t31
      t12524 = 96*t11347*t12231*sqq
      t12527 = 256*t11318*t12231*BCsThlq
      t12529 = 64*t12250*t12164
      t12530 = t12162*t32
      t12532 = 32*t12162*t370
      t12533 = t12192*t32
      t12535 = 32*t12192*t370
      t12537 = 64*t12250*t12171
      t12538 = t12503-t12507+t12510-t12514+t12518+t12521-t12524-t12527+t
     #12529-t12530+t12532+t12533-t12535-t12537
      t12544 = t12215*t151*t11110
      t12549 = t12200*t594*t598*t12201*t11111
      A1e0c6p1t1n2 = 24*t11283*(t12263+t12299+t12338+t12377+t12419+t1245
     #6+t12501+t12538)*t12544*t12549
      A1e0c7p1t1n2 = 12*t596*(t31-t11098-t12159+t11093+t12158+t12160-t11
     #094)*t12224*t12229
      A1e0c8p1t1n2 = 0
      t12555 = t12234-t12236-t12238-t12240+t12243-t12244+t12246-t12249-t
     #12252-t12254-t12256+t12259-t12262
      t12556 = -t12266+t12269+t12272+t12274+t12276+t12278+t12280+t12282-
     #t12284-t12286-t12289-t12292+t12295+t12298
      t12558 = -t12303-t12306+t12310-t12312-t12314-t12316+t12318+t12320+
     #t12322-t12325+t12328+t12331-t12334+t12337
      t12559 = t12341-t12344-t12347-t12350-t12353+t12356+t12360+t12362-t
     #12364+t12366+t12368-t12371-t12373+t12376
      t12562 = t12381-t12384+t12387+t12390-t12393+t12396+t12399-t12403-t
     #12406+t12409+t12412+t12415+t12418
      t12563 = -t12422-t12425+t12428+t12431+t12435+t12437-t12439+t12442-
     #t12444+t12446-t12448-t12450-t12452+t12455
      t12565 = -t12460-t12463-t12466+t12469+t12472-t12475-t12479+t12483+
     #t12487-t12490+t12493+t12495-t12497-t12500
      t12566 = -t12503+t12507-t12510-t12514+t12518+t12521-t12524-t12527+
     #t12529-t12530+t12532+t12533-t12535-t12537
      A1e0c9p1t1n2 = 24*t11283*(t12555+t12556+t12558+t12559+t12562+t1256
     #3+t12565+t12566)*t12544*t12549
      t12573 = 8*t12209
      t12575 = 16*t11086*t603
      t12577 = 8*t11080*t603
      t12579 = 16*t11086*Mzsq
      t12581 = 8*t11080*Mzsq
      t12583 = 4*t11531*t1799
      t12585 = 4*t2609*t11531
      t12588 = -t11530+t11546+t964+t435+t11505-t11555-t11556-16*t11181+t
     #11506+t11557-t11518
      t12591 = t11561*t11110
      t12594 = t12202*t12199*t11062
      A1e0c10p1t1n2 = 12*(-t12573+t816+t12213-t12575+t12577-t12579+t1258
     #1-t200+t12583+t12585+t12588)*t12194*t12591*t12594
      A1e0c11p1t1n2 = 0
      t12597 = 4*t11617*t33
      t12599 = 32*t675*t11531
      t12601 = 4*t11609*t816
      t12603 = 32*t11086*t656
      t12604 = -t3842+t11571-t11783-t12597+t12599-t12601+t12603-t11579-t
     #11587-t943+t11602+t11861+t11868+t11872-t11873-t11875
      t12606 = 64*t11234*t603
      t12608 = 32*t11240*t603
      t12609 = 64*t11907
      t12611 = 32*t1642*t11079
      t12613 = 128*t11234*sqq
      t12615 = 64*t6533*t11079
      t12617 = 8*t11609*t32
      t12619 = 64*t11086*t31
      t12620 = t11876+t11877-t1093-t4087-t12606+t12608-t12609+t12611+t12
     #613-t12615+t12617-t12619+t11751-t11948-t11949+t11950+t11951
      A1e0c12p1t1n2 = 24*(t12604+t12620)*t12194*BtQlq*t4620*t12215*t1178
     #0*t12228
      A1e0c1p1t2n2 = A1e0c1p1t1n2
      A1e0c2p1t2n2 = 0
      A1e0c3p1t2n2 = A1e0c3p1t1n2
      A1e0c4p1t2n2 = A1e0c4p1t1n2
      A1e0c5p1t2n2 = 0
      A1e0c6p1t2n2 = A1e0c6p1t1n2
      A1e0c7p1t2n2 = A1e0c7p1t1n2
      A1e0c8p1t2n2 = 0
      A1e0c9p1t2n2 = A1e0c9p1t1n2
      A1e0c10p1t2n2 = A1e0c10p1t1n2
      A1e0c11p1t2n2 = 0
      A1e0c12p1t2n2 = A1e0c12p1t1n2
      t12628 = 4*t11609
      A2e0c1p1t1n2 = 48*t964*(-t11961+sqq+t12628)*t12224*t12229
      A2e0c2p1t1n2 = 0
      t12633 = Er*t12171
      t12635 = 2*t403*t12633
      t12636 = BCsThlq*t43
      t12638 = 32*t12636*t12633
      t12640 = 16*t396*t12633
      t12641 = t31*t12191
      t12642 = j2*CsThr
      t12643 = t12641*t12642
      t12644 = t2*t12191
      t12646 = 16*t12644*t12642
      t12648 = 8*t865*t12192
      t12650 = 16*t11086*t12192
      t12651 = Er*t12164
      t12653 = 16*t396*t12651
      t12655 = 2*t403*t12651
      t12657 = 32*t12636*t12651
      t12659 = 32*t44*t12167
      t12661 = 16*t11086*t12162
      t12662 = t12635-t12638-t12640-t12643-t12646-t12648-t12650+t12653-t
     #12655+t12657+t12659+t12661
      t12664 = 8*t865*t12162
      t12666 = 64*t11284*t12184
      t12668 = 64*t44*t12178
      t12670 = 32*t44*t12174
      t12674 = 4*t11609*sqq*j2*j4
      t12678 = 4*t31*t11073*t53*j4
      t12679 = t11073*Er
      t12680 = t12679*j4
      t12682 = 32*t396*t12680
      t12686 = 4*t11609*sqq*t12191*j2
      t12688 = 16*t11086*t12193
      t12690 = 8*t865*t12193
      t12691 = t12679*t12177
      t12693 = 64*t12636*t12691
      t12695 = 4*t403*t12691
      t12697 = 64*t44*t12181
      t12698 = t12664+t12666+t12668-t12670-t12674+t12678-t12682+t12686+t
     #12688+t12690-t12693+t12695-t12697
      t12701 = 32*t396*t12691
      t12703 = 64*t12636*t12680
      t12705 = 4*t12157*t12193
      t12706 = BCsThlq*t19
      t12707 = t12706*sqq
      t12708 = t12651*CsThr
      t12710 = 8*t12707*t12708
      t12711 = t12633*CsThr
      t12713 = 8*t12707*t12711
      t12716 = 64*t12636*Er*t12181
      t12718 = 4*t141*t12181
      t12720 = 32*t143*t12181
      t12722 = 16*t12707*t12680
      t12724 = 16*t12707*t12691
      t12726 = 16*t396*t12711
      t12729 = 8*t12706*t317*t12171
      t12730 = -t12701-t12703-t12705+t12710-t12713+t12716-t12718+t12720+
     #t12722+t12724+t12726+t12729
      t12732 = 32*t12636*t12711
      t12734 = 2*t403*t12711
      t12737 = 16*t12707*t12679*t12231
      t12740 = 8*t12706*t317*t12164
      t12742 = 32*t12636*t12708
      t12744 = 2*t403*t12708
      t12746 = 16*t396*t12708
      t12747 = 32*t12432
      t12749 = 16*t11320*j4
      t12750 = t11463*j4
      t12752 = 16*t12644*j2
      t12753 = t12641*j2
      t12754 = 32*t12440
      t12755 = t12732-t12734-t12737-t12740-t12742+t12744-t12746-t12747-t
     #12749-t12750+t12752+t12753+t12754
      A2e0c3p1t1n2 = 96*t11969*(t12662+t12698+t12730+t12755)*t12215*t151
     #*t12229
      t12762 = 8*t12157
      t12764 = 4*t11609*t603
      t12766 = 4*t2609*BCsThlq
      A2e0c4p1t1n2 = 12*sqq*(-t12762+t656+t12764+t12766+t596+t12048-t120
     #49-t12050-t12054)*t12194*t12591*t12594
      A2e0c5p1t1n2 = 0
      t12772 = -t12583+t12575+t12579-t12585+t12573-t12213+t12144+t12145-
     #t964+t11518-t12106-t12107-t12129-t816+t12130
      A2e0c6p1t1n2 = 24*t12772*t12194*BtQlq*t12544*t12594
      A2e0c1p1t2n2 = A2e0c1p1t1n2
      A2e0c2p1t2n2 = 0
      A2e0c3p1t2n2 = A2e0c3p1t1n2
      A2e0c4p1t2n2 = A2e0c4p1t1n2
      A2e0c5p1t2n2 = 0
      A2e0c6p1t2n2 = A2e0c6p1t1n2
      t12777 = -t11094+t11098-t12158+t11093+t31+t12159+t12160
      t12778 = -t12162-t12166-t12169+t12173+t12176+t12180+t12183+t12186+
     #t12192+t12193
      t12782 = 1/(1+CsThr)
      t12783 = t12782*t149
      A1e0c1p1t1n3 = -12*t12777*t12778*t12196*t12783*t12203
      A1e0c2p1t1n3 = 0
      t12787 = t12778*BtQlq
      t12788 = t371-t12213+t12208-t11182+t12211+t12210-t32-t39+t540-t436
     #+t435
      t12789 = 1/t11272
      t12793 = t151*t12782
      t12794 = t12793*t149
      t12795 = t12794*t12203
      A1e0c3p1t1n3 = 24*t12787*t12788*t12789*t363*t12795
      t12798 = t12787*t151
      t12801 = t11986*t12227
      t12802 = t7894*t12782*t12801
      A1e0c4p1t1n3 = 12*t596*t12777*t12798*t12802
      A1e0c5p1t1n3 = 0
      t12804 = -t12234+t12236+t12238-t12240+t12243-t12244-t12246+t12249-
     #t12252+t12254+t12256-t12259+t12262
      t12805 = t12266+t12269-t12272+t12274+t12276+t12278-t12280-t12282-t
     #12284+t12286+t12289-t12292-t12295-t12298
      t12807 = -t12303-t12306+t12310-t12312+t12314-t12316+t12318+t12320+
     #t12322-t12325+t12328+t12331-t12334-t12337
      t12808 = t12341-t12344+t12347+t12350+t12353-t12356+t12360-t12362-t
     #12364+t12366-t12368-t12371+t12373-t12376
      t12811 = t12381+t12384-t12387+t12390+t12393-t12396+t12399-t12403+t
     #12406-t12409+t12412+t12415+t12418
      t12812 = -t12422-t12425+t12428-t12431+t12435-t12437+t12439-t12442+
     #t12444-t12446+t12448-t12450+t12452-t12455
      t12814 = t12460-t12463+t12466+t12469+t12472-t12475-t12479-t12483+t
     #12487-t12490+t12493-t12495-t12497+t12500
      t12815 = t12503+t12507-t12510-t12514-t12518+t12521-t12524+t12527-t
     #12529+t12530-t12532-t12533+t12535+t12537
      t12824 = t12783*t598*t594*t12789*t151
      A1e0c6p1t1n3 = -24*t11283*(t12804+t12805+t12807+t12808+t12811+t128
     #12+t12814+t12815)*t12203*t12824
      A1e0c7p1t1n3 = 12*t596*(t31-t11098+t12159+t11093-t12158+t12160-t11
     #094)*t12798*t12802
      A1e0c8p1t1n3 = 0
      t12831 = -t12234+t12236+t12238+t12240+t12243-t12244-t12246+t12249-
     #t12252+t12254+t12256-t12259+t12262
      t12832 = t12266+t12269-t12272+t12274+t12276-t12278-t12280-t12282-t
     #12284+t12286+t12289-t12292-t12295-t12298
      t12834 = -t12303-t12306-t12310-t12312+t12314+t12316+t12318+t12320-
     #t12322+t12325+t12328+t12331-t12334-t12337
      t12835 = t12341-t12344+t12347+t12350+t12353-t12356+t12360+t12362-t
     #12364+t12366-t12368-t12371+t12373-t12376
      t12838 = t12381+t12384-t12387+t12390+t12393-t12396+t12399+t12403+t
     #12406-t12409+t12412-t12415+t12418
      t12839 = -t12422-t12425+t12428+t12431+t12435-t12437+t12439-t12442+
     #t12444-t12446-t12448-t12450-t12452+t12455
      t12841 = t12460-t12463-t12466+t12469+t12472-t12475-t12479-t12483-t
     #12487+t12490+t12493+t12495-t12497+t12500
      t12842 = -t12503-t12507+t12510-t12514-t12518+t12521-t12524+t12527-
     #t12529+t12530-t12532-t12533+t12535+t12537
      A1e0c9p1t1n3 = -24*t11283*(t12831+t12832+t12834+t12835+t12838+t128
     #39+t12841+t12842)*t12203*t12824
      t12856 = t594*t12782*t149*t12203
      A1e0c10p1t1n3 = 12*(t12573+t816-t12213+t12575+t12577+t12579+t12581
     #-t200-t12583-t12585+t12588)*t12778*t11562*t12856
      A1e0c11p1t1n3 = 0
      t12858 = -t3842+t11571-t11783+t12597-t12599+t12601-t12603-t11579-t
     #11587-t943+t11602+t11861+t11868+t11872-t11873-t11875
      t12859 = t11876+t11877-t1093-t4087+t12606+t12608+t12609+t12611-t12
     #613-t12615-t12617+t12619+t11751-t11948-t11949+t11950+t11951
      t12863 = t4620*t12789
      A1e0c12p1t1n3 = -24*(t12858+t12859)*t12778*BtQlq*t12863*t4619*t127
     #82*t12801
      A1e0c1p1t2n3 = A1e0c1p1t1n3
      A1e0c2p1t2n3 = 0
      A1e0c3p1t2n3 = A1e0c3p1t1n3
      A1e0c4p1t2n3 = A1e0c4p1t1n3
      A1e0c5p1t2n3 = 0
      A1e0c6p1t2n3 = A1e0c6p1t1n3
      A1e0c7p1t2n3 = A1e0c7p1t1n3
      A1e0c8p1t2n3 = 0
      A1e0c9p1t2n3 = A1e0c9p1t1n3
      A1e0c10p1t2n3 = A1e0c10p1t1n3
      A1e0c11p1t2n3 = 0
      A1e0c12p1t2n3 = A1e0c12p1t1n3
      A2e0c1p1t1n3 = -48*t964*(t11961-sqq+t12628)*t12798*t12802
      A2e0c2p1t1n3 = 0
      t12874 = t12635+t12638-t12640+t12643+t12646-t12648+t12650+t12653-t
     #12655-t12657-t12659-t12661
      t12875 = t12664+t12666+t12668+t12670+t12674+t12678-t12682-t12686+t
     #12688-t12690+t12693+t12695+t12697
      t12877 = -t12701+t12703-t12705+t12710-t12713+t12716+t12718-t12720-
     #t12722-t12724-t12726-t12729
      t12878 = t12732+t12734-t12737+t12740-t12742-t12744+t12746-t12747-t
     #12749-t12750+t12752+t12753+t12754
      A2e0c3p1t1n3 = 96*t11969*(t12874+t12875+t12877+t12878)*t11110*t122
     #01*t11986*t12782*t7894*t12789*t151
      A2e0c4p1t1n3 = -12*sqq*(-t12762-t656+t12764+t12766-t596-t12048+t12
     #049+t12050+t12054)*t12778*t11562*t12856
      A2e0c5p1t1n3 = 0
      t12895 = t12579+t12575-t12213-t12583-t12585+t12573+t12129-t12130+t
     #816-t11518-t12145+t964+t12106+t12107-t12144
      A2e0c6p1t1n3 = 24*t12895*t12778*BtQlq*t12863*t12795
      A2e0c1p1t2n3 = A2e0c1p1t1n3
      A2e0c2p1t2n3 = 0
      A2e0c3p1t2n3 = A2e0c3p1t1n3
      A2e0c4p1t2n3 = A2e0c4p1t1n3
      A2e0c5p1t2n3 = 0
      A2e0c6p1t2n3 = A2e0c6p1t1n3
      t12900 = t12191*CsThr
      t12904 = log(DCMPLX(Eq/(Eq-Er)))
      t12905 = t12904*CsThr
      t12906 = j4-t12191+t12900-t12904+t12905
      t12907 = 4*t59
      t12908 = 64*t1344
      t12910 = 4*t75*t1
      t12911 = 16*t137
      t12912 = 12*t145
      t12913 = 8*t258
      t12914 = t34-t36+t56-t12907-t12908-t12910+t12911-t2248-t12912+t142
     #-t144+t12913
      t12915 = 32*t557
      t12916 = 8*t6863
      t12917 = 32*t1308
      t12918 = 8*t6666
      t12919 = t220*sqq
      t12920 = 16*t12919
      t12921 = 2*t5494
      t12922 = 8*t7949
      t12923 = 16*t41
      t12924 = -t12915+t348-t32-t39-t12916+t12917+t12918+t12920+t12921-t
     #12922+t12923+t45
      t12925 = t12914+t12924
      t12928 = t12199*t12201
      A1e0c1p1t1n4 = -12*t12906*t12925*t153*t571*t12928
      A1e0c2p1t1n4 = 0
      t12932 = 4*t55
      t12933 = 32*t1344
      t12934 = 8*t143
      t12935 = 16*t2856
      t12936 = t47*t370
      t12937 = 32*t12936
      t12938 = t435-t436+t12932-t12907+t12933-t12910-t355+t354-t12934+t1
     #2913-t12935+t12937
      t12939 = t370*BCsThl
      t12940 = 64*t12939
      t12941 = 16*t557
      t12942 = 16*t6666
      t12943 = 32*t12919
      t12944 = 16*t80
      t12945 = -t12940+t12941-t348+t371-t32-t39+t540-t12942+t12943+t1292
     #1-t12922+t12944
      t12948 = 1/t7885
      t12951 = t567*t149
      A1e0c3p1t1n4 = 24*t12906*(t12938+t12945)*t12948*t569*t12951*t151*t
     #12199*t12201
      t12962 = 1/(-t1233+t2469+t584*sll-t584*Mzsq)
      t12964 = t151*t12962*t12928
      A1e0c4p1t1n4 = 12*t12906*Mzsq*t12925*t149*t12964
      A1e0c5p1t1n4 = 0
      t12967 = 64*t1642*t1635
      t12969 = 32*t12936*t656
      t12971 = 128*t1642*t796
      t12973 = 32*t12936*t804
      t12975 = 4*t796*t3041
      t12977 = 8*t504*t3041
      t12979 = 32*t800*t5394
      t12981 = 64*t1642*t504
      t12983 = 16*t208*t67
      t12985 = 64*t1642*t2519
      t12987 = 4*t800*t3041
      t12989 = 64*t12936*t191
      t12991 = 64*t220*t208
      t12993 = 128*t12939*t191
      t12995 = 8*t67*t816
      t12997 = 64*t12939*t656
      t12999 = 32*t220*t816
      t13000 = t665+t12967-t12969-t12971-t12973+t12975-t12977-t12979+t12
     #981-t12983+t12985+t12987+t12989+t12991-t12993+t12995+t12997-t12999
      t13002 = 64*t12939*t804
      t13004 = 8*t675*t246
      t13006 = 32*t1713*t608
      t13008 = 32*t1642*t58
      t13009 = 16*t869
      t13010 = 16*t878
      t13011 = t2519*t19
      t13013 = 16*t13011*t4277
      t13015 = 8*t2509*t1328
      t13016 = t75*Er
      t13019 = 2*t13016*t5657*Mzsq
      t13021 = 4*t1847*t7871
      t13022 = t306*t603
      t13024 = 8*t13022*t318
      t13026 = 8*t1774*t7871
      t13027 = t181*t603
      t13029 = 16*t13027*t318
      t13032 = 16*t711*t47*t4357
      t13033 = t13002+t13004-t13006-t13008-t13009-t13010+t922-t945-t946+
     #t959+t13013-t13015-t13019-t13021-t13024+t13026+t13029+t13032
      t13036 = 2*t1891*t141
      t13039 = 32*t711*BCsThl*t4357
      t13041 = 4*t1862*t141
      t13042 = 2*t1094
      t13043 = 8*t1098
      t13044 = 32*t1101
      t13045 = 4*t1103
      t13046 = 4*t2129
      t13047 = 32*t1226
      t13048 = 16*t1228
      t13049 = -t13036-t13039+t13041-t13042+t13043+t13044+t13045+t1106-t
     #1108+t1148-t1150-t13046+t1152-t1170-t1172-t1189+t13047-t13048
      t13050 = 8*t1240
      t13051 = 4*t1242
      t13052 = 16*t1244
      t13053 = 2*t1249
      t13054 = 16*t1303
      t13055 = 4*t1434
      t13056 = 8*t1450
      t13058 = 64*t764*t9222
      t13060 = 8*t306*t8465
      t13062 = 32*t5531*t4357
      t13064 = 4*t796*t1331
      t13066 = 16*t675*t1639
      t13068 = 16*t711*t13016
      t13070 = 2*t800*t4591
      t13071 = -t13050+t13051-t13052-t13053-t13054+t1351-t1370+t1392+t14
     #06+t13055+t13056-t13058-t13060+t13062+t13064-t13066+t13068-t13070
      t13075 = 4*t504*t4591
      t13077 = 32*t243*t143
      t13080 = 16*t711*t243*Er
      t13082 = 16*t831*t143
      t13083 = t796*Er
      t13085 = 64*t764*t13083
      t13087 = 32*t1343*t4341
      t13089 = 32*t1343*t4357
      t13091 = 32*t220*t1609
      t13093 = 32*t1343*t4277
      t13095 = 32*t44*t1276
      t13097 = 32*t220*t1328
      t13099 = 8*t67*t1328
      t13101 = 8*t67*t1609
      t13103 = 32*t181*t4248
      t13105 = 16*t306*t4248
      t13106 = t246*Er
      t13108 = 8*t711*t13106
      t13110 = 64*t5531*t318
      t13111 = t58*Er
      t13113 = 32*t764*t13111
      t13114 = t13075-t13077+t13080+t13082+t13085-t13087-t13089-t13091-t
     #13093+t13095-t13097+t13099+t13101-t13103+t13105-t13108-t13110+t131
     #13
      t13116 = 16*t328*t4341
      t13118 = 16*t328*t4357
      t13120 = 16*t181*t8465
      t13122 = 32*t5531*t4341
      t13124 = 2*t796*t4591
      t13126 = 4*t663*t141
      t13128 = 64*t1343*t318
      t13129 = t2519*Er
      t13131 = 32*t764*t13129
      t13133 = 8*t800*t1754
      t13134 = t663*Er
      t13136 = 32*t711*t13134
      t13138 = 64*t675*t1888
      t13140 = 4*t75*t9215
      t13141 = 8*t2169
      t13142 = 4*t2172
      t13143 = 4*t1618
      t13144 = 8*t1621
      t13145 = t13116+t13118+t13120+t13122-t13124-t13126+t13128-t13131-t
     #13133-t13136+t13138+t13140+t13141-t13142-t1591+t1611+t13143-t13144
      t13159 = 32*t1633-16*t1637-16*t1640-t1669+4*t1673-16*t1675+4*t1689
     #-8*t1698+t1705+4*t1784+t1794-t1798+t1802-t1804-8*t1806+4*t1808+4*t
     #1814-16*t1824
      t13163 = -16*t1826+32*t1828+4*t1971+t1994+t1995-t1996-t1997-t1998-
     #t1999+t2001+t2002-t2004+t2005+t2006+t2007-t2009-t2010-t2012
      t13164 = t13159+t13163
      t13169 = t151*t2052*t2054
      t13173 = t2057*t12948*t12200*t12201
      A1e0c6p1t1n4 = 12*(t13000+t13033+t13049+t13071+t13114+t13145+t1316
     #4)*t12906*t13169*t13173
      t13177 = 32*t2625*BCsThl*t26
      t13179 = 32*t2640*t5272
      t13181 = 64*t2640*t5325
      t13183 = 16*t2625*t2519
      t13185 = 8*t2625*t58
      t13187 = 4*t2365*t33
      t13189 = 32*t2625*t133
      t13190 = 64*t2132
      t13191 = 128*t2134
      t13192 = 16*t2137
      t13193 = 12*t2139
      t13194 = 16*t2142
      t13196 = 8*t2365*t145
      t13198 = 8*t2613*t13129
      t13199 = t2066-t2067+t13177+t13179-t13181+t13183-t13185-t13187-t13
     #189-t2088+t2130-t2131+t13190-t13191+t13192-t13193+t13194-t13196+t1
     #3198
      t13201 = 8*t2613*t13111
      t13203 = 2*t2365*t141
      t13204 = 4*t2150
      t13205 = 32*t2194
      t13206 = 96*t2196
      t13207 = 48*t2198
      t13208 = 16*t2207
      t13209 = 16*t2187
      t13210 = t13201+t13203-t13204-t2154-t2163+t2170-t2173-t13205+t1320
     #6-t13207+t13208+t2210+t2212+t2214+t2215+t2216-t2217+t13209-t2220
      t13213 = t2227*t151
      A1e0c7p1t1n4 = 12*t12906*(t13199+t13210)*t13213*t149*t12962*t12928
      A1e0c8p1t1n4 = 0
      t13218 = 4*t2280
      t13219 = 12*t2414
      t13221 = 4*t2365*t3041
      t13223 = 64*t3355*t504
      t13225 = 32*t2365*t5394
      t13227 = 64*t3355*t2519
      t13229 = 8*t2625*t246
      t13231 = 32*t3355*t58
      t13233 = 32*t1713*t2206
      t13235 = 128*t3355*t796
      t13237 = 16*t13027*t10266
      t13238 = t13218+t13219+t13221+t13223-t13225+t13227+t13229-t13231-t
     #13233-t13235+t13237
      t13241 = 32*t2613*BCsThl*t4357
      t13243 = 4*t2769*t141
      t13245 = 2*t6938*t1358
      t13246 = t1343*t603
      t13248 = 32*t13246*t4317
      t13249 = t1343*sll
      t13251 = 64*t13249*t4317
      t13253 = 4*t6318*t1358
      t13255 = 32*t13246*t4277
      t13256 = t5531*sll
      t13258 = 64*t13256*t4317
      t13259 = t306*sll
      t13261 = 16*t13259*t10266
      t13262 = t181*sll
      t13264 = 32*t13262*t10266
      t13266 = 16*t13011*t5884
      t13267 = -t13241+t13243-t13245-t13248+t13251+t13253-t13255-t13258+
     #t13261-t13264+t13266
      t13271 = 2*t13016*t5657*t2063
      t13273 = 8*t2509*t2893
      t13275 = 8*t10805*t1328
      t13278 = 16*t2613*t47*t4357
      t13280 = 8*t13022*t10266
      t13281 = t5531*t603
      t13283 = 32*t13281*t4317
      t13285 = 4*t1891*t3548
      t13286 = t220*t603
      t13288 = 32*t13286*t1328
      t13290 = 32*t13281*t4277
      t13292 = 64*t3355*t1635
      t13293 = -t13271-t13273+t13275+t13278-t13280+t13283+t13285-t13288+
     #t13290+t13292+t2751
      t13299 = t2776-t2780+t2834-8*t2836-20*t2838-t2858+80*t2859+4*t2877
     #-t2882-48*t2883-t2888+t2928
      t13308 = -8*t2931+4*t2933-16*t2935-16*t2937+4*t2943-t2951+4*t2952+
     #t2972+t3036+t3038-t3040
      t13309 = t3043+t3045-t3047+t3049+t3051+t3058-t3060-t3062+t3064+t30
     #66+t3068
      t13310 = t13308+t13309
      t13314 = t3070-t3072+t3075+t3077+t3079-t3081-t3083-t3085-16*t3086+
     #4*t3088+8*t3090
      t13315 = 32*t3097
      t13316 = 32*t3101
      t13317 = 8*t3105
      t13319 = 32*t220*t1672
      t13321 = 4*t796*t944
      t13323 = 32*t12936*t1677
      t13325 = 16*t2625*t1639
      t13327 = 32*t12936*t1620
      t13329 = 64*t2640*t9222
      t13331 = 64*t12939*t1620
      t13333 = 64*t2625*t1888
      t13334 = t13315-t3100+t13316-t13317-t13319+t13321-t13323-t13325-t1
     #3327-t13329+t13331+t13333
      t13339 = 8*t67*t1672
      t13341 = 64*t12939*t1677
      t13343 = 8*t2365*t1754
      t13345 = 64*t220*t1589
      t13347 = 128*t12939*t1328
      t13349 = 64*t12936*t1328
      t13352 = 16*t208*t67*Mzsq
      t13354 = 32*t2640*t13111
      t13356 = 8*t2613*t13106
      t13358 = 16*t2613*t13016
      t13360 = 8*t504*t944
      t13361 = t13339+t13341-t13343+t13345-t13347+t13349-t13352+t13354-t
     #13356+t13358-t13360
      t13363 = 2*t2365*t4591
      t13365 = 16*t13027*t4115
      t13366 = 8*t3361
      t13367 = 4*t3363
      t13368 = 32*t3365
      t13369 = 16*t3372
      t13370 = 4*t3400
      t13371 = -t13363+t13365-t3360-t13366+t13367+t13368-t13369+t3375+t3
     #399+t13370+t3403
      t13373 = 32*t3404
      t13374 = 8*t3406
      t13375 = 16*t3408
      t13376 = 4*t3413
      t13377 = 16*t3418
      t13378 = 16*t3425
      t13380 = 8*t13022*t4115
      t13382 = 2*t1891*t3968
      t13383 = 16*t3499
      t13384 = t13373-t13374-t13375-t3411+t13376-t13377-t13378-t13380-t1
     #3382+t3094-t13383
      t13386 = 32*t2640*t13129
      t13389 = 4*t75*t262*t2063
      t13391 = 8*t67*t2893
      t13393 = 32*t1343*t5884
      t13395 = 32*t44*t3015
      t13397 = 32*t220*t2893
      t13399 = 64*t2640*t13083
      t13400 = 8*t3507
      t13401 = 32*t3524
      t13402 = 2*t3526
      t13403 = -t13386+t13389+t13391-t13393+t13395-t13397+t13399-t3506+t
     #13400+t3523+t13401-t13402
      t13408 = t3528+t3530+4*t3531-t3534+t3541-16*t3542-t3550-t3552-t355
     #6-t3558+t3560
      t13410 = -t3563+t3565+t3567-t3569-t3571+t3573+t3575-t3586-t3588+t3
     #596-16*t3597
      t13417 = -16*t3599+t3603-t3605-16*t3606-t3609-t3612+4*t3613+4*t361
     #5+t3625-t3733-2*t3736
      t13418 = -t3764-t3766+t3767-t3769-t3781-t3783-t3785+t3787-t3789-t3
     #790+t3808+t3820
      t13420 = t13408+t13410+t13417+t13418
      A1e0c9p1t1n4 = 12*(t13238+t13267+t13293+t13299+t13310+t13314+t1333
     #4+t13361+t13371+t13384+t13403+t13420)*t12906*t151*t3838*t13173
      t13428 = 8*t67*t656
      t13430 = 64*t1343*t4154
      t13432 = 4*t75*t1799
      t13434 = 32*t1635*t80
      t13436 = 32*t796*t80
      t13438 = 16*t220*t656
      t13439 = 64*t9386
      t13441 = 64*t67*t4573
      t13443 = 32*t800*t80
      t13445 = 32*t220*t191
      t13447 = 32*t608*t80
      t13448 = 32*t3933
      t13449 = 24*t3935
      t13450 = -t3842-t13428-t13430-t13432-t13434+t13436+t13438-t13439+t
     #13441+t13443-t13445-t13447-t3925-t943+t4407-t13448+t13449+t3938
      t13451 = 16*t3941
      t13452 = 16*t3943
      t13453 = 16*t3952
      t13454 = 32*t3954
      t13455 = 16*t3956
      t13456 = 12*t3964
      t13457 = 12*t3970
      t13458 = t181*t4317
      t13459 = 8*t13458
      t13461 = 2*t75*t1358
      t13463 = 16*t306*t318
      t13464 = -t3940-t13451+t13452+t3946-t3948+t4413+t4415-t13453+t1345
     #4+t13455-t13456-t3967+t3969-t13457+t4005+t13459+t13461-t13463
      t13466 = t181*t318
      t13467 = 16*t13466
      t13468 = t2519*t336
      t13469 = 16*t13468
      t13470 = t796*t145
      t13471 = 8*t13470
      t13472 = t181*t4357
      t13473 = 8*t13472
      t13475 = 2*t75*t2854
      t13476 = t800*t145
      t13477 = 8*t13476
      t13479 = 4*t75*t6908
      t13480 = 8*t4472
      t13481 = 128*t9729
      t13483 = 4*t75*t2609
      t13485 = 64*t1343*t4276
      t13487 = 8*t75*t262
      t13489 = 64*t5531*t4114
      t13491 = 16*t2519*t675
      t13493 = 8*t58*t675
      t13494 = 4*t4008
      t13495 = -t13467+t13469-t13471+t13473+t13475-t13477-t13479+t13480+
     #t13481-t13483-t13485+t13487-t13489+t13491-t13493-t13494+t4022+t403
     #2
      t13502 = -4*t4033+8*t4061-t4066+8*t4503-t1093-t4087-32*t4080-t4090
     #+16*t4047+16*t4084+t4094-t4095-t4096+t4097+t4098+t4099+t4100+t4101
      t13510 = t569/(t152-t582-t804+t791)*t12928
      A1e0c10p1t1n4 = 12*t12906*(t13450+t13464+t13495+t13502)*t571*t1351
     #0
      A1e0c11p1t1n4 = 0
      t13513 = 16*t1211*t6159
      t13515 = 128*t12936*t1436
      t13517 = 8*t815*t2751
      t13520 = 4*t815*t791*t1436
      t13523 = 2*t815*t4759*t208
      t13525 = t815*t2300*t208
      t13526 = 2*t13525
      t13533 = t5616-12*t5748-8*t8378-64*t5777-t5837+t5839+t5841+16*t585
     #9+16*t5861-t5865-t8521
      t13536 = 6*t5895
      t13539 = 6*t5905
      t13548 = 2*t5979
      t13549 = 4*t5981
      t13553 = 48*t5949-32*t5951-12*t5953+48*t5955+16*t5959+t13548-t1354
     #9-8*t5984+16*t5986+16*t5989-t5993
      t13554 = -32*t5887+t13536+16*t5897-t5902-24*t5903+t13539+t5908-8*t
     #5914-t5921+8*t5947+t13553
      t13556 = 4*t6030
      t13568 = t6053+t6055-32*t6056-4*t6058+16*t6060+48*t6062+8*t6064+48
     #*t6066-12*t6069-t6082+t6084
      t13569 = t5995+t8531-t13556-32*t6033+16*t6035-24*t6037+t6040-t6042
     #-t6044+t6048+t13568
      t13573 = t6086-t6088-4*t6089+8*t6092-12*t6095+t6098-t6101+t6106-t6
     #108-t6110
      t13574 = 16*t6211
      t13575 = 8*t6213
      t13576 = 32*t6275
      t13577 = BCsThl*t905
      t13580 = 12*t582*t13577*t1
      t13582 = 96*t1713*t1221
      t13584 = 96*t12936*t1298
      t13586 = t2469*t10710*t1
      t13587 = 8*t13586
      t13589 = 64*t1713*t2149
      t13591 = 12*t1211*t944
      t13592 = -t6210+t13574-t13575+t6219+t13576-t13580+t13582+t13584+t1
     #3587-t13589-t13591
      t13598 = 128*t1713*t4869*sqq
      t13600 = 64*t1713*t2769
      t13602 = 192*t5871*t1888
      t13604 = 128*t5871*t1862
      t13605 = t206*t6938
      t13606 = 96*t13605
      t13608 = t206*t4869*t31
      t13609 = 64*t13608
      t13611 = 64*t1713*t4870
      t13613 = 24*t4994*t1688
      t13614 = t2365*t1652
      t13615 = 4*t13614
      t13617 = 32*t6533*t2851
      t13620 = 24*t206*t6944
      t13622 = 4*t10710*t6251
      t13624 = 192*t1713*t1694
      t13626 = 16*t206*t7380
      t13627 = t4994*t1694
      t13628 = 96*t13627
      t13630 = 64*t5178*t8270
      t13632 = 24*t4994*t1221
      t13633 = t5272*t32
      t13635 = 32*t5178*t13633
      t13636 = t4994*t1218
      t13637 = 96*t13636
      t13639 = 192*t1713*t1218
      t13641 = 128*t1713*t3123
      t13642 = -t13620-t13622-t13624+t13626+t13628+t13630-t13632-t13635+
     #t13637-t13639+t13641
      t13644 = t5394*t4870
      t13645 = 32*t13644
      t13646 = t6517*t1862
      t13647 = 64*t13646
      t13649 = 8*t796*t6160
      t13651 = 64*t12936*t1672
      t13653 = 32*t6533*t6898
      t13655 = 8*t5394*t2851
      t13657 = 16*t5192*t8273
      t13658 = t5394*t2769
      t13659 = 32*t13658
      t13661 = 16*t206*t3107
      t13663 = 24*t206*t7127
      t13665 = t206*t2956
      t13666 = 64*t13665
      t13667 = t206*t6364
      t13668 = 96*t13667
      t13669 = t5325*t32
      t13671 = 64*t5178*t13669
      t13673 = 192*t5871*t1314
      t13675 = 128*t4853*t13016
      t13677 = 128*t4853*t13106
      t13680 = 12*t504*t959*t1
      t13682 = 96*t5871*t1639
      t13685 = 8*t504*t6514*t1
      t13687 = 12*t10710*t6195
      t13688 = -t13666+t13668+t13671-t13673-t13675+t13677-t6398-t13680+t
     #13682+t13685-t13687
      t13692 = 64*t6533*t4870
      t13694 = 64*t6533*t2769
      t13696 = 8*t5394*t6898
      t13698 = 16*t6517*t1632
      t13701 = 8*t315*t5899*t905
      t13703 = 128*t6481*t1862
      t13705 = 32*t1713*t6898
      t13707 = 32*t1713*t2851
      t13709 = 96*t1713*t1688
      t13718 = -64*t1713*t5677*sqq-64*t12936*t2830-t6470-t6472+16*t6473-
     #64*t6475-t6478-t6480-t6483-t6485+t6487
      t13721 = t6509+t6511-t6513+t6516+t6519-t6521+t6523+t6525-t6527+t65
     #29+t6535
      t13722 = t6489-t6491-t6494+t6496+t6498-t6500-t6501+t6503-t6505+t65
     #07+t13721
      t13726 = 32*t6540
      t13727 = 8*t6542
      t13728 = t280*Mzsq
      t13730 = 48*t13728*t8270
      t13731 = t280*t2063
      t13733 = 32*t13731*t533
      t13734 = t5178*t37
      t13736 = 48*t13734*t4317
      t13737 = t4853*BCsThl
      t13739 = 96*t13737*t4317
      t13741 = 256*t5871*t75
      t13744 = 6*t815*t3768
      t13746 = 6*t815*t9563
      t13747 = t815*t3808
      t13748 = 4*t13747
      t13750 = 96*t13256*t10266
      t13752 = 64*t13256*t10270
      t13754 = 96*t13249*t10266
      t13755 = t13249*t10270
      t13756 = 64*t13755
      t13757 = t5178*t47
      t13759 = 24*t13757*t10266
      t13761 = 32*t4994*t4364
      t13762 = t4994*t10710
      t13763 = 128*t13762
      t13764 = t13744+t13746-t13748+t6660-t13750+t13752+t13754-t13756+t1
     #3759+t13761-t13763
      t13767 = 2*t815*t7309
      t13768 = t815*t3780
      t13769 = 2*t13768
      t13772 = 6*t815*t804*t905
      t13775 = 6*t815*t1095*t32
      t13778 = 4*t815*t4940*t32
      t13780 = 4*t815*t3820
      t13782 = 48*t13728*t13669
      t13784 = 32*t13731*t13016
      t13787 = 24*t13728*t13633
      t13789 = 16*t13731*t13106
      t13790 = t4853*t47
      t13792 = 96*t13790*t4317
      t13793 = t5178*BCsThl
      t13795 = 48*t13793*t10266
      t13797 = 4*t1215*t6251
      t13798 = t5120*t1672
      t13799 = 128*t13798
      t13800 = t5120*t2761
      t13801 = 32*t13800
      t13803 = 8*t326*t2761
      t13804 = t5120*t4838
      t13805 = 32*t13804
      t13806 = t5657*t603
      t13808 = 12*t13806*t8273
      t13809 = t5657*t4628
      t13811 = 8*t13809*t8270
      t13812 = t13787-t13789-t13792-t13795-t13797-t13799+t13801-t13803+t
     #13805+t13808-t13811
      t13815 = t4154*t32
      t13817 = 24*t13259*t13815
      t13818 = t5912*t31
      t13820 = 16*t13259*t13818
      t13822 = 16*t5644*t13083
      t13823 = t5268*Mzsq
      t13825 = 16*t13823*t9222
      t13827 = 8*t5644*t9222
      t13828 = t6620*BCsThl
      t13829 = t13828*t5983
      t13830 = 32*t13829
      t13834 = 4*t403*t4628*t532*Mzsq
      t13836 = 8*t1691*t3548
      t13839 = 128*t1713*t47*t1620
      t13841 = 16*t1891*t944
      t13843 = t3120*t1754
      t13844 = 4*t13843
      t13846 = 32*t13828*t5913
      t13847 = t6620*t47
      t13849 = 32*t13847*t5913
      t13852 = 16*t396*t37*t5913
      t13855 = 64*t5112*BCsThl*t5980
      t13858 = 64*t5112*t47*t5980
      t13859 = t689*BCsThl
      t13861 = 16*t13859*t5913
      t13864 = 16*t396*t2063*t13134
      t13866 = 32*t13847*t5983
      t13868 = 32*t13823*t13083
      t13871 = 16*t4994*t47*t1620
      t13872 = -t13844+t13846-t13849-t13852-t13855+t13858-t13861-t13864-
     #t13866+t13868+t13871
      t13874 = t5123*t1620
      t13875 = 64*t13874
      t13877 = 48*t13734*t4357
      t13879 = 32*t13734*t5980
      t13881 = 96*t13790*t4357
      t13883 = 64*t13790*t5980
      t13884 = t280*t37
      t13886 = 16*t13884*t5913
      t13889 = 8*t689*t47*t5913
      t13892 = 32*t689*t37*t5980
      t13895 = 12*t5192*t37*t8465
      t13897 = 32*t326*t1672
      t13900 = 32*t13256*t5983
      t13902 = 96*t13249*t8465
      t13904 = 48*t13884*t8465
      t13906 = 48*t13793*t8465
      t13908 = 24*t13757*t8465
      t13910 = 16*t13731*t13134
      t13912 = 96*t13737*t4357
      t13914 = 256*t5131*t1620
      t13915 = t5912*sqq
      t13917 = 32*t13884*t13915
      t13919 = 32*t13256*t5913
      t13921 = 64*t13256*t13915
      t13922 = -t13900+t13902-t13904-t13906+t13908-t13910+t13912+t13914+
     #t13917-t13919+t13921
      t13927 = 64*t13737*t5980
      t13929 = 32*t13249*t5913
      t13931 = 64*t13249*t13915
      t13933 = 8*t326*t4838
      t13934 = t13249*t5983
      t13935 = 32*t13934
      t13937 = 96*t13256*t8465
      t13939 = 48*t13262*t13815
      t13941 = 32*t13262*t13818
      t13942 = 4*t8727
      t13943 = 8*t8730
      t13948 = 4*t6933
      t13950 = -t6824+t6871-t6876-t8903+16*t6892+16*t6899-32*t6902+t6910
     #-t6928-t13948+8*t6983
      t13960 = -64*t6985+8*t6987-4*t6989-32*t6991+t7014+16*t7045+t8974-1
     #2*t7079+48*t7133-12*t7139
      t13961 = 48*t7147
      t13962 = 4*t7193
      t13963 = 32*t7220
      t13964 = 4*t7239
      t13965 = 8*t7241
      t13968 = 64*t13728*t5325*t816
      t13971 = 16*t13731*t5325*t656
      t13974 = 32*t13728*t5272*t816
      t13977 = 8*t13731*t5272*t656
      t13981 = 4*t5657*t37*t5405*t656
      t13984 = 4*t13809*t532*t596
      t13985 = t13961-t13962-t13963-t13964+t13965+t13968-t13971-t13974+t
     #13977+t13981+t13984
      t13991 = t5912*t596
      t13994 = t4154*t596
      t14012 = -t7415+48*t7421-4*t7423-12*t7425-t7431+48*t7438-4*t7446+8
     #*t7448+16*t7472+16*t7511-t7572
      t14013 = -16*t13806*t532*t964+8*t13259*t13991+64*t13884*t13994-64*
     #t7296-8*t9025+12*t9027+12*t9029-12*t7315+16*t7319+t7384+t14012
      t14015 = 128*t13249*t13994
      t14017 = 128*t13256*t13994
      t14019 = 16*t13262*t13991
      t14021 = 16*t13757*t13994
      t14023 = 32*t13793*t13994
      t14026 = -t7701-t7703-t7705+t7707-t7711+t7713-t7714-t7723-32*t7731
     #-t7736+t7766
      t14034 = t4620*t569
      t14038 = t12928*t567*t571*t4108
      s2 = -12
      s5 = t6548+t13604+t13606+t13598-t13600-t13602+t13592+t13872+t13836
     #-t13839+t13841+t13830+t13834+t13573+t13817-t13820-t13822-t13825+t1
     #3827+t13569+t13812+t13554+t13526+t13533+t13772+t13775-t13778-t1378
     #0+t13523-t13782+t13784+t13513-t13515-t13517-t13520+t13764+t13767+t
     #13769+t13741+t13739+t14026-t14021+t14023-t14015+t14017-t14019+t140
     #13+t13985+t13960+t13950+t13935-t13937-t13939+t13941+t13942-t13943-
     #t13933+t13922-t13927
      s4 = t13929-t13931+t13897+t13895+t13889+t13892-t13886-t13877+t1387
     #9-t13881+t13883-t13875+t6698-t5389+t4743-t4736-t13736+t13733-t1373
     #0+t13726-t13727+t13722+t13709+t13718+t13705+t13707-t13701+t13703+t
     #6539-t13696+t13698-t6537-t13694-t13692+t13688+t13661-t13663+t7700+
     #t7697+t7679-t7694-t7695+t13659-t13655-t13657-t13651+t13653+t13649+
     #t13645-t13647+t13642+t6434-t13609-t13611-t13613-t13615+t13617+t488
     #0+s5+t6683
      s5 = t12906
      s3 = s4*s5
      s1 = s2*s3
      s2 = t12948*t14034*t14038
      A1e0c12p1t1n4 = s1*s2
      A1e0c1p1t2n4 = A1e0c1p1t1n4
      A1e0c2p1t2n4 = 0
      A1e0c3p1t2n4 = A1e0c3p1t1n4
      A1e0c4p1t2n4 = A1e0c4p1t1n4
      A1e0c5p1t2n4 = 0
      A1e0c6p1t2n4 = A1e0c6p1t1n4
      A1e0c7p1t2n4 = A1e0c7p1t1n4
      A1e0c8p1t2n4 = 0
      A1e0c9p1t2n4 = A1e0c9p1t1n4
      A1e0c10p1t2n4 = A1e0c10p1t1n4
      A1e0c11p1t2n4 = 0
      A1e0c12p1t2n4 = A1e0c12p1t1n4
      t14041 = 24*t6666
      t14042 = 16*t6863
      t14043 = 32*t1352
      t14044 = -t7840+t7837+t7838+t32-t142-t34-t14041+t14042-t12921-t140
     #43+t12917+t12910
      A2e0c1p1t1n4 = -48*t12906*t14044*Mzsq*t149*t12964
      A2e0c2p1t1n4 = 0
      t14051 = 32*t75*t8012
      t14052 = t7898*BCsThl
      t14054 = 8*t14052*t3041
      t14055 = t1215*t8776
      t14056 = 32*t14055
      t14057 = 4*t13525
      t14058 = -t14051+t14054-t14056-t8054-t14057+t8088+t8091-t8093-t810
     #2+t8104+t8106+t8109+t8111-t8118+t8122-t8140-t8144-t8152+t8172
      t14061 = -t8177+t8186-t8235-t8240-4*t8248+t8277-t8327-t8332-t8334+
     #t8337+t8347+t8357+t8371+t8394-t8396-t8460-t8464+8*t8466-t8469-t847
     #3
      t14063 = 8*t5885
      t14064 = 4*t6029
      t14067 = 64*t8012*t2519*t26
      t14068 = t32*t26
      t14071 = 8*t220*t14068*Mzsq
      t14073 = t905*t26*t1
      t14075 = 2*t504*t14073
      t14078 = 16*t504*t6517*t26
      t14080 = 32*t7907*t1218
      t14082 = 8*t14052*t1331
      t14084 = 4*t14052*t4591
      t14085 = -t8477+t14063-t8522+t8523+t5953-t8526+t8529-t14064-t6057-
     #t8532+t8533+t6069-t14067-t14071-t14075+t14078+t14080+t14082-t14084
      t14087 = t220*t1436*t26
      t14088 = 16*t14087
      t14091 = 32*t12939*t208*t26
      t14092 = t939*t26
      t14094 = 32*t12939*t14092
      t14096 = 4*t41*t6364
      t14097 = t1211*t4962
      t14098 = 2*t14097
      t14099 = t31*t26
      t14100 = t14099*Mzsq
      t14102 = 16*t12939*t14100
      t14106 = 16*t12939*t152*t26*Mzsq
      t14109 = 64*t8012*t504*t26
      t14110 = t816*t26
      t14112 = 16*t12939*t14110
      t14115 = 8*t220*t1093*t26
      t14117 = 64*t8107*t13129
      t14120 = 64*t8107*t504*Er
      t14122 = 4*t41*t1589
      t14124 = 64*t7929*t1314
      t14125 = t7928*BCsThl
      t14126 = t3958*Mzsq
      t14128 = 64*t14125*t14126
      t14130 = 64*t7929*t1888
      t14132 = 16*t14052*t1754
      t14134 = 128*t7929*t1862
      t14135 = t796*t14073
      t14138 = 128*t8012*t796*t26
      t14139 = -t14088-t14091+t14094-t14096+t14098+t14102-t14106-t14109+
     #t14112-t14115+t14117+t14120+t14122-t14124-t14128-t14130-t14132+t14
     #134+t14135+t14138
      t14144 = 16*t12939*t1095*t26
      t14147 = 8*t504*t5268*t40
      t14150 = 64*t1343*t208*t40
      t14152 = 64*t8012*t1888
      t14154 = 4*t41*t1667
      t14156 = 128*t8107*t13083
      t14158 = 32*t7907*t1694
      t14160 = 4*t41*t6938
      t14162 = 64*t7929*t1891
      t14164 = 16*t14052*t4589
      t14165 = t13577*t1257
      t14167 = 32*t8012*t1270
      t14170 = 128*t14125*t191*t40
      t14171 = 16*t13586
      t14172 = 8*t13614
      t14173 = 2*t6473
      t14174 = -t14144+t14147-t14150+t14152+t14154-t14156+t14158-t14160-
     #t14162+t14164+t14165-t14167+t14170-t14171+t14172-t13659+t13666+t85
     #86-t14173
      t14175 = t5123*t3973
      t14176 = 8*t14175
      t14177 = t5120*t2115
      t14178 = 128*t14177
      t14179 = t14052*t603
      t14181 = 8*t14179*t141
      t14183 = t1211*t31*t1257
      t14184 = t14052*sll
      t14186 = 16*t14184*t145
      t14187 = t6938*t1257
      t14188 = 2*t14187
      t14189 = t12939*t603
      t14191 = 32*t14189*t3973
      t14193 = 64*t7940*t1620
      t14195 = 4*t1352*t1667
      t14196 = t13286*t14092
      t14197 = 8*t14196
      t14199 = 32*t14179*t145
      t14201 = t13286*t4573*t2063
      t14202 = 64*t14201
      t14203 = t13286*t2115
      t14204 = 64*t14203
      t14206 = 32*t14189*t4574
      t14207 = t8590+t8597+t8598-t8599+t8600-t8601+t14176+t14178+t14181-
     #t14183+t14186-t14188-t14191-t14193-t14195+t14197-t14199-t14202-t14
     #204+t14206
      t14209 = t13286*t14100
      t14210 = 16*t14209
      t14211 = t1691*t4589
      t14213 = 32*t7940*t1609
      t14214 = t1694*t8776
      t14215 = 16*t14214
      t14217 = 8*t1352*t1672
      t14219 = 4*t14184*t141
      t14221 = 4*t13859*t14126
      t14223 = 32*t13828*t14126
      t14224 = 8*t13747
      t14225 = 128*t13755
      t14226 = 4*t13768
      t14227 = 64*t13829
      t14228 = 8*t13843
      t14229 = 64*t13934
      t14230 = 2*t6892
      t14231 = t14210-t14211+t14213+t14215+t14217-t14219-t14221+t14223-t
     #8665+t14224-t14225-t14226-t13801+t14227+t14228+t8724+t14229-t14230
     #+t8906
      t14232 = t1891*t1242
      t14233 = 2*t14232
      t14234 = t1891*t2139
      t14235 = 16*t14234
      t14237 = t13286*t191*t3972
      t14238 = 16*t14237
      t14239 = t8938+t8965-t9039-t9049+t9058-t9060+t9061-t9062-t9064+t90
     #65+t9069-t9071+t9073-t9088+t9090+t9092+t14233+t14235-t14238-t9095
      A2e0c3p1t1n4 = 24*(t14058+t14061+t14085+t14139+t14174+t14207+t1423
     #1+t14239)*t12906*t2052*t9104*t12948*t149*t151*t9106*t12199*t12201
      t14252 = 64*t5773
      t14254 = 32*t181*t3961
      t14256 = 32*t181*t3930
      t14257 = t14252-t14254-t14256+t13432+t13436-t13439+t13443-t9134+t9
     #43-t9136-t9137+t3940-t3946+t9138+t9139+t9140-t3969+t9143
      t14258 = 24*t13458
      t14259 = 48*t13466
      t14260 = 32*t13468
      t14261 = 16*t13470
      t14262 = 24*t13472
      t14263 = 16*t13476
      t14264 = -t4005-t14258-t13461+t14259-t14260+t14261-t14262-t13475+t
     #14263+t13479+t13483-t13487+t1093+t4096+t9151-t4097-t4099-t9153
      A2e0c4p1t1n4 = -12*t12906*(t14257+t14264)*t571*t13510
      A2e0c5p1t1n4 = 0
      t14271 = 128*t4853*t8963
      t14272 = 48*t4735
      t14273 = 16*t14055
      t14276 = 96*t4819*t8465
      t14278 = 32*t4819*t5913
      t14280 = 64*t4819*t13915
      t14282 = 96*t4853*t4357
      t14284 = 64*t4853*t5980
      t14286 = 32*t6620*t5913
      t14288 = 64*t5112*t5980
      t14289 = 64*t8276
      t14291 = 96*t5120*t14100
      t14292 = t14068*BCsThl
      t14294 = 24*t9215*t14292
      t14296 = 24*t6190*t14292
      t14297 = -t14276-t14278+t14280-t14282+t14284-t14286+t14288+t14289-
     #t14291+t14294+t14296
      t14301 = 16*t6248*t14099*BCsThl
      t14302 = t4628*t26
      t14305 = 64*t5120*t14302*sqq
      t14306 = t14302*Mzsq
      t14308 = 32*t5120*t14306
      t14311 = 128*t4819*Mzsq*t4357
      t14313 = 96*t5120*t14110
      t14316 = 64*t5123*t4524*Mzsq
      t14318 = 32*t12919*t14306
      t14322 = 8*t33*t4628*t3972*BCsThl
      t14323 = 32*t8468
      t14324 = 32*t8472
      t14335 = -96*t4853*t4317-96*t4819*t10266+t8521-t13536-32*t5897+48*
     #t5901+48*t5903-t13539-t5908+16*t5914+16*t5920
      t14343 = -t13548+t13549+16*t5984-32*t5986-32*t5989+48*t5992-t5995-
     #t8531+t13556+48*t6037
      t14344 = 16*t6081
      t14345 = 32*t6085
      t14346 = 64*t6097
      t14347 = 128*t14087
      t14348 = 32*t14097
      t14349 = 16*t6209
      t14350 = 32*t6211
      t14351 = t14344-t14345+t6088-t14346-t6106+t6108+t14347-t14348+t143
     #49-t14350+t13575
      t14353 = 64*t6275
      t14354 = 48*t13605
      t14355 = 32*t13608
      t14356 = 48*t13627
      t14357 = 48*t13636
      t14359 = 16*t13644
      t14360 = 32*t13646
      t14361 = 16*t13658
      t14362 = 32*t13665
      t14363 = 48*t13667
      t14364 = 16*t6397
      t14365 = -t14359+t14360-t13649-t14361+t14362-t14363-t13675+t14364+
     #t13680-t13685+t13687
      t14369 = 48*t6471
      t14370 = 16*t6484
      t14371 = 48*t6493
      t14372 = 16*t6499
      t14379 = -32*t6506-t6511+48*t6512-t6516-32*t6518+t6521-t6529-t6539
     #-64*t6540+t13727-32*t6547
      t14391 = -96*t14175+64*t14177+24*t14183-16*t14187-96*t14196-32*t14
     #201-32*t14203+64*t14209+24*t14211+8*t14214
      t14392 = 64*t13762
      t14393 = t13739-t13744-t13746+t13748+t13754-t13756+t14392-t6683-t1
     #3767-t13769-t13772
      t14399 = -t13775+t13778+t13780+t13797+64*t13798-16*t13800-16*t1380
     #4+t13830-t13836-t13841
      t14400 = 32*t13874
      t14401 = t13844+t13846-t13855+t14400+t13902+t13912-t13927+t13929-t
     #13931+t13935-t6871
      t14403 = 64*t7696
      t14408 = 8*t6248*t26*sqq*t800
      t14410 = 16*t7058*t801
      t14411 = 32*t14232
      t14412 = 8*t14234
      t14413 = 128*t14237
      t14414 = t7703+t7705+t7711-t7713+t7714+t7723+t14408-t14410-t14411+
     #t14412+t14413
      A2e0c6p1t1n4 = 12*(t13591-t13587+t13580+t6928-t6910+t6876-t13526-t
     #13523-t13513+t13517+t13520-t14015-t14356-t14357+t13948-t14353-t143
     #54+t14355+t14351+t14343+t14335-t14323-t14324-t14318+t14322-t14313+
     #t14316-t14308+t14311-t14301+t14305+t14297+t14271+t14272-t14273-t73
     #84-t4743+t13701+t6505-t6503+t6501+t7701-t7700+t7694-t6496+t13615+t
     #13622-t4880-t6489+t14414+t14401-t14403+t14399+t14393+t14391+t14379
     #+t14371+t14372+t14369+t14370+t14365)*t12906*t12948*t14034*t14038
      A2e0c1p1t2n4 = A2e0c1p1t1n4
      A2e0c2p1t2n4 = 0
      A2e0c3p1t2n4 = A2e0c3p1t1n4
      A2e0c4p1t2n4 = A2e0c4p1t1n4
      A2e0c5p1t2n4 = 0
      A2e0c6p1t2n4 = A2e0c6p1t1n4
      t14423 = -j4+t12191+t12900+t12904+t12905
      t14424 = t34-t36+t56-t12907+t12908+t12910+t12911-t2248-t12912+t142
     #-t144+t12913
      t14425 = -t12915+t348-t32-t39+t12916-t12917-t12918-t12920-t12921-t
     #12922+t12923+t45
      t14426 = t14424+t14425
      t14429 = t12782*t12201
      A1e0c1p1t1n5 = -12*t14423*t14426*t153*t571*t14429
      A1e0c2p1t1n5 = 0
      t14433 = t435-t436+t12932-t12907-t12933+t12910-t355+t354-t12934+t1
     #2913-t12935+t12937
      t14434 = t12940+t12941-t348+t371-t32-t39+t540+t12942-t12943-t12921
     #-t12922+t12944
      t14437 = 1/t7890
      A1e0c3p1t1n5 = -24*t14423*(t14433+t14434)*t14437*t569*t12951*t1279
     #3*t12201
      t14445 = t14423*t580
      t14448 = t586*t12201*t12794
      A1e0c4p1t1n5 = 12*t14426*Mzsq*t14445*t14448
      A1e0c5p1t1n5 = 0
      t14450 = t665+t12967-t12969+t12971-t12973-t12975+t12977+t12979-t12
     #981-t12983-t12985-t12987+t12989-t12991+t12993+t12995-t12997+t12999
      t14451 = -t13002+t13004-t13006-t13008-t13009-t13010+t922-t945-t946
     #+t959-t13013-t13015+t13019-t13021-t13024+t13026-t13029+t13032
      t14453 = t13036+t13039-t13041-t13042+t13043+t13044+t13045+t1106-t1
     #108+t1148-t1150-t13046+t1152-t1170-t1172-t1189+t13047-t13048
      t14454 = -t13050+t13051-t13052-t13053-t13054+t1351-t1370+t1392+t14
     #06+t13055+t13056-t13058-t13060+t13062-t13064-t13066-t13068+t13070
      t14457 = -t13075-t13077+t13080+t13082-t13085+t13087+t13089+t13091+
     #t13093+t13095+t13097+t13099+t13101+t13103+t13105-t13108-t13110+t13
     #113
      t14458 = t13116+t13118-t13120+t13122+t13124-t13126-t13128+t13131+t
     #13133-t13136-t13138-t13140+t13141-t13142-t1591+t1611+t13143-t13144
      t14466 = t2057*t14437*t12783*t12201
      A1e0c6p1t1n5 = -12*(t14450+t14451+t14453+t14454+t14457+t14458+t131
     #64)*t14423*t13169*t14466
      t14469 = t2066-t2067-t13177+t13179+t13181-t13183-t13185+t13187-t13
     #189-t2088+t2130-t2131+t13190-t13191+t13192-t13193+t13194+t13196-t1
     #3198
      t14470 = t13201-t13203-t13204-t2154-t2163+t2170-t2173-t13205+t1320
     #6-t13207+t13208+t2210+t2212+t2214+t2215+t2216-t2217+t13209-t2220
      A1e0c7p1t1n5 = 12*(t14469+t14470)*t14423*t2226*t14429*t13213*t149
      A1e0c8p1t1n5 = 0
      t14477 = t13218+t13219-t13221-t13223+t13225-t13227+t13229-t13231-t
     #13233+t13235-t13237
      t14478 = t13241-t13243+t13245+t13248-t13251-t13253+t13255-t13258+t
     #13261+t13264-t13266
      t14480 = t13271-t13273+t13275+t13278-t13280+t13283-t13285+t13288+t
     #13290+t13292+t2751
      t14483 = t13315-t3100+t13316-t13317+t13319-t13321-t13323-t13325-t1
     #3327-t13329-t13331-t13333
      t14487 = t13339-t13341+t13343-t13345+t13347+t13349-t13352+t13354-t
     #13356-t13358+t13360
      t14488 = t13363-t13365-t3360-t13366+t13367+t13368-t13369+t3375+t33
     #99+t13370+t3403
      t14490 = t13373-t13374-t13375-t3411+t13376-t13377-t13378-t13380+t1
     #3382+t3094-t13383
      t14491 = t13386-t13389+t13391+t13393+t13395+t13397-t13399-t3506+t1
     #3400+t3523+t13401-t13402
      A1e0c9p1t1n5 = -12*(t14477+t14478+t14480+t13299+t13310+t13314+t144
     #83+t14487+t14488+t14490+t14491+t13420)*t14423*t151*t3838*t14466
      t14501 = -t3842-t13428+t13430+t13432-t13434-t13436-t13438+t13439+t
     #13441-t13443+t13445-t13447-t3925-t943+t4407-t13448+t13449+t3938
      t14502 = -t3940-t13451+t13452+t3946-t3948+t4413+t4415-t13453+t1345
     #4+t13455-t13456-t3967+t3969-t13457+t4005-t13459-t13461-t13463
      t14504 = t13467-t13469+t13471-t13473-t13475+t13477+t13479+t13480-t
     #13481+t13483+t13485-t13487-t13489-t13491-t13493-t13494+t4022+t4032
      t14510 = t14429*t571*t569
      A1e0c10p1t1n5 = 12*(t14501+t14502+t14504+t13502)*t14423*t4108*t145
     #10
      A1e0c11p1t1n5 = 0
      t14515 = -t6210+t13574-t13575+t6219+t13576+t13580+t13582+t13584-t1
     #3587-t13589+t13591
      t14520 = -t13620+t13622+t13624+t13626-t13628+t13630-t13632-t13635-
     #t13637+t13639-t13641
      t14523 = t13666-t13668-t13671+t13673+t13675+t13677-t6398+t13680+t1
     #3682-t13685+t13687
      t14532 = -t13744-t13746+t13748+t6660-t13750+t13752-t13754+t13756+t
     #13759+t13761+t13763
      t14535 = t13787-t13789-t13792+t13795+t13797+t13799-t13801-t13803-t
     #13805+t13808-t13811
      t14539 = t13844-t13846-t13849-t13852+t13855+t13858+t13861-t13864-t
     #13866-t13868+t13871
      t14542 = -t13900-t13902-t13904+t13906+t13908-t13910-t13912-t13914+
     #t13917-t13919+t13921
      t14548 = t13961-t13962-t13963-t13964+t13965-t13968+t13971-t13974+t
     #13977+t13981+t13984
      t14560 = t2052*t569*t2054
      t14565 = t580*t14437*t151*t12951*t14429
      s2 = 12
      s5 = t6548-t13604-t13606-t13598+t13600+t13602-t13836-t13839-t13841
     #-t13830+t13834+t13573+t13817-t13820+t13822-t13825+t13827+t13569+t1
     #3554-t13526+t13533-t13772-t13775+t13778+t13780-t13523+t13782-t1378
     #4-t13513-t13515+t13517+t13520-t13767-t13769-t13741-t13739+t14026-t
     #14021-t14023+t14015+t14017+t14019+t14013+t13960+t13950-t13935-t139
     #37+t13939-t13941+t13942-t13943-t13933+t13927-t13929+t13931+t13897+
     #t13895+t13889+t13892
      s4 = -t13886-t13877+t13879-t13881+t13883+t13875+t14548+t6698+t1454
     #2+t14539+t14532+t14535+t14515+t14520+t14523-t5389+t4743-t4736-t137
     #36+t13733-t13730+t13726-t13727+t13722+t13709+t13718+t13705+t13707+
     #t13701-t13703+t6539-t13696+t13698-t6537+t13694+t13692+t13661-t1366
     #3+t7700+t7697+t7679-t7694-t7695-t13659-t13655-t13657-t13651+t13653
     #-t13649-t13645+t13647+t6434+t13609+t13611-t13613+t13615+t13617+t48
     #80+s5+t6683
      s5 = t14423
      s3 = s4*s5
      s1 = s2*s3
      s2 = t594*t14560*t14565
      A1e0c12p1t1n5 = s1*s2
      A1e0c1p1t2n5 = A1e0c1p1t1n5
      A1e0c2p1t2n5 = 0
      A1e0c3p1t2n5 = A1e0c3p1t1n5
      A1e0c4p1t2n5 = A1e0c4p1t1n5
      A1e0c5p1t2n5 = 0
      A1e0c6p1t2n5 = A1e0c6p1t1n5
      A1e0c7p1t2n5 = A1e0c7p1t1n5
      A1e0c8p1t2n5 = 0
      A1e0c9p1t2n5 = A1e0c9p1t1n5
      A1e0c10p1t2n5 = A1e0c10p1t1n5
      A1e0c11p1t2n5 = 0
      A1e0c12p1t2n5 = A1e0c12p1t1n5
      t14567 = -t14043+t14042+t7840+t12910-t32+t34+t12917-t7837-t12921+t
     #142-t7838-t14041
      A2e0c1p1t1n5 = 48*t14567*Mzsq*t14445*t14448
      A2e0c2p1t1n5 = 0
      t14571 = -t14051+t14054-t14056+t8054-t14057-t8088-t8091+t8093+t810
     #2-t8104-t8106-t8109-t8111+t8118-t8122+t8140+t8144+t8152-t8172
      t14574 = t8477-t14063+t8522-t8523-t5953+t8526-t8529+t14064+t6057+t
     #8532-t8533-t6069-t14067-t14071-t14075+t14078+t14080+t14082-t14084
      t14575 = -t14088-t14091+t14094-t14096+t14098+t14102-t14106-t14109+
     #t14112-t14115+t14117+t14120-t14122-t14124-t14128-t14130-t14132+t14
     #134+t14135+t14138
      t14578 = -t14144+t14147-t14150+t14152-t14154-t14156+t14158-t14160-
     #t14162+t14164+t14165-t14167+t14170-t14171+t14172-t13659+t13666-t85
     #86+t14173
      t14579 = -t8590-t8597-t8598+t8599-t8600+t8601+t14176+t14178+t14181
     #-t14183+t14186-t14188-t14191-t14193-t14195+t14197-t14199-t14202-t1
     #4204+t14206
      t14581 = t14210-t14211+t14213+t14215+t14217-t14219-t14221+t14223+t
     #8665+t14224-t14225-t14226-t13801+t14227+t14228-t8724+t14229+t14230
     #-t8906
      t14582 = -t8938-t8965+t9039+t9049-t9058+t9060-t9061+t9062+t9064-t9
     #065-t9069+t9071-t9073+t9088-t9090-t9092+t14233+t14235-t14238+t9095
      A2e0c3p1t1n5 = 24*(t14571-t14061+t14574+t14575+t14578+t14579+t1458
     #1+t14582)*t14423*t2052*t9104*t14437*t151*t149*t9106*t12782*t12201
      t14595 = t14252-t14254-t14256+t13432+t13436-t13439+t13443+t9134-t9
     #43+t9136+t9137-t3940+t3946-t9138-t9139-t9140+t3969-t9143
      t14596 = t4005-t14258-t13461+t14259-t14260+t14261-t14262-t13475+t1
     #4263+t13479+t13483-t13487-t1093-t4096-t9151+t4097+t4099+t9153
      A2e0c4p1t1n5 = 12*(t14595+t14596)*t14423*t4108*t14510
      A2e0c5p1t1n5 = 0
      t14602 = t14276+t14278-t14280+t14282-t14284+t14286-t14288-t14289-t
     #14291+t14294+t14296
      t14609 = -t14344+t14345-t6088+t14346+t6106-t6108+t14347-t14348-t14
     #349+t14350-t13575
      t14612 = -t14359+t14360-t13649-t14361+t14362-t14363-t13675-t14364+
     #t13680-t13685+t13687
      t14619 = t13739-t13744-t13746+t13748+t13754-t13756+t14392+t6683-t1
     #3767-t13769-t13772
      t14622 = t13844+t13846-t13855+t14400+t13902+t13912-t13927+t13929-t
     #13931+t13935+t6871
      t14625 = -t7703-t7705-t7711+t7713-t7714-t7723+t14408-t14410-t14411
     #+t14412+t14413
      A2e0c6p1t1n5 = 12*(t13591-t13587+t13580-t6928+t6910-t6876-t13526-t
     #13523-t13513+t13517+t13520-t14015+t14609+t14612+t14619+t14622+t146
     #25-t14356-t14357-t13948+t14602+t14353-t14354+t14355-t14343-t14335+
     #t14323+t14324-t14318+t14322-t14313+t14316-t14308-t14311-t14301+t14
     #305-t14271-t14272-t14273+t7384+t4743+t13701-t6505+t6503-t6501-t770
     #1+t7700-t7694+t6496+t13615+t13622+t4880+t6489+t14403+t14399+t14391
     #-t14379-t14371-t14372-t14369-t14370)*t14423*t594*t14560*t14565
      A2e0c1p1t2n5 = A2e0c1p1t1n5
      A2e0c2p1t2n5 = 0
      A2e0c3p1t2n5 = A2e0c3p1t1n5
      A2e0c4p1t2n5 = A2e0c4p1t1n5
      A2e0c5p1t2n5 = 0
      A2e0c6p1t2n5 = A2e0c6p1t1n5

      RETURN
      END
