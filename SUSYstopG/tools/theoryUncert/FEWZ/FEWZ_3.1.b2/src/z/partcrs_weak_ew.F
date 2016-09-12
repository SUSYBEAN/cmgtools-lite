
      SUBROUTINE InitConstw()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      DOUBLE COMPLEX spf10
      DOUBLE COMPLEX spf11
      DOUBLE COMPLEX spf12
      DOUBLE COMPLEX spf13
      DOUBLE COMPLEX spf14
      DOUBLE COMPLEX spf15
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
      DOUBLE COMPLEX spf38
      DOUBLE COMPLEX spf39
      DOUBLE COMPLEX spf40
      DOUBLE COMPLEX spf41
      DOUBLE COMPLEX spf8
      DOUBLE COMPLEX spf9

C      t1 = mur**2
	if (murflag.gt.0d0) then 
		t1 = dble(Mz**2)
	else
		t1 = mur**2
	endif
      t2 = 1/Mzsq
      t3 = t1*t2
      t4 = 1/cw2
      t6 = log(DCMPLX(t3*t4))
      t7 = cw2*t6
      t8 = log(DCMPLX(t3))
      t10 = cw2-1
      t11 = 1/t10
      wC1 = -(-t7+t8-cw2+1)*t11/6
      t14 = Mzsq*cw2
      t15 = t14*t6
      t16 = Mh**2
      t17 = 1/t16
      t19 = log(DCMPLX(t1*t17))
      t20 = t16*t19
      t23 = 1/(-t14+t16)
      wC2 = (-t15+t20-t14+t16)*t23/6
      wC3 = (-Mzsq*t8+t20-Mzsq+t16)/(-Mzsq+t16)/6
      t30 = Mt**2
      t33 = log(DCMPLX(t1/t30))
      wC4 = (t30*t33-t15+t30-t14)/(t30-t14)/6
      spf15 = log(DCMPLX(-t1/(-t14+Ieps)))
      t42 = t2*t4
      t43 = t42*Ieps
c      spf18 = log(DCMPLX(-1-t43))
	spf18 = log(-1-Ieps)
      wC5 = 1-spf18/2+spf15/2
      spf13 = log(DCMPLX(-t43+t4))
      t48 = 4*t14
      t53 = sqrt(Mzsq*(-t48+Mzsq+4*cw2*Ieps))
      t54 = Mzsq-t53
      spf20 = log(DCMPLX(-t54*t2*t4))
      t57 = Mzsq+t53
      spf21 = log(DCMPLX(-t57*t2*t4))
      t60 = 2*t14
      t61 = -Mzsq-t53+t60
      spf24 = log(DCMPLX(-t61*t2*t4))
      t64 = -Mzsq+t53+t60
      spf25 = log(DCMPLX(-t64*t2*t4))
      spf32 = log(DCMPLX(-t61/t57))
      spf33 = log(DCMPLX(-t64/t54))
      t78 = 2*spf15*Mzsq*cw2
      t89 = -2*t14*spf24+t48+spf32*t53+spf33*Mzsq-spf33*t53+t78+spf32*Mz
     #sq-2*t14*spf25-2*spf13*Mzsq*cw2+2*t14*spf21+2*t14*spf20
      wC6 = t89*t2*t4/48
      t92 = t16*t2
      spf14 = log(DCMPLX(-t43+t92*t4))
      t95 = -t48+t16
      t96 = sqrt(t95)
      t97 = t16*Mzsq
      t100 = t16**2
      t107 = sqrt((-4*t97*cw2+t100+4*t14*Ieps)*t17/t95)
      t108 = t96*t107
      t109 = Mh-t108
      spf34 = log(DCMPLX(-Mh*t109*t42))
      t112 = Mh+t108
      spf35 = log(DCMPLX(-Mh*t112*t42))
      t116 = t96*Mh*t107
      t117 = t16-t116-t60
      spf36 = log(DCMPLX(t117*t2*t4))
      t120 = t16+t116-t60
      spf37 = log(DCMPLX(t120*t2*t4))
      t123 = 1/Mh
      spf40 = log(DCMPLX(t117*t123/t109))
      spf41 = log(DCMPLX(t120*t123/t112))
      t147 = t16*spf40+t48-2*spf14*Mzsq*cw2+t78+2*t14*spf35+2*t14*spf34-
     #2*t14*spf37+Mh*spf41*t108-Mh*spf40*t108+t16*spf41-2*t14*spf36
      wC7 = t147*t2*t4/48
      spf19 = log(DCMPLX((-t30+t14+Ieps)/(-t30+Ieps)))
      t157 = t14*t33
      wC8 = (t60+spf19*t30-spf19*Mzsq*cw2+t157)*t2*t4/12
      spf8 = log(DCMPLX(-t1/(-Mzsq+Ieps)))
      t164 = t2*Ieps
c      spf9 = log(DCMPLX(-1-t164))
	spf9 = log(-1-Ieps)
      wC9 = 1-spf9/2+spf8/2
      spf16 = log(DCMPLX(-t164+cw2))
      t169 = 4*Ieps
      t172 = sqrt((Mzsq-t48+t169)*t2)
      t173 = -1+t172
      t174 = 1+t172
      spf26 = log(DCMPLX(t173/t174))
      spf27 = log(DCMPLX(t174/t173))
      t181 = spf8/4
      wC10 = 1.D0/2.D0-spf16/4+spf26/8+t181+spf26*t172/8+spf27/8-spf27*t
     #172/8
      spf10 = log(DCMPLX(-t164+t92))
      t188 = 4*Mzsq
      t189 = -t188+t16
      t190 = sqrt(t189)
      t198 = sqrt((-4*t97+t100+4*Mzsq*Ieps)*t17/t189)
      t199 = t190*t198
      t200 = Mh-t199
      spf28 = log(DCMPLX(-Mh*t200*t2))
      t203 = Mh+t199
      spf29 = log(DCMPLX(-Mh*t203*t2))
      t207 = t190*Mh*t198
      t208 = 2*Mzsq
      t209 = t16-t207-t208
      spf30 = log(DCMPLX(t209*t2))
      t211 = t16+t207-t208
      spf31 = log(DCMPLX(t211*t2))
      spf38 = log(DCMPLX(t209*t123/t200))
      spf39 = log(DCMPLX(t211*t123/t203))
      t237 = -2*Mzsq*spf31+t188+2*Mzsq*spf29+2*Mzsq*spf28+Mh*spf39*t199-
     #2*Mzsq*spf30+t16*spf38+2*spf8*Mzsq-Mh*spf38*t199+t16*spf39-2*spf10
     #*Mzsq
      wC11 = -t237*t2/48
      t240 = t30*t2
      spf11 = log(DCMPLX(-t164+t240))
      t245 = sqrt((Mzsq-4*t30+t169)*t2)
      t246 = -1+t245
      t247 = 1+t245
      spf22 = log(DCMPLX(t246/t247))
      spf23 = log(DCMPLX(t247/t246))
      wC12 = 1.D0/2.D0-spf23*t245/8+spf22*t245/8+t181+spf23/8-spf11/4+sp
     #f22/8
      t259 = dble(Mzsq)
      spf12 = log(DCMPLX(-t1/(-t259+Ieps)))
      spf17 = log(DCMPLX(-1-1/t259*Ieps))
      wC13 = 1-spf17/2+spf12/2
      t268 = t30*Mzsq
      t270 = t268*wC12*cw2
      t272 = cw2**2
      t273 = wC12*t272
      t274 = t268*t273
      t276 = t33*cw2
      t277 = t268*t276
      t279 = t33*t272
      t280 = t268*t279
      t288 = Mzsq**2
      t289 = t288*wC12
      t290 = t289*cw2
      t292 = t288*t272
      t293 = t292*wC10
      t296 = t288*cw2
      t302 = t288*wC10
      t306 = 320*t270-256*t274-80*t277+64*t280-28*t268*wC12-3*t97*t19-18
     #*t97*wC3+160*t290+408*t293+64*t292-76*t296+36*t100*wC11+9*t100*wC3
     #-6*t302+432*t288*wC11
      t307 = t288*wC9
      t314 = t307*t272
      t318 = t296*wC10
      t320 = t272*cw2
      t321 = t288*t320
      t322 = t321*wC10
      t324 = t289*t272
      t326 = t307*cw2
      t332 = t292*t6
      t334 = t321*t6
      t337 = -206*t307-68*t289+9*t288*wC3-3*t288*t8-320*t314-144*t97*wC1
     #1-96*t318+288*t322-128*t324+400*t326+34*t268*t33-6*t296*t6+24*t332
     #-72*t334+38*t288
      wC14 = (t306+t337)*t11*t42/72
      t344 = dble(wC13)
      t348 = t16*t30
      t349 = t272**2
      t350 = Mzsq*t349
      t354 = Mzsq*t320
      t358 = Mzsq*t272
      t379 = t16*t288
      t383 = t349*cw2
      t387 = t100*Mzsq
      t388 = t272*wC3
      t391 = t272*wC11
      t394 = t320*wC3
      t397 = t320*wC11
      t400 = t30**2
      t401 = t400*Mzsq
      t402 = t272*wC8
      t408 = -1024*t348*t350*wC12+1280*t348*t354*wC12-112*t348*t358*wC12
     #-864*t348*t358*wC8+432*t348*t14*wC8+256*t348*t350*t33-320*t348*t35
     #4*t33+100*t348*t358*t33-36*t348*t157+1632*t379*t349*wC10+1152*t379
     #*t383*wC10-72*t387*t388-576*t387*t391-36*t387*t394-144*t387*t397+8
     #64*t401*t402-12672*t379*t320*wC6
      t409 = t272*wC7
      t412 = t272*wC2
      t421 = t30*t288
      t422 = t320*wC8
      t439 = t320*wC12
      t442 = t349*wC12
      t460 = -1440*t387*t409+288*t387*t412+2736*t379*t320*wC1-13824*t379
     #*t349*wC6+864*t421*t422+1728*t379*t422+4608*t379*t320*wC7+864*t379
     #*t320*wC5-648*t379*t320*wC2-384*t379*t320*wC10+112*t421*t439-1280*
     #t421*t442-824*t379*t272*wC9+36*t379*t388-272*t379*t273+1728*t379*t
     #391-24*t379*t272*wC10-1280*t379*t349*wC9
      t469 = t30*t272
      t493 = t16*t400
      t494 = cw2*wC8
      t507 = -512*t379*t442+1024*t421*t383*wC12-54*t348*t14+54*t97*t469-
     #2016*t379*cw2*wC6+504*t379*cw2*wC1+14400*t379*t272*wC6-3276*t379*t
     #272*wC1+1600*t379*t320*wC9+72*t379*t394+640*t379*t439+576*t379*t39
     #7-864*t493*t494-432*t421*t402+720*t387*cw2*wC7-180*t387*cw2*wC2-86
     #4*t379*t402
      t526 = log(DCMPLX(cw2))
      t543 = t320*t33
      t554 = log(DCMPLX(t14*t17))
      t558 = -2304*t379*t409-432*t379*t272*wC5+540*t379*t412-432*t401*t4
     #94+576*t379*t349*spf15-864*t379*t320*spf15+288*t379*t272*spf15+114
     #*t379*t272*t526+192*t379*t320*t526-256*t421*t383*t33+320*t421*t349
     #*t33-72*t401*t279+72*t493*t276-100*t421*t543+36*t401*t276+36*t421*
     #t279+96*t379*t320*t8-6*t379*t320*t554
      t578 = t288*Mzsq
      t579 = t578*t383
      t582 = t578*t349
      t585 = t578*t320
      t589 = t578*t272
      t594 = t100*t16
      t601 = 6*t379*t272*t554-144*t379*t320*t6+96*t379*t349*t6-288*t379*
     #t383*t6+240*t379*t272*t6+3*t379*cw2+13824*t579*wC6-192*t582*t526+8
     #24*t585*wC9-1267*t585-3*t589+2678*t582-1408*t579-36*t493-144*wC7*t
     #594+36*wC2*t594-3456*t582*wC7
      t602 = t594*t272
      t616 = t578*t349*t272
      t635 = t469*t288
      t637 = t401*cw2
      t639 = t578*cw2
      t642 = 36*t602*wC3-240*t585*t6+272*t585*wC12+1728*t585*wC7+2016*t5
     #89*wC6+1280*t579*wC9-1152*t616*wC10-2736*t582*wC1-114*t585*t526-16
     #00*t582*wC9+144*t602*wC11-504*t589*wC1+512*t579*wC12-2675*t379*t32
     #0+1264*t379*t272+54*t635+36*t637+144*t639*wC6
      t662 = t594*cw2
      t669 = t401*t272
      t679 = -36*t639*wC1+864*t585*wC8+432*t585*wC5-396*t585*wC2+432*t49
     #3*wC8+3*t387*cw2-144*t379*wC6+36*t379*wC1+72*t493*cw2+288*t662*wC7
     #-72*t662*wC2+1408*t379*t349-72*t669-36*t585*wC3-1728*t585*wC11+24*
     #t585*wC10+12672*t582*wC6
      t697 = t320*t30*t288
      t703 = t582*t8
      t717 = -640*t582*wC12+384*t582*wC10-1632*t579*wC10-14400*t585*wC6+
     #3276*t585*wC1-1728*t582*wC8-864*t582*wC5+432*t582*wC2-54*t697-3*t3
     #87*t272-36*t493*t33-96*t703-288*t585*spf15-576*t579*spf15+864*t582
     #*spf15-96*t579*t6+288*t616*t6+144*t582*t6
      t722 = t10**2
      t723 = 1/t722
      t724 = 1/t288
      if (whichscheme .eq. 2) then
      wC16 = -7.D0/8.D0*t6+2.D0/9.D0*t33-49.D0/108.D0+20.D0/9.D0*t344
      else if (whichscheme .eq. 1) then
      wC16 = -(t408+t460+t507+t558+t601+t642+t679+t717)*t23*t723*t724/t2
     #72/576
      else
      wC16 = 0
      end if
      t755 = -28*t270+320*t274+16*t277-80*t280-144*t14*t16*wC7+18*t14*t1
     #6*wC2-108*t14*t30*wC8-144*t14*t16*wC11-18*t14*t16*wC3-256*t268*t43
     #9+64*t268*t543-68*t290
      t760 = t288*wC6
      t775 = -96*t293-258*t292+38*t296+220*t321+36*t760-9*t288*wC1+36*t1
     #00*wC7-9*t100*wC2-108*t400*wC8+9*t400*t33+400*t314-6*t318+408*t322
      t799 = 160*t324-206*t326+12*t332+36*t334+576*t296*wC6-54*t296*wC1+
     #432*t296*wC11+9*t296*wC3-1728*t760*t320-2448*t760*t272+288*t302*t3
     #49-72*t321*wC1+135*t292*wC1
      t812 = cw2*t100
      t827 = 432*t292*wC7-9*t292*wC2+216*t292*wC8+108*t292*wC5-128*t289*
     #t320-320*t307*t320+36*t812*wC11+9*t812*wC3+24*t292*t8+72*t321*spf1
     #5-72*t292*spf15-72*t288*t349*t6+9*t400
      wC18 = -t724*t4*(t755+t775+t799+t827)*t723/144
      t838 = t400*t30
      t849 = -58*t585+2*t421*cw2+71*t635-t589-22*t637+9*t838+8*t585*t8-t
     #401+108*t838*wC4+2*t589*t8+216*t585*wC4-4*t582
      t850 = t272*t8
      t853 = t30*t8
      t858 = t8*cw2
      t865 = log(DCMPLX(t240*t4))
      t879 = 8*t401*t850-16*t321*t853-16*t292*t853+8*t401*t858-324*t292*
     #t30*wC4-18*t401*cw2*t865-36*t421*t272*t865-4*t296*t853+2*t401*t8-4
     #*t669+8*t703+8*t697
      wC20 = t4*(t849+t879)*t2/(-t292+t321+2*t268*cw2-2*t268*t272-t400+t
     #400*cw2)/288
      t891 = 2*t8
      t892 = 36*t7
      t894 = 8*t858
      t895 = 8*t850
      t896 = 4*t272
      wC22 = (t891+t892-22*cw2+t894+t895-1-t896)*t4*t11/288
      wC24 = (2*cw2+t895-t894+t891-t896-1+4*t7)*t4*t11/32
      wC26 = (-10*cw2-1+t892+32*t850-16*t858+t891-16*t272)*t4*t11/288
      wC28 = (-t891+2*t858+1-cw2)*t4/72

      RETURN
      END

      SUBROUTINE CrsecSeqw()

      IMPLICIT DOUBLE COMPLEX (t, s)

      include 'constcommon_ew.f'

      include 'eventcommon_ew.f'

      include 'crseccommon_ew.f'

      DOUBLE COMPLEX spf101
      DOUBLE COMPLEX spf106
      DOUBLE COMPLEX spf107
      DOUBLE COMPLEX spf108
      DOUBLE COMPLEX spf109
      DOUBLE COMPLEX spf112
      DOUBLE COMPLEX spf113
      DOUBLE COMPLEX spf114
      DOUBLE COMPLEX spf117
      DOUBLE COMPLEX spf118
      DOUBLE COMPLEX spf119
      DOUBLE COMPLEX spf120
      DOUBLE COMPLEX spf121
      DOUBLE COMPLEX spf122
      DOUBLE COMPLEX spf123
      DOUBLE COMPLEX spf124
      DOUBLE COMPLEX spf125
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
      DOUBLE COMPLEX spf150
      DOUBLE COMPLEX spf151
      DOUBLE COMPLEX spf152
      DOUBLE COMPLEX spf153
      DOUBLE COMPLEX spf156
      DOUBLE COMPLEX spf157
      DOUBLE COMPLEX spf158
      DOUBLE COMPLEX spf159
      DOUBLE COMPLEX spf162
      DOUBLE COMPLEX spf163
      DOUBLE COMPLEX spf166
      DOUBLE COMPLEX spf167
      DOUBLE COMPLEX spf168
      DOUBLE COMPLEX spf169
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
      DOUBLE COMPLEX spf180
      DOUBLE COMPLEX spf181
      DOUBLE COMPLEX spf182
      DOUBLE COMPLEX spf183
      DOUBLE COMPLEX spf184
      DOUBLE COMPLEX spf185
      DOUBLE COMPLEX spf186
      DOUBLE COMPLEX spf187
      DOUBLE COMPLEX spf188
      DOUBLE COMPLEX spf189
      DOUBLE COMPLEX spf190
      DOUBLE COMPLEX spf191
      DOUBLE COMPLEX spf194
      DOUBLE COMPLEX spf195
      DOUBLE COMPLEX spf196
      DOUBLE COMPLEX spf197
      DOUBLE COMPLEX spf198
      DOUBLE COMPLEX spf20
      DOUBLE COMPLEX spf201
      DOUBLE COMPLEX spf202
      DOUBLE COMPLEX spf203
      DOUBLE COMPLEX spf204
      DOUBLE COMPLEX spf205
      DOUBLE COMPLEX spf206
      DOUBLE COMPLEX spf207
      DOUBLE COMPLEX spf208
      DOUBLE COMPLEX spf21
      DOUBLE COMPLEX spf210
      DOUBLE COMPLEX spf211
      DOUBLE COMPLEX spf212
      DOUBLE COMPLEX spf214
      DOUBLE COMPLEX spf215
      DOUBLE COMPLEX spf216
      DOUBLE COMPLEX spf217
      DOUBLE COMPLEX spf218
      DOUBLE COMPLEX spf219
      DOUBLE COMPLEX spf22
      DOUBLE COMPLEX spf220
      DOUBLE COMPLEX spf221
      DOUBLE COMPLEX spf222
      DOUBLE COMPLEX spf223
      DOUBLE COMPLEX spf224
      DOUBLE COMPLEX spf225
      DOUBLE COMPLEX spf226
      DOUBLE COMPLEX spf227
      DOUBLE COMPLEX spf228
      DOUBLE COMPLEX spf229
      DOUBLE COMPLEX spf230
      DOUBLE COMPLEX spf231
      DOUBLE COMPLEX spf232
      DOUBLE COMPLEX spf233
      DOUBLE COMPLEX spf234
      DOUBLE COMPLEX spf235
      DOUBLE COMPLEX spf236
      DOUBLE COMPLEX spf237
      DOUBLE COMPLEX spf238
      DOUBLE COMPLEX spf239
      DOUBLE COMPLEX spf24
      DOUBLE COMPLEX spf240
      DOUBLE COMPLEX spf241
      DOUBLE COMPLEX spf242
      DOUBLE COMPLEX spf243
      DOUBLE COMPLEX spf244
      DOUBLE COMPLEX spf245
      DOUBLE COMPLEX spf246
      DOUBLE COMPLEX spf247
      DOUBLE COMPLEX spf248
      DOUBLE COMPLEX spf249
      DOUBLE COMPLEX spf250
      DOUBLE COMPLEX spf251
      DOUBLE COMPLEX spf252
      DOUBLE COMPLEX spf253
      DOUBLE COMPLEX spf254
      DOUBLE COMPLEX spf255
      DOUBLE COMPLEX spf256
      DOUBLE COMPLEX spf257
      DOUBLE COMPLEX spf258
      DOUBLE COMPLEX spf259
      DOUBLE COMPLEX spf26
      DOUBLE COMPLEX spf260
      DOUBLE COMPLEX spf261
      DOUBLE COMPLEX spf262
      DOUBLE COMPLEX spf263
      DOUBLE COMPLEX spf264
      DOUBLE COMPLEX spf265
      DOUBLE COMPLEX spf266
      DOUBLE COMPLEX spf267
      DOUBLE COMPLEX spf268
      DOUBLE COMPLEX spf269
      DOUBLE COMPLEX spf27
      DOUBLE COMPLEX spf270
      DOUBLE COMPLEX spf271
      DOUBLE COMPLEX spf272
      DOUBLE COMPLEX spf273
      DOUBLE COMPLEX spf274
      DOUBLE COMPLEX spf275
      DOUBLE COMPLEX spf276
      DOUBLE COMPLEX spf277
      DOUBLE COMPLEX spf278
      DOUBLE COMPLEX spf279
      DOUBLE COMPLEX spf28
      DOUBLE COMPLEX spf280
      DOUBLE COMPLEX spf281
      DOUBLE COMPLEX spf282
      DOUBLE COMPLEX spf283
      DOUBLE COMPLEX spf284
      DOUBLE COMPLEX spf285
      DOUBLE COMPLEX spf286
      DOUBLE COMPLEX spf287
      DOUBLE COMPLEX spf288
      DOUBLE COMPLEX spf289
      DOUBLE COMPLEX spf29
      DOUBLE COMPLEX spf290
      DOUBLE COMPLEX spf291
      DOUBLE COMPLEX spf292
      DOUBLE COMPLEX spf293
      DOUBLE COMPLEX spf294
      DOUBLE COMPLEX spf295
      DOUBLE COMPLEX spf296
      DOUBLE COMPLEX spf297
      DOUBLE COMPLEX spf298
      DOUBLE COMPLEX spf299
      DOUBLE COMPLEX spf30
      DOUBLE COMPLEX spf300
      DOUBLE COMPLEX spf301
      DOUBLE COMPLEX spf302
      DOUBLE COMPLEX spf303
      DOUBLE COMPLEX spf304
      DOUBLE COMPLEX spf305
      DOUBLE COMPLEX spf306
      DOUBLE COMPLEX spf307
      DOUBLE COMPLEX spf308
      DOUBLE COMPLEX spf309
      DOUBLE COMPLEX spf310
      DOUBLE COMPLEX spf311
      DOUBLE COMPLEX spf312
      DOUBLE COMPLEX spf313
      DOUBLE COMPLEX spf314
      DOUBLE COMPLEX spf315
      DOUBLE COMPLEX spf316
      DOUBLE COMPLEX spf317
      DOUBLE COMPLEX spf318
      DOUBLE COMPLEX spf319
      DOUBLE COMPLEX spf32
      DOUBLE COMPLEX spf320
      DOUBLE COMPLEX spf321
      DOUBLE COMPLEX spf322
      DOUBLE COMPLEX spf323
      DOUBLE COMPLEX spf324
      DOUBLE COMPLEX spf325
      DOUBLE COMPLEX spf326
      DOUBLE COMPLEX spf327
      DOUBLE COMPLEX spf328
      DOUBLE COMPLEX spf329
      DOUBLE COMPLEX spf330
      DOUBLE COMPLEX spf331
      DOUBLE COMPLEX spf332
      DOUBLE COMPLEX spf333
      DOUBLE COMPLEX spf334
      DOUBLE COMPLEX spf335
      DOUBLE COMPLEX spf336
      DOUBLE COMPLEX spf337
      DOUBLE COMPLEX spf338
      DOUBLE COMPLEX spf339
      DOUBLE COMPLEX spf34
      DOUBLE COMPLEX spf340
      DOUBLE COMPLEX spf341
      DOUBLE COMPLEX spf342
      DOUBLE COMPLEX spf343
      DOUBLE COMPLEX spf344
      DOUBLE COMPLEX spf345
      DOUBLE COMPLEX spf346
      DOUBLE COMPLEX spf347
      DOUBLE COMPLEX spf348
      DOUBLE COMPLEX spf349
      DOUBLE COMPLEX spf35
      DOUBLE COMPLEX spf350
      DOUBLE COMPLEX spf351
      DOUBLE COMPLEX spf352
      DOUBLE COMPLEX spf353
      DOUBLE COMPLEX spf354
      DOUBLE COMPLEX spf355
      DOUBLE COMPLEX spf36
      DOUBLE COMPLEX spf37
      DOUBLE COMPLEX spf38
      DOUBLE COMPLEX spf43
      DOUBLE COMPLEX spf44
      DOUBLE COMPLEX spf45
      DOUBLE COMPLEX spf46
      DOUBLE COMPLEX spf50
      DOUBLE COMPLEX spf51
      DOUBLE COMPLEX spf57
      DOUBLE COMPLEX spf58
      DOUBLE COMPLEX spf59
      DOUBLE COMPLEX spf60
      DOUBLE COMPLEX spf61
      DOUBLE COMPLEX spf62
      DOUBLE COMPLEX spf63
      DOUBLE COMPLEX spf64
      DOUBLE COMPLEX spf65
      DOUBLE COMPLEX spf67
      DOUBLE COMPLEX spf68
      DOUBLE COMPLEX spf69
      DOUBLE COMPLEX spf70
      DOUBLE COMPLEX spf71
      DOUBLE COMPLEX spf72
      DOUBLE COMPLEX spf73
      DOUBLE COMPLEX spf78
      DOUBLE COMPLEX spf79
      DOUBLE COMPLEX spf80
      DOUBLE COMPLEX spf81
      DOUBLE COMPLEX spf83
      DOUBLE COMPLEX spf84
      DOUBLE COMPLEX spf85
      DOUBLE COMPLEX spf86
      DOUBLE COMPLEX spf87
      DOUBLE COMPLEX spf88
      DOUBLE COMPLEX spf89
      DOUBLE COMPLEX spf90
      DOUBLE COMPLEX spf93
      DOUBLE COMPLEX spf94
      DOUBLE COMPLEX spf96
      DOUBLE COMPLEX spf97
      DOUBLE COMPLEX spf98
      DOUBLE COMPLEX wI100pc1
      DOUBLE COMPLEX wI101pc1
      DOUBLE COMPLEX wI102pc1
      DOUBLE COMPLEX wI103pc1
      DOUBLE COMPLEX wI104pc1
      DOUBLE COMPLEX wI105pc1
      DOUBLE COMPLEX wI106pc1
      DOUBLE COMPLEX wI107pc1
      DOUBLE COMPLEX wI108pc1
      DOUBLE COMPLEX wI109pc1
      DOUBLE COMPLEX wI110pc1
      DOUBLE COMPLEX wI111pc1
      DOUBLE COMPLEX wI112pc1
      DOUBLE COMPLEX wI113pc1
      DOUBLE COMPLEX wI114pc1
      DOUBLE COMPLEX wI115pc1
      DOUBLE COMPLEX wI116pc1
      DOUBLE COMPLEX wI117pc1
      DOUBLE COMPLEX wI118pc1
      DOUBLE COMPLEX wI119pc1
      DOUBLE COMPLEX wI120pc1
      DOUBLE COMPLEX wI121pc1
      DOUBLE COMPLEX wI122pc1
      DOUBLE COMPLEX wI123pc1
      DOUBLE COMPLEX wI124pc1
      DOUBLE COMPLEX wI125pc1
      DOUBLE COMPLEX wI126pc1
      DOUBLE COMPLEX wI127pc1
      DOUBLE COMPLEX wI128pc1
      DOUBLE COMPLEX wI129pc1
      DOUBLE COMPLEX wI130pc1
      DOUBLE COMPLEX wI131pc1
      DOUBLE COMPLEX wI132pc1
      DOUBLE COMPLEX wI133pc1
      DOUBLE COMPLEX wI134pc1
      DOUBLE COMPLEX wI135pc1
      DOUBLE COMPLEX wI136pc1
      DOUBLE COMPLEX wI137pc1
      DOUBLE COMPLEX wI138pc1
      DOUBLE COMPLEX wI139pc1
      DOUBLE COMPLEX wI140pc1
      DOUBLE COMPLEX wI141pc1
      DOUBLE COMPLEX wI142pc1
      DOUBLE COMPLEX wI143pc1
      DOUBLE COMPLEX wI144pc1
      DOUBLE COMPLEX wI145pc1
      DOUBLE COMPLEX wI146pc1
      DOUBLE COMPLEX wI147pc1
      DOUBLE COMPLEX wI148pc1
      DOUBLE COMPLEX wI149pc1
      DOUBLE COMPLEX wI150pc1
      DOUBLE COMPLEX wI151pc1
      DOUBLE COMPLEX wI152pc1
      DOUBLE COMPLEX wI153pc1
      DOUBLE COMPLEX wI154pc1
      DOUBLE COMPLEX wI155pc1
      DOUBLE COMPLEX wI156pc1
      DOUBLE COMPLEX wI157pc1
      DOUBLE COMPLEX wI158pc1
      DOUBLE COMPLEX wI159pc1
      DOUBLE COMPLEX wI160pc1
      DOUBLE COMPLEX wI161pc1
      DOUBLE COMPLEX wI162pc1
      DOUBLE COMPLEX wI163pc1
      DOUBLE COMPLEX wI164pc1
      DOUBLE COMPLEX wI165pc1
      DOUBLE COMPLEX wI166pc1
      DOUBLE COMPLEX wI29pc1
      DOUBLE COMPLEX wI30pc1
      DOUBLE COMPLEX wI31pc1
      DOUBLE COMPLEX wI32pc1
      DOUBLE COMPLEX wI33pc1
      DOUBLE COMPLEX wI34pc1
      DOUBLE COMPLEX wI35pc1
      DOUBLE COMPLEX wI36pc1
      DOUBLE COMPLEX wI37pc1
      DOUBLE COMPLEX wI38pc1
      DOUBLE COMPLEX wI39pc1
      DOUBLE COMPLEX wI40pc1
      DOUBLE COMPLEX wI41pc1
      DOUBLE COMPLEX wI42pc1
      DOUBLE COMPLEX wI43pc1
      DOUBLE COMPLEX wI44pc1
      DOUBLE COMPLEX wI45pc1
      DOUBLE COMPLEX wI46pc1
      DOUBLE COMPLEX wI47pc1
      DOUBLE COMPLEX wI48pc1
      DOUBLE COMPLEX wI49pc1
      DOUBLE COMPLEX wI50pc1
      DOUBLE COMPLEX wI51pc1
      DOUBLE COMPLEX wI52pc1
      DOUBLE COMPLEX wI53pc1
      DOUBLE COMPLEX wI54pc1
      DOUBLE COMPLEX wI55pc1
      DOUBLE COMPLEX wI56pc1
      DOUBLE COMPLEX wI57pc1
      DOUBLE COMPLEX wI58pc1
      DOUBLE COMPLEX wI59pc1
      DOUBLE COMPLEX wI60pc1
      DOUBLE COMPLEX wI61pc1
      DOUBLE COMPLEX wI62pc1
      DOUBLE COMPLEX wI63pc1
      DOUBLE COMPLEX wI64pc1
      DOUBLE COMPLEX wI65pc1
      DOUBLE COMPLEX wI66pc1
      DOUBLE COMPLEX wI67pc1
      DOUBLE COMPLEX wI68pc1
      DOUBLE COMPLEX wI69pc1
      DOUBLE COMPLEX wI70pc1
      DOUBLE COMPLEX wI71pc1
      DOUBLE COMPLEX wI72pc1
      DOUBLE COMPLEX wI73pc1
      DOUBLE COMPLEX wI74pc1
      DOUBLE COMPLEX wI75pc1
      DOUBLE COMPLEX wI76pc1
      DOUBLE COMPLEX wI77pc1
      DOUBLE COMPLEX wI78pc1
      DOUBLE COMPLEX wI79pc1
      DOUBLE COMPLEX wI80pc1
      DOUBLE COMPLEX wI81pc1
      DOUBLE COMPLEX wI82pc1
      DOUBLE COMPLEX wI83pc1
      DOUBLE COMPLEX wI84pc1
      DOUBLE COMPLEX wI85pc1
      DOUBLE COMPLEX wI86pc1
      DOUBLE COMPLEX wI87pc1
      DOUBLE COMPLEX wI88pc1
      DOUBLE COMPLEX wI91pc1
      DOUBLE COMPLEX wI92pc1
      DOUBLE COMPLEX wI93pc1
      DOUBLE COMPLEX wI94pc1
      DOUBLE COMPLEX wI95pc1
      DOUBLE COMPLEX wI96pc1
      DOUBLE COMPLEX wI97pc1
      DOUBLE COMPLEX wI98pc1
      DOUBLE COMPLEX wI99pc1
      DOUBLE COMPLEX wV1
      DOUBLE COMPLEX wV10
      DOUBLE COMPLEX wV11
      DOUBLE COMPLEX wV12
      DOUBLE COMPLEX wV13
      DOUBLE COMPLEX wV14
      DOUBLE COMPLEX wV15
      DOUBLE COMPLEX wV16
      DOUBLE COMPLEX wV17
      DOUBLE COMPLEX wV18
      DOUBLE COMPLEX wV19
      DOUBLE COMPLEX wV2
      DOUBLE COMPLEX wV20
      DOUBLE COMPLEX wV25
      DOUBLE COMPLEX wV26
      DOUBLE COMPLEX wV27
      DOUBLE COMPLEX wV28
      DOUBLE COMPLEX wV3
      DOUBLE COMPLEX wV31
      DOUBLE COMPLEX wV32
      DOUBLE COMPLEX wV33
      DOUBLE COMPLEX wV34
      DOUBLE COMPLEX wV35
      DOUBLE COMPLEX wV36
      DOUBLE COMPLEX wV4
      DOUBLE COMPLEX wV5
      DOUBLE COMPLEX wV6
      DOUBLE COMPLEX wV7
      DOUBLE COMPLEX wV8
      DOUBLE COMPLEX wV9
      DOUBLE COMPLEX DENSPENCE
      DOUBLE COMPLEX IM
      PARAMETER (IM=(0.D0,1.D0))
      DOUBLE PRECISION spcZip
      PARAMETER (spcZip=0d0)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

C      t1 = mur**2
	if (murflag.gt.0d0) then 
		t1 = dble(Mz**2)
	else
		t1 = mur**2
	endif
      spf24 = log(DCMPLX(-t1/(-s+Ieps)))
      t5 = 1/s
      t6 = t5*Ieps
      t7 = Mzsq*cw2
      t8 = t7*t5
      spf38 = log(DCMPLX(-t6+t8))
      t10 = sqrt(s)
      t11 = 4*t7
      t12 = 4*Ieps
      t13 = s-t11+t12
      t14 = sqrt(t13)
      t15 = t10*t14
      t16 = -s+t15
      t17 = s+t15
      spf166 = log(DCMPLX(t16/t17))
      spf167 = log(DCMPLX(t17/t16))
      t22 = spf24/4
      t26 = sqrt(t13*t5)
      wV1 = 1.D0/2.D0+t22+spf166/8+spf167/8+spf166*t26/8-spf167*t26/8-sp
     #f38/4
      t32 = Mzsq*t5
      spf27 = log(DCMPLX(-t6+t32))
      t34 = 4*Mzsq
      t35 = s-t34+t12
      t36 = sqrt(t35)
      t37 = t10*t36
      t38 = -s+t37
      t39 = s+t37
      spf142 = log(DCMPLX(t38/t39))
      spf143 = log(DCMPLX(t39/t38))
      t46 = sqrt(t35*t5)
      wV2 = 1.D0/2.D0+t22+spf142/8+spf142*t46/8-spf143*t46/8+spf143/8-sp
     #f27/4
      t53 = Mh**2
      t54 = s**2
      t57 = s*Mzsq
      t58 = 2*t57
      t59 = t53**2
      t60 = t53*Mzsq
      t62 = Mzsq**2
      t63 = s*Ieps
      t66 = sqrt(t54-2*t53*s-t58+t59-2*t60+t62+4*t63)
      t67 = -s+t53-Mzsq-t66
      spf184 = log(DCMPLX(t67*t5))
      t69 = -s+t53-Mzsq+t66
      spf185 = log(DCMPLX(t69*t5))
      t71 = s+t53-Mzsq-t66
      spf194 = log(DCMPLX(-t71*t5))
      t73 = s+t53-Mzsq+t66
      spf195 = log(DCMPLX(-t73*t5))
      spf246 = log(DCMPLX(t71/t67))
      spf247 = log(DCMPLX(t73/t69))
      t100 = -2*s*spf184-spf247*Mzsq+spf247*t66+2*spf27*s-2*s*spf185-2*s
     #pf24*s+2*s*spf195+spf246*t53+2*s*spf194-spf246*t66-spf246*s-4*s-sp
     #f246*Mzsq-spf247*s+spf247*t53
      wV3 = -t100*t5/48
      t103 = Mt**2
      t104 = 4*t103
      t105 = s-t104+t12
      t106 = sqrt(t105)
      t107 = t10*t106
      t108 = s+t107
      t109 = t108*t5
      spf62 = log(DCMPLX(-t109))
      t110 = -s+t107
      t111 = t110*t5
      spf63 = log(DCMPLX(t111))
      spf145 = log(DCMPLX(t110/t108))
      spf146 = log(DCMPLX(t108/t110))
      t121 = sqrt(t105*t5)
      t124 = log(DCMPLX(2.D0))
      wV4 = 1.D0/2.D0-spf63/4-spf62/4+t22+spf145/8+spf146/8+spf145*t121/
     #8+t124/2-spf146*t121/8
      t128 = 2*t103
      t129 = s*CsTh
      t130 = 2*Ieps
      t131 = t128+s+t129-t130
      t132 = -t103+Ieps
      spf107 = log(DCMPLX(-t131/t132))
      t135 = 2*s
      t136 = 2*t129
      t141 = s*t124
      t142 = t141*CsTh
      t143 = s*spf107
      t145 = 1/t103
      t147 = log(DCMPLX(t1*t145))
      t148 = t147*s
      t152 = 1+CsTh
      t153 = 1/t152
      wV5 = (t135+t136+2*t103*t124-2*t103*spf107+t141+t142-t143-t143*CsT
     #h+t148+t148*CsTh)*t5*t153/12
      spf22 = log(DCMPLX(-1-t6))
      t156 = 1/t1
      spf32 = log(DCMPLX(Ieps-s*t156))
      t159 = -t7+Ieps
      t160 = 1/t159
      spf37 = log(DCMPLX(-s*t160))
      t162 = s+t7
      t164 = t7/t162
      spf50 = DENSPENCE(DCMPLX(t164),spcZip)
      t165 = s+Ieps
      spf61 = log(DCMPLX(t165*t160))
      t167 = 1/Mzsq
      t168 = 1/cw2
      t169 = t167*t168
      spf67 = log(DCMPLX(1-t169*Ieps))
      spf68 = DENSPENCE(DCMPLX(-t7*t160),spcZip)
      t173 = t162*t160
      spf94 = log(DCMPLX(-t173))
      spf96 = DENSPENCE(DCMPLX(-t173),spcZip)
      t174 = t156*Ieps
      t175 = s*cw2
      t176 = t175*Mzsq
      t177 = cw2**2
      t178 = t62*t177
      spf113 = log(DCMPLX(t174+(-t176-t178)*t5*t156))
      t184 = 1/(-t176-t178+t63)
      t185 = t178*t184
      spf120 = log(DCMPLX(-t185))
      spf121 = DENSPENCE(DCMPLX(-t185),spcZip)
      spf125 = log(DCMPLX(-t7*s*t184))
      spf139 = log(DCMPLX(s*t159*t184))
      spf156 = DENSPENCE(DCMPLX(-t7*t162*t184),spcZip)
      t192 = log(DCMPLX(t164))
      t194 = log(DCMPLX(-t7*t156))
      t200 = log(DCMPLX(-t7*t162*t5*t156))
      t211 = PI**2
      t214 = 6*spf68
      t227 = 6*t192*t194-6*t192*t200+6*spf156-6*spf121-6*spf120*spf139+6
     #*spf120*spf67+6*spf120*spf125-t211+6*spf50-6*spf96+t214-6*spf94*sp
     #f61+6*spf94*spf22+6*spf94*spf37-6*t192*spf32+6*t192*spf113-6*t192*
     #spf38
      wV6 = t227*t5/36
      t229 = Mzsq+s
      t231 = Mzsq/t229
      spf20 = DENSPENCE(DCMPLX(t231),spcZip)
      t232 = -Mzsq+Ieps
      t233 = 1/t232
      spf26 = log(DCMPLX(-s*t233))
      spf29 = log(DCMPLX(1-t167*Ieps))
      spf30 = DENSPENCE(DCMPLX(-Mzsq*t233),spcZip)
      spf34 = log(DCMPLX(t165*t233))
      t239 = t229*t233
      spf45 = log(DCMPLX(-t239))
      spf46 = DENSPENCE(DCMPLX(-t239),spcZip)
      spf86 = log(DCMPLX(t174+(-t62-t57)*t5*t156))
      t245 = 1/(-t62-t57+t63)
      t246 = t62*t245
      spf87 = log(DCMPLX(-t246))
      spf88 = DENSPENCE(DCMPLX(-t246),spcZip)
      spf98 = log(DCMPLX(-t57*t245))
      spf112 = log(DCMPLX(s*t232*t245))
      t250 = Mzsq*t229
      spf117 = DENSPENCE(DCMPLX(-t250*t245),spcZip)
      t252 = log(DCMPLX(t231))
      t254 = log(DCMPLX(-Mzsq*t156))
      t257 = t5*t156
      t259 = log(DCMPLX(-t250*t257))
      t272 = 6*spf30
      t285 = -6*t252*t254+6*t252*t259-6*spf117+6*spf88+6*spf87*spf112-6*
     #spf87*spf29-6*spf87*spf98+t211-6*spf20+6*spf46-t272+6*spf45*spf34-
     #6*spf45*spf22-6*spf45*spf26+6*t252*spf32-6*t252*spf86+6*t252*spf27
      wV7 = -t285*t5/36
      t288 = -s+t7
      t290 = t288*t167*t168
      spf51 = DENSPENCE(DCMPLX(t290),spcZip)
      t291 = 2*t7
      t292 = s-t11
      t293 = sqrt(t292)
      t294 = t10*t293
      t296 = 1/(-s+t291-t294)
      spf140 = DENSPENCE(DCMPLX(2*t7*t296),spcZip)
      t300 = 1/(-s+t291+t294)
      spf141 = DENSPENCE(DCMPLX(2*t7*t300),spcZip)
      spf162 = DENSPENCE(DCMPLX(2*t288*t296),spcZip)
      spf163 = DENSPENCE(DCMPLX(2*t288*t300),spcZip)
      t309 = log(DCMPLX(t7/t288))
      t313 = log(DCMPLX(t290))
      t322 = -t292
      t324 = sqrt(-t322*t5)
      t325 = -1-t324
      t326 = log(DCMPLX(t325))
      t327 = t124*t326
      t331 = -t325
      t332 = s*t331
      t333 = t324*s
      t334 = -s+t291-t333
      t335 = 1/t334
      t337 = log(DCMPLX(-t332*t335))
      t341 = log(DCMPLX(t7*t335))
      t344 = 1-t324
      t345 = log(DCMPLX(t344))
      t346 = t124*t345
      t348 = log(DCMPLX(t331))
      t349 = t124*t348
      t351 = s*t167
      t352 = t351*t168
      t353 = log(DCMPLX(t352))
      t356 = -t344
      t357 = s*t356
      t359 = log(DCMPLX(-t357*t335))
      t363 = log(DCMPLX(t288*t335))
      t366 = -6*IM*t309*PI-6*IM*t313*PI-6*spf163-6*spf162+6*spf141+6*spf
     #140+6*spf51-t211+6*t327+6*t309*t326+6*t337*t124+6*t337*t341-6*t346
     #-6*t349+6*t313*t353-6*t359*t124-6*t359*t363
      t367 = log(DCMPLX(-t352))
      t369 = log(DCMPLX(-t8))
      t371 = -s+t291+t333
      t372 = 1/t371
      t374 = log(DCMPLX(t288*t372))
      t376 = log(DCMPLX(-s*t372))
      t380 = log(DCMPLX(-s*t335))
      t382 = log(DCMPLX(t356))
      t385 = t124*t382
      t387 = log(DCMPLX(t357*t372))
      t390 = log(DCMPLX(t7*t372))
      t393 = log(DCMPLX(-t371*t5))
      t400 = log(DCMPLX(-t334*t5))
      t403 = log(DCMPLX(t332*t372))
      t406 = -t313*t367+t309*t369+t374*t376+t374*t326+t363*t380+t363*t38
     #2+t309*t382+t385+t387*t124+t387*t390-t309*t393-t341*t348-t341*t380
     #-t390*t376-t390*t345-t309*t400-t403*t374-t403*t124
      wV8 = (t366+6*t406)*t5/36
      t410 = -t103+t7+Ieps
      t411 = 1/t410
      t412 = t7*t411
      spf85 = log(DCMPLX(-t412))
      t413 = t132*t411
      spf93 = log(DCMPLX(t413))
      t414 = t7-t103+Ieps-s
      spf124 = log(DCMPLX(t410/t414))
      t418 = 1/(-s+t291-t128+t130-t294)
      spf149 = log(DCMPLX(-t418*s))
      t421 = 1/(-s+t291-t128+t130+t294)
      spf150 = log(DCMPLX(-t421*s))
      t423 = 2*t6
      spf151 = log(DCMPLX(-t423-(-s+t291-t128-t294)*t5))
      spf152 = log(DCMPLX(-t423-(-s+t291-t128+t294)*t5))
      t430 = t410*t418
      spf182 = log(DCMPLX(t430))
      t431 = t410*t421
      spf183 = log(DCMPLX(t431))
      t432 = t414*t418
      spf190 = log(DCMPLX(t432))
      t433 = t414*t421
      spf191 = log(DCMPLX(t433))
      spf196 = DENSPENCE(DCMPLX(2*t430),spcZip)
      spf197 = DENSPENCE(DCMPLX(2*t431),spcZip)
      spf201 = DENSPENCE(DCMPLX(2*t432),spcZip)
      spf202 = DENSPENCE(DCMPLX(2*t433),spcZip)
      t438 = -s+t294
      spf210 = log(DCMPLX(t438*t421))
      t440 = s+t294
      spf211 = log(DCMPLX(t440*t421))
      t442 = Ieps**2
      t445 = t103*Mzsq
      t446 = t445*cw2
      t447 = 2*t446
      t448 = s*t103
      t449 = t103**2
      t450 = t442+(-t128+t291-s)*Ieps-t447+t178+t448+t449
      t451 = 1/t450
      spf212 = log(DCMPLX(t7*s*t451))
      spf214 = log(DCMPLX(-t438*t418))
      spf215 = log(DCMPLX(-t440*t418))
      t456 = s*t132
      spf218 = log(DCMPLX(-t456*t451))
      t458 = -t103+t7
      t460 = 2*t458*Ieps
      t461 = t458**2
      t462 = t442+t460+t461
      t463 = t462*t451
      spf224 = log(DCMPLX(t463))
      spf226 = DENSPENCE(DCMPLX(t463),spcZip)
      t464 = s*t410
      spf227 = log(DCMPLX(-t464*t451))
      spf228 = log(DCMPLX(-t450*t5*t411))
      t469 = t414*t410*t451
      spf238 = log(DCMPLX(t469))
      spf239 = DENSPENCE(DCMPLX(t469),spcZip)
      t481 = -spf226+spf224*spf227+spf215*spf182-spf124*spf152+spf190*t3
     #82+spf190*spf149-spf183*t345-spf183*spf150-spf182*t348-spf182*spf1
     #49-spf124*spf151+spf215*t124-spf202-spf201+spf197+spf196+t327-t346
     #-t349
      t499 = t385+spf124*spf228-spf211*spf191-spf211*t124+spf238*spf212-
     #spf238*spf85-spf238*spf227+spf210*spf183+spf210*t124+spf191*spf150
     #+spf191*t326-spf224*spf218+spf224*spf93+spf239+spf124*t326-spf214*
     #t124+spf124*t382-spf124*spf85-spf214*spf190
      wV9 = (t481+t499)*t5/4
      t502 = -s+Mzsq
      t503 = t502*t167
      spf21 = DENSPENCE(DCMPLX(t503),spcZip)
      t504 = 2*Mzsq
      t505 = s-t34
      t506 = sqrt(t505)
      t507 = t10*t506
      t509 = 1/(-s+t504-t507)
      spf108 = DENSPENCE(DCMPLX(2*Mzsq*t509),spcZip)
      t513 = 1/(-s+t504+t507)
      spf109 = DENSPENCE(DCMPLX(2*Mzsq*t513),spcZip)
      spf118 = DENSPENCE(DCMPLX(2*t502*t509),spcZip)
      spf119 = DENSPENCE(DCMPLX(2*t502*t513),spcZip)
      t520 = -t505
      t522 = sqrt(-t520*t5)
      t523 = 1+t522
      t524 = s*t523
      t525 = t522*s
      t526 = -s+t504-t525
      t527 = 1/t526
      t529 = log(DCMPLX(-t524*t527))
      t531 = log(DCMPLX(Mzsq*t527))
      t533 = log(DCMPLX(t523))
      t535 = -s+t504+t525
      t536 = 1/t535
      t538 = log(DCMPLX(t524*t536))
      t540 = log(DCMPLX(t502*t536))
      t543 = 1-t522
      t544 = log(DCMPLX(t543))
      t546 = -t543
      t547 = s*t546
      t549 = log(DCMPLX(t547*t536))
      t552 = log(DCMPLX(Mzsq*t536))
      t555 = log(DCMPLX(t502*t527))
      t556 = log(DCMPLX(t546))
      t559 = log(DCMPLX(-s*t527))
      t564 = log(DCMPLX(-t547*t527))
      t570 = log(DCMPLX(Mzsq/t502))
      t574 = -t529*t531+t531*t533+t538*t540+t538*t124+t124*t544-t549*t12
     #4-t549*t552-t555*t556-t555*t559-t529*t124+t531*t559+t564*t124+t564
     #*t555+t552*t544-t570*t556-t124*t556+t124*t533
      t577 = log(DCMPLX(-t523))
      t582 = log(DCMPLX(-t32))
      t586 = log(DCMPLX(-t535*t5))
      t590 = log(DCMPLX(-s*t536))
      t594 = log(DCMPLX(-t526*t5))
      t601 = log(DCMPLX(t503))
      t602 = log(DCMPLX(t351))
      t605 = log(DCMPLX(-t351))
      t619 = -6*t570*t577-6*t124*t577-6*t570*t582+6*t570*t586+6*t552*t59
     #0+6*t570*t594-6*t540*t590-6*t540*t577-6*t601*t602+6*t601*t605+6*IM
     #*t601*PI+6*IM*t570*PI+6*spf119+6*spf118-6*spf109-6*spf108-6*spf21+
     #t211
      wV10 = -(6*t574+t619)*t5/36
      spf64 = log(DCMPLX(-t111))
      spf65 = log(DCMPLX(t109))
      spf78 = log(DCMPLX(t412))
      spf101 = log(DCMPLX(-t413))
      t623 = -t103+Ieps+t7+s
      spf123 = log(DCMPLX(t410/t623))
      t627 = 1/(s-t128+t130+t291-t107)
      spf135 = log(DCMPLX(s*t627))
      t630 = 1/(s-t128+t130+t291+t107)
      spf136 = log(DCMPLX(s*t630))
      spf137 = log(DCMPLX(t423+(s-t128+t291-t107)*t5))
      spf138 = log(DCMPLX(t423+(s-t128+t291+t107)*t5))
      t638 = t410*t627
      spf173 = log(DCMPLX(t638))
      t639 = t410*t630
      spf174 = log(DCMPLX(t639))
      t640 = t623*t627
      spf176 = log(DCMPLX(t640))
      t641 = t623*t630
      spf177 = log(DCMPLX(t641))
      spf180 = DENSPENCE(DCMPLX(2*t638),spcZip)
      spf181 = DENSPENCE(DCMPLX(2*t639),spcZip)
      spf186 = DENSPENCE(DCMPLX(2*t640),spcZip)
      spf187 = DENSPENCE(DCMPLX(2*t641),spcZip)
      spf203 = log(DCMPLX(t110*t630))
      spf204 = log(DCMPLX(t108*t630))
      t648 = t442+t460+t449-t447+t176+t178
      t649 = 1/t648
      spf205 = log(DCMPLX(t7*s*t649))
      spf206 = log(DCMPLX(-t110*t627))
      spf207 = log(DCMPLX(-t108*t627))
      spf216 = log(DCMPLX(-t456*t649))
      t655 = t462*t649
      spf219 = log(DCMPLX(t655))
      spf220 = DENSPENCE(DCMPLX(t655),spcZip)
      spf222 = log(DCMPLX(t464*t649))
      spf223 = log(DCMPLX(t648*t5*t411))
      t660 = t623*t410*t649
      spf234 = log(DCMPLX(t660))
      spf235 = DENSPENCE(DCMPLX(t660),spcZip)
      t680 = -t124*spf64+t124*spf63-spf173*spf64-spf173*spf135+spf123*sp
     #f223+t124*spf62-t124*spf65+spf234*spf216-spf203*spf177-spf203*t124
     #+spf176*spf62+spf176*spf135-spf234*spf101-spf234*spf222-spf219*spf
     #205+spf219*spf78+spf219*spf222-spf123*spf138+spf177*spf63
      t694 = spf204*spf174+spf204*t124+spf206*t124+spf206*spf173-spf174*
     #spf136-spf174*spf65-spf123*spf101+spf123*spf63+spf123*spf62+spf177
     #*spf136-spf207*spf176-spf207*t124-spf123*spf137+spf235-spf220+spf1
     #80+spf181-spf186-spf187
      wV11 = (t680+t694)*t5/4
      t697 = s+t129
      spf69 = log(DCMPLX(t130+t697*t156))
      t700 = s*t152
      spf70 = log(DCMPLX(t700*t160))
      t702 = t5*t153
      t703 = t702*Ieps
      t704 = t702*t7
      spf71 = log(DCMPLX(t703-t704))
      t707 = t697*t5
      spf84 = log(DCMPLX(2*t703-t707*t153))
      t710 = -t291+s+t129
      t711 = 1/t710
      t712 = t7*t711
      spf97 = DENSPENCE(DCMPLX(-2*t712),spcZip)
      t714 = s+t129-t130
      spf106 = log(DCMPLX(-t714*t160))
      t716 = t710*t160
      spf114 = log(DCMPLX(t716))
      spf122 = DENSPENCE(DCMPLX(t716/2),spcZip)
      t718 = t697*Ieps
      t719 = 2*t178
      t720 = Mzsq*CsTh
      t721 = t175*t720
      t723 = 1/(t718+t719-t176-t721)
      t724 = t178*t723
      spf172 = log(DCMPLX(t724))
      spf175 = DENSPENCE(DCMPLX(2*t724),spcZip)
      t726 = t153*t156
      t727 = t726*Ieps
      t728 = t707*t727
      spf198 = log(DCMPLX(t728+(t719-t176-t721)*t5*t726))
      spf208 = log(DCMPLX(-t700*t7*t723))
      spf217 = log(DCMPLX(t700*t159*t723))
      spf221 = DENSPENCE(DCMPLX(-t7*t710*t723),spcZip)
      t746 = t702*t156
      t748 = log(DCMPLX(-t7*t710*t746))
      t751 = log(DCMPLX(-t712))
      t764 = 6*spf114*spf70-6*spf114*spf106+6*spf84*spf114-6*t748*t124-6
     #*t748*t751-6*spf217*t124-6*spf217*spf172-6*t124*spf70+t214+6*spf22
     #1-6*spf175-6*spf122+6*spf97-t211
      t768 = spf84*t124
      t769 = spf67*t124
      t770 = t124*spf69
      t777 = t194*t124
      t779 = spf198*t124+spf198*t751+spf67*spf172-t768+t769-t770+t194*t7
     #51-t751*spf71-t751*spf69+t124*spf106+spf208*t124+spf208*spf172+t77
     #7-t124*spf71
      wV12 = -(t764+6*t779)*t5*t153/18
      t785 = -t291+s+t129+t128-t130
      t787 = t410/t785
      spf147 = log(DCMPLX(-t787))
      spf148 = DENSPENCE(DCMPLX(-2*t787),spcZip)
      t789 = t32*t153
      spf153 = log(DCMPLX(2*t789*cw2*t156*Ieps-t7*(-t291+s+t129+t128)*t7
     #46))
      t801 = t129*t103
      t803 = -2*t442+(t136+t135+t104)*Ieps-t801-2*t449-t448
      t806 = t153/t131
      spf225 = log(DCMPLX(t803*t5*t806))
      t808 = t291+s+t129
      spf229 = log(DCMPLX(t808*t5*t727+(-t447+t719-t176-t721)*t5*t726))
      t815 = t808*Ieps
      t817 = 1/(t815-t447+t719-t176-t721)
      spf230 = log(DCMPLX(-t700*t7*t817))
      t820 = t815-t447-t176-t721
      t821 = 1/t820
      spf231 = log(DCMPLX(t700*t131*t821))
      spf232 = log(DCMPLX(t820*t5*t806))
      spf233 = log(DCMPLX(t700*t159*t817))
      t830 = t410*Mzsq*cw2*t817
      spf236 = log(DCMPLX(t830))
      spf237 = DENSPENCE(DCMPLX(2*t830),spcZip)
      t838 = t54*CsTh
      t839 = 2*t838
      t840 = 4*t448
      t841 = 4*t801
      t842 = CsTh**2
      t843 = t54*t842
      t844 = 4*t449
      t845 = 4*t442+(-6*s-6*t129-8*t103)*Ieps+t54+t839+t840+t841+t843+t8
     #44
      spf240 = log(DCMPLX(-t845*t5*t806))
      t849 = t410*t131*t821
      spf241 = log(DCMPLX(-t849))
      spf242 = DENSPENCE(DCMPLX(-t849),spcZip)
      spf243 = DENSPENCE(DCMPLX(-t7*t785*t817),spcZip)
      t853 = t785*t131*t821
      spf244 = log(DCMPLX(t853))
      spf245 = DENSPENCE(DCMPLX(t853/2),spcZip)
      spf248 = log(DCMPLX(t803*t821))
      spf249 = log(DCMPLX(-t845*t821))
      t884 = -6*t124*spf231+6*spf230*t124-6*spf244*spf249-6*spf153*spf14
     #7-6*spf153*t124-6*spf240*t124+6*spf229*spf147+6*spf229*t124+6*spf1
     #48-t211+6*spf241*spf248+6*spf230*spf236-6*spf233*spf236+6*spf67*sp
     #f236-6*spf147*spf232
      t894 = -spf241*spf231-spf241*spf225-t124*spf232-spf147*spf69-spf23
     #3*t124+t124*spf249+t194*spf147+spf240*spf244+spf244*spf231-spf245+
     #t769-t770+spf243+spf242+t777-spf237
      wV13 = -(t884+6*t894)*t5*t153/18
      spf43 = log(DCMPLX(t700*t233))
      spf44 = log(DCMPLX(t703-t789))
      t902 = -t504+s+t129
      t903 = 1/t902
      t904 = Mzsq*t903
      spf59 = DENSPENCE(DCMPLX(-2*t904),spcZip)
      t906 = t902*t233
      spf72 = log(DCMPLX(t906))
      spf81 = log(DCMPLX(-t714*t233))
      spf89 = DENSPENCE(DCMPLX(t906/2),spcZip)
      t909 = 2*t62
      t910 = t57*CsTh
      t912 = 1/(t718+t909-t57-t910)
      t913 = t62*t912
      spf144 = log(DCMPLX(t913))
      spf159 = DENSPENCE(DCMPLX(2*t913),spcZip)
      spf169 = log(DCMPLX(t728+(t909-t57-t910)*t5*t726))
      spf170 = log(DCMPLX(-t57*t152*t912))
      spf179 = log(DCMPLX(t700*t232*t912))
      t923 = Mzsq*t902
      spf188 = DENSPENCE(DCMPLX(-t923*t912),spcZip)
      t925 = log(DCMPLX(-t904))
      t935 = log(DCMPLX(-t923*t746))
      t944 = t254*t124
      t952 = -t272+t211-6*t254*t925+6*t925*spf44-6*spf169*t925-6*spf84*s
     #pf72+6*t935*t925+6*t124*spf44-6*spf169*t124+6*t124*spf43-6*t944+6*
     #t935*t124-6*spf170*spf144-6*spf170*t124
      t960 = spf29*t124
      t961 = -t124*spf81-spf72*spf43+spf179*spf144+spf179*t124+spf72*spf
     #81-spf188+t925*spf69-spf29*spf144-t960+spf159+spf89-spf59+t768+t77
     #0
      wV14 = (t952+6*t961)*t5*t153/18
      t966 = -1+CsTh
      t967 = t966*s
      spf57 = log(DCMPLX(-t967*t233))
      t969 = 1/t966
      t970 = t969*t5
      t971 = t970*Ieps
      t972 = t32*t969
      spf58 = log(DCMPLX(-t971+t972))
      t974 = t504-s+t129
      t975 = 1/t974
      t976 = Mzsq*t975
      spf60 = DENSPENCE(DCMPLX(2*t976),spcZip)
      t978 = -s+t129+t130
      spf73 = log(DCMPLX(t978*t233))
      t980 = -s+t129
      spf79 = log(DCMPLX(t130-t980*t156))
      t983 = t974*t233
      spf80 = log(DCMPLX(-t983))
      t985 = t980*t969
      spf83 = log(DCMPLX(-2*t971-t985*t5))
      spf90 = DENSPENCE(DCMPLX(-t983/2),spcZip)
      t991 = 1/(t980*Ieps-t909+t57-t910)
      t992 = t62*t991
      spf157 = log(DCMPLX(-t992))
      spf158 = DENSPENCE(DCMPLX(-2*t992),spcZip)
      spf168 = log(DCMPLX(t985*t257*Ieps+(-t909+t57-t910)*t5*t969*t156))
      spf171 = log(DCMPLX(-t967*Mzsq*t991))
      spf178 = log(DCMPLX(t967*t232*t991))
      t1005 = Mzsq*t974
      spf189 = DENSPENCE(DCMPLX(-t1005*t991),spcZip)
      t1014 = log(DCMPLX(t976))
      t1025 = log(DCMPLX(-t1005*t970*t156))
      t1034 = -t272-6*spf189+6*t124*spf57+t211+6*spf178*spf157-6*spf29*s
     #pf157-6*t254*t1014+6*t124*spf79+6*t1014*spf79+6*spf178*t124+6*t102
     #5*t1014+6*spf80*spf73-6*spf80*spf57-6*spf168*t124
      t1044 = -spf168*t1014+t1025*t124+t1014*spf58-spf171*spf157-t124*sp
     #f73+spf83*t124-spf83*spf80-spf171*t124+t124*spf58+spf158+spf90-spf
     #60-t944-t960
      wV15 = -(t1034+6*t1044)*t969*t5/18
      t1050 = sqrt(t152)
      t1051 = 2*t721
      t1052 = 4*t178
      t1053 = 2*t176
      t1056 = 4*t176
      t1057 = 4*t721
      t1058 = 8*t178
      t1059 = t54+t838-t1056-t1057+t1058
      t1060 = sqrt(t1059)
      t1061 = 1/t1060
      t1063 = t711*t167
      t1065 = t1063*t168*Ieps
      t1066 = t1050*(t1051-t1052+t1053)*t1061*t1065
      t1067 = t54*cw2
      t1068 = t1067*t720
      t1069 = 4*t1068
      t1070 = t54*s
      t1071 = t1070*CsTh
      t1072 = s*t62
      t1073 = t1072*t177
      t1074 = 8*t1073
      t1075 = t54*t1050
      t1076 = t1075*t1060
      t1077 = t54*Mzsq
      t1078 = t1077*cw2
      t1079 = 4*t1078
      t1083 = t1063*t168
      spf256 = DENSPENCE(DCMPLX(-t1066/2-t1050*(-t1069+t1071+t1074-t1076
     #-t1079+t1070)*t1061*t1083/2),spcZip)
      spf257 = DENSPENCE(DCMPLX(t1066/2+t1050*(-t1069+t1071+t1074+t1076-
     #t1079+t1070)*t1061*t1083/2),spcZip)
      t1093 = CsTh*t62
      t1094 = t1093*t177
      t1095 = 4*t1094
      t1096 = Mzsq*t842
      t1097 = t175*t1096
      t1098 = 2*t1097
      t1100 = 1/t1050
      t1102 = (-t1052-t1095+t1053+t1057+t1098)*t1100*t1061
      t1103 = t1102*t1065
      t1104 = t1050*t1060
      t1106 = 4*t1104*t178
      t1107 = s*t1050
      t1111 = 2*t1107*t1060*Mzsq*cw2
      t1113 = t7*CsTh
      t1115 = 2*t1107*t1060*t1113
      t1117 = t1075*t1060*CsTh
      t1118 = 2*t1071
      t1119 = t1070*t842
      t1120 = 8*t1068
      t1121 = t843*t7
      t1122 = 4*t1121
      t1123 = t129*t178
      t1124 = 8*t1123
      t1125 = -t1106+t1111+t1115-t1076-t1117+t1070+t1118+t1119-t1079-t11
     #20-t1122+t1074+t1124
      t1126 = t1125*t1100
      spf287 = log(DCMPLX(t1103+t1126*t1061*t1083))
      t1130 = t1106-t1111-t1115+t1076+t1117+t1070+t1118+t1119-t1079-t112
     #0-t1122+t1074+t1124
      t1131 = t1130*t1100
      spf289 = log(DCMPLX(-t1103-t1131*t1061*t1083))
      t1137 = t1102*t711*t300*Ieps
      t1140 = cw2*t10*t293
      t1142 = 2*t1104*Mzsq*t1140
      t1143 = t10*s
      t1144 = t1143*t1050
      t1145 = t1060*t293
      t1146 = t1144*t1145
      t1148 = t1144*t1145*CsTh
      t1149 = t1111-t1142+t1146+t1148+t1070+t1118+t1119-t1079-t1120-t112
     #2+t1074+t1124
      t1151 = t1061*t711
      t1152 = t1151*t300
      t1154 = t1137+t1149*t1100*t1152
      spf292 = log(DCMPLX(t1154))
      t1157 = t1102*t711*t296*Ieps
      t1158 = t1111+t1142-t1146-t1148+t1070+t1118+t1119-t1079-t1120-t112
     #2+t1074+t1124
      t1160 = t1151*t296
      t1162 = t1157+t1158*t1100*t1160
      spf293 = log(DCMPLX(t1162))
      spf296 = DENSPENCE(DCMPLX(t1154),spcZip)
      spf297 = DENSPENCE(DCMPLX(t1162),spcZip)
      t1163 = t1070+t1118+t1119-t1079-t1120-t1122+t1074+t1124+t1148-t111
     #1-t1142+t1146
      t1166 = -t1157-t1163*t1100*t1160
      spf298 = log(DCMPLX(t1166))
      t1167 = t1070+t1118+t1119-t1079-t1120-t1122+t1074+t1124-t1148-t111
     #1+t1142-t1146
      t1170 = -t1137-t1167*t1100*t1152
      spf299 = log(DCMPLX(t1170))
      spf302 = DENSPENCE(DCMPLX(t1166),spcZip)
      spf303 = DENSPENCE(DCMPLX(t1170),spcZip)
      spf308 = log(DCMPLX(t1157+t1126*t1160))
      spf309 = log(DCMPLX(t1137+t1126*t1152))
      spf312 = log(DCMPLX(-t1157-t1131*t1160))
      spf313 = log(DCMPLX(-t1137-t1131*t1152))
      t1181 = sqrt(-s*t322)
      t1183 = 1/(t291-s-t1181)
      t1185 = log(DCMPLX(t7*t1183))
      t1191 = log(DCMPLX(t7/(t291-s+t1181)))
      t1194 = log(DCMPLX(t704))
      t1201 = -spf298*spf287-spf298*t1185-spf299*spf287-spf299*t1191+spf
     #287*t124+spf287*t1194+spf293*spf289+spf293*t1185+spf292*spf289+spf
     #292*t1191-spf289*t124
      t1209 = -spf289*t1194+spf298*spf308+spf299*spf309-spf293*spf312-sp
     #f292*spf313-2*spf256+spf302+spf303+2*spf257-spf297-spf296
      t1211 = t54*t152
      t1214 = sqrt(t1211*t1059)
      wV16 = -(t1201+t1209)/t1214/2
      spf28 = -t458
      t1218 = 1/t10
      t1219 = t7*t1218
      t1224 = CsTh*t103
      t1225 = t1224*t7
      t1227 = -t1058-8*t1094+t1056+8*t721+4*t1097+8*t446+8*t1225-t844-t8
     #40-t841-t54-t839-t843
      t1229 = sqrt(t1227*t296)
      t1230 = 1/t1229
      t1231 = t128+s+t129
      t1232 = 1/t1231
      t1233 = t1230*t1232
      t1236 = 16*t1219*t1233*Ieps
      t1237 = 2*t1119
      t1238 = 16*t1073
      t1239 = 16*t1123
      t1240 = t175*t445
      t1241 = 16*t1240
      t1242 = t129*t446
      t1243 = 16*t1242
      t1244 = 8*t1078
      t1245 = 16*t1068
      t1246 = 8*t1121
      t1247 = t449*s
      t1248 = 8*t1247
      t1249 = t54*t103
      t1250 = 8*t1249
      t1251 = t1249*CsTh
      t1252 = 8*t1251
      t1253 = 2*t1070
      t1254 = 4*t1071
      t1255 = sqrt(2.D0)
      t1256 = t10*t54
      t1257 = t1255*t1256
      t1258 = t1257*t1229
      t1259 = t1255*s
      t1260 = t1229*t103
      t1263 = 2*t1259*t1260*t293
      t1264 = t1255*t54
      t1265 = t1229*t293
      t1266 = t1264*t1265
      t1268 = t1264*t1265*CsTh
      t1269 = t1229*CsTh
      t1270 = t1257*t1269
      t1271 = t1255*t1143
      t1273 = 2*t1271*t1260
      t1274 = -t1237-t1238-t1239+t1241+t1243+t1244+t1245+t1246-t1248-t12
     #50-t1252-t1253-t1254+t1258+t1263+t1266+t1268+t1270+t1273
      t1275 = t1274*t1218
      t1277 = t169*t1232
      spf314 = log(DCMPLX(-t1236+t1275*t1230*t1277))
      t1280 = -t1238-t1239-t1237+t1244+t1245+t1246+t1241+t1243-t1263-t12
     #48-t1266-t1268-t1273-t1253-t1254-t1252-t1258-t1270-t1250
      t1281 = t1280*t1218
      spf315 = log(DCMPLX(t1236-t1281*t1230*t1277))
      t1286 = t296*Ieps
      t1287 = t1233*t1286
      t1289 = 16*t178*t1218*t1287
      t1290 = t1233*t296
      spf316 = log(DCMPLX(-t1289+t1281*t1290))
      spf317 = log(DCMPLX(t1289-t1275*t1290))
      t1295 = t1255*Mzsq
      t1297 = t1218*t1230
      t1301 = 2*t1295*cw2*t1297*t1232*Ieps
      t1302 = t1259*CsTh
      t1304 = 8*t1302*t446
      t1305 = t54*t1229
      t1306 = t293*CsTh
      t1307 = t1305*t1306
      t1309 = 4*t1264*t1224
      t1311 = 4*t1264*t7
      t1312 = t1255*t1070
      t1314 = 2*t1312*CsTh
      t1316 = 8*t1259*t178
      t1317 = t1255*cw2
      t1318 = t448*Mzsq
      t1320 = 8*t1317*t1318
      t1322 = 8*t1264*t1113
      t1323 = t842*cw2
      t1324 = t1323*Mzsq
      t1326 = 4*t1264*t1324
      t1328 = 8*t1259*t1094
      t1329 = t1256*t1229
      t1330 = s*t1229
      t1331 = t103*t293
      t1333 = 2*t1330*t1331
      t1334 = t1305*t293
      t1335 = t1143*t1229
      t1337 = 2*t1335*t103
      t1338 = t1329*CsTh
      t1339 = t1255*t449
      t1341 = 4*t1339*s
      t1342 = t1312*t842
      t1344 = 4*t1264*t103
      t1345 = t1304-t1307-t1309+t1311-t1314-t1316+t1320+t1322+t1326-t132
     #8-t1312-t1329-t1333-t1334-t1337-t1338-t1341-t1342-t1344
      t1346 = t1345*t1218
      spf318 = DENSPENCE(DCMPLX(t1301-t1346*t1230*t1277/4),spcZip)
      t1351 = t1304+t1307-t1309+t1311-t1314-t1316+t1320+t1322+t1326-t132
     #8-t1312+t1329+t1333+t1334+t1337+t1338-t1341-t1342-t1344
      t1352 = t1351*t1218
      spf319 = DENSPENCE(DCMPLX(-t1301+t1352*t1230*t1277/4),spcZip)
      t1357 = t1255*t62
      t1361 = 4*t1357*t177*t1218*t1287
      spf320 = DENSPENCE(DCMPLX(-t1361+t1346*t1290/2),spcZip)
      spf321 = DENSPENCE(DCMPLX(t1361-t1352*t1290/2),spcZip)
      t1369 = 16*t1297*Ieps
      t1371 = t1255*t10*t1229
      t1373 = 8*t1371*t446
      t1377 = 4*t1271*t1229*Mzsq*cw2
      t1378 = t1271*t1229
      t1380 = 4*t1378*t1113
      t1381 = -t1238-t1239-t1253+t1373+t1377+t1380+t1244+t1245+t1246+t12
     #41+t1243
      t1382 = -t1266-t1268-t1248-t1263-t1270-t1254-t1237-t1258-t1250-t12
     #52-t1273
      t1383 = t1381+t1382
      t1384 = t1383*t1218
      t1385 = 1/t62
      t1387 = 1/t177
      t1388 = t1230*t1385*t1387
      spf322 = log(DCMPLX(-t1369+t1384*t1388))
      t1391 = -t1238-t1239+t1273-t1377-t1380+t1241+t1243+t1244+t1245+t12
     #46-t1373
      t1392 = t1258+t1270+t1266+t1268-t1253-t1250-t1252-t1248-t1254-t123
     #7+t1263
      t1393 = t1391+t1392
      t1394 = t1393*t1218
      spf323 = log(DCMPLX(t1369-t1394*t1388))
      t1397 = t1384*t1230
      spf324 = log(DCMPLX(-t1236+t1397*t1277))
      t1400 = t1394*t1230
      spf325 = log(DCMPLX(t1236-t1400*t1277))
      t1406 = 16*t1219*t1230*t296*Ieps
      t1407 = t169*t296
      spf326 = log(DCMPLX(-t1406+t1397*t1407))
      spf327 = log(DCMPLX(t1406-t1400*t1407))
      spf328 = log(DCMPLX(-t1289+t1384*t1290))
      spf329 = log(DCMPLX(t1289-t1394*t1290))
      t1416 = t7*spf28
      t1417 = t103*spf28
      t1418 = spf28**2
      t1419 = SQRT(t1418)
      t1420 = t1419*t103
      t1422 = t1419*Mzsq*cw2
      t1424 = 1/(t1416+t1417+t1420-t1422)
      t1426 = t1297*t1424*Ieps
      t1428 = 16*t1416*t1426
      t1430 = t1383*spf28*t1218
      t1431 = t1230*t1424
      t1432 = t1431*t169
      spf330 = log(DCMPLX(-t1428+t1430*t1432))
      t1436 = 1/(-t1416-t1417+t1420-t1422)
      t1438 = t1297*t1436*Ieps
      t1440 = 16*t1416*t1438
      t1442 = t1393*spf28*t1218
      t1443 = t1230*t1436
      t1444 = t1443*t169
      spf331 = log(DCMPLX(-t1440+t1442*t1444))
      spf332 = log(DCMPLX(t1440-t1430*t1444))
      spf333 = log(DCMPLX(t1428-t1442*t1432))
      t1452 = t178*spf28*t1218
      t1455 = 16*t1452*t1431*t1286
      t1456 = t1431*t296
      spf334 = log(DCMPLX(-t1455+t1430*t1456))
      t1461 = 16*t1452*t1443*t1286
      t1462 = t1443*t296
      spf335 = log(DCMPLX(-t1461+t1442*t1462))
      spf336 = log(DCMPLX(t1461-t1430*t1462))
      spf337 = log(DCMPLX(t1455-t1442*t1456))
      t1471 = 4*t1371*t7*t1420
      t1472 = CsTh*spf28
      t1475 = 4*t1378*t7*t1472
      t1476 = cw2*spf28
      t1477 = t445*t1476
      t1479 = 4*t1371*t1477
      t1481 = 2*t1119*spf28
      t1483 = t1257*t1229*spf28
      t1485 = 8*t1249*t1472
      t1487 = 8*t1077*t1476
      t1488 = t177*spf28
      t1490 = 16*t1072*t1488
      t1491 = t178*t1419
      t1493 = 4*t1371*t1491
      t1494 = t178*spf28
      t1496 = 4*t1371*t1494
      t1498 = 16*t801*t1416
      t1500 = t1306*spf28
      t1501 = t1264*t1229*t1500
      t1502 = t1471+t1475+t1479-t1481-t1483-t1485+t1487-t1490-t1493-t149
     #6+t1498-t1501
      t1503 = t1259*t1229
      t1504 = t1331*spf28
      t1506 = 2*t1503*t1504
      t1507 = t1378*t1416
      t1508 = 4*t1507
      t1510 = 8*t1247*spf28
      t1512 = 16*t129*t1494
      t1515 = 2*t1271*t1260*spf28
      t1517 = t1264*t1265*spf28
      t1519 = t1257*t1269*spf28
      t1521 = 8*t843*t1416
      t1522 = t720*spf28
      t1524 = 16*t1067*t1522
      t1525 = t445*spf28
      t1527 = 16*t175*t1525
      t1529 = 2*t1070*spf28
      t1531 = 8*t1249*spf28
      t1533 = 4*t1071*spf28
      t1534 = -t1506+t1508-t1510-t1512-t1515-t1517-t1519+t1521+t1524+t15
     #27-t1529-t1531-t1533
      t1539 = t1436*t167*t168
      spf338 = log(DCMPLX(-t1440+(t1502+t1534)*t1218*t1230*t1539))
      t1542 = t1471-t1475-t1479-t1481+t1483-t1485+t1487-t1490-t1493+t149
     #6+t1498+t1501
      t1543 = t1506-t1508-t1510-t1512+t1515+t1517+t1519+t1521+t1524+t152
     #7-t1529-t1531-t1533
      t1548 = t1424*t167*t168
      spf339 = log(DCMPLX(-t1428+(t1542+t1543)*t1218*t1230*t1548))
      t1551 = -t1471+t1475+t1479-t1481-t1483-t1485+t1487-t1490+t1493-t14
     #96+t1498-t1501
      spf340 = log(DCMPLX(t1428-(t1551+t1534)*t1218*t1230*t1548))
      t1557 = -t1471-t1475-t1479-t1481+t1483-t1485+t1487-t1490+t1493+t14
     #96+t1498+t1501
      spf341 = log(DCMPLX(t1440-(t1557+t1543)*t1218*t1230*t1539))
      t1563 = t1295*t1476
      t1565 = 2*t1563*t1438
      t1566 = t1312*spf28
      t1568 = 8*t1302*t1477
      t1569 = t293*spf28
      t1570 = t1305*t1569
      t1572 = 4*t1264*t1417
      t1574 = t1312*t842*spf28
      t1576 = 2*t1312*t1472
      t1579 = 4*t1339*s*spf28
      t1581 = 2*t1335*t1417
      t1582 = t1329*t1472
      t1583 = t10*t1229
      t1586 = 4*t1583*t103*t1416
      t1589 = 4*t1583*t1419*t446
      t1591 = 8*t1302*t1494
      t1592 = -t1566+t1568-t1570-t1572-t1574-t1576-t1579-t1581-t1582+t15
     #86+t1589-t1591
      t1595 = 4*t1264*t842*t1416
      t1598 = 8*t1264*cw2*t1522
      t1601 = 8*t1317*s*t1525
      t1603 = cw2*CsTh
      t1606 = 4*t1335*Mzsq*t1603*spf28
      t1607 = t1329*spf28
      t1610 = 4*t1264*t1224*spf28
      t1611 = t1305*t1500
      t1613 = 2*t1330*t1504
      t1615 = 4*t1583*t1491
      t1617 = 4*t1583*t1494
      t1619 = 8*t1259*t1494
      t1621 = 4*t1264*t1416
      t1622 = t1335*t1416
      t1623 = 4*t1622
      t1624 = t1595+t1598+t1601+t1606-t1607-t1610-t1611-t1613-t1615-t161
     #7-t1619+t1621+t1623
      spf342 = DENSPENCE(DCMPLX(-t1565+(t1592+t1624)*t1218*t1230*t1539/4
     #),spcZip)
      t1632 = 2*t1563*t1426
      t1633 = -t1566+t1568-t1570-t1572-t1574-t1576-t1579-t1581-t1582+t15
     #86-t1589-t1591
      t1634 = t1595+t1598+t1601+t1606-t1607-t1610-t1611-t1613+t1615-t161
     #7-t1619+t1621+t1623
      spf343 = DENSPENCE(DCMPLX(t1632-(t1633+t1634)*t1218*t1230*t1548/4)
     #,spcZip)
      t1641 = -t1566+t1568+t1570-t1572-t1574-t1576-t1579+t1581+t1582-t15
     #86+t1589-t1591
      t1642 = t1595+t1598+t1601-t1606+t1607-t1610+t1611+t1613-t1615+t161
     #7-t1619+t1621-t1623
      spf344 = DENSPENCE(DCMPLX(-t1632+(t1641+t1642)*t1218*t1230*t1548/4
     #),spcZip)
      t1649 = -t1566+t1568+t1570-t1572-t1574-t1576-t1579+t1581+t1582-t15
     #86-t1589-t1591
      t1650 = t1595+t1598+t1601-t1606+t1607-t1610+t1611+t1613+t1615+t161
     #7-t1619+t1621-t1623
      spf345 = DENSPENCE(DCMPLX(t1565-(t1649+t1650)*t1218*t1230*t1539/4)
     #,spcZip)
      t1657 = t1357*t1488
      t1662 = 4*t1657*t1297*t1436*t296*Ieps
      t1667 = 2*t1330*Mzsq*cw2*t293*spf28
      t1671 = 2*t1330*t1419*t7*t293
      t1672 = -t1566+t1568-t1570-t1572-t1574-t1576-t1579-t1582+t1667+t16
     #71+t1586+t1589-t1591+t1595
      t1673 = 6*t1622
      t1675 = 2*t1335*t1422
      t1676 = t1420*t293
      t1678 = 2*t1330*t1676
      t1680 = 2*t1335*t1420
      t1681 = t1598+t1601+t1606-t1607-t1610-t1611-t1615-t1617-t1619+t162
     #1+t1673+t1675-t1678-t1680
      spf346 = DENSPENCE(DCMPLX(-t1662+(t1672+t1681)*t1218*t1462/2),spcZ
     #ip)
      t1691 = 4*t1657*t1297*t1424*t296*Ieps
      t1692 = -t1566+t1568-t1570-t1572-t1574-t1576-t1579-t1582+t1667-t16
     #71+t1586-t1589-t1591+t1595
      t1693 = t1598+t1601+t1606-t1607-t1610-t1611+t1615-t1617-t1619+t162
     #1+t1673-t1675+t1678+t1680
      spf347 = DENSPENCE(DCMPLX(t1691-(t1692+t1693)*t1218*t1456/2),spcZi
     #p)
      t1699 = -t1566+t1568+t1570-t1572-t1574-t1576-t1579+t1582-t1667+t16
     #71-t1586+t1589-t1591+t1595
      t1700 = t1598+t1601-t1606+t1607-t1610+t1611-t1615+t1617-t1619+t162
     #1-t1673+t1675-t1678-t1680
      spf348 = DENSPENCE(DCMPLX(-t1691+(t1699+t1700)*t1218*t1456/2),spcZ
     #ip)
      t1706 = -t1566+t1568+t1570-t1572-t1574-t1576-t1579+t1582-t1667-t16
     #71-t1586-t1589-t1591+t1595
      t1707 = t1598+t1601-t1606+t1607-t1610+t1611+t1615+t1617-t1619+t162
     #1-t1673-t1675+t1678+t1680
      spf349 = DENSPENCE(DCMPLX(t1662-(t1706+t1707)*t1218*t1462/2),spcZi
     #p)
      t1715 = t1229*t1419
      t1718 = 2*s*t293*t1255*t1715*t7
      t1721 = 2*t1503*t7*t1569
      t1722 = t1471+t1475+t1479-t1481+t1718+t1721-t1483-t1485+t1487-t149
     #0-t1493-t1496+t1498-t1501
      t1723 = 6*t1507
      t1725 = 2*t1378*t1422
      t1727 = 2*t1503*t1676
      t1730 = 2*t1271*t1715*t103
      t1731 = t1723+t1725-t1727-t1510-t1512-t1517-t1519+t1521+t1524+t152
     #7-t1730-t1529-t1531-t1533
      spf350 = log(DCMPLX(-t1461+(t1722+t1731)*t1218*t1462))
      t1736 = t1471-t1475-t1479-t1481+t1718-t1721+t1483-t1485+t1487-t149
     #0-t1493+t1496+t1498+t1501
      t1737 = -t1723+t1725-t1727-t1510-t1512+t1517+t1519+t1521+t1524+t15
     #27-t1730-t1529-t1531-t1533
      spf351 = log(DCMPLX(-t1455+(t1736+t1737)*t1218*t1456))
      t1742 = -t1471+t1475+t1479-t1481-t1718+t1721-t1483-t1485+t1487-t14
     #90+t1493-t1496+t1498-t1501
      t1743 = t1723-t1725+t1727-t1510-t1512-t1517-t1519+t1521+t1524+t152
     #7+t1730-t1529-t1531-t1533
      spf352 = log(DCMPLX(t1455-(t1742+t1743)*t1218*t1456))
      t1748 = -t1471-t1475-t1479-t1481-t1718-t1721+t1483-t1485+t1487-t14
     #90+t1493+t1496+t1498+t1501
      t1749 = -t1723-t1725+t1727-t1510-t1512+t1517+t1519+t1521+t1524+t15
     #27+t1730-t1529-t1531-t1533
      spf353 = log(DCMPLX(t1461-(t1748+t1749)*t1218*t1462))
      t1754 = t54**2
      t1755 = t1754*t1070
      t1756 = t1755*t1181
      t1757 = t842**2
      t1758 = t1756*t1757
      t1759 = t1754*s
      t1760 = t1759*t1181
      t1762 = 24*t1760*t449
      t1763 = t449*t103
      t1764 = t1763*t1759
      t1766 = 32*t1764*CsTh
      t1767 = t1754*t54
      t1768 = t449*t1767
      t1770 = 48*t1768*CsTh
      t1772 = 6*t1756*t842
      t1773 = t842*CsTh
      t1775 = 4*t1756*t1773
      t1778 = sqrt(s*t1227*t1183)
      t1779 = t1255*t1778
      t1780 = t1779*t1763
      t1781 = cw2*t1181
      t1782 = t1077*t1781
      t1784 = 32*t1780*t1782
      t1785 = t1755*t1255
      t1786 = t1785*t1778
      t1787 = t1779*t1070
      t1788 = t1181*t103
      t1791 = 122*t1787*t1788*t178
      t1792 = t178*t1255
      t1793 = t1778*t54
      t1794 = t449*t1181
      t1797 = 116*t1792*t1793*t1794
      t1798 = t177*t1181
      t1799 = t1072*t1798
      t1801 = 24*t1780*t1799
      t1802 = t1779*t1754
      t1803 = t449*CsTh
      t1804 = t7*t1803
      t1806 = 92*t1802*t1804
      t1807 = t1778*t1070
      t1810 = 220*t1792*t1807*t1803
      t1811 = t62*Mzsq
      t1812 = t177*cw2
      t1813 = t1811*t1812
      t1814 = t1813*t1255
      t1817 = 160*t1814*t1793*t1788
      t1820 = 68*t1787*t1794*t7
      t1821 = t1181*Mzsq
      t1822 = cw2*t103
      t1825 = 44*t1802*t1821*t1822
      t1826 = t1758+t1762+t1766+t1770+t1772+t1775+t1784-t1786-t1791-t179
     #7-t1801+t1806-t1810+t1817+t1820+t1825
      t1827 = t1181*t62
      t1828 = t177*CsTh
      t1829 = t1827*t1828
      t1832 = t62**2
      t1833 = t177**2
      t1834 = t1832*t1833
      t1835 = t1834*t1255
      t1836 = t1181*t842
      t1840 = t1181*CsTh
      t1844 = t1832*Mzsq
      t1845 = t1833*cw2
      t1846 = t1844*t1845
      t1847 = t1846*t1255
      t1848 = t1778*s
      t1849 = t1848*t1840
      t1852 = t103*t62
      t1853 = t1852*t1828
      t1859 = t1778*t449
      t1860 = s*t1181
      t1864 = t1779*t1759
      t1865 = t103*t842
      t1866 = t7*t1865
      t1871 = cw2*t1773
      t1878 = t1821*t1603
      t1881 = t177*t1773
      t1885 = t177*t842
      t1889 = t1848*t1788
      t1892 = t1852*t1885
      t1895 = -95*t1802*t1829-56*t1835*t1793*t1836-128*t1835*t1793*t1840
     #+32*t1847*t1849-372*t1802*t1853+176*t1814*t1793*t1803+64*t1814*t18
     #59*t1860+56*t1864*t1866+112*t1864*t1225+9*t1864*t1821*t1871+27*t18
     #64*t1821*t1323+27*t1864*t1878-25*t1802*t1827*t1881-85*t1802*t1827*
     #t1885-80*t1835*t1889-174*t1802*t1892
      t1899 = 196*t1787*t1813*t1865
      t1902 = 536*t1787*t1813*t1224
      t1905 = 32*t1835*t1848*t1803
      t1906 = t1779*t54
      t1909 = 48*t1906*t1834*t1865
      t1912 = 320*t1906*t1834*t1224
      t1916 = 20*t1787*t1813*t1181*t1773
      t1919 = 112*t1787*t1813*t1836
      t1920 = t1813*t1840
      t1922 = 164*t1787*t1920
      t1923 = t1778*t103
      t1926 = 64*t1847*t1923*t129
      t1927 = t1754*t1181
      t1928 = t1779*t1927
      t1930 = 44*t1928*t1866
      t1932 = t1785*t1778*t1773
      t1934 = 118*t1927*t1813
      t1935 = t1070*t1181
      t1937 = 168*t1935*t1834
      t1940 = 3*t1785*t1778*t842
      t1943 = 3*t1785*t1778*CsTh
      t1944 = t1767*t1181
      t1947 = 8*t1944*t103*t1773
      t1948 = t1899+t1902-t1905-t1909-t1912+t1916+t1919+t1922+t1926+t193
     #0-t1932-t1934+t1937-t1940-t1943+t1947
      t1950 = 24*t1944*t1865
      t1952 = 24*t1944*t1224
      t1953 = t449*t842
      t1955 = 24*t1760*t1953
      t1957 = 48*t1760*t1803
      t1958 = t1763*CsTh
      t1960 = 32*t1927*t1958
      t1961 = t1779*t1944
      t1962 = t1181*t54
      t1964 = 112*t1846*t1962
      t1966 = 10*t1944*t7
      t1968 = 46*t1760*t178
      t1969 = t1759*t1255
      t1971 = 12*t1969*t1859
      t1972 = t1767*t1255
      t1974 = 6*t1972*t1923
      t1975 = t1759*t62
      t1976 = t177*t103
      t1978 = 346*t1975*t1976
      t1979 = t1754*t62
      t1980 = t177*t449
      t1982 = 636*t1979*t1980
      t1983 = t1070*t103
      t1985 = 848*t1983*t1834
      t1986 = t1767*Mzsq
      t1988 = 82*t1986*t1822
      t1989 = t1759*Mzsq
      t1990 = cw2*t449
      t1992 = 204*t1989*t1990
      t1993 = t1950+t1952+t1955+t1957+t1960-t1961-t1964-t1966+t1968-t197
     #1-t1974+t1978+t1982+t1985-t1988-t1992
      t1996 = t1754*t1832
      t1999 = 1032*t1996*t1833*t842
      t2002 = 1016*t1996*t1833*CsTh
      t2003 = t1754*t1763
      t2005 = 216*t2003*t7
      t2006 = t7*t449
      t2008 = 92*t1802*t2006
      t2009 = t178*t1181
      t2011 = 35*t1802*t2009
      t2012 = t1759*t1811
      t2015 = 822*t2012*t1812*t842
      t2016 = t1812*CsTh
      t2018 = 662*t2012*t2016
      t2019 = t1935*t449
      t2020 = t178*t842
      t2022 = 300*t2019*t2020
      t2024 = 672*t2019*t1094
      t2025 = t1927*t62
      t2026 = t1976*t1773
      t2028 = 170*t2025*t2026
      t2029 = t1976*t842
      t2031 = 570*t2025*t2029
      t2032 = t1976*CsTh
      t2034 = 630*t2025*t2032
      t2035 = t1935*t1811
      t2036 = t1812*t103
      t2037 = t2036*t842
      t2039 = 676*t2035*t2037
      t2040 = t2036*CsTh
      t2042 = 944*t2035*t2040
      t2043 = t1935*t1763
      t2045 = 152*t2043*t1113
      t2046 = t1760*Mzsq
      t2047 = t1822*t1773
      t2049 = 66*t2046*t2047
      t2050 = t1999+t2002-t2005+t2008-t2011-t2015-t2018+t2022+t2024+t202
     #8+t2031+t2034-t2039-t2042-t2045-t2049
      t2051 = t1822*t842
      t2053 = 198*t2046*t2051
      t2054 = t1822*CsTh
      t2056 = 198*t2046*t2054
      t2057 = t1779*t1767
      t2058 = t7*t1773
      t2060 = 11*t2057*t2058
      t2062 = 33*t2057*t1324
      t2064 = 33*t2057*t1113
      t2065 = t178*t1773
      t2067 = 41*t1864*t2065
      t2069 = 133*t1864*t2020
      t2071 = 143*t1864*t1094
      t2072 = t1778*t1754
      t2075 = 56*t1814*t2072*t1773
      t2078 = 240*t1814*t2072*t842
      t2081 = 312*t1814*t2072*CsTh
      t2084 = 188*t1835*t1807*t842
      t2087 = 344*t1835*t1807*CsTh
      t2090 = 6*t1864*t1788*t842
      t2091 = t1788*CsTh
      t2093 = 12*t1864*t2091
      t2096 = 12*t1802*t1794*CsTh
      t2097 = -t2053-t2056+t2060+t2062+t2064-t2067-t2069-t2071+t2075+t20
     #78+t2081-t2084-t2087-t2090-t2093-t2096
      t2101 = 48*t1835*t1848*t449
      t2104 = 228*t1792*t1807*t449
      t2105 = t1070*cw2
      t2106 = t2105*Mzsq
      t2108 = 48*t1780*t2106
      t2109 = t54*t62
      t2110 = t2109*t177
      t2112 = 72*t1780*t2110
      t2115 = 208*t1814*t1793*t449
      t2116 = t1763*t1811
      t2117 = t2116*t1812
      t2118 = t1779*s
      t2120 = 16*t2117*t2118
      t2123 = 48*t1847*t1793*t842
      t2124 = t1793*CsTh
      t2126 = 160*t1847*t2124
      t2127 = t1832*t62
      t2128 = t1833*t177
      t2129 = t2127*t2128
      t2131 = t1848*CsTh
      t2133 = 32*t2129*t1255*t2131
      t2134 = t1763*t62
      t2136 = t1962*CsTh
      t2138 = 184*t2134*t177*t2136
      t2139 = t1927*Mzsq
      t2140 = t1990*t842
      t2142 = 156*t2139*t2140
      t2143 = t1990*CsTh
      t2145 = 312*t2139*t2143
      t2146 = t1813*t1181
      t2147 = t54*t449
      t2150 = 160*t2146*t2147*t842
      t2153 = 520*t2146*t2147*CsTh
      t2154 = t103*t1832
      t2155 = t2154*t1833
      t2157 = 576*t2155*t2136
      t2160 = 16*t1835*t1807*t1773
      t2161 = -t2101-t2104+t2108-t2112+t2115+t2120+t2123+t2126-t2133+t21
     #38-t2142-t2145-t2150-t2153+t2157-t2160
      t2162 = t2036*t1773
      t2164 = 136*t2035*t2162
      t2165 = t1962*t842
      t2167 = 272*t2155*t2165
      t2168 = t1860*CsTh
      t2170 = 32*t2117*t2168
      t2171 = t103*t1844
      t2174 = 32*t2171*t1845*t2168
      t2175 = t449*t1832
      t2178 = 64*t2175*t1833*t2168
      t2179 = t1813*t103
      t2181 = 340*t1787*t2179
      t2183 = 272*t1906*t2155
      t2185 = 72*t1787*t2146
      t2188 = 64*t1847*t1923*s
      t2191 = 72*t1835*t1793*t1181
      t2194 = 392*t1996*t1833*t1773
      t2195 = t1848*t1181
      t2197 = 32*t1847*t2195
      t2198 = t1976*t62
      t2200 = 198*t1802*t2198
      t2202 = 56*t1864*t446
      t2203 = t1821*cw2
      t2205 = 9*t1864*t2203
      t2206 = t1767*t62
      t2208 = 342*t2206*t1885
      t2210 = 242*t2206*t1828
      t2211 = -t2164+t2167-t2170-t2174+t2178+t2181-t2183+t2185+t2188-t21
     #91+t2194+t2197-t2200+t2202+t2205+t2208+t2210
      t2226 = t1755*Mzsq
      t2239 = t1070*t1844
      t2249 = t54*t2127
      t2258 = t1070*t449
      t2261 = -80*t2012*t1812*t1757-434*t2012*t1812*t1773+214*t2206*t188
     #1+50*t2206*t177*t1757-48*t2226*t1871-72*t2226*t1323-48*t2226*t1603
     #-12*t2226*cw2*t1757+32*t1996*t1833*t1757-128*t2239*t1845*t1773-592
     #*t2239*t1845*t842-800*t2239*t1845*CsTh+128*t2249*t2128*t842+256*t2
     #249*t2128*CsTh+448*t2147*t1834-888*t2258*t1813
      t2262 = t1754*t1811
      t2264 = 764*t2262*t2036
      t2265 = t449**2
      t2268 = 80*t7*t2265*t1070
      t2270 = 352*t1249*t1846
      t2273 = 80*t2265*t54*t178
      t2274 = t177*t1070
      t2276 = 440*t2134*t2274
      t2277 = t54*t1812
      t2279 = 272*t2116*t2277
      t2283 = 8*t1754*t1255*t1778*t1763
      t2285 = 32*t1927*t1763
      t2286 = t1779*t1935
      t2288 = 220*t2286*t1853
      t2292 = 108*t178*t1779*t2147*t1840
      t2294 = 68*t2286*t1804
      t2296 = 98*t2286*t1892
      t2297 = t1813*t1779
      t2298 = t1249*t1840
      t2300 = 224*t2297*t2298
      t2302 = 88*t1928*t1225
      t2305 = 48*t2297*t1247*t1840
      t2308 = 64*t2297*t1249*t1836
      t2309 = -t2264-t2268-t2270+t2273+t2276-t2279-t2283+t2285-t2288-t22
     #92+t2294-t2296+t2300+t2302+t2305+t2308
      t2314 = 80*t1834*t1779*t448*t1840
      t2315 = t1754**2
      t2317 = 194*t2012*t1812
      t2319 = 16*t1935*t2265
      t2320 = t1755*t103
      t2322 = 8*t2320*t1773
      t2324 = 336*t2239*t1845
      t2326 = 128*t2249*t2128
      t2328 = 12*t2226*cw2
      t2330 = 64*t2206*t177
      t2332 = 344*t1996*t1833
      t2334 = t1779*t1944*t1773
      t2335 = t2315*t1757
      t2337 = 4*t2315*CsTh
      t2339 = 6*t2315*t842
      t2341 = 4*t2315*t1773
      t2342 = -t2314+t2315-t2317+t2319+t2322-t2324+t2326-t2328+t2330+t23
     #32-t2334+t2335+t1756+t2337+t2339+t2341
      t2344 = 16*t2265*t1754
      t2345 = 32*t1764
      t2346 = 24*t1768
      t2347 = 8*t2320
      t2358 = (-64*t1846*t1181+32*t1834*t54+32*t1834*t1860-128*t1846*s+6
     #4*t2129)*t442
      t2360 = 96*t1834*t1247
      t2362 = 40*t1835*t1849
      t2364 = 40*t1906*t1920
      t2366 = 8*t1787*t1829
      t2367 = t1779*t103
      t2370 = 16*t2367*t2109*t1798
      t2372 = 64*t1814*t1889
      t2375 = 160*t1834*t1181*t801
      t2376 = t178*t1070
      t2378 = 32*t2376*t2091
      t2380 = 8*t1787*t2009
      t2382 = 40*t1906*t2146
      t2384 = 144*t2155*t2118
      t2387 = 96*t1814*t1923*t54
      t2389 = t1181*t1255*t1778
      t2391 = 48*t2155*t2389
      t2393 = 16*t2367*t2376
      t2395 = 56*t1835*t2195
      t2397 = 104*t1835*t2124
      t2398 = t2360-t2362+t2364-t2366-t2370+t2372+t2375+t2378-t2380+t238
     #2-t2384+t2387-t2391-t2393-t2395-t2397
      t2401 = 144*t1813*t54*t2091
      t2403 = 8*t1802*t1094
      t2404 = t1813*CsTh
      t2406 = 56*t1787*t2404
      t2408 = 32*t1847*t2131
      t2410 = 52*t1813*t1935
      t2412 = 144*t1834*t1962
      t2414 = 8*t178*t1927
      t2416 = 160*t1846*t1860
      t2419 = 8*t178*t1759*t842
      t2422 = 16*t178*t1759*CsTh
      t2424 = 32*t2129*t1779
      t2425 = t1754*t103
      t2427 = 32*t178*t2425
      t2429 = 384*t1834*t1249
      t2431 = 208*t1813*t1983
      t2433 = 192*t1846*t448
      t2435 = 96*t1846*t843
      t2437 = 480*t1846*t838
      t2438 = -t2401-t2403+t2406+t2408-t2410+t2412+t2414-t2416+t2419+t24
     #22-t2424+t2427+t2429-t2431-t2433-t2435-t2437
      t2441 = 192*t2129*t129
      t2443 = 168*t1834*t1119
      t2445 = 400*t1834*t1071
      t2446 = t1754*t842
      t2448 = 68*t1813*t2446
      t2449 = t1754*CsTh
      t2451 = 136*t1813*t2449
      t2453 = 32*t178*t2258
      t2455 = 144*t1813*t2147
      t2456 = t1247*t1181
      t2458 = 80*t1813*t2456
      t2460 = 32*t1846*t2389
      t2462 = 120*t1834*t1906
      t2463 = t1935*t842
      t2465 = 52*t1813*t2463
      t2467 = 96*t1846*t2118
      t2468 = t1935*CsTh
      t2470 = 104*t1813*t2468
      t2472 = 224*t1834*t2136
      t2473 = t178*t1759
      t2474 = 8*t2473
      t2476 = 68*t1813*t1754
      t2478 = 384*t1846*t54
      t2479 = t2441+t2443+t2445-t2448-t2451+t2453-t2455-t2458+t2460-t246
     #2-t2465+t2467-t2470+t2472+t2474-t2476-t2478
      t2481 = 192*t2129*s
      t2483 = 232*t1834*t1070
      t2486 = 32*t178*t1983*t1181
      t2488 = 160*t1846*t2168
      t2490 = 80*t1834*t2165
      t2493 = 16*t178*t1927*CsTh
      t2496 = 8*t178*t1927*t842
      t2500 = 32*t2171*t1845*t1255*t1778
      t2501 = t1070*t1811
      t2504 = 56*t1779*t2501*t1812
      t2507 = 8*t1779*t1979*t177
      t2508 = t1249*t1181
      t2510 = 144*t1813*t2508
      t2513 = 32*t178*t2425*CsTh
      t2515 = 192*t1846*t801
      t2516 = t1983*CsTh
      t2518 = 208*t1813*t2516
      t2520 = 384*t1834*t1251
      t2521 = t1788*s
      t2523 = 160*t1834*t2521
      t2524 = t2147*t1181
      t2526 = 32*t178*t2524
      t2527 = t2481+t2483+t2486-t2488+t2490+t2493+t2496+t2500+t2504-t250
     #7-t2510+t2513-t2515-t2518+t2520+t2523+t2526
      t2532 = 918*t1975*t2029
      t2534 = 978*t1975*t2032
      t2537 = 1200*t1979*t1980*CsTh
      t2540 = 424*t2134*t2274*CsTh
      t2543 = 224*t2116*t2277*CsTh
      t2544 = t1834*t1773
      t2546 = 128*t1983*t2544
      t2547 = t1834*t842
      t2549 = 1072*t1983*t2547
      t2550 = t1834*CsTh
      t2552 = 1792*t1983*t2550
      t2554 = 408*t1989*t2143
      t2556 = 204*t1989*t2140
      t2557 = t2344+t2345+t2346+t2347+t2358+(t2398+t2438+t2479+t2527)*Ie
     #ps+t2532+t2534+t2537+t2540-t2543+t2546+t2549+t2552-t2554-t2556
      t2586 = t1813*t842
      t2596 = -216*t2003*t1113-82*t1986*t2047-246*t1986*t2051-246*t1986*
     #t2054+408*t1935*t2547+456*t1935*t2550+564*t1979*t1980*t842+286*t19
     #75*t2026+142*t1760*t2065+234*t1760*t2020+170*t1760*t1094-214*t1927
     #*t1813*t1773-450*t1927*t2586-386*t1927*t2404+32*t1760*t178*t1757-4
     #0*t1944*t1113
      t2632 = -60*t1944*t1324-112*t1846*t2165-224*t1846*t2136-10*t1944*t
     #7*t1757-40*t1944*t2058-376*t2262*t2162-1516*t2262*t2037-1904*t2262
     #*t2040-1432*t2258*t2404-544*t2258*t2586+24*t2320*t842+576*t2147*t2
     #550+128*t2147*t2547-32*t1927*t1813*t1757+120*t1935*t2544-256*t1249
     #*t1846*t842
      t2633 = t2596+t2632
      t2636 = 608*t1249*t1846*CsTh
      t2639 = 3*t1779*t1944*t842
      t2642 = 3*t1779*t1944*CsTh
      t2646 = 64*t2175*t1833*t1181*s
      t2649 = 12*t1969*t1859*CsTh
      t2652 = 6*t1972*t1923*t842
      t2655 = 12*t1972*t1923*CsTh
      t2659 = 304*t2154*t1833*t54*t1181
      t2663 = 32*t2171*t1845*t1181*s
      t2665 = 360*t1813*t2524
      t2667 = 156*t1927*t2006
      t2670 = 6*t1779*t1760*t103
      t2673 = 12*t1779*t1927*t449
      t2675 = 112*t1846*t1906
      t2677 = 32*t2129*t2118
      t2679 = 172*t1834*t1787
      t2680 = -t2636-t2639-t2642+t2646-t2649-t2652-t2655+t2659-t2663-t26
     #65-t2667-t2670-t2673+t2675-t2677-t2679
      t2682 = 128*t1813*t1802
      t2684 = 51*t1779*t2473
      t2687 = 11*t1779*t1986*cw2
      t2689 = 66*t1760*t446
      t2693 = 152*t1935*t1763*Mzsq*cw2
      t2695 = 404*t1935*t2179
      t2697 = 230*t1927*t2198
      t2699 = 24*t2320*CsTh
      t2701 = 24*t1768*t842
      t2703 = 4*t1756*CsTh
      t2705 = 8*t1944*t103
      t2709 = 372*t1935*t449*t62*t177
      t2710 = t177*t54
      t2713 = 200*t2134*t2710*t1181
      t2717 = 48*t7*t1181*t2265*t54
      t2721 = 16*t2265*t62*t1798*s
      t2725 = 48*t2116*t1812*t1181*s
      t2727 = 8*t1779*t2043
      t2728 = t2682-t2684+t2687-t2689-t2693-t2695+t2697+t2699+t2701+t270
     #3+t2705+t2709+t2713-t2717+t2721-t2725-t2727
      t2733 = 1/t1778
      t2735 = t2118*t1878
      t2737 = 2*t2367*t2203
      t2739 = 6*t2367*t176
      t2741 = 3*t1906*t1113
      t2742 = t1781*CsTh
      t2744 = 8*t1318*t2742
      t2745 = t2118*t2203
      t2746 = t1779*t1962
      t2747 = t1779*t1071
      t2749 = 24*t838*t178
      t2751 = 6*t1119*t7
      t2753 = 12*t2105*t720
      t2754 = t2735+t2737+t2739+t2741-t2744+t2745-t2746-t2747+t2749-t275
     #1-t2753
      t2756 = 16*t1077*t1822
      t2758 = 8*t2109*t1885
      t2759 = s*t1811
      t2761 = 16*t2759*t2016
      t2762 = 4*t2298
      t2764 = 16*t448*t178
      t2766 = 8*t1247*t7
      t2768 = 2*t1779*t1249
      t2769 = 4*t1782
      t2770 = 8*t1799
      t2771 = -t2756+t2758-t2761+t2762+t2764-t2766-t2768-t2769+t2770+t17
     #54+t2463
      t2773 = 2*t2468
      t2774 = 4*t2508
      t2775 = 4*t2516
      t2777 = 16*t2759*t1812
      t2778 = 16*t2110
      t2779 = 6*t2106
      t2780 = 4*t2456
      t2782 = 2*t1779*t2521
      t2784 = 8*t1077*t2742
      t2787 = 4*t1077*t1781*t842
      t2788 = t2773+t2774+t2775-t2777+t2778-t2779-t1787+t2780-t2782-t278
     #4-t2787
      t2790 = 8*t448*t2203
      t2792 = 3*t1779*t1078
      t2794 = 16*t1077*t2054
      t2796 = 16*t448*t1094
      t2799 = 8*t1072*t1798*CsTh
      t2800 = t1779*t2136
      t2801 = 2*t2449
      t2802 = 4*t2147
      t2803 = 4*t1983
      t2807 = (-16*t1813+t1074+8*t2009)*Ieps
      t2808 = -t2790+t2792-t2794+t2796+t2799-t2800+t2446+t1935+t2801+t28
     #02+t2803+t2807
      t2814 = t1385*t1387*t1183
      spf354 = log(DCMPLX((t1826+t1895+t1948+t1993+t2050+t2097+t2161+t22
     #11+t2261+t2309+t2342+t2557+t2633+t2680+t2728)*t2733/(t2754+t2771+t
     #2788+t2808)*t2814))
      t2816 = t1758+t1762+t1766+t1770+t1772+t1775-t1784+t1786+t1791+t179
     #7+t1801-t1806+t1810-t1817-t1820-t1825
      t2819 = -t1899-t1902+t1905+t1909+t1912-t1916-t1919-t1922-t1926-t19
     #30+t1932-t1934+t1937+t1940+t1943+t1947
      t2820 = t1950+t1952+t1955+t1957+t1960+t1961-t1964-t1966+t1968+t197
     #1+t1974+t1978+t1982+t1985-t1988-t1992
      t2823 = t1999+t2002-t2005-t2008+t2011-t2015-t2018+t2022+t2024+t202
     #8+t2031+t2034-t2039-t2042-t2045-t2049
      t2824 = -t2053-t2056-t2060-t2062-t2064+t2067+t2069+t2071-t2075-t20
     #78-t2081+t2084+t2087+t2090+t2093+t2096
      t2826 = t2101+t2104-t2108+t2112-t2115-t2120-t2123-t2126+t2133+t213
     #8-t2142-t2145-t2150-t2153+t2157+t2160
      t2827 = -t2164+t2167-t2170-t2174+t2178-t2181+t2183-t2185-t2188+t21
     #91+t2194-t2197+t2200-t2202-t2205+t2208+t2210
      t2831 = -t2264-t2268-t2270+t2273+t2276-t2279+t2283+t2285+t2288+t22
     #92-t2294+t2296-t2300-t2302-t2305-t2308
      t2833 = t2314+t2315-t2317+t2319+t2322-t2324+t2326-t2328+t2330+t233
     #2+t2334+t2335+t1756+t2337+t2339+t2341
      t2834 = t2360+t2362-t2364+t2366+t2370-t2372+t2375+t2378+t2380-t238
     #2+t2384-t2387+t2391+t2393+t2395+t2397
      t2835 = -t2401+t2403-t2406-t2408-t2410+t2412+t2414-t2416+t2419+t24
     #22+t2424+t2427+t2429-t2431-t2433-t2435-t2437
      t2837 = t2441+t2443+t2445-t2448-t2451+t2453-t2455-t2458-t2460+t246
     #2-t2465-t2467-t2470+t2472+t2474-t2476-t2478
      t2838 = t2481+t2483+t2486-t2488+t2490+t2493+t2496-t2500-t2504+t250
     #7-t2510+t2513-t2515-t2518+t2520+t2523+t2526
      t2842 = t2344+t2345+t2346+t2347+t2358+(t2834+t2835+t2837+t2838)*Ie
     #ps+t2532+t2534+t2537+t2540-t2543+t2546+t2549+t2552-t2554-t2556
      t2845 = -t2636+t2639+t2642+t2646+t2649+t2652+t2655+t2659-t2663-t26
     #65-t2667+t2670+t2673-t2675+t2677+t2679
      t2846 = -t2682+t2684-t2687-t2689-t2693-t2695+t2697+t2699+t2701+t27
     #03+t2705+t2709+t2713-t2717+t2721-t2725+t2727
      t2852 = -t2735-t2737-t2739-t2741-t2744-t2745+t2746+t2747+t2749-t27
     #51-t2753
      t2853 = -t2756+t2758-t2761+t2762+t2764-t2766+t2768-t2769+t2770+t17
     #54+t2463
      t2855 = t2773+t2774+t2775-t2777+t2778-t2779+t1787+t2780+t2782-t278
     #4-t2787
      t2856 = -t2790-t2792-t2794+t2796+t2799+t2800+t2446+t1935+t2801+t28
     #02+t2803+t2807
      spf355 = log(DCMPLX(-(t2816-t1895+t2819+t2820+t2823+t2824+t2826+t2
     #827+t2261+t2831+t2833+t2842+t2633+t2845+t2846)*t2733/(t2852+t2853+
     #t2855+t2856)*t2814))
      t2864 = t1419/spf28
      t2865 = t2864*t103
      t2866 = t2864*t7
      t2870 = log(DCMPLX(t7/(t7+t103+t2865-t2866)))
      t2873 = log(DCMPLX(t7*t1232))
      t2878 = log(DCMPLX(-t7/(-t7-t103+t2865-t2866)))
      t2882 = log(DCMPLX(t1231*t145))
      t2884 = -spf349-spf348+spf347+spf346+spf345+spf344-spf343-spf342-s
     #pf321+spf320-spf319+spf318+spf352*spf334-t2870*spf352+t2873*spf314
     #+t2878*spf353-spf353*spf335-t2882*spf322
      t2922 = -2*spf315*spf322-2*t124*spf322-5*t124*spf331+2*t2870*spf35
     #1-2*spf351*spf337+2*spf351*spf327+3*t124*spf335-2*spf327*spf354+2*
     #spf323*spf354+2*t1185*spf354-2*spf317*spf328+3*t124*spf328-2*spf33
     #8*spf332+2*t2878*spf338+2*spf338*spf322+3*t124*spf337-2*spf352*spf
     #326-2*spf314*spf325
      t2960 = -2*spf339*spf323+2*spf339*spf333-2*t2870*spf339+2*spf340*s
     #pf322-2*spf322*spf355+2*spf326*spf355-2*t1185*spf355-2*t2873*spf31
     #5-5*t124*spf333+5*t124*spf332-3*t124*spf336-2*spf340*spf330+2*spf3
     #16*spf329-3*t124*spf329+5*t124*spf325+5*t124*spf330+2*spf353*spf32
     #7-2*spf350*spf326
      t2997 = 2*spf350*spf336-2*t2878*spf350+2*t2882*spf323+2*t124*spf32
     #3+2*spf314*spf323+2*t2870*spf340+2*t2873*spf317+2*spf315*spf324-5*
     #t124*spf324-2*t2882*spf327-2*spf316*spf327-2*t2873*spf316+2*spf317
     #*spf326+2*t2882*spf326+2*spf341*spf331-2*t2878*spf341-3*t124*spf33
     #4-2*spf341*spf323
      wV17 = -(2*t2884+t2922+t2960+t2997)*t1255*t2733/(s+t1181)/188
      t3006 = 2*t910
      t3007 = 4*t62
      t3010 = 4*t57
      t3011 = 4*t910
      t3012 = 8*t62
      t3013 = t54+t838-t3010-t3011+t3012
      t3014 = sqrt(t3013)
      t3015 = 1/t3014
      t3018 = t167*t903*Ieps
      t3019 = t1050*(t3006-t3007+t58)*t3015*t3018
      t3020 = t1077*CsTh
      t3021 = 4*t3020
      t3022 = 8*t1072
      t3023 = t1075*t3014
      t3024 = 4*t1077
      t3028 = t3015*t167*t903
      spf252 = DENSPENCE(DCMPLX(-t3019/2-t1050*(t1071-t3021+t3022-t3023-
     #t3024+t1070)*t3028/2),spcZip)
      spf253 = DENSPENCE(DCMPLX(t3019/2+t1050*(t1071-t3021+t3022+t3023-t
     #3024+t1070)*t3028/2),spcZip)
      t3037 = 4*t1093
      t3039 = 2*t57*t842
      t3042 = (-t3007-t3037+t58+t3011+t3039)*t1100*t3015
      t3043 = t3042*t3018
      t3044 = t1050*t3014
      t3046 = 4*t3044*t62
      t3047 = t3014*Mzsq
      t3049 = 2*t1107*t3047
      t3052 = 2*t1107*t3047*CsTh
      t3054 = t1075*t3014*CsTh
      t3055 = 8*t3020
      t3056 = t1077*t842
      t3057 = 4*t3056
      t3058 = t1072*CsTh
      t3059 = 8*t3058
      t3060 = -t3046+t3049+t3052-t3023-t3054+t1070+t1118+t1119-t3024-t30
     #55-t3057+t3022+t3059
      t3061 = t3060*t1100
      spf259 = log(DCMPLX(t3043+t3061*t3028))
      t3064 = t3046-t3049-t3052+t3023+t3054+t1070+t1118+t1119-t3024-t305
     #5-t3057+t3022+t3059
      t3065 = t3064*t1100
      spf261 = log(DCMPLX(-t3043-t3065*t3028))
      t3070 = t3042*t903*t513*Ieps
      t3072 = Mzsq*t10*t506
      t3074 = 2*t3044*t3072
      t3075 = t3014*t506
      t3076 = t1144*t3075
      t3078 = t1144*t3075*CsTh
      t3079 = t3049-t3074+t3076+t3078+t1070+t1118+t1119-t3024-t3055-t305
     #7+t3022+t3059
      t3081 = t3015*t903
      t3082 = t3081*t513
      t3084 = t3070+t3079*t1100*t3082
      spf264 = log(DCMPLX(t3084))
      t3087 = t3042*t903*t509*Ieps
      t3088 = t3049+t3074-t3076-t3078+t1070+t1118+t1119-t3024-t3055-t305
     #7+t3022+t3059
      t3090 = t3081*t509
      t3092 = t3087+t3088*t1100*t3090
      spf265 = log(DCMPLX(t3092))
      spf268 = DENSPENCE(DCMPLX(t3084),spcZip)
      spf269 = DENSPENCE(DCMPLX(t3092),spcZip)
      t3093 = t1070+t1118+t1119-t3024-t3055-t3057+t3022+t3059+t3078-t304
     #9-t3074+t3076
      t3096 = -t3087-t3093*t1100*t3090
      spf270 = log(DCMPLX(t3096))
      t3097 = t1070+t1118+t1119-t3024-t3055-t3057+t3022+t3059-t3078-t304
     #9+t3074-t3076
      t3100 = -t3070-t3097*t1100*t3082
      spf271 = log(DCMPLX(t3100))
      spf274 = DENSPENCE(DCMPLX(t3096),spcZip)
      spf275 = DENSPENCE(DCMPLX(t3100),spcZip)
      spf280 = log(DCMPLX(t3087+t3061*t3090))
      spf281 = log(DCMPLX(t3070+t3061*t3082))
      spf284 = log(DCMPLX(-t3087-t3065*t3090))
      spf285 = log(DCMPLX(-t3070-t3065*t3082))
      t3111 = sqrt(-s*t520)
      t3115 = log(DCMPLX(Mzsq/(t504-s-t3111)))
      t3121 = log(DCMPLX(Mzsq/(t504-s+t3111)))
      t3124 = log(DCMPLX(t789))
      t3131 = spf270*spf259+spf270*t3115+spf271*spf259+spf271*t3121-spf2
     #59*t124-spf259*t3124-spf265*spf261-spf265*t3115-spf264*spf261-spf2
     #64*t3121+spf261*t124
      t3139 = spf261*t3124-spf270*spf280-spf271*spf281+spf265*spf284+spf
     #264*spf285+spf268+2*spf252-spf274-spf275-2*spf253+spf269
      t3143 = sqrt(t1211*t3013)
      wV18 = (t3131+t3139)/t3143/2
      t3148 = sqrt(-t966)
      t3150 = t54-t838-t1056+t1057+t1058
      t3151 = sqrt(t3150)
      t3152 = 1/t3151
      t3155 = 1/(t291-s+t129)
      t3157 = t169*t3155*Ieps
      t3158 = (t1052+t1051-t1053)*t3148*t3152*t3157
      t3159 = t54*t3148
      t3160 = t3159*t3151
      t3164 = t169*t3155
      spf254 = DENSPENCE(DCMPLX(t3158/2+(-t1069+t1071-t1074-t1070-t3160+
     #t1079)*t3148*t3152*t3164/2),spcZip)
      spf255 = DENSPENCE(DCMPLX(-t3158/2-(-t1069+t1071-t1074-t1070+t3160
     #+t1079)*t3148*t3152*t3164/2),spcZip)
      t3175 = 1/t3148
      t3177 = (-t1052+t1095+t1053-t1057+t1098)*t3175*t3152
      t3178 = t3177*t3157
      t3179 = t3148*t3151
      t3181 = 4*t3179*t178
      t3182 = s*t3148
      t3186 = 2*t3182*t3151*Mzsq*cw2
      t3189 = 2*t3182*t3151*t1113
      t3191 = t3159*t3151*CsTh
      t3192 = t3181-t3186+t3189+t3160-t3191+t1070-t1118+t1119-t1079+t112
     #0-t1122+t1074-t1124
      t3193 = t3192*t3175
      spf286 = log(DCMPLX(t3178+t3193*t3152*t3164))
      t3197 = -t3181+t3186-t3189-t3160+t3191+t1070-t1118+t1119-t1079+t11
     #20-t1122+t1074-t1124
      t3198 = t3197*t3175
      spf288 = log(DCMPLX(-t3178-t3198*t3152*t3164))
      t3204 = t3177*t3155*t296*Ieps
      t3205 = t1143*t3148
      t3206 = t3151*t293
      t3208 = t3205*t3206*CsTh
      t3211 = 2*t3179*Mzsq*t1140
      t3212 = t3205*t3206
      t3213 = t1070-t1118+t1119-t1079+t1120-t1122+t1074-t1124-t3208-t318
     #6-t3211+t3212
      t3215 = t3152*t3155
      t3216 = t3215*t296
      t3218 = t3204+t3213*t3175*t3216
      spf290 = log(DCMPLX(t3218))
      t3221 = t3177*t3155*t300*Ieps
      t3222 = t1070-t1118+t1119-t1079+t1120-t1122+t1074-t1124+t3208-t318
     #6+t3211-t3212
      t3224 = t3215*t300
      t3226 = t3221+t3222*t3175*t3224
      spf291 = log(DCMPLX(t3226))
      spf294 = DENSPENCE(DCMPLX(t3218),spcZip)
      spf295 = DENSPENCE(DCMPLX(t3226),spcZip)
      t3227 = t3186-t3211+t3212-t3208+t1070-t1118+t1119-t1079+t1120-t112
     #2+t1074-t1124
      t3230 = -t3221-t3227*t3175*t3224
      spf300 = log(DCMPLX(t3230))
      t3231 = t3186+t3211-t3212+t3208+t1070-t1118+t1119-t1079+t1120-t112
     #2+t1074-t1124
      t3234 = -t3204-t3231*t3175*t3216
      spf301 = log(DCMPLX(t3234))
      spf304 = DENSPENCE(DCMPLX(t3230),spcZip)
      spf305 = DENSPENCE(DCMPLX(t3234),spcZip)
      spf306 = log(DCMPLX(t3204+t3193*t3216))
      spf307 = log(DCMPLX(t3221+t3193*t3224))
      spf310 = log(DCMPLX(-t3204-t3198*t3216))
      spf311 = log(DCMPLX(-t3221-t3198*t3224))
      t3252 = log(DCMPLX(-t7*t970))
      t3256 = spf290*spf310+spf291*spf311-spf301*spf306-spf290*spf288-sp
     #f290*t1185-spf291*spf288-spf291*t1191+spf288*t124+spf288*t3252+spf
     #301*spf286+spf301*t1185
      t3264 = spf300*spf286+spf300*t1191-spf286*t124-spf286*t3252-2*spf2
     #55+spf294+spf295+2*spf254-spf305-spf304-spf300*spf307
      t3266 = t54*t966
      t3268 = sqrt(-t3266*t3150)
      wV19 = -(t3256+t3264)/t3268/2
      t3274 = t54-t838-t3010+t3011+t3012
      t3275 = sqrt(t3274)
      t3276 = 1/t3275
      t3279 = t167*t975*Ieps
      t3280 = (t3007+t3006-t58)*t3148*t3276*t3279
      t3281 = t3159*t3275
      t3285 = t3276*t167*t975
      spf250 = DENSPENCE(DCMPLX(t3280/2+(t1071-t3021-t3022-t1070-t3281+t
     #3024)*t3148*t3285/2),spcZip)
      spf251 = DENSPENCE(DCMPLX(-t3280/2-(t1071-t3021-t3022-t1070+t3281+
     #t3024)*t3148*t3285/2),spcZip)
      t3296 = (-t3007+t3037+t58-t3011+t3039)*t3175*t3276
      t3297 = t3296*t3279
      t3298 = t3148*t3275
      t3300 = 4*t3298*t62
      t3301 = t3275*Mzsq
      t3303 = 2*t3182*t3301
      t3306 = 2*t3182*t3301*CsTh
      t3308 = t3159*t3275*CsTh
      t3309 = t3300-t3303+t3306+t3281-t3308+t1070-t1118+t1119-t3024+t305
     #5-t3057+t3022-t3059
      t3310 = t3309*t3175
      spf258 = log(DCMPLX(t3297+t3310*t3285))
      t3313 = -t3300+t3303-t3306-t3281+t3308+t1070-t1118+t1119-t3024+t30
     #55-t3057+t3022-t3059
      t3314 = t3313*t3175
      spf260 = log(DCMPLX(-t3297-t3314*t3285))
      t3319 = t3296*t975*t509*Ieps
      t3320 = t3275*t506
      t3322 = t3205*t3320*CsTh
      t3324 = 2*t3298*t3072
      t3325 = t3205*t3320
      t3326 = t1070-t1118+t1119-t3024+t3055-t3057+t3022-t3059-t3322-t330
     #3-t3324+t3325
      t3328 = t3276*t975
      t3329 = t3328*t509
      t3331 = t3319+t3326*t3175*t3329
      spf262 = log(DCMPLX(t3331))
      t3334 = t3296*t975*t513*Ieps
      t3335 = t1070-t1118+t1119-t3024+t3055-t3057+t3022-t3059+t3322-t330
     #3+t3324-t3325
      t3337 = t3328*t513
      t3339 = t3334+t3335*t3175*t3337
      spf263 = log(DCMPLX(t3339))
      spf266 = DENSPENCE(DCMPLX(t3331),spcZip)
      spf267 = DENSPENCE(DCMPLX(t3339),spcZip)
      t3340 = t3303-t3324+t3325-t3322+t1070-t1118+t1119-t3024+t3055-t305
     #7+t3022-t3059
      t3343 = -t3334-t3340*t3175*t3337
      spf272 = log(DCMPLX(t3343))
      t3344 = t3303+t3324-t3325+t3322+t1070-t1118+t1119-t3024+t3055-t305
     #7+t3022-t3059
      t3347 = -t3319-t3344*t3175*t3329
      spf273 = log(DCMPLX(t3347))
      spf276 = DENSPENCE(DCMPLX(t3343),spcZip)
      spf277 = DENSPENCE(DCMPLX(t3347),spcZip)
      spf278 = log(DCMPLX(t3319+t3310*t3329))
      spf279 = log(DCMPLX(t3334+t3310*t3337))
      spf282 = log(DCMPLX(-t3319-t3314*t3329))
      spf283 = log(DCMPLX(-t3334-t3314*t3337))
      t3365 = log(DCMPLX(-t972))
      t3368 = -spf262*spf282-spf263*spf283+spf273*spf278+spf272*spf279+s
     #pf262*spf260+spf262*t3115+spf263*spf260+spf263*t3121-spf260*t124-s
     #pf260*t3365-spf273*spf258
      t3376 = -spf273*t3115-spf272*spf258-spf272*t3121+spf258*t124+spf25
     #8*t3365+2*spf251-spf266-spf267-2*spf250+spf277+spf276
      t3379 = sqrt(-t3266*t3274)
      wV20 = (t3368+t3376)/t3379/2
      t3382 = 5*t57
      t3386 = t62*cw2
      t3388 = s*t177
      t3389 = t3388*Mzsq
      t3392 = 14*t54
      t3393 = t62*t1812
      t3396 = t62*t1833
      t3399 = s*t1812*Mzsq
      t3402 = s*t1833*Mzsq
      t3404 = t3382+148*t2710-44*t1067+508*t178-40*t3386-244*t3389-68*t1
     #76+t3392-932*t3393-52*t2277+464*t3396+488*t3399-208*t3402
      t3406 = conjg(cw2)
      t3407 = t3406*s
      t3409 = conjg(Mzsq)
      t3410 = t3406*t3409
      t3412 = t3406**2
      t3413 = t3412*t3409
      t3416 = t966**2
      t3418 = cw2-1
      t3419 = t3418**2
      t3420 = 1/t3419
      t3424 = 1/(-s+t3409)
      t3425 = t1387*t3424
      t3426 = 1/t3406
      t3428 = 1/(-1+t3406)
      t3429 = t3426*t3428
      t3430 = t502**2
      t3431 = 1/t3430
      t3433 = t3425*t3429*t3431
      wV25 = -t211*t3404*(2*t3407-8*t3410+8*t3413+s)*t3416*t3420*t3433/2
     #7648
      t3436 = t152**2
      t3437 = 4*t3410
      t3440 = 37*t57
      t3441 = 52*t2710
      t3452 = 1/t3418
      t3453 = t3452*t1387
      t3454 = t3424*t3426
      t3456 = t3453*t3454*t3431
      wV26 = -t3436*(t3437-s)*(-t3440+t3441-24*t1067-612*t178+148*t3386+
     #192*t3389+104*t176-t54+464*t3393-208*t3399)*t211*t3456/13824
      t3459 = 2*t3410
      t3461 = t3436*(t3459-s)
      t3463 = 12*t1067
      t3467 = 114*t176
      wV27 = -t3461*(-25*t57+t3441-t3463-264*t178+50*t3386-4*t3389+t3467
     #-13*t54+232*t3393-104*t3399)*t211*t3456/1728
      t3480 = 38*t176
      t3481 = 26*t3386
      t3486 = (t3410-s)*t3416
      t3489 = t1387*t3431*t3454
      wV28 = -t211*(-26*t3389+58*t178-7*t54+26*t1067-t3480-t3481+26*t57)
     #*t3486*t3489/54
      t3493 = 20*t1067
      t3495 = 20*t3386
      t3502 = -t3382+140*t2710+t3493+284*t178-t3495-380*t3389+t3480-t339
     #2-472*t3393-80*t2277+208*t3396+400*t3399-80*t3402
      t3507 = (4*t3407-t3437+4*t3413-s)*t3416
      wV31 = t211*t3502*t3507*t3420*t3433/15552
      t3512 = t3436*(t3459+s)
      t3513 = 40*t2710
      t3518 = 208*t3393
      t3519 = 80*t3399
      wV32 = -t3512*(t3440-t3513+t3463-336*t178+74*t3386+288*t3389-86*t1
     #76+t54+t3518-t3519)*t211*t3456/6912
      wV33 = -t3461*(-13*t57+t3513+t3463-216*t178+t3481-40*t3389+t3467-2
     #5*t54+t3518-t3519)*t211*t3456/6912
      wV34 = -t211*(-20*t3389+52*t178-t54+t3493-t3480-t3495+20*t57)*t348
     #6*t3489/216
      t3542 = t1811*t177
      t3544 = 63*t1318
      t3546 = 63*t1249
      t3549 = t1852*cw2
      t3551 = t1811*cw2
      t3556 = t3386*s
      t3559 = t1811*t1833
      t3565 = t3393*s
      t3567 = t2710*Mzsq
      t3573 = -1136*t3542+t3544+144*t1240-t3546+36*t1067*t103-180*t3549+
     #80*t3551-180*t1976*t57+56*t1077-152*t3556+1888*t1813-832*t3559+180
     #*t2198+20*t1072-80*t1078+1520*t1073-1600*t3565-560*t3567+320*t2277
     #*Mzsq+320*t3396*s
      t3577 = t3420*t3424
      t3578 = t3577*t3426
      t3579 = t3428*t167
      wV35 = -t211*t3573*t3507*t1387*t3578*t3579*t3431/248832
      t3588 = 2*t1077
      t3596 = 672*t3542-t3544-90*t1240+t3546+90*t3549-148*t3551-t3588+17
     #2*t3556-416*t1813-74*t1072-24*t1078-576*t1073+160*t3565+80*t3567
      t3601 = t3426*t167
      wV36 = t3512*t3596*t211*t1387*t3452*t3424*t3601*t3431/27648
      t3605 = t211*wV10
      t3606 = t3420*t3428
      t3607 = t3605*t3606
      t3608 = t3406*t3424
      t3609 = Mzsq*t3409
      t3611 = t3607*t3608*t3609
      t3613 = t3420*cw2
      t3614 = t3605*t3613
      t3615 = t3428*t3424
      t3616 = t3615*t3609
      t3619 = t3420*t177
      t3620 = t3605*t3619
      t3621 = t3620*t3616
      t3623 = t3428*t3406
      t3625 = t3424*Mzsq*t3409
      t3626 = t3623*t3625
      t3627 = t3620*t3626
      t3631 = t3615*t838
      t3634 = t3420*t1387
      t3635 = t3605*t3634
      t3636 = t3635*t3616
      t3638 = s*t3409
      t3639 = t3638*CsTh
      t3640 = t3615*t3639
      t3641 = t3635*t3640
      t3643 = t3420*t168
      t3644 = t3605*t3643
      t3645 = t3424*s
      t3646 = t3645*Mzsq
      t3647 = t3429*t3646
      t3650 = t3634*t3428
      t3652 = t3608*t3639
      t3653 = t3605*t3650*t3652
      t3657 = t3424*t54
      t3658 = t3657*CsTh
      t3659 = t3429*t3658
      t3660 = t3635*t3659
      t3662 = t3620*t3659
      t3668 = t3454*t57
      t3669 = t3607*t3668
      t3671 = t3620*t3647
      t3673 = t3635*t3626
      t3677 = t3619*t3428
      t3679 = t3605*t3677*t3652
      t3681 = -13.D0/50.D0*t3611-12.D0/25.D0*t3614*t3616+8.D0/25.D0*t362
     #1-8.D0/25.D0*t3627+12.D0/25.D0*t3614*t3626-3.D0/50.D0*t3614*t3631+
     #t3636/200-t3641/400+3.D0/400.D0*t3644*t3647+t3653/400-3.D0/50.D0*t
     #3644*t3616+t3660/3200+t3662/50+3.D0/50.D0*t3614*t3647+3.D0/50.D0*t
     #3644*t3626-13.D0/400.D0*t3669-t3671/25-t3673/200-3.D0/100.D0*t3614
     #*t3659+4.D0/25.D0*t3679
      t3682 = t3613*t3428
      t3686 = t3643*t3428
      t3690 = t3607*t3652
      t3692 = t3409*CsTh
      t3693 = t3645*t3692
      t3694 = t3607*t3693
      t3696 = t3635*t3631
      t3700 = t3620*t3640
      t3706 = t3454*t838
      t3707 = t3607*t3706
      t3709 = t3615*t57
      t3710 = t3635*t3709
      t3712 = t3620*t3631
      t3716 = t3635*t3647
      t3720 = t3620*t3709
      t3722 = t3605*t3420
      t3723 = t3722*t3631
      t3727 = t3722*t3616
      t3729 = t3722*t3709
      t3731 = -6.D0/25.D0*t3605*t3682*t3652-3.D0/100.D0*t3605*t3686*t365
     #2+13.D0/100.D0*t3690-13.D0/100.D0*t3694+t3696/1600+3.D0/100.D0*t36
     #44*t3640-4.D0/25.D0*t3700+6.D0/25.D0*t3614*t3640-3.D0/800.D0*t3644
     #*t3659+13.D0/800.D0*t3707-t3710/800+t3712/25-3.D0/400.D0*t3644*t36
     #31-t3716/1600+3.D0/200.D0*t3644*t3709-2.D0/25.D0*t3720+13.D0/400.D
     #0*t3723+3.D0/25.D0*t3614*t3709+13.D0/50.D0*t3727-13.D0/200.D0*t372
     #9
      wI29pc1 = t3681+t3731
      t3732 = t211*wV14
      t3733 = t3732*t3615
      t3734 = t3613*t3639
      t3737 = t3638*t842
      t3738 = t3613*t3737
      t3742 = t3733*t3634*t843
      t3745 = t3733*t3634*t838
      t3747 = t3732*t3454
      t3748 = t3606*t843
      t3749 = t3747*t3748
      t3751 = t3428*t168
      t3752 = t3420*t54
      t3753 = t3752*CsTh
      t3757 = t3606*t838
      t3758 = t3747*t3757
      t3761 = t3733*t3634*t57
      t3763 = t3752*t842
      t3767 = t3420*CsTh
      t3769 = t3733*t3767*t57
      t3771 = t3643*t910
      t3774 = t3619*t910
      t3775 = t3733*t3774
      t3777 = t3613*t910
      t3780 = t3428*t177
      t3782 = t3747*t3780*t3753
      t3785 = t3747*t3780*t3763
      t3794 = t3732*t3608*t3428
      t3795 = t720*t3409
      t3796 = t3643*t3795
      t3804 = -3.D0/25.D0*t3733*t3734-3.D0/25.D0*t3733*t3738-t3742/3200-
     #t3745/3200-13.D0/1600.D0*t3749+3.D0/1600.D0*t3747*t3751*t3753-13.D
     #0/1600.D0*t3758+t3761/1600+3.D0/1600.D0*t3747*t3751*t3763+13.D0/40
     #0.D0*t3769-3.D0/400.D0*t3733*t3771+t3775/25-3.D0/50.D0*t3733*t3777
     #-t3782/100-t3785/100+3.D0/800.D0*t3733*t3643*t843+3.D0/800.D0*t373
     #3*t3643*t838-3.D0/100.D0*t3794*t3796+3.D0/100.D0*t3733*t3796-3.D0/
     #400.D0*t3733*t3643*t57
      t3806 = t3733*t3767*t3609
      t3808 = t3732*t3424
      t3809 = t3808*t3748
      t3811 = t3732*t3608
      t3812 = t3606*t3795
      t3813 = t3811*t3812
      t3815 = t3613*t3795
      t3818 = t3606*t3609
      t3819 = t3808*t3818
      t3823 = t3808*t3757
      t3825 = t3619*t3795
      t3826 = t3733*t3825
      t3828 = t3794*t3825
      t3830 = t3606*t57
      t3831 = t3808*t3830
      t3834 = t3733*t3619*t838
      t3837 = t3733*t3619*t843
      t3840 = t3733*t3619*t57
      t3842 = t3428*t1387
      t3843 = t3420*s
      t3844 = t3843*Mzsq
      t3846 = t3747*t3842*t3844
      t3851 = t3634*t3639
      t3852 = t3733*t3851
      t3854 = t3634*t3737
      t3855 = t3733*t3854
      t3857 = t3794*t3851
      t3859 = t3794*t3854
      t3862 = t3747*t3606*t910
      t3864 = -13.D0/100.D0*t3806-13.D0/800.D0*t3809+13.D0/100.D0*t3813+
     #6.D0/25.D0*t3733*t3815-13.D0/100.D0*t3819-6.D0/25.D0*t3794*t3815-1
     #3.D0/800.D0*t3823-4.D0/25.D0*t3826+4.D0/25.D0*t3828+13.D0/400.D0*t
     #3831-t3834/50-t3837/50+t3840/25+t3846/3200-3.D0/50.D0*t3733*t3613*
     #t57+t3852/800+t3855/800-t3857/800-t3859/800+13.D0/800.D0*t3862
      t3867 = t3732*t3454*t3428
      t3870 = t3634*t3795
      t3871 = t3794*t3870
      t3873 = t3867*t3774
      t3875 = t3733*t3870
      t3877 = t3634*t910
      t3878 = t3733*t3877
      t3883 = t3420*Mzsq
      t3884 = t3883*t3409
      t3885 = t3751*t3884
      t3888 = t3811*t3818
      t3893 = t3428*cw2
      t3894 = t3893*t3884
      t3897 = t3780*t3884
      t3898 = t3811*t3897
      t3901 = t3733*t3619*t3609
      t3914 = t3867*t3877
      t3917 = t3747*t3842*t3763
      t3920 = t3747*t3842*t3753
      t3922 = t3747*t3830
      t3924 = -3.D0/100.D0*t3867*t3777+t3871/400+t3873/50-t3875/400+t387
     #8/1600+3.D0/100.D0*t3733*t3643*t3609-3.D0/100.D0*t3811*t3885+13.D0
     #/100.D0*t3888+6.D0/25.D0*t3733*t3613*t3609-6.D0/25.D0*t3811*t3894+
     #4.D0/25.D0*t3898-4.D0/25.D0*t3901+3.D0/100.D0*t3733*t3613*t838+3.D
     #0/100.D0*t3733*t3613*t843-3.D0/800.D0*t3867*t3771-3.D0/800.D0*t374
     #7*t3751*t3844+t3914/3200-t3917/6400-t3920/6400+13.D0/800.D0*t3922
      t3925 = t3842*t3884
      t3926 = t3811*t3925
      t3929 = t3747*t3780*t3844
      t3941 = t3733*t3634*t3609
      t3943 = t3619*t3639
      t3944 = t3794*t3943
      t3946 = t3619*t3737
      t3947 = t3794*t3946
      t3954 = t3811*t3606*t3639
      t3957 = t3811*t3606*t3737
      t3959 = t3643*t3639
      t3962 = t3643*t3737
      t3966 = t3733*t3843*t3692
      t3970 = t3409*t842
      t3972 = t3733*t3843*t3970
      t3976 = t3733*t3943
      t3978 = t3733*t3946
      t3980 = t3926/400+t3929/50-3.D0/100.D0*t3747*t3893*t3844+3.D0/200.
     #D0*t3747*t3893*t3753+3.D0/200.D0*t3747*t3893*t3763-t3941/400-2.D0/
     #25.D0*t3944-2.D0/25.D0*t3947+3.D0/25.D0*t3794*t3734+3.D0/25.D0*t37
     #94*t3738-13.D0/200.D0*t3954-13.D0/200.D0*t3957+3.D0/200.D0*t3794*t
     #3959+3.D0/200.D0*t3794*t3962+13.D0/200.D0*t3966-3.D0/200.D0*t3733*
     #t3959+13.D0/200.D0*t3972-3.D0/200.D0*t3733*t3962+2.D0/25.D0*t3976+
     #2.D0/25.D0*t3978
      wI30pc1 = t3804+t3864+t3924+t3980
      t3982 = t211*wV18
      t3983 = t3982*t3634
      t3985 = t3429*t3657*Mzsq
      t3986 = t3983*t3985
      t3988 = t3982*t3643
      t3989 = t3424*t1070
      t3990 = t3429*t3989
      t3993 = t3982*t3650
      t3994 = t57*t3409
      t3995 = t3608*t3994
      t3996 = t3993*t3995
      t3998 = t3982*t3606
      t3999 = t62*t3409
      t4001 = t3998*t3608*t3999
      t4003 = t3615*t3999
      t4006 = t3998*t3995
      t4008 = t3645*t3609
      t4009 = t3998*t4008
      t4011 = t3615*t3994
      t4014 = t3982*t3686
      t4017 = t3645*t3795
      t4018 = t3993*t4017
      t4020 = t3608*s
      t4021 = t4020*t3795
      t4022 = t3993*t4021
      t4024 = t3986/12800+3.D0/12800.D0*t3988*t3990+t3996/1600-39.D0/400
     #.D0*t4001-9.D0/400.D0*t3988*t4003+13.D0/400.D0*t4006-13.D0/400.D0*
     #t4009+3.D0/400.D0*t3988*t4011-3.D0/400.D0*t4014*t3995-t4018/800+t4
     #022/800
      t4025 = t3609*t842
      t4026 = t3645*t4025
      t4027 = t3993*t4026
      t4029 = t3424*t62
      t4030 = t4029*t3692
      t4031 = t3998*t4030
      t4033 = t4020*t4025
      t4034 = t3993*t4033
      t4036 = t3982*t3420
      t4038 = t4036*t3780*t3989
      t4040 = t3999*CsTh
      t4041 = t3608*t4040
      t4042 = t3998*t4041
      t4044 = t3615*t4040
      t4052 = t3982*t3619
      t4053 = t4052*t4044
      t4055 = t3982*t3613
      t4061 = -t4027/1600+13.D0/400.D0*t4031+t4034/1600-t4038/400-13.D0/
     #400.D0*t4042-3.D0/400.D0*t3988*t4044+3.D0/6400.D0*t4036*t3751*t398
     #9+3.D0/400.D0*t4014*t4041+t4053/25-3.D0/50.D0*t4055*t4044+3.D0/800
     #.D0*t4036*t3893*t3989
      t4063 = t3982*t3677
      t4064 = t4063*t4041
      t4066 = t3615*t1119
      t4067 = t4036*t4066
      t4069 = t3982*t3682
      t4072 = t1070*t1773
      t4073 = t3615*t4072
      t4074 = t4036*t4073
      t4078 = t4052*t4011
      t4080 = t3454*t3058
      t4081 = t3998*t4080
      t4083 = t4063*t3995
      t4085 = t3615*t1071
      t4086 = t4036*t4085
      t4088 = t3615*t1077
      t4089 = t4036*t4088
      t4091 = t3615*t1072
      t4092 = t4036*t4091
      t4094 = -t4064/25-13.D0/6400.D0*t4067+3.D0/50.D0*t4069*t4041-13.D0
     #/6400.D0*t4074-3.D0/50.D0*t4069*t3995-t4078/25-13.D0/3200.D0*t4081
     #+t4083/25-13.D0/6400.D0*t4086+13.D0/1600.D0*t4089-39.D0/1600.D0*t4
     #092
      t4095 = t3454*t3020
      t4098 = t3454*t3056
      t4102 = t3429*t3989*CsTh
      t4103 = t3983*t4102
      t4106 = t3429*t3989*t842
      t4107 = t3983*t4106
      t4110 = t3429*t4029*s
      t4115 = t3454*t1119
      t4116 = t3998*t4115
      t4118 = t4052*t4110
      t4120 = t54*t3409
      t4121 = t4120*CsTh
      t4122 = t3615*t4121
      t4125 = t3454*t1071
      t4126 = t3998*t4125
      t4132 = -3.D0/1600.D0*t4014*t4095-3.D0/3200.D0*t4014*t4098-t4103/5
     #1200-t4107/51200+9.D0/400.D0*t4055*t4110-3.D0/400.D0*t4055*t3985-1
     #3.D0/12800.D0*t4116-3.D0/200.D0*t4118-3.D0/1600.D0*t3988*t4122-13.
     #D0/12800.D0*t4126+3.D0/12800.D0*t3988*t4102+3.D0/12800.D0*t3988*t4
     #106
      t4135 = t4052*t3985
      t4137 = t4052*t4102
      t4139 = t4052*t4106
      t4141 = t3615*t4120
      t4142 = t4036*t4141
      t4144 = t3608*t4121
      t4147 = t3657*t3692
      t4148 = t3998*t4147
      t4150 = t3998*t4144
      t4152 = t3615*t3058
      t4153 = t4052*t4152
      t4158 = t3998*t3645*t1093
      t4161 = t3998*t3657*t720
      t4163 = t4135/200-t4137/800-t4139/800+13.D0/1600.D0*t4142+3.D0/160
     #0.D0*t4014*t4144+13.D0/1600.D0*t4148-13.D0/1600.D0*t4150-t4153/100
     #+3.D0/200.D0*t4055*t4152-13.D0/1600.D0*t4158+13.D0/800.D0*t4161
      t4166 = t3615*t1070
      t4167 = t4036*t4166
      t4170 = t3429*t3989*t1773
      t4171 = t4052*t4170
      t4173 = t3454*t4072
      t4174 = t3998*t4173
      t4180 = t3615*t3056
      t4183 = t3615*t3020
      t4186 = t4120*t1773
      t4187 = t3608*t4186
      t4188 = t4063*t4187
      t4190 = t3615*t4186
      t4191 = t4052*t4190
      t4195 = 3.D0/1600.D0*t3988*t4152-13.D0/6400.D0*t4167-t4171/800-13.
     #D0/12800.D0*t4174+3.D0/1600.D0*t4055*t4170+3.D0/12800.D0*t3988*t41
     #70-3.D0/200.D0*t4055*t4180-3.D0/100.D0*t4055*t4183-t4188/100+t4191
     #/100+3.D0/200.D0*t4069*t4187
      t4197 = t3606*t3406
      t4198 = t3982*t4197
      t4199 = t4198*t4017
      t4201 = t3998*t4017
      t4207 = t3998*t4026
      t4215 = t4198*t4026
      t4220 = t4036*t3842*t3989
      t4222 = 13.D0/200.D0*t4199-13.D0/200.D0*t4201+3.D0/200.D0*t4014*t4
     #017-3.D0/200.D0*t4014*t4021-13.D0/400.D0*t4207-3.D0/800.D0*t3988*t
     #4183-3.D0/200.D0*t4055*t4122+3.D0/400.D0*t4014*t4026+13.D0/400.D0*
     #t4215-3.D0/400.D0*t4014*t4033-t4220/25600
      t4226 = t3998*t3657*t1096
      t4228 = t3983*t4003
      t4234 = t3657*t3970
      t4235 = t3998*t4234
      t4239 = t3983*t4170
      t4241 = t4063*t4017
      t4243 = t4063*t4021
      t4249 = 3.D0/200.D0*t4069*t4144+13.D0/1600.D0*t4226+3.D0/1600.D0*t
     #4228-3.D0/1600.D0*t3988*t4180+3.D0/25.D0*t4069*t4017+13.D0/1600.D0
     #*t4235-3.D0/25.D0*t4069*t4021-t4239/51200-2.D0/25.D0*t4241+2.D0/25
     #.D0*t4243+3.D0/50.D0*t4069*t4026-3.D0/50.D0*t4069*t4033
      t4253 = t4120*t842
      t4254 = t3608*t4253
      t4255 = t3998*t4254
      t4257 = t3615*t4253
      t4266 = t4063*t4033
      t4268 = t3657*t3409
      t4269 = t3623*t4268
      t4274 = t3409*t1773
      t4275 = t3657*t4274
      t4276 = t3998*t4275
      t4278 = t4036*t4003
      t4280 = t3998*t4187
      t4282 = -13.D0/1600.D0*t4255-3.D0/1600.D0*t3988*t4257-3.D0/1600.D0
     #*t3988*t4141+3.D0/6400.D0*t3988*t4073+3.D0/6400.D0*t3988*t4066+t42
     #66/25+3.D0/1600.D0*t3988*t4269+3.D0/1600.D0*t4014*t4254+13.D0/1600
     #.D0*t4276+39.D0/400.D0*t4278-13.D0/1600.D0*t4280
      t4285 = t4052*t4183
      t4287 = t4036*t3990
      t4291 = t3454*t1077
      t4292 = t3998*t4291
      t4294 = t4052*t4180
      t4296 = t3993*t4144
      t4298 = t3983*t4122
      t4304 = t4052*t4091
      t4306 = -3.D0/1600.D0*t3988*t4190+t4285/50-13.D0/12800.D0*t4287+3.
     #D0/1600.D0*t4014*t4187+13.D0/3200.D0*t4292+t4294/100-t4296/6400+t4
     #298/6400+9.D0/1600.D0*t3988*t4091+9.D0/200.D0*t4055*t4091-3.D0/100
     #.D0*t4304
      t4309 = t3998*t3608*t4120
      t4313 = t3983*t4152
      t4317 = t3983*t4183
      t4319 = t3983*t4180
      t4323 = t3983*t4257
      t4325 = t3993*t4254
      t4327 = t3983*t4190
      t4329 = t3993*t4187
      t4331 = -13.D0/1600.D0*t4309-3.D0/1600.D0*t3988*t4088-t4313/6400-3
     #.D0/200.D0*t4055*t4141+t4317/3200+t4319/6400+3.D0/200.D0*t4055*t42
     #69+t4323/6400-t4325/6400+t4327/6400-t4329/6400
      t4336 = t3998*t4095
      t4338 = t3998*t4098
      t4340 = t4052*t4088
      t4344 = t4063*t4095
      t4346 = t4063*t4098
      t4356 = 3.D0/800.D0*t4055*t4085+3.D0/6400.D0*t3988*t4085+13.D0/160
     #0.D0*t4336+13.D0/3200.D0*t4338+t4340/100-3.D0/200.D0*t4055*t4088+t
     #4344/100+t4346/200-3.D0/200.D0*t4069*t4095-3.D0/400.D0*t4069*t4098
     #+3.D0/400.D0*t4069*t4080+3.D0/3200.D0*t4014*t4080
      t4359 = t4052*t4141
      t4361 = t4052*t4073
      t4363 = t4052*t4269
      t4367 = t4052*t4085
      t4369 = t4052*t4066
      t4373 = t4063*t4080
      t4379 = t4063*t4144
      t4381 = t4359/100-t4361/400-t4363/100+3.D0/800.D0*t4055*t4073-t436
     #7/400-t4369/400+3.D0/800.D0*t4055*t4066-t4373/200-3.D0/200.D0*t405
     #5*t4257+3.D0/200.D0*t4069*t4254-t4379/100
      t4382 = t4063*t4254
      t4384 = t4052*t4122
      t4386 = t4052*t4257
      t4390 = t3993*t4095
      t4392 = t3983*t4085
      t4394 = t3983*t4269
      t4396 = t3993*t4098
      t4398 = t3993*t4080
      t4400 = t3983*t4066
      t4402 = t3983*t4073
      t4404 = -t4382/100+t4384/100+t4386/100-3.D0/200.D0*t4055*t4190+t43
     #90/6400-t4392/25600-t4394/6400+t4396/12800-t4398/12800-t4400/25600
     #-t4402/25600
      t4406 = t4052*t3990
      t4410 = t4029*t3409
      t4411 = t3623*t4410
      t4412 = t3983*t4411
      t4418 = t3983*t4011
      t4422 = t3454*t1072
      t4423 = t3998*t4422
      t4427 = t4063*t4026
      t4429 = t3983*t3990
      t4431 = -t4406/800+3.D0/1600.D0*t4055*t3990-3.D0/1600.D0*t4412+3.D
     #0/1600.D0*t4055*t4102+3.D0/1600.D0*t4055*t4106-t4418/1600-3.D0/320
     #0.D0*t3988*t3985-39.D0/3200.D0*t4423+9.D0/3200.D0*t3988*t4110-t442
     #7/25-t4429/51200
      t4432 = t3983*t4088
      t4434 = t3983*t4091
      t4436 = t3983*t4141
      t4438 = t3983*t4110
      t4442 = t3983*t4044
      t4444 = t4052*t4003
      t4446 = t3993*t4041
      t4448 = t4052*t4411
      t4456 = t4432/6400-3.D0/6400.D0*t4434+t4436/6400-3.D0/12800.D0*t44
     #38+9.D0/400.D0*t3988*t4411+t4442/1600+3.D0/25.D0*t4444-t4446/1600-
     #3.D0/25.D0*t4448+9.D0/50.D0*t4055*t4411-9.D0/50.D0*t4055*t4003+3.D
     #0/50.D0*t4055*t4011
      wI31pc1 = t4024+t4061+t4094+t4132+t4163+t4195+t4222+t4249+t4282+t4
     #306+t4331+t4356+t4381+t4404+t4431+t4456
      t4460 = t1070*t211
      t4461 = t4460*wV19
      t4462 = t3606*t3454
      t4465 = t3606*t3424
      t4468 = t54*t211
      t4470 = t3424*t3409
      t4471 = t3606*t4470
      t4474 = wV19*t3420
      t4475 = t4468*t4474
      t4476 = t3623*t4470
      t4479 = t4460*t4474
      t4480 = t3424*CsTh
      t4481 = t3429*t4480
      t4484 = t3606*t4480
      t4487 = t3615*t3692
      t4490 = t4480*t3409
      t4491 = t3623*t4490
      t4494 = t3424*t842
      t4495 = t3429*t4494
      t4498 = t3606*t4494
      t4501 = t3615*t3970
      t4504 = t4494*t3409
      t4505 = t3623*t4504
      t4508 = t3424*t1773
      t4512 = t3606*t4508
      t4518 = t4508*t3409
      wI32pc1 = -t4461*t4462/80-t4461*t4465/40+t4468*wV19*t4471/10-t4475
     #*t4476/10+3.D0/80.D0*t4479*t4481+3.D0/40.D0*t4461*t4484-3.D0/10.D0
     #*t4475*t4487+3.D0/10.D0*t4475*t4491-3.D0/80.D0*t4479*t4495-3.D0/40
     #.D0*t4461*t4498+3.D0/10.D0*t4475*t4501-3.D0/10.D0*t4475*t4505+t447
     #9*t3429*t4508/80+t4461*t4512/40-t4475*t3615*t4274/10+t4475*t3623*t
     #4518/10
      t4522 = t211*wV2
      t4523 = t4522*t3424
      t4524 = t3893*t3843
      t4527 = t3420*t3409
      t4528 = t3751*t4527
      t4531 = t3780*t3843
      t4532 = t4523*t4531
      t4534 = t3606*t3409
      t4535 = t4523*t4534
      t4537 = t4522*t3454
      t4538 = t3842*t3843
      t4539 = t4537*t4538
      t4541 = t4522*t3615
      t4542 = t3634*t129
      t4543 = t4541*t4542
      t4545 = t4522*t3608
      t4546 = t3780*t4527
      t4547 = t4545*t4546
      t4549 = t3893*t4527
      t4554 = t3751*t3843
      t4557 = t3843*CsTh
      t4559 = t4537*t3842*t4557
      t4561 = t3619*t129
      t4562 = t4541*t4561
      t4564 = t3643*t129
      t4568 = t4537*t3780*t4557
      t4570 = t3767*t3409
      t4571 = t3842*t4570
      t4572 = t4545*t4571
      t4574 = t3634*t3692
      t4575 = t4541*t4574
      t4577 = t3613*t129
      t4582 = t3842*t4527
      t4583 = t4545*t4582
      t4585 = t4537*t4531
      t4587 = -3.D0/50.D0*t4523*t4524+3.D0/100.D0*t4523*t4528+t4532/25-1
     #3.D0/100.D0*t4535+t4539/3200-t4543/1600+4.D0/25.D0*t4547-6.D0/25.D
     #0*t4545*t4549-3.D0/100.D0*t4545*t4528-3.D0/800.D0*t4537*t4554-t455
     #9/3200-t4562/25+3.D0/400.D0*t4541*t4564-t4568/50-t4572/400+t4575/4
     #00+3.D0/50.D0*t4541*t4577-3.D0/100.D0*t4537*t4524+t4583/400+t4585/
     #50
      t4594 = t3606*t129
      t4595 = t4537*t4594
      t4597 = t3643*t3692
      t4600 = t3751*t4570
      t4603 = t3613*t3692
      t4606 = t3606*t3692
      t4607 = t4545*t4606
      t4609 = t3893*t4570
      t4612 = t3619*t3692
      t4613 = t4541*t4612
      t4615 = t3780*t4570
      t4616 = t4545*t4615
      t4618 = t4523*t4606
      t4621 = t4523*t3429*t3843
      t4623 = t4523*t4538
      t4625 = t4523*t4582
      t4627 = t4523*t4594
      t4629 = t4523*t4546
      t4633 = t3623*t4527
      t4634 = t4523*t4633
      t4638 = t3606*s
      t4639 = t4523*t4638
      t4641 = 3.D0/800.D0*t4537*t3751*t4557+3.D0/100.D0*t4537*t3893*t455
     #7-13.D0/800.D0*t4595-3.D0/100.D0*t4541*t4597+3.D0/100.D0*t4545*t46
     #00-6.D0/25.D0*t4541*t4603-13.D0/100.D0*t4607+6.D0/25.D0*t4545*t460
     #9+4.D0/25.D0*t4613-4.D0/25.D0*t4616+13.D0/100.D0*t4618+13.D0/800.D
     #0*t4621+t4623/1600-t4625/400-13.D0/400.D0*t4627-4.D0/25.D0*t4629+6
     #.D0/25.D0*t4523*t4549+13.D0/100.D0*t4634-3.D0/400.D0*t4523*t4554+1
     #3.D0/400.D0*t4639
      wI33pc1 = t4587+t4641
      t4642 = t4460*wV20
      t4643 = t4642*t4465
      t4645 = wV20*t3424
      t4646 = t4460*t4645
      t4648 = t3429*t3613*t842
      t4653 = t4646*t3429*t3619*t842
      t4663 = t3420*t1773
      t4665 = t4646*t3429*t4663
      t4669 = t4646*t3429*t3619*t1773
      t4672 = t4468*t4645*t3406
      t4673 = t4672*t4571
      t4675 = t4468*t4645
      t4676 = t4675*t4571
      t4678 = t3420*t842
      t4679 = t4678*t3409
      t4680 = t3842*t4679
      t4681 = t4675*t4680
      t4683 = t3842*t4678
      t4684 = t4646*t4683
      t4687 = t4646*t3842*t4663
      t4689 = t4672*t4680
      t4694 = t4663*t3409
      t4695 = t3842*t4694
      t4696 = t4675*t4695
      t4698 = t4672*t4695
      t4702 = t4646*t3429*t3634*t1773
      t4705 = t3429*t3634*t842
      t4706 = t4646*t4705
      t4708 = t3751*t3767
      t4713 = t4646*t3429*t3619*CsTh
      t4715 = 13.D0/1600.D0*t4643-9.D0/400.D0*t4646*t4648+3.D0/200.D0*t4
     #653+3.D0/3200.D0*t4646*t3429*t3643*t1773+3.D0/400.D0*t4646*t3429*t
     #3613*t1773-13.D0/3200.D0*t4665-t4669/200-3.D0/1600.D0*t4673+3.D0/1
     #600.D0*t4676-3.D0/1600.D0*t4681+3.D0/6400.D0*t4684-t4687/6400+3.D0
     #/1600.D0*t4689-3.D0/1600.D0*t4642*t3615*t3643+t4696/1600-t4698/160
     #0-t4702/12800+3.D0/12800.D0*t4706+9.D0/1600.D0*t4646*t4708-3.D0/20
     #0.D0*t4713
      t4717 = t4646*t3780*t3767
      t4719 = t3893*t3767
      t4722 = t3429*t3613
      t4726 = t4646*t3429*t3619
      t4728 = t4675*t4582
      t4730 = t3634*t3409
      t4731 = t3623*t4730
      t4732 = t4675*t4731
      t4735 = t3429*t3613*CsTh
      t4739 = t3429*t3643*CsTh
      t4742 = t3429*t3767
      t4743 = t4646*t4742
      t4747 = t4675*t4606
      t4751 = t3623*t4570
      t4752 = t4675*t4751
      t4758 = t4675*t4615
      t4760 = t4672*t4615
      t4762 = t3429*t3634
      t4763 = t4646*t4762
      t4767 = t4675*t4633
      t4769 = -3.D0/100.D0*t4717+9.D0/200.D0*t4646*t4719-3.D0/400.D0*t46
     #46*t4722+t4726/200-t4728/1600+t4732/1600+9.D0/400.D0*t4646*t4735+9
     #.D0/3200.D0*t4646*t4739-39.D0/3200.D0*t4743-9.D0/400.D0*t4675*t460
     #0+39.D0/400.D0*t4747+9.D0/400.D0*t4672*t4600-39.D0/400.D0*t4752-9.
     #D0/50.D0*t4675*t4609+9.D0/50.D0*t4672*t4609+3.D0/25.D0*t4758-3.D0/
     #25.D0*t4760+t4763/12800+3.D0/400.D0*t4675*t4528+13.D0/400.D0*t4767
      t4771 = t3643*t3409
      t4772 = t3623*t4771
      t4775 = t3893*t4679
      t4778 = t3623*t4679
      t4779 = t4675*t4778
      t4781 = t3606*t3970
      t4782 = t4675*t4781
      t4784 = t3751*t4679
      t4787 = t3893*t4678
      t4793 = t4646*t3780*t4678
      t4795 = t3751*t4678
      t4798 = t3842*t3767
      t4799 = t4646*t4798
      t4801 = t3619*t3409
      t4803 = t4675*t3623*t4801
      t4805 = t4675*t4546
      t4807 = t3613*t3409
      t4808 = t3623*t4807
      t4815 = t3780*t4679
      t4816 = t4675*t4815
      t4818 = t4672*t4815
      t4821 = t4646*t3780*t4663
      t4826 = t3751*t4694
      t4829 = -3.D0/400.D0*t4675*t4772+9.D0/50.D0*t4675*t4775+39.D0/400.
     #D0*t4779-39.D0/400.D0*t4782-9.D0/400.D0*t4672*t4784-9.D0/200.D0*t4
     #646*t4787+9.D0/400.D0*t4675*t4784+3.D0/100.D0*t4793-9.D0/1600.D0*t
     #4646*t4795-3.D0/6400.D0*t4799+t4803/25-t4805/25-3.D0/50.D0*t4675*t
     #4808+3.D0/50.D0*t4675*t4549-9.D0/50.D0*t4672*t4775-3.D0/25.D0*t481
     #6+3.D0/25.D0*t4818-t4821/100+3.D0/1600.D0*t4646*t3751*t4663-3.D0/4
     #00.D0*t4675*t4826
      t4834 = t4675*t3606*t4274
      t4839 = t4675*t3623*t4694
      t4841 = t3893*t4694
      t4846 = t3780*t4694
      t4847 = t4672*t4846
      t4849 = t4675*t4846
      t4852 = t3429*t3634*CsTh
      t4853 = t4646*t4852
      t4855 = t3429*t3643
      t4858 = t4642*t4462
      t4861 = t4642*t3615*t3634
      t4863 = t4642*t4484
      t4865 = t4642*t4512
      t4867 = t4642*t4498
      t4869 = t4468*wV20
      t4870 = t4869*t4471
      t4876 = t4642*t3615*t3619
      t4879 = t3429*t3643*t842
      t4882 = t3429*t4678
      t4883 = t4646*t4882
      t4885 = 3.D0/200.D0*t4646*t3893*t4663+13.D0/400.D0*t4834+3.D0/400.
     #D0*t4672*t4826-13.D0/400.D0*t4839-3.D0/50.D0*t4675*t4841+3.D0/50.D
     #0*t4672*t4841-t4847/25+t4849/25-3.D0/12800.D0*t4853-3.D0/3200.D0*t
     #4646*t4855+13.D0/3200.D0*t4858+t4861/6400-39.D0/1600.D0*t4863-13.D
     #0/1600.D0*t4865+39.D0/1600.D0*t4867-13.D0/400.D0*t4870-3.D0/200.D0
     #*t4642*t3615*t3613+t4876/100-9.D0/3200.D0*t4646*t4879+39.D0/3200.D
     #0*t4883
      wI34pc1 = t4715+t4769+t4829+t4885
      t4887 = t4468*wV8
      t4888 = t4887*t4462
      t4890 = t4887*t4465
      t4892 = s*t211
      t4893 = t4892*wV8
      t4894 = t4893*t4471
      t4896 = wV8*t3420
      t4897 = t4892*t4896
      t4898 = t4897*t4476
      t4900 = t4468*t4896
      t4909 = t4900*t4495
      t4911 = t4887*t4498
      t4913 = t4897*t4501
      t4915 = t4897*t4505
      wI35pc1 = t4888/40+t4890/20-t4894/5+t4898/5-t4900*t4481/20-t4887*t
     #4484/10+2.D0/5.D0*t4897*t4487-2.D0/5.D0*t4897*t4491+t4909/40+t4911
     #/20-t4913/5+t4915/5
      spf36 = log(DCMPLX(t1/t714))
      t4919 = t211*t3424
      t4920 = t4919*t3623
      t4921 = t3409*spf36
      t4922 = t3643*t4921
      t4925 = t4919*t3429
      t4926 = t129*spf36
      t4928 = t4925*t3634*t4926
      t4931 = t4925*t3619*t4926
      t4933 = s*spf36
      t4934 = t3634*t4933
      t4935 = t4925*t4934
      t4937 = t3619*t4921
      t4938 = t4920*t4937
      t4940 = t4919*t3893
      t4941 = t3767*t4933
      t4944 = t3619*t4933
      t4945 = t4925*t4944
      t4947 = t4919*t3751
      t4950 = t3767*t4921
      t4953 = t3692*spf36
      t4957 = t4920*t4950
      t4961 = t4919*t3780
      t4962 = t4961*t4941
      t4967 = t4919*t3428
      t4968 = t3767*t141
      t4969 = t4967*t4968
      t4971 = -3.D0/200.D0*t4920*t4922-t4928/6400-t4931/100+t4935/6400+2
     #.D0/25.D0*t4938+3.D0/100.D0*t4940*t4941+t4945/100+3.D0/800.D0*t494
     #7*t4941-3.D0/200.D0*t4947*t4950+3.D0/25.D0*t4920*t3613*t4953-13.D0
     #/200.D0*t4957-3.D0/25.D0*t4940*t4950-t4962/50+3.D0/200.D0*t4920*t3
     #643*t4953-13.D0/800.D0*t4969
      t4972 = t3409*t124
      t4973 = t3619*t4972
      t4974 = t4967*t4973
      t4976 = t3643*t4972
      t4979 = t4919*t3406
      t4981 = t4979*t3606*t4972
      t4983 = t3613*t4972
      t4986 = t3619*t141
      t4987 = t4967*t4986
      t4993 = t3767*t4972
      t4996 = t4961*t4968
      t4998 = t3692*t124
      t5002 = t4920*t4993
      t5009 = t4961*t4993
      t5012 = t4920*t3619*t4998
      t5014 = -2.D0/25.D0*t4974+3.D0/200.D0*t4967*t4976+13.D0/200.D0*t49
     #81+3.D0/25.D0*t4967*t4983+t4987/50+3.D0/100.D0*t4940*t4968+3.D0/80
     #0.D0*t4947*t4968-3.D0/25.D0*t4940*t4993-t4996/50+3.D0/25.D0*t4920*
     #t3613*t4998-13.D0/200.D0*t5002-3.D0/200.D0*t4947*t4993+3.D0/200.D0
     #*t4920*t3643*t4998+2.D0/25.D0*t5009-2.D0/25.D0*t5012
      t5016 = t3613*t141
      t5021 = t3634*t4972
      t5022 = t4920*t5021
      t5024 = t3643*t141
      t5027 = t4925*t4986
      t5031 = t4967*t4934
      t5034 = t4920*t3634*t4998
      t5036 = t4919*t3842
      t5037 = t5036*t4993
      t5039 = t3634*t4921
      t5040 = t4967*t5039
      t5042 = t4919*t3426
      t5044 = t5042*t3606*t4933
      t5049 = t4925*t4968
      t5055 = t4925*t3619*t142
      t5057 = -3.D0/100.D0*t4967*t5016-3.D0/200.D0*t4925*t5016+t5022/800
     #-3.D0/800.D0*t4967*t5024+t5027/100-3.D0/1600.D0*t4925*t5024+t5031/
     #3200-t5034/800+t5037/800-t5040/800+13.D0/1600.D0*t5044+3.D0/200.D0
     #*t4925*t3613*t142-13.D0/1600.D0*t5049+3.D0/1600.D0*t4925*t3643*t14
     #2-t5055/100
      t5058 = t5036*t4968
      t5060 = t3613*t4921
      t5064 = t4979*t3606*t4921
      t5068 = t4967*t4937
      t5070 = t4967*t4941
      t5072 = t4967*t4950
      t5077 = t4925*t3634*t142
      t5079 = t4967*t4944
      t5081 = t3643*t4933
      t5084 = t3613*t4933
      t5087 = t4967*t4612
      t5091 = t4979*t4606
      t5093 = -t5058/3200+3.D0/25.D0*t4967*t5060+13.D0/200.D0*t5064+3.D0
     #/200.D0*t4967*t4922-2.D0/25.D0*t5068-13.D0/800.D0*t5070+13.D0/200.
     #D0*t5072-3.D0/25.D0*t4920*t4983-t5077/6400+t5079/50-3.D0/800.D0*t4
     #967*t5081-3.D0/100.D0*t4967*t5084+4.D0/25.D0*t5087-6.D0/25.D0*t496
     #7*t4603-13.D0/100.D0*t5091
      t5096 = t5042*t4594
      t5100 = t4979*t4582
      t5102 = t5042*t4531
      t5104 = t4919*t4534
      t5108 = t4919*t4638
      t5112 = t4979*t4546
      t5118 = t5042*t4538
      t5122 = t4967*t4561
      t5126 = -13.D0/800.D0*t5096-3.D0/100.D0*t4967*t4597+t5100/400+t510
     #2/50-13.D0/100.D0*t5104-3.D0/100.D0*t5042*t4524+13.D0/400.D0*t5108
     #-3.D0/800.D0*t5042*t4554+4.D0/25.D0*t5112-3.D0/100.D0*t4979*t4528-
     #6.D0/25.D0*t4979*t4549+t5118/3200+3.D0/50.D0*t4967*t4577-t5122/25+
     #3.D0/400.D0*t4967*t4564
      t5127 = t4967*t4542
      t5129 = t4967*t4574
      t5132 = t5042*t3606*t141
      t5134 = t3634*t141
      t5135 = t4967*t5134
      t5137 = t4967*t4993
      t5139 = t4967*t5021
      t5142 = t4967*t4527*t124
      t5145 = t4967*t3843*t124
      t5149 = t4920*t4973
      t5151 = t4920*t4574
      t5153 = t4925*t4542
      t5155 = t4925*t5134
      t5157 = t4925*t4561
      t5161 = -t5127/1600+t5129/400+13.D0/1600.D0*t5132+t5135/3200+13.D0
     #/200.D0*t5137-t5139/800-13.D0/200.D0*t5142+13.D0/800.D0*t5145-3.D0
     #/200.D0*t4920*t4976+2.D0/25.D0*t5149-t5151/400-t5153/3200+t5155/64
     #00-t5157/50+3.D0/800.D0*t4925*t4564
      t5164 = t4967*t4527*spf36
      t5167 = t4967*t3843*spf36
      t5169 = t4967*t4570
      t5173 = t5042*t4638
      t5175 = t3634*s
      t5176 = t4967*t5175
      t5178 = t4967*t4730
      t5180 = t4967*t4801
      t5188 = t4979*t4534
      t5190 = t3613*s
      t5193 = t4920*t4612
      t5197 = -13.D0/200.D0*t5164+13.D0/800.D0*t5167+13.D0/100.D0*t5169+
     #3.D0/100.D0*t4925*t4577+13.D0/800.D0*t5173+t5176/1600-t5178/400-4.
     #D0/25.D0*t5180+3.D0/100.D0*t4920*t4597+6.D0/25.D0*t4920*t4603+6.D0
     #/25.D0*t4967*t4807+13.D0/100.D0*t5188-3.D0/50.D0*t4967*t5190-4.D0/
     #25.D0*t5193+3.D0/100.D0*t4967*t4771
      t5198 = t4961*t4950
      t5200 = t4920*t5039
      t5202 = t3619*s
      t5203 = t4967*t5202
      t5206 = t4920*t3619*t4953
      t5208 = t3643*s
      t5213 = t4967*t4557
      t5217 = t5036*t4950
      t5221 = t4925*t4941
      t5224 = t4920*t3634*t4953
      t5232 = t5036*t4941
      t5234 = 2.D0/25.D0*t5198+t5200/800+t5203/25-2.D0/25.D0*t5206-3.D0/
     #400.D0*t4967*t5208-3.D0/200.D0*t4925*t5084-13.D0/400.D0*t5213-3.D0
     #/1600.D0*t4925*t5081+t5217/800-3.D0/25.D0*t4920*t5060-13.D0/1600.D
     #0*t5221-t5224/800+3.D0/200.D0*t4925*t3613*t4926+3.D0/1600.D0*t4925
     #*t3643*t4926-t5232/3200
      wI36pc1 = t4971+t5014+t5057+t5093+t5126+t5161+t5197+t5234
      t5237 = t3605*t1387
      t5238 = t5237*t3668
      t5240 = t5237*t3693
      t5242 = t3605*t168
      t5243 = t5242*t3625
      t5245 = t3605*t3625
      t5247 = t5242*t3668
      t5249 = t5237*t3706
      t5251 = t3605*t3424
      t5252 = t3426*s
      t5253 = t5252*Mzsq
      t5254 = t5251*t5253
      t5256 = t5237*t3625
      t5258 = t5251*t3639
      t5260 = t5242*t3693
      t5262 = t5242*t3706
      t5264 = t3426*t54
      t5266 = t5251*t5264*CsTh
      wI37pc1 = t5238/96-t5240/48+t5243/3-2.D0/3.D0*t5245-t5247/12+t5249
     #/192+t5254/6-t5256/24-t5258/3+t5260/6-t5262/24+t5266/12
      t5268 = t211*wV15
      t5269 = t5268*t3426
      t5271 = t5269*t3425*t57
      t5273 = t5268*t3424
      t5274 = t1387*s
      t5276 = t5273*t5274*t3692
      t5278 = t168*Mzsq
      t5280 = t5273*t5278*t3409
      t5282 = t5268*t3625
      t5284 = t3424*t168
      t5286 = t5269*t5284*t57
      t5289 = t5269*t3425*t838
      t5291 = t5269*t3646
      t5293 = t1387*Mzsq
      t5295 = t5273*t5293*t3409
      t5297 = t5273*t3639
      t5299 = t168*s
      t5301 = t5273*t5299*t3692
      t5304 = t5269*t5284*t838
      t5306 = t5269*t3658
      t5308 = -t5271/192+t5276/96-t5280/6+t5282/3+t5286/24-t5289/384-t52
     #91/12+t5295/48+t5297/6-t5301/12+t5304/48-t5306/24
      t5309 = t5268*t3454
      t5311 = t5309*t5274*t720
      t5314 = t5273*t5274*t3970
      t5316 = t168*CsTh
      t5318 = t5273*t5316*t3609
      t5320 = t5273*t3795
      t5323 = t5309*t5299*t720
      t5326 = t5269*t3425*t843
      t5329 = t5269*t4480*t57
      t5331 = t1387*CsTh
      t5333 = t5273*t5331*t3609
      t5335 = t5273*t3737
      t5338 = t5273*t5299*t3970
      t5341 = t5269*t5284*t843
      t5344 = t5269*t3657*t842
      t5346 = t5311/192-t5314/96+t5318/6-t5320/3-t5323/24+t5326/384+t532
     #9/12-t5333/48-t5335/6+t5338/12-t5341/48+t5344/24
      wI38pc1 = t5308+t5346
      t5347 = t4460*wV18
      t5349 = t5347*t3454*t1387
      t5351 = t4468*wV18
      t5352 = t5284*t3409
      t5353 = t5351*t5352
      t5357 = t4468*wV18*t3424*t3409
      t5360 = t5347*t3454*t168
      t5362 = t3425*t3409
      t5363 = t5351*t5362
      t5367 = t4460*wV18*t3426*t3424
      t5370 = t5347*t3454*t5331
      t5373 = t5351*t5284*t3692
      t5375 = t5351*t4490
      t5378 = t5347*t3454*t5316
      t5381 = t5351*t3425*t3692
      t5384 = t5347*t3454*CsTh
      t5386 = t5349/768+t5353/24-t5357/12-t5360/96-t5363/192+t5367/48+t5
     #370/256+t5373/8-t5375/4-t5378/32-t5381/64+t5384/16
      t5387 = t1387*t842
      t5389 = t5347*t3454*t5387
      t5392 = t5351*t5284*t3970
      t5394 = t5351*t4504
      t5396 = t168*t842
      t5398 = t5347*t3454*t5396
      t5401 = t5351*t3425*t3970
      t5404 = t5347*t3454*t842
      t5406 = t1387*t1773
      t5408 = t5347*t3454*t5406
      t5411 = t5351*t5284*t4274
      t5413 = t5351*t4518
      t5415 = t168*t1773
      t5417 = t5347*t3454*t5415
      t5420 = t5351*t3425*t4274
      t5423 = t5347*t3454*t1773
      t5425 = t5389/256+t5392/8-t5394/4-t5398/32-t5401/64+t5404/16+t5408
     #/768+t5411/24-t5413/12-t5417/96-t5420/192+t5423/48
      wI39pc1 = t5386+t5425
      t5426 = t4522*t5362
      t5428 = t4522*t1387
      t5429 = t3454*s
      t5430 = t5428*t5429
      t5432 = t4522*t5352
      t5434 = t4522*t4470
      t5436 = t4522*t168
      t5437 = t5436*t5429
      t5439 = t4522*t5429
      t5441 = t5428*t4490
      t5443 = t3454*t129
      t5444 = t5428*t5443
      t5446 = t5436*t4490
      t5448 = t4522*t4490
      t5450 = t5436*t5443
      t5452 = t3426*CsTh
      t5453 = t5452*s
      t5454 = t4523*t5453
      wI40pc1 = t5426/48-t5430/192-t5432/6+t5434/3+t5437/24-t5439/12+t54
     #41/48-t5444/192-t5446/6+t5448/3+t5450/24-t5454/12
      t5456 = t211*wV20
      t5457 = t5456*t3424
      t5458 = t3426*t1070
      t5460 = t5457*t5458*t842
      t5462 = t5456*t168
      t5463 = t5462*t4410
      t5465 = t5456*t1387
      t5466 = t5465*t4410
      t5469 = t5457*t5458*t1773
      t5472 = t5457*t5264*Mzsq
      t5475 = t5457*t5252*t62
      t5477 = t3454*t1070
      t5478 = t5465*t5477
      t5480 = t5462*t5477
      t5483 = t5457*t5458*CsTh
      t5485 = t5462*t4268
      t5487 = t5457*t4253
      t5489 = t5457*t4186
      t5491 = t5457*t4040
      t5493 = -t5460/192-t5463/8+t5466/64+t5469/192+t5472/48-t5475/16-t5
     #478/3072+t5480/384+t5483/192-t5485/96+t5487/48-t5489/48-t5491/12
      t5494 = t5465*t4268
      t5496 = t5457*t3994
      t5498 = t5457*t4121
      t5500 = t5465*t4008
      t5502 = t5462*t4275
      t5504 = t5462*t4422
      t5506 = t5465*t4275
      t5508 = t5462*t4147
      t5510 = t5462*t4008
      t5512 = t5462*t4030
      t5514 = t5465*t4030
      t5516 = t5465*t4422
      t5518 = t5462*t4291
      t5520 = t5462*t4173
      t5522 = t5494/768-t5496/12-t5498/48-t5500/192+t5502/96+t5504/32-t5
     #506/768+t5508/96+t5510/24+t5512/24-t5514/192-t5516/256-t5518/96-t5
     #520/384
      t5524 = t5465*t4125
      t5526 = t5465*t4115
      t5528 = t5465*t4291
      t5530 = t5465*t4173
      t5532 = t5264*t1096
      t5533 = t5457*t5532
      t5535 = t5264*t720
      t5536 = t5457*t5535
      t5538 = t5462*t4115
      t5540 = t5462*t4125
      t5542 = t5465*t4234
      t5544 = t5465*t4147
      t5546 = t5252*t1093
      t5547 = t5457*t5546
      t5549 = t57*t3692
      t5550 = t5457*t5549
      t5552 = t57*t3970
      t5553 = t5457*t5552
      t5555 = t5524/3072-t5526/3072+t5528/768+t5530/3072+t5533/48-t5536/
     #24+t5538/384-t5540/384+t5542/768-t5544/768+t5547/48+t5550/6-t5553/
     #12
      t5556 = t5462*t4234
      t5558 = t5456*t5284
      t5559 = t5558*t5532
      t5561 = t5558*t5535
      t5563 = t5456*t3425
      t5564 = t5563*t5549
      t5566 = t5563*t5535
      t5568 = t5563*t5552
      t5570 = t5563*t5546
      t5572 = t5563*t5532
      t5574 = t5558*t5546
      t5576 = t5558*t5549
      t5578 = t5558*t5552
      t5580 = t5456*t4268
      t5582 = t5456*t4410
      t5584 = t5456*t5477
      t5586 = -t5556/96-t5559/96+t5561/48+t5564/96-t5566/384-t5568/192+t
     #5570/768+t5572/768-t5574/96-t5576/12+t5578/24+t5580/48+t5582/4-t55
     #84/192
      wI41pc1 = t5493+t5522+t5555+t5586
      spf35 = log(DCMPLX(-t1/t978))
      t5590 = t211*t168
      t5591 = t5590*t3424
      t5592 = t3692*spf35
      t5593 = t5591*t5592
      t5595 = t5252*spf35
      t5596 = t5591*t5595
      t5598 = t211*t1387
      t5599 = t5598*t3424
      t5600 = t5599*t5595
      t5602 = t5591*t5453
      t5604 = t5599*t5453
      t5606 = t5599*t5592
      t5608 = t5252*t124
      t5609 = t5591*t5608
      t5611 = t5599*t5608
      t5613 = t5042*t142
      t5615 = t5599*t4998
      t5617 = t5591*t4998
      t5620 = t5042*t129*spf35
      t5622 = t5452*t141
      t5623 = t5599*t5622
      t5626 = t5452*s*spf35
      t5627 = t5599*t5626
      t5629 = t5591*t5626
      t5631 = t5591*t5622
      t5633 = t4919*t3409
      t5635 = t4919*t4972
      t5637 = t5593/12-t5596/48+t5600/384-t5602/24+t5604/192-t5606/96-t5
     #609/48+t5611/384+t5613/24-t5615/96+t5617/12+t5620/24+t5623/384+t56
     #27/384-t5629/48-t5631/48-t5633/3-t5635/6
      t5639 = t4919*t3409*spf35
      t5641 = t4919*t5252
      t5643 = t5590*t4470
      t5645 = t5598*t4470
      t5647 = t4919*t3692
      t5649 = t4470*spf35
      t5650 = t5590*t5649
      t5652 = t4919*t5595
      t5654 = t4919*t5592
      t5656 = t4470*t124
      t5657 = t5598*t5656
      t5659 = t5598*t5649
      t5661 = t4919*t5453
      t5663 = t4919*t5608
      t5665 = t4919*t4998
      t5667 = t5590*t5656
      t5669 = t5598*t4490
      t5671 = t5590*t4490
      t5673 = t5598*t5429
      t5675 = t5590*t5429
      t5677 = -t5639/6+t5641/12+t5643/6-t5645/48-t5647/3+t5650/12+t5652/
     #24-t5654/6-t5657/96-t5659/96+t5661/12+t5663/24-t5665/6+t5667/12-t5
     #669/48+t5671/6+t5673/192-t5675/24
      wI42pc1 = t5637+t5677
      wI43pc1 = t5238/8-t5247/2-t5258/2+t5260/2+t5249/16+t5254/2+t5243-t
     #5245-t5256/4-t5240/8-t5262/4+t5266/4
      t5700 = -t5271/16+t5286/4+t5297/4-t5301/4-t5289/32-t5291/4-t5280/2
     #+t5282/2+t5295/8+t5276/16+t5304/8-t5306/8
      t5713 = t5311/16-t5323/4-t5335/4+t5338/4+t5326/32+t5329/4+t5318/2-
     #t5320/2-t5333/8-t5314/16-t5341/8+t5344/8
      wI44pc1 = t5700+t5713
      t5726 = -t5363/32+t5349/64+t5353/8-t5357/8-t5360/16+t5367/16-3.D0/
     #32.D0*t5381+3.D0/64.D0*t5370+3.D0/8.D0*t5373-3.D0/8.D0*t5375-3.D0/
     #16.D0*t5378+3.D0/16.D0*t5384
      t5739 = -3.D0/32.D0*t5401+3.D0/64.D0*t5389+3.D0/8.D0*t5392-3.D0/8.
     #D0*t5394-3.D0/16.D0*t5398+3.D0/16.D0*t5404-t5420/32+t5408/64+t5411
     #/8-t5413/8-t5417/16+t5423/16
      wI45pc1 = t5726+t5739
      wI46pc1 = t5426/8-t5430/16-t5432/2+t5434/2+t5437/4-t5439/4+t5441/8
     #-t5444/16-t5446/2+t5448/2+t5450/4-t5454/4
      t5765 = -t5460/64-3.D0/8.D0*t5463+3.D0/32.D0*t5466+t5469/64+t5472/
     #16-3.D0/16.D0*t5475-t5478/256+t5480/64+t5483/64-t5485/32+t5487/32-
     #t5489/32-t5491/8
      t5780 = t5494/128-t5496/8-t5498/32-t5500/32+t5502/32+3.D0/16.D0*t5
     #504-t5506/128+t5508/32+t5510/8+t5512/8-t5514/32-3.D0/64.D0*t5516-t
     #5518/16-t5520/64
      t5795 = t5524/256-t5526/256+t5528/64+t5530/256+t5533/16-t5536/8+t5
     #538/64-t5540/64+t5542/128-t5544/128+t5547/16+t5550/4-t5553/8
      t5810 = -t5556/32-t5559/16+t5561/8+t5564/16-t5566/32-t5568/32+t557
     #0/64+t5572/64-t5574/16-t5576/4+t5578/8+t5580/32+3.D0/8.D0*t5582-t5
     #584/64
      wI47pc1 = t5765+t5780+t5795+t5810
      t5830 = t5593/4-t5596/8+t5600/32-t5602/4+t5604/16-t5606/16-t5609/8
     #+t5611/32+t5613/8-t5615/16+t5617/4+t5620/8+t5623/32+t5627/32-t5629
     #/8-t5631/8-t5633/2-t5635/4
      t5849 = -t5639/4+t5641/4+t5643/2-t5645/8-t5647/2+t5650/4+t5652/8-t
     #5654/4-t5657/16-t5659/16+t5661/4+t5663/8-t5665/4+t5667/4-t5669/8+t
     #5671/2+t5673/16-t5675/4
      wI48pc1 = t5830+t5849
      t5850 = abs(DCMPLX(Mzsq))
      t5851 = t5850**2
      t5852 = t57-t5851
      t5853 = t5852**2
      t5854 = 1/t5853
      t5855 = t5854*t62
      t5857 = s*t5851
      t5858 = abs(DCMPLX(t3418))
      t5859 = t5858**2
      t5860 = abs(DCMPLX(cw2))
      t5861 = t5860**2
      t5862 = 1/t5861
      t5863 = t5859*t5862
      t5864 = t5857*t5863
      t5867 = t5854*t168
      t5868 = t3605*t5867
      t5869 = t5851*CsTh
      t5870 = t1077*t5869
      t5872 = t5854*t54
      t5874 = Mzsq*t5851
      t5875 = t5863*CsTh
      t5876 = t5874*t5875
      t5879 = t54*t1811
      t5880 = t5879*t5863
      t5883 = t5867*t54
      t5887 = t5854*t1387
      t5891 = t3605*t5854
      t5894 = t5851**2
      t5895 = t5894*CsTh
      t5896 = t5274*t5895
      t5899 = t1387*t1070
      t5900 = t5899*t1093
      t5903 = t62*t5851
      t5904 = t5299*t5903
      t5907 = t5299*t5895
      t5910 = t168*t1070
      t5911 = t5910*t1093
      t5914 = t5274*t5903
      t5916 = -t3605*t5855*t5864/2+t5868*t5870+t3605*t5872*t5876/4-t5868
     #*t5880/2-t3605*t5883*t5876/4-t3605*t5887*t5870/2-t5891*t5870/2+t58
     #91*t5896/4+t5891*t5900/4-2*t5891*t5904-t5891*t5907/2-t5891*t5911/2
     #+t5891*t5914
      t5919 = t5854*Mzsq
      t5920 = t5919*t5894
      t5923 = t5872*t1811
      t5926 = t1072*t5851
      t5928 = t5278*t5894
      t5930 = t168*t54
      t5931 = t5930*t1811
      t5933 = t5293*t5894
      t5936 = t1387*t54
      t5937 = t5936*t1811
      t5940 = s*t5894
      t5941 = t5940*CsTh
      t5944 = t1070*t62
      t5945 = t5944*CsTh
      t5948 = t5944*t5875
      t5951 = t5851*t5859
      t5952 = t5951*t5862
      t5953 = t1072*t5952
      t5956 = t5854*t1070
      t5958 = t62*t5859
      t5959 = t5862*CsTh
      t5960 = t5958*t5959
      t5963 = t5891*t5880/2-t3605*t5920/2-t3605*t5923/2+t5891*t5926+t589
     #1*t5928+t5891*t5931-t5891*t5933/2-t5891*t5937/2+t5891*t5941/4+t589
     #1*t5945/4+t5868*t5948/4+t5868*t5953/2-t3605*t5956*t5960/4
      wI49pc1 = t5916+t5963
      t5964 = t5854*t842
      t5965 = t3732*t5964
      t5969 = t5854*CsTh
      t5970 = t3732*t5969
      t5976 = t3732*t5854
      t5980 = Mzsq*t5894
      t5992 = t62*t842
      t6002 = -t5965*t1077*t5952/8-t5970*t5880/4-t3732*t5872*t5876/8+t59
     #76*t5904-t5976*t5880/4-t5976*t5316*t5980/2+t5976*t5396*t5940/4-t59
     #76*t5316*t5879/2-t5976*t5914/2-t5976*t5899*t5992/8-t5976*t5900/8-t
     #5976*t5896/8+t5976*t5911/4
      t6011 = t5894*t842
      t6018 = t129*t5903
      t6021 = t54*t5851
      t6022 = t1096*t6021
      t6027 = t3732*t5883
      t6034 = t5863*t842
      t6035 = t5874*t6034
      t6038 = t3732*t5867
      t6042 = t5976*t5907/4+t5976*t5910*t5992/4+t5976*t5331*t5980/4-t597
     #6*t5274*t6011/8+t5976*t5331*t5879/4-t5976*t6018/2+t5976*t6022/4+t5
     #976*t5870/4+t6027*t5876/8+t3732*t5923/4+t3732*t5920/4+t6027*t6035/
     #8-t6038*t5944*t6034/8
      t6046 = t842*s
      t6047 = t6046*t5894
      t6059 = t720*t5894
      t6076 = t5976*t5933/4-t5976*t6047/8-t5976*t1119*t62/8-t5976*t5941/
     #8-t5976*t5945/8-t5976*t5926/2+t5976*t6059/4+t5976*t838*t1811/4-t59
     #76*t5928/2+t5965*t5944*t5863/8-t5976*t5931/2+t5976*t5937/4+t5970*t
     #5953/4
      t6077 = t5867*CsTh
      t6084 = t3732*t5887
      t6093 = t1811*t5859
      t6094 = t6093*t5862
      t6110 = -t3732*t6077*t5953/4+t6038*t6018-t6038*t5870/2+t6084*t6022
     #/4+t6084*t5870/4-t6038*t6022/2-t6084*t6018/2+t6038*t838*t6094/4-t6
     #038*t5948/8+t6038*t5880/4+t3732*t5855*t5864/4-t6038*t5953/4+t3732*
     #t5956*t5960/8
      wI50pc1 = t6002+t6042+t6076+t6110
      t6112 = t3982*t5854
      t6113 = t54*t5894
      t6114 = t6113*CsTh
      t6117 = t3982*t5867
      t6118 = t1070*Mzsq
      t6119 = t6118*t5869
      t6122 = t57*t6011
      t6125 = t57*t5895
      t6128 = t2759*t5869
      t6131 = t2109*t5869
      t6134 = t5851*t842
      t6135 = t6118*t6134
      t6138 = t2109*t6134
      t6142 = t6118*t5851*t1773
      t6145 = t3982*t5887
      t6156 = -t6112*t6114/64-t6117*t6119/16-t6117*t6122/8-t6117*t6125/4
     #-t6117*t6128/4+t6117*t6131/2-t6117*t6135/16+t6117*t6138/4-t6117*t6
     #142/16+t6145*t6125/8+t6145*t6119/32+t6145*t6128/8+t6145*t6122/16+t
     #6145*t6142/32
      t6172 = t5863*t1773
      t6179 = t54*t1832
      t6192 = t3982*t5956
      t6193 = t5862*t842
      t6200 = t5854*t1754
      t6201 = t3982*t6200
      t6204 = t6145*t6135/32-t6145*t6131/4-t6145*t6138/8+t6117*t2501*t58
     #75/8+t6117*t2501*t6034/16-t6117*t1979*t5875/64-t6117*t1979*t6172/6
     #4-t6117*t1979*t6034/64-t6117*t6179*t5875/16+3.D0/16.D0*t6117*t2759
     #*t5952-t6117*t2109*t5952/16+t6117*t6118*t5952/64-t6192*t6093*t6193
     #/16-t6192*t6093*t5959/8+t6201*t5960/64
      t6207 = t5958*t5862*t1773
      t6210 = t5958*t6193
      t6214 = t1832*t5859
      t6218 = t5894*t1773
      t6222 = t5936*t6011
      t6225 = t5936*t5895
      t6244 = t6201*t6207/64+t6201*t6210/64+t3982*t5872*t6214*t5959/16-t
     #6112*t5936*t6218/64-t6112*t6222/64-t6112*t6225/64-t6112*t6131/4+t6
     #112*t6128/8+t6112*t6122/16+t6112*t6125/8+t6112*t6119/32-t6112*t613
     #8/8+t6112*t6142/32+t6112*t6135/32
      t6245 = t2501*t5863
      t6248 = t1979*t5863
      t6251 = t6179*t5863
      t6254 = t1811*t842
      t6258 = t1811*CsTh
      t6262 = t168*t1754
      t6266 = t1832*CsTh
      t6276 = t62*t1773
      t6283 = t168*t1811
      t6287 = t168*t62
      t6294 = t1387*t1754
      t6298 = -t6112*t6245/16+t6112*t6248/64+3.D0/16.D0*t6112*t6251-t611
     #2*t5910*t6254/8-t6112*t5910*t6258/4+t6112*t6262*t1093/32+t6112*t59
     #30*t6266/8-t6112*t5278*t5940/8-t6112*t5910*t5874/16+t6112*t6262*t6
     #276/32+t6112*t6262*t5992/32-3.D0/4.D0*t6112*t6283*t5857+t6112*t628
     #7*t5895/8+t6112*t5930*t5903/4-t6112*t6294*t1093/64
      t6325 = t1387*t62
      t6326 = t6325*t5895
      t6329 = t5936*t5903
      t6332 = t1387*t1811
      t6333 = t6332*t5857
      t6336 = t5293*t5940
      t6339 = t5899*t5874
      t6343 = t3982*t5854*s
      t6344 = t1811*t5851
      t6348 = -t6112*t5936*t6266/16+t6112*t5930*t6218/32+t6112*t5930*t60
     #11/32+t6112*t5930*t5895/32+t6112*t5899*t6254/16+t6112*t5899*t6258/
     #8-t6112*t6294*t6276/64-t6112*t6294*t5992/64-t6112*t6326/16-t6112*t
     #6329/8+3.D0/8.D0*t6112*t6333+t6112*t6336/16+t6112*t6339/32-3.D0/16
     #.D0*t6343*t6344*t5863
      t6349 = t5874*t5863
      t6356 = t3982*t5855
      t6364 = t6344*t5875
      t6368 = t3982*t5867*t1070
      t6371 = t5874*t6172
      t6376 = t6021*t6034
      t6379 = t5867*t62
      t6380 = t3982*t6379
      t6381 = t6021*t5875
      t6392 = -t6192*t6349/64-t6117*t6248/64-3.D0/16.D0*t6117*t6251+t635
     #6*t6021*t5863/16+t6117*t6245/16+t3982*t5867*s*t6364/16+t6368*t5876
     #/64-t6192*t6371/64-t6192*t6035/64+t6356*t6376/16-t6380*t6381/8-t63
     #80*t6376/16-t6343*t6364/16+t6356*t6381/8-t6192*t5876/64
      t6425 = t6113*t842
      t6437 = t6368*t6035/64+t6368*t6371/64-t3982*t5872*t5894/64-3.D0/16
     #.D0*t3982*t5872*t1832-t3982*t6200*t62/64+t3982*t5956*t1811/16-3.D0
     #/16.D0*t3982*t5855*t5894-3.D0/16.D0*t6112*t5936*t1832-t6112*t5936*
     #t5894/64-3.D0/16.D0*t6112*t6325*t5894-t6112*t6113*t1773/64-t6112*t
     #6425/64-t6112*t6294*t62/64+3.D0/8.D0*t6112*t5930*t1832+t6112*t5899
     #*t1811/16
      t6438 = t6179*CsTh
      t6453 = t2501*t842
      t6456 = t2501*CsTh
      t6459 = t1979*CsTh
      t6462 = t2109*t5851
      t6465 = t2759*t5851
      t6468 = t62*t5894
      t6469 = t6468*CsTh
      t6472 = t57*t5894
      t6475 = t6118*t5851
      t6481 = t1979*t842
      t6484 = -t6112*t6438/16+t6112*t5930*t5894/32+3.D0/8.D0*t6112*t6287
     #*t5894-t6112*t5910*t1811/8+t6112*t6262*t62/32+t6112*t6453/16+t6112
     #*t6456/8-t6112*t6459/64-t6112*t6462/8+3.D0/8.D0*t6112*t6465-t6112*
     #t6469/16+t6112*t6472/16+t6112*t6475/32-t6112*t1979*t1773/64-t6112*
     #t6481/64
      wI51pc1 = t6156+t6204+t6244+t6298+t6348+t6392+t6437+t6484
      t6487 = t5867*t5894
      t6490 = t5887*t5894
      t6496 = t5969*t5894
      t6499 = t4522*t5854
      t6509 = t720*t5857
      t6514 = t2109*t5863
      t6517 = t4522*t5867
      t6520 = t5874*t5863*s
      t6526 = -t4522*t6487/2+t4522*t6490/4+t4522*t5872*t62/4-t4522*t6496
     #/4-t6499*t5293*t5857/2+t6499*t6287*t838/2-t6499*t6325*t838/4+t6499
     #*t6509/2+t6499*t5278*t5857-t6499*t6514/4-t6517*t6509-t6517*t6520/4
     #+t4522*t5887*t6509/2
      t6530 = t5958*t5862
      t6536 = t4522*t5969
      t6547 = t6325*t54
      t6550 = t838*t62
      t6553 = t5331*t5894
      t6556 = t6287*t54
      t6559 = t5316*t5894
      t6562 = t5857*Mzsq
      t6565 = t4522*t5919*t5864/4-t6517*t838*t6530/4+t6517*t6514/4+t6536
     #*t6514/4+t4522*t5854*t5894/4+t4522*t6077*t6520/4-t6536*t6520/4+t64
     #99*t6547/4-t6499*t6550/4-t6499*t6553/4-t6499*t6556/2+t6499*t6559/2
     #-t6499*t6562/2
      wI52pc1 = t6526+t6565
      t6566 = wV20*t5854
      t6570 = t1754*t211
      t6571 = t6570*wV20
      t6608 = t4468*t6566*t5894/16-3.D0/16.D0*t6571*t5887*t1093+3.D0/8.D
     #0*t4869*t5867*t5895-t4642*t5887*t5874/8+t4642*t5867*t5874/4-3.D0/8
     #.D0*t6571*t5867*t5992+3.D0/8.D0*t6571*t5867*t1093-t4869*t5887*t621
     #8/16+3.D0/16.D0*t4869*t5887*t6011-3.D0/16.D0*t4869*t5887*t5895-t65
     #71*t5887*t6276/16+3.D0/16.D0*t6571*t5887*t5992+t4869*t5867*t6218/8
      t6624 = t5854*t1773
      t6631 = t6570*t6566
      t6637 = t4460*t6566
      t6643 = t4460*t6566*t168
      t6649 = -3.D0/8.D0*t4869*t5867*t6011-3.D0/8.D0*t4642*t5964*t5874+3
     #.D0/8.D0*t4642*t5969*t5874+t6571*t5867*t6276/8-t6571*t5855*t5863/1
     #6+t4642*t6624*t5874/8+t6570*t6566*t62/16-t6631*t5415*t6530/16+t663
     #1*t6207/16-t6637*t6371/16+3.D0/16.D0*t6637*t6035+3.D0/16.D0*t6643*
     #t5876+t6631*t6287*t5863/16
      t6685 = -3.D0/16.D0*t6643*t6035+t6643*t6371/16+t6637*t5406*t5874/8
     #-3.D0/8.D0*t6637*t5387*t5874-t6637*t5415*t5874/4+3.D0/4.D0*t6637*t
     #5396*t5874-3.D0/4.D0*t6637*t5316*t5874+3.D0/8.D0*t6637*t5331*t5874
     #-3.D0/16.D0*t6637*t5876-t6637*t5278*t5952/16+t4869*t6490/16+t6571*
     #t5887*t62/16-t6571*t6379/8
      t6720 = -t4869*t6487/8-3.D0/16.D0*t6571*t5969*t62+3.D0/16.D0*t6571
     #*t5855*t842-t4642*t5919*t5851/8-3.D0/16.D0*t4869*t6496+3.D0/16.D0*
     #t4869*t5964*t5894-t4869*t6624*t5894/16-t6571*t6624*t62/16-3.D0/16.
     #D0*t6631*t5316*t6530+t6637*t6349/16+3.D0/16.D0*t6631*t5960-3.D0/16
     #.D0*t6631*t6210+3.D0/16.D0*t6631*t5396*t6530
      wI53pc1 = t6608+t6649+t6685+t6720
      t6722 = t211*t5854
      t6727 = t168*t5894
      t6731 = t2109*spf36
      t6734 = t5895*spf36
      t6739 = t1387*t5894
      t6743 = t5895*t124
      t6746 = t2109*t124
      t6752 = t6722*t168
      t6755 = t6722*CsTh
      t6758 = t6722*Mzsq
      t6761 = t5874*t141
      t6764 = t6722*t1387
      t6768 = t6722*t54
      t6770 = t5958*t5862*t124
      t6774 = t1093*t54*t124
      t6779 = -t6722*t6556/2+t6722*t6559/2-t6722*t6727*spf36/4+t6722*t67
     #31/8-t6722*t6734/8+t6722*t6547/4+t6722*t6739*spf36/8-t6722*t6743/8
     #+t6722*t6746/8+t6722*t5894*t124/8+t6752*t6514/4+t6755*t6514/4+t675
     #8*t5864/4+t6752*t6761/2+t6764*t6509/2-t6752*t6509-t6768*t6770/8+t6
     #752*t6774/4-t6764*t6774/8
      t6783 = t1093*t54*spf36
      t6788 = t5874*t4933
      t6798 = t5958*t5862*spf36
      t6828 = t6755*t6761/4+t6752*t6783/4-t6764*t6783/8+t6752*t6788/2+t6
     #755*t6788/4-t6764*t6761/4-t6764*t6788/4-t6768*t6798/8+t6722*t5894*
     #spf36/8+t6722*t6739*t124/8-t6722*t6727*t124/4-t6722*t6562/2-t6722*
     #t6550/4-t6722*t6553/4+t6752*t6562-t6768*t6530/4+t6752*t6734/4-t675
     #5*t6731/8-t6758*t5857*spf36/4-t6752*t6731/4
      t6857 = t6722*t5316
      t6860 = t6722*t5278
      t6865 = t6722*t720
      t6870 = t5316*Mzsq
      t6871 = t6722*t6870
      t6872 = t5862*s
      t6874 = t5951*t6872*t124
      t6877 = t6764*t6731/8-t6764*t6550/4-t6755*t6746/8+t6752*t6550/2-t6
     #764*t6562/2-t6764*t6743/8+t6755*t6562/2-t6722*t5895/4+t6764*t6746/
     #8-t6752*t6746/4+t6722*t2109/4+t6722*t5894/4-t6758*t5857*t124/4+t68
     #57*t6520/4-t6860*t5864/4-t6857*t6514/4-t6865*t5864/4+t6752*t6743/4
     #+t6871*t6874/8
      t6880 = t6722*t5874
      t6890 = t6722*t5930
      t6893 = t6722*t5331
      t6896 = t6722*t838
      t6902 = t5951*t6872*spf36
      t6930 = -t6865*t6874/8+t6880*t5863*t141/8-t6857*t2109*t5863*t124/8
     #-t6860*t6874/8+t6890*t6770/8+t6893*t6761/4+t6896*t6770/8-t6857*t67
     #61/2-t6865*t6902/8+t6871*t6902/8-t6764*t6734/8-t6857*t2109*t5863*s
     #pf36/8-t6860*t6902/8+t6890*t6798/8+t6896*t6798/8-t6722*t6727/2+t68
     #80*t5863*t4933/8-t6857*t6788/2+t6893*t6788/4+t6722*t6739/4
      wI54pc1 = t6779+t6828+t6877+t6930
      t6932 = t211*wV1
      t6933 = t3420*t1812
      t6934 = t6932*t6933
      t6937 = 1/(-t1070+t4120+t3588-2*t3994-t1072+t3999)
      t6938 = t6937*t1811
      t6939 = t3429*t6938
      t6940 = t6934*t6939
      t6942 = t6932*t3613
      t6943 = t3428*t6937
      t6944 = t6943*t4040
      t6945 = t6942*t6944
      t6947 = t1070*t3420
      t6948 = t6932*t6947
      t6949 = t6937*t842
      t6950 = t3842*t6949
      t6951 = t6948*t6950
      t6953 = t6937*CsTh
      t6954 = t3429*t6953
      t6955 = t6948*t6954
      t6957 = t3842*t6953
      t6958 = t6948*t6957
      t6960 = t6932*t3752
      t6961 = t3406*t6937
      t6962 = t6961*t3409
      t6963 = t3842*t6962
      t6964 = t6960*t6963
      t6966 = t6932*t3843
      t6967 = t3426*t6937
      t6968 = t6967*t62
      t6969 = t3780*t6968
      t6970 = t6966*t6969
      t6972 = t6932*t3619
      t6973 = t6972*t6939
      t6975 = t6932*t5190
      t6976 = t6943*t3795
      t6977 = t6975*t6976
      t6979 = t6932*t4638
      t6980 = t6961*t3795
      t6981 = t6979*t6980
      t6983 = t6938*CsTh
      t6984 = t3429*t6983
      t6987 = t6937*t62
      t6988 = t6987*t842
      t6989 = t3429*t6988
      t6990 = t6966*t6989
      t6992 = t3893*t6968
      t6993 = t6966*t6992
      t6995 = t6932*t5208
      t6996 = t3623*t6937
      t6997 = t6996*t3795
      t6998 = t6995*t6997
      t7000 = t6932*t3682
      t7001 = t6961*t4040
      t7002 = t7000*t7001
      t7004 = t6966*t6976
      t7006 = t6995*t6976
      t7008 = t6932*t3606
      t7009 = t7008*t7001
      t7011 = t6940/40-113.D0/120.D0*t6945-t6951/7680-t6955/160+t6958/38
     #40-t6964/1920+49.D0/960.D0*t6970-11.D0/320.D0*t6973+17.D0/60.D0*t6
     #977+13.D0/80.D0*t6981-3.D0/160.D0*t6942*t6984+11.D0/640.D0*t6990-1
     #13.D0/1920.D0*t6993+t6998/120+113.D0/120.D0*t7002-13.D0/80.D0*t700
     #4-t7006/120-11.D0/40.D0*t7009
      t7012 = t6975*t6997
      t7014 = t6972*t6944
      t7016 = t6932*t3677
      t7017 = t7016*t7001
      t7019 = t6953*t3999
      t7020 = t7008*t7019
      t7022 = t5992*t3409
      t7023 = t6943*t7022
      t7024 = t6942*t7023
      t7026 = t6943*t4025
      t7027 = t6995*t7026
      t7029 = t6966*t7026
      t7031 = t6961*t7022
      t7032 = t7000*t7031
      t7034 = t6961*t4025
      t7035 = t6979*t7034
      t7037 = t6996*t4025
      t7038 = t6995*t7037
      t7040 = t6975*t7026
      t7042 = t7016*t7031
      t7044 = t6972*t7023
      t7046 = t7008*t7031
      t7048 = t6975*t7037
      t7050 = t6943*t3970
      t7051 = t6960*t7050
      t7053 = t3752*t168
      t7054 = t6932*t7053
      t7055 = t6949*t3409
      t7056 = t3623*t7055
      t7057 = t7054*t7056
      t7059 = t6949*t3999
      t7060 = t7008*t7059
      t7062 = -17.D0/60.D0*t7012+49.D0/60.D0*t7014-49.D0/60.D0*t7017+11.
     #D0/40.D0*t7020+113.D0/240.D0*t7024+t7027/240+13.D0/160.D0*t7029-11
     #3.D0/240.D0*t7032-13.D0/160.D0*t7035-t7038/240-17.D0/120.D0*t7040+
     #49.D0/120.D0*t7042-49.D0/120.D0*t7044+11.D0/80.D0*t7046+17.D0/120.
     #D0*t7048-t7051/40-t7057/960-11.D0/80.D0*t7060
      t7064 = t5*t3420
      t7066 = t6932*t7064*t1812
      t7067 = t6254*t3409
      t7068 = t6996*t7067
      t7069 = t7066*t7068
      t7071 = t6943*t7067
      t7072 = t7066*t7071
      t7074 = t6960*t7056
      t7076 = t6943*t6254
      t7077 = t6934*t7076
      t7079 = t6943*t5992
      t7080 = t6966*t7079
      t7082 = t6949*Mzsq
      t7084 = t6960*t3893*t7082
      t7088 = t7064*cw2
      t7089 = t6932*t7088
      t7092 = t3751*t7055
      t7093 = t6960*t7092
      t7095 = t6932*t7064
      t7096 = t6938*t3409
      t7100 = t6943*t3999
      t7101 = t6942*t7100
      t7103 = t6937*t3409
      t7104 = t7103*Mzsq
      t7105 = t3751*t7104
      t7106 = t6966*t7105
      t7108 = t3751*t7103
      t7109 = t6960*t7108
      t7111 = t6987*t3409
      t7112 = t3623*t7111
      t7113 = t6942*t7112
      t7115 = t6943*t3609
      t7116 = t6966*t7115
      t7118 = t3623*t7104
      t7119 = t6966*t7118
      t7121 = t6995*t7118
      t7123 = t6975*t7118
      t7125 = t7069/5-t7072/5+t7074/40+t7077/20+11.D0/320.D0*t7080+17.D0
     #/480.D0*t7084+3.D0/160.D0*t6942*t7076-3.D0/40.D0*t7089*t7071+t7093
     #/960-3.D0/40.D0*t7095*t3893*t7096+113.D0/240.D0*t7101+t7106/240+t7
     #109/960-113.D0/240.D0*t7113+13.D0/160.D0*t7116-13.D0/160.D0*t7119-
     #t7121/240+17.D0/120.D0*t7123
      t7127 = t6966*t3893*t7104
      t7129 = t6932*t54
      t7130 = t3606*t7103
      t7131 = t7129*t7130
      t7133 = t6932*t3420
      t7135 = t7133*t3780*t6938
      t7137 = t6937*Mzsq
      t7138 = t3606*t7137
      t7139 = t7129*t7138
      t7141 = t6932*t1070
      t7142 = t3606*t6967
      t7143 = t7141*t7142
      t7148 = t6932*s
      t7149 = t3606*t6987
      t7150 = t7148*t7149
      t7152 = t3643*t6943
      t7153 = t7141*t7152
      t7155 = t3634*t6943
      t7156 = t7141*t7155
      t7158 = t3606*t6949
      t7159 = t7141*t7158
      t7161 = t3606*t6953
      t7162 = t7141*t7161
      t7164 = t1812*t3428
      t7166 = t7133*t7164*t6938
      t7168 = t7133*t7100
      t7170 = t3623*t7096
      t7171 = t7066*t7170
      t7174 = t6960*t3893*t7137
      t7176 = t6953*Mzsq
      t7177 = t3429*t7176
      t7178 = t6960*t7177
      t7180 = t3751*t6953
      t7181 = t6948*t7180
      t7183 = t6953*t62
      t7184 = t3429*t7183
      t7185 = t6975*t7184
      t7187 = -17.D0/120.D0*t7127-t7131/40-11.D0/160.D0*t7135-13.D0/640.
     #D0*t7139+t7143/320+3.D0/160.D0*t7133*t3893*t6938+11.D0/320.D0*t715
     #0-t7153/3840-t7156/7680+t7159/160-t7162/80+t7166/20-11.D0/80.D0*t7
     #168+t7171/5+17.D0/480.D0*t7174+13.D0/640.D0*t7178+t7181/1920+113.D
     #0/960.D0*t7185
      t7191 = t6932*t3752*cw2
      t7192 = t7191*t7177
      t7194 = t6961*t3999
      t7195 = t7008*t7194
      t7197 = t6972*t7100
      t7199 = t6972*t7112
      t7201 = t3751*t7176
      t7202 = t6960*t7201
      t7204 = t3893*t7183
      t7205 = t6966*t7204
      t7207 = t6943*t720
      t7208 = t6960*t7207
      t7210 = t3780*t7183
      t7211 = t6966*t7210
      t7213 = t6943*t6258
      t7214 = t6972*t7213
      t7216 = t6258*t3409
      t7217 = t6996*t7216
      t7220 = t7064*t177
      t7221 = t6932*t7220
      t7222 = t7221*t7217
      t7224 = t6943*t7216
      t7225 = t7221*t7224
      t7227 = t6953*t3409
      t7228 = t3623*t7227
      t7229 = t7054*t7228
      t7231 = t6934*t6984
      t7233 = t6943*t3692
      t7234 = t6960*t7233
      t7236 = t7066*t7224
      t7238 = t6960*t7228
      t7240 = t7066*t7217
      t7242 = -17.D0/480.D0*t7192+11.D0/80.D0*t7195-49.D0/120.D0*t7197+4
     #9.D0/120.D0*t7199+t7202/480+113.D0/480.D0*t7205+13.D0/320.D0*t7208
     #-49.D0/240.D0*t7211+11.D0/80.D0*t7214-3.D0/20.D0*t7089*t7217+11.D0
     #/20.D0*t7222-11.D0/20.D0*t7225+t7229/480-t7231/20+t7234/20+2.D0/5.
     #D0*t7236-t7238/20-2.D0/5.D0*t7240
      t7243 = t6934*t7213
      t7246 = t6960*t3893*t7176
      t7248 = t6943*t1093
      t7249 = t6966*t7248
      t7251 = t3751*t7082
      t7252 = t6960*t7251
      t7254 = t3752*t1387
      t7255 = t6932*t7254
      t7256 = t7255*t7228
      t7258 = t3751*t7227
      t7259 = t6960*t7258
      t7265 = t6943*t1096
      t7266 = t6960*t7265
      t7268 = t3893*t6988
      t7269 = t6966*t7268
      t7271 = t6972*t6984
      t7273 = t3780*t6988
      t7274 = t6966*t7273
      t7276 = t6972*t7076
      t7278 = t3842*t7227
      t7279 = t6960*t7278
      t7283 = t7221*t7071
      t7285 = t6932*t5202
      t7286 = t7285*t7184
      t7288 = t3751*t6967
      t7289 = t6948*t7288
      t7291 = -t7243/10-17.D0/240.D0*t7246-11.D0/160.D0*t7249-t7252/960+
     #t7256/960-t7259/480+3.D0/20.D0*t7089*t7224-3.D0/80.D0*t6942*t7213-
     #13.D0/640.D0*t7266-113.D0/960.D0*t7269+11.D0/160.D0*t7271+49.D0/48
     #0.D0*t7274-11.D0/160.D0*t7276-t7279/960+3.D0/40.D0*t7089*t7068+11.
     #D0/40.D0*t7283-49.D0/480.D0*t7286-t7289/7680
      t7293 = t7221*t7068
      t7295 = t6967*t842
      t7296 = t3751*t7295
      t7297 = t6948*t7296
      t7299 = t6967*CsTh
      t7300 = t3751*t7299
      t7301 = t6948*t7300
      t7303 = t3842*t6967
      t7304 = t6948*t7303
      t7306 = t3429*t6949
      t7307 = t6948*t7306
      t7309 = t3842*t7103
      t7310 = t6960*t7309
      t7312 = t3751*t7137
      t7313 = t6960*t7312
      t7315 = t3893*t6987
      t7316 = t6966*t7315
      t7318 = t3780*t6987
      t7319 = t6966*t7318
      t7324 = t7095*t3780*t7096
      t7326 = t7221*t7170
      t7328 = t3751*t6962
      t7329 = t6960*t7328
      t7331 = t3623*t7103
      t7332 = t6960*t7331
      t7334 = t3429*t7082
      t7335 = t6960*t7334
      t7338 = t7095*t7164*t7096
      t7340 = t6949*t1811
      t7341 = t3429*t7340
      t7344 = t6975*t6989
      t7346 = -11.D0/40.D0*t7293-t7297/7680+t7301/3840-t7304/15360+t7307
     #/320+t7310/1920-t7313/960-113.D0/960.D0*t7316+49.D0/480.D0*t7319+3
     #.D0/40.D0*t7089*t7170+11.D0/40.D0*t7324-11.D0/40.D0*t7326-t7329/96
     #0+t7332/40-13.D0/1280.D0*t7335-t7338/5+3.D0/320.D0*t6942*t7341-113
     #.D0/1920.D0*t7344
      t7347 = t3751*t6949
      t7348 = t6948*t7347
      t7350 = t7191*t7334
      t7352 = t6972*t7341
      t7354 = t6934*t7341
      t7356 = t7255*t7056
      t7358 = t7285*t6989
      t7360 = t3842*t7295
      t7361 = t6948*t7360
      t7363 = t3842*t7055
      t7364 = t6960*t7363
      t7366 = t3429*t7137
      t7367 = t6960*t7366
      t7369 = t6967*Mzsq
      t7370 = t3751*t7369
      t7371 = t6960*t7370
      t7373 = t3429*t6987
      t7374 = t6966*t7373
      t7376 = t3606*t6937
      t7377 = t7141*t7376
      t7379 = t7054*t7177
      t7381 = t7054*t7334
      t7383 = t6966*t7184
      t7385 = t3842*t7299
      t7386 = t6948*t7385
      t7391 = t6960*t3893*t7369
      t7393 = -t7348/3840+17.D0/960.D0*t7350-11.D0/320.D0*t7352+t7354/40
     #-t7356/1920+49.D0/960.D0*t7358-t7361/15360+t7364/1920-13.D0/1280.D
     #0*t7367-t7371/1920+11.D0/640.D0*t7374+t7377/160+t7379/960-t7381/19
     #20-11.D0/320.D0*t7383+t7386/7680+3.D0/320.D0*t6942*t6939+17.D0/960
     #.D0*t7391
      wI55pc1 = t7011+t7062+t7125+t7187+t7242+t7291+t7346+t7393
      wI56pc1 = wV25
      t7396 = t211*wV3
      t7397 = t7396*t3634
      t7398 = t53*t54
      t7399 = t6943*t7398
      t7400 = t7397*t7399
      t7402 = t7396*t3606
      t7403 = t59*s
      t7405 = t7402*t6967*t7403
      t7407 = t7396*t3643
      t7408 = t6937*t59
      t7410 = t3429*t7408*s
      t7413 = t6943*t7403
      t7414 = t7397*t7413
      t7416 = t7103*t54
      t7417 = t3623*t7416
      t7418 = t7397*t7417
      t7421 = t7402*t6967*t7398
      t7423 = t6937*t53
      t7425 = t3429*t7423*t54
      t7428 = t6943*t1077
      t7429 = t7397*t7428
      t7431 = t6943*t1072
      t7432 = t7397*t7431
      t7434 = t6953*t1070
      t7435 = t3429*t7434
      t7440 = t7396*t3650
      t7441 = t6967*t3020
      t7442 = t7440*t7441
      t7444 = t7396*t3686
      t7447 = t7400/1000-t7405/500+3.D0/2000.D0*t7407*t7410-t7414/2000-t
     #7418/500+t7421/250-3.D0/1000.D0*t7407*t7425-t7429/200-t7432/2000-3
     #.D0/1000.D0*t7407*t7435+3.D0/100.D0*t7407*t7428+t7442/200-3.D0/100
     #.D0*t7444*t7441
      t7448 = t7402*t7019
      t7450 = t7397*t7410
      t7452 = t7397*t7425
      t7454 = t59*t3409
      t7456 = t7402*t6953*t7454
      t7462 = t6967*t3056
      t7465 = CsTh*t59
      t7466 = t7465*t3409
      t7467 = t6943*t7466
      t7470 = t6961*t7466
      t7471 = t7402*t7470
      t7473 = t60*s
      t7474 = t6943*t7473
      t7479 = t6961*CsTh
      t7480 = t53*t3409
      t7481 = t7480*s
      t7482 = t7479*t7481
      t7483 = t7440*t7482
      t7485 = t60*t3409
      t7486 = t7479*t7485
      t7487 = t7440*t7486
      t7489 = t6943*t7485
      t7490 = t7397*t7489
      t7492 = -4.D0/125.D0*t7448-t7450/4000+t7452/2000-4.D0/125.D0*t7456
     #-3.D0/125.D0*t7444*t7001+3.D0/125.D0*t7407*t6944+3.D0/200.D0*t7444
     #*t7462+3.D0/125.D0*t7407*t7467+4.D0/125.D0*t7471-3.D0/500.D0*t7407
     #*t7474-3.D0/125.D0*t7444*t7470-t7483/125-t7487/125-t7490/250
      t7494 = t6943*t7481
      t7495 = t7397*t7494
      t7497 = t7396*t3420
      t7498 = t6937*t1070
      t7500 = t7497*t3842*t7498
      t7502 = t7497*t7399
      t7504 = t7479*t3994
      t7505 = t7440*t7504
      t7507 = t7497*t7431
      t7509 = t3429*t7498
      t7510 = t7497*t7509
      t7515 = t6943*t1119
      t7516 = t7497*t7515
      t7518 = t6953*t3994
      t7519 = t7440*t7518
      t7521 = CsTh*t53
      t7522 = t7521*t54
      t7523 = t6943*t7522
      t7527 = t7402*t6953*t7398
      t7529 = t7465*s
      t7530 = t6943*t7529
      t7533 = t6943*t3058
      t7536 = -t7495/250-t7500/2000+t7502/125+t7505/25-t7507/250-t7510/5
     #00+3.D0/1000.D0*t7497*t3751*t7498-t7516/250-t7519/25+3.D0/250.D0*t
     #7407*t7523-2.D0/125.D0*t7527-3.D0/500.D0*t7407*t7530-3.D0/500.D0*t
     #7407*t7533
      t7537 = t6953*t7473
      t7538 = t7440*t7537
      t7541 = t7402*t6953*t7403
      t7543 = t3606*t3426
      t7544 = t7396*t7543
      t7545 = t7544*t7537
      t7547 = t6953*t4120
      t7548 = t7402*t7547
      t7550 = t6943*t4121
      t7553 = t7299*t7473
      t7556 = t7402*t7059
      t7558 = t6961*t4121
      t7561 = t7402*t7031
      t7565 = t7402*t7558
      t7570 = t7402*t6949*t7454
      t7572 = t6943*t1070
      t7573 = t7497*t7572
      t7575 = -t7538/500+t7541/125-t7545/125-4.D0/125.D0*t7548+3.D0/125.
     #D0*t7407*t7550+3.D0/500.D0*t7444*t7553+2.D0/125.D0*t7556-3.D0/125.
     #D0*t7444*t7558-2.D0/125.D0*t7561-3.D0/250.D0*t7407*t7023+4.D0/125.
     #D0*t7565+3.D0/250.D0*t7444*t7031+2.D0/125.D0*t7570-t7573/250
      t7578 = t842*t59
      t7579 = t7578*t3409
      t7580 = t6961*t7579
      t7581 = t7402*t7580
      t7583 = t6943*t7579
      t7588 = t6961*t842
      t7589 = t7588*t7481
      t7590 = t7440*t7589
      t7592 = t7588*t7485
      t7593 = t7440*t7592
      t7595 = t6949*t4120
      t7596 = t7402*t7595
      t7598 = t7588*t3994
      t7599 = t7440*t7598
      t7601 = t6943*t4253
      t7604 = t6961*t4253
      t7607 = t7402*t7604
      t7610 = t7402*t6967*t1077
      t7612 = t6949*t1072
      t7613 = t7402*t7612
      t7616 = t7402*t6967*t1072
      t7618 = -2.D0/125.D0*t7581-3.D0/250.D0*t7407*t7583+3.D0/250.D0*t74
     #44*t7580+t7590/250+t7593/250+2.D0/125.D0*t7596-t7599/50-3.D0/250.D
     #0*t7407*t7601+3.D0/250.D0*t7444*t7604-2.D0/125.D0*t7607-t7610/50-t
     #7613/250-t7616/500
      t7619 = t7397*t7509
      t7621 = t6987*s
      t7622 = t3429*t7621
      t7625 = t6949*t1077
      t7626 = t7402*t7625
      t7628 = t6949*t3994
      t7629 = t7440*t7628
      t7631 = t6949*t7485
      t7632 = t7440*t7631
      t7634 = t7578*s
      t7635 = t6943*t7634
      t7638 = t6943*t3056
      t7641 = t5992*s
      t7642 = t6943*t7641
      t7646 = t7402*t6949*t7398
      t7648 = t842*t53
      t7649 = t7648*t54
      t7650 = t6943*t7649
      t7654 = t7402*t6949*t7403
      t7656 = t6949*t7473
      t7657 = t7544*t7656
      t7659 = t7440*t7656
      t7661 = t7295*t7473
      t7664 = -t7619/4000+3.D0/2000.D0*t7407*t7622-t7626/25+t7629/50-t76
     #32/250+3.D0/1000.D0*t7407*t7635+3.D0/100.D0*t7407*t7638+3.D0/1000.
     #D0*t7407*t7642+t7646/125-3.D0/500.D0*t7407*t7650-t7654/250+t7657/2
     #50+t7659/1000-3.D0/1000.D0*t7444*t7661
      t7666 = t6967*t3058
      t7667 = t7402*t7666
      t7669 = t6949*t7481
      t7670 = t7440*t7669
      t7672 = t7397*t7435
      t7674 = t7498*t842
      t7675 = t3429*t7674
      t7676 = t7397*t7675
      t7678 = t7402*t7537
      t7680 = t6967*t7473
      t7681 = t7440*t7680
      t7683 = t6953*t7481
      t7688 = t6953*t7485
      t7695 = t7402*t7683
      t7699 = t7667/250-t7670/250+t7672/2000-t7676/4000-2.D0/125.D0*t767
     #8+t7681/2000-6.D0/125.D0*t7444*t7683+3.D0/250.D0*t7444*t7537-6.D0/
     #125.D0*t7444*t7688+6.D0/125.D0*t7444*t7482+6.D0/125.D0*t7444*t7486
     #+8.D0/125.D0*t7695-6.D0/25.D0*t7444*t7504
      t7708 = t7397*t7112
      t7710 = t7402*t7656
      t7712 = t7396*t4197
      t7713 = t7712*t7688
      t7715 = t7402*t7688
      t7717 = t7712*t7518
      t7719 = t7402*t7518
      t7721 = t7712*t7683
      t7723 = t6943*t4120
      t7724 = t7397*t7723
      t7730 = 6.D0/25.D0*t7444*t7518-3.D0/125.D0*t7444*t7589+3.D0/125.D0
     #*t7444*t7669-3.D0/500.D0*t7444*t7656-t7708/500+t7710/125-8.D0/125.
     #D0*t7713+8.D0/125.D0*t7715+8.D0/25.D0*t7717-8.D0/25.D0*t7719-8.D0/
     #125.D0*t7721+t7724/500+3.D0/125.D0*t7444*t7631-3.D0/125.D0*t7444*t
     #7592
      t7739 = t3623*t7408*t3409
      t7740 = t7397*t7739
      t7742 = t7402*t7669
      t7744 = t7402*t7441
      t7746 = t7402*t7628
      t7748 = t7712*t7669
      t7752 = t7712*t7628
      t7754 = t7397*t7523
      t7756 = t7402*t7001
      t7758 = t7712*t7631
      t7760 = t7402*t7631
      t7762 = -3.D0/25.D0*t7444*t7628+3.D0/25.D0*t7444*t7598-t7740/500-4
     #.D0/125.D0*t7742+t7744/25+4.D0/25.D0*t7746+4.D0/125.D0*t7748-3.D0/
     #1000.D0*t7444*t7666-4.D0/25.D0*t7752-t7754/500+4.D0/125.D0*t7756+4
     #.D0/125.D0*t7758-4.D0/125.D0*t7760
      t7763 = t7440*t7666
      t7765 = t7397*t7530
      t7767 = t6967*t7529
      t7770 = t6943*t3994
      t7771 = t7397*t7770
      t7773 = t6961*t7481
      t7778 = t7402*t7767
      t7782 = t6961*t7485
      t7785 = t6961*t3994
      t7790 = t7402*t7773
      t7793 = t7402*t7423*t3638
      t7795 = t7402*t7785
      t7797 = t7763/2000+t7765/1000-3.D0/1000.D0*t7444*t7767+t7771/50-3.
     #D0/125.D0*t7444*t7773+3.D0/125.D0*t7407*t7494+t7778/250+3.D0/125.D
     #0*t7407*t7489-3.D0/125.D0*t7444*t7782+3.D0/25.D0*t7444*t7785-3.D0/
     #25.D0*t7407*t7770+4.D0/125.D0*t7790-4.D0/125.D0*t7793-4.D0/25.D0*t
     #7795
      t7799 = t7137*t3638
      t7800 = t7402*t7799
      t7803 = t7402*t7423*t3609
      t7805 = t7402*t7782
      t7807 = t7440*t7558
      t7815 = t7397*t7474
      t7817 = t7497*t7100
      t7819 = t7497*t7428
      t7821 = t7497*t7723
      t7823 = t6943*t7454
      t7824 = t7497*t7823
      t7826 = t7497*t7413
      t7828 = 4.D0/25.D0*t7800-4.D0/125.D0*t7803+4.D0/125.D0*t7805+t7807
     #/250-3.D0/500.D0*t7407*t7399+3.D0/1000.D0*t7407*t7413+3.D0/1000.D0
     #*t7407*t7431+t7815/1000+2.D0/125.D0*t7817-t7819/25+2.D0/125.D0*t78
     #21+2.D0/125.D0*t7824-t7826/250
      t7829 = t6943*t1071
      t7830 = t7497*t7829
      t7832 = t7397*t7100
      t7834 = t7402*t7680
      t7836 = t7397*t7823
      t7838 = t7402*t7194
      t7845 = t7402*t6961*t7454
      t7853 = t7440*t7773
      t7855 = t7397*t7622
      t7857 = t7440*t7785
      t7859 = t7830/125+t7832/500+t7834/250+t7836/500-2.D0/125.D0*t7838-
     #3.D0/1000.D0*t7444*t7680-3.D0/250.D0*t7407*t7100-2.D0/125.D0*t7845
     #+3.D0/250.D0*t7407*t7112+3.D0/250.D0*t7407*t7739-3.D0/250.D0*t7407
     #*t7823+t7853/250-t7855/4000-t7857/50
      t7862 = t7440*t7782
      t7864 = t7402*t7462
      t7866 = t6967*t7641
      t7867 = t7402*t7866
      t7871 = t6967*t7634
      t7872 = t7402*t7871
      t7876 = t7397*t7635
      t7878 = t7397*t7650
      t7880 = t6967*t7649
      t7881 = t7402*t7880
      t7885 = t7440*t7604
      t7887 = t6967*t7522
      t7890 = t7402*t7887
      t7892 = t7862/250-t7864/50-t7867/500+3.D0/2000.D0*t7444*t7866-t787
     #2/500+3.D0/2000.D0*t7444*t7871-t7876/2000+t7878/1000+t7881/250-3.D
     #0/1000.D0*t7444*t7880-t7885/500+3.D0/500.D0*t7444*t7887-t7890/125
      t7893 = t6943*t3020
      t7894 = t7397*t7893
      t7896 = t7397*t7533
      t7898 = t7440*t7553
      t7900 = t6953*t1072
      t7901 = t7402*t7900
      t7903 = t6953*t1077
      t7904 = t7402*t7903
      t7908 = t7397*t7550
      t7910 = t7397*t7638
      t7912 = t7397*t7601
      t7914 = t7440*t7661
      t7916 = t7397*t7642
      t7918 = t7440*t7031
      t7920 = t7440*t7683
      t7924 = t7894/100+t7896/1000-t7898/1000+t7901/125+2.D0/25.D0*t7904
     #-3.D0/50.D0*t7407*t7893-t7908/250-t7910/200+t7912/500+t7914/2000-t
     #7916/2000-t7918/500+t7920/125+3.D0/1000.D0*t7407*t7515
      t7926 = t7440*t7580
      t7928 = t7397*t7023
      t7930 = t7397*t7583
      t7933 = t7402*t7423*t57
      t7935 = t7440*t7001
      t7937 = t7440*t7470
      t7939 = t7440*t7688
      t7941 = t6937*t54
      t7942 = t7941*Mzsq
      t7943 = t3429*t7942
      t7946 = t7397*t6944
      t7950 = t7397*t7467
      t7952 = t7397*t7943
      t7956 = -t7926/500+t7928/500+t7930/500+t7933/125+t7935/250+t7937/2
     #50+t7939/125+3.D0/200.D0*t7407*t7943-t7946/250-3.D0/500.D0*t7407*t
     #7829-t7950/250-t7952/400+3.D0/2000.D0*t7407*t7509
      t7958 = t7402*t6967*t1071
      t7960 = t7440*t7767
      t7962 = t7440*t7887
      t7964 = t7397*t7829
      t7966 = t7440*t7866
      t7968 = t7440*t7462
      t7971 = t7402*t6967*t1119
      t7975 = t7440*t7880
      t7977 = t7440*t7871
      t7979 = t7397*t7515
      t7986 = t7402*t6961*t4120
      t7988 = t7958/250+t7960/2000-t7962/1000+t7964/1000-t7966/4000-t796
     #8/400-t7971/500+3.D0/2000.D0*t7407*t7675+t7975/2000-t7977/4000-t79
     #79/2000-3.D0/250.D0*t7407*t7723+3.D0/250.D0*t7407*t7417-2.D0/125.D
     #0*t7986
      wI57pc1 = t7447+t7492+t7536+t7575+t7618+t7664+t7699+t7730+t7762+t7
     #797+t7828+t7859+t7892+t7924+t7956+t7988
      t7992 = t211*wV4
      t7993 = t7992*t3606
      t7994 = t6937*t103
      t7995 = t7994*t3609
      t7996 = t7993*t7995
      t7998 = t7992*t3752
      t7999 = t7998*t7363
      t8001 = t7992*t5202
      t8002 = t8001*t6989
      t8004 = t7992*t6947
      t8005 = t8004*t7347
      t8007 = t7992*t7254
      t8008 = t8007*t7056
      t8010 = t7992*t5190
      t8011 = t8010*t6989
      t8013 = t7998*t7334
      t8015 = t7998*t7278
      t8017 = t8001*t7184
      t8019 = t8004*t6957
      t8021 = t8004*t7300
      t8023 = t8004*t7296
      t8025 = 3.D0/64.D0*t7996+17.D0/5120.D0*t7999-t8002/160+11.D0/10240
     #.D0*t8005-17.D0/5120.D0*t8008+t8011/80-3.D0/1024.D0*t8013-17.D0/25
     #60.D0*t8015+t8017/80+17.D0/10240.D0*t8019-11.D0/10240.D0*t8021+11.
     #D0/20480.D0*t8023
      t8026 = t8004*t6950
      t8028 = t8004*t6954
      t8030 = t7992*t3686
      t8031 = t445*t3409
      t8032 = t6961*t8031
      t8033 = t8030*t8032
      t8035 = t7993*t8032
      t8038 = t7998*t3751*t7994
      t8040 = t7992*t7053
      t8041 = t8040*t7056
      t8043 = t8010*t7184
      t8045 = t7998*t7370
      t8047 = t7992*t3843
      t8048 = t8047*t7115
      t8050 = t8004*t7360
      t8053 = t7993*t6967*t1852
      t8055 = t7994*Mzsq
      t8056 = t3429*t8055
      t8057 = t8047*t8056
      t8059 = -17.D0/20480.D0*t8026+t8028/640+3.D0/64.D0*t8033-3.D0/64.D
     #0*t8035-59.D0/10240.D0*t8038+11.D0/2560.D0*t8041-t8043/40+3.D0/102
     #4.D0*t8045+3.D0/128.D0*t8048-17.D0/40960.D0*t8050-t8053/80-3.D0/51
     #2.D0*t8057
      t8061 = t3429*t7994
      t8062 = t7998*t8061
      t8064 = t8004*t7306
      t8066 = t8047*t7210
      t8068 = t7998*t7108
      t8070 = t7994*t3409
      t8072 = t8047*t3842*t8070
      t8074 = t8047*t7204
      t8076 = t7998*t7207
      t8078 = t7998*t7092
      t8080 = t7998*t7228
      t8082 = t7992*t3643
      t8083 = t6943*t8031
      t8084 = t8082*t8083
      t8086 = t7992*t5208
      t8087 = t3623*t8070
      t8088 = t8086*t8087
      t8090 = t8047*t8087
      t8092 = 19.D0/5120.D0*t8062-t8064/1280+t8066/40-11.D0/2560.D0*t806
     #8+7.D0/5120.D0*t8072-t8074/20+3.D0/256.D0*t8076-11.D0/2560.D0*t807
     #8+t8080/80-3.D0/64.D0*t8084-59.D0/2560.D0*t8088+19.D0/640.D0*t8090
      t8093 = t103*t3409
      t8094 = t6943*t8093
      t8095 = t8047*t8094
      t8097 = t8047*t6992
      t8099 = t6967*t103
      t8101 = t7998*t3751*t8099
      t8103 = t7998*t7258
      t8105 = t7998*t7177
      t8107 = t7998*t7201
      t8109 = t8047*t7248
      t8111 = t8040*t7228
      t8113 = t7998*t7233
      t8115 = t7992*t3619
      t8116 = t6943*t1852
      t8117 = t8115*t8116
      t8119 = t6943*t445
      t8120 = t8047*t8119
      t8122 = t7992*t3613
      t8123 = t8122*t8116
      t8125 = -19.D0/640.D0*t8095+t8097/80-59.D0/20480.D0*t8101+11.D0/12
     #80.D0*t8103+3.D0/512.D0*t8105-3.D0/256.D0*t8107+t8109/40-11.D0/128
     #0.D0*t8111-t8113/80-t8117/40-3.D0/256.D0*t8120+t8123/20
      t8128 = t7992*t7220
      t8129 = t1852*t3409
      t8130 = t6943*t8129
      t8131 = t8128*t8130
      t8133 = t7992*t7088
      t8134 = t8133*t8130
      t8136 = t6996*t8129
      t8137 = t8128*t8136
      t8139 = t8047*t7105
      t8141 = t8122*t7100
      t8143 = t7992*t4638
      t8144 = t8143*t7034
      t8146 = t8086*t7037
      t8148 = t7992*t3682
      t8149 = t8148*t7031
      t8151 = t7993*t7031
      t8153 = t7992*t3677
      t8154 = t8153*t7031
      t8156 = t7993*t7059
      t8158 = t8115*t7023
      t8160 = t8131/10-t8134/5-t8137/10-3.D0/128.D0*t8139-t8141/10-3.D0/
     #128.D0*t8144+3.D0/128.D0*t8146+t8149/10-t8151/20-t8154/20+t8156/20
     #+t8158/20
      t8161 = t8122*t7023
      t8163 = t8086*t7026
      t8165 = t8133*t8136
      t8167 = t7998*t7312
      t8169 = t6949*t103
      t8170 = t3429*t8169
      t8171 = t8007*t8170
      t8173 = t7993*t7194
      t8175 = t8115*t7112
      t8177 = t8115*t7100
      t8179 = t7064*t3428
      t8180 = t7992*t8179
      t8182 = t8180*t6961*t8129
      t8184 = t7998*t6963
      t8186 = t1224*t3409
      t8187 = t6943*t8186
      t8188 = t8086*t8187
      t8190 = t6953*t8129
      t8191 = t8180*t8190
      t8193 = -t8161/10-3.D0/128.D0*t8163+t8165/5+3.D0/512.D0*t8167-7.D0
     #/40960.D0*t8171-t8173/20-t8175/20+t8177/20-t8182/10-17.D0/5120.D0*
     #t8184-59.D0/1280.D0*t8188-t8191/5
      t8195 = t6953*t1852
      t8196 = t7993*t8195
      t8198 = t1224*Mzsq
      t8199 = t6943*t8198
      t8200 = t8086*t8199
      t8202 = t8047*t8199
      t8204 = t6996*t8186
      t8205 = t8086*t8204
      t8208 = t8143*t6961*t8186
      t8210 = t8047*t8187
      t8212 = t1224*t62
      t8213 = t6943*t8212
      t8214 = t8115*t8213
      t8216 = t6943*CsTh
      t8217 = t8216*t8129
      t8218 = t8128*t8217
      t8220 = t8122*t8213
      t8222 = t7992*t5
      t8223 = t8222*t3677
      t8224 = t7479*t8129
      t8225 = t8223*t8224
      t8227 = t8133*t8217
      t8229 = t7992*t1070
      t8230 = t8229*t7158
      t8232 = t8196/20-3.D0/128.D0*t8200+3.D0/128.D0*t8202+59.D0/1280.D0
     #*t8205-19.D0/320.D0*t8208+19.D0/320.D0*t8210+t8214/20-t8218/5-t822
     #0/10+t8225/5+2.D0/5.D0*t8227-t8230/640
      t8233 = t8229*t7142
      t8235 = t8229*t7161
      t8237 = t7992*t3420
      t8238 = t8237*t8116
      t8240 = t8222*t3682
      t8241 = t8240*t8224
      t8243 = t8180*t8224
      t8245 = t6949*t8129
      t8246 = t8180*t8245
      t8248 = t6949*t1852
      t8249 = t7993*t8248
      t8251 = t1865*Mzsq
      t8252 = t6943*t8251
      t8253 = t8086*t8252
      t8255 = t8047*t6969
      t8257 = t1865*t3409
      t8258 = t6943*t8257
      t8259 = t8086*t8258
      t8261 = t6996*t8257
      t8262 = t8086*t8261
      t8265 = t8143*t6961*t8257
      t8267 = -t8233/1280+t8235/320-t8238/40-2.D0/5.D0*t8241+t8243/5+t82
     #46/10-t8249/40+3.D0/256.D0*t8253-t8255/160+59.D0/2560.D0*t8259-59.
     #D0/2560.D0*t8262+19.D0/640.D0*t8265
      t8271 = t8047*t8258
      t8273 = t8047*t8252
      t8275 = t1865*t62
      t8276 = t6943*t8275
      t8277 = t8115*t8276
      t8279 = t8122*t8276
      t8281 = t6943*t842
      t8282 = t8281*t8129
      t8283 = t8128*t8282
      t8285 = t7588*t8129
      t8286 = t8223*t8285
      t8288 = t7998*t7331
      t8290 = t8133*t8282
      t8292 = t8180*t8285
      t8294 = t8240*t8285
      t8296 = t7998*t7328
      t8298 = t8047*t7315
      t8300 = -19.D0/640.D0*t8271-3.D0/256.D0*t8273-t8277/40+t8279/20+t8
     #283/10-t8286/10-t8288/160-t8290/5-t8292/10+t8294/5+11.D0/2560.D0*t
     #8296+t8298/40
      t8301 = t8229*t7155
      t8303 = t8047*t7318
      t8305 = t7992*s
      t8306 = t8305*t7149
      t8308 = t6953*t103
      t8309 = t3429*t8308
      t8310 = t8040*t8309
      t8313 = t7998*t3842*t8308
      t8315 = t7998*t8309
      t8317 = t7998*t7309
      t8319 = t8047*t7026
      t8321 = t8143*t6980
      t8323 = t8086*t6997
      t8325 = t7998*t7050
      t8328 = t7998*t3842*t8169
      t8330 = -17.D0/20480.D0*t8301-t8303/80-t8306/80+59.D0/10240.D0*t83
     #10+7.D0/10240.D0*t8313-19.D0/2560.D0*t8315+17.D0/5120.D0*t8317+3.D
     #0/128.D0*t8319+3.D0/64.D0*t8321-3.D0/64.D0*t8323+t8325/160-7.D0/20
     #480.D0*t8328
      t8332 = t7998*t8170
      t8334 = t8040*t8170
      t8336 = t7992*t5175
      t8337 = t8336*t8087
      t8339 = t6943*t1224
      t8340 = t7998*t8339
      t8342 = t7992*t54
      t8343 = t3606*t7994
      t8344 = t8342*t8343
      t8346 = t8336*t8204
      t8348 = t6967*t8212
      t8349 = t8148*t8348
      t8351 = t8047*t7184
      t8353 = t8153*t8348
      t8355 = t8148*t7001
      t8357 = t6967*t8251
      t8358 = t8143*t8357
      t8360 = t8004*t7385
      t8362 = 19.D0/5120.D0*t8332-59.D0/20480.D0*t8334-7.D0/5120.D0*t833
     #7-19.D0/1280.D0*t8340+19.D0/2560.D0*t8344+7.D0/2560.D0*t8346-t8349
     #/20+t8351/80+t8353/40-t8355/5-3.D0/512.D0*t8358+17.D0/20480.D0*t83
     #60
      t8363 = t3429*t6937
      t8364 = t8363*t8251
      t8365 = t8086*t8364
      t8367 = t6967*t8275
      t8368 = t7993*t8367
      t8371 = t7998*t3751*t8169
      t8373 = t8336*t8258
      t8375 = t6943*t1865
      t8376 = t7998*t8375
      t8378 = t8342*t7138
      t8380 = t8336*t8261
      t8382 = t7993*t7001
      t8384 = t8153*t7001
      t8386 = t7998*t7056
      t8388 = t7993*t7019
      t8390 = t8115*t6944
      t8392 = 3.D0/512.D0*t8365-t8368/80-59.D0/10240.D0*t8371+7.D0/5120.
     #D0*t8373+19.D0/2560.D0*t8376-3.D0/512.D0*t8378-7.D0/5120.D0*t8380+
     #t8382/10+t8384/10-t8386/160-t8388/10-t8390/10
      t8395 = t8086*t6976
      t8397 = t8122*t6944
      t8399 = t8047*t6976
      t8402 = t8047*t3751*t8070
      t8405 = t8047*t3751*t8055
      t8407 = t8047*t7273
      t8409 = t7992*t7064
      t8410 = t8409*t8130
      t8412 = t6953*t8031
      t8413 = t7993*t8412
      t8415 = t8047*t7268
      t8417 = t8342*t7130
      t8419 = t7992*t4197
      t8420 = t8419*t8412
      t8422 = t8004*t7288
      t8424 = 3.D0/64.D0*t8395+t8397/5-3.D0/64.D0*t8399+59.D0/2560.D0*t8
     #402+3.D0/256.D0*t8405-t8407/80+t8410/10-3.D0/32.D0*t8413+t8415/40+
     #t8417/160+3.D0/32.D0*t8420+11.D0/20480.D0*t8422
      t8426 = t8030*t7479*t8031
      t8428 = t7998*t7265
      t8430 = t8363*t8198
      t8431 = t8086*t8430
      t8433 = t7993*t8348
      t8435 = t8007*t8309
      t8437 = t8030*t8412
      t8439 = t8047*t6989
      t8441 = t6949*t8031
      t8442 = t7993*t8441
      t8444 = t8086*t8056
      t8447 = t7998*t3842*t8099
      t8449 = t8229*t7152
      t8452 = t8030*t7588*t8031
      t8454 = -3.D0/32.D0*t8426-3.D0/512.D0*t8428-3.D0/256.D0*t8431+t843
     #3/40+7.D0/20480.D0*t8435+3.D0/32.D0*t8437-t8439/160+3.D0/64.D0*t84
     #42+3.D0/512.D0*t8444-7.D0/40960.D0*t8447+11.D0/10240.D0*t8449+3.D0
     #/64.D0*t8452
      t8456 = t8047*t7079
      t8458 = t8237*t7100
      t8460 = t8419*t8441
      t8462 = t8030*t8441
      t8464 = t7998*t7366
      t8467 = t7998*t3842*t7994
      t8469 = t7998*t7251
      t8471 = t8047*t7373
      t8473 = t8040*t7177
      t8475 = t8007*t7228
      t8477 = t8040*t7334
      t8479 = t8004*t7180
      t8481 = -t8456/80+t8458/20-3.D0/64.D0*t8460-3.D0/64.D0*t8462-3.D0/
     #1024.D0*t8464-7.D0/20480.D0*t8467+3.D0/512.D0*t8469-t8471/160-3.D0
     #/512.D0*t8473+17.D0/2560.D0*t8475+3.D0/1024.D0*t8477-11.D0/5120.D0
     #*t8479
      t8482 = t8229*t7376
      t8484 = t7994*t62
      t8485 = t3429*t8484
      t8486 = t8115*t8485
      t8488 = t8148*t8367
      t8490 = t8153*t8367
      t8492 = t8122*t8485
      t8494 = t8122*t7112
      t8496 = t8086*t7118
      t8498 = t8047*t7118
      t8500 = t6967*t8198
      t8501 = t8143*t8500
      t8503 = t8004*t7303
      t8506 = t7998*t3751*t8308
      t8508 = t8336*t8187
      t8510 = -t8482/640-t8486/80+t8488/40-t8490/80+t8492/40+t8494/10+3.
     #D0/128.D0*t8496-3.D0/128.D0*t8498+3.D0/256.D0*t8501-17.D0/40960.D0
     #*t8503+59.D0/5120.D0*t8506-7.D0/2560.D0*t8508
      wI58pc1 = t8025+t8059+t8092+t8125+t8160+t8193+t8232+t8267+t8300+t8
     #330+t8362+t8392+t8424+t8454+t8481+t8510
      t8514 = t211*wV6
      t8515 = t8514*t177
      t8517 = 1/(-t57+t3609+t54-t3638)
      t8518 = t8517*t3428
      t8519 = t3420*t1811
      t8520 = t8518*t8519
      t8521 = t8515*t8520
      t8523 = t8514*t5910
      t8524 = t8518*t3767
      t8525 = t8523*t8524
      t8527 = t8517*t3426
      t8528 = t8527*t3606
      t8529 = t8523*t8528
      t8531 = t8518*t4678
      t8532 = t8523*t8531
      t8534 = t1070*t8517
      t8535 = t8514*t8534
      t8536 = t8535*t4742
      t8538 = t8535*t4882
      t8541 = t8527*t3606*t842
      t8542 = t8523*t8541
      t8544 = t8514*t175
      t8545 = t8518*t3884
      t8546 = t8544*t8545
      t8548 = t177*t5
      t8549 = t8514*t8548
      t8550 = t8519*t3409
      t8551 = t8518*t8550
      t8552 = t8549*t8551
      t8554 = s*t8517
      t8555 = t8514*t8554
      t8556 = t8555*t3818
      t8558 = t3420*t62
      t8559 = t8518*t8558
      t8560 = t8544*t8559
      t8562 = t177*t8517
      t8563 = t8514*t8562
      t8564 = t3606*t3999
      t8565 = t8563*t8564
      t8567 = t8521/60-t8525/60+t8529/240+t8532/120+t8536/30-t8538/60+t8
     #542/240+t8546/3-t8552/15-2.D0/15.D0*t8556+t8560/24+4.D0/15.D0*t856
     #5
      t8568 = t8554*t3406
      t8569 = t8514*t8568
      t8570 = t3606*t4025
      t8571 = t8569*t8570
      t8573 = t175*t8517
      t8574 = t8514*t8573
      t8575 = t4197*t4025
      t8576 = t8574*t8575
      t8578 = cw2*t8517
      t8579 = t8514*t8578
      t8580 = t3606*t7022
      t8581 = t8579*t8580
      t8584 = t8514*t8578*t3406
      t8585 = t8584*t8580
      t8588 = t8514*t8562*t3406
      t8589 = t8588*t8580
      t8591 = t8563*t8580
      t8593 = t8555*t8570
      t8595 = t8574*t8570
      t8597 = t4197*t3795
      t8598 = t8574*t8597
      t8600 = t8569*t3812
      t8602 = t3606*t4040
      t8603 = t8588*t8602
      t8605 = t8584*t8602
      t8607 = 2.D0/15.D0*t8571-t8576/3-t8581/6+t8585/6-4.D0/15.D0*t8589+
     #4.D0/15.D0*t8591-2.D0/15.D0*t8593+t8595/3+2.D0/3.D0*t8598-4.D0/15.
     #D0*t8600+8.D0/15.D0*t8603-t8605/3
      t8609 = t8514*t1070
      t8610 = t8609*t8528
      t8612 = t8609*t8531
      t8614 = t8609*t8524
      t8616 = t8579*t8602
      t8619 = t8514*t1812*t8520
      t8621 = t8514*t54
      t8622 = t8518*t3883
      t8623 = t8621*t8622
      t8625 = t8518*t4527
      t8626 = t8621*t8625
      t8628 = t8563*t8602
      t8630 = t8555*t3812
      t8632 = t8514*t168
      t8634 = t8632*t8534*t3606
      t8636 = t8574*t3812
      t8638 = t8514*t1067
      t8639 = t3606*Mzsq
      t8640 = t8527*t8639
      t8641 = t8638*t8640
      t8643 = -t8610/60-t8612/30+t8614/15+t8616/3-t8619/60+t8623/30+2.D0
     #/15.D0*t8626-8.D0/15.D0*t8628+4.D0/15.D0*t8630+t8634/120-2.D0/3.D0
     #*t8636-t8641/24
      t8644 = t54*t8517
      t8645 = t8514*t8644
      t8646 = t8645*t4778
      t8648 = t8579*t8564
      t8650 = t3429*t3883
      t8651 = t8645*t8650
      t8653 = t3606*t62
      t8654 = t8527*t8653
      t8655 = t8544*t8654
      t8657 = t3883*CsTh
      t8658 = t3429*t8657
      t8659 = t8645*t8658
      t8661 = t3767*t62
      t8662 = t3429*t8661
      t8663 = t8574*t8662
      t8665 = t4678*Mzsq
      t8666 = t3429*t8665
      t8667 = t8645*t8666
      t8670 = t8527*t3606*CsTh
      t8671 = t8523*t8670
      t8673 = t1812*t5
      t8674 = t8673*t8517
      t8675 = t8514*t8674
      t8676 = t4197*t7067
      t8677 = t8675*t8676
      t8679 = t1812*t8517
      t8680 = t8514*t8679
      t8681 = t3606*t6254
      t8682 = t8680*t8681
      t8684 = t8548*t8517
      t8685 = t8514*t8684
      t8686 = t8685*t8676
      t8688 = t3606*t7067
      t8689 = t8675*t8688
      t8691 = -2.D0/15.D0*t8646-t8648/6+t8651/60+t8655/48-t8659/30-t8663
     #/24+t8667/60-t8671/120-t8677/15-t8682/60+t8686/15+t8689/15
      t8694 = t8514*t3388
      t8695 = t4678*t62
      t8696 = t8518*t8695
      t8697 = t8694*t8696
      t8699 = t8518*t8665
      t8700 = t8638*t8699
      t8702 = t8685*t8688
      t8704 = t3429*t8695
      t8705 = t8574*t8704
      t8707 = t8519*CsTh
      t8708 = t3429*t8707
      t8709 = t8563*t8708
      t8711 = t3429*t8519
      t8712 = t8563*t8711
      t8714 = t8563*t8681
      t8716 = t3606*t1096
      t8717 = t8645*t8716
      t8719 = t8544*t8696
      t8721 = t3388*t8517
      t8722 = t8514*t8721
      t8723 = t8722*t8662
      t8725 = t5930*t8517
      t8726 = t8514*t8725
      t8727 = t8726*t4751
      t8729 = t8680*t8708
      t8731 = -t8697/15-t8700/12-t8702/15+t8705/48-t8709/60+t8712/120+t8
     #714/60+t8717/30+t8719/24+t8723/15-t8727/15+t8729/60
      t8732 = t8514*t5930
      t8733 = t8518*t4570
      t8734 = t8732*t8733
      t8736 = t8726*t4778
      t8738 = t8722*t8704
      t8740 = t4678*t1811
      t8741 = t3429*t8740
      t8742 = t8563*t8741
      t8744 = t1067*t8517
      t8745 = t8514*t8744
      t8746 = t8745*t8666
      t8748 = t3606*t6258
      t8749 = t8680*t8748
      t8751 = t4197*t7216
      t8752 = t8685*t8751
      t8754 = t8518*t4679
      t8755 = t8732*t8754
      t8757 = t8680*t8741
      t8759 = t8518*t8661
      t8760 = t8694*t8759
      t8762 = t8518*t8657
      t8763 = t8638*t8762
      t8765 = t8645*t4781
      t8767 = t8734/15+t8736/30-t8738/30+t8742/120-t8746/24+t8749/30-2.D
     #0/15.D0*t8752-t8755/30-t8757/120+2.D0/15.D0*t8760+t8763/6+2.D0/15.
     #D0*t8765
      t8769 = t8675*t8751
      t8771 = t8645*t4751
      t8773 = t3606*t7216
      t8774 = t8675*t8773
      t8776 = t8563*t8748
      t8778 = t3606*t720
      t8779 = t8645*t8778
      t8781 = t8544*t8759
      t8783 = t8645*t4606
      t8785 = t8745*t8658
      t8787 = t8558*t3409
      t8788 = t3623*t8787
      t8789 = t8579*t8788
      t8791 = t8694*t8654
      t8793 = t8685*t8773
      t8795 = t8517*t3406
      t8796 = t8795*t4534
      t8797 = t8732*t8796
      t8799 = 2.D0/15.D0*t8769+4.D0/15.D0*t8771-2.D0/15.D0*t8774-t8776/3
     #0-t8779/15-t8781/12-4.D0/15.D0*t8783+t8785/12+t8789/6-t8791/30+2.D
     #0/15.D0*t8793+t8797/30
      t8800 = t8680*t8711
      t8803 = t8514*t8673*t8551
      t8805 = t3623*t8550
      t8806 = t8675*t8805
      t8808 = t8645*t4633
      t8810 = t8694*t8559
      t8812 = t8685*t8805
      t8814 = t8638*t8622
      t8816 = t8732*t8625
      t8818 = t8563*t8788
      t8820 = t3623*t3884
      t8821 = t8574*t8820
      t8823 = t8555*t8820
      t8825 = t8518*t3420
      t8826 = t8609*t8825
      t8828 = -t8800/120+t8803/15-t8806/15-2.D0/15.D0*t8808-t8810/15+t88
     #12/15-t8814/12-t8816/30-4.D0/15.D0*t8818-t8821/3+2.D0/15.D0*t8823-
     #t8826/30
      wI59pc1 = t8567+t8607+t8643+t8691+t8731+t8767+t8799+t8828
      t8831 = t211*wV7
      t8833 = t8831*t8554*t3426
      t8834 = t3842*t8695
      t8835 = t8833*t8834
      t8837 = t8831*t8534
      t8838 = t8837*t4798
      t8840 = t8831*t8527
      t8841 = t3893*t8519
      t8842 = t8840*t8841
      t8845 = t8840*t3893*t8707
      t8847 = t8831*t8554
      t8848 = t8847*t8704
      t8850 = t8837*t4735
      t8852 = t8831*t1070
      t8853 = t8852*t8531
      t8855 = t8847*t8834
      t8857 = t8837*t4648
      t8859 = t8831*t8517
      t8860 = t8859*t8564
      t8862 = t8840*t8748
      t8864 = t8840*t8681
      t8868 = t8847*t3429*t3643*t62
      t8870 = -t8835/8000+t8838/2000+6.D0/625.D0*t8842-12.D0/625.D0*t884
     #5-207.D0/20000.D0*t8848+13.D0/5000.D0*t8850+9.D0/10000.D0*t8853-t8
     #855/4000-13.D0/10000.D0*t8857+207.D0/2500.D0*t8860+27.D0/2500.D0*t
     #8862-27.D0/5000.D0*t8864+23.D0/10000.D0*t8868
      t8871 = t3751*t8519
      t8872 = t8840*t8871
      t8874 = t8831*t8518
      t8876 = t8874*t3619*t3999
      t8878 = t8831*t8644
      t8879 = t8878*t8658
      t8882 = t8852*t8518*t3613
      t8884 = t8831*t8795
      t8886 = t8884*t3780*t8787
      t8888 = t8878*t4600
      t8890 = t8831*s
      t8891 = t8890*t8559
      t8895 = t8847*t3429*t3613*t62
      t8897 = t8859*t8841
      t8899 = t8847*t3894
      t8901 = t5*t8517
      t8902 = t8901*t3428
      t8903 = t8831*t8902
      t8905 = t8903*t3613*t7216
      t8907 = t8831*t8568
      t8908 = t8907*t3897
      t8910 = t8847*t3897
      t8913 = t8852*t8518*t3634
      t8915 = t8872/1000+52.D0/625.D0*t8876+9.D0/1000.D0*t8879-13.D0/500
     #0.D0*t8882-52.D0/625.D0*t8886+3.D0/625.D0*t8888-207.D0/10000.D0*t8
     #891+179.D0/10000.D0*t8895+12.D0/625.D0*t8897-7.D0/125.D0*t8899+96.
     #D0/625.D0*t8905-26.D0/625.D0*t8908+26.D0/625.D0*t8910-t8913/4000
      t8917 = t8907*t3894
      t8919 = t8907*t3925
      t8922 = t8884*t3842*t8787
      t8924 = t8859*t8681
      t8926 = t8884*t8564
      t8929 = t8878*t3842*t3883
      t8932 = t8874*t3643*t3999
      t8935 = t8884*t3751*t8787
      t8938 = t8847*t3751*t8558
      t8940 = t8837*t4855
      t8943 = t8831*t8795*t3428
      t8944 = t3634*t4040
      t8945 = t8943*t8944
      t8948 = t8878*t3751*t3883
      t8952 = t8878*t3429*t3613*Mzsq
      t8954 = 7.D0/125.D0*t8917-t8919/500-t8922/1000-27.D0/2500.D0*t8924
     #-207.D0/2500.D0*t8926-t8929/2000-23.D0/1250.D0*t8932+23.D0/1250.D0
     #*t8935+23.D0/5000.D0*t8938+3.D0/10000.D0*t8940+t8945/500+2.D0/625.
     #D0*t8948+7.D0/1000.D0*t8952
      t8955 = t8878*t4731
      t8957 = t3780*t8519
      t8958 = t8840*t8957
      t8962 = t8847*t3429*t3619*t62
      t8965 = t8907*t3650*t3795
      t8967 = t3634*t7022
      t8968 = t8943*t8967
      t8971 = t8878*t3893*t3883
      t8973 = t8878*t4633
      t8976 = t8907*t3650*t4025
      t8978 = t8878*t4772
      t8980 = t8554*t3428
      t8981 = t8831*t8980
      t8982 = t8981*t3870
      t8985 = t8981*t3634*t4025
      t8988 = t8831*t8901*t3406
      t8989 = t3780*t8550
      t8990 = t8988*t8989
      t8992 = t8831*t8901
      t8993 = t8992*t8989
      t8996 = t8847*t3780*t8558
      t8998 = -t8955/1000-13.D0/2500.D0*t8958-13.D0/1250.D0*t8962+t8965/
     #250-t8968/1000+7.D0/500.D0*t8971+9.D0/2500.D0*t8973-t8976/500+3.D0
     #/1250.D0*t8978-t8982/250+t8985/500-26.D0/625.D0*t8990+26.D0/625.D0
     #*t8993-13.D0/625.D0*t8996
      t9002 = t8847*t3893*t8558
      t9004 = t3893*t8550
      t9005 = t8988*t9004
      t9007 = t8874*t8944
      t9009 = t8859*t8748
      t9011 = t8874*t8967
      t9013 = t8878*t4808
      t9015 = t8878*t4582
      t9017 = t8878*t4549
      t9020 = t8988*t3686*t7216
      t9023 = t8847*t3842*t8558
      t9025 = t8837*t4722
      t9027 = t3751*t8550
      t9028 = t8992*t9027
      t9031 = t8988*t3686*t7067
      t9033 = 179.D0/5000.D0*t9002+48.D0/625.D0*t9005-t9007/500+27.D0/12
     #50.D0*t9009+t9011/1000-13.D0/1250.D0*t9013+t9015/1000+13.D0/1250.D
     #0*t9017-2.D0/125.D0*t9020-t9023/4000-13.D0/10000.D0*t9025-t9028/12
     #5+t9031/125
      t9034 = t8981*t3825
      t9036 = t8988*t8688
      t9038 = t8831*t54
      t9039 = t9038*t8622
      t9041 = t8988*t8773
      t9044 = t8981*t3619*t4025
      t9046 = t8992*t8688
      t9048 = t8992*t8773
      t9050 = t8859*t8711
      t9053 = t8907*t3677*t3795
      t9056 = t8874*t4678*t3999
      t9058 = t3613*t7022
      t9059 = t8874*t9058
      t9062 = t8874*t3643*t6258
      t9064 = t8847*t3925
      t9067 = t8981*t3643*t4025
      t9069 = -52.D0/625.D0*t9034-27.D0/625.D0*t9036-9.D0/1000.D0*t9039+
     #54.D0/625.D0*t9041+26.D0/625.D0*t9044+27.D0/625.D0*t9046-54.D0/625
     #.D0*t9048-27.D0/5000.D0*t9050+52.D0/625.D0*t9053+207.D0/2500.D0*t9
     #056-179.D0/1250.D0*t9059-t9062/250+t9064/500-8.D0/625.D0*t9067
      t9072 = t8907*t3677*t4025
      t9074 = t3643*t4040
      t9075 = t8874*t9074
      t9078 = t8874*t3643*t6254
      t9081 = t8874*t3634*t3999
      t9083 = t8943*t9074
      t9085 = t3643*t7022
      t9086 = t8874*t9085
      t9088 = t8847*t8570
      t9090 = t8907*t3812
      t9092 = t8988*t9027
      t9094 = t8859*t8957
      t9096 = t8943*t9085
      t9098 = t1811*t3409
      t9100 = t8992*t3606*t9098
      t9102 = t8837*t4787
      t9104 = -26.D0/625.D0*t9072+23.D0/625.D0*t9075+t9078/500+t9081/100
     #0-23.D0/625.D0*t9083-23.D0/1250.D0*t9086+9.D0/250.D0*t9088+9.D0/12
     #5.D0*t9090+t9092/125-13.D0/1250.D0*t9094+23.D0/1250.D0*t9096+27.D0
     #/625.D0*t9100-13.D0/5000.D0*t9102
      t9105 = t8992*t8805
      t9107 = t8837*t4719
      t9110 = t8907*t3686*t3795
      t9112 = t3613*t4040
      t9113 = t8943*t9112
      t9115 = t8884*t8602
      t9118 = t8878*t3780*t3883
      t9120 = t3619*t4040
      t9121 = t8943*t9120
      t9123 = t8981*t3815
      t9126 = t8907*t3682*t3795
      t9129 = t8884*t3893*t8787
      t9131 = t8874*t9112
      t9133 = t8874*t9120
      t9136 = t8874*t3767*t3999
      t9139 = t8831*t8644*t3426
      t9140 = t3780*t8657
      t9141 = t9139*t9140
      t9143 = -27.D0/625.D0*t9105+13.D0/2500.D0*t9107-16.D0/625.D0*t9110
     #-179.D0/625.D0*t9113+207.D0/1250.D0*t9115-13.D0/1250.D0*t9118+104.
     #D0/625.D0*t9121+14.D0/125.D0*t9123-14.D0/125.D0*t9126+179.D0/1250.
     #D0*t9129+179.D0/625.D0*t9131-104.D0/625.D0*t9133-207.D0/1250.D0*t9
     #136+13.D0/1250.D0*t9141
      t9147 = t3751*t8661
      t9148 = t8833*t9147
      t9150 = t8907*t3885
      t9152 = t3780*t8665
      t9153 = t9139*t9152
      t9155 = t8847*t8820
      t9157 = t3751*t8695
      t9158 = t8833*t9157
      t9160 = t8837*t4762
      t9162 = t8837*t4742
      t9164 = t8837*t4683
      t9166 = t8837*t4739
      t9169 = t8840*t3751*t8707
      t9171 = t9038*t8625
      t9173 = t8837*t4879
      t9175 = t8837*t4882
      t9177 = -23.D0/5000.D0*t9148+8.D0/625.D0*t9150-13.D0/2500.D0*t9153
     #-9.D0/250.D0*t9155+23.D0/10000.D0*t9158-t9160/8000-9.D0/10000.D0*t
     #9162-t9164/4000-3.D0/5000.D0*t9166-t9169/500-9.D0/2500.D0*t9171+3.
     #D0/10000.D0*t9173+9.D0/20000.D0*t9175
      t9178 = t8847*t3818
      t9181 = t8840*t3751*t8740
      t9183 = t8518*t3643
      t9184 = t8852*t9183
      t9186 = t8847*t9147
      t9188 = t8847*t9157
      t9190 = t8878*t4609
      t9192 = t8878*t4775
      t9195 = t8831*t8644*t3406
      t9196 = t9195*t4609
      t9198 = t9195*t4775
      t9201 = t8903*t3643*t7216
      t9204 = t8903*t3643*t7067
      t9206 = t3842*t8657
      t9207 = t8878*t9206
      t9209 = t8852*t8528
      t9211 = t8878*t9140
      t9213 = 9.D0/250.D0*t9178+t9181/1000+3.D0/5000.D0*t9184-23.D0/2500
     #.D0*t9186+23.D0/5000.D0*t9188-13.D0/625.D0*t9190+13.D0/1250.D0*t91
     #92+13.D0/625.D0*t9196-13.D0/1250.D0*t9198+2.D0/125.D0*t9201-t9204/
     #125+t9207/1000+9.D0/20000.D0*t9209+13.D0/625.D0*t9211
      t9215 = t8878*t9152
      t9217 = t8837*t4705
      t9219 = t3842*t8665
      t9220 = t8878*t9219
      t9222 = t8859*t8871
      t9224 = t3842*t8661
      t9225 = t8847*t9224
      t9227 = t8852*t8524
      t9229 = t8878*t4528
      t9231 = t8847*t3885
      t9233 = t8981*t3796
      t9235 = t9195*t4784
      t9237 = t8878*t8650
      t9241 = t8878*t3429*t3619*Mzsq
      t9243 = t8847*t3812
      t9245 = -13.D0/1250.D0*t9215-t9217/8000-t9220/2000+t9222/500+t9225
     #/2000-9.D0/5000.D0*t9227-3.D0/1250.D0*t9229-8.D0/625.D0*t9231+16.D
     #0/625.D0*t9233+3.D0/1250.D0*t9235-9.D0/2000.D0*t9237-13.D0/2500.D0
     #*t9241-9.D0/125.D0*t9243
      t9248 = t8878*t3429*t3643*Mzsq
      t9250 = t8847*t8662
      t9252 = t3751*t8665
      t9253 = t9139*t9252
      t9255 = t3751*t8657
      t9256 = t9139*t9255
      t9259 = t8847*t3429*t8558
      t9261 = t8837*t4852
      t9263 = t8878*t4778
      t9265 = t8878*t4781
      t9268 = t8988*t3677*t7067
      t9270 = t8852*t8825
      t9272 = t3780*t8695
      t9273 = t8847*t9272
      t9276 = t8988*t3682*t7067
      t9279 = t8859*t3606*t1811
      t9282 = t8903*t3619*t7067
      t9284 = t9248/625+207.D0/10000.D0*t9250+t9253/625-2.D0/625.D0*t925
     #6-207.D0/20000.D0*t9259+t9261/4000+9.D0/2500.D0*t9263-9.D0/2500.D0
     #*t9265-26.D0/625.D0*t9268+9.D0/10000.D0*t9270-13.D0/625.D0*t9273+4
     #8.D0/625.D0*t9276-27.D0/2500.D0*t9279+26.D0/625.D0*t9282
      t9287 = t3893*t8665
      t9288 = t8878*t9287
      t9290 = t8837*t4708
      t9292 = t3893*t8661
      t9293 = t8833*t9292
      t9295 = t9195*t4571
      t9297 = t3893*t8657
      t9298 = t9139*t9297
      t9301 = t8840*t3780*t8707
      t9303 = t8992*t9004
      t9306 = t8874*t3613*t3999
      t9308 = t8907*t8570
      t9311 = t8907*t3686*t4025
      t9313 = t8943*t9058
      t9316 = t8907*t3682*t4025
      t9318 = t8884*t8580
      t9320 = 7.D0/500.D0*t9288-3.D0/2500.D0*t9290-179.D0/5000.D0*t9293+
     #t9295/500-7.D0/500.D0*t9298+13.D0/1250.D0*t9301-48.D0/625.D0*t9303
     #-179.D0/1250.D0*t9306-9.D0/250.D0*t9308+8.D0/625.D0*t9311+179.D0/1
     #250.D0*t9313+7.D0/125.D0*t9316-207.D0/2500.D0*t9318
      t9321 = t3619*t7022
      t9322 = t8943*t9321
      t9325 = t8981*t3613*t4025
      t9327 = t8874*t9321
      t9329 = t8878*t4571
      t9331 = t8878*t8666
      t9334 = t8840*t3893*t8740
      t9336 = t9195*t4680
      t9338 = t3893*t8695
      t9339 = t8833*t9338
      t9341 = t9139*t9287
      t9343 = t8837*t4795
      t9346 = t8840*t3780*t8740
      t9348 = t8847*t9338
      t9350 = t8878*t8716
      t9352 = t3780*t8661
      t9353 = t8833*t9352
      t9355 = -52.D0/625.D0*t9322-7.D0/125.D0*t9325+52.D0/625.D0*t9327-t
     #9329/500-9.D0/2000.D0*t9331+6.D0/625.D0*t9334-t9336/1000+179.D0/10
     #000.D0*t9339+7.D0/1000.D0*t9341+3.D0/5000.D0*t9343-13.D0/2500.D0*t
     #9346+179.D0/5000.D0*t9348-9.D0/1000.D0*t9350+13.D0/625.D0*t9353
      t9358 = t8874*t3619*t6254
      t9360 = t9139*t9206
      t9362 = t9139*t9219
      t9364 = t8833*t9224
      t9368 = t8878*t3429*t3634*Mzsq
      t9371 = t8847*t3606*t5992
      t9373 = t8878*t4680
      t9375 = t8833*t9272
      t9378 = t8874*t3619*t6258
      t9380 = t8847*t9352
      t9383 = t8988*t3682*t7216
      t9385 = t8847*t9292
      t9387 = t8878*t8778
      t9389 = -13.D0/1250.D0*t9358+t9360/2000-t9362/4000+t9364/4000-t936
     #8/4000-207.D0/10000.D0*t9371+t9373/1000-13.D0/1250.D0*t9375+13.D0/
     #625.D0*t9378+26.D0/625.D0*t9380-96.D0/625.D0*t9383-179.D0/2500.D0*
     #t9385+9.D0/500.D0*t9387
      t9392 = t8847*t3429*t3634*t62
      t9394 = t8878*t4751
      t9396 = t8878*t4784
      t9398 = t8878*t4606
      t9400 = t9195*t4600
      t9403 = t8903*t3613*t7067
      t9406 = t8874*t3613*t6254
      t9409 = t8903*t3619*t7216
      t9412 = t8988*t3677*t7216
      t9415 = t8874*t3613*t6258
      t9417 = t8878*t9252
      t9420 = t8847*t3606*t1093
      t9422 = t8878*t9255
      t9424 = t8878*t9297
      t9426 = -t9392/8000-9.D0/1250.D0*t9394-3.D0/1250.D0*t9396+9.D0/125
     #0.D0*t9398-3.D0/625.D0*t9400-48.D0/625.D0*t9403+12.D0/625.D0*t9406
     #-52.D0/625.D0*t9409+52.D0/625.D0*t9412-24.D0/625.D0*t9415+2.D0/625
     #.D0*t9417+207.D0/5000.D0*t9420-4.D0/625.D0*t9422-7.D0/250.D0*t9424
      wI60pc1 = t8870+t8915+t8954+t8998+t9033+t9069+t9104+t9143+t9177+t9
     #213+t9245+t9284+t9320+t9355+t9389+t9426
      t9430 = t211*wV8
      t9431 = t1813*t9430
      t9433 = 1/(t1072-t5874-t1077+t5857)
      t9434 = s*t9433
      t9436 = 1/(-cw2+t5861)
      t9437 = t9436*t3420
      t9438 = t9437*CsTh
      t9439 = t9434*t9438
      t9440 = t9431*t9439
      t9442 = t1834*t9430
      t9443 = t9433*t5862
      t9444 = t9437*t842
      t9445 = t9443*t9444
      t9446 = t9442*t9445
      t9448 = t3559*t9430
      t9450 = t5862*t9436*t3420
      t9451 = t9434*t9450
      t9452 = t9448*t9451
      t9454 = t3393*t9430
      t9455 = t9433*t9436
      t9456 = t3767*t5851
      t9457 = t9455*t9456
      t9458 = t9454*t9457
      t9460 = t4678*t5851
      t9461 = t9455*t9460
      t9462 = t9454*t9461
      t9464 = t1832*t1812
      t9465 = t9464*t9430
      t9466 = t9455*t3767
      t9467 = t9465*t9466
      t9469 = t9455*t4678
      t9470 = t9465*t9469
      t9472 = t7*t4893
      t9475 = t9455*t3420*t5861*t5851
      t9476 = t9472*t9475
      t9478 = t9430*t5
      t9479 = t3542*t9478
      t9480 = t9455*t3420
      t9482 = t5861*CsTh*t5851
      t9483 = t9480*t9482
      t9484 = t9479*t9483
      t9486 = t3542*t9430
      t9487 = t9486*t9439
      t9489 = t9434*t9444
      t9490 = t9486*t9489
      t9493 = t5861*t842*t5851
      t9494 = t9480*t9493
      t9495 = t9479*t9494
      t9497 = t9442*t9466
      t9499 = t1813*t9478
      t9500 = t9499*t9483
      t9502 = t9499*t9494
      t9504 = -13.D0/100.D0*t9440-t9446/80+13.D0/400.D0*t9452+13.D0/25.D
     #0*t9458-13.D0/50.D0*t9462+t9467/20-t9470/40+3.D0/25.D0*t9476+t9484
     #/5+t9487/10-t9490/20-t9495/10-t9497/20-t9500/5+t9502/10
      t9505 = t178*t9430
      t9506 = t9505*t9475
      t9508 = t9430*t9433
      t9509 = t3386*t9508
      t9510 = t9437*t9482
      t9511 = t9509*t9510
      t9513 = t9431*t9451
      t9516 = t1832*t1845*t9430
      t9517 = t9443*t9438
      t9518 = t9516*t9517
      t9520 = t9516*t9445
      t9522 = t9442*t9469
      t9524 = t9437*t9493
      t9525 = t9509*t9524
      t9527 = t54*t9433
      t9529 = t9505*t9527*t9437
      t9531 = t9505*t9461
      t9533 = t9472*t9494
      t9535 = t3386*t9430
      t9536 = t9535*t9475
      t9538 = t9472*t9483
      t9540 = t9479*t9475
      t9542 = t9499*t9475
      t9544 = t177*t211
      t9545 = wV8*s
      t9547 = t9544*t9545*t9433
      t9548 = Mzsq*t9436
      t9550 = t9547*t9548*t9456
      t9552 = 13.D0/50.D0*t9506+2.D0/5.D0*t9511-t9513/40-t9518/40+t9520/
     #80+t9522/40-t9525/5+3.D0/100.D0*t9529+t9531/5+3.D0/25.D0*t9533-t95
     #36/5-6.D0/25.D0*t9538-t9540/10+t9542/10+6.D0/25.D0*t9550
      t9554 = t178*t9508
      t9555 = t9554*t9510
      t9557 = t9554*t9524
      t9560 = t9547*t9548*t9460
      t9563 = t9437*t5851
      t9564 = t5*t9433*t9563
      t9565 = t9431*t9564
      t9567 = t1813*t4893
      t9568 = t9567*t9445
      t9570 = t3393*t4887
      t9571 = t9570*t9517
      t9573 = t9431*t9489
      t9575 = t9499*t9461
      t9577 = t9570*t9445
      t9579 = t9567*t9517
      t9581 = t9443*t9437
      t9582 = t9442*t9581
      t9584 = t9505*t9457
      t9586 = t3559*t4893
      t9587 = t9586*t9517
      t9589 = t9586*t9445
      t9591 = t9499*t9457
      t9593 = -13.D0/25.D0*t9555+13.D0/50.D0*t9557-3.D0/25.D0*t9560+t956
     #5/10-t9568/40-3.D0/100.D0*t9571+13.D0/200.D0*t9573+t9575/10+3.D0/2
     #00.D0*t9577+t9579/20-t9582/80-2.D0/5.D0*t9584-13.D0/200.D0*t9587+1
     #3.D0/400.D0*t9589-t9591/5
      t9595 = t9455*t3420*t5851
      t9596 = t9454*t9595
      t9599 = t9505*t9527*t9444
      t9601 = t3559*t9478
      t9602 = t9601*t9457
      t9605 = t9454*t9527*t9450
      t9607 = t9601*t9461
      t9609 = t9516*t9581
      t9612 = t9433*Mzsq
      t9614 = t9544*t9545*t9612*t9563
      t9616 = t9464*t211
      t9618 = wV8*t9433*t9437
      t9619 = t9616*t9618
      t9622 = t1834*t211*t9618
      t9625 = t9505*t9527*t9438
      t9627 = t9434*t9437
      t9628 = t9431*t9627
      t9630 = t9448*t9564
      t9632 = t9505*t9595
      t9634 = t9442*t9517
      t9636 = t9486*t9627
      t9638 = -13.D0/50.D0*t9596+3.D0/100.D0*t9599+t9602/5+3.D0/200.D0*t
     #9605-t9607/10+t9609/80-3.D0/25.D0*t9614-t9619/40+t9622/40-3.D0/50.
     #D0*t9625+13.D0/200.D0*t9628-t9630/10+t9632/5+t9634/40-t9636/20
      wI61pc1 = t9504+t9552+t9593+t9638
      spf14 = log(DCMPLX(-t1/t165))
      t9642 = t4892*t3613
      t9643 = t1093*spf14
      t9644 = t6943*t9643
      t9645 = t9642*t9644
      t9647 = t211*t3420
      t9648 = t9647*t3842
      t9649 = t7454*wC3
      t9650 = t6953*t9649
      t9651 = t9648*t9650
      t9653 = t4892*t3643
      t9654 = t60*wC3
      t9655 = t6943*t9654
      t9658 = t59*wC3
      t9659 = t6943*t9658
      t9663 = t9647*t3842*t3406
      t9664 = t9663*t9650
      t9666 = t9647*t3780
      t9667 = t6967*t6258
      t9668 = t9666*t9667
      t9670 = t6961*t9649
      t9671 = t9648*t9670
      t9673 = t4892*t3634
      t9674 = t7465*wC3
      t9675 = t6943*t9674
      t9676 = t9673*t9675
      t9678 = t4892*t3686
      t9679 = t6967*t9674
      t9682 = t211*t5
      t9683 = t9682*t3613
      t9684 = t6943*t9098
      t9685 = t9683*t9684
      t9687 = t3999*wC3
      t9688 = t6961*t9687
      t9689 = t9648*t9688
      t9691 = t9647*t3751
      t9692 = t3609*wC3
      t9693 = t7423*t9692
      t9696 = t9682*t3643
      t9697 = t9696*t7170
      t9700 = t3429*t7408*wC3
      t9701 = t9673*t9700
      t9703 = t9682*t3606
      t9705 = t9703*t6961*t9098
      t9707 = -163.D0/5100.D0*t9645+t9651/1700+3.D0/3400.D0*t9653*t9655-
     #3.D0/6800.D0*t9653*t9659-t9664/1700+t9668/7650+t9671/3400-t9676/68
     #00+3.D0/6800.D0*t9678*t9679+14.D0/1275.D0*t9685+t9689/3400-3.D0/85
     #0.D0*t9691*t9693-2.D0/765.D0*t9697+t9701/27200+6.D0/425.D0*t9705
      t9709 = t3429*t6987*wC3
      t9710 = t9673*t9709
      t9712 = t9691*t9667
      t9714 = t9647*t168
      t9715 = t6943*t1811
      t9716 = t9714*t9715
      t9718 = t4892*t3606
      t9719 = Mzsq*spf14
      t9720 = t7103*t9719
      t9721 = t9718*t9720
      t9723 = t3609*spf14
      t9724 = t6943*t9723
      t9725 = t9642*t9724
      t9727 = t4892*t3677
      t9728 = t6961*t9723
      t9729 = t9727*t9728
      t9731 = t4892*t3619
      t9732 = t9731*t9724
      t9734 = t3999*spf14
      t9735 = t6961*t9734
      t9736 = t9691*t9735
      t9738 = t9647*t3623
      t9739 = t3409*spf14
      t9740 = t6987*t9739
      t9741 = t9738*t9740
      t9743 = t3426*t3420
      t9744 = t5590*t9743
      t9745 = t9744*t7638
      t9747 = t4892*t3682
      t9748 = t9747*t9728
      t9750 = t9718*t9728
      t9752 = t9653*t9644
      t9754 = t6953*t9734
      t9755 = t9648*t9754
      t9757 = t62*spf14
      t9758 = t6953*t9757
      t9759 = t9718*t9758
      t9761 = t4468*t3606
      t9762 = t6953*t9739
      t9763 = t9761*t9762
      t9765 = t9710/27200+t9712/1530-t9716/1530+271.D0/5100.D0*t9721-182
     #.D0/3825.D0*t9725-26.D0/1275.D0*t9729+26.D0/1275.D0*t9732-8.D0/382
     #5.D0*t9736-239.D0/7650.D0*t9741+329.D0/61200.D0*t9745+182.D0/3825.
     #D0*t9748-271.D0/5100.D0*t9750+4.D0/3825.D0*t9752+t9755/1530+239.D0
     #/15300.D0*t9759-t9763/450
      t9767 = t9682*t3686
      t9768 = t9098*spf14
      t9769 = t6953*t9768
      t9770 = t9767*t9769
      t9772 = t4468*t3643
      t9773 = t3692*spf14
      t9774 = t6943*t9773
      t9775 = t9772*t9774
      t9778 = t9647*t3751*t3406
      t9781 = t9647*t6943
      t9782 = t3409*wC3
      t9783 = t1093*t9782
      t9784 = t9781*t9783
      t9786 = t6967*t9654
      t9787 = t9718*t9786
      t9789 = t62*wC3
      t9790 = t6943*t9789
      t9793 = t1*t167
      t9794 = log(DCMPLX(t9793))
      t9795 = t3692*t9794
      t9796 = t6961*t9795
      t9797 = t9761*t9796
      t9799 = t9682*t3677
      t9800 = t6961*t7067
      t9801 = t9799*t9800
      t9803 = t9682*t3619
      t9804 = t9803*t7071
      t9806 = t9098*t9794
      t9807 = t6953*t9806
      t9808 = t9767*t9807
      t9810 = t1811*t9794
      t9811 = t6949*t9810
      t9812 = t9691*t9811
      t9814 = t4468*t3686
      t9815 = t9814*t9796
      t9817 = t62*t9794
      t9818 = t6949*t9817
      t9821 = t3999*t9794
      t9822 = t6949*t9821
      t9823 = t9648*t9822
      t9825 = t1096*t9794
      t9826 = t6943*t9825
      t9827 = t9772*t9826
      t9829 = t6943*t9795
      t9830 = t9772*t9829
      t9832 = -4.D0/765.D0*t9770+2.D0/153.D0*t9775+3.D0/850.D0*t9778*t96
     #50+2.D0/425.D0*t9784-t9787/1700-3.D0/6800.D0*t9653*t9790+2.D0/425.
     #D0*t9797-2.D0/3825.D0*t9801+2.D0/3825.D0*t9804+4.D0/765.D0*t9808+t
     #9812/1530+4.D0/1275.D0*t9815-11.D0/3400.D0*t9718*t9818+t9823/3060-
     #59.D0/61200.D0*t9827-4.D0/1275.D0*t9830
      t9833 = t4892*t3650
      t9834 = t3409*t9794
      t9835 = t9834*Mzsq
      t9836 = t6953*t9835
      t9837 = t9833*t9836
      t9839 = Mzsq*t9794
      t9840 = t6949*t9839
      t9841 = t9761*t9840
      t9843 = t9663*t9822
      t9845 = t5992*t9794
      t9846 = t6943*t9845
      t9847 = t9653*t9846
      t9849 = t4468*t3619
      t9850 = t9849*t9826
      t9852 = t4468*t3613
      t9853 = t9852*t9826
      t9855 = t9691*t7019
      t9857 = t9647*t3893
      t9858 = t9857*t9822
      t9860 = t3893*t3406
      t9861 = t9647*t9860
      t9862 = t9861*t9822
      t9864 = t9731*t9846
      t9866 = t9642*t9846
      t9868 = t9691*t9822
      t9870 = t6949*t9835
      t9871 = t9718*t9870
      t9873 = t9691*t7903
      t9877 = t5992*t9834
      t9880 = -13.D0/15300.D0*t9837+73.D0/10200.D0*t9841-t9843/3060+13.D
     #0/15300.D0*t9847+26.D0/3825.D0*t9850-41.D0/3060.D0*t9853-2.D0/1275
     #.D0*t9855-83.D0/3825.D0*t9858+83.D0/3825.D0*t9862-13.D0/3825.D0*t9
     #864+83.D0/15300.D0*t9866-13.D0/3825.D0*t9868-73.D0/2550.D0*t9871-3
     #29.D0/15300.D0*t9873-11.D0/850.D0*t9738*t9822+11.D0/850.D0*t9781*t
     #9877
      t9883 = t9738*t7547
      t9885 = t9691*t7900
      t9887 = t3406*t3420
      t9888 = t5590*t9887
      t9889 = t9888*t7550
      t9891 = t9778*t9822
      t9893 = t9666*t9822
      t9895 = t9857*t7903
      t9897 = t3780*t3406
      t9898 = t9647*t9897
      t9899 = t9898*t9822
      t9901 = t9747*t9870
      t9903 = t9666*t9740
      t9905 = t9682*t3682
      t9906 = t9905*t9769
      t9908 = t9678*t9728
      t9910 = t7588*t9835
      t9911 = t9747*t9910
      t9913 = t9682*t4197
      t9914 = t9913*t9769
      t9916 = t9703*t9769
      t9918 = t7479*t9768
      t9919 = t9905*t9918
      t9921 = -t9883/3825+t9885/2550-106.D0/3825.D0*t9889+13.D0/3825.D0*
     #t9891+52.D0/3825.D0*t9893-8.D0/255.D0*t9895-52.D0/3825.D0*t9899+41
     #.D0/765.D0*t9901+5.D0/153.D0*t9903-58.D0/1275.D0*t9906+127.D0/5100
     #.D0*t9908-41.D0/765.D0*t9911-12.D0/425.D0*t9914+12.D0/425.D0*t9916
     #+58.D0/1275.D0*t9919
      t9922 = t9678*t9870
      t9924 = t4892*t4197
      t9925 = t9924*t9870
      t9927 = t9666*t9811
      t9929 = t9857*t9811
      t9931 = t9678*t9910
      t9933 = t1811*spf14
      t9934 = t6953*t9933
      t9935 = t9666*t9934
      t9937 = t9857*t9934
      t9939 = t9767*t9918
      t9941 = t9727*t9870
      t9943 = t9799*t9769
      t9945 = t9663*t9754
      t9947 = t9799*t9918
      t9949 = t9727*t9910
      t9951 = t4460*t3643
      t9953 = t3429*t6949*spf14
      t9954 = t9951*t9953
      t9956 = t4460*t3634
      t9957 = t6937*spf14
      t9958 = t3429*t9957
      t9959 = t9956*t9958
      t9961 = t7588*t9806
      t9962 = t9799*t9961
      t9964 = 59.D0/15300.D0*t9922+73.D0/2550.D0*t9925-13.D0/3825.D0*t99
     #27+8.D0/1275.D0*t9929-59.D0/15300.D0*t9931-43.D0/7650.D0*t9935+29.
     #D0/2550.D0*t9937+4.D0/765.D0*t9939-104.D0/3825.D0*t9941+86.D0/3825
     #.D0*t9943-t9945/1530-86.D0/3825.D0*t9947+104.D0/3825.D0*t9949+t995
     #4/1224-59.D0/122400.D0*t9959-52.D0/3825.D0*t9962
      t9966 = t6949*t9806
      t9967 = t9799*t9966
      t9969 = t9905*t9966
      t9971 = t9691*t9934
      t9973 = t6953*t9723
      t9974 = t9833*t9973
      t9976 = t6943*t9768
      t9977 = t9683*t9976
      t9979 = t6961*t9768
      t9980 = t9703*t9979
      t9982 = t9905*t9979
      t9984 = t6938*t9739
      t9985 = t9703*t9984
      t9987 = t9905*t9961
      t9989 = t720*spf14
      t9990 = t6943*t9989
      t9991 = t9772*t9990
      t9993 = t9803*t9976
      t9995 = t9767*t9979
      t9997 = t9799*t9979
      t9999 = t9767*t9961
      t10001 = t9703*t9966
      t10003 = t9913*t9966
      t10005 = 52.D0/3825.D0*t9967-32.D0/1275.D0*t9969+t9971/765+t9974/3
     #060+29.D0/1275.D0*t9977+6.D0/425.D0*t9980-29.D0/1275.D0*t9982-6.D0
     #/425.D0*t9985+32.D0/1275.D0*t9987-127.D0/10200.D0*t9991-43.D0/3825
     #.D0*t9993-2.D0/765.D0*t9995+43.D0/3825.D0*t9997+2.D0/765.D0*t9999+
     #6.D0/425.D0*t10001-6.D0/425.D0*t10003
      t10006 = t3970*t9794
      t10007 = t6943*t10006
      t10008 = t9852*t10007
      t10010 = t9833*t9910
      t10012 = t6949*t9834
      t10013 = t9761*t10012
      t10015 = t6953*t9719
      t10016 = t9761*t10015
      t10018 = t6949*t9768
      t10019 = t9905*t10018
      t10021 = t4468*t3682
      t10022 = t6961*t10006
      t10023 = t10021*t10022
      t10025 = t9761*t10022
      t10027 = t9647*t6996
      t10028 = t7521*t9692
      t10029 = t10027*t10028
      t10031 = t6949*t9654
      t10032 = t9718*t10031
      t10034 = t9648*t9735
      t10036 = t6943*t9739
      t10037 = t9852*t10036
      t10039 = t9814*t10022
      t10043 = log(DCMPLX(t1/t53))
      t10044 = t3409*t10043
      t10045 = t7423*t10044
      t10046 = t9718*t10045
      t10048 = t7588*t9768
      t10049 = t9905*t10048
      t10051 = t9703*t10018
      t10053 = t7480*t10043
      t10054 = t6953*t10053
      t10057 = -26.D0/3825.D0*t10008-13.D0/30600.D0*t10010+t10013/425+27
     #1.D0/10200.D0*t10016+29.D0/1275.D0*t10019+26.D0/3825.D0*t10023-t10
     #025/425+4.D0/425.D0*t10029-t10032/850+t10034/3060+13.D0/2550.D0*t1
     #0037-2.D0/1275.D0*t10039+t10046/1275-29.D0/1275.D0*t10049-6.D0/425
     #.D0*t10051+t9678*t10054/850
      t10061 = t9718*t10054
      t10063 = t9924*t10054
      t10065 = t211*cw2
      t10066 = t10065*t3420
      t10067 = t10066*t7431
      t10070 = t9647*t3842*t6937
      t10071 = t7648*t9692
      t10072 = t10070*t10071
      t10074 = t7578*t9782
      t10075 = t9781*t10074
      t10078 = t6949*t53*t9692
      t10079 = t9663*t10078
      t10081 = t6949*t9649
      t10084 = t5992*t9782
      t10085 = t9781*t10084
      t10089 = t9738*t10081
      t10091 = t7479*t10053
      t10094 = t6961*t10053
      t10095 = t9718*t10094
      t10097 = t6949*t9687
      t10100 = t9738*t10097
      t10104 = -2.D0/1275.D0*t10061+2.D0/1275.D0*t10063+91.D0/3060.D0*t1
     #0067+t10072/1700-t10075/425-t10079/1700-3.D0/1700.D0*t9778*t10081-
     #t10085/425+3.D0/1700.D0*t9691*t10081+t10089/425-t9678*t10091/850-t
     #10095/1275-3.D0/1700.D0*t9778*t10097+t10100/425+3.D0/1700.D0*t9691
     #*t10097
      t10106 = t9647*t3751*t6937
      t10109 = t9781*t10071
      t10113 = t9647*t3428
      t10114 = t10113*t7059
      t10116 = t6943*t10053
      t10121 = t53*t10043
      t10122 = t6967*t10121
      t10123 = t9761*t10122
      t10125 = t9682*t3420
      t10126 = t10125*t9897
      t10127 = t3999*t147
      t10128 = t8169*t10127
      t10129 = t10126*t10128
      t10131 = t6949*t10053
      t10132 = t9833*t10131
      t10135 = t3429*t7423*t10043
      t10138 = t10027*t10071
      t10140 = t9696*t9684
      t10142 = t9663*t10097
      t10144 = t4460*t3606
      t10145 = CsTh*t9794
      t10146 = t6967*t10145
      t10147 = t10144*t10146
      t10150 = t3429*t6953*t9794
      t10151 = t9951*t10150
      t10153 = t9905*t10128
      t10155 = -3.D0/850.D0*t10106*t10071+2.D0/425.D0*t10109+3.D0/850.D0
     #*t9778*t10078+13.D0/225.D0*t10114-t9653*t10116/1700+t9678*t10094/1
     #700-t10123/10200+64.D0/3825.D0*t10129+t10132/10200+t9772*t10135/13
     #600-2.D0/425.D0*t10138+2.D0/765.D0*t10140+t10142/3400+t10147/1700+
     #t10151/2550+128.D0/3825.D0*t10153
      t10157 = t9648*t10097
      t10159 = t9913*t10128
      t10161 = t10125*t9860
      t10162 = t10161*t10128
      t10164 = t10066*t7100
      t10166 = t6943*t10145
      t10167 = t9956*t10166
      t10169 = t4468*t3650
      t10170 = t6967*t9825
      t10171 = t10169*t10170
      t10173 = t9913*t10018
      t10175 = t6949*t9933
      t10176 = t9666*t10175
      t10178 = t9767*t10048
      t10180 = t4460*t3420
      t10182 = log(DCMPLX(t9793*t168))
      t10184 = t10180*t6943*t10182
      t10186 = t6967*t9845
      t10187 = t9833*t10186
      t10189 = t4460*t3613
      t10191 = t3429*t6949*t9794
      t10192 = t10189*t10191
      t10194 = t842*t9794
      t10195 = t6967*t10194
      t10196 = t10144*t10195
      t10198 = t9951*t10191
      t10200 = t6967*t6254
      t10201 = t9857*t10200
      t10203 = t9666*t10200
      t10205 = -t10157/3400+64.D0/3825.D0*t10159-128.D0/3825.D0*t10162-9
     #1.D0/765.D0*t10164-t10167/3060-13.D0/244800.D0*t10171+6.D0/425.D0*
     #t10173+43.D0/15300.D0*t10176-2.D0/765.D0*t10178+t10184/850-t10187/
     #24480+13.D0/15300.D0*t10192-t10196/3400-t10198/5100-7.D0/5100.D0*t
     #10201-t10203/15300
      t10206 = t6943*t10194
      t10207 = t9956*t10206
      t10209 = t10113*t6983
      t10211 = t4468*t3677
      t10212 = t10211*t10170
      t10214 = t4468*t3634
      t10215 = t10214*t9826
      t10217 = t6937*t9794
      t10218 = t10217*Mzsq
      t10219 = t3429*t10218
      t10220 = t9849*t10219
      t10222 = t720*t9794
      t10223 = t6967*t10222
      t10224 = t9814*t10223
      t10226 = t6943*t9839
      t10227 = t10214*t10226
      t10229 = t3429*t10217
      t10230 = t9956*t10229
      t10232 = t211*t3406
      t10233 = t10232*t3420
      t10234 = t10233*t7100
      t10236 = t10169*t10022
      t10238 = t9761*t10170
      t10240 = t6967*t1811
      t10241 = t10113*t10240
      t10243 = t9951*t10166
      t10245 = t6258*t9794
      t10246 = t6967*t10245
      t10247 = t9691*t10246
      t10249 = t9647*t177
      t10250 = t10249*t9715
      t10252 = t1093*t9794
      t10253 = t6943*t10252
      t10254 = t9673*t10253
      t10256 = t10207/6120-3.D0/425.D0*t10209+13.D0/3825.D0*t10212-13.D0
     #/122400.D0*t10215+13.D0/3825.D0*t10220+59.D0/61200.D0*t10224-13.D0
     #/122400.D0*t10227+t10230/12240-13.D0/225.D0*t10234+t10236/1530+73.
     #D0/20400.D0*t10238+3.D0/1700.D0*t10241+t10243/1275-t10247/1530-t10
     #250/7650+t10254/6120
      t10259 = t10214*t9829
      t10261 = t9647*t3429
      t10262 = t10261*t9811
      t10264 = t9727*t10186
      t10266 = t9666*t10246
      t10268 = t9857*t10246
      t10270 = t9747*t10186
      t10272 = t5598*t3420
      t10273 = t10272*t7572
      t10275 = t9714*t7572
      t10277 = t211*t3426
      t10278 = t10277*t3420
      t10279 = t10278*t7572
      t10281 = t9761*t10223
      t10283 = t10211*t10223
      t10285 = t10169*t9796
      t10287 = t10066*t9715
      t10289 = t10189*t10166
      t10291 = t6943*t10222
      t10292 = t10214*t10291
      t10294 = t10259/765-3.D0/1700.D0*t10262-13.D0/7650.D0*t10264+13.D0
     #/3825.D0*t10266-8.D0/1275.D0*t10268+83.D0/30600.D0*t10270-3.D0/170
     #0.D0*t10273+53.D0/15300.D0*t10275+t10279/61200-73.D0/10200.D0*t102
     #81-26.D0/3825.D0*t10283-t10285/765-7.D0/2550.D0*t10287-13.D0/3825.
     #D0*t10289+13.D0/61200.D0*t10292
      t10295 = t6967*t10252
      t10296 = t9727*t10295
      t10298 = t6953*t9810
      t10299 = t10261*t10298
      t10301 = t10113*t7340
      t10303 = t6953*t9687
      t10306 = t9673*t9790
      t10308 = t9857*t10175
      t10310 = t7521*t10043
      t10311 = t6967*t10310
      t10312 = t9761*t10311
      t10314 = t10066*t7723
      t10316 = t6943*t10310
      t10317 = t10214*t10316
      t10319 = t7103*spf14
      t10320 = t3623*t10319
      t10321 = t9852*t10320
      t10323 = t9738*t9693
      t10325 = t1093*t10182
      t10326 = t6967*t10325
      t10327 = t9718*t10326
      t10332 = t6967*t6258*t10182
      t10333 = t9666*t10332
      t10335 = t842*t10182
      t10336 = t6943*t10335
      t10337 = t9951*t10336
      t10339 = t9647*t7164
      t10340 = t10339*t10332
      t10342 = 13.D0/3825.D0*t10296+3.D0/850.D0*t10299+3.D0/850.D0*t1030
     #1+3.D0/850.D0*t9778*t10303+t10306/13600-29.D0/5100.D0*t10308+t1031
     #2/5100+13.D0/3825.D0*t10314+t10317/20400-13.D0/2550.D0*t10321-2.D0
     #/425.D0*t10323-t10327/200-t9814*t10311/6800-3.D0/425.D0*t10333+t10
     #337/1700+2.D0/425.D0*t10340
      t10344 = t9098*t10182
      t10345 = t6943*t10344
      t10346 = t9683*t10345
      t10348 = t3409*t10182
      t10349 = t1093*t10348
      t10350 = t9781*t10349
      t10352 = t6987*t10182
      t10353 = t3429*t10352
      t10354 = t9731*t10353
      t10356 = t6943*t10348
      t10357 = t9772*t10356
      t10359 = t5992*t10182
      t10360 = t6967*t10359
      t10361 = t9747*t10360
      t10363 = t62*t10182
      t10364 = t6943*t10363
      t10365 = t9642*t10364
      t10368 = t6967*t6254*t10182
      t10369 = t9857*t10368
      t10371 = t1811*t10182
      t10372 = t6967*t10371
      t10373 = t9857*t10372
      t10375 = t10348*Mzsq
      t10376 = t6953*t10375
      t10377 = t9678*t10376
      t10379 = Mzsq*t10182
      t10380 = t6943*t10379
      t10381 = t9772*t10380
      t10385 = t9772*t3623*t7103*t10182
      t10387 = t9852*t10380
      t10389 = t9642*t10353
      t10391 = t9731*t10364
      t10393 = t6953*t10371
      t10394 = t10339*t10393
      t10396 = t9803*t10345
      t10398 = 4.D0/425.D0*t10346+t10350/25-t10354/1020-t10357/425-t1036
     #1/510-t10365/255-t10369/850-t10373/850-13.D0/2550.D0*t10377-13.D0/
     #20400.D0*t10381+t10385/425+7.D0/1700.D0*t10387-t10389/510-t10391/5
     #10+4.D0/425.D0*t10394-12.D0/425.D0*t10396
      t10400 = t9682*t6933*t3428
      t10401 = t6961*t10344
      t10402 = t10400*t10401
      t10405 = t9682*t6933*t10345
      t10407 = t9799*t10401
      t10409 = t9905*t10401
      t10412 = t9761*t6961*t10348
      t10414 = t7479*t10344
      t10415 = t10400*t10414
      t10417 = t9905*t10414
      t10419 = t9799*t10414
      t10421 = t6953*t10344
      t10422 = t9799*t10421
      t10424 = t10400*t10421
      t10426 = t6943*t10325
      t10427 = t9731*t10426
      t10429 = t6953*t10363
      t10430 = t9718*t10429
      t10432 = t10066*t7428
      t10434 = t720*t10182
      t10435 = t6943*t10434
      t10436 = t9852*t10435
      t10438 = t3692*t10182
      t10439 = t6961*t10438
      t10440 = t9814*t10439
      t10442 = t9799*t10018
      t10444 = -8.D0/425.D0*t10402+8.D0/425.D0*t10405+12.D0/425.D0*t1040
     #7-4.D0/425.D0*t10409+2.D0/425.D0*t10412+16.D0/425.D0*t10415+8.D0/4
     #25.D0*t10417-24.D0/425.D0*t10419+24.D0/425.D0*t10422-16.D0/425.D0*
     #t10424+t10427/255-t10430/100+4.D0/255.D0*t10432-7.D0/850.D0*t10436
     #-2.D0/425.D0*t10440-43.D0/3825.D0*t10442
      t10449 = t9799*t10048
      t10451 = t3970*spf14
      t10452 = t6943*t10451
      t10453 = t9852*t10452
      t10455 = t6949*t9734
      t10456 = t9663*t10455
      t10458 = t6961*t10451
      t10459 = t10021*t10458
      t10461 = t9761*t10458
      t10463 = t10233*t7723
      t10465 = t10066*t7076
      t10467 = t10066*t6939
      t10469 = t1093*wC3
      t10470 = t6967*t10469
      t10471 = t9718*t10470
      t10476 = t10180*t3751*t10217
      t10478 = t10180*t10166
      t10481 = t10180*t3893*t10217
      t10484 = t10180*t3842*t10217
      t10486 = t10180*t10206
      t10488 = 43.D0/3825.D0*t10449+13.D0/2550.D0*t10453+t10456/3060-13.
     #D0/2550.D0*t10459-t10461/900+t10463/7650-7.D0/2550.D0*t10465-7.D0/
     #5100.D0*t10467-t10471/1700+3.D0/6800.D0*t9678*t10470-t10476/2550+t
     #10478/850+13.D0/7650.D0*t10481+t10484/6120-t10486/1700
      t10489 = t6987*t9782
      t10490 = t10113*t10489
      t10492 = t4892*t3420
      t10494 = t10492*t3842*t6987
      t10496 = t7408*t9782
      t10497 = t10113*t10496
      t10499 = t6953*t9654
      t10502 = t6943*t9817
      t10505 = t6943*t9810
      t10506 = t9714*t10505
      t10508 = t4468*t3420
      t10509 = t10508*t10226
      t10511 = t9767*t9966
      t10513 = t9714*t7431
      t10515 = t10249*t10505
      t10517 = t7588*t10053
      t10518 = t9833*t10517
      t10520 = t10113*t10298
      t10522 = t6943*t9834
      t10523 = t10508*t10522
      t10525 = t7648*t10043
      t10526 = t6967*t10525
      t10527 = t9761*t10526
      t10529 = t6943*t10121
      t10530 = t10214*t10529
      t10532 = t9833*t9870
      t10534 = -t10490/425+t10494/12240-t10497/425-3.D0/1700.D0*t9678*t1
     #0499-11.D0/3400.D0*t10492*t10502+t10506/1530+73.D0/10200.D0*t10509
     #-2.D0/765.D0*t10511-t10513/5100-13.D0/3825.D0*t10515-t10518/10200+
     #3.D0/425.D0*t10520+t10523/425-t10527/10200-t10530/40800+13.D0/3060
     #0.D0*t10532
      t10536 = t10066*t10505
      t10538 = t10113*t9811
      t10540 = t6967*t9810
      t10541 = t10113*t10540
      t10543 = t10272*t7515
      t10545 = t6987*t9834
      t10550 = t9666*t9754
      t10552 = t842*spf14
      t10553 = t6943*t10552
      t10554 = t10180*t10553
      t10556 = CsTh*spf14
      t10557 = t6943*t10556
      t10558 = t10180*t10557
      t10560 = t10113*t10175
      t10562 = t10113*t9934
      t10564 = t10180*t10229
      t10566 = t6943*t9933
      t10567 = t10249*t10566
      t10569 = t10113*t9740
      t10571 = t9714*t10566
      t10573 = t9898*t10455
      t10575 = 8.D0/1275.D0*t10536-3.D0/850.D0*t10538-3.D0/1700.D0*t1054
     #1-3.D0/1700.D0*t10543+11.D0/850.D0*t10113*t10545+t9814*t10526/1360
     #0-10.D0/153.D0*t10550-t10554/3600+t10558/1800+3.D0/850.D0*t10560-3
     #.D0/425.D0*t10562-t10564/3400+43.D0/15300.D0*t10567+239.D0/7650.D0
     #*t10569-t10571/1530-5.D0/153.D0*t10573
      t10576 = t10508*t10036
      t10578 = t6943*t9757
      t10579 = t10492*t10578
      t10581 = t10066*t10566
      t10583 = t6937*t10182
      t10584 = t3429*t10583
      t10585 = t10180*t10584
      t10587 = t10180*t10336
      t10589 = t6987*t10348
      t10590 = t10113*t10589
      t10593 = t10180*t3751*t9957
      t10596 = t10180*t3893*t9957
      t10598 = t6967*t9933
      t10599 = t10113*t10598
      t10602 = t10180*t3842*t9957
      t10604 = t10180*t9958
      t10606 = t6943*t9719
      t10607 = t10508*t10606
      t10610 = t10180*t3751*t10583
      t10612 = CsTh*t10182
      t10613 = t6943*t10612
      t10614 = t10180*t10613
      t10616 = t6943*t10371
      t10617 = t10249*t10616
      t10619 = t9666*t10393
      t10621 = t10576/900-239.D0/30600.D0*t10579-29.D0/5100.D0*t10581+t1
     #0585/1700+t10587/850-t10590/50+t10593/612-13.D0/10200.D0*t10596+3.
     #D0/1700.D0*t10599-59.D0/61200.D0*t10602-t10604/7200-271.D0/20400.D
     #0*t10607+t10610/1700-t10614/425+3.D0/425.D0*t10617-6.D0/425.D0*t10
     #619
      t10625 = t9647*t1812*t10616
      t10627 = t10508*t10380
      t10629 = t10113*t7019
      t10631 = t10113*t7612
      t10633 = t9857*t10393
      t10635 = t10113*t7595
      t10637 = t10508*t10356
      t10639 = t10492*t10364
      t10641 = t10066*t10616
      t10643 = t10066*t7829
      t10645 = t10113*t7799
      t10647 = t9738*t10496
      t10649 = t10113*t7900
      t10651 = t10113*t7547
      t10653 = t10113*t7903
      t10655 = -2.D0/425.D0*t10625-29.D0/5100.D0*t10627-26.D0/225.D0*t10
     #629-13.D0/900.D0*t10631+2.D0/425.D0*t10633-t10635/7650-2.D0/425.D0
     #*t10637+t10639/200-t10641/425+13.D0/7650.D0*t10643+253.D0/2550.D0*
     #t10645+t10647/425+13.D0/450.D0*t10649+t10651/3825+253.D0/5100.D0*t
     #10653
      t10657 = t9772*t6943*t10438
      t10659 = t6953*t10379
      t10660 = t9761*t10659
      t10662 = t9714*t7829
      t10664 = t3999*t10182
      t10665 = t6953*t10664
      t10666 = t9857*t10665
      t10668 = t6953*t10348
      t10669 = t9761*t10668
      t10671 = t5598*t3426
      t10672 = t3606*t7498
      t10673 = t10671*t10672
      t10675 = t9738*t10665
      t10677 = t6949*t10371
      t10678 = t9666*t10677
      t10680 = t9772*t10435
      t10682 = t9857*t10677
      t10684 = t9714*t7428
      t10686 = t6949*t10121
      t10687 = t9761*t10686
      t10689 = t6943*t10525
      t10692 = t9714*t7100
      t10694 = t9714*t7723
      t10696 = t7479*t10375
      t10697 = t9747*t10696
      t10699 = 2.D0/425.D0*t10657+29.D0/2550.D0*t10660-53.D0/7650.D0*t10
     #662-8.D0/255.D0*t10666+4.D0/425.D0*t10669-3.D0/3400.D0*t10673-t106
     #75/25+3.D0/425.D0*t10678+13.D0/10200.D0*t10680-t10682/425+329.D0/3
     #0600.D0*t10684-t10687/5100+t9772*t10689/6800+t10692/1275-53.D0/382
     #5.D0*t10694-14.D0/425.D0*t10697
      t10701 = t9747*t10376
      t10703 = t10065*t3426
      t10704 = t10703*t10672
      t10706 = t9924*t10376
      t10708 = t9718*t10376
      t10710 = t9666*t10372
      t10712 = t5590*t3426
      t10713 = t10712*t10672
      t10715 = t9678*t10696
      t10717 = t9861*t10665
      t10719 = t9857*t10455
      t10721 = t9691*t9740
      t10723 = t5992*spf14
      t10724 = t6967*t10723
      t10725 = t9727*t10724
      t10727 = t1096*spf14
      t10728 = t6967*t10727
      t10729 = t10021*t10728
      t10731 = t9898*t10665
      t10733 = t9666*t10665
      t10735 = t10339*t10372
      t10737 = t10278*t7829
      t10739 = 14.D0/425.D0*t10701-13.D0/30600.D0*t10704+58.D0/1275.D0*t
     #10706-58.D0/1275.D0*t10708+3.D0/850.D0*t10710+53.D0/30600.D0*t1071
     #3+13.D0/2550.D0*t10715+8.D0/255.D0*t10717-163.D0/2550.D0*t10719+8.
     #D0/3825.D0*t10721-5.D0/1224.D0*t10725+91.D0/15300.D0*t10729+4.D0/2
     #55.D0*t10731-4.D0/255.D0*t10733-t10735/425-t10737/30600
      t10741 = t6953*t53*t9692
      t10742 = t9663*t10741
      t10744 = t10278*t7431
      t10746 = t10278*t7428
      t10748 = t6943*t10469
      t10751 = t6953*t9789
      t10752 = t9718*t10751
      t10754 = t9761*t10728
      t10756 = t9747*t10724
      t10758 = t10066*t7213
      t10760 = t10249*t7213
      t10762 = t10272*t7829
      t10764 = t10125*t9684
      t10766 = t10272*t7100
      t10768 = t9714*t7076
      t10770 = t9714*t7213
      t10772 = t10249*t7076
      t10774 = t10492*t9659
      t10776 = t10742/850-13.D0/1800.D0*t10744-253.D0/20400.D0*t10746+3.
     #D0/3400.D0*t9653*t10748-t10752/850-271.D0/40800.D0*t10754+163.D0/2
     #0400.D0*t10756+7.D0/1275.D0*t10758+t10760/3825+3.D0/850.D0*t10762-
     #6.D0/425.D0*t10764-t10766/3060-t10768/1530+t10770/765-t10772/7650+
     #t10774/1700
      t10780 = t10492*t9790
      t10782 = t10249*t6939
      t10785 = t9718*t6967*t10363
      t10788 = t3429*t7137*t10182
      t10789 = t9852*t10788
      t10791 = t6967*t10379
      t10792 = t9761*t10791
      t10794 = t9714*t6939
      t10796 = t9718*t10131
      t10798 = t10113*t10200
      t10800 = t10113*t9667
      t10802 = t9718*t10724
      t10804 = t6943*t10723
      t10805 = t9673*t10804
      t10807 = t9861*t10455
      t10809 = t9857*t9740
      t10815 = t10780/1700-t10782/15300+t10785/400+7.D0/3400.D0*t10789-2
     #9.D0/10200.D0*t10792-t10794/3060+t10796/1275+3.D0/1700.D0*t10798-3
     #.D0/850.D0*t10800-239.D0/61200.D0*t10802+t10805/12240+163.D0/2550.
     #D0*t10807-163.D0/2550.D0*t10809+3.D0/1700.D0*t9691*t10496-3.D0/170
     #0.D0*t9691*t9670
      t10816 = t9761*t10439
      t10818 = t9772*t10788
      t10820 = t9663*t9693
      t10822 = t9642*t10426
      t10824 = t1096*t10182
      t10825 = t6943*t10824
      t10826 = t9772*t10825
      t10828 = t6949*t10379
      t10829 = t9761*t10828
      t10831 = t9905*t10421
      t10833 = t9691*t10455
      t10835 = t6949*t9719
      t10836 = t9761*t10835
      t10838 = t9924*t10131
      t10844 = t10339*t10677
      t10846 = t3970*t10182
      t10847 = t6961*t10846
      t10848 = t9814*t10847
      t10851 = t9772*t6943*t10846
      t10853 = t9852*t10825
      t10855 = -4.D0/425.D0*t10816-13.D0/40800.D0*t10818-t10820/1700+2.D
     #0/255.D0*t10822-13.D0/20400.D0*t10826-29.D0/5100.D0*t10829-8.D0/42
     #5.D0*t10831+8.D0/3825.D0*t10833-271.D0/20400.D0*t10836-t10838/1275
     #+t9678*t10517/1700-3.D0/13600.D0*t9653*t9700-2.D0/425.D0*t10844+t1
     #0848/425-t10851/425+7.D0/1700.D0*t10853
      t10857 = t6949*t9723
      t10858 = t9678*t10857
      t10860 = t3970*t9757
      t10861 = t9781*t10860
      t10863 = t9683*t7071
      t10865 = t9814*t10728
      t10867 = t9951*t10553
      t10869 = t6949*t10375
      t10870 = t9678*t10869
      t10872 = t6949*t10363
      t10873 = t9718*t10872
      t10875 = t5992*t10348
      t10876 = t9781*t10875
      t10878 = t6949*t10664
      t10879 = t9738*t10878
      t10883 = t7588*t10375
      t10884 = t9678*t10883
      t10886 = t9861*t10878
      t10888 = t9857*t10878
      t10890 = t9738*t10455
      t10892 = t9691*t10175
      t10894 = t9718*t10857
      t10896 = -127.D0/5100.D0*t10858+239.D0/7650.D0*t10861+14.D0/1275.D
     #0*t10863+127.D0/40800.D0*t10865+t10867/612+13.D0/5100.D0*t10870+t1
     #0873/200-t10876/50+t10879/50-t9678*t10131/1700-13.D0/5100.D0*t1088
     #4-4.D0/255.D0*t10886+4.D0/255.D0*t10888-239.D0/7650.D0*t10890-t108
     #92/1530+271.D0/5100.D0*t10894
      t10897 = t9718*t10869
      t10899 = t6967*t9789
      t10900 = t9718*t10899
      t10904 = t6943*t10375
      t10905 = t9653*t10904
      t10909 = t9951*t3429*t6953*t10182
      t10911 = t9924*t10869
      t10913 = t9747*t10883
      t10915 = t9747*t10869
      t10917 = t6943*t10359
      t10918 = t9731*t10917
      t10920 = t6949*t10348
      t10921 = t9761*t10920
      t10923 = t9898*t10878
      t10925 = t9666*t10878
      t10927 = t9778*t10455
      t10929 = t10278*t7515
      t10931 = t6949*t10344
      t10932 = t10400*t10931
      t10934 = t5598*t9743
      t10935 = t10934*t7829
      t10937 = 29.D0/1275.D0*t10897+t10900/3400-3.D0/13600.D0*t9653*t970
     #9+13.D0/5100.D0*t10905-t10909/1700-29.D0/1275.D0*t10911+7.D0/425.D
     #0*t10913-7.D0/425.D0*t10915-t10918/510-2.D0/425.D0*t10921-2.D0/255
     #.D0*t10923+2.D0/255.D0*t10925-8.D0/3825.D0*t10927+t10929/61200+8.D
     #0/425.D0*t10932+3.D0/1700.D0*t10935
      t10940 = t10261*t7900
      t10943 = t10144*t6967*t10612
      t10945 = t9727*t10857
      t10947 = t6938*spf14
      t10948 = t10113*t10947
      t10952 = t9951*t3429*t6949*t10182
      t10955 = t10144*t6967*t10335
      t10957 = t7588*t10344
      t10958 = t10400*t10957
      t10960 = t7588*t9723
      t10961 = t9727*t10960
      t10964 = t10180*t6943*t9794
      t10966 = t9678*t10724
      t10968 = t6967*t9658
      t10969 = t9718*t10968
      t10971 = t10214*t10689
      t10973 = t6254*spf14
      t10974 = t6967*t10973
      t10975 = t9691*t10974
      t10977 = t10214*t10452
      t10979 = t9738*t10489
      t10981 = 13.D0/900.D0*t10940-t10943/850+26.D0/1275.D0*t10945+3.D0/
     #850.D0*t10948+t10952/3400+t10955/1700-8.D0/425.D0*t10958-26.D0/127
     #5.D0*t10961-t10964/1700-t10966/3825+t10969/3400-t10971/40800-t1097
     #5/3060+59.D0/15300.D0*t10977+t10979/425
      t10984 = t9905*t10957
      t10986 = t9799*t10957
      t10988 = t9799*t10931
      t10990 = t6967*t9839
      t10991 = t9761*t10990
      t10993 = t8093*t147
      t10994 = t6949*t10993
      t10995 = t9718*t10994
      t10997 = t9924*t10994
      t10999 = t9905*t10931
      t11001 = t9951*t10584
      t11003 = t9642*t10917
      t11005 = t9761*t10847
      t11007 = t7588*t10993
      t11008 = t9678*t11007
      t11010 = t9747*t10857
      t11012 = t6938*t9794
      t11013 = t10113*t11012
      t11015 = t7994*t10127
      t11016 = t9703*t11015
      t11018 = t445*t147
      t11019 = t6949*t11018
      t11020 = t9718*t11019
      t11022 = 3.D0/1700.D0*t9691*t10489-4.D0/425.D0*t10984+12.D0/425.D0
     #*t10986-12.D0/425.D0*t10988+73.D0/20400.D0*t10991-8.D0/3825.D0*t10
     #995+8.D0/3825.D0*t10997+4.D0/425.D0*t10999+t11001/3400-t11003/255+
     #2.D0/425.D0*t11005-11.D0/7650.D0*t11008-182.D0/3825.D0*t11010-3.D0
     #/850.D0*t11013-64.D0/3825.D0*t11016+t11020/510
      t11024 = t3609*t147
      t11025 = t1865*t11024
      t11026 = t9781*t11025
      t11028 = t9861*t9754
      t11030 = t6943*t10727
      t11031 = t10214*t11030
      t11033 = t9898*t9754
      t11035 = t7994*t11024
      t11038 = t6943*t11018
      t11041 = t6943*t10993
      t11042 = t9653*t11041
      t11046 = t9727*t10360
      t11048 = t10934*t7515
      t11050 = t10261*t7903
      t11052 = t8169*t11024
      t11055 = t1093*t9834
      t11058 = t1224*t147
      t11059 = t6967*t11058
      t11060 = t10169*t11059
      t11062 = t1865*t147
      t11063 = t6967*t11062
      t11064 = t10169*t11063
      t11066 = t9747*t10960
      t11068 = -2.D0/255.D0*t11026-163.D0/1275.D0*t11028+t11031/24480+10
     #.D0/153.D0*t11033+2.D0/255.D0*t9691*t11035-t9653*t11038/510+11.D0/
     #7650.D0*t11042+2.D0/255.D0*t10106*t11025-t11046/1020-3.D0/3400.D0*
     #t11048+253.D0/10200.D0*t11050-2.D0/255.D0*t9778*t11052-11.D0/425.D
     #0*t9781*t11055-t11060/3600+t11064/7200+182.D0/3825.D0*t11066
      t11069 = t10027*t11025
      t11071 = t9691*t10298
      t11074 = t9647*t3780*t3426
      t11075 = t1852*t147
      t11076 = t6949*t11075
      t11077 = t11074*t11076
      t11079 = t6953*t9839
      t11080 = t9761*t11079
      t11082 = t9857*t9754
      t11084 = t9691*t9754
      t11086 = t7299*t9654
      t11087 = t9833*t11086
      t11089 = t10113*t7625
      t11091 = t9747*t10295
      t11093 = t10021*t10223
      t11095 = t10065*t9887
      t11096 = t11095*t6944
      t11099 = t9647*t3893*t3426
      t11100 = t11099*t11076
      t11102 = t6953*t9817
      t11107 = t9678*t10186
      t11109 = t9678*t9973
      t11111 = 2.D0/255.D0*t11069-t11071/765+8.D0/3825.D0*t11077-73.D0/5
     #100.D0*t11080+163.D0/1275.D0*t11082-16.D0/3825.D0*t11084+t11087/68
     #00-253.D0/10200.D0*t11089-83.D0/15300.D0*t11091+41.D0/3060.D0*t110
     #93-182.D0/765.D0*t11096-16.D0/3825.D0*t11100+11.D0/1700.D0*t9718*t
     #11102-11.D0/6800.D0*t9718*t10186+13.D0/30600.D0*t11107+127.D0/2550
     #.D0*t11109
      t11119 = t9678*t10295
      t11121 = t9727*t9836
      t11123 = t10021*t10170
      t11125 = t62*t147
      t11126 = t1865*t11125
      t11127 = t9781*t11126
      t11129 = t9951*t10206
      t11131 = t10261*t11076
      t11133 = t7479*t9835
      t11134 = t9727*t11133
      t11136 = t9778*t9754
      t11138 = t9727*t9973
      t11140 = t3692*t9757
      t11141 = t9781*t11140
      t11143 = t9738*t9754
      t11145 = t9718*t9973
      t11147 = t10066*t7572
      t11149 = t6967*t10824
      t11150 = t10021*t11149
      t11152 = 11.D0/3400.D0*t9718*t10295-13.D0/15300.D0*t11119+208.D0/3
     #825.D0*t11121-41.D0/6120.D0*t11123+16.D0/3825.D0*t11127-t11129/255
     #0+8.D0/3825.D0*t11131-208.D0/3825.D0*t11134+16.D0/3825.D0*t11136-5
     #2.D0/1275.D0*t11138-239.D0/3825.D0*t11141+239.D0/3825.D0*t11143-27
     #1.D0/2550.D0*t11145-13.D0/15300.D0*t11147+7.D0/3400.D0*t11150
      t11153 = t9761*t11149
      t11155 = t9857*t10540
      t11157 = t9924*t10857
      t11159 = t10113*t7621
      t11161 = t9814*t10170
      t11163 = t6254*t9794
      t11164 = t6967*t11163
      t11165 = t9691*t11164
      t11167 = t10113*t7942
      t11169 = t9673*t9846
      t11171 = t10214*t10007
      t11173 = t10113*t7416
      t11175 = t9678*t10960
      t11177 = t103*t147
      t11178 = t6943*t11177
      t11179 = t10508*t11178
      t11181 = t9857*t7900
      t11183 = t9666*t7903
      t11185 = t10249*t7431
      t11187 = t9814*t11059
      t11189 = -29.D0/10200.D0*t11153+4.D0/1275.D0*t11155-271.D0/5100.D0
     #*t11157-13.D0/900.D0*t11159-59.D0/122400.D0*t11161+t11165/3060-253
     #.D0/10200.D0*t11167-t11169/12240-t11171/1530-t11173/7650+127.D0/51
     #00.D0*t11175+2.D0/3825.D0*t11179-91.D0/1530.D0*t11181+26.D0/3825.D
     #0*t11183-133.D0/7650.D0*t11185+11.D0/30600.D0*t11187
      t11191 = t6943*t11058
      t11192 = t10214*t11191
      t11194 = t9666*t11164
      t11196 = t9857*t11164
      t11198 = t9814*t11149
      t11200 = t9718*t10360
      t11202 = t9814*t11063
      t11204 = t6943*t11062
      t11205 = t10214*t11204
      t11207 = t9761*t11059
      t11209 = t9696*t7224
      t11211 = t9696*t7071
      t11213 = t9852*t10291
      t11215 = t9666*t7900
      t11218 = t3429*t7994*t147
      t11219 = t10214*t11218
      t11221 = t9849*t10291
      t11223 = t7479*t9723
      t11224 = t9727*t11223
      t11226 = t7103*t9794
      t11227 = t3623*t11226
      t11228 = t10214*t11227
      t11230 = -t11192/1800-13.D0/7650.D0*t11194+4.D0/1275.D0*t11196-13.
     #D0/40800.D0*t11198+t11200/400-11.D0/61200.D0*t11202+t11205/3600-2.
     #D0/3825.D0*t11207-4.D0/765.D0*t11209+2.D0/765.D0*t11211+41.D0/1530
     #.D0*t11213+133.D0/3825.D0*t11215+t11219/7200-52.D0/3825.D0*t11221+
     #52.D0/1275.D0*t11224+t11228/1530
      t11231 = t6943*t9806
      t11232 = t9696*t11231
      t11234 = t9666*t7019
      t11236 = t9544*t9887
      t11237 = t11236*t6944
      t11239 = t10189*t10206
      t11241 = t9647*t9715
      t11243 = t9731*t10253
      t11247 = t3409*t147
      t11248 = t445*t11247
      t11249 = t9781*t11248
      t11251 = t9647*t7572
      t11253 = t9747*t9973
      t11255 = t9673*t10116
      t11257 = t9738*t11035
      t11260 = t9761*t6961*t9739
      t11262 = t9772*t10522
      t11264 = t9772*t10219
      t11266 = t6943*t9835
      t11267 = t9673*t11266
      t11269 = -2.D0/765.D0*t11232-532.D0/3825.D0*t11234+532.D0/3825.D0*
     #t11237+13.D0/7650.D0*t11239+3.D0/850.D0*t11241+26.D0/3825.D0*t1124
     #3-2.D0/255.D0*t9778*t11035-2.D0/255.D0*t11249+t11251/30600+364.D0/
     #3825.D0*t11253+t11255/10200+2.D0/255.D0*t11257-t11260/900+2.D0/127
     #5.D0*t11262-59.D0/122400.D0*t11264+13.D0/30600.D0*t11267
      t11272 = t9642*t10253
      t11274 = t6953*t9821
      t11275 = t9648*t11274
      t11277 = t9738*t10303
      t11279 = t9857*t7019
      t11281 = t9691*t7059
      t11283 = t10249*t7428
      t11287 = t9648*t10303
      t11289 = t9691*t7595
      t11291 = t9747*t11223
      t11293 = t6953*t11018
      t11294 = t9718*t11293
      t11296 = t9772*t10320
      t11298 = t7994*t11125
      t11299 = t10113*t11298
      t11301 = t9666*t10455
      t11303 = t9833*t10857
      t11305 = -83.D0/7650.D0*t11272-t11275/1530-2.D0/425.D0*t11277+182.
     #D0/765.D0*t11279+t11281/1275-13.D0/3825.D0*t11283-3.D0/1700.D0*t96
     #91*t9688+t11287/1700-53.D0/3825.D0*t11289-364.D0/3825.D0*t11291-t1
     #1294/255+t11296/153+16.D0/3825.D0*t11299+5.D0/153.D0*t11301-t11303
     #/6120
      t11306 = t9731*t10578
      t11308 = t9833*t9728
      t11310 = t9924*t9973
      t11312 = t6953*t10993
      t11313 = t9718*t11312
      t11315 = t9663*t11274
      t11317 = t9653*t10253
      t11319 = t9663*t10303
      t11321 = t9691*t10540
      t11323 = t9772*t10291
      t11325 = t9772*t11030
      t11327 = t9691*t7625
      t11329 = t9849*t10606
      t11331 = t9738*t7595
      t11333 = t9691*t7612
      t11335 = t9888*t7601
      t11337 = t9648*t9740
      t11339 = -5.D0/612.D0*t11306+t11308/6120+271.D0/2550.D0*t11310+16.
     #D0/3825.D0*t11313+t11315/1530-13.D0/7650.D0*t11317-t11319/1700+t11
     #321/3060+59.D0/30600.D0*t11323+127.D0/20400.D0*t11325+329.D0/30600
     #.D0*t11327-13.D0/2550.D0*t11329+t11331/7650-t11333/5100+53.D0/3825
     #.D0*t11335-t11337/3060
      t11341 = t9852*t10606
      t11343 = t9642*t10578
      t11345 = t9956*t10150
      t11347 = t9956*t10191
      t11349 = t10113*t7111
      t11351 = t9666*t9735
      t11353 = t9857*t9735
      t11355 = t8308*t10127
      t11356 = t10161*t11355
      t11358 = t10249*t7100
      t11360 = t9905*t11355
      t11362 = t10272*t7723
      t11364 = t9678*t11223
      t11366 = t9666*t10540
      t11368 = t10113*t7674
      t11370 = t9648*t10545
      t11372 = t9673*t10502
      t11374 = 91.D0/7650.D0*t11341+163.D0/10200.D0*t11343-t11345/6120+t
     #11347/12240+13.D0/225.D0*t11349-5.D0/153.D0*t11351+163.D0/2550.D0*
     #t11353+256.D0/3825.D0*t11356+266.D0/3825.D0*t11358-256.D0/3825.D0*
     #t11360+3.D0/425.D0*t11362-127.D0/2550.D0*t11364-13.D0/7650.D0*t113
     #66+t11368/30600+t11370/3060-t11372/12240
      t11375 = t9857*t11274
      t11377 = t9861*t11274
      t11379 = t9653*t10578
      t11381 = t9696*t9976
      t11383 = t9857*t7625
      t11385 = t9857*t7595
      t11387 = t9653*t9724
      t11389 = t9772*t10036
      t11391 = t10339*t10368
      t11393 = t9772*t10226
      t11395 = t10214*t10522
      t11397 = t6961*t9821
      t11398 = t9648*t11397
      t11400 = t9653*t10502
      t11402 = t9673*t9724
      t11404 = t9772*t10606
      t11406 = t9673*t7184
      t11408 = 166.D0/3825.D0*t11375-166.D0/3825.D0*t11377-2.D0/3825.D0*
     #t11379+2.D0/765.D0*t11381+4.D0/255.D0*t11383+13.D0/3825.D0*t11385-
     #127.D0/5100.D0*t11387-t11389/153-t11391/425-59.D0/61200.D0*t11393-
     #t11395/1530-t11398/3060+13.D0/15300.D0*t11400-t11402/6120+127.D0/2
     #0400.D0*t11404-t11406/12240
      t11414 = t9673*t9655
      t11416 = t11095*t7601
      t11418 = t9803*t9684
      t11420 = t7137*spf14
      t11421 = t3429*t11420
      t11422 = t10214*t11421
      t11424 = t6987*spf14
      t11425 = t3429*t11424
      t11426 = t9673*t11425
      t11428 = t9888*t7023
      t11430 = t9738*t7059
      t11432 = t10189*t9958
      t11434 = t9913*t11355
      t11438 = t9642*t10502
      t11440 = t9803*t11231
      t11442 = t9731*t10502
      t11444 = t9849*t10226
      t11446 = 3.D0/6800.D0*t9678*t9786-t11414/6800-13.D0/3825.D0*t11416
     #+2.D0/3825.D0*t11418+t11422/48960+t11426/24480-t11428/1275-13.D0/2
     #25.D0*t11430-13.D0/20400.D0*t11432-128.D0/3825.D0*t11434+t9772*t10
     #529/6800+83.D0/15300.D0*t11438+52.D0/3825.D0*t11440-13.D0/3825.D0*
     #t11442+26.D0/3825.D0*t11444
      t11447 = t9852*t10226
      t11449 = t9951*t9958
      t11451 = t6961*t9806
      t11452 = t9799*t11451
      t11454 = t9683*t11231
      t11456 = t9666*t10368
      t11458 = t9761*t11063
      t11460 = t9666*t11298
      t11462 = t10214*t10320
      t11464 = t6961*t10664
      t11465 = t9857*t11464
      t11467 = t9744*t7893
      t11469 = t9691*t11274
      t11471 = t9849*t11421
      t11473 = t9857*t10598
      t11475 = t9738*t10589
      t11477 = t9857*t10589
      t11481 = -41.D0/3060.D0*t11447+t11449/1224-52.D0/3825.D0*t11452-32
     #.D0/1275.D0*t11454+3.D0/850.D0*t11456+t11458/3825+16.D0/3825.D0*t1
     #1460-59.D0/15300.D0*t11462-4.D0/255.D0*t11465-329.D0/30600.D0*t114
     #67+26.D0/3825.D0*t11469-13.D0/5100.D0*t11471-29.D0/10200.D0*t11473
     #+t11475/50+4.D0/255.D0*t11477+11.D0/425.D0*t9738*t11274
      t11483 = t5598*t9887
      t11484 = t11483*t7723
      t11487 = t6961*t103*t10127
      t11488 = t9905*t11487
      t11490 = t9905*t11015
      t11492 = t9799*t11487
      t11494 = t9799*t11015
      t11496 = t9857*t11298
      t11498 = t9778*t11274
      t11500 = t9666*t10545
      t11502 = t9666*t11274
      t11504 = t9898*t11274
      t11506 = t6961*t10993
      t11507 = t9678*t11506
      t11509 = t9913*t11015
      t11511 = t6961*t10375
      t11512 = t9678*t11511
      t11514 = t7103*t10379
      t11515 = t9718*t11514
      t11517 = t7479*t10993
      t11518 = t9678*t11517
      t11520 = t9673*t10748
      t11522 = -3.D0/425.D0*t11484-128.D0/3825.D0*t11488+128.D0/3825.D0*
     #t11490+64.D0/3825.D0*t11492-64.D0/3825.D0*t11494-32.D0/3825.D0*t11
     #496-26.D0/3825.D0*t11498+52.D0/3825.D0*t11500-104.D0/3825.D0*t1150
     #2+104.D0/3825.D0*t11504-11.D0/7650.D0*t11507+64.D0/3825.D0*t11509-
     #13.D0/5100.D0*t11512+29.D0/1275.D0*t11515+11.D0/3825.D0*t11518-t11
     #520/6800
      t11523 = t9648*t10489
      t11527 = t7994*t11247
      t11528 = t9718*t11527
      t11530 = t10070*t10028
      t11532 = t9718*t10499
      t11534 = t9666*t10598
      t11536 = t9731*t11425
      t11538 = t9718*t11506
      t11541 = t7994*Mzsq*t147
      t11542 = t9718*t11541
      t11544 = t9833*t10094
      t11546 = t4892*t7543
      t11547 = t11546*t11293
      t11549 = t9905*t11451
      t11551 = t6938*t9834
      t11552 = t9703*t11551
      t11554 = t9852*t11421
      t11556 = t9718*t9836
      t11558 = t9718*t11511
      t11560 = -t11523/3400+3.D0/3400.D0*t9678*t10031-8.D0/3825.D0*t1152
     #8-t11530/850+t11532/425+43.D0/30600.D0*t11534-5.D0/1224.D0*t11536+
     #8.D0/3825.D0*t11538+t11542/510-t11544/10200-t11547/510+32.D0/1275.
     #D0*t11549+6.D0/425.D0*t11552+91.D0/15300.D0*t11554+73.D0/1275.D0*t
     #11556-29.D0/1275.D0*t11558
      t11563 = t9747*t11133
      t11565 = t9691*t10545
      t11567 = t9747*t9836
      t11569 = t9703*t11451
      t11571 = t9767*t11451
      t11573 = t6967*t9719
      t11574 = t9761*t11573
      t11576 = t9642*t11425
      t11581 = t9833*t11312
      t11583 = t6967*t9757
      t11584 = t9718*t11583
      t11586 = t6961*t9835
      t11587 = t9833*t11586
      t11589 = t9673*t10578
      t11591 = t9772*t11421
      t11595 = t9852*t10522
      t11597 = 82.D0/765.D0*t11563-13.D0/3825.D0*t11565-82.D0/765.D0*t11
     #567-6.D0/425.D0*t11569+2.D0/765.D0*t11571-271.D0/40800.D0*t11574+1
     #63.D0/20400.D0*t11576+t9678*t7299*t11018/510+t11581/450-239.D0/612
     #00.D0*t11584-13.D0/30600.D0*t11587+t11589/12240+127.D0/40800.D0*t1
     #1591+t9678*t11293/255-26.D0/3825.D0*t11595
      t11598 = t9852*t11227
      t11600 = t9642*t10904
      t11602 = t9747*t11511
      t11604 = t9799*t11355
      t11606 = t6953*t11075
      t11607 = t9857*t11606
      t11609 = t9666*t11606
      t11611 = t9678*t9836
      t11613 = t9924*t9836
      t11615 = t9857*t10545
      t11617 = t9653*t11425
      t11620 = t9761*t6961*t9834
      t11622 = t10126*t11355
      t11624 = t9772*t11227
      t11626 = t9703*t9807
      t11628 = t9678*t11133
      t11632 = 26.D0/3825.D0*t11598-7.D0/425.D0*t11600+7.D0/425.D0*t1160
     #2+128.D0/3825.D0*t11604+64.D0/3825.D0*t11607-32.D0/3825.D0*t11609-
     #59.D0/7650.D0*t11611-73.D0/1275.D0*t11613-83.D0/3825.D0*t11615-t11
     #617/3825-t11620/425-128.D0/3825.D0*t11622-2.D0/1275.D0*t11624-12.D
     #0/425.D0*t11626+59.D0/7650.D0*t11628-11.D0/850.D0*t9738*t10545
      t11634 = t9924*t11312
      t11636 = t9833*t10091
      t11638 = t6953*t10121
      t11639 = t9761*t11638
      t11641 = t9666*t10589
      t11643 = t9691*t11397
      t11645 = t9857*t11397
      t11649 = t6987*t9794
      t11650 = t3429*t11649
      t11651 = t9731*t11650
      t11653 = t9666*t11397
      t11655 = t9642*t11650
      t11657 = t6967*t9817
      t11660 = t9852*t10219
      t11662 = t6967*t11177
      t11663 = t9761*t11662
      t11665 = t9653*t11650
      t11667 = t10261*t11606
      t11669 = t9833*t10994
      t11671 = -16.D0/3825.D0*t11634+t11636/5100+t11639/2550+2.D0/255.D0
     #*t11641+13.D0/3825.D0*t11643+83.D0/3825.D0*t11645-3.D0/3400.D0*t96
     #78*t11086-13.D0/7650.D0*t11651-52.D0/3825.D0*t11653+83.D0/30600.D0
     #*t11655-11.D0/6800.D0*t9718*t11657-41.D0/6120.D0*t11660+t11663/382
     #5+13.D0/30600.D0*t11665-16.D0/3825.D0*t11667-t11669/900
      t11672 = t11099*t11606
      t11674 = t6949*t11177
      t11675 = t9761*t11674
      t11677 = t11074*t11606
      t11682 = t9772*t11204
      t11684 = t9833*t11007
      t11686 = t8308*t11024
      t11689 = t1224*t11024
      t11692 = t11546*t11019
      t11694 = t9666*t11464
      t11696 = t9951*t10613
      t11700 = t9833*t10054
      t11702 = t9833*t10295
      t11704 = t9857*t10332
      t11706 = t9747*t10326
      t11708 = 32.D0/3825.D0*t11672+2.D0/3825.D0*t11675-16.D0/3825.D0*t1
     #1677-t9678*t7295*t11018/1020-11.D0/30600.D0*t11682+t11684/900+4.D0
     #/255.D0*t9778*t11686-4.D0/255.D0*t10106*t11689+t11692/1020-2.D0/25
     #5.D0*t11694-t11696/850-t9772*t10316/3400-t11700/5100+t11702/12240+
     #t11704/425+t11706/255
      t11713 = t9727*t10326
      t11715 = t6967*t10434
      t11716 = t10021*t11715
      t11718 = t9905*t9800
      t11720 = t9857*t9667
      t11722 = t9761*t11715
      t11724 = t9814*t11715
      t11726 = t9648*t10081
      t11728 = t9673*t7079
      t11732 = t10169*t10223
      t11734 = t6967*t9643
      t11735 = t9833*t11734
      t11737 = t9956*t10557
      t11739 = t11546*t10031
      t11741 = t6949*t9789
      t11742 = t9718*t11741
      t11744 = t9683*t7170
      t11746 = t11713/510-7.D0/1700.D0*t11716-14.D0/1275.D0*t11718+7.D0/
     #2550.D0*t11720+29.D0/5100.D0*t11722+13.D0/20400.D0*t11724-t11726/3
     #400+t11728/12240-3.D0/850.D0*t9691*t10303+13.D0/122400.D0*t11732-t
     #11735/12240+59.D0/30600.D0*t11737-t11739/1700+t11742/1700-14.D0/12
     #75.D0*t11744
      t11747 = t9648*t9693
      t11749 = t10065*t9743
      t11750 = t11749*t7428
      t11752 = t9544*t9743
      t11753 = t11752*t7431
      t11755 = t9833*t10031
      t11757 = t11749*t7431
      t11759 = t7295*t9654
      t11762 = t7578*wC3
      t11763 = t6943*t11762
      t11766 = t6967*t10556
      t11767 = t10144*t11766
      t11770 = t3429*t6953*spf14
      t11771 = t10189*t11770
      t11773 = t6967*t9989
      t11774 = t10169*t11773
      t11776 = t10169*t10311
      t11778 = t9744*t7829
      t11780 = t9781*t5549
      t11782 = t10214*t10606
      t11784 = t9663*t10081
      t11786 = t6949*t9658
      t11787 = t9718*t11786
      t11789 = t11747/1700+2.D0/255.D0*t11750-133.D0/15300.D0*t11753-t11
     #755/6800+91.D0/6120.D0*t11757+3.D0/6800.D0*t9678*t11759-3.D0/6800.
     #D0*t9653*t11763+t11767/3600+13.D0/10200.D0*t11771-t11774/24480+t11
     #776/40800-53.D0/15300.D0*t11778-253.D0/1275.D0*t11780+t11782/24480
     #+t11784/3400+t11787/1700
      t11791 = t9691*t10598
      t11793 = t10214*t10036
      t11795 = t11752*t7893
      t11797 = t10169*t10526
      t11799 = t9781*t5552
      t11801 = t9857*t7628
      t11803 = t9857*t7547
      t11805 = t5992*wC3
      t11806 = t6943*t11805
      t11809 = t9691*t7628
      t11811 = t9951*t11770
      t11813 = t10189*t9953
      t11815 = t10169*t10728
      t11817 = t10189*t10150
      t11819 = t9673*t7373
      t11821 = t6967*t10552
      t11822 = t10144*t11821
      t11824 = t9673*t11763
      t11826 = -t11791/3060+59.D0/15300.D0*t11793+13.D0/3825.D0*t11795-t
     #11797/81600+253.D0/2550.D0*t11799-16.D0/255.D0*t11801-26.D0/3825.D
     #0*t11803-3.D0/6800.D0*t9653*t11806-329.D0/7650.D0*t11809-t11811/61
     #2-13.D0/20400.D0*t11813+t11815/48960-13.D0/7650.D0*t11817+t11819/2
     #4480-t11822/7200+t11824/13600
      t11827 = t9803*t7170
      t11829 = t9833*t10724
      t11831 = t9956*t10553
      t11833 = t9778*t7628
      t11835 = t9738*t7628
      t11837 = t9861*t7628
      t11839 = t9781*t11689
      t11842 = t7423*Mzsq*wC3
      t11843 = t9718*t11842
      t11845 = t10027*t11689
      t11847 = t9772*t11178
      t11849 = t9682*t7376
      t11850 = t1865*t10127
      t11851 = t11849*t11850
      t11853 = t1224*t11125
      t11854 = t9781*t11853
      t11856 = t9852*t9774
      t11858 = t9678*t10994
      t11860 = t10261*t11298
      t11862 = t9799*t10128
      t11864 = -2.D0/3825.D0*t11827+t11829/24480-59.D0/61200.D0*t11831+3
     #29.D0/7650.D0*t11833-253.D0/2550.D0*t11835+16.D0/255.D0*t11837+4.D
     #0/255.D0*t11839-t11843/850-4.D0/255.D0*t11845-11.D0/30600.D0*t1184
     #7-64.D0/3825.D0*t11851-32.D0/3825.D0*t11854-13.D0/1275.D0*t11856+1
     #1.D0/7650.D0*t11858+8.D0/3825.D0*t11860-64.D0/3825.D0*t11862
      t11867 = t9857*t11076
      t11871 = t9666*t11076
      t11873 = t6961*t9773
      t11874 = t10021*t11873
      t11876 = t9761*t11873
      t11878 = t9666*t7628
      t11880 = t9898*t7628
      t11882 = t9678*t11312
      t11884 = t9772*t11218
      t11886 = t1224*t10127
      t11887 = t11849*t11886
      t11889 = t10214*t11178
      t11891 = t9833*t11506
      t11893 = t6953*t11177
      t11894 = t9761*t11893
      t11896 = t6967*t11018
      t11897 = t9718*t11896
      t11899 = t9673*t11041
      t11901 = -32.D0/3825.D0*t11867-t9678*t11019/510+16.D0/3825.D0*t118
     #71+13.D0/1275.D0*t11874+t11876/450+52.D0/3825.D0*t11878-52.D0/3825
     #.D0*t11880-11.D0/3825.D0*t11882-11.D0/61200.D0*t11884+128.D0/3825.
     #D0*t11887+t11889/3600+t11891/900-4.D0/3825.D0*t11894+t11897/1020-t
     #11899/900
      t11904 = t6967*t11075
      t11905 = t9666*t11904
      t11907 = t11749*t7829
      t11909 = t9814*t11873
      t11911 = t9731*t9644
      t11913 = t9857*t11904
      t11915 = t9833*t11223
      t11917 = t9849*t9990
      t11919 = t9772*t11191
      t11921 = t9833*t11517
      t11923 = t9814*t10458
      t11925 = t9731*t10804
      t11927 = t9833*t10960
      t11929 = t9744*t7515
      t11931 = t11749*t7515
      t11933 = t6949*t9757
      t11934 = t9718*t11933
      t11936 = -t9678*t11896/1020+8.D0/3825.D0*t11905+13.D0/15300.D0*t11
     #907-2.D0/153.D0*t11909+5.D0/306.D0*t11911-16.D0/3825.D0*t11913-t11
     #915/3060+13.D0/1275.D0*t11917+11.D0/15300.D0*t11919-t11921/450+t11
     #923/153-5.D0/612.D0*t11925+t11927/6120+53.D0/30600.D0*t11929-13.D0
     #/30600.D0*t11931-239.D0/30600.D0*t11934
      t11938 = t9653*t10804
      t11940 = t9648*t10455
      t11942 = t9852*t11030
      t11944 = t9642*t10804
      t11946 = t9849*t11030
      t11948 = t9857*t7799
      t11950 = t6258*spf14
      t11951 = t6967*t11950
      t11952 = t9857*t11951
      t11954 = t10211*t11773
      t11956 = t10261*t9934
      t11958 = t9767*t10018
      t11960 = t9772*t10452
      t11962 = t6949*t9739
      t11963 = t9761*t11962
      t11965 = t9666*t11951
      t11967 = t10169*t11873
      t11969 = t10189*t10557
      t11971 = t9888*t7723
      t11973 = -2.D0/3825.D0*t11938-t11940/3060+91.D0/7650.D0*t11942+163
     #.D0/10200.D0*t11944-13.D0/2550.D0*t11946-16.D0/255.D0*t11948+29.D0
     #/5100.D0*t11952+13.D0/2550.D0*t11954-3.D0/850.D0*t11956+2.D0/765.D
     #0*t11958-t11960/153+t11963/900-43.D0/15300.D0*t11965+59.D0/7650.D0
     #*t11967+13.D0/5100.D0*t11969+53.D0/3825.D0*t11971
      t11974 = t11095*t7723
      t11976 = t9888*t7100
      t11978 = t9852*t9990
      t11980 = t6953*t9658
      t11981 = t9718*t11980
      t11983 = t10021*t11773
      t11985 = t11546*t10499
      t11987 = t7465*t9782
      t11988 = t9781*t11987
      t11990 = t9727*t11734
      t11992 = t9691*t7799
      t11994 = t11095*t7100
      t11996 = t11236*t7100
      t11998 = t9888*t7770
      t12000 = t9738*t7799
      t12002 = t9761*t11773
      t12004 = t9747*t11734
      t12006 = t11095*t7770
      t12008 = -13.D0/3825.D0*t11974-t11976/1275-91.D0/3825.D0*t11978-t1
     #1981/850-91.D0/7650.D0*t11983+t11985/850+2.D0/425.D0*t11988+5.D0/6
     #12.D0*t11990-329.D0/7650.D0*t11992+91.D0/765.D0*t11994-266.D0/3825
     #.D0*t11996+329.D0/7650.D0*t11998-253.D0/2550.D0*t12000+271.D0/2040
     #0.D0*t12002-163.D0/10200.D0*t12004+16.D0/255.D0*t12006
      t12012 = t9718*t11734
      t12014 = t9673*t9644
      t12016 = t9814*t11773
      t12018 = t9951*t10557
      t12020 = t9666*t7799
      t12022 = t11236*t7770
      t12024 = t9678*t11734
      t12026 = t9691*t11951
      t12028 = t10214*t9774
      t12030 = t11483*t7550
      t12032 = t11749*t7533
      t12034 = t9744*t7533
      t12036 = t9956*t11770
      t12038 = t10214*t9990
      t12040 = t11749*t7893
      t12042 = 239.D0/30600.D0*t12012-t12014/6120-127.D0/20400.D0*t12016
     #-t12018/306+52.D0/3825.D0*t12020-52.D0/3825.D0*t12022+2.D0/3825.D0
     #*t12024+t12026/1530-59.D0/7650.D0*t12028+6.D0/425.D0*t12030-91.D0/
     #3060.D0*t12032+t12034/5100+59.D0/61200.D0*t12036-t12038/12240-4.D0
     #/255.D0*t12040
      t12043 = t9905*t9807
      t12045 = t11752*t7533
      t12047 = t9648*t7547
      t12051 = t10261*t7625
      t12053 = t11752*t7638
      t12055 = t10261*t7612
      t12057 = t9956*t9953
      t12061 = t7479*t9806
      t12062 = t9905*t12061
      t12064 = t9714*t7515
      t12066 = t9648*t10496
      t12068 = t9857*t10298
      t12070 = t9799*t12061
      t12072 = t9744*t7642
      t12074 = t10214*t10135
      t12076 = 64.D0/1275.D0*t12043+133.D0/7650.D0*t12045-6.D0/425.D0*t1
     #2047-3.D0/850.D0*t9691*t9650-253.D0/20400.D0*t12051-13.D0/7650.D0*
     #t12053-13.D0/1800.D0*t12055-59.D0/122400.D0*t12057+3.D0/425.D0*t10
     #106*t10028-64.D0/1275.D0*t12062+53.D0/15300.D0*t12064-t12066/3400-
     #16.D0/1275.D0*t12068+104.D0/3825.D0*t12070-t12072/10200-t12074/816
     #00
      t12078 = t9799*t9807
      t12080 = t9666*t10298
      t12082 = t9673*t9659
      t12084 = t11483*t7601
      t12086 = t11749*t7642
      t12088 = t9772*t10007
      t12090 = t9673*t11650
      t12092 = t11749*t7638
      t12094 = t11752*t7642
      t12096 = t9648*t7595
      t12098 = t10211*t10728
      t12100 = t10189*t10229
      t12102 = t10214*t10219
      t12106 = t9857*t7059
      t12108 = t9781*t10028
      t12110 = -104.D0/3825.D0*t12078+26.D0/3825.D0*t12080+t12082/13600-
     #3.D0/425.D0*t12084+91.D0/6120.D0*t12086+2.D0/1275.D0*t12088-t12090
     #/24480+2.D0/255.D0*t12092-133.D0/15300.D0*t12094+3.D0/425.D0*t1209
     #6-13.D0/5100.D0*t12098+13.D0/15300.D0*t12100-13.D0/244800.D0*t1210
     #2-3.D0/425.D0*t9778*t10741-91.D0/765.D0*t12106-4.D0/425.D0*t12108
      t12111 = t9951*t10229
      t12113 = t9642*t11266
      t12115 = t11095*t7023
      t12117 = t9857*t7612
      t12119 = t7103*t9839
      t12120 = t9718*t12119
      t12122 = t9718*t11586
      t12124 = t9747*t11586
      t12126 = t9666*t7625
      t12128 = t9666*t7612
      t12130 = t9666*t7059
      t12132 = t9678*t11586
      t12134 = t9913*t9807
      t12136 = t9653*t11266
      t12138 = t11236*t7023
      t12140 = t9857*t7518
      t12142 = t9691*t7518
      t12144 = -t12111/5100+41.D0/765.D0*t12113+91.D0/765.D0*t12115+91.D
     #0/3060.D0*t12117-73.D0/2550.D0*t12120+73.D0/2550.D0*t12122-41.D0/7
     #65.D0*t12124-13.D0/3825.D0*t12126-133.D0/7650.D0*t12128+266.D0/382
     #5.D0*t12130-59.D0/15300.D0*t12132+12.D0/425.D0*t12134+59.D0/15300.
     #D0*t12136-266.D0/3825.D0*t12138+32.D0/255.D0*t12140+329.D0/3825.D0
     #*t12142
      t12147 = t9731*t11266
      t12149 = t9727*t11586
      t12151 = t10508*t10529
      t12153 = t9778*t7518
      t12155 = t9833*t11133
      t12157 = t9852*t9829
      t12159 = t9767*t12061
      t12161 = t10066*t7515
      t12163 = t10021*t9796
      t12165 = t9738*t7518
      t12167 = t9861*t7518
      t12169 = t6953*t9834
      t12170 = t9761*t12169
      t12172 = t9666*t7518
      t12174 = t9898*t7518
      t12176 = t11095*t7550
      t12178 = -104.D0/3825.D0*t12147+104.D0/3825.D0*t12149-t12151/5100-
     #329.D0/3825.D0*t12153+13.D0/15300.D0*t12155+52.D0/3825.D0*t12157-4
     #.D0/765.D0*t12159-13.D0/15300.D0*t12161-52.D0/3825.D0*t12163+253.D
     #0/1275.D0*t12165-32.D0/255.D0*t12167-2.D0/425.D0*t12170-104.D0/382
     #5.D0*t12172+104.D0/3825.D0*t12174+26.D0/3825.D0*t12176
      t12179 = t9738*t9650
      t12181 = t9888*t6944
      t12183 = t9738*t7019
      t12185 = t9744*t7428
      t12187 = t6961*t7216
      t12188 = t9905*t12187
      t12190 = t11752*t7428
      t12192 = t9799*t12187
      t12194 = t9803*t7224
      t12196 = t9833*t10470
      t12198 = t6967*t11762
      t12199 = t9833*t12198
      t12201 = t6967*t11805
      t12202 = t9833*t12201
      t12204 = t9833*t9679
      t12206 = t9683*t7224
      t12208 = t9673*t7248
      t12210 = t9673*t6989
      t12214 = -2.D0/425.D0*t12179+2.D0/1275.D0*t12181+26.D0/225.D0*t121
     #83+329.D0/61200.D0*t12185+28.D0/1275.D0*t12188-13.D0/7650.D0*t1219
     #0+4.D0/3825.D0*t12192-4.D0/3825.D0*t12194-t12196/13600+t12199/2720
     #0+t12202/27200-t12204/13600-28.D0/1275.D0*t12206-t12208/6120+t1221
     #0/24480-3.D0/13600.D0*t9678*t12198
      t12216 = t9833*t11759
      t12218 = t9673*t11806
      t12220 = t6938*t3970
      t12221 = t9703*t12220
      t12225 = t9833*t10499
      t12227 = t9744*t7431
      t12229 = t9691*t7547
      t12231 = t9648*t7059
      t12233 = t9718*t12198
      t12235 = t9648*t7194
      t12239 = t9718*t12201
      t12241 = t6938*t3692
      t12242 = t9703*t12241
      t12244 = t9767*t12187
      t12246 = t9648*t7031
      t12248 = t9648*t7001
      t12250 = -t12216/13600+t12218/13600-6.D0/425.D0*t12221+3.D0/3400.D
     #0*t9653*t9675+t12225/3400-t12227/10200+106.D0/3825.D0*t12229-t1223
     #1/3060+t12233/3400+t12235/3060-3.D0/13600.D0*t9678*t12201+t12239/3
     #400+12.D0/425.D0*t12242+4.D0/765.D0*t12244+t12246/3060-t12248/1530
      t12251 = t9648*t7019
      t12253 = t9703*t12187
      t12255 = t9767*t9800
      t12257 = t9703*t9800
      t12259 = t9691*t10200
      t12263 = t60*t9782
      t12264 = t9781*t12263
      t12266 = t10113*t7434
      t12269 = t10180*t6943*spf14
      t12271 = t9857*t10974
      t12273 = t10189*t10553
      t12275 = t10169*t10458
      t12277 = t9666*t10974
      t12279 = t10261*t10175
      t12281 = t9718*t9679
      t12283 = t9833*t9786
      t12285 = t12251/1530-12.D0/425.D0*t12253-2.D0/765.D0*t12255+6.D0/4
     #25.D0*t12257-t12259/3060+3.D0/850.D0*t9778*t9693+2.D0/425.D0*t1226
     #4-t12266/15300-t12269/3600-29.D0/10200.D0*t12271-13.D0/10200.D0*t1
     #2273-59.D0/15300.D0*t12275+43.D0/30600.D0*t12277+3.D0/1700.D0*t122
     #79-t12281/1700-t12283/13600
      wI62pc1 = t10534+t10937+t10896+t10855+t10699+t11022+t10981+t10488+
     #t11152+t11111+t10655+t11068+t10444+t11269+t11189+t10621+t11230+t11
     #374+t11339+t11305+t10398+t11481+t11446+t11408+t11632+t11597+t11560
     #+t11522+t11746+t11708+t11671+t9707+t11864+t11826+t11789+t11973+t11
     #936+t11901+t9832+t9765+t12076+t12042+t12214+t12008+t12178+t9921+t1
     #2144+t12110+t9880+t12285+t12250+t10005+t9964+t10155+t10104+t10057+
     #t10815+t10575+t10256+t10205+t10776+t10342+t10739+t10294
      t12292 = t3426*t3452
      t12293 = t12292*t6938
      t12294 = t6932*t177*t12293
      t12296 = t3452*t6937
      t12297 = t12296*t3609
      t12298 = t7148*t12297
      t12300 = t6932*t5930
      t12301 = t12296*t3692
      t12302 = t12300*t12301
      t12304 = t6932*t5899
      t12305 = t12292*t6949
      t12306 = t12304*t12305
      t12308 = t5*cw2
      t12309 = t6932*t12308
      t12310 = t12296*t7216
      t12311 = t12309*t12310
      t12313 = t12292*t7137
      t12314 = t12300*t12313
      t12316 = t12292*t6953
      t12317 = t12304*t12316
      t12319 = cw2*t3426
      t12320 = t6932*t12319
      t12321 = t12296*t6258
      t12322 = t12320*t12321
      t12324 = t6932*t5252
      t12325 = t12296*t5992
      t12326 = t12324*t12325
      t12328 = t12296*t1093
      t12329 = t12324*t12328
      t12331 = t6932*t3452
      t12332 = t12331*t7019
      t12334 = t12292*t7176
      t12335 = t12300*t12334
      t12337 = -t12294/16+7.D0/96.D0*t12298+t12302/192-t12306/3072-5.D0/
     #16.D0*t12311-5.D0/1536.D0*t12314-t12317/1536+5.D0/64.D0*t12322+47.
     #D0/768.D0*t12326+47.D0/384.D0*t12329-47.D0/96.D0*t12332-5.D0/768.D
     #0*t12335
      t12338 = t12292*t7082
      t12339 = t12300*t12338
      t12341 = cw2*t3452
      t12342 = t6932*t12341
      t12343 = t12342*t7059
      t12345 = s*t3452
      t12346 = t6932*t12345
      t12347 = t6949*t3609
      t12348 = t12346*t12347
      t12350 = t12331*t7059
      t12352 = t6932*t5299
      t12353 = t12296*t4025
      t12354 = t12352*t12353
      t12356 = t6953*t3609
      t12357 = t12346*t12356
      t12359 = t6932*cw2
      t12360 = t12359*t12293
      t12362 = t12342*t7019
      t12364 = t12296*t9098
      t12365 = t12309*t12364
      t12367 = t12292*t6987
      t12368 = t7148*t12367
      t12370 = t12296*t3795
      t12371 = t12352*t12370
      t12373 = t6932*t8548
      t12374 = t12296*t7067
      t12375 = t12373*t12374
      t12376 = t12375/4
      t12377 = -5.D0/1536.D0*t12339+43.D0/96.D0*t12343+7.D0/96.D0*t12348
     #-47.D0/192.D0*t12350+5.D0/384.D0*t12354+7.D0/48.D0*t12357+5.D0/128
     #.D0*t12360+43.D0/48.D0*t12362-5.D0/32.D0*t12365+47.D0/768.D0*t1236
     #8+5.D0/192.D0*t12371+t12376
      t12379 = t6932*t5910
      t12380 = t12379*t12316
      t12382 = t12379*t12305
      t12384 = t6932*t175
      t12385 = t12384*t12367
      t12387 = t6932*t5264
      t12388 = t12296*t720
      t12389 = t12387*t12388
      t12391 = t6932*t5936
      t12392 = t12391*t12301
      t12394 = t7129*t12313
      t12396 = t12292*t7183
      t12397 = t12384*t12396
      t12400 = t6932*t177*t3426
      t12401 = t12400*t12321
      t12403 = t168*t3452
      t12404 = t12403*t7103
      t12405 = t7129*t12404
      t12407 = t12296*t3970
      t12408 = t12391*t12407
      t12410 = t3453*t7103
      t12411 = t7129*t12410
      t12413 = t12296*t1096
      t12414 = t12387*t12413
      t12416 = -t12380/768-t12382/1536-43.D0/384.D0*t12385-7.D0/192.D0*t
     #12389+t12392/384-7.D0/384.D0*t12394-43.D0/192.D0*t12397-t12401/8+t
     #12405/384+t12408/768+t12411/768-7.D0/384.D0*t12414
      t12417 = t12296*t6254
      t12418 = t12320*t12417
      t12420 = t12400*t12417
      t12422 = t12292*t6988
      t12423 = t12384*t12422
      t12425 = t12373*t12364
      t12426 = t12425/4
      t12427 = t12331*t7111
      t12429 = t12373*t12310
      t12430 = t12429/2
      t12431 = t12300*t12407
      t12433 = t12309*t12374
      t12435 = t12352*t12297
      t12437 = t168*t3426
      t12438 = t12437*t12296
      t12439 = t7141*t12438
      t12441 = t1387*t3426
      t12442 = t12441*t12296
      t12443 = t7141*t12442
      t12445 = t12296*t3999
      t12446 = t12359*t12445
      t12448 = 5.D0/128.D0*t12418-t12420/16-43.D0/384.D0*t12423+t12426-4
     #7.D0/192.D0*t12427+t12430+t12431/384-5.D0/32.D0*t12433+5.D0/384.D0
     #*t12435-t12439/1536-t12443/3072+43.D0/96.D0*t12446
      wI63pc1 = t12337+t12377+t12416+t12448
      wI64pc1 = wV26
      t12450 = t5862*t168
      t12451 = t7396*t12450
      t12452 = 1/t5852
      t12456 = 1/(t6118-t6021-2*t2109+2*t6562+t2759-t5903)
      t12457 = t12452*t12456
      t12461 = t12451*t12457*t3452*t1832*t54
      t12463 = t5862*t12452
      t12465 = t7396*t12463*t12456
      t12466 = t3452*Mzsq
      t12467 = t7403*t5851
      t12468 = t12466*t12467
      t12469 = t12465*t12468
      t12471 = t1119*t5851
      t12472 = t12466*t12471
      t12473 = t12465*t12472
      t12476 = t7396*t12450*t12452
      t12477 = t12456*t3452
      t12478 = t12477*t842
      t12480 = t12476*t12478*t6462
      t12482 = t12477*t6475
      t12483 = t12476*t12482
      t12485 = t3452*t842
      t12486 = t12485*t6462
      t12487 = t12465*t12486
      t12489 = t12477*Mzsq
      t12491 = t12476*t12489*t12471
      t12493 = t12477*s
      t12494 = t6254*t5851
      t12496 = t12476*t12493*t12494
      t12498 = t12477*t6462
      t12499 = t12476*t12498
      t12501 = t168*t12452
      t12502 = t7396*t12501
      t12503 = t12477*t6468
      t12504 = t12502*t12503
      t12505 = t12504/40
      t12506 = t1387*t12452
      t12508 = t7396*t12506*t12456
      t12509 = t3452*t1811
      t12510 = t6134*s
      t12511 = t12509*t12510
      t12512 = t12508*t12511
      t12514 = t12477*t6465
      t12515 = t12502*t12514
      t12516 = t12515/40
      t12517 = t12502*t12482
      t12518 = t12517/40
      t12519 = t5869*s
      t12520 = t12509*t12519
      t12521 = t12508*t12520
      t12524 = t53*t5851*t54
      t12525 = t12466*t12524
      t12526 = t12508*t12525
      t12528 = t12501*t12456
      t12529 = t7396*t12528
      t12531 = t12485*t59*t6562
      t12532 = t12529*t12531
      t12533 = t12532/40
      t12534 = t12502*t12498
      t12535 = t12534/4
      t12537 = t12485*t53*t5926
      t12538 = t12529*t12537
      t12539 = t12538/20
      t12540 = t12461/640+t12469/160+t12473/160-t12480/64-t12483/640+t12
     #487/16-t12491/640-t12496/640-t12499/64-t12505-t12512/160+t12516+t1
     #2518-t12521/80+t12526/80+t12533+t12535-t12539
      t12541 = t3452*CsTh
      t12543 = t12541*t53*t5926
      t12544 = t12529*t12543
      t12545 = t12544/10
      t12546 = t12477*t6113
      t12547 = t12502*t12546
      t12548 = t12547/40
      t12549 = t3452*t62
      t12550 = t5869*t54
      t12551 = t12549*t12550
      t12552 = t12508*t12551
      t12554 = t12508*t12486
      t12557 = t3452*t53*t5926
      t12558 = t12508*t12557
      t12560 = t12541*t6475
      t12561 = t12508*t12560
      t12564 = t12477*t7465*t5894
      t12565 = t12502*t12564
      t12566 = t12565/20
      t12567 = t12508*t12472
      t12569 = t12508*t12468
      t12572 = t12466*t842*t12524
      t12573 = t12529*t12572
      t12574 = t12573/20
      t12576 = t12477*t7578*t5894
      t12577 = t12502*t12576
      t12578 = t12577/40
      t12580 = t12466*CsTh*t12524
      t12581 = t12529*t12580
      t12582 = t12581/10
      t12583 = t7396*t12506
      t12584 = t12583*t12546
      t12586 = t12583*t12503
      t12589 = t12477*t60*t5894
      t12590 = t12502*t12589
      t12591 = t12590/20
      t12592 = t12477*t6472
      t12593 = t12502*t12592
      t12594 = t12593/4
      t12595 = t12477*t6469
      t12596 = t12502*t12595
      t12597 = t12596/20
      t12599 = t12477*t5992*t5894
      t12600 = t12502*t12599
      t12601 = t12600/40
      t12602 = -t12545-t12548-t12552/8-t12554/16+t12558/80-t12561/80-t12
     #566-t12567/160-t12569/160-t12574-t12578-t12582+t12584/160+t12586/1
     #60+t12591-t12594-t12597-t12601
      t12604 = t7521*t5894
      t12605 = t12466*t12604
      t12606 = t12508*t12605
      t12608 = t7648*t5894
      t12609 = t12466*t12608
      t12610 = t12508*t12609
      t12612 = t12466*t6047
      t12613 = t12508*t12612
      t12615 = t12466*t5941
      t12616 = t12508*t12615
      t12618 = t7396*t12463
      t12619 = t12618*t12514
      t12621 = t12618*t12498
      t12623 = t12477*t6114
      t12624 = t12502*t12623
      t12625 = t12624/20
      t12626 = t53*t5894
      t12627 = t12626*s
      t12628 = t12541*t12627
      t12629 = t12508*t12628
      t12631 = t12477*t12627
      t12632 = t12502*t12631
      t12633 = t12632/20
      t12634 = t12618*t12482
      t12636 = t12583*t12564
      t12638 = t12583*t12589
      t12640 = t12477*t62
      t12641 = t7521*t1070
      t12643 = t12476*t12640*t12641
      t12645 = t12583*t12576
      t12647 = t12485*t12627
      t12648 = t12508*t12647
      t12650 = t59*t5894
      t12651 = t12477*t12650
      t12652 = t12502*t12651
      t12653 = t12652/40
      t12654 = t12477*t6425
      t12655 = t12502*t12654
      t12656 = t12655/40
      t12657 = t7648*t1070
      t12659 = t12476*t12640*t12657
      t12661 = -t12606/40-t12610/80+t12613/16+t12616/8+t12619/160+t12621
     #/16-t12625-t12629/40+t12633+t12634/160+t12636/80-t12638/80-t12643/
     #160+t12645/160-t12648/80-t12653-t12656-t12659/320
      t12662 = t12583*t12592
      t12665 = t12465*t12549*t12641
      t12667 = t12583*t12595
      t12669 = t12583*t12599
      t12671 = t12583*t12514
      t12674 = t12465*t12549*t12657
      t12676 = t12476*t12514
      t12678 = t12583*t12631
      t12680 = t12583*t12623
      t12683 = t7396*t5862*t12528
      t12684 = t12683*t12580
      t12686 = t12583*t12654
      t12688 = t12583*t12498
      t12690 = t12583*t12482
      t12692 = t12683*t12543
      t12694 = t12508*t12531
      t12696 = t12508*t12572
      t12698 = t12529*t12560
      t12699 = t12698/20
      t12700 = t12508*t12580
      t12702 = t12662/16+t12665/40+t12667/80+t12669/160-t12671/160+t1267
     #4/80-t12676/640-t12678/80+t12680/80+t12684/160+t12686/160-t12688/1
     #6-t12690/160+t12692/160-t12694/160+t12696/80+t12699+t12700/40
      t12705 = t12529*t12472
      t12706 = t12705/40
      t12707 = t12529*t12551
      t12708 = t12707/2
      t12709 = t12529*t12486
      t12710 = t12709/4
      t12711 = t12529*t12468
      t12712 = t12711/40
      t12714 = t12541*t59*t6562
      t12715 = t12508*t12714
      t12717 = t12529*t12525
      t12718 = t12717/20
      t12720 = t59*t62*t54
      t12722 = t12465*t12485*t12720
      t12724 = t12529*t12557
      t12725 = t12724/20
      t12726 = t12508*t12543
      t12728 = t12529*t12511
      t12729 = t12728/40
      t12730 = t12508*t12537
      t12732 = t12529*t12520
      t12733 = t12732/20
      t12734 = t12529*t12605
      t12735 = t12734/10
      t12736 = t12529*t12615
      t12737 = t12736/2
      t12738 = t12529*t12609
      t12739 = t12738/20
      t12740 = t12477*t12720
      t12741 = t12618*t12740
      t12743 = t12529*t12612
      t12744 = t12743/4
      t12745 = t12529*t12628
      t12746 = t12745/10
      t12747 = t12706+t12708+t12710+t12712-t12715/80-t12718-t12722/160-t
     #12725+t12726/40+t12729+t12730/80+t12733+t12735-t12737+t12739-t1274
     #1/160-t12744+t12746
      t12748 = t12529*t12647
      t12749 = t12748/20
      t12750 = t12529*t12714
      t12751 = t12750/20
      t12752 = t12465*t12580
      t12754 = t12583*t12651
      t12756 = t842*t1832
      t12757 = t12756*t54
      t12758 = t12477*t12757
      t12759 = t12618*t12758
      t12761 = t12465*t12543
      t12763 = t12477*t6481
      t12764 = t12618*t12763
      t12766 = t12683*t12714
      t12768 = t12477*t6456
      t12769 = t12476*t12768
      t12771 = t12465*t12531
      t12774 = t12476*t12489*t12524
      t12776 = t12465*t12560
      t12778 = t12683*t12572
      t12780 = t12476*t12763
      t12782 = t12465*t12714
      t12784 = t12477*t6459
      t12785 = t12476*t12784
      t12787 = t12477*t6453
      t12788 = t12476*t12787
      t12790 = t12683*t12537
      t12792 = t12749+t12751-t12752/40+t12754/160-t12759/160-t12761/40-t
     #12764/160-t12766/320+t12769/32+t12771/160+t12774/320+t12776/80+t12
     #778/320+t12780/640+t12782/80+t12785/320+t12788/64+t12790/320
      t12794 = t53*t62
      t12795 = t12794*t5851
      t12797 = t12476*t12493*t12795
      t12799 = t12477*CsTh
      t12801 = t12476*t12799*t6475
      t12804 = t12476*t12799*t6462
      t12807 = t12476*t12799*t6465
      t12809 = t12465*t12551
      t12811 = t12465*t12520
      t12813 = t12465*t12525
      t12815 = t12465*t12557
      t12818 = t12476*t12489*t12467
      t12821 = t12618*t12477*t1979
      t12825 = t12451*t12457*t12509*t1070
      t12829 = t12451*t12457*t12549*t1754
      t12832 = t12618*t12477*t6179
      t12834 = t12618*t12784
      t12837 = t53*t1811*t54
      t12838 = t12477*t12837
      t12839 = t12476*t12838
      t12841 = t12477*t6438
      t12842 = t12476*t12841
      t12844 = t12794*t1070
      t12845 = t12477*t12844
      t12846 = t12476*t12845
      t12848 = t12618*t12768
      t12850 = t12797/320-t12801/320-t12804/32-t12807/320+t12809/8+t1281
     #1/80-t12813/80-t12815/80-t12818/640-t12821/160+t12825/64+t12829/64
     #0-t12832/160-t12834/80-t12839/320+t12842/320-t12846/320-t12848/8
      t12851 = t12618*t12787
      t12853 = t12618*t12845
      t12855 = t12618*t12838
      t12857 = t12476*t12740
      t12859 = t12465*t12537
      t12861 = t12465*t12572
      t12864 = t12618*t12477*t2501
      t12866 = t12618*t12841
      t12868 = t12476*t12758
      t12870 = t12683*t12531
      t12872 = t12465*t12511
      t12875 = t12476*t12799*t12837
      t12878 = t12465*t12541*t12837
      t12881 = t12465*t12485*t12837
      t12884 = t12476*t12478*t12837
      t12887 = t12465*t12541*t12720
      t12890 = t12476*t12799*t12720
      t12893 = t12476*t12478*t12720
      t12895 = -t12851/16+t12853/80+t12855/80+t12857/640-t12859/80-t1286
     #1/80-t12864/16-t12866/80+t12868/640-t12870/640+t12872/160-t12875/1
     #60+t12878/40+t12881/80-t12884/320-t12887/80+t12890/320+t12893/640
      wI65pc1 = t12540+t12602+t12661+t12702+t12747+t12792+t12850+t12895
      t12898 = t7992*t5936
      t12899 = t12292*t7994
      t12900 = t12898*t12899
      t12902 = t7992*t5274
      t12903 = t12296*t8186
      t12904 = t12902*t12903
      t12906 = t7992*t5930
      t12907 = t12906*t12301
      t12909 = t8305*t12297
      t12911 = t7992*t5252
      t12912 = t12296*t8198
      t12913 = t12911*t12912
      t12915 = t7992*t12319
      t12917 = t12915*t12296*t8212
      t12919 = t12296*t8257
      t12920 = t12902*t12919
      t12922 = t12906*t12407
      t12924 = t12296*t8251
      t12925 = t12911*t12924
      t12927 = t8342*t12313
      t12929 = t8229*t12442
      t12931 = t8342*t12404
      t12933 = t8342*t12410
      t12935 = t8229*t12438
      t12937 = t7992*t12292
      t12938 = t12937*t8248
      t12941 = t7992*t5299*t3426
      t12942 = t12941*t12924
      t12945 = t12915*t12296*t8275
      t12947 = t12296*t1852
      t12948 = t12915*t12947
      t12950 = -7.D0/8192.D0*t12900+7.D0/1024.D0*t12904-t12907/128+t1290
     #9/32-t12913/32+t12917/16+7.D0/2048.D0*t12920-t12922/256-t12925/64-
     #t12927/128-17.D0/8192.D0*t12929-t12931/256+17.D0/2048.D0*t12933+t1
     #2935/1024-t12938/32+25.D0/1024.D0*t12942+t12945/32+t12948/32
      t12951 = t12906*t12899
      t12953 = t7992*t175
      t12954 = t12953*t12367
      t12956 = t8305*t12367
      t12958 = t7992*t5299
      t12959 = t12296*t8093
      t12960 = t12958*t12959
      t12963 = t7992*t5*t3452
      t12964 = t7994*t3999
      t12965 = t12963*t12964
      t12966 = t12965/8
      t12967 = t12958*t12903
      t12970 = t7992*t12308*t3452
      t12971 = t12970*t8190
      t12972 = t12971/4
      t12973 = t12963*t8190
      t12974 = t12973/4
      t12975 = t12292*t8169
      t12976 = t12898*t12975
      t12978 = t12958*t12297
      t12980 = t7992*t5910
      t12981 = t12980*t12316
      t12983 = t12970*t8245
      t12984 = t12983/8
      t12985 = t7992*t12345
      t12986 = t12985*t12356
      t12988 = t7992*t3426
      t12989 = t12988*t12947
      t12991 = t7992*t12341
      t12992 = t12991*t7019
      t12993 = t12992/8
      t12994 = t7992*cw2
      t12995 = t12994*t12445
      t12996 = t12995/16
      t12997 = t12958*t12370
      t12999 = t12958*t12919
      t13001 = -23.D0/2048.D0*t12951+t12954/64-t12956/64+23.D0/512.D0*t1
     #2960+t12966+23.D0/256.D0*t12967-t12972+t12974-7.D0/8192.D0*t12976-
     #25.D0/512.D0*t12978+t12981/512-t12984+t12986/16-t12989/32-t12993-t
     #12996-25.D0/256.D0*t12997+23.D0/512.D0*t12999
      t13003 = t12963*t8245
      t13004 = t13003/8
      t13005 = t12958*t12353
      t13007 = t12991*t7059
      t13008 = t13007/16
      t13009 = t12985*t12347
      t13011 = t12953*t12396
      t13013 = t12902*t12959
      t13015 = t12898*t12301
      t13017 = t12292*t8308
      t13018 = t12906*t13017
      t13020 = t12953*t12422
      t13022 = t7992*t5264
      t13023 = t13022*t12388
      t13025 = t7992*t3452
      t13026 = t13025*t7059
      t13027 = t13026/16
      t13028 = t12906*t12975
      t13030 = t12898*t12407
      t13032 = t13022*t12413
      t13034 = t12296*t445
      t13035 = t12911*t13034
      t13037 = t7992*t5899
      t13038 = t13037*t12305
      t13040 = t7992*t12403
      t13041 = t13040*t7995
      t13043 = t12980*t12305
      t13045 = t13004-25.D0/512.D0*t13005-t13008+t13009/32+t13011/32+7.D
     #0/2048.D0*t13013+17.D0/1024.D0*t13015-23.D0/1024.D0*t13018+t13020/
     #64-t13023/64+t13027-23.D0/2048.D0*t13028+17.D0/2048.D0*t13030-t130
     #32/128-t13035/64-17.D0/8192.D0*t13038-25.D0/256.D0*t13041+t13043/1
     #024
      t13046 = t12906*t12313
      t13048 = t13040*t8412
      t13050 = t13025*t7019
      t13051 = t13050/8
      t13052 = t13040*t8441
      t13054 = t7992*t12296
      t13055 = t1865*t3609
      t13056 = t13054*t13055
      t13058 = t12911*t12328
      t13060 = t1224*t3609
      t13061 = t13054*t13060
      t13063 = t13025*t7995
      t13065 = t12906*t12334
      t13067 = t7992*t12308
      t13069 = t13067*t12296*t8129
      t13070 = t13069/8
      t13071 = t13037*t12316
      t13073 = t12937*t8195
      t13075 = t12906*t12338
      t13077 = t12941*t12912
      t13079 = t12898*t13017
      t13081 = t13025*t7111
      t13082 = t13081/16
      t13083 = t12911*t12325
      t13085 = t12292*t8055
      t13086 = t12958*t13085
      t13088 = 25.D0/2048.D0*t13046-25.D0/128.D0*t13048+t13051-25.D0/256
     #.D0*t13052+t13056/16-t13058/32+t13061/8+t13063/16+25.D0/1024.D0*t1
     #3065-t13070-17.D0/4096.D0*t13071-t13073/16+25.D0/2048.D0*t13075+25
     #.D0/512.D0*t13077-7.D0/4096.D0*t13079+t13082-t13083/64+25.D0/1024.
     #D0*t13086
      wI66pc1 = t12950+t13001+t13045+t13088
      t13091 = t8632*t8534*t12292
      t13093 = t8517*t3452
      t13094 = t13093*t5452
      t13095 = t8523*t13094
      t13097 = t3426*t842
      t13098 = t13093*t13097
      t13099 = t8523*t13098
      t13101 = t3452*t3409
      t13103 = t8632*t8644*t13101
      t13105 = t13093*t3692
      t13106 = t8732*t13105
      t13108 = t13093*t3970
      t13109 = t8732*t13108
      t13111 = t3426*Mzsq
      t13112 = t13093*t13111
      t13113 = t8621*t13112
      t13115 = t12292*t720
      t13116 = t8645*t13115
      t13118 = t12292*t1096
      t13119 = t8645*t13118
      t13121 = t3426*t62
      t13122 = t13093*t13121
      t13123 = t8544*t13122
      t13126 = t13093*t5452*t62
      t13127 = t8544*t13126
      t13130 = t13093*t13121*t842
      t13131 = t8544*t13130
      t13133 = t13091/64+t13095/32+t13099/64-t13103/16-t13106/8-t13109/1
     #6+t13113/16+t13116/8+t13119/16+5.D0/64.D0*t13123+5.D0/32.D0*t13127
     #+5.D0/64.D0*t13131
      t13134 = t3426*t1811
      t13135 = t13093*t13134
      t13136 = t8515*t13135
      t13138 = t12292*t6258
      t13139 = t8563*t13138
      t13141 = t12292*t6254
      t13142 = t8563*t13141
      t13145 = t13093*t3999
      t13146 = t8514*cw2*t13145
      t13148 = t12541*t3999
      t13149 = t8579*t13148
      t13151 = t12485*t3999
      t13152 = t8579*t13151
      t13155 = t13093*t3609
      t13156 = t8514*s*t13155
      t13157 = t13156/4
      t13158 = t12541*t3609
      t13159 = t8555*t13158
      t13160 = t13159/2
      t13161 = t12485*t3609
      t13162 = t8555*t13161
      t13163 = t13162/4
      t13164 = t13093*t9098
      t13165 = t8549*t13164
      t13167 = t13093*t7216
      t13168 = t8549*t13167
      t13170 = t13093*t7067
      t13171 = t8549*t13170
      t13173 = t13136/32+t13139/16+t13142/32-5.D0/16.D0*t13146-5.D0/8.D0
     #*t13149-5.D0/16.D0*t13152-t13157-t13160-t13163-t13165/8-t13168/4-t
     #13171/8
      wI67pc1 = t13133+t13173
      t13175 = t8847*t12292*t5992
      t13178 = t8878*t12292*t5293
      t13181 = t8847*t12292*t6325
      t13184 = t8878*t12403*t3692
      t13187 = t8878*t12403*t3970
      t13190 = t8878*t12292*t1113
      t13192 = t8831*t13093
      t13194 = t13192*t1603*t3999
      t13196 = t1387*t3409
      t13198 = t13192*t13196*t1093
      t13201 = t8992*t12341*t9098
      t13203 = t8847*t13158
      t13207 = t8847*t12403*t3609
      t13210 = t8837*t12292*t5316
      t13213 = t8878*t12292*t1324
      t13216 = t13192*t12319*t6254
      t13219 = t8847*t12292*t6287
      t13222 = t8847*t12403*t3795
      t13226 = t8847*t12292*t5387*t62
      t13229 = t13192*t13196*t5992
      t13232 = t8847*t3453*t3795
      t13235 = t8847*t3453*t4025
      t13238 = t8847*t12403*t4025
      t13240 = 11.D0/72.D0*t13207-5.D0/252.D0*t13210-13.D0/504.D0*t13213
     #-13.D0/504.D0*t13216-47.D0/1008.D0*t13219+11.D0/36.D0*t13222+37.D0
     #/8064.D0*t13226-37.D0/2016.D0*t13229-37.D0/504.D0*t13232-37.D0/100
     #8.D0*t13235+11.D0/72.D0*t13238
      t13243 = t8992*t12403*t7067
      t13246 = t8992*t12509*t3692
      t13249 = t8992*t12509*t3970
      t13252 = t8992*t12403*t7216
      t13255 = t13192*t1323*t3999
      t13257 = t8847*t13161
      t13260 = t8847*t12341*t3609
      t13263 = t8992*t12403*t9098
      t13266 = t8847*t12341*t3795
      t13269 = t8847*t12341*t4025
      t13273 = t8992*t12341*t7067
      t13276 = t13192*t6287*t3692
      t13279 = t13192*t6287*t3970
      t13282 = t8992*t12341*t7216
      t13286 = t8847*t12292*t5396*t62
      t13289 = t13192*t12437*t6258
      t13292 = t8859*t3453*t3999
      t13294 = t13093*t12437
      t13295 = t8852*t13294
      t13298 = t8852*t13093*t12441
      t13300 = t8831*t5
      t13301 = t13300*t13164
      t13304 = t8859*t12292*t3551
      t13306 = 13.D0/126.D0*t13273+47.D0/126.D0*t13276+47.D0/252.D0*t132
     #79+13.D0/63.D0*t13282-47.D0/1008.D0*t13286-37.D0/1008.D0*t13289-37
     #.D0/2016.D0*t13292-5.D0/504.D0*t13295+37.D0/8064.D0*t13298-10.D0/6
     #3.D0*t13301-13.D0/504.D0*t13304
      t13309 = t8890*t13122
      t13311 = t9038*t13112
      t13314 = t9038*t13093*t13196
      t13316 = t8852*t13098
      t13318 = t8852*t13094
      t13320 = t9038*t13105
      t13322 = t8859*t13151
      t13324 = t8859*t13148
      t13327 = t8859*t12403*t3999
      t13330 = t8859*t12292*t6283
      t13333 = t8859*t13141
      t13335 = t8859*t13138
      t13338 = t13093*t168*t3409
      t13339 = t9038*t13338
      t13341 = t8890*t13155
      t13344 = t8859*t12341*t3999
      t13346 = t9038*t13108
      t13349 = t8852*t13093*t3426
      t13352 = t8859*t12292*t1811
      t13355 = t8859*t13101*t62
      t13358 = t9038*t13093*t3409
      t13361 = t8847*t3453*t3609
      t13363 = 5.D0/126.D0*t13333+5.D0/63.D0*t13335+5.D0/126.D0*t13339-5
     #3.D0/252.D0*t13341+13.D0/63.D0*t13344-13.D0/504.D0*t13346+13.D0/20
     #16.D0*t13349+5.D0/126.D0*t13352-173.D0/504.D0*t13355-13.D0/504.D0*
     #t13358-37.D0/1008.D0*t13361
      t13367 = t8847*t12292*t5316*t62
      t13371 = t8847*t12292*t1603*t62
      t13374 = t8878*t3453*t3692
      t13376 = t8878*t13115
      t13380 = t8847*t12292*t1323*t62
      t13383 = t8878*t3453*t3970
      t13385 = t8878*t13118
      t13388 = t8878*t12292*t7
      t13391 = t8847*t12292*t3386
      t13394 = t13192*t12437*t6254
      t13398 = t8837*t12292*t5387
      t13401 = t8837*t12292*t5396
      t13404 = t8878*t12292*t5278
      t13407 = t8847*t12292*t1093
      t13410 = t8878*t12292*t6870
      t13413 = t8837*t12292*t5331
      t13415 = t5396*Mzsq
      t13417 = t8878*t12292*t13415
      t13420 = t13192*t12319*t6258
      t13424 = t8878*t12292*t5293*CsTh
      t13428 = t8878*t12292*t5293*t842
      t13432 = t8847*t12292*t6325*CsTh
      t13434 = 37.D0/8064.D0*t13398-5.D0/504.D0*t13401-11.D0/288.D0*t134
     #04+173.D0/1008.D0*t13407-11.D0/144.D0*t13410+37.D0/4032.D0*t13413-
     #11.D0/288.D0*t13417-13.D0/252.D0*t13420+37.D0/2016.D0*t13424+37.D0
     #/4032.D0*t13428+37.D0/4032.D0*t13432
      wI68pc1 = -37.D0/2016.D0*t13330-173.D0/252.D0*t13324+47.D0/252.D0*
     #t13327-173.D0/504.D0*t13322+13.D0/1008.D0*t13318-13.D0/252.D0*t133
     #20+13.D0/2016.D0*t13316+173.D0/2016.D0*t13309+53.D0/1008.D0*t13311
     #-37.D0/2016.D0*t13314+t13306+13.D0/126.D0*t13269+13.D0/63.D0*t1326
     #6+13.D0/126.D0*t13260+37.D0/504.D0*t13263-10.D0/63.D0*t13249+37.D0
     #/252.D0*t13252+13.D0/63.D0*t13255-53.D0/252.D0*t13257-20.D0/63.D0*
     #t13246+37.D0/504.D0*t13243+t13240-37.D0/1008.D0*t13198+13.D0/126.D
     #0*t13201-53.D0/126.D0*t13203+26.D0/63.D0*t13194-13.D0/252.D0*t1319
     #0+37.D0/8064.D0*t13181+5.D0/63.D0*t13184+5.D0/126.D0*t13187+37.D0/
     #4032.D0*t13178+173.D0/2016.D0*t13175+t13434-37.D0/2016.D0*t13394-1
     #3.D0/252.D0*t13391-37.D0/1008.D0*t13374+53.D0/504.D0*t13376-13.D0/
     #252.D0*t13380-37.D0/2016.D0*t13383+53.D0/1008.D0*t13385-13.D0/504.
     #D0*t13388-13.D0/126.D0*t13371-47.D0/504.D0*t13367+t13363
      t13437 = t3452*t9433
      t13438 = t5*t5851
      t13440 = t9486*t13437*t13438
      t13441 = t13440/8
      t13444 = t9486*t13437*t13438*CsTh
      t13445 = t13444/4
      t13448 = t9486*t13437*t13438*t842
      t13449 = t13448/8
      t13450 = t5862*t3452
      t13452 = t9486*t13450*t9434
      t13456 = t9486*t13450*t9434*CsTh
      t13460 = t9486*t13450*t9434*t842
      t13464 = t9616*wV8*t5862*t13437
      t13468 = t9465*t13450*t9433*CsTh
      t13472 = t9465*t13450*t9433*t842
      t13474 = t3386*t211
      t13478 = t13474*wV8*t3452*t9433*t5851
      t13479 = t13478/4
      t13481 = t9535*t13437*t5869
      t13482 = t13481/2
      t13484 = t9535*t13437*t6134
      t13485 = t13484/4
      wI69pc1 = t13441+t13445+t13449-t13452/16-t13456/8-t13460/16-t13464
     #/32-t13468/16-t13472/32+t13479+t13482+t13485
      t13486 = t4468*t12441
      t13487 = t12296*t9825
      t13488 = t13486*t13487
      t13490 = t10065*t12292
      t13491 = t13490*t10298
      t13493 = t4892*t12341
      t13494 = t13493*t10857
      t13495 = 39.D0/1384.D0*t13494
      t13496 = t4468*t12403
      t13497 = t13496*t10668
      t13499 = t10065*t12296
      t13500 = t13499*t10860
      t13502 = t5590*t12296
      t13503 = t13502*t10860
      t13505 = t4892*t12403
      t13506 = t13505*t10857
      t13508 = t13505*t10869
      t13510 = t4892*t12296
      t13511 = t1096*t10348
      t13512 = t13510*t13511
      t13514 = t13505*t11514
      t13516 = t9682*t12403
      t13517 = t13516*t9984
      t13519 = -5.D0/5536.D0*t13488+13.D0/1384.D0*t13491-t13495-9.D0/138
     #4.D0*t13497-125.D0/2768.D0*t13500+77.D0/5536.D0*t13503-709.D0/1107
     #2.D0*t13506+21.D0/5536.D0*t13508+69.D0/5536.D0*t13512+21.D0/5536.D
     #0*t13514+37.D0/2768.D0*t13517
      t13521 = t9682*t177*t3452
      t13522 = t13521*t10931
      t13523 = 9.D0/346.D0*t13522
      t13524 = t4460*t12437
      t13526 = t13524*t12296*t10612
      t13528 = t9682*t12341
      t13529 = t13528*t10931
      t13530 = 27.D0/1384.D0*t13529
      t13532 = t13524*t12296*t10335
      t13534 = t13528*t10421
      t13535 = 27.D0/692.D0*t13534
      t13536 = t13521*t10421
      t13537 = 9.D0/173.D0*t13536
      t13538 = t13496*t10920
      t13540 = t4468*t12437
      t13541 = t12296*t10525
      t13542 = t13540*t13541
      t13544 = t12296*t10121
      t13545 = t13486*t13544
      t13547 = t4468*t12319
      t13548 = t12296*t10222
      t13549 = t13547*t13548
      t13551 = t4468*t12292
      t13552 = t13551*t11079
      t13554 = t211*t3452
      t13555 = t13554*t3426
      t13556 = t13555*t7621
      t13558 = -t13523+9.D0/5536.D0*t13526+t13530+9.D0/11072.D0*t13532+t
     #13535-t13537-9.D0/2768.D0*t13538+3.D0/11072.D0*t13542-3.D0/44288.D
     #0*t13545-13.D0/692.D0*t13549+147.D0/5536.D0*t13552-479.D0/22144.D0
     #*t13556
      t13560 = t13555*t7674
      t13562 = t13555*t7942
      t13564 = t4460*t3426
      t13566 = t13564*t12296*spf14
      t13568 = t13554*t6937
      t13569 = t13568*t9734
      t13571 = t6967*t1070
      t13572 = t13554*t13571
      t13574 = t4468*t3452
      t13575 = t13574*t10319
      t13577 = t4892*t12292
      t13578 = t13577*t11102
      t13580 = t4892*t12437
      t13581 = t12296*t10252
      t13582 = t13580*t13581
      t13584 = t4892*t12319
      t13585 = t13584*t13581
      t13587 = t13540*t13487
      t13589 = t4468*t3453
      t13590 = t13589*t10012
      t13592 = -13.D0/11072.D0*t13560-615.D0/22144.D0*t13562-39.D0/22144
     #.D0*t13566+89.D0/2768.D0*t13569-13.D0/11072.D0*t13572+39.D0/5536.D
     #0*t13575-93.D0/5536.D0*t13578+57.D0/5536.D0*t13582+13.D0/1384.D0*t
     #13585-51.D0/11072.D0*t13587-37.D0/5536.D0*t13590
      t13593 = t13547*t13487
      t13595 = t5590*t12292
      t13596 = t13595*t9811
      t13598 = t13490*t9811
      t13600 = t12296*t11062
      t13601 = t13540*t13600
      t13603 = t12296*t11177
      t13604 = t13486*t13603
      t13606 = t13502*t11248
      t13608 = t12296*t11058
      t13609 = t13540*t13608
      t13611 = t13540*t13603
      t13613 = t9682*t12296
      t13614 = t1852*t11247
      t13615 = t13613*t13614
      t13616 = 4.D0/173.D0*t13615
      t13618 = t13580*t12296*t11018
      t13620 = t13502*t11689
      t13623 = t13554*t6953*t11248
      t13625 = -13.D0/1384.D0*t13593+37.D0/11072.D0*t13596+13.D0/2768.D0
     #*t13598-t13601/2768+17.D0/22144.D0*t13604+25.D0/1384.D0*t13606-t13
     #609/1384-t13611/2768-t13616-25.D0/5536.D0*t13618+25.D0/692.D0*t136
     #20-4.D0/173.D0*t13623
      t13628 = t13528*t11015
      t13629 = 4.D0/173.D0*t13628
      t13630 = t4892*t3453
      t13631 = t13630*t11527
      t13633 = t13502*t11025
      t13636 = t13554*t6949*t11248
      t13638 = t13505*t10131
      t13639 = 3.D0/2768.D0*t13638
      t13640 = t13630*t10054
      t13642 = t13630*t10131
      t13644 = t13555*t10947
      t13646 = t13554*t1387
      t13647 = t13646*t7416
      t13649 = t13568*t11248
      t13651 = t13555*t11298
      t13653 = t13629-17.D0/5536.D0*t13631+25.D0/1384.D0*t13633-2.D0/173
     #.D0*t13636-t13639+3.D0/5536.D0*t13640+3.D0/11072.D0*t13642+5.D0/69
     #2.D0*t13644+243.D0/11072.D0*t13647-2.D0/173.D0*t13649+t13651/173
      t13654 = t4468*t3426
      t13655 = t12296*t9719
      t13656 = t13654*t13655
      t13658 = t13574*t11962
      t13660 = t13574*t9762
      t13662 = t4892*t3426
      t13663 = t12296*t9757
      t13664 = t13662*t13663
      t13666 = t12296*t9933
      t13667 = t10703*t13666
      t13669 = t13568*t11140
      t13671 = t12296*t10552
      t13672 = t13564*t13671
      t13674 = t12296*t10556
      t13675 = t13564*t13674
      t13677 = t4468*t1387
      t13678 = t12296*t9739
      t13679 = t13677*t13678
      t13681 = t10712*t13666
      t13683 = t12296*t9810
      t13684 = t10703*t13683
      t13686 = t10712*t13683
      t13688 = -423.D0/22144.D0*t13656+39.D0/5536.D0*t13658+39.D0/2768.D
     #0*t13660-89.D0/11072.D0*t13664-43.D0/11072.D0*t13667+89.D0/1384.D0
     #*t13669-39.D0/22144.D0*t13672-39.D0/11072.D0*t13675+317.D0/22144.D
     #0*t13679-37.D0/11072.D0*t13681+13.D0/2768.D0*t13684+37.D0/11072.D0
     #*t13686
      t13690 = t12296*t9834
      t13691 = t13677*t13690
      t13693 = t4468*t168
      t13694 = t13693*t13678
      t13696 = t13646*t9740
      t13698 = t13568*t10860
      t13700 = t10065*t3452
      t13701 = t13700*t10545
      t13702 = 13.D0/692.D0*t13701
      t13703 = t12296*t9817
      t13704 = t13662*t13703
      t13706 = t12296*t9839
      t13707 = t13654*t13706
      t13709 = t9682*t3452
      t13710 = t13709*t11551
      t13712 = t12296*t10194
      t13713 = t13564*t13712
      t13715 = t13646*t10545
      t13717 = t13693*t13690
      t13719 = t13568*t9877
      t13721 = -37.D0/5536.D0*t13691-205.D0/11072.D0*t13694-37.D0/11072.
     #D0*t13696+89.D0/2768.D0*t13698-t13702-93.D0/11072.D0*t13704+147.D0
     #/11072.D0*t13707+5.D0/173.D0*t13710+13.D0/5536.D0*t13713+37.D0/110
     #72.D0*t13715+5.D0/346.D0*t13717+93.D0/2768.D0*t13719
      t13722 = t13568*t11055
      t13724 = t4460*t168
      t13725 = t12292*t10217
      t13726 = t13724*t13725
      t13728 = t168*t6937
      t13730 = t13555*t13728*t1072
      t13733 = t13555*t13728*t1077
      t13735 = t13574*t10012
      t13737 = t13574*t12169
      t13739 = t5590*t3452
      t13740 = t13739*t10545
      t13742 = t4892*t3452
      t13743 = t13742*t12119
      t13745 = t13700*t7799
      t13746 = 13.D0/692.D0*t13745
      t13747 = t13568*t5549
      t13749 = t13555*t7903
      t13751 = t13555*t7900
      t13753 = 93.D0/1384.D0*t13722-5.D0/1384.D0*t13726-5.D0/2768.D0*t13
     #730+583.D0/22144.D0*t13733-13.D0/1384.D0*t13735-13.D0/692.D0*t1373
     #7-57.D0/2768.D0*t13740-147.D0/2768.D0*t13743-t13746+615.D0/2768.D0
     #*t13747-615.D0/11072.D0*t13749-479.D0/11072.D0*t13751
      t13757 = t13739*t7059
      t13759 = t13739*t7595
      t13761 = t1387*t6937
      t13763 = t13555*t13761*t1071
      t13765 = t13646*t7595
      t13767 = t13540*t13706
      t13769 = t13555*t7612
      t13771 = t13555*t7625
      t13773 = t13646*t7547
      t13775 = t13739*t7799
      t13777 = t13568*t5552
      t13779 = t13700*t7059
      t13781 = 5.D0/692.D0*t13757-143.D0/5536.D0*t13759-243.D0/22144.D0*
     #t13763+243.D0/11072.D0*t13765-51.D0/11072.D0*t13767-479.D0/22144.D
     #0*t13769-615.D0/22144.D0*t13771+243.D0/5536.D0*t13773-583.D0/5536.
     #D0*t13775+615.D0/5536.D0*t13777-71.D0/692.D0*t13779
      t13782 = cw2*t6937
      t13784 = t13555*t13782*t1077
      t13787 = t13555*t13782*t1072
      t13790 = t13555*t13728*t1071
      t13793 = t13555*t13761*t1119
      t13795 = t13739*t7547
      t13797 = t13739*t7019
      t13799 = t13700*t7019
      t13802 = t13555*t13728*t1119
      t13804 = t10712*t12321
      t13806 = t10703*t12321
      t13808 = t13739*t10496
      t13809 = 9.D0/2768.D0*t13808
      t13810 = t12296*t10310
      t13811 = t13486*t13810
      t13813 = 13.D0/2768.D0*t13784+71.D0/2768.D0*t13787+143.D0/11072.D0
     #*t13790-243.D0/44288.D0*t13793-143.D0/2768.D0*t13795+5.D0/346.D0*t
     #13797-71.D0/346.D0*t13799+143.D0/22144.D0*t13802-37.D0/5536.D0*t13
     #804-t13806/346+t13809-3.D0/22144.D0*t13811
      t13815 = t13595*t7903
      t13817 = t13595*t7625
      t13819 = t10712*t12417
      t13821 = t13646*t7059
      t13823 = t13646*t7019
      t13825 = t13580*t13703
      t13827 = t13493*t12119
      t13828 = 13.D0/346.D0*t13827
      t13829 = t13646*t10496
      t13831 = t13646*t10489
      t13833 = t13547*t13706
      t13835 = t13584*t13703
      t13837 = 583.D0/11072.D0*t13815+583.D0/22144.D0*t13817-37.D0/11072
     #.D0*t13819-37.D0/11072.D0*t13821-37.D0/5536.D0*t13823+57.D0/11072.
     #D0*t13825+t13828-9.D0/11072.D0*t13829-9.D0/11072.D0*t13831-13.D0/1
     #384.D0*t13833+13.D0/2768.D0*t13835
      t13838 = t13739*t10489
      t13839 = 9.D0/2768.D0*t13838
      t13840 = t9682*cw2
      t13841 = t13840*t12364
      t13843 = t10703*t12417
      t13845 = t4892*t1387
      t13846 = t13845*t12367
      t13848 = t9682*t168
      t13849 = t13848*t12364
      t13851 = t13709*t12220
      t13853 = t13709*t12241
      t13855 = t13630*t12119
      t13857 = t13516*t11551
      t13859 = t13630*t9836
      t13861 = t13555*t9934
      t13863 = t13709*t9984
      t13865 = t13839+t13841/173-t13843/692+37.D0/44288.D0*t13846+37.D0/
     #2768.D0*t13849-5.D0/173.D0*t13851-10.D0/173.D0*t13853+5.D0/1384.D0
     #*t13855-37.D0/2768.D0*t13857+5.D0/692.D0*t13859+5.D0/346.D0*t13861
     #-5.D0/173.D0*t13863
      t13868 = t13595*t7900
      t13870 = t13490*t7900
      t13872 = t12296*t10145
      t13873 = t13564*t13872
      t13875 = t4460*t1387
      t13876 = t12292*t9957
      t13877 = t13875*t13876
      t13879 = t13724*t13876
      t13881 = t13739*t9740
      t13883 = t13700*t9740
      t13885 = t13742*t9720
      t13887 = t13555*t10175
      t13889 = t13490*t7903
      t13891 = t13595*t7612
      t13893 = -5.D0/1384.D0*t13868+71.D0/1384.D0*t13870+13.D0/2768.D0*t
     #13873-317.D0/88576.D0*t13877+205.D0/44288.D0*t13879+77.D0/5536.D0*
     #t13881-125.D0/2768.D0*t13883+423.D0/5536.D0*t13885+5.D0/692.D0*t13
     #887+13.D0/1384.D0*t13889-5.D0/2768.D0*t13891
      t13894 = t13555*t9811
      t13896 = t13555*t10298
      t13899 = t13693*t12296*t10348
      t13901 = t12296*t10371
      t13902 = t10703*t13901
      t13904 = t12296*t10379
      t13905 = t13654*t13904
      t13907 = t12296*t10363
      t13908 = t13662*t13907
      t13910 = t4460*t12441
      t13911 = t13910*t13712
      t13913 = t13700*t10589
      t13915 = t13524*t13712
      t13917 = t13528*t11551
      t13918 = 13.D0/692.D0*t13917
      t13919 = t13496*t12169
      t13922 = t9544*t3426*t13901
      t13924 = -5.D0/692.D0*t13894-5.D0/346.D0*t13896-9.D0/2768.D0*t1389
     #9-27.D0/5536.D0*t13902-69.D0/22144.D0*t13905+93.D0/22144.D0*t13908
     #+37.D0/22144.D0*t13911-3.D0/692.D0*t13913-5.D0/1384.D0*t13915-t139
     #18+5.D0/173.D0*t13919+9.D0/1384.D0*t13922
      t13926 = t13875*t13725
      t13928 = t13742*t11514
      t13930 = t13568*t10875
      t13932 = t13493*t9870
      t13933 = 13.D0/346.D0*t13932
      t13934 = t13505*t9836
      t13936 = t3692*t9839
      t13937 = t13510*t13936
      t13939 = t13502*t11055
      t13941 = t13490*t7612
      t13943 = t13490*t7625
      t13945 = t13568*t10349
      t13948 = t13724*t12292*t10583
      t13950 = t13499*t9877
      t13951 = 13.D0/692.D0*t13950
      t13952 = 37.D0/22144.D0*t13926+69.D0/5536.D0*t13928-93.D0/5536.D0*
     #t13930+t13933+51.D0/1384.D0*t13934-147.D0/1384.D0*t13937-57.D0/138
     #4.D0*t13939+71.D0/2768.D0*t13941+13.D0/2768.D0*t13943-93.D0/2768.D
     #0*t13945+9.D0/11072.D0*t13948-t13951
      t13953 = t13505*t9870
      t13955 = t3970*t9839
      t13956 = t13510*t13955
      t13958 = t13502*t9877
      t13960 = t13499*t5549
      t13961 = 13.D0/346.D0*t13960
      t13962 = t13502*t5549
      t13964 = t13505*t12119
      t13966 = t5598*t12296
      t13967 = t13966*t11055
      t13969 = t13568*t4253
      t13971 = t13528*t9807
      t13972 = 13.D0/346.D0*t13971
      t13973 = t13516*t9807
      t13975 = t6258*t9834
      t13976 = t13613*t13975
      t13978 = t13630*t9870
      t13980 = 51.D0/2768.D0*t13953-147.D0/2768.D0*t13956-57.D0/2768.D0*
     #t13958-t13961-583.D0/2768.D0*t13962+51.D0/2768.D0*t13964+37.D0/553
     #6.D0*t13967+13.D0/2768.D0*t13969-t13972-37.D0/1384.D0*t13973+10.D0
     #/173.D0*t13976+5.D0/1384.D0*t13978
      t13985 = t13496*t10012
      t13987 = t13966*t9877
      t13989 = t13496*t9762
      t13991 = t13966*t11140
      t13993 = t13516*t9769
      t13995 = t13516*t9966
      t13997 = t13499*t5552
      t13998 = 13.D0/692.D0*t13997
      t13999 = t13502*t5552
      t14002 = t9682*t12341*t6937
      t14003 = t14002*t11850
      t14004 = 4.D0/173.D0*t14003
      t14006 = t13564*t12296*t9794
      t14008 = t13577*t11541
      t14010 = 5.D0/346.D0*t13985+37.D0/11072.D0*t13987-205.D0/5536.D0*t
     #13989-37.D0/5536.D0*t13991+37.D0/1384.D0*t13993-37.D0/2768.D0*t139
     #95-t13998-583.D0/5536.D0*t13999+t14004+13.D0/5536.D0*t14006+t14008
     #/346
      t14011 = t13490*t11298
      t14013 = t13568*t7022
      t14015 = t13613*t11850
      t14016 = 4.D0/173.D0*t14015
      t14017 = t13630*t9973
      t14019 = t6258*t9739
      t14020 = t13613*t14019
      t14022 = t13966*t10860
      t14024 = t13505*t10994
      t14026 = t10277*t12296
      t14027 = t14026*t11853
      t14029 = t13505*t11527
      t14031 = t13568*t10664
      t14033 = t13568*t3994
      t14035 = t12296*t1811
      t14036 = t10277*t14035
      t14038 = -t14011/173+479.D0/5536.D0*t14013-t14016-37.D0/11072.D0*t
     #14017-10.D0/173.D0*t14020-37.D0/11072.D0*t14022+t14024/692+2.D0/17
     #3.D0*t14027+t14029/692-93.D0/5536.D0*t14031+615.D0/5536.D0*t14033+
     #5.D0/692.D0*t14036
      t14040 = t13568*t4121
      t14042 = t13630*t11312
      t14044 = t4892*t12441
      t14045 = t14044*t12325
      t14047 = t13739*t7111
      t14049 = t13568*t9821
      t14051 = t13516*t10018
      t14053 = t13528*t9769
      t14055 = t13496*t11962
      t14057 = t13630*t10994
      t14059 = t13502*t9783
      t14060 = 9.D0/1384.D0*t14059
      t14061 = t13577*t11293
      t14063 = 13.D0/1384.D0*t14040-17.D0/2768.D0*t14042+37.D0/44288.D0*
     #t14045+5.D0/692.D0*t14047+93.D0/2768.D0*t14049+37.D0/2768.D0*t1405
     #1+43.D0/1384.D0*t14053-205.D0/11072.D0*t14055-17.D0/5536.D0*t14057
     #+t14060+t14061/173
      t14064 = t13490*t11606
      t14066 = t14026*t11126
      t14069 = t4892*t12437*t3452
      t14070 = t14069*t11293
      t14072 = t14069*t11019
      t14074 = t14002*t11886
      t14075 = 8.D0/173.D0*t14074
      t14076 = t13577*t11019
      t14078 = t13490*t11076
      t14080 = t13502*t11987
      t14081 = 9.D0/1384.D0*t14080
      t14082 = t13568*t4040
      t14084 = t12296*t10469
      t14085 = t13580*t14084
      t14087 = t13574*t11226
      t14089 = t13486*t13600
      t14091 = -2.D0/173.D0*t14064+t14066/173-25.D0/2768.D0*t14070-25.D0
     #/5536.D0*t14072+t14075+t14076/346-t14078/173+t14081+479.D0/2768.D0
     #*t14082-9.D0/5536.D0*t14085-13.D0/1384.D0*t14087+17.D0/22144.D0*t1
     #4089
      t14094 = t12296*t9674
      t14095 = t13580*t14094
      t14097 = t13502*t12263
      t14098 = 9.D0/1384.D0*t14097
      t14099 = t13486*t13608
      t14101 = t13554*t7111
      t14103 = t13613*t11886
      t14104 = 8.D0/173.D0*t14103
      t14105 = t13505*t11312
      t14107 = t12296*t9789
      t14108 = t13580*t14107
      t14111 = t4892*t12441*t3452
      t14112 = t14111*t10499
      t14114 = t12296*t11762
      t14115 = t13580*t14114
      t14117 = t12296*t11805
      t14118 = t13580*t14117
      t14120 = t13555*t11012
      t14122 = -9.D0/5536.D0*t14095-t14098+17.D0/11072.D0*t14099+479.D0/
     #5536.D0*t14101-t14104+t14105/346-9.D0/11072.D0*t14108-9.D0/11072.D
     #0*t14112-9.D0/11072.D0*t14115-9.D0/11072.D0*t14118-5.D0/692.D0*t14
     #120
      t14123 = t14111*t10031
      t14125 = t12296*t9658
      t14126 = t13580*t14125
      t14128 = t13516*t12241
      t14130 = t13516*t12220
      t14132 = t13502*t10028
      t14133 = 9.D0/692.D0*t14132
      t14134 = t14044*t14117
      t14136 = t13528*t12241
      t14138 = t12296*t9654
      t14139 = t14044*t14138
      t14141 = t14044*t14114
      t14143 = t13528*t12220
      t14145 = t13502*t10074
      t14146 = 9.D0/2768.D0*t14145
      t14147 = t13502*t10071
      t14148 = 9.D0/1384.D0*t14147
      t14149 = -9.D0/22144.D0*t14123-9.D0/11072.D0*t14126+37.D0/1384.D0*
     #t14128+37.D0/2768.D0*t14130-t14133+9.D0/44288.D0*t14134+2.D0/173.D
     #0*t14136-9.D0/22144.D0*t14139+9.D0/44288.D0*t14141+t14143/173+t141
     #46-t14148
      t14151 = t13580*t14138
      t14153 = t10712*t14035
      t14155 = t13555*t7434
      t14157 = t13502*t10084
      t14158 = 9.D0/2768.D0*t14157
      t14159 = t13966*t10028
      t14161 = t13966*t10071
      t14163 = t14044*t14107
      t14165 = t14044*t14125
      t14167 = t13966*t12263
      t14169 = t10703*t14035
      t14171 = t13966*t9783
      t14173 = t13966*t11987
      t14175 = 9.D0/5536.D0*t14151-37.D0/11072.D0*t14153-13.D0/5536.D0*t
     #14155+t14158+9.D0/2768.D0*t14159+9.D0/5536.D0*t14161+9.D0/44288.D0
     #*t14163+9.D0/44288.D0*t14165+9.D0/5536.D0*t14167-t14169/692-9.D0/5
     #536.D0*t14171-9.D0/5536.D0*t14173
      t14176 = t13554*t7416
      t14178 = t13630*t10857
      t14180 = t6254*t9739
      t14181 = t13613*t14180
      t14183 = t13528*t10018
      t14185 = t12296*t10727
      t14186 = t13486*t14185
      t14188 = t12296*t9643
      t14189 = t14044*t14188
      t14191 = t13524*t13674
      t14193 = t13524*t13671
      t14195 = t12296*t9989
      t14196 = t13486*t14195
      t14198 = t13540*t13655
      t14200 = t13580*t13663
      t14202 = t13584*t13663
      t14204 = 13.D0/2768.D0*t14176-37.D0/22144.D0*t14178-5.D0/173.D0*t1
     #4181+43.D0/2768.D0*t14183+37.D0/88576.D0*t14186+37.D0/22144.D0*t14
     #189+205.D0/22144.D0*t14191+205.D0/44288.D0*t14193+37.D0/44288.D0*t
     #14196+709.D0/44288.D0*t14198-77.D0/22144.D0*t14200+125.D0/11072.D0
     #*t14202
      t14208 = t13547*t13655
      t14210 = t13739*t7416
      t14212 = t12296*t10723
      t14213 = t14044*t14212
      t14215 = t3692*t9719
      t14216 = t13510*t14215
      t14218 = t13493*t9973
      t14219 = 39.D0/692.D0*t14218
      t14220 = t13499*t11140
      t14222 = t13502*t11140
      t14224 = t13505*t9973
      t14226 = t13700*t7111
      t14228 = t3970*t9719
      t14229 = t13510*t14228
      t14231 = t14069*t10499
      t14233 = 39.D0/5536.D0*t14208-143.D0/5536.D0*t14210+37.D0/44288.D0
     #*t14213+423.D0/2768.D0*t14216-t14219-125.D0/1384.D0*t14220+77.D0/2
     #768.D0*t14222-709.D0/5536.D0*t14224-71.D0/692.D0*t14226+423.D0/553
     #6.D0*t14229+9.D0/2768.D0*t14231
      t14234 = t13966*t10084
      t14236 = t13966*t10074
      t14238 = t14069*t10031
      t14240 = t13555*t7340
      t14242 = t14044*t14084
      t14244 = t14044*t14094
      t14246 = t14044*t12328
      t14248 = t13555*t6983
      t14251 = t13584*t12296*t10359
      t14252 = 3.D0/2768.D0*t14251
      t14253 = t13490*t10677
      t14255 = t13577*t10872
      t14257 = t13551*t10828
      t14259 = -9.D0/11072.D0*t14234-9.D0/11072.D0*t14236+9.D0/5536.D0*t
     #14238+5.D0/692.D0*t14240+9.D0/22144.D0*t14242+9.D0/22144.D0*t14244
     #+37.D0/22144.D0*t14246+5.D0/346.D0*t14248+t14252-27.D0/5536.D0*t14
     #253+93.D0/22144.D0*t14255-69.D0/22144.D0*t14257
      t14261 = t13486*t13541
      t14263 = t13540*t13810
      t14265 = t9544*t12292
      t14266 = t14265*t10393
      t14269 = t13540*t12296*t10824
      t14271 = t14265*t10677
      t14274 = t13584*t12296*t10325
      t14275 = 3.D0/1384.D0*t14274
      t14276 = t13490*t10393
      t14278 = t13577*t10429
      t14280 = t13551*t10659
      t14282 = t13540*t13904
      t14284 = t13499*t10349
      t14286 = -3.D0/44288.D0*t14261+3.D0/5536.D0*t14263+9.D0/692.D0*t14
     #266-21.D0/22144.D0*t14269+9.D0/1384.D0*t14271+t14275-27.D0/2768.D0
     #*t14276+93.D0/11072.D0*t14278-69.D0/11072.D0*t14280-21.D0/22144.D0
     #*t14282-3.D0/346.D0*t14284
      t14288 = t13540*t12296*t10434
      t14290 = t13584*t13907
      t14291 = 3.D0/2768.D0*t14290
      t14292 = t6938*t10348
      t14293 = t13528*t14292
      t14294 = 27.D0/1384.D0*t14293
      t14295 = t13521*t14292
      t14296 = 9.D0/346.D0*t14295
      t14297 = t13499*t10875
      t14299 = t13505*t10376
      t14301 = t720*t10348
      t14302 = t13510*t14301
      t14304 = t13709*t7096
      t14307 = t13555*t13728*t1070
      t14310 = t13555*t13761*t1070
      t14312 = t13540*t13544
      t14314 = t13505*t10045
      t14315 = 3.D0/2768.D0*t14314
      t14316 = -21.D0/11072.D0*t14288+t14291+t14294-t14296-3.D0/692.D0*t
     #14297+21.D0/2768.D0*t14299+69.D0/2768.D0*t14302-5.D0/173.D0*t14304
     #+143.D0/22144.D0*t14307-243.D0/44288.D0*t14310+3.D0/11072.D0*t1431
     #2-t14315
      t14319 = t13630*t9720
      t14321 = t13528*t9984
      t14323 = t13584*t14188
      t14325 = t13551*t9840
      t14327 = t13505*t10054
      t14328 = 3.D0/1384.D0*t14327
      t14329 = t13630*t10045
      t14331 = t12296*t9845
      t14332 = t13584*t14331
      t14334 = t13577*t9818
      t14336 = t13547*t14195
      t14338 = t13551*t10015
      t14340 = t13490*t9934
      t14342 = -37.D0/22144.D0*t14319+43.D0/2768.D0*t14321+125.D0/5536.D
     #0*t14323+147.D0/11072.D0*t14325-t14328+3.D0/11072.D0*t14329+13.D0/
     #2768.D0*t14332-93.D0/11072.D0*t14334+39.D0/2768.D0*t14336-423.D0/1
     #1072.D0*t14338-43.D0/5536.D0*t14340
      t14343 = t13540*t14195
      t14345 = t13577*t9758
      t14347 = t13595*t9934
      t14349 = t13580*t14331
      t14351 = t14044*t13703
      t14353 = t13486*t13706
      t14355 = t13589*t9762
      t14357 = t13584*t14212
      t14359 = t13580*t14188
      t14361 = t13547*t14185
      t14363 = t13910*t13872
      t14365 = t13595*t10175
      t14367 = 709.D0/22144.D0*t14343-89.D0/5536.D0*t14345-37.D0/5536.D0
     #*t14347+57.D0/11072.D0*t14349-37.D0/44288.D0*t14351-5.D0/5536.D0*t
     #14353+317.D0/11072.D0*t14355+125.D0/11072.D0*t14357-77.D0/11072.D0
     #*t14359+39.D0/5536.D0*t14361+37.D0/11072.D0*t14363-37.D0/11072.D0*
     #t14365
      t14369 = t13589*t11962
      t14371 = t13580*t14212
      t14373 = t13577*t11933
      t14375 = t13540*t14185
      t14377 = t13551*t10835
      t14379 = t13490*t10175
      t14381 = t13493*t9720
      t14382 = 39.D0/1384.D0*t14381
      t14383 = t13505*t9720
      t14385 = t13540*t13548
      t14387 = t13589*t12169
      t14389 = t13595*t10298
      t14391 = t13486*t13655
      t14393 = 317.D0/22144.D0*t14369-77.D0/22144.D0*t14371-89.D0/11072.
     #D0*t14373+709.D0/44288.D0*t14375-423.D0/22144.D0*t14377-43.D0/1107
     #2.D0*t14379-t14382-709.D0/11072.D0*t14383-51.D0/5536.D0*t14385-37.
     #D0/2768.D0*t14387+37.D0/5536.D0*t14389+37.D0/88576.D0*t14391
      t14394 = t14044*t13663
      t14396 = t13910*t13674
      t14398 = t13910*t13671
      t14400 = t13528*t9966
      t14401 = 13.D0/692.D0*t14400
      t14402 = t13499*t11055
      t14403 = 13.D0/346.D0*t14402
      t14404 = t14044*t13581
      t14406 = t13493*t9836
      t14407 = 13.D0/173.D0*t14406
      t14408 = t6254*t9834
      t14409 = t13613*t14408
      t14411 = t13646*t7111
      t14413 = t13486*t13548
      t14415 = t13524*t13872
      t14417 = t14044*t14331
      t14419 = 37.D0/44288.D0*t14394-317.D0/44288.D0*t14396-317.D0/88576
     #.D0*t14398-t14401-t14403-37.D0/22144.D0*t14404+t14407+5.D0/173.D0*
     #t14409-37.D0/11072.D0*t14411-5.D0/2768.D0*t14413-5.D0/692.D0*t1441
     #5-37.D0/44288.D0*t14417
      wI70pc1 = t13519+t13952+t13924+t13653+t13893+t13625+t14091+t14063+
     #t13592+t14038+t14010+t14393+t14367+t14342+t14286+t14316+t13753+t14
     #259+t13721+t13688+t14233+t14419+t14204+t14175+t13865+t14149+t13837
     #+t13558+t13813+t14122+t13781+t13980
      t14446 = -7.D0/192.D0*t12335-7.D0/384.D0*t12339+5.D0/12.D0*t12343+
     #t12348/48-9.D0/32.D0*t12350+7.D0/192.D0*t12354+t12357/24+3.D0/32.D
     #0*t12360+5.D0/6.D0*t12362-3.D0/16.D0*t12365+9.D0/64.D0*t12368
      t14467 = -t12414/96+3.D0/32.D0*t12418-t12420/8-5.D0/24.D0*t12423+t
     #12426-9.D0/32.D0*t12427+t12430-3.D0/16.D0*t12433+7.D0/192.D0*t1243
     #5-t12443/768+5.D0/12.D0*t12446
      wI71pc1 = -t12294/8+t12298/48-t12306/768-3.D0/8.D0*t12311-7.D0/384
     #.D0*t12314-t12317/384+3.D0/16.D0*t12322+9.D0/64.D0*t12326+9.D0/32.
     #D0*t12329-9.D0/16.D0*t12332+t14446+7.D0/96.D0*t12371+t12376-5.D0/2
     #4.D0*t12385-t12389/48+t12392/192-t12394/96-5.D0/12.D0*t12397-t1240
     #1/4+t12408/384+t12411/384+t14467
      wI72pc1 = wV27
      t14481 = t12461/160+t12469/80+t12473/80-t12480/16-t12483/160+t1248
     #7/8-t12491/160-t12496/160-t12499/16-t12505-t12512/80+t12516+t12518
     #-t12521/40+t12526/40+t12533+t12535-t12539
      t14490 = -t12545-t12548-t12552/4-t12554/8+t12558/40-t12561/40-t125
     #66-t12567/80-t12569/80-t12574-t12578-t12582+t12584/80+t12586/80+t1
     #2591-t12594-t12597-t12601
      t14506 = -t12606/20-t12610/40+t12613/8+t12616/4+t12619/80+t12621/8
     #-t12625-t12629/20+t12633+t12634/80+t12636/40-t12638/40-t12643/40+t
     #12645/80-t12648/40-t12653-t12656-t12659/80
      t14524 = t12662/8+t12665/20+t12667/40+t12669/80-t12671/80+t12674/4
     #0-t12676/160-t12678/40+t12680/40+t12684/40+t12686/80-t12688/8-t126
     #90/80+t12692/40-t12694/80+t12696/40+t12699+t12700/20
      t14532 = t12706+t12708+t12710+t12712-t12715/40-t12718-t12722/80-t1
     #2725+t12726/20+t12729+t12730/40+t12733+t12735-t12737+t12739-t12741
     #/80-t12744+t12746
      t14549 = t12749+t12751-t12752/20+t12754/80-t12759/80-t12761/20-t12
     #764/80-t12766/80+t12769/8+t12771/80+t12774/80+t12776/40+t12778/80+
     #t12780/160+t12782/40+t12785/80+t12788/16+t12790/80
      t14569 = t12797/80-t12801/80-t12804/8-t12807/80+t12809/4+t12811/40
     #-t12813/40-t12815/40-t12818/160-t12821/80+t12825/16+t12829/160-t12
     #832/80-t12834/40-t12839/80+t12842/80-t12846/80-t12848/4
      t14588 = -t12851/8+t12853/40+t12855/40+t12857/160-t12859/40-t12861
     #/40-t12864/8-t12866/40+t12868/160-t12870/160+t12872/80-t12875/40+t
     #12878/20+t12881/40-t12884/80-t12887/40+t12890/80+t12893/160
      wI73pc1 = t14481+t14490+t14506+t14524+t14532+t14549+t14569+t14588
      t14609 = -7.D0/2048.D0*t12900+7.D0/512.D0*t12904-7.D0/256.D0*t1290
     #7+3.D0/64.D0*t12909-3.D0/32.D0*t12913+t12917/8+7.D0/1024.D0*t12920
     #-7.D0/512.D0*t12922-3.D0/64.D0*t12925-3.D0/128.D0*t12927-17.D0/204
     #8.D0*t12929-7.D0/512.D0*t12931+17.D0/1024.D0*t12933+7.D0/1024.D0*t
     #12935-t12938/16+15.D0/256.D0*t12942+t12945/16+t12948/16
      t14622 = -13.D0/1024.D0*t12951+t12954/32-t12956/32+13.D0/512.D0*t1
     #2960+t12966+13.D0/256.D0*t12967-t12972+t12974-7.D0/2048.D0*t12976-
     #15.D0/256.D0*t12978+7.D0/512.D0*t12981-t12984+3.D0/32.D0*t12986-t1
     #2989/16-t12993-t12996-15.D0/128.D0*t12997+13.D0/512.D0*t12999
      t14639 = t13004-15.D0/256.D0*t13005-t13008+3.D0/64.D0*t13009+t1301
     #1/16+7.D0/1024.D0*t13013+17.D0/512.D0*t13015-13.D0/512.D0*t13018+t
     #13020/32-3.D0/64.D0*t13023+t13027-13.D0/1024.D0*t13028+17.D0/1024.
     #D0*t13030-3.D0/128.D0*t13032-3.D0/64.D0*t13035-17.D0/2048.D0*t1303
     #8-15.D0/128.D0*t13041+7.D0/1024.D0*t13043
      t14655 = 15.D0/512.D0*t13046-15.D0/64.D0*t13048+t13051-15.D0/128.D
     #0*t13052+3.D0/32.D0*t13056-t13058/16+3.D0/16.D0*t13061+3.D0/32.D0*
     #t13063+15.D0/256.D0*t13065-t13070-17.D0/1024.D0*t13071-t13073/8+15
     #.D0/512.D0*t13075+15.D0/128.D0*t13077-7.D0/1024.D0*t13079+t13082-t
     #13083/32+15.D0/256.D0*t13086
      wI74pc1 = t14609+t14622+t14639+t14655
      wI75pc1 = t13113/8+t13116/4+t13119/8+t13123/16+t13127/8+t13131/16+
     #t13091/16+t13095/8+t13099/16-t13146/8-t13149/4-t13152/8-t13103/8-t
     #13106/4-t13109/8-t13157-t13160-t13163
      t14694 = 3.D0/16.D0*t13207-17.D0/248.D0*t13210-13.D0/248.D0*t13213
     #-13.D0/248.D0*t13216-21.D0/248.D0*t13219+3.D0/8.D0*t13222+25.D0/19
     #84.D0*t13226-25.D0/992.D0*t13229-25.D0/248.D0*t13232-25.D0/496.D0*
     #t13235+3.D0/16.D0*t13238
      t14718 = 13.D0/124.D0*t13273+21.D0/62.D0*t13276+21.D0/124.D0*t1327
     #9+13.D0/62.D0*t13282-21.D0/248.D0*t13286-25.D0/496.D0*t13289-25.D0
     #/992.D0*t13292-17.D0/496.D0*t13295+25.D0/1984.D0*t13298-17.D0/124.
     #D0*t13301-13.D0/248.D0*t13304
      t14743 = 17.D0/248.D0*t13333+17.D0/124.D0*t13335+17.D0/248.D0*t133
     #39-15.D0/62.D0*t13341+13.D0/62.D0*t13344-13.D0/248.D0*t13346+13.D0
     #/496.D0*t13349+17.D0/248.D0*t13352-81.D0/248.D0*t13355-13.D0/248.D
     #0*t13358-25.D0/496.D0*t13361
      t14767 = 25.D0/1984.D0*t13398-17.D0/496.D0*t13401-3.D0/32.D0*t1340
     #4+81.D0/248.D0*t13407-3.D0/16.D0*t13410+25.D0/992.D0*t13413-3.D0/3
     #2.D0*t13417-13.D0/124.D0*t13420+25.D0/496.D0*t13424+25.D0/992.D0*t
     #13428+25.D0/992.D0*t13432
      wI76pc1 = t14743+t14718-25.D0/992.D0*t13330-81.D0/124.D0*t13324+21
     #.D0/124.D0*t13327-81.D0/248.D0*t13322+13.D0/248.D0*t13318-13.D0/12
     #4.D0*t13320+13.D0/496.D0*t13316+81.D0/496.D0*t13309+15.D0/124.D0*t
     #13311-25.D0/992.D0*t13314+13.D0/124.D0*t13269+13.D0/62.D0*t13266+1
     #3.D0/124.D0*t13260+25.D0/496.D0*t13263-17.D0/124.D0*t13249+25.D0/2
     #48.D0*t13252+13.D0/62.D0*t13255-15.D0/62.D0*t13257-17.D0/62.D0*t13
     #246+25.D0/496.D0*t13243-25.D0/496.D0*t13198+13.D0/124.D0*t13201-15
     #.D0/31.D0*t13203+13.D0/31.D0*t13194-13.D0/124.D0*t13190+25.D0/1984
     #.D0*t13181+17.D0/124.D0*t13184+17.D0/248.D0*t13187+25.D0/992.D0*t1
     #3178+81.D0/496.D0*t13175-25.D0/992.D0*t13394-13.D0/124.D0*t13391-2
     #5.D0/496.D0*t13374+15.D0/62.D0*t13376-13.D0/124.D0*t13380-25.D0/99
     #2.D0*t13383+15.D0/124.D0*t13385-13.D0/248.D0*t13388-13.D0/62.D0*t1
     #3371-21.D0/124.D0*t13367+t14767+t14694
      wI77pc1 = -t13464/16-t13468/8-t13472/16-t13452/8-t13456/4-t13460/8
     #+t13441+t13445+t13449+t13479+t13482+t13485
      t14786 = -7.D0/2768.D0*t13488+13.D0/692.D0*t13491-t13495-9.D0/692.
     #D0*t13497-67.D0/1384.D0*t13500+93.D0/5536.D0*t13503-23.D0/346.D0*t
     #13506+21.D0/2768.D0*t13508+27.D0/2768.D0*t13512+21.D0/2768.D0*t135
     #14+25.D0/2768.D0*t13517
      t14795 = -t13523+9.D0/1384.D0*t13526+t13530+9.D0/2768.D0*t13532+t1
     #3535-t13537-9.D0/1384.D0*t13538+3.D0/5536.D0*t13542-3.D0/11072.D0*
     #t13545-13.D0/346.D0*t13549+55.D0/1384.D0*t13552-249.D0/5536.D0*t13
     #556
      t14808 = -13.D0/2768.D0*t13560-353.D0/5536.D0*t13562-39.D0/5536.D0
     #*t13566+91.D0/2768.D0*t13569-13.D0/2768.D0*t13572+39.D0/2768.D0*t1
     #3575-47.D0/1384.D0*t13578+59.D0/2768.D0*t13582+13.D0/692.D0*t13585
     #-13.D0/5536.D0*t13587-25.D0/2768.D0*t13590
      t14820 = -13.D0/692.D0*t13593+25.D0/5536.D0*t13596+13.D0/1384.D0*t
     #13598-7.D0/2768.D0*t13601+17.D0/5536.D0*t13604+15.D0/692.D0*t13606
     #-7.D0/1384.D0*t13609-7.D0/2768.D0*t13611-t13616-15.D0/1384.D0*t136
     #18+15.D0/346.D0*t13620-6.D0/173.D0*t13623
      t14832 = t13629-17.D0/2768.D0*t13631+15.D0/692.D0*t13633-3.D0/173.
     #D0*t13636-t13639+3.D0/2768.D0*t13640+3.D0/5536.D0*t13642+17.D0/138
     #4.D0*t13644+231.D0/5536.D0*t13647-3.D0/173.D0*t13649+2.D0/173.D0*t
     #13651
      t14845 = -109.D0/2768.D0*t13656+39.D0/2768.D0*t13658+39.D0/1384.D0
     #*t13660-91.D0/5536.D0*t13664-13.D0/1384.D0*t13667+91.D0/1384.D0*t1
     #3669-39.D0/5536.D0*t13672-39.D0/2768.D0*t13675+281.D0/11072.D0*t13
     #679-25.D0/5536.D0*t13681+13.D0/1384.D0*t13684+25.D0/5536.D0*t13686
      t14858 = -25.D0/2768.D0*t13691-181.D0/5536.D0*t13694-25.D0/5536.D0
     #*t13696+91.D0/2768.D0*t13698-t13702-47.D0/2768.D0*t13704+55.D0/276
     #8.D0*t13707+17.D0/692.D0*t13710+13.D0/1384.D0*t13713+25.D0/5536.D0
     #*t13715+17.D0/692.D0*t13717+47.D0/1384.D0*t13719
      t14870 = 47.D0/692.D0*t13722-17.D0/1384.D0*t13726-17.D0/2768.D0*t1
     #3730+167.D0/2768.D0*t13733-13.D0/692.D0*t13735-13.D0/346.D0*t13737
     #-59.D0/2768.D0*t13740-55.D0/1384.D0*t13743-t13746+353.D0/1384.D0*t
     #13747-353.D0/2768.D0*t13749-249.D0/2768.D0*t13751
      t14885 = 17.D0/1384.D0*t13757-131.D0/2768.D0*t13759-231.D0/5536.D0
     #*t13763+231.D0/5536.D0*t13765-13.D0/5536.D0*t13767-249.D0/5536.D0*
     #t13769-353.D0/5536.D0*t13771+231.D0/2768.D0*t13773-167.D0/1384.D0*
     #t13775+353.D0/2768.D0*t13777-151.D0/1384.D0*t13779
      t14897 = 13.D0/1384.D0*t13784+151.D0/2768.D0*t13787+131.D0/2768.D0
     #*t13790-231.D0/11072.D0*t13793-131.D0/1384.D0*t13795+17.D0/692.D0*
     #t13797-151.D0/692.D0*t13799+131.D0/5536.D0*t13802-25.D0/2768.D0*t1
     #3804-17.D0/1384.D0*t13806+t13809-3.D0/5536.D0*t13811
      t14909 = 167.D0/1384.D0*t13815+167.D0/2768.D0*t13817-25.D0/5536.D0
     #*t13819-25.D0/5536.D0*t13821-25.D0/2768.D0*t13823+59.D0/5536.D0*t1
     #3825+t13828-9.D0/5536.D0*t13829-9.D0/5536.D0*t13831-13.D0/692.D0*t
     #13833+13.D0/1384.D0*t13835
      t14921 = t13839+17.D0/1384.D0*t13841-17.D0/2768.D0*t13843+25.D0/11
     #072.D0*t13846+25.D0/2768.D0*t13849-17.D0/692.D0*t13851-17.D0/346.D
     #0*t13853+7.D0/1384.D0*t13855-25.D0/2768.D0*t13857+7.D0/692.D0*t138
     #59+17.D0/692.D0*t13861-17.D0/692.D0*t13863
      t14935 = -17.D0/1384.D0*t13868+151.D0/1384.D0*t13870+13.D0/692.D0*
     #t13873-281.D0/22144.D0*t13877+181.D0/11072.D0*t13879+93.D0/5536.D0
     #*t13881-67.D0/1384.D0*t13883+109.D0/1384.D0*t13885+17.D0/1384.D0*t
     #13887+13.D0/692.D0*t13889-17.D0/2768.D0*t13891
      t14947 = -17.D0/1384.D0*t13894-17.D0/692.D0*t13896-9.D0/1384.D0*t1
     #3899-27.D0/2768.D0*t13902-27.D0/5536.D0*t13905+45.D0/5536.D0*t1390
     #8+25.D0/5536.D0*t13911+3.D0/1384.D0*t13913-17.D0/1384.D0*t13915-t1
     #3918+17.D0/346.D0*t13919+9.D0/692.D0*t13922
      t14959 = 25.D0/5536.D0*t13926+27.D0/2768.D0*t13928-45.D0/2768.D0*t
     #13930+t13933+13.D0/1384.D0*t13934-55.D0/692.D0*t13937-59.D0/1384.D
     #0*t13939+151.D0/2768.D0*t13941+13.D0/1384.D0*t13943-45.D0/1384.D0*
     #t13945+9.D0/2768.D0*t13948-t13951
      t14970 = 13.D0/2768.D0*t13953-55.D0/1384.D0*t13956-59.D0/2768.D0*t
     #13958-t13961-167.D0/692.D0*t13962+13.D0/2768.D0*t13964+25.D0/2768.
     #D0*t13967+13.D0/1384.D0*t13969-t13972-25.D0/1384.D0*t13973+17.D0/3
     #46.D0*t13976+7.D0/1384.D0*t13978
      t14984 = 17.D0/692.D0*t13985+25.D0/5536.D0*t13987-181.D0/2768.D0*t
     #13989-25.D0/2768.D0*t13991+25.D0/1384.D0*t13993-25.D0/2768.D0*t139
     #95-t13998-167.D0/1384.D0*t13999+t14004+13.D0/1384.D0*t14006+3.D0/3
     #46.D0*t14008
      t14996 = -2.D0/173.D0*t14011+249.D0/2768.D0*t14013-t14016-25.D0/55
     #36.D0*t14017-17.D0/346.D0*t14020-25.D0/5536.D0*t14022+7.D0/1384.D0
     #*t14024+4.D0/173.D0*t14027+7.D0/1384.D0*t14029-45.D0/2768.D0*t1403
     #1+353.D0/2768.D0*t14033+17.D0/1384.D0*t14036
      t15008 = 13.D0/692.D0*t14040-17.D0/1384.D0*t14042+25.D0/11072.D0*t
     #14045+17.D0/1384.D0*t14047+47.D0/1384.D0*t14049+25.D0/2768.D0*t140
     #51+13.D0/346.D0*t14053-181.D0/5536.D0*t14055-17.D0/2768.D0*t14057+
     #t14060+3.D0/173.D0*t14061
      t15019 = -4.D0/173.D0*t14064+2.D0/173.D0*t14066-15.D0/692.D0*t1407
     #0-15.D0/1384.D0*t14072+t14075+3.D0/346.D0*t14076-2.D0/173.D0*t1407
     #8+t14081+249.D0/1384.D0*t14082-9.D0/2768.D0*t14085-13.D0/692.D0*t1
     #4087+17.D0/5536.D0*t14089
      t15031 = -9.D0/2768.D0*t14095-t14098+17.D0/2768.D0*t14099+249.D0/2
     #768.D0*t14101-t14104+7.D0/692.D0*t14105-9.D0/5536.D0*t14108-9.D0/2
     #768.D0*t14112-9.D0/5536.D0*t14115-9.D0/5536.D0*t14118-17.D0/1384.D
     #0*t14120
      t15041 = -9.D0/5536.D0*t14123-9.D0/5536.D0*t14126+25.D0/1384.D0*t1
     #4128+25.D0/2768.D0*t14130-t14133+9.D0/11072.D0*t14134+17.D0/692.D0
     #*t14136-9.D0/5536.D0*t14139+9.D0/11072.D0*t14141+17.D0/1384.D0*t14
     #143+t14146-t14148
      t15054 = 9.D0/2768.D0*t14151-25.D0/5536.D0*t14153-13.D0/1384.D0*t1
     #4155+t14158+9.D0/1384.D0*t14159+9.D0/2768.D0*t14161+9.D0/11072.D0*
     #t14163+9.D0/11072.D0*t14165+9.D0/2768.D0*t14167-17.D0/2768.D0*t141
     #69-9.D0/2768.D0*t14171-9.D0/2768.D0*t14173
      t15067 = 13.D0/1384.D0*t14176-25.D0/11072.D0*t14178-17.D0/692.D0*t
     #14181+13.D0/692.D0*t14183+25.D0/22144.D0*t14186+25.D0/5536.D0*t141
     #89+181.D0/5536.D0*t14191+181.D0/11072.D0*t14193+25.D0/11072.D0*t14
     #196+23.D0/692.D0*t14198-93.D0/11072.D0*t14200+67.D0/2768.D0*t14202
      t15081 = 39.D0/2768.D0*t14208-131.D0/2768.D0*t14210+25.D0/11072.D0
     #*t14213+109.D0/692.D0*t14216-t14219-67.D0/692.D0*t14220+93.D0/2768
     #.D0*t14222-23.D0/173.D0*t14224-151.D0/1384.D0*t14226+109.D0/1384.D
     #0*t14229+9.D0/1384.D0*t14231
      t15093 = -9.D0/5536.D0*t14234-9.D0/5536.D0*t14236+9.D0/2768.D0*t14
     #238+17.D0/1384.D0*t14240+9.D0/5536.D0*t14242+9.D0/5536.D0*t14244+2
     #5.D0/5536.D0*t14246+17.D0/692.D0*t14248-t14252-27.D0/2768.D0*t1425
     #3+45.D0/5536.D0*t14255-27.D0/5536.D0*t14257
      t15105 = -3.D0/11072.D0*t14261+3.D0/2768.D0*t14263+9.D0/346.D0*t14
     #266-21.D0/5536.D0*t14269+9.D0/692.D0*t14271-t14275-27.D0/1384.D0*t
     #14276+45.D0/2768.D0*t14278-27.D0/2768.D0*t14280-21.D0/5536.D0*t142
     #82+3.D0/692.D0*t14284
      t15114 = -21.D0/2768.D0*t14288-t14291+t14294-t14296+3.D0/1384.D0*t
     #14297+21.D0/1384.D0*t14299+27.D0/1384.D0*t14302-17.D0/692.D0*t1430
     #4+131.D0/5536.D0*t14307-231.D0/11072.D0*t14310+3.D0/5536.D0*t14312
     #-t14315
      t15127 = -25.D0/11072.D0*t14319+13.D0/692.D0*t14321+67.D0/1384.D0*
     #t14323+55.D0/2768.D0*t14325-t14328+3.D0/5536.D0*t14329+13.D0/1384.
     #D0*t14332-47.D0/2768.D0*t14334+39.D0/1384.D0*t14336-109.D0/1384.D0
     #*t14338-13.D0/692.D0*t14340
      t15140 = 23.D0/346.D0*t14343-91.D0/2768.D0*t14345-25.D0/2768.D0*t1
     #4347+59.D0/5536.D0*t14349-25.D0/11072.D0*t14351-7.D0/2768.D0*t1435
     #3+281.D0/5536.D0*t14355+67.D0/2768.D0*t14357-93.D0/5536.D0*t14359+
     #39.D0/2768.D0*t14361+25.D0/2768.D0*t14363-25.D0/5536.D0*t14365
      t15153 = 281.D0/11072.D0*t14369-93.D0/11072.D0*t14371-91.D0/5536.D
     #0*t14373+23.D0/692.D0*t14375-109.D0/2768.D0*t14377-13.D0/1384.D0*t
     #14379-t14382-23.D0/346.D0*t14383-13.D0/2768.D0*t14385-25.D0/1384.D
     #0*t14387+25.D0/2768.D0*t14389+25.D0/22144.D0*t14391
      t15163 = 25.D0/11072.D0*t14394-281.D0/11072.D0*t14396-281.D0/22144
     #.D0*t14398-t14401-t14403-25.D0/5536.D0*t14404+t14407+17.D0/692.D0*
     #t14409-25.D0/5536.D0*t14411-7.D0/1384.D0*t14413-17.D0/692.D0*t1441
     #5-25.D0/11072.D0*t14417
      wI78pc1 = t14996+t14897+t14870+t14858+t14845+t14832+t14820+t14808+
     #t14959+t15153+t15163+t15140+t15127+t15114+t15105+t15093+t15081+t14
     #970+t15067+t15054+t14947+t15041+t14935+t15031+t15019+t14921+t14795
     #+t14786+t14909+t15008+t14984+t14885
      t15168 = t12437*t7137
      t15175 = t13761*t3409
      t15178 = t12441*t6937
      t15187 = t13728*t3609
      t15190 = t6967*t720
      t15195 = t13761*t3692
      t15198 = t12441*t6953
      t15201 = t7129*t15168/12+t6932*t5*t13782*t9098/4-t7129*t15175/192+
     #t7141*t15178/192+17.D0/48.D0*t6932*t7111-t12359*t10240/4-17.D0/48.
     #D0*t7148*t6968-t7148*t15187/12-t12300*t15190/6-t12309*t12241/2+t71
     #29*t15195/96-t7141*t15198/96
      t15202 = t6932*t6937
      t15207 = t6967*t1093
      t15212 = t6967*t1096
      t15217 = t13761*t3970
      t15220 = t12441*t6949
      t15227 = t6967*t5992
      t15232 = -17.D0/24.D0*t15202*t4040+t12359*t9667/2+17.D0/24.D0*t714
     #8*t15207+t12352*t12356/6+t12300*t15212/12+t12309*t12220/4-t7129*t1
     #5217/192+t7141*t15220/192+17.D0/48.D0*t15202*t7022-t12359*t10200/4
     #-17.D0/48.D0*t7148*t15227-t12352*t12347/12
      wI79pc1 = t15201+t15232
      wI80pc1 = wV28
      t15233 = t12456*t168
      t15292 = t12618*t15233*t6459/20+t12618*t15233*t6456/2+t12465*t5396
     #*t6462/4+t12465*t5299*t12494/40-t12618*t15233*t6481/40+t12465*t687
     #0*t12524/10-t12465*t13415*t12524/20+t12465*t5316*s*t12795/10-t1246
     #5*t6870*t12467/20-t12465*t5299*t842*t12795/20+t12465*t5278*s*t7578
     #*t5851/40-t12618*t15233*t6453/4+t12618*t15233*t12844/20-t12618*t15
     #233*t1979/40-t12618*t15233*t2501/4-t12618*t15233*t6179/40+t12618*t
     #15233*t6462/4+t12618*t15233*t6465/40
      t15296 = t6021*Mzsq
      t15309 = t7396*t12457
      t15331 = t1387*t53
      t15341 = t5869*Mzsq
      t15345 = t6134*Mzsq
      t15353 = -t12508*t7648*t5926/20-t12508*t7648*t15296/20+t12508*t752
     #1*t5926/10+t12508*t7578*t6562/40-t12508*t7465*t6562/20-t15309*t533
     #1*t12627/10+t12508*t7521*t15296/10+t15309*t5387*t12627/20-t15309*t
     #5293*t6047/4+t15309*t5293*t12608/20+t15309*t5293*t5941/2-t15309*t5
     #293*t12604/10-t15309*t15331*t15296/20-t15309*t15331*t5926/20-t1530
     #9*t6325*t12550/2-t15309*t5899*t15341/20+t15309*t5899*t15345/40+t15
     #309*t1387*t59*t6562/40
      t15355 = t5992*t5851
      t15385 = t7396*t12452
      t15386 = t12456*t1387
      t15405 = t15309*t5936*t15355/4+t15309*t6332*t12510/40-t15309*t6332
     #*t12519/20+t15309*t15331*t5940/20+t15309*t6225/20-t15309*t6222/40-
     #t15309*t5387*t12650/40+t15309*t5331*t12650/20+t15309*t6326/20-t153
     #09*t6336/4+t15309*t5293*t12626/20-t15385*t15386*t6468/40-t15309*t6
     #325*t6011/40+t15309*t6333/40-t15385*t15386*t12650/40+t15309*t6329/
     #4+t15309*t6339/40-t15385*t15386*t6113/40
      t15460 = -t12465*t5316*t12837/10+t12618*t15233*t12837/20+t12618*t1
     #5233*t6438/20-t12618*t15233*t12720/40-t12618*t15233*t12757/40+t124
     #65*t5316*t12720/20+t12465*t5396*t12837/20-t12465*t5396*t12720/40+t
     #12465*t6287*t12657/20-t12465*t6287*t12641/10-t12465*t5278*t12524/2
     #0-t12465*t5299*t12795/20-t12465*t5316*t6475/20+t12618*t15233*t6475
     #/40-t12465*t5316*t6465/20+t12465*t5278*t12467/40-t12465*t5316*t646
     #2/2+t12465*t5278*t12471/40
      wI81pc1 = t15292+t15353+t15405+t15460
      t15466 = t6949*t8093
      t15469 = t6953*t8093
      t15491 = t7992*t168
      t15504 = -t7992*t7111/16-t8222*t12964/8-7.D0/512.D0*t12902*t15466+
     #7.D0/256.D0*t12902*t15469-t12988*t8195/4-t8305*t15207/8+5.D0/64.D0
     #*t8305*t15187-7.D0/512.D0*t8305*t13761*t8093+t8305*t15227/16-5.D0/
     #64.D0*t8342*t15168+17.D0/512.D0*t8229*t15220+17.D0/256.D0*t8342*t1
     #5195+t12988*t8248/8+5.D0/32.D0*t15491*t7995-17.D0/512.D0*t8342*t15
     #217+7.D0/512.D0*t8342*t12441*t7994-17.D0/256.D0*t8229*t15198+7.D0/
     #512.D0*t12898*t6967*t1865
      t15511 = t7992*t6937
      t15532 = t7992*t5*t6937
      t15544 = t7992*t13728
      t15549 = -5.D0/32.D0*t12958*t8357+5.D0/32.D0*t12906*t15190-5.D0/64
     #.D0*t12906*t15212+t15511*t4040/8-t15511*t7022/16+17.D0/512.D0*t822
     #9*t15178+5.D0/16.D0*t12958*t8500-7.D0/256.D0*t12898*t6967*t1224+t1
     #2988*t8484/8-17.D0/512.D0*t8342*t15175+t8305*t6968/16-5.D0/32.D0*t
     #12958*t12356-t15532*t3970*t1852/8-5.D0/32.D0*t12958*t6967*t445+t15
     #532*t3692*t1852/4+5.D0/64.D0*t12958*t12347-5.D0/16.D0*t15544*t1306
     #0+5.D0/32.D0*t15544*t13055
      wI82pc1 = t15504+t15549
      t15550 = s*t12452
      t15551 = t8831*t15550
      t15552 = t9433*t1811
      t15556 = t8831*t12501
      t15557 = t9433*t62
      t15558 = t15557*t5895
      t15561 = t15557*t6011
      t15564 = t54*t12452
      t15565 = t8831*t15564
      t15569 = t8831*t5299
      t15570 = t12452*t9433
      t15572 = t15570*t5980*t842
      t15575 = t8831*t5910
      t15576 = t15570*t15341
      t15582 = t15570*t6059
      t15591 = t15570*t5980
      t15594 = t8831*t5930
      t15595 = t15570*t5895
      t15598 = t9433*t1832
      t15599 = t15598*t5869
      t15603 = t15570*t6258*t5851
      t15606 = t8831*t5936
      t15609 = t15551*t15552*t6134/4-3.D0/8.D0*t15556*t15558+3.D0/16.D0*
     #t15556*t15561-t15565*t15557*t5869/4+3.D0/16.D0*t15569*t15572+t1557
     #5*t15576/4-t15551*t9612*t6011/16-3.D0/8.D0*t15569*t15582+t15565*t1
     #5557*t6134/8+t15551*t9612*t5895/8+3.D0/16.D0*t15569*t15591-t15594*
     #t15595/8+t15556*t15599/4+3.D0/4.D0*t15569*t15603+t15606*t15595/8
      t15610 = t15570*t6011
      t15614 = t8831*t5*t12452
      t15618 = t5*t168
      t15619 = t8831*t15618
      t15620 = t1811*t5894
      t15630 = t8831*t6262
      t15631 = t15570*t1093
      t15635 = t8831*t5930*t12452
      t15636 = t15557*t5952
      t15640 = t8831*t5299*t12452
      t15641 = t15552*t5952
      t15644 = t9433*t1844
      t15649 = t8831*t5910*t12452
      t15656 = t15598*t5875
      t15659 = t15598*t6034
      t15663 = t8831*t1070*t12452
      t15664 = t15552*t5875
      t15667 = t15552*t6034
      t15672 = -t15606*t15610/16-t15614*t15552*t6011/16-t15619*t15570*t1
     #5620*CsTh/8+t15594*t15610/16+t15614*t15552*t5895/8-t15630*t15631/8
     #+t15635*t15636/8+t15640*t15641/16-t15551*t15644*t5875/8+t15649*t96
     #12*t5952/16+t15551*t15644*t6034/16-t15565*t15656/4+t15565*t15659/8
     #-t15663*t15664/8+t15663*t15667/16+t15635*t15656/8
      t15675 = t8831*t6262*t12452
      t15679 = t8831*t15570
      t15680 = t1832*t5851
      t15697 = t8831*t15550*t9433
      t15703 = t15598*t6134
      t15707 = t15570*t1093*t5851
      t15714 = t8831*t15564*t9433
      t15727 = t15675*t15557*t5875/8+t15679*t15680*t5875/8-t15635*t15659
     #/16+t15649*t15664/4-t15649*t15667/8-t15675*t15557*t6034/16-t15679*
     #t15680*t6034/16+t15697*t6364/4-t15697*t6344*t6034/8-t15556*t15703/
     #8+3.D0/4.D0*t15594*t15707-t15551*t15552*t5869/2+t15714*t5903*t5875
     #/8-t15714*t5903*t6034/16+t15551*t15644*t5863/16+t15663*t15552*t586
     #3/16
      t15741 = t15557*t5851
      t15748 = t9612*t5851
      t15755 = t15552*t5851
      t15765 = t8831*t12506
      t15778 = t15565*t15598*t5863/8-t15594*t15570*t6214*t5862/16-t15575
     #*t15570*t6094/8-t15630*t15570*t6530/16-t15635*t15741*t5875/4+t1563
     #5*t15741*t6034/8-t15649*t15748*t5875/8+t15649*t15748*t6034/16+t156
     #40*t15755*t6034/16-t15679*t15680*t5863/16-t15640*t15755*t5875/8+t1
     #5765*t15558/8-t15551*t15641/8-t15565*t15636/16+t15619*t15570*t1562
     #0*t842/16-t15765*t15561/16
      t15783 = t8831*t5899
      t15784 = t15570*t15345
      t15787 = t8831*t5274
      t15790 = t15570*t15355
      t15795 = t15570*t12494
      t15803 = t15570*t1844*t842
      t15806 = t15570*t12756
      t15809 = t15570*t6254
      t15812 = t15570*t5992
      t15815 = t15570*t6258
      t15818 = t15570*t15620
      t15821 = t15570*t5903
      t15824 = -t15606*t15707/2+t15783*t15784/8-t15787*t15603/4+t15606*t
     #15790/4-t15783*t15576/4+t15787*t15795/8+t15787*t15582/4-t15787*t15
     #572/8+t15569*t15803/16+3.D0/16.D0*t15594*t15806+3.D0/16.D0*t15575*
     #t15809+t15630*t15812/16-3.D0/8.D0*t15575*t15815+t15619*t15818/16+t
     #15606*t15821/4
      t15826 = t15570*t1844*CsTh
      t15829 = t15570*t6266
      t15832 = t15570*t5874
      t15835 = t15570*t6344
      t15846 = t8831*t6294
      t15863 = -t15569*t15826/8-3.D0/8.D0*t15594*t15829+t15783*t15832/8+
     #t15787*t15835/8-t15787*t15591/8-t15575*t15832/8-3.D0/8.D0*t15594*t
     #15821-3.D0/8.D0*t15569*t15835+t15846*t15631/8-t15606*t15806/16-t15
     #783*t15809/8-3.D0/8.D0*t15569*t15795-t15846*t15812/16+t15783*t1581
     #5/4-t15575*t15784/8+t15606*t15829/8
      t15868 = t15570*t6468
      t15873 = t8831*t12452
      t15876 = t8831*t168
      t15902 = -3.D0/8.D0*t15594*t15790-t8831*t1387*t15868/16-t13300*t15
     #818/16+t15873*t15558/4-t15876*t15570*t15680/8+3.D0/16.D0*t15876*t1
     #5868+t8890*t15835/4+t15873*t15703/8-t15873*t15599/4-t15873*t15561/
     #8-t8890*t15591/16+t9038*t15821/8-t9038*t15806/8+t9038*t15829/4-t88
     #52*t15809/16+t8852*t15815/8
      t15910 = t8831*t1754
      t15920 = t9433*t5894
      t15951 = -t8890*t15803/16+t8890*t15826/8-t8852*t12506*t15552/8-t15
     #910*t12506*t15557/16-t9038*t12506*t15598/16+t15910*t12501*t15557/1
     #6+t9038*t12501*t15920/16+3.D0/16.D0*t8852*t12501*t15552-t9038*t125
     #06*t15920/16+3.D0/16.D0*t9038*t12501*t15598+t8890*t12501*t15644/16
     #-t15873*t15557*t5894/8+t15873*t15598*t5851/8-t9038*t15570*t1832/8-
     #t8852*t15570*t1811/16-t8890*t15570*t1844/16
      wI83pc1 = t15609+t15672+t15727+t15778+t15824+t15863+t15902+t15951
      t15954 = t9682*t7096
      t15956 = t10277*t10947
      t15958 = t13486*t11893
      t15960 = t4892*t13761
      t15962 = t15960*t1224*t11247
      t15964 = t15960*t13955
      t15966 = t13580*t11102
      t15968 = t13540*t9840
      t15970 = t13486*t11674
      t15972 = t5590*t6949
      t15973 = t15972*t11248
      t15975 = t5590*t6953
      t15976 = t15975*t11248
      t15978 = t13580*t11293
      t15980 = t13580*t11019
      t15983 = t15960*t1865*t11247
      t15985 = t13580*t11541
      t15987 = 13.D0/692.D0*t15954-13.D0/692.D0*t15956+17.D0/692.D0*t159
     #58-17.D0/692.D0*t15962-55.D0/2768.D0*t15964+13.D0/173.D0*t15966+13
     #.D0/692.D0*t15968-17.D0/1384.D0*t15970-5.D0/173.D0*t15973+10.D0/17
     #3.D0*t15976-10.D0/173.D0*t15978+5.D0/173.D0*t15980+17.D0/1384.D0*t
     #15983+5.D0/173.D0*t15985
      t15989 = t9682*t6953*t13614
      t15992 = t9682*t6949*t13614
      t15994 = t14044*t10751
      t15996 = t14044*t11980
      t15998 = t14044*t11741
      t16000 = t5598*t7111
      t16002 = t14044*t11786
      t16004 = t14044*t11842
      t16006 = t14044*t10499
      t16008 = t14044*t10031
      t16010 = t5598*t6953
      t16011 = t16010*t12263
      t16013 = t5598*t6949
      t16014 = t16013*t12263
      t16016 = t4892*t13728
      t16017 = t16016*t14215
      t16019 = t16016*t14228
      t16021 = t9682*t13728
      t16022 = t16021*t14019
      t16024 = -8.D0/173.D0*t15989+4.D0/173.D0*t15992+9.D0/1384.D0*t1599
     #4+9.D0/1384.D0*t15996-9.D0/2768.D0*t15998+13.D0/692.D0*t16000-9.D0
     #/2768.D0*t16002+9.D0/1384.D0*t16004-9.D0/692.D0*t16006+9.D0/1384.D
     #0*t16008+9.D0/692.D0*t16011-9.D0/1384.D0*t16014-63.D0/346.D0*t1601
     #7+63.D0/692.D0*t16019+13.D0/346.D0*t16022
      t16026 = t15960*t14215
      t16028 = t16021*t14180
      t16030 = t15960*t14228
      t16032 = t13540*t10015
      t16034 = t13580*t9758
      t16036 = t13540*t10835
      t16038 = t13580*t11933
      t16040 = t5590*t7416
      t16042 = t10277*t6938
      t16044 = t211*t6937
      t16045 = t16044*t3999
      t16047 = t14044*t9758
      t16049 = t13486*t10835
      t16051 = t14044*t11933
      t16053 = t13486*t10015
      t16055 = -13.D0/692.D0*t16026-13.D0/692.D0*t16028+13.D0/1384.D0*t1
     #6030+63.D0/346.D0*t16032-39.D0/692.D0*t16034-63.D0/692.D0*t16036+3
     #9.D0/1384.D0*t16038+13.D0/692.D0*t16040-13.D0/692.D0*t16042-20.D0/
     #173.D0*t16045+13.D0/346.D0*t16047-13.D0/1384.D0*t16049-13.D0/692.D
     #0*t16051+13.D0/692.D0*t16053
      t16056 = t9682*t13782
      t16058 = t16056*t6258*t10348
      t16061 = t16056*t6254*t10348
      t16063 = t13540*t10659
      t16065 = t13540*t10828
      t16067 = t13486*t10686
      t16070 = t15960*t7521*t10044
      t16073 = t15960*t7648*t10044
      t16075 = t16016*t14301
      t16077 = t16016*t13511
      t16079 = t13486*t11638
      t16081 = t10277*t7621
      t16083 = t5590*t7111
      t16085 = t10277*t7340
      t16087 = t16044*t7022
      t16089 = t10277*t11012
      t16091 = 9.D0/173.D0*t16058-9.D0/346.D0*t16061-3.D0/692.D0*t16063+
     #3.D0/1384.D0*t16065+3.D0/2768.D0*t16067+3.D0/1384.D0*t16070-3.D0/2
     #768.D0*t16073+3.D0/692.D0*t16075-3.D0/1384.D0*t16077-3.D0/1384.D0*
     #t16079+20.D0/173.D0*t16081-13.D0/692.D0*t16083-13.D0/692.D0*t16085
     #-20.D0/173.D0*t16087+13.D0/692.D0*t16089
      t16094 = t16044*t9821
      t16096 = t16044*t9734
      t16098 = t5598*t7416
      t16100 = t13677*t11662
      t16102 = t9682*t6937
      t16103 = t16102*t13614
      t16105 = t13677*t11962
      t16107 = t10712*t10175
      t16109 = t13662*t11933
      t16111 = t13840*t14292
      t16113 = t13662*t10429
      t16115 = t10703*t10393
      t16117 = t4892*t168
      t16118 = t16117*t11514
      t16120 = t13662*t10872
      t16122 = t10703*t10677
      t16124 = -13.D0/692.D0*t16094-67.D0/1384.D0*t16096-129.D0/1384.D0*
     #t16098-17.D0/1384.D0*t16100+4.D0/173.D0*t16103-181.D0/2768.D0*t161
     #05+13.D0/692.D0*t16107+67.D0/1384.D0*t16109-9.D0/346.D0*t16111+3.D
     #0/173.D0*t16113-9.D0/173.D0*t16115-3.D0/1384.D0*t16118-3.D0/346.D0
     #*t16120+9.D0/346.D0*t16122
      t16125 = t5590*t6937
      t16126 = t16125*t5549
      t16128 = t16125*t11248
      t16130 = t13845*t10045
      t16132 = t13677*t10122
      t16134 = t13693*t10791
      t16136 = t13845*t11527
      t16138 = t10277*t6937
      t16139 = t16138*t11126
      t16141 = t16138*t11853
      t16143 = t10712*t7612
      t16145 = t10712*t7900
      t16147 = t10712*t7625
      t16149 = t10712*t7903
      t16151 = t16125*t5552
      t16153 = t13845*t15207
      t16155 = t5598*t6937
      t16156 = t16155*t12263
      t16158 = -113.D0/346.D0*t16126-5.D0/173.D0*t16128-3.D0/2768.D0*t16
     #130+3.D0/2768.D0*t16132+3.D0/1384.D0*t16134+17.D0/1384.D0*t16136-4
     #.D0/173.D0*t16139+8.D0/173.D0*t16141+13.D0/692.D0*t16143-13.D0/346
     #.D0*t16145-113.D0/692.D0*t16147+113.D0/346.D0*t16149+113.D0/692.D0
     #*t16151+13.D0/346.D0*t16153-9.D0/1384.D0*t16156
      t16160 = t13848*t12220
      t16162 = t13848*t12241
      t16164 = t16155*t11987
      t16166 = t16155*t9783
      t16168 = t13845*t10968
      t16170 = t13845*t10899
      t16172 = t16155*t10074
      t16174 = t16155*t10084
      t16176 = t13662*t11102
      t16178 = t13654*t11079
      t16180 = t13677*t12169
      t16182 = t13693*t10990
      t16184 = t13845*t15227
      t16186 = t13677*t10012
      t16188 = t10712*t9811
      t16190 = -13.D0/692.D0*t16160+13.D0/346.D0*t16162-9.D0/1384.D0*t16
     #164-9.D0/1384.D0*t16166-9.D0/2768.D0*t16168-9.D0/2768.D0*t16170+9.
     #D0/2768.D0*t16172+9.D0/2768.D0*t16174-13.D0/346.D0*t16176+13.D0/17
     #3.D0*t16178-13.D0/173.D0*t16180+13.D0/692.D0*t16182-13.D0/692.D0*t
     #16184+13.D0/346.D0*t16186-13.D0/692.D0*t16188
      t16191 = t13662*t9818
      t16193 = t13654*t9840
      t16195 = t13724*t10146
      t16197 = t13848*t11551
      t16199 = t16117*t12119
      t16201 = t13724*t10195
      t16203 = t13677*t10990
      t16205 = t13875*t10195
      t16207 = t13875*t10146
      t16209 = t13845*t11657
      t16211 = t16102*t14408
      t16213 = t13693*t10012
      t16215 = t16155*t9877
      t16217 = t16102*t13975
      t16219 = t13693*t12169
      t16221 = 13.D0/692.D0*t16191-13.D0/346.D0*t16193-13.D0/173.D0*t161
     #95+13.D0/692.D0*t16197-13.D0/692.D0*t16199+13.D0/346.D0*t16201+55.
     #D0/2768.D0*t16203-13.D0/346.D0*t16205+13.D0/173.D0*t16207+13.D0/69
     #2.D0*t16209-13.D0/692.D0*t16211-13.D0/346.D0*t16213-13.D0/692.D0*t
     #16215+13.D0/346.D0*t16217+13.D0/173.D0*t16219
      t16225 = t10712*t10298
      t16227 = t16155*t11055
      t16229 = t16125*t9877
      t16231 = t4892*t6937
      t16232 = t16231*t13955
      t16234 = t16125*t11055
      t16236 = t16231*t13936
      t16238 = t16117*t11657
      t16240 = t13693*t11573
      t16242 = t13845*t12119
      t16244 = t13724*t11766
      t16246 = t16117*t11583
      t16248 = t13845*t11583
      t16250 = t13677*t11573
      t16252 = t13724*t11821
      t16254 = 13.D0/346.D0*t16225+13.D0/346.D0*t16227+13.D0/346.D0*t162
     #29+13.D0/346.D0*t16232-13.D0/173.D0*t16234-13.D0/173.D0*t16236-13.
     #D0/346.D0*t16238-63.D0/692.D0*t16240-55.D0/2768.D0*t16242+39.D0/69
     #2.D0*t16244+39.D0/1384.D0*t16246-13.D0/692.D0*t16248-13.D0/1384.D0
     #*t16250-39.D0/1384.D0*t16252
      t16255 = t16125*t10860
      t16257 = t16231*t14215
      t16259 = t16125*t11140
      t16261 = t13875*t11821
      t16263 = t13875*t11766
      t16265 = t16117*t9720
      t16267 = t16231*t14228
      t16269 = t16102*t14019
      t16271 = t13693*t9762
      t16273 = t16155*t11140
      t16275 = t13693*t11962
      t16277 = t16155*t10860
      t16279 = t13848*t9984
      t16281 = t16102*t14180
      t16283 = t13845*t9720
      t16285 = -39.D0/1384.D0*t16255+39.D0/692.D0*t16257+39.D0/692.D0*t1
     #6259+181.D0/2768.D0*t16261-181.D0/1384.D0*t16263+63.D0/692.D0*t162
     #65-39.D0/1384.D0*t16267-13.D0/346.D0*t16269-39.D0/692.D0*t16271-13
     #.D0/346.D0*t16273+39.D0/1384.D0*t16275+13.D0/692.D0*t16277-13.D0/6
     #92.D0*t16279+13.D0/692.D0*t16281+13.D0/1384.D0*t16283
      t16287 = t13654*t10015
      t16289 = t13662*t9758
      t16291 = t16044*t10664
      t16293 = t13654*t10835
      t16295 = t13677*t9762
      t16297 = t10712*t9934
      t16299 = t10277*t6983
      t16301 = t16155*t9649
      t16303 = t16155*t4040
      t16305 = t16102*t7216
      t16307 = t16155*t7022
      t16309 = t16102*t7067
      t16311 = t13845*t6968
      t16313 = t10712*t7340
      t16315 = -39.D0/692.D0*t16287-67.D0/692.D0*t16289+3.D0/346.D0*t162
     #91+39.D0/1384.D0*t16293+181.D0/1384.D0*t16295-13.D0/346.D0*t16297+
     #13.D0/346.D0*t16299+9.D0/2768.D0*t16301-13.D0/346.D0*t16303-13.D0/
     #346.D0*t16305+13.D0/692.D0*t16307+13.D0/692.D0*t16309-13.D0/692.D0
     #*t16311+13.D0/692.D0*t16313
      t16316 = t16138*t11075
      t16319 = t10703*t6938*t10182
      t16321 = t13662*t10352
      t16323 = t16044*CsTh
      t16324 = t16323*t10664
      t16326 = t16044*t842
      t16327 = t16326*t10664
      t16329 = t5590*t10240
      t16331 = t10671*t7674
      t16333 = t10712*t7434
      t16335 = t10712*t7674
      t16337 = t10712*t6983
      t16339 = t13848*t7096
      t16341 = t16155*t9687
      t16343 = t16138*t3056
      t16345 = t16138*t7641
      t16347 = t16155*t4253
      t16349 = -4.D0/173.D0*t16316+9.D0/346.D0*t16319-3.D0/346.D0*t16321
     #-3.D0/173.D0*t16324+3.D0/346.D0*t16327+13.D0/692.D0*t16329+129.D0/
     #1384.D0*t16331+13.D0/346.D0*t16333-13.D0/692.D0*t16335-13.D0/346.D
     #0*t16337-13.D0/692.D0*t16339+9.D0/2768.D0*t16341+13.D0/692.D0*t163
     #43+20.D0/173.D0*t16345-129.D0/1384.D0*t16347
      t16352 = t10712*t7942
      t16354 = t10712*t7621
      t16356 = t10671*t7434
      t16358 = t16125*t4253
      t16360 = t16125*t7022
      t16362 = t16138*t3058
      t16364 = t16138*t3020
      t16366 = t16125*t3994
      t16368 = t16155*t4121
      t16370 = t6967*spf14
      t16371 = t13875*t16370
      t16373 = t13724*t16370
      t16375 = t16323*t3994
      t16377 = t16326*t3994
      t16379 = t16125*t4121
      t16381 = -113.D0/692.D0*t16352+13.D0/692.D0*t16354-129.D0/692.D0*t
     #16356+13.D0/692.D0*t16358-13.D0/692.D0*t16360-40.D0/173.D0*t16362-
     #13.D0/346.D0*t16364+113.D0/692.D0*t16366+129.D0/692.D0*t16368+181.
     #D0/2768.D0*t16371-39.D0/1384.D0*t16373+13.D0/346.D0*t16375-13.D0/6
     #92.D0*t16377-13.D0/346.D0*t16379
      t16382 = t16125*t4040
      t16384 = t16125*t9734
      t16386 = t16231*t9723
      t16388 = t16155*t9734
      t16390 = t16102*t9768
      t16392 = t16138*t11950
      t16394 = t16138*t10973
      t16396 = t13654*t11420
      t16398 = t13662*t11424
      t16400 = t10712*t10947
      t16402 = t13677*t10319
      t16404 = t16323*t9734
      t16406 = t16326*t9734
      t16408 = t16155*t9821
      t16410 = t16125*t9821
      t16412 = 13.D0/346.D0*t16382-39.D0/1384.D0*t16384-39.D0/1384.D0*t1
     #6386+13.D0/692.D0*t16388+13.D0/692.D0*t16390+13.D0/346.D0*t16392-1
     #3.D0/692.D0*t16394+39.D0/1384.D0*t16396+67.D0/1384.D0*t16398+13.D0
     #/692.D0*t16400-181.D0/2768.D0*t16402+67.D0/692.D0*t16404-67.D0/138
     #4.D0*t16406-13.D0/692.D0*t16408+13.D0/346.D0*t16410
      t16414 = t6967*t9794
      t16415 = t13724*t16414
      t16417 = t13875*t16414
      t16419 = t10712*t11012
      t16421 = t13693*t10319
      t16423 = t16326*t9821
      t16425 = t16102*t9806
      t16427 = t13693*t11226
      t16429 = t16138*t10245
      t16431 = t16138*t11163
      t16433 = t16231*t9835
      t16435 = t13662*t11649
      t16437 = t16323*t9821
      t16439 = t5598*t13571
      t16441 = t5590*t13571
      t16443 = t16044*t4040
      t16445 = 13.D0/346.D0*t16415-13.D0/346.D0*t16417-13.D0/692.D0*t164
     #19+39.D0/1384.D0*t16421-13.D0/692.D0*t16423-13.D0/692.D0*t16425-13
     #.D0/346.D0*t16427-13.D0/346.D0*t16429+13.D0/692.D0*t16431+13.D0/34
     #6.D0*t16433+13.D0/692.D0*t16435+13.D0/346.D0*t16437+129.D0/1384.D0
     #*t16439-13.D0/692.D0*t16441+40.D0/173.D0*t16443
      t16446 = t16044*t3994
      t16448 = t14044*t11102
      t16450 = t10277*t7942
      t16452 = t13486*t11079
      t16454 = t14044*t9818
      t16456 = t13486*t9840
      t16458 = t16021*t13975
      t16460 = t13540*t11079
      t16462 = t15960*t13936
      t16464 = t16021*t14408
      t16466 = t13580*t9818
      t16468 = t16016*t13936
      t16470 = t16016*t13955
      t16472 = t13677*t11226
      t16474 = t13654*t10218
      t16476 = -13.D0/692.D0*t16446-13.D0/346.D0*t16448+13.D0/692.D0*t16
     #450-55.D0/1384.D0*t16452+13.D0/692.D0*t16454+55.D0/2768.D0*t16456-
     #13.D0/346.D0*t16458-13.D0/346.D0*t16460+55.D0/1384.D0*t16462+13.D0
     #/692.D0*t16464-13.D0/346.D0*t16466+13.D0/346.D0*t16468-13.D0/692.D
     #0*t16470+13.D0/346.D0*t16472-13.D0/346.D0*t16474
      wI84pc1 = t15987+t16024+t16055+t16091+t16124+t16158+t16190+t16221+
     #t16254+t16285+t16315+t16349+t16381+t16412+t16445+t16476
      t16480 = t13554*t3428
      t16481 = t54*wC26
      t16482 = t7137*t16481
      t16485 = t3406*t3452
      t16486 = t16485*t3428
      t16487 = t5590*t16486
      t16488 = t3638*wC16
      t16489 = t7176*t16488
      t16490 = t16487*t16489
      t16492 = t3638*wC18
      t16493 = t7082*t16492
      t16494 = t16487*t16493
      t16496 = t13554*t6943
      t16497 = t1096*t16492
      t16498 = t16496*t16497
      t16500 = t3452*t3428
      t16501 = t10065*t16500
      t16502 = t3638*wC14
      t16503 = t6953*t16502
      t16504 = t16501*t16503
      t16506 = t10065*t16486
      t16507 = t16506*t16503
      t16509 = s*wC14
      t16510 = t3692*t16509
      t16513 = t10232*t16500
      t16516 = t16500*t6937
      t16517 = t10232*t16516
      t16518 = t16517*t16497
      t16520 = t3999*wC24
      t16521 = t6953*t16520
      t16522 = t16501*t16521
      t16524 = t16506*t16521
      t16526 = -t16480*t16482/64+t16490/40+t16494/40-t16498/40+t16504/10
     #-t16507/10-3.D0/40.D0*t16496*t16510+3.D0/40.D0*t16513*t16503+t1651
     #8/40-t16522/10+t16524/10
      t16527 = t16487*t16503
      t16529 = t16513*t16521
      t16531 = t1077*wC26
      t16532 = t6943*t16531
      t16535 = t5598*t12292
      t16536 = t1070*wC18
      t16537 = t6943*t16536
      t16538 = t16535*t16537
      t16540 = t10277*t16500
      t16541 = t6937*s
      t16543 = t16541*t62*wC24
      t16544 = t16540*t16543
      t16546 = t1077*wC18
      t16547 = t6943*t16546
      t16548 = t13595*t16547
      t16550 = t1077*wC24
      t16551 = t6943*t16550
      t16552 = t13595*t16551
      t16555 = t7941*Mzsq*wC24
      t16556 = t16540*t16555
      t16558 = s*wC26
      t16559 = t1093*t16558
      t16562 = t720*t16481
      t16565 = t6987*t16558
      t16568 = -t16527/80-t16529/10-t13595*t16532/1280+t16538/640+t16544
     #/160+t16548/320-t16552/1280-t16556/128-t16496*t16559/40+t16496*t16
     #562/32+t16480*t16565/80
      t16570 = t5590*t16500
      t16571 = t4120*wC24
      t16572 = t6953*t16571
      t16573 = t16570*t16572
      t16575 = t4120*wC18
      t16576 = t6953*t16575
      t16577 = t16570*t16576
      t16579 = t12292*t3428
      t16580 = t10065*t16579
      t16581 = t1072*wC24
      t16582 = t6949*t16581
      t16583 = t16580*t16582
      t16585 = t1070*wC16
      t16586 = t6949*t16585
      t16587 = t16570*t16586
      t16589 = t5598*t16500
      t16590 = t6949*t16546
      t16591 = t16589*t16590
      t16593 = t6949*t16536
      t16594 = t16570*t16593
      t16596 = t54*wC14
      t16597 = t6949*t16596
      t16598 = t16570*t16597
      t16600 = t1070*wC26
      t16601 = t6949*t16600
      t16604 = t5598*t16486
      t16605 = t6949*t16575
      t16606 = t16604*t16605
      t16608 = t16501*t16597
      t16612 = t13595*t6943*t1071*wC24
      t16614 = t16573/80-t16577/20-t16583/160+t16587/320-t16591/320-t165
     #94/160+t16598/640+t16570*t16601/640+t16606/80+t16608/80-t16612/640
      t16619 = t838*wC14
      t16620 = t6943*t16619
      t16621 = t13595*t16620
      t16625 = t13595*t6943*t1071*wC16
      t16627 = t1070*wC24
      t16628 = t6953*t16627
      t16629 = t16540*t16628
      t16631 = t6953*t16536
      t16632 = t16589*t16631
      t16635 = t6943*t1071*wC18
      t16636 = t13595*t16635
      t16638 = t6953*t16585
      t16639 = t16540*t16638
      t16641 = t6953*t16600
      t16644 = t6953*t16596
      t16647 = t16540*t16631
      t16649 = t13490*t16620
      t16651 = t5598*t16579
      t16652 = t16651*t16590
      t16654 = -t13595*t6943*t1071*wC26/640-t16621/640-t16625/320-t16629
     #/320-t16632/160+t16636/160-t16639/160-t16540*t16641/320+3.D0/320.D
     #0*t16540*t16644+t16647/160-t16649/80-t16652/640
      t16659 = t13595*t6943*t1119*wC24
      t16665 = t843*wC14
      t16666 = t6943*t16665
      t16667 = t13595*t16666
      t16671 = t13595*t6943*t1119*wC16
      t16673 = t3638*wC24
      t16674 = t7082*t16673
      t16675 = t16487*t16674
      t16677 = t7137*t16492
      t16678 = t16487*t16677
      t16681 = t6943*t1119*wC18
      t16682 = t13595*t16681
      t16685 = t3609*s*wC18
      t16686 = t16496*t16685
      t16688 = t16513*t16677
      t16690 = t16570*t16628
      t16692 = t16570*t16677
      t16694 = t16659/1280+t13595*t6943*t1119*wC26/1280+t16667/1280+t166
     #71/640-t16675/160+t16678/40-t16682/320-t16686/40+t16688/40-t16690/
     #320-t16692/40
      t16695 = t1072*wC26
      t16696 = t6953*t16695
      t16699 = t1077*wC16
      t16700 = t6953*t16699
      t16701 = t16540*t16700
      t16703 = t6953*t16531
      t16706 = t1072*wC16
      t16707 = t6953*t16706
      t16708 = t16540*t16707
      t16710 = t5590*t16579
      t16711 = t16710*t16700
      t16713 = t6953*t16546
      t16714 = t16710*t16713
      t16716 = t6953*t16550
      t16717 = t16710*t16716
      t16721 = t6953*t16581
      t16722 = t16540*t16721
      t16724 = t16540*t16716
      t16726 = t16589*t16576
      t16728 = -t16540*t16696/80+t16701/32+t16540*t16703/64-t16708/40+t1
     #6711/320-t16714/160+t16717/640+t16710*t16703/640-t16722/80+t16724/
     #64+t16726/40
      t16730 = t16580*t16707
      t16732 = t16580*t16700
      t16734 = t16540*t16713
      t16736 = t16580*t16716
      t16742 = t16570*t16631
      t16744 = t16580*t16721
      t16746 = t16570*t16638
      t16748 = t16589*t16713
      t16750 = t16570*t16644
      t16752 = t16730/40-t16732/40-t16734/160-t16736/80-t16580*t16703/80
     #+t16580*t16696/80+t16742/80+t16744/80-t16746/160+t16748/160-t16750
     #/320
      t16755 = t16604*t16576
      t16757 = t16501*t16644
      t16759 = t6949*t16695
      t16762 = t6949*t16699
      t16763 = t16540*t16762
      t16765 = t6949*t16531
      t16768 = t16501*t16762
      t16770 = t6949*t16571
      t16771 = t16487*t16770
      t16773 = t16487*t16605
      t16775 = t6949*t16706
      t16776 = t16540*t16775
      t16778 = t16710*t16762
      t16780 = t16710*t16590
      t16782 = -t16570*t16641/320-t16755/40-t16757/40+t16540*t16759/160-
     #t16763/64-t16540*t16765/128+t16768/40+t16771/160-t16773/40+t16776/
     #80-t16778/640+t16780/320
      t16786 = t6949*t16550
      t16787 = t16710*t16786
      t16789 = t7082*t16488
      t16790 = t16506*t16789
      t16792 = t1096*t16488
      t16794 = t16496*t16792/8
      t16797 = t10065*t16516
      t16798 = t16797*t16792
      t16801 = t16517*t16792/8
      t16802 = t16487*t16789
      t16804 = t16540*t16786
      t16806 = t16589*t16605
      t16808 = t16540*t16590
      t16810 = t16580*t16786
      t16812 = -t16787/1280+t16790/10+t16794-t16710*t16765/1280-t16798/1
     #0-t16801-t16802/80-t16804/128-t16806/80+t16808/320+t16810/160
      t16813 = t16580*t16775
      t16815 = t16580*t16762
      t16821 = t16540*t16586
      t16825 = t6949*t16627
      t16826 = t16540*t16825
      t16828 = t16589*t16593
      t16832 = t16540*t16593
      t16834 = t13490*t16666
      t16836 = -t16813/80+t16815/80+t16580*t16765/160-t16580*t16759/160+
     #t16821/320+t16540*t16601/640+t16826/640+t16828/320-3.D0/640.D0*t16
     #540*t16597-t16832/320+t16834/160
      t16838 = t16570*t16555
      t16840 = Mzsq*wC18
      t16841 = t7941*t16840
      t16842 = t16570*t16841
      t16846 = t16501*t16543
      t16848 = s*wC16
      t16849 = t6987*t16848
      t16850 = t16501*t16849
      t16852 = t54*wC24
      t16853 = t7103*t16852
      t16854 = t16513*t16853
      t16856 = t54*wC16
      t16857 = t7103*t16856
      t16858 = t16513*t16857
      t16860 = t7103*t16481
      t16863 = Mzsq*wC16
      t16864 = t7941*t16863
      t16865 = t16570*t16864
      t16867 = t54*wC18
      t16868 = t7103*t16867
      t16869 = t16589*t16868
      t16871 = t5590*t16485
      t16872 = t4120*wC16
      t16874 = t16871*t6943*t16872
      t16876 = -t16838/640+t16842/160-t16570*t16482/640-t16846/80-t16850
     #/40+t16854/80+t16858/40+t16513*t16860/80-t16865/320-t16869/80+t168
     #74/80
      t16877 = t16501*t16555
      t16883 = t16501*t16864
      t16886 = t16871*t6943*t16571
      t16888 = t6943*t16575
      t16889 = t16871*t16888
      t16891 = t16513*t16868
      t16893 = t4120*wC26
      t16897 = t3409*wC26
      t16898 = t6987*t16897
      t16901 = t16604*t16677
      t16905 = t10065*t16485
      t16906 = t3999*wC26
      t16910 = t16877/80-t16501*t16565/80+t16501*t16482/80+t16883/40+t16
     #886/160-t16889/40-t16891/40+t16871*t6943*t16893/160+t16513*t16898/
     #20-t16901/80+t16501*t16898/20-t16905*t6943*t16906/20
      t16913 = t7103*t16509
      t16914 = t16501*t16913
      t16916 = t6943*t16502
      t16917 = t16905*t16916
      t16921 = t3409*wC24
      t16922 = t6987*t16921
      t16923 = t16513*t16922
      t16925 = t16871*t16916
      t16927 = t16501*t16922
      t16929 = t3970*t16481
      t16932 = s*wC24
      t16933 = t5992*t16932
      t16934 = t16496*t16933
      t16936 = t5992*t16848
      t16937 = t16496*t16936
      t16939 = t6949*t16872
      t16940 = t16570*t16939
      t16942 = t1096*t16852
      t16943 = t16496*t16942
      t16945 = -t16914/20+t16917/20-3.D0/80.D0*t16513*t16913+t16923/20+t
     #16925/160+t16927/20-t16496*t16929/80+t16934/80+t16937/40-t16940/80
     #-t16943/64
      t16946 = t3970*t16867
      t16947 = t16496*t16946
      t16949 = t16513*t16605
      t16951 = t6949*t16893
      t16954 = t6949*t16906
      t16961 = t16604*t16493
      t16963 = t6949*t16502
      t16964 = t16501*t16963
      t16966 = t16506*t16963
      t16968 = t3970*t16509
      t16971 = t7941*wC14
      t16974 = t16947/40-t16949/40+t16487*t16951/160+t16513*t16954/20+t1
     #6501*t16954/20-t16506*t16954/20-t16961/80-t16964/20+t16966/20+3.D0
     #/80.D0*t16496*t16968-3.D0/320.D0*t16480*t16971
      t16978 = t16487*t16963
      t16980 = t6949*t16520
      t16981 = t16513*t16980
      t16983 = t16501*t16980
      t16985 = t7498*wC18
      t16986 = t16480*t16985
      t16988 = t3638*wC26
      t16989 = t7082*t16988
      t16992 = t1096*t16988
      t16999 = t16506*t16674
      t17001 = t16480*t16555
      t17003 = -3.D0/80.D0*t16513*t16963+t16978/160+t16981/20+t16983/20-
     #t16986/160+t16506*t16989/20+t16496*t16992/16-t16517*t16992/16-t164
     #87*t16989/160+t16999/20-t17001/64
      t17004 = t16480*t16849
      t17006 = t16480*t16543
      t17010 = t16480*t16857
      t17012 = t16480*t16853
      t17018 = t16480*t16631
      t17020 = t16480*t16638
      t17022 = t16480*t16628
      t17024 = t6943*t16627
      t17025 = t13555*t17024
      t17027 = t6943*t16600
      t17030 = t17004/40+t17006/80-t16480*t16860/80-t17010/40-t17012/80-
     #t16480*t16641/160+3.D0/160.D0*t16480*t16644+t17018/80-t17020/80-t1
     #7022/160+t17025/640+t13555*t17027/640
      t17035 = t6943*t16585
      t17036 = t13555*t17035
      t17038 = t13555*t16537
      t17040 = t6943*t16596
      t17045 = t16480*t16868
      t17047 = t16480*t16864
      t17049 = t16480*t16841
      t17051 = t1096*t16673
      t17053 = t16517*t17051/16
      t17054 = t13646*t16537
      t17060 = t17036/320-t17038/320-3.D0/640.D0*t13555*t17040+3.D0/80.D
     #0*t16480*t16913+t17045/40-t17047/32+t17049/160-t17053+t17054/320+t
     #16480*t16601/320-3.D0/320.D0*t16480*t16597
      t17061 = t16480*t16593
      t17063 = t16480*t16586
      t17065 = t16480*t16825
      t17067 = t13700*t17040
      t17071 = t13739*t17040
      t17073 = t16797*t17051
      t17076 = t16496*t17051/16
      t17077 = t13739*t16537
      t17079 = t13739*t17035
      t17081 = t13739*t17024
      t17083 = -t17061/160+t17063/160+t17065/320+t17067/80+t13739*t17027
     #/640+t17071/640-t17073/20+t17076-t17077/160+t17079/320+t17081/640
      t17085 = t3409*wC16
      t17086 = t6987*t17085
      t17087 = t16480*t17086
      t17089 = t16480*t16922
      t17093 = t6953*t16893
      t17099 = t5598*t16516
      t17100 = t720*t16492
      t17101 = t17099*t17100
      t17103 = t16540*t16582
      t17105 = t16570*t16825
      t17107 = t5992*t16558
      t17110 = t1096*t16481
      t17113 = t16570*t16770
      t17115 = -t17087/10-t17089/20-t16480*t16898/20+t16570*t17093/80+t1
     #6496*t1093*t16897/10-t17101/40+t17103/160+t17105/640+t16496*t17107
     #/80-t16496*t17110/64-t17113/160
      t17116 = t16570*t16605
      t17118 = t16570*t16503
      t17121 = t16496*t1093*t16921
      t17124 = t16496*t1093*t17085
      t17132 = t16496*t5992*t16921
      t17134 = t17099*t16497
      t17136 = t16570*t16963
      t17139 = t16496*t5992*t17085
      t17141 = t16651*t16713
      t17144 = t16905*t6943*t16520
      t17146 = t17116/40+t17118/80+t17121/10+t17124/5-t16570*t16951/160-
     #t16496*t5992*t16897/20-t17132/20+t17134/80-t17136/160-t17139/10+t1
     #7141/320-t17144/20
      t17149 = t16513*t17086
      t17151 = t16501*t17086
      t17153 = t3999*wC16
      t17155 = t16905*t6943*t17153
      t17157 = t16535*t16635
      t17159 = t16535*t16681
      t17161 = t7498*wC26
      t17164 = t7137*t16988
      t17167 = t7137*t16673
      t17168 = t16570*t17167
      t17170 = t16570*t16786
      t17172 = t13595*t17040
      t17174 = t7137*t16488
      t17175 = t16570*t17174
      t17177 = t17149/10+t17151/10-t17155/10-t17157/320+t17159/640+t1648
     #0*t17161/320+t16570*t17164/160+t17168/160-t17170/640+t17172/1280+t
     #17175/80
      t17178 = t16570*t16716
      t17180 = t16570*t16590
      t17184 = t843*t16840
      t17185 = t16496*t17184
      t17189 = t843*t16863
      t17190 = t16496*t17189
      t17192 = t16570*t16762
      t17194 = t16513*t16770
      t17196 = t16513*t16939
      t17198 = t16501*t16582
      t17200 = t16501*t16775
      t17202 = t17178/320+t17180/160-t16570*t16765/640+t17185/160+t16513
     #*t16951/80-t17190/32-t17192/320+t17194/80+t17196/40-t17198/80-t172
     #00/40
      t17206 = t16570*t16713
      t17210 = t16487*t16939
      t17212 = t16501*t16786
      t17214 = t16570*t16700
      t17216 = t838*t16840
      t17217 = t16496*t17216
      t17219 = t838*t16863
      t17220 = t16496*t17219
      t17222 = t16513*t16572
      t17224 = t6953*t16872
      t17225 = t16513*t17224
      t17227 = t16501*t16721
      t17229 = t13595*t17027/1280-t17206/80+t16570*t16703/320+t17210/80+
     #t17212/80+t17214/160-t17217/80+t17220/16-t17222/40-t17225/20+t1722
     #7/40
      t17230 = t16501*t16707
      t17234 = t16506*t16980
      t17236 = t6949*t17153
      t17237 = t16513*t17236
      t17239 = t16501*t17236
      t17241 = t16506*t17236
      t17243 = t3970*t16852
      t17244 = t16496*t17243
      t17246 = t3970*t16856
      t17247 = t16496*t17246
      t17249 = t16570*t16857
      t17251 = t16570*t16853
      t17253 = t16570*t16868
      t17257 = t17230/20-t16513*t17093/40-t17234/20+t17237/10+t17239/10-
     #t17241/10-t17244/80-t17247/40-t17249/80-t17251/160+t17253/40-t1657
     #0*t16860/160
      t17261 = t16589*t16677
      t17263 = t16570*t16913
      t17266 = t13490*t6943*t16706
      t17268 = t6943*t16699
      t17269 = t13490*t17268
      t17271 = t16540*t16841
      t17273 = t16589*t16841
      t17275 = t13490*t16551
      t17279 = t7498*wC16
      t17280 = t16480*t17279
      t17286 = t13490*t6943*t16581
      t17288 = t17261/80-t17263/160-t17266/80+t17269/80+t17271/320-t1727
     #3/320+t17275/160+t13490*t16532/160+t17280/160-t13490*t6943*t16695/
     #160-t17286/160
      t17289 = t16487*t17224
      t17291 = t16501*t16716
      t17294 = t5598*t16485*t16888
      t17296 = t13490*t17040
      t17298 = t3609*t16558
      t17301 = t16501*t16700
      t17303 = t16487*t16572
      t17313 = -t17289/40-t17291/40+t17294/80+t17296/160+t16496*t17298/1
     #6-t17301/20-t17303/80+t16501*t16696/40-t16501*t16703/40-t16501*t17
     #164/20+t16506*t17164/20
      t17315 = t16487*t16576
      t17317 = t3692*t16867
      t17318 = t16496*t17317
      t17322 = t16501*t17167
      t17326 = t16506*t17167
      t17328 = t3609*t16932
      t17330 = t16496*t17328/16
      t17332 = t16513*t17167/16
      t17333 = t16487*t17167
      t17335 = t16513*t16576
      t17339 = t17315/20-t17318/20-t16487*t17164/160-t17322/20-t16513*t1
     #7164/16+t17326/20+t17330-t17332-t17333/160+t17335/20-t16487*t17093
     #/80
      t17340 = t16501*t17174
      t17342 = t16506*t17174
      t17344 = t3609*t16848
      t17346 = t16496*t17344/8
      t17348 = t16513*t17174/8
      t17349 = t7498*wC24
      t17350 = t16480*t17349
      t17352 = t16496*t17100
      t17354 = t16517*t17100
      t17356 = t16487*t17174
      t17358 = t7176*t16492
      t17359 = t16487*t17358
      t17361 = t5590*t16516
      t17362 = t720*t16673
      t17363 = t17361*t17362
      t17365 = t720*t16488
      t17366 = t17361*t17365
      t17368 = t17361*t17100
      t17370 = -t17340/10+t17342/10+t17346-t17348+t17350/320+t17352/20-t
     #17354/20-t17356/80-t17359/20-t17363/80-t17366/40+t17368/20
      t17373 = t720*t16988
      t17378 = t7176*t16988
      t17383 = t16797*t17362
      t17385 = t7176*t16673
      t17386 = t16506*t17385
      t17393 = t16496*t17362/8
      t17395 = t16517*t17362/8
      t17396 = t16487*t17385
      t17398 = -t17361*t17373/80+t16797*t17373/10-t16506*t17378/10+t1648
     #7*t17378/80+t17383/10-t17386/10-t16496*t17373/8+t16517*t17373/8-t1
     #7393+t17395+t17396/80
      t17399 = t16797*t17365
      t17405 = t13595*t17024
      t17407 = t16506*t16489
      t17410 = t16496*t17365/4
      t17411 = t6953*t16906
      t17418 = t16604*t17358
      t17421 = t16517*t17365/4
      t17422 = t17399/5-t16501*t16759/80+t16501*t16765/80+t17405/1280-t1
     #7407/5-t17410-t16513*t17411/10-t16501*t17411/10+t16506*t17411/10+t
     #17418/40+t17421
      t17424 = t6953*t17153
      t17425 = t16513*t17424
      t17427 = t16501*t17424
      t17429 = t16535*t16547
      t17431 = t17361*t17051
      t17433 = t17361*t16792
      t17435 = t16506*t17424
      t17437 = t3692*t16852
      t17438 = t16496*t17437
      t17440 = t3692*t16856
      t17441 = t16496*t17440
      t17443 = t17361*t16497
      t17449 = -t17425/5-t17427/5-t17429/640+t17431/160+t17433/80+t17435
     #/5+t17438/40+t17441/20-t17443/40+t17361*t16992/160-t16797*t16992/2
     #0
      t17450 = t13595*t17035
      t17452 = t13595*t16537
      t17454 = t3692*t16481
      t17457 = t1093*t16932
      t17458 = t16496*t17457
      t17462 = t16540*t16864
      t17466 = t16540*t16849
      t17468 = t13595*t17268
      t17470 = t1093*t16848
      t17471 = t16496*t17470
      t17473 = t16570*t17224
      t17475 = t720*t16852
      t17476 = t16496*t17475
      t17478 = t17450/640-t17452/320+t16496*t17454/40-t17458/40+t16540*t
     #16565/160-t17462/64-t16540*t16482/128+t17466/80-t17468/640-t17471/
     #20+t17473/40+t17476/32
      wI85pc1 = t17257+t17478+t17398+t16614+t16945+t17449+t17339+t16910+
     #t16876+t17146+t16836+t17422+t17313+t17370+t17288+t16812+t16782+t17
     #115+t17060+t16752+t17030+t16728+t17003+t16974+t17083+t16568+t16694
     #+t17229+t16526+t17202+t16654+t17177
      t17483 = t3999*wC28
      t17484 = t16044*t17483
      t17485 = t17484/8
      t17487 = t16138*t1077*wC28
      t17488 = t17487/32
      t17491 = t16044*t17153
      t17493 = t5590*t6967
      t17494 = t17493*t17219
      t17496 = t17493*t17216
      t17500 = t17493*t17189
      t17502 = t17493*t17184
      t17504 = t10712*t16841
      t17507 = t3609*s*wC28
      t17508 = t15975*t17507
      t17509 = t17508/16
      t17510 = t15972*t17507
      t17511 = t17510/32
      t17512 = Mzsq*wC28
      t17513 = t838*t17512
      t17514 = t17493*t17513
      t17516 = t17485+t17488-t16138*t16695/288+t17491/36+t17494/288-t174
     #96/288+t17493*t16562/576+t17500/576-t17502/576-t17504/576-t17509-t
     #17511+t17514/64
      t17519 = t5598*t6967
      t17520 = t17519*t17184
      t17522 = t10712*t16864
      t17526 = t15975*t17344
      t17528 = t17519*t17216
      t17530 = t16010*t16685
      t17532 = t15975*t16685
      t17534 = t843*t17512
      t17535 = t17493*t17534
      t17537 = t16013*t16685
      t17539 = t15972*t16685
      t17543 = t15972*t17344
      t17545 = t10277*t16971
      t17547 = t17493*t17110/1152+t17520/576+t17522/576-t15975*t17298/14
     #4-t17526/72+t17528/288-t17530/72+t17532/72+t17535/128-t17537/144+t
     #17539/144-t15972*t17298/288-t17543/144+t17545/288
      t17549 = t16138*t16619
      t17551 = t16138*t16706
      t17554 = t10712*t7941*t17512
      t17556 = t10671*t16841
      t17559 = t16138*t1072*wC28
      t17560 = t17559/32
      t17563 = t16125*t16685
      t17565 = t16138*t17219
      t17567 = t1070*wC28
      t17569 = t10712*t6949*t17567
      t17572 = t10712*t6953*t17567
      t17574 = t10671*t16593
      t17576 = t10671*t16631
      t17578 = t10712*t16593
      t17580 = t17549/144-t17551/144+t17554/128+t17556/576-t17560-t16138
     #*t16559/144+t17563/144+t17565/72-t17569/128-t17572/64-t17574/576-t
     #17576/288+t17578/576
      t17581 = t10712*t16586
      t17583 = t10712*t16597
      t17587 = t10712*t16631
      t17589 = t10712*t16638
      t17591 = t10712*t16644
      t17599 = t16138*t17470
      t17601 = t16138*t17189
      t17605 = t16155*t17317
      t17607 = t16323*t17507
      t17608 = t17607/4
      t17609 = -t17581/576-t17583/1152-t10712*t16601/1152+t17587/288-t17
     #589/288-t17591/576-t10712*t16641/576+t10712*t16482/1152+t16138*t16
     #562/144-t17599/72+t17601/144-t16138*t17107/288+t17605/72-t17608
      t17612 = t16155*t16946
      t17614 = t16138*t16936
      t17618 = t16138*t17513
      t17619 = t17618/16
      t17620 = t16326*t17507
      t17621 = t17620/8
      t17622 = t16125*t16946
      t17624 = t16125*t17246
      t17626 = t62*wC28
      t17628 = t16138*t6046*t17626
      t17629 = t17628/32
      t17630 = t16138*t17534
      t17631 = t17630/32
      t17633 = t16138*t129*t17626
      t17634 = t17633/16
      t17635 = t16125*t17317
      t17637 = t16125*t17440
      t17639 = t16125*t17344
      t17641 = t17612/144-t17614/144+t16138*t17110/288+t17619-t17621-t17
     #622/144+t17624/144-t17629+t17631-t17634-t17635/72+t17637/72-t17639
     #/144
      t17644 = t16125*t16968
      t17650 = t16155*t16685
      t17652 = t16125*t16510
      t17656 = t16323*t17344
      t17660 = t54*wC28
      t17662 = t16125*t3970*t17660
      t17663 = t17662/32
      t17665 = t16125*t3692*t17660
      t17666 = t17665/16
      t17667 = t16125*t17507
      t17668 = t17667/32
      t17669 = t16326*t17344
      t17671 = t16326*t16502
      t17673 = -t16125*t17298/288+t17644/288+t16125*t16929/288+t16125*t1
     #7454/144-t17650/144+t17652/144-t16326*t17298/72-t17656/18-t16323*t
     #17298/36+t17663+t17666-t17668-t17669/36-t17671/72
      t17675 = t16044*Mzsq
      t17677 = t17675*t3638*wC28
      t17678 = t17677/8
      t17681 = t17675*t16488
      t17683 = t16323*t17483
      t17684 = t17683/4
      t17685 = t16326*t17483
      t17686 = t17685/8
      t17689 = t10712*t16971
      t17691 = t10712*t17279
      t17693 = t10712*t16985
      t17696 = t10712*t7498*wC28
      t17698 = t10671*t16985
      t17700 = t16138*t16665
      t17703 = t16125*t4120*wC28
      t17704 = t17703/32
      t17705 = -t17678-t17675*t16988/72-t17681/36+t17684+t17686-t10712*t
     #17161/1152-t17689/1152-t17691/576+t17693/576-t17696/128-t17698/576
     #+t17700/288+t17704
      t17706 = t16155*t16575
      t17708 = t16125*t16872
      t17710 = t16125*t16575
      t17714 = t16125*t16502
      t17718 = t16326*t17153
      t17720 = t16323*t16502
      t17724 = t16323*t17153
      t17728 = t16138*t16699
      t17732 = t16044*t16502
      t17734 = t17706/144+t17708/144-t17710/144+t16125*t16893/288+t17714
     #/288+t16326*t16906/72+t17718/36-t17720/36+t16323*t16906/36+t17724/
     #18+t16138*t16531/288+t17728/144+t16044*t16906/72-t17732/72
      wI86pc1 = t17516+t17547+t17580+t17609+t17641+t17673+t17705+t17734
      t17749 = t17485+t17487/16+t17491/16+t17494/32-t17496/32+t17500/64-
     #t17502/64-t17504/64-t17508/8-t17510/16+t17514/16+t17520/64+t17522/
     #64
      t17754 = t16044*t16520
      t17765 = -t17526/16+t17528/32-t17530/16+t17532/16+t17754/32+t17535
     #/32-t17537/32+t17539/32-t17543/32+t17545/64+t17549/32-t17551/32+t1
     #7554/32+t17556/64
      t17780 = -t17559/16+t17563/32+t17565/16-t17569/32-t17572/16-t17574
     #/64-t17576/32+t17578/64-t17581/64-t17583/128+t17587/32-t17589/32-t
     #17591/64
      t17781 = t15972*t17328
      t17783 = t17493*t17475
      t17785 = t17493*t16942
      t17787 = t15975*t17328
      t17797 = -t17781/64+t17783/64+t17785/128-t17787/32-t17599/16+t1760
     #1/32+t17605/16-t17608+t17612/32-t17614/32+t17618/8-t17621-t17622/3
     #2+t17624/32
      t17800 = t16125*t16571
      t17802 = t10712*t17349
      t17808 = t16326*t16520
      t17810 = t16323*t16520
      t17812 = t17675*t16673
      t17814 = t16138*t16581
      t17816 = t16138*t16550
      t17820 = t17800/64-t17802/128-t17628/16+t17630/16-t17633/8-t17635/
     #16+t17808/32+t17810/16-t17812/32-t17814/64+t17816/64+t17637/16-t17
     #639/32
      t17832 = t17644/64-t17650/32+t17652/32-t17656/8+t17662/16+t17665/8
     #-t17667/16-t17669/16-t17671/32-t17678-t17681/16+t17684+t17686-t176
     #89/128
      t17847 = -t17691/64+t17693/64-t17696/32-t17698/64+t17700/64+t17703
     #/16+t17706/32+t17708/32-t17710/32+t17714/64+t17718/16-t17720/16+t1
     #7724/8
      t17848 = t16125*t17243
      t17850 = t16138*t16942
      t17852 = t16138*t16933
      t17854 = t16138*t17475
      t17856 = t16138*t17457
      t17858 = t16125*t17437
      t17860 = t16125*t17328
      t17862 = t16326*t17328
      t17864 = t16323*t17328
      t17866 = t10712*t16628
      t17868 = t10712*t16555
      t17870 = t10712*t16825
      t17874 = t17848/64+t17850/64-t17852/64+t17854/32-t17856/32+t17858/
     #32-t17860/64-t17862/32-t17864/16-t17866/64+t17868/128-t17870/128+t
     #17728/32-t17732/32
      wI87pc1 = t17749+t17765+t17780+t17797+t17820+t17832+t17847+t17874
      t17877 = t17487/8
      t17884 = t17508/4
      t17885 = t17510/8
      t17887 = t17514/4
      t17894 = t17535/8
      t17898 = -t17887+t17520/52+t17522/52+t17526/26-t17528/26+t17530/26
     #-t17532/26+t17894-t17537/52+t17539/52-t17543/52
      t17903 = t17554/8
      t17905 = t17559/8
      t17908 = t17569/8
      t17909 = t17572/4
      t17922 = -t17574/52+t17576/26+t17578/52-t17581/52-t17583/104-t1758
     #7/26+t17589/26+t17591/52+t17599/26+t17601/52-t17605/26
      t17927 = t17618/4
      t17930 = t17628/8
      t17931 = t17630/8
      t17932 = t17633/4
      t17941 = t17662/8
      t17942 = t17665/4
      t17943 = t17667/8
      t17945 = t17635/26-t17637/26-t17639/52+t17644/104-t17650/52-t17652
     #/52+t17656/26+t17941-t17942-t17943-t17669/52
      t17952 = t17696/8
      t17956 = t17703/8
      t17966 = t17700/104+t17956+t17706/52+t17708/52-t17710/52+t17714/10
     #4+t17718/52+t17720/52-t17724/26+t17728/52-t17732/104
      wI88pc1 = -t17504/52+t17500/52-t17502/52+t17931+t17932+t17496/26+t
     #17491/52-t17494/26+t17485-t17927-t17930+t17909+t17922-t17905-t1790
     #8+t17903+t17898-t17885+t17877+t17884+t17686-t17678-t17681/52-t1768
     #4-t17671/104-t17622/52+t17624/52-t17621-t17614/52+t17612/52+t17608
     #+t17563/52-t17565/26+t17556/52-t17549/52-t17551/52-t17698/52+t1769
     #3/52-t17689/104-t17691/52+t17545/104+t17966-t17952+t17945
      t17969 = t3428*s
      t17970 = t3454*t17969
      t17973 = t3615*s
      t17976 = t3615*t3409
      t17981 = t6932*t9743
      t17982 = t3615*t129
      t17989 = t6932*t9887
      wI91pc1 = t7133*t17970/16-t7133*t17973/4+t7133*t17976/4-t7133*t447
     #6/4-t17981*t17982/16+t7133*t17982/4-t7133*t4487/4+t17989*t4487/4
      t18004 = -t3716/704+t3723/44-t3712/22+t3720/11-t3729/22-t3707/176-
     #t3696/352+t3662/88+t3710/176+t3690/44-t3679/22+t3660/1408
      t18017 = t3669/88-t3671/44+t3673/176-t3636/176-t3694/44+t3700/22+t
     #3641/352-t3653/352+t3727/22-t3611/22-t3621/11+t3627/11
      wI92pc1 = t18004+t18017
      t18018 = t211*wV12
      t18019 = t18018*t3577
      t18020 = t17969*t7
      t18023 = t18018*t9743
      t18024 = t3615*t176
      t18027 = t18018*t3420
      t18032 = t17969*t3692
      t18035 = t18018*t9887
      t18038 = t3428*Mzsq
      t18039 = t3409*cw2
      t18040 = t18038*t18039
      t18043 = t3609*cw2
      t18044 = t3615*t18043
      t18047 = t3428*CsTh
      t18048 = t18047*t176
      t18054 = t3615*t843
      t18059 = t17969*t3970
      t18065 = t18047*t18043
      t18068 = t9887*t3424
      wI93pc1 = -t18019*t18020/4+t18023*t18024/16+t18027*t3631/8-t18023*
     #t3631/32-t18019*t18032/8+t18035*t3640/8+t18019*t18040/4-t18035*t18
     #044/4-t18019*t18048/4+t18018*t3578*t18048/16+t18027*t18054/8-t1802
     #3*t18054/32-t18019*t18059/8+t18035*t3615*t3737/8+t18019*t18065/4-t
     #18018*t18068*t18065/4
      t18084 = t3742/704+t3745/704+t3749/352+t3758/352-t3761/352+t3769/4
     #4-t3775/22-t3782/176-t3785/176-t3806/44-t3809/88+t3813/44
      t18097 = -t3819/44-t3823/88+t3826/22-t3828/22+t3831/44+t3834/44+t3
     #837/44-t3840/22+t3846/1408-t3852/704-t3855/704+t3857/704
      t18111 = t3859/704-t3862/176-t3871/352+t3873/88+t3875/352-t3878/35
     #2+t3888/44-t3898/22+t3901/22+t3914/1408-t3917/2816-t3920/2816
      t18124 = -t3922/176-t3926/352+t3929/88+t3941/352+t3944/44+t3947/44
     #-t3954/88-t3957/88+t3966/88+t3972/88-t3976/44-t3978/44
      wI94pc1 = t18084+t18097+t18111+t18124
      t18126 = t211*wV16
      t18127 = t18126*t3420
      t18130 = t18126*t9743
      t18141 = t18126*t3577
      t18142 = t3428*t54
      t18143 = t18142*t7
      t18146 = t17969*t178
      t18149 = t18126*t9887
      t18150 = t178*t3409
      t18154 = t18142*t1113
      t18157 = t18126*t3578
      t18160 = t18142*t1324
      t18163 = t17969*t1094
      t18172 = t18126*t18068
      t18173 = t3428*t3409
      t18174 = t18173*t176
      t18177 = t3428*t62
      t18178 = t177*t3409
      t18180 = t18177*t18178*CsTh
      t18183 = t18126*t4465
      t18184 = t3638*t1324
      t18187 = t18127*t4166/64-t18130*t4085/256+t18130*t3615*t1078/64-t1
     #8130*t4066/256-3.D0/64.D0*t18130*t3615*t1073-t18141*t18143/16+3.D0
     #/16.D0*t18141*t18146+3.D0/16.D0*t18149*t3615*t18150-t18141*t18154/
     #8+t18157*t18154/32+t18157*t18160/64-t18157*t18163/64-t18141*t18160
     #/16+t18141*t18163/16-t18130*t4073/256-t18172*t18174/16-t18141*t181
     #80/16+t18183*t18184/16
      t18188 = t18173*s
      t18189 = t18188*t1324
      t18198 = t3638*t1113
      t18201 = t18188*t1113
      t18204 = t18177*t18178
      t18207 = t18142*t3692
      t18212 = t18142*t3970
      t18217 = t18142*t4274
      t18232 = -t18172*t18189/16+t18172*t18180/16+t18149*t4141/64+t18141
     #*t18174/16+t18183*t18198/8-t18172*t18201/8-3.D0/16.D0*t18141*t1820
     #4-t18141*t18207/64+t18149*t4122/64-t18141*t18212/64+t18149*t4257/6
     #4-t18141*t18217/64+t18149*t4190/64+t18127*t4066/64+t18127*t4085/64
     #-t18127*t4141/64+t18127*t4073/64-t18127*t3990/256
      wI95pc1 = t18187+t18232
      t18246 = t3986/5632-t3996/1408-3.D0/176.D0*t4001+t4006/176-t4009/1
     #76+t4018/704-t4022/704+t4027/1408+t4031/176-t4034/1408+t4038/352-t
     #4042/176-t4053/88
      t18261 = t4064/88-t4067/704-t4074/704+t4078/88+t4081/704-t4083/88-
     #t4086/704+t4089/176-3.D0/176.D0*t4092-t4103/22528-t4107/22528+t411
     #6/2816-3.D0/352.D0*t4118+t4126/2816
      t18276 = t4135/352-t4137/1408-t4139/1408+t4142/704+t4148/704-t4150
     #/704+t4153/88-t4158/176+t4161/88-t4167/704-t4171/1408+t4174/2816+t
     #4188/352
      t18291 = -t4191/352+t4199/88-t4201/88-t4207/176+t4215/176+t4220/56
     #32+t4226/176-3.D0/1408.D0*t4228+t4235/704-t4239/22528+t4241/44-t42
     #43/44-t4255/704-t4266/88
      t18307 = t4276/704+3.D0/176.D0*t4278-t4280/704-t4285/44+t4287/2816
     #-t4292/704-t4294/88+t4296/5632-t4298/5632+3.D0/88.D0*t4304-t4309/7
     #04+t4313/1408-t4317/704
      t18322 = -t4319/1408-t4323/5632+t4325/5632-t4327/5632+t4329/5632-t
     #4336/352-t4338/704-t4340/88+t4344/176+t4346/352-t4359/352+t4361/35
     #2+t4363/352+t4367/352
      t18337 = t4369/352-t4373/352+t4379/352+t4382/352-t4384/352-t4386/3
     #52+t4390/2816+t4392/5632+t4394/5632+t4396/5632-t4398/5632+t4400/56
     #32+t4402/5632
      t18352 = -t4406/1408+3.D0/1408.D0*t4412+t4418/1408+3.D0/704.D0*t44
     #23+t4427/88-t4429/22528-t4432/1408+3.D0/1408.D0*t4434-t4436/5632-3
     #.D0/5632.D0*t4438-t4442/1408-3.D0/88.D0*t4444+t4446/1408+3.D0/88.D
     #0*t4448
      wI96pc1 = t18246+t18261+t18276+t18291+t18307+t18322+t18337+t18352
      t18367 = t4539/1408+t4625/352-t4623/352-t4535/44+t4634/44+t4629/22
     #-t4547/22-t4583/352-t4532/22+t4639/44-t4621/176+t4585/88
      t18380 = -t4559/1408-t4575/352+t4543/352+t4618/44-t4607/44-t4613/2
     #2+t4616/22+t4572/352+t4562/22-t4627/44+t4595/176-t4568/88
      wI97pc1 = t18367+t18380
      t18393 = t4643/176+3.D0/352.D0*t4653+t4665/704-t4669/352+3.D0/1408
     #.D0*t4673-3.D0/1408.D0*t4676+3.D0/1408.D0*t4681-3.D0/1408.D0*t4684
     #+t4687/1408-3.D0/1408.D0*t4689-t4696/1408+t4698/1408
      t18406 = -t4702/5632+3.D0/5632.D0*t4706-3.D0/352.D0*t4713+3.D0/88.
     #D0*t4717+t4726/352+t4728/1408-t4732/1408+3.D0/704.D0*t4743+3.D0/17
     #6.D0*t4747-3.D0/176.D0*t4752-3.D0/88.D0*t4758+3.D0/88.D0*t4760
      t18420 = t4763/5632+t4767/176+3.D0/176.D0*t4779-3.D0/176.D0*t4782-
     #3.D0/88.D0*t4793+3.D0/1408.D0*t4799-t4803/88+t4805/88+3.D0/88.D0*t
     #4816-3.D0/88.D0*t4818+t4821/88+t4834/176
      t18433 = -t4839/176+t4847/88-t4849/88-3.D0/5632.D0*t4853-t4858/704
     #-t4861/1408-3.D0/176.D0*t4863-t4865/176+3.D0/176.D0*t4867-t4870/17
     #6-t4876/88-3.D0/704.D0*t4883
      wI98pc1 = t18393+t18406+t18420+t18433
      t18435 = t4888/32
      t18436 = t4890/8
      t18437 = t9430*t3577
      t18439 = t18437*t18020/2
      t18442 = t9430*t9743*t18024/8
      t18443 = t4909/32
      t18444 = t4911/8
      t18445 = t4894/8
      t18446 = t4898/8
      t18447 = t4913/8
      t18448 = t4915/8
      t18450 = t18437*t18040/2
      t18453 = t9430*t9887*t18044/2
      wI99pc1 = t18435-t18436+t18439-t18442+t18443-t18444+t18445-t18446+
     #t18447-t18448-t18450+t18453
      t18454 = t4928/2816
      t18455 = t4931/176
      t18456 = t4935/2816
      t18457 = t4938/44
      t18458 = t4945/176
      t18460 = t4962/44
      t18462 = t4974/44
      t18464 = t4987/44
      t18465 = t4996/44
      t18467 = t5009/44
      t18468 = t5012/44
      t18469 = t5022/704
      t18470 = t5027/176
      t18471 = t5031/704
      t18472 = -t18454-t18455+t18456-t18457+t18458+7.D0/176.D0*t4957+t18
     #460+7.D0/176.D0*t4969+t18462-7.D0/176.D0*t4981-t18464+t18465+7.D0/
     #176.D0*t5002-t18467+t18468-t18469+t18470-t18471
      t18473 = t5034/704
      t18474 = t5037/704
      t18475 = t5040/704
      t18478 = t5055/176
      t18479 = t5058/704
      t18481 = t5068/44
      t18484 = t5077/2816
      t18485 = t5079/44
      t18486 = t5087/22
      t18489 = t5100/352
      t18490 = t5102/88
      t18491 = t18473-t18474+t18475+7.D0/704.D0*t5044-7.D0/704.D0*t5049-
     #t18478+t18479-7.D0/176.D0*t5064+t18481+7.D0/176.D0*t5070-7.D0/176.
     #D0*t5072-t18484-t18485-t18486+7.D0/88.D0*t5091-7.D0/352.D0*t5096-t
     #18489+t18490
      t18495 = t5112/22
      t18496 = t5118/1408
      t18497 = t5122/22
      t18498 = t5127/352
      t18499 = t5129/352
      t18501 = t5135/704
      t18503 = t5139/704
      t18506 = t5149/44
      t18507 = t5151/352
      t18508 = t5153/1408
      t18509 = t5155/2816
      t18510 = t5157/88
      t18511 = 7.D0/88.D0*t5104-7.D0/88.D0*t5108-t18495+t18496+t18497+t1
     #8498-t18499+7.D0/704.D0*t5132-t18501-7.D0/176.D0*t5137+t18503+7.D0
     #/176.D0*t5142-7.D0/176.D0*t5145-t18506+t18507-t18508+t18509-t18510
      t18516 = t5176/352
      t18517 = t5178/352
      t18518 = t5180/22
      t18520 = t5193/22
      t18521 = t5198/44
      t18522 = t5200/704
      t18523 = t5203/22
      t18524 = t5206/44
      t18526 = t5217/704
      t18528 = t5224/704
      t18529 = t5232/704
      t18530 = 7.D0/176.D0*t5164-7.D0/176.D0*t5167-7.D0/88.D0*t5169+7.D0
     #/352.D0*t5173-t18516+t18517+t18518-7.D0/88.D0*t5188+t18520-t18521-
     #t18522-t18523+t18524+7.D0/88.D0*t5213-t18526-7.D0/704.D0*t5221+t18
     #528+t18529
      wI100pc1 = t18472+t18491+t18511+t18530
      wI101pc1 = -t5238/36-t5247/9-t5266/18-t5262/18-t5249/72-t5254/9-t5
     #260/9-t5256/18-t5258/9-t5240/36-2.D0/9.D0*t5243-2.D0/9.D0*t5245
      t18556 = t5271/72+t5286/18+t5306/36+t5304/36+t5289/144+t5291/18+t5
     #301/18+t5295/36+t5297/18+t5276/72+t5280/9+t5282/9
      t18569 = -t5311/72-t5323/18-t5344/36-t5341/36-t5326/144-t5329/18-t
     #5338/18-t5333/36-t5335/18-t5314/72-t5318/9-t5320/9
      wI102pc1 = t18556+t18569
      t18582 = t5349/288+t5363/144+t5357/36+t5360/72+t5353/36+t5367/72+t
     #5370/96+t5381/48+t5375/12+t5378/24+t5373/12+t5384/24
      t18595 = t5389/96+t5401/48+t5394/12+t5398/24+t5392/12+t5404/24+t54
     #08/288+t5420/144+t5413/36+t5417/72+t5411/36+t5423/72
      wI103pc1 = t18582+t18595
      wI104pc1 = t5426/36+t5430/72+t5432/9+t5434/9+t5437/18+t5439/18+t54
     #41/36+t5444/72+t5446/9+t5448/9+t5450/18+t5454/18
      t18621 = -t5460/288-t5463/12-t5466/48+t5469/288+t5472/72-t5475/24-
     #t5478/1152-t5480/288+t5483/288-t5485/144-t5487/144+t5489/144+t5491
     #/36
      t18636 = -t5494/576+t5496/36+t5498/144+t5500/144+t5502/144-t5504/2
     #4+t5506/576+t5508/144+t5510/36+t5512/36+t5514/144-t5516/96+t5518/7
     #2+t5520/288
      t18651 = t5524/1152-t5526/1152+t5528/288+t5530/1152+t5533/72-t5536
     #/36-t5538/288+t5540/288-t5542/576+t5544/576+t5547/72-t5550/18+t555
     #3/36
      t18666 = -t5556/144+t5559/72-t5561/36-t5564/72-t5566/144+t5568/144
     #+t5570/288+t5572/288+t5574/72-t5576/18+t5578/36-t5580/144-t5582/12
     #-t5584/288
      wI105pc1 = t18621+t18636+t18651+t18666
      t18686 = t5593/18+t5596/36+t5600/144+t5602/18+t5604/72+t5606/72+t5
     #609/36+t5611/144+t5613/36+t5615/72+t5617/18+t5620/36+t5623/144+t56
     #27/144+t5629/36+t5631/36+t5633/9+t5635/18
      t18705 = t5639/18+t5641/18+t5643/9+t5645/36+t5647/9+t5650/18+t5652
     #/36+t5654/18+t5657/72+t5659/72+t5661/18+t5663/36+t5665/18+t5667/18
     #+t5669/36+t5671/9+t5673/72+t5675/18
      wI106pc1 = t18686+t18705
      wI107pc1 = -wI43pc1
      wI108pc1 = -wI44pc1
      wI109pc1 = -wI45pc1
      wI110pc1 = -wI46pc1
      wI111pc1 = -wI47pc1
      wI112pc1 = -wI48pc1
      t18706 = t6940/32
      t18707 = 23.D0/48.D0*t6945
      t18708 = t6951/1536
      t18709 = t6955/64
      t18710 = t6958/768
      t18711 = t6964/1536
      t18712 = 29.D0/384.D0*t6970
      t18713 = t6973/32
      t18714 = t6977/3
      t18715 = 7.D0/32.D0*t6981
      t18716 = t6990/128
      t18717 = 23.D0/384.D0*t6993
      t18718 = 5.D0/192.D0*t6998
      t18719 = 23.D0/48.D0*t7002
      t18720 = 7.D0/32.D0*t7004
      t18721 = 5.D0/192.D0*t7006
      t18722 = -t18706-t18707+t18708+t18709-t18710+t18711-t18712+t18713+
     #t18714+t18715-t18716+t18717-t18718+t18719-t18720+t18721
      t18723 = t7009/16
      t18724 = t7012/3
      t18725 = 29.D0/48.D0*t7014
      t18726 = 29.D0/48.D0*t7017
      t18727 = t7020/16
      t18728 = 23.D0/96.D0*t7024
      t18729 = 5.D0/384.D0*t7027
      t18730 = 7.D0/64.D0*t7029
      t18731 = 23.D0/96.D0*t7032
      t18732 = 7.D0/64.D0*t7035
      t18733 = 5.D0/384.D0*t7038
      t18734 = t7040/6
      t18735 = 29.D0/96.D0*t7042
      t18736 = 29.D0/96.D0*t7044
      t18737 = t7046/32
      t18738 = t7048/6
      t18739 = t7051/32
      t18740 = -t18723-t18724+t18725-t18726+t18727+t18728-t18729+t18730-
     #t18731-t18732+t18733-t18734+t18735-t18736+t18737+t18738-t18739
      t18742 = t7057/384
      t18743 = t7060/32
      t18744 = t7069/8
      t18745 = t7072/8
      t18746 = t7074/32
      t18747 = t7077/8
      t18748 = t7080/32
      t18749 = t7084/6
      t18750 = t7093/384
      t18751 = 23.D0/96.D0*t7101
      t18752 = 5.D0/384.D0*t7106
      t18753 = t7109/384
      t18754 = 23.D0/96.D0*t7113
      t18755 = 7.D0/64.D0*t7116
      t18756 = 7.D0/64.D0*t7119
      t18757 = 5.D0/384.D0*t7121
      t18758 = -t18742-t18743+t18744-t18745+t18746+t18747+t18748+t18749+
     #t18750+t18751-t18752+t18753-t18754+t18755-t18756+t18757
      t18759 = t7123/6
      t18760 = t7127/6
      t18761 = t7131/32
      t18762 = t7135/8
      t18763 = 7.D0/64.D0*t7139
      t18764 = t7143/128
      t18765 = t7150/32
      t18766 = t7153/384
      t18767 = t7156/1536
      t18768 = t7159/32
      t18769 = t7162/16
      t18770 = t7166/8
      t18771 = t7168/32
      t18772 = t7171/8
      t18773 = t7174/6
      t18774 = 7.D0/128.D0*t7178
      t18775 = t7181/192
      t18776 = t18759-t18760-t18761-t18762-t18763-t18764+t18765-t18766+t
     #18767+t18768-t18769+t18770-t18771+t18772+t18773-t18774+t18775
      t18779 = 23.D0/192.D0*t7185
      t18780 = t7192/12
      t18781 = t7195/32
      t18782 = 29.D0/96.D0*t7197
      t18783 = 29.D0/96.D0*t7199
      t18784 = 5.D0/192.D0*t7202
      t18785 = 23.D0/48.D0*t7205
      t18786 = 7.D0/32.D0*t7208
      t18787 = 29.D0/48.D0*t7211
      t18788 = t7214/4
      t18789 = t7222/4
      t18790 = t7225/4
      t18791 = t7229/192
      t18792 = t7231/16
      t18793 = t7234/16
      t18794 = t7236/4
      t18795 = -t18779+t18780+t18781-t18782+t18783-t18784+t18785+t18786-
     #t18787+t18788+t18789-t18790+t18791+t18792+t18793+t18794
      t18796 = t7238/16
      t18797 = t7240/4
      t18798 = t7243/4
      t18799 = t7246/3
      t18800 = t7249/16
      t18801 = 5.D0/384.D0*t7252
      t18802 = t7256/768
      t18803 = t7259/192
      t18804 = 7.D0/64.D0*t7266
      t18805 = 23.D0/96.D0*t7269
      t18806 = t7271/16
      t18807 = 29.D0/96.D0*t7274
      t18808 = t7276/8
      t18809 = t7279/768
      t18810 = t7283/8
      t18811 = 29.D0/192.D0*t7286
      t18812 = t7289/1536
      t18813 = -t18796-t18797-t18798-t18799-t18800+t18801-t18802-t18803-
     #t18804-t18805-t18806+t18807-t18808+t18809+t18810+t18811+t18812
      t18815 = t7293/8
      t18816 = t7297/1536
      t18817 = t7301/768
      t18818 = t7304/6144
      t18819 = t7307/128
      t18820 = t7310/1536
      t18821 = 5.D0/384.D0*t7313
      t18822 = 23.D0/96.D0*t7316
      t18823 = 29.D0/96.D0*t7319
      t18824 = t7324/8
      t18825 = t7326/8
      t18826 = t7329/384
      t18827 = t7332/32
      t18828 = 7.D0/256.D0*t7335
      t18829 = t7338/8
      t18830 = 23.D0/384.D0*t7344
      t18831 = -t18815+t18816-t18817-t18818-t18819-t18820+t18821-t18822+
     #t18823+t18824-t18825-t18826+t18827+t18828-t18829+t18830
      t18832 = t7348/384
      t18833 = t7350/24
      t18834 = t7352/32
      t18835 = t7354/32
      t18836 = t7356/1536
      t18837 = 29.D0/384.D0*t7358
      t18838 = t7361/6144
      t18839 = t7364/1536
      t18840 = 7.D0/256.D0*t7367
      t18841 = 5.D0/1536.D0*t7371
      t18842 = t7374/128
      t18843 = t7377/32
      t18844 = 5.D0/768.D0*t7379
      t18845 = 5.D0/1536.D0*t7381
      t18846 = t7383/64
      t18847 = t7386/3072
      t18848 = t7391/24
      t18849 = -t18832-t18833+t18834-t18835+t18836-t18837-t18838-t18839+
     #t18840-t18841-t18842+t18843+t18844-t18845+t18846+t18847-t18848
      wI113pc1 = t18722+t18740+t18758+t18776+t18795+t18813+t18831+t18849
      t18870 = -t7400/320+t7405/640+t7414/640+t7418/640-t7421/320+t7429/
     #64+t7432/640+t7442/128-t7448/80-t7450/2560+t7452/1280-t7456/80+t74
     #71/80+t7483/160+t7487/160+t7490/320+t7495/320+t7500/640
      t18889 = t7502/80-t7505/32-t7507/160+t7510/640-t7516/160+t7519/32-
     #t7527/40+t7538/160+t7541/80+t7545/160-t7548/80+t7556/160-t7561/160
     #+t7565/80+t7570/160-t7573/160-t7581/160-t7590/320
      t18909 = -t7593/320+t7596/160+t7599/64-t7607/160+t7610/64-t7613/16
     #0+t7616/640-t7619/2560-t7626/16-t7629/64+t7632/320+t7646/80-t7654/
     #160-t7657/320-t7659/320-t7667/320+t7670/320+t7672/1280
      t18928 = -t7676/2560-t7678/40+t7681/1280+t7695/40+t7708/640+t7710/
     #80-t7713/40+t7715/40+t7717/8-t7719/8-t7721/40-t7724/640+t7740/640-
     #t7742/80-t7744/32+t7746/16+t7748/80-t7752/16
      t18949 = t7754/160+t7756/80+t7758/80-t7760/80+t7763/1280-t7765/320
     #-t7771/64-t7778/320+t7790/80-t7793/80-t7795/16+t7800/16-t7803/80+t
     #7805/80-t7807/320-t7815/320+t7817/160-t7819/16
      t18968 = t7821/160+t7824/160-t7826/160+t7830/80-t7832/640-t7834/32
     #0-t7836/640-t7838/160-t7845/160-t7853/320-t7855/2560+t7857/64-t786
     #2/320+t7864/64+t7867/640+t7872/640+t7876/640-t7878/320
      t18988 = -t7881/320+t7885/640+t7890/160-t7894/32-t7896/320-t7898/6
     #40+t7901/80+t7904/8+t7908/320+t7910/64-t7912/640+t7914/1280+t7916/
     #640+t7918/640-t7920/160+t7926/640-t7928/640-t7930/640
      t19007 = t7933/80-t7935/320-t7937/320-t7939/160+t7946/320+t7950/32
     #0-t7952/256-t7958/320+t7960/1280-t7962/640-t7964/320-t7966/2560-t7
     #968/256+t7971/640+t7975/1280-t7977/2560+t7979/640-t7986/160
      wI114pc1 = t18870+t18889+t18909+t18928+t18949+t18968+t18988+t19007
      wI115pc1 = wV31
      t19011 = 17.D0/5632.D0*t7999
      t19014 = 17.D0/5632.D0*t8008
      t19017 = 17.D0/2816.D0*t8015
      t19019 = 17.D0/2816.D0*t8019
      t19022 = -3.D0/88.D0*t7996-t19011+t8002/176-5.D0/704.D0*t8005+t190
     #14-t8011/88-3.D0/704.D0*t8013+t19017-t8017/88-t19019-5.D0/1408.D0*
     #t8021+5.D0/2816.D0*t8023
      t19023 = 17.D0/5632.D0*t8026
      t19030 = 17.D0/22528.D0*t8050
      t19035 = t19023-t8028/2816+3.D0/88.D0*t8035-5.D0/352.D0*t8038-5.D0
     #/704.D0*t8041+t8043/44-3.D0/176.D0*t8048-t19030+t8053/88-3.D0/352.
     #D0*t8057-25.D0/5632.D0*t8062+t8064/5632
      t19049 = t8066/22+5.D0/704.D0*t8068-7.D0/5632.D0*t8072-t8074/11-3.
     #D0/88.D0*t8076+5.D0/704.D0*t8078+t8080/704-5.D0/352.D0*t8088+25.D0
     #/1408.D0*t8090-25.D0/1408.D0*t8095-t8097/88+5.D0/1408.D0*t8101
      t19055 = t8117/22
      t19058 = t8131/22
      t19060 = t8137/22
      t19062 = -5.D0/352.D0*t8103+3.D0/352.D0*t8105+t8109/22+5.D0/352.D0
     #*t8111-t8113/704-t19055+3.D0/88.D0*t8120+t8123/11+t19058-t8134/11-
     #t19060-t8141/22
      t19077 = 3.D0/176.D0*t8144+t8149/22-t8151/44-t8154/44+t8156/44+t81
     #58/44-t8161/22+t8165/11-7.D0/22528.D0*t8171-t8173/44-t8175/44+t817
     #7/44
      t19079 = 17.D0/5632.D0*t8184
      t19087 = t8214/11
      t19088 = t8218/11
      t19090 = -t8182/22+t19079-5.D0/176.D0*t8188-t8191/11+t8196/11-3.D0
     #/44.D0*t8202+5.D0/176.D0*t8205-25.D0/704.D0*t8208+25.D0/704.D0*t82
     #10+t19087-t19088-2.D0/11.D0*t8220
      t19092 = t8225/11
      t19104 = t19092+2.D0/11.D0*t8227-t8230/1408+t8233/5632+t8235/704-t
     #8238/22-2.D0/11.D0*t8241+t8243/11+t8246/22-t8249/22+t8255/176+5.D0
     #/352.D0*t8259
      t19109 = t8277/22
      t19111 = t8283/22
      t19112 = t8286/22
      t19117 = -5.D0/352.D0*t8262+25.D0/1408.D0*t8265-25.D0/1408.D0*t827
     #1+3.D0/88.D0*t8273-t19109+t8279/11+t19111-t19112-t8288/1408-t8290/
     #11-t8292/22+t8294/11
      t19123 = 17.D0/5632.D0*t8301
      t19129 = 17.D0/5632.D0*t8317
      t19133 = -5.D0/704.D0*t8296+t8298/22+t19123-t8303/44-t8306/44-5.D0
     #/704.D0*t8310-7.D0/2816.D0*t8313+25.D0/2816.D0*t8315-t19129-3.D0/1
     #76.D0*t8319-3.D0/88.D0*t8321+t8325/1408
      t19143 = t8353/44
      t19146 = 7.D0/5632.D0*t8328-25.D0/5632.D0*t8332+5.D0/1408.D0*t8334
     #+7.D0/5632.D0*t8337-25.D0/704.D0*t8340+25.D0/1408.D0*t8344-7.D0/28
     #16.D0*t8346+t8349/22-t8351/88-t19143-t8355/11-3.D0/352.D0*t8358
      t19148 = 17.D0/11264.D0*t8360
      t19160 = t19148+t8368/88-5.D0/352.D0*t8371-7.D0/5632.D0*t8373+25.D
     #0/1408.D0*t8376+3.D0/176.D0*t8378+7.D0/5632.D0*t8380+t8382/22+t838
     #4/22-t8386/1408-t8388/22-t8390/22
      t19165 = t7992*t1067
      t19166 = t3606*t7176
      t19167 = t19165*t19166
      t19176 = t7992*t175*t3406
      t19177 = t3606*t7104
      t19178 = t19176*t19177
      t19180 = t8397/11+3.D0/88.D0*t8399+5.D0/352.D0*t8402-t8407/44+3.D0
     #/88.D0*t19167+t8410/22+3.D0/44.D0*t8413+t8415/22+t8417/1408-3.D0/4
     #4.D0*t8420+5.D0/2816.D0*t8422-3.D0/176.D0*t19178
      t19184 = t3606*t7082
      t19185 = t19165*t19184
      t19188 = t7992*t1067*t3426
      t19189 = t19188*t19166
      t19191 = t19188*t19184
      t19196 = t8010*t7026
      t19199 = t19176*t7376*t4025
      t19204 = 3.D0/176.D0*t8428-3.D0/176.D0*t19185-3.D0/352.D0*t19189+3
     #.D0/704.D0*t19191-t8433/44+7.D0/11264.D0*t8435+t8439/176+3.D0/176.
     #D0*t19196-3.D0/176.D0*t19199-3.D0/88.D0*t8442-7.D0/22528.D0*t8447-
     #5.D0/704.D0*t8449
      t19205 = t6943*Mzsq
      t19207 = t19165*t9743*t19205
      t19209 = t3606*t8055
      t19210 = t12953*t19209
      t19212 = t8010*t6976
      t19218 = t19176*t7376*t3795
      t19220 = t19165*t7138
      t19225 = t8122*t8083
      t19227 = 3.D0/704.D0*t19207-3.D0/88.D0*t19210-3.D0/88.D0*t19212-t8
     #456/44+t8458/44+3.D0/88.D0*t8460+3.D0/88.D0*t19218-3.D0/176.D0*t19
     #220-3.D0/704.D0*t8464+7.D0/5632.D0*t8467+t8471/176+3.D0/88.D0*t192
     #25
      t19231 = t7992*cw2*t3406*t3420
      t19232 = t19231*t8083
      t19235 = t7992*t175*t3426
      t19236 = t19235*t19209
      t19238 = t8148*t8412
      t19240 = t8216*t8031
      t19241 = t19231*t19240
      t19243 = t8148*t8441
      t19245 = t8281*t8031
      t19246 = t19231*t19245
      t19248 = t8010*t8199
      t19250 = 17.D0/2816.D0*t8475
      t19253 = t8010*t8252
      t19256 = t19235*t7376*t8198
      t19258 = -3.D0/88.D0*t19232+3.D0/352.D0*t19236-3.D0/44.D0*t19238+3
     #.D0/44.D0*t19241+3.D0/88.D0*t19243-3.D0/88.D0*t19246+3.D0/44.D0*t1
     #9248-t19250+5.D0/352.D0*t8479-t8482/1408-3.D0/88.D0*t19253-3.D0/17
     #6.D0*t19256
      t19260 = t19235*t7376*t8251
      t19262 = t12953*t19177
      t19264 = t8486/88
      t19266 = t8490/88
      t19271 = 17.D0/22528.D0*t8503
      t19274 = 3.D0/352.D0*t19260+3.D0/176.D0*t19262+t19264-t8488/44+t19
     #266-t8492/44+t8494/22+3.D0/176.D0*t8498+3.D0/176.D0*t8501-t19271+5
     #.D0/176.D0*t8506+7.D0/2816.D0*t8508
      wI116pc1 = t19022+t19035+t19049+t19062+t19077+t19090+t19104+t19117
     #+t19133+t19146+t19160+t19180+t19204+t19227+t19258+t19274
      t19290 = t8521/16+t8525/32+t8529/256-t8532/64-t8536/64+t8538/128+t
     #8542/256+t8546/8-t8552/16-t8556/32+7.D0/64.D0*t8560+5.D0/32.D0*t85
     #65
      t19303 = t8571/32-t8576/8-7.D0/64.D0*t8581+7.D0/64.D0*t8585-5.D0/3
     #2.D0*t8589+5.D0/32.D0*t8591-t8593/32+t8595/8+t8598/4-t8600/16+5.D0
     #/16.D0*t8603-7.D0/32.D0*t8605
      t19317 = t8610/128-t8612/32+t8614/16+7.D0/32.D0*t8616-t8619/16+t86
     #23/32+t8626/32-5.D0/16.D0*t8628+t8630/16-t8634/64-t8636/4+t8641/32
      t19330 = -t8646/32-7.D0/64.D0*t8648-t8651/128-7.D0/256.D0*t8655+t8
     #659/64+7.D0/128.D0*t8663-t8667/128-t8671/128-t8677/16-t8682/16+t86
     #86/16+t8689/16
      t19345 = -5.D0/32.D0*t8697-t8700/8-t8702/16-7.D0/256.D0*t8705+t870
     #9/32-t8712/64+t8714/16+t8717/32+7.D0/64.D0*t8719-5.D0/64.D0*t8723+
     #t8727/32-t8729/32
      t19358 = -t8734/32-t8736/64+5.D0/128.D0*t8738-t8742/64+t8746/32+t8
     #749/8-t8752/8+t8755/64+t8757/64+5.D0/16.D0*t8760+t8763/4+t8765/32
      t19372 = t8769/8+t8771/16-t8774/8-t8776/8-t8779/16-7.D0/32.D0*t878
     #1-t8783/16-t8785/16+7.D0/64.D0*t8789+5.D0/128.D0*t8791+t8793/8-t87
     #97/64
      t19385 = t8800/64+t8803/16-t8806/16-t8808/32-5.D0/32.D0*t8810+t881
     #2/16-t8814/8+t8816/64-5.D0/32.D0*t8818-t8821/8+t8823/32-t8826/32
      wI117pc1 = t19290+t19303+t19317+t19330+t19345+t19358+t19372+t19385
      t19401 = 5.D0/18432.D0*t8835+5.D0/2304.D0*t8838+t8842/128-t8845/64
     #-7.D0/1536.D0*t8848+5.D0/576.D0*t8850-7.D0/384.D0*t8853-5.D0/4608.
     #D0*t8855-5.D0/1152.D0*t8857-7.D0/384.D0*t8860+7.D0/768.D0*t8862-7.
     #D0/1536.D0*t8864+t8868/4608
      t19416 = 5.D0/4608.D0*t8872-5.D0/144.D0*t8876-7.D0/768.D0*t8879+5.
     #D0/288.D0*t8882+5.D0/144.D0*t8886+t8888/64+7.D0/384.D0*t8891+13.D0
     #/1152.D0*t8895-t8897/32-t8899/288-t8905/16+5.D0/288.D0*t8908-5.D0/
     #288.D0*t8910-5.D0/4608.D0*t8913
      t19431 = t8917/288-5.D0/2304.D0*t8919-5.D0/4608.D0*t8922+7.D0/384.
     #D0*t8924+7.D0/384.D0*t8926-5.D0/2304.D0*t8929+t8932/1152-t8935/115
     #2-t8938/1152-t8940/512+5.D0/2304.D0*t8945+13.D0/1152.D0*t8948-t895
     #2/1152
      t19446 = -5.D0/4608.D0*t8955-5.D0/1152.D0*t8958-5.D0/576.D0*t8962+
     #5.D0/1152.D0*t8965-5.D0/4608.D0*t8968+t8971/288-7.D0/384.D0*t8973-
     #5.D0/2304.D0*t8976+t8978/128-5.D0/1152.D0*t8982+5.D0/2304.D0*t8985
     #+5.D0/288.D0*t8990-5.D0/288.D0*t8993+5.D0/144.D0*t8996
      t19462 = -13.D0/288.D0*t9002-t9005/32-5.D0/2304.D0*t9007-7.D0/192.
     #D0*t9009+5.D0/4608.D0*t9011+5.D0/288.D0*t9013+5.D0/4608.D0*t9015-5
     #.D0/288.D0*t9017+5.D0/576.D0*t9020-5.D0/4608.D0*t9023-5.D0/1152.D0
     #*t9025+5.D0/1152.D0*t9028-5.D0/1152.D0*t9031
      t19477 = 5.D0/144.D0*t9034+7.D0/384.D0*t9036-7.D0/384.D0*t9039-7.D
     #0/192.D0*t9041-5.D0/288.D0*t9044-7.D0/384.D0*t9046+7.D0/192.D0*t90
     #48-7.D0/1536.D0*t9050-5.D0/144.D0*t9053-7.D0/384.D0*t9056+13.D0/28
     #8.D0*t9059+5.D0/576.D0*t9062+5.D0/2304.D0*t9064-13.D0/1152.D0*t906
     #7
      t19492 = 5.D0/288.D0*t9072-t9075/576-5.D0/1152.D0*t9078+5.D0/4608.
     #D0*t9081+t9083/576+t9086/1152+7.D0/384.D0*t9088+7.D0/192.D0*t9090-
     #5.D0/1152.D0*t9092+5.D0/288.D0*t9094-t9096/1152-7.D0/384.D0*t9100+
     #5.D0/288.D0*t9102
      t19507 = 7.D0/384.D0*t9105-5.D0/144.D0*t9107-13.D0/576.D0*t9110+13
     #.D0/144.D0*t9113-7.D0/192.D0*t9115+5.D0/288.D0*t9118-5.D0/72.D0*t9
     #121+t9123/144-t9126/144-13.D0/288.D0*t9129-13.D0/144.D0*t9131+5.D0
     #/72.D0*t9133+7.D0/192.D0*t9136+5.D0/576.D0*t9141
      t19524 = -t9148/2304+13.D0/1152.D0*t9150-5.D0/1152.D0*t9153-7.D0/3
     #84.D0*t9155+t9158/4608+5.D0/18432.D0*t9160-7.D0/768.D0*t9162-5.D0/
     #4608.D0*t9164+t9166/256-5.D0/2304.D0*t9169+7.D0/384.D0*t9171-t9173
     #/512+7.D0/1536.D0*t9175
      t19539 = 7.D0/384.D0*t9178+5.D0/4608.D0*t9181+t9184/128+t9186/576-
     #t9188/1152+5.D0/144.D0*t9190-5.D0/288.D0*t9192-5.D0/144.D0*t9196+5
     #.D0/288.D0*t9198-5.D0/576.D0*t9201+5.D0/1152.D0*t9204+5.D0/1152.D0
     #*t9207+7.D0/1536.D0*t9209-5.D0/144.D0*t9211
      t19554 = 5.D0/288.D0*t9215+5.D0/18432.D0*t9217-5.D0/2304.D0*t9220-
     #5.D0/1152.D0*t9222+5.D0/2304.D0*t9225+7.D0/192.D0*t9227-t9229/128-
     #13.D0/1152.D0*t9231+13.D0/576.D0*t9233+t9235/128+7.D0/1536.D0*t923
     #7-5.D0/1152.D0*t9241-7.D0/192.D0*t9243
      t19569 = -13.D0/4608.D0*t9248+7.D0/768.D0*t9250-13.D0/4608.D0*t925
     #3+13.D0/2304.D0*t9256-7.D0/1536.D0*t9259-5.D0/9216.D0*t9261-7.D0/3
     #84.D0*t9263+7.D0/384.D0*t9265+5.D0/288.D0*t9268-7.D0/384.D0*t9270+
     #5.D0/144.D0*t9273-t9276/32+7.D0/384.D0*t9279-5.D0/288.D0*t9282
      t19585 = t9288/288-t9290/64-13.D0/576.D0*t9293+5.D0/2304.D0*t9295+
     #t9298/576+5.D0/576.D0*t9301+t9303/32+13.D0/288.D0*t9306-7.D0/384.D
     #0*t9308+13.D0/1152.D0*t9311-13.D0/288.D0*t9313+t9316/288+7.D0/384.
     #D0*t9318
      t19600 = 5.D0/144.D0*t9322-t9325/288-5.D0/144.D0*t9327-5.D0/2304.D
     #0*t9329+7.D0/1536.D0*t9331+t9334/128-5.D0/4608.D0*t9336+13.D0/1152
     #.D0*t9339-t9341/1152+t9343/128-5.D0/1152.D0*t9346-13.D0/288.D0*t93
     #48-7.D0/384.D0*t9350+5.D0/288.D0*t9353
      t19615 = 5.D0/288.D0*t9358-5.D0/4608.D0*t9360+5.D0/9216.D0*t9362-5
     #.D0/9216.D0*t9364+5.D0/9216.D0*t9368+7.D0/384.D0*t9371+5.D0/4608.D
     #0*t9373-5.D0/576.D0*t9375-5.D0/144.D0*t9378-5.D0/72.D0*t9380+t9383
     #/16+13.D0/144.D0*t9385+7.D0/192.D0*t9387
      t19630 = 5.D0/18432.D0*t9392+7.D0/192.D0*t9394-t9396/128-7.D0/192.
     #D0*t9398-t9400/64+t9403/32-t9406/32+5.D0/144.D0*t9409-5.D0/144.D0*
     #t9412+t9415/16+13.D0/1152.D0*t9417-7.D0/192.D0*t9420-13.D0/576.D0*
     #t9422-t9424/144
      wI118pc1 = t19401+t19416+t19431+t19446+t19462+t19477+t19492+t19507
     #+t19524+t19539+t19554+t19569+t19585+t19600+t19615+t19630
      t19649 = -11.D0/32.D0*t9440+t9446/64-11.D0/256.D0*t9452+11.D0/32.D
     #0*t9458-11.D0/64.D0*t9462+t9467/8-t9470/16+3.D0/32.D0*t9476+t9484/
     #8+t9487/4-t9490/8-t9495/16-t9497/8-t9500/8+t9502/16
      t19665 = 11.D0/64.D0*t9506+t9511/4+t9513/32+t9518/32-t9520/64+t952
     #2/16-t9525/8+3.D0/32.D0*t9529+t9531/8+3.D0/32.D0*t9533-t9536/8-3.D
     #0/16.D0*t9538-t9540/16+t9542/16+3.D0/16.D0*t9550
      t19682 = -11.D0/32.D0*t9555+11.D0/64.D0*t9557-3.D0/32.D0*t9560+t95
     #65/16+t9568/32+3.D0/64.D0*t9571+11.D0/64.D0*t9573+t9575/16-3.D0/12
     #8.D0*t9577-t9579/16+t9582/64-t9584/4+11.D0/128.D0*t9587-11.D0/256.
     #D0*t9589-t9591/8
      t19698 = -11.D0/64.D0*t9596+3.D0/32.D0*t9599+t9602/8-3.D0/128.D0*t
     #9605-t9607/16-t9609/64-3.D0/32.D0*t9614-t9619/16+t9622/16-3.D0/16.
     #D0*t9625+11.D0/64.D0*t9628-t9630/16+t9632/8-t9634/32-t9636/8
      wI119pc1 = t19649+t19665+t19682+t19698
      t19701 = t9651/1792
      t19702 = t9664/1792
      t19704 = t9671/3584
      t19705 = t9676/1792
      t19707 = t9689/3584
      t19708 = 5.D0/4032.D0*t9697
      t19709 = t9701/14336
      t19710 = t9705/192
      t19711 = t9710/14336
      t19712 = 5.D0/8064.D0*t9712
      t19713 = 5.D0/4032.D0*t9716
      t19714 = -3.D0/56.D0*t9645-t19701+t19702-17.D0/4032.D0*t9668-t1970
     #4+t19705-t9685/224-t19707-t19708+t19709+t19710+t19711-t19712-t1971
     #3
      t19717 = 5.D0/672.D0*t9729
      t19718 = 5.D0/672.D0*t9732
      t19719 = 13.D0/8064.D0*t9736
      t19724 = 13.D0/4032.D0*t9752
      t19725 = 5.D0/8064.D0*t9755
      t19728 = 5.D0/2016.D0*t9770
      t19730 = t9721/224-5.D0/504.D0*t9725-t19717+t19718+t19719-17.D0/10
     #08.D0*t9741-t9745/32256+5.D0/504.D0*t9748-t9750/224-t19724-t19725+
     #17.D0/504.D0*t9759+67.D0/4032.D0*t9763-t19728-127.D0/8064.D0*t9775
      t19732 = t9784/448
      t19733 = t9787/1792
      t19734 = t9797/48
      t19737 = 5.D0/2016.D0*t9808
      t19738 = 5.D0/4032.D0*t9812
      t19739 = t9815/112
      t19740 = 5.D0/16128.D0*t9823
      t19741 = 19.D0/4032.D0*t9827
      t19742 = t9830/112
      t19743 = 13.D0/16128.D0*t9837
      t19744 = 41.D0/2688.D0*t9841
      t19745 = 5.D0/16128.D0*t9843
      t19746 = t9847/1008
      t19747 = t19732+t19733+t19734-17.D0/2016.D0*t9801+17.D0/2016.D0*t9
     #804+t19737+t19738-t19739-t19740-t19741+t19742+t19743+t19744+t19745
     #-t19746
      t19748 = 5.D0/504.D0*t9850
      t19749 = 23.D0/1008.D0*t9853
      t19750 = t9855/224
      t19751 = t9858/252
      t19752 = t9862/252
      t19753 = 5.D0/1008.D0*t9864
      t19754 = t9866/252
      t19755 = t9868/1008
      t19756 = 41.D0/2688.D0*t9871
      t19759 = t9885/224
      t19761 = t9891/1008
      t19762 = 5.D0/1008.D0*t9893
      t19763 = t19748-t19749+t19750-t19751+t19752-t19753+t19754+t19755-t
     #19756-t9873/4032-61.D0/2016.D0*t9883-t19759+109.D0/4032.D0*t9889-t
     #19761+t19762
      t19767 = 5.D0/1008.D0*t9899
      t19768 = 23.D0/1008.D0*t9901
      t19772 = 23.D0/1008.D0*t9911
      t19773 = t9914/96
      t19774 = t9916/96
      t19776 = 19.D0/4032.D0*t9922
      t19777 = 41.D0/2688.D0*t9925
      t19778 = 5.D0/1008.D0*t9927
      t19779 = t9929/112
      t19780 = -t9895/48-t19767+t19768+13.D0/1008.D0*t9903-t9906/112+t99
     #08/768-t19772-t19773+t19774+t9919/112+t19776+t19777-t19778+t19779
      t19781 = 19.D0/4032.D0*t9931
      t19784 = 5.D0/2016.D0*t9939
      t19785 = 5.D0/504.D0*t9941
      t19787 = 5.D0/8064.D0*t9945
      t19789 = 5.D0/504.D0*t9949
      t19791 = 59.D0/64512.D0*t9959
      t19792 = 5.D0/1008.D0*t9962
      t19793 = 5.D0/1008.D0*t9967
      t19794 = t9969/112
      t19795 = 5.D0/2016.D0*t9971
      t19796 = -t19781-t9935/1008+t9937/112+t19784-t19785+t9943/1008+t19
     #787-t9947/1008+t19789+127.D0/64512.D0*t9954-t19791-t19792+t19793-t
     #19794+t19795
      t19798 = 5.D0/16128.D0*t9974
      t19800 = t9980/192
      t19802 = t9985/192
      t19803 = t9987/112
      t19806 = 5.D0/4032.D0*t9995
      t19808 = 5.D0/4032.D0*t9999
      t19809 = t10001/192
      t19810 = t10003/192
      t19811 = 5.D0/504.D0*t10008
      t19812 = 13.D0/32256.D0*t10010
      t19813 = -t19798+t9977/224+t19800-t9982/224-t19802+t19803-t9991/38
     #4-t9993/2016-t19806+t9997/2016+t19808+t19809-t19810-t19811+t19812
      t19814 = t10013/96
      t19817 = 5.D0/504.D0*t10023
      t19818 = t10025/96
      t19819 = t10029/224
      t19820 = t10032/448
      t19821 = 5.D0/16128.D0*t10034
      t19822 = 5.D0/672.D0*t10037
      t19823 = t10039/224
      t19824 = t10046/2688
      t19826 = t10051/192
      t19827 = t10061/1344
      t19828 = t10063/1344
      t19829 = t19814+t10016/112+t10019/224+t19817-t19818+t19819-t19820-
     #t19821+t19822+t19823+t19824-t10049/224-t19826-t19827+t19828
      t19834 = t10072/1792
      t19835 = t10075/896
      t19836 = t10079/1792
      t19837 = t10085/896
      t19838 = t10089/896
      t19839 = t10095/2688
      t19840 = t10100/896
      t19841 = t10109/448
      t19843 = t10123/10752
      t19844 = t10129/126
      t19845 = t10132/10752
      t19846 = t10138/448
      t19847 = 5.D0/126.D0*t10067-t19834-t19835+t19836-t19837+t19838-t19
     #839+t19840+t19841+191.D0/8064.D0*t10114+t19843+t19844-t19845-t1984
     #6
      t19848 = 5.D0/4032.D0*t10140
      t19849 = t10142/3584
      t19850 = t10147/192
      t19851 = t10151/448
      t19852 = t10153/63
      t19853 = t10157/3584
      t19854 = t10159/126
      t19855 = t10162/63
      t19857 = 5.D0/4032.D0*t10167
      t19858 = 13.D0/129024.D0*t10171
      t19859 = t10173/192
      t19861 = 5.D0/4032.D0*t10178
      t19863 = t19848-t19849-t19850+t19851+t19852+t19853+t19854-t19855-5
     #.D0/126.D0*t10164+t19857-t19858+t19859+t10176/2016-t19861+t10184/1
     #12
      t19865 = 5.D0/64512.D0*t10187
      t19866 = 5.D0/2016.D0*t10192
      t19867 = t10196/384
      t19868 = t10198/896
      t19871 = 5.D0/8064.D0*t10207
      t19872 = t10209/96
      t19873 = 5.D0/2016.D0*t10212
      t19874 = 13.D0/32256.D0*t10215
      t19875 = 5.D0/2016.D0*t10220
      t19876 = 19.D0/8064.D0*t10224
      t19877 = 13.D0/32256.D0*t10227
      t19878 = 5.D0/32256.D0*t10230
      t19880 = -t19865-t19866+t19867-t19868-t10201/896+17.D0/8064.D0*t10
     #203-t19871-t19872-t19873+t19874-t19875-t19876+t19877+t19878-191.D0
     #/8064.D0*t10234
      t19881 = 5.D0/8064.D0*t10236
      t19882 = 41.D0/10752.D0*t10238
      t19883 = t10241/768
      t19884 = t10243/112
      t19885 = 5.D0/8064.D0*t10247
      t19887 = 5.D0/8064.D0*t10254
      t19888 = 5.D0/4032.D0*t10259
      t19889 = t10262/768
      t19890 = 5.D0/4032.D0*t10264
      t19891 = 5.D0/2016.D0*t10266
      t19892 = t10268/224
      t19893 = t10270/1008
      t19894 = 3.D0/448.D0*t10273
      t19896 = -t19881-t19882-t19883-t19884+t19885-17.D0/2016.D0*t10250-
     #t19887-t19888+t19889+t19890-t19891+t19892-t19893+t19894-109.D0/806
     #4.D0*t10275
      t19900 = 41.D0/5376.D0*t10281
      t19901 = 5.D0/1008.D0*t10283
      t19902 = 5.D0/4032.D0*t10285
      t19904 = 5.D0/252.D0*t10289
      t19905 = 13.D0/16128.D0*t10292
      t19906 = 5.D0/2016.D0*t10296
      t19907 = t10299/384
      t19908 = t10301/192
      t19909 = t10306/3584
      t19911 = t10312/5376
      t19912 = 5.D0/1008.D0*t10314
      t19913 = t10317/5376
      t19914 = -61.D0/16128.D0*t10279+t19900+t19901+t19902+t10287/224-t1
     #9904-t19905-t19906-t19907+t19908-t19909-t10308/224-t19911+t19912-t
     #19913
      t19915 = 5.D0/672.D0*t10321
      t19916 = t10323/448
      t19920 = t10340/224
      t19930 = -t19915-t19916+9.D0/1792.D0*t10327+3.D0/448.D0*t10333-t10
     #337/448-t19920+t10346/224+9.D0/448.D0*t10350+t10354/384+t10357/448
     #+t10361/1344-t10365/336+t10369/896+t10373/896+13.D0/2688.D0*t10377
      t19937 = t10394/56
      t19939 = t10402/112
      t19940 = t10405/112
      t19944 = t10415/56
      t19947 = 13.D0/5376.D0*t10381-t10385/448+t10387/112+t10389/1344-t1
     #0391/96+t19937-3.D0/224.D0*t10396-t19939+t19940+3.D0/224.D0*t10407
     #-t10409/224+t10412/112+t19944+t10417/112-3.D0/112.D0*t10419
      t19949 = t10424/56
      t19957 = 5.D0/672.D0*t10453
      t19958 = 5.D0/16128.D0*t10456
      t19959 = 5.D0/672.D0*t10459
      t19963 = 3.D0/112.D0*t10422-t19949+t10427/48-9.D0/448.D0*t10430+t1
     #0432/96-t10436/56+t10440/224-t10442/2016+t10449/2016+t19957-t19958
     #-t19959+67.D0/8064.D0*t10461+61.D0/4032.D0*t10463+t10465/224
      t19969 = t10471/1792
      t19970 = t10476/224
      t19971 = t10478/48
      t19972 = 5.D0/504.D0*t10481
      t19973 = 5.D0/8064.D0*t10484
      t19974 = t10486/96
      t19975 = t10490/896
      t19976 = 5.D0/16128.D0*t10494
      t19977 = t10497/896
      t19978 = 5.D0/4032.D0*t10506
      t19979 = 41.D0/2688.D0*t10509
      t19980 = 5.D0/4032.D0*t10511
      t19981 = t10513/448
      t19982 = -t10467/896+t19969+t19970+t19971+t19972-t19973-t19974-t19
     #975-t19976-t19977+t19978+t19979-t19980+t19981
      t19983 = 5.D0/1008.D0*t10515
      t19984 = t10518/10752
      t19985 = t10520/96
      t19986 = t10523/96
      t19987 = t10527/10752
      t19988 = t10530/10752
      t19989 = 13.D0/32256.D0*t10532
      t19990 = t10536/112
      t19991 = t10538/192
      t19992 = t10541/768
      t19993 = 3.D0/448.D0*t10543
      t19997 = t10560/192
      t19998 = -t19983+t19984+t19985+t19986+t19987+t19988-t19989+t19990-
     #t19991+t19992+t19993-13.D0/504.D0*t10550+67.D0/8064.D0*t10554-67.D
     #0/4032.D0*t10558+t19997
      t20000 = t10562/96
      t20001 = t10564/384
      t20004 = 5.D0/4032.D0*t10571
      t20013 = 5.D0/672.D0*t10596
      t20014 = t10599/768
      t20015 = -t20000+t20001+t10567/2016+17.D0/1008.D0*t10569-t20004-13
     #.D0/1008.D0*t10573-67.D0/8064.D0*t10576-17.D0/1008.D0*t10579-t1058
     #1/224-t10585/448+t10587/112-9.D0/896.D0*t10590-127.D0/16128.D0*t10
     #593-t20013-t20014
      t20016 = 59.D0/16128.D0*t10602
      t20023 = t10625/112
      t20031 = t20016-67.D0/32256.D0*t10604-t10607/224-t10610/448-t10614
     #/56+3.D0/224.D0*t10617-3.D0/112.D0*t10619-t20023-53.D0/2688.D0*t10
     #627-191.D0/4032.D0*t10629-191.D0/8064.D0*t10631+t10633/112-61.D0/4
     #032.D0*t10635-t10637/112+9.D0/896.D0*t10639
      t20035 = 5.D0/504.D0*t10643
      t20037 = t10647/896
      t20046 = 3.D0/1792.D0*t10673
      t20048 = -t10641/224+t20035+11.D0/896.D0*t10645+t20037+191.D0/4032
     #.D0*t10649+61.D0/2016.D0*t10651+11.D0/448.D0*t10653-t10657/224+53.
     #D0/1344.D0*t10660+109.D0/4032.D0*t10662-t10666/168+t10669/56-t2004
     #6-9.D0/448.D0*t10675
      t20053 = t10687/2688
      t20054 = t10692/448
      t20058 = 5.D0/4032.D0*t10704
      t20064 = 3.D0/224.D0*t10678-13.D0/2688.D0*t10680-t10682/224+t10684
     #/8064-t20053-t20054+109.D0/8064.D0*t10694-t10697/56+t10701/56+t200
     #58+53.D0/1344.D0*t10706-53.D0/1344.D0*t10708-3.D0/896.D0*t10710+10
     #9.D0/32256.D0*t10713-13.D0/2688.D0*t10715
      t20068 = 13.D0/8064.D0*t10721
      t20073 = t10735/448
      t20075 = t10742/896
      t20078 = t10752/448
      t20081 = t10717/168-3.D0/112.D0*t10719-t20068+13.D0/4032.D0*t10725
     #-5.D0/2016.D0*t10729+t10731/48-t10733/48+t20073+61.D0/8064.D0*t107
     #37-t20075+191.D0/32256.D0*t10744+11.D0/3584.D0*t10746-t20078+t1075
     #4/896-3.D0/448.D0*t10756
      t20084 = 3.D0/224.D0*t10762
      t20085 = t10764/192
      t20086 = 5.D0/16128.D0*t10766
      t20087 = 5.D0/4032.D0*t10768
      t20088 = 5.D0/2016.D0*t10770
      t20090 = t10774/896
      t20091 = t10780/896
      t20096 = 5.D0/16128.D0*t10794
      t20097 = -t10758/112+17.D0/1008.D0*t10760-t20084-t20085+t20086-t20
     #087+t20088-17.D0/2016.D0*t10772+t20090+t20091+17.D0/8064.D0*t10782
     #-9.D0/3584.D0*t10785-t10789/448+53.D0/10752.D0*t10792+t20096
      t20101 = t10796/2688
      t20102 = t10798/768
      t20103 = t10800/384
      t20105 = 5.D0/16128.D0*t10805
      t20110 = t10820/1792
      t20115 = t20101-t20102+t20103+17.D0/4032.D0*t10802-t20105+3.D0/112
     #.D0*t10807-3.D0/112.D0*t10809-t10816/56-13.D0/21504.D0*t10818+t201
     #10+t10822/168+13.D0/5376.D0*t10826-53.D0/2688.D0*t10829-t10831/112
      t20116 = 13.D0/8064.D0*t10833
      t20118 = t10838/2688
      t20119 = t10844/112
      t20131 = -t20116-t10836/224-t20118-t20119-t10848/448+t10851/448+t1
     #0853/112-t10858/768+17.D0/1008.D0*t10861-t10863/224-t10865/3072-12
     #7.D0/16128.D0*t10867-13.D0/5376.D0*t10870+9.D0/896.D0*t10873-9.D0/
     #896.D0*t10876
      t20138 = 5.D0/4032.D0*t10892
      t20141 = t10900/3584
      t20148 = 9.D0/896.D0*t10879+13.D0/5376.D0*t10884-t10886/336+t10888
     #/336-17.D0/1008.D0*t10890-t20138+t10894/224+53.D0/2688.D0*t10897-t
     #20141-13.D0/5376.D0*t10905-t10909/896-53.D0/2688.D0*t10911+t10913/
     #112-t10915/112-t10918/96
      t20152 = 13.D0/8064.D0*t10927
      t20154 = t10932/112
      t20155 = 3.D0/896.D0*t10935
      t20158 = 5.D0/672.D0*t10945
      t20159 = t10948/192
      t20162 = t10958/112
      t20163 = 5.D0/672.D0*t10961
      t20164 = -t10921/112-t10923/96+t10925/96+t20152-61.D0/16128.D0*t10
     #929+t20154+t20155-191.D0/16128.D0*t10940+t10943/224+t20158+t20159+
     #t10952/1792-t10955/448-t20162-t20163
      t20167 = t10964/96
      t20168 = 13.D0/32256.D0*t10966
      t20169 = t10969/3584
      t20170 = t10971/10752
      t20171 = 5.D0/16128.D0*t10975
      t20172 = 59.D0/16128.D0*t10977
      t20173 = t10979/896
      t20177 = 41.D0/10752.D0*t10991
      t20178 = t10995/4032
      t20179 = t10997/4032
      t20182 = -t20167-t20168-t20169+t20170+t20171-t20172+t20173-t10984/
     #224+3.D0/224.D0*t10986-3.D0/224.D0*t10988-t20177-t20178+t20179+t10
     #999/224+t11001/1792
      t20185 = 5.D0/2016.D0*t11008
      t20187 = t11013/192
      t20188 = t11016/126
      t20189 = t11020/168
      t20190 = t11026/168
      t20192 = 5.D0/32256.D0*t11031
      t20194 = 5.D0/2016.D0*t11042
      t20196 = 3.D0/1792.D0*t11048
      t20198 = -t11003/336+t11005/112+t20185-5.D0/504.D0*t11010-t20187-t
     #20188-t20189+t20190-3.D0/56.D0*t11028-t20192+13.D0/504.D0*t11033-t
     #20194+t11046/384-t20196-11.D0/1792.D0*t11050
      t20200 = 17.D0/32256.D0*t11060
      t20201 = 17.D0/64512.D0*t11064
      t20203 = t11069/168
      t20204 = 5.D0/2016.D0*t11071
      t20205 = t11077/504
      t20206 = 41.D0/1344.D0*t11080
      t20208 = 13.D0/4032.D0*t11084
      t20209 = t11087/3584
      t20211 = t11091/504
      t20212 = 23.D0/2016.D0*t11093
      t20214 = t11100/252
      t20215 = -t20200+t20201+5.D0/504.D0*t11066-t20203-t20204-t20205-t2
     #0206+3.D0/56.D0*t11082+t20208+t20209-11.D0/896.D0*t11089+t20211-t2
     #0212-5.D0/63.D0*t11096+t20214
      t20216 = t11107/4032
      t20218 = t11119/2016
      t20220 = t9857*t11035/168
      t20221 = 5.D0/252.D0*t11121
      t20222 = 23.D0/4032.D0*t11123
      t20223 = t11127/126
      t20224 = t11129/224
      t20225 = t11131/504
      t20226 = 5.D0/252.D0*t11134
      t20227 = 13.D0/4032.D0*t11136
      t20228 = 5.D0/336.D0*t11138
      t20232 = t20216+t11109/384-t20218-t20220+t20221+t20222+t20223+t202
     #24-t20225-t20226-t20227-t20228-17.D0/504.D0*t11141+17.D0/504.D0*t1
     #1143-t11145/112
      t20238 = 5.D0/1008.D0*t11147
      t20241 = t11155/448
      t20244 = 19.D0/16128.D0*t11161
      t20245 = 5.D0/16128.D0*t11165
      t20247 = 5.D0/16128.D0*t11169
      t20248 = 5.D0/8064.D0*t11171
      t20251 = t11179/4032
      t20252 = -t20238-t11150/448+53.D0/10752.D0*t11153-t20241-t11157/22
     #4-191.D0/8064.D0*t11159+t20244-t20245-11.D0/896.D0*t11167+t20247+t
     #20248-61.D0/4032.D0*t11173+t11175/768+t20251
      t20254 = 5.D0/504.D0*t11183
      t20256 = 5.D0/4032.D0*t11187
      t20257 = 17.D0/8064.D0*t11192
      t20258 = 5.D0/4032.D0*t11194
      t20259 = t11196/448
      t20262 = 5.D0/8064.D0*t11202
      t20263 = 17.D0/16128.D0*t11205
      t20264 = t11207/8064
      t20265 = 5.D0/2016.D0*t11209
      t20266 = 5.D0/4032.D0*t11211
      t20267 = 23.D0/504.D0*t11213
      t20268 = -5.D0/63.D0*t11181+t20254-53.D0/2016.D0*t11185+t20256+t20
     #257+t20258-t20259-13.D0/21504.D0*t11198-9.D0/3584.D0*t11200-t20262
     #-t20263+t20264-t20265+t20266+t20267
      t20271 = 17.D0/64512.D0*t11219
      t20272 = 5.D0/252.D0*t11221
      t20273 = 5.D0/336.D0*t11224
      t20274 = 5.D0/8064.D0*t11228
      t20275 = 5.D0/4032.D0*t11232
      t20278 = 5.D0/504.D0*t11239
      t20279 = t11241/192
      t20280 = 5.D0/504.D0*t11243
      t20281 = t11249/168
      t20284 = t11255/10752
      t20285 = 53.D0/1008.D0*t11215+t20271-t20272+t20273-t20274-t20275-5
     #3.D0/1008.D0*t11234+53.D0/1008.D0*t11237+t20278+t20279+t20280+t202
     #81+61.D0/4032.D0*t11251+5.D0/252.D0*t11253-t20284
      t20286 = t11257/168
      t20288 = t11262/224
      t20289 = 19.D0/16128.D0*t11264
      t20290 = 13.D0/32256.D0*t11267
      t20291 = t11272/126
      t20292 = 5.D0/8064.D0*t11275
      t20293 = t11277/448
      t20295 = t11281/448
      t20296 = 5.D0/1008.D0*t11283
      t20297 = t11287/1792
      t20300 = t11294/84
      t20301 = -t20286+67.D0/8064.D0*t11260-t20288+t20289-t20290-t20291+
     #t20292-t20293+5.D0/63.D0*t11279-t20295-t20296-t20297+109.D0/8064.D
     #0*t11289-5.D0/252.D0*t11291+t20300
      t20305 = t11299/126
      t20307 = 5.D0/32256.D0*t11303
      t20309 = 5.D0/32256.D0*t11308
      t20311 = t11313/2016
      t20312 = 5.D0/8064.D0*t11315
      t20313 = t11317/504
      t20314 = t11319/1792
      t20315 = 5.D0/16128.D0*t11321
      t20316 = 19.D0/2016.D0*t11323
      t20318 = -127.D0/16128.D0*t11296+t20305+13.D0/1008.D0*t11301+t2030
     #7-13.D0/1008.D0*t11306-t20309+t11310/112+t20311-t20312+t20313+t203
     #14-t20315+t20316+t11325/768
      t20320 = 5.D0/672.D0*t11329
      t20322 = t9861*t11035/168
      t20324 = t4892*t12319*t3420
      t20327 = t20324*t8216*t11018/336
      t20329 = t11333/448
      t20331 = 5.D0/16128.D0*t11337
      t20334 = 5.D0/16128.D0*t11345
      t20335 = 5.D0/32256.D0*t11347
      t20339 = t9747*t11293/84
      t20340 = t11327/8064-t20320+t20322+t20327+61.D0/4032.D0*t11331+t20
     #329-109.D0/8064.D0*t11335+t20331+5.D0/504.D0*t11341+3.D0/112.D0*t1
     #1343-t20334+t20335+191.D0/8064.D0*t11349-13.D0/1008.D0*t11351-t203
     #39
      t20343 = t9747*t11019/168
      t20346 = t20324*t8281*t11018/672
      t20348 = 2.D0/63.D0*t11356
      t20350 = 2.D0/63.D0*t11360
      t20351 = 3.D0/448.D0*t11362
      t20353 = 5.D0/4032.D0*t11366
      t20355 = 5.D0/16128.D0*t11370
      t20356 = 5.D0/16128.D0*t11372
      t20357 = t11375/126
      t20358 = t11377/126
      t20359 = 13.D0/8064.D0*t11379
      t20360 = t20343-t20346+3.D0/112.D0*t11353+t20348+53.D0/2016.D0*t11
     #358-t20350-t20351-t11364/384+t20353+61.D0/4032.D0*t11368-t20355+t2
     #0356+t20357-t20358+t20359
      t20361 = 5.D0/4032.D0*t11381
      t20363 = 5.D0/1008.D0*t11385
      t20366 = t9861*t11052/168
      t20367 = t10065*t7376
      t20369 = t20367*t11025/168
      t20372 = t9642*t11038/168
      t20374 = t9861*t11686/84
      t20376 = t20324*t11038/672
      t20377 = t11391/448
      t20378 = 19.D0/4032.D0*t11393
      t20380 = t20367*t11689/84
      t20381 = 5.D0/8064.D0*t11395
      t20382 = 5.D0/16128.D0*t11398
      t20383 = t20361+t11383/96+t20363-t11387/768+t20366-t20369+127.D0/1
     #6128.D0*t11389+t20372-t20374-t20376+t20377-t20378+t20380+t20381+t2
     #0382
      t20387 = t11400/1008
      t20388 = 5.D0/32256.D0*t11402
      t20390 = 5.D0/32256.D0*t11406
      t20391 = t11414/1792
      t20392 = 5.D0/1008.D0*t11416
      t20394 = 5.D0/129024.D0*t11422
      t20395 = 5.D0/64512.D0*t11426
      t20396 = t11428/448
      t20398 = 5.D0/2688.D0*t11432
      t20399 = t11434/63
      t20400 = t11438/252
      t20401 = -t20387+t20388+t11404/768-t20390+t20391-t20392+17.D0/2016
     #.D0*t11418+t20394+t20395+t20396-191.D0/8064.D0*t11430+t20398-t2039
     #9+t20400
      t20402 = 5.D0/1008.D0*t11440
      t20403 = 5.D0/1008.D0*t11442
      t20404 = 5.D0/504.D0*t11444
      t20405 = 23.D0/1008.D0*t11447
      t20407 = 5.D0/1008.D0*t11452
      t20408 = t11454/112
      t20410 = t11458/16128
      t20411 = t11460/126
      t20412 = 59.D0/16128.D0*t11462
      t20415 = t11469/504
      t20416 = 5.D0/2688.D0*t11471
      t20417 = t20402-t20403+t20404-t20405+127.D0/64512.D0*t11449-t20407
     #-t20408-3.D0/896.D0*t11456-t20410+t20411+t20412-t11465/336+t11467/
     #16128-t20415+t20416
      t20422 = 3.D0/448.D0*t11484
      t20423 = t11488/63
      t20424 = t11490/63
      t20425 = t11492/126
      t20426 = t11494/126
      t20427 = t11496/63
      t20428 = t11498/504
      t20429 = 5.D0/1008.D0*t11500
      t20430 = 5.D0/504.D0*t11502
      t20431 = 5.D0/504.D0*t11504
      t20432 = 5.D0/2016.D0*t11507
      t20433 = t11509/126
      t20434 = t11473/896+9.D0/896.D0*t11475+t11477/336+t20422-t20423+t2
     #0424+t20425-t20426-t20427+t20428+t20429-t20430+t20431+t20432+t2043
     #3
      t20437 = 5.D0/1008.D0*t11518
      t20438 = t11520/1792
      t20439 = t11523/3584
      t20440 = t11528/4032
      t20441 = t11530/896
      t20442 = t11532/224
      t20445 = t11538/4032
      t20446 = t11542/168
      t20447 = t11544/10752
      t20448 = t11547/336
      t20449 = t11549/112
      t20450 = 13.D0/5376.D0*t11512+53.D0/2688.D0*t11515-t20437+t20438+t
     #20439-t20440+t20441+t20442-t11534/8064+13.D0/4032.D0*t11536+t20445
     #-t20446+t20447-t20448+t20449
      t20453 = t11552/192
      t20455 = 41.D0/1344.D0*t11556
      t20457 = 23.D0/504.D0*t11563
      t20458 = t11565/1008
      t20459 = 23.D0/504.D0*t11567
      t20460 = t11569/192
      t20461 = 5.D0/4032.D0*t11571
      t20464 = 17.D0/8064.D0*t11581
      t20466 = 13.D0/32256.D0*t11587
      t20467 = 5.D0/16128.D0*t11589
      t20468 = t20453-5.D0/2016.D0*t11554+t20455-53.D0/2688.D0*t11558+t2
     #0457+t20458-t20459-t20460+t20461+t11574/896-3.D0/448.D0*t11576-t20
     #464+17.D0/4032.D0*t11584+t20466-t20467
      t20470 = 5.D0/504.D0*t11595
      t20471 = 5.D0/504.D0*t11598
      t20474 = t11604/63
      t20475 = 2.D0/63.D0*t11607
      t20476 = t11609/63
      t20477 = 19.D0/2016.D0*t11611
      t20478 = 41.D0/1344.D0*t11613
      t20479 = t11615/252
      t20480 = 13.D0/32256.D0*t11617
      t20481 = t11620/96
      t20482 = t11622/63
      t20483 = t11624/224
      t20484 = -t11591/3072-t20470+t20471-t11600/112+t11602/112+t20474+t
     #20475-t20476-t20477-t20478-t20479-t20480-t20481-t20482+t20483
      t20486 = t11626/96
      t20487 = 19.D0/2016.D0*t11628
      t20488 = t11634/2016
      t20489 = t11636/5376
      t20490 = t11639/1344
      t20492 = t11643/1008
      t20493 = t11645/252
      t20494 = 5.D0/4032.D0*t11651
      t20495 = 5.D0/1008.D0*t11653
      t20496 = t11655/1008
      t20497 = 23.D0/4032.D0*t11660
      t20498 = t11663/16128
      t20499 = t11665/4032
      t20500 = t11667/252
      t20501 = -t20486+t20487-t20488-t20489+t20490+t11641/96-t20492+t204
     #93+t20494-t20495-t20496+t20497-t20498+t20499+t20500
      t20502 = 17.D0/16128.D0*t11669
      t20503 = t11672/126
      t20504 = t11675/4032
      t20505 = t11677/252
      t20506 = 5.D0/2016.D0*t11682
      t20507 = 17.D0/16128.D0*t11684
      t20508 = t11692/672
      t20511 = t11700/5376
      t20512 = 5.D0/32256.D0*t11702
      t20517 = t20502-t20503+t20504+t20505+t20506-t20507+t20508-t11694/9
     #6+t11696/224+t20511+t20512-t11704/448-t11706/672-t11713/192+t11716
     #/224
      t20526 = t11726/3584
      t20527 = 5.D0/16128.D0*t11728
      t20528 = 13.D0/64512.D0*t11732
      t20529 = 5.D0/32256.D0*t11735
      t20530 = 59.D0/8064.D0*t11737
      t20531 = t11739/1792
      t20532 = t11742/896
      t20534 = t11747/1792
      t20536 = t11718/224+t11720/448-53.D0/5376.D0*t11722+13.D0/10752.D0
     #*t11724+t20526-t20527+t20528-t20529-t20530+t20531+t20532+t11744/22
     #4-t20534-t11750/384
      t20538 = t11755/1792
      t20541 = 5.D0/1344.D0*t11771
      t20542 = 5.D0/64512.D0*t11774
      t20543 = t11776/21504
      t20546 = 5.D0/32256.D0*t11782
      t20547 = t11784/3584
      t20548 = t11787/896
      t20549 = 5.D0/16128.D0*t11791
      t20550 = 59.D0/16128.D0*t11793
      t20551 = 5.D0/2016.D0*t11795
      t20552 = 53.D0/8064.D0*t11753+t20538-5.D0/504.D0*t11757+67.D0/1612
     #8.D0*t11767-t20541-t20542+t20543-109.D0/16128.D0*t11778-11.D0/448.
     #D0*t11780-t20546-t20547+t20548+t20549-t20550-t20551
      t20554 = t11797/43008
      t20557 = 5.D0/504.D0*t11803
      t20560 = 5.D0/2688.D0*t11813
      t20561 = 5.D0/129024.D0*t11815
      t20562 = 5.D0/1008.D0*t11817
      t20563 = 5.D0/64512.D0*t11819
      t20565 = t11824/3584
      t20567 = 5.D0/64512.D0*t11829
      t20568 = 59.D0/16128.D0*t11831
      t20569 = -t20554+11.D0/896.D0*t11799-t11801/96-t20557-t11809/8064-
     #127.D0/32256.D0*t11811+t20560+t20561+t20562+t20563-67.D0/32256.D0*
     #t11822-t20565-17.D0/2016.D0*t11827+t20567+t20568
      t20573 = t11839/84
      t20574 = t11843/448
      t20575 = t11845/84
      t20576 = 5.D0/2016.D0*t11847
      t20577 = t11851/126
      t20578 = t11854/63
      t20579 = 5.D0/336.D0*t11856
      t20580 = 5.D0/2016.D0*t11858
      t20581 = t11860/504
      t20582 = t11862/126
      t20583 = t11867/63
      t20584 = t11871/126
      t20585 = t11833/8064-11.D0/896.D0*t11835+t11837/96-t20573-t20574+t
     #20575+t20576-t20577-t20578-t20579-t20580-t20581-t20582-t20583+t205
     #84
      t20588 = 5.D0/336.D0*t11874
      t20590 = 5.D0/1008.D0*t11878
      t20591 = 5.D0/1008.D0*t11880
      t20592 = 5.D0/1008.D0*t11882
      t20593 = 5.D0/8064.D0*t11884
      t20594 = t11887/63
      t20595 = 17.D0/16128.D0*t11889
      t20596 = 17.D0/16128.D0*t11891
      t20597 = t11894/2016
      t20598 = t11897/672
      t20599 = 17.D0/16128.D0*t11899
      t20600 = t11905/504
      t20601 = 5.D0/2016.D0*t11907
      t20602 = t20588-67.D0/4032.D0*t11876+t20590-t20591+t20592-t20593+t
     #20594-t20595-t20596-t20597+t20598+t20599-t20600-t20601
      t20605 = t11913/252
      t20606 = 5.D0/16128.D0*t11915
      t20607 = 5.D0/336.D0*t11917
      t20608 = 5.D0/1008.D0*t11919
      t20609 = 17.D0/8064.D0*t11921
      t20612 = 5.D0/32256.D0*t11927
      t20614 = 5.D0/4032.D0*t11931
      t20616 = 13.D0/8064.D0*t11938
      t20617 = 5.D0/16128.D0*t11940
      t20618 = 127.D0/8064.D0*t11909+13.D0/504.D0*t11911+t20605+t20606+t
     #20607-t20608+t20609-127.D0/16128.D0*t11923-13.D0/1008.D0*t11925-t2
     #0612+109.D0/32256.D0*t11929+t20614-17.D0/1008.D0*t11934+t20616+t20
     #617
      t20622 = 5.D0/672.D0*t11946
      t20625 = 5.D0/1344.D0*t11954
      t20626 = t11956/384
      t20627 = 5.D0/4032.D0*t11958
      t20631 = 59.D0/8064.D0*t11967
      t20632 = 5.D0/336.D0*t11969
      t20634 = 5.D0/1008.D0*t11974
      t20635 = 5.D0/504.D0*t11942+3.D0/112.D0*t11944-t20622-t11948/96-t1
     #1952/448-t20625+t20626+t20627+127.D0/16128.D0*t11960-67.D0/8064.D0
     #*t11963+t11965/4032-t20631+t20632-109.D0/8064.D0*t11971-t20634
      t20636 = t11976/448
      t20638 = t11981/448
      t20640 = t11985/896
      t20641 = t11988/448
      t20651 = t20636-5.D0/252.D0*t11978-t20638+5.D0/1008.D0*t11983-t206
     #40+t20641-13.D0/2016.D0*t11990-t11992/8064+5.D0/126.D0*t11994-53.D
     #0/2016.D0*t11996+t11998/8064-11.D0/896.D0*t12000-t12002/448+3.D0/2
     #24.D0*t12004+t12006/96
      t20656 = 5.D0/8064.D0*t12014
      t20659 = 5.D0/1008.D0*t12020
      t20660 = 5.D0/1008.D0*t12022
      t20661 = 13.D0/16128.D0*t12024
      t20662 = 5.D0/8064.D0*t12026
      t20663 = 59.D0/8064.D0*t12028
      t20664 = 3.D0/224.D0*t12030
      t20666 = t12034/896
      t20667 = 59.D0/32256.D0*t12036
      t20668 = 5.D0/16128.D0*t12038
      t20669 = -17.D0/2016.D0*t12012+t20656+t12016/1536+127.D0/8064.D0*t
     #12018+t20659-t20660+t20661-t20662+t20663-t20664+5.D0/252.D0*t12032
     #+t20666+t20667+t20668
      t20671 = t12043/56
      t20673 = 3.D0/224.D0*t12047
      t20675 = 5.D0/4032.D0*t12053
      t20677 = 59.D0/64512.D0*t12057
      t20678 = t12062/56
      t20680 = t12066/3584
      t20681 = t12068/56
      t20682 = 5.D0/504.D0*t12070
      t20683 = t12072/1792
      t20684 = t12074/43008
      t20685 = t12040/192+t20671-53.D0/4032.D0*t12045+t20673+11.D0/3584.
     #D0*t12051+t20675+191.D0/32256.D0*t12055-t20677-t20678-109.D0/8064.
     #D0*t12064+t20680-t20681+t20682-t20683-t20684
      t20687 = 5.D0/504.D0*t12078
      t20688 = 5.D0/504.D0*t12080
      t20689 = t12082/3584
      t20690 = 3.D0/448.D0*t12084
      t20692 = t12088/224
      t20693 = 5.D0/64512.D0*t12090
      t20696 = 3.D0/448.D0*t12096
      t20697 = 5.D0/2688.D0*t12098
      t20698 = 5.D0/2016.D0*t12100
      t20699 = 13.D0/129024.D0*t12102
      t20701 = t12108/224
      t20702 = -t20687+t20688-t20689+t20690-5.D0/504.D0*t12086-t20692-t2
     #0693-t12092/384+53.D0/8064.D0*t12094-t20696+t20697-t20698-t20699-5
     #.D0/126.D0*t12106-t20701
      t20703 = t12111/896
      t20704 = 23.D0/1008.D0*t12113
      t20707 = 41.D0/2688.D0*t12120
      t20708 = 41.D0/2688.D0*t12122
      t20709 = 23.D0/1008.D0*t12124
      t20710 = 5.D0/1008.D0*t12126
      t20713 = 19.D0/4032.D0*t12132
      t20714 = t12134/96
      t20715 = 19.D0/4032.D0*t12136
      t20718 = -t20703+t20704+5.D0/126.D0*t12115+5.D0/126.D0*t12117-t207
     #07+t20708-t20709-t20710-53.D0/2016.D0*t12128+53.D0/2016.D0*t12130-
     #t20713+t20714+t20715-53.D0/2016.D0*t12138+t12140/48
      t20722 = 5.D0/504.D0*t12147
      t20723 = 5.D0/504.D0*t12149
      t20724 = t12151/2688
      t20726 = 13.D0/16128.D0*t12155
      t20727 = 5.D0/252.D0*t12157
      t20728 = 5.D0/2016.D0*t12159
      t20729 = 5.D0/1008.D0*t12161
      t20730 = 5.D0/252.D0*t12163
      t20733 = t12170/48
      t20734 = 5.D0/504.D0*t12172
      t20735 = 5.D0/504.D0*t12174
      t20736 = t12142/4032-t20722+t20723-t20724-t12153/4032-t20726+t2072
     #7-t20728-t20729-t20730+11.D0/448.D0*t12165-t12167/48-t20733-t20734
     #+t20735
      t20737 = 5.D0/504.D0*t12176
      t20738 = t12179/448
      t20739 = t12181/224
      t20743 = 5.D0/4032.D0*t12190
      t20746 = t12196/7168
      t20747 = t12199/14336
      t20748 = t12202/14336
      t20749 = t12204/7168
      t20751 = 5.D0/8064.D0*t12208
      t20752 = t20737-t20738-t20739+191.D0/4032.D0*t12183-t12185/32256-t
     #12188/112+t20743+17.D0/1008.D0*t12192-17.D0/1008.D0*t12194-t20746+
     #t20747+t20748-t20749+t12206/112+t20751
      t20754 = 5.D0/64512.D0*t12210
      t20755 = t12216/7168
      t20756 = t12218/3584
      t20757 = t12221/192
      t20758 = t12225/896
      t20759 = t12227/1792
      t20761 = 5.D0/16128.D0*t12231
      t20762 = t12233/3584
      t20763 = 5.D0/16128.D0*t12235
      t20764 = t12239/3584
      t20765 = t12242/96
      t20766 = 5.D0/2016.D0*t12244
      t20767 = 5.D0/16128.D0*t12246
      t20768 = 5.D0/8064.D0*t12248
      t20769 = t20754-t20755-t20756-t20757-t20758-t20759-109.D0/4032.D0*
     #t12229+t20761-t20762-t20763-t20764+t20765+t20766-t20767+t20768
      t20770 = 5.D0/8064.D0*t12251
      t20771 = t12253/96
      t20772 = 5.D0/4032.D0*t12255
      t20773 = t12257/192
      t20774 = 5.D0/16128.D0*t12259
      t20775 = t12264/448
      t20779 = 5.D0/672.D0*t12273
      t20780 = 59.D0/16128.D0*t12275
      t20782 = t12279/768
      t20783 = t12281/1792
      t20784 = t12283/7168
      t20785 = -t20770-t20771-t20772+t20773+t20774+t20775-61.D0/2016.D0*
     #t12266+67.D0/8064.D0*t12269+t12271/896-t20779+t20780-t12277/8064-t
     #20782+t20783-t20784
      wI120pc1 = t20769+t19847+t20434+t20718+t19813+t20702+t20164+t20285
     #+t20585+t20685+t20148+t20268+t20131+t20115+t20097+t20252+t20232+t2
     #0215+t20569+t20552+t20536+t20517+t19796+t19963+t19947+t20501+t2040
     #1+t19930+t19914+t20383+t19896+t19780+t20360+t20340+t19763+t20318+t
     #20301+t19714+t19747+t19880+t20081+t20064+t20048+t20031+t20015+t199
     #98+t20752+t20484+t19863+t19982+t20785+t20669+t20468+t20450+t19730+
     #t20651+t20635+t20417+t20618+t20736+t20198+t20602+t20182+t19829
      t20791 = t12294/12
      t20792 = 11.D0/72.D0*t12298
      t20793 = t12302/72
      t20794 = t12306/1152
      t20795 = t12311/12
      t20796 = 11.D0/576.D0*t12314
      t20797 = t12317/576
      t20798 = t12322/24
      t20799 = 13.D0/288.D0*t12326
      t20800 = 13.D0/144.D0*t12329
      t20801 = 13.D0/72.D0*t12332
      t20802 = 11.D0/288.D0*t12335
      t20803 = t20791+t20792+t20793-t20794-t20795-t20796-t20797-t20798-t
     #20799-t20800-t20801-t20802
      t20804 = 11.D0/576.D0*t12339
      t20805 = 13.D0/36.D0*t12343
      t20806 = 11.D0/72.D0*t12348
      t20807 = 13.D0/144.D0*t12350
      t20808 = 11.D0/288.D0*t12354
      t20809 = 11.D0/36.D0*t12357
      t20810 = t12360/48
      t20811 = 13.D0/18.D0*t12362
      t20812 = t12365/24
      t20813 = 13.D0/288.D0*t12368
      t20814 = 11.D0/144.D0*t12371
      t20815 = t12375/6
      t20816 = -t20804+t20805+t20806-t20807-t20808+t20809-t20810+t20811-
     #t20812-t20813-t20814+t20815
      t20818 = t12380/144
      t20819 = t12382/288
      t20820 = 13.D0/72.D0*t12385
      t20821 = 11.D0/72.D0*t12389
      t20822 = t12392/288
      t20823 = 11.D0/144.D0*t12394
      t20824 = 13.D0/36.D0*t12397
      t20825 = t12401/6
      t20826 = t12405/144
      t20827 = t12408/576
      t20828 = t12411/576
      t20829 = 11.D0/144.D0*t12414
      t20830 = t20818+t20819+t20820+t20821-t20822+t20823+t20824+t20825+t
     #20826-t20827-t20828+t20829
      t20831 = t12418/48
      t20832 = t12420/12
      t20833 = 13.D0/72.D0*t12423
      t20834 = t12425/6
      t20835 = 13.D0/144.D0*t12427
      t20836 = t12429/3
      t20837 = t12431/144
      t20838 = t12433/24
      t20839 = 11.D0/288.D0*t12435
      t20840 = t12439/288
      t20841 = t12443/1152
      t20842 = 13.D0/36.D0*t12446
      t20843 = -t20831+t20832+t20833+t20834-t20835+t20836+t20837-t20838-
     #t20839+t20840-t20841+t20842
      wI121pc1 = t20803+t20816+t20830+t20843
      t20863 = t12461/360-t12469/180-t12473/180-t12480/36-t12483/360-t12
     #487/18-t12491/360-t12496/360-t12499/36-t12504/90+t12512/180+t12515
     #/90+t12517/90+t12521/90-t12526/90+t12532/90+t12534/9-t12538/45
      t20882 = -2.D0/45.D0*t12544-t12547/90+t12552/9+t12554/18-t12558/90
     #+t12561/90-t12565/45+t12567/180+t12569/180-t12573/45-t12577/90-2.D
     #0/45.D0*t12581-t12584/180-t12586/180+t12590/45-t12593/9-t12596/45-
     #t12600/90
      t20902 = t12606/45+t12610/90-t12613/18-t12616/9-t12619/180-t12621/
     #18-t12624/45+t12629/45+t12632/45-t12634/180-t12636/90+t12638/90-t1
     #2643/90-t12645/180+t12648/90-t12652/90-t12655/90-t12659/180
      t20921 = -t12662/18-t12665/45-t12667/90-t12669/180+t12671/180-t126
     #74/90-t12676/360+t12678/90-t12680/90+t12684/90-t12686/180+t12688/1
     #8+t12690/180+t12692/90+t12694/180-t12696/90+t12698/45-t12700/45
      t20942 = t12705/90+2.D0/9.D0*t12707+t12709/9+t12711/90+t12715/90-t
     #12717/45+t12722/180-t12724/45-t12726/45+t12728/90-t12730/90+t12732
     #/45+2.D0/45.D0*t12734-2.D0/9.D0*t12736+t12738/45+t12741/180-t12743
     #/9+2.D0/45.D0*t12745
      t20961 = t12748/45+t12750/45+t12752/45-t12754/180+t12759/180+t1276
     #1/45+t12764/180-t12766/180+t12769/18-t12771/180+t12774/180-t12776/
     #90+t12778/180+t12780/360-t12782/90+t12785/180+t12788/36+t12790/180
      t20981 = t12797/180-t12801/180-t12804/18-t12807/180-t12809/9-t1281
     #1/90+t12813/90+t12815/90-t12818/360+t12821/180+t12825/36+t12829/36
     #0+t12832/180+t12834/90-t12839/180+t12842/180-t12846/180+t12848/9
      t21000 = t12851/18-t12853/90-t12855/90+t12857/360+t12859/90+t12861
     #/90+t12864/18+t12866/90+t12868/360-t12870/360-t12872/180-t12875/90
     #-t12878/45-t12881/90-t12884/180+t12887/90+t12890/180+t12893/360
      wI122pc1 = t20863+t20882+t20902+t20921+t20942+t20961+t20981+t21000
      wI123pc1 = wV32
      t21008 = 2.D0/27.D0*t12917
      t21013 = 17.D0/3456.D0*t12929
      t21015 = 17.D0/1728.D0*t12933
      t21019 = t12945/27
      t21020 = t12948/27
      t21021 = -7.D0/3456.D0*t12900-7.D0/864.D0*t12904+13.D0/432.D0*t129
     #07-t12909/108-t12913/54-t21008-7.D0/1728.D0*t12920+13.D0/864.D0*t1
     #2922-t12925/108-t12927/216-t21013+13.D0/864.D0*t12931-t21015+13.D0
     #/1728.D0*t12935+t12938/27-5.D0/432.D0*t12942-t21019-t21020
      t21028 = 4.D0/27.D0*t12971
      t21033 = 2.D0/27.D0*t12983
      t21040 = 53.D0/1728.D0*t12951-t12954/54+t12956/54+53.D0/864.D0*t12
     #960+2.D0/27.D0*t12965+53.D0/432.D0*t12967-t21028+4.D0/27.D0*t12973
     #-7.D0/3456.D0*t12976-5.D0/432.D0*t12978+13.D0/864.D0*t12981-t21033
     #-t12986/54+t12989/27-2.D0/27.D0*t12992-t12995/27-5.D0/216.D0*t1299
     #7+53.D0/864.D0*t12999
      t21048 = 17.D0/864.D0*t13015
      t21054 = 17.D0/1728.D0*t13030
      t21057 = 17.D0/3456.D0*t13038
      t21060 = 2.D0/27.D0*t13003-5.D0/432.D0*t13005-t13007/27-t13009/108
     #-t13011/27-7.D0/1728.D0*t13013-t21048+53.D0/864.D0*t13018-t13020/5
     #4-t13023/108+t13026/27+53.D0/1728.D0*t13028-t21054-t13032/216-t130
     #35/108-t21057-5.D0/216.D0*t13041+13.D0/1728.D0*t13043
      t21070 = 2.D0/27.D0*t13069
      t21071 = 17.D0/1728.D0*t13071
      t21079 = -5.D0/864.D0*t13046-5.D0/108.D0*t13048+2.D0/27.D0*t13050-
     #5.D0/216.D0*t13052-t13056/54+t13058/27-t13061/27-t13063/54-5.D0/43
     #2.D0*t13065-t21070-t21071+2.D0/27.D0*t13073-5.D0/864.D0*t13075-5.D
     #0/216.D0*t13077-7.D0/1728.D0*t13079+t13081/27+t13083/54-5.D0/432.D
     #0*t13086
      wI124pc1 = t21021+t21040+t21060+t21079
      t21093 = t13103/48+t13106/24+t13109/48+t13091/96+t13095/48+t13099/
     #96-t13113/48-t13116/24-t13119/48-7.D0/96.D0*t13123-7.D0/48.D0*t131
     #27-7.D0/96.D0*t13131
      t21106 = -t13136/24-t13139/12-t13142/24-t13165/12-t13168/6-t13171/
     #12-7.D0/48.D0*t13146-7.D0/24.D0*t13149-7.D0/48.D0*t13152-t13156/24
     #-t13159/12-t13162/24
      wI125pc1 = t21093+t21106
      t21129 = 31.D0/864.D0*t13207+17.D0/432.D0*t13210-5.D0/216.D0*t1321
     #3-5.D0/216.D0*t13216-5.D0/216.D0*t13219+31.D0/432.D0*t13222-37.D0/
     #3456.D0*t13226-37.D0/1728.D0*t13229-37.D0/432.D0*t13232-37.D0/864.
     #D0*t13235+31.D0/864.D0*t13238
      t21153 = -5.D0/108.D0*t13273-5.D0/54.D0*t13276-5.D0/108.D0*t13279-
     #5.D0/54.D0*t13282-5.D0/216.D0*t13286-37.D0/864.D0*t13289-37.D0/172
     #8.D0*t13292+17.D0/864.D0*t13295-37.D0/3456.D0*t13298+17.D0/216.D0*
     #t13301-5.D0/216.D0*t13304
      t21178 = 17.D0/432.D0*t13333+17.D0/216.D0*t13335+17.D0/432.D0*t133
     #39+7.D0/216.D0*t13341-5.D0/54.D0*t13344-5.D0/216.D0*t13346-5.D0/43
     #2.D0*t13349+17.D0/432.D0*t13352+29.D0/216.D0*t13355-5.D0/216.D0*t1
     #3358-37.D0/864.D0*t13361
      t21202 = -37.D0/3456.D0*t13398+17.D0/864.D0*t13401+31.D0/1728.D0*t
     #13404+29.D0/216.D0*t13407+31.D0/864.D0*t13410-37.D0/1728.D0*t13413
     #+31.D0/1728.D0*t13417-5.D0/108.D0*t13420-37.D0/864.D0*t13424-37.D0
     #/1728.D0*t13428-37.D0/1728.D0*t13432
      wI126pc1 = t21202+t21178+t21153+t21129-37.D0/1728.D0*t13330+29.D0/
     #108.D0*t13324-5.D0/108.D0*t13327+29.D0/216.D0*t13322-5.D0/216.D0*t
     #13318-5.D0/108.D0*t13320-5.D0/432.D0*t13316+29.D0/432.D0*t13309+7.
     #D0/432.D0*t13311-37.D0/1728.D0*t13314-5.D0/108.D0*t13269-5.D0/54.D
     #0*t13266-5.D0/108.D0*t13260-37.D0/864.D0*t13263+17.D0/216.D0*t1324
     #9-37.D0/432.D0*t13252-5.D0/54.D0*t13255+7.D0/216.D0*t13257+17.D0/1
     #08.D0*t13246-37.D0/864.D0*t13243-37.D0/864.D0*t13198-5.D0/108.D0*t
     #13201+7.D0/108.D0*t13203-5.D0/27.D0*t13194-5.D0/108.D0*t13190-37.D
     #0/3456.D0*t13181+17.D0/216.D0*t13184+17.D0/432.D0*t13187-37.D0/172
     #8.D0*t13178+29.D0/432.D0*t13175-37.D0/1728.D0*t13394-5.D0/108.D0*t
     #13391-37.D0/864.D0*t13374+7.D0/216.D0*t13376-5.D0/108.D0*t13380-37
     #.D0/1728.D0*t13383+7.D0/432.D0*t13385-5.D0/216.D0*t13388-5.D0/54.D
     #0*t13371-5.D0/108.D0*t13367
      wI127pc1 = t13464/24+t13468/12+t13472/24+t13452/12+t13456/6+t13460
     #/12+t13440/12+t13444/6+t13448/12+t13478/6+t13481/3+t13484/6
      t21228 = -t13488/477-4.D0/477.D0*t13491-2.D0/159.D0*t13494+2.D0/26
     #5.D0*t13497-52.D0/2385.D0*t13500-31.D0/9540.D0*t13503-109.D0/4770.
     #D0*t13506-7.D0/1590.D0*t13508+19.D0/1590.D0*t13512-7.D0/1590.D0*t1
     #3514+37.D0/4770.D0*t13517
      t21241 = -4.D0/265.D0*t13522+t13526/265+3.D0/265.D0*t13529+t13532/
     #530+6.D0/265.D0*t13534-8.D0/265.D0*t13536+t13538/265-t13542/3180-t
     #13545/6360+8.D0/477.D0*t13549-26.D0/795.D0*t13552+41.D0/1908.D0*t1
     #3556
      t21254 = -t13560/477+41.D0/3180.D0*t13562-t13566/318+101.D0/4770.D
     #0*t13569-t13572/477-t13575/159+8.D0/795.D0*t13578-t13582/1590-4.D0
     #/477.D0*t13585+7.D0/636.D0*t13587+37.D0/4770.D0*t13590
      t21267 = 4.D0/477.D0*t13593-37.D0/9540.D0*t13596-2.D0/477.D0*t1359
     #8-13.D0/4770.D0*t13601+17.D0/9540.D0*t13604+2.D0/477.D0*t13606-13.
     #D0/2385.D0*t13609-13.D0/4770.D0*t13611-32.D0/2385.D0*t13615+t13618
     #/477+4.D0/477.D0*t13620+16.D0/2385.D0*t13623
      t21281 = 32.D0/2385.D0*t13628+17.D0/4770.D0*t13631+2.D0/477.D0*t13
     #633+8.D0/2385.D0*t13636-t13638/1590-t13640/1590-t13642/3180-17.D0/
     #2385.D0*t13644-27.D0/1060.D0*t13647+8.D0/2385.D0*t13649-16.D0/2385
     #.D0*t13651
      t21294 = 13.D0/1060.D0*t13656-t13658/159-2.D0/159.D0*t13660+101.D0
     #/9540.D0*t13664+t13667/2385+101.D0/2385.D0*t13669-t13672/318-t1367
     #5/159-317.D0/19080.D0*t13679+37.D0/9540.D0*t13681-2.D0/477.D0*t136
     #84-37.D0/9540.D0*t13686
      t21308 = 37.D0/4770.D0*t13691+169.D0/9540.D0*t13694+37.D0/9540.D0*
     #t13696+101.D0/4770.D0*t13698-4.D0/477.D0*t13701+4.D0/795.D0*t13704
     #-13.D0/795.D0*t13707+34.D0/2385.D0*t13710+2.D0/477.D0*t13713-37.D0
     #/9540.D0*t13715-34.D0/2385.D0*t13717+8.D0/795.D0*t13719
      t21321 = 16.D0/795.D0*t13722-17.D0/2385.D0*t13726-17.D0/4770.D0*t1
     #3730-32.D0/2385.D0*t13733+4.D0/477.D0*t13735+8.D0/477.D0*t13737-t1
     #3740/1590-26.D0/795.D0*t13743-4.D0/477.D0*t13745+41.D0/795.D0*t137
     #47+41.D0/1590.D0*t13749+41.D0/954.D0*t13751
      t21336 = -17.D0/2385.D0*t13757+119.D0/4770.D0*t13759-27.D0/1060.D0
     #*t13763-27.D0/1060.D0*t13765+7.D0/636.D0*t13767+41.D0/1908.D0*t137
     #69+41.D0/3180.D0*t13771-27.D0/530.D0*t13773-64.D0/2385.D0*t13775+4
     #1.D0/1590.D0*t13777-23.D0/477.D0*t13779
      t21349 = -2.D0/477.D0*t13784-23.D0/954.D0*t13787+119.D0/4770.D0*t1
     #3790-27.D0/2120.D0*t13793+119.D0/2385.D0*t13795-34.D0/2385.D0*t137
     #97-46.D0/477.D0*t13799+119.D0/9540.D0*t13802+37.D0/4770.D0*t13804-
     #5.D0/477.D0*t13806+t13808/530-t13811/3180
      t21362 = -64.D0/2385.D0*t13815-32.D0/2385.D0*t13817+37.D0/9540.D0*
     #t13819+37.D0/9540.D0*t13821+37.D0/4770.D0*t13823-t13825/3180+8.D0/
     #477.D0*t13827+t13829/1060+t13831/1060+4.D0/477.D0*t13833-2.D0/477.
     #D0*t13835
      t21375 = t13838/530-5.D0/477.D0*t13841-5.D0/954.D0*t13843+37.D0/19
     #080.D0*t13846+37.D0/4770.D0*t13849-34.D0/2385.D0*t13851-68.D0/2385
     #.D0*t13853-2.D0/477.D0*t13855-37.D0/4770.D0*t13857-4.D0/477.D0*t13
     #859-34.D0/2385.D0*t13861-34.D0/2385.D0*t13863
      t21389 = -17.D0/2385.D0*t13868-23.D0/477.D0*t13870+4.D0/477.D0*t13
     #873-317.D0/38160.D0*t13877+169.D0/19080.D0*t13879-31.D0/9540.D0*t1
     #3881-52.D0/2385.D0*t13883+13.D0/530.D0*t13885-17.D0/2385.D0*t13887
     #-4.D0/477.D0*t13889-17.D0/4770.D0*t13891
      t21402 = 17.D0/2385.D0*t13894+34.D0/2385.D0*t13896+t13899/265+3.D0
     #/530.D0*t13902+19.D0/3180.D0*t13905-17.D0/3180.D0*t13908+37.D0/954
     #0.D0*t13911-11.D0/795.D0*t13913-17.D0/2385.D0*t13915-4.D0/477.D0*t
     #13917-68.D0/2385.D0*t13919-2.D0/265.D0*t13922
      t21416 = 37.D0/9540.D0*t13926+19.D0/1590.D0*t13928-17.D0/1590.D0*t
     #13930+8.D0/477.D0*t13932+7.D0/159.D0*t13934-52.D0/795.D0*t13937-t1
     #3939/795-23.D0/954.D0*t13941-2.D0/477.D0*t13943-17.D0/795.D0*t1394
     #5+t13948/530-4.D0/477.D0*t13950
      t21429 = 7.D0/318.D0*t13953-26.D0/795.D0*t13956-t13958/1590-8.D0/4
     #77.D0*t13960-128.D0/2385.D0*t13962+7.D0/318.D0*t13964-37.D0/4770.D
     #0*t13967-2.D0/477.D0*t13969-8.D0/477.D0*t13971-37.D0/2385.D0*t1397
     #3+68.D0/2385.D0*t13976-2.D0/477.D0*t13978
      t21445 = -34.D0/2385.D0*t13985-37.D0/9540.D0*t13987+169.D0/4770.D0
     #*t13989+37.D0/4770.D0*t13991+37.D0/2385.D0*t13993-37.D0/4770.D0*t1
     #3995-4.D0/477.D0*t13997-64.D0/2385.D0*t13999+32.D0/2385.D0*t14003+
     #2.D0/477.D0*t14006+4.D0/2385.D0*t14008
      t21458 = 16.D0/2385.D0*t14011+41.D0/954.D0*t14013-32.D0/2385.D0*t1
     #4015+37.D0/9540.D0*t14017-68.D0/2385.D0*t14020+37.D0/9540.D0*t1402
     #2-13.D0/2385.D0*t14024-32.D0/2385.D0*t14027-13.D0/2385.D0*t14029-1
     #7.D0/1590.D0*t14031+41.D0/1590.D0*t14033-17.D0/2385.D0*t14036
      t21471 = -4.D0/477.D0*t14040+17.D0/2385.D0*t14042+37.D0/19080.D0*t
     #14045-17.D0/2385.D0*t14047+8.D0/795.D0*t14049+37.D0/4770.D0*t14051
     #+4.D0/2385.D0*t14053+169.D0/9540.D0*t14055+17.D0/4770.D0*t14057+t1
     #4059/265+8.D0/2385.D0*t14061
      t21484 = 32.D0/2385.D0*t14064-16.D0/2385.D0*t14066+2.D0/477.D0*t14
     #070+t14072/477+64.D0/2385.D0*t14074+4.D0/2385.D0*t14076+16.D0/2385
     #.D0*t14078+t14080/265+41.D0/477.D0*t14082+t14085/530+4.D0/477.D0*t
     #14087+17.D0/9540.D0*t14089
      t21498 = t14095/530-t14097/265+17.D0/4770.D0*t14099+41.D0/954.D0*t
     #14101-64.D0/2385.D0*t14103-26.D0/2385.D0*t14105+t14108/1060-t14112
     #/530+t14115/1060+t14118/1060+17.D0/2385.D0*t14120
      t21511 = -t14123/1060+t14126/1060+37.D0/2385.D0*t14128+37.D0/4770.
     #D0*t14130-2.D0/265.D0*t14132+t14134/2120-10.D0/477.D0*t14136-t1413
     #9/1060+t14141/2120-5.D0/477.D0*t14143+t14145/530-t14147/265
      t21525 = -t14151/530+37.D0/9540.D0*t14153-2.D0/477.D0*t14155+t1415
     #7/530-t14159/265-t14161/530+t14163/2120+t14165/2120-t14167/530-5.D
     #0/954.D0*t14169+t14171/530+t14173/530
      t21538 = -2.D0/477.D0*t14176+37.D0/19080.D0*t14178-34.D0/2385.D0*t
     #14181+2.D0/2385.D0*t14183+37.D0/38160.D0*t14186+37.D0/9540.D0*t141
     #89+169.D0/9540.D0*t14191+169.D0/19080.D0*t14193+37.D0/19080.D0*t14
     #196-109.D0/9540.D0*t14198-31.D0/19080.D0*t14200-26.D0/2385.D0*t142
     #02
      t21553 = -t14208/159+119.D0/4770.D0*t14210+37.D0/19080.D0*t14213+1
     #3.D0/265.D0*t14216-4.D0/159.D0*t14218-104.D0/2385.D0*t14220-31.D0/
     #4770.D0*t14222-109.D0/2385.D0*t14224-23.D0/477.D0*t14226+13.D0/530
     #.D0*t14229-t14231/265
      t21566 = t14234/1060+t14236/1060-t14238/530-17.D0/2385.D0*t14240+t
     #14242/1060+t14244/1060+37.D0/9540.D0*t14246-34.D0/2385.D0*t14248-1
     #1.D0/1590.D0*t14251+3.D0/530.D0*t14253-17.D0/3180.D0*t14255+19.D0/
     #3180.D0*t14257
      t21579 = -t14261/6360-t14263/1590-4.D0/265.D0*t14266-7.D0/3180.D0*
     #t14269-2.D0/265.D0*t14271-11.D0/795.D0*t14274+3.D0/265.D0*t14276-1
     #7.D0/1590.D0*t14278+19.D0/1590.D0*t14280-7.D0/3180.D0*t14282-22.D0
     #/795.D0*t14284
      t21592 = -7.D0/1590.D0*t14288-11.D0/1590.D0*t14290+3.D0/265.D0*t14
     #293-4.D0/265.D0*t14295-11.D0/795.D0*t14297-7.D0/795.D0*t14299+19.D
     #0/795.D0*t14302-34.D0/2385.D0*t14304+119.D0/9540.D0*t14307-27.D0/2
     #120.D0*t14310-t14312/3180-t14314/1590
      t21606 = 37.D0/19080.D0*t14319+2.D0/2385.D0*t14321-52.D0/2385.D0*t
     #14323-13.D0/795.D0*t14325-t14327/795-t14329/3180-2.D0/477.D0*t1433
     #2+4.D0/795.D0*t14334-2.D0/159.D0*t14336+13.D0/530.D0*t14338+2.D0/2
     #385.D0*t14340
      t21619 = -109.D0/4770.D0*t14343+101.D0/4770.D0*t14345+37.D0/4770.D
     #0*t14347-t14349/3180-37.D0/19080.D0*t14351-t14353/477-317.D0/9540.
     #D0*t14355-26.D0/2385.D0*t14357-31.D0/9540.D0*t14359-t14361/159+37.
     #D0/4770.D0*t14363+37.D0/9540.D0*t14365
      t21633 = -317.D0/19080.D0*t14369-31.D0/19080.D0*t14371+101.D0/9540
     #.D0*t14373-109.D0/9540.D0*t14375+13.D0/1060.D0*t14377+t14379/2385-
     #2.D0/159.D0*t14381-109.D0/4770.D0*t14383+7.D0/318.D0*t14385+37.D0/
     #2385.D0*t14387-37.D0/4770.D0*t14389+37.D0/38160.D0*t14391
      t21646 = 37.D0/19080.D0*t14394-317.D0/19080.D0*t14396-317.D0/38160
     #.D0*t14398-4.D0/477.D0*t14400-8.D0/477.D0*t14402-37.D0/9540.D0*t14
     #404+16.D0/477.D0*t14406+34.D0/2385.D0*t14409+37.D0/9540.D0*t14411-
     #2.D0/477.D0*t14413-34.D0/2385.D0*t14415-37.D0/19080.D0*t14417
      wI128pc1 = t21228+t21375+t21389+t21362+t21349+t21336+t21321+t21308
     #+t21294+t21281+t21267+t21241+t21254+t21592+t21579+t21566+t21553+t2
     #1538+t21525+t21511+t21498+t21484+t21471+t21458+t21445+t21429+t2141
     #6+t21402+t21646+t21633+t21619+t21606
      wI129pc1 = wI71pc1
      wI130pc1 = wI73pc1
      wI131pc1 = wV33
      wI132pc1 = wI74pc1
      t21673 = 57.D0/448.D0*t13207-11.D0/224.D0*t13210-5.D0/112.D0*t1321
     #3-5.D0/112.D0*t13216-3.D0/56.D0*t13219+57.D0/224.D0*t13222+13.D0/1
     #792.D0*t13226-13.D0/896.D0*t13229-13.D0/224.D0*t13232-13.D0/448.D0
     #*t13235+57.D0/448.D0*t13238
      t21697 = 5.D0/56.D0*t13273+3.D0/14.D0*t13276+3.D0/28.D0*t13279+5.D
     #0/28.D0*t13282-3.D0/56.D0*t13286-13.D0/448.D0*t13289-13.D0/896.D0*
     #t13292-11.D0/448.D0*t13295+13.D0/1792.D0*t13298-11.D0/112.D0*t1330
     #1-5.D0/112.D0*t13304
      t21722 = 11.D0/224.D0*t13333+11.D0/112.D0*t13335+11.D0/224.D0*t133
     #39-3.D0/16.D0*t13341+5.D0/28.D0*t13344-5.D0/112.D0*t13346+5.D0/224
     #.D0*t13349+11.D0/224.D0*t13352-27.D0/112.D0*t13355-5.D0/112.D0*t13
     #358-13.D0/448.D0*t13361
      t21746 = 13.D0/1792.D0*t13398-11.D0/448.D0*t13401-57.D0/896.D0*t13
     #404+27.D0/112.D0*t13407-57.D0/448.D0*t13410+13.D0/896.D0*t13413-57
     #.D0/896.D0*t13417-5.D0/56.D0*t13420+13.D0/448.D0*t13424+13.D0/896.
     #D0*t13428+13.D0/896.D0*t13432
      wI133pc1 = t21746+t21722-13.D0/896.D0*t13330-27.D0/56.D0*t13324+3.
     #D0/28.D0*t13327-27.D0/112.D0*t13322+5.D0/112.D0*t13318-5.D0/56.D0*
     #t13320+5.D0/224.D0*t13316+27.D0/224.D0*t13309+3.D0/32.D0*t13311-13
     #.D0/896.D0*t13314+5.D0/56.D0*t13269+5.D0/28.D0*t13266+5.D0/56.D0*t
     #13260+13.D0/448.D0*t13263-11.D0/112.D0*t13249+13.D0/224.D0*t13252+
     #5.D0/28.D0*t13255-3.D0/16.D0*t13257-11.D0/56.D0*t13246+13.D0/448.D
     #0*t13243-13.D0/448.D0*t13198+5.D0/56.D0*t13201-3.D0/8.D0*t13203+5.
     #D0/14.D0*t13194-5.D0/56.D0*t13190+13.D0/1792.D0*t13181+11.D0/112.D
     #0*t13184+11.D0/224.D0*t13187+13.D0/896.D0*t13178+27.D0/224.D0*t131
     #75-13.D0/896.D0*t13394-5.D0/56.D0*t13391-13.D0/448.D0*t13374+3.D0/
     #16.D0*t13376-5.D0/56.D0*t13380-13.D0/896.D0*t13383+3.D0/32.D0*t133
     #85-5.D0/112.D0*t13388-5.D0/28.D0*t13371-3.D0/28.D0*t13367+t21697+t
     #21673
      t21760 = -t13488/680+t13491/68-3.D0/136.D0*t13494-9.D0/680.D0*t134
     #97-7.D0/136.D0*t13500+57.D0/5440.D0*t13503-43.D0/680.D0*t13506+21.
     #D0/2720.D0*t13508+27.D0/2720.D0*t13512+21.D0/2720.D0*t13514+13.D0/
     #2720.D0*t13517
      t21773 = -9.D0/340.D0*t13522+9.D0/1360.D0*t13526+27.D0/1360.D0*t13
     #529+9.D0/2720.D0*t13532+27.D0/680.D0*t13534-9.D0/170.D0*t13536-9.D
     #0/1360.D0*t13538+3.D0/5440.D0*t13542-3.D0/10880.D0*t13545-t13549/3
     #4+t13552/40-3.D0/64.D0*t13556
      t21786 = -t13560/272-329.D0/5440.D0*t13562-3.D0/544.D0*t13566+109.
     #D0/2720.D0*t13569-t13572/272+3.D0/272.D0*t13575-2.D0/85.D0*t13578+
     #7.D0/544.D0*t13582+t13585/68-t13587/5440-13.D0/2720.D0*t13590
      t21799 = -t13593/68+13.D0/5440.D0*t13596+t13598/136-7.D0/2720.D0*t
     #13601+t13604/320+3.D0/136.D0*t13606-7.D0/1360.D0*t13609-7.D0/2720.
     #D0*t13611-2.D0/85.D0*t13615-3.D0/272.D0*t13618+3.D0/68.D0*t13620-3
     #.D0/85.D0*t13623
      t21813 = 2.D0/85.D0*t13628-t13631/160+3.D0/136.D0*t13633-3.D0/170.
     #D0*t13636-3.D0/2720.D0*t13638+3.D0/2720.D0*t13640+3.D0/5440.D0*t13
     #642+11.D0/1360.D0*t13644+219.D0/5440.D0*t13647-3.D0/170.D0*t13649+
     #t13651/85
      t21826 = -37.D0/1088.D0*t13656+3.D0/272.D0*t13658+3.D0/136.D0*t136
     #60-109.D0/5440.D0*t13664-t13667/136+109.D0/1360.D0*t13669-3.D0/544
     #.D0*t13672-3.D0/272.D0*t13675+49.D0/2176.D0*t13679-13.D0/5440.D0*t
     #13681+t13684/136+13.D0/5440.D0*t13686
      t21840 = -13.D0/2720.D0*t13691-29.D0/1088.D0*t13694-13.D0/5440.D0*
     #t13696+109.D0/2720.D0*t13698-t13701/68-t13704/85+t13707/80+11.D0/6
     #80.D0*t13710+t13713/136+13.D0/5440.D0*t13715+11.D0/680.D0*t13717+2
     #.D0/85.D0*t13719
      t21853 = 4.D0/85.D0*t13722-11.D0/1360.D0*t13726-11.D0/2720.D0*t137
     #30+161.D0/2720.D0*t13733-t13735/68-t13737/34-7.D0/544.D0*t13740-t1
     #3743/40-t13745/68+329.D0/1360.D0*t13747-329.D0/2720.D0*t13749-3.D0
     #/32.D0*t13751
      t21868 = 11.D0/1360.D0*t13757-7.D0/160.D0*t13759-219.D0/5440.D0*t1
     #3763+219.D0/5440.D0*t13765-t13767/5440-3.D0/64.D0*t13769-329.D0/54
     #40.D0*t13771+219.D0/2720.D0*t13773-161.D0/1360.D0*t13775+329.D0/27
     #20.D0*t13777-151.D0/1360.D0*t13779
      t21881 = t13784/136+151.D0/2720.D0*t13787+7.D0/160.D0*t13790-219.D
     #0/10880.D0*t13793-7.D0/80.D0*t13795+11.D0/680.D0*t13797-151.D0/680
     #.D0*t13799+7.D0/320.D0*t13802-13.D0/2720.D0*t13804-11.D0/1360.D0*t
     #13806+9.D0/2720.D0*t13808-3.D0/5440.D0*t13811
      t21894 = 161.D0/1360.D0*t13815+161.D0/2720.D0*t13817-13.D0/5440.D0
     #*t13819-13.D0/5440.D0*t13821-13.D0/2720.D0*t13823+7.D0/1088.D0*t13
     #825+t13827/34-9.D0/5440.D0*t13829-9.D0/5440.D0*t13831-t13833/68+t1
     #3835/136
      t21907 = 9.D0/2720.D0*t13838+11.D0/1360.D0*t13841-11.D0/2720.D0*t1
     #3843+13.D0/10880.D0*t13846+13.D0/2720.D0*t13849-11.D0/680.D0*t1385
     #1-11.D0/340.D0*t13853+t13855/340-13.D0/2720.D0*t13857+t13859/170+1
     #1.D0/680.D0*t13861-11.D0/680.D0*t13863
      t21921 = -11.D0/1360.D0*t13868+151.D0/1360.D0*t13870+t13873/68-49.
     #D0/4352.D0*t13877+29.D0/2176.D0*t13879+57.D0/5440.D0*t13881-7.D0/1
     #36.D0*t13883+37.D0/544.D0*t13885+11.D0/1360.D0*t13887+t13889/68-11
     #.D0/2720.D0*t13891
      t21934 = -11.D0/1360.D0*t13894-11.D0/680.D0*t13896-9.D0/1360.D0*t1
     #3899-27.D0/2720.D0*t13902-27.D0/5440.D0*t13905+9.D0/1088.D0*t13908
     #+13.D0/5440.D0*t13911+3.D0/1360.D0*t13913-11.D0/1360.D0*t13915-t13
     #917/68+11.D0/340.D0*t13919+9.D0/680.D0*t13922
      t21948 = 13.D0/5440.D0*t13926+27.D0/2720.D0*t13928-9.D0/544.D0*t13
     #930+t13932/34+t13934/1360-t13937/20-7.D0/272.D0*t13939+151.D0/2720
     #.D0*t13941+t13943/136-9.D0/272.D0*t13945+9.D0/2720.D0*t13948-t1395
     #0/68
      t21961 = t13953/2720-t13956/40-7.D0/544.D0*t13958-t13960/34-161.D0
     #/680.D0*t13962+t13964/2720+13.D0/2720.D0*t13967+t13969/136-t13971/
     #34-13.D0/1360.D0*t13973+11.D0/340.D0*t13976+t13978/340
      t21977 = 11.D0/680.D0*t13985+13.D0/5440.D0*t13987-29.D0/544.D0*t13
     #989-13.D0/2720.D0*t13991+13.D0/1360.D0*t13993-13.D0/2720.D0*t13995
     #-t13997/68-161.D0/1360.D0*t13999+2.D0/85.D0*t14003+t14006/136+3.D0
     #/340.D0*t14008
      t21990 = -t14011/85+3.D0/32.D0*t14013-2.D0/85.D0*t14015-13.D0/5440
     #.D0*t14017-11.D0/340.D0*t14020-13.D0/5440.D0*t14022+7.D0/1360.D0*t
     #14024+2.D0/85.D0*t14027+7.D0/1360.D0*t14029-9.D0/544.D0*t14031+329
     #.D0/2720.D0*t14033+11.D0/1360.D0*t14036
      t22003 = t14040/68-t14042/80+13.D0/10880.D0*t14045+11.D0/1360.D0*t
     #14047+2.D0/85.D0*t14049+13.D0/2720.D0*t14051+t14053/34-29.D0/1088.
     #D0*t14055-t14057/160+9.D0/1360.D0*t14059+3.D0/170.D0*t14061
      t22016 = -2.D0/85.D0*t14064+t14066/85-3.D0/136.D0*t14070-3.D0/272.
     #D0*t14072+4.D0/85.D0*t14074+3.D0/340.D0*t14076-t14078/85+9.D0/1360
     #.D0*t14080+3.D0/16.D0*t14082-9.D0/2720.D0*t14085-t14087/68+t14089/
     #320
      t22030 = -9.D0/2720.D0*t14095-9.D0/1360.D0*t14097+t14099/160+3.D0/
     #32.D0*t14101-4.D0/85.D0*t14103+7.D0/680.D0*t14105-9.D0/5440.D0*t14
     #108-9.D0/2720.D0*t14112-9.D0/5440.D0*t14115-9.D0/5440.D0*t14118-11
     #.D0/1360.D0*t14120
      t22043 = -9.D0/5440.D0*t14123-9.D0/5440.D0*t14126+13.D0/1360.D0*t1
     #4128+13.D0/2720.D0*t14130-9.D0/680.D0*t14132+9.D0/10880.D0*t14134+
     #11.D0/680.D0*t14136-9.D0/5440.D0*t14139+9.D0/10880.D0*t14141+11.D0
     #/1360.D0*t14143+9.D0/2720.D0*t14145-9.D0/1360.D0*t14147
      t22057 = 9.D0/2720.D0*t14151-13.D0/5440.D0*t14153-t14155/136+9.D0/
     #2720.D0*t14157+9.D0/1360.D0*t14159+9.D0/2720.D0*t14161+9.D0/10880.
     #D0*t14163+9.D0/10880.D0*t14165+9.D0/2720.D0*t14167-11.D0/2720.D0*t
     #14169-9.D0/2720.D0*t14171-9.D0/2720.D0*t14173
      t22070 = t14176/136-13.D0/10880.D0*t14178-11.D0/680.D0*t14181+t141
     #83/68+13.D0/21760.D0*t14186+13.D0/5440.D0*t14189+29.D0/1088.D0*t14
     #191+29.D0/2176.D0*t14193+13.D0/10880.D0*t14196+43.D0/1360.D0*t1419
     #8-57.D0/10880.D0*t14200+7.D0/272.D0*t14202
      t22085 = 3.D0/272.D0*t14208-7.D0/160.D0*t14210+13.D0/10880.D0*t142
     #13+37.D0/272.D0*t14216-3.D0/68.D0*t14218-7.D0/68.D0*t14220+57.D0/2
     #720.D0*t14222-43.D0/340.D0*t14224-151.D0/1360.D0*t14226+37.D0/544.
     #D0*t14229+9.D0/1360.D0*t14231
      t22098 = -9.D0/5440.D0*t14234-9.D0/5440.D0*t14236+9.D0/2720.D0*t14
     #238+11.D0/1360.D0*t14240+9.D0/5440.D0*t14242+9.D0/5440.D0*t14244+1
     #3.D0/5440.D0*t14246+11.D0/680.D0*t14248-3.D0/2720.D0*t14251-27.D0/
     #2720.D0*t14253+9.D0/1088.D0*t14255-27.D0/5440.D0*t14257
      t22111 = -3.D0/10880.D0*t14261+3.D0/2720.D0*t14263+9.D0/340.D0*t14
     #266-21.D0/5440.D0*t14269+9.D0/680.D0*t14271-3.D0/1360.D0*t14274-27
     #.D0/1360.D0*t14276+9.D0/544.D0*t14278-27.D0/2720.D0*t14280-21.D0/5
     #440.D0*t14282+3.D0/680.D0*t14284
      t22124 = -21.D0/2720.D0*t14288-3.D0/2720.D0*t14290+27.D0/1360.D0*t
     #14293-9.D0/340.D0*t14295+3.D0/1360.D0*t14297+21.D0/1360.D0*t14299+
     #27.D0/1360.D0*t14302-11.D0/680.D0*t14304+7.D0/320.D0*t14307-219.D0
     #/10880.D0*t14310+3.D0/5440.D0*t14312-3.D0/2720.D0*t14314
      t22138 = -13.D0/10880.D0*t14319+t14321/68+7.D0/136.D0*t14323+t1432
     #5/80-3.D0/1360.D0*t14327+3.D0/5440.D0*t14329+t14332/136-t14334/85+
     #3.D0/136.D0*t14336-37.D0/544.D0*t14338-t14340/68
      t22151 = 43.D0/680.D0*t14343-109.D0/2720.D0*t14345-13.D0/2720.D0*t
     #14347+7.D0/1088.D0*t14349-13.D0/10880.D0*t14351-t14353/680+49.D0/1
     #088.D0*t14355+7.D0/272.D0*t14357-57.D0/5440.D0*t14359+3.D0/272.D0*
     #t14361+13.D0/2720.D0*t14363-13.D0/5440.D0*t14365
      t22165 = 49.D0/2176.D0*t14369-57.D0/10880.D0*t14371-109.D0/5440.D0
     #*t14373+43.D0/1360.D0*t14375-37.D0/1088.D0*t14377-t14379/136-3.D0/
     #136.D0*t14381-43.D0/680.D0*t14383-t14385/2720-13.D0/1360.D0*t14387
     #+13.D0/2720.D0*t14389+13.D0/21760.D0*t14391
      t22178 = 13.D0/10880.D0*t14394-49.D0/2176.D0*t14396-49.D0/4352.D0*
     #t14398-t14400/68-t14402/34-13.D0/5440.D0*t14404+t14406/17+11.D0/68
     #0.D0*t14409-13.D0/5440.D0*t14411-t14413/340-11.D0/680.D0*t14415-13
     #.D0/10880.D0*t14417
      wI134pc1 = t22070+t22124+t22057+t22043+t22030+t22016+t22003+t21990
     #+t21977+t21961+t21948+t21934+t21921+t21907+t21894+t21881+t21868+t2
     #1853+t21840+t21826+t21813+t21799+t21786+t22178+t21773+t21760+t2213
     #8+t22165+t22151+t22111+t22098+t22085
      wI135pc1 = wV34
      t22197 = t15954/68-t15956/68+t15958/40-t15962/40-43.D0/2720.D0*t15
     #964+t15966/17+t15968/68-t15970/80-t15973/34+t15976/17-t15978/17+t1
     #5980/34+t15983/80+t15985/34
      t22213 = -4.D0/85.D0*t15989+2.D0/85.D0*t15992+9.D0/1360.D0*t15994+
     #9.D0/1360.D0*t15996-9.D0/2720.D0*t15998+t16000/68-9.D0/2720.D0*t16
     #002+9.D0/1360.D0*t16004-9.D0/680.D0*t16006+9.D0/1360.D0*t16008+9.D
     #0/680.D0*t16011-9.D0/1360.D0*t16014-3.D0/17.D0*t16017+3.D0/34.D0*t
     #16019+t16022/34
      t22229 = -t16026/68-t16028/68+t16030/136+3.D0/17.D0*t16032-3.D0/68
     #.D0*t16034-3.D0/34.D0*t16036+3.D0/136.D0*t16038+t16040/68-t16042/6
     #8-2.D0/17.D0*t16045+t16047/34-t16049/136-t16051/68+t16053/68
      t22245 = 9.D0/170.D0*t16058-9.D0/340.D0*t16061-3.D0/680.D0*t16063+
     #3.D0/1360.D0*t16065+3.D0/2720.D0*t16067+3.D0/1360.D0*t16070-3.D0/2
     #720.D0*t16073+3.D0/680.D0*t16075-3.D0/1360.D0*t16077-3.D0/1360.D0*
     #t16079+2.D0/17.D0*t16081-t16083/68-t16085/68-2.D0/17.D0*t16087+t16
     #089/68
      t22262 = -t16094/68-7.D0/136.D0*t16096-123.D0/1360.D0*t16098-t1610
     #0/80+2.D0/85.D0*t16103-163.D0/2720.D0*t16105+t16107/68+7.D0/136.D0
     #*t16109-9.D0/340.D0*t16111+3.D0/170.D0*t16113-9.D0/170.D0*t16115-3
     #.D0/1360.D0*t16118-3.D0/340.D0*t16120+9.D0/340.D0*t16122
      t22278 = -11.D0/34.D0*t16126-t16128/34-3.D0/2720.D0*t16130+3.D0/27
     #20.D0*t16132+3.D0/1360.D0*t16134+t16136/80-2.D0/85.D0*t16139+4.D0/
     #85.D0*t16141+t16143/68-t16145/34-11.D0/68.D0*t16147+11.D0/34.D0*t1
     #6149+11.D0/68.D0*t16151+t16153/34-9.D0/1360.D0*t16156
      t22295 = -t16160/68+t16162/34-9.D0/1360.D0*t16164-9.D0/1360.D0*t16
     #166-9.D0/2720.D0*t16168-9.D0/2720.D0*t16170+9.D0/2720.D0*t16172+9.
     #D0/2720.D0*t16174-t16176/34+t16178/17-t16180/17+t16182/68-t16184/6
     #8+t16186/34-t16188/68
      t22311 = t16191/68-t16193/34-t16195/17+t16197/68-t16199/68+t16201/
     #34+43.D0/2720.D0*t16203-t16205/34+t16207/17+t16209/68-t16211/68-t1
     #6213/34-t16215/68+t16217/34+t16219/17
      t22329 = t16225/34+t16227/34+t16229/34+t16232/34-t16234/17-t16236/
     #17-t16238/34-3.D0/34.D0*t16240-43.D0/2720.D0*t16242+3.D0/68.D0*t16
     #244+3.D0/136.D0*t16246-t16248/68-t16250/136-3.D0/136.D0*t16252
      t22345 = -3.D0/136.D0*t16255+3.D0/68.D0*t16257+3.D0/68.D0*t16259+1
     #63.D0/2720.D0*t16261-163.D0/1360.D0*t16263+3.D0/34.D0*t16265-3.D0/
     #136.D0*t16267-t16269/34-3.D0/68.D0*t16271-t16273/34+3.D0/136.D0*t1
     #6275+t16277/68-t16279/68+t16281/68+t16283/136
      t22361 = -3.D0/68.D0*t16287-7.D0/68.D0*t16289+3.D0/340.D0*t16291+3
     #.D0/136.D0*t16293+163.D0/1360.D0*t16295-t16297/34+t16299/34+9.D0/2
     #720.D0*t16301-t16303/34-t16305/34+t16307/68+t16309/68-t16311/68+t1
     #6313/68
      t22377 = -2.D0/85.D0*t16316+9.D0/340.D0*t16319-3.D0/340.D0*t16321-
     #3.D0/170.D0*t16324+3.D0/340.D0*t16327+t16329/68+123.D0/1360.D0*t16
     #331+t16333/34-t16335/68-t16337/34-t16339/68+9.D0/2720.D0*t16341+t1
     #6343/68+2.D0/17.D0*t16345-123.D0/1360.D0*t16347
      t22394 = -11.D0/68.D0*t16352+t16354/68-123.D0/680.D0*t16356+t16358
     #/68-t16360/68-4.D0/17.D0*t16362-t16364/34+11.D0/68.D0*t16366+123.D
     #0/680.D0*t16368+163.D0/2720.D0*t16371-3.D0/136.D0*t16373+t16375/34
     #-t16377/68-t16379/34
      t22410 = t16382/34-3.D0/136.D0*t16384-3.D0/136.D0*t16386+t16388/68
     #+t16390/68+t16392/34-t16394/68+3.D0/136.D0*t16396+7.D0/136.D0*t163
     #98+t16400/68-163.D0/2720.D0*t16402+7.D0/68.D0*t16404-7.D0/136.D0*t
     #16406-t16408/68+t16410/34
      t22427 = t16415/34-t16417/34-t16419/68+3.D0/136.D0*t16421-t16423/6
     #8-t16425/68-t16427/34-t16429/34+t16431/68+t16433/34+t16435/68+t164
     #37/34+123.D0/1360.D0*t16439-t16441/68+4.D0/17.D0*t16443
      t22443 = -t16446/68-t16448/34+t16450/68-43.D0/1360.D0*t16452+t1645
     #4/68+43.D0/2720.D0*t16456-t16458/34-t16460/34+43.D0/1360.D0*t16462
     #+t16464/68-t16466/34+t16468/34-t16470/68+t16472/34-t16474/34
      wI136pc1 = t22197+t22213+t22229+t22245+t22262+t22278+t22295+t22311
     #+t22329+t22345+t22361+t22377+t22394+t22410+t22427+t22443
      t22447 = t1072*wC22
      t22448 = t6953*t22447
      t22451 = t1077*wC22
      t22452 = t6949*t22451
      t22455 = t6949*t22447
      t22460 = t1070*wC22
      t22461 = t6949*t22460
      t22466 = t6953*t22451
      t22469 = t6953*t22460
      t22474 = t16490/32
      t22476 = t16494/32
      t22477 = t16498/16
      t22478 = t16504/16
      t22479 = t16507/16
      t22480 = t3999*wC22
      t22481 = t6949*t22480
      t22488 = t6943*t22451
      t22491 = t6943*t22460
      t22494 = t4120*wC22
      t22495 = t6949*t22494
      t22500 = t22476+t22477-t22478+t22479-t16513*t22481/32-t16501*t2248
     #1/32+t16506*t22481/32+t13595*t22488/512-t13595*t22491/512-t16513*t
     #22495/32-t16501*t22452/32
      t22506 = t16518/16
      t22507 = t16522/16
      t22508 = t16524/16
      t22509 = t16527/64
      t22510 = t16529/16
      t22511 = t16538/256
      t22512 = t16544/128
      t22513 = t16548/128
      t22514 = t16552/512
      t22515 = t16501*t22455/32+t16487*t22495/128-t22506+t22507-t22508-t
     #22509+t22510-t22511+t22512-t22513+t22514
      t22520 = t54*wC22
      t22521 = t7137*t22520
      t22524 = t16556/64
      t22525 = t7103*t22520
      t22528 = t16573/64
      t22529 = t16577/16
      t22530 = t16583/128
      t22531 = t16587/64
      t22532 = t16591/64
      t22533 = t16594/32
      t22534 = t16580*t22452/128-t16580*t22455/128-t16570*t22521/128-t22
     #524-t16570*t22525/128+t22528-t22529-t22530+t22531-t22532-t22533
      t22537 = t16598/128
      t22538 = t16606/64
      t22539 = t16608/32
      t22546 = s*wC22
      t22547 = t6987*t22546
      t22550 = t3409*wC22
      t22551 = t6987*t22550
      t22556 = t16612/256
      t22559 = t22537+t22538-t22539-t16480*t22461/32+t16480*t22525/32+t1
     #6480*t22521/16-t16480*t22547/32+t16480*t22551/32+t13739*t22491/128
     #+t22556+t16480*t22469/16
      t22560 = t16621/256
      t22561 = t16625/128
      t22562 = t16629/64
      t22563 = t16632/32
      t22564 = t16636/64
      t22565 = t16639/32
      t22566 = t16647/32
      t22567 = t16649/64
      t22568 = t16652/256
      t22569 = t16659/512
      t22570 = t16667/512
      t22571 = t22560+t22561-t22562-t22563-t22564-t22565+t22566-t22567+t
     #22568-t22569-t22570
      t22573 = t16671/256
      t22574 = t16675/128
      t22575 = t16678/32
      t22576 = t16682/128
      t22577 = t16686/16
      t22578 = t16688/16
      t22579 = t16690/64
      t22580 = t16692/32
      t22581 = t16701/16
      t22582 = t16708/32
      t22583 = t16711/128
      t22584 = -t22573-t22574+t22575+t22576+t22577-t22578-t22579-t22580+
     #t22581-t22582-t22583
      t22585 = t16714/64
      t22586 = t16717/256
      t22587 = t16722/64
      t22588 = t16724/32
      t22589 = t16726/32
      t22590 = t16730/32
      t22591 = t16732/32
      t22592 = t16734/32
      t22593 = t16736/64
      t22594 = t3638*wC22
      t22595 = t720*t22594
      t22598 = t7176*t22594
      t22601 = t22585-t22586-t22587+t22588+t22589+t22590-t22591-t22592-t
     #22593-t16797*t22595/16+t16506*t22598/16
      t22613 = t1096*t22594
      t22616 = t16742/16
      t22617 = t16744/64
      t22618 = t16746/32
      t22619 = t16748/32
      t22620 = t16750/64
      t22622 = t16755/32
      t22623 = t16757/16
      t22626 = t7082*t22594
      t22629 = t16763/32
      t22630 = t16768/16
      t22631 = t16771/128
      t22632 = t16773/32
      t22633 = t16776/64
      t22634 = t16778/256
      t22635 = t16780/128
      t22636 = -t22622+t22623+t13555*t22491/128-t16506*t22626/32-t22629-
     #t22630+t22631-t22632+t22633+t22634-t22635
      t22638 = t16787/512
      t22645 = t16790/16
      t22646 = t16798/16
      t22647 = t16802/64
      t22650 = t7137*t22594
      t22653 = t22638-t16496*t22613/16+t16517*t22613/16-t16487*t22626/12
     #8-t22645-t16794+t22646+t16801-t22647+t17361*t22613/128+t16501*t226
     #50/32
      t22656 = t3609*t22546
      t22661 = t16804/64
      t22662 = t16806/64
      t22672 = t3692*t22520
      t22675 = t720*t22520
      t22678 = -t16506*t22650/32-t16496*t22656/16+t16513*t22650/16-t2266
     #1-t22662-t16487*t22650/128+t13490*t22488/128-t13490*t6943*t22447/1
     #28+t16570*t22650/128-t16496*t22672/16-t16496*t22675/8
      t22683 = t1093*t22546
      t22686 = t16808/64
      t22687 = t16810/128
      t22688 = t16813/64
      t22689 = t16815/64
      t22690 = t6953*t22494
      t22700 = t16821/64
      t22701 = t16570*t22466/64+t16496*t22683/16+t22686+t22687-t22688+t2
     #2689+t16570*t22690/64-t16496*t1093*t22550/16+t16513*t22690/16+t165
     #01*t22466/16+t22700
      t22702 = t16826/128
      t22703 = t16828/64
      t22704 = t16832/64
      t22705 = t16834/128
      t22706 = t16838/128
      t22719 = t16842/32
      t22720 = t16846/32
      t22721 = t22702+t22703-t22704+t22705-t22706+t13595*t6943*t1071*wC2
     #2/256-t16540*t22469/64-t13595*t6943*t1119*wC22/512+t16540*t22461/1
     #28+t22719+t22720
      t22723 = t16850/16
      t22724 = t16854/32
      t22725 = t16858/16
      t22726 = t16865/64
      t22727 = t16869/64
      t22728 = t16874/64
      t22729 = t16877/32
      t22730 = t16883/16
      t22737 = t22723-t22724-t22725-t22726-t22727+t22728-t22729-t22730-t
     #16513*t22525/32-t16501*t22521/32+t16501*t22547/32
      t22743 = t16886/128
      t22744 = t16889/32
      t22745 = t16891/16
      t22755 = t16901/64
      t22756 = t16914/32
      t22757 = t16871*t6943*t22494/128-t16513*t22551/32+t22743-t22744+t2
     #2745-t16501*t22551/32+t16905*t6943*t22480/32-t16501*t22448/16-t164
     #87*t22690/64-t22755+t22756
      t22762 = t16917/32
      t22763 = t16923/32
      t22764 = t16925/128
      t22765 = t16927/32
      t22766 = t6953*t22480
      t22773 = t16934/32
      t22774 = t16937/16
      t22775 = t16940/64
      t22777 = t16943/16
      t22778 = t16947/16
      t22779 = t16949/16
      t22780 = t16961/64
      t22781 = t16964/32
      t22782 = t16966/32
      t22783 = t16978/128
      t22784 = t16981/32
      t22785 = t16983/32
      t22786 = t16986/16
      t22787 = t16999/32
      t22788 = t22777-t22778+t22779-t22780+t22781-t22782+t22783-t22784-t
     #22785+t22786-t22787
      t22790 = t17001/16
      t22791 = t17004/16
      t22792 = t17006/32
      t22793 = t17010/16
      t22794 = t17012/32
      t22795 = t17018/8
      t22796 = t17020/8
      t22797 = t17022/16
      t22798 = t17025/128
      t22799 = t17036/64
      t22800 = t17038/64
      t22801 = t22790-t22791-t22792+t22793+t22794-t22795+t22796+t22797+t
     #22798+t22799-t22800
      t22802 = t17045/16
      t22803 = t17047/8
      t22804 = t17049/16
      t22805 = t17054/64
      t22806 = t17061/16
      t22807 = t17063/16
      t22808 = t17065/32
      t22809 = t17067/32
      t22810 = t17071/128
      t22811 = t17073/32
      t22812 = -t22802+t22803-t22804+t17053+t22805+t22806-t22807-t22808-
     #t22809+t22810+t22811
      t22815 = t17077/32
      t22816 = t17079/64
      t22817 = t17081/128
      t22818 = t17087/16
      t22819 = t17089/32
      t22820 = t17101/32
      t22821 = t17103/128
      t22822 = t17105/128
      t22823 = t17113/128
      t22824 = t17116/32
      t22825 = -t17076-t22815+t22816+t22817+t22818+t22819-t22820+t22821+
     #t22822-t22823+t22824
      t22826 = t17118/64
      t22827 = t17121/16
      t22828 = t17124/8
      t22829 = t17132/32
      t22830 = t17134/64
      t22831 = t17136/128
      t22832 = t17139/16
      t22833 = t17141/128
      t22834 = t17144/32
      t22835 = t17149/16
      t22836 = t17151/16
      t22837 = t22826-t22827-t22828+t22829+t22830-t22831+t22832-t22833+t
     #22834-t22835-t22836
      t22839 = t17155/16
      t22840 = t17157/128
      t22841 = t17159/256
      t22842 = t3970*t22520
      t22845 = t1096*t22520
      t22848 = t17168/128
      t22849 = t17170/128
      t22850 = t17172/512
      t22851 = t17175/64
      t22852 = t17178/64
      t22853 = t17180/32
      t22854 = t22839+t22840-t22841+t16496*t22842/32+t16496*t22845/16+t2
     #2848-t22849-t22850+t22851+t22852+t22853
      t22855 = t17185/16
      t22856 = t17190/8
      t22857 = t17192/64
      t22858 = t17194/32
      t22859 = t17196/16
      t22860 = t17198/32
      t22861 = t17200/16
      t22862 = t17206/16
      t22863 = t17210/64
      t22864 = t17212/32
      t22865 = t17214/32
      t22866 = -t22855+t22856-t22857-t22858-t22859+t22860+t22861-t22862+
     #t22863-t22864+t22865
      t22870 = t17217/8
      t22871 = t17220/4
      t22872 = t17222/16
      t22873 = t17225/8
      t22874 = t17227/16
      t22875 = t17230/8
      t22876 = t17234/32
      t22877 = t17237/16
      t22878 = t17239/16
      t22879 = t17241/16
      t22881 = t17244/32
      t22882 = t17247/16
      t22883 = t17249/64
      t22884 = t17251/128
      t22885 = t17253/32
      t22886 = t17261/64
      t22887 = t17263/128
      t22888 = t17266/64
      t22889 = t17269/64
      t22890 = t17271/64
      t22891 = t17273/64
      t22892 = t22881+t22882-t22883-t22884+t22885+t22886-t22887-t22888+t
     #22889+t22890-t22891
      t22894 = t17275/128
      t22895 = t17280/16
      t22896 = t17286/128
      t22897 = t17289/32
      t22898 = t17291/16
      t22899 = t17294/64
      t22900 = t17296/128
      t22901 = t17301/8
      t22902 = t17303/64
      t22903 = t17315/16
      t22904 = t17318/8
      t22905 = t22894-t22895-t22896-t22897+t22898+t22899+t22900+t22901-t
     #22902+t22903+t22904
      t22908 = t5992*t22546
      t22916 = t17322/32
      t22917 = t17326/32
      t22918 = t17333/128
      t22919 = t17335/8
      t22920 = t17340/16
      t22921 = -t16570*t22452/128-t16496*t22908/32-t16570*t22495/128+t16
     #496*t5992*t22550/32+t22916-t22917-t17330+t17332-t22918-t22919+t229
     #20
      t22924 = t17342/16
      t22925 = t17350/32
      t22926 = t17352/8
      t22927 = t17354/8
      t22928 = t17356/64
      t22929 = t17359/16
      t22930 = t17363/64
      t22931 = t17366/32
      t22932 = t17368/16
      t22933 = -t22924-t17346+t17348-t22925-t22926+t22927-t22928-t22929-
     #t22930-t22931+t22932
      t22934 = t17383/16
      t22935 = t17386/16
      t22936 = t17396/64
      t22937 = t17399/8
      t22938 = t7498*wC22
      t22941 = t17405/512
      t22942 = t17407/8
      t22943 = t17418/32
      t22944 = -t22934+t22935+t17393-t17395+t22936-t22937-t16480*t22938/
     #32-t22941+t22942+t17410+t22943
      t22946 = t17425/8
      t22947 = t17427/8
      t22948 = t17429/256
      t22949 = t17431/128
      t22950 = t17433/64
      t22951 = t17435/8
      t22952 = t17438/16
      t22953 = t17441/8
      t22954 = t17443/32
      t22955 = t17450/256
      t22956 = -t17421+t22946+t22947+t22948+t22949+t22950-t22951-t22952-
     #t22953-t22954-t22955
      t22957 = t17452/128
      t22958 = t17458/16
      t22959 = t17462/32
      t22960 = t17466/64
      t22961 = t17468/256
      t22962 = t17471/8
      t22963 = t17473/32
      t22964 = t17476/8
      t22971 = t22957+t22958-t22959+t22960+t22961+t22962+t22963-t22964+t
     #16540*t22547/128-t16540*t22521/64-t16710*t22466/256
      wI137pc1 = -t22762+t22825+t22701+t22812+t22801+t22788+t22653-t2276
     #3+t22764-t22765+t22515+t22500+t22474+t22944+t22678+t22933+t22971+t
     #22737+t22921+t22956+t22636+t22905+t22892-t22774-t22775-t22878+t228
     #79+t16513*t22766/16-t22875+t22876-t22877+t22873-t22874+t16540*t224
     #66/32-t16570*t22469/64-t16580*t22466/64-t16540*t22448/64+t16710*t2
     #2452/512+t16540*t22455/128-t16540*t22452/64+t16570*t22461/128+t165
     #80*t22448/64-t16506*t22766/16+t16501*t22766/16-t17361*t22595/64+t1
     #6797*t22613/32+t16496*t22595/8-t16517*t22595/8+t16487*t22598/64-t2
     #2871+t22872+t22866+t22870-t22773-t22620-t22618+t22619+t22616+t2261
     #7+t22601+t22854+t22584+t22571+t22559+t22721+t22837+t22534+t22757
      t22989 = -t17484/12+t17487/24-t17491/54-t17494/108+t17496/108-t175
     #00/216+t17502/216+t17504/216-t17508/12-t17510/24-t17514/24-t17520/
     #216-t17522/216
      t22990 = t17526/54
      t22991 = t17528/108
      t22992 = t17530/54
      t22993 = t17532/54
      t22994 = t17535/48
      t22995 = t17537/108
      t22996 = t17539/108
      t22997 = t17543/108
      t22998 = t17545/216
      t22999 = t17549/108
      t23000 = t17551/108
      t23001 = t17554/48
      t23002 = t17556/216
      t23003 = t17559/24
      t23004 = -t22990-t22991-t22992+t22993-t22994-t22995+t22996-t22997+
     #t22998+t22999-t23000-t23001-t23002-t23003
      t23020 = t17563/108
      t23021 = t17565/54
      t23022 = t17569/48
      t23023 = t17572/24
      t23024 = t17574/216
      t23025 = t17576/108
      t23026 = -t16323*t22480/54-t16326*t22480/108-t16138*t22447/216+t16
     #138*t22451/216+t16125*t22494/216+t10712*t22938/432+t17675*t22594/1
     #08+t23020+t23021+t23022+t23023+t23024+t23025
      t23027 = t17578/216
      t23028 = t17581/216
      t23029 = t17583/432
      t23030 = t17587/108
      t23031 = t17589/108
      t23032 = t17591/216
      t23033 = t17599/54
      t23034 = t17601/108
      t23035 = t17605/54
      t23036 = t17607/6
      t23037 = t17612/108
      t23038 = t17614/108
      t23039 = t17618/12
      t23040 = t17620/12
      t23041 = -t23027+t23028+t23029-t23030+t23031+t23032-t23033+t23034+
     #t23035+t23036+t23037-t23038+t23039+t23040
      t23044 = t17622/108
      t23045 = t17624/108
      t23046 = t17628/24
      t23047 = t17630/24
      t23048 = t17633/12
      t23049 = t17635/54
      t23050 = t17637/54
      t23051 = t17639/108
      t23052 = t17644/216
      t23053 = t17650/108
      t23054 = t17652/108
      t23055 = t17656/27
      t23056 = t17662/24
      t23057 = -t23044+t23045-t23046+t23047-t23048-t23049+t23050-t23051+
     #t23052-t23053+t23054+t23055+t23056
      t23058 = t17665/12
      t23059 = t17667/24
      t23060 = t17669/54
      t23061 = t17671/108
      t23062 = t17677/12
      t23063 = t17681/54
      t23064 = t17683/6
      t23065 = t17685/12
      t23066 = t17689/432
      t23067 = t17691/216
      t23068 = t17693/216
      t23069 = t17696/48
      t23070 = t17698/216
      t23071 = t17700/216
      t23072 = t23058-t23059+t23060+t23061+t23062+t23063-t23064-t23065+t
     #23066+t23067-t23068+t23069+t23070+t23071
      t23074 = t17703/24
      t23075 = t17706/108
      t23076 = t17708/108
      t23077 = t17710/108
      t23078 = t17714/216
      t23079 = t17718/54
      t23080 = t17720/54
      t23081 = t17724/27
      t23084 = t17728/108
      t23085 = t17732/108
      t23090 = t23074+t23075+t23076-t23077+t23078-t23079+t23080-t23081-t
     #16044*t22480/108+t23084+t23085-t17493*t22675/216-t17493*t22845/432
      t23119 = -t15975*t22656/108-t15972*t22656/216-t16125*t22656/216+t1
     #6125*t22672/108+t16125*t22842/216-t10712*t22521/432+t10712*t22469/
     #216+t10712*t22461/432-t16138*t22683/108+t16138*t22675/108-t16138*t
     #22908/216+t16138*t22845/216+t16323*t22656/54+t16326*t22656/108
      wI138pc1 = t22989+t23004+t23026+t23041+t23057+t23072+t23090+t23119
      t23132 = t17484/16+t17488+t17491/8+t17494/16-t17496/16+t17500/32-t
     #17502/32-t17504/32-t17509-t17511+t17514/32+t17520/32+t17522/32
      t23147 = -t17526/8+t17528/16-t17530/8+t17532/8+t17754/16+t17535/64
     #-t17537/16+t17539/16-t17543/16+t17545/32+t17549/16-t17551/16+t1755
     #4/64+t17556/32
      t23161 = -t17560+t17563/16+t17565/8-t17569/64-t17572/32-t17574/32-
     #t17576/16+t17578/32-t17581/32-t17583/64+t17587/16-t17589/16-t17591
     #/32
      t23175 = -t17781/32+t17783/32+t17785/64-t17787/16-t17599/8+t17601/
     #16+t17605/8-t17607/8+t17612/16-t17614/16+t17619-t17620/16-t17622/1
     #6+t17624/16
      t23188 = t17800/32-t17802/64-t17629+t17631-t17634-t17635/8+t17808/
     #16+t17810/8-t17812/16-t17814/32+t17816/32+t17637/8-t17639/16
      t23200 = t17644/32-t17650/16+t17652/16-t17656/4+t17663+t17666-t176
     #68-t17669/8-t17671/16-t17677/16-t17681/8+t17683/8+t17685/16-t17689
     #/64
      t23214 = -t17691/32+t17693/32-t17696/64-t17698/32+t17700/32+t17704
     #+t17706/16+t17708/16-t17710/16+t17714/32+t17718/8-t17720/8+t17724/
     #4
      t23229 = t17848/32+t17850/32-t17852/32+t17854/16-t17856/16+t17858/
     #16-t17860/32-t17862/16-t17864/8-t17866/32+t17868/64-t17870/64+t177
     #28/16-t17732/16
      wI139pc1 = t23132+t23147+t23161+t23175+t23188+t23200+t23214+t23229
      t23248 = -t17887+t17520/40+t17522/40+t17526/20-t17528/20+t17530/20
     #-t17532/20+t17894-t17537/40+t17539/40-t17543/40
      t23268 = -t17574/40+t17576/20+t17578/40-t17581/40-t17583/80-t17587
     #/20+t17589/20+t17591/40+t17599/20+t17601/40-t17605/20
      t23284 = t17635/20-t17637/20-t17639/40+t17644/80-t17650/40-t17652/
     #40+t17656/20+t17941-t17942-t17943-t17669/40
      t23303 = t17700/80+t17956+t17706/40+t17708/40-t17710/40+t17714/80+
     #t17718/40+t17720/40-t17724/20+t17728/40-t17732/80
      wI140pc1 = -t17504/40+t17500/40-t17502/40+t17931+t17932+t17496/20+
     #t17491/40-t17494/20+t17485-t17927-t17930+t17909-t17905-t17908+t179
     #03-t17885+t17877+t17884+t17686-t17678-t17681/40-t17684-t17671/80-t
     #17622/40+t17624/40-t17621-t17614/40+t17612/40+t17608+t17563/40-t17
     #565/20+t17556/40-t17549/40-t17551/40+t23284+t23268+t23248+t23303-t
     #17698/40+t17693/40-t17689/80-t17691/40+t17545/80-t17952
      t23306 = t211*wV13
      t23307 = t23306*t3577
      t23308 = t3429*t176
      t23311 = t3429*t448
      t23321 = t23306*t3420
      t23322 = t3615*t448
      t23325 = t3615*t8093
      t23328 = t3623*t18043
      t23338 = t3623*t8093
      t23341 = t23307*t23308/16-t23307*t23311/32-t23307*t3429*t838/32-t2
     #3307*t18020/4+t23307*t18040/4+t23321*t23322/8-t23321*t23325/8-t233
     #07*t23328/4-t23307*t18032/8+t23307*t3623*t3639/8+t23321*t3631/8+t2
     #3307*t23338/8
      t23348 = t3429*t843
      t23366 = t3623*t3737
      t23374 = t23306*t3578*t18048/16-t23307*t3429*t801/32-t23307*t23348
     #/32-t23307*t18048/4+t23307*t18065/4+t23307*t18047*t448/8-t23307*t1
     #8047*t8093/8-t23306*t18068*t18065/4-t23307*t18059/8+t23307*t23366/
     #8+t23321*t18054/8+t23307*t3623*t8186/8
      wI141pc1 = t23341+t23374
      t23375 = t211*wV17
      t23376 = t23375*t3578
      t23387 = t23375*t18068
      t23390 = t23375*t4465
      t23397 = t23375*t3577
      t23408 = t17969*t8186
      t23416 = t23376*t18154/32+t23376*t18160/64-t23376*t18163/64+t23376
     #*t17969*CsTh*t446/64-t23387*t18189/16+t23390*t18184/16-t23387*t182
     #01/8+t23390*t18198/8-t23397*t3429*t4072/256-t23397*t18212/64+t2339
     #7*t18142*t1865/32-t23397*t18143/16-t23397*t23408/32+t23397*t18142*
     #t1224/32+t23397*t18174/16
      t23428 = t3428*t103
      t23430 = t18039*CsTh
      t23453 = t17969*t8257
      t23458 = -t23390*t1242/16+3.D0/16.D0*t23397*t18146-t23397*t18217/6
     #4+t23397*t3623*t4186/64-3.D0/16.D0*t23397*t18204-t23387*t23428*Mzs
     #q*t23430/16+t23397*t3429*t1078/64-t23397*t17969*t8093/16+t23390*t4
     #45*t23430/16-t23397*t3429*t1249/64-t23397*t3429*t1247/64+t23387*t2
     #3408/32-t23387*t18174/16-t23397*t23453/32+t23387*t23453/32
      t23464 = t23428*t18043
      t23478 = t3893*t1318
      t23497 = -t23397*t18160/16-t23397*t18207/64-3.D0/16.D0*t23387*t234
     #64-t23397*t18154/8+3.D0/16.D0*t23397*t23464-t23397*t3429*t1119/256
     #-t23397*t18180/16+t23387*t18180/16+3.D0/64.D0*t23376*t23478-3.D0/1
     #6.D0*t23397*t23478+t23397*t18163/16+t23397*t3623*t4120/64-t23397*t
     #3429*t1251/128-t23397*t3429*t1071/256+t23397*t3623*t4121/64
      t23505 = t449*t3409
      t23519 = t23375*t3420
      t23541 = t23397*t3623*t4253/64-t23397*t3429*t1249*t842/128+t23397*
     #t3623*t23505/16-3.D0/64.D0*t23397*t3429*t1073+t23397*t3623*t448*t3
     #409/16+3.D0/16.D0*t23397*t3623*t18150+t23519*t4166/64-t23519*t3615
     #*t23505/16-t23519*t4141/64+t23519*t4085/64-t23519*t3990/256+t23519
     #*t3615*t1247/16+t23519*t3615*t1249/16+t23519*t4073/64+t23519*t4066
     #/64
      wI142pc1 = t23416+t23458+t23497+t23541
      t23543 = t211*wV5
      t23544 = t23543*t3420
      t23553 = t23543*t3577
      wI143pc1 = t23544*t17970/16-t23544*t17973/4+t23544*t17976/4-t23544
     #*t4476/4-t23553*t3429*t129/16+t23544*t17982/4-t23544*t4487/4+t2355
     #3*t3623*t3692/4
      t23568 = t9430*t3420
      wI144pc1 = t18435-t18450+t18445-t18442+t18437*t23311/16+t18443-t18
     #437*t23338/4-t18446+t18439-t23568*t23322/4-t18444+t18453+t23568*t2
     #3325/4+t18447-t18448-t18436
      t23573 = t211*wV9
      t23574 = t23573*t3420
      t23578 = t23573*t3577
      wI145pc1 = t23574*t3454*t18142/32-t23578*t18040/2+t23574*t3615*t36
     #38/8-t23578*t23308/8+t23578*t23311/16+t23578*t23348/32-t23578*t233
     #38/4-t23578*t3623*t3638/8+t23578*t18020/2-t23574*t23322/4-t23574*t
     #18054/8+t23578*t23328/2+t23574*t23325/4+t23578*t18059/8-t23578*t23
     #366/8-t23574*t3615*t54/8
      t23616 = -t18454-t18455+t18456-t18457+t18458-t4957/88+t18460-t4969
     #/88+t18462+t4981/88-t18464+t18465-t5002/88-t18467+t18468-t18469+t1
     #8470-t18471
      t23624 = t18473-t18474+t18475-t5044/352+t5049/352-t18478+t18479+t5
     #064/88+t18481-t5070/88+t5072/88-t18484-t18485-t18486-t5091/44+t509
     #6/176-t18489+t18490
      t23632 = -t5104/44+t5108/44-t18495+t18496+t18497+t18498-t18499-t51
     #32/352-t18501+t5137/88+t18503-t5142/88+t5145/88-t18506+t18507-t185
     #08+t18509-t18510
      t23640 = -t5164/88+t5167/88+t5169/44-t5173/176-t18516+t18517+t1851
     #8+t5188/44+t18520-t18521-t18522-t18523+t18524-t5213/44-t18526+t522
     #1/352+t18528+t18529
      wI146pc1 = t23616+t23624+t23632+t23640
      t23644 = t3843*t3406
      t23645 = t6932*t23644
      t23648 = t6932*t4197
      t23653 = t167*t6937
      t23654 = t23653*t103
      t23655 = t3606*t23654
      t23658 = t7064*t3406
      t23659 = t6932*t23658
      t23664 = t12359*t23658
      t23667 = wV1*s
      t23669 = t5598*t23667*t3420
      t23672 = -t18706-t18707+t18708+t18709-t18710+t18711-t18712+t18713+
     #9.D0/64.D0*t6966*t8187-9.D0/64.D0*t23645*t8187+15.D0/64.D0*t23648*
     #t8412+3.D0/32.D0*t7000*t8412+t18714+t18715+3.D0/512.D0*t12304*t236
     #55+3.D0/64.D0*t23659*t8282+3.D0/64.D0*t7089*t8282-3.D0/64.D0*t2366
     #4*t8282+3.D0/512.D0*t23669*t8258
      t23673 = wV1*t3420
      t23674 = t5598*t23673
      t23675 = t1953*t3409
      t23676 = t6943*t23675
      t23679 = wV1*t5
      t23681 = t5590*t23679*t3420
      t23682 = t19205*t23675
      t23685 = t6949*t449
      t23686 = t3579*t23685
      t23689 = wV1*t1070
      t23691 = t5598*t23689*t3420
      t23692 = t3579*t8169
      t23698 = t9743*t6943*t103
      t23701 = t3429*t23685
      t23706 = t5598*t23667
      t23708 = t9743*t6943*t449
      t23711 = t6937*t449
      t23712 = t3606*t23711
      t23715 = t5590*wV1
      t23716 = t23715*t23658
      t23719 = 3.D0/256.D0*t23674*t23676-3.D0/64.D0*t23681*t23682-t18716
     #+t18717-t18718+t18719-t18720+t18721-t18723-t18724+t18725+3.D0/256.
     #D0*t7255*t23686+3.D0/512.D0*t23691*t23692+3.D0/2048.D0*t7255*t8170
     #+3.D0/2048.D0*t12391*t23698+3.D0/1024.D0*t23669*t23701-9.D0/256.D0
     #*t12300*t8343+3.D0/1024.D0*t23706*t23708-3.D0/64.D0*t12352*t23712+
     #3.D0/64.D0*t23716*t23682
      t23721 = t6932*t8179
      t23722 = t7137*t23675
      t23727 = t23711*t3409
      t23728 = t3623*t23727
      t23731 = 3.D0/64.D0*t23721*t23722-t18726+t18727+t18728-t18729+t187
     #30-t18731-t18732+t18733-t18734+t18735-t18736+t18737+t18738-t18739-
     #t18742-t18743-3.D0/512.D0*t23669*t8087-3.D0/256.D0*t23674*t23728
      t23732 = t3606*t8070
      t23737 = t5590*t23673
      t23738 = t6943*t23505
      t23753 = t23673*t3406
      t23754 = t10065*t23753
      t23761 = t5590*t23753
      t23767 = t5590*t23673*t3428
      t23772 = 9.D0/256.D0*t12352*t23732-9.D0/256.D0*t6995*t8087+t18744-
     #t18745+t18746+t18747+3.D0/64.D0*t23737*t23738-3.D0/64.D0*t23737*t2
     #3728+3.D0/64.D0*t17989*t23738-9.D0/128.D0*t6966*t8094+9.D0/128.D0*
     #t6966*t8087-15.D0/128.D0*t17989*t8083-3.D0/64.D0*t6942*t8083+3.D0/
     #64.D0*t23754*t8083-3.D0/256.D0*t23669*t8187-3.D0/32.D0*t23754*t192
     #40-9.D0/128.D0*t23761*t19240-15.D0/64.D0*t7008*t8412+9.D0/128.D0*t
     #23767*t8412+9.D0/256.D0*t6995*t8258
      t23775 = t1803*t3409
      t23776 = t6943*t23775
      t23779 = t19205*t23775
      t23788 = t3606*t8308
      t23791 = t6953*t449
      t23792 = t3606*t23791
      t23797 = t5598*t23753
      t23806 = t6949*t23505
      t23821 = -3.D0/128.D0*t23674*t23776+3.D0/32.D0*t23681*t23779+3.D0/
     #32.D0*t23721*t8190-9.D0/64.D0*t6960*t8339+15.D0/64.D0*t6966*t8199+
     #9.D0/128.D0*t12300*t23788+3.D0/32.D0*t12352*t23792+3.D0/256.D0*t23
     #669*t8204+3.D0/128.D0*t23797*t23776+3.D0/64.D0*t7008*t8248-9.D0/25
     #6.D0*t6995*t8261+3.D0/64.D0*t23737*t23676-3.D0/64.D0*t7008*t23806+
     #3.D0/64.D0*t17989*t23676-3.D0/64.D0*t23761*t23676-9.D0/128.D0*t696
     #6*t8258+9.D0/128.D0*t23645*t8258-15.D0/128.D0*t23648*t8441-3.D0/64
     #.D0*t7000*t8441
      t23830 = t23653*t449
      t23831 = t3429*t23830
      t23837 = t3429*t7137*t449
      t23840 = t3579*t23727
      t23843 = 3.D0/64.D0*t23754*t19245+9.D0/256.D0*t23761*t19245+15.D0/
     #128.D0*t7008*t8441-9.D0/256.D0*t23767*t8441-3.D0/1024.D0*t7255*t23
     #831+t18748+t18749+t18750+t18751-t18752+t18753-t18754+t18755-t18756
     #+t18757+t18759-t18760-3.D0/256.D0*t17981*t8116-3.D0/256.D0*t23737*
     #t23837-3.D0/256.D0*t23669*t23840
      t23851 = t6943*t1803
      t23856 = t3579*t8070
      t23861 = t3623*t167
      t23862 = t23861*t23727
      t23865 = t23861*t8070
      t23870 = 3.D0/256.D0*t12352*t23708-3.D0/32.D0*t7008*t8195+3.D0/32.
     #D0*t6942*t8213-3.D0/32.D0*t6966*t23851-3.D0/32.D0*t6975*t8199-3.D0
     #/512.D0*t7255*t23856+9.D0/256.D0*t23761*t8083+3.D0/256.D0*t23669*t
     #23862+3.D0/512.D0*t7255*t23865+9.D0/128.D0*t7129*t8343-t18761-t187
     #62-t18763-t18764+t18765-t18766+t18767+t18768-t18769
      t23879 = t6943*t1953
      t23886 = t18038*t23685
      t23889 = t7137*t1953
      t23900 = t23673*t3426
      t23901 = t10065*t23900
      t23904 = t18038*t23791
      t23907 = t3843*t3426
      t23908 = t6932*t23907
      t23911 = t18770-t18771-9.D0/1024.D0*t6995*t8056+9.D0/1024.D0*t1230
     #0*t23698-3.D0/512.D0*t12391*t8343-3.D0/64.D0*t6942*t8276+3.D0/64.D
     #0*t6966*t23879+t18772+3.D0/64.D0*t6975*t8252+9.D0/256.D0*t6995*t82
     #52+3.D0/64.D0*t23737*t23886-3.D0/64.D0*t7008*t23889+9.D0/512.D0*t6
     #995*t8430-9.D0/512.D0*t7054*t8309+3.D0/256.D0*t12391*t23788+3.D0/1
     #28.D0*t23706*t23792+t18773-3.D0/128.D0*t23901*t8213-3.D0/128.D0*t1
     #7981*t23904-15.D0/256.D0*t23908*t8199
      t23917 = t3429*t23791
      t23922 = t5590*t23900
      t23925 = t3579*t6937
      t23926 = t23925*t23675
      t23929 = t23925*t8257
      t23934 = t5598*wV1
      t23935 = t23934*t23644
      t23938 = t3752*t3406
      t23939 = t23934*t23938
      t23948 = t3606*t8169
      t23951 = t3606*t23685
      t23954 = 3.D0/128.D0*t6975*t8430+3.D0/128.D0*t6966*t23917+9.D0/256
     #.D0*t6960*t8309-3.D0/256.D0*t23922*t23886-3.D0/256.D0*t23669*t2392
     #6-3.D0/512.D0*t7255*t23929-t18774+t18775-t18779+t18780-3.D0/256.D0
     #*t17981*t8276+3.D0/256.D0*t23935*t23926+3.D0/512.D0*t23939*t23929+
     #3.D0/256.D0*t6995*t23701+t18781-9.D0/1024.D0*t6995*t8364+9.D0/1024
     #.D0*t7054*t8170-3.D0/512.D0*t12391*t23948-3.D0/256.D0*t23706*t2395
     #1
      t23967 = t18038*t23711
      t23972 = -t18782+t18783-t18784+t18785+t18786-t18787+3.D0/256.D0*t2
     #3901*t8276+3.D0/256.D0*t17981*t23886+15.D0/512.D0*t23908*t8252+t18
     #788-3.D0/256.D0*t6975*t8364-3.D0/256.D0*t6966*t23701-9.D0/512.D0*t
     #6960*t8170+3.D0/256.D0*t17981*t23967+15.D0/512.D0*t6966*t8056+t187
     #89-t18790+t18791+t18792+t18793
      t23978 = t3429*t23711
      t23983 = 15.D0/128.D0*t7008*t7995-3.D0/256.D0*t6975*t8056-3.D0/256
     #.D0*t6966*t23978-9.D0/512.D0*t6960*t8061+t18794-t18796-t18797-t187
     #98-t18799-t18800+t18801-t18802-t18803-t18804-t18805-t18806+t18807-
     #t18808+t18809
      t23984 = t18810+t18811+t18812-t18815+t18816-t18817-t18818-t18819-t
     #18820+t18821-t18822+t18823+t18824-t18825-t18826+t18827+t18828-t188
     #29+t18830-t18832
      t23997 = t7137*t1803
      t24002 = t7137*t23775
      t24013 = -t18833+t18834-t18835+t18836+3.D0/64.D0*t7148*t23712-3.D0
     #/64.D0*t7133*t23967-t18837+3.D0/64.D0*t7133*t8116-t18838-t18839-9.
     #D0/128.D0*t6995*t8199-3.D0/32.D0*t23737*t23904+3.D0/32.D0*t7008*t2
     #3997-3.D0/32.D0*t23716*t23779-3.D0/32.D0*t23721*t24002+3.D0/32.D0*
     #t23659*t23779-3.D0/32.D0*t23659*t8217-3.D0/32.D0*t7089*t8217+3.D0/
     #32.D0*t23664*t8217
      t24014 = t3429*t23654
      t24017 = t3429*t167
      t24018 = t24017*t23791
      t24021 = t24017*t8308
      t24024 = t24017*t23685
      t24027 = t24017*t8169
      t24040 = t18038*t23727
      t24045 = t18840-t18841-3.D0/2048.D0*t23691*t24014+3.D0/512.D0*t725
     #5*t24018+3.D0/1024.D0*t23691*t24021-t18842-3.D0/1024.D0*t7255*t240
     #24-3.D0/2048.D0*t23691*t24027+3.D0/64.D0*t23659*t8130+3.D0/64.D0*t
     #7089*t8130-15.D0/128.D0*t6966*t8119+t18843+t18844-t18845+t18846+t1
     #8847+3.D0/512.D0*t23706*t23732+3.D0/256.D0*t23674*t23738-3.D0/64.D
     #0*t23681*t24040-9.D0/256.D0*t23737*t8083
      t24069 = t3606*t23830
      t24074 = t23925*t23775
      t24077 = t23925*t8186
      t24088 = -3.D0/64.D0*t7095*t8130-9.D0/128.D0*t6995*t8187+9.D0/128.
     #D0*t6995*t8204-3.D0/32.D0*t23737*t23776+3.D0/32.D0*t23761*t23776-3
     #.D0/64.D0*t7089*t8136+3.D0/128.D0*t17981*t8213+3.D0/128.D0*t23922*
     #t23904-3.D0/64.D0*t23659*t23682-3.D0/256.D0*t23706*t23712+3.D0/256
     #.D0*t6942*t8485+3.D0/256.D0*t12391*t24069+9.D0/128.D0*t6960*t8375+
     #3.D0/128.D0*t23669*t24074+3.D0/256.D0*t7255*t24077-3.D0/128.D0*t23
     #935*t24074-3.D0/256.D0*t23939*t24077-3.D0/64.D0*t23721*t8245-15.D0
     #/128.D0*t6966*t8252
      t24107 = t3623*Mzsq
      t24108 = t24107*t23727
      t24111 = t3579*t23791
      t24114 = t3579*t8308
      t24125 = t6953*t23505
      t24132 = -9.D0/256.D0*t12300*t23948-t18848-3.D0/64.D0*t12352*t2395
     #1-3.D0/512.D0*t23669*t8261-3.D0/256.D0*t23797*t23676-3.D0/64.D0*t6
     #942*t8116-3.D0/128.D0*t6995*t23917+3.D0/64.D0*t12384*t19209+9.D0/2
     #56.D0*t12352*t19209+3.D0/64.D0*t23737*t23967+3.D0/64.D0*t23681*t24
     #108-3.D0/128.D0*t7255*t24111-3.D0/256.D0*t23691*t24114-3.D0/1024.D
     #0*t7255*t8309-3.D0/512.D0*t23669*t23917+3.D0/64.D0*t7095*t24040-3.
     #D0/64.D0*t23659*t24040+3.D0/32.D0*t7008*t24125-3.D0/32.D0*t17989*t
     #23776-3.D0/64.D0*t7133*t23738
      wI147pc1 = t23672+t23719+t23731+t23772+t23821+t23843+t23870+t23911
     #+t23954+t23972+t23983+t23984+t24013+t24045+t24088+t24132
      t24136 = t211*wV11
      t24137 = t24136*t5930
      t24138 = t8517*t167
      t24140 = t24138*t3606*t449
      t24143 = t24136*t3388
      t24146 = t24136*t8548
      t24149 = t24136*t8578
      t24150 = t9743*t103
      t24151 = t18177*t24150
      t24154 = t24136*t175
      t24157 = t24136*t1070
      t24160 = t24136*t8684
      t24165 = t24136*t8562
      t24166 = t3420*t103
      t24167 = t18177*t24166
      t24170 = t8578*Mzsq
      t24171 = t24136*t24170
      t24172 = t3606*t8257
      t24175 = t4197*t8257
      t24178 = t8517*Mzsq
      t24179 = t24136*t24178
      t24182 = t5936*t8517
      t24183 = t24136*t24182
      t24184 = t9743*t449
      t24185 = t3579*t24184
      t24188 = t24136*t8517
      t24189 = t3606*t23505
      t24194 = t24136*t8679
      t24197 = -t24137*t24140/160-9.D0/80.D0*t24143*t8559-t24146*t8551/2
     #0+3.D0/160.D0*t24149*t24151+3.D0/32.D0*t24154*t8559-t24157*t8825/8
     #0+t24160*t8805/20+9.D0/320.D0*t24143*t8654+3.D0/40.D0*t24165*t2416
     #7-7.D0/80.D0*t24171*t24172+7.D0/80.D0*t24171*t24175+t24179*t24172/
     #20+t24183*t24185/1280+t24188*t24189/16-3.D0/128.D0*t24154*t8654+t2
     #4194*t8711/80
      t24200 = t3767*t449
      t24201 = t24017*t24200
      t24204 = t4678*t449
      t24205 = t24017*t24204
      t24208 = t12308*t8517
      t24209 = t24136*t24208
      t24210 = t24166*t3409
      t24211 = t18177*t24210
      t24215 = t18177*t3406*t24210
      t24218 = t5274*t8517
      t24219 = t24136*t24218
      t24220 = t3579*t3420
      t24221 = t24220*t23775
      t24228 = t24136*t8674
      t24231 = t24136*t1067
      t24234 = t24136*t8644
      t24235 = t3767*t103
      t24236 = t3429*t24235
      t24239 = t8578*t62
      t24240 = t24136*t24239
      t24243 = t8554*Mzsq
      t24244 = t24136*t24243
      t24247 = t24136*t8573
      t24248 = t18038*t3426
      t24249 = t24248*t24235
      t24252 = t24136*t8554
      t24253 = t3429*t24200
      t24260 = -t24165*t8711/80-t24183*t24201/640+t24183*t24205/1280+3.D
     #0/40.D0*t24209*t24211-3.D0/40.D0*t24209*t24215-t24219*t24221/160+t
     #24194*t8748/10-t24165*t8748/10-t24228*t8773/10+3.D0/20.D0*t24231*t
     #8762+t24234*t24236/40-3.D0/80.D0*t24240*t24236-t24244*t24236/40+7.
     #D0/160.D0*t24247*t24249-t24252*t24253/32-3.D0/16.D0*t24154*t8759+9
     #.D0/40.D0*t24143*t8759
      t24266 = t24136*s
      t24267 = t1387*t8517
      t24268 = t24267*t167
      t24269 = t24266*t24268
      t24270 = t4197*t23775
      t24273 = t5299*t8517
      t24274 = t24136*t24273
      t24277 = t24136*t8725
      t24280 = t24136*t5910
      t24287 = t24220*t23675
      t24296 = t3429*t24204
      t24299 = t4678*t103
      t24300 = t3429*t24299
      t24307 = t24160*t8773/10-t24160*t8751/10+t24269*t24270/160+7.D0/32
     #0.D0*t24274*t24253-t24277*t24236/160+t24280*t8524/80+t24228*t8751/
     #10-3.D0/40.D0*t24234*t8778+t24219*t24287/320-t24194*t8681/20+t2416
     #5*t8681/20-t24280*t8531/160-7.D0/640.D0*t24274*t24296+t24277*t2430
     #0/320-9.D0/80.D0*t24143*t8696-t24160*t8688/20
      t24314 = t4197*t23675
      t24327 = t24248*t24299
      t24334 = t18038*t24150
      t24337 = t3606*t8093
      t24342 = t3420*t449
      t24343 = t3429*t24342
      t24346 = t3429*t24166
      t24349 = t24228*t8688/20-t24228*t8676/20-3.D0/40.D0*t24231*t8699-t
     #24269*t24314/320+3.D0/32.D0*t24154*t8696+t24160*t8676/20+3.D0/160.
     #D0*t24240*t24300+t24244*t24300/80+t24252*t24296/64-7.D0/320.D0*t24
     #247*t24327-t24234*t24300/80+3.D0/80.D0*t24234*t8716+t24252*t24334/
     #80+t24179*t24337/20-7.D0/320.D0*t24247*t24334+t24252*t24343/64-t24
     #234*t24346/80
      t24352 = t168*t8517
      t24353 = t24136*t24352
      t24354 = t3606*t23775
      t24358 = t3579*t9743*t1763
      t24361 = t3606*t8186
      t24364 = t24352*t3428
      t24365 = t24136*t24364
      t24366 = t9887*t23775
      t24369 = t4197*t8186
      t24372 = t3606*t1803
      t24375 = t18038*t24235
      t24380 = t3606*t1763
      t24393 = t24352*t167
      t24394 = t24136*t24393
      t24395 = t3420*t1763
      t24396 = t24395*t3409
      t24397 = t3623*t24396
      t24400 = t8653*t8257
      t24403 = 7.D0/80.D0*t24353*t24354+t24219*t24358/320-t24274*t24361/
     #40-7.D0/80.D0*t24365*t24366+t24274*t24369/40+t24252*t24372/8-7.D0/
     #40.D0*t24247*t24375-3.D0/16.D0*t24240*t4751-t24188*t24380/40-3.D0/
     #40.D0*t24244*t4751+t24277*t4751/80-3.D0/20.D0*t24247*t3812+3.D0/20
     #.D0*t24247*t8597-3.D0/40.D0*t24160*t24211+t24394*t24397/80+3.D0/40
     #.D0*t24209*t24400
      t24406 = t24136*t5
      t24407 = t24406*t24239
      t24410 = t3606*t1865
      t24413 = t18038*t24299
      t24416 = t8518*t24299
      t24421 = t24136*t5299
      t24422 = t8518*t24204
      t24427 = t24136*t8721
      t24430 = t3579*t24200
      t24439 = t24136*t8534
      t24442 = t24136*t8744
      t24445 = t24136*t8902
      t24446 = t1763*t3409
      t24447 = t24446*t842
      t24448 = t9887*t24447
      t24453 = 3.D0/40.D0*t24160*t24215-3.D0/40.D0*t24407*t24175+t24234*
     #t24410/20-t24252*t24413/20-t24137*t24416/80-t24280*t8670/320+7.D0/
     #160.D0*t24421*t24422+3.D0/160.D0*t24234*t8658-9.D0/160.D0*t24427*t
     #8662+t24183*t24430/160+3.D0/64.D0*t24247*t8662-t24194*t8708/40+t24
     #165*t8708/40-t24439*t4742/160-3.D0/80.D0*t24442*t8658-t24445*t2444
     #8/40+t24280*t8541/640
      t24459 = t3579*t24204
      t24462 = t24266*t24393
      t24467 = t3606*t1958
      t24470 = t1763*t842
      t24471 = t3606*t24470
      t24474 = t24446*CsTh
      t24475 = t3606*t24474
      t24478 = t3606*t24447
      t24481 = t15618*t8517
      t24482 = t24136*t24481
      t24489 = t24136*t8901
      t24490 = t3606*t24446
      t24495 = t24395*CsTh
      t24496 = t3579*t24495
      t24499 = t8518*t24342
      t24502 = -3.D0/320.D0*t24234*t8666+9.D0/320.D0*t24427*t8704-t24183
     #*t24459/320+t24462*t24270/80-t24462*t24314/160-t24353*t24467/20+t2
     #4353*t24471/40+t24394*t24475/40-t24394*t24478/80+t24482*t24475/20-
     #t24482*t24478/40+t24482*t24397/40+t24489*t24490/40-t24489*t24397/4
     #0-t24274*t24496/40-t24266*t24499/16
      t24507 = t24136*t168
      t24508 = t8518*t24395
      t24513 = t24136*t54
      t24514 = t8518*t24166
      t24524 = t24136*t177
      t24527 = t3429*t24395
      t24538 = t24136*t5936
      t24544 = t24178*t3606*t103
      t24547 = t24157*t8524/40-t24157*t8531/80+t24507*t24508/40+t24157*t
     #8528/320+t24513*t24514/20+3.D0/80.D0*t24513*t8622+t24513*t8625/80-
     #t24136*t1812*t8520/20+t24524*t8520/20+t24188*t24527/160-t24157*t91
     #83/160+t24188*t24467/20-t24188*t24471/40-t24234*t4633/80-t24538*t2
     #4140/320-3.D0/40.D0*t24149*t24167+7.D0/80.D0*t24154*t24544
      t24553 = t24342*t3409
      t24554 = t23861*t24553
      t24567 = t8562*t62
      t24568 = t24136*t24567
      t24571 = t8653*t8186
      t24576 = t3606*t1224
      t24581 = t8518*t24235
      t24584 = t8518*t24200
      t24591 = -9.D0/40.D0*t24165*t8602-t24274*t24554/160+3.D0/32.D0*t24
     #149*t8788+3.D0/80.D0*t24252*t8820-t24137*t8796/160+3.D0/40.D0*t241
     #54*t8545+9.D0/80.D0*t24165*t8564+9.D0/40.D0*t24568*t4751-3.D0/20.D
     #0*t24209*t24571+3.D0/20.D0*t24407*t24369-t24234*t24576/10+t24252*t
     #24375/10+t24137*t24581/40-7.D0/80.D0*t24421*t24584+t24234*t4781/80
     #-3.D0/160.D0*t24568*t24300
      t24592 = t24136*t24268
      t24593 = t4197*t24474
      t24596 = t4197*t24447
      t24615 = t24136*t24267
      t24616 = t3579*t24396
      t24623 = t18177*t24235
      t24626 = t18177*t24299
      t24633 = t24592*t24593/40-t24592*t24596/80+t24277*t24430/80-t24277
     #*t24459/160-t24592*t24475/40-3.D0/32.D0*t24149*t8580-3.D0/80.D0*t2
     #4252*t8570+t24137*t8754/160-t24234*t4778/80+t24592*t24478/80+t2461
     #5*t24616/80-t24274*t24221/80+t24274*t24287/160-3.D0/20.D0*t24165*t
     #24623+3.D0/40.D0*t24165*t24626+3.D0/20.D0*t24160*t24571-3.D0/40.D0
     #*t24160*t24400
      t24637 = t24406*t24567
      t24652 = t24136*t5274
      t24653 = t24138*t24380
      t24656 = t24017*t24495
      t24661 = t3606*t1953
      t24672 = -t24394*t24593/40-3.D0/20.D0*t24637*t24369+3.D0/40.D0*t24
     #637*t24175-3.D0/128.D0*t24247*t8704+t24439*t4882/320+t24194*t8741/
     #80-t24165*t8741/80+3.D0/160.D0*t24442*t8666-t24652*t24653/80-t2421
     #9*t24656/160-3.D0/40.D0*t24149*t24626-t24252*t24661/16+7.D0/80.D0*
     #t24247*t24413+3.D0/32.D0*t24240*t4778+t24394*t24596/80-t24482*t245
     #93/20
      t24679 = t9887*t24474
      t24686 = t24136*t8518
      t24687 = t3767*t23505
      t24692 = t24395*t842
      t24693 = t24017*t24692
      t24698 = t24136*t15618
      t24699 = t8518*t24396
      t24704 = t18038*t24166
      t24711 = t8518*t24210
      t24714 = t3623*t24210
      t24717 = t24482*t24596/40-t24489*t24475/20+t24489*t24478/40+t24445
     #*t24679/20-3.D0/40.D0*t24247*t8820-9.D0/80.D0*t24165*t8788-t24686*
     #t24687/8+t24686*t24366/8+t24219*t24693/320-t24353*t24616/80-t24698
     #*t24699/40+3.D0/160.D0*t24231*t8640-t24252*t24704/20-t24137*t24514
     #/80+7.D0/160.D0*t24421*t24499+t24421*t24711/80-t24274*t24714/80
      t24730 = t3623*t24553
      t24739 = t9887*t23505
      t24742 = t18038*t24210
      t24747 = t3579*t24553
      t24756 = 3.D0/80.D0*t24244*t4778-t24277*t4778/160+3.D0/40.D0*t2424
     #7*t8570-3.D0/40.D0*t24247*t8575-7.D0/160.D0*t24353*t24189+7.D0/160
     #.D0*t24353*t24730-t24252*t24337/20+t24252*t24714/20-t24179*t24714/
     #20-t24686*t24739/16-7.D0/80.D0*t24149*t24742+7.D0/80.D0*t24171*t24
     #714+t24274*t24747/160-3.D0/32.D0*t24149*t8564+t24252*t24361/10-3.D
     #0/80.D0*t24252*t3818
      t24763 = t24136*t8980
      t24764 = t9887*t8186
      t24767 = t24178*t3428
      t24768 = t24136*t24767
      t24780 = t8518*t24184
      t24797 = t24137*t8625/160+t24219*t24747/320+t24421*t24653/80-t2476
     #3*t24764/10+t24768*t24764/10+7.D0/40.D0*t24171*t24361-t24179*t2436
     #1/10+t24136*t8673*t8551/20-t24274*t24358/320-7.D0/640.D0*t24421*t2
     #4780-t24234*t4606/40+3.D0/16.D0*t24149*t8602+3.D0/40.D0*t24252*t38
     #12-t24228*t8805/20-3.D0/40.D0*t24231*t8622-t24219*t24554/320+t2427
     #7*t24185/640
      t24827 = t3579*t24692
      t24830 = t8518*t24150
      t24833 = -t24353*t24527/160+9.D0/80.D0*t24165*t8580-9.D0/80.D0*t24
     #568*t4778+t24280*t8528/640-3.D0/320.D0*t24234*t8650-3.D0/160.D0*t2
     #4165*t24151-t24592*t24397/80-t24277*t24201/320+t24277*t24205/640-7
     #.D0/40.D0*t24171*t24369+t24274*t24172/80+t24274*t24656/160-t24274*
     #t24693/320+t24219*t24496/40-t24219*t24827/80+t24137*t24830/320
      t24836 = t3606*t23675
      t24839 = t9887*t23675
      t24842 = t4678*t23505
      t24857 = t9887*t8257
      t24862 = t3429*t24495
      t24865 = t3429*t24692
      t24868 = t9743*t1958
      t24871 = t9743*t24470
      t24876 = -t24274*t24175/80-7.D0/160.D0*t24353*t24836+7.D0/160.D0*t
     #24365*t24839+t24686*t24842/16-t24137*t8733/80+t24234*t4751/40+3.D0
     #/20.D0*t24149*t24623+t24274*t24827/80-t24686*t24839/16-t24252*t241
     #72/20+t24763*t24857/20-t24768*t24857/20+t24353*t24862/80-t24353*t2
     #4865/160-t24686*t24868/80+t24686*t24871/160+3.D0/80.D0*t24568*t242
     #36
      wI148pc1 = t24197+t24260+t24307+t24349+t24403+t24453+t24502+t24547
     #+t24591+t24633+t24672+t24717+t24756+t24797+t24833+t24876
      wI149pc1 = wV35
      t24897 = -13.D0/11264.D0*t8021+13.D0/22528.D0*t8023+t19023+25.D0/5
     #632.D0*t8028-27.D0/2816.D0*t8033+93.D0/1408.D0*t8035-t8038/128-13.
     #D0/5632.D0*t8041-35.D0/2816.D0*t8043+27.D0/22528.D0*t8045+39.D0/14
     #08.D0*t8048
      t24910 = -t19030+41.D0/2816.D0*t8053-93.D0/5632.D0*t8057+t8062/281
     #6-25.D0/11264.D0*t8064-t8066/176+13.D0/5632.D0*t8068-25.D0/5632.D0
     #*t8072+35.D0/704.D0*t8074-9.D0/176.D0*t8133*t7217+39.D0/704.D0*t80
     #76
      t24922 = 13.D0/5632.D0*t8078-25.D0/1408.D0*t8080+27.D0/2816.D0*t80
     #84-t8088/128-t8090/704+t8095/704+35.D0/5632.D0*t8097+t8101/512-13.
     #D0/2816.D0*t8103-39.D0/2816.D0*t8105+27.D0/2816.D0*t8107
      t24932 = wV4*t1070
      t24933 = t24932*t3420
      t24934 = t5590*t24933
      t24937 = -35.D0/1408.D0*t8109+13.D0/2816.D0*t8111+25.D0/1408.D0*t8
     #113-t19055+93.D0/1408.D0*t8120+73.D0/704.D0*t8123+t19058-73.D0/704
     #.D0*t8134-t19060+27.D0/5632.D0*t8139+9.D0/2816.D0*t24934*t23692
      t24952 = 9.D0/352.D0*t8133*t7068+9.D0/1408.D0*t8040*t23686+9.D0/70
     #4.D0*t8086*t24074+35.D0/1408.D0*t8141-39.D0/1408.D0*t8144-27.D0/56
     #32.D0*t8146-35.D0/1408.D0*t8149+35.D0/2816.D0*t8151+t8154/352-35.D
     #0/2816.D0*t8156-t8158/352
      t24962 = t15491*t23644
      t24975 = 35.D0/1408.D0*t8161+27.D0/5632.D0*t8163-9.D0/704.D0*t8040
     #*t24111-9.D0/704.D0*t8082*t23967-9.D0/1408.D0*t8086*t23926-9.D0/70
     #4.D0*t24962*t24074+9.D0/1408.D0*t24962*t23926+9.D0/2816.D0*t8040*t
     #23865-9.D0/1408.D0*t12898*t24069-9.D0/2816.D0*t13037*t23655+9.D0/1
     #408.D0*t8086*t23862
      t24978 = wV4*t5
      t24980 = t5590*t24978*t3420
      t24983 = t7992*t23658
      t24988 = t7992*t9887
      t24993 = wV4*t3420
      t24994 = t24993*t3406
      t24995 = t5590*t24994
      t25006 = -9.D0/704.D0*t8409*t24040-9.D0/704.D0*t24980*t24108+9.D0/
     #704.D0*t24983*t24040-9.D0/352.D0*t7993*t24125+9.D0/352.D0*t24988*t
     #23776-9.D0/1408.D0*t8082*t23676+9.D0/1408.D0*t24995*t23676+9.D0/70
     #4.D0*t7993*t23806-9.D0/704.D0*t24988*t23676-9.D0/1408.D0*t8082*t23
     #738+9.D0/1408.D0*t8082*t23728
      t25016 = t5598*t24933
      t25024 = t5590*t24993*t3426
      t25044 = -9.D0/704.D0*t8336*t24074-27.D0/5632.D0*t8167-25.D0/22528
     #.D0*t8171+35.D0/2816.D0*t8173+t8175/352-t8177/352-41.D0/704.D0*t81
     #82+t19079-t8188/64-41.D0/352.D0*t8191+41.D0/352.D0*t8196
      t25055 = t5598*wV4
      t25056 = t25055*t23644
      t25059 = t25055*t23938
      t25064 = 27.D0/1408.D0*t8200-93.D0/704.D0*t8202+t8205/64+t8208/352
     #-t8210/352+t19087-9.D0/1408.D0*t8007*t24077+9.D0/176.D0*t8115*t721
     #3+9.D0/704.D0*t25056*t24074+9.D0/1408.D0*t25059*t24077-9.D0/704.D0
     #*t12902*t23792
      t25075 = t7992*t9743
      t25082 = -t19088-73.D0/352.D0*t8220+t19092+73.D0/352.D0*t8227+9.D0
     #/2816.D0*t8086*t23917-9.D0/352.D0*t8115*t7076-9.D0/176.D0*t8128*t7
     #224+9.D0/176.D0*t8128*t7217+9.D0/1408.D0*t25075*t23904-9.D0/1408.D
     #0*t8047*t23917+9.D0/2816.D0*t25024*t23886
      t25095 = t5598*t24994
      t25108 = 9.D0/1408.D0*t8336*t23926+9.D0/2816.D0*t8007*t23929+9.D0/
     #1408.D0*t12902*t23951-9.D0/352.D0*t24980*t23779-9.D0/704.D0*t12958
     #*t23792-9.D0/704.D0*t25095*t23776-9.D0/704.D0*t8047*t23879-9.D0/70
     #4.D0*t8082*t23886-9.D0/1408.D0*t25056*t23926-9.D0/2816.D0*t25059*t
     #23929-9.D0/5632.D0*t8086*t23701
      t25115 = t15491*t23658
      t25132 = 9.D0/352.D0*t8128*t7071-9.D0/352.D0*t8128*t7068+9.D0/704.
     #D0*t7993*t23889-9.D0/704.D0*t25115*t23682-9.D0/704.D0*t8180*t23722
     #+9.D0/2816.D0*t8040*t24018-9.D0/5632.D0*t8040*t24024+9.D0/1408.D0*
     #t12906*t24069+9.D0/2816.D0*t12980*t23655+9.D0/352.D0*t8133*t7170-9
     #.D0/704.D0*t24988*t23738
      t25154 = -9.D0/2816.D0*t8040*t23856-9.D0/1408.D0*t8086*t23840+9.D0
     #/2816.D0*t8082*t23837-9.D0/5632.D0*t12958*t23708+25.D0/2816.D0*t82
     #30-25.D0/11264.D0*t8233+9.D0/1408.D0*t8336*t23840+9.D0/2816.D0*t80
     #07*t23856-9.D0/1408.D0*t8336*t23862-9.D0/5632.D0*t8040*t23831-9.D0
     #/2816.D0*t8007*t23865
      t25157 = t15491*t23938
      t25163 = t3606*t7096
      t25180 = 9.D0/1408.D0*t8040*t24077-9.D0/1408.D0*t25157*t24077+9.D0
     #/704.D0*t8122*t6984+9.D0/352.D0*t9544*t24978*t25163-9.D0/352.D0*t8
     #128*t7170+9.D0/1408.D0*t12902*t23712-9.D0/1408.D0*t8122*t6939+9.D0
     #/176.D0*t8133*t7224-9.D0/176.D0*t8122*t7213-9.D0/2816.D0*t8040*t23
     #929+9.D0/2816.D0*t25157*t23929
      t25191 = t5598*t24993
      t25220 = -9.D0/2816.D0*t25075*t23967+9.D0/2816.D0*t8047*t23978-9.D
     #0/1408.D0*t25191*t23738-73.D0/352.D0*t8241+41.D0/352.D0*t8243+41.D
     #0/704.D0*t8246-41.D0/704.D0*t8249-27.D0/2816.D0*t8253-t8255/1408+t
     #8259/128+9.D0/704.D0*t24980*t24040
      t25237 = 9.D0/704.D0*t8082*t23776-9.D0/704.D0*t24995*t23776+9.D0/7
     #04.D0*t24983*t23682-t8262/128-t8265/704+t8271/704+93.D0/1408.D0*t8
     #273-t19109+73.D0/704.D0*t8279-9.D0/1408.D0*t25191*t23676+9.D0/704.
     #D0*t24980*t23682
      t25253 = 9.D0/1408.D0*t12958*t23951+9.D0/1408.D0*t25095*t23676+t19
     #111-t19112+25.D0/2816.D0*t8288-73.D0/704.D0*t8290-41.D0/704.D0*t82
     #92+9.D0/704.D0*t8007*t24111+9.D0/1408.D0*t25016*t24114-9.D0/704.D0
     #*t8115*t6984+9.D0/2816.D0*t8336*t23917
      t25266 = 73.D0/704.D0*t8294-13.D0/5632.D0*t8296-35.D0/1408.D0*t829
     #8+t19123+t8303/352+35.D0/2816.D0*t8306-t8310/256-25.D0/2816.D0*t83
     #13-t8315/1408-t19129+9.D0/704.D0*t8237*t23738
      t25267 = t3606*t6938
      t25291 = 9.D0/352.D0*t12994*t25267-9.D0/352.D0*t9544*wV4*t25267-9.
     #D0/704.D0*t8305*t23712+9.D0/704.D0*t8237*t23967-9.D0/1408.D0*t8007
     #*t23686-9.D0/2816.D0*t25016*t23692+9.D0/1408.D0*t8115*t7341-9.D0/5
     #632.D0*t8336*t23701+9.D0/1408.D0*t12958*t23712-9.D0/5632.D0*t12902
     #*t23708+9.D0/1408.D0*t25191*t23728
      t25309 = -9.D0/352.D0*t13067*t25163+9.D0/352.D0*t8047*t23851+9.D0/
     #352.D0*t8082*t23904-9.D0/352.D0*t7993*t23997+9.D0/352.D0*t25115*t2
     #3779+39.D0/1408.D0*t8319+39.D0/704.D0*t8321+27.D0/2816.D0*t8323-25
     #.D0/2816.D0*t8325+25.D0/5632.D0*t8328+t8332/2816
      t25319 = t8334/512+25.D0/5632.D0*t8337+t8340/352-t8344/704-25.D0/2
     #816.D0*t8346+73.D0/1408.D0*t8349+35.D0/5632.D0*t8351-t19143+35.D0/
     #704.D0*t8355-93.D0/5632.D0*t8358+t19148
      t25345 = 35.D0/1408.D0*t8388+t8390/176-27.D0/2816.D0*t8395-35.D0/7
     #04.D0*t8397-39.D0/704.D0*t8399+t8402/128-27.D0/2816.D0*t8405+t8407
     #/352-39.D0/704.D0*t19167+41.D0/704.D0*t8410+93.D0/704.D0*t8413
      t25358 = -35.D0/1408.D0*t8415-25.D0/2816.D0*t8417-93.D0/704.D0*t84
     #20+13.D0/22528.D0*t8422+39.D0/1408.D0*t19178+27.D0/1408.D0*t8426-3
     #9.D0/1408.D0*t8428-27.D0/5632.D0*t8431+39.D0/1408.D0*t19185+39.D0/
     #2816.D0*t19189-39.D0/5632.D0*t19191
      t25370 = -41.D0/1408.D0*t8433+25.D0/11264.D0*t8435-27.D0/1408.D0*t
     #8437-35.D0/11264.D0*t8439-39.D0/1408.D0*t19196+39.D0/1408.D0*t1919
     #9-93.D0/1408.D0*t8442+27.D0/11264.D0*t8444-25.D0/22528.D0*t8447-13
     #.D0/5632.D0*t8449-39.D0/5632.D0*t19207
      t25384 = -3.D0/64.D0*t19210+39.D0/704.D0*t19212-27.D0/2816.D0*t845
     #2+35.D0/2816.D0*t8456-35.D0/2816.D0*t8458+93.D0/1408.D0*t8460-39.D
     #0/704.D0*t19218+39.D0/1408.D0*t19220+27.D0/2816.D0*t8462+39.D0/563
     #2.D0*t8464+25.D0/5632.D0*t8467
      t25396 = -27.D0/5632.D0*t8469-35.D0/11264.D0*t8471+3.D0/64.D0*t192
     #25-3.D0/64.D0*t19232+3.D0/256.D0*t19236-3.D0/32.D0*t19238+3.D0/32.
     #D0*t19241-27.D0/11264.D0*t8473+3.D0/64.D0*t19243-3.D0/64.D0*t19246
     #+3.D0/32.D0*t19248
      t25411 = -t19250-9.D0/1408.D0*t24934*t24114-9.D0/352.D0*t8133*t707
     #1+9.D0/352.D0*t8122*t7076+27.D0/22528.D0*t8477+13.D0/2816.D0*t8479
     #+25.D0/2816.D0*t8482-3.D0/64.D0*t19253-3.D0/128.D0*t19256+3.D0/256
     #.D0*t19260-39.D0/1408.D0*t19262
      t25420 = t19264-73.D0/2816.D0*t8488+t19266-73.D0/2816.D0*t8492-35.
     #D0/1408.D0*t8494-27.D0/5632.D0*t8496-39.D0/1408.D0*t8498+93.D0/281
     #6.D0*t8501-t19271+t8506/64+25.D0/2816.D0*t8508
      s1 = 9.D0/11264.D0*t25016*t24027-9.D0/2816.D0*t25075*t23886+9.D0/2
     #816.D0*t8047*t23701+9.D0/5632.D0*t24934*t24021+9.D0/352.D0*t8180*t
     #24002-9.D0/352.D0*t24983*t23779+9.D0/704.D0*t25191*t23776-9.D0/140
     #8.D0*t8122*t7341-9.D0/11264.D0*t24934*t24027+9.D0/5632.D0*t8007*t2
     #3831-9.D0/11264.D0*t24934*t24014+9.D0/1408.D0*t8115*t6939+9.D0/112
     #64.D0*t25016*t24014-9.D0/2816.D0*t8007*t24018-9.D0/5632.D0*t25016*
     #t24021-9.D0/1408.D0*t25024*t23904+9.D0/5632.D0*t8007*t24024+t24937
     #+t24922+t24910+t24897+t25420+t25411+t25396-41.D0/704.D0*t8238-25.D
     #0/1408.D0*t8235+t25358+73.D0/704.D0*t8165+t25345+t25384+t25319+t25
     #370+t25309+t25291
      wI150pc1 = s1+t25266+t25253+t25237+t25220-t8376/704-39.D0/1408.D0*
     #t8378+25.D0/5632.D0*t8380-25.D0/5632.D0*t8373+41.D0/2816.D0*t8368-
     #t8371/128+t25180+27.D0/11264.D0*t8365+t25154+t25132+t25108+t25082-
     #t19011+t25064-t19019-35.D0/1408.D0*t8382-t8384/176+25.D0/2816.D0*t
     #8386+t19014+t19017+t25044+t25006+t8017/704+39.D0/5632.D0*t8013+35.
     #D0/5632.D0*t8011-13.D0/5632.D0*t8005-t8002/1408+t24975+t24952-93.D
     #0/1408.D0*t7996
      t25443 = -t8525/96-t8529/768+t8532/192+t8536/192-t8538/384-t8542/7
     #68-t8546/48-t8556/96+t8560/192-t8565/96+t8571/96+t8576/48-t8581/19
     #2+t8585/192+t8589/96-t8591/96-t8593/96-t8595/48
      t25462 = -t8598/24-t8600/48-t8603/48-t8605/96-t8610/384+t8612/96-t
     #8614/48+t8616/96+t8623/96-t8626/96+t8628/48+t8630/48+t8634/192+t86
     #36/24-t8641/192+t8646/96-t8648/192-t8651/384
      t25482 = -t8655/768+t8659/192+t8663/384-t8667/384+t8671/384+t8697/
     #96+t8700/48-t8705/768+t8717/96+t8719/192+t8723/192-t8727/96+t8734/
     #96+t8736/192-t8738/384-t8746/192-t8755/192-t8760/48
      t25501 = -t8763/24-t8765/96-t8771/48-t8779/48-t8781/96+t8783/48+t8
     #785/96+t8789/192-t8791/384+t8797/192+t8808/96+t8810/96+t8814/48-t8
     #816/192+t8818/96+t8821/48+t8823/96+t8826/96
      wI151pc1 = t25443+t25462+t25482+t25501
      t25518 = 7.D0/24.D0*t9440-t9446/96+7.D0/192.D0*t9452-7.D0/24.D0*t9
     #458+7.D0/48.D0*t9462-t9467/12+t9470/24-t9476/8-t9484/12-t9487/6+t9
     #490/12+t9495/24+t9497/12+t9500/12-t9502/24
      t25534 = -7.D0/48.D0*t9506-t9511/6-t9513/48-t9518/48+t9520/96-t952
     #2/24+t9525/12-t9529/8-t9531/12-t9533/8+t9536/12+t9538/4+t9540/24-t
     #9542/24-t9550/4
      t25551 = 7.D0/24.D0*t9555-7.D0/48.D0*t9557+t9560/8-t9565/24-t9568/
     #48-t9571/16-7.D0/48.D0*t9573-t9575/24+t9577/32+t9579/24-t9582/96+t
     #9584/6-7.D0/96.D0*t9587+7.D0/192.D0*t9589+t9591/12
      t25567 = 7.D0/48.D0*t9596-t9599/8-t9602/12+t9605/32+t9607/24+t9609
     #/96+t9614/8+t9619/24-t9622/24+t9625/4-7.D0/48.D0*t9628+t9630/24-t9
     #632/12+t9634/48+t9636/12
      wI152pc1 = t25518+t25534+t25551+t25567
      t25569 = t23573*t8517
      t25574 = t23573*t8573
      t25577 = t23573*t8578
      t25580 = t23573*t24481
      t25583 = t23573*t24239
      t25588 = t23573*t8562
      t25593 = t23573*t8684
      t25598 = t23573*t177
      t25599 = t8901*t62
      t25600 = t25598*t25599
      t25603 = t25569*t24471/40+t25569*t24380/40+t25574*t24327/40+t25577
     #*t8580/10-t25580*t24475/20-3.D0/160.D0*t25583*t24300+t25580*t24478
     #/40+3.D0/20.D0*t25588*t24623-3.D0/40.D0*t25588*t24626-3.D0/20.D0*t
     #25593*t24571+3.D0/40.D0*t25593*t24400+3.D0/20.D0*t25600*t24369
      t25604 = t23573*t8554
      t25631 = -t25604*t24372/10+3.D0/40.D0*t25577*t24626-3.D0/40.D0*t25
     #600*t24175-3.D0/20.D0*t25577*t24623+t25604*t24661/20+t25574*t24375
     #/5-t25574*t24413/10+t25583*t4751/5-t25583*t4778/10+t25574*t3812/10
     #-t25574*t8570/20-t25574*t8597/10+t25574*t8575/20
      t25637 = t23573*t24567
      t25642 = t23573*t24208
      t25647 = t23573*cw2
      t25648 = t25647*t25599
      t25653 = t23573*t24352
      t25658 = t23573*t8901
      t25665 = t25588*t8602/4-t25588*t8580/8-t25637*t4751/4+t25637*t4778
     #/8+3.D0/20.D0*t25642*t24571-3.D0/40.D0*t25642*t24400-3.D0/20.D0*t2
     #5648*t24369+3.D0/40.D0*t25648*t24175+t25653*t24467/20-t25653*t2447
     #1/40+t25658*t24475/20-t25658*t24478/40+t25580*t24593/20
      t25668 = t23573*t8902
      t25671 = t23573*t5930
      t25674 = t23573*t5299
      t25679 = t23573*t24243
      t25684 = t23573*t8644
      t25697 = -t25580*t24596/40-t25668*t24679/20-t25671*t24416/40+3.D0/
     #40.D0*t25674*t24584-3.D0/80.D0*t25674*t24422+9.D0/320.D0*t25679*t2
     #4236+t25668*t24448/40-t25684*t24576/20+t25684*t24410/40-9.D0/80.D0
     #*t25604*t24375+9.D0/160.D0*t25604*t24413+t25671*t24581/20-t25574*t
     #24249/20
      t25708 = t23573*t24218
      t25711 = t23573*t24273
      t25716 = t23573*t8725
      t25723 = t23573*t175
      t25726 = t23573*t3388
      t25729 = t25604*t24253/40+t25684*t24236/80+3.D0/320.D0*t25674*t247
     #80+t25653*t24527/160-t25708*t24693/640-3.D0/160.D0*t25711*t24253-3
     #.D0/80.D0*t25637*t24236-t25716*t24236/80+3.D0/80.D0*t25583*t24236+
     #3.D0/160.D0*t25588*t24151-t25723*t8559/10+t25726*t8559/8
      t25732 = t23573*t8548
      t25739 = t23573*t15618
      t25754 = t23573*t24364
      t25759 = t25671*t24830/160+t25732*t8551/20-t25593*t8805/20-3.D0/16
     #0.D0*t25577*t24151+t25739*t24699/40-t25711*t24361/20+t25711*t24172
     #/40+t25711*t24369/20-t25711*t24175/40-3.D0/40.D0*t25653*t24354+3.D
     #0/80.D0*t25653*t24836+3.D0/40.D0*t25754*t24366-3.D0/80.D0*t25754*t
     #24839
      t25761 = t23573*t8518
      t25776 = t23573*t8980
      t25781 = t23573*t24767
      t25786 = t23573*t24170
      t25789 = t23573*t168
      t25792 = t25761*t24687/10-t25761*t24842/20-t25761*t24366/10+t25761
     #*t24839/20+t25604*t24361/20-t25604*t24172/40-t25569*t24467/20-t257
     #76*t24764/20+t25776*t24857/40-9.D0/80.D0*t25781*t24764+9.D0/160.D0
     #*t25781*t24857-t25786*t24361/5-t25789*t24508/40
      t25811 = t23573*t24178
      t25822 = -t25569*t24189/20-t25598*t8520/20+t23573*t54*t24514/40-t2
     #5569*t24527/160+t23573*t1812*t8520/20+t25786*t24172/10+t25786*t243
     #69/5-t25786*t24175/10+9.D0/80.D0*t25811*t24361-t25723*t24544/10+t2
     #5577*t8564/10-t25708*t24747/160-9.D0/160.D0*t25811*t24172
      t25826 = t23573*t24182
      t25829 = t23573*t24267
      t25832 = t23573*t1067
      t25837 = t23573*t5274
      t25848 = t23573*t8679
      t25855 = t25826*t24201/320-t25829*t24616/160+t25832*t8622/20+t2570
     #8*t24554/160+t25837*t24653/160-t25577*t8788/10-t25723*t8545/20+t25
     #574*t8820/20-t25588*t8564/8-t25848*t8711/80-t25826*t24205/640-t256
     #04*t24296/80
      t25860 = t23573*t8721
      t25881 = t23573*t1387*t8554*t167
      t25886 = t25588*t8741/80+t25588*t8788/8+t25860*t8662/16-t25574*t86
     #62/20-t25588*t8708/40+t25653*t24865/160+t25723*t8759/5-t25723*t869
     #6/10-t25726*t8759/4+t25726*t8696/8+t25832*t8699/20-t25881*t24270/8
     #0+t25881*t24314/160
      t25914 = -t25653*t24862/80-9.D0/160.D0*t25811*t24337-t25860*t8704/
     #32-t25588*t8681/20-t25832*t8762/10+3.D0/40.D0*t25642*t24215-t25658
     #*t24490/40-t25580*t24397/40+t25658*t24397/40-9.D0/640.D0*t25604*t2
     #4334+t25574*t8704/40-t25708*t24496/80+t25708*t24827/160
      t25919 = t23573*t8744
      t25942 = t25848*t8708/40-t25848*t8741/80+t25919*t8658/40-t25826*t2
     #4185/640+t25723*t8654/40+t25588*t8711/80-3.D0/40.D0*t25642*t24211-
     #t25826*t24430/80+t25826*t24459/160+t25574*t24334/40-t25604*t24343/
     #80-t25593*t8773/10+t25593*t8688/20
      t25949 = t23573*t24268
      t25962 = t23573*t8674
      t25971 = t25593*t8751/10-t25593*t8676/20-t25949*t24593/80+t25949*t
     #24596/160+t25761*t24868/80-t25761*t24871/160-t25848*t8748/10+t2584
     #8*t8681/20+t25962*t8773/10-t25962*t8688/20-t25962*t8751/10-t25726*
     #t8654/32
      t25978 = t23573*t5936
      t25999 = -3.D0/40.D0*t25588*t24167+3.D0/40.D0*t25593*t24211-3.D0/4
     #0.D0*t25593*t24215+t25978*t24140/160+3.D0/40.D0*t25577*t24167-t256
     #84*t24346/160-t25708*t24358/640-t25919*t8666/80-t25832*t8640/80+9.
     #D0/160.D0*t25604*t24704-t25671*t24514/40-3.D0/80.D0*t25674*t24499+
     #t25708*t24656/320
      t26029 = t25949*t24397/160-t23573*t8673*t8551/20+t25962*t8805/20+t
     #25588*t8748/10+t25708*t24221/80-t25684*t24300/160+t23573*s*t24499/
     #20+t25962*t8676/20+t25674*t24711/40-t25711*t24714/40+3.D0/80.D0*t2
     #5653*t24189-t25708*t24287/160-3.D0/80.D0*t25653*t24730
      t26056 = t25761*t24739/20+t25949*t24475/80-t25949*t24478/160+3.D0/
     #320.D0*t25711*t24296+3.D0/160.D0*t25637*t24300+t25716*t24300/160-t
     #25604*t24337/40+t25604*t24714/40+9.D0/160.D0*t25811*t24714+t25577*
     #t24742/10-t25786*t24714/10-9.D0/640.D0*t25679*t24300-t25577*t8602/
     #5
      wI153pc1 = t25603+t25631+t25665+t25697+t25729+t25759+t25792+t25822
     #+t25855+t25886+t25914+t25942+t25971+t25999+t26029+t26056
      t26065 = -5.D0/64.D0*t9645-t19701+t19702+11.D0/8064.D0*t9668-t1970
     #4+t19705+3.D0/448.D0*t9685-t19707-t19708+t19709+t19710+t19711-t197
     #12-t19713-3.D0/896.D0*t9721-17.D0/8064.D0*t9725-t19717+t19718+t197
     #19
      t26071 = t3409*wC4
      t26072 = t23711*t26071
      t26075 = t7994*t26071
      t26078 = t9647*t18038
      t26084 = t6953*wC4
      t26085 = t18177*t26084
      t26090 = t4892*t9743
      t26091 = t1953*wC4
      t26092 = t6943*t26091
      t26095 = t7137*t26091
      t26098 = t1865*wC4
      t26099 = t7137*t26098
      t26102 = t13845*t4197
      t26104 = t23505*wC4
      t26105 = t23653*t842*t26104
      t26108 = -53.D0/2304.D0*t9741-11.D0/64512.D0*t9745+17.D0/8064.D0*t
     #9748+3.D0/896.D0*t9750-t19724-t19725+53.D0/1152.D0*t9759-3.D0/224.
     #D0*t9691*t26072+3.D0/448.D0*t9718*t26075-3.D0/224.D0*t26078*t26075
     #+51.D0/896.D0*t9718*t7137*t26071-3.D0/112.D0*t9731*t26085+107.D0/8
     #064.D0*t9763-t19728-281.D0/16128.D0*t9775+3.D0/1792.D0*t26090*t260
     #92-3.D0/1792.D0*t10261*t26095-3.D0/896.D0*t11546*t26099+3.D0/1792.
     #D0*t26102*t26105
      t26110 = t4468*t9743
      t26111 = t6943*t26098
      t26114 = t6949*wC4
      t26115 = t18177*t26114
      t26120 = t3970*wC4
      t26121 = t6938*t26120
      t26126 = t3692*wC4
      t26127 = t6938*t26126
      t26130 = t10065*t23658
      t26131 = t9715*t26126
      t26134 = t9544*t23658
      t26137 = t9715*t26120
      t26140 = t5590*t8639
      t26141 = t8093*wC4
      t26142 = t6953*t26141
      t26145 = t6987*t26098
      t26148 = t842*wC4
      t26149 = t6943*t26148
      t26152 = t5590*t7543
      t26155 = t18038*t26114
      t26165 = t19732+3.D0/1792.D0*t26110*t26111+3.D0/224.D0*t9731*t2611
     #5-3.D0/56.D0*t9642*t26115+3.D0/224.D0*t9799*t26121+t19733+3.D0/28.
     #D0*t9642*t26085-3.D0/112.D0*t9799*t26127-3.D0/112.D0*t26130*t26131
     #+3.D0/112.D0*t26134*t26131+3.D0/224.D0*t26130*t26137-3.D0/896.D0*t
     #26140*t26142+3.D0/1792.D0*t10261*t26145-3.D0/448.D0*t9951*t26149+3
     #.D0/1792.D0*t26152*t26095+3.D0/448.D0*t9772*t26155+9.D0/224.D0*t98
     #52*t26155-51.D0/896.D0*t9761*t7137*t26148-3.D0/1792.D0*t9833*t2610
     #5
      t26166 = t6938*t26148
      t26171 = t5590*t23907
      t26174 = t449*wC4
      t26175 = t6943*t26174
      t26180 = t3752*t3426
      t26181 = t5590*t26180
      t26184 = t7994*wC4
      t26185 = t18038*t26184
      t26193 = 3.D0/224.D0*t9857*t26166-3.D0/224.D0*t9666*t26166-3.D0/89
     #6.D0*t26171*t26092-3.D0/1792.D0*t9673*t26175+t19734+11.D0/4032.D0*
     #t9801-11.D0/4032.D0*t9804+t19737+t19738-t19739-t19740-t19741+t1974
     #2+t19743+t19744-3.D0/7168.D0*t26181*t26111+3.D0/7168.D0*t26171*t26
     #185-3.D0/7168.D0*t9772*t3429*t26184+3.D0/1792.D0*t20324*t26185+t19
     #745
      t26200 = -t19746+t19748-t19749+t19750-t19751+t19752-t19753+t19754+
     #t19755-t19756-11.D0/8064.D0*t9873-59.D0/4032.D0*t9883-t19759+227.D
     #0/8064.D0*t9889-t19761+t19762+t9895/42-t19767+t19768
      t26208 = t23711*wC4
      t26209 = t3579*t26208
      t26215 = t6987*t26141
      t26219 = t3429*t7137*wC4
      t26223 = t3429*t6987*wC4
      t26227 = -3.D0/224.D0*t26134*t26137+3.D0/896.D0*t9678*t26142+5.D0/
     #288.D0*t9903-t9906/56+23.D0/10752.D0*t9908-t19772-3.D0/448.D0*t977
     #2*t26209-3.D0/1792.D0*t11749*t18177*t26184+3.D0/448.D0*t9703*t2621
     #5-9.D0/896.D0*t9852*t26219-3.D0/896.D0*t9731*t26223-t19773+t19774+
     #t9919/56+t19776+t19777-t19778+t19779-t19781
      t26236 = -5.D0/504.D0*t9935+t9937/56+t19784-t19785+5.D0/504.D0*t99
     #43+t19787-5.D0/504.D0*t9947+t19789+281.D0/129024.D0*t9954-t19791-t
     #19792+t19793-t19794+t19795-t19798+t9977/112+t19800-t9982/112-t1980
     #2
      t26240 = t3428*t1811
      t26241 = t6937*wC4
      t26242 = t26240*t26241
      t26245 = t3429*t26208
      t26248 = t7103*wC4
      t26249 = t26240*t26248
      t26252 = t6943*t26104
      t26259 = t19803-23.D0/5376.D0*t9991-5.D0/1008.D0*t9993-t19806+5.D0
     #/1008.D0*t9997+t19808+t19809-t19810+3.D0/896.D0*t11752*t26242+3.D0
     #/7168.D0*t9673*t26245-3.D0/224.D0*t9683*t26249+3.D0/224.D0*t9888*t
     #26252-3.D0/448.D0*t9738*t26072-t19811+t19812+t19814-3.D0/448.D0*t1
     #0016+t10019/112+t19817-t19818
      t26265 = t19819-t19820-t19821+t19822+t19823+t19824-t10049/112-t198
     #26-t19827+t19828+295.D0/4032.D0*t10067-t19834-t19835+t19836-t19837
     #+t19838-t19839+t19840+t19841
      t26267 = t4892*t9887
      t26268 = t6943*t26141
      t26273 = 661.D0/16128.D0*t10114-3.D0/448.D0*t26267*t26268+t19843+t
     #19844-t19845-t19846+t19848-t19849-t19850+t19851+t19852+t19853+t198
     #54-t19855-295.D0/4032.D0*t10164+t19857-t19858+t19859+5.D0/1008.D0*
     #t10176
      t26279 = -t19861+t10184/224-t19865-t19866+t19867-t19868+3.D0/1792.
     #D0*t10201-11.D0/16128.D0*t10203-t19871-t19872-t19873+t19874-t19875
     #-t19876+t19877+t19878-661.D0/16128.D0*t10234-t19881-t19882
      t26284 = -t19883-t19884+t19885+11.D0/4032.D0*t10250-t19887-t19888+
     #t19889+t19890-t19891+t19892-t19893+t19894-227.D0/16128.D0*t10275-5
     #9.D0/32256.D0*t10279+t19900+t19901+t19902-3.D0/448.D0*t10287-t1990
     #4-t19905
      t26291 = t3429*t26241
      t26294 = t18177*t26241
      t26297 = t6987*t26071
      t26302 = t6943*t26071
      t26307 = CsTh*wC4
      t26308 = t6943*t26307
      t26314 = t18038*t26241
      t26317 = t103*wC4
      t26318 = t6987*t26317
      t26321 = t6943*t26317
      t26327 = -t19906-t19907+t19908-t19909+3.D0/224.D0*t10066*t26242-3.
     #D0/224.D0*t10249*t26242-3.D0/448.D0*t10180*t26291+27.D0/896.D0*t10
     #492*t26294-27.D0/896.D0*t10113*t26297+3.D0/448.D0*t10113*t26072-3.
     #D0/112.D0*t10508*t26302+3.D0/112.D0*t10180*t26149-3.D0/56.D0*t1018
     #0*t26308-3.D0/448.D0*t13724*t3606*t26241-51.D0/896.D0*t10508*t2631
     #4-3.D0/448.D0*t10113*t26318-3.D0/448.D0*t10508*t26321-3.D0/448.D0*
     #t10492*t26175-t10308/112
      t26342 = -t19911+t19912-t19913-t19915-t19916+9.D0/3584.D0*t10327+t
     #10333/128-t10337/896-t19920+3.D0/448.D0*t10346+9.D0/896.D0*t10350+
     #17.D0/5376.D0*t10354+t10357/896-t10361/672+t10365/168+3.D0/1792.D0
     #*t10369+3.D0/1792.D0*t10373+t10377/384+t10381/768-t10385/896
      t26349 = t10065*t8639
      t26363 = t10387/448+27.D0/896.D0*t9738*t26297+3.D0/56.D0*t9857*t26
     #297-3.D0/224.D0*t26349*t26142+3.D0/1792.D0*t26110*t26321-t10389/67
     #2-17.D0/1344.D0*t10391+t19937-t10396/64-t19939+t19940+t10407/64-3.
     #D0/448.D0*t10409+t10412/224+t19944+3.D0/224.D0*t10417-t10419/32+t1
     #0422/32-t19949
      t26378 = 17.D0/672.D0*t10427-9.D0/896.D0*t10430-t10432/84-t10436/2
     #24+t10440/448-5.D0/1008.D0*t10442+5.D0/1008.D0*t10449+t19957+3.D0/
     #112.D0*t10180*t6943*wC4-t19958-t19959+107.D0/16128.D0*t10461+59.D0
     #/8064.D0*t10463-3.D0/448.D0*t10465+3.D0/1792.D0*t10467+t19969+t199
     #70+t19971+t19972-t19973
      t26383 = t10065*t23644
      t26384 = t19205*t26120
      t26387 = t5598*t23907
      t26388 = t1803*wC4
      t26389 = t6943*t26388
      t26392 = 9.D0/224.D0*t26383*t26384-3.D0/3584.D0*t26387*t26389-t199
     #74-t19975-t19976-t19977+t19978+t19979-t19980+t19981-t19983+t19984+
     #t19985+t19986+t19987+t19988-t19989+t19990-t19991
      t26406 = t19992+t19993-5.D0/144.D0*t10550+107.D0/16128.D0*t10554-1
     #07.D0/8064.D0*t10558+t19997-t20000+t20001+5.D0/1008.D0*t10567+53.D
     #0/2304.D0*t10569-t20004-5.D0/288.D0*t10573-107.D0/16128.D0*t10576-
     #53.D0/2304.D0*t10579-t10581/112-t10585/896+t10587/224-9.D0/1792.D0
     #*t10590-281.D0/32256.D0*t10593
      t26422 = -t20013-t20014+t20016-107.D0/64512.D0*t10604+3.D0/896.D0*
     #t10607-t10610/896-t10614/112+t10617/64-t10619/32-t20023-55.D0/5376
     #.D0*t10627-661.D0/8064.D0*t10629-661.D0/16128.D0*t10631+3.D0/224.D
     #0*t10633-59.D0/8064.D0*t10635-t10637/224+9.D0/1792.D0*t10639-3.D0/
     #448.D0*t10641+t20035
      t26439 = -23.D0/1792.D0*t10645+t20037+661.D0/8064.D0*t10649+59.D0/
     #4032.D0*t10651-23.D0/896.D0*t10653-t10657/448+55.D0/2688.D0*t10660
     #+227.D0/8064.D0*t10662+t10666/84+t10669/112-t20046-9.D0/896.D0*t10
     #675+t10678/64-t10680/384-3.D0/448.D0*t10682+11.D0/16128.D0*t10684-
     #t20053-t20054+227.D0/16128.D0*t10694-t10697/224
      t26457 = t10701/224+t20058+55.D0/2688.D0*t10706-55.D0/2688.D0*t107
     #08-t10710/256+227.D0/64512.D0*t10713-t10715/384-t10717/84-5.D0/128
     #.D0*t10719-t20068+5.D0/1152.D0*t10725-17.D0/32256.D0*t10729+17.D0/
     #672.D0*t10731-17.D0/672.D0*t10733+t20073+59.D0/16128.D0*t10737-t20
     #075+661.D0/64512.D0*t10744-23.D0/7168.D0*t10746
      t26467 = -t20078-3.D0/3584.D0*t10754-5.D0/512.D0*t10756+3.D0/224.D
     #0*t10758-11.D0/2016.D0*t10760-t20084-t20085+t20086-t20087+t20088+1
     #1.D0/4032.D0*t10772+t20090+t20091-11.D0/16128.D0*t10782-9.D0/7168.
     #D0*t10785-t10789/1792+55.D0/21504.D0*t10792+t20096+t20101
      t26480 = t23925*t26091
      t26483 = -t20102+t20103+53.D0/9216.D0*t10802-t20105+5.D0/128.D0*t1
     #0807-5.D0/128.D0*t10809-t10816/112-t10818/3072+t20110-t10822/84+t1
     #0826/768-55.D0/5376.D0*t10829-3.D0/224.D0*t10831-t20116+3.D0/896.D
     #0*t10836-t20118-t20119-t10848/896+3.D0/1792.D0*t26181*t26480
      t26484 = t26240*t26114
      t26487 = t23653*t26091
      t26491 = t23653*CsTh*t26104
      t26494 = t1224*wC4
      t26495 = t6987*t26494
      t26513 = -3.D0/896.D0*t11749*t26484+3.D0/1792.D0*t10169*t26487+3.D
     #0/896.D0*t9833*t26491-3.D0/896.D0*t10261*t26495+t10851/896+t10853/
     #448-23.D0/10752.D0*t10858+53.D0/2304.D0*t10861+3.D0/448.D0*t10863-
     #23.D0/43008.D0*t10865-281.D0/32256.D0*t10867-t10870/768+9.D0/1792.
     #D0*t10873-9.D0/1792.D0*t10876+9.D0/1792.D0*t10879+t10884/768+t1088
     #6/168-t10888/168-53.D0/2304.D0*t10890-t20138
      t26531 = -3.D0/896.D0*t10894+55.D0/5376.D0*t10897-t20141-t10905/76
     #8-t10909/1792-55.D0/5376.D0*t10911+t10913/448-t10915/448-17.D0/134
     #4.D0*t10918-t10921/224-17.D0/1344.D0*t10923+17.D0/1344.D0*t10925+t
     #20152-59.D0/32256.D0*t10929+t20154+t20155-661.D0/32256.D0*t10940+t
     #10943/448+t20158
      t26538 = t20159+t10952/3584-t10955/896-t20162-t20163-t20167-t20168
     #-t20169+t20170+t20171-t20172+t20173-3.D0/448.D0*t10984+t10986/64-t
     #10988/64-t20177-t20178+t20179+3.D0/448.D0*t10999
      t26547 = t7137*t26388
      t26550 = t18038*t26084
      t26555 = t6938*t26307
      t26560 = t11001/3584+t11003/168+t11005/224+t20185-17.D0/8064.D0*t1
     #1010-t20187-t20188-t20189+t20190-5.D0/64.D0*t11028+3.D0/224.D0*t99
     #51*t26308-3.D0/896.D0*t26152*t26547-3.D0/224.D0*t9772*t26550-9.D0/
     #112.D0*t9852*t26550-3.D0/112.D0*t9857*t26555-t20192+5.D0/144.D0*t1
     #1033-t20194+17.D0/5376.D0*t11046
      t26563 = t7137*t26141
      t26576 = t7137*t26174
      t26579 = -t20196+23.D0/3584.D0*t11050-t20200+t20201+17.D0/8064.D0*
     #t11066-t20203-t20204-t20205-t20206+3.D0/448.D0*t9857*t26563+5.D0/6
     #4.D0*t11082+t20208+t20209+23.D0/1792.D0*t11089+3.D0/224.D0*t9738*t
     #26563-3.D0/448.D0*t9861*t26563+3.D0/448.D0*t9772*t26302+3.D0/1792.
     #D0*t26090*t26175-3.D0/1792.D0*t10261*t26576+t20211
      t26589 = t23653*t26388
      t26599 = -t20212-295.D0/2016.D0*t11096+t20214+t20216+23.D0/5376.D0
     #*t11109-t20218-t20220+t20221+t20222+t20223-27.D0/448.D0*t9718*t698
     #7*t26307+3.D0/112.D0*t9905*t26127+3.D0/224.D0*t9814*t26589+3.D0/89
     #6.D0*t11099*t26495+3.D0/896.D0*t9673*t26389+t20224-t20225-t20226+5
     #1.D0/448.D0*t9761*t7137*t26307
      t26611 = t10065*t26180
      t26614 = t4460*t9743
      t26617 = t9544*t23907
      t26628 = 3.D0/1792.D0*t26181*t26209+3.D0/7168.D0*t26387*t26092+51.
     #D0/3584.D0*t26110*t26155-t20227-t20228-53.D0/1152.D0*t11141+53.D0/
     #1152.D0*t11143+3.D0/448.D0*t11145-t20238-t11150/1792+55.D0/21504.D
     #0*t11153-9.D0/896.D0*t26611*t26155-3.D0/448.D0*t26614*t26149-3.D0/
     #896.D0*t26617*t26115+3.D0/896.D0*t11752*t26484-27.D0/3584.D0*t2609
     #0*t26115+51.D0/3584.D0*t26110*t26314-t20241+3.D0/896.D0*t11157-661
     #.D0/16128.D0*t11159
      t26630 = t1803*t26071
      t26633 = t6953*t26104
      t26645 = t20244-t20245-3.D0/224.D0*t9781*t26630-3.D0/112.D0*t9778*
     #t26633+3.D0/224.D0*t9738*t26633+3.D0/224.D0*t9924*t26142+23.D0/179
     #2.D0*t11167+t20247+t20248-59.D0/8064.D0*t11173+23.D0/10752.D0*t111
     #75+t20251-295.D0/2016.D0*t11181+t20254-151.D0/4032.D0*t11185+t2025
     #6+t20257+t20258-t20259
      t26660 = -t11198/3072-9.D0/7168.D0*t11200-t20262-t20263+t20264-t20
     #265+t20266+t20267+151.D0/2016.D0*t11215-3.D0/1792.D0*t26181*t26155
     #+3.D0/1792.D0*t9951*t3429*t26114+3.D0/224.D0*t20324*t26115+3.D0/22
     #4.D0*t9642*t26223+t20271-t20272+t20273-t20274-t20275-151.D0/2016.D
     #0*t11234+151.D0/2016.D0*t11237
      t26670 = t20278+t20279+t20280+t20281+59.D0/8064.D0*t11251+17.D0/40
     #32.D0*t11253-t20284-t20286+107.D0/16128.D0*t11260-t20288+t20289-t2
     #0290-t20291+t20292-t20293+295.D0/2016.D0*t11279-t20295-t20296-t202
     #97
      t26672 = t23925*t26388
      t26675 = t26240*t26084
      t26692 = t6949*t26317
      t26706 = 227.D0/16128.D0*t11289-3.D0/896.D0*t26181*t26672+3.D0/448
     #.D0*t11749*t26675-3.D0/896.D0*t10169*t26589+27.D0/1792.D0*t26090*t
     #26085-17.D0/4032.D0*t11291+t20300-281.D0/32256.D0*t11296+3.D0/224.
     #D0*t26614*t26308+3.D0/448.D0*t26617*t26085-3.D0/448.D0*t11752*t266
     #75+3.D0/448.D0*t9857*t26145-3.D0/448.D0*t9761*t26692+3.D0/896.D0*t
     #26181*t26550-3.D0/896.D0*t9951*t3429*t26084-3.D0/112.D0*t20324*t26
     #085-51.D0/1792.D0*t26110*t26550+9.D0/448.D0*t26611*t26550+t20305
      t26713 = t6949*t26174
      t26724 = 5.D0/288.D0*t11301+t20307-5.D0/288.D0*t11306-t20309-3.D0/
     #448.D0*t11310+t20311-t20312+t20313+t20314-t20315+3.D0/1792.D0*t977
     #2*t26111-3.D0/448.D0*t9718*t26713+3.D0/224.D0*t9653*t26092-3.D0/44
     #8.D0*t9691*t26095+3.D0/448.D0*t26078*t26713+t20316+23.D0/10752.D0*
     #t11325+11.D0/16128.D0*t11327-t20320
      t26729 = t5590*t23644
      t26732 = t9682*t8639
      t26733 = t6949*t26104
      t26736 = t13848*t4197
      t26737 = t7082*t26104
      t26745 = t20322+t20327+59.D0/8064.D0*t11331+t20329-227.D0/16128.D0
     #*t11335+t20331+17.D0/8064.D0*t11341+5.D0/128.D0*t11343-t20334+t203
     #35+3.D0/448.D0*t26729*t26384-3.D0/448.D0*t26732*t26733-3.D0/448.D0
     #*t26736*t26737+3.D0/448.D0*t9913*t26737+661.D0/16128.D0*t11349-5.D
     #0/288.D0*t11351-t20339+t20343-t20346+5.D0/128.D0*t11353
      t26754 = t20348+151.D0/4032.D0*t11358-t20350-t20351-23.D0/5376.D0*
     #t11364+t20353+59.D0/8064.D0*t11368-t20355+t20356+t20357-t20358+t20
     #359+t20361-t11383/84+t20363-23.D0/10752.D0*t11387+t20366-t20369+28
     #1.D0/32256.D0*t11389
      t26755 = t9682*t8653
      t26756 = t6949*t26141
      t26759 = t9647*t18177
      t26760 = t6949*t26071
      t26763 = t6943*t26120
      t26766 = t6987*t26120
      t26772 = 3.D0/448.D0*t26755*t26756-27.D0/896.D0*t26759*t26760+3.D0
     #/448.D0*t9772*t26763+27.D0/896.D0*t9738*t26766-3.D0/112.D0*t9761*t
     #26760+t20372-t20374-t20376+t20377-t20378+t20380+t20381+t20382-t203
     #87+t20388+23.D0/10752.D0*t11404-t20390+t20391-t20392
      t26780 = -11.D0/4032.D0*t11418+t20394+t20395+t20396-661.D0/16128.D
     #0*t11430+t20398-t20399+t20400+t20402-t20403+t20404-t20405+281.D0/1
     #29024.D0*t11449+3.D0/56.D0*t9857*t26766-t20407-t20408-t11456/256-t
     #20410+t20411
      t26786 = t20412+t11465/168+11.D0/32256.D0*t11467-t20415+t20416+t11
     #473/448+9.D0/1792.D0*t11475-t11477/168+t20422-t20423+t20424+t20425
     #-t20426-t20427+t20428+t20429-t20430+t20431+t20432+t20433
      t26794 = t7137*t26120
      t26798 = t11512/768+55.D0/5376.D0*t11515-t20437+t20438+t20439-t204
     #40+t20441+t20442-5.D0/4032.D0*t11534+5.D0/1152.D0*t11536+t20445-t2
     #0446+51.D0/896.D0*t9718*t26794+t20447-t20448+t20449+t20453-17.D0/3
     #2256.D0*t11554+t20455
      t26806 = -55.D0/5376.D0*t11558+t20457+t20458-t20459-t20460+t20461-
     #3.D0/3584.D0*t11574-5.D0/512.D0*t11576-t20464+53.D0/9216.D0*t11584
     #+t20466-t20467-23.D0/43008.D0*t11591-t20470+t20471-t11600/448+t116
     #02/448+t20474+t20475
      t26808 = t19205*t26494
      t26811 = t7176*t26141
      t26826 = t6943*t26494
      t26829 = -t20476-t20477-t20478-t20479-t20480-3.D0/896.D0*t20324*t2
     #6808+3.D0/224.D0*t9861*t26811+3.D0/112.D0*t9691*t26633-3.D0/224.D0
     #*t9718*t26142+3.D0/112.D0*t26078*t26142-t20481-t20482+3.D0/112.D0*
     #t9666*t26555+3.D0/448.D0*t26171*t26389-3.D0/3584.D0*t26171*t26808+
     #3.D0/3584.D0*t26181*t26826+t20483-t20486+t20487
      t26830 = t4468*t9887
      t26831 = t6943*t26126
      t26834 = t7137*t26126
      t26839 = t19205*t26126
      t26843 = t16117*t4197
      t26848 = t7176*t26104
      t26851 = t7183*t26141
      t26854 = -t20488-3.D0/56.D0*t26830*t26831+9.D0/112.D0*t9747*t26834
     #+51.D0/448.D0*t9924*t26834-9.D0/112.D0*t26383*t26839-t20489+t20490
     #+17.D0/1344.D0*t11641+3.D0/224.D0*t26843*t26491+3.D0/224.D0*t9678*
     #t26834-3.D0/224.D0*t9767*t26848+3.D0/224.D0*t9905*t26851-t20492+t2
     #0493+t20494-t20495-t20496+t20497-t20498+t20499
      t26877 = t20500+t20502-t20503+t20504+t20505+t20506-t20507+3.D0/224
     #.D0*t9731*t26294-3.D0/56.D0*t9642*t26294+3.D0/224.D0*t9803*t26249+
     #3.D0/224.D0*t26130*t26249-3.D0/224.D0*t26134*t26249+t20508-17.D0/1
     #344.D0*t11694+t11696/448-3.D0/1792.D0*t11483*t26252+3.D0/1792.D0*t
     #9648*t26072+3.D0/1792.D0*t9691*t26563-3.D0/1792.D0*t9653*t26268
      t26878 = t1953*t26071
      t26893 = t5590*t23938
      t26896 = t7137*t26494
      t26903 = t6953*t26317
      t26911 = t9647*t24107
      t26914 = t7082*t26141
      t26919 = 3.D0/448.D0*t9781*t26878+3.D0/224.D0*t9778*t26733-3.D0/44
     #8.D0*t9738*t26733-3.D0/448.D0*t9924*t26756+3.D0/448.D0*t26349*t267
     #56+3.D0/224.D0*t9913*t26851-3.D0/224.D0*t10161*t26851+3.D0/224.D0*
     #t26893*t26831+3.D0/224.D0*t9747*t26896+t20511+t20512-3.D0/112.D0*t
     #9718*t26896+3.D0/896.D0*t9678*t26896+3.D0/224.D0*t26759*t26903-3.D
     #0/224.D0*t9857*t26495-3.D0/896.D0*t9772*t26826-3.D0/896.D0*t11704+
     #3.D0/224.D0*t26911*t26756-3.D0/448.D0*t9861*t26914-3.D0/224.D0*t96
     #91*t26733
      t26938 = t23653*t26104
      t26943 = 3.D0/448.D0*t9718*t26756-3.D0/224.D0*t26078*t26756+3.D0/2
     #24.D0*t9761*t26903+t11706/336-17.D0/2688.D0*t11713+t11716/896-3.D0
     #/448.D0*t11718-3.D0/896.D0*t11720-55.D0/10752.D0*t11722+t11724/153
     #6+t20526+3.D0/1792.D0*t26140*t26756-3.D0/1792.D0*t9678*t26756+3.D0
     #/448.D0*t9678*t26938-3.D0/448.D0*t26843*t26105-t20527+t20528-t2052
     #9-t20530
      t26948 = t6988*t26141
      t26962 = t20531+t20532-3.D0/448.D0*t11744-t20534+t11750/336+151.D0
     #/16128.D0*t11753+t20538-295.D0/16128.D0*t11757+3.D0/448.D0*t10161*
     #t26948-3.D0/448.D0*t26893*t26763-3.D0/448.D0*t9747*t26099+3.D0/224
     #.D0*t9718*t26099-3.D0/448.D0*t26759*t26692+107.D0/32256.D0*t11767-
     #t20541-t20542+t20543-227.D0/32256.D0*t11778+23.D0/896.D0*t11780-t2
     #0546
      t26975 = t6953*t26174
      t26988 = -t20547+t20548+t20549-t20550-t20551-3.D0/448.D0*t9678*t26
     #794+3.D0/448.D0*t9767*t26737-3.D0/448.D0*t9905*t26948-3.D0/448.D0*
     #t9913*t26948-t20554-3.D0/224.D0*t26078*t26975-3.D0/224.D0*t26729*t
     #26839-3.D0/896.D0*t26110*t26826-3.D0/896.D0*t26090*t26389-23.D0/17
     #92.D0*t11799+t11801/84-t20557-11.D0/16128.D0*t11809-281.D0/64512.D
     #0*t11811
      t26994 = t20560+t20561+t20562+t20563-107.D0/64512.D0*t11822-t20565
     #+11.D0/4032.D0*t11827+t20567+t20568+11.D0/16128.D0*t11833+23.D0/17
     #92.D0*t11835-t11837/84-t20573-t20574+t20575+t20576-t20577-t20578-t
     #20579
      t26997 = -t20580-t20581-t20582-t20583+t20584+t20588-107.D0/8064.D0
     #*t11876+t20590-t20591+t20592-t20593+t20594-t20595-t20596-t20597+t2
     #0598+t20599-t20600-t20601
      t27008 = 281.D0/16128.D0*t11909+5.D0/144.D0*t11911+t20605+t20606+t
     #20607-t20608+t20609-281.D0/32256.D0*t11923-5.D0/288.D0*t11925-t206
     #12+227.D0/64512.D0*t11929+t20614-53.D0/2304.D0*t11934+t20616+t2061
     #7+17.D0/8064.D0*t11942+5.D0/128.D0*t11944-t20622+t11948/84-t11952/
     #224
      t27028 = t19205*t26098
      t27031 = -t20625+t20626+t20627+281.D0/32256.D0*t11960-107.D0/16128
     #.D0*t11963+5.D0/2016.D0*t11965-t20631+t20632-227.D0/16128.D0*t1197
     #1-t20634+t20636-17.D0/4032.D0*t11978+3.D0/896.D0*t10261*t26547+3.D
     #0/448.D0*t11546*t26896-3.D0/896.D0*t26102*t26491-3.D0/112.D0*t2691
     #1*t26142-3.D0/896.D0*t26090*t26185-3.D0/1792.D0*t9678*t26099+3.D0/
     #7168.D0*t26171*t27028
      t27047 = t7137*t26104
      t27064 = t6953*t26071
      t27070 = 3.D0/1792.D0*t20324*t27028-3.D0/448.D0*t9814*t26487-3.D0/
     #1792.D0*t11099*t26145-3.D0/1792.D0*t9673*t26092+27.D0/896.D0*t9718
     #*t6987*t26148-3.D0/224.D0*t9905*t26121+3.D0/224.D0*t9718*t26975-3.
     #D0/448.D0*t9703*t27047-3.D0/448.D0*t5590*t23658*t19205*t26104+3.D0
     #/448.D0*t9913*t27047-3.D0/896.D0*t11749*t26242-3.D0/224.D0*t9678*t
     #26491-t20638+3.D0/896.D0*t9663*t26633-3.D0/896.D0*t9648*t26633+3.D
     #0/56.D0*t9761*t27064-51.D0/448.D0*t9718*t26834+17.D0/16128.D0*t119
     #83-t20640
      t27077 = t23925*t26104
      t27080 = t6987*t26126
      t27087 = t5598*t26180
      t27100 = t20641-5.D0/576.D0*t11990-11.D0/16128.D0*t11992+295.D0/40
     #32.D0*t11994-151.D0/4032.D0*t11996+3.D0/1792.D0*t5598*t23644*t2707
     #7+3.D0/112.D0*t9666*t27080+3.D0/28.D0*t9861*t27080-3.D0/112.D0*t98
     #98*t27080-3.D0/7168.D0*t27087*t26480+11.D0/16128.D0*t11998+23.D0/1
     #792.D0*t12000+3.D0/1792.D0*t12002+5.D0/256.D0*t12004-t12006/84-53.
     #D0/4608.D0*t12012+t20656-3.D0/1792.D0*t9653*t26185+3.D0/1792.D0*t1
     #0214*t26209
      t27130 = 3.D0/448.D0*t9857*t26318+3.D0/1792.D0*t9772*t26321+3.D0/2
     #24.D0*t9653*t26175+23.D0/21504.D0*t12016+281.D0/16128.D0*t12018+3.
     #D0/448.D0*t26130*t18177*t6937*t26141-3.D0/448.D0*t9772*t3623*t2624
     #8-3.D0/1792.D0*t9833*t26938-3.D0/448.D0*t9642*t26185+3.D0/448.D0*t
     #9767*t27047-3.D0/448.D0*t9905*t26215-3.D0/448.D0*t9913*t26215+3.D0
     #/224.D0*t9718*t7137*t26317+t20659-t20660+t20661-t20662+t20663-t206
     #64+295.D0/8064.D0*t12032
      t27157 = t20666+t20667+t20668+3.D0/896.D0*t9778*t26811-3.D0/896.D0
     #*t26729*t8216*t26141+3.D0/3584.D0*t27087*t26672-27.D0/3584.D0*t260
     #90*t26294-3.D0/1792.D0*t9772*t26219-t12040/168+t20671-151.D0/8064.
     #D0*t12045+t20673+3.D0/1792.D0*t9951*t26291-3.D0/1792.D0*t9778*t265
     #63+3.D0/1792.D0*t26729*t26268-3.D0/896.D0*t9653*t26245-23.D0/7168.
     #D0*t12051+t20675+661.D0/64512.D0*t12055
      t27170 = -t20677-t20678-227.D0/16128.D0*t12064-3.D0/224.D0*t9772*t
     #26831-27.D0/448.D0*t9738*t27080-3.D0/28.D0*t9857*t27080+3.D0/224.D
     #0*t26732*t26633+3.D0/224.D0*t26736*t26848+t20680-t20681+t20682-t20
     #683-t20684-t20687+t20688-t20689+t20690-295.D0/16128.D0*t12086-t206
     #92
      t27185 = -t20693-3.D0/224.D0*t9913*t26848-3.D0/224.D0*t26755*t2614
     #2+27.D0/448.D0*t26759*t27064-3.D0/112.D0*t9653*t26389+t12092/336+1
     #51.D0/16128.D0*t12094-t20696+t20697-t20698-t20699-295.D0/4032.D0*t
     #12106-t20701-t20703+t20704+295.D0/4032.D0*t12115+295.D0/4032.D0*t1
     #2117-t20707+t20708
      t27194 = -t20709-t20710-151.D0/4032.D0*t12128+151.D0/4032.D0*t1213
     #0-t20713+t20714+t20715-151.D0/4032.D0*t12138-t12140/42+11.D0/8064.
     #D0*t12142-t20722+t20723+3.D0/224.D0*t9691*t26547-t20724-11.D0/8064
     #.D0*t12153-t20726+t20727-t20728-t20729-t20730
      t27201 = t18038*t26248
      t27219 = -23.D0/896.D0*t12165+t12167/42-t20733-t20734+t20735+t2073
     #7-3.D0/448.D0*t9691*t26576+3.D0/448.D0*t26729*t27201-3.D0/224.D0*t
     #9666*t26766-3.D0/56.D0*t9861*t26766+3.D0/224.D0*t9898*t26766-t2073
     #8-3.D0/1792.D0*t9663*t26733+3.D0/1792.D0*t9648*t26733+3.D0/448.D0*
     #t9678*t26105-t20739+661.D0/8064.D0*t12183-11.D0/64512.D0*t12185+3.
     #D0/224.D0*t12188
      t27232 = t20743+3.D0/112.D0*t26830*t26763-9.D0/224.D0*t9747*t26794
     #-51.D0/896.D0*t9924*t26794-11.D0/2016.D0*t12192+11.D0/2016.D0*t121
     #94-t20746+t20747+t20748-t20749-3.D0/224.D0*t12206+t20751+t20754+3.
     #D0/448.D0*t10113*t26576-t20755-t20756-t20757-t20758-t20759-227.D0/
     #8064.D0*t12229
      t27243 = t18177*t26248
      t27246 = t20761-t20762-t20763-t20764+t20765+t20766-t20767-3.D0/179
     #2.D0*t9778*t26914+3.D0/1792.D0*t26729*t8281*t26141-3.D0/7168.D0*t2
     #7087*t26209-3.D0/224.D0*t9666*t26297+t20768-t20770-t20771-t20772-3
     #.D0/56.D0*t11095*t27243+t20773+t20774+t20775
      t27274 = 3.D0/224.D0*t11236*t27243+3.D0/112.D0*t26830*t26302-9.D0/
     #224.D0*t9642*t27201-51.D0/896.D0*t26267*t27201+9.D0/224.D0*t26383*
     #t27201-59.D0/4032.D0*t12266-3.D0/448.D0*t26729*t27077+3.D0/1792.D0
     #*t10261*t26318+107.D0/16128.D0*t12269+3.D0/1792.D0*t9744*t18038*t2
     #6208-3.D0/448.D0*t9653*t27201+3.D0/448.D0*t9772*t26314+9.D0/224.D0
     #*t9852*t26314+t12271/448-t20779+t20780-5.D0/4032.D0*t12277-t20782+
     #t20783-t20784
      wI154pc1 = t27194+t27274+t27246+t27070+t27031+t26531+t26919+t27232
     #+t26997+t27008+t26994+t26877+t26513+t26854+t26988+t26962+t26829+t2
     #6806+t26798+t26786+t26780+t26772+t26754+t26745+t26724+t26943+t2648
     #3+t26273+t26467+t26706+t26457+t26670+t26439+t26422+t26265+t26392+t
     #26406+t26327+t26342+t26259+t26660+t26363+t26378+t26645+t26108+t260
     #65+t26236+t26227+t26279+t26284+t26200+t26628+t26193+t26599+t26579+
     #t26165+t26560+t26538+t27219+t27185+t27170+t27157+t27130+t27100
      t27281 = t5590*t23679*Mzsq
      t27282 = t12296*t23675
      t27286 = t6932*t5*t62
      t27289 = wV1*t54
      t27291 = t5598*t27289*t3426
      t27292 = t167*t3452
      t27293 = t27292*t23791
      t27297 = t3601*t12296*t449
      t27300 = t6932*t13121
      t27305 = t5590*wV1*t3426
      t27306 = t12466*t23791
      t27310 = t5590*t23667*t3426
      t27315 = t12292*t23685
      t27321 = t5598*t27289*t167
      t27325 = t5598*t23667*t167
      t27331 = -t27281*t27282/16-t27286*t12919/16-t27291*t27293/32-t1239
     #1*t27297/64-t27300*t12296*t1224/16-t27305*t27306/16-t27310*t12912/
     #16+t12300*t13017/16+t12352*t27315/32+t12324*t12924/32-t27321*t1291
     #9/64-t27325*t27282/32-t27300*t12296*t1865/32
      t27332 = t12466*t23685
      t27340 = t5598*t23689*t3426
      t27341 = t27292*t8308
      t27344 = t27292*t23685
      t27347 = t27292*t8169
      t27350 = -t27305*t27332/32+t20791+t20792-t27310*t12924/32+t12300*t
     #12975/32-t27340*t27341/64-t27291*t27344/64-t27340*t27347/128+t2079
     #3-t20794-t20795-t20796-t20797-t20798
      t27352 = -t20799-t20800-t20801-t20802-t20804+t20805+t20806-t20807-
     #t20808+t20809-t20810+t20811-t20812
      t27353 = -t20813-t20814+t20815+t20818+t20819+t20820+t20821-t20822+
     #t20823+t20824+t20825+t20826-t20827-t20828
      t27359 = t12296*t23505
      t27364 = t20829-t20831+t20832+t20833+t20834-t20835+t20836+t20837-t
     #20838-t20839+t6932*Mzsq*t12959/16+t23715*t27359/16+t23934*t27359/3
     #2
      t27372 = t12292*t23791
      t27382 = t5590*wV1*Mzsq
      t27389 = wV1*t3452
      t27390 = t5598*t27389
      t27393 = -t6932*t3426*t12947/32+t20840-t20841+t20842+t12324*t13034
     #/32+t12391*t13017/64+t23706*t27372/32+t12391*t12975/128+t23706*t27
     #315/64+t12352*t12959/16-t27382*t12959/16+t23706*t12959/64+t23706*t
     #12903/32+t27390*t24125/16
      t27395 = t12296*t23775
      t27401 = t12466*t23727
      t27408 = t12292*t23711
      t27411 = t12466*t23711
      t27414 = t27292*t8070
      t27417 = t27292*t23727
      t27425 = t3601*t12296*t103
      t27430 = -t27281*t27395/8-t27286*t12903/8-t5590*t23679*t27401/16-t
     #27286*t12959/16+t12391*t12899/128+t23706*t27408/64-t27305*t27411/3
     #2-t12391*t27414/64-t23706*t27417/32-t12352*t13085/32+t12300*t12899
     #/32-t12304*t27425/128-t27382*t12903/8
      t27431 = t6932*t12466
      t27448 = t5590*t27389
      t27461 = t27431*t15466/16+t12352*t27408/32+t27431*t15469/8+t12352*
     #t27372/16+t12324*t12912/16-t27321*t12903/32-t27325*t27395/16+t1235
     #2*t12919/16+t27448*t23806/16-t27382*t12919/16+t23706*t12919/64+t27
     #390*t23806/32+t12352*t12903/8+t27448*t24125/8
      wI155pc1 = t27331+t27350+t27352+t27353+t27364+t27393+t27430+t27461
      t27464 = t3426*t103
      t27465 = t12549*t27464
      t27476 = t5452*t449
      t27477 = t13093*t27476
      t27482 = t5452*t103
      t27483 = t12466*t27482
      t27486 = t27292*t23775
      t27489 = t13093*t27482
      t27492 = t13097*t449
      t27493 = t13093*t27492
      t27498 = t13097*t103
      t27499 = t12466*t27498
      t27502 = t27292*t23675
      t27505 = t13093*t27498
      t27508 = t27292*t27492
      t27511 = t24149*t27465/16+t24280*t13098/96-t24234*t13118/48-7.D0/9
     #6.D0*t24154*t13130-t24165*t13141/24-5.D0/48.D0*t24421*t27477-t2414
     #6*t13167/6+t24252*t27483/12+t24219*t27486/48+t24137*t27489/24-5.D0
     #/96.D0*t24421*t27493-t24146*t13170/12+t24252*t27499/24+t24219*t275
     #02/96+t24137*t27505/48+t24183*t27508/192
      t27513 = t13093*t3426*t449
      t27516 = t27292*t24446
      t27521 = t3426*t1763
      t27522 = t13093*t27521
      t27525 = t12466*t8093
      t27534 = t13093*t23505
      t27539 = t12549*t27482
      t27542 = t12549*t27498
      t27545 = t12292*t1958
      t27548 = t12292*t24470
      t27551 = t27292*t24474
      t27554 = t27292*t24447
      t27560 = -5.D0/96.D0*t24421*t27513+t24615*t27516/24+t24513*t13338/
     #48-t24507*t27522/48+t24188*t27525/12+t24157*t13294/96-t24513*t1311
     #2/48-t24524*t13135/24-5.D0/48.D0*t24507*t27534-t24266*t13155/24+t2
     #4149*t27539/8+t24149*t27542/16-t24353*t27545/24-t24353*t27548/48+t
     #24615*t27551/12+t24615*t27554/24-7.D0/48.D0*t24136*cw2*t13145
      t27562 = t24136*t13093
      t27565 = t13093*t8186
      t27568 = t12541*t23505
      t27573 = t13093*t8257
      t27576 = t12485*t23505
      t27587 = t13093*t24474
      t27590 = t13093*t24447
      t27593 = t12549*t8186
      t27596 = t12549*t8257
      t27603 = t27562*t13060/6+t24421*t27565/12-5.D0/24.D0*t24353*t27568
     #+t27562*t13055/12+t24421*t27573/24-5.D0/48.D0*t24353*t27576-t24146
     #*t13164/12-7.D0/48.D0*t24149*t13151-t24252*t13161/24+t24137*t13108
     #/48-t24698*t27587/12-t24698*t27590/24+t24209*t27593/4+t24209*t2759
     #6/8+t24280*t13094/48-t24234*t13115/24
      t27608 = t27292*t27476
      t27612 = t13093*t3601*t449
      t27617 = t12466*t27464
      t27621 = t13093*t3601*t1763
      t27624 = t13093*t24446
      t27627 = t13093*t8093
      t27637 = t27292*t27521*CsTh
      t27641 = t27292*t27521*t842
      t27645 = t13093*t8129
      t27650 = t13093*t167*t3409*t449
      t27653 = t13093*t27464
      t27656 = -7.D0/48.D0*t24154*t13126-t24165*t13138/12+t24183*t27608/
     #96+t24538*t27612/192-7.D0/96.D0*t24154*t13122+t24252*t27617/24+t24
     #652*t27621/48-t24698*t27624/24+t24421*t27627/24-7.D0/24.D0*t24149*
     #t13148-t24252*t13158/12+t24137*t13105/24+t24219*t27637/24+t24219*t
     #27641/48+t24136*t12308*t27645/8+t24652*t27650/96+t24137*t27653/48
      wI156pc1 = t27511+t27560+t27603+t27656
      wI157pc1 = wV36
      t27659 = t5598*wV4*t3452
      t27663 = t5590*t24978*Mzsq
      t27676 = t5590*wV4*t3426
      t27689 = -t27659*t24125/24+t27663*t27395/24+t5590*t24978*t27401/48
     #-t12902*t27408/96+t12898*t27414/96+t12902*t27417/48+t27676*t27411/
     #96+t13037*t27425/192-t12994*t12293/48-t15491*t27359/48-t25055*t273
     #59/48-25.D0/3456.D0*t12900-25.D0/864.D0*t12904
      t27702 = -t12907/864+37.D0/864.D0*t12909-17.D0/432.D0*t12913-t2100
     #8-25.D0/1728.D0*t12920-t12922/1728-17.D0/864.D0*t12925+37.D0/1728.
     #D0*t12927-t21013-t12931/1728-t21015-t12935/3456+41.D0/864.D0*t1293
     #8-t12915*t12417/48
      t27710 = wV4*t54
      t27712 = t5598*t27710*t167
      t27717 = t5598*wV4*s*t167
      t27729 = -t13067*t12364/24-t13067*t12310/12-t12958*t27372/48+t2771
     #2*t12903/48+t27717*t27395/24+t27676*t27306/48-t13067*t12374/24-t12
     #942/864-t21019-t21020+35.D0/1728.D0*t12951+t12954/432-t12958*t2731
     #5/96
      t27739 = t5598*t24932*t3426
      t27743 = t5598*t27710*t3426
      t27753 = t27712*t12919/96-t12958*t27408/96+t27717*t27282/48+t27676
     #*t27332/96+t27739*t27341/96+t27743*t27344/96-13.D0/1728.D0*t12956+
     #35.D0/864.D0*t12960+41.D0/432.D0*t12965+35.D0/432.D0*t12967-t21028
     #+41.D0/216.D0*t12973-25.D0/3456.D0*t12976+7.D0/1728.D0*t12978
      t27773 = -t12981/1728+t27739*t27347/192-t13040*t24125/24-t13040*t2
     #3806/48-t27659*t23806/48+t27663*t27282/48-t21033+37.D0/432.D0*t129
     #86+41.D0/864.D0*t12989+t12992/108+t12995/216+7.D0/864.D0*t12997+35
     #.D0/864.D0*t12999
      t27792 = -t12915*t12321/24+t27743*t27293/48+t12898*t27297/96-t1290
     #2*t27372/48-t12902*t27315/96+41.D0/432.D0*t13003+7.D0/1728.D0*t130
     #05+t13007/216+37.D0/864.D0*t13009+t13011/216-25.D0/1728.D0*t13013-
     #t21048+35.D0/864.D0*t13018+t13020/432
      t27805 = 37.D0/864.D0*t13023-13.D0/864.D0*t13026+35.D0/1728.D0*t13
     #028-t21054+37.D0/1728.D0*t13032-17.D0/864.D0*t13035-t21057-t13041/
     #432-t13043/3456+7.D0/3456.D0*t13046-t13048/216-13.D0/432.D0*t13050
     #-t13052/432
      t27818 = -17.D0/432.D0*t13056-13.D0/864.D0*t13058-17.D0/216.D0*t13
     #061-17.D0/432.D0*t13063+7.D0/1728.D0*t13065-t21070-t21071+41.D0/43
     #2.D0*t13073+7.D0/3456.D0*t13075-t13077/432-25.D0/1728.D0*t13079-13
     #.D0/864.D0*t13081-13.D0/1728.D0*t13083-t13086/864
      wI158pc1 = t27689+t27702+t27729+t27753+t27773+t27792+t27805+t27818
      wI159pc1 = wI126pc1
      t27837 = t23573*t13093
      t27848 = -t25674*t27513/24+t25829*t27516/48+t25577*t27539/8+t25577
     #*t27542/16-t25653*t27545/24-t25653*t27548/48+t25829*t27551/24+t258
     #29*t27554/48+5.D0/24.D0*t27837*t13060-t25674*t27565/6-t25653*t2756
     #8/6+5.D0/48.D0*t27837*t13055-t25674*t27573/12
      t27877 = -t25653*t27576/12-t25732*t13164/12-t25577*t13151/6-t25739
     #*t27587/12-t25739*t27590/24+t25642*t27593/4+t25642*t27596/8-t25723
     #*t13126/6-t25588*t13138/12+t25826*t27608/48+t25978*t27612/96-t2572
     #3*t13122/12+5.D0/96.D0*t25604*t27617+t25837*t27621/96
      t27906 = -t25739*t27624/24-t25674*t27627/12-t25577*t13148/3+t25708
     #*t27637/48+t25708*t27641/96+t23573*t12308*t27645/8+t25837*t27650/4
     #8-t25671*t27653/24+t25577*t27465/16-t25723*t13130/12-t25588*t13141
     #/24-t25732*t13167/6+5.D0/48.D0*t25604*t27483
      t27935 = -t25674*t27477/12-t25789*t27522/48+5.D0/48.D0*t25569*t275
     #25-t25598*t13135/24-t25789*t27534/12+t25708*t27486/24-t25671*t2748
     #9/12-t25674*t27493/24-t25647*t13145/6-t25732*t13170/12+5.D0/96.D0*
     #t25604*t27499+t25708*t27502/48-t25671*t27505/24+t25826*t27508/96
      wI160pc1 = t27848+t27877+t27906+t27935
      t27951 = -5.D0/2376.D0*t13488-5.D0/594.D0*t13491-5.D0/396.D0*t1349
     #4-35.D0/1188.D0*t13500-31.D0/9504.D0*t13503-245.D0/9504.D0*t13506-
     #t13508/1584+t13512/396-t13514/1584+37.D0/4752.D0*t13517-t13522/66+
     #t13529/66+t13534/33-t13536/33
      t27966 = -t13542/3168-t13545/6336+5.D0/297.D0*t13549-13.D0/396.D0*
     #t13552+85.D0/2376.D0*t13556-5.D0/2376.D0*t13560-t13562/792-5.D0/15
     #84.D0*t13566+73.D0/2376.D0*t13569-5.D0/2376.D0*t13572-5.D0/792.D0*
     #t13575+t13578/99-t13582/1584-5.D0/594.D0*t13585
      t27973 = t211*Mzsq
      t27974 = t27973*t12296
      t27978 = t5590*t12466
      t27983 = t12549*t26241
      t27991 = 35.D0/3168.D0*t13587+37.D0/4752.D0*t13590+5.D0/594.D0*t13
     #593-37.D0/9504.D0*t13596-5.D0/1188.D0*t13598-t27974*t1865*t26071/8
     #8+t27978*t26756/176-t13505*t26756/176+t13584*t27983/88-13.D0/4752.
     #D0*t13601+17.D0/9504.D0*t13604+5.D0/1188.D0*t13606-13.D0/2376.D0*t
     #13609-13.D0/4752.D0*t13611
      t27998 = t5590*t13111
      t27999 = t12296*t26091
      t28002 = t10065*t13134
      t28003 = t12296*t26307
      t28006 = t4892*t13121
      t28014 = t12296*t26388
      t28017 = t5590*t5253
      t28018 = t12296*t26494
      t28023 = -4.D0/297.D0*t13615+5.D0/2376.D0*t13618+5.D0/594.D0*t1362
     #0+2.D0/297.D0*t13623+4.D0/297.D0*t13628+17.D0/4752.D0*t13631+t2799
     #8*t27999/176-t28002*t28003/44-5.D0/88.D0*t28006*t28003-t13540*t124
     #66*t26084/44+t13524*t28003/44-t13580*t28014/88+t28017*t28018/176-t
     #13540*t28018/176
      t28026 = t4892*t13111
      t28029 = t10277*t12549
      t28035 = t4468*t13111
      t28039 = t5598*t5264*t167
      t28042 = t12296*t26148
      t28056 = -t28026*t28018/88+t28029*t26692/176+t13584*t12549*t26084/
     #44+5.D0/88.D0*t28035*t28003-t28039*t28014/176-5.D0/176.D0*t28006*t
     #28042-t13540*t12466*t26114/88+t13524*t28042/88+5.D0/1188.D0*t13633
     #+t13636/297-t13638/1584-t13640/1584-t13642/3168-17.D0/2376.D0*t136
     #44
      t28065 = t12296*t26317
      t28068 = t9682*t12549
      t28077 = t12296*t26174
      t28081 = t10065*t5*t1811
      t28082 = t12296*t26071
      t28088 = t5598*t351*t3452
      t28092 = t13584*t12549*t26114/88+5.D0/176.D0*t28035*t28042+t13580*
     #t12466*t26184/352-t13540*t28065/352+t28068*t26075/88-9.D0/352.D0*t
     #13647+t13649/297+t14044*t28014/176-t28002*t28042/88+t14044*t28077/
     #352-t28081*t28082/44-t28026*t28065/176-t28088*t26733/176-2.D0/297.
     #D0*t13651
      t28108 = t13656/132-5.D0/792.D0*t13658-5.D0/396.D0*t13660+73.D0/47
     #52.D0*t13664+5.D0/1188.D0*t13667+73.D0/1188.D0*t13669-5.D0/1584.D0
     #*t13672-5.D0/792.D0*t13675-317.D0/19008.D0*t13679+37.D0/9504.D0*t1
     #3681-5.D0/1188.D0*t13684-37.D0/9504.D0*t13686+37.D0/4752.D0*t13691
     #+49.D0/2376.D0*t13694
      t28123 = 37.D0/9504.D0*t13696+73.D0/2376.D0*t13698-5.D0/594.D0*t13
     #701+t13704/198-13.D0/792.D0*t13707+17.D0/1188.D0*t13710+5.D0/1188.
     #D0*t13713-37.D0/9504.D0*t13715-17.D0/1188.D0*t13717+t13719/99+2.D0
     #/99.D0*t13722-17.D0/2376.D0*t13726-17.D0/4752.D0*t13730-137.D0/950
     #4.D0*t13733
      t28141 = 5.D0/594.D0*t13735+5.D0/297.D0*t13737-t13740/1584-13.D0/3
     #96.D0*t13743-5.D0/594.D0*t13745-t13747/198-t13749/396+85.D0/1188.D
     #0*t13751-17.D0/2376.D0*t13757+8.D0/297.D0*t13759-9.D0/352.D0*t1376
     #3-9.D0/352.D0*t13765+35.D0/3168.D0*t13767+85.D0/2376.D0*t13769
      t28156 = -t13771/792-9.D0/176.D0*t13773-137.D0/4752.D0*t13775-t137
     #77/396-20.D0/297.D0*t13779-5.D0/1188.D0*t13784-10.D0/297.D0*t13787
     #+8.D0/297.D0*t13790-9.D0/704.D0*t13793+16.D0/297.D0*t13795-17.D0/1
     #188.D0*t13797-40.D0/297.D0*t13799+4.D0/297.D0*t13802+37.D0/4752.D0
     #*t13804
      t28172 = 5.D0/594.D0*t13806+t13808/528-t13811/3168-137.D0/4752.D0*
     #t13815-137.D0/9504.D0*t13817+37.D0/9504.D0*t13819+37.D0/9504.D0*t1
     #3821+37.D0/4752.D0*t13823-t13825/3168+5.D0/297.D0*t13827+t13829/10
     #56+t13831/1056+5.D0/594.D0*t13833-5.D0/1188.D0*t13835
      t28188 = t13838/528+5.D0/594.D0*t13841+5.D0/1188.D0*t13843+37.D0/1
     #9008.D0*t13846+37.D0/4752.D0*t13849-17.D0/1188.D0*t13851-17.D0/594
     #.D0*t13853+t28029*t26903/88-5.D0/1188.D0*t13855-37.D0/4752.D0*t138
     #57-5.D0/594.D0*t13859-17.D0/1188.D0*t13861-17.D0/1188.D0*t13863-17
     #.D0/2376.D0*t13868
      t28205 = -20.D0/297.D0*t13870+5.D0/594.D0*t13873-317.D0/38016.D0*t
     #13877+49.D0/4752.D0*t13879-31.D0/9504.D0*t13881-35.D0/1188.D0*t138
     #83+t13885/66-17.D0/2376.D0*t13887-5.D0/594.D0*t13889-17.D0/4752.D0
     #*t13891+17.D0/2376.D0*t13894+17.D0/1188.D0*t13896+t13902/132+t1390
     #5/792
      t28220 = -t13908/1584+37.D0/9504.D0*t13911-7.D0/396.D0*t13913-17.D
     #0/2376.D0*t13915-5.D0/594.D0*t13917-17.D0/594.D0*t13919-t13922/132
     #+37.D0/9504.D0*t13926+t13928/396-t13930/792+5.D0/297.D0*t13932+35.
     #D0/792.D0*t13934-13.D0/198.D0*t13937-t13939/792
      t28236 = -10.D0/297.D0*t13941-5.D0/1188.D0*t13943-t13945/396-5.D0/
     #594.D0*t13950+35.D0/1584.D0*t13953-13.D0/396.D0*t13956-t13958/1584
     #-5.D0/297.D0*t13960-137.D0/2376.D0*t13962+35.D0/1584.D0*t13964-37.
     #D0/4752.D0*t13967-5.D0/1188.D0*t13969-5.D0/297.D0*t13971-37.D0/237
     #6.D0*t13973
      t28238 = t12296*t26104
      t28243 = t12296*t26098
      t28249 = t12466*t26241
      t28263 = -t5598*t351*t28238/176+t27998*t28014/88-t28026*t28243/176
     #-t13486*t27292*t26208/352-t13540*t28249/88-t13580*t28077/176+17.D0
     #/594.D0*t13976-5.D0/1188.D0*t13978-17.D0/1188.D0*t13985-37.D0/9504
     #.D0*t13987+49.D0/1188.D0*t13989+37.D0/4752.D0*t13991+37.D0/2376.D0
     #*t13993-37.D0/4752.D0*t13995-5.D0/594.D0*t13997
      t28282 = -137.D0/4752.D0*t13999+4.D0/297.D0*t14003+5.D0/1188.D0*t1
     #4006+t14008/594+2.D0/297.D0*t14011+85.D0/1188.D0*t14013-4.D0/297.D
     #0*t14015+37.D0/9504.D0*t14017-17.D0/594.D0*t14020+37.D0/9504.D0*t1
     #4022-13.D0/2376.D0*t14024-4.D0/297.D0*t14027-13.D0/2376.D0*t14029-
     #t14031/792
      t28297 = -t14033/396-17.D0/2376.D0*t14036-5.D0/594.D0*t14040+17.D0
     #/2376.D0*t14042+37.D0/19008.D0*t14045-17.D0/2376.D0*t14047+t14049/
     #99+37.D0/4752.D0*t14051+5.D0/297.D0*t14053+49.D0/2376.D0*t14055+17
     #.D0/4752.D0*t14057+t14059/264+t14061/297+4.D0/297.D0*t14064
      t28308 = t12296*t26126
      t28319 = -2.D0/297.D0*t14066+5.D0/1188.D0*t14070+5.D0/2376.D0*t140
     #72+8.D0/297.D0*t14074+t14076/594+2.D0/297.D0*t14078+t14080/264-t28
     #088*t26633/88-t28081*t28308/22-t13580*t27999/176+t28017*t28243/352
     #-t13540*t28243/352+85.D0/594.D0*t14082+t14085/528
      t28334 = 5.D0/594.D0*t14087+17.D0/9504.D0*t14089+t14095/528-t14097
     #/264+17.D0/4752.D0*t14099+85.D0/1188.D0*t14101-8.D0/297.D0*t14103-
     #13.D0/1188.D0*t14105+t14108/1056-t14112/528+t14115/1056+t14118/105
     #6+17.D0/2376.D0*t14120-t14123/1056
      t28357 = t14126/1056-t13502*t26630/44-t27974*t1224*t26071/44+t2797
     #8*t26142/88-t13505*t26142/88-t13502*t26878/88+37.D0/2376.D0*t14128
     #+37.D0/4752.D0*t14130-t14132/132+t14134/2112+5.D0/297.D0*t14136-t1
     #4139/1056+t14141/2112+5.D0/594.D0*t14143
      t28372 = t14145/528-t14147/264-t14151/528+37.D0/9504.D0*t14153-5.D
     #0/1188.D0*t14155+t14157/528-t14159/264-t14161/528+t14163/2112+t141
     #65/2112-t14167/528+5.D0/1188.D0*t14169+t14171/528+t14173/528
      t28378 = t12296*t26120
      t28391 = -5.D0/1188.D0*t14176+37.D0/19008.D0*t14178-17.D0/1188.D0*
     #t14181+5.D0/594.D0*t14183-t28081*t28378/44+t14044*t27999/352+37.D0
     #/38016.D0*t14186+37.D0/9504.D0*t14189+49.D0/2376.D0*t14191+49.D0/4
     #752.D0*t14193+37.D0/19008.D0*t14196-245.D0/19008.D0*t14198-31.D0/1
     #9008.D0*t14200-35.D0/2376.D0*t14202
      t28402 = t5590*t57
      t28406 = t5590*t32*t3452
      t28413 = -5.D0/792.D0*t14208+8.D0/297.D0*t14210+37.D0/19008.D0*t14
     #213+t14216/33-5.D0/198.D0*t14218-35.D0/594.D0*t14220-31.D0/4752.D0
     #*t14222-245.D0/4752.D0*t14224-20.D0/297.D0*t14226+t14229/66-t28402
     #*t28308/22+t28406*t26633/44+t28068*t26142/44+t13496*t27064/22
      t28417 = t10065*t12549
      t28420 = t4892*t12466
      t28446 = t28417*t27064/22+5.D0/44.D0*t28420*t27064+t13966*t26630/8
     #8-t28402*t28378/44+t28406*t26733/88+t28068*t26756/88+t13496*t26760
     #/44+t28417*t26760/44+5.D0/88.D0*t28420*t26760+t13966*t26878/176+t2
     #7978*t26075/176-t13505*t26075/176+t27998*t28077/176-t14231/264
      t28465 = t14234/1056-t28402*t28082/44-t28039*t27999/352+t5590*t32*
     #t28238/88+t14236/1056-t14238/528-17.D0/2376.D0*t14240+t14242/1056+
     #t14244/1056+37.D0/9504.D0*t14246-17.D0/1188.D0*t14248-7.D0/792.D0*
     #t14251+t14253/132-t14255/1584
      t28472 = t211*t62*t3452
      t28484 = t14257/792-t14261/6336-t14263/1584-t14266/66-5.D0/88.D0*t
     #28472*t26248-t14269/3168-t14271/132-7.D0/396.D0*t14274+t14276/66-t
     #14278/792+t14280/396-t14282/3168-7.D0/198.D0*t14284-t14288/1584
      t28499 = -7.D0/792.D0*t14290+t14293/66-t14295/66-7.D0/396.D0*t1429
     #7-t14299/792+t14302/198-17.D0/1188.D0*t14304+4.D0/297.D0*t14307-9.
     #D0/704.D0*t14310-t14312/3168-t14314/1584+37.D0/19008.D0*t14319+5.D
     #0/594.D0*t14321-35.D0/1188.D0*t14323
      t28516 = -13.D0/792.D0*t14325-t14327/792-t14329/3168-5.D0/1188.D0*
     #t14332+t14334/198-5.D0/396.D0*t14336+t14338/66+5.D0/594.D0*t14340-
     #245.D0/9504.D0*t14343+73.D0/2376.D0*t14345+37.D0/4752.D0*t14347-t1
     #4349/3168-37.D0/19008.D0*t14351-5.D0/2376.D0*t14353
      t28531 = -317.D0/9504.D0*t14355-35.D0/2376.D0*t14357-31.D0/9504.D0
     #*t14359-5.D0/792.D0*t14361+37.D0/4752.D0*t14363+37.D0/9504.D0*t143
     #65-317.D0/19008.D0*t14369-31.D0/19008.D0*t14371+73.D0/4752.D0*t143
     #73-245.D0/19008.D0*t14375+t14377/132+5.D0/1188.D0*t14379-5.D0/396.
     #D0*t14381-245.D0/9504.D0*t14383
      t28555 = 35.D0/1584.D0*t14385+37.D0/2376.D0*t14387-37.D0/4752.D0*t
     #14389+37.D0/38016.D0*t14391+37.D0/19008.D0*t14394-5.D0/88.D0*t2847
     #2*t26760+t13474*t28082/44-t13739*t26072/88-t27973*t3452*t26075/88+
     #5.D0/88.D0*t4892*Mzsq*t28082+t13646*t26072/176-317.D0/19008.D0*t14
     #396-317.D0/38016.D0*t14398-5.D0/594.D0*t14400
      t28581 = -5.D0/297.D0*t14402-37.D0/9504.D0*t14404+10.D0/297.D0*t14
     #406+17.D0/1188.D0*t14409+t10277*t62*t28065/176-t10703*t12509*t2624
     #1/88-5.D0/176.D0*t13662*t27983+t13724*t12292*t26241/88+5.D0/176.D0
     #*t13654*t28249+t13693*t28082/44-5.D0/44.D0*t28472*t27064+37.D0/950
     #4.D0*t14411-5.D0/1188.D0*t14413-17.D0/1188.D0*t14415-37.D0/19008.D
     #0*t14417
      wI161pc1 = t28372+t28141+t28123+t28108+t28092+t28056+t28023+t27951
     #+t27991+t27966+t28581+t28555+t28531+t28516+t28499+t28484+t28357+t2
     #8465+t28334+t28446+t28319+t28297+t28282+t28263+t28413+t28236+t2822
     #0+t28391+t28205+t28188+t28172+t28156
      wI162pc1 = wI136pc1
      t28587 = -t22511+t22512-t22513+t22514-t22524+t22528-t22529-t22530+
     #t22531-t22532-t22533
      t28589 = t22537+t22538-t22539+t22556+t22560+t22561-t22562-t22563-t
     #22564-t22565+t22566
      t28590 = -t22567+t22568-t22569-t22570-t22573-t22574+t22575+t22576+
     #t22577-t22578-t22579
      t28593 = -t22580+t22581-t22582-t22583+t22585-t22586-t22587+t22588+
     #t22589+t22590-t22591
      t28594 = -t22592-t22593+t22616+t22617-t22618+t22619-t22620-t22622+
     #t22623-t22629-t22630
      t28596 = t22631-t22632+t22633+t22634-t22635+t22638-t22645-t16794+t
     #22646+t16801-t22647
      t28597 = t54*wC20
      t28598 = t7103*t28597
      t28601 = t3638*wC20
      t28602 = t7137*t28601
      t28605 = Mzsq*wC20
      t28606 = t7941*t28605
      t28609 = -t22661-t22662+t22686+t22687-t22688+t22689-t16570*t28598/
     #128+t16501*t28602/32-t16501*t28606/32+t22700+t22702
      t28616 = t1072*wC20
      t28620 = t62*wC20
      t28621 = t16541*t28620
      t28624 = t4120*wC20
      t28630 = -t22727+t22728-t22729-t22730+t22743-t22744-t16506*t28602/
     #32-t13490*t6943*t28616/128+t16501*t28621/32+t16871*t6943*t28624/12
     #8-t16513*t28598/32
      t28632 = t3409*wC20
      t28633 = t6987*t28632
      t28638 = t3999*wC20
      t28642 = t22745-t22755+t22756-t22762-t22763+t22764-t22765-t16513*t
     #28633/32-t16501*t28633/32+t16905*t6943*t28638/32-t22773
      t28643 = t6949*t28638
      t28652 = t1070*wC20
      t28653 = t6953*t28652
      t28660 = t6949*t28652
      t28663 = t720*t28601
      t28668 = -t22774-t22775+t22777-t16501*t28643/32+t16506*t28643/32+t
     #13595*t6943*t1071*wC20/256-t16540*t28653/64-t13595*t6943*t1119*wC2
     #0/512+t16540*t28660/128-t17361*t28663/64-t16797*t28663/16
      t28671 = t7176*t28601
      t28674 = t1096*t28601
      t28681 = t7082*t28601
      t28688 = t1077*wC20
      t28689 = t6943*t28688
      t28694 = -t22778+t22779+t16506*t28671/16+t17361*t28674/128-t16496*
     #t28674/16+t16517*t28674/16-t16487*t28681/128+t16797*t28674/32-t165
     #06*t28681/32+t13490*t28689/128+t16570*t28602/128
      t28695 = t6953*t28688
      t28698 = t838*t28605
      t28701 = t6953*t28624
      t28709 = t129*t28620
      t28712 = t6953*t28616
      t28717 = t3692*t28597
      t28720 = t16570*t28695/64-t16496*t28698/8+t16570*t28701/64-t22780+
     #t22781-t16496*t1093*t28632/16+t16501*t28695/16+t16496*t28709/16-t1
     #6501*t28712/16-t16487*t28701/64-t16496*t28717/16
      t28722 = t843*t28605
      t28725 = t6949*t28624
      t28732 = t3609*s*wC20
      t28739 = -t22782+t22783-t22784-t22785+t22786+t16496*t28722/16-t165
     #70*t28725/128+t16496*t5992*t28632/32-t16496*t28732/16+t16513*t2860
     #2/16-t16487*t28602/128
      t28746 = t6953*t28638
      t28753 = -t16710*t28695/256+t16540*t28695/32-t22787+t16513*t28701/
     #16+t16513*t28746/16+t16501*t28746/16-t16506*t28746/16+t22790-t2279
     #1-t22792+t22793
      t28759 = t17053+t22805+t22806-t22807-t22808-t22809+t22810+t22811-t
     #17076-t22815+t22816
      t28761 = t7498*wC20
      t28764 = t6949*t28688
      t28773 = t22817+t22818+t22819-t16480*t28761/32-t22820+t22821+t2282
     #2-t16570*t28764/128+t16496*t28663/8-t16517*t28663/8+t16487*t28671/
     #64
      t28774 = -t22823+t22824+t22826-t22827-t22828+t22829+t22830-t22831+
     #t22832-t22833+t22834
      t28777 = -t22835-t22836+t22839+t22840-t22841+t22848-t22849-t22850+
     #t22851+t22852+t22853
      t28784 = -t22855+t22856-t22857-t22858-t22859+t22860+t22861-t16540*
     #t28712/64-t16570*t28653/64+t16570*t28660/128-t22862
      t28794 = t6949*t28616
      t28799 = t22863-t22864-t16580*t28695/64+t16580*t28712/64+t16710*t2
     #8764/512-t16540*t28764/64+t16540*t28794/128+t22865+t22870-t22871+t
     #16580*t28764/128
      t28802 = t6943*t28652
      t28807 = -t16580*t28794/128-t13595*t28802/512+t13595*t28689/512+t2
     #2872+t22873-t22874-t22875+t22876-t22877-t22878+t22879
      t28819 = t6046*t28620
      t28825 = t16501*t28794/32+t22886-t22887-t22888+t22889+t22890-t2289
     #1+t22894-t22895-t22896-t22897
      t28829 = t3970*t28597
      t28836 = t22898+t16487*t28725/128+t16496*t28829/32+t22899+t22900+t
     #22901-t22902-t16513*t28725/32-t16513*t28643/32+t22903+t22904
      t28837 = t22916-t22917-t17330+t17332-t22918-t22919+t22920-t22924-t
     #17346+t17348-t22925
      t28840 = -t22926+t22927-t22928-t22929-t22930-t22931+t22932-t22934+
     #t22935+t17393-t17395
      t28841 = t22936-t22937-t22941+t22942+t17410+t22943-t17421+t22946+t
     #22947+t22948+t22949
      t28845 = t22950-t22951-t22952-t22953-t22954-t22955+t22957+t22958-t
     #22959+t16480*t28598/32+t22960
      t28860 = t22961+t16480*t28653/16-t16480*t28660/32+t13555*t28802/12
     #8+t16480*t28606/16+t13739*t28802/128+t16480*t28633/32-t16480*t2862
     #1/32+t22962+t22963-t22964
      wI163pc1 = t28841-t22725-t22726+t28837+t28840+t28836+t28825+t28807
     #+t28799+t28774+t28777+t28784+t28773+t28759+t22803-t22804+t28753-t2
     #2802-t22800+t28739+t22798+t22799+t22796+t22797+t22794-t22795+t2872
     #0+t28694+t28668+t28642+t28630+t28609+t22510-t22508-t22509-t22506+t
     #22507+t22479+t22477-t22478+t22474+t22476+t22719+t22720+t28593+t285
     #94+t28596+t28587+t28589+t28590+t22705-t22706+t22885-t22883-t22884-
     #t16496*t28819/32-t16570*t28606/128-t16501*t28764/32-t16540*t28606/
     #64+t16540*t28621/128+t22881+t22882+t28860+t22723-t22724+t28845+t22
     #703-t22704
      t28869 = -t22990-t22991-t22992+t22993-t22994-t22995+t22996-t22997+
     #t22998+t22999-t23000-t23001-t16326*t28638/108+t17675*t28601/108
      t28889 = -t23002-t23003+t16138*t28688/216-t16138*t28616/216+t16125
     #*t28624/216-t16323*t28638/54+t10712*t28761/432-t17493*t28698/216-t
     #17493*t28722/432-t15975*t28732/108-t15972*t28732/216+t23020+t23021
      t28890 = t23022+t23023+t23024+t23025-t23027+t23028+t23029-t23030+t
     #23031+t23032-t23033+t23034+t23035+t23036
      t28893 = t23037-t23038+t23039+t23040-t23044+t23045-t23046+t23047-t
     #23048-t23049+t23050-t23051+t23052
      t28894 = -t23053+t23054+t23055+t23056+t23058-t23059+t23060+t23061+
     #t23062+t23063-t23064-t23065+t23066+t23067
      t28896 = -t23068+t23069+t23070+t23071+t23074+t23075+t23076-t23077+
     #t23078-t23079+t23080-t23081+t23084
      t28923 = -t16044*t28638/108+t23085-t16138*t28819/216+t16323*t28732
     #/54+t16326*t28732/108-t16125*t28732/216+t16125*t28717/108+t16125*t
     #28829/216-t10712*t28606/432+t10712*t28653/216+t10712*t28660/432+t1
     #6138*t28698/108-t16138*t28709/108+t16138*t28722/216
      wI164pc1 = t22989+t28869+t28889+t28890+t28893+t28894+t28896+t28923
      wI165pc1 = wI139pc1
      wI166pc1 = wI140pc1
      t28926 = alf**2
      t28927 = t28926*alf
      t28928 = 1/PI
      t28929 = t28927*t28928
      A1e0c1p1t1w = 1024.D0/3.D0*t28929*wI163pc1
      A1e0c2p1t1w = -1152*t28929*wI164pc1
      A1e0c3p1t1w = 512.D0/3.D0*t28929*wI165pc1
      A1e0c4p1t1w = 2560.D0/3.D0*t28929*wI166pc1
      A1e0c5p1t1w = 128.D0/27.D0*t28927*(72*wI147pc1+90*wI148pc1+180*wI1
     #14pc1+243*wV35+176*wI150pc1-162*wI151pc1+432*wI118pc1-81*wI152pc1-
     #135*wI153pc1+252*wI154pc1)*t28928
      A1e0c6p1t1w = 32*t28927*(8*wI155pc1+8*wI156pc1+15*wI122pc1-16*wV36
     #+16*wI158pc1-24*wI159pc1+12*wI160pc1+22*wI161pc1)*t28928
      t28959 = 72*wI71pc1
      t28960 = 90*wI73pc1
      t28961 = 108*wV33
      t28962 = 128*wI74pc1
      t28963 = 27*wI75pc1
      t28964 = 168*wI133pc1
      t28965 = 54*wI77pc1
      t28966 = 170*wI134pc1
      A1e0c7p1t1w = -64.D0/27.D0*t28927*(-t28959-t28960+t28961-t28962-t2
     #8963-t28964-t28965-t28966)*t28928
      t28971 = 72*wI79pc1
      t28972 = 90*wI81pc1
      t28973 = 27*wV34
      t28974 = 128*wI82pc1
      t28975 = 240*wI83pc1
      A1e0c8p1t1w = -64.D0/27.D0*t28927*(-t28971-t28972+t28973-t28974+t2
     #8975-170*wI162pc1)*t28928
      t28993 = 18*wI91pc1+132*wI92pc1+54*wI141pc1+132*wI94pc1-3384*wI142
     #pc1-176*wI96pc1+44*wI97pc1-44*wI98pc1-54*wI143pc1+27*wI144pc1+18*w
     #I145pc1-22*wI146pc1
      A1e0c9p1t1w = -64.D0/9.D0*t28927*t28993*t28928
      A1e0c10p1t1w = -64*t28927*(6*wI101pc1+6*wI102pc1-2*wI103pc1+2*wI10
     #4pc1-8*wI105pc1-wI106pc1)*t28928
      A1e0c11p1t1w = 128.D0/9.D0*t28927*(6*wI107pc1+6*wI108pc1+2*wI109pc
     #1+2*wI110pc1+8*wI111pc1+wI112pc1)*t28928
      t29021 = t28927*(6*wI49pc1+6*wI50pc1-8*wI51pc1+2*wI52pc1-2*wI53pc1
     #-wI54pc1)*t28928
      A1e0c12p1t1w = 256.D0/9.D0*t29021
      A1e0c13p1t1w = 1024.D0/3.D0*t28929*wI137pc1
      A1e0c14p1t1w = -1152*t28929*wI138pc1
      A1e0c15p1t1w = 512.D0/3.D0*t28929*wI139pc1
      A1e0c16p1t1w = 2560.D0/3.D0*t28929*wI140pc1
      A1e0c17p1t1w = 32.D0/27.D0*t28927*(288*wI113pc1+720*wI114pc1+243*w
     #V31+704*wI116pc1+432*wI117pc1+1728*wI118pc1+864*wI119pc1+1008*wI12
     #0pc1)*t28928
      A1e0c18p1t1w = 8.D0/3.D0*t28927*(96*wI121pc1+180*wI122pc1-96*wV32+
     #192*wI124pc1+144*wI125pc1-288*wI126pc1+216*wI127pc1+265*wI128pc1)*
     #t28928
      A1e0c19p1t1w = 64.D0/27.D0*t28927*(72*wI129pc1+90*wI130pc1-t28961+
     #128*wI132pc1+t28963+t28964+t28965+t28966)*t28928
      A1e0c20p1t1w = -64.D0/27.D0*t28927*(-t28971-t28972+t28973-t28974+t
     #28975-170*wI136pc1)*t28928
      A1e0c21p1t1w = 128.D0/9.D0*t28927*(-9*wI91pc1-66*wI92pc1-27*wI93pc
     #1-66*wI94pc1+36*wI95pc1+88*wI96pc1-22*wI97pc1+22*wI98pc1-27*wI99pc
     #1+11*wI100pc1)*t28928
      A1e0c22p1t1w = -2560.D0/3.D0*t28929*wI85pc1
      A1e0c23p1t1w = 3072*t28929*wI86pc1
      A1e0c24p1t1w = 4096.D0/3.D0*t28929*wI87pc1
      A1e0c25p1t1w = 13312.D0/3.D0*t28929*wI88pc1
      A1e0c26p1t1w = 16.D0/27.D0*t28927*(1440*wI55pc1-864*wV25+2250*wI57
     #pc1+2560*wI58pc1+810*wI59pc1-7500*wI60pc1+2700*wI61pc1+3825*wI62pc
     #1)*t28928
      A1e0c27p1t1w = 256.D0/27.D0*t28927*(72*wI63pc1-54*wV26+90*wI65pc1+
     #128*wI66pc1+27*wI67pc1+189*wI68pc1+81*wI69pc1+173*wI70pc1)*t28928
      A1e0c28p1t1w = 256.D0/27.D0*t28927*(t28959-27*wV27+t28960+t28962+t
     #28963+186*wI76pc1+t28965+173*wI78pc1)*t28928
      A1e0c29p1t1w = 64.D0/27.D0*t28927*(288*wI79pc1-27*wV28+360*wI81pc1
     #+512*wI82pc1-1248*wI83pc1+692*wI84pc1)*t28928
      A1e0c30p1t1w = -160.D0/9.D0*t28927*(120*wI29pc1+120*wI30pc1-160*wI
     #31pc1+9*wI32pc1+40*wI33pc1-40*wI34pc1-27*wI35pc1-20*wI36pc1)*t2892
     #8
      A1e0c31p1t1w = 512.D0/3.D0*t28927*(6*wI37pc1+6*wI38pc1+2*wI39pc1+2
     #*wI40pc1+8*wI41pc1+wI42pc1)*t28928
      A1e0c32p1t1w = 1024.D0/9.D0*t28927*(6*wI43pc1+6*wI44pc1+2*wI45pc1+
     #2*wI46pc1+8*wI47pc1+wI48pc1)*t28928
      A1e0c33p1t1w = -2048.D0/9.D0*t29021
      A1e0c1p1t2w = A1e0c1p1t1w
      A1e0c2p1t2w = A1e0c2p1t1w
      A1e0c3p1t2w = A1e0c3p1t1w
      A1e0c4p1t2w = A1e0c4p1t1w
      A1e0c5p1t2w = A1e0c5p1t1w
      A1e0c6p1t2w = A1e0c6p1t1w
      A1e0c7p1t2w = A1e0c7p1t1w
      A1e0c8p1t2w = A1e0c8p1t1w
      A1e0c9p1t2w = A1e0c9p1t1w
      A1e0c10p1t2w = A1e0c10p1t1w
      A1e0c11p1t2w = A1e0c11p1t1w
      A1e0c12p1t2w = A1e0c12p1t1w
      A1e0c13p1t2w = A1e0c13p1t1w
      A1e0c14p1t2w = A1e0c14p1t1w
      A1e0c15p1t2w = A1e0c15p1t1w
      A1e0c16p1t2w = A1e0c16p1t1w
      A1e0c17p1t2w = A1e0c17p1t1w
      A1e0c18p1t2w = A1e0c18p1t1w
      A1e0c19p1t2w = A1e0c19p1t1w
      A1e0c20p1t2w = A1e0c20p1t1w
      A1e0c21p1t2w = A1e0c21p1t1w
      A1e0c22p1t2w = A1e0c22p1t1w
      A1e0c23p1t2w = A1e0c23p1t1w
      A1e0c24p1t2w = A1e0c24p1t1w
      A1e0c25p1t2w = A1e0c25p1t1w
      A1e0c26p1t2w = A1e0c26p1t1w
      A1e0c27p1t2w = A1e0c27p1t1w
      A1e0c28p1t2w = A1e0c28p1t1w
      A1e0c29p1t2w = A1e0c29p1t1w
      A1e0c30p1t2w = A1e0c30p1t1w
      A1e0c31p1t2w = A1e0c31p1t1w
      A1e0c32p1t2w = A1e0c32p1t1w
      A1e0c33p1t2w = A1e0c33p1t1w
      A1e0c1p1t3w = A1e0c1p1t2w
      A1e0c2p1t3w = A1e0c2p1t2w
      A1e0c3p1t3w = A1e0c3p1t2w
      A1e0c4p1t3w = A1e0c4p1t2w
      A1e0c5p1t3w = A1e0c5p1t2w
      A1e0c6p1t3w = A1e0c6p1t2w
      A1e0c7p1t3w = A1e0c7p1t2w
      A1e0c8p1t3w = A1e0c8p1t2w
      A1e0c9p1t3w = A1e0c9p1t2w
      A1e0c10p1t3w = A1e0c10p1t2w
      A1e0c11p1t3w = A1e0c11p1t2w
      A1e0c12p1t3w = A1e0c12p1t2w
      A1e0c13p1t3w = A1e0c13p1t2w
      A1e0c14p1t3w = A1e0c14p1t2w
      A1e0c15p1t3w = A1e0c15p1t2w
      A1e0c16p1t3w = A1e0c16p1t2w
      A1e0c17p1t3w = A1e0c17p1t2w
      A1e0c18p1t3w = A1e0c18p1t2w
      A1e0c19p1t3w = A1e0c19p1t2w
      A1e0c20p1t3w = A1e0c20p1t2w
      A1e0c21p1t3w = A1e0c21p1t2w
      A1e0c22p1t3w = A1e0c22p1t2w
      A1e0c23p1t3w = A1e0c23p1t2w
      A1e0c24p1t3w = A1e0c24p1t2w
      A1e0c25p1t3w = A1e0c25p1t2w
      A1e0c26p1t3w = A1e0c26p1t2w
      A1e0c27p1t3w = A1e0c27p1t2w
      A1e0c28p1t3w = A1e0c28p1t2w
      A1e0c29p1t3w = A1e0c29p1t2w
      A1e0c30p1t3w = A1e0c30p1t2w
      A1e0c31p1t3w = A1e0c31p1t2w
      A1e0c32p1t3w = A1e0c32p1t2w
      A1e0c33p1t3w = A1e0c33p1t2w

      RETURN
      END

