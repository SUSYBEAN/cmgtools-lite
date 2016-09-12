
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Auxiliary functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function s13q(z,l1,l2,l3,l4)
      double precision z,l1,l2,l3,l4
      double precision u,r,t,s34,s34p,s34m,en,k,pt,
     & ctheta,stheta,cphi,beta,gamma,s34s,
     & term1,term2,term3

      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s34m = (r-t)*(1-r*t)/r
      s34p = (r+t)*(1+r*t)/r
      s34 = l1*s34m
      s34s = l1*l3*(1-l3)*(1+u)**2/r/(r+t)/(1+r*t) 
      en = 0.5d0*(1-s34+z)
      k = (1-u)/(1+u)*en
      pt = dsqrt(dabs(en**2-k**2-z))
      
      ctheta = 2.0d0*l2-1.0d0
      stheta = dsqrt(dabs(1-ctheta**2))
      cphi = 2.0d0*l4-1.0d0
      
      beta = dsqrt(dabs(pt**2+k**2))/(1-en)
      gamma = (1-en)/dsqrt(s34)

      term1 = (1-en)/4.0d0*(1-beta*ctheta)
      term2 = (1-en)/4.0d0*(beta-ctheta)*en/(1+u)*(1-2.0d0*l3)/
     &     (1-l3+l3*z)*2.0d0/dsqrt(dabs(1+s34s*s34-2*s34s-2*z*s34s))
      term3 = stheta*cphi*dsqrt(dabs(s34))/4.0d0*r*dsqrt(l3*(1-l3)
     &     *(1-l1)*(1-l1*(r-t)*(1-r*t)/(r+t)/(1+r*t)))*2.0d0/
     &     dsqrt(dabs(1+s34s*s34-2*s34s-2*z*s34s))

      s13q = -2.0d0*(term1+term2+term3)      
      
      RETURN
      END

      double precision function s23q(z,l1,l2,l3,l4)
      double precision z,l1,l2,l3,l4
      double precision u,r,t,s34,s34p,s34m,en,k,pt,
     & ctheta,stheta,cphi,beta,gamma,s34s,
     & term1,term2,term3

CCC   Just flip the sign of term2,term3 in s13

      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s34m = (r-t)*(1-r*t)/r
      s34p = (r+t)*(1+r*t)/r
      s34 = l1*s34m
      s34s = l1*l3*(1-l3)*(1+u)**2/r/(r+t)/(1+r*t) 
      en = 0.5d0*(1-s34+z)
      k = (1-u)/(1+u)*en
      pt = dsqrt(dabs(en**2-k**2-z))
      
      ctheta = 2.0d0*l2-1.0d0
      stheta = dsqrt(dabs(1-ctheta**2))
      cphi = 2.0d0*l4-1.0d0
      
      beta = dsqrt(dabs(pt**2+k**2))/(1-en)
      gamma = (1-en)/dsqrt(s34)

      term1 = (1-en)/4.0d0*(1-beta*ctheta)
      term2 = (1-en)/4.0d0*(beta-ctheta)*en/(1+u)*(1-2.0d0*l3)/
     &     (1-l3+l3*z)*2.0d0/dsqrt(dabs(1+s34s*s34-2*s34s-2*z*s34s))
      term3 = stheta*cphi*dsqrt(dabs(s34))/4.0d0*r*dsqrt(l3*(1-l3)
     &     *(1-l1)*(1-l1*(r-t)*(1-r*t)/(r+t)/(1+r*t)))*2.0d0/
     &     dsqrt(dabs(1+s34s*s34-2*s34s-2*z*s34s))

      s23q = -2.0d0*(term1-term2-term3)      
      
      RETURN
      END

CCC   Spin correlation dot products:
CCC   Copy these over from the FEHiP photon phase-space; note that 
CCC   the l5,l6,l7 variables are already generated in the denom 
CCC   funtions; l6 is generated as sin^2(Pi/2*var(8))
CCC   Also need to return the c3l1 angle, for 1/s24^2 terms

      subroutine fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,
     & l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,
     & qflag)
      IMPLICIT double precision(a-z)
      double precision z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,
     & l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision B1,B2,Ome,FFF,Discr,s2anorm,s2bnorm,s1anorm,
     & s1bnorm,qflag
      double precision pi,pt3,en3,pz3,pt4,en4,pz4,pt1,en1,pz1,
     & pt2,en2,pz2,pT,en,pz,phi3l1,lbuf,lbuf2,u,r,t,c3v,cl1v,cphi
     
      pi = 3.1415926536d0

      cphi = 2.0d0*l8-1.0d0

      B1=z-s1v
      B2=z-s2v
      Ome = dabs(B1*B2-z)
      FFF = dsqrt(z*Ome*l5*dabs(1.0d0-l5))
      Discr = FFF*(2.0d0*l6-1d0)
      s2anorm = l5
      s2bnorm = (1.0d0-l5)
      s1anorm = l5*Ome + z*(1.0d0-l5)+2.0d0*Discr
      s1bnorm = (1.0d0-l5)*Ome + z*l5-2.0d0*Discr

      s1l1 = -s1anorm/B2
      s2l1 = -s2anorm*B2
      s1l2 = -s1bnorm/B2
      s2l2 = -s2bnorm*B2

      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)

      pt3 = dsqrt(dabs(s13*s23))
      en3 = -0.5d0*(s13+s23)
      pz3 = -0.5d0*(s23-s13)

      pt4 = dsqrt(dabs(s14*s24))
      en4 = -0.5d0*(s14+s24)
      pz4 = -0.5d0*(s24-s14)

      pt1 = dsqrt(dabs(s1l1*s2l1))
      en1 = -0.5d0*(s1l1+s2l1)
      pz1 = -0.5d0*(s2l1-s1l1)

      pt2 = dsqrt(dabs(s1l2*s2l2))
      en2 = -0.5d0*(s1l2+s2l2)
      pz2 = -0.5d0*(s2l2-s1l2)

      pT  = dsqrt(dabs(((s1v-z)*(s2v-z)-z)))
      en = -0.5d0*((s1v-z)+(s2v-z))
      pz = -0.5d0*((s2v-z)-(s1v-z))

CCC   Very little sensitivity to ptbuf or lbuf
      ptbuf = 0.0001d0
      lbuf = 0.0001d0
      lbuf2 = 0.9999d0

      if (qflag.lt.(0.0d0)) then
         if ((pt3.lt.ptbuf).or.(pt1.lt.ptbuf)) then
            s3l1 = 2.0d0*(en3*en1-pz3*pz1)
            c3l1 = 0.0d0
         else
            if (pT.lt.ptbuf) then
                c3l1=2.0d0*l6-1.0d0
            else
               c3v = (en3*en-pz3*pz-(s3v-z)/2.0d0)/pt3/pT
               if (c3v.lt.(-1d0)) then
                  c3v=-1d0
               elseif (c3v.gt.(1d0)) then
                  c3v=1d0
               endif
               
               cl1v = (en1*en-pz1*pz-z/2.0d0)/pt1/pT
               if (cl1v.lt.(-1d0)) then
                  cl1v=-1d0
               elseif (cl1v.gt.(1d0)) then
                  cl1v=1d0
               endif
               if (l7.lt.(0.5d0)) then
                  c3l1=c3v*cl1v+cphi*dsqrt(dabs(1d0-c3v**2)*
     &                 (1d0-cl1v**2))
               else
                  c3l1=c3v*cl1v-cphi*dsqrt(dabs(1d0-c3v**2)*
     &                 (1d0-cl1v**2))
               endif
            endif
            s3l1 = 2.0d0*(en3*en1-pz3*pz1-pt3*pt1*c3l1)
         endif
      else

CCC   If l3=0 or l3=1, then sqrt(s13*s23/l2) vanishes, and we can ignore c3l1
CCC   Therefore, we don't need to worry about the l3->0,1 limit of c3v

      if ((l2.lt.lbuf).or.(l1.gt.(lbuf2))) then
         c3v = 1.0d0-2.0d0*l4
      elseif (l1.lt.lbuf) then
         c3v = -1.0d0
      elseif (l3.lt.lbuf) then
         t1 = dsqrt(2.D0)
         t2 = l1*l2
         t8 = dsqrt(1-l2)
         t11 = dsqrt(1-l1)
         t13 = dsqrt(l2)
         t14 = dsqrt(l1)
         t16 = t**2
         t18 = dsqrt(1+t16)
         t20 = t13*t14*t18*t1
         t23 = t8*t11
         t26 = t14*t18*t1
         t48 = dsqrt(l1*t16-3*t2*t16+2*l2*t16-2*t26*t13*t11*t8*t+4
     &        *t20*t23*t*l4-t2+l1)
         c3v = t1*(2*t2*t-2*t*l2-2*l4*t8*t11*t20+t23*t13*t26)
     &        /t11/t13/t48/2.0d0
      elseif (l3.gt.(lbuf2)) then
         t1 = dsqrt(l2)
         t4 = dsqrt(1-l1)
         t7 = dsqrt(2.D0)
         t8 = t**2
         t10 = dsqrt(1+t8)
         t11 = t7*t10
         t12 = dsqrt(l1)
         t14 = dsqrt(1-l2)
         t15 = t12*t14
         t21 = 2**(1.D0/4.D0)
         t22 = t21**2
         t27 = t7*l2
         t49 = dsqrt(t7*t8*l1+2*t27*t8-3*t27*t8*l1-4*t*t14*t4*t1*t12
     &        *t10+8*l4*t*t14*t4*t1*t12*t10-t27*l1+t7*l1)
         c3v = -(2*t1*t*t4+2*t11*t15*l4-t11*t15)*t22*t21/t49/2   
      else

         c3v = (en3*en-pz3*pz-(s3v-z)/2.0d0)/pt3/pT
      endif

      if ((pt1.lt.ptbuf).or.(pT.lt.ptbuf)) then
         cl1v = -1.0d0+2.0d0*l6
      else
         cl1v = (en1*en-pz1*pz-z/2.0d0)/pt1/pT
      endif

      if (c3v.lt.(-1.0d0)) then
         c3v=-1.0d0
      elseif (c3v.gt.(1.0d0)) then
         c3v=1.0d0
      endif
     
      if (cl1v.lt.(-1.0d0)) then
         cl1v=-1.0d0
      elseif (cl1v.gt.(1.0d0)) then
         cl1v=1.0d0
      endif

      if (l7.gt.(0.5d0)) then
         c3l1=c3v*cl1v+cphi*dsqrt(dabs(1d0-c3v**2)*(1d0-cl1v**2))
      else
         c3l1=c3v*cl1v-cphi*dsqrt(dabs(1d0-c3v**2)*(1d0-cl1v**2))
      endif
      
      if (c3l1.lt.(-1.0d0)) then
         c3l1=-1.0d0
      elseif (c3l1.gt.(1.0d0)) then
         c3l1=1.0d0
      endif

      s3l1 = 2.0d0*(en3*en1-pz3*pz1-pt3*pt1*c3l1)
      endif

      s1l2 = s1v-z-s1l1
      s2l2 = s2v-z-s2l1 
      s3l2 = -z-s3l1+s3v 
      s4l1 = -z-s1l1-s2l1-s3l1
      s4l2 = s4v+s1l1+s2l1+s3l1

      RETURN
      END

CCC   This subroutine returns s3l1,s4l1,s3l2,s4l2 with 1-z extracted

      subroutine zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      IMPLICIT double precision(a-z)
      double precision z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,s4l2p,c3l1,
     & qflag
      double precision B1,B2,Ome,FFF,Discr,s2anorm,s2bnorm,s1anorm,
     & s1bnorm
      double precision pi,pt3,en3,pz3,pt4,en4,pz4,pt1,en1,pz1,
     & pt2,en2,pz2,pT,en,pz
      double precision s13,s23,s14,s24,s34,s1v,s2v,pt3p,en3p,pz3p,
     & pt4p,en4p,pz4p,s1l1p,s2l1p,zbuf,ptbuf,lbuf,lbuf2,u,r,t,cphi
     
      cphi = 2.0d0*l8-1.0d0
 
      pi = 3.1415926536d0

      B1=z-s1vp*(1.0d0-z)
      B2=z-s2vp*(1.0d0-z)
      Ome = dabs(B1*B2-z)
      FFF = dsqrt(z*Ome*l5*dabs(1d0-l5))
      Discr = FFF*(2d0*l6-1d0)
      s2anorm = l5
      s2bnorm = (1d0-l5)
      s1anorm = l5*Ome + z*(1d0-l5)+2d0*Discr
      s1bnorm = (1d0-l5)*Ome + z*l5-2d0*Discr

      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)

      s1l1 = -s1anorm/B2
      s2l1 = -s2anorm*B2
      s1l2 = -s1bnorm/B2
      s2l2 = -s2bnorm*B2

      s13=(1.0d0-z)*s13p
      s23=(1.0d0-z)*s23p
      s14=(1.0d0-z)*s14p
      s24=(1.0d0-z)*s24p
      s34=(1.0d0-z)*s34p
      s1v=(1.0d0-z)*s1vp
      s2v=(1.0d0-z)*s2vp

      pt3 = dsqrt(dabs(s13*s23))
      en3 = -0.5d0*(s13+s23)
      pz3 = -0.5d0*(s23-s13)

      pt4 = dsqrt(dabs(s14*s24))
      en4 = -0.5d0*(s14+s24)
      pz4 = -0.5d0*(s24-s14)

CCC   These are defined with 1-z extracted

      pt3p = dsqrt(dabs(s13p*s23p))
      en3p = -0.5d0*(s13p+s23p)
      pz3p = -0.5d0*(s23p-s13p)

      pt4p = dsqrt(dabs(s14p*s24p))
      en4p = -0.5d0*(s14p+s24p)
      pz4p = -0.5d0*(s24p-s14p)

      pt1 = dsqrt(dabs(s1l1*s2l1))
      en1 = -0.5d0*(s1l1+s2l1)
      pz1 = -0.5d0*(s2l1-s1l1)

      pt2 = dsqrt(dabs(s1l2*s2l2))
      en2 = -0.5d0*(s1l2+s2l2)
      pz2 = -0.5d0*(s2l2-s1l2)

      pT  = dsqrt(dabs(((s1v-z)*(s2v-z)-z)))
      en = -0.5d0*((s1v-z)+(s2v-z))
      pz = -0.5d0*((s2v-z)-(s1v-z))

CCC   Very little sensitivity to ptbuf or lbuf
      ptbuf = 0.0001d0
      lbuf = 0.0001d0
      lbuf2 = 0.9999d0

      if (qflag.lt.(0.0d0)) then
         if ((pt3.lt.ptbuf).or.(pt1.lt.ptbuf)) then
            s3l1p = 2.0d0*(en3p*en1-pz3p*pz1)
            c3l1 = 0.0d0
         else
            if (pT.lt.ptbuf) then
                c3l1=2.0d0*l6-1.0d0
            else
               c3v = (en3*en-pz3*pz-(s3v-z)/2.0d0)/pt3/pT
               if (c3v.lt.(-1d0)) then
                  c3v=-1d0
               elseif (c3v.gt.(1d0)) then
                  c3v=1d0
               endif
               
               cl1v = (en1*en-pz1*pz-z/2.0d0)/pt1/pT
               if (cl1v.lt.(-1d0)) then
                  cl1v=-1d0
               elseif (cl1v.gt.(1d0)) then
                  cl1v=1d0
               endif
               if (l7.lt.(0.5d0)) then
                  c3l1=c3v*cl1v+cphi*dsqrt(dabs(1d0-c3v**2)*
     &                 (1d0-cl1v**2))
               else
                  c3l1=c3v*cl1v-cphi*dsqrt(dabs(1d0-c3v**2)*
     &                 (1d0-cl1v**2))
               endif
            endif
            s3l1p = 2.0d0*(en3p*en1-pz3p*pz1-pt3p*pt1*c3l1)
         endif
      else

CCC   If l3=0 or l3=1, then sqrt(s13*s23/l2) vanishes, and we can ignore c3l1
CCC   Therefore, we don't need to worry about the l3->0,1 limit of c3v

      if ((l2.lt.lbuf).or.(l1.gt.(lbuf2))) then
         c3v = 1.0d0-2.0d0*l4
      elseif (l1.lt.lbuf) then
         c3v = -1.0d0
      elseif (l3.lt.lbuf) then
         t1 = dsqrt(2.D0)
         t2 = l1*l2
         t8 = dsqrt(1-l2)
         t11 = dsqrt(1-l1)
         t13 = dsqrt(l2)
         t14 = dsqrt(l1)
         t16 = t**2
         t18 = dsqrt(1+t16)
         t20 = t13*t14*t18*t1
         t23 = t8*t11
         t26 = t14*t18*t1
         t48 = dsqrt(l1*t16-3*t2*t16+2*l2*t16-2*t26*t13*t11*t8*t+4
     &        *t20*t23*t*l4-t2+l1)
         c3v = t1*(2*t2*t-2*t*l2-2*l4*t8*t11*t20+t23*t13*t26)
     &        /t11/t13/t48/2.0d0
      elseif (l3.gt.(lbuf2)) then
         t1 = dsqrt(l2)
         t4 = dsqrt(1-l1)
         t7 = dsqrt(2.D0)
         t8 = t**2
         t10 = dsqrt(1+t8)
         t11 = t7*t10
         t12 = dsqrt(l1)
         t14 = dsqrt(1-l2)
         t15 = t12*t14
         t21 = 2**(1.D0/4.D0)
         t22 = t21**2
         t27 = t7*l2
         t49 = dsqrt(t7*t8*l1+2*t27*t8-3*t27*t8*l1-4*t*t14*t4*t1*t12
     &        *t10+8*l4*t*t14*t4*t1*t12*t10-t27*l1+t7*l1)
         c3v = -(2*t1*t*t4+2*t11*t15*l4-t11*t15)*t22*t21/t49/2   
      else

         c3v = (en3*en-pz3*pz-(s3v-z)/2.0d0)/pt3/pT
      endif

      if ((pt1.lt.ptbuf).or.(pT.lt.ptbuf)) then
         cl1v = -1.0d0+2.0d0*l6
      else
         cl1v = (en1*en-pz1*pz-z/2.0d0)/pt1/pT
      endif

      if (c3v.lt.(-1.0d0)) then
         c3v=-1.0d0
      elseif (c3v.gt.(1.0d0)) then
         c3v=1.0d0
      endif
     
      if (cl1v.lt.(-1.0d0)) then
         cl1v=-1.0d0
      elseif (cl1v.gt.(1.0d0)) then
         cl1v=1.0d0
      endif

      if (l7.gt.(0.5d0)) then
         c3l1=c3v*cl1v+cphi*dsqrt(dabs(1d0-c3v**2)*(1d0-cl1v**2))
      else
         c3l1=c3v*cl1v-cphi*dsqrt(dabs(1d0-c3v**2)*(1d0-cl1v**2))
      endif
      
      if (c3l1.lt.(-1.0d0)) then
         c3l1=-1.0d0
      elseif (c3l1.gt.(1.0d0)) then
         c3l1=1.0d0
      endif

      s3l1p = 2.0d0*(en3p*en1-pz3p*pz1-pt3p*pt1*c3l1)
      endif

CCC   What I want to do is the following: expand
CCC   s2l1=-l5+(1-z)*s2l1p
CCC   s1l1=-(1-l5)+(1-z)*s1l1p
CCC   To get s1l1p, what I will do is take (s1l1(z)-s1l1(1))/(1-z)
CCC   I will do this numerically; for z less than some buffer
CCC   value, go to the derivative

      s2l1p = l5*(1.0d0+s2vp)
      zbuf = 0.99999d0
      if (z.lt.zbuf) then
         s1l1p = (s1l1+(1.0d0-l5))/(1.0d0-z)
      else
         s1l1p = -(1.0d0-l3)*(1.0d0-l5)-2.0d0*(2.0d0*l6-1)*dsqrt(l3
     & *(1.0d0-l3)*(1.0d0-l1)*l5*(1.0d0-l5))
      endif

      s1l2 = s1v-z-s1l1
      s2l2 = s2v-z-s2l1
      s3l2p = -s3l1p-s13p-s23p-s34p
      s4l1p = 1.0d0-s1l1p-s2l1p-s3l1p
      s4l2p = s13p+s23p+s1l1p+s2l1p+s3l1p

      RETURN
      END

CCC   This routine gets the regulatory factors for the 1/s34^2 topology

      double precision function s34spinc(z,l1,l2,l3,l4,l5,l6,l7,l8)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,l5,l6,l7,l8,l1buf,cphi,u,r,t

      if (l7.gt.(0.5d0)) then
         cphi = 2.0d0*l8-1.0d0
      else
         cphi = -(2.0d0*l8-1.0d0)
      endif

      l1buf = 0.999d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)

      if (l1.lt.l1buf) then

      t1 = r**2
      t2 = t1*l1
      t3 = t2*t
      t4 = l1*r
      t7 = -1+l1
      t8 = r*t
      t10 = t7*(-1+t8)
      t12 = r-t
      t14 = t**2
      t15 = r*t14
      t16 = l1*t
      t17 = t4*t14
      t18 = r+t1*t+t+t15-t4+t3+t16-t17
      t21 = l5-1
      t25 = dsqrt(dabs(l5*l6*t10*t12*t18*(-1+l6)*t21))
      t31 = l2-1
      t36 = dsqrt(dabs(-l2*l4*t7*t18*(l4-1)*t31/r))
      t42 = 1/t18/(r+t)
      t43 = 1/t7
      t48 = t1**2
      t49 = t48*r
      t50 = t49*l2
      t51 = l1**2
      t52 = t51*l1
      t53 = dsqrt(l1)
      t54 = t53*t52
      t55 = t14**2
      t56 = t55*t14
      t57 = t54*t56
      t64 = dsqrt(dabs(-t10*t12*t18*t21*l5))
      t65 = t49*t64
      t66 = t65*l4
      t71 = dsqrt(dabs(l2*t7*t31*r*t18))
      t72 = t71*t
      t73 = t72*l1
      t76 = t1*r
      t77 = t76*l2
      t78 = t53*t51
      t79 = t78*t14
      t80 = t79*l5
      t83 = t49*l5
      t84 = t51**2
      t85 = t53*t84
      t88 = t48*l1
      t89 = t72*l5
      t92 = t48**2
      t93 = t92*t1
      t95 = t53*t84*l1
      t97 = t55*t
      t98 = l5*t97
      t100 = t53*l1
      t102 = l5*t14
      t105 = t78*t55
      t108 = t48*t76
      t109 = t108*t64
      t111 = t14*t
      t112 = t71*t111
      t113 = t112*t51
      t116 = t92*r
      t117 = t116*l2
      t118 = t56*l5
      t119 = t118*t53
      t121 = l2*t14
      t122 = t121*t95
      t125 = t49*t84
      t127 = l4*l5
      t128 = t127*t55
      t131 = t92*l5
      t132 = l2*t111
      t133 = t132*t85
      t136 = -26*t50*t57-8*t66*t73-27*t77*t80-5*t83*t85-8*t88*t89-t93*t9
     #5*t98-2*t76*t100*t102-2*t50*t105+4*t109*l4*t113+t117*t119-5*t83*t1
     #22-24*t125*t71*t128-40*t131*t133
      t137 = t48*t1
      t138 = t137*t54
      t141 = t76*l6
      t142 = t141*t64
      t143 = l4*t71
      t145 = t143*t*t51
      t148 = t93*l2
      t149 = t98*t53
      t152 = t78*t56
      t155 = t108*l2
      t156 = t100*t97
      t157 = t156*t64
      t160 = t48*l6
      t161 = t160*t64
      t163 = t143*t55*t51
      t166 = t93*t97
      t167 = l2*t53
      t170 = r*l2
      t173 = r*t64
      t176 = t49*t51
      t178 = t143*l5
      t181 = t108*l1
      t185 = t92*t111
      t186 = t78*l5
      t190 = l2*t100
      t193 = t116*t14
      t194 = t167*l5
      t197 = -2*t138*t132+16*t142*t145+4*t148*t149+8*t83*t152-6*t155*t15
     #7-32*t161*t163-2*t166*t167+2*t170*t80-2*t173*t113+12*t176*t56*t178
     #+4*t181*t56*t178-6*t185*t186+2*t93*t111*t190-3*t193*t194
      t199 = t76*l5
      t202 = t137*l2
      t203 = t54*t97
      t204 = t203*l5
      t207 = t55*t111
      t208 = t207*t100
      t210 = t76*t78
      t211 = l5*t55
      t214 = t49*t55
      t215 = t100*l5
      t217 = t85*t55
      t220 = t48*l2
      t221 = t78*t97
      t224 = t137*t111
      t230 = t108*t71
      t231 = t127*t14
      t234 = t49*t111
      t237 = t54*t64*l6
      t242 = t100*t56
      t245 = 12*t199*t122+54*t202*t204+t131*t208+12*t210*t211-t214*t215+
     #6*t155*t217-6*t220*t221+2*t224*t167-84*t224*t52*t178+4*t230*t231+2
     #4*t234*l2*t237+4*t166*t190+4*t117*t242
      t246 = t116*t78
      t247 = t56*l2
      t248 = t247*l5
      t251 = l2*t
      t252 = t251*t85
      t257 = t207*l5
      t258 = t257*t53
      t262 = t1*l6*t64
      t263 = t71*l1
      t264 = t263*t14
      t267 = t108*t100
      t268 = t55**2
      t270 = t268*l2*l5
      t273 = t77*t78
      t275 = t111*t64*l6
      t278 = t1*l5
      t279 = t132*t95
      t282 = t48*l5
      t283 = t251*t95
      t287 = t100*t*t64
      t290 = t132*t78
      t293 = t1*l2
      t295 = t54*t111*l5
      t300 = t54*l5
      t301 = t300*t55
      t304 = 11*t246*t248-4*t131*t252+t93*t54*t98-6*t202*t258+4*t262*t26
     #4+3*t267*t270-12*t273*t275-14*t278*t279-2*t282*t283+2*t50*t287+10*
     #t65*t290-18*t293*t295-2*t77*t287-102*t155*t301
      t307 = t49*t85
      t310 = t76*t64
      t311 = t72*t51
      t314 = t137*t64
      t315 = t314*l4
      t316 = t71*t14
      t317 = t316*t51
      t320 = t55*t71
      t321 = t320*l5
      t324 = t71*t51
      t325 = t324*t97
      t328 = t137*l5
      t329 = t328*l4
      t330 = t112*t84
      t333 = t108*t53
      t334 = t251*t64
      t337 = t48*t14
      t341 = t92*t76
      t342 = t341*l2
      t343 = t85*l5
      t344 = t343*t55
      t346 = t78*t111
      t347 = t346*l5
      t350 = t137*t
      t351 = t52*t71
      t352 = t351*l5
      t355 = t76*t51
      t361 = -42*t307*t248+4*t310*t311-8*t315*t317-2*t176*t321+4*t65*t32
     #5+48*t329*t330-2*t333*t334-16*t337*l2*t237+t342*t344+6*t148*t347+2
     #4*t350*t352+16*t355*t71*t231+8*t185*t352
      t362 = t137*t78
      t364 = l2*t207*l5
      t367 = t132*t53
      t370 = t49*l6
      t371 = t370*t64
      t375 = t78*t*t64
      t380 = t167*t64
      t383 = t108*t111
      t387 = l2*t95
      t388 = t387*t97
      t391 = t93*l5
      t392 = l2*t85
      t393 = t392*t97
      t396 = t49*t52
      t397 = t396*t71
      t400 = t76*t52
      t401 = t71*l5
      t402 = t401*t14
      t405 = t100*t111
      t406 = t405*l5
      t413 = 4*t362*t364+4*t142*t367-8*t371*t325+2*t50*t375-58*t220*t295
     #-4*t337*t380+12*t383*l2*t237-10*t278*t388-2*t391*t393+60*t397*t231
     #+2*t400*t402+6*t293*t406-4*t293*t347+4*t77*t375
      t416 = t54*t*l5
      t424 = t108*t52
      t425 = t424*t71
      t428 = t108*t54
      t431 = t137*t14
      t435 = t49*l1
      t441 = t137*l6*t64
      t442 = l2*t55
      t443 = t442*t54
      t446 = t71*t97
      t447 = t446*l1
      t450 = t108*t14
      t451 = l5*l2
      t452 = t451*t85
      t455 = t100*t55
      t456 = t455*t64
      t459 = t48*t111
      t462 = -6*t293*t416+5*t328*t283-4*t396*t56*t178+48*t425*t231-29*t4
     #28*t248-12*t431*l2*t237+16*t435*t178-8*t66*t113-16*t441*t443+4*t37
     #1*t447+50*t450*t452+6*t220*t456-8*t459*t178
      t463 = t108*t85
      t466 = t242*l5
      t472 = t108*t84
      t475 = t48*t78
      t478 = t92*l2
      t479 = t215*t
      t487 = t97*t85
      t491 = t76*t53
      t493 = t48*t
      t494 = t95*l5
      t497 = t54*t14
      t500 = 24*t463*t248-15*t155*t466+4*t155*t287+t391*t388+16*t472*t32
     #1+44*t475*t364+4*t478*t479-8*t214*l1*t178-7*t170*t344+26*t131*t487
     #+t333*t118-t491*t102-3*t493*t494-2*t77*t497
      t506 = t48*t85
      t509 = t78*t207
      t514 = t310*l4
      t519 = t76*t111
      t523 = t48*t64
      t524 = t523*l4
      t525 = t320*t51
      t530 = t92*t55
      t533 = l2*t64*l6
      t536 = t72*t52
      t539 = t48*t52
      t540 = t446*l5
      t545 = -12*t199*t57-4*t506*t257-12*t220*t509+2*t199*t242-8*t514*t3
     #11+4*t214*t178+12*t519*l2*t237+16*t524*t525-2*t310*t73-8*t530*t53*
     #t533+4*t514*t536-6*t539*t540+5*t148*t204
      t549 = t56*t71*l5
      t554 = t55*t95
      t557 = t53*l5
      t560 = t1*t111
      t570 = t111*t85
      t573 = t53*t56
      t581 = -84*t328*t133+2*t472*t549-8*t523*t317+30*t83*t554-4*t459*t5
     #57+5*t560*t300-2*t65*t112+2*t83*t71*t84-2*t65*t72-2*t278*t570-t77*
     #t573+4*t142*t536+2*t396*t549-2*t214*t401
      t586 = t92*t84
      t590 = t76*t
      t594 = t49*t78
      t601 = t324*l5
      t608 = t95*t207*l5
      t617 = t316*t52
      t620 = -4*t2*t111*t178+8*t586*t97*t178+4*t590*l2*t237-22*t594*t270
     #-6*t148*t406-4*t523*t316-24*t350*t601+24*t88*t111*t178+3*t478*t608
     #+52*t425*t128-35*t594*t248-24*t176*t178-12*t161*t617
      t621 = t221*l5
      t624 = t85*t268
      t626 = t54*t207
      t633 = t156*l5
      t638 = t1*t64
      t639 = t263*t55
      t650 = t316*t84
      t653 = t48*t71
      t654 = l5*t
      t657 = t49*t53
      t659 = -10*t293*t621+t50*t624+12*t282*t626+2*t293*t221+4*t142*t73+
     #31*t202*t633-20*t282*t388-4*t638*t639-30*t328*t252+12*t88*t97*t178
     #+4*t310*t113-32*t83*l4*t650+4*t653*t654-t657*l5
      t664 = t112*t52
      t690 = t49*t56
      t693 = -10*t594*l5+12*t66*t664+12*t176*t401-2*t310*t112+28*t278*t3
     #93-14*t220*t479-8*t371*t145+48*t350*t51*t178-12*t594*t97*t533-28*t
     #176*t71*t231-28*t224*t601-4*t586*t540-2*t690*t401
      t696 = t108*l5
      t704 = t49*t14
      t707 = t76*t97
      t708 = t71*t64
      t711 = t215*t14
      t716 = t137*t85
      t719 = t143*t14
      t722 = t108*t78
      t723 = t268*l5
      t733 = -4*t15*t452-24*t696*l4*t650+5*t199*t152-8*t371*t113-2*t704*
     #t557-2*t707*t708+5*t77*t711-8*t142*t311+14*t716*t257-16*t161*t719-
     #4*t722*t723-2*t478*t156-19*t117*t301-4*t333*t97*t533
      t735 = t49*t54
      t738 = t48*t100
      t744 = t203*t64
      t747 = t116*t54
      t751 = t143*t111*t52
      t754 = t300*t14
      t760 = t54*t55
      t765 = t137*l1
      t773 = 33*t735*t211+11*t738*t98+t83*t57-4*t262*t617-2*t155*t744-3*
     #t747*t211-24*t371*t751-59*t50*t754-4*t83*t143*t84+6*t117*t760-2*t6
     #38*t264-16*t765*t97*t178+4*t371*t367-8*t142*t447
      t774 = t116*l5
      t780 = t108*t55
      t783 = r*l5
      t785 = t92*t85
      t788 = t143*t111
      t795 = t442*t95
      t804 = t92*t54
      t807 = 5*t774*t554+24*t282*t570+t170*t217+2*t780*t557+t783*t217-6*
     #t785*t364-8*t371*t788+4*t478*t149+8*t314*t443+5*t783*t795-8*t66*t3
     #25-12*t396*t321+12*t220*t626-11*t804*t98
      t814 = t351*t55
      t827 = t79*t64
      t833 = t116*t84
      t837 = t638*l4
      t842 = 2*t774*t152-4*t441*t814+4*t310*t447+12*t696*t650-78*t202*t3
     #47-114*t220*t204+14*t176*t402+4*t478*t827-55*t138*t98+t77*t79-4*t8
     #33*t71*t128+8*t837*t639-4*t315*t814
      t844 = t143*t111*t51
      t854 = t155*t100
      t856 = t97*t64*l6
      t865 = t708*t111
      t868 = t143*t
      t880 = 16*t142*t844-30*t131*t393+3*t116*t95*t118+25*t199*t795+12*t
     #854*t856+14*t478*t221-2*t310*t536+76*t155*t344+4*t141*t865+4*t65*t
     #868+3*t783*t122-5*t220*t608-2*t396*t446*t64+6*t435*t549
      t884 = t92*l1
      t897 = l2*t54
      t898 = t111*l5
      t908 = 4*t478*t626+4*t884*t97*t178-2*t884*t540-3*t117*t105+3*t246*
     #t211-t116*t55*t167+t193*t167+t897*t898+2*t117*t57-3*t155*t79+t307*
     #t121+t117*t79+16*t202*t626
      t911 = t137*t97
      t914 = t48*t54
      t927 = t48*t97
      t939 = -10*t77*t760-16*t911*t392+21*t914*t98+t155*t343+6*t523*t814
     #-4*t463*t723+t77*t215-26*t424*t321+6*t914*t132-2*t927*t557+35*t50*
     #t80+36*t224*t186+2*t109*t367-8*t653*t127*t
      t956 = t49*t95
      t965 = t211*t53
      t970 = -6*t459*t392-17*t210*t248-71*t224*t300+14*t131*t388-8*t355*
     #t402-18*t594*t102+4*t293*t203+10*t956*t118-10*t493*t300-3*t155*t15
     #2-8*t50*t344+2*t117*t965-12*t161*t639
      t971 = t108*t51
      t978 = t14*t85
      t999 = t54*t268
      t1002 = t263*l5
      t1008 = -12*t971*t56*t178-28*t765*t111*t178-30*t696*t978-12*t524*t
     #814-8*t371*t73-32*t459*t51*t178+9*t155*t711-4*t472*t56*t178-13*t73
     #8*t364-8*t148*t633-26*t220*t633-4*t50*t999-2*t780*t1002-8*t185*t51
     #*t178
      t1010 = t112*l5
      t1024 = t14*t64*l6
      t1030 = t143*t51*t97
      t1037 = t*t85
      t1042 = t137*t55
      t1047 = 2*t2*t1010-2*t310*t367+14*t50*t217+16*t131*l4*t330+14*t765
     #*t1010+4*t293*t100*t1024+38*t478*t204+16*t371*t1030+4*t371*t311+76
     #*t328*t570+10*t282*t1037+16*t83*t650+4*t1042*t380+30*t696*t554
      t1062 = t405*t64
      t1068 = t143*t111*l1
      t1079 = 15*t117*t344-2*t181*t549+14*t307*t247+33*t722*t248-4*t155*
     #t999-4*t155*t215+2*t638*t617+2*t155*t1062+88*t282*t393-8*t142*t106
     #8+12*t478*t416-4*t155*t300+16*t371*t844+4*t315*t639
      t1083 = t396*t97
      t1090 = t92*t52
      t1102 = t400*t71
      t1105 = t708*l6
      t1112 = t72*t84
      t1118 = -8*t1083*t143*t64*l6-2*t50*t744-12*t1090*t97*t178-77*t77*t
     #344-64*t696*t217-t170*t754-28*t971*t55*t178-12*t1102*t128+4*t707*t
     #1105+8*t765*t89+12*t514*t664-4*t282*t1112-4*t202*t100*t1024
      t1120 = t116*t85
      t1130 = t100*t14*t64
      t1135 = t*t64*l6
      t1140 = t137*t51
      t1149 = t95*t14
      t1154 = t459*t300+4*t1120*t248-t560*t494-16*t185*t52*t178+4*t690*t
     #178-4*t478*t1130+4*t657*l2*t1135-8*t273*t1135+24*t1140*t97*t178-12
     #*t1140*t540-t780*t215+32*t696*t57+10*t696*t1149+16*t329*t1112
      t1161 = r*l6*t64
      t1175 = t108*t95
      t1177 = t112*l1
      t1185 = 6*t193*t452+4*t214*t1002-8*t1161*t844+15*t328*t388-6*t50*t
     #157+12*t539*t97*t178-2*t230*t102-10*t696*t122-t1175*t270-2*t310*t1
     #177+t478*t406-12*t50*t711+18*t735*t270
      t1189 = t76*t71
      t1192 = t220*t100
      t1194 = t55*t64*l6
      t1200 = t76*l1
      t1203 = t76*t14
      t1209 = t49*t71
      t1214 = t78*t268
      t1217 = t49*t100
      t1225 = t956*t270+6*t523*t443+2*t1189*t102-12*t1192*t1194+12*t355*
     #t55*t178+2*t1200*t402-32*t1203*t452+4*t173*l4*t113-4*t1209*t127-4*
     #t199*t978+6*t50*t1214-4*t1217*t723-13*t50*t215-10*t48*t95*t98
      t1232 = t100*t268
      t1240 = t143*t97*l1
      t1249 = t897*t64
      t1252 = t116*t51
      t1256 = t387*t56
      t1259 = t186*t
      t1262 = -2*t220*t346+4*t459*t401-4*t50*t1232-2*t560*t352-8*t262*t5
     #25-8*t371*t1240+2*t1209*t102-8*t142*t788+4*t142*t325-6*t519*t1249+
     #4*t1252*t55*t178+5*t83*t1256-44*t202*t1259
      t1263 = t76*t54
      t1274 = t1*t54
      t1277 = t105*l5
      t1295 = -5*t1263*t102-22*t478*t621-12*t88*t1010-3*t267*t102-2*t314
     #*t639+2*t1274*t132+6*t50*t1277+30*t202*t406+2*t199*t650+16*t142*t1
     #240-28*t328*t393+16*t161*t317-2*t48*t53*t654+8*t337*t1249
      t1299 = t202*t78
      t1324 = -3*t50*t79+8*t1299*t1024+16*t262*t163+58*t202*t416+122*t20
     #2*t295+t1120*t247-t83*t387+2*t657*t118-2*t155*t455+4*t514*t1177-7*
     #t185*t300+8*t362*t257+9*t463*t247+3*t77*t105
      t1333 = t50*t100
      t1336 = t116*t64
      t1343 = t92*t100
      t1355 = 6*t155*t1214-4*t77*t186-2*t185*t190-2*t224*t190+12*t1333*t
     #856-2*t1336*t290+t657*t723+t333*t451+12*t125*t321+4*t1343*t98+2*t1
     #209*l5-4*t1189*t231+4*t435*t71*t231-6*t185*t194
      t1370 = t143*l1*t14
      t1388 = t48*t51
      t1391 = 4*t370*t865-6*t131*t626+3*t1343*t898+16*t524*t317-16*t441*
     #t1370+2*t92*t64*t443+6*t478*t346+4*t220*t208-2*t202*t156-8*t854*t1
     #135-8*t435*t401-2*t109*t113+12*t1388*t540
      t1394 = t455*l5
      t1398 = t282*l4
      t1420 = -3*t50*t119+t77*t1394+2*t1336*t367+16*t1398*t330+4*t220*t1
     #130-2*t435*t402+49*t428*t211-19*t1263*t211+7*t50*t343+t1175*t723-8
     #*t396*t401+4*t50*t1062+10*t199*t1256-2*t590*t1249
      t1422 = t1*t14
      t1431 = t137*t100
      t1442 = t77*t100
      t1445 = t137*t95
      t1452 = -4*t1422*l2*t237+5*t1217*t102-t774*t57+22*t50*t186-5*t1431
     #*t257+6*t155*t242-16*t722*t211-2*t214*t190-8*t514*t447-4*t1442*t27
     #5-5*t1445*t257+6*t523*t264+4*t65*t788
      t1466 = t1*t97
      t1481 = 2*t716*t364-8*t220*t827+42*t224*t352+5*t199*t554+64*t328*t
     #487+t450*t167+20*t131*t570-2*t1466*t190+9*t210*t102-18*t459*t186-4
     #*t854*t275+4*t780*l1*t178+t199*t387-2*t65*t311
      t1487 = t108*l6*t64
      t1512 = -36*t463*t118-4*t1487*t367-8*t282*t330+24*t1299*t1194+45*t
     #478*t295+2*t833*t321-16*t927*t392-2*t491*t211+10*t735*l5+8*t1398*t
     #1112+2*t314*t814+6*t220*t149+4*t1161*t113
      t1543 = 2*t804*t364+16*t161*t525-8*t371*t868+8*t523*t719+54*t220*t
     #347-4*t657*t97*t533-3*t774*t1256-15*t293*t204-6*t310*t664+16*t88*t
     #*t178+53*t1263*t248+2*t282*t252-24*t328*t330+4*t514*t325
      t1549 = t76*t55
      t1552 = t49*t
      t1579 = -6*t383*t1249+2*t202*t1130-2*t1549*t190+4*t1552*l2*t237+8*
     #t459*t194+4*t142*t1177-12*t77*t57+4*t314*t317-2*t65*t447+4*t638*t5
     #25+4*t1442*t1135+6*t155*t186+4*t176*t55*t178-8*t142*t143*t*t52
      t1607 = 14*t971*t321+24*t397*t128+t307*t723+6*t493*t194-8*t514*t11
     #3-t170*t105-r*t78*t211-8*t478*t78*t1024+4*t66*t1177+4*t66*t311+4*t
     #1083*t1105-8*t523*t525-6*t65*t664+6*t400*t321
      t1637 = 16*t396*t178-6*t478*t456-25*t696*t795-t333*t270-12*t117*t7
     #54+24*t161*t1370-4*t314*t264+39*t77*t754-12*t138*t364+50*t155*t127
     #7-12*t524*t264-24*t424*t402-12*t722*t97*t533
      t1652 = t268*t53
      t1659 = t143*t52*t55
      t1671 = -12*t234*t1249+8*t493*t352+2*t927*t190+4*t780*t178-16*t493
     #*t52*t178+4*t50*t54*t856+t50*t1652-4*t391*t133-8*t837*t525+24*t161
     #*t1659-12*t435*t56*t178+27*t155*t80-6*t202*t149+3*t657*t451
      t1673 = t76*t85
      t1699 = 6*t1673*t247+3*t696*t242-14*t328*t626+14*t202*t203+9*t77*t
     #217+6*t1673*t118+4*t1274*t98+t155*t624-4*t77*t343+6*t594*t723-4*t7
     #35*t723+6*t202*t221-24*t475*t98
      t1702 = t92*t97
      t1709 = t155*t78
      t1730 = t143*t16
      t1733 = -6*t478*t203-4*t1702*t392-5*t83*t242+4*t185*t601+4*t1709*t
     #1135-2*t657*t334-6*t424*t549-8*t1042*t53*t533-42*t1673*t248-4*t110
     #2*t231+56*t224*t51*t178+9*t202*t479+4*t310*t788-8*t142*t1730
      t1739 = t97*l2*t64
      t1766 = 13*t1217*t270+6*t722*t1739+2*t155*t497+5*t170*t1277+3*t463
     #*t270-4*t116*l6*t64*t367+7*t155*t1394+5*t738*t654-12*t478*t1259-4*
     #t1120*t118+2*t1422*t1249+2*t333*t1739+5*t117*t1277
      t1790 = t220*t78
      t1797 = -4*t1466*t392+6*t117*t711-8*t1192*t1024-4*t278*t487-16*t30
     #7*t118+2*t185*t167+8*t765*t540+4*t117*t78*t275+6*t431*t1249-8*t441
     #*t317-4*t155*t1232+24*t1790*t1194+2*t109*t290-4*t224*t401
      t1807 = t143*t64
      t1816 = t137*t52
      t1828 = t143*l1*t55
      t1831 = t442*t78
      t1834 = 2*t1549*t401-4*t523*t320+66*t220*t621-8*t1487*t844+4*t1083
     #*t1807+4*t66*t447-4*t1333*t1135-12*t524*t617+16*t1816*t540+8*t337*
     #t53*t533+4*t707*t1807+8*t160*t708*t14+24*t161*t1828-12*t314*t1831
      t1844 = t137*t84
      t1864 = -2*t50*t301+6*t523*t639-4*t202*t827+t293*t479+3*t77*t119+1
     #6*t1844*t97*t178-18*t50*t300-56*t83*t217-4*t431*t380+3*t117*t80+16
     #*t371*t1730-8*t1333*t275+4*t530*t380-5*t450*t194
      t1872 = t143*t14*t52
      t1890 = t143*t55
      t1900 = -5*t774*t795+8*t262*t1872+2*t350*t194-8*t142*t1030-16*t765
     #*t*t178-3*t170*t301+2*t342*t965-4*t1209*t231+8*t441*t1659-16*t161*
     #t1890-13*t1431*t898+4*t155*t54*t856-12*t523*t1831
      t1932 = 2*t202*t608+20*t328*t279-32*t1816*t97*t178+4*t333*l2*t1135
     #+6*t77*t300-54*t914*t364-6*t804*t132+8*t431*t53*t533-5*t747*t248-1
     #5*t282*t279-13*t1431*t98-12*t696*t1256-2*t428*t270+4*t837*t264
      t1963 = t956*l5-5*t1431*t654+4*t199*t217+2*t1*t78*t98+6*t428*t723+
     #7*t293*t633-10*t478*t347-4*t92*t78*t98-12*t524*t639-71*t155*t754-8
     #*t181*t71*t231-4*t199*l4*t650-8*t262*t1370
      t1992 = t92*t95
      t1995 = 6*t50*t466+6*t310*t290-4*t478*t208-4*t1200*t71*t231+4*t514
     #*t73+6*t202*t346+4*t50*t242+12*t738*t898+5*t1217*l5-32*t472*t71*t1
     #28-2*t155*t375-8*t1844*t540+8*t224*t178-15*t1992*t98
      t2025 = t143*t14*t51
      t2028 = 24*t161*t1872+8*t262*t639+2*t77*t1062-66*t202*t621+4*t371*
     #t1177+36*t704*t452-20*t971*t71*t231-6*t355*t321-10*t185*t494+8*t88
     #4*t111*t178+16*t1790*t1024+6*t971*t549-32*t161*t2025
      t2050 = t53*t207
      t2059 = -2*t293*t1130-24*t1388*t97*t178+15*t83*t1149+10*t1431*t364
     #+10*t131*t279-2*t117*t455+9*t50*t152+2*t704*t190-15*t459*t494+8*t7
     #85*t257+2*t478*t2050-3*t1992*t257+2*t459*t190+22*t220*t203
      t2089 = -10*t155*t105-6*t594*t211-2*t193*t190+3*t1203*t194-2*t65*t
     #1177-2*t478*t258-6*t88*t540-2*t1702*t167-2*t65*t367+4*t65*t113-8*t
     #117*t466-7*t77*t466+5*t704*t194+10*t971*t402
      t2090 = t85*t207
      t2095 = t50*t78
      t2124 = -2*t478*t2090-4*t884*t1010-20*t2095*t275+4*t1549*l1*t178-2
     #*t1252*t321-8*t328*t1112-3*t657*t270+69*t735*t248-t155*t119+4*t370
     #*t708*t+78*t77*t301+12*t371*t664-8*t707*t71*l4*t64*l6+8*t278*t252
      t2148 = t1*t100
      t2155 = -3*t278*t283-7*t307*t270+8*t911*t178-2*t722*t270-10*t224*t
     #194+4*t77*t965-8*t441*t1828-8*t371*t1068-2*t1552*t1249-10*t50*t760
     #-2*t2148*t98-2*t2148*t898-2*t560*t392
      t2183 = 2*t137*t207*t557-2*t148*t346+29*t428*t102+t333*t102-2*t148
     #*t203-t451*t111*t95+t342*t105-t341*t55*t167-12*t155*t57+20*t350*t1
     #86-30*t350*t300+2*t202*t2050+4*t911*t557+2*t224*t557
      t2210 = t267*t723-44*t83*t978-15*t696*t152+30*t362*t98-30*t1445*t9
     #8-30*t224*t494+6*t155*t760+15*t1175*t118+2*t77*t242-7*t722*t102-4*
     #t224*t392-12*t282*t509-4*t911*t401
      t2233 = t478*t100
      t2240 = -31*t220*t406+2*t117*t1394+8*t523*t1890-4*t1709*t275-2*t31
     #0*t325-8*t131*t330-31*t77*t1277-4*t2095*t1135+44*t735*t102-3*t342*
     #t1277+30*t278*t133+8*t2233*t1024+6*t523*t617-2*t911*t167
      t2266 = 5*t77*t152-6*t202*t2090+t155*t1652-t491*t118+4*t282*t208-5
     #*t350*t494-12*t202*t509-2*t927*t167+t1673*t121+4*t282*t487+3*t199*
     #t1149+2*t657*t1739+16*t459*t601
      t2297 = 4*t65*t73-30*t396*t402-48*t350*t52*t178+8*t441*t264-4*t220
     #*t2090+20*t328*t1037+3*t1343*t364-12*t161*t814-7*t117*t152-12*t161
     #*t443+4*t560*t52*t178-6*t155*t965-4*t837*t617+4*t148*t295
      t2328 = 12*t424*t56*t178+8*t315*t264-t50*t573+16*t441*t2025+6*t109
     #0*t540-8*t142*t113-16*t262*t1828+8*t160*t708*t55-24*t142*t751+6*t5
     #0*t1394+8*t220*t1259+6*t594*t1739+12*t2233*t1194+4*t441*t639
      t2359 = -12*t161*t264-2*t50*t965-4*t478*t633-2*t780*t401-5*t774*t2
     #17-6*t176*t549+42*t282*t133-2*t1549*t1002+12*t142*t664+4*t181*t402
     #-4*t92*l6*t64*t443+4*t1487*t113+28*t506*t364-4*t1549*t178
      s34spinc = 8*(t3+r-t4+t)*t25*t*t36*t42*t43*cphi+(t1391+t545+t2089+
     #t1079+t1420+t1118+t1864+t1452+t1047+t1543+t304+t1481+t2297+t1008+t
     #1900+t1733+t1671+t2028+t970+t939+t500+t1579+t2210+t1932+t908+t1963
     #+t1154+t2124+t842+t2359+t462+t733+t1797+t807+t1995+t245+t1607+t163
     #7+t773+t1355+t1699+t2266+t2328+t1295+t2183+t2155+t693+t1324+t197+t
     #620+t361+t659+t880+t413+t1225+t1185+t136+t1262+t1834+t2059+t1512+t
     #1766+t581+t2240)/t1/(t3+t15-t17+r-t4+t16)*t42*t43/(t1+t8+l1-t4*t)
      else
      t1 = r**2
      t3 = t1*t+t
      t5 = r*t
      t8 = (1-t5)*(r-t)
      t10 = 1+t1
      t13 = 1-l5
      t17 = dsqrt(dabs(l5*l6*t8*t*t10*(1-l6)*t13))
      t24 = l2-1
      t29 = dsqrt(dabs(-l2*l4*t*t10*(1-l4)*t24/r))
      t31 = 1/t3/2
      t34 = 1/(r+t)
      t39 = t1*r
      t41 = dsqrt(2.D0)
      t42 = t39*l6*t41
      t47 = dsqrt(dabs(t8*t*t10*t13*l5))
      t53 = dsqrt(dabs(-2*l2*t24*t5*t10))
      t54 = t**2
      t55 = t54*t
      t57 = t47*l4*t53*t55
      t60 = t1**2
      t61 = t60*r
      t62 = t61*t41
      t64 = t47*t53*t55
      t69 = t60**2
      t70 = t69*r
      t71 = t70*l2
      t72 = t54**2
      t75 = t69*t1
      t76 = t72*t
      t92 = t60*t1
      t97 = t69*t39*t72
      t100 = l2*t55
      t103 = t69*l2
      t106 = l2*t72
      t109 = -24*t42*t57-6*t62*t64+12*t42*t64-10*t71*t72+2*t75*t76*l5-8*
     #t39*l5*t72-20*t72*t61*l2-2*t70*l5*t72+8*t55*t60*l2+12*t55*t92*l2-2
     #*t97*l2+2*t100*t1+8*t103*t76-2*t106*r
      t110 = t76*t60
      t118 = t39*t41
      t121 = l4*t53*t55
      t124 = t60*t39
      t125 = t124*l2
      t128 = r*l5
      t137 = l5*t1
      t140 = l2*t76
      t143 = t75*l2
      t146 = t61*l5
      t149 = l5*l2
      t152 = t92*l5
      t155 = t92*l2
      t158 = 8*t110*l2-10*t72*t39*l2+8*t110*l5+12*t118*t47*t121-20*t125*
     #t72-2*t128*t72+8*t69*t55*l2+8*t69*t76*l5+2*t137*t76+2*t140*t1+2*t1
     #43*t76-12*t146*t72+2*t149*t55+12*t152*t76+12*t155*t76
      t166 = t61*l6*t41
      t170 = t124*l6*t41
      t176 = r*l6*t41
      t185 = t76*l5
      t197 = t60*l5
      t200 = 2*t143*t55-8*t124*l5*t72-24*t166*t57+4*t170*t64+6*t100*t137
     #-8*t176*t57-6*t118*t64+4*t176*t64+2*t128*t106-4*t143*t185-2*t143*t
     #55*l5-24*t155*t185-4*t140*t137+2*t97*t149+4*t197*t100
      t208 = t72*l5
      t220 = r*t41
      t228 = t124*t41
      t239 = -16*t197*t140-4*t152*t100-6*t69*l5*t100+10*t71*t208+20*t146
     #*t106+10*l2*t39*t208-16*t103*t185+20*t125*t208+4*t220*t47*t121-8*t
     #170*t57+12*t166*t64-2*t228*t64+12*t62*t47*t121-2*t220*t64+4*t228*t
     #47*t121
      s34spinc = -16*t3*t17*t*t29*t31*t34*cphi-(t109+t158+t200+t239)/t1/
     #t3*t31*t34/t10

      endif

      RETURN
      END

      double precision function s34spinc2(z,l1,l2,l3,l4,l5,l6,l7,l8)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,l5,l6,l7,l8,l1buf,cphi,u,r,t
      
      if (l7.gt.(0.5d0)) then
         cphi = 2.0d0*l8-1.0d0
      else
         cphi = -(2.0d0*l8-1.0d0)
      endif

      l1buf = 0.999d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
     
      if (l1.lt.l1buf) then
     
      t2 = l1-1
      t4 = r**2
      t6 = t**2
      t7 = r*t6
      t8 = l1*r
      t9 = t4*l1
      t10 = t9*t
      t11 = l1*t
      t12 = t8*t6
      t13 = r+t4*t+t+t7-t8+t10+t11-t12
      t16 = l2-1
      t17 = 1/r
      t21 = dsqrt(dabs(-l2*l4*t2*t13*(l4-1)*t16*t17))
      t23 = r*t
      t25 = t2*(-1+t23)
      t27 = r-t
      t30 = l5-1
      t34 = dsqrt(dabs(l5*l6*t25*t27*t13*(l6-1)*t30))
      t36 = t8*t
      t43 = 1/t13/(1+t23)
      t44 = 1/t2
      t49 = t4**2
      t50 = t49*l5
      t51 = l1**2
      t52 = t51**2
      t53 = dsqrt(l1)
      t54 = t53*t52
      t55 = l2*t54
      t56 = t6*t
      t57 = t55*t56
      t60 = t4*r
      t61 = t60*l6
      t66 = dsqrt(dabs(-t25*t27*t13*t30*l5))
      t67 = t61*t66
      t72 = dsqrt(dabs(l2*t2*t16*r*t13))
      t73 = t72*t51
      t74 = t73*t
      t77 = t49*r
      t78 = t77*l6
      t79 = t78*t66
      t80 = l4*t72
      t81 = t51*t56
      t82 = t80*t81
      t85 = t77*l5
      t86 = t53*l1
      t87 = l2*t86
      t90 = t77*t86
      t91 = t6**2
      t92 = l2*t91
      t93 = t92*l5
      t96 = r*t72
      t97 = t96*l4
      t98 = t66*t56
      t99 = t98*l1
      t102 = t51*l1
      t103 = t53*t102
      t104 = t60*t103
      t107 = t49*t4
      t108 = t107*l2
      t109 = t103*t56
      t112 = t60*l5
      t113 = t80*t102
      t116 = t49**2
      t117 = t116*l6
      t118 = t66*t53
      t119 = t118*t6
      t122 = t49*t60
      t123 = t122*l6
      t124 = t53*t51
      t125 = t66*t124
      t126 = t125*t56
      t129 = t107*l5
      t131 = t53*t52*l1
      t132 = t131*t
      t135 = t112*l4
      t136 = t72*l1
      t137 = t136*t6
      t140 = l2*t131
      t141 = t140*t
      t144 = t107*t66
      t145 = l2*t124
      t146 = t145*t6
      t149 = 140*t50*t57-4*t67*t74+24*t79*t82+18*t85*t87+24*t90*t93+8*t9
     #7*t99+6*t104*t91+12*t108*t109-12*t112*t113+4*t117*t119+4*t123*t126
     #+4*t129*t132+12*t135*t137-8*t50*t141-12*t144*t146
      t150 = t4*l5
      t151 = t140*t56
      t154 = t91*t6
      t156 = t154*l2*l5
      t159 = t107*t72
      t160 = t159*l4
      t161 = t66*l6
      t162 = t51*t6
      t163 = t161*t162
      t166 = t145*t
      t169 = t49*l6
      t170 = t125*t6
      t173 = t122*t91
      t178 = t77*t6
      t179 = t53*l2
      t182 = t169*t66
      t183 = t72*t102
      t184 = t183*t91
      t187 = t122*l5
      t188 = t183*t6
      t191 = t116*r
      t192 = t191*t53
      t195 = t60*t72
      t196 = t195*l4
      t197 = t161*t11
      t200 = t122*l2
      t201 = t200*t66
      t202 = l6*t124
      t203 = t202*t
      t206 = t49*t124
      t207 = t91*t56
      t209 = l2*t207*l5
      t212 = t191*l5
      t213 = t179*t6
      t216 = -8*t150*t151+18*t104*t156-24*t160*t163+32*t50*t166-12*t169*
     #t170-2*t173*t72+44*t129*t166-2*t178*t179+8*t182*t184+14*t187*t188+
     #4*t192*t93-24*t196*t197-8*t201*t203+20*t206*t209+4*t212*t213
      t220 = t50*l4
      t221 = t73*t56
      t226 = t122*t86
      t232 = t122*t124
      t235 = t136*t91
      t238 = t87*t56
      t241 = t49*t72
      t242 = t241*l4
      t243 = t66*t102
      t244 = t243*t6
      t247 = t4*l6
      t248 = t247*t66
      t251 = t77*t72
      t252 = l4*t66
      t253 = t252*t
      t256 = t140*t91
      t259 = t77*t154
      t264 = 4*t61*t126-60*t220*t221+6*t90*t154+3*t226*t154-4*t187*t86+t
     #187*t53-3*t232*t154+12*t135*t235+32*t79*t238+8*t242*t244+8*t248*t1
     #84-8*t251*t253+32*t187*t256+4*t259*t179+2*t226*t156
      t267 = t4*t66
      t270 = t191*t124
      t273 = t66*t86
      t274 = t273*t56
      t279 = t77*t103
      t280 = t91**2
      t282 = t280*l2*l5
      t285 = t80*t
      t288 = t107*t53
      t289 = t91*t
      t290 = l2*t289
      t293 = t77*t280
      t297 = t124*t6
      t300 = t77*l2
      t301 = t103*t91
      t304 = t4*t72
      t305 = t304*l4
      t307 = t161*t51*t91
      t310 = t72*t52
      t311 = t310*t
      t314 = t179*t56
      t319 = 12*t150*t57+4*t267*t213+6*t270*t93-8*t61*t274+6*t187*t213+2
     #0*t279*t282+16*t79*t285+6*t288*t290-5*t293*l5*t86+43*t85*t297-20*t
     #300*t301+24*t305*t307-36*t220*t311+6*t50*t314-4*t226*t91
      t322 = t122*t103
      t324 = t77*t66
      t325 = t179*t
      t328 = t49*t53
      t329 = t290*l5
      t332 = t72*t6
      t335 = t122*t51
      t336 = t91*t72
      t341 = t60*t91
      t344 = t60*t102
      t347 = t60*t154
      t354 = t60*t66
      t355 = t56*t53
      t358 = t4*t289
      t365 = t322*t154+8*t324*t325-8*t328*t329+2*t112*t332-2*t335*t336+4
     #*t187*t336-4*t341*t80-2*t344*t332+6*t347*t73+2*t112*t136-6*t112*t7
     #3-2*t354*t355-4*t358*t80+2*t267*t297+4*t50*t285
      t366 = t80*t6
      t369 = t273*t
      t372 = t55*t289
      t375 = t183*t
      t380 = t122*t53
      t384 = t60*l2
      t385 = t384*t66
      t395 = t310*t6
      t400 = t60*t52
      t401 = t72*l5
      t402 = t401*t91
      t405 = -4*t112*t366-8*t123*t369+126*t50*t372-12*t150*t375-16*t300*
     #t274+2*t380*t6+t187*t54-8*t385*t203+4*t150*t221+t322*t6+10*t279*t9
     #1-10*t288*t329-42*t85*t395+6*t187*t124-30*t400*t402
      t407 = t49*t66
      t412 = t4*t124
      t417 = t4*t103
      t424 = t60*t124
      t427 = t55*t6
      t430 = t87*t
      t436 = t118*t
      t439 = r*t51
      t446 = -12*t407*t146-r*t86*t91+2*t412*t56+3*t412*t289-2*t417*t289+
     #r*t103*t91+2*t251*t154-3*t424*t154-124*t85*t427-16*t324*t430-8*t19
     #6*t161*t289+4*t61*t436+2*t439*t402-12*t324*l4*t221-t380*t154
      t450 = t77*t53
      t454 = r*t52
      t463 = r*l6*t66
      t464 = t136*t
      t467 = t60*t51
      t468 = t336*l4
      t471 = t80*t51
      t476 = t107*t289
      t479 = t116*l5
      t482 = t252*t6
      t488 = t107*t56
      t493 = 2*t85*t53+t328*t56+t450*t6-48*t220*t74-6*t454*t402-8*t80*t6
     #6*l6*t51*t6+4*t463*t464-20*t467*t468+12*t112*t471+2*t358*t72-8*t47
     #6*t113-8*t479*t151-12*t241*t482-6*t191*t54*t93-4*t488*t72-6*t232*t
     #93
      t497 = t66*t51
      t498 = t497*t
      t503 = l6*t53
      t504 = t503*t
      t507 = t54*t154
      t510 = t150*l4
      t513 = t497*t6
      t516 = t49*t289
      t519 = l2*t103
      t520 = t519*t56
      t523 = t107*t124
      t526 = t77*t102
      t532 = t161*t56*l1
      t535 = t116*t53
      t539 = t519*t6
      t542 = -4*t97*t498-8*t354*t430-8*t201*t504+16*t187*t507-12*t510*t3
     #11+8*t305*t513+4*t516*t72-172*t129*t520+66*t523*t329+20*t526*t468-
     #2*t450*t154-16*t97*t532+3*t535*t56-t85*t131+66*t187*t539
      t545 = t60*t86
      t548 = t251*l4
      t551 = r*l5
      t552 = t551*l4
      t555 = t116*t66
      t558 = t77*t124
      t565 = t136*t154
      t571 = t116*t103
      t574 = t116*t124
      t580 = 16*t85*t124+3*t545*t154-48*t548*t532-16*t552*t184-4*t555*t1
     #46-3*t558*t6-10*t450*t156+6*t322*t91-16*t135*t565-t535*t289-3*t288
     #*t289-2*t571*t289+3*t574*t289-t574*t56-2*t571*t56
      t582 = t77*t51
      t585 = t140*t154
      t588 = t122*t66
      t589 = t86*t
      t598 = t73*t6
      t602 = t80*l5
      t605 = t77*t54
      t614 = t131*t154
      t617 = t73*t289
      t622 = -16*t582*t468+12*t112*t585+4*t588*t589+4*t467*t332-4*t304*t
     #81+6*t112*t183+24*t135*t598+24*t467*t154*t602-102*t605*t156+4*t526
     #*t366-16*t78*t369+6*t212*t146-6*t187*t614+12*t67*t617+4*t112*t213
      t625 = t86*t56
      t632 = t60*l1
      t633 = t80*t154
      t636 = t300*t66
      t638 = l6*t56*t53
      t641 = t80*t56
      t646 = t131*t91
      t648 = r*l2
      t649 = t91*t86
      t652 = t116*t4
      t653 = t652*l5
      t654 = t54*t56
      t656 = t191*l2
      t657 = t86*t6
      t662 = t54*t289
      t666 = t66*l1*t
      t671 = -6*t212*t427+4*t354*t625-10*t341*t183-204*t605*t93+12*t632*
     #t633-16*t636*t638+16*t79*t641+14*t479*t109-t551*t646+2*t648*t649-t
     #653*t654+4*t656*t657-6*t112*t614-8*t479*t662+4*t97*t666+4*t653*t31
     #4
      t676 = t519*t91
      t685 = t187*l4
      t692 = t85*l4
      t695 = t72*t154
      t698 = t131*t207
      t701 = t91*t53
      t704 = t49*t56
      t707 = t207*t103
      t712 = -12*t259*t471+76*t187*t676+4*t344*t366+4*t335*t468+8*t248*t
     #188+28*t685*t137+12*t169*t119+4*t344*t633-32*t692*t565-2*t112*t695
     #+4*t50*t698-2*t112*t701-2*t704*t179+20*t50*t707+9*t112*t657
      t713 = t54*t91
      t716 = t154*t86
      t719 = l2*l5
      t724 = t53*t6
      t731 = t72*t
      t736 = t207*t53
      t738 = t289*t124
      t741 = t49*l2
      t750 = 33*t112*t713+3*t112*t716-12*t605*t719-6*t112*t137-2*t267*t7
     #24+6*t407*t297+28*t85*t519-2*t129*t731+4*t150*t314+t129*t736-33*t1
     #29*t738+12*t741*t109-27*t50*t738+6*t328*t290+t293*l5*t53
      t754 = t252*t91
      t767 = t107*l6
      t768 = t767*t66
      t771 = t55*t
      t783 = t87*t6
      t788 = -66*t129*t662+4*t304*t754+46*t112*t188-4*t195*t253-36*t85*t
     #646-28*t129*t375-20*t558*t282+16*t768*t366+42*t129*t771+8*t324*t16
     #6-16*t516*t113+2*t77*t131*t282-92*t135*t188-18*t112*t783-8*t385*t6
     #38
      t789 = t183*t56
      t795 = t122*t52
      t806 = t72*t289
      t821 = t122*t6
      t824 = t122*t154
      t827 = 56*t510*t789-4*t439*t91*t602-12*t795*t402+2*t212*t256+4*t48
     #8*t183+8*t488*t80-2*t50*t731-6*t50*t806-2*t150*t355+8*t187*t87-16*
     #t112*t646+6*t187*t297-18*t85*t716-8*t85*t657-4*t821*t179+2*t824*t1
     #79
      t832 = t4*l2
      t833 = t103*t289
      t836 = t207*t124
      t843 = t53*t
      t845 = t103*t
      t849 = t49*t54
      t850 = l5*t207
      t858 = t243*t56
      t863 = t741*t66
      t864 = t503*t6
      t871 = 4*t832*t833-10*t50*t836-6*t488*t179-3*t50*t355+t479*t843+36
     #*t129*t845+t551*t301-15*t849*t850-2*r*t91*t183-4*t247*t170-4*t548*
     #t858+86*t85*t184-24*t863*t864+8*t212*t539+36*t50*t771
      t872 = t519*t
      t878 = t719*t
      t883 = t98*t53
      t894 = t125*t
      t901 = t116*t56
      t908 = -16*t479*t872+4*t195*t252*t289-4*t328*t878+8*t248*t598+4*t6
     #1*t883+2*t96*t498+8*t305*t244-4*t97*t858+4*t150*t520+4*t123*t894+2
     #*t476*t72-2*t704*t72-2*t901*t72+2*t195*t154+2*t341*t72
      t915 = t161*t102*t91
      t918 = t107*t54
      t925 = t124*t56
      t935 = t131*t289
      t940 = t77*t52
      t942 = l4*l5
      t943 = t942*t91
      t946 = t695*l5
      t951 = -2*t821*t72-12*t552*t188-16*t305*t915+30*t918*t209+2*t96*t8
     #58+2*t267*t332-21*t50*t925+62*t50*t833-2*t653*t355-2*t652*t56*t179
     #+24*t50*t935-88*t50*t789+84*t940*t72*t943-24*t582*t946+36*t50*t617
      t953 = t161*t102*t6
      t956 = t479*l4
      t961 = t73*t91
      t968 = t66*t72
      t969 = t968*t6
      t987 = l2*t6
      t992 = -16*t305*t953+16*t956*t221+6*t479*t771-12*t182*t961+24*t67*
     #t80*l1*t289-8*t767*t969-8*t107*t102*t641+12*t552*t395+20*t79*t464+
     #15*t112*t507-3*t212*t297-2*t656*t301+6*t212*t649+6*t558*t987-12*t1
     #69*t969
      t995 = t77*l1
      t1004 = t354*l4
      t1005 = t136*t289
      t1012 = t80*t289
      t1022 = t80*l1
      t1031 = -4*t995*t366-28*t685*t188+2*t8*t336+54*t112*t184-12*t1004*
     #t1005-18*t187*t783-2*t112*t310+8*t129*t1012-2*t526*t332+6*t259*t73
     #-16*t49*t102*t641-4*t112*t1022+8*t324*t314+18*t129*t74-2*t551*t137
      t1032 = t124*t
      t1041 = t140*t289
      t1046 = t122*l1
      t1049 = t108*t66
      t1050 = t202*t6
      t1058 = t124*t91
      t1063 = t519*t289
      t1070 = -2*t354*t1032+4*t354*t166-4*t8*t468-4*t195*t99-8*t479*t104
     #1+6*t656*t649-12*t1046*t468+24*t1049*t1050+48*t582*t154*t602-12*t4
     #67*t946-9*t112*t1058+60*t135*t395+4*t150*t1063+2*t251*t858-72*t220
     #*t617
      t1082 = t60*t53
      t1090 = t161*t51*t
      t1103 = t968*t
      t1106 = -8*t123*t274-2*t588*t1032-4*t324*t1032-2*t335*t332-4*t85*t
     #310+2*t1082*t92+24*t129*t617-16*t196*t532+8*t196*t1090-4*t85*t695+
     #2*t354*t731-t212*t646+2*t116*t72*t81-2*t380*t156-4*t61*t1103
      t1107 = t131*t56
      t1110 = t131*t6
      t1113 = t129*l4
      t1121 = t183*t289
      t1126 = t310*t56
      t1133 = t243*t91
      t1147 = 4*t150*t1107-6*t112*t1110-36*t1113*t74-60*t129*t789-2*t191
     #*t6*t179-4*t1004*t1121-12*t548*t498-76*t1113*t1126+4*t942*t221+6*t
     #150*t311+8*t305*t1133-8*t1046*t366+12*t768*t598+4*t588*t325+96*t58
     #2*t91*t602-2*t648*t301
      t1176 = t116*l2
      t1177 = t1176*t66
      t1184 = 10*t341*t145-6*t347*t87-48*t50*t1041-2*t178*t72-2*t85*t72+
     #4*t122*t102*t468-24*t226*t93-3*t232*t91-4*t942*t1126-54*t187*t427-
     #8*t182*t468+2*t450*t92-8*t1177*t864+8*t78*t436+4*t901*t80
      t1204 = t154*t103
      t1210 = t103*t6
      t1219 = -2*t632*t336-16*t242*t953+4*t476*t183-2*t344*t695+6*t407*t
     #332+8*t516*t183-16*t85*t614+51*t85*t507+10*t85*t280*t124-62*t85*t1
     #204+5*t85*t280*t54-4*t300*t1210-4*t79*t789+26*t112*t539-3*t479*t73
     #8
      t1233 = t122*t54
      t1249 = t72*t56
      t1250 = t1249*l1
      t1253 = -6*t479*t1032+4*t1176*t109+2*t535*t290+24*t129*t935+2*t117
     #6*t925-6*t112*t297-84*t1233*t93+10*t150*t372+3*t212*t713+10*t90*t2
     #82+t545*t6-4*t241*t1133-18*t467*t402+12*t50*t1012+24*t79*t1250
      t1277 = l6*t86*t
      t1284 = t154*t53
      t1290 = 6*t212*t657+36*t467*t91*t602-32*t1113*t221-8*t358*t113-2*t
     #200*t1210+4*t588*t314-12*t187*t395+6*t200*t657-2*t191*t91*t145+20*
     #t129*t707+16*t201*t1277+2*t551*t256-28*t479*t520+t187*t1284-13*t18
     #7*t1204-10*t251*t666
      t1295 = t497*t91
      t1306 = t145*t56
      t1326 = 12*t112*t649+6*t304*t1295-8*t479*t221-12*t248*t961-2*t384*
     #t1210-2*t270*t987+4*t588*t1306+176*t220*t789+2*t551*t598-3*t212*t3
     #01-4*t212*t1210-9*t112*t1204+16*t479*t372+57*t85*t1058-8*t4*t102*t
     #641
      t1331 = t4*t53
      t1334 = t124*t154
      t1337 = t54*t
      t1360 = -16*t187*t646+8*t187*t519+2*t1331*t290+6*t384*t1334-18*t50
     #*t1337+28*t50*t845-18*t108*t738+12*t741*t833+2*t85*t140-18*t112*t3
     #01+4*t129*t698-2*t112*t724+3*t653*t925-2*t150*t109+2*t347*t179
      t1388 = t310*t289
      t1393 = -8*t201*t638+102*t85*t713+4*t80*t513-16*t78*t274+3*t85*t70
     #1-32*t685*t184-16*t636*t504-2*t384*t657-4*t450*t719-2*t212*t724-8*
     #t467*t366+6*t324*t221+8*t78*t126+18*t50*t1388+3*t479*t925
      t1402 = t243*t
      t1409 = t60*t54
      t1429 = -2*t551*t676-24*t1113*t311-6*t1176*t738-8*t78*t1103+2*t195
     #*t1402+10*t129*t314-2*t85*t213-30*t1409*t156+4*t632*t366-16*t242*t
     #915+4*t4*t56*t183+70*t50*t109+2*t50*t843-63*t50*t662-5*t129*t355+4
     #*t832*t109
      t1435 = t832*t66
      t1464 = t107*t103
      t1467 = 8*t1435*t1050+6*t574*t329+32*t85*t585-t1082*t91-72*t129*t8
     #72+12*t112*t146+12*t454*t72*t943-12*t692*t137+16*t526*t946+4*t1176
     #*t833-6*t195*t666-3*t479*t1337-6*t479*t355+8*t479*t845-6*t1464*t28
     #9
      t1475 = t191*t86
      t1496 = 3*t523*t56+t104*t154-t652*t124*t56+t191*t103*t91-3*t1475*t
     #91-2*t1475*t6-2*t424*t6+t270*t6+t192*t6-5*t424*t91+2*t279*t6+4*t47
     #9*t1107+4*t61*t894+84*t692*t598-2*t328*t209
      t1521 = t140*t6
      t1531 = -44*t510*t1126-12*t384*t301-6*t901*t179+12*t150*t617-4*t26
     #7*t146+4*t354*t314+16*t123*t66*t238-2*t380*t92-6*t832*t738+20*t344
     #*t468+8*t248*t366+32*t85*t1521-8*t187*t468+4*t632*t468+8*t4*l4*t22
     #1
      t1562 = -4*t67*t1121-16*t196*t858+24*t50*t1107-42*t85*t598+4*t247*
     #t119+6*t85*t54+2*t279*t154+6*t212*t676-2*t354*t925-2*t632*t332+96*
     #t220*t1121+18*t424*t93+24*t548*t1090-6*t354*t617-2*t450*t282
      t1569 = t310*t154
      t1596 = -6*t653*t1306+12*t67*t464+8*t692*t1569+12*t407*t213-2*t212
     #*t701-2*t96*t666-4*t300*t1204-6*t545*t156-6*t288*t878+16*t85*t783+
     #86*t129*t109-114*t85*t301+26*t322*t156-2*t150*t833-8*t588*t430
      t1618 = t202*t56
      t1631 = -2*t288*t209+24*t129*t1107-6*t479*t1306-32*t526*t154*t602-
     #6*t450*t93+72*t85*t256+8*t67*t80*t102*t289-3*t187*t724-24*t1113*t1
     #388-8*t385*t1618+30*t479*t57-20*t85*t1022-8*t50*t641-38*t187*t301-
     #48*t1113*t617
      t1663 = 6*t380*t93+6*t354*t1005+3*t187*t1334+4*t335*t366-9*t558*t9
     #1-6*t1464*t56+16*t129*t221-114*t558*t93-3*t112*t1334-2*t384*t1204-
     #5*t150*t662-12*t79*t74-6*t150*t654-66*t1409*t93-3*t328*t289
      t1666 = t54*t6
      t1697 = 8*t97*t1090+62*t85*t1666+6*t150*t925-t85*t280*t131+12*t995
     #*t633+6*t85*t137+56*t1113*t375+18*t50*t311+2*t652*l2*t925+42*t187*
     #t713+4*r*t102*t468+t380*t91-8*t200*t274-t1082*t154+5*t206*t56-14*t
     #85*t103
      t1729 = 9*t206*t289-2*t535*t878+t104*t6-12*t692*t235+4*t135*t1569+
     #84*t692*t395+3*t187*t1058-86*t85*t146-40*t50*t375+80*t220*t375-12*
     #t129*t641-6*t824*t87-4*t832*t925+9*t187*t657-t187*t716
      t1747 = t4*t51
      t1762 = -6*t187*t1110+6*t173*t145-12*t767*t170-8*t385*t504+120*t11
     #13*t789+4*t479*t1126+8*t551*t184+60*t400*t72*t943+12*t1747*t1012-6
     #*t551*t395-2*t354*t843-48*t129*t151-2*t1233*t719+3*t212*t1666+6*t4
     #24*t156
      t1772 = l4*t154
      t1783 = t77*t91
      t1796 = 6*t150*t74-4*t150*t641+18*t112*t1666-21*t129*t1337-4*t251*
     #t1772+4*t85*t80+10*t85*t136-2*t200*t1204+8*t324*t625-10*t1783*t183
     #-33*t187*t1210+2*t995*t336-2*t173*t183-12*t232*t719+4*t588*t625
      t1807 = t49*t51
      t1830 = 4*t173*t80+6*t1046*t336-6*t995*t695+12*t108*t833-3*t212*t1
     #058+24*t1807*t1012+4*t424*t987+12*t187*t649+8*t200*t649-12*t200*t3
     #01-27*t129*t925+70*t129*t833-12*t259*t87-78*t129*t654+27*t187*t166
     #6-82*t85*t1210
      t1841 = t289*t53
      t1863 = -10*t741*t925-22*t129*t1032+6*t200*t1334+6*t479*t833+t479*
     #t1841+4*t479*t935-15*t479*t654+24*t50*t74+24*t510*t375+24*t182*t36
     #6+4*t161*t598-4*t463*t221-4*t995*t468-12*t212*t783-2*t150*t806
      t1888 = t968*t91
      t1895 = 156*t129*t57+2*t555*t297-108*t135*t184+20*t523*t209+12*t12
     #9*t1388+20*t548*t666-28*t150*t789-32*t129*t1121-12*t510*t1388+36*t
     #90*t156+8*t196*t99+3*t150*t738+4*t247*t1888+12*t479*t166+4*t552*t1
     #37
      t1908 = t161*t102*t56
      t1930 = -8*t97*t197-16*t85*t1110+5*t129*t1841+40*t85*t1334+6*t144*
     #t297+8*t97*t1908+16*t385*t1277+30*t50*t221+5*t85*t1284+12*t196*t66
     #6+4*t169*t1888-16*t344*t154*t602+132*t129*t372-2*t112*t1569-6*t112
     #*t235
      t1964 = 24*t863*t1050-10*t85*t280*t103+t85*t724+12*t300*t1334+54*t
     #206*t329-70*t50*t654+2*t653*t57+16*t67*t238+12*t49*l4*t221-8*t248*
     #t468+8*t187*t598-4*t304*t244+6*t551*t188+4*t241*t754+4*t354*t325+8
     #*t324*t1306
      t1976 = t49*t103
      t1996 = -12*t67*t1005-t1331*t289-2*t417*t56+3*t288*t56-6*t1976*t28
     #9-6*t1976*t56-6*t558*t154+9*t523*t289+t192*t91+t270*t91-3*t226*t6-
     #4*t187*t103-t450*t91+8*t112*t565+12*t112*t1521
      t2028 = -36*t220*t1388-4*t116*l4*t221-16*t150*t1121+78*t85*t188-32
     #*t1233*t156+36*t85*t471+30*t849*t209+4*t129*t285+4*t354*t1306+24*t
     #548*t99-8*t1435*t864+4*t150*t1012-4*t552*t598-8*t510*t221+32*t685*
     #t235
      t2060 = 64*t1113*t1121+228*t85*t676-4*t196*t498+18*t1783*t145-28*t
     #85*t113-4*t96*t99-10*t129*t836-12*t79*t221-12*t242*t1295-6*t108*t9
     #25+14*t85*t183+3*t129*t843+2*t582*t332-4*t117*t170-8*t201*t1618
      t2088 = t140*t207
      t2091 = -15*t918*t850-8*t241*t513-8*t187*t366+4*t50*t1841+t50*t736
     #+4*t50*t132-16*t50*t1032-172*t692*t184+124*t279*t156-12*t251*t99-1
     #6*t67*t789+36*t112*t676-6*t632*t695+6*t129*t1249-8*t50*t2088
      t2124 = -12*t85*t649-40*t548*t197+2*t195*t498-16*t636*t1618-8*t479
     #*t641-8*t384*t274+8*t1177*t1050+8*t195*t858-6*t232*t156+8*t67*t285
     #-4*t196*t1402-156*t692*t188+4*t123*t436+8*t582*t336+4*t178*t80
      t2155 = -2*t85*t336-6*t241*t81-6*t144*t724+4*t187*t332-2*t588*t843
     #+4*t479*t1249+24*t685*t395+8*t182*t188-4*t67*t375-8*t956*t1126+8*t
     #78*t883+12*t479*t314-24*t1049*t864+8*t552*t235+24*t242*t307
      t2188 = 12*t187*t585+4*t588*t166+4*t61*t968*t289-40*t1976*t209-42*
     #t940*t402+54*t129*t1306-8*t129*t141-12*t187*t146-2*t192*t92+8*t344
     #*t946+4*t85*t468-4*t304*t482-4*t247*t969-8*t159*t482+10*t467*t336
      t2202 = t107*t51
      t2223 = -2*t195*t66*t289-4*t476*t80-80*t558*t156-2*t968*t162-18*t7
     #41*t738-4*t195*t1772-6*t2202*t806-4*t324*t843-4*t324*t925+2*t995*t
     #332-4*t324*t355+8*t324*t589+8*t1004*t221+4*t555*t213+8*t67*t221-14
     #0*t50*t520
      t2239 = t80*t52
      t2258 = -4*t551*t235-12*t112*t598-10*t605*t282-4*t463*t789+12*t144
     #*t213-4*t582*t366+4*t112*t2239+32*t112*t256-12*t150*t1306+12*t1004
     #*t617-2*t267*t336-18*t85*t73-2*t588*t355-2*t588*t925+4*t821*t80
      t2290 = 4*t1046*t332+4*t358*t183-4*t304*t513-4*t85*t1569+4*t1082*t
     #93+16*t85*t565-8*t251*t252*t56+6*t150*t1388+4*t354*t589-12*t1807*t
     #806-4*t129*t806+4*t144*t332+4*t50*t1249+4*t324*t731+4*t324*t1249
      t2323 = -2*t526*t695-12*t1475*t93-16*t335*t91*t602-40*t1464*t209-2
     #*t407*t336+12*t2202*t1012+16*t242*t513+38*t129*t1126+6*t251*t498-8
     #*t129*t2088-4*t241*t244-2*t401*t81-14*t187*t137-6*t1747*t806-16*t1
     #87*t235
      t2358 = 8*t67*t1250+32*t636*t1277-3*t187*t701-48*t129*t1041+32*t19
     #6*t1908+32*t510*t1121-16*t305*t163+12*t767*t119-2*t380*t719-32*t24
     #2*t163+4*t150*t935+2*t354*t1121+8*t196*t161*t102*t+8*t463*t82+12*t
     #129*t311+8*t242*t1133
      t2394 = -8*t78*t968*t56-32*t558*t719-13*t112*t1210+8*t78*t894-48*t
     #582*t402-12*t347*t471+t652*t53*t56-9*t85*t86+2*t150*t1249-48*t50*t
     #151-6*t159*t513-12*t510*t74-24*t67*t80*t51*t289+4*t112*t633+8*t463
     #*t1250
      t2415 = r*t66
      t2426 = 22*t150*t1126-30*t112*t395+8*t704*t183-12*t479*t1063-16*t6
     #85*t598-4*t354*t221+6*t241*t1295+4*t526*t633+8*t335*t402-6*t412*t3
     #29+2*t2415*t221-124*t50*t1063+16*t187*t184-108*t220*t1126+12*t187*
     #t1521
      t2459 = -48*t50*t1121+8*t85*t2239+8*t85*t633+54*t50*t1126+16*t182*
     #t598-4*t2415*l4*t221+164*t85*t539-4*t463*t74-6*t407*t724-2*t535*t3
     #29-36*t112*t427-56*t50*t872-140*t129*t1063+12*t160*t513+4*t123*t88
     #3
      t2494 = -16*t636*t203+4*t704*t80+6*t85*t235-8*t150*t1041+42*t50*t1
     #306-16*t67*t82-8*t516*t80-24*t510*t617+24*t795*t72*t943+2*t401*t52
     #*t56-8*t61*t369-2*t555*t724-24*t545*t93-12*t305*t1295-4*t304*t1133
     #+8*t548*t1908
      s34spinc2 = 8*t21*t34*t*(t9+t23-t36+1)*t43*t44*cphi+(t1729+t1070+t
     #1253+t1106+t319+t2494+t1147+t264+t1697+t1393+t405+t1290+t1467+t196
     #4+t446+t493+t1930+t1219+t580+t1360+t542+t712+t750+t1762+t149+t871+
     #t1496+t216+t951+t1031+t1531+t1796+t1830+t2155+t1663+t2223+t1184+t2
     #426+t2459+t2358+t1562+t1326+t365+t2323+t788+t2188+t2028+t1631+t142
     #9+t2091+t827+t992+t671+t1863+t622+t2290+t1996+t2394+t908+t2258+t15
     #96+t1895+t2060+t2124)*t17/(t10+t7-t12+r-t8+t11)*t43*t44/(t4+t23+l1
     #-t36)
      else
      t3 = r**2
      t4 = 1+t3
      t7 = l2-1
      t8 = 1/r
      t12 = dsqrt(dabs(-l2*l4*t*t4*(1-l4)*t7*t8))
      t14 = r*t
      t17 = (1-t14)*(r-t)
      t21 = 1-l5
      t25 = dsqrt(dabs(l5*l6*t17*t*t4*(1-l6)*t21))
      t29 = t3*t+t
      t31 = 1/t29/2
      t34 = 1/(1+t14)
      t39 = t3**2
      t40 = t**2
      t41 = t40*t
      t42 = t39*t41
      t48 = dsqrt(dabs(-2*l2*t7*t14*t4))
      t49 = t48*l4
      t50 = dsqrt(2.D0)
      t56 = dsqrt(dabs(t17*t*t4*t21*l5))
      t64 = t56*t48*t40
      t69 = t40**2
      t70 = t69*r
      t72 = t39*t3
      t73 = t72*t48
      t74 = t73*l4
      t75 = t50*t56
      t76 = t75*t40
      t83 = t3*t48
      t84 = t83*l4
      t94 = t39**2
      t95 = t94*l5
      t98 = t3*r
      t99 = t39*t98
      t100 = t99*l2
      t103 = t94*t3
      t108 = 8*t42-8*t49*t50*t56*l6*t40+12*t39*l6*t50*t64+2*t41*t3-2*t70
     #+4*t74*t76-2*t48*t50*t56*t40+12*t84*t76+4*l6*t50*t64+4*t72*l6*t50*
     #t64-8*t95*t41+16*t100*t69+2*t103*t41+4*t70*l2
      t109 = t94*r
      t113 = t69*t98
      t116 = t39*r
      t117 = t69*t116
      t122 = t41*t72
      t125 = l2*t41
      t128 = t39*t48
      t131 = t109*l5
      t140 = r*l5
      t146 = t72*l5
      t149 = t3*l5
      t152 = 4*t109*l2*t69+16*t113*l2+24*t117*l2-16*t42*l2-24*t122*l2-4*
     #t125*t3-6*t128*t76+2*t131*t69-4*t103*l2*t41-8*t39*l5*t41+2*t140*t6
     #9+8*t98*l5*t69-12*t146*t41-2*t149*t41
      t157 = t103*l5
      t167 = t75*l6*t40
      t172 = t128*l4
      t187 = l2*t69
      t190 = 8*t99*l5*t69-2*t157*t41+12*t116*l5*t69-16*t94*l2*t41-8*t74*
     #t167-24*t84*t167-24*t172*t167+4*t49*t76-6*t83*t76-2*t109*t69+12*t1
     #72*t76+8*t94*t41+4*t149*t125-4*t131*t187
      t216 = t69*l5
      t223 = 24*t146*t125-12*t117+12*l6*t3*t50*t64+12*t122-8*t113-2*t73*
     #t76+16*t39*l2*l5*t41-4*t140*t187+4*t157*t125+16*t95*t125-16*t113*l
     #2*l5-24*l2*t116*t216-16*t100*t216-8*t99*t69
      s34spinc2 = -16*t12*t25*t*t4*t31*t34*cphi-(t108+t152+t190+t223)*t8
     #/t29*t31*t34/t4
      endif
      RETURN
      END

	


CCC   These are the ln(1-z)/z, polylog(2,x)/x, etc. that are needed 
CCC   in the RV pieces

      double precision function s13e0(x)
      IMPLICIT double precision(s,t)
      double precision x,xbuf

      xbuf = 0.1d0
      if (x.gt.xbuf) then
         s13e0 = dlog(1.0d0-x)/x
      else
      t2 = x**2
      t3 = t2*x
      t4 = t2**2
      t5 = t4*t3
      t6 = t4**2
      t7 = t6*t5
      t8 = t6**2
      t9 = t8**2
      t17 = t4*x
      t19 = t4*t2
      t23 = t6*x
      t25 = -1-x/2-t9*t7/48-t9*t8/49-t2/3-t3/4-t4/5-t17/6-t19/7-t5/8-t6/
     #9-t23/10
      t26 = t6*t2
      t28 = t6*t3
      t30 = t6*t4
      t32 = t6*t17
      t34 = t6*t19
      t48 = -t26/11-t28/12-t30/13-t32/14-t34/15-t7/16-t8/17-t8*x/18-t8*t
     #2/19-t8*t3/20-t8*t4/21-t8*t17/22
      t73 = -t8*t19/23-t8*t5/24-t8*t6/25-t8*t23/26-t8*t26/27-t8*t28/28-t
     #8*t30/29-t8*t32/30-t8*t34/31-t8*t7/32-t9/33-t9*x/34
      t100 = -t9*t2/35-t9*t3/36-t9*t4/37-t9*t17/38-t9*t19/39-t9*t5/40-t9
     #*t6/41-t9*t23/42-t9*t26/43-t9*t28/44-t9*t30/45-t9*t32/46-t9*t34/47
      s13e0 = t25+t48+t73+t100
   
      endif
      END

      
      double precision function s13e1(x)
      IMPLICIT double precision(s,t)
      double precision x

      t2 = x**2
      t3 = t2**2
      t4 = t3**2
      t5 = t4*t2
      t6 = t4**2
      t7 = t6**2
      t10 = t2*x
      t11 = t4*t10
      t14 = t4*t3
      t17 = t3*x
      t18 = t4*t17
      t21 = t3*t2
      t22 = t4*t21
      t27 = t3*t10
      t28 = t4*t27
      t34 = 1+x/4+t7*t5/1849+t7*t11/1936+t7*t14/2025+t7*t18/2116+t7*t22/
     #2209+t7*t6/2401+t7*t28/2304+t2/9+t10/16+t3/25
      t39 = t4*x
      t48 = t17/36+t21/49+t27/64+t4/81+t39/100+t5/121+t11/144+t14/169+t1
     #8/196+t22/225+t28/256+t6/289
      t74 = t6*x/324+t6*t2/361+t6*t10/400+t6*t3/441+t6*t17/484+t6*t21/52
     #9+t6*t27/576+t6*t4/625+t6*t39/676+t6*t5/729+t6*t11/784+t6*t14/841
      t100 = t6*t18/900+t6*t22/961+t6*t28/1024+t7/1089+t7*x/1156+t7*t2/1
     #225+t7*t10/1296+t7*t3/1369+t7*t17/1444+t7*t21/1521+t7*t27/1600+t7*
     #t4/1681+t7*t39/1764
      s13e1 = t34+t48+t74+t100

      END
      
         
      double precision function s13e2(x)
      IMPLICIT double precision(s,t)
      double precision x

      t1 = x**2
      t2 = t1**2
      t3 = t2**2
      t4 = t3*x
      t6 = t3*t1
      t9 = t1*x
      t10 = t3*t9
      t12 = t3*t2
      t14 = t2*x
      t15 = t3*t14
      t17 = t2*t9
      t18 = t3*t17
      t20 = t3**2
      t22 = t2*t1
      t23 = t3*t22
      t27 = 1+t4/1000+t6/1331+x/8+t10/1728+t12/2197+t15/2744+t18/4096+t2
     #0/4913+t23/3375+t1/27+t9/64
      t47 = t2/125+t14/216+t22/343+t17/512+t3/729+t20*x/5832+t20*t1/6859
     #+t20*t9/8000+t20*t2/9261+t20*t14/10648+t20*t22/12167+t20*t17/13824
      t65 = t20**2
      t73 = t20*t3/15625+t20*t4/17576+t20*t6/19683+t20*t10/21952+t20*t12
     #/24389+t20*t15/27000+t20*t23/29791+t20*t18/32768+t65/35937+t65*x/3
     #9304+t65*t1/42875+t65*t9/46656
      t100 = t65*t2/50653+t65*t14/54872+t65*t22/59319+t65*t17/64000+t65*
     #t3/68921+t65*t4/74088+t65*t6/79507+t65*t10/85184+t65*t12/91125+t65
     #*t15/97336+t65*t23/103823+t65*t18/110592+t65*t20/117649
      s13e2 = t27+t47+t73+t100

      END

      double precision function s13e3(x)
      IMPLICIT double precision(s,t)
      double precision x

      t2 = x**2
      t4 = t2*x
      t6 = t2**2
      t8 = t6*x
      t10 = t6*t2
      t12 = t6*t4
      t14 = t6**2
      t16 = t14*x
      t18 = t14*t2
      t20 = t14*t4
      t22 = 1+x/16+t2/81+t4/256+t6/625+t8/1296+t10/2401+t12/4096+t14/656
     #1+t16/10000+t18/14641+t20/20736
      t23 = t14*t6
      t25 = t14*t8
      t27 = t14*t10
      t29 = t14*t12
      t31 = t14**2
      t47 = t23/28561+t25/38416+t27/50625+t29/65536+t31/83521+t31*x/1049
     #76+t31*t2/130321+t31*t4/160000+t31*t6/194481+t31*t8/234256+t31*t10
     #/279841+t31*t12/331776
      t65 = t31**2
      t73 = t31*t14/390625+t31*t16/456976+t31*t18/531441+t31*t20/614656+
     #t31*t23/707281+t31*t25/810000+t31*t27/923521+t31*t29/1048576+t65/1
     #185921+t65*x/1336336+t65*t2/1500625+t65*t4/1679616
      t100 = t65*t6/1874161+t65*t8/2085136+t65*t10/2313441+t65*t12/25600
     #00+t65*t14/2825761+t65*t16/3111696+t65*t18/3418801+t65*t20/3748096
     #+t65*t23/4100625+t65*t25/4477456+t65*t27/4879681+t65*t29/5308416+t
     #65*t31/5764801
      s13e3 = t22+t47+t73+t100

      END

CCC   These are more RV terms needed

      double precision function sbqe0(l1)
      IMPLICIT double precision(s,t)
      double precision l1,l1buf

      l1buf=0.1d0
      if (l1.gt.l1buf) then
         sbqe0=dlog(1.0d0-l1)/l1
      else
      t2 = l1**2
      t4 = t2*l1
      t6 = t2**2
      t8 = t6*l1
      t10 = t6*t2
      t12 = t6*t4
      t14 = t6**2
      t31 = t14**2
      t41 = -t14*t2/11-t14*t4/12-t14*t6/13-t14*t8/14-t14*t10/15-t14*t12/
     #16-t31/17-t31*l1/18-t31*t2/19-t31*t4/20-t31*t6/21
      sbqe0 = -1-l1/2-t2/3-t4/4-t6/5-t8/6-t10/7-t12/8-t14/9-t14*l1/10+t4
     #1
      endif
      
      END

      double precision function sbqe1(l1)
      IMPLICIT double precision(s,t)
      double precision l1,l1buf

      l1buf=0.1d0
      if (l1.gt.l1buf) then
         sbqe1=-0.5d0*(dlog(1.0d0-l1))**2/l1
      else
      t2 = l1**2
      t4 = t2*l1
      t6 = t2**2
      t8 = t6*l1
      t10 = t6*t2
      t12 = t6*t4
      t14 = t6**2
      t30 = t14**2
      sbqe1 = -l1/2-t2/2-11.D0/24.D0*t4-5.D0/12.D0*t6-137.D0/360.D0*t8-7
     #.D0/20.D0*t10-363.D0/1120.D0*t12-761.D0/2520.D0*t14-7129.D0/25200.
     #D0*t14*l1-671.D0/2520.D0*t14*t2-83711.D0/332640.D0*t14*t4-6617.D0/
     #27720.D0*t14*t6-1145993.D0/5045040.D0*t14*t8-1171733.D0/5405400.D0
     #*t14*t10-1195757.D0/5765760.D0*t14*t12-143327.D0/720720.D0*t30-421
     #42223.D0/220540320.D0*t30*l1-751279.D0/4084080.D0*t30*t2-275295799
     #.D0/1551950400.D0*t30*t4-55835135.D0/325909584.D0*t30*t6
      endif

      END

      double precision function sbqe2(l1)
      IMPLICIT double precision(s,t)
      double precision l1,l1buf

      l1buf=0.1d0
      if (l1.gt.l1buf) then
         sbqe2=1.0d0/6.0d0*(dlog(1.0d0-l1))**3/l1
      else
      t1 = l1**2
      t3 = t1*l1
      t5 = t1**2
      t13 = t5**2
      sbqe2 = -t1/6-t3/4-7.D0/24.D0*t5-5.D0/16.D0*l1*t5-29.D0/90.D0*t1*t
     #5-469.D0/1440.D0*t5*t3-29531.D0/90720.D0*t13-1303.D0/4032.D0*t13*l
     #1-16103.D0/50400.D0*t13*t1
      endif

      END
      
      double precision function zqe0(l1,z,l3)
      IMPLICIT double precision(s,t)
      double precision l1,z,l3

      if (l1.lt.(0.5d0)) then
      t2 = l1**2
      t3 = t2**2
      t4 = t3**2
      t5 = t4*l1
      t6 = z**2
      t7 = t6**2
      t8 = t7*t6
      t12 = t3*t2
      t13 = t6*z
      t19 = t2*l1
      t20 = t3*t19
      t23 = t7*z
      t28 = t7*t13
      t37 = -2520-21168*t5*t8-280*t4+7200*t12*t13-504*t3-504*t3*t7-11025
     #*t20*t7+6615*t20*t23+2205*t20*z+315*t20*t28-2205*t20*t8+11025*t20*
     #t13+1680*t2*z
      t52 = t7**2
      t67 = -840*t2*t6-19600*t4*t7+2240*t4*z-7840*t4*t8+2240*t4*t28-7840
     #*t4*t6+15680*t4*t23-280*t4*t52+15680*t4*t13+1890*t19*z+630*t19*t13
     #-1890*t19*t6+2016*t3*z+2016*t3*t13
      t71 = t3*l1
      t99 = -3024*t3*t6+420*t71*t23+2100*t71*z+4200*t71*t13-4200*t71*t6-
     #2100*t71*t7+31752*t5*t23+21168*t5*t13+252*t5*t52*z-2268*t5*t52-317
     #52*t5*t7+9072*t5*t28+2268*t5*z-9072*t5*t6
      t121 = 1260*z*l1+2160*t12*t23-5400*t12*t7-5400*t12*t6+2160*t12*z-3
     #60*t12*t8-360*t12-315*t20-1260*l1-840*t2-6615*t20*t6-252*t5-420*t7
     #1-630*t19
      zqe0 = (-1+z)*(t37+t67+t99+t121)/2520    
      else
      t4 = dlog(1-l1+z*l1)
      zqe0 = -1/l1*t4
      endif

      END

      double precision function zqe1(l1,z,l3)
      IMPLICIT double precision(s,t)
      double precision l1,z,l3

      if (l1.lt.(0.5d0)) then
      t1 = l1**2
      t2 = t1**2
      t3 = t2*t1
      t4 = z**2
      t5 = t3*t4
      t7 = t2**2
      t8 = t7*z
      t10 = t4*z
      t11 = t4**2
      t12 = t11*t10
      t13 = t12*t7
      t17 = t2*l1
      t19 = t7*l1
      t20 = t19*t10
      t22 = t11*t4
      t23 = t7*t22
      t25 = t19*t22
      t29 = t1*l1
      t30 = t2*t29
      t32 = t11**2
      t35 = t17*z
      t37 = t19*t11
      t39 = -57.D0/20.D0*t5+481.D0/280.D0*t8+t13/2-t3*t12/49-137.D0/360.
     #D0*t17+2509.D0/210.D0*t20-5.D0/2.D0*t23-z-191.D0/24.D0*t25+t17*t22
     #/36-363.D0/1120.D0*t30+t30*t32/64+77.D0/60.D0*t35-1879.D0/120.D0*t
     #37
      t41 = t30*t22
      t43 = t2*t4
      t45 = t4*t7
      t47 = t19*t12
      t49 = t17*t10
      t51 = t19*z
      t53 = t11*z
      t54 = t7*t53
      t56 = t30*z
      t59 = t19*t32
      t61 = t30*t10
      t63 = t32*z
      t66 = t17*t4
      t68 = -l1/2-t41/2-7.D0/6.D0*t43-341.D0/70.D0*t45+17.D0/6.D0*t47+3.
     #D0/2.D0*t49+4609.D0/2520.D0*t51+73.D0/12.D0*t54+223.D0/140.D0*t56-
     #7.D0/20.D0*t3-t59/2+319.D0/60.D0*t61-t7*t63/81-47.D0/24.D0*t66
      t70 = t2*t10
      t72 = t30*t11
      t77 = t7*t10
      t95 = -2520-5400*t5+630*t29*t10+315*t30*t12+2240*t8+2240*t13-420*t
     #17+21168*t20-7840*t23-21168*t25-504*t2*t11-315*t30+2100*t35
      t113 = -31752*t37-1260*l1-2205*t41+1260*z*l1-3024*t43-7840*t45+907
     #2*t47+4200*t49+2268*t51+252*t19*t63+15680*t54+2205*t56-840*t1*t4-3
     #60*t3
      t123 = t17*t11
      t126 = t29*z
      t128 = t7*t11
      t130 = t30*t4
      t134 = t19*t4
      t136 = -2268*t59+11025*t61-4200*t66+2016*t70+420*t17*t53-11025*t72
     #+15680*t77-2100*t123-630*t29+1890*t126-19600*t128-6615*t130-280*t7
     #*t32-9072*t134
      t137 = t30*t53
      t139 = t29*t4
      t143 = t3*t10
      t149 = t19*t53
      t151 = t3*z
      t153 = t2*z
      t155 = t3*t11
      t157 = t3*t53
      t159 = t1*z
      t161 = 6615*t137-1890*t139-280*t7-840*t1+7200*t143-360*t3*t22-504*
     #t2-252*t19+31752*t149+2160*t151+2016*t153-5400*t155+2160*t157+1680
     #*t159
      t165 = dlog(z)
      t180 = t70/2-107.D0/24.D0*t72+t19*t32*t4/100+743.D0/90.D0*t77-(-1+
     #z)*(t95+t113+t136+t161)*t165/2520-t123/2+t4*l1/4-11.D0/24.D0*t29+5
     #.D0/6.D0*t126-533.D0/60.D0*t128-t2*t53/25-153.D0/40.D0*t130-3349.D
     #0/560.D0*t134-t1*t10/9
      t196 = 13.D0/6.D0*t137-t139/2-761.D0/2520.D0*t7-t1/2+37.D0/12.D0*t
     #143+t29*t11/16-5.D0/12.D0*t2-7129.D0/25200.D0*t19+55.D0/4.D0*t149+
     #29.D0/20.D0*t151+13.D0/12.D0*t153-11.D0/6.D0*t155+t157/2+t159/2
      zqe1 = t39+t68+t180+t196
      else
      t2 = 1/l1
      t3 = 1/l3*t2
      t4 = z*l1
      t7 = dlog(1-l1+t4-t4*l3)
      t10 = dlog(1-l1+t4)
      t12 = dlog(z)
      t15 = t10**2
      zqe1 = t3*t7-t3*t10+t2*t12*t10-t2*t15/2
      endif
      
      END

      double precision function zqe2(l1,z,l3)
      IMPLICIT double precision(s,t)
      double precision l1,z,l3

      if (l1.lt.(0.5d0)) then
      t1 = l1**2
      t2 = t1*l1
      t3 = t1**2
      t4 = t3*t2
      t5 = z**2
      t6 = t5**2
      t7 = t4*t6
      t9 = t6*z
      t10 = t4*t9
      t12 = t4*z
      t14 = t6**2
      t15 = t4*t14
      t17 = t5*z
      t18 = t17*t4
      t20 = t4*t5
      t22 = t1*t17
      t24 = t3*t1
      t26 = t3**2
      t28 = t26*l1
      t30 = t26*t6
      t32 = t26*z
      t34 = 3.D0/4.D0*t7-t10/6-319.D0/315.D0*t12+t15/512-179.D0/120.D0*t
     #18+235.D0/144.D0*t20-t22/27+29.D0/90.D0*t24+29531.D0/90720.D0*t26+
     #1303.D0/4032.D0*t28+251.D0/120.D0*t30-349.D0/288.D0*t32
      t35 = t6*t5
      t36 = t26*t35
      t53 = -2520-11025*t7+6615*t10+2205*t12+11025*t18-6615*t20-360*t24-
     #280*t26-252*t28-19600*t30+2240*t32-280*t26*t14-7840*t36
      t55 = t3*l1
      t60 = t24*t9
      t62 = t24*t6
      t64 = t24*t5
      t66 = t24*z
      t70 = t24*t17
      t72 = t26*t5
      t74 = t26*t9
      t78 = t6*t17
      t79 = t28*t78
      t81 = -504*t3+420*t55*t9+1260*z*l1+2160*t60-5400*t62-5400*t64+2160
     #*t66-360*t24*t35+7200*t70-7840*t72+15680*t74-1260*l1-315*t4+9072*t
     #79
      t83 = t28*z
      t85 = t28*t5
      t89 = t4*t35
      t91 = t1*z
      t95 = t26*t78
      t99 = t2*t5
      t103 = t3*t17
      t105 = t55*t6
      t107 = t14*z
      t110 = t28*t14
      t112 = 2268*t83-9072*t85+315*t4*t78-2205*t89+1680*t91-840*t1*t5+22
     #40*t95+630*t2*t17-1890*t99-504*t3*t6+2016*t103-2100*t105+252*t28*t
     #107-2268*t110
      t113 = t26*t17
      t115 = t2*z
      t117 = t3*z
      t119 = t3*t5
      t121 = t55*z
      t123 = t55*t17
      t125 = t55*t5
      t127 = t28*t9
      t129 = t28*t17
      t131 = t28*t6
      t133 = t28*t35
      t138 = 15680*t113+1890*t115+2016*t117-3024*t119+2100*t121+4200*t12
     #3-4200*t125+31752*t127+21168*t129-31752*t131-21168*t133-840*t1-630
     #*t2-420*t55
      t142 = dlog(z)
      t143 = t142**2
      t163 = 107.D0/24.D0*t7-13.D0/6.D0*t10-223.D0/140.D0*t12-t15/64-319
     #.D0/60.D0*t18+153.D0/40.D0*t20+t22/9+7.D0/20.D0*t24+761.D0/2520.D0
     #*t26+7129.D0/25200.D0*t28+533.D0/60.D0*t30-481.D0/280.D0*t32+5.D0/
     #2.D0*t36
      t177 = 5.D0/12.D0*t3-t60/2+11.D0/6.D0*t62+57.D0/20.D0*t64-29.D0/20
     #.D0*t66-37.D0/12.D0*t70+341.D0/70.D0*t72-73.D0/12.D0*t74+z+l1/2+36
     #3.D0/1120.D0*t4-17.D0/6.D0*t79-4609.D0/2520.D0*t83+3349.D0/560.D0*
     #t85
      t187 = t28*t14*t5
      t189 = t26*t107
      t193 = t2*t6
      t195 = t3*t9
      t198 = t89/2-t91/2-t95/2+t99/2-t103/2+t105/2+t110/2-t187/100+t189/
     #81-743.D0/90.D0*t113-5.D0/6.D0*t115-t193/16+t195/25-13.D0/12.D0*t1
     #17
      t200 = t55*t35
      t210 = l1*t5
      t212 = t24*t78
      t216 = 7.D0/6.D0*t119-t200/36-77.D0/60.D0*t121-3.D0/2.D0*t123+47.D
     #0/24.D0*t125-55.D0/4.D0*t127-2509.D0/210.D0*t129+1879.D0/120.D0*t1
     #31+191.D0/24.D0*t133+t1/2-t210/4+t212/49+11.D0/24.D0*t2+137.D0/360
     #.D0*t55
      t223 = t36/6+7.D0/24.D0*t3+(-1+z)*(t53+t81+t112+t138)*t143/5040+t6
     #2/6+119.D0/120.D0*t64-29.D0/36.D0*t66-5.D0/8.D0*t70+(t163+t177+t19
     #8+t216)*t142+6077.D0/2520.D0*t72-7.D0/8.D0*t74-z+469.D0/1440.D0*t4
      t237 = -t79/6-1571.D0/1120.D0*t83+1063.D0/320.D0*t85+t187/1000-t18
     #9/729-23.D0/8.D0*t113-t115/6+t193/64-t195/125-3.D0/8.D0*t117+t119/
     #6+t200/216
      t249 = -71.D0/120.D0*t121-t123/6+t125/2-67.D0/24.D0*t127-12287.D0/
     #2520.D0*t129+665.D0/144.D0*t131+t133+t1/6+t210/8-t212/343+t2/4+5.D
     #0/16.D0*t55
      zqe2 = t34+t223+t237+t249
      else
      t1 = 1/l3
      t2 = 1/l1
      t3 = t1*t2
      t4 = z*l1
      t7 = dlog(1-l1+t4-t4*l3)
      t8 = t7**2
      t11 = dlog(l3)
      t12 = t11*t2
      t14 = dlog(z)
      t18 = t14**2
      t21 = dlog(1-l1+t4)
      t24 = t21**2
      zqe2 = t3*t8/2+(-t12*t1-t3*t14)*t7-t2*t18*t21/2+(t2*t24/2+t3*t21)*
     #t14-t2*t24*t21/6-t3*t24/2+t12*t1*t21
      endif

      END

      double precision function sqe1(l1)
      IMPLICIT double precision(s,t)
      double precision l1

      if (l1.lt.(0.5d0)) then
      t2 = l1**2
      t4 = t2*l1
      t6 = t2**2
      t14 = t6**2
      sqe1 = -1-l1/4-t2/9-t4/16-t6/25-t6*l1/36-t6*t2/49-t6*t4/64-t14/81-
     #t14*l1/100-t14*t2/121
      else
      t2 = dlog(1-l1)
      t3 = l1**2
      t4 = t3*l1
      t5 = t3**2
      t6 = t5*t4
      t9 = t5**2
      t10 = t9*t3
      t14 = t9*l1
      t16 = t5*l1
      t18 = t5*t3
      t20 = 0.3141592653589793D1**2
      t33 = 3586319.D0/254016.D0+40123.D0/84.D0*t6*t2-7381.D0/2520.D0*t2
     #*t10+12012823.D0/141120.D0*t9-108194591.D0/6350400.D0*t14-1435891.
     #D0/2016.D0*t16+15411317.D0/30240.D0*t18+77*t20*t16-74471.D0/168.D0
     #*t2*t3-11393.D0/12.D0*t2*t18+1968329.D0/1270080.D0*t10-26739449.D0
     #/52920.D0*t4+15344789.D0/21600.D0*t5-103767707.D0/1270080.D0*l1-55
     #991.D0/2520.D0*t2-13501741.D0/52920.D0*t6
      t66 = -11.D0/6.D0*t20+55.D0/2.D0*t20*t6+69851.D0/504.D0*t2*l1+55.D
     #0/6.D0*t20*l1-55.D0/2.D0*t20*t3+76781.D0/84.D0*t2*t4+55*t20*t4-781
     #67.D0/60.D0*t2*t5-8959.D0/56.D0*t2*t9-55.D0/6.D0*t20*t9+80939.D0/2
     #520.D0*t14*t2-77*t20*t5+79091.D0/60.D0*t2*t16-55*t20*t18+11.D0/6.D
     #0*t20*t14-t20*t10/6+3928917.D0/15680.D0*t3
      sqe1 = t33+t66
      endif
      
      END

      double precision function sqe2(l1)
      IMPLICIT double precision(s,t)
      double precision l1

      if (l1.lt.(0.5d0)) then
      t2 = l1**2
      t4 = t2*l1
      t6 = t2**2
      t14 = t6**2
      sqe2 = -1+l1/8+7.D0/54.D0*t2+19.D0/192.D0*t4+113.D0/1500.D0*t6+127
     #.D0/2160.D0*t6*l1+323.D0/6860.D0*t6*t2+691.D0/17920.D0*t6*t4+6569.
     #D0/204120.D0*t14+6877.D0/252000.D0*t14*l1+78671.D0/3354120.D0*t14*
     #t2
      else
      t1 = l1**2
      t2 = t1**2
      t3 = t2**2
      t5 = 0.3141592653589793D1**2
      t6 = t2*l1
      t10 = dlog(1-l1)
      t11 = t2*t1
      t16 = t10**2
      t23 = t3*t1
      t28 = t1*l1
      t29 = t2*t28
      t32 = -103374326407.D0/3200601600.D0-29353051921.D0/118540800.D0*t
     #3-79091.D0/360.D0*t6*t5+49039331.D0/30240.D0*t10*t11+11393.D0/72.D
     #0*t5*t11-78167.D0/120.D0*t16*t2+79091.D0/120.D0*t16*t6-8959.D0/112
     #.D0*t16*t3-7381.D0/5040.D0*t16*t23-11393.D0/24.D0*t16*t11+40123.D0
     #/168.D0*t16*t29
      t44 = t3*l1
      t48 = 197010683477.D0/266716800.D0*t29+665654530691.D0/3200601600.
     #D0*l1+55991.D0/15120.D0*t5-55991.D0/5040.D0*t16-74471.D0/336.D0*t1
     #6*t1-80017096283.D0/118540800.D0*t1+212667113.D0/6350400.D0*t10+37
     #3705386829.D0/266716800.D0*t28+76781.D0/168.D0*t16*t28+79612018463
     #9.D0/16003008000.D0*t44-351904421.D0/6350400.D0*t10*t44
      t70 = 74471.D0/1008.D0*t1*t5-69851.D0/3024.D0*l1*t5-55972381.D0/25
     #4016.D0*t10*l1-80939.D0/15120.D0*t5*t44+8959.D0/336.D0*t5*t3+10254
     #3859.D0/141120.D0*t10*t1-111751769417.D0/76204800.D0*t11-109101219
     #881.D0/54432000.D0*t2-10084132.D0/6615.D0*t10*t28-76781.D0/504.D0*
     #t5*t28+7381.D0/15120.D0*t5*t23
      t91 = 32160403.D0/6350400.D0*t10*t23+47446259.D0/21600.D0*t10*t2-4
     #0123.D0/504.D0*t5*t29+78167.D0/360.D0*t5*t2-21677681.D0/26460.D0*t
     #10*t29+1147469347.D0/564480.D0*t6-37620007.D0/16800.D0*t10*t6+7769
     #893.D0/28224.D0*t10*t3+80939.D0/5040.D0*t16*t44-14528236349.D0/320
     #0601600.D0*t23+69851.D0/1008.D0*t16*l1
      sqe2 = t32+t48+t70+t91
      endif
      
      END

      double precision function Diquad71(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+l1*(1-r*t)/2.0d0/r/(r+t)
      A1 = l1*l2*(1+u)**2/2.0d0
      A2 = (1-l1/2.0d0)*(1-l2)*Kp*r*(1-l1*Km/2.0d0/Kp)
      Diquad71 = l2*s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END

      double precision function Diquad721(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+(1-l1*l2/2.0d0)*(1-r*t)/r/(r+t)
      A1 = (1-l1*l2/2.0d0)*(1+u)**2
      A2 = l1/2.0d0*(1-l2)*Kp*r*(1-(1-l1*l2/2.0d0)*Km/Kp)
      Diquad721 = s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END
            
      double precision function Diquad722(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+(1-l1/2.0d0)*(1-r*t)/r/(r+t)
      A1 = (1-l1/2.0d0)*l2*(1+u)**2      
      A2 = 1/2.0d0*(1-l1*l2)*Kp*r*(1-(1-l1/2.0d0)*Km/Kp)
      Diquad722 = l2*s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END

      double precision function Dsquar311(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+(1-l1*l2)*(1-r*t)/r/(r+t)
      A1 = (1-l1*l2)*(1+u)**2
      A2 = l1*(1-l2)*Kp*r*(1-(1-l1*l2)*Km/Kp)
      Dsquar311 = s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END

      double precision function Dsquar312(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+(1-l1)*(1-r*t)/r/(r+t)
      A1 = l2*(1-l1)*(1+u)**2
      A2 = (1-l1*l2)*Kp*r*(1-(1-l1)*Km/Kp)
      Dsquar312 = l2*s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END

      double precision function Dsquar321(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+l1*l2*(1-r*t)/r/(r+t)
      A1 = l1*(1-l2)*(1+u)**2
      A2 = (1-l1*l2)*Kp*r*(1-l1*l2*Km/Kp)
      Dsquar321 = s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END

      double precision function Dsquar322(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,r,t,u,A1,A2,Kp,Km
      u = (l3-l3*z+z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      Kp = (r+t)*(1+r*t)
      Km = (r-t)*(1-r*t)
      s1 = 1+l1*(1-r*t)/r/(r+t)
      A1 = (1-l1*l2)*(1+u)**2
      A2 = l2*(1-l1)*Kp*r*(1-l1*Km/Kp)
      Dsquar322 = l2*s1+(A1+A2+2.0d0*(2.0d0*l4-1)*dsqrt(A1*A2))
     & /Kp/r/s1
      RETURN
      END
  
      double precision function k1(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4
      t3 = l2*l3
      t11 = dsqrt(t3*(1-l2)*(1-l3))
      k1 = 1/(1-(1-z)*l1*(l2+l3-2*t3-2*(2*l4-1)*t11))
      RETURN
      END

      double precision function s13f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4
      s13f = -l1*(1-z)*(1-l2)
      RETURN
      END

      double precision function s23f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4
      s23f = -l1*l2*(1-z)
      RETURN
      END

      double precision function s14f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,k1
      external k1
      t1 = k1(z,l1,l2,l3,l4)
      s14f = -t1*(1-z)*(1-l1)*(1-l3)
      RETURN
      END

      double precision function s24f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,k1
      external k1
      t1 = k1(z,l1,l2,l3,l4)
      s24f = -t1*(1-z)*(1-l1)*l3
      RETURN
      END

      
      double precision function zs13f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4
      zs13f = -l1*(1-l2)
      RETURN
      END

      double precision function zs23f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4
      zs23f = -l1*l2
      RETURN
      END

      double precision function zs14f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,k1
      external k1
      t1 = k1(z,l1,l2,l3,l4)
      zs14f = -t1*(1-l1)*(1-l3)
      RETURN
      END

      double precision function zs24f(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,k1
      external k1
      t1 = k1(z,l1,l2,l3,l4)
      zs24f = -t1*(1-l1)*l3
      RETURN
      END

      double precision function s3vf(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,s14f,s24f
      external s14f,s24f
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      s3vf = 1+t1+t2
      RETURN
      END

      double precision function s4vf(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,s13f,s23f
      external s13f,s23f
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      s4vf = 1+t1+t2
      RETURN
      END

      double precision function f1(y,z,l1)
      IMPLICIT double precision(s,t)
      double precision y,z,l1,r,t,u
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t3 = 1/r
      t7 = r*t
      f1 = t2*t3/(r+t)/(1+t7)*(1-l1*(r-t)*(1-t7)*t3/(1+z))
      RETURN
      END

      double precision function f2(y,z)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t6 = (1+u)**2
      f2 = r*(r+t)*(1+r*t)/t6
      RETURN
      END

      double precision function f3(y,z)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = r*t
      f3 = (r-t)*(1-t2)/(r+t)/(1+t2)
      RETURN
      END

      double precision function s34r2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      s34r2 = t2/r/(r+t)/(1+r*t)
      RETURN
      END

      double precision function s23r2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s23r2 = 1+l1*(1-r*t)/r/(r+t)
      RETURN
      END

      double precision function s1vr2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s1vr2 = 1-l1*r*(1-r*t)/(r+t)
      RETURN
      END

      double precision function s2vr2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s2vr2 = 1-l1*(r-t)/r/(1+r*t)
      RETURN
      END

      double precision function s23rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s23rf = -l2*(1+l1*(1-r*t)/r/(r+t))
      RETURN
      END

      double precision function s34rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      s34rf = l1*t2/r/(r+t)/(1+r*t)
      RETURN
      END

      double precision function s1vrf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      s1vrf = -1+l1*r*(1-r*t)/(r+t)
      RETURN
      END

      double precision function s13rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = r+t
      t4 = 1/r/t2
      t5 = r*t
      t6 = 1+t5
      t8 = 1-t5
      t14 = 1-l2
      t17 = (1+u)**2
      t19 = 1-l1
      t27 = r*t2*t6-l1*r*(r-t)*t8
      t36 = dsqrt(l1*l2*t19*t14*t17*t27)
      s13rf = -t4/t6/(1+l1*t8*t4)*(l1*t14*t17+l2*t19*t27+2*(2*l4-1)*t36)
      RETURN
      END

      double precision function s3vrf(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,s13rf,s23rf,s34rf
      external s13rf,s23rf,s34rf
      t1 = 1-l3
      t2 = 1-z
      t4 = s13rf(z,l1,l2,l3,l4)
      t7 = s23rf(z,l1,l2,l3,l4)
      t10 = t2**2
      t11 = s34rf(z,l1,l2,l3,l4)
      s3vrf = z-t1*t2*t4-l3*t2*t7-l3*t1*t10*t11
      RETURN
      END

      double precision function s4vrf(z,l1,l2,l3,l4)
      IMPLICIT double precision(s,t)
      double precision z,l1,l2,l3,l4,s13rf,s23rf
      external s13rf,s23rf
      t2 = 1-z
      t4 = s13rf(z,l1,l2,l3,l4)
      t7 = s23rf(z,l1,l2,l3,l4)
      s4vrf = 1+(1-l3)*t2*t4+l3*t2*t7
      RETURN
      END

      double precision function p1a(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t3 = r-t
      t5 = r*t
      t6 = 1-t5
      t7 = 1/r
      t16 = r+t
      p1a = -t2*(1-l1*t3*t6*t7/(1+z))*(1+l1*t6*t7/t16)/(l1*t2+(1-l1)*(r*
     #t16*(1+t5)-l1*r*t3*t6))
      RETURN
      END

      double precision function p2a(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t6 = (1+u)**2
      t7 = t6**2
      p2a = r*(r+t)*(1+r*t)/t7
      RETURN
      END

      double precision function f13a(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t4 = (1+u)**2
      t9 = r*t
      f13a = l1*(1-l2)*t4+(1-l1)*(r*(r+t)*(1+t9)-l1*r*(r-t)*(1-t9))
      RETURN
      END

      double precision function g13a(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      a = l1*(1-l2)*(1+u)**2+(1-l1)*(r*(r+t)*(1+r*t)-l1*r*(r-t)*(1-r*t))
      b = l2*(1-l1)*(r*(r+t)*(1+r*t)-l1*r*(r-t)*(1-r*t))
      t4 = dsqrt(a*b)
      g13a = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s23I13a(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t3 = r*t
      t4 = 1-t3
      t7 = r+t
      s23I13a = t2*(1+l1*t4/r/t7)/(l1*t2+(1-l1)*(r*t7*(1+t3)-l1*r*(r-t)*
     #t4))
      RETURN
      END

      double precision function s23I13arf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t3 = (1+u)**2
      t4 = r*t
      t5 = 1-t4
      t8 = r+t
      s23I13arf = -l1*l2*t3*(1+l1*t5/r/t8)/(l1*t3+(1-l1)*(r*t8*(1+t4)-l1
     #*r*(r-t)*t5))
      RETURN
      END

      double precision function s13I13arf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,num,g13a
      external g13a
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = r+t
      t4 = 1/r/t2
      t5 = r*t
      t6 = 1+t5
      t8 = 1-t5
      t16 = (1+u)**2
      t17 = l1*t16
      t19 = (1-l2)**2
      num = -t4/t6/(1+l1*t8*t4)*t17*t19*(t17+(1-l1)*(r*t2*t6-l1*r*(r-t)*
     #t8))
      s13I13arf = num/g13a(z,l1,l2,y,l4)
      RETURN
      END

      double precision function s24I13a(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t4 = (1+u)**2
      t7 = r+t
      t9 = r*t
      t14 = 1-t9
      t18 = (1-l1)*(r*t7*(1+t9)-l1*r*(r-t)*t14)
      s24I13a = (l1*(1-l2)*t4+t18)*(1+l1*t14/r/t7)/(l1*t4+t18)
      RETURN
      END

      double precision function f13a4(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t3 = 1-l1*l2/2
      t5 = (1+u)**2
      t9 = r*t
      f13a4 = t3*t5/2+l1*(r*(r+t)*(1+t9)-t3*r*(r-t)*(1-t9))/2
      RETURN
      END

      double precision function g13a4(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t3 = 1-l1*l2/2
      t5 = (1+u)**2
      t9 = r*t
      a = t3*t5/2+l1*(r*(r+t)*(1+t9)-t3*r*(r-t)*(1-t9))/2
      t6 = r*t
      b = l1*(1-l2/2)*(r*(r+t)*(1+t6)-(1-l1*l2/2)*r*(r-t)*(1-t6))/2
      t4 = dsqrt(a*b)
      g13a4 = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s24I13a4(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = l1*l2
      t3 = 1-t1/2
      t5 = (1+u)**2
      t6 = t3*t5
      t7 = r+t
      t9 = r*t
      t14 = 1-t9
      t17 = r*t7*(1+t9)-t3*r*(r-t)*t14
      s24I13a4 = (t6/2+l1*t17/2)*(1+t3*t14/r/t7)/(t6+t1*t17/2)
      RETURN
      END

      double precision function s13I13a4rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,den13,pref,g13a4
      external g13a4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = l1*l2
      t3 = 1-t1/2
      t5 = (1+u)**2
      t9 = r*t
      den13 = t3*t5+t1*(r*(r+t)*(1+t9)-t3*r*(r-t)*(1-t9))/2
      t3 = 1-l1*l2/2
      t5 = 1/r
      t8 = 1/(r+t)
      t9 = r*t
      t20 = (1+u)**2
      pref = -t3*l2*t5*t8/(1+t9)/(1+t3*(1-t9)*t5*t8)*t20/4
      s13I13a4rf = den13*pref/g13a4(z,l1,l2,y,l4)
      RETURN
      END

      double precision function f13a5(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = 1-l1/2
      t5 = (1+u)**2
      t9 = r*t
      f13a5 = l2*t2*t5/2+r*(r+t)*(1+t9)/2-t2*r*(r-t)*(1-t9)/2
      RETURN
      END

      double precision function g13a5(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = 1-l1/2
      t5 = (1+u)**2
      t9 = r*t
      a = l2*t2*t5/2+r*(r+t)*(1+t9)/2-t2*r*(r-t)*(1-t9)/2
      t6 = r*t
      b = (1-l1*l2/2)*(r*(r+t)*(1+t6)-(1-l1/2)*r*(r-t)*(1-t6))/2
      t4 = dsqrt(a*b)
      g13a5 = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s24I13a5(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = 1-l1/2
      t5 = (1+u)**2
      t7 = r+t
      t9 = r*t
      t11 = r*t7*(1+t9)
      t14 = 1-t9
      t16 = t2*r*(r-t)*t14
      s24I13a5 = (l2*t2*t5/2+t11/2-t16/2)*(1+t2*t14/r/t7)/(t2*t5+l1*(t11
     #-t16)/2)
      RETURN
      END

      double precision function s13I13a5rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,den13,pref,g13a5
      external g13a5
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = 1-l1/2
      t4 = (1+u)**2
      t8 = r*t
      den13 = t2*t4+l1*(r*(r+t)*(1+t8)-t2*r*(r-t)*(1-t8))/2
      t2 = 1-l1/2
      t4 = l2**2
      t5 = 1/r
      t9 = 1/(r+t)
      t10 = r*t
      t21 = (1+u)**2
      pref = -t2*l1*t4*t5*t9/(1+t10)/(1+t2*(1-t10)*t5*t9)*t21/4
      s13I13a5rf = den13*pref/g13a5(z,l1,l2,y,l4)
      RETURN
      END

      double precision function p2b(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      p2b = 1/t2
      RETURN
      END

      double precision function f13b(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t8 = r*t
      f13b = l1*t2+l2*(1-l1)*(r*(r+t)*(1+t8)-l1*r*(r-t)*(1-t8))
      RETURN
      END

      double precision function g13b(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      a = l1*(1+u)**2+l2*(1-l1)*(r*(r+t)*(1+r*t)-l1*r*(r-t)*(1-r*t))
      b = l1*(1-l2)*(1+u)**2
      t4 = dsqrt(a*b)
      g13b = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s23I13b(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t3 = l1*t2
      t4 = 1-l1
      t6 = r+t
      t8 = r*t
      t13 = 1-t8
      t16 = r*t6*(1+t8)-l1*r*(r-t)*t13
      s23I13b = (t3+l2*t4*t16)*(1+l1*t13/r/t6)/(t3+t4*t16)
      RETURN
      END

      double precision function s24I13b(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = r+t
      t3 = r*t
      t8 = 1-t3
      t11 = r*t1*(1+t3)-l1*r*(r-t)*t8
      t20 = (1+u)**2
      s24I13b = t11*(1+l1*t8/r/t1)/(l1*t20+(1-l1)*t11)
      RETURN
      END

      double precision function s13I13brf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,num,g13b
      external g13b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = r+t
      t4 = 1/r/t2
      t5 = r*t
      t6 = 1+t5
      t8 = 1-t5
      t15 = l2**2
      t16 = 1-l1
      t24 = r*t2*t6-l1*r*(r-t)*t8
      t26 = (1+u)**2
      num = -t4/t6/(1+l1*t8*t4)*t15*t16*t24*(l1*t26+t16*t24)
      s13I13brf = num/g13b(z,l1,l2,y,l4)
      RETURN
      END

      double precision function f13b1(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = l1*l2
      t6 = r*t
      t18 = (1+u)**2
      f13b1 = (1-t1/2)*(r*(r+t)*(1+t6)-t1*r*(r-t)*(1-t6)/2)/2+l1*t18/2
      RETURN
      END

      double precision function g13b1(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = l1*l2
      t6 = r*t
      t18 = (1+u)**2
      a = (1-t1/2)*(r*(r+t)*(1+t6)-t1*r*(r-t)*(1-t6)/2)/2+l1*t18/2
      t5 = (1+u)**2
      b = l1*(1-l2/2)*t5/2
      t4 = dsqrt(a*b)
      g13b1 = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s23I13b1(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,f13b1
      external f13b1
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = f13b1(z,l1,l2,y)
      t2 = l1*l2
      t3 = r*t
      t4 = 1-t3
      t7 = r+t
      t15 = (1+u)**2
      s23I13b1 = t1*(1+t2*t4/r/t7/2)/(t2*t15/2+(1-t2/2)*(r*t7*(1+t3)-t2*
     #r*(r-t)*t4/2))
      RETURN
      END

      double precision function s13I13b1rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,den13,pref,g13b1
      external g13b1
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = l1*l2
      t3 = (1+u)**2
      t10 = r*t
      den13 = t1*t3/2+(1-t1/2)*(r*(r+t)*(1+t10)-t1*r*(r-t)*(1-t10)/2)
      t1 = l1*l2
      t5 = r+t
      t7 = r*t
      t8 = 1+t7
      t12 = 1-t7
      t18 = 1/r
      t19 = 1/t5
      pref = -l2*(1-t1/2)*(r*t5*t8-t1*r*(r-t)*t12/2)*t18*t19/t8/(1+t1*t1
     #2*t18*t19/2)/4
      s13I13b1rf = den13*pref/g13b1(z,l1,l2,y,l4)
      RETURN
      END

      double precision function f13b2(z,l1,l2,y)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t6 = r*t
      t18 = (1+u)**2
      f13b2 = l2*(1-l1/2)*(r*(r+t)*(1+t6)-l1*r*(r-t)*(1-t6)/2)/2+t18/2
      RETURN
      END

      double precision function g13b2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,a,b
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t6 = r*t
      t18 = (1+u)**2
      a = l2*(1-l1/2)*(r*(r+t)*(1+t6)-l1*r*(r-t)*(1-t6)/2)/2+t18/2
      t5 = (1+u)**2
      b = (1-l1*l2/2)*t5/2
      t4 = dsqrt(a*b)
      g13b2 = a+b+2*(2*l4-1)*t4
      RETURN
      END

      double precision function s23I13b2(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,f13b2
      external f13b2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t1 = f13b2(z,l1,l2,y)
      t2 = r*t
      t3 = 1-t2
      t6 = r+t
      t14 = (1+u)**2
      s23I13b2 = t1*(1+l1*t3/r/t6/2)/(l1*t14/2+(1-l1/2)*(r*t6*(1+t2)-l1*
     #r*(r-t)*t3/2))
      RETURN
      END

      double precision function s13I13b2rf(z,l1,l2,y,l4)
      IMPLICIT double precision(s,t)
      double precision y,z,r,t,u,l1,l2,l4,den13,pref,g13b2
      external g13b2
      u = (y+z-y*z)/(1-y+y*z)
      r = dsqrt(u)
      t = dsqrt(z)
      t2 = (1+u)**2
      t9 = r*t
      den13 = l1*t2/2+(1-l1/2)*(r*(r+t)*(1+t9)-l1*r*(r-t)*(1-t9)/2)
      t1 = l2**2
      t5 = r+t
      t7 = r*t
      t8 = 1+t7
      t12 = 1-t7
      t21 = 1/r/t5
      pref = -l1*t1*(1-l1/2)*(r*t5*t8-l1*r*(r-t)*t12/2)*t21/t8/(1+l1*t12
     #*t21/2)/4
      s13I13b2rf = den13*pref/g13b2(z,l1,l2,y,l4)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   POLYLOGS (Babis's code)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C STEP 1: series expansions 

      double precision function sli2(x)
      double precision fin, x
      integer i 
      fin=0 
      do i=1, 13, 1
        fin = fin + (x**i)/dble(i**2)       
      end do
      sli2=fin
      end 

      double precision function sli3(x)
      double precision fin, x
      integer i 
      fin=0 
      do i=1, 13, 1
        fin = fin + (x**i)/dble(i**3)       
      end do
      sli3=fin
      end 
      
      double precision function ss12(x)
      double precision fin, x
      integer i,j,k
      fin=0
      do k=2, 13, 1 
        do i=1, k, 1
          do j=1, k, 1
            if ((i+j).eq.(k)) then 
              fin= fin+x**(i+j)/dble(i)/dble(j)/dble(i+j)/2.0d0
            end if            
          end do 
        end do
      end do 
      ss12=fin
      end 

C STEP 2: Analytic continuations 


      double precision function Li2(x) 
      double precision x, pi, zeta2
      double precision sli2
      external sli2

      pi=3.1415926535897932385D0
      zeta2=Pi**2/6.0D0
      
      if ((x.lt.(0.5)).or.(x.eq.(0.5)).AND.(x.gt.0)) then     
        Li2=sli2(x)
      end if
      
      if ((x.gt.(0.5)).AND.(x.lt.1)) THEN 
       Li2=-dlog(x)*dlog(1.0d0-x)+pi**2/6.0d0-sli2(1.0d0-x)             
      end if
      
      if (x.lt.0) then 
       Li2 = -0.5d0*dlog(1.0d0-x)**2-sli2(x/(x-1.0d0))
      end if
      
      if  (x.eq.1) then 
        Li2=zeta2
      end if
    
      if (x.eq.0) then
        Li2=0.0D0
      end if
      end


      double precision function Li3(x) 
      double precision x, pi, zeta2, zeta3
      double precision sli2, sli3, ss12
      external sli2, sli3, ss12

      pi=3.1415926535897932385D0
      zeta2=Pi**2/6.0D0
      zeta3=1.2020569031595942854D0
      if ((x.lt.(0.5)).or.(x.eq.(0.5)).AND.(x.gt.0)) then     
        Li3=sli3(x)
      end if
      
      if ((x.gt.(0.5)).AND.(x.lt.1)) THEN 
       Li3=-ss12(1.0d0-x)-0.5D0*dlog(x)**2*dlog(1.0d0-x)-dlog(x)
     & *sli2(1.0d0-x)+zeta2*dlog(x)+zeta3          
      end if
      
      if (x.lt.0) then 
       Li3 = -sli3(x/(x-1))+ss12(x/(x-1.0d0))-dlog(1.0d0-x)
     & *sli2(x/(x-1.0d0))-1.0D0/6.0D0*dlog(1.0d0-x)**3
      end if
      
      if  (x.eq.1) then 
        Li3=zeta3
      end if
      if (x.eq.0) then
        Li3=0.0D0
      end if
      end

      double precision function S12f(x) 
      double precision x, pi, zeta2, zeta3
      double precision sli2, sli3, ss12
      external sli2, sli3, ss12

      pi=3.1415926535897932385D0
      zeta2=pi**2/6.0D0
      zeta3=1.2020569031595942854D0
      if ((x.lt.(0.5)).or.(x.eq.(0.5)).AND.(x.gt.0)) then     
        S12f=ss12(x)
      end if
      
      if ((x.gt.(0.5)).AND.(x.lt.1)) THEN 
       S12f=-sli3(1.0d0-x)+0.5D0*dlog(x)*dlog(1.0d0-x)**2
     &      +dlog(1.0d0-x)*sli2(1.0d0-x)+zeta3 

      end if     
      if (x.lt.0) then 
       S12f=ss12(x/(x-1.0d0))+1.0D0/6.0D0*dlog(1.0d0-x)**3 
      end if      
      if  (x.eq.1) then 
        S12f=zeta3
      end if
      if (x.eq.0) then
        S12f=0.0D0
      end if
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

     
