
      double precision invmassmin,invmassmax,tmassmin,tmassmax,
     & ZpTmin,ZpTmax,ZYmin,ZYmax,
     & leppTmin,leppTmax,aleppTmin,aleppTmax,softleppTmin,softleppTmax,
     & hardleppTmin,hardleppTmax,lepetamin,lepetamax,
     & alepetamin,alepetamax,softlepetamin,softlepetamax,
     & hardlepetamin,hardlepetamax,
     & jet1ptmin,
     & lldRmin,lldPhimin,lldPhimax,ljdRmin,jjdRmin,Rsep,
     & minptj4obs,maxetaj4obs,zY_CScutoff

      double precision minptr4obs,maxetar4obs,
     & lrdRmin,lrdRrecb
C     & epsiso,delRiso

      integer lepabseta,alepabseta,softlepabseta,hardlepabseta
      integer minnumjet,maxnumjet
      integer minnumpht,maxnumpht

      character*4 jetalgo

      double precision ptbuf

      COMMON/cutcontrol/ invmassmin,invmassmax,tmassmin,tmassmax,
     & ZpTmin,ZpTmax,ZYmin,ZYmax,
     & leppTmin,leppTmax,aleppTmin,aleppTmax,
     & softleppTmin,softleppTmax,hardleppTmin,hardleppTmax,
     & lepabseta,lepetamin,lepetamax,alepetamin,alepetamax,alepabseta,
     & softlepabseta,softlepetamin,softlepetamax,
     & hardlepabseta,hardlepetamin,hardlepetamax,
     & jet1ptmin,minptj4obs,maxetaj4obs,
     & lldRmin,lldPhimin,lldPhimax,ljdRmin,jjdRmin,
     & zY_CScutoff,Rsep,
     & minptr4obs,maxetar4obs,
     & lrdRmin,lrdRrecb,
C     & epsiso,delRiso
     & ptbuf,
     & jetalgo,
     & minnumjet,maxnumjet,
     & minnumpht,maxnumpht
