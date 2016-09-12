      DOUBLE PRECISION FUNCTION CS0w(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum,blum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
         blum = bbbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
         blum = bbbmlum(1,iset)
      endif
      dlum = dlum-blum
      CS0w = ulum*(A1e0c22p1t1w+A1e0c23p1t1w+A1e0c24p1t1w+
     #A1e0c25p1t1w+A1e0c26p1t1w+A1e0c27p1t1w+A1e
     #0c28p1t1w+A1e0c29p1t1w+A1e0c30p1t1w+A1e0c31p1t1w
     #+A1e0c32p1t1w+A1e0c33p1t1w)
     #       + dlum*(A1e0c10p1t2w+A1e0c11p1t2w+A1e0c12p1t2w+
     #A1e0c13p1t2w+A1e0c14p1t2w+A1e0c15p1t2w+A1e0c16p1
     #t2w+A1e0c17p1t2w+A1e0c18p1t2w+A1e0c19p1t2w+
     #A1e0c20p1t2w+A1e0c21p1t2w)
     #       + blum*(A1e0c1p1t3w+A1e0c2p1t3w+A1e0c3p1t3w+A1e0c
     #4p1t3w+A1e0c5p1t3w+A1e0c6p1t3w+A1e0c7p1t3w+
     #A1e0c8p1t3w+A1e0c9p1t3w+A1e0c10p1t3w+A1e
     #0c11p1t3w+A1e0c12p1t3w)
      RETURN
      END
