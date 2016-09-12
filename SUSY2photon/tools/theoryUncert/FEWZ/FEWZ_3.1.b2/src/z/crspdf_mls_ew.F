
      DOUBLE PRECISION FUNCTION CS0z_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      double precision weakadj
c      common/weakadjst/ weakadj
      weakadj = 1d0
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0z_mls =
     #+ulum*(A1e0c1p1t1z*geu**2*gel**2+A1e0c2p1t1z*(gvusq*gvlsq
     #+gvusq*galsq+gausq*gvlsq+gausq*galsq)*weakadj**2+A1e0c3p1t1z*geu*gel*
     #gvu*gvl)+dlum*(A1e0c1p1t2z*ged**2*gel**2+A1e0c2p1t
     #2z*(gvdsq*gvlsq+gvdsq*galsq+gadsq*gvlsq+gadsq*galsq)*weakadj**2+A1e0c
     #3p1t2z*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1z*gvugau*gvlgal*weakadj**2+A2e0c2p1t1z*geu*gel*ga
     #u*gal)+dlum*(A2e0c1p1t2z*gvdgad*gvlgal*weakadj**2+A2e0c2p1t
     #2z*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0s_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0s_mls =
     #+ulum*(A1e0c1p1t1s*(GeuSq*geu**2*gel**2+geu**2*GelSq*gel**2)+A1e0c2p1
     #t1s*GeuGel*geu**2*gel**2+A1e0c3p1t1s*(gvusq*gvlsq*GeuSq+gvusq*gvlsq*GelSq
     #+gvusq*galsq*GeuSq+gvusq*galsq*GelSq+gausq*gvlsq*GeuSq
     #+gausq*gvlsq*GelSq+gausq*galsq*GeuSq+gausq*galsq*GelSq
     #)+A1e0c4p1t1s*(gvusq*gvlsq*GeuGel+gvusq*galsq*GeuGel+gausq*
     #gvlsq*GeuGel+gausq*galsq*GeuGel)+A1e0c5p1t1s*(GeuSq*geu*gel*gvu*
     #gvl+GelSq*geu*gel*gvu*gvl)+A1e0c6p1t1s*GeuGel*geu*gel*gvu*gvl)+dlum
     #*(A1e0c1p1t2s*(GedSq*ged**2*gel**2+ged**2*GelSq*gel**2)+A1e0c2p1t
     #2s*GedGel*ged**2*gel**2+A1e0c3p1t2s*(gvdsq*gvlsq*GedSq+gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GedSq+gvdsq*galsq*GelSq+gadsq*gvlsq*GedSq
     #+gadsq*gvlsq*GelSq+gadsq*galsq*GedSq+gadsq*galsq*GelSq)
     #+A1e0c4p1t2s*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+gadsq*g
     #vlsq*GedGel+gadsq*galsq*GedGel)+A1e0c5p1t2s*(GedSq*ged*gel*gvd*g
     #vl+GelSq*ged*gel*gvd*gvl)+A1e0c6p1t2s*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1s*(gvugau*gvlgal*GeuSq+gvugau*gvlgal*GelSq
     #)+A2e0c2p1t1s*gvugau*gvlgal*GeuGel+A2e0c3p1t1s*(GeuSq*geu*gel
     #*gau*gal+GelSq*geu*gel*gau*gal)+A2e0c4p1t1s*GeuGel*geu*gel*gau*gal)
     #+dlum*(A2e0c1p1t2s*(gvdgad*gvlgal*GedSq+gvdgad
     #*gvlgal*GelSq)+A2e0c2p1t2s*gvdgad*gvlgal*GedGel+A2e0c3p1t2s*(
     #GedSq*ged*gel*gad*gal+GelSq*ged*gel*gad*gal)+A2e0c4p1t2s*GedGel*ged*gel*g
     #ad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0c_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum,ulum1,ulum2,dlum1,dlum2
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
         ulum1 = uubplum(2,iset)
         ulum2 = uubplum(3,iset)
         dlum1 = ddbplum(2,iset)
         dlum2 = ddbplum(3,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
         ulum1 = uubmlum(2,iset)
         ulum2 = uubmlum(3,iset)
         dlum1 = ddbmlum(2,iset)
         dlum2 = ddbmlum(3,iset)
      endif
      CS0c_mls =
     #+ulum*(A1e0c1p1t1c*GeuSq*geu**2*gel**2+A1e0c2p1t1c*geu**2*GelSq*gel**2
     #+A1e0c3p1t1c*(GeuGel*geu**2*gel**2+gvusq*gvlsq*GeuGel+gvusq*galsq*GeuGel
     #+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel+GeuGel*geu*gel*g
     #vu*gvl)+A1e0c4p1t1c*(gvusq*gvlsq*GeuSq+gvusq*galsq*GeuSq+gau
     #sq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1e0c5p1t1c*(gvusq*gvlsq*GelSq
     #+gvusq*galsq*GelSq+gausq*gvlsq*GelSq+gausq*galsq*GelSq
     #)+A1e0c6p1t1c*GeuSq*geu*gel*gvu*gvl+A1e0c7p1t1c*GelSq*geu*gel*gvu*gv
     #l)+dlum*(A1e0c1p1t2c*GedSq*ged**2*gel**2+A1e0c2p1t2c*ged
     #**2*GelSq*gel**2+A1e0c3p1t2c*(GedGel*ged**2*gel**2+gvdsq*gvlsq*GedGel+gvdsq
     #*galsq*GedGel+gadsq*gvlsq*GedGel+gadsq*galsq*GedGel+GedGel*ged*gel*g
     #vd*gvl)+A1e0c4p1t2c*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq
     #+gadsq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2c*(gvdsq
     #*gvlsq*GelSq+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*
     #galsq*GelSq)+A1e0c6p1t2c*GedSq*ged*gel*gvd*gvl+A1e0c7p1t2c*GelSq*ged*gel
     #*gvd*gvl)
     #+ulum2*(A1e0c1p2t1c*GeuSq*geu**2*gel**2+A1e0c2p2t1c*geu**2*GelSq*gel**2
     #+A1e0c3p2t1c*(GeuGel*geu**2*gel**2+gvusq*gvlsq*GeuGel+gvusq*galsq*GeuGel
     #+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel+GeuGel*geu*gel*g
     #vu*gvl)+A1e0c4p2t1c*(gvusq*gvlsq*GeuSq+gvusq*galsq*GeuSq+g
     #ausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1e0c5p2t1c*(gvusq*gvlsq*GelSq
     #+gvusq*galsq*GelSq+gausq*gvlsq*GelSq+gausq*galsq*GelSq
     #)+A1e0c6p2t1c*GeuSq*geu*gel*gvu*gvl+A1e0c7p2t1c*GelSq*geu*gel*gvu*
     #gvl)+dlum2*(A1e0c1p2t2c*GedSq*ged**2*gel**2+A1e0c2p2t2
     #c*ged**2*GelSq*gel**2+A1e0c3p2t2c*(GedGel*ged**2*gel**2+gvdsq*gvlsq*GedGel+g
     #vdsq*galsq*GedGel+gadsq*gvlsq*GedGel+gadsq*galsq*GedGel+GedGel*ged*gel*g
     #vd*gvl)+A1e0c4p2t2c*(gvdsq*gvlsq*GedSq+gvdsq*gal
     #sq*GedSq+gadsq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p2t2c*
     #(gvdsq*gvlsq*GelSq+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+ga
     #dsq*galsq*GelSq)+A1e0c6p2t2c*GedSq*ged*gel*gvd*gvl+A1e0c7p2t2c*GelSq*ged*gel
     #*gvd*gvl)
     #+ulum1*(A1e0c1p3t1c*GeuSq*geu**2*gel**2+A1e0c2p3t1c*geu**2*GelSq*gel**2
     #+A1e0c3p3t1c*(GeuGel*geu**2*gel**2+gvusq*gvlsq*GeuGel+gvusq*galsq*GeuGel
     #+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel+GeuGel*geu*gel*g
     #vu*gvl)+A1e0c4p3t1c*(gvusq*gvlsq*GeuSq+gvusq*galsq*GeuSq+g
     #ausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1e0c5p3t1c*(gvusq*gvlsq*GelSq
     #+gvusq*galsq*GelSq+gausq*gvlsq*GelSq+gausq*galsq*GelSq
     #)+A1e0c6p3t1c*GeuSq*geu*gel*gvu*gvl+A1e0c7p3t1c*GelSq*geu*gel*gvu*
     #gvl)+dlum1*(A1e0c1p3t2c*GedSq*ged**2*gel**2+A1e0c2p3t2
     #c*ged**2*GelSq*gel**2+A1e0c3p3t2c*(GedGel*ged**2*gel**2+gvdsq*gvlsq*GedGel+g
     #vdsq*galsq*GedGel+gadsq*gvlsq*GedGel+gadsq*galsq*GedGel+GedGel*ged*gel*g
     #vd*gvl)+A1e0c4p3t2c*(gvdsq*gvlsq*GedSq+gvdsq*gal
     #sq*GedSq+gadsq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p3t2c*
     #(gvdsq*gvlsq*GelSq+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+ga
     #dsq*galsq*GelSq)+A1e0c6p3t2c*GedSq*ged*gel*gvd*gvl+A1e0c7p3t2c*GelSq*ged*gel
     #*gvd*gvl)
     #+ulum*(A2e0c1p1t1c*gvugau*gvlgal*GeuSq+A2e0c2p1t1c*gvu
     #gau*gvlgal*GelSq+A2e0c3p1t1c*(gvugau*gvlgal*GeuGel+GeuGel*geu*gel*g
     #au*gal)+A2e0c4p1t1c*GeuSq*geu*gel*gau*gal+A2e0c5p1t1c*GelSq*geu*gel
     #*gau*gal)+dlum*(A2e0c1p1t2c*gvdgad*gvlgal*GedSq
     #+A2e0c2p1t2c*gvdgad*gvlgal*GelSq+A2e0c3p1t2c*(gvdgad*gvlgal*GedGel
     #+GedGel*ged*gel*gad*gal)+A2e0c4p1t2c*GedSq*ged*gel*gad*gal+A2e0
     #c5p1t2c*GelSq*ged*gel*gad*gal)
     #+ulum2*(A2e0c1p2t1c*gvugau*gvlgal*GeuSq+A2e0c2p2t1c*g
     #vugau*gvlgal*GelSq+A2e0c3p2t1c*(gvugau*gvlgal*GeuGel+GeuGel*geu*gel*g
     #au*gal)+A2e0c4p2t1c*GeuSq*geu*gel*gau*gal+A2e0c5p2t1c*GelSq*geu*gel
     #*gau*gal)+dlum2*(A2e0c1p2t2c*gvdgad*gvlgal*GedSq
     #+A2e0c2p2t2c*gvdgad*gvlgal*GelSq+A2e0c3p2t2c*(gvdgad*gvlgal*GedGel
     #+GedGel*ged*gel*gad*gal)+A2e0c4p2t2c*GedSq*ged*gel*gad*gal
     #+A2e0c5p2t2c*GelSq*ged*gel*gad*gal)
     #+ulum1*(A2e0c1p3t1c*gvugau*gvlgal*GeuSq+A2e0c2p3t1c*g
     #vugau*gvlgal*GelSq+A2e0c3p3t1c*(gvugau*gvlgal*GeuGel+GeuGel*geu*gel*g
     #au*gal)+A2e0c4p3t1c*GeuSq*geu*gel*gau*gal+A2e0c5p3t1c*GelSq*geu*gel
     #*gau*gal)+dlum1*(A2e0c1p3t2c*gvdgad*gvlgal*GedSq
     #+A2e0c2p3t2c*gvdgad*gvlgal*GelSq+A2e0c3p3t2c*(gvdgad*gvlgal*GedGel
     #+GedGel*ged*gel*gad*gal)+A2e0c4p3t2c*GedSq*ged*gel*gad*gal
     #+A2e0c5p3t2c*GelSq*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0v_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0v_mls =
     #+ulum*(A1e0c1p1t1v*(GeuSq*geu**2*gel**2+geu**2*GelSq*gel**2)+A1e0c2p1
     #t1v*GeuGel*geu**2*gel**2+A1e0c3p1t1v*(gvusq*gvlsq*GeuSq+gvusq*gvlsq*GelSq
     #+gvusq*galsq*GeuSq+gvusq*galsq*GelSq+gausq*gvlsq*GeuSq
     #+gausq*gvlsq*GelSq+gausq*galsq*GeuSq+gausq*galsq*GelSq
     #)+A1e0c4p1t1v*(gvusq*gvlsq*GeuGel+gvusq*galsq*GeuGel+gausq*
     #gvlsq*GeuGel+gausq*galsq*GeuGel)+A1e0c5p1t1v*(GeuSq*geu*gel*gvu*
     #gvl+GelSq*geu*gel*gvu*gvl)+A1e0c6p1t1v*GeuGel*geu*gel*gvu*gvl)
     #+dlum*(A1e0c1p1t2v*(GedSq*ged**2*gel**2+ged**2*GelSq*gel**2
     #)+A1e0c2p1t2v*GedGel*ged**2*gel**2+A1e0c3p1t2v*(gvdsq*gvlsq*GedSq+gvd
     #sq*gvlsq*GelSq+gvdsq*galsq*GedSq+gvdsq*galsq*GelSq+gadsq
     #*gvlsq*GedSq+gadsq*gvlsq*GelSq+gadsq*galsq*GedSq+gadsq*ga
     #lsq*GelSq)+A1e0c4p1t2v*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+gadsq*g
     #vlsq*GedGel+gadsq*galsq*GedGel)+A1e0c5p1t2v*(GedSq*ged*gel*gvd*g
     #vl+GelSq*ged*gel*gvd*gvl)+A1e0c6p1t2v*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1v*(gvugau*gvlgal*GeuSq+gvugau*gvlgal*GelSq
     #)+A2e0c2p1t1v*gvugau*gvlgal*GeuGel+A2e0c3p1t1v*(GeuSq*geu*gel
     #*gau*gal+GelSq*geu*gel*gau*gal)+A2e0c4p1t1v*GeuGel*geu*gel*gau*gal)
     #+dlum*(A2e0c1p1t2v*(gvdgad*gvlgal*GedSq+gvdgad
     #*gvlgal*GelSq)+A2e0c2p1t2v*gvdgad*gvlgal*GedGel+A2e0c3p1t2v*(
     #GedSq*ged*gel*gad*gal+GelSq*ged*gel*gad*gal)+A2e0c4p1t2v*GedGel*ged*gel*g
     #ad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n_mls =
     #+ulum*(A1e0c1p1t1n*GeuSq*geu**2*gel**2+A1e0c2p1t1n*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n*GeuGel*geu**2*gel**2+A1e0c4p1t1n*(gvusq*gvlsq*GeuSq+gvusq
     #*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1e0c5
     #p1t1n*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n*(gvusq*gvlsq*GeuGel+gvusq
     #*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel)+A1e0
     #c7p1t1n*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n*GelSq*geu*gel*gvu*gvl+A1e0c9p1
     #t1n*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1t2n*GedSq*ged**2
     #*gel**2+A1e0c2p1t2n*ged**2*GelSq*gel**2+A1e0c3p1t2n*GedGel*ged**2*gel**2+A1e0c
     #4p1t2n*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gadsq*gvlsq*GedSq
     #+gadsq*galsq*GedSq)+A1e0c5p1t2n*(gvdsq*gvlsq*GelSq+gvdsq
     #*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq)+A1e0c6p
     #1t2n*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+gadsq*gvlsq*GedGel
     #+gadsq*galsq*GedGel)+A1e0c7p1t2n*GedSq*ged*gel*gvd*gvl+A1e0c8
     #p1t2n*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n*gvugau*gvlgal*GeuSq+A2e0c2p1t1n*gvu
     #gau*gvlgal*GelSq+A2e0c3p1t1n*gvugau*gvlgal*GeuGel+A2e0c4p1t1
     #n*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n*GelSq*geu*gel*gau*gal+A2e0c6p1t1n*GeuGel*geu*gel*g
     #au*gal)+dlum*(A2e0c1p1t2n*gvdgad*gvlgal*GedSq
     #+A2e0c2p1t2n*gvdgad*gvlgal*GelSq+A2e0c3p1t2n*gvdgad
     #*gvlgal*GedGel+A2e0c4p1t2n*GedSq*ged*gel*gad*gal+A2e0c5p1t2n*GelSq*ged*gel
     #*gad*gal+A2e0c6p1t2n*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n1_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n1_mls =
     #+ulum*(A1e0c1p1t1n1*GeuSq*geu**2*gel**2+A1e0c2p1t1n1*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n1*GeuGel*geu**2*gel**2+A1e0c4p1t1n1*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n1*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n1*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n1*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n1*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n1*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n1*GedSq*ged**2*gel**2+A1e0c2p1t2n1*ged**2*GelSq*gel**2+A1e0c3p1t2n1*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n1*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n1*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n1*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n1*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n1*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n1*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n1*gvugau*gvlgal*GeuSq+A2e0c2p1t1n1*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n1*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n1*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n1*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n1*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n1*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n1*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #1*gvdgad*gvlgal*GedGel+A2e0c4p1t2n1*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n1*GelSq*ged*gel*gad*gal+A2e0c6p1t2n1*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n2_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n2_mls =
     #+ulum*(A1e0c1p1t1n2*GeuSq*geu**2*gel**2+A1e0c2p1t1n2*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n2*GeuGel*geu**2*gel**2+A1e0c4p1t1n2*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n2*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n2*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n2*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n2*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n2*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n2*GedSq*ged**2*gel**2+A1e0c2p1t2n2*ged**2*GelSq*gel**2+A1e0c3p1t2n2*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n2*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n2*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n2*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n2*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n2*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n2*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n2*gvugau*gvlgal*GeuSq+A2e0c2p1t1n2*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n2*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n2*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n2*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n2*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n2*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n2*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #2*gvdgad*gvlgal*GedGel+A2e0c4p1t2n2*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n2*GelSq*ged*gel*gad*gal+A2e0c6p1t2n2*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n3_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n3_mls =
     #+ulum*(A1e0c1p1t1n3*GeuSq*geu**2*gel**2+A1e0c2p1t1n3*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n3*GeuGel*geu**2*gel**2+A1e0c4p1t1n3*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n3*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n3*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n3*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n3*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n3*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n3*GedSq*ged**2*gel**2+A1e0c2p1t2n3*ged**2*GelSq*gel**2+A1e0c3p1t2n3*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n3*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n3*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n3*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n3*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n3*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n3*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n3*gvugau*gvlgal*GeuSq+A2e0c2p1t1n3*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n3*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n3*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n3*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n3*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n3*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n3*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #3*gvdgad*gvlgal*GedGel+A2e0c4p1t2n3*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n3*GelSq*ged*gel*gad*gal+A2e0c6p1t2n3*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n4_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n4_mls =
     #+ulum*(A1e0c1p1t1n4*GeuSq*geu**2*gel**2+A1e0c2p1t1n4*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n4*GeuGel*geu**2*gel**2+A1e0c4p1t1n4*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n4*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n4*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n4*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n4*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n4*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n4*GedSq*ged**2*gel**2+A1e0c2p1t2n4*ged**2*GelSq*gel**2+A1e0c3p1t2n4*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n4*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n4*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n4*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n4*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n4*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n4*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n4*gvugau*gvlgal*GeuSq+A2e0c2p1t1n4*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n4*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n4*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n4*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n4*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n4*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n4*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #4*gvdgad*gvlgal*GedGel+A2e0c4p1t2n4*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n4*GelSq*ged*gel*gad*gal+A2e0c6p1t2n4*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n5_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n5_mls =
     #+ulum*(A1e0c1p1t1n5*GeuSq*geu**2*gel**2+A1e0c2p1t1n5*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n5*GeuGel*geu**2*gel**2+A1e0c4p1t1n5*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n5*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n5*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n5*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n5*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n5*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n5*GedSq*ged**2*gel**2+A1e0c2p1t2n5*ged**2*GelSq*gel**2+A1e0c3p1t2n5*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n5*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n5*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n5*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n5*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n5*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n5*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n5*gvugau*gvlgal*GeuSq+A2e0c2p1t1n5*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n5*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n5*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n5*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n5*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n5*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n5*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #5*gvdgad*gvlgal*GedGel+A2e0c4p1t2n5*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n5*GelSq*ged*gel*gad*gal+A2e0c6p1t2n5*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n6_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n6_mls =
     #+ulum*(A1e0c1p1t1n6*GeuSq*geu**2*gel**2+A1e0c2p1t1n6*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n6*GeuGel*geu**2*gel**2+A1e0c4p1t1n6*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n6*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n6*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n6*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n6*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n6*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n6*GedSq*ged**2*gel**2+A1e0c2p1t2n6*ged**2*GelSq*gel**2+A1e0c3p1t2n6*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n6*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n6*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n6*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n6*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n6*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n6*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n6*gvugau*gvlgal*GeuSq+A2e0c2p1t1n6*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n6*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n6*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n6*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n6*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n6*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n6*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #6*gvdgad*gvlgal*GedGel+A2e0c4p1t2n6*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n6*GelSq*ged*gel*gad*gal+A2e0c6p1t2n6*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n7_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n7_mls =
     #+ulum*(A1e0c1p1t1n7*GeuSq*geu**2*gel**2+A1e0c2p1t1n7*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n7*GeuGel*geu**2*gel**2+A1e0c4p1t1n7*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n7*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n7*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n7*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n7*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n7*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n7*GedSq*ged**2*gel**2+A1e0c2p1t2n7*ged**2*GelSq*gel**2+A1e0c3p1t2n7*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n7*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n7*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n7*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n7*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n7*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n7*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n7*gvugau*gvlgal*GeuSq+A2e0c2p1t1n7*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n7*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n7*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n7*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n7*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n7*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n7*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #7*gvdgad*gvlgal*GedGel+A2e0c4p1t2n7*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n7*GelSq*ged*gel*gad*gal+A2e0c6p1t2n7*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n8_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n8_mls =
     #+ulum*(A1e0c1p1t1n8*GeuSq*geu**2*gel**2+A1e0c2p1t1n8*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n8*GeuGel*geu**2*gel**2+A1e0c4p1t1n8*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n8*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n8*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n8*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n8*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n8*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n8*GedSq*ged**2*gel**2+A1e0c2p1t2n8*ged**2*GelSq*gel**2+A1e0c3p1t2n8*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n8*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n8*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n8*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n8*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n8*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n8*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n8*gvugau*gvlgal*GeuSq+A2e0c2p1t1n8*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n8*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n8*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n8*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n8*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n8*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n8*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #8*gvdgad*gvlgal*GedGel+A2e0c4p1t2n8*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n8*GelSq*ged*gel*gad*gal+A2e0c6p1t2n8*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n9_mls(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0n9_mls =
     #+ulum*(A1e0c1p1t1n9*GeuSq*geu**2*gel**2+A1e0c2p1t1n9*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n9*GeuGel*geu**2*gel**2+A1e0c4p1t1n9*(gvusq*gvlsq*GeuSq+
     #gvusq*galsq*GeuSq+gausq*gvlsq*GeuSq+gausq*galsq*GeuSq)+A1
     #e0c5p1t1n9*(gvusq*gvlsq*GelSq+gvusq*galsq*GelSq+gausq*gvlsq*GelSq
     #+gausq*galsq*GelSq)+A1e0c6p1t1n9*(gvusq*gvlsq*GeuGel
     #+gvusq*galsq*GeuGel+gausq*gvlsq*GeuGel+gausq*galsq*GeuGel
     #)+A1e0c7p1t1n9*GeuSq*geu*gel*gvu*gvl+A1e0c8p1t1n9*GelSq*geu*gel*gvu*gvl+
     #A1e0c9p1t1n9*GeuGel*geu*gel*gvu*gvl)+dlum*(A1e0c1p1
     #t2n9*GedSq*ged**2*gel**2+A1e0c2p1t2n9*ged**2*GelSq*gel**2+A1e0c3p1t2n9*GedGel*ged**2*gel**2
     #+A1e0c4p1t2n9*(gvdsq*gvlsq*GedSq+gvdsq*galsq*GedSq+gad
     #sq*gvlsq*GedSq+gadsq*galsq*GedSq)+A1e0c5p1t2n9*(gvdsq*gvlsq*GelSq
     #+gvdsq*galsq*GelSq+gadsq*gvlsq*GelSq+gadsq*galsq*GelSq
     #)+A1e0c6p1t2n9*(gvdsq*gvlsq*GedGel+gvdsq*galsq*GedGel+g
     #adsq*gvlsq*GedGel+gadsq*galsq*GedGel)+A1e0c7p1t2n9*GedSq*ged*gel
     #*gvd*gvl+A1e0c8p1t2n9*GelSq*ged*gel*gvd*gvl+A1e0c9p1t2n9*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n9*gvugau*gvlgal*GeuSq+A2e0c2p1t1n9*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n9*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n9*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n9*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n9*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n9*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n9*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #9*gvdgad*gvlgal*GedGel+A2e0c4p1t2n9*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n9*GelSq*ged*gel*gad*gal+A2e0c6p1t2n9*GedGel*ged*gel*gad*gal)
      RETURN
      END
