
      DOUBLE PRECISION FUNCTION CS0z_msv(whichlum,iset)
      IMPLICIT NONE
      INTEGER whichlum,iset
      DOUBLE PRECISION ulum,dlum
      include 'constcommon_ew.f'
      include 'crseccommon_ew.f'
      double complex weakadj
c      common/weakadjst/ weakadj
      weakadj = 1d0
      if (whichlum.eq.1) then
         ulum = uubplum(1,iset)
         dlum = ddbplum(1,iset)
      else
         ulum = uubmlum(1,iset)
         dlum = ddbmlum(1,iset)
      endif
      CS0z_msv =
     #+ulum*(A1e0c1p1t1z*geu**2*gel**2+A1e0c2p1t1z*(gvusq*gvlsq
     #+gausq*gvlsq)*weakadj**2+A1e0c3p1t1z*(gvusq*galsq+gausq*galsq)*weakadj**2+A1e0c4
     #p1t1z*geu*gel*gvu*gvl)+dlum*(A1e0c1p1t2z*ged**2*ge
     #l**2+A1e0c2p1t2z*(gvdsq*gvlsq+gadsq*gvlsq)*weakadj**2+A1e0c3p1t2z*(gvdsq
     #*galsq+gadsq*galsq)*weakadj**2+A1e0c4p1t2z*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1z*gvugau*gvlgal*weakadj**2+A2e0c2p1t1z*geu*gel*ga
     #u*gal)+dlum*(A2e0c1p1t2z*gvdgad*gvlgal*weakadj**2+A2e0c2p1t
     #2z*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0s_msv(whichlum,iset)
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
      CS0s_msv =
     #+ulum*(A1e0c1p1t1s*GeuSq*geu**2*gel**2+A1e0c2p1t1s*geu**2*GelSq*gel**2
     #+A1e0c3p1t1s*GeuGel*geu**2*gel**2+A1e0c4p1t1s*(gvusq*gvlsq*GeuSq+gausq
     #*gvlsq*GeuSq)+A1e0c5p1t1s*(gvusq*gvlsq*GelSq+gausq*gvlsq*GelSq
     #)+A1e0c6p1t1s*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel)+
     #A1e0c7p1t1s*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e0c8p1t1
     #s*gvusq*galsq*GelSq+A1e0c9p1t1s*(gvusq*galsq*GeuGel+gausq*g
     #alsq*GeuGel)+A1e0c10p1t1s*gausq*galsq*GelSq+A1e0c11p1t1s*GeuSq*geu*gel
     #*gvu*gvl+A1e0c12p1t1s*GelSq*geu*gel*gvu*gvl+A1e0c13p1t1s*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2s*GedSq*ged**2*gel**2+A
     #1e0c2p1t2s*ged**2*GelSq*gel**2+A1e0c3p1t2s*GedGel*ged**2*gel**2+A1e0c4p1t2s*(gv
     #dsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2s*(gvdsq*gvlsq*GelSq
     #+gadsq*gvlsq*GelSq)+A1e0c6p1t2s*(gvdsq*gvlsq*GedGel+
     #gadsq*gvlsq*GedGel)+A1e0c7p1t2s*(gvdsq*galsq*GedSq+gadsq*ga
     #lsq*GedSq)+A1e0c8p1t2s*gvdsq*galsq*GelSq+A1e0c9p1t2s*(gvdsq*
     #galsq*GedGel+gadsq*galsq*GedGel)+A1e0c10p1t2s*gadsq*galsq*GelSq
     #+A1e0c11p1t2s*GedSq*ged*gel*gvd*gvl+A1e0c12p1t2s*GelSq*ged*gel*gvd*
     #gvl+A1e0c13p1t2s*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1s*gvugau*gvlgal*GeuSq+A2e0c2p1t1s*gvu
     #gau*gvlgal*GelSq+A2e0c3p1t1s*gvugau*gvlgal*GeuGel+A2e0c4p1t1
     #s*GeuSq*geu*gel*gau*gal+A2e0c5p1t1s*GelSq*geu*gel*gau*gal+A2e0c6p1t1s*GeuGel*geu*gel*g
     #au*gal)+dlum*(A2e0c1p1t2s*gvdgad*gvlgal*GedSq
     #+A2e0c2p1t2s*gvdgad*gvlgal*GelSq+A2e0c3p1t2s*gvdgad
     #*gvlgal*GedGel+A2e0c4p1t2s*GedSq*ged*gel*gad*gal+A2e0c5p1t2s*GelSq*ged*gel
     #*gad*gal+A2e0c6p1t2s*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0c_msv(whichlum,iset)
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
      CS0c_msv =
     #+ulum*(A1e0c1p1t1c*GeuSq*geu**2*gel**2+A1e0c2p1t1c*(geu**2*GelSq*gel**2
     #+GeuGel*geu**2*gel**2+gvusq*gvlsq*GelSq+gvusq*gvlsq*GeuGel+gvusq
     #*galsq*GelSq+gvusq*galsq*GeuGel+gausq*gvlsq*GelSq+gausq*g
     #vlsq*GeuGel+gausq*galsq*GelSq+gausq*galsq*GeuGel+GelSq*geu*gel
     #*gvu*gvl+GeuGel*geu*gel*gvu*gvl)+A1e0c3p1t1c*(gvusq*gvlsq*GeuSq
     #+gausq*gvlsq*GeuSq)+A1e0c4p1t1c*(gvusq*galsq*GeuSq+gausq*ga
     #lsq*GeuSq)+A1e0c5p1t1c*GeuSq*geu*gel*gvu*gvl)+dlum*
     #(A1e0c1p1t2c*GedSq*ged**2*gel**2+A1e0c2p1t2c*(ged**2*GelSq*gel**2+GedSq*GelSq*ged*gel
     #+gvdsq*gvlsq*GelSq+gvdsq*gvlsq*GedGel+gvdsq*galsq*GelSq+
     #gvdsq*galsq*GedGel+gadsq*gvlsq*GelSq+gadsq*gvlsq*GedGel+g
     #adsq*galsq*GelSq+gadsq*galsq*GedGel+GelSq*ged*gel*gvd*gvl+GedGel*ged*gel*g
     #vd*gvl)+A1e0c3p1t2c*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq
     #)+A1e0c4p1t2c*(gvdsq*galsq*GedSq+gadsq*galsq*GedSq)+A1
     #e0c5p1t2c*GedSq*ged*gel*gvd*gvl)
     #+ulum2*(A1e0c1p2t1c*GeuSq*geu**2*gel**2+A1e0c2p2t1c*(geu**2*GelSq*gel**2
     #+GeuGel*geu**2*gel**2+gvusq*gvlsq*GelSq+gvusq*gvlsq*GeuGel+gvusq
     #*galsq*GelSq+gvusq*galsq*GeuGel+gausq*gvlsq*GelSq+gausq
     #*gvlsq*GeuGel+gausq*galsq*GelSq+gausq*galsq*GeuGel+GelSq*geu*gel
     #*gvu*gvl+GeuGel*geu*gel*gvu*gvl)+A1e0c3p2t1c*(gvusq*gvlsq*GeuSq
     #+gausq*gvlsq*GeuSq)+A1e0c4p2t1c*(gvusq*galsq*GeuSq+gausq*
     #galsq*GeuSq)+A1e0c5p2t1c*GeuSq*geu*gel*gvu*gvl)+dlum2
     #*(A1e0c1p2t2c*GedSq*ged**2*gel**2+A1e0c2p2t2c*(ged**2*GelSq*gel**2+GedGel*ged**2*gel**2
     #+gvdsq*gvlsq*GelSq+gvdsq*gvlsq*GedGel+gvdsq*galsq*GelSq
     #+gvdsq*galsq*GedGel+gadsq*gvlsq*GelSq+gadsq*gvlsq*GedGel
     #+gadsq*galsq*GelSq+gadsq*galsq*GedGel+GelSq*ged*gel*gvd*gvl+
     #GedGel*ged*gel*gvd*gvl)+A1e0c3p2t2c*(gvdsq*gvlsq*GedSq+gadsq*gv
     #lsq*GedSq)+A1e0c4p2t2c*(gvdsq*galsq*GedSq+gadsq*galsq*GedSq
     #)+A1e0c5p2t2c*GedSq*ged*gel*gvd*gvl)
     #+ulum1*(A1e0c1p3t1c*GeuSq*geu**2*gel**2+A1e0c2p3t1c*(geu**2*GelSq*gel**2
     #+GeuGel*geu**2*gel**2+gvusq*gvlsq*GelSq+gvusq*gvlsq*GeuGel+gvusq
     #*galsq*GelSq+gvusq*galsq*GeuGel+gausq*gvlsq*GelSq+gausq
     #*gvlsq*GeuGel+gausq*galsq*GelSq+gausq*galsq*GeuGel+GelSq*geu*gel
     #*gvu*gvl+GeuGel*geu*gel*gvu*gvl)+A1e0c3p3t1c*(gvusq*gvlsq*GeuSq
     #+gausq*gvlsq*GeuSq)+A1e0c4p3t1c*(gvusq*galsq*GeuSq+gausq*
     #galsq*GeuSq)+A1e0c5p3t1c*GeuSq*geu*gel*gvu*gvl)+dlum1
     #*(A1e0c1p3t2c*GedSq*ged**2*gel**2+A1e0c2p3t2c*(ged**2*GelSq*gel**2+GedGel*ged**2*gel**2
     #+gvdsq*gvlsq*GelSq+gvdsq*gvlsq*GedGel+gvdsq*galsq*GelSq
     #+gvdsq*galsq*GedGel+gadsq*gvlsq*GelSq+gadsq*gvlsq*GedGel
     #+gadsq*galsq*GelSq+gadsq*galsq*GedGel+GelSq*ged*gel*gvd*gvl+
     #GedGel*ged*gel*gvd*gvl)+A1e0c3p3t2c*(gvdsq*gvlsq*GedSq+gadsq*gv
     #lsq*GedSq)+A1e0c4p3t2c*(gvdsq*galsq*GedSq+gadsq*galsq*GedSq
     #)+A1e0c5p3t2c*GedSq*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1c*gvugau*gvlgal*GeuSq+A2e0c2p1t1c*(gv
     #ugau*gvlgal*GelSq+gvugau*gvlgal*GeuGel+GelSq*geu*gel*gau*gal+GeuGel*geu*gel*g
     #au*gal)+A2e0c3p1t1c*GeuSq*geu*gel*gau*gal)+dlum
     #*(A2e0c1p1t2c*gvdgad*gvlgal*GedSq+A2e0c2p1t2c*(gvdgad*g
     #vlgal*GelSq+gvdgad*gvlgal*GedGel+GelSq*ged*gel*gad*gal+GedGel*ged*gel*g
     #ad*gal)+A2e0c3p1t2c*GedSq*ged*gel*gad*gal)
     #+ulum2*(A2e0c1p2t1c*gvugau*gvlgal*GeuSq+A2e0c2p2t1c*(
     #gvugau*gvlgal*GelSq+gvugau*gvlgal*GeuGel+GelSq*geu*gel*gau*gal+
     #GeuGel*geu*gel*gau*gal)+A2e0c3p2t1c*GeuSq*geu*gel*gau*gal)+dlum2
     #*(A2e0c1p2t2c*gvdgad*gvlgal*GedSq+A2e0c2p2t2c*(gvdga
     #d*gvlgal*GelSq+gvdgad*gvlgal*GedGel+GelSq*ged*gel*gad*gal+GedGel*ged*gel*g
     #ad*gal)+A2e0c3p2t2c*GedSq*ged*gel*gad*gal)
     #+ulum1*(A2e0c1p3t1c*gvugau*gvlgal*GeuSq+A2e0c2p3t1c*(
     #gvugau*gvlgal*GelSq+gvugau*gvlgal*GeuGel+GelSq*geu*gel*gau*gal+
     #GeuGel*geu*gel*gau*gal)+A2e0c3p3t1c*GeuSq*geu*gel*gau*gal)+dlum1
     #*(A2e0c1p3t2c*gvdgad*gvlgal*GedSq+A2e0c2p3t2c*(gvdga
     #d*gvlgal*GelSq+gvdgad*gvlgal*GedGel+GelSq*ged*gel*gad*gal+GedGel*ged*gel*g
     #ad*gal)+A2e0c3p3t2c*GedSq*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0v_msv(whichlum,iset)
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
      CS0v_msv =
     #+ulum*(A1e0c1p1t1v*GeuSq*geu**2*gel**2+A1e0c2p1t1v*geu**2*GelSq*gel**2
     #+A1e0c3p1t1v*GeuGel*geu**2*gel**2+A1e0c4p1t1v*(gvusq*gvlsq*GeuSq+gausq
     #*gvlsq*GeuSq)+A1e0c5p1t1v*(gvusq*gvlsq*GelSq+gausq*gvlsq*GelSq
     #)+A1e0c6p1t1v*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel)+
     #A1e0c7p1t1v*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e0c8p1t1
     #v*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1v*(gvusq*
     #galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1v*GeuSq*geu*gel*gvu*
     #gvl+A1e0c11p1t1v*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1v*GeuGel*geu*gel*gvu
     #*gvl)+dlum*(A1e0c1p1t2v*GedSq*ged**2*gel**2+A1e0c2p1t2v*
     #ged**2*GelSq*gel**2+A1e0c3p1t2v*GedGel*ged**2*gel**2+A1e0c4p1t2v*(gvdsq*gvlsq*GedSq
     #+gadsq*gvlsq*GedSq)+A1e0c5p1t2v*(gvdsq*gvlsq*GelSq+ga
     #dsq*gvlsq*GelSq)+A1e0c6p1t2v*(gvdsq*gvlsq*GedGel+gadsq*gvlsq*GedGel
     #)+A1e0c7p1t2v*(gvdsq*galsq*GedSq+gadsq*galsq*GedSq
     #)+A1e0c8p1t2v*(gvdsq*galsq*GelSq+gadsq*galsq*GelSq)+A1e0c9p1
     #t2v*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel)+A1e0c10p1t2v*GedSq*ged*gel
     #*gvd*gvl+A1e0c11p1t2v*GelSq*ged*gel*gvd*gvl+A1e0c12p1t2v*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1v*gvugau*gvlgal*GeuSq+A2e0c2p1t1v*gvu
     #gau*gvlgal*GelSq+A2e0c3p1t1v*gvugau*gvlgal*GeuGel+A2e0c4p1t1
     #v*GeuSq*geu*gel*gau*gal+A2e0c5p1t1v*GelSq*geu*gel*gau*gal+A2e0c6p1t1v*GeuGel*geu*gel*g
     #au*gal)+dlum*(A2e0c1p1t2v*gvdgad*gvlgal*GedSq
     #+A2e0c2p1t2v*gvdgad*gvlgal*GelSq+A2e0c3p1t2v*gvdgad
     #*gvlgal*GedGel+A2e0c4p1t2v*GedSq*ged*gel*gad*gal+A2e0c5p1t2v*GelSq*ged*gel
     #*gad*gal+A2e0c6p1t2v*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n_msv(whichlum,iset)
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
      CS0n_msv =
     #+ulum*(A1e0c1p1t1n*GeuSq*geu**2*gel**2+A1e0c2p1t1n*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n*GeuGel*geu**2*gel**2+A1e0c4p1t1n*(gvusq*gvlsq*GeuSq+gausq
     #*gvlsq*GeuSq)+A1e0c5p1t1n*(gvusq*gvlsq*GelSq+gausq*gvlsq*GelSq
     #)+A1e0c6p1t1n*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel)+
     #A1e0c7p1t1n*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e0c8p1t1
     #n*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n*(gvusq*
     #galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n*GeuSq*geu*gel*gvu*
     #gvl+A1e0c11p1t1n*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n*GeuGel*geu*gel*gvu
     #*gvl)+dlum*(A1e0c1p1t2n*GedSq*ged**2*gel**2+A1e0c2p1t2n*
     #ged**2*GelSq*gel**2+A1e0c3p1t2n*GedGel*ged**2*gel**2+A1e0c4p1t2n*(gvdsq*gvlsq*GedSq
     #+gadsq*gvlsq*GedSq)+A1e0c5p1t2n*(gvdsq*gvlsq*GelSq+ga
     #dsq*gvlsq*GelSq)+A1e0c6p1t2n*(gvdsq*gvlsq*GedGel+gadsq*gvlsq*GedGel
     #)+A1e0c7p1t2n*(gvdsq*galsq*GedSq+gadsq*galsq*GedSq
     #)+A1e0c8p1t2n*(gvdsq*galsq*GelSq+gadsq*galsq*GelSq)+A1e0c9p1
     #t2n*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel)+A1e0c10p1t2n*GedSq*ged*gel
     #*gvd*gvl+A1e0c11p1t2n*GelSq*ged*gel*gvd*gvl+A1e0c12p1t2n*GedGel*ged*gel*g
     #vd*gvl)
     #+ulum*(A2e0c1p1t1n*gvugau*gvlgal*GeuSq+A2e0c2p1t1n*gvu
     #gau*gvlgal*GelSq+A2e0c3p1t1n*gvugau*gvlgal*GeuGel+A2e0c4p1t1
     #n*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n*GelSq*geu*gel*gau*gal+A2e0c6p1t1n*GeuGel*geu*gel*g
     #au*gal)+dlum*(A2e0c1p1t2n*gvdgad*gvlgal*GedSq
     #+A2e0c2p1t2n*gvdgad*gvlgal*GelSq+A2e0c3p1t2n*gvdgad
     #*gvlgal*GedGel+A2e0c4p1t2n*GedSq*ged*gel*gad*gal+A2e0c5p1t2n*GelSq*ged*gel
     #*gad*gal+A2e0c6p1t2n*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n1_msv(whichlum,iset)
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
      CS0n1_msv =
     #+ulum*(A1e0c1p1t1n1*GeuSq*geu**2*gel**2+A1e0c2p1t1n1*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n1*GeuGel*geu**2*gel**2+A1e0c4p1t1n1*(gvusq*gvlsq*GeuSq+
     #gausq*gvlsq*GeuSq)+A1e0c5p1t1n1*(gvusq*gvlsq*GelSq+gausq*gv
     #lsq*GelSq)+A1e0c6p1t1n1*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel
     #)+A1e0c7p1t1n1*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e
     #0c8p1t1n1*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n1
     #*(gvusq*galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n1*GeuSq*geu*gel
     #*gvu*gvl+A1e0c11p1t1n1*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n1*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2n1*GedSq*ged**2*gel**2
     #+A1e0c2p1t2n1*ged**2*GelSq*gel**2+A1e0c3p1t2n1*GedGel*ged**2*gel**2+A1e0c4p1t2n
     #1*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2n1*(gvdsq
     #*gvlsq*GelSq+gadsq*gvlsq*GelSq)+A1e0c6p1t2n1*(gvdsq*gvlsq*GedGel
     #+gadsq*gvlsq*GedGel)+A1e0c7p1t2n1*(gvdsq*galsq*GedSq+g
     #adsq*galsq*GedSq)+A1e0c8p1t2n1*(gvdsq*galsq*GelSq+gadsq*gal
     #sq*GelSq)+A1e0c9p1t2n1*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel
     #)+A1e0c10p1t2n1*GedSq*ged*gel*gvd*gvl+A1e0c11p1t2n1*GelSq*ged*gel*gvd
     #*gvl+A1e0c12p1t2n1*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n1*gvugau*gvlgal*GeuSq+A2e0c2p1t1n1*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n1*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n1*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n1*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n1*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n1*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n1*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #1*gvdgad*gvlgal*GedGel+A2e0c4p1t2n1*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n1*GelSq*ged*gel*gad*gal+A2e0c6p1t2n1*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n2_msv(whichlum,iset)
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
      CS0n2_msv =
     #+ulum*(A1e0c1p1t1n2*GeuSq*geu**2*gel**2+A1e0c2p1t1n2*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n2*GeuGel*geu**2*gel**2+A1e0c4p1t1n2*(gvusq*gvlsq*GeuSq+
     #gausq*gvlsq*GeuSq)+A1e0c5p1t1n2*(gvusq*gvlsq*GelSq+gausq*gv
     #lsq*GelSq)+A1e0c6p1t1n2*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel
     #)+A1e0c7p1t1n2*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e
     #0c8p1t1n2*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n2
     #*(gvusq*galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n2*GeuSq*geu*gel
     #*gvu*gvl+A1e0c11p1t1n2*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n2*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2n2*GedSq*ged**2*gel**2
     #+A1e0c2p1t2n2*ged**2*GelSq*gel**2+A1e0c3p1t2n2*GedGel*ged**2*gel**2+A1e0c4p1t2n
     #2*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2n2*(gvdsq
     #*gvlsq*GelSq+gadsq*gvlsq*GelSq)+A1e0c6p1t2n2*(gvdsq*gvlsq*GedGel
     #+gadsq*gvlsq*GedGel)+A1e0c7p1t2n2*(gvdsq*galsq*GedSq+g
     #adsq*galsq*GedSq)+A1e0c8p1t2n2*(gvdsq*galsq*GelSq+gadsq*gal
     #sq*GelSq)+A1e0c9p1t2n2*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel
     #)+A1e0c10p1t2n2*GedSq*ged*gel*gvd*gvl+A1e0c11p1t2n2*GelSq*ged*gel*gvd
     #*gvl+A1e0c12p1t2n2*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n2*gvugau*gvlgal*GeuSq+A2e0c2p1t1n2*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n2*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n2*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n2*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n2*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n2*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n2*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #2*gvdgad*gvlgal*GedGel+A2e0c4p1t2n2*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n2*GelSq*ged*gel*gad*gal+A2e0c6p1t2n2*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n3_msv(whichlum,iset)
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
      CS0n3_msv =
     #+ulum*(A1e0c1p1t1n3*GeuSq*geu**2*gel**2+A1e0c2p1t1n3*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n3*GeuGel*geu**2*gel**2+A1e0c4p1t1n3*(gvusq*gvlsq*GeuSq+
     #gausq*gvlsq*GeuSq)+A1e0c5p1t1n3*(gvusq*gvlsq*GelSq+gausq*gv
     #lsq*GelSq)+A1e0c6p1t1n3*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel
     #)+A1e0c7p1t1n3*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e
     #0c8p1t1n3*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n3
     #*(gvusq*galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n3*GeuSq*geu*gel
     #*gvu*gvl+A1e0c11p1t1n3*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n3*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2n3*GedSq*ged**2*gel**2
     #+A1e0c2p1t2n3*ged**2*GelSq*gel**2+A1e0c3p1t2n3*GedGel*ged**2*gel**2+A1e0c4p1t2n
     #3*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2n3*(gvdsq
     #*gvlsq*GelSq+gadsq*gvlsq*GelSq)+A1e0c6p1t2n3*(gvdsq*gvlsq*GedGel
     #+gadsq*gvlsq*GedGel)+A1e0c7p1t2n3*(gvdsq*galsq*GedSq+g
     #adsq*galsq*GedSq)+A1e0c8p1t2n3*(gvdsq*galsq*GelSq+gadsq*gal
     #sq*GelSq)+A1e0c9p1t2n3*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel
     #)+A1e0c10p1t2n3*GedSq*ged*gel*gvd*gvl+A1e0c11p1t2n3*GelSq*ged*gel*gvd
     #*gvl+A1e0c12p1t2n3*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n3*gvugau*gvlgal*GeuSq+A2e0c2p1t1n3*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n3*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n3*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n3*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n3*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n3*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n3*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #3*gvdgad*gvlgal*GedGel+A2e0c4p1t2n3*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n3*GelSq*ged*gel*gad*gal+A2e0c6p1t2n3*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n4_msv(whichlum,iset)
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
      CS0n4_msv =
     #+ulum*(A1e0c1p1t1n4*GeuSq*geu**2*gel**2+A1e0c2p1t1n4*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n4*GeuGel*geu**2*gel**2+A1e0c4p1t1n4*(gvusq*gvlsq*GeuSq+
     #gausq*gvlsq*GeuSq)+A1e0c5p1t1n4*(gvusq*gvlsq*GelSq+gausq*gv
     #lsq*GelSq)+A1e0c6p1t1n4*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel
     #)+A1e0c7p1t1n4*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e
     #0c8p1t1n4*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n4
     #*(gvusq*galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n4*GeuSq*geu*gel
     #*gvu*gvl+A1e0c11p1t1n4*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n4*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2n4*GedSq*ged**2*gel**2
     #+A1e0c2p1t2n4*ged**2*GelSq*gel**2+A1e0c3p1t2n4*GedGel*ged**2*gel**2+A1e0c4p1t2n
     #4*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2n4*(gvdsq
     #*gvlsq*GelSq+gadsq*gvlsq*GelSq)+A1e0c6p1t2n4*(gvdsq*gvlsq*GedGel
     #+gadsq*gvlsq*GedGel)+A1e0c7p1t2n4*(gvdsq*galsq*GedSq+g
     #adsq*galsq*GedSq)+A1e0c8p1t2n4*(gvdsq*galsq*GelSq+gadsq*gal
     #sq*GelSq)+A1e0c9p1t2n4*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel
     #)+A1e0c10p1t2n4*GedSq*ged*gel*gvd*gvl+A1e0c11p1t2n4*GelSq*ged*gel*gvd
     #*gvl+A1e0c12p1t2n4*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n4*gvugau*gvlgal*GeuSq+A2e0c2p1t1n4*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n4*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n4*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n4*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n4*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n4*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n4*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #4*gvdgad*gvlgal*GedGel+A2e0c4p1t2n4*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n4*GelSq*ged*gel*gad*gal+A2e0c6p1t2n4*GedGel*ged*gel*gad*gal)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CS0n5_msv(whichlum,iset)
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
      CS0n5_msv =
     #+ulum*(A1e0c1p1t1n5*GeuSq*geu**2*gel**2+A1e0c2p1t1n5*geu**2*GelSq*gel**2
     #+A1e0c3p1t1n5*GeuGel*geu**2*gel**2+A1e0c4p1t1n5*(gvusq*gvlsq*GeuSq+
     #gausq*gvlsq*GeuSq)+A1e0c5p1t1n5*(gvusq*gvlsq*GelSq+gausq*gv
     #lsq*GelSq)+A1e0c6p1t1n5*(gvusq*gvlsq*GeuGel+gausq*gvlsq*GeuGel
     #)+A1e0c7p1t1n5*(gvusq*galsq*GeuSq+gausq*galsq*GeuSq)+A1e
     #0c8p1t1n5*(gvusq*galsq*GelSq+gausq*galsq*GelSq)+A1e0c9p1t1n5
     #*(gvusq*galsq*GeuGel+gausq*galsq*GeuGel)+A1e0c10p1t1n5*GeuSq*geu*gel
     #*gvu*gvl+A1e0c11p1t1n5*GelSq*geu*gel*gvu*gvl+A1e0c12p1t1n5*GeuGel*geu*gel*g
     #vu*gvl)+dlum*(A1e0c1p1t2n5*GedSq*ged**2*gel**2
     #+A1e0c2p1t2n5*ged**2*GelSq*gel**2+A1e0c3p1t2n5*GedGel*ged**2*gel**2+A1e0c4p1t2n
     #5*(gvdsq*gvlsq*GedSq+gadsq*gvlsq*GedSq)+A1e0c5p1t2n5*(gvdsq
     #*gvlsq*GelSq+gadsq*gvlsq*GelSq)+A1e0c6p1t2n5*(gvdsq*gvlsq*GedGel
     #+gadsq*gvlsq*GedGel)+A1e0c7p1t2n5*(gvdsq*galsq*GedSq+g
     #adsq*galsq*GedSq)+A1e0c8p1t2n5*(gvdsq*galsq*GelSq+gadsq*gal
     #sq*GelSq)+A1e0c9p1t2n5*(gvdsq*galsq*GedGel+gadsq*galsq*GedGel
     #)+A1e0c10p1t2n5*GedSq*ged*gel*gvd*gvl+A1e0c11p1t2n5*GelSq*ged*gel*gvd
     #*gvl+A1e0c12p1t2n5*GedGel*ged*gel*gvd*gvl)
     #+ulum*(A2e0c1p1t1n5*gvugau*gvlgal*GeuSq+A2e0c2p1t1n5*g
     #vugau*gvlgal*GelSq+A2e0c3p1t1n5*gvugau*gvlgal*GeuGel+A2e0c4p
     #1t1n5*GeuSq*geu*gel*gau*gal+A2e0c5p1t1n5*GelSq*geu*gel*gau*gal+A2e0c6p1t
     #1n5*GeuGel*geu*gel*gau*gal)+dlum*(A2e0c1p1t2n5*gvdga
     #d*gvlgal*GedSq+A2e0c2p1t2n5*gvdgad*gvlgal*GelSq+A2e0c3p1t2n
     #5*gvdgad*gvlgal*GedGel+A2e0c4p1t2n5*GedSq*ged*gel*gad*gal+A2e0c5p1
     #t2n5*GelSq*ged*gel*gad*gal+A2e0c6p1t2n5*GedGel*ged*gel*gad*gal)
      RETURN
      END
