CCC part of pdfcontrol.f
      integer whichas, pdfCL, NPDFerror
      double precision fitscale
      COMMON/pdfcontrol/ fitscale, whichas, pdfCL, NPDFerror
      integer NSMAX
      parameter(NSMAX=204)

      double precision uubplum(3,0:NSMAX)
      double precision ddbplum(3,0:NSMAX)
      double precision uubmlum(3,0:NSMAX)
      double precision ddbmlum(3,0:NSMAX)
      double precision bbbplum(1,0:NSMAX)
      double precision bbbmlum(1,0:NSMAX)
      double precision phphplum(1,0:NSMAX)
      double precision phphmlum(1,0:NSMAX)
      common/partlum/ uubplum,ddbplum,uubmlum,ddbmlum,
     &                bbbplum,bbbmlum,
     &                phphplum,phphmlum

