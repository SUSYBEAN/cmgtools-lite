      subroutine fitpdf(fscale,iset,Np)
      implicit none

      integer Np,iset
      double precision fscale
      integer i,i1,i2,i3
      
      include 'pdfcontrol.f'
      
      double precision xval,pdfval(1:10)
      double precision xx,yphot,yglu,yupv,ydnv,yusea,ydsea,ysv,yssea,ychm,
     &   ybot
      double precision hh,xa1,xa2,xa3,ya1,ya2,ya3,a,b,c

      dimension xx(Np) 
      dimension yphot(Np), yglu(Np), yupv(Np), ydnv(Np)
      dimension yusea(Np), ydsea(Np), ysv(Np), yssea(Np)
      dimension ychm(Np), ybot(Np)

CCCC    ================================================================ 

      do i=1,Np
         xx(i) = tau + (1-tau)*i*(i-1d0)/Np/(Np-1)
         
         xval = xx(i)
         call ObtainPDFs(xval,fscale,iset,pdfval,0) !call PDF set

         yupv(i) = pdfval(1)
         ydnv(i) = pdfval(2)
         yusea(i) = pdfval(3)
         ydsea(i) = pdfval(4)
         ysv(i) = pdfval(5)
         yssea(i) = pdfval(6)
         ychm(i) = pdfval(7)
         ybot(i) = pdfval(8)
         yglu(i) = pdfval(9)
         yphot(i) = pdfval(10)
      enddo
      
      hh = Np*(Np-1)/(1-tau)
      
      hhph = hh
      hhg = hh
      hhuv = hh
      hhdv = hh
      hhus = hh
      hhds = hh
      hhsv = hh
      hhss = hh
      hhchm = hh
      hhbtm = hh

      do i=1,Np 
         if (i==1) then
            i1=i
            i2=i+1
            i3=i+2
         elseif (i==Np) then 
            i1=i-2
            i2=i-1
            i3=i
         else
            i1=i-1
            i2=i
            i3=i+1
         endif 
         
         xa1 = xx(i1)
         xa2 = xx(i2)
         xa3 = xx(i3)
         
         ya1 =yphot(i1)
         ya2 =yphot(i2)
         ya3 =yphot(i3)
        
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
        
         xxph(iset,i) = xx(i)
         AAPH(iset,i) = a
         BBPH(iset,i) = b
         CCPH(iset,i) = c

         ya1 =yglu(i1)
         ya2 =yglu(i2)
         ya3 =yglu(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxg(iset,i) = xx(i)
         AAG(iset,i) = a
         BBG(iset,i) = b
         CCG(iset,i) = c
   
         ya1 =yupv(i1)
         ya2 =yupv(i2)
         ya3 =yupv(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxuv(iset,i) = xx(i)
         AAUV(iset,i) = a
         BBUV(iset,i) = b
         CCUV(iset,i) = c
         
         ya1 =ydnv(i1)
         ya2 =ydnv(i2)
         ya3 =ydnv(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxdv(iset,i) = xx(i)
         AADV(iset,i) = a
         BBDV(iset,i) = b
         CCDV(iset,i) = c
         
         ya1 =yusea(i1)
         ya2 =yusea(i2)
         ya3 =yusea(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxus(iset,i) = xx(i)
         AAUS(iset,i) = a
         BBUS(iset,i) = b
         CCUS(iset,i) = c

         ya1 =ydsea(i1)
         ya2 =ydsea(i2)
         ya3 =ydsea(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxds(iset,i) = xx(i)
         AADS(iset,i) = a
         BBDS(iset,i) = b
         CCDS(iset,i) = c
         
         ya1 =ysv(i1)
         ya2 =ysv(i2)
         ya3 =ysv(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxsv(iset,i) = xx(i)
         AASV(iset,i) = a
         BBSV(iset,i) = b
         CCSV(iset,i) = c

         ya1 =yssea(i1)
         ya2 =yssea(i2)
         ya3 =yssea(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxss(iset,i) = xx(i)
         AASS(iset,i) = a
         BBSS(iset,i) = b
         CCSS(iset,i) = c
        
         ya1 =ychm(i1)
         ya2 =ychm(i2)
         ya3 =ychm(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxchm(iset,i) = xx(i)
         AACHM(iset,i) = a
         BBCHM(iset,i) = b
         CCCHM(iset,i) = c
         
         ya1 =ybot(i1)
         ya2 =ybot(i2)
         ya3 =ybot(i3)
         
         call givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
         
         xxbtm(iset,i) = xx(i)
         AABTM(iset,i) = a
         BBBTM(iset,i) = b
         CCBTM(iset,i) = c
         
      enddo

      return
      end
      
      
      subroutine givecoeff(xa1,xa2,xa3,ya1,ya2,ya3, a, b, c)
      implicit none
      double precision xa1,xa2,xa3,ya1,ya2,ya3,a,b,c
      double precision den

      den=(xa1-xa2)*(xa2-xa3)*(xa3-xa1)
      a = -ya1*(xa2-xa3)/den
      a = a-ya2*(xa3-xa1)/den
      a = a-ya3*(xa1-xa2)/den
      
      b =  ya1*(xa2*xa2-xa3*xa3)/den
      b = b+ya2*(xa3*xa3-xa1*xa1)/den
      b = b+ya3*(xa1*xa1-xa2*xa2)/den
      
      c = -ya1*(xa2-xa3)*xa2*xa3/den
      c = c-ya2*(xa3-xa1)*xa1*xa3/den
      c = c-ya3*(xa1-xa2)*xa1*xa2/den
      return
      end
