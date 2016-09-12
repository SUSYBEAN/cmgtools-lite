CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC Running Coupling Constant CCCCCCCCCCCCCCCCCCCCCCCC

      double precision function SSC(x,asZ,modeas)
      implicit real*8(a-h,o-z)
      double precision asZ
      real*8 Mz
      Integer modeas,i

CCC   asZ should be consisten with PDF fit value

      if (modeas.eq.0) then
        
CCC   LO alpha_s/Pi

      b0= 1.916666667

      as0=asZ/3.1415926D0
      Mz=91.187D0
      step=dlog(x**2/Mz**2)/100
                                                                                                                                            
      do i=1,100
         as=as0 - step*b0*as0**2*(1.D0)
         as0 = as
      enddo
      SSC=as

      elseif (modeas.eq.1) then
 
CCC   NLO alpha_s/Pi
 
      b0= 1.916666667
      b1= 1.260869565
                                                                                                                                            
      as0=asZ/3.1415926D0
      Mz=91.187D0
      step=dlog(x**2/Mz**2)/100
                                                                                                                                            
      do i=1,100
         as=as0 - step*b0*as0**2*(1.D0+b1*as0)
         as0 = as
      enddo
      SSC=as
      
      elseif (modeas.eq.2) then

CC    NNLO alphas/pi                                                                                                                                            
      b0= 1.916666667
      b1=1.260869565
      b2=1.474788647
                                                            
      as0=asZ/3.1415926D0
      Mz=91.187D0
      step=dlog(x**2/Mz**2)/100D0
 
      do i=1,100
 
         as=as0 - step*b0*as0**2*(1.D0+b1*as0 + b2*as0**2 )
         as0 = as
      enddo
 
      SSC=as
 
      else
 
      SSC=0
 
      endif
 
      RETURN
      end

