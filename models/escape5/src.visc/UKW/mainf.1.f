
c----------------------------------------------------------------
c              UrQMD adapted by Klaus WERNER  April 2010 
c----------------------------------------------------------------

      integer function geteposcode(n)
      geteposcode = idtrafo('pdg','nxs',n)
      end

      subroutine urqmdmain
      include 'aaa.h'
      laproj = 79 
      maproj = 197 
      latarg = 79 
      matarg = 197
      engy = 200    !GeV
      nptlpt=0   ! update (Klaus) 23.04.2010
c      nptlbd=390   ! update (Klaus) 23.04.2010
      call uKW
      end
      
c-------------- get initial particles ---------------------------------
      subroutine getInitialParticles(nn)
      include 'aaa.h'
      integer ninpart

      call cxxninit(ninpart) ;
      do nn=1,ninpart
      istptl(nn)=0
      ityptl(nn)=0   ! type: target, etc
      call cxxinit(nn,idptl(nn),xorptl(1,nn),xorptl(2,nn),
     . xorptl(3,nn),xorptl(4,nn),pptl(1,nn),pptl(2,nn),pptl(3,nn),
     . pptl(4,nn),pptl(5,nn))
      idptl(nn)=idtrafo('pdg','nxs', idptl(nn))
      end do
      nn=ninpart

      print*,'(info) ',nn,' particles read from file:'
c      do i=1,nn
c      print*, idptl(i),xorptl(1,i),xorptl(2,i),
c     . xorptl(3,i),xorptl(4,i),pptl(1,i),pptl(2,i),pptl(3,i),
c     . pptl(4,i),pptl(5,i)
c      end do
      nptlbd=nn   ! update (Klaus) 23.04.2010
      end

c----------------------- store initial particles to file -------------
      subroutine storepart(nn)
      include 'aaa.h'
      call cxxninit(ninpart) ;
      do nn=1,ninpart
      istptl(nn)=0
      ityptl(nn)=0   ! type: target, etc
      call cxxinit(nn,idptl(nn),xorptl(1,nn),xorptl(2,nn),
     . xorptl(3,nn),xorptl(4,nn),pptl(1,nn),pptl(2,nn),pptl(3,nn),
     . pptl(4,nn),pptl(5,nn))
      idptl(nn)=idtrafo('pdg','nxs', idptl(nn))
      end do
      nn=ninpart

      print*,'(info) ',nn,' particles written:'
      open (unit=20,file="uuu.data",action="write",status="replace")
      do i=1,nn
      write(20,*) istptl(i),ityptl(i),idptl(i),xorptl(1,i),xorptl(2,i),
     . xorptl(3,i),xorptl(4,i),pptl(1,i),pptl(2,i),pptl(3,i),
     . pptl(4,i),pptl(5,i)
      end do
      end

c---------------------- treat final particles -------------------------
      subroutine treatFinalParticles(nptl1,nptl2)
      include 'aaa.h'
      call cxxnfinal(nptl2-nptl1+1)
      print*, '---------final particles-------------'
      do i=nptl1,nptl2
      idptl(i)=idtrafo('nxs','pdg', idptl(i))
c      print*, idptl(i),xorptl(1,i),xorptl(2,i),
c     . xorptl(3,i),xorptl(4,i),pptl(1,i),pptl(2,i),pptl(3,i),
c     . pptl(4,i),pptl(5,i)
      call cxxfinal(i-nptl1+1,idptl(i),xorptl(1,i),xorptl(2,i),
     . xorptl(3,i),xorptl(4,i),pptl(1,i),pptl(2,i),pptl(3,i),
     . pptl(4,i),pptl(5,i))
      enddo
      end
      
      
      
