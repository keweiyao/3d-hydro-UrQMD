
c----------------------------------------------------------------
c              UrQMD adapted by Klaus WERNER  April 2010 
c----------------------------------------------------------------

      subroutine uKW
      include 'aaa.h'
      character*500 edir
      common/cdir/edir
      call idresi
      edir='./ '
      ihacas=3
      call uinitial
      call uepos
      call uexit
      end
      subroutine uinitial
      implicit none
c     debug and validity range
      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 40000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn
      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer
      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
      real*8 
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8
     
      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      real*8 spx(nspl), spPauy(nspl), outPau(nspl), 
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)
      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww
      real*8 
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax), 
     +     dectime(nmax), tform(nmax), xtotfac(nmax)
      
      
      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid
     
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz
      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep
c spectator arrays
	integer smax
	parameter(smax=500)  ! maximum number of spectators
	real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +	       p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +	       sfmass(smax)
	
        integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)
	integer nspec
	common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass
	common /sisys/ sspin, scharge, sityp, siso3, suid
	common /ssys/ nspec
        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)
        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
	integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
	integer   CTOption(numcto)
	real*8    CTParam(numctp)
        common /options/CTOption,CTParam
	real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz 
	      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes
      parameter (minnuc=1) ! lowest baryon particle ID 
      parameter (minmes=100) ! lowest meson particle ID
      parameter (maxmes=132) ! hightest meson particle ID
c number of resonances of a kind
      parameter (numnuc=16) ! number of nucleon resonances
      parameter (numdel=10) ! number of delta resonances
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1) ! highest nucleon ID
      parameter (mindel=minnuc+maxnuc)   ! lowest delta ID
      parameter (maxdel=mindel+numdel-1) ! highest delta ID
c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1) ! lowest baryon resonance ID
      parameter (maxres=maxdel)   ! highest (nonstrange) baryon 
                                  ! resonance ID
c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13) ! number of lambda states
      parameter (numsig=9)  ! number of sigma states
      parameter (numcas=6)  ! number of cascade states
      parameter (numome=1)  ! number of omega states
      parameter (minlam=mindel+numdel)   ! ID of lowest lambda state
      parameter (maxlam=minlam+numlam-1) ! ID of highest lambda state
      parameter (minsig=minlam+numlam)   ! ID of lowest sigma state
      parameter (maxsig=minsig+numsig-1) ! ID of highest sigma state
      parameter (mincas=minsig+numsig)   ! ID of lowest cascade state
      parameter (maxcas=mincas+numcas-1) ! ID of highest cascade state
      parameter (minome=mincas+numcas)   ! ID of lowest omega state
      parameter (maxome=minome+numome-1) ! ID of highest omega state
c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc) ! ID of lowest baryon state
      parameter (maxbar=maxome) ! ID of highest baryon state
      parameter (offmeson=minmes) ! offset between zero and lowest 
                                  ! meson state
      parameter (maxmeson=maxmes) ! ID of highest meson state
c... these variables are in principal obsolete and should be exchanged 
c           were referenced 
c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)   ! ID of kaon
      parameter (itrho=104)    ! ID of rho meson 
      parameter (itome=103)    ! ID of omega meson
      parameter (iteta=102)    ! ID of eta
      parameter (itphi=109)    ! ID of phi
      parameter (itetapr=107)  ! ID of eta'
      parameter (pimeson=101)  ! ID of $\pi$
      parameter (nucleon=minnuc) ! ID of nucleon
      integer itmin,itmax
      parameter (itmin=minnuc)  ! lowest defined ID
      parameter (itmax=maxmes)  ! highest defined ID
      parameter (maxbra=11)  ! decay channels for $s=0$ baryon resonances
      parameter (maxbrm=25) ! decay channels for meson resonances
      parameter (maxbrs1=10)! decay channels for $s=1$ baryon resonances
      parameter (maxbrs2=3) ! decay channels for $s=2$ baryon resonances
       integer mlt2it(maxmes-minmes) ! meson IDs sorted by multipletts
      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)      ! offset for mass generation
      parameter (mresmin=1.0765d0)  ! minimum baryon resonance mass
      parameter (mresmax=5d0)       ! maximum baryon resonance mass
      character*45 versiontag
      common /versioning/ versiontag
      real*8 massres(minbar:maxbar),widres(minbar:maxbar)
      real*8 branmes(0:maxbrm,minmes+1:maxmes)
      real*8 branres(0:maxbra,minnuc+1:maxdel)
      real*8 branbs1(0:maxbrs1,minlam+1:maxsig)
      real*8 branbs2(0:maxbrs2,mincas+1:maxcas)
      integer Jres(minbar:maxbar)
      integer Jmes(minmes:maxmes)
      integer pares(minbar:maxbar),pames(minmes:maxmes)
      integer Isores(minbar:maxbar), Isomes(minmes:maxmes)
      integer brtype(4,0:maxbra),bmtype(4,0:maxbrm)
      integer bs1type(4,0:maxbrs1),bs2type(4,0:maxbrs2)
      real*8 massmes(minmes:maxmes)
      real*8 mmesmn(minmes:maxmes)
      real*8 widmes(minmes:maxmes)
      integer strres(minbar:maxbar),strmes(minmes:maxmes)
      integer lbr(0:maxbra,minnuc+1:maxdel)
      integer lbs1(0:maxbrs1,minlam+1:maxsig)
      integer lbs2(0:maxbrs2,mincas+1:maxcas)
      integer lbm(0:maxbrm,minmes+1:maxmes)
      common /resonances/ massres,widres,massmes,widmes,mmesmn,
     ,                    branres,branmes,branbs1,branbs2,
     ,                    bs1type,bs2type,lbs1,lbs2,lbm,
     ,                    jres,jmes,lbr,brtype,pares,pames,
     ,                    bmtype,
     ,                    Isores,Isomes,strres,strmes,mlt2it
      !----------------------------------------------------------------------------     
      !epos common blocks for particle list
      !----------------------------------------------------------------------------     
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer iappl,model
      common/appli/iappl,model
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list
      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)
      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
      integer nptl
      integer nevt,kolevt,koievt,npjevt,ntgevt,npnevt,nppevt,
     *        ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt,
     *        nglevt,minfra,maxfra
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt   
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
      integer laproj,maproj,latarg,matarg,nptlbd,nptlpt
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/c4ptl/nptlbd,nptlpt     
      integer istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      integer ihacas
      common/chacas/ihacas
      !---------------------------------------------------------------         
      integer i,nn,idtmp,ityptmp,iso3tmp,itmp
      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      real*8 dectim
      external dectim
      real*8 mintime,eb
      integer j,k,icount,npold
      integer strcount
      common /inewpartxx/ strcount
      real*8 etot
      integer igeteost,igethyt,ireadhyt
      common/hydr2/igeteost,igethyt,ireadhyt
      integer nrap,ntau,nrad,id
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      integer nupa(0:10)
      common/cnupa/nupa
      real*8 deluutau
      integer nuutaumx
      common /cf15/deluutau,nuutaumx
      integer nui,nw,nwi
      data nui/0/
      save nui
      !integer nptest,ia,jca(3),nflav
      !parameter (nflav=6)
      !integer jc(nflav,2),ic(2)
      integer itypart(nmax),iorpart(nmax)
      common /city/itypart /corpart/iorpart

c##############################     subroutine uinitial

      nui=nui+1
      if(nui.eq.1)then
        write(ifmt,'(a,50a1)')' (info)',('u',i=1,50)
        do j=-1,40
          do i=-10,10
            dcoll(i,j)=0d0
          enddo
          do i=-10,10
            dpart(i,j)=0d0
          enddo
          do i=0,20
            rcoll(i,j)=0
          enddo
        enddo
        zevents=0  
        do i=0,10
        nupa(i)=0
        enddo
      endif
      zevents=zevents+1
      
      do nn=1,nmax
        itypart(nn)=0
        iorpart(nn)=0
      enddo
       
      call uinit(0)
      call osc_header
      call osc99_header

      npart = 0
      npold = 0
      nbar=0
      nmes=0
      uid_cnt=0
c reset counters
c     all collisions/decays
      ctag  = 0
c     all decays
      dectag = 0
c     number of prod. hard resonances
      NHardRes=0
c     number of prod. soft resonances
      NSoftRes=0
c     number of prod. resonances via decay
      NDecRes=0
c     number of blocked collisions
      NBlColl=0
c     number of elastic collisions
      NElColl=0
c     number of strings
      strcount=1
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
	spin(j)  = 0
	ncoll(j) = 0
	lstcoll(j)=0
	r0(j) = 0.0
	rx(j)	 = 0.0
	ry(j)	 = 0.0
	rz(j)	 = 0.0
	p0(j)	 = 0.0
	px(j)	 = 0.0
	py(j)	 = 0.0
	pz(j)	 = 0.0
	frr0(j) = 0.0
	frrx(j)    = 0.0
	frry(j)    = 0.0
	frrz(j)    = 0.0
	frp0(j)    = 0.0
	frpx(j)    = 0.0
	frpy(j)    = 0.0
	frpz(j)    = 0.0
	fmass(j) = 0.0
	charge(j)= 0
	iso3(j)  = 0
	ityp(j)  = 0
	dectime(j)= 0.0
	origin(j)=0
	tform(j)=0.0
	xtotfac(j)=1.0
	strid(j)=0
	uid(j)=0
	 ffermpx(j) = 0.0
	 ffermpy(j) = 0.0
	 ffermpz(j) = 0.0

	 do 21 k=1,2
	    p0td(k,j)=0.d0
	    pxtd(k,j)=0.d0
	    pytd(k,j)=0.d0
	    pztd(k,j)=0.d0
	    fmasstd(k,j)=0.d0
	    ityptd(k,j)=0
	   iso3td(k,j)=0
 21	 continue
 20   continue


c epos event info to u event info      
      bimp = bimevt
      
c initialise
      npart=0  !number of particles transferred to u
      mintime = 1d2 !the minimum formation time

c energy of spectators
      etot=0 
      do nn=1,nptlpt      
        if(istptl(nn).eq.0)etot=etot+pptl(4,nn)
      enddo

c ptls from testFile

      !-----------------------------------
      ! 1st run with optns:  set ihacas 2 
      ! 2nd run with optns:  set ihacas 3
      !-----------------------------------
      if(ihacas.eq.3)then
        nptlpt=0
        call getInitialParticles(nn)
        nptlbd=nn
      endif  

c decay unknown particles

      nw=nptlbd
      nwi=nw+1
      do nn=nptlpt+1,nptlbd       
        if(istptl(nn).eq.0)then
          call idcorrect(ityptl(nn),idptl(nn),idtmp)
          call pdg2id(ityptmp,iso3tmp,idtmp)
          if(abs(ityptmp).gt.1000)then  !unkown hadrons
            nw=nw+1
            !write(ifch,*)' ~~~~~~~ unknown hadron',idtmp
            !.      ,ityptmp,iso3tmp,nn,nw
            istptl(nw)   = istptl(nn)
            idptl(nw)    = idptl(nn)
            ityptl(nw)   = ityptl(nn)
            do i=1,4
            xorptl(i,nw) = xorptl(i,nn)
            enddo
            do i=1,5
            pptl(i,nw)   = pptl(i,nn)
            enddo
            ibptl(1,nw)  = ibptl(1,nn)
            iorptl(nw)   = nn
            do i=1,2
            tivptl(i,nw)  = tivptl(i,nn)
            enddo
            istptl(nn)=3
            ifrptl(1,nn)=nw
            ifrptl(2,nn)=nw
          endif
        endif
      enddo       
      if(nw.ge.nwi)then
        nptl=nw
        call decayall(nwi,999)
        nptlbd=nptl
        !call alist(' ~~~~~~~ list from unknown hadrons&',nwi,nptl)
        do nn=nwi,nptl  
          if(istptl(nn).eq.0)then
            call idcorrect(ityptl(nn),idptl(nn),idtmp)
          endif
          if(istptl(nn).eq.1)istptl(nn)=4
        enddo     
      endif

c fill in the baryons first     
      !write(ifmt,'(a/)')' '
      !nptest=0
      !do i=1,3
      !jca(i)=0
      !enddo
      if(ihacas.eq.2)call openNewTestFile
      if(ish.ge.2)call alist('(u) initial baryons&',0,0)
      if(ish.ge.2)write(ifch,*)'check between ',nptlpt+1,' and ',nptlbd 
      if(ish.ge.2)write(ifch,*)' '
      nbar = 0
      do nn=nptlpt+1,nptlbd       
        if(istptl(nn).eq.0)then
        if(ihacas.eq.2)call writeParticleTestFile(nn)
           if(ityptl(nn).eq.61)then
           idtmp=idptl(nn)
           else
           idtmp=idtrafo('nxs','pdg',idptl(nn))
           endif
           call pdg2id(ityptmp,iso3tmp,idtmp)
           if(abs(ityptmp).gt.1000)then
              write(ifch,*)' ~~~~~~~ ERROR: unkown hadron'
     .         ,ityptl(nn),idtmp,ityptmp,iso3tmp,nn
           endif
           if(abs(ityptmp).le.maxbar) then
               if(ish.ge.2)call alist('&',nn,nn)
               !ia=abs(idptl(nn))
               !if(ia.eq.2130.or.ia.eq.1230.or.ia.eq.1130.or.ia.eq.2230
               !.         .or.ia.eq.1131.or.ia.eq.1231.or.ia.eq.2231)
               !.         !nptest=nptest+1
               !call idtr4(idptl(nn),ic)
               !call iddeco(ic,jc)
               !do i=1,3
               !jca(i)=jca(i)+jc(i,1)-jc(i,2)
               !enddo
               !if(jc(3,1)-jc(3,2).ne.0)write(ifmt,'(i6,$)')idptl(nn)
               etot=etot+pptl(4,nn)
               nbar=nbar+1
               call DelayFormation(nn)
               r0(nbar)=xorptl(4,nn)
               rx(nbar)=xorptl(1,nn)
               ry(nbar)=xorptl(2,nn)
               rz(nbar)=xorptl(3,nn)
               p0(nbar)=pptl(4,nn)
               px(nbar)=pptl(1,nn)
               py(nbar)=pptl(2,nn)
               pz(nbar)=pptl(3,nn)
               fmass(nbar)=pptl(5,nn)
               ityp(nbar)=ityptmp
               iso3(nbar)=iso3tmp
               if(ityptl(nn).eq.61)then
               charge(nbar)=ibptl(1,nn)
               else
               charge(nbar)=fchg(iso3(nbar),ityp(nbar))
               endif
               itypart(nbar)=ityptl(nn)
               iorpart(nbar)=0
               lstcoll(nbar)=0
               ncoll(nbar)=0
               origin(nbar)=0
               tform(nbar)=r0(nbar)
               dectime(nbar)=dectim(nbar,1)+tform(nbar)
               xtotfac(nbar)=0d0
               if(r0(nbar).lt.mintime) mintime = r0(nbar)
           endif
        endif
      enddo
      
c then fill in the mesons
      if(ish.ge.2)call alist('(u) initial mesons&',0,0)
      if(ish.ge.2)write(ifch,*)'check between ',nptlpt+1,' and ',nptlbd 
      if(ish.ge.2)write(ifch,*)' '
      nmes = 0
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.0)then
           if(ityptl(nn).eq.61)then
           idtmp=idptl(nn)
           else
           idtmp=idtrafo('nxs','pdg',idptl(nn))
           endif
           call pdg2id(ityptmp,iso3tmp,idtmp)
           !if(abs(ityptmp).gt.1000)then
           !   write(ifch,*)' ~~~~~~~ ERROR: unkown hadron'
           !.       ,ityptl(nn),idtmp,ityptmp,iso3tmp,nn
           !endif
           if(abs(ityptmp).ge.minmes) then
               !if()then
               !print*,'+++++ ids: ',nn,idptl(nn),idtmp,ityptmp
               !endif
               if(ish.ge.2)call alist('&',nn,nn)
               nmes=nmes+1
               !call idtr4(idptl(nn),ic)
               !call iddeco(ic,jc)
               !do i=1,3
               !jca(i)=jca(i)+jc(i,1)-jc(i,2)
               !enddo
               !if(jc(3,1)-jc(3,2).ne.0)write(ifmt,'(i6,$)')idptl(nn)
               etot=etot+pptl(4,nn)
               itmp=nbar+nmes
               call DelayFormation(nn)
               r0(itmp)=xorptl(4,nn)
               rx(itmp)=xorptl(1,nn)
               ry(itmp)=xorptl(2,nn)
               rz(itmp)=xorptl(3,nn)
               p0(itmp)=pptl(4,nn)
               px(itmp)=pptl(1,nn)
               py(itmp)=pptl(2,nn)
               pz(itmp)=pptl(3,nn)
               fmass(itmp)=pptl(5,nn)
               ityp(itmp)=ityptmp
               iso3(itmp)=iso3tmp
               if(ityptl(nn).eq.61)then
               charge(itmp)=ibptl(1,nn)
               else
               charge(itmp)=fchg(iso3(itmp),ityp(itmp))
               endif
               itypart(itmp)=ityptl(nn)
               iorpart(itmp)=0
               lstcoll(itmp)=0
               ncoll(itmp)=0
               origin(itmp)=0
               tform(itmp)=r0(itmp)
               dectime(itmp)=dectim(itmp,1)+tform(itmp)
               xtotfac(itmp)=0d0
               if(r0(itmp).lt.mintime) mintime = r0(itmp)
           endif
        endif
      enddo

      if(ihacas.eq.2)call writeEndEventTestFile

      !print*,'+++++ etot/eini=',etot/(197*200.)

      !write(ifmt,'(a,i5,3x,$)')'+++++ Nptest:',nptest
      !write(ifmt,'(a,3i5,3x,$)')'+++++ Flavor:',jca

      npart = nbar + nmes
      if(ish.ge.2)write(ifch,*)' '
      if(ish.ge.2)write(ifch,*)'nr of ptls: ',npart

      if(npart.eq.0)then
        write(ifmt,*)'********** uinitial: no particles'
        write(ifmt,*)'check between ',nptlpt+1,' and ',nptlbd 
        do nn=nptlpt+1,nptlbd     
           if(ityptl(nn).eq.61)then
           idtmp=idptl(nn)
           else
           idtmp=idtrafo('nxs','pdg',idptl(nn))
           endif
           call pdg2id(ityptmp,iso3tmp,idtmp)
           write(ifmt,*)abs(ityptmp).ge.minmes,
     .     idptl(nn),istptl(nn) ,ityptl(nn)   
        enddo  
        stop'091130'
      endif 

c back to the same starting time
      do i = 1, npart
         call getind(p0(i),pz(i),r0(i),rx(i),ry(i),rz(i)
     .    ,nrap,ntau,nrad,ityp(i),iso3(i),id)
          !if(ish.ge.2)write(ifch,'(4x,i6,i7,4x,3i3,3x,2f7.2,3x,2f7.2)')
          !.    i,id,ntau,nrad,nrap,sqrt(r0(i)**2-rz(i)**2)
          !.      ,sqrt(rx(i)**2+ry(i)**2)
         dpart(nrap,ntau)=dpart(nrap,ntau)+1
         dpart(nrap,nuutaumx+1)=dpart(nrap,nuutaumx+1)+1
         !save freeze-out configuration, in case of no further
         !rescatterings
         frr0(i) = r0(i)
         frrx(i) = rx(i)
         frry(i) = ry(i)
         frrz(i) = rz(i)
         frp0(i) = p0(i)
         frpx(i) = px(i)
         frpy(i) = py(i)
         frpz(i) = pz(i)
         rx(i)=rx(i)-px(i)/p0(i)*(r0(i)-mintime)
         ry(i)=ry(i)-py(i)/p0(i)*(r0(i)-mintime)
         rz(i)=rz(i)-pz(i)/p0(i)*(r0(i)-mintime)
         r0(i)=mintime
      enddo

      acttime=mintime

c      write(*,*)'DEBUG INFO (epos.f): ',mintime,npart,istmax,nbar,nmes

c keep  epos ptl list  
      nptl=nptlpt
      !the following for analysis, may be commented
      if(ish.ge.2)call alist('(u) kept EPOS list&',0,0)
      do nn=nptlpt+1,nptlbd  
        if(istptl(nn).eq.0)istptl(nn)=3
        nptl=nptl+1
        if(ish.ge.2)call alist('&',nn,nn)
      enddo
      
      if(nui.eq.1)write(ifmt,'(a,50a1)')' (info)',('u',i=1,50)
      return
      end

      subroutine DelayFormation(j)
      include 'aaa.h'
      if(tauhac.le.0.)return
      if(ityptl(j).ne.60)return
      r=rangen()
      tauran=-tauhac*alog(r)*pptl(4,j)/pptl(5,j)
      do j=1,4
        xorptl(j,j)=xorptl(j,j)
     &   + pptl(j,j) / pptl(4,j) * tauran
      enddo
      end   
      subroutine idcorrect(ity,id,idtmp)
      if(ity.eq.61)then
        if(id.eq.310)id=311
        if(id.eq.130)id=-311
        idtmp=id
      else
        if(id.eq.20)id=230
        if(id.eq.-20)id=-230
        idtmp=idtrafo('nxs','pdg',id)
      endif
      end !~~~~~~~~~~~~~~~~



                subroutine uexit      !transfer -> EPOS


      implicit none
c     debug and validity range
      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 40000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn
      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer
      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
      real*8 
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8
      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      real*8 spx(nspl), spPauy(nspl), outPau(nspl), 
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)
      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww
      real*8 
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax), 
     +     dectime(nmax), tform(nmax), xtotfac(nmax)
      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz
      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep
c spectator arrays
	integer smax
	parameter(smax=500)  ! maximum number of spectators
	real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +	       p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +	       sfmass(smax)
       integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)
	integer nspec
	common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass
	common /sisys/ sspin, scharge, sityp, siso3, suid
	common /ssys/ nspec
        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)
        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
	integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
	integer   CTOption(numcto)
	real*8    CTParam(numctp)
        common /options/CTOption,CTParam
	real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz 
	      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes
      parameter (minnuc=1) ! lowest baryon particle ID 
      parameter (minmes=100) ! lowest meson particle ID
      parameter (maxmes=132) ! hightest meson particle ID
c number of resonances of a kind
      parameter (numnuc=16) ! number of nucleon resonances
      parameter (numdel=10) ! number of delta resonances
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1) ! highest nucleon ID
      parameter (mindel=minnuc+maxnuc)   ! lowest delta ID
      parameter (maxdel=mindel+numdel-1) ! highest delta ID
c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1) ! lowest baryon resonance ID
      parameter (maxres=maxdel)   ! highest (nonstrange) baryon 
                                  ! resonance ID
c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13) ! number of lambda states
      parameter (numsig=9)  ! number of sigma states
      parameter (numcas=6)  ! number of cascade states
      parameter (numome=1)  ! number of omega states
      parameter (minlam=mindel+numdel)   ! ID of lowest lambda state
      parameter (maxlam=minlam+numlam-1) ! ID of highest lambda state
      parameter (minsig=minlam+numlam)   ! ID of lowest sigma state
      parameter (maxsig=minsig+numsig-1) ! ID of highest sigma state
      parameter (mincas=minsig+numsig)   ! ID of lowest cascade state
      parameter (maxcas=mincas+numcas-1) ! ID of highest cascade state
      parameter (minome=mincas+numcas)   ! ID of lowest omega state
      parameter (maxome=minome+numome-1) ! ID of highest omega state
c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc) ! ID of lowest baryon state
      parameter (maxbar=maxome) ! ID of highest baryon state
      parameter (offmeson=minmes) ! offset between zero and lowest 
                                  ! meson state
      parameter (maxmeson=maxmes) ! ID of highest meson state
c... these variables are in principal obsolete and should be exchanged 
c were referenced 
c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)   ! ID of kaon
      parameter (itrho=104)    ! ID of rho meson 
      parameter (itome=103)    ! ID of omega meson
      parameter (iteta=102)    ! ID of eta
      parameter (itphi=109)    ! ID of phi
      parameter (itetapr=107)  ! ID of eta'
      parameter (pimeson=101)  ! ID of $\pi$
      parameter (nucleon=minnuc) ! ID of nucleon
      integer itmin,itmax
      parameter (itmin=minnuc)  ! lowest defined ID
      parameter (itmax=maxmes)  ! highest defined ID
      parameter (maxbra=11)  ! decay channels for $s=0$ baryon resonances
      parameter (maxbrm=25) ! decay channels for meson resonances
      parameter (maxbrs1=10)! decay channels for $s=1$ baryon resonances
      parameter (maxbrs2=3) ! decay channels for $s=2$ baryon resonances
       integer mlt2it(maxmes-minmes) ! meson IDs sorted by multipletts
      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)      ! offset for mass generation
      parameter (mresmin=1.0765d0)  ! minimum baryon resonance mass
      parameter (mresmax=5d0)       ! maximum baryon resonance mass
      character*45 versiontag
      common /versioning/ versiontag
      real*8 massres(minbar:maxbar),widres(minbar:maxbar)
      real*8 branmes(0:maxbrm,minmes+1:maxmes)
      real*8 branres(0:maxbra,minnuc+1:maxdel)
      real*8 branbs1(0:maxbrs1,minlam+1:maxsig)
      real*8 branbs2(0:maxbrs2,mincas+1:maxcas)
      integer Jres(minbar:maxbar)
      integer Jmes(minmes:maxmes)
      integer pares(minbar:maxbar),pames(minmes:maxmes)
      integer Isores(minbar:maxbar), Isomes(minmes:maxmes)
      integer brtype(4,0:maxbra),bmtype(4,0:maxbrm)
      integer bs1type(4,0:maxbrs1),bs2type(4,0:maxbrs2)
      real*8 massmes(minmes:maxmes)
      real*8 mmesmn(minmes:maxmes)
      real*8 widmes(minmes:maxmes)
      integer strres(minbar:maxbar),strmes(minmes:maxmes)
      integer lbr(0:maxbra,minnuc+1:maxdel)
      integer lbs1(0:maxbrs1,minlam+1:maxsig)
      integer lbs2(0:maxbrs2,mincas+1:maxcas)
      integer lbm(0:maxbrm,minmes+1:maxmes)
      common /resonances/ massres,widres,massmes,widmes,mmesmn,
     ,                    branres,branmes,branbs1,branbs2,
     ,                    bs1type,bs2type,lbs1,lbs2,lbm,
     ,                    jres,jmes,lbr,brtype,pares,pames,
     ,                    bmtype,
     ,                    Isores,Isomes,strres,strmes,mlt2it
      !epos common blocks for particle list
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list
      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)
      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
      integer inbxxx
      real rinptl(mxptl),vrad,taugm,rangen
      common/c6ptl/rinptl,vrad,inbxxx
      integer nptl
      integer nevt,kolevt,koievt,npjevt,ntgevt,npnevt,nppevt,
     *        ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt,
     *        nglevt,minfra,maxfra
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
      integer istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer nn,idpdgg,idepos  ,nptlep
      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      integer pdgid
      external pdgid
      real*8 dectim
      external dectim
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      !integer nptest, ia
      !integer jc(nflav,2),ic(2),jca(3),i,nflav
      !parameter (nflav=6)
      integer itypart(nmax),iorpart(nmax)
      common /city/itypart /corpart/iorpart

c######################################   subroutine uexit      !transfer -> EPOS

      if(ish.ge.2)call alist('(u) final ptls&',0,0)
      if(ish.ge.2)write(ifch,*)' '
      
      !write(ifmt,'(a)')' '
      !nptest=0
      !do i=1,3
      !jca(i)=0
      !enddo
      nptlep=nptl
      do nn=1,npart  
	   idpdgg=pdgid(ityp(nn),iso3(nn))
           idepos=idtrafo('pdg','nxs',idpdgg)
              
	       nptl=nptl+1
               
               iorptl(nptl)=0
               jorptl(nptl)=0
               istptl(nptl)=0
               ifrptl(1,nptl)=0
               ifrptl(2,nptl)=0
               xorptl(4,nptl)=frr0(nn)
               xorptl(1,nptl)=frrx(nn)
               xorptl(2,nptl)=frry(nn)
               xorptl(3,nptl)=frrz(nn)
               pptl(4,nptl)=frp0(nn)
               pptl(1,nptl)=frpx(nn)
               pptl(2,nptl)=frpy(nn)
               pptl(3,nptl)=frpz(nn)
               pptl(5,nptl)=fmass(nn)
               idptl(nptl)=idepos
               istptl(nptl)=0
               ityptl(nptl)=itypart(nn)
               iorptl(nptl)=iorpart(nn)
               if(ityptl(nptl).eq.61)ityptl(nptl)=60
               if(ityptl(nptl).eq.62)ityptl(nptl)=60
               rinptl(nptl)=-9999
               tivptl(1,nptl)=xorptl(4,nptl)
               call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
               tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
               if(ish.ge.2)call alist('&',nptl,nptl)
      enddo
      !write(ifmt,'(i5)')nptest
      !write(ifmt,'(3i5)')jca
      
      write(ifmt,'(a,i6,a,i6,a)')                 
     . 'return ptls',nptlep+1,'  -',nptl
       call treatFinalParticles(nptlep+1,nptl)
     
      end


         subroutine  uepos


cc---------------------------------------------------------------------
c   main modul
cc---------------------------------------------------------------------

      implicit none
      integer idpdgg,idepos,pdgid,idtrafo
      real amepos
      include 'coms.f'
      include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx,iret
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer nuj
      data nuj/0/
      save nuj
      integer iuskip
      common /ciuskip/iuskip
      integer iutime,kk
      double precision ddt
      COMMON /NPTLC/ddt(ncollmax,4), iutime(5)
      logical ktime         

c#######################################    subroutine  uepos

      if(iuskip.eq.1)goto 777
      if(iuskip.eq.2)goto 778

      nuj=nuj+1

      mc=0
      mp=0
      noc=0
      kk=0
      ktime=.false.
     
      time = 0.0  !time is the system time at the BEGINNING of every timestep

      !initialize random number generator
      !call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

      !old time if an old fort.14 is used 
      if(CTOption(40).eq.1)time=acttime
      if(CTOption(40).eq.3)time=acttime

      !write headers to file
      call output(13)
      !call output(14)
      call output(15)
      call output(16)
      !if(event.eq.1)call output(17)
      call osc99_event(-1)

      !for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)
      !participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

      otime = outsteps*dtimestep  !compute time of output

      steps = 0  !reset time step counter

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! loop over all timesteps
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 20  steps=1,nsteps  

         if (eos.ne.0) then
           do j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
           enddo
         end if

         !we are at the beginning of the timestep, set current time (acttime)
         acttime = time

         if(CTOption(16).ne.0) goto 103  !option for MD without collision term

         call colload  ! Load collision table with next collisions in current timestep

         ! check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then
 101        continue              !entry-point for collision loop in case  
            k = 0                 !      of full colload after every coll
 100        continue              !normal entry-point for collision loop 
            kk=kk+1 
            if(ktime)call timer(iutime)  
            if(ktime)ddt(kk,1)=iutime(3)+0.001*iutime(4)
            call getnext(k)       !get next collision
            if (k.eq.0) goto 102  !exit collision loop if no collisions are left

            !propagate all particles to next collision time
            !store actual time in acttime, propagation time st=cttime(k)-acttime
	    st=cttime(k)-acttime
            if(ktime)call timer(iutime)  
            if(ktime)ddt(kk,2)=iutime(3)+0.001*iutime(4)
            call cascstep(acttime,st)
            acttime = cttime(k)   !new actual time (for upcoming collision)

            !perform collision 

            if(cti2(k).gt.0.and.
     .           abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(ifmt,*)' ***(E) wrong collision update (col) ***'
               write(ifmt,*)cti1(k),cti2(k),
     .              ctsqrts(k),sqrts(cti1(k),cti2(k))
            else if(cti2(k).eq.0.and.
     .              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
              write(ifmt,*)' *** main(W) wrong collision update (decay)'
              write(ifmt,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     .              fmass(cti1(k)),ctsqrts(k)
            endif

            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

            !store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))
            if(cti2(k).eq.0)then !~~~~~~~to avoid rare crash~~~~~~~~
              call checkdecay(ctsqrts(k),it1,iso3(cti1(k)),iret)
              if(iret.eq.1)then
	        idpdgg=pdgid(ityp(cti1(k)),iso3(cti1(k)))
                idepos=idtrafo('pdg','nxs',idpdgg)
                call idmass(idepos,amepos)
                ctsqrts(k)=amepos
                fmass(cti1(k))=amepos
                !print*,'(info) checkdecay: iret=1 --> take epos mass'
                call checkdecay(ctsqrts(k),it1,iso3(cti1(k)),iret)
                if(iret.eq.1)stop'checkdecay: still iret=1\n\n'
              endif
            endif !~~~~~~~~~~~~~~~~
            !increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter
               mc=mc+1
            endif
            !perform scattering/decay
            cti1sav = cti1(k)                       
            cti2sav = cti2(k)   
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))
            
            !update collision table 

            !normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
                 !new collision partners for pauli-blocked states (nexit=0)
                 if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                   cti1(k) = cti1sav 
                   cti2(k) = cti2sav 
                 endif

                 if(ktime)call timer(iutime)  
                 if(ktime)ddt(kk,3)=iutime(3)+0.001*iutime(4)
                 call collupd(cti1(k),1)
                 if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                 ncharge=0
                 !new collision partners for scattered/produced particles (nexit><0)
                 if(ktime)call timer(iutime)  
                 if(ktime)ddt(kk,3)=iutime(3)+0.001*iutime(4)
                 do i=1,nexit
                   !ncharge is used for charge conservation check
                   ncharge=ncharge+charge(inew(i))
                   call collupd(inew(i),1)
                 enddo
               endif
               !update collisions for partners of annihilated particles
               do ii=1,nsav
                  call collupd(ctsav(ii),1)
               enddo
               nsav=0
            else ! (CTOption(17).ne.0)
              !full collision load
              call colload
            endif
            if(ktime)call timer(iutime)  
            if(ktime)ddt(kk,4)=iutime(3)+0.001*iutime(4) 
            if(ktime)print*,'*****',kk
     .       ,nint((ddt(kk,2)-ddt(kk,1))*1000)
     .       ,nint((ddt(kk,3)-ddt(kk,2))*1000)
     .       ,nint((ddt(kk,4)-ddt(kk,3))*1000)

            if (CTOption(17).eq.0) goto 100
            goto 101

            !this is the point to jump to after all collisions in the timestep
            !have been taken care of
 102        continue

         endif ! (nct.gt.0)


         !After all collisions in the timestep are done, propagate to end of 
         !the timestep.

         !point to jump to in case of MD without collision term
 103     continue

         time = time+dtimestep  !increment timestep

         !After all collisions in the timestep are done, propagate to end of 
         !the timestep.
         call cascstep(acttime,time-acttime)

         !in case of potential interaction, do MD propagation step
         if (eos.ne.0) then
            ! set initial conditions for MD propagation-step
            do j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
            enddo
            !now molecular dynamics trajectories
            call proprk(time,dtimestep)
         endif ! (eos.ne.0)

         !perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
            if(CTOption(28).eq.2)call spectrans(otime)
            call file14outx(steps/outsteps)
         endif ! output handling

 20   continue ! time step loop

      acttime=time
      
      if(npart.eq.0)stop
     . '\n \n    ***** STOP in uepos: no particles (2) *****  \n\n'
      
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! optional decay of all unstable 
  !  particles before final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  777 continue
      
      !DANGER: pauli-blocked decays are not performed !!!
      if(CTOption(18).eq.0.and.CTOption(51).eq.0
     .   .or.iuskip.eq.1 ) then
         !print*,'(info) npart=',npart,' before final decay' !----------
         !no do-loop is used because npart changes in loop-structure
         i=0
         nct=0
         actcol=0
         CTOption(10)=1  !disable Pauli-Blocker for final decays
 40      continue  !decay loop structure starts here
         i=i+1
         if(dectime(i).lt.1.d30) then !if particle unstable
 41         continue
            isstable = .false.
            do stidx=1,nstable
               if (ityp(i).eq.stabvec(stidx)) then
                  !write (6,*) 'no decay of particle ',ityp(i)
                  isstable = .true.
               endif
            enddo
            if (.not.isstable) then
               !~~~~~~~to avoid rare crash~~~~~~~~
               call checkdecay(fmass(i),ityp(i),iso3(i),iret)
               if(iret.eq.1)then
	         idpdgg=pdgid(ityp(i),iso3(i))
                 idepos=idtrafo('pdg','nxs',idpdgg)
                 call idmass(idepos,amepos)
                 fmass(i)=amepos
                 print*,'(info) checkdecay: iret=1 --> '
     .           ,'take epos mass (2)'
                 call checkdecay(fmass(i),ityp(i),iso3(i),iret)
                 if(iret.eq.1)stop'checkdecay: still iret=1 (2) \n\n'
               endif
               !~~~~~~~~~~~~~~~~
               call scatter(i,0,0.d0,fmass(i),xdummy) !perform decay
               !backtracing if decay-product is unstable itself
               if(dectime(i).lt.1.d30) goto 41
            endif
         endif
         !check next particle
         if(i.lt.npart) goto 40
         !print*,'(info) npart=',npart,' after final decay'  !---------------
      endif ! final decay

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !     final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(CTOption(28).eq.2)call spectrans(otime)

      call file13out(nsteps)
      !call file14out(nsteps)
      call file16out
      call osc_event
      call osc99_event(1)
      call osc99_eoe
      
      mp=mp+npart
      if(ctag.eq.0)then
         !!!!!write(*,*)'(W) No collision in event ',event
         noc=noc+1
      endif

      if(nuj.eq.1)write(ifmt,'(a,50a1)')' (info)',('u',i=1,50)
 
  778 continue
  
      end


c###################################################################################
c###################################################################################

           subroutine input(io)

c###################################################################################
c###################################################################################

c-----------------------------------------------------------------------
c     This subroutine reads the UQMD input file (unit=9) 
c
c input : for ({\\tt io=0} input-file will be processed, 
c         else default values assumed 
c output: information in common-block coms.f
c
c-----------------------------------------------------------------------

      implicit none

      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'inputs.f'
      include 'boxinc.f'
      
      integer laproj,maproj,latarg,matarg
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      integer icinpu
      real engy,elepti,elepto,angmue
      common/lept1/engy,elepti,elepto,angmue,icinpu

      character*3 flag
      character*77 inputstr,file9,fheader,file14,file15,file16,file17
      character*77 file13,file10,file19,file20
      integer line,proflg,tarflg,impflg,beamflg,inx,ival,partid
      integer eosflg,i,io
      real*8 rval,caltim,outtim
      logical dtflag,bret

      character CTOStrng(numcto)*60
      character CTPStrg(numctp)*60

c setting of internal parameters values:
      real*8 valint(1)
      common /values/ valint
      
      logical infu

      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      
      integer ncnt
      data ncnt /0/

      save

c##############################################################################
c######################################  subroutine input(io)
c##############################################################################

      ncnt=ncnt+1
      infu=info
      if(ncnt.gt.1)infu=.false.


      valint(1)=0.d0    

      bret=io.ne.0
      goto 108

      entry inpini  
c  called by some test programs

      bret=.true.
 108  continue
  
c initialize counters
      line=0
      boxflag=0
      mbflag=0
      edens=0.d0
      para=0
      solid=0
      mbox=0

c the following flags check, wether all necessary input is given 
c projectile
      proflg=0
      prspflg=0
c target
      tarflg=0
      trspflg=0
c impact parameter
      impflg=0
c incident beam energy
      beamflg=0
      srtflag=0
      firstev=0
c equation of state
      eosflg=0
c excitation function
      nsrt=1
	 npb=1
      efuncflag=0
c default number of events
      nevents=1
c default seed for random number generator
      ranseed=0
c default number of timesteps
      nsteps=1000
c use standard time-step
      dtflag=.false.
c skip conditions on unit 14, 15, 16 & 18
      bf13=.false.
      bf14=.false.
      bf15=.false.
      bf16=.false.
      bf18=.false.
      bf19=.false.
      bf20=.false.
      do 111 i=1,numcto
         CTOdc(i)='  '
 111  continue
      do 112 i=1,numctp
         CTPdc(i)='  '
 112  continue
      do 113 i=1,maxstables
         stabvec(i)=0
 113  continue
      nstable = 0

c default settings for CTParam and CTOption cccccccccccccccccccccccccccccc
      CTParam(1)=1.d0  
      CTPStrg(1)='scaling factor for decay-width'
      CTParam(2)=0.52d0 
      CTPStrg(2)='used for minimal stringmass & el/inel cut in makestr'
      CTParam(3)=2d0 
      CTPStrg(3)='velocity exponent for modified AQM'  
      CTParam(4)=0.3d0 
      CTPStrg(4)='transverse pion mass, used in make22 & strexct'
      CTParam(5)=0d0 
      CTPStrg(5)='probabil. for quark rearrangement in cluster'
      CTParam(6)=0.4d0  
      CTPstrg(6)='strangeness probability'
      CTParam(7)=0.d0 
      CTPStrg(7)='charm probability (not yet implemented in UQMD)'
      CTParam(8)=0.093d0 
      CTPStrg(8)='probability to create a diquark'
      CTParam(9)=0.35d0 
      CTPStrg(9)='kinetic energy cut off for last string break'
      CTParam(10)=0.25d0 
      CTPStrg(10)='min. kinetic energy for hadron in string'
      CTParam(11)=0.d0 
      CTPStrg(11)='fraction of non groundstate resonances'
      CTParam(12)=.5d0  
      CTPStrg(12)='probability for rho 770 in String'
      CTParam(13)=.27d0 
      CTPStrg(13)='probability for rho 1450 (rest->rho1700)'
      CTParam(14)=.49d0 
      CTPStrg(14)='probability for omega 782'
      CTParam(15)=.27d0 
      CTPStrg(15)='probability for omega 1420(rest->om1600)'
      CTParam(16)=1.0d0 
      CTPStrg(16)='mass cut betw. rho770 and rho 1450'
      CTParam(17)=1.6d0 
      CTPSTRG(17)='mass cut betw. rho1450 and rho1700'
      CTParam(18)=.85d0 
      CTPStrg(18)='mass cut betw. om 782 and om1420'
      CTParam(19)=1.55d0
      CTPStrg(19)='mass cut betw. om1420 and om1600'
      CTParam(20)=0.0d0
      CTPStrg(20)=' distance for second projectile'
      CTParam(21)=0.0d0
      CTPStrg(21)=' deformation parameter'
      CTParam(25)=.9d0 
      CTPStrg(25)=' probability for diquark not to break'
      CTParam(26)=50d0 
      CTPStrg(26)=' maximum trials to get string masses'
      CTParam(27)=1d0 
      CTPStrg(27)=' scaling factor for xmin in string excitation'
      CTParam(28)=1d0 
      CTPStrg(28)=' scaling factor for transverse fermi motion'
      CTParam(29)=0.4 
      CTPStrg(29)=' single strange di-quark suppression factor '
      CTParam(30)=1.5 
      CTPStrg(30)=' radius offset for initialisation  '
      CTParam(31)=1.6d0 
      CTPStrg(31)=' sigma of gaussian for tranverse momentum tranfer '
      CTParam(32)=0d0
      CTPStrg(32)=' alpha-1 for valence quark distribution  '
      CTParam(33)=2.5d0
      CTPStrg(33)=' betav for valence quark distribution  (DPM)'
      CTParam(34)=0.1
      CTPStrg(34)=' minimal x multiplied with ecm  '
      CTParam(35)=3.0
      CTPStrg(35)=' offset for cut for the FSM '
      CTParam(36)=0.275d0
      CTPStrg(36)=' fragmentation function parameter a  '
      CTParam(37)=0.42d0
      CTPStrg(37)=' fragmentation function parameter b  '
      CTParam(38)=1.08d0
      CTPStrg(38)=' diquark pt scaling factor '
      CTParam(39)=0.8d0
      CTPStrg(39)=' strange quark pt scaling factor '
      CTParam(40)=0.5d0
      CTPStrg(40)=' betas-1 for valence quark distribution (LEM)'
      CTParam(41)=0.0
      CTPStrg(41)=' distance of initialisation'
      CTParam(42)=0.55d0
      CTPStrg(42)=' width of gaussian -> pt in string-fragmentation '
      CTParam(43)=5.d0
      CTPStrg(43)=' maximum kinetic energy in mesonic clustr '
      CTParam(44)=.8d0
      CTPStrg(44)=' prob. of double vs. single excitation for AQM inel.'
      CTParam(45)=0.5
      CTPStrg(45)=' offset for minimal mass generation of strings'
      CTParam(46)=800000
      CTPStrg(46)=' maximal number of rejections for initialisation'
      CTParam(47)=1.0
      CTPStrg(47)=' field feynman fragmentation funct. param. a'
      CTParam(48)=2.0
      CTPStrg(48)=' field feynman fragmentation funct. param. b'
      CTParam(49)=50.5
      CTPStrg(49)=' Energy cut-off for exclusive PYTHIA use'
      CTParam(50)=1d0 
      CTPStrg(50)=' enhancement factor for 0- mesons'
      CTParam(51)=1d0 
      CTPStrg(51)=' enhancement factor for 1- mesons'
      CTParam(52)=1d0
      CTPStrg(52)=' enhancement factor for 0+ mesons'
      CTParam(53)=1d0
      CTPStrg(53)=' enhancement factor for 1+ mesons'   
      CTParam(54)=1d0 
      CTPStrg(54)=' enhancement factor for 2+ mesons'   
      CTParam(55)=1d0
      CTPStrg(55)=' enhancement factor for 1+-mesons'   
      CTParam(56)=1d0
      CTPStrg(56)=' enhancement factor for 1-*mesons'   
      CTParam(57)=1d0
      CTPStrg(57)=' enhancement factor for 1-*mesons'    
      CTParam(58)=1.d0
      CTPStrg(58)=' scaling factor for DP time-delay'

cc
      CTOption(1)=0  
      CTOStrng(1)=' resonance widths are mass dependent '
      CTOption(2)=0
      CTOStrng(2)=' conservation of scattering plane'
      CTOption(3)=0  
      CTOStrng(3)=' use modified detailed balance'
      CTOption(4)=0  
      CTOStrng(4)=' no initial conf. output '
      CTOption(5)=0  
      CTOStrng(5)=' fixed impact parameter'
      CTOption(6)=0  
      CTOStrng(6)=' no first collisions inside proj/target'
      CTOption(7)=0  
      CTOStrng(7)=' elastic cross section enabled (<>0:total=inelast)'
      CTOption(8)=0  
      CTOStrng(8)=' extrapolate branching ratios '
      CToption(9)=0  
      CTOStrng(9)=' use tabulated pp cross sections ' 
      CTOption(10)=0 
      CTOStrng(10)=' enable Pauli Blocker'
      CTOption(11)=0 
      CTOStrng(11)=' mass reduction for cascade initialisation' 
      CTOption(12)=0 
      CTOStrng(12)=' string condition =0 (.ne.0 no strings)'
      CTOption(13)=0 
      CTOStrng(13)=' enhanced file16 output '
      CTOption(14)=0 
      CTOStrng(14)=' cos(the) is distributet between -1..1 '
      CTOption(15)=0 
      CTOStrng(15)=' allow mm&mb-scattering'
      CTOption(16)=0 
      CTOStrng(16)=' propagate without collisions'
      CTOption(17)=0 
      CTOStrng(17)=' colload after every timestep '
      CTOption(18)=0 
      CTOStrng(18)=' final decay of unstable particles'
      CTOption(19)=0  
      CTOStrng(19)=' allow bbar annihilaion'
      CTOption(20)=0
      CTOStrng(20)=' dont generate e+e- instead of bbar'
      CTOption(21)=0
      CTOStrng(21)=' use field feynman frgm. function'
      CTOption(22)=1
      CTOStrng(22)=' use lund excitation function'
      CTOption(23)=0
      CTOStrng(23)=' lorentz contraction of projectile & targed'
      CTOption(24)=1
      CTOStrng(24)=' Wood-Saxon initialization'
      CTOption(25)=0
      CTOStrng(25)=' phase space corrections for resonance mass'
      CTOption(26)=0
      CTOStrng(26)=' use z -> 1-z for diquark-pairs'
      CTOption(27)=0 
      CTOStrng(27)=' reference frame (1=target, 2=projectile, else=cms)'
      CTOption(28)=0
      CTOStrng(28)=' propagate spectators also '
      CTOption(29)=2
      CTOStrng(29)=' no transverse momentum in clustr '
      CTOption(30)=1
      CTOStrng(30)=' frozen fermi motion '
      CTOption(31)=0
      CTOStrng(31)=' reduced mass spectrum in string'
      CTOption(32)=0
      CTOStrng(32)=' masses are distributed acc. to m-dep. widths'
      CTOption(33)=0
      CTOStrng(33)=' use tables & m-dep. for pmean in fprwdt & fwidth'
      CTOption(34)=1
      CTOStrng(34)=' lifetme according to m-dep. width'
      CTOption(35)=1
      CTOStrng(35)=' generate high precision tables'
      CTOption(36)=0
      CTOStrng(36)=' normalize Breit-Wigners with m.dep. widths '
      CTOption(37)=1
      CTOStrng(37)=' heavy quarks form di-quark clusters'
      CTOption(38)=0
      CTOStrng(38)=' scale p-pbar to b-bbar with equal p_lab '
      CTOption(39)=0
      CTOStrng(39)=' dont call pauliblocker'
      CTOption(40)=0
      CTOStrng(40)=' read old fort.14 file '
      CTOption(41)=0
      CTOStrng(41)=' generate extended output for cto40'
      CTOption(42)=0
      CTOStrng(42)=' hadrons now have color fluctuations'
      CTOption(43)=0
      CTOStrng(43)=" don't generate dimuon intead of dielectron output"
      CTOption(44)=0
      CTOStrng(44)=' PYTHIA for hard scatterings also at low energies'
      CTOption(45)=0
      CTOStrng(45)=' not used at the moment'

      CTOption(50)=0
      CTOStrng(50)=' two-body rescatterings are not allowed'
      CTOption(51)=0
      CTOStrng(51)=' decays are not allowed'
      
      if(bret)return

c initialize arrays for special PRO/TAR combinations
      do 10 i=1,2
         spityp(i)=0
         spiso3(i)=0
 10   continue
c header for output files
      fheader=' this is the default uqmd-fileheader'

c open fortran-unit 9 for input
c and units 14, 15 for output
      call getenv('ftn09',file9)
      call getenv('ftn10',file10)
      call getenv('ftn13',file13)
      call getenv('ftn14',file14)
      call getenv('ftn15',file15)
      call getenv('ftn16',file16)
      call getenv('ftn17',file17)
      call getenv('ftn19',file19)
      call getenv('ftn20',file20)
ckw      if (file9(1:4).ne.'    ') then
ckw         OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')
ckw      endif
      if (file10(1:4).ne.'    ') then
         OPEN(UNIT=10,FILE=file10,STATUS='OLD',FORM='FORMATTED')
         CTOption(40)=1
         nevents=100000
      endif
      if (file13(1:4).ne.'    ') then
         OPEN(UNIT=13,FILE=file13,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file14(1:4).ne.'    ') then
         OPEN(UNIT=14,FILE=file14,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file15(1:4).ne.'    ') then
         OPEN(UNIT=15,FILE=file15,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file16(1:4).ne.'    ') then
         OPEN(UNIT=16,FILE=file16,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file17(1:4).ne.'    ') then
         OPEN(UNIT=17,FILE=file17,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file19(1:4).ne.'    ') then
         OPEN(UNIT=19,FILE=file19,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file20(1:4).ne.'    ') then
         OPEN(UNIT=20,FILE=file20,STATUS='unknown',FORM='FORMATTED')
      endif
c
c stop input if old event is read in
      if(CTOption(40).eq.1) return


c this entry is used to read cto,ctp and tim statements
c in case of old event readin
      entry getparams
      
ckw      close(9)
ckw      OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')
 
c read input lines
 1    continue
      line=line+1
ckw      read(9,99) flag,inputstr
       
      inputstr( 1:40)='                                        '
      inputstr(41:77)='                                     '
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !    settings
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(line.eq.1)then
        flag='pro'                         !~~~~projectile (Ap, Zp)
        write(inputstr(1:8),'(2i4)')maproj,laproj
      elseif(line.eq.2)then              
        flag='tar'                         !~~~~target (Ap, Zp)
        write(inputstr(1:8),'(2i4)')matarg,latarg
      elseif(line.eq.3)then     
        flag='nev'                         !~~~~number of events
        inputstr(1:2)=' 1'
      elseif(line.eq.4)then 
        flag='tim'    !~~~~propagation time, output time step (in fm/c)
        inputstr(1:8)=' 400 400'
      elseif(line.eq.5)then 
        flag='ecm'                         !~~~~cms energy in AGeV
        write(inputstr(1:5),'(i5)')nint(engy)
      elseif(line.eq.6)then 
        flag='imp'                         !~~~~impact parameter (in fm)
        inputstr(1:2)=' 0'
      elseif(line.eq.7)then 
        flag='   '                         !~~~~random number seed
        !inputstr(1:11)=' 1134570653'
      elseif(line.eq.8)then 
        flag='eos'                         !~~~~equation of state
        inputstr(1:2)=' 0'
      elseif(line.eq.9)then 
        flag='cto'                     
        inputstr(1:4)=' 5 1'
      elseif(line.eq.10)then 
        flag='cto'                     
        inputstr(1:5)=' 40 3'
      elseif(line.eq.11)then 
        flag='f13'                     
        inputstr(1:1)=' '
      elseif(line.eq.12)then 
        flag='   '   !'f14'             !~~~~~~~~~suppress output                 
        inputstr(1:1)=' '
      elseif(line.eq.13)then 
        flag='f15'                     
        inputstr(1:1)=' '
      elseif(line.eq.14)then 
        flag='f16'                     
        inputstr(1:1)=' '
      elseif(line.eq.15)then 
        flag='f19'                     
        inputstr(1:1)=' '
      elseif(line.eq.16)then 
        flag='f20'                     
c      elseif(line.eq.17)then !~~~~Pauli blocker 
c        flag='cto'                     
c        inputstr(1:5)=' 10 1'
c      elseif(line.eq.??)then  !~~~~~for fast run comment the following flags!!!!
c        flag='cdt'                 !~~~~~~~~~  timestep        
c        inputstr(1:1)='1'
c      elseif(line.eq.??)then 
c        flag='tim'                 !~~~~~~~~~ caltim, outtim    
c        inputstr(1:6)='400 20'
      else
        goto2
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c 3    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c select action according to flag:
c #  : treat line as a comment
      if(flag(1:1).eq.'#') goto 1
c blanks are comments, too
      if(flag(1:1).eq.' ') goto 1
c xxx: treat line as end of input marker
      if(flag.eq.'xxx'.or.flag.eq.'end') then
         goto 2
cc cal: header for output-files
c      if(flag.eq.'cal') then
c         fheader=inputstr
c pro: define projectile
      elseif(flag.eq.'pro') then
         proflg=proflg+1
         read(inputstr,fmt=*,err=88,end=88) Ap,Zp
         if(proflg.gt.1) then
            write(ifmt,*)'multiple definitions for projectile system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c PRO: define special projectile
      elseif(flag.eq.'PRO') then
         proflg=proflg+1
         prspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(1),spiso3(1)
         Ap=1
         if(proflg.gt.1) then
            write(ifmt,*)'multiple definitions for projectile system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c tar: define target
      elseif(flag.eq.'tar') then
         tarflg=tarflg+1
         read(inputstr,fmt=*,err=88,end=88) At,Zt
         if(tarflg.gt.1) then
            write(ifmt,*)'multiple definitions for target system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c TAR: define special target
      elseif(flag.eq.'TAR') then
         tarflg=tarflg+1
         trspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(2),spiso3(2)
         At=1
         if(tarflg.gt.1) then
            write(ifmt,*)'multiple definitions for target system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c box: define a box with a length in fm
c	parameters: 2: energie
c		    3: 1 =solid		
c		    4: 1 = walls
        elseif(flag.eq.'box') then          
           boxflag=boxflag+1                     
          read(inputstr,fmt=*,err=88,end=88) lbox,edens,solid,para
	   if (edens.gt.0.d0) edensflag=1
		
           if (lbox.le.0) then
              write(ifmt,*) 'Error, lenght<=0'
              stop                          
           endif                            
           lboxhalbe=lbox/2.d0
           lboxd=lbox*2.d0

 	   if (edens.lt.0.d0) then 
	      write(ifmt,*) 'Error, a negativ energy '
	      stop
	   endif
           
           if(boxflag.gt.1) then            
            write(ifmt,*)'multiple boxes are defined'
            stop                            
        endif                  
c bpt: define particles in the box
c parameters: ityp, iso3, mpart, pmax
	elseif(flag.eq.'bpt') then
	   if (edens.gt.0.d0) then
	      write(ifmt,*) 'Error, energie is already defined'
	      stop
	   endif
           mbox=mbox+1
           read(inputstr,fmt=*,err=88,end=88) 
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox),bptpmax(mbox)
	   edensflag=0 
 	   if (bptpart(mbox).le.0) then 
	      write(ifmt,*) 'Error, a negativ particle number'
	      stop
	   endif
           if(boxflag.lt.1) then
            write(ifmt,*)'no box is defined'          
	    stop
        endif
c bpe: define particles in the box with a given energy
c parameters: ityp, iso3, mpart, 
	elseif(flag.eq.'bpe') then
	   if (edens.le.0) then
	      write(ifmt,*) 'Error, no energie is defined'
	      stop
	   endif
           mbox=mbox+1
	   read(inputstr,fmt=*,err=88,end=88) 
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox)
           if(boxflag.lt.1) then
            write(ifmt,*)'no box is defined'          
	    stop
        endif
c ene: beam energy (lab-system)
      elseif(flag.eq.'ene'.or.flag.eq.'elb') then
         beamflg=beamflg+1
         read(inputstr,fmt=*,err=88,end=88) ebeam 
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if (ebeam.le.200) then
           write(ifmt,*)'Calculation at ebeam.le.200 A GeV:'
           write(ifmt,*)'parameter nmax in coms.f may be decreased!'
         endif
c plb: beam momentum (lab-system)
      elseif(flag.eq.'plb') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbeam 
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
       if (pbeam.le.200) then
            write(ifmt,*)'Calculation at pbeam.le.200 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
       endif
c PLB: beam momentum ( LAb-system, excitation function possible)
      elseif(flag.eq.'PLB'.or.flag.eq.'PLG') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbmin,pbmax,npb 
         pbeam=pbmin
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if(npb.gt.1.and.flag.eq.'PLB') efuncflag=1
         if(npb.gt.1.and.flag.eq.'PLG') efuncflag=2
         if(abs(pbmax-pbmin).le.1.d-6) then
            npb=1
            efuncflag=0
         endif
         if (pbmax.le.200) then
            write(ifmt,*)'Calculations at pbmax.le.200 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
         endif
c ecm:  c.m.energy 
      elseif(flag.eq.'ecm') then
         beamflg=beamflg+1
         srtflag=1
         read(inputstr,fmt=*,err=88,end=88) ecm 
         srtmin=ecm
         srtmax=ecm
         nsrt=1
         efuncflag=0 
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if (ecm.le.20) then 
          if(infu)then
        write(ifmt,*)'(info) Calculation at sroot.le.20 A GeV:'
        write(ifmt,*)'(info) parameter nmax in coms.f may be decreased!'
          endif
         endif
c ENE: beam energy (sqrt(s): CM-system, excitation function possible)
      elseif(flag.eq.'ENE'.or.flag.eq.'ELG') then
         beamflg=beamflg+1
         srtflag=1
         read(inputstr,fmt=*,err=88,end=88) srtmin,srtmax,nsrt 
         ecm=srtmin
c        if(flag.eq.'ELG')ecm=1d1**dlog10(srtmin)
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if(nsrt.gt.1.and.flag.eq.'ENE') efuncflag=1
         if(nsrt.gt.1.and.flag.eq.'ELG') efuncflag=2
         if(abs(srtmax-srtmin).le.1.d-6) then
            nsrt=1
            efuncflag=0
         endif
         if (srtmax.le.20) then
            write(ifmt,*)'Calculations at srootmax.le.20 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
         endif
c imp: impact parameter
      elseif(flag.eq.'imp') then
         bmin=0.d0
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bdist 
         if(bdist.lt.0d0)then
           CTOption(5)=1
           bdist=abs(bdist)
           write(ifmt,*)'randomly choosen impact parameter:',
     ,             ' CTOption(5) is set to 1'
         end if
         if(impflg.gt.1) then
            write(ifmt,*)'multiple definitions for impact parameter:'
            write(ifmt,*)'-> last entry will be used'
         endif
c IMP: impact parameter
      elseif(flag.eq.'IMP') then
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bmin,bdist 
         CTOption(5)=1
         if(impflg.gt.1) then
            write(ifmt,*)'multiple definitions for impact parameter:'
            write(ifmt,*)'-> last entry will be used'
         endif
c eos: impact parameter
      elseif(flag.eq.'eos') then
         eosflg=eosflg+1
         read(inputstr,fmt=*,err=88,end=88) eos 
         if(eosflg.gt.1) then
            write(ifmt,*)'multiple definitions for equation of state:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if (eos.ne.0) then
            CTOption(24)=0
         endif
c nev: number of events
      elseif(flag.eq.'nev') then
         read(inputstr,fmt=*,err=88,end=88) nevents 
c rsd: 
      elseif(flag.eq.'rsd') then
         read(inputstr,fmt=*,err=88,end=88) ranseed
c cdt: collision time step
      elseif(flag.eq.'cdt') then
         read(inputstr,fmt=*,err=88,end=88) dtimestep
         dtflag=.true.
c tim: time of propatation
      elseif(flag.eq.'tim') then
         read(inputstr,fmt=*,err=88,end=88) caltim, outtim 
c stb: keep particle stable
      elseif(flag.eq.'stb') then
         read(inputstr,fmt=*,err=88,end=88) partid
         if (nstable.lt.maxstables) then
            nstable = nstable + 1
            stabvec(nstable) = partid
         else
            write(ifmt,*) 'Warning: too many stable particles defined!'
         endif
c cto: collision term options
      elseif(flag.eq.'cto') then
         read(inputstr,fmt=*,err=88,end=88) inx,ival
         if(ncnt.eq.1)
     &    write(ifmt,*)'(info) CTOption(',inx,')=',CTOption(inx)
     &        ,CTOStrng(inx)(1:index(CTOStrng(inx),' '))
     &      ,'is changed to',ival
         CTOption(inx)=ival
         CTOdc(inx)=' *'
c ctp: collision term parameter
      elseif(flag.eq.'ctp') then
         read(inputstr,fmt=*,err=88,end=88) inx,rval
         CTParam(inx)=rval
         CTPdc(inx)=' *'
         write(ifmt,*)'CTParam(',inx,'):   ',CTPStrg(inx)
     ,             ,'is changed to',rval
      elseif (flag.eq.'f13') then
         bf13=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 13'
      elseif (flag.eq.'f14') then
         bf14=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 14'
      elseif (flag.eq.'f15') then
         bf15=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 15'
      elseif (flag.eq.'iou') then
         read(inputstr,fmt=*,err=88,end=88) inx,ival
         call uounit(inx,ival)
         write(ifmt,*)'file',inx,'will be written on unit',ival
      elseif (flag.eq.'f16') then
         bf16=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 16'
      elseif (flag.eq.'f18') then
          bf18=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 18'
      elseif (flag.eq.'f19') then
          bf19=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 19'
      elseif (flag.eq.'f20') then
          bf20=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 20'
      else
         write(ifmt,*)'undefined opcode in input-file on line',line
         stop
      endif
      goto 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue


c fast CASCADE mode
      if(.not.dtflag.and.eos.eq.0) dtimestep=outtim
c
      nsteps=int(0.01+caltim/dtimestep)
      outsteps=int(0.01+outtim/dtimestep)
      if(infu)write(ifmt,*)'(info) nsteps, outsteps :'
     . ,nsteps, outsteps
      if(infu)write(ifmt,*)'(info) dtimestep :', dtimestep

c stop input if old event is read in
      if(CTOption(40).eq.1) return


c here some validity checks of the input should be performed
	if (boxflag.eq.1.and.mbox.eq.0) then
	    write(ifmt,*) 'Error: no particles in the box.'
	    stop
	ElseIf (boxflag.eq.0) then
      if(proflg.eq.0) then
         write(ifmt,*)'Error: no projectile specified in input.'
         stop
      elseif(tarflg.eq.0) then
         write(ifmt,*)'Error: no target specified in input.'
         stop
      elseif((impflg.eq.0)) then
         write(ifmt,*)'Error: no impact parameter in input.'
         stop
      elseif(beamflg.eq.0.and.prspflg.eq.0) then
         write(ifmt,*)'Error: no incident beam energy specified.'
         stop
      endif
c EndIf for the Box
	EndIf      
      if (efuncflag.ne.0.and.
     &    mod(nevents,max(nsrt,npb)).ne.0) then
         write(ifmt,*)'INPUT: the number of events divided by the ',
     ,   'number of energies requested is no integer.'
      end if      
c
c constraints for skyrme pots:
      if(eos.ne.0.and.((srtflag.eq.0.and.ebeam.gt.4d0)
     &             .or.(srtflag.eq.1.and.srtmax.gt.3.3d0)
     &             .or.(srtflag.eq.2.and.pbeam.gt.4.9))) then
         write(ifmt,*)'***(W) I switched off the potentials'
         eos=0
      end if
      if(eos.ne.0) then
         CTOption(11)=1
         CTOption(28)=0
         CTOption(30)=0
      endif
c

c now print the selected analysis

c...some input combinations should be avoided and/or commented
      if(CTOption(7).ne.0.and.At*Ap.ne.1)then
        write(ifmt,*)'Warning: CTOption(7)=',CTOption(7), 
     ,  ' no elastic collisions in NN',
     ,  ' should not be used for serious calculations!'
      end if

      if(CTOption(18).ne.0)then
        write(ifmt,*)'Warning: CTOption(18)=',CToption(18),': ',
     ,  'unstable particles will not decay after propagation.'
      end if


      if(CTOption(31).ne.0)then
        write(ifmt,*)'Warning: CTOption(31)=',CToption(31),': ',
     ,  "Not yet completly implemented. Don't use for serious",
     ,  'calculations (not yet..).' 
      end if

      if(CTParam(28).lt.0d0.or.CTParam(28).gt.1d0)then
        write(ifmt,*)'Warning: CTParam(28)=',CTParam(28),': ',
     ,  'should be between 0 and 1. it will be corrected.'
        CTParam(28)=min(1d0,max(0d0,CTParam(28)))
      end if
      
      return

 88   write(ifmt,*) 'syntax-error in input-file on line ',line
     .   ,'   flag ',flag
      write(ifmt,*)inputstr
      stop
      end


c######################################################################
c######################################################################
c##################
c##################        output
c##################
c######################################################################
c######################################################################


c--------------------------------------------------------------------------------
                  subroutine file14outx(itime)
c--------------------------------------------------------------------------------
      implicit none
      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'
      integer i,itotcoll,iinelcoll,ii,ix,iy,itime,ncnt
      integer nexpart,idpdgg,idepos,idtrafo,pdgid,neta,ij,nrap,nrapid
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      real*8 sigmatot,t,tf,z,x,y,p1,p2,p3,p4
      logical go
      common /outco2/sigmatot
      include 'outcom.f'
      data ncnt /0/
      save 

      if(bf14)return
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      ! print*,'(file14outx)',ttime,npart
      !@ ,itotcoll,NElColl,iinelcoll,NBlColl,dectag,
      !@     NHardRes,NSoftRes,NDecRes
       !---------------------------------------------------
       !  r0(i), rx(i), ry(i), rz(i)   ................... x4
       !  p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i)
       !                    ,pz(i)+ffermpz(i),fmass(i) ... p5
       !  ityp(i)  ..... particle id 
       !  iso3(i) ...... 2 times the isospin of a particle
       !  charge(i) .... charge of the particle
       !  lstcoll(i) ... index of last collision partner
       !  ncoll(i) ..... number of collisions
       !  origin(i) ....
       !  dectime(i) ...
       !  tform(i) ..... formation time
       !  xtotfac(i) ... cross section 
       !                 (zero if the particle is not yet formed)
       !  uid(i) ......
       !-----------------------------------------------------
       nexpart=0
       ncnt=ncnt+1
       if(ncnt.eq.1)then
       do neta=-6,6
        zetasum(neta,itime)=0
       enddo
       do nrap=-10,10
        zrapsum(nrap,itime)=0
       enddo
       do ii=1,2
       do ij=-10,10
       do neta=-10,10
        zeta(ii,ij,neta,itime)=0
       enddo
       enddo
       enddo
       endif
       do i=1,npart
         t=r0(i)
         tf=tform(i)
         if(t.gt.tf)then
          nexpart=nexpart+1
	  idpdgg=pdgid(ityp(i),iso3(i))
          idepos=idtrafo('pdg','nxs',idpdgg)
          ! some codes like 40323 are not in list -> idepos=0
          z=rz(i)
          x=rx(i)
          y=ry(i)
          ix=-10+(x+10.5)
          iy=-10+(y+10.5)
          p4=p0(i)
          p1=px(i)+ffermpx(i)
          p2=py(i)+ffermpy(i)
          p3=pz(i)+ffermpz(i)
          nrap=nrapid(p3,p4)
          zrapsum(nrap,itime)=zrapsum(nrap,itime)+1
          neta=nrapid(z,t)
          zetasum(neta,itime)=zetasum(neta,itime)+1
          do ii=1,2
           go=.false.
           if(ii.eq.1.and.abs(y).le.1.)go=.true.
           if(ii.eq.2.and.abs(x).le.1.)go=.true.
           if(go)then
            if(ii.eq.1)ij=ix
            if(ii.eq.2)ij=iy
            if(ij.ge.-10.and.ij.le.10)then
             zeta(ii,ij,neta,itime)=zeta(ii,ij,neta,itime)+1
            endif
           endif
          enddo
         endif  
       enddo
       
       !write(ifmt,'(a,2i3,i6)')'(file14outx)',itime,nsteps,nexpart
       
       !do ii=1,2
       !print*,' '
       !do neta=-6,6
       !write(ifmt,'(3x,13i4)')(zeta(ii,ij,neta,itime),ij=-5,5)
       !enddo
       !enddo

      end


c----------------------------------------------------------------------
            subroutine uplot
c----------------------------------------------------------------------

      include 'coms.f'
      !if(nsteps.gt.1)call uplot3
      end
      
      !------------------------
      subroutine uplot3
      !------------------------
      include 'aaa.h'
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!-------------------'
      write(ifhi,'(a,i3)') '!   uplot2     '
      write(ifhi,'(a)')    '!-------------------'
      print*,'nevent=',nevent
      do itime=1,31,2
       call getchar(itime,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -5 5'
       write(ifhi,'(a)')    'txt  "xaxis y "'              
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis dn/dy "'
       write(ifhi,'(a)')'array 2'
       do nrap=-5,5
        x=nrap
        y=zrapsum(nrap,itime)
        write(ifhi,'(2e13.5)')x,y/nevent
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      end 

      !------------------------
      subroutine uplot2
      !------------------------
      include 'aaa.h'
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!--------------------------'
      write(ifhi,'(a,i3)') '!   uplot2     '
      write(ifhi,'(a)')    '!--------------------------'
      print*,'nevent=',nevent
      do itime=1,31,2
       call getchar(itime,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -5 5'
       write(ifhi,'(a)')    'txt  "xaxis [c] "'              
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis dn/d[c] "'
       write(ifhi,'(a)')'array 2'
       do neta=-5,5
        x=neta
        y=zetasum(neta,itime)
        write(ifhi,'(2e13.5)')x,y/nevent
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      end 
      
      !------------------------
      subroutine uplot1
      !------------------------
      include 'aaa.h'
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!-----------------------'
      write(ifhi,'(a,i3)') '!   uplot1    '
      write(ifhi,'(a)')    '!-----------------------'
      np=0
      do neta=-4,4,2
      do ii=1,2
      do itime=1,31,2
       np=np+1
       call getchar(np,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -10 10'
       if(ii.eq.1)write(ifhi,'(a)')    'txt  "xaxis x (fm)"'
       if(ii.eq.2)write(ifhi,'(a)')    'txt  "xaxis y (fm)"'
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.1 0.9 "  [c]=',neta,'"'
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis ptl density "'
       write(ifhi,'(a)')'array 2'
       do ij=-10,10
        x=ij
        y=zeta(ii,ij,neta,itime)
        write(ifhi,'(2e13.5)')x,y/nevent/2.
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      enddo
      enddo
      end 


c---------------------------------------------------------------------------------
      subroutine checkdecay(mm1,ii1,iiz1,iret)
c---------------------------------------------------------------------------------

cinput mm1   : mass of  particle 
cinput ii1   : ID of  particle 
cinput iiz1  : $2\cdot I_3$ of particle 


      implicit none
      integer i1,iz1,ii1,iiz1,is,iret
      real*8 m1,mm1
      include 'comres.f'
      include 'options.f'
      integer strit

      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio

      iret=0
      i1=ii1
      iz1=iiz1
      m1=mm1
      
      is=iabs(strit(i1))

      if(iabs(i1).ge.minmes)then ! meson dec. 
         call anndexx(m1,i1,iz1,
     .        maxbrm ,minmes+1,maxmes,bmtype,branmes,iret)
      else if(is.eq.0)then   ! n*,d,d*
         call anndexx(m1,i1,iz1,
     .        maxbra,minnuc+1,maxdel,brtype,branres,iret)
      else if(is.eq.1)then   ! 
         call anndexx(m1,i1,iz1,
     .         maxbrs1,minlam+1,maxsig,bs1type,branbs1,iret)
      else if(is.eq.2)then
         call anndexx(m1,i1,iz1,
     .        maxbrs2,mincas+1,maxcas,bs2type,branbs2,iret)
      else
         write(ifmt,*)'make22(anndexx): s=',is,'not included'
         stop
      end if
      return
      end
      
      subroutine anndexx(m1,i1,iz1,
     &            maxbr,mini,maxi,btype,branch,iret)

      implicit none

      include 'comres.f'
      include 'comwid.f'
      include 'newpart.f'
      include 'options.f'

      real*8 pi,cc
      parameter(pi=3.1415927,cc=0.38937966)
      integer maxbr,mini,maxi,btype(4,0:maxbr)
      real*8 branch(0:maxbr,mini:maxi)
      integer icnt,iret
      integer i,i1,iz1    !,j
      real*8 m1,prob(0:100),sum
      real*8 mminit,fbrancx
      integer isoit
  
      do 3 i=0,maxbr
         if(isoit(btype(1,i))+isoit(btype(2,i))+isoit(btype(3,i))+
     &      isoit(btype(4,i)).lt.iabs(iz1).or.
     &        m1.lt.mminit(btype(1,i))+mminit(btype(2,i))
     &             +mminit(btype(3,i))+mminit(btype(4,i)) )then
            prob(i)=0.d0
         else
            prob(i)=fbrancx(i,iabs(i1),iz1,m1,branch(i,iabs(i1)),
     &           btype(1,i),btype(2,i),btype(3,i),btype(4,i))
         endif
 3    continue

      icnt=0
      call getbran(prob,0,100,sum,0,maxbr,i)

      if(i.gt.maxbr)then
         iret=1
         !write(ifmt,*)'anndexx(dec): no final state found for:',i1,m1,iz1
         !write(ifmt,*)'please check minimal masses: m1,m1min,m2min'
         !write(ifmt,*)'and iso3 of decaying particle'
         !write(ifmt,*)(prob(j),j=0,maxbr)
         !stop
      end if


      return
      end

c----------------------------------------------------------------------------------
                   subroutine f15outchxx(colldens)
c----------------------------------------------------------------------------------

      implicit none
      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'
      include 'outcom.f'
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer ninit,id,nupa(0:10)
      common/cnupa/nupa
      common /cninit/ninit
      integer nrap,ntau,i,ii
      integer nuutaumx,nrad,nin
      real*8 sigmatot,colldens,deluutau
      common /outco2/sigmatot
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      common /cf15/deluutau,nuutaumx
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      real*8 colldensdummy
      colldensdummy=colldens
      
      nin=ninit
ccc  if (bf15) return
      !if(ish.ge.2)write(ifch,*)'-----------'
      do i=1,nin
      call getind(tp0(i),tpz(i),tr0(i),trx(i),try(i),trz(i)
     . ,nrap,ntau,nrad,tityp(i),tiso3(i),id)
      !if(ish.ge.2)write(ifch,'(4x,i6,i5,3x,3i3)')tind(i),id,ntau,nrad,nrap
      !write(ifch,'(i6,$)')id
       !         istr=strit(tityp(i))
       !         ich = fchg(tiso3(i),tityp(i))
       !         
       !         write(15,*) tind(i),tr0(i),trx(i),try(i),trz(i),
       !     @                   tp0(i),tpx(i),tpy(i),tpz(i),tm(i),
       !     @                   tityp(i),tiso3(i),ich,tlcoll(i),
       !     @                   tcoll(i),istr,torigin(i)
      enddo
      !if(ish.ge.2)write(ifch,*)'     --------->'   
      !write(ifch,'(a,$)')'     --------->  '   
      if(nexit.le.10)nupa(nexit)=nupa(nexit)+1   
      do ii=1,nexit
       i=inew(ii)
       call getind(p0(i),pz(i),r0(i),rx(i),ry(i),rz(i)
     .  ,nrap,ntau,nrad,ityp(i),iso3(i),id)
       dcoll(nrap,ntau)=dcoll(nrap,ntau)+1
       rcoll(nrad,ntau)=rcoll(nrad,ntau)+1
       dcoll(nrap,nuutaumx+1)=dcoll(nrap,nuutaumx+1)+1
       rcoll(nrad,nuutaumx+1)=rcoll(nrad,nuutaumx+1)+1
      !if(ish.ge.2)write(ifch,'(4x,i6,i5,3x,3i3)')i,id,ntau,nrad,nrap
      !write(ifch,'(i6,$)')id
      enddo
      !write(ifch,'(a)')' '
      end

      subroutine xf15
      integer nupa(0:10)
      common/cnupa/nupa
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer nuutaumx
      real*8 deluutau,avpa,sum
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      common /cf15/deluutau,nuutaumx
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      if(ish.ge.2)write(ifch,*)'Initial particles(y,tau)'
      do j=-1,nuutaumx/2
      if(ish.ge.2)write(ifch,'(5x,11i5)')
     . (nint(dpart(i,j)/zevents),i=-5,5)
      enddo
      if(ish.ge.2)write(ifch,'(5x,11i5)')
     . (nint(dpart(i,nuutaumx+1)/zevents),i=-5,5)
      if(ish.ge.2)write(ifch,*)'ParticlesFromCollisions(y,tau)'
      do j=-1,nuutaumx+1
      if(ish.ge.2)write(ifch,'(5x,11i5)')
     . (nint(dcoll(i,j)/zevents),i=-5,5)
      enddo
      avpa=0
      sum=0
      do i=0,10
        sum=sum+nupa(i)
        avpa=avpa+nupa(i)*i
      enddo
      avpa=avpa/sum
      if(ish.ge.2)write(ifch,*)'nexit distribution:'
      if(ish.ge.2)write(ifch,'(5x,7i5,3x,a,f5.2)')
     . (nint(nupa(j)/zevents),j=0,6)
     . ,'mean=',avpa
      if(ish.ge.2)write(ifch,*)
     . 'ParticlesFromCollisions(y,tau)/mean_nexit'
      if(ish.ge.2)write(ifch,'(5x,11i5)')
     . (nint(dcoll(i,nuutaumx+1)/zevents/avpa),i=-5,5)
      if(ish.ge.2)write(ifch,*)
     . 'Collisions/Initial_particle(y,tau)'
      if(ish.ge.2)write(ifch,'(5x,11f5.2)')
     . (dcoll(i,nuutaumx+1)/avpa/dpart(i,nuutaumx+1),i=-5,5)
      !if(ish.ge.2)write(ifch,*)'Collisions(rad,tau)'
      !do j=-1,nuutaumx
      !if(ish.ge.2)write(ifch,'(16i4)')(rcoll(i,j),i=1,16)
      !enddo
      end

      subroutine getind(p0,pz,r0,rx,ry,rz,nrap,ntau,nrad,ityp,iso3,id)
      implicit none
      integer itab
      common/citab/itab
      real*8 pz,p0,r0,rx,ry,rz,tau,deluutau,rad
      integer nrap,nrapid,nuutaumx,ntau,nrad,ityp,iso3
     . ,id,idpdgg,pdgid,idtrafo
      common /cf15/deluutau,nuutaumx
      deluutau=1
      nuutaumx=20
      nrap=nrapid(pz,p0)
      ntau=-1
      tau=r0**2-rz**2
      if(tau.gt.0d0)then
       tau=sqrt(tau) 
       if(tau+deluutau/2.gt.(nuutaumx+1)*deluutau)then
        ntau=nuutaumx+1
       else
        ntau=(tau+deluutau/2)/deluutau
        ntau=min(ntau,nuutaumx+1)
       endif
      endif
      rad=-1
      rad=rx**2+ry**2
      nrad=1
      if(rad.gt.0d0)then
       rad=sqrt(rad)
       if(rad.gt.20)then
         nrad=20
       else
        nrad=1+rad
        nrad=min(nrad,20)
       endif
      endif
      idpdgg=pdgid(ityp,iso3)
      id=idpdgg
      if(itab.eq.1)id=idtrafo('pdg','nxs',idpdgg)
      end 

c------------------------------------------------------------------------
      subroutine printReaction(id1,id2,istra,kmin,kmax)
c------------------------------------------------------------------------
      !id1,id2: indx incoming
      !kk=kmin,kmax
      !kk=1: lstr computed if istra=1, may be used  via /clstr/
      !kk=2: full printout
      !kk=1 or 2: ity
      !-------------------------------------------------------
      !id of out-particles determined in subroutine anndex
      !branch i from call getbran(prob,0,100,sum,0,maxbr,i)
      !prob contains the branching ratios
      !-------------------------------------------------------
      include 'newpart.f'
      include 'outcom.f'
      include 'coms.f'
      integer ninit,nin,idpdgg,pdgid,idtrafo,id,ii,kk
      common /cninit/ninit
      logical lstrange,lstr
      common/clstr/lstr
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer itypart(nmax),iorpart(nmax),ity(2)
      common /city/itypart /corpart/iorpart
      ity(1)=itypart(id1)
      if(id2.ne.0)then
        ity(2)=itypart(id2)
      else
        ity(2)=0
      endif
      ipdg=1   !ipdg may  be set 0 or 1
      lstr=.true.
      if(istra.eq.1)lstr=.false.
      nin=ninit
      do kk=kmin,kmax
        if(kk.eq.2.and.lstr)write(ifch,'(a,$)')'    '   
        do i=1,nin
          idpdgg=pdgid(tityp(i),tiso3(i))
          if(ipdg.eq.1)then
          id=idpdgg
          else
          id=idtrafo('pdg','nxs',idpdgg)
          endif
          if(kk.eq.1.and.lstrange(id))lstr=.true.
          if(kk.eq.2.and.lstr)write(ifch,'(2i6,$)')id,ity(i)
        enddo
        if(kk.eq.2.and.lstr)write(ifch,'(a,$)')'  ----->  '   
        do ii=1,nexit
          i=inew(ii)
          idpdgg=pdgid(ityp(i),iso3(i))
          if(ipdg.eq.1)then
          id=idpdgg
          else
          id=idtrafo('pdg','nxs',idpdgg)
          endif
          if(nin.eq.1)then
            itypart(i)=ity(1)
            iorpart(i)=-999
          else
            iorpart(i)=0
            if(nexit.eq.2)then
            itypart(i)=ity(ii)
            else
            itypart(i)=62  
            endif
          endif
          if(kk.eq.1.and.lstrange(id))lstr=.true.
          if(kk.eq.2.and.lstr)write(ifch,'(2i6,$)')id,itypart(i)
        enddo
        if(kk.eq.2.and.lstr)write(ifch,'(a)')'   '   
      enddo
      end

c----------------------------------------------------------------
      logical function lstrange(id)
c----------------------------------------------------------------
      include 'coms.f'
      integer id,ia
      logical lstr
      lstr=.false.
      ia=abs(id)
      if(ia.eq.2130.or.ia.eq.1230.or.ia.eq.1130.or.ia.eq.2230
     .         .or.ia.eq.1131.or.ia.eq.1231.or.ia.eq.2231)
     . lstr=.true.
      lstrange=lstr
      end

cc------------------------------------------------------------------------
      subroutine getchar(np,ch)
cc------------------------------------------------------------------------
      character*4 ch
      ch='    '
      if(np.le.9)then
       write(ch,'(a,i1)')np
      elseif(np.le.99)then
       write(ch,'(a,i2)')np
      elseif(np.le.999)then
       write(ch,'(a,i3)')np
      else
       ch='????'
      endif 
      end

cc---------------------------------------------------------------------
          integer function nrapid(p3,p4)
cc---------------------------------------------------------------------
      real*8 p3,p4
      if(p4-p3.le.0.)then
      nrapid=1000
      elseif(p4+p3.le.0.)then
      nrapid=-1000
      else
      rap=0.5*log((p4+p3)/(p4-p3))
      nrapid=-19+int(rap+19.5)
      endif
      nrapid=max(nrapid,-10)
      nrapid=min(nrapid, 10)
      end

cc---------------------------------------------------------------------
          subroutine writeParticleTestFile(nn)
cc---------------------------------------------------------------------
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list
      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)
      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
      open(unit=89,file='uuu.data',status='old',access='append')
      write(89,*)'ptl'
      write(89,*)
     .   istptl(nn)
     .  ,ityptl(nn)
     .  ,idptl(nn)
     .  ,xorptl(1,nn)
     .  ,xorptl(2,nn)
     .  ,xorptl(3,nn)
     .  ,xorptl(4,nn)
     .  ,pptl(1,nn)
     .  ,pptl(2,nn)
     .  ,pptl(3,nn)
     .  ,pptl(4,nn)
     .  ,pptl(5,nn)
      close(89)
      end

      
      subroutine openNewTestFile
      data ncnttf/0/
      save ncnttf
      ncnttf=ncnttf+1
      if(ncnttf.eq.1)then
      open(unit=89,file='uuu.data',status='new')
      write(89,*)'test'
      close(89)
      endif
      end
      
      subroutine writeEndEventTestFile
      open(unit=89,file='uuu.data',status='old',access='append')
      write(89,*)'end'
      close(89)
      end

     
c###################################################################################
c###################################################################################
c#######
c#######
c#######       random number stuff
c#######
c#######
c###################################################################################
c###################################################################################


      double precision function ranff(idummy)
      jdummy=idummy
      ranff = rangen()
      return
      end
      subroutine sseed(ranseed)
      integer ranseed
      dummyseed=ranseed
      return
      end


c###################################################################################
c###################################################################################
c#######
c#######
c#######       tables.dat stuff
c#######
c#######
c###################################################################################
c###################################################################################

c----------------------------------------------------------
      subroutine loadwtab (io)
c----------------------------------------------------------
c   output   : information in common-block comwid.f
c
c     load the tabulated branching ratios from disk
c 
c----------------------------------------------------------

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer ios, nsp, io, ver
      character*50 deftab
      character*8 defexe
      logical b
      character*500 edir
      common/cdir/edir

c set the defaultname of the file, containing the table
      parameter (deftab='tables.dat ', defexe='uqmd.exe')

      b=io.eq.1

c set the name of the table

      tabname=edir(1:index(edir,' ')-1)
     .//deftab(1:index(deftab,' '))

      if(b)write (6,*) 'Looking for the tabulated decay width...'
      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
     .      status='old')
      if (ios.eq.0) then
         if(b)write (6,'(a)') ' (info) reading tables.dat'
         read (75,*,err=99) ver, nsp, tabx, fbtaby, pbtaby,
     .   fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .	 tabxnd, frrtaby
         if(b)write (6,*) 'O.K.'
         wtabflg=3
         if (ver.eq.tabver) then
            if(b)write (6,*) 'tabver=',ver,'  O.K.'
         else
            write (6,*) '(info) wrong tables.dat!'
            write (6,*) '(info) tabver should be',tabver,'
     .        ,instead of',ver
            stop
         endif
         if (nsp.eq.widnsp) then
            if(b)write (6,*) 'widnsp=',nsp,'  O.K.'
         else
            write (6,*) '(info) wrong table!'
            write (6,*) '(info) widnsp should be',widnsp
     .       ,', instead of',nsp
   	    stop
	 endif
         close (unit=75, status='keep')
      else 
        write (6,*) 'No file ',tabname(1:index(tabname,' ')-1)
        close (unit=75, status='delete')
        call mkwtab
      endif

      return

  99  continue
      write (6,*) '(info) read error, reading from '
     .,tabname(1:index(tabname,' ')-1)
      stop
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine savewtab
c
c     Revision : 1.0
c
c     save the tabulated branching ratios to disk
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer ios

      write (6,*) 'Writing new table...'
c try to generate a new file 
      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
     .      status='new')
c if it succedds ...
      if (ios.eq.0) then
c write the tables into the file
         write (75,*) tabver, widnsp, tabx, fbtaby, 
     .        pbtaby, fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .	    tabxnd, frrtaby
c otherwise complain
      else
         write (6,*) 'Error: ',tabname(1:index(tabname,' ')-1)
     .    ,'exists!'
      endif
c close the file
      close (unit=75, status='keep')
      
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkwtab
c
c     Revision : 1.0
c
coutput   : information in common-block comwid.f
c
c     tabulate the mass dependent branching ratios 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      real*8 fwidth,m,first,last,delta,abl0,abln,mir,mminit,fbrancx
      real*8 massit,bran,smass,bwnorm,fppfit
      integer i,bchan,itp,isoit,cmin,cmax,i1,i2,i3,i4,ii1

      write (6,*) 'Generating table...'
c this indicates, that all tables are still empty
      wtabflg=0

c high precision splines from mintab to maxtab1
c lower precicision between maxtab1 and maxtab2

c now fill the x-values
c start with 'mintab'
      first=mintab
c 66 % of all fixpoints between mintab and maxtab1
c calculate the steps
      delta=(maxtab1-mintab)/((widnsp-1d0)*2d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab1>mintab in comwid'
         stop
      endif
c store the values into 'tabx'
      do 10 i=1,int(widnsp*2./3.)
         m=first+(i-1)*delta
         tabx(i)=m
 10   continue
c 33 % of all fixpoints with larger delta between maxtab1 and maxtab2
      delta=(maxtab2-maxtab1)/((widnsp-1d0)*1d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab2>maxtab1 in comwid'
         stop
      endif
c store the values into 'tabx'
	do 11 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabx(i)=m
 11   continue

c now fill the y-values of the full branching ratios

c these are the first derivatives at the first an the last point
c of the interpolating function. a value greater than 1E30 signals the 
c 'spline' routine to set the boundary condition for a natural spline
c with zero second derivative 
      abl0=2D30
      abln=2D30

c loop over all baryons
      do 20 itp=minbar,maxbar
c loop over all x-values
         do 21 i=1,widnsp
c store the values ...
            fbtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 21      continue
c calculate the second derivate and store it in 'fbtaby(,,2)'
         call spline (tabx(1),fbtaby(1,itp,1),widnsp,abl0,abln,
     .                fbtaby(1,itp,2))
 20   continue
      write (6,*) '(1/7) ready.'

c loop over all mesons
      do 30 itp=minmes,maxmes
c loop over all x-values
         do 31 i=1,widnsp
c store the values ...
            fmtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 31      continue
c calculate the second derivate and store it in 'fmtaby(,,2)'
         call spline (tabx(1),fmtaby(1,itp,1),widnsp,abl0,abln,
     .                fmtaby(1,itp,2))
 30   continue
      write (6,*) '(2/7) ready.'

c the flag indicates, that now all full widths are tabulated 
      wtabflg=1

c now fill the y-values of the partial branching ratios

c loop over all baryons
      do 40 itp=minbar,maxbar
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels
         call brange (itp, cmin, cmax)
c check, if there are any decay channels         
         if (cmax.gt.0) then
c loop over all decay channels
            do 41 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type (itp,bchan,bran,i1,i2,i3,i4)
c check, if decay is allowed

               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values               
                  do 42 i=1,widnsp
c store the values
                     pbtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 42               continue
c calculate the second derivate and store it in 'pbtaby(,2,,)'
                  call spline (tabx(1),pbtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pbtaby(1,2,itp,bchan))
               end if
 41         continue
         end if
 40   continue
      write (6,*) '(3/7) ready.'

c loop over all mesons
      do 50 itp=minmes,maxmes
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels 
         call brange (itp, cmin, cmax)
c check, if there are any decay channels         
         if (cmax.gt.0) then
            do 51 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type(itp,bchan,bran,i1,i2,i3,i4)
c!!!
               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values               
                  do 52 i=1,widnsp
                     pmtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 52               continue
c calculate the second derivate and store it in 'pmtaby(,2,,)'
                  call spline (tabx(1),pmtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pmtaby(1,2,itp,bchan))
               end if
 51         continue
         end if
 50   continue

      write (6,*) '(4/7) ready.'


c calculate the norm integral of the Breit-Wigner functions
c   with mass dependent widths

c..baryons
	do 60 i=minbar,maxbar
	   bwbarnorm(i)=bwnorm(i)
60	continue
      write (6,*) '(5/7) ready.'

c.. mesons
	do 61 i=minmes,maxmes
	   bwmesnorm(i)=bwnorm(i)
61	continue
      write (6,*) '(6/7) ready.'

c now all branching ratios and BW-integrals are tabulated 
      wtabflg=2

ce tabulate fppfit
c fill the x-values
c range of tabulated cross sections
      first=2d0*massit(nucleon)+massit(pimeson)
	last=maxtab1
c calculate the steps
c the energies are weighted quadratically
      delta=(last-first)/((widnsp-1)*2./3.)**2
c store the values into 'tabx'
c 66 % of all fixpoints between mintab and maxtab1
      do 69 i=1,int(widnsp*2./3.)
         m=first+(i-1)**2*delta
         tabxnd(i)=m
 69   continue
c 33 % of all fixpoints with larger, constant delta between maxtab1 and maxtab2
	delta=(maxtab2-last)/((widnsp-1)*1./3.)
	do 70 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabxnd(i)=m
 70   continue


c.. all pp-exit channels
c loop over first out-particle N & D
	do 81 ii1=1,2
	  if(ii1.eq.1)i1=minnuc
	  if(ii1.eq.2)i1=mindel
c loop over second out-particle N(1440)..maxdel
	  do 82 i2=minnuc+1,maxdel
c loop over all x-values
          do 83 i=1,widnsp
c store the values ...
	      frrtaby(i,1,ii1,i2)=fppfit(99,tabxnd(i),i1,i2)
83	    continue
c calculate the second derivate and store it in 'frrtaby(,,2)'
          call spline (tabxnd(1),frrtaby(1,1,ii1,i2),widnsp,abl0,abln,
     .                frrtaby(1,2,ii1,i2))
82	  continue
81	continue


      write (6,*) '(7/7) ready.'

c pp cross sections are now tabulated
	wtabflg=3

c save the table on disk
      call savewtab

      return
      end


c######################################################################
c######################################################################
c#########
c#########                    utilities
c#########
c######################################################################
c######################################################################

     
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splint (xa,ya,y2a,n,x)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c     Date     : 03/07/96
c     Revision : 1.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop '\n\n STOP in splint: bad xa input\n\n'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .            (b*b*b-b)*y2a(khi))*(h*h)/6d0
      splint=y

      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splintth (xa,ya,y2a,n,x,th)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c                modified my H. Weber
c     Date     : 03/07/96
c     Revision : 1.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     split routine with nice threshold behaviour for cross sections
c

      implicit none

      include 'comres.f'
      include 'comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h,th

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop '\n\n STOP in splint: bad xa input \n\n'
      if (xa(khi).lt.(th+2*h)) then
c linear approximation close to threshold (within 2h)
         splintth=ya(khi)*(x-th)/(xa(khi)-th)
      else
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .        (b*b*b-b)*y2a(khi))*(h*h)/6d0
         splintth=y
      endif
      
      return
      end



