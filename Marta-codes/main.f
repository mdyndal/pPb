**********************************************************
*       MCGEN ver. 0.1 - KRAKOW - LOUVAIN-LA-NEUVE       *
* Program for analyzing u+u- pair in pp->pp+u+u- process *
* Created 28 October 2013                                *
* Authors: Antoni.Szczurek@ifj.edu.pl                    *
*          Marta.Luszczak@ur.edu.pl                      *
*          Wolfgang.Schafer@ifj.edu.pl                   *
*          Gustavo.Dasilveira@uclouvain.be               *
*          Laurent.Forthomme@uclouvain.be                *
*                                                        *
**********************************************************
      program main

      integer nev,nprt,file
      integer ilepton,proc,ngen
      double precision pl(4,9)
      double precision sigma,errsigma
      double precision fmap_mstw_F2,fmap_mstw_FL  

      parameter (nev=10000)             ! number of events NEV
      parameter (nprt=nev/10)           ! printing period Nprt

      common/output/sigma,errsigma,pl
      common/settings/ilepton,proc,ngen

      common/transfer1/fmap_mstw_F2(120,120),
     2                fmap_mstw_FL(120,120)

      save /transfer1/

      file=21

      call gmuini
      call gmucha
      call gmubeg
c      open(unit=file,file='events.dat',status='unknown')
c      do iev=1,ngen
c         call gmugna
c         print *,iev,'passed gmugna'
c         call gmufil
c         print *,iev,'passed gmufil'
c         do ipart=1,9
c            print *,'Particle',ipart
c            do ind=1,4
c               print *,'p(',ind,')=',pl(ind,ipart)
c            enddo
c         enddo
c         if(mod(iev,nprt).eq.0) print *,' Event nr = ',iev
c      enddo
c      close(unit=file)

      end program main
      
      
      
      
