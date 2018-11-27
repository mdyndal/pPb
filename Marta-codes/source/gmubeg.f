      subroutine gmubeg

      implicit none
      
      integer proc,ilepton
      integer iend,ndim,npoin,nprin,ntreat,aitmx,ncvg
      common/vggen/iend,ndim,npoin,nprin,ntreat,aitmx,ncvg
c      data ndim,npoin/7,100/
      data npoin/100/
      double precision sigma
      external sigma
      double precision osigma,oerrsigma,pl(4,9)
c
c   makes default assignments for Vegas

      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds,ngen
      double precision alph
      double precision avgi, sd, chi2a
      double precision xl,xu,acc,si,swgt,schi,xi

      double precision fmap_mstw_F2,fmap_mstw_FL

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds
      common/settings/ilepton,proc,ngen
 
ccc      common/transfer1/fmap_mstw_F2(120,120),
ccc     >                fmap_mstw_FL(120,120)


      common/output/osigma,oerrsigma,pl

c      data acc/-1./,itmx/5/
      data ncall/5000/,itmx/5/,nprn/5/,acc/-1./,
     1     xl/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./,
     2     xu/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./,
     3     alph/1.5/,ndmx/50/,mds/1/,ndev/6/,
     4     ndo/1/,xi/550*1./,it/0/,si,swgt,schi/3*0./
     
cc       save/transfer1/

      if(proc.eq.1)then
         ndim=8
      elseif((proc.eq.2).or.(proc.eq.3))then
         ndim=9
      elseif(proc.eq.4)then
         ndim=10
      else
         print *,'GMUBEG: ERROR check the process number!'
         stop
      endif
      
c-----Probing of the phase space
c      nprn=-1
      nprn=nprin
      ncall=1000000
      itmx=10
c      itmx=2
      print 1000, ndim,ncall,itmx
      call vegas(ndim,sigma,avgi,sd,chi2a)
      print *,'Total cross-section =',avgi,' +/-',sd,'nb'
      
c-----Vegas integration
      nprn=nprin
      ncall=ncvg
      itmx=aitmx
      print 1001, ndim,ncall,itmx
      call vegas1(ndim,sigma,avgi,sd,chi2a)
      print *,'Total cross-section =',avgi,' +/-',sd,'nb'
      osigma=avgi
      oerrsigma=sd
      
C-----Preparation for events generation if needed
      if (iend.lt.2) then
         print *,
     &        ' GMUBEG : Program stops without SETGEN call and'//
     &        ' events generation , IEND < 2 '
         stop
      else if (nprin.gt.1) then
         print *,'GMUBEG : ===> SETGEN is operative... '
      endif
      
      if (iend.lt.3) then
         print *,
     &        ' GMUBEG : Program stops without events generation,'//
     &        ' IEND < 3'
         stop
      endif
c      call setgen(sigma,ndim,npoin,nprin,ntreat)
      call setgen(sigma,ndim,npoin,1,ntreat)

 1000 format('*** Vegas warm-up ******************************'/
     +       '* Number of dimensions to integrate:',i10,' *'/
     +       '* Number of function calls         :',i10,' *'/
     +       '* Maximal number of iterations     :',i10,' *'/
     +       '************************************************'
     +)
 1001 format('*** Vegas integration **************************'/
     +       '* Number of dimensions to integrate:',i10,' *'/
     +       '* Number of function calls         :',i10,' *'/
     +       '* Maximal number of iterations     :',i10,' *'/
     +       '************************************************'
     +)

      return
      end
