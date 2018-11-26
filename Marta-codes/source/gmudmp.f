      subroutine gmudmp

*************************************************************************
*                                                                       *
*  Subroutine to dump all the input parameters                          *
*                                                                       *
*************************************************************************

      implicit none

      integer iend,ndim,npoin,nprin,ntreat,itmx,ncvg
      integer ilepton,proc,ngen
      double precision inp1,inp2,ptdiff_min,ptdiff_max,y_min,y_max
      double precision q1t_min,q1t_max,q2t_min,q2t_max
      double precision pt_min,pt_max
      double precision amx_min,amx_max,amy_min,amy_max

      common/vggen/iend,ndim,npoin,nprin,ntreat,itmx,ncvg
      common/cuts/inp1,inp2,
     +     ptdiff_min,ptdiff_max,
     +     y_min,y_max,
     +     pt_min,pt_max,
     +     q1t_min,q1t_max,q2t_min,q2t_max,
     +     amx_min,amx_max,amy_min,amy_max
      common/settings/ilepton,proc,ngen

      print *,'*** Parameters dump ************************************'
      print 1000,'Production mode',iend
      print 1000,'Printout mode',nprin
      print 1000,'Number of fct calls in Vegas',ncvg
      print 1000,'Number of Vegas iterations',itmx
      print 1000,'Use Treat ?',ntreat
      print *,'********************************************************'
      print 1000,'Number of events to generate',ngen
      print 1000,'Leptons pair',ilepton
      print 1000,'Process type',proc
      if(proc.eq.1) print 1003,' Elastic'
      if(proc.eq.2) print 1003,' Elastic-Inelastic'
      if(proc.eq.3) print 1003,' Inelastic-Elastic'
      if(proc.eq.4) print 1003,' Inelastic-Inelastic'
      print 1001,'First incoming protons pz',inp1
      print 1001,'Second incoming protons pz',inp2
      print 1002,'Leptons pT difference',ptdiff_min,ptdiff_max
      print 1002,'Leptons rapididy',y_min,y_max
      print 1002,'Leptons pT',pt_min,pt_max
      print 1002,'Q1T',q1t_min,q1t_max
      print 1002,'Q2T',q2t_min,q2t_max
      print 1002,'Proton 1 remnant mass',amx_min,amx_max
      print 1002,'Proton 2 remnant mass',amy_min,amy_max
      print *,'********************************************************'

 1000 format('* ',a37,':',i15,' *')
 1001 format('* ',a37,':',f15.2,' *')
 1002 format('* ',a37,': [',f5.1,',',f6.1,'] *')
 1003 format('* -------->',a44,' *')
      return
      end
