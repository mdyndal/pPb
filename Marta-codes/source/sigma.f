c     =================================================================
c
c     Integration function
c     -----------------------------------------------------------------  
      function sigma(x)
      implicit none
      double precision sigma
      double precision x(11),wgt,jac,tmp
      double precision inp1,inp2
      double precision q1t_min, q1t_max, q2t_min, q2t_max
      double precision y_min, y_max
      double precision pt_min, pt_max
      double precision ptdiff_min, ptdiff_max
      double precision q1t,q2t,phiq1t,phiq2t,ptdiff,phiptdiff,pi 
      double precision lq1t_min,lq1t_max,lq1t,lq2t_min,lq2t_max,lq2t
      double precision amx_min,amx_max,amx
      double precision amy_min,amy_max,amy
      double precision am_p
      double precision fmap_mstw_F2,fmap_mstw_FL
      double precision integrand
c     -----------------------------------------------------------------
      integer itmx,nprn,it,ndo
      integer ndim,ncvg,igraph,npoin,nprin,ntreat
      double precision si,si2,swgt,schi,xi,scalls,d,di
      double precision pt1x,pt1y,y1,pt2x,pt2y,y2,am_l
      double precision q10,q1x,q1y,q1z,q20,q2x,q2y,q2z
      double precision ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z
      common/vgb2/ndo,it,si,si2,swgt,schi,xi(50,11),scalls,
     +     d(50,11),di(50,11)
      common/vegpar/ndim,ncvg,itmx,nprn,igraph,npoin,nprin,ntreat
      common/kinematics/pt1x,pt1y,y1,pt2x,pt2y,y2,am_l,
     +     ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z,
     +     q10,q1x,q1y,q1z,q20,q2x,q2y,q2z
      common/cuts/inp1,inp2,
     +     ptdiff_min,ptdiff_max,
     +     y_min,y_max,
     +     pt_min,pt_max,
     +     q1t_min,q1t_max,q2t_min,q2t_max,
     +     amx_min,amx_max,amy_min,amy_max
c      common/vgb3/wgt
      integer proc,ilepton,ngen
      common/settings/ilepton,proc,ngen

ccc      common/transfer1/fmap_mstw_F2(120,120),
ccc     2                fmap_mstw_FL(120,120)

ccc        save/transfer1/
c     =================================================================
c     PHASE SPACE (range of integration)
c     =================================================================
c
      am_p = 0.93827203d0

      pi = 4.d0*datan(1.d0)

c     -----------------------------------------------------------------

c      q1t = q1t_min + (q1t_max-q1t_min)*x(1)
c      q2t = q2t_min + (q2t_max-q2t_min)*x(2)
      lq1t_max = dlog(q1t_max)
      lq1t_min = -10.d0         ! FIXME sufficiently low ??
      lq2t_max = dlog(q2t_max)
      lq2t_min = -10.d0
      lq1t = lq1t_min + (lq1t_max-lq1t_min)*x(1)
      lq2t = lq2t_min + (lq2t_max-lq2t_min)*x(2)
      q1t = dexp(lq1t)
      q2t = dexp(lq2t)

      phiq1t = 2.d0*pi*x(3)
      phiq2t = 2.d0*pi*x(4)
      y1 = y_min+(y_max-y_min)*x(5)
      y2 = y_min+(y_max-y_min)*x(6)
      ptdiff = ptdiff_min + (ptdiff_max-ptdiff_min)*x(7)
      phiptdiff = 2.d0*pi*x(8)

c     =================================================================
c     jacobian: x(n) ----> phase space
c     =================================================================
c      jac = (q1t_max-q1t_min) * (q2t_max-q2t_min) 
c     2    * (y_max-y_min)**2 * (2.d0*pi)**3 * (ptdiff_max-ptdiff_min)
      jac = (lq1t_max-lq1t_min) * q1t * (lq2t_max-lq2t_min) * q2t
     2    * (y_max-y_min)**2 * (2.d0*pi)**3 * (ptdiff_max-ptdiff_min)

      if(proc.ge.2) then
         amx = amx_min + (amx_max - amx_min)*x(9)
         jac = jac * (amx_max - amx_min) * 2.d0 * amx 
         amy = am_p
         if(proc.eq.2) then
		amy = amx
		amx = am_p
	 endif
         if(proc.eq.4) then
            amy = amy_min + (amy_max - amy_min)*x(10)
            jac = jac * (amy_max - amy_min) * 2.d0 * amy
         endif
      endif

c     =================================================================
c     calling subroutine with aintegrand
c     =================================================================

      tmp=wgt*jac/itmx

c      call INCqqbar(integrand,proc,q1t,q2t,phiq1t,phiq2t,y1,y2,
c     2              ptdiff,phiptdiff,amx,amy,tmp)

      call INCqqbar(integrand,proc,q1t,q2t,phiq1t,phiq2t,y1,y2,
     2              ptdiff,phiptdiff,amx,amy)



c     =================================================================
c     sigma total
c     =================================================================
      sigma = jac*integrand

c     =================================================================
      return
      end
