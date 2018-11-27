      subroutine gmufil

**********************************************************
*                                                        *
*  Fill the output file with the kinematic quantities    *
*  for the in- and outgoing particles. The convention    *
*  is the following :                                    *
*   (1) = incoming proton 1                              *
*   (2) = incoming proton 2                              *
*   (3) = outgoing proton (or proton remnant) 1          *
*   (4) = central two-photon system                      *
*   (5) = outgoing proton (or proton remnant) 2          *
*   (6) = outgoing lepton 1                              *
*   (7) = outgoing lepton 2                              *
*                                                        *
**********************************************************

      implicit none

      integer nlimax
      double precision pi2
      double precision pt1x,pt1y,pt2x,pt2y,y1,y2
      double precision am_l
      double precision pt1,pt2,et1,et2
      double precision sigma,errsigma
      double precision px_0,px_x,px_y,px_z,py_0,py_x,py_y,py_z
      double precision ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z

      parameter(pi2=2.*3.14159265)
      parameter(nlimax=7) ! maximal number of particles in the event

      common/kinematics/pt1x,pt1y,y1,pt2x,pt2y,y2,am_l,
     +     ak10,ak1x,ak1y,ak1z,ak20,ak2x,ak2y,ak2z,
     +     px_0,px_x,px_y,px_z,py_0,py_x,py_y,py_z
      common/output/sigma,errsigma,pl

*---- Convention used here : (E,px,py,pz)
      double precision pl(4,nlimax)
      
      pt1=dsqrt(pt1x**2+pt1y**2)
      et1=dsqrt(pt1**2+am_l**2)
      pt2=dsqrt(pt2x**2+pt2y**2)
      et2=dsqrt(pt2**2+am_l**2)

*---- First incoming proton
      pl(1,1)=ak10
      pl(2,1)=ak1x
      pl(3,1)=ak1y
      pl(4,1)=ak1z
*---- Second incoming proton
      pl(1,2)=ak20
      pl(2,2)=ak2x
      pl(3,2)=ak2y
      pl(4,2)=ak2z
*---- First outgoing proton (or remnant)
      pl(1,3)=px_0
      pl(2,3)=px_x
      pl(3,3)=px_y
      pl(4,3)=px_z
*---- Central two-photon system
      pl(1,4)=0.d0
      pl(2,4)=0.d0
      pl(3,4)=0.d0
      pl(4,4)=0.d0
*---- Second outgoing proton (or remnant)
      pl(1,5)=py_0
      pl(2,5)=py_x
      pl(3,5)=py_y
      pl(4,5)=py_z
*---- First outgoing lepton
      pl(1,6)=et1*dcosh(y1)
      pl(2,6)=pt1x
      pl(3,6)=pt1y
      pl(4,6)=et1*dsinh(y1)
*---- Second outgoing lepton
      pl(1,7)=et2*dcosh(y2)
      pl(2,7)=pt2x
      pl(3,7)=pt2y
      pl(4,7)=et2*dsinh(y2)
      
      end
