      module f2_grid
                implicit none
                real*8, public :: fmap_MSTW_F2(120,120)
                real*8, public :: fmap_MSTW_FL(120,120)

        contains

      subroutine READ_GRID(fmap_MSTW_F2,fmap_MSTW_FL)

                real*8 :: fmap_MSTW_F2(120,120)
                real*8 :: fmap_MSTW_FL(120,120)
        real*8 :: rx,rx_min,drx        
        real*8 :: rmu2,rmu2_min,drmu2        
        real*8 :: alog10_x,alog10_mu2,F2,FL
        integer :: i,j,irx,irmu2,nrx,nrmu2,nmax


      open(unit=11,file='F2_MSTW_ij.dat',
     &     status='unknown')

      open(unit=12,file='FL_MSTW_ij.dat',
     &     status='unknown') 

c     =================================================================

      While(i<120) do

      read(11,*) i, j, fmap_mstw_f2(i,j)
      read(12,*) i, j, fmap_mstw_fl(i,j)

      enddo


      close(unit=11)

      end subroutine READ_GRID

      end module f2_grid 
