      subroutine gmucha

*************************************************************************
*                                                                       *
*  Subroutine to fetch all the input parameters from an external        *
*  card. This input card is fetched at the argument to this program     *
*  if one is specified, or in a default file placed in this directory,  *
*   'input.card'                                                        *
*                                                                       *
*************************************************************************

      implicit none

      integer iend,ndim,npoin,nprin,ntreat,itmx,ncvg
      double precision inp1,inp2
      double precision q1t_min,q1t_max,q2t_min,q2t_max
      double precision y_min,y_max
      double precision ptdiff_min,ptdiff_max
      double precision amx_min,amx_max,amy_min,amy_max
      double precision pt_min,pt_max
      character(len=32) file
      integer i,lun,maxln,ngen
      integer ilepton,proc
      character(len=8) key
      character(len=10) value
      logical fexst

      common/vggen/iend,ndim,npoin,nprin,ntreat,itmx,ncvg
      common/cuts/inp1,inp2,
     +     ptdiff_min,ptdiff_max,y_min,y_max,
     +     pt_min,pt_max,
     +     q1t_min,q1t_max,q2t_min,q2t_max,
     +     amx_min,amx_max,amy_min,amy_max
      common/settings/ilepton,proc,ngen

      lun=15
      maxln=25

      if (iargc().gt.0) then
         call getarg(1,file)
      else
         file='input.card'
      endif

*---- Makes sure the file exists      
      inquire(file=file,exist=fexst)
      if (fexst.eqv..false.) then
         print *,'GMUCHA: ERROR! Input card does not exist!'
         stop
      endif

*---- Read the parameter card using key/value pairs
      open(lun,file=file,status='old')
      do i=1,maxln
         read(lun,1000,end=10) key,value
c         print *,'=',trim(key),'=',value
         if (trim(key).eq."IEND")   read(value,*) iend
         if (trim(key).eq."NGEN")   read(value,*) ngen
         if (trim(key).eq."PROC")   read(value,*) proc
         if (trim(key).eq."NPRIN")  read(value,*) nprin
         if (trim(key).eq."NCVG")   read(value,*) ncvg
         if (trim(key).eq."NTREAT") read(value,*) ntreat
         if (trim(key).eq."ITMX")   read(value,*) itmx
         if (trim(key).eq."ILEPTON")read(value,*) ilepton
         if (trim(key).eq."INP1")   read(value,*) inp1
         if (trim(key).eq."INP2")   read(value,*) inp2
         if (trim(key).eq."Q1TMIN") read(value,*) q1t_min
         if (trim(key).eq."Q1TMAX") read(value,*) q1t_max
         if (trim(key).eq."Q2TMIN") read(value,*) q2t_min
         if (trim(key).eq."Q2TMAX") read(value,*) q2t_max
         if (trim(key).eq."PTMIN")  read(value,*) pt_min
         if (trim(key).eq."PTMAX")  read(value,*) pt_max
         if (trim(key).eq."PTDFMIN")read(value,*) ptdiff_min
         if (trim(key).eq."PTDFMAX")read(value,*) ptdiff_max
         if (trim(key).eq."YMIN")   read(value,*) y_min
         if (trim(key).eq."YMAX")   read(value,*) y_max
         if (trim(key).eq."MXMIN")  read(value,*) amx_min
         if (trim(key).eq."MXMAX")  read(value,*) amx_max
         if (trim(key).eq."MYMIN")  read(value,*) amy_min
         if (trim(key).eq."MYMAX")  read(value,*) amy_max
      enddo
 10   continue
      close(lun)

c     Checks for the q1t/q2t integration boundaries
c---- FIXME FIXME move me somewhere else !
      if (((proc.eq.3).or.(proc.eq.4)).and.(q1t_max.lt.50.0)) then
         q1t_max = 50.0
      endif
      if (((proc.eq.2).or.(proc.eq.4)).and.(q2t_max.lt.50.0)) then
         q2t_max = 50.0
      endif

      call gmudmp
*     
 1000 format(a8,a10)
      return
      end
