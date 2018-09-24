!
!                                      **************************************
!                                      *     module pgnc                    *
!                                      *     R. J. Purser                   *
!                                      *     NOAA/NCEP/EMC September  2018  *
!                                      *     jim.purser@noaa.gov            *
!                                      **************************************
!
! A suite of routines for generating the coefficients needed to define the
! gnomonic cubed-sphere transformations and their derivatives.
!
! DIRECT DEPENDENCIES
! Libraries[their Modules]: pmat[pmat, pmat4, pmat5]
! Additional Modules      : pietc, pkind
!
!=============================================================================
module pgnc
!=============================================================================
! Generic smooth cubed-sphere gnomonic grid transformations that include
! equidistant, edge-uniform, equiangular and others in this continuous family,
! are accommodated with a parameter setting, n=0, and a non-negative real
! parameter, alpha, that specifies which member of this continuous family
! is intended. Alpha=0 implies it is the equidistant version; alpha=0.5 implies
! it is the edge-uniform grid presently used in the FV3; alpha=1.0 implies
! that it is the equiangular version. For alpha<1, the grid spacing along
! a tile-median gets finer towards the edge; alpha>1 is also allowed and
! causes the grid spacing to get coarser towards the edge.
! Alternatively, we can use integer n>0 to indicate that the intended 
! grid is one of the Mobius net kind, where this n specifies the intended
! degree of continuity of the transformation across the edges of the zones
! runiing along the edges of the cube where the Mobius net attribute of the
! grid pertains (three-way inetrsections among the three grid line families).
! The Mobius net grids are not perfectly smooth at these transitions, but
! the mappring and its first n derivatives are continuous there. In these
! cases, the angular halfwidth of each Mobius zone (measured in degrees along
! the perpendicular tile median) is the re-purposed real parameter, alpha.
! The degree of continuity, n, can go from 1 up to nn=8.
!
! Each of the cubic grids defined as above can also be further transformed
! by a Schmidt conformal mapping transformaton that enhances the resolution
! at the center of tile-1 of the protoframe by a factor, S, when S>1, or
! degrades the resolution there when S<1. (S=1 corresponds to no enhancement.)
! The Schmidt transformation is optional and, if it is not initialized, is
! omitted by default.
!
! For each of these, there is complete freedom to choose how the indices of the
! six cubic panels or "tiles" are arranged and their local coordinates
! oriented relative to a standard cubic "protoframe", whose own arrangement
! of panels, when opened up ("developed") into a connected planar map, looks
! like the sketch:
!
!  +-------+
!  |       |
!  |   6   |
!  !  -->  |
!  +-------+-------+-------+-------+
!  |       |       |       |       |
!  |   2   |   3   |   4   |   5   |
!  |  -->  |  -->  |  -->  |  -->  |
!  +-------+-------+-------+-------+
!  |       |
!  |   1   |
!  |  -->  |
!  +-------+
!
! Here, the short arrow indicate the direction of the local "x" coordinate
! of the right-handed pair, (x,y) of local map coordinates. Generally, it is
! assumed that the local coordinates have their origins at the tile centers
! and both x and y extend between -1 and +1 from edge to edge. The particular
! arrangement of the user's choice of panel indexing is specified by a 
! permutation of the six integers 1--6 that label the panels, so that the
! permutation, [6,1,2,4,5,3] puts user-defined label "6" in the position of
! the above protoframe's tile-1, and so on. Each permutation defined in this
! way can be lexicographically ordered, and therefore indexed from 0 to 719
! (the number of possible permuations of six objects is 6! = 720). The 
! index corresponding to [1,2,3,4,5,6] is 0; that corresponding to [6,1,2,4,5,3]
! is 603, that corresponding to [6,5,4,3,2,1] is 719. The index of any given
! permutation of the six indices can be evaluated conveniently by a call to
! subroutine, indexofperm. 
!
! The tile-1 of the protoframe (which is NOT necessarily the user's defined
! tile-1) is always the one that defines the orientation of the cubic frame
! in space, so it is not necessary to separately define the orientation of
! the local (x,y) map coordinates of that special tile. But the orientations
! of the five remaining tiles, 2--6, of the protoframe can each be oriented
! in any one of the four possible orientations, measured in right-angles
! counterclockwise from the protoframe standard. Thus, the series of
! rotation indices corresponding to the protoframe tiles 2--6, each index being
! between 0 and 3, defines a 5-digit "quaternary", or base-4, numbers. If we 
! take protoframe tile-6 to correspond to the "units", tile-5 to the "fours"
! and so on (tile-1 is always initialized with a "0"-digit of the
! "1024s", which is why the rotation index is <1024). Evaluating this single
! equivalent number "quat" between 0 and 1023, can be conveniently done
! for any explicit rotations code, by a call to subroutine, indexofquat.
!
! The orientation in space of the cube as a whole is done with complete
! generality by: (1) specifying the longitude (degrees positive in the easterly
! sense) and the latitude (degrees positive in the northerly sense) of the
! center of the protoframe's tile-1; (2) specifying the azimuthal angle of 
! twist (degrees positive in the counterclockwise sense) of the actual 
! local x-axis direction relative to the local geographical East, or the
! equivalent direction in the special limiting cases where this tile
! center is either the North Pole or the South Pole (lat = +90 or -90).
!
! To sum up, the user-defined cubic arrangement (the "topology") of six map
! panels, and the orientation of the whole in space, can be specified
! by providing to an initialization subroutine, which is here called 
! "ini_orientation", the permutation index (in [0:719]), rotations index 
! (in [0:1023]), the protoframe tile-1 longitude, latitude and azimuth (all 
! in degrees as defined above), or, since the azimuth is very frequently
! zero, a call to the same routine that omits azimuth, in which case it
! defaults to zero.
! This basic topology and orientation can be
! done independently for up to NNTRAN separate transformations, indexed by
! ITRAN from 1 to NNTRAN, as the first argument to ini_orientation.
!
! Once a transformation indexed by ITRAN has been sufficiently initialized
! we can conduct transformations between the map coordinate pair, XM=(x,y) of
! a specied user-defined tile or map panel, IPAN, to the earth-centered
! cartesian 3-vector, XC (and thence to latitide and longitude if desired):
! call xmtoxc(XM,XC,XCD,IPAN,ITRAN)
! and the 3*2 matrix, XCD, is also provided as the Jacobian, d(XC)/d(XM) of
! this transformation, where XC is always a unit vector (i.e., we take the
! radius of the spherical earth to be 1). We can omit ITRAN if, by default, we
! assume it to be ITRAN=1, as we might very often do when we are only
! concerned with a single cubed-sphere transformation).
!
! Conversely, we might wish to know which map panel a given Cartesian 3-vector
! XC corresponds to, and which coordinate pair XM on this map panel (or tile)
! is defined within the standard range [-1.,1.] for each component. In this 
! case, we initial set IPAN=0, and call the routine to invert the mapping:
! call xctoxm(XC,XM,XMD,IPAN,ITRAN).
! As expected, XMD is a kind of Jacobian, except that concept is not formally
! defined when going from a 3-vector to a 2-vector, so instead, it is the
! usual generalized inverse of the earlier defined XCD.
!
! If we want to FORCE the tile, or panel, index IPAN, to be a specified one,
! regardless of whether the coordinates (x,y)=XM actually will end up lying
! within the proper range, [-1.,1.], then we set IPAN to be the desired
! index in the range 1--6 and:
! call xctoxm(XC,XM,XMD,IPAN,ITRAN)
! otherwsie exactly as before. This is very useful during interpolations
! when it becomes necessary to know the coordinates of a point (even a grid 
! point) inside panel IPAN_A with respect to the map coorinates of the
! neighboring panel, say IPAN_B.
!=============================================================================
use pkind,      only: dp                     ! double precision real kind
use pietc,      only: T,F,u0,u1,u2,u4,pi,dtor! Useful constants
implicit none
private
public:: fin_generic,ini_orientation,ini_combination,     &
         ini_gridtype,ini_schmidt,query_cogs,             &
         xmtoxc,xctoxm,                                   &
         permofindex,indexofperm,quatofindex,indexofquat, &
         recpq,getpqawrtb,pqawrtb

integer,parameter             :: m=6      ! number of tiles in the cube
integer,parameter             :: nntran=4 ! Allowance of distict transformations
integer,parameter             :: nn=8     ! Max degree of continuity
real(dp),parameter            :: piq=pi/4,u4opi=u4/pi
real(dp),dimension(3,3,nntran):: rot6
real(dp),dimension(3,3,6)     :: rotp
real(dp),dimension(nn,nntran) :: gtv
real(dp),dimension(nntran)    :: gtkay,gtphit,gtat,gtcb,gtatcb,gtschm
integer,dimension(6,nntran)   :: kgperm,kgpermi
integer,dimension(6,nntran)   :: kgrot
integer,dimension(nntran)     :: gtcont
logical,dimension(nntran)     :: linic,linio,linig,linis
logical                       :: lini
data lini/.false./

interface ini_generic;    module procedure ini_generic;          end interface
interface fin_generic;    module procedure fin_generic;          end interface
interface ini_orientation
   module procedure ini_o_rr,ini_o_rri,ini_o_rrr,ini_o_rrri;     end interface
interface ini_combination;module procedure ini_c_ii,ini_c_iii;   end interface
interface ini_gridtype
   module procedure ini_g_r,ini_g_ri,ini_g_ir,ini_g_iri;         end interface
interface ini_schmidt;    module procedure ini_s_r,ini_s_ri;     end interface
interface query_cogs
   module procedure query_cogs,query_cogs_i;                     end interface
interface xmtoxc;         module procedure xmtoxc,xmtoxc_i;      end interface
interface xctoxm;         module procedure xctoxm,xctoxm_i;      end interface
interface eatomn;         module procedure eatomn;               end interface
interface mntoea;         module procedure mntoea;               end interface
interface eatogg;         module procedure eatogg;               end interface
interface ggtoea;         module procedure ggtoea;               end interface
interface phitoa;         module procedure phitoa;               end interface
interface atophi;         module procedure atophi;               end interface
interface getbs;          module procedure getbs;                end interface
interface gety;           module procedure gety;                 end interface
interface getphipe;       module procedure getphipe;             end interface
interface getphipo;       module procedure getphipo;             end interface

interface permofindex;    module procedure permofindex;          end interface
interface indexofperm;    module procedure indexofperm;          end interface
interface quatofindex;    module procedure quatofindex;          end interface
interface indexofquat;    module procedure indexofquat;          end interface
interface recpq;          module procedure recpq;                end interface
interface getpqawrtb;     module procedure getpqawrtb;           end interface
interface pqawrtb;        module procedure pqawrtb;              end interface

contains

!==============================================================================
subroutine ini_generic!                                           [ini_generic]
!==============================================================================
! Initialize the array rotp used by all the gnomonic transformations, but
! set all the particular status flags, except lini, to "uninitialized".
!==============================================================================
rotp(:,1,1)=(/ u0, u1,u0/);rotp(:,2,1)=(/ u1,u0,u0/);rotp(:,3,1)=(/ u0, u0,-u1/)
rotp(:,1,2)=(/ u0, u1,u0/);rotp(:,2,2)=(/ u0,u0,u1/);rotp(:,3,2)=(/ u1, u0, u0/)
rotp(:,1,3)=(/-u1, u0,u0/);rotp(:,2,3)=(/ u0,u0,u1/);rotp(:,3,3)=(/ u0, u1, u0/)
rotp(:,1,4)=(/ u0,-u1,u0/);rotp(:,2,4)=(/ u0,u0,u1/);rotp(:,3,4)=(/-u1, u0, u0/)
rotp(:,1,5)=(/ u1, u0,u0/);rotp(:,2,5)=(/ u0,u0,u1/);rotp(:,3,5)=(/ u0,-u1, u0/)
rotp(:,1,6)=(/ u0, u1,u0/);rotp(:,2,6)=(/-u1,u0,u0/);rotp(:,3,6)=(/ u0, u0, u1/)
lini=T
linic=F ! combinatorics not initialized for any transformation
linio=F ! orientation of principal tile not initialized for any transformation
linig=F ! grid type not initialized for any transformation
linis=F ! Schmidt factor not initialized for any transformation
end subroutine ini_generic

!=============================================================================
subroutine fin_generic!                                          [fin_generic]
!=============================================================================
! Finalization -- reset all status flags, lini and lormode 
!  to "uninitialized".
!============================================================================
lini=F
linic=F
linio=F
linig=F
linis=F
end subroutine fin_generic

!=============================================================================
subroutine ini_c_ii(gperm,grot)!                             [ini_combination]
!=============================================================================
integer, intent(in ):: gperm,grot
call ini_c_iii(gperm,grot,1)
end subroutine ini_c_ii
!=============================================================================
subroutine ini_c_iii(gperm,grot,itran)!                      [ini_combination]
!=============================================================================
! Specify the combinatorial arrangements of user-defined panels indices in
! the standard protoframe (via the permutation-of-6 index, gperm) and the
! relative rotations of the tiles in protoframe faces, 2--5, via the
! equivalent, grot, of the 5-digit quaternary number, grot.
!=============================================================================
integer, intent(in ):: gperm,grot,itran
if(.not.lini)call ini_generic
if(itran<1 .or. itran>nntran)stop 'In ini_combination; itran out of range'
if(gperm<0 .or. gperm>719)   stop 'In ini_combination; gperm out of range'
if(grot<0  .or. grot>1023)   stop 'In ini_combination; grot out of range'
if(.not.lini)call ini_generic
call permofindex(gperm,kgperm(:,itran))
call quatofindex(grot ,kgrot (:,itran))
call recpq(kgperm(:,itran), kgpermi(:,itran))
linic(itran)=T
end subroutine ini_c_iii

!=============================================================================
subroutine ini_o_rr(lon,lat)!                                [ini_orientation]
!=============================================================================
real(dp),intent(in ):: lon,lat
call ini_o_rrri(lon,lat,u0,1)
end subroutine ini_o_rr
!=============================================================================
subroutine ini_o_rri(lon,lat,itran)!                         [ini_orientation]
!=============================================================================
real(dp),intent(in ):: lon,lat
integer ,intent(in ):: itran
call ini_o_rrri(lon,lat,u0,itran)
end subroutine ini_o_rri
!=============================================================================
subroutine ini_o_rrr(lon,lat,az)!                            [ini_orientation]
!=============================================================================
real(dp),intent(in ):: lon,lat,az
call ini_o_rrri(lon,lat,az,1)
end subroutine ini_o_rrr
!=============================================================================
subroutine ini_o_rrri(lon,lat,az,itran)!                     [ini_orientation]
!=============================================================================
! Initialize the orientation in earth-centered cartesian space of the 
! protoframe in which the cubed-sphere grid geometry if referred. The center
! of the protoframe's face-1 is positioned at (degrees) longitide and
! latitude, (lon,lat) and its local map x-axis is rotated counterclockwise
! relative to local East by (degrees) angle, az. As usual, itran is the index
! of the stored transformation (assumed by default to be =1 if omitted).
!============================================================================= 
real(dp),intent(in ):: lon,lat,az
integer ,intent(in ):: itran
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: twist
real(dp)               :: rlon,clon,slon,rlat,clat,slat,raz,caz,saz
!=============================================================================
if(itran<1 .or. itran>nntran)stop 'In ini_orientation; itran out of range'
if(.not.lini)call ini_generic
if(.not.linic(itran))stop &
     'In ini_orientation; a call to ini_combination must precede this call'
rlon=lon*dtor; clon=cos(rlon); slon=sin(rlon)
rlat=lat*dtor; clat=cos(rlat); slat=sin(rlat)
raz =az *dtor; caz =cos(raz ); saz =sin(raz )
rot6(:,1,itran)=(/-slat*clon,-slat*slon,clat/)
rot6(:,2,itran)=(/-slon,clon,u0/)
rot6(:,3,itran)=(/-clat*clon,-clat*slon,-slat/)
twist(:,1)=(/caz,-saz,u0/);twist(:,2)=(/saz,caz,u0/);twist(:,3)=(/u0,u0,u1/)
rot6(:,:,itran)=matmul(rot6(:,:,itran),twist)
linio(itran)=T
linig(itran)=F
linis(itran)=F
end subroutine ini_o_rrri

!=============================================================================
subroutine ini_g_r(alpha)!                                      [ini_gridtype]
!=============================================================================
real(dp),intent(in ):: alpha
call ini_g_ri(alpha,1)
end subroutine ini_g_r
!=============================================================================
subroutine ini_g_ri(alpha,itran)!                               [ini_gridtype]
!=============================================================================
real(dp),intent(in ):: alpha
integer, intent(in ):: itran
!=============================================================================
if(.not.lini)call ini_generic
if(itran<1 .or. itran>nntran)stop 'In ini_gridtype; itran out of range'
if(.not.linio(itran))stop &
     'In ini_gridtype; a call to ini_orientation must precede this call'
if(alpha<u0)stop 'In ini_gridtype; parameter alpha must not be negative'
gtcont(itran)=0
gtcb(itran)=sqrt(alpha)
gtatcb(itran)=atan(gtcb(itran))
print '("Grid transformation,",i2," initialized as generic gnomonic type,")',&
     itran
print '("Stretching-profile parameter, alpha = [cos[beta]]^2 = ",f12.8)',alpha
linig(itran)=T
linis(itran)=F
end subroutine ini_g_ri
!=============================================================================
subroutine ini_g_ir(n,alpha)!                                    [ini_gridtype]
!=============================================================================
integer, intent(in ):: n
real(dp),intent(in ):: alpha
if(n==0)then; call ini_g_ri(alpha,1)
else
              call ini_g_iri(n,alpha,1)
endif
end subroutine ini_g_ir
!=============================================================================
subroutine ini_g_iri(n,alpha,itran)!                            [ini_gridtype]
!=============================================================================
integer ,intent(in ):: n
real(dp),intent(in ):: alpha
integer, intent(in ):: itran
!=============================================================================
if(.not.lini)call ini_generic
if(itran<1 .or. itran>nntran)stop 'In ini_gridtype; itran out of range'
if(.not.linio(itran))stop &
     'In ini_gridtype; a call to ini_orientation must precede this call'
if(n==0)then; call ini_g_ri(alpha,itran)
elseif(n<0 .or. n>8)then
   stop 'In ini_gridtype; n is out of bounds'
else
! Treat the case of the Mobius net grid of order of continuity, n, and with
! a half-width angle of the Mobius zone of alpha degrees.
if(alpha<=u0 .or. alpha>= 45)stop &
     'In ini_gridtype; half-width alpha of Mobius zone out of bounds'
gtcont(itran)=n
gtphit(itran)=(45-alpha)*dtor
call getbs(n,gtphit(itran), &
     gtat(itran),gtkay(itran),gtv(1:n,itran)); gtv(n+1:nn,itran)=0
print '("Grid transformation,",i2," initialized as Mobius net type,")',itran
print '("Order of continuity, n =",t35,i2)',n
print '("Mobius zone halfwidth [degrees]:"t35,f12.8)',alpha
print '("Transition index variable:",t35,f12.8)',gtat(itran)
endif
linig(itran)=T
linis(itran)=F
end subroutine ini_g_iri

!=============================================================================
subroutine ini_s_r(s)!                                           [ini_schmidt]
!=============================================================================
real(dp),intent(in ):: s
call ini_s_ri(s,1)
end subroutine ini_s_r
!=============================================================================
subroutine ini_s_ri(s,itran)!                                    [ini_schmidt]
!=============================================================================
! Initialize the Schmidt enhancement factor of transformation ITRAN to be S
! (This initialization is optional; the default, equivalent to S=1, is to
! not use the Schmidt transformation.)
!=============================================================================
real(dp),intent(in ):: s
integer, intent(in ):: itran
!=============================================================================
if(itran<1 .or. itran>nntran)stop 'In ini_schmidt; itran out of range'
if(.not.lini)call ini_generic
if(.not.linig(itran))stop &
     'In ini_schmidt; a call to ini_gridtype must precede this call'
gtschm(itran)=s
linis(itran)=T
end subroutine ini_s_ri

!=============================================================================
subroutine query_cogs(li,lic,lio,lig,lis,& !                      [query_cogs]
     gperm,grot,grot6,n,alpha,s)                                
!=============================================================================
logical,intent(out)                :: li,lic,lio,lig,lis
integer,dimension(6),   intent(out):: gperm,grot
real(dp),dimension(3,3),intent(out):: grot6
integer,intent(out)                :: n
real(dp),intent(out)               :: alpha,s
call query_cogs_i(li,lic,lio,lig,lis, gperm,grot,grot6,n,alpha,s, 1)
end subroutine query_cogs
!=============================================================================
subroutine query_cogs_i(li,lic,lio,lig,lis,& !                    [query_cogs]
     gperm,grot,grot6,n,alpha,s,itran)                                
!=============================================================================
! Return the parameters that describe the status of transformation, itran.
logical,intent(out)                :: li,lic,lio,lig,lis
integer,dimension(6),   intent(out):: gperm,grot
real(dp),dimension(3,3),intent(out):: grot6
integer,intent(out)                :: n
real(dp),intent(out)               :: alpha,s
integer ,intent(in )               :: itran
!=============================================================================
if(itran<1 .or. itran>nntran)stop 'In query_cogs; itran is out of bounds'
li=lini
if(.not.li)return
lic=linic(itran)
lio=linio(itran)
lig=linig(itran)
lis=linis(itran)
if(.not.lic)return
gperm=kgperm(:,itran)
grot =kgrot (:,itran)
if(.not.lio)return
grot6=rot6(:,:,itran)
if(.not.lig)return
n=gtcont(itran)! If n>0 (Mobius) it is interpreted as the degree of continuity
if(n==0)then
   alpha=gtcb(itran)**2 ! Generic smooth gnomonic type
else
   alpha=45-gtphit(itran)/dtor ! Mobius net type, alpha = halfwidth of zone
endif
if(.not.lis)return
s=gtschm(itran)
end subroutine query_cogs_i

!=============================================================================
subroutine xmtoxc(xm,xc,xcd,ipan)!                                    [xmtoxc]
!=============================================================================
! Take the transformation index to be its default, itran=1.
!=============================================================================
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
integer,                intent(in ):: ipan
!=============================================================================
call xmtoxc_i(xm,xc,xcd,ipan,1)
end subroutine xmtoxc
!=============================================================================
subroutine xmtoxc_i(xm,xc,xcd,ipan,itran)!                            [xmtoxc]
!=============================================================================
! Trasform from the 2-vector map coordinates, xm, of user-defined panel 
! (or tile), ipan, using the specification of stored transformation, itran,
! and return the earth-centered cartesian unit 3-vector, xc, and its 3*2
! jacobian matrix (derivative wrt xm), xcd.
!=============================================================================
use pmat5,only: ctoc_schm
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
integer,                intent(in ):: ipan,itran
!-----------------------------------------------------------------------------
real(dp),dimension(3,3)    :: xcd1
real(dp),dimension(2,2)    :: xm1d
real(dp),dimension(2,2)    :: twist
real(dp),dimension(3)      :: xc1
real(dp),dimension(2)      :: xm1
real(dp),dimension(nn)     :: v
real(dp)                   :: a,aa,aap,b,bb,bbp,z,zzz,aapzzz,bbpzzz,ab,&
                              cb,atcb,phit,at,kay
integer, dimension(2,2,0:3):: twists
integer                    :: i,irot,jpan,n
data twists/1,0,0,1,   0,1,-1,0,  -1,0,0,-1,  0,-1,1,0/
!=============================================================================
if(itran<1.or.itran>nntran)stop 'In xmtoxc; transformation index out of bounds'
if(.not.linig(itran))then
   print '("In xmtoxc; transformation, ",i2," has not been initialized")',itran
   stop
endif
n=gtcont(itran)
if(n==0)then
   cb=gtcb(itran)
   atcb=gtatcb(itran)
   xm1d=0;do i=1,2;call ggtoea(cb,atcb,xm(i),xm1(i),xm1d(i,i));enddo
else
   phit=gtphit(itran)
   kay=gtkay(itran)
   at=gtat(itran)
   v(1:n)=gtv(1:n,itran)
   xm1d=0
   do i=1,2;call mntoea(n,phit,at,kay,v(1:n),xm(i),xm1(i),xm1d(i,i));enddo
endif
jpan=kgpermi(ipan,itran)
irot=kgrot(jpan,itran)
twist=twists(:,:,irot)
xm1=matmul(twist,xm1)
a=tan(xm1(1)*piq); aa=a*a; aap=aa+u1
b=tan(xm1(2)*piq); bb=b*b; bbp=bb+u1
ab=a*b
z=u1/sqrt(u1+aa+bb); zzz=z**3; aapzzz=aap*zzz; bbpzzz=bbp*zzz
xc=(/a*z,b*z,z/)
xcd(:,1)=piq*(/bbp*aapzzz,  -ab*aapzzz,  -a*aapzzz/)
xcd(:,2)=piq*(/-ab*bbpzzz,  aap*bbpzzz,  -b*bbpzzz/)
xcd=matmul(xcd,twist)
xc =matmul(rotp(:,:,jpan),xc)
xcd=matmul(rotp(:,:,jpan),xcd)
xcd=matmul(xcd,xm1d)
if(linis(itran))then
   call ctoc_schm(gtschm(itran),xc,xc1,xcd1)
   xc=xc1
   xcd=matmul(xcd1,xcd)
endif
xc =matmul(rot6(:,:,itran),xc)
xcd=matmul(rot6(:,:,itran),xcd)
end subroutine xmtoxc_i

!=============================================================================
subroutine xctoxm(xc,xm,xmd,ipan)!                                    [xctoxm]
!=============================================================================
real(dp),dimension(3),  intent(in   ):: xc
real(dp),dimension(2),  intent(  out):: xm
real(dp),dimension(2,3),intent(  out):: xmd
integer,                intent(inout):: ipan
!=============================================================================
call xctoxm_i(xc,xm,xmd,ipan,1)
end subroutine xctoxm
!=============================================================================
subroutine xctoxm_i(xc,xm,xmd,ipan,itran)!                            [xctoxm]
!=============================================================================
! The inverse transformation of xmtoxc, with the pseudo-jacobian, 2*3 matrix
! xmd, being the pseudo-inverse of xcd=d(xc)/x(xm).
!=============================================================================
use pmat,         only: inv
use pmat4,        only: normalized
use pmat5,        only: ctoc_schm
real(dp),dimension(3),  intent(in   ):: xc
real(dp),dimension(2),  intent(  out):: xm
real(dp),dimension(2,3),intent(  out):: xmd
integer,                intent(inout):: ipan
integer,                intent(in   ):: itran
!-----------------------------------------------------------------------------
real(dp),dimension(6)     :: ax
real(dp),dimension(3)     :: xct,xc1
real(dp),dimension(3,2)   :: xcd
real(dp),dimension(3,3)   :: rot,xccd
real(dp),dimension(2,2)   :: em
real(dp),dimension(2)     :: xm1
real(dp),dimension(2,2)   :: twist
real(dp),dimension(nn)    :: v
real(dp)                  :: cb,atcb,kay,phit
integer, dimension(1)     :: ii
integer,dimension(2,2,0:3):: twists
integer                   :: i,irot,jpan,n
data twists/1,0,0,1,   0,1,-1,0,  -1,0,0,-1,  0,-1,1,0/
!=============================================================================
if(itran<1.or.itran>nntran)stop 'In xctoxm; transformation index out of bounds'
if(.not.linig(itran))then
   print '("In xctoxm; transformation, ",i2," has not been initialized")',itran
   stop
endif
if(ipan<0.or.ipan>6)stop 'In xctoxm; invalid ipan specified'
xccd=transpose(rot6(:,:,itran))
xct=normalized(xc)
xct=matmul(xccd,xct)
if(linis(itran))then
   call ctoc_schm(u1/gtschm(itran),xct,xc1)
   xct=xc1
endif
if(ipan==0)then
   ax=(/-xct(3),xct(1),xct(2),-xct(1),-xct(2),xct(3)/);ii=maxloc(ax); jpan=ii(1)
   ipan=kgperm(jpan,itran)
else
   jpan=kgpermi(ipan,itran)
endif
irot=kgrot(jpan,itran)
twist=twists(:,:,irot)

! Temporarily rotate the normalized cartesian vector so that it belongs 
! to tile 2 of the protoframe (whose center is at [1,0,0])
rot=matmul(rotp(:,:,2),transpose(rotp(:,:,jpan)))
xct=matmul(rot,xct)
xm=(/atan(xct(2)/xct(1)),atan(xct(3)/xct(1))/)/piq
xm=matmul(xm,twist)
n=gtcont(itran)
if(n==0)then
   cb=gtcb(itran)
   atcb=gtatcb(itran)
   em=0;do i=1,2;call eatogg(cb,atcb,xm(i),xm1(i),em(i,i));enddo
   xm=xm1
else
   phit=gtphit(itran)
   kay =gtkay(itran)
   v(1:n)=gtv(1:n,itran)
   em=0;do i=1,2;call eatomn(n,phit,kay,v(1:n),xm(i),xm1(i),em(i,i));enddo
   xm=xm1
endif
call xmtoxc(xm,xct,xcd,ipan,itran)
! Construct the generalized inverse of xcd and call it xmd:
em=matmul(transpose(xcd),xcd)
call inv(em)
xmd=matmul(em,transpose(xcd))
end subroutine xctoxm_i

!==============================================================================
subroutine eatomn(n,phit,kay,v,xea,xmn,xmnd)!                          [eatomn]
!==============================================================================
! Equiangular (EA) to Mobius-net (MN) gnomonic map coordinates, together
! with the jacobian, xmnd=d(xmn)/d(xea), of this transformation.
!==============================================================================
integer,                intent(in ):: n
real(dp),               intent(in ):: phit,kay
real(dp),dimension(n),  intent(in ):: v
real(dp),               intent(in ):: xea
real(dp),               intent(out):: xmn,xmnd
!------------------------------------------------------------------------------
real(dp):: phi,da,dda
!==============================================================================
phi=xea*piq
call phitoa(n,phit,kay,v, phi,xmn,da,dda)
xmnd=da*piq
end subroutine eatomn

!==============================================================================
subroutine mntoea(n,phit,at,kay,v,xmn,xea,xead)!                       [mntoea]
!==============================================================================
! Mobius-net (MN) to equiangular (EA) gnomonic map coordinates, together
! with the jacobian, xead=d(xea)/d(xmn) of this transformation.
!==============================================================================
integer,                intent(in ):: n
real(dp),               intent(in ):: phit,at,kay
real(dp),dimension(n),  intent(in ):: v
real(dp),               intent(in ):: xmn
real(dp),               intent(out):: xea,xead
!-----------------------------------------------------------------------------
real(dp):: phi,dphi,ddphi
!==============================================================================
call atophi(n,phit,at,kay,v, xmn, phi,dphi,ddphi)
xea =phi *u4opi
xead=dphi*u4opi
end subroutine mntoea

!=============================================================================
subroutine eatogg(cb,atcb,xea,xgg,xggd)!                              [eatogg]
!=============================================================================
! Equiangular (EA) to generic gnomonic (GG) map coordinates, together 
! with the jacobian, xggd=d(xgg)/d(xea) of this transformation.
!==============================================================================
real(dp),intent(in ):: cb,atcb
real(dp),intent(in ):: xea
real(dp),intent(out):: xgg,xggd
!-----------------------------------------------------------------------------
real(dp):: tea
!=============================================================================
tea=tan(piq*xea)
xgg=atan(cb*tea)/atcb
xggd=(piq*cb/atcb)*(1+tea**2)/(1+cb**2*tea**2)
end subroutine eatogg

!=============================================================================
subroutine ggtoea(cb,atcb,xgg,xea,xead)!                              [ggtoea]
!=============================================================================
! Generic gnomonic (GG) to equiangular (EA) gnomonic map coordinates, together
! with the jacobian, xead=d(xea)/d(xgg) of this transformation.
!=============================================================================
real(dp),intent(in ):: cb,atcb
real(dp),intent(in ):: xgg
real(dp),intent(out):: xea,xead
!-----------------------------------------------------------------------------
real(dp):: tgg
!=============================================================================
tgg=tan(atcb*xgg)
xea =u4opi*atan(tgg/cb)
xead=u4opi*cb*atcb*(u1+tgg**2)/(cb**2+tgg**2)
end subroutine ggtoea

!=============================================================================
subroutine phitoa(n,phit,kay,v, phi,a,da,dda)!                        [phitoa]
!=============================================================================
! For the Mobius net grids, with the order of continuity, n, and the
! transitional angle, phit, from the tile median marking the onset of the 
! Mobius zone, and given a general angle from the median, phi, return the
! corresponding index variable, a, as well as its first and second derivatives,
! da and dda, wrt phi.
!=============================================================================
integer,              intent(in ):: n
real(dp),             intent(in ):: phit,kay
real(dp),dimension(n),intent(in ):: v
real(dp),             intent(in ):: phi
real(dp),             intent(out):: a,da,dda
!-----------------------------------------------------------------------------
real(dp),dimension(n):: phipo,phipe
real(dp)             :: aphi,tphi,phic2,s2r,c2r
!=============================================================================
aphi=abs(phi)
if(aphi>phit)then
   tphi=tan(aphi)
   a=u1+kay*log(tphi)
   phic2=(piq-aphi)*2; s2r=sin(phic2); c2r=cos(phic2)
   da=kay*u2/c2r
   dda=-kay*u4*s2r/c2r**2
else
   call getphipo(n,aphi,phipo)! odd powers
   call getphipe(n,aphi,phipe)! even powers
   a =dot_product(v,phipo)
   da=dot_product(v,phipe)
   if(n==1)then
      dda=0
   else
      dda=dot_product(v(2:),phipo(1:n-1))
   endif
endif
if(aphi*phi<0)then
   a  =-a
   dda=-dda
endif
end subroutine phitoa

!==============================================================================
subroutine atophi(n,phit,at,kay,v,a,phi,dphi,ddphi)!                   [atophi]
!==============================================================================
! Given the order of continuity, n, at the edge of the Mobius zone,
! the Mobius zone transition angle, phit, and its corresponding index 
! variable, at, the calibration constant, kay, the n spline coefficients,
! v, and a given particular value of index variable, a, evaluate the
! corresponding phi, and its first and second derivatives, dphi and ddphi,
!  wrt a. 
!==============================================================================
integer,              intent(in ):: n
real(dp),             intent(in ):: phit,at,kay
real(dp),dimension(n),intent(in ):: v
real(dp),             intent(in ):: a
real(dp),             intent(out):: phi,dphi,ddphi
!------------------------------------------------------------------------------
integer, parameter:: nit=20 ! <- Newton iterations allowance
real(dp),parameter:: eps=1.e-12! <- Newton convergence criterion
real(dp)          :: aa,ac,e,ee,eep,eec,ares,ar,da,dda
integer           :: i
!==============================================================================
if(a==0)then
   phi=0;    dphi=u1/v(1); ddphi=0; return
endif
aa=abs(a)
if(aa>=at)then
   ac=(u1-aa)/kay
   e=exp(-ac)
   ee=e*e
   eep=kay*(u1+ee)
   eec=u1-ee
   phi  =atan(e)
   dphi =e/eep
   ddphi=e*eec/eep**2
else
   phi=aa/v(1)
   ares=1
   do i=1,nit
      call phitoa(n,phit,kay,v, phi,ar,da,dda)
      if(abs(ares)<eps)exit
      ares=ar-aa
      phi=phi-ares/da
   enddo
   if(i>nit)print '("In atophi; WARNING: Newton iterations not converging")'
   dphi = u1/da
   ddphi=-dda/da**3
endif
if(a*aa<0)then ! Flip some signs for the negative angle cases:
   phi  =-phi
   ddphi=-ddphi
endif
end subroutine atophi

!==============================================================================
subroutine getphipe(n,phi,phipe)!                                    [getphipe]
!==============================================================================
! Get n-vector of the first n even powers of phi divided by the factorial
! of each exponent.
!==============================================================================
integer,              intent(in ):: n
real(dp),             intent(in ):: phi
real(dp),dimension(n),intent(out):: phipe
!------------------------------------------------------------------------
integer:: i,k
!==============================================================================
phipe(1)=u1
do i=2,n
   k=(i*2-2)*(i*2-3)
   phipe(i)=phipe(i-1)*phi**2/k
enddo
end subroutine getphipe

!==============================================================================
subroutine getphipo(n,phi,phipo)!                                    [getphipo]
!==============================================================================
! Get n-vector of the first n odd powers of phi divided by the factorial
! of each exponent.
!==============================================================================
integer,              intent(in ):: n
real(dp),             intent(in ):: phi
real(dp),dimension(n),intent(out):: phipo
!------------------------------------------------------------------------
integer:: i,k
!==============================================================================
phipo(1)=phi
do i=2,n
   k=(i*2-1)*(i*2-2)
   phipo(i)=phipo(i-1)*phi**2/k
enddo
end subroutine getphipo

!==============================================================================
subroutine getbs(n,phit,at,kay,v)!                                      [getbs]
!==============================================================================
! For a given degree of continuity, n, with 1 <= n <= nn, at the transition
! angle, phit, with 0 <= phit <= pi/4, return the normalizing constant, kay,
! for calibrating the Mobius net part of the index variable function a(phi)
! (where phi>phit) and the n successive odd-degree coefficients, v, for the
! polynomial spline a(phi) for phi < phit. Note that the value, and first
! n derivatives, for the spline and Mobius portions of a(phi), match at phit
! itself. As a convenience, also output the transitional index variable, at,
! that corresponds to phit.
!==========================================================================
use pmat, only: inv
integer,              intent(in ):: n
real(dp),             intent(in ):: phit
real(dp),             intent(out):: at,kay
real(dp),dimension(n),intent(out):: v
!--------------------------------------------------------------------------
real(dp),dimension(0:n,0:n):: em
real(dp),dimension(0:nn)   :: y
real(dp),dimension(0:n)    :: vs
real(dp),dimension(0:n*2-1):: phitp
integer                    :: i,j,k,n2m
!===========================================================================
if(n>nn)stop 'In getbs; n must not exceed its max allowed value, nn'
n2m=n*2-1
phitp(0)=1
do k=1,n2m
   phitp(k)=phitp(k-1)*phit/k
enddo
call gety(phit,y)

em=0
do j=1,n
do i=0,n
   k=2*j-i-1
   if(k<0)cycle
   em(i,j)=phitp(k)
enddo
enddo
em(:,0)=-y(0:n)
call inv(em)
vs=0; vs(0)=1
vs=matmul(em,vs)
kay=vs(0); v=vs(1:n)
at=u1+kay*y(0)
end subroutine getbs

!============================================================================
subroutine gety(x,y)!                                                  [gety]
!============================================================================
! Return the value, and the first 8 derivatives, of the function y=ln(tan(x))
!============================================================================
real(dp),               intent(in ):: x
real(dp),dimension(0:8),intent(out):: y
!----------------------------------------------------------------------------
real(dp):: x2,cx2,sx2,q,p2,p4,p6
!============================================================================
x2=x*2; cx2=cos(x2); sx2=sin(x2); p2=sx2*sx2; p4=p2*p2; p6=p2*p4
y(0)=log(tan(x))
q=-2/sx2
y(1)=-q;                            q=-q*y(1)
y(2)=-q*cx2;                        q=-q*y(1)
y(3)=-q*(   2-     p2);             q=-q*y(1)
y(4)=-q*(   6-     p2)*cx2;         q=-q*y(1)
y(5)=-q*(  24-  20*p2+    p4);      q=-q*y(1)
y(6)=-q*( 120-  60*p2+    p4)*cx2;  q=-q*y(1)
y(7)=-q*( 720- 840*p2+182*p4-p6);   q=-q*y(1)
y(8)=-q*(5040-4200*p2+546*p4-p6)*cx2
end subroutine gety

!-----------------------------
! Routines that relate to the combinatorial aspects of specifying the
! arrangement of user-defined tiles and their orienations relative to the
! standard geometry of what we refer to as the cubic protoframe:

!=============================================================================
subroutine permofindex(i,perm)!                                  [permofindex]
!=============================================================================
! Deliver the Ith from the list of lexicographically-ordered permutations
! of the first m integers, where the list starts at index 0 with the
! default permutation, {1,2,...,M}, where M=6 for the cube.
! The index is first represented in mixed-base {1,2,3..}
! form, MIXED,  and PERM initially the default permutation. Then elements of 
! MIXED dictate how each element in turn of PERM is chosen from what
! numbers remain. 
!=============================================================================
integer,             intent(in ):: i
integer,dimension(m),intent(out):: perm
!-----------------------------------------------------------------------------
integer,dimension(m):: mixed
integer             :: j,k,L,Lp
!=============================================================================
if(i<0)stop 'In permofindex; index must not be negative'
L=i
mixed(1)=0; perm(1)=1
do j=2,m; perm(j)=j; Lp=L/j; mixed(j)=L-Lp*j; L=Lp; enddo
if(L>0)stop 'In permofindex; index is too large'
do j=1,m
   k=j+mixed(m+1-j); L=perm(k); perm(j+1:k)=perm(j:k-1); perm(j)=L
enddo
end subroutine permofindex

!=============================================================================
subroutine indexofperm(perm,i)!                                  [indexofperm]
!=============================================================================
! Return the index, I, that this permutation, PERM, of the first M=6 integers
! would have in a lexicographically-ordering of them, where the first 
! permutation, [1,2,..] has index 0. E.g., when M=3, lexicographically ordered
! list, I going from 0 to 5, is [1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]
!=============================================================================
integer,dimension(m),intent(in ):: perm
integer,             intent(out):: i
!-----------------------------------------------------------------------------
integer:: j,jc,k,p,ff
!=============================================================================
ff=1; i=0
do j=1,m-1
   ff=ff*j! factorial
   jc=m-j; p=perm(jc); do k=jc+1,m; if(p>perm(k))i=i+ff; enddo
enddo
end subroutine indexofperm

!=============================================================================
subroutine quatofindex(i,quat)!                                  [quatofindex]
!=============================================================================
! With M=6, and freedom to rotate tile maps, 2 -- m, of the protoframe 
! Return the quaternary (M-1)-digit, QUAT (ordered such that the units are in
! the last position, M=6) that corresponds to index I.
!=============================================================================
integer,             intent(in ):: i
integer,dimension(m),intent(out):: quat
!-----------------------------------------------------------------------------
integer:: ii,j
!=============================================================================
ii=i
quat(1)=0
do j=m,2,-1
   quat(j)=mod(ii,4)
   ii=ii/4
enddo
if(ii>0)stop 'In quatofindex; i is too large for only (m-1) digits'
end subroutine quatofindex

!=============================================================================
subroutine indexofquat(quat,i)!                                  [indexofquat]
!=============================================================================
! Return the index, I, that corresponds to the M quaternary digits, QUAT,
! when ordered such that the units are in position M, the 4's in position
! M-1, and so on.
!=============================================================================
integer,dimension(m),intent(in ):: quat
integer,             intent(out):: i
!-----------------------------------------------------------------------------
integer:: j
!==============================================================================
i=0
do j=2,m ! <- Ignore digit quat(1) corresponding to "1024's"; it should be 0.
   i=4*i+quat(j)
enddo
end subroutine indexofquat

!=============================================================================
subroutine recpq(pa,recpa)!                                            [recpq]
!=============================================================================
! For a relative permutation, pa, 
! return the corresponding codes for the reciprocal relationship, recpa.
!=============================================================================
integer,dimension(m),  intent(in ):: pa
integer,dimension(m),  intent(out):: recpa
!-----------------------------------------------------------------------------
integer:: i,p
!=============================================================================
do i=1,m; p=pa(i); recpa(p)=i; enddo
end subroutine recpq

!=============================================================================
subroutine getpqawrtb(pa,qa,pb,qb, pawrtb,qawrtb)!                [getpqawrtb]
!=============================================================================
! Given the explicit 6-index permutation, pa, and rotation code, qa, from
! the standard protoframe, of frame A, and the corresponding permutation, pb,
! and rotation code, qb, of frame B, return the relative permutation 
! and rotation code of A with respect to frame B.
!=============================================================================
integer,dimension(m),intent(in ):: pa,pb
integer,dimension(m),intent(in ):: qa,qb
integer,dimension(m),intent(out):: pawrtb
integer,dimension(m),intent(out):: qawrtb
!-----------------------------------------------------------------------------
integer:: a,b,i
!=============================================================================
do i=1,m
   a=pa(i)
   b=pb(i)
   pawrtb(b)=a
   qawrtb(b)=modulo(qa(i)-qb(i),4)
enddo
end subroutine getpqawrtb

!=============================================================================
subroutine pqawrtb(ib,jb,pb,pawrtb,qawrtb,ia,ja,pa,qa)!              [pqawrtb]
!=============================================================================
! Given the permutation, pawrtb, defining the cubic tile indices of frame A
! with respect to frame B, and the relative rotations (in right angles), qawrtb,
! of A's map coordinates relative to the 2nd through 6th of B's map coordinates,
! then, for tile-centered integer grid coordinate, (ib,jb), in tile pb of 
! frame B, return the corresponding coordinates, (ia,ja), of this grid
! coordinates in frame A together with the tile index, pa, there, and relative
! map rotation (in right angles) of frame A with respect to frame B.
!============================================================================= 
integer,             intent(in ):: ib,jb,pb
integer,dimension(m),intent(in ):: pawrtb,qawrtb
integer,             intent(out):: ia,ja,pa,qa
!=============================================================================
pa=pawrtb(pb)
qa=qawrtb(pb)
select case(qa)
case(0);   ia= ib;   ja= jb
case(1);   ia= jb;   ja=-ib
case(2);   ia=-ib;   ja=-jb
case(3);   ia=-jb;   ja= ib
end select
end subroutine pqawrtb

end module pgnc




