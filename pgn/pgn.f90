!
!                                      ***************************************
!                                      *      module pgn                     *
!                                      *      R. J. Purser                   *
!                                      *      NOAA/NCEP/EMC July  2018       *
!                                      *      jim.purser@noaa.gov            *
!                                      ***************************************
!
! A suite of routines for generating the coefficients needed to define the
! gnomonic cubed-sphere transformations and their derivatives.
!
! DIRECT DEPENDENCIES
! Libraries[their Modules]: pmat[pmat, peuc, pgeo]
! Additional Modules      : pgn_tables, pietc, pkind
!
!=============================================================================
module pgn
!=============================================================================
! Three types of cubed-sphere gnomonic grid transformations, and the 
! corresponding Schmidt refinements of these transformations, are accommodated
! by this suite: (1) the Equiangular (EA); (2) the Edge-uniform (ED); and (3)
! the Mobius net (MN) types.
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
! in which the short arrow indicate the direction of the local "x" coordinate
! of the right-handed pair, (x,y) of local map coordinates. Generally, it is
! assumed that the local coordinates have their origins at the tile centers
! and both x and y extend between -1 and +1 from edge to edge. The particular
! arrangement of the user's choice of panel indexing is specified by a 
! permutation of the six integers 1--6 that label the panels, so that the
! permutation, [6,1,2,4,5,3] puts user-defined label "6" in the position of
! the above protoframe's tile-1, and so on. Each permutation defined in this
! way can be lexicographically ordered, and therefore indexed from 0 to 719
! (the number of possible permuations of six objects is 6!=720). The 
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
! in any one of the four possible orientations, measured in rightangles
! counterclockwise from the protoframe standard. Thus, the series of
! rotation indices corresponding to the protoframe tiles 2--6, each index being
! between 0 and 3, defines a 5-digit "quaternary", or base-4, numbers. If we 
! take protoframe tile-6 to correspond to the "units", tile-5 to the "fours"
! and so on (tile-1 is always initialized with a "0"-digit of the
! "1024s", which is why the rotation index is <1024). Evaluating this single
! equivalent number "quat" between 0 and 1023, can be conveniently done
! for any explicit rotations code, by a call to subroutine indexofquat.
!
! The orientation in space of the cube as a whole is done with complete
! generality by specifying the longitude (degrees positive in the easterly
! sense) and the latitude (degrees positive in the northerly sense) of the
! center of the protoframe's tile-1, together with an azimuthal angle of 
! twist (degrees positive in the counterclockwise sense) of the actual 
! local x-axis direction relative to the local geographical East, or the
! equivalent direction in the special limiting cases where this tile
! center is either the North Pole or the South Pole, when the defining
! latitude --> +90 or -90 with longitude kept fixed in some continuous 
! progression.
!
! To sum up, the user-defined cubic arrangement (the "topology") of six map
! panels, and the orientation of the whole in space, can be specified
! by providing to an initialization subroutine, which is here called 
! "ini_orientation", the permutation index (in [0:719]), rotations index 
! (in [0:1023]), the protoframe tile-1 longitude, latitude and azimuth (all 
! in degrees as defined above), or, since the azimuth is very frequently
! zero, a call to the same routine that omits this argument to stil obtain
! this zero azimuth default. This basic topology and orientation can be
! done independently for up to NNTRAN separate transformations, indexed by
! ITRAN from 1 to NNTRAN, as the first argument to ini_orientation.
!
! For each transformation, once the topology/orientation is initialized
! (not before!) we must next define which type of gnomonic grid and 
! associated map, or "grid index coordinate" we wish this transformation 
! to correspond to. We do this by a call to one of the three subroutines
! as follows:
! call ini_eqcube(ITRAN) -- make transformation ITRAN of the equiangular type
! call ini_edcube(ITRAN) -- make of the edge-unform type (used in the FV3)
! call ini_mncube(ITRAN,MOBZONE) -- make it of the Mobius Net type with a
!                           half-width angle (degrees) of its Mobius zones
!                           of MOBZONE.
! An optional final initialization routine allows a further SCHMIDT
! conformal refinement transformation to be applied:
! call ini_schmidt(ITRAN,S)
! where S is the refinement factor (S=1 mean no refinement). The "focus" of
! the Schmidt refinement is alway the center of the protoframe's tile-1.
!
! Once a transformation indexed by ITRAN has been sufficiently initialized
! we can conduct transformations between the map coordinate pair, XM=(x,y) of
! a specied user-defined tile or map panel, IPAN, to the earth-centered
! cartesian 3-vector, XC (and thence to latitide and longitude is desired):
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
use pgn_tables, only: nntran,rotp,rot6,linigtran,lormode,kgmode
implicit none
integer,parameter:: m=6 ! <- Number of permutable tiles (faces) of the cube
private
public:: fin_generic,ini_orientation,ini_eacube,ini_edcube,ini_mncube, &
         ini_schmidt,  xmtoxc,xctoxm,                                  &
         permofindex,indexofperm,quatofindex,indexofquat,              &
         recpq,getpqawrtb,pqawrtb

interface ini_generic;    module procedure ini_generic;          end interface
interface fin_generic;    module procedure fin_generic;          end interface
interface ini_orientation;module procedure ini_or_rr,ini_or_rrr; end interface
interface ini_eacube;     module procedure ini_eacube;           end interface
interface ini_edcube;     module procedure ini_edcube;           end interface
interface ini_mncube;     module procedure ini_mncube;           end interface
interface ini_schmidt;    module procedure ini_schmidt;          end interface

interface xmtoxc;         module procedure xmtoxc,xmtoxc_i;      end interface
interface xctoxm;         module procedure xctoxm,xctoxm_i;      end interface
interface mobnetcoefs;    module procedure mobnetcoefs;          end interface
interface eatomn;         module procedure eatomn;               end interface
interface mntoea;         module procedure mntoea;               end interface
interface eatoed;         module procedure eatoed;               end interface
interface edtoea;         module procedure edtoea;               end interface
interface phitoa;         module procedure phitoa;               end interface
interface atophi;         module procedure atophi;               end interface

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
! set all the particular status flags, except linigtran, to "uninitialized".
!==============================================================================
use pgn_tables,only: linigtran,lormode,rotp
!==============================================================================
rotp(:,1,1)=(/ u0, u1,u0/);rotp(:,2,1)=(/ u1,u0,u0/);rotp(:,3,1)=(/ u0, u0,-u1/)
rotp(:,1,2)=(/ u0, u1,u0/);rotp(:,2,2)=(/ u0,u0,u1/);rotp(:,3,2)=(/ u1, u0, u0/)
rotp(:,1,3)=(/-u1, u0,u0/);rotp(:,2,3)=(/ u0,u0,u1/);rotp(:,3,3)=(/ u0, u1, u0/)
rotp(:,1,4)=(/ u0,-u1,u0/);rotp(:,2,4)=(/ u0,u0,u1/);rotp(:,3,4)=(/-u1, u0, u0/)
rotp(:,1,5)=(/ u1, u0,u0/);rotp(:,2,5)=(/ u0,u0,u1/);rotp(:,3,5)=(/ u0,-u1, u0/)
rotp(:,1,6)=(/ u0, u1,u0/);rotp(:,2,6)=(/-u1,u0,u0/);rotp(:,3,6)=(/ u0, u0, u1/)
linigtran=T
lormode  =F
end subroutine ini_generic

!=============================================================================
subroutine fin_generic!                                          [fin_generic]
!=============================================================================
! Finalization -- reset all status flags, linigtran and lormode 
!  to "uninitialized".
!============================================================================
use pgn_tables,only: linigtran,lormode
linigtran=F
lormode  =F
end subroutine fin_generic

!=============================================================================
subroutine ini_or_rr(itran,gperm,grot,lon,lat)!              [ini_orientation]
!=============================================================================
integer, intent(in ):: itran,gperm,grot
real(dp),intent(in ):: lon,lat
!=============================================================================
call ini_or_rrr(itran,gperm,grot,lon,lat,u0)
end subroutine ini_or_rr

!=============================================================================
subroutine ini_or_rrr(itran,gperm,grot,lon,lat,az)!          [ini_orientation]
!=============================================================================
use pgn_tables,   only: nntran,kgperm,kgpermi,kgrot,kgroti, &
                        linigtran,lormode,kgmode,rot6
integer, intent(in ):: itran,gperm,grot
real(dp),intent(in ):: lon,lat,az
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: twist
real(dp)               :: rlon,clon,slon,rlat,clat,slat,raz,caz,saz
!=============================================================================
if(itran<1 .or. itran>nntran)stop 'In ini_orientation; itran out of range'
if(gperm<0 .or. gperm>719)stop 'In orient_gntran; gperm out of range'
if(grot<0  .or. grot>1023)stop 'In orient_gntran; grot out of range'
if(.not.linigtran)call ini_generic
call permofindex(gperm,kgperm(:,itran))
call quatofindex(grot ,kgrot (:,itran))
call recpq(kgperm(:,itran),kgrot(:,itran), kgpermi(:,itran),kgroti(:,itran))
rlon=lon*dtor; clon=cos(rlon); slon=sin(rlon)
rlat=lat*dtor; clat=cos(rlat); slat=sin(rlat)
raz =az *dtor; caz =cos(raz ); saz =sin(raz )
rot6(:,1,itran)=(/-slat*clon,-slat*slon,clat/)
rot6(:,2,itran)=(/-slon,clon,u0/)
rot6(:,3,itran)=(/-clat*clon,-clat*slon,-slat/)
twist(:,1)=(/caz,-saz,u0/);twist(:,2)=(/saz,caz,u0/);twist(:,3)=(/u0,u0,u1/)
rot6(:,:,itran)=matmul(rot6(:,:,itran),twist)
lormode(itran)=T
kgmode(itran)=0
end subroutine ini_or_rrr

!=============================================================================
subroutine ini_eacube(itran)!                                     [ini_eacube]
!=============================================================================
! Initialize the mode of transformation ITRAN to be EQUIANGULAR
!=============================================================================
use pgn_tables, only: nntran,lormode,kgmode,lschm
integer,intent(in ):: itran
if(itran<1 .or. itran>nntran)stop 'In ini_eacube; itran out of range'
if(.not.lormode(itran))then
   print'("In ini_eacube; Orientation of the cube of transformation",i2)',itran
   print'("has not been initialized.")'
   stop
endif
kgmode(itran)=1
lschm(itran)=F
end subroutine ini_eacube

!=============================================================================
subroutine ini_edcube(itran)!                                     [ini_edcube]
!=============================================================================
! Initialize the mode of transformation ITRAN to be EDGE-uniform gnomonic
!=============================================================================
use pgn_tables, only: nntran,lormode,kgmode,lschm
integer,intent(in ):: itran
if(itran<1 .or. itran>nntran)stop 'In ini_edcube; itran out of range'
if(.not.lormode(itran))then
   print'("In ini_edcube; Orientation of the cube of transformation",i2)',itran
   print'("has not been initialized.")'
   stop
endif
kgmode(itran)=2
lschm(itran)=F
end subroutine ini_edcube

!=============================================================================
subroutine ini_mncube(itran,mobzone)!                             [ini_mncube]
!=============================================================================
! Initialize the mode of transformation ITRAN to be MOBIUS-NET gnomonic
!=============================================================================
use pgn_tables, only: nntran,lormode,kgmode,lschm,mg_alpha,mg_at,mg_ae, &
     mg_c1,mg_c3,mg_cs
integer, intent(in ):: itran
real(dp),intent(in ):: mobzone
!-----------------------------------------------------------------------------
if(itran<1 .or. itran>nntran)stop 'In ini_mncube; itran out of range'
if(.not.lormode(itran))then
   print'("In ini_mncube; Orientation of the cube of transformation",i2)',itran
   print'("has not been initialized.")'
   stop
endif
mg_alpha(itran)=mobzone*dtor
call mobnetcoefs(mg_alpha(itran),mg_at(itran),mg_ae(itran), &
     mg_c1(itran),mg_c3(itran),mg_cs(itran))
kgmode(itran)=3
lschm (itran)=F
end subroutine ini_mncube

!=============================================================================
subroutine ini_schmidt(itran,s)!                                 [ini_schmidt]
!=============================================================================
! Initialize the Schmidt enhancement factor of transformation ITRAN to be S
!=============================================================================
use pgn_tables, only: nntran,kgmode,lschm,sschm
integer, intent(in ):: itran
real(dp),intent(in ):: s
!=============================================================================
if(itran<1 .or. itran>nntran)stop 'In ini_schmidt; itran out of range'
if(kgmode(itran)==0)then
   print'("In ini_schmidt; Gnomonic grid type for transformation",i2)',itran
   print'("has not been initialized")'
   stop
endif
sschm(itran)=s
lschm(itran)=T
end subroutine ini_schmidt

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
use pgn_tables,only: mg_alpha,mg_at,mg_ae,mg_c1,mg_c3,mg_cs,&
                       sschm,kgpermi,kgrot,kgmode,lschm
use pgeo,      only: ctoc_schm
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
integer,                intent(in ):: ipan,itran
!-----------------------------------------------------------------------------
real(dp),parameter         :: piq=pi/4
real(dp),dimension(3,3)    :: xcd1
real(dp),dimension(2,2)    :: xegd
real(dp),dimension(2,2)    :: twist
real(dp),dimension(3)      :: xc1
real(dp),dimension(2)      :: xm1
real(dp)                   :: a,aa,aap,b,bb,bbp,z,zzz,aapzzz,bbpzzz,ab
integer, dimension(2,2,0:3):: twists
integer                    :: irot,jpan
data twists/1,0,0,1,   0,1,-1,0,  -1,0,0,-1,  0,-1,1,0/
!=============================================================================
if(itran<1.or.itran>nntran)stop 'In xmtoxc; transformation index out of bounds'
select case(kgmode(itran))
case(1)
   xm1=xm
   xegd=u0; xegd(1,1)=u1; xegd(2,2)=u1
case(2)
   call edtoea(xm,xm1,xegd)
case(3)
   call mntoea(mg_alpha(itran),mg_at(itran),mg_ae(itran),&
        mg_c1(itran),mg_c3(itran),mg_cs(itran),xm,xm1,xegd) 
case default
   print'("In xmtoxc; transformation, ",i2," has not been initialized")',itran
   stop
end select
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
xcd=matmul(xcd,xegd)
if(lschm(itran))then
   call ctoc_schm(sschm(itran),xc,xc1,xcd1)
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
use pgn_tables, only: mg_alpha,mg_at,mg_ae,mg_c1,mg_c3,mg_cs,&
                        lschm,sschm,kgperm,kgpermi,kgrot
use pmat,         only: inv
use peuc,         only: normalized
use pgeo,         only: ctoc_schm
real(dp),dimension(3),  intent(in   ):: xc
real(dp),dimension(2),  intent(  out):: xm
real(dp),dimension(2,3),intent(  out):: xmd
integer,                intent(inout):: ipan
integer,                intent(in   ):: itran
!-----------------------------------------------------------------------------
real(dp),parameter        :: piq=pi/4
real(dp),dimension(6)     :: ax
real(dp),dimension(3)     :: xct,xc1
real(dp),dimension(3,2)   :: xcd
real(dp),dimension(3,3)   :: rot,xccd
real(dp),dimension(2,2)   :: em
real(dp),dimension(2)     :: xm1
real(dp),dimension(2,2)   :: twist
integer, dimension(1)     :: ii
integer,dimension(2,2,0:3):: twists
integer                   :: irot,jpan
data twists/1,0,0,1,   0,1,-1,0,  -1,0,0,-1,  0,-1,1,0/
!=============================================================================
if(itran<1.or.itran>nntran)stop 'In xctoxm; transformation index out of bounds'
if(.not.lormode(itran))then
   print'("In xctoxm; transformation, ",i2," has not been initialized")',itran
   stop
endif
if(ipan<0.or.ipan>6)stop 'In xctoxm; invalid ipan specified'
xccd=transpose(rot6(:,:,itran))
xct=normalized(xc)
xct=matmul(xccd,xct)
if(lschm(itran))then
   call ctoc_schm(u1/sschm(itran),xct,xc1)
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
select case(kgmode(itran))
case(2)
   call eatoed(xm,xm1,em)
   xm=xm1
case(3)
   call eatomn(mg_alpha(itran),mg_ae(itran),&
        mg_c1(itran),mg_c3(itran),mg_cs(itran),xm,xm1,em)
   xm=xm1
end select
call xmtoxc(xm,xct,xcd,ipan,itran)
! Construct the generalized inverse of xcd and call it xmd:
em=matmul(transpose(xcd),xcd)
call inv(em)
xmd=matmul(em,transpose(xcd))
end subroutine xctoxm_i

!==============================================================================
subroutine mobnetcoefs(alpha,at,ae,c1,c3,cs)!                     [mobnetcoefs]
!==============================================================================
! From the given Mobius zone width, alpha (in radians), return the 
! coefficients needed to evaluate the transformation and related quantities
! for the construction of a gnomonic cubic grid possessing the Mobius-net
! three-way intersection property within the hexagons, centered on each
! cube-corner, formed by the set-intersections of the three radiating Mobius
! zones there. These Mobius zones extend out either side of the great circles
! in which the edges lie, and to and angle of ALPHA on either side. The index
! coordinate, a, of the modified coordinate lines is assumed to be +1 or -1 at
! cube edges, and 0 at the medians. The transition index coordinate, AT,
! is the demarcation between the pure Mobius zone, a>AT, and the spline-
! interpolated remainder, a<AT (AT is usually assumed to be only slightly
! smaller than unity). The transformations between index variable and the
! angle within the Mobius zones are essentially the Gudermannian and
! inverse-Gudermannian functions.
! The spline that performs the interpolation is of the
! "tensioned" kind and is therefore partially composed of hyperbolic (or
! exponential) functions, and the same zone half-width, ALPHA, is adopted to
! define the inherent e-folding scale of these splines. The splines are 
! defined in angle coordinates, phi, or a scaled version, p=phi/ALPHA,
! and ensure that, as a function of p, the index variable, a, and its first
! three derivatives wrt p remain continuous. The functional form of a(p)
! is: a(p) = c1*p - c3*p**3/6 + cs*2*sinh(p)/exp(pe), where pe=(pi/4)/alpha,
! the p at the cube edge. The hyperbolic function is written this way to
! avoid overflow problems when alpha is very small, since we can express that
! term (and its even derivatives) equivalently, 
! cs*2*sinh(p)/exp(pe) = cs*(exp(p-pe)+exp(-p-pe))*tanh(p),
! and its odd-derivatives:
! d/dp ... = cs*2*cosh(p)/exp(pe) = cs*(exp(p-pe)+exp(-p-pe))
! where the exponentials in the expanded forms are smaller than 1.
!=============================================================================
real(dp),intent(in ):: alpha
real(dp),intent(out):: at,ae,c1,c3,cs
!-----------------------------------------------------------------------------
real(dp),parameter:: u8=8,piq=pi/4
real(dp)          :: alpha2,c2r,s2r,pe,da,dda,ddda,pt,tt,&
                     ptmtt,dett,c2oe,s2oe,aemat
!==============================================================================
alpha2=alpha*2! Twice the angular width of, phi_e-phi_t (the full Mobius zone)
c2r=cos(alpha2); s2r=sin(alpha2)
! To make the exponential functions easier to manipulate, we use a scaled
! angular coordinate, p = phi/alpha, when fitting the tensioned spline.
!
! Relations involving a and its p-derivatives at the transition point, phi_t:
aemat=u2*atanh(tan(alpha))            ! a_e-a_t
da  =alpha   *(u2/c2r)                ! da/dp|_t
dda =alpha**2*(-u4*s2r/c2r**2)        ! d^2 a/dp^2|_t
ddda=alpha**3*(u8*(u1+s2r**2)/c2r**3) ! d^3 a/dp^3|_t
pe=piq/alpha ! Edge, in p units
pt=pe-u1     ! Transition point, in p units
tt=tanh(pt); c2oe=exp(-u1)+exp(-pe-pt); s2oe=tt*c2oe
ptmtt=pt-tt
dett=c2oe*ptmtt ! In 2*2 matrix form, this is found to be the determinant.
c3=(-dda+tt*ddda)/ptmtt
cs=(-dda+pt*ddda)/dett
c1= da  +c3*pt**2/2 - cs*c2oe
at=c1*pt-c3*pt**3/6 + cs*s2oe
ae=at+aemat
! Rescale all coefficients by ae so that parameter a becomes one at the edge:
ae=u1/ae
c1=c1*ae
c3=c3*ae
cs=cs*ae
at=at*ae
end subroutine mobnetcoefs

!==============================================================================
subroutine eatomn(alpha,ae,c1,c3,cs,xea,xmn,xmnd)!                     [eatomn]
!==============================================================================
! Equiangular (EA) to Mobius-net (MN) gnomonic map coordinates
!==============================================================================
real(dp),               intent(in ):: alpha,ae,c1,c3,cs
real(dp),dimension(2)  ,intent(in ):: xea
real(dp),dimension(2)  ,intent(out):: xmn
real(dp),dimension(2,2),intent(out):: xmnd
!------------------------------------------------------------------------------
real(dp),parameter:: piq=pi/4
real(dp)          :: phi,a,da,dda
integer           :: i
!==============================================================================
xmnd=0
do i=1,2
   phi=xea(i)*piq
   call phitoa(alpha,ae,c1,c3,cs,phi,a,da,dda)
   xmn(i)=a
   xmnd(i,i)=da*piq
enddo
end subroutine eatomn

!==============================================================================
subroutine mntoea(alpha,at,ae,c1,c3,cs,xmn,xea,xead)!                  [mntoea]
!==============================================================================
! Mobius-net (MN) to equiangular (EA) gnomonic map coordinates
!==============================================================================
real(dp),               intent(in ):: alpha,at,ae,c1,c3,cs
real(dp),dimension(2)  ,intent(in ):: xmn
real(dp),dimension(2)  ,intent(out):: xea
real(dp),dimension(2,2),intent(out):: xead
!-----------------------------------------------------------------------------
real(dp),parameter:: u4opi=u4/pi
real(dp)          :: a,phi,dphi,ddphi
integer           :: i
!==============================================================================
xead=0
do i=1,2
   a=xmn(i)
   call atophi(alpha,at,ae,c1,c3,cs, a, phi,dphi,ddphi)
   xea(i)=phi*u4opi
   xead(i,i)=dphi*u4opi
enddo
end subroutine mntoea

!=============================================================================
subroutine edtoea(xed,xea,xead)!                                      [edtoea]
!=============================================================================
! Edge-type (ED) to equiangular (EA) gnomonic map coordinates
!=============================================================================
use pietc, only: u1,r2,or2,pi
real(dp),dimension(2),  intent(in ):: xed
real(dp),dimension(2),  intent(out):: xea
real(dp),dimension(2,2),intent(out):: xead
!-----------------------------------------------------------------------------
real(dp),parameter:: u4opi=u4/pi
real(dp)          :: rcorn,ted
integer           :: i
!=============================================================================
rcorn=atan(or2)! Radians corner latitude (= 35.26 degrees)
xead(1,2)=0; xead(2,1)=0
do i=1,2
   ted=tan(rcorn*xed(i))
   xea(i)=u4opi*atan(r2*ted)
   xead(i,i)=u4opi*r2*rcorn*(u1+ted**2)/(u1+2*ted**2)
enddo
end subroutine edtoea

!=============================================================================
subroutine eatoed(xea,xed,xedd)!                                      [eatoed]
!=============================================================================
! Equiangular (EA) to Edge-type (ED) gnomonic map coordinates
!=============================================================================
use pietc, only: u1,r2,or2,pi
real(dp),dimension(2),  intent(in ):: xea
real(dp),dimension(2),  intent(out):: xed
real(dp),dimension(2,2),intent(out):: xedd
!-----------------------------------------------------------------------------
real(dp),parameter:: piq=pi/4
real(dp)          :: rcorn,tea
integer           :: i
!=============================================================================
rcorn=atan(or2)! Radians corner latitude (= 35.26 degrees)
xedd(1,2)=0; xedd(2,1)=0
do i=1,2
   tea=tan(piq*xea(i))
   xed(i)=atan(tea*or2)/rcorn
   xedd(i,i)=piq*or2*(u1+tea**2)/(rcorn*(u1+tea**2/2))
enddo
end subroutine eatoed

!=============================================================================
subroutine phitoa(alpha,ae,c1,c3,cs,phi,a,da,dda)!                    [phitoa]
!=============================================================================
! Given Mobius zone width, alpha (in radians), and the tensioned-spline
! interpolation parameters, at,ae,c1,c3,cs, and given a particular value
! angle phi from the tile-median, evaluate the corresponding Mobius-net
! index function, a, and its derivatives wrt phi.
!=============================================================================
real(dp),intent(in ):: alpha,ae,c1,c3,cs, phi
real(dp),intent(out):: a,da,dda
!------------------------------------------------------------------------------
real(dp),parameter:: piq=pi/4
real(dp)          :: p,pe,pt,pp,aphi,phic2,s2r,c2r,c2oe,s2oe,tt,oalpha
!==============================================================================
if(alpha==0)then
   a  =phi/piq
   da =u1/piq
   dda=0
   return
endif
oalpha=u1/alpha
aphi=abs(phi)
p=aphi*oalpha
pe=piq*oalpha
pt=pe-u1
if(p>=pt)then
   a=u1+ae*log(tan(aphi))
   phic2=(piq-aphi)*2; s2r=sin(phic2); c2r=cos(phic2)
   da  =ae*u2/c2r
   dda=-ae*u4*s2r/c2r**2
else
   c2oe=exp((p-pe))+exp(-(p+pe)); tt=tanh(p); s2oe=tt*c2oe
   pp=p*p
   a  =(c1*p-c3*p*pp/6 +cs*s2oe)
   da =(c1  -c3*pp/2   +cs*c2oe)*oalpha
   dda=(    -c3*p      +cs*s2oe)*oalpha**2
endif
if(phi*aphi<0)then
   a=-a
   dda=-dda
endif
end subroutine phitoa

!==============================================================================
subroutine atophi(alpha,at,ae,c1,c3,cs, a, phi,dphi,ddphi)!            [atophi]
!==============================================================================
! Given Mobius zone width, alpha (in radians), and the tensioned-spline
! interpolation parameters, at,ae,c1,c3,cs, and given a particular value
! of a, evaluate the corresponding phi, and its derivatives wrt a. 
!==============================================================================
real(dp),intent(in ):: alpha,at,ae,c1,c3,cs, a
real(dp),intent(out):: phi,dphi,ddphi
!------------------------------------------------------------------------------
integer, parameter:: nit=20 ! <- Newton iterations allowance
real(dp),parameter:: piq=pi/4,eps=1.e-12! <- Newton convergence criterion
real(dp)          :: p,pe,da,pt,tt,c2oe,s2oe,ac,e,ee,eep,eec,ap,aa,pp
integer           :: i
!==============================================================================
if(alpha==0)then
   phi  =a*piq
   dphi =piq
   ddphi=0
   return
endif
pe=piq/alpha ! Edge, in p units
pt=pe-u1     ! Transition point, in p units
if(a==0)then
   phi=0;    dphi=alpha/(c1+cs*2*exp(-pe)); ddphi=0; return
endif
aa=abs(a)
if(aa>=at)then
   ac=(u1-aa)/ae
   e=exp(-ac)
   ee=e*e
   eep=ae*(u1+ee)
   eec=u1-ee
   phi  =atan(e)
   dphi =e/eep
   ddphi=e*eec/eep**2
else
   p=aa*pt/at
   ap=2*eps
   do i=1,nit
      c2oe=exp((p-pe))+exp(-(p+pe)); tt=tanh(p); s2oe=tt*c2oe
      pp=p*p
      da=c1-c3*pp/2+cs*c2oe ! da/dp
      if(abs(ap)<eps)exit
      ap=c1*p-c3*p*pp/6 +cs*s2oe-aa
      p=p-ap/da
   enddo
   if(i>nit)print'("In atophi; WARNING: Newton iterations not converging")'
   phi  = alpha*p
   dphi = alpha/da
   ddphi=-alpha*(-c3*p+cs*s2oe)/da**3
endif
if(a*aa<0)then ! Flip some signs for the negative angle cases:
   phi  =-phi
   ddphi=-ddphi
endif
end subroutine atophi

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
subroutine recpq(pa,qa,recpa,recqa)!                                   [recpq]
!=============================================================================
! For a relative permutation, pa, and rotation quaternary digits, qa,
! return the corresponding codes for the reciprocal relationship,
! recpa, recqa.
!=============================================================================
integer,dimension(m),  intent(in ):: pa,qa
integer,dimension(m),  intent(out):: recpa,recqa
!-----------------------------------------------------------------------------
integer:: i,p
!=============================================================================
do i=1,m
   p=pa(i)
   recpa(p)=i
   recqa(p)=modulo(-qa(i),4)
enddo
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

end module pgn




