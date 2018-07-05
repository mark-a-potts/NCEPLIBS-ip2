!
!                                      ***************************************
!                                      *      module pgn_tables            *
!                                      *      R. J. Purser                   *
!                                      *      NOAA/NCEP/EMC July 2018       *
!                                      *      jim.purser@noaa.gov            *
!                                      ***************************************
!
! Tables of coefficients and flags relating to as many as nntran distinct
! gnomonic-cube transoformations between Earth-centered normalized 
! cartesian 3-vectors, xc, and tile-indexed map coordinates, xm. 
! 
! The conditions {-1. < xm(1) < +1., -1. < xm(2) < +1.} define the interior of
! each of the six tiles indexed, ipan={1,..,6}. Each cube configuration is
! defined in relation to the orientation of what we shall call the "protoframe"
! whose own six tiles, and their accompanying right-handed (x,y) map 
! coordinates are arranged, relative to each other, such that a "development"
! of them on a single flat sheet might look like:
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
! where the arrows point in the direction of the local map's x-axis (the y-axis
! being 90-degrees counterclockwise from it, always). The orientation of this
! frame on the Earth is fixed by defining the latitude and longitude of the
! "focus" at the center of protoframe's tile-1, together with an optional
! azimuth angle by which the x-axis at this focus is twisted in the counter-
! clockwise sense relative to local East, all angles being specified in
! degrees, latitudes positive towards the North, longitudes positive towards
! the East. In the special cases where the focal latitude is at either -90
! of +90, the orientation is taken to be that which obtains when this
! polar latitude is taken as the continuous limit, the longutude and azimuth
! remaining fixed. For a rigid cube, the locations and orientations of the
! other tiles follows unambiguously from the aforementioned three angle
! specifications. But, in the case where a Schmidt conformal transformation
! is also applied, it is always applied with the focus (the center of 
! protoframe's tile-1)  as the point where the specified extremal enhancement 
! factor (the default value, 1, being no enhancement, values bigger than 1 
! being positive enhancement) applies.
!
! To accommodate the full range of possible choices of tile-topologies that
! a user might draw from, we allow for each tile of the rigidly rotated
! protoframe to be renumbered by the user through the specification of a
! general permutation of the tile indices, {1,2,3,4,5,6}. This is done by
! listing for each prototile in order, the corresponding user-defined
! tile index, such as "[6,1,2,4,5,3]", which implies that the user intends
! the center of his/her "tile-6" to be the orientation-defining focus. Such
! a permutation can always be reduced to a single integer index in [0--719]
! by simply indexing all the possible permutations in lexicographic order.
! Thus, the trivial identity permutation, [1,2,3,4,5,6] has index 0, while
! the nontrivial example [6,1,2,4,5,3] (which is what the FV3 uses) 
! corresponds to the lexicographically-indexed permutation, 603. But we 
! must also accommodate the user's choice of the orientations of each of
! the map coordinates in each of the protoframe tiles except the first
! (since that tile's absolute orientation is already defined by the
! three angles, including the azimuth angle, which makes any further twisting 
! of this tile's coordinate frame redundant. Also, since the right-handed
! local map coodinates each have only four possible orientations (multiples
! 0, 1, 2, or 3 right-angles counterclockwise relative to the standard
! local orientation of that tile of the protoframe), we can express the set of 
! five possible local coordinate orientations of the tiles 2--6 of the 
! protoframe by five digits of a base-4, or  "quaternary", number. The
! convention we adopt is to make "units" of the quaternary expansion of the
! "rotations index" correspond to the tile-6 (of the protoframe), the "fours"
! to tile-5, and so on, giving all possible options for the five nonfocal 
! tiles' coordinate orientations distilled into a single rotations index
! in the range [0--1023].
!     
! Once the orientation, optional Schmidt stretching, tile permutation
! and rotations codes have been decided, there is also the choice of which
! flavor of gnomonic grid construction to use. The simplest and most
! natural choice is the "equiangular" gnomonic which we denote using the
! abbreviation "ea". The angles subtended by the planes defining the 
! great-circle grid lines and the plane of the tile median (the "zero
! coordinate" of this local family) are exactly equally spaced and 
! continue seamlessly, beyond each cube-edge, to correspond exactly to the 
! equivalent family of grid lines of the neighboring tile, so that 
! interpolation between the points of one tile's grid and the smooth extension 
! of the grid of the neighboring tile entails only single one-dimensional
! operations in each instance. The FV3 model uses another distribution
! of angles between successive planes of each family of grid lines which
! are not equally-spaced, but which mark out points along the 
! projected edges of the cube that are uniformly-spaced on the sphere.
! We can refer to this as the "edge grid", abbreviated to "ed". There is
! another possibility that might be convenient for applications (such as
! data assimilation) where frequent accurate and smooth interpolations need
! to be performed between the overlapping extensions of the grids of
! neighboring tiles, and that is the "Mobius-net" grid, abbreviated to "mn".
! In this grid, not only are the grids line families that include a cube-edge
! as a member smoothly continuous across that edge (though not precisely 
! equally-spaced in angle), so that, like the EA grid, the interchange of data
! between neighboring tiles involves at worst only one-dimensional interpolation
! stencils, but in the vicinity of the eight cube-corners, all three families
! of grid lines shared in common by the three neighboring tiles that meet
! there, intersect in perfect three-way intersections to form what is known
! as a "Mobius-net". This, in a patch containing each given corner, there
! is no necessity for any nontrivial interpolation to be performed at all.
! Farther from the corner, outside a well-defined "Mobius zone" the grid
! lines are smoothly continued by a blending procedure based on integrated
! tensioned splines that make the transition virtually seamless, and such
! that, where the one-dimensional stencils ARE needed, there is never any
! confusion, or need for additional computations, to ensure that the source 
! points of those stencils are unambiguously defined, even when the stencils
! are relatively broad (more than four points, for example). Compared to
! the EA grid, this can make interpolations computationally simpler between
! the grids of neighboring tiles.
!
! A code is specified to determines which of these flavors of gnomonic
! grids is intended:
! kgmode=1 ==> "Equiangular" or  "EA" grid
! kgmode=2 ==> "Edge" or "ED" grid
! kgmode=3 ==> "Mobius-net" or "MN" grid.
!
! In the case of the Mobius-net grid, an additional real angular parameter
! (in degrees) specifies how far from the cube edge the transition is between
! the Mobius zone (where perfect three-way intersections occur) and the
! spline-blended remaining grid, and this same parameter also determines
! the characteristic "e-folding" scale inherent in the definition of the
! tensioned splines that peform this blending.
!
!=============================================================================
module pgn_tables
!=============================================================================
use pkind, only: dp
implicit none
integer,parameter             :: nntran=4
real(dp),dimension(nntran)    :: sschm,mg_alpha,mg_at,mg_ae,mg_c1,mg_c3,mg_cs
real(dp),dimension(3,3,nntran):: rot6
real(dp),dimension(3,3,6)     :: rotp
integer,dimension(  6,nntran) :: kgperm,kgpermi
integer,dimension(2:6,nntran) :: kgrot,kgroti
integer,dimension(nntran)     :: kgmode
logical,dimension(nntran)     :: lormode,lschm
logical                       :: linigtran
data linigtran/.false./
end module pgn_tables

