!                                             *******************************
!                                             *     module permute_cube     *
!                                             *     R. J. Purser 2018       *
!                                             *     jim.purser@noaa.gov     *
!                                             *******************************
!
!=============================================================================
! Utility routines relating permutation indices to explicit permutations
! of 6 cube tiles (faces) and 5-digit quaternary-number rotation-code to the 
! relative orientations of tile grid map coordinates (in right-angles). Also,
! code to generate the relative permutations and relative map orientation codes
! connecting pairs of cubic frames that are specified by their permutations
! and map-rotations with respect to the standard "protoframe".
!=============================================================================
module permute_cube
!=============================================================================
implicit none
integer,parameter:: m=6 ! <- Number of permutable tiles (faces) of the cube
private
public:: permofindex,indexofperm,quatofindex,indexofquat, &
         recpq,getpqawrtb,pqawrtb
interface permofindex; module procedure permofindex;        end interface
interface indexofperm; module procedure indexofperm;        end interface
interface quatofindex; module procedure quatofindex;        end interface
interface indexofquat; module procedure indexofquat;        end interface
interface recpq;       module procedure recpq;              end interface
interface getpqawrtb;  module procedure getpqawrtb;         end interface
interface pqawrtb;     module procedure pqawrtb;            end interface

contains

!=============================================================================
subroutine permofindex(i,perm)
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
if(i<0)stop'In permofindex; index must not be negative'
L=i
mixed(1)=0; perm(1)=1
do j=2,m; perm(j)=j; Lp=L/j; mixed(j)=L-Lp*j; L=Lp; enddo
if(L>0)stop'In permofindex; index is too large'
do j=1,m
   k=j+mixed(m+1-j); L=perm(k); perm(j+1:k)=perm(j:k-1); perm(j)=L
enddo
end subroutine permofindex

!=============================================================================
subroutine indexofperm(perm,i)
!=============================================================================
! Return the index, I, that this permutation, PERM, of the first M=6 integers
! would have in a lexicographically-ordering of them, where the first 
! permutation, [1,2,..] has index 0. E.g., when M=3, lexicographically ordered
! list, I going from 0 to 5, is [1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]
!=============================================================================
integer,dimension(m),intent(in ):: perm
integer,             intent(out):: i
!-----------------------------------------------------------------------------
integer:: j,jc,k,p,f
!=============================================================================
f=1; i=0
do j=1,m-1
   f=f*j! factorial
   jc=m-j; p=perm(jc); do k=jc+1,m; if(p>perm(k))i=i+f; enddo
enddo
end subroutine indexofperm

!=============================================================================
subroutine quatofindex(i,quat)
!=============================================================================
! With M=6, and freedom to rotate tiles maps, 2 -- m, of the protoframe 
! Return the quaternary (M-1)-digit, QUAT (ordered such that the units are in
! the last position, M=6) that corresponds to index I.
!=============================================================================
integer,               intent(in ):: i
integer,dimension(2:m),intent(out):: quat
!-----------------------------------------------------------------------------
integer:: ii,j
!=============================================================================
ii=i
do j=m,2,-1
   quat(j)=mod(ii,4)
   ii=ii/4
enddo
if(ii>0)stop'In quatofindex; i is too large for only (m-1) digits'
end subroutine quatofindex

!=============================================================================
subroutine indexofquat(quat,i)
!=============================================================================
! Return the index, I, that corresponds to the M quaternary digits, QUAT,
! when ordered such that the units are in position M, the 4's in position
! M-1, and so on.
!=============================================================================
integer,dimension(2:m),intent(in ):: quat
integer,               intent(out):: i
!-----------------------------------------------------------------------------
integer:: j,p
!==============================================================================
i=0
do j=2,m
   i=4*i+quat(j)
enddo
end subroutine indexofquat

!=============================================================================
subroutine recpq(pa,qa,recpa,recqa)
!=============================================================================
! For a relative ptime-permuation, pa, and 2:6 rotation quaternary digits, qa,
! return the corresponding codes for the reciprocal relationship,
! recpa, recqa.
!=============================================================================
integer,dimension(m),  intent(in ):: pa
integer,dimension(2:m),intent(in ):: qa
integer,dimension(m),  intent(out):: recpa
integer,dimension(2:m),intent(out):: recqa
!-----------------------------------------------------------------------------
integer:: i,p
!=============================================================================
do i=1,m
   p=pa(i)
   recpa(p)=i
   if(i>1)recqa(p)=modulo(-qa(i),4)
enddo
end subroutine recpq

!=============================================================================
subroutine getpqawrtb(pa,qa,pb,qb, pawrtb,qawrtb)
!=============================================================================
! Given the explicit 6-index permutation, pa, and rotation code, qa, from
! the standard protoframe, of frame A, and the corresponding permutation, pb,
! and rotation code, qb, of frame B, return the reltative permuation 
! and rotation code of A with respect to frame B.
!=============================================================================
integer,dimension(  m),intent(in ):: pa,pb
integer,dimension(2:m),intent(in ):: qa,qb
integer,dimension(  m),intent(out):: pawrtb
integer,dimension(2:m),intent(out):: qawrtb
!-----------------------------------------------------------------------------
integer:: a,b,i
!=============================================================================
do i=1,m
   a=pa(i)
   b=pb(i)
   pawrtb(b)=a
   if(i>1)qawrtb(b)=modulo(qa(i)-qb(i),4)
enddo
end subroutine getpqawrtb

!=============================================================================
subroutine pqawrtb(ib,jb,pb,pawrtb,qawrtb,ia,ja,pa,qa)
!=============================================================================
! Given the permutation, pawrtb, defining the cubic tile indices of frame A
! with respect to frame B, and the relative rotations (in right angles), qawrtb,
! of A's map ccordinates relative to the 2nd through 6th of B's map coordinates,
! then, for tile-centered integer grid coordinate, (ib,jb), in tile pb of 
! frame B, return the corresponding coordinates, (ia,ja), of this grid
! coordinates in frame A together with the tile index, pa, there, and relative
! map rotation (in right angles) of frame A with respect to frame B.
!============================================================================= 
integer,               intent(in ):: ib,jb,pb
integer,dimension(  m),intent(in ):: pawrtb
integer,dimension(2:m),intent(in ):: qawrtb
integer,               intent(out):: ia,ja,pa,qa
!=============================================================================
pa=pawrtb(pb)
qa=qawrtb(pb)
select case(qa)
case(0)
   ia= ib
   ja= jb
case(1)
   ia= jb
   ja=-ib
case(2)
   ia=-ib
   ja=-jb
case(3)
   ia=-jb
   ja= ib
end select
end subroutine pqawrtb

end module permute_cube 
