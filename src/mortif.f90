! From https://github.com/szaghi/MORTIF/blob/master/src/lib/mortif.f90
!
module mortif
!-----------------------------------------------------------------------------------------------------------------------------------
! MORTIF, MORTon Indexer (Z-order) Fortran environment.
!
!  A library to encode/decode integer indexes into Morton's (Z-order) ordering.
!  Morton's code (Z-order) is a scheme to map multi-dimensional arrays onto to a linear with a great deal of spatial locality.
!
! #### References
!  [1] *A Computer Oriented Geodetic Data Base and a New Technique in File Sequencing*, Morton G.M., technical report, IBM, 1966.
!  [2] *On Spatial Orders and Location Codes*, Stocco, LJ and Schrack, G, IEEE Transaction on Computers, vol 58, n 3, March 2009.
!  [3] *Out-of-Core Construction of Sparse Voxel Octrees*, J. Baert, A. Lagae and Ph. Dutre, Proceedings of the Fifth ACM
!  SIGGRAPH/Eurographics conference on High-Performance Graphics, 2013.
!-----------------------------------------------------------------------------------------------------------------------------------

  use iso_fortran_env

  implicit none


  integer, parameter :: I8P = INT64
  integer, parameter :: I4P = INT32
  integer, parameter :: I2P = INT16
  integer, parameter :: I1P = INT8

  integer, parameter :: R4P = REAL32

  save
  private
  public :: morton2D, demorton2D, morton3D, demorton3D
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Binary masks: used into the bits dilating and contracting algorithms.
!> 0000000000000000000000000000000011111111111111111111111111111111.
integer(I8P) , parameter :: mask32_32=int(Z'FFFFFFFF',         I8P)
!> 0000000000000000000000000000000000000000000000001111111111111111.
integer(I8P) , parameter :: mask16_48=int(Z'FFFF',             I8P)
!> 1111111111111111000000000000000000000000000000001111111111111111.
integer(I8P) , parameter :: mask16_32=int(Z'FFFF00000000FFFF', I8P)
!> 0000000000000000111111111111111100000000000000001111111111111111.
integer(I8P) , parameter :: mask16_16=int(Z'FFFF0000FFFF',     I8P)
!> 0000000000000000000000000000000000000000000000000000000011111111.
integer(I8P) , parameter :: mask8_56 =int(Z'FF',               I8P)
!> 0000000011111111000000000000000011111111000000000000000011111111.
integer(I8P) , parameter :: mask8_16 =int(Z'FF0000FF0000FF',   I8P)
!> 0000000011111111000000001111111100000000111111110000000011111111.
integer(I8P) , parameter :: mask8_8  =int(Z'FF00FF00FF00FF',   I8P)
!> 0000000000000000000000000000000000000000000000000000000000001111.
integer(I8P) , parameter :: mask4_60 =int(Z'F',                I8P)
!> 1111000000001111000000001111000000001111000000001111000000001111.
integer(I8P) , parameter :: mask4_8  =int(Z'F00F00F00F00F00F', I8P)
!> 0000111100001111000011110000111100001111000011110000111100001111.
integer(I8P) , parameter :: mask4_4  =int(Z'F0F0F0F0F0F0F0F',  I8P)
!> 0000000000000000000000000000000000000000000000000000000000000011.
integer(I8P) , parameter :: mask2_62 =int(Z'3',                I8P)
!> 0011000011000011000011000011000011000011000011000011000011000011.
integer(I8P) , parameter :: mask2_4  =int(z'30C30C30C30C30C3', I8P)
!> 0011001100110011001100110011001100110011001100110011001100110011.
integer(I8P) , parameter :: mask2_2  =int(Z'3333333333333333', I8P)
!> 1001001001001001001001001001001001001001001001001001001001001001.
integer(I8P) , parameter :: mask1_2  =int(Z'9249249249249249', I8P)
!> 0101010101010101010101010101010101010101010101010101010101010101.
integer(I8P) , parameter :: mask1_1  =int(Z'5555555555555555', I8P)

integer(I8P), parameter :: signif(1:5) = [mask2_62,  &
                                          mask4_60,  &
                                          mask8_56,  &
                                          mask16_48, &
                                          mask32_32] !< Binary mask for selecting significant bits.

integer(I8P), parameter :: mask(1:6,1:2) = reshape([mask1_1,mask2_2,mask4_4,mask8_8, mask16_16,mask32_32,  & ! 1 bit interleaving.
                                                    mask1_2,mask2_4,mask4_8,mask8_16,mask16_32,mask32_32], & ! 2 bits interleaving.
                                                   [6,2]) !< Binary mask for perfoming significant bits shifting.

integer(I1P), parameter :: shft(1:6,1:2) = reshape([1_I1P,2_I1P,4_I1P,8_I1P, 16_I1P,32_I1P,  & ! 1 bit interleaving.
                                                    2_I1P,4_I1P,8_I1P,16_I1P,32_I1P,64_I1P], & ! 2 bits interleaving.
                                                   [6,2]) !< Shift number array.

real(R4P), parameter :: log10_2_inv = 1._R4P/log10(2._R4P) !< Real parameter for computing the number of shifts (Ns).
!<
!< The number of shifts is computed by \(2^{Ns}=b\) where `b` is the significant bits. As a consequence
!< \(Ns=\frac{log(b)}{log2}\) therefore it is convenient to store the value of \(\frac{1}{log2}\).
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function dilatate(i, b, z) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Dilatate integer of 32 bits to integer of 64 bits.
  !<
  !< See *On Spatial Orders and Location Codes*, Stocco, LJ and Schrack, G, IEEE Transaction on Computers, vol 58, n 3, March 2009.
  !< The resulting integer has 64 bits; it has only `b` significant bits interleaved by `z` zeros: `bb/zx0/bb-1/zx0../b1/zx0/b0`;
  !< e.g. for `(b=4, z=1)`: `b3/b2/b1/b0 => b3/0/b2/0/b1/0/b0`.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in) :: i !< Input integer.
  integer(I2P), intent(in) :: b !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I1P), intent(in) :: z !< Number of zero 'i' (1/2).
  integer(I8P)             :: d !< Dilated integer.
  integer(I1P)             :: l !< Counter.
  integer(I1P)             :: m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = int(log10(b*1._R4P)*log10_2_inv, I1P)
  d = int(i, I8P)
  d = iand(d, signif(l))
  do m=l, 1_I1P, -1_I1P !(5/4/3/2/1,1,-1)
    d = iand(ior(d, ishft(d, shft(m,z))), mask(m,z))
  enddo
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dilatate

  elemental subroutine contract(i, b, z, c)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Contract integer of 64 bits into integer of 32 bits.
  !<
  !< See *On Spatial Orders and Location Codes*, Stocco, LJ and Schrack, G, IEEE Transaction on Computers, vol 58, n 3, March 2009.
  !< The resulting integer(int8/int16/int32) has only `b' significant bits obtained by the following contraction:
  !< if `bb/zx0/bb-1/zx0../b1/zx0/b0 => bb/bb-1/.../b1/b0`; e.g. for `(b=4,z=1)`: `b3/0/b2/0/b1/0/b0 => b3/b2/b1/b0`.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)  :: i !< Input integer.
  integer(I2P), intent(in)  :: b !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I1P), intent(in)  :: z !< Number of zero 'i' (1/2).
  integer(I4P), intent(out) :: c !< Contracted integer.
  integer(I8P)              :: d !< Temporary dilated integer.
  integer(I1P)              :: m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = iand(i, mask(1,z))
  do m=1_I1P, int(log10(b*1._R4P)*log10_2_inv, I1P), 1_I1P !(1,5/4/3/2/1)
    d = iand(ieor(d, ishft(d, -shft(m,z))), mask(m+1,z))
  enddo
  c = int(d, I4P)
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine contract

  elemental function morton2D(i, j, b) result(code)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Encode 2 integer (32 bits) indexes into 1 integer (64 bits) Morton's code.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)           :: i    !< I index.
  integer(I4P), intent(in)           :: j    !< J index.
  integer(I2P), intent(in), optional :: b    !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I8P)                       :: code !< Morton's code.
  integer(I8P)                       :: di   !< Dilated indexe.
  integer(I8P)                       :: dj   !< Dilated indexe.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dilatate(i=i, b=b, z=1_I1P)
    dj = dilatate(i=j, b=b, z=1_I1P)
  else
    di = dilatate(i=i, b=32_I2P, z=1_I1P)
    dj = dilatate(i=j, b=32_I2P, z=1_I1P)
  endif
  code = ishft(dj,1) + di
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton2D

  elemental subroutine demorton2D(code, i, j, b)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Decode 1 integer (64 bits) Morton's code into 2 integer (32 bits) indexes.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)           :: code !< Morton's code.
  integer(I4P), intent(inout)        :: i    !< I index.
  integer(I4P), intent(inout)        :: j    !< J index.
  integer(I2P), intent(in), optional :: b    !< Number of significant bits of 'i' (2/4/8/16/32).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    call contract(i=code,           b=b, z=1_I1P, c=i)
    call contract(i=ishft(code,-1), b=b, z=1_I1P, c=j)
  else
    call contract(i=code,           b=32_I2P, z=1_I1P, c=i)
    call contract(i=ishft(code,-1), b=32_I2P, z=1_I1P, c=j)
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton2D

  elemental function morton3D(i, j, k, b) result(code)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Encode 3 integer (32 bits) indexes into 1 integer (64 bits) Morton's code.
  !<
  !< @note Due to 64 bits limit of the Morton's code, the 3 allowed-side of indexes is limited to 21 bits.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(in)           :: i    !< I index.
  integer(I4P), intent(in)           :: j    !< J index.
  integer(I4P), intent(in)           :: k    !< K index.
  integer(I2P), intent(in), optional :: b    !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I8P)                       :: code !< Morton's code.
  integer(I8P)                       :: di   !< Dilated indexes.
  integer(I8P)                       :: dj   !< Dilated indexes.
  integer(I8P)                       :: dk   !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dilatate(i=i, b=b, z=2_I1P)
    dj = dilatate(i=j, b=b, z=2_I1P)
    dk = dilatate(i=k, b=b, z=2_I1P)
  else
    di = dilatate(i=i, b=32_I2P, z=2_I1P)
    dj = dilatate(i=j, b=32_I2P, z=2_I1P)
    dk = dilatate(i=k, b=32_I2P, z=2_I1P)
  endif
  code = ishft(dk,2) + ishft(dj,1) + di
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton3D

  elemental subroutine demorton3D(code, i, j, k, b)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Decode 1 integer (64 bits) Morton's code into 3 integer (16 bits) indexes.
  !<
  !< @note Due to 64 bits limit of the Morton's code, the 3 allowed-side of indexes is limited to 21 bits.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(in)           :: code !< Morton's code.
  integer(I4P), intent(inout)        :: i    !< I index.
  integer(I4P), intent(inout)        :: j    !< J index.
  integer(I4P), intent(inout)        :: k    !< K index.
  integer(I2P), intent(in), optional :: b    !< Number of significant bits of 'i' (2/4/8/16).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    call contract(i=code,           b=b, z=2_I1P, c=i)
    call contract(i=ishft(code,-1), b=b, z=2_I1P, c=j)
    call contract(i=ishft(code,-2), b=b, z=2_I1P, c=k)
  else
    call contract(i=code,           b=32_I2P, z=2_I1P, c=i)
    call contract(i=ishft(code,-1), b=32_I2P, z=2_I1P, c=j)
    call contract(i=ishft(code,-2), b=32_I2P, z=2_I1P, c=k)
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton3D
endmodule mortif
