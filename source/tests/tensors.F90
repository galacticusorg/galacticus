!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{RST
Contains a program to test tensor functionality.
!!}

program Test_Tensors
  !!{RST
  Tests of coordinate system functions.
  !!}
  use :: Display          , only : displayVerbositySet           , verbosityLevelStandard
  use :: Tensors          , only : assignment(=)                 , operator(*)                   , max                  , tensorIdentityR2D3Sym, &
          &                        tensorNullR2D3Sym             , tensorRank2Dimension3Symmetric                       , tensorUnitR2D3Sym
  use :: Unit_Tests       , only : Assert                        , Unit_Tests_Begin_Group        , Unit_Tests_End_Group , Unit_Tests_Finish
  type            (tensorRank2Dimension3Symmetric) :: tensorR2D3Sym,resultTensor,otherTensor
  double precision                                 :: serialized(6)
  integer                                          :: fileHandle

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Tensors")

  ! Test rank 2, 3D, symmetric tensors.
  call Unit_Tests_Begin_Group("Rank 2, 3D, symmetric tensors")
  tensorR2D3Sym=reshape([1.0d0,2.0d0,3.0d0,2.0d0,4.0d0,5.0d0,3.0d0,5.0d0,6.0d0],[3,3])
  resultTensor =tensorR2D3Sym+tensorR2D3Sym
  call Assert(                                                                           &
       &      'addition'                                                               , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([2.0d0,4.0d0,6.0d0,4.0d0,8.0d0,10.0d0,6.0d0,10.0d0,12.0d0],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  resultTensor =tensorR2D3Sym-tensorR2D3Sym
  call Assert(                                                                           &
       &      'subtraction'                                                            , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  resultTensor =2.0d0*tensorR2D3Sym
  call Assert(                                                                           &
       &      'scalar multiplication (scalar first)'                                   , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([2.0d0,4.0d0,6.0d0,4.0d0,8.0d0,10.0d0,6.0d0,10.0d0,12.0d0],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  resultTensor =tensorR2D3Sym*2.0d0
  call Assert(                                                                           &
       &      'scalar multiplication (tensor first)'                                   , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([2.0d0,4.0d0,6.0d0,4.0d0,8.0d0,10.0d0,6.0d0,10.0d0,12.0d0],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  resultTensor =tensorR2D3Sym/2.0d0
  call Assert(                                                                           &
       &      'scalar division'                                                        , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([0.5d0,1.0d0,1.5d0,1.0d0,2.0d0,2.5d0,1.5d0,2.5d0,3.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call resultTensor%reset()
  call Assert(                                                                           &
       &      'set to zero'                                                            , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call resultTensor%setToUnity()
  call Assert(                                                                           &
       &      'set to unity'                                                           , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call resultTensor%setToIdentity()
  call Assert(                                                                           &
       &      'set to identity'                                                        , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call Assert(                                                                           &
       &      'supplied identity'                                                      , &
       &      tensorIdentityR2D3Sym%toMatrix()                                         , &
       &      reshape([1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call Assert(                                                                           &
       &      'supplied zero'                                                          , &
       &      tensorNullR2D3Sym%toMatrix()                                             , &
       &      reshape([0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  call Assert(                                                                           &
       &      'contraction'                                                            , &
       &      tensorR2D3Sym%contract()                                                 , &
       &      11.0d0                                                                   , &
       &      absTol=1.0d-6                                                              &
       &     )
  call Assert(                                                                           &
       &      'double contraction'                                                     , &
       &      tensorR2D3Sym%doubleContract(tensorR2D3Sym)                              , &
       &      129.0d0                                                                  , &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Construction from six components (public interface) and conversion to a matrix.
  resultTensor=tensorRank2Dimension3Symmetric(1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0)
  call Assert(                                                                           &
       &      'construct from components'                                              , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([1.0d0,2.0d0,3.0d0,2.0d0,4.0d0,5.0d0,3.0d0,5.0d0,6.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Null constructor sets all components to zero.
  resultTensor=tensorRank2Dimension3Symmetric()
  call Assert(                                                                           &
       &      'null constructor'                                                       , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Supplied unit tensor (all elements unity).
  call Assert(                                                                           &
       &      'supplied unit'                                                          , &
       &      tensorUnitR2D3Sym%toMatrix()                                             , &
       &      reshape([1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Element accessor, using the (i,j) enumeration (0-based indices), including the
  ! symmetric off-diagonal reflection.
  call Assert('element (0,0)'           ,tensorR2D3Sym%element(0,0),1.0d0,absTol=1.0d-6)
  call Assert('element (0,1)'           ,tensorR2D3Sym%element(0,1),2.0d0,absTol=1.0d-6)
  call Assert('element (0,2)'           ,tensorR2D3Sym%element(0,2),3.0d0,absTol=1.0d-6)
  call Assert('element (1,1)'           ,tensorR2D3Sym%element(1,1),4.0d0,absTol=1.0d-6)
  call Assert('element (1,2)'           ,tensorR2D3Sym%element(1,2),5.0d0,absTol=1.0d-6)
  call Assert('element (2,2)'           ,tensorR2D3Sym%element(2,2),6.0d0,absTol=1.0d-6)
  call Assert('element (2,0) [symmetry]',tensorR2D3Sym%element(2,0),3.0d0,absTol=1.0d-6)
  call Assert('element (2,1) [symmetry]',tensorR2D3Sym%element(2,1),5.0d0,absTol=1.0d-6)
  ! isZero predicate.
  call Assert('isZero (non-zero tensor)',tensorR2D3Sym    %isZero(),.false.)
  call Assert('isZero (null tensor)'    ,tensorNullR2D3Sym%isZero(),.true. )
  ! Increment (in-place addition).
  resultTensor=tensorR2D3Sym
  call resultTensor%increment(tensorR2D3Sym)
  call Assert(                                                                           &
       &      'increment'                                                              , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([2.0d0,4.0d0,6.0d0,4.0d0,8.0d0,10.0d0,6.0d0,10.0d0,12.0d0],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Element-by-element max().
  otherTensor =tensorRank2Dimension3Symmetric(0.0d0,3.0d0,0.0d0,5.0d0,0.0d0,7.0d0)
  resultTensor=max(tensorR2D3Sym,otherTensor)
  call Assert(                                                                           &
       &      'element-wise max'                                                       , &
       &      resultTensor%toMatrix()                                                  , &
       &      reshape([1.0d0,3.0d0,3.0d0,3.0d0,5.0d0,5.0d0,3.0d0,5.0d0,7.0d0   ],[3,3]), &
       &      absTol=1.0d-6                                                              &
       &     )
  ! Tensor/matrix equality operator.
  call Assert(                                                                                         &
       &      'matrix equality (equal)'                                                              , &
       &      tensorR2D3Sym == reshape([1.0d0,2.0d0,3.0d0,2.0d0,4.0d0,5.0d0,3.0d0,5.0d0,6.0d0],[3,3]), &
       &      .true.                                                                                   &
       &     )
  call Assert(                                                                                         &
       &      'matrix equality (unequal)'                                                            , &
       &      tensorR2D3Sym == reshape([1.0d0,2.0d0,3.0d0,2.0d0,4.0d0,5.0d0,3.0d0,5.0d0,9.0d0],[3,3]), &
       &      .false.                                                                                  &
       &     )
  ! Vector projection: x . A . x with x=(1,1,1) is the sum of all matrix elements.
  call Assert(                                                                                 &
       &      'vector projection'                                                            , &
       &      tensorR2D3Sym%vectorProject([1.0d0,1.0d0,1.0d0])                               , &
       &      31.0d0                                                                         , &
       &      absTol=1.0d-6                                                                    &
       &     )
  ! Serialize to an array (packed upper triangle).
  call tensorR2D3Sym%serialize(serialized)
  call Assert(                                                                                 &
       &      'serialize'                                                                    , &
       &      serialized                                                                     , &
       &      [1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0]                                          , &
       &      absTol=1.0d-6                                                                    &
       &     )
  ! Deserialize from an array (round-trips through the same packing).
  call resultTensor%deserialize([10.0d0,20.0d0,30.0d0,40.0d0,50.0d0,60.0d0])
  call Assert(                                                                                 &
       &      'deserialize'                                                                  , &
       &      resultTensor%toMatrix()                                                        , &
       &      reshape([10.0d0,20.0d0,30.0d0,20.0d0,40.0d0,50.0d0,30.0d0,50.0d0,60.0d0],[3,3]), &
       &      absTol=1.0d-6                                                                    &
       &     )
  ! Serialized property count.
  call Assert('serialize count',tensorR2D3Sym%serializeCount(),6)
  ! Binary dump/read round-trip pins the on-disk storage contract.
  open(newunit=fileHandle,status='scratch',form='unformatted',action='readwrite')
  call tensorR2D3Sym%dumpRaw(fileHandle)
  rewind(fileHandle)
  call resultTensor%readRaw(fileHandle)
  close(fileHandle)
  call Assert(                                                                                &
       &      'binary dump/read round-trip'                                                 , &
       &      resultTensor%toMatrix()                                                       , &
       &      tensorR2D3Sym%toMatrix()                                                      , &
       &      absTol=1.0d-6                                                                   &
       &     )
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Tensors
