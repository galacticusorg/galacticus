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

!!{
Contains a program to test tensor functionality.
!!}

program Test_Tensors
  !!{
  Tests of coordinate system functions.
  !!}
  use :: Display   , only : displayVerbositySet           , verbosityLevelStandard
  use :: Tensors   , only : assignment(=)                 , operator(*)           , tensorIdentityR2D3Sym, tensorNullR2D3Sym, &
          &                 tensorRank2Dimension3Symmetric
  use :: Unit_Tests, only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group , Unit_Tests_Finish
  type(tensorRank2Dimension3Symmetric) :: tensorR2D3Sym,resultTensor

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
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Tensors
