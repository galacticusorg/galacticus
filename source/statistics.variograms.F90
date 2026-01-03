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
Contains a module which implements variogram models for Gaussian Process regression.
!!}

module Statistics_Variograms
  !!{
  Implements variogram models for Gaussian Process regression.
  !!}
  private

  !![
  <functionClass>
   <name>variogram</name>
   <descriptiveName>Variograms</descriptiveName>
   <description>Class providing variogram models for Gaussian Process regression.</description>
   <default>spherical</default>
   <data>double precision :: separationNormalization, semiVarianceNormalization</data>
   <method name="fit" >
     <description>Fit the variogram model to provided data.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ), dimension(:) :: separations, semiVariances</argument>
   </method>
   <method name="fitGeneric" >
     <description>Fit a generic variogram model to provided data.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type            (enumerationVariogramFitOptionType), intent(in   )                            :: variogramFitOption               </argument>
     <argument>double precision                                   , intent(in   )             , dimension(:) :: separations       , semiVariances</argument>
     <argument>double precision                                   , intent(  out), allocatable, dimension(:) :: C                                </argument>
     <code>
       call variogramFitGeneric_(self,variogramFitOption,separations,semiVariances,C)
     </code>
   </method>
   <method name="countParameters" >
     <description>Return the number of parameters in the variogram model.</description>
     <type>integer(c_size_t)</type>
     <pass>yes</pass>
   </method>
   <method name="modelInitialGuess" >
     <description>Provide an initial guess for the parameters, $C$, of the variogram model.</description>
     <type>double precision, allocatable, dimension(:)</type>
     <argument>double precision, intent(in   ), dimension(:) :: separations, semiVariances</argument>
     <pass>yes</pass>
   </method>
   <method name="modelF" >
     <description>Evaluate the loss function, $f$, of the variogram model for the given parameters, $C$, and separations and semi-variances.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ), dimension(:) :: C, separations, semivariances</argument>
   </method>
   <method name="modelDF" >
     <description>Evaluate the gradients of the loss function, $\partial f/\partial C$, of the variogram model for the given parameters, $C$, and separations and semi-variances.</description>
     <type>double precision, allocatable, dimension(:)</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ), dimension(:) :: C, separations, semivariances</argument>
   </method>
   <method name="modelFDF" >
     <description>Evaluate the loss function, $f$, and its gradients, $\partial f/\partial C$, of the variogram model for the given parameters, $C$, and separations and semi-variances.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   )             , dimension(:) :: C   , separations, semivariances</argument>
     <argument>double precision, intent(  out)                            :: f                               </argument>
     <argument>double precision, intent(  out), allocatable, dimension(:) :: dfdC                            </argument>
   </method>
   <method name="variogram" >
     <description>Returns the variogram evaluated at the given separation. If no separation is provided the result for infinite separation is returned.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ), optional :: separation</argument>
   </method>
   <method name="correlation" >
     <description>Return the correlation evaluated at the given separation.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: separation</argument>
   </method>
  </functionClass>
  !!]

  class           (variogramClass), pointer                   :: self_
  integer                                                     :: binCount
  double precision                , allocatable, dimension(:) :: separationsBinned, semiVariancesBinned
  !$omp threadprivate(self_,separationsBinned,semiVariancesBinned,binCount)

  ! Enumeration of variogram fitting options
  !![
  <enumeration>
   <name>variogramFitOption</name>
   <description>Specifies options for fitting variograms.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="mean"   />
   <entry label="median" />
   <entry label="maximum"/>
  </enumeration>
  !!]

contains
  
  subroutine variogramFitGeneric_(self,variogramFitOption,separations,semiVariances,C)
    !!{
    Compute best fit coefficients for the variogram model.
    !!}
    use            :: Error                     , only : Error_Report
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Interface_GSL             , only : GSL_Continue   , GSL_ENoProg    , GSL_Success
    use            :: Multidimensional_Minimizer, only : multiDMinimizer
    use            :: Sorting                   , only : sortIndex
    implicit none
    class           (variogramClass                   ), intent(inout), target                     :: self
    type            (enumerationVariogramFitOptionType), intent(in   )                             :: variogramFitOption
    double precision                                   , intent(in   ), dimension(:)               :: separations                 , semiVariances
    double precision                                   , intent(  out), allocatable , dimension(:) :: C
    double precision                                                  , allocatable , dimension(:) :: separationsNormalized       , semiVariancesNormalized
    integer         (c_size_t                         )               , allocatable , dimension(:) :: rank
    integer                                            , parameter                                 :: iterationsMaximum    =10000
    double precision                                   , parameter                                 :: gradientTolerance    =1.0d-2
    double precision                                   , parameter                                 :: binWidthMaximum      =2.0d0
    integer                                            , parameter                                 :: binCountMaximum      =100000
    type            (multiDMinimizer                  )                                            :: minimizer_
    integer                                                                                        :: status                      , iteration
    integer         (c_size_t                         )                                            :: k                           , j                        , &
         &                                                                                            i
    double precision                                                                               :: currentMinimum
    logical                                                                                        :: converged
    
    ! Allocate workspace.
    self_ => self
    allocate(C                      (self%countParameters()))
    allocate(separationsNormalized  (size(separations     )))
    allocate(semiVariancesNormalized(size(separations     )))
    allocate(semiVariancesBinned    (size(separations     )))
    allocate(separationsBinned      (size(separations     )))
    allocate(rank                   (size(separations     )))
    ! Compute normalized separations and variances.
    self%separationNormalization  =sum(separations  )/dble(size(separations))
    self%semiVarianceNormalization=sum(semiVariances)/dble(size(separations))
    separationsNormalized         =separations       /self%separationNormalization
    semiVariancesNormalized       =semiVariances     /self%semiVarianceNormalization
    ! Get rank ordering by separation.
    rank=sortIndex(separationsNormalized)
    ! Compute binned estimates of the mean semi-variances.
    binCount=0
    j       =0
    k       =1
    do while (.true.)
       j=j+1
       if (j-k == binCountMaximum .or. separationsNormalized(rank(j)) > binWidthMaximum*separationsNormalized(rank(k)) .or. j == size(separations)) then
          binCount=binCount+1
          select case (variogramFitOption%ID)
          case (variogramFitOptionMean   %ID)
             separationsBinned  (binCount)=sum   (separationsNormalized  (rank(k:j)))/dble(j-k+1)
             semiVariancesBinned(binCount)=sum   (semiVariancesNormalized(rank(k:j)))/dble(j-k+1)
          case (variogramFitOptionMedian %ID)
             i=(j+k)/2
             separationsBinned  (binCount)=       separationsNormalized  (rank(  i))
             semiVariancesBinned(binCount)=       semiVariancesNormalized(rank(  i))
          case (variogramFitOptionMaximum%ID)
             i=(j+k)/2
             separationsBinned  (binCount)=       separationsNormalized  (rank(  i))
             semiVariancesBinned(binCount)=maxval(semiVariancesNormalized(rank(k:j)))
          end select
          k=j+1
          if (j == size(separations)) exit
       end if
    end do
    ! Build the minimizer.
    minimizer_=multiDMinimizer(self%countParameters(),variogramModelF,variogramModeldF,variogramModelFdF)
    ! Initial guess for the parameters.
    C=self%modelInitialGuess(separationsBinned(1:binCount),semiVariancesBinned(1:binCount))
    call minimizer_%set(C,stepSize=1.0d0,tolerance=1.0d-2)
    ! Iterate the minimizer until a sufficiently good solution is found.
    converged=.false.
    iteration=0
    do while (                                &
         &     .not.converged                 &
         &    .and.                           &
         &     iteration <  iterationsMaximum &
         &   )
       iteration=iteration+1
       call minimizer_%iterate(status)
       currentMinimum=minimizer_%minimum     (                                                  )
       converged     =minimizer_%testGradient(toleranceAbsolute=gradientTolerance*currentMinimum)
       if (status == GSL_ENoProg) exit
       if (status /= GSL_Success) call Error_Report('failed to iterate minimizer'//{introspection:location})
    end do
    ! Extract the best fit parameters.
    C=minimizer_%x()
    ! Clean up.
    deallocate(separationsNormalized  )
    deallocate(semiVariancesNormalized)
    deallocate(separationsBinned      )
    deallocate(semiVariancesBinned    )
    deallocate(rank                   )
    return
  end subroutine variogramFitGeneric_

  double precision function variogramModelF(C)
    !!{
    Function to be minimized when fitting the variogram.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: C

    variogramModelF=self_%modelF(C,separationsBinned(1:binCount),semiVariancesBinned(1:binCount))
    return
  end function variogramModelF

  function variogramModelDF(C) result(dfdC)
    !!{
    Derivatives of the function to be minimized when fitting the variogram.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: C
    double precision               , dimension(size(C)) :: dfdC

    dfdC=self_%modelDF(C,separationsBinned(1:binCount),semiVariancesBinned(1:binCount))
    return
  end function variogramModelDF

  subroutine variogramModelFDF(C,f,dfdC)
    !!{
    Computes both function and derivatives to be minimized when fitting the variogram.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: C
    double precision, intent(  out)                     :: f
    double precision, intent(  out), dimension(size(C)) :: dfdC
    double precision, allocatable  , dimension(     : ) :: dfdC_
    
    call self_%modelFdF(C,separationsBinned(1:binCount),semiVariancesBinned(1:binCount),f,dfdC_)
    dfdC=dfdC_
    return
  end subroutine variogramModelFDF
  
end module Statistics_Variograms
