!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  A primordial power spectrum class which provides a piecewise power-law power spectrum.
  !!}

  !![
  <powerSpectrumPrimordial name="powerSpectrumPrimordialPiecewisePowerLaw">
   <description>
    Implements a piecewise power-law primordial power spectrum, possibly with a running index. The primordial power spectrum has the
    form:
    \begin{equation}
     P(k) \propto A_\mathrm{i} k^{n_\mathrm{eff,i}(k)},
    \end{equation}
    where
    \begin{equation}
     n_\mathrm{eff,i}(k) = n_\mathrm{s,i} + {1\over 2}{\d n \over \d \ln k}_i \ln \left( {k \over k_\mathrm{ref,i}} \right) + {1\over 6}{\d^2 n \over \d \ln k^2}_i \left[ \ln \left( {k \over k_\mathrm{ref}} \right) \right]^2,
    \end{equation}
    where $n_\mathrm{s,i}=${\normalfont \ttfamily [index]} is the power spectrum index at wavenumber
    $k_\mathrm{ref,i}=${\normalfont \ttfamily [wavenumberReference]}, $\d n / \d \ln k_i=${\normalfont \ttfamily [running]}, and $\d^2 n / \d \ln k^2_i=${\normalfont \ttfamily [runningRunning]}
    describes the running of this index with wavenumber. The subscript ``i'', which runs from $1$ to $N$ refers to each interval of the piecewise power-law. Note that $k_\mathrm{ref,i}$ is defined only for $i\ge 2$. For the first ($i=1$) interval, the wavenumber ranges from $0$ to $k_\mathrm{ref,2}$.

    The amplitudes, $A_i$, are chosen to make the power spectrum continuous.
   </description>
  </powerSpectrumPrimordial>
  !!]
  type, extends(powerSpectrumPrimordialClass) :: powerSpectrumPrimordialPiecewisePowerLaw
     !!{
     A power-law primordial power spectrum class.
     !!}
     private
     double precision, allocatable, dimension(:) :: index_        , running            , &
          &                                         runningRunning, wavenumberReference, &
          &                                         normalization
   contains
     !![
     <methods>
       <method method="indexEffective" description="Compute the local effective index of the power spectrum."  />
       <method method="indices"        description="Compute the indices corresponding to the given wavenumber."/>
     </methods>
     !!]
     procedure :: power                 => piecewisePowerLawPower
     procedure :: logarithmicDerivative => piecewisePowerLawLogarithmicDerivative
     procedure :: indexEffective        => piecewisePowerLawIndexEffective
     procedure :: indices               => piecewisePowerLawIndices
  end type powerSpectrumPrimordialPiecewisePowerLaw

  interface powerSpectrumPrimordialPiecewisePowerLaw
     !!{
     Constructors for the \refClass{powerSpectrumPrimordialPiecewisePowerLaw} primordial power spectrum class.
     !!}
     module procedure piecewisePowerLawConstructorParameters
     module procedure piecewisePowerLawConstructorInternal
  end interface powerSpectrumPrimordialPiecewisePowerLaw

contains

  function piecewisePowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumPrimordialPiecewisePowerLaw} primordial power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumPrimordialPiecewisePowerLaw)                              :: self
    type            (inputParameters                         ), intent(inout)               :: parameters
    double precision                                          , allocatable  , dimension(:) :: index_ , wavenumberReference, &
         &                                                                                     running, runningRunning

    if (parameters%isPresent('index'              )) then
       allocate(index_             (parameters%count('index'              )))
    else
       allocate(index_             (2                                      ))
       end if
    if (parameters%isPresent('running'            )) then
       allocate(running            (parameters%count('running'            )))
    else
       allocate(running            (2                                      ))
       end if
    if (parameters%isPresent('runningRunning'     )) then
       allocate(runningRunning     (parameters%count('runningRunning'     )))
    else
       allocate(runningRunning     (2                                      ))
       end if
    if (parameters%isPresent('wavenumberReference')) then
       allocate(wavenumberReference(parameters%count('wavenumberReference')))
    else
       allocate(wavenumberReference(1                                      ))
       end if
    !![
    <inputParameter>
      <name>index</name>
      <variable>index_</variable>
      <source>parameters</source>
      <defaultValue>[0.9649d0,0.9649d0]</defaultValue>
      <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
      <description>The index of the power-law primordial power spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>running</name>
      <source>parameters</source>
      <defaultValue>[0.0d0,0.0d0]</defaultValue>
      <description>The running, $\d n_\mathrm{s} / \d \ln k$, of the power spectrum index.</description>
    </inputParameter>
    <inputParameter>
      <name>runningRunning</name>
      <source>parameters</source>
      <defaultValue>[0.0d0,0.0d0]</defaultValue>
      <description>The running-of-the-running, $\d^2 n_\mathrm{s} / \d \ln k^2$, of the power spectrum index.</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberReference</name>
      <source>parameters</source>
      <defaultValue>[1.0d0]</defaultValue>
      <description>When a running power spectrum index is used, this is the wavenumber, $k_\mathrm{ref}$, at which the index is equal to {\normalfont \ttfamily [index]}.</description>
    </inputParameter>
    !!]
    self=powerSpectrumPrimordialPiecewisePowerLaw(index_,running,runningRunning,wavenumberReference)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function piecewisePowerLawConstructorParameters

  function piecewisePowerLawConstructorInternal(index_,running,runningRunning,wavenumberReference) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumPrimordialPiecewisePowerLaw} primordial power spectrum class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (powerSpectrumPrimordialPiecewisePowerLaw)                              :: self
    double precision                                          , intent(in   ), dimension(:) :: index_ , wavenumberReference, &
         &                                                                                     running, runningRunning
    integer                                                                                 :: i
    !![
    <constructorAssign variables="index_, wavenumberReference, running, runningRunning"/>
    !!]

    ! Validate.
    if (size(index_) <                            2) call Error_Report("'index' must have at least two elements"                        //{introspection:location})
    if (size(index_) /= size(running            )  ) call Error_Report("'running' must have same size as 'index'"                       //{introspection:location})
    if (size(index_) /= size(runningRunning     )  ) call Error_Report("'runningRunning' must have same size as 'index'"                //{introspection:location})
    if (size(index_) /= size(wavenumberReference)+1) call Error_Report("'wavenumberReference' must have one fewer elements than 'index'"//{introspection:location})
    ! Compute normalizations.
    allocate(self%normalization(size(self%index_)))
    self%normalization(1:2)=1.0d0
    if (size(self%index_) > 2) then
       do i=3,size(self%index_)
          self%normalization(i)=+   self%normalization      (                              i-1    ) &
               &                *(                                                                  &
               &                  + self%wavenumberReference(                              i-1    ) &
               &                  / self%wavenumberReference(                                  i-2) &
               &                 )**self%indexEffective     (self%wavenumberReference(i-1),i-1,i-2)
       end do
    end if
    return
  end function piecewisePowerLawConstructorInternal

  double precision function piecewisePowerLawPower(self,wavenumber)
    !!{
    Return the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPiecewisePowerLaw), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber
    integer                                                                   :: i         , j

    call                       self%indices            (wavenumber,i,j)
    piecewisePowerLawPower=+   self%normalization      (           i  ) &
         &                 *(                                           &
         &                   +      wavenumber                          &
         &                   / self%wavenumberReference(             j) &
         &                  )**self%indexEffective     (wavenumber,i,j)
    return
  end function piecewisePowerLawPower

  double precision function piecewisePowerLawLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPiecewisePowerLaw), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber
    integer                                                                   :: i         , j

    call                                   self%indices       (wavenumber,i,j)
    piecewisePowerLawLogarithmicDerivative=self%indexEffective(wavenumber,i,j)
    return
  end function piecewisePowerLawLogarithmicDerivative

  double precision function piecewisePowerLawIndexEffective(self,wavenumber,i,j)
    !!{
    Compute the local effective index of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPiecewisePowerLaw), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber
    integer                                                   , intent(in   ) :: i         , j

    piecewisePowerLawIndexEffective=+self%index_(i)                &
         &                       +1.0d0/2.0d0                      &
         &                       *self%running(i)                  &
         &                       *log(                             &
         &                            +     wavenumber             &
         &                            /self%wavenumberReference(j) &
         &                           )                             &
         &                       +1.0d0/6.0d0                      &
         &                       *self%runningRunning(i)           &
         &                       *log(                             &
         &                            +     wavenumber             &
         &                            /self%wavenumberReference(j) &
         &                           )**2
    return
  end function piecewisePowerLawIndexEffective
  
  subroutine piecewisePowerLawIndices(self,wavenumber,i,j)
    !!{
    Compute the indices corresponding to the given wavenumber.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPiecewisePowerLaw), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber
    integer                                                   , intent(  out) :: i         , j

    if (wavenumber < self%wavenumberReference(1)) then
       i=1
       j=1
    else if (wavenumber >= self%wavenumberReference(size(self%wavenumberReference))) then
       i=size(self%index_)
       j=size(self%wavenumberReference)
    else
       j=1
       do while (wavenumber > self%wavenumberReference(j+1))
          j=j+1
       end do
       i=j+1
    end if
    return
  end subroutine piecewisePowerLawIndices

