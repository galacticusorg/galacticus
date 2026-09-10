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
  Implementation of a timescale for star formation which limits the timescale to be above a given multiple of the component dynamical time.
  !!}

  !![
  <starFormationTimescale name="starFormationTimescaleLowerLimited" docformat="rst">
   <description>
   A timescale for star formation which limits the timescale to be above a given multiple of the component dynamical time.
   </description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleLowerLimited
     !!{RST
     Implementation of a timescale for star formation which limits the timescale to be above a given multiple of the component dynamical time.
     !!}
     private
     class           (starFormationTimescaleClass), pointer :: starFormationTimescale_ => null()
     double precision                                       :: timescaleMinimum
     logical                                                :: diskSupported                    , spheroidSupported
   contains
     procedure :: timescale => lowerLimitedTimescale
  end type starFormationTimescaleLowerLimited

  interface starFormationTimescaleLowerLimited
     !!{RST
     Constructors for the :galacticus-class:`starFormationTimescaleLowerLimited` timescale for star formation class.
     !!}
     module procedure lowerLimitedConstructorParameters
     module procedure lowerLimitedConstructorInternal
  end interface starFormationTimescaleLowerLimited

contains

  function lowerLimitedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`starFormationTimescaleLowerLimited` timescale for star formation class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleLowerLimited)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (starFormationTimescaleClass       ), pointer       :: starFormationTimescale_
    double precision                                                    :: timescaleMinimum

    !![
    <inputParameter docformat="rst">
      <name>timescaleMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>
      The minimum timescale for star formation in units of the dynamical time.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationTimescale"  name="starFormationTimescale_"  source="parameters"/>
    !!]
    self=starFormationTimescaleLowerLimited(timescaleMinimum,starFormationTimescale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationTimescale_"/>
    !!]
    return
  end function lowerLimitedConstructorParameters

  function lowerLimitedConstructorInternal(timescaleMinimum,starFormationTimescale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`starFormationTimescaleLowerLimited` timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleLowerLimited)                        :: self
    double precision                                    , intent(in   )         :: timescaleMinimum
    class           (starFormationTimescaleClass       ), intent(in   ), target :: starFormationTimescale_
    !![
    <constructorAssign variables="timescaleMinimum, *starFormationTimescale_"/>
    <componentPropertyAssert class="disk"     properties="velocity radius" require="gettable" assignTo="self%diskSupported"    />
    <componentPropertyAssert class="spheroid" properties="velocity radius" require="gettable" assignTo="self%spheroidSupported"/>
    !!]
    return
  end function lowerLimitedConstructorInternal

  double precision function lowerLimitedTimescale(self,component)
    !!{RST
    Return a star formation rate for the given ``component`` which is limited to be no smaller than a given multiple of the dynamical time.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponent    , nodeComponentDisk, nodeComponentSpheroid
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (starFormationTimescaleLowerLimited), intent(inout) :: self
    class           (nodeComponent                     ), intent(inout) :: component
    double precision                                                    :: timeDynamical, velocity, &
         &                                                                 radius

    select type (component)
    class is (nodeComponentDisk    )
       !![
       <componentPropertyAssert class="disk"     properties="velocity radius" require="gettable" condition="self%diskSupported"    />
       !!]
       velocity=component%velocity()
       radius  =component%radius  ()
    class is (nodeComponentSpheroid)
       !![
       <componentPropertyAssert class="spheroid" properties="velocity radius" require="gettable" condition="self%spheroidSupported"/>
       !!]
       velocity=component%velocity()
       radius  =component%radius  ()
    class default
       velocity=0.0d0
       radius  =0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Check for zero velocity.
    if (velocity <= 0.0d0) then
       ! No well defined answer in this case.
       lowerLimitedTimescale=0.0d0
    else
       ! Get the dynamical time in Gyr.
       timeDynamical=+MpcPerKmPerSToGyr &
            &        *radius            &
            &        /velocity
       ! Compute the star formation timescale using a simple scaling factor.
       lowerLimitedTimescale=max(                                                           &
            &                     self%starFormationTimescale_%timescale       (component), &
            &                    +self                        %timescaleMinimum             &
            &                    *                             timeDynamical                &
            &                   )
    end if
    return
  end function lowerLimitedTimescale
