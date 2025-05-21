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
  Implementation of the \cite{baugh_can_2005} timescale for star formation.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <starFormationTimescale name="starFormationTimescaleBaugh2005">
   <description>
    A star formation timescale class which adopts the star formation rate given by a modified version of the
    \cite{baugh_can_2005} prescription:
    \begin{equation}
    \tau_\star = \tau_0 (V_\mathrm{disk}/V_0)^\alpha a^\beta
    \end{equation}
    where $\tau_0=${\normalfont \ttfamily [timescale]}, $\alpha=${\normalfont \ttfamily [exponentVelocity]},
    $\beta=${\normalfont \ttfamily [exponentExpansionFactor]}, and $V_0=${\normalfont \ttfamily [velocityNormalization]}.
   </description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleBaugh2005
     !!{
     Implementation of the \cite{baugh_can_2005} timescale for star formation.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_    => null()
     double precision                                    :: timescale_                      , exponentVelocity     , &
          &                                                 exponentExpansionFactor         , velocityNormalization
   contains
     final     ::              baugh2005Destructor
     procedure :: timescale => baugh2005Timescale
  end type starFormationTimescaleBaugh2005

  interface starFormationTimescaleBaugh2005
     !!{
     Constructors for the \refClass{starFormationTimescaleBaugh2005} timescale for star formation class.
     !!}
     module procedure baugh2005ConstructorParameters
     module procedure baugh2005ConstructorInternal
  end interface starFormationTimescaleBaugh2005

contains

  function baugh2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationTimescaleBaugh2005} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleBaugh2005)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    double precision                                                 :: timescale              , exponentVelocity     , &
         &                                                              exponentExpansionFactor, velocityNormalization

    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>8.0d0</defaultValue>
      <description>The timescale (in Gyr) for star formation in the \cite{baugh_can_2005} prescription.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <defaultValue>-3.0d0</defaultValue>
      <description>The exponent for velocity in the \cite{baugh_can_2005} prescription for star formation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentExpansionFactor</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent for expansion factor in the \cite{baugh_can_2005} prescription for star formation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityNormalization</name>
      <defaultValue>200.0d0</defaultValue>
      <description>The normalization velocity, $V_0$..</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=starFormationTimescaleBaugh2005(timescale,exponentVelocity,exponentExpansionFactor,velocityNormalization,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function baugh2005ConstructorParameters

  function baugh2005ConstructorInternal(timescale,exponentVelocity,exponentExpansionFactor,velocityNormalization,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationTimescaleBaugh2005} timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleBaugh2005)                        :: self
    double precision                                 , intent(in   )         :: timescale              , exponentVelocity     , &
         &                                                                      exponentExpansionFactor, velocityNormalization
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="exponentVelocity, exponentExpansionFactor, velocityNormalization, *cosmologyFunctions_"/>
    !!]

    self%timescale_=timescale
    return
  end function baugh2005ConstructorInternal

  subroutine baugh2005Destructor(self)
    !!{
    Destructor for the \refClass{starFormationTimescaleBaugh2005} timescale for star formation class.
    !!}
    implicit none
    type(starFormationTimescaleBaugh2005), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine baugh2005Destructor

  double precision function baugh2005Timescale(self,component)
    !!{
    Returns the timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component} in the {\normalfont
    \ttfamily baugh2005} timescale model.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponent, nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class           (starFormationTimescaleBaugh2005), intent(inout) :: self
    class           (nodeComponent                  ), intent(inout) :: component
    class           (nodeComponentBasic             ), pointer       :: basic
    double precision                                                 :: expansionFactor, velocity, &
         &                                                              time

    select type (component)
    class is (nodeComponentDisk    )
       velocity=component%velocity()
    class is (nodeComponentSpheroid)
       velocity=component%velocity()
    class default
       velocity=0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    if (velocity <= 0.0) then
       baugh2005Timescale=0.0d0
    else
       basic              =>  component%hostNode%basic                              (    )
       time               =  +basic             %time                               (    )
       expansionFactor    =  +self              %cosmologyFunctions_%expansionFactor(time)
       baugh2005Timescale =  +                                               self%timescale_             &
            &                *(velocity       /self%velocityNormalization)**self%exponentVelocity        &
            &                * expansionFactor                            **self%exponentExpansionFactor
    end if
    return
  end function baugh2005Timescale
