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
  Implements a supernovae type Ia class which operates through the differential delay time distribution.
  !!}

  !![
  <supernovaeTypeIa name="supernovaeTypeIaDifferentialDTD" abstract="yes">
   <description>
    A supernovae type Ia class which operates through the differential delay time distribution.
   </description>
  </supernovaeTypeIa>
  !!]
  type, abstract, extends(supernovaeTypeIaMassIndependentDTD) :: supernovaeTypeIaDifferentialDTD
     !!{
     A supernovae type Ia class which operates through the differential delay time distribution.
     !!}
     private
   contains
     !![
     <methods>
       <method method="numberDifferential" description="Return the number of type Ia supernovae per Solar mass of stars formed per unit time at a given population age and metallicity."/>
       <method method="timeDelayMinimum"   description="Return the minimum time in the type Ia supernovae delay time distribution."                                                     />
       <method method="timeDelayMaximum"   description="Return the maximum time in the type Ia supernovae delay time distribution."                                                     />
     </methods>
     !!]
     procedure                                       :: numberCumulative   => differentialDTDNumberCumulative
     procedure(numberDifferentialTemplate), deferred :: numberDifferential
     procedure                                       :: timeDelayMinimum   => differentialDTTimeDelayMinimum
     procedure                                       :: timeDelayMaximum   => differentialDTTimeDelayMaximum
  end type supernovaeTypeIaDifferentialDTD

  abstract interface
     double precision function numberDifferentialTemplate(self,age,metallicity)
       !!{
       Interface for differential number of Type Ia SNe.
       !!}
       import supernovaeTypeIaDifferentialDTD
       class           (supernovaeTypeIaDifferentialDTD), intent(inout) :: self
       double precision                                 , intent(in   ) :: age , metallicity
     end function numberDifferentialTemplate
  end interface

  ! Sub-module-scope variables used in integration.
  class           (supernovaeTypeIaDifferentialDTD), pointer :: self_
  double precision                                           :: metallicity_
  !$omp threadprivate(self_,metallicity_)
  
contains
  
  double precision function differentialDTTimeDelayMinimum(self) result(time)
    !!{
    Compute the minimum time in the type Ia supernovae delay time distribution.
    !!}
    implicit none
    class(supernovaeTypeIaDifferentialDTD), intent(inout) :: self
    !$GLC attributes unused :: self

    time=0.0d0
    return
  end function differentialDTTimeDelayMinimum
  
  double precision function differentialDTTimeDelayMaximum(self) result(time)
    !!{
    Compute the maximum time in the type Ia supernovae delay time distribution.
    !!}
    implicit none
    class(supernovaeTypeIaDifferentialDTD), intent(inout) :: self
    !$GLC attributes unused :: self

    time=huge(0.0d0)
    return
  end function differentialDTTimeDelayMaximum
  
  double precision function differentialDTDNumberCumulative(self,age,metallicity) result(number)
    !!{
    Compute the cumulative number of Type Ia supernovae originating per unit interval of secondary star mass with given
    {\normalfont \ttfamily initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. Here we
    assume that the total number of Type Ias is specified independent of secondary star mass.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (supernovaeTypeIaDifferentialDTD), intent(inout), target :: self
    double precision                                 , intent(in   )         :: age        , metallicity
    type            (integrator                     ), allocatable  , save   :: integrator_
    !$omp threadprivate(integrator_)

    if (age < self%timeDelayMinimum()) then
       ! The minimum time in the delay time distribution has not yet been reached.
       number=0.0d0
    else
       ! Integrate to find the cumulative number of type Ia supernovae at this time.
       if (.not.allocated(integrator_)) then
          allocate(integrator_)
          integrator_=integrator(timeDelayDistributionDifferential,toleranceRelative=1.0d-3)
       end if
       self_ => self
       metallicity_=metallicity
       number=integrator_%integrate(self%timeDelayMinimum(),min(age,self%timeDelayMaximum()))       
    end if
    return
  end function differentialDTDNumberCumulative
  
  double precision function timeDelayDistributionDifferential(time)
    !!{
    Integrand used to compute the cumulative number of type Ia supernovae occuring by a given time.
    !!}
    implicit none
    double precision, intent(in   ) :: time
    
    timeDelayDistributionDifferential=self_%numberDifferential(time,metallicity_)
    return
  end function timeDelayDistributionDifferential

