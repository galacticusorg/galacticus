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
Implements an N-body data operator which computes the total mass of particles.
!!}

  !![
  <nbodyOperator name="nbodyOperatorMassTotal">
   <description>An N-body data operator which computes the total mass of particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorMassTotal
     !!{
     An N-body data operator which computes the total mass of particles.
     !!}
     private
   contains
     procedure :: operate => convexHullMassTotalOperate
  end type nbodyOperatorMassTotal

  interface nbodyOperatorMassTotal
     !!{
     Constructors for the {\normalfont \ttfamily massTotal} N-body operator class.
     !!}
     module procedure convexHullMassTotalConstructorParameters
  end interface nbodyOperatorMassTotal

contains

  function convexHullMassTotalConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily convexHullMassTotal} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorMassTotal)                :: self
    type(inputParameters       ), intent(inout) :: parameters

    self=nbodyOperatorMassTotal()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function convexHullMassTotalConstructorParameters

  subroutine convexHullMassTotalOperate(self,simulations)
    !!{
    Compute the massTotal of the points.
    !!}
    use :: Error  , only : Error_Report
    use :: Display, only : displayIndent, displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorMassTotal), intent(inout)               :: self
    type            (nBodyData             ), intent(inout), dimension(:) :: simulations
    integer         (c_size_t              ), pointer      , dimension(:) :: idParticle
    double precision                        , pointer      , dimension(:) :: massParticle
    integer                                                               :: i
    double precision                                                      :: massTotal

    call displayIndent('compute total mass',verbosityLevelStandard)
    do i=1,size(simulations)
       if      (simulations(i)%propertiesReal%exists('massParticle')) then
          massParticle => simulations(i)%propertiesReal   %value('massParticle')
          massTotal=sum(massParticle)
          nullify(massParticle)
       else if (simulations(i)%attributesReal%exists('massParticle')) then
          idParticle   => simulations(i)%propertiesInteger%value('idParticle'  )
          massTotal=simulations(i)%attributesReal%value('massParticle')*dble(size(idParticle))
          nullify(idParticle)
       else
          massTotal=0.0d0
          call Error_Report('particle masses are not known'//{introspection:location})
       end if
       call simulations(i)%attributesReal%set           (keyCH        ='massTotal',value         =massTotal)
       call simulations(i)%analysis      %writeAttribute(attributeName='massTotal',attributeValue=massTotal)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine convexHullMassTotalOperate
