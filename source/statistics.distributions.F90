!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements a class of distributions.

module Statistics_Distributions
  use Pseudo_Random
  private
  public :: distributionNew

  ! Define the basic distribution type.
  type, public :: distribution
  end type distribution

  ! Define a 1D distribution.
  type, abstract, public, extends(distribution) :: distribution1D
     type(pseudoRandom) :: randomNumberGenerator
   contains
     !@ <objectMethods>
     !@   <object>distribution1D</object>
     !@   <objectMethod>
     !@     <method>density</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero x\argin</arguments>
     !@     <description>Return the probability density at {\tt x}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>cumulative</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero x\argin</arguments>
     !@     <description>Return the cumulative probability at {\tt x}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inverse</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero p\argin</arguments>
     !@     <description>Return the value of the independent variable corresponding to cumulative probability {\tt p}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>sample</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return a random deviate from the distribution.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(distribution1DDensity   ), deferred :: density
     procedure(distribution1DCumulative), deferred :: cumulative
     procedure                                     :: inverse    => distribution1DInverse
     procedure                                     :: sample     => distribution1DSample
  end type distribution1D
  
  abstract interface
     double precision function distribution1DDensity(self,x)
       import :: distribution1D
       class           (distribution1D), intent(in   ) :: self
       double precision                , intent(in   ) :: x
     end function distribution1DDensity
  end interface
  
  abstract interface
     double precision function distribution1DCumulative(self,x)
       import :: distribution1D
       class           (distribution1D), intent(in   ) :: self
       double precision                , intent(in   ) :: x
     end function distribution1DCumulative
  end interface
  
  ! Define a list of distributions.
  type, public :: distributionList
     class(distribution), pointer :: thisDistribution
  end type distributionList

  ! Include all distribution types.
  include 'statistics.distributions.uniform.type.inc'
  include 'statistics.distributions.loguniform.type.inc'
  include 'statistics.distributions.normal.type.inc'
  include 'statistics.distributions.Cauchy.type.inc'
  include 'statistics.distributions.Student-t.type.inc'
  include 'statistics.distributions.Gamma.type.inc'
  include 'statistics.distributions.Voight.type.inc'

contains

  function distributionNew(definition) result (newDistribution)
    !% Create a new distribution from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class           (distribution), pointer                :: newDistribution
    type            (node        ), pointer, intent(in   ) :: definition
    type            (node        ), pointer                :: typeElement
    double precision                                       :: distributionMinimum         , distributionMaximum      , &
         &                                                    distributionMean            , distributionVariance     , &
         &                                                    distributionScale           , distributionMedian       , &
         &                                                    distributionDegreesOfFreedom, distributionShape        , &
         &                                                    distributionRate            , distributionWidth
    logical                                                :: distributionMinimumExists   , distributionMaximumExists

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("uniform")
       allocate(distributionUniform :: newDistribution)
       select type (newDistribution)
       type is (distributionUniform)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"minimum"),distributionMinimum)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"maximum"),distributionMaximum)
          newDistribution=distributionUniform(distributionMinimum,distributionMaximum)
       end select
    case ("logUniform")
       allocate(distributionLogUniform :: newDistribution)
       select type (newDistribution)
       type is (distributionLogUniform)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"minimum"),distributionMinimum)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"maximum"),distributionMaximum)
          newDistribution=distributionLogUniform(distributionMinimum,distributionMaximum)
       end select
    case ("normal")
       allocate(distributionNormal :: newDistribution)
       select type (newDistribution)
       type is (distributionNormal)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"mean"    ),distributionMean    )
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"variance"),distributionVariance)
          distributionMinimumExists=XML_Path_Exists(definition,"minimum")
          distributionMaximumExists=XML_Path_Exists(definition,"maximum")
          if (distributionMinimumExists) call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"minimum"),distributionMinimum)
          if (distributionMaximumExists) call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"maximum"),distributionMaximum)
          if      (     distributionMinimumExists.and.     distributionMaximumExists) then
             newDistribution=distributionNormal(&
                  &                             distributionMean                        , &
                  &                             distributionVariance                    , &
                  &                             limitLower          =distributionMinimum, &
                  &                             limitUpper          =distributionMaximum  &
                  &                            )
          else if (     distributionMinimumExists.and..not.distributionMaximumExists) then
             newDistribution=distributionNormal(                                          &
                  &                             distributionMean                        , &
                  &                             distributionVariance                    , &
                  &                             limitLower          =distributionMinimum  &
                  &                            )
          else if (.not.distributionMinimumExists.and.     distributionMaximumExists) then
             newDistribution=distributionNormal(                                          &
                  &                             distributionMean                        , &
                  &                             distributionVariance                    , &
                  &                             limitUpper          =distributionMaximum  &
                  &                            )
          else if (.not.distributionMinimumExists.and..not.distributionMaximumExists) then
             newDistribution=distributionNormal(                                          &
                  &                             distributionMean                        , &
                  &                             distributionVariance                      &
                  &                            )
          end if
       end select
    case ("Gamma")
       allocate(distributionGamma :: newDistribution)
       select type (newDistribution)
       type is (distributionGamma)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"rate" ),distributionRate )
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"shape"),distributionShape)
          distributionMinimumExists=XML_Path_Exists(definition,"minimum")
          distributionMaximumExists=XML_Path_Exists(definition,"maximum")
          if (distributionMinimumExists) call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"minimum"),distributionMinimum)
          if (distributionMaximumExists) call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"maximum"),distributionMaximum)
          if      (     distributionMinimumExists.and.     distributionMaximumExists) then
             newDistribution=distributionGamma(                                           &
                  &                             distributionShape                       , &
                  &                             distributionRate                        , &
                  &                             limitLower          =distributionMinimum, &
                  &                             limitUpper          =distributionMaximum  &
                  &                            )
          else if (     distributionMinimumExists.and..not.distributionMaximumExists) then
             newDistribution=distributionGamma(                                           &
                  &                             distributionShape                       , &
                  &                             distributionRate                        , &
                  &                             limitLower          =distributionMinimum  &
                  &                            )
          else if (.not.distributionMinimumExists.and.     distributionMaximumExists) then
             newDistribution=distributionGamma(                                           &
                  &                             distributionShape                       , &
                  &                             distributionRate                        , &
                  &                             limitUpper          =distributionMaximum  &
                  &                            )
          else if (.not.distributionMinimumExists.and..not.distributionMaximumExists) then
             newDistribution=distributionGamma(                                           &
                  &                             distributionShape                       , &
                  &                             distributionRate                          &
                  &                            )
          end if
       end select
    case ("Cauchy")
       allocate(distributionCauchy :: newDistribution)
       select type (newDistribution)
       type is (distributionCauchy)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"median"),distributionMedian)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"scale" ),distributionScale )
          newDistribution=distributionCauchy(distributionMedian,distributionScale)
       end select
    case ("StudentT")
       allocate(distributionStudentT :: newDistribution)
       select type (newDistribution)
       type is (distributionStudentT)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"degreesOfFreedom"),distributionDegreesOfFreedom)
          newDistribution=distributionStudentT(distributionDegreesOfFreedom)
       end select
    case ("Voight")
       allocate(distributionVoight :: newDistribution)
       select type (newDistribution)
       type is (distributionVoight)
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"scale"   ),distributionScale   )
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"median"  ),distributionMedian  )
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(definition,"variance"),distributionVariance)
          newDistribution=distributionVoight(distributionScale,distributionMedian,sqrt(distributionVariance))
       end select
    case default
       call Galacticus_Error_Report('distributionNew','distribution type is unrecognized')
    end select
    return
  end function distributionNew

  double precision function distribution1DSample(self)
    !% Sample from a 1D distribution.
    implicit none
    class(distribution1D), intent(inout) :: self
 
    ! Draw a random number uniformly from 0 to 1 and use the inverse of our self to get the
    ! corresponding random variate.
    distribution1DSample=self%inverse(                                                          &
         &                            self%randomNumberGenerator%sample(                        &
         &                                                              ompThreadOffset=.true., &
         &                                                              mpiRankOffset  =.true.  &
         &                                                             )                        &
         &                           )
    return
  end function distribution1DSample

  double precision function distribution1DInverse(self,p)
    !% Null implementation of 1D inversion.
    use Galacticus_Error
    class           (distribution1D), intent(in   ) :: self
    double precision                , intent(in   ) :: p

    call Galacticus_Error_Report('distribution1DInverse', &
         &                       'not implemented'        &
         &                      )
    return
  end function distribution1DInverse

  ! Include all distribution methods.
  include 'statistics.distributions.uniform.methods.inc'
  include 'statistics.distributions.loguniform.methods.inc'
  include 'statistics.distributions.normal.methods.inc'
  include 'statistics.distributions.Cauchy.methods.inc'
  include 'statistics.distributions.Student-t.methods.inc'
  include 'statistics.distributions.Gamma.methods.inc'
  include 'statistics.distributions.Voight.methods.inc'

end module Statistics_Distributions
