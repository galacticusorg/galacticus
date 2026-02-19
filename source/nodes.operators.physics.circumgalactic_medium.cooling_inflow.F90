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
  Implements a node operator class that implements inflow of gas from the \gls{cgm} due to cooling.
  !!}

  use :: Cooling_Rates             , only : coolingRateClass
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType

  !![
  <nodeOperator name="nodeOperatorCGMCoolingInflow">
   <description>
    A node operator class that implements inflow of gas from the \gls{cgm} due to cooling.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMCoolingInflow
     !!{
     A node operator class that implements inflow of gas from the \gls{cgm} due to cooling.
     !!}
     private
     class(coolingRateClass            ), pointer :: coolingRate_ => null()
     type (enumerationComponentTypeType)          :: component
   contains
     final     ::                          cgmCoolingInflowDestructor
     procedure :: differentialEvolution => cgmCoolingInflowDifferentialEvolution
  end type nodeOperatorCGMCoolingInflow
  
  interface nodeOperatorCGMCoolingInflow
     !!{
     Constructors for the \refClass{nodeOperatorCGMCoolingInflow} node operator class.
     !!}
     module procedure cgmCoolingInflowConstructorParameters
     module procedure cgmCoolingInflowConstructorInternal
  end interface nodeOperatorCGMCoolingInflow
  
contains
  
  function cgmCoolingInflowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMCoolingInflow} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type (nodeOperatorCGMCoolingInflow)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(coolingRateClass            ), pointer       :: coolingRate_
    type (varying_string              )                :: component

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component to which cooling gas should be directed.</description>
    </inputParameter>
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=nodeOperatorCGMCoolingInflow(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function cgmCoolingInflowConstructorParameters

  function cgmCoolingInflowConstructorInternal(component,coolingRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMCoolingInflow} node operator class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    use :: Error                     , only : Error_Report
    implicit none
    type (nodeOperatorCGMCoolingInflow)                        :: self
    type (enumerationComponentTypeType), intent(in   )         :: component
    class(coolingRateClass            ), intent(in   ), target :: coolingRate_
    !![
    <constructorAssign variables="component, *coolingRate_"/>
    !!]

    if     (                                                                                                     &
         &   component /= componentTypeDisk                                                                      &
         &  .and.                                                                                                &
         &   component /= componentTypeSpheroid                                                                  &
         & ) call Error_Report("only 'disk', and 'spheroid' components are supported"//{introspection:location})
    return
  end function cgmCoolingInflowConstructorInternal

  subroutine cgmCoolingInflowDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMCoolingInflow} node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMCoolingInflow), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine cgmCoolingInflowDestructor

  subroutine cgmCoolingInflowDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform inflow of \gls{cgm} gas due to cooling.
    !!}
    use :: Abundances_Structure      , only : abundances          , operator(*)
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo, nodeComponentDisk    , nodeComponentSpheroid
    use :: Galactic_Structure_Options, only : componentTypeDisk   , componentTypeSpheroid
    implicit none
    class           (nodeOperatorCGMCoolingInflow), intent(inout), target  :: self
    type            (treeNode                    ), intent(inout), target  :: node
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: functionInterrupt
    integer                                       , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    type            (abundances                  )                         :: rateMassAbundancesCooling
    double precision                                                       :: rateMassCooling
    !$GLC attributes unused :: propertyType

    hotHalo => node%hotHalo()
    if (hotHalo%mass() <= 0.0d0) return
    rateMassCooling          =+self%coolingRate_%rate           (node)
    rateMassAbundancesCooling=+                  rateMassCooling       &
         &                    *     hotHalo     %abundances     (    ) &
         &                    /     hotHalo     %mass           (    )
    call hotHalo    %massRate         (-rateMassCooling                                      )
    call hotHalo    %abundancesRate   (-rateMassAbundancesCooling                            )
    if (self%component == componentTypeDisk) then
       disk     => node%disk    ()
       call disk    %massGasRate      (+rateMassCooling          ,interrupt,functionInterrupt)
       call disk    %abundancesGasRate(+rateMassAbundancesCooling,interrupt,functionInterrupt)
    else if (self%component == componentTypeSpheroid) then
       spheroid => node%spheroid()
       call spheroid%massGasRate      (+rateMassCooling          ,interrupt,functionInterrupt)
       call spheroid%abundancesGasRate(+rateMassAbundancesCooling,interrupt,functionInterrupt)
    else
       call Error_Report('unexpected component - this should not happen'//{introspection:location})
    end if
    return
  end subroutine cgmCoolingInflowDifferentialEvolution
  
