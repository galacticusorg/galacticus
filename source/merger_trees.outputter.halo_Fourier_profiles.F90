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
  Implements a merger tree outputter class that outputs $k$-space density profiles as needed for halo model calculations.
  !!}
  
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Galactic_Filters        , only : galacticFilterClass
  use :: IO_HDF5                 , only : hdf5Object

  !![
  <mergerTreeOutputter name="mergerTreeOutputterHaloFourierProfiles">
   <description>
    A merger tree outputter class which outputs $k$-space density profiles as needed for halo model calculations. A
    ``{\normalfont \ttfamily haloModel}'' group is created in the \glc\ output file. This group contains the following:
  
    \begin{description}
  
     \item [{\normalfont \ttfamily wavenumber}] A dataset giving the wavenumbers (in units of Mpc$^{-1}$) at which all output
     power spectra are tabulated. The minimum and maximum wavenumbers to tabulate are determined by the {\normalfont \ttfamily
     [haloModelWavenumberMinimum]} and {\normalfont \ttfamily [haloModelWavenumberMaximum]} parameters respectively, while the
     number of points to tabulate in each decade of wavenumber is determined by the {\normalfont \ttfamily
     [haloModelWavenumberPointsPerDecade]} parameter.
  
     \item [{\normalfont \ttfamily powerSpectrum}] A dataset giving the linear theory power spectrum (in units of Mpc$^3$
     normalized to $z=0$ at each wavenumber specified in the {\normalfont \ttfamily wavenumber} dataset.
  
     \item [{\normalfont \ttfamily Output\{i\}/mergerTree\{j\}/fourierProfile\{k\}}] A dataset giving the Fourier transform of
     the dark matter halo density profile (dimensionless and normalized to unity at small wavenumber) for the node with index
     {\normalfont \ttfamily k} in merger tree with index {\normalfont \ttfamily j} at output number {\normalfont \ttfamily
     i}. Profiles are written only for nodes which are isolated, and are tabulated at the wavenumbers given in the {\normalfont
     \ttfamily wavenumber} group. Note that wavenumbers are assumed to be comoving.
  
    \end{description}
  
    Finally, each numbered output group is given two additional attributes, {\normalfont \ttfamily linearGrowthFactor} and
    {\normalfont \ttfamily linearGrowthFactorLogDerivative} which give the growth factor, $D$, and its logarithmic derivative,
    $\d \ln D / \d \ln a$ at the output time.
   </description>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterHaloFourierProfiles
     !!{
     Implementation of a merger tree outputter class that outputs $k$-space density profiles as needed for halo model calculations.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer                   :: cosmologyFunctions_       => null()
     class           (darkMatterHaloScaleClass ), pointer                   :: darkMatterHaloScale_      => null()
     class           (darkMatterProfileDMOClass), pointer                   :: darkMatterProfileDMO_     => null()
     class           (galacticFilterClass      ), pointer                   :: galacticFilter_           => null()
     integer                                                                :: wavenumberPointsPerDecade          , wavenumberCount
     double precision                                                       :: wavenumberMaximum                  , wavenumberMinimum
     double precision                           , allocatable, dimension(:) :: wavenumber
     type            (hdf5Object               )                            :: outputGroup
   contains
     final     ::               haloFourierProfilesDestructor
     procedure :: outputTree => haloFourierProfilesOutputTree
     procedure :: outputNode => haloFourierProfilesOutputNode
     procedure :: finalize   => haloFourierProfilesFinalize
  end type mergerTreeOutputterHaloFourierProfiles

  interface mergerTreeOutputterHaloFourierProfiles
     !!{
     Constructors for the \refClass{mergerTreeOutputterHaloFourierProfiles} merger tree outputter.
     !!}
     module procedure haloFourierProfilesConstructorParameters
     module procedure haloFourierProfilesConstructorInternal
  end interface mergerTreeOutputterHaloFourierProfiles
  
contains
  
  function haloFourierProfilesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterHaloFourierProfiles} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeOutputterHaloFourierProfiles)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    class           (galacticFilterClass                   ), pointer       :: galacticFilter_
    double precision                                                        :: wavenumberMinimum        , wavenumberMaximum
    integer                                                                 :: wavenumberPointsPerDecade

    !![
    <inputParameter>
      <name>wavenumberPointsPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of points per decade in wavenumber at which to tabulate power spectra for the halo model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The minimum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMaximum</name>
      <defaultValue>1.0d4</defaultValue>
      <description>The maximum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticFilter"       name="galacticFilter_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=mergerTreeOutputterHaloFourierProfiles(wavenumberPointsPerDecade,wavenumberMinimum,wavenumberMaximum,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"      />
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function haloFourierProfilesConstructorParameters

  function haloFourierProfilesConstructorInternal(wavenumberPointsPerDecade,wavenumberMinimum,wavenumberMaximum,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOutputterHaloFourierProfiles} merger tree outputter class.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (mergerTreeOutputterHaloFourierProfiles)                        :: self
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass             ), intent(in   ), target :: darkMatterProfileDMO_
    class           (galacticFilterClass                   ), intent(in   ), target :: galacticFilter_
    double precision                                        , intent(in   )         :: wavenumberMinimum        , wavenumberMaximum
    integer                                                 , intent(in   )         :: wavenumberPointsPerDecade
    !![
    <constructorAssign variables="wavenumberPointsPerDecade, wavenumberMinimum, wavenumberMaximum, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *galacticFilter_"/>
    !!]
    
    ! Build a grid of wavenumbers.
    self%wavenumberCount=int(log10(self%wavenumberMaximum/self%wavenumberMinimum)*dble(self%wavenumberPointsPerDecade))+1
    allocate(self%wavenumber(self%wavenumberCount))
    self%wavenumber=Make_Range(self%wavenumberMinimum,self%wavenumberMaximum,self%wavenumberCount,rangeType=rangeTypeLogarithmic)
    return
  end function haloFourierProfilesConstructorInternal

  subroutine haloFourierProfilesDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOutputterHaloFourierProfiles} merger tree outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterHaloFourierProfiles), intent(inout) :: self
    
    call self%finalize()
    !![
    <objectDestructor name="self%galacticFilter_"      />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine haloFourierProfilesDestructor

  subroutine haloFourierProfilesFinalize(self)
    !!{
    Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    !!}
    !$ use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeOutputterHaloFourierProfiles), intent(inout) :: self

    !$ call hdf5Access%set  ()
    if (self%outputGroup%isOpen()) call self%outputGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine haloFourierProfilesFinalize
  
  subroutine haloFourierProfilesOutputTree(self,tree,indexOutput,time)
    !!{
    Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    !!}
    use    :: Output_HDF5                     , only : outputFile
    use    :: Galacticus_Nodes                , only : treeNode                , nodeComponentBasic
    !$ use :: HDF5_Access                     , only : hdf5Access
    use    :: ISO_Varying_String              , only : var_str
    use    :: Mass_Distributions              , only : massDistributionClass
    use    :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    use    :: Numerical_Constants_Astronomical, only : megaParsec
    use    :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerTreeOutputterHaloFourierProfiles), intent(inout)               :: self
    type            (mergerTree                            ), intent(inout), target       :: tree
    integer         (c_size_t                              ), intent(in   )               :: indexOutput
    double precision                                        , intent(in   )               :: time
    type            (treeNode                              )               , pointer      :: node
    class           (nodeComponentBasic                    )               , pointer      :: basic
    class           (massDistributionClass                 )               , pointer      :: massDistribution_
    double precision                                        , allocatable  , dimension(:) :: fourierProfile
    type            (mergerTreeWalkerAllNodes              )                              :: treeWalker
    type            (hdf5Object                            )                              :: outputGroup      , treeGroup   , &
         &                                                                                   dataset
    integer         (c_size_t                              )                              :: treeIndexPrevious
    double precision                                                                      :: expansionFactor  , radiusVirial
    integer                                                                               :: i
    !$GLC attributes unused :: time
    
    allocate(fourierProfile(self%wavenumberCount))
    !$ call hdf5Access%set  ()
    if (.not.self%outputGroup%isOpen()) then
       self%outputGroup=outputFile%openGroup("haloFourierProfiles","Halo model data.")
       call self   %outputGroup%writeDataset  (self%wavenumber ,'wavenumber','Wavenumber at which Fourier transform of density profile is tabulated [Mpc⁻¹].',datasetReturned=dataset)
       call dataset            %writeAttribute(1.0d0/megaParsec,'unitsInSI'                                                                                                          )
       call dataset            %close         (                                                                                                                                      )
    end if
    outputGroup=self%outputGroup%openGroup(char(var_str('output')//indexOutput),char(var_str("Fourier space density profiles of halos for all trees at output number ")//indexOutput//"."))
    !$ call hdf5Access%unset()
    treeIndexPrevious=-huge(0_c_size_t)
    treeWalker       =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       if (.not.self%galacticFilter_%passes(node)) cycle       
       if (node%hostTree%index /= treeIndexPrevious) then
          treeIndexPrevious=node%hostTree%index
          !$ call hdf5Access%set  ()
          if (treeGroup%isOpen()) call treeGroup%close()
          treeGroup=outputGroup%openGroup(char(var_str('tree')//node%hostTree%index),"Fourier space density profiles of halos for each tree.")
          !$ call hdf5Access%unset()
       end if
       basic           => node%basic                              (            )
       expansionFactor =  self%cosmologyFunctions_%expansionFactor(basic%time())
       ! Construct profile. (Our wavenumbers are comoving, so we must convert them to physical coordinates before passing them to
       ! the dark matter profile k-space routine.)
       massDistribution_ => self%darkMatterProfileDMO_%get         (node)
       radiusVirial      =  self%darkMatterHaloScale_ %radiusVirial(node)
       do i=1,self%waveNumberCount
          fourierProfile(i)=massDistribution_%fourierTransform(radiusVirial,self%wavenumber(i)/expansionFactor)
       end do
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       !$ call hdf5Access%set  ()
       call treeGroup%writeDataset(fourierProfile,char(var_str('node')//node%index()),"The Fourier-space density profile.")
       !$ call hdf5Access%unset()
    end do
    !$ call hdf5Access%set  ()
    if (treeGroup%isOpen()) call treeGroup  %close()
    call                         outputGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine haloFourierProfilesOutputTree

  subroutine haloFourierProfilesOutputNode(self,node,indexOutput)
    !!{
    Perform no output.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (mergerTreeOutputterHaloFourierProfiles), intent(inout) :: self
    type   (treeNode                              ), intent(inout) :: node
    integer(c_size_t                              ), intent(in   ) :: indexOutput
    !$GLC attributes unused :: self, node, indexOutput

    call Error_Report('output of single nodes is not supported'//{introspection:location})
    return
  end subroutine haloFourierProfilesOutputNode
