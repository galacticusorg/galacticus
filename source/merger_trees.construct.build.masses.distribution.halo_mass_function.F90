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
  Implementation of a merger tree halo mass function sampling class in which the sampling rate is proportional to the halo mass function.
  !!}

  use :: Cosmological_Density_Field, only : haloEnvironmentClass
  use :: Halo_Mass_Functions       , only : haloMassFunctionClass

  !![
  <mergerTreeBuildMassDistribution name="mergerTreeBuildMassDistributionHaloMassFunction">
   <description>
    A merger tree build mass distribution class in which the sampling density function equal to the halo mass function,
    \begin{equation}
    \gamma(M) = \hbox{minmax}(\phi_\mathrm{min},\phi_\mathrm{max},\mathrm{d} n(M)/\mathrm{d}\log M [1 + p_1
    \log_{10}(M/10^{13}M_\odot) + p_2 \log_{10}^2(M/10^{13}M_\odot)),
    \end{equation}
     where $\phi_\mathrm{min}=${\normalfont \ttfamily [abundanceMinimum]}, $\phi_\mathrm{max}=${\normalfont \ttfamily
     [abundanceMaximum]}, $p_1=${\normalfont \ttfamily [modifier1]}, $p_2=${\normalfont \ttfamily [modifier2]}, and
    \begin{equation}
     \hbox{minmax}(a,b,x) = \left\{\begin{array}{ll}a &amp; \hbox{ if } x &lt; a \\ x &amp; \hbox{ if } a \leq x \leq b\\ b &amp; \hbox{ if }
     x &gt; b,\end{array}\right.
    \end{equation}
    resulting in a sample of halos representative of a volume of space.
   </description>
  </mergerTreeBuildMassDistribution>
  !!]
  type, extends(mergerTreeBuildMassDistributionClass) :: mergerTreeBuildMassDistributionHaloMassFunction
     !!{
     Implementation of merger tree halo mass function sampling class in which the sampling rate is proportional to the halo mass function.
     !!}
     private
     class           (haloMassFunctionClass), pointer :: haloMassFunction_ => null()
     class           (haloEnvironmentClass ), pointer :: haloEnvironment_  => null()
     double precision                                 :: abundanceMinimum           , abundanceMaximum, &
          &                                              modifier1                  , modifier2
   contains
     final     ::           haloMassFunctionDestructor
     procedure :: sample => haloMassFunctionSample
  end type mergerTreeBuildMassDistributionHaloMassFunction

  interface mergerTreeBuildMassDistributionHaloMassFunction
     !!{
     Constructors for the \refClass{mergerTreeBuildMassDistributionHaloMassFunction} merger tree halo mass function sampling class.
     !!}
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface mergerTreeBuildMassDistributionHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the haloMassFunction merger tree halo mass function sampling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildMassDistributionHaloMassFunction)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    class           (haloMassFunctionClass                          ), pointer       :: haloMassFunction_
    class           (haloEnvironmentClass                           ), pointer       :: haloEnvironment_
    double precision                                                                 :: abundanceMinimum , abundanceMaximum, &
         &                                                                              modifier1        , modifier2

    !![
    <inputParameter>
      <name>abundanceMinimum</name>
      <defaultValue>-1.0d0</defaultValue>
      <description>The abundance (in units of Mpc$^{-3}$) below which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>abundanceMaximum</name>
      <defaultValue>-1.0d0</defaultValue>
      <description>The abundance (in units of Mpc$^{-3}$) above which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>modifier1</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Coefficient of the polynomial modifier applied to the halo mass function when sampling halo masses for tree construction.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>modifier2</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Coefficient of the polynomial modifier applied to the halo mass function when sampling halo masses for tree construction.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="parameters"/>
    <objectBuilder class="haloEnvironment" name="haloEnvironment_"   source="parameters"/>
    !!]
    self=mergerTreeBuildMassDistributionHaloMassFunction(abundanceMinimum,abundanceMaximum,modifier1,modifier2,haloMassFunction_,haloEnvironment_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloMassFunction_"/>
    <objectDestructor name="haloEnvironment_" />
    !!]
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(abundanceMinimum,abundanceMaximum,modifier1,modifier2,haloMassFunction_,haloEnvironment_) result(self)
    !!{
    Internal constructor for the haloMassFunction merger tree halo mass function sampling class.
    !!}
    implicit none
    type            (mergerTreeBuildMassDistributionHaloMassFunction)                        :: self
    double precision                                                 , intent(in   )         :: abundanceMinimum, abundanceMaximum, &
         &                                                                                      modifier1       , modifier2
    class           (haloMassFunctionClass                          ), intent(in   ), target :: haloMassFunction_
    class           (haloEnvironmentClass                           ), intent(in   ), target :: haloEnvironment_
    !![
    <constructorAssign variables="abundanceMinimum, abundanceMaximum, modifier1, modifier2, *haloMassFunction_, *haloEnvironment_"/>
    !!]

    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassDistributionHaloMassFunction} merger tree halo mass sampling class.
    !!}
    implicit none
    type(mergerTreeBuildMassDistributionHaloMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%haloMassFunction_"/>
    <objectDestructor name="self%haloEnvironment_" />
    !!]
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionSample(self,mass,time,massMinimum,massMaximum)
    !!{
    Computes the halo mass function sampling rate using a volume-limited sampling.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, nodeComponentBasic, treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class           (mergerTreeBuildMassDistributionHaloMassFunction), intent(inout) :: self
    double precision                                                 , intent(in   ) :: mass                , massMaximum, &
         &                                                                              massMinimum         , time
    double precision                                                 , parameter     :: massZeroPoint=1.0d13
    class           (nodeComponentBasic                             ), pointer       :: basic
    type            (mergerTree                                     ), target        :: tree
    !$GLC attributes unused :: massMinimum, massMaximum

    ! Create a work node.
    tree %nodeBase          => treeNode               (-1_kind_int8,tree)
    tree %nodeBase%hostTree => tree
    basic                   => tree    %nodeBase%basic(autoCreate=.true.)
    call tree %properties         %initialize          (                               )
    call self %haloEnvironment_   %overdensityLinearSet(tree%nodeBase,overdensity=0.0d0)
    call basic%timeSet                                 (time                           )
    call basic%timeLastIsolatedSet                     (time                           )
    call basic%massSet                                 (mass                           )
    ! Construct sampling rate.
    haloMassFunctionSample=+                                         mass                &
         &                 *self%haloMassFunction_%differential(time,mass,tree%nodeBase) &
         &                 *10.0d0**(                                                    &
         &                           +self%modifier1*log10(mass/massZeroPoint)           &
         &                           +self%modifier2*log10(mass/massZeroPoint)**2        &
         &                          )
    call tree%nodeBase%destroy()
    deallocate(tree%nodeBase)
    ! Limit sampling rate.
    if (self%abundanceMinimum > 0.0d0)  &
         & haloMassFunctionSample       &
         & =max(                        &
         &      haloMassFunctionSample, &
         &      self%abundanceMinimum   &
         &     )
    if (self%abundanceMaximum > 0.0d0)  &
         & haloMassFunctionSample       &
         & =min(                        &
         &      haloMassFunctionSample, &
         &      self%abundanceMaximum   &
         &     )
    return
  end function haloMassFunctionSample
