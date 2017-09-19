!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!+ Contributions to this file made by: Daniel McAndrew.

  !% An implementation of accretion from the \gls{igm} onto halos using filtering mass of the \gls{igm}
  !% calculated from an equation from \cite{naoz_formation_2007}.

  !# <accretionHalo name="accretionHaloNaozBarkana2007">
  !#  <description>Accretion onto halos using filtering mass of the \gls{igm} calculated from an equation from \cite{naoz_formation_2007}.</description>
  !# </accretionHalo>

  type, extends(accretionHaloSimple) :: accretionHaloNaozBarkana2007
     !% A halo accretion class using filtering mass of the \gls{igm} calculated from an equation from \cite{naoz_formation_2007}.
     private
     double precision                 :: rateAdjust                  , massMinimum             , &
          &                              filteredFractionRateStored  , filteredFractionStored
     logical                          :: filteredFractionRateComputed, filteredFractionComputed
     integer         (kind=kind_int8) :: lastUniqueID
   contains
     !@ <objectMethods>
     !@   <object>accretionHaloNaozBarkana2007</object>
     !@   <objectMethod>
     !@     <method>filteredFraction</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout</arguments>
     !@     <description>Returns the fraction of potential accretion onto a halo from the \gls{igm} which succeeded.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>filteredFractionRate</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout</arguments>
     !@     <description>Returns the fraction of potential accretion rate onto a halo from the \gls{igm} which succeeds.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>filteredFractionCompute</method>
     !@     <type>double precision</type>
     !@     <arguments>\doublezero\ massHalo\argin, \doublezero\ massFiltering\argin</arguments>
     !@     <description>Returns the fraction of potential accretion onto a halo from the \gls{igm} which succeeded given the halo and filtering masses.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: calculationReset       => naozBarkana2007CalculationReset
     procedure :: branchHasBaryons       => naozBarkana2007BranchHasBaryons
     procedure :: accretionRate          => naozBarkana2007AccretionRate
     procedure :: accretedMass           => naozBarkana2007AccretedMass
     procedure :: failedAccretionRate    => naozBarkana2007FailedAccretionRate
     procedure :: failedAccretedMass     => naozBarkana2007FailedAccretedMass
     procedure :: filteredFraction       => naozBarkana2007FilteredFraction
     procedure :: filteredFractionRate   => naozBarkana2007FilteredFractionRate
     procedure :: filteredFractionCompute => naozBarkana2007FilteredFractionCompute
  end type accretionHaloNaozBarkana2007

  interface accretionHaloNaozBarkana2007
     !% Constructors for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
     module procedure naozBarkana2007Constructor
     module procedure naozBarkana2007DefaultConstructor
  end interface accretionHaloNaozBarkana2007

  interface assignment(=)
     module procedure naozBarkana2007FromSimple
  end interface assignment(=)

  ! Initialization state.
  logical                     :: naozBarkana2007DefaultInitialized   =.false.

  ! Parameters controlling the Naoz & Barkana (2007) accretion rate model.
  double precision            :: naozBarkana2007RateAdjust                   , naozBarkana2007MassMinimum

  ! Virial density contrast definition used by Gnedin (2000) to define halos, and therefore used in the filtering mass fitting functions.
  double precision, parameter :: naozBarkana2007VirialDensityContrast=200.0d0
  
contains

  subroutine naozBarkana2007FromSimple(naozBarkana2007,simple)
    !% Assign a {\normalfont \ttfamily simple} halo accretion object to a {\normalfont \ttfamily naozBarkana2007} halo accretion object.
    implicit none
    type(accretionHaloNaozBarkana2007), intent(inout) :: naozBarkana2007
    type(accretionHaloSimple         ), intent(in   ) :: simple
    
    naozBarkana2007%reionizationSuppressionTime    =simple%reionizationSuppressionTime
    naozBarkana2007%reionizationSuppressionVelocity=simple%reionizationSuppressionVelocity
    naozBarkana2007%negativeAccretionAllowed       =simple%negativeAccretionAllowed
    naozBarkana2007%accreteNewGrowthOnly           =simple%accreteNewGrowthOnly
    naozBarkana2007%radiation                      =simple%radiation
    return
  end subroutine naozBarkana2007FromSimple
  
  function naozBarkana2007DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
    use Input_Parameters
    implicit none
    type(accretionHaloNaozBarkana2007), target :: naozBarkana2007DefaultConstructor

    ! Get default parameters.
    if (.not.naozBarkana2007DefaultInitialized) then
       !$omp critical(accretionHaloNaozBarkana2007DefaultInitialize)
       if (.not.naozBarkana2007DefaultInitialized) then
          !# <inputParameter>
          !#   <name>accretionHaloNaozBarkana2007RateAdjust</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.3d0</defaultValue>
          !#   <description>The dimensionless multiplier for the rate at which the halo gas content adjusts to changes in the filtering mass.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>naozBarkana2007RateAdjust</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>accretionHaloNaozBarkana2007MassMinimum</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.0d0</defaultValue>
          !#   <description>The minimum mass of gas accreted into a halo below which the mass is truncated to zero.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>naozBarkana2007MassMinimum</variable>
          !# </inputParameter>
          ! Record that class is now initialized.
          naozBarkana2007DefaultInitialized=.true.
       end if
       !$omp end critical(accretionHaloNaozBarkana2007DefaultInitialize)
    end if
    naozBarkana2007DefaultConstructor            =accretionHaloSimple()
    naozBarkana2007DefaultConstructor%rateAdjust =naozBarkana2007RateAdjust
    naozBarkana2007DefaultConstructor%massMinimum=naozBarkana2007MassMinimum
    return
  end function naozBarkana2007DefaultConstructor
       
  function naozBarkana2007Constructor(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly,rateAdjust,massMinimum)
    !% Default constructor for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
    implicit none
    type            (accretionHaloNaozBarkana2007), target        :: naozBarkana2007Constructor
    double precision                              , intent(in   ) :: reionizationSuppressionTime, reionizationSuppressionVelocity, &
         &                                                           rateAdjust                 , massMinimum
    logical                                                       :: negativeAccretionAllowed   , accreteNewGrowthOnly
    
    naozBarkana2007Constructor            =accretionHaloSimple(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    naozBarkana2007Constructor%rateAdjust =rateAdjust
    naozBarkana2007Constructor%massMinimum=massMinimum
    return
  end function naozBarkana2007Constructor

  subroutine naozBarkana2007CalculationReset(self,node)
    !% Reset the accretion rate calculation.
    implicit none
    class(accretionHaloNaozBarkana2007), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    self%filteredFractionComputed    =.false.
    self%filteredFractionRateComputed=.false.
    self%lastUniqueID                =node%uniqueID()
    return
  end subroutine naozBarkana2007CalculationReset

  logical function naozBarkana2007BranchHasBaryons(self,node)
    !% Returns true if this branch can accrete any baryons.
    use Galacticus_Nodes
    use Accretion_Halo_Totals
    use Cosmology_Parameters
    implicit none
    class           (accretionHaloNaozBarkana2007), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    type            (treeNode                    )               , pointer :: branchNode
    class           (cosmologyParametersClass    )               , pointer :: cosmologyParameters_
    class           (accretionHaloTotalClass     )               , pointer :: accretionHaloTotal_
    class           (nodeComponentBasic          )               , pointer :: basic
    double precision                                                       :: fractionBaryons     , massHaloMinimum

    cosmologyParameters_            =>  cosmologyParameters             ()
    accretionHaloTotal_             =>  accretionHaloTotal              ()
    fractionBaryons                 =  +cosmologyParameters_%OmegaBaryon() &
         &                             /cosmologyParameters_%OmegaMatter()
    massHaloMinimum                 =  +self                %massMinimum   &
         &                             /fractionBaryons
    naozBarkana2007BranchHasBaryons =   .false.
    branchNode                      =>  node
    do while (associated(branchNode))
       basic => branchnode%basic()
       if (accretionHaloTotal_%accretedMass(branchNode)*self%filteredFraction(branchNode) >= massHaloMinimum) then
          naozBarkana2007BranchHasBaryons=.true.
          exit
       end if
       branchNode => branchNode%walkBranch(node)
    end do
    return
  end function naozBarkana2007BranchHasBaryons

  double precision function naozBarkana2007FilteredFraction(self,node)
    !% Returns the baryonic mass fraction in a halo after the effects of the filtering mass.
    use Galacticus_Nodes
    use Intergalactic_Medium_State
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    class           (accretionHaloNaozBarkana2007 ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (nodeComponentBasic           ), pointer       :: basic
    class           (intergalacticMediumStateClass), pointer       :: igmState_
    double precision                                               :: massFiltering, massHalo

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Evaluate the filtering mass suppression fitting formula as defined by Naoz & Barkana (2007;
    ! http://adsabs.harvard.edu/abs/2007MNRAS.377..667N). We use a halo mass in this formula defined in the same way (∆=200) as in
    ! the original work by Gnedin (2000; http://adsabs.harvard.edu/abs/2000ApJ...542..535G) based on the discussion of halo
    ! definition in Naoz, Yoshida, & Gnedin (2013; http://adsabs.harvard.edu/abs/2013ApJ...763...27N).
    if (.not.self%filteredFractionComputed) then
       igmState_                     => intergalacticMediumState                         (                                                  )
       basic                         => node                               %basic        (                                                  )
       massFiltering                 =  igmState_                          %filteringMass(basic%time()                                      )
       massHalo                      =  Dark_Matter_Profile_Mass_Definition              (node   ,naozBarkana2007VirialDensityContrast)
       self%filteredFractionStored   =  self%filteredFractionCompute(massHalo,massFiltering)
       self%filteredFractionComputed =  .true.
    end if
    naozBarkana2007FilteredFraction=self%filteredFractionStored
    return
  end function naozBarkana2007FilteredFraction

  double precision function naozBarkana2007FilteredFractionRate(self,node)
    !% Returns the baryonic mass accretion rate fraction in a halo after the effects of the filtering mass.
    use Galacticus_Nodes
    use Intergalactic_Medium_State
    use Dark_Matter_Profile_Mass_Definitions
    use Math_Exponentiation
    implicit none
    class           (accretionHaloNaozBarkana2007 ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (nodeComponentBasic           ), pointer       :: basic
    class           (intergalacticMediumStateClass), pointer       :: igmState_
    double precision                                               :: massFiltering, massHalo

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Evaluate the rate of change of the filtering mass suppression fitting formula as defined by Naoz & Barkana (2007;
    ! http://adsabs.harvard.edu/abs/2007MNRAS.377..667N). We use a halo mass in this formula defined in the same way (∆=200) as in
    ! the original work by Gnedin (2000; http://adsabs.harvard.edu/abs/2000ApJ...542..535G) based on the discussion of halo
    ! definition in Naoz, Yoshida, & Gnedin (2013; http://adsabs.harvard.edu/abs/2013ApJ...763...27N). The rate of change here
    ! assumes that the filtering mass is constant in time.
    if (.not.self%filteredFractionRateComputed) then
       igmState_                        => intergalacticMediumState                         (                                                  )
       basic                            => node                               %basic        (                                                  )
       massFiltering                    =  igmState_                          %filteringMass(basic%time()                                      )
       massHalo                         =  Dark_Matter_Profile_Mass_Definition              (node         ,naozBarkana2007VirialDensityContrast)
       if (.not.self%filteredFractionComputed) then
          self%filteredFractionStored   =  self%filteredFractionCompute(massHalo,massFiltering)
          self%filteredFractionComputed =  .true.
       end if
       self%filteredFractionRateStored  = +           self%filteredFractionStored  &
            &                             *(+1.0d0                                 &
            &                               +cubeRoot(self%filteredFractionStored) &
            &                               *(                                     &
            &                                 +2.0d0**(1.0d0/3.0d0)                &
            &                                 -1.0d0                               &
            &                                )                                     &
            &                               *24.0d0                                &
            &                               *massFiltering                         &
            &                               /massHalo                              &
            &                              )
       self%filteredFractionRateComputed=.true.
    end if
    naozBarkana2007FilteredFractionRate=self%filteredFractionRateStored
    return
  end function naozBarkana2007FilteredFractionRate

  double precision function naozBarkana2007FilteredFractionCompute(self,massHalo,massFiltering)
    !% Compute the filtered fraction.
    implicit none
    double precision                               , intent(in   ) :: massHalo, massFiltering
    class           (accretionHaloNaozBarkana2007 ), intent(inout) :: self
    !GCC$ attributes unused :: self

    naozBarkana2007FilteredFractionCompute=+1.0d0                     &
         &                                 /(                         &
         &                                    +1.0d0                  &
         &                                    +(                      &
         &                                      +2.0d0**(1.0d0/3.0d0) &
         &                                      -1.0d0                &
         &                                     )                      &
         &                                    *8.0d0                  &
         &                                    *massFiltering          &
         &                                    /massHalo               &
         &                                   )**3
    return
  end function naozBarkana2007FilteredFractionCompute

  double precision function naozBarkana2007AccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Accretion_Halo_Totals
    use Dark_Matter_Halo_Scales
    implicit none
    class           (accretionHaloNaozBarkana2007), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer                                       , intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (nodeComponentHotHalo        ), pointer       :: hotHalo
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class           (accretionHaloTotalClass     ), pointer       :: accretionHaloTotal_
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: growthRate          , filteredFraction, &
         &                                                           filteredFractionRate, fractionAccreted

    naozBarkana2007AccretionRate=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    accretionHaloTotal_  => accretionHaloTotal ()
    darkMatterHaloScale_ => darkMatterHaloScale()
    basic                => node%basic         ()
    hotHalo              => node%hotHalo       ()
    ! Find the post-filtering accretion rate fraction.
    filteredFractionRate=self%filteredFractionRate(node)
    ! Compute the mass accretion rate onto the halo.
    naozBarkana2007AccretionRate=+cosmologyParameters_%OmegaBaryon  (    ) &
         &                       /cosmologyParameters_%OmegaMatter  (    ) &
         &                       *accretionHaloTotal_ %accretionRate(node) &
         &                       *filteredFractionRate
    ! Test for negative accretion.
    if (.not.self%negativeAccretionAllowed.and.accretionHaloTotal_%accretionRate(node) < 0.0d0) then
       ! Accretion rate is negative, and not allowed. Return zero accretion rate.
       naozBarkana2007AccretionRate=+0.0d0
    else
       ! Adjust the rate to allow mass to flow back-and-forth from accreted to unaccreted reservoirs if the current mass fraction
       ! differs from that expected given the filtering mass.
       growthRate                  =+self                %rateAdjust               &
            &                       /darkMatterHaloScale_%dynamicalTimescale(node)
       filteredFraction            =+self                %filteredFraction  (node)
       fractionAccreted            =+  hotHalo           %          mass    (    ) &
            &                       /(                                             &
            &                         +hotHalo           %          mass    (    ) &
            &                         +hotHalo           %unaccretedMass    (    ) &
            &                        )
       naozBarkana2007AccretionRate=+naozBarkana2007AccretionRate                  &
            &                       -(                                             &
            &                         +hotHalo           %          mass    (    ) &
            &                         +hotHalo           %unaccretedMass    (    ) &
            &                        )                                             &
            &                       *(                                             &
            &                         +fractionAccreted                            &
            &                         -filteredFraction                            &
            &                        )                                             &
            &                       *growthRate
    end if
    ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
    if (self%accreteNewGrowthOnly .and. accretionHaloTotal_%accretedMass(node) < basic%massMaximum()) naozBarkana2007AccretionRate=0.0d0
    return
  end function naozBarkana2007AccretionRate

  double precision function naozBarkana2007AccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Accretion_Halo_Totals
    implicit none
    class           (accretionHaloNaozBarkana2007), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer                                       , intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class           (accretionHaloTotalClass     ), pointer       :: accretionHaloTotal_
    double precision                                              :: filteredFraction

    naozBarkana2007AccretedMass=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get required objects.
    basic                => node               %basic           (    )
    cosmologyParameters_ => cosmologyParameters                 (    )
    accretionHaloTotal_  => accretionHaloTotal                  (    )
    ! Get the filtered mass fraction.
    filteredFraction     =self                 %filteredFraction(node)
    ! Get the default cosmology.
    naozBarkana2007AccretedMass=+cosmologyParameters_%OmegaBaryon (    ) &
         &                      /cosmologyParameters_%OmegaMatter (    ) &
         &                      *accretionHaloTotal_ %accretedMass(node) &
         &                      *filteredFraction
         return
  end function naozBarkana2007AccretedMass

  double precision function naozBarkana2007FailedAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Accretion_Halo_Totals
    implicit none
    class           (accretionHaloNaozBarkana2007), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer                                       , intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (nodeComponentHotHalo        ), pointer       :: hotHalo
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class           (accretionHaloTotalClass     ), pointer       :: accretionHaloTotal_
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: growthRate          , filteredFraction, &
         &                                                           filteredFractionRate, fractionAccreted

    naozBarkana2007FailedAccretionRate=0.0d0
    if (accretionMode               == accretionModeCold) return
    if (node         %isSatellite()                     ) return
    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters                     (    )
    accretionHaloTotal_  => accretionHaloTotal                      (    )
    darkMatterHaloScale_ => darkMatterHaloScale                     (    )
    basic                => node               %basic               (    )
    hotHalo              => node               %hotHalo             (    )
    ! Get the post-filtering accretion rate fraction.
    filteredFractionRate =  self               %filteredFractionRate(node)
    ! Test for negative accretion.
    if (.not.self%negativeAccretionAllowed.and.accretionHaloTotal_%accretionRate(node) < 0.0d0) then
       naozBarkana2007FailedAccretionRate=+cosmologyParameters_%OmegaBaryon  (    ) &
            &                             /cosmologyParameters_%OmegaMatter  (    ) &
            &                             *accretionHaloTotal_ %accretionRate(node)
    else
       ! Compute the rate of failed accretion.
       naozBarkana2007FailedAccretionRate=+cosmologyParameters_%OmegaBaryon  (    ) &
            &                             /cosmologyParameters_%OmegaMatter  (    ) &
            &                             *accretionHaloTotal_ %accretionRate(node) &
            &                             *(                                        &
            &                               +1.0d0                                  &
            &                               -filteredFractionRate                   &
            &                              )
       ! Adjust the rate to allow mass to flow back-and-forth from accreted to unaccreted reservoirs if the current mass fraction
       ! differs from that expected given the filtering mass.
       growthRate                        =+self                %rateAdjust               &
            &                             /darkMatterHaloScale_%dynamicalTimescale(node)
       filteredFraction                  =+self                %filteredFraction  (node)
       fractionAccreted                  =+  hotHalo           %          mass    (    ) &
            &                             /(                                             &
            &                               +hotHalo           %          mass    (    ) &
            &                               +hotHalo           %unaccretedMass    (    ) &
            &                              )
       naozBarkana2007FailedAccretionRate=+naozBarkana2007FailedAccretionRate            &
            &                             +(                                             &
            &                               +hotHalo           %          mass    (    ) &
            &                               +hotHalo           %unaccretedMass    (    ) &
            &                              )                                             &
            &                             *(                                             &
            &                               +fractionAccreted                            &
            &                               -filteredFraction                            &
            &                              )                                             &
            &                             *growthRate
    end if
    ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
    if (self%accreteNewGrowthOnly .and. accretionHaloTotal_%accretedMass(node) < basic%massMaximum()) naozBarkana2007FailedAccretionRate=0.0d0
    return
  end function naozBarkana2007FailedAccretionRate

  double precision function naozBarkana2007FailedAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Accretion_Halo_Totals
    implicit none
    class           (accretionHaloNaozBarkana2007), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer                                       , intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class           (accretionHaloTotalClass     ), pointer       :: accretionHaloTotal_
    double precision                                              :: filteredFraction

    naozBarkana2007FailedAccretedMass=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get required objects.
    basic                             => node                %basic           (    )
    cosmologyParameters_              => cosmologyParameters                  (    )
    accretionHaloTotal_               => accretionHaloTotal                   (    )
    ! Get the failed fraction.
    filteredFraction                  =  self                %filteredFraction(node)
    ! Get the default cosmology.
    naozBarkana2007FailedAccretedMass = +cosmologyParameters_%OmegaBaryon     (    ) &
         &                              /cosmologyParameters_%OmegaMatter     (    ) &
         &                              *accretionHaloTotal_ %accretedMass    (node) &
         &                              *(                                           &
         &                                +1.0d0                                     &
         &                                -filteredFraction                          &
         &                               )
    return
  end function naozBarkana2007FailedAccretedMass
  
