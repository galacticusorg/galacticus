!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo mass accretion histories using the \cite{wechsler_concentrations_2002} algorithm.

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryWechsler2002">
  !#  <description>Dark matter halo mass accretion histories using the \cite{wechsler_concentrations_2002} algorithm.</description>
  !# </darkMatterHaloMassAccretionHistory>

  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryWechsler2002
     !% A dark matter halo mass accretion historiy class using the \cite{wechsler_concentrations_2002} algorithm.
     private
     logical          :: formationRedshiftCompute
     double precision :: formationRedshift
   contains
     procedure :: time => wechsler2002Time
  end type darkMatterHaloMassAccretionHistoryWechsler2002
  
  interface darkMatterHaloMassAccretionHistoryWechsler2002
     !% Constructors for the {\normalfont \ttfamily wechsler2002} dark matter halo mass accretion history class.
     module procedure wechsler2002DefaultConstructor
     module procedure wechsler2002Constructor
  end interface darkMatterHaloMassAccretionHistoryWechsler2002

  ! Default parameters.
  logical          :: accretionHistoryWechslerFormationRedshiftCompute
  double precision :: accretionHistoryWechslerFormationRedshift
  
  ! Initialization status.
  logical          :: wechsler2002Initialized                          =.false.
  
contains

  function wechsler2002DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily wechsler2002} dark matter halo mass accretion history class.
    use Input_Parameters
    implicit none
    type(darkMatterHaloMassAccretionHistoryWechsler2002), target  :: wechsler2002DefaultConstructor
    
    if (.not.wechsler2002Initialized) then
       !$omp critical(wechsler2002Initialize)
       if (.not.wechsler2002Initialized) then       
          !@ <inputParameter>
          !@   <name>accretionHistoryWechslerFormationRedshiftCompute</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Compute formation redshift automatically for \cite{wechsler_concentrations_2002} halo mass accretion histories?
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('accretionHistoryWechslerFormationRedshiftCompute',accretionHistoryWechslerFormationRedshiftCompute,defaultValue=.true.)
          if (.not.accretionHistoryWechslerFormationRedshiftCompute) then
             ! In this case, read the formation redshift.
             !@ <inputParameter>
             !@   <name>accretionHistoryWechslerFormationRedshift</name>
             !@   <defaultValue>0.4</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The formation redshift to use in \cite{wechsler_concentrations_2002} halo mass accretion histories.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('accretionHistoryWechslerFormationRedshift',accretionHistoryWechslerFormationRedshift,defaultValue=0.4d0)             
          end if
          ! Record that module is now initialized.
          wechsler2002Initialized=.true.
       end if
       !$omp end critical(wechsler2002Initialize)
    end if
    wechsler2002DefaultConstructor=wechsler2002Constructor(accretionHistoryWechslerFormationRedshiftCompute,accretionHistoryWechslerFormationRedshift)
    return
  end function wechsler2002DefaultConstructor

  function wechsler2002Constructor(formationRedshiftCompute,formationRedshift)
    !% Generic constructor for the {\normalfont \ttfamily wechsler2002} dark matter halo mass accretion history class.
    use Galacticus_Error
    implicit none
    type            (darkMatterHaloMassAccretionHistoryWechsler2002), target                  :: wechsler2002Constructor
    logical                                                         , intent(in   )           :: formationRedshiftCompute
    double precision                                                , intent(in   ), optional :: formationRedshift

    wechsler2002Constructor%formationRedshiftCompute=formationRedshiftCompute
    if (.not.formationRedshiftCompute) then
       if (present(formationRedshift)) then
          wechsler2002Constructor%formationRedshift=formationRedshift
       else
          call Galacticus_Error_Report('wechsler2002Constructor','formation redshift must be provided')
       end if
    end if
    return
  end function wechsler2002Constructor
    
  double precision function wechsler2002Time(self,node,mass)
    !% Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily thisNode} using the algorithm of
    !% \cite{wechsler_concentrations_2002}.
    use Cosmology_Functions
    implicit none
    class           (darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout)          :: self
    type            (treeNode                                      ), intent(inout), pointer :: node
    double precision                                                , intent(in   )          :: mass
    class           (nodeComponentBasic                            )               , pointer :: baseBasic
    class           (cosmologyFunctionsClass                       )               , pointer :: cosmologyFunctions_
    double precision                                                                         :: expansionFactor                   , expansionFactorBase, &
         &                                                                                      mergerTreeFormationExpansionFactor

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions      ()
    baseBasic           => node              %basic()
    select case (self%formationRedshiftCompute)
    case (.true.)
       ! Compute the expansion factor at formation.
       mergerTreeFormationExpansionFactor=expansionFactorAtFormation(baseBasic%mass())
    case (.false.)
       ! Use the specified formation redshift.
       mergerTreeFormationExpansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(self%formationRedshift)
    end select
    ! Get the expansion factor at the tree base.
    expansionFactorBase=cosmologyFunctions_%expansionFactor(baseBasic%time())
    ! Compute the expansion factor for the current node.
    expansionFactor=expansionFactorBase/(1.0d0-0.5d0*log(mass/baseBasic%mass())/mergerTreeFormationExpansionFactor)
    ! Find the time corresponding to this expansion factor.
    wechsler2002Time=cosmologyFunctions_%cosmicTime(expansionFactor)
    return

  contains
    
    double precision function expansionFactorAtFormation(haloMass)
      !% Computes the expansion factor at formation using the simple model of \cite{bullock_profiles_2001}.
      use Power_Spectra
      use Critical_Overdensity
      implicit none
      double precision, intent(in   ) :: haloMass
      double precision, parameter     :: haloMassFraction   =0.015d0                            ! Wechsler et al. (2002;  Astrophysical Journal, 568:52-70).
      double precision                :: formationTime              , haloMassCharacteristic, &
           &                             sigmaCharacteristic
     
      ! Compute the characteristic mass at formation time.
      haloMassCharacteristic=haloMassFraction*haloMass
      ! Compute the corresponding rms fluctuation in the density field (i.e. sigma(M)).
      sigmaCharacteristic=Cosmological_Mass_Root_Variance(haloMassCharacteristic)
      ! Get the time at which this equals the critical overdensity for collapse.
      formationTime=Time_of_Collapse(criticalOverdensity=sigmaCharacteristic,mass=haloMass)
      ! Get the corresponding expansion factor.
      expansionFactorAtFormation=cosmologyFunctions_%expansionFactor(formationTime)
      return
    end function expansionFactorAtFormation

  end function wechsler2002Time
