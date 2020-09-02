!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.

  use :: Kind_Numbers, only : kind_int8

  !# <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksBlitz2006">
  !#  <description>The \cite{blitz_role_2006} star formation rate surface density law for galactic disks.</description>
  !# </starFormationRateSurfaceDensityDisks>
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksBlitz2006
     !% Implementation of the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.
     private
     integer         (kind_int8) :: lastUniqueID
     logical                     :: factorsComputed
     double precision            :: heightToRadialScaleDisk  , pressureCharacteristic             , &
          &                         pressureExponent         , starFormationFrequencyNormalization, &
          &                         surfaceDensityCritical   , surfaceDensityExponent             , &
          &                         velocityDispersionDiskGas, radiusDisk                         , &
          &                         massGas                  , hydrogenMassFraction               , &
          &                         massStellar
   contains
     !@ <objectMethods>
     !@   <object>starFormationRateSurfaceDensityDisksBlitz2006</object>
     !@   <objectMethod>
     !@     <method>calculationReset</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} node\arginout</arguments>
     !@     <description>Reset memoized calculations.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                     blitz2006Destructor
     procedure :: autoHook         => blitz2006AutoHook
     procedure :: calculationReset => blitz2006CalculationReset
     procedure :: rate             => blitz2006Rate
  end type starFormationRateSurfaceDensityDisksBlitz2006

  interface starFormationRateSurfaceDensityDisksBlitz2006
     !% Constructors for the {\normalfont \ttfamily blitz2006} star formation surface density rate in disks class.
     module procedure blitz2006ConstructorParameters
     module procedure blitz2006ConstructorInternal
  end interface starFormationRateSurfaceDensityDisksBlitz2006

contains

  function blitz2006ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily blitz2006} star formation surface density rate in disks class which takes a parameter set as input.
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent

    !# <inputParameter>
    !#   <name>velocityDispersionDiskGas</name>
    !#   <defaultSource>\citep{leroy_star_2008}</defaultSource>
    !#   <defaultValue>10.0d0</defaultValue>
    !#   <description>The velocity dispersion of gas in disks.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>heightToRadialScaleDisk</name>
    !#   <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
    !#   <defaultValue>0.137d0</defaultValue>
    !#   <description>The ratio of scale height to scale radius for disks in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>surfaceDensityCritical</name>
    !#   <defaultSource>\citep{bigiel_star_2008}</defaultSource>
    !#   <defaultValue>200.0d0</defaultValue>
    !#   <description>The surface density (in units of $M_\odot$ pc$^{-2}$) in the ``Blitz-Rosolowsky2006'' star formation timescale calculation at which low-density truncation begins.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>surfaceDensityExponent</name>
    !#   <defaultSource>\citep{bigiel_star_2008}</defaultSource>
    !#   <defaultValue>0.4d0</defaultValue>
    !#   <description>The exponent for surface density in the ``Blitz-Rosolowsky2006'' star formation timescale calculation at in the high density regime.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>starFormationFrequencyNormalization</name>
    !#   <defaultSource>\citep{leroy_star_2008}</defaultSource>
    !#   <defaultValue>5.25d-10</defaultValue>
    !#   <description>The star formation frequency (in the low-density limit and in units of yr$^{-1}$) in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>pressureCharacteristic</name>
    !#   <defaultSource>\citep{blitz_role_2006}</defaultSource>
    !#   <defaultValue>4.54d0</defaultValue>
    !#   <description>The characteristic pressure (given as $P_0/k_\mathrm{B}$ in units of K cm$^{-3}$) in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>pressureExponent</name>
    !#   <defaultSource>\citep{blitz_role_2006}</defaultSource>
    !#   <defaultValue>0.92d0</defaultValue>
    !#   <description>The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=starFormationRateSurfaceDensityDisksBlitz2006(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent)
    !# <inputParametersValidate source="parameters"/>
    return
  end function blitz2006ConstructorParameters

  function blitz2006ConstructorInternal(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent) result(self)
    !% Internal constructor for the {\normalfont \ttfamily blitz2006} star formation surface density rate from disks class.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : giga                   , hecto     , kilo, mega
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    double precision                                               , intent(in   ) :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent
     !# <constructorAssign variables="velocityDispersionDiskGas, heightToRadialScaleDisk, surfaceDensityCritical, surfaceDensityExponent, starFormationFrequencyNormalization, pressureCharacteristic, pressureExponent"/>

    self%lastUniqueID   =-1_kind_int8
    self%factorsComputed=.false.
    ! Validate
    if (pressureExponent < 0.0d0) call Galacticus_Error_Report('pressureExponent < 0 violates assumptions'//{introspection:location})
    ! Convert parameters to internal units.
    self%surfaceDensityCritical             =self%surfaceDensityCritical*(mega**2)                                                    ! Convert to M☉/Mpc².
    self%starFormationFrequencyNormalization=self%starFormationFrequencyNormalization*giga                                            ! Convert to Gyr⁻¹.
    self%pressureCharacteristic             =self%pressureCharacteristic*boltzmannsConstant*((hecto*megaParsec)**3)/massSolar/kilo**2 ! Convert to M☉(km/s)²/Mpc.
    return
  end function blitz2006ConstructorInternal

  subroutine blitz2006AutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    call calculationResetEvent%attach(self,blitz2006CalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine blitz2006AutoHook

  subroutine blitz2006Destructor(self)
    !% Destructor for the blitz2006 cooling radius class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    call calculationResetEvent%detach(self,blitz2006CalculationReset)
    return
  end subroutine blitz2006Destructor

  subroutine blitz2006CalculationReset(self,node)
    !% Reset the Kennicutt-Schmidt relation calculation.
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node

    self%factorsComputed=.false.
    self%lastUniqueID   =node%uniqueID()
    return
  end subroutine blitz2006CalculationReset

  double precision function blitz2006Rate(self,node,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    !% in the galactic disk of {\normalfont \ttfamily node}. The disk is assumed to obey the
    !% \cite{blitz_role_2006} star formation rule.
    use :: Abundances_Structure                , only : abundances
    use :: Galactic_Structure_Options          , only : componentTypeDisk                 , coordinateSystemCylindrical, massTypeGaseous, massTypeStellar
    use :: Galactic_Structure_Surface_Densities, only : Galactic_Structure_Surface_Density
    use :: Galacticus_Nodes                    , only : nodeComponentDisk                 , treeNode
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Constants_Astronomical        , only : gravitationalConstantGalacticus
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: radius
    class           (nodeComponentDisk                            ), pointer       :: disk
    type            (abundances                                   ), save          :: abundancesFuel
    !$omp threadprivate(abundancesFuel)
    double precision                                                               :: molecularFraction, pressureRatio        , &
         &                                                                            surfaceDensityGas, surfaceDensityStellar

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if factors have been precomputed.
    if (.not.self%factorsComputed) then
       ! Get the disk properties.
       disk             => node%disk       ()
       self%massGas     =  disk%massGas    ()
       self%massStellar =  disk%massStellar()
       self%radiusDisk  =  disk%radius     ()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       abundancesFuel=disk%abundancesGas()
       call abundancesFuel%massToMassFraction(self%massGas)
       self%hydrogenMassFraction=abundancesFuel%hydrogenMassFraction()
       ! Record that factors have now been computed.
       self%factorsComputed=.true.
    end if
    ! Return zero rate for non-positive radius or mass.
    if (self%massGas <= 0.0d0 .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
       blitz2006Rate=0.0d0
       return
    end if
    ! Get gas and stellar surface densities.
    surfaceDensityGas    =Galactic_Structure_Surface_Density(node,[radius,0.0d0,0.0d0],coordinateSystem=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    surfaceDensityStellar=Galactic_Structure_Surface_Density(node,[radius,0.0d0,0.0d0],coordinateSystem=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeStellar)
    ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.
    pressureRatio=+0.5d0                                   &
         &        *Pi                                      &
         &        *gravitationalConstantGalacticus         &
         &        *surfaceDensityGas                       &
         &        *(                                       &
         &          +surfaceDensityGas                     &
         &          +self%velocityDispersionDiskGas        &
         &          *sqrt(                                 &
         &                +surfaceDensityStellar           &
         &                /Pi                              &
         &                /gravitationalConstantGalacticus &
         &                /self%heightToRadialScaleDisk    &
         &                /self%radiusDisk                 &
         &               )                                 &
         &         )                                       &
         &        /self%pressureCharacteristic
    ! Compute the molecular fraction, limited to 100% molecular.
    if (pressureRatio >= 1.0d0) then
       molecularFraction=                                         1.0d0
    else
       molecularFraction=min(pressureRatio**self%pressureExponent,1.0d0)
    end if
    ! Compute the star formation rate surface density.
    blitz2006Rate=+surfaceDensityGas                        &
         &        *self%hydrogenMassFraction                &
         &        *molecularFraction                        &
         &        *self%starFormationFrequencyNormalization &
         &        *(                                        &
         &          +1.0d0                                  &
         &          +(                                      &
         &            +self%hydrogenMassFraction            &
         &            *surfaceDensityGas                    &
         &            /self%surfaceDensityCritical          &
         &           )**self%surfaceDensityExponent         &
         &         )
    return
  end function blitz2006Rate
