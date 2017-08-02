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

  !% An implementation of virial orbits which assumes fixed orbital parameters.

  use Dark_Matter_Halo_Scales
  
  !# <virialOrbit name="virialOrbitFixed">
  !#  <description>Virial orbits assuming fixed orbital parameters.</description>
  !# </virialOrbit>
  type, extends(virialOrbitClass) :: virialOrbitFixed
     !% A virial orbit class that assumes fixed orbital parameters.
     private
     double precision                                      :: velocityRadial                  , velocityTangential
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_   => null()
   contains
     final     ::                              fixedDestructor
     procedure :: orbit                     => fixedOrbit
     procedure :: densityContrastDefinition => fixedDensityContrastDefinition
  end type virialOrbitFixed
  
  interface virialOrbitFixed
     !% Constructors for the {\normalfont \ttfamily fixed} virial orbit class.
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface virialOrbitFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fixed} satellite virial orbit class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (virialOrbitFixed          )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    class           (darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    double precision                                            :: velocityRadial        , velocityTangential

    !# <inputParameter>
    !#   <name>velocityRadial</name>
    !#   <defaultValue>-0.90d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>velocityTangential</name>
    !#   <defaultValue>0.75d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="virialDensityContrast"  name="virialDensityContrast_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"    name="darkMatterHaloScale_"    source="parameters"/>
    self=virialOrbitFixed(velocityRadial,velocityTangential,virialDensityContrast_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(velocityRadial,velocityTangential,virialDensityContrast_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily fixed} virial orbits class.
    use Galacticus_Error
    implicit none
    type            (virialOrbitFixed          )                        :: self
    class           (virialDensityContrastClass), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterHaloScaleClass  ), intent(in   ), target :: darkMatterHaloScale_
    double precision                            , intent(in   )         :: velocityRadial        , velocityTangential
    !# <constructorAssign variables="velocityRadial, velocityTangential, *virialDensityContrast_, *darkMatterHaloScale_"/>

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily fixed} virial orbits class.
    implicit none
    type(virialOrbitFixed), intent(inout) :: self

    !# <objectDestructor name="self%virialDensityContrast_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"  />
    return
  end subroutine fixedDestructor

  function fixedOrbit(self,node,host,acceptUnboundOrbits)
    !% Return fixed orbital parameters for a satellite.
    use Galacticus_Error
    use Galacticus_Display
    use Dark_Matter_Profile_Mass_Definitions
    use ISO_Varying_String
    implicit none
    type            (keplerOrbit       )                        :: fixedOrbit
    class           (virialOrbitFixed  ), intent(inout), target :: self
    type            (treeNode          ), intent(inout)         :: host               , node
    logical                             , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic), pointer               :: hostBasic          , basic
    double precision                                            :: velocityHost       , massHost     , &
         &                                                         radiusHost         , massSatellite
    type            (varying_string    )                        :: message
    character       (len=12            )                        :: label
    !GCC$ attributes unused :: acceptUnboundOrbits
    
    ! Reset the orbit.
    call fixedOrbit%reset()
    ! Get required objects.
    basic     => node%basic()
    hostBasic => host%basic()
    ! Find mass, radius, and velocity in the host corresponding to the our virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%time()),radiusHost,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,self%virialDensityContrast_%densityContrast(    basic%mass(),    basic%time())                        )
    ! Set basic properties of the orbit - do not allow the satellite mass to exceed the host mass.
    call fixedOrbit%massesSet(min(massSatellite,massHost),massHost)
    call fixedOrbit%radiusSet(radiusHost)
    call fixedOrbit%velocityRadialSet    (self%velocityRadial    *velocityHost)
    call fixedOrbit%velocityTangentialSet(self%velocityTangential*velocityHost)
    ! Propagate the orbit to the virial radius under the default density contrast definition.
    radiusHost=self%darkMatterHaloScale_%virialRadius(host)
    if (fixedOrbit%radiusApocenter() >= radiusHost .and. fixedOrbit%radiusPericenter() <= radiusHost) then
       call fixedOrbit%propagate(radiusHost,infalling=.true.)
       call fixedOrbit%massesSet(min(basic%mass(),hostBasic%mass()),hostBasic%mass())
    else
       call Galacticus_Verbosity_Level_Set(verbosityStandard)
       call Galacticus_Display_Indent('Satellite node')
       call node%serializeASCII()
       call Galacticus_Display_Indent('Host node'     )
       call host%serializeASCII()
       call Galacticus_Display_Indent('Host node'     )
       message="orbit does not reach halo radius"               //char(10)
       write (label,'(e12.6)') massSatellite
       message=message//"      satellite mass = "//label//" M☉" //char(10)
       write (label,'(e12.6)') massHost
       message=message//"           host mass = "//label//" M☉" //char(10)
       write (label,'(e12.6)') radiusHost
       message=message//"         host radius = "//label//" Mpc"//char(10)
       write (label,'(e12.6)') fixedOrbit%radiusPericenter()
       message=message//"    orbit pericenter = "//label//" Mpc"//char(10)
       write (label,'(e12.6)') fixedOrbit%radiusApocenter()
       message=message//"     orbit apocenter = "//label//" Mpc"
       call Galacticus_Error_Report('fixedOrbit',message)
    end if
    return
  end function fixedOrbit

  function fixedDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of fixed virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: fixedDensityContrastDefinition
    class(virialOrbitFixed          ), intent(inout) :: self

    fixedDensityContrastDefinition => self%virialDensityContrast_
    return
  end function fixedDensityContrastDefinition

