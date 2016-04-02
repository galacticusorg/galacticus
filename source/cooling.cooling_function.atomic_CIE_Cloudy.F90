!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Implements a cooling function class which utilizes the {\normalfont \scshape Cloudy} code to
  !% compute cooling in collisional ionization equilibrium.
  
  !# <coolingFunction name="coolingFunctionAtomicCIECloudy" defaultThreadPrivate="yes">
  !#  <description>
  !#   Class providing cooling function by utilizing the {\normalfont \scshape Cloudy} code to
  !#   compute cooling in collisional ionization equilibrium. {\normalfont \scshape Cloudy} will
  !#   be downloaded, compiled and run automatically if necessary\footnote{{\normalfont \scshape
  !#   Cloudy} is used to generate a file which contains a tabulation of the cooling function
  !#   suitable for reading by the {\normalfont \ttfamily CIE from file} method. Generation of
  !#   the tabulation typically takes several hours, but only needs to be done once as the
  !#   stored table is simply read back in on later runs.}.  
  !#  </description>
  !# </coolingFunction>
  type, extends(coolingFunctionCIEFile) :: coolingFunctionAtomicCIECloudy
     !% A cooling function class which utilizes the {\normalfont \scshape Cloudy} code to compute
     !% the cooling function in collisional ionization equilibrium.
     private
     logical :: initialized
   contains
     !@ <objectMethods>
     !@   <object>coolingFunctionAtomicCIECloudy</object>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless type(abundances)\textgreater} gasAbunances\argin</arguments>
     !@     <description>Run {\normalfont \scshape Cloudy} to tabulate the cooling function as necessary.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                       atomicCIECloudyDestructor
     procedure :: tabulate                           => atomicCIECloudyTabulate
     procedure :: coolingFunction                    => atomicCIECloudyCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => atomicCIECloudyCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => atomicCIECloudyCoolingFunctionDensityLogSlope
     procedure :: descriptor                         => atomicCIECloudyDescriptor
  end type coolingFunctionAtomicCIECloudy

  interface coolingFunctionAtomicCIECloudy
     !% Constructors for the ``atomic CIE Cloudy'' cooling function class.
     module procedure atomicCIECloudyConstructorParameters
     module procedure atomicCIECloudyConstructorInternal
  end interface coolingFunctionAtomicCIECloudy

  ! File names for the cooling function and chemical state data.
  character       (len=55)           :: atomicCIECloudyCoolingFunctionFileName='data/cooling/cooling_function_Atomic_CIE_Cloudy.xml'
  character       (len=55)           :: atomicCIECloudyChemicalStateFileName  ='data/chemicalState/chemical_state_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision       , parameter :: atomicCIECloudyMetallicityMaximumDefault=30.0d0 ! Thirty times Solar.

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision       , parameter :: atomicCIECloudyMetallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

  ! Factor by which metallicity must exceed currently tabulated maximum before we retabulate.
  double precision       , parameter :: atomicCIECloudyMetallicityTolerance     = 0.1d0

contains

  function atomicCIECloudyConstructorParameters(parameters)
    !% Constructor for the ``atomic CIE Cloudy'' cooling function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(coolingFunctionAtomicCIECloudy)                :: atomicCIECloudyConstructorParameters
    type(inputParameters               ), intent(inout) :: parameters

    atomicCIECloudyConstructorParameters=atomicCIECloudyConstructorInternal()
    return
  end function atomicCIECloudyConstructorParameters
  
  function atomicCIECloudyConstructorInternal()
    !% Internal constructor for the ``atomic CIE Cloudy'' cooling function class.
    implicit none
    type(coolingFunctionAtomicCIECloudy) :: atomicCIECloudyConstructorInternal

    ! Initialize.
    atomicCIECloudyConstructorInternal%temperaturePrevious     =-1.0d0
    atomicCIECloudyConstructorInternal%metallicityPrevious     =-1.0d0
    atomicCIECloudyConstructorInternal%temperatureSlopePrevious=-1.0d0
    atomicCIECloudyConstructorInternal%metallicitySlopePrevious=-1.0d0
    atomicCIECloudyConstructorInternal%resetMetallicity        =.true.
    atomicCIECloudyConstructorInternal%resetTemperature        =.true.
    atomicCIECloudyConstructorInternal%initialized             =.false.
   return
  end function atomicCIECloudyConstructorInternal
  
  subroutine atomicCIECloudyDestructor(self)
    !% Destructor for the ``atomic CIE Cloudy'' cooling function class.
    use Numerical_Interpolation
    implicit none
    type(coolingFunctionAtomicCIECloudy), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine atomicCIECloudyDestructor

  subroutine atomicCIECloudyTabulate(self,gasAbundances)
    !% Create the chemical state.
    use System_Command
    use Galacticus_Input_Paths
    use String_Handling
    use Galacticus_Display
    implicit none
    class    (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type     (abundances                    ), intent(in   ) :: gasAbundances
    logical                                                  :: makeFile
    character(len=32                        )                :: metallicityLabel
    type     (varying_string                )                :: command
    integer                                                  :: status
    
    ! Determine if we need to retabulate the chemical state.
    if (.not.self%initialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) > 0.0d0) then
          makeFile=(                                                                      &
               &     min(                                                                 &
               &         Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar)    , &
               &         atomicCIECloudyMetallicityMaximumLimit                           &
               &        )                                                                 &
               &    >                                                                     &
               &     self%metallicityMaximum*(1.0d0+atomicCIECloudyMetallicityTolerance)  &
               &   )
       else
          makeFile=.false.
       end if
       if (makeFile) then
          ! Remove the cooling function file so that a new one will be created.
          command='rm -f '//char(Galacticus_Input_Path())//trim(atomicCIECloudyCoolingFunctionFileName)
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the metallicity is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) > 0.0d0) then
          self%metallicityMaximum=min(                                                                 &
               &                      max(                                                             &
               &                           atomicCIECloudyMetallicityMaximumDefault                  , &
               &                          +3.0d0                                                       &
               &                          *Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) &
               &                         )                                                           , &
               &                           atomicCIECloudyMetallicityMaximumLimit                      &
               &                     )
       else
          self%metallicityMaximum=atomicCIECloudyMetallicityMaximumDefault
       end if
       write (metallicityLabel,'(e12.6)') log10(self%metallicityMaximum)
       ! Test if we can compile the Cloudy driver script.
       command='perl -c '//char(Galacticus_Input_Path())//'scripts/aux/Atomic_CIE_Cloudy_Driver.pl &> /dev/null'
       call System_Command_Do(command,status)
       if (status == 0) then       
          ! Run Atomic_CIE_Cloudy wrapper script.
          command=char(Galacticus_Input_Path())//'scripts/aux/Atomic_CIE_Cloudy_Driver.pl'   //' '// &
               &  metallicityLabel                                                           //' '// &
               &  char(Galacticus_Input_Path())//trim(atomicCIECloudyCoolingFunctionFileName)//' '// &
               &  char(Galacticus_Input_Path())//trim(atomicCIECloudyChemicalStateFileName  )
          command=command//" "//cieFileFormatVersionCurrent
          call System_Command_Do(command)
       end if
       ! Call routine to read in the tabulated data.
       call self%readFile(char(Galacticus_Input_Path()//trim(atomicCIECloudyCoolingFunctionFileName)))
       ! Flag that cooling function is now initialized.
       self%initialized=.true.
    end if
    return
  end subroutine atomicCIECloudyTabulate
  
  double precision function atomicCIECloudyCoolingFunction(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function for collisional ionization equilibrium as computed by
    !% {\normalfont \scshape Cloudy}.
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    type            (radiationStructure            ), intent(in   ) :: radiation

    call                           self%tabulate                              (                                  gasAbundances                            )
    atomicCIECloudyCoolingFunction=self%coolingFunctionCIEFile%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunction

  double precision function atomicCIECloudyCoolingFunctionTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling function with respect to temperature for
    !% collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.  read
    !% from a file.
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    type            (radiationStructure            ), intent(in   ) :: radiation

    call                                              self%tabulate                                                 (                                  gasAbundances                            )
    atomicCIECloudyCoolingFunctionTemperatureLogSlope=self%coolingFunctionCIEFile%coolingFunctionTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunctionTemperatureLogSlope

  double precision function atomicCIECloudyCoolingFunctionDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling function with respect to density for
    !% collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    type            (radiationStructure            ), intent(in   ) :: radiation

    call                                          self%tabulate                                             (                                  gasAbundances                            )
    atomicCIECloudyCoolingFunctionDensityLogSlope=self%coolingFunctionCIEFile%coolingFunctionDensityLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunctionDensityLogSlope

  subroutine atomicCIECloudyDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type (inputParameters               ), intent(inout) :: descriptor
 
    call descriptor%addParameter("coolingFunctionMethod","atomicCIECloudy")
    return
  end subroutine atomicCIECloudyDescriptor
