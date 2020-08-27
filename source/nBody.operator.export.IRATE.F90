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

!% Contains a module which implements an N-body data operator which exports N-body data to IRATE format.

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Cosmology_Functions , only : cosmologyFunctionsClass

  !# <nbodyOperator name="nbodyOperatorExportIRATE">
  !#  <description>An N-body data operator which exports N-body data to IRATE format.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorExportIRATE
     !% An N-body data operator which exports N-body data to IRATE format.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     type            (varying_string          )          :: fileName
     integer                                             :: snapshot
     double precision                                    :: redshift
   contains
     final     ::            exportIRATEDestructor
     procedure :: operate => exportIRATEOperate
  end type nbodyOperatorExportIRATE

  interface nbodyOperatorExportIRATE
     !% Constructors for the {\normalfont \ttfamily exportIRATE} N-body operator class.
     module procedure exportIRATEConstructorParameters
     module procedure exportIRATEConstructorInternal
  end interface nbodyOperatorExportIRATE
    
contains
  
  function exportIRATEConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily exportIRATE} N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorExportIRATE)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    type            (varying_string          )                :: fileName
    integer                                                   :: snapshot
    double precision                                          :: redshift

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file to which data should be exported.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>snapshot</name>
    !#   <source>parameters</source>
    !#   <description>The snapshot index of the data.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <description>The redshift of the data.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=nbodyOperatorExportIRATE(fileName,snapshot,redshift,cosmologyParameters_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    return
  end function exportIRATEConstructorParameters

  function exportIRATEConstructorInternal(fileName,snapshot,redshift,cosmologyParameters_,cosmologyFunctions_) result (self)
    !% Internal constructor for the {\normalfont \ttfamily exportIRATE} N-body operator class.
    implicit none
    type            (nbodyOperatorExportIRATE)                        :: self
    type            (varying_string          ), intent(in   )         :: fileName
    integer                                   , intent(in   )         :: snapshot
    double precision                          , intent(in   )         :: redshift
    class           (cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="fileName, snapshot, redshift, *cosmologyParameters_, *cosmologyFunctions_"/>

    return
  end function exportIRATEConstructorInternal

  subroutine exportIRATEDestructor(self)
    !% Destructor for {\normalfont \ttfamily exportIRATE} importer class.
    implicit none
    type(nbodyOperatorExportIRATE), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine exportIRATEDestructor

  subroutine exportIRATEOperate(self,simulations)
    !% Output simulation data to an IRATE-format file.
    use :: Galacticus_Display              , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, verbosityStandard
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: IO_HDF5                         , only : hdf5Object, hdf5Access
    use :: IO_IRATE                        , only : irate
    use :: ISO_Varying_String              , only : char
    use :: Numerical_Constants_Astronomical, only : massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (nbodyOperatorExportIRATE), intent(inout)               :: self
    type            (nBodyData               ), intent(inout), dimension(:) :: simulations
    type            (irate                   )                              :: irate_
    character       (len=13                  )                              :: snapshotLabel
    type            (hdf5Object              )                              :: irateFile         , snapshotGroup, &
         &                                                                     halosGroup        , dataset
    integer                                                                 :: i
    type            (varying_string          )                              :: datasetDescription, unitName
    double precision                                         , dimension(3) :: unitscgs

    call Galacticus_Display_Indent('export simulation to IRATE file',verbosityStandard)
    if (size(simulations) /= 1) call Galacticus_Error_Report('precisely 1 simulation should be supplied'//{introspection:location})
    irate_=irate(char(self%fileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call irate_%writeHalos(                                                &
         &                                     self      %snapshot       , &
         &                                     self      %redshift       , &
         &                 center             =simulations(1)%position   , &
         &                 velocity           =simulations(1)%velocity   , &
         &                 IDs                =simulations(1)%particleIDs, &
         &                 overwrite          =.true.                    , &
         &                 objectsOverwritable=.true.                      &
         &                )
    ! Write any additional properties to the file.
    if     (                                                 &
         &      simulations(1)%propertiesInteger%size()  > 0 &
         &  .or.                                             &
         &      simulations(1)%propertiesReal   %size()  > 0 &
         & ) then       
       write (snapshotLabel,'(a,i5.5)') 'Snapshot',self%snapshot
       !$ call hdf5Access%set()
       call irateFile%openFile(char(self%fileName),readOnly=.false.)
       snapshotGroup=irateFile    %openGroup(snapshotLabel)
       halosGroup   =snapshotGroup%openGroup('HaloCatalog')
       do i=1,simulations(1)%propertiesInteger%size()
          select case (char(simulations(1)%propertiesInteger%key(i)))
          case ('particleID'               )
             datasetDescription="Unique ID for each particle."
          case ('descendentID'             )
             datasetDescription="Unique ID of descendent."
          case ('progenitorCount'          )
             datasetDescription="Number of progenitors."
          case ('hostID'                   )
             datasetDescription="ID of immediate host."
          case ('hostRootID'               )
             datasetDescription="ID of root host."
          case ('descendentHostID'         )
             datasetDescription="ID of descendent's immediate host"
          case ('isPhantom'                )
             datasetDescription="Zero (0) for real particles, non-zero for phantom particles."
          case default
             datasetDescription="Unknown property."
          end select
          call halosGroup%writeDataset  (simulations(1)%propertiesInteger%value(i),char(simulations(1)%propertiesInteger%key(i)),char(datasetDescription)                        )
       end do
       do i=1,simulations(1)%propertiesReal   %size()
          select case (char(simulations(1)%propertiesReal   %key(i)))
          case ('expansionFactor'          )
             datasetDescription="Cosmological expansion factor."
             unitName="dimensionless"
             unitscgs=0.0d0
          case ('descendentExpansionFactor')
             datasetDescription="Cosmological expansion factor of descendent."
             unitName="dimensionless"
             unitscgs=0.0d0
          case ('massVirial'               )
             datasetDescription="Halo virial mass."
             unitName="massSolar"
             unitscgs=[massSolar*kilo,0.0d0,0.0d0]
          case default
             datasetDescription="Unknown property."
          end select
          call halosGroup%writeDataset  (simulations(1)%propertiesReal   %value(i),char(simulations(1)%propertiesReal   %key(i)),char(datasetDescription),datasetReturned=dataset)
          call dataset   %writeAttribute(char(unitName)                           ,'unitname'                                                                                    )
          call dataset   %writeAttribute(     unitscgs                            ,'unitscgs'                                                                                    )
          call dataset   %close         (                                                                                                                                        )
       end do
       call halosGroup   %close()
       call snapshotGroup%close()
       call irateFile    %close()
       !$ call hdf5Access%unset()
    end if
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine exportIRATEOperate
