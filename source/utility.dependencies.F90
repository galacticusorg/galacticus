!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements dependency versioning.
!!}

module Dependencies
  !!{
  Implements dependency versioning.
  !!}
  use :: Hashes, only : varyingStringHash
  implicit none
  private
  public :: dependencyVersion

  type   (varyingStringHash) :: dependencies_
  logical                    :: initialized =.false.
  
contains

  function dependencyVersion(dependency,majorOnly)
    !!{
    Return the version number to use for a named dependency.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath     , pathTypeExec
    use :: ISO_Varying_String, only : varying_string, trim        , assignment(=), char
    implicit none
    type     (varying_string)                          :: dependencyVersion
    character(len=*         ), intent(in   )           :: dependency
    logical                  , intent(in   ), optional :: majorOnly
    integer                                            :: fileUnit
    character(len=256       )                          :: line             , dependency_
    type     (varying_string)                          :: version_
    integer                                            :: indexSeparator   , status
    !![
    <optionalArgument name="majorOnly" defaultsTo=".false."/>
    !!]
    
    !$omp critical(dependenciesInitialize)
    if (.not.initialized) then
       call dependencies_%initialize()
       open(newUnit=fileUnit,file=char(inputPath(pathTypeExec))//'aux/dependencies.yml',status='old',form='formatted',iostat=status)
       do while (status == 0)
          read (fileUnit,'(a)',iostat=status) line
          if (status /= 0) exit
          indexSeparator=index(line,":")
          if (indexSeparator == 0) then
             call Error_Report('badly-formed YAML'//{introspection:location})
          else
             dependency_=line(1:indexSeparator-1 )
             version_   =line(  indexSeparator+2:)
             call dependencies_%set(trim(dependency_),trim(version_))
          end if
       end do
       close(fileUnit)
       initialized=.true.
    end if
    !$omp end critical(dependenciesInitialize)
    if (dependencies_%exists(trim(dependency))) then
       dependencyVersion=dependencies_%value(trim(dependency))
       if (majorOnly_ .and. index(dependencyVersion,".") > 1) &
            & dependencyVersion=extract(dependencyVersion,1,index(dependencyVersion,".")-1)
    else
       call Error_Report('dependency not found'//{introspection:location})
    end if
    return
  end function dependencyVersion
  
end module Dependencies
