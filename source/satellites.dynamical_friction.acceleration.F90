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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module that implements calculations of the acceleration due to dynamical friction for satellites.

module Satellite_Dynamical_Friction
  !% Implements calculations of dynamical friction for satellites.

  use Galacticus_Nodes, only : treeNode
  
  !# <functionClass>
  !#  <name>satelliteDynamicalFriction</name>
  !#  <descriptiveName>Dynamical friction models.</descriptiveName>
  !#  <description>Class providing models of dynamical friction.</description>
  !#  <default>chandrasekhar1943</default>
  !#  <method name="acceleration" >
  !#   <description>Returns the satellite acceleration due to dynamical friction for {\normalfont \ttfamily node} (in units of km/s/Gyr).</description>
  !#   <type>double precision, dimension(3)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !# </functionClass>

end module Satellite_Dynamical_Friction
