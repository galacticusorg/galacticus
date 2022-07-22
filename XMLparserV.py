import xml.etree.ElementTree as ET
import subprocess

def changeParameters(tau, v, heat, mloss, halo_mass):
    # change parameters
    tree = ET.parse('isolatedHaloParametersDDM.xml')
    root = tree.getroot()
    particle = root.find('darkMatterParticle')
    lifetime = particle.find('lifetime')
    massSplitting = particle.find('massSplitting')
    heating = particle.find('heating')
    massLoss = particle.find('massLoss')
    lifetime.set('value', tau)
    massSplitting.set('value', str(v/300000))
    heating.set('value', heat)
    massLoss.set('value', mloss)
    file_name = 'results/T' + tau + 'V' + str(v) + 'H' + heat + 'ML' + mloss + 'HM' + halo_mass + '.hdf5'
    name = root.find('outputFileName')
    name.set('value', file_name)
    tree.write('isolatedHaloParametersDDM.xml')
    # change halo tree
    tree = ET.parse('isolatedHaloTree.xml')
    root = tree.getroot()
    for node in root.iter('node'):
        basic = node.find('basic')
        mass = basic.find('mass')
        mass.text = halo_mass
    tree.write('isolatedHaloTree.xml')

if __name__ == '__main__':
    lifetimes = ['100', '50', '10', '1', '0.1']
    vel = [10, 50, 100, 200, 500]
    for lt in lifetimes:
        for v in vel:
            for h in ['true', 'false']:
                for ml in ['true', 'false']:
                    changeParameters(lt, v, h, ml, '1.0e12')
                    subprocess.run(["./Galacticus.exe", "isolatedHaloParametersDDM.xml"])