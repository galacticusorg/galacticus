import xml.etree.ElementTree as ET
import subprocess

def changeParameters(tau, eps, heat, mloss, halo_mass):
    # change parameters
    tree = ET.parse('isolatedHaloParametersDDM.xml')
    root = tree.getroot()
    particle = root.find('darkMatterParticle')
    lifetime = particle.find('lifetime')
    massSplitting = particle.find('massSplitting')
    heating = particle.find('heating')
    massLoss = particle.find('massLoss')
    lifetime.set('value', tau)
    massSplitting.set('value', eps)
    heating.set('value', heat)
    massLoss.set('value', mloss)
    eps = eps[2:]
    file_name = 'results/T' + tau + 'MS' + eps + 'H' + heat + 'ML' + mloss + 'HM' + halo_mass + '.hdf5'
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
    lifetimes = ['10', '20']
    massSplitting = ['0.0001', '0.001', '0.01', '0.1']
    halo_mass = ['1.0e8', '1.0e9', '1.0e10', '1.0e11', '1.0e12']
    for hm in halo_mass:
        for lt in lifetimes:
            for ms in massSplitting:
                changeParameters(lt, ms, 'false', 'true', hm)
                subprocess.run(["./Galacticus.exe", "isolatedHaloParametersDDM.xml"])
