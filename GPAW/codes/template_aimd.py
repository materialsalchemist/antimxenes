import ase
import ase.md
from gpaw import GPAW, PW
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.md.nose_hoover_chain


atoms = ase.io.read(file)

# Describe the interatomic interactions with the PBE
atoms.calc = GPAW(mode=PW(500),xc='PBE',kpts=(1,1,1))

dyn = ase.md.nose_hoover_chain.NoseHooverChainNVT(atoms = atoms, timestep= 1, temperature_K=300, tdamp=100)


# def printenergy(a=atoms):  # store a reference to atoms in the definition.
#     """Function to print the potential, kinetic and total energy."""
#     epot = a.get_potential_energy() / len(a)
#     ekin = a.get_kinetic_energy() / len(a)
#     print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
#           'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

dyn.run(5000)
