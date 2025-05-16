from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
import ase
from ase.spectrum.band_structure import BandStructure

# Perform standard ground state calculation (with plane wave basis)
si = ase.io.read(file)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(10, 10, 1),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01))
si.calc = calc
si.get_potential_energy()
ef = calc.get_fermi_level()
calc.write(file[:-5]+'.gpw')

lat = si.cell.get_bravais_lattice()
pathlist = list(lat.get_special_points())
l = len(pathlist)
npoints = l * 20
path = ''
for i in range(len(pathlist)):
    path += pathlist[i]
# Restart from ground state and fix potential:
calc = GPAW(file[:-5]+'.gpw').fixed_density(
    nbands=60,
    symmetry='off',
    kpts={'path': path, 'npoints': npoints})
ef = calc.get_fermi_level()
bs = calc.band_structure()
energies = bs.energies - ef
path = bs.path
bss = BandStructure(path = path, energies=energies)



# bs.plot(filename=file[:-5]+'_bs1.png', emax = ef +2, emin = ef -2)
from gpaw import GPAW, restart
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 5)


slab, calc = restart(file[:-5]+'.gpw')
e, dos = calc.get_dos(spin=0, npts=2001, width=0.2)
ef = calc.get_fermi_level()

ed = []
dosd = []
for i in range(len(dos)):
    if e[i] > ef - 2 and e[i] < ef + 2:
        ed.append(e[i])
        dosd.append(dos[i])



f, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[5,1], layout = 'compressed')
plt.subplot(1,2,1)
ax1.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='minor', labelsize=15)

bss.plot(ax = ax1, emax = 2, emin = -2)
ax1.set_ylabel('$E - E_F$ (eV)', fontsize = 20)

plt.subplot(1,2,2)
ax2.plot(dosd, ed - ef)
ax2.set_ylim(-2,2)
ax2.axhline(y=0, color='k', linestyle='dotted')
# ax2.plot(dos_reduce, np.ones(len(dos_reduce)) * e_f, 'k.')
ax2.set_xlabel('DOS')
# plt.ylabel('Energy (eV)')
ax2.set_xticks([])
ax2.set_yticks([])
# ax2.set_aspect(3)
f.savefig('images/bs_dos_'+file[11:-5]+'.png')

