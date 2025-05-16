#optim
from ase.io import read
from gpaw import GPAW, PW
from ase.optimize import BFGS
import ase
import gpaw
import os
from ase.filters import FrechetCellFilter



def run_gpaw(input_cif, planewave, kpts, output_folder):

    name = ''
    for i in range(1,len(input_cif)):
        name = name + input_cif[-i]
        if input_cif[-i] == '/':
            break
    name = input_cif[:-4]
    monx = ase.io.read(input_cif)
    calc_PBE = GPAW(mode=PW(planewave),xc='PBE',kpts=kpts)
    monx.set_calculator(calc_PBE)
    box = FrechetCellFilter(monx, mask = [1,1,0,1,1,1],exp_cell_factor=20)

    dyn = BFGS(box)
    dyn.run(fmax=0.05)
    monx.write(output_folder +'/'+ name)
    monx.write(output_folder +'/'+ name[:-4]+'.json')
    return

def run_if_error(input_cif, planewave, kpts, output_folder):
    try: 
        run_gpaw(input_cif, planewave[0], kpts, output_folder)
    except (ValueError, gpaw.KohnShamConvergenceError):
        try:
            run_gpaw(input_cif, planewave[1], kpts, output_folder)
        except (ValueError, gpaw.KohnShamConvergenceError):
            try: 
                run_gpaw(input_cif, planewave[2], kpts, output_folder)
            except (ValueError, gpaw.KohnShamConvergenceError):
                pass

#=======================================================================================================
planewave = [900,800,700]
kpts = [10,10,1]

run_if_error(file, planewave, kpts, './results')
