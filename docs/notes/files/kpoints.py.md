# kpoints.py
The `kpoints.py` file can be used to generate kpoints for all possible types of calculations that we do in the Marom group. To use it, copy the code block below and put it in `~/bin/kpoints.py`, then run the following command to add a <a href="https://en.wikipedia.org/wiki/Shebang_(Unix)" target="_blank">Shebang</a> line to the top of the `kpoints.py` file which tells your system which python interpreter to use when you make the file executable. 

```bash
sed -i '1s,.*,'"#\!$(which python)"',' ~/bin/kpoints.py
```

This should append something like the following to the top of the `kpoints.py` file:

```txt
#!/global/homes/d/ddardzin/.local/mambaforge/bin/python
```

Now you can make the `kpoints.py` file and executable file by running the following command.

```bash
chmod +x ~/bin/incar.py
```


```python
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.ase import AseAtomsAdaptor
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.io import read
import os
import argparse
import warnings

def parse_arguments():
    parser = argparse.ArgumentParser(description='param')
    parser.add_argument(
        '-k',
        '--kpoints',
        help="File path to KPOINTS file",
        type=str,
    )
    parser.add_argument(
        '-i',
        '--ibzkpt',
        help="File path to IBZKPT file",
        default='IBZKPT',
        type=str,
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Output file name",
        default='KPOINTS',
        type=str,
    )
    parser.add_argument(
        '-p',
        '--poscar',
        help="POSCAR file name, used for finding the band structure k-path",
        default='POSCAR',
        type=str,
    )
    parser.add_argument(
        '-e',
        '--hse',
        help="Generates a KPOINTS file for an HSE band calculation",
        action='store_true'
    )
    parser.add_argument(
        '-b',
        '--band',
        help="Generates a general Line-mode KPOINTS file",
        action='store_true'
    )
    parser.add_argument(
        '-g',
        '--grid',
        help="Generates a general mesh grid KPOINTS file",
        action='store_true'
    )
    parser.add_argument(
        '-ad',
        '--autodensity',
        help="k-point density to be converted to a Monkhorst-Pack grid",
        type=float
    )
    parser.add_argument(
        '-d',
        '--density',
        nargs='+',
        help="Grid density of a grid KPOINTS file (e.g. -d 7 7 7 or --density 7 7 7)",
        default=[7, 7, 7],
        type=int
    )
    parser.add_argument(
        '-n',
        '--segments',
        help="Number of segements between each high symmetry point for a band path",
        default=50,
        type=int
    )
    parser.add_argument(
        '-c',
        '--coords',
        help="String of high symmetry points for the band structure (e.g. -c GXWLGK)",
        type=str
    )
    parser.add_argument(
        '-l',
        '--list',
        help="This option will print out a list of the high symmetry points and their labels of a given POSCAR",
        action='store_true'
    )
    return parser.parse_args()

args = parse_arguments()
output_path = args.output
kpoints_path = args.kpoints
ibzkpt_path = args.ibzkpt
hse = args.hse
band = args.band
grid = args.grid
density = args.density
autodensity = args.autodensity
poscar = args.poscar
segments = args.segments
coords = args.coords
list_coords = args.list


def get_high_symm_points(poscar_path, coords):
    if os.path.isfile(poscar):
        st = Structure.from_file(poscar)
        path_types = ["latimer_munro", "setyawan_curtarolo", "hinuma"]
        kp = HighSymmKpath(st, path_type=path_types[1])
        print(kp)
        kpoints = Kpoints.automatic_linemode(segments, kp)
        kpoints.labels = ';'.join(kpoints.labels).replace('\\Gamma', 'G').split(';')
           
        if coords is not None:
            coords = coords.upper()
            if all(c in kpoints.labels for c in coords):
                kpath = [c for c in coords]
                
                for i in reversed(range(len(kpath))):
                    if i < len(kpath) - 1 and i > 0:
                        kpath.insert(i+1, kpath[i])
                        
                info_dict = {l: c for l, c in zip(
                    kpoints.labels, kpoints.kpts)}
                    
                kpoints.labels = kpath
                kpoints.kpts = [info_dict[l] for l in kpath]
            else:
                warnings.warn(
                    f'The high-symmetry points ({coords}) were not found in the generated BZ, please check the output KPOINTS file for the correct labels.', stacklevel=3)
                    
        return kpoints
    else:
        raise BaseException(
            'The POSCAR file does not exist, please make sure the file name is correct.')


def interpolate(point1, point2, n):
    xs = np.linspace(point1[0], point2[0], n)
    ys = np.linspace(point1[1], point2[1], n)
    zs = np.linspace(point1[2], point2[2], n)
    
    return np.c_[xs, ys, zs]
    
if band:
    if hse:
        if os.path.isfile(ibzkpt_path):
            with open(ibzkpt_path, 'r') as ibzkpt_file:
                ibzkpt = ibzkpt_file.read()
                
            split_ibzkpt = ibzkpt.split('\n')
            
            ibzkpt_nb = int(split_ibzkpt[1].strip())
        else:
            raise BaseException(
                f'The IBZKPT file cannot be found at {ibzkpt_path}')
                
        kpoints = get_high_symm_points(
            poscar_path=poscar,
            coords=coords,
        )
        
        num_kpts = kpoints.num_kpts
        kpts = kpoints.kpts
         
        all_kpoints = np.vstack(
            [interpolate(kpts[i], kpts[i+1], num_kpts)
             for i in range(len(kpts) - 1) if not i % 2],
        )
        
        total_kpts = ibzkpt_nb + all_kpoints.shape[0]
        split_ibzkpt[1] = f'\t{int(total_kpts)}'
        
        strings = [
            '\t'.join([''] + ['%.14f' % x for x in line] + ['\t\t 0']) for line in all_kpoints
        ]
        
        output = '\n'.join(split_ibzkpt) + '\n'.join(strings)
    else:
        try:
            kpoints = get_high_symm_points(
                poscar_path=poscar,
                coords=coords,
            )
            output = kpoints.__str__()
        except BaseException as x:
            raise
            
if grid:
    if autodensity:
        atoms = read(poscar)
        kpts = kptdensity2monkhorstpack(
            atoms,
            kptdensity=autodensity,
            even=False
        )
        d = f'{kpts[0]} '
        for i in kpts[1:]:
            d += f'{i} '
    else:
        if len(density) != 3 and len(density) != 1:
            print('The grid density should be 3 elements long, automatically using a density of 7 7 7. See kpoints.py -h for an example')
            d = '7 7 7 '
        elif len(density) == 3:
            d = f'{density[0]} '
            for i in density[1:]:
                d += f'{i} '
                
    output = f"Automatic kpoint scheme\n0\nGamma\n{d}\n"
    
if band or grid or hse:
    with open(output_path, 'w') as out:
        out.write(output)
elif list_coords:
    kpoints = get_high_symm_points(
        poscar_path=poscar,
        coords=coords,
    )
    labels = np.array(kpoints.labels)
    kpts = np.array(kpoints.kpts)   
    unique_labels, inds = np.unique(labels, return_index=True)
    unique_kpts = kpts[inds]
    
    for l, k in zip(unique_labels, unique_kpts):
        strings = ', '.join(['%.3f' % x for x in k])
        print(f"{l} ---> [{strings}]")
else:
    raise BaseException(
        "NO FILE WRITTEN: Please make sure to select an option (hse, band, grid, list), runkpoints.py -h for more information.")
(
```