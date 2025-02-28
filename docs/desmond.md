
## Usage

``` py
from mdscribe.desmond import Multisim

# read template .msj and .cfg
md_msj = Multisim(template="desmond-md.msj")
md_cfg = Multisim(template="desmond-md.cfg")

with open(msj_file,"w") as msj:
    # modify desmond msj template
    md_msj.dot.simulate[-1].cfg_file = cfg_file_basename
    
    # Setting up restraints using the restraints keyword:
    # https://www.schrodinger.com/kb/332119
    if args.posres_force > 0.0:
        # print the restraints in the multisim log file
        md_msj.dot.simulate[-1].print_restraint = 'true'

        # add the new terms defined in "restraints.new" to existing restraints.
        # The default is restraints.existing = ignore which will 
        # delete existing terms before adding any new ones.
        # md_msj.dot.simulate[-1].restraints.existing = 'retain'

        md_msj.dot.simulate[-1].restraints.new = [
            {
                'name'              : 'posre_harm',
                'atoms'             : [ f'"{args.posres}"' ],
                'force_constants'   : [ args.posres_force, ] * 3,
            }
            ]
        # force constants in the x, y, and z direction

    # write to a .msj file
    md_msj.write(msj)

with open(cfg_file,"w") as cfg:
    # modify desmond cfg template
    md_cfg.dot.randomize_velocity.seed = random.randint(1000, 9999)
    md_cfg.dot.time = total_simulation_time
    md_cfg.dot.temperature = t_schedule
    md_cfg.dot.trajectory.interval = args.interval
    
    # write to a .cfg file
    md_cfg.write(cfg)
```

## Metadynamics

### Collective variable types

| type    | Description                  |  Default Width | Wall | Floor | Atom sites |
| ----    | ---------------------------- | ------------- | ---- | ----- | ---------- |
| dist    |    Distance                               | 0.05 Å     |     yes  |   yes  |  2 |
| angle   |    Angle                                  | 2.5°       |    yes   |  yes   |  3 |
| dihedral  |   Dihedral                              |  5.0°      |      no  |    no  |  4 |
| rgyr      |  Radius of gyration                     | 0.1 Å      |     no   |   no   |  1 |
| rgyr_mass |  Mass-weighted radius of gyration       | 0.1 Å      |     no   |   no   |  1 |
| rmsd      |  RMSD from aligned starting structure   | 0.1 Å      |     no   |   no   |  1 |
| rmsd_symm |  Symmetry aware RMSD                    | 0.1 Å      |     no   |   no   |  1 |
| zdist     |  Distance along the z axis              | 0.05 Å     |     yes  |   yes  |  1 |
| zdist0    |  Absolute distance along the z axis     | 0.1 Å      |     yes  |   yes  |  1 |
| whim1     |  WHIM1 - first principal moment [35]    | 0.5 Å2     |     no   |   no   |  1 |
| whim2     |  WHIM2 - second principal moment [35]   | 0.25 Å2    |     no   |   no   |  1 |

### Choice of the bias factor (kTemp)

In well-tempered metadynamics, the height of the deployed Gaussians are rescaled (decreased) during simulation time by:
omega_0 * exp(-V/(kB * ΔT))
Where omega_0 is the initial hill height, V is the bias potential and the denominator in the exponential, kB * ΔT,  
is the bias factor (kTemp). During well-tempered metadynamics, the dynamics of the system are effectively accelerated 
(without heating up the system) up to T + ΔT, where T is the chosen MD simulation temperature. 
The choice of the bias factor value is guided by the highest barrier in the simulation system 
which the well-tempered metadynamics run should overcome. Here are some suggestions, 
assuming that the initial hill height ω0 (height parameter in the Metadynamics panel) has been set to 0.3 kcal/mol:

| Max. barrier height (kcal/mol) | kTemp (kcal/mol) |
| ------------------------------ | ---------------- |
| 3                              | 1.7              |
| 6                              | 3.4              |
| 10                             | 5.6              |
| 15                             | 8.4              |
| 20                             | 11.2             |

### Examples


``` py
distance:
    cv = [
        {atom = ["res. UNK" "res. 849" ]
        type = dist
        wall = 40
        floor = 10
        width = 0.05
        }
    ]

dihedral:
    cv = [
        {atom = [404 406 407 415 ]
        type = dihedral
        width = 5.0
        }
        {atom = [406 407 415 417 ]
        type = dihedral
        width = 5.0
        }
    ]

zdist(membrane):
    cv = [
        {atom = ["res. UNK"]
        type = zdist
        width = 0.05
        wall = 20
        floor = 5
        }
        ]

rmsd:
    cv = [
        {atom = ["res. UNK" "res. 849" ]
        type = rmsd
        width = 0.1
        }
    ]
```
    
### Extending Simulations

See https://www.schrodinger.com/kb/788642 for extending metadynamics simulation.


## Command-Line Interface (CLI)

`desmond.Multisim` class are used to build below CLI.

| Command-line interface | Description |
| ---------------------- | ----------- |
| `mdinfo`                 | Display a running Desmond MD simulation information in the Amber-style |
| `batch-desmond-setup`    | Setup/parametrize Desmond MD simulations |
| `batch-desmond-min`      | Run Desmond energy minimizations |
| `batch-desmond-md`       | Run Desmond MD simulations |
| `batch-desmond-metad`    | Run metadynamics MD simulations |
| `batch-desmond-pli`      | Analyze protein-ligand interactions |
| `batch-desmond-extend`   | Extend Desmond MD simulations |
| `batch-desmond-report`   | Generate reports for Desmond MD simulations |
| `batch-desmond-dihedral` | Analyze dihedral angles from Desmond MD trajectories |
| `batch-desmond-ligrmsd`  | Analyze ligand rmsd from Desmond MD trajectories |
| `batch-desmond-rg`       | Analyze radius of gyration from Desmond MD trajectories |

For more helps, `$ command --help`.

::: mdscribe.desmond.cli
::: mdscribe.desmond