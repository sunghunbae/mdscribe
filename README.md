Python Utilities for Molecular Dynamics
=======================================

Introduction
============

MDScribe is a set of command-line scripts and a library written to make 
setting up the molecular dynamics simulations easier for Desmond, OpenMM, and Amber.

Install
===============

```shell
$ pip install mdscribe
```

Usage for Desmond
=================

```python
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

    md_msj.write(msj)

with open(cfg_file,"w") as cfg:
    # modify desmond cfg template
    md_cfg.dot.randomize_velocity.seed = random.randint(1000, 9999)
    md_cfg.dot.time = total_simulation_time
    md_cfg.dot.temperature = t_schedule
    md_cfg.dot.trajectory.interval = args.interval
    md_cfg.write(cfg)
```