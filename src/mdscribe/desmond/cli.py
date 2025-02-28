import os
import sys
import re
import copy
import random
import glob
import argparse
import subprocess
import pandas as pd

from datetime import timedelta, datetime
from importlib.resources import files
script_path = files('mdscribe.desmond.script')

from mdscribe import desmond

SCHRODINGER = os.environ["SCHRODINGER"]
USER=os.environ["USER"]

schrodinger_hosts = "{}/schrodinger.hosts".format(SCHRODINGER)
schrodinger_run = "{}/run".format(SCHRODINGER)
multisim = "{}/utilities/multisim".format(SCHRODINGER)



def batch_report() -> None:
    subprocess.run([schrodinger_run, script_path.joinpath('batch-desmond-report.py')])


def batch_dihedral() -> None:
    subprocess.run([schrodinger_run, script_path.joinpath('batch-desmond-dihedral.py')])


def batch_distance() -> None:
    subprocess.run([schrodinger_run, script_path.joinpath('batch-desmond-distance.py')])


def batch_ligrmsd() -> None:
    subprocess.run([schrodinger_run, script_path.joinpath('batch-desmond-ligrmsd.py')])


def batch_rg() -> None:
    subprocess.run([schrodinger_run, script_path.joinpath('batch-desmond-rg.py')])


def batch_extend() -> None:
    # *-in.cms
    # read .cfg file and get last_time
    # /home/shbae/local/schrodinger2020-2/desmond \
    # -JOBNAME ${j} -HOST localhost -gpu -restore ${j}.cpt -in ${j}-in.cms -cfg mdsim.last_time=500000 -WAIT
    raise NotImplementedError


def mdinfo():
    simulation_time = re.compile(r'last_time = "(?P<t>[.0-9]+)"')
    #    last_time = "100000.0"
    progress = re.compile(r'Chemical time:\s+(?P<finished>[.0-9]+) ps, Step:\s+[0-9]+, ns/day:\s+(?P<rate>[.0-9]+)')
    #Chemical time:         33507.6000 ps, Step: 5584600, ns/day:      296.231

    parser = argparse.ArgumentParser(description="desmond MD info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tmpdir', dest='tmpdir', default=None, help='desmond temporary directory')
    parser.add_argument('-l', '--logfile', dest='logfile', default=None)
    args = parser.parse_args()

    if args.tmpdir is None:
        with open(schrodinger_hosts, "r") as f:
            for line in f:
                if line.strip().startswith("#"):
                    continue
                c = line.strip().split(":")
                if len(c) != 2: 
                    continue
                k = c[0].strip()
                v = c[1].strip()
                if k == 'tmpdir':
                    args.tmpdir = os.path.join(v, USER)
            if not os.path.exists(args.tmpdir):
                sys.stderr.write("Use -t or --tmpdir to give temporary directory\n")
                sys.exit(1)


    if args.logfile is None:
        args.logfile = args.tmpdir+"/desmond_*/*.log"
        logfile = None
        for f in glob.glob(args.logfile):
            if not "multisim" in f:
                logfile = f
                break
        if not logfile:
            sys.stderr.write(f"desmond log not found in {args.logfile}\n")
            sys.stderr.write("Use -l or --logfile to give log file path\n")
            sys.stderr.write("Use -t or --tmpdir to give temporary directory\n")
            sys.exit(2)


    total = None
    rate = []
    with open(logfile,"r") as f:
        for line in f:
            m = simulation_time.search(line)
            n = progress.search(line)
            if not total and m and m.group("t"):
                total = float(m.group("t"))*0.001
            if n and n.group("finished") and n.group("rate"):
                rate.append(float(n.group("rate")))
                finished = float(n.group("finished"))*0.001

        n = len(rate)
        print("-"*80)
        print(f"| SCHRODINGER= {SCHRODINGER}")
        print(f"| tmpdir: {args.tmpdir}")
        print(f"|")
        print(f"| Desmond MD Timing Info")
        print(f"| ----------------------")
        print(f"|")
        print(f"| {os.path.dirname(logfile)}/")
        print(f"| {os.path.basename(logfile)}")
        print(f"|")

        if n == 0:
            if total:
                print(f"| Total     {total:9.2f} ns")
            print(f"| Timing data not available yet")
            print("-"*80)
            sys.exit(0)

        remaining = total-finished
        avg_rate = sum(rate[-n:])/n
        eta = 24.0*remaining/avg_rate # hours
        eta_time = datetime.now() + timedelta(hours=eta)
        print(f"| Total     {total:9.2f} ns")
        print(f"| Completed {finished:9.2f} ns ({100.0*finished/total:5.1f}%)")
        print(f"| Remaining {remaining:9.2f} ns")
        print(f"| ")
        print(f"| Average timings")
        print(f"|    ns/day = {avg_rate:9.2f}")
        print(f"| ")
        print(f"| Estimated time")
        print(f"|    remaining = {eta:.2f} hours")
        print(f"|    completed = {eta_time.ctime()}")
        print("-"*80)


def batch_setup():
    setup_msj = desmond.Multisim(template="desmond-setup.msj")
    
    parser = argparse.ArgumentParser(description="batch gdesmond md system setup",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--conc', dest="conc", type=float, default=0.15, help="salt concentration in M")
    parser.add_argument('-d','--dist', dest="dist", type=float, default=10.0, help="buffer distance in A")
    parser.add_argument('-l','--lipid',  dest="lipid", default="", help="lipid bilayer")
    parser.add_argument('--cpp', dest="cpp", default=False, action="store_true", help="CPP simulation setup")
    parser.add_argument('-s','--solvent',  dest="solvent", default="TIP3P", help="solvent model")
    parser.add_argument('-i','--counterion',  dest="counterion", default="Na", help="neutralizing ion")
    parser.add_argument('-n','--negative', dest="neg", default="Cl", help="negative salt ion")
    parser.add_argument('-p','--positive', dest="pos", default="Na", help="positive salt ion")
    parser.add_argument('-f','--forcefield', dest="forcefield", default="S-OPLS", help="forcefield")
    parser.add_argument('-j','--jobfile', dest="job_file", default="desmond_setup_1.sh", help="job filename")
    parser.add_argument('-m','--msjfile', dest="msj_file", default="desmond_setup_1.msj", help="msj filename")
    parser.add_argument('-a','--appendix', dest="appendix", default="", help="job name appendix")
    parser.add_argument('mae', nargs="+", help="desmond mae file")
    args = parser.parse_args()

    if args.appendix:
        msj_file = args.msj_file[:-4] + "_" + args.appendix + ".msj"
    else:
        msj_file = args.msj_file

    if args.appendix:
        job_file = args.job_file[:-3] + "_" + args.appendix + ".sh"
    else:
        job_file = args.job_file

    # job file (.sh) and msj file (.msj) should match
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"

    while os.path.exists(msj_file):
        splited = msj_file.replace(".msj","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        msj_file = "_".join(splited) + ".msj"

    with open(msj_file, "w") as msj:
        setup_msj.dot.build_geometry.solvent = str(args.solvent)
        setup_msj.dot.build_geometry.add_counterion.ion = str(args.counterion)
        setup_msj.dot.build_geometry.salt.concentration = str(args.conc)
        setup_msj.dot.build_geometry.box.size = f"[ {args.dist} {args.dist} {args.dist} ]"
        setup_msj.dot.build_geometry.salt.negative_ion = str(args.neg)
        setup_msj.dot.build_geometry.salt.positive_ion = str(args.pos)
        if args.cpp:
            setup_msj.dot.build_geometry.rezero_system = "False"
            args.lipid = "POPC"
        if args.lipid:
            setup_msj.dot.build_geometry.membrane_box.lipid = str(args.lipid)
        else:
            # remove membrane_box block
            setup_msj.dot.build_geometry.pop("membrane_box")
        setup_msj.dot.build_geometry.override_forcefield = str(args.forcefield)
        setup_msj.dot.assign_forcefield.forcefield = str(args.forcefield)
        setup_msj.dot.assign_forcefield.water = str(args.solvent)
        setup_msj.write(msj)


    with open("README","a") as readme, open(job_file,"w") as job:
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"Force Field   = {args.forcefield}\n")
        readme.write(f"Solvent       = {args.solvent}\n")
        readme.write(f"Counter Ion   = {args.counterion}\n")
        readme.write(f"Positive Ion  = {args.pos}\n")
        readme.write(f"Negative Ion  = {args.neg}\n")
        readme.write(f"Concentration = {args.conc} (M)\n")
        readme.write(f"Size          = {args.dist} (A)\n")
        readme.write(f"Lipid         = {args.lipid}\n")
        readme.write(f"msjfile       = {msj_file}\n")
        readme.write(f"Jobfile       = {job_file}\n")
        readme.write(f"Input structure(s):\n")

        for i, infile in enumerate(args.mae):
            prefix = os.path.basename(infile).split(".")[0]
            job_name = f"desmond_setup-{prefix}"
            if args.appendix:
                job_name += f"-{args.appendix}"
            cms_file = f"{job_name}-out.cms"
            job.write(f"if [ ! -f {cms_file} ]\n")
            job.write(f"then\n")
            job.write(f"{multisim} \\\n")
            job.write(f"  -JOBNAME {job_name} \\\n")
            job.write(f"  -m {msj_file} {os.path.abspath(infile)} \\\n") 
            job.write(f"  -o {cms_file} \\\n")
            job.write(f"  -HOST localhost:20 -maxjob 20 -WAIT\n")
            job.write(f"fi\n")
            readme.write(f"[{i+1:02d}] {infile}\n")
        readme.write("\n\n")

    os.chmod(job_file, 0o777)


def batch_md():
    md_msj = desmond.Multisim(template="desmond-md.msj")
    md_cfg = desmond.Multisim(template="desmond-md.cfg")

    parser = argparse.ArgumentParser(description="batch gdesmond md jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, 
                        metavar="gpu_device", help="gpu device id")
    parser.add_argument('-T', dest="temp", nargs="+", type=float, default=[300.0,], 
                        metavar="temperature", help="temperature in K")
    parser.add_argument('-R', dest="ramp", nargs="+", type=float, default=[10.0,],
                        metavar="tempramp", help="heat and cool ramps in ps/K")
    parser.add_argument('-t', dest="simulation_time", nargs="+", type=float, default=[100.0,],
                        metavar="simulation_time", help="simulation time in ns")
    parser.add_argument('-i', dest="interval", type=float, default=100.0,
                        metavar="interval", help="frame interval in ps")
    parser.add_argument('--posres', dest="posres", default="res. UNK",
                        metavar="posres", help="ASL for positional restraints during prod. sim.")
    parser.add_argument('--posres-force', dest="posres_force", type=float, default=0.0,
                        metavar="posres_force", help="positional restraints force constant (kcal/mol/A**2)")
    parser.add_argument('-p', dest="prefix", default="r", help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, help="number of repeats")
    parser.add_argument('-j', dest="job_file", default="desmond_md_job_1.sh", help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)

    opt  = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella '
    opt += '-set stage[1].set_family.md.jlaunch_opt=["-gpu"] -lic "DESMOND_GPGPU:16"'

    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"

    if len(args.temp) > 1:
        try:
            assert len(args.temp) == len(args.simulation_time)
            assert len(args.temp) == (len(args.ramp) + 1)
        except:
            print("For a variable temperature simulaton, the number of temperatures and simulation times ")
            print("should match and temperature ramp(s) (ps/K) should be given between temperatures.")
            print("Please check -T, -t and -R options")
            sys.exit(0)
    else:
        args.ramp = [] # not used

    # Note: if not specified, 
    # all times are in the unit of ps and 
    # energy is in the unit of kcal/mol.

    with open("README","a") as readme, open(job_file,"w") as job:
        print(f"Job file = {job_file}")

        dirs = [ f"{args.prefix}{n:02d}" for n in range(args.start, args.start+args.repeat) ]

        t_schedule = [ [ args.temp[0], 0.0 ], ]

        if len(args.temp) > 1 :

            elapsed = 0
            prev_temp = args.temp[0]
            idx = 0

            for (temp,simulation_time) in zip(args.temp, args.simulation_time):


                deltaT = abs(temp - prev_temp)

                if deltaT > 0.001: # ramp
                    elapsed += args.ramp[idx] * deltaT # ramp: ps/K
                    t_schedule.append([temp, elapsed])
                    idx += 1

                elapsed += simulation_time * 1000. # ns -> ps
                t_schedule.append([temp, elapsed]) # ns -> ps
                prev_temp = temp
            total_simulation_time = elapsed # ps
        else:
            total_simulation_time = sum(args.simulation_time) * 1000.0 # ps

        if not args.interval:
            # args.simulation_time in ns and args.interval in ps
            # default: make 1000 frames
            args.interval = total_simulation_time / 1000.
        
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"GPU device              = {args.gpu_device}\n")
        readme.write(f"Temperature (K)         = {args.temp}\n")
        readme.write(f"Temperature Ramp (ps/K) = {args.ramp}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Temperature schedule    = {str(t_schedule)}\n")
        readme.write(f"Total Sim. Time (ns)    = {total_simulation_time/1000.0}\n")
        readme.write(f"Trajectory interval(ps) = {args.interval}\n")
        readme.write(f"Repeat                  = {args.repeat}\n")
        readme.write( "Directory               = %s\n" % " ".join(dirs))
        readme.write(f"Jobfile                 = {job_file}\n\n")
        
        for i, infile in enumerate(cms_files):
            info = f"[{i+1:02d}] {infile}"
            print(info)
            readme.write(info+"\n")
        readme.write("\n\n")

        job.write(f'export CUDA_VISIBLE_DEVICES="{args.gpu_device}"\n\n')

        for n in range(args.start, args.start+args.repeat): 
            outdir = f"{args.prefix}{n:02d}"
            outdir_abspath = os.path.abspath(outdir)
            job.write(f"cd {outdir_abspath}/\n\n")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            msj_file = f"{outdir}/desmond_md_job_{n:02d}.msj"
            cfg_file = f"{outdir}/desmond_md_job_{n:02d}.cfg"
            cfg_file_basename= os.path.basename(cfg_file)
            
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

            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                
                job_name = f"desmond_md_job_{n:02d}_{prefix}"
                job.write('if [ ! -f {}/{} ]\n'.format( outdir_abspath, f"{job_name}-out.cms",))
                job.write('then\n')
                job.write('{} -JOBNAME {} -m {} -c {} -description "{}" {} {} -o {} -WAIT\n'.format(
                    multisim, 
                    job_name, 
                    os.path.basename(msj_file),
                    os.path.basename(cfg_file),
                    "GPU desmond MD",
                    opt,
                    os.path.join("..",infile),
                    f"{job_name}-out.cms",
                ))
                job.write('fi\n\n')

    os.chmod(job_file, 0o777)


def batch_min():
    min_msj = desmond.Multisim(template="desmond-min.msj")
    min_cfg = desmond.Multisim(template="desmond-min.cfg")

    opt = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella '
    opt += '-lic "DESMOND_GPGPU:16" '
    opt += '-description "minimization"'

    parser = argparse.ArgumentParser(description="batch gdesmond minimization jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, 
                        help="gpu device id")
    parser.add_argument('-t', dest="simulation_time", type=float, default=100.0, 
                        help="simulation time in ps")
    parser.add_argument('-p', dest="prefix", type=str, default="r", 
                        help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, 
                        help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, 
                        help="number of repeats")
    parser.add_argument('-j', dest="job_file", type=str, default="desmond_min_job_1.sh", 
                        help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)


    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"


    with open("README","a") as readme, open(job_file,"w") as job:
        print("\n" + job_file + "\n")
        outdir_nums = list(range(args.start, args.start+args.repeat))
        outdirs = [f"{args.prefix}{num:02d}" for num in outdir_nums]
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"GPU device              = {args.gpu_device}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Repeat                  = {args.repeat}\n")
        readme.write( "Directory               = %s\n" % " ".join(outdirs))
        readme.write(f"Jobfile                 = {job_file}\n\n")

        job.write(f'export CUDA_VISIBLE_DEVICES="{args.gpu_device}"\n\n')

        for i, infile in enumerate(cms_files):
            info = f"[{i+1}] {infile}"
            print(info)
            readme.write(info+"\n")
        print()
        readme.write("\n")

        for (num, outdir) in zip(outdir_nums, outdirs):
            outdir_abspath = os.path.abspath(outdir)
            job.write(f"cd {outdir_abspath}/\n\n")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            cfg_file = f"{outdir}/desmond_min_job_{num:02d}.cfg"
            msj_file = f"{outdir}/desmond_min_job_{num:02d}.msj"
            cfg_file_basename= os.path.basename(cfg_file)

            with open(cfg_file,"w") as cfg, open(msj_file,"w") as msj:
                min_cfg.dot.randomize_velocity.seed = str(random.randint(1000,9999))
                min_cfg.dot.time = str(args.simulation_time) # (ps)
                min_cfg.write(cfg)

                min_msj.dot.simulate.cfg_file = cfg_file_basename
                min_msj.write(msj)
            
            # append to job file
            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                job_name = f"desmond_min_job_{num:02d}_{prefix}"
                job.write('{} -JOBNAME {} -m {} -c {} {} {} -o {} -WAIT\n\n'.format(
                    multisim, 
                    job_name,
                    os.path.basename(msj_file),
                    os.path.basename(cfg_file),
                    opt,
                    infile,
                    f"{job_name}-out.cms",
                ))

    os.chmod(job_file, 0o777)


def index_or_asl(expr:str) -> str:
    if len(expr.split()) == 1:
        return expr
    else:
        return f'"{expr}"'
        

def batch_metad():
    """Metadynamics Simulations.
    """

    metad_msj = desmond.Multisim(template="desmond-metad.msj")
    metad_cfg = desmond.Multisim(template="desmond-metad.cfg")

    opt = '-HOST localhost -maxjob 1 -cpu 1 -mode umbrella -lic "DESMOND_GPGPU:16" '
    opt += '-description "metadynamics"'

    parser = argparse.ArgumentParser(description="batch gdesmond metadynamics jobs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--meta-height', dest="meta_height", type=float, default=0.03, help="Height of the repulsive Gaussian potential in kcal/mol (0.03)")
    parser.add_argument('--meta-interval', dest="meta_interval", type=float, default=0.09, help="Interval in ps at which Gaussian potentials are added (0.09)")
    parser.add_argument('--meta-first', dest="meta_first", type=float, default=0.0, help="Time in ps at which the Gaussian potentials are first added (0)")
    parser.add_argument('--meta-last', dest="meta_last", type=float, default=-1.0, help="Time in ps at which the Gaussian potentials are last added (simulation time)")
    parser.add_argument('--meta-kTemp', dest="meta_kTemp", type=float, default=2.4, help="Perform the well-tempered metadynamics (2.4)")
    parser.add_argument('--cv-dist', dest='cv_dist', nargs=4, action='append', default=[], help="atom1 atom2 width(0.05) wall")
    parser.add_argument('--cv-angle', dest='cv_angle', nargs=4, action='append', default=[], help="atom1 atom2 atom3 width(2.5)")
    parser.add_argument('--cv-dihedral', dest='cv_dihedral', nargs=5, action='append', default=[], help="atom1 atom2 atom3 atom4 width(5.0)")
    parser.add_argument('--cv-rgyr', dest='cv_rgyr', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rgyr-mass', dest='cv_rgyr_mass', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rmsd', dest='cv_rmsd', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-rmsd-symm', dest='cv_rmsd_symm', nargs=2, action='append', default=[], help="atom width(0.1)")
    parser.add_argument('--cv-zdist', dest='cv_zdist', nargs=3, action='append', default=[], help="atom width wall(0.05)")
    parser.add_argument('--cv-zdist0', dest='cv_zdist0', nargs=3, action='append', default=[], help="atom width wall(0.1)")
    parser.add_argument('-g', dest="gpu_device", type=int, default=0, help="gpu device id")
    parser.add_argument('-T', dest="temperature", type=float, default=300.0, help="temperature in Kelvin")
    parser.add_argument('-t', dest="simulation_time", type=float, default=40.0, help="simulation time in ns")
    parser.add_argument('-i', dest="interval", type=float, default=40.0, help="frame interval in ps")
    parser.add_argument('-p', dest="prefix", type=str, default="m", help="directory prefix")
    parser.add_argument('-s', dest="start", type=int, default=1, help="directory start")
    parser.add_argument('-r', dest="repeat", type=int, default=1, help="number of repeats")
    parser.add_argument('-j', dest="job_file", type=str, default="desmond_metadynamics_job_1.sh", help="job filename")
    parser.add_argument('cms', nargs="+", help="desmond cms file")
    args = parser.parse_args()

    try:
        cms_files = [os.path.abspath(f) for f in args.cms]
        assert(len(cms_files) > 0)
    except:
        print(".cms file(s) not found")
        sys.exit(0)

    if not args.interval:
        args.interval = args.simulation_time # it will save 1000 frames

    job_file = args.job_file
    while os.path.exists(job_file):
        splited = job_file.replace(".sh","").split("_")
        splited[-1] = str(int(splited[-1]) + 1)
        job_file = "_".join(splited) + ".sh"


    with open("README","a") as readme, open(job_file,"w") as job:
        print("\n" + job_file + "\n")
        outdir_nums = list(range(args.start, args.start+args.repeat))
        outdirs = [f"{args.prefix}{num:02d}" for num in outdir_nums]
        cmd_echo = ""
        for argv in sys.argv:
            if cmd_echo:
                cmd_echo += " "
            if " " in argv:
                cmd_echo += f'"{argv}"'
            else:
                cmd_echo += f'{argv}'
        readme.write(f"{cmd_echo}\n\n")
        readme.write(f"GPU device              = {args.gpu_device}\n")
        readme.write(f"Temperature (K)         = {args.temperature}\n")
        readme.write(f"Simulation Time (ns)    = {args.simulation_time}\n")
        readme.write(f"Trajectory interval(ps) = {args.interval}\n")
        readme.write(f"Repeat                  = {args.repeat}\n")
        readme.write( "Directory               = %s\n" % " ".join(outdirs))
        readme.write(f"Jobfile                 = {job_file}\n\n")

        job.write(f'export CUDA_VISIBLE_DEVICES="{args.gpu_device}"\n\n')

        for i, infile in enumerate(cms_files):
            info = f"[{i+1}] {infile}"
            print(info)
            readme.write(info+"\n")
        print()
        readme.write("\n")

        for (num, outdir) in zip(outdir_nums, outdirs):
            outdir_abspath = os.path.abspath(outdir)
            job.write(f"cd {outdir_abspath}/\n\n")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            cfg_file = f"{outdir}/desmond_metadynamics_job_{num:02d}.cfg"
            msj_file = f"{outdir}/desmond_metadynamics_job_{num:02d}.msj"
            cfg_file_basename= os.path.basename(cfg_file)

            with open(cfg_file,"w") as cfg, open(msj_file,"w") as msj:
                metad_cfg.dot.randomize_velocity.seed = str(random.randint(1000,9999))
                metad_cfg.dot.time = str(args.simulation_time*1000.0)
                metad_cfg.dot.trajectory.interval = str(args.interval)
                metad_cfg.dot.temperature = f"[ [{args.temperature} 0] ]"
                metad_cfg.write(cfg)

                metad_msj.dot.simulate[-1].cfg_file = cfg_file_basename
                
                # metadynamics
                cv_list = []


                for atom1, atom2, width, wall in args.cv_dist:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    cv.atom = f'[{atom1} {atom2}]'
                    cv.type = 'dist'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)
                
                for atom1, atom2, atom3, width in args.cv_angle:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    atom3 = index_or_asl(atom3)
                    cv.atom = f'[{atom1} {atom2} {atom3}]'
                    cv.type = 'angle'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, atom2, atom3, atom4, width in args.cv_dihedral:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    atom2 = index_or_asl(atom2)
                    atom3 = index_or_asl(atom3)
                    atom4 = index_or_asl(atom4)
                    cv.atom = f'[{atom1} {atom2} {atom3} {atom4}]'
                    cv.type = 'dihedral'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rmsd:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rmsd'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rmsd_symm:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rmsd_symm'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rgyr:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rgyr'
                    cv.width = str(width)
                    cv_list.append(cv)
                
                for atom1, width in args.cv_rgyr_mass:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'rgyr_mass'
                    cv.width = str(width)
                    cv_list.append(cv)

                for atom1, width, wall in args.cv_zdist:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'zdist'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)
                
                for atom1, width, wall in args.cv_zdist0:
                    cv = copy.deepcopy(metad_msj.dot.simulate[-1].meta.cv[0])
                    atom1 = index_or_asl(atom1)
                    cv.atom = f'[{atom1}]'
                    cv.type = 'zdist0'
                    cv.width = str(width)
                    if not wall.startswith('-'):
                        cv.wall = str(wall)
                    cv_list.append(cv)

                metad_msj.dot.simulate[-1].meta.cv = cv_list
                metad_msj.dot.simulate[-1].meta.height = str(args.meta_height)
                metad_msj.dot.simulate[-1].meta.interval = str(args.meta_interval)
                metad_msj.dot.simulate[-1].meta.first = str(args.meta_first)

                if args.meta_last > 0.0:
                    metad_msj.dot.simulate[-1].meta.last = str(args.meta_last)
                
                # well-tempered metadynamics
                if args.meta_kTemp > 0.0:
                    metad_msj.dot.simulate[-1].meta.kTemp = str(args.meta_kTemp)
                
                metad_msj.write(msj)
            
            # append to job file
            for infile in cms_files:
                prefix_ = re.sub(r'desmond_setup[-_]','', os.path.basename(infile))
                prefix  = prefix_.replace("-out.cms","")
                job_name = f"desmond_metadynamics_job_{num:02d}_{prefix}"
                job.write('{} -JOBNAME {} -m {} -c {} {} {} -o {} -WAIT\n\n'.format(
                    multisim, 
                    job_name,
                    os.path.basename(msj_file),
                    os.path.basename(cfg_file),
                    opt,
                    infile,
                    f"{job_name}-out.cms",
                ))

    os.chmod(job_file, 0o777)


def read_eaf(filename, verbose=False):
    HBond = {}
    Hydrophobic = {}
    WaterBridge = {}
    Polar = {}
    HalogenBond = {}
    LigWat = {}
    Metal = {}
    PiCat = {}
    PiPi = {}

    result = desmond.expr.parse_file(filename)
    d = result.as_dict()
    desmond.traverse_dict(d)
    dot = desmond.DotMap(d)

    for section in dot.Keywords:
        try:
            assert section.ProtLigInter.HBondResult
            num_frames = len(section.ProtLigInter.HBondResult)
            for frame in section.ProtLigInter.HBondResult:
                # [[3 "_:ARG_143:HH22" d-s "L-FRAG_0:N6" ]]
                for (frameno, prot, hbond_type, lig) in frame:
                    prot = prot.strip('\"')
                    (_, resid, atom) = prot.split(":")
                    (resName, resSeq) = resid.split("_")
                    resSeq = int(resSeq)
                    if resSeq in HBond:
                        HBond[resSeq]['count'] += 1
                    else:
                        HBond[resSeq] = {'resName': resName, 'count':1 }
            for resSeq in sorted(HBond):
                fraction = float(HBond[resSeq]['count'])/num_frames
                if verbose:
                    print(f"HBond {HBond[resSeq]['resName']}_{resSeq} {fraction:5.3f} {num_frames}")
        except:
            pass

        try:
            assert section.ProtLigInter.HydrophobicResult
            num_frames = len(section.ProtLigInter.HydrophobicResult)
            for frame in section.ProtLigInter.HydrophobicResult:
                # [[0 "_:PHE_223" L-FRAG_0 ] [0 "_:ALA_241" L-FRAG_0 ]]
                for (frameno, prot, lig) in frame:
                    prot = prot.strip('\"')
                    (_, resid) = prot.split(":")
                    (resName, resSeq) = resid.split("_")
                    resSeq = int(resSeq)
                    if resSeq in Hydrophobic:
                        Hydrophobic[resSeq]['count'] += 1
                    else:
                        Hydrophobic[resSeq] = {'resName': resName, 'count':1 }
            for resSeq in sorted(Hydrophobic):
                fraction = float(Hydrophobic[resSeq]['count'])/num_frames
                if verbose:
                    print(f"Hydrophobic {Hydrophobic[resSeq]['resName']}_{resSeq} {fraction:5.3f} {num_frames}")
        except:
            pass
        
        try:
            assert section.ProtLigInter.PolarResult
            num_frames = len(section.ProtLigInter.PolarResult)
            for frame in section.ProtLigInter.PolarResult:
                # [[1 "_:GLU_216:OE2" b "L-FRAG_1:N3" 4.45 ]]
                for (frameno, prot, _, lig, _) in frame:
                    prot = prot.strip('\"')
                    (_, resid, atom) = prot.split(":")
                    (resName, resSeq) = resid.split("_")
                    resSeq = int(resSeq)
                    if resSeq in Polar:
                        Polar[resSeq]['count'] += 1
                    else:
                        Polar[resSeq] = {'resName': resName, 'count':1 }
            for resSeq in sorted(Polar):
                fraction = float(Polar[resSeq]['count'])/num_frames
                if verbose:
                    print(f"Polar {Polar[resSeq]['resName']}_{resSeq} {fraction:5.3f} {num_frames}")
        except:
            pass

        try:
            assert section.ProtLigInter.WaterBridgeResult
            num_frames = len(section.ProtLigInter.WaterBridgeResult)
            for frame in section.ProtLigInter.WaterBridgeResult:
                # [[3 "_:GLU_216:OE2" a "L-FRAG_0:N2" a 2431 ]]
                for (frameno, prot, _, lig, _, _) in frame:
                    prot = prot.strip('\"')
                    (_, resid, atom) = prot.split(":")
                    (resName, resSeq) = resid.split("_")
                    resSeq = int(resSeq)
                    if resSeq in WaterBridge:
                        WaterBridge[resSeq]['count'] += 1
                    else:
                        WaterBridge[resSeq] = {'resName': resName, 'count':1 }
            for resSeq in sorted(WaterBridge):
                fraction = float(WaterBridge[resSeq]['count'])/num_frames
                if verbose:
                    print(f"WaterBridge {WaterBridge[resSeq]['resName']}_{resSeq} {fraction:5.3f} {num_frames}")
        except:
            pass

    return num_frames, HBond, Hydrophobic, Polar, WaterBridge, HalogenBond, LigWat, Metal, PiCat, PiPi


def batch_pli():
    teststring ="""Keywords = [
    {RMSD = {
        ASL = "((mol. 1 and backbone) and not (atom.ele H) and not (mol. 2))"
        Frame = 0
        Panel = pl_interact_survey
        Result = [0.0 1.161 1.286 1.331 1.176 1.195 ]
        SelectionType = Backbone
        Tab = pl_rmsd_tab
        Type = ASL
        Unit = Angstrom
        }
    }
    {RMSD = {
        ASL = "(mol. 1 and sidechain and not (mol. 2))"
        Frame = 0
        Panel = pl_interact_survey
        Result = [0.0 1.161 1.286 1.331 1.176 1.195 ]
        SelectionType = "Side chains"
        Tab = pl_rmsd_tab
        Type = ASL
        Unit = Angstrom
        }
    }
    ]"""
    # result = expr.parse_string(teststring)
    # d = result.as_dict()
    # traverse_dict(d)
    # dot = DotMap(d)
    # try:
    #     assert dot.Keywords[1].RMSD.SelectionType == '"Side chains"'
    #     print("ok")
    # except:
    #     print("error")

    parser = argparse.ArgumentParser(description="Average Protein-Ligand Interactions",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--out', dest='out', default='mean-PLIF', help="output basename")
    parser.add_argument('eaf', nargs='+', default=[], help='input -out.eaf filename(s)')
    args = parser.parse_args()

    if len(args.eaf) == 0:
        argparse.print_help()
        sys.exit(0)

    total_num_frames = 0
    data = {}
    for filename in args.eaf:
        (num_frames, HBond, Hydrophobic, Polar, WaterBridge, 
        HalogenBond, LigWat, Metal, PiCat, PiPi) = read_eaf(filename, verbose=False)
        
        total_num_frames += num_frames
        
        print(f"{filename}  {num_frames} frames")
        
        for resSeq in HBond:
            if resSeq in data:
                data[resSeq]['hbond'] += HBond[resSeq]['count']
            else:
                data[resSeq] = {'resName': HBond[resSeq]['resName'], 
                                'hbond': HBond[resSeq]['count'],
                                'hydrophobic': 0,
                                'polar': 0,
                                'waterbridge': 0,
                                }
        for resSeq in Hydrophobic:
            if resSeq in data:
                data[resSeq]['hydrophobic'] += Hydrophobic[resSeq]['count']
            else:
                data[resSeq] = {'resName': Hydrophobic[resSeq]['resName'], 
                                'hbond': 0,
                                'hydrophobic': Hydrophobic[resSeq]['count'],
                                'polar': 0,
                                'waterbridge': 0,
                                }
        for resSeq in Polar:
            if resSeq in data:
                data[resSeq]['polar'] += Polar[resSeq]['count']
            else:
                data[resSeq] = {'resName': Polar[resSeq]['resName'], 
                                'hbond': 0,
                                'hydrophobic': 0,
                                'polar': Polar[resSeq]['count'],
                                'waterbridge': 0,
                                }
        for resSeq in WaterBridge:
            if resSeq in data:
                data[resSeq]['waterbridge'] += WaterBridge[resSeq]['count']
            else:
                data[resSeq] = {'resName': WaterBridge[resSeq]['resName'],
                                'hbond': 0,
                                'hydrophobic' : 0,
                                'polar': 0,
                                'waterbridge': WaterBridge[resSeq]['count'],
                                } 

    csvdata = {'resid':[], 
            'resSeq':[], 
            'resName':[], 
            'hbond':[], 
            'hydrophobic':[],
            'polar': [],
            'waterbridge': [],
            }
    
    for resSeq in sorted(data):
        csvdata['resSeq'].append(resSeq)
        csvdata['resName'].append(data[resSeq]['resName'])
        csvdata['resid'].append(f"{data[resSeq]['resName']}_{resSeq}")
        csvdata['hbond'].append(float(data[resSeq]['hbond'])/total_num_frames)
        csvdata['hydrophobic'].append(float(data[resSeq]['hydrophobic'])/total_num_frames)
        csvdata['polar'].append(float(data[resSeq]['polar'])/total_num_frames)
        csvdata['waterbridge'].append(float(data[resSeq]['waterbridge'])/total_num_frames)

    df = pd.DataFrame(csvdata)
    df.to_csv(args.out + '.csv', index=False, float_format='%.4f')
    g = df.loc[:, ~df.columns.isin(['resSeq','resName'])].plot.bar(
        x="resid", 
        stacked=True,
        title="Protein-Ligand Interactions", 
        xlabel='Residue', 
        ylabel='Fraction of MD trajectory', 
        figsize=(8,3), 
        fontsize=8,
        )
    fig = g.get_figure()
    fig.savefig(args.out + '.pdf', bbox_inches="tight", pad_inches=0.2, dpi=150)