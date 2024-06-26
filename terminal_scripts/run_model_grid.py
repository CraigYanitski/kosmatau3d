# python3.8

import argparse
import numpy as np
import os

from kosmatau3d import models

# from pandas.compat.pyarrow import pa
from pprint import pprint


# ARGUMENTS
# , fontsize=16 Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "-f",
    "--folder",
    type=str,
    default="/mnt/Milky-Way_models/yanitski/projects/pdr/KT3_history", #"/mnt/yanitski_backup/yanitski/projects/pdr/KT3_history",
    help="Folder containing kosmatau3d models",
)
parser.add_argument(
    "-g",
    "--grid",
    type=str,
    default="convergence",
    choices=[
        "convergence",
        "f_cm-cm",
        "f_icm-icm",
        "cm-icm",
        "cm-f_fuv",
        "f_cm-f_icm",
        "r_cmz-f_fuv",
        "f_hi-f_fuv",
        "f_cl-f_icl-f_n-f_fuv",
        "fuv_cl",
        "fuv_icl",
        "const_fuv",
        "vox_disp",
        "f_fuv_gc",
        "disp_gc",
        "fuv_gc-disp_gc",
        "fuv_gc-mass_gc",
        "wnm",
    ],
    help="Parameters varied in grid",
)
parser.add_argument(
    "-r", "--resolution", type=int, default=400, help="Voxel size in the model"
)
parser.add_argument(
    "-m",
    "--mp",
    type=int,
    default=8,
    help="Number of threads to open during RT \n  set to 0 for debugging",
)
parser.add_argument(
    "-e",
    "--exec",
    type=str,
    default="all",
    choices=["all", "model", "RT", "HI"],
    help="Which portion of the code to evaluate",
)
parser.add_argument(
    "-o",
    "--overwrite",
    type=str,
    default="false",
    choices=["true", "false"],
    help="Overwrite the pre-computed model",
)
parser.add_argument(
    "--test",
    type=str,
    default="false",
    choices=["true", "false"],
    help="test grid with the first model in the 'convergence' grid",
)
parser.add_argument(
    "-i",
    "--index",
    type=int,
    default=1,
    help="Index from which to begin grid. Useful if a grid is interrupted.",
)
parser.add_argument(
    "-d", "--debug", type=bool, default=False, help="Use debugging mode"
)
parser.add_argument("-t", "--timed", type=bool, default=False, help="Use timing mode.")
parser.add_argument(
    "-v", "--verbose", type=bool, default=False, help="Use verbose mode."
)

args = parser.parse_args()

if args.test == "true":
    args.grid = "convergence"
    test = True
else:
    test = False

# convergence
if args.grid == "convergence":
    if test:
        folder = "r{}_convergence_test/"
        res = (400,)
    else:
        folder = "r{}_convergence/"
        res = (400, 200, 100)
    diameter = list(36 - _ / 1e3 for _ in res)
    grid_flag = True
    param_keys = ("resolution", "x", "y", "new_grid")
    params = [res, diameter, diameter]
    param_folders = list(folder.format(*_) for _ in zip(*params))
    params.append([grid_flag] * len(params[0]))
# f_cm-cm
elif args.grid == "f_cm-cm":
    folder = f"r{args.resolution}" + "_fcm{}_cm{}/"
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    param_keys = ("h2_mass_factor", "clump_mass_range", "clump_mass_number")
    fmass1_mesh, i_cl_mesh = list(
        _.flatten() for _ in np.meshgrid(fmass1, np.arange(len(mass_cl)))
    )
    mass_bin_mesh = []
    for _ in range(i_cl_mesh.size):
        mass_bin_mesh.append([mass_cl[i_cl_mesh[_]], [-2]])
    mass_num_mesh = list(
        list(int(np.max(mbin) - np.min(mbin)) + 1 for mbin in _) for _ in mass_bin_mesh
    )
    params = (fmass1_mesh, mass_bin_mesh, mass_num_mesh)
    mass_cl_strs = list("_".join([str(_) for _ in mass[0]]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip(*(params[0], mass_cl_strs)))
# f_icm-icm
elif args.grid == "f_icm-icm":
    folder = f"r{args.resolution}" + "_ficm{}_icm{}/"
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    param_keys = ("hi_mass_factor", "clump_mass_range", "clump_mass_number")
    fmass2_mesh, i_icl_mesh = list(
        _.flatten() for _ in np.meshgrid(fmass2, np.arange(len(mass_icl)))
    )
    mass_bin_mesh = []
    for _ in range(i_icl_mesh.size):
        mass_bin_mesh.append([[0, 2], mass_icl[i_icl_mesh[_]]])
    mass_num_mesh = list(
        list(int(np.max(mbin) - np.min(mbin)) + 1 for mbin in _) for _ in mass_bin_mesh
    )
    params = (fmass2_mesh, mass_bin_mesh, mass_num_mesh)
    mass_icl_strs = list("_".join([str(_) for _ in mass[1]]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip(*(params[0], mass_icl_strs)))
# cm-icm
elif args.grid == "cm-icm":
    folder = f"r{args.resolution}" + "_cm{}_icm{}/"
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    param_keys = ("clump_mass_range", "clump_mass_number")
    i_cl_mesh, i_icl_mesh = list(
        _.flatten()
        for _ in np.meshgrid(np.arange(len(mass_cl)), np.arange(len(mass_icl)))
    )
    mass_bin_mesh = []
    for _ in range(i_cl_mesh.size):
        mass_bin_mesh.append([mass_cl[i_cl_mesh[_]], mass_icl[i_icl_mesh[_]]])
    mass_num_mesh = list(
        list(int(np.max(mbin) - np.min(mbin)) + 1 for mbin in _) for _ in mass_bin_mesh
    )
    params = (mass_bin_mesh, mass_num_mesh)
    mass_cl_strs = list("_".join([str(_) for _ in mass[0]]) for mass in params[0])
    mass_icl_strs = list("_".join([str(_) for _ in mass[1]]) for mass in params[0])
    param_folders = list(folder.format(*_) for _ in zip(*(mass_cl_strs, mass_icl_strs)))
# cm-icm
elif args.grid == "cm-f_fuv":
    folder = f"r{args.resolution}" + "_cm{}_f_fuv{:.2f}/"
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    # mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    f_fuv = 10 ** np.linspace(0, 4, num=9)
    param_keys = ("clump_mass_range", "clump_mass_number", "fuv_factor")
    i_cl_mesh, f_fuv_mesh = list(
        _.flatten() for _ in np.meshgrid(np.arange(len(mass_cl)), f_fuv)
    )
    mass_bin_mesh = []
    for _ in range(i_cl_mesh.size):
        mass_bin_mesh.append([mass_cl[i_cl_mesh[_]], [-3]])
    mass_num_mesh = list(
        list(int(np.max(mbin) - np.min(mbin)) + 1 for mbin in _) for _ in mass_bin_mesh
    )
    params = (mass_bin_mesh, mass_num_mesh, f_fuv_mesh)
    # pprint(params)
    mass_cl_strs = list("_".join([str(_) for _ in mass[0]]) for mass in params[0])
    mass_icl_strs = list("_".join([str(_) for _ in mass[1]]) for mass in params[0])
    param_folders = list(folder.format(*_) for _ in zip(*(mass_cl_strs, f_fuv_mesh)))
    # pprint(param_folders)
# f_cm-f_icm
elif args.grid == "f_cm-f_icm":
    folder = f"r{args.resolution}" + "_f_cm{}_f_icm{:}/"
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    param_keys = ("h2_mass_factor", "hi_mass_factor")
    params = list(_.flatten() for _ in np.meshgrid(fmass1, fmass2))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# r_cmz-f_fuv
elif args.grid == "r_cmz-f_fuv":
    folder = f"r{args.resolution}" + "_rcmz{}_f_fuv{:.1f}/"
    r_cmz = (200, 600, 1000, 1400)  # np.arange(0, 1001, 50, dtype=int)
    f_fuv = 10 ** np.linspace(0, 4, num=9)
    param_keys = ("r_cmz", "fuv_factor")
    params = list(_.flatten() for _ in np.meshgrid(r_cmz, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# f_hi-f_fuv
elif args.grid == "f_hi-f_fuv":
    folder = f"r{args.resolution}" + "_fhi{}_fuv{:.1f}/"
    f_hi = [1.0, 0.8, 0.6, 0.4]  # [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
    f_fuv = 10 ** np.linspace(0, 4, num=9)
    param_keys = ("interclump_hi_ratio", "fuv_factor")
    params = list(_.flatten() for _ in np.meshgrid(f_hi, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# f_cl-f_icl-f_n-f_fuv
elif args.grid == "f_cl-f_icl-f_n-f_fuv":
    folder = f"r{args.resolution}" + "_f_cm{:.1f}_f_icm{:.1f}_f_n{:.0f}_f_fuv{:.0f}/"
    f_cl = 10 ** np.array([-1, 0, 1, 2], dtype=float)
    f_icl = 10 ** np.array([-1, 0, 1, 2], dtype=float)
    f_d = 10 ** np.array([0, 1, 2], dtype=float)
    f_fuv = 10 ** np.array([0, 1, 2], dtype=float)
    param_keys = ("h2_mass_factor", "hi_mass_factor", "density_factor", "fuv_factor")
    params = list(_.flatten() for _ in np.meshgrid(f_cl, f_icl, f_d, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# f_fuv_gc
elif args.grid == "f_fuv_gc":
    folder = f"r{args.resolution}" + "_f_fuv_gc{:.2f}_r_gc{}/"
    fuv_gc = 10 ** np.linspace(0, 4, num=9)
    r_gc = (200, 600, 1000, 1400)
    param_keys = ("scale_gc", "r_gc")
    params = list(_.flatten() for _ in np.meshgrid(fuv_gc, r_gc))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# disp_gc
elif args.grid == "disp_gc":
    folder = f"r{args.resolution}" + "_r_gc{:.2f}_disp_gc{:.2f}/"
    r_gc = (200, 600)
    disp_gc = (200, 300)
    param_keys = (
        "disp_core",
        "r_core",
    )
    params = list(_.flatten() for _ in np.meshgrid(r_gc, disp_gc))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# fuv_gc-disp_gc
elif args.grid == "fuv_gc-disp_gc":
    folder = f"r{args.resolution}" + "_f_fuv_gc{:.2f}_r_fuv_gc{}_disp_gc{}_r_disp_gc{}/"
    r_gc = (200, 600, 1000, 1400)
    # fuv_gc = (10**np.array([0, 0.5, 1, 1.5, 2]), )
    fuv_gc = (*10 ** np.array([0.5, 1, 1.5, 2, 2.5, 3]),)
    r_core = (200,)
    disp_gc = (
        150,
        175,
        200,
    )
    grid_flag = (True,)
    param_keys = ("scale_gc", "r_gc", "disp_core", "r_core", "new_grid")
    params = list(
        _.flatten() for _ in np.meshgrid(fuv_gc, r_gc, disp_gc, r_core, grid_flag)
    )
    param_folders = list(folder.format(*_) for _ in zip(*params))
# fuv_gc-mass_gc
elif args.grid == "fuv_gc-mass_gc":
    folder = (
        f"r{args.resolution}" + "_r_gc{}_f_fuv_gc{:.2f}_f_mhi_gc{:.2f}_f_mh2_gc{:.2f}/"
    )
    r_gc = (200, 600, 1000, 1400)
    # fuv_gc = (10**np.array([0, 0.5, 1, 1.5, 2]), )
    fuv_gc = (*10 ** np.array([0, 0.5, 1, 1.5, 2, 2.5, 3.0]),)
    mhi_gc = (*10 ** np.array([0, 0.5, 1, 1.5]),)
    mh2_gc = (*10 ** np.array([0, 0.5, 1, 1.5]),)
    grid_flag = (True,)
    param_keys = ("r_gc", "scale_gc", "mhi_gc", "mh2_gc", "new_grid")
    params = list(
        _.flatten() for _ in np.meshgrid(r_gc, fuv_gc, mhi_gc, mh2_gc, grid_flag)
    )
    param_folders = list(folder.format(*_) for _ in zip(*params))
# fuv_cl
elif args.grid == "fuv_cl":
    folder = f"r{args.resolution}" + "_log_fuv_cl{:.2f}_f_fuv{:.2f}/"
    log_fuv = np.arange(0.1, 1.0, 0.2)  # 0.1, 2.6, 0.2
    f_fuv = 10 ** np.arange(0.25, 3.1, 0.25)
    param_keys = ("clump_log_fuv", "fuv_factor")
    params = list(_.flatten() for _ in np.meshgrid(log_fuv, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# fuv_icl
elif args.grid == "fuv_icl":
    folder = f"r{args.resolution}" + "_log_fuv_icl{:.2f}/"
    log_fuv = np.arange(1.1, 2.6, 0.1)
    param_keys = ("interclump_log_fuv",)
    params = list(_.flatten() for _ in np.meshgrid(log_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# const_fuv
elif args.grid == "const_fuv":
    folder = f"r{args.resolution}" + "_log_fuv_cl{:.2f}_log_fuv_icl{:.2f}/"
    log_fuv = np.arange(1.1, 2.6, 0.1)
    param_keys = (
        "clump_log_fuv",
        "interclump_log_fuv",
    )
    params = list(_.flatten() for _ in np.meshgrid(log_fuv, log_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# vox_disp
elif args.grid == "vox_disp":
    folder = f"r{args.resolution}" + "_vox_disp{:.2f}/"
    vox_disp = [0.001, 1.1, 1.86, 2.42, 4.47, 5.82, 10.72]
    param_keys = ("disp_gmc",)
    params = list(_.flatten() for _ in np.meshgrid(vox_disp))
    param_folders = list(folder.format(*_) for _ in zip(*params))
# wnm
elif args.grid == "wnm":
    folder = f"r{args.resolution}" + "_f_wnm{:.1f}_log_fuv_wnm{:.1f}/"
    f_wnm = [0.2, 0.3, 0.4, 0.5]
    log_fuv = [2, 3, 4, 5]
    interclump_idx = [False, True, True]
    interclump_wnm_idx = [False, False, True]
    clump_mass_range = [[0, 2], [-3], [-3]]
    clump_mass_number = [3, 1, 1]
    clump_n_max = [1, 100, 100]
    ensemble_mass_factor = [1, 1, 1]
    param_keys = (
        "interclump_wnm_ratio",
        "interclump_wnm_log_fuv",
        "interclump_idx",
        "interclump_wnm_idx",
        "clump_mass_range",
        "clump_mass_number",
        "clump_n_max",
        "ensemble_mass_factor",
    )
    params = list(_.flatten() for _ in np.meshgrid(f_wnm, log_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
    for _ in (
        interclump_idx,
        interclump_wnm_idx,
        clump_mass_range,
        clump_mass_number,
        clump_n_max,
        ensemble_mass_factor,
    ):
        params.append([_] * len(params[0]))
# unknown
else:
    params = None
    param_keys = None
    param_folders = None
    exit("Invalid grid queried; exiting...")

if args.exec == "model":
    run_model = True
    run_rt = False
    hi = False
elif args.exec == "RT":
    run_model = False
    run_rt = True
    hi = False
elif args.exec == "HI":
    print("HI calculation")
    run_model = False
    run_rt = False
    hi = True
elif args.exec == "all":
    run_model = True
    run_rt = True
    hi = False
else:
    run_model = False
    run_rt = False
    hi = False

if args.overwrite == "true":
    overwrite = True
else:
    overwrite = False

index = args.index - 1

debug = args.debug


# KOSMATAU3D PARAMETERS
# Edit these parameters according to the model you want to produce.
parameters = {
    # Data files
    "clump_tau_grid_file": "clump_tau_LineCenter.dat",
    "clump_tb_grid_file": "clump_Tmb_LineCenter.dat",
    "clump_taufuv_grid_file": "RhoMassAFUV.dat",
    "clump_column_density_file": "clumpMeanCols.dat",
    "clump_temperature_file": "clumpTemperatures_filled.dat",
    "interclump_tau_grid_file": "interclumpTauLineCenter.dat",
    "interclump_tb_grid_file": "interclumpTmbLineCenter.dat",
    "interclump_dust_tau_grid_file": "interclumpDustTau.dat",
    "interclump_dust_tb_grid_file": "interclumpDustSED.dat",
    # 'interclump_tau_grid_file': 'clumpTauLineCenter_reduced.dat',
    # 'interclump_tb_grid_file': 'clumpTmbLineCenter_reduced.dat',
    # 'interclump_dust_tau_grid_file': 'clumpDustTau_reduced.dat',
    # 'interclump_dust_tb_grid_file': 'clumpDustTmb_reduced.dat',
    "interclump_taufuv_grid_file": "interclumpTauFUV.dat",
    "interclump_column_density_file": "interclumpMeanCols.dat",
    "interclump_temperature_file": "interclumpTemperatures_filled.dat",
    "h2_surface_density_file": "h2_surface-density_marasco-bacchini.dat",
    "hi_surface_density_file": "hi_surface-density_marasco-bacchini.dat",
    "h2_scale_height_file": "h2_scale-height_bacchini.dat",
    "hi_scale_height_file": "hi_scale-height_bacchini.dat",
    "like_clumps": False,
    "all_full": True,
    "h_number_density_file": "h_number-density_cubick.dat",
    "fuv_file": "galactic_FUV_complete.dat",
    "average_fuv": False,
    "l_range": (912, 2066),
    "scale_gc": 1.0,
    "mhi_gc": 1.0,
    "mh2_gc": 1.0,
    "new_grid": True,
    # Model information
    "history_path": args.folder,
    "directory": r"/MilkyWay",
    "folder": "",
    "x": 35.6,
    "y": 35.6,
    "z": 0,
    "model_type": "disk",
    # Model parameters
    "resolution": args.resolution,
    # 'molecules': 'all',
    "transitions": [
        "C+ 1",
        "C 1",
        "C 2",
        "CO 1",
        "CO 2",
        "CO 3",
        "CO 4",
        "CO 5",
        "CO 6",
        "CO 7",
        "CO 8",  #'CO 9', 'CO 10',
        "13C+ 1",
        #'13C 1', '13C 2',
        "13CO 1",
        "13CO 2",
    ],  #'13CO 3', '13CO 4', '13CO 5',
    #'13CO 6', '13CO 7', '13CO 8', '13CO 9', '13CO 10',
    #'HCO+ 1', 'HCO+ 2', 'HCO+ 3', 'HCO+ 4', 'HCO+ 5',
    #'H13CO+ 1', 'H13CO+ 2', 'H13CO+ 3', 'H13CO+ 4', 'H13CO+ 5',
    #'H3O+ 1', 'H3O+ 2', 'H3O+ 3', 'H3O+ 4', 'H3O+ 5'],
    # 'dust': 'PAH',
    # "dust": ["240um", "550um"],
    "dust": ["20.um", "30.24um", "70.79um", "158.5um", "177.8um", "240um", "550um", "1.1mm", "1.3mm", "2.4mm", "1.8mm", "3.1mm"],
    "abundances": ["C+", "C", "CO", "13C+", "13C", "13CO"],
    "clump_mass_range": [[0, 2], [-3]],  # [-3]
    "clump_mass_number": [3, 1],  # 1
    "clump_n_max": [1, 100],  # 100
    "clump_log_fuv": None,
    "interclump_log_fuv": None,
    "interclump_wnm_log_fuv": 2,
    "interclump_wnm_ratio": 0,
    "interclump_idx": (False, True),  # True
    "interclump_wnm_idx": (False, False),  # True
    "interclump_density": 19.11,
    "disp_gmc": 0.001,  # 1.1*40**0.38, #0.001,
    "velocity_range": [-350, 350],
    "velocity_number": 701,
    # Flags
    "suggested_calc": True,
    "dilled": False,
    # Property factors
    "hi_mass_factor": 1,
    "h2_mass_factor": 1,
    "ensemble_mass_factor": [1, 1],  # 1
    "fuv_factor": 1,
    "density_factor": 1,
    "interclump_hi_ratio": 1,
    "r_cmz": 0,
    # Logging
    "timed": args.timed,
    "debug": args.debug,
    "verbose": args.verbose,
}

f_mass1 = 1.0
f_mass2 = 1.0
f_m_ens = [1, 1]
f_hi = 1
f_density = 1.0
f_uv = 1
mass_cl = [0, 2]
mass_icl = [-2]
r_cmz = 0
models.constants.fuv_ism = 1
fuv = 1.8

i_max = len(params[0])
print(i_max)

for i, param in enumerate(list(zip(*params))[index:]):

    # Set the required parameters
    for _, p in enumerate(param_keys):
        parameters[p] = param[_]
    parameters["folder"] = param_folders[index:][i]

    # Set model folder directory
    directory = (
        parameters["history_path"]
        + parameters["directory"]
        + "/"
        + parameters["folder"]
    )

    # Check if data exists and will be overwritten
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Print details
    print(f"\n\n   ==> Model {i+index+1} of {i_max}")
    print(f"       {parameters['folder']}")
    print(f"       {'-'*len(parameters['folder'])}")
    if args.verbose:
        print()
        pprint(parameters)
    print()
    print("resolution: ".rjust(30) + f"{parameters['resolution']}")
    print("clump mass range: ".rjust(30) + f"{parameters['clump_mass_range']}")
    print("clump mass number: ".rjust(30) + f"{parameters['clump_mass_number']}")
    print()
    print("R_CMZ: ".rjust(30) + f"{parameters['r_cmz']}")
    print("ensemble mass factor: ".rjust(30) + f"{parameters['ensemble_mass_factor']}")
    print("H2 mass factor: ".rjust(30) + f"{parameters['h2_mass_factor']}")
    print("HI mass factor: ".rjust(30) + f"{parameters['hi_mass_factor']}")
    print("HI interclump ratio: ".rjust(30) + f"{parameters['interclump_hi_ratio']}")
    print("clump density factor: ".rjust(30) + f"{parameters['density_factor']}")
    print("FUV factor: ".rjust(30) + f"{parameters['fuv_factor']}")
    print()
    print("    -- parameters changed")
    for _, p in enumerate(param_keys):
        print(f"  {p} ->".rjust(30) + f" {param[_]}")
    print()

    # Initialise model instance
    kosma = models.Model(**parameters)

    # Run model if selected
    if run_model and (
        (not "species_emissivity.fits" in os.listdir(directory)) or overwrite
    ):
        kosma.calculateModel(kind="linear", multiprocessing=0)
    else:
        print("Will not overwrite model.")

    # Run radiative transfer if selected
    if run_rt and (
        not ("synthetic_intensity.fits" in os.listdir(directory)) or overwrite
    ):
        models.radiativeTransfer.calculateObservation(
            directory=directory,
            dim="spherical",
            multiprocessing=args.mp,
            hi=False,
            vel_pool=False,
            pencil_beam=True,
            slRange=[(-np.pi, np.pi), (-2 * np.pi / 180, 2 * np.pi / 180)],
            nsl=[721, 9],
            terminal=True,
            debug=False,
        )
    else:
        print("Will not overwrite transition datacubes.")
    if (run_rt or hi) and (
        not ("synthetic_hi_intensity.fits" in os.listdir(directory)) or overwrite
    ):
        models.radiativeTransfer.calculateObservation(
            directory=directory,
            dim="spherical",
            multiprocessing=args.mp,
            hi=True,
            vel_pool=False,
            pencil_beam=True,
            slRange=[(-np.pi, np.pi), (-2 * np.pi / 180, 2 * np.pi / 180)],
            nsl=[721, 9],
            terminal=True,
            debug=False,
        )
    else:
        print("Will not overwrite HI datacube.")

    del kosma
