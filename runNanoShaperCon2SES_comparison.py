import os
import subprocess
import time
import statistics
import parmed as pmd
from pathlib import Path
import shutil

################### Settings ########################

# To create the input folders please clone https://github.com/yxwu21/Con2SES/
INPUT_FOLDER = "./Con2SES/datasets/protein_complex_test/"
#INPUT_FOLDER = "./Con2SES/datasets/protein_test/"
#INPUT_FOLDER = "./Con2SES/datasets/nucleic_acid_test/"

# Temporary folder were to store NS results
WORK_DIR = "work"

# Path to NanoShaper binary 
NANOSHAPER_BIN = "NanoShaper"  

# Output file name in the current directory
OUTPUT_TIMES_FILE = "nanoshaper_timings_protein_complex_noclean.txt"

# if you want to skip memory clean up, below in the NS config file please set
#
# Skip_Mem_CleanUp = true
#

#####################################################

os.makedirs(WORK_DIR, exist_ok=True)

times = []

for pdb_file in Path(INPUT_FOLDER).glob("*.pdb"):
    pdb_name = pdb_file.stem
    pdb_path = os.path.abspath(pdb_file)
    run_dir = Path(WORK_DIR) / pdb_name
    run_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing {pdb_name}")
    
    if ("nucleic" in INPUT_FOLDER):
        # Step 4: Convert to XYZR
        xyzr_path = run_dir / "mol.xyzr"
        with open(pdb_path) as f_in, open(xyzr_path, "w") as f_out:
                for line in f_in:
                    if line.startswith(("ATOM", "HETATM")):
                        v = line.split()        
                        x = float(v[5])
                        y = float(v[6])
                        z = float(v[7])
                        r = float(v[8])
                        f_out.write(f"{x:.3f} {y:.3f} {z:.3f} {r:.3f}\n")
    else:
        # Step 1: Write tleap input
        tleap_input = run_dir / "tleap.in"
        tleap_input.write_text(f"""
        source leaprc.protein.ff14SB
        mol = loadpdb {pdb_path}
        saveamberparm mol mol.prmtop mol.inpcrd
        quit
        """)

        # Step 2: Run tleap
        subprocess.run(["tleap", "-f", str(tleap_input.resolve())], cwd=run_dir, check=False)

        # Step 3: Convert to PQR using ParmEd
        parm = pmd.load_file(str(run_dir / "mol.prmtop"), str(run_dir / "mol.inpcrd"))
        pqr_path = run_dir / "mol.pqr"
        parm.save(str(pqr_path), format="pqr", overwrite=True)

        # Step 4: Convert to XYZR
        xyzr_path = run_dir / "mol.xyzr"
        with open(pqr_path) as f_in, open(xyzr_path, "w") as f_out:
                for line in f_in:
                    if line.startswith(("ATOM", "HETATM")):
                        v = line.split()        
                        x = float(v[5])
                        y = float(v[6])
                        z = float(v[7])
                        r = float(v[9])
                        f_out.write(f"{x:.3f} {y:.3f} {z:.3f} {r:.3f}\n")

    # Step 5: Write NanoShaper config
    config_path = run_dir / "conf.prm"
    config_path.write_text(f"""

Compute_Vertex_Normals = true
Save_Mesh_MSMS_Format = false
Load_Balancing = true
Print_Available_Surfaces = false
Grid_scale = 2.0
Grid_perfil = 80.0 
XYZR_FileName = {xyzr_path.resolve()}   
Build_epsilon_maps  = false
Build_status_map = true
Tri2Balls = false
Surface = ses
Smooth_Mesh = true
Number_thread = 8
Skin_Surface_Parameter = 0.45
Blobbyness = -2.5
Skip_Mem_CleanUp = true
Patch_Based_Algorithm = true
Analytical_Ray_Vs_Torus_Intersection = true
Force_Serial_Build = false
Max_Num_Atoms = -1
Domain_Shrinkage = 1.0
Optimize_Grids = true
# Enable or disable cavity detection together with the volume conditional
# filling of voids and cavities
Cavity_Detection_Filling = false

# It is the value of the minimal volume of a cavity to get filled if 
# cavity detection is enabled. 
# The default value is an approximation of the volume of the water molecule 
# default value is 11.4, this is the approximate volume of a water molecule 
# in Angstrom
Conditional_Volume_Filling_Value = 11.4

# If this flag is true, cavities where a sphere of Probe_Radius cannot fit, 
# are removed.
# Use this feature when cavity detection is enabled to filter out bad shaped 
# cavities whose
# volume is higher than Conditional_Volume_Filling_Value.
Keep_Water_Shaped_Cavities = false

# The radius of the sphere that represents a water molecule in Angstrom
# default value is 1.4 Angstrom
Probe_Radius = 1.4
Max_Probes_Self_Intersections = 100
Self_Intersections_Grid_Coefficient = 1.5

# Enable accurate triangulation: if accurate triangulation is enable all points 
# are sampled from the original surface. If disabled the points are not 
# analytically sampled and an high memory saving can be obtained
# together with a 3x speed-up on ray casting if both epsmap is disabled
# MC phase will be slower because vertices are calculated on the fly
Accurate_Triangulation = true

# Perform triangulation using a single ray-casting process. Vertex data is 
# inferred.
Triangulation = true

# Check duplicated vertices when reading
Check_duplicated_vertices = true

# If true save the status map. Enable this for cavity detection and 
# visualization of the coloured FD grid
Save_Status_map = false

# Save Skin/SES in a PovRay file for ray-tracing. 
# This is a purely graphics representation because the surface is not in
# a left handed system as it should in Pov-Ray
Save_PovRay = false
""")

    # Step 6: Run NanoShaper and time it
    start_time = time.time()
    subprocess.run([NANOSHAPER_BIN, str(config_path.resolve())], cwd=run_dir, check=True)
    end_time = time.time()

    elapsed = end_time - start_time
    times.append((pdb_name, elapsed))
    print(f"{pdb_name}: {elapsed:.2f} seconds")

# Step 7: Save results
mean_time = statistics.mean(t for _, t in times)
std_dev = statistics.stdev(t for _, t in times) if len(times) > 1 else 0.0

with open(OUTPUT_TIMES_FILE, "w") as out_file:
    for name, t in times:
        out_file.write(f"{name}\t{t:.4f} s\n")
    out_file.write(f"\nAverage: {mean_time:.4f} s\n")
    out_file.write(f"Std Dev: {std_dev:.4f} s\n")

print(f"\nSaved timing results to: {OUTPUT_TIMES_FILE}")
