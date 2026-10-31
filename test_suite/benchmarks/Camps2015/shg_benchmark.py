#!/usr/bin/env python3
"""Stochastically heated grains benchmark (Camps et al. 2015, A&A 580, A87).

Runs mcfost -benchmark_SHG for each grain type (silicate, graphite, PAH) and
input field (Mathis ISRF x U, diluted blackbodies), builds the emissivity
lambda.eps_lambda [W sr-1 per H] from the temperatures and temperature
probabilities written by mcfost, and compares it with the codes of the paper.

  shg_benchmark.py run     [-w workdir] [-m mcfost] [-f fields] [-g grains] [-t n_T] [-j jobs]
  shg_benchmark.py compare [-w workdir] [-f fields] [-g grains]

The dust model and the results of the paper are extracted once into the data
directory next to this script, from the copy in the MCFOST utils
($MCFOST_UTILS/benchmarks/Camps2015), or else downloaded from
https://shg.ugent.be.
"""

import argparse
import concurrent.futures
import io
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile

import numpy as np
from astropy.io import fits

URL = "https://shg.ugent.be/downloads/"
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# constants of the benchmark (shg.ugent.be/html/_benchmark.html)
h = 6.62606957e-34
c = 2.99792458e8
k = 1.3806488e-23
Rsun_to_au = 6.957e8 / 1.495978707e11  # as mcfost

# input fields: U x J_Mathis and dilution factors of the blackbodies
U_MATHIS = [f"{u:.0e}" for u in 10.0 ** np.arange(-4, 7)]
BLACKBODIES = {3000: 8.28e-12, 6000: 2.23e-13, 9000: 2.99e-14, 12000: 7.23e-15, 15000: 2.36e-15, 18000: 9.42e-16}
FIELDS = [f"Mathis_U_{u}" for u in U_MATHIS] + [f"BlackBody_T_{T:05d}" for T in BLACKBODIES]

# grain type: optical file, size distribution file, calorimetry file, column in the result
# files, number of mcfost sizes (as in the 2015 MCFOST results), size range [um]
GRAINS = {
    "Sil": ("suvSil_121_1201.dat", "ZDA_BARE_GR_S_SzDist_Sil.dat", "Silicate_Calorimetry_1000.dat", 2, 64, (3.55e-4, 0.37)),
    "Gra": ("Gra_121_1201.dat", "ZDA_BARE_GR_S_SzDist_Gra.dat", "Graphitic_Calorimetry_1000.dat", 3, 62, (3.55e-4, 0.33)),
    "PAH": ("PAH_28_1201_neu.dat", "ZDA_BARE_GR_S_SzDist_PAH.dat", "Graphitic_Calorimetry_1000.dat", 4, 24, (3.5e-4, 5e-3)),
}

REFERENCE_CODES = ["DustEM", "SKIRT", "CRT", "DIRTY", "TRADING", "DARTRAY"]

# wavelength range of the comparison [um], and lowest emissivity compared (fraction of the peak)
WL_RANGE = (3.0, 1000.0)
EPS_MIN = 1e-3


def B_lambda(wl, T):
    """Planck function [W m-3 sr-1], wl in m"""
    x = np.minimum(h * c / (wl * k * T), 700.0)
    return 2 * h * c**2 / wl**5 / np.expm1(x)


def download():
    utils = os.environ.get("MCFOST_UTILS", os.path.join(os.path.dirname(DATA), "../../../utils"))
    for name in ("DustModel", "Results_FullSolution"):
        if not os.path.isdir(os.path.join(DATA, name)):
            archive = os.path.join(utils, "benchmarks", "Camps2015", name + ".zip")
            if os.path.isfile(archive):
                print(f"Extracting {archive}")
                zipfile.ZipFile(archive).extractall(DATA)
            else:
                print(f"Downloading {URL}{name}.zip")
                with urllib.request.urlopen(URL + name + ".zip") as r:
                    zipfile.ZipFile(io.BytesIO(r.read())).extractall(DATA)


def read_optical(grain):
    """Q_abs(a, lambda) = Q_ext - Q_sca from the Misselt file (as mcfost), sizes and wavelengths [um]"""
    with open(os.path.join(DATA, "DustModel/GrainInputs", GRAINS[grain][0])) as f:
        lines = f.read().split("\n")
    a, Q, wl = [], [], None
    i = 0
    while i < len(lines):
        if "radius" in lines[i]:
            a.append(float(lines[i].split()[0]))
            block = np.array([l.split() for l in lines[i + 2:i + 2 + 1201]], float)
            Q.append(block[:, 4] - block[:, 3])
            wl = block[:, 1]
            i += 1203
        else:
            i += 1
    return np.array(a), np.array(Q), wl


def mathis_field(wl_mum):
    """J_Mat(lambda) [W m-3 sr-1] (Camps et al. 2015, Eq. 11)"""
    wl = wl_mum * 1e-6
    J = 1e-14 * B_lambda(wl, 7500) + 1e-13 * B_lambda(wl, 4000) + 4e-13 * B_lambda(wl, 3000)
    J = np.where(wl_mum < 0.250, 0.0566 * wl_mum**-1.6678, J)
    J = np.where(wl_mum < 0.134, 1.627, J)
    J = np.where(wl_mum < 0.110, 3069 * wl_mum**3.4172, J)
    return np.where(wl_mum < 0.0912, 0.0, J)


def input_field(field, wl_mum):
    if field.startswith("Mathis"):
        return float(field.split("_")[-1]) * mathis_field(wl_mum)
    T_star = int(field[-5:])
    return BLACKBODIES[T_star] * B_lambda(wl_mum * 1e-6, T_star)


def para_file(grain, field):
    """mcfost parameter file: one thin spherical shell around a 1 Rsun star.
    The field is imposed by -benchmark_SHG: the Mathis field scaled by -Mathis_field,
    or the star diluted at the distance of the shell (dilution factor (R_star/2R)^2)"""
    optical, _, _, _, n_grains, (a_min, a_max) = GRAINS[grain]
    if field.startswith("Mathis"):
        T_star, R = 3000.0, 1000.0
    else:
        T_star = float(field[-5:])
        R = 0.5 / np.sqrt(BLACKBODIES[int(T_star)]) * Rsun_to_au
    # temperature step on the benchmark grid (1201 log points from 1e-3 to 1e4 um) from 0.02 um:
    # mcfost bin edges are the geometric midpoints. The stars have no flux at shorter
    # wavelengths, and the input fields are negligible there
    r = 10 ** (7 / 1200)
    i_min = int(np.ceil(1200 * np.log10(0.02 / 1e-3) / 7))
    n_lambda, wl_min, wl_max = 1201 - i_min, 1e-3 * r ** (i_min - 0.5), 1e4 * r ** 0.5
    return f"""4.1                      mcfost version

#Number of photon packages
  1.28e3                  nbr_photons_eq_th  : T computation
  1.28e3                  nbr_photons_lambda : SED computation
  1.28e3                  nbr_photons_image  : images computation

#Wavelength
  {n_lambda}  {wl_min:.10g} {wl_max:.10g}  n_lambda, lambda_min, lambda_max [mum]
  T F T                   compute temperature?, compute sed?, use default wavelength grid for output ?
  IMLup.lambda            wavelength file (if previous parameter is F)
  F F                     separation of different contributions?, stokes parameters?

#Grid geometry and size
  2                       1 = cylindrical, 2 = spherical
  1 1 1 1                 n_rad (log distribution), nz (or n_theta), n_az, n_rad_in

#Maps
  11 11 3000.             grid (nx,ny), size [AU]
  45.  45.  1  F          RT: imin, imax, n_incl, centered ?
  0    0.   1             RT: az_min, az_max, n_az angles
  140.0                   distance (pc)
  0.                      disk PA

#Scattering method
  2                       1=exact phase function, 2=hg function with same g (2 implies the loss of polarizarion)

#Symmetries
  T                       image symmetry
  T                       central symmetry
  T                       axial symmetry (important only if N_phi > 1)

#Disk physics
  0     0.50  1.0         dust_settling (0=no settling, 1=parametric, 2=Dubrulle, 3=Fromang), exp_strat, a_strat (for parametric settling)
  F                       dust radial migration
  F                       sublimate dust
  F                       hydostatic equilibrium
  F  1e-5                 viscous heating, alpha_viscosity

#Number of zones : 1 zone = 1 density structure + corresponding grain properties
  1

#Density structure
  3                       zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
  1.e-10    100.          dust mass,  gas-to-dust mass ratio
  10.  100.0  2           scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)
  {R * (1 - 1e-3):.8g}  0.0  {R * (1 + 1e-3):.8g}  100.  Rin, edge, Rout, Rc (AU)
  1.125                   flaring exponent, unused for envelope
  0.0  0.0                surface density exponent (or -gamma for tappered-edge disk or volume density for envelope)

#Grain properties
  1  Number of species
  Mie  1 2  0.0  1.0  0.9 Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),  porosity, mass fraction, Vmax (for DHS)
  Misselt_{optical}  1.0  Optical indices file, volume fraction
  3                       Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
  {a_min:.6g} {a_max:.6g} 3.5 {n_grains}   amin, amax [mum], aexp, n_grains (log distribution)

#Molecular RT settings
  T T T                   lpop, laccurate_pop, LTE
  0.05 cs                 Turbulence velocity, unit km/s or cs
  1                       Number of molecules
  co.dat 6                molecular data filename, level max up to which NLTE populations are calculated
  T 1.e-4 abundance.fits.gz   cst molecule abundance ?, abundance, abundance file
  T  2                    ray tracing ?,  number of lines in ray-tracing
  2 3                     transition numbers
  -10.0 10.0 41           vmin and vmax [km/s], number of velocity bins betwen vmin and vmax

#Atoms settings / share some informations with molecules
 1              number of atoms
 H_6.atom       all levels treated in details at the moment
 F              non-LTE ?
 0              initial solution, 0 LTE, 1 from file
 1000 101       vmax (km/s), n_points for ray-traced images and total flux
 T 1            images (T) or total flux (F) ? Number of lines for images
 3 2            upper level -> lower level (Atomic model dependent)

#Star properties
 1 Number of stars
 {T_star:.1f}  1.0  1.0  0.0  0.0  0.0  T Temp, radius (solar radius),M (solar mass),x,y,z (AU), automatic spectrum?
 lte4000-3.5.NextGen.fits.gz
 0.0  2.2  fUV, slope_fUV
"""


def prepare_inputs(workdir):
    """Dust files in the mcfost format: Misselt file with the calorimetry file on the first line"""
    inputs = os.path.join(workdir, "inputs")
    os.makedirs(inputs, exist_ok=True)
    grain_dir = os.path.join(DATA, "DustModel/GrainInputs")
    for optical, _, calorimetry, _, _, _ in GRAINS.values():
        with open(os.path.join(inputs, "Misselt_" + optical), "w") as out, open(os.path.join(grain_dir, optical)) as f:
            out.write(calorimetry + "\n" + f.read())
        shutil.copy(os.path.join(grain_dir, calorimetry), inputs)
    return inputs


def run_one(mcfost, workdir, inputs, grain, field, n_T):
    run_dir = os.path.join(workdir, f"{grain}_{field}")
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(run_dir)
    for f in os.listdir(inputs):
        os.symlink(os.path.join(inputs, f), os.path.join(run_dir, f))
    with open(os.path.join(run_dir, "SHG.para"), "w") as f:
        f.write(para_file(grain, field))
    cmd = [mcfost, "SHG.para", "+disk_struct", "-nT", str(n_T), "-star_bb", "-benchmark_SHG"]
    if field.startswith("Mathis"):
        cmd += ["-Mathis_field", field.split("_")[-1]]
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"), MCFOST_AUTO_UPDATE="0")
    with open(os.path.join(run_dir, "mcfost.log"), "w") as log:
        r = subprocess.run(cmd, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT, env=env)
    return grain, field, r.returncode


def run(mcfost, workdir, grains, fields, n_T=300, jobs=None):
    """Runs mcfost for every grain type and field, in parallel; returns the failed runs"""
    download()
    inputs = prepare_inputs(workdir)
    failed = []
    with concurrent.futures.ThreadPoolExecutor(jobs or os.cpu_count()) as ex:
        futures = [ex.submit(run_one, mcfost, workdir, inputs, g, f, n_T) for f in fields for g in grains]
        for future in concurrent.futures.as_completed(futures):
            grain, field, rc = future.result()
            if rc != 0:
                failed.append(f"{grain}_{field}")
    return failed


def emissivity(workdir, grain, field):
    """lambda.eps_lambda [W sr-1 per H] on the benchmark grid, and the ratio of the emitted to
    the absorbed power for each mcfost grain size"""
    run_dir = os.path.join(workdir, f"{grain}_{field}")
    with fits.open(os.path.join(run_dir, "data_th/Temperature_nRE.fits.gz")) as f:
        T_eq = np.atleast_1d(np.squeeze(f[0].data)).astype(float)
        is_eq = np.atleast_1d(np.squeeze(f[1].data)).astype(bool)
        tab_T = np.squeeze(f[2].data).astype(float)
        P = np.atleast_2d(np.squeeze(f[3].data)).astype(float)
    a = np.atleast_1d(np.squeeze(fits.getdata(os.path.join(run_dir, "data_disk/grain_sizes.fits.gz"))))
    a_lo = np.atleast_1d(np.squeeze(fits.getdata(os.path.join(run_dir, "data_disk/grain_sizes_min.fits.gz"))))
    a_hi = np.atleast_1d(np.squeeze(fits.getdata(os.path.join(run_dir, "data_disk/grain_sizes_max.fits.gz"))))
    if P.shape[0] != a.size:
        P = P.T

    # cross-sections as mcfost: linear in Q, weights in log(a)
    a_file, Q_file, wl_mum = read_optical(grain)
    Q = np.array([np.interp(np.log(a), np.log(a_file), Q_file[:, j]) for j in range(wl_mum.size)])
    C_abs = np.pi * (a * 1e-6) ** 2 * Q  # (n_wl, n_grains) [m2]

    # number of grains per H in each mcfost size bin
    sd = np.loadtxt(os.path.join(DATA, "DustModel/SizeInputs", GRAINS[grain][1]), skiprows=4)
    f_a = np.exp(np.interp(np.log(a), np.log(sd[:, 0]), np.log(sd[:, 1]))) * 1e2  # [m-1 H-1]
    f_a[(a < sd[0, 0] * 0.999) | (a > sd[-1, 0] * 1.001)] = 0.0
    n_H = f_a * (a_hi - a_lo) * 1e-6

    wl = wl_mum * 1e-6
    B_grain = np.where(is_eq[None, :], B_lambda(wl[:, None], np.maximum(T_eq, 1.0)[None, :]),
                       B_lambda(wl[:, None], tab_T[None, :]) @ P.T)
    eps = (C_abs * B_grain) @ n_H * wl

    # energy balance with the bin widths of mcfost (the Mathis field is discontinuous)
    dwl = np.gradient(wl)[:, None]
    absorbed = (C_abs * input_field(field, wl_mum)[:, None] * dwl).sum(axis=0)
    emitted = (C_abs * B_grain * dwl).sum(axis=0)
    return wl_mum, eps, emitted / np.maximum(absorbed, 1e-300)


def reference(code, field):
    path = os.path.join(DATA, "Results_FullSolution", f"{code}_SHG_{field}.dat")
    return np.loadtxt(path, comments="#") if os.path.exists(path) else None


def compare(workdir, grains, fields, codes=REFERENCE_CODES[:3], out=sys.stdout):
    """mcfost / reference code over WL_RANGE where the reference is above EPS_MIN of its peak,
    and emitted / absorbed power per grain size. Returns {(field, grain): {name: ratios}}"""
    results = {}
    fmt = lambda r: f"{np.median(r):6.3f} [{r.min():5.2f},{r.max():5.2f}]"
    print(f"{'field':20s} {'grain':5s} {'E_em/E_abs':>20s} " + " ".join(f"{c:>20s}" for c in codes), file=out)
    for field in fields:
        for grain in grains:
            wl, eps, balance = emissivity(workdir, grain, field)
            ratios = {"E_em/E_abs": balance}
            col = GRAINS[grain][3]
            for code in codes:
                ref = reference(code, field)
                if ref is not None:
                    ok = (wl >= WL_RANGE[0]) & (wl <= WL_RANGE[1]) & (ref[:, col] > ref[:, col].max() * EPS_MIN)
                    ratios[code] = eps[ok] / ref[ok, col]
            results[(field, grain)] = ratios
            print(f"{field:20s} {grain:5s} " + " ".join(fmt(ratios[c]) if c in ratios else f"{'-':>20s}"
                                                   for c in ["E_em/E_abs"] + list(codes)), file=out)
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("action", choices=["run", "compare"])
    p.add_argument("-w", "--workdir", default="Camps2015_SHG")
    p.add_argument("-m", "--mcfost", default=shutil.which("mcfost"))
    p.add_argument("-f", "--fields", default="all", help="comma-separated, e.g. Mathis_U_1e+00,BlackBody_T_03000")
    p.add_argument("-g", "--grains", default="Sil,Gra,PAH")
    p.add_argument("-t", "--n_T", type=int, default=300, help="number of temperature bins (mcfost -nT)")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="mcfost runs in parallel")
    args = p.parse_args()

    workdir = os.path.abspath(args.workdir)
    fields = FIELDS if args.fields == "all" else args.fields.split(",")
    grains = args.grains.split(",")
    if args.action == "run":
        failed = run(os.path.abspath(args.mcfost), workdir, grains, fields, args.n_T, args.jobs)
        if failed:
            print("FAILED (see mcfost.log):", " ".join(failed))
    download()
    compare(workdir, grains, fields)


if __name__ == "__main__":
    main()
