import astropy.io.fits as fits
import numpy as np
import subprocess
import pytest
import shutil
import glob
import os

_mcfost_bin = "../src/mcfost"

# Get list of models using directory names (filtering out hidden and AppleDouble files)
model_list = [
    d for d in os.listdir("test_data/")
    if not d.startswith(".") and os.path.isdir(os.path.join("test_data", d))
]

# If running on CI, only run some of the tests
#if os.environ.get('CI', None) == 'true':
#    model_list = ["ref3.0","ref3.0_multi","debris","discF_00500"]

#model_list = ["ref3.0"]#,"debris","discF_00500"]

wl_list = ["1.0","10","100","1000"]
wl_list_pola = ["1.0","1000"]
wl_list_contrib = ["1.0","100","1000"]


# flag to skip the calculations and set-up the thresholds in a easy way
compute_models = True

# change this to test results downloaded from other machines
test_dir = "."

if test_dir != ".":
    compute_models = False

def mcfost(filename,opt=""):
    cmd = _mcfost_bin+" "+filename+" "+opt
    result = subprocess.call(cmd.split())

def clean_results(model_name):
    cmd = ["rm","-rf",model_name," *.tmp"]
    result = subprocess.call(cmd)

def all_almost_equal(x,y,threshold=0.01):
    # test if all the values of two arrays are almost equal
    return (abs((x-y)) < threshold * x).all()

def MC_similar(x,y,threshold=0.01,mask_threshold=1e-24):
    # test if two arrays have the same at the 75% percentile
    # ignoring values that are very small as they are very noisy

    #mask = (abs(x) < mask_threshold) | (abs(y) < mask_threshold) # Unit: W.m-2 or W.m-2/pixel
    mask = abs(x) < mask_threshold # Unit: W.m-2 or W.m-2/pixel : only using reference for the mask
    x_ma = np.ma.masked_where(mask, x)
    y_ma = np.ma.masked_where(mask, y)

    #return (abs((x_ma-y_ma)/x_ma).mean() < threshold)
    print("MC relative difference is:", np.percentile(abs((x_ma-y_ma)/x_ma).compressed(), 75))
    return ( np.percentile(abs((x_ma-y_ma)/x_ma).compressed(), 75)  < threshold )

def test_mcfost_bin():
    # We first test if the mcfost binary actually exists and runs
    assert os.path.isfile(_mcfost_bin)
    try:
        subprocess.call(_mcfost_bin)
    except OSError as e:
        assert not e.errno == os.errno.ENOENT

@pytest.mark.parametrize("model_name", model_list)
def test_Temperature(model_name):
    if compute_models:
        clean_results(model_name) # removing all previous calculations

        # Run the mcfost model
        filename = "test_data/"+model_name+"/"+model_name+".para"
        if (model_name == "discF_00500"):
            opt=" -phantom test_data/"+model_name+"/"+model_name+" -not_random_Voronoi"
        else:
            opt=""
        mcfost(filename,opt="-mol -root_dir "+model_name+opt)

    # Read the results
    T_name = model_name+"/data_th/Temperature.fits.gz"
    T = fits.getdata(test_dir+"/"+T_name)
    T_ref = fits.getdata("test_data/"+T_name)

    print("Maximum T difference", (abs(T-T_ref)/(T_ref+1e-30)).max())
    print("Mean T difference   ", (abs(T-T_ref)/(T_ref+1e-30)).mean())

    assert MC_similar(T_ref,T,threshold=0.05)

@pytest.mark.parametrize("model_name", model_list)
def test_SED(model_name):
    # Re-use the previous mcfost model

    # Read the results
    SED_name = model_name+"/data_th/sed_rt.fits.gz"
    if (not os.path.isfile(SED_name)):
        pytest.skip("No SED found")
    SED = fits.getdata(test_dir+"/"+SED_name)[0,:,:,:]
    SED_ref = fits.getdata("test_data/"+SED_name)[0,:,:,:]

    if model_name == "ref3.0_multi":
        threshold=0.12
    else:
        threshold=0.1

    print("Maximum SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=threshold)

@pytest.mark.parametrize("model_name", model_list)
def test_SED_mc(model_name):
    # Re-use the previous mcfost model
    # Step-2 MC SED: catches step-1 photons leaking into step 2 (lsed_complete)

    # Read the results
    SED_name = model_name+"/data_th/sed_mc.fits.gz"
    if (not os.path.isfile("test_data/"+SED_name)):
        pytest.skip("No reference MC SED found")
    SED = fits.getdata(test_dir+"/"+SED_name)[0,:,:,:]
    SED_ref = fits.getdata("test_data/"+SED_name)[0,:,:,:]

    print("Maximum MC SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean MC SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=0.15)

@pytest.mark.parametrize("model_name", model_list)
def test_SED_th(model_name):
    # Re-use the previous mcfost model
    # Step-1 SED must not be written empty (e.g. after hydrostatic iterations)

    SED_name = model_name+"/data_th/.sed_th.fits.gz"
    if (not os.path.isfile(test_dir+"/"+SED_name)):
        pytest.skip("No step-1 SED found")
    SED = fits.getdata(test_dir+"/"+SED_name)

    assert np.isfinite(SED).all()
    assert SED.max() > 0.0

@pytest.mark.parametrize("model_name", model_list)
def test_SED_contrib(model_name):
    # Re-use the previous mcfost model

    # Read the results
    SED_name = model_name+"/data_th/sed_rt.fits.gz"
    if (not os.path.isfile(SED_name)):
        pytest.skip("No SED found")
    SED = fits.getdata(test_dir+"/"+SED_name)[1:,:,:,:]
    SED_ref = fits.getdata("test_data/"+SED_name)[1:,:,:,:]

    print("Maximum SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=0.15)

@pytest.mark.parametrize("model_name", model_list)
def test_mol_map(model_name):
    # Re-use the previous mcfost model

    # Read the results
    image_name = model_name+"/data_CO/lines.fits.gz"
    image = fits.getdata(test_dir+"/"+image_name)
    image_ref = fits.getdata("test_data/"+image_name)

    print("Maximum mol map difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean mol map difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.1)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list)
def test_image(model_name, wl):

    if os.environ.get('CI', None) == 'true':
        if (model_name == "discF_00500"):
            # This requires too much memory on CI
            pytest.skip("Skipping images on phantom dump on CI")

    if compute_models:
        # Run the mcfost model
        filename = "test_data/"+model_name+"/"+model_name+".para"
        if (model_name == "discF_00500"):
            opt=" -phantom test_data/"+model_name+"/"+model_name+" -not_random_Voronoi"
        else:
            opt=""
        mcfost(filename,opt="-img "+wl+" -root_dir "+model_name+opt)

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"

    hdr = fits.getheader("test_data/"+image_name)
    n_incl = hdr['NAXIS3']

    threshold=0.1
    if (model_name == "ref3.0") and (wl == "100"): # weird difference on linux, ifort, openmp=no, release=no
        threshold=0.11

    if (model_name == "ref3.0_multi") and (wl == "1.0"): # weird difference on linux, ifort, openmp=yes, release=no
        threshold=0.11

    if (model_name == "ref4.1_nLTE") and (wl == "100"):
        threshold=0.13

    if (model_name == "ref4.1_nLTE") and (wl == "10"):
        threshold=0.12

    # We just keep intensity
    for i in range(n_incl):
        image = np.nan_to_num(fits.getdata(test_dir+"/"+image_name))
        image_ref = fits.getdata("test_data/"+image_name)

        print(model_name, " wl=", wl)
        print("-- i=", i)
        print(image.shape)
        print(image_ref.shape)

        image = image[0,:,i,:,:]
        image_ref = image_ref[0,:,i,:,:]

        print("Min ref", image_ref.min(), image_ref.max())
        print("Min    ", image.min(), image.max())

        print("Maximum image difference", i, (abs(image-image_ref)/(image_ref+1e-30)).max())
        print("Mean image difference   ", i, (abs(image-image_ref)/(image_ref+1e-30)).mean())

        print("Threshold=", threshold)
        print("MC similar", MC_similar(image_ref,image,threshold=threshold))
        print("-------------")

    image = fits.getdata(test_dir+"/"+image_name)
    image_ref = fits.getdata("test_data/"+image_name)

    image = image[0,:,:,:,:]
    image_ref = image_ref[0,:,:,:,:]

    assert MC_similar(image_ref,image,threshold=threshold)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list_pola)
def test_pola(model_name, wl):
    # Re-use previous calculation

    if os.environ.get('CI', None) == 'true':
        if (model_name == "discF_00500"):
            # This requires too much memory on CI
            pytest.skip("Skipping images on phantom dump on CI")

    if (model_name == "ref3.0_multi"):
        pytest.skip("Skipping pola for ref3.0_multi for now")

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"
    image = fits.getdata(test_dir+"/"+image_name)
    if ((image.shape[0] != 4) & (image.shape[0] != 8)):
        pytest.skip("No polarisation found")
    image_ref = fits.getdata("test_data/"+image_name)

    # We just keep Stokes Q, U
    image = image[[1,2],:,:,:,:]
    image_ref = image_ref[[1,2],:,:,:,:]

    print("Maximum pola difference", (abs(image-image_ref)/(image_ref+1e-36)).max())
    print("Mean pola difference   ", (abs(image-image_ref)/(image_ref+1e-36)).mean())

    if model_name == "debris":
        mask_threshold = 1e-32
    elif model_name == "ref4.1_PAH":
        mask_threshold = 1e-28
    elif model_name == "ref4.1_hydrostatic":
        mask_threshold = 1e-28
    else:
        mask_threshold = 1e-21

    threshold=0.1
    if (model_name == "ref4.1_nLTE") and (wl == "1000"):
        threshold=0.18

    assert MC_similar(image_ref,image,threshold=threshold,mask_threshold=mask_threshold)

@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list_contrib)
def test_contrib(model_name, wl):
    # Re-use previous calculation

    if os.environ.get('CI', None) == 'true':
        if (model_name == "discF_00500"):
            # This requires too much memory on CI
            pytest.skip("Skipping images on phantom dump on CI")

    if (model_name == "ref3.0_multi") and (wl == "1.0"):
        pytest.skip("Skipping contrib at 1.0 for ref3.0_multi for now")

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"
    image = fits.getdata(test_dir+"/"+image_name)
    if ((image.shape[0] != 5) & (image.shape[0] != 8)):
        pytest.skip("No contrib found")
    image_ref = fits.getdata("test_data/"+image_name)

    # We just keep separate contributions
    image = image[[4,5,6,7],:,:,:,:]
    image_ref = image_ref[[4,5,6,7],:,:,:,:]

    print("Maximum contrib difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean contrib difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    if model_name == "ref3.0_multi":
        mask_threshold=1e-19
    elif (model_name == "ref4.1_PAH"):
        mask_threshold=1e-22
    else:
        mask_threshold=1e-23

    threshold=0.1
    if (model_name == "ref3.0") and (wl == "100"): # weird difference on linux, ifort, openmp=no, release=no
        threshold=0.11

    if (model_name == "ref4.1_nLTE") and (wl == "100"):
        threshold=0.13

    assert MC_similar(image_ref,image,threshold=threshold,mask_threshold=mask_threshold)


def test_prodimo():
    # ProDiMo mode must run through step 2 (xN_abs is owned by setup_ProDiMo)
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_prodimo"
    clean_results(root_dir)
    filename = "test_data/ref3.0/ref3.0.para"
    result = subprocess.call((_mcfost_bin+" "+filename+" -prodimo -root_dir "+root_dir).split())
    assert result == 0
    assert len(glob.glob(root_dir+"/data_ProDiMo*/forProDiMo.fits.gz")) > 0


def test_hydrostatic_restart():
    # Restarting a hydrostatic model past n_iter_physics must still produce a
    # temperature and step-1 SED, instead of writing empty arrays
    if not compute_models:
        pytest.skip("Only run when computing models")
    model_name = "ref4.1_hydrostatic"
    if (not os.path.isfile(model_name+"/data_th/Temperature_iter0001.fits.gz")):
        pytest.skip("No hydrostatic run to restart from")
    root_dir = model_name+"_restart"
    clean_results(root_dir)
    shutil.copytree(model_name, root_dir)
    filename = "test_data/"+model_name+"/"+model_name+".para"
    result = subprocess.call((_mcfost_bin+" "+filename+" -n_iter_physics 1 -root_dir "+root_dir).split())
    assert result == 0

    SED = fits.getdata(root_dir+"/data_th/.sed_th.fits.gz")
    assert SED.max() > 0.0

    T = fits.getdata(root_dir+"/data_th/Temperature.fits.gz")
    T_ref = fits.getdata("test_data/"+model_name+"/data_th/Temperature.fits.gz")
    assert MC_similar(T_ref,T,threshold=0.05)


def test_mol_spectrum_circular_pixels():
    # Molecular spectrum with circular pixels (1x1 map, RT_line_method = 1) and a
    # symmetric image: same spectrum with 1 thread and with all threads (F5), and
    # symmetric in velocity. The symmetrisation used to write outside the spectrum
    # array (old velocity indexing), giving garbage fluxes and crashes.
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_mol_spectrum"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/spectrum.para"
    make_para("test_data/ref3.0/ref3.0.para", para, {"grid (nx,ny)": set_first_value("1 1")})

    spectra = []
    for threads, run_dir in [(1, root_dir+"/t1"), (None, root_dir+"/tN")]:
        # temperature of the reference model
        os.makedirs(run_dir+"/data_th")
        shutil.copy("test_data/ref3.0/data_th/Temperature.fits.gz", run_dir+"/data_th")
        run_mcfost(para, run_dir, "-mol -no_T", threads=threads)
        spectra.append(fits.getdata(run_dir+"/data_CO/lines.fits.gz").astype(float))

    s1, sN = spectra
    assert np.all(np.isfinite(s1)) and np.all(s1 > 0)
    # the stellar flux is computed with random rays: small differences between thread counts
    assert np.abs(sN - s1).max() < 1e-3 * s1.max()
    # axes (azimuth, inclination, transition, velocity, y, x): symmetric in velocity
    assert np.allclose(s1, np.flip(s1, axis=3), rtol=1e-5)


def test_mol_tau_surface():
    # Molecular tau = 1 surface map (-tau_surface) with several inclinations: every
    # inclination has surface pixels and the map is identical with 1 thread and with all
    # threads. The map used to be reset inside the parallel region at each inclination
    # (earlier inclinations erased, race between threads, F6), and was kept with the
    # 1x1 pixel dimensions of the dust SED ray tracing (out-of-bounds writes).
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_mol_tau_surface"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/tau_surface.para"
    make_para("test_data/ref3.0/ref3.0.para", para, {"grid (nx,ny)": set_first_values("21 21")})

    maps = []
    for threads, run_dir in [(1, root_dir+"/t1"), (None, root_dir+"/tN")]:
        # temperature of the reference model
        os.makedirs(run_dir+"/data_th")
        shutil.copy("test_data/ref3.0/data_th/Temperature.fits.gz", run_dir+"/data_th")
        run_mcfost(para, run_dir, "-mol -no_T -tau_surface 1", threads=threads)
        maps.append(fits.getdata(run_dir+"/data_CO/tau=1_surface.fits.gz"))

    # axes (x y z, azimuth, inclination, y, x)
    for i in range(maps[0].shape[2]):
        assert (maps[0][:,0,i] != 0).any(), f"inclination {i}: empty tau = 1 surface"
    assert np.array_equal(maps[0], maps[1])


def test_mol_nlte():
    # Molecular NLTE populations (ref3.0 with LTE = F on a coarse grid, without the
    # accurate-population step): each step of the iterations converges in a few
    # iterations. The convergence criterion used to keep the maximum of all previous
    # iterations (maxdiff reset inside the parallel region with a max reduction, F1) and
    # never converged. The dust opacity in J_mol_loc was read out of bounds (F14): only
    # detected by a bounds-checked build.
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_mol_nlte"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/nlte.para"
    make_para("test_data/ref3.0/ref3.0.para", para,
              {"n_rad (log distribution)": set_first_values("20 10 1 4"),
               "lpop, laccurate_pop, LTE": set_first_values("T F F")})

    stdout = run_mcfost(para, root_dir, "-mol", timeout=300)

    # "Step i Iteration n" ... "Relative difference = x" ... "Threshold = y"
    steps = {}
    for line in stdout.splitlines():
        words = line.split()
        if line.startswith(" Step") and "Iteration" in line:
            step, iteration = int(words[1]), int(words[3])
        elif "Relative difference" in line:
            diff = float(words[-1])
        elif "Threshold" in line:
            steps[step] = (iteration, diff, float(words[-1]))
    assert sorted(steps) == [1, 2]
    for step, (n_iter, diff, threshold) in steps.items():
        assert n_iter <= 10, f"step {step}: {n_iter} iterations"
        assert diff < threshold, f"step {step}: relative difference {diff}"

    lines = fits.getdata(root_dir+"/data_CO/lines.fits.gz")
    assert np.all(np.isfinite(lines)) and lines.max() > 0


def test_seed_image():
    # With -seed, results are written in root_dir/seed=<n>: an image run with the same
    # seed must read the temperature computed there (J5: Tfile ignored seed_dir)
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_seed"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/seed.para"
    make_para("test_data/ref3.0/ref3.0.para", para, {})

    run_mcfost(para, root_dir, "-seed 7", timeout=300)
    assert os.path.isfile(root_dir+"/seed=7/data_th/Temperature.fits.gz")

    stdout = run_mcfost(para, root_dir, "-seed 7 -img 1.0", timeout=300)
    assert "seed=7/data_th/Temperature.fits.gz" in stdout
    assert os.path.isfile(root_dir+"/seed=7/data_1.0/RT.fits.gz")


@pytest.mark.parametrize("args", ["0 1e-8", "150 1e-8", "2 1e-8", "x 1e-8", "1 abc"])
def test_mdot_option_invalid(args):
    # -Mdot <star> <Mdot>: the star number used to index a fixed-size array without
    # check (H6: -Mdot 0 or -Mdot 150 wrote outside the array); must stop with an error
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_mdot"
    clean_results(root_dir)
    result = subprocess.run((os.path.abspath(_mcfost_bin)+" test_data/ref3.0/ref3.0.para -op 1.0 -root_dir "
                             +root_dir+" -Mdot "+args).split(), capture_output=True, text=True, timeout=60)
    assert result.returncode != 0
    assert "-Mdot" in result.stdout


def test_mdot_option():
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_mdot"
    clean_results(root_dir)
    stdout = run_mcfost("test_data/ref3.0/ref3.0.para", root_dir, "-op 1.0 -Mdot 1 1e-8", timeout=60)
    assert "Mdot=" in stdout


# Pascucci et al. 2004 benchmark: SEDs and midplane temperatures compared with the
# RADMC and RADICAL results of the paper (distributed with the MCFOST utils).
# Optically thin cases run for every pull request, optically thick ones with --run-slow.
# Largest differences with MCFOST 4.1.16 (October 2026), SED (all wavelengths, >= 1 um
# for tau = 100 where some RADMC optical fluxes are 0) and midplane temperature: 3.9 % and
# 2.1 % (tau = 0.1), 8.1 % and 2.7 % (tau = 1), 15.6 % and 3.9 % (tau = 10), 15.6 % and
# 9.8 % (tau = 100); RADMC and RADICAL differ by up to 8, 23 and 7 % for tau = 1, 10, 100.
_pascucci_cases = [
    pytest.param("0.1", "01tau", 0.06, 0.0, 0.05),
    pytest.param("1", "1tau", 0.10, 0.0, 0.05),
    pytest.param("10", "10tau", 0.20, 0.0, 0.06, marks=pytest.mark.slow),
    pytest.param("100", "100tau", 0.20, 1.0, 0.12, marks=pytest.mark.slow),
]


def _pascucci_results_dir():
    utils = os.environ.get("MCFOST_UTILS", os.path.join("..", "utils"))
    return os.path.join(utils, "benchmarks", "Pascucci2004")


def _read_pascucci_seds(filename):
    # lambda [um], inclination [deg], lambda.L_lambda [erg/s]; 2 header lines (n_lambda, n_incl)
    rows = [l.split() for l in open(filename) if l.strip() and not l.lstrip().startswith(";")][2:]
    return np.array([[float(x) for x in r[:3]] for r in rows if len(r) >= 3])


@pytest.mark.parametrize("tau, ref_dir, sed_tol, lambda_min, temp_tol", _pascucci_cases)
def test_pascucci_benchmark(tau, ref_dir, sed_tol, lambda_min, temp_tol):
    if not compute_models:
        pytest.skip("Only run when computing models")
    results = _pascucci_results_dir()
    if not os.path.isdir(results):
        pytest.skip("Pascucci 2004 results not found in the MCFOST utils (mcfost -update_utils)")

    root_dir = "Pascucci2004_tau="+tau
    clean_results(root_dir)
    # the benchmark mode (isotropic scattering, tabulated cross-sections) is selected by the
    # name of the parameter file; +disk_struct writes the grid for the temperature comparison
    run_mcfost("benchmarks/Pascucci2004/Pascucci_tau="+tau+".para", root_dir, "-star_bb +disk_struct",
               timeout=3600)

    sed = fits.getdata(root_dir+"/data_th/sed_rt.fits.gz").astype(np.float64)
    temperature = fits.getdata(root_dir+"/data_th/Temperature.fits.gz").astype(np.float64)
    r_midplane = fits.getdata(root_dir+"/data_disk/grid.fits.gz")[0,0,0,:] # au

    au_to_cm = 1.495978707e13
    distance = 140 * 3.0856775814913673e18 # cm
    for code in ["RADMC", "RADICAL"]:
        ref = _read_pascucci_seds(results+"/"+code+"/"+ref_dir+"/allseds.dat")
        for k, incl in enumerate([12.5, 77.5]):
            ref_incl = ref[ref[:,1] == incl]
            # lambda.F_lambda [W/m2] at 140 pc -> lambda.L_lambda [erg/s]
            l_lambda = sed[0,0,k,:] * 1e3 * 4*np.pi * distance**2
            ratio = l_lambda / ref_incl[:,2]
            mask = ref_incl[:,0] >= lambda_min
            print(code, "i =", incl, "SED ratio:", ratio[mask].min(), ratio[mask].max())
            assert np.all(abs(ratio[mask] - 1) < sed_tol), \
                f"SED differs from {code} at i = {incl} deg by more than {sed_tol}"

        # midplane temperature: (r [cm], theta [rad], T [K]), theta closest to pi/2
        dtemp = np.loadtxt(results+"/"+code+"/"+ref_dir+"/dtemp.dat", skiprows=1)
        theta = np.unique(dtemp[:,1])
        theta_mid = theta[np.argmin(abs(theta - np.pi/2))]
        mid = dtemp[dtemp[:,1] == theta_mid]
        r = mid[:,0] / au_to_cm
        T_mcfost = np.interp(np.log(r), np.log(r_midplane), temperature[0,:])
        ratio = T_mcfost / mid[:,2]
        print(code, "midplane temperature ratio:", ratio.min(), ratio.max())
        assert np.all(abs(ratio - 1) < temp_tol), \
            f"midplane temperature differs from {code} by more than {temp_tol}"


# Pinte et al. 2009 benchmark (very optically thick disc, anisotropic scattering): SEDs
# compared with MCMax and MCFOST 2009, midplane temperatures with MCMax, TORUS and
# MCFOST 2009, vertical temperatures at 0.2 and 200 au with MCMax and MCFOST 2009
# (results distributed with the MCFOST utils). Long tests (--run-slow).
# Largest differences with MCFOST 4.1.16 over several runs (October 2026):
# - SED per wavelength (points above 1e-3 of the maximum): 16 % (tau = 1e3, i <= 76 deg),
#   32 % (tau = 1e6, i <= 70 deg, short wavelengths at 18 deg); closer to edge-on the short
#   wavelengths differ by up to a factor 3 between all codes, only the median over
#   wavelength is compared there (8 % at most, at 81.4 deg where the star is seen through
#   the disc surface);
# - temperatures inside the disc: 7.6 % (tau = 1e3), 19 % (tau = 1e6, inner au, where
#   MCMax and TORUS differ by a similar amount).
# tau = 1e6 also compares the 1 um images (stellar pixels excluded) with MCMax and MCFOST 2009
# at the 10 inclinations, and TORUS at 70 and 87 deg: scattered flux 7.5 % (MCMax), cuts
# along the major axis 15 % (i <= 64 deg) and medians of the major and minor axis cuts 14 %,
# flux-weighted linear polarisation degree 24 % (MCMax is 10-20 % lower, MCFOST 2009 before
# the multiple scattering fix of 2026 about 8 % higher).
_pinte2009_cases = [
    pytest.param("1e3", 76.0, 0.25, 0.12, "", marks=pytest.mark.slow),
    pytest.param("1e6", 70.0, 0.40, 0.25, "-diff_approx", marks=pytest.mark.slow),
]
_pinte2009_inclinations = ["18.2", "31.8", "41.4", "49.5", "56.6", "63.3", "69.5", "75.5", "81.4", "87.1"]


@pytest.mark.parametrize("tau, incl_max, sed_tol, temp_tol, opt", _pinte2009_cases)
def test_pinte2009_benchmark(tau, incl_max, sed_tol, temp_tol, opt):
    if not compute_models:
        pytest.skip("Only run when computing models")
    utils = os.environ.get("MCFOST_UTILS", os.path.join("..", "utils"))
    results = os.path.join(utils, "benchmarks", "Pinte2009")
    if not os.path.isdir(results):
        pytest.skip("Pinte 2009 results not found in the MCFOST utils (mcfost -update_utils)")

    root_dir = "Pinte2009_tau="+tau
    clean_results(root_dir)
    # -diff_approx (tau = 1e6): temperature of the optically thick midplane, not reached by
    # the packets, from the diffusion approximation (1 K otherwise)
    run_mcfost("benchmarks/Pinte2009/Pinte2009_tau="+tau+".para", root_dir,
               "-star_bb +disk_struct "+opt, timeout=3600)

    sed_file = fits.open(root_dir+"/data_th/sed_rt.fits.gz")
    sed = sed_file[0].data.astype(np.float64)        # lambda.F_lambda [W/m2] at 140 pc
    wl = sed_file[1].data.astype(np.float64)         # micron
    for k, incl in enumerate(_pinte2009_inclinations):
        for code, filename in [("MCMax", "mcmax/tau="+tau+"_aniso/SED/SED"+incl+".dat"),
                               ("MCFOST 2009", "mcfost/tau="+tau+"_aniso/SED/tau="+tau+"_aniso_inc"
                                +"{:03d}".format(round(np.degrees(np.arccos(0.95 - 0.1*k))))+".dat")]:
            ref = np.loadtxt(results+"/"+filename, comments="#")
            sed_interp = np.exp(np.interp(np.log(ref[:,0]), np.log(wl), np.log(np.maximum(sed[0,0,k,:], 1e-300))))
            mask = (ref[:,1] > 1e-3 * ref[:,1].max()) & (ref[:,0] >= wl.min()) & (ref[:,0] <= wl.max())
            ratio = sed_interp[mask] / ref[mask,1]
            print(code, "i =", incl, "SED ratio:", ratio.min(), ratio.max(), "median", np.median(ratio))
            assert abs(np.median(ratio) - 1) < 0.12, f"SED differs from {code} at i = {incl} deg"
            if float(incl) <= incl_max:
                assert np.all(abs(ratio - 1) < sed_tol), f"SED differs from {code} at i = {incl} deg"

    temperature = root_dir+"/data_th/Temperature_DA.fits.gz"
    if not os.path.isfile(temperature):
        temperature = root_dir+"/data_th/Temperature.fits.gz"
    T = fits.getdata(temperature).astype(np.float64) # (z, r)
    grid = fits.getdata(root_dir+"/data_disk/grid.fits.gz")
    r = grid[0,0,0,:]

    # midplane: radius [au], T [K]
    for code, filename in [("MCMax", "mcmax/tau="+tau+"_aniso/T/midplane_temp_tau"+tau+".dat"),
                           ("TORUS", "torus/tau="+tau+"_aniso/T/tau="+tau+"_aniso_radial.dat"),
                           ("MCFOST 2009", "mcfost/tau="+tau+"_aniso/T/Radial_temp_tau="+tau+".txt")]:
        ref = np.loadtxt(results+"/"+filename, comments="#")
        ref = ref[(ref[:,0] >= r.min()) & (ref[:,0] <= r.max())] # some codes extend beyond the disc
        ratio = np.interp(np.log(ref[:,0]), np.log(r), T[0,:]) / ref[:,1]
        print(code, "midplane temperature ratio:", ratio.min(), ratio.max())
        assert np.all(abs(ratio - 1) < temp_tol), f"midplane temperature differs from {code}"

    # vertical cuts: z [au], T [K], at R = 0.2 and 200 au (z of the cells interpolated in log r)
    for R, label in [(0.2, "0.2AU"), (200., "200AU")]:
        i = np.searchsorted(r, R)
        f = np.log(R/r[i-1]) / np.log(r[i]/r[i-1])
        z = np.exp((1-f) * np.log(grid[1,0,:,i-1]) + f * np.log(grid[1,0,:,i]))
        T_R = (1-f) * T[:,i-1] + f * T[:,i]
        for code, filename in [("MCMax", "mcmax/tau="+tau+"_aniso/T/vertical_temp_"+label+"_tau"+tau+".dat"),
                               ("MCFOST 2009", "mcfost/tau="+tau+"_aniso/T/Vertical_temp_R="+label+"_tau="+tau+".txt")]:
            ref = np.loadtxt(results+"/"+filename, comments="#")
            z_ref = abs(ref[:,0])
            mask = (z_ref >= z.min()) & (z_ref <= z.max())
            ratio = np.interp(np.log(z_ref[mask]), np.log(z), T_R) / ref[mask,1]
            print(code, "temperature ratio at R =", R, "au:", ratio.min(), ratio.max())
            assert np.all(abs(ratio - 1) < temp_tol), f"temperature at R = {R} au differs from {code}"

    if tau == "1e6":
        _check_pinte2009_images(root_dir, results)


def _check_pinte2009_images(root_dir, results):
    # 1 um images and linear polarisation degree of the tau = 1e6 model
    run_mcfost("benchmarks/Pinte2009/Pinte2009_tau=1e6.para", root_dir, "-star_bb -img 1.0", timeout=3600)
    image = fits.getdata(root_dir+"/data_1.0/RT.fits.gz").astype(np.float64) # (Stokes, az, incl, y, x)
    c = image.shape[-1] // 2

    def properties(I, P):
        # stellar pixels excluded; cuts averaged over 5 pixels along the major (x) and minor (y) axes
        disc = I.copy()
        disc[c-3:c+4,c-3:c+4] = 0.0
        flux = disc.sum()
        return flux, disc[c-2:c+3,:].mean(axis=0), disc[:,c-2:c+3].mean(axis=1), np.sum(P * disc) / flux

    for k, incl in enumerate(_pinte2009_inclinations):
        n = round(np.degrees(np.arccos(0.95 - 0.1*k)))
        I = image[0,0,k]
        P = np.sqrt(image[1,0,k]**2 + image[2,0,k]**2) / np.maximum(I, 1e-300)
        flux, major, minor, P_mean = properties(I, P)
        refs = [("MCFOST 2009", "mcfost/tau=1e6_aniso/images/intensity_map_inc{:03d}.fits.gz".format(n),
                 "mcfost/tau=1e6_aniso/pola/linear_polarization_map_inc{:03d}.fits.gz".format(n), 1.0),
                ("MCMax", "mcmax/tau=1e6_aniso/images/ImageF_i"+incl+"_l00001.00_fov00900.0.fits",
                 "mcmax/tau=1e6_aniso/pola/ImageP_i"+incl+"_l00001.00_fov00900.0.fits", 1.0)]
        if n in (70, 87): # TORUS: polarisation degree in %
            refs.append(("TORUS", "torus/tau=1e6_aniso/images/tau=1e6_aniso_inc{:03d}_natural_image001.fits".format(n),
                         "torus/tau=1e6_aniso/pola/tau=1e6_aniso_inc{:03d}_natural_pol001.fits".format(n), 0.01))
        for code, image_file, pola_file, P_unit in refs:
            I_ref = fits.getdata(results+"/"+image_file).astype(np.float64)
            P_ref = fits.getdata(results+"/"+pola_file).astype(np.float64) * P_unit
            flux_ref, major_ref, minor_ref, P_mean_ref = properties(I_ref, P_ref)
            print(code, "i =", incl, "scattered flux ratio:", flux/flux_ref, "polarisation:", P_mean, P_mean_ref)
            assert abs(flux/flux_ref - 1) < 0.10, f"1 um scattered flux differs from {code} at i = {incl} deg"
            if code == "TORUS" and n == 87: # cuts and polarisation not consistent with the other codes
                continue
            for name, cut, cut_ref in [("major", major, major_ref), ("minor", minor, minor_ref)]:
                mask = cut_ref > 1e-2 * cut_ref.max()
                ratio = cut[mask] / cut_ref[mask]
                print(code, "i =", incl, name, "axis cut ratio:", ratio.min(), ratio.max(), "median", np.median(ratio))
                assert abs(np.median(ratio) - 1) < 0.20, f"{name} axis cut differs from {code} at i = {incl} deg"
                if name == "major" and float(incl) <= 64:
                    assert np.all(abs(ratio - 1) < 0.15), f"major axis cut differs from {code} at i = {incl} deg"
            assert abs(P_mean/P_mean_ref - 1) < 0.30, f"polarisation degree differs from {code} at i = {incl} deg"


# Pinte et al. 2009 benchmark, temperatures with isotropic scattering (compared with RADMC,
# ProDiMo, TORUS and MCFOST 2009) and without scattering (RADICAL, ProDiMo, MCFOST 2009),
# tau = 1e3. Largest differences with MCFOST 4.1.16 (October 2026): medians 3 %, points
# 14 % (isotropic, ProDiMo at 200 au) and 18 % (no scattering, ProDiMo and RADICAL near
# the inner and outer edges, where their grids are coarser).
def _read_radmc_struct(filename):
    # RADMC / RADICAL structure file: spherical grid (r [cm], theta from the pole [rad]) and
    # temperature array (ir, itheta); returns the midplane temperature profile (r [au], T [K])
    blocks = {}
    key = None
    for line in open(filename):
        if line.startswith("#"):
            key = line.strip("# \n")
            blocks[key] = []
        else:
            blocks[key].extend(float(x) for x in line.split())
    nr, nt = (int(x) for x in blocks["NR,NT"])
    r = np.array(blocks["Spherical radius grid [cm] (middle of cell)"]) / 1.495978707e13
    theta = np.array(blocks["Theta grid [rad, from pole] (middle of cell)"])
    T = np.array(blocks["Temperature array (for ir=0,nr-1 do for it=0,nt-1 do ...)"]).reshape(nr, nt)
    return r, T[:, np.argmin(abs(theta - np.pi/2))]


@pytest.mark.slow
@pytest.mark.parametrize("case, opt", [("iso", "-isotropic"), ("no_scatt", "-no_scattering")])
def test_pinte2009_temperature(case, opt):
    if not compute_models:
        pytest.skip("Only run when computing models")
    utils = os.environ.get("MCFOST_UTILS", os.path.join("..", "utils"))
    results = os.path.join(utils, "benchmarks", "Pinte2009")
    if not os.path.isdir(results):
        pytest.skip("Pinte 2009 results not found in the MCFOST utils (mcfost -update_utils)")

    root_dir = "Pinte2009_tau=1e3_"+case
    clean_results(root_dir)
    run_mcfost("benchmarks/Pinte2009/Pinte2009_tau=1e3.para", root_dir, "-star_bb +disk_struct "+opt, timeout=3600)
    T = fits.getdata(root_dir+"/data_th/Temperature.fits.gz").astype(np.float64) # (z, r)
    grid = fits.getdata(root_dir+"/data_disk/grid.fits.gz")
    r = grid[0,0,0,:]

    def load(filename, skip=0):
        return np.loadtxt(results+"/"+filename, comments="#", skiprows=skip)

    prodimo = "prodimo/tau=1e3_"+case+"/T/ProDiMo_ltau3_"
    mcfost09 = "mcfost/tau=1e3_"+case+"/T/"
    midplane = {"ProDiMo": load(prodimo+"z0.0.dat")[:,[0,2]], "MCFOST 2009": load(mcfost09+"Radial_temp_tau=1e3.txt")}
    if case == "iso":
        midplane["RADMC"] = np.transpose(_read_radmc_struct(results+"/radmc/tau=1e3_iso/T/struct_3E-8.out"))
        midplane["TORUS"] = load("torus/tau=1e3_iso/T/tau=1e3_iso_radial.dat")
    else:
        midplane["RADICAL"] = np.transpose(_read_radmc_struct(results+"/radical/tau=1e3_no_scatt/T/struct_3E-8_noscat_radical.out"))
    # vertical cuts (z [au], T [K]); ProDiMo files have 2 header lines without "#"
    vertical = {(0.2, "ProDiMo"): load(prodimo+"r0.2.dat", 2)[:,[1,2]],
                (200., "ProDiMo"): load(prodimo+"r200.dat", 2)[:,[1,2]],
                (0.2, "MCFOST 2009"): load(mcfost09+"Vertical_temp_R=0.2AU_tau=1e3.txt"),
                (200., "MCFOST 2009"): load(mcfost09+"Vertical_temp_R=200AU_tau=1e3.txt")}

    for code, ref in midplane.items():
        ref = ref[(ref[:,0] >= r.min()) & (ref[:,0] <= r.max())]
        ratio = np.interp(np.log(ref[:,0]), np.log(r), T[0,:]) / ref[:,1]
        print(code, "midplane temperature ratio:", ratio.min(), ratio.max(), "median", np.median(ratio))
        assert abs(np.median(ratio) - 1) < 0.05, f"midplane temperature differs from {code}"
        assert np.all(abs(ratio - 1) < 0.25), f"midplane temperature differs from {code}"

    for (R, code), ref in vertical.items():
        i = np.searchsorted(r, R)
        f = np.log(R/r[i-1]) / np.log(r[i]/r[i-1])
        z = np.exp((1-f) * np.log(grid[1,0,:,i-1]) + f * np.log(grid[1,0,:,i]))
        T_R = (1-f) * T[:,i-1] + f * T[:,i]
        z_ref = abs(ref[:,0])
        mask = (z_ref >= z.min()) & (z_ref <= z.max())
        ratio = np.interp(np.log(z_ref[mask]), np.log(z), T_R) / ref[mask,1]
        print(code, "temperature ratio at R =", R, "au:", ratio.min(), ratio.max(), "median", np.median(ratio))
        assert abs(np.median(ratio) - 1) < 0.05, f"temperature at R = {R} au differs from {code}"
        assert np.all(abs(ratio - 1) < 0.25), f"temperature at R = {R} au differs from {code}"


def test_sublimation_image():
    # An image run after a temperature run with dust sublimation must read back
    # the sublimation radius and use the same inner radius as the temperature run
    if not compute_models:
        pytest.skip("Only run when computing models")
    # ref4.1_hydrostatic with hydrostatic equilibrium turned off: a single-zone 4.1 disk
    root_dir = "ref4.1_sublimation"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/sublimation.para"
    with open("test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para") as f_in, open(para, "w") as f_out:
        for line in f_in:
            if "sublimate dust" in line:
                line = line.replace("F", "T", 1)
            elif "hydostatic equilibrium" in line:
                line = line.replace("T", "F", 1)
            f_out.write(line)

    # -no_diff_approx: with the default grid, the first cells are optically thick once rmin moves inwards
    def run(opt):
        result = subprocess.run((_mcfost_bin+" "+para+" -no_diff_approx -root_dir "+root_dir+" "+opt).split(),
                                capture_output=True, text=True)
        print(result.stdout)
        assert result.returncode == 0, result.stdout[-2000:]
        return [l for l in result.stdout.splitlines() if "New minimum radius" in l]

    rmin_th = run("")
    rmin_img = run("-img 1.0")
    assert len(rmin_th) == 1 and len(rmin_img) == 1
    rmin_th = float(rmin_th[0].split()[-2])
    rmin_img = float(rmin_img[0].split()[-2])
    assert rmin_th > 0.0
    assert abs(rmin_img - rmin_th) < 1e-6 * rmin_th
    assert os.path.isfile(root_dir+"/data_1.0/RT.fits.gz")

    # The radial grid must follow the inner radius moved by sublimation: otherwise the
    # first radial cell spans from the sublimation radius to the original inner radius
    # and reaches unphysical temperatures (~2e4 K). Silicates sublimate at 1500 K here.
    T = fits.getdata(root_dir+"/data_th/Temperature.fits.gz")
    assert T.max() < 2000.0


def make_para(src, dst, edits):
    # Writes a copy of the parameter file src to dst. edits maps a text found on a
    # line (e.g. the comment of a parameter) to a function returning the new line.
    with open(src) as f_in, open(dst, "w") as f_out:
        for line in f_in:
            for key, edit in edits.items():
                if key in line:
                    line = edit(line)
            f_out.write(line)


def set_first_value(value):
    # replaces the first value of a parameter line, keeping the rest of the line
    return lambda line: "  "+value+" "+" ".join(line.split()[1:])+"\n"


def set_first_values(values):
    # replaces the first values of a parameter line (as many as given, separated by
    # spaces), keeping the rest of the line
    n = len(values.split())
    return lambda line: "  "+values+" "+" ".join(line.split()[n:])+"\n"


def run_mcfost(para, root_dir, opt="", cwd=None, threads=None, timeout=None):
    # cwd: directory to run from (para and root_dir are then relative to it), e.g. to
    # pick up dust files placed in the run directory, which is searched first
    # threads: number of OpenMP threads (default: environment)
    # timeout: maximum run time [s], the test fails if it is exceeded
    env = None
    if threads is not None:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    result = subprocess.run((os.path.abspath(_mcfost_bin)+" "+para+" -root_dir "+root_dir+" "+opt).split(),
                            capture_output=True, text=True, cwd=cwd, env=env, timeout=timeout)
    print(result.stdout)
    assert result.returncode == 0, result.stdout[-3000:]
    return result.stdout


def read_sublimation(root_dir):
    # Temperature (one value per cell) and sublimation mask (n_cells, n_grains or n_zones),
    # both in the internal cell order
    T = fits.getdata(root_dir+"/data_th/Temperature.fits.gz").ravel()
    mask = fits.getdata(root_dir+"/data_th/sublimation_mask.fits.gz") > 0.5
    assert mask.shape[0] == T.size
    return T, mask


def utils_T_sub(dust_file):
    # constant sublimation temperature of a dust file of the MCFOST utils (2nd value of
    # the first non-comment line), used with -constant_sublimation_temperature
    with open(os.path.join(os.environ["MCFOST_UTILS"], "Dust", dust_file)) as f:
        for line in f:
            if not line.lstrip().startswith("#"):
                return float(line.split()[1])


def flux_conservation(stdout):
    lines = [l for l in stdout.splitlines() if "Flux conservation" in l]
    assert len(lines) == 1
    return abs(float(lines[0].split()[4]))


def test_sublimation_front():
    # Single-zone 4.1 disk with its inner radius inside the sublimation radius: the
    # sublimation front is found by iterating the temperature calculation. At
    # convergence, the remaining dust is below T_sub and the sublimated cells (where a
    # trace of dust gives the test-particle temperature) are above it, within the
    # Monte Carlo hysteresis band. The constant sublimation temperature of the dust
    # file is used (the density-dependent law is tested in test_sublimation_power_law).
    if not compute_models:
        pytest.skip("Only run when computing models")
    T_sub = utils_T_sub("Draine_Si_sUV.dat")
    root_dir = "ref4.1_sublimation_front"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/front.para"
    make_para("test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para", para,
              {"sublimate dust": set_first_value("T"),
               "hydostatic equilibrium": set_first_value("F"),
               "Rin, edge, Rout": set_first_value("0.03")})

    stdout = run_mcfost(para, root_dir, "-no_diff_approx -constant_sublimation_temperature")
    assert "Dust sublimation converged" in stdout
    assert flux_conservation(stdout) < 2.0

    T, mask = read_sublimation(root_dir)
    sublimated = mask.any(axis=1)
    assert sublimated.sum() > 100
    assert T[~sublimated].max() < 1.1 * T_sub
    assert T[sublimated].min() > 0.9 * T_sub


def test_sublimation_multi_population():
    # Two zones and three populations (ref3.0_multi in LTE): silicates sublimate in
    # the inner zone, the ice-coated grains of the outer zone beyond the snow line, and
    # PAHs never sublimate. Each population sublimates independently (variable dust).
    # Only the innermost zone is extended to the sublimation radius (no overlapping
    # grid regions). The outer zone starts at 10 au, where the ice would not reach the
    # sublimation temperature of the utils (165 K): a copy of ice_opct.dat with
    # T_sub = 100 K is placed in the run directory, which is searched first.
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_multi_sublimation"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/multi.para"
    make_para("test_data/ref3.0_multi/ref3.0_multi.para", para,
              {"sublimate dust": set_first_value("T"),
               "Heating method": set_first_value("1"),
               "nbr_photons_eq_th": set_first_value("1.28e5"),
               # inner zone only (the outer zone starts at 10 au)
               "Rin, edge, Rout": lambda line: set_first_value("0.03")(line) if line.split()[0] == "1." else line})
    T_sub_ice = 100.0
    with open(os.path.join(os.environ["MCFOST_UTILS"], "Dust", "ice_opct.dat")) as f_in, \
         open(root_dir+"/ice_opct.dat", "w") as f_out:
        data_line = False
        for line in f_in:
            if not data_line and not line.lstrip().startswith("#"):
                line = line.split()[0]+" "+str(T_sub_ice)+"\n"
                data_line = True
            f_out.write(line)

    stdout = run_mcfost("multi.para", ".", "-no_diff_approx -no_SED -constant_sublimation_temperature",
                        cwd=root_dir)
    assert "Dust sublimation converged" in stdout

    T, mask = read_sublimation(root_dir)
    # grain sizes : 50 silicates (inner zone), 50 silicate + ice (outer zone), 1 PAH
    assert mask.shape[1] == 101
    silicates, ices, pah = mask[:,0:50], mask[:,50:100], mask[:,100]
    assert silicates.any()
    assert ices.any()
    assert not pah.any()
    assert T[silicates.any(axis=1)].min() > 0.7 * utils_T_sub("Draine_Si_sUV.dat")
    assert T[ices.any(axis=1)].min() > 0.7 * T_sub_ice


def test_sublimation_hydrostatic():
    # Hydrostatic equilibrium and sublimation iterated together: the gas structure,
    # dust redistribution and sublimation mask are updated in the same physics loop
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref4.1_sublimation_hydrostatic"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/hydro.para"
    make_para("test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para", para,
              {"sublimate dust": set_first_value("T"),
               "Rin, edge, Rout": set_first_value("0.03"),
               "nbr_photons_eq_th": set_first_value("1.28e5")})

    stdout = run_mcfost(para, root_dir, "-no_diff_approx -n_iter_physics 4")
    assert "Updating hydrostatic equilibrium" in stdout
    assert "Updating dust sublimation" in stdout
    assert os.path.isfile(root_dir+"/data_disk/gas_density.fits.gz")

    T, mask = read_sublimation(root_dir)
    assert np.isfinite(T).all()
    assert mask.any()

    # Image: the structure of the temperature run is restored in the same order (grid
    # extended to the sublimation radius, then hydrostatic density, then sublimation
    # mask), otherwise the grid redefinition overwrites the hydrostatic density
    stdout = run_mcfost(para, root_dir, "-img 1.0")
    lines = stdout.splitlines()
    def first_line(text):
        return next(i for i, l in enumerate(lines) if text in l)
    assert first_line("Reading sublimation file") < first_line("Restoring hydrostatic density") \
        < first_line("Reading sublimation mask")
    assert os.path.isfile(root_dir+"/data_1.0/RT.fits.gz")


def test_sublimation_sed_only():
    # SED computed without the temperature calculation (-no_T), from the temperature,
    # sublimation radius and sublimation mask of a previous run: same SED as the one
    # computed after the temperature calculation (without the mask, the sublimated
    # cells keep their dust at the test-particle temperature: ~2.5 x the 1-5 micron flux)
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref4.1_sublimation_sed"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/sed.para"
    make_para("test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para", para,
              {"sublimate dust": set_first_value("T"),
               "hydostatic equilibrium": set_first_value("F"),
               "Rin, edge, Rout": set_first_value("0.03")})
    opt = "-no_diff_approx -constant_sublimation_temperature"

    run_mcfost(para, root_dir, opt)
    sed_full = fits.getdata(root_dir+"/data_th/sed_rt.fits.gz")
    # the SED files are not overwritten
    for f in glob.glob(root_dir+"/data_th/sed_*.fits.gz"):
        os.remove(f)

    stdout = run_mcfost(para, root_dir, opt+" -no_T")
    assert "Reading sublimation mask" in stdout
    sed = fits.getdata(root_dir+"/data_th/sed_rt.fits.gz")
    # total flux (first contribution), all inclinations and wavelengths
    ratio = sed[0] / np.where(sed_full[0] > 0, sed_full[0], 1.0)
    assert abs(np.median(ratio) - 1.0) < 0.02
    assert ratio.max() < 1.3


def test_sublimation_per_grain_size():
    # nLTE grains have one temperature per grain size: the small (hotter) ice-coated
    # grains of ref3.0_multi sublimate over a larger region than the large ones, and
    # some cells are only partially sublimated. PAHs (nRE) never sublimate.
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref3.0_multi_sublimation_nlte"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/multi_nlte.para"
    make_para("test_data/ref3.0_multi/ref3.0_multi.para", para,
              {"sublimate dust": set_first_value("T"),
               "nbr_photons_eq_th": set_first_value("1.28e5"),
               "Rin, edge, Rout": lambda line: set_first_value("0.03")(line) if line.split()[0] == "1." else line})

    stdout = run_mcfost(para, root_dir, "-no_diff_approx -no_SED")
    assert "Dust sublimation converged" in stdout

    T, mask = read_sublimation(root_dir)
    # grain sizes : 50 silicates (LTE), 50 silicate + ice (nLTE, T_sub = 100 K), 1 PAH (nRE)
    ices, pah = mask[:,50:100], mask[:,100]
    n_cells_per_size = ices.sum(axis=0)
    assert n_cells_per_size[0] > 0
    assert n_cells_per_size[0] > n_cells_per_size[-1]
    assert (ices.any(axis=1) & ~ices.all(axis=1)).sum() > 0
    assert not pah.any()


def test_sublimation_power_law():
    # Density-dependent sublimation temperature, read from a .sub file next to the
    # optical indices file: T_sub = T0 (rho_gas / rho0)^gamma (Pollack et al. 1994).
    # The .sub file is placed in the run directory (searched before the utils dust
    # directory) so that it does not affect the other tests. At convergence, each cell
    # is on the correct side of its own T_sub(rho), within the hysteresis band, and the
    # sublimation temperature varies along the curved front (dense midplane hotter).
    # -constant_sublimation_temperature ignores the .sub file.
    if not compute_models:
        pytest.skip("Only run when computing models")
    T0, rho0, gamma = 2000.0, 1.0, 0.0195
    root_dir = "ref4.1_sublimation_power_law"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    make_para("test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para", root_dir+"/power_law.para",
              {"sublimate dust": set_first_value("T"),
               "hydostatic equilibrium": set_first_value("F"),
               "Rin, edge, Rout": set_first_value("0.03")})
    with open(root_dir+"/Draine_Si_sUV.sub", "w") as f:
        f.write("# test sublimation law\n"
                "format_version = 1\n"
                "law = power_law\n"
                f"T0 = {T0}\n"
                f"rho0 = {rho0}\n"
                f"gamma = {gamma}\n")

    # the constant law is used when forced, even with a .sub file
    stdout = run_mcfost("power_law.para", "constant", "-constant_sublimation_temperature -disk_struct", cwd=root_dir)
    assert "Sublimation law from" not in stdout

    # +disk_struct writes the gas density on the grid moved to the sublimation radius
    stdout = run_mcfost("power_law.para", ".", "-no_diff_approx -no_SED +disk_struct", cwd=root_dir)
    assert "Sublimation law from Draine_Si_sUV.sub" in stdout
    assert "Dust sublimation converged" in stdout

    T, mask = read_sublimation(root_dir)
    sublimated = mask.any(axis=1)
    rho = fits.getdata(root_dir+"/data_disk/gas_density.fits.gz").ravel() # g/cm3
    assert rho.size == T.size
    T_sub = T0 * (np.maximum(rho, 1e-30) / rho0)**gamma
    assert sublimated.sum() > 100
    assert (T / T_sub)[~sublimated].max() < 1.1
    assert (T / T_sub)[sublimated].min() > 0.9
    assert T_sub[sublimated].max() > 1.3 * T_sub[sublimated].min()


def nre_absorbed_fractions(stdout):
    # fraction of the energy absorbed by non-equilibrium grains, first nRE pass of each
    # temperature calculation
    values = [float(l.split()[-4]) for l in stdout.splitlines() if "Non equilibrium grains have absorbed" in l]
    return values[0::2]


def check_hydrostatic_temperature_stable(root_dir, tolerance=0.05):
    T1 = fits.getdata(root_dir+"/data_th/Temperature_iter0001.fits.gz")
    T2 = fits.getdata(root_dir+"/data_th/Temperature_iter0002.fits.gz")
    mask = T1 > 2.0
    assert abs(np.median(T2[mask]/T1[mask]) - 1.0) < tolerance


def test_hydrostatic_nre():
    # Hydrostatic iterations with non-equilibrium grains (PAHs): the packet budget,
    # emissivities and equilibrium status modified during a temperature calculation are
    # restored before the next one, and the split between equilibrium and
    # non-equilibrium absorption is rebuilt with the new density. The energy absorbed
    # by the PAHs and the temperature must stay consistent between iterations.
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref4.1_PAH_hydrostatic"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/pah_hydro.para"
    make_para("test_data/ref4.1_PAH/ref4.1_PAH.para", para,
              {"hydostatic equilibrium": set_first_value("T")})

    stdout = run_mcfost(para, root_dir, "-n_iter_physics 3")
    fractions = nre_absorbed_fractions(stdout)
    assert len(fractions) == 3
    for f in fractions[1:]:
        assert abs(f/fractions[0] - 1.0) < 0.3
    check_hydrostatic_temperature_stable(root_dir)


def test_hydrostatic_nlte():
    # Hydrostatic iterations with nLTE grains (one temperature per grain size)
    if not compute_models:
        pytest.skip("Only run when computing models")
    root_dir = "ref4.1_nLTE_hydrostatic"
    clean_results(root_dir)
    os.makedirs(root_dir, exist_ok=True)
    para = root_dir+"/nlte_hydro.para"
    make_para("test_data/ref4.1_nLTE/ref4.1_nLTE.para", para,
              {"hydostatic equilibrium": set_first_value("T")})

    run_mcfost(para, root_dir, "-n_iter_physics 3")
    check_hydrostatic_temperature_stable(root_dir)
