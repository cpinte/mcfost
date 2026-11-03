import os
import re
import shutil
import subprocess
import numpy as np
import astropy.io.fits as fits
import pytest

_test_dir = os.path.dirname(os.path.abspath(__file__))
_mcfost_bin = os.path.join(_test_dir, "..", "src", "mcfost")
_ref3_para = os.path.join(_test_dir, "test_data", "ref3.0", "ref3.0.para")


def MC_similar(x, y, threshold=0.03, mask_threshold=1e-24):
    mask = abs(x) < mask_threshold
    x_ma = np.ma.masked_where(mask, x)
    y_ma = np.ma.masked_where(mask, y)
    p75 = np.percentile(abs((x_ma - y_ma) / x_ma).compressed(), 75)
    return p75 < threshold


def run_mcfost_cmd(args, timeout=None, env=None):
    # No time limit by default, as run_mcfost in test_mcfost.py: the ref3.0 temperature runs
    # take more than 120 s on the macOS Intel runners
    cmd = [_mcfost_bin] + args
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    return result


def test_mrw_cli_flags(tmp_path):
    """Verify MRW command-line flags and input validation."""
    assert os.path.isfile(_mcfost_bin), f"mcfost binary not found at {_mcfost_bin}"

    # -help should list -mrw options
    res = run_mcfost_cmd(["-help"])
    assert "-mrw" in res.stdout
    assert "-no_mrw" in res.stdout
    assert "-gamma_mrw" in res.stdout

    # Test invalid gamma < 2.0 exits with error
    res_err = run_mcfost_cmd([_ref3_para, "-gamma_mrw", "1.5"])
    assert res_err.returncode != 0
    assert "gamma_MRW must be" in res_err.stdout or "gamma_MRW must be" in res_err.stderr

    # Test valid custom gamma via -mrw 3.0 with fast -op mode
    dir_test = str(tmp_path / "test_cli")
    res_gamma = run_mcfost_cmd([_ref3_para, "-root_dir", dir_test, "-op", "0.81", "-mrw", "3.0"])
    assert res_gamma.returncode == 0
    # list-directed output: the number of digits depends on the compiler (gfortran, ifort)
    m = re.search(r"MRW enabled with gamma =\s*([0-9.Ee+-]+)", res_gamma.stdout)
    assert m is not None and float(m.group(1)) == 3.0


def test_mrw_temperature_and_sed(tmp_path):
    """Verify that MRW produces temperatures and thermal SEDs consistent with standard MC."""
    dir_std = str(tmp_path / "run_std")
    dir_mrw = str(tmp_path / "run_mrw")

    # Run Standard MC
    res_std = run_mcfost_cmd([_ref3_para, "-root_dir", dir_std, "-no_mrw"])
    assert res_std.returncode == 0, f"Standard MC failed:\n{res_std.stdout[-2000:]}"

    # Run MRW
    res_mrw = run_mcfost_cmd([_ref3_para, "-root_dir", dir_mrw, "-mrw"])
    assert res_mrw.returncode == 0, f"MRW run failed:\n{res_mrw.stdout[-2000:]}"
    assert "MRW enabled" in res_mrw.stdout

    # Compare Temperature
    T_std = fits.getdata(os.path.join(dir_std, "data_th", "Temperature.fits.gz")).astype(np.float64)
    T_mrw = fits.getdata(os.path.join(dir_mrw, "data_th", "Temperature.fits.gz")).astype(np.float64)

    assert T_std.shape == T_mrw.shape
    assert MC_similar(T_std, T_mrw, threshold=0.03), "MRW temperature differs from Standard MC by > 3%"

    med_rel_diff = np.median(abs(T_mrw - T_std) / np.maximum(T_std, 2.7)) * 100
    print(f"MRW vs Standard MC median T difference: {med_rel_diff:.3f}%")
    assert med_rel_diff < 0.5, f"MRW median T difference {med_rel_diff:.3f}% exceeds 0.5%"

    # Compare ray-traced thermal SED (component 3, directly produced by the computed temperature field)
    sed_file_std = os.path.join(dir_std, "data_th", "sed_rt.fits.gz")
    sed_file_mrw = os.path.join(dir_mrw, "data_th", "sed_rt.fits.gz")

    if os.path.exists(sed_file_std) and os.path.exists(sed_file_mrw):
        sed_std = fits.getdata(sed_file_std).astype(np.float64)
        sed_mrw = fits.getdata(sed_file_mrw).astype(np.float64)

        # Thermal SED component (index 3)
        thermal_std = sed_std[3, 0, :, :]
        thermal_mrw = sed_mrw[3, 0, :, :]

        for i_inc in range(thermal_std.shape[0]):
            mask = thermal_std[i_inc, :] > 1e-4 * thermal_std[i_inc, :].max()
            assert MC_similar(thermal_std[i_inc, mask], thermal_mrw[i_inc, mask], threshold=0.02), (
                f"Thermal SED differs by > 2% at inclination index {i_inc}"
            )

        # Total SED component (index 0): verify median difference across all inclinations is < 3%
        total_std = sed_std[0, 0, :, :]
        total_mrw = sed_mrw[0, 0, :, :]
        mask_tot = total_std > 1e-4 * total_std.max()
        med_sed_diff = np.median(abs(total_mrw[mask_tot] - total_std[mask_tot]) / total_std[mask_tot]) * 100
        print(f"MRW vs Standard MC median Total SED difference: {med_sed_diff:.3f}%")
        assert med_sed_diff < 3.0, f"Total SED median difference {med_sed_diff:.3f}% exceeds 3%"


@pytest.mark.slow
def test_mrw_optically_thick_benchmark(tmp_path):
    """Verify MRW speedup and accuracy on Pinte et al. 2009 tau=1e4 optically thick disk."""
    bench_para = os.path.join(_test_dir, "benchmarks", "Pinte2009", "Pinte2009_tau=1e4.para")
    if not os.path.exists(bench_para):
        pytest.skip("Pinte2009_tau=1e4.para not found")

    # Read and configure 1.28e5 photons for temperature only
    test_para_path = str(tmp_path / "pinte_1e4_test.para")
    with open(bench_para) as f:
        lines = f.readlines()
    lines[3] = "  1.28e5                  nbr_photons_eq_th  : T computation\n"
    lines[4] = "  0                       nbr_photons_lambda : SED computation\n"
    lines[9] = "  T F F                   compute temperature?, compute sed?, use default wavelength grid for output ?\n"
    with open(test_para_path, "w") as f:
        f.writelines(lines)

    dir_std = str(tmp_path / "run_std_1e4")
    dir_mrw = str(tmp_path / "run_mrw_1e4")

    # Fixed seed and 1 thread: reproducible, the check tests the MRW bias and not the noise
    env = dict(os.environ, OMP_NUM_THREADS="1")
    res_std = run_mcfost_cmd([test_para_path, "-star_bb", "-root_dir", dir_std, "-no_mrw", "-seed", "1"], env=env)
    assert res_std.returncode == 0

    res_mrw = run_mcfost_cmd([test_para_path, "-star_bb", "-root_dir", dir_mrw, "-mrw", "-seed", "1"], env=env)
    assert res_mrw.returncode == 0

    T_std = fits.getdata(os.path.join(dir_std, "seed=1", "data_th", "Temperature.fits.gz")).astype(np.float64)
    T_mrw = fits.getdata(os.path.join(dir_mrw, "seed=1", "data_th", "Temperature.fits.gz")).astype(np.float64)

    # Midplane check
    valid = T_std[0, :] > 2.7
    rel_diff_mid = abs(T_mrw[0, valid] - T_std[0, valid]) / T_std[0, valid]
    assert np.median(rel_diff_mid) < 0.05, f"Midplane median difference {np.median(rel_diff_mid):.3f} >= 5%"
