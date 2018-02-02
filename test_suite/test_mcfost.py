import astropy.io.fits as fits
import numpy as np
import subprocess
import shutil
import os

_mcfost_bin = "../src/mcfost"

def mcfost(filename,opt=""):
    cmd = _mcfost_bin+" "+filename+" "+opt
    result = subprocess.call(cmd.split())

def clean_results():
    result = subprocess.call("rm","-rf data_*")

def all_almost_equal(x,y,threshold=0.01):
    # test if all the values of two arrays are almost equal
    return (abs((x-y)) < threshold * x).all()

def MC_similar(x,y,threshold=0.01):
    # test if two arrays have the same 75% percentile
    # ignoring values that are very small as they are very noisy

    mask = (abs(x) < 1e-300) | (abs(y) < 1e-300)
    x_ma = np.ma.masked_where(mask, x)
    y_ma = np.ma.masked_where(mask, y)

    #return (abs((x_ma-y_ma)/x_ma).mean() < threshold)
    return ( np.percentile(abs((x_ma-y_ma)/x_ma), 75)  < threshold )


def test_Temperature():
    # Run the mcfost model
    filename = "test_data/ref3.0/ref3.0.para"
    if (os.path.isdir("data_th")):
        shutil.rmtree("data_th")
    mcfost(filename,opt="-mol")

    # Read the results
    T = fits.getdata("data_th/Temperature.fits.gz")
    T_ref = fits.getdata("test_data/ref3.0/data_th/Temperature.fits.gz")

    print("Maximum T difference", (abs(T-T_ref)/(T_ref+1e-30)).max())
    print("Mean T difference   ", (abs(T-T_ref)/(T_ref+1e-30)).mean())

    assert MC_similar(T_ref,T,threshold=0.05)

def test_SED():
    # Re-use the previous mcfost model

    # Read the results
    SED = fits.getdata("data_th/sed_rt.fits.gz")
    SED_ref = fits.getdata("test_data/ref3.0/data_th/sed_rt.fits.gz")

    print("Maximum SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=0.05)


def test_mol_map():
    # Re-use the previous mcfost model

    # Read the results
    image = fits.getdata("data_CO/lines.fits.gz")
    image_ref = fits.getdata("test_data/ref3.0/data_CO/lines.fits.gz")

    print("Maximum mol map difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean mol map difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)

def test_image():
    # Run the mcfost model
    filename = "test_data/ref3.0/ref3.0.para"
    if (os.path.isdir("data_1.0")):
        shutil.rmtree("data_1.0")
    mcfost(filename,opt="-img 1.0")

    # Read the results
    image = fits.getdata("data_1.0/RT.fits.gz")
    image_ref = fits.getdata("test_data/ref3.0/data_1.0/RT.fits.gz")

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)
