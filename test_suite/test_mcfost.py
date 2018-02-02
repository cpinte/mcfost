import astropy.io.fits as fits
import numpy as np
import subprocess
import pytest
import shutil
import glob
import os

# Reference tests are computed with v3.0.28

_mcfost_bin = "../src/mcfost"

model_list = glob.glob1("test_data/","*")
wl_list = ["1.0","10","100","1000"]
wl_list_pola = ["1.0","1000"]
wl_list_contrib = ["1.0","100","1000"]


def mcfost(filename,opt=""):
    cmd = _mcfost_bin+" "+filename+" "+opt
    result = subprocess.call(cmd.split())

def clean_results():
    cmd = ["rm","-rf"]+glob.glob("data_*")
    result = subprocess.call(cmd)

def all_almost_equal(x,y,threshold=0.01):
    # test if all the values of two arrays are almost equal
    return (abs((x-y)) < threshold * x).all()

def MC_similar(x,y,threshold=0.01,mask_threshold=1e-20):
    # test if two arrays have the same at the 75% percentile
    # ignoring values that are very small as they are very noisy

    mask = (abs(x) < mask_threshold) | (abs(y) < mask_threshold) # Unit: W.m-2 or W.m-2/pixel
    x_ma = np.ma.masked_where(mask, x)
    y_ma = np.ma.masked_where(mask, y)

    #return (abs((x_ma-y_ma)/x_ma).mean() < threshold)
    return ( np.percentile(abs((x_ma-y_ma)/x_ma), 75)  < threshold )

def test_mcfost_bin():
    # We first test if the mcfost binary actually exists
    return os.path.isfile(_mcfost_bin)

@pytest.mark.parametrize("model_name", model_list)
def test_Temperature(model_name):
    # Run the mcfost model
    filename = "test_data/"+model_name+"/"+model_name+".para"

    clean_results() # removing all previous calculations
    mcfost(filename,opt="-mol")

    # Read the results
    T = fits.getdata("data_th/Temperature.fits.gz")
    T_ref = fits.getdata("test_data/ref3.0/data_th/Temperature.fits.gz")

    print("Maximum T difference", (abs(T-T_ref)/(T_ref+1e-30)).max())
    print("Mean T difference   ", (abs(T-T_ref)/(T_ref+1e-30)).mean())

    assert MC_similar(T_ref,T,threshold=0.05)

@pytest.mark.parametrize("model_name", model_list)
def test_SED(model_name):
    # Re-use the previous mcfost model

    # Read the results
    SED = fits.getdata("data_th/sed_rt.fits.gz")
    SED_ref = fits.getdata("test_data/"+model_name+"/data_th/sed_rt.fits.gz")

    print("Maximum SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=0.05)


@pytest.mark.parametrize("model_name", model_list)
def test_mol_map(model_name):
    # Re-use the previous mcfost model

    # Read the results
    image = fits.getdata("data_CO/lines.fits.gz")
    image_ref = fits.getdata("test_data/"+model_name+"/data_CO/lines.fits.gz")

    print("Maximum mol map difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean mol map difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list)
def test_image(model_name, wl):
    # Run the mcfost model
    filename = "test_data/"+model_name+"/"+model_name+".para"
    mcfost(filename,opt="-img "+wl)

    # Read the results
    image = fits.getdata("data_"+wl+"/RT.fits.gz")
    image_ref = fits.getdata("test_data/"+model_name+"/data_"+wl+"/RT.fits.gz")

    # We just keep intensity
    image = image[0,:,:,:,:]
    image_ref = image_ref[0,:,:,:,:]

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list_pola)
def test_pola(model_name, wl):
    # Re-use previous calculation

    # Read the results
    image = fits.getdata("data_"+wl+"/RT.fits.gz")
    image_ref = fits.getdata("test_data/"+model_name+"/data_"+wl+"/RT.fits.gz")

    # We just keep Stokes Q, U
    image = image[[1,2],:,:,:,:]
    image_ref = image_ref[[1,2],:,:,:,:]

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)

@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list_contrib)
def test_contrib(model_name, wl):
    # Re-use previous calculation

    # Read the results
    image = fits.getdata("data_"+wl+"/RT.fits.gz")
    image_ref = fits.getdata("test_data/"+model_name+"/data_"+wl+"/RT.fits.gz")

    # We just keep separate contributions
    image = image[[4,5,6,7],:,:,:,:]
    image_ref = image_ref[[4,5,6,7],:,:,:,:]

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)
