import astropy.io.fits as fits
import numpy as np
import subprocess
import pytest
import shutil
import glob
import os

# Reference tests are computed with db94244c41cad9c1517363b1607f8596af1c55c0
_mcfost_bin = "../src/mcfost"

# Get list of models using directory names
model_list = glob.glob1("test_data/","*")

# If running on CI, only test ref3.0
if os.environ.get('CI', None) == 'true':
    model_list = ["ref3.0"]

wl_list = ["1.0","10","100","1000"]
wl_list_pola = ["1.0","1000"]
wl_list_contrib = ["1.0","100","1000"]


def mcfost(filename,opt=""):
    cmd = _mcfost_bin+" "+filename+" "+opt
    result = subprocess.call(cmd.split())

def clean_results(model_name):
    #cmd = ["rm","-rf"]+glob.glob("data_*")
    cmd = ["rm","-rf",model_name]
    result = subprocess.call(cmd)

def all_almost_equal(x,y,threshold=0.01):
    # test if all the values of two arrays are almost equal
    return (abs((x-y)) < threshold * x).all()

def MC_similar(x,y,threshold=0.01,mask_threshold=1e-25):
    # test if two arrays have the same at the 75% percentile
    # ignoring values that are very small as they are very noisy

    #mask = (abs(x) < mask_threshold) | (abs(y) < mask_threshold) # Unit: W.m-2 or W.m-2/pixel
    mask = abs(x) < mask_threshold # Unit: W.m-2 or W.m-2/pixel : only using reference for the mask
    x_ma = np.ma.masked_where(mask, x)
    y_ma = np.ma.masked_where(mask, y)

    #return (abs((x_ma-y_ma)/x_ma).mean() < threshold)
    return ( np.percentile(abs((x_ma-y_ma)/x_ma), 75)  < threshold )

def test_mcfost_bin():
    # We first test if the mcfost binary actually exists and runs
    assert os.path.isfile(_mcfost_bin)
    try:
        subprocess.call(_mcfost_bin)
    except OSError as e:
        assert not e.errno == os.errno.ENOENT

@pytest.mark.parametrize("model_name", model_list)
def test_Temperature(model_name):
    clean_results(model_name) # removing all previous calculations

    # Run the mcfost model
    filename = "test_data/"+model_name+"/"+model_name+".para"
    if (model_name == "discF_00500"):
        opt=" -phantom test_data/"+model_name+"/"+model_name
    else:
        opt=""
    mcfost(filename,opt="-mol -root_dir "+model_name+opt)

    # Read the results
    T_name = model_name+"/data_th/Temperature.fits.gz"
    T = fits.getdata(T_name)
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
        pytest.skip("No SED")
    SED = fits.getdata(SED_name)
    SED_ref = fits.getdata("test_data/"+SED_name)

    print("Maximum SED difference", (abs(SED-SED_ref)/(SED_ref+1e-30)).max())
    print("Mean SED difference   ", (abs(SED-SED_ref)/(SED_ref+1e-30)).mean())

    assert MC_similar(SED_ref,SED,threshold=0.1)


@pytest.mark.parametrize("model_name", model_list)
def test_mol_map(model_name):
    # Re-use the previous mcfost model

    # Read the results
    image_name = model_name+"/data_CO/lines.fits.gz"
    image = fits.getdata(image_name)
    image_ref = fits.getdata("test_data/"+image_name)

    print("Maximum mol map difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean mol map difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list)
def test_image(model_name, wl):
    # Run the mcfost model
    filename = "test_data/"+model_name+"/"+model_name+".para"
    if (model_name == "discF_00500"):
        opt=" -phantom test_data/"+model_name+"/"+model_name
    else:
        opt=""
    mcfost(filename,opt="-img "+wl+" -root_dir "+model_name+opt)

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"
    image = fits.getdata(image_name)
    image_ref = fits.getdata("test_data/"+image_name)

    # We just keep intensity
    image = image[0,:,:,:,:]
    image_ref = image_ref[0,:,:,:,:]

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.1)


@pytest.mark.parametrize("model_name", model_list)
@pytest.mark.parametrize("wl", wl_list_pola)
def test_pola(model_name, wl):
    # Re-use previous calculation

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"
    image = fits.getdata(image_name)
    if ((image.shape[0] != 4) & (image.shape[0] != 8)):
        pytest.skip("No pola")
    image_ref = fits.getdata("test_data/"+image_name)

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

    # Skip test on CI because it currently fails
    # TODO: fix this test failure
    if os.environ.get('CI', None) == 'true':
        pytest.skip("CI")

    # Read the results
    image_name = model_name+"/data_"+wl+"/RT.fits.gz"
    image = fits.getdata(image_name)
    if ((image.shape[0] != 5) & (image.shape[0] != 8)):
        pytest.skip("No contrib")
    image_ref = fits.getdata("test_data/"+image_name)

    # We just keep separate contributions
    image = image[[4,5,6,7],:,:,:,:]
    image_ref = image_ref[[4,5,6,7],:,:,:,:]

    print("Maximum image difference", (abs(image-image_ref)/(image_ref+1e-30)).max())
    print("Mean image difference   ", (abs(image-image_ref)/(image_ref+1e-30)).mean())

    assert MC_similar(image_ref,image,threshold=0.05)
