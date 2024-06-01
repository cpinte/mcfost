import astropy.io.fits as fits
import numpy as np
import subprocess
import pytest
import shutil
import glob
import os

_mcfost_bin = "../src/mcfost"

# Get list of models using directory names
model_list = glob.glob1("test_data/","*")

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

    # We just keep intensity
    for i in range(n_incl):
        image = np.nan_to_num(fits.getdata(test_dir+"/"+image_name))
        image_ref = fits.getdata("test_data/"+image_name)

        print("-- i=", i)
        print(image.shape)
        print(image_ref.shape)

        image = image[0,:,i,:,:]
        image_ref = image_ref[0,:,i,:,:]

        print("Min ref", image_ref.min(), image_ref.max())
        print("Min    ", image.min(), image.max())

        print("Maximum image difference", i, (abs(image-image_ref)/(image_ref+1e-30)).max())
        print("Mean image difference   ", i, (abs(image-image_ref)/(image_ref+1e-30)).mean())

        print("MC similar", MC_similar(image_ref,image,threshold=0.1))
        print("-------------")

    image = fits.getdata(test_dir+"/"+image_name)
    image_ref = fits.getdata("test_data/"+image_name)

    image = image[0,:,:,:,:]
    image_ref = image_ref[0,:,:,:,:]

    threshold=0.1
    if (model_name == "ref3.0") and (wl == "100"): # weird difference on linux, ifort, openmp=no, release=no
        threshold=0.11

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
    else:
        mask_threshold = 1e-21

    assert MC_similar(image_ref,image,threshold=0.1,mask_threshold=mask_threshold)

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
    else:
        mask_threshold=1e-23

    threshold=0.1
    if (model_name == "ref3.0") and (wl == "100"): # weird difference on linux, ifort, openmp=no, release=no
        threshold=0.11

    assert MC_similar(image_ref,image,threshold=threshold,mask_threshold=mask_threshold)
