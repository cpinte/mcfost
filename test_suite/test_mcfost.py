import astropy.io.fits as fits
import subprocess
import shutil
import os

def mcfost(filename,opt=""):
    args = filename+" "+opt
    result = subprocess.call(["../src/mcfost",args])

def clean_results():
    result = subprocess.call("rm","-rf data_*")

def all_almost_equal(x,y,threshold=0.01):
    # test if all the values of two arrays are almost equal
    return (abs((x-y)) < threshold * x).all()

def similar(x,y,threshold=0.01):
    # test if two arrays have the same average
    return (abs((x-y)/x).mean() < threshold)


def test_SED():
    # Run the mcfost model
    filename = "test_data/ref3.0/ref3.0.para"
    if (os.path.isdir("data_th")):
        shutil.rmtree("data_th")
    mcfost(filename)

    # Read the results
    T = fits.getdata("data_th/Temperature.fits.gz")
    T_ref = fits.getdata("test_data/ref3.0/data_th/Temperature.fits.gz")

    print("Maximum T difference", (abs(T-T_ref)/T_ref).max())
    print("Mean T difference   ", (abs(T-T_ref)/T_ref).mean())

    assert similar(T_ref,T,threshold=0.05)
