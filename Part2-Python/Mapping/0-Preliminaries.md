---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.6.0
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Datasets

We need to download some data for the examples. Some files are bundled with the notebooks because I am not smart enough to understand how to download them through the various web / ftp interfaces directly or I needed to filter or compress the information for the purposes of the class. Please go to the original sites if you need to use these data for anything other than the demonstrations in these notebooks.

Most datasets are too large for the repository (and generally, that's not a place to keep anything which is not the primary target of revision management. The bundled data are compressed and will be unpacked (copied) to the `../../Data/Resources` directory by the "download" functions in this notebook. If you mess them up, just run the download again. Anything that is undamaged will just be be skipped anyway.

## Global Magnetic Data

Magnetic intensity data from [geomag.org](http://geomag.org/models/EMAG2/EMAG2_V2.tif)

## Topography data

ETOPO1 images are from NOAA - we use their geotifs which are subsampled from the original (enormous) dataset but 

## NASA Blue marble images

Winter and Summer images for the Earth are grabbed for plotting examples. The winter ones (June) are used by default as these have less ice in the N. Hemisphere. 

## Earthquake hypocentres

Are grabbed from the [NEIC](http://earthquake.usgs.gov/earthquakes/search/) - they are in the geoJSON format since that is very simple to read with python. The downloads are limited at 20k events so the time and magnitude range is whatever it takes to get just under this limit. The filenames give clues about that, but, so does the catalogue itself once it is read in.

## Global age grid 

Taken from Earthbyte and reduced in size by throwing away the grid information and saving in compressed numpy format. 

## Global strain rate

I grabbed the second invariant of the strain rate from the [global strain rate map](http://gsrm.unavco.org/intro) project through the 'Jules Vernes' portal. 


```{code-cell} ipython3
:hide_input: false

%%sh

ls -l ../../Data/Reference/
ls -l ../../Data/Resources/

cp ../../Data/Reference/* ../../Data/Resources/
```

+++ {"solution": "shown", "solution_first": true}

## The datasets that we use for this course are kept in a central location. 

These can be downloaded on demand but the smaller files are stored away in the `Reference` directory and copied to the `Resources` directory. 

   - Why do I do that ? 

So if you happen to delete or break one of the datasets, we can get a new copy. __Expand the cell__ to see the details of the download and caching mechanism. 

   - `md5` checksums
   - the `requests` module for handling urls (properly)
   - exception handling (!)

``` python

import requests
import os
import time
import sys
import shutil

def download_file(url, local_filename, expected_size):
    
    if (os.path.isfile(url)):
        shutil.copy(url, local_filename)
        
# and so on ...         
``` 

```{code-cell} ipython3
:hide_input: false
:solution: shown

## Here are some functions to download or install files. 


import requests
import os
import time
import sys
import shutil

def download_file(url, local_filename, expected_size):
    
    
# We might want to bundle some files if they are small / compressed or not readily available for download

    if (os.path.isfile(url)):
        shutil.copy(url, local_filename)
    
    else:
        r = requests.get(url, stream=True)

        start_time = time.time()
        last_time = start_time
        datasize = 0

        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=10000000): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
                    datasize += len(chunk)
                    if (time.time() - last_time) > 2.5:
                        print ("{:.1f} Mb in {:.1f}s / {}".format(datasize / 1.0e6, time.time() - start_time, expected_size))
                        last_time = time.time()

    return 
```

```{code-cell} ipython3
:hide_input: false
:solution: shown

import hashlib
from functools import partial

def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 4096), b''):
            d.update(buf)
    return d.hexdigest()




def download_cached_file(location_url, local_file, expected_md5, expected_size="Unknown"):

    try:
        assert md5sum(local_file) == expected_md5
        print ("Using cached file - {}".format(local_file))
        return(2)
    
    except (IOError, AssertionError) as error_info:
        # No file or the wrong file ... best go download it
        # print "Assertion failed - {}".format(sys.exc_info())
        
        try:
            data_file = download_file(location_url, local_file, expected_size)
            print ("Downloaded from {}".format( location_url ))
            return(1)

        except:
            print ("Unable to download {} [{}] ".format( location_url, sys.exc_info() ))
            return(0)


def report_cached_file(local_file, expected_md5):
    
    import os.path
    
    if not os.path.isfile(local_file):
        print ("Local file {} does not exist".format(local_file))
    else:
        if len(expected_md5) == 0 or md5sum(local_file) == expected_md5:
            print ("Cached file {} is valid".format(local_file))
        else:
            print ("Cached file {} failed, checksum {}".format(local_file, md5sum(local_file)))
        
        
```

+++ {"solution": "hidden", "solution_first": true}

## The resources "database"

The files are managed by keeping a table of the download point and a checksum to validate the file.

Although it is a little laborious to set this up, it is a useful way to keep track of the data sources. For very large collections, you would probably do better with a database not a hand-built dictionary

The database takes the form of a dictionary like this:

```python

resource_list = [
# Global magnetic data 
{
 "local_file":"../../Data/Resources/EMAG2_image_V2.tif",
 "md5":"c4944c27ed22ddc89225e506936b016b",
 "url":"http://geomag.org/models/EMAG2/EMAG2_image_V2.tif",
 "expected_size":"133Mb"
},

# and so on ... 

```

```{code-cell} ipython3
:hide_input: false
:solution: hidden

# ../../Data/Resources required for the tutorial
# The checksum is to make sure any cached file is the right one.
# If the checksum is wrong the file gets downloaded again.

# Some of these files are local, stored in the repo in compressed form and installed into the
# ../../Data/Resources directory. 

resource_list = [
# Global magnetic data 
{
 "local_file":"../../Data/Resources/EMAG2_image_V2.tif",
 "md5":"c4944c27ed22ddc89225e506936b016b",
 "url":"http://geomag.org/models/EMAG2/EMAG2_image_V2.tif",
 "expected_size":"133Mb"
},
    
# Australian Mag Data from AuScope portal (built by hand)
{
 "local_file":"../../Data/Resources/AusMagAll.tiff.zip",
 "md5":'',
 "url":"../../Data/Reference/AusMagAll.tiff.zip",
 "expected_size":"10Mb"
},
    
# Blue Marble Image Geotiff June 2004 (11Mb)  
{
 "local_file": "../../Data/Resources/BlueMarbleNG-TB_2004-06-01_rgb_3600x1800.TIFF",
 "md5": '55e1399674d43a26353d84c114d7ff80',
 "url":"http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=526294&cs=rgb&format=TIFF&width=3600&height=1800",
 "expected_size":"11Mb"
},
    
    
# Blue Marble Image with Bathymetry png - December 2004 
{
 "local_file": "../../Data/Resources/BlueMarbleNG-TB_2004-12-01_rgb_3600x1800.TIFF",
 "md5": '825da4b417ae19e24d2ea6db6cf8ad21',
 "url":"http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=526311&cs=rgb&format=TIFF&width=3600&height=1800",
 "expected_size":"11Mb"
},
    
# Etopo1 - color image 
{
 "local_file": "../../Data/Resources/color_etopo1_ice_low.tif.zip",
 "md5": '3f53d72e85393deb28874db0c76fbfcb',
 "url":"https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/color_etopo1_ice_low.tif.zip",
 "expected_size":"35Mb"
},
    
# Etopo1 - bw shaded relief image 
{
 "local_file": "../../Data/Resources/etopo1_grayscale_hillshade.tif.zip",
 "md5": 'f38b8dc6cd971d72e77edaa837d5b85c',
 "url":"https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/image/etopo1_grayscale_hillshade.tif.zip",
 "expected_size":"115Mb"
},
    
    
# ETOPO_c_geotiff_10800 ... cell centred data for ETOPO1 
{
 "local_file": "../../Data/Resources/ETOPO_c_geotiff_10800.tif.zip",
 "md5":"fbd0478d0974ca433cdd5cb3530d0d48",
 "url":"https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip",
 "expected_size":"320Mb"
},    
    
# Natural Earth Hypsometry Images (same grid size as ETOPO1 image above)
{
 "local_file":"../../Data/Resources/HYP_50M_SR_W.zip",
 "md5":"f3323bc6f38ddec98cff4936fb324ec5",
 "url":"http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/HYP_50M_SR_W.zip",
 "expected_size":"102Mb"
},   
 
# Natural Earth Ocean Bottom Images (same grid size as ETOPO1 image above)
{
 "local_file":"../../Data/Resources/OB_50M.zip",
 "md5":"69a2ed6ff8cf30b2dedd2d80461f9e43",
 "url":"http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/raster/OB_50M.zip",
 "expected_size":"50Mb"
},   
            
# Global age grid  [LOCAL]
{
 "local_file":"../../Data/Resources/global_age_data.3.6.z.npz",
 "md5":'',
 "url":"../../Data/Reference/global_age_data.3.6.z.npz",
 "expected_size":"22Mb"
},
 
# Global age grid (raw) [LOCAL]
{
 "local_file":"../../Data/Resources/global_age_data.3.6.xyz.zip",
 "md5":'',
 "url":"../../Data/Reference/global_age_data.3.6.xyz",
 "expected_size":"22Mb"
},


# Global second invariant of strain rate [LOCAL]
{
 "local_file":"../../Data/Resources/sec_invariant_strain_0.2.dat.zip",
 "md5":'',
 "url":"../../Data/Reference/sec_invariant_strain_0.2.dat.zip",
 "expected_size":"4Mb"
},
  
    
# Earthquake data for various regions as json files [LOCAL]
{
 "local_file":"../../Data/Resources/Earthquake_data.zip",
 "md5":'',
 "url":"../../Data/Reference/Earthquake_data.zip",
 "expected_size":"4Mb"
},
    
    
# Global surface velocity vector fields in various reference frames:
# http://gsrm.unavco.org/model/
# There are others ... 
    
{
 "local_file":"../../Data/Resources/velocity_AU.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_AU.nc",
 "expected_size":"1Mb"
},   
  
{
 "local_file":"../../Data/Resources/velocity_EU.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_EU.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_IN.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_IN.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_NA.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_NA.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_OK.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_OK.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_NNR.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_NNR.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_PA.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_PA.nc",
 "expected_size":"1Mb"
},   


{
 "local_file":"../../Data/Resources/velocity_TA.nc",
 "md5":'',
 "url":"../../Data/Reference/velocity_TA.nc",
 "expected_size":"1Mb"
},   

       
]

```

```{code-cell} ipython3
for resource in resource_list:
    print ("\nChecking {:s} ({})".format(resource["local_file"], resource["md5"]))
    report_cached_file(resource["local_file"],resource["md5"])
```

```{code-cell} ipython3
# md5sum("../../Data/Resources/EMAG2_image_V2.tif")
  

filename="Readme.md"

import hashlib 
from functools import partial

def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 4096), b''):
            d.update(buf)
    return d.hexdigest()

print(md5sum(filename))
```

```{code-cell} ipython3
for resource in resource_list:
    print ("\nDownloading {:s}".format(resource["local_file"]))
    download_cached_file(resource["url"], resource["local_file"], resource["md5"], resource["expected_size"] )
```

```{code-cell} ipython3
# Extract any files that were downloaded as zip archives (keep the originals to avoid re-downloading)

import zipfile
import glob

for zipped in glob.glob("../../Data/Resources/*.zip"):
    with zipfile.ZipFile(zipped) as zf:
        zf.extractall("../../Data/Resources")
        print ("Unzipped {}".format(zipped))
```

```{code-cell} ipython3
# if you add a new url / file to the download list above, add the md5 checksum as well !
# unless it is small enough to have a local copy. I would suggest still providing a url
# in the comments just in case.

md5sum("../../Data/Resources/etopo1_grayscale_hillshade.tif.zip") 
```

```{code-cell} ipython3

```
