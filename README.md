# Spatial Augmented Multiomics Interface (SAMI) 
©July 11, 2023 University of Florida Research Foundation, Inc. All Rights Reserved.

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

We present a method for the simultaneous analysis of spatial metabolome, lipidome, and glycome from a single tissue section using mass spectrometry imaging. Our workflow includes a computational pipeline called __Spatial Augmented Multiomics Interface (SAMI)__ that offers multiomics integration, high dimensionality clustering, spatial anatomical mapping with matched multiomics features, and metabolic pathway enrichment to providing unprecedented insights into the spatial distribution and interaction of these biomolecules in mammalian tissue biology.
![Main Figure](https://github.com/XinBiostats/SAMI/blob/main/figures/main.png)

## Implement
SAMI can be run through two different ways:

### 1. Using Docker (Recommended):
We have pre-configured the environment for you using Docker, which ensures a consistent and reliable environment and make it easy to get started.

#### Steps:
- Clone SAMI from Github Repository:
```bash
git clone https://github.com/XinBiostats/SAMI
```
- Download [testing data](https://www.dropbox.com/scl/fo/qjdk94golwij84xfii15b/h?rlkey=etrdydm1iw86ntcprbem2wivn&dl=1) from Dropbox and put it in "./SAMI/datasets/".  
- Download [libaries](https://www.dropbox.com/scl/fo/d1ssah39qkfr9er640dai/APjGcNoiOX9zW84QvhCJtuw?rlkey=joosfjb7vjs72t3t4pmyxrrol&st=uqcbacfy&dl=1) for pathway enrichment analysis from Dropbox and put them in "./SAMI/lib/".  
- Download __Docker__ desktop from [Docker website](https://www.docker.com), and install it on your machine.  
- Open the Terminal or PowerShell(Windows), then run below command with your SAMI path:
```bash
docker run -it --rm --user root -e GRANT_SUDO=yes -p 8888:8888 -v "YOUR_SAMI_PATH:/home/jovyan/work" xinbiostats/sami:latest

example: docker run -it --rm --user root -e GRANT_SUDO=yes -p 8888:8888 -v "/Users/xin.ma/Desktop/SAMI:/home/jovyan/work" xinbiostats/sami:latest
```
- Find the highlighted link in your teminal and copy it to your browser. The link will not be the same, but will show up at same place.
![docker_link](https://github.com/XinBiostats/SAMI/blob/main/figures/docker_link.png)
- All set! You can play with our Demo now. ([demo](https://github.com/XinBiostats/SAMI/blob/main/demo)) 

### 2. Using Conda:
Create your own environment for SAMI.(Due to potential incompatibility issues caused by different operating systems and versions, it is recommended to use Docker.）

#### Steps:
- Clone SAMI from Github Repository:
```bash
git clone https://github.com/XinBiostats/SAMI
```
- Download [testing data](https://www.dropbox.com/scl/fo/qjdk94golwij84xfii15b/h?rlkey=etrdydm1iw86ntcprbem2wivn&dl=1) from Dropbox and put it in "./SAMI/datasets/".  
- Download [libaries](https://www.dropbox.com/scl/fo/d1ssah39qkfr9er640dai/APjGcNoiOX9zW84QvhCJtuw?rlkey=joosfjb7vjs72t3t4pmyxrrol&st=uqcbacfy&dl=1) for pathway enrichment analysis from Dropbox and put them in "./SAMI/lib/".  
- Open the Terminal or PowerShell(Windows), then install requirements:
```bash
conda env create -f environment.yml
```
-  Activate SAMI environment, find your R installation's home directory.
```bash
conda activate SAMI
   
R RHOME
```
For example:
```bash
/Users/xin.ma/anaconda3/envs/SAMI/lib/R
```
Open [./SAMI/pathway.py](https://github.com/XinBiostats/SAMI/blob/main/SAMI/pathway.py) scripts, update os.environ['R_HOME'] using your R home directory.
```python
os.environ['R_HOME'] = "/Users/xin.ma/anaconda3/envs/SAMI/lib/R"
```
- All set! You can play with our Demo now. ([demo](https://github.com/XinBiostats/SAMI/blob/main/demo)) 
