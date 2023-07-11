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

## Installation
1. Download SAMI:
```bash
git clone https://github.com/XinBiostats/SAMI
```
2. Requirements: SAMI is implemented by both Python and R. To install requirements
```bash
pip install -r requirements.txt
```
## Usage
1. Download data from [Dropbox](https://www.dropbox.com/scl/fo/9ntjdocvj3rlopjrjnrxu/h?dl=0&rlkey=3ipjglxydmiioxika44bw2lzv) and put downloaded data into datasets folder.
2. We created a demo ([demo.ipynb](https://github.com/XinBiostats/SAMI/blob/main/demo/demo.ipynb)) to demonstrate how to use SAMI. The results will be displayed inline and stored in [results](https://github.com/XinBiostats/SAMI/tree/main/results) folder.

