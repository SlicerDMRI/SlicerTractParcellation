# TractParcellation
Automatic parcellation of white matter tractography using [White matter analysis (WMA)](https://whitematteranalysis.readthedocs.io/en/latest/index.html) method and [O’Donnell Research Group (ORG) atlas](https://dmri.slicer.org/atlases/) .


<div align="center">
  <img src="https://github.com/JoshuaKening/SlicerTractParcellation/assets/129246293/a1a48305-d42e-4bdc-8787-bb5b52752d90" alt="Tract Parcellation" style="width:40%;">
  <br>
  <em>Fig 1: Visual representation of brain tract parcellation.</em>
</div>

# Tutorial
The Tutorial for the Powpoint version can be found in [releases](https://github.com/JoshuaKening/SlicerTractParcellation/releases).
* Start 3D Slicer
* Go to ```TractParcellation``` module
* Install ```WMA```
* Download ```WM atlas``` (optional)
* Select input:
  * Select one of the three input methods
    * ```From Slicer```: Select from ```Input Fiberbundle```
    * ```From File```: Choose file path from ```Input File```
    * ```From Directory```: Choose folder path from ```Input Folder```
* Select output path from ```Output Folder```
* Select the registration mode: ```affine``` or ```affine + nonlinear```
* Decide whether to clear unnecessary intermediate results
* Select the number of threads involved in the computation
* Click ```Apply```


# User interface
### Message Box
- **Functionality**: Displays the installation status of the WMA toolkit and WM atlas.

### Installation Module
- **Functionality**: Provides an interface for automatically downloading the WMA toolkits and data sets.

### IO Module
- **Functionality**: Provides selection for input and output paths and types.

### Advanced Parameters Module
- **Functionality**:
  - Offers two registration methods.
  - Ability to automatically clear intermediate results.
  - Option to select the number of processors involved in the calculation.
# License
The contents of this repository are released under a [Slicer license](https://github.com/Slicer/Slicer/blob/main/License.txt).
