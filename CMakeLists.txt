cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerWMA)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://dmri.slicer.org/atlases")
set(EXTENSION_CATEGORY "Diffusion")
set(EXTENSION_CONTRIBUTORS "Fan Zhang (UESTC; BWH; HMS), Kening Zhang (UESTC), Steve Pieper (Isomics), Lauren O'Donnell (BWH; HMS)")
set(EXTENSION_DESCRIPTION "The Anatomical Tract Parcellation module segments fiber clusters and white matter tracts from input whole brain tractography.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/SlicerDMRI/SlicerTractParcellation/main/SlicerTractParcellation.png")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(AnatomicalTractParcellation)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
