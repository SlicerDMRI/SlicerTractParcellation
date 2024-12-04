cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerWMA)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/SlicerDMRI/SlicerTractParcellation")
set(EXTENSION_CATEGORY "Diffusion")
set(EXTENSION_CONTRIBUTORS "Fan Zhang (UESTC; BWH; HMS), Kening Zhang (UESTC), Lauren O'Donnell (BWH; HMS)")
set(EXTENSION_DESCRIPTION "This is a Slicer extension to provide access to the whitematteranalysis package.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/JoshuaKening/SlicerTractParcellation%20/main/image.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/JoshuaKening/SlicerTractParcellation%20/main/Screenshot.png")
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
