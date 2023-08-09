import SimpleITK as sitk
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description='Simple script to binarize signed-distance nrrd file')
parser.add_argument("-i", "--input", help = "File path of input nrrd file", required=True)
parser.add_argument("-o", "--output", help = "Filename of output file without extension", required=True)
args = vars(parser.parse_args())
input_file = args["input"]
input_dir = os.path.dirname(input_file)
output_file = os.path.join(input_dir,args["output"]+'.mhd')

image = sitk.ReadImage(input_file)
image_dim = image.GetSize()

if len(image_dim) == 3: # it's a single-channel 3D nrrd file
    image_direction = image.GetDirection()
    image_bin_spacing = image.GetSpacing()
    image_bin_origin = image.GetOrigin()
    image_bin_direction = [image_bin_spacing[0], 0., 0., 0., image_bin_spacing[1], 0., 0., 0., image_bin_spacing[2]]

    image_array = sitk.GetArrayFromImage(image).swapaxes(0,2) # need to perform swapaxes as SimpleITK uses x-y-z and Numpy uses z-y-x
    mask_array = np.zeros([image_dim[0], image_dim[1], image_dim[2]], dtype="uint8")
    mask_array[image_array < 0] = 1
    image_bin = sitk.GetImageFromArray(mask_array.swapaxes(0,2)) # need to perform swapaxes as SimpleITK uses x-y-z and Numpy uses z-y-x
    image_bin.SetSpacing(image_bin_spacing)
    image_bin.SetOrigin(image_bin_origin)
    image_bin.SetDirection(image_bin_direction)
    sitk.WriteImage(image_bin, output_file)

else: # it's a multi-channel 3D nrrd file 
    image_direction = image.GetDirection()
    image_bin_spacing = image.GetSpacing()[:3]
    image_bin_origin = image.GetOrigin()[:3]
    image_bin_direction = [image_bin_spacing[0], 0., 0., 0., image_bin_spacing[1], 0., 0., 0., image_bin_spacing[2]]
    num_comp = image_dim[3]
    analyte = {0:"CALC",
               1:"LRNC",
               2:"IPH",
               3:"PVAT",
               4:"MATX"}
    for eachcomp in range(num_comp):
        image_array = sitk.GetArrayFromImage(image[:, :, :, eachcomp]).swapaxes(0,2) # need to perform swapaxes as SimpleITK uses x-y-z and Numpy uses z-y-x
        mask_array = np.zeros([image_dim[0], image_dim[1], image_dim[2]], dtype="uint8")
        mask_array[image_array < 0] = 1
        image_bin = sitk.GetImageFromArray(mask_array.swapaxes(0,2)) # need to perform swapaxes as SimpleITK uses x-y-z and Numpy uses z-y-x)
        image_bin.SetSpacing(image_bin_spacing)
        image_bin.SetOrigin(image_bin_origin)
        image_bin.SetDirection(image_bin_direction)
        sitk.WriteImage(image_bin, output_file.replace('.mhd',f'_{analyte[eachcomp]}.mhd'))
