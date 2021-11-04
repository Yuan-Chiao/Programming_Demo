import numpy as np
from dataio import load_single_image
from scipy.ndimage import zoom
import SimpleITK as sitk
from dataio import import_data_filename, write_nii

from utils import get_nii_data

def image_rcv(path, crop_image_size, y_pred):

    x0 = get_nii_data(path)
    y0 = get_nii_data(path[:-10]+'mask.nii.gz')

    # crop based on the labels
    tmp = y0 - np.min(y0)
    I0 = np.sum(tmp, axis=(1, 2)) # Dimension at Z
    index0 = I0>0
    result_idx0=np.array(range(len(index0)))
    crop0=result_idx0[index0]

    I1 = np.sum(tmp, axis=(0, 2)) # Dimension at Y
    index1 = I1>0
    result_idx1=np.array(range(len(index1)))
    crop1=result_idx1[index1]

    I2 = np.sum(tmp, axis=(0, 1)) # Dimension at X
    index2 = I2>0
    result_idx2=np.array(range(len(index2)))
    crop2=result_idx2[index2]

    new_shape = (crop0[-1]-crop0[0]+1, crop1[-1]-crop1[0]+1, crop2[-1]-crop2[0]+1)

    y = np.around(zoom(y_pred, np.divide(new_shape, crop_image_size), order=1))

    y0[crop0[0]:crop0[-1]+1, crop1[0]:crop1[-1]+1, crop2[0]:crop2[-1]+1]=y
    return y0


def image_rcv_raw(path, crop_image_size, y_pred):

    x0 = get_nii_data(path)
    y0 = get_nii_data(path[:-10]+'mask.nii.gz')
    new_shape = y0.shape
    y = np.around(zoom(y_pred, np.divide(new_shape, crop_image_size), order=1))
    return y


def load_image_correct_oritation(filename):
    image_obj = sitk.ReadImage(filename)
    direction = image_obj.GetDirection()
    origin = np.asarray(image_obj.GetOrigin())
    spacing = np.asarray(image_obj.GetSpacing())
    affine = SimpleRot(direction)
    data = sitk.GetArrayFromImage(image_obj).astype(np.float32)
    image_size = np.asarray(data.shape)
    center = (image_size/2-1)*spacing
    affine.SetCenter([center[0], center[1], center[2]])
    image_obj.SetOrigin([0,0,0])
    image_obj.SetDirection([1,0,0,0,1,0,0,0,1])
    data = sitk.GetArrayFromImage(image_obj).astype(np.float32)

    # newimage = resample(image_obj, affine) # resample is not working for placenta

    # print('testfromhere2')
    # data = sitk.GetArrayFromImage(newimage).astype(np.float32)
    # print(np.max(data))
    # print(np.min(data))
    # print(data.shape)

    # data = sitk.GetArrayFromImage(newimage).astype(np.float32)
    return data


def resample(image, transform):
    # Output image Origin, Spacing, Size, Direction are taken from the reference
    # image in this call to Resample
    reference_image = image
    interpolator = sitk.sitkCosineWindowedSinc
    default_value = 0
    # print(transform)
    return sitk.Resample(image, reference_image, transform,
                         interpolator, default_value)


def SimpleRot(matrix):
    dimension = 3
    affine = sitk.AffineTransform(3)
    matrix = np.array(matrix).reshape((dimension, dimension))
    target = np.array([[1,0,0],
                       [0,1,0],
                       [0,0,-1]])
    transform_matrix = target @ np.linalg.inv(matrix)
    affine.SetMatrix(transform_matrix.ravel())
    return affine


def normlize_mean_std(tmp):
    tmp_std = np.std(tmp) + 0.0001
    tmp_mean = np.mean(tmp)
    tmp = (tmp - tmp_mean) / tmp_std
    return tmp


def normlize_min_max(tmp):
    tmp_max = np.amax(tmp)
    tmp_min = np.amin(tmp)
    tmp = (tmp - tmp_min) / (tmp_max-tmp_min)
    return tmp
