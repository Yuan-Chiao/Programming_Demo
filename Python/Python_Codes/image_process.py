import numpy as np
from dataio import load_single_image
from scipy.ndimage import zoom
import SimpleITK as sitk
from dataio import import_data_filename, write_nii

from utils import get_nii_data

def crop_pad3D(x, target_size, shift=[0, 0, 0]):
    'crop or zero-pad the 3D volume to the target size'
    small = 0
    y = np.ones(target_size, dtype=np.float32)*small
    current_size = x.shape
    pad_size = [0, 0, 0]

    for dim in range(3):
        if current_size[dim] > target_size[dim]:
            pad_size[dim] = 0
        else:
            pad_size[dim] = int(np.ceil((target_size[dim] - current_size[dim])/2.0))
    # pad first
    x1 = np.pad(x, [[pad_size[0], pad_size[0]], [pad_size[1], pad_size[1]], [pad_size[2], pad_size[2]]], 'constant', constant_values=small)
    # crop on x1
    start_pos = np.ceil((np.asarray(x1.shape) - np.asarray(target_size))/2.0)
    start_pos = start_pos.astype(int)
    y = x1[(shift[0]+start_pos[0]):(shift[0]+start_pos[0]+target_size[0]),
           (shift[1]+start_pos[1]):(shift[1]+start_pos[1]+target_size[1]),
           (shift[2]+start_pos[2]):(shift[2]+start_pos[2]+target_size[2])]
    return y


def crop_edge3D(x, target_size):
    # print(x.shape)
    small = 0
    y = np.ones(target_size, dtype=np.float32)*small

    tmp = x - np.min(x)
    # print(np.sum(tmp, axis=(1,2)).shape)
    I0 = np.sum(tmp, axis=(1, 2))
    index0 = I0>0
    I1 = np.sum(tmp, axis=(0, 2))
    index1 = I1>0
    I2 = np.sum(tmp, axis=(0, 1))
    index2 = I2>0
    new_shape = (sum(index0), sum(index1), sum(index2))
    image_cropped = x[index0, ...]
    image_cropped = image_cropped[:, index1, :]
    image_cropped = image_cropped[..., index2]
    # print(new_shape)
    # print('crop image: '+str(image_cropped.shape))

    y = zoom(image_cropped, (target_size[0]/len(index0), target_size[1]/len(index1), target_size[2]/len(index2)))
    return y


def crop_edge_pair(path, target_size,n):
    # small = 0 # No Need Actually (Yuan)
    # y = np.ones(target_size, dtype=np.float32)*small  # No Need Actually (Yuan)

    # x0 = load_image_correct_oritation(path)    #  image
    # # print(x0.shape) # Brain (126, 159, 117) # Pla (52, 256, 256)
    # # y0 = load_image_correct_oritation((path[:-12]+'brain_tissue_labels.nii.gz'))  # label for Brain
    # y0 = load_image_correct_oritation(path[:-10]+'mask.nii.gz')  # label for Placenta
    # # print(y0.shape) # Brain (126, 159, 117) # Pla (52, 256, 256)

    x0 = get_nii_data(path)
    y0 = get_nii_data(path[:-10]+'mask.nii.gz')

    # x=np.flip(x0,0)
    # print(x.shape)
    # write_nii(x, 'X'+str(n)) # +str(path[-39:-28])
    # y=np.flip(y0,0)
    # print(y0.shape)
    # write_nii(y, 'Y'+str(n))
    # # print(np.min(x0)) # 0
    # # print(np.max(x0)) # 9

    # crop based on the labels
    tmp = y0 - np.min(y0)
    # print(np.min(y0)) # 0
    # print(np.max(y0)) # 9
    # print(np.sum(tmp, axis=(1,2)).shape) # Brain (126,) # Pla (52,)
    I0 = np.sum(tmp, axis=(1, 2)) # Dimension at Z
    # print(I0)
    # print(I0.shape) # (126,)
    index0 = I0>0
    result_idx0=np.array(range(len(index0)))
    crop0=result_idx0[index0]

    I1 = np.sum(tmp, axis=(0, 2)) # Dimension at Y
    # print(I1)
    # print(I1.shape) # (159,)
    index1 = I1>0
    result_idx1=np.array(range(len(index1)))
    crop1=result_idx1[index1]

    I2 = np.sum(tmp, axis=(0, 1)) # Dimension at X
    # print(I2)
    # print(I2.shape) # (117,)
    index2 = I2>0
    result_idx2=np.array(range(len(index2)))
    crop2=result_idx2[index2]

    new_shape = (crop0[-1]-crop0[0]+1, crop1[-1]-crop1[0]+1, crop2[-1]-crop2[0]+1)
    print('Placenta segmentation dimension for subject ID=', n, ':', new_shape)

    y_cropped = y0[crop0[0]:crop0[-1]+1, crop1[0]:crop1[-1]+1, crop2[0]:crop2[-1]+1]
    x_cropped = x0[crop0[0]:crop0[-1]+1, crop1[0]:crop1[-1]+1, crop2[0]:crop2[-1]+1]

    # print('crop image: '+str(image_cropped.shape))
    # # print(np.divide(target_size, new_shape))
    x = zoom(x_cropped, np.divide(target_size, new_shape), order=1)
    # print(new_shape)
    # print(np.divide(target_size, new_shape)) # [0.98765432 1.08910891 1.05882353]
    # print(x.shape) # Brain (80, 110, 90) Pla (100, 100, 100)
    # print(np.max(x)) # 1663.5151 for Placenta
    # print(np.min(x)) # 5.939394 for Placenta
    y = np.around(zoom(y_cropped, np.divide(target_size, new_shape), order=0))
    # print(y.shape) # Brain (80, 110, 90) Pla (100, 100, 100)
    # print(np.max(y)) # 1 for Placenta
    # print(np.min(y)) # 0 for Placenta

    return x, y


def crop_edge_pair_raw(path, target_size,n):

    x0 = get_nii_data(path)
    y0 = get_nii_data(path[:-10]+'mask.nii.gz')
    new_shape = y0.shape
    print(new_shape)
    x = zoom(x0, np.divide(target_size, new_shape), order=1)
    y = np.around(zoom(y0, np.divide(target_size, new_shape), order=0))
    return x, y



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


def resample(image, transform):  # What is this for? (Yuan)
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
