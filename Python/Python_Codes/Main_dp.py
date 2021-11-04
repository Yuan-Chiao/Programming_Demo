import os
from model import unet3d
from tensorflow.keras import optimizers
from tensorflow.keras.utils import multi_gpu_model
from tensorflow.keras.callbacks import ModelCheckpoint
from generator_array import Generator
import numpy as np
from prediction_my import evaluate
from image_process import crop_edge_pair
from dataio import import_data_filename, write_nii
from one_hot_label import redefine_label
import metrics
from viewer import view3plane
# Set up environ parameters
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0,1"
# check the gpu connection
from tensorflow.python.client import device_lib
# print(device_lib.list_local_devices())

import tensorflow as tf
from scipy.ndimage import zoom
from image_process import load_image_correct_oritation
from utils import get_nii_data
from image_recovered import image_rcv
from image_process import crop_edge_pair_raw
from prediction_my import evaluate_raw
from utils import get_nii_affine
from utils import save_image

def main():
    subject_index = range(369)   # 369 is the total number of subjects
    n_subject = len(subject_index)
    shift_pixel = [0, 0, 0]
    zoom_scale = [1, 1, 1]
    crop_image_size = [100,100,100]  # [superior-inferior, anterior-posterior, left-right]
    image_size = [int(crop_image_size[0]*zoom_scale[0]), int(crop_image_size[1]*zoom_scale[1]), int(crop_image_size[2]*zoom_scale[2])]
    patch_size = [64, 64, 64]
    labels = [0, 1]  # 0: no placenta; 1: placenta

    MainFolder='/home/yuan/Desktop/DeepLearning/Placenta_seg/' # Put all codes in this folder
    DataFolder = MainFolder+'Data_Placenta/' # Put all placenta image data in this folder

    print('********************   '+'Read Data'+'   ********************')
    subject_list = import_data_filename(DataFolder, '_placenta_ana.nii.gz') # '_placenta_ana.nii.gz' is the T2-weighted anatomical MR images of placenta without segmentation

    X = np.zeros((n_subject, image_size[0], image_size[1], image_size[2]), dtype=np.float32) # For T2 images
    Y = np.zeros((n_subject, image_size[0], image_size[1], image_size[2]), dtype=np.float32) # For segmentation images
    affines = np.zeros((n_subject, 4, 4), dtype=np.float32) # For affines

    k = 0
    for n in range(0, n_subject):
        X[k], Y[k] = crop_edge_pair(subject_list[n], crop_image_size,n)  # X: For T2 images, Y: For segmentation images
        affines[k] = get_nii_affine(subject_list[n])
        k += 1
    print('..The number of subjects for training and image size:' + str(X.shape))

    ID = 0
    with open('test_id.txt', 'w') as patient_txt:
        for n in range(ID, n_subject):
            Folder_Crop = MainFolder+'Data_Placenta_Crop/' # This folder will save the cropped T2 and segmentation images
            if not os.path.exists(Folder_Crop):
                os.makedirs(Folder_Crop)
            save_image(X[n], affines[n], Folder_Crop + 'x_'+ str(subject_list[n][66:70])+'_'+str(subject_list[n][76:77]) + '.nii')
            save_image(Y[n], affines[n], Folder_Crop + 'y_'+ str(subject_list[n][66:70])+'_'+str(subject_list[n][76:77]) + '.nii')
            patient_txt.write(subject_list[n][66:77]+'\n')
    
    # Separate the training set into training and validaiton groups for DL model
    # I used 240 (80%) subjects for training and 60 (20%) subjects for validation
    training_generator = Generator(X[0:240, ...], Y[0:240, ...], batch_size=2, patch_size=patch_size, labels=labels)
    print('..training_n_subjecr: ' + str(training_generator.n_subject))
    validation_generator = Generator(X[240:300, ...], Y[240:300, ...], batch_size=2, patch_size=patch_size, labels=labels)
    print('..validation_n_subjecr: ' + str(validation_generator.n_subject))

    print('********************   '+'Build Model'+'   ********************')
    single_model = unet3d(patch_size+[1], len(labels))
    model = multi_gpu_model(single_model, gpus=2)
    optimizer = optimizers.Adam(lr=1e-4)
    model.compile(optimizer=optimizer,
              loss='binary_crossentropy', # I used 'binary_crossentropy' for loss function. Other option: 'categorical_crossentropy'
              metrics=[metrics.dice_multi, 'binary_crossentropy'])

    print('********************   '+'Model Fitting'+'   ********************')
    check_pointer = ModelCheckpoint(filepath='mymodel_weights.h5',verbose=1, save_best_only=True, save_weights_only=True)
    tensor_board_callback = tf.keras.callbacks.TensorBoard(log_dir='./log')
    model.load_weights('mymodel_weights.h5')   # Enable this if you want to use pre-trained model
    hist = model.fit_generator(generator=training_generator,
                               steps_per_epoch=128,
                               epochs=50,
                               validation_data=validation_generator,
                               validation_steps=8,
                               verbose=2,
                               callbacks=[check_pointer]
                               )
    model.save_weights('mymodel_weights.h5')

    print('********************   '+'Prediction'+'   ********************')
    #  This will export the predicted segmentation based on the trained DL model and export Dice, Precision, Sensitivity, Specificity
    model.load_weights('mymodel_weights.h5')

    print('**********   '+'Training Data'+'   **********')
    Initial_Idx=0
    Sbj_Len=240
    evaluate(X[Initial_Idx:Sbj_Len], Y[Initial_Idx:Sbj_Len], single_model, image_size, patch_size, labels,MainFolder,subject_list,Initial_Idx,crop_image_size,affines[Initial_Idx:Sbj_Len])

    print('**********   '+'Validation Data'+'   **********')
    Initial_Idx=240
    Sbj_Len=300
    evaluate(X[Initial_Idx:Sbj_Len], Y[Initial_Idx:Sbj_Len], single_model, image_size, patch_size, labels,MainFolder,subject_list,Initial_Idx,crop_image_size,affines[Initial_Idx:Sbj_Len])

    print('**********   '+'Testing Data'+'   **********')
    Initial_Idx=300
    Sbj_Len=369
    evaluate(X[Initial_Idx:Sbj_Len], Y[Initial_Idx:Sbj_Len], single_model, image_size, patch_size, labels,MainFolder,subject_list,Initial_Idx,crop_image_size,affines[Initial_Idx:Sbj_Len])

if __name__ == '__main__':
    main()

print('end of code')