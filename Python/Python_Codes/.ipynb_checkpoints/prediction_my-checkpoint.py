import numpy as np
from patch3d import patch
from image_process import normlize_mean_std, crop_pad3D
import metrics
from dataio import write_nii
from one_hot_label import restore_labels

from image_recovered import image_rcv
from image_recovered import image_rcv_raw
from utils import save_image

import metrics2
import os

class predict(object):
    'run the model on test data'

    def __init__(self, model, image_size, patch_size=(32, 32, 32), labels=[1]):
        'Initialization'
        self.patch_size = patch_size
        # self.data_file = data_file
        self.labels = labels
        self.stride = [8, 8, 8]
        self.image_size = np.asarray(image_size)
        self.model = model
        self.loc_patch = patch(self.image_size, patch_size, self.stride)

    def __run__(self, X):
        'test on one image each time'
        output_size = np.append((self.loc_patch.size_after_pad),len(self.labels))
        Y0 = np.zeros(output_size,  dtype=np.float32)
        batch_size = 2
        X0 = np.zeros([batch_size]+self.patch_size+[1])

        for index in range(np.ceil(self.loc_patch.n_patch/batch_size).astype(int)):
            batch_of_patch_index = np.asarray(range(batch_size)) + batch_size*index
            batch_of_patch_index[batch_of_patch_index>=self.loc_patch.n_patch]=0
            # get one batch_size
            for n, selected_patch in enumerate(batch_of_patch_index):
                patch = self.loc_patch.__get_single_patch__(X, selected_patch)
                patch = normlize_mean_std(patch)
                X0[n] = patch[np.newaxis, ..., np.newaxis]
            # predict this patches
            prediction = self.model.predict(X0)
            prediction_flip = self.model.predict(X0[::-1, ...])
            prediction = prediction + prediction_flip[::-1, ...]
            # put the label back
            for n, selected_patch in enumerate(batch_of_patch_index):
                Y0 = self.loc_patch.__put_single_patch__(Y0, np.squeeze(prediction[n]), selected_patch)

        Y = restore_labels(Y0, self.labels)

        result = crop_pad3D(Y, self.image_size)
        return result


def test(x, model, image_size, patch_size, labels):
    prediction = np.zeros(x.shape)
    predictor = predict(model, image_size, patch_size, labels)
    for n in range(x[0]):
        prediction[n] = predictor.__run__(x[n])
    return prediction


def evaluate(x, y_true, model, image_size, patch_size, labels,MainFolder,subject_list,Initial_Idx,crop_image_size,affines):
    n_subject = x.shape[0]
    predictor = predict(model, image_size, patch_size, labels)
    dice = 0.0
    dice_label_0 = 0.0
    dice_label_1 = 0.0
    preci = 0.0
    sensi = 0.0
    specif = 0.0
    np.set_printoptions(precision=3)
    for n in range(n_subject):
        y_pred = predictor.__run__(x[n])

        tmp_dice = metrics.dice_multi_array(y_true[n], y_pred, labels)
        print(str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77])+': Dice='+str(tmp_dice))
        dice=np.append(dice,np.mean(tmp_dice))
        dice_label_0 = np.append(dice_label_0,tmp_dice[0]) # for label "0": non-placenta
        dice_label_1= np.append(dice_label_1,tmp_dice[1])  # for label "1": placenta

        tmp_preci = metrics2.precision(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77])+': Precision='+str(tmp_preci))
        preci=np.append(preci,tmp_preci)

        tmp_sensi = metrics2.sensitivity(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77])+': Sensitivity='+str(tmp_sensi))
        sensi=np.append(sensi,tmp_sensi)

        tmp_specif = metrics2.specificity(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77])+': Specificity='+str(tmp_specif))
        specif=np.append(specif,tmp_specif)

        Folder_Pred = MainFolder+'Data_Placenta_Pred/'
        if not os.path.exists(Folder_Pred):
            os.makedirs(Folder_Pred)
        save_image(y_pred, affines[n], MainFolder+'Data_Placenta_Crop/'+'y_pred_'+ str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77]) + '.nii')
        y_pred_rcv=image_rcv(subject_list[Initial_Idx+n], crop_image_size, y_pred)
        save_image(y_pred_rcv, affines[n], Folder_Pred + 'y_pred_'+ str(subject_list[Initial_Idx+n][66:70])+'_'+str(subject_list[Initial_Idx+n][76:77]) + '.nii')

    dice=dice[1:]
    dice_label_0=dice_label_0[1:]
    dice_label_1=dice_label_1[1:]
    preci=preci[1:]
    sensi=sensi[1:]
    specif=specif[1:]

    print('*****   '+'Dice'+'   *****')
    print(dice)
    print('Dice mean:' + str(np.mean(dice)))
    print('Dice std:' + str(np.std(dice)))
    print(dice_label_0)
    print('Dice label 0 mean:' + str(np.mean(dice_label_0))) # Dice average for label "0": non-placenta
    print('Dice label 0 std:' + str(np.std(dice_label_0)))
    print(dice_label_1)
    print('Dice label 1 mean:' + str(np.mean(dice_label_1))) # Dice average for label "1": placenta
    print('Dice label 1 std:' + str(np.std(dice_label_1)))

    print('*****   '+'Precision'+'   *****')
    print(preci)
    print('Precision mean:' + str(np.mean(preci)))
    print('Precision std:' + str(np.std(preci)))

    print('*****   '+'Sensitivity'+'   *****')
    print(sensi)
    print('Sensitivity mean:' + str(np.mean(sensi)))
    print('Sensitivity std:' + str(np.std(sensi)))

    print('*****   '+'Specificity'+'   *****')
    print(specif)
    print('Specificity mean:' + str(np.mean(specif)))
    print('Specificity std:' + str(np.std(specif)))

def evaluate_raw(x, y_true, model, image_size, patch_size, labels,MainFolder,subject_list,Initial_Idx,crop_image_size,affines):
    n_subject = x.shape[0]
    predictor = predict(model, image_size, patch_size, labels)
    metric = 0.0
    metric2 = 0.0
    metric3 = 0.0
    metric4 = 0.0
    np.set_printoptions(precision=3)
    for n in range(n_subject):
        y_pred = predictor.__run__(x[n])
        tmp = metrics.dice_multi_array(y_true[n], y_pred, labels)
        print(str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76])+': '+str(tmp))
        metric += tmp

        tmp2 = metrics2.precision(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76])+': '+str(tmp2))
        metric2 += tmp2

        tmp3 = metrics2.sensitivity(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76])+': '+str(tmp3))
        metric3 += tmp3

        tmp4 = metrics2.specificity(y_pred, y_true[n])
        print(str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76])+': '+str(tmp4))
        metric4 += tmp4

        save_image(y_pred, affines[n], MainFolder+'Data_Placenta_Crop_raw/'+'y_pred_'+ str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76]) + '.nii')
        y_pred_rcv=image_rcv_raw(subject_list[Initial_Idx+n], crop_image_size, y_pred)
        save_image(y_pred_rcv, affines[n], MainFolder+'Data_Placenta_Pred_raw/' + 'y_pred_'+ str(subject_list[Initial_Idx+n][65:69])+'_'+str(subject_list[Initial_Idx+n][75:76]) + '.nii')

    return metric/n_subject, metric2/n_subject, metric3/n_subject, metric4/n_subject
