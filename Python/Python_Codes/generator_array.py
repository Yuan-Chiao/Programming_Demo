import numpy as np
import tensorflow.keras
from one_hot_label import multi_class_labels
from dataio import write_label_nii, write_nii
from patch3d import patch
from image_process import normlize_mean_std


class Generator(tensorflow.keras.utils.Sequence):
    'Generates data for Keras, based on array data X and Y'

    def __init__(self, X, Y, batch_size=32, patch_size=[32, 32, 32], labels=[1]):
        'Initialization'
        self.batch_size = batch_size
        self.patch_size = np.asarray(patch_size)
        self.X = X
        self.Y = Y
        self.labels = labels
        self.image_size = np.asarray(X.shape[1:])
        self.n_subject = X.shape[0]
        self.loc_patch = patch(np.asarray(X.shape[1:]), patch_size, stride=[2, 2, 2])
        tmp = np.indices((self.n_subject, self.loc_patch.n_patch))
        self.indices = np.column_stack((np.ndarray.flatten(tmp[0]), np.ndarray.flatten(tmp[1])))
        self.on_epoch_end()

        # print('n_subject: '+str(self.n_subject))
        # self.loc_patch.__info__()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(self.n_subject*self.loc_patch.n_patch/self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'

        # Generate indexes of the batch, batch means how many patch is used once
        X = np.zeros((self.batch_size, self.patch_size[0], self.patch_size[1], self.patch_size[2], 1), dtype=np.float32)
        Y = np.zeros((self.batch_size, self.patch_size[0], self.patch_size[1], self.patch_size[2], len(self.labels)), dtype=np.float32)  # channel_last by default
        # Generate data
        # print('generator index: '+str(index))
        # print('generator len: '+str(self.__len__()))
        sub_index = self.indices[np.arange(self.batch_size) + index*self.batch_size]
        for batch_count in range(self.batch_size):
            image_index = sub_index[batch_count, 0]
            patch_index = sub_index[batch_count, 1]

            label_data = self.Y[image_index]
            label = self.loc_patch.__get_single_patch__(label_data, patch_index)
            Y[batch_count] = multi_class_labels(label, self.labels)

            image_data = self.X[image_index]
            image = self.loc_patch.__get_single_patch__(image_data, patch_index)
            image = normlize_mean_std(image)
            X[batch_count, :, :, :, 0] = image


        # TODO augmentation
        if np.random.uniform() > 0.5:
            X = X[:,::-1,:]
            Y = Y[:,::-1,:]
        if np.random.uniform() > 0.5:
            X = X[::-1,:,:]
            Y = Y[::-1,:,:]
        if np.random.uniform() > 0.5:
            X = X[:,:,::-1]
            Y = Y[:,:,::-1]

        # X = X0[..., np.newaxis]
        # change to one-hot-label
        # for label_index in range(len(self.labels)):
        #     Y[Y0 == self.labels[label_index], label_index] = 1

        # print('labels ' + str(Y.shape))
        # print('image ' + str(image.shape))
        # print(np.mean(image_data))
        # write_nii(X[0,:,:,:,0], 'image'+ str(index)+'.nii')
        # write_label_nii(Y[0,...],'generator')
        # write_nii(np.argmax(Y[1,...],3),'generatorss.nii')
        return X, Y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        np.random.shuffle(self.indices)
        # np.random.shuffle(self.patch_indices)
        # self.indexes = np.arange(len(self.list_IDs))
