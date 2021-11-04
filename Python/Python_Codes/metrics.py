# metrics
import tensorflow as tf
from tensorflow.keras import backend as K
import numpy as np


def dice_arrary(y_true, y_pred):
    epsilon = K.epsilon()
    intersection = np.sum(np.logical_and(y_true, y_pred))
    return (2. * intersection + epsilon) / (np.sum(y_true) +
                                            np.sum(y_pred) +
                                            epsilon)


def dice_tensor(y_true, y_pred):
    ' calc the dice on tf.tensor object'
    tmp = tf.logical_and(y_true, y_pred)
    intersection = tf.reduce_sum(tf.cast(tmp, dtype=tf.float32))
    dice = (2. * intersection + K.epsilon()) / (tf.reduce_sum(tf.cast(y_true, dtype=tf.float32))
                                                + tf.reduce_sum(tf.cast(y_pred, dtype=tf.float32))
                                                + K.epsilon())
    return dice


# def dice_coefficient_loss(y_true, y_pred):
#     return 1-dice_tensor(y_true, y_pred)


def dice_multi(y_true, y_pred):
    'work on the tensor'
    dice_value = 0.0
    n_labels = y_pred.get_shape().as_list()[-1]
    prediction = tf.argmax(y_pred, -1)
    for i in range(n_labels):
        yi_true = tf.slice(y_true, [0, 0, 0, 0, i], [-1, -1, -1, -1, 1])
        yi_true = tf.cast(yi_true[..., 0], dtype=tf.bool)
        yi_pred = tf.equal(prediction, i)
        dice_value += dice_tensor(yi_true, yi_pred)
    return dice_value/n_labels


def dice_multi_array(y_true, y_pred, labels):
    n_labels = len(labels)
    dice_value = np.zeros(n_labels, dtype=np.float32)
    for i in range(n_labels):
        yi_true = (y_true == labels[i])
        yi_pred = (y_pred == labels[i])
        # print(yi_true.shape)
        # print(yi_true.dtype)
        dice_value[i] = dice_arrary(yi_true, yi_pred)
    return dice_value


def dice_multi_loss(y_true, y_pred):
    return 1-dice_multi(y_true, y_pred)


def jaccard_distance_loss(y_true, y_pred, smooth=1e-5):
    """
    Jaccard = (|X & Y|)/ (|X|+ |Y| - |X & Y|)
            = sum(|A*B|)/(sum(|A|)+sum(|B|)-sum(|A*B|))

    The jaccard distance loss is usefull for unbalanced datasets. This has been
    shifted so it converges on 0 and is smoothed to avoid exploding or disapearing
    gradient.

    Ref: https://en.wikipedia.org/wiki/Jaccard_index

    @url: https://gist.github.com/wassname/f1452b748efcbeb4cb9b1d059dce6f96
    @author: wassname
    """
    intersection = K.sum(K.abs(y_true * y_pred), axis=-1)
    sum_ = K.sum(K.abs(y_true) + K.abs(y_pred), axis=-1)
    jac = (intersection + smooth) / (sum_ - intersection + smooth)
    return (1 - jac) * smooth
