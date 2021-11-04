# U-Net 3D model core
from tensorflow.keras.layers import Conv3D, MaxPool3D, concatenate, Input, Dropout, PReLU, Conv3DTranspose, BatchNormalization
from tensorflow.keras.models import Model
# from keras_contrib.layers import CRF


def unet_core(x, filter_size=8, kernel_size=(3, 3, 3)):
    x = Conv3D(filters=filter_size,
               kernel_size=kernel_size,
               padding='same',
               kernel_initializer='he_normal')(x)
    x = BatchNormalization()(x)
    x = PReLU()(x)
    x = Conv3D(filters=filter_size,
               kernel_size=kernel_size,
               padding='same',
               kernel_initializer='he_normal')(x)
    x = BatchNormalization()(x)
    x = PReLU()(x)
    return x


def unet3d(patch_size, n_label):
    input_layer = Input(shape=patch_size)
    d1 = unet_core(input_layer, filter_size=96, kernel_size=(3, 3, 3))
    l = MaxPool3D(strides=(2, 2, 2))(d1)
    d2 = unet_core(l, filter_size=192, kernel_size=(3, 3, 3))
    l = MaxPool3D(strides=(2, 2, 2))(d2)
    d3 = unet_core(l, filter_size=384, kernel_size=(3, 3, 3))
    l = MaxPool3D(strides=(2, 2, 2))(d3)

    b = unet_core(l, filter_size=768, kernel_size=(3, 3, 3))
    Dropout(0.5)

    l = Conv3DTranspose(filters=384, kernel_size=(2, 2, 2),  padding='same', strides=2, kernel_initializer='he_normal')(b)
    l = concatenate([l, d3], axis=-1)
    u3 = unet_core(l, filter_size=384, kernel_size=(3, 3, 3))
    l = Conv3DTranspose(filters=192, kernel_size=(2, 2, 2),  padding='same', strides=2, kernel_initializer='he_normal')(u3)
    l = concatenate([l, d2], axis=-1)
    u2 = unet_core(l, filter_size=192, kernel_size=(3, 3, 3))
    l = Conv3DTranspose(filters=96, kernel_size=(2, 2, 2),  padding='same', strides=2, kernel_initializer='he_normal')(u2)
    l = concatenate([l, d1], axis=-1)
    u1 = unet_core(l, filter_size=96, kernel_size=(3, 3, 3))

    output_layer = Conv3D(filters=n_label, kernel_size=(1, 1, 1), activation='sigmoid')(u1)
    # output_layer = CRF(n_label)
    model = Model(input_layer, output_layer)
    return model
