import matplotlib.pyplot as plt
import numpy as np


def view3plane(I, figurename=None, pos=None):
    if not figurename:
        figurename = 'temp'
    if not pos:
        image_size = I.shape
        pos = np.asarray(image_size)//2
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].imshow(I[pos[0], :, :], cmap='gray'), ax[0].title.set_text('dim=0')
    ax[1].imshow(I[:, pos[1], :], cmap='gray'), ax[1].title.set_text('dim=1')
    ax[2].imshow(I[:, :, pos[2]], cmap='gray'), ax[2].title.set_text('dim=2')
    plt.tight_layout()
    plt.savefig(figurename+'.jpg')
    plt.close(fig)
