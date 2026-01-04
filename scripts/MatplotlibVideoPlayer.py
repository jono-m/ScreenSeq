import matplotlib.pyplot as plt
from typing import *

import numpy as np


class MatplotlibVideoPlayer:
    def __init__(self, axis: plt.Axes):
        self.axis = axis
        self.lastSize = None
        self.img = None

    def Play(self, images: List[np.ndarray], **kwargs):
        for img in images:
            self.ShowFrame(img, **kwargs)

    def ShowFrame(self, img: np.ndarray, delay=0.01, **kwargs):
        if self.lastSize != img.shape:
            self.axis.clear()
            self.img = self.axis.imshow(img, **kwargs)
        else:
            self.img.set_data(img)
        plt.pause(delay)
