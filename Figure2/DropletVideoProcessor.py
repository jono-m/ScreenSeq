import cv2
from PIL import Image, ImageEnhance
import numpy as np
import skimage
import scipy
from Util import PrintUpdate


class DropletVideoProcessor:
    def __init__(self):
        self.images = []
        self.droplets = []
        self.segmented = []
        self.fps = -1

    def Read(self, path: str, fps, startTime=None, endTime=None):
        print("Reading...", end="")
        video = cv2.VideoCapture(path)
        i = 0
        self.images = []
        self.fps = fps
        endFrame = None if endTime is None else int(endTime*fps)
        startFrame = None if startTime is None else int(startTime * fps)
        while video.isOpened():
            ret, frame = video.read()
            if not ret or (endFrame is not None and i >= endFrame):
                break

            PrintUpdate(str(i))
            i += 1

            if startFrame is not None and i < startFrame:
                continue

            image = Image.fromarray(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
            self.images.append(image)
        print("Done")
        PrintUpdate.Done()
        video.release()

    def Preprocess(self, rotation, crop, brightness):
        print("Preprocessing...")
        for i, image in enumerate(self.images):
            PrintUpdate(str(i))
            image = image.rotate(rotation, Image.Resampling.BILINEAR).crop(crop)
            brightener = ImageEnhance.Brightness(image)
            image = brightener.enhance(brightness)
            self.images[i] = image
        PrintUpdate.Done()
        print("Done")

    def Process(self):
        print("Converting...")
        grayscale = np.stack([skimage.color.rgb2gray(i) for i in self.images])
        foregrounds = self._RemoveBackground(grayscale)
        print("Segmenting...", end="")
        for i, image in enumerate(self.images):
            rps, segmented = self._Segment(foregrounds[i])
            self.segmented.append(segmented)
            self.droplets.append(rps)
            PrintUpdate(str(i))
        PrintUpdate.Done()
        print("Done")

    def _RemoveBackground(self, grayscale):
        print("Background removal...")
        background = scipy.stats.mode(grayscale, axis=0).mode
        foregrounds = np.abs(grayscale - background)
        return foregrounds

    def _Segment(self, grayscale):
        f = skimage.filters.gaussian(grayscale, 3) ** 0.5
        binarized = f > skimage.filters.threshold_otsu(f)
        clean = skimage.morphology.binary_closing(binarized, skimage.morphology.disk(10))
        clean = scipy.ndimage.binary_fill_holes(clean)
        clean = skimage.segmentation.clear_border(clean)
        rps = skimage.measure.regionprops(skimage.measure.label(clean))
        return (rps, clean)
