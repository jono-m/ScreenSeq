import sys
from typing import *
from multiprocessing import Pool
import cv2
import numpy as np
import skimage
import scipy
import matplotlib.pyplot as plt
from scripts.PrintUpdate import *
from scripts.MatplotlibVideoPlayer import MatplotlibVideoPlayer

CropSettings = Tuple[float, float, float, float]


class Droplet:
    def __init__(self):
        self.image = None
        self.rp = None
        self.time = 0

def DetectDropletsInVideo(paths: str,
                          crop: CropSettings = None, nBackgroundFrames=100,
                          minDropletSize: int = None, maxDropletSize: int = None,
                          maxEccentricity: float = 0.7,
                          stride: int = 1,
                          debugMode=False,
                          parallel=False) -> List[Droplet]:

    if debugMode:
        fig, axs = plt.subplots(1, 3)
        fig.canvas.mpl_connect("close_event", sys.exit)
        videoPlayer = MatplotlibVideoPlayer(axs[0])

    dropletVideo = DropletVideo(paths[0], crop, stride)
    print("Identifying background...")
    backgroundFrames = [skimage.color.rgb2gray(frame) for frame in dropletVideo]
    dropletVideo.release()
    background = scipy.stats.mode(np.stack(backgroundFrames), axis=0)[0]
    droplets: List[Droplet] = []

    print("Detecting droplets...", end="")
    nFrames = 0

    if parallel and not debugMode:
        chunkSize = 1000
        pool = Pool()
        for chunkDroplets in pool.imap(DetectDropletsParallel, ParallelIterator(dropletVideo, background,
                                                                                minDropletSize, maxDropletSize,
                                                                                maxEccentricity), chunkSize):
            droplets += chunkDroplets
            nFrames += 1
            if nFrames % chunkSize == 0:
                PrintUpdate(str(nFrames) + " (" + str(len(droplets)) + " droplets found)")
        pool.close()
        return droplets

    for path in paths:
        dropletVideo = DropletVideo(path, crop, stride)
        for image in dropletVideo:
            nFrames += 1
            newDroplets = DetectDropletsInImage(image, background, minDropletSize, maxDropletSize, maxEccentricity, debugMode)
            for droplet in newDroplets:
                droplet.time = nFrames
            droplets += newDroplets
            PrintUpdate(str(nFrames) + " (" + str(len(droplets)) + " droplets found)")
            if debugMode:
                for droplet in newDroplets:
                    rp = droplet.rp
                    image[rp.coords[:, 0], rp.coords[:, 1]] = [255, 0, 0]
                videoPlayer.ShowFrame(image)
                axs[1].clear()
                axs[1].hist([droplet.rp.area for droplet in droplets])
                axs[1].set_xlabel("Area")
                axs[2].clear()
                axs[2].hist([droplet.rp.eccentricity for droplet in droplets])
                axs[2].set_xlabel("Eccentricity")
        dropletVideo.release()
    print("Done")
    PrintUpdate.Done()
    dropletVideo.release()
    return droplets


_disk = skimage.morphology.disk(3)


class DropletVideo:
    def __init__(self, path: str, crop: CropSettings, stride: int):
        self.video = cv2.VideoCapture(path)
        self.crop = crop
        self.stride = stride

    def __iter__(self):
        return self

    def __next__(self):
        for i in range(self.stride):
            self.video.read()
        ret, frame = self.video.read()
        if not ret:
            raise StopIteration
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        if self.crop is not None:
            frame = CropImage(frame, self.crop)
        return frame

    def release(self):
        self.video.release()


class ParallelIterator:
    def __init__(self, video: DropletVideo, background, minDropletSize, maxDropletSize, maxEccentricity):
        self.video = video
        self.params = [background, minDropletSize, maxDropletSize, maxEccentricity, False]
        self.i = 0

    def __iter__(self):
        return self

    def __next__(self):
        return [next(self.video), self.params]

def DetectDropletsParallel(parallel):
    return DetectDropletsInImage(parallel[0], *parallel[1])

def DetectDropletsInImage(image: np.ndarray, background: np.ndarray, minDropletSize: int,
                          maxDropletSize: int, maxEccentricity: float, debugMode: bool):
    grayscale = skimage.color.rgb2gray(image)
    foreground = np.abs(grayscale - background)
    segmented = Segment(foreground)
    rps = skimage.measure.regionprops(segmented)
    droplets = []
    for rp in rps:
        if ((minDropletSize is not None and rp.area < minDropletSize) or
                (maxDropletSize is not None and rp.area > maxDropletSize)):
            continue
        if maxEccentricity is not None and rp.eccentricity > maxEccentricity:
            continue
        droplet = Droplet()
        droplet.image = image[rp.bbox[0]:rp.bbox[2], rp.bbox[1]:rp.bbox[3]]
        droplet.rp = rp
        droplets.append(droplet)
    return droplets


def Segment(foreground: np.ndarray) -> np.ndarray:
    foreground = skimage.filters.gaussian(foreground, 3) ** 0.5
    segmented = foreground > skimage.filters.threshold_triangle(foreground)
    segmented = skimage.morphology.binary_closing(segmented, _disk)
    segmented = scipy.ndimage.binary_fill_holes(segmented)
    segmented = skimage.segmentation.clear_border(segmented)
    return skimage.measure.label(segmented)


def CropImage(image: np.ndarray, cropFraction: CropSettings):
    height, width = image.shape[:2]

    crop = [cropFraction[0] * width, cropFraction[1] * height,
            (1 - cropFraction[2]) * width, (1 - cropFraction[3]) * height]
    crop = [round(c) for c in crop]
    image = image[crop[1]:crop[3], crop[0]:crop[2]]
    return image
