from scripts.DetectDropletsInVideo import DetectDropletsInVideo, DropletVideo
from DropletFluorescentImageProcessor import ComputeDropletIntensities
import numpy as np
import nd2reader
from PIL import  Image


class Analysis:
    def __init__(self):
        self.droplets = DetectDropletsInVideo(r"../data/ColorMicroscopy/dropletVideo.mp4",
                                              minDropletSize=None,
                                              maxDropletSize=None,
                                              maxEccentricity=.5,
                                              crop=(0.35, 0.3, 0, 0.3),
                                              debugMode=False,
                                              stride=0,
                                              parallel=False)

        self.fps = 29.97
        self.images = DropletVideo(r"../data/ColorMicroscopy/dropletVideo.mp4",crop=(0.35, 0.3, 0, 0.3),stride=0)
        self.images = [Image.fromarray(image) for image in self.images]
        self.colors = [d.meanColor for d in self.droplets]
        self.sizes = [d.propDict['area'] for d in self.droplets]

        # Read in the two images
        nd2file = nd2reader.ND2Reader(r"../data/FluorMicroscopy/Inlet1.nd2")
        nd2file.iter_axes = 'c'
        self.noSwitchFluorescence = nd2file.get_frame(0)
        self.noSwitchBF = nd2file.get_frame(1)
        nd2file.close()

        nd2file = nd2reader.ND2Reader(r"../data/FluorMicroscopy/Inlet2.nd2")
        nd2file.iter_axes = 'c'
        self.switchFluorescence = nd2file.get_frame(0)
        self.switchBF = nd2file.get_frame(1)
        nd2file.close()

        self.intensitiesNoSwitch = ComputeDropletIntensities(self.noSwitchFluorescence, self.noSwitchBF)
        self.intensitiesSwitch = ComputeDropletIntensities(self.switchFluorescence, self.switchBF)

if __name__ == "__main__":
    Analysis()