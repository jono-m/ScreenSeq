from DropletVideoProcessor import DropletVideoProcessor
from DropletFluorescentImageProcessor import ComputeDropletIntensities
import numpy as np
import nd2reader


class Analysis:
    def __init__(self):
        self.dvp = DropletVideoProcessor()
        self.dvp.Read(r"../data/ColorMicroscopy/dropletVideo.mp4",
                      29.73,
                      endTime=25)
        self.dvp.PreprocessAll(-0.5, (210, 190, 512, 350), 1.2)
        self.dvp.ProcessAll()

        self.sizes = []
        self.colors = []

        for image, droplets in zip(self.dvp.images, self.dvp.droplets):
            size = np.mean([droplet.axis_major_length for droplet in droplets])
            coords = np.concatenate([np.asarray(droplet.coords) for droplet in droplets])
            color = np.asarray(image)[coords[:, 0], coords[:, 1]].mean(axis=0)

            self.colors.append(color)
            self.sizes.append(size)

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
