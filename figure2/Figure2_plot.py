import matplotlib.pyplot as plt
import matplotlib.transforms
import matplotlib.patches
import pandas
import typing
import numpy as np
import seaborn

from Figure2_analysis import Analysis
from PIL import Image

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 4
plt.rcParams['font.family'] = 'arial'
plt.rcParams['axes.linewidth'] = 0.7
plt.rcParams['xtick.major.pad'] = 1
plt.rcParams['ytick.major.pad'] = 1


class Plot:
    def __init__(self, analysis: Analysis):
        self.a = analysis

        plt.close(plt.gcf())

        self.fig = plt.figure(figsize=(4, 3), layout='none', dpi=200)
        self.axes = {}

        hr = [3, 2, 3]
        gs = plt.GridSpec(nrows=3, ncols=2, hspace=0.1, wspace=0.1, height_ratios=hr)
        self.axes['a'] = self.fig.add_subplot(gs[0, 0])
        self.axes['b'] = self.fig.add_subplot(gs[0, 1])
        self.axes['c'] = self.fig.add_subplot(gs[1, :])
        bottomgs = plt.GridSpec(nrows=3, ncols=3, hspace=0.4, wspace=0.2, width_ratios=[2, 2, 3], height_ratios=hr)
        self.axes['d'] = self.fig.add_subplot(bottomgs[2, 0])
        self.axes['e'] = self.fig.add_subplot(bottomgs[2, 1])
        bottomgs2 = plt.GridSpec(nrows=3, ncols=3, hspace=0.4, wspace=0.4, width_ratios=[2, 2, 3], height_ratios=hr)
        self.axes['f'] = self.fig.add_subplot(bottomgs2[2, 2])

        self._DropletPlot(self.axes['c'])

        self._DropletPopulationPlot(self.axes['f'])

        self._Label()
        self._Image(r"../data/Photos/chipPhoto.png", self.axes["a"])
        self._Image(r"../data/ColorMicroscopy/junctionImage.png", self.axes["b"])

        self._Image(Image.fromarray(self.a.noSwitchFluorescence), self.axes["d"])
        self._Image(Image.fromarray(self.a.switchFluorescence), self.axes["e"])
        self.axes["d"].set_xlabel("Without switching")
        self.axes["e"].set_xlabel("With switching")
        plt.show()

    def _Label(self):
        for label, ax in self.axes.items():
            ax.text(0.0, 1.0, label, transform=(
                    ax.transAxes + matplotlib.transforms.ScaledTranslation(-0.05, 0.08, self.fig.dpi_scale_trans)),
                    fontsize=8, va="top", ha="right", fontfamily="Arial", fontweight="bold")

    def _Image(self, image: typing.Union[Image.Image, str], ax: plt.Axes):
        image = Image.open(image) if isinstance(image, str) else image
        ax.imshow(np.asarray(image), aspect="equal")
        ax.set_xticks([])
        ax.set_yticks([])
        w, h = ax.get_window_extent().size
        asp = w / h
        iasp = image.width / image.height
        if asp > iasp:
            gap = (image.height - image.width / asp) / 2
            ax.set_ylim([image.height - gap, gap])
        else:
            gap = (image.width - image.height * asp) / 2
            ax.set_xlim([gap, image.width - gap])

    def _DropletPopulationPlot(self, ax: plt.Axes):
        data = pandas.DataFrame()
        data['Fluorescence intensity'] = self.a.intensitiesNoSwitch + self.a.intensitiesSwitch
        data['Switching'] = ['No switching'] * len(self.a.intensitiesNoSwitch) + ['Switching'] * len(self.a.intensitiesSwitch)
        data.replace([np.inf, -np.inf], np.nan, inplace=True)
        seaborn.histplot(data=data, x='Fluorescence intensity', hue="Switching", ax=ax, common_norm=False, kde=True)
        ax.set_ylabel("Droplet count", labelpad=1)
        ax.set_xlabel("Fluorescence intensity (au)")

        ax.legend(ax.lines, [t.get_text() for t in ax.get_legend().texts], fontsize=4, frameon=False)

    def _DropletPlot(self, ax: plt.Axes):
        insetTimes = [0, 5, 10, 15, 20]
        insetFrames = [min(int(x * self.a.fps), len(self.a.sizes) - 1) for x in insetTimes]

        # Compute limits
        heightRatios = [2, 1, 1]
        sizes = np.asarray(self.a.sizes) * 100 / 75

        # Limits
        tick = 15
        bottom = int(sizes.min())
        bottom = bottom - ((bottom - 1) % tick + 1)
        top = int(sizes.max())
        top = top + (tick - (top % tick))
        ticks = np.arange(bottom, top + 1, tick)

        sizeHeight = (top - bottom) * 1.05
        colorHeight = sizeHeight * (heightRatios[1] / heightRatios[2])
        insetHeight = sizeHeight * (heightRatios[0] / heightRatios[2])
        totalHeight = sizeHeight + colorHeight + insetHeight

        # Make plot
        ax.plot(sizes, c="black", linewidth=0.5)
        ax.set_xlim([0, sizes.size - 1])
        ax.set_ylim([bottom, bottom + totalHeight])
        ax.set_ylabel("Droplet length ($\mu m$)")
        ax.yaxis.set_label_coords(-0.06, (sizeHeight / 2) / totalHeight)
        ax.set_yticks(ticks)
        ax.set_xlabel("Time (s)", labelpad=0)
        ax.set_xticks(insetFrames, insetTimes)

        ax.spines['left'].set_bounds(bottom, bottom + sizeHeight + colorHeight)
        ax.spines['right'].set_bounds(bottom, bottom + sizeHeight + colorHeight)
        ax.spines['top'].set_visible(False)

        lw = ax.spines['left'].get_linewidth()

        # add color plot
        colors = np.asarray(self.a.colors)[None, :, :]
        ax.imshow(colors / 255,
                  extent=(0,
                          sizes.size - 1,
                          bottom + sizeHeight,
                          bottom + sizeHeight + colorHeight),
                  aspect="auto")
        ax.axhline(bottom + sizeHeight, linewidth=lw, color="black")
        ax.axhline(bottom + sizeHeight + colorHeight, linewidth=lw, color="black")

        images = self.a.images

        imagePadding = 10

        imageHeight = insetHeight - imagePadding
        maxWidth = imageHeight * images[0].width / images[0].height

        # with those images at the height and packed evenly, what width is each?
        imageSpacing = 60
        availableWidth = sizes.size - imageSpacing * (len(insetFrames) - 1)
        imageWidth = (availableWidth / len(insetFrames))
        imageWidth = min(imageWidth, maxWidth)
        imageAspect = imageWidth / imageHeight
        images = [i.crop((i.width - i.height * imageAspect, 0, i.width, i.height)) for i in images]
        w, h = ax.get_window_extent().size
        pixelAspect = (h / totalHeight) * (sizes.size / w)
        imageWidth = pixelAspect * imageWidth

        imageSpacing = (sizes.size - imageWidth * len(insetFrames)) / (len(insetFrames) - 1)
        for i, t in enumerate(insetFrames):
            l = ax.axvline(t, (sizeHeight) / totalHeight, 1, linewidth=lw, color="black", linestyle="solid")
            l.set_clip_on(False)

            top = totalHeight + bottom
            bot = top - imageHeight
            left = i * (imageWidth + imageSpacing)
            right = left + imageWidth
            extents = [left, right, bot, top]
            ax.imshow(np.asarray(images[t]),
                      aspect='auto',
                      extent=extents, zorder=2)
            patch = matplotlib.patches.Rectangle((left, bot), imageWidth-1, imageHeight,
                                                 edgecolor=ax.spines['left'].get_edgecolor(), facecolor="none",
                                                 zorder=5,
                                                 linewidth=lw)
            patch.set_clip_on(False)
            ax.add_patch(patch)
