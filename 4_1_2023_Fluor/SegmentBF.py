import numpy as np
import scipy.ndimage
import seaborn
from scipy.ndimage import binary_fill_holes
import skimage
from PIL import Image
from pathlib import Path
import matplotlib.pyplot as plt

rootPath = Path(r"C:\Users\jonoj\Documents\Fluorescence2\BF")


def Segment(i):
    # Segmentation algorithm is slightly different depending on if we are dealing with a droplet bilayer.
    segmentOverlap = [False,
                      True,
                      False,
                      True,
                      False,
                      False,
                      True,
                      True]
    toSegment = list(rootPath.iterdir())
    SegmentMonolayer(toSegment[i])
    # if segmentOverlap[i]:
    #     SegmentBilayer(toSegment[i])
    # else:
    #     SegmentMonolayer(toSegment[i])


# def SegmentBilayer(path):
#     image = np.asarray(Image.open(path).convert(mode="L"))
#     thresholded = skimage.feature.canny(image, 2, low_threshold=1, high_threshold=10)
#     thresholded2 = skimage.morphology.area_closing(thresholded, 50)
#     thresholded = np.where(thresholded, False, thresholded2)
#
#     segmented = np.zeros_like(thresholded)
#     regions = skimage.measure.regionprops(skimage.measure.label(thresholded))
#     goodregions = []
#     for region in regions:
#         if region.axis_minor_length > 0 and region.axis_major_length / region.axis_minor_length < 1.5:
#             coords_np = tuple(np.asarray(region.coords).T)
#             segmented[coords_np] = True
#             goodregions.append(region)
#     Preview(image, segmented)


def SegmentMonolayer(path):
    image = np.asarray(Image.open(path).convert(mode="L"))
    thresholded = skimage.feature.canny(image, 2, low_threshold=0, high_threshold=20)
    radius = 10
    hough = skimage.transform.hough_circle(thresholded, radius)
    accums, cxs, cys, radii = skimage.transform.hough_circle_peaks(hough, [radius], threshold=0.5*np.max(hough), min_xdistance=int(radius*1.5),
                                                                   min_ydistance=int(radius*1.5))

    segmented = np.zeros_like(thresholded)
    for cx, cy in zip(cxs, cys):
        radius = 6
        if cx-radius < 0 or cy-radius < 0 or cx+radius >= segmented.shape[0] or cy+radius >= segmented.shape[0]:
            continue
        segmented[skimage.draw.circle_perimeter(cy, cx, radius)] = True

    # regions = skimage.measure.regionprops(skimage.measure.label(thresholded))
    # goodregions = regions
    # for region in regions:
    #     if region.axis_minor_length > 0 and region.axis_major_length / region.axis_minor_length < 1.5:
    #         goodregions.append(region)

    # quartiles = np.percentile([r.area for r in goodregions], [25, 75])
    # iqr = quartiles[1]-quartiles[0]
    # threshLower = quartiles[0] - 1.5*iqr
    # threshUpper = quartiles[1] + 1.5 * iqr
    # goodregions = [region for region in goodregions if threshLower <= region.area <= threshUpper]
    # for region in goodregions:
    #     coords_np = tuple(np.asarray(region.coords).T)
    #     segmented[coords_np] = True
    Preview(image, segmented)


fig, axsAll = plt.subplots(2, 3, sharex='all', sharey='all')
axs = (x for x in axsAll.flatten())


def Preview(imageNP, segmentedNP):
    overlay = np.stack([(imageNP.astype(float) * 0.5).astype(np.uint8)] * 3, axis=-1)
    overlay[segmentedNP, 1] = 255
    plt.sca(next(axs))
    plt.imshow(imageNP, cmap="gray", interpolation="nearest")
    plt.sca(next(axs))
    plt.imshow(segmentedNP, cmap="gray", interpolation="nearest")
    plt.sca(next(axs))
    plt.imshow(overlay)


Segment(0)
Segment(2)

for ax in axsAll.flatten():
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis('off')
    ax.set_aspect('equal')
fig.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
plt.show()
