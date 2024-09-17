import numpy as np
import scipy.ndimage
import seaborn
from scipy.ndimage import binary_fill_holes
import skimage
from PIL import Image
from pathlib import Path
import matplotlib.pyplot as plt

rootPath = Path(r"BF")


def DoSegment(bfImage, radius):
    blurred = skimage.filters.gaussian(bfImage, 1)
    grad = skimage.filters.scharr(blurred)
    edges = grad > skimage.filters.threshold_mean(grad)
    edges = skimage.morphology.remove_small_objects(edges, 100)
    distance = 10
    radii = np.arange(8, 15)
    hough = skimage.transform.hough_circle(edges, radii)
    _, circleXs, circleYs, circleRadii = skimage.transform.hough_circle_peaks(hough, radii,
                                                                              min_xdistance=distance,
                                                                              min_ydistance=distance,
                                                                              normalize=False)

    segmented = np.zeros_like(edges)
    for circleX, circleY, circleRadius in zip(circleXs, circleYs, circleRadii):
        if circleX - circleRadius < 0 or circleY - circleRadius < 0 or \
                circleX + circleRadius >= segmented.shape[0] or circleY + circleRadius >= \
                segmented.shape[0]:
            continue
        segmented[skimage.draw.circle_perimeter(circleY, circleX, circleRadius)] = True

    return segmented


def DoSegmentRinse(preview=False):
    path = rootPath / "Rinse2001.png"
    bfImage = np.asarray(Image.open(path).convert(mode="L"))

    seg = skimage.morphology.remove_small_objects(bfImage > 40, 100)

    regions = skimage.measure.regionprops(skimage.measure.label(seg))
    segmented = np.zeros_like(seg)
    for region in regions:
        if region.eccentricity > 0.5:
            continue
        segmented[skimage.draw.disk(region.centroid, 2)] = True

    Image.fromarray(segmented).convert(mode="1").save(
        rootPath / "Segmented" / (path.stem + ".png"))
    if preview:
        Preview(bfImage, segmented)


def DoSegmentCrosstalk(preview=False):
    path = rootPath / "Image20h001.png"
    bfImage = np.asarray(Image.open(path).convert(mode="L"))

    seg = skimage.morphology.remove_small_objects(bfImage > 45, 100)

    regions = skimage.measure.regionprops(skimage.measure.label(seg))
    segmented = np.zeros_like(seg)
    for region in regions:
        if region.eccentricity > 0.5 or region.area > 5000:
            continue
        segmented[skimage.draw.disk(region.centroid, 2)] = True

    Image.fromarray(segmented).convert(mode="1").save(
        rootPath / "Segmented" / (path.stem + ".png"))
    if preview:
        Preview(bfImage, segmented)


def DoSegmentInlets(preview=False):
    path = rootPath / "Inlet8.png"
    bfImage = np.asarray(Image.open(path).convert(mode="L"))

    seg = skimage.morphology.remove_small_objects(bfImage > 45, 100)

    regions = skimage.measure.regionprops(skimage.measure.label(seg))
    segmented = np.zeros_like(seg)
    for region in regions:
        if region.eccentricity > 0.5 or region.area > 1000:
            continue
        segmented[skimage.draw.disk(region.centroid, 2)] = True

    Image.fromarray(segmented).convert(mode="1").save(
        rootPath / "Segmented" / (path.stem + ".png"))
    if preview:
        Preview(bfImage, segmented)


def Preview(imageNP, segmentedNP):
    fig, axsAll = plt.subplots(1, 3, sharex='all', sharey='all')
    axs = (x for x in axsAll.flatten())
    overlay = np.stack([(imageNP.astype(float) * 0.75).astype(np.uint8)] * 3, axis=-1)
    overlay[segmentedNP, 1] = 255
    plt.sca(next(axs))
    plt.imshow(imageNP, cmap="gray", interpolation="nearest")
    plt.sca(next(axs))
    plt.imshow(segmentedNP, cmap="gray", interpolation="nearest")
    plt.sca(next(axs))
    plt.imshow(overlay)

    for ax in axsAll.flatten():
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axis('off')
        ax.set_aspect('equal')
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()


DoSegmentInlets(preview=True)
plt.show()
