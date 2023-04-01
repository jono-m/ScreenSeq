import numpy as np
from PIL import Image
import skimage


def FindDroplets(image, area_min, area_max, circularity_min, show=False):
    seg = image > skimage.filters.threshold_otsu(image)
    seg = skimage.morphology.erosion(seg)
    seg = skimage.morphology.area_closing(seg)
    regions = skimage.measure.regionprops(skimage.measure.label(seg))
    droplets = []
    for r in regions:
        circularity = r.area * 4 * np.pi / (r.perimeter_crofton ** 2)
        if circularity >= circularity_min and area_min <= r.area <= area_max:
            r.circularity = circularity
            r.coords_np = tuple(np.asarray(r.coords).T)
            droplets.append(r)

    if show:
        image = np.round(Scale01(image)*255).astype(np.uint8)
        i = np.stack([image, image, image], axis=-1)
        for droplet in droplets:
            i[:, :, 1][droplet.coords_np] = 255

        Image.fromarray(i).show()

    return droplets


def Scale01(i):
    i = i.astype(float)
    return (i - np.min(i)) / (np.max(i) - np.min(i))


def Show(i):
    return Image.fromarray(np.floor(Scale01(i) * 255).astype(np.uint8)).show()
