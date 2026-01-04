from scripts.DetectDropletsInVideo import DetectDropletsInVideo
import pickle
if __name__ == "__main__":
    droplets = DetectDropletsInVideo(
        r"../data/ColorMicroscopy/screen_3_6_2025.mp4",
        minDropletSize=8000,
        maxDropletSize=None,
        maxEccentricity=0.7,
        crop=(0, .3, 0, .15),
        debugMode=False,
        stride=10,
        parallel=True)

    with open("droplets.pkl", "wb") as file:
        pickle.dump(droplets, file)