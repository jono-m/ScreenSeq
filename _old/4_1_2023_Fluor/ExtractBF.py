import numpy as np
from pathlib import Path
from nd2reader import ND2Reader
from PIL import Image

for path in Path(r"C:\Users\jonoj\Documents\Fluorescence2").iterdir():
    if path.suffix != ".nd2":
        continue
    outputPath = path.parent / "BF"
    outputPath.mkdir(exist_ok=True)
    with ND2Reader(str(path.absolute())) as nd2file:
        bfChannel = nd2file.metadata["channels"].index("Brightfield")
        nd2file.iter_axis = "c"
        bf = nd2file.get_frame(bfChannel)
        bf = (bf - np.min(bf)) / (np.max(bf) - np.min(bf))
        bf = (bf * 255).astype(np.uint8)
        Image.fromarray(bf).save(outputPath / (path.stem + ".png"))
