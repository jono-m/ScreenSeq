import importlib
import DropletFluorescentImageProcessor


def Run():
    importlib.reload(DropletFluorescentImageProcessor)
    DropletFluorescentImageProcessor.Process()
