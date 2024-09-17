import matplotlib.pyplot as plt

class MPLVideoPlayer:
    def __init__(self):
        plt.ion()
        self.img = None

    def Play(self, imgs, **kwargs):
        for img in imgs:
            self.Frame(img, **kwargs)

    def Frame(self, img, **kwargs):
        if self.img is None:
            self.img = plt.imshow(img, **kwargs)
            plt.show()
        else:
            img.set_data(img)
        plt.draw()
        plt.pause(0.01)


class PrintUpdate:
    lastText = ""

    def __init__(self, text):
        print("\b" * len(PrintUpdate.lastText), end="")
        PrintUpdate.lastText = text
        print(text, end="")

    @staticmethod
    def Done():
        PrintUpdate.lastText = ""
        print("")
