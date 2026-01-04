class PrintUpdate:
    lastText = ""

    def __init__(self, text: str):
        print("\b" * len(PrintUpdate.lastText), end="")
        PrintUpdate.lastText = text
        print(text, end="")

    @staticmethod
    def Done():
        PrintUpdate.lastText = ""
        print("")
