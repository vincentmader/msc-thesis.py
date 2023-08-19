from typing import Optional


class GridspecSubplot():
    __slots__ = ["title", "xlabel", "ylabel"]

    def __init__(
        self, 
        title:      Optional[str]   = None,
        xlabel:     Optional[str]   = None,
        ylabel:     Optional[str]   = None,
    ):
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel

    def draw(self, _):
        raise Exception("Not implemented.")
