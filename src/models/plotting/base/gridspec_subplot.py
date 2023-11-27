from typing import Optional


class GridspecSubplot():
    __slots__ = ["title", "xlabel", "ylabel"]

    def __init__(
        self, 
        title:      Optional[str]   = None,
        xlabel:     Optional[str]   = None,
        ylabel:     Optional[str]   = None,
    ):
        self.title  = title  if title  is not None else ""
        self.xlabel = xlabel if xlabel is not None else ""
        self.ylabel = ylabel if ylabel is not None else ""

    def draw(self, _):
        raise Exception("Not implemented.")
