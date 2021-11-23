from pydantic import BaseModel
from typing import Dict, List, Union
from pathlib import Path
import math
from matplotlib.figure import Figure
from shapely.geometry import Polygon
import numpy as np


class Point3D(BaseModel):
    l: float
    x: float
    y: float
    z: float


class Crosssection(BaseModel):
    id: str = ""
    points: List[Point3D] = []

    @classmethod
    def from_dam_data(
        cls, surfaceline: str, locations: Dict
    ) -> Union[None, "Crosssection"]:
        result = Crosssection()
        args = surfaceline.split(";")
        result.id = args[0]

        if not result.id in locations.keys():
            print(f"Could not find location for crosssection '{result.id}'")
            return None

        xref = locations[result.id]["x"]
        yref = locations[result.id]["y"]

        xs = [float(args[i]) for i in range(1, len(args), 3)]
        ys = [float(args[i]) for i in range(2, len(args), 3)]
        zs = [float(args[i]) for i in range(3, len(args), 3)]
        ls = []
        # if the location is correct we expect values decreasing in distance (-)
        # and later increasing in distance (+)
        past_reference_point = False
        for i, p in enumerate(zip(xs, ys)):
            dl = math.hypot(p[0] - xref, p[1] - yref)

            if i == 0:
                ls.append(-dl)
            else:
                past_reference_point = abs(ls[-1]) > dl
                if past_reference_point:
                    dl *= -1
                ls.append(dl)

        # check if we have a reference point (if there is no negative value in ls we have a problem)
        if min(ls) >= 0.0:
            print(f"Could not find reference point from crosssection {result.id}")
            return None

        result.points = [
            Point3D(x=xs[i], y=ys[i], z=zs[i], l=round(ls[i], 3))
            for i in range(len(xs))
        ]
        return result

    @property
    def zmin(self) -> float:
        if len(self.points) <= 0:
            raise ValueError("Trying to find zmin on a crosssection without points")
        return min([p.z for p in self.points])

    @property
    def zmax(self) -> float:
        if len(self.points) <= 0:
            raise ValueError("Trying to find zmax on a crosssection without points")
        return max([p.z for p in self.points])

    def to_surfaceline(self) -> str:
        line = f"{self.id};"
        for p in self.points:
            line += f"{p.x:.03f};{p.y:.03f};{p.z:.03f};"
        return line[:-1]  # remove last semi colon


class CrosssectionList(BaseModel):
    crosssections: List[Crosssection] = []

    @classmethod
    def from_dam_csv(
        cls,
        filename_surfacelines: Union[Path, str],
        filename_locations: Union[Path, str],
    ):
        result = CrosssectionList()
        surfacelines = open(filename_surfacelines, "r").readlines()
        locationlines = open(filename_locations, "r").readlines()

        # first define the locations
        locations = {}
        for line in locationlines[1:]:
            args = line.split(";")
            locations[args[0]] = {"x": float(args[1]), "y": float(args[2])}

        # now generate the crosssections
        for line in surfacelines[1:]:
            crosssection = Crosssection.from_dam_data(line, locations)
            if crosssection is not None:
                result.crosssections.append(crosssection)

        return result

    @property
    def zmin(self) -> float:
        return min([crs.zmin for crs in self.crosssections])

    def _handle_crosssection(
        self,
        crosssection: Crosssection,
        zmin: float,
        weight_length: float,
        offset_from_referenceline: float,
    ):

        result = 0.0
        zmax = crosssection.zmax + 0.1
        ls = np.linspace(
            offset_from_referenceline,
            offset_from_referenceline + int(weight_length),
            int(weight_length),
        )
        crspoints = [[p.l, p.z] for p in crosssection.points]
        crspoints += [[crspoints[-1][0], zmin - 1.0], [crspoints[0][0], zmin - 1.0]]
        crspolygon = Polygon(crspoints)
        for i in range(1, len(ls)):
            lmin = ls[i - 1]
            lmax = ls[i]
            rect = Polygon([(lmin, zmax), (lmax, zmax), (lmax, zmin), (lmin, zmin)])
            ipolygon = rect.intersection(crspolygon)
            result += pow(ipolygon.area, 3)  # * (lmid - lmin) / LENGTH_FOR_WEIGHTS

        return result

    def get_n_normative(
        self,
        output_path: Union[Path, str],
        num_results: int = 3,
        weight_length: float = 20.0,
        offset_from_referenceline: float = 5.0,
    ):
        selection_weighted = []
        fig = Figure(figsize=(15, 6))
        ax = fig.add_subplot()

        f = open(Path(output_path) / "surfacelines_normative.csv", "w")
        f.write("LOCATIONID;X1;Y1;Z1;.....;Xn;Yn;Zn;(Profiel)\n")

        zmin = self.zmin
        for crs in self.crosssections:
            weight = self._handle_crosssection(
                crs, zmin, weight_length, offset_from_referenceline
            )
            xs = [p.l for p in crs.points]
            ys = [p.z for p in crs.points]
            ax.plot(xs, ys, "k:")
            selection_weighted.append((weight, crs))

        # sort on weight
        selection_weighted = sorted(selection_weighted, key=lambda x: x[0])

        for _, crs in selection_weighted[:num_results]:
            xs = [p.l for p in crs.points]
            ys = [p.z for p in crs.points]
            ax.plot(xs, ys, label=f"{crs.id}")
            f.write(f"{crs.to_surfaceline()}\n")

        f.close()
        ax.set_title("normative crosssections")
        fig.legend()
        fig.savefig(Path(output_path) / f"normative_crosssections_result.png")


if __name__ == "__main__":
    cl = CrosssectionList.from_dam_csv(
        "./testdata/surfacelines.csv", "./testdata/locations.csv"
    )
    cl.get_n_normative("./testdata/output")
