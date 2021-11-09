from typing import List, Tuple
import math

from point import Point

class PointVector:
  def __init__(self, arg):
    self.points = []
    if isinstance(arg, int):
      for i in range(0, arg):
        self.points.append(Point(0, 0))
    else:
      for point in arg:
        self.points.append(point)

  def add(self, ps: PointVector) -> PointVector: return PointVector([p.add(ps[i]) for i, p in enumerate(self.points)])
  def sub(self, ps: PointVector) -> PointVector: return PointVector([p.sub(ps[i]) for i, p in enumerate(self.points)]) 
  def dot(self, ps: PointVector) -> List[float]: return [p.dot(ps[i]) for i, p in enumerate(self.points)]
  def mul(self, s: float) -> PointVector: return PointVector([p.mul(s) for p in self.points])
  def div(self, s: float) -> PointVector: return PointVector([p.div(s) for p in self.points])
  def normalize(self) -> PointVector: return PointVector([p.normalize() for p in self.points])

  def __len__(self):
    return len(self.points)

  def __iter__(self):
    return self.points.__iter__()

  def __str__(self) -> str: 
    return "[" + ", ".join(["Point(" + str(p.x) + ", " + str(p.y) + ")" for p in self.points]) + "]"

  def __getitem__(self, x):
    return self.points[x]