from typing import List, Tuple
import math

class Point:
  pass

class Point:
  def __init__(self, x: float, y:float):
    self.x: float = x
    self.y: float = y

  def add(self, p: Point) -> Point: return Point(self.x + p.x, self.y + p.y)
  def sub(self, p: Point) -> Point: return Point(self.x - p.x, self.y - p.y)
  def mul(self, s: float) -> Point: return Point(self.x * s, self.y * s)
  def div(self, s: float) -> Point: return Point(self.x / s, self.y / s)
  def normalize(self) -> float: 
    if self.len() == 0:
      return self
    else:
      return self.div(self.len())
  def dot(self, p: Point) -> float: return self.x * p.x + self.y * p.y
  def len(self) -> float: return (self.x**2 + self.y**2)**0.5
  def angle(self) -> float: 
    x = self.x
    y = self.y if self.y != 0 else 1e-60
    return math.atan(self.y/self.x)

  @staticmethod 
  def fromAngle(alpha: float) -> Point:
    return Point(math.cos(alpha), math.sin(alpha))

  @staticmethod 
  def xy_list(points: List[Point]) -> Tuple[List[float], List[float]]:
    x_list: List[float] = [p.x for p in points]
    y_list: List[float] = [p.y for p in points]
    return (x_list, y_list)

  def __str__(self) -> str: return "Point(" + str(self.x) + ", " + str(self.y) + ")"

  def __mul__(self, other):
    if isinstance(other, float) or isinstance(other, int):
      return self.mul(other)

  def __truediv__(self, other):
    if isinstance(other, float) or isinstance(other, int):
      return self.div(other)

  def __add__(self, other):
    if isinstance(other, float) or isinstance(other, int):
      return Point(self.x + other, self.y + other)
    if isinstance(other, Point):
      return self.add(other)

  def __sub__(self, other):
    if isinstance(other, float) or isinstance(other, int):
      return Point(self.x - other, self.y - other)
    if isinstance(other, Point):
      return self.sub(other)