from typing import List, Tuple, Any, Optional

import random
import math
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation

from point import Point
from point_vector import PointVector

PointList = List[Point]

fig = plt.figure()
subplt = plt.subplot()
subplt.axis('scaled')
subplt.set_xlim([-6, 6])
subplt.set_ylim([-4, 14])
lineplot, = subplt.plot([], [], 'k-')

class Bubble:
  def __init__(self, 
      npoints: int, pos: Point, 
      k: float = 1.0, mass: float = 1.0, field_strength=900, radius: float = 1.0, g: float = 20.0, gamma: float = 0.1):
    self.k = k # spring constant
    self.g = g # grav constant
    self.gamma = gamma # damping factor
    self.mass = mass
    self.point_mass = mass / npoints
    self.points: PointList = []
    connections = list(range(1, npoints)) + [0]
    self.connections: List[int] = connections
    self.velocities: PointList = [Point(0, 0) for n in range(0, npoints)]
    self.accels: PointList = [Point(0, 0) for n in range(0, npoints)]
    self.pres_accels: PointList = [Point(0, 0) for n in range(0, npoints)]
    self.field_strength = field_strength*self.point_mass

    phi = 2*np.pi/(npoints*2)
    A = np.pi * radius ** 2
    self.internal_energy = A/(np.cos(phi)/np.sin(phi/2)*1/(cV*4)/k)

    angles = np.arange(0, 2*np.pi, 2*np.pi/npoints)

    for alpha in angles:
      point = Point.fromAngle(alpha).mul(radius).add(pos)
      self.points.append(point)

  def get_connection_xy_list(self) -> Tuple[List[float], List[float]]:
    x_list = []
    y_list = []
    for i, j in enumerate(self.connections):
      p1= self.points[i]
      p2 = self.points[j]
      x_list.append(p1.x); x_list.append(p2.x);
      y_list.append(p1.y); y_list.append(p2.y);
    return (x_list, y_list)

  def get_dE_as_position(self, p_acc, vel, acc):
    p_F = p_acc.mul(self.point_mass)
    dEvector = p_F.dot( vel.mul(dt).add( acc.mul(dt*dt*0.5) ) )
    return -1*sum(dEvector)

  def verlet_step(self):
    def PV(x): return PointVector(x)
    E_0 = self.internal_energy
    pos_0 = PV(self.points)
    vel_0 = PV(self.velocities)
    acc_0 = PV(self.accels)
    pres_acc_0 = PV(self.pres_accels)

    new_pos = pos_0.add( vel_0.mul(dt) ).add( acc_0.mul(dt*dt*0.5) )
    dE_pos = self.get_dE_as_position(pres_acc_0, vel_0, acc_0)
    new_E = E_0 + dE_pos

    force = PV(self.point_forces(new_pos, vel_0, new_E))
    pres_force = PV(self.pressure_forces(new_pos, new_E))

    new_acc = force.div(self.point_mass)
    new_pres_acc = pres_force.div(self.point_mass)

    new_vel = vel_0.add( (acc_0.add(new_acc)).mul(dt*0.5) )

    self.points = new_pos
    self.accels = new_acc
    self.pres_accels = new_pres_acc
    self.velocities = new_vel
    self.internal_energy = new_E

  def rk_prim(self, S):
    E_index = 0
    X_index = 1
    V_index = len(self.points) + 1
    E = S[E_index].x
    pos = S[X_index:V_index]
    vel = S[V_index:]
    forces = PointVector(self.point_forces(pos, vel, E))
    pressure_forces = PointVector(self.pressure_forces(pos, E))
    dE = -1*sum(pressure_forces.dot(vel))
    acc = forces.div(self.point_mass)
    return PointVector([Point(dE, 0)] + vel[:] + acc[:])

  def runge_kutta_step(self):
    E_index = 0
    X_index = 1
    V_index = len(self.points) + 1

    E = self.internal_energy
    pos = self.points
    vel = self.velocities

    S = PointVector([Point(E, 0)] + pos[:] + vel[:])

    K1 = self.rk_prim(S                   ).mul(dt)
    K2 = self.rk_prim(S.add( K1.mul(0.5)) ).mul(dt)
    K3 = self.rk_prim(S.add( K2.mul(0.5)) ).mul(dt)
    K4 = self.rk_prim(S.add( K3         ) ).mul(dt)

    dS = K1.add( K2.mul(2) ).add( K3.mul(2) ).add(K4).div(6)
    new_S = S.add(dS)

    self.internal_energy = new_S[E_index].x
    self.points = new_S[X_index:V_index]
    self.velocities = new_S[V_index:]

  def point_forces(self, points: PointList, velocities: PointList, IE: float) -> PointList:
    spring_forces = self.spring_forces(points)
    pressure_forces = self.pressure_forces(points, IE)
    gravity_forces = self.gravity_forces(points)
    contact_forces = self.contact_forces(points, velocities)
    damping_forces = self.damping_forces(velocities)
    forces = [spring_forces[i]
        .add(pressure_forces[i])
        .add(gravity_forces[i])
        .add(damping_forces[i])
        .add(contact_forces[i])
        for i in range(0, len(points))]
    return forces

  def damping_forces(self, velocities: PointList) -> PointList:
    center_vel = self.center_velocity(velocities)
    center_damp = center_vel.mul(0) #-self.gamma*0.2)
    forces: PointList = [v.sub(center_vel).mul(-self.gamma*self.point_mass).add(center_damp) for v in velocities]
    return forces

  def gravity_forces(self, points: PointList) -> PointList:
    forces: PointList = [Point(0, -self.g * self.point_mass) for p in points]
    return forces

  def contact_forces(self, points: PointList, vels: PointList) -> PointList:
    forces: PointList = [Point(0, 0) for p in points]
    for i, p in enumerate(points):
      if p.y < 0:
        forces[i] = Point(0, self.field_strength) #Point(0, 10 * (1 - np.e**p.y))
    return forces

  def spring_forces(self, points: PointList) -> PointList:
    forces: PointList = [Point(0, 0) for p in points]
    for i, j in enumerate(self.connections):
      p_i: Point = points[i]
      p_j: Point = points[j]
      distance_vector: Point = p_j.sub(p_i)
      forces[i] = forces[i].add(distance_vector.mul(self.k))
      forces[j] = forces[j].add(distance_vector.mul(-self.k))
    return forces

  def pressure_forces(self, points: PointList, IE: float) -> PointList:
    center = self.center(points)
    forces: PointList = [Point(0, 0) for p in points]
    area = self.area(points)
    for i, j in enumerate(self.connections):
      p_i: Point = points[i]
      p_j: Point = points[j]
      distance_vec: Point = p_j.sub(p_i)
      distance: float = p_i.sub(p_j).len()
      midpoint: Point = p_i.add(p_j).div(2)
      normal_vector: Point = Point(distance_vec.y, -1*distance_vec.x).normalize()
      NkbT = IE / cV
      force = normal_vector.mul(distance * NkbT / area)
      forces[i] = forces[i].add(force.div(2))
      forces[j] = forces[j].add(force.div(2))
    return forces

  def area(self, points = None) -> float:
    if points is None: points = self.points
    x_list, y_list = Point.xy_list(points)
    return 0.5*np.abs(np.dot(x_list,np.roll(y_list,1))-np.dot(y_list,np.roll(x_list,1)))

  def center(self, points = None) -> Point:
    if points is None: points = self.points
    center: Point = Point(0, 0)
    for p in points: center = center.add(p)
    return center.div(len(self.points))

  def center_velocity(self, velocities = None) -> Point:
    if velocities is None: velocities = self.velocities
    center_vel: Point = Point(0, 0)
    for v in velocities: center_vel = center_vel.add(v)
    return center_vel.div(len(self.points))

  def kinetic_energy(self, vels = None) -> float:
    if vels == None: vels = self.velocities
    total_kE: float = 0
    for v in vels:
      total_kE += 0.5 * self.point_mass * v.len()**2
    return total_kE

  def spring_energy(self, points = None) -> float:
    if points == None: points = self.points
    total_pE: float = 0
    for i, j in enumerate(self.connections):
      p_i: Point = points[i]
      p_j: Point = points[j]
      distance: float = p_j.sub(p_i).len()
      total_pE += 0.5 * self.k * distance**2
    return total_pE

  def grav_energy(self, points = None) -> float:
    if points == None: points = self.points
    total_gE:float = 0
    for p in points:
      total_gE += self.point_mass * self.g * p.y
      total_gE += -p.y * self.field_strength if p.y < 0 else 0
    return total_gE

dt: float = 0.001
cV = 5/2
t = 0

steps_per_frame = 10
numsteps = 3000
numframes = int(numsteps/steps_per_frame)
frame = 0

GRAVITY_DROP = True
IMPART_VELOCITY = False
gravity_frame = 0.1 if GRAVITY_DROP else 2
velocity_frame = 0.1 if IMPART_VELOCITY else 2
start_frame = 0.1

def make_bubble():
  return Bubble(30, Point(0, 10), k = 1250, radius = 2, g = 0, gamma = 15)#1/4)

bubble = make_bubble()

plt.title("Creating water balloons...")

plt.axhline(y = 0, linestyle = "--", c="gray")

def animate(framenr):
  print(framenr*steps_per_frame)
  global t, ts, Es, frame
  for i in range(0, steps_per_frame):
    bubble.runge_kutta_step()
    #bubble.verlet_step()
    if frame >= int(numframes * start_frame) + 1:
      t += dt
  if frame == int(numframes * gravity_frame):
    bubble.g = 20
    #bubble2.g = 20
    plt.title("Adding gravity...")
  if frame == int(numframes * velocity_frame):
    c = bubble.center()
    for i, p in enumerate(bubble.points):
      dir = p.sub(c).normalize()
      if isinstance(bubble.velocities, PointVector):
        bubble.velocities = bubble.velocities.points
      bubble.velocities[i] = dir.mul(10)
    plt.title("Adding gravity...")

  frame += 1

  x_list, y_list = bubble.get_connection_xy_list()
  lineplot.set_data(x_list, y_list)

  if frame < int(numframes * start_frame): return

def init():
  x_list, y_list = bubble.get_connection_xy_list()
  lineplot.set_data(x_list, y_list)
  #x_list, y_list = bubble2.get_connection_xy_list()
  #lineplot2.set_data(x_list, y_list)
  return lineplot,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=numframes, interval=10, blit=False, repeat=False)

plt.rcParams.update({'font.size': 20})
plt.show()
plt.clf()