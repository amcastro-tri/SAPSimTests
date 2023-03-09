"""Example simulation of a walking Cassie."""

import argparse
from dataclasses import dataclass
from math import sqrt
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from pydrake.common import FindResourceOrThrow
from pydrake.geometry import DrakeVisualizer
from pydrake.geometry import ProximityProperties
from pydrake.geometry import (Box, Cylinder, HalfSpace, Sphere)
from pydrake.math import RigidTransform, RotationMatrix, RollPitchYaw

# Multibody
from pydrake.multibody.meshcat import JointSliders
from pydrake.multibody.plant import (
    AddMultibodyPlantSceneGraph,
    CoulombFriction,
    DiscreteContactSolver,
)
from pydrake.multibody.parsing import Parser
from pydrake.multibody.tree import (
    PrismaticJoint,
    RotationalInertia,
    SpatialInertia,
    UnitInertia,    
)

# Systems
from pydrake.systems.framework import DiagramBuilder
from pydrake.systems.analysis import (PrintSimulatorStatistics, Simulator)
from pydrake.systems.primitives import (
    MatrixGain,    
    ConstantVectorSource,
    Multiplexer)


# Misc.
from pydrake.all import (MeshcatVisualizer, StartMeshcat)
from pydrake.all import (AbstractValue,
                         LeafSystem,
                         Rgba)

def xyz_rpy_deg(xyz, rpy_deg):
    """Shorthand for defining a pose."""
    rpy_deg = np.asarray(rpy_deg)
    return RigidTransform(RollPitchYaw(rpy_deg * np.pi / 180), xyz)

@dataclass
class ContactProperties:
    stiffness: float
    relaxation_time: float
    friction: float
    barrier_distance: float = 1.0
    """Barrier distance in meters. Default corresponds to a large value for
    which the contact potential reduces to the compliant model.
    The barrier of two bodies is the sum of each body's barrier distance. """    

def sphere_volume(r):
    return 4./3.* math.pi * r * r *r

def AddGround(contact_properties, plant):
    properties = ProximityProperties()
    properties.AddProperty(
        "material", "point_contact_stiffness", contact_properties.stiffness)
    properties.AddProperty(
        "material", "relaxation_time", contact_properties.relaxation_time)
    properties.AddProperty(
        "material", "barrier_distance", contact_properties.barrier_distance)        
    properties.AddProperty(
        "material", "coulomb_friction",
        CoulombFriction(
            contact_properties.friction,
            contact_properties.friction))
    plant.RegisterCollisionGeometry(
      plant.world_body(), RigidTransform.Identity(),
      HalfSpace(), "ground_collision", properties)
    plant.RegisterVisualGeometry(
      plant.world_body(), RigidTransform([0, 0, -0.05]),
      Box(5, 5, 0.1), "ground_visual", [0, 0, 1, 1])


def AddSphere(name, mass, radius, contact_properties, color, plant):
    I_Bo = mass * UnitInertia.SolidSphere(radius)
    M_Bo = SpatialInertia.MakeFromCentralInertia(
        mass, np.array([0, 0, 0]), I_Bo)
    body = plant.AddRigidBody(name, M_Bo)

    # Contact properties.
    properties = ProximityProperties()
    properties.AddProperty(
        "material", "point_contact_stiffness", contact_properties.stiffness)
    properties.AddProperty(
        "material", "relaxation_time", contact_properties.relaxation_time)
    properties.AddProperty(
        "material", "barrier_distance", contact_properties.barrier_distance)        
    properties.AddProperty(
        "material", "coulomb_friction",
        CoulombFriction(
            contact_properties.friction,
            contact_properties.friction))
    
    plant.RegisterCollisionGeometry(
      body, RigidTransform.Identity(),
      Sphere(radius), name + "_collision", properties)
    plant.RegisterVisualGeometry(
      body, RigidTransform.Identity(),
      Sphere(radius), name + "_visual", color)
    return body

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target_realtime_rate", type=float, default=1.0,
        help="Desired rate relative to real time.  See documentation for "
             "Simulator::set_target_realtime_rate() for details.")
    parser.add_argument(
        "--simulation_time", type=float, default=1.0,
        help="Desired duration of the simulation in seconds.")
    parser.add_argument(
        "--time_step", type=float, default=0.01,
        help="If greater than zero, the plant is modeled as a system with "
             "discrete updates and period equal to this time_step. "
             "If 0, the plant is modeled as a continuous system.")
    args = parser.parse_args()

    # Start the visualizer.
    meshcat = StartMeshcat()
    meshcat.Set2dRenderMode(xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5)

    # Build model.
    builder = DiagramBuilder()
    plant, scene_graph = AddMultibodyPlantSceneGraph(
        builder=builder, time_step=args.time_step)
    plant.set_discrete_contact_solver(DiscreteContactSolver.kSap)

    g = np.linalg.norm(plant.gravity_field().gravity_vector())
    print(f"g={g}")
    contact_properties = ContactProperties(
        stiffness=1.0e5, relaxation_time=0.01, friction=1.0)

    AddGround(contact_properties, plant)

    # Sphere 1
    density = 1000.0    
    radius1 = 0.02
    mass1 = sphere_volume(radius1) * density
    color = [1,0,0,1]
    print(f"m1 = {mass1}. w1 = {mass1*g}")
    body1 = AddSphere("body1", mass1, radius1, contact_properties, color, plant)

    # Sphere 2
    density = 1000.0    
    radius2 = 0.1
    mass = sphere_volume(radius2) * density
    color = [0,1,0,1]
    print(f"m2 = {mass}. w2 = {mass*g}")
    body2 = AddSphere("body2", mass, radius2, contact_properties, color, plant)

    plant.Finalize()

    # Add viz.
    visualizer = MeshcatVisualizer.AddToBuilder(builder, scene_graph, meshcat)

    # If we wanted the Drake viz.
    #DrakeVisualizer.AddToBuilder(builder, scene_graph)

    # Done defining the diagram.
    diagram = builder.Build()

    # Create context and set initial condition.
    context = diagram.CreateDefaultContext()
    plant_context = plant.GetMyMutableContextFromRoot(context)
    plant.SetFreeBodyPose(
        plant_context, body1, RigidTransform([0, 0, radius1]))
    plant.SetFreeBodyPose(
        plant_context, body2, RigidTransform([-0.000001, 0, 2*radius1+radius2]))

    # Publish initial visualization
    simulator = Simulator(diagram, context)
    simulator.set_target_realtime_rate(args.target_realtime_rate)
    simulator.Initialize()

    input("Press Enter to continue...")

    # Run sim.    
    simulator.AdvanceTo(args.simulation_time)
    PrintSimulatorStatistics(simulator)


if __name__ == "__main__":
    main()
