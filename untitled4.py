# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 22:16:04 2021

@author: JOSE
"""

import oommfc as oc
import discretisedfield as df
import micromagneticmodel as mm

region = df.Region(p1=(-40e-9, -40e-9, 0), p2=(60e-9, 40e-9, 1e-9))
mesh = df.Mesh(region=region, cell=(1e-9, 1e-9, 1e-9))

gamma0 = 2.211e5  # gyromagnetic ratio (m/As)
alpha = 0.3  # Gilbert damping

mesh.k3d()

system = mm.System(name='skyrmion')

system.energy = (mm.Exchange(A=1.5e-11)
               + mm.DMI(D=3e-3, crystalclass='Cnv')
               + mm.UniaxialAnisotropy(K=0.8e6, u=(0, 0, 1))
               + mm.Demag() )
system.dynamics = mm.Precession(gamma0=2.211e5) + mm.Damping(alpha=alpha)


Ms = 580e3


def Ms_fun(pos):
    """Function to set magnitude of magnetisation: zero outside cylindric shape,
    Ms inside cylinder.

    Cylinder radius is 50nm.

    """
    x, y, z = pos
    if (x**2 + y**2)**0.5 < 40e-9:
        return Ms
    else:
        return 0


def m_init(pos):
    """Function to set initial magnetisation direction:
    -z inside cylinder (r=10nm),
    +z outside cylinder.
    y-component to break symmetry.

    """
    x, y, z = pos
    if (x**2 + y**2)**0.5 < 8e-9:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


# create system with above geometry and initial magnetisation
system.m = df.Field(mesh, dim=3, value=m_init, norm=Ms)

system.m.norm.k3d_nonzero()

system.m.plane('z').mpl()

# minimize the energy
md = oc.MinDriver()
md.drive(system)

# Plot relaxed configuration: vectors in z-plane
system.m.plane('z').mpl()


ux = -39.91924  # velocity in x-direction (m/s)
beta = 0.6  # non-adiabatic STT parameter

system.dynamics += mm.ZhangLi(u=ux, beta=beta)  # please notice the use of `+=` operator
td = oc.TimeDriver()
td.drive(system, t=0.6e-9, n=100)

system.m.z.plane('z').k3d_scalar()


