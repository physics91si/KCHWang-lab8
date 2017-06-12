#!/usr/bin/python3

# Physics 91SI
# molecule 2017
# Lab 8

# Modules you won't need
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Modules you will need
import numpy as np
from particle import Particle, Molecule

# TODO: Implement this function
def init_molecule():
    """Create Particles p1 and p2 inside boundaries and return a molecule
    connecting them"""
    p1 = Particle(np.array([0.2, 0.2]), 1)
    p2 = Particle(np.array([0.8, 0.8]), 2)
    spring_cos = 1
    equil_len = 0.5
    return Molecule(p1.pos, p2.pos, p1.m, p2.m, spring_cos, equil_len)


# TODO: Implement this function
def time_step(dt, mol):
    """Sets new positions and velocities of the particles attached to mol"""
    F = mol.get_force()
    delta_v_1 = F / mol.p1.m * dt
    delta_v_2 = -F / mol.p2.m * dt
    delta_pos_1 = mol.p1.vel * dt
    delta_pos_2 = mol.p2.vel * dt
    
    mol.p1.vel += delta_v_1
    mol.p2.vel += delta_v_2
    mol.p1.pos += delta_pos_1
    mol.p2.pos += delta_pos_2
    
    print('p1_pos:', mol.p1.pos[0], mol.p1.pos[1], 'p2_pos:', mol.p2.pos[0], mol.p2.pos[1])


#############################################
# The rest of the file is already implemented
#############################################

def run_dynamics(n, dt, xlim=(0, 1), ylim=(0, 1)):
    """Calculate each successive time step and animate it"""
    mol = init_molecule()

    # Animation stuff
    fig, ax = plt.subplots()
    line, = ax.plot((mol.p1.pos[0], mol.p2.pos[0]), (mol.p1.pos[1], mol.p2.pos[1]), '-o')
    
    ax.clear()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.title('Dynamics simulation')
    dynamic_ani = animation.FuncAnimation(fig, update_anim, n,
            fargs=(dt, mol,line), interval=50, blit=False)
    plt.show()

def update_anim(i,dt, mol,line):
    """Update and draw the molecule. Called by FuncAnimation"""
    time_step(dt, mol)
    line.set_data([(mol.p1.pos[0], mol.p2.pos[0]),
                   (mol.p1.pos[1], mol.p2.pos[1])])
    return line,

if __name__ == '__main__':
    # Set the number of iterations and time step size
    n = 10
    dt = .1
    run_dynamics(n, dt)
