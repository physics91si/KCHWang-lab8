# Physics 91SI
# Spring 2017
# Lab 8

import numpy as np

class Particle:
    """Stores information about a particle with mass, position, and velocity."""
    
    def __init__(self, Position, M):
        """Create a particle with position (numpy array of len 2) and mass."""
        self.pos = Position   # Sets x position
        self.m = M  # Sets mass
        # Initial velocity and acceleration set to be zero
        self.vel = np.zeros((2,))
        self.acc = np.zeros((2,))
        
        
class Molecule:
    """Contains information about the entire two-particle system."""
    
    def __init__(self, init_pos_1, init_pos_2, M_1, M_2, spring_cos, equi_len):
        """Creates two particles with position and mass. 
        Sets up spring constant and equilibrium length as k and L0."""
        self.p1 = Particle(init_pos_1, M_1)
        self.p2 = Particle(init_pos_2, M_2)
        self.k = spring_cos
        self.L0 = equi_len
        
    def get_displ(self):
        """Returns the vector pointing from particle p1 to particle p2"""
        return self.p2.pos - self.p1.pos
    
    def get_force(self):
        """Returns the force (vector) due to the spring on particle p1"""
        displ = self.get_displ()
        equil = displ / np.linalg.norm(displ) * self.L0
        return self.k * (displ - equil)
        
        
        
        
if __name__ == '__main__':
    p = Particle(0, 1)
    print('initial position:', p.pos)
    print('initial mass:', p.m)
    p.pos = -1
    p.m = 2
    print('updated position:', p.pos)
    print('updated mass:', p.m)
    
    p1 = Particle(np.array([0.2, 0.2]), 1)
    p2 = Particle(np.array([0.8, 0.8]), 1)
    spring_cos = 5
    equil_len = 0.5
    molecule = Molecule(p1.pos, p2.pos, p1.m, p2.m, spring_cos, equil_len)
    print('equil_len:', equil_len, 'force:', molecule.get_force())
    
    equil_len = 1
    molecule = Molecule(p1.pos, p2.pos, p1.m, p2.m, spring_cos, equil_len)
    print('equil_len:', equil_len, 'force:', molecule.get_force())
    
    equil_len = 1.5
    molecule = Molecule(p1.pos, p2.pos, p1.m, p2.m, spring_cos, equil_len)
    print('equil_len:', equil_len, 'force:', molecule.get_force())
    