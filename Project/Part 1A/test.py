import numpy as np
from numba import jit
import sys 
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d.axes3d as p3
#from ast2000solarsystem import AST2000SolarSystem
import random as random


seed = 1

H2_molar_mass = 2.01588

gravitational_constant = 6.67408e-11

avogadro = 6.02214179e23
boltzmann = 1.380648e-23


satellite_mass = 1000


#solar_system = AST2000SolarSystem(seed)


def escape_velocity(M, R):
    v_escape = np.sqrt((2*G*M)/R)
    return v_escape 


class rocket_engine():
    def __init__(self, dimensions = 3, temperature = 10000, N = 1E5, mass = 1.67e-27, length = 1e-6):
        self.k = 1.38064852e-23
        self.T = temperature
        self.N = N
        self.m = mass
        self.L = length 
        self.dim = dimensions
        self.sigma = np.sqrt((self.k*self.T)/self.m)
        self.x = self.position()
        self.v = self.velocities()
        self.steps = 1000
        self.time = 10E-9
        self.dt  = self.time/self.steps
        
    
    def velocities(self):
        return np.random.normal(0,self.sigma, size=(int(self.N),self.dim))

    def position(self):
        return np.random.uniform(0,self.L,   size=(int(self.N),self.dim))


    
    def meanvel(self):
        self.v_s = 0
        for i in range(int(self.N)):
            self.v_s += np.sqrt(self.v[i,0]**2+self.v[i,1]**2+self.v[i,2]**2)
        return self.v_s

    def meankin(self):
        m = self.m
        vel = 0
        for i in range(int(self.N)):
            vel += self.v[i,0]**2 + self.v[i,1]**2 + self.v[i,2]**2
        return 0.5*m*vel

    def test_mean(self):
        """
        making a test function that runs meankin() and meanvel() and checks the
        computed velocity and kinetic energy and the relative error between them
        anything below 1% is perfectly acceptable
        """
        m = self.m
        analytical_mean = 1.5*self.T*self.k
        computed_mean   = 0 
        for j in self.v:
            computed_mean += self.meankin() 
            computed_mean  = computed_mean/self.N
            relative_error =    abs(analytical_mean - computed_mean)/analytical_mean
            print("----------Kinetic energy----------")
            print("{:<20}{:g}".format("Computed mean:", computed_mean))
            print("{:<20}{:g}".format("Analytical mean:", analytical_mean))
            print("{:<20}{:.2f}%".format("Relative error:", relative_error * 100))
            print("-----------------------------")
            break
        assert relative_error < 0.02, "the mean kinetic energy is off"

        print("----------Velocity----------")


        analytical_vel = np.sqrt(8*self.k*self.T/(np.pi*m))
        computed_vel   = 0
        for i in self.v: 
            computed_vel += self.meanvel()
            computed_vel = computed_vel/self.N 
            relative_error = abs(analytical_vel - computed_vel)/analytical_vel
            print("{:<20}{:g}".format("Computed velocity:", computed_vel))
            print("{:<20}{:g}".format("Analytical velocity:", analytical_vel))
            print("{:<20}{:.2f}%".format("Relative error:", relative_error *100))
            print("-----------------------------")
            break
        assert relative_error < 0.02, "the mean velocity is off"


    def box_simulation(self):
        p = 0
        momentum = 0
        escapies = 0
        colliding = 0
        short_side = 0.25*self.L 
        long_side = 0.75*self.L 
        exiting = np.zeros_like(self.x, dtype = bool)
        for t in range(int(self.steps)):
            self.x += self.v * self.dt
            exiting_velocities = abs(self.v[:,2])
            case_1 = np.logical_and(np.greater_equal(self.x, self.L), np.greater(self.v, 0)) # 
            case_2 = np.logical_and(np.less_equal(self.x, 0), np.less(self.v, 0))
            collision = np.logical_and(case_1, case_2)

            case_3 = np.logical_and(np.less_equal(self.x[:,2], 0), np.less(self.v[:,2], 0))
            case_4 = np.logical_and(np.greater_equal(self.x[:,0], short_side), np.less_equal(self.x[:,0],long_side))
            case_5 = np.logical_and(np.greater_equal(self.x[:,1], short_side), np.less_equal(self.x[:,1],long_side))
            collision_2 = np.logical_and(case_3, case_4, case_5)
            exiting[:,0]= collision_2
            exiting[:,1] =collision_2
            exiting[:,2] = collision_2

 
            escapies      += np.sum(collision_2)
            p += np.sum(collision_2*exiting_velocities)

            variable = exiting.astype(np.int)
            spin1 = np.logical_and(collision, np.logical_not(exiting)).astype(np.int)


            colliding  += np.sum(spin1)
            spin2 = spin1.astype(np.int)*2-1
            self.v = -self.v*spin2

        p+= (2 * self.m* p * escapies)/self.dt #force on top wall
        particle_per  = escapies/self.time         #The box force averaged over all time steps
        mean_force = p/self.steps 
        box_mass = particle_per * self.m
        pressure = mean_force/self.L**2

        print("----")
        result1 = print('There are {:g} particles exiting the gas box per second.'\
                        .format(particle_per))
        result2 = print("The average {:g} Thrust per second."\
                        .format(mean_force))
        result3 = print('We would need {:g} Boxes'\
                        .format(box_mass))
        result4 = print("Nummerical pressure on the wall {:g}".format(pressure))
        result5 = print("Analytical pressure on the wall {:g}".format(self.N/self.L**3*self.k*self.T))
        return result1, result2, result3, result4, result5, 


    def test_pressure(self):



    """
    def plot(self):
        T, N, dt, dim, steps = self.T, self.N, self.dt, self.dim, self.steps
        def update_lines( num, dataLines, lines):
            for line, data in zip(lines, dataLines):
                line.set_data(data[0:2, num - 1:num])
                line.set_3d_properties(data[2, num - 1:num])
            return lines

    # Attach 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    frames = 100

    # Run the actual simulation
    x_momentum, datax = self.escape()
    lines = []
    data = []
    for i in range(N):
        data.append([datax[i]])  # wrap data inside another layer of [], needed for animation!
        lines.append([ax.plot(data[i][0][0, 0:1], data[i][0][1, 0:1], data[i][0][2, 0:1], 'o')[0]])

    # Set the axes properties
    ax.set_xlim3d([0.0, dim])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, dim])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, dim])
    ax.set_zlabel('Z')

    ax.set_title('Particle Animation')

    hole = patches.Circle((dim / 2, dim/ 2), dim / 4)
    hole.set_facecolor('black')
    ax.add_patch(hole)
    art3d.pathpatch_2d_to_3d(hole, z=0, zdir="z")

    ani = [i for i in xrange(number_of_particles)]
    for i in ani:
        ani[i] = animation.FuncAnimation(fig, update_lines, frames, fargs=(data[i], lines[i]),
                                         interval=50, blit=False)

    plt.show()



    def debog(self):
        exiting, col_wall, f = escape()
        part_s = exiting/self.time
        test = print(part_s)
        return test
    """
    def maxwell(self, x):
        sigma = np.sqrt(self.k * self.T/self.m)
        exponent = -x**2/(2*sigma**2)
        return 1/(sigma*np.sqrt(2*np.pi))*np.exp(exponent)


"""
x1 = np.linspace(-25000, 25000, 51)
x_label = ["v_x", "v_y", "v_z"]
gas = rocket_engine()
for i, label in enumerate(x_label):
    plt.style.use("classic")
    plt.grid()
    plt.hist(gas.v[:,i], bins=31, density = True, histtype = "step")
    plt.plot(x1, gas.maxwell(x1), "r-")
    plt.show()

    def escape(self):
        exiting = 0
        f = 0
        col_wall = 0
        p_zdirection = 0
        short_side = 0.25*L 
        long_side = 0.75*L 
        condition_exit = np.zeros_like(x, dtype = bool)
        for i in range(int(steps)):
            x += v * dt 
            condition1 = np.logical_and(np.greater_equal(self.x, self.L), np.greater(self.v, 0)) # 
            out_side_box2 = np.logical_and(np.less_equal(self.x, 0), np.less(self.v, 0))
            condition_collision = np.logical_or(condition1,condition2)

            condition3 = np.logical_and(np.less_equal(self.x[:,2], 0), np.less(self.v[:,2], 0))
            condition4 = np.logical_and(np.greater_equal(self.x[:,0], short_side), np.less_equal(self.x[:,0],long_side))
            condition5 = np.logical_and(np.greater_equal(self.x[:,1], short_side), np.less_equal(self.x[:,1],long_side))

            temp = np.logical_and(condition3, condition4, condition5)
            condition_exit[:,0] = temp
            condition_exit[:,1] = temp
            condition_exit[:,2] = temp
            exiting      += np.sum(temp)
            p_zdirection += np.sum(temp*abs(v[:,2]))

            refill_mat = condition_exit.astype(np.int8)
            flip_bool = np.logical_and(condition_collision, np.logical_not(condition_exit)).astype(np.int8)

            x[:,2] = x[:,2] + refill_mat[:,0]*L

            col_wall = col_wall + np.sum(flip_bool)
            #flip = flip_bool.astype(np.int8)*2-1
            v = -v

        f+= (2*m*p_zdirection*exiting)/dt #force on top wall
        particle_per  = exiting/time         #The box force averaged over all time steps
        mean_force = f/steps 
        box_mass = particle_per * m

        print("----")
        result1 = print('There are {:g} particles exiting the gas box per second.'\
                        .format(particle_per))
        result2 = print("The average {:g} Thrust per second."\
                        .format(mean_force))
        result3 = print('We would need {:g} Boxes'\
                        .format(box_mass))
        return result1, result2, result3  
"""

if __name__ == "__main__":
    A = rocket_engine()
    #result2 = A.box_escape()
    #result3 = A.test_mean()
    result3 = A.box_simulation()
    
    

    

