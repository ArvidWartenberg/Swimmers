import numpy as np
import matplotlib.pyplot as plt
import sys
import random

class Particle:
    def __init__(self, x=0, y=0, phi=0, v=1):
        self.x = x
        self.y = y
        self.phi = phi
        self.v = v
        self.X = [x]
        self.Y = [y]
        self.PHI = [phi]
        self.V = [v]

    def update(self, dx, dy, dphi, dv=0):
        self.x += dx
        self.y += dy
        self.phi += dphi
        self.v += dv
        self.X.append(self.x)
        self.Y.append(self.y)
        self.PHI.append(self.phi)
        self.V.append(self.v)



class System:
    def __init__(self, DT, DR, dt, N, V=None):
        self.DT = DT
        self.DR = DR
        self.dt = dt
        self.N = N
        self.particles = []
        self.make_particles(V)

    def make_particles(self, V=None):
        for i in range(self.N):
            if V == None:
                self.particles.append(Particle(phi=np.random.uniform(0,2*np.pi)))
            else:
                self.particles.append(Particle(phi=np.random.uniform(0,2*np.pi), v=V[i]))


    def time_step(self):
        for p in self.particles:
            dx = (p.v*np.cos(p.phi) + np.sqrt(2*self.DT)*np.random.normal(0, 1))*self.dt
            dy = (p.v*np.sin(p.phi) + np.sqrt(2*self.DT)*np.random.normal(0, 1))*self.dt
            dphi = np.sqrt(2*self.DR)*np.random.normal(0, 1)*self.dt
            p.update(dx=dx,dy=dy,dphi=dphi)


    def iterate(self, T):

        t = 0
        while t < T:
            self.time_step()
            t += 1


    def plot(self):

        C_arr = []
        plt.subplot(1,2,1)
        for p in self.particles:
            r = random.random()
            b = random.random()
            g = random.random()
            color = np.array([r, g, b])
            c = color.reshape(1, -1)
            plt.scatter(p.X,p.Y, s=4, alpha=.5, c=c)
            plt.plot(p.X, p.Y, linewidth=1, alpha=.5, color=c.reshape(3))
            C_arr.append(c.reshape(3))
            plt.scatter(p.x,p.y, s=100, alpha=1, c=c)

        Tau = np.arange(1, 1000, 10)
        MSD = system.MSD(Tau)
        plt.tight_layout(True)
        plt.title('Particle trajectories', fontsize=15)
        plt.subplot(1,2,2)
        ax = plt.gca()
        plt.plot(Tau * dt, MSD[0, :], label='$v=$%.1f [$\\mu$ms$^{-1}$]' % system.particles[0].v, color=C_arr[0])
        plt.plot(Tau * dt, MSD[1, :], label='$v=$%.1f [$\\mu$ms$^{-1}$]' % system.particles[1].v, color=C_arr[1])
        plt.plot(Tau * dt, MSD[2, :], label='$v=$%.1f [$\\mu$ms$^{-1}$]' % system.particles[2].v, color=C_arr[2])
        plt.plot(Tau * dt, MSD[3, :], label='$v=$%.1f [$\\mu$ms$^{-1}$]' % system.particles[3].v, color=C_arr[3])
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('$\\tau$ [s]', fontsize=13)
        plt.ylabel('MSD($\\tau$) [$\\mu$m$^2$]', fontsize=13)
        plt.title('MSD vs. $\\tau$ for different $v$', fontsize=15)
        plt.legend()
        plt.tight_layout(True)
        plt.show()

    def MSD(self, Tau):
        MSD = np.zeros((len(self.particles), len(Tau)))
        for p_ix in range(len(self.particles)):
            p = self.particles[p_ix]
            for i in range(len(Tau)):
                tau = Tau[i]
                msd = 0
                T = len(p.X)-tau-1
                for t in range(T):
                    msd += (p.X[t+tau]-p.X[t])**2 + (p.Y[t+tau]-p.Y[t])**2
                MSD[p_ix,i] += msd/T
        return MSD

DT, DR, dt = .22, .16, .1
N=4
V=[0,1,2,3]
system = System(N=N,DT=DT,DR=DR,dt=dt,V=V)
system.iterate(2000)
system.plot()
