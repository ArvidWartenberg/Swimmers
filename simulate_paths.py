import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from collections import deque
import sys
import random

class Particle:
    def __init__(self, L, v, maxlen=0):
        self.r = np.array([np.random.uniform(0,L),np.random.uniform(0,L),0])
        self.r_ = deque(maxlen=maxlen)
        for i in range(10):
            self.r_.append(self.r)
        self.theta = np.random.uniform(0,2*np.pi)
        self.v = v



class System:
    def __init__(self, N, M, L, R, eta, v, T_0, r_c, dir):
        self.L = L
        self.N = N
        self.M = M
        self.R = R
        self.eta = eta
        self.v = v
        self.T_0 = T_0
        self.r_c = r_c
        self.P_A = []
        self.P_P = []
        self.dir=dir
        self.make_particles()
        self.sampled_t = []
        self.fix_overlap(reps=2)

    def make_particles(self):
        for i in range(self.N):
            self.P_A.append(Particle(L=self.L, v=self.v, maxlen=10000))
        for i in range(self.M):
            self.P_P.append(Particle(L=self.L, v=self.v))

    def diff_mat(self):
        mat = np.zeros((self.N,self.N,3))
        for i in range(self.N):
            for j in range(self.N): # Comp-time can be halved here
                mat[i,j,:] = self.P_A[i].r-self.P_A[j].r

    def periodic_distance(self, r_ni_):
        r_ni = np.zeros(3)
        for j in [0, 1]:
            if r_ni_[j] > self.L / 2:
                r_ni[j] = r_ni_[j] - self.L
            elif r_ni_[j] < -self.L / 2:
                r_ni[j] = self.L + r_ni_[j]
            else:
                r_ni[j] = r_ni_[j]
        return r_ni

    def update_p_P(self, n):
        r = self.P_P[n].r
        theta = self.P_P[n].theta
        xi = np.random.uniform(-np.pi,np.pi)
        theta += xi

        v = self.P_P[n].v
        self.P_P[n].r += 0.1*v*np.array([np.cos(theta), np.sin(theta), 0])
        self.P_P[n].r = self.P_P[n].r%(self.L)

        self.P_P[n].r = self.P_P[n].r % (self.L)
        self.P_P[n].r_.append(np.copy(self.P_P[n].r))
        self.P_P[n].theta = theta


    def update_p_A(self, n):
        theta = self.P_A[n].theta
        v_vec = np.array([np.cos(theta), np.sin(theta), 0])
        e_z = np.array([0,0,1])
        T_n = 0
        r = self.P_A[n].r
        if self.T_0 != 0:
            for i in range(self.N): # Do mod L somewhere here
                if i != n:
                    r_ni_ = (self.P_A[i].r-r)
                    r_ni = self.periodic_distance(r_ni_)

                    if np.linalg.norm(r_ni) < self.r_c:
                        T_n += self.T_0 * np.dot(v_vec,r_ni)/np.linalg.norm(r_ni)**2*np.dot(np.cross(v_vec,r_ni), e_z)

            for i in range(self.M): # Do mod L somewhere here
                r_ni_ = (self.P_P[i].r-r)
                r_ni = self.periodic_distance(r_ni_)

                if np.linalg.norm(r_ni) < self.r_c:
                    T_n -= self.T_0 * np.dot(v_vec,r_ni)/np.linalg.norm(r_ni)**2*np.dot(np.cross(v_vec,r_ni), e_z)

        theta = self.P_A[n].theta
        xi = np.random.uniform(-self.eta,self.eta)/2
        theta += T_n + xi

        v = self.P_A[n].v
        self.P_A[n].r += v*np.array([np.cos(theta), np.sin(theta), 0])
        self.P_A[n].r = self.P_A[n].r%(self.L)

        self.P_A[n].r = self.P_A[n].r % (self.L)
        self.P_A[n].r_.append(np.copy(self.P_A[n].r))
        self.P_A[n].theta = theta





    def fix_overlap(self,reps):
        for n in range(self.N):
            for i in range(self.N):  # Do mod L somewhere here
                if i != n:
                    delta_cm = self.periodic_distance(self.P_A[n].r-self.P_A[i].r)
                    delta_cm_norm = np.linalg.norm(delta_cm)
                    if delta_cm_norm < 2*self.R:
                        delta_cm_normed = delta_cm/delta_cm_norm
                        overlap = 2*self.R - delta_cm_norm
                        self.P_A[n].r += delta_cm_normed * 1 / 2 * overlap
                        self.P_A[i].r -= delta_cm_normed * 1 / 2 * overlap

                        self.P_A[n].r_[-1] = self.P_A[n].r
                        self.P_A[i].r_[-1] = self.P_A[i].r

            for i in range(self.M):  # Do mod L somewhere here
                    delta_cm = self.periodic_distance(self.P_A[n].r - self.P_P[i].r)
                    delta_cm_norm = np.linalg.norm(delta_cm)
                    if delta_cm_norm < 2 * self.R:
                        delta_cm_normed = delta_cm / delta_cm_norm
                        overlap = 2 * self.R - delta_cm_norm
                        self.P_A[n].r += delta_cm_normed * 1 / 2 * overlap
                        self.P_P[i].r -= delta_cm_normed * 1 / 2 * overlap

                        self.P_A[n].r_[-1] = self.P_A[n].r

        for m in range(self.M):
            for i in range(self.M):  # Do mod L somewhere here
                if i != m:
                    delta_cm = self.periodic_distance(self.P_P[m].r - self.P_P[i].r)
                    delta_cm_norm = np.linalg.norm(delta_cm)
                    if delta_cm_norm < 2 * self.R:
                        delta_cm_normed = delta_cm / delta_cm_norm
                        overlap = 2 * self.R - delta_cm_norm
                        self.P_P[m].r += delta_cm_normed * 1 / 2 * overlap
                        self.P_P[i].r -= delta_cm_normed * 1 / 2 * overlap

    def time_step(self):
        for i in range(self.N):
            self.update_p_A(i)

        for i in range(self.M):
            self.update_p_P(i)


    def iterate(self, T, plot_freq):
        t = 0
        while t < T:
            self.time_step()
            if t%plot_freq == 0:
                self.fix_overlap(reps=2)
                self.gen_state_fig(t)
                print('Fin iter %i'%t)
            t += 1

    def gen_state_fig(self, t):
        dir = self.dir + '/t%i'%t
        self.sampled_t.append(t)
        plt.figure()
        ax = plt.gca()
        c = 'gray'
        fill = True
        for p_ix in range(len(self.P_P)):
            p = self.P_P[p_ix]
            circle = plt.Circle((p.r[0],p.r[1]), 1, color=c, fill=True, alpha=0.2)
            ax.add_artist(circle)
            circle2 = plt.Circle((p.r[0],p.r[1]), 1, color=c, fill=False, alpha=1)
            ax.add_artist(circle2)


        c = 'r'
        fill = False
        for p_ix in range(len(self.P_A)):
            p = self.P_A[p_ix]
            for i in range(len(p.r_)-1,1,-2):
                r_f = p.r_[i]
                r_i = p.r_[i-2]
                if np.linalg.norm(r_f-r_i) > self.L/2:
                    continue
                plt.plot([r_f[0],r_i[0]],[r_f[1],r_i[1]],color='b', alpha=.03, linewidth=3)
            circle = plt.Circle((p.r[0],p.r[1]), 1, color=c, fill=True)
            ax.add_artist(circle)

        plt.axis([0, self.L, 0, self.L])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.axis('off')
        plt.savefig(dir, dpi=100)
        plt.close()

    def plot_grid_time(self, dir, eta):
        fig = plt.figure()

        while True:
            plt.ion()

            for t in self.sampled_t:
                file = self.dir + '/t' + str(t) + '.png'
                plt.clf()
                image = img.imread(file)
                plt.imshow(image)
                plt.axis('off')
                plt.title('$\\eta$=%.3f, t=%i'%(eta,t), fontsize=15)
                plt.pause(0.005)
                plt.draw()
N=10
M=400
L=50
R=1
eta=2*np.pi
v=.1
T_0=1
r_c=4
dir = 'path_figs'

system = System(N=N, M=M, L=L,R=R,eta=eta,v=v,T_0=T_0, r_c=r_c, dir=dir)
system.iterate(T=100, plot_freq=10)
system.plot_grid_time(dir, eta)

print()
 