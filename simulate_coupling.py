import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from collections import deque
import sys
import random

class Particle:
    def __init__(self, L, v):
        self.r = np.array([np.random.uniform(0,L),np.random.uniform(0,L),0])
        self.r_ = deque(maxlen=50)
        for i in range(10):
            self.r_.append(self.r)
        self.theta = np.random.uniform(0,2*np.pi)
        self.v = v



class System:
    def __init__(self, N, L, R, eta, v, T_0, r_c, dir):
        self.L = L
        self.N = N
        self.R = R
        self.eta = eta
        self.v = v
        self.T_0 = T_0
        self.r_c = r_c
        self.P = []
        self.dir=dir
        self.make_particles()
        self.sampled_t = []

    def make_particles(self):
        for i in range(self.N):
            self.P.append(Particle(L=self.L, v=self.v))

    def diff_mat(self):
        mat = np.zeros((self.N,self.N,3))
        for i in range(self.N):
            for j in range(self.N): # Comp-time can be halved here
                mat[i,j,:] = self.P[i].r-self.P[j].r

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

    def update_p(self, n):
        theta = self.P[n].theta
        v_vec = np.array([np.cos(theta), np.sin(theta), 0])
        e_z = np.array([0,0,1])
        T_n = 0
        r = self.P[n].r
        if self.T_0 != 0:
            for i in range(self.N): # Do mod L somewhere here
                if i != n:
                    r_ni_ = (self.P[i].r-r)
                    '''r_ni = np.zeros(3)
                    for j in [0,1]:
                        if r_ni_[j] > self.L/2:
                            r_ni[j] = r_ni_[j]-self.L
                        elif r_ni_[j] < -self.L/2:
                            r_ni[j] = self.L + r_ni_[j]
                        else:
                            r_ni[j] = r_ni_[j]'''
                    r_ni = self.periodic_distance(r_ni_)

                    if np.linalg.norm(r_ni) < self.r_c:
                        T_n += self.T_0 * np.dot(v_vec,r_ni)/np.linalg.norm(r_ni)**2*np.dot(np.cross(v_vec,r_ni), e_z)


        theta = self.P[n].theta
        xi = np.random.uniform(-self.eta,self.eta)/2
        theta += T_n + xi

        v = self.P[n].v
        self.P[n].r += v*np.array([np.cos(theta), np.sin(theta), 0])
        self.P[n].r = self.P[n].r%(self.L)

        for i in range(self.N):  # Do mod L somewhere here
            if i != n:
                delta_cm = self.periodic_distance(self.P[n].r-self.P[i].r)
                delta_cm_norm = np.linalg.norm(delta_cm)
                if delta_cm_norm < 2*self.R:
                    delta_cm_normed = delta_cm/delta_cm_norm
                    overlap = 2*self.R - delta_cm_norm
                    self.P[n].r += delta_cm_normed * 1 / 2 * overlap
                    self.P[i].r -= delta_cm_normed * 1 / 2 * overlap

                    self.P[n].r_[-1] = self.P[n].r
                    self.P[i].r_[-1] = self.P[i].r


        self.P[n].r = self.P[n].r % (self.L)
        self.P[n].r_.append(np.copy(self.P[n].r))
        self.P[n].theta = theta








    def time_step(self):
        for i in range(self.N):
            self.update_p(i)

    def iterate(self, T, plot_freq):
        t = 0
        while t < T:
            self.time_step()
            if t%plot_freq == 0:
                self.gen_state_fig(t)
                print('Fin iter %i'%t)
            t += 1

    def gen_state_fig(self, t):

        edges = {}
        for i in range(self.N):
            first = True
            for j in range(0,self.N):   #(i+1,self.N):
                if np.linalg.norm(self.periodic_distance(self.P[i].r - self.P[j].r)) < 2.2*self.R:
                    if first:
                        first = False
                        edges[i] = set()
                        edges[i].add(j)
                    edges[i].add(j)

        assignment = {}

        def recursion(key, cluster):
            if assignment.__contains__(key): return
            assignment[key] = cluster
            if not edges.__contains__(key): return
            connections = edges.pop(key)
            for con in list(connections):
                recursion(con, cluster)

        cluster = 0

        while bool(edges):
            k = list(edges.keys())[0]
            recursion(k, cluster)
            cluster += 1

        clusters = {}
        for key in list(assignment):
            if not clusters.__contains__(assignment[key]):
                clusters[assignment[key]] = 1
            else:
                clusters[assignment[key]] += 1

        cc = {2:'red', 3:'orange', 4:'lime', 5:'cyan', 6:'magenta', 7:'darkviolet', 8:'lightpink'}
        for j in range(9,100):
            cc[j] = 'black'

        dir = self.dir + '/t%i'%t
        self.sampled_t.append(t)
        plt.figure()
        ax = plt.gca()
        for p_ix in range(len(self.P)):
            p = self.P[p_ix]
            c = 'r'
            fill = False
            if assignment.__contains__(p_ix):
                if clusters[assignment[p_ix]] > 1:
                    c = cc[clusters[assignment[p_ix]]]
                    fill = True
            circle = plt.Circle((p.r[0],p.r[1]), 1, color=c, fill=fill)
            ax.add_artist(circle)
            for i in range(len(p.r_)-1,1,-2):
                r_f = p.r_[i]
                r_i = p.r_[i-2]
                if np.linalg.norm(r_f-r_i) > self.L/2:
                    break
                plt.plot([r_f[0],r_i[0]],[r_f[1],r_i[1]],color=c)

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

def main():
    N=100
    L=100
    R=1
    eta=2*np.pi
    v=.05
    T_0=1
    r_c=6
    dir = 'coupling_figs'

    system = System(N=N,L=L,R=R,eta=eta,v=v,T_0=T_0, r_c=r_c, dir=dir)
    system.iterate(T=100, plot_freq=20)
    system.plot_grid_time(dir, eta)

if __name__ == "__main__":
    main()