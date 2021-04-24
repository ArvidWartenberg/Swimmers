import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from collections import deque
import sys
import random


class Particle:
    def __init__(self, L, v):
        self.r = np.array([np.random.uniform(0, L), np.random.uniform(0, L), 0])
        self.r_ = deque(maxlen=30)
        for i in range(10):
            self.r_.append(self.r)
        self.theta = np.random.uniform(0, 2 * np.pi)
        self.v = v


class System:
    def __init__(self, N, L, R, eta, v, T_0, r_c):
        self.L = L
        self.N = N
        self.R = R
        self.eta = eta
        self.v = v
        self.T_0 = T_0
        self.r_c = r_c
        self.P = []
        self.assignment = {}
        self.clusters = {}
        self.n_clusters = 0
        self.make_particles()
        self.sampled_t = []

    def make_particles(self):
        for i in range(self.N):
            self.P.append(Particle(L=self.L, v=self.v))

    def diff_mat(self):
        mat = np.zeros((self.N, self.N, 3))
        for i in range(self.N):
            for j in range(self.N):  # Comp-time can be halved here
                mat[i, j, :] = self.P[i].r - self.P[j].r

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
        e_z = np.array([0, 0, 1])
        T_n = 0
        r = self.P[n].r
        if self.T_0 != 0:
            for i in range(self.N):  # Do mod L somewhere here
                if i != n:
                    r_ni_ = (self.P[i].r - r)
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
                        T_n += self.T_0 * np.dot(v_vec, r_ni) / np.linalg.norm(r_ni) ** 2 * np.dot(
                            np.cross(v_vec, r_ni), e_z)

        theta = self.P[n].theta
        xi = np.random.uniform(-self.eta, self.eta) / 2
        theta += T_n + xi

        v = self.P[n].v
        self.P[n].r += v * np.array([np.cos(theta), np.sin(theta), 0])
        self.P[n].r = self.P[n].r % (self.L)

        for i in range(self.N):  # Do mod L somewhere here
            if i != n:
                delta_cm = self.periodic_distance(self.P[n].r - self.P[i].r)
                delta_cm_norm = np.linalg.norm(delta_cm)
                if delta_cm_norm < 2 * self.R:
                    delta_cm_normed = delta_cm / delta_cm_norm
                    overlap = 2 * self.R - delta_cm_norm
                    self.P[n].r += delta_cm_normed * 1 / 2 * overlap
                    self.P[i].r -= delta_cm_normed * 1 / 2 * overlap

                    self.P[n].r_[-1] = self.P[n].r
                    self.P[i].r_[-1] = self.P[i].r
                    if self.assignment.__contains__(i) and self.assignment.__contains__(n):

                    elif self.assignment.__contains__(i):
                        self.assignment[n] = self.assignment[i]

                    elif self.assignment.__contains__(n):
                        self.assignment[i] = self.assignment[n]
                    else:
                        self.assignment[i] = self.n_clusters
                        self.assignment[n] = self.n_clusters
                        self.n_clusters += 1

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
            if t % plot_freq == 0:
                self.gen_state_fig(t)
                print('Fin iter %i' % t)
            t += 1

    def gen_state_fig(self, t):
        dir = 'plots1/t%i' % t
        self.sampled_t.append(t)
        plt.figure()
        ax = plt.gca()
        for p in self.P:
            # plt.scatter(p.r[0],p.r[1], color="none", edgecolor="blue")
            circle = plt.Circle((p.r[0], p.r[1]), 1, color='blue', fill=False)
            ax.add_artist(circle)
            for i in range(len(p.r_) - 1, 1, -2):
                r_f = p.r_[i]
                r_i = p.r_[i - 2]
                if np.linalg.norm(r_f - r_i) > self.L / 2:
                    break
                plt.plot([r_f[0], r_i[0]], [r_f[1], r_i[1]], color='b')

            # plt.plot([p.r[0],p.r_[0],p.r__[0]],[p.r[1],p.r_[1],p.r__[1]],c='b')
        plt.axis([0, self.L, 0, self.L])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.axis('off')
        plt.savefig(dir, dpi=100)
        plt.close()

    def plot_grid_time(self, dir):
        fig = plt.figure()

        while True:
            plt.ion()

            for t in self.sampled_t:
                file = dir + '/t' + str(t) + '.png'
                plt.clf()
                image = img.imread(file)
                plt.imshow(image)
                plt.axis('off')
                plt.title('t=%i' % t, fontsize=15)
                plt.pause(0.005)
                plt.draw()


system = System(N=30, L=100, R=1, eta=.002 * np.pi, v=.05, T_0=1, r_c=5)
system.iterate(T=1000, plot_freq=10)
system.plot_grid_time('plots1')

print()
