from __future__ import print_function, division

import abc
from collections import namedtuple

import numpy as np


class Element:
    __metaclass__ = abc.ABCMeta
    gauss_points = None
    gauss_weights = None
    dofs = None
    strains_components = None
    
    def __init__(self, nodes):
        self.xe = np.zeros((len(nodes), 3))
        for i in range(3):
            self.xe[:, i] = [n.coordinates[i] for n in nodes]
        self.node_labels = [n.label for n in nodes]

    def J(self, xi, eta, zeta):  # noqa
        return np.dot(self.d(xi, eta, zeta), self.xe)

    def volume(self):
        return np.sum([np.linalg.det(self.J(*gp)*w) for gp, w in zip(self.gauss_points, self.gauss_weights)])

    @abc.abstractmethod
    def d(self, xi, eta, zeta):
        pass

    @abc.abstractmethod
    def B(self, xi, eta, zeta):  # noqa
        pass

    @abc.abstractmethod
    def N(self, xi, eta, zeta):  # noqa
        pass

    def gp_volume(self, i):
        gp = self.gauss_points[i, :]
        return np.linalg.det(self.J(*gp))*self.gauss_weights[i]

    def gauss_point_coordinates(self, i):
        gp = self.gauss_points[i, :]
        N_vec = self.N(*gp)  # noqa
        N = np.zeros((3, 3*len(N_vec)))  # noqa
        for i in range(3):
            N[i, i::3] = N_vec
        return np.dot(N, self.xe.flatten())


class C3D8(Element):
    dofs = 24
    strains_components = 6*8
    local_nodal_pos = np.array([[-1, -1, -1],
                                [1, -1, -1],
                                [1, 1, -1],
                                [-1, 1, -1],
                                [-1, -1, 1],
                                [1, -1, 1],
                                [1, 1, 1],
                                [-1, 1, 1]])
    gauss_points = np.zeros((8, 3))
    counter = 0
    gauss_weights = np.ones(8)
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                gauss_points[counter, :] = k/np.sqrt(3), j/np.sqrt(3), i/np.sqrt(3)
                counter += 1

    def __init__(self, nodes):
        super(C3D8, self).__init__(nodes)

    def d(self, xi, eta, zeta):
        """
        Returns the a matrix 3x8 matrix with derivatives of the shape functions with respect to x y and z
        :param xi:     xi-coordinate
        :param eta:    eta-coordinate
        :param zeta:   zeta-coordinate
        :return:       3x8 matrix with derivatives
        """

        d_matrix = np.zeros((3, 8))
        for i in range(8):
            d_matrix[0, i] = (1. + eta*self.local_nodal_pos[i, 1]) * \
                             (1. + zeta*self.local_nodal_pos[i, 2])*self.local_nodal_pos[i, 0]/8
            d_matrix[1, i] = (1. + xi*self.local_nodal_pos[i, 0]) * \
                             (1. + zeta*self.local_nodal_pos[i, 2])*self.local_nodal_pos[i, 1]/8
            d_matrix[2, i] = (1. + xi*self.local_nodal_pos[i, 0]) * \
                             (1. + eta*self.local_nodal_pos[i, 1])*self.local_nodal_pos[i, 2]/8
        return d_matrix

    def N(self, xi, eta, zeta):
        return ((1 + self.local_nodal_pos[:, 0]*xi)
                * (1 + self.local_nodal_pos[:, 1]*eta)
                * (1 + self.local_nodal_pos[:, 2]*zeta)/8)

    def B(self, xi, eta, zeta):
        B = np.zeros((6, 24))  # noqa
        jacobian = self.J(xi, eta, zeta)
        d = self.d(xi, eta, zeta)
        for i in range(8):
            dx_avg = [np.linalg.solve(self.J(*gp), self.d(*gp)[:, i])*np.linalg.det(self.J(*gp))
                      for gp in self.gauss_points]
            dx_avg = sum(dx_avg)/self.volume()

            for j in range(3):
                for k in range(3):
                    B[j, 3*i + k] += dx_avg[k]/3

            dx = np.linalg.solve(jacobian, d[:, i])
            for j in range(3):
                for k in range(3):
                    if j == k:
                        B[j, 3*i + k] += 2*dx[k]/3
                    else:
                        B[j, 3*i + k] += -dx[k]/3

            B[3, 3*i] += dx[1]
            B[3, 3*i + 1] += dx[0]

            B[4, 3*i] += dx[2]
            B[4, 3*i + 2] += dx[0]

            B[5, 3*i + 1] += dx[2]
            B[5, 3*i + 2] += dx[1]
        return B


if __name__ == '__main__':
    Node = namedtuple('Node', ['coordinates', 'label'])
    node_list = [Node(coordinates=[0, 0, 0], label=1),
                 Node(coordinates=[2, 0, 0], label=2),
                 Node(coordinates=[2, 2, 0], label=3),
                 Node(coordinates=[0, 2, 0], label=4),
                 Node(coordinates=[0, 0, 1], label=5),
                 Node(coordinates=[2, 0, 1], label=6),
                 Node(coordinates=[2, 2, 1], label=7),
                 Node(coordinates=[0, 2, 1], label=8)]
    element = C3D8(node_list)
    element.B(0, 0, 0)
    print(element.gauss_point_coordinates(4))


elements = {'C3D8': C3D8}
