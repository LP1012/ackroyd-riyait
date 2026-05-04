import meshio
import numpy as np
import scipy

from scipy.integrate import dblquad

from linear_triangles import *


def mu(theta, phi):
    return np.sin(phi) * np.cos(theta)


def eta(theta, phi):
    return np.sin(phi) * np.sin(theta)


def K_mu_mu(vertices, j, n):
    result = dblquad(
        K_mu_mu_integrand, 0, 1, lambda xi: 0, lambda xi: 1 - xi, args=(j, n, vertices)
    )
    return result


def K_mu_mu_integrand(eta, xi, j, n, vertices):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    mu_sqrd = mu(theta, phi) ** 2
    det_jacobian = DeterminantJacobian(vertices)
    return (
        mu_sqrd
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


spatial_mesh = meshio.read("../gmsh/square-source-coarse.msh")
angle_mesh = meshio.read("../gmsh/solid-angle.msh")
