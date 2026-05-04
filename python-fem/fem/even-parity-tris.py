import meshio
import numpy as np
import scipy

from scipy.integrate import dblquad

from linear_triangles import *


def Mu(theta, phi):
    return np.sin(phi) * np.cos(theta)


def Eta(theta, phi):
    return np.sin(phi) * np.sin(theta)


def K_mu_mu(vertices, j, n, det_jacobian):
    result = dblquad(
        K_mu_mu_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def K_mu_mu_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    mu_sqrd = Mu(theta, phi) ** 2
    return (
        mu_sqrd
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


def K_eta_eta(eta, xi, j, n, vertices, det_jacobian):
    result = dblquad(
        K_eta_eta_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def K_eta_eta_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    eta_sqrd = Eta(theta, phi) ** 2
    return (
        eta_sqrd
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


def K_mu_eta(eta, xi, j, n, vertices, det_jacobian):
    result = dblquad(
        K_mu_eta_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def K_mu_eta_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    eta_mu = Eta(theta, phi) * Mu(theta, phi)
    return (
        eta_mu
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


def M_xx(x_gradients, i, m, triangle_area):
    return x_gradients[i] * x_gradients[m] * triangle_area


def M_yy(y_gradients, i, m, triangle_area):
    return y_gradients[i] * y_gradients[m] * triangle_area


def M_xy(x_gradients, y_gradients, i, m, triangle_area):
    return x_gradients[i] * y_gradients[m] * triangle_area


def M_yx(x_gradients, y_gradients, i, m, triangle_area):
    return x_gradients[i] * x_gradients[m] * triangle_area


def B_mu(j, n, vertices, det_jacobian):
    result = dblquad(
        B_mu_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def B_mu_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    mu = Mu(theta, phi)
    return (
        mu
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


def B_eta(j, n, vertices, det_jacobian):
    result = dblquad(
        B_eta_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def B_eta_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    eta = Eta(theta, phi)
    return (
        eta
        * LinearBasis(xi, eta, j)
        * LinearBasis(xi, eta, n)
        * np.sin(phi)
        * det_jacobian
    )


def K(j, n, vertices, det_jacobian):
    result = dblquad(
        K_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )
    return result


def K_integrand(xi, eta, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    return (
        LinearBasis(xi, eta, j) * LinearBasis(xi, eta, n) * np.sin(phi) * det_jacobian
    )


def M_matrix(verices):
    area = TriangleArea(verices)
    return (area / 6) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=float)


def A(n, vertices, det_jacobian):
    result = dblquad(
        A_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(n, vertices, det_jacobian),
    )
    return result


def A_integrand(xi, eta, n, vertices, det_jacobian):
    _, phi = ReferenceToReal(vertices, xi, eta)
    return LinearBasis(xi, eta, n) * det_jacobian * np.sin(phi)


spatial_mesh = meshio.read("../gmsh/square-source.msh")
angle_mesh = meshio.read("../gmsh/solid-angle.msh")

print(spatial_mesh)
line_cells = next(cb.data for cb in spatial_mesh.cells if cb.type == "line")
print(line_cells)
