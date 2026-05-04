import meshio
import numpy as np
import scipy

from scipy.integrate import dblquad, quad

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
    )[0]
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


def K_eta_eta(xi, eta, j, n, vertices, det_jacobian):
    result = dblquad(
        K_eta_eta_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )[0]
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


def K_mu_eta(xi, eta, j, n, vertices, det_jacobian):
    result = dblquad(
        K_mu_eta_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )[0]
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
    return y_gradients[i] * x_gradients[m] * triangle_area


def B_mu(j, n, vertices, det_jacobian):
    result = dblquad(
        B_mu_integrand,
        0,
        1,
        lambda xi: 0,
        lambda xi: 1 - xi,
        args=(j, n, vertices, det_jacobian),
    )[0]
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
    )[0]
    return result


def B_eta_integrand(eta, xi, j, n, vertices, det_jacobian):
    theta, phi = ReferenceToReal(vertices, xi, eta)
    y_cosine = Eta(theta, phi)
    return (
        y_cosine
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
    )[0]
    return result


def K_integrand(eta, xi, j, n, vertices, det_jacobian):
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
    )[0]
    return result


def A_integrand(eta, xi, n, vertices, det_jacobian):
    _, phi = ReferenceToReal(vertices, xi, eta)
    return LinearBasis(xi, eta, n) * det_jacobian * np.sin(phi)


def S(xi, eta, vertices, det_jacobian, scattering_cross_sections):
    xs = scattering_cross_sections / 4 / np.pi
    return xs * M_matrix(vertices)


def Q(m, vertices, det_jacobian, Qs):
    return (
        Qs[m]
        * dblquad(
            Q_integrand,
            0,
            1,
            lambda xi: 0,
            lambda xi: 1 - xi,
            args=(m, vertices, det_jacobian),
        )[0]
    )


def Q_integrand(eta, xi, m, vertices, det_jacobian):
    return LinearBasis(xi, eta, m) * det_jacobian


# to do:
# - write methods for local matrix creation
# - write methods to evaluate boundary integrals (must always be positive)
# - global matrix assembly
# - material properties (total cross section, scattering cross section, source)
# - scalar flux computation

spatial_mesh = meshio.read("gmsh/square-source.msh")
angle_mesh = meshio.read("gmsh/solid-angle.msh")

triangles = []
triangle_phys = []

lines = []
line_phys = []

for block, phys in zip(spatial_mesh.cells, spatial_mesh.cell_data["gmsh:physical"]):
    if block.type == "triangle":
        triangles.append(block.data)
        triangle_phys.append(phys)
    elif block.type == "line":
        lines.append(block.data)
        line_phys.append(phys)


triangles = np.vstack(triangles)
triangle_phys = np.concatenate(triangle_phys)

lines = np.vstack(lines)
line_phys = np.concatenate(line_phys)

west_edges = lines[line_phys == 13]
north_edges = lines[line_phys == 14]
east_edges = lines[line_phys == 15]
south_edges = lines[line_phys == 16]

print("Triangles:", len(triangles))
print("Lines:", len(lines))

print("Unique triangle regions:", np.unique(triangle_phys))
print("Unique boundary tags:", np.unique(line_phys))

Sigma_t = np.zeros(len(triangles))
Sigma_s = np.zeros(len(triangles))
Qs = np.zeros(len(triangles))

# check the numbering on these once you run... I totally guessed :-)
Sigma_t[triangle_phys == 1] = 0.2
Sigma_t[triangle_phys == 2] = 1e-5
Sigma_t[triangle_phys == 3] = 0.2

Sigma_s[triangle_phys == 1] = 0
Sigma_s[triangle_phys == 2] = 0
Sigma_s[triangle_phys == 3] = 0

Qs[triangle_phys == 3] = 6  # fix value later...
