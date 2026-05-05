import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from scipy.sparse import lil_matrix, kron
from scipy.sparse.linalg import gmres

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


def K_eta_eta(vertices, j, n, det_jacobian):
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


def K_mu_eta(vertices, j, n, det_jacobian):
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


points = spatial_mesh.points[:, :2]
n_nodes = len(points)

# --- Global matrices ---
M_global = lil_matrix((n_nodes, n_nodes))
M_xx_global = lil_matrix((n_nodes, n_nodes))
M_yy_global = lil_matrix((n_nodes, n_nodes))
M_xy_global = lil_matrix((n_nodes, n_nodes))
M_yx_global = lil_matrix((n_nodes, n_nodes))

B_x_global = lil_matrix((n_nodes, n_nodes))
B_y_global = lil_matrix((n_nodes, n_nodes))

b_global = np.zeros(n_nodes)

# --- Element assembly ---
for e, tri in enumerate(triangles):
    verts = points[tri]  # (3,2)

    area = TriangleArea(verts)
    x_grads, y_grads = BuildBasisDxDy(verts)

    dx, dy = BuildBasisDxDy(verts)

    M_loc = M_matrix(verts)

    Mxx_loc = np.zeros((3, 3))
    Myy_loc = np.zeros((3, 3))
    Mxy_loc = np.zeros((3, 3))
    Myx_loc = np.zeros((3, 3))

    for i in range(3):
        for m in range(3):
            Mxx_loc[i, m] = M_xx(x_grads, i, m, area)
            Myy_loc[i, m] = M_yy(y_grads, i, m, area)
            Mxy_loc[i, m] = M_xy(x_grads, y_grads, i, m, area)
            Myx_loc[i, m] = M_yx(x_grads, y_grads, i, m, area)

    # --- Assemble into global ---
    for i_local, i_global in enumerate(tri):
        for m_local, m_global in enumerate(tri):
            M_global[i_global, m_global] += M_loc[i_local, m_local]
            M_xx_global[i_global, m_global] += Mxx_loc[i_local, m_local]
            M_yy_global[i_global, m_global] += Myy_loc[i_local, m_local]
            M_xy_global[i_global, m_global] += Mxy_loc[i_local, m_local]
            M_yx_global[i_global, m_global] += Myx_loc[i_local, m_local]

    # --- Source vector ---
    # ∫ T_i dτ = area / 3
    for i_local, i_global in enumerate(tri):
        b_global[i_global] += Qs[e] * (area / 3) / (4 * np.pi)


# --- Boundary assembly ---
for edge, tag in zip(lines, line_phys):
    n1, n2 = edge
    p1 = points[n1]
    p2 = points[n2]

    edge_vec = p2 - p1
    length = np.linalg.norm(edge_vec)

    # outward normal (rotate edge vector)
    normal = np.array([edge_vec[1], -edge_vec[0]])
    normal = normal / np.linalg.norm(normal)

    nx, ny = normal

    # local edge mass matrix
    B_edge = (length / 6.0) * np.array([[2, 1], [1, 2]])

    nodes = [n1, n2]

    for i_local, i_global in enumerate(nodes):
        for m_local, m_global in enumerate(nodes):
            B_x_global[i_global, m_global] += B_edge[i_local, m_local] * nx
            B_y_global[i_global, m_global] += B_edge[i_local, m_local] * ny


# --- Convert to CSR for solving later ---
M_global = M_global.tocsr()
M_xx_global = M_xx_global.tocsr()
M_yy_global = M_yy_global.tocsr()
M_xy_global = M_xy_global.tocsr()
M_yx_global = M_yx_global.tocsr()
B_x_global = B_x_global.tocsr()
B_y_global = B_y_global.tocsr()

angle_points = angle_mesh.points[:, :2]

angle_tris = []
for block in angle_mesh.cells:
    if block.type == "triangle":
        angle_tris.append(block.data)

angle_tris = np.vstack(angle_tris)

n_angle = len(angle_points)

from scipy.sparse import lil_matrix

K_mu_mu_g = lil_matrix((n_angle, n_angle))
K_eta_eta_g = lil_matrix((n_angle, n_angle))
K_mu_eta_g = lil_matrix((n_angle, n_angle))
K_g = lil_matrix((n_angle, n_angle))

B_mu_g = lil_matrix((n_angle, n_angle))
B_eta_g = lil_matrix((n_angle, n_angle))

A_vec = np.zeros(n_angle)

for tri in angle_tris:
    verts = angle_points[tri]
    detJ = DeterminantJacobian(verts)

    for j_local, j_global in enumerate(tri):
        for n_local, n_global in enumerate(tri):

            K_mu_mu_g[j_global, n_global] += K_mu_mu(verts, j_local, n_local, detJ)
            K_eta_eta_g[j_global, n_global] += K_eta_eta(verts, j_local, n_local, detJ)
            K_mu_eta_g[j_global, n_global] += K_mu_eta(verts, j_local, n_local, detJ)
            K_g[j_global, n_global] += K(j_local, n_local, verts, detJ)

            B_mu_g[j_global, n_global] += B_mu(j_local, n_local, verts, detJ)
            B_eta_g[j_global, n_global] += B_eta(j_local, n_local, verts, detJ)

    for n_local, n_global in enumerate(tri):
        A_vec[n_global] += A(n_local, verts, detJ)


K_mu_mu_g = K_mu_mu_g.tocsr()
K_eta_eta_g = K_eta_eta_g.tocsr()
K_mu_eta_g = K_mu_eta_g.tocsr()
K_g = K_g.tocsr()
B_mu_g = B_mu_g.tocsr()
B_eta_g = B_eta_g.tocsr()

A_global = (
    kron(B_mu_g, B_x_global)
    + kron(B_eta_g, B_y_global)
    + (1 / Sigma_t.mean())
    * (
        kron(K_mu_mu_g, M_xx_global)
        + kron(K_eta_eta_g, M_yy_global)
        + kron(K_mu_eta_g, M_xy_global + M_yx_global)
    )
    + kron(K_g, M_global)
)

b_space = np.zeros(len(points))
b_full = np.kron(A_vec, b_space)

psi, info = gmres(A_global, b_full, atol=1e-10)

n_space = len(points)

psi = psi.reshape((n_space, n_angle))


phi = psi @ A_vec

# --- spatial coordinates ---
points = spatial_mesh.points[:, :2]
x = points[:, 0]
y = points[:, 1]

# --- scalar flux (already computed) ---
# phi shape: (n_nodes,)
z = phi

# --- build triangulation ---
triang = mtri.Triangulation(x, y, triangles=triangles)

# --- plot ---
plt.figure(figsize=(6, 5))
tpc = plt.tricontourf(triang, z, levels=50, cmap="viridis")
plt.colorbar(tpc, label="Scalar Flux φ")

plt.triplot(triang, linewidth=0.3, color="k", alpha=0.4)

plt.title("Scalar Flux on Spatial Mesh Nodes")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")

plt.show()
