import numpy as np


def CreateJacobian(vertices):
    x1, y1 = vertices[0]
    x2, y2 = vertices[1]
    x3, y3 = vertices[2]

    return np.array([[x2 - x1, x3 - x1], [y2 - y1, y3 - y1]])


def DeterminantJacobian(vertices):
    return np.abs(np.linalg.det(CreateJacobian(vertices)))


def TriangleArea(vertices):
    return DeterminantJacobian(vertices) / 2


def ReferenceToReal(vertices, xi, eta):
    x1, y1 = vertices[0]
    x2, y2 = vertices[1]
    x3, y3 = vertices[2]

    x = x1 + (x2 - x1) * xi + (x3 - x1) * eta
    y = y1 + (y2 - y1) * xi + (y3 - y1) * eta

    return x, y


def LinearBasis(xi, eta, i):
    match i:
        case 0:
            return 1 - xi - eta
        case 1:
            return xi
        case 2:
            return eta
        case _:
            raise ValueError(f"i must be 0, 1, or 2. i={i}")


def LinearBasisXiDerivative(i):
    match i:
        case 0:
            return -1
        case 1:
            return 1
        case 2:
            return 0
        case _:
            raise ValueError(f"i must be 0, 1, or 2. i={i}")


def LinearBasisEtaDerivative(i):
    match i:
        case 0:
            return -1
        case 1:
            return 0
        case 2:
            return 1
        case _:
            raise ValueError(f"i must be 0, 1, or 2. i={i}")


def RealGradients(vertices, i):
    jacobian = CreateJacobian(vertices)
    xi_derivative = LinearBasisXiDerivative(i)
    eta_derivative = LinearBasisEtaDerivative(i)
    reference_gradients = np.array([xi_derivative, eta_derivative])
    real_gradients = np.linalg.solve(jacobian.T, reference_gradients)
    return real_gradients


def BuildBasisDxDy(vertices):
    dxs = np.zeros(3)
    dys = np.zeros(3)
    for i in range(3):
        gradients = RealGradients(vertices, i)
        dxs[i] = gradients[0]  # dT_i/dx
        dys[i] = gradients[1]  # dT_i/dy
    return dxs, dys
