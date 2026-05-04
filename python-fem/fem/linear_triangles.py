import numpy as np


def CreateJacobian(vertices):
    x1, y1 = vertices[0]
    x2, y2 = vertices[1]
    x3, y3 = vertices[2]

    return np.array([[x2 - x1, x3 - x1], [y2 - y1, y3 - y1]])


def DeterminantJacobian(vertices):
    return np.abs(np.linalg.det(CreateJacobian(vertices)))


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
