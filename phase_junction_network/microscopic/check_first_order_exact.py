#!/usr/bin/env python3
"""Exact rational coefficient check for the first-order frame/connection identity."""
from __future__ import annotations
import argparse
import itertools
import json
from fractions import Fraction as F


def h_basis():
    out = []
    for component in range(6):
        h = [[F(0) for _ in range(3)] for _ in range(3)]
        if component < 3:
            h[component][component] = F(1)
        else:
            i, j = ((0, 1), (0, 2), (1, 2))[component - 3]
            h[i][j] = h[j][i] = F(1)
        out.append(h)
    return out


def c_basis():
    out = []
    for a in range(3):
        for b in range(3):
            for c in range(b, 3):
                connection = [
                    [[F(0) for _ in range(3)] for _ in range(3)]
                    for _ in range(3)
                ]
                connection[a][b][c] = F(1)
                connection[a][c][b] = F(1)
                out.append(connection)
    return out


HB = h_basis()
CB = c_basis()


def add_tensor(left, right):
    if isinstance(left[0], list):
        return [add_tensor(a, b) for a, b in zip(left, right)]
    return [a + b for a, b in zip(left, right)]


def q_connection(connection):
    value = F(0)
    for i, a, b in itertools.product(range(3), repeat=3):
        value += (
            connection[a][b][i] * connection[b][a][i]
            - connection[a][i][i] * connection[b][a][b]
        )
    return value


def quadratic_matrix(basis, value_function):
    size = len(basis)
    matrix = [[F(0) for _ in range(size)] for _ in range(size)]
    diagonal = [value_function(element) for element in basis]
    for i in range(size):
        matrix[i][i] = diagonal[i]
        for j in range(i):
            matrix[i][j] = matrix[j][i] = (
                value_function(add_tensor(basis[i], basis[j]))
                - diagonal[i]
                - diagonal[j]
            ) / 2
    return matrix


J = quadratic_matrix(CB, q_connection)


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(left, right):
    return [
        [
            sum(
                (left[i][inner] * right[inner][j] for inner in range(len(right))),
                F(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def scale(scalar, matrix):
    return [[scalar * value for value in row] for row in matrix]


def gamma(h, k):
    connection = [
        [[F(0) for _ in range(3)] for _ in range(3)]
        for _ in range(3)
    ]
    for a, b, c in itertools.product(range(3), repeat=3):
        connection[a][b][c] = (
            k[b] * h[a][c] + k[c] * h[a][b] - k[a] * h[b][c]
        ) / 2
    return connection


def connection_coordinates(connection):
    return [
        connection[a][b][c]
        for a in range(3)
        for b in range(3)
        for c in range(b, 3)
    ]


def gamma_map(k):
    columns = [connection_coordinates(gamma(h, k)) for h in HB]
    return [list(row) for row in zip(*columns)]


def fp_value(h, k):
    k_squared = sum((value * value for value in k), F(0))
    h_squared = sum(
        (h[i][j] * h[i][j] for i in range(3) for j in range(3)),
        F(0),
    )
    k_dot_h = [
        sum((k[i] * h[i][j] for i in range(3)), F(0))
        for j in range(3)
    ]
    trace = sum((h[i][i] for i in range(3)), F(0))
    return (
        k_squared * h_squared
        - 2 * sum((value * value for value in k_dot_h), F(0))
        + 2 * sum((k_dot_h[j] * k[j] for j in range(3)), F(0)) * trace
        - k_squared * trace * trace
    )


def fp_matrix(k):
    return quadratic_matrix(HB, lambda h: fp_value(h, k))


def schur_matrix(k):
    mapping = gamma_map(k)
    return scale(F(-4), matmul(transpose(mapping), matmul(J, mapping)))


def h_coordinates(h):
    return [h[0][0], h[1][1], h[2][2], h[0][1], h[0][2], h[1][2]]


def gauge_matrix(k):
    columns = []
    for direction in range(3):
        xi = [F(0), F(0), F(0)]
        xi[direction] = F(1)
        h = [
            [k[i] * xi[j] + xi[i] * k[j] for j in range(3)]
            for i in range(3)
        ]
        columns.append(h_coordinates(h))
    return [list(row) for row in zip(*columns)]


def is_zero(matrix):
    return all(value == 0 for row in matrix for value in row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()

    coefficient_points = [
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
        (F(1), F(1), F(0)),
        (F(1), F(0), F(1)),
        (F(0), F(1), F(1)),
    ]
    coefficient_failures = [
        tuple(int(x) for x in k)
        for k in coefficient_points
        if schur_matrix(k) != fp_matrix(k)
    ]

    gauge_points = [
        (F(2), F(-1), F(3)),
        (F(1), F(2), F(4)),
        (F(-3), F(5), F(2)),
        (F(7), F(-4), F(1)),
    ]
    gauge_failures = [
        tuple(int(x) for x in k)
        for k in gauge_points
        if not is_zero(matmul(fp_matrix(k), gauge_matrix(k)))
    ]

    report = {
        "status": "exact rational coefficient verification",
        "connection_components": len(CB),
        "frame_components": len(HB),
        "homogeneous_quadratic_coefficient_points": [
            [int(x) for x in k] for k in coefficient_points
        ],
        "schur_equals_fierz_pauli_coefficient_failures": coefficient_failures,
        "exact_gauge_null_test_points": [
            [int(x) for x in k] for k in gauge_points
        ],
        "exact_gauge_null_failures": gauge_failures,
        "interpretation": (
            "A homogeneous quadratic 6x6 kernel in three momentum variables is "
            "fixed by the three axis and three pair-sum evaluations. Equality at "
            "all six points proves coefficient equality over the rationals."
        ),
    }
    assert not coefficient_failures
    assert not gauge_failures
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
