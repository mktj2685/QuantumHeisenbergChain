from typing import Dict
import argparse
from math import sqrt
import numpy as np
from numpy import linalg as LA


def parse_args():
    parser = argparse.ArgumentParser(
        description='Exact diagonalization of quantum Heisenberg chain.')
    parser.add_argument('--S', type=float, required=True,
                        help='Consider spin-S chain.')
    parser.add_argument('--N', type=int, required=True,
                        help='Number of spins.')

    return parser.parse_args()

def operators(S: float) -> Dict[str, np.ndarray]:
    # operator matrix dimension is 2S+1.
    dim = int(2*S+1)

    # Kronecker delta function.
    def delta(x, y): return 1 if x == y else 0

    # Initialize dictionary of spin operators.
    ops = {
        'Sx': np.zeros((dim, dim), dtype=np.float64),
        'Sy': np.zeros((dim, dim), dtype=np.complex128),
        'Sz': np.zeros((dim, dim), dtype=np.float64),
        'S+': np.zeros((dim, dim), dtype=np.float64),
        'S-': np.zeros((dim, dim), dtype=np.float64),
        'S2': np.zeros((dim, dim), dtype=np.float64),
    }

    # Calculate elements of each operator.
    for i in range(dim):
        for j in range(dim):
            m1 = S - i
            m2 = S - j
            ops['Sx'][i, j] = (delta(m1, m2+1)+delta(m1+1, m2)
                               )*(sqrt(S*(S+1)-m1*m2))/2
            ops['Sy'][i, j] = -1j * \
                (delta(m1, m2+1)-delta(m1+1, m2))*(sqrt(S*(S+1)-m1*m2))/2
            ops['Sz'][i, j] = delta(m1, m2)*m1
            ops['S+'][i, j] = delta(m1, m2+1)*(sqrt(S*(S+1)-m1*m2))
            ops['S-'][i, j] = delta(m1+1, m2)*(sqrt(S*(S+1)-m1*m2))
            ops['S2'][i, j] = delta(m1, m2)*(S*(S+1))

    return ops


def hamiltonian(S: float, N: int, bc: str = 'periodic') -> np.ndarray:
    """
    Calculate matrix representation, eignenvalues and eigenvectors of below Hamiltonian.
        H = Σ^n_i Hi
        Hi = (S+(i)S-(i+1) + S-(i)S+(i+1))/2 + Sz(i)Sz(i+1)

        where
            n : number of spins.
            S(i) : spin operator on i-th site.
    """
    dim1 = int(2*S+1)   # 1-spin state space
    dimN = dim1**N      # N-spin state

    H = np.zeros((dimN, dimN), dtype=np.float64)    # Hamiltonian
    # Identity matrix witch operates single spin.
    I = np.identity(dim1, dtype=np.float64)

    ops = operators(S)  # spin operators.

    # Calculate Hamiltonian elements (bulk part).
    for i in range(N-1):
        Hi = 0.5 * (np.kron(ops['S+'], ops['S-']) + np.kron(ops['S-'],
                    ops['S+'])) + np.kron(ops['Sz'], ops['Sz'])
        for _ in range(0, i):
            Hi = np.kron(I, Hi)
        for _ in range(i+2, N):
            Hi = np.kron(Hi, I)
        H += Hi

    # Calculate Hamiltonian elements (boundary part).
    if bc == 'periodic':
        # Identity matrix which operates bulk spins.
        I = np.identity(dim1**(N-2), dtype=np.float64)
        Hi = 0.5 * (np.kron(ops['S+'], np.kron(I, ops['S-'])) + np.kron(
            ops['S-'], np.kron(I, ops['S+']))) + np.kron(ops['Sz'], np.kron(I, ops['Sz']))
        H += Hi

    return H


if __name__ == '__main__':
    # Get arguments.
    args = parse_args()

    # Get matrix representation.
    H = hamiltonian(args.S, args.N)

    # Diagonalize.
    w, v = LA.eigh(H)

    # print eigenenergy and eigenstate.
    for i in range(4):
        print(f'E{i} = {w[i]}')
        print(f'|ψ{i}> = {v[:,i]}')
