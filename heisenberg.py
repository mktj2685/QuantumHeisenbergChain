import numpy as np
from numpy import linalg as LA


def Hamiltonian(n:int, bc:str='open'):
    """
    Calculate matrix representation, eignenvalues and eigenvectors of below Hamiltonian.

        H = Σ^n_i Hi
        Hi = (S+(i)S-(i+1) + S-(i)S+(i+1))/2 + Sz(i)Sz(i+1)
        
        where
            n : number of spins.
            S(i) : spin operator on i-th site.
    """
    # Set spin operators.
    Sp = np.array([[0, 1], [0, 0]], dtype=np.float64)
    Sm = np.array([[0, 0], [1, 0]], dtype=np.float64)
    Sz = 0.5 * np.array([[1, 0], [0, -1]], dtype=np.float64)

    # Set identity matrix.
    I2 = np.identity(2, dtype=np.float64)

    # Calculate matrix representation of Hamiltonian.
    H = np.zeros((2**n, 2**n), dtype=np.float64)

    ## Bulk part.
    for i in range(n-1):
        Hi = 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp)) + np.kron(Sz, Sz)
        for _ in range(0, i):
            Hi = np.kron(I2, Hi)
        for _ in range(i+2, n):
            Hi = np.kron(Hi, I2)
        H += Hi

    ## Boundary part
    # NOTE Plese comment out follow lines if you want open boudary.
    if bc == 'periodic':
        In2 = np.identity(2**(n-2), dtype=np.float64)
        Hi = 0.5 * (np.kron(Sp, np.kron(In2, Sm)) + np.kron(Sm, np.kron(In2, Sp))) + np.kron(Sz, np.kron(In2, Sz))
        H += Hi

    # Calculate eigenvalues and eigenvectors.
    w, v = LA.eigh(H)

    # Set return object.
    ret = {
        'Hamiltonian': H,
        'Eigenvalues': w,
        'Eigenvectors': v
    }

    return ret


if __name__ == '__main__':
    ret = Hamiltonian(12)
    print(ret['Hamiltonian'])
    print(ret['Eigenvalues'])
    print(ret['Eigenvectors'])