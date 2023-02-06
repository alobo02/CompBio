#!/usr/bin/env python3

# TODO: Module DocString
"""_summary_
"""

import sys
import numpy as np


class GlobalSeqeuenceAlignment:
    """_summary_
    """

    def __init__(
        self,
        x: str,
        y: str,
        M: float = 4,
        m: float = -2,
        g: float = -2
    ) -> None:
        """_summary_

        Parameters
        ----------
        x : str
            _description_
        y : str
            _description_
        M : float, optional
            _description_, by default 4
        m : float, optional
            _description_, by default -2
        g : float, optional
            _description_, by default -2
        """

        # Strings
        self.x = '_' + x
        self.y = '_' + y
        # Scoring Scheme
        self.M = M
        self.m = m
        self.g = g
        # Obtain dimension of matrix
        m_dim, n_dim = len(self.x), len(self.y)
        self.m_dim = m_dim
        self.n_dim = n_dim
        # Initialize alignment (A) matrix
        A = np.zeros((m_dim, n_dim))
        A[:, 0] = np.arange(0, g, g*m_dim)
        A[0, :] = np.arange(0, g, g*n_dim)
        self.A = A
        # Initialize pointer (P) matrix
        P = np.zeros((m_dim, n_dim))
        self.P = P
        # Compute matrix scores
        self._compute_matrix()

    def _compute_matrix(self) -> None:
        """_summary_
        """
        for i in range(1, self.m_dim):
            for j in range(1, self.n_dim):
                self._update_matrix(i, j)

    def _update_matrix(self, i: int, j: int) -> None:
        """_summary_

        Parameters
        ----------
        i : int
            _description_
        j : int
            _description_
        """
        cases = [
            self.A[i-1, j-1] + self._score(self.x[i], self.y[j]),
            self.A[i-1, j] + self.g,
            self.A[i, j-1] + self.g
        ]
        self.A[i, j] = max(cases)
        self.P[i, j] = np.argmax(cases) + 1

    def _score(self, x_i: str, y_j: str) -> float:
        """_summary_

        Parameters
        ----------
        x_i : str
            _description_
        y_j : str
            _description_

        Returns
        -------
        float
            _description_
        """
        if x_i == y_j:
            return self.M
        else:
            return self.m

    def get_optimal_alignment(self) -> tuple[str, str]:
        """_summary_

        Returns
        -------
        tuple[str, str]
            _description_
        """
        i = self.m_dim-1
        j = self.n_dim-1

        x_stack = []
        y_stack = []

        while i > 0 and j > 0:
            pointer = self.P[i, j]
            if pointer == 1:
                x_stack.append(self.x[i])
                y_stack.append(self.y[j])
                i -= 1
                j -= 1
            elif pointer == 2:
                x_stack.append(self.x[i])
                y_stack.append('-')
                i -= 1
            elif pointer == 3:
                x_stack.append('-')
                y_stack.append(self.y[j])
                j -= 1
        x_stack.extend(['-'] * i)
        y_stack.extend(['-'] * j)

        x_aligned = ''.join(x_stack[::-1])
        y_aligned = ''.join(y_stack[::-1])

        return (x_aligned, y_aligned)


def parse_FASTA_file(fasta_fh) -> tuple[str, str]:
    """_summary_

    Parameters
    ----------
    fasta_fh : _type_
        _description_

    Returns
    -------
    tuple[str, str]
        _description_
    """
    x_y = [[], []]
    i = -1
    for line in fasta_fh.readlines():
        if line.startswith('>'):
            i += 1
        else:
            x_y[i].append(line.rstrip())
    x_y = tuple(map(lambda z: ''.join(z), x_y))
    return x_y


if __name__ == "__main__":
    x, y = parse_FASTA_file(sys.stdin)
    align = GlobalSeqeuenceAlignment(x, y)
    print(*align.get_optimal_alignment(), sep='\n')
