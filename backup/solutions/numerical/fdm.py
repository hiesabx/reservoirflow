"""Finite Difference Method (FDM) Module"""
import numpy as np
import scipy.sparse as ss
import sympy as sym

from reservoirflow.utils.helpers import _lru_cache


class FDM:
    """Finite Difference Method (FDM) Class"""

    def __init__(
        self,
        model,
    ):
        self.model = model
        self.dtype = model.dtype
        self.grid = model.grid
        self.fluid = model.fluid

    @_lru_cache(maxsize=None)
    def get_cell_trans(
        self,
        cell_id=None,
        cell_coord=None,
        boundary: bool = False,
    ):
        """Returns transmissibility (T) at all cell faces.

        Parameters
        ----------
        cell_id : int, iterable of int
            cell id based on natural order as int.
        cell_coord : iterable of int
            cell coordinate (i,j,k) as a tuple of int.
        boundary : bool, optional
            include boundary cells.

        Returns
        -------
        dict
            dictionary with cell neighbors' id as keys and
            transmissibility as values.
        """
        cell_G = self.grid.get_cell_G(cell_id, cell_coord, boundary)
        mu_B = self.fluid.mu * self.fluid.B
        return {k: v / mu_B for k, v in cell_G.items()}

    def get_cells_trans(
        self,
        boundary: bool = False,
        sparse=False,
        vectorize=True,
    ):
        if vectorize:
            mu_B = self.fluid.mu * self.fluid.B
            return self.grid.get_cells_G(boundary, sparse) / mu_B

        n = self.grid.get_n(boundary)
        if sparse:
            cells_T = ss.lil_matrix((n, n), dtype=self.dtype)
        else:
            cells_T = np.zeros((n, n), dtype=self.dtype)

        if boundary:
            for cell_id in self.grid.cells_id:
                T = self.get_cell_trans(cell_id, None, False)
                for cell_n_id in T.keys():
                    cells_T[cell_id, cell_n_id] = T[cell_n_id]
        else:
            for cell_id in self.grid.cells_id:
                cell_i = self.cells_i_dict[cell_id]
                T = self.get_cell_trans(cell_id, None, False)
                cells_n_i = [self.cells_i_dict[x] for x in T.keys()]
                for cell_n_id, cell_n_i in zip(T.keys(), cells_n_i):
                    cells_T[cell_i, cell_n_i] = T[cell_n_id]

        return cells_T
