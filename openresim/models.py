# %% 1. Import Statements:
# from tabnanny import verbose
import time
from openresim import base, grids, fluids, wells, plots
import numpy as np
import sympy as sym
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import matplotlib.pyplot as plt
from tqdm import tqdm
# from functools import lru_cache
# import threading
# from concurrent import futures


# %% 2. Model Class:
class Model(base.Base):
    '''
    Model class to create a model.
    '''
    name = 'Model'

    def __init__(self, grid: grids.Grid, fluid: fluids.SinglePhaseFluid, well: wells.Well = None, pi: int = None, dt: int = 1, dtype: str = 'double', unit='field'):
        '''
        dd
        '''
        # super().__init__(unit)
        self.grid = grid
        self.fluid = fluid
        self.dtype = dtype
        self.dt = dt
        self.flow_equations_terms_dict = {}
        
        # Initial Step: 
        self.A = None
        self.nsteps = 1
        self.tstep = 0

        # Set properties: compressibility and RHS
        self.set_properties()

        # Pressure array: shape = (nsteps, nx+2)
        self.pressures = np.ones((1, self.grid.nx+2), dtype=self.dtype) * np.nan

        # Initial aressure: required in case of compressible fluid
        self.pi = pi
        if pi is not None:
            self.pressures[0][1:-1] = pi

        # Rate array:
        self.rates = np.zeros((1, self.grid.nx+2), dtype=self.dtype)

        # Well dict:
        self.wells = {}
        if well is not None:
            self.set_well(well)

    def set_transmissibility(self):
        self.transmissibility = self.trans = self.grid.G / \
            (self.fluid.mu * self.fluid.B)
    set_trans = set_transmissibility

    def get_well_G(self, i):
        G_n = 2*np.pi * \
            self.factors['transmissibility conversion'] * \
            self.grid.k[i]*self.grid.dz[i]
        G_d = np.log(self.wells[i]['r_eq']/self.wells[i]['r']*12)
        if 's' in self.wells[i].keys():
            G_d += self.wells[i]['s']
        return G_n/G_d

    def get_well_r_eq(self, i):
        return 0.14*(self.grid.dx[i]**2 + self.grid.dy[i]**2)**0.5

    def set_well(self, i=None, well=None, q=None, pwf=None, r=None, s=None):
        if well is not None:
            if i is None:
                i = well.i
            self.wells[i] = vars(well)
        else:
            assert i is not None, "i must be defined"
            if not i in self.wells.keys():
                self.wells[i] = {}
            if q is not None:
                self.wells[i]['q'] = q
            if pwf is not None:
                self.wells[i]['pwf'] = pwf
            if r is not None:
                self.wells[i]['r'] = r
            if s is not None:
                self.wells[i]['s'] = s
        if 'q' in self.wells[i]:
            self.rates[0][i] = self.wells[i]['q']
        self.wells[i]['r_eq'] = self.get_well_r_eq(i)
        self.wells[i]['G'] = self.get_well_G(i)

    def set_boundary(self, i: int, cond: str, v: float):
        """Set model boundary condition

        Args:
            i (int): boundary grid index
                - can be obtained from Grid class using get_boundaries() method. 
            cond (str): boundary condition
                - contant pressure: 'pressure', 'press', 'p'
                - constant rate: 'rate', 'q'
                - constant pressure gradient: 'pressure gradient', 'press grad', 'gradient', 'grad', 'pg', 'g'
            v (float): value
        """
        if cond in ['rate', 'q']:
            self.rates[0][i] = v
        if cond in ['pressure gradient', 'press grad', 'gradient', 'grad', 'pg', 'g']:
            self.rates[0][i] = self.trans[i] * self.grid.dx[i] * v
        if cond in ['pressure', 'press', 'p']:
            self.pressures[0][i] = v

    def set_boundaries(self, b_dict: dict):
        """

        Args:
            b_dict (dict): _description_
        """
        self.b_dict = b_dict
        for i in b_dict.keys():
            [(cond, v)] = b_dict[i].items()
            self.set_boundary(i, cond, v)

    def __set_RHS(self):
        if self.compressibility_type == 'incompressible':
            self.RHS = np.zeros(self.grid.nx + 2, dtype=self.dtype)
        elif self.compressibility_type == 'compressible':
            self.RHS = (self.grid.volume * self.grid.phi * self.compressibility) / \
                (self.factors['volume conversion'] * self.fluid.B * self.dt)

    def set_properties(self):
        """

        """
        self.set_transmissibility()
        if self.fluid.compressibility_type == self.grid.compressibility_type == 'incompressible':
            self.set_compressibility(0)
        else:
            self.set_compressibility(
                self.fluid.compressibility + self.grid.compressibility)
        self.__set_RHS()
    set_props = set_properties

    def get_i_flow_equation(self, i, verbose=False):
        
        assert i > 0 and i <= self.grid.nx, 'grid index is out of range.'
        exec(f"p{i}=sym.Symbol('p{i}')")

        if i not in self.flow_equations_terms_dict:

            i_neighbors = self.grid.get_neighbors(i)
            i_boundaries = self.grid.get_boundaries(i)

            if verbose:
                print(f'i: {i} - Neighbors: {i_neighbors} - Boundaries: {i_boundaries}')

            # exec(f"p{i}=sym.Symbol('p{i}')")
            # ToDo: keep pressure constant at specific cell (requires A adjust)
            # if not np.isnan(self.pressures[self.tstep][i]):
            #     exec(f"p{i} = {self.pressures[self.tstep][i]}")
            terms = []

            # 1. Flow from grid neighbors:
            for neighbor in i_neighbors:
                exec(f"p{neighbor} = sym.Symbol('p{neighbor}')")
                # To Do: keep pressure constant at specific cell (requires A adjust)
                # if not np.isnan(self.pressures[self.tstep][neighbor]):
                #     exec(f"p{neighbor} = {self.pressures[self.tstep][neighbor]}")
                exec(f"""n_term = self.trans[{min(neighbor,i)}] * ((p{neighbor} - p{i})
                                - (self.fluid.g * (self.grid.z[{neighbor}] - self.grid.z[{i}])))
                """)
                terms.append(locals()['n_term'])

            # 2. Flow from grid boundaries:
            for boundary in i_boundaries:
                exec(f"p{boundary}=sym.Symbol('p{boundary}')")
                if not np.isnan(self.pressures[self.tstep][boundary]):
                    exec(f"""b_term = self.trans[{min(boundary,i)}] * 2 * ((p{boundary} - p{i}) 
                                    - (self.fluid.g * (self.grid.z[{boundary}]-self.grid.z[{i}])))
                    """)
                    exec(f"b_term = b_term.subs(p{boundary}, {self.pressures[self.tstep][boundary]})")
                else:
                    exec(f"b_term = self.rates[self.tstep][{boundary}]")
                terms.append(locals()['b_term'])

            # 3. Flow from grid well:
            if i in self.wells:
                if 'q' in self.wells[i]:
                    terms.append(self.wells[i]['q'])
                else:
                    exec(f"""w_term = - self.wells[{i}]['G'] / (self.fluid.B*self.fluid.mu) * (p{i} - self.wells[{i}]['pwf'])""")
                    terms.append(locals()['w_term'])
                    
            if verbose:
                print('terms:', terms)

            self.flow_equations_terms_dict[i] = terms

        else:
            terms = self.flow_equations_terms_dict[i]

        # 4. Accumulation term:
        if self.RHS[i] == 0:
            a_term = 0
        else:
            try:
                exec(f"accumulation = self.RHS[{i}] * (p{i} - {self.pressures[self.tstep][i]})")
            except:
                raise Exception("Initial pressure (pi) must be specified")
            a_term = locals()['accumulation']

        # 4. Overall grid flow equation:
        i_flow_equation = sym.Eq(sum(terms), a_term)
        if i_flow_equation.lhs.as_coefficients_dict()[1] != 0 or \
            i_flow_equation.rhs.as_coefficients_dict()[1] != 0:
            i_flow_equation = i_flow_equation.simplify()

        # 5. Find lhs and rhs:
        i_lhs = dict(sorted(
            i_flow_equation.lhs.as_coefficients_dict().items(), key=lambda x: str(x[0])))
        i_rhs = i_flow_equation.rhs
        return i_lhs, i_rhs

    # @lru_cache(maxsize=None)
    def get_flow_equations(self, verbose=False):
        for i in self.grid.order[self.grid.i_blocks.astype('bool')]:
            i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
            print(f'Grid {i}: {i_lhs}, {i_rhs}')

    # @lru_cache(maxsize=None)
    def update_matrix(self, i, sparse=True, verbose=False):
        '''
        Update coefficient matrix (A) and result vector (d). 
        Note: arrays are passed by reference.
        '''
        i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
        i_lhs = list(i_lhs.values())
        self.d[i-1] = np.array(i_rhs).astype(self.dtype)
        
        if self.tstep == 0:
            if i == 1: 
                self.A[i-1, i-1:i+1] = i_lhs
            elif i == self.grid.nx:
                self.A[i-1, i-2:i] = i_lhs
            else:
                self.A[i-1, i-2:i+1] = i_lhs


    def __get_matrix(self, sparse=False, verbose=False): 
        '''Create coefficient matrix (A) and result vector (d).
        
        Very slow and used only for testing purposes.
        '''
        if sparse: 
            self.d = ss.lil_matrix((self.grid.nx, 1), dtype=self.dtype)
            self.A = ss.lil_matrix((self.grid.nx, self.grid.nx), dtype=self.dtype)
        else: 
            self.d = np.zeros((self.grid.nx, 1), dtype=self.dtype)
            self.A = np.zeros((self.grid.nx, self.grid.nx), dtype=self.dtype)
        
        for i in range(1, self.grid.nx+1):
            i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
            self.d[i-1] = np.array(i_rhs).astype(self.dtype)
            i_lhs = np.array(list(i_lhs.values())).astype(self.dtype)
            if i == 1: 
                self.A[i-1, i-1:i+1] = i_lhs
            elif i == self.grid.nx:
                self.A[i-1, i-2:i] = i_lhs
            else:
                self.A[i-1, i-2:i+1] = i_lhs
                
        if verbose:      
            print('- A:\n', self.A)
            print('- d:\n', self.d)
            
        
    # @lru_cache(maxsize=None)
    def get_matrix(self, sparse=True, verbose=False):
        '''Create coefficient matrix (A) and result vector (d).
        '''
        if self.grid.D == 1:

            # Construct d vector:
            if all(self.RHS == 0):
                self.d = ss.lil_matrix((self.grid.nx, 1), dtype=self.dtype)
            else:
                try:
                    self.d = ss.lil_matrix((-self.RHS[1:-1]*self.pressures[self.tstep][1:-1]).reshape(-1, 1))
                except:
                    raise Exception("Initial pressure (pi) must be specified")
            if not sparse:
                self.d = self.d.toarray() # ss.lil_matrix(self.d, dtype=self.dtype)

            # Construct A matrix:
            if self.tstep == 0:
                self.A = ss.diags([-self.trans[1:]-self.trans[:-1]-self.RHS[1:-1], 
                                    self.trans[1:-1], # East trans for interior blocks
                                    self.trans[1:-1]], # West trans for interior blocks
                                    [0, 1, -1],
                                    shape=(self.grid.nx, self.grid.nx),
                                    format='lil',  # “dia”, “csr”, “csc”, “lil”
                                    dtype=self.dtype)
                self.A[ 0, 0] = - self.trans[ 1] - self.RHS[ 1]
                self.A[-1,-1] = - self.trans[-2] - self.RHS[-1]
                if not sparse:
                    self.A = self.A.toarray()

            # Update matrix if there is pressure or flow in 'west' boundary:
            if not np.isnan(self.pressures[self.tstep][0]) or self.rates[self.tstep][0] != 0:
                self.update_matrix(1, sparse, verbose)

            # Update matrix if there is pressure or flow in 'east' boundary:
            if not np.isnan(self.pressures[self.tstep][-1]) or self.rates[self.tstep][-1] != 0:
                # at last grid: self.grid.nx or -2
                self.update_matrix(self.grid.nx, sparse, verbose)

            # Update matrix in wells i_blocks:
            for i in self.wells.keys():
                self.update_matrix(i, sparse, verbose)

            if verbose:
                if sparse:
                    print('- A:\n', self.A.toarray())
                    print('- d:\n', self.d.toarray())
                else:
                    print('- A:\n', self.A)
                    print('- d:\n', self.d)

            return self.A, self.d

    # @lru_cache(maxsize=None)
    def solve(self, sparse=True, check_MB=True, update=True, verbose=False):

        self.get_matrix(sparse, verbose)
        # self.__get_matrix(sparse, verbose)
        
        if sparse:
            pressures = ssl.spsolve(self.A.tocsc(), self.d)
        else:
            pressures = np.linalg.solve(self.A, self.d).flatten()
            # same as: np.dot(np.linalg.inv(self.A), self.d)

        if self.grid.D == 1:
            # Update pressures:
            if update:
                self.tstep += 1
                self.pressures = np.vstack(
                    [self.pressures, self.pressures[-1]])
                # self.pressures[self.tstep] = self.pressures[self.tstep-1].copy()
                self.pressures[self.tstep][1:-1] = pressures
                self.rates = np.vstack([self.rates, self.rates[-1]])
                # self.rates[self.tstep] = self.rates[self.tstep-1].copy()

                # Update wells:
                for i in self.wells.keys():
                    if 'q' in self.wells[i]:
                        self.wells[i]['pwf'] = self.pressures[self.tstep][i] + \
                            (self.wells[i]['q']*self.fluid.B *
                             self.fluid.mu/self.wells[i]['G'])
                    if 'pwf' in self.wells[i]:
                        self.wells[i]['q'] = - self.wells[i]['G'] / (self.fluid.B*self.fluid.mu) * \
                            (self.pressures[self.tstep]
                             [i] - self.wells[i]['pwf'])
                        self.rates[self.tstep][i] = self.wells[i]['q']

                # Update boundaries:
                for boundary in self.grid.boundaries:
                    i = 1 if boundary == 0 else boundary-1
                    if not np.isnan(self.pressures[self.tstep][boundary]):
                        self.rates[self.tstep][boundary] = self.trans[min(i, boundary)] * 2 * (
                            (self.pressures[self.tstep][boundary] - self.pressures[self.tstep][i]) - (self.fluid.g * (self.grid.z[boundary]-self.grid.z[i])))

                if check_MB:
                    self.check_MB(verbose)

        if verbose:
            print('- Pressures:\n', self.pressures[self.tstep])
            print('- rates:\n', self.rates[self.tstep])

        return pressures

    def run(self, nsteps=10, sparse=True, check_MB=True, verbose=False):
        self.nsteps += nsteps
        start_time = time.time()
        for _ in tqdm(range(1, nsteps+1), unit='steps', colour='green', position=0, leave=True):
            self.solve(sparse, check_MB, update=True, verbose=verbose)
        duration = round(time.time() - start_time, 2)
        print(f'Simulation run of {nsteps} steps is finished in {duration} seconds.')

    def check_MB(self, verbose=False, error_threshold=0.1):
        """Material Balance Check

        """
        if verbose:
            print(f'Error in step {self.tstep}')
        if self.comp_type == 'incompressible':
            self.error = self.rates[self.tstep].sum()  # must add up to 0
            if verbose:
                print(f"    - Error: {self.error}")
        elif self.comp_type == 'compressible':
            # Check MB error over a time step:
            self.incremental_error = (self.RHS[1:-1] * (self.pressures[self.tstep][1:-1] -
                                      self.pressures[self.tstep-1][1:-1])).sum() / self.rates[self.tstep].sum()
            # Check MB error from initial state to current time step: (less accurate)
            self.cumulative_error = (self.RHS[1:-1] * self.dt * (self.pressures[self.tstep]
                                     [1:-1] - self.pressures[0][1:-1])).sum() / (self.dt * self.tstep * self.rates.sum())
            self.error = abs(self.incremental_error - 1)
            if verbose:
                print(f"    - Incremental Error: {self.incremental_error}")
                print(f"    -  Cumulative Error: {self.cumulative_error}")
                print(f"    -       Total Error: {self.incremental_error+self.cumulative_error}")       
            
        assert abs(self.error) < error_threshold, f"""
        Material balance error ({self.error}) higher than the allowed error ({error_threshold})."""

    def plot(self, property:str='pressures', i:int=None, tstep:int=None):      
        
        if tstep is None: 
            tstep = self.tstep
        
        if i is not None:
            exec(f"plt.plot(self.{property}[:, i].flatten())")
            plt.xlabel('Days')
        elif tstep is not None:
            exec(f"plt.plot(self.{property}[tstep, :].flatten())")
            plt.xlabel('Grid (i)')
            plt.xticks(ticks=range(0, self.grid.nx+2))
        plt.grid()
        plt.show()

    def plot_grid(self, property:str='pressures', tstep:int=None):
        if tstep is None: 
            tstep = self.tstep
        exec(f"plt.imshow(self.{property}[tstep][1:-1][np.newaxis, :])")
        plt.colorbar(label=f'{property.capitalize()} ({self.units[property[:-1]]})')
        plt.title(f'{property.capitalize()} Distribution')
        plt.yticks([])
        plt.xlabel('Grid (i)')
        plt.xticks(ticks=range(0, 4), labels=range(1, 5))
        plt.show()

    def show_grid(self, property:str, show_centers=True, show_boundary=False, show_bounds=False):
        plots.show_grid(self, property, show_centers,
                        show_boundary, show_bounds)

    def copy(self):
        """Copy model (under development)

        Returns:
            _type_: _description_
        """
    # https://stackoverflow.com/questions/48338847/how-to-copy-a-python-class-instance-if-deepcopy-does-not-work
        copy_model = Model(grid=self.grid, fluid=self.fluid, pi=self.pi,
                           dt=self.dt, dtype=self.dtype, unit=self.unit)
        # for w in self.wells:
        #     well = wells.Well(self.wells[w])
        #     copy_model.set_well(well)
        # copy_model.set_boundaries(self.b_dict)
        return copy_model


if __name__ == '__main__':
    grid = grids.CartGrid(nx=4, ny=1, nz=1, dx=300, dy=350,
                        dz=40, phi=0.27, kx=270, comp=1*10**-6, dtype='single')
    fluid = fluids.SinglePhaseFluid(
        mu=0.5, B=1, rho=50, comp=1*10**-5, dtype='single')
    model = Model(grid, fluid, pi=4000, dtype='single')
    model.set_well(i=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
    model.run(nsteps=6, sparse=False, check_MB=True, verbose=False)
    print(model.pressures)
    print(model.rates)