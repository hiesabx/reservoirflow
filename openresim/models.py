#%% 1. Import Statements:
from tabnanny import verbose
from openresim import base, grids, fluids, wells, plots
import numpy as np
import sympy as sym
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import matplotlib.pyplot as plt
# import copy


#%% 2. Model Class: 
class Model(base.Base):
    '''
    Model class to create a model.
    '''
    name = 'Model'
    def __init__(self, grid, fluid, well=None, pi=None, dt=1, dtype='double', unit='field'):
        '''
        dd
        '''
        #super().__init__(unit)
        self.dt = dt
        self.dtype = dtype # np.single, np.double
        self.grid = grid # composition
        self.fluid = fluid # composition
        self.pressures = (self.grid.blocks * np.nan).astype(self.dtype)
        self.pi = pi
        if pi != None:
            self.pressures[1:-1] = pi
        self.rates = np.zeros(grid.nx + 2, dtype=self.dtype)
        self.wells = {}
        if well != None:
            self.set_well(well)
        self.set_properties()


    def set_transmissibility(self):
        self.transmissibility = self.trans = self.grid.G / (self.fluid.mu * self.fluid.B)
    set_trans = set_transmissibility


    def get_well_G(self, i):
        G_n =  2*np.pi*self.factors['transmissibility conversion']*self.grid.k[i]*self.grid.dz[i]
        G_d = np.log(self.wells[i]['r_eq']/self.wells[i]['r']*12)
        if 's' in self.wells[i].keys():
            G_d += self.wells[i]['s']
        return G_n/G_d


    def get_well_r_eq(self, i):
        return 0.14*(self.grid.dx[i]**2 + self.grid.dy[i]**2)**0.5


    def set_well(self, i=None, well=None, q=None, pwf=None, r=None, s=None):
        if well != None:
            if i == None:
                i = well.i
            self.wells[i] = vars(well)
        else:
            assert i != None, "i must be defined"
            if not i in self.wells.keys():
                self.wells[i] = {}
            if q != None:
                self.wells[i]['q'] = q
            if pwf != None:
                self.wells[i]['pwf'] = pwf
            if r != None:
                self.wells[i]['r'] = r
            if s != None:
                self.wells[i]['s'] = s
        if 'q' in self.wells[i]:
            self.rates[i] = self.wells[i]['q']
        self.wells[i]['r_eq'] = self.get_well_r_eq(i)
        self.wells[i]['G'] = self.get_well_G(i)


    def set_boundary(self, i, cond, v):
        if cond in ['rate', 'q']:
            self.rates[i] = v
        if cond in ['pressure gradient', 'press grad', 'gradient', 'grad', 'pg', 'g']:
            self.rates[i] = self.trans[i] * self.grid.dx[i] * v
        if cond in ['pressure', 'press', 'p']:
            self.pressures[i] = v


    def set_boundaries(self, b_dict):
        '''
        
        '''
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
            self.set_compressibility(self.fluid.compressibility + self.grid.compressibility)
        self.__set_RHS()
    set_props = set_properties


    def get_i_flow_equation(self, i, verbose=False):
        # ToDo: look at trans!
        # if i < 0:
        #     i = self.grid.nx + 2 + i
        assert i > 0 and i <= self.grid.nx, 'grid index is out of range.'

        i_neighbors = self.grid.get_neighbors(i)
        i_boundaries = self.grid.get_boundaries(i)

        if verbose: print(f'i: {i} - Neighbors: {i_neighbors} - Boundaries: {i_boundaries}')
        
        exec(f"p{i}=sym.Symbol('p{i}')")
        # ToDo: keep pressure constant at specific cell (requires A adjust)
        # if not np.isnan(self.pressures[i]):
        #     exec(f"p{i} = {self.pressures[i]}")
        terms = []

        # 1. Flow from grid neighbors:
        for neighbor in i_neighbors:
            exec(f"p{neighbor} = sym.Symbol('p{neighbor}')")
            # To Do: keep pressure constant at specific cell (requires A adjust)
            # if not np.isnan(self.pressures[neighbor]):
            #     exec(f"p{neighbor} = {self.pressures[neighbor]}")
            exec(f"n_term = self.trans[{neighbor}] * ((p{neighbor} - p{i}) - (self.fluid.g * (self.grid.z[{neighbor}]-self.grid.z[{i}])))")
            terms.append(locals()['n_term'])

        # 2. Flow from grid boundaries:
        for boundary in i_boundaries:
            exec(f"p{boundary}=sym.Symbol('p{boundary}')")
            if not np.isnan(self.pressures[boundary]):
                exec(f"b_term = self.trans[{min(boundary,i)}] * 2 * ((p{boundary} - p{i}) - (self.fluid.g * (self.grid.z[{boundary}]-self.grid.z[{i}])))")
                exec(f"b_term = b_term.subs(p{boundary}, {self.pressures[boundary]})")
            else:
                exec(f"b_term = self.rates[{boundary}]")
            terms.append(locals()['b_term'])

        # 3. Flow from grid well:
        if i in self.wells:
            if 'q' in self.wells[i]:
                terms.append(self.wells[i]['q'])
            else:
                exec(f"w_term = - self.wells[{i}]['G'] / (self.fluid.B*self.fluid.mu) * (p{i} - self.wells[{i}]['pwf'])")
                terms.append(locals()['w_term'])
        if verbose: print('terms:', terms)

        # 4. Accumulation term:
        if self.RHS[i] == 0:
            exec(f"accumulation = 0")
        else:
            try:
                exec(f"accumulation = self.RHS[{i}] * (p{i} - {self.pressures[i]})")
            except:
                raise Exception("Initial pressure (pi) must be specified")
        a_term = locals()['accumulation']

        # 4. Overall grid flow equation:
        i_flow_equation = sym.Eq(sum(terms), a_term)
        if i_flow_equation.lhs.as_coefficients_dict()[1] != 0:
            i_flow_equation = i_flow_equation.simplify()

        # 5. Find lhs and rhs:
        i_lhs = dict(sorted(i_flow_equation.lhs.as_coefficients_dict().items(), key=lambda x: str(x[0])))
        i_rhs = i_flow_equation.rhs

        return i_lhs, i_rhs


    def get_flow_equations(self, verbose=False):
        for i in self.grid.order[self.grid.i_blocks.astype('bool')]:
            i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
            print(f'Grid {i}: {i_lhs}, {i_rhs}')


    def update_matrix(self, i, A, d, sparse, verbose=False):
        # arrays are passed by reference
        n = self.grid.nx
        i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
        i_lhs = list(i_lhs.values())
        d[i-1] = np.array(i_rhs).astype(self.dtype)
        
        if i < 3:
            start = 0
        else: 
            start = i - 2
        
        if i+len(i_lhs) >= n: 
            end = n
        else:
            end = start + len(i_lhs)
            
        if sparse:
            array = A[i-1,start:end].toarray()[0]
        else:
            array = A[i-1,start:end]
        
        assert array.shape[0] == len(i_lhs), 'lhs of grid {} does not match the coefficients matrix'.format(i)
        A[i-1,start:end] = i_lhs # pass by reference
    

    def get_matrix(self, sparse, plot=False, verbose=False):
        
        if self.grid.D == 1:
            n = self.grid.nx
            if all(self.RHS == 0):
                d = np.zeros((n,1), dtype=self.dtype)
            else: 
                try:
                    d = (-self.RHS[1:-1] * self.pressures[1:-1]).reshape(-1, 1)
                except:
                    raise Exception("Initial pressure (pi) must be specified")
            if sparse:
                A = ss.diags([-self.trans*2-self.RHS[:-1], self.trans[1:-1], self.trans[:-2]],
                             [0, 1, -1], 
                             shape=(n,n),
                             format='lil', #“dia”, “csr”, “csc”, “lil”
                             dtype=self.dtype) 
                # d = ss.lil_matrix((n,1), dtype=self.dtype) #ss.csc_matrix(d)
                d = ss.lil_matrix(d, dtype=self.dtype)
            else:
                A = np.zeros((n,n), dtype=self.dtype)
                # d = np.zeros((n,1), dtype=self.dtype)
                # fill coefficient matrix tri-diagonals with trans values:
                np.fill_diagonal(A, -self.trans*2-self.RHS[:-1])
                rng = np.arange(n-1)
                A[rng, rng+1] = self.trans[1:-1] # east trans for interior blocks To Do: confirm if this is right
                A[rng+1, rng] = self.trans[:-2] # west trans for interior blocks except last.To Do: confirm if this is right

            # Update matrix if there is pressure or flow in 'west' boundary:
            if not np.isnan(self.pressures[0]) or self.rates[0] != 0:
                self.update_matrix(1, A, d, sparse, verbose)
            
            # Update matrix if there is pressure or flow in 'east' boundary:
            if not np.isnan(self.pressures[-1]) or self.rates[-1] != 0:
                self.update_matrix(self.grid.nx, A, d, sparse, verbose) # at last grid: self.grid.nx or -2

            # Update matrix in wells i_blocks:
            for i in self.wells.keys():
                self.update_matrix(i, A, d, sparse, verbose)
            
            if plot: 
                if sparse:
                    plt.imshow(A.toarray())
                else:
                    plt.imshow(A)
                plt.show()
            
            if verbose:
                if sparse:
                    print('- A:\n', A.toarray()), print('- d:\n', d.toarray())
                else:
                    print('- A:\n', A), print('- d:\n', d)

            return A, d


    def solve(self, sparse=True, check_MB=True, update=True, verbose=False):
        self.A, self.d = self.get_matrix(sparse, verbose=verbose)

        if sparse:
            pressures = ssl.spsolve(self.A.tocsc(), self.d).flatten()
        else:
            pressures = np.linalg.solve(self.A, self.d).flatten() # same as: np.dot(np.linalg.inv(A),d)
        
        if self.grid.D == 1:
            # Update pressures:
            if update:
                self.pressures[1:-1] = pressures

                # Update wells:
                for i in self.wells.keys():
                    if 'q' in self.wells[i]:
                        self.wells[i]['pwf'] = self.pressures[i] + \
                                                (self.wells[i]['q']*self.fluid.B*self.fluid.mu/self.wells[i]['G'])
                    if 'pwf' in self.wells[i]:
                        self.wells[i]['q'] =  - self.wells[i]['G'] / (self.fluid.B*self.fluid.mu) * \
                                                (self.pressures[i] - self.wells[i]['pwf'])
                        self.rates[i] = self.wells[i]['q']
                # Update boundaries:
                for boundary in self.grid.boundaries:
                    i = 1 if boundary==0 else boundary-1
                    if not np.isnan(self.pressures[boundary]):
                        self.rates[boundary] = self.trans[min(i,boundary)] * 2 * ((self.pressures[boundary] - self.pressures[i]) - (self.fluid.g * (self.grid.z[boundary]-self.grid.z[i])))
                    else:
                        pass
        
        if check_MB:
            self.check_MB(verbose)
        
        if verbose: print('- Pressures:\n', self.pressures), print('- rates:\n', self.rates)

        return pressures
    
    
    def check_MB(self, verbose, error=0.1):
        """Material Balance Check
        
        """
        if self.comp_type == 'incompressible':
            self.error = self.rates.sum() # must add up to 0

        if self.comp_type == 'compressible':
            # Check MB error over a time step: 
            self.incremental_error = self.RHS.sum() / self.rates.sum()
            # Check MB error from initial state to current time step: (less accurate)
            # self.cumulative_error = self.RHS.sum() * self.dt / (self.rates.sum())
            self.error = abs(self.incremental_error - 1)

        assert abs(self.error) < error, 'Material balance error ({}) higher than the allowed error ({}).'.format(self.error, error)
        if verbose: print(f"- Solved with an acceptable MB Error of {self.error} {self.units['error']}")


    def plot(self, property:str):
        exec(f"plt.imshow(self.{property}[1:-1][np.newaxis, :])")
        plt.colorbar(label=f'{property} ({self.units[property[:-1]]})')
        plt.yticks([])
        plt.xlabel('Grid (i)')
        plt.xticks(ticks=range(0, 4), labels=range(1, 5))
        plt.show()


    def show_grid(self, property:str, show_centers=True, show_boundary=False, show_bounds=False):
        plots.show_grid(self, property, show_centers, show_boundary, show_bounds)

    # https://stackoverflow.com/questions/48338847/how-to-copy-a-python-class-instance-if-deepcopy-does-not-work
    def copy(self):
        copy_model = Model(grid=self.grid, fluid=self.fluid, pi=self.pi, 
                            dt=self.dt, dtype=self.dtype, unit=self.unit)
        # for w in self.wells:
        #     well = wells.Well(self.wells[w])
        #     copy_model.set_well(well)
        # copy_model.set_boundaries(self.b_dict)
        return copy_model


if __name__ == '__main__':
    grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=400, phi=0.27, k=270, dtype='double')
    grid.dx[2]=600
    #grid.get_pv_grid()
    fluid = fluids.Fluid(mu=0.5 , B=1, dtype='double')
    model = Model(grid, fluid, dtype='single')
    model.set_well(i=4, q=-600, s=1.5, r=3.5) # model.set_well(i=4, wells.Well(i=4, q=-600, s=0, r=3.5))
    model.set_well(i=1, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
    # print(model.pressures)
    # A, d = model.get_matrix(sparse=True)
    # print(A.toarray())
    # print(d.toarray())
    # A, d = model.get_matrix(sparse=False)
    # print(A)
    # print(d)
    # model.get_flow_equations()
    # print(model.get_matrix(sparse=False)[0])
    # print(model.get_matrix(sparse=True)[0].toarray())
    model.solve(sparse=False, check_MB=True, verbose=False)
    #model.plot('pressures')
    # model.plot_grid()
    plots.show_grid(model, property='pressures', show_centers=True, show_boundary=False)
    # Comp
    # model.set_compressibility(10) # to do: this method should not be exposed!
    
    # Units: 
    # model.set_units('metric')

    # Reporting:
    # print(model.__doc__)
    # print('- repr', repr(model))
    # model.report()
