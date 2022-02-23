from tabulate import tabulate

class Base():
    name = None
    units_dict = {
        'field':
            {
                'transmissibility': 'STB/D-psi',
                'error': 'STB/D',
                'pressure': 'psia',
            },
        'metric':
            {
                'transmissibility': 'M3/D-bar',
                'error': 'M3/D'
            }
    }
    factors_dict = {
        'field':
            {
                'transmissibility conversion': 0.001127,
                'gravity conversion': 0.21584*10**-3,
                'gravitational acceleration': 32.174,
            },
        'metric':
            {
                'transmissibility conversion': 0,
                'gravity conversion': 0,
                'gravitational acceleration': 0,
            }
        
    }
    
    unit = 'field'
    units = units_dict[unit]
    factors = factors_dict[unit]

    def __init__(self, unit='field'):
        self.set_units(unit)


    def set_units(self, unit='field'):
        if unit in self.units_dict.keys():
            self.unit = unit
            self.units = self.units_dict[unit]
            self.factors = self.factors_dict[unit]
        else:
            raise ValueError(f"The selected unit system ({unit}) is unknown!")


    def set_compressibility(self, comp):
        """
        
        """
        self.compressibility = self.comp = comp
        if comp == 0:
            self.compressibility_type = self.comp_type = 'incompressible'
            #self.__set_compressibility_type('incompressible')
        elif comp > 0:
            self.compressibility_type = self.comp_type = 'compressible'
            #self.__set_compressibility_type('compressible')
        else: 
            # self.__set_compressibility_type('unknown')
            raise ValueError('Compressibility type is unknown!')
    set_comp = set_compressibility


    # def __set_compressibility_type(self, comp_type: str):
    #     """
        
    #     """
    #     if 'incomp' in comp_type:
    #         self.compressibility_type = self.comp_type = 'incompressible'
    #     elif 'comp' in comp_type:
    #         self.compressibility_type = self.comp_type = 'compressible'
    #     else:
    #         raise ValueError('Compressibility type is unknown!')
    # __set_comp_type = __set_compressibility_type


    def report(self, prop=None, ifmt=0):
        props = vars(self)
        tablefmt = ["pipe", "plain", "simple", "fancy_grid", "presto", "orgtbl"]
        if prop == None:
            print(f'{self.name} Information: \n')
            table = tabulate([(str(k),str(v)) for k, v in props.items()], 
                        headers=['Property', 'Value'], 
                        showindex="always", 
                        tablefmt=tablefmt[ifmt]
                    )
            print(table)
            print(' ')
        elif prop in props.keys():
                print(f' - {prop}: {props[prop]}')
        else:
            print('Unknown property! Use prop=None to print all available properties.')

            
    def __str__(self):
        self.report()
        return ''
    

    def __repr__(self): # print(repr(class))
        return str(vars(self))
        

if __name__ == '__main__':
    b = Base('metric')
    b.name = 'b'
    print(repr(b))
    k = Base()
    k.name = 'k'
    b.report()
    k.report()
