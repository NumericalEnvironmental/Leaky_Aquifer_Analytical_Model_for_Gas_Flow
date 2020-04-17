###################################################################################
#
# radial_gas_ss.py
# version 1
#
# by Walt McNab (2020)
#
# leaky aquifer-inspired analytical solution for (linearized) ideal gas extraction
# DeGlee (1930, 1951) solution
#
###################################################################################

from numpy import *
from scipy.special import *
from scipy.spatial import distance
from pandas import *
from scipy.optimize import fsolve
from scipy.integrate import quad


# constants

g = 9.807                   # gravitational acceleration
R = 8.314                   # universal gas constant (J K-1 mol-1)
delta = 0.02                # distance offset for gradient calculation (m)
rw = 0.08                   # example well radius (m)


### classes

class Reservoir:

    def __init__(self):
        # aquifer properties
        in_file = open('reservoir.txt','r')
        i = 0
        for line in in_file:
            line_input = line.split()
            if line_input[0] == 'aquifer_h_perm':
                self.k = float(line_input[1])
            elif line_input[0] == 'aquitard_v_perm':
                kc = float(line_input[1])
            elif line_input[0] == 'aquifer_thickness':
                self.b = float(line_input[1])
            elif line_input[0] == 'aquitard_thickness':
                bc = float(line_input[1])                
            else:              # ambient gas pressure
                self.P0 = float(line_input[1])
            i += 1
        in_file.close()
        self.B = sqrt(self.k*self.b*bc/kc)
        
    def RadialFlow(self, r, Qm, gas):
        # leaky steady-sate solution for square of gas pressure perturbation
        dP2 = k0(r/self.B) * gas.u*Qm/(pi*self.k*self.b*gas.beta)
        return dP2

    def PressureDiff(self, x, y, wp, Qm, gas):
        # estimate pressure perturbation across set of points (x & y arrays) from a particular well
        gp = transpose([x, y])          
        rc = distance.cdist(gp, wp)             # distance array w.r.t. to x & y arrays and well location (wp)
        rt = transpose(rc)[0]
        dP2 = self.RadialFlow(rt, Qm, gas)      # square of the pressure associated with well
        P = sqrt(self.P0**2 + dP2)
        return P - self.P0                      # computed pressure minus background (i.e., delta P)

    def Flux(self, grad, gas):
        # compute gas mass fluxes per unit area
        jx = -grad[0] * self.k  * gas.beta / (2*gas.u)
        jy = -grad[1] * self.k  * gas.beta / (2*gas.u)                  
        return jx, jy

    def FluxGrid(self, grid, reservoir, gas):
        # find gas mass fluxes at grid points
        grad_x = ((reservoir.P0+grid.dPx)**2 - (reservoir.P0+grid.dP)**2) / delta
        grad_y = ((reservoir.P0+grid.dPy)**2 - (reservoir.P0+grid.dP)**2) / delta
        grid.jx, grid.jy = self.Flux([grad_x, grad_y], gas)
        return grid


class Polygon:

    def __init__(self):
        # closed polygon to calculate integrated fluxes toward well array
        points = read_csv('polygon.csv')
        self.x = array(points['x'])
        self.y = array(points['y'])
        self.m = []             # slope of polygon side
        self.n = []             # normal vector to polygon side
        self.L = []             # length of polygon side
        for i in range(len(self.x)-1):
            dx = self.x[i+1] - self.x[i]            
            dy = self.y[i+1] - self.y[i]
            self.m.append(dy/dx)
            L = sqrt(dx**2 + dy**2)
            self.L.append(L)
            self.n.append([-dy/L, dx/L])

    def Qn(self, x, well, aiw, reservoir, grid, gas, indx, intake):
        # normal flux vector (kg m^-2 sec^-1) along polygon side      
        y = self.y[indx] + self.m[indx]*(x-self.x[indx])
        dPw = reservoir.PressureDiff(array(well['x']), array(well['y']), array([[x, y]]), array(well['Q']), gas)        
        dPwx = reservoir.PressureDiff(array(well['x']), array(well['y']), array([[x+delta, y]]), array(well['Q']), gas)         
        dPwy = reservoir.PressureDiff(array(well['x']), array(well['y']), array([[x, y+delta]]), array(well['Q']), gas) 
        dP = sum(dPw)
        dPx = sum(dPwx)
        dPy = sum(dPwy)        
        if intake:
            dPaiw = reservoir.PressureDiff(array(aiw['x']), array(aiw['y']), array([[x, y]]), array(aiw['Q']), gas)        
            dPaiwx = reservoir.PressureDiff(array(aiw['x']), array(aiw['y']), array([[x+delta, y]]), array(aiw['Q']), gas)         
            dPaiwy = reservoir.PressureDiff(array(aiw['x']), array(aiw['y']), array([[x, y+delta]]), array(aiw['Q']), gas)
            dP += sum(dPaiw)
            dPx += sum(dPaiwx)
            dPy += sum(dPaiwy)
        grad_x = ((reservoir.P0+dPx)**2 - (reservoir.P0+dP)**2) / delta
        grad_y = ((reservoir.P0+dPy)**2 - (reservoir.P0+dP)**2) / delta                  
        jx, jy = reservoir.Flux([grad_x, grad_y], gas) 
        return -dot([jx, jy], self.n[indx])    
   
    def QTotal(self, well, aiw, reservoir, grid, gas, intake):
        QTot = 0.
        for i in range(len(self.x)-1):
            x0 = min([self.x[i], self.x[i+1]])
            xf = max([self.x[i], self.x[i+1]]) 
            QTot += (self.L[i]/abs(xf-x0)) * reservoir.b * quad(self.Qn, x0, xf, args=(well, aiw, reservoir, grid, gas, i, intake))[0]
        return QTot    
    
    
class Gas:

    def __init__(self):
        # gas properties
        in_file = open('gas.txt','r')
        i = 0
        for line in in_file:
            line_input = line.split()
            if line_input[0] == 'temperature':
                T = float(line_input[1])
            elif line_input[0] == 'viscosity':
                self.u = float(line_input[1])
            else:               # molecular_weight':
                M = float(line_input[1])
            i += 1
        in_file.close()
        self.beta = M/(R*T)    # compressibility factor

        
class Grid:

    def __init__(self):
        # create a uniform grid with supplied specs and convert to a n x 2 arrays
        in_file = open('grid.txt','r')
        i = 0
        for line in in_file:
            line_input = line.split()
            if line_input[0] == 'x0':
                x0 = float(line_input[1])
            elif line_input[0] == 'xf':
                xf = float(line_input[1])
            elif line_input[0] == 'y0':
                y0 = float(line_input[1])
            elif line_input[0] == 'yf':
                yf = float(line_input[1])                
            elif line_input[0] == 'num_x':
                num_x = int(line_input[1])                
            else:               # grid discretization along y-direction
                num_y = int(line_input[1])
            i += 1
        in_file.close()        
        xGrid = linspace(x0, xf, num_x+1)
        yGrid = linspace(y0, yf, num_y+1)
        X, Y = meshgrid(xGrid,yGrid)
        self.x = X.flatten()
        self.y = Y.flatten()
        self.dP = zeros(len(self.x), float)
        self.dPx = zeros(len(self.x), float)   
        self.dPy = zeros(len(self.x), float)           
        self.jx = zeros(len(self.x), float)
        self.jy = zeros(len(self.x), float)        


### utility functions   

def AirOpt(Qc, Pw, aiw, reservoir, gas, xaiw, yaiw):
    # optimize flow rates in air injection wells to reduce local vacuum to zero
    Pc = zeros(len(Pw), float)      # Pc = calibrated pressure
    for i in range(len(aiw)):    
        xw = aiw['x'].iloc[i]
        yw = aiw['y'].iloc[i]
        wp = array([[xw, yw]])
        dP = reservoir.PressureDiff(xaiw, yaiw, wp, Qc[i], gas)
        Pc += dP
    return Pc + Pw   
 
    
def AiwFluxes(aiw, well, reservoir, gas):
    # air intake well fluxes
    xaiw = array(aiw['x']) + rw
    yaiw = array(aiw['y'])
    Qc0 = array(aiw['Q'])
    Pw = zeros(len(xaiw), float)
    for i in range(len(well)):    # initial pressures at air intake well locations
        xw = well['x'].iloc[i]
        yw = well['y'].iloc[i]
        Qm = well['Q'].iloc[i]
        wp = array([[xw, yw]])
        dP = reservoir.PressureDiff(xaiw, yaiw, wp, Qm, gas)
        Pw += dP
    Qc = fsolve(AirOpt, Qc0, args=(Pw, aiw, reservoir, gas, xaiw, yaiw))
    return Qc    


def WellPressures(intake, well, aiw, reservoir, gas):
    ### compute total delta pressure at each well in model
    if intake:
        namePts = concatenate((array(well['name']), array(aiw['name'])))
        xPts = concatenate((array(well['x']), array(aiw['x']))) + rw
        yPts = concatenate((array(well['y']), array(aiw['y'])))
    else:
        namePts = array(well['name'])
        xPts = array(well['x']) + rw
        yPts = array(well['y'])        
    Pr = zeros(len(namePts), float)    
    for i in range(len(well)):    
        xw = well['x'].iloc[i]
        yw = well['y'].iloc[i]
        Qm = well['Q'].iloc[i]
        wp = array([[xw, yw]])
        dP = reservoir.PressureDiff(xPts, yPts, wp, Qm, gas)
        Pr += dP
    if intake:
        for i in range(len(aiw)):    
            xw = aiw['x'].iloc[i]
            yw = aiw['y'].iloc[i]
            Qm = aiw['Q'].iloc[i]
            wp = array([[xw, yw]])
            dP = reservoir.PressureDiff(xPts, yPts, wp, Qm, gas)
            Pr += dP
    return namePts, Pr


def PressureGrid(wellSet, grid, reservoir, gas):
    # update grid values for total vacuum and to inform flux calculations
    for i in range(len(wellSet)):    
        xw = wellSet['x'].iloc[i]
        yw = wellSet['y'].iloc[i]
        Qm = wellSet['Q'].iloc[i]
        wp = array([[xw, yw]])
        dP = reservoir.PressureDiff(grid.x, grid.y, wp, Qm, gas)
        grid.dP += dP
        dPx = reservoir.PressureDiff(grid.x+delta, grid.y, wp, Qm, gas)
        grid.dPx += dPx
        dPy = reservoir.PressureDiff(grid.x, grid.y+delta, wp, Qm, gas)
        grid.dPy += dPy    
    return grid


### main script

def GasFlow(intake, budget):

    gas = Gas()
    print('Read gas data.')
    
    reservoir = Reservoir()
    print('Read reservoir data.')
    print('Note:', 'B = ', reservoir.B, '3*b = ', 3*reservoir.b)
    print('\t', 'Validity = ', reservoir.B > 3*reservoir.b)     # condition required of DeGlee model
    
    well = read_csv('wells.csv')
    print('Read SVE well data.')
    if intake:
        aiw = read_csv('air intake wells.csv')
        print('Read air injection well data.')    
    
    grid = Grid()
    print('Set up grid.')

    if budget:
        polygon = Polygon()
        print('Read surrounding polygon.')

    if intake:
        Q = AiwFluxes(aiw, well, reservoir, gas)
        aiw['Q'] = Q        
        print('Computed fluxes in air intake wells.')
    else: aiw = []

    # vacuum summary
    namePts, Pr = WellPressures(intake, well, aiw, reservoir, gas)
    print('Well', 'Vacuum (Pa)')
    for i in range(len(namePts)):
        print(namePts[i], Pr[i])
    refP = min(Pr)                  # largest vacuum in model

    ### populate grid for plotting
    grid = PressureGrid(well, grid, reservoir, gas)
    if intake: grid = PressureGrid(aiw, grid, reservoir, gas)
    print('Populated grid points.')

    grid = reservoir.FluxGrid(grid, reservoir, gas)     # add local gas fluxes to grid
    ratio = grid.dP/refP                        # calculate ratio of local vacuum to largest vacuum
    output = DataFrame(data = {'x':grid.x, 'y':grid.y, 'vacuum':-grid.dP, 'jx':grid.jx, 'jy':grid.jy, 'ratio':ratio})
    output.to_csv('output.csv', index=False)

    if budget:
        Qb = polygon.QTotal(well, aiw, reservoir, grid, gas, intake)
        print('Analyzed flow budget:')
        print('\t', 'Bounding polygon = ', Qb)
        if intake:
            Qaiw = aiw['Q'].sum()
            print('\t', 'Air inflow wells = ', Qaiw)            
        else: Qaiw = 0.
        Qwell = well['Q'].sum()
        print('\t', 'SVE wells = ', Qwell)
        QL = -Qwell - Qb - Qaiw
        print('\t', 'Leakage = ', QL)
        
    print('Done.')


### run script

intake = True            # intake = True --> air intake wells considered; intake == False --> not considered
budget = True           # budget = True --> integrate mass flux along boundary polygon and sum fluxes
GasFlow(intake, budget)
