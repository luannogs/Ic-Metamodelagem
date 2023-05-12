from pysis.flowsheet import Simulation
import time
import numpy
from scipy.stats import _qmc
import pysis.flowsheet as sim_f

fs = Simulation(path=f'C:\\Users\\LuanNogueira\Downloads\\Ic-Metamodelagem\DWSIMSIMS\\HidrogenaçãoCO2A.hsc')
fs.update_flowsheet()
fs.set_visible(visibility=1)
operations=fs.Operations
materialStreams=fs.MatStreams
distillationColumn=sim_f.DistillationColumn(fs.Operations['T-100'].COMObject)
Pfr=sim_f.PFR(fs.Operations['PFR-100'].COMObject)
fs.Solver.CanSolve=True


def simulated_values() ->list: 
    # Variáveis de interesse: 
    # Vazões de utilidade/MP
    # Custos Equipamentos 
    # Metanol Produzido 
    
   
    return OperationPower('Rotating Equipment')+[
        
    materialStreams['9'].get_massflow('kg/h'),

    materialStreams['DestL'].get_massflow('kg/h'),

    Pfr.get_volume('m3')

    
    ]

def set_dof(x:list):
    
    fs.Solver.CanSolve=False

    materialStreams['13'].set_temperature(x[0],'C')

    materialStreams['9'].set_molarflow(x[1],'kgmole/h')
    
    fs.Operations['TEE-100'].SplitsValue=(x[2],1-x[2])

    distillationColumn.set_specifications(spec='Reflux Ratio', value=x[3])
    
    Pfr.set_volume(x[4])

    materialStreams['10'].set_pressure(x[5],'kpa')
    materialStreams['8'].set_pressure(x[5],'kpa')
    materialStreams['REC'].set_pressure(x[5],'kpa')
    fs.Solver.CanSolve = True
def RReciclo():
    materialStreams['REC'].set_molarflow(5.204e4)
    materialStreams['REC'].set_compmolarflow(materialStreams['REC-2'].get_compmolarflow())

def run(x:list):
    start=time.time()
    y=[]
    for _ in x:
        set_dof(_)
        trials=0
        is_converged=False
        is_valid=False
        while(is_converged==False or is_valid==False) and trials<=1:
            set_dof(_)
            if trials>0:
                RReciclo()
            is_valid=all([i>-3000 for i in simulated_values()])

            is_converged=distillationColumn.get_convergence()
            
            trials+=1
        y.append(simulated_values()+[is_valid, is_converged])
    sim_f.Simulation.save(fs)
    end=time.time()
    print(end-start)
    return y
    # fs.close()

sampler=_qmc.LatinHypercube(d=6)
sample=sampler.random(n=2)
_qmc.discrepancy(sample)

l_bounds=[180,4500,0.001,1.25,35,4500]

u_bounds=[240,6500,0.05,1.8,55,5500]

x=list(_qmc.scale(sample,l_bounds, u_bounds))

def OperationPower(operation :str):
    PowerList=[]
    for _ in fs.EnerStreams.keys(): 
        for ops in fs.EnerStreams[_].get_connections().values():
            for l in ops:
                if(fs.Operations[str(l)].classification == operation):
                    PowerList.append(fs.EnerStreams[_].get_power())
    return (PowerList)

operations['E-100'].COMObject.AreaValue()

