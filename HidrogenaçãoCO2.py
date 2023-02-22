from pysis.flowsheet import Simulation

#Simulação Principal com erro
fs = Simulation(path=r'C:\Users\LuanNogueira\Downloads\Ic-Metamodelagem\DWSIMSIMS\HidrogenaçãoCO2A.hsc')

#Simulação de teste funcionando 
#fs = Simulation(path=r'C:\Users\LuanNogueira\Downloads\Ic-Metamodelagem\DWSIMSIMS\Test1.hsc')

n_s='1'
propToRead={
    'Pressure':'Bar'
}
Read_Properties = fs.MatStreams[n_s].get_properties(propToRead)
print(Read_Properties)