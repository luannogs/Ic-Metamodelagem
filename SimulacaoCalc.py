import pysis
import numpy
from scipy.stats import _qmc
import json
import PyPinch.src.PyPinch as pinch
import csv
from pina import PinchAnalyzer, make_stream


def LHS(n):
    sampler = _qmc.LatinHypercube(d=6)
    sample = sampler.random(n)
    _qmc.discrepancy(sample)
    l_bounds = [180, 4500, 0.001, 1.25, 35, 4500]
    u_bounds = [240, 6500, 0.05, 1.8, 55, 5500]
    x = (_qmc.scale(sample, l_bounds, u_bounds).tolist())
    
    lhs= dict(zip([int(i) for i in range(0,n)],x))
    return lhs


class SimulationOps:
    def __init__(self, path):
        self.path = path
        self.fs = pysis.Simulation(path)
        # self.Operations = self.fs.Operations
        # self.fs.MatStreams = self.fs.MatStreams
        # distillationColumn = 0
        self.path = path
        self.fs.set_visible(1)
        pass

    def sampling_sim(self, samples:dict, iterations_to_reset:int,file_y:str):
        y = {}
        results={}     
        f=0
        
        try:
            with open(file_y,'r') as file:
                d=(json.load(file))
                file.close()
            keys_off=[str(i) for i in list(samples.keys()) if int(i) not in [int(i) for i in list(d.keys())]]
            samples_d=dict(zip(keys_off,[samples[i] for i in keys_off]))
        except:
            samples_d=samples

        for keys,lists in samples_d.items():

            try:
                self.set_dof(x=list(lists))
                trials = 0
                is_converged = pysis.flowsheet.DistillationColumn(
                    self.fs.Operations["T-100"].COMObject
                ).get_convergence()
                stream_convergeance = self.stream_convergeance()

                while (
                    is_converged == False or stream_convergeance == False
                ) and trials <= 3:
                    self.set_dof(x=lists)
                    if trials > 0:
                        try:
                            self.RReciclo()
                            
                        except:
                            try:
                                self.reset_stream()
                            except:
                                None
                            

                        else:
                            print("Resetando Reciclo e correntes")
                            self.set_dof(x=list(lists))

                    stream_convergeance = self.stream_convergeance()
                    is_converged = pysis.flowsheet.DistillationColumn(
                        self.fs.Operations["T-100"].COMObject
                    ).get_convergence()
                    trials += 1
                

                self.pinch_analysis()
                results={}
                
                results.update(self.UCEquipamentos())
                results.update(self.DesignVariables())
                results.update(self.CustoMP())
                results.update(self.CustoUtilidades())
                results.update({'Methanol':self.fs.MatStreams["DestL"].get_compmassflow()["Methanol"]})
                results.update({"Lucro":self.profit()})
                results.update({'ConvergênciaTD':pysis.flowsheet.DistillationColumn(
                        self.fs.Operations["T-100"].COMObject
                    ).get_convergence()})
                results.update(self.QntdUtilidades())
                results.update({'Convergencia_correntes':self.stream_convergeance()})
                
                y=({keys:results})

                try:
                    with open(file_y,'r') as json_file:
                        data=dict(json.load(json_file))
                        json_file.close()
                    data.update(y)
                    with open(file_y,'w') as json_file:
                        json.dump(data,json_file)
                        json_file.close()
                except:
                    with open(file_y,'w') as json_file:
                        json.dump(y,json_file)
                        json_file.close()

                f+=1
                print("OK: {}".format([str(keys)]))    
                if f%iterations_to_reset==0:
                    try:
                        print('reset_simulação')
                        self.fs.close()
                        self.fs = 0
                        self.fs = pysis.Simulation(self.path)
                        self.fs.set_visible(1)
                    except:
                        print('Falha Reset Simulação')
                        None
                
            except:
                print("Erro: {}".format([(keys)]))
                try:
                    self.fs.close()
                    self.fs = 0
                    self.fs = pysis.Simulation(self.path)
                    self.fs.set_visible(1)
                except:
                    self.fs = pysis.Simulation(self.path)
                    self.fs.set_visible(1)
                    print('Rerun')
                    
                    
                    self.sampling_sim(samples,iterations_to_reset,file_y)



    def DesignVariables(self):
        
        return {
            "VolumeReator": pysis.flowsheet.PFR(
                self.fs.Operations["PFR-100"].COMObject
            ).get_volume(),
            "AreaTC": [
                self.fs.Operations[u].COMObject.UAValue * 1000 / 1700
                for u in self.fs.Operations.keys()
                if self.fs.Operations[u].classification == "Heat Transfer Equipment"
                and any(
                    [
                        x == "HeatTransferArea"
                        for x in self.fs.Operations[u].COMObject.__dir__()
                    ]
                )
            ],
            "VazaoWasteWater": self.fs.MatStreams[
                "LIQ"
            ].COMObject.IdealLiquidVolumeFlowValue
            * 1000,
            "CaldeiraVazao": self.fs.MatStreams["UtilidadeQuente"].get_massflow(),
            "VolumePackingTD": self.fs.Operations["CustoDestilação"]
            .COMObject.Cell("A1")
            .CellValue,
            "PotenciaCompressores": [
                self.fs.Operations[u].COMObject.EnergyValue
                for u in self.fs.Operations.keys()
                if self.fs.Operations[u].classification == "Rotating Equipment" and u != 'K-105'
            ],
            "MassaProcessVessel": [
                self.fs.Operations["ProcessVessel"].COMObject.Cell("B7").CellValue,
                self.fs.Operations["ProcessVessel"].COMObject.Cell("B8").CellValue,
            ],
        }
    


    
    def pinch_analysis_plot(self):
        
        ps1=make_stream(self.fs.Operations['E-100'].COMObject.DutyValue,self.fs.MatStreams["2"].get_temperature('C'),self.fs.MatStreams["3"].get_temperature('C'))
        ps2=make_stream(self.fs.Operations['E-101'].COMObject.DutyValue,self.fs.MatStreams["4"].get_temperature('C'),self.fs.MatStreams["5"].get_temperature('C'))
        ps3=make_stream(self.fs.Operations['E-102'].COMObject.DutyValue,self.fs.MatStreams["6"].get_temperature('C'),self.fs.MatStreams["7"].get_temperature('C'))
        ps4=make_stream(-self.fs.Operations['E-103'].COMObject.DutyValue,self.fs.MatStreams["12"].get_temperature('C'),self.fs.MatStreams["13"].get_temperature('C'))
        ps5=make_stream(self.fs.Operations['E-104'].COMObject.DutyValue,self.fs.MatStreams["14"].get_temperature('C'),self.fs.MatStreams["15"].get_temperature('C'))
        ps6=make_stream(self.fs.Operations['E-105'].COMObject.DutyValue,self.fs.MatStreams["22"].get_temperature('C'),self.fs.MatStreams["23"].get_temperature('C'))
        # pscond=make_stream(self.fs.Operations['Condensador'].COMObject.DutyValue,self.fs.MatStreams["35"].get_temperature('C'),self.fs.MatStreams["40"].get_temperature('C'))
        # psref=make_stream(-self.fs.Operations['Refervedor'].COMObject.DutyValue,self.fs.MatStreams["43"].get_temperature('C'),self.fs.MatStreams["44"].get_temperature('C'))

        min_temp_diff=10
        temp_shift=min_temp_diff/2

        analyzer=PinchAnalyzer(temp_shift)
        analyzer.add_streams(ps1,ps2,ps3,ps4,ps5,ps6)

        return analyzer.grand_composite_curve

    
    def listDesignVariables(self):
        return numpy.array(list(self.DesignVariables().values())).flatten()
    
    def UCEquipamentos(self):
        design_variables = self.DesignVariables()

        return {
            "CustoReator": (4.00
            * (61500 + 32500 * (design_variables["VolumeReator"]) ** 0.80))*813/499.6,
            # "CustoTrocadoresDeCalor": (sum(
            #     [1.0 * (28000 + 54 * u**1.20) for u in design_variables["AreaTC"]]
            # ))*813/499.6,
            "CustoTorreResfriamento": (2.5
            * (170000 + 1500 * design_variables["VazaoWasteWater"] ** 0.90))*813/499.6,
            "CustoDestilação": (4000 * design_variables["VolumePackingTD"]+4.0*(17400+79*(20*3.842*15e-3*7890)**0.85))*813/499.6,
            "CustoComp":( sum(
                [
                    (8400 + 3100 * (u) ** 0.6)
                    for u in design_variables["PotenciaCompressores"]
                ]
            ))*813/499.6,
            "CustoProcessVessel": (sum(
                [
                    4.0 * (17400 + 79 * u** 0.85) 
                    for u in design_variables["MassaProcessVessel"]
                ]
            ))*813/499.6,
        }

    def CustoEquipamentos(self):
        return sum(list(self.UCEquipamentos().values()))

    def QntdUtilidades(self):
        # self.pinch_analysis()
        return {
            "UtilidadeFria":self.fs.MatStreams['UtilidadeFria-2'].COMObject.IdealLiquidVolumeFlowValue*3600 ,
            "UtilidadeQuente":self.fs.EnerStreams['UQE-2'].get_power()/1e6*3600,
            "VazaoWasteWater": self.fs.MatStreams[
                "LIQ"
            ].COMObject.IdealLiquidVolumeFLowValue
            * 3600,
            "PotenciaEletricidade": self.fs.EnerStreams["CustoEletricidade"].get_power()
            / 1e3,
        }

    def pinch_analysis(self):
        ps1=make_stream(self.fs.Operations['E-100'].COMObject.DutyValue,self.fs.MatStreams["2"].get_temperature('C'),self.fs.MatStreams["3"].get_temperature('C'))
        ps2=make_stream(self.fs.Operations['E-101'].COMObject.DutyValue,self.fs.MatStreams["4"].get_temperature('C'),self.fs.MatStreams["5"].get_temperature('C'))
        ps3=make_stream(self.fs.Operations['E-102'].COMObject.DutyValue,self.fs.MatStreams["6"].get_temperature('C'),self.fs.MatStreams["7"].get_temperature('C'))
        ps4=make_stream(-self.fs.Operations['E-103'].COMObject.DutyValue,self.fs.MatStreams["12"].get_temperature('C'),self.fs.MatStreams["13"].get_temperature('C'))
        ps5=make_stream(self.fs.Operations['E-104'].COMObject.DutyValue,self.fs.MatStreams["14"].get_temperature('C'),self.fs.MatStreams["15"].get_temperature('C'))
        ps6=make_stream(self.fs.Operations['E-105'].COMObject.DutyValue,self.fs.MatStreams["22"].get_temperature('C'),self.fs.MatStreams["23"].get_temperature('C'))
        # pscond=make_stream(self.fs.Operations['Condensador'].COMObject.DutyValue,self.fs.MatStreams["35"].get_temperature('C'),self.fs.MatStreams["40"].get_temperature('C'))
        # psref=make_stream(-self.fs.Operations['Refervedor'].COMObject.DutyValue,self.fs.MatStreams["43"].get_temperature('C'),self.fs.MatStreams["44"].get_temperature('C'))

        min_temp_diff=10
        temp_shift=min_temp_diff/2

        analyzer=PinchAnalyzer(temp_shift)
        analyzer.add_streams(ps1,ps2,ps3,ps4,ps5,ps6)
        
        self.fs.Operations['E-109-2'].COMObject.DutyValue=float(analyzer.cold_utility_target)
        self.fs.Operations['E-108-2'].COMObject.DutyValue=float(analyzer.hot_utility_target)
        pass  

    def pinch_analysis_comparer(self):
        ps1=make_stream(self.fs.Operations['E-100'].COMObject.DutyValue,self.fs.MatStreams["2"].get_temperature('C'),self.fs.MatStreams["3"].get_temperature('C'))
        ps2=make_stream(self.fs.Operations['E-101'].COMObject.DutyValue,self.fs.MatStreams["4"].get_temperature('C'),self.fs.MatStreams["5"].get_temperature('C'))
        ps3=make_stream(self.fs.Operations['E-102'].COMObject.DutyValue,self.fs.MatStreams["6"].get_temperature('C'),self.fs.MatStreams["7"].get_temperature('C'))
        ps4=make_stream(-self.fs.Operations['E-103'].COMObject.DutyValue,self.fs.MatStreams["12"].get_temperature('C'),self.fs.MatStreams["13"].get_temperature('C'))
        ps5=make_stream(self.fs.Operations['E-104'].COMObject.DutyValue,self.fs.MatStreams["14"].get_temperature('C'),self.fs.MatStreams["15"].get_temperature('C'))
        ps6=make_stream(self.fs.Operations['E-105'].COMObject.DutyValue,self.fs.MatStreams["22"].get_temperature('C'),self.fs.MatStreams["23"].get_temperature('C'))
        # pscond=make_stream(self.fs.Operations['Condensador'].COMObject.DutyValue,self.fs.MatStreams["35"].get_temperature('C'),self.fs.MatStreams["40"].get_temperature('C'))
        # psref=make_stream(-self.fs.Operations['Refervedor'].COMObject.DutyValue,self.fs.MatStreams["43"].get_temperature('C'),self.fs.MatStreams["44"].get_temperature('C'))

        min_temp_diff=10
        temp_shift=min_temp_diff/2

        analyzer=PinchAnalyzer(temp_shift)
        analyzer.add_streams(ps1,ps2,ps3,ps4,ps5,ps6)
        
        return float(analyzer.cold_utility_target),float(analyzer.hot_utility_target)
        

    def CustoUtilidades(self):
        qntdutilidades = self.QntdUtilidades()
        return {
            "CustoResfriamento": qntdutilidades['UtilidadeFria'] * 0.03,
            "CustoAquecimento": qntdutilidades['UtilidadeQuente'] * 12.1,
            "CustoWasteWater": qntdutilidades["VazaoWasteWater"] * 1.5,
            "CustoEletricidade": qntdutilidades["PotenciaEletricidade"] * 94.5,
        }

    def CustoTotalUtilidades(self):
        return sum(list(self.CustoUtilidades().values()))

    def VazaoMP(self):
        return {
            "QntdH2": self.fs.MatStreams["9"].COMObject.IdealLiquidVolumeFlowValue
            * 3600,
            # ton
            "QntdCO2": self.fs.MatStreams["1"].get_massflow() / 1000,
        }

    def CustoMP(self):
        return {
            "CustoH2": self.VazaoMP()["QntdH2"] * 5.240,
            "CustoCO2": self.VazaoMP()["QntdCO2"] * 95.50,
        }

    def CustoTotalMP(self):
        return sum(list(self.CustoMP().values()))

    def obj_f(self):
        return [
            self.CustoEquipamentos(),
            self.CustoTotalMP(),
            self.CustoTotalUtilidades(),
            self.fs.MatStreams["DestL"].get_compmolarflow()["Methanol"],
        ]

    def set_dof_u(self, x:list):

        try:
            self.set_dof(x)

            
            
            trials = 0
            is_converged = pysis.flowsheet.DistillationColumn(
                self.fs.Operations["T-100"].COMObject
            ).get_convergence()
            stream_convergeance = self.stream_convergeance()

            trials=0

            while (
                is_converged == False or stream_convergeance == False
            ) and trials <= 3:
                
                self.set_dof(x)
                
                if trials > 0:
                    try:
                        self.RReciclo()
                        
                    except:
                        try:
                            self.reset_stream()
                        except:
                            None
                    else:
                        print("Resetando Reciclo e correntes")
                        self.set_dof(x)

                stream_convergeance = self.stream_convergeance()
                is_converged = pysis.flowsheet.DistillationColumn(
                    self.fs.Operations["T-100"].COMObject
                ).get_convergence()
                trials += 1
            if (
                is_converged == False or stream_convergeance == False
            ):
                print('Não Convergência')
            # self.pinch_analysis()
        except:
            self.reset_simulation()
    def set_dof(self, x):
        self.fs.Solver.CanSolve = False

        self.fs.MatStreams["13"].set_temperature(x[0], "C")

        self.fs.MatStreams["9"].set_molarflow(x[1], "kgmole/h")

        self.fs.Operations["TEE-100"].SplitsValue = (x[2], 1 - x[2])

        pysis.flowsheet.DistillationColumn(
            self.fs.Operations["T-100"].COMObject
        ).set_specifications(spec="Reflux Ratio", value=x[3])

        pysis.flowsheet.PFR(self.fs.Operations["PFR-100"].COMObject).set_volume(x[4])

        self.fs.MatStreams["10"].set_pressure(x[5], "kpa")
        self.fs.MatStreams["8"].set_pressure(x[5], "kpa")
        self.fs.MatStreams["20"].set_pressure(x[5], "kpa")

        self.fs.MatStreams["REC-2"].set_pressure(x[5], "kpa")

        
        self.fs.Solver.CanSolve = True
        
    def func(self, x:list):

        try:
            self.set_dof(x)

            
            
            trials = 0
            is_converged = pysis.flowsheet.DistillationColumn(
                self.fs.Operations["T-100"].COMObject
            ).get_convergence()
            stream_convergeance = self.stream_convergeance()

            trials=0

            while (
                is_converged == False or stream_convergeance == False
            ) and trials <= 3:
                
                self.set_dof(x)
                
                if trials > 0:
                    try:
                        self.RReciclo()
                        
                    except:
                        try:
                            self.reset_stream()
                        except:
                            None
                    else:
                        print("Resetando Reciclo e correntes")
                        self.set_dof(x)

                stream_convergeance = self.stream_convergeance()
                is_converged = pysis.flowsheet.DistillationColumn(
                    self.fs.Operations["T-100"].COMObject
                ).get_convergence()
                trials += 1
            if (
                is_converged == False or stream_convergeance == False
            ):
                print('Não Convergência')
            # self.pinch_analysis()
        except:
            self.reset_simulation()

        return self.profit()

    def RReciclo(self):
        self.fs.MatStreams["REC"].set_molarflow(
            self.fs.MatStreams["REC-2"].get_molarflow()
        )
        self.fs.MatStreams["REC"].set_compmolarflow(
            self.fs.MatStreams["REC-2"].get_compmolarflow()
        )
        

    def reset_stream(self):
        unconverged_streams = [
            i
            for i in list(self.fs.MatStreams.keys())
            if self.fs.MatStreams[i].COMObject.MassFlow() < 0
            or self.fs.MatStreams[i].COMObject.MolarEntropy() < 0
        ]
        streams_pressure = dict(
            zip(
                list(self.fs.MatStreams.keys()),
                [
                    self.fs.MatStreams[u].get_pressure()
                    for u in list(self.fs.MatStreams.keys())
                ],
            )
        )
        unconverged_streams_bool = any(
            [
                self.fs.MatStreams[i].COMObject.MassFlow() < 0
                or self.fs.MatStreams[i].COMObject.MolarEntropy() < 0
                for i in list(self.fs.MatStreams.keys())
            ]
        )
        for i in unconverged_streams:
            if self.fs.MatStreams[i].COMObject.CanModifyStream == True:
                self.fs.MatStreams[i].set_pressure(100)
                self.fs.MatStreams[i].set_pressure(streams_pressure[i])

    def reset_simulation(self):
        self.fs.close()
        self.fs = 0
        self.fs = pysis.Simulation(self.path)
        self.fs.set_visible(1)
        self.fs.Solver.CanSolve = True
        return self

    def stream_convergeance(self):
        stream_convergeance = all(
            [
                i >= 0
                for i in [f.get_massflow() for f in self.fs.MatStreams.values()]
                + [f.get_power() for f in self.fs.EnerStreams.values()]
            ]
        )
        return stream_convergeance

    def ColumnConvergeance(self):
        is_converged = pysis.flowsheet.DistillationColumn(
            self.fs.Operations["T-100"].COMObject
        ).get_convergence()
        return is_converged

    def profit(self):
        # self.pinch_analysis()
        opex = (self.CustoTotalUtilidades() + self.CustoTotalMP()) * 24 * 30 * 12
        capex = (
            self.CustoEquipamentos()
            + self.CustoEquipamentos() * 0.35
            + (self.CustoEquipamentos() + self.CustoEquipamentos() * 0.35) * 0.20
            + (self.CustoEquipamentos() + self.CustoEquipamentos() * 0.35) * 0.30
            + (self.CustoEquipamentos() + self.CustoEquipamentos() * 0.35) * 0.15
        )

        tac = (capex * 0.20 + opex) / 12 / 30 / 24
        methanol_sell = float(
            self.fs.MatStreams["DestL"].get_compmassflow()["Methanol"] * 575 / 1000
        )
        return (methanol_sell - tac)
