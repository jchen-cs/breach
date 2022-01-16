% Simple test program that just checks a spec
clear;
mdl = 'Heat'; % A system with a linear function

Sys = BreachSimulinkSystem(mdl); 
Sys.PrintSignals(); % Should get us an output named "ModelOut"

SysFalsify = Sys.copy();

outdoor_gen = fixed_cp_signal_gen({'OutdoorTemp'}); % Just a constant signal
heat_gen = fixed_cp_signal_gen({'HeaterTemp'});      

InputGen = BreachSignalGen({outdoor_gen, heat_gen});
InputGen.SetParam({'OutdoorTemp_u0', 'HeaterTemp_u0'}, [0, 50]);

SysFalsify.SetInputGen(InputGen);

ComfortableTemp = STL_Formula('ComfortableTemp', 'alw((RoomTemp[t] >= 50) & (RoomTemp[t] <= 80))');

global iterations;
iterations = 0;

SysFalsify.SetParamRanges({'OutdoorTemp_u0', 'HeaterTemp_u0'}, [0 30; 40 50]);
req = BreachRequirement(ComfortableTemp);
req.SetParam('t_start', 1);

falsify = FalsificationProblem(SysFalsify, req);
%falsify.solver = 'simulated_annealing'
%falsify.startSample = []
falsify.solve();
fprintf('%d iterations\n', round(iterations))

cex = falsify.GetFalse();
cex.PlotSignals({'OutdoorTemp', 'HeaterTemp', 'RoomTemp'});