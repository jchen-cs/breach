% Simple test program that just checks a spec
clear;
mdl = 'TestSystem'; % A system with a constant 1

Sys = BreachSimulinkSystem(mdl); 
Sys.PrintSignals(); % Should get us an output named "ModelOut"

AboveZero = STL_Formula('AboveZero', 'alw(ModelOut[t] > 0)');

Time = 0:0.05:10;
Sys.Sim(Time);
Sys.CheckSpec(AboveZero);

