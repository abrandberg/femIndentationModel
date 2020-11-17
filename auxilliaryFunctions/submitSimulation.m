function A = submitSimulation(simLoc,simName,submitLine,ctrl)

cd(simLoc)
system(submitLine)
cd(ctrl.workDir)

% import results
A = importFodi([simLoc filesep simName '.fodi']);

A = [-1000 4*1e3 1].*A;   % Factor 4 to expand double symmetry

