function createDatFile(simLoc,simName,simInputs)
%createDatFile(simLoc,simName,simInputs) generates a .dat file at the specified
% location by merging the pre-defined APDL header, the inputs supplied by the user
% and the pre-defined APDL body.
%
% INPUTS :  simLoc    - String specifiying the absolute or relative path to the 
%                       location where the .dat file should be placed.
%
%           simName   - String with the name of the .dat file to be generated.
%           
%           simInputs - Struct containing the user inputs.
%           .friction - Friction coefficient between indenter and half space.
%           .nlgeom   - Nonlinear solution or not [0,1].
%           .FMax     - Applied max load, supplied as nano-Newton.
%           .iRadius  - Indenter radius, supplied in nano-meter.
%           .EL       - Longitudinal Young's modulus, supplied in MPa.
%           .ET       - Transverse Young's modulus, supplied in MPa.    
%           .nu       - Poisson's ratio, currently both of them.
%           .GLT      - Shear modulus, supplied in MPa.
%           .rAZ      - Rotation around Z from X toward Y
%           .rAX      - Rotation around X from Y toward Z
%           .rAY      - Rotation around Y from Z toward X
%           .visc     - Viscosity included or not [0,1]
%           .vFac     - Viscous coefficient in the prony series.
%
% OUTPUTS : NONE
%
% ABOUT : 
% The main advantage of using this method is that all inputs specified by the user
% are located in a single place, and if something goes wrong, the code does not run.
% Alternative approaches using STRFIND tend suffer from an inability to tell how
% many replacement operations were performed, and tend to not raise exceptions
% if for some reason the code did not change the variables.
%
% created by : August Brandberg augustbr at kth dot se
% date: 2020-11-14
%

% Prepare simulation : generate strings
fricCond = ['fricCond = ' num2str(simInputs.friction)];
nlCond   = ['nlCond = '   num2str(simInputs.nlgeom)];
FMax     = ['Fmax = '     num2str(simInputs.FMax/1e3)];
iRadius  = ['iRadius = '  num2str(simInputs.iRadius/1000)];
eLong    = ['eLong = '    num2str(simInputs.EL)];
eTran    = ['eTran = '    num2str(simInputs.ET)];
nu       = ['nu = '       num2str(simInputs.nu)];
GLT      = ['GLT = '      num2str(simInputs.GLT)];
rAZ      = ['rAZ = '      num2str(simInputs.rAZ)];
rAX      = ['rAX = '      num2str(simInputs.rAX)];
rAY      = ['rAY = '      num2str(simInputs.rAY)];
visc     = ['visco = '    num2str(simInputs.visc)];
vFac     = ['vFac = '     num2str(simInputs.vFac)];



% Prepare simulation : import header and body
part1SimFile = 'part1Sim.txt';
fileID = fopen(part1SimFile,'r');
part1Sim = fread(fileID,'*char');
fclose(fileID);
    
part2SimFile = 'part2Sim_v1.txt';
fileID = fopen(part2SimFile,'r');
part2Sim = fread(fileID,'*char');
fclose(fileID);



% Write the real simulation file
fileID = fopen(horzcat(simLoc,filesep,simName,'.dat'),'w');
fprintf(fileID,'%s\n',part1Sim,fricCond,nlCond,FMax,iRadius,eLong,eTran,nu,GLT,rAZ,rAX,rAY,visc,vFac,part2Sim);
fclose(fileID);

end % End of function createDatFile.m