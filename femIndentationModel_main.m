%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% femIndentationModel_main
%
% Generates and runs ANSYS models.
%
% See README.md for documentation.
%
% created by : August Brandberg augustbr at kth dot se
% date : 2020-11-14
%
clear; close all; clc
format compact

% Preliminary instructions
ctrl.workDir = cd;
ctrl.verbose = 0;
ctrl.formatSpec = '%20s %40s %40s \n';
addpath([ctrl.workDir filesep 'auxilliaryFunctions'])

% File and folder names and paths
simLoc = [ctrl.workDir filesep 'runFolder']
simName = 'file'


% Specify ANSYS call according to your installation
ansysExecPath = '"C:\Program Files\ANSYS Inc\v202\ansys\bin\winx64\MAPDL.exe"';
ansysLic      = '-p aa_t_a';
ansysNP       = '-np 8';
ansysMisc     = '-dis -mpi INTELMPI -s read -l en-us -b';
ansysJobname  = ['-j ' simName];
ansysRunDir   = ['-dir ' simLoc];
ansysInput    = ['-i ' simLoc filesep simName '.dat'];
ansysOutput   = ['-o ' simLoc filesep simName '.out'];

submitLine = strjoin({ansysExecPath , ansysLic , ansysNP , ansysMisc , ansysJobname , ansysRunDir , ansysInput , ansysOutput},' ');




% The function modulusFitter.m fits the FEM output to determine the
% indentation modulus. These are inputs, created to match the conditions we
% compare against. See modulusFitter.m for complete documentation.
hyperParameters.epsilon              = 0.75;
hyperParameters.sampleRate           = 2000;
hyperParameters.thermpnt             = 2000;
hyperParameters.unloadingFitFunction = 'Oliver-Pharr';
hyperParameters.constrainHead        = 0;
hyperParameters.constrainTail        = 0;
hyperParameters.unloadingFitRange    = 1900;
hyperParameters.compensateCreep      = 1;



% Generate comparison between FEM and Jäger et al.
elPlotJager(ctrl,hyperParameters,simLoc,simName,submitLine,'ELParaResume.mat')

% Generate comparison between FEM and Jäger et al.
etPlotJager(ctrl,hyperParameters,simLoc,simName,submitLine,'ETParaResume.mat')


% Generate comparison between FEM and Jäger et al.
gltPlotJager(ctrl,hyperParameters,simLoc,simName,submitLine,'GLTParaResume.mat')



% Check elasticity assumption
calibrateVisco(ctrl,hyperParameters,simLoc,simName,submitLine,'ViscoCalibration.mat')

% Check model assumptions
checkAssumptions(ctrl,hyperParameters,simLoc,simName,submitLine,'CheckAssumptions.mat')




