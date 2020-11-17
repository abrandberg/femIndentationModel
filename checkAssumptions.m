function checkAssumptions(ctrl,hyperParameters,simLoc,simName,submitLine,savePath)


hyperParameters.unloadingFitFunction = 'Feng';%'Oliver-Pharr';



% Preparare simulation : Switches
simInputsA(1).friction = 0 ;                                        % Frictional coefficient.
simInputsA(1).nlgeom = 1 ;                                          % Account for large deformations.
simInputsA(1).FMax = 20000/4;                                       % Peak force, nano-Newton. Note, only 1/4th of real value!
simInputsA(1).iRadius = 300;                                        % Indenter radius, nano-meter.
simInputsA(1).EL = 11.2e3;                                          % Elastic longitudinal modulus, MPa.
simInputsA(1).ET = 1.85e3;                                          % Elastic transverse modulus, MPa.
simInputsA(1).GLT = 2.54e3;                                         % Elastic transverse modulus, MPa.
simInputsA(1).nu = 0.25*simInputsA(1).ET/simInputsA(1).EL;          
simInputsA(1).rAZ = 90;                                             
simInputsA(1).rAX = 0;                                              
simInputsA(1).rAY = 0;                                              
simInputsA(1).visc = 1;
simInputsA(1).vFac = 1.4;

simInputsB(1).friction = 0.0;
simInputsB(1).nlgeom = 1 ;
simInputsB(1).FMax = 20000/4;
simInputsB(1).iRadius = 300;
simInputsB(1).EL = 11.2e3;
simInputsB(1).ET = 1.85e3;
simInputsB(1).GLT = 2.54e3;
simInputsB(1).nu = 0.25*simInputsB(1).ET/simInputsB(1).EL;
simInputsB(1).rAZ = 0;
simInputsB(1).rAX = 0;
simInputsB(1).rAY = 0;
simInputsB(1).visc = 1;
simInputsB(1).vFac = 1.4;




% Preparare simulation : NLGEOM
simInputsA(2).friction = 0;
simInputsA(2).nlgeom = 0;
simInputsA(2).FMax = 20000/4;
simInputsA(2).iRadius = 300;
simInputsA(2).EL = 11.2e3;
simInputsA(2).ET = 1.85e3;
simInputsA(2).GLT = 2.54e3;
simInputsA(2).nu = 0.25*simInputsA(2).ET/simInputsA(2).EL;
simInputsA(2).rAZ = 90;
simInputsA(2).rAX = 0;
simInputsA(2).rAY = 0;
simInputsA(2).visc = 1;
simInputsA(2).vFac = 1.4;

simInputsB(2).friction = 0.0;
simInputsB(2).nlgeom = 0;
simInputsB(2).FMax = 20000/4;
simInputsB(2).iRadius = 300;
simInputsB(2).EL = 11.2e3;
simInputsB(2).ET = 1.85e3;
simInputsB(2).GLT = 2.54e3;
simInputsB(2).nu = 0.25*simInputsB(2).ET/simInputsB(2).EL;
simInputsB(2).rAZ = 0;
simInputsB(2).rAX = 0;
simInputsB(2).rAY = 0;
simInputsB(2).visc = 1;
simInputsB(2).vFac = 1.4;



% Preparare simulation : FRICTION
simInputsA(3).friction = 0.5;
simInputsA(3).nlgeom = 1;
simInputsA(3).FMax = 20000/4;
simInputsA(3).iRadius = 300;
simInputsA(3).EL = 11.2e3;
simInputsA(3).ET = 1.85e3;
simInputsA(3).GLT = 2.54e3;
simInputsA(3).nu = 0.25*simInputsA(2).ET/simInputsA(2).EL;
simInputsA(3).rAZ = 90;
simInputsA(3).rAX = 0;
simInputsA(3).rAY = 0;
simInputsA(3).visc = 1;
simInputsA(3).vFac = 1.4;

simInputsB(3).friction = 0.5;
simInputsB(3).nlgeom = 1;
simInputsB(3).FMax = 20000/4;
simInputsB(3).iRadius = 300;
simInputsB(3).EL = 11.2e3;
simInputsB(3).ET = 1.85e3;
simInputsB(3).GLT = 2.54e3;
simInputsB(3).nu = 0.25*simInputsB(2).ET/simInputsB(2).EL;
simInputsB(3).rAZ = 0;
simInputsB(3).rAX = 0;
simInputsB(3).rAY = 0;
simInputsB(3).visc = 1;
simInputsB(3).vFac = 1.4;




    indentationSet(1).relativeHumidity = 0;
    indentationSet(1).indenterType = 'hemispherical';
    indentationSet(1).indentationNormal = 'L';
    indentationSet(1).springConstant = nan;
    indentationSet(1).areaFile = nan;
    indentationSet(1).targetDir = nan;
    






figure;
subplot(1,4,3)
plot(8.45,2.92,'r+')
hold on

subplot(1,4,4)
plot(59,79,'r+')
hold on

x0 = 1;

for cLoop = 1:numel(simInputsA)
    
    tic
    createDatFile(simLoc,simName,simInputsA(cLoop))
    A = submitSimulation(simLoc,simName,submitLine,ctrl);
    A = A(1:end-1,:);
    F = @(x,xdata)1e9*(4*sqrt(simInputsA(cLoop).iRadius/1e9)*simInputsA(cLoop).EL*1e6/(3)/(1-(simInputsA(cLoop).nu)^2)*(xdata*1e-9).^x(1));
    
    [x,~,~,~,~] = lsqcurvefit(F,x0,A(:,1),A(:,2));

    % Visualization
    subplot(1,4,1)
    plot(A(:,1),A(:,2),'o-','DisplayName','FEM')
    hold on
    plot(A(:,1),F(x,A(:,1)),'k--','DisplayName',['Fit $F(z) = 4E \sqrt{R}/(3(1-\nu^2))' ' \cdot z^{' num2str(round(x(1),3)) '}$'])
    legend('location','best','interpreter','latex')
    xlabel('$z$ [nm]','interpreter','latex')
    ylabel('$F$ [nN]','interpreter','latex')
    pause(0.5)
        
    % Up-sample data
    resultFile = interp1(A(:,3),A,linspace(A(1,3),A(end,3),2000*A(end,3)));
    indentationSet(cLoop).designatedName = num2str(cLoop);
    [Er,~,~,~,~] = modulusFitter(indentationSet,ctrl,hyperParameters,resultFile(:,[1 2]))
    
    
    ErSaveA(cLoop) = Er;
    
    outputSave(cLoop).fodiA = A;
    
% B sim
    createDatFile(simLoc,simName,simInputsB(cLoop))
    B = submitSimulation(simLoc,simName,submitLine,ctrl);
    B = B(1:end-1,:);
    F = @(x,xdata)1e9*(4*sqrt(simInputsB(cLoop).iRadius/1e9)*simInputsB(cLoop).EL*1e6/(3)/(1-(simInputsB(cLoop).nu)^2)*(xdata*1e-9).^x(1));
    [x,~,~,~,~] = lsqcurvefit(F,x0,B(:,1),B(:,2));


    % Visualization
    subplot(1,4,2)
    plot(B(:,1),B(:,2),'o-','DisplayName','FEM')
    hold on
    plot(B(:,1),F(x,B(:,1)),'k--','DisplayName',['Fit $F(z) = 4E \sqrt{R}/(3(1-\nu^2))' ' \cdot z^{' num2str(round(x(1),3)) '}$'])
    legend('location','best','interpreter','latex')
    xlabel('$z$ [nm]','interpreter','latex')
    ylabel('$F$ [nN]','interpreter','latex')
    pause(0.5)
        
    % Up-sample data
    resultFile = interp1(B(:,3),B,linspace(B(1,3),B(end,3),2000*B(end,3)));
    indentationSet(cLoop).designatedName = num2str(cLoop);
    [Er,~,~,~,~] = modulusFitter(indentationSet,ctrl,hyperParameters,resultFile(:,[1 2]))
    
    
    ErSaveB(cLoop) = Er;
    outputSave(cLoop).fodiB = B;
    
    % plot the indentation modulus
    subplot(1,4,3)
    plot(ErSaveA(cLoop),ErSaveB(cLoop),'sk')
    hold on
    xlabel('E_L'); ylabel('E_T')
    
    
    
    % Plot the indentation depth
    subplot(1,4,4)
    plot(max(A(:,1)),max(B(:,1)),'or')
    hold on
    xlabel('x0_L'); ylabel('x0_T')
    
    pause(0.5)
    toc
end



save(savePath,'-v7.3')
