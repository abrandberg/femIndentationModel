function etPlotJager(ctrl,hyperParameters,simLoc,simName,submitLine,savePath)

x0 = 1;

El = 55*1e3;
Et = 10*1e3;
Glt = 3*1e3;
nutl = 0.25;
nutt = 0.25;
ETTry = 1e3*[5:15 10];

simInputs(1).friction = 0 ;                                                         % Frictional coefficient.
simInputs(1).nlgeom = 0 ;                                                           % Account for large deformations.
simInputs(1).FMax = 20000/4;                                                        % Peak force, nano-Newton. Note, only 1/4th of real value!
simInputs(1).iRadius = 300;                                                         % Indenter radius, nano-meter.
simInputs(1).EL = El;                                                               % Elastic longitudinal modulus, MPa.
simInputs(1).ET = ETTry(1);                                                         % Elastic transverse modulus, MPa.
simInputs(1).GLT = Glt;                                                             %
simInputs(1).nu = nutl*simInputs(1).ET/simInputs(1).EL;                             %
simInputs(1).rAZ = 90;                                                              %
simInputs(1).rAX = 0;                                                               %
simInputs(1).rAY = 0;                                                               %
simInputs(1).visc = 0;                                                              %
simInputs(1).vFac = 0; 

% Copy the settings to every simulation in the sweep, add appropriate ET value.
for bLoop = 2:length(ETTry)                                                     
   simInputs(bLoop) = simInputs(1);
   simInputs(bLoop).ET = ETTry(bLoop);
   simInputs(bLoop).nu = nutl*simInputs(bLoop).ET/simInputs(bLoop).EL;
end


[~,~,ETJager0,ETJager90,~,~,~,~] = importJagerFigure4();
figure;
subplot(1,2,2)
plot(ETJager0(:,1),ETJager0(:,2),'bs-','DisplayName','J0')
hold on
plot(ETJager90(:,1),ETJager90(:,2),'bs-','DisplayName','J90')
clear ErSave
for cLoop = 1:length(ETTry)
    tic
    createDatFile(simLoc,simName,simInputs(cLoop))
    A = submitSimulation(simLoc,simName,submitLine,ctrl);
    A = A(1:end-1,:);
    F = @(x,xdata)1e9*(4*sqrt(simInputs(cLoop).iRadius/1e9)*simInputs(cLoop).EL*1e6/(3)/(1-(simInputs(cLoop).nu)^2)*(xdata*1e-9).^x(1));

    [x,~,~,~,~] = lsqcurvefit(F,x0,A(:,1),A(:,2));

    % Visualization
    subplot(1,2,1)
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
    
    
    simOutputs(cLoop).Er = Er;
    ErSave(cLoop) = Er;
    subplot(1,2,2)
    plot(simInputs(cLoop).ET/1000,Er,'sk')
    pause(0.5)
    toc
end

% Change direction from along EL direction to ET direction
for dLoop = 1:length(ETTry)   
   simInputs(dLoop).rAZ = 0;
end


for eLoop = 1:length(ETTry)
    tic
    createDatFile(simLoc,simName,simInputs(eLoop))
    A = submitSimulation(simLoc,simName,submitLine,ctrl);
    A = A(1:end-1,:);
    F = @(x,xdata)1e9*(4*sqrt(simInputs(eLoop).iRadius/1e9)*simInputs(eLoop).EL*1e6/(3)/(1-(simInputs(eLoop).nu)^2)*(xdata*1e-9).^x(1));
    
    [x,~,~,~,~] = lsqcurvefit(F,x0,A(:,1),A(:,2));

    % Visualization
    subplot(1,2,1)
    plot(A(:,1),A(:,2),'o-','DisplayName','FEM')
    hold on
    plot(A(:,1),F(x,A(:,1)),'k--','DisplayName',['Fit $F(z) = 4E \sqrt{R}/(3(1-\nu^2))' ' \cdot z^{' num2str(round(x(1),3)) '}$'])
    legend('location','best','interpreter','latex')
    xlabel('$z$ [nm]','interpreter','latex')
    ylabel('$F$ [nN]','interpreter','latex')
    pause(0.5)
    
    
    % Up-sample data
    resultFile = interp1(A(:,3),A,linspace(A(1,3),A(end,3),2000*A(end,3)));
    indentationSet(eLoop).designatedName = num2str(eLoop);
    [Er,~,~,~,~] = modulusFitter(indentationSet,ctrl,hyperParameters,resultFile(:,[1 2]))
    
    
    simOutputs(eLoop).Er = Er;
    ErSave2(eLoop) = Er;
    subplot(1,2,2)
    plot(simInputs(eLoop).ET/1000,Er,'sk')
    pause(0.5)
    toc
end

save(savePath,'-v7.3')






