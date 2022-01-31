function TAP_combinedModel

% (c) 2020 Adrián Lahuerta Lavieja and Martin Johansson;
% e-mail: adrian.lahuerta@kuleuven.be and martin.n.johansson@ericsson.com

% The code (or parts of it) may be used for non-profit purposes as long
% the copyright notice is included and [A] is credited and cited.

% [A] A. Lahuerta-Lavieja, M. Johansson, C. Larsson, U. Gustavsson, and 
% G. A. E. Vandenbosch, "Computationally efficient millimeter-wave 
% backscattering models: A combined blockage and backscattering model,"
% submitted to IEEE Antennas Wireless Propag. Lett., 2022.

%% Scenario for testing blockage and backscattering
scenario.number = 3; % Scenario from Fig. 3 in [A]

%% Print scenario
fprintf(['\n Scenario ' num2str(scenario.number) ' \n \n'])

%% Implemented model: string
models.List = {'3DFresnel'; 'erf'; 'METIS'; 'mmMAGIC'; 'ITUFresnel'};
% models.List = {'3DFresnel'};

%% PO data
PO.available = 1;

%% MoM data
MoM.available = 1;

%% Main code
addpath('Fresnel_Integrals')
scenario                     = getScenarioInfo(scenario);
[txSurf0, rxSurf0, scenario] = getAntennas(scenario);
[scattSurf0, scenario]       = getScatteringSurfaces(scenario);
scenario                     = getSweepVariable(scenario);
PO                           = getPOresults(scenario, PO);
MoM                          = getMoMresults(scenario, MoM);
[scenario, models]           = calculateModels(scenario, models, scattSurf0, txSurf0, rxSurf0);
getErrorMetric(MoM, models);
plotResults(scenario, PO, MoM, models);
end

%% Frequency
function scenario = getScenarioInfo(scenario)
scenario.f      = 28e9;
scenario.c0     = 299792458;
scenario.z0     = 377;
scenario.lambda = scenario.c0/scenario.f;
scenario.k0     = 2*pi/scenario.lambda;
end

%% Antennas
function [txSurf0, rxSurf0, scenario] = getAntennas(scenario)
switch scenario.number
    case 3
        yAntennaWallDistance = 0;
        xAntennaSeparation   = 0;
        zAntHeight           = 0;
        TxPhiAntToWall = pi/2;                         % Point along positive y-axis
        RxPhiAntToWall = TxPhiAntToWall;
        
        txSurf0 = isotropic_TAP();
        rxSurf0 = nardaV637_TAP(scenario.f);
        
        txSurf0.geomT{end+1} = [-xAntennaSeparation/2 -yAntennaWallDistance zAntHeight]; % Move x meters along negative y-axis
        txSurf0.geomT{end+1} = [0 0 1 TxPhiAntToWall]; % Phi rotation
        txSurf0.geomT{end+1} = [-1 0 0];               % Move 1 meter along negative y-axis
        
        rxSurf0.geomT{end+1} = [xAntennaSeparation/2 -yAntennaWallDistance zAntHeight]; % Move x meters along negative y-axis
        rxSurf0.geomT{end+1} = [0 0 1 RxPhiAntToWall]; % Phi rotation
        
    otherwise
        error('Scenario not defined!')
end

[~, xyzTx0] = reduce_geomT(txSurf0.geomT);
[~, xyzRx0] = reduce_geomT(rxSurf0.geomT);

scenario.xyzTx0 = xyzTx0;
scenario.xyzRx0 = xyzRx0;

end

%% Scattering surfaces
function [scattSurf0, scenario] = getScatteringSurfaces(scenario)
nWidth    = 2;
nHeight   = 2;
nSamples  = [nWidth nHeight];

switch scenario.number
    case 3
        scattSurf0 = smallRectangleForBlockage(nSamples, scenario.f);
end

%% Gamma
switch scenario.number
    case 3
        gamma = -1; % PEC
end

scenario.gamma = gamma;
end

%% Sweep variable
function scenario = getSweepVariable(scenario)

resolution = getResolution(scenario);
interval   = getInterval(scenario);

nSamples  = abs((interval(2) - interval(1)))/resolution + 1;
sweepVals = linspace(interval(1), interval(2), nSamples);

idxPO = 1:numel(sweepVals);

sweepVals = sweepVals*pi/180;

scenario.sweepVals = sweepVals;
scenario.idxPO     = idxPO;

end

function resolution = getResolution(scenario)
switch scenario.number
    case 3 % Rx antenna position in the circle
        resolution = 0.1;  % [degree]
    otherwise
        resolution = [];
end
end

function interval = getInterval(scenario)
switch scenario.number
    case 3 % Rx antenna position in the circle
        interval = [-90 270];                  % [degree]
    otherwise
        interval = [];
end
end

%% Models calculation
function [scenario, models] = calculateModels(scenario, models, scattSurf0, txSurf0, rxSurf0)

% Variable preallocation for speed
models.resultsTh = zeros(numel(scenario.sweepVals), 1);
models.resultsPh = zeros(numel(scenario.sweepVals), 1);
gammaFresnel = zeros(numel(scenario.sweepVals),1);
[gammaFresnelSingleReflection] = deal(gammaFresnel);

for k = 1:numel(scenario.sweepVals)  % Move/rotate of surface(s) or antenna (s);
    
    % Move/rotate surface(s)
    switch scenario.number
        case 3
            scattSurf{1} = scattSurf0{1};
            R = 1; % 1 meter radius circle
            x_position        = R*cos(scenario.sweepVals(k));
            y_position        = R*sin(scenario.sweepVals(k));       % y0 = -1 is the starting position
            rxSurf0.geomT{4}  = [x_position y_position 0];
            scenario.xyzRx0   = rxSurf0.geomT{4};
            rxAntennaPointing = atan2(y_position, x_position) - pi; % Point to the center of coordinates
            rxSurf0.geomT{5}  = [0 0 1 rxAntennaPointing];
    end
    
    scattSurfOrdered = scattSurf;
    
    % Coordinates for surfaces needed for Fresnel integrals
    for z = 1:numel(scattSurfOrdered)
        iteration.xyz0(z,:) = scattSurfOrdered{z}.centroid;         % Center of surface
        iteration.xyz1(z,:) = scattSurfOrdered{z}.vertices(1,:);    % Corner 1 of surface
        iteration.t1(z,:)   = scattSurfOrdered{z}.vertices(2,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.t2(z,:)   = scattSurfOrdered{z}.vertices(3,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.n(z,:)    = scattSurfOrdered{z}.normalVectors(1,:); % Surface normal
        iteration.area(z,:) = scattSurfOrdered{z}.area;
    end
    
    
    iteration.xyzTx(1,:) = scenario.xyzTx0;
    iteration.xyzRx(1,:) = scenario.xyzRx0;
    [iteration.xyzRP(1,:), iteration.xyzTxImage(1,:), scenario.blockedPath(1)] = getSpecularReflectionPoint(scenario.xyzTx0, scenario.xyzRx0, iteration.xyz1(1,:), iteration.n(1,:));
    iteration.xyzRxImage(1,:) = getImageLocation(scenario.xyzRx0, iteration.xyz1(1,:), iteration.n(1,:));
    
    % Get distance through specular reflection point
    R0reflex        = norm(scenario.xyzTx0-iteration.xyzRP(1,:)) + norm(scenario.xyzRx0-iteration.xyzRP(1,:));
    scenario.R(k,1) = norm(scenario.xyzRx0-iteration.xyzTxImage(1,:));
    
    % Find point used for free space path loss and antenna gain calculation
    % Point on surface closest to specular reflection point
    [iteration.xyzAmp,s,t] = getPointOnSurfaceClosestToRP(iteration.xyzRP(1,:), iteration.xyz1(1,:), iteration.t1(1,:), iteration.t2(1,:));
    
    % Find end points for Fresnel integrals
    % Project the reflection point on (the extensions of) the four
    % rectangle edges and use these for the calculation "excess distance".
    iteration.xyzEdges          = getProjectionsOnEdges(iteration.xyzRP(1,:), iteration.xyz1(1,:), iteration.t1(1,:), iteration.t2(1,:));
    iteration.xyzEdgesMidpoints = getEdgesMidpoint(scattSurf, iteration.t1(1,:), iteration.t2(1,:));
    
    % Get the gain toward the single ray scattering point
    [iteration.eThTx, iteration.ePhTx, iteration.eThRx, iteration.ePhRx] = ...
        getAntennasTowardSingleRayScatteringPoint(iteration.xyzAmp, scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
    
    % LoS component
    scenario.rLOS(k) = max(0,norm(scenario.xyzRx0-scenario.xyzTx0));
    gffLOS(k)        = exp(-1j*scenario.k0*scenario.rLOS(k))/scenario.rLOS(k) * scenario.lambda/(4*pi); %#ok<*AGROW>
    [eThTx(k), ePhTx(k), eThRx(k), ePhRx(k)] = getAntennaTowardAntenna(scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
    
    % Signs for summing of integrals
    iteration.signs = 1 - 2*[s>0 s<1 t>0 t<1];
    
    % Excess path lengths for Fresnel integrals
    iteration.deltaR(:,1)  = max(0,sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzTx0),2)) + sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzRx0),2)) - R0reflex);
    
    % The phase reference is the specular reflection point.
    rMiddle = norm(iteration.xyzRP(1,:)-scenario.xyzTx0) + norm(iteration.xyzRP(1,:)-scenario.xyzRx0);
    scenario.rAmp(k,1) = norm(iteration.xyzAmp-scenario.xyzTx0) + norm(iteration.xyzAmp-scenario.xyzRx0);
    
    % Green's function times "Friis factor"
    iteration.gff    = exp(-1j*scenario.k0*rMiddle)/rMiddle * scenario.lambda/(4*pi);
    
    % Calculate the different models
    if any(strcmp(models.List, '3DFresnel'))
        gammaFresnelSingleReflection(k,1) = ourFresnelModel(scenario, iteration, 1);
    end
    if any(strcmp(models.List, 'erf'))
        gammaErfSingleReflection(k,1)     = ourErfModel(scenario, iteration, 1);
    end
    if any(strcmp(models.List, 'METIS'))
        gammaMetisSingleReflection(k,1)   = metisModel(scenario, iteration, 1);
    end
    if any(strcmp(models.List, 'METISRCS'))
        gammaMetisRCSSingleReflection(k,1)= metisRCSmodel(scenario, iteration, 1);
    end
    if any(strcmp(models.List, 'mmMAGIC'))
        gammaMMMagicSingleReflection(k,1) = mmMAGICModel(scenario, iteration, k, 1);
    end
    if any(strcmp(models.List, 'ITUFresnel'))
        gammaITUFresSingleReflection(k,1) = ITUfresnelModel(scenario, iteration, 1);
    end
    if any(strcmp(models.List, 'ITUEmpirical'))
        gammaITUEmpSingleReflection(k,1)  = ITUsemiEmpModel(scenario, iteration, k, 1);
    end
    
    if scenario.blockedPath
        if any(strcmp(models.List, '3DFresnel'))
            gammaFresnel(k,1) = -gammaFresnelSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'erf'))
            gammaErf(k,1) = -gammaErfSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'METIS'))
            gammaMetis(k,1) = -gammaMetisSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'METISRCS'))
            gammaMetisRCS(k,1) = -gammaMetisRCSSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'mmMAGIC'))
            gammaMMMagic(k,1) = -gammaMMMagicSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'ITUFresnel'))
            gammaITUFres(k,1) = -gammaITUFresSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'ITUEmpirical'))
            gammaITUEmp(k,1) = -gammaITUEmpSingleReflection(k,1);
        end
    else
        if any(strcmp(models.List, '3DFresnel'))
            gammaFresnel(k,1) = scenario.gamma(1)*gammaFresnelSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'erf'))
            gammaErf(k,1) = scenario.gamma(1)*gammaErfSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'METIS'))
            gammaMetis(k,1) = scenario.gamma(1)*gammaMetisSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'METISRCS'))
            gammaMetisRCS(k,1) = scenario.gamma(1)*gammaMetisRCSSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'mmMAGIC'))
            gammaMMMagic(k,1) = scenario.gamma(1)*gammaMMMagicSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'ITUFresnel'))
            gammaITUFres(k,1) = scenario.gamma(1)*gammaITUFresSingleReflection(k,1);
        end
        if any(strcmp(models.List, 'ITUEmpirical'))
            gammaITUEmp(k,1) = scenario.gamma(1)*gammaITUEmpSingleReflection(k,1);
        end
    end
    
    counter = 1; % Variable to count the number of methods that have been calculated at the moment
    if any(strcmp(models.List, '3DFresnel'))
        modelsTmp(counter, k) = gammaFresnel(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'erf'))
        modelsTmp(counter, k) = gammaErf(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'METIS'))
        modelsTmp(counter, k) = gammaMetis(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'mmMAGIC'))
        modelsTmp(counter, k) = gammaMMMagic(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'ITUFresnel'))
        modelsTmp(counter, k) = gammaITUFres(k, 1);
    end
    modelsTh(:,k) = modelsTmp(:,k)*iteration.eThTx*iteration.eThRx*iteration.gff;
    modelsPh(:,k) = modelsTmp(:,k)*iteration.ePhTx*iteration.ePhRx*iteration.gff;
end

% LoS component
LOScomponentTh = eThTx.*eThRx.*gffLOS;
LOScomponentPh = ePhTx.*ePhRx.*gffLOS;

models.resultsTh = modelsTh + LOScomponentTh;
models.resultsPh = modelsPh + LOScomponentPh;
end

function PO = getPOresults(scenario, PO)
tempPO    = loadPOdata(scenario.number);
POresults = tempPO(:);

PO.results = POresults;
end

function MoM = getMoMresults(scenario, MoM)

if MoM.available == 1
    tempMoM = loadMoMdata(scenario.number);
    MoMresults = tempMoM(:);
else
    MoMresults = [];
end

MoM.results = MoMresults;
end

function [scenario, PO, MoM, models] = plotResults(scenario, PO, MoM, models)


scenario.sweepVals = scenario.sweepVals/pi*180;

% Amplitude plot
[models.tikzAmp, PO.tikzAmp, MoM.tikzAmp]       = plotAmplitudResults(scenario, models, PO, MoM);
% Phase plot
[models.tikzPhase, PO.tikzPhase, MoM.tikzPhase] = plotPhaseResults(scenario, models, PO, MoM);

if PO.available && not(any(strcmp(models.List, 'PO'))) % Add PO to the legend before display
    models.List = [models.List; 'PO'];
end
if MoM.available && not(any(strcmp(models.List, 'MoM'))) % Add MoM to the legend before display
    models.List = [models.List; 'MoM'];
end

models.legend = models.List;


figure(999);  legend(models.List)
figure(1000); legend(models.List)
end

function [modelsTikzAmp, POtikzAmp, MoMtikzAmp] = plotAmplitudResults(scenario, models, PO, MoM)

figure(999);
clf

results = models.resultsTh.';

xValue = scenario.sweepVals;
modelsTikzAmp = zeros(size(results));
index = 1:numel(xValue);

for iii = 1:numel(models.List)
    modelsTikzAmp(:, iii) = 10*log10(abs2(results(index,iii)));
end

if PO.available
    POtikzAmp = 10*log10(abs2(PO.results(index,1)));
else
    POtikzAmp = [];
end

if MoM.available
    MoMtikzAmp = 10*log10(abs2(MoM.results(index,1)));
else
    MoMtikzAmp = [];
end

plotScenario3Amplitude(scenario, xValue, modelsTikzAmp, POtikzAmp, MoMtikzAmp);
end

function [modelsTikzPhase, POtikzPhase, MoMtikzPhase] = plotPhaseResults(scenario, models, PO, MoM)

figure(1000);
clf

results = models.resultsTh.';

xValue = scenario.sweepVals;
modelsTikzPhase = zeros(size(results));

index = 1:numel(xValue);

for iii = 1:numel(models.List)
    
            phase_offset = MoM.results(index,1);
    
    modelsTikzPhase(:, iii) = unwrap(angle(results(index,iii)./phase_offset(index,1)))*180/pi;
    
end

if PO.available
    POtikzPhase = unwrap(angle(PO.results(index,1)./phase_offset(index,1)))*180/pi;
else
    POtikzPhase = [];
end

if MoM.available
    MoMtikzPhase = unwrap(angle(MoM.results(index,1)./phase_offset(index,1)))*180/pi;
else
    MoMtikzPhase = [];
end

plotScenario3Phase(scenario, xValue, modelsTikzPhase, POtikzPhase, MoMtikzPhase);
end

function plotScenario3Amplitude(scenario, x, models, PO, MoM)
subplot(1,3,1)
index = find(x >= -80 & x <= -70);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([-80 -70])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on

subplot(1,3,2)
index = find(x >= -5 & x <= 5);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([-5 5])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on

subplot(1,3,3)
index = find(x >= 80 & x <= 90);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([80 90])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on
end

function plotScenario3Phase(scenario, x, models, PO, MoM)
subplot(1,3,1)
index = find(x >= -80 & x <= -70);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([-80 -70])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on

subplot(1,3,2)
index = find(x >= -5 & x <= 5);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([-5 5])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on

subplot(1,3,3)
index = find(x >= 80 & x <= 90);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([80 90])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on
end

function getProperXaxisLabel(scenario)
switch scenario.number
    case 3
        xlabel('\alpha (deg)')
    otherwise
        xlabel('Undefined label. Please update')
end
end

function getProperAxisLimitsAmplitude(scenario)
switch scenario.number
    case 3
        ylim([-90 -45]);
    otherwise
end
end

function getProperAxisLimitsPhase(scenario)
switch scenario.number
    case 3
        xlim([-90 90]);
        ylim([-90 90]);
    otherwise
end
end

function [POresults] = loadPOdata(scenario)
load(['SimulationData\POscenario' num2str(scenario)], 'POresults');
end

function [MoMresults] = loadMoMdata(scenario)
load(['SimulationData\MoMscenario' num2str(scenario)], 'MoMresults');
end

function xyzTxImage  = getImageLocation(xyzTx,xyzOnSurface,n)
% getImageLocation  Image source in plane containing
% 'xyzOnSurface' with normal 'n'.
xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
end

function [xyzRP, xyzTxImage, blockedPath]  = getSpecularReflectionPoint(xyzTx, xyzRx, xyzOnSurface, n)
% getSpecularReflectionPoint  Specular reflection point on plane containing
% 'xyzOnSurface' with normal 'n'.

blockedPath = false;

% First check whether the Tx and Rx are in different sides of the plane
% defined by the scattering surface

xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
u         = (xyzRx-xyzTxImage)/norm(xyzRx-xyzTxImage);
w         = xyzTxImage-xyzOnSurface;
s         = -(n*w.')/(n*u.');
xyzRP     = xyzTxImage + s*u;

% Check if the segment u1 intersects with the scattering surface. If it
% does, it's a reflection situation; otherwise, a blockage situation
u1         = (xyzRx-xyzTxImage);
s1         = -(n*w.')/(n*u1.');

if (s1 < 0 || s1 > 1) && not(any(ismembertol([0 1], s1))) && not(isinf(s1))
    
    blockedPath = true;
    
    u         = (xyzRx-xyzTx)/norm(xyzRx-xyzTx);
    w         = xyzTx-xyzOnSurface;
    s         = -(n*w.')/(n*u.');
    xyzRP     = xyzTx + s*u;
    xyzTxImage = xyzTx;
    
    return
end
end

function [xyz0,s0,t0] = getPointOnSurfaceClosestToRP(xyzRp,xyz1moved,t1,t2)
u  = t1/norm(t1); % Unit vector along first edge from first to second corner
v  = t2/norm(t2); % Unit vector along second edge from first to third corner
s0 = (xyzRp-xyz1moved)*u.'/norm(t1); % Projection of vector from rectangle center to reflection point on 'u'
t0 = (xyzRp-xyz1moved)*v.'/norm(t2); % Projection of vector from rectangle center to reflection point on 'v'
s  = s0;
if s<0
    s = 0; % Move to first edge; RP outside rectangle
elseif s>1
    s = 1; % Move to edge opposite of first edge; RP outside rectangle
end
t = t0;
if t<0
    t = 0; % Move to second edge
elseif t>1
    t = 1; % Move to edge opposite of second edge; RP outside rectangle
end
xyz0 = xyz1moved + s*t1 + t*t2;
end

function xyzEdges = getProjectionsOnEdges(xyzRp,xyz1moved,t1,t2)
% Points on line including xyzRp and parallel to to t1 and t2
u = t1/norm(t1); % Unit vector along first edge from first to second corner
v = t2/norm(t2); % Unit vector along second edge from first to third corner

% Project the vector xyz1moved-xyzRp on 'u' to get 'pu', and then
% compute xyzRp + pu and xyzRp + pu + t1 which are the points used for
% the Fresnel integral in the first "dimension" (of the surface). And
% correspondingly for the second dimension.
pu = ((xyz1moved-xyzRp)*u(:))*u; % Projection of xyz1moved-xyzRp on u
pv = ((xyz1moved-xyzRp)*v(:))*v; % Projection of xyz1moved-xyzRp on v
xyzEdges = xyzRp + [pu;pu+t1;pv;pv+t2];
end

function xyzEdgesMidpoints = getEdgesMidpoint(scattSurf,t1,t2)
xyzEdgesMidpoints(1, :) = scattSurf{1}.vertices(1,:) + t2/2;
xyzEdgesMidpoints(2, :) = scattSurf{1}.vertices(2,:) + t2/2;
xyzEdgesMidpoints(3, :) = scattSurf{1}.vertices(1,:) + t1/2;
xyzEdgesMidpoints(4, :) = scattSurf{1}.vertices(3,:) + t1/2;
end

function res = smallRectangleForBlockage(nSamples, freq)
% Test case for comparison with MoM

surfaceWidth  = 0.3;
surfaceHeight = 0.5;
surfaceOffset = [0 0 0];
surfScatt     = TAP_CreateRectangle(surfaceWidth,surfaceHeight,nSamples,freq);
surfScatt     = TAP_rotateObject(surfScatt, [1 0 0], pi/2); % Scattering currents in xz-plane
surfScatt     = TAP_moveObject(surfScatt,surfaceOffset);
res           = {surfScatt};

end

function [eThTx, ePhTx, eThRx, ePhRx] = getAntennasTowardSingleRayScatteringPoint(xyzAmp, xyzTx0, xyzRx0, txSurf0, rxSurf0)
% Antenna gains toward antenna gain reference point.

thTxMiddle = atan2(sqrt(sum(abs2(xyzAmp(:,1:2)-xyzTx0(1:2)),2)),xyzAmp(:,3)-xyzTx0(3));
phTxMiddle = atan2(xyzAmp(2)-xyzTx0(2),xyzAmp(1)-xyzTx0(1));
thRxMiddle = atan2(sqrt(sum(abs2(xyzAmp(:,1:2)-xyzRx0(1:2)),2)),xyzAmp(:,3)-xyzRx0(3));
phRxMiddle = atan2(xyzAmp(2)-xyzRx0(2),xyzAmp(1)-xyzRx0(1));

rHat_x  = sin(thTxMiddle)*cos(phTxMiddle);
rHat_y  = sin(thTxMiddle)*sin(phTxMiddle);
rHat_z  = cos(thTxMiddle);
R       = reduce_geomT(txSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;

[eThTx, ePhTx] = txSurf0.function_handle(theta, phi);

rHat_x  = sin(thRxMiddle)*cos(phRxMiddle);
rHat_y  = sin(thRxMiddle)*sin(phRxMiddle);
rHat_z  = cos(thRxMiddle);
R       = reduce_geomT(rxSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThRx, ePhRx] = rxSurf0.function_handle(theta, phi);

end

function [eThTx, ePhTx, eThRx, ePhRx] = getAntennaTowardAntenna(xyzTx0, xyzRx0, txSurf0, rxSurf0)

thTxToRx   = atan2(sqrt(sum(abs2(xyzRx0(:,1:2)-xyzTx0(1:2)),2)),xyzRx0(:,3)-xyzTx0(3));
phTxToRx   = atan2(xyzRx0(2)-xyzTx0(2),xyzRx0(1)-xyzTx0(1));

thRxToTx   = atan2(sqrt(sum(abs2(xyzTx0(:,1:2)-xyzRx0(1:2)),2)),xyzTx0(:,3)-xyzRx0(3));
phRxToTx   = atan2(xyzTx0(2)-xyzRx0(2),xyzTx0(1)-xyzRx0(1));

rHat_x  = sin(thTxToRx)*cos(phTxToRx);
rHat_y  = sin(thTxToRx)*sin(phTxToRx);
rHat_z  = cos(thTxToRx);
R       = reduce_geomT(txSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThTx, ePhTx] = txSurf0.function_handle(theta, phi);

rHat_x  = sin(thRxToTx)*cos(phRxToTx);
rHat_y  = sin(thRxToTx)*sin(phRxToTx);
rHat_z  = cos(thRxToTx);
R       = reduce_geomT(rxSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThRx, ePhRx] = rxSurf0.function_handle(theta, phi);
end

function gammaFresnel = ourFresnelModel(scenario, iteration, m)
% 3D Fresnel model
% [1] A. Lahuerta-Lavieja, M. Johansson, U. Gustavsson, T. A. H. Bressner, and G. A. E. Vandenbosch,
% "Computationally-efficient millimeter-wave back-scattering models," to be published in IEEE Trans. Antennas Propag., 2020.

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
fresnelInts = fresnelC(sqrt(4/lambda*deltaR)) - 1j*fresnelS(sqrt(4/lambda*deltaR));
gammaFresnel = 0.5*1j*(((fresnelInts(1)*signs(1)+fresnelInts(2)*signs(2))*(fresnelInts(3)*signs(3)+fresnelInts(4)*signs(4))));
end

function gammaErf = ourErfModel(scenario, iteration, m)
% erf model [1]
% [1] A. Lahuerta-Lavieja, M. Johansson, U. Gustavsson, T. A. H. Bressner, and G. A. E. Vandenbosch,
% "Computationally-efficient millimeter-wave back-scattering models," to be published in IEEE Trans. Antennas Propag., 2020.

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
erfInts = erf(sqrt(4/lambda*deltaR));
gammaErf = ((erfInts(1)*signs(1)+erfInts(2)*signs(2))*(erfInts(3)*signs(3)+erfInts(4)*signs(4)))*0.25;
end

function gammaMetis = metisModel(scenario, iteration, m)
% M-METIS model
% Based on METIS D1.4 Section C.1.4: Shadowing objects
% url: https://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;
gammaMetis = ((metisTerms(1)*signs(1)+metisTerms(2)*signs(2))*(metisTerms(3)*signs(3)+metisTerms(4)*signs(4)));
end

function [gammaMetisRCS] =  metisRCSmodel(scenario, iteration, m)
% METIS RCS scattering model
% Based on METIS D1.4 Section C.1.5 Scattering objects
% url: https://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
R1     = iteration.R1;
R2     = iteration.R2;
area   = iteration.area(:,m);

radius = sqrt(area/pi);
RCS    = pi*radius^2;
metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;
gammaMetisRCS0 = ((metisTerms(1)*signs(1)+metisTerms(2)*signs(2))*(metisTerms(3)*signs(3)+metisTerms(4)*signs(4)));
Ssc             = 1/(4*pi*(R1)^2);
correction = (4*pi*area)/lambda^2;
gammaMetisRCS =  sqrt(correction*Ssc*RCS*(lambda/(4*pi*R2))^2*(1-gammaMetisRCS0)^2);
end

function gammaMMMagic = mmMAGICModel(scenario, iteration, k, m)
% M-mmMAGIC model
% Based on mmMAGIC D2.2 Section 4.6.3 Blockage
% url: https://bscw.5g-mmmagic.eu/pub/bscw.cgi/d202656/mmMAGIC_D2-2.pdf

lambda     = scenario.lambda;
k0         = scenario.k0;
xyzRx0     = scenario.xyzRx0;
R          = scenario.R(k);
deltaR     = iteration.deltaR(:,m);
signs      = iteration.signs;
xyzEdges   = iteration.xyzEdges;
xyzTxImage = iteration.xyzTxImage(m,:);

metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;

edgesVectorsTx = xyzEdges - xyzTxImage;
edgesVectorsRx = xyzRx0 - xyzEdges;
cos_mmMAGIC_phi = dot(edgesVectorsTx, edgesVectorsRx, 2)./(sqrt(sum(abs2(edgesVectorsTx),2)).*sqrt(sum(abs2(edgesVectorsRx),2)));

mmMAGIC_phase_projection = exp(-1j*k0*(sqrt(sum(abs2(xyzEdges-xyzTxImage),2)) + sqrt(sum(abs2(xyzEdges-xyzRx0),2))));
mmMAGIC_phase_R = exp(-1j*k0*R);
mmMAGIC_F = cos_mmMAGIC_phi.*(0.5-metisTerms);
mmMAGICTerms = (0.5 - mmMAGIC_phase_projection/mmMAGIC_phase_R.*mmMAGIC_F);

gammaMMMagic = ((mmMAGICTerms(1)*signs(1)+mmMAGICTerms(2)*signs(2))*(mmMAGICTerms(3)*signs(3)+mmMAGICTerms(4)*signs(4)));
end

function gammaITUFres = ITUfresnelModel(scenario, iteration, m)
% Based on ITU Recommendation P.526-14 Section 5.2.1.1: Fresnel integral method
% url: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.526-14-201801-I!!PDF-E.pdf

lambda     = scenario.lambda;
xyzRx0     = scenario.xyzRx0;
xyzEdges   = iteration.xyzEdges;
xyzTxImage = iteration.xyzTxImage(m,:);
xyzRP      = iteration.xyzRP(m,:);
xyzTx0     = xyzTxImage;

t1 = xyzEdges(2,:) - xyzEdges(1,:);
t1 = t1/norm(t1);
t2 = xyzEdges(4,:) - xyzEdges(3,:);
t2 = t2/norm(t2);
n  = cross(t1,t2);
n  = abs(n/norm(n));

zt = dot(xyzRP-xyzTx0, n);
zr = dot(xyzRx0-xyzRP, n);
xy = xyzEdges - xyzRP;

rho_r1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzRx0.* t1, 2);
rho_t1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzTx0.* t1, 2);
rho_r2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzRx0.* t2, 2);
rho_t2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzTx0.* t2, 2);

phi = zeros(4,1);
phi(1:2) = atan2(rho_r1, zr) + atan2(rho_t1, zt);
phi(3:4) = atan2(rho_r2, zr) + atan2(rho_t2, zt);

x = dot(xy, t1(ones(size(xy,1),1),:), 2);
y = dot(xy, t2(ones(size(xy,1),1),:), 2);
xy = [x(1:2); y(3:4)];

v_ITU = sign(xy).*sqrt(2/lambda.*abs(xy).^1.18*abs(1/zr + 1/zt)^0.18.*abs(phi).^0.82);

Cx = fresnelC(v_ITU(2)) - fresnelC(v_ITU(1));
Cy = fresnelC(v_ITU(4)) - fresnelC(v_ITU(3));
Sx = fresnelS(v_ITU(2)) - fresnelS(v_ITU(1));
Sy = fresnelS(v_ITU(4)) - fresnelS(v_ITU(3));

gammaITUFres = conj((Cx*Sy + Sx*Cy) + 1j*(Sx*Sy - Cx*Cy))*0.5;
end

function result = abs2(x)
result = real(x).^2 + imag(x).^2;
end

function [surface, nVerts] = TAP_CreateRectangle(width,height,nVerts, frequency)
x     = linspace(0,width,nVerts(1));
x     = x - mean(x);
y     = linspace(0,height,nVerts(2));
y     = y - mean(y);
[x,y] = ndgrid(x,y);
surface.vertices = [x(:) y(:) zeros(size(x(:)))];
nVertsTot        = nVerts(1)*nVerts(2);
faces            = (1:(nVertsTot-nVerts(1))).' + [0 1 nVerts(1)+1 nVerts(1)];
faces(nVerts(1):nVerts(1):end,:) = [];
surface.faces    = faces;
surface.tag = 'rectangleSurface';
surface.frequency=frequency;
surface.type = 'rectangular surface';
siz = size(surface.faces);
surface.centroid = reshape(mean(reshape(surface.vertices(surface.faces.',:),[siz(2) siz(1) 3])),[siz(1) 3]);

faces       = surface.faces;
xV = reshape(surface.vertices(faces,1),siz);
yV = reshape(surface.vertices(faces,2),siz);
zV = reshape(surface.vertices(faces,3),siz);

if size(xV,2)==4
    v1 = [xV(:,3) - xV(:,2) yV(:,3) - yV(:,2) zV(:,3) - zV(:,2)];
    v2 = [xV(:,4) - xV(:,3) yV(:,4) - yV(:,3) zV(:,4) - zV(:,3)];
elseif size(xV,2)==3
    v1 = [xV(:,1) - xV(:,1) yV(:,2) - yV(:,1) zV(:,2) - zV(:,1)];
    v2 = [xV(:,3) - xV(:,2) yV(:,3) - yV(:,2) zV(:,3) - zV(:,2)];
end
xp = [v1(:,2).*v2(:,3)-v2(:,2).*v1(:,3) v1(:,3).*v2(:,1)-v2(:,3).*v1(:,1) v1(:,1).*v2(:,2)-v2(:,1).*v1(:,2)];
nV = xp./sqrt(sum(xp.^2,2));
surface.normalVectors = nV;

crossProd           = zeros([size(surface.faces,1),3,size(surface.faces,2)-2]);
numVerts              = size(surface.faces,2);
for k = 1:numVerts
    k2 = mod(k,numVerts)+1;
    crossProd(:,:,k) = cross([xV(:,k ) yV(:,k ) zV(:,k )],...
        [xV(:,k2) yV(:,k2) zV(:,k2)],2);
end
crossProdSum = sum(crossProd,3);
surface.area  = abs(0.5*sum(surface.normalVectors.*crossProdSum,2));

end

function obj = TAP_moveObject(obj, dist)
switch obj.type
    case 'antenna'
        obj.geomT = {dist obj.geomT{:}};
    otherwise
        obj.vertices = bsxfun(@plus, obj.vertices, dist);
        obj.centroid = bsxfun(@plus, obj.centroid, dist);
end
end

function obj = TAP_rotateObject(obj, rotAxis, alpha)
if alpha==0
    rotMatrixR =eye(3);
else
    rotAxis=rotAxis/sqrt(rotAxis*rotAxis');
    A=[0 -rotAxis(3) rotAxis(2);rotAxis(3) 0 -rotAxis(1);-rotAxis(2) rotAxis(1) 0];
    rotMatrixR = eye(3) + A^2 + A*sin(alpha) - A^2*cos(alpha);
    rotMatrixR = rotMatrixR';
end
switch obj.type
    case 'antenna'
        obj.geomT = {[rotAxis alpha] obj.geomT{:}};
    otherwise
        obj.vertices = obj.vertices*rotMatrixR;
        if isfield(obj, 'centroid')
            obj.centroid = obj.centroid*rotMatrixR;
        end
        if isfield(obj, 'normalVectors')
            obj.normalVectors = obj.normalVectors*rotMatrixR;
        end
end
end

function [R,T] = reduce_geomT(geomT)
affineT = eye(4);
for m = 1:numel(geomT)
    thisT = eye(4);
    if numel(geomT{m}) == 3
        thisT(1:3,4) = geomT{m};
    elseif numel(geomT{m}) == 4
        if geomT{m}(4)==0
            thisT(1:3,1:3)=eye(3);
        else
            n=geomT{m}(1:3);
            n=n/sqrt(n*n');
            A=[0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
            thisT(1:3,1:3) = eye(3) + A^2 + A*sin(geomT{m}(4)) - A^2*cos(geomT{m}(4));
        end
    end
    affineT = affineT*thisT;
end
R = affineT(1:3,1:3);
T = affineT(1:3,4)';
end

function source = nardaV637_TAP(f,varargin)
% nardaV637_TAP
% Leuven Rx horn: Narda V637
% https://nardamiteq.com/product-spec/waveguidehornantennas_Standard_Gain_Horns_2.60_to_40_GHz.pdf
f_scaling               = f/28e9;
para.f                  = f;                  % Frequency (Hz)
% Values and estimates from data sheet and data sheet picture
para.a                  = 0.007112/f_scaling; % Feeding wave guide width (m); WR28
para.b                  = 0.003556/f_scaling; % Feeding wave guide height (m); "-
para.width              = 0.02292/f_scaling;  % Aperture width (m); with wall thickness removed (2x2 mm)
para.height             = 0.01683/f_scaling;  % Aperture height (m); with wall thickness removed (2x2 mm)
para.length             = 0.033724/f_scaling; % Length (along the wave guide feed) of the flaring section
para.scaleToDirectivity = [];
[~,~,para]              = horn_pattern(pi/2,0,para);
fhandle = eval([' @(theta,phi)(',func2str(@horn_pattern),'(theta,phi,para));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'nardav637';
source.geomT{1}      = [0 0 0];
source.geomT{2}      = [1 0 0 0];
if nargin>1 && not(isempty(varargin{1}))
    source.geomT{3} = [0 0 1 varargin{1}];
end
end

function source = schwarzbeck9170_TAP(f,varargin)
% Schwarzbeck9170
% Leuven Tx horn: Schwarzbeck 9170
% http://schwarzbeck.de/Datenblatt/k9170.pdf
f_scaling               = f/28e9;
para.f                  = f;     % Frequency (Hz)
% Values and estimates from data sheet and data sheet picture
para.a                  = 0.010668/f_scaling; % Feeding wave guide width (m); assumed WR42
para.b                  = 0.004318/f_scaling; % Feeding wave guide height (m); "-
para.width              = 0.059/f_scaling;    % Aperture width (m); wall thickness removed
para.height             = 0.044/f_scaling;    % Aperture height (m); wall thickness removed
para.length             = 0.062/f_scaling;    % Length (along the wave guide feed) of the flaring section
para.scaleToDirectivity = [];
[~,~,para]              = horn_pattern(pi/2,0,para);
fhandle = eval([' @(theta,phi)(',func2str(@horn_pattern),'(theta,phi,para));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'schwarzbeck9170';
source.geomT{1} = [0 0 0];
source.geomT{2} = [1 0 0 0];
if nargin>1 && not(isempty(varargin{1}))
    source.geomT{3} = [0 0 1 varargin{1}];
end
end

function [E_theta,E_phi,modelParameters] = horn_pattern(theta,phi,varargin)

if not(isempty(varargin{1}))
    f = varargin{1}.f;
    b = varargin{1}.b;
    a = varargin{1}.a;
    height = varargin{1}.height;
    width = varargin{1}.width;
    length = varargin{1}.length;
    scaleToDirectivity = varargin{1}.scaleToDirectivity;
end

modelParameters.f = f;
modelParameters.b = b;
modelParameters.a = a;
modelParameters.height = height;
modelParameters.width = width;
modelParameters.length = length;
modelParameters.scaleToDirectivity = scaleToDirectivity;

if isempty(scaleToDirectivity)
    tmp = modelParameters;
    tmp.scaleToDirectivity = 1;
    th = linspace(0,pi,361);
    ph = linspace(-pi,pi,721);
    dO = dOmega_TAP(th,ph,[0 pi -pi pi]);
    [th,ph]     = ndgrid(th,ph);
    [E_th,E_ph] = horn_pattern(th,ph,tmp);
    average     = dO*(abs(E_th(:)).^2+abs(E_ph(:)).^2)/4/pi;
    scaleToDirectivity = sqrt(1./average);
    modelParameters.scaleToDirectivity = scaleToDirectivity;
end

[E_theta,E_phi] = alg_horn(theta,phi,f,a,b,width,height,length,scaleToDirectivity);
end

function [E_theta,E_phi] = alg_horn(theta,phi,f,a,b,width,height,length,scaleToDirectivity)
c_0 = 299792458;
k  = 2*pi*f/c_0;
ky = k*sin(theta).*sin(phi);
kz = k*cos(theta);

apexE  = height*sqrt((length/(height-b))^2+.25);
sqrt_factorE = 1./sqrt(pi*apexE*k);

t1      = sqrt_factorE*(-.5*k*height-kz*apexE);
t2      = sqrt_factorE*(.5*k*height-kz*apexE);
fresnel = conj((fresnelC(t2)+1j*fresnelS(t2))-(fresnelC(t1)+1j*fresnelS(t1)));
I2      = sqrt(pi*apexE./k).*exp(.5j*kz.^2.*apexE./k).*fresnel;

apexH  = width*sqrt((length/(width-a))^2+.25);
sqrt_factorH = 1./sqrt(pi*apexH*k);

kyprim = ky+pi/width;
t1prim = sqrt_factorH*(-.5*k*width-kyprim*apexH);
t2prim = sqrt_factorH*(.5*k*width-kyprim*apexH);
fresnelprim = conj((fresnelC(t2prim)+1j*fresnelS(t2prim))-(fresnelC(t1prim)+1j*fresnelS(t1prim)));

kybis  = ky-pi/width;
t1bis = sqrt_factorH*(-.5*k*width-kybis*apexH);
t2bis = sqrt_factorH*(.5*k*width-kybis*apexH);
fresnelbis = conj((fresnelC(t2bis)+1j*fresnelS(t2bis))-(fresnelC(t1bis)+1j*fresnelS(t1bis)));
I1 = .5*sqrt(pi*apexH./k).*(...
    exp(.5j*kyprim.^2.*apexH./k).*fresnelprim + ...
    exp(.5j*kybis.^2.*apexH./k).*fresnelbis);

Amp     = I1.*I2*scaleToDirectivity;
E_theta = -Amp.*(cos(phi) + sin(theta));
E_phi   = Amp.*cos(theta).*sin(phi);
end

function dO=dOmega_TAP(theta,phi,varargin)
limits = varargin{1};
w=(2*mod(floor(theta/pi),2)-1).*cos(theta)+2*floor(theta/pi);
w_av=conv(w,[.5 .5]);
w_av(1)=(2*mod(floor(limits(1)/pi),2)-1).*cos(limits(1))+2*floor(limits(1)/pi);
w_av(end)=(2*mod(floor(limits(2)/pi),2)-1).*cos(limits(2))+2*floor(limits(2)/pi);
dw=abs(diff(w_av));
phi_av=conv(phi,[.5 .5]);
phi_av(1)=limits(3);
phi_av(end)=limits(4);
dphi=abs(diff(phi_av));
[dw,dphi]=ndgrid(dw,dphi);
dO=(dw(:).*dphi(:))';
end

function source = isotropic_TAP(varargin)
parameters.pol = [1 0];
fhandle = eval([' @(theta,phi)(',func2str(@isotropic_pattern),'(theta,phi,parameters));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'isotropic';
source.geomT{1}         = [0 0 0];
source.geomT{2}         = [1 0 0 0];
if nargin>0 && not(isempty(varargin{1}))
    source.geomT{3}     = [0 0 1 varargin{1}];
end
end

function [E_theta, E_phi, modelParameters] = isotropic_pattern(theta, phi, varargin)
pol  = [1 0];
modelParameters.pol = pol;
E_theta=pol(1)*ones(size(theta+phi));
E_phi  =pol(2)*ones(size(theta+phi));
end

function [errorTable, extraError] = getErrorMetric(y, models)
% y is the MoM reference
% xEst is the estimation

xEst = models.resultsTh.';

if y.available
    errorTable = zeros(numel(models.List),1);
    [~, b] = size(xEst);   
    idx = 1:size(y.results,1);
    % Remove parts of data where PO or model is non-applicable        
    idx = idx(11:end-10); % For 0.1 deg resolution, remove 1 deg at beginning and end
    
    for i = 1:b
        errorTable(i,1) = 10*log10(NMSE(y.results(idx), xEst(idx, i)));
    end
    metric  = 'momNMSEdB';    
    
    format bank
    T = table(models.List, errorTable,  'VariableNames', {'Model', metric, });
    disp(T)
    format short
else
    errorTable = [];
    extraError = [];
end
end

function NMSE = NMSE(x, xEst)
idx = find(not(isnan(xEst)));
x    = x(idx);
xEst = xEst(idx);
NMSE  = sum(abs(x - xEst).^2)/sum(abs(x).^2);
end


