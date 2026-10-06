% voltage in kV
V = [2.0 2.5 3.0 3.5 4.0 4.5 5.0]';

% for the ring measurements each row in the matrix is associated with a
% voltage the top is 2.0 and the bottom is 5.0

% inner ring measurements in cm
innerMin = [
    3.00 2.80 3.08 3.08 3.10;
    2.75 2.59 2.50 2.62 2.58;
    2.47 2.38 2.42 2.49 2.39;
    2.19 2.26 2.34 2.30 2.22;
    2.09 2.14 1.93 1.88 2.05;
    1.89 1.76 1.82 1.83 1.71;
    1.80 1.61 1.55 1.61 1.77
];

innerMax = [
    3.35 3.41 3.47 3.48 3.60;
    3.00 3.05 3.01 3.04 3.06;
    2.96 2.90 2.86 2.92 3.05;
    2.60 2.69 2.55 2.65 2.50;
    2.55 2.60 2.39 2.41 2.49;
    2.40 2.25 2.29 2.34 2.18;
    2.02 2.19 2.09 2.11 2.26
];

% outer ring measurements in cm
outerMin = [
    5.41 5.39 5.48 5.41 5.30;
    4.82 4.89 4.88 4.88 4.79;
    4.50 4.48 4.48 4.55 4.46;
    4.08 4.19 4.10 3.95 4.14;
    3.82 3.69 3.71 3.65 3.68;
    3.32 3.29 3.10 3.12 3.19;
    3.26 3.25 3.39 3.28 3.37
];

outerMax = [
    5.96 5.79 5.92 5.85 5.96;
    5.41 5.20 5.29 5.30 5.31;
    5.08 4.99 4.90 4.83 5.04;
    4.46 4.62 4.49 4.47 4.71;
    4.38 4.29 4.21 4.26 4.20;
    3.90 4.09 4.01 3.95 3.90;
    3.53 3.67 3.50 3.61 3.62
];

% constants
h = 6.62607015e-34;
me = 9.1093837e-31;
e = 1.602176634e-19;
L = 0.125;      % target to screen, m
dL = 0.002;     % uncertainty in L, m
R = 0.065;      % radius of bulb, m

% midpoint of the two edges
innerTrials = (innerMin + innerMax)/2;
outerTrials = (outerMin + outerMax)/2;

% mean of the 5 trials
Dinner_cm = mean(innerTrials, 2);
Douter_cm = mean(outerTrials, 2);

% standard deviation of inner and outer
innerRandom = std(innerTrials, 0, 2);
outerRandom = std(outerTrials, 0, 2);

% thickness of the rings 
innerWidth = mean((innerMax-innerMin)/2, 2);
outerWidth = mean((outerMax-outerMin)/2, 2);

% caliper uncertainty
edgeError = 0.05;
instrumentError = edgeError/sqrt(2)

% total uncertainty in measured D for inner and outer calculated in
% quadrature
DinnerErr_cm = sqrt(innerRandom.^2 + innerWidth.^2 + instrumentError.^2)
DouterErr_cm = sqrt(outerRandom.^2 + outerWidth.^2 + instrumentError.^2)

% convert to meters
Dinner = Dinner_cm/100;
Douter = Douter_cm/100;
DinnerErr = DinnerErr_cm/100;
DouterErr = DouterErr_cm/100;

% correct each measured diameter for the curved bulb
[DinnerCorr, DinnerCorrErr] = correctDiameter(Dinner, DinnerErr, L, R);
[DouterCorr, DouterCorrErr] = correctDiameter(Douter, DouterErr, L, R);

% convert corrected diameters back to cm for graph
DinnerCorr_cm = DinnerCorr*100;
DouterCorr_cm = DouterCorr*100;
DinnerCorrErr_cm = DinnerCorrErr*100;
DouterCorrErr_cm = DouterCorrErr*100;

x = 1 ./ sqrt(V);

% power supply had 4 marks between each whole kV on the meter so there are 5 spaces per kV
dV = (1/5)/2;

% uncertainty in x = V^(-1/2)
xErr = 0.5 * dV ./ (V.^(3/2))

[mInner, mInnerErr] = weightedFit(x, DinnerCorr_cm, DinnerCorrErr_cm);
[mOuter, mOuterErr] = weightedFit(x, DouterCorr_cm, DouterCorrErr_cm);
mInnerErr
mOuterErr
% constant when V is in kV
lambda = h / sqrt(2*me*e*1000);

% slope is in cm*sqrt(kV), convert to meters
mInner_m = mInner/100;
mOuter_m = mOuter/100;

% D = (2*L*lambda/d)*(1/sqrt(V))
dInner = (2*L*lambda)/mInner_m;
dOuter = (2*L*lambda)/mOuter_m;

% convert to nm
dInner = dInner*1e9;
dOuter = dOuter*1e9;

dInnerFitErr = dInner*(mInnerErr/mInner);
dOuterFitErr = dOuter*(mOuterErr/mOuter);
dInnerLErr = dInner*(dL/L);
dOuterLErr = dOuter*(dL/L);

% combine length and fit error in quadrature
dInnerErr = sqrt(dInnerFitErr^2 + dInnerLErr^2);
dOuterErr = sqrt(dOuterFitErr^2 + dOuterLErr^2);

acceptedInner = 0.213;
acceptedOuter = 0.123;

innerPercentError = abs(dInner-acceptedInner)/acceptedInner*100;
outerPercentError = abs(dOuter-acceptedOuter)/acceptedOuter*100;

innerPercentUncertainty = dInnerErr/dInner*100;
outerPercentUncertainty = dOuterErr/dOuter*100;

fprintf('\nInner ring:\n')
fprintf('slope = %.4f cm*sqrt(kV)\n', mInner)
fprintf('d = %.4f +/- %.4f nm\n', dInner, dInnerErr)
fprintf('percent uncertainty = %.2f %%\n', innerPercentUncertainty)
fprintf('percent error = %.2f %%\n', innerPercentError)

fprintf('\nOuter ring:\n')
fprintf('slope = %.4f cm*sqrt(kV)\n', mOuter)
fprintf('d = %.4f +/- %.4f nm\n', dOuter, dOuterErr)
fprintf('percent uncertainty = %.2f %%\n', outerPercentUncertainty)
fprintf('percent error = %.2f %%\n', outerPercentError)

%plotting with error bars
figure
errorbar(x, DinnerCorr_cm, DinnerCorrErr_cm, DinnerCorrErr_cm, xErr, xErr, 'o', 'LineStyle', 'none')
hold on
errorbar(x, DouterCorr_cm, DouterCorrErr_cm, DouterCorrErr_cm, xErr, xErr, 's', 'LineStyle', 'none')

xfit = linspace(0,max(x)*1.05,100);
plot(xfit, mInner*xfit, 'LineWidth', 1.5)
plot(xfit, mOuter*xfit, 'LineWidth', 1.5)

xlabel('$\frac{1}{\sqrt{V_a}} \quad (\frac{1}{\sqrt{kV}})$', 'Interpreter', 'latex')
ylabel('Corrected ring diameter D (cm)', 'Interpreter', 'latex')
title('Electron Diffraction Experiment')
legend('Inner Ring Error Bars', 'Outer Ring Error Bars', 'Inner Ring Fit', 'Outer Ring Fit', 'Location', 'northwest')

function [Dcorr, DcorrErr] = correctDiameter(D, DErr, L, R)
    % distance from center of bulb to ring position
    q = sqrt(R^2-(D/2).^2);
    % amount the curved screen moves toward the target
    correction = R-q;
    % actual target-to-ring horizontal distance
    Leff = L-correction;
    % diameter corrected back to distance L
    Dcorr = D .* (L./Leff);
    % propagate D uncertainty numerically
    Dplus = D + DErr;
    Dminus = D - DErr;
    qPlus = sqrt(R^2-(Dplus/2).^2);
    qMinus = sqrt(R^2-(Dminus/2).^2);
    Lplus = L-(R-qPlus);
    Lminus = L-(R-qMinus);
    corrPlus = Dplus .* (L./Lplus);
    corrMinus = Dminus .* (L./Lminus);
    DcorrErr = abs(corrPlus-corrMinus)/2;
end

function [m, mErr] = weightedFit(x, y, yErr)
    w = 1./yErr.^2;
    % weighted fit through the origin
    m = sum(w.*x.*y) / sum(w.*x.^2);
    mErr = 1/sqrt(sum(w.*x.^2));
end