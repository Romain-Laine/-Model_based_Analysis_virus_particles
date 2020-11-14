% Try lsqnonline solution for radVP; sigPos
%
% This script now converges for the test case - it gets a reasonable radVP
% which (at 60.4 nm) is a much better radVP than the mean value or the peak
% position. 
% But the sigLoc2D value moves only slightly toward its true value (of 25)
% - is this because I'm missing a sqrt(2) 
% - yes, possibly: The 2D Gaussian filter may be 15 nm as the 1D sigma - 
%  but the scatterplot simulation uses it as the 2D value. 
%  - 1-D value in scatterplot sim is more like 25/1.4 ~~ 11 nm. 

% % Define inputs:
obsMeanRad = 62.2888;
obsStdRad  =  15.4338;
% (These were from radVP = 60, sigLoc2D = 20, radLink - 10)

% obsMeanRad = 62.706;
% obsStdRad  =  15.4;
% % (These are from radVP = 60, sigLoc2D = 25)

initRadVP = obsMeanRad;
initSigPos = obsStdRad;

obsMeas = [obsMeanRad,obsStdRad];

x0(1) = initRadVP;
x0(2)= initSigPos;


% x = lsqnonlin( circ2dGaussVirusGetParams,inputs )


% Write an iterative loop:
count = 0;
hXm = 1;   % Step size in radVP, nm, for calculating required step
hXs = 1;

lam = 0.1; % Amount of convergence to attempt.

cXm = initRadVP;  % Inialise current radVP guess
cXs = initSigPos;
dXm = 1; % To start the loop OK
cXlist = zeros(100);
while ( abs(dXm) > 0.01 && count < 100)
    
   tXm = cXm + hXm;
   tXs = cXs + hXs;
   
   [estMeasCurrRad,estMeasCurrStd] = circ2dGaussVirusGetParams([cXm, cXs]);
   [estMeasAltMRad,estMeasAltMStd] = circ2dGaussVirusGetParams([tXm, cXs]);
   [estMeasAltSRad,estMeasAltSStd] = circ2dGaussVirusGetParams([tXm, tXs]);
   
   gradWithM = mean(estMeasAltMRad - estMeasCurrRad)/hXm;
   gradWithS = mean(estMeasAltSStd - estMeasCurrStd)/hXs;
   
   dXm = -( mean(estMeasCurrRad - obsMeas(1)) )/gradWithM;
   dXs = -( mean(estMeasCurrStd - obsMeas(2)) )/gradWithS;
   
   cXm = cXm + lam*dXm;
   cXs = cXs + lam*dXs;
   
   % PUT IN TEST - move if solution is better - solve slow convergence!
   [cXm,cXs]
   
   count = count+1;
   cXlist(count) = cXm;
end

if count<100
  cXlist(count+1:end) = [];
end

figure(4)
plot(cXlist, 'lineWidth',3)
set(gca,'fontSize',14)
xlabel('Iteration','fontSize',14)
ylabel('Virus Protein radius, nm','fontSize',14)
title('Inverse simulation of Virus Protein radius','fontSize',14)

estRadMean
estRadStd