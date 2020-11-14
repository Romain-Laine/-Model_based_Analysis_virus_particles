% Virus STORM data simulation
% 
% 2013 November, EJR, Version 1-1
%
% PURPOSE:
%    Simulate localisation microscopy results for some virus geometries
%    Present and analyse results, and see if input geometry is recoverable 
%
% METHOD 1:
%    Simulate 1 virus in the XY plane only (flag2DcircleVirus = 1)
%    Let Virus Protein (VP) lie on a circular ring, (radius radVP).
%    Let there be numProtPerVir labelled proteins per virus,
%    And the position of these proteins on the circle is random
%     (Could be regularly spaced; the above may describe partial labelling)
%    Let a number (numFluorPerVir) of fluorophores be attached per protein
%    Constant distance, random orientation of fluorophore attachment
%      (Attachment distance can be set to zero to simplify the case)
%    Let localisation error be from a 2D, circularly-symmetric Gaussian 
%      Note the "fixed distance" and "Gaussian distance" may represent 
%      more factors than just "Attachment" and "Localisation Error"
%      In fact, the Gaussian term could be used alone, for all factors.
%    Simulate in units of nm. 
%    (Could convert nm to pixel widths for rainSTORM visualisation input.)
%    Let localisations be generated:
%       Randomly, 
%       Homogenously (all dyes equally likely, with equal error probs),
%       Free from mislocalisations / overlapping image artifacts.
%
%    OBSERVATIONS (Method 1)
%     For radVP    = 60;  
%         radLink  = 10;            
%         sigThomp = 20;
%      mean(radLocs) from known origin is larger than radVP, by about
%      2.3 nm (plus or minus 1 nm, by repeated simulation).
%      This is about a 4% error from method bias (as expected - EJR).
%      Thus VP radius could not be obtained to better than 2 nm by this 
%      analysis, but the data look precise enough to do better...
%     And for:
%         numLocs = 1000; numProt = 500, 1 dye per protein
%       mean(positLocs) is typically 2-3 nm away from the (0,0) origin
%       which indicates the error that would be added by estimating 
%       the displacements from the Centre of Mass of Localisations. 
%       Note that in this simulation the centre is known to lie at (0,0). 
%       A single colour-correct capsid position would have to be very 
%       precise to out-perform 2-3 nm precision, but it may be possible
%       since the capsid can be bright (as it contains many molecules).
%       What about the average error from capsids in multiple viruses ???
%     How  can the ~4% radius-increasing bias be improved upon?
%       (A) By modelling this case, and fitting observed mean / std / kurt 
%       of the radLocs histogram to get parameters (radVP, sigThomp). 
%       (B) By incorporating (measurable) Thompson precision data
%       (C) By incorporating (esimated) antibody size to the analysis
%     And how can this model be improved upon?
%        By making it 3-D, to simulate a spherical shell of virus protein
%        By introducting Elliptical or one-side bias in VP labelling
%          (although simulating a small number of Virus proteins would 
%           automatically generate 1-sided bias in the labelling density).
%        By contraining fluorophores to lie outside the VP radius
%          (which seems to be a physically likely case - but is not 
%           significant if the localisation error is much bigger than the 
%           fluorophore attachment distance, as may be likely).
%
%     Beware dimensions in this model - see "sigLoc2D" math and comments


%
% 1. INPUTS and FLOW CONTROL
%
flag2DcircleVirus = 1;    % Simulate in XY-plane, a ring of virus protein
flagScatterPlot   = 1;    

numLocs        = 50000;   % Number of localisations sampled from virus
numProtPerVir  = 20000;     % Number of proteins on a virus
numFluorPerPro = 1;       % Number of fluorophores per protein
numFluors = numProtPerVir*numFluorPerPro;
radVP    = 60;            % Virus protein radius (from centre of virus)
radLink  = 10;          % Fluorophore attachment radius
sigLoc2D = 20;         % 2d Localisation error (sqrt 2 * Thompson formula)


% 
% 2. SIMULATION OF DATA
% 
if (flag2DcircleVirus) % For a METHOD 1 simulation
 angleProteins = 2*pi *rand(numProtPerVir,1);
 positProteins = radVP * [cos(angleProteins),sin(angleProteins)];
 
 angleAttach = 2*pi * rand(numFluors, 1);
 dispAttach  = radLink * [cos(angleAttach),sin(angleAttach)];
 positProRep = zeros(numFluors,2); % Pre-allocate mem
 
 for lpFl = 1:numFluorPerPro
     positProRep( ((lpFl-1)*numProtPerVir+1):(lpFl*numProtPerVir),: ) =...
        positProteins;
 end
 
 positFluor = positProRep + dispAttach;
 
 % Simulate localisations 1 at a time (simple, but not efficient):
 positLocs = zeros(numLocs,2);          % Allocate memory
 angleLocs = 2*pi * rand(numLocs, 1);   % Angle of localisation error
 disOfLocs = sigLoc2D*randn(numLocs,1); % Distance of localisation error
 dispLocs  = [disOfLocs.*cos(angleLocs), disOfLocs.*sin(angleLocs)];
 for lpLoc = 1:numLocs
     rFN              = ceil( rand * numFluors); % A random fluorophore
     positLocs(lpLoc,:) = positFluor(rFN,:) + dispLocs(lpLoc,:);
 end
 
end

% Determine radial positions of fluorophores and localisations
radFluors = sqrt( positFluor(:,1).^2 + positFluor(:,2).^2 );
radLocs   = sqrt( positLocs(:,1).^2 + positLocs(:,2).^2 );

% Now try fitting the parameters of a model to radLocs
% Could try NLLS error; or max likelyhood... 



%
% 3. Visualisation of results
%
if (flagScatterPlot)
 figure(1)
 scatter(positProteins(:,1),positProteins(:,2),'b');
 axis equal
 
 hold on
 scatter(positFluor(:,1),positFluor(:,2),'gx');
 scatter(positLocs(:,1),positLocs(:,2),'rx');
 title('Positions','fontSize',14)
 legend('Protein','Fluorophore','Localisation');
 xlabel('X, nm','fontSize',14)
 ylabel('Y, nm','fontSize',14)
 set(gca,'fontSize',14)
 
 hold off
end
 
figure(2)
hist(radFluors, (radVP-radLink-1:1:radVP+radLink+1) );
title('Fluorophore radial distance','fontSize',14)
xlabel('Radial Distance, nm','fontSize',14)
ylabel('Number','fontSize',14)
set(gca,'fontSize',14)

figure(3)
hist(radLocs,50);
title('Localisation radial distance','fontSize',14)
xlabel('Radial Distance, nm','fontSize',14)
ylabel('Number','fontSize',14)
set(gca,'fontSize',14)

mean(radLocs)
std(radLocs)
% mean(positLocs) % Bias of centre position is determined by CoM(locs)