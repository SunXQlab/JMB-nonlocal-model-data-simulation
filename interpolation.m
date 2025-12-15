clear; clc; close all;

set(0, 'defaultfigurecolor', 'w');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'defaultTextInterpreter', 'latex');


% Import SPOT  data
csv1 = importdata('spot.csv');
spotX = csv1.data(:, 2);  % X coordinates of SPOTs
spotY = csv1.data(:, 3);  % Y coordinates of SPOTs

% Import cell type proportion data
csv2 = importdata('st_ct_proportion.csv');

% Extract phenotype data
phenotypeData = csv1.data(:, 1);  % HIF1A expression values

% Extract cell type proportions 
malignantData = csv2.data(:, 1) ;     % Malignant cell proportion
endothelialData = csv2.data(:, 5);      % Endothelial cell proportion

%% Create grid for interpolation
% Define grid from 0 to 600 with 6-unit spacing
gridStep = 6;
gridMax = 600;
gridX = 0:gridStep:gridMax;
gridY = 0:gridStep:gridMax;

% Create meshgrid for interpolation
[X, Y] = meshgrid(gridX, gridY);
gridPointsX = X(:);  % Flatten to column vectors
gridPointsY = Y(:);

%% Identify grid points close to SPOT locations
% Find grid points within 80 units of any SPOT
distanceThreshold = 60;
closeToSpot = false(size(gridPointsX));

% Vectorized distance calculation for efficiency
for i = 1:length(gridPointsX)
    distances = sqrt((gridPointsX(i) - spotX).^2 + (gridPointsY(i) - spotY).^2);
    if any(distances < distanceThreshold)
        closeToSpot(i) = true;
    end
end

% Get indices of grid points far from SPOTs
farGridIndices = find(~closeToSpot);
farGridX = gridPointsX(farGridIndices);
farGridY = gridPointsY(farGridIndices);

% Combine SPOT points with distant grid points
allX = [spotX; farGridX];
allY = [spotY; farGridY];

% Create zero values for distant grid points
zeroValues = zeros(size(farGridX));

%% Prepare data arrays for interpolation
% Combine actual SPOT data with zeros for distant grid points
phenotypeAll = [phenotypeData; zeroValues];
malignantAll = [malignantData; zeroValues];
endothelialAll = [endothelialData; zeroValues];

%% Perform scattered data interpolation
% Interpolate phenotype data using natural neighbor interpolation
F1 = scatteredInterpolant(allX, allY, phenotypeAll, 'natural');
W1 = F1(X, Y);
W1 = W1 / max(W1(:));  % Normalize to [0, 1]

% Interpolate malignant cell proportion
F2 = scatteredInterpolant(allX, allY, malignantAll, 'natural');
W2 = F2(X, Y);

% Interpolate endothelial cell proportion
F3 = scatteredInterpolant(allX, allY, endothelialAll, 'natural');
W3 = F3(X, Y);

%% Convert pixel coordinates to physical units
% Assuming each pixel corresponds to 0.005 physical units
scaleFactor = 0.005 / gridStep;
X_scaled = X * scaleFactor;
Y_scaled = Y * scaleFactor;

%% Visualize interpolated data
% Figure 1: Phenotypic state heatmap
figure(1);
pcolor(X_scaled, Y_scaled, W1);
shading flat;
colorbar
 xticks([0 0.5])
 yticks([0 0.5])
 xlabel('$r_1$');
 ylabel('$r_2$');
axis square 
title('phenotypic state');

% Figure 2: Malignant cell  heatmap
figure(2);
pcolor(X_scaled, Y_scaled, W2 * 1e8);
shading flat;
colorbar
 xticks([0 0.5])
 yticks([0 0.5])
 xlabel('$r_1$');
 ylabel('$r_2$');
axis square 
title('malignant cell ');

% Figure 3: Endothelial cell proportion heatmap
figure(3);
pcolor(X_scaled, Y_scaled, W3);
shading flat;
colorbar
 xticks([0 0.5])
 yticks([0 0.5])
 xlabel('$r_1$');
 ylabel('$r_2$');
axis square 
title('endothelial cell');



% Save phenotype data
dlmwrite('phenotype.csv', W1, 'precision', '%.6f', 'delimiter', ',');
    
% Save malignant cell data
dlmwrite('malignant.csv', W2 * 1e8 , 'precision', '%.6f', 'delimiter', ',');
    
% Save endothelial cell data
dlmwrite('endothelial.csv', W3, 'precision', '%.6f', 'delimiter', ',');
    
% Identify vascular points (endothelial > 0.05)
vascularMask = W3 > 0.05;
vascularX = X(vascularMask) / gridStep;
vascularY = Y(vascularMask) / gridStep;
    
% Save vascular positions
vascularPos = [vascularX(:), vascularY(:)];
dlmwrite('vascular_pos.csv', vascularPos', 'precision', '%.6f', 'delimiter', ',');
    
% Plot vascular points
figure(5);
plot(vascularX, vascularY, 'b.', 'MarkerSize', 10);
axis square;
grid on;
xlabel('X (Physical Units)');
ylabel('Y (Physical Units)');
title('Vascular Point Locations');


%% Refresh display
drawnow;
