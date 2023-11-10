clc; close all; clear;
addpath(genpath('utils/'))
set(0,'DefaultFigureWindowStyle','docked')

%% some parameters
% fmap size: k2-by-k1
k1 = 100;
k2 = 100;
% params to compute WKS descriptors
% #WKS = length(1:numSkip:numTimes)
numTimes = 100;
numSkip = 1;
% relative weights of different terms to optimize a functional map
para.a = 2e-1;
para.b = 1e-2;
para.c = 8e-4;  % weight for the Laplacian mask term!
para.alpha = 1e-1;
% params to pre-process the meshes
meshOptions = {'IfComputeGeoDist',false,'IfComputeLB',true,'IfComputeNormals',true,'numEigs',100};
% params to visualize the maps
plotOptions = {'IfShowCoverage',false,'OverlayAxis','y','cameraPos',[0,90]};
%% read the mesh
mesh_dir = '../data/PAIR_003/';
s1_name = 'target.off';
s2_name = 'source.off';
S1 = MESH.MESH_IO.read_shape([mesh_dir, s1_name]);
S2 = MESH.MESH_IO.read_shape([mesh_dir, s2_name]);

%% center data
S1.surface.VERT = S1.surface.VERT - mean(S1.surface.VERT);
S2.surface.VERT = S2.surface.VERT - mean(S2.surface.VERT);

%% preprocess the meshes
S1 = MESH.preprocess(S1,meshOptions{:});
S2 = MESH.preprocess(S2,meshOptions{:});

%% compute the WKS descriptors
B1 = S1.evecs(:,1:k1); B2 = S2.evecs(:,1:k2);
Ev1 = S1.evals(1:k1); Ev2 = S2.evals(1:k2); 

fct1_all = heatKernelSignature(B1, Ev1, S1.A, numTimes); %fct1_all = fMAP.waveKernelSignature(B1, Ev1, S1.A, numTimes);
fct2_all = heatKernelSignature(B2, Ev2, S2.A, numTimes);    %fct2_all = fMAP.waveKernelSignature(B2, Ev2, S2.A, numTimes);
fct1 = fct1_all(:,1:numSkip:end);
fct2 = fct2_all(:,1:numSkip:end);

%% Plot WKS descriptors on shapes 
numPlots = 50; %numPlots = size(fct1,2);
numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots / numRows);

% Plotting the target shape with all selected WKS functions as subplots
figure('name', 'WKS target');
verticesTarget = S1.surface.VERT;
facesTarget = S1.surface.TRIV;

for iPlot=1:numPlots
    i = iPlot;
    subplot(numRows, numCols,iPlot);
    h = trisurf(facesTarget, verticesTarget(:,1), verticesTarget(:,2), verticesTarget(:,3), fct1(:,i), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    title(['WKS function ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

% Plotting the source shape with all selected WKS functions as subplots
figure('name', 'WKS source');
verticesSource = S2.surface.VERT;
facesSource = S2.surface.TRIV;

for iPlot=1:numPlots
    i = iPlot;
    subplot(numRows, numCols,iPlot);
    h = trisurf(facesSource, verticesSource(:,1), verticesSource(:,2), verticesSource(:,3), fct2(:,i), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    title(['WKS function ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

%% Select a file to load landmarks based descriptors
use_lm = false;

if(use_lm)
[file,path] = uigetfile('*.txt','Select the landmarks based descriptors file for source shape');
lm_fct_2 = load([path,file]);

[file,path] = uigetfile('*.txt','Select the landmarks based descriptors file for target shape');
lm_fct_1 = load([path,file]);

%% Select only landmarks based descriptors 1,2,6,7,11,12,16,17,21,22,26,27
lm_fct_1 = lm_fct_1(:,[1,6,11,16,21,26]);
lm_fct_2 = lm_fct_2(:,[1,6,11,16,21,26]);

%lm_fct_1 = lm_fct_1(:,[1,2,6,7,11,12,16,17,21,22,26,27]);
%lm_fct_2 = lm_fct_2(:,[1,2,6,7,11,12,16,17,21,22,26,27]);


%% Plot landmarks based descriptors on shapes
% Plotting the target shape with all selected landmarks based functions as subplots
figure('name', 'Landmarks based target');
numPlots = size(lm_fct_1,2);
numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots / numRows);

for i=1:numPlots
    subplot(numRows, numCols,i);
    h = trisurf(facesTarget, verticesTarget(:,1), verticesTarget(:,2), verticesTarget(:,3), lm_fct_1(:,i), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    title(['Landmarks based function ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

% Plotting the source shape with all selected landmarks based functions as subplots
figure('name', 'Landmarks based source');

for i=1:numPlots
    subplot(numRows, numCols,i);
    h = trisurf(facesSource, verticesSource(:,1), verticesSource(:,2), verticesSource(:,3), lm_fct_2(:,i), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    title(['Landmarks based function ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end



%% Append the landmarks based descriptors to the WKS descriptors
fct1 = [lm_fct_1]; %fct1 = [fct1,lm_fct_1];
fct2 = [lm_fct_2]; %fct2 = [fct2,lm_fct_2];

end

%% optimize the functional map using the standard or the complex resolvent Laplacian term
[C12, M_old] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'standard');
[C12_slant, M_slant] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'slant');
[C12_new, M_new] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'complRes');
T21 = fMAP.fMap2pMap(B1,B2,C12);
T21_slant = fMAP.fMap2pMap(B1,B2,C12_slant);
T21_new = fMAP.fMap2pMap(B1,B2,C12_new);

%%
% visualize the computed maps
figure();
subplot(1,3,1);
MESH.PLOT.visualize_map_colors(S2,S1,T21,plotOptions{:}); title('standard Mask');
subplot(1,3,2);
MESH.PLOT.visualize_map_colors(S2,S1,T21_slant,plotOptions{:}); title('slanted Mask');
subplot(1,3,3);
MESH.PLOT.visualize_map_colors(S2,S1,T21_new,plotOptions{:}); title('complex resolvent Mask');

% visualize the mask 
figure();
subplot(1,3,1); imagesc(M_old); axis square;  title('Standard Laplacian Mask');
subplot(1,3,2); imagesc(M_slant); axis square;  title('Standard Laplacian Mask');
subplot(1,3,3); imagesc(M_new); axis square;  title('Complex Resolvent Laplacian Mask');
