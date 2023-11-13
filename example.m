clc; clear;
% close all
addpath(genpath('utils/'))
set(0,'DefaultFigureWindowStyle','docked')
dbstop if error

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
mesh_dir = '../data/PAIR_TEST_TEMP/';
shapeTarget_name = 'target.off';
shapeSource_name = 'source.off';
shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);

%% center data
shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

%% preprocess the meshes
shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});
shapeSource = MESH.preprocess(shapeSource,meshOptions{:});

%% compute the WKS descriptors
BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2); 

fctTarget_all = heatKernelSignature(BTarget, EvTarget, shapeTarget.A, numTimes); %fctTarget_all = fMAP.waveKernelSignature(BTarget, EvTarget, shapeTarget.A, numTimes);
fctSource_all = heatKernelSignature(BSource, EvSource, shapeSource.A, numTimes);    %fctSource_all = fMAP.waveKernelSignature(BSource, EvSource, shapeSource.A, numTimes);
fctTarget = fctTarget_all(:,1:numSkip:end);
fctSource = fctSource_all(:,1:numSkip:end);

%% Plot WKS descriptors on shapes 
numPlots = 50; %numPlots = size(fctTarget,2);
numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots / numRows);

% Plotting the target shape with all selected WKS functions as subplots
figure('name', 'WKS target');
verticesTarget = shapeTarget.surface.VERT;
facesTarget = shapeTarget.surface.TRIV;

for iPlot=1:numPlots
    i = iPlot;
    subplot(numRows, numCols,iPlot);
    h = trisurf(facesTarget, verticesTarget(:,1), verticesTarget(:,2), verticesTarget(:,3), fctTarget(:,i), 'FaceColor', 'interp');
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
verticesSource = shapeSource.surface.VERT;
facesSource = shapeSource.surface.TRIV;

for iPlot=1:numPlots
    i = iPlot;
    subplot(numRows, numCols,iPlot);
    h = trisurf(facesSource, verticesSource(:,1), verticesSource(:,2), verticesSource(:,3), fctSource(:,i), 'FaceColor', 'interp');
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
use_lm = true;
load_lm = false;

if(use_lm)
    if(load_lm)
[file,path] = uigetfile('*.txt','Select the landmarks based descriptors file for source shape');
lm_fct_Source = load([path,file]);

[file,path] = uigetfile('*.txt','Select the landmarks based descriptors file for target shape');
lm_fct_Target = load([path,file]);

%% Select only landmarks based descriptors 1,2,6,7,11,12,16,17,21,22,26,27
lm_fct_Target = lm_fct_Target(:,[1,6,11,16,21,26]);
lm_fct_Source = lm_fct_Source(:,[1,6,11,16,21,26]);

%lm_fct_Target = lm_fct_Target(:,[1,2,6,7,11,12,16,17,21,22,26,27]);
%lm_fct_Source = lm_fct_Source(:,[1,2,6,7,11,12,16,17,21,22,26,27]);

    else
        % Select the file containing the landmarks indices
        [file,path] = uigetfile('*.txt','Select the landmarks indices file');
        lm_idx = load([path,file]);

        % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
        lm_idx_Target = lm_idx(:,1);
        lm_idx_Source = lm_idx(:,2);

        % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
        timesteps_lm = 200;

        lm_fct_Target = fMAP.compute_descriptors_with_landmarks(shapeTarget,100,lm_idx_Target,timesteps_lm,1);
        lm_fct_Source = fMAP.compute_descriptors_with_landmarks(shapeSource,100,lm_idx_Source,timesteps_lm,1);

        % Keep values at index 1,2,3,4,5,6 and then regularly skip values
        nb_lm = size(lm_idx,1);
        skip_timestep_lm = 50;
        total_skip_lm = skip_timestep_lm * nb_lm;
        initial_idx = 1:nb_lm;
        idx = initial_idx;

        for i=1:(timesteps_lm/skip_timestep_lm-1)
            idx = [idx, initial_idx + i*total_skip_lm];
        end

        lm_fct_Target = lm_fct_Target(:,idx);
        lm_fct_Source = lm_fct_Source(:,idx);
        normc(lm_fct_Target), normc(lm_fct_Source);

    end


%% Plot shapes with landmarks as points
figure('name', 'Landmarks target');
% loop over all the landmarks and plot them as points in subplots
numPlot = length(lm_idx);
numRows = ceil(sqrt(numPlot));
numCols = ceil(numPlot / numRows);
for i=1:numPlot
    subplot(numRows, numCols,i);
    h = trisurf(facesTarget, verticesTarget(:,1), verticesTarget(:,2), verticesTarget(:,3), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    hold on;
    scatter3(verticesTarget(lm_idx_Target(i),1),verticesTarget(lm_idx_Target(i),2),verticesTarget(lm_idx_Target(i),3),100,'filled');
    title(['Landmark ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

figure('name', 'Landmarks source');
% loop over all the landmarks and plot them as points in subplots
for i=1:numPlot
    subplot(numRows, numCols,i);
    h = trisurf(facesSource, verticesSource(:,1), verticesSource(:,2), verticesSource(:,3), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    hold on;
    scatter3(verticesSource(lm_idx_Source(i),1),verticesSource(lm_idx_Source(i),2),verticesSource(lm_idx_Source(i),3),100,'filled');
    title(['Landmark ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

%% Plot landmarks based descriptors on shapes
% Plotting the target shape with all selected landmarks based functions as subplots
figure('name', 'Landmarks based target');
numPlots = size(lm_fct_Target,2);
numRows = ceil(sqrt(numPlots));
numCols = ceil(numPlots / numRows);

for i=1:numPlots
    subplot(numRows, numCols,i);
    h = trisurf(facesTarget, verticesTarget(:,1), verticesTarget(:,2), verticesTarget(:,3), lm_fct_Target(:,i), 'FaceColor', 'interp');
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
    h = trisurf(facesSource, verticesSource(:,1), verticesSource(:,2), verticesSource(:,3), lm_fct_Source(:,i), 'FaceColor', 'interp');
    set(h,'edgecolor','none');
    title(['Landmarks based function ', num2str(i)]);
    view(-200,-30);
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %colorbar; % Display the color bar indicating the color scale
end

drawnow;

%% Append the landmarks based descriptors to the WKS descriptors
fctTarget = [lm_fct_Target]; %fctTarget = [fctTarget,lm_fct_Target];
fctSource = [lm_fct_Source]; %fctSource = [fctSource,lm_fct_Source];

end

%% optimize the functional map using the standard or the complex resolvent Laplacian term
[C_target2source, M_old] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');
[C_target2source_slant, M_slant] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');
[C_target2source_new, M_new] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'complRes');
T_source2target = fMAP.fMap2pMap(BTarget,BSource,C_target2source);
T_source2target_slant = fMAP.fMap2pMap(BTarget,BSource,C_target2source_slant);
T_source2target_new = fMAP.fMap2pMap(BTarget,BSource,C_target2source_new);

%%
% visualize the computed maps
figure();
subplot(1,3,1);
MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title('standard Mask');
subplot(1,3,2);
MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant,plotOptions{:}); title('slanted Mask');
subplot(1,3,3);
MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('complex resolvent Mask');

% visualize the mask 
figure();
subplot(1,3,1); imagesc(M_old); axis square;  title('Standard Laplacian Mask');
subplot(1,3,2); imagesc(M_slant); axis square;  title('Standard Laplacian Mask');
subplot(1,3,3); imagesc(M_new); axis square;  title('Complex Resolvent Laplacian Mask');
