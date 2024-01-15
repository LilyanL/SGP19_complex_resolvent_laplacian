%% Let's setup everything nicely
clc; clear all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode

% create folders to export results and figures
mkdir results;
destinationFolder = '.\results\';

% Create a subfolder to store the maps
maps_dir = [destinationFolder 'maps\'];
mkdir(maps_dir);

% Initialization to store data
pairs_array = {};

%% some parameters
% fmap size: k2-by-k1
k1 = 100;
k2 = 100;

% params to compute WKS or descriptors
numTimesGlobalDescriptors = 200;
numSkipGlobalDescriptors = 20;

% relative weights of different terms to optimize a functional map
para.a = 2e-1;
para.b = 1e-2;
para.c = 8e-4;  % weight for the Laplacian mask term!
para.alpha = 1e-1;

% params to pre-process the meshes
meshOptions = {'IfComputeGeoDist',false,'IfComputeLB',true,'IfComputeNormals',true,'numEigs',100};

% params to visualize the maps
plotOptions = {'IfShowCoverage',false,'OverlayAxis','y','cameraPos',[0,90]};

% list of folders contained in ..\data\ from which we load the meshes and the landmarks to compute the maps
listFolders = {'PAIR_001', 'PAIR_002'};% 'PAIR_003', 'PAIR_004', 'PAIR_005'};

% list of methods to compute the descriptors
listMethods = {'WKS', 'HKS'}; %, 'SIHKS', 'EKS', 'WKS+SIHKS', 'WKS+EKS', 'SIHKS+EKS', 'WKS+SIHKS+EKS'};

%% Display options
displayDescriptorsLocal = false;
displayDescriptorsGlobal = false;
displayBasisFunctions = false;

%% Display all the pairs of meshes and landmarks and pre-process the meshes
% Create cell array of absolute paths to the meshes
foldersPaths = cell(1,length(listFolders));
for i = 1:length(listFolders)
    foldersPaths{i} = [pwd '\..\data\' listFolders{i} '\'];
end

pairs_array = CollectionLoadShapes(foldersPaths, meshOptions, true, true);

%% Compute and display the basis functions for each shape
if(displayBasisFunctions)
    CollectionDisplayBasisFunctions(pairs_array, k1, k2);
end

%% For each method, compute the global descriptors for each shape


% Call the function
pairs_array = CollectionComputeGlobalDescriptors(pairs_array, listMethods, k1, k2, numTimesGlobalDescriptors, numSkipGlobalDescriptors, displayDescriptorsGlobal);

%% For each method, compute and plot the local descriptors for each shape
num_skip = 15;
timesteps_lm = 100;
pairs_array = CollectionComputeLocalDescriptors(pairs_array, listMethods, k1, k2, timesteps_lm, num_skip, displayDescriptorsLocal);

%% Compute the functional maps using different descriptors
for i = 1:length(pairs_array)
    % load the pair of shapes from the array
    curPairShapes = pairs_array{i};

    shapeSource = curPairShapes.shape_source;
    shapeTarget = curPairShapes.shape_target;
    curFolderName = curPairShapes.pair_folder_name;

    % for each descriptor combination method, compute the functional map
    for nbMethod = 1:length(listMethodsMaps)
        method = listMethodsMaps{nbMethod};

        % load the descriptors from the PairShapes object
        fctTarget = []; fctSource = []; methodString = [];
        for j = 1:length(method)
            methodString = [methodString method{j}{1} ' ' method{j}{2} ' + '];
            if strcmp(method{j}{2}, 'global')
                fctTarget = [fctTarget curPairShapes.descriptors_global_target{strcmp(curPairShapes.descriptors_global_target_labels, method{j}{1})}];
                fctSource = [fctSource curPairShapes.descriptors_global_source{strcmp(curPairShapes.descriptors_global_source_labels, method{j}{1})}];
            elseif strcmp(method{j}{2}, 'local')
                fctTarget = [fctTarget curPairShapes.descriptors_local_target{strcmp(curPairShapes.descriptors_local_target_labels, method{j}{1})}];
                fctSource = [fctSource curPairShapes.descriptors_local_source{strcmp(curPairShapes.descriptors_local_source_labels, method{j}{1})}];
            end
        end
        methodString = methodString(1:end-3); % remove the last ' + '


        BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
        EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2);

        fprintf('Computing the functional map using %s descriptors...\n', methodString);

        % optimize the functional map using the standard or the complex resolvent Laplacian term
        fprintf('Computing the functional map using %s descriptors and the standard Laplacian term...\n', methodString);
        [C_target2source, M_old] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');

        fprintf('Computing the functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
        [C_target2source_slant, M_slant] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');

        fprintf('Computing the functional map using %s descriptors and the complex resolvent Laplacian term...\n', methodString);
        [C_target2source_new, M_new] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'complRes');

        T_source2target = fMAP.fMap2pMap(BTarget,BSource,C_target2source);
        T_source2target_slant = fMAP.fMap2pMap(BTarget,BSource,C_target2source_slant);
        T_source2target_new = fMAP.fMap2pMap(BTarget,BSource,C_target2source_new);

        % visualize the computed maps
        figure('Name', ['Folder ' curFolderName ' - FM using descriptors ' methodString],'NumberTitle','off');
        subplot(1,3,1);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title('standard Mask');
        subplot(1,3,2);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant,plotOptions{:}); title('slanted Mask');
        subplot(1,3,3);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('complex resolvent Mask');

        % Store the mappings in the PairShapes object
        curPairShapes.mappings{nbMethod} = [T_source2target; T_source2target_slant; T_source2target_new];
        curPairShapes.mappings_Labels{nbMethod} ={[methodString ' - standard'], [methodString ' - slant'], [methodString ' - complRes']};

        % refine the mapping with zoomOut (final_map = refineZoomOut(initial_matches, initial_dim, S1, S2)
        T_source2target_new = refineZoomOut(T_source2target_new, size(C_target2source_new,1), shapeSource, shapeTarget, ['Folder ' curFolderName ' - FM using descriptors ' methodString ' + zoomOut']);
        curPairShapes.mappings{nbMethod} = [curPairShapes.mappings{nbMethod}; T_source2target_new];
        curPairShapes.mappings_Labels{nbMethod} = [curPairShapes.mappings_Labels{nbMethod} ' - zoomOut'];


        % on a new figure, visualize the mapping with complex resolvent Laplacian term and the drilling paths
        figure('Name', ['Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (direct)'],'NumberTitle','off');

        %Display both shapes with the drilling paths
        %plotName = ['Shapes with drilling paths - Folder ' listFolders{i} ' - FM using descriptors ' methodString ' (direct')];
        sourceTitle = ['Source shape (' curFolderName ')'];
        targetTitle = ['Target shape (' curFolderName ')'];

        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'direct');

        figure('Name', ['Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (connex)'],'NumberTitle','off');
        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'connex', [7 7]);

        % export the maps to text files for later use
        map_name = [maps_dir 'map_' curFolderName '_' methodString '_standard.txt'];
        dlmwrite(map_name, T_source2target, 'delimiter', ' ');
        map_name = [maps_dir 'map_' curFolderName '_' methodString '_slant.txt'];
        dlmwrite(map_name, T_source2target_slant, 'delimiter', ' ');
        map_name = [maps_dir 'map_' curFolderName '_' methodString '_complRes.txt'];
        dlmwrite(map_name, T_source2target_new, 'delimiter', ' ');

    end
end


%% Export all figures
saveAllFigures(destinationFolder, 'png');

%% Check memory
memory