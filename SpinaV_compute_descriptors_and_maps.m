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
%listFolders = {'PAIR_001','PAIR_002','PAIR_003', 'PAIR_004', 'PAIR_005'};% };
listFolders = {'PAIR_002_lowres'};

% list of methods to compute the descriptors
listMethods = {'WKS', 'HKS'}; %, 'SIHKS', 'EKS', 'WKS+SIHKS', 'WKS+EKS', 'SIHKS+EKS', 'WKS+SIHKS+EKS'};

% list of methods to compute the maps stored as cell arrays of strings, each cell array indicating the method and the locality of the descriptors
listMethodsMaps = {
    %{{'WKS', 'global'}};
    %{{'WKS', 'local'}};
    %{{'WKS', 'global'}, {'WKS', 'local'}};
    %
    %{{'HKS', 'global'}};
    %{{'HKS', 'local'}};
    %{{'HKS', 'global'}, {'HKS', 'local'}};
%
    %{{'WKS', 'local'}, {'HKS', 'local'}};
    %{{'WKS', 'global'}, {'HKS', 'global'}, {'WKS', 'local'}, {'HKS', 'local'}};
    {{'WKS', 'global'}, {'HKS', 'local'}};
};

%% Display options
displayShapePairs = true;
displayShapePairsWithPaths = true;
displayBasisFunctions = true;
displayDescriptorsGlobal = true;
displayDescriptorsLocal = true;

%% Display all the pairs of meshes and landmarks and pre-process the meshes
% Create cell array of absolute paths to the meshes
foldersPaths = cell(1,length(listFolders));
for i = 1:length(listFolders)
    foldersPaths{i} = [pwd '\..\data\' listFolders{i} '\'];
end

loadOption = 'default'; % 'default' or 'source_only' or 'target_only
pairs_array = CollectionLoadShapes(foldersPaths, meshOptions, displayShapePairs, displayShapePairsWithPaths, loadOption);

%% Duplicate the pairs of shapes to apply random noise and transforms to the copies and store them in a temporary cell array
% For each pair of shapes, create nbCopies copies with random noise and transforms using function transformShape (limits are provided as parameters)
% transformShape(mesh, noiseMagnitude, rotationMin, rotationMax, translationMin, translationMax)
nbCopies = 2;
pairs_array_tmp = cell(1, length(pairs_array)*nbCopies);
noiseMagnitude = 0; % 0.5;
rotationMin = -pi/4; %-pi/4;
rotationMax = pi/4; %pi/4;
translationMin = -300; %-300;
translationMax = 300; %300;

for i = 1:length(pairs_array)
    display(['Creating duplicates for pair ' num2str(i) ' out of ' num2str(length(pairs_array))]);
    curPairShapes = pairs_array{i};
    for j = 1:nbCopies
        curPairShapesCopy = curPairShapes;

        % if it is the first copy, no noise or transform is applied to have a clean copy for comparison
        if(j==1)
            pairs_array_tmp{(i-1)*nbCopies+j} = curPairShapesCopy;
            continue;
        end

        % add noise and transform the source shape
        curPairShapesCopy.shape_source = transformShape(curPairShapesCopy.shape_source, noiseMagnitude, translationMin, translationMax, rotationMin, rotationMax);
        [curPairShapesCopy.shape_source, noiseVectorSource, transformSource] = transformShape(curPairShapesCopy.shape_source, noiseMagnitude, rotationMin, rotationMax, translationMin, translationMax);

        % add noise and transform the target shape
        curPairShapesCopy.shape_target = transformShape(curPairShapesCopy.shape_target, noiseMagnitude, translationMin, translationMax, rotationMin, rotationMax);
        [curPairShapesCopy.shape_target, noiseVectorTarget, transformTarget] = transformShape(curPairShapesCopy.shape_target, noiseMagnitude, rotationMin, rotationMax, translationMin, translationMax);

        % store the noise vectors and the transforms in the PairShapes object
        curPairShapesCopy.noise_vector_source = noiseVectorSource;
        curPairShapesCopy.noise_vector_target = noiseVectorTarget;
        curPairShapesCopy.transform_source = transformSource;
        curPairShapesCopy.transform_target = transformTarget;

        pairs_array_tmp{(i-1)*nbCopies+j} = curPairShapesCopy;
    end
end

pairs_array = pairs_array_tmp;
clear pairs_array_tmp;

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
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString],'NumberTitle','off');
        subplot(1,3,1);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title('standard Mask');
        subplot(1,3,2);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant,plotOptions{:}); title('slanted Mask');
        subplot(1,3,3);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('complex resolvent Mask');

        % Store the mappings in the PairShapes object
        curPairShapes.mappings{nbMethod} = [T_source2target; T_source2target_slant; T_source2target_new];
        curPairShapes.mappings_Labels{nbMethod} ={[methodString ' - standard'], [methodString ' - slant'], [methodString ' - complRes']};

        % compute the functional and geometric map from source to target
        fprintf('Computing the reverse functional map using %s descriptors and the standard Laplacian term...\n', methodString);
        [C_source2target, M_old] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'standard');

        fprintf('Computing the reverse functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
        [C_source2target_slant, M_slant] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'slant');

        fprintf('Computing the reverse functional map using %s descriptors and the complex resolvent Laplacian term...\n', methodString);
        [C_source2target_new, M_new] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'complRes');

        T_target2source = fMAP.fMap2pMap(BSource,BTarget,C_source2target);
        T_target2source_slant = fMAP.fMap2pMap(BSource,BTarget,C_source2target_slant);
        T_target2source_new = fMAP.fMap2pMap(BSource,BTarget,C_source2target_new);

        % TODO: store the maps in the PairShapes object

        % refine the mapping with BCICP ([T21, T12] = bcicp_refine(S1,S2,B1,B2,T21_ini, T12_ini,num_iter))
        [T_target2source_bcicp, T_source2target_bcicp] = bcicp_refine(shapeSource, shapeTarget, BSource, BTarget, T_target2source, T_source2target,  10);
        [T_target2source_slant_bcicp, T_source2target_slant_bcicp] = bcicp_refine(shapeSource, shapeTarget, BSource, BTarget, T_target2source_slant, T_source2target_slant, 10);
        [T_target2source_new_bcicp, T_source2target_new_bcicp] = bcicp_refine(shapeSource, shapeTarget, BSource, BTarget, T_target2source_new, T_source2target_new, 10);

        % visualize the computed maps with BCICP
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and BCICP'],'NumberTitle','off');
        subplot(1,3,1);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_bcicp,plotOptions{:}); title('standard Mask');
        subplot(1,3,2);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant_bcicp,plotOptions{:}); title('slanted Mask');
        subplot(1,3,3);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new_bcicp,plotOptions{:}); title('complex resolvent Mask');

        % refine the mapping with zoomOut (final_map = refineZoomOut(initial_matches, initial_dim, S1, S2)
        T_source2target_new = refineZoomOut(T_source2target_new, size(C_target2source_new,1), shapeSource, shapeTarget, [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' + zoomOut']);
        curPairShapes.mappings{nbMethod} = [curPairShapes.mappings{nbMethod}; T_source2target_new];
        curPairShapes.mappings_Labels{nbMethod} = [curPairShapes.mappings_Labels{nbMethod} ' - zoomOut'];


        % on a new figure, visualize the mapping with complex resolvent Laplacian term and the drilling paths
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (direct)'],'NumberTitle','off');

        %Display both shapes with the drilling paths
        %plotName = ['Shapes with drilling paths - Folder ' listFolders{i} ' - FM using descriptors ' methodString ' (direct')];
        sourceTitle = ['Source shape (' curFolderName ')'];
        targetTitle = ['Target shape (' curFolderName ')'];

        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'direct');

        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (connex)'],'NumberTitle','off');
        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'connex', [7 7]);

        % export the maps to text files for later use
        map_name = [maps_dir num2str(i) ' - map_' curFolderName '_' methodString '_standard.txt'];
        dlmwrite(map_name, T_source2target, 'delimiter', ' ');
        map_name = [maps_dir num2str(i) ' - map_' curFolderName '_' methodString '_slant.txt'];
        dlmwrite(map_name, T_source2target_slant, 'delimiter', ' ');
        map_name = [maps_dir num2str(i) ' - map_' curFolderName '_' methodString '_complRes.txt'];
        dlmwrite(map_name, T_source2target_new, 'delimiter', ' ');

    end
end


%% Export all figures
saveAllFigures(destinationFolder, 'png');

%% Save all data
save("results\all_data.mat");
%% Check memory
memory