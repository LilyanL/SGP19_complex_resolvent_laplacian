%% Let's setup everything nicely
clc; clear all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode

% create folders to export results and figures
mkdir results;
destinationFolder = './results/';

% Create a subfolder to store the maps
maps_dir = [destinationFolder 'maps/'];
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

% list of folders contained in ../data/ from which we load the meshes and the landmarks to compute the maps
listFolders = {'PAIR_001', 'PAIR_002'};% 'PAIR_003', 'PAIR_004', 'PAIR_005'};

% list of methods to compute the descriptors
listMethods = {'WKS', 'HKS'}; %, 'SIHKS', 'EKS', 'WKS+SIHKS', 'WKS+EKS', 'SIHKS+EKS', 'WKS+SIHKS+EKS'};

%% Display all the pairs of meshes and landmarks
for i = 1:length(listFolders)

    curPairShapes = PairShapes;

    %% Load the meshes
    mesh_dir = ['../data/' listFolders{i} '/'];
    shapeTarget_name = 'target.off';
    shapeSource_name = 'source.off';
    shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
    shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);

    %% center data
    shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
    shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

    % preprocess the meshes
    shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});
    shapeSource = MESH.preprocess(shapeSource,meshOptions{:});

    landmarks_file = ['../data/' listFolders{i} '/landmarks.txt'];
    lm_idx = load(landmarks_file)+1; % +1 because matlab starts indexing at 1

    % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
    lm_idx_Target = lm_idx(:,1);
    lm_idx_Source = lm_idx(:,2);

    % Assign a color to each landmark using the colormap 'jet'
    colorMap = jet(length(lm_idx_Target));

    %Display both shapes with their landmarks
    plotName = ['Shapes with landmarks - Folder ' listFolders{i}];
    figure('Name', plotName,'NumberTitle','off');
    subplot(1,2,1);
    display_shape(shapeSource);
    scatter3(shapeSource.surface.VERT(lm_idx_Source,1), shapeSource.surface.VERT(lm_idx_Source,2), shapeSource.surface.VERT(lm_idx_Source,3), 100, colorMap, 'filled');
    title(['Source shape (' listFolders{i} ')']);

    % Add labels to the landmarks
    for j = 1:length(lm_idx_Source)
        text(shapeSource.surface.VERT(lm_idx_Source(j),1), shapeSource.surface.VERT(lm_idx_Source(j),2), shapeSource.surface.VERT(lm_idx_Source(j),3), num2str(j), 'FontSize', 14);
    end

    % Add a text box under the plot with the number of landmarks
    dim = [.2 .05 .3 .3];
    str = ['Number of landmarks: ' num2str(length(lm_idx_Target))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');

    subplot(1,2,2);
    display_shape(shapeTarget);
    scatter3(shapeTarget.surface.VERT(lm_idx_Target,1), shapeTarget.surface.VERT(lm_idx_Target,2), shapeTarget.surface.VERT(lm_idx_Target,3), 100, colorMap, 'filled');
    title(['Target shape (' listFolders{i} ')']);

    % Add labels to the landmarks
    for j = 1:length(lm_idx_Target)
        text(shapeTarget.surface.VERT(lm_idx_Target(j),1), shapeTarget.surface.VERT(lm_idx_Target(j),2), shapeTarget.surface.VERT(lm_idx_Target(j),3), num2str(j), 'FontSize', 14);
    end

    % If '/drilling_paths.txt' exists, load the drilling paths
    readPaths = [];
    drillingIndex_entry = [];
    drillingIndex_exit = [];
    drilling_paths_file = [mesh_dir 'drilling_paths_source.txt'];
    if exist(drilling_paths_file, 'file') == 2
        readPaths = readmatrix([mesh_dir 'drilling_paths_source.txt']);
        drillingIndex_entry = readPaths(:,1);
        drillingIndex_exit = readPaths(:,2);
    end

    % Plot the drilling paths if they exist
    subplot(1,2,2);
    colorMap = jet(length(drillingIndex_entry));
    for j = 1:length(drillingIndex_entry)
        entryPoint = shapeTarget.surface.VERT(drillingIndex_entry(j),:);
        endPoint = shapeTarget.surface.VERT(drillingIndex_exit(j),:);
        plotTrajectory(entryPoint,endPoint, colorMap(j,:), 3, 0.5);
    end

    % Store all useful data in the PairShapes object
    curPairShapes.shape_source = shapeSource;
    curPairShapes.shape_target = shapeTarget;

    curPairShapes.landmarks_source = lm_idx_Source;
    curPairShapes.landmarks_target = lm_idx_Target;

    curPairShapes.trajectories_source = readPaths;

    % Store the PairShapes object in the array
    pairs_array{i} = curPairShapes;

end

%% Compute and display the basis functions for each shape
% for i = 1:length(listFolders)
%     mesh_dir = ['../data/' listFolders{i} '/'];
%     shapeTarget_name = 'target.off';
%     shapeSource_name = 'source.off';
%     shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
%     shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);
%
%     % center data
%     shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
%     shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);
%
%     % preprocess the meshes
%     shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});
%     shapeSource = MESH.preprocess(shapeSource,meshOptions{:});
%
%     % compute descriptors
%     BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
%     EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2);
%
%     % plot the descriptors
%     numPlots = k1;
%     numRows = ceil(sqrt(numPlots));
%     numCols = ceil(numPlots / numRows);
%
%     figure('Name', ['Folder ' listFolders{i}  '- Basis functions '  ' - Source'],'NumberTitle','off');
%     for j = 1:numPlots
%         subplot(numRows, numCols, j);
%         h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), BSource(:,j), 'FaceColor', 'interp');
%         set(h, 'edgecolor', 'none');
%         axis equal; axis off; hold on;
%         title(['Basis function ' num2str(j)]);
%     end
%
%     figure('Name', ['Folder ' listFolders{i}  '- Basis functions ' ' - Target'],'NumberTitle','off');
%     for j = 1:numPlots
%         subplot(numRows, numCols, j);
%         h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), BTarget(:,j), 'FaceColor', 'interp');
%         set(h, 'edgecolor', 'none');
%         axis equal; axis off; hold on;
%         title(['Basis function ' num2str(j)]);
%     end
%
% end

%% For each method, compute the global descriptors for each shape
for nbMethod = 1:length(listMethods)
    method = listMethods{nbMethod};
    fprintf('Computing global %s descriptors for each shape...\n', method);

    for i = 1:length(listFolders)
        % load the pair of shapes from the array
        curPairShapes = pairs_array{i};
        shapeSource = curPairShapes.shape_source;
        shapeTarget = curPairShapes.shape_target;

        % compute descriptors
        BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
        EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2);

        if strcmp(method, 'WKS')
            fctTarget_all = fMAP.waveKernelSignature(BTarget, EvTarget, shapeTarget.A, numTimesGlobalDescriptors);
            fctSource_all = fMAP.waveKernelSignature(BSource, EvSource, shapeSource.A, numTimesGlobalDescriptors);
        elseif strcmp(method, 'HKS')
            fctTarget_all = heatKernelSignature(BTarget, EvTarget, shapeTarget.A, numTimesGlobalDescriptors);
            fctSource_all = heatKernelSignature(BSource, EvSource, shapeSource.A, numTimesGlobalDescriptors);
        end

        % ignore some descriptors
        fctTarget = fctTarget_all(:,1:numSkipGlobalDescriptors:end);
        fctSource = fctSource_all(:,1:numSkipGlobalDescriptors:end);

        % plot the descriptors
        numPlots = size(fctTarget,2);
        numRows = ceil(sqrt(numPlots));
        numCols = ceil(numPlots / numRows);

        figure('Name', ['Folder ' listFolders{i} ' - Global descriptors ' method  ' - Source'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = (j-1)*numSkipGlobalDescriptors+1; % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            display_shape(shapeSource, fctSource(:,j));
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure('Name', ['Folder ' listFolders{i} ' - Global descriptors ' method  ' - Target'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = (j-1)*numSkipGlobalDescriptors+1; % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            display_shape(shapeTarget, fctTarget(:,j));
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        % Store the descriptors in the PairShapes object
        curPairShapes.descriptors_global_source{nbMethod} = fctSource;
        curPairShapes.descriptors_global_target{nbMethod} = fctTarget;

        % Store the descriptors labels in the PairShapes object (WKS or HKS)
        curPairShapes.descriptors_global_source_labels{nbMethod} = method;
        curPairShapes.descriptors_global_target_labels{nbMethod} = method;

        % Update the pair of shapes in the array
        pairs_array{i} = curPairShapes;

    end
end

%% For each method, compute and plot the local descriptors for each shape
for nbMethod = 1:length(listMethods)
    method = listMethods{nbMethod};
    fprintf('Computing local %s descriptors for each shape...\n', method);
    for i = 1:length(listFolders)
        % load the pair of shapes from the array
        curPairShapes = pairs_array{i};
        shapeSource = curPairShapes.shape_source;
        shapeTarget = curPairShapes.shape_target;

        % compute descriptors
        BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
        EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2);

        % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
        lm_idx_Target = curPairShapes.landmarks_target;
        lm_idx_Source = curPairShapes.landmarks_source;

        % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
        timesteps_lm = 100;

        %compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
        numEigs = 100;
        lm_fct_Target = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeTarget,numEigs,lm_idx_Target,timesteps_lm,1, method);
        lm_fct_Source = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx_Source,timesteps_lm,1, method);

        num_skip = 15;
        firstTimeSteps = 1+(1- 1:(size(lm_idx,1)))*timesteps_lm;
        idx=[];
        for lm_nb = 1:size(lm_idx,1)
            added_idx = firstTimeSteps(lm_nb):num_skip:(firstTimeSteps(lm_nb+1)-1);
            idx = [idx added_idx];
        end

        lm_fct_Target = lm_fct_Target(:,idx);
        lm_fct_Source = lm_fct_Source(:,idx);
        normc(lm_fct_Target); normc(lm_fct_Source);

        % plot the descriptors
        numPlots = size(lm_fct_Target,2);
        numRows = ceil(sqrt(numPlots));
        numCols = ceil(numPlots / numRows);

        figure('Name', ['Folder ' listFolders{i} ' - Local descriptors ' method  ' - Source'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = idx(j); % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            display_shape(shapeSource, lm_fct_Source(:,j));
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure('Name', ['Folder ' listFolders{i} ' - Local descriptors ' method  ' - Target'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = idx(j); % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            display_shape(shapeTarget, lm_fct_Target(:,j));
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        % Store the descriptors in the PairShapes object
        curPairShapes.descriptors_local_source{nbMethod} = lm_fct_Source;
        curPairShapes.descriptors_local_target{nbMethod} = lm_fct_Target;

        % Store the descriptors labels in the PairShapes object (WKS or HKS)
        curPairShapes.descriptors_local_source_labels{nbMethod} = method;
        curPairShapes.descriptors_local_target_labels{nbMethod} = method;

        % Update the pair of shapes in the array
        pairs_array{i} = curPairShapes;
    end

    %% Compute the functional maps using different descriptors
    for i = 1:length(listFolders)
        % load the pair of shapes from the array
        curPairShapes = pairs_array{i};

        shapeSource = curPairShapes.shape_source;
        shapeTarget = curPairShapes.shape_target;

        % load the descriptors from the PairShapes object
        fctTarget = curPairShapes.descriptors_local_target{nbMethod};
        fctSource = curPairShapes.descriptors_local_source{nbMethod};

        % optimize the functional map using the standard or the complex resolvent Laplacian term

        fprintf('Computing the functional map using %s descriptors and the standard Laplacian term...\n', method);
        [C_target2source, M_old] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');

        fprintf('Computing the functional map using %s descriptors and the slanted Laplacian term...\n', method);
        [C_target2source_slant, M_slant] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');

        fprintf('Computing the functional map using %s descriptors and the complex resolvent Laplacian term...\n', method);
        [C_target2source_new, M_new] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'complRes');

        T_source2target = fMAP.fMap2pMap(BTarget,BSource,C_target2source);
        T_source2target_slant = fMAP.fMap2pMap(BTarget,BSource,C_target2source_slant);
        T_source2target_new = fMAP.fMap2pMap(BTarget,BSource,C_target2source_new);

        % Store the mappings in the PairShapes object
        curPairShapes.mappings{nbMethod} = [T_source2target; T_source2target_slant; T_source2target_new];
        curPairShapes.mappings_Labels{nbMethod} ={[method ' - standard'], [method ' - slant'], [method ' - complRes']};

        % visualize the computed maps
        figure('Name', ['Folder ' listFolders{i} ' - FM using local descriptors ' method],'NumberTitle','off');
        subplot(1,3,1);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title('standard Mask');
        subplot(1,3,2);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant,plotOptions{:}); title('slanted Mask');
        subplot(1,3,3);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('complex resolvent Mask');

        % on a new figure, visualize the mapping with complex resolvent Laplacian term and the drilling paths
        figure('Name', ['Folder ' listFolders{i} ' - FM using local descriptors ' method ' and drilling paths'],'NumberTitle','off');

        %Display both shapes with their landmarks
        plotName = ['Shapes with landmarks - Folder ' listFolders{i}];
        subplot(1,2,1);
        display_shape(shapeSource);
        title(['Source shape (' listFolders{i} ')']);

        subplot(1,2,2);
        display_shape(shapeTarget);
        title(['Target shape (' listFolders{i} ')']);

        % Get the drilling paths from the PairShapes object for the target shape
        drillingIndex_entry_source = curPairShapes.trajectories_source(:,1);
        drillingIndex_exit_source = curPairShapes.trajectories_source(:,2);

        % Compute the entry and exit points of the drilling paths on the source shape using the mapping
        drillingIndex_entry_target = T_source2target_new(drillingIndex_entry_source);
        drillingIndex_exit_target = T_source2target_new(drillingIndex_exit_source);

        
        % Plot the drilling paths on the source shape
        subplot(1,2,1);
        colorMap = jet(length(drillingIndex_entry_source));
        for j = 1:length(drillingIndex_entry_source)
            entryPoint = shapeSource.surface.VERT(drillingIndex_entry_source(j),:);
            endPoint = shapeSource.surface.VERT(drillingIndex_exit_source(j),:);
            plotTrajectory(entryPoint,endPoint, colorMap(j,:), 3, 0.5);
        end

        % Plot the drilling paths on the target shape
        subplot(1,2,2);
        colorMap = jet(length(drillingIndex_entry_target));
        for j = 1:length(drillingIndex_entry_target)
            entryPoint = shapeTarget.surface.VERT(drillingIndex_entry_target(j),:);
            endPoint = shapeTarget.surface.VERT(drillingIndex_exit_target(j),:);
            plotTrajectory(entryPoint,endPoint, colorMap(j,:), 3, 0.5);
        end

        % export the maps to text files for later use
        map_name = [maps_dir 'map_' listFolders{i} '_' method '_standard.txt'];
        dlmwrite(map_name, T_source2target, 'delimiter', ' ');
        map_name = [maps_dir 'map_' listFolders{i} '_' method '_slant.txt'];
        dlmwrite(map_name, T_source2target_slant, 'delimiter', ' ');
        map_name = [maps_dir 'map_' listFolders{i} '_' method '_complRes.txt'];
        dlmwrite(map_name, T_source2target_new, 'delimiter', ' ');

    end
end


%% Export all figures
saveAllFigures(destinationFolder, 'png');

%% Check memory
memory