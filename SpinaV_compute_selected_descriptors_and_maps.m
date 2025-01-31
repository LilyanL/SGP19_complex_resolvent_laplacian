%% Let's setup everything nicely
clc; clear all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode

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
listFolders = {'PAIR_006'}; %{'PAIR_001', 'PAIR_002', 'PAIR_003'}; %, 'PAIR_004', 'PAIR_005'};

% list of methods to compute the descriptors
listMethods = {'WKS', 'HKS'}; %, 'SIHKS', 'EKS', 'WKS+SIHKS', 'WKS+EKS', 'SIHKS+EKS', 'WKS+SIHKS+EKS'};

%% Display all the pairs of meshes and landmarks
for i = 1:length(listFolders)
    mesh_dir = ['../data/' listFolders{i} '/'];
    shapeTarget_name = 'target.off';
    shapeSource_name = 'source.off';
    shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
    shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);
    
    %% center data
    shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
    shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

    landmarks_file = ['../data/' listFolders{i} '/landmarks.txt'];
    lm_idx = load(landmarks_file)+1;

    % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
    lm_idx_Target = lm_idx(:,1);
    lm_idx_Source = lm_idx(:,2);

    % Assign a color to each landmark using the colormap 'jet'
    colorMap = jet(length(lm_idx_Target));

    %Display both shapes with their landmarks
    plotName = ['Shapes with landmarks - Folder ' listFolders{i}];
    figure('Name', plotName,'NumberTitle','off');
    subplot(1,2,1);
    h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
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
    h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
    scatter3(shapeTarget.surface.VERT(lm_idx_Target,1), shapeTarget.surface.VERT(lm_idx_Target,2), shapeTarget.surface.VERT(lm_idx_Target,3), 100, colorMap, 'filled');
    title(['Target shape (' listFolders{i} ')']);

    % Add labels to the landmarks
    for j = 1:length(lm_idx_Target)
        text(shapeTarget.surface.VERT(lm_idx_Target(j),1), shapeTarget.surface.VERT(lm_idx_Target(j),2), shapeTarget.surface.VERT(lm_idx_Target(j),3), num2str(j), 'FontSize', 14);
    end

end

%% For each method, compute the global descriptors for each shape
for nbMethod = 1:length(listMethods)
    method = listMethods{nbMethod};
    fprintf('Computing global %s descriptors for each shape...\n', method);
    for i = 1:length(listFolders)
        mesh_dir = ['../data/' listFolders{i} '/'];
        shapeTarget_name = 'target.off';
        shapeSource_name = 'source.off';
        shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
        shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);

        % center data
        shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
        shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

        % preprocess the meshes
        shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});
        shapeSource = MESH.preprocess(shapeSource,meshOptions{:});

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
            h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), fctSource(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure('Name', ['Folder ' listFolders{i} ' - Global descriptors ' method  ' - Target'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = (j-1)*numSkipGlobalDescriptors+1; % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), fctTarget(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

    end
end

%% For each method, compute and plot the local descriptors for each shape. For HKS, keep the first descriptor only. For WKS, keep the selected descritor if different from zero
descriptor_for_computation_source = [];
descriptor_for_computation_target = [];

for nbMethod = 1:length(listMethods)
    method = listMethods{nbMethod};
    fprintf('Computing local %s descriptors for each shape...\n', method);
    for i = 1:length(listFolders)
        mesh_dir = ['../data/' listFolders{i} '/'];
        shapeTarget_name = 'target.off';
        shapeSource_name = 'source.off';
        shapeTarget = MESH.MESH_IO.read_shape([mesh_dir, shapeTarget_name]);
        shapeSource = MESH.MESH_IO.read_shape([mesh_dir, shapeSource_name]);
        
        % center data
        shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);
        shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

        % preprocess the meshes
        shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});
        shapeSource = MESH.preprocess(shapeSource,meshOptions{:});

        BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
        EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2); 

        % load landmarks indices
        landmarks_file = ['../data/' listFolders{i} '/landmarks.txt'];
        lm_idx = load(landmarks_file)+1;

        % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
        lm_idx_Target = lm_idx(:,1);
        lm_idx_Source = lm_idx(:,2);

        % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
        timesteps_lm = 100;

        %compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
        numEigs = 100;

        for j = 1:length(lm_idx_Target)
            lm_fct_Target = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeTarget,numEigs,lm_idx_Target(j),timesteps_lm,1, method);
            lm_fct_Source = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx_Source(j),timesteps_lm,1, method);

            if strcmp(method, 'WKS')
                descriptor_for_computation_source = [descriptor_for_computation_source, lm_fct_Source(:,33)];
                descriptor_for_computation_target = [descriptor_for_computation_target, lm_fct_Target(:,1)];
            elseif strcmp(method, 'HKS')
                descriptor_for_computation_source = [descriptor_for_computation_source, lm_fct_Source(:,1)];
                descriptor_for_computation_target = [descriptor_for_computation_target, lm_fct_Target(:,1)];
            end
        end
        
    end
end

descriptor_for_computation_source = normc(descriptor_for_computation_source);
descriptor_for_computation_target = normc(descriptor_for_computation_target);

%% Plot all the descriptors for each shape
 % plot the descriptors
 numPlots = size(descriptor_for_computation_source,2);
 numRows = ceil(sqrt(numPlots));
 numCols = ceil(numPlots / numRows);

 figure('Name', ['Folder ' listFolders{i} ' - Selected descriptors '  ' - Source'],'NumberTitle','off');
 for j = 1:numPlots
     numberOfDescriptor = j; % number of descriptor, taking into account the ignored ones
     subplot(numRows, numCols, j);
     h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), descriptor_for_computation_source(:,j), 'FaceColor', 'interp');
     set(h, 'edgecolor', 'none');
     axis equal; axis off; hold on;
     title(['Descriptor ' num2str(numberOfDescriptor)]);
 end

 figure('Name', ['Folder ' listFolders{i} ' - Selected descriptors '  ' - Target'],'NumberTitle','off');
 for j = 1:numPlots
     numberOfDescriptor = j; % number of descriptor, taking into account the ignored ones
     subplot(numRows, numCols, j);
     h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), descriptor_for_computation_target(:,j), 'FaceColor', 'interp');
     set(h, 'edgecolor', 'none');
     axis equal; axis off; hold on;
     title(['Descriptor ' num2str(numberOfDescriptor)]);
 end


%% Choose the descriptors to use for the functional map
        fctTarget = [descriptor_for_computation_target]; %fctTarget = [fctTarget,lm_fct_Target];
        fctSource = [descriptor_for_computation_source]; %fctSource = [fctSource,lm_fct_Source];

        %% optimize the functional map using the standard or the complex resolvent Laplacian term
        [C_target2source, M_old] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');
        [C_target2source_slant, M_slant] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');
        [C_target2source_new, M_new] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'complRes');
        T_source2target = fMAP.fMap2pMap(BTarget,BSource,C_target2source);
        T_source2target_slant = fMAP.fMap2pMap(BTarget,BSource,C_target2source_slant);
        T_source2target_new = fMAP.fMap2pMap(BTarget,BSource,C_target2source_new);

        % visualize the computed maps
        figure('Name', ['Folder ' listFolders{i} ' - FM using local descriptors ' method],'NumberTitle','off');
        subplot(1,3,1);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title('standard Mask');
        subplot(1,3,2);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_slant,plotOptions{:}); title('slanted Mask');
        subplot(1,3,3);
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('complex resolvent Mask');


%% Export all figures to png
mkdir results;
destinationFolder = './results/';
saveAllFigures(destinationFolder, 'png');

%% Check memory
memory

        