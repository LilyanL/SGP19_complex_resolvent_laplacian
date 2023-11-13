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
numSkipGlobalDescriptors = 10;

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
%     figure('Name', ['Basis functions - Folder ' listFolders{i} '- Source'],'NumberTitle','off');
%     for j = 1:numPlots
%         subplot(numRows, numCols, j);
%         h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), BSource(:,j), 'FaceColor', 'interp');
%         set(h, 'edgecolor', 'none');
%         axis equal; axis off; hold on;
%         title(['Basis function ' num2str(j)]);
%     end
% 
%     figure('Name', ['Basis functions - Folder ' listFolders{i} '- Target'],'NumberTitle','off');
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
        
        figure('Name', ['Global descriptors ' method ' - Folder ' listFolders{i} '- Source'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = (j-1)*numSkipGlobalDescriptors+1; % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), fctSource(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure('Name', ['Global descriptors ' method ' - Folder ' listFolders{i} '- Target'],'NumberTitle','off');
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

%% For each method, compute and plot the local descriptors for each shape
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

        % load landmarks indices
        landmarks_file = ['../data/' listFolders{i} '/landmarks.txt'];
        lm_idx = load(landmarks_file)+1;

        % Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
        lm_idx_Target = lm_idx(:,1);
        lm_idx_Source = lm_idx(:,2);

        % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
        timesteps_lm = 10;

        %compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
        numEigs = 100;
        lm_fct_Target = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeTarget,numEigs,lm_idx_Target,timesteps_lm,1, method);
        lm_fct_Source = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx_Source,timesteps_lm,1, method);

        num_skip = 1;
        % % Keep first values and then regularly skip values
        % nb_lm = size(lm_idx,1);
        % skip_timestep_lm = 1;
        % total_skip_lm = skip_timestep_lm * nb_lm;
        % initial_idx = 1:nb_lm;
        % idx = initial_idx;
        % 
        % for skipStep=1:(timesteps_lm/skip_timestep_lm-1)
        %     idx = [idx, initial_idx + skipStep*total_skip_lm];
        % end

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

        figure('Name', ['Local descriptors ' method ' - Folder ' listFolders{i} '- Source'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = idx(j); % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_Source(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure('Name', ['Local descriptors ' method ' - Folder ' listFolders{i} '- Target'],'NumberTitle','off');
        for j = 1:numPlots
            numberOfDescriptor = idx(j); % number of descriptor, taking into account the ignored ones
            subplot(numRows, numCols, j);
            h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), lm_fct_Target(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

    end
end

        