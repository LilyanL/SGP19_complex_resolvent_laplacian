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
k1 = 50;
k2 = 50;

% params to compute WKS or descriptors
numTimesGlobalDescriptors = 100;
numSkipGlobalDescriptors = 1;

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
listFolders = {'PAIR_008'};
%listFolders = {'PAIR_001','PAIR_002'};

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
    {{'WKS', 'global'}, {'HKS', 'local'}, {'complexResolvent'}, {'BCICP'}, {'zoomOut'}};
    %{{'WKS', 'global'}, {'HKS', 'local'}, {'complexResolvent'}, {'zoomOut'}};
    };

%% Display options
displayShapePairs = true;
displayShapePairsWithPaths = true;
displayBasisFunctions = false;
displayDescriptorsGlobal = false;
displayDescriptorsLocal = false;

%% Figures used to display the final results
figErrors = figure('Name', 'Errors', 'NumberTitle', 'off');
mean_errors_array = {};

%% Display all the pairs of meshes and landmarks and pre-process the meshes
% Create cell array of absolute paths to the meshes
fprintf(' \n------------------------------ ');
fprintf('Loading and displaying the shapes and landmarks...');
fprintf(' ------------------------------ \n\n');

foldersPaths = cell(1,length(listFolders));
for i = 1:length(listFolders)
    foldersPaths{i} = [pwd '\..\data\' listFolders{i} '\'];
end

loadOption = 'default'; % 'default' or 'source_only' or 'target_only
pairs_array = CollectionLoadShapes(foldersPaths, meshOptions, displayShapePairs, displayShapePairsWithPaths, loadOption);

drawnow;
%% Duplicate the pairs of shapes to apply random noise and transforms to the copies and store them in a temporary cell array
% For each pair of shapes, create nbCopies copies with random noise and transforms using function transformShape (limits are provided as parameters)
% transformShape(mesh, noiseMagnitude, rotationMin, rotationMax, translationMin, translationMax)

fprintf(' \n------------------------------ ');
fprintf('Copying the shapes and applying random noise and transforms...');
fprintf(' ------------------------------ \n\n');

nbCopies = 3;
pairs_array_tmp = cell(1, length(pairs_array)*nbCopies);
noiseMagnitude = 5; %5;% 0; % 0.5;
rotationMin    = -30;% -pi/4; %-pi/4;
rotationMax    = 30;% pi/4; %pi/4;
translationMin = -360;% -300; %-300;
translationMax = 360;% 300; %300;

% If a file containing transforms to apply to the shapes is provided, load it. Else, create the transforms and store them in a file.
% Transform file format: each line contains the 6 parameters of a transform (3 for translation, 3 for rotation)

transformsParameters = zeros(nbCopies*length(pairs_array), 6);
transformParametersFilePath = [pwd '\..\data\TRANSFORMS_001\'];

if exist(transformParametersFilePath, 'dir')
    transformParametersFilePath = [transformParametersFilePath 'transformsParameters.txt'];
    if exist(transformParametersFilePath, 'file')
        transformsParameters = dlmread(transformParametersFilePath);
        if(size(transformsParameters,1) < nbCopies*length(pairs_array))
            error(['The number of transforms in the file is stricly inferior to the number of copies and pairs of shapes (' num2str(nbCopies*length(pairs_array)) '). Change folder or randomly generate transforms in the code.']);
        end
    else
        transformsParameters = zeros(length(pairs_array), 6);
        for i = 1:(nbCopies*length(pairs_array))
            if i ==1 || mod(i, nbCopies) == 1
                transformsParameters(i,:) = [0, 0, 0, 0, 0, 0]; % no noise or transform is applied to have a clean copy for comparison
                continue;
            end
            transformsParameters(i,:) = [randi([translationMin, translationMax]), randi([translationMin, translationMax]), randi([translationMin, translationMax]), randi([rotationMin, rotationMax]), randi([rotationMin, rotationMax]), randi([rotationMin, rotationMax])];
        end
        dlmwrite([pwd '\transforms.txt'], transformsParameters, 'delimiter', ' ');
    end
end

for i = 1:length(pairs_array)
    display(['Creating duplicates for pair ' num2str(i) ' out of ' num2str(length(pairs_array))]);
    curPairShapes = pairs_array{i};
    for j = 1:nbCopies
        curPairShapesCopy = curPairShapes;

        % for each segment, apply a (nbSegment-1)*5 degree of rotation around the z axis before applying the random noise and global transform to the whole shape
        for nbSegment=2:size(curPairShapesCopy.segmentations_target,2)
            % Retrieve indices of current vertebra
            curIndicesToExtract = curPairShapesCopy.segmentations_target{nbSegment};

            % Compute new vertices for current vertebra
            curArray= curPairShapesCopy.shape_target.surface.VERT(curIndicesToExtract,:);
            curPointCloud = pointCloud(curArray);
            curTransform = rigidtform3d(makehgtform('zrotate', deg2rad((nbSegment-1)*5))); % rotate around z axis by (i-1)*5 degrees
            curPointCloud = pctransform(curPointCloud, curTransform);
            curArray = curPointCloud.Location;

            % Update the vertices of the current vertebra
            curPairShapesCopy.shape_target.surface.VERT(curIndicesToExtract,:) = curArray;

            % Apply same rigid transform to every component of the shape surface
            curPairShapesCopy.shape_target.surface.X(curIndicesToExtract,:) = curArray(:,1);
            curPairShapesCopy.shape_target.surface.Y(curIndicesToExtract,:) = curArray(:,2);
            curPairShapesCopy.shape_target.surface.Z(curIndicesToExtract,:) = curArray(:,3);

            % Multiply transform_target for the current segment with curTransform
            curPairShapesCopy.transform_target{nbSegment} = rigidtform3d(curTransform.A * curPairShapesCopy.transform_target{nbSegment}.A);
    
        end

        % Create a figure to display the target shape before and after the rotation (curPairShapesCopy vs curPairShapes)
        figure('Name', ['Shape ' num2str(i) ' before and after rotation'],'NumberTitle','off');
        subplot(1,2,2);
        title('After rotation');
        display_shape(curPairShapesCopy.shape_target);
        hold on;

        subplot(1,2,1);
        title('Before rotation');
        display_shape(curPairShapes.shape_target);
        hold on;
    


        % TBD: check if commented code can be deleted
        % add noise and transform the source shape %mesh, noise, translation, rotation
        %noise = noiseMagnitude * randn(size(curPairShapesCopy.shape_source.surface.VERT));
        %[curPairShapesCopy.shape_source, transformSource] = transformShape(curPairShapesCopy.shape_source, noise, transformsParameters((i-1)*nbCopies+j, 1:3), transformsParameters((i-1)*nbCopies+j, 4:6));

        % if j ==1 || mod(j, nbCopies) == 1
        %     noise = noiseMagnitude * zeros(size(curPairShapesCopy.shape_target.surface.VERT));
        %     transformTarget = rigidtform3d;
        % else
            % add noise and transform the target shape
            noise = noiseMagnitude * randn(size(curPairShapesCopy.shape_target.surface.VERT));
            [curPairShapesCopy.shape_target, transformTarget] = transformShape(curPairShapesCopy.shape_target, noise, transformsParameters((i-1)*nbCopies+j, 1:3), transformsParameters((i-1)*nbCopies+j, 4:6));

        % end
        % store the noise vectors and the transforms in the PairShapes object
        curPairShapesCopy.noise_vector_source = [];
        curPairShapesCopy.noise_vector_target = noise;

        % for each segment, store the transform
        for nbSegment=1:size(curPairShapesCopy.segmentations_source,2)
            % Update the transform for the current segment by multiplying it with the transformTarget
            curPairShapesCopy.transform_target{nbSegment} = rigidtform3d(transformTarget.A * curPairShapesCopy.transform_target{nbSegment}.A);
        end

        curPairShapesCopy.transform_source{1} = rigidtform3d;
        %curPairShapesCopy.transform_target{1} = transformTarget;

        pairs_array_tmp{(i-1)*nbCopies+j} = curPairShapesCopy;
    end
end

pairs_array = pairs_array_tmp;
clear pairs_array_tmp;
drawnow;
%% Compute and display the basis functions for each shape
if(displayBasisFunctions)
    fprintf('\n===============================')
    fprintf('Computing and displaying the basis functions...');
    fprintf('=============================== \n\n');

    CollectionDisplayBasisFunctions(pairs_array, k1, k2);
end
drawnow;
%% For each method, compute the global descriptors for each shape
fprintf(' \n=============================== ');
fprintf('Computing the global descriptors for each shape...');
fprintf(' =============================== \n\n');
pairs_array = CollectionComputeGlobalDescriptors(pairs_array, listMethods, k1, k2, numTimesGlobalDescriptors, numSkipGlobalDescriptors, displayDescriptorsGlobal);
drawnow;
%% For each method, compute and plot the local descriptors for each shape
fprintf(' \n=============================== ');
fprintf('Computing the local descriptors for each shape...');
fprintf(' =============================== \n\n');

num_skip = 15;
timesteps_lm = 100;
pairs_array = CollectionComputeLocalDescriptors(pairs_array, listMethods, k1, k2, timesteps_lm, num_skip, displayDescriptorsLocal);
drawnow;
%% Compute the functional maps using different descriptors
fprintf(' \n=============================== ');
fprintf('Computing the functional maps...');
fprintf(' =============================== \n\n');

for i = 1:length(pairs_array)
    % load the pair of shapes from the array
    curPairShapes = pairs_array{i};

    shapeSource = curPairShapes.shape_source;
    shapeTarget = curPairShapes.shape_target;
    curFolderName = curPairShapes.pair_folder_name;

    % for each descriptor combination method, compute the functional map
    for nbMethod = 1:length(listMethodsMaps)
        method = listMethodsMaps{nbMethod};
        maskMethodName = '';

        % Extract all the first elements of the cell array to find the mask method name
        firstElements = cellfun(@(x) x{1}, method, 'UniformOutput', false);

        if(any(strcmp(firstElements, 'standard')))
            maskMethodName = 'standard';
        elseif(any(strcmp(firstElements, 'slant')))
            maskMethodName = 'slant';
        elseif(any(strcmp(firstElements, 'complexResolvent')))
            maskMethodName = 'complexResolvent';
        end

        % load the descriptors from the PairShapes object
        fctTarget = []; fctSource = []; methodString = [maskMethodName ' + '];
        for j = 1:length(method)
            if(length(method{j}) == 2)
                methodString = [methodString method{j}{1} ' ' method{j}{2} ' + ' ];

                if strcmp(method{j}{2}, 'global')
                    fctTarget = [fctTarget curPairShapes.descriptors_global_target{strcmp(curPairShapes.descriptors_global_target_labels, method{j}{1})}];
                    fctSource = [fctSource curPairShapes.descriptors_global_source{strcmp(curPairShapes.descriptors_global_source_labels, method{j}{1})}];
                elseif strcmp(method{j}{2}, 'local')
                    fctTarget = [fctTarget curPairShapes.descriptors_local_target{strcmp(curPairShapes.descriptors_local_target_labels, method{j}{1})}];
                    fctSource = [fctSource curPairShapes.descriptors_local_source{strcmp(curPairShapes.descriptors_local_source_labels, method{j}{1})}];
                end


            elseif(length(method{j}) == 1)
                methodString = [methodString method{j}{1} ' + '];
            end
        end
        methodString = methodString(1:end-3); % remove the last ' + '


        BTarget = shapeTarget.evecs(:,1:k1); BSource = shapeSource.evecs(:,1:k2);
        EvTarget = shapeTarget.evals(1:k1); EvSource = shapeSource.evals(1:k2);

        fprintf(' \n------------------------------ ');
        fprintf(['Computing the functional map for vertebra number ' num2str(i) ' out of ' num2str(length(pairs_array))]);
        fprintf(' ------------------------------ \n');

        fprintf('Computing the functional map using %s descriptors...\n', methodString);
        mapFigureTitle = '';

        % If the standard Laplacian term is used, compute the standard Laplacian matrices
        if any(strcmp(maskMethodName, 'standard'))
            fprintf('Computing the functional map using %s descriptors and the standard Laplacian term...\n', methodString);
            [C_target2source, ~] = compute_fMap(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');
            mapFigureTitle = 'standard Mask';
        end

        % If the slanted Laplacian term is used, compute the slanted Laplacian matrices
        if any(strcmp(maskMethodName, 'slant'))
            fprintf('Computing the functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
            [C_target2source, ~] = compute_fMap(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');
            mapFigureTitle = 'slanted Mask';
        end

        % If the complex resolvent Laplacian term is used, compute the complex resolvent Laplacian matrices
        if any(strcmp(maskMethodName, 'complexResolvent'))
            fprintf('Computing the functional map using %s descriptors and the complex resolvent Laplacian term...\n', methodString);
            [C_target2source, ~] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'complRes');
            mapFigureTitle = 'complex resolvent Mask';
        end

        T_source2target = fMAP.fMap2pMap(BTarget,BSource,C_target2source);

        % visualize the computed maps
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString],'NumberTitle','off');
        MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target,plotOptions{:}); title(mapFigureTitle);

        % Store the mappings in the PairShapes object
        curPairShapes.mappings{nbMethod} = [T_source2target];
        curPairShapes.mappings_Labels{nbMethod} ={methodString };

        if(any(strcmp(firstElements, 'BCICP')))


            fprintf(' \n------------------------------ ');
            fprintf('Computing the reverse maps for BCICP...');
            fprintf(' ------------------------------\n');
            drawnow;

            if any(strcmp(maskMethodName, 'standard'))
                fprintf('Computing the reverse functional map using %s descriptors and the standard Laplacian term...\n', methodString);
                [C_source2target, ~] = compute_fMap(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'standard');
            end

            if any(strcmp(maskMethodName, 'slant'))
                fprintf('Computing the reverse functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
                [C_source2target, ~] = compute_fMap(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'slant');
            end

            if any(strcmp(maskMethodName, 'complexResolvent'))
                fprintf('Computing the reverse functional map using %s descriptors and the complex resolvent Laplacian term...\n', methodString);
                [C_source2target, ~] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'complRes');
            end

            T_target2source = fMAP.fMap2pMap(BSource,BTarget,C_source2target);

            % TODO: store the maps in the PairShapes object

            % refine the mapping with BCICP ([T21, T12] = bcicp_refine(S1,S2,B1,B2,T21_ini, T12_ini,num_iter))

            fprintf(' \n------------------------------ ');
            fprintf('Refining the maps with BCICP...');
            fprintf(' ------------------------------\n');
            drawnow;

            tic;
            [T_target2source_bcicp, T_source2target_new] = bcicp_refine(shapeSource, shapeTarget, BSource, BTarget, T_target2source, T_source2target,  1);
            fprintf('Time needed to compute the BCICP map: %f seconds\n', toc);

            % visualize the computed maps with BCICP
            figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and BCICP'],'NumberTitle','off');
            MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('standard Mask');

        end

        fprintf(' ------------------------------ ');
        fprintf('Refining the maps with ZoomOut...');
        fprintf(' ------------------------------\n');

        if(any(strcmp(firstElements, 'zoomOut')))
            % refine the mapping with zoomOut (final_map = refineZoomOut(initial_matches, initial_dim, S1, S2)
            tic;
            T_source2target_new = refineZoomOut(T_source2target, size(C_target2source,1), shapeSource, shapeTarget, [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' + zoomOut']);
            fprintf('Time needed to compute the zoomOut map: %f seconds\n', toc);
        end

        curPairShapes.mappings{nbMethod} = T_source2target_new;
        curPairShapes.mappings_Labels{nbMethod} = methodString;

        fprintf(' \n------------------------------ ');
        fprintf('Visualizing the maps...');
        fprintf(' ------------------------------\n');

        %Display both shapes with the drilling paths
        %plotName = ['Shapes with drilling paths - Folder ' listFolders{i} ' - FM using descriptors ' methodString ' (direct')];
        sourceTitle = ['Source shape (' curFolderName ')'];
        targetTitle = ['Target shape (' curFolderName ')'];

        % on a new figure, visualize the mapping with complex resolvent Laplacian term and the drilling paths
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (direct)'],'NumberTitle','off');
        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'direct');

        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (connex)'],'NumberTitle','off');
        display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'connex', [7 7]);


        fprintf(' ------------------------------ ');
        fprintf('Computing the rigid transforms between the shapes using all the maps obtained so far...');
        fprintf(' ------------------------------\n');

        % compute the rigid transforms between the shapes using all the maps obtained so far
        %all_Transforms is a cell array containing the transforms for each method
        %The first dimension is the method used
        %The second dimension is the transform for the whole shape or for each vertebra if a segmentation is provided
        % function used: all_Transforms = computeCorrespondances(M, N, all_matches, all_matches_names, manual_vertebra_segmentation, transformComputationMethods, partial_manual_landmarks_vertices);

        % Point cloud on partial vertebra
        %curPoints =  vertebrae_segmentation{i};

        % Transform computations
        source_segmentation = curPairShapes.segmentations_source;
        source_segmentations_labels = curPairShapes.segmentations_source_labels;

        for nbSegment=1:size(source_segmentation,2)
            % Create a temp figure to display the segment of the source shape being registered and the matching points on the target shape
            % figure('Name', [num2str(i) ' - Segment ' num2str(nbSegment) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (direct)'],'NumberTitle','off');
            % hold on;
            % %subplot(1,2,1);
            % title('Segment of source shape');
            % scatter3(shapeSource.surface.VERT(source_segmentation{nbSegment},1), shapeSource.surface.VERT(source_segmentation{nbSegment},2), shapeSource.surface.VERT(source_segmentation{nbSegment},3), 5, 'filled');
            % hold on;
            % 
            % % Display the matching points on the target shape for the current segment
            % %subplot(1,2,2);
            % %title('Segment of target shape');
            % matchingIndices = T_source2target_new(source_segmentation{nbSegment});
            % scatter3(shapeTarget.surface.VERT(matchingIndices,1), shapeTarget.surface.VERT(matchingIndices,2), shapeTarget.surface.VERT(matchingIndices,3), 5, 'filled');

            curPointsIndices =  source_segmentation{nbSegment};
            disp(['==== Computing rigid transform for vertebra #', num2str(nbSegment)]);
            [source_T_target, ~, ~] = computeTransformBetweenShapes(shapeTarget,shapeSource, curPointsIndices, T_source2target_new, 'ransac');
            disp('====  Rigid transform computed');
            curPairShapes.registrations_source_to_target{nbMethod}{nbSegment} = source_T_target;
        end


        % PLOT REGISTERED SHAPES
        % This displays the vertebrae of the source shape  mapped to the target shape
        % The source shape is displayed as a black mesh and the target shape is displayed as a full color mesh
        % The registration is done using the vertebrae segmentation and breakpoints can be added to the display loop
        % in the function plotTwoShapes
        figure('Name', [num2str(i) ' - Folder ' curFolderName ' - Moving source to target'],'NumberTitle','off');
        hold on;
        title('Shape source moved to shape target !');

        manual_vertebra_segmentation_colors_per_vertex = {ones(shapeSource.nv, 3)};
        transformToPlot = curPairShapes.registrations_source_to_target{nbMethod};
        plotTwoShapes(shapeTarget,shapeSource, curPairShapes.segmentations_source, transformToPlot, manual_vertebra_segmentation_colors_per_vertex);


        trajOnSource = curPairShapes.trajectories_source;
        indexEntryPointsSource = trajOnSource(:,1);
        indexEndPointsSource = trajOnSource(:,2);

        entryPointsSource = shapeSource.surface.VERT(indexEntryPointsSource', :);
        endPointsSource = shapeSource.surface.VERT(indexEndPointsSource', :);

        trajColors = hsv(size(entryPointsSource, 1));

        % For each drilling path, identify the corresponding segmentation and apply the corresponding transform
        for curTrajNb = 1:size(entryPointsSource, 1)

            % determine to which segment the current trajectory belongs, i.e. what is the index of the segment containing the entry and end points
            for nbSegment=1:size( curPairShapes.segmentations_source,2)
                if any(ismember( curPairShapes.segmentations_source{nbSegment}, indexEntryPointsSource(curTrajNb))) && any(ismember( curPairShapes.segmentations_source{nbSegment}, indexEndPointsSource(curTrajNb)))
                    break;
                end
            end

            curTransform = transformToPlot{1,nbSegment};%TBD: genericity for multiple vertebrae!

            curEntryPointsSource = entryPointsSource(curTrajNb, :);
            curEndPointsSource = endPointsSource(curTrajNb, :);

            pointCloudentryPointsSource = pointCloud(curEntryPointsSource);
            pointCloudEndPointsSource = pointCloud(curEndPointsSource);

            %Compute corresponding points using the transform
            entryPointsTargetVertices = pctransform(pointCloudentryPointsSource, curTransform); %Here we invert the transform because we want to go from M to N
            endPointsTargetVertices = pctransform(pointCloudEndPointsSource, curTransform);

            entryPointsTargetVertices = entryPointsTargetVertices.Location;
            endPointsTargetVertices = endPointsTargetVertices.Location;

            curColor = trajColors(curTrajNb, :);
            [mainplotTarget, entryExtensionPlotTarget, endExtensionPlotTarget] = plotTrajectory(entryPointsTargetVertices, endPointsTargetVertices, curColor, 5, 1);
        end



        fprintf(' ------------------------------ ');
        fprintf('Exporting the maps...');
        fprintf(' ------------------------------\n');

        % export the maps to text files for later use
        map_name = [maps_dir num2str(i) ' - map_' curFolderName '_' methodString '_' maskMethodName '.txt'];
    end

    pairs_array{i} = curPairShapes;

    registration_errors = curPairShapes.computeRegistrationErrors();
    average_errors = curPairShapes.computeAverageRegistrationErrors(registration_errors);
    mean_errors_array{end+1} = average_errors;

    % Update the final error figure by plotting all the errors stored in the mean_errors_array
    figure(figErrors);

    % Plot the errors
    % Plot error of X component for each pair of shapes and each method

    % Error of X component for each pair of shapes for method 1
    subplot(2,3,1);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end

    hold on;
    title('X component error');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');

    % Plot error of Y component for each pair of shapes and each method
    subplot(2,3,2);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end
    title('Y component error');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');

    % Plot error of Z component for each pair of shapes and each method
    subplot(2,3,3);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+2*length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end
    hold on;
    title('Z component error');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');

    % Plot error of rotation R for each pair of shapes and each method
    subplot(2,3,4);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+3*length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end
    hold on;
    title('Rotation error (R)');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');

    % Plot error of rotation P for each pair of shapes and each method

    % Error of translation for each pair of shapes for method 1
    subplot(2,3,5);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+4*length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end
    hold on;
    title('Rotation error (P)');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');

    % Plot error of rotation Yaw for each pair of shapes and each method

    % Error of scaling for each pair of shapes for method 1
    subplot(2,3,6);
    hold off;
    for nbMethod = 1:length(listMethodsMaps)
        plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+5*length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
        hold on;
    end
    hold on;
    title('Rotation error (Y)');
    xlabel('Pair of shapes');
    ylabel('Error');
    legend('show');
end


%% Export all figures
saveAllFigures(destinationFolder, 'png');

%% Save all data
save("results\all_data.mat");
%% Check memory
memory