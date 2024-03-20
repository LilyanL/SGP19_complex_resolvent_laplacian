%% Let's setup everything nicely
clc; clear all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode

% create folders to export results and figures
mkdir results;
destinationFolder_root = '.\results\';

% Create a subfolder to store the maps
maps_dir = [destinationFolder_root 'maps\'];
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
listFolders = {'PAIR_009'};
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
    %{{'WKS', 'global'}, {'HKS', 'local'}, {'slant'}, {'BCICP'}, {'zoomOut'}};
    %{{'WKS', 'global'}, {'HKS', 'local'}, {'complexResolvent'}, {'zoomOut'}};
    %{{'WKS', 'global'}, {'HKS', 'local'}, {'standard'}, {'BCICP'}, {'zoomOut'}};
    };

%% Display options
displayShapePairs = true;
displayShapePairsWithPaths = true;
displayBasisFunctions = false;
displayDescriptorsGlobal = false;
displayDescriptorsLocal = false;

%% Data used to display the final results
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

drawnow;
%% Duplicate the pairs of shapes to apply random noise and transforms to the copies and store them in a temporary cell array
% For each pair of shapes, create nbCopies copies with random noise and transforms using function transformShape (limits are provided as parameters)
% transformShape(mesh, noiseMagnitude, rotationMin, rotationMax, translationMin, translationMax)

fprintf(' \n============================== ');
fprintf('Copying the shapes and applying random noise and transforms...');
fprintf(' ============================== \n');

nbCopies = 5;
noiseMagnitudeVec = [0.5 ]; %5;% 0; % 0.5;
angleTorsionVertebra = 25;
rotationMin    = -30;% -pi/4; %-pi/4;
rotationMax    = 30;% pi/4; %pi/4;
translationMin = -360;% -300; %-300;
translationMax = 360;% 300; %300;

centerData = true;
centerDataOption = 'partial';

 extensionLengthRatio = 0.2;
 trajectoryWidthMillimeters = 5;
 trajectoryWidth = trajectoryWidthMillimeters/0.352806/2; %LineWidth: 1 point is 0,352806 mm so 6 mm is 17 points

for curNoiseValue = noiseMagnitudeVec
    noiseMagnitude = curNoiseValue;
    destinationFolder = [destinationFolder_root 'noise_' num2str(noiseMagnitude) '\'];

    close all;
    figErrors = figure('Name', 'Errors', 'NumberTitle', 'off');
    figErrorsMAE = figure('Name', 'Mean Absolute Errors', 'NumberTitle', 'off');
    figErrorsRMSE = figure('Name', 'Root Mean Squared Errors', 'NumberTitle', 'off');

    if ~exist(destinationFolder, 'dir')
        mkdir(destinationFolder);
    end

    pairs_array = CollectionLoadShapes(foldersPaths, meshOptions, displayShapePairs, displayShapePairsWithPaths, loadOption, centerData, centerDataOption);
    pairs_array_tmp = cell(1, length(pairs_array)*nbCopies);


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

    firstCaseDisplayed = false;
    for i = 1:length(pairs_array)
        display(['Creating duplicates for pair ' num2str(i) ' out of ' num2str(length(pairs_array))]);
        
        curPairShapes = pairs_array{i};
        for j = 1:nbCopies
            numCurrentFig = 1;
            curPairShapesCopy = curPairShapes;

            % for each segment, apply a (nbSegment-1)*5 degree of rotation around the z axis before applying the random noise and global transform to the whole shape
            for nbSegment=2:size(curPairShapesCopy.segmentations_target,2)
                % Retrieve indices of current vertebra
                curIndicesToExtract = curPairShapesCopy.segmentations_target{nbSegment};

                % Compute new vertices for current vertebra
                curArray= curPairShapesCopy.shape_target.surface.VERT(curIndicesToExtract,:);
                curPointCloud = pointCloud(curArray);
                curTransform = rigidtform3d(makehgtform('zrotate', deg2rad((nbSegment-1)*angleTorsionVertebra))); % rotate around z axis by (i-1)*5 degrees
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

            if (~firstCaseDisplayed)
                % Create a figure to display the target shape before and after the rotation (curPairShapesCopy vs curPairShapes)
                figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Shapes before (color) and after (black) rotation'],'NumberTitle','off');
                numCurrentFig = numCurrentFig+1;
                title('Before and after rotation same figure');
                display_shape(curPairShapesCopy.shape_target);
                hold on;
                h = display_shape(curPairShapes.shape_target);
                set(h, 'edgecolor', 'black');
                set(h, 'FaceColor', 'none');
                hold on;

                firstCaseDisplayed = true;
            end

            noise = noiseMagnitude * randn(size(curPairShapesCopy.shape_target.surface.VERT));
            [curPairShapesCopy.shape_target, transformTarget] = transformShape(curPairShapesCopy.shape_target, noise, transformsParameters((i-1)*nbCopies+j, 1:3), transformsParameters((i-1)*nbCopies+j, 4:6));

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
        fprintf('=============================== \n');

        CollectionDisplayBasisFunctions(pairs_array, k1, k2);
    end
    drawnow;
    %% For each method, compute the global descriptors for each shape
    fprintf(' \n=============================== ');
    fprintf('Computing the global descriptors for each shape...');
    fprintf(' =============================== \n');
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

        % Display the two shapes on the same plot in a new figure
        figure('Name', [num2str(i) '.' num2str(numCurrentFig) '  - Folder ' curFolderName ' - Shapes on same plot' ],'NumberTitle','off');
        numCurrentFig = numCurrentFig + 1;
        display_shape(shapeSource);
        hold on;
        display_shape(shapeTarget);
        title(['Shapes ' num2str(i) ' - ' curFolderName]);



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
            fprintf(['Computing the functional map for pair of shapes #' num2str(i) ' out of ' num2str(length(pairs_array))]);
            fprintf(' ------------------------------ \n');

            fprintf('Computing the functional map using %s descriptors...\n', methodString);
            mapFigureTitle = '';

            % If the standard Laplacian term is used, compute the standard Laplacian matrices
            if any(strcmp(maskMethodName, 'standard'))
                fprintf('Computing the functional map using %s descriptors and the standard Laplacian term...\n', methodString);
                [C_target2source, ~] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'standard');
                mapFigureTitle = 'standard Mask';
            end

            % If the slanted Laplacian term is used, compute the slanted Laplacian matrices
            if any(strcmp(maskMethodName, 'slant'))
                fprintf('Computing the functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
                [C_target2source, ~] = compute_fMap_complRes(shapeTarget,shapeSource,BTarget,BSource,EvTarget,EvSource,fctTarget,fctSource,para, 'slant');
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
            figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - FM using descriptors ' methodString],'NumberTitle','off');
            numCurrentFig = numCurrentFig +1;
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
                    [C_source2target, ~] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'standard');
                end

                if any(strcmp(maskMethodName, 'slant'))
                    fprintf('Computing the reverse functional map using %s descriptors and the slanted Laplacian term...\n', methodString);
                    [C_source2target, ~] = compute_fMap_complRes(shapeSource,shapeTarget,BSource,BTarget,EvSource,EvTarget,fctSource,fctTarget,para, 'slant');
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
                [T_target2source_bcicp, T_source2target_new] = bcicp_refine(shapeSource, shapeTarget, BSource, BTarget, T_target2source, T_source2target,  5);
                T_source2target_bcicp = T_source2target_new;
                fprintf('Time needed to compute the BCICP map: %f seconds\n', toc);

                % visualize the computed maps with BCICP
                figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and BCICP'],'NumberTitle','off');
                numCurrentFig = numCurrentFig + 1;

                MESH.PLOT.visualize_map_colors(shapeSource,shapeTarget,T_source2target_new,plotOptions{:}); title('After BCICP');
            end

            fprintf(' ------------------------------ ');
            fprintf('Refining the maps with ZoomOut...');
            fprintf(' ------------------------------\n');

            if(any(strcmp(firstElements, 'zoomOut')))
                % refine the mapping with zoomOut (final_map = refineZoomOut(initial_matches, initial_dim, S1, S2)
                tic;
                T_source2target_new = refineZoomOut(T_source2target_new, size(C_target2source,1), shapeSource, shapeTarget, [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' + zoomOut']);
                numCurrentFig = numCurrentFig + 1;
                fprintf('Time needed to compute the zoomOut map: %f seconds\n', toc);
            end
    
            %% 
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
            figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (direct)'],'NumberTitle','off');
            numCurrentFig = numCurrentFig + 1;
            display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'direct');

            figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - FM using descriptors ' methodString ' and drilling paths (connex)'],'NumberTitle','off');
            numCurrentFig = numCurrentFig + 1;
            display_pair_shapes_and_paths(shapeSource, shapeTarget, sourceTitle, targetTitle, curPairShapes.trajectories_source, T_source2target_new, 'connex', [7 7]);


            %%
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

                % If segmentation_source_upper_part is not empty, use it to keep only indices corresponding to the upper part of the vertebra
                if ~isempty(curPairShapes.segmentation_source_upper_part)
                    curPointsIndices = intersect(curPointsIndices, curPairShapes.segmentation_source_upper_part);
                    % display the percentage of points kept
                    disp(['Percentage of points kept for vertebra #', num2str(nbSegment), ' after upper part segmentation: ', num2str(length(curPointsIndices)/length(source_segmentation{nbSegment})*100), '%']);
                end

                disp(['==== Computing rigid transform for vertebra #', num2str(nbSegment)]);
                [source_T_target, ~, ~] = computeTransformBetweenShapes(shapeTarget,shapeSource, curPointsIndices, T_source2target_new, 'ransac');
                disp('====  Rigid transform computed');
                curPairShapes.registrations_source_to_target{nbMethod}{nbSegment} = source_T_target;
            end

            %% PLOT REGISTERED SHAPES
            % This displays the vertebrae of the source shape  mapped to the target shape
            % The source shape is displayed as a black mesh and the target shape is displayed as a full color mesh
            % The registration is done using the vertebrae segmentation and breakpoints can be added to the display loop
            % in the function plotTwoShapes
            figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - Moving source to target'],'NumberTitle','off');
            numCurrentFig = numCurrentFig + 1;
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
                [mainplotTarget, entryExtensionPlotTarget, endExtensionPlotTarget] = plotTrajectory(entryPointsTargetVertices, endPointsTargetVertices, curColor, trajectoryWidth, extensionLengthRatio);
            end



            
        end

        pairs_array{i} = curPairShapes;

        registration_errors = curPairShapes.computeRegistrationErrors();
        average_errors = curPairShapes.computeAverageRegistrationErrors(registration_errors);
        mean_errors_array{end+1} = average_errors;

        %% Display the errors on each component of the registration for each pair of shapes and each method

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
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('X component error');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        % Plot error of Y component for each pair of shapes and each method
        subplot(2,3,2);
        hold off;
        for nbMethod = 1:length(listMethodsMaps)
            plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
            hold on;
        end

        hold on;
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('Y component error');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        % Plot error of Z component for each pair of shapes and each method
        subplot(2,3,3);
        hold off;
        for nbMethod = 1:length(listMethodsMaps)
            plot(1:length(mean_errors_array), cellfun(@(x) x(nbMethod+2*length(listMethodsMaps)), mean_errors_array), '-o','DisplayName', ['Method ' num2str(nbMethod)]);
            hold on;
        end
        hold on;
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('Z component error');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        % Plot error of rotation R for each pair of shapes and each method
        subplot(2,3,4);
        hold off;
        for nbMethod = 1:length(listMethodsMaps)
            y = cellfun(@(x) x(nbMethod+3*length(listMethodsMaps)), mean_errors_array);
            y = rad2deg(y);
            plot(1:length(mean_errors_array), y, '-o','DisplayName', ['Method ' num2str(nbMethod)]);
            hold on;
        end
        hold on;
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('Rotation error (Rx)');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        % Plot error of rotation P for each pair of shapes and each method

        % Error of translation for each pair of shapes for method 1
        subplot(2,3,5);
        hold off;
        for nbMethod = 1:length(listMethodsMaps)
            y = cellfun(@(x) x(nbMethod+4*length(listMethodsMaps)), mean_errors_array);
            y = rad2deg(y);
            plot(1:length(mean_errors_array), y, '-o','DisplayName', ['Method ' num2str(nbMethod)]);
            hold on;
        end
        hold on;
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('Rotation error (Ry)');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        % Plot error of rotation Yaw for each pair of shapes and each method

        % Error of scaling for each pair of shapes for method 1
        subplot(2,3,6);
        hold off;
        for nbMethod = 1:length(listMethodsMaps)
            y = cellfun(@(x) x(nbMethod+5*length(listMethodsMaps)), mean_errors_array);
            y = rad2deg(y);
            plot(1:length(mean_errors_array), y, '-o','DisplayName', ['Method ' num2str(nbMethod)]);
            hold on;
        end
        hold on;
        yline(0, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '0');
        title('Rotation error (Rz)');
        xlabel('Pair of shapes #');
        ylabel('Error (mm)');
        legend('show');

        drawnow

    

    %% On a new figure, take the first pair of shapes and display the following figures:
    % Figure 1: the source shape with the drilling paths. All the drilling paths are displayed in the same color (blue)
    % Figure 2: the target shape with the drilling paths. The target shape is transformed by applied the inverse of the ground truth transform to bring it back to the source shape.
    %            All the drilling paths are displayed in the same color (red)
    % Figure 3: the target shape with both the ground truth drilling paths and the computed drilling paths. The ground truth drilling paths are displayed in blue and the computed drilling paths are displayed in red.
    %            Again, the target shape is transformed by applied the inverse of the ground truth transform to bring it back to the source shape.
    % Overall, this allows to visualize with the same camera view the source shape, the target shape and the drilling paths on both shapes and to easily compare the ground truth drilling paths and the computed drilling paths.

   
    groundTruthTrajColor = [0 0 1];
    computedTrajColor = [1 0 0];

    % Take the first pair of shapes
    curPairShapes = pairs_array{i};

    % ==============================
    % ==============================
    fig1 = figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - Drilling paths comparison 1 - Plan on source shape'],'NumberTitle','off');
    % Display the source shape with the drilling paths
  
    view(0,-50);
    hold on;
    title('Plan on source shape');
    display_shape(curPairShapes.shape_source);
    hold on;
    % Display the drilling paths on the source shape
    
    for curTrajNb = 1:size(curPairShapes.trajectories_source, 1)
        curTraj = curPairShapes.trajectories_source(curTrajNb, :);
        curEntry = curPairShapes.shape_source.surface.VERT(curTraj(1), :);
        curEnd = curPairShapes.shape_source.surface.VERT(curTraj(2), :);

        %plotTrajectory(entryPoint,endPoint, color, radius, extensionLengthRatio, parentFigure)
        
        plotTrajectory(curEntry, curEnd, groundTruthTrajColor, trajectoryWidth, extensionLengthRatio); %LineWidth: 1 point is 0,352806 mm so 6 mm is 17 points
        %plot3([curEntry(1) curEnd(1)], [curEntry(2) curEnd(2)], [curEntry(3) curEnd(3)], 'Color', sourceTrajColor, 'LineWidth', 2);
    end

    % ==============================
    % ==============================
    % Display the target shape with the drilling paths
    % Transform the target shape by applying the inverse of the ground truth transform to bring it back to the source shape
    % Transform the drilling paths by the inverse of the ground truth transform
    %fig2 = figure('Name', 'Drilling paths comparison 2 - Plan on target shape','NumberTitle','off');
    fig2 = figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - Drilling paths comparison 2 - Plan on target shape'],'NumberTitle','off');
    view(10,-15);
    hold on;
    
    targetShape = curPairShapes.shape_target;
    display_shape(targetShape);
    title('Target shape with perfectly transferred drilling paths');

    % Display the drilling paths on the target shape
    

    for curTrajNb = 1:size(curPairShapes.trajectories_source, 1)
        curTraj = curPairShapes.trajectories_source(curTrajNb, :);

        % Identify to which segment the current trajectory belongs
        for nbSegment=1:size(curPairShapes.segmentations_source,2)
            if any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(1))) && any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(2)))
                break;
            end
        end

        % Retrieve the transform for the current segment
        curTransform = curPairShapes.registrations_source_to_target{1}{nbSegment};
        curEntry = curPairShapes.shape_source.surface.VERT(curTraj(1), :);
        curEnd = curPairShapes.shape_source.surface.VERT(curTraj(2), :);

        % Apply the transform to the source drilling path to bring it to the target shape
        curEntry = pctransform(pointCloud(curEntry), curTransform).Location;
        curEnd = pctransform(pointCloud(curEnd), curTransform).Location;

        plotTrajectory(curEntry, curEnd, groundTruthTrajColor, trajectoryWidth, extensionLengthRatio); %LineWidth: 1 point is 0,352806 mm so 6 mm is 17 points
        %plot3([curEntry(1) curEnd(1)], [curEntry(2) curEnd(2)], [curEntry(3) curEnd(3)], 'Color', groundTruthTrajColor, 'LineWidth', 2);
    end

    % ==============================
    % ==============================
    % Display the target shape with both the ground truth drilling paths and the computed drilling paths
    % Transform the target shape by applying the inverse of the ground truth transform to bring it back to the source shape
    % Transform the drilling paths by the inverse of the ground truth transform
    %fig3 = figure('Name', 'Drilling paths comparison 3 - Transfer to target shape','NumberTitle','off');
    fig3 = figure('Name', [num2str(i) '.' num2str(numCurrentFig) ' - Folder ' curFolderName ' - Drilling paths comparison 3 - Transfer to target shape'],'NumberTitle','off');
    view(10,-15);
    hold on;
    title('Target shape with both ground truth and computed drilling paths');
    display_shape(targetShape);
    hold on;

    % Display the ground truth drilling paths on the target shape
    for curTrajNb = 1:size(curPairShapes.trajectories_source, 1)
        curTraj = curPairShapes.trajectories_source(curTrajNb, :);

        % Identify to which segment the current trajectory belongs
        for nbSegment=1:size(curPairShapes.segmentations_source,2)
            if any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(1))) && any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(2)))
                break;
            end
        end

        % Compute the transform for the current segment
        curTransform = curPairShapes.transform_target{nbSegment};

        curEntry = curPairShapes.shape_source.surface.VERT(curTraj(1), :);
        curEnd = curPairShapes.shape_source.surface.VERT(curTraj(2), :);

        % Apply the transform to the source drilling path to bring it to the target shape
        curEntry = pctransform(pointCloud(curEntry), curTransform).Location;
        curEnd = pctransform(pointCloud(curEnd), curTransform).Location;

        plotTrajectory(curEntry, curEnd, groundTruthTrajColor, trajectoryWidth, extensionLengthRatio); %LineWidth: 1 point is 0,352806 mm so 6 mm is 17 points
        %plot3([curEntry(1) curEnd(1)], [curEntry(2) curEnd(2)], [curEntry(3) curEnd(3)], 'Color', groundTruthTrajColor, 'LineWidth', 2);
    end

    % Display the computed drilling paths on the target shape
    for curTrajNb = 1:size(curPairShapes.trajectories_source, 1)
        curTraj = curPairShapes.trajectories_source(curTrajNb, :);

        % Find to which segment the current trajectory belongs
        for nbSegment=1:size(curPairShapes.segmentations_source,2)
            if any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(1))) && any(ismember(curPairShapes.segmentations_source{nbSegment}, curTraj(2)))
                break;
            end
        end

        % Get the computed transform for the current segment
        computedTransform = curPairShapes.registrations_source_to_target{1}{nbSegment};%TBD: fix this! Poser les calculs, l'ordre de multiplcation des matrices pose problème
        inverseComputedTransform = invert(computedTransform);

        % Apply the ground truth transform to the source drilling path to bring it to the target shape
        curEntry = curPairShapes.shape_source.surface.VERT(curTraj(1), :);
        curEnd = curPairShapes.shape_source.surface.VERT(curTraj(2), :);

        %curEntry = pctransform(pointCloud(curEntry), groundTruthTransform).Location;
        %curEnd = pctransform(pointCloud(curEnd), groundTruthTransform).Location;

        % Apply the computed transform to the source drilling path to bring it back to the source shape
        curEntry = pctransform(pointCloud(curEntry), computedTransform).Location;
        curEnd = pctransform(pointCloud(curEnd), computedTransform).Location;

        % Compute the transform for the current segment
        %curRotation = rigidtform3d(makehgtform('zrotate', deg2rad((nbSegment-1)*angleTorsionVertebra))); % rotate around z axis by (i-1)*5 degrees

        %curEntry = pctransform(pointCloud(curEntry), invert(curRotation)).Location;
        %curEnd = pctransform(pointCloud(curEnd), invert(curRotation)).Location;

        % Plot the transformed drilling path
        plotTrajectory(curEntry, curEnd, computedTrajColor, trajectoryWidth, extensionLengthRatio); %LineWidth: 1 point is 0,352806 mm so 6 mm is 17 points
        % plot3([curEntry(1) curEnd(1)], [curEntry(2) curEnd(2)], [curEntry(3) curEnd(3)], 'Color', computedTrajColor, 'LineWidth', 2);
    end

    %Export the 3 figures to a single folder named 'drilling_paths_comparison' with the noise value combined in the name
    folderName = ['.\results\drilling_paths_comparison_noise_' num2str(curNoiseValue)];
    mkdir(folderName);
    figureList = [fig1 fig2 fig3];
    for iFig = 1:length(figureList)
        FigHandle = figureList(iFig);
        FigName   = [FigHandle.Name];
        set(0, 'CurrentFigure', FigHandle);
    
        saveas(FigHandle, fullfile(folderName,  [FigName '.png']));
        savefig(fullfile(folderName, [FigName '.fig']));
    end

    end

    %% On new figures, compute the mean absolute errors, the root mean squared errors for each vertebra and each method
    % More precisely, for each method, compute the mean absolute error, the root mean squared error and the R2 score for each vertebra
    % over all the copies of the shapes and store the results as a matrix

    % Initialize the matrices to store the errors
    errorsRMSETranslation = zeros(length(listMethodsMaps), size(curPairShapesCopy.segmentations_source,2));
    errorsMAETranslation = zeros(length(listMethodsMaps), size(curPairShapesCopy.segmentations_source,2));
    errorsRMSERotation = zeros(length(listMethodsMaps), size(curPairShapesCopy.segmentations_source,2));
    errorsMAERotation = zeros(length(listMethodsMaps), size(curPairShapesCopy.segmentations_source,2));

    % Loop over all the pairs of shapes to compute the errors
    %for nbSegment=1:size(curPairShapesCopy.segmentations_source,2)
    for nbShapes = 1:length(pairs_array)
        curPairShapes = pairs_array{nbShapes};
        registration_errors = curPairShapes.computeRegistrationErrors();
        for nbMethod = 1:length(listMethodsMaps)
            for nbSegment=1:size(curPairShapes.segmentations_source,2)
                % Compute the errors for the current pair of shapes and the current method and the current segment

                %curErrorRMSE is the root mean squared error for the current pair of shapes, the current method and the current segment.
                %It must be explicitly computed
                curError = registration_errors{nbMethod, nbSegment};
                %curErrorRMSETranslation = curError.

                errorsRMSETranslation(nbMethod, nbSegment) = errorsRMSETranslation(nbMethod, nbSegment) + (curError.error_X^2 + curError.error_Y^2 + curError.error_Z^2);
                errorsMAETranslation(nbMethod, nbSegment) = errorsMAETranslation(nbMethod, nbSegment) + (abs(curError.error_X) + abs(curError.error_Y) + abs(curError.error_Z));
                errorsRMSERotation(nbMethod, nbSegment) = errorsRMSERotation(nbMethod, nbSegment) + (curError.error_R^2 + curError.error_P^2 + curError.error_Yaw^2);
                errorsMAERotation(nbMethod, nbSegment) = errorsMAERotation(nbMethod, nbSegment) + (abs(curError.error_R) + abs(curError.error_P) + abs(curError.error_Yaw));
            end
        end
    end

    % Compute the mean errors for each method and each segment
    errorsRMSETranslation = errorsRMSETranslation / length(pairs_array);
    errorsMAETranslation = errorsMAETranslation / length(pairs_array);
    errorsRMSERotation = errorsRMSERotation / length(pairs_array);
    errorsMAERotation = errorsMAERotation / length(pairs_array);

    % Apply square root to the mean squared errors to obtain the root mean squared errors
    errorsRMSETranslation = sqrt(errorsRMSETranslation);
    errorsRMSERotation = sqrt(errorsRMSERotation);

    % Plot the mean absolute errors, the root mean squared errors for each vertebra and each method
    x_segments = 1:size(curPairShapesCopy.segmentations_source,2);
    figure(figErrorsMAE);
    subplot(2,1,1);
    for curMethod = 1:length(listMethodsMaps)
        plot(x_segments, errorsMAETranslation(curMethod,:), '-o', 'DisplayName', ['Method ' num2str(curMethod)]);
        hold on;
    end

    title('Mean Absolute Errors - Translation');
    xlabel('Vertebra #');
    ylabel('Error (mm)');
    legend('show');

    subplot(2,1,2);
    for curMethod = 1:length(listMethodsMaps)
        plot(x_segments, rad2deg(errorsMAERotation(curMethod,:)), '-o', 'DisplayName', ['Method ' num2str(curMethod)]);
        hold on;
    end

    title('Mean Absolute Errors - Rotation');
    xlabel('Vertebra #');
    ylabel('Error (degrees)');
    legend('show');

    figure(figErrorsRMSE);
    subplot(2,1,1);
    for curMethod = 1:length(listMethodsMaps)
        plot(x_segments, errorsRMSETranslation(curMethod,:), '-o', 'DisplayName', ['Method ' num2str(curMethod)]);
        hold on;
    end

    title('Root Mean Squared Errors - Translation');
    xlabel('Vertebra #');
    ylabel('Error (mm)');
    legend('show');

    subplot(2,1,2);
    for curMethod = 1:length(listMethodsMaps)
        plot(x_segments, rad2deg(errorsRMSERotation(curMethod,:)), '-o', 'DisplayName', ['Method ' num2str(curMethod)]);
        hold on;
    end

    title('Root Mean Squared Errors - Rotation');
    xlabel('Vertebra #');
    ylabel('Error (degrees)');
    legend('show');

    drawnow;

    %% Open a CSV file to store the errors for the current noise value
    csvFileName = ['results\errors_noise_' num2str(curNoiseValue) '.csv'];
    fid = fopen(csvFileName, 'w');
    % Open a CSV file to store the errors for the current noise value
    csvFileName = ['results\errors_noise_' num2str(curNoiseValue) '.csv'];
    fid = fopen(csvFileName, 'w');

    % Write separator
    fprintf(fid, 'sep=,\n');

    % Write the header row
    header = 'Method,Segment,Error_RMSE_Translation,Error_MAE_Translation,Error_RMSE_Rotation,Error_MAE_Rotation\n';
    fprintf(fid, header);

    % Write the errors for each method and each segment
    for nbMethod = 1:length(listMethodsMaps)
        for nbSegment = 1:size(curPairShapesCopy.segmentations_source,2)
            % Get the errors for the current method and segment
            error_RMSE_Translation = errorsRMSETranslation(nbMethod, nbSegment);
            error_MAE_Translation = errorsMAETranslation(nbMethod, nbSegment);
            error_RMSE_Rotation = rad2deg(errorsRMSERotation(nbMethod, nbSegment));
            error_MAE_Rotation = rad2deg(errorsMAERotation(nbMethod, nbSegment));

            % Write the errors to the CSV file
            row = [num2str(nbMethod), ',', num2str(nbSegment), ',', num2str(error_RMSE_Translation), ',', num2str(error_MAE_Translation), ',', num2str(error_RMSE_Rotation), ',', num2str(error_MAE_Rotation), '\n'];
            fprintf(fid, row);
        end
    end

    % Close the CSV file
    fclose(fid);
    fclose("all");  %dirty


    %% Export all figures
    saveAllFigures(destinationFolder, 'png');
    saveAllFigures(destinationFolder, 'fig');

    %% Save all data
    save(['results\all_data_noise_' num2str(curNoiseValue) '.mat']);

    end
%% Check memory
memory