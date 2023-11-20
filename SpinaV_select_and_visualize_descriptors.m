%% Let's setup everything nicely
clc; clear all, close all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode


%% some parameters
k1 = 100;
% params to pre-process the meshes
meshOptions = {'IfComputeGeoDist',false,'IfComputeLB',true,'IfComputeNormals',true,'numEigs',100};

% params to visualize the maps
plotOptions = {'IfShowCoverage',false,'OverlayAxis','y','cameraPos',[0,90]};

%% Select and load the source meshe using the GUI
isLoadedSource = false;
[file,path] = uigetfile('*.off','Select the source shape to load');
if isequal(file,0)
    disp('User selected Cancel');
    return;
else
    isLoadedSource = true;
    disp(['User selected ', fullfile(path,file)]);
end

% load the mesh
shapeSource = MESH.MESH_IO.read_shape([path, file]);

% center data
shapeSource.surface.VERT = shapeSource.surface.VERT - mean(shapeSource.surface.VERT);

% preprocess the mesh
shapeSource = MESH.preprocess(shapeSource,meshOptions{:});

% compute eigenvecs and eigenvalues
B = shapeSource.evecs(:,1:k1); 
Ev = shapeSource.evals(1:k1);

%% Select a 3D point on the mesh using the GUI and compute the descriptors
% select a point on the mesh
selectionFigure = figure('Name', 'Select the entry point');
WKSFigure = figure(2);
HKSFigure = figure(3);

figure(1);
title('Select the entry point');
subplot(2,2, [1 3]);
h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), 'FaceColor', 'interp', 'EdgeColor', 'none');
alpha(0.5);
hold on;

pointCloud = shapeSource.surface.VERT';
% set the callback, pass pointCloud to the callback function
%set(selectionFigure, 'WindowButtonDownFcn', {@callbackClickA3DPoint, pointCloud}); 

% wait for user to close the figure
while (ishandle(selectionFigure) || ishandle(WKSFigure)||ishandle(HKSFigure))

    % wait until the user has clicked on the mesh
    while waitforbuttonpress ~= 0
        pause 1;
    end

    point = get(gca, 'CurrentPoint'); % mouse click position
    camPos = get(gca, 'CameraPosition'); % camera position
    camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to

    camDir = camPos - camTgt; % camera direction
    camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector

    % build an orthonormal frame based on the viewing direction and the 
    % up vector (the "view frame")
    zAxis = camDir/norm(camDir);    
    upAxis = camUpVect/norm(camUpVect); 
    xAxis = cross(upAxis, zAxis);
    yAxis = cross(zAxis, xAxis);

    rot = [xAxis; yAxis; zAxis]; % view rotation 

    % the point cloud represented in the view frame
    rotatedPointCloud = rot * pointCloud; 

    % the clicked point represented in the view frame
    rotatedPointFront = rot * point' ;

    % find the nearest neighbour to the clicked point 
    pointCloudIndex = dsearchn(rotatedPointCloud(1:2,:)', ... 
        rotatedPointFront(1:2));

    h = findobj(gca,'Tag','pt'); % try to find the old point
    selectedPoint = pointCloud(:, pointCloudIndex); 

    if isempty(h) % if it's the first click (i.e. no previous point to delete)
        
        % highlight the selected point
        h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
            selectedPoint(3,:), 'r.', 'MarkerSize', 20); 
        set(h,'Tag','pt'); % set its Tag property for later use   

    else % if it is not the first click

        delete(h); % delete the previously selected point
        
        % highlight the newly selected point
        h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
            selectedPoint(3,:), 'r.', 'MarkerSize', 20);  
        set(h,'Tag','pt');  % set its Tag property for later use

    end

    fprintf('you clicked on point number %d\n', pointCloudIndex);

    %% Compute an visualize the descriptors
    lm_idx = pointCloudIndex;

    % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
    timesteps_lm = 100;

    %compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
    numEigs = 100;
    method = 'WKS';
    lm_fct_WKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx,timesteps_lm,1, method);
    method = 'HKS';
    lm_fct_HKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx,timesteps_lm,1, method);

    % ignore some of the descriptors
    skip_lm_step = 4;
    lm_fct_WKS = lm_fct_WKS(:,1:skip_lm_step:end);
    lm_fct_HKS = lm_fct_HKS(:,1:skip_lm_step:end);

    % normalize the descriptors
    lm_fct_WKS = normc(lm_fct_WKS);
    lm_fct_HKS = normc(lm_fct_HKS);

    %% Display the descriptors
    subplot(2,2, 2);
    h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_WKS(:,1), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
    title(['First descriptor WKS']);
    colorbar;

    subplot(2,2, 4);
    h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_HKS(:,1), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
    title(['First descriptor HKS']);
    colorbar;

    % Display all the descriptors with WKS in the same figure
    numPlots = size(lm_fct_WKS,2);
    numRows = ceil(sqrt(numPlots));
    numCols = ceil(numPlots / numRows);
    figure(2);
        for j = 1:numPlots
            numberOfDescriptor = 1+(j-1)*skip_lm_step;
            subplot(numRows, numCols, j);
            h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_WKS(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

        figure(3);
        for j = 1:numPlots
            numberOfDescriptor = 1+(j-1)*skip_lm_step;
            subplot(numRows, numCols, j);
            h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_HKS(:,j), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['Descriptor ' num2str(numberOfDescriptor)]);
        end

    figure(1);
    fprintf('ready for next point');

end



%% load landmarks indices
%landmarks_file = ['../data/' listFolders{i} '/landmarks.txt'];
%lm_idx = load(landmarks_file)+1;
%% Column 1 corresponds to the landmarks indices for the target shape and column 2 for the source shape
%lm_idx_Target = lm_idx(:,1);
%lm_idx_Source = lm_idx(:,2);
%% Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
%timesteps_lm = 100;
%%compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
%numEigs = 100;