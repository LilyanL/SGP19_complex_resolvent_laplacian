%% Let's setup everything nicely
clc; clear all, close all; % close all;
set(0,'DefaultFigureWindowStyle','docked') % docked or normal
dbstop if error % debug mode

%% some initialization
textboxIndexSource = [];
textboxIndexTarget = [];

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

%% Select and load the target meshe using the GUI
isLoadedTarget = false;
[file,path] = uigetfile('*.off','Select the target shape to load');
if isequal(file,0)
    disp('User selected Cancel');
    return;
else
    isLoadedTarget = true;
    disp(['User selected ', fullfile(path,file)]);
end

% load the mesh
shapeTarget = MESH.MESH_IO.read_shape([path, file]);

% center data
shapeTarget.surface.VERT = shapeTarget.surface.VERT - mean(shapeTarget.surface.VERT);

% preprocess the mesh
shapeTarget = MESH.preprocess(shapeTarget,meshOptions{:});

% compute eigenvecs and eigenvalues
A = shapeTarget.evecs(:,1:k1);
Ew = shapeTarget.evals(1:k1);

%% Prepare the figures with some useful attributes
% Figure 1: select the entry point and display specific descriptors
selectionFigure = figure('Name', 'Select the entry point');

%Source: shape
subplot(4,2, [1 3], 'ButtonDownFcn', {@Callback, 1});

%Source: descriptor 1
subplot(4,2, [2], 'ButtonDownFcn', {@Callback, 2});

%Source: descriptor 2
subplot(4,2, [4], 'ButtonDownFcn', {@Callback, 3});

%Target
subplot(4,2, [5 7], 'ButtonDownFcn', {@Callback, 4});

%Target: descriptor 1
subplot(4,2, [6], 'ButtonDownFcn', {@Callback, 5});

%Target: descriptor 2
subplot(4,2, [8], 'ButtonDownFcn', {@Callback, 6});

% Figure 2: display all the descriptors with WKS in the same figure
WKSFigure = figure(2);

% Figure 3: display all the descriptors with HKS in the same figure
HKSFigure = figure(3);

%% Display the source and target shapes
figure(1);
title('Select the entry point');
subplot(4,2, [1 3]);
trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), 'FaceColor', 'interp', 'EdgeColor', 'none', 'ButtonDownFcn', {@Callback, 1});
alpha(0.5);
view(-75, 5);
hold on;

subplot(4,2, [5 7], 'ButtonDownFcn', {@Callback, 4});
trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), 'FaceColor', 'interp', 'EdgeColor', 'none', 'ButtonDownFcn', {@Callback, 4});
alpha(0.5);
view(-75, 5);
hold on;

%% Select a 3D point on the mesh using the GUI and compute the descriptors
pointCloudSource = shapeSource.surface.VERT';
pointCloudTarget = shapeTarget.surface.VERT';
% set the callback, pass pointCloudSource to the callback function
%set(selectionFigure, 'WindowButtonDownFcn', {@callbackClickA3DPoint, pointCloudSource}); 

% wait for user to close the figure
while (ishandle(selectionFigure) || ishandle(WKSFigure)||ishandle(HKSFigure))

    % wait until the user has clicked on the mesh
    while waitforbuttonpress ~= 0
        pause 1;
    end

    FigH = ancestor(gca, 'figure');
    % disp(Index)
    % UserData = [get(FigH, 'UserData'), Index];
    % set(FigH, 'Userdata', UserData);
    % 
    C = get(selectionFigure, 'UserData')

    
    if isempty(C) || (C(1) ~= 1 && C(1)~=4)
        continue;
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
    rotatedPointCloudSource = rot * pointCloudSource; 
    rotatedPointCloudTarget = rot * pointCloudTarget;

    % the clicked point represented in the view frame
    rotatedPointFront = rot * point' ;

    % find the nearest neighbour to the clicked point 
    pointCloudIndex = 0;
    if(C == 1)
        pointCloudIndex = dsearchn(rotatedPointCloudSource(1:2,:)', ... 
            rotatedPointFront(1:2));
    else 
        if (C == 4)
            pointCloudIndex = dsearchn(rotatedPointCloudTarget(1:2,:)', ... 
                rotatedPointFront(1:2));
        end
    end

    h = findobj(gca,'Tag','pt'); % try to find the old point

    selectedPoint = [];
    if(C == 1)
        selectedPoint = pointCloudSource(:, pointCloudIndex); 
    else
        if(C == 4)
            selectedPoint = pointCloudTarget(:, pointCloudIndex); 
        end
    end

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

    if (C ==1)
        fprintf('you clicked on point number %d on the source shape\n', pointCloudIndex);
        %Add a text box under the plot with the point index number
        if(~isempty(textboxIndexSource))
            delete(textboxIndexSource);
        end

        str = num2str(pointCloudIndex);
        textboxIndexSource = annotation('textbox', [0.2, 0.6, 0.1, 0.1], 'String', str);       

    else
        if (C == 4)
            fprintf('you clicked on point number %d on the target shape\n', pointCloudIndex);
            %Add a text box under the plot with the point index number
            if(~isempty(textboxIndexTarget))
                delete(textboxIndexTarget);
            end

            str = num2str(pointCloudIndex);
            textboxIndexTarget = annotation('textbox', [0.2, 0.1, 0.1, 0.1], 'String', str);  
        end
    end


    %% Compute an visualize the descriptors
    lm_idx = pointCloudIndex;

    % Compute the landmarks based descriptors using compute_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip)
    timesteps_lm = 100;

    %compute_chosen_local_descriptors_with_landmarks(S,numEigs,landmarks,t,num_skip, method)
    lm_fct_WKS = []; lm_fct_HKS=[];
    numEigs = 100;
    method = 'WKS';
    if(C == 1)
        lm_fct_WKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx,timesteps_lm,1, method);
    else
        if(C == 4)
            lm_fct_WKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeTarget,numEigs,lm_idx,timesteps_lm,1, method);
        end
    end

    method = 'HKS';
    if(C == 1)
        lm_fct_HKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeSource,numEigs,lm_idx,timesteps_lm,1, method);
    else
        if(C == 4)
            lm_fct_HKS = fMAP.compute_chosen_local_descriptors_with_landmarks(shapeTarget,numEigs,lm_idx,timesteps_lm,1, method);
        end
    end
    
    % ignore some of the descriptors
    skip_lm_step = 4;
    lm_fct_WKS = lm_fct_WKS(:,1:skip_lm_step:end);
    lm_fct_HKS = lm_fct_HKS(:,1:skip_lm_step:end);

    % normalize the descriptors
    lm_fct_WKS = normc(lm_fct_WKS);
    lm_fct_HKS = normc(lm_fct_HKS);

    %% Display the descriptors
    if (C == 1)
    subplot(4,2, 2);
    h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_WKS(:,1), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
    title(['First descriptor WKS']);
    colorbar;

    subplot(4,2, 4);
    h = trisurf(shapeSource.surface.TRIV, shapeSource.surface.VERT(:,1), shapeSource.surface.VERT(:,2), shapeSource.surface.VERT(:,3), lm_fct_HKS(:,1), 'FaceColor', 'interp');
    set(h, 'edgecolor', 'none');
    axis equal; axis off; hold on;
    title(['First descriptor HKS']);
    colorbar;

    else
        if (C == 4)
            subplot(4,2, 6);
            h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), lm_fct_WKS(:,1), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['First descriptor WKS']);
            colorbar;

            subplot(4,2, 8);
            h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), lm_fct_HKS(:,1), 'FaceColor', 'interp');
            set(h, 'edgecolor', 'none');
            axis equal; axis off; hold on;
            title(['First descriptor HKS']);
            colorbar;
        end
    end

    % Display all the descriptors with WKS in the same figure
    numPlots = size(lm_fct_WKS,2);
    numRows = ceil(sqrt(numPlots));
    numCols = ceil(numPlots / numRows);

    if(C == 1)
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

    else
        if(C==4)

            figure(4);
            for j = 1:numPlots
                numberOfDescriptor = 1+(j-1)*skip_lm_step;
                subplot(numRows, numCols, j);
                h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), lm_fct_WKS(:,j), 'FaceColor', 'interp');
                set(h, 'edgecolor', 'none');
                axis equal; axis off; hold on;
                title(['Descriptor ' num2str(numberOfDescriptor)]);
            end

            figure(5);
            for j = 1:numPlots
                numberOfDescriptor = 1+(j-1)*skip_lm_step;
                subplot(numRows, numCols, j);
                h = trisurf(shapeTarget.surface.TRIV, shapeTarget.surface.VERT(:,1), shapeTarget.surface.VERT(:,2), shapeTarget.surface.VERT(:,3), lm_fct_HKS(:,j), 'FaceColor', 'interp');
                set(h, 'edgecolor', 'none');
                axis equal; axis off; hold on;
                title(['Descriptor ' num2str(numberOfDescriptor)]);
            end
        end
    end

    figure(1);
    fprintf('ready for next point');

end

function Callback(ObjectH, EventData, Index)
FigH = ancestor(ObjectH, 'figure');
disp(Index)
UserData = [Index];%UserData = [get(FigH, 'UserData'), Index];
set(FigH, 'Userdata', UserData);

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