classdef PairShapes
    properties
        shape_source
        shape_target

        pair_folder_name

        landmarks_source
        landmarks_target

        noise_vector_source %if non was applied, should be 0
        noise_vector_target %if non was applied, should be 0
        transform_source %cell array of registrations. Each cell contains cells of registrations (one rigid3dtransform per segment). If non was applied, should be identity
        transform_target %cell array of registrations. Each cell contains cells of registrations (one rigid3dtransform per segment). If non was applied, should be identity

        mappings %cell array of mappings. Each cell contains cells of mappings (one cell per segment)
        mappings_Labels %cell array of mappings names
        registrations_source_to_target %cell array of registrations. Each cell contains cells of registrations (one rigid3dtransform per segment)

        descriptors_global_source
        descriptors_global_target
        descriptors_local_source
        descriptors_local_target

        descriptors_global_source_labels
        descriptors_global_target_labels
        descriptors_local_source_labels
        descriptors_local_target_labels

        trajectories_source %mapping can directly be applied
        trajectories_target %mapping is trickier to apply (reverse search)

        

    end
    methods

        % computeRegistrationErrors(obj)
        % This function compares the transforms between source and target shapes for different methods and segments.
        % It calculates the registration errors for translation and rotation in each dimension.
        %
        % Input:
        % - obj: The object containing the mappings and transforms
        %
        % Output:
        % - registration_errors: A cell array containing the registration errors for each method and segment
        function registration_errors = computeRegistrationErrors(obj)
            num_methods = numel(obj.mappings_Labels);
            num_segments = numel(obj.transform_source{1});
            registration_errors = cell(num_methods, num_segments);
            
            for i = 1:num_methods
                for j = 1:num_segments
                    registration_error = RegistrationError();
                    registration_error.error_X = obj.transform_target{j}.Translation(1) - obj.registrations_source_to_target{i}{j}.Translation(1);
                    registration_error.error_Y = obj.transform_target{j}.Translation(2) - obj.registrations_source_to_target{i}{j}.Translation(2);
                    registration_error.error_Z = obj.transform_target{j}.Translation(3) - obj.registrations_source_to_target{i}{j}.Translation(3);
                    registration_error.error_R = obj.transform_target{j}.Rotation(1) - obj.registrations_source_to_target{i}{j}.Rotation(1);
                    registration_error.error_P = obj.transform_target{j}.Rotation(2) - obj.registrations_source_to_target{i}{j}.Rotation(2);
                    registration_error.error_Yaw = obj.transform_target{j}.Rotation(3) - obj.registrations_source_to_target{i}{j}.Rotation(3);
                    
                    registration_errors{i, j} = registration_error;
                end
            end
        end

       
        % Input:
        %   - registration_errors: a cell array containing the registration errors for each method and segment.
        % 
        % Output:
        %   - average_errors: a matrix containing the average registration errors for each method and component (X, Y, Z, R, P, Yaw).
       
        function average_errors = computeAverageRegistrationErrors(obj, registration_errors)
            num_methods = size(registration_errors, 1);
            num_segments = size(registration_errors, 2);
            num_components = 6; % X, Y, Z, R, P, Yaw
            
            average_errors = zeros(num_methods, num_components);
            
            for i = 1:num_methods
                for j = 1:num_segments
                    error = registration_errors{i, j};
                    average_errors(i, 1) = average_errors(i, 1) + error.error_X;
                    average_errors(i, 2) = average_errors(i, 2) + error.error_Y;
                    average_errors(i, 3) = average_errors(i, 3) + error.error_Z;
                    average_errors(i, 4) = average_errors(i, 4) + error.error_R;
                    average_errors(i, 5) = average_errors(i, 5) + error.error_P;
                    average_errors(i, 6) = average_errors(i, 6) + error.error_Yaw;
                end
            end
            
            average_errors = average_errors / num_segments;
        end

    end

end