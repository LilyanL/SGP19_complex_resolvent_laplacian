classdef PairShapes
    properties
        shape_source
        shape_target

        pair_folder_name

        landmarks_source
        landmarks_target

        mappings
        mappings_Labels

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
       
    end
 end