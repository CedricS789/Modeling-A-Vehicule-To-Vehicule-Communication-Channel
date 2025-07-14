function valid_rays_found = findReflectedRaysRecursive(current_tx_pos, RX_pos, walls, num_reflections_remaining, processed_wall_indices, all_images_pos, params)
% Recursively finds reflected rays using the image method.
%
% INPUTS:
%   current_tx_pos              - Position of the current source: real TX or an image.
%   RX_pos                      - Position of the final receiver which is constant.
%   walls                       - Struct array of wall definitions.
%   num_reflections_remaining   - The number of additional reflections to find.
%   processed_wall_indices      - The sequence of wall indices used so far.
%   all_images_pos              - The sequence of image source positions generated so far.
%   params                      - Struct with simulation parameters.
%
% OUTPUTS:
%   valid_rays_found    - A cell array of completed, valid ray data structs
%                         found from this recursive branch.

    valid_rays_found = {};   % Initialize a cell array to store results from this recursive path.

    if num_reflections_remaining == 0
        ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos); 
        if ~isempty(ray_coordinates)
            ray_data = calculatePhysicalProperties(ray_coordinates, processed_wall_indices, walls, params);  % If the path is valid, compute its physical characteristics.
            valid_rays_found = {ray_data};                                                                   % Package the complete, valid ray data to be returned up the call stack.
        end
        return;                                                                                              % Terminate this branch of recursion and return the findings.
    end

    last_wall_index = 0;                                                                        % Will hold the index of the wall used in the previous reflection step.
    if ~isempty(processed_wall_indices)
        last_wall_index = processed_wall_indices(end);                                          % Get the index of the most recent wall to avoid immediate re-reflection.
    end
    
    for i = 1:length(walls)
        if i == last_wall_index
            continue;
        end
        
        new_image_source = findSymmetricAcrossLine(current_tx_pos, walls(i).coordinates);
        
        % Make the recursive call
        rays_from_branch = findReflectedRaysRecursive(new_image_source, RX_pos, walls, ...           % Pass the new virtual source as the starting point for the next level.
            num_reflections_remaining - 1, [processed_wall_indices, i], ...                          % Decrement remaining reflections and append the current wall to the path history.
            [all_images_pos, {new_image_source}], params);                                           % Append the new virtual source to the list of image positions.
        

        valid_rays_found = [valid_rays_found, rays_from_branch];
    end
end