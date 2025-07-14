function valid_rays_found = findReflectedRaysRecursive(current_tx_pos, RX_pos, walls, num_reflections_remaining, processed_wall_indices, all_images_pos, params)
% findReflectedRaysRecursive - Recursively finds reflected rays using the image method.
%
% This function is the core algorithm of the ray tracer. It builds a reflection tree
% where each node is an image source created by reflecting a previous source
% across a wall.
% INPUTS:
%   current_tx_pos              - Position of the current source: real TX or an image.
%   RX_pos                      - Position of the final receiver (constant).
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

    % --- BASE CASE: No more reflections needed ---
    % When the desired number of reflections is reached, the full path from the
    % final image source to the receiver can be traced and validated.
    if num_reflections_remaining == 0
        % The full reflection sequence is now hypothesized; validate its geometric path.
        ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos);
        
        % A non-empty result from validateRayPath means the path is unobstructed
        % and physically possible.
        if ~isempty(ray_coordinates)
            ray_data = calculatePhysicalProperties(ray_coordinates, processed_wall_indices, walls, params);  % If the path is valid, compute its physical characteristics.
            valid_rays_found = {ray_data};                                                                   % Package the complete, valid ray data to be returned up the call stack.
        end
        return;                                                                                              % Terminate this branch of recursion and return the findings.
    end

    % --- RECURSIVE STEP: Find the next reflection ---
    last_wall_index = 0;                                                                  % Will hold the index of the wall used in the previous reflection step.
    if ~isempty(processed_wall_indices)
        last_wall_index = processed_wall_indices(end);                                         % Get the index of the most recent wall to avoid immediate re-reflection.
    end
    
    % Iterate through all possible walls to serve as the next reflecting surface.
    for i = 1:length(walls)
        % This check prevents the ray from immediately reflecting back and forth off
        % the same wall, which is a physically redundant path.
        if i == last_wall_index
            continue;
        end
        
        % This is the core of the Method of Images: create a new virtual source
        % by reflecting the current source across the plane of the wall.
        new_image_source = findSymmetricAcrossLine(current_tx_pos, walls(i).coordinates);
        
        % --- Make the recursive call ---
        % Descend one level deeper into the reflection tree with the new state.
        rays_from_branch = findReflectedRaysRecursive(new_image_source, RX_pos, walls, ...           % Pass the new virtual source as the starting point for the next level.
            num_reflections_remaining - 1, [processed_wall_indices, i], ...                          % Decrement remaining reflections and append the current wall to the path history.
            [all_images_pos, {new_image_source}], params);                                           % Append the new virtual source to the list of image positions.
        
        % Collect any valid rays that were found in the deeper recursive calls
        % and add them to this level's list of findings.
        valid_rays_found = [valid_rays_found, rays_from_branch];
    end
end