function valid_rays_found = findReflectedRaysRecursive(current_source_pos, rx_pos, walls, reflections_to_go, wall_indices_path, image_source_path, params)
% findReflectedRaysRecursive - Recursively finds reflected rays using the image method.
%
% This function is the core of the ray tracer. It builds a reflection tree
% where each node is an image source created by reflecting a previous source
% across a wall.
%
% INPUTS:
%   current_source_pos  - Position of the current source (real TX or an image).
%   rx_pos              - Position of the final receiver (constant).
%   walls               - Struct array of wall definitions.
%   reflections_to_go   - The number of additional reflections to find.
%   wall_indices_path   - The sequence of wall indices used so far.
%   image_source_path   - The sequence of image source positions generated so far.
%   params          - Struct with simulation parameters.
%
% OUTPUTS:
%   valid_rays_found    - A cell array of completed, valid ray data structs
%                         found from this recursive branch.

    valid_rays_found = {};

    % --- BASE CASE: No more reflections needed ---
    % If we have reached the desired reflection order, we can trace the final
    % path from the last image source to the receiver and validate it.
    if reflections_to_go == 0
        % Validate the geometry of the potential ray path.
        ray_coords = validateRayPath(rx_pos, walls, wall_indices_path, image_source_path);
        
        % If the path is geometrically valid (unobstructed), calculate its
        % physical properties (distance, gain, etc.).
        if ~isempty(ray_coords)
            ray_data = calculatePhysicalProperties(ray_coords, wall_indices_path, walls, params);
            valid_rays_found = {ray_data}; % Return the completed ray data.
        end
        return; % End this branch of the recursion.
    end

    % --- RECURSIVE STEP: Find the next reflection ---
    last_wall_index = 0;
    if ~isempty(wall_indices_path)
        last_wall_index = wall_indices_path(end);
    end
    
    % Iterate through all possible walls for the next reflection.
    for i = 1:length(walls)
        % To avoid trivial reflections (e.g., back and forth on the same wall),
        % do not reflect across the same wall consecutively.
        if i == last_wall_index
            continue;
        end
        
        % Create a new image source by reflecting the current source across wall 'i'.
        new_image_source = reflectPointAcrossLine(current_source_pos, walls(i).coords);
        
        % --- Make the recursive call ---
        % Descend one level deeper into the reflection tree.
        rays_from_branch = findReflectedRaysRecursive(new_image_source, rx_pos, walls, ...
            reflections_to_go - 1, [wall_indices_path, i], ...
            [image_source_path, {new_image_source}], params);
        
        % Collect all valid rays found from this deeper branch.
        valid_rays_found = [valid_rays_found, rays_from_branch];
    end
end
