function valid_rays_from_recursion = findReflectedRaysRecursive(current_source_position, receiver_position, walls, reflections_remaining, wall_indices_sequence, image_source_sequence, simulation_parameters)
% findReflectedraysRecursive - Recursively finds reflection rays using the image method.
%
% This is the heart of the ray tracer. It builds a "tree" of reflections.
% Each call to this function represents a node in the tree.
%
% INPUTS:
%   current_source_position - Position of the current source (can be the real TX or an image).
%   receiver_position       - Position of the final receiver (this never changes).
%   walls                   - Struct array of wall definitions.
%   reflections_remaining   - The number of additional reflections needed for this branch.
%   wall_indices_sequence   - The sequence of wall indices used to get to this point.
%   image_source_sequence   - The sequence of image source positions generated so far.
%   simulation_parameters   - Struct with simulation parameters.
%
% OUTPUTS:
%   valid_rays_from_recursion - A cell array of completed, valid ray data structs
%                                found from this node and its children.

    valid_rays_from_recursion = {};

    % --- BASE CASE ---
    if reflections_remaining == 0
        ray_coordinates = validateRayGeometry(receiver_position, walls, wall_indices_sequence, image_source_sequence);
        if ~isempty(ray_coordinates)
            ray_data_struct = calculatePhysicalProperties(ray_coordinates, wall_indices_sequence, walls, simulation_parameters);
            valid_rays_from_recursion = {ray_data_struct};
        end
        return;
    end

    % --- RECURSIVE STEP ---
    last_wall_index = 0;
    if ~isempty(wall_indices_sequence)
        last_wall_index = wall_indices_sequence(end);
    end
    
    for i = 1:length(walls)
        if i == last_wall_index, continue; end
        
        new_image_source = reflectPointAcrossWall(current_source_position, walls(i).coords);
        
        rays_from_this_branch = findReflectedRaysRecursive(new_image_source, receiver_position, walls, ...
            reflections_remaining - 1, [wall_indices_sequence, i], ...
            [image_source_sequence, {new_image_source}], simulation_parameters);
        
        valid_rays_from_recursion = [valid_rays_from_recursion, rays_from_this_branch];
    end
end
