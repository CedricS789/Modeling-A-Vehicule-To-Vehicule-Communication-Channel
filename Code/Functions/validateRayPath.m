function final_ray_coords = validateRayPath(rx_pos, walls, wall_indices_path, image_source_path)
% validateRayPath - Traces and validates the geometry of a potential ray path.
%
% Given a sequence of image sources, this function determines if a physically
% valid, unobstructed ray exists. It works by tracing the path backwards
% from the receiver to the transmitter, finding reflection points and
% checking each segment for obstructions.
%
% INPUTS:
%   rx_pos              - Position of the final receiver [x, y].
%   walls               - Struct array of wall definitions.
%   wall_indices_path   - The sequence of wall indices for this ray (e.g., [4 1 2]).
%   image_source_path   - The cell array of image source positions for this ray.
%
% OUTPUTS:
%   final_ray_coords    - An (N+2)x2 matrix of ray coordinates [TX; R1;..; RX]
%                         if the ray is valid. Otherwise, an empty array [].

    final_ray_coords = [];
    num_reflections = length(image_source_path) - 1;
    if num_reflections < 1, return; end

    % Pre-allocate matrix for the ray's vertices.
    % Path: [TX, Reflection_1, ..., Reflection_N, RX]
    potential_ray_coords = zeros(num_reflections + 2, 2);
    potential_ray_coords(1,:) = image_source_path{1}; % Start with the true TX
    potential_ray_coords(end,:) = rx_pos;             % End with the RX

    is_path_valid = true;
    current_target_point = rx_pos;

    % --- Trace backwards from RX to TX ---
    for i = num_reflections:-1:1
        % The line segment for this trace is from the current target point
        % (RX or a reflection point) to the image source that created it.
        image_source_for_segment = image_source_path{i+1};
        reflecting_wall_index = wall_indices_path(i);
        reflecting_wall_coords = walls(reflecting_wall_index).coords;
        
        % The reflection point must lie on the actual wall segment.
        reflection_point = findSegmentIntersection(...
            image_source_for_segment, current_target_point, ...
            reflecting_wall_coords(1,:), reflecting_wall_coords(2,:));

        % If there's no intersection, the geometry is invalid.
        if isempty(reflection_point)
            is_path_valid = false;
            break;
        end
        potential_ray_coords(i+1,:) = reflection_point;

        % --- Check for obstructions for the current ray segment ---
        % The segment is from reflection_point to current_target_point.
        for j = 1:length(walls)
            % A segment cannot be obstructed by its own reflecting wall.
            % Also, it cannot be obstructed by the wall of the *next* reflection
            % point, as it is supposed to touch that wall.
            is_this_segments_wall = (j == reflecting_wall_index);
            is_next_segments_wall = (i < num_reflections && j == wall_indices_path(i+1));
            
            if ~is_this_segments_wall && ~is_next_segments_wall
                % Check for intersection with any other wall in the environment.
                if ~isempty(findSegmentIntersection(reflection_point, current_target_point, walls(j).coords(1,:), walls(j).coords(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
        if ~is_path_valid, break; end
        
        % For the next iteration, the new target is the reflection point we just found.
        current_target_point = reflection_point;
    end

    % --- Final check for the first segment (TX to first reflection) ---
    if is_path_valid
        first_reflection_wall_idx = wall_indices_path(1);
        tx_pos = potential_ray_coords(1,:);
        first_reflection_point = potential_ray_coords(2,:);
        
        for j = 1:length(walls)
            % This segment should only touch the first reflecting wall.
            if j ~= first_reflection_wall_idx
                if ~isempty(findSegmentIntersection(tx_pos, first_reflection_point, walls(j).coords(1,:), walls(j).coords(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
    end

    % If the path has remained valid through all checks, store it.
    if is_path_valid
        final_ray_coords = potential_ray_coords;
    end
end
