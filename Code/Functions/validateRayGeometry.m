function ray_coordinates = validateRayGeometry(receiver_position, walls, wall_indices_sequence, image_source_sequence)
% validaterayGeometry - Traces and validates the geometry of a potential ray.
%
% Given a complete sequence of image sources, this function determines if a
% physically valid, unobstructed ray exists. It works by tracing the ray
% backwards from the receiver to the transmitter.
%
% INPUTS:
%   receiver_position     - Position of the final receiver [x, y].
%   walls                 - Struct array of wall definitions.
%   wall_indices_sequence - The sequence of wall indices for this ray (e.g., [4 1 2]).
%   image_source_sequence - The cell array of image source positions for this ray.
%
% OUTPUTS:
%   ray_coordinates      - An (N+2)x2 matrix of ray coordinates [TX; R1;..; RX]
%                           if the ray is valid. Otherwise, an empty array [].

    ray_coordinates = [];
    num_reflections = length(image_source_sequence) - 1;
    if num_reflections < 1, return; end

    potential_ray_points = zeros(num_reflections + 2, 2);
    potential_ray_points(1,:) = image_source_sequence{1};
    potential_ray_points(end,:) = receiver_position;

    is_ray_geometry_valid = true;
    traceback_target_point = receiver_position;

    for i = num_reflections:-1:1
        image_source_for_segment = image_source_sequence{i+1};
        reflecting_wall_index = wall_indices_sequence(i);
        reflecting_wall_coords = walls(reflecting_wall_index).coords;
        
        reflection_point = findSegmentIntersection(image_source_for_segment, traceback_target_point, reflecting_wall_coords(1,:), reflecting_wall_coords(2,:));

        if isempty(reflection_point)
            is_ray_geometry_valid = false; break;
        end
        potential_ray_points(i+1,:) = reflection_point;

        for j = 1:length(walls)
            is_reflecting_wall_of_this_segment = (j == reflecting_wall_index);
            is_reflecting_wall_of_next_segment = (i < num_reflections && j == wall_indices_sequence(i+1));
            
            if ~is_reflecting_wall_of_this_segment && ~is_reflecting_wall_of_next_segment
                if ~isempty(findSegmentIntersection(reflection_point, traceback_target_point, walls(j).coords(1,:), walls(j).coords(2,:)))
                    is_ray_geometry_valid = false; break;
                end
            end
        end
        if ~is_ray_geometry_valid, break; end
        traceback_target_point = reflection_point;
    end

    if is_ray_geometry_valid
        first_reflection_wall_idx = wall_indices_sequence(1);
        for j = 1:length(walls)
            if j ~= first_reflection_wall_idx
                if ~isempty(findSegmentIntersection(potential_ray_points(1,:), potential_ray_points(2,:), walls(j).coords(1,:), walls(j).coords(2,:)))
                    is_ray_geometry_valid = false; break;
                end
            end
        end
    end

    if is_ray_geometry_valid
        ray_coordinates = potential_ray_points;
    end
end
