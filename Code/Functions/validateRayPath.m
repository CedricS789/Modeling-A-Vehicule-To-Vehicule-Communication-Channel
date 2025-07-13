function final_ray_coordinates = validateRayPath(rx_pos, walls, processed_wall_indices, all_tx_images_pos)
% validateRayPath Traces and validates the geometry of a potential ray path.
%
% This function determines if a physically valid, unobstructed ray exists
% by tracing the path backwards from the receiver to the transmitter,
% identifying reflection points, and checking each segment for obstructions.
%
% INPUTS:
%   rx_pos                  - Position of the final receiver [x, y].
%   walls                   - Struct array of wall definitions. Each struct
%                           contains 'coordinates' (start and end points).
%   processed_wall_indices  - Sequence of wall indices for this ray, indicating
%                           reflection order from transmitter to receiver.
%   all_tx_images_pos       - Cell array of image source positions.
%                           all_tx_images_pos{1} is the true transmitter (TX) position.
%                           all_tx_images_pos{i+1} is the image source for the i-th reflection.
%
% OUTPUTS:
%   final_ray_coordinates - An (N+2)x2 matrix of ray coordinates from original
%                           image going through reflections and to the receptor [TX; R1;..; RN; RX]
%                           if the ray is valid (N is the number of reflections).
%                           Returns an empty array [] if the ray path is invalid
%                           (e.g., no intersection or obstruction).
    final_ray_coordinates = [];                                                           % Default to an invalid path; will be populated only if all checks pass.
    K = length(all_tx_images_pos) - 1;                                      % The number of reflections is one less than the total number of sources (original + images).
    % LOS ray (K < 1) is not handled by this function.
    if K < 1                                                                % This function only handles reflected paths.
        return;                                                                           % Exit if no reflections are specified.
    end
    % Pre-allocate matrix for all ray vertices: TX, N reflection points, RX.
    potential_ray_coordinates = zeros(K + 2, 2);                            % Allocate memory for the path vertices: TX, N reflection points, and RX.
    original_tx_pos = all_tx_images_pos{1};                                               % The true transmitter.
    potential_ray_coordinates(1,:) = original_tx_pos;                                     % The path must start at the true transmitter (TX).
    potential_ray_coordinates(end,:) = rx_pos;                                            % The path must terminate at the receiver (RX).
    is_path_valid = true;                                                                 % Assume the path is valid until a check fails.
    current_target_point = rx_pos;                                                        % Begin the backward trace from the receiver's location.

    % --- Trace backwards from RX to TX and check for obstructions ---
    for i = K:-1:1                                                          % Trace path segments backwards, from the final reflection toward the first.
        image_source_for_segment = all_tx_images_pos{i+1};                                % The ray segment appears to originate from the corresponding image source.
        reflecting_wall_index = processed_wall_indices(i);                                % Identify the wall associated with this step of the reflection sequence.
        reflecting_wall_coordinates = walls(reflecting_wall_index).coordinates;           % Get the start and end coordinates of this reflecting wall.

        % Find the reflection point on the current reflecting wall.
        % This is the intersection of the line from the current target point
        % to the image source with the reflecting wall segment.
        reflection_point = findSegmentIntersection(...                                     % Calculate the exact point of reflection on the wall.
            image_source_for_segment, current_target_point, ...                            % A line is cast from the image source to the target of the previous segment.
            reflecting_wall_coordinates(1,:), reflecting_wall_coordinates(2,:));           % The intersection with the finite wall segment is the reflection point.
        if isempty(reflection_point)                                                      % If the line of sight to the image misses the physical wall segment.
            is_path_valid = false;                                                        % the geometric path is impossible.
            break;                                                                        % Stop processing, as the entire ray path is invalid.
        end
        potential_ray_coordinates(i+1,:) = reflection_point;                              % Store this valid reflection point as a vertex in the ray path.

        % Check for obstructions on the current ray segment (reflection_point to current_target_point).
        for j = 1:length(walls)                                                           % Iterate through all walls in the environment to check for blockages.
            % Skip checking the current reflecting wall and the next reflecting wall
            % (if applicable) as they are part of the intended path.
            is_this_segments_wall = (j == reflecting_wall_index);                         % The ray is allowed to hit the wall it's reflecting from.
            is_next_segments_wall = (i < K && j == processed_wall_indices(i+1)); % And the wall of the next reflection (the previous in this backward trace).
            if ~is_this_segments_wall && ~is_next_segments_wall                           % For all other walls...
                % If an intersection with any other wall is found, the segment is obstructed.
                if ~isempty(findSegmentIntersection(reflection_point, current_target_point, ...
                        walls(j).coordinates(1,:), walls(j).coordinates(2,:)))            % check if the current ray segment is blocked.
                    is_path_valid = false;                                                % If it is blocked, the path is invalid.
                    break;                                                                % Exit the inner loop; no need to check other walls.
                end
            end
        end
        if ~is_path_valid                                                                 % If the inner loop found an obstruction.
            break;                                                                        % exit the main loop as well.
        end
        current_target_point = reflection_point;                                          % For the next backward step, this reflection point is the new target.
    end

    % --- Final check for the first segment (TX to first reflection point) ---
    if is_path_valid                                                                      % If all previous segments are valid, check the first one.
        first_reflection_wall_idx = processed_wall_indices(1);                            % Get the wall corresponding to the first reflection after the TX.
        tx_pos = potential_ray_coordinates(1,:);                                          % Get the starting position of the ray (the transmitter).
        first_reflection_point = potential_ray_coordinates(2,:);                          % Get the first reflection point calculated in the loop above.
        for j = 1:length(walls)                                                           % Iterate through all walls to check for obstructions.
            % The first segment should only interact with the first reflecting wall.
            % Any intersection with another wall indicates an obstruction.
            if j ~= first_reflection_wall_idx                                             % The ray is allowed to hit its target reflecting wall.
                if ~isempty(findSegmentIntersection(tx_pos, first_reflection_point, walls(j).coordinates(1,:), walls(j).coordinates(2,:)))  % Check for intersections with any OTHER wall.
                    is_path_valid = false;                                                % If it hits another wall first, it's obstructed.
                    break;                                                                % Exit the loop, the path is invalid.
                end
            end
        end
    end

    % If all checks pass, assign the calculated coordinates to the output.
    if is_path_valid                                                                      % If the path has survived all geometric and obstruction checks.
        final_ray_coordinates = potential_ray_coordinates;                                % commit the calculated vertices to the output.
    end
end
