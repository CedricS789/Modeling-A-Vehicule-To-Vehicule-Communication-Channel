function final_ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos)
% Validates a ray path by checking for obstructions.

    final_ray_coordinates = [];
    M = length(all_images_pos) - 1;
    if M < 1
        return;
    end
    potential_ray_coordinates = zeros(M + 2, 2);
    original_tx_pos = all_images_pos{1};
    potential_ray_coordinates(1,:) = original_tx_pos;
    potential_ray_coordinates(end,:) = RX_pos;
    is_path_valid = true;
    current_target_point = RX_pos;

    for k = M:-1:1
        image_source_for_segment = all_images_pos{k+1};
        reflecting_wall_index = processed_wall_indices(k);
        reflecting_wall_coordinates = walls(reflecting_wall_index).coordinates;

        reflection_point = findSegmentIntersection(...
            image_source_for_segment, current_target_point, ...
            reflecting_wall_coordinates(1,:), reflecting_wall_coordinates(2,:));
        if isempty(reflection_point)
            is_path_valid = false;
            break;
        end
        potential_ray_coordinates(k+1,:) = reflection_point;

        for j = 1:length(walls)
            is_this_segments_wall = (j == reflecting_wall_index);
            is_next_segments_wall = (k < M && j == processed_wall_indices(k+1));
            if ~is_this_segments_wall && ~is_next_segments_wall
                if ~isempty(findSegmentIntersection(reflection_point, current_target_point, ...
                        walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
        if ~is_path_valid
            break;
        end
        current_target_point = reflection_point;
    end

    if is_path_valid
        first_reflection_wall_idx = processed_wall_indices(1);
        tx_pos = potential_ray_coordinates(1,:);
        first_reflection_point = potential_ray_coordinates(2,:);
        for j = 1:length(walls)
            if j ~= first_reflection_wall_idx
                if ~isempty(findSegmentIntersection(tx_pos, first_reflection_point, ...
                        walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
    end

    if is_path_valid
        final_ray_coordinates = potential_ray_coordinates;
    end
end