function valid_rays_found = findReflectedRaysRecursive(current_tx_pos, RX_pos, walls, num_reflections_remaining, processed_wall_indices, all_images_pos, params)
% Recursively finds reflected rays using the image method.

    valid_rays_found = {};

    if num_reflections_remaining == 0
        ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos);
        if ~isempty(ray_coordinates)
            ray_data = calculatePhysicalProperties(ray_coordinates, processed_wall_indices, walls, params);
            valid_rays_found = {ray_data};
        end
        return;
    end

    last_wall_index = 0;
    if ~isempty(processed_wall_indices)
        last_wall_index = processed_wall_indices(end);
    end
    
    for i = 1:length(walls)
        if i == last_wall_index
            continue;
        end
        
        new_image_source = findSymmetricAcrossLine(current_tx_pos, walls(i).coordinates);
        
        rays_from_branch = findReflectedRaysRecursive(new_image_source, RX_pos, walls, ...
            num_reflections_remaining - 1, [processed_wall_indices, i], ...
            [all_images_pos, {new_image_source}], params);
        
        valid_rays_found = [valid_rays_found, rays_from_branch];
    end
end