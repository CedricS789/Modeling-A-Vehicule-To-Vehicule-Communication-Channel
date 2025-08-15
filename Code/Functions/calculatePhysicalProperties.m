function ray_data = calculatePhysicalProperties(ray_coordinates, hit_walls_indices, walls, params)
% Calculates physical properties for a valid ray path.
    
    c = params.c;
    M = size(ray_coordinates, 1) - 2;
    d_tot = 0;
    gamma_tot_n = 1;

    for m = 1:(M + 1)
        start_segment = ray_coordinates(m,:);
        end_segment = ray_coordinates(m+1,:);
        
        d_tot = d_tot + norm(end_segment - start_segment);
        tau_n = d_tot/c;
        
        is_reflection_segment = (m <= M);
        if is_reflection_segment
            incident_vector = end_segment - start_segment;
            reflecting_wall = walls(hit_walls_indices(m));
            
            wall_vector = reflecting_wall.coordinates(2,:) - reflecting_wall.coordinates(1,:);
            normal_vector = [wall_vector(2), -wall_vector(1)];
            
            cos_theta_n = abs(dot(incident_vector, normal_vector)) / (norm(incident_vector) * norm(normal_vector));
            sin_theta_n = 1 - cos_theta_n^2;
            epsilon_r = reflecting_wall.eps_r;
            
            gamma_num = cos_theta_n - sqrt(epsilon_r - sin_theta_n);
            gamma_den = cos_theta_n + sqrt(epsilon_r - sin_theta_n);
            gamma_m = gamma_num / gamma_den;
            
            gamma_tot_n = gamma_tot_n * gamma_m;
        end
    end
    
    ray_data.coordinates = ray_coordinates;
    ray_data.type = sprintf('%d-Refl', M);
    ray_data.distance_total = d_tot;
    ray_data.tau_n = tau_n;
    ray_data.theta_n = acosd(cos_theta_n);
    ray_data.gamma_tot_n = gamma_tot_n;
    ray_data.alpha_n = calculateAlpha_n(ray_data, params);
end