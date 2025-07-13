function ray_data = calculatePhysicalProperties(ray_coords, wall_indices_path, walls, params)
% calculatePhysicalProperties - Computes physical properties for a valid ray.
%
% For a geometrically valid ray path, this function calculates the total
% distance, the product of reflection coefficients, and the complex channel gain.
%
% INPUTS:
%   ray_coords          - (N+2)x2 matrix of ray vertices [TX; R1;..; RX].
%   wall_indices_path   - The sequence of wall indices for this ray.
%   walls               - Struct array of wall definitions.
%   params          - Struct with simulation parameters.
%
% OUTPUTS:
%   ray_data            - A struct containing all calculated physical properties.

    num_reflections = size(ray_coords, 1) - 2;
    total_distance = 0;
    cumulative_reflection_coeff = 1;

    % Iterate over each segment of the ray path.
    for i = 1:(num_reflections + 1)
        start_segment = ray_coords(i,:);
        end_segment = ray_coords(i+1,:);
        
        % Accumulate total path distance.
        total_distance = total_distance + norm(end_segment - start_segment);
        
        % If this segment ends in a reflection, calculate the coefficient.
        is_reflection_segment = (i <= num_reflections);
        if is_reflection_segment
            incident_vector = end_segment - start_segment;
            reflecting_wall = walls(wall_indices_path(i));
            
            wall_vector = reflecting_wall.coords(2,:) - reflecting_wall.coords(1,:);
            % The normal vector is perpendicular to the wall's direction vector.
            normal_vector = [wall_vector(2), -wall_vector(1)];
            
            % Calculate the cosine of the angle of incidence (theta_i).
            cos_theta_i = abs(dot(incident_vector, normal_vector)) / (norm(incident_vector) * norm(normal_vector));
            
            % Using sin^2 + cos^2 = 1
            sin_theta_i_sq = 1 - cos_theta_i^2;
            
            epsilon_r = reflecting_wall.eps_r;
            
            % Fresnel reflection coefficient for perpendicular (E-field) polarization.
            % This is a standard formula in electromagnetics.
            gamma_num = cos_theta_i - sqrt(epsilon_r - sin_theta_i_sq);
            gamma_den = cos_theta_i + sqrt(epsilon_r - sin_theta_i_sq);
            reflection_coeff = gamma_num / gamma_den;
            
            cumulative_reflection_coeff = cumulative_reflection_coeff * reflection_coeff;
        end
    end
    
    % Store all calculated data in a struct.
    ray_data.path = ray_coords;
    ray_data.type = sprintf('%d-Refl', num_reflections);
    ray_data.dist = total_distance;
    ray_data.gamma_prod = cumulative_reflection_coeff;
    
    % Calculate the final complex gain for this ray.
    ray_data.alpha_n = calculateAlpha_n(ray_data, params);
end
