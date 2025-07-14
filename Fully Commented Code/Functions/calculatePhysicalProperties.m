function ray_data = calculatePhysicalProperties(ray_coordinates, hit_walls_indices, walls, params)
% calculatePhysicalProperties - Computes physical properties for a valid ray.
%
% For a geometrically valid ray path, this function calculates the total
% distance, the product of reflection coefficients, and the complex channel gain.
%
% INPUTS:
%   ray_coordinates     - (N+2)x2 matrix of ray vertices [TX; R1;..; RX].
%   hit_walls_indices   - The sequence of wall indices for this ray.
%   walls               - Struct array of wall definitions.
%   params              - Struct with simulation parameters.
%
% OUTPUTS:
%   ray_data            - A struct containing all calculated physical properties: gamma_prod, dist, alpha_n

    K = size(ray_coordinates, 1) - 2;                                    % The number of reflections is the vertex count minus the TX and RX.
    d_tot = 0;                                                           % This will accumulate the total travel length of the ray path.
    gamma_prod = 1;                                                      % This will accumulate the product of reflection losses at each bounce.

    % Iterate over each straight-line segment that makes up the full ray path.
    for n = 1:(K + 1)
        start_segment = ray_coordinates(n,:);
        end_segment = ray_coordinates(n+1,:);
        
        d_tot = d_tot + norm(end_segment - start_segment);                          % Add the length of the current segment to the ray's total path length.
        
        is_reflection_segment = (n <= K);                                           % True for all segments except the final one that lands at the receiver.
        if is_reflection_segment
            incident_vector = end_segment - start_segment;                          % The vector representing the ray's direction as it strikes the wall.
            reflecting_wall = walls(hit_walls_indices(n));
            
            wall_vector = reflecting_wall.coordinates(2,:) - reflecting_wall.coordinates(1,:);
            normal_vector = [wall_vector(2), -wall_vector(1)];                                      % A vector perpendicular to the wall's surface, essential for angle calculations.
            
            % Use the dot product to find the cosine of the angle between the
            % incident ray and the line normal to the wall's surface.
            cos_theta_n = abs(dot(incident_vector, normal_vector)) / (norm(incident_vector) * norm(normal_vector));
            sin_theta_n = 1 - cos_theta_n^2;
            epsilon_r = reflecting_wall.eps_r;                                      % The wall's relative permittivity, a material property governing reflection.
            
            % The Fresnel reflection coefficient models how much of the wave's
            % energy is reflected, based on angle, polarization, and material.
            gamma_num = cos_theta_n - sqrt(epsilon_r - sin_theta_n);            % Numerator of the Fresnel formula for perpendicular polarization.
            gamma_den = cos_theta_n + sqrt(epsilon_r - sin_theta_n);            % Denominator of the Fresnel formula for perpendicular polarization.
            reflection_coeff = gamma_num / gamma_den;
            
            gamma_prod = gamma_prod * reflection_coeff; % The total reflection loss is the product of losses at each surface.
        end
    end
    
    % Package all calculated ray characteristics into a single struct for output.
    ray_data.coordinates = ray_coordinates;                                                 % The geometric path vertices from TX to RX.
    ray_data.type = sprintf('%d-Refl', K);                                                  % A string label identifying the type of ray, example: "2-Refl".
    ray_data.distance_total = d_tot;                                                        % The total distance the ray travels.
    ray_data.gamma_prod = gamma_prod;                                                       % The combined reflection coefficient from all bounces.
    
    % Combine path loss and reflection effects to find the final complex gain.
    ray_data.alpha_n = calculateAlpha_n(ray_data, params);
end