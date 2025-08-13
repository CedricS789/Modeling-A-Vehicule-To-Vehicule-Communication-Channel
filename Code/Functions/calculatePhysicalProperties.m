function ray_data = calculatePhysicalProperties(ray_coordinates, hit_walls_indices, walls, params)
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
%   ray_data            - A struct containing all calculated physical properties: gamma_tot_n, dist, alpha_n
    
    c = params.c;
    M = size(ray_coordinates, 1) - 2;                                    % The number of reflections is the vertex count minus the TX and RX.
    d_tot = 0;                                                           % This will accumulate the total travel length of the ray path.
    gamma_tot_n = 1;                                                      % This will accumulate the product of reflection losses at each bounce.

    for m = 1:(M + 1)
        start_segment = ray_coordinates(m,:);
        end_segment = ray_coordinates(m+1,:);
        
        d_tot = d_tot + norm(end_segment - start_segment);                          % Add the length of the current segment to the ray's total path length.
        tau_n = d_tot/c;
        
        is_reflection_segment = (m <= M);                                           % True for all segments except the final one that lands at the receiver.
        if is_reflection_segment
            incident_vector = end_segment - start_segment;                          % The vector representing the ray's direction as it strikes the wall.
            reflecting_wall = walls(hit_walls_indices(m));
            
            wall_vector = reflecting_wall.coordinates(2,:) - reflecting_wall.coordinates(1,:);
            normal_vector = [wall_vector(2), -wall_vector(1)];                                      % A vector perpendicular to the wall's surface, essential for angle calculations.
            
            
            cos_theta_n = abs(dot(incident_vector, normal_vector)) / (norm(incident_vector) * norm(normal_vector));
            sin_theta_n = 1 - cos_theta_n^2;
            epsilon_r = reflecting_wall.eps_r;
            
           
            gamma_num = cos_theta_n - sqrt(epsilon_r - sin_theta_n);            % Numerator of the Fresnel formula for perpendicular polarization.
            gamma_den = cos_theta_n + sqrt(epsilon_r - sin_theta_n);            % Denominator of the Fresnel formula for perpendicular polarization.
            gamma_m = gamma_num / gamma_den;
            
            gamma_tot_n = gamma_tot_n * gamma_m;                         % The total reflection loss is the product of losses at each surface.
        end
    end
    
    ray_data.coordinates = ray_coordinates;                                                 % The geometric path vertices from TX to RX.
    ray_data.type = sprintf('%d-Refl', M);                                                  % A string label identifying the type of ray.
    ray_data.distance_total = d_tot;                                                        % The total distance the ray travels.
    ray_data.tau_n = tau_n;
    ray_data.theta_n = acosd(cos_theta_n);
    ray_data.gamma_tot_n = gamma_tot_n;                                                       % The combined reflection coefficient from all bounces.
    ray_data.alpha_n = calculateAlpha_n(ray_data, params);
end