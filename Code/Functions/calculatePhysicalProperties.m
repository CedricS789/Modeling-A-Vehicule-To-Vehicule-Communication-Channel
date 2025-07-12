function ray_data = calculatePhysicalProperties(ray_coordinates, wall_indices_sequence, walls, simulation_parameters)
% calculatePhysicalProperties - Calculates physical properties for a valid ray.
%
% This function takes a set of geometrically valid ray coordinates and
% computes the essential physical parameters of that propagation ray.
%
% INPUTS:
%   ray_coordinates      - (N+2)x2 matrix of ray vertices [TX; R1;..; RX].
%   wall_indices_sequence - The sequence of wall indices for this ray.
%   walls                 - Struct array of wall definitions.
%   simulation_parameters - Struct with simulation parameters.
%
% OUTPUTS:
%   ray_data             - A struct containing all the calculated physical
%                           properties for the ray.

    num_reflections = size(ray_coordinates, 1) - 2;
    total_distance = 0;
    cumulative_gamma = 1;

    for i = 1:(num_reflections + 1)
        start_segment = ray_coordinates(i,:);
        end_segment = ray_coordinates(i+1,:);
        
        total_distance = total_distance + norm(end_segment - start_segment);
        
        is_reflection_segment = (i <= num_reflections);
        if is_reflection_segment
            incident_vector = end_segment - start_segment;
            reflecting_wall = walls(wall_indices_sequence(i));
            
            wall_vector = reflecting_wall.coords(2,:) - reflecting_wall.coords(1,:);
            normal_vector = [wall_vector(2), -wall_vector(1)];
            
            cos_theta_i = abs(dot(incident_vector/norm(incident_vector), normal_vector/norm(normal_vector)));
            sin_theta_i_sq = 1 - cos_theta_i^2;
            
            epsilon_r = reflecting_wall.eps_r;
            
            gamma_numerator = cos_theta_i - sqrt(epsilon_r - sin_theta_i_sq);
            gamma_denominator = cos_theta_i + sqrt(epsilon_r - sin_theta_i_sq);
            reflection_coefficient_perp = gamma_numerator / gamma_denominator;
            
            cumulative_gamma = cumulative_gamma * reflection_coefficient_perp;
        end
    end
    
    ray_data.coordinates = ray_coordinates;
    ray_data.type = sprintf('%d-Refl', num_reflections);
    ray_data.dist = total_distance;
    ray_data.gamma_prod = cumulative_gamma;
    
    ray_data.gain = calculateAlpha_n(ray_data, simulation_parameters);
end
