function symmetric_point = findSymmetricAcrossLine(point, line_coordinates)
% Reflects a 2D point across a line segment.
%
% INPUTS:
%   point       - [x, y] of the point to reflect.
%   line_coordinates - [x1 y1; x2 y2] defining the line of reflection.
%
% OUTPUTS:
%   symmetric_point - 1x2 vector [x, y] of the resulting image point.

    line_point_1 = line_coordinates(1,:);                                               % Select one point on the line to act as a reference origin.
    
    line_direction_vector = line_coordinates(2,:) - line_point_1;                       % Calculate the vector that defines the direction of the line.

    normal_vector = [line_direction_vector(2), -line_direction_vector(1)];              % Find a vector that is normal to the reflection line.
    
    vector_to_point = point - line_point_1;                                             % Define a vector from the reference point on the line to the point being symmetric.
    
    scale = 2 * dot(vector_to_point, normal_vector) / dot(normal_vector, normal_vector);
    
    symmetric_point = point - scale * normal_vector;                                    % Subtracting the scaled normal vector moves the original point to its mirror image.
end