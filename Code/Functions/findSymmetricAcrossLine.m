function symmetric_point = findSymmetricAcrossLine(point, line_coords)
% findSymmetricAcrossLine - Reflects a 2D point across a line segment.
%
% This function calculates the image of a point when symmetric across the
% infinite line defined by the line segment's two endpoints. This is a core
% geometric operation for the Method of Images.
%
% Formula: P' = P - 2 * projection_of_v_onto_n
% where P' is the symmetric point, P is the original point, v is a vector
% from the line to P, and n is the normal vector to the line.
%
% INPUTS:
%   point       - 1x2 vector [x, y] of the point to reflect.
%   line_coords - 2x2 matrix [x1 y1; x2 y2] defining the line of reflection.
%
% OUTPUTS:
%   symmetric_point - 1x2 vector [x, y] of the resulting image point.

    line_point_1 = line_coords(1,:);                                              % Select one point on the line to act as a reference origin.
    
    line_direction_vector = line_coords(2,:) - line_point_1;                      % Calculate the vector that defines the direction of the line.
    
    % To reflect a point, we need to move it along a line perpendicular
    % to the line of reflection. For a 2D vector [dx, dy], its normal is [dy, -dx].
    normal_vector = [line_direction_vector(2), -line_direction_vector(1)];        % Find a vector that is perpendicular (normal) to the reflection line.
    
    vector_to_point = point - line_point_1;                                       % Define a vector from the reference point on the line to the point being symmetric.
    
    % This scale factor is derived from the vector projection formula. It
    % determines how far to move along the normal to get to the symmetric position.
    % The dot product projects our vector onto the normal, and the factor of 2
    % moves it twice that distance (once to the line, once past it).
    scale = 2 * dot(vector_to_point, normal_vector) / dot(normal_vector, normal_vector);
    
    symmetric_point = point - scale * normal_vector;                              % Subtracting the scaled normal vector moves the original point to its mirror image.
end