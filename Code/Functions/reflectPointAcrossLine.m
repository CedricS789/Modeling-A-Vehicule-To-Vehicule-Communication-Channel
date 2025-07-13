function reflected_point = reflectPointAcrossLine(point, line_coords)
% reflectPointAcrossLine - Reflects a 2D point across a line segment.
%
% This function calculates the image of a point when reflected across the
% infinite line defined by the line segment's two endpoints. This is a core
% geometric operation for the Method of Images.
%
% Formula: P' = P - 2 * projection_of_v_onto_n
% where P' is the reflected point, P is the original point, v is a vector
% from the line to P, and n is the normal vector to the line.
%
% INPUTS:
%   point       - 1x2 vector [x, y] of the point to reflect.
%   line_coords - 2x2 matrix [x1 y1; x2 y2] defining the line of reflection.
%
% OUTPUTS:
%   reflected_point - 1x2 vector [x, y] of the resulting image point.

    % A point on the line of reflection.
    line_point_1 = line_coords(1,:);
    
    % Vector defining the direction of the line.
    line_direction_vector = line_coords(2,:) - line_point_1;
    
    % A vector perpendicular (normal) to the line's direction.
    % For a 2D vector [dx, dy], a normal is [dy, -dx].
    normal_vector = [line_direction_vector(2), -line_direction_vector(1)];
    
    % Vector from a point on the line to the point we want to reflect.
    vector_to_point = point - line_point_1;
    
    % Apply the vector reflection formula.
    % The dot product term calculates the projection length. We travel
    % twice this distance along the normal vector from the original point.
    scale = 2 * dot(vector_to_point, normal_vector) / dot(normal_vector, normal_vector);
    
    reflected_point = point - scale * normal_vector;
end
