function reflected_point = reflectPointAcrossWall(point_to_reflect, wall_coordinates)
% REFLECTPOINTACROSSWALL - Reflects a 2D point across a line segment.
%
% This function calculates the image of a point when reflected across the
% infinite line defined by a wall's two endpoints. This is a core
% geometric operation for the Method of Images.
%
% The formula for reflection of a point P across a line defined by point L0
% and normal vector n is: P' = P - 2 * dot(P - L0, n) / |n|^2 * n.
%
% INPUTS:
%   point_to_reflect - 1x2 vector [x, y] of the point.
%   wall_coordinates - 2x2 matrix [x1 y1; x2 y2] defining the line of reflection.
%
% OUTPUTS:
%   reflected_point  - 1x2 vector [x, y] of the resulting image point.

    % A point on the line of reflection.
    line_point = wall_coordinates(1,:);
    
    % Vector defining the direction of the wall.
    wall_direction_vector = wall_coordinates(2,:) - line_point;
    
    % A vector perpendicular (normal) to the wall's direction.
    % For a 2D vector [x, y], the perpendicular vector is [y, -x] or [-y, x].
    normal_vector = [wall_direction_vector(2), -wall_direction_vector(1)];
    
    % Vector from the point on the line to the point we want to reflect.
    point_vector = point_to_reflect - line_point;
    
    % Apply the vector reflection formula.
    % The term dot(point_vector, normal_vector) is the projection of point_vector onto the normal.
    % We travel twice this projected distance along the normal vector from the original point.
    reflected_point = point_to_reflect - 2 * dot(point_vector, normal_vector) / dot(normal_vector, normal_vector) * normal_vector;
end
