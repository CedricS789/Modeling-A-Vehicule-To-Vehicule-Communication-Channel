function symmetric_point = findSymmetricAcrossLine(point, line_coordinates)
% Reflects a 2D point across a line segment.
% INPUTS:
%   point - [x, y] of the point to reflect.
%   line_coordinates - [x1 y1; x2 y2] defining the line.
% OUTPUTS:
%   symmetric_point - [x, y] of the resulting image point.

    line_point_1 = line_coordinates(1,:);
    line_direction_vector = line_coordinates(2,:) - line_point_1;
    normal_vector = [line_direction_vector(2), -line_direction_vector(1)];
    vector_to_point = point - line_point_1;
    scale = 2 * dot(vector_to_point, normal_vector) / dot(normal_vector, normal_vector);
    symmetric_point = point - scale * normal_vector;
end