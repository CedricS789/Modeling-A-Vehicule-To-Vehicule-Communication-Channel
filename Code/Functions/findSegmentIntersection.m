function intersection_point = findSegmentIntersection(segment1_point1, segment1_point2, segment2_point1, segment2_point2)
% FINDSEGMENTINTERSECTION - Finds the intersection point of two line segments.
%
% This function determines if two line segments, defined by their endpoints,
% intersect. It uses a standard vector cross-product method to solve the
% system of linear equations for the intersection.
%
% The method represents each segment as a parametric equation P = P0 + t*v,
% where 'v' is the direction vector and 't' is a parameter. An intersection
% exists only if the solution yields 't' values between 0 and 1 for BOTH segments.
%
% INPUTS:
%   segment1_point1, segment1_point2 - 1x2 vectors for the endpoints of segment 1.
%   segment2_point1, segment2_point2 - 1x2 vectors for the endpoints of segment 2.
%
% OUTPUTS:
%   intersection_point - 1x2 vector [x, y] of the intersection, or an
%                        empty array [] if they do not intersect on the segments.

    % Define the direction vectors for the two line segments.
    segment1_vector = segment1_point2 - segment1_point1;
    segment2_vector = segment2_point2 - segment2_point1;

    % Calculate the 2D "cross product" of the direction vectors.
    % If this is zero, the lines are parallel.
    v1_cross_v2 = segment1_vector(1)*segment2_vector(2) - segment1_vector(2)*segment2_vector(1);

    % If the cross product is very close to zero, the lines are parallel
    % and are considered not to intersect for this application.
    if abs(v1_cross_v2) < 1e-10
        intersection_point = [];
        return;
    end

    % Vector between the starting points of the two segments.
    start_points_vector = segment2_point1 - segment1_point1;
    
    % Solve for the parameters t and u for the intersection point.
    % Intersection Point = segment1_point1 + t * segment1_vector
    % Intersection Point = segment2_point1 + u * segment2_vector
    t_param = (start_points_vector(1)*segment2_vector(2) - start_points_vector(2)*segment2_vector(1)) / v1_cross_v2;
    u_param = (start_points_vector(1)*segment1_vector(2) - start_points_vector(2)*segment1_vector(1)) / v1_cross_v2;

    % Check if the intersection point lies within BOTH segments.
    % The parameters 't' and 'u' must be between 0 and 1.
    % A small tolerance (1e-9) is used to robustly include endpoints.
    if (t_param >= -1e-9 && t_param <= 1+1e-9) && (u_param >= -1e-9 && u_param <= 1+1e-9)
        % If the conditions are met, calculate the intersection point.
        intersection_point = segment1_point1 + t_param * segment1_vector;
    else
        % Otherwise, the infinite lines intersect, but the segments do not.
        intersection_point = [];
    end
end
