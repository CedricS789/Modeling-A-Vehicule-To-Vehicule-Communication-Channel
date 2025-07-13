function intersection_point = findSegmentIntersection(p1, p2, p3, p4)
% findSegmentIntersection - Finds the intersection point of two line segments.
%
% This function determines if two line segments, L1 (from p1 to p2) and
% L2 (from p3 to p4), intersect. It uses a standard vector cross-product
% method to solve the system of linear equations.
%
% INPUTS:
%   p1, p2 - 1x2 vectors for the endpoints of segment 1.
%   p3, p4 - 1x2 vectors for the endpoints of segment 2.
%
% OUTPUTS:
%   intersection_point - 1x2 vector [x, y] of the intersection, or an
%                        empty array [] if they do not intersect on the segments.

    % Represent segments as P = P_start + t * direction_vector
    v1 = p2 - p1; % Direction vector of segment 1
    v2 = p4 - p3; % Direction vector of segment 2

    % Calculate the 2D "cross product" of the direction vectors.
    % If this is zero, the lines are parallel or collinear.
    v1_cross_v2 = v1(1)*v2(2) - v1(2)*v2(1);

    % Use a small tolerance to handle floating-point inaccuracies.
    if abs(v1_cross_v2) < 1e-10
        intersection_point = []; % Lines are parallel, no unique intersection.
        return;
    end

    % Vector between the starting points of the two segments.
    p3_minus_p1 = p3 - p1;
    
    % Solve for the parameters 't' and 'u'.
    % Intersection Point = p1 + t * v1
    % Intersection Point = p3 + u * v2
    % 't' is the fractional distance along segment 1 where intersection occurs.
    % 'u' is the fractional distance along segment 2 where intersection occurs.
    t = (p3_minus_p1(1)*v2(2) - p3_minus_p1(2)*v2(1)) / v1_cross_v2;
    u = (p3_minus_p1(1)*v1(2) - p3_minus_p1(2)*v1(1)) / v1_cross_v2;

    % For the segments to intersect, both 't' and 'u' must be in [0, 1].
    % A small tolerance is used to robustly include cases where the
    % intersection is exactly at an endpoint.
    if (t >= -1e-9 && t <= 1+1e-9) && (u >= -1e-9 && u <= 1+1e-9)
        % Intersection is valid, calculate the point.
        intersection_point = p1 + t * v1;
    else
        % The infinite lines intersect, but the segments themselves do not.
        intersection_point = [];
    end
end
