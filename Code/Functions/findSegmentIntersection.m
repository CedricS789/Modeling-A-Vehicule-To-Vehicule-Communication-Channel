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

    % Define segments using parametric form: P = P_start + t * direction_vector
    v1 = p2 - p1;                                                                 % Vector representing the direction and magnitude of segment 1.
    v2 = p4 - p3;                                                                 % Vector representing the direction and magnitude of segment 2.

    % The system is solved using Cramer's rule. The denominator is the
    % determinant of the matrix formed by the direction vectors.
    v1_cross_v2 = v1(1)*v2(2) - v1(2)*v2(1);                                       % This is the 2D equivalent of a cross-product's magnitude.

    % If the determinant is near zero, the vectors are parallel (or collinear),
    % meaning the lines don't have a single, unique intersection point.
    if abs(v1_cross_v2) < 1e-10                                                   % Use a small tolerance to account for floating-point errors.
        intersection_point = [];                                                  % Return empty as parallel lines do not intersect.
        return;
    end

    % Create a vector from the start point of segment 1 to the start of segment 2.
    p3_minus_p1 = p3 - p1;
    
    % Solve for parameters 't' and 'u' which represent the fractional distance
    % along each segment where the intersection occurs.
    % Intersection point I satisfies: I = p1 + t*v1 and I = p3 + u*v2.
    t = (p3_minus_p1(1)*v2(2) - p3_minus_p1(2)*v2(1)) / v1_cross_v2;               % 't' is the parameter for segment 1 (p1 to p2).
    u = (p3_minus_p1(1)*v1(2) - p3_minus_p1(2)*v1(1)) / v1_cross_v2;               % 'u' is the parameter for segment 2 (p3 to p4).

    % An intersection exists on the segments only if both parameters are between
    % 0 and 1, inclusive. This confirms the intersection lies within the bounds
    % of both finite line segments.
    if (t >= -1e-9 && t <= 1+1e-9) && (u >= -1e-9 && u <= 1+1e-9)                  % Use tolerance to robustly include intersections at the endpoints.
        % The intersection is valid and lies on both segments.
        intersection_point = p1 + t * v1;                                         % Calculate the coordinates using the parameter t and segment 1's equation.
    else
        % The segments do not intersect, even if their infinite line extensions do.
        intersection_point = [];
    end
end