function intersection_point = findSegmentIntersection(p1, p2, p3, p4)
% Determines if two line segments, L1 (p1 to p2) and L2 (p3 to p4), intersect.
% INPUTS:
%   p1, p2 - Endpoints of segment 1.
%   p3, p4 - Endpoints of segment 2.
% OUTPUTS:
%   intersection_point - Coordinates of intersection, or empty if no intersection.

    v1 = p2 - p1;
    v2 = p4 - p3;

    v1_cross_v2 = v1(1)*v2(2) - v1(2)*v2(1);

    if abs(v1_cross_v2) < 1e-10
        intersection_point = [];
        return;
    end

    p3_minus_p1 = p3 - p1;
    t = (p3_minus_p1(1)*v2(2) - p3_minus_p1(2)*v2(1)) / v1_cross_v2;
    u = (p3_minus_p1(1)*v1(2) - p3_minus_p1(2)*v1(1)) / v1_cross_v2;

    if (t >= -1e-9 && t <= 1+1e-9) && (u >= -1e-9 && u <= 1+1e-9)
        intersection_point = p1 + t * v1;
    else
        intersection_point = [];
    end
end