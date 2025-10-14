function [X, Y, Z] = mycylinder(startpoint, endpoint, r)
%function [X, Y, Z] = mycylinder(startpoint, endpoint, r)
%function obj = mycylinder(startpoint, endpoint, r)
% Создёт цилиндр
%   startpoint - начальная точка оси, 3x1
%   endpoint   - конечная точка оси, 3x1
%   r          - радиус цилиндра
% чтобы нарисовать цилинд, вызывать mesh(X, Y, Z);
    [X, Y, Z] = cylinder(r);
    Z = Z*norm(endpoint-startpoint);
    e1 = [1; 0; 0];
    e2 = [0; 1; 0];
    e3 = [0; 0; 1];
    dir = endpoint - startpoint;
    dir = dir/norm(dir);
    [n, b] = getnormals(dir);
    dcm = [e1, e2, e3]'*[n, b, dir];
    for n = 1 : numel(X)
        p = dcm*[X(n); Y(n); Z(n)];
        X(n) = p(1); Y(n) = p(2); Z(n) = p(3);
    end
    if nargout < 3
        X = mesh(X, Y, Z);
    end
end
