function [points, tri] = plane(normal, edgesize)
    [tau1, tau2] = getnormals(normal); %две образующих в плоскости
    DevToSatMatrix = [normal, tau1, tau2];
    points = edgesize/2*[0, 1, -1; 0, 1, 1; 0, -1, 1; 0, -1, -1]';
    tri = [1, 2, 3; 1, 3, 4];    
    for n = 1 : size(points, 2)
        points(:, n) = DevToSatMatrix*points(:, n);
    end
end


