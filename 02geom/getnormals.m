function [n, b] = getnormals(v)
%function [n, b] = getnormals(v)
%������� ��� ������� � ���������� �������
%[v, n, b] - ������ ������ ��������

    %������� �����-������ �����������, �� ������������ v
    assert(numel(find(isnan(v))) == 0, 'getnormals: vector contains NaN');
    assert(numel(find(isinf(v))) == 0, 'getnormals: vector contains Inf');
    assert(sum(abs(v)) > 0, 'getnormals: vector should not be zero');
    n = ones(size(v));
    for k = 1 : numel(v)
        if(v(k) ~= 0)
            n(k) = 0;
            numerator = dot(n, v);
            n(k) = -numerator/v(k);
            break;
        end
    end
    n = n/norm(n);
    %������� ���������
    b = cross(v, n);
    b = b/norm(b);    
end
