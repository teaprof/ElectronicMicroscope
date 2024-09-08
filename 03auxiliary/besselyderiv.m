function y = besselyderiv(m, r)
%function y = besselyderiv(m, r)
%производная функции Неймана
    y = - bessely(m + 1, r);
    %при r = 0 получаем 0/0, который заменяем просто на 0
    idx = r~= 0;
    y(idx) = y(idx) + m*bessely(m, r(idx))./r(idx);
end
