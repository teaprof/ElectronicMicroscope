function y = besseljderiv(m, r)
%function y = besseljderiv(m, r)
%производная функции Бесселя
    y = - besselj(m + 1, r);
    idx = r ~= 0;
    %при r = 0 получаем 0/0, который заменяем просто на 0
    y(idx) = y(idx) + m*besselj(m, r(idx))./r(idx);
end
