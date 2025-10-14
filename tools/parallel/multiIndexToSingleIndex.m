function IDX = multiIndexToSingleIndex(MultiIndex, MaxValues)
%
%function IDX = multiIndexToSingleIndex(MultiIndex, MaxValues)
%отображает множество индексов (i1, i2, ..., ik) на один индекс
%используется для организации одного parfor-цикла вместо множества
%вложенных циклов
%MultiIndex - нумерация начинается с 1, а не с нуля
%MaxValues - количество элементов по каждой размерности
%
%Функция полностью совместима с sub2ind, но более удобна в использовании,
%если размерность известна только в момент выполнения программы.
%
% Вызов
%   IDX = multiIndexToSingleIndex(MultiIndex, MaxValues)
% полностью эквивалентен вызову
%   subs = num2cell(MultiIndex);
%   IDX = sub2ind(MaxValues, subs{:});
% но не использует при этом cell arrays.
%
%see also: 
%   singleIndexToMultiIndex
%   sub2ind ind2sub
    IDX = 0;
    M = 1;
    for k = 1 : numel(MultiIndex)
        IDX = IDX + (MultiIndex(k) - 1)*M;
        M = M*MaxValues(k);
    end
    IDX = IDX + 1;
end