function varargout = SingleIndexToMultiIndex(Index, MaxValues)
%function MultiIndex = SingleIndexToMultiIndex(Index, MaxValues)
%преобразует Index из отрезка [1 prod(MaxValues)] в массив MultiIndex, 
%такой что 1<=MultiIndex(i)<=MaxValues(i)
%
%Вызов
%   MultiIndex = SingleIndexToMultiIndex(Index, MaxValues)
%эквивалентен вызову
%   dim = numel(MaxValues)
%   c = cell([1 dim]);
%   [c{:}] = ind2sub(MaxValues, SingleIndex);
%   MultiIndex = [c{:}];
%но не использует cell arrays.
%
%example:
%   MaxValues = [10 15];
%   res = zeros(1, prod(MaxValues));
%   parfor k = 1 : prod(MaxValues)
%       MultiIndex = SingleIndexToMultiIndex(k, MaxValues);
%       res(k) = sin(MultiIndex(1)/MaxValues(1)*2*pi + MultiIndex(2)/MaxValues(2)*2*pi);
%   end
%   
%   z = zeros(MaxValues);
%   for k = 1 : prod(MaxValues)
%       MultiIndex = SingleIndexToMultiIndex(k, MaxValues)
%       z(MultiIndex(1), MultiIndex(2)) = res(k);
%   end
%   surf(z)
%
%
%see also: MultiIndexToSingleIndex, ind2sub, sub2ind
    MultiIndex = zeros(size(MaxValues));
    Index = Index - 1;
    for i = 1 : numel(MaxValues)
            CurMax = MaxValues(i);
            MultiIndex(i) = mod(Index, CurMax) + 1;
            Index = floor(Index/CurMax);
    end
    assert(numel(find(MultiIndex > MaxValues)) == 0);
    if(nargout == 1)
        varargout{1} = MultiIndex;
    else
        %assert(numel(MaxValues) == numel(nargout));
        varargout = num2cell(MultiIndex);
    end
end
