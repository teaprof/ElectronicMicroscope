function Tasks = PrepareParallelTasks(varargin)
%function Tasks = PrepareParallelTasks(varargin)
%Создает массив Tasks, каждая строка которого - комбинация вида (n1, n2, ..., nk), 
%такая что ni берется из массива varargin(i)
%
%example:
%   X = [0:0.1:10];
%   Y = [0:0.1:10];
%   MaxValues = [numel(X) numel(Y)];
%   Tasks = PrepareParallelTasks(X, Y);
%   parfor k = 1 : size(Tasks, 1)
%       x = Tasks(k, 1);
%       y = Tasks(k, 2);
%       res(k) = sin(x + y);
%   end
%   
%   z = zeros(MaxValues);
%   for k = 1 : prod(MaxValues)
%       MultiIndex = SingleIndexToMultiIndex(k, MaxValues);
%       z(MultiIndex(1), MultiIndex(2)) = res(k);
%   end
%   surf(z)
%
%see also: MultiIndexToSingleIndex, SingleIndexToMultiIndex
    MaxIndexes = zeros(1, numel(varargin));
    for k = 1 : numel(varargin)
        MaxIndexes(k) = numel(varargin{k});
    end 
    Tasks = zeros(prod(MaxIndexes), numel(varargin));
    for k = 1 : prod(MaxIndexes)
        M = SingleIndexToMultiIndex(k, MaxIndexes);
        for m = 1 : numel(MaxIndexes)
            Tasks(k, m) = varargin{m}(M(m));
        end
    end
end