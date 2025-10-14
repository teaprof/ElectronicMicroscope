function [first, last] = getSubJob(Njobs, curWorker, totalWorkers)
%function [first, last] = getSubJob(N, curWorker, totalWorkers)
    first = floor(Njobs*(curWorker-1)/totalWorkers) + 1;
    last = floor(Njobs*curWorker/totalWorkers);
    if curWorker == totalWorkers
        last = Njobs;
    end
end
