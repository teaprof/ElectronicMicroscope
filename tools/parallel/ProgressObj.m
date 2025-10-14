classdef ProgressObj < handle
%Usage:
%     P = ProgressObj(N);
%     parfor n = 1 : N
%         x = rand;
%         yy(n) = f(x);
%         P.increase(1);
%     end
%     P.done();
    
    properties
        ntasks
        finished
        queue
        tstart
        ymin
%         listener %can not be perfectly passed to workers
    end
    methods
        function obj = ProgressObj(nTasks)
            if nargin == 0
                nTasks = 0;
            end
            obj.ntasks = nTasks;
            obj.queue = parallel.pool.DataQueue;
            obj.finished = 0;
%             obj.listener = afterEach(obj.queue, @(data)addprogress(obj, data));
            afterEach(obj.queue, @(data)addprogress(obj, data));            
            obj.tstart = tic;
            obj.progressBar(0);
            obj.ymin = []; %no value yet
        end
        function done(obj)
            fprintf('\n');
        end
        function obj = start(obj)
            obj.finished = 0;
        end
        function obj = setNTasks(obj, nTasks)
            %change local copy of obj
            obj.ntasks = nTasks;
            %send reset command to master copy
            send(obj.queue, [1 nTasks]);
        end
        function obj = addprogress(obj, msg)
            switch(msg(1))
                case 1
                    obj.ntasks = msg(2);
                case 2
                    obj.finished = obj.finished + msg(2);
                    if numel(msg) >= 3
                        if numel(obj.ymin) == 0
                            obj.ymin = msg(3);
                        else
                            obj.ymin = min(obj.ymin, msg(3));
                        end
                    end
                    obj.progressBar(1);
            end
        end
        function obj = increase(obj, amount, ymin)
            if nargin < 2
                amount = 1;
            end
            if nargin < 3
                ymin = [];
            end
            send(obj.queue, [2 amount ymin]);
        end
        function percents = get(obj)
            if(obj.finished == 0 || obj.ntasks == 0)
                percents = 0;
            else
                percents = 100.0*obj.finished/obj.ntasks;
            end
        end
        function progressBar(obj, deletePrev)           
            barlen = 40;
            
            percents = obj.get;
            curlen = floor((barlen-2)*percents/100);
            barstr = repmat('=', 1, curlen);
            while(length(barstr) < barlen - 2)
                barstr = [barstr ' '];
            end
            percentstr = sprintf('%8.3f %%', percents);
            
            elapsed = toc(obj.tstart);
            if percents > 0
                total = elapsed*100/percents;
            else
                total = Inf;
            end
            remaining = total - elapsed;
            
            statusstr = sprintf('%7.1f elasped, %7.1f total, %7.1f remaining\n', elapsed, total, remaining);            
            
            
%             if(deletePrev)
%                 backslashes = repmat('\b', 1, 170);
%                 fprintf(backslashes);
%             end
            fprintf('[%s]%s\n%s', barstr, percentstr, statusstr);
            
            if(numel(obj.ymin) > 0)
                fprintf('Minimum found: %g\n', obj.ymin);
            end
            fprintf('\n');
            
%             fprintf(str);
        end
    end
end
