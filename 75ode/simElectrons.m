function traj = simElectrons(r0, v0, EMfieldSolution, tspan, electronsPerParticle)
    x0 = PhaseToX(r0, v0);
    N = size(r0, 2);    
%1. Interaction with retarded potentials    
%     sol = ode_retarded(@(t, x, thistory, xhistory) rightSideRetarded(t, x, N, xhistory, thistory, EMfieldSolution, electronsPerParticle), tspan, x0, 'acceleration');

%2. Interaction with classic potentials
    sol = ode45(@(t, x) rightSideClassic(t, x, N, EMfieldSolution, electronsPerParticle), tspan, x0);

%3. No interaction 
%    sol = ode45(@(t, x) rightSideNoInteraction(t, x, N, EMfieldSolution), tspan, x0);
    
    traj.t = sol.x;
    traj.rv = sol.y;
    traj.Nparticles = N;
    [X, Y, Z] = odesolToPhase(sol);
    traj.X = X; traj.Y = Y; traj.Z = Z;
%     checkEnergyLaw(traj);
end

function sol = ode_retarded(fcn, time, x0, stepalgorithm)
    %инициализируем историю - надо две точки, чтобы можно было сторить
    %интерполяцию. Моменты времени выбираем t=-1 c (берём с запасом) и t = time(1),
    %предполагая, что time(1) > -1
    x(:, 1) = x0;
    x(:, 2) = x0;
    t(1) = -1;
    t(2) = time(1);
    curx = x0;    
    curt = time(1);
    curstep = (time(2) - time(1))/20;
    finish = false;
    tstart = tic;
    while(~finish)
        switch stepalgorithm
            case 'adaptive' %эта ветка не тестировалась после исправлений
                [curstep, finish] = choosestep(2*curstep, curt, time(2));
                if(curstep < 0)
                    break;
                end
                GoodStep = false;
                while(GoodStep == false)
                    %Один большой шаг
                    nextx = rkstep_retarded(fcn, curt, curx, curstep, t, x);
                    
                    %Два маленьких шага
                    x_intermediate = rkstep_retarded(fcn, curt, curx, curstep/2, t, x);
                    t_intermediate = curt + curstep/2;
                    nextx_2 = rkstep_retarded(fcn, t_intermediate, x_intermediate, curstep/2, [t t_intermediate], [x x_intermediate]);
                    
                    xt_abs_err = abs(nextx - nextx_2);
                    idx = find(abs(nextx_2) > 1e-20);
                    eps = abs(xt_abs_err(idx)./nextx_2(idx));
                    
                    GoodStep = max(xt_abs_err)<1e-4 && max(eps) < 1e-6;
                    if(~GoodStep)
                        curstep = curstep/2;
                    end
                end
                [~, finish] = choosestep(curstep, curt, time(2));
                curx = nextx_2;
                curt = curt + curstep;
            case 'acceleration'
                xt = fcn(curt, curx, t, x);
                N = size(x0, 1)/6;
                [vt, at] = XtoPhase(xt, N);
                curstep = min(1e4/max(abs(at), [], 'all'), 1e-13);
%                 curstep = max(min(curstep, time(2) - curt), 1e-15);
                [curstep, finish] = choosestep(curstep, curt, time(2));
                curx = rkstep_retarded(fcn, curt, curx, curstep, t, x);                
                curt = curt + curstep;
            otherwise
                error('stepalgorithm is unknown');
        end
        x(:, end+1) = curx;
        t(:, end+1) = curt;
        fprintf('curt = %10.6g, curstep = %12.8g\n', curt, curstep);
        telapsed = toc(tstart);
        testimated = telapsed/(curt-time(1))*(time(2)-time(1));
        tremaining = testimated - telapsed;
        fprintf('Elapsed %f secs, estimated total = %f secs, remaining = %f secs\n', telapsed, testimated, tremaining);
    end
    
    sol.x = t(2:end);
    sol.y = x(:, 2:end);
end

function [step, finish] = choosestep(curstep, curt, tend)
%function [step, finish] = choosestep(curstep, curt, tend)
%finish - истина, если этот шаг последний
    finish = false;
    maxstep = tend - curt;
    step = min(curstep, maxstep);
    nextt = curt + step;
    assert(tend ~= 0);
    if((tend - nextt)/tend < 10*eps)
        step = maxstep;
        finish = true;
    end        
end

function x = rkstep_retarded(fcn, t0, x0, step, thistory, xhistory)
    k1 = fcn(t0, x0, thistory, xhistory);
    t2 = t0 + step/2;
    x2 = x0 + step/2*k1;
    thist_ext = [thistory t2];
    xhist_ext = [xhistory x2];
    k2 = fcn(t2, x2, thist_ext, xhist_ext);
    t3 = t2;
    x3 = x0 + step/2*k2;
    thist_ext(end) = t3;
    xhist_ext(:, end) = x3;
    k3 = fcn(t3, x3, thist_ext, xhist_ext);
    t4 = t0 + step;
    x4 = x0 + step*k3;    
    thist_ext(end) = t4;
    xhist_ext(:, end) = x4;
    k4 = fcn(t4, x4, thist_ext, xhist_ext);
    x = x0 + step/6*(k1 + 2*k2 + 2*k3 + k4);
end

function [xhist, vhist] = PhaseToHistory(xhistory, particlen)
    N = size(xhistory, 1)/6;
    idx = (particlen-1)*3 + 1;
    xhist = xhistory(idx:idx+2, :);
    idx = 3*N + idx;
    vhist = xhistory(idx:idx+2, :);
end


