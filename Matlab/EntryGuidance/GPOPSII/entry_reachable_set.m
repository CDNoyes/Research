function entry_reachable_set()

dtr = pi/180;

min_alt = 8e3;

Nper = 15;

% System Models:
mars = Mars();
vm = VehicleModel();

% auxdata.ref = ref;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.delta.CD = 0;
% Initial states:
%[radius long lat velocity fpa heading]
for N = [ 550, 600, 675, 800, 1000] % 400
    switch N
        case 0
            x0 = [3540e3; 0*dtr; 0*dtr; 5505; -14.15*dtr; 0*dtr]';
        case 400
            x0 = [ 3.44798709e+06  1.23845619e-01  0  5.48299238e+03 -1.69371202e-01 0];
        case 550
            x0 = [ 3.42606332e+06  1.68402816e-01  0  4.55018727e+03 -1.00352125e-01 0];
        case 600
            x0 = [ 3.42262762e+06  1.80735324e-01  0  3.91773878e+03 -6.07013321e-02 0];
        case 675
            x0 = [ 3.41971007e+06  1.95905383e-01  0  3.04810497e+03 -5.40734510e-02 0];
        case 800
            x0 = [ 3.41704669e+06  2.14174154e-01  0  2.04726686e+03 -5.15120036e-02 0];
            
        case 1000
            x0 = [ 3.41443052e+06  2.32450788e-01  0  1.19474151e+03 -1.86866763e-02 0];
    end

    % x0 = cell2mat(x0);
    
    sols = {};
    
    for iter = 1:(Nper+4)
        % Target info:
        if iter > 4
            lons = linspace(sols{2}.state(end,2), sols{3}.state(end,2)*0.99, Nper);
            target_lon = lons(iter-4);
        end
        
        
        
        % Bounds
        t0 = 0;
        tfl = 1;
        tfu = 500;
        
        bounds.phase.initialtime.lower = t0;
        bounds.phase.initialtime.upper = t0;
        bounds.phase.finaltime.lower = tfl;
        bounds.phase.finaltime.upper = tfu;
        
        lb = [mars.radiusEquatorial+min_alt
            -pi/2 * 0
            -pi/2
            300
            -45*pi/180 %-25 degrees, FPA
            -pi/2      % Heading
            -pi/2]';   % Bank angle
        
        ub = [x0(1)+100
            pi
            pi/2
            x0(4)+500
            10*dtr % FPA constraint degrees
            pi
            pi/2]';
        

        vel_max = 470;
        
        
            
        if iter == 1 || iter == 4 % Zero crossrange 
            xfl = [lb(1), lb(2), lb(3)*0, lb(4), lb(5), lb(6), lb(7)];
            xfu = [ub(1), ub(2), ub(3)*0, vel_max, ub(5), ub(6), ub(7)];
            
        elseif iter > 4% Fixed lon i.e. downrange
            xfl = [lb(1), target_lon, lb(3), lb(4), lb(5), lb(6), lb(7)];
            xfu = [ub(1), target_lon, ub(3), vel_max, ub(5), ub(6), ub(7)];    
            
        else % Free lat/lon
            xfl = lb;
            xfu = [ub(1), ub(2:3), vel_max, ub(5), ub(6), ub(7)];
        end
        
        bounds.phase.initialstate.lower =[x0, lb(7)*0];
        bounds.phase.initialstate.upper = [x0, ub(7)*0];
        bounds.phase.state.lower = [lb];
        bounds.phase.state.upper = [ub];
        bounds.phase.finalstate.lower = [xfl];
        bounds.phase.finalstate.upper = [xfu];
        
        bank_rate = 20*dtr;
        bounds.phase.control.lower = -bank_rate;
        bounds.phase.control.upper = bank_rate;
        bounds.phase.integral.lower = 0;
        bounds.phase.integral.upper = 5000;
        
        % Initial Guess
        tGuess = [t0;0.5*(tfu+tfl)];
        guess.phase.time = tGuess;
        guess.phase.state = [[x0,0];0.5*(xfu+xfl)];
        guess.phase.control = zeros(2,1);
        guess.phase.integral = 0;
        
        % Initial Mesh
        meshphase.colpoints = 4*ones(1,10);
        meshphase.fraction = .1*ones(1,10);
        setup.name = 'EntryReachableSet';
        setup.functions.continuous = @EntryConstraints;
        if iter == 1 || iter == 2
            setup.functions.endpoint = @MinDR;
        elseif iter == 3 || iter == 4
            setup.functions.endpoint = @MaxDR;
        else
            setup.functions.endpoint = @MaxCR;
        end
        setup.auxdata = auxdata;
        setup.mesh.phase = meshphase;
        setup.bounds = bounds;
        
        setup.guess = guess;
        setup.nlp.solver = 'snopt';
        setup.derivatives.supplier = 'sparseCD';
        setup.derivatives.derivativelevel = 'first';
        setup.scales.method = 'automatic-guessUpdate';
        setup.method = 'RPM-Differentiation'; % Differentiation
        setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao'
        setup.mesh.tolerance = 1e-5; % default 1e-3
        setup.mesh.colpointsmin = 3;
        setup.mesh.colpointsmax = 10;
        setup.mesh.maxiterations = 20;
        setup.displaylevel = 1; % default 2
        sol = gpops2(setup);
        Sol = sol.result.solution.phase;
        Sol = rmfield(Sol, 'controlRadau');
        Sol = rmfield(Sol, 'integral');
        Sol = rmfield(Sol, 'costate');
        Sol = rmfield(Sol, 'timeRadau');
        Sol = rmfield(Sol, 'integralmultipliers');
        sols{end+1} = Sol;
        save(['./results/',num2str(N),'/rs_',num2str(iter),'.mat'], '-struct', 'Sol')
        
    end
end
end
function output = MaxDR( input )
xf = input.phase.finalstate;
output.objective = -xf(2);
end

function output = MinDR( input )
xf = input.phase.finalstate;
output.objective = xf(2);
end

function output = MaxCR( input )
xf = input.phase.finalstate;
output.objective = -xf(3);
end