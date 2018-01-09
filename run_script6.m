%Simplifying McNabb's run_script.m

clear
close all

%% Run Script Settings
Save_mat            = 1;
Plot_VFs            = false;
Plot_Values         = true;
Display_Loop        = false;
Show_Results_Table  = true;
period				= 6;
ant_scale			= 3.6; % num_ants = ceil(num_FOBs/ant_scale)
cl_scale			= 5;   % cl_size  = ceil(num_FOBs/cl_scale)

%% Begin Procedure
load Baseline4
rng('default')

global	DISTANCES       ...
        PROBABILITIES   ...
		NUM_BLUE        ...
		FOB             ...
        VEHICLE         ...
		PARAMETERS      ...
		ACO_PARAMETERS	...
		HEURISTICS

DISTANCES                       = MapBlue.ShortestPaths.Distances;	% Get shortest paths
PROBABILITIES                   = MapBlue.Probabilities(1).Prob;	% Get route success probabilities
NUM_BLUE                        = size(DISTANCES,1);                % Get number customers + depot

FOB.NUM                         = NUM_BLUE - 1;                                 % Get number FOBs 
FOB.CAPACITY                    = 8000;                                         % Set custonmer max capacity
FOB.INITIAL_INVENTORY           = MapBlue.InitInventory'/2;                       % Get initial inventory levels
FOB.INITIAL_INVENTORY(NUM_BLUE) = 0;                                            % Returns depot demand to 0
FOB.THREAT                      = MapBlue.FOBThreatVec(1).ThreatVec';           % Get FOB threat vectors
FOB.SERVICE_TIMES				= zeros(1,NUM_BLUE);                            % Initialize service times
FOB.TIME_WINDOWS				= [zeros(1,NUM_BLUE);ones(1,NUM_BLUE)*period];  % Set time windows to entire period

VEHICLE.NUM_MIN                 = 1;                           % Set number of vehicles
VEHICLE.NUM_MAX                 = 10;                           % Set number of vehicles
VEHICLE.CAPACITY                = GeneralSettings.UAVCapacity; % Set Vehicle capacity
VEHICLE.DIST_LIM                = 494;
VEHICLE.SPEED                   = 148;

PARAMETERS.VFCS					= [0 20000 -2]; %[0 24000 -2];			% Set value function parameters [0 1 0];  [0 20000 -2]; 
PARAMETERS.PERIOD				= period;				% Length of single time period
PARAMETERS.SCALING_FACTOR		= 500;					% Unit-izes delivery increments
PARAMETERS.DISPLAY_LOOP			= Display_Loop;			% Display control message in command window
PARAMETERS.PLOT_VFS				= Plot_VFs;				% Plot value functions
PARAMETERS.NON_IMPROV_LIM		= 15;					% Num non-improving moves allowed for altering demand
PARAMETERS.FMINCONOPTIONS		= optimset('Display','off','Algorithm','active-set');
PARAMETERS.ValueFuncMax			= fmincon(@(x) value_function(x, PARAMETERS.VFCS, -1), 0, 0, 0, 0, 0, 0, FOB.CAPACITY, [], PARAMETERS.FMINCONOPTIONS);
PARAMETERS.ValueFuncMaxScaled	= round( PARAMETERS.ValueFuncMax/PARAMETERS.SCALING_FACTOR, 0);
PARAMETERS.VFCS_SCALED			= [	PARAMETERS.VFCS(1), ...
									PARAMETERS.VFCS(2)/PARAMETERS.SCALING_FACTOR, ...
									PARAMETERS.VFCS(2)/(-2 * PARAMETERS.SCALING_FACTOR * PARAMETERS.ValueFuncMaxScaled) ...
												*(PARAMETERS.VFCS(3) ~= 0) ...		%  ~= 0 => function is linear so theta(3) = 0
												*(-1)^(PARAMETERS.VFCS(3) > 0)];	%  > 0  => function is concave up

ACO_PARAMETERS.ALPHA            = 1;
ACO_PARAMETERS.BETA             = 3;
ACO_PARAMETERS.RUN_TIME         = 0;
ACO_PARAMETERS.PBEST            = 0.95;
ACO_PARAMETERS.Q0               = 0.9;
ACO_PARAMETERS.RHO              = 0.7;
ACO_PARAMETERS.GP_RATIO         = 0;
ACO_PARAMETERS.NUM_ANTS         = ceil(FOB.NUM/ant_scale);
ACO_PARAMETERS.PHI              = 0;
ACO_PARAMETERS.ITERATIONS       = 15;
ACO_PARAMETERS.CL_SIZE			= ceil(FOB.NUM/cl_scale);

HEURISTICS.DISTANCES			= [];
HEURISTICS.TIMES				= [];
HEURISTICS.LOG_PROB_SURVIVE		= [];
HEURISTICS.LOG_PROB_DESTRUCT	= [];

%% Plot value functions
if PARAMETERS.PLOT_VFS
	PlotValueFunction(FOB,PARAMETERS)
	pause(1)
end


Solution = Controlv2;

N		 = VEHICLE.NUM_MAX - VEHICLE.NUM_MIN + 1;
results	 = zeros(5,N);



%% Show results summary

if Show_Results_Table
	DmdVal = zeros(N,FOB.NUM);
	for i = 1:N
        results(1,i) = size(Solution.value_solution{i,1}.sbest,1);
        results(3,i) = Solution.value_solution{i,1}.TotValue;
        results(5,i) = sum(sum(Solution.value_solution{i,1}.Loads));
		X = Solution.demand(i,2:end);
		for j = 1:FOB.NUM
			DmdVal(i,j) = value_function(X(j), PARAMETERS.VFCS_SCALED, 1);
		end
	end
	results (2,:) = sum(DmdVal,2);
    results (4,:) = sum(Solution.demand,2)'*PARAMETERS.SCALING_FACTOR;
    array2table(results,'RowNames',{'# Vehicles','Demand Value','Soln Value','Cum Demand','Tot Delivered'})
end



%% Show plot of loads and values per iteration

if Plot_Values
    flrsqrt = floor( sqrt( N)); 
    clgsqrt = ceil( sqrt( N));
    prootN  = flrsqrt * clgsqrt;
    if prootN < N
        rows = clgsqrt;
        cols = rows;
    else
        rows = flrsqrt;
        cols = clgsqrt;
    end
    
    MinYL = min( cell2mat( Solution.valuei));
	MaxYL = max( cell2mat( Solution.values));
    if MinYL == MaxYL
        MinYL = MinYL - abs( 0 - MaxYL) * 0.05;
        MaxYL = MaxYL + abs( 0 - MaxYL) * 0.05;
    else
        MinYL = MinYL - 0.05 * ( MaxYL - MinYL);
        MaxYL = MaxYL + 0.05 * ( MaxYL - MinYL); 
    end
    
    MinYR = min( 0, min( cell2mat( Solution.loads)))  * PARAMETERS.SCALING_FACTOR;
	MaxYR = max( cell2mat( Solution.loads)) * PARAMETERS.SCALING_FACTOR; 
    MaxYR = MaxYR + 0.06 * ( MaxYR - MinYR); 
    
    
    %%
	figure; hold on
    for i = 1:N
        MaxX = length(Solution.values{i});
        X    = 1:1:MaxX;
        YL1  = Solution.values{i}; % best value so far
        YL2  = Solution.valuei{i}; % iteration value
        YR   = Solution.loads{i}*PARAMETERS.SCALING_FACTOR; % iteration loads
        subplot(rows,cols,i)
            yyaxis left
            plot(X,YL1,X,YL2,':')
            axis([-inf MaxX MinYL MaxYL])
            if mod(i,cols) == 1, ylabel('Route Value'), end
        
            yyaxis right
            plot(YR)
            axis([-inf MaxX MinYR MaxYR])
            if mod(i,cols) == 0, ylabel('Total Delivered'), end
        
            if i > (rows - 1)*cols, xlabel('Iteration'), end
			if i == N 
				legend('Best so far','Iteration value','Load','Location','southeast')
				legend('boxoff')
			end
    end
	hold off
	pause(1)
end


%%
disp(['Run Time ' num2str(Solution.run_time) ' sec'])

clearvars -except	GeneralSettings ...
					CodeRunTime ...
					MapGrid ...
					MapBlue ...
					MapRed ...
					Options ...
					Solution ...
                    solution ...
					Save_mat ...
                    VEHICLE         ...
                    FOB             ...
                    PARAMETERS      ...
                    ACO_PARAMETERS  ...
                    HEURISTICS

switch Save_mat
    case 1
        clear Save_mat 
        save('Solution10.mat', 'Solution', '-mat');
    case 2
        clear Save_mat
        save(uiputfile('*.mat','Save Workspace As'));
    otherwise
        clear Save_mat
        warning('Solution.mat not saved.')
end



%
%