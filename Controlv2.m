function full_solution = Controlv2
% DESCRIPTION: secondary function containing main ACO functions.
% --------------------------------------------------------------------------------------------------
% ----Outputs---------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
% .value_solution	Cell array with all output data; further detailed below.
%
% .demand			Customer demands.  Columns correspond to customers, i; rows correspond to number
%					of vehicles (i.e., for vehicle 1,...,N, customer j should receive the specified 
%					delivery quantity.
%
% .values			Table showing how the solution value changed as the algorithm proceeded. Rows 
%					correspond to # of vehicles.  So looking at row i, you can see what the value is
%					for each iteration of the algorithm as the demands are altered. This info should 
%					give you an idea of how many non-improving moves to allow.
%
% .runtime			Total run time of the routing portion of the algorithm in seconds.

% --------------------------------------------------------------------------------------------------
% ----Notes-----------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------- 
% Note 1: This scaling factor will automatically adjust everything. I.e., if you have an 8000 lb 
% limit but only want to consider changes in 500 lb increments, set the scaling factor to 500 and 
% the routing algorithm will then run on a scale of 0-16 in 1 unit increments, converting back to 
% the original scale at the end.
%
% Note 2: Define time as whatever (hours, minutes, days) but make sure that the veh_speed is 
% input as distance/period. If dist_limit is in miles and period is in hours, make sure speed is 
% miles per hour.
%
% Note 3: Specifies max num non-improving moves when altering demand.

CntlTime = tic;
global	DISTANCES       ...
        PROBABILITIES   ...
		HEURISTICS		...
		NUM_BLUE        ...
        VEHICLE         ...
		FOB             ...
		PARAMETERS      ...
		ACO_PARAMETERS

% Set basic parameters
sf					= PARAMETERS.SCALING_FACTOR; % Note 1
nimlim				= PARAMETERS.NON_IMPROV_LIM; % Note 3
Theta				= PARAMETERS.VFCS;
MaxV_scaled			= PARAMETERS.ValueFuncMaxScaled;
theta				= PARAMETERS.VFCS_SCALED;

% FOB Specs
d_dist              = MoveDepotData(DISTANCES);
p_dist              = MoveDepotData(PROBABILITIES);
threat              = MoveDepotData(FOB.THREAT);
current_inventory   = MoveDepotData(FOB.INITIAL_INVENTORY);
service				= MoveDepotData(FOB.SERVICE_TIMES); 
tw					= MoveDepotData(FOB.TIME_WINDOWS); 

% Vehicle specs
Vehicles			= VEHICLE.NUM_MIN:1:VEHICLE.NUM_MAX;
VehRange			= VEHICLE.NUM_MAX - VEHICLE.NUM_MIN + 1;
veh_speed			= VEHICLE.SPEED;

% Heuristic arrays
HEURISTICS.DISTANCES		 = d_dist;
HEURISTICS.TIMES			 = d_dist./veh_speed;
HEURISTICS.LOG_PROB_SURVIVE	 = -log(p_dist);
HEURISTICS.LOG_PROB_DESTRUCT = -log(1 - p_dist);

% Scale vehicle and FOB inventory capacities
Q_scaled            = VEHICLE.CAPACITY/sf;
C_scaled            = FOB.CAPACITY/sf;
current_inventory	= current_inventory./sf;


%% Loop 1
% Calculate the total value of the new solution.

%Initializations.
value_solution		= cell(VehRange,1);
values				= [];
valuei				= [];
loads				= [];
routes				= [];

for i = 1:VehRange 
	Veh = Vehicles(i);
    if PARAMETERS.DISPLAY_LOOP, disp(['Loop 1: Recalc Solution ',num2str(Veh),' Vehicles']), end
	
    % Define an initial demand. This defines a proportional demand for each customer who is below their max_value.
    demand(i,:)  = DefineInitialDemand(Veh,current_inventory,threat,MaxV_scaled,Q_scaled,C_scaled,NUM_BLUE); %#ok
	
    % Run ACO to get an initial solution.							
    ACO_solution = ACO_value(	Veh, demand(i,:), tw, service, Q_scaled, theta, current_inventory, threat, ...
								ACO_PARAMETERS, VEHICLE, HEURISTICS, NUM_BLUE);
	
	values{i}		= ACO_solution.TotValue; %#ok
    valuei{i}		= values{i};  %#ok
    loads{i}		= sum(sum(ACO_solution.Loads));   %#ok
    routes{i}{1}	= ACO_solution.sbest; %#ok
	
    % Assign the solution to the overall highest value solution.
    value_solution{i} = ACO_solution;
end



%% Loop 2
% Tweak Demands: For each number of vehicles, we want to start tweaking the demands to see if we can increase the overall value of the solution.

for i = 1:VehRange
	Veh = Vehicles(i);
    if PARAMETERS.DISPLAY_LOOP, disp(['Loop 2: Demand Tweaks ',num2str(Veh),' Vehicles']), end
    
	trigger		   = 0;								% Initialize while-loop terminator
    current_value  = value_solution{i}.TotValue;	% Store the value of the current best solution.
	
	order_increase = 2:1:NUM_BLUE;					% Build a list of all customers for demand increase.
    order_decrease = order_increase;				% Build a list of all customers for demand decrease.

    % Only choose a demand increase from those that aren't already full.
    % Remove customers who are full from the list of increase eligibles.
    for j = FOB.NUM:-1:1
        if (demand(i,order_increase(j)) + current_inventory(order_increase(j)))== C_scaled % >= C_scaled
            order_increase(j) = [];
        end
    end

    % Only choose a demand decrease from those that aren't already empty.
    % Remove customers who aren't receiving a delivery from the list of decrease eligibles.
    for j = size(order_decrease,2):-1:1
        if demand(i,order_decrease(j)) == 0
            order_decrease(j) = [];
        end
    end

    
    if isempty(order_increase)									% If no eligible increases, force a decrease.
        choice1 = 2;
    elseif isempty(order_decrease)								% If no eligible decreases, force an increase.
        choice1 = 1;
    elseif isempty(order_increase) && isempty(order_decrease)	% If no eligible increases or decreases, then we're done. Iterate the for loop.
        continue;
    else
        choice1 = randperm(2,1);
    end
    
    counter = 2;

    % Repeat the following until we've reached the non-improving iteration limit.
	while trigger < nimlim
		
		if trigger > 0
			%Build a list of all customers.
			order_increase = 2:1:NUM_BLUE;
			order_decrease = order_increase;
			
			% Only choose a demand increase from those that aren't already full.
			% Remove customers who are full from the list of increase eligibles.
            for j = size(order_increase,2):-1:1
                if (demand(i,order_increase(j)) + current_inventory(order_increase(j))) == C_scaled %>= C_scaled
                    order_increase(j) = [];
                end
            end
            
			% Only choose a demand decrease from those that aren't already empty.
			% Remove customers who aren't receiving a delivery from the list of decrease eligibles.
            for j = size(order_decrease,2):-1:1
                if demand(i,order_decrease(j)) == 0
                    order_decrease(j) = [];
                end
            end		
			
			if isempty(order_increase)									% If no eligible increases, force a decrease.	
				choice1 = 2;   
			elseif isempty(order_decrease)								% If no eligible decreases, force an increase.
				choice1 = 1;    
			elseif isempty(order_increase) && isempty(order_decrease)	% If no eligible increases or decreases, then iterate the loop.
				continue;       
			else														% o.w. randomly choose one
				choice1 = randperm(2,1); 
			end
		end
		
		temp_demand = demand(i,:); % Store the demand in a temporary array so we can alter it to see if it yields an improved solution.
        
		%Increase a demand.
		if choice1 == 1
			% Of the non-full demands, choose one.
			choice2								 = randperm(length(order_increase),1);
			temp_demand(order_increase(choice2)) = temp_demand(order_increase(choice2)) + 1;
			if (temp_demand(order_increase(choice2)) + current_inventory(order_increase(choice2))) > C_scaled
				trigger = trigger + 1;
				continue;
			end
		
		% Decrease a demand.
		elseif choice1 == 2
			% Of the non-full demands, choose one.
			choice2 = randperm(size(order_decrease,2),1);
			temp_demand(order_decrease(choice2)) = temp_demand(order_decrease(choice2)) - 1;
			if temp_demand(order_decrease(choice2)) < 0
				trigger = trigger + 1;
				continue;
			end
		end
		
		% Run another check. If a FOB is not allowing split deliveries, then do not allow its demand to exceed the capacity of a single vehicle, Q.
        if choice1 == 1
            if (threat(order_increase(choice2)) == 1) && (demand(i,order_increase(choice2)) > Q_scaled)
                demand(i,order_increase(choice2)) = Q_scaled;
            end
        end
        
		% Run ACO to get new solution.		
        ACO_solution2		= ACO_value(Veh, demand(i,:), tw, service, Q_scaled, theta, current_inventory, threat, ...
										ACO_PARAMETERS, VEHICLE, HEURISTICS, NUM_BLUE);
							
		% Store the value of the new solution.
		new_value			= ACO_solution2.TotValue;
		
		% Compare the new solution to the previous best solution.
		if (new_value > current_value)
			% Reset trigger to 0 if new solution is better and copy this new solution.
			trigger				= 0;
			value_solution{i}	= ACO_solution2;
			current_value		= new_value;
            demand(i,:) = temp_demand;
		else
			% Increment the trigger to show that we had a non-improving move.
			trigger = trigger + 1;
        end
        
        % Record iteration results
        values{i}(counter)   = current_value; %#ok
        valuei{i}(counter)   = new_value;
        loads{i}(counter)    = sum(sum(ACO_solution2.Loads) );
        routes{i}{counter}	 = ACO_solution2.sbest;
        counter				 = counter + 1;
	end
end

% Re-scale the demands and loads arrays back to the original unit scale using scaling factor.
for i = 1:VehRange
    value_solution{i}.Loads  = (value_solution{i}.Loads) * sf;
end

% Output the solution.
full_solution = struct(	'value_solution',	{value_solution}, ...
						'demand',			{demand}, ...
						'values',			{values}, ...
						'valuei',			{valuei}, ...
						'loads',			{loads}, ...
						'routes',			{routes}, ...
						'run_time',			{toc(CntlTime)});
end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
function NewData = MoveDepotData(Data)
[numR,numC]		 = size(Data);		% Determine dimensions of Data
if numR == numC
	temprow		 = Data(numR,:);	% Put depot row values into temp variable
	Data(numR,:) = [];				% Clear depot row from Data
	Data		 = [temprow;Data];	% Make depot row the first entry
	
	tempcol		 = Data(:,numC);	% Put depot column values into temp variable
	Data(:,numC) = [];				% Clear depot column from Data
	Data		 = [tempcol,Data];	% Make depot column the first entry
	NewData		 = Data;
elseif numC == 1	
	temprow		 = Data(numR,:);	% Put depot row values into temp variable
	Data(numR,:) = [];				% Clear depot row from Data
	Data		 = [temprow;Data];	% Make depot row the first entry
	NewData		 = Data;
else
	tempcol		 = Data(:,numC);	% Put depot row values into temp variable
	Data(:,numC) = [];				% Clear depot row from Data
	Data		 = [tempcol,Data];	% Make depot row the first entry
	NewData		 = Data;
end
end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
function demand = DefineInitialDemand(	AvailableVehicles, current_inventory, threat, ...
										MaxV_scaled, Q_scaled, C_scaled,NUM_BLUE)
% Define an initial demand. This defines a proportional demand for each customer who is below their max_value.
% global NUM_BLUE

demand = zeros(1,NUM_BLUE);

for cust = 2:NUM_BLUE
    if current_inventory(cust) < MaxV_scaled
        demand(cust) = MaxV_scaled - current_inventory(cust);
    end
end

demand = ceil(demand.*((AvailableVehicles*Q_scaled)/(sum(demand))));


for cust = 2:NUM_BLUE			% Check #1: Ensure no FOBs capacity is exceeded.
    if demand(cust) > (C_scaled - current_inventory(cust))
        demand(cust) = (C_scaled - current_inventory(cust));
    end
end

for cust = 2:NUM_BLUE			% Check #2: If a customer is not allowing split deliveries, limit it's demand to exceed  a single-vehicle's capacity.
    if (threat(cust) == 1) && (demand(cust) > Q_scaled)
        demand(cust) = Q_scaled;
    end
end
end

%
%