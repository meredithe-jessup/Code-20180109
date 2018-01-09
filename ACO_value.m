function full_solution = ACO_value(	num_vehicles, demand, tw, service, Q, ...
									theta, current_inventory, threat, ACO_PARAMETERS, VEHICLE, HEURISTICS, NUM_BLUE)
%This code is for a simple VRP with normal constraints. Inputs are:
%heuristic matrix (nxn), demand array (1xn), time window array (2xn),
%service time array (1xn), max vehicle load, heuristic importance,
%pheromone importance, the number of iterations, and allotted run time.
%(Note that only of the last two is necessary for the code to run, so put
%in a non-negative value for one you'd like to use and zero for the other.
%If you put in non-negative values for both, the algorithm will default to
%run time).

%Inputs: 
% heuristic = matrix of distances/costs between nodes (symmetric)
% demand = array of demands for each node 
% tw = "2xn" array with beginning (1st row) and end (2nd row) of time window 
% service = array of service times for each node 
% Q = max vehicle load 
% alpha = relative importance of heuristic information 
% beta = relative importance of pheromone information
% ro = pheromone evaporation rate (now calculated) 
% min = min pheromone level
% max = max pheromone level

tstart		= tic;
hx          = NUM_BLUE;
dist_limit  = VEHICLE.DIST_LIM;
heuristic1  = HEURISTICS.TIMES;
heuristic2  = HEURISTICS.DISTANCES;
heuristic3  = HEURISTICS.LOG_PROB_DESTRUCT;
p_dist2     = HEURISTICS.LOG_PROB_SURVIVE; 
alpha       = ACO_PARAMETERS.ALPHA;
beta        = ACO_PARAMETERS.BETA;
run_time    = ACO_PARAMETERS.RUN_TIME; 
pbest       = ACO_PARAMETERS.PBEST; 
qo          = ACO_PARAMETERS.Q0;        % acceptance criteria for city
rho         = ACO_PARAMETERS.RHO;
gp_ratio    = ACO_PARAMETERS.GP_RATIO; 
num_ants    = ACO_PARAMETERS.NUM_ANTS; 
phi         = ACO_PARAMETERS.PHI; 
iterations  = ACO_PARAMETERS.ITERATIONS;
cl_size		= ACO_PARAMETERS.CL_SIZE;

%Calculate rho by Stutzle's rule based on number of iterations.  If iteration count isn't given, estimate based on allotted run_time.
% if iterations>0
%     sbar = iterations;
% else
%     sbar = run_time*5;
% end
% rho = ((1-(.05)^(1/(hx-1)))/(((hx-1)/2-1)*(.05)^(1/(hx-1))))^(1/sbar);


%Form the candidate lists;
candidate_list		= zeros(hx,cl_size);
candidate_list(1,:) = 1:1:cl_size;

for i = 2:hx
    [~,temp] = sort(heuristic1(i,2:hx));
    candidate_list(i,1:cl_size) = temp(2:(cl_size + 1)) + 1;
end

%Initialize global best solution, cost, and array of all costs
sbest								= [];
antcostcounter						= 0;
antcosts							= zeros(1,max((ceil(sqrt(hx-1))*iterations),1));
loads_bestofiterationroute			= [];
loads_best							= [];
start_times_bestofiterationroute	= [];
start_times_best					= [];
route_cost_bestofiterationroute		= [];
route_value_best					= [];
s_value_best						= [];
s_value_bestofiteration				= [];
route_value_bestofiteration			= [];
route_dist_bestofiterationroute		= [];
all_value_bestofiteration			= [];


%Run a single iteration with a beta of zero, which is equivalent to a probabilistic nearest neighbor solution.
ant1		= ant_value(	heuristic1, heuristic2, candidate_list, cl_size, ones(hx), demand, tw, service, Q, ...
							alpha, 0, 0, 1, hx, qo, phi, dist_limit, threat, theta, current_inventory);
						
mx			= 1/(1-rho)*1/((hx-1)*ant1.cost);	%Use this initial solution to estimate the maximum pheromone amount.
pheromone	= ones(hx).*mx;						%Form pheromone matrix, setting all pheromones to max value.

%Using Stutzle's rule, construct min
pdec		= pbest^(1/(hx-2));
mn			= (2*mx*(1-pdec))/((hx-1)*pdec);

%Define my limit as either time or iterations
if run_time > 0
    limit	= run_time;
    counter = toc;
else
    limit	= iterations;
    counter = 0;
end

iteration_counter = 0;
pheromone_counter = 0;

%Main ACO function
while counter < limit
    
    iteration_counter = iteration_counter + 1;
    pheromone_counter = pheromone_counter + 1;
    
    %Set iteration best cost to big M, in this case 10x the probabilistic nearest neighbor solution calculated above.
    bestofiterationcost     = ant1.cost*10;
    bestofiterationsolution = [];
    
    %# of ants is equal to the number of customers
	for i = 1:num_ants
        
        %Subrouting constructs an entire solution for 1 ant
		ant = ant_value(	heuristic1, heuristic2, candidate_list, cl_size, pheromone, demand, tw, service, Q, ...
										alpha, beta, mn, mx, hx, qo, phi, dist_limit, threat, theta, current_inventory);
    
        %Calculate the values of the given routes.
        all_value = route_values(ant.s, ant.loads, current_inventory, theta);
        if num_vehicles < size(all_value,1)
            route_value = all_value(1:num_vehicles,:);
        else
            route_value = all_value;
        end

        if size(ant.s,2) > size(ant.loads,2)
            temp		= zeros(size(ant.loads,1),(size(ant.s,2)-size(ant.loads,2)));
            ant.loads	= [ant.loads,temp];
        end
        
        if size(route_value,1) < num_vehicles
            temp		= zeros(num_vehicles - size(route_value,1),size(route_value,2));
            route_value = [route_value;temp]; %#ok
        end
        
        %Calculate the total value of a solution.
        s_value = sum(sum(route_value));
        
        %Calculate probabilities.
        %Initialize to 0, except first column to 1 (we're assuming the base doesn't get bombed, so probability of survival is 1 prior to departure).
        prob	= zeros(size(route_cost_bestofiterationroute));
		
        %Compute route survival probabilities.
		[S1,S2] = size(ant.s);
        for i31 = 1:S1
            for i32 = 2:S2
                %Check to see if we're at the end of a route.
                if ant.s(i31,i32) == 0
                    break;
				end
				
                %If not, add probability of survival of reaching that node (i.e., probability of survival associated with previous leg).
                prob(i31,i32) = heuristic3(ant.s(i31,i32-1),ant.s(i31,i32));
            end
        end
        
        %If the first run, then store the solution as the best.
        if i == 1
            s_value_bestofiteration				= s_value;
            route_value_bestofiteration			= route_value;
            bestofiterationcost					= ant.cost;
            bestofiterationsolution				= ant.s;
            start_times_bestofiterationroute	= ant.start_times;
            loads_bestofiterationroute			= ant.loads;
            route_cost_bestofiterationroute		= ant.route_cost;
            route_dist_bestofiterationroute		= ant.route_dist;
            all_value_bestofiteration			= all_value;

            if counter == 0
                s_value_best		= s_value;
                route_value_best	= route_value;
                costbest			= ant.cost;
                sbest				= ant.s;
                start_times_best	= ant.start_times;
                loads_best			= ant.loads;
                route_dist_best		= ant.route_dist;
                prob_best			= prob;
            end
        end
            
       
        %Update global best
        if s_value > s_value_best
            s_value_best		= s_value;
            route_value_best	= route_value;
            costbest			= ant.cost;
            sbest				= ant.s;
            start_times_best	= ant.start_times;
            loads_best			= ant.loads;
            route_dist_best		= ant.route_dist;
            prob_best			= prob;
        end

        %Store current cost in array
        antcostcounter				= antcostcounter + 1;
        antcosts(antcostcounter)	= ant.cost;
        
        if s_value > s_value_bestofiteration
            s_value_bestofiteration				= s_value;
            route_value_bestofiteration			= route_value;
            bestofiterationcost					= ant.cost;
            bestofiterationsolution				= ant.s;
            start_times_bestofiterationroute	= ant.start_times;
            loads_bestofiterationroute			= ant.loads;
            route_cost_bestofiterationroute		= ant.route_cost;
            route_dist_bestofiterationroute		= ant.route_dist;
            all_value_bestofiteration			= all_value;
        end
        
	end

   
%% Local Search based on value function

    feeder_solution     = bestofiterationsolution;
    feeder_cost         = bestofiterationcost;
    feeder_route_cost   = route_cost_bestofiterationroute;
    feeder_start_times  = start_times_bestofiterationroute;
    feeder_loads        = loads_bestofiterationroute;
    feeder_route_dist   = route_dist_bestofiterationroute;
    feeder_s_value      = s_value_bestofiteration;
    feeder_route_value  = route_value_bestofiteration;
    feeder_all_value    = all_value_bestofiteration;
    lowest_num_veh      = size(feeder_solution,1);
	
    change3 = 1;
    trigger = 1;

	while trigger == 1
        
        trigger = 0;
		
		% change1 = 0 temporarily disable 2-opt*
        vale = min(num_vehicles, size(feeder_solution,1));
        if (vale > 1) && (size(feeder_solution,1) > num_vehicles)
            change1 = 1;
            change2 = 1;
			
        else %(size(feeder_solution,1)-num_vehicles) <= 0
            change1 = 0;
            change2 = 0;        
        end

		while change1 == 1
			
			if (size(feeder_solution,1) - num_vehicles) <= 0
				break;
			end
			
			[	feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, ...
				feeder_route_dist, feeder_s_value, feeder_route_value, feeder_all_value, change1] ...
							= two_opt_star_value(	feeder_solution, feeder_cost, feeder_route_cost, heuristic1, heuristic2,  tw, service, Q, ...
													feeder_start_times, feeder_loads, dist_limit, feeder_route_dist, num_vehicles, feeder_s_value, ...
													feeder_route_value, feeder_all_value, current_inventory, theta);
			if change1 == 1
				trigger = 1;
				change2 = 1;
				change3 = 1;
			end
		end
		
		% change2 = 0 temporarily disable Or-opt        
        while change2 == 1

            if (size(feeder_solution,1) - num_vehicles) <= 0
                break;
            end
            
            [	feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, ...
				feeder_route_dist, feeder_s_value, feeder_route_value, feeder_all_value, change2] ...
						= Or_opt_v2_value(		feeder_solution, feeder_cost, feeder_route_cost, heuristic1, heuristic2, tw, service, Q, ...
												feeder_start_times, feeder_loads, dist_limit, num_vehicles, feeder_route_dist, feeder_s_value, ...
												feeder_route_value, feeder_all_value);
            if change2 == 1
                trigger = 1;
                change1 = 1;
                change3 = 1;
            end
        end

        while change3 == 1

            [	feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, ...
				feeder_route_value, feeder_all_value, change3] ...
						= relocate_value(		feeder_solution, feeder_cost, feeder_route_cost, heuristic1, tw, service, Q, ...
												feeder_start_times, feeder_loads, feeder_route_value, feeder_all_value);
            if change3 == 1
                if num_vehicles > 1
                    trigger = 1;
                    change1 = 1;
                    change2 = 1;
                end
            end
        end
	end
    
    
    %Update global best
    if (feeder_s_value > s_value_best) || ((feeder_s_value == s_value_best) && (size(feeder_all_value,1) < lowest_num_veh))
        s_value_best        = feeder_s_value;
        route_value_best    = feeder_route_value;
        costbest            = feeder_cost;
        sbest               = feeder_solution;
        start_times_best    = feeder_start_times;
        loads_best          = feeder_loads;
        route_dist_best     = feeder_route_dist;
    end
    

%% Local Search to optimize distances, run only on the appropriate number of routes.

    %Trim and recalculate as appropriate.
	feeder_solution		= TrimFeeder(num_vehicles, feeder_solution);
	feeder_route_cost	= TrimFeeder(num_vehicles, feeder_route_cost);
	feeder_start_times	= TrimFeeder(num_vehicles, feeder_start_times);
	feeder_loads		= TrimFeeder(num_vehicles, feeder_loads);
	feeder_route_dist	= TrimFeeder(num_vehicles, feeder_route_dist);
	feeder_cost         = sum(feeder_route_cost);
	
    % Initialize to 0, then initialize the first column to 1. We're assume  base doesn't get bombed, so route initial probability of survival is 1.
    feeder_prob = zeros(size(feeder_route_cost));

    %Compute route survival probabilities.
    for i31 = 1:size(feeder_solution,1)
        for i32 = 2:size(feeder_solution,2)
            %Check to see if we're at the end of a route.
            if feeder_solution(i31,i32) == 0
                break;
            end
            
            %If not, add probability of survival of reaching that node (i.e., probability of survival associated with previous leg).
            feeder_prob(i31,i32) = heuristic3(feeder_solution(i31,i32-1),feeder_solution(i31,i32));
        end
    end

    change1 = 1;
    change2 = 1;
    change3 = 1;
    trigger = 1;

    while trigger == 1
 
        trigger = 0;
      
        while change1 == 1
            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, feeder_route_dist, change1] ...
						= two_opt_star(		feeder_solution, feeder_cost, feeder_prob, heuristic1, heuristic2, heuristic3, tw, service, Q, ...
											feeder_start_times, feeder_loads, dist_limit, feeder_route_dist);
            if change1 == 1
                trigger = 1;
                change2 = 1;
                change3 = 1;
            end
        end
        
% change2=0;%temporarily disable Or-opt
        while change2 == 1
            
            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, feeder_route_dist, change2] ...
						= Or_opt_v2(		feeder_solution, feeder_cost, feeder_prob, heuristic1, heuristic2, heuristic3, tw, service, Q, ...
											feeder_start_times, feeder_loads, dist_limit, feeder_route_dist);
            if change2 == 1
                trigger = 1;
                change1 = 1;
                change3 = 1;
            end
        end

        while change3 == 1

            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, change3] ...
						= relocate(			feeder_solution, feeder_cost, feeder_prob, heuristic1, heuristic3, tw, service, Q, ...
											feeder_start_times, feeder_loads);
            if change3 == 1
                trigger = 1;
                change1 = 1;
                change2 = 1;
            end
        end
    end

    % Recalculate the route_value. Total value is unchanged but customers (and their associated values) were shifted around.
    route_value     = route_values(feeder_solution, feeder_loads, current_inventory, theta);
    vale            = min(size(route_value,1),num_vehicles);
    route_value     = route_value(1:vale,:);
    
    ls_solution     = feeder_solution;
    ls_start_times  = feeder_start_times;
    ls_loads        = feeder_loads;
    ls_dist         = feeder_route_dist;
    ls_value        = feeder_s_value;
    ls_route_value  = route_value;
    ls_prob         = feeder_prob;
    
%__________________________________________________________________________
    
    %Update global best
	route_dist_best     = ls_dist;
	sbest               = ls_solution;
	start_times_best    = ls_start_times;
	loads_best          = ls_loads;
	s_value_best        = ls_value;
	route_value_best    = ls_route_value;
	prob_best           = ls_prob;
    
    %Local pheromone update.  Decrease all by rho.
    pheromone = pheromone.*(1-rho);
    
    %Check against min values.
    for x1 = size(pheromone,1)
        for x2 = size(pheromone,2)
            if pheromone(x1,x2) < mn
                pheromone(x1,x2) = mn;
            end
        end
    end
    
    %Global pheromone update using global best
    for k=1:size(sbest,1)-1
        for j=1:size(sbest,2)-1
            
            %Only update if the entry is an actual path in a route
            if (sbest(k,j) ~= 0) && (sbest(k,j+1) ~= 0)
                pheromone(sbest(k,j),sbest(k,j+1)) = pheromone(sbest(k,j),sbest(k,j+1)) + mx/3*gp_ratio;%(1/costbest)*gp_ratio;
                
                %Check pheromones against min/max values.  Only need to check against max because we're just adding positive pheromone at this step.
                if pheromone(sbest(k,j),sbest(k,j+1)) > mx
                    pheromone(sbest(k,j),sbest(k,j+1)) = mx;
                end
               
                %Maintain symmetry
                pheromone(sbest(k,j+1),sbest(k,j)) = pheromone(sbest(k,j),sbest(k,j+1));
                
            end
        end
    
    end
    
    %Global pheromone update using iteration best
    for k=1:size(ls_solution,1)-1
        for j=1:size(ls_solution,2)-1
            
            %Only update if the entry is an actual path in a route
            if (ls_solution(k,j) ~= 0) && (ls_solution(k,j+1) ~= 0)
                pheromone(ls_solution(k,j),ls_solution(k,j+1)) = pheromone(ls_solution(k,j),ls_solution(k,j+1)) + mx/3*gp_ratio;
                
                %Check pheromones against min/max values
                if pheromone(ls_solution(k,j),ls_solution(k,j+1)) > mx
                    pheromone(ls_solution(k,j),ls_solution(k,j+1)) = mx;
                end
 
                %Maintain symmetry
                pheromone(ls_solution(k,j+1),ls_solution(k,j)) = pheromone(ls_solution(k,j),ls_solution(k,j+1));

            end
        end
    end
    
    
    if run_time > 0
        counter = toc;
    else
        counter=counter+1;
    end
end

%Fix sbest to match Ian's numbering (depot is last, not first).
for j = 1:size(sbest,2)
    for i = 1:size(sbest,1)
        if sbest(i,j) == 1
            sbest(i,j) = hx; %#ok
        elseif sbest(i,j) > 1
            sbest(i,j) = sbest(i,j) - 1; %#ok
        end
    end
end

%Compute the cumulative probabilities.
cumul_prob = zeros(size(sbest));
for i = 1:size(sbest,1)
    for j = 2:size(sbest,2)
        if sbest(i,j) == 0
            break;
        end
        
        cumul_prob(i,j) = cumul_prob(i,j-1) + p_dist2(sbest(i,j-1),sbest(i,j));
    end
end

%Convert probabilities back from log format to actual probabilities.
prob_best		= 1-exp(-prob_best);
prob_best(:,1)	= 1;
cumul_prob		= exp(-cumul_prob);


%Trim all of the parts of the solution to make sure they're the same size, meaning delete any dummy columns.

%Sbest
x = size(sbest,2);
y = sum(sbest,1);
for k = x:-1:1
    if y(k) == 0
        sbest(:,k) = [];
    else
        break;
    end
end

%Start_times
y = size(start_times_best,2);
if y > x
    for k = y:-1:x+1
        start_times_best(:,k) = [];
    end
end

%loads
y = size(loads_best,2);
if y > x
    for k = y:-1:x+1
        loads_best(:,k) = [];
    end
end

%route_value
y = size(route_value_best,2);
if y > x
    for k = y:-1:x+1
        route_value_best(:,k) = [];
    end
end

%prob
y = size(prob_best,2);
if y > x
    for k = y:-1:x+1
        prob_best(:,k) = [];
    end
end

%cumul_prob
y = size(cumul_prob,2);
if y > x
    for k = y:-1:x+1
        cumul_prob(:,k) = [];
    end
end

%Output solution
full_solution = struct(	'sbest',		{sbest},					...	%Best solution
						'TotDist',		{sum(sum(ant.route_dist))},	...	%This is the total distance required to fulfill the solution.
						'Start_Times',	{start_times_best},			...	%Total time for each route. The last row value > 0 is vehicle's RTB time.
						'Loads',		{loads_best},				...	%This is the load amount for each delivery.
						'Run_Time1',	{toc(tstart)},				...	%This is the total time it took to generate all of this.
						'Route_Dists',	{route_dist_best},			...	%This is the distance covered by each route.
						'TotValue',		{s_value_best},				...	%This is the total value added by executing this solution.
						'Route_Values',	{route_value_best},			...	%These are the individual values of each delivery.
						'Prob_Surv',	{prob_best},				...	%These are the probabilities of survival for each leg of a route.
						'Cumu_Prob',	{cumul_prob});					%These are the cumulative probabilites for each route.
end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
function [solution, cost, route_cost, start_times, loads, route_dist,	change] ...
					= Or_opt_v2(			solution, cost, route_cost, heuristic, heuristic2,	heuristic3,	tw, service, Q,	start_times, loads, ...
										dist_limit,					route_dist)
%Given an initial solution, this program will run Or-opt to try to find an improved solution. This differs from Or-opt function in that this one will
%try to move larger segments first, iteratively decreasing the size until it finds a feasible improving move.

%Initialization
change = 0;

%Randomly re-order an index of routes.  This will tell us what order we're going to look at routes.
routes		= transpose(randperm(size(solution,1)));
num_routes	= size(routes,1);

%For each route, see if Or-opt* will yield an improved solution. 
%This loop controls the first route.
for i = 1:(size(solution,1)-1)
    
    %Assign the first route.
    route1 = routes(i);
    
    %Check each arc of the route in question, except the depot return arc.
    %This loop controls the arc from the first route.
    customer1 = 1;
    while solution(route1,customer1+1) > 0
        
        %Select a route with which to perform the switch. Note that we will only run 2-opt* on the routes occuring after the current route in
        %this array, avoiding unnecessary double work. This loop controls the second route.
        for j = (i+1):num_routes
			
            route2		= routes(j); %Assign the second route.
            
            %Select first arc from second route. Don't select the depot return route because then we have no more arcs from which to
            %select the second arc. This loop contols the first arc from the second route.
            customer21	= 1;
			
            while solution(route2,customer21+1) > 1

%--------------------------------------------------------------------------
                %Initialize customer22 to the last customer prior to the depot
                customer22 = customer21 + 1;
                for check = (customer22+1):(size(solution,2))
                    if solution(route2,check) == 1
                        break;
                    end
                    customer22 = check;
                end
                
                %Need to select an additional arc from the second route. This will tell us which sequence of customers we're going to attempt to move
				%to route1. This loop contols the second arc from the second route.
                while customer22 > customer21
                    
                    insertion_size = customer22 - customer21;
                    
%--------------------------------------------------------------------------
                    %Now we have a valid arc from the first route, which indicates the insertion point, and two valid arcs from the second route, 
					%which indicates the sequence of customers that we'll move. Build the first new route.
                    temp1				= [solution(route1,:)	,	zeros(1,insertion_size)];
                    temp1_loads			= [loads(route1,:)		,	zeros(1,insertion_size)];
                    temp1_start_times	= [start_times(route1,:),	zeros(1,insertion_size)];

                    %Prior to customer1, route1 is unchanged.  So zero out everything after that.  Reason for zeroing out is in case route ends up
					%shorter than the original, we don't have any vestigial customers from the original route.
                    for k = (customer1+1):(size(temp1,2))
                        temp1(1,k)				= 0;
                        temp1_loads(1,k)		= 0;
                        temp1_start_times(1,k)	= 0;
                    end

                    %Now, add in the piece of the second route.
                    for k = 1:insertion_size
                        %Add the route2 customer to route1.
                        temp1(customer1 + k) = solution(route2,customer21 + k);

                        %Add the route2 customer's to route1 loads.
                        temp1_loads(customer1 + k) = loads(route2,customer21 + k);

                        %Add the route2 customer's start time to route1 start  times.
						serv1 = service(temp1(customer1 + k - 1));
						heur1 = heuristic(temp1(customer1 + k - 1),temp1(customer1 + k));
                        temp1_start_times(customer1 + k) = temp1_start_times(customer1 + k - 1) + serv1 + heur1;
						
                        %If we arrive before the opening of a time window, don't start service until the time window opens.
                        if temp1_start_times(customer1 + k) < tw(1,temp1(customer1 + k))
                            temp1_start_times(customer1 + k) = tw(1,temp1(customer1 + k));
                        end
                    end
                    
                    %Now add the end of the first route back.
                    k = 1;
                    while solution(route1,customer1 + k) > 0
                        
                        %Add the remainder of the route1 customers.
                        temp1(customer1 + insertion_size + k) = solution(route1,customer1 + k);
                        
                        %Add the remainder of the route1 customers' loads.
                        temp1_loads(customer1 + insertion_size + k) = loads(route1,customer1 + k);

                        %Add the remainder of the route1 customers' start times.
						serv2 = service(temp1(customer1 + insertion_size + k - 1));
						heur2 = heuristic(temp1(customer1 + insertion_size + k - 1),temp1(customer1 + insertion_size + k));
                        temp1_start_times(customer1 + insertion_size + k) = temp1_start_times(customer1 + insertion_size + k - 1) + serv2 + heur2;
						
                        %If we arrive before the opening of a time window, don't start service until the time window opens.
                        if temp1_start_times(customer1 + insertion_size + k) < tw(1,temp1(customer1 + insertion_size + k))
                            temp1_start_times(customer1 + insertion_size + k) = tw(1,temp1(customer1 + insertion_size + k));
                        end
                        
                        k = k+1;
                        %If the increment puts us beyond the end of the solution array, end the loop.
                        if (customer1 + k) > size(solution,2)
                            break;
                        end
                    end


%--------------------------------------------------------------------------
                    %First new route is now built, so start check load and time windows for feasibility.

                    %Check load of new route. If it exceeds vehicle capacity, then iterate to the next customer22 of route2.
                    if sum(temp1_loads) > Q
                        customer22 = customer22 - 1;
                        continue;
                    end

                    %Check TWs of new route. If any deadlines violated, move to new arc in route. Start at k=2 b/c don't need to check depot times.
                    breaker = 0;
                    k = 2;
                    while temp1(k) > 1
                        if temp1_start_times(k) > tw(2,temp1(k))
                            breaker = 1;
                            break;
                        end
                        k=k+1;
                    end

                    %If breaker == 1, then this route busts a time window. So iterate to the next customer22.
                    if breaker == 1
                        customer22 = customer22 - 1;
                        continue;
                    end

    %--------------------------------------------------------------------------
                    %Build the second new route.
                    temp2 = solution(route2,:);
                    temp2_loads = loads(route2,:);
                    temp2_start_times = start_times(route2,:);

                    %Prior to customer21, the route is unchanged.  Zero out everything after that.  Reason for zeroing out is in case the route ends 
					%up shorter than the original, we don't have any vestigial customers from the original route.
                    for k = (customer21+1):(size(temp2,2))
                        temp2(1,k) = 0;
                        temp2_loads(1,k) = 0;
                        temp2_start_times(1,k) = 0;
                    end

                    %Now add shift everything from customer22+1 to  position starting at customer21+1 to account for the sequence moved to route1.
                    k = 1;
                    while solution(route2,customer22 + k) > 0
                        
                        %Add the remainder of the route2 customers.
                        temp2(customer21 + k) = solution(route2,customer22 + k);
                        
                        %Add the remainder of the route1 customers' loads.
                        temp2_loads(customer21 + k) = loads(route2,customer22 + k);

                        %Add the remainder of the route1 customers' start times.
						serv3 = service(temp2(customer21 + k - 1));
						heur3 = heuristic(temp2(customer21 + k - 1),temp2(customer21 + k));
                        temp2_start_times(customer21 + k) = temp2_start_times(customer21 + k - 1) + serv3 + heur3;
						
                        %If we arrive before the opening of a time window, don't start service until the time window opens.
                        if temp2_start_times(customer21 + k) < tw(1,temp2(customer21 + k))
                            temp2_start_times(customer21 + k) = tw(1,temp2(customer21 + k));
                        end
                        
                        k = k+1;
                        %If the increment puts us beyond the end of the solution array, end the loop.
                        if (customer22 + k) > size(solution,2)
                            break;
                        end
                    end

                    %Don't need to check load of route2 because all we've done is remove loads.

                    %Don't need to check the time windows of the route2 because we've shifted everything to earlier times (or remained the same).

    %--------------------------------------------------------------------------
                    %New solution is feasible wrt to both vehicle loads and time windows, so check for improvement.
                    
                    %First, check route distances against distance limit. Compute the dist of the first route.
                    k = 1;
                    temp1_dist = 0;
                    while temp1(k+1) > 0
                        temp1_dist = temp1_dist + heuristic2(temp1(k),temp1(k+1));
                        k = k+1;
                        if k+1 > size(temp1,2)
                            break;
                        end
                    end

                    %Compute the dist of the second route.
                    k = 1;
                    temp2_dist = 0;
                    while temp2(k+1) > 0
                        temp2_dist = temp2_dist + heuristic2(temp2(k),temp2(k+1));
                        k = k+1;
                        if k+1 > size(temp2,2)
                            break;
                        end
                    end

                    %If distance limit is busted, then don't implement. Move on to next arc.
                    if temp1_dist > dist_limit || temp2_dist > dist_limit
                        customer22 = customer22 - 1;
                        continue;
                    end                    

                    %Compute the cost of the first route.
                    k = 1;
                    temp1_cost = 0;
                    while temp1(k+1) > 0
                        temp1_cost = temp1_cost + heuristic3(temp1(k),temp1(k+1));
                        k = k+1;
                        if k+1 > size(temp1,2)
                            break;
                        end
                    end

                    %Compute the cost of the second route.
                    k = 1;
                    temp2_cost = 0;
                    while temp2(k+1) > 0
                        temp2_cost = temp2_cost + heuristic3(temp2(k),temp2(k+1));
                        k = k+1;
                        if k+1 > size(temp2,2)
                            break;
                        end
                    end

                    %Update the costs associated with the solution.
                    new_cost = cost - route_cost(route1) - route_cost(route2) + temp1_cost + temp2_cost;

                    %If no cost improvement, then don't implement. Move on to next arc.
                    if new_cost >= cost
                        customer22 = customer22 - 1;
                        continue;
                    end

    %--------------------------------------------------------------------------
                    %At this point, we've found a cheaper solution that is feasible wrt to both vehicle loads and time windows. check solution for
					%improvement. implement the solution and return to the calling function.

                    %This is a binary indicator. Set to 1 because we have found a better solution.
                    change = 1;

                    %Update the costs associated with the solution.
                    cost = new_cost;
                    route_cost(route1) = temp1_cost;
                    route_cost(route2) = temp2_cost;
                    
                    %Update route distances.
                    route_dist(route1) = temp1_dist;
                    route_dist(route2) = temp2_dist;
          
                    %Update route1 in the solution.
                    ssize = size(solution,2);
                    for k = 1:(size(temp1,2))
                        solution(route1,k) = temp1(k);
                    end
                    if size(temp1,2) < ssize
                        for k = (size(temp1,2)+1):ssize
                            solution(route1,k) = 0;
                        end
                    end

                    %Update route2 in the solution.
                    for k = 1:(size(temp2,2))
                        solution(route2,k) = temp2(k);
                    end
                    if size(temp2,2) < ssize
                        for k = (size(temp2,2)+1):ssize
                            solution(route2,k) = 0;
                        end
                    end

                    %Update route1 start_times.
                    for k = 1:(size(temp1_start_times,2))
                        start_times(route1,k) = temp1_start_times(k);
                    end
                    if size(temp1_start_times,2) < ssize
                        for k = (size(temp1_start_times,2)+1):ssize
                            solution(route1,k) = 0;
                        end
                    end

                    %Update route2 start_times.
                    for k = 1:(size(temp2_start_times,2))
                        start_times(route2,k) = temp2_start_times(k);
                    end

                    %Update route1 loads.
                    for k = 1:(size(temp1_loads,2))
                        loads(route1,k) = temp1_loads(k);
                    end

                    %Update route2 loads.
                    for k = 1:(size(temp2_loads,2))
                        loads(route2,k) = temp2_loads(k);
                    end

                    %If either of the routes is shorter than the one it's replacing, zero out any vestigial customers/loads/start times.
                    if size(temp1,2) < ssize
                        for k = (size(temp1,2)+1):ssize
                            solution(route1,k) = 0;
                            loads(route1,k) = 0;
                            start_times(route1,k) = 0;
                        end
                    end
                    
                    if size(temp2,2) < ssize
                        for k = (size(temp2,2)+1):ssize
                            solution(route2,k) = 0;
                            loads(route2,k) = 0;
                            start_times(route2,k) = 0;
                        end
                    end
                    
                    %Check the second route and if it's empty (no customers), then delete it.
                    if solution(route2,2) == 1
                        solution(route2,:) = [];
                        loads(route2,:) = [];
                        start_times(route2,:) = [];
                        route_cost(route2) = [];
                        route_dist(route2) = [];
                    end

                    %Run the fix_routes subroutine to make sure we aren't visiting the same customer twice within a single route.
                    [solution, cost, route_cost, start_times, loads, route_dist] ...
						= fix_routes(solution, cost, route_cost, heuristic, heuristic2, heuristic3, tw, service, start_times, loads, route_dist);

                    %Delete any zero columns from the solution, plus the corresponding columns from the loads and start times arrays.
                    index = size(solution,2);
                    for k = 0:(index-1)
                        y = sum(solution,1);
                        if y(index - k) == 0
                            solution(:,(index - k)) = [];
                            loads(:,(index - k)) = [];
                            start_times(:,(index - k)) = [];
                        else
                            break;
                        end
                    end
                    
                    %We've found an improved solution. Return to the calling function.
                    return;
                end
                
                %If we reach this line, then we've either exhausted  customer22 possibilities or we've ruled out any solutions using the current 
				%customer1/customer21 pairing. In either case, move on to the next customer21.
                customer21 = customer21 + 1;
                
            end
        end
        
        %Go to the next customer in the first route.
        customer1 = customer1 + 1;
        if (customer1+1) > size(solution,2)
            break
        end
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------

function [solution, cost, route_cost, start_times, loads, route_dist, value, route_value, all_value, change] ...
					= Or_opt_v2_value(		solution,	cost, route_cost, heuristic, heuristic2,			tw, service, Q, start_times, loads, ...
										dist_limit, num_vehicles,	route_dist, value, route_value, all_value)
%Given an initial solution, this program will run Or-opt to try to find an improved solution. This differs from Or-opt function in that this one will
%try to move larger segments first, iteratively decreasing the size until it finds a feasible improving move.

%Initialization
change = 0;

if size(route_value,1) < num_vehicles
    temp		= zeros(num_vehicles - size(route_value,1),size(route_value,2));
    route_value = [route_value;temp];
end

%Randomly re-order an index of routes.  This will tell us what order we're going to look at routes.
good_routes = transpose(randperm(num_vehicles,1));

bad_routes = transpose(randperm(size(solution,1)-num_vehicles,1));

%For each route, see if 2-opt* will yield an improved solution. 
%This loop controls the first route.
for i = 1:size(good_routes)
    
    %Assign the first route.
    route1 = good_routes(i);
    
    %Check each arc of the route in question, except the depot return arc.
    %This loop controls the arc from the first route.
    customer1 = 1;
    while solution(route1,customer1+1) > 0
        
        %Select a route with which to perform the switch. Note that we will only run 2-opt* on the routes occuring after the current route in
        %this array, avoiding unnecessary double work. This loop controls the second route.
        for j = 1:size(bad_routes)

            %Assign the second route.
            route2 = bad_routes(j)+num_vehicles;
            
            %Select first arc from second route. Don't select depot return route because then we have no more arcs from which to select second arc.
            %This loop contols the first arc from the second route.
            customer21 = 1;
            while solution(route2,customer21+1) > 1

%--------------------------------------------------------------------------
                %Initialize customer22 to the last customer prior to the depot
                customer22 = customer21 + 1;
                for check = (customer22+1):(size(solution,2))
                    if solution(route2,check) == 1
                        break;
                    end
                    customer22 = check;
                end
                
                %Need to select additional arc from second route. This determines which customer sequence we're going to attempt to move to route1.
                %This loop contols the second arc from the second route.
                while customer22 > customer21
                    
                    insertion_size = customer22 - customer21;
                    
%--------------------------------------------------------------------------
                    %Now we have a valid arc from the first route, which indicates the insertion point, and two valid arcs from the second route, 
					%which indicates the sequence of customers that we'll move.
					
                    %Build the first new route.
                    temp1				= [solution(route1,:),zeros(1,insertion_size)];
                    temp1_loads			= [loads(route1,:),zeros(1,insertion_size)];
					
%                    all_value
                    temp1_rv			= [all_value(route1,:),zeros(1,insertion_size)];
                    temp1_start_times	= [start_times(route1,:),zeros(1,insertion_size)];

                    %Prior to customer1, route1 is unchanged.  So zero out everything after that.  Reason for zeroing out is in case the route ends 
					%up shorter than the original, we don't have any vestigial customers from the original route.
                    for k = (customer1+1):(size(temp1,2))
                        temp1(1,k)				= 0;
                        temp1_loads(1,k)		= 0;
                        temp1_rv(1,k)			= 0;
                        temp1_start_times(1,k)	= 0;
                    end

                    %Now, add in the piece of the second route.
                    for k = 1:insertion_size
                        %Add the route2 customer to route1.
                        temp1(customer1 + k)		= solution(route2,customer21 + k);

                        %Add the route2 customer's to route1 loads.
                        temp1_loads(customer1 + k)	= loads(route2,customer21 + k);

                        %Add the route2 customer's to route1 values.
                        temp1_rv(customer1 + k)		= all_value(route2,customer21 + k);

                        %Add the route2 customer's start time to route1 start times.
						serv1 = service(temp1(customer1 + k - 1));
						heur1 = heuristic(temp1(customer1 + k - 1),temp1(customer1 + k));
                        temp1_start_times(customer1 + k) = temp1_start_times(customer1 + k - 1) + serv1 + heur1;
						
                        %If we arrive before the opening of a time window don't start service until the time window opens.
                        if temp1_start_times(customer1 + k) < tw(1,temp1(customer1 + k))
                            temp1_start_times(customer1 + k) = tw(1,temp1(customer1 + k));
                        end
                    end

                    %Now add the end of the first route back.
                    k = 1;
                    while solution(route1,customer1 + k) > 0
                        
                        %Add the remainder of the route1 customers.
                        temp1(customer1 + insertion_size + k)		= solution(route1,customer1 + k);
                        
                        %Add the remainder of the route1 customers' loads.
                        temp1_loads(customer1 + insertion_size + k) = loads(route1,customer1 + k);

                        %Add the remainder of the route1 customers' value.
                        temp1_rv(customer1 + insertion_size + k)	= all_value(route1,customer1 + k);

                        %Add the remainder of the route1 customers' start times.
						serv2 = service(temp1(customer1 + insertion_size + k - 1));
						heur2 = heuristic(temp1(customer1 + insertion_size + k - 1),temp1(customer1 + insertion_size + k)); 
                        temp1_start_times(customer1 + insertion_size + k) = temp1_start_times(customer1 + insertion_size + k - 1) + serv2 + heur2;
						
                        %If we arrive before the opening of a time window, don't start service until the time window opens.
                        if temp1_start_times(customer1 + insertion_size + k) < tw(1,temp1(customer1 + insertion_size + k))
                            temp1_start_times(customer1 + insertion_size + k) = tw(1,temp1(customer1 + insertion_size + k));
                        end
                        
                        k = k+1;
                        %If the increment puts us beyond the end of the solution array, end the loop.
                        if (customer1 + k) > size(solution,2)
                            break;
                        end
                    end

%all_value
%--------------------------------------------------------------------------
                    %First new route is now built, so start check load and time windows for feasibility.

                    %Check load of new route. If it exceeds vehicle capacity, then iterate to the next customer22 of route2.
                    if sum(temp1_loads) > Q
                        customer22 = customer22 - 1;
                        continue;
                    end

                    %Check TWs of new route. If any deadlines violated, move to new arc in route. Start at k=2 b/c don't need to check depot times.
                    breaker = 0;
                    k = 2;
                    while temp1(k) > 1
                        if temp1_start_times(k) > tw(2,temp1(k))
                            breaker = 1;
                            break;
                        end
                        k=k+1;
                    end

                    %If breaker == 1, then this route busts a time window. So iterate to the next customer22.
                    if breaker == 1
                        customer22 = customer22 - 1;
                        continue;
                    end
                    
                    %If no value improvement, then don't implement. Move on to next arc.
                    if sum(temp1_rv) <= sum(route_value(route1,:))
                        customer22 = customer22 - 1;
                        continue;
                    end
    %--------------------------------------------------------------------------
                    %Build the second new route.
                    temp2				= solution(route2,:);
                    temp2_loads			= loads(route2,:);
                    temp2_rv			= all_value(route2,:);
                    temp2_start_times	= start_times(route2,:);

                    %Prior to customer21, the route is unchanged.  Zero out everything after that.  Reason for zeroing out is in case the route ends 
					%up shorter than the original, we don't have any vestigial customers from the original route.
                    for k = (customer21+1):(size(temp2,2))
                        temp2(1,k)				= 0;
                        temp2_loads(1,k)		= 0;
                        temp2_rv(1,k)			= 0;
                        temp2_start_times(1,k)	= 0;
                    end

                    %Now add shift everything from customer22+1 to position starting at customer21+1 to account for sequence moved to route1.
                    k = 1;
                    while solution(route2,customer22 + k) > 0
                        
                        %Add the remainder of the route2 customers.
                        temp2(customer21 + k)		= solution(route2,customer22 + k);
                        
                        %Add the remainder of the route1 customers' loads.
                        temp2_loads(customer21 + k) = loads(route2,customer22 + k);

                        %Add the remainder of the route1 customers' loads.
                        temp2_rv(customer21 + k)	= all_value(route2,customer22 + k);

                        %Add the remainder of the route1 customers' start times.
						serv3 = service(temp2(customer21 + k - 1));
						heur3 = heuristic(temp2(customer21 + k - 1),temp2(customer21 + k));
                        temp2_start_times(customer21 + k) = temp2_start_times(customer21 + k - 1) + serv3 + heur3;
						
                        %If we arrive before the opening of a time window, don't start service until the time window opens.
                        if temp2_start_times(customer21 + k) < tw(1,temp2(customer21 + k))
                            temp2_start_times(customer21 + k) = tw(1,temp2(customer21 + k));
                        end
                        
                        k = k+1;
                        %If the increment puts us beyond the end of solution array, end the loop.
                        if (customer22 + k) > size(solution,2)
                            break;
                        end
                    end

                    %Don't need to check load of route2 because all we've done is remove loads.

                    %Don't need to check the time windows of the route2 because we've shifted everything to earlier times (or remained the same).

    %--------------------------------------------------------------------------
                    %New solution is feasible wrt to both vehicle loads time windows, so check for improvement.
                    
                    %First, check route distances against distance limit. Compute the dist of the first route.
                    k = 1;
                    temp1_dist = 0;
                    while temp1(k+1) > 0
                        temp1_dist = temp1_dist + heuristic2(temp1(k),temp1(k+1));
                        k = k+1;
                        if k+1 > size(temp1,2)
                            break;
                        end
                    end

                    %Compute the dist of the second route.
                    k = 1;
                    temp2_dist = 0;
                    while temp2(k+1) > 0
                        temp2_dist = temp2_dist + heuristic2(temp2(k),temp2(k+1));
                        k = k+1;
                        if k+1 > size(temp2,2)
                            break;
                        end
                    end

                    %If distance limit is busted, then don't implement. Move to next arc.
                    if temp1_dist > dist_limit || temp2_dist > dist_limit
                        customer22 = customer22 - 1;
                        continue;
                    end                    

                    %Compute the cost of the first route.
                    k = 1;
                    temp1_cost = 0;
                    while temp1(k+1) > 0
                        temp1_cost = temp1_cost + heuristic(temp1(k),temp1(k+1));
                        k = k+1;
                        if k+1 > size(temp1,2)
                            break;
                        end
                    end

                    %Compute the cost of the second route.
                    k = 1;
                    temp2_cost = 0;
                    while temp2(k+1) > 0
                        temp2_cost = temp2_cost + heuristic(temp2(k),temp2(k+1));
                        k = k+1;
                        if k+1 > size(temp2,2)
                            break;
                        end
                    end

                    %Update the costs associated with the solution.
                    new_cost = cost - route_cost(route1) - route_cost(route2) + temp1_cost + temp2_cost;

    %--------------------------------------------------------------------------
                    %At this point, we've found a cheaper solution that feasible wrt to both vehicle loads and time windows. check the solution for 
					%improvement.implement the solution and return to the calling function.

                    %This is a binary indicator. Set to 1 because we have a better solution.
                    change = 1;

                    %Update the costs associated with the solution.
                    cost				= new_cost;
                    route_cost(route1)	= temp1_cost;
                    route_cost(route2)	= temp2_cost;
                    
                    %Update route distances.
                    route_dist(route1)	= temp1_dist;
                    route_dist(route2)	= temp2_dist;

                    %Update route values.
                    for k = 1:(size(temp1_rv,2))
                        route_value(route1,k)	= temp1_rv(k);
                        all_value(route1,k)		= temp1_rv(k);
                    end

                    %Update route2 loads.
                    for k = 1:(size(temp2_rv,2))
                        all_value(route2,k) = temp2_rv(k);
                    end
                    
                    %Update route1 in the solution.
                    ssize = size(solution,2);
                    for k = 1:(size(temp1,2))
                        solution(route1,k) = temp1(k);
                    end
                    if size(temp1,2) < ssize
                        for k = (size(temp1,2)+1):ssize
                            solution(route1,k) = 0;
                        end
                    end

                    %Update route2 in the solution.
                    for k = 1:(size(temp2,2))
                        solution(route2,k) = temp2(k);
                    end
                    if size(temp2,2) < ssize
                        for k = (size(temp2,2)+1):ssize
                            solution(route2,k) = 0;
                        end
                    end

                    %Update route1 start_times.
                    for k = 1:(size(temp1_start_times,2))
                        start_times(route1,k) = temp1_start_times(k);
                    end
                    if size(temp1_start_times,2) < ssize
                        for k = (size(temp1_start_times,2)+1):ssize
                            solution(route1,k) = 0;
                        end
                    end

                    %Update route2 start_times.
                    for k = 1:(size(temp2_start_times,2))
                        start_times(route2,k) = temp2_start_times(k);
                    end

                    %Update route1 loads.
                    for k = 1:(size(temp1_loads,2))
                        loads(route1,k) = temp1_loads(k);
                    end

                    %Update route2 loads.
                    for k = 1:(size(temp2_loads,2))
                        loads(route2,k) = temp2_loads(k);
                    end

                    %If either of the routes is shorter than the one it's replacing, zero out any vestigial customers/loads/start times.
                    if size(temp1,2) < ssize
                        for k = (size(temp1,2)+1):ssize
                            solution(route1,k)		= 0;
                            loads(route1,k)			= 0;
                            start_times(route1,k)	= 0;
                            all_value(route1,k)		= 0;
                            route_value(route1,k)	= 0;
                        end
                    end
                    
                    if size(temp2,2) < ssize
                        for k = (size(temp2,2)+1):ssize
                            solution(route2,k)		= 0;
                            loads(route2,k)			= 0;
                            start_times(route2,k)	= 0;
                            all_value(route2,k)		= 0;
                        end
                    end
                    
                    %Check the second route and if it's empty (no customers), then delete it.
                    if solution(route2,2) == 1
                        solution(route2,:)		= [];
                        loads(route2,:)			= [];
                        start_times(route2,:)	= [];
                        route_cost(route2)		= [];
                        route_dist(route2)		= [];
                        all_value(route2,:)		= [];
                    end

                    %Run the fix_routes subroutine to make sure we aren't visiting the same customer twice within a single route.
                    [solution, cost, route_cost, start_times, loads, route_dist, route_value, all_value] ...
								= fix_routes_value(	solution, cost, route_cost, heuristic, heuristic2, tw, service, start_times, loads, ...
													route_dist, num_vehicles, route_value, all_value);
                    value = sum(sum(route_value));

                    %Delete any zero columns from the solution, plus the corresponding columns from the loads and start times arrays.
                    index = size(solution,2);
                    for k = 0:(index-1)
                        y = sum(solution,1);
                        if y(index - k) == 0
                            solution(:,(index - k))		= [];
                            loads(:,(index - k))		= [];
                            start_times(:,(index - k))	= [];
                        else
                            break;
                        end
                    end
                    if size(route_value,2) > size(solution,2)
                        index = size(solution,2) - size(route_value,2);
                        for k = 1:index
                            route_value(:,size(route_value,2)) = [];
                        end
                    end                    
                    %We've found an improved solution. Return to the calling function.
                    return;
                end
                
                %If we reach this line, then we've either exhausted customer22 possibilities or we've ruled out any feasible solutions using the 
				%current customer1/customer21 pairing. In either case, move on to the next customer21.
                customer21 = customer21 + 1;
            end
        end
        
        %Go to the next customer in the first route.
        customer1 = customer1 + 1;
        if (customer1+1) > size(solution,2)
            break
        end
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%This function is an intra-route operator, greedily moving customers within a route if a better feasible solution is found.
function [solution, cost, route_cost, start_times, loads,				change] ...
				= relocate(					solution,	cost, route_cost, heuristic,			heuristic3,	tw, service, ~, start_times, loads)

change	= 0;
r		= size(solution,1); % total number of routes

%Run for each route
for route = 1:r
    
    %Try to move each node. Start at 2 and end NLT second-to-last spot because we don't want to move the depot.
    node = 2;
    while node < size(solution,2)
        
        %If we've reached the depot, move on to the next route.
        if solution(route,node) == 1
            break;
        end
        
        %Otherwise, move it to a new spot and calculate costs and time windows.
        new_spot = 2;
        while new_spot < size(solution,2)
            
            %Don't try to move it to the same spot it currently occupies.
            if node == new_spot
                new_spot = new_spot + 1;
                continue;
            end
            
            %If we've reached the depot, move on to next prospective node.
            if solution(route,new_spot) == 1
                break;
            end
            
            %Otherwise, do move it.
            temp		= solution(route,:);
            temp_loads	= loads(route,:);
            temp_st		= zeros(1,size(start_times,2));
            temp_rc		= 0;
            
            %If we're moving the customer to an earlier spot in the route.
            if new_spot < node
                
                %Calculate new solution.
                temp((new_spot+1):node)			= solution(route,new_spot:(node-1));
                temp(new_spot)					= solution(route,node);
                
                %Calculate new loads.
                temp_loads((new_spot+1):node)	= loads(route,new_spot:(node-1));
                temp_loads(new_spot)			= loads(route,node);
                
                %Calculate new costs.
                for i = 2:(size(temp,2))
                    
                    %if we're at the depot, end the loop.
                    if temp(i) == 0
                        break;
                    end
                    
                    temp_rc			= temp_rc + heuristic3(temp(i-1),temp(i));							%Otherwise, add arc to cost.
                    temp_st(i)		= temp_st(i-1) + service(temp(i)) + heuristic(temp(i-1),temp(i));	%And calculate start time.
                    
                    %Check time windows.
                    if temp_st(i) < tw(1,temp(i))
                        temp_st(i)	= tw(1,temp(i));
                    end
                end
                
            %If we're moving the customer to a later spot in the route.
            else
                
                %Calculate new solution.
                temp((node):new_spot-1)			= solution(route,(node+1):new_spot);
                temp(new_spot)					= solution(route,node);

                %Calculate new solution.
                temp_loads((node):new_spot-1)	= loads(route,(node+1):new_spot);
                temp_loads(new_spot)			= loads(route,node);
                
                %Calculate new costs.
                for i = 2:(size(temp,2))
                    
                    %if we're past the depot, end the loop.
                    if temp(i) == 0
                        break;
                    end
                    
                    temp_rc			= temp_rc + heuristic3(temp(i-1),temp(i));							%Otherwise, add arc to cost.
                    temp_st(i)		= temp_st(i-1) + service(temp(i)) + heuristic(temp(i-1),temp(i));	%And calculate start time.
                    
                    %Check time windows.
                    if temp_st(i) < tw(1,temp(i))
                        temp_st(i)	= tw(1,temp(i));
                    end
                end
            end
            
            %If no cost improvement, then don't move.
            if temp_rc >= route_cost(route)
                new_spot = new_spot + 1;
                continue;
            end
            
            %If time windows infeasible, then don't move. Don't check depot times.
            breaker = 0;
            for i = 2:(size(temp_st,2)-1)
                
                %If we've reached the depot, stop.
                if temp(i) == 1
                    break;
                end
                
                %If time window busted, then we can't move to this spot. So kill this loop.
                if temp_st(i) > tw(2,temp(i))
                    breaker = 1;
                    break;
                end
            end
            
            %And iterate to the next spot.
            if breaker == 1
                new_spot = new_spot + 1;
                continue;
            end
            
%-------------------------------------------------------------------------%
            %If we reach this point, then we have a feasible solution with an improved cost, so implement the solution.
            
            change					= 1;				%Set binary indicator to 1.
            solution(route,:)		= temp;				%Change solution.
            loads(route,:)			= temp_loads;		%Change loads.
            start_times(route,:)	= temp_st;			%Change start times.
            route_cost(route)		= temp_rc;			%Change route costs.
            cost					= sum(route_cost);	%Change cost.
            node					= 2;				%Since we've now relocated a node, start over with this route.
            break;
        end
        
        node = node + 1;
                
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%This function is an intra-route operator, greedily moving customers within a route if a better feasible solution is found.
function [solution, cost, route_cost, start_times, loads,						value, route_value,				change] ...
				= relocate_value(			solution,	cost, route_cost, heuristic,						tw, service, ~, start_times, loads, ...
																				value, route_value)

change = 0;

%r = total number of routes
%r = size(solution,1);
r = size(route_value,1);%min(num_vehicles,size(solution,1));

%Run for each route
for route = 1:r
    
    %Try to move each node. Start at 2 and end NLT second-to-last spot because we don't want to move the depot.
    node = 2;
    while node < size(solution,2)
        
        %If we've reached the depot, move on to the next route.
        if solution(route,node) == 1
            break;
        end
        
        %Otherwise, move it to a new spot and calculate costs and time windows.
        new_spot = 2;
        while new_spot < size(solution,2)
            
            %Don't try to move it to the same spot it currently occupies.
            if node == new_spot
                new_spot = new_spot + 1;
                continue;
            end
            
            %If we've reached the depot, move on to next prospective node.
            if solution(route,new_spot) == 1
                break;
            end
            
            %Otherwise, do move it.
            temp		= solution(route,:);
            temp_loads	= loads(route,:);
            temp_st		= zeros(1,size(start_times,2));
            temp_rc		= 0;
            temp_rv		= route_value(route,:);
            
            %If we're moving the customer to an earlier spot in the route.
            if new_spot < node
                
                %Calculate new solution.
                temp((new_spot+1):node)			= solution(route,new_spot:(node-1));
                temp(new_spot)					= solution(route,node);
                
                %Calculate new loads.
                temp_loads((new_spot+1):node)	= loads(route,new_spot:(node-1));
                temp_loads(new_spot)			= loads(route,node);
                
                %Calculate new values.
                temp_rv((new_spot+1):node)		= route_value(route,new_spot:(node-1));
                temp_rv(new_spot)				= route_value(route,node);

                %Calculate new costs.
                for i = 2:(size(temp,2))
                    
                    %if we're at the depot, end the loop.
                    if temp(i) == 0
                        break;
                    end
                    
                    temp_rc		= temp_rc + heuristic(temp(i-1),temp(i));							% Otherwise, add arc to cost.
                    temp_st(i)	= temp_st(i-1) + service(temp(i)) + heuristic(temp(i-1),temp(i));	% And calculate start time.
                    
                    %Check time windows.
                    if temp_st(i) < tw(1,temp(i))
                        temp_st(i) = tw(1,temp(i));
                    end
                end
                
            %If we're moving the customer to a later spot in the route.
            else
                
                %Calculate new solution.
                temp((node):new_spot-1)			= solution(route,(node+1):new_spot);
                temp(new_spot)					= solution(route,node);

                %Calculate new loads.
                temp_loads((node):new_spot-1)	= loads(route,(node+1):new_spot);
                temp_loads(new_spot)			= loads(route,node);
                
                %Calculate new loads.
                temp_rv((node):new_spot-1)		= route_value(route,(node+1):new_spot);
                temp_rv(new_spot)				= route_value(route,node);

                %Calculate new costs.
                for i = 2:(size(temp,2))
                    
                    %if we're past the depot, end the loop.
                    if temp(i) == 0
                        break;
                    end
                    
                    %Otherwise, add arc to cost.
                    temp_rc = temp_rc + heuristic(temp(i-1),temp(i));
                    
                    %And calculate start time.
                    temp_st(i) = temp_st(i-1) + service(temp(i)) + heuristic(temp(i-1),temp(i));
                    
                    %Check time windows.
                    if temp_st(i) < tw(1,temp(i))
                        temp_st(i) = tw(1,temp(i));
                    end
                end
            end
            
            %If no cost improvement, then don't move.
            if temp_rc >= route_cost(route)
                new_spot = new_spot + 1;
                continue;
            end
            
            %If time windows infeasible, then don't move. Don't check depot times.
            breaker = 0;
            for i = 2:(size(temp_st,2)-1)
                
                %If we've reached the depot, stop.
                if temp(i) == 1
                    break;
                end
                
                %If time window busted, then we can't move to this spot. So kill this loop.
                if temp_st(i) > tw(2,temp(i))
                    breaker = 1;
                    break;
                end
            end
            
            %And iterate to the next spot.
            if breaker == 1
                new_spot = new_spot + 1;
                continue;
            end
            
%-------------------------------------------------------------------------%
            %If we reach this point, then we have a feasible solution with an improved cost, so implement the solution.
            change				= 1;			%Set binary indicator to 1.
            solution(route,:)	= temp;			%Change solution.
            loads(route,:)		= temp_loads;	%Change loads.

            %Change start times.
            for j = 1:size(start_times,2)
                if j > size(temp_st,2)
                    start_times(route,j) = 0;
                else
                    start_times(route,j) = temp_st(j);
                end
            end
            
            %Change route costs.
            route_cost(route) = temp_rc;
            
            %Change cost.
            cost = sum(route_cost);
            
            %Move the two values to correspond with the new order of delivery.
            route_value(route,:) = temp_rv;
            
            %Since we've now relocated a node, start over with this route.
            node = 2;
            break;
        end
        
        node = node + 1;
                
    end
end

end




% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%This function will calculate the value of each delivery given a solution.
function val	= route_values(routes, loads, current_inventory, theta)

%Initialize
temp_inventory = current_inventory;
val = routes.*0;

%Using current inventory and total demand, calculate the increase in value by delivering a load to each customer.
for i = 1:size(routes,1)
    for j = 2:size(routes,2)
        %If we've reached the end of a route, break the j loop and go to the next iteration of the i loop.
        if routes(i,j) == 1
            break;

        end

        %Calculate the value of a delivery.
		vf1		 = value_function( loads( i,j) + temp_inventory( routes( i,j))	, theta, 1);
		vf2		 = value_function( temp_inventory( routes( i,j))				, theta, 1);
        val(i,j) = vf1 - vf2;
        
        %Update the inventory to reflect that delivery so that future value calculations are accurate.
        temp_inventory(routes(i,j)) = temp_inventory(routes(i,j)) + loads(i,j);
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%Given an initial solution, this program will run 2-opt* to try to find an improved solution.
function [solution, cost, route_cost, start_times, loads, route_dist,											change] ...
				= two_opt_star(				solution,	cost, route_cost, heuristic, heuristic2, heuristic3, tw, service, Q, start_times, loads, ...
										dist_limit, route_dist)

change	  = 0; %Initialization

if size(solution,2) > size(loads,2)
    temp  = zeros(size(loads,1),(size(solution,2)-size(loads,2)));
    loads = [loads,temp];
end


%Randomly re-order an index of routes.  This will tell us what order we're going to look at routes.
routes = transpose(randperm(size(solution,1)));

%For each route, see if 2-opt* will yield an improved solution. This loop controls the first route.
for i = 1:(size(solution,1)-1)
    
    %Assign the first route, then delete it from the list of routes to check.
    route1		= routes(i);
    num_routes	= size(routes,1);
    
    %Check each arc of the route in question. This loop controls the arc from the first route.
    customer1	= 1;
    while solution(route1,customer1+1) > 0
        
        %Select a route with which to perform the switch. Note that we will only run 2-opt* on the routes occuring after the current route in
        %this array, avoiding unnecessary double work. This loop controls the second route.
        for j = (i+1):num_routes

            %Assign the second route.
            route2 = routes(j);
            
            %When checking  first arc in first route, don't try to swap with first arc from second route. This swap will yield an identical solution.
            if customer1 == 1
                customer2 = 2;
            else
                customer2 = 1;
            end
            
            %Select arc from second route. This loop contols the arc from the second route.
            while solution(route2,customer2+1) > 0

                %Check the arcs to make sure we're making a valid exchange. First check is to make sure that we're not just swapping depot returns.
                if (solution(route2,customer2+1) == 1) && (solution(route1,customer1+1) == 1)
                    %If we are, that means we must be at the end of both routes.  So end the current route2 iteration and go to next second route.
                    break;
                    
                    %Next check ensure depot isn't self-connected (i.e., first node of the second route to last node of the first route).
                elseif (solution(route2,customer2) == 1) && (solution(route1,customer1+1) == 1)
                    %If this is true, then go to the next customer in the second route.
                    customer2 = customer2 + 1;
                    continue;
                    
                    %Final, re-check if depot is self-connected (i.e., first node of the first route to last node of the second route).
                elseif (solution(route2,customer2+1) == 1) && (solution(route1,customer1) == 1)
                    %If this is true, then we've reached the end of this route2.  So go to the next second route.
                    break;
                end
                
%--------------------------------------------------------------------------
                %Now we have two valid arcs from two valid routes. Build the first new route.
                temp1 = solution(route1,:);
                temp1_loads = loads(route1,:);
                temp1_start_times = start_times(route1,:);
                
                %Prior to customer1, route1 is unchanged.  So zero out everything after that.  Reason for zeroing out is in case the route ends up 
				%shorter than the original, we don't have any vestigial customers from the original route.
                for k = (customer1+1):(size(temp1,2))
                    temp1(1,k) = 0;
                    temp1_loads(1,k) = 0;
                    temp1_start_times(1,k) = 0;
                end

                %Now, add in the piece of the second route.
                k = 1;
				
                %Stop when we reach the depot return for route2.
                while solution(route2,customer2 + k) > 0
                    if customer2 + k > size(loads,2)
                        break;
                    end
                    
                    temp1(customer1 + k) = solution(route2,customer2 + k); %Add the route2 customer to route1.
                    temp1_loads(customer1 + k) = loads(route2,customer2 + k); %Add the route2 customer's to route1 loads.
                    
                    %Add the route2 customer's start time to route1 start times.
					tmp1								= temp1_start_times(customer1 + k - 1);
					serv								= service(temp1(customer1 + k - 1));
					heur								= heuristic(temp1(customer1 + k - 1),temp1(customer1 + k));
					temp1_start_times(customer1 + k)	=  tmp1 + serv + heur;
                    
					%If we arrive before the opening of a time window, don't start service until the time window opens.
                    if temp1_start_times(customer1 + k) < tw(1,temp1(customer1 + k))
                        temp1_start_times(customer1 + k) = tw(1,temp1(customer1 + k));
                    end
                    
                    k = k+1;
                    %If the increment puts us beyond the end of the solution array, end the loop.
                    if (customer2 + k) > size(solution,2)
                        break;
                    end
                end

                %First new route is now built, so start check load and time windows for feasibility.

                %Check load of new route. If it exceeds vehicle capacity, then iterate to the next customer of route2.
                if sum(temp1_loads) > Q
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
                %Check TWs of new route. If any deadlines are violated, move to new arc in route. Start at k=2 b/c don't need to check depot times.
                breaker = 0;
                k = 2;

                while temp1(k) > 1
                    if temp1_start_times(k) > tw(2,temp1(k))
                        breaker = 1;
                        break;
                    end
                    k=k+1;
                end
                
                %If breaker == 1, then this route busts a time window.
                if breaker == 1
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
%--------------------------------------------------------------------------
                %Build the second new route.
                temp2 = solution(route2,:);
                temp2_loads = loads(route2,:);
                temp2_start_times = start_times(route2,:);
                
                %Prior to customer2, the route is unchanged.  Zero out everything after that.  Reason for zeroing out is in case the route ends up 
				%shorter than the original, we don't have any vestigial customers from the original route.
                for k = (customer2+1):(size(temp2,2))
                    temp2(1,k)				= 0;
                    temp2_loads(1,k)		= 0;
                    temp2_start_times(1,k)	= 0;
                end
                
                %Now, add in the piece of the second route.
                k = 1;
                %Stop when we reach the depot return for route2.
                while solution(route1,customer1 + k) > 0
                    
                    temp2(customer2 + k)				= solution(route1,customer1 + k);	%Add the route1 customer to route2.
                    temp2_loads(customer2 + k)			= loads(route1,customer1 + k);		%Add the route1 customer's load to route2 loads.

                    %Add the route2 customer's start time to route1 start times.
					tmp2								= temp2_start_times(customer2 + k - 1);
					serv								= service(temp2(customer2 + k - 1));
					heur								= heuristic(temp2(customer2 + k - 1),temp2(customer2 + k));
                    temp2_start_times(customer2 + k)	= tmp2 + serv + heur;
					
                    %If we arrive before the opening of a time window, don't start service until the time window opens.
                    if temp2_start_times(customer2 + k) < tw(1,temp2(customer2 + k))
                        temp2_start_times(customer2 + k) = tw(1,temp2(customer2 + k));
                    end
                    
                    k = k+1;
                    %If the increment puts us beyond the end of the solution array, end the loop.
                    if (customer1 + k) > size(solution,2)
                        break;
                    end
                end
                
                %Check load of new route. If it exceeds vehicle capacity, then move on to a new route2. We don't need to iterate to the next customer
				%in the current route2 because that would just the same piece of route1 along with another route2 customer.
                if sum(temp2_loads) > Q
                    break;
                end
                
                %Check the time windows of the new route. If any deadlines are violated, move on.
                k = 2;
                while temp2(k) > 1
                    if temp2_start_times(k) > tw(2,temp2(k))
                        breaker = 1;
                        break;
                    end
                    k = k+1;
                end

                %Moving to the next arc in route2 will just add another customer to route2. Since we satisfy the triangle inequality, then adding 
				%another customer to route2 can't possibly satisfy the time windows. So break the arc selection loop and move on to the next route2.
                if breaker == 1
                    break;
                end
                
%--------------------------------------------------------------------------
                %New solution is feasible wrt to both vehicle loads and time windows, so check for improvement.
                %Compute the dist of the first route.
                k			= 1;
                temp1_dist	= 0;
                while temp1(k+1) > 0
                    temp1_dist	= temp1_dist + heuristic2(temp1(k),temp1(k+1));
                    k			= k+1;
                    if k+1 > size(temp1,2)
                        break;
                    end
                end
                
                %Compute the dist of the second route.
                k			= 1;
                temp2_dist	= 0;
                while temp2(k+1) > 0
                    temp2_dist	= temp2_dist + heuristic2(temp2(k),temp2(k+1));
                    k			= k+1;
                    if k+1 > size(temp2,2)
                        break;
                    end
                end
                
                %If distance limit is busted, then don't implement. Move on to next arc.
                if temp1_dist > dist_limit || temp2_dist > dist_limit
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
                %Compute the cost of the first route.
                k			= 1;
                temp1_cost	= 0;
                while temp1(k+1) > 0
                    temp1_cost	= temp1_cost + heuristic3(temp1(k),temp1(k+1));
                    k			= k+1;
                    if k+1 > size(temp1,2)
                        break;
                    end
                end
                
                %Compute the cost of the second route.
                k			= 1;
                temp2_cost	= 0;
                while temp2(k+1) > 0
                    temp2_cost = temp2_cost + heuristic3(temp2(k),temp2(k+1));
                    k = k+1;
                    if k+1 > size(temp2,2)
                        break;
                    end
                end

                %Update the costs associated with the solution.
                new_cost	= cost - route_cost(route1) - route_cost(route2) + temp1_cost + temp2_cost;
                
                %If no cost improvement, then don't implement. Move on to next arc.
                if new_cost >= cost
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end

%--------------------------------------------------------------------------
                %At this point, we've found a cheaper solution that is feasible wrt to both vehicle loads and time windows. So check the solution for
				%improvement. Implement the solution and return to the calling function.
                
                %This is a binary indicator. Set to 1 because we have found a better solution.
                change = 1;

                %Update the costs associated with the solution.
                cost				= new_cost;
                route_cost(route1)	= temp1_cost;
                route_cost(route2)	= temp2_cost;
                
                %Update route distances.
                route_dist(route1)	= temp1_dist;
                route_dist(route2)	= temp2_dist;
                                
                %Update route1 in the solution.
                for k = 1:(size(temp1,2))
                    solution(route1,k) = temp1(k);
                end
                
                %Update route2 in the solution.
                for k = 1:(size(temp2,2))
                    solution(route2,k) = temp2(k);
                end
                
                %Update route1 start_times.
                for k = 1:(size(temp1_start_times,2))
                    start_times(route1,k) = temp1_start_times(k);
                end
                
                %Update route2 start_times.
                for k = 1:(size(temp2_start_times,2))
                    start_times(route2,k) = temp2_start_times(k);
                end
                
                %Update route1 loads.
                for k = 1:(size(temp1_loads,2))
                    loads(route1,k) = temp1_loads(k);
                end
                
                %Update route2 loads.
                for k = 1:(size(temp2_loads,2))
                    loads(route2,k) = temp2_loads(k);
                end
                
                %Run the fix_routes subroutine to make sure we aren't visiting the same customer twice within a single route.
                [solution, cost, route_cost, start_times, loads, route_dist] ...
					= fix_routes(solution, cost, route_cost, heuristic, heuristic2, heuristic3, tw, service, start_times, loads, route_dist);

                %Delete any zero columns from the solution, plus the
                %corresponding columns from the loads and start times
                %arrays.
                index = size(solution,2);
                for k = 0:(index-1)
                    y = sum(solution,1);
                    if y(index - k) == 0
                        solution(:,(index - k)) = [];
                        loads(:,(index - k)) = [];
                        start_times(:,(index - k)) = [];
                    else
                        break;
                    end
                end
                
                %We've found an improved solution. Return to the calling
                %function.
                return;
            end
        end
        
        %Go to the next customer in the first route.
        customer1 = customer1 + 1;
        
        if (customer1 + 1) > (size(solution,2))
            break;
        end
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%Given an initial solution, this program will run 2-opt* to try to find an improved solution.
function [solution, cost, route_cost, start_times, loads, route_dist,			value, route_value, all_value,	change] ...
				= two_opt_star_value(		solution,	cost, route_cost, heuristic, heuristic2,			tw, service, Q, start_times, loads, ...
										dist_limit, route_dist, num_vehicles,	value, route_value, all_value, current_inventory, theta)

						
%Initialization
change = 0;


if size(route_value,1) < num_vehicles
    temp		= zeros(num_vehicles - size(route_value,1),size(route_value,2));
    route_value = [route_value;temp];
end

%Randomly re-order an index of routes.  This will tell us what order we're going to look at routes.
good_routes = transpose(randperm(num_vehicles,1));
bad_routes	= transpose(randperm(size(solution,1)-num_vehicles,1));

%For each route, see if 2-opt* will yield an improved solution. 
%This loop controls the first route.
% while numel(routes) > 1
for i = 1:size(good_routes)
    
    %Assign the first route, then delete it from the list of routes to check.
    route1 = good_routes(i);
%     num_routes = size(routes,1);
    
    %Check each arc of the route in question. This loop controls the arc from the first route.
    customer1 = 1;
    while solution(route1,customer1+1) > 0
        
        %Select a route with which to perform the switch. Note that we will only run 2-opt* on the routes occuring after the current route in
        %this array, avoiding unnecessary double work. This loop controls the second route.
        for j = 1:size(bad_routes)

            %Assign the second route.
            route2 = bad_routes(j)+num_vehicles;
            
            %if we're checking the first arc in the first route, then don't try to swap with the first arc from the second route. 
			%This swap will yield an identical solution.
            if customer1 == 1
                customer2 = 2;
            else
                customer2 = 1;
            end
            
            %Select arc from second route. This loop contols the arc from the second route.
            while solution(route2,customer2+1) > 0

                %Check the arcs to make sure we're making a valid exchange. First check is to make sure that we're not just swapping depot returns.
                if (solution(route2,customer2+1) == 1) && (solution(route1,customer1+1) == 1)
                    %If we are, we must be at the end of both routes.  So end the current route2 iteration and go to the next second route.
                    break;
                end
                
%--------------------------------------------------------------------------
                %Now we have two valid arcs from two valid routes.
                %Build the first new route.
                temp1 = solution(route1,:);
                temp1_loads = loads(route1,:);
                temp1_start_times = start_times(route1,:);
                temp1_rv = all_value(route1,:);
                
                %Prior to customer1, route1 is unchanged.  So zero out
                %everything after that.  Reason for zeroing out is in case
                %the route ends up shorter than the original, we don't have
                %any vestigial customers from the original route.
                for k = (customer1+1):(size(temp1,2))
                    temp1(1,k) = 0;
                    temp1_loads(1,k) = 0;
                    temp1_start_times(1,k) = 0;
                    temp1_rv(1,k) = 0;
                end

                %Now, add in the piece of the second route.
                k = 1;
                %Stop when we reach the depot return for route2.
                while solution(route2,customer2 + k) > 0
                    
                    temp1(customer1 + k)		= solution( route2,customer2 + k); %Add the route2 customer to route1.
                    temp1_loads(customer1 + k)	= loads(	route2,customer2 + k); %Add the route2 customer's to route1 loads.
                    temp1_rv(customer1 + k)		= all_value(route2,customer2 + k); %Add the route2 customer's to route1 value.
                    
                    %Add the route2 customer's start time to route1 start times.
					serv								= service(temp1(customer1 + k - 1));
					heur								= heuristic(temp1(customer1 + k - 1),temp1(customer1 + k));
                    temp1_start_times(customer1 + k)	= temp1_start_times(customer1 + k - 1) + serv + heur;
					
                    %If we arrive before the opening of a time window, don't start service until the time window opens.
                    if temp1_start_times(customer1 + k) < tw(1,temp1(customer1 + k))
                        temp1_start_times(customer1 + k) = tw(1,temp1(customer1 + k));
                    end
                    
                    k = k+1;
                    %If the increment puts us beyond the end of the solution array, end the loop.
                    if (customer2 + k) > size(solution,2)
                        break;
                    end
                end

                %First new route is now built, so start check load and time windows for feasibility.
                %Check load of new route. If it exceeds vehicle capacity, then iterate to the next customer of route2.
                if sum(temp1_loads) > Q
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
                %Check TWs of new route. If any deadlines are violated, move to new arc in route. Start at k=2 b/c don't need to check depot times.
                breaker = 0;
                k = 2;
                while temp1(k) > 1
                    if temp1_start_times(k) > tw(2,temp1(k))
                        breaker = 1;
                        break;
                    end
                    k=k+1;
                end
                
                %If breaker == 1, then this route busts a time window.
                if breaker == 1
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
                %If no value improvement, then don't implement. Move on to
                %next arc.
%                 temp1_rv = route_values(temp1, temp1_loads, current_inventory, theta);
                
                if sum(temp1_rv) <= sum(route_value(route1,:))
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
%--------------------------------------------------------------------------
                %Build the second new route.
                temp2 = solution(route2,:);
                temp2_loads = loads(route2,:);
                temp2_start_times = start_times(route2,:);
                temp2_rv = all_value(route2,:);
                
                %Prior to customer2, the route is unchanged.  Zero out everything after that.  Reason for zeroing out is in case the route ends up 
				%shorter than the original, we don't have any vestigial customers from the original route.
                for k = (customer2+1):(size(temp2,2))
                    temp2(1,k) = 0;
                    temp2_loads(1,k) = 0;
                    temp2_start_times(1,k) = 0;
                    temp2_rv(1,k) = 0;
                end
                
                %Now, add in the piece of the second route.
                k = 1;
                %Stop when we reach the depot return for route2.
                while solution(route1,customer1 + k) > 0
                    
                    temp2(		customer2 + k)	= solution(route1,customer1 + k);	%Add the route1 customer to route2.
                    temp2_loads(customer2 + k)	= loads(route1,customer1 + k);		%Add the route1 customer's load to route2 loads.
                    temp2_rv(	customer2 + k)	= all_value(route1,customer1 + k);	%Add the route1 customer's to route2 value.

                    %Add the route2 customer's start time to route1 start times.
					serv								= service(temp2(customer2 + k - 1));
					heur								= heuristic(temp2(customer2 + k - 1),temp2(customer2 + k));
                    temp2_start_times(customer2 + k)	= temp2_start_times(customer2 + k - 1) + serv + heur;
                    
					%If we arrive before the opening of a time window, don't start service until the time window opens.
                    if temp2_start_times(customer2 + k) < tw(1,temp2(customer2 + k))
                        temp2_start_times(customer2 + k) = tw(1,temp2(customer2 + k));
                    end
                    
                    k = k+1;
                    %If the increment puts us beyond the end of the solution array, end the loop.
                    if (customer1 + k) > size(solution,2)
                        break;
                    end
                end
                  
                %Check the time windows of the new route. If any deadlines are violated, move on.
                k = 2;
                while temp2(k) > 1
                    if temp2_start_times(k) > tw(2,temp2(k))
                        breaker = 1;
                        break;
                    end
                    k=k+1;
                end

                %Moving to the next arc in route2 will just add another customer to route2. Since we satisfy the triangle inequality, then adding 
				%another customer to route2 can't possibly satisfy the time windows. So break the arc selection loop and move on to the next route2.
                if breaker == 1
                    break;
                end
                
                
%--------------------------------------------------------------------------
                %New solution is feasible wrt to both vehicle loads and time windows, so check for improvement.
                %Compute the dist of the first route.
                k = 1;
                temp1_dist = 0;
                while temp1(k+1) > 0
                    temp1_dist = temp1_dist + heuristic2(temp1(k),temp1(k+1));
                    k = k+1;
                    if k+1 > size(temp1,2)
                        break;
                    end
                end
                
                %Compute the dist of the second route.
                k = 1;
                temp2_dist = 0;
                while temp2(k+1) > 0
                    temp2_dist = temp2_dist + heuristic2(temp2(k),temp2(k+1));
                    k = k+1;
                    if k+1 > size(temp2,2)
                        break;
                    end
                end
                
                %If distance limit is busted, then don't implement. Move on to next arc.
                if temp1_dist > dist_limit %|| temp2_dist > dist_limit
                    customer2 = customer2 + 1;
                    if (customer2 + 1) > (size(solution,2))
                        break;
                    end
                    continue;
                end
                
                %Compute the cost of the first route.
                k = 1;
                temp1_cost = 0;
                while temp1(k+1) > 0
                    temp1_cost = temp1_cost + heuristic(temp1(k),temp1(k+1));
                    k = k+1;
                    if k+1 > size(temp1,2)
                        break;
                    end
                end
                
                %Compute the cost of the second route.
                k = 1;
                temp2_cost = 0;
                while temp2(k+1) > 0
                    temp2_cost = temp2_cost + heuristic(temp2(k),temp2(k+1));
                    k = k+1;
                    if k+1 > size(temp2,2)
                        break;
                    end
                end

                %Update the costs associated with the solution.
                new_cost = cost - route_cost(route1) - route_cost(route2) + temp1_cost + temp2_cost;
                
%--------------------------------------------------------------------------
                %At this point, we've found a cheaper solution that is feasible wrt to both vehicle loads and time windows. So check the solution for
				%improvement.implement the solution and return to the calling function.
                
                %This is a binary indicator. Set to 1 because we have found a better solution.
                change = 1;

                %Update the costs associated with the solution.
                cost						= new_cost;
                route_cost(route1)			= temp1_cost;
                route_cost(route2)			= temp2_cost;
                
                %Update route distances.
                route_dist(route1)			= temp1_dist;
                route_dist(route2)			= temp2_dist;
                
                %Update route values.
                for k = 1:(size(temp1_rv,2))
                    route_value(route1,k)	= temp1_rv(k);
                    all_value(route1,k)		= temp1_rv(k);
                end                
                for k = 1:(size(temp2_rv,2))
                    all_value(route2,k)		= temp2_rv(k);
                end                
                                
                %Update route1 in the solution.
                for k = 1:(size(temp1,2))
                    solution(route1,k)		= temp1(k);
                end
                
                %Update route2 in the solution.
                for k = 1:(size(temp2,2))
                    solution(route2,k)		= temp2(k);
                end
                
                %Update route1 start_times.
                for k = 1:(size(temp1_start_times,2))
                    start_times(route1,k)	= temp1_start_times(k);
                end
                
                %Update route2 start_times.
                for k = 1:(size(temp2_start_times,2))
                    start_times(route2,k)	= temp2_start_times(k);
                end
                
                %Update route1 loads.
                for k = 1:(size(temp1_loads,2))
                    loads(route1,k)			= temp1_loads(k);
                end
                
                %Update route2 loads.
                for k = 1:(size(temp2_loads,2))
                    loads(route2,k)			= temp2_loads(k);
                end
                
                %Run the fix_routes subroutine to make sure we aren't visiting the same customer twice within a single route.
                [solution, cost, route_cost, start_times, loads, route_dist, route_value, all_value] ...
					= fix_routes_value(solution, cost, route_cost, heuristic, heuristic2, tw, service, ...
										start_times, loads, route_dist, num_vehicles, route_value, all_value);
                value = sum(sum(route_value));
                
                %Delete any rows that depict a depot-to-depot route.
                for inc = size(solution,1):-1:1
                    if solution(inc,2) == 1
                        solution(inc,:)			= [];
                        loads(	 inc,:)			= [];
                        all_value(inc,:)		= [];
                        if inc <= num_vehicles
                            route_value(inc,:)	= [];
                        end
                        start_times(inc,:)		= [];
                        route_dist( inc,:)		= [];
                        route_cost( inc,:)		= [];
                    end
                end

                %Delete any zero columns from the solution, plus the corresponding columns from the loads and start times arrays.
                index	= size(solution,2);
                y		= sum(solution,1);
                for k = 0:(index-1)
                    if y(index - k) == 0
                        solution(:,(index - k))		= [];
                        loads(:,(index - k))		= [];
                        start_times(:,(index - k))	= [];
                        all_value(:,(index - k))	= [];
                    else
                        break;
                    end
                end
                if size(route_value,2) > size(solution,2)
                    index = size(solution,2) - size(route_value,2);
                    for k = 1:index
                        route_value(:,size(route_value,2)) = [];
                    end
                end
                
                %We've found an improved solution. Return to the calling function.
                return;
            end
        end
        
        %Go to the next customer in the first route.
        customer1 = customer1 + 1;
        
        if (customer1 + 1) > (size(solution,2))
            break;
        end
    end
end

end


% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%Trims unneeded 
function v = TrimFeeder(num_vehicles, vect)
	t  = min(num_vehicles, size(vect,1));
	v  = vect(1:t,:);
end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%Used by Or_opt_v2 and two_opt_star
%Ensures a customer isn't visited multiple times within a route. If so, the multiple visits are combined into a single visit at the first location. 
function [solution, cost, route_cost, start_times, loads, route_dist] ...
					= fix_routes(			solution, cost, route_cost, heuristic, heuristic2, heuristic3,	tw, service,	start_times, loads, ...
																	route_dist)


rows	= size(solution,1);

%columns = total size of matrix
columns = size(solution,2);

%Iterate through each route
for route = 1:rows
    
    %For each route, check each of the customers to see if they occur twice. Note that we don't need to check the depot.
    for index1 = 2:columns

        %If the customer directly after customer1 is the depot, then we're done checking this route.
        if solution(route,index1+1) <= 1
            break;
        elseif (index1+1) > size(solution,2)
            break;
        end
        
        %For each customer1, check each customer that occurs after it.
        for index2 = index1+1:columns
        
            %If customer2 is the depot, then we're done checking customer1. Move on to next customer1.
            if solution(route,index2) <= 1
               break;
            end

            %Assign the customers.
            customer1 = solution(route,index1);
            customer2 = solution(route,index2);

            %Check to see if the customers are the same.
            if customer1 == customer2
                
                %If they are the same, then combine them into customer1's location. Update solution, cost, route_cost, start_times, and loads.
                temp				= solution(route,:);
                temp_start_times	= start_times(route,:);
                temp_loads			= loads(route,:);
                for k = index2:(size(temp,2))
                    temp_start_times(1,k) = 0;
                end
                
                %Update solution.
                temp(index2:(size(temp,2)-1))	= temp((index2 + 1):size(temp,2));
                temp(size(temp,2))				= 0;
                
                %Update costs.
				heur31				= heuristic3(solution(route,index2 - 1),	solution(route,index2 + 1));
				heur32				= heuristic3(solution(route,index2 - 1),	solution(route,index2));
				heur33				= heuristic3(solution(route,index2),		solution(route,index2 + 1));
                route_change		= heur31 - heur32 - heur33;
				
                cost				= cost + route_change;
                route_cost(route)	= route_cost(route) + route_change;
                
                %Update distance.
				heur21				= heuristic2(solution(route,index2 - 1),	solution(route,index2 + 1));
				heur22				= heuristic2(solution(route,index2 - 1),	solution(route,index2));
				heur23				= heuristic2(solution(route,index2),		solution(route,index2 + 1));
                route_change2		= heur21 - heur22 - heur23; 
				
                route_dist(route)	= route_dist(route) + route_change2;                
                solution(route,:)	= temp;

                %Update start times. We've updated the solution, so index2 is no longer a duplicate. Nothing prior to index2 has changed wrt start 
				%times, so just update everything from there on.
                i = index2;
                while solution(route,i) > 0
                
                    %Update times.
                    temp_start_times(i) = temp_start_times(i - 1) + service(solution(route,i - 1)) + heuristic(solution(route,i - 1),solution(route,i));
                    
                    %Check time windows. We're shifting everything to an earlier time, so only need to check beginning of time windows.
                    if temp_start_times(i) < tw(1,solution(route,i))
                        temp_start_times(i) = tw(1,solution(route,i));
                    end
                    
                    %Increment
                    i = i+1;
                    
                    %Make sure we don't overrun the matrix.
                    if i > size(solution,2)
                        break;
                    end
                end
                start_times(route,:) = temp_start_times;
                
                %Update loads.
                temp_loads(index1)	= temp_loads(index1) + temp_loads(index2);
                temp_loads(index2)	= [];
                temp_loads			= [temp_loads,0];
                loads(route,:)		= temp_loads;
            end
        end
    end
end


end


% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
%Used by Or_opt_v2_value and two_opt_star_value
%Ensures a customer isn't visited multiple times within a route. If so, the multiple visits are combined into a single visit at the first location.
function [solution, cost, route_cost, start_times, loads, route_dist, route_value, all_value] ...
					= fix_routes_value(		solution, cost, route_cost, heuristic, heuristic2,				tw, service,    start_times, loads, ...
																	route_dist, num_vehicles, route_value, all_value)


[rows, columns] = size(solution); % rows = total number of routes, columns = total size of matrix

%Iterate through each route
for route = 1:rows
	
	%For each route, check each of the customers to see if they occur twice. Note that we don't need to check the depot.
	for index1 = 2:columns
		
		%If the customer directly after customer1 is the depot, then we're done checking this route.
		if solution(route,index1+1) <= 1, break; end
		
		%For each customer1, check each customer that occurs after it.
		for index2 = index1+1:columns
			
			%If customer2 is the depot, then we're done checking the customer1. Move on to next customer1.
			if solution(route,index2) <= 1, break; end
			
			%Assign the customers.
			customer1 = solution(route,index1);
			customer2 = solution(route,index2);
			
			%Check to see if the customers are the same.
			if customer1 == customer2
				
				%If they are the same, then combine them into customer1's location. Update solution, cost, route_cost, start_times, and loads.
				temp				= solution(		route, :);
				temp_start_times	= start_times(	route, :);
				temp_loads			= loads(		route, :);
				temp_rv				= all_value(	route, :);
				for k = index2:(size(temp,2))
					temp_start_times(1,k) = 0;
				end
				
				%Update solution.
				temp(index2:(size(temp,2)-1)) = temp((index2 + 1):size(temp,2));
				temp(size(temp,2)) = 0;
				
				%Update costs.
				heur1				= heuristic( solution( route, index2 - 1),		solution( route, index2 + 1));
				heur2				= heuristic( solution( route, index2 - 1),		solution( route, index2));
				heur3				= heuristic( solution( route, index2),			solution( route, index2 + 1));
				route_change		= heur1 - heur2 - heur3;
				cost				= cost				+ route_change;
				route_cost(route)	= route_cost(route) + route_change;
				
				%Update distance.
				heur21				= heuristic2(solution(route,index2 - 1),	solution(route,index2 + 1));
				heur22				= heuristic2(solution(route,index2 - 1),	solution(route,index2));
				heur23				= heuristic2(solution(route,index2),		solution(route,index2 + 1));
				route_change2		= heur21 - heur22 - heur23;
				route_dist(route)	= route_dist(route) + route_change2;
				solution(route,:)	= temp;
				
				%Update start times. We've updated the solution, so index2 is no longer a duplicate. Nothing prior to index2 has changed wrt start
				%times, so just update everything from there on.
				i = index2;
				while solution(route,i) > 0
					
					%Update times.
					tempst				= temp_start_times( i - 1);
					serv				= service(solution( route, i - 1));
					heur				= heuristic(solution( route, i - 1), solution( route,i));
					temp_start_times(i) = tempst + serv + heur;
					
					%Check time windows. We're shifting everything to an earlier time, so only need to check beginning of time windows.
					if temp_start_times(i) < tw(1,solution(route,i))
						temp_start_times(i) = tw(1,solution(route,i));
					end
					
					%Increment
					i = i+1;
					
					%Make sure we don't overrun the matrix.
					if i > size(solution,2)
						break;
					end
				end
				start_times(route,:) = temp_start_times;
				
				%Update loads.
				temp_loads(index1)	= temp_loads(index1) + temp_loads(index2);
				temp_loads(index2)	= [];
				temp_loads			= [temp_loads,0]; %#ok
				loads(route,:)		= temp_loads;
				
				%Update values.
				temp_rv(index1)		= temp_rv(index1) + temp_rv(index2);
				temp_rv(index2)		= [];
				temp_rv				= [temp_rv,0]; %#ok
				all_value(route,:)	= temp_rv;
				if route <= num_vehicles
					for i1 = size(temp_rv,2):-1:1
						if temp_rv(i1) == 0
							temp_rv(i1) = [];
						else
							break;
						end
					end
					route_value(route,:) = zeros(1,size(route_value,2));
					
					for i2 = 1:size(temp_rv,2)
						route_value(route,i2) = temp_rv(i2);
					end
				end
			end
		end
	end
end


end



%
%
