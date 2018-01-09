%This code is for a simple VRP with normal constraints. Inputs are:
%heuristic matrix (square), demand array (1xn), max vehicle load,
%heuristic importance, pheromone importance, pheromone evaporation rate,
%min pheromone level, and max pheromone level.

%Inputs:
%heuristic	= matrix of distances/costs between nodes (symmetric)
%demand		= array of demands for each node
%Q			= max vehicle load
%alpha		= relative importance of heuristic information
%beta		= relative importance of pheromone information
%ro			= pheromone evaporation rate
%min		= min pheromone level
%max		= max pheromone level

%Potential inputs not currently included
%q0			= acceptance criteria for city selection (Ant System modification)
%ants		= number of ants (to include just introduce as input and use as
%       index for outermost for loop
% function [cost, s, route_cost, pheromone, start_times, loads, route_dist]...
% 	= ant_value(heuristic, heuristic2, candidate_list, cl_size, pheromone, demand, tw, service, Q, ...
% 				alpha, beta, min, max, hx, qo, phi, dist_limit, threat, theta, current_inventory)
function ant = ant_value(heuristic, heuristic2, candidate_list, cl_size, pheromone, demand, tw, service, Q, ...
				alpha, beta, min, max, hx, qo, phi, dist_limit, threat, theta, current_inventory)

%Initialize unvisited_array to entire customer list
unvisited = ones(hx-1,1);
for j = 1:hx-1
    unvisited(j,1) = j+1;
end

%Remove any customers with 0 demand from the unvisited array.
for j = size(unvisited,1):-1:1
    if demand(1,unvisited(j)) == 0
        unvisited(j) = [];
    end
end

s_demand	= demand;	%Initialize demand array for this solution
s			= 1;		%Initialize current solution.
v			= 1;		%Initialize vehicle counter
vcap		= Q;		%Initialize remaining vehicle capacity
r			= 1;		%Initialize number of nodes in current route
loc			= 1;		%Initialize current location to depot
cost		= 0;		%Initialize total cost

%Initialize time required by current trip (i.e., current departure from depot).  This value resets every time the depot is visited to simulate
%having multiple vehicles.  It is the cost incurred by travel distances plus waiting times.
time			= 0;

%Initialize cost for current route.  Note that time and route_cost are the same except time includes any waiting periods whereas route_cost doesn't.
route_cost		= 0;
route_dist		= 0;
route_feasible	= unvisited; % Initialize list of customers that can be feasibly added to current route
start_times		= [];		 % Initialize the start times for service at each stop.
loads			= [];		 % Initialize the load that each truck delivers to each customer.

%Build routes
while numel(unvisited) > 0
    
    
    ploc  = loc; %Remember previous location
    r     = r+1; %Add another node to current route

    %Subroutine to update list of feasible locations for current route
    temp  = zeros(hx,1);
    place = 0;
	
    %Reason for the if statement: if r=2, then we are at the depot and the entire set is feasible so we don't need to update feasibility.
    if (r > 2)  && (numel(route_feasible) > 0)
        for j = 1:size(route_feasible)
            %check to ensure time window compatibility
            if time + service(ploc) + heuristic(ploc,route_feasible(j)) <= tw(2,route_feasible(j))
                %check to ensure distance limit compatibility
                if route_dist + heuristic2(ploc,route_feasible(j)) + heuristic2(route_feasible(j),1) <= dist_limit
                    place		= place + 1;
                    temp(place) = route_feasible(j);
                end
            end
        end
        
        route_feasible = temp((1:place),1);
    end
    
    %Get next location from empirical_selection subroutine
    loc			= empirical_selection_value(ploc, candidate_list, cl_size, route_feasible, heuristic, pheromone, ...
											alpha, beta, qo, loads, threat, s_demand, Q, theta, current_inventory);

    %Run construct_route subroutine to add next location to route
    [unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, pheromone, start_times, loads, loc, route_dist] ...
				= construct_route_value(	loc, ploc, unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, heuristic, ...
											heuristic2, tw, service, pheromone, min, max, Q, start_times, loads, phi, route_dist);

    %If depot is the only unvisited location, add edge from current location back to depot to complete routing
    if (numel(unvisited) == 0) && (loc ~= 1)
        s(v,r+1)		= 1;
        cost			= cost + heuristic(loc,1);				%Update cost
        route_cost(v)	= route_cost(v) + heuristic(loc,1);		%Update route_cost
        route_dist(v)	= route_dist(v) + heuristic2(loc,1);	%Update route_cost 
        
	end
end

ant = struct(	'cost',			cost,		...
				's',			s,			...		
				'route_cost',	route_cost, ...
				'pheromone',	pheromone,	...
				'start_times',	start_times,...
				'loads',		loads,		...
				'route_dist',	route_dist);
end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
function [unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, pheromone, start_times, loads, loc, route_dist] ...
						= construct_route_value(loc, ploc, unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, heuristic, ... 
												heuristic2, tw, service, pheromone, min, max, Q, start_times, loads, phi, route_dist)

%If depot is the next location, then proceed to depot and start a new route.
if loc == 1

    s(v,r)				= 1;									%Add depot to current route
    cost				= cost + heuristic(ploc,1);				%Update cost to reflect return to depot
    route_cost(v)		= route_cost(v) + heuristic(ploc,1);	%Update route costs.
    route_dist(v)		= route_dist(v) + heuristic2(ploc,1);	%Update route distances.
    start_times(v,r)	= time + heuristic(ploc,1);				%Update the service start times.
    loads(v,r)			= 0;									%Update loads.
    v					= v + 1;								%Add another vehicle
    
    %Add new row to s, start_times, and loads for the new vehicle
    s					= [s;zeros(1,size(s,2))];
    start_times			= [start_times;zeros(1,size(start_times,2))];
    loads				= [loads;zeros(1,size(loads,2))];
    route_cost			= [route_cost;0];
    route_dist			= [route_dist;0];
    
    s(v,1)				= 1;							%Set first location in new route as depot
    start_times(v,1)	= 0;							%Set the start time as zero
    loads(v,1)			= 0;							%Set the load as zero
    r					= 1;							%Reset # of nodes on current route to 1.
    time				= 0;							%Reset trip time to zero since we're back at the depot.
    route_feasible		= unvisited;					%Reset list of feasible locations to all unvisited locations
    vcap				= Q;							%Reset vehicle capacity to Q
    pheromone(loc,ploc) = pheromone(loc,ploc)*(1-phi);	%Local pheromone update
    
    %Check pheromones against min/max
    if pheromone(loc,ploc) < min
        pheromone(loc,ploc) = min;
        
    elseif pheromone(loc,ploc) > max
        pheromone(loc,ploc) = max;
        
    end
    
    pheromone(ploc,loc) = pheromone(loc,ploc);			%Maintain symmetry of pheromone matrix
    
%If next location is not the depot, then attempt to add it to the route.
else
    %If demand can fit on current vehicle, add all demand to current route
    if s_demand(loc) < vcap

        s(v,r) = loc; %Add current location to current route
        
        %Remove loc from list of unvisited cities
        for k=1:size(unvisited)
            if unvisited(k) == loc
                unvisited(k) = [];
                break;
            end
        end
        
        %Remove loc from list of feasible locations
        for k=1:size(route_feasible)
            if route_feasible(k) == loc
                route_feasible(k) = [];
                break;
            end
        end
        
        vcap				= vcap - s_demand(loc);					 %Update vehicle capacity
        cost				= cost + heuristic(ploc,loc);			 %Update total routing cost
        
        %Update trip time.
        time				= time + heuristic(ploc,loc);			 %Add route time.
		
        if (time < tw(1,loc)), time = tw(1,loc); end				 %Add waiting time, if applicable.

        start_times(v,r)	= time;									 %Update the service start times.
        time				= time + service(1,loc);				 %Add service time to total time.
        loads(v,r)			= s_demand(loc);						 %Update loads
        route_cost(v)		= route_cost(v) + heuristic(ploc,loc);	 %Update cost of current route
        route_dist(v)		= route_dist(v) + heuristic2(ploc,loc);  %Update dist of current route
        pheromone(loc,ploc) = pheromone(loc,ploc)*(1-phi);			 %Local pheromone update
        
        if pheromone(loc,ploc) < min, pheromone(loc,ploc) = min; end %Check pheromones against min/max
        
        pheromone(ploc,loc) = pheromone(loc, ploc);					 %Maintain symmetry of pheromone matrix
        
    %If demand is exactly the remaining capacity, add it to the current route, return to depot and start a new route.
    elseif s_demand(loc) == vcap
        
        s(v,r)		= loc;								%Add current location to route
        loads(v,r)	= s_demand(loc);					%Update loads
        
        %Update trip time.
        time		= time + heuristic(ploc,loc);		%Add route time.
		
        if time < tw(1,loc), time = tw(1,loc);end		%Add waiting time, if applicable.
        
        start_times(v,r)	= time;						%Update the service start times.
        time				= time + service(1,loc);	%Add service time to total time.
        s(v,r+1)			= 1;						%Add depot to current route because vehicle is capacitated
        loads(v,r+1)		= 0;						%Update loads.
    
        %Update trip time.
        time				= time + heuristic(loc,1);	%Add route time.
        start_times(v,r+1)	= time;						%Update the service start times.
        
        %Remove loc from list of unvisited cities
        for k=1:size(unvisited)
            if unvisited(k) == loc
                unvisited(k) = [];
                break
            end
        end
        
        cost			= cost + heuristic(ploc,loc) + heuristic(loc,1);			%Update cost to reflect return to depot
        time			= 0;														%Reset trip time to zero since we're back at the depot.
        route_cost(v)	= route_cost(v) + heuristic(ploc,loc) + heuristic(loc,1);	%Update route cost to zero since we're back at the depot
        route_dist(v)	= route_dist(v) + heuristic2(ploc,loc) + heuristic2(loc,1); %Update route dist to zero since we're back at the depot
        route_feasible	= unvisited;												%Update feasible locations to all unvisited

        if numel(unvisited) > 0
            
            v = v + 1; %Add another vehicle

            %Add new row to s, start_times, and loads for the new vehicle
            s			= [s;zeros(1,size(s,2))];
            start_times = [start_times;zeros(1,size(start_times,2))];
            loads		= [loads;zeros(1,size(loads,2))];
            route_cost	= [route_cost;0];
            route_dist	= [route_dist;0];

            s(v,1)				= 1; %Make first location the depot
            start_times(v,1)	= 0; %Set the start time as zero
            loads(v,1)			= 0; %Set the load as zero
            r					= 1; %Update number of nodes in current route
            vcap				= Q; %Update vehicle capacity
        end
        
        %Local pheromone update
        pheromone(loc,ploc)		= pheromone(loc,ploc)*(1-phi);
        pheromone(loc,1)		= pheromone(loc,1)*(1-phi);
        
        %Check pheromones against min/max. Since 0 <= phi < 1, only need to check min because pheromone can only decrease.
        if pheromone(loc,ploc)	< min, pheromone(loc,ploc)	= min; end
        if pheromone(loc,1)		< min, pheromone(loc,1)		= min; end
        
        %Maintain symmetry of pheromone matrix
        pheromone(ploc,loc)		= pheromone(loc, ploc);
        pheromone(1,loc)		= pheromone(loc, 1);
        
        loc						= 1; %Set current location to depot
    
    %If demand cannot fit on current vehicle, fit as much as possible, update demand array, return to depot and start a new route.  Note that the 
	%location is left on the unvisited list because its demand is not completely filled.
    else     
        
        s(v,r)			= loc;					%Add current location to route
        s(v,r+1)		= 1;					%Add depot to current route because vehicle is capacitated
        loads(v,r)		= vcap;					%Update loads
        s_demand(loc)	= s_demand(loc)-vcap;	%Update demand at current location to reflect partial fulfillment
        
        
        cost			= cost + heuristic(ploc,loc)		   + heuristic(loc,1);	%Update cost to reflect return to depot
        route_cost(v)	= route_cost(v) + heuristic(ploc,loc)  + heuristic(loc,1);	%Update route cost to zero since we're back at the depot
        route_dist(v)	= route_dist(v) + heuristic2(ploc,loc) + heuristic2(loc,1); %Update route dist to zero since we're back at the depot
        loads(v,r+1)	= 0;														%Update loads.
    
        %Update trip time.
        time						= time + heuristic(ploc,loc);	%Add route time.
        if time < tw(1,loc), time	= tw(1,loc); end				%Add waiting time, if applicable.
        start_times(v,r)			= time;							%Update the service start times.
        time						= time + service(1,loc);		%Add service time to total time.
        s(v,r+1)					= 1;							%Add depot to current route because vehicle is capacitated

        %Update trip time.
        time						= time + heuristic(loc,1);		%Add route time.
        start_times(v,r+1)			= time;							%Update the service start times.
        time						= 0;							%Reset trip time to zero since we're back at the depot.
        route_feasible				= unvisited;					%Update feasible locations to all unvisited
        v							= v + 1;						%Add another vehicle
        
        %Add new row to s, start_times, and loads for the new vehicle
        s					= [s;zeros(1,size(s,2))];
        start_times			= [start_times;zeros(1,size(start_times,2))];
        loads				= [loads;zeros(1,size(loads,2))];
        route_cost			= [route_cost;0];
        route_dist			= [route_dist;0];
        
        s(v,1)				= 1; %Make first location the depot
        start_times(v,1)	= 0; %Set the start time as zero
        loads(v,1)			= 0; %Set the load as zero
        r					= 1; %Update number of nodes in current route
        vcap				= Q; %Update vehicle capacity
        
        %Local pheromone update
        pheromone(loc,ploc) = pheromone( loc,ploc)*(1-phi);
        pheromone(loc,1)	= pheromone( loc,1)*(1-phi);
        
        %Check pheromones against min/max
        if pheromone(loc,ploc)	< min, pheromone(loc,ploc)	= min; end
        if pheromone(loc,1)		< min, pheromone(loc,1)		= min; end
        
        %Maintain symmetry of pheromone matrix
        pheromone(ploc,loc) = pheromone(loc, ploc);
        pheromone(1,loc)	= pheromone(loc, 1);
        
        loc					= 1; %Set current location to depot
    
    end
end

end



% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
function loc = empirical_selection_value(	loc, candidate_list, cl_size, route_feasible, heuristic, pheromone, alpha, beta, qo, loads, threat, ...
											demand, Q, theta, current_inventory)

%Construct an updated feasible list consisting of the intersection of the feasible list and the candidate list for the current location.
feasible			= zeros(size(heuristic,1),1);
feasible_counter	= 0;
q					= rand(1,1);
ploc				= loc;

%current location is the depot, then everything remaining is feasible
if loc == 1
    feasible			= route_feasible;
    feasible_counter	= 1;

%otherwise, build the candidate list with only eligible candidates with a positive remaining load
elseif numel(route_feasible) > 0
    for k = 1:cl_size
        for j = 1:size(route_feasible,1)
            if candidate_list(loc,k) == route_feasible(j)
                feasible_counter			= feasible_counter + 1;
                feasible(feasible_counter)	= route_feasible(j);
                break;
            end
        end
    end
end

valid_loc = 0; %Binary variable saying that we have not yet chosen a valid location.

while valid_loc == 0
    
    if feasible_counter > 0, feasible = feasible(1:feasible_counter); end %???

    %If depot is the only unvisited node, send route to depot
    if feasible_counter == 0
        loc = 1;

    %If only 1 feasible location, go there.
    elseif (numel(feasible) == 1)
        loc = feasible(1);

    %Check q0 criterion. If q < q0, make a greedy move.
    elseif q < qo

        %Get cells corresponding to current city and unvisited cities from heuristic and pheromone matrices and multiply them together
        emp_dist = [feasible,zeros(size(feasible,1),1),zeros(size(feasible,1),1)];

        %Make an empirical distribution based on heuristic and pheromone info.
        for k = 1:size(feasible,1)
			%***
			vf1a			= value_function( demand( feasible( k)),			theta, 1);
			vf1b			= value_function( current_inventory( feasible( k)), theta, 1);
			tau1			= pheromone( loc,feasible( k));
            emp_dist(k,2)	= ((vf1a - vf1b)^beta)*(tau1^alpha);
        end

        %Find the maximum.
        position = 1;
        for i = 2:size(emp_dist,1)
            if emp_dist(i,2) > emp_dist(position(1),2)
                position = i;
            elseif emp_dist(i,2) == emp_dist(position(1),2)
                position = [position,i];
            end
        end

        %if a single maximum value, choose it as the next customer.
        if size(position, 2) == 1
            loc  = emp_dist(position,1);

        %if multiple optimals, randomly select one.
        else
            perm = randperm(size(position,2));
            loc  = emp_dist(position(perm(1)),1);
        end

    %Otherwise, empirically choose next location
    else
        %Get cells corresponding to current city and unvisited cities from heuristic and pheromone matrices and multiply them together
        emp_dist	= [feasible,zeros(size(feasible,1),1),zeros(size(feasible,1),1)];

        for k = 1:size(feasible,1)
			%***
			vf2a			= value_function( demand( feasible( k)) + current_inventory( feasible( k)), theta, 1);
			vf2b			= value_function( current_inventory( feasible( k)),						theta, 1);
			tau2			= pheromone( loc, feasible( k));
            emp_dist(k,2)	= (vf2a - vf2b)^beta * ( tau2^alpha);
        end

        %Rescale the second column to lie between 0 and 1.
        lo			= min(emp_dist(:,2));
        hi			= max(emp_dist(:,2));
		
        if hi == lo
            %check to make sure we won't get a "/0" error. In this case, give all equal probability.
            for k = 1:size(emp_dist,1)
                emp_dist(k,2) = 1/size(emp_dist,1);
            end
        else            
            %Rescale.
            for k = 1:size(emp_dist,1)
                emp_dist(k,2) = (emp_dist(k,2)-lo)/(hi-lo);
			end
			
            %Rescale again such that the sum of the second column is 1.
            total = sum(emp_dist(:,2));
            for k = 1:size(emp_dist,1)
                emp_dist(k,2) = emp_dist(k,2)/total;
            end
        end
        
        %Add a cumulative column
        emp_dist(1,3) = emp_dist(1,2);
        for k = 2:size(emp_dist,1)
            emp_dist(k,3)=emp_dist(k-1,3)+emp_dist(k,2);
        end
        
        p		= rand(1,1);
        breaker = 0;
        k		= 0;
        
        while breaker == 0
            k = k+1;
                if p<emp_dist(k,3)
                    loc		= emp_dist(k,1);
                    breaker = 1;
                end
        end
    end
    
    %check to make sure we aren't forcing a split delivery in a hi threat environment. If we are not, then return the current location. If we are, 
	%then remove the current location from the list of feasible locations and re-draw.
    if threat(loc) == 0
        valid_loc = 1;
    else
        if ploc == 1 || loc == 1
            valid_loc = 1;
        elseif sum(loads(size(loads,1),:)) + demand(loc) <= Q
            valid_loc = 1;
        else
            feasible_counter = feasible_counter - 1;
            for i = 1:size(feasible,2)
                if feasible(i) == loc, feasible(i) = []; end
            end
        end
    end
end

end