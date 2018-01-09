
function [unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, pheromone, start_times, loads, loc, route_dist] ...
						= construct_route_value(loc, ploc, unvisited, route_feasible, s, cost, time, route_cost, vcap, r, v, s_demand, heuristic, ... 
												heuristic2, tw, service, pheromone, min, max, Q, rho, start_times, loads, phi, route_dist, threat)

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