function [solution, cost, route_cost, start_times, loads, route_dist, value, route_value, all_value, change] ...
	= two_opt_star_value(solution, cost, route_cost, heuristic, heuristic2, tw, service, ...
							Q, start_times, loads, dist_limit, route_dist, num_vehicles, value, route_value, all_value, current_inventory, theta)
%Given an initial solution, this program will run 2-opt* to try to find an improved solution.

						
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