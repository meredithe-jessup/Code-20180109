
function [solution, cost, route_cost, start_times, loads, route_dist, value, route_value, all_value, change] ...
					= Or_opt_v2_value(	solution, cost, route_cost, heuristic, heuristic2, tw, service, Q, start_times, ...
										loads, dist_limit, num_vehicles, route_dist, value, route_value, all_value)
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