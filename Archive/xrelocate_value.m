%This function is an intra-route operator, greedily moving customers within a route if a better feasible solution is found.
function [solution, cost, route_cost, start_times, loads, value, route_value, change] ...
				= relocate_value(solution, cost, route_cost, heuristic, tw, service, ~, start_times, loads, value, route_value, num_vehicles)

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