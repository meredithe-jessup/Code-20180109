function [solution, cost, route_cost, start_times, loads, route_dist] ...
		= fix_routes(solution, cost, route_cost, heuristic, heuristic2, heuristic3, tw, service, start_times, loads, route_dist)
%Fix ensures that the same customer isn't visited multiple times within a single route. If so, the multiple visits are combined into a single visit
%at the first location.

rows = size(solution,1);

%columns = total size of matrix
columns = size(solution,2);

%Iterate through each route
for route = 1:rows
    
    %For each route, check each of the customers to see if they occur
    %twice. Note that we don't need to check the depot.
    for index1 = 2:columns

        %If the customer directly after customer1 is the depot, then we're
        %done checking this route.
        if solution(route,index1+1) <= 1
            break;
        elseif (index1+1) > size(solution,2)
            break;
        end
        
        %For each customer1, check each customer that occurs after it.
        for index2 = index1+1:columns
        
            %If customer2 is the depot, then we're done checking the
            %customer1. Move on to next customer1.
            if solution(route,index2) <= 1
               break;
            end

            %Assign the customers.
            customer1 = solution(route,index1);
            customer2 = solution(route,index2);

            %Check to see if the customers are the same.
            if customer1 == customer2
                
                %If they are the same, then combine them into customer1's
                %location. Update solution, cost, route_cost, start_times,
                %and loads.
                temp = solution(route,:);
                temp_start_times = start_times(route,:);
                temp_loads = loads(route,:);
                for k = index2:(size(temp,2))
                    temp_start_times(1,k) = 0;
                end
                
                %Update solution.
                temp(index2:(size(temp,2)-1)) = temp((index2 + 1):size(temp,2));
                temp(size(temp,2)) = 0;
                
                %Update costs.
                route_change = heuristic3(solution(route,index2 - 1),solution(route,index2 + 1)) - heuristic3(solution(route,index2 - 1),solution(route,index2)) - heuristic3(solution(route,index2),solution(route,index2 + 1)); 
                cost = cost + route_change;
                route_cost(route) = route_cost(route) + route_change;
                
                %Update distance.
                route_change2 = heuristic2(solution(route,index2 - 1),solution(route,index2 + 1)) - heuristic2(solution(route,index2 - 1),solution(route,index2)) - heuristic2(solution(route,index2),solution(route,index2 + 1)); 
                route_dist(route) = route_dist(route) + route_change2;                

                solution(route,:) = temp;

                %Update start times. We've updated the solution, so index2
                %is no longer a duplicate. Nothing prior to index2 has
                %changed wrt start times, so just update everything from
                %there on.
                i = index2;
                while solution(route,i) > 0
                
                    %Update times.
                    temp_start_times(i) = temp_start_times(i - 1) + service(solution(route,i - 1)) + heuristic(solution(route,i - 1),solution(route,i));
                    
                    %Check time windows. We're shifting everything to an
                    %earlier time, so only need to check beginning of time
                    %windows.
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
                temp_loads(index1) = temp_loads(index1) + temp_loads(index2);
%                 temp_loads(index2:(size(temp,2)-1)) = temp((index2 + 1):size(temp,2));
%                 temp_loads(size(temp,2)) = 0;
                temp_loads(index2) = [];
                temp_loads = [temp_loads,0];
                loads(route,:) = temp_loads;
            end
        end
    end
end


end