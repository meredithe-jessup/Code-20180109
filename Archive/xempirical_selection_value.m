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