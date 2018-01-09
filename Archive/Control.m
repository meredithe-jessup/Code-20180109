function [full_solution]=Control(d_dist, p_dist, threat, num_vehicles, current_inventory, theta, DispLoop)


%--------------------------------------------------------------------------
%-----Decoder ring for the input-----%
% All inputs are assuming that the depot is listed last (i.e., row/column n
% for each entry is the depot).
%
% d_dist = a matrix of the physical distances between all of the nodes. 
% This is used to ensure that the vehicle does not exceed its range. Also 
% used in conjunction with veh_speed parameter to determine time info.
%
% p_dist = a matrix with the probability of survival when traversing
% between each of the nodes.
%
% threat = 1xn vector giving the threat at each customer and the depot.
% 
% num_vehicles = the total number of vehicles allotted for the current time
% step. The algorithm will solve the problem for 1 vehicle, 2 vehicles, ...
% up to the num_vehicles.
%
% current_inventory = 1xn vector, the inventory of each customer at the 
% beginning of the current time step.
%
% theta = 1x3 vector, the parameters for the quadratic reward curve.
% 
% DispLoop = toggle Loop counter output
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%-----Decoder ring for the output-----%
% When you open the solution, you'll see 4 entries.
%  Entry 1: A cell array with all of the output data layered within it.  
%  This entry is further detailed below.
%
%  Entry 2: This is the demand for the customers.  Columns correspond to
%  customers (i.e., column j = customer j) and rows correspond the number
%  of vehicles (i.e., row i = i vehicles).  So entry i,j says if I allot i
%  vehicles, then customer j should receive this delivery quantity.
%
%  Entry 3: This table shows how the value of the solution changed as the
%  algorithm proceeded.  Rows correspond to # of vehicles.  So looking at
%  row i, you can see what the value is for each iteration of the algorithm
%  as the demands are altered. This info should give you an idea of how
%  many non-improving moves to allow.
%
%  Entry 4: This is the total run time of the routing portion of the
%  algorithm in seconds.
%
%
%--------------------------------------------------------------------------
% Now, I'll explain Entry 1 from above in further detail.  It is a nx1 cell
% array where row i corresponds to an allotment of i vehicles.  Inside of
% each of these cells is another cell array, 1x11, that contains all of the
% data for the i vehicle solution.
%
%
% Cell 1 = The solution.  In other words, given i vehicles, these are the
% routes that you should give to your vehicles.
%
% Cell 2 = The total distance required to fulfill the solution.
%
% Cell 3 = The time that service begins for each of the customers. The 
% entries correspond to Cell 1, the solution.  In other words, for each
% customer from the solution in Cell 1, there is a corresponding start time
% in this matrix.  This cell will give you the total time for each route. 
% The last non-zero value in each row is the time the vehicle arrives back at the depot.
% 
% Cell 4 = This is the delivery amount corresponding to each customer in
% the solution given in Cell 1.
%
% Cell 5 = This is the amount of time the algorithm required to solve the
% problem given i vehicles.
%
% Cell 6 = Same as Cell 5, just calculated slightly differently. Should be 
% within a couple seconds of each other.
% 
% Cell 7 = The distance covered by each route given in Cell 1.
% 
% Cell 8 = The total value added by executing this solution.
%
% Cell 9 = The individual values of each delivery.
%
% Cell 10 = The probabilities of survival for each leg of a route.  I.e.,
% cell i,j gives the probability of arriving at the customer given in the
% i,j cell of the solution from Cell 1 after departing the customer from
% the i,j-1 cell of the solution from Cell 1.
%
% Cell 11 = The cumulative probabilites for each route. I.e., cell i,j is
% simply the product of cells i,1 , i,2 , ... , i,j-1.
%--------------------------------------------------------------------------


%Start the timers.
tic;
timer_start = cputime;

if nargin < 7
    DispLoop = false;
end


%--------------------------------------------------------------------------
% This section are all the parameters that you can change to suit your
% problem.  You can also easily change any of these to an input if you'd
% like.  I hard-coded them for simplicity for now.

%The algorithm will treat deliveries as single units. If we want to use 
%something else, this scaling factor will automatically adjust everything. 
%I.e., if you have an 8000 lb limit but only want to consider changes in
%500 lb increments, set the scaling factor to 500 and the routing algorithm
%will then run on a scale of 0-16 in 1 unit increments, converting back to
%the original scale at the end.
scaling_factor = 500;

%This tells us how many non-improving moves when we're altering the demand
%before we call it done.
non_improv_limit = 15; 

%1xn vector, with service time for each customer (including a 0 for depot)
service = zeros(1,size(d_dist,1));

%Vehicle capacity.
q = 8000;

%The range of the vehicle.
dist_limit = 494;

%Speed of the vehicle. Ensure units are consistent.
veh_speed = 148;

%Capacity of each customer. Currently assumed homogeneous capacities.
Customer_cap = 8000;

%Length of the time period for the problem.  Can define as whatever (hours,
%minutes, days) but make sure that the veh_speed is input as
%distance/period.  I.e., if dist_limit is in miles and period is in hours,
%make sure speed is miles per hour.  If dist_limit is in parsecs and period
%is in years, make sure veh_speed is parsecs per year.
period = 6;

%If not using time windows, define a dummy time window vector that spans
%the entire period.
tw = zeros(2,size(d_dist,1));
for i = 1:size(tw,2)
    tw(2,i) = period;%set to period length
end
%--------------------------------------------------------------------------


%Fix inputs to move the depot info to the front. My code uses customer1 as
%the depot. Ian's uses customer37 as the depot.
%Fix distances.
x = size(d_dist,1);
temp = d_dist(x,:);
d_dist(x,:) = [];
d_dist = [temp;d_dist];
temp = d_dist(:,x);
d_dist(:,x) = [];
d_dist = [temp,d_dist];

x = size(p_dist,1);
temp = p_dist(x,:);
p_dist(x,:) = [];
p_dist = [temp;p_dist];
temp = p_dist(:,x);
p_dist(:,x) = [];
p_dist = [temp,p_dist];


%Fix threat
if size(threat,2) < x
    threat = [0,threat];
else
    temp = threat(1,x);
    threat(x) = [];
    threat = [temp,threat];
end

%Fix current_inventory
if size(current_inventory,2) < x
    current_inventory = [0,current_inventory];
else
    temp = current_inventory(1,x);
    current_inventory(x) = [];
    current_inventory = [temp,current_inventory];
end

%Initialize.
full_solution = cell(3,1);
value_solution = cell(num_vehicles,2);

%Incorporate scaling factor.
Q = q/scaling_factor;
customer_capacity = Customer_cap/scaling_factor;

%ACO parameters.
beta     = 5;
run_time = 0;
pbest    = 0.95;
qo       = 0.9;
rho      = 0.7;
gp_ratio = 0;
num_ants = 36;
phi      = 0;

%Solve for the maximum point on the reward curve defined by the thetas
%input.
global CONSTANTS
CONSTANTS.Theta = theta;
myOptions=optimset('Display','off','Algorithm','active-set');
max_value = fmincon(@myfun,0,0,0,0,0,0,customer_capacity*scaling_factor,[],myOptions);
max_value = round(max_value/scaling_factor,0);

%Need to rescale the current inventory and thetas to match the scaling
%factor.
current_inventory = current_inventory./scaling_factor;
theta(2) = theta(2)/scaling_factor;
theta(3) = theta(2)/(-2*max_value);

%Initialize.
demand = zeros(num_vehicles,size(d_dist,1));
values = [];
    
%Calculate the total value of the new solution.
for i = 1:num_vehicles
    if DispLoop
        display('Loop 1')
        i
    end
    %Define an initial demand. This defines a proportional demand for each
    %customer who is below their max_value.
    for cust = 2:size(demand,2)
        if current_inventory(cust) < max_value
            demand(i,cust) = max_value - current_inventory(cust);
        end
    end
    demand(i,:) = ceil(demand(i,:).*((i*Q)/(sum(demand(i,:)))));
            
    %Run a check to make sure we won't exceed anyone's demand.
    for cust = 2:size(demand,2)
        if demand(i,cust) > (customer_capacity - current_inventory(cust))
            demand(i,cust) = (customer_capacity - current_inventory(cust));
        end
    end
    
    %Run another check. If a customer is not allowing split deliveries, then
    %do not allow its demand to exceed the capacity of a single vehicle, Q.
    for cust = 2:size(demand,2)
        if (threat(cust) == 1) && (demand(i,cust) > Q)
            demand(i,cust) = Q;
        end
    end    
    
    %Run ACO to get an initial solution.
    ACO_solution = Copy_of_ACO_value(d_dist, p_dist, demand(i,:), tw, service, Q, beta, run_time, pbest, qo, rho, gp_ratio, num_ants, phi, dist_limit, veh_speed, threat, i, theta, current_inventory);

    values(i,1) = ACO_solution{8};

    %Assign the solution to the overall highest value solution.
    value_solution{i,1} = ACO_solution;
end


%For each number of vehicles, we want to start tweaking the demands to see
%if we can increase the overall value of the solution.
for i = 1:num_vehicles
    if DispLoop
        display('Loop 2')
        i
    end
    %Store the value of the current best solution.
    current_value = value_solution{i,1}{1,8};

    %Set up a while loop to go until we meet some ending condition
    trigger = 0;
        
    %Build a list of all customers.
    order_increase = zeros(1,size(demand,2)-1);
    for j = 1:(size(demand,2)-1)
        order_increase(1,j) = j+1;
    end
    order_decrease = order_increase;

    %Only choose a demand increase from those that aren't already full.
    %Remove customers who are full from the list of increase eligibles.
    for j = size(order_increase,2):-1:1
        if (demand(i,order_increase(j)) + current_inventory(order_increase(j))) == customer_capacity
            order_increase(j) = [];
        end
    end

    %Only choose a demand decrease from those that aren't already empty.
    %Remove customers who aren't receiving a delivery from the list of 
    %decrease eligibles.
    for j = size(order_decrease,2):-1:1
        if demand(i,order_decrease(j)) == 0
            order_decrease(j) = [];
        end
    end

    %If no eligible increases, force a decrease.
    if size(order_increase,2) == 0
        choice1 = 2;
    %If no eligible decreases, force an increase.
    elseif size(order_decrease,2) == 0
        choice1 = 1;
    %If no eligible increases or decreases, then we're done. Iterate the
    %for loop.
    elseif (size(order_increase,2)==0) && (size(order_decrease,2) == 0)
        continue;
    else
        choice1 = randperm(2,1);
    end
    
    counter = 2;

    %Repeat the following until we've reached the non-improving iteration
    %limit.
    while trigger < non_improv_limit

        if trigger > 0
            %Build a list of all customers.
            order_increase = zeros(1,size(demand,2)-1);
            for j = 1:(size(demand,2)-1)
                order_increase(1,j) = j+1;
            end
            order_decrease = order_increase;

            %Only choose a demand increase from those that aren't already full.
            %Remove customers who are full from the list of increase eligibles.
            for j = size(order_increase,2):-1:1
                if (demand(i,order_increase(j)) + current_inventory(order_increase(j))) == customer_capacity
                    order_increase(j) = [];
                end
            end

            %Only choose a demand decrease from those that aren't already empty.
            %Remove customers who aren't receiving a delivery from the list of 
            %decrease eligibles.
            for j = size(order_decrease,2):-1:1
                if demand(i,order_decrease(j)) == 0
                    order_decrease(j) = [];
                end
            end

            %If no eligible increases, force a decrease.
            if size(order_increase,2) == 0
                choice1 = 2;
            %If no eligible decreases, force an increase.
            elseif size(order_decrease,2) == 0
                choice1 = 1;
            %If no eligible increases or decreases, then we're done. Iterate the
            %for loop.
            elseif (size(order_increase,2)==0) && (size(order_decrease,2) == 0)
                continue;
            else
                choice1 = randperm(2,1);
            end
        end

        %Store the demand in a temporary array so we can alter it to see if
        %it yields an improved solution.
        temp_demand = demand(i,:);

        %Increase a demand.
        if choice1 == 1
            %Of the non-full demands, choose one.
            choice2 = randperm(size(order_increase,2),1);
            temp_demand(1,order_increase(choice2)) = temp_demand(1,order_increase(choice2)) + 1;
            if (temp_demand(1,order_increase(choice2)) + current_inventory(order_increase(choice2))) > customer_capacity
                trigger = trigger + 1;
                continue;
            end
        
        %Decrease a demand.
        elseif choice1 == 2
            %Of the non-full demands, choose one.
            choice2 = randperm(size(order_decrease,2),1);
            temp_demand(1,order_decrease(choice2)) = temp_demand(1,order_decrease(choice2)) - 1;
            if temp_demand(1,order_decrease(choice2)) < 0
                trigger = trigger + 1;
                continue;
            end
        end
        
        %Run another check. If a depot is not allowing split deliveries, then
        %do not allow its demand to exceed the capacity of a single vehicle, Q.
        if choice1 == 1
            if (threat(1,order_increase(choice2)) == 1) && (demand(i,order_increase(choice2)) > Q)
                demand(i,order_increase(choice2)) = Q;
            end
        end

        %Run ACO to get new solution.
        ACO_solution2 = Copy_of_ACO_value(d_dist, p_dist, temp_demand, tw, service, Q, beta, run_time, pbest, qo, rho, gp_ratio, num_ants, phi, dist_limit, veh_speed, threat, i, theta, current_inventory);

        %Store the value of the new solution.
        new_value = ACO_solution2{8};
        values(i,counter) = current_value;
        counter = counter + 1;
        
        %Compare the new solution to the previous best solution.
        if (new_value > current_value) %|| ((new_value == current_value) && (choice1 == 2))
            %Reset trigger to 0 if new solution is better and copy this new
            %solution.
            trigger = 0;
            value_solution{i,1} = ACO_solution2;
            current_value = new_value;
            demand(i,:) = temp_demand;
        else
            %Increment the trigger to show that we had a
            %non-improving move.
            trigger = trigger + 1;
        end

    end
    
    
%Note: thought this was a good idea, couldn't get it to work. I was merely
%trying to ensure that each vehicle that departed the depot was full.
% %--------------------------------------------------------------------------
% %New addition to code. If something's broken, it's probably in here.
% 
% %Add in code to augment any solutions that don't use the full allotment
%     %of vehicles, if such an augmentation exists.
%     %Find the customer with the lowest inventory and add him to a new
%     %route.
%     s = value_solution{i,1};
%     loads = value_solution{i,4};
%     rv = value_solution{i,14};
%     st = value_solution{i,3};
%     rdist = value_solution{i,12};
%     if size(s,1) < i
%         
%         %Add as many routes as necessary, or as many as we can.
%         for veh = 1:(i-size(s,1))
%             %Determine which customers have some unfilled demand.
%             filled_demand = zeros(1,size(demand,1));
%             for i1 = 1:size(s,1)
%                 for i2 = 2:size(s,2)
%                     if s1(i1,i2) == 1
%                         break;
%                     end
%                     filled_demand(1,s(i1,i2)) = filled_demand(1,s(i1,i2)) + loads(i1,i2);
%                 end
%             end
% 
%             %Calculate the remaining demand.
%             remaining_demand = demand(i,:) - filled_demand;
% 
%             %Now, pick the customer with the largest value.
%             %For now, pick customer with lowest inventory and delivery his
%             %remaining demand.
%             
%             %Initialize new route.
%             addon_s = zeros(1,size(s,2));
%             addon_s(1) = 1;
%             
%             %Initialize new route_value.
%             addon_rv = zeros(1,size(rv,2));
%             
%             %Initialize new loads.
%             addon_loads = zeros(1,size(loads,2));
%             
%             %Initialize new start times.
%             addon_st = zeros(1,size(st,2));
%             
%             %Initialize new route distances.
%             addon_rdist = 0;
%             
%             condition = 0;
%             while condition == 0
% 
%                 %Find the most remaining demand.
%                 [x,y] = min(remaining_demand);
% 
%             end
%         end
%     end
% %--------------------------------------------------------------------------

end

%Re-scale the demands and loads arrays back to the original unit scale 
%using scaling factor.
for i = 1:num_vehicles
    temp = value_solution{i,1}{1,4};
    temp = temp.*scaling_factor;
    value_solution{i,1}{1,4} = temp;
end
demand = demand.*scaling_factor;

%Output the solution.
full_solution{1} = value_solution;
full_solution{2} = demand;
full_solution{3} = values;
full_solution{4} = cputime - timer_start;

%--------------------------------------------------------------------------
%-----Decoder ring for the input-----%
% All inputs are assuming that the depot is listed last (i.e., row/column n
% for each entry is the depot).
%
% d_dist = a matrix of the physical distances between all of the nodes. 
% This is used to ensure that the vehicle does not exceed its range. Also 
% used in conjunction with veh_speed parameter to determine time info.
%
% p_dist = a matrix with the probability of survival when traversing
% between each of the nodes.
%
% threat = 1xn vector giving the threat at each customer and the depot.
% 
% num_vehicles = the total number of vehicles allotted for the current time
% step. The algorithm will solve the problem for 1 vehicle, 2 vehicles, ...
% up to the num_vehicles.
%
% current_inventory = 1xn vector, the inventory of each customer at the 
% beginning of the current time step.
%
% theta = 1x3 vector, the parameters for the quadratic reward curve.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%-----Decoder ring for the output-----%
% When you open the solution, you'll see 4 entries.
%  Entry 1: A cell array with all of the output data layered within it.  This
%  entry is further detailed below.
%
%  Entry 2: This is the demand for the customers.  Columns correspond to
%  customers (i.e., column j = customer j) and rows correspond the number
%  of vehicles (i.e., row i = i vehicles).  So entry i,j says if I allot i
%  vehicles, then customer j should receive this delivery quantity.
%
%  Entry 3: This table shows how the value of the solution changed as the
%  algorithm proceeded.  Rows correspond to # of vehicles.  So looking at
%  row i, you can see what the value is for each iteration of the algorithm
%  as the demands are altered. This info should give you an idea of how
%  many non-improving moves to allow.
%
%  Entry 4: This is the total run time of the routing portion of the
%  algorithm in seconds.
%
%
%--------------------------------------------------------------------------
% Now, I'll explain Entry 1 from above in further detail.  It is a nx1 cell
% array where row i corresponds to an allotment of i vehicles.  Inside of
% each of these cells is another cell array, 1x11, that contains all of the
% data for the i vehicle solution.
%
%
% Cell 1 = The solution.  In other words, given i vehicles, these are the
% routes that you should give to your vehicles.
%
% Cell 2 = The total distance required to fulfill the solution.
%
% Cell 3 = The time that service begins for each of the customers. The 
% entries correspond to Cell 1, the solution.  In other words, for each
% customer from the solution in Cell 1, there is a corresponding start time
% in this matrix.  This cell will give you the total time for each route. 
% The last non-zero value in each row is the time the vehicle arrives back at the depot.
% 
% Cell 4 = This is the delivery amount corresponding to each customer in
% the solution given in Cell 1.
%
% Cell 5 = This is the amount of time the algorithm required to solve the
% problem given i vehicles.
%
% Cell 6 = Same as Cell 5, just calculated slightly differently. Should be 
% within a couple seconds of each other.
% 
% Cell 7 = The distance covered by each route given in Cell 1.
% 
% Cell 8 = The total value added by executing this solution.
%
% Cell 9 = The individual values of each delivery.
%
% Cell 10 = The probabilities of survival for each leg of a route.  I.e.,
% cell i,j gives the probability of arriving at the customer given in the
% i,j cell of the solution from Cell 1 after departing the customer from
% the i,j-1 cell of the solution from Cell 1.
%
% Cell 11 = The cumulative probabilites for each route. I.e., cell i,j is
% simply the product of cells i,1 , i,2 , ... , i,j-1.
%--------------------------------------------------------------------------

end