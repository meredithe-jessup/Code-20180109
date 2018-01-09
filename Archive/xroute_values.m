%This function will calculate the value of each delivery given a solution.

function val = route_values(routes, loads, current_inventory, theta)

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