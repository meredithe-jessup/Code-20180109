load 'Baseline4'

for j = 1:37
    for i = 1:37
        if i > j
            temp1 = MapBlue.ShortestPaths.Paths{i,j};
            temp2 = fliplr(temp1);
            MapBlue.ShortestPaths.Paths{i,j} = temp2;
        end
    end
end

clear SPP i j ans temp temp1 temp2