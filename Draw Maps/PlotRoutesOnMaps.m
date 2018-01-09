% Main function to generate Maps
clear
close all
check1 = [];
check2 = [];

try
    check1 = MapBlue.DepCoord;
catch
    load 'Baseline5.mat'
end

try
    check2 = Solution.run_time;
catch
    load 'Solution2.mat'
end

% maps = [1]; %#ok

XFOB = MapBlue.FOBCoord(:,1);
YFOB = MapBlue.FOBCoord(:,2);

XDep = MapBlue.DepCoord(1);
YDep = MapBlue.DepCoord(2);

SPPs = MapBlue.ShortestPaths.Paths;

if GeneralSettings.GridSize >= 30
	sz = [];
elseif GeneralSettings.GridSize >= 10
	sz = 200;
else
	sz = 300;
end

for v = 1:4
    
	% Set figure layout
    MapName		  = ['Map ', num2str(v)];
	figure( 'Name',MapName,     'NumberTitle','off',    'units','inches',	 'position',[0,0,8,6.5]);   hold on
	
    AdjustPlotPosition;          
    
	MapGrid     = DrawHexMap(GeneralSettings.GridSize,  GeneralSettings.HexSize,    Options.Labeling,   [], Options.TogTicks);
	plotDep     = scatter(XDep,YDep,    sz,     'sk',   'MarkerFaceColor','g');
    plotFOB     = scatter(XFOB,YFOB,    sz/4,   'k',    'MarkerFaceColor','b');
    
    RedPatches  = MapRed(1).RedPatches;
    pHandle     = MapGrid.Handles;
    NumRed      = length(RedPatches);
    
    for i = 1:NumRed
        ptch = RedPatches(i);
        set(pHandle(ptch),  'FaceColor','flat',     'FaceVertexCData',[1 0.3 0.1],     'CDataMapping','scaled')
    end
    
    R = size(Solution.value_solution{v}.sbest,1 );
    
    for r = 1:R
        cmap              = colormap(winter(R*3) );
        
        Order             = Solution.value_solution{v}.sbest(r,:);
        Order(Order == 0) = [];
        Nlinks            = length(Order)-1;
        RoutePath         = [];
        
        for i = 1:Nlinks
            tmpPath       = []; %#ok
            a             = Order(i);
            b             = Order(i+1);
            
            if i > 1
                tmpPath    = cell2mat(SPPs(a,b) );
                tmpPath(1) = []; 
            else
                tmpPath    = cell2mat(SPPs(a,b) );
            end
            
            RoutePath      = [RoutePath tmpPath]; %#ok
        end
    end
    
    [Mdim,Ndim] = size(Solution.value_solution{v}.Loads);
    
    fs = reshape(Solution.value_solution{v}.sbest',1,Mdim*Ndim);
    ls = reshape(Solution.value_solution{v}.Loads',1,Mdim*Ndim);
    
    [bdx, mdx, ndx] = unique(fs);
    result = accumarray(ndx, ls);
    
    Fs = bdx;
    Ls = string(result);
    
    text(XFOB(Fs') + 1,YFOB(Fs,1) + 0.5,Ls);
    for t = 1:Nlinks+1
        plot(MapGrid.XY(RoutePath,1),MapGrid.XY(RoutePath,2),'Color',cmap(r*3,:),'LineWidth', 2);
    end
        
	hold off
end



