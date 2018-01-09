clear
close all

saves = false; % 'BaseGrid', 'Blue', 'Red', true = 'All'

if ~saves
	warning('Figures will not be saved')
end

Inputs.datafile    = 1;
Inputs.sendtext	   = false;
Inputs.GridSize	   = 25;
Inputs.HexSize	   = 4;
Inputs.Labeling    = 'Labels'; % Array of labels for each cell--'Labels' or 'IDs'
Inputs.Orientation = 1;
Inputs.TogTicks    = 1;
Inputs.redlvl      = 0.1;
Inputs.NumMaps     = 4;
Inputs.Sync        = true;


%%
switch Inputs.datafile
	case 1
		load 'R1_coord.mat';
		Inputs.data = R1_coord;
	case 2
		load 'C1_coord.mat';
		Inputs.data = C1_coord;
	case 3
		load 'RC1_coord.mat';
		Inputs.data = RC1_coord;
end

Depot       = Inputs.data(end,:);
FOBs		= Inputs.data(1:(end-1),:);



%% Set random number stream
% If Labeling is not specified, default to no labels
if exist('Inputs.Sync','var') == 0 || isempty(Inputs.Sync); Inputs.Sync = true; end

if Inputs.Sync
	Inputs.RNGState = rng('default'); % Sync random streams
end



%% Create Maps
for mp = 0:Inputs.NumMaps
    
	% Set figure layout
    MapName    = ['Map ',num2str(mp)];
	figs(mp+1) = figure('Name',MapName,'NumberTitle','off','units','inches','position',[0,0,8,7.5]); hold on
    title(['t = ',num2str(mp)])
    
    
    % Set and save grid layout
	MapGrid = DrawHexMap(Inputs.GridSize,Inputs.HexSize,Inputs.Labeling,[],Inputs.TogTicks);
	if mp == 0 && strcmp(saves,'All') || strcmp(saves,'BaseGrid')
		saveas(gcf, strcat('Grid',num2str(Inputs.GridSize)),'fig');
		saveas(gcf, strcat('Grid',num2str(Inputs.GridSize)),'png');
    end
	
    
	% Plot blue forces (depot and FOBs)
	MapBlue = PlotBlue(Inputs.GridSize,MapGrid,Depot,FOBs);
	if mp == 0 && strcmp(saves,'All') || strcmp(saves,'Blue')
		MapNameFOB = strcat('GridFOBs',num2str(Inputs.GridSize));
		saveas(gcf,MapNameFOB,'fig');
		saveas(gcf,MapNameFOB,'png');
	end
	hold off
	
    
	if mp == 1      % Plot initial threats
		MapRed  = PlotRed(MapGrid,Inputs.redlvl,MapBlue.DepLabel,Inputs.Sync);
		if  strcmp(saves,'All') || strcmp(saves,'Red')
			HexMapName = NameHexMap(Inputs.GridSize,36,50,Inputs.redlvl,1,0,0,0,0);
			saveas(gcf, strcat(HexMapName,'_',num2str(mp)),'png');
			saveas(gcf, HexMapName,'fig');
		end
		hold off
        
	elseif mp > 1   % Create Alternate Map(s)
		RedShiftMap = PlotRedShift(   ...
                    MapGrid.NumCells, ...
                    MapRed.RedLabels, ...
                    MapGrid.Patches,  ...
                    MapGrid.XY,       ...
                    MapGrid.Handles,  ...
                    MapGrid.GridNbrsP,...
                    MapGrid.GridNbrsL,...
                    0.9, 0.1);
		if strcmp(saves,'All') || strcmp(saves,'Red')
			HexMapName = NameHexMap(Inputs.GridSize,36,50,Inputs.redlvl,1,0,0,0,0);
			saveas(gcf, strcat(HexMapName,'_',num2str(mp)),'png');
			saveas(gcf, HexMapName,'fig');
		end
		hold off
	end
end


%% Send text message when program completes

if Inputs.sendtext
	Send_Text;
	openfig('figure1.fig');
end