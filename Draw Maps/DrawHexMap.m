function HexMap = DrawHexMap(GridSize,HexSize,Labeling,Orientation,TogTicks)
% --------------------------------------------------------------------------------------------------
% DrawHexMap draws and generates labels for a hex map of a specified number of cells and
% and cell diameter
% --------------------------------------------------------------------------------------------------
% INPUTS:
%  GridSize				Specify size of grid
%  HexSize				Specify diameter of hexes
%  labeling				Label all regions
%  Orientation			Specify orients of cells: 0 = edges horizontal, 1 = edges vertical
%  TogTicks				Turn axis tick marks on or off
%
% OUTPUTS:
%  HexMap.Labels		Array of labels for each cell--'Labels' or 'IDs'
%		 .PatchIDs		Array of patch IDs for hexagonal cells
%		 .NumCells		Total number of cells in grid
%		 .AdjGph		Graph of grid area
%		 .AdjMat		Adjacency matrix for grid area
%		 .GridNbrs		Array of neighbor nodes for each cell
%		 .Handles		Array of polygon handles for cell
%		 .MaxXY			Max X,Y grid coordinates 
% ==================================================================================================


%% Check variables

% Check GridSize
if exist('GridSize','var')==0 || isempty(GridSize)
	% If GridSize is not specified, default 25
	GridSize = 25; 	
elseif mod(GridSize,2)== 0
	% If the specified grid size is even, make it odd so that meshgrid will generate properly
	GridSize = GridSize + 1; 
end
%
if GridSize > 50
    warning('Larger grids may take a long time to compute.')
elseif GridSize > 110
    error('Input Error: Grids much greater than 100x100 will cause Matlab instability.')
end

% If HexSize is not specified, default 1
if exist('HexSize','var')==0 || isempty(HexSize); HexSize = 1; end

% If Labeling is not specified, default to no labels
if exist('Labeling','var')==0 || isempty(Labeling); Labeling = ''; end

% If orientation is not specified, default to edges vertical
if exist('Orientation','var')==0 || isempty(Orientation); Orientation = 1; end

% If tickmarks are not specified, default to ticks off
if exist('TogTicks','var')==0 || isempty(TogTicks); TogTicks = 1; end


%% Generate hexagonal grid

Rad3Over2 = sqrt(3) / 2;

switch Orientation
	case  0
		[X0,Y0] = meshgrid(0:1:GridSize);
		n		= size(X0,1);
		s		= numel(X0);
		X		= Rad3Over2 * X0 * HexSize;
		Y		= (Y0 + repmat([0 0.5],[n,n/2])) * HexSize;
	case 1
		[Y0,X0] = meshgrid(0:1:GridSize);
		n		= size(Y0,1);
		s		= numel(Y0);
		X		= (X0 + repmat([0 0.5],[n,n/2]))' * HexSize;
		Y		= (Rad3Over2 * Y0)' * HexSize;
end
XY    = [reshape(X,1,s)' reshape(Y,1,s)'];
[v,p] = voronoin(XY);

% To store handles to the generated patches
pHandle = nan(s,1);  

% Draw patches and store patch handles
if TogTicks
    set(gca,'Xtick',[],'Ytick',[]);
end
for i = 1:s
    pHandle(i) = patch(v(p{i},1), v(p{i},2),NaN,'EdgeColor',[0.75 0.75 0.75]);
end


%% Generates Patch IDs and labels for all cells

lbl = 0;
for i = 1:s
	if numel(p{i})>=5 && all(p{i}~=1)
		lbl		   = lbl + 1;
		lblid(lbl) = i;
		
		HexMap.Labels{i}   = num2str(lbl);
		HexMap.PatchIDs{i} = num2str(i);
		
	else
		HexMap.Labels{i}   = '';
		HexMap.PatchIDs{i} = '';
	end
end

% Store labeling values
HexMap.Labels	= HexMap.Labels';
HexMap.PatchIDs = HexMap.PatchIDs';

% if Labeling is true then labels all cells	
if strcmp(Labeling,'IDs') || strcmp(Labeling,'Labels')
    if GridSize > 52
        warning('Grid is too large to completely label.')
	else
        % Ensure label are a legible size
        if GridSize < 22
            labelsize = 12;
        elseif GridSize < 42
            labelsize = 8;
        else
            labelsize = 4;
        end
        
        switch Labeling
            case 'Labels'
                tgrd = text(XY(:,1), XY(:,2), ...
					HexMap.Labels,...
					'HorizontalAlignment','center', ...
					'FontSize',labelsize, ...
					'Color',[0.75 0.75 0.75]);
            case 'IDs'
                tptc = text(XY(:,1), XY(:,2), ...
					HexMap.PatchIDs, ...
					'HorizontalAlignment','center', ...
					'FontSize',labelsize, ...
					'Color',[0.75 0.75 0.75]);
        end
    end
end


%% Find grid neighbors

nbr  = 0;						% initialize counter for grid neighbors array
nmat2  = 0;						% initialize
AdjMat = zeros(length(lblid));	% initialize adjacency matrix

NNs   = rangesearch(XY,XY,HexSize+0.1); % find all nearest neighbors

% Find neighbors within region bounds 
for grd = 1:s
    if ~isempty(HexMap.PatchIDs{grd})
        nbrmat                 = cell2mat(NNs(grd));
        nbrmat                 = intersect(lblid,nbrmat);
        nbrmat(nbrmat(:)==grd) = [];		% Remove self from neighbor list
        nbr                    = nbr + 1;
		Patches(nbr)		   = grd;		% Preserve patch ID
        Labels2(nbr)           = nbr;		% Preserve grid label
        GridNbrsP(nbr)         = {nbrmat};	% Store neighbors
		
		nbrL = zeros(1,length(nbrmat));
		for i = 1:length(nbrmat)
			nbrL(i) = find(lblid==nbrmat(i));
		end
		GridNbrsL(nbr)         = {nbrL};	% Store neighbors
		XY2(nbr,1)			   = XY(grd,1);	% Store Coordinate
		XY2(nbr,2)			   = XY(grd,2);	% Store Coordinate
        
        for ctr = 1:length(nbrmat)
            nmat       = nbrmat(ctr);
            nmat2(ctr) = find(lblid==nmat);
		end
		
        AdjMat(nbr,nmat2) = 1;
        nmat2 = [];
    end
end

AdjGph				= graph(AdjMat);				% Creat Graph
AdjGph.Edges.Weight = AdjGph.Edges.Weight*HexSize;	% Change arc weights
AdjMat				= adjacency(AdjGph);			% Create sparse adjacency matrix
GridNbrsP			= GridNbrsP';					% Flip neighbor list to vertical
ArcList				= AdjGph.Edges{:,1};			% Put Edges in an array (for convenience)


%% Make outputs

HexMap.NumCells  = lbl;
HexMap.AdjGph    = AdjGph;
HexMap.AdjMat	 = AdjMat;
HexMap.GridNbrsP = GridNbrsP;
HexMap.GridNbrsL = GridNbrsL;
HexMap.Patches   = Patches;
HexMap.Labels2   = Labels2;
HexMap.Handles   = pHandle;
HexMap.ArcList   = ArcList;
HexMap.XY		 = XY2;
HexMap.MaxXY	 = max(XY2);


% Resize axes
xyM = max(XY2);
xym = min(XY2);

switch  Orientation
	case 0
		axis equal, axis([-inf inf 0 GridSize*HexSize]), zoom on
	case 1
% 		axis equal, axis([0.5*HexSize GridSize*HexSize -inf inf ]), zoom on
		axis equal, axis([xym(1) xyM(1) xym(2) xyM(2)]), zoom on
end
end

%
%