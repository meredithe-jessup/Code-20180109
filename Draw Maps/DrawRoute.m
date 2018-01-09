function Route = DrawRoute(DepCoord,FOBCoords,FOBNAMES,ORDER)

for s = 1:length(FOBNAMES)
    StrLn(s) = length(num2str(FOBNAMES(s)));
    NAMES{s} = FOBNAMES(s);
end

DepX = DepCoord(1);
DepY = DepCoord(2);

FOBX = FOBCoords(:,1);
FOBY = FOBCoords(:,2);

Route.pth  = plot([DepX;FOBX(ORDER);DepX],[DepY;FOBY(ORDER);DepY],'k');
Route.fobs = scatter(FOBX,FOBY,100,'c','filled','MarkerEdgeColor','k');
Route.txt  = text(FOBX+1.5,FOBY,NAMES,'Color','b','FontName','FixedWidth');

end