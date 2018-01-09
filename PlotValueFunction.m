function PlotValueFunction(FOB,PARAMETERS)

    X0 = roots(PARAMETERS.VFCS);
    X1 = linspace(min(X0),FOB.CAPACITY);
    X2 = X1/PARAMETERS.SCALING_FACTOR;
	
    Y1 = PARAMETERS.VFCS(1)			+ PARAMETERS.VFCS(2)*X1			+ PARAMETERS.VFCS(3)*(X1.^2);
    Y2 = PARAMETERS.VFCS_SCALED(1)	+ PARAMETERS.VFCS_SCALED(2)*X2	+ PARAMETERS.VFCS_SCALED(3)*(X2.^2);

	figure; hold on
	subplot(1,2,1)
		plot(X1,Y1);
		axis square
		axis([-inf max(X1) min(Y1) max(Y1)*1.1]);
	subplot(1,2,2)
		plot(X2,Y2,'r');
		axis square
		axis([-inf max(X2) min(Y2) max(Y2)*1.1]);
	hold off
end