clc
close all
clear all

%%
filename = 'Ub=Uc_PeakFrequencies_N50.mat';
load(filename);

%%
Excludes = cell(1,3);
Excludes{1} = [1:5];
Excludes{2} = [6:11];
Excludes{3} = [1:5];

Data = cell(1,3);
Data{1} = FreqLow;
Data{2} = FreqMedium;
Data{3} = FreqHigh;

Names = cell(1,3);
Names{1} = 'FreqLow';
Names{2} = 'FreqMedium';
Names{3} = 'FreqHigh';
%%

for i = 1:3
[xData, yData] = prepareCurveData(U, Data{i} );

ft = fittype( 'poly1' );
ex = excludedata( xData, yData, 'Indices', Excludes{i} );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = ex;

[fitresult, gof] = fit( xData, yData, ft, opts );
linfit_str = [Names{i},' = ' ,num2str(fitresult.p1),'*U + ', num2str(fitresult.p2)];

figure( 'Name', Names{i} );
h = plot( fitresult, xData, yData, ex );
legend( h, Names{i}, ['Excluded ',Names{i}], linfit_str, 'Location', 'NorthEast' );

xlabel( 'U' );
ylabel( Names{i} );
end
