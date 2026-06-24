function OptimizationResultsLoadData = OptimizationResultsLoad(T,A)

TMin = 4;
TMax = 8;
Tstep = 1;

AMin = 0.5;
AMax = 2.0;
Astep = 0.25;

[TMatrix,AMatrix] = meshgrid(TMin:Tstep:TMax,AMin:Astep:AMax);
KdampingMatrix = load('KposKdamping\KdampingMatrix.mat');
KposMatrix = load('KposKdamping\KposMatrix.mat');

Kpos = interp2(TMatrix,AMatrix,KposMatrix.Kpos_M,T,A);
Kdamping = interp2(TMatrix,AMatrix,KdampingMatrix.Kdamping_M,T,A);

Data = [];
Data.A = A;
Data.T = T;

if isnan(Kdamping)
    Kdamping = -500; % set to min damping if T or A are not in range
end

if isnan(Kpos)
    Kpos = 0;
end

Data.Kdamping = Kdamping;
Data.Kpos = Kpos;


OptimizationResultsLoadData = Data;%load the file



end