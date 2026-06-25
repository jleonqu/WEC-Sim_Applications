nS_T = 5;
nS_A = 7;

T_vec = linspace(4,8,nS_T);
A_vec = linspace(0.5,2.0,nS_A);

T_M = zeros(nS_T,nS_A);
A_M = zeros(nS_T,nS_A);
Kdamping_M = zeros(nS_T,nS_A);
Kpos_M = zeros(nS_T,nS_A);


for i = 1:1:nS_T
    for j = 1:1:nS_A
        WaveConditionNumber = (i-1)*nS_A+j;%WCN
        WCNString = num2str(WaveConditionNumber);
        %fileName = strcat('TimeDomainResultsV1\TimeDomainSim_',WCNString,'.mat');
        fileName = strcat('Kpos_Kdamping_',WCNString,'.mat');
        ResultsData = load(fileName);

        T_M(i,j) = ResultsData.Data.T;
        A_M(i,j) = ResultsData.Data.A;
        Kdamping_M(i,j) = ResultsData.Data.Kdamping;
        Kpos_M(i,j) = ResultsData.Data.Kpos;
    end
end

figure(1)
h1 = heatmap(T_vec, A_vec,T_M');
%h1.Colormap = parula;
h1.Title = 'Wave Period [s]';
h1.YLabel = 'A (m)';
h1.XLabel = 'T (s)';

figure(2)
h2 = heatmap(T_vec, A_vec,A_M');
%h1.Colormap = parula;
h2.Title = 'Wave Amplitude [m]';
h2.YLabel = 'A (m)';
h2.XLabel = 'T (s)';

figure(8)
h8 = heatmap(T_vec, A_vec,Kdamping_M');
h8.Title = 'Kdamping';
h8.YLabel = 'A (m)';
h8.XLabel = 'T (s)';

figure(9)
h9 = heatmap(T_vec, A_vec,Kpos_M');
h9.Title = 'Kpos';
h9.YLabel = 'A (m)';
h9.XLabel = 'T (s)';
%h9.XLabel = 'T (s)';
%h9.XLabel = 'T (s)';
% 
%f = figure(9); exportgraphics(f,'KposWaveBot2X.png','Resolution',300)​