% This function generates the input data for the performance tests of the
% motors. The output of this function is a .mat file.
% The inputs of this function are:
%%%% wShaft: it's a vector with the values of shaft speed that will be
%    tested: wShaft = [w1 w2 w3 ... wn]
%    wShaft -----> rpm
%%%% torqueShaft: it's a vector with the values of torque that will be
%    tested: torqueShaft = [torque1 torque2 torque3 ... torquen]
%    torqueShaft -----> Nm
%%%% freq: it's the data frequency in Hz
%    freq -----> Hz
%%%% tTest: time for each test [wi,torquei]
%    tTest -----> seconds

function experimentDataOut = experimentDataFunc(wShaft,torqueShaft,freq,tTest)

count = 0;
dataOutput = [];
for i=1:1:length(wShaft)
    wShaft_i = wShaft(i);
    for j=1:1:length(torqueShaft)
        count = length(torqueShaft)*(i-1)+j;
        torqueShaft_j = torqueShaft(j);
        dataPreliminary = ones(freq*tTest,3);
        dataPreliminary(:,1) = dataPreliminary(:,1)*wShaft_i;
        dataPreliminary(:,2) = dataPreliminary(:,2)*torqueShaft_j;
        dataPreliminary(:,3) = dataPreliminary(:,3)*count;
        dataOutput = [dataOutput;dataPreliminary];
    end
end

experimentDataOut = [];
experimentDataOut.wShaft_rpm = dataOutput(:,1);
experimentDataOut.torqueShaft_Nm = dataOutput(:,2);
experimentDataOut.count = dataOutput(:,3);
experimentDataOut.freq_Hz = freq;
experimentDataOut.tTest_s = tTest;

end