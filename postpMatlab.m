clc; clear;
root = pwd;
path = root + "/run/data/";


%% Store de 1d BL data
filesData = dir(path+"*1d_bl*.out");

for i=1:length(filesData)
    filename = filesData(i).name;
    data = readtable(path+filename,"FileType","text");
    datTemp = table2array(data);
    wss(:,i) = datTemp(:,1);
    redelta(:,i) = datTemp(:,2);
    retheta(:,i) = datTemp(:,3);
end


n = linspace(1,128,128);
figure()
plot(n,redelta(:,100)./retheta(:,100))
%%
filesData = dir(path+"*uPlus*.out");
for i=1:length(filesData)
    filename = filesData(i).name;
    data =(readtable(path+filename,"FileType","text"));

    cellMask = varfun(@iscell, data, 'OutputFormat', 'uniform');
    for col = find(cellMask)
        stringIdx = cellfun(@ischar, data{:,col}) | cellfun(@isstring, data{:,col});
        numericCells = cellfun(@str2double, data{:,col}, 'UniformOutput', false);
        data.(data.Properties.VariableNames{col}) = cell2mat(numericCells);
    end

    uPlus(:,:,i) = table2array(data);

end

filesData = dir(path+"*yPlus*.out");
for i=1:length(filesData)
    filename = filesData(i).name;
    data =(readtable(path+filename,"FileType","text"));

    cellMask = varfun(@iscell, data, 'OutputFormat', 'uniform');
    for col = find(cellMask)
        stringIdx = cellfun(@ischar, data{:,col}) | cellfun(@isstring, data{:,col});
        numericCells = cellfun(@str2double, data{:,col}, 'UniformOutput', false);
        data.(data.Properties.VariableNames{col}) = cell2mat(numericCells);
    end

    yPlus(:,:,i) = table2array(data);
    % data =flip(data,1)
end
figure()
semilogx(yPlus(:,100,3),uPlus(:,100,3))