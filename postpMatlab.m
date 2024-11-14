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


n = linspace(1,100,100);
hold on 
for i=1:100
    t = i*ones(100);
    plot3(n,t,wss(1:100,i))
end 
hold off
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
    % data =flip(data,1)
    contourf(ddata)
    colorbar
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
    contourf(ddata)
    colorbar
end
