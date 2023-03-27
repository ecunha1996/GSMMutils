function data = readExpMatrix(filename)

[~,sheet_name]=xlsfinfo(filename);

data = {30};

for k=1:length(sheet_name)
    % Import the data
    data{k} = readtable(filename, 'Sheet', sheet_name{k}, 'VariableNamingRule', 'preserve');
end
A = table2array(data{32}(:,{'Trial'}));
data{32}.Properties.RowNames = A;
end