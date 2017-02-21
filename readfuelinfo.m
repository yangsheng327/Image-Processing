function [excelfilename, stanwinfilename]=readfuelinfo(Fileinfo)

%cin=fopen('N:\Experiment\Experiment Analysis\Stanwin Files\fuels.csv','r');

cin=fopen('../Stanwin Files/fuels.csv','r');

%cin=fopen('fuels.csv','r');
i=1;
line=fgetl(cin); %Gets rid of headers

while(~feof(cin))
    line=fgetl(cin);
    %Fuel info is stored as: fuel, oxidizer/inert, temperature,
    %stanwinfile, resultsfile
    [fuelinfo(i,1) fuelinfo(i,2) temp(i) fuelinfo(i,3) fuelinfo(i,4)]=strread(line,'%s %s %u %s %s','delimiter',',');
    i=i+1;
end


temperature=Fileinfo.temperature;
inert=Fileinfo.oxygeninert;
fuel=Fileinfo.fuel;


% strcmpi(fuel,char(fuelinfo(117,1)))
% strcmpi(inert,char(fuelinfo(117,2)))
% temperature
% temp(117)

for j=1:i-1
    if( strcmpi(fuel,char(fuelinfo(j,1))) && strcmpi(inert,char(fuelinfo(j,2))) && temperature==temp(j) )
        excelfilename=char(fuelinfo(j,4));
        stanwinfilename=char(fuelinfo(j,3));
    end
end

fclose(cin);

return