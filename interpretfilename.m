function Fileinfo = interpretfilename(ExperimentFileName)
%Reads file name, extracts useful information and fills in text boxes on
%control panel.
%Example file name:
%ch4_air_1atm_phi120_ppcm938_sg075_s095_ss010_8k_512by512_jan30_2006andrew.csv

A=strread(char(ExperimentFileName),'%s','delimiter','_');

Fileinfo.fuel=char(A(1));
Fileinfo.oxygeninert=char(A(2));


Fileinfo.temperature = 298; %default
Fileinfo.ppcm = 944/10;  %default

for i=3:length(A)
    current=char(A(i));
    
    if(strcmp(current(length(current)),'k')) %Check for frame rate
        Fileinfo.framerate=1000*str2num(current(1:length(current)-1));
        %set(handles.framespersecbox,'String',num2str(1000*framerate));
    end
    
    if(strcmp(current(1),'T')) %Check for Temperature
        Fileinfo.temperature=str2num(current(2:length(current)));
    end
    
    if(length(current)>3)
        if(strcmp(current(length(current)-2:length(current)),'atm')) %Check for pressure
            Fileinfo.pressure=str2num(char(current(1:length(current)-3)));
        elseif(strcmp(current(1:3),'phi')) %Check for phi
            Fileinfo.phi=str2num(char(current(4:length(current))))/100;
        end
    end
    
    if(length(current)>4) %Check for ppcm
        if(strcmp(current(1:4),'ppcm'))
            Fileinfo.ppcm=str2num(current(5:length(current)))/10;
        elseif(strcmp(current(1:3),'fps'))
            Fileinfo.framerate=str2num(current(4:length(current)));
        elseif(strcmp(current(length(current)-3):current(length(current)),'ppcm'))
            Fileinfo.ppcm=str2num(current(1:length(current)-4))/10;
        end
    end

end


