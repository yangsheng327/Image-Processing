function [rhoratio, rhou]=readstanwin(handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Opens appropriate stanwin file for specified fuel and calculates density
%ratio for specified equivalence ratio and pressure

equivalence=handles.phi;
pressure=handles.pressure;
filename=handles.stanwinfilename;
stanin=fopen([pwd,'\..',char(filename)],'r');
i=1;
line=fgetl(stanin);
titles(i,:)=strread(line,'%s',-1,'delimiter',',');

%Figure out which column is used for densities and given pressures and phis
i=1;
while ~strcmp('Phi',titles(i))&& i<=length(titles)
    i=i+1;
end
Phiindex=i;

while ~strcmp('P[atm]',titles(i))&& i<=length(titles)
    i=i+1;
end
Pressureindex=i;

while ~strcmp('V1(cm3/gm)',titles(i))&& i<=length(titles)
    i=i+1;
end
V1index=i;

while ~strcmp('V2(cm3/gm)',titles(i))&& i<=length(titles)
    i=i+1;
end
V2index=i;


philist = [];
preslist = [];
V1list = [];
V2list = [];

j=1;
%save for relevant data
while(~feof(stanin))
    line=fgetl(stanin);
    values=strread(line,'%f',-1,'delimiter',',');
	philist(j) = values(Phiindex);
	preslist(j) = values(Pressureindex);
	V1list(j) = values(V1index);
	V2list(j) = values(V2index);
	j = j+1;

end
fclose(stanin);


%Search for relevant data
flag=0;
    for i=1:length(philist)
        if(philist(i)==equivalence && preslist(i)==pressure)
        rhoratio=V2list(i)/V1list(i);
        rhou=1/V1list(i);
		flag=1;
        break;
        end
    end
	
	if flag==0
	
	           phiclose = [];
			   V1close = [];
			   V2close = [];
	           j = 1;
	           for i=1:length(philist)
	               if(preslist(i)==pressure)
			          phiclose(j)=philist(i);
					  V1close(j)=V1list(i);
					  V2close(j)=V2list(i);
				      j = j+1;
	               end
			   end

			[~, index1] = min(abs(equivalence-phiclose));  
			if phiclose(index1) > equivalence 
			    index2 = index1-1;
			else
				index2 = index1+1;
			end
			%phiclose(index1)
			%phiclose(index2)
			%V1close(index2)
			%V1close(index1)
			%V2close(index2)
			%V2close(index1)
			V1want = V1close(index1)+(equivalence-phiclose(index1))*(V1close(index2)-V1close(index1))/(phiclose(index2)-phiclose(index1));
			V2want = V2close(index1)+(equivalence-phiclose(index1))*(V2close(index2)-V2close(index1))/(phiclose(index2)-phiclose(index1));
			%V2want
			%V1want
	        rhoratio=V2want/V1want;
            rhou=1/V1want;	
	end


return