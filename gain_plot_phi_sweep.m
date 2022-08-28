clc;
clear all
close all;



%% Searching for probe order

promt='Introduce the number of angles measured: '; %Ask the user for the number of angles measured
degrees=zeros(1,input(promt)); %Create a vector with the obtained length

fid=fopen('modelo4_phi_sweep_dat.txt','r'); %Open the file to be read
aux2=0;
for i=1:length(degrees)
    linenum=2+1004*aux2; %Jump all the numerical data
    var(i)=textscan(fid,'#"Frequency / GHz"	"E_Field (Farfield) (Spherical) (90 %d 1000)(Abs) [1] [Magnitude]"',1,'Delimiter','\n','HeaderLines',linenum-1); %Find the order in the data headers
    aux2=aux2+1; %Update the auxiliar variable
    fseek(fid,0,"bof"); %Move to the specified position in the file
end
fclose(fid); %Close the file
var=cell2mat(var); %Convert a cell-type variable in a int32

degrees=var; %The vector contains the order given by CST
degrees1=-(length(degrees)-1)/2*10:10:(length(degrees)-1)/2*10; %The vector contains the order needed for the plot



%% Importing data

pos=1;
vmax=1;
vmin=200;
for i=0:length(degrees)-1
    data(:,1)=dlmread('modelo4_phi_sweep_dat.txt','',[3+1004*i 0 1003+1004*i 0]); %Read the frequency samples
    data(:,2)=dlmread('modelo4_phi_sweep_dat.txt','',[3+1004*i 1 1003+1004*i 1]); %Read the electric field values
    nm=length(data(:,2)); %Number of samples index
    aux=1;
    for k=1:1:nm
        f(pos)=data(k,1); %Create a vector with the frequency samples
        v(pos)=data(k,2); %create a vector with the electric field values
        pos=pos+1; %Update the position variable
    end
end



%% Ordering the data in the needed order


pos=0;
newv=zeros(1,length(v)); %Create the vector for the reordered data
for i=1:length(degrees)
    ind=Encontrar(degrees1,degrees(i)); %Use the created function "Encontrar" to obtain the index where the following sample chunk belongs
    for j=1:nm
             newv(pos+j)=v(((ind-1)*nm)+j); %Reorder the data in the new vector
    end
    pos=pos+nm; %Jump all the chunk
end

v=newv; %Return the data to the vector v
nf=length(v); %Total number of samples index (in all the angles)



%% Importing the reflection coefficient


s=importdata('modelo4_s11.txt'); %Import the file into a struct
for i=1:length(s.data(:,1))
    s11(i)=s.data(i,2); %Create a vector with the values of the reflection coefficient in dB
    s11(i)=10^(s11(i)/20); %Convert the reflection coefficient into linear values
end



%% Converting the electric field into gain


pos=1;
for i=1:1:nf
    if f(i)==f(1)
        pos=1; %Restart the position for s11 parameters every angle
    end
    v(i)=10^(v(i)/20); %Turn the electric field in dB(V/m) into linear values
    v(i)=(abs(v(i))^2); %Get the square of the electric field
    gain(i)=v(i)/(1-s11(pos)^2); %Distinguish the operations applied to the absolut gain and the realized gain
    rlzdgain(i)=v(i)/30; %Conversion to realized gain
    gain(i)=gain(i)/30; %Conversion to absolut gain
    gain(i)=10*log10(gain(i)); %Turn it into dB
    rlzdgain(i)=10*log10(rlzdgain(i)); %Turn it into dB
    pos=pos+1; %Update the positioning variable
end

pos=1;
agmax=1;
agmin=200;
rgmax=1;
rgmin=200;
for i=1:nf
    if rlzdgain(i)>=rgmax
        rgmax=rlzdgain(i); %Obtain the realized gain value
    end
    if rlzdgain(i)<=rgmin
         rgmin=rlzdgain(i); %Obtain the minimum realized gain value
    end
    if i==nm
        fmax=f(i); %Obtain the maximum frequency value
    end
end
fmin=f(1); %Obtain the minimum frequency value

for i=1:nf
    if gain(i)>=agmax
        agmax=gain(i); %Obtain the maximum gain value
    end
    if gain(i)<=agmin
         agmin=gain(i); %Obtain the minimum gain value
    end
end

freq=zeros(1,nm); %Create a number of samples-length vector
i=1;
while (i<(length(f)/length(degrees))+1)
        freq(i)=f(i); %Obtain the frequency values just once, not as many times as angles we have
        i=i+1; %Update the loop variable
end



%% Creating the matrix with the gain value for each angle (row) and frequency point (column)


C1=zeros(length(degrees), nf/length(degrees)); %Create the absolut gain matrix that will be plotted
C2=zeros(length(degrees), nf/length(degrees)); %Create the realized gain matrix that will be plotted
aux=0;
j=1;
for k=1:1:length(degrees)
    for i=aux+1:1:nf
        if f(i)==fmax
            j=1;
            aux=i;
            break
        end
        C1(k,j)=gain(i); %Fill the absolut gain matrix
        C2(k,j)=rlzdgain(i); %Fill the realized gain matrix
        j=j+1;
    end
end

%% Relaxing the resolution (in order to achieve more visibility of gain values in the plot)

M1=zeros(length(degrees), 334); %Create a absolut gain matrix to get a lower resolution
M2=zeros(length(degrees), 334); %Create a realized gain matrix to get a lower resolution
aux=1;
for i=1:3:length(freq) %Sampling the frequency vector (334 samples instead of 1001)
    frequency(aux)=freq(i); %Fill the vector
    aux=aux+1; %Update the positioning variable
end
aux=1;
for k=1:length(degrees) %Sampling each row of C such that we obtain a new vector M[13x334]
    aux=1;
    for i=1:3:length(freq)-3
        M1(k,aux)=C1(k,i); %Fill the absolut gain matrix
        M2(k,aux)=C2(k,i); %Fill the realized gain matrix
        aux=aux+1; %Update the positioning variable
    end
end

%% Plotting



figure;
surf(frequency,degrees1,M1);
title('Absolut Gain in E-plane');
xlabel('Frequency');
ylabel('Phi');
zlabel('Absolut Gain');
axis([fmin fmax-0.036 degrees1(1) degrees1(length(degrees1)) agmin agmax]); %Variable depending on the antena operating frequency
c=colorbar;

figure;

surf(frequency,degrees1,M2); %Plot
title('Realized Gain in E-plane');
xlabel('Frequency');
ylabel('Phi');
zlabel('Realized Gain');
axis([fmin fmax-0.036 degrees1(1) degrees1(length(degrees1)) rgmin rgmax]); %Variable depending on the antena operating frequency
c=colorbar;

%PLOT ADVISE: View -> Property inspector -> Colour and transparency maps ->
%colormap=Turbo.
%PLOT ADVISE: View -> Property inspector -> Colour and transparency maps ->
%CLin=[lb,ub] where lb and ub are the lower bound and the upper bound of
%our gain values.


%% Functions



%This function returns the position of the element Y in the vector x
function[ind]=Encontrar(x,Y)
    n=length(x);
    for i=1:n
        if x(i)==Y
            ind=i;
            break;
        end
    end
end