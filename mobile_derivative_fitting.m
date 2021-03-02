filename = '8V\B_Tauc.txt';             %file to be read

headerslineIn = 2;                      %file formatting info
delimiterIn = '\t';                     %

ptsder = 3;                             %pts to use in the derivative calculation - don't change
ptsfit = 8;                             %pts to use in the horizontal fitting

eVmin = 3;                              %energy boundaries
eVmax = 4;                              %

%DATA ACQUISITION%
imported = importdata(filename,delimiterIn,headerslineIn);  %importing

rowscols = size(imported.data);         %initialising data:
eV    = zeros(rowscols(1),1);           % energy
cody  = zeros(rowscols(1),1);           % imaginary refractive index
eVder = zeros(rowscols(1),1);           % energy for the derivative plot
der   = zeros(rowscols(1),1);           % derivative

for i = (1:rowscols(1))                 %assigning data
    eV(i) = imported.data(i,1);         %
    eVder(i) = imported.data(i,1);      %
    cody(i) = imported.data(i,2);       %
end                                     %

%DERIVATIVE CALCULATION%
delta = floor(ptsder/2);                %pts below/above the "central one" for the calculation
i = rowscols(1)-delta;                  %we go right-to-left
while i > delta                         %derivative calculation
    der(i) = ( cody(i+delta) - cody(i-delta) ) / ( eV(i+delta) - eV(i-delta) );
    i = i-ptsder;                       %
end                                     %

i = rowscols(1);                        %removing unassigned elements of the arrays
while i > 0                             %
    if der(i) == 0                      %
        eVder(i)=[];                    %
        der(i)=[];                      %
    end                                 %
    i = i-1;                            %
end                                     %

%HORIZONTAL FIT CALCULATION%
rowscolsder = size(eVder);              %derivative dataset dimension

A = 1;                                  %finding energy boundaries
while eVder(A) < eVmin                  %
    A = 1+A;                            %
end                                     %
                                        %
B = rowscolsder(1)-1;                   %
while eVder(B) > eVmax                  %
    B = B-1;                            %
end                                     %

index = 1;                              %serial index of the fits
for i = (A:B)                           % number of fits to be performed
    for j = ((i+ptsfit-1):(B+1))        %
        index = 1+index;                %
    end                                 %
end                                     %

fitnumber = index-1;                    %actual number of performed fits
lower = zeros((index-1),1);             %lower and-
upper = zeros((index-1),1);             %-upper boundaries of the fit
avg   = zeros((index-1),1);             %the horizontal fit is actually the mean value
chi2  = zeros((index-1),1);             %reduced chi squared

index = 1;                              %horizontal fit and reduced-chi^2 calculation
for i = (A:B)                           %
    for j = ((i+ptsfit-1):(B+1))        %
        lower(index) = i;               %
        upper(index) = j;               %
        avg(index) = mean(der(i:j));    %
        chi2(index) = 0;                %
        for k = (i:j)                   %
            chi2(index) = ( ( (der(k)-avg(index)) * (der(k)-avg(index)) ) / avg(index) ) + chi2(index);
        end                             %
        chi2(index) = chi2(index) / (ptsfit-2);
        index = 1+index;                %
    end                                 %
end                                     %

%FINDING BEST FIT PARAMETERS%
i = fitnumber;
while i > 0                             %removing negative chi squared
    if chi2(i) < 0                      %
        chi2(i) = [];                   %
        avg(i) = [];                    %
        lower(i) = [];                  %
        upper(i) = [];                  %
    end                                 % 
    i = i - 1;                          %
end                                     %

[minchi2, minindex] = min(chi2);
fitavg = avg(minindex);
fitlower = eVder(lower(minindex));
fitupper = eVder(upper(minindex));

%PLOTTING%
plottitle = filename(4);              
if length(plottitle) > 1
        plottitle(2) = '†';
        plottitle(3) = '-';
end

plot(eVder,der,'x')
title(plottitle)                    
xlim([eVmin eVmax])                                                    
xlabel('eV')    
ylabel('Tauc/Cody Pseudo-Derivative')
line([eVmin, eVmax],[fitavg,fitavg])
txt = ['Boundaries: ' num2str(fitlower) '-' num2str(fitupper) ' eV; Chi^2 = ' num2str(minchi2) char(10) 'Pts used for derivative: ' num2str(ptsder) '; Pts used for fitting: ' num2str(ptsfit)];
text(3.35, 200, txt)