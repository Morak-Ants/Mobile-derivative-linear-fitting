filename = '8V\D_Tauc.txt';             %file to be read

headerslineIn = 2;                      %file formatting info
delimiterIn = '\t';                     %

ptsder = 3;                             %pts to use in the derivative calculation

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
B = rowscolsder(1);                     %

fitpts = B - A + 1;                     %max number of points in fits
lower = zeros(fitpts,1);                %lower and-
upper = zeros(fitpts,1);                %-upper boundaries of the fit
avg   = zeros(fitpts,1);                %the horizontal fit is actually the mean value
chi2  = zeros(fitpts,1);                %reduced chi squared

i = fitpts;                             %horizontal fit and reduced-chi^2 calculation
while i > delta                         %
    lower(i) = A-i+fitpts;              %
    upper(i) = B;                       %    
    avg(i) = mean(der(lower(i):B));     %
    chi2(i) = 0;                        %
    for k = (lower(i):B)                %
        chi2(i) = ( ( (der(k)-avg(i)) * (der(k)-avg(i)) ) / avg(i) ) + chi2(i);
    end                                 %
    chi2(i) = chi2(i) / (i-2);          %
    i = i-1;                            %
end                                     %

%FINDING BEST FIT PARAMETERS%
i = fitpts;                             % the search starts from the highest # of pts for the fit
ptsfit = 0;                             % initialising the best number of pts for fitting
while i > 0                             %searching for the first chi2<1
    if chi2(i) < 1                      %
        ptsfit = i;                     %
        break                           %
    end                                 % 
    if  chi2(i) < 2 && chi2(i-1) > chi2(i) && chi2(i-2) > chi2(i-1)
        ptsfit = i;                     %searching for minima @ chi2<2
        break                           %
    end                                 %
    i = i - 1;                          %
end                                     %

if ptsfit < ptsder                      %if no chi2<1 or minimum is achieved, we choose the minimum chi2
    [~, ptsfit] = min(chi2(ptsder:fitpts));
    ptsfit = ptsfit + ptsder - 1;       %
end                                     %

minchi2 = chi2(ptsfit);
fitavg = avg(ptsfit);
fitlower = eVder(lower(ptsfit));
fitupper = eVder(upper(ptsfit));

%PLOTTING%
plottitle = filename(4);              

plot(eVder,der,'x')
title(plottitle)                    
xlim([eVmin eVmax])                                                    
xlabel('eV')    
ylabel('Tauc/Cody Pseudo-Derivative')
line([eVmin, eVmax],[fitavg,fitavg])
txt = ['Boundaries: ' num2str(fitlower) '-' num2str(fitupper) ' eV; Chi^2 = ' num2str(minchi2) char(10) 'Pts used for derivative: ' num2str(ptsder) '; Pts used for fitting: ' num2str(ptsfit)];
text(3.35, 750, txt)