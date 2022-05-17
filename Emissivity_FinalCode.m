%% Notes on data
%%%%% Pressure and temperature pairing (for Earth)
%%%%% From https://www.pdas.com/index.html 
% 1/e atm (37275 Pa) --> 239.5 K ~ 240K
% 1/e^2 atm (13713 Pa) --> 216.6 K ~ 217K
% 0.5 atm (50662 Pa) --> 252.4 K ~ 252K
% 1 atm (101325 Pa) --> 288.1K ~ 288K
%                    --> 1600K if simulating runaway greehouse

% Wavenumber range : 500 cm^-1 to 11501 cm^-1
% 500 half-widths
% cutoff on intensity : 1.0*10^(-28)

%% Constants
g       = 9.8               ; % Gravitational constant
h       = 6.626*10^(-34)    ; % Planck constant [J/Hz]
c       = 299792458         ; % Speed of light [m/s]
kb      = 1.38*10^(-23)     ; % BoltzmanSen constant [m^2*kg*s^-2*K^-1]
sigma   = 3.67*10^(-8)      ; % Stefan-Boltzman constant

%% User inputs 
%%% User input values: Temperature of the planet, Temperature higher
%%% altitude (if doing a higher altitude calculation), Pressure, Concentration of gases
%%% in atmosphere (mixing ratio). If using water vapor, enter the precipitable water vapor
%%% content in m.

Temperature_planet = 288; % This is the planet's SURFACE temperature
Pressure_atm = 1; % in atm
Pressure = Pressure_atm * 101325; %in Pa
Conc_1 = 0.216; % Here using Conc_1 as water vapor so it is a precipitable water vapor content rather and mixing ratio
Conc_2 = 0.00000182;
Conc_3 = 0.0004;

%% upload data
%%% SPECTRA downloads are in the form of textfiles which have a column of
%%% wavelengths and a column of cross-section data. The wavelengths should be
%%% equivalently spaced for each gas (if the selections for the
%%% gas's cross-sections data were the same) so you can use 1 gas's
%%% wavelength data for everything. 

% Define path directory for data
path = '/Users/isabellesangha/Downloads/Honours_research_project/';

% Specify folder for data
folder = "Data";

% Specify which dataset is which gas
Gas_1 = 'H2O';
Gas_2 = 'CH4';
Gas_3 = 'CO2';

% Specify Pressure and temperature (as strings). Must match what you entered above for atmospheric
% pressure.

T = "288";
P = "1";

% fill in the file name of the textfiles from SPECTRA
Filename_1 = sprintf("%s%s/%s_%sK_%satm.txt",path, folder, Gas_1,T,P);
Filename_2 = sprintf("%s%s/%s_%sK_%satm.txt",path,folder, Gas_2,T,P);
Filename_3 = sprintf("%s%s/%s_%sK_%satm.txt",path, folder, Gas_3,T,P);

%Read the textfile
FID_1 = fopen(Filename_1);
FID_2 = fopen(Filename_2);
FID_3 = fopen(Filename_3);

%%% The first few lines of the textfile from SPECTRA contain information about
%%% the datatype, the wavenumber interval, the resolution, etc. Using '%f%f'
%%% skips the first few lines of text directly pulls the data for wavenumber,
%%% cross-section, and absorption. They are ordered in a vector such that it
%%% goes: wavenumber, cross-section value, absorption,... 

% wavenumber data is in cm^{-1}, wavelength data is in micrometers and xsec
% data is in cm^2/molecule

% for loop runs through the first 6 lines of the data that give information
% about data format 
for i = 1:6
    text = fgetl(FID_1);
end
data_1 = fscanf(FID_1, '%f%f');
    wavenumber = data_1(1:3:end);
    wavelength = 10000000./wavenumber*0.001; %Converting from wavenumber in cm^(-1) to wavelength in micrometer
    xsec_1 = data_1(2:3:end); 
    
for i = 1:6
    text = fgetl(FID_2);
end
data_2 = fscanf(FID_2, '%f%f');
    xsec_2 = data_2(2:3:end);
    
for i = 1:6
    text = fgetl(FID_3);
end
data_3 = fscanf(FID_3, '%f%f');
    xsec_3 = data_3(2:3:end);

%% Calculating absorption coef, optical depth, transmissivity
%%% From the cross-section data we want to determine the absorption
%%% coefficient (kappa). This is done using the relationship between the
%%% cross-section value and the molecular weight of the gas

% Molar Mass of each gas
MolarMass_1 = 18.015;
MolarMass_2 = 44.01;
MolarMass_3 = 16.04;

% Molecular weight of each molecule in kg/molecule
mu_1 = MolarMass_1/1000/(6.022*10^(23));
mu_2 = MolarMass_2/1000/(6.022*10^(23));
mu_3 = MolarMass_3/1000/(6.022*10^(23));

% Absorption coefficent in m^2/kg
k_1 = xsec_1*0.0001/mu_1;
k_2 = xsec_2*0.0001/mu_2;
k_3 = xsec_3*0.0001/mu_3;

%%% Now calculate the weighting coefficient of each gas to be used when
%%% determine the average absorption. Slightly different method if water.
%%% Here using gas_1 as water vapor. 

density_H2O = 997; % Water vapor density

WeightingCoef_1 = Conc_1*density_H2O; % Using precipitable water vapor to get water vapor weighting coefficient
WeightingCoef_2 = Pressure*Conc_2/g; % Weighting coefficient calculation for well-mixed gases
WeightingCoef_3 = Pressure*Conc_3/g;

%%% Using the values calculated above we can now calculate the optical depth,
%%% transmittance, and absorbance for each wavelength for each gas

% Optical depth
Optical_depth_1 = k_1*WeightingCoef_1;
Optical_depth_2 = k_2*WeightingCoef_1;
Optical_depth_3 = k_3*WeightingCoef_1;

% Transmittance
Transmittance_1 = exp(Optical_depth_1*(-1)); 
Transmittance_2 = exp(Optical_depth_2*(-1));
Transmittance_3 = exp(Optical_depth_3*(-1));

% Absorbance
Absorbance_1 = 1 - Transmittance_1;
Absorbance_2 = 1 - Transmittance_2;
Absorbance_3 = 1 - Transmittance_3;

%%% To calculate the total absorption coefficient (for an atmosphere that is
%%% a mixture of the gases) we must sum the optical depths for each gas
%%% (these have already been appropriately weighted according to the
%%% concentration of the gas)

Optical_depth_tot = Optical_depth_1 + Optical_depth_2 + Optical_depth_3;
k_tot             = Optical_depth_tot/(Pressure/g);
Transmittance_tot = exp(Optical_depth_tot*(-1));
Absorbance_tot    = 1 - Transmittance_tot; 

%% Figure 1 : plot of absorption cross-sections (cm^2/molecule) for each gas
fig1 = figure(1);
t = tiledlayout(3,1);
t.YLabel.String = 'Absorption cross-section (cm^2/molecule)';

%top 
nexttile 
semilogx(wavelength, xsec_1)
title(sprintf("%s Absorption cross-section",Gas_1));

%middle (1) 
nexttile 
semilogx(wavelength, xsec_2)
title(sprintf("%s Absorption cross-section",Gas_2));

%middle (2) 
nexttile 
semilogx(wavelength, xsec_3)
title(sprintf("%s Absorption cross-section",Gas_3));
xlabel('Wavelength(\mum)');
%% Earth's emission spectrum
% To calculate the radiation emitted by Earth ( in [W*m^-2*sr^-1*Hz^-1]) use Planck's law

Spectral_radiance = (2*h*c^2)./(wavelength.*10^(-6)).^5 .* (1./(exp(h*c./((wavelength.*10^(-6)).* (kb*Temperature_planet)))-1)); 

%% Calculating mean absorbance (emissivity)
% Calculating wavelength interval in meters
dlambda = wavelength.^2;

%%%%%%% Mean Transmittance %%%%%%
%%% Using a weighted mean to calculated the transmittance averaged over the
%%% longwave wavelengths. The transmittance is weighted by the Planck
%%% function 


% Initialize the numerator and denominator for both the average weighted by
% the Planck function

Num = 0;
Den = 0;


%%% for loop runs through an iterative sum that weights each tranmittance
%%% value at each wavelength by the corresponding value given by plancks law
%%% or its derivativec at that wavelength. The numerator is the area under
%%% the weighting function (either planck or partial) and is found by summing
%%% the value at each wavelength

Start = 1; 
End = length(wavelength);

for i = Start : End 
    
    this_num = Transmittance_tot(i)*Spectral_radiance(i)*dlambda(i);
    this_den = Spectral_radiance(i)*dlambda(i);
    
    Num = Num + this_num; 
    Den = Den + this_den;

end

Mean_transmittance = Num/Den;

% Using the mean transmittance we can calculate the emissivity (mean
% absorbance) of the atmospehre

Emissivity = 1 - Mean_transmittance;

%% Figure 2 : Plotting Planck function for the planet along with the extinction coefficient
%%% Here we normalize and scale the Planck function so that the maximum is
%%% the maximum of the emission function

fig2 = figure(2);
t = tiledlayout(3,1);
t.YLabel.String = "Mass Extinction Coefficietn (m__2/kg)";
t.Title.String = "Plots of Molecule's Extinction Coefficients with the Planetary Planck function)";

% top 
nexttile
semilogx(wavelength*10^(6), k_1)
hold on
semilogx(wavelength*10^(6), Spectral_radiance./Den/max(Spectral_radiance./Den)*max(k_1))
title(sprintf("%s Mass Extinction Coefficient", Gas_1))

% middle (1)
nexttile
semilogx(wavelength*10^(6), k_2)
hold on
semilogx(wavelength*10^(6), Spectral_radiance./Den/max(Spectral_radiance./Den)*max(k_2))
title(sprintf("%s Mass Extinction Coefficient", Gas_2))

% bottom
nexttile
semilogx(wavelength*10^(6), k_3)
hold on
semilogx(wavelength*10^(6), Spectral_radiance./Den/max(Spectral_radiance./Den)*max(k_3))
title(sprintf("%s Mass Extinction Coefficient", Gas_3))
%% Emissivity as function of molecule concentration
%%% Here the emissivity is plotted as a function of the concentration of a
%%% greenhouse molecule. You can specify the background concentration of
%%% gases and run multiple different background gas compositions 

%%% Define concentration of the gases that will not have varying
%%% concentrations during the run. Here defining Conc_1_pt2 as water vapor
%%% concentration and Conc_2_pt2 as CH4 concentration since we want CO2 to be
%%% the gas with varying concentration. 

Conc_1_pt2 = 0 ;
Conc_2_pt2 = 0 ;

% Use same method as above to calculate the optical depth of the two
% constant concentration gases 
WeightingCoef_1_pt2 = Conc_1_pt2*density_H2O; 
WeightingCoef_2_pt2 = Pressure*Conc_2_pt2/g; 

Optical_depth_1_pt2 = k_1*WeightingCoef_1_pt2;
Optical_depth_2_pt2 = k_2*WeightingCoef_1_pt2;

%%% For the gas you will be plotting emissivity against, specify the lower
%%% bound and upper bound for the concentration range as well as the
%%% interval you want between the concentration

lowBound    = 0         ;
highBound   = 0.002     ;
interval    = 0.000005  ;

Conc_3_pt2 = lowBound: interval : highBound;

% Create an empty vector of the same length for emissivity values
EmissivtyVect = zeros(length(Conc_3_pt2),1);

% Double for-loop will run through the same method as above for calculating
% the mean atmospheric transmittance but this time will calculation for
% each concentration of the gas in the range defined. 
for i = 1 : length(Conc_3_pt2)
    Weighting_coef = Pressure * Conc_3_pt2(i)/g;
    Optical_depth_3_pt2 = k_3 * Weighting_coef;
    Optical_depth_tot_pt2 = Optical_depth_1_pt2 + Optical_depth_2_pt2 + Optical_depth_3_pt2;
    k_tot_pt2 = Optical_depth_tot_pt2/(Pressure/g);
    Transmittance_tot_pt2 = exp(Optical_depth_tot_pt2 * (-1));
    
    num_pt2 = 0;
    den_pt2 = 0;
    
    for j = 1 : length(wavelength)
        this_num_pt2 = Transmittance_tot_pt2(j)*Spectral_radiance(j)*dlambda(j);
        this_den_pt2 = Spectral_radiance(j)*dlambda(j);
        
        num_pt2 = num_pt2 + this_num_pt2;
        den_pt2 = den_pt2 + this_den_pt2;
    end
    
    Mean_transmittance_pt2 = num_pt2/den_pt2;
    EmissivityVect(i) = 1 - Mean_transmittance_pt2;
end

%% Figure 3: Plot of emissivity as a function of gas concentration (in ppm)     
fig3 = figure(3);
plot(Conc_3_pt2*10^6, EmissivityVect, "Linewidth",2)
hold on
% add lines for pre-industrial CO2 levels (if plotting as function of CO2,
% otherwise comment out this section
xline(280, '--', 'Pre-industrial')
xline(280*2, '-','Double Pre-industrial')
title(sprintf("Atmospheric Emissivity as a functino of %s at 5km", Gas_3), "fontsize", 13)
xlabel(sprintf("%s concentration (ppmv)", Gas_3), "fontsize", 12)
ylabel("Emissvity", "fontsize",12)

%% Figure 4: Plot of emissivity as a function of gas concentration (ppm) on a loglog scale
fig4 = figure(4);
loglog(Conc_3_pt2*10^6, EmissivityVect, "Linewidth",2)
hold on
% add lines for pre-industrial CO2 levels (if plotting as function of CO2,
% otherwise comment out this section
xline(280, '--', 'Pre-industrial')
xline(280*2, '-','Double Pre-industrial')
title(sprintf("Atmospheric Emissivity as a functino of %s at 5km", Gas_3), "fontsize", 13)
xlabel(sprintf("%s concentration (ppmv)", Gas_3), "fontsize", 12)
ylabel("Emissvity", "fontsize",12)
