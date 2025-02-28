%% PLOT PARAMETERS
% Options
print_format = 'png';
print_font = 'Arial';
print_fontsize = 15;
print_size = [5.5 3];
% Plot Parameters
markersize = 5;
fontsize = print_fontsize;
linewidth = 2;
fontname = print_font;
% Set defaul plot parameters
set(groot,'defaultLineLineWidth',linewidth);
set(groot,'defaultLineMarkerSize',markersize);
set(groot,'defaultAxesFontSize',fontsize);
set(groot,'defaultAxesTitleFontSizeMultiplier',1)
set(groot,'defaulttextfontsize',fontsize);
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultAxesBox','on')
set(groot,'DefaultFigureWindowStyle','normal');
% Level Fuel Consumption contour plot [kg/h]
Level_FC{1} = 0:1:5;
Level_FC{2} = 7;
Level_FC{3} = 10:5:30;
% Level BSFC contour plot [g/kWh]
Level_BSFC{1} = 200:5:220;
Level_BSFC{2} = 220:10:240;
Level_BSFC{3} = 240:20:300;
Level_BSFC{4} = 300:50:400;
% Level EM efficiency contour plot [-]
Level_Eff{1} = 0.7:0.05:0.9;
Level_Eff{2} = 0.9:0.01:0.95;
% Level Regeneration contour plot [-]
Level_Reg{1} = 0:0.1:0.8;
Level_Reg{2} = 0.8:0.05:0.95;