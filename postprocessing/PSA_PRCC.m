% Function to run the model with parameter snesitivity analysis cases
%
% Inputs: model        -- simbio data from PSA
%         params       -- object containing LHS values of model parameters
%         varargin     -- flag that indicates if this analysis is done on
%                         plausible patients or all the patients
%
% Outputs: simDataPSA  -- results from the LHS simulations


function PSA_PRCC(params_in,params_out,varargin)

% find if all the simulations are used or just the plausible ones
if (nargin == 3)
    if strcmp(varargin{1},'plausible')
        n_PSA = length(params_out.iPatientPlaus);
        index = params_out.iPatientPlaus;
    elseif strcmp(varargin{1},'patient')
        n_PSA = length(params_out.iPatient);
        index = params_out.iPatient;
    else
        n_PSA = length(params_out.iPatient);
        index = params_out.iPatient;
    end
else
    n_PSA = length(params_out.iPatient);
    index = params_out.iPatient;
end

% load parameter labels
for i = 1:length(params_in.names)
    namesIn{i} = params_in.(params_in.names{i}).ScreenName;
end

for i = 1:length(params_out.names)
    namesOut{i} = params_out.(params_out.names{i}).ScreenName;
end

% calculate PRCCs
alpha = 0.05/size(params_in.all(index,:), 2);
[rho,pval,~] = PRCC(params_in.all(index,:), params_out.post(index,:), 1, namesIn, alpha);
% by Simeone Marino, May 29 2007
% modified by Marissa Renardy, October 8 2020

% replace NaNs to zero to be able to heatmap
k = find(isnan(rho))';
rho(k) = 0;

figure
h = imagesc(rho);
hold on

[rows,cols] = size(rho);
for i = 1:rows
    for j = 1:cols
        if pval(i,j) <= 0.001 / size(params_in.all(index,:), 2)
            textHandles(i,j) = text(j,i,'***',...
                'horizontalAlignment','center', 'FontSize',14);
        elseif pval(i,j) <= 0.01 / size(params_in.all(index,:), 2)
            textHandles(i,j) = text(j,i,'**',...
                'horizontalAlignment','center', 'FontSize',14);
        elseif pval(i,j) <= 0.05 / size(params_in.all(index,:), 2)
            textHandles(i,j) = text(j,i,'*',...
                'horizontalAlignment','center', 'FontSize',14);
        end
    end
end

ax = gca;
ax.XTick = 1:length(namesIn);
ax.XTickLabel = namesIn;
ax.YTick = 1:length(namesOut);
ax.YTickLabel = namesOut;
ax.TickLength = [0 0];

set(gcf, 'Position', [100 200 1665 830]);
set(gca,'YTickLabelRotation',35)
set(gca,'XTickLabelRotation',35)
set(gca,'FontSize',12, 'Position', [.27 .55 .60 .40]);

mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
% - Author:   Víctor Martínez-Cagigal
% - Date:     19/11/2018
% - Version:  1.0
% - E-mail:   vicmarcag (at) gmail (dot) com
% - Biomedical Engineering Group (University of Valladolid), Spain

% colorbar('southoutside');
colormap(mycolormap);
hBar = colorbar;
caxis([-1 1])
set(hBar, 'Position', [.90 .55 .02 .4]);
set(gca,'FontSize',14);
