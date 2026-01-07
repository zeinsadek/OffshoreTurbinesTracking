%% Demo for plotBarStackGroups(stackData, groupLabels)

clear; close all; clc
% addpath('/Users/zeinsadek/Documents/MATLAB/MatlabFunctions/plotBarStackGroups')

% Dimensions expected by your function:
% stackData(i,j,k) = (Group, Stack, StackElement)
NumGroups = 5;   % number of x-axis groups
NumStacks = 4;   % number of bars per group (side-by-side)
NumElems  = 3;   % number of stacked pieces per bar


% Make toy data (positive values so stacks look nice)
% rng(0)
% stackData = 0.5 + 3*rand(NumGroups, NumStacks, NumElems);

% Structured fake data to troubleshoot
stack_values = [1,2,3];
for g = 1:NumGroups
    for s = 1:NumStacks
        for e = 1:NumElems
            stackData(g,s,e) = g * stack_values(e);
        end
    end
end


% Color control
barColors   = lines(NumStacks);   % one color per BAR
stackAlpha  = linspace(0.4,1,size(stackData,3));  % opacity per STACK

% Labels for groups
% pretend these are your x values (intentionally unsorted)
xVals = [50 45 40 35 30];
% xVals = fliplr(xVals);
groupLabels = arrayfun(@(x) sprintf("S%g", x), xVals, 'UniformOutput', false);


% Sort by xVals (or use your own custom order)
[xSorted, idx] = sort(xVals, 'ascend');

stackData2  = stackData(idx,:,:);
labels2     = groupLabels(idx);

% Sanity printout (this is the "confidence" step)
disp(table((1:NumGroups)', xSorted(:), string(labels2(:)), ...
    'VariableNames', {'RowInPlot','x','Label'}))


% Call your function
figure('color', 'white')
hold on
plotBarStackGroups(stackData2, barColors, stackAlpha, labels2);
title('Grouped stacked bars using plotBarStackGroups()')
ylabel('Value')



% Create dummy legend handles (patches) that won't affect the plot
dummy = gobjects(NumStacks,1);
for i = 1:NumStacks
    dummy(i) = patch(NaN, NaN, barColors(i,:), ...
        'EdgeColor','none', ...
        'FaceAlpha',1);   % main color only
end
hold off

% Create legend using ONLY these handles
legend(dummy, {'Case 1','Case 2','Case 3','Case 4'}, 'Location','northwest', 'box', 'off');




%%

function [] = plotBarStackGroups(stackData, barColors, stackAlpha, groupLabels)
% Plot a set of stacked bars, but group them according to labels provided.
%
% Params: 
%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)

NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
figure('color', 'white')
hold on
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked', 'HandleVisibility', 'off');
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    
    % ---- color stacks within THIS bar ----
    for k = 1:size(Y,2)   % stack elements
        h(i,k).FaceColor = barColors(i,:);
        h(i,k).FaceAlpha = stackAlpha(k);
        h(i,k).EdgeColor = 'none';
    end

end
hold off
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end 
