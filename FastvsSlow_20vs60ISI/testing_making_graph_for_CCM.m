
load('zbrain3D.mat','brain3D','BrainRegions3D');

figure;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%%%% copy the average matrix of the CCMs from the "getting the short data
%%%% together.xls" file and put it in an variable named R.

n=8;

names={'pallium_fasthab','thalamus_fasthab','tectum_fasthab','hindbrain_fasthab','pallium_nonhab','tectum_nonhab','habenula_nonhab','tectum_slopehab'};

% set the source of the lines:
% set the target of the lines:
[row,col] = find(R); %% this only works if I have 0s where there were no connections.

%weights = nonzeros(triu(R,1));
weights = R(find(R));
% create the graph object:
G = graph(row,col,weights,names);
% mark the lines to remove from the graph:
threshold = 0.001; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*10;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = [1134,998,800,640,1124,880,1070,820];
y = [309,332,400,340,390,400,370,450];
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

view(-90,90)
