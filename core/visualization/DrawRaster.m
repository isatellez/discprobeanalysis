function [varargout] = DrawRaster(spkcell,mrkcell,groupingvar,rowcolors)
% Draw a raster plot of neural spiking data for a single unit
% Syntax:
%       [H, hspk, hbhv, hzro] = DrawRaster(spkcell,mrkcell,groupingvar,rowcolors)
%
% Inputs:
%    spkcell - cell array of length==trials, each cell contains spike times
%    mrkcell - cell array of length==trials, each cell contains marker times
%    groupingvar(optional) - vector of IDs to group trials
%    rowcolors(optional) - Ntrials x 3 RGB for per-row backgrounds
%
% By RAMON BARTOLO [2016], lightly adjusted for rowcolors.

if nargin==1 && ~iscell(spkcell) && strcmp(spkcell,'test')
    [spkcell,mrkcell] = genTestData; %#ok<NASGU>
end

if nargin<3
    groupingvar(:,1:numel(spkcell)) = 1;
else
    if numel(spkcell) ~= numel(groupingvar)
        error('groupingvar must have as many elements as cells with spike times');
    end
    if numel(spkcell) ~= numel(mrkcell)
        error('markers must have as many cells as the spike times array');
    end
end

useColors = (nargin>=4) && ~isempty(rowcolors);
if useColors
    rowcolors = double(rowcolors);

    if max(rowcolors,[],'all') > 1
        rowcolors = rowcolors ./ 255;
    end

    if size(rowcolors,2) ~= 3
        error('rowcolors must be Ntrials x 3');
    end
end

maxargout = 4;

group = unique(groupingvar);
ngrps = length(group(:,1));

TrIdx = 0;
xvec = [];
yvec = [];
bhvvec = [];
row_order = zeros(numel(spkcell),1);

for gr = 1:ngrps
    idx = find(groupingvar==group(gr));
    gsc = spkcell(idx);
    gmc = mrkcell(idx);
    grtrs = numel(gsc);

    for tr = 1:grtrs
        if ~isempty(gsc{tr})
            txvec(1,:) = gsc{tr};
            tyvec(1,:) = zeros(size(txvec)) + TrIdx + 0.5;
            xvec = cat(2,xvec,txvec);
            yvec = cat(2,yvec,tyvec);
        end

        tbhvvec(1,:) = gmc{tr};
        tbhvvec(tbhvvec==0) = [];
        tbhvvec(2,:) = TrIdx+0.5;
        bhvvec = cat(2,bhvvec,tbhvvec);

        TrIdx = TrIdx + 1;
        row_order(TrIdx) = idx(tr);
        clear txvec tyvec tbhvvec
    end

    if ngrps>1 && gr~=ngrps
        divliney(gr,1:2) = [TrIdx, TrIdx]; %#ok<AGROW>
    end

    clear gsc gmc grtrs idx
end

h(1) = figure;
hold on;
xlabel("Time (s)");
ylabel("Trials");

% per-row background in the first 200 ms
if useColors
    xs = [];
    if ~isempty(xvec),   xs = [xs xvec];        end
    if ~isempty(bhvvec), xs = [xs bhvvec(1,:)]; end
    if isempty(xs), xs = [0 1]; end
    xl = [min(xs) max(xs)];

    colorWin = [0 0.20];
    x1 = max(colorWin(1), xl(1));
    x2 = min(colorWin(2), xl(2));

    if x2 > x1
        for r = 1:TrIdx
            c = rowcolors(row_order(r),:);
            patch([x1 x2 x2 x1], [r-0.5 r-0.5 r+0.5 r+0.5], c, ...
                  'EdgeColor','none','FaceAlpha',0.25);
        end
    end
end

h(2) = plot(xvec,yvec,'|','Color','k','MarkerSize',1);
if ~isempty(bhvvec)
    h(3) = plot(bhvvec(1,:),bhvvec(2,:),'.','Color',[0.9,0.2,0.1],'MarkerSize',1);
end
axis('tight');
h(4) = line([0,0],[0,TrIdx],'Color',[0.9,0.2,0.2],'LineWidth',1);

if ngrps>1
    A = get(gca,'XLim');
    divlinex = [ones(ngrps-1,1)*A(1), ones(ngrps-1,1)*A(2)];
    if exist('divliney','var') && ~isempty(divlinex)
        line(divlinex',divliney','Color','k','Linewidth',1);
    end
end

set(gca,'YTick',[1,TrIdx],'TickDir','out');
hold off;

na = nargout;
if na > 0
    if na>maxargout
        warning('too many output arguments, ignoring the extra arguments');
        na = maxargout;
    end
    for k = 1:na
        varargout{k} = h(k);
    end
end

end
