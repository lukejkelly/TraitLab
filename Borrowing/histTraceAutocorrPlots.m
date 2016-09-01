function histTraceAutocorrPlots(outputInfo, destInfo, sampInds)
% Function to produce summaries of MCMC output.

% Setting up.
GlobalSwitches; GlobalValues;

% Variable indices.
lInd = 2;   rInd = 3;   mInd = 4;   kInd = 8;   cInd = 10;  bInd = 12;

% Gathering output.
nOut  = size({outputInfo.path}, 2);
stats = cell(nOut, 1);
for i = 1:nOut
    stats{i} = getfield(readoutput(outputInfo(i).path), 'stats')';
end

% Populating output arrays.
[xLike, xRoot, xMu, xBM, xKap, xDur, xCat] = deal(zeros(length(sampInds), nOut));
for i = 1:nOut
    xLike(:, i) = stats{i}(sampInds, lInd);
    xRoot(:, i) = stats{i}(sampInds, rInd);
    xMu  (:, i) = exp(stats{i}(sampInds, mInd));
    xBM  (:, i) = exp(diff(stats{i}(sampInds, [mInd, bInd]), [], 2));
    xKap (:, i) = stats{i}(sampInds, kInd);
    xCat (:, i) = stats{i}(sampInds, cInd);
    xDur (:, i) = -log(1 - xKap(:, i)) ./ xMu(:, i) .* xCat(:, i);    
end

% Creating plots.

% Log-likelihood.
ESSLike = ESSFun(xLike);

histFun(xLike, 'Log-likelihood')
printVarieties(gcf, 'histLike', outputInfo, destInfo, ESSLike)

traceFun(xLike, 'Log-likelihood', sampInds)
printVarieties(gcf, 'traceLike', outputInfo, destInfo, ESSLike)

autocorrFun(ESSLike, 'Log-likelihood')
printVarieties(gcf, 'autocorrLike', outputInfo, destInfo, ESSLike)

% Root age.
ESSRoot = ESSFun(xRoot);

histFun(xRoot, 'Root age $ -t_1 $')
printVarieties(gcf, 'histRoot', outputInfo, destInfo, ESSRoot)

traceFun(xRoot, 'Root age $ -t_1 $', sampInds)
printVarieties(gcf, 'traceRoot', outputInfo, destInfo, ESSRoot)

autocorrFun(ESSRoot, 'Root age $ -t_1 $')
printVarieties(gcf, 'autocorrRoot', outputInfo, destInfo, ESSRoot)

% Death rate.
ESSMu = ESSFun(xMu);

histFun(xMu, 'Death rate $ \mu $')
printVarieties(gcf, 'histMu', outputInfo, destInfo, ESSMu)

traceFun(xMu, 'Death rate $ \mu $', sampInds)
printVarieties(gcf, 'traceMu', outputInfo, destInfo, ESSMu)

autocorrFun(ESSMu, 'Death rate $ \mu $')
printVarieties(gcf, 'autocorrMu', outputInfo, destInfo, ESSMu)

% Relative transfer rate.

% If any SD samples then set corresponding column of xBM to 'NaN'.
ESSBM = ESSFun(xBM);

histFun(xBM, 'Relative transfer rate $ \beta / \mu $')
printVarieties(gcf, 'histBM', outputInfo, destInfo, ESSBM)

traceFun(xBM, 'Relative transfer rate $ \beta / \mu $', sampInds)
printVarieties(gcf, 'traceBM', outputInfo, destInfo, ESSBM)

autocorrFun(ESSBM, 'Relative transfer rate $ \beta / \mu $')
printVarieties(gcf, 'autocorrBM', outputInfo, destInfo, ESSBM)

% Catastrophe severity.
ESSKap = ESSFun(xKap);

histFun(xKap, 'Catastrophe severity $ \kappa $')
printVarieties(gcf, 'histKappa', outputInfo, destInfo, ESSKap)

traceFun(xKap, 'Catastrophe severity $ \kappa $', sampInds)
printVarieties(gcf, 'traceKappa', outputInfo, destInfo, ESSKap)

autocorrFun(ESSKap, 'Catastrophe severity $ \kappa $')
printVarieties(gcf, 'autocorrKappa', outputInfo, destInfo, ESSKap)

% Catastrophe duration.
ESSDur = ESSFun(xDur);

histFun(xDur, 'Total catastrophe length $ \delta * |C| $')
printVarieties(gcf, 'histDur', outputInfo, destInfo, ESSDur)

traceFun(xDur, 'Total catastrophe length $ \delta * |C| $', sampInds)
printVarieties(gcf, 'traceDur', outputInfo, destInfo, ESSDur)

autocorrFun(ESSDur, 'Total catastrophe length $ \delta * |C| $')
printVarieties(gcf, 'autocorrDur', outputInfo, destInfo, ESSDur)

% Number of catastrophes.
ESSCat = ESSFun(xCat);

histFun(xCat, 'Number of catastrophes $ |C| $')
printVarieties(gcf, 'histCat', outputInfo, destInfo, ESSCat)

traceFun(xCat, 'Number of catastrophes $ |C| $', sampInds)
printVarieties(gcf, 'traceCat', outputInfo, destInfo, ESSCat)

autocorrFun(ESSCat, 'Number of catastrophes $ |C| $')
printVarieties(gcf, 'autocorrCat', outputInfo, destInfo, ESSCat)

end


% ESS function.
function [ESS] = ESSFun(X)

% Number of sequences and lags.
nSeq = size(X, 2);
numLags = round(length(X) / 3); % roundn(length(X) / 3, 3);

% ESS struct inputs.
[acf, dispNegInd, ess] = deal(cell(nSeq, 1));

for i = 1:nSeq
    acf{i}        = autocorr(X(:, i), numLags);
    essNI         = essNegInd(acf{i});
    dispNegInd{i} = min([2 * essNI, numLags]);
    ess{i}        = round(length(X) / (1 + 2 * sum(acf{i}(2:(essNI + 1)))));
end
    
% Creating output struct.
ESS = struct('acf', acf, 'negInd', dispNegInd, 'ess', ess);

end

% Function to obtain index at which ACF becomes negligible.
function [negInd] = essNegInd(gamma)

% Using method of Geyer (1992).
Gamma = gamma(1:2:(end - 1)) + gamma(2:2:end);

% First index that Gamma is negative.
zeroInd = find(Gamma < 0, 1);

% Index after which acf entries are 'negligible'.
negInd = 2 * (zeroInd - 1) + 1;

end

% Histogram function.
function [] = histFun(Y, xLabel)

% Setting up.
sf = 1.01;

% Binning data and creating bar plot.
if all(rem(Y, 1) == 0)
    [N, X] = hist(Y, 0:(max(Y(:))), 'Visible', 'off');
    xLims  = [0, max(Y(:)) + 1];
else
    [N, X] = hist(Y, 20, 'Visible', 'off');
    xLims  = [min(Y(:)), max(Y(:))];
end
bar(X, bsxfun(@rdivide, N, sum(N) * diff(X(1:2))), 'BarWidth', 1, 'LineWidth', 0.25);

% Plotting window parameters.
xlim(xLims);
ylim([0, max(N(:)) / (length(Y) * diff(X(1:2))) * sf]);
xlabel(xLabel, 'Interpreter', 'LaTeX');
ylabel('Relative frequency', 'Interpreter', 'LaTeX');

end

% Traceplot function.
function [] = traceFun(X, yLabel, sampInds)

% Subsample to save space, if necessary.
% sampInds = sampInds(100:100:end);
% X        = X(100:100:end, :);

yLims = [min(X(:)), max(X(:))];

% Jitter output if X is integral.
if all(mod(X, 1) == 0)
    X = X - 0.025 + rand(size(X)) / 20; ...
    yLims = [floor(min(X(:))), ceil(max(X(:)))];
end

% Binning data and creating bar plot.
plot(sampInds, X, 'MarkerSize', 6);

% Plotting window parameters.
xlim(sampInds([1, end]));
ylim(yLims);
xlabel('Iteration', 'Interpreter', 'LaTeX');
ylabel(yLabel, 'Interpreter', 'LaTeX');
xTickLocs = roundn(sampInds(1), 1):roundn(sampInds(end) / 4, 3):...
    roundn(sampInds(end), 1);
set(gca, 'XTick', xTickLocs);
set(gca, 'XTickLabel', sprintf('%2.1e\n', xTickLocs));

end

% Autocorrelation plot function.
function [] = autocorrFun(ESS, yLabel)

% Finding index after which both entries become negligible.
lagInds = 1:max([1, 10^(floor(log10(max([ESS.negInd]))) - 1)]):max([ESS.negInd]);

% Creating stem plot.
Y = [ESS.acf];
plot(lagInds - 1, Y(lagInds, :));

% Plotting window parameters.
xlim(lagInds([1, end]) - 1);
ylim([min(min([ESS(:).acf])), 1]);
xlabel('Lag', 'Interpreter', 'LaTeX');
ylabel([yLabel, ' sample ACF'], 'Interpreter', 'LaTeX');

end


% Add legend and print colour/b&w and square/rectangle varieties of plot.
function [] = printVarieties(fig, plotName, outputInfo, destInfo, ESS)

% Plotting window parameters.
fName   = {'PaperPosition', 'PaperSize'};
fSquare = {[0 0 10 10], [10 10]};
set(fig, fName, fSquare); 

% Legend entries.
labs = cell(size({outputInfo.label}, 2), 1);
for i = 1:size(labs, 1)
    labs{i} = sprintf('%s (%d)', outputInfo(i).label, ESS(i).ess);
end

% Legend.
legend(labs, 'location', 'NorthEast');
set(legend, 'Interpreter', 'LaTeX');

% Printing.
pause(0.5); print(fig, '-dpdf', [destInfo, plotName, '.pdf']);

% Clear figure in preparation for next round.
clf;

end
