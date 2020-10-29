% Uniform distribution parameters
a = 0.5;
b = 2;

% Scaling uniform distributions
xc = 0.9;
yc = 0.5;

% Setting up
n = 1e5; nn = 1:n;

rp = @(v) v * (a + rand * (b - a));
dp = @(u, v) (v * a <= u && u <= v * (a + (b - a))) / (v * (b - a));
ldp = @(u, v) log(1 * (v * a <= u && u <= v * (a + (b - a)))) - log(v * (b - a));


type = 'indirect';
[mc, mcl, mcus] = deal(nan(n, 2));
switch type
case 'direct'
    % Sampling x = xc * U(a, b)
    rx = @() rp(xc); dx = @(x) dp(x, xc); ldx = @(x) ldp(x, xc);
    ry = @() rp(yc); dy = @(y) dp(y, yc); ldy = @(y) ldp(y, yc);

    oracle = a * [xc, yc] + (b - a) * rand(n, 2) .* [xc, yc];
    [mc(:, 1), mc(:, 2)] = arrayfun(@(~) maximalCoupling(rx, dx, ry, dy), nn);
    [mcl(:, 1), mcl(:, 2)] = arrayfun(@(~) maximalCouplingLog(rx, ldx, ry, ldy), nn);
    [mcus(:, 1), mcus(:, 2)] = arrayfun(@(~) maximalCouplingUniformScaling(xc, yc, a, b), nn);

    ranges = [a, b]' * [xc, yc];
case 'indirect'
    % Sampling x = 1 - (1 - xc) * U(a, b)
    xc = max(0.01, min(0.99, xc));
    yc = max(0.01, min(0.99, yc));

    rx = @() 1 - rp(1 - xc); dx = @(x) dp(1 - x, 1 - xc); ldx = @(x) ldp(1 - x, 1 - xc);
    ry = @() 1 - rp(1 - yc); dy = @(y) dp(1 - y, 1 - yc); ldy = @(y) ldp(1 - y, 1 - yc);

    oracle = 1 - b * (1 - [xc, yc]) + rand(n, 2) .* (1 - [xc, yc]) * (b - a);
    [mc(:, 1), mc(:, 2)] = arrayfun(@(~) maximalCoupling(rx, dx, ry, dy), nn);
    [mcl(:, 1), mcl(:, 2)] = arrayfun(@(~) maximalCouplingLog(rx, ldx, ry, ldy), nn);

    [x_t, y_t] = arrayfun(@(~) maximalCouplingUniformScaling(1 - xc, 1 - yc, a, b), nn);
    mcus =  1 - [x_t', y_t'];

    ranges = 1 - [b, a]' * (1 - [xc, yc]);
end

figure(1);
for j = 1:2
    subplot(1, 2, j);
    [n1, e1] = histcounts(oracle(:, j), 101, 'Normalization', 'cdf');
    [n2, e2] = histcounts(mc(:, j), 101, 'Normalization', 'cdf');
    [n3, e3] = histcounts(mcl(:, j), 101, 'Normalization', 'cdf');
    [n4, e4] = histcounts(mcus(:, j), 101, 'Normalization', 'cdf');

    plot([e1; e2; e3; e4]', [zeros(1, 4); [n1; n2; n3; n4]'], ':', ...
        'LineWidth', 2); hold on;
    plot(ranges(:, j), [0, 1], '--', 'LineWidth', 2); hold off
    axis('tight');
end
legend

if isequal(type, 'direct') && ~(b * xc < a * yc || b * yc < a * xc)
    if xc <= yc
        ol = (b * xc - a * yc) / (yc * (b - a));
    else
        ol = (b * yc - a * xc) / (xc * (b - a));
    end
elseif isequal(type, 'indirect') && (a * (1 - xc) <= b * (1 - yc) ...
                                     || (a * (1 - yc) <= b * (1 - xc)))
    if xc <= yc
        ol = (b * (1 - yc) - a * (1 - xc)) / ((1 - xc) * (b - a));
    else
        ol = (b * (1 - xc) - a * (1 - yc)) / ((1 - yc) * (b - a));
    end
else
    ol = 0;
end
[ol, mean([mc(:, 1) == mc(:, 2), mcl(:, 1) == mcl(:, 2), mcus(:, 1) == mcus(:, 2)])]

[min(oracle); min(mc); min(mcl); min(mcus)]'
[max(oracle); max(mc); max(mcl); max(mcus)]'

% cov(oracle), cov(mc), cov(mcl), cov(mcus)

figure(2);
subplot(2, 2, 1); plot(oracle(:, 1), oracle(:, 2), 'o');
subplot(2, 2, 2); plot(mc(:, 1), mc(:, 2), 'o');
subplot(2, 2, 3); plot(mcl(:, 1), mcl(:, 2), 'o');
subplot(2, 2, 4); plot(mcus(:, 1), mcus(:, 2), 'o');

figure(3);
subplot(3, 1, 1); plot(sort(oracle), sort(mc) - sort(oracle))
subplot(3, 1, 2); plot(sort(oracle), sort(mcl) - sort(oracle))
subplot(3, 1, 3); plot(sort(oracle), sort(mcus) - sort(oracle))

figure(4);
subplot(2, 1, 1); plot(repmat(sort(oracle), 1, 3), [sort(mc), sort(mcl), sort(mcus)])
subplot(2, 1, 2); plot(repmat(sort(oracle), 1, 2), [sort(mc) - sort(mcl), sort(mc) - sort(mcus)])

figure(5);
plot(sort(repmat(sort(oracle), 1, 3) - [sort(mc), sort(mcl), sort(mcus)]))
