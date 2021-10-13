clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

% generate grid
[R, C] = GenerateFermatSpiral(400);
R = R'; C = C';

% nearest-neighbor preconditioning
tic
% route_precon = precondition_route(R,C);
route_precon = [R,C];
toc

% 2opt algorithm
tic
route2Opt = two_opt(route_precon(:,1),route_precon(:,2));
toc

%% compare cost
cost1 = round(sqrt(cost([R,C])));
cost2 = round(sqrt(cost(route_precon)));
cost3 = round(sqrt(cost(route2Opt)));

%% show results
axesRange = max(abs([R(:);C(:)]));
figure(1)
subplot(1,3,1)
plot(R,C,'ko-','MarkerfaceColor','k')
hold on
plot(R(1),C(1),'ro-','MarkerfaceColor','r')
hold off, axis square, title({'original',['distance: ',num2str(cost1)]})
axis(axesRange*[-1 1 -1 1])
subplot(1,3,2)
plot(route_precon(:,1),route_precon(:,2),'ko-','MarkerfaceColor','k')
hold on
plot(route_precon(1,1),route_precon(1,2),'ro-','MarkerfaceColor','r')
hold off, axis square, title({'preconditioned',['distance: ',num2str(cost2)]})
axis(axesRange*[-1 1 -1 1])
subplot(1,3,3)
plot(route2Opt(:,1),route2Opt(:,2),'ko-','MarkerfaceColor','k')
hold on
plot(route2Opt(1,1),route2Opt(1,2),'ro-','MarkerfaceColor','r')
hold off, axis square, title({'2-opt result',['distance: ',num2str(cost3)]})
axis(axesRange*[-1 1 -1 1])

%% PRECONDITIONING
function route = precondition_route(R,C)
numPos = size(R,1);
idx = 2:numPos-1;
idx_nearest = [1;zeros(numPos-2, 1);numPos];
for k=2:numPos-1
    distances = (R(idx_nearest(k-1))-R(idx)).^2 + (C(idx_nearest(k-1))-C(idx)).^2;
    next_idx = find(distances == min(distances));
    idx_nearest(k) = idx(next_idx(1));
    idx = setdiff(idx, idx(next_idx(1)));
end
R = R(idx_nearest);
C = C(idx_nearest);
route = [R,C];
end
%% 2OPT
function route = two_opt(R,C)

route = [R,C];
best_route = route;
cost_best = cost(best_route);
progress_made = true;
while progress_made
    progress_made = false;
    for k = 1:length(route)-1
        for l=k+1:length(route)
            % reverse potential crossing
            new_route = [route(1:k,:); ...
                        flipud(route(k+1:l-1,:));...
                        route(l:end,:)];
            % update best route (if lowest cost achieved)
            cost_new = cost(new_route);
            if cost_new < cost_best
                best_route = new_route;
                cost_best = cost_new;
                progress_made = true;
            end
            
        end
    end
    route = best_route;
end
end
%% COST
function r = cost(route)
r = sum(diff(route(:,1)).^2 + diff(route(:,2)).^2);
end