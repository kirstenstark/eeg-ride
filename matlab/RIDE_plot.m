function f = RIDE_plot(results,name,channel,varargin)



color = {[0,0,0],[0,0,1],[1,0,0],[0,0.5,0],[0,1,0]};

h = gca;
if ~isempty(varargin) h = varargin{1}; end


t_axis = linspace(results.cfg.epoch_twd(1),results.cfg.epoch_twd(2),size(results.erp,1));

if length(name) == 1 && strcmpi(name{1},'r') t_axis = t_axis-mean(results.latency0{end});end

hold off;
% figure;
for j = 1:length(name)
    eval(['plot(h,t_axis,results.',name{j},'(:,channel),''color'',color{j});']);
    hold on;
end

for j = 1:length(name) name{j}(name{j}=='_')=' ';end
legend(h,name);
axis tight;xlabel('time after stimulus (ms)');ylabel('potential (\muV)');
if length(name) == 1 && strcmpi(name{1},'r') xlabel('time after RT (ms)');end




    