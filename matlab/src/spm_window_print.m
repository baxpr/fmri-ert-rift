function spm_window_print(p_filename)

% Get SPM's "Graphics" figure and print it some useful way

% Append index to filename to make it novel
[p,n,e] = fileparts(p_filename);
k = 0;
while 1
	k = k + 1;
	p_filename = fullfile(p,sprintf('%s_%03d%s',n,k,e));
	if ~exist(p_filename,'file')
		break
	end
end

% Print
print( ...
	spm_figure('FindWin','Graphics'), ...
	['-d' e(2:end)], ...
	p_filename ...
	);
