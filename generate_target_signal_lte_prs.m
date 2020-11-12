function [tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length)
% generate target signal(lte prs)
% this function is used in tdoa simulation program
%
% [input]
% - ndlrb: number of downlink resource block
% - nprsrb: number of position reference signal resource block
% - subframe_length: subframe length. number of prs signal is proportional to subframe length
%

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC, 'R.5'
% enb = lteRMCDL('R.5', 'FDD', 20);   % Get configuration based on RMC, 'R.5'
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = nprsrb;        % Bandwidth of PRS in resource blocks
% enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks, default = 2
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes
enb.NDLRB = ndlrb;  % default = 15
% enb.CyclicPrefix = 'Extended';
% enb.DuplexMode = 'TDD';
% Warning: Using default value for parameter field TDDConfig (0) 
% Warning: Using default value for parameter field SSC (0) 
enb;

% get ofdm modulation related information
info = lteOFDMInfo(enb);
info;

fs = info.SamplingRate;
nfft = info.Nfft;

resource_grid = lteDLResourceGrid(enb); % Empty resource grid
initial_grid_length = numel(resource_grid);
% size(resource_grid);

% concatenate resource grid
resource_grid = repmat(resource_grid, 1, subframe_length);
prs_symbol = ltePRS(enb); % Map PRS REs
prs_index = ltePRSIndices(enb);

% map prs symbol to resource grid
for n = 1 : subframe_length
    resource_grid(prs_index + initial_grid_length * (n - 1)) = prs_symbol;
end

tx_signal = lteOFDMModulate(enb, resource_grid);        % OFDM modulate
size(resource_grid);
tx_signal_length = size(tx_signal); % column vector

end


