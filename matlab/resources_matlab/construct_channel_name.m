function channel_name = construct_channel_name(channels, references, separators)

channel_name = cellfun(@(s) cellfun(@(c) [c, s], channels, 'Un', 0), separators, 'Un', 0);
channel_name = cellfun(@(r) cellfun(@(c) [c, r], [channel_name{:}], 'Un', 0), references, 'Un', 0);
channel_name = [channel_name{:}];
end