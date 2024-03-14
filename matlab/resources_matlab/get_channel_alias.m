function channel_name = get_channel_alias(channel_name)


prefix = {'', 'EEG ', 'EOG '};
suffix = {'', '-Ref', '-Gnd'};

if ~iscell(channel_name)
    channel_name = {channel_name};
end

channel_name = cellfun(@(p) cellfun(@(c) [p, c], channel_name, 'Un', 0), prefix, 'Un', 0);
channel_name = cellfun(@(s) cellfun(@(c) [c, s], [channel_name{:}], 'Un', 0), suffix, 'Un', 0);
channel_name = [channel_name{:}];
channel_name = reduce_channel_name(channel_name);

end