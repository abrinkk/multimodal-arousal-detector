function channel_name = reduce_channel_name(channel_name)

separators = {':', ';', '.', ',', '*', '+', '?', '#', ...
    '(', ')', '[', ']','-', '_', '/' ,'\'};

if ~iscell(channel_name)
    channel_name = {channel_name};
end

channel_name = replace(channel_name, separators, ' ');
channel_name = regexprep(channel_name,'\s+',' ');
channel_name = lower(cellfun(@(c) strtrim(c), channel_name, 'Un', 0));
channel_name = erase(channel_name,[" ","-"]); % test
end