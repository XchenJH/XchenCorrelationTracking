function i2cc = load_channel_matrix(features)

if any(strcmpi('grayCC', features))
    temp = load('i2ccs');
    i2cc = temp.i2ccs;
elseif any(strcmpi('grayCCr', features))
    temp = load('i2ccs_rot');
    i2cc = temp.i2ccs_rot;
else
    i2cc = [];
end

end