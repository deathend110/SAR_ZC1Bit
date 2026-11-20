function zc1bit_check_parameters(saved, current, context)
%ZC1BIT_CHECK_PARAMETERS 仅核对实际计算参数；旧记录缺失新字段时不阻断恢复。

if nargin < 3 || strlength(string(context)) == 0
    context = "恢复参数";
end

saved = zc1bit_parameters(saved);
current = zc1bit_parameters(current);

if istable(saved) && istable(current)
    saved_vars = string(saved.Properties.VariableNames);
    curr_vars = string(current.Properties.VariableNames);
    for col = curr_vars
        if ~ismember(col, saved_vars)
            continue;
        end
        if ~isequaln(saved.(col), current.(col))
            mismatch(context + "." + col, saved.(col), current.(col));
        end
    end
    return;
end

if ~isstruct(saved) || ~isstruct(current)
    return;
end

names = string(fieldnames(current)).';
for name = names
    if ~isfield(saved, name)
        continue;
    end
    a = saved.(name);
    b = current.(name);
    if (isstruct(a) && isstruct(b)) || (istable(a) && istable(b))
        zc1bit_check_parameters(a, b, context + "." + name);
    elseif ~isequaln(a, b)
        mismatch(context + "." + name, a, b);
    end
end
end

function mismatch(name, old_val, new_val)
error("zc1bit:ParameterMismatch", ...
    "参数 %s 不一致：原值=%s；当前值=%s。", ...
    name, describe_value(old_val), describe_value(new_val));
end

function text = describe_value(value)
if isnumeric(value) || islogical(value)
    if numel(value) <= 16
        text = string(mat2str(value));
    else
        text = sprintf("[%s 大小 %s]", class(value), mat2str(size(value)));
    end
elseif isstring(value) || ischar(value)
    vals = string(value);
    if numel(vals) <= 8
        text = join(vals(:).', ", ");
    else
        text = sprintf("[string 数量 %d]", numel(vals));
    end
else
    text = string(class(value));
end
end
