function state = zc1bit_load_state(path_value, initial)
%ZC1BIT_LOAD_STATE 仅读取检查点 state，不校验版本、路径或时间戳。

state = initial;
if isfile(path_value)
    loaded = load(path_value, "state");
    if ~isfield(loaded, "state")
        error("zc1bit:CheckpointState", "检查点缺少 state 字段：%s", path_value);
    end
    state = loaded.state;
end
end
