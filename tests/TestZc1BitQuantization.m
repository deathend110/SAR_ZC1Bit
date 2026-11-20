classdef TestZc1BitQuantization < matlab.unittest.TestCase
    % 1-bit 非对称符号量化与伪 1.5-bit 过零置零规则（B0, R1..R9, L1）真值表单元测试。

    methods (TestClassSetup)
        function addRepoPaths(t)
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(repo));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo, "src", "quantization")));
        end
    end

    methods (Test)
        function asymmetricSgnPlusMapsZeroToPlusOne(t)
            S = complex(single([0, -0.5, 0.5, -1e-6]), single([0, 1, -1, 0]));
            U = complex(single(zeros(size(S))));
            q = quantize_1bit_with_u(S, U);
            t.verifyClass(q, "single");
            t.verifyEqual(real(q), single([1, -1, 1, -1]));
            t.verifyEqual(imag(q), single([1, 1, -1, 1]));

            % 加阈值约定 S + U：零点 S+U==0 仍映射为 +1
            S2 = complex([1, -1], [0.5, -0.5]);
            U2 = complex([-1, 0.5], [-0.5, 0.5]);
            q2 = quantize_1bit_with_u(S2, U2);
            t.verifyEqual(real(q2), [1, -1]);
            t.verifyEqual(imag(q2), [1, 1]);

            t.verifyError(@() quantize_1bit_with_u(S, zeros(1, 2)), "zc1bit:ThresholdSizeMismatch");
        end

        function truthTableAllElevenRulesRangeAndAzimuth(t)
            b = [1, 1, -1, -1, 1, -1, 1, 1].';
            S1 = complex(b, -b);

            expected = struct( ...
                "B0", [1,  1, -1, -1,  1, -1,  1,  1].', ...
                "R1", [1,  0,  0,  0,  0,  0,  0,  1].', ...
                "R2", [1,  0, -1, -1,  0, -1,  1,  1].', ...
                "R3", [1,  1,  0, -1,  1,  0,  1,  1].', ...
                "R4", [1,  1, -1, -1,  0, -1,  0,  1].', ...
                "R5", [1,  1, -1,  0,  1,  0,  1,  1].', ...
                "R6", [1,  0, -1, -1,  0, -1,  0,  1].', ...
                "R7", [1,  1,  0,  0,  1,  0,  1,  1].', ...
                "L1", [1,  0, -1,  0,  0,  0,  1,  1].', ...
                "R8", [1,  0,  0, -1,  0,  0,  1,  1].', ...
                "R9", [1,  1, -1,  0,  0,  0,  0,  1].');

            rule_names = string(fieldnames(expected)).';
            for rkey = rule_names
                % 距离向（dim 1，列向量）
                out_r = apply_zc_pseudo15bit(S1, rkey, struct("Axis", "range", "DecisionMode", "parallel"));
                t.verifySize(out_r, size(S1));
                t.verifyEqual(real(out_r), expected.(rkey), "Rule " + rkey + " real range mismatch");
                % 虚部取 -b，对应反向预期值
                exp_im = -expected.(rkey);
                if ismember(rkey, ["R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9"])
                    % 反号通道等价于对称规则置换，逐项核对虚部独立正确性
                    out_im_direct = apply_zc_pseudo15bit(complex(-b, -b), rkey, ...
                        struct("Axis", "range", "DecisionMode", "parallel"));
                    t.verifyEqual(imag(out_r), real(out_im_direct));
                else
                    t.verifyEqual(imag(out_r), exp_im);
                end

                % 方位向（dim 2，行向量）
                out_a = apply_zc_pseudo15bit(S1.', rkey, struct("Axis", "azimuth", "DecisionMode", "parallel"));
                t.verifySize(out_a, size(S1.'));
                t.verifyEqual(out_a, out_r.');
            end
        end

        function parallelVsSerialAndL1Equivalence(t)
            % 序列 [+1, -1, +1, +1]：R1 并行得 [0,0,0,1]，串行得 [0,0,1,1]
            b = [1; -1; 1; 1];
            S1 = complex(b, ones(size(b)));
            y_par = apply_zc_pseudo15bit(S1, "R1", struct("Axis", "range", "DecisionMode", "parallel"));
            y_ser = apply_zc_pseudo15bit(S1, "R1", struct("Axis", "range", "DecisionMode", "serial"));
            t.verifyEqual(real(y_par), [0; 0; 0; 1]);
            t.verifyEqual(real(y_ser), [0; 0; 1; 1]);
            t.verifyEqual(imag(y_par), [1; 1; 1; 1]);

            % L1 定义等价性：y[n]=(b[n]+b[n+1])/2 (n<N), y[N]=b[N]
            rng(11);
            raw = complex(2 * randi([0, 1], 17, 13) - 1, 2 * randi([0, 1], 17, 13) - 1);
            y_l1_par = apply_zc_pseudo15bit(raw, "L1", struct("Axis", "range", "DecisionMode", "parallel"));
            y_l1_ser = apply_zc_pseudo15bit(raw, "L1", struct("Axis", "range", "DecisionMode", "serial"));
            expected_re = [(real(raw(1:end - 1, :)) + real(raw(2:end, :))) / 2; real(raw(end, :))];
            expected_im = [(imag(raw(1:end - 1, :)) + imag(raw(2:end, :))) / 2; imag(raw(end, :))];
            t.verifyEqual(real(y_l1_par), expected_re);
            t.verifyEqual(imag(y_l1_par), expected_im);
            t.verifyEqual(y_l1_ser, y_l1_par);
        end

        function dcRemovalZeroRatioBiasAnd2DJoint(t)
            % 周期平衡序列上 R1/R8/R9 严格零偏置；L1 仅在端点不等时存在 (b[N]-b[1])/(2N) 边界项，端点相等时严格零偏置
            cyc = repmat([1; 1; -1; -1], 6, 4);
            S1 = complex(cyc, -cyc);
            for rkey = ["R1", "R8", "R9"]
                [~, d] = apply_zc_pseudo15bit(S1, rkey, struct("Axis", "range", "RemoveDC", false));
                t.verifyEqual(d.MeanBiasReal, 0, "AbsTol", 1e-12);
                t.verifyEqual(d.MeanBiasImag, 0, "AbsTol", 1e-12);
            end
            [~, d_l1_asym] = apply_zc_pseudo15bit(S1, "L1", struct("Axis", "range", "RemoveDC", false));
            t.verifyEqual(d_l1_asym.MeanBiasReal, -1 / size(cyc, 1), "AbsTol", 1e-12);
            t.verifyEqual(d_l1_asym.MeanBiasImag, 1 / size(cyc, 1), "AbsTol", 1e-12);

            cyc_sym = repmat([1; -1; -1; 1], 6, 4);
            [~, d_l1_sym] = apply_zc_pseudo15bit(complex(cyc_sym, -cyc_sym), "L1", struct("Axis", "range", "RemoveDC", false));
            t.verifyEqual(d_l1_sym.MeanBiasReal, 0, "AbsTol", 1e-12);
            t.verifyEqual(d_l1_sym.MeanBiasImag, 0, "AbsTol", 1e-12);

            % R2 单侧置零含负偏置，开启 RemoveDC 后严格零均值且诊断保留原始偏置
            [y_dc, d2] = apply_zc_pseudo15bit(single(S1), "R2", struct("Axis", "range", "RemoveDC", true));
            t.verifyClass(y_dc, "single");
            t.verifyLessThan(d2.MeanBiasReal, 0);
            t.verifyGreaterThan(d2.ZeroRatio, 0);
            t.verifyEqual(mean(double(real(y_dc(:)))), 0, "AbsTol", 1e-6);
            t.verifyEqual(mean(double(imag(y_dc(:)))), 0, "AbsTol", 1e-6);

            % 长度 1 边界不变性与非法输入拦截
            one = complex(1, -1);
            t.verifyEqual(apply_zc_pseudo15bit(one, "R1", struct("Axis", "range")), one);
            t.verifyError(@() apply_zc_pseudo15bit(complex(0, 1), "R1"), "zc1bit:Invalid1BitInput");

            % 二维联合置零取并集：结果保持 {-1,0,+1} 且尺寸不变
            M = complex([1, -1; 1, 1], [-1, -1; 1, -1]);
            y_2d = apply_zc_pseudo15bit(M, "R1", struct("Axis", "both_union"));
            t.verifySize(y_2d, size(M));
            t.verifyEqual(real(y_2d), [0, 0; 1, 0]);
        end

        function serialOptimizedMatchesNaiveReference(t)
            rng(99);
            x = complex(2 * randi([0, 1], 19, 17) - 1, 2 * randi([0, 1], 19, 17) - 1);
            rules = ["B0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "L1", "R8", "R9"];
            for rkey = rules
                y_r = apply_zc_pseudo15bit(x, rkey, struct("Axis", "range", "DecisionMode", "serial"));
                ref_r = naive_serial_1d_dim1(x, rkey);
                t.verifyEqual(y_r, ref_r, "Serial range mismatch for rule " + rkey);

                y_a = apply_zc_pseudo15bit(x, rkey, struct("Axis", "azimuth", "DecisionMode", "serial"));
                ref_a = naive_serial_1d_dim1(x.', rkey).';
                t.verifyEqual(y_a, ref_a, "Serial azimuth mismatch for rule " + rkey);
            end
        end
    end
end

function y = naive_serial_1d_dim1(x, rule_key)
    re = naive_serial_channel_dim1(real(x), rule_key);
    im = naive_serial_channel_dim1(imag(x), rule_key);
    y = complex(re, im);
end

function v = naive_serial_channel_dim1(x_ch, rule_key)
    switch string(rule_key)
        case "B0"
            flags = [false, false, false, false];
        case "R1"
            flags = [true, true, true, true];
        case "R2"
            flags = [true, false, false, false];
        case "R3"
            flags = [false, true, false, false];
        case "R4"
            flags = [false, false, false, true];
        case "R5"
            flags = [false, false, true, false];
        case "R6"
            flags = [true, false, false, true];
        case "R7"
            flags = [false, true, true, false];
        case "L1"
            flags = [true, false, true, false];
        case "R8"
            flags = [true, true, false, false];
        case "R9"
            flags = [false, false, true, true];
        otherwise
            error("Unexpected rule key: %s", string(rule_key));
    end
    v = x_ch;
    [n_rows, n_cols] = size(v);
    for col = 1:n_cols
        for row = 1:(n_rows - 1)
            a = v(row, col);
            b = v(row + 1, col);
            if a == 1 && b == -1
                if flags(1), v(row, col) = 0; end
                if flags(2), v(row + 1, col) = 0; end
            elseif a == -1 && b == 1
                if flags(3), v(row, col) = 0; end
                if flags(4), v(row + 1, col) = 0; end
            end
        end
    end
end
