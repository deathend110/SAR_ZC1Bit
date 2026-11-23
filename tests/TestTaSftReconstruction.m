classdef TestTaSftReconstruction < matlab.unittest.TestCase
    % 阈值感知 SFT 1-bit 接收端重构（TA-SFT）单元测试。
    % 验证合规输入校验、边界退化等价、5 电平量化契约、逆 Mills 比稳定与端到端断点恢复。

    methods (TestClassSetup)
        function addRepoPaths(t)
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(repo));
            for name = ["workflow", "quantization", "threshold", "sampling", "sarproc", "metrics", "compat"]
                t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo, "src", name)));
            end
        end
    end

    methods (Test)
        function rejectsInvalid1BitInputAndCarrierMismatch(t)
            c1 = complex([1, -1; -1, 1], [1, 1; -1, -1]);
            carrier = complex([0.8, 0.6; -0.6, 0.8], [0.6, -0.8; 0.8, 0.6]);

            % 含 0 或非 ±1 值必须拒绝
            bad_in = complex([1, 0; -1, 1], [1, 1; -1, -1]);
            t.verifyError(@() apply_ta_sft_reconstruct(bad_in, carrier, "TA1_gate_zc"), "zc1bit:Invalid1BitInput");

            % 缺载波或载波尺寸不匹配、含 NaN 必须拒绝
            t.verifyError(@() apply_ta_sft_reconstruct(c1, [], "TA1_gate_zc"), "zc1bit:InvalidCarrier");
            t.verifyError(@() apply_ta_sft_reconstruct(c1, ones(1, 2), "TA2_anchor_zc"), "zc1bit:InvalidCarrier");
            bad_carrier = carrier;
            bad_carrier(1, 1) = NaN;
            t.verifyError(@() apply_ta_sft_reconstruct(c1, bad_carrier, "TA3_mmse_est"), "zc1bit:InvalidCarrier");

            % 非法方法名与负超参必须拒绝
            t.verifyError(@() apply_ta_sft_reconstruct(c1, carrier, "UNKNOWN"), "zc1bit:UnknownTaSftMethod");
            t.verifyError(@() apply_ta_sft_reconstruct(c1, carrier, "TA1_gate_zc", struct("Tau", -0.1)), "zc1bit:InvalidHyperparam");
            t.verifyError(@() apply_ta_sft_reconstruct(c1, carrier, "TA2_anchor_q5", struct("Alpha", NaN)), "zc1bit:InvalidHyperparam");
            t.verifyError(@() apply_ta_sft_reconstruct(c1, carrier, "TA3_mmse_fixed", struct("Gamma", -1)), "zc1bit:InvalidHyperparam");
        end

        function ta1GateZcBoundaryThreeLevelAndAxisSymmetry(t)
            b_col = [1; -1; -1; 1; -1; -1];
            c1 = complex(single(b_col), single(-b_col));
            % 第 1-2 边中点载波小（|u_mid|=0.1<=0.3），第 3-4 与 4-5 边中点载波大（|u_mid|=0.8>0.3）
            u_re = [0.1; 0.1; 0.8; 0.8; 0.8; 0.2];
            u_im = [0.9; 0.9; 0.2; 0.2; 0.2; 0.9];
            carrier = complex(single(u_re), single(u_im));

            % tau = 0 严格退化为 B0
            [y0, d0] = apply_ta_sft_reconstruct(c1, carrier, "TA1_gate_zc", struct("Axis", "range", "Tau", 0));
            t.verifyClass(y0, "single");
            t.verifyEqual(y0, c1);
            t.verifyEqual(d0.ZeroRatio, 0);

            % tau = 1 严格退化为 R1（由于 |u_mid| <= 1）
            y1 = apply_ta_sft_reconstruct(c1, carrier, "TA1_gate_zc", struct("Axis", "range", "Tau", 1));
            yr1 = apply_ta_sft_reconstruct(c1, [], "R1_par_raw", struct("Axis", "range"));
            t.verifyEqual(y1, yr1);

            % tau = 0.3：实部仅 1-2 边置 0，虚部仅 3-4 与 4-5 边置 0
            [yg, dg] = apply_ta_sft_reconstruct(c1, carrier, "TA1_gate_zc", struct("Axis", "range", "Tau", 0.3));
            t.verifyEqual(real(yg), single([0; 0; -1; 1; -1; -1]));
            t.verifyEqual(imag(yg), single([-1; 1; 0; 0; 0; 1]));
            t.verifyTrue(all(ismember(real(yg(:)), single([-1, 0, 1]))));
            t.verifyTrue(all(ismember(imag(yg(:)), single([-1, 0, 1]))));
            t.verifyGreaterThan(dg.ZeroRatio, 0);

            % 方位向转置严格对称
            yg_az = apply_ta_sft_reconstruct(c1.', carrier.', "TA1_gate_zc", struct("Axis", "azimuth", "Tau", 0.3));
            t.verifyEqual(yg_az, yg.');
        end

        function ta2AnchorContinuousAndQ5FiveLevelContract(t)
            % 序列 [+1, -1, +1, +1]：节点 2 同时被边 (1,2) 与边 (2,3) 命中，验证先平均后量化
            b = [1; -1; 1; 1];
            c1 = complex(b, -b);
            % 边 (1,2) 中点 u_mid = 0.7 -> v_1 = -0.7; 边 (2,3) 中点 u_mid = -0.3 -> v_2 = +0.3
            % 节点 1 得 -0.7；节点 2 得 (-0.7 + 0.3)/2 = -0.2；节点 3 得 +0.3；节点 4 保持 +1
            u_re = [0.6; 0.8; -1.4; 0.0];
            u_im = [-0.8; -0.6; 0.2; 0.0];
            carrier = complex(u_re, u_im);

            % alpha = 0 严格退化为 R1
            yr1 = apply_ta_sft_reconstruct(c1, [], "R1_par_raw", struct("Axis", "range"));
            y_a0 = apply_ta_sft_reconstruct(c1, carrier, "TA2_anchor_zc", struct("Axis", "range", "Alpha", 0));
            yq_a0 = apply_ta_sft_reconstruct(c1, carrier, "TA2_anchor_q5", struct("Axis", "range", "Alpha", 0));
            t.verifyEqual(y_a0, yr1);
            t.verifyEqual(yq_a0, yr1);

            % 连续版精确数值检查（alpha = 1.0）
            yc = apply_ta_sft_reconstruct(c1, carrier, "TA2_anchor_zc", struct("Axis", "range", "Alpha", 1.0));
            t.verifyEqual(real(yc), [-0.7; -0.2; 0.3; 1.0], "AbsTol", 1e-12);

            % 5 电平量化版：(-0.7 -> -0.5), (-0.2 -> 0), (+0.3 -> +0.5), (+1 -> +1)
            yq5 = apply_ta_sft_reconstruct(c1, carrier, "TA2_anchor_q5", struct("Axis", "range", "Alpha", 1.0));
            t.verifyEqual(real(yq5), [-0.5; 0.0; 0.5; 1.0], "AbsTol", 1e-12);
            allowed_q5 = [-1, -0.5, 0, 0.5, 1];
            t.verifyTrue(all(ismember(real(yq5(:)), allowed_q5)));
            t.verifyTrue(all(ismember(imag(yq5(:)), allowed_q5)));

            % 随机大矩阵压力测试：绝不能出现 ±0.25 等非 5 电平值，且方位转置等价
            rng(2026);
            rand_1bit = complex(2 * randi([0, 1], 32, 24) - 1, 2 * randi([0, 1], 32, 24) - 1);
            theta = 2 * pi * rand(32, 24);
            rand_car = exp(1i * theta);
            yq5_big = apply_ta_sft_reconstruct(rand_1bit, rand_car, "TA2_anchor_q5", struct("Axis", "range", "Alpha", 2.0));
            t.verifyTrue(all(ismember(real(yq5_big(:)), allowed_q5)));
            t.verifyTrue(all(ismember(imag(yq5_big(:)), allowed_q5)));
            yq5_az = apply_ta_sft_reconstruct(rand_1bit.', rand_car.', "TA2_anchor_q5", struct("Axis", "azimuth", "Alpha", 2.0));
            t.verifyEqual(yq5_az, yq5_big.');
        end

        function ta3MmseStabilityAndBlindGammaEstimation(t)
            % 1. gamma = 0 退化为 B0
            c1 = complex([1, -1; 1, -1], [-1, 1; 1, -1]);
            car = exp(1i * [0.1, 1.2; 2.3, 3.4]);
            y0 = apply_ta_sft_reconstruct(c1, car, "TA3_mmse_fixed", struct("Gamma", 0));
            t.verifyEqual(y0, c1);

            % 2. 极端大 gamma=25 下必须全有限且满足大负 z 逆 Mills 比渐近 |imr(z) + z| < 0.1
            extreme_c1 = complex([1, -1], [-1, 1]);
            extreme_car = complex([-1, 1], [1, -1]); % z = b * gamma * u = -25
            y_ext = apply_ta_sft_reconstruct(extreme_c1, extreme_car, "TA3_mmse_fixed", struct("Gamma", 25));
            t.verifyTrue(all(isfinite(real(y_ext)), "all") && all(isfinite(imag(y_ext)), "all"));
            t.verifyEqual(real(y_ext), [25, -25], "AbsTol", 0.2);

            % 3. 合成圆对称复高斯信号，验证 TA3_mmse_est 无偏恢复真实 gamma（无任何 sigma_hat 传入）
            rng(42);
            nr = 256;
            nc = 256;
            S = complex(randn(nr, nc), randn(nr, nc)); % 单路标准差 sigma = 1
            [c_idx, r_idx] = meshgrid(1:nc, 1:nr);
            theta = 2 * pi * (0.27 * r_idx + 0.13 * c_idx);
            carrier = exp(1i * theta);

            for gamma_true = [0.5, 1.0, 1.8]
                U = gamma_true * carrier;
                b_1bit = quantize_1bit_with_u(S, U);
                [y_est, d_est] = apply_ta_sft_reconstruct(b_1bit, carrier, "TA3_mmse_est");
                t.verifyLessThan(abs(d_est.EstGamma - gamma_true), 0.08);
                t.verifyTrue(all(isfinite(real(y_est)), "all") && all(isfinite(imag(y_est)), "all"));
            end
        end

        function endToEndWorkflowAndParameterOnlyResume(t)
            folder = TestTaSftReconstruction.createFixtureDir("ta_sft_e2e");
            cleanup = onCleanup(@() TestTaSftReconstruction.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_test.mat");
            save(param_file, "-struct", "S60");

            data_dir = fullfile(folder, "raw", "SceneA");
            mkdir(data_dir);
            rng(7);
            raw = complex(single(randn(S60.nrn, 24)), single(randn(S60.nrn, 24)));
            save(fullfile(data_dir, "rstart 1.mat"), "raw", "-v7.3");

            cfg = zc1bit_smoke_config(fullfile(folder, "raw"), param_file);
            cfg.results_root = fullfile(folder, "results");
            cfg.dataset_names = "SceneA";
            cfg.enl_dataset_names = "SceneA";
            cfg.num_samples_per_dataset = 4; % 2 calibration (1,3) / 2 test (2,4)
            cfg.raw_decimation_step = 1;
            cfg.ta_sft.range_allocations = [1, 1];
            cfg.ta_sft.azimuth_allocations = [1, 2];
            cfg.ta_sft.tau_grid = [0.2, 0.6];
            cfg.ta_sft.alpha_grid = [0.5, 1.5];
            cfg.ta_sft.gamma_grid = [0.6, 1.4];
            cfg.ta_sft.phase_mismatch_deg = [0, 30];
            cfg.ta_sft.phase_mismatch_units = "ab_R1A2_azimuth_asft";
            cfg.ta_sft.leakage_groups = ["R1A1", "R1A2"];
            cfg.ta_sft.b0_locked_table = table( ...
                ["rb_R1A1_range_rsft"; "ab_R1A2_azimuth_asft"], ...
                ["B0_par_raw"; "B0_par_raw"], ...
                ["s0_fr0p2_fa0"; "sm3_fr0_fa0p6"], ...
                ["rb_R1A1_range_rsft_s0_fr0p2_fa0"; "ab_R1A2_azimuth_asft_sm3_fr0_fa0p6"], ...
                [0; -3], [0.2; 0], [0; 0.6], [0; 0], ...
                'VariableNames', {'UnitKey', 'VariantKey', 'ParamTag', 'ParamKey', 'STR_dB', 'FrOverBr', 'FaOverBa', 'Phi0'});

            out = run_zc1bit_ta_sft(cfg, "all");
            t.verifyEqual(height(out.search.units), 2);
            t.verifyEqual(height(out.search.locked_params), 14);
            t.verifyEqual(height(out.evaluate.summary_table), 14);
            t.verifyEqual(height(out.evaluate.wilcoxon_table), 22); % (6 vs_B0 + 5 vs_R1) * 2 units
            t.verifyEqual(height(out.evaluate.robustness_table), 14); % 2 angles * 7 methods * 1 unit
            t.verifyTrue(all(ismember(out.evaluate.detail_table.LocalSampleIdx, [2, 4])));
            t.verifyTrue(all(isfinite(out.evaluate.summary_table.SSIM_Mean)));
            t.verifyTrue(all(isfinite(out.evaluate.wilcoxon_table.DeltaSSIM_Mean)));
            t.verifyGreaterThan(height(out.evaluate.leakage_table), 0);

            % 验证路径变更 + NaN 覆盖回波后 100% 命中断点零重算
            moved_raw = fullfile(folder, "raw_moved");
            copyfile(fullfile(folder, "raw"), moved_raw);
            raw = complex(nan(S60.nrn, 24), nan(S60.nrn, 24)); %#ok<NASGU>
            save(fullfile(moved_raw, "SceneA", "rstart 1.mat"), "raw", "-v7.3");

            cfg2 = cfg;
            cfg2.data_root = moved_raw;
            cfg2.runtime.cache_gb = 2;
            out2 = run_zc1bit_ta_sft(cfg2, "all");
            t.verifyEqual(out2.search.locked_params.Calibration_SSIM_Mean, out.search.locked_params.Calibration_SSIM_Mean);
            t.verifyEqual(out2.evaluate.summary_table.SSIM_Mean, out.evaluate.summary_table.SSIM_Mean);
            t.verifyEqual(out2.evaluate.robustness_table.SSIM_Mean, out.evaluate.robustness_table.SSIM_Mean);
        end
    end

    methods (Static)
        function root = workspaceTestTmpRoot()
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            root = fullfile(repo, "tmp", "tests");
        end

        function folder = createFixtureDir(tag)
            root = TestTaSftReconstruction.workspaceTestTmpRoot();
            if ~isfolder(root)
                mkdir(root);
            end
            [~, uniq] = fileparts(tempname(char(root)));
            folder = fullfile(root, string(tag) + "_" + string(uniq));
            mkdir(folder);
        end

        function removeFixture(folder)
            resolved = string(java.io.File(char(folder)).getCanonicalPath());
            root = string(java.io.File(char(TestTaSftReconstruction.workspaceTestTmpRoot())).getCanonicalPath());
            if resolved == root || ~startsWith(lower(resolved), lower(root + string(filesep)))
                error("zc1bit:UnsafeCleanup", "拒绝清理工作区测试临时目录以外的路径：%s", resolved);
            end
            if isfolder(resolved)
                rmdir(resolved, "s");
            end
        end
    end
end

function S60 = smallImagingParameters()
S60 = struct();
S60.nrn = 12;
S60.nan = 10;
S60.C = 3e8;
S60.R0 = 1000;
S60.Fs = 1e6;
S60.Tp = 4e-6;
S60.gama = 1e10;
S60.lambda = 0.03;
S60.v = 100;
S60.prf = 1000;
S60.Ta = 0.004;
S60.fc = 10e9;
S60.B = 1e6;
S60.Bd = 200;
S60.R_total = 6;
S60.A_num = 6;
S60.tnrn = 2 * S60.R0 / S60.C + ((0:S60.nrn - 1).' - floor(S60.nrn / 2)) / S60.Fs;
S60.tnan = ((0:S60.nan - 1) - floor(S60.nan / 2)) / S60.prf;
S60.fnrn = ((0:S60.nrn - 1).' - floor(S60.nrn / 2)) * S60.Fs / S60.nrn;
S60.fnan = ((0:S60.nan - 1).' - floor(S60.nan / 2)) * S60.prf / S60.nan;
end
