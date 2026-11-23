classdef TestSftP15Reconstruction < matlab.unittest.TestCase
    % 严格伪 1.5-bit（{+1, 0, -1}）SFT 重构算法契约、真值表与工作流断点测试。

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
        function rejectsInvalidInputsAndEnforcesThreeLevelContract(t)
            rng(2026);
            re_in = single(2 * randi([0, 1], 16, 12) - 1);
            im_in = single(2 * randi([0, 1], 16, 12) - 1);
            c1 = complex(re_in, im_in);
            ph = 2 * pi * rand(16, 12);
            carrier = exp(1i * ph);

            % 1. 非 {+1,-1} 输入拒绝
            bad_in = c1;
            bad_in(1, 1) = 0;
            t.verifyError(@() apply_sft_p15_reconstruct(bad_in, carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 0.5)), "zc1bit:Invalid1BitInput");

            % 2. 缺载波/错尺寸/非有限载波校验
            t.verifyError(@() apply_sft_p15_reconstruct(c1, [], "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 0.5)), "zc1bit:InvalidCarrier");
            t.verifyError(@() apply_sft_p15_reconstruct(c1, ones(2, 2), "P15_3_arcsin_lin", ...
                struct("Axis", "range", "TauLin", 0.5)), "zc1bit:InvalidCarrier");
            bad_car = carrier;
            bad_car(3, 3) = NaN;
            t.verifyError(@() apply_sft_p15_reconstruct(c1, bad_car, "P15_4_dual_domain", ...
                struct("Axis", "range", "Tau", 0.8, "MinRunLength", 3)), "zc1bit:InvalidCarrier");

            % 3. 全盲方案允许空载波
            y_p2_blind = apply_sft_p15_reconstruct(c1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 1, "MinRunLength", 2));
            t.verifyClass(real(y_p2_blind), "single");

            % 4. 非法方法名与超参拒绝
            t.verifyError(@() apply_sft_p15_reconstruct(c1, carrier, "P15_unknown", struct("Axis", "range")), ...
                "zc1bit:UnknownMethod");
            t.verifyError(@() apply_sft_p15_reconstruct(c1, carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", -0.1)), "zc1bit:InvalidHyperparam");
            t.verifyError(@() apply_sft_p15_reconstruct(c1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 2, "MinRunLength", 1)), "zc1bit:InvalidHyperparam");

            % 5. 全部 6 方法输出严格属于 {-1, 0, +1} 且保留 single 精度
            specs = { ...
                "B0_par_raw", struct("Axis", "range"); ...
                "R1_par_raw", struct("Axis", "range"); ...
                "P15_1_carrier_agree", struct("Axis", "range", "Tau", 0.5); ...
                "P15_2_rlp_she", struct("Axis", "range", "HalfWidth", 2, "MinRunLength", 5); ...
                "P15_3_arcsin_lin", struct("Axis", "range", "TauLin", 0.4); ...
                "P15_4_dual_domain", struct("Axis", "range", "Tau", 0.7, "MinRunLength", 3)};
            for k = 1:size(specs, 1)
                out = apply_sft_p15_reconstruct(c1, carrier, specs{k, 1}, specs{k, 2});
                t.verifyClass(real(out), "single");
                t.verifySize(out, size(c1));
                t.verifyTrue(all(real(out) == -1 | real(out) == 0 | real(out) == 1, "all"));
                t.verifyTrue(all(imag(out) == -1 | imag(out) == 0 | imag(out) == 1, "all"));
            end
        end

        function p15_1CarrierAgreeDegeneracyAndOddSymmetry(t)
            b = [1; 1; -1; -1];
            u = [0.8; 0.3; -0.8; 0.3];
            c1 = complex(b, -b);
            carrier = complex(u, -u);

            % tau = 1.0 严格退化为 B0
            y_b0 = apply_sft_p15_reconstruct(c1, carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 1.0));
            t.verifyEqual(y_b0, c1);

            % 手算真值：b.*u = [0.8, 0.3, 0.8, -0.3] > 0.5 => [0, 1, 0, -1]
            [y, diag] = apply_sft_p15_reconstruct(c1, carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 0.5));
            t.verifyEqual(real(y), [0; 1; 0; -1]);
            t.verifyEqual(imag(y), [0; -1; 0; 1]);
            t.verifyEqual(diag.ZeroRatio, 0.5, "AbsTol", 1e-12);

            % 奇对称校验：(-b, -u) -> -y
            y_neg = apply_sft_p15_reconstruct(-c1, -carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 0.5));
            t.verifyEqual(y_neg, -y);
        end

        function p15_2RlpSheTruthTableHalfWidthAndAxisSymmetry(t)
            % 游程依次为：3, 1, 1, 2, 6, 5；仅第 5 条边 (idx 13->14) 左右均 >= 3 且 >= 5
            col = [1; 1; 1; -1; 1; -1; -1; 1; 1; 1; 1; 1; 1; -1; -1; -1; -1; -1];
            S1 = complex([col, -col], [col, -col]);

            % 1. Lmin = Inf 严格退化为 B0
            [y_inf, d_inf] = apply_sft_p15_reconstruct(S1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 1, "MinRunLength", Inf));
            t.verifyEqual(y_inf, S1);
            t.verifyEqual(d_inf.ProtectedEdgeRatio, 1.0, "AbsTol", 1e-12);

            % 2. (w=1, Lmin=1) 严格退化为 R1_par_raw
            y_l1 = apply_sft_p15_reconstruct(S1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 1, "MinRunLength", 1));
            y_r1 = apply_sft_p15_reconstruct(S1, [], "R1_par_raw", struct("Axis", "range"));
            t.verifyEqual(y_l1, y_r1);

            % 3. (w=1, Lmin=3)：5 条边仅边 5 置零，其余 4 条受保护 (4/5 = 0.8)
            [y_w1, d_w1] = apply_sft_p15_reconstruct(S1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 1, "MinRunLength", 3));
            exp_col = col;
            exp_col([13, 14]) = 0;
            t.verifyEqual(real(y_w1(:, 1)), exp_col);
            t.verifyEqual(real(y_w1(:, 2)), -exp_col);
            t.verifyEqual(d_w1.ProtectedEdgeRatio, 0.8, "AbsTol", 1e-12);

            % 4. (w=2, Lmin=5)：边 5 两侧各置 2 点 (idx 12, 13, 14, 15)
            y_w2 = apply_sft_p15_reconstruct(S1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 2, "MinRunLength", 5));
            exp_w2 = col;
            exp_w2([12, 13, 14, 15]) = 0;
            t.verifyEqual(real(y_w2(:, 1)), exp_w2);
            t.verifyEqual(real(y_w2(:, 2)), -exp_w2);

            % 5. 距离/方位转置等价性（含多列跨列切断检验）
            y_az = apply_sft_p15_reconstruct(S1.', [], "P15_2_rlp_she", ...
                struct("Axis", "azimuth", "HalfWidth", 2, "MinRunLength", 5));
            t.verifyEqual(y_az, y_w2.');
        end

        function p15_3ArcsinLinDirectionAndOddSymmetry(t)
            b = [1; -1; 1; -1; 1];
            u = [-0.8; -0.8; 0.8; 0.8; 0.2];
            c1 = complex(b, -b);
            carrier = complex(u, -u);

            % TauLin = 1.0 严格退化为 B0
            y_b0 = apply_sft_p15_reconstruct(c1, carrier, "P15_3_arcsin_lin", ...
                struct("Axis", "range", "TauLin", 1.0));
            t.verifyEqual(y_b0, c1);

            % 边 1 (u_mid=-0.8) 削 +1 侧 (idx 1); 边 3 (u_mid=+0.8) 削 -1 侧 (idx 4)
            y = apply_sft_p15_reconstruct(c1, carrier, "P15_3_arcsin_lin", ...
                struct("Axis", "range", "TauLin", 0.5));
            t.verifyEqual(real(y), [0; -1; 1; 0; 1]);

            % 奇对称与方位转置等价
            y_neg = apply_sft_p15_reconstruct(-c1, -carrier, "P15_3_arcsin_lin", ...
                struct("Axis", "range", "TauLin", 0.5));
            t.verifyEqual(y_neg, -y);
            y_az = apply_sft_p15_reconstruct(c1.', carrier.', "P15_3_arcsin_lin", ...
                struct("Axis", "azimuth", "TauLin", 0.5));
            t.verifyEqual(y_az, y.');
        end

        function p15_4DualDomainUnionAndDegeneracy(t)
            rng(99);
            S1 = complex(2 * randi([0, 1], 24, 6) - 1, 2 * randi([0, 1], 24, 6) - 1);
            carrier = exp(1i * 2 * pi * rand(24, 6));

            % (1.0, Inf) -> B0
            y_b0 = apply_sft_p15_reconstruct(S1, carrier, "P15_4_dual_domain", ...
                struct("Axis", "range", "Tau", 1.0, "MinRunLength", Inf));
            t.verifyEqual(y_b0, S1);

            % (0.6, Inf) -> P15_1
            y_p1 = apply_sft_p15_reconstruct(S1, carrier, "P15_1_carrier_agree", ...
                struct("Axis", "range", "Tau", 0.6));
            y_4p1 = apply_sft_p15_reconstruct(S1, carrier, "P15_4_dual_domain", ...
                struct("Axis", "range", "Tau", 0.6, "MinRunLength", Inf));
            t.verifyEqual(y_4p1, y_p1);

            % (1.0, 3) -> P15_2 (w=1, lmin=3)
            y_p2 = apply_sft_p15_reconstruct(S1, [], "P15_2_rlp_she", ...
                struct("Axis", "range", "HalfWidth", 1, "MinRunLength", 3));
            y_4p2 = apply_sft_p15_reconstruct(S1, carrier, "P15_4_dual_domain", ...
                struct("Axis", "range", "Tau", 1.0, "MinRunLength", 3));
            t.verifyEqual(y_4p2, y_p2);

            % 并集契约：任意 (0.6, 3) 置零点 == P15_1 置零点 | P15_2 置零点
            y_4 = apply_sft_p15_reconstruct(S1, carrier, "P15_4_dual_domain", ...
                struct("Axis", "range", "Tau", 0.6, "MinRunLength", 3));
            t.verifyEqual(real(y_4) == 0, (real(y_p1) == 0) | (real(y_p2) == 0));
            t.verifyEqual(imag(y_4) == 0, (imag(y_p1) == 0) | (imag(y_p2) == 0));
        end

        function endToEndSmokeAndParameterOnlyResume(t)
            folder = TestSftP15Reconstruction.createFixtureDir("sft_p15_smoke");
            cleanup = onCleanup(@() TestSftP15Reconstruction.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_test.mat");
            save(param_file, "-struct", "S60");

            data_dir = fullfile(folder, "raw", "SceneA");
            mkdir(data_dir);
            rng(2026);
            raw = single(complex(randn(S60.nrn, 24), randn(S60.nrn, 24)));
            save(fullfile(data_dir, "rstart 1.mat"), "raw");

            cfg = zc1bit_smoke_config(fullfile(folder, "raw"), param_file);
            cfg.results_root = fullfile(folder, "results");
            cfg.dataset_names = "SceneA";
            cfg.enl_dataset_names = "SceneA";
            cfg.enl_window_size = 4;
            cfg.enl_stride = 2;
            cfg.num_samples_per_dataset = 4; % 2 calibration (1,3), 2 test (2,4)

            % 直注 B0 SFT 锁参表与精简网格（含 Inf/边界点）
            cfg.sft_p15.range_allocations = [1, 1];
            cfg.sft_p15.azimuth_allocations = [1, 2];
            cfg.sft_p15.p1_tau_grid = [0.6, 1.0];
            cfg.sft_p15.p2_w1_lmin_grid = [2, inf];
            cfg.sft_p15.p2_w2_lmin_grid = 5;
            cfg.sft_p15.p3_tau_grid = [0.5, 1.0];
            cfg.sft_p15.p4_tau_grid = [0.8, 1.0];
            cfg.sft_p15.p4_lmin_grid = [3, inf];
            cfg.sft_p15.phase_mismatch_deg = [0, 30, 90];
            cfg.sft_p15.phase_mismatch_units = ["rb_R1A1_range_rsft", "ab_R1A2_azimuth_asft"];
            cfg.sft_p15.leakage_groups = ["R1A1", "R1A2"];
            cfg.sft_p15.b0_locked_table = table( ...
                ["rb_R1A1_range_rsft"; "ab_R1A2_azimuth_asft"], ...
                ["B0_par_raw"; "B0_par_raw"], ...
                ["s0_fr0p3_fa0"; "sm2_fr0_fa0p6"], ...
                ["s0_fr0p3_fa0"; "sm2_fr0_fa0p6"], ...
                [0; -2], [0.3; 0], [0; 0.6], [0; 0], ...
                'VariableNames', {'UnitKey', 'VariantKey', 'ParamTag', 'ParamKey', 'STR_dB', 'FrOverBr', 'FaOverBa', 'Phi0'});

            out = run_zc1bit_sft_p15(cfg, "all");
            t.verifyEqual(height(out.search.locked_params), 2 * 6);
            t.verifyEqual(height(out.evaluate.summary_table), 2 * 6);
            t.verifyEqual(height(out.evaluate.wilcoxon_table), 2 * (5 + 4)); % 5 vs_B0 + 4 vs_R1
            t.verifyTrue(all(ismember(out.evaluate.detail_table.LocalSampleIdx, [2; 4])));
            t.verifyFalse(any(ismember(out.evaluate.detail_table.LocalSampleIdx, [1; 3])));

            % 验证载波基频抑制字段存在且有限
            t.verifyTrue(all(isfinite(out.evaluate.leakage_table.CarrierToneRatio)));
            t.verifyTrue(all(isfinite(out.evaluate.leakage_table.CarrierSuppression_dB)));

            % 验证全盲 P15_2_rlp_she 在全部相偏角度下指标严格一致
            rob = out.evaluate.robustness_table;
            p2_rob = rob(rob.UnitKey == "ab_R1A2_azimuth_asft" & rob.MethodKey == "P15_2_rlp_she", :);
            t.verifyEqual(height(p2_rob), 3);
            t.verifyLessThan(max(abs(p2_rob.SSIM_Mean - p2_rob.SSIM_Mean(1))), 1e-12);

            % 验证修改非计算参数与 NaN 覆写回波后 100% 断点命中零重算
            raw_nan = complex(nan(S60.nrn, 24), nan(S60.nrn, 24));
            save(fullfile(data_dir, "rstart 1.mat"), "raw_nan", "-v7.3");
            cfg.runtime.cache_gb = 2;
            out_resumed = run_zc1bit_sft_p15(cfg, "all");
            t.verifyEqual(out_resumed.evaluate.summary_table.SSIM_Mean, out.evaluate.summary_table.SSIM_Mean);
        end
    end

    methods (Static)
        function root = workspaceTestTmpRoot()
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            root = fullfile(repo, "tmp", "tests");
        end

        function folder = createFixtureDir(tag)
            root = TestSftP15Reconstruction.workspaceTestTmpRoot();
            if ~isfolder(root)
                mkdir(root);
            end
            [~, uniq] = fileparts(tempname(char(root)));
            folder = fullfile(root, string(tag) + "_" + string(uniq));
            mkdir(folder);
        end

        function removeFixture(folder)
            resolved = string(java.io.File(char(folder)).getCanonicalPath());
            root = string(java.io.File(char(TestSftP15Reconstruction.workspaceTestTmpRoot())).getCanonicalPath());
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
