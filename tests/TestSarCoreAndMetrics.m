classdef TestSarCoreAndMetrics < matlab.unittest.TestCase
    % 快路径成像等价性、ENL/ROI/Entropy/泄漏比、Wilcoxon 检验与加噪函数单元测试。

    methods (TestClassSetup)
        function addRepoPaths(t)
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(repo));
            for name = ["sarproc", "sampling", "threshold", "quantization", "metrics", "compat"]
                t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo, "src", name)));
            end
        end
    end

    methods (Test)
        function fastFocusMatchesStandardPathAcrossAllAllocations(t)
            S60 = smallImagingParameters();
            ctx = build_fast_focus_context(S60, [1, 2, 3], [1, 2, 3]);
            rng(41);
            S = complex(randn(S60.nrn, S60.nan), randn(S60.nrn, S60.nan));

            img_gt_std = generate_gt_image(S, S60);
            img_gt_fast = fast_focus_complex_roi(S, ctx, 1, 1);
            t.verifyEqual(img_gt_fast, img_gt_std, "AbsTol", 1e-10);

            allocs = [1, 1; 2, 1; 1, 2; 2, 2; 3, 3];
            for k = 1:size(allocs, 1)
                rq = allocs(k, 1);
                aq = allocs(k, 2);
                S_up = two_dim_upsample_fft(S, rq, aq);
                rs = RandStream("mt19937ar", "Seed", 100 + k);
                U = build_splitrt_threshold(S_up, 0.6, rs);
                c1 = quantize_1bit_with_u(S_up, U);
                cp15 = apply_zc_pseudo15bit(c1, "R1", struct("Axis", "range"));

                img_std = focus_upsampled_channel(cp15, S60, rq, aq);
                img_fast = fast_focus_complex_roi(cp15, ctx, rq, aq);
                t.verifyEqual(img_fast, img_std, "AbsTol", 1e-10);
            end
        end

        function metricsEnlEntropyWilcoxonAndGaussian(t)
            % ENL 数学公式与常量图 Inf
            intensity = [1, 2; 3, 4];
            img = sqrt(intensity);
            expected_enl = mean(intensity(:))^2 / var(intensity(:), 0);
            t.verifyEqual(compute_enl(img, 1, 1, 2), expected_enl, "AbsTol", 1e-12);
            t.verifyEqual(compute_enl(ones(4), 1, 1, 4), Inf);

            % Entropy 256-bin
            t.verifyEqual(compute_entropy(zeros(8), 256), 0, "AbsTol", 1e-12);
            t.verifyEqual(compute_entropy([0, 1; 0, 1], 256), 1, "AbsTol", 1e-12);

            % 固定均匀区 ROI 选择确定性与 build_enl_regions 过滤
            rng(7);
            g1 = 0.5 + 0.1 * randn(96);
            g1(33:64, 33:64) = 0.4 + 0.001 * randn(32);
            manifest = table((1:2).', ["SAR_Dataset_filed"; "SAR_Dataset_city1_histeq"], ...
                'VariableNames', {'SampleID', 'Dataset'});
            cfg = struct("enl_dataset_names", ["SAR_Dataset_filed", "SAR_Dataset_port"], ...
                "enl_window_size", 32, "enl_stride", 8);
            [is_enl, tops, lefts, ~, roi_tbl] = build_enl_regions({g1; g1}, manifest, cfg);
            t.verifyEqual(is_enl, [true; false]);
            t.verifyEqual(height(roi_tbl), 1);
            [top0, left0] = select_uniform_roi(g1, 32, 8);
            t.verifyEqual([tops(1), lefts(1)], [top0, left0]);

            % 配对 Wilcoxon（行列向量一致 & 零差值返回 p=1）
            a = [1.2, 2.5, 3.1, 4.8, 5.2, 6.0];
            b = [1.0, 2.1, 2.9, 4.0, 5.0, 5.5];
            w_row = compute_paired_wilcoxon(a, b);
            w_col = compute_paired_wilcoxon(a.', b.');
            t.verifyEqual(w_row, w_col);
            t.verifyGreaterThan(w_row.MeanDiff, 0);
            t.verifyTrue(isfinite(w_row.PValue));
            w_zero = compute_paired_wilcoxon(a, a);
            t.verifyEqual(w_zero.MeanDiff, 0);
            t.verifyEqual(w_zero.PValue, 1);

            % 频谱带外泄漏比非负且有界
            ref = randn(16, 16) + 1i * randn(16, 16);
            [off_r, rr, ar] = compute_leakage_metrics(ref, ref, 0.35);
            t.verifyGreaterThanOrEqual(off_r, 0);
            t.verifyLessThanOrEqual(off_r, 1);
            t.verifyGreaterThanOrEqual(rr, 0);
            t.verifyGreaterThanOrEqual(ar, 0);

            % gaussian SNR 口径与零均值
            rng(23);
            sig = complex(randn(400, 300), randn(400, 300));
            [noise, gstats] = gaussian(sig, 5, false);
            t.verifySize(noise, size(sig));
            t.verifyLessThan(abs(gstats.ActualSNRdB - 5), 0.05);
            t.verifyLessThan(abs(gstats.NoiseMean), 1e-10);
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
