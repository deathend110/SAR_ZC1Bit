classdef TestSamplingAndThreshold < matlab.unittest.TestCase
    % 上采样/裁剪网格解析与 ZT/SplitRT/RSFT/ASFT/2D-SFT 阈值构建单元测试。

    methods (TestClassSetup)
        function addRepoPaths(t)
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(repo));
            for name = ["sampling", "threshold", "quantization", "compat"]
                t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo, "src", name)));
            end
        end
    end

    methods (Test)
        function samplingGridUpsampleAndCrop(t)
            S60 = struct("nrn", 12, "nan", 10, "Fs", 60e6, "prf", 240, ...
                "B", 45e6, "Bd", 200, "R0", 1e4, "C", 3e8, "tnrn", zeros(12, 1));
            g_frac = resolve_sampling_grid(sqrt(3)^2, 1.5, [12, 10], S60);
            t.verifyEqual(g_frac.q_range, 3);
            t.verifyEqual(g_frac.upsampled_size, [36, 15]);
            t.verifyEqual(g_frac.Fs_up, 3 * S60.Fs);
            t.verifyEqual(g_frac.PRF_up, 1.5 * S60.prf);

            rng(5);
            S = complex(randn(12, 10), randn(12, 10));
            S_up = two_dim_upsample_fft(S, 3, 2);
            t.verifySize(S_up, [36, 20]);
            S_down = two_dim_downsample_fft(S_up, 3, 2, S60);
            t.verifySize(S_down, size(S));

            % 整数倍率常量信号先插值再裁频谱与理论倍率一致
            const_sig = complex(ones(12, 10), zeros(12, 10));
            c_up = two_dim_upsample_fft(const_sig, 2, 3);
            c_down = two_dim_downsample_fft(c_up, 2, 3, S60);
            t.verifyEqual(c_down, const_sig * 6, "AbsTol", 1e-12);
        end

        function thresholdFamiliesNyquistAndPhaseSum(t)
            S60 = struct("nrn", 16, "nan", 12, "Fs", 60e6, "prf", 240, ...
                "B", 45e6, "Bd", 200, "R0", 1e4, "C", 3e8);
            rng(21);
            S = complex(randn(16, 12), randn(16, 12));
            S_up = two_dim_upsample_fft(S, 4, 2);
            grid = resolve_sampling_grid(4, 2, size(S), S60);

            % ZT 全零且尺寸一致
            U_zt = build_threshold(S_up, grid, struct("ThresholdType", "ZT"), S60);
            t.verifySize(U_zt, size(S_up));
            t.verifyEqual(U_zt, complex(zeros(size(S_up))));

            % SplitRT：本地 RandStream 严格可复现且幅值为 As * sigma_hat
            s1 = RandStream("mt19937ar", "Seed", 12345);
            s2 = RandStream("mt19937ar", "Seed", 12345);
            U_rt1 = build_threshold(S_up, grid, struct("ThresholdType", "RT", "As", 0.6), ...
                S60, struct("rand_stream", s1));
            U_rt2 = build_threshold(S_up, grid, struct("ThresholdType", "RT", "As", 0.6), ...
                S60, struct("rand_stream", s2));
            t.verifyEqual(U_rt1, U_rt2);
            expected_amp = 0.6 * sqrt(2 / pi) * mean(abs(S_up(:)));
            t.verifyEqual(abs(U_rt1), repmat(expected_amp, size(S_up)), "AbsTol", 1e-12);

            % RSFT / ASFT 作为 2D-SFT 单轴特例严格等价，且幅值满足 sigma_hat / 10^(STR/20)
            p_rsft = struct("ThresholdType", "RSFT", "STR_dB", -4, "fr_over_Br", 1.2, "fa_over_Ba", 0, "phi0", 0);
            p_2d_r = struct("ThresholdType", "2D-SFT", "STR_dB", -4, "fr_over_Br", 1.2, "fa_over_Ba", 0, "phi0", 0);
            [U_rsft, meta_r] = build_threshold(S_up, grid, p_rsft, S60);
            U_2d_r = build_threshold(S_up, grid, p_2d_r, S60);
            t.verifyEqual(U_rsft, U_2d_r, "AbsTol", 1e-12);
            t.verifyEqual(meta_r.Amplitude, estimate_sigma_hat(S_up) / 10^(-4 / 20), "AbsTol", 1e-12);

            p_asft = struct("ThresholdType", "ASFT", "STR_dB", 2, "fr_over_Br", 0, "fa_over_Ba", 0.8, "phi0", 0);
            p_2d_a = struct("ThresholdType", "2D-SFT", "STR_dB", 2, "fr_over_Br", 0, "fa_over_Ba", 0.8, "phi0", 0);
            t.verifyEqual(build_threshold(S_up, grid, p_asft, S60), ...
                build_threshold(S_up, grid, p_2d_a, S60), "AbsTol", 1e-12);

            % 2D-SFT 相位先相加再取复指数，幅值恒模
            p_2d = struct("ThresholdType", "2D-SFT", "STR_dB", 0, "fr_over_Br", 1.0, "fa_over_Ba", 0.5, "phi0", 0);
            [U_2d, meta_2d] = build_threshold(S_up, grid, p_2d, S60);
            t.verifyEqual(abs(U_2d), repmat(meta_2d.Amplitude, size(S_up)), "AbsTol", 1e-12);

            % Nyquist 越界必须抛出错误
            fr_lim = 0.98 * grid.Fs_up / (2 * S60.B);
            p_bad = struct("ThresholdType", "RSFT", "STR_dB", 0, "fr_over_Br", fr_lim, "fa_over_Ba", 0, "phi0", 0);
            t.verifyError(@() build_threshold(S_up, grid, p_bad, S60), "zc1bit:SFTNyquistViolation");
        end
    end
end
