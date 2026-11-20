classdef TestZc1BitWorkflow < matlab.unittest.TestCase
    % 路径守卫、35/35 清单防泄漏、仅实际计算参数断点恢复与阶段一冒烟闭环测试。

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
        function pathGuardIsolatesForeignDirectoryAndRestoresOnError(t)
            cfg = zc1bit_config();
            folder = string(tempname);
            mkdir(folder);
            cleanup = onCleanup(@() TestZc1BitWorkflow.removeFixture(folder));

            shadow = fullfile(folder, "Range_Compress.m");
            fid = fopen(shadow, "w", "n", "UTF-8");
            fprintf(fid, "function varargout = Range_Compress(varargin)\nerror('foreign:Shadow','外部冲突函数');\nend\n");
            fclose(fid);

            prev = pwd;
            dir_cleanup = onCleanup(@() cd(prev));
            cd(folder);
            t.verifyEqual(string(which("Range_Compress")), shadow);

            guard = setup_zc1bit_paths(cfg); %#ok<NASGU>
            t.verifyEqual(string(pwd), cfg.repo_root);
            t.verifyEqual(string(which("Range_Compress")), fullfile(cfg.repo_root, "src", "sarproc", "Range_Compress.m"));
            clear guard
            t.verifyEqual(string(pwd), folder);

            % 报错退出时也必须恢复原始目录
            bad_cfg = cfg;
            bad_cfg.stage = "unknown_stage";
            t.verifyError(@() run_zc1bit_workflow(bad_cfg), "zc1bit:UnknownStage");
            t.verifyEqual(string(pwd), folder);
            clear dir_cleanup
        end

        function manifestSplitAndParameterOnlyCheckpointResume(t)
            folder = string(tempname);
            mkdir(folder);
            cleanup = onCleanup(@() TestZc1BitWorkflow.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_test.mat");
            save(param_file, "-struct", "S60");

            data_dir = fullfile(folder, "raw", "SceneA");
            mkdir(data_dir);
            raw = sparse(1, 1, 1i, S60.nrn, 40);
            save(fullfile(data_dir, "rstart 1.mat"), "raw");
            save(fullfile(data_dir, "rstart 2.mat"), "raw");

            cfg = zc1bit_smoke_config(fullfile(folder, "raw"), param_file);
            cfg.results_root = fullfile(folder, "results");
            cfg.dataset_names = "SceneA";
            cfg.num_samples_per_dataset = 10;
            out1 = run_zc1bit_workflow(cfg);
            m = out1.record.manifest;
            t.verifyEqual(sum(m.Split == "calibration"), 5);
            t.verifyEqual(sum(m.Split == "test"), 5);
            t.verifyEqual(m.LocalSampleIdx(m.Split == "calibration"), [1; 3; 5; 7; 9]);
            t.verifyEqual(m.LocalSampleIdx(m.Split == "test"), [2; 4; 6; 8; 10]);

            % 改变 data_root、results_root、workers 等非计算参数，检查点仍应匹配通过
            moved_raw = fullfile(folder, "raw_moved");
            copyfile(fullfile(folder, "raw"), moved_raw);
            cfg2 = cfg;
            cfg2.data_root = moved_raw;
            cfg2.screen.workers = 4;
            cfg2.screen.use_fast_path = false;
            out2 = run_zc1bit_workflow(cfg2);
            t.verifyEqual(out2.record.manifest.CStart, m.CStart);
            t.verifyTrue(all(startsWith(out2.record.manifest.FilePath, moved_raw)));

            % 原子保存/读取检查点，修改实际计算参数则报错 ParameterMismatch
            cp_file = fullfile(folder, "cp.mat");
            saved_state = struct("parameters", zc1bit_parameters(cfg), "completed", true);
            zc1bit_atomic_save(cp_file, struct("state", saved_state));
            loaded_state = zc1bit_load_state(cp_file, struct());
            zc1bit_check_parameters(loaded_state.parameters, cfg2);

            bad_cfg = cfg;
            bad_cfg.screen.rt_as_values = [0.5, 0.7];
            t.verifyError(@() zc1bit_check_parameters(loaded_state.parameters, bad_cfg), "zc1bit:ParameterMismatch");
        end

        function stage1SmokeClosureAndTestIsolation(t)
            folder = string(tempname);
            mkdir(folder);
            cleanup = onCleanup(@() TestZc1BitWorkflow.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_smoke.mat");
            save(param_file, "-struct", "S60");

            raw_root = fullfile(folder, "raw");
            scene_dir = fullfile(raw_root, "SmokeScene");
            mkdir(scene_dir);
            rng(33);
            % 构造 24 列，其中后半段包含测试窗口；验证阶段一绝不依赖/运行测试样本
            raw = complex(randn(S60.nrn, 24), randn(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw", "-v7.3");

            cfg = zc1bit_smoke_config(raw_root, param_file);
            cfg.results_root = fullfile(folder, "results");
            cfg.stage = "screen";
            out = run_zc1bit_workflow(cfg);

            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_detail.csv")));
            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_summary.csv")));
            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_rule_ranking.csv")));
            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_cal_enl_roi_manifest.csv")));
            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_results.mat")));
            t.verifyTrue(isfile(fullfile(out.paths.screen, "screen_metadata.txt")));

            % 明细表仅包含 1 个校准样本（LocalSampleIdx=1），绝不含测试样本（LocalSampleIdx=2）
            detail = out.screen.detail_table;
            t.verifyTrue(all(detail.LocalSampleIdx == 1));
            t.verifyFalse(any(detail.LocalSampleIdx == 2));

            % 验证 R1A1 在 range_branch 与 azimuth_branch 的 Key 互不冲突
            [tgroups, ~] = zc1bit_build_screen_tasks(cfg, S60);
            t.verifyEqual(numel(unique(tgroups.ThresholdKey)), height(tgroups));

            % 用全 NaN 覆盖原始回波内容并修改 workers 后重跑 screen，必须 100% 命中断点零重算
            raw_nan = complex(nan(S60.nrn, 24), nan(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw_nan", "-v7.3");
            cfg.screen.workers = 2;
            out_resumed = run_zc1bit_screen(cfg);
            t.verifyEqual(out_resumed.summary_table.PSNR_Mean, out.screen.summary_table.PSNR_Mean);

            % 验证正式配置在本地 FS60_params.mat 下可无冲突构建全部阶段一任务与 Nyquist 检查
            real_cfg = zc1bit_config();
            if isfile(real_cfg.parameter_file_60)
                S60_real = load(real_cfg.parameter_file_60);
                [real_tg, real_vs] = zc1bit_build_screen_tasks(real_cfg, S60_real);
                t.verifyEqual(numel(unique(real_tg.ThresholdKey)), height(real_tg));
                t.verifyEqual(numel(unique(real_vs.VariantKey)), height(real_vs));
                t.verifyGreaterThan(height(real_tg), 0);
                t.verifyGreaterThan(height(real_vs), 0);
            end

            % 阶段二未启用前必须明确拒绝对外调用
            t.verifyError(@() run_zc1bit_search(cfg), "zc1bit:StageNotImplemented");
            t.verifyError(@() run_zc1bit_confirm(cfg), "zc1bit:StageNotImplemented");
            t.verifyError(@() run_zc1bit_evaluate(cfg), "zc1bit:StageNotImplemented");
        end
    end

    methods (Static)
        function removeFixture(folder)
            resolved = string(java.io.File(char(folder)).getCanonicalPath());
            root = string(java.io.File(tempdir).getCanonicalPath());
            if resolved == root || ~startsWith(lower(resolved), lower(root + string(filesep)))
                error("zc1bit:UnsafeCleanup", "拒绝清理临时目录以外的路径：%s", resolved);
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
