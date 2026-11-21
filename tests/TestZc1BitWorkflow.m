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
            folder = TestZc1BitWorkflow.createFixtureDir("pathguard");
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
            folder = TestZc1BitWorkflow.createFixtureDir("manifest");
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

            base_cfg = zc1bit_config();
            t.verifyEqual(base_cfg.tmp_root, fullfile(base_cfg.repo_root, "tmp"));
            t.verifyTrue(startsWith(lower(folder), lower(fullfile(base_cfg.repo_root, "tmp", "tests") + string(filesep))));
            t.verifyEqual(out1.paths.tmp, base_cfg.tmp_root);
            t.verifyError(@() TestZc1BitWorkflow.removeFixture(base_cfg.repo_root), "zc1bit:UnsafeCleanup");
            t.verifyError(@() TestZc1BitWorkflow.removeFixture(TestZc1BitWorkflow.workspaceTestTmpRoot()), "zc1bit:UnsafeCleanup");

            % 改变 data_root、results_root、tmp_root、runtime 内存参数等非计算参数，检查点仍应匹配通过
            moved_raw = fullfile(folder, "raw_moved");
            copyfile(fullfile(folder, "raw"), moved_raw);
            cfg2 = cfg;
            cfg2.data_root = moved_raw;
            cfg2.tmp_root = fullfile(folder, "tmp_moved");
            cfg2.runtime.cache_gb = 2;
            cfg2.runtime.preload_signals = false;
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
            folder = TestZc1BitWorkflow.createFixtureDir("smoke");
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

            % 验证扁平断点目录：无子目录，且包含 gt_cal_cache.mat 与全部 ThresholdKey .mat
            cp_screen_dir = fullfile(out.paths.checkpoints, "screen");
            cp_entries = dir(cp_screen_dir);
            cp_names = string({cp_entries.name});
            cp_is_subdir = [cp_entries.isdir] & ~ismember(cp_names, [".", ".."]);
            t.verifyFalse(any(cp_is_subdir));
            expected_cp_files = ["gt_cal_cache.mat"; tgroups.ThresholdKey + ".mat"];
            t.verifyTrue(all(ismember(expected_cp_files, cp_names)));

            % 用全 NaN 覆盖原始回波内容并修改 runtime 内存参数后重跑 screen，必须 100% 命中断点零重算
            raw_nan = complex(nan(S60.nrn, 24), nan(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw_nan", "-v7.3");
            cfg.runtime.cache_gb = 2;
            cfg.runtime.preload_signals = false;
            out_resumed = run_zc1bit_screen(cfg);
            t.verifyEqual(out_resumed.summary_table.PSNR_Mean, out.screen.summary_table.PSNR_Mean);

            % 验证正式配置在本地 FS60_params.mat 下可无冲突构建全部 36 个单向四阈值阶段一任务
            real_cfg = zc1bit_config();
            if isfile(real_cfg.parameter_file_60)
                S60_real = load(real_cfg.parameter_file_60);
                [real_tg, real_vs] = zc1bit_build_screen_tasks(real_cfg, S60_real);
                t.verifyEqual(numel(unique(real_tg.ThresholdKey)), height(real_tg));
                t.verifyEqual(numel(unique(real_vs.VariantKey)), height(real_vs));
                t.verifyEqual(height(real_tg), 36);
                t.verifyFalse(any(real_tg.ThresholdType == "2D-SFT"));
                t.verifyTrue(all(real_tg.Range_q == 1 | real_tg.Azimuth_q == 1));
                t.verifyGreaterThan(height(real_vs), 0);
            end

            % 未完成前序阶段时必须明确拒绝跳过锁参直接 confirm / evaluate
            t.verifyError(@() run_zc1bit_confirm(cfg), "zc1bit:MissingSearchResults");
            t.verifyError(@() run_zc1bit_evaluate(cfg), "zc1bit:MissingLockedConfirm");
        end

        function cacheModeMatchesStreamingFallback(t)
            folder = TestZc1BitWorkflow.createFixtureDir("cache_vs_stream");
            cleanup = onCleanup(@() TestZc1BitWorkflow.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_smoke.mat");
            save(param_file, "-struct", "S60");

            raw_root = fullfile(folder, "raw");
            scene_dir = fullfile(raw_root, "SmokeScene");
            mkdir(scene_dir);
            rng(33);
            raw = complex(randn(S60.nrn, 24), randn(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw", "-v7.3");

            cfg_cached = zc1bit_smoke_config(raw_root, param_file);
            cfg_cached.results_root = fullfile(folder, "results_cached");
            cfg_cached.stage = "screen";
            cfg_cached.runtime.cache_gb = 1;
            cfg_cached.runtime.preload_signals = true;
            out_cached = run_zc1bit_workflow(cfg_cached);

            cfg_stream = zc1bit_smoke_config(raw_root, param_file);
            cfg_stream.results_root = fullfile(folder, "results_stream");
            cfg_stream.stage = "screen";
            cfg_stream.runtime.cache_gb = 1e-9;
            cfg_stream.runtime.preload_signals = false;
            out_stream = run_zc1bit_workflow(cfg_stream);

            s1 = out_cached.screen.summary_table;
            s2 = out_stream.screen.summary_table;
            t.verifyEqual(s1.PSNR_Mean, s2.PSNR_Mean);
            t.verifyEqual(s1.SSIM_Mean, s2.SSIM_Mean);
            t.verifyEqual(s1.ZeroRatio_Mean, s2.ZeroRatio_Mean);
            t.verifyEqual(s1.ENL_Mean, s2.ENL_Mean);

            % 验证 zc1bit_atomic_save 默认 -v7 往返读写一致性
            v7_file = fullfile(folder, "atomic_v7_roundtrip.mat");
            zc1bit_atomic_save(v7_file, struct("summary_table", s1, "tag", "v7_ok"));
            loaded = load(v7_file);
            t.verifyEqual(loaded.summary_table, s1);
            t.verifyEqual(loaded.tag, "v7_ok");
        end

        function stage2AutoTop5SearchConfirmEvaluateClosure(t)
            folder = TestZc1BitWorkflow.createFixtureDir("stage2_closure");
            cleanup = onCleanup(@() TestZc1BitWorkflow.removeFixture(folder));

            S60 = smallImagingParameters();
            param_file = fullfile(folder, "S60_smoke.mat");
            save(param_file, "-struct", "S60");

            raw_root = fullfile(folder, "raw");
            scene_dir = fullfile(raw_root, "SmokeScene");
            mkdir(scene_dir);
            rng(33);
            raw = complex(randn(S60.nrn, 24), randn(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw", "-v7.3");

            cfg = zc1bit_smoke_config(raw_root, param_file);
            cfg.results_root = fullfile(folder, "results");
            cfg.screen.rules = ["B0", "R1", "R3", "R4", "R7", "L1", "R9"];
            cfg.stage2.top_k = 5;
            cfg.stage2.sft_refine = true;

            % 未运行阶段一时，自动 Top 5 解析与阶段二执行顺序守卫必须报错
            t.verifyError(@() zc1bit_resolve_stage2_variants(cfg), "zc1bit:MissingScreenRanking");
            t.verifyError(@() run_zc1bit_confirm(cfg), "zc1bit:MissingSearchResults");
            t.verifyError(@() run_zc1bit_evaluate(cfg), "zc1bit:MissingLockedConfirm");

            % 先运行阶段一生成新版 screen_rule_ranking.csv（含过采样指标列）
            cfg.stage = "screen";
            out_s1 = run_zc1bit_workflow(cfg);
            rank_tbl = out_s1.screen.ranking_table;
            t.verifyTrue(all(ismember(["UpsampledConfigCount", "UpsampledMeanSSIM", "UpsampledMeanPSNR"], ...
                string(rank_tbl.Properties.VariableNames))));

            % 验证优化后自动读取 Top 5 非 B0 变体 + B0_par_raw（共 6 变体）及去冗余/去重开关
            sel_all = zc1bit_resolve_stage2_variants(cfg, out_s1.paths);
            t.verifyEqual(height(sel_all), 6);
            t.verifyEqual(sel_all.VariantKey(1), "B0_par_raw");
            t.verifyEqual(sel_all.VariantRole, ["baseline"; "top1"; "top2"; "top3"; "top4"; "top5"]);
            t.verifyFalse(any(sel_all.RuleKey(2:6) == "B0"));

            cfg_nodup = cfg;
            cfg_nodup.stage2.dedup_equivalent = false;
            sel_nodup = zc1bit_resolve_stage2_variants(cfg_nodup, out_s1.paths);
            expected_non_b0 = rank_tbl.VariantKey(rank_tbl.RuleKey ~= "B0");
            t.verifyEqual(sel_nodup.VariantKey(2:6), expected_non_b0(1:5));

            cfg_distinct = cfg;
            cfg_distinct.stage2.distinct_rules_only = true;
            cfg_distinct.stage2.top_k = 3;
            sel_distinct = zc1bit_resolve_stage2_variants(cfg_distinct, out_s1.paths);
            cand_rules = sel_distinct.RuleKey(2:end);
            t.verifyEqual(numel(unique(cand_rules)), numel(cand_rules));

            % 串行执行阶段二 search -> confirm -> evaluate
            out_search = run_zc1bit_search(cfg);
            t.verifyTrue(isfile(fullfile(out_search.output_dir, "search_selected_variants.csv")));
            t.verifyTrue(isfile(fullfile(out_search.output_dir, "search_summary.csv")));
            t.verifyTrue(isfile(fullfile(out_search.output_dir, "search_best_params.csv")));
            t.verifyTrue(isfile(fullfile(out_search.output_dir, "search_results.mat")));
            t.verifyTrue(isfile(fullfile(out_search.output_dir, "search_metadata.txt")));
            t.verifyTrue(all(out_search.cal_manifest.Split == "calibration"));
            t.verifyEqual(height(out_search.best_params_table), height(out_search.units) * 6);

            out_confirm = run_zc1bit_confirm(cfg);
            t.verifyTrue(isfile(fullfile(out_confirm.output_dir, "locked_params.csv")));
            t.verifyTrue(isfile(fullfile(out_confirm.output_dir, "locked_confirm_plan.csv")));
            t.verifyTrue(isfile(fullfile(out_confirm.output_dir, "confirm_results.mat")));
            t.verifyTrue(isfile(fullfile(out_confirm.output_dir, "confirm_metadata.txt")));
            t.verifyTrue(any(out_confirm.locked_confirm_plan.Track == "independent_recalibrated"));
            t.verifyTrue(any(out_confirm.locked_confirm_plan.Track == "b0_transferred"));

            out_eval = run_zc1bit_evaluate(cfg);
            t.verifyTrue(isfile(fullfile(out_eval.output_dir, "evaluate_detail.csv")));
            t.verifyTrue(isfile(fullfile(out_eval.output_dir, "evaluate_summary.csv")));
            t.verifyTrue(isfile(fullfile(out_eval.output_dir, "evaluate_leakage_summary.csv")));
            t.verifyTrue(isfile(fullfile(out_eval.output_dir, "evaluate_test_enl_roi_manifest.csv")));
            t.verifyTrue(isfile(fullfile(out_eval.output_dir, "evaluate_results.mat")));
            t.verifyTrue(isfile(fullfile(out_eval.stats_dir, "paired_wilcoxon_vs_b0.csv")));

            % 测试集严格隔离：评估明细只含测试样本（LocalSampleIdx == 2），绝不含校准样本
            t.verifyTrue(all(out_eval.detail_table.LocalSampleIdx == 2));
            t.verifyFalse(any(out_eval.detail_table.LocalSampleIdx == 1));
            t.verifyGreaterThan(height(out_eval.wilcoxon_table), 0);
            t.verifyTrue(all(isfinite(out_eval.wilcoxon_table.DeltaSSIM_Mean)));
            t.verifyGreaterThan(height(out_eval.leakage_table), 0);
            t.verifyTrue(all(out_eval.leakage_table.OffSupportRatio >= 0 & out_eval.leakage_table.OffSupportRatio <= 1));

            % 验证扁平断点结构（无子目录）与全 NaN 覆盖回波后的 100% 断点命中零重算
            for cp_sub = ["search", "evaluate"]
                cp_entries = dir(fullfile(out_s1.paths.checkpoints, cp_sub));
                cp_names = string({cp_entries.name});
                t.verifyFalse(any([cp_entries.isdir] & ~ismember(cp_names, [".", ".."])));
            end
            raw_nan = complex(nan(S60.nrn, 24), nan(S60.nrn, 24));
            save(fullfile(scene_dir, "rstart 1.mat"), "raw_nan", "-v7.3");
            cfg.runtime.cache_gb = 2;
            out_search_resumed = run_zc1bit_search(cfg);
            out_eval_resumed = run_zc1bit_evaluate(cfg);
            t.verifyEqual(out_search_resumed.best_params_table.Calibration_SSIM_Mean, out_search.best_params_table.Calibration_SSIM_Mean);
            t.verifyEqual(out_eval_resumed.summary_table.SSIM_Mean, out_eval.summary_table.SSIM_Mean);

            % 验证正式配置在本地 FS60_params.mat 下构建全部 54 个阶段二单元（无 BARU、无 2D-SFT）
            real_cfg = zc1bit_config();
            if isfile(real_cfg.parameter_file_60)
                S60_real = load(real_cfg.parameter_file_60);
                real_units = zc1bit_build_stage2_units(real_cfg, S60_real);
                t.verifyEqual(height(real_units), 54);
                t.verifyEqual(numel(unique(real_units.UnitKey)), 54);
                t.verifyFalse(any(real_units.ThresholdType == "2D-SFT"));
                t.verifyTrue(all(real_units.Range_q == 1 | real_units.Azimuth_q == 1));
            end
        end
    end

    methods (Static)
        function root = workspaceTestTmpRoot()
            repo = string(fileparts(fileparts(mfilename("fullpath"))));
            root = fullfile(repo, "tmp", "tests");
        end

        function folder = createFixtureDir(tag)
            root = TestZc1BitWorkflow.workspaceTestTmpRoot();
            if ~isfolder(root)
                mkdir(root);
            end
            [~, uniq] = fileparts(tempname(char(root)));
            folder = fullfile(root, string(tag) + "_" + string(uniq));
            mkdir(folder);
        end

        function removeFixture(folder)
            resolved = string(java.io.File(char(folder)).getCanonicalPath());
            root = string(java.io.File(char(TestZc1BitWorkflow.workspaceTestTmpRoot())).getCanonicalPath());
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
