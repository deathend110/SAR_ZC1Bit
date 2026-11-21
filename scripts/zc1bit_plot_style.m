function style = zc1bit_plot_style()
%ZC1BIT_PLOT_STYLE 绘图脚本统一视觉参数；不改变任何数据口径。

style = struct();
style.font_name = "Arial";
style.font_size = 11;
style.title_size = 12;
style.label_size = 12;
style.line_width = 1.8;
style.marker_size = 5;
style.baseline_color = [0, 0, 0];
style.baseline_style = "--";
style.pseudo_colors = lines(8);
style.figure_position = [100, 100, 980, 540];
style.export_resolution = 300;
style.sig_marker_size = 7;
end
