function Export2Folder(figure_name,folder_name) %exports figure to articale figure folder
    set(gcf, 'Color', 'w');

    export_fig(figure_name)  %cool function - check internet for support

    movefile(figure_name, fullfile(folder_name,figure_name));

end    