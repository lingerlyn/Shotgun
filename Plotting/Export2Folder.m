function Export2Folder(figure_name,folder_name) %exports figure to articale figure folder
    set(gcf, 'Color', 'w');
%     set(findall(gcf,'type','text'),'fontWeight','bold')
    set(findall(gcf,'type','axe'),'box','off')  %remove all boxes
    export_fig(figure_name, '-rgb', '-r600')  %cool function - check internet for support

    movefile(figure_name, fullfile(folder_name,figure_name));

end    