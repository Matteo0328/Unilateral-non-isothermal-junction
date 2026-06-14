function exportPDF(filesAndFigs, opts)
%EXPORTPDF Exports figures to PDF handling themes, margins, A4 scale, and folders.
%
%   exportPDF(filesAndFigs, Name, Value) exports one or multiple figures 
%   to PDF files, with advanced control over layout, scaling, and themes.
%
%   Inputs:
%       filesAndFigs - Can be a single filename (uses current figure), or
%                      a comma-separated list of ('filename', figure_handle) pairs.
%
%   Name-Value Options:
%       Folder       - Target subfolder for saved PDFs [string/char]. Default: ''
%       Margin       - Border thickness in centimeters [double]. Default: 0
%       A4Scale      - Scales the figure to fit an A4 sheet [logical]. Default: false
%       Orientation  - 'auto', 'landscape', or 'portrait'. Default: 'auto'
%       KeepAspect   - Maintains screen aspect ratio on A4 [logical]. Default: true
%
%   Examples:
%       exportPDF('my_plot')
%       exportPDF('file1', fig1, 'file2', fig2)
%       exportPDF('my_plot', gcf, Folder="Results/Report")
%       exportPDF('my_plot', gcf, Margin=1.5)
%       exportPDF('my_plot', gcf, A4Scale=true, Orientation="landscape")

    % Function argument validation and defaults
    arguments (Repeating)
        filesAndFigs 
    end
    
    arguments
        opts.Folder (1,:) char = '' 
        opts.Margin (1,1) double = 0
        opts.A4Scale (1,1) logical = false
        opts.Orientation (1,:) char {mustBeMember(opts.Orientation, {'auto', 'landscape', 'portrait'})} = 'auto'
        opts.KeepAspect (1,1) logical = true
    end
    
    margin_cm = opts.Margin;
    a4_scale = opts.A4Scale;
    orientation = opts.Orientation;
    keep_aspect = opts.KeepAspect;
    out_folder = opts.Folder;
    
    % Handle destination folder creation
    if ~isempty(out_folder)
        if ~exist(out_folder, 'dir')
            [status, msg] = mkdir(out_folder);
            if ~status
                error('Cannot create folder "%s". Reason: %s', out_folder, msg);
            end
            fprintf(' -> Created new folder: %s\n', out_folder);
        end
    end
    
    % Parse inputs into filenames and figure handles
    fileNames = {};
    figHandles = gobjects(0);
    
    if length(filesAndFigs) == 1 && (ischar(filesAndFigs{1}) || isstring(filesAndFigs{1}))
        fileNames = {char(filesAndFigs{1})};
        figHandles = gcf;
    else
        for i = 1:2:length(filesAndFigs)-1
            if (ischar(filesAndFigs{i}) || isstring(filesAndFigs{i})) && isgraphics(filesAndFigs{i+1})
                fileNames{end+1} = char(filesAndFigs{i});
                figHandles(end+1) = filesAndFigs{i+1};
            else
                error('Arguments must be passed as pairs: ''filename'', figure_handle.');
            end
        end
    end
    
    if isempty(fileNames)
        error('No valid files or figures found to export.');
    end
    
    % Process and export each figure
    for i = 1:length(figHandles)
        f = figHandles(i);
        fname = fileNames{i};
        
        if ~endsWith(fname, '.pdf', 'IgnoreCase', true)
            fname = [fname, '.pdf'];
        end
        
        if ~isempty(out_folder)
            fname = fullfile(out_folder, fname);
        end
        
        try
            % Temporarily force light theme for clean PDF export
            hasTheme = isprop(f, 'Theme');
            if hasTheme
                oldTheme = f.Theme;
                f.Theme = 'light';
                drawnow; 
            else
                oldColor = f.Color;
                oldInvert = f.InvertHardcopy;
                f.Color = 'w';
                f.InvertHardcopy = 'on';
            end
            
            % Export logic
            if margin_cm == 0 && ~a4_scale
                % Standard tight export
                if isa(f, 'matlab.ui.Figure')
                    set(f, 'PaperPositionMode', 'auto');
                end
                exportgraphics(f, fname, 'ContentType', 'vector', 'BackgroundColor', 'white');
            else
                % Custom layout export (Margins or A4)
                if ~isa(f, 'matlab.ui.Figure')
                    error('Layout management requires a Figure object.');
                end
                
                oldUnits = f.Units;
                oldPaperUnits = f.PaperUnits;
                oldPaperPosMode = f.PaperPositionMode;
                
                f.Units = 'centimeters';
                f.PaperUnits = 'centimeters';
                
                figPos = f.Position;
                W_fig = figPos(3);
                H_fig = figPos(4);
                
                if a4_scale
                    if strcmpi(orientation, 'landscape')
                        target_W = 29.7; target_H = 21.0;
                    elseif strcmpi(orientation, 'portrait')
                        target_W = 21.0; target_H = 29.7;
                    else
                        if W_fig > H_fig
                            target_W = 29.7; target_H = 21.0;
                        else
                            target_W = 21.0; target_H = 29.7;
                        end
                    end
                    
                    avail_W = target_W - 2*margin_cm;
                    avail_H = target_H - 2*margin_cm;
                    
                    if avail_W <= 0 || avail_H <= 0
                        error('Margin is too large for an A4 sheet.');
                    end
                    
                    if keep_aspect
                        scale_factor = min(avail_W / W_fig, avail_H / H_fig);
                        W_plot = W_fig * scale_factor;
                        H_plot = H_fig * scale_factor;
                    else
                        W_plot = avail_W;
                        H_plot = avail_H;
                    end
                    
                    pos_x = (target_W - W_plot) / 2;
                    pos_y = (target_H - H_plot) / 2;
                    
                    f.PaperSize = [target_W, target_H];
                    f.PaperPosition = [pos_x, pos_y, W_plot, H_plot];
                else
                    f.PaperSize = [W_fig + 2*margin_cm, H_fig + 2*margin_cm];
                    f.PaperPosition = [margin_cm, margin_cm, W_fig, H_fig];
                end
                
                f.PaperPositionMode = 'manual';
                print(f, fname, '-dpdf', '-vector');
                
                % Restore original dimension parameters
                f.Units = oldUnits;
                f.PaperUnits = oldPaperUnits;
                f.PaperPositionMode = oldPaperPosMode;
            end
            
            % Restore original visual theme
            if hasTheme
                f.Theme = oldTheme;
                drawnow;
            else
                f.Color = oldColor;
                f.InvertHardcopy = oldInvert;
            end
            
            % Print success status to terminal
            if a4_scale
                orient_str = 'Auto';
                if strcmpi(orientation, 'landscape'), orient_str = 'Landscape'; end
                if strcmpi(orientation, 'portrait'), orient_str = 'Portrait'; end
                
                aspect_str = 'Proportional';
                if ~keep_aspect, aspect_str = 'Stretched'; end
                
                fprintf(' -> Saved on A4 (%s, %s): %s\n', orient_str, aspect_str, fname);
            else
                fprintf(' -> Saved: %s\n', fname);
            end
            
        catch ME
            % Safety catch to ensure theme is restored even if export fails
            if exist('hasTheme','var') && hasTheme && exist('oldTheme','var')
                f.Theme = oldTheme;
            elseif exist('oldColor','var')
                f.Color = oldColor;
                f.InvertHardcopy = oldInvert;
            end
            warning('Failed to save %s. Reason: %s', fname, ME.message);
        end
    end
    
    fprintf('Export completed.\n');
end