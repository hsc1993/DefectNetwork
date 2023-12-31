clear all;close all;clc


case_str = 'bes';  
case_str = 'orthogonal';  
case_str = 'orthogonal_merged';  


flag_figure = 1;
% flag_figure = 0;

number_seeds_list = [20];
cellsize_matrix_list = [];
gnd_matrix_list = [];


for i = 1:size(number_seeds_list,2)
number_seeds = number_seeds_list(i);
n_trial = 1;



%% configuration
x = 300;
y = 300;
z = 300;
volume = x*y*z;

gnd_matrix = [];
gnd_old_matrix = [];
cellsize_matrix = [];


for ii = 1:n_trial
    dislocationLength = 0;
    %% general data import
    fid_general = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_config',num2str(ii),'.dat'),'r');
    formatSpec = '%s';
    general_info = textscan(fid_general,'%s', 'Delimiter', '\n');
    fclose(fid_general);
    
    %% data import
    % voro volume data import 
    fid_vorovolume = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_voro_volume',num2str(ii),'.dat'),'r');
    formatSpec = '%f';
    voro_volume = fscanf(fid_vorovolume,formatSpec);
    fclose(fid_vorovolume);
    
    cellsize = 0;
    for i = 1:size(voro_volume)
        cellsize = cellsize+voro_volume(i)^(1/3);
    end
    cellsize = cellsize/size(voro_volume,1);
    cellsize_matrix = [cellsize_matrix,cellsize];

    % voro gnd data import 
    fid_vorognd = fopen(strcat(case_str,num2str(number_seeds),'//', num2str(number_seeds),'_plot_voro_gnd',num2str(ii),'.dat'),'r');
    formatSpec = '%f';
    voro_gnd = fscanf(fid_vorognd,formatSpec);
    fclose(fid_vorognd);
    
    %% data processing
    voro_gnd_pervol = zeros(1,size(voro_gnd,1));
    total_l1 = 0;
    num_nonzero_gnd = 0;
    voro_volume_nonzero = 0;
    for i = 1:size(voro_gnd)
        total_l1 = total_l1 + voro_gnd(i); % A
        voro_gnd_pervol(i) = voro_gnd(i)/voro_volume(i);
        if voro_gnd(i) ~= 0
            num_nonzero_gnd = num_nonzero_gnd+1;
            voro_volume_nonzero = voro_volume_nonzero+voro_volume(i);
        end
    end
    
    gnd_among_nonzerovolume = total_l1/voro_volume_nonzero*1e20;
    gnd = total_l1/volume*1E20;  % m^-2 
    gnd_matrix = [gnd_matrix,gnd_among_nonzerovolume];
    gnd_old_matrix = [gnd_old_matrix,gnd];

    % dislocation start data import 
    formatSpec = '%f %f %f';
    fid_seg_start = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_dislocation_seg_start',num2str(ii),'.dat'),'rt');
    line = fgetl(fid_seg_start);
    seg_start_array = [str2num(line)];

    while ischar(line)
        line = fgetl(fid_seg_start);
        if line==-1
            break;
        end
        line_num = str2num(line);
        seg_start_array = [seg_start_array; line_num];
    end
    fclose(fid_seg_start);


    % dislocation end data import 
    formatSpec = '%f %f %f';
    fid_seg_end = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_dislocation_seg_end',num2str(ii),'.dat'),'rt');
    line = fgetl(fid_seg_end);
    seg_end_array = [str2num(line)];
    
    while ischar(line)
        line = fgetl(fid_seg_end);
        if line==-1
            break;
        end
        line_num = str2num(line);
        seg_end_array = [seg_end_array; line_num];
    end
    fclose(fid_seg_end);
    for segcount = 1:size(seg_start_array,1)
        dislocationLength = dislocationLength+norm(seg_end_array(segcount,:)-seg_start_array(segcount,:));
    end

    % dislocation burger data import 
    formatSpec = '%f %f %f';
    fid_seg_end = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_dislocation_seg_burger',num2str(ii),'.dat'),'rt');
    line = fgetl(fid_seg_end);
    burger_array = [str2num(line)];
    
    while ischar(line)
        line = fgetl(fid_seg_end);
        if line==-1
            break;
        end
        line_num = str2num(line);
        burger_array = [burger_array; line_num];
    end
    fclose(fid_seg_end);

    % calculate length total length for each burgers type
    length_burger_array = [];
    for i = 1:size(burger_array,1)
        length_burger_array = [length_burger_array;[burger_array(i,:) norm(seg_end_array(i,:)-seg_start_array(i,:))]];
    end

    % collapse lengths into a stand alone matrix
    collapse_length_burger_array = [[1.4300   1.4300    -1.4300 0];
                                    [-1.4300   -1.4300    -1.4300 0];
                                    [1.4300   -1.4300    1.4300 0];
                                    [-1.4300   1.4300    1.4300 0];
                                    ];
    
    for i = 1:size(length_burger_array,1)
        currentburger = length_burger_array(i,1:3);
        currentlength = length_burger_array(i,4);

            
        for j = 1:size(collapse_length_burger_array,1)
            if currentburger(1) == collapse_length_burger_array(j,1) && currentburger(2) == collapse_length_burger_array(j,2) && currentburger(3) == collapse_length_burger_array(j,3)
                collapse_length_burger_array(j,4) = collapse_length_burger_array(j,4)+currentlength;
            end
        end
 
    end
    collapse_length_burger_array(:,4)
        
    
    % dislocation belong data import 
    formatSpec = '%f %f %f';
    fid_seg_belong = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_dislocation_seg_cell_belong',num2str(ii),'.dat'),'rt');
    line = fgetl(fid_seg_belong);
    seg_cell_belong_array = [str2num(line)];
    
    while ischar(line)
        line = fgetl(fid_seg_belong);
        if line==-1
            break;
        end
        line_num = str2num(line);
        seg_cell_belong_array = [seg_cell_belong_array; line_num];
    end
    fclose(fid_seg_belong);


    % dislocation vorvx data import 
    formatSpec = '%f %f %f';
    fid_vorvx_scaled = fopen(strcat(case_str,num2str(number_seeds),'//',num2str(number_seeds),'_plot_vorvx_scaled',num2str(ii),'.dat'),'rt');
    line = fgetl(fid_vorvx_scaled);

    vorvx_scaled_array = [];
    vorvx_scaled_array_temp = [];
    while ischar(line)
        line = fgetl(fid_vorvx_scaled);
        if line==-1
            break;
        end
        line_num = str2num(line);
        vorvx_scaled_array_temp = [vorvx_scaled_array_temp; line_num];

        if strcmp(line,'new_seg')
            vorvx_scaled_array_temp_cell = {vorvx_scaled_array_temp};
            vorvx_scaled_array = [vorvx_scaled_array,vorvx_scaled_array_temp_cell];
            vorvx_scaled_array_temp = [];
        end
    end
    fclose(fid_vorvx_scaled);


    
    %% PLOT
    if flag_figure == 1
        for i = 1:size(vorvx_scaled_array,2)
            col(i,:) = [voro_gnd_pervol(i),voro_gnd_pervol(i),voro_gnd_pervol(i)];    
        end

        figure('position',[0 0 600 600],'Color',[1 1 1]);
        for i = 1:size(vorvx_scaled_array,2)
            K = convhulln(vorvx_scaled_array{i});
            if col(i,1) == 0
                continue
            end
            C_voro = voro_gnd_pervol(i)*1e20;
            trisurf(K,vorvx_scaled_array{i}(:,1),vorvx_scaled_array{i}(:,2),vorvx_scaled_array{i}(:,3),C_voro,'FaceAlpha',0.2,'EdgeAlpha',0)
            hold on;
        end

        for i = 1:size(seg_start_array(:,1))
            rng(seg_cell_belong_array(i));
%             C = rand(3,1);
            b = burger_array(i,:);
            if norm((norm(b)-norm(b(1))-norm(b(2))-norm(b(3))))<1e-3
                C = 'b';
            else
                C = 'r';
            end
            if i < size(seg_start_array(:,1),1)
                seg_x = [seg_start_array(i,1),seg_end_array(i,1)];
                seg_y = [seg_start_array(i,2),seg_end_array(i,2)];
                seg_z = [seg_start_array(i,3),seg_end_array(i,3)];
                p = plot3(seg_x,seg_y,seg_z,'Color',C);
                p.LineWidth = 3;
            end
            
        end
            
        axis('equal')
        xlabel('X');ylabel('Y');zlabel('Z');
        c = colorbar('Ticks',[0.5e17,1e17,1.5e17,2e17,2.5e17,3e17,3.5e17,4e17,4.5e17,5e17],...
         'TickLabels',{'0.5\times10^{17}','1\times10^{17}','1.5\times10^{17}','2\times10^{17}','2.5\times10^{17}','3\times10^{17}','3.5\times10^{17}','4\times10^{17}','4.5\times10^{17}','5\times10^{17}'});
        pos = c.Position;
        set(c,'Position',[pos(1)-0.07 pos(2) pos(3)*1.2 pos(4)])% To change size

        w = c.LineWidth;
        c.LineWidth = 2;
        c.FontSize = 70;

        set(gca,'FontSize',20)
        view(130,30)
        view(85,40)
        set(gca,'visible','off')

    end

end

output_matrix = [cellsize_matrix;gnd_matrix];
fid_output = fopen(strcat(case_str,num2str(number_seeds),'//',case_str,num2str(number_seeds),'_cellsize_gnd.txt'),'w');

for i = 1:2
    for j = 1:size(output_matrix,2)
        if i == 1
            formatSpec = '%8.4f \t';
            fprintf(fid_output,formatSpec,output_matrix(i,j));
        end
        if i == 2
            formatSpec = '%8.4e \t';
            fprintf(fid_output,formatSpec,output_matrix(i,j));
        end
    end
    fprintf(fid_output,'\n');
end



cellsize_matrix_list = [cellsize_matrix_list, cellsize_matrix]
gnd_matrix_list = [gnd_matrix_list, gnd_matrix]


%% output voro colors for 3d rendering
fid_output = fopen(strcat(case_str,num2str(number_seeds),'//',case_str,num2str(number_seeds),'_voro_color.txt'),'w');

formatSpec = '%8.6f,';
for i = 1:size(voro_gnd_pervol,2)
    fprintf(fid_output,formatSpec,voro_gnd_pervol(i)*1e20);
end

fprintf(fid_output,'\n');



end

