%Code for post processing
clc,clear;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Data Input Section %%%%%%%%%%
% i = Spanwise, j = Streamwise, k = Wall-normal 
    data_set = "example_data.mat";%AEM full field, containing U,V,W,zpos and xpos
    data.Reynolds_number = 3200; % Reynolds numbers of the data source
%%%% End of Data Input Section %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(data_set);
    data.zpos = zpos_delta;% Wall-normal grid postions
    correction_para = 0;% 加入修正量以匹配对数律
    U = U - correction_para;
    data.dx = mean(diff(xpos_delta)); % Streamwise data point spacing

%%%% Create Data Folder %%%%
    data_fold = ['data_',num2str(data.Reynolds_number)];
    currentFolder = pwd;
    data_fold_path = fullfile(currentFolder, data_fold);
    if ~exist(data_fold_path, 'dir')
        mkdir(data_fold_path);
        disp(['Folder created: ', data_fold_path]);
    else
        disp(['Folder already exists: ', data_fold_path]);
    end
%%%%%% End of Creating %%%%%%%

%%%%%%%% Calculating Mean Velocity Profile %%%%%%%%
    data.mean_U = squeeze(mean(mean(U,2),1));%save mean velocity profile
    figure;
        semilogx(data.zpos,data.mean_U,'-x',linewidth=2,DisplayName='Raw Data');
        hold on
        semilogx(data.zpos,log(data.zpos*data.Reynolds_number)/0.41+5,linewidth=1,DisplayName='Log-law');
        hold off
        xlabel('Wall-normal Distance $z^+$',Interpreter='latex');
        ylabel('Streamwise Velocity $U^+$',Interpreter='latex');
        legend(Location="northwest");
        filename = fullfile(data_fold_path,['Mean Velocity Profile of Re_tau = ',num2str(data.Reynolds_number),'.pdf']);
        saveas(gca,filename);
        close(gcf);
%%%% End of Calculating Mean Velocity Profile %%%%%


%%%%%%% Calculating Critical Droplet Size %%%%%%%%
    left_bound = 1e-5;%lower bound of initial droplet size
    right_bound = 1e-3;%upper bound of initial droplet size


    data.ratio = [0.01,0.02,0.05,0.10];
    for i = 1:length(data.ratio)
        [data.critical_value(i),data.physical_duration{i},data.physical_threshold{i},data.physical_tau{i},data.alternative_density{i}] = ...
            critical_value(U,V,data.dx,data.zpos,data.Reynolds_number,data.ratio(i),left_bound,right_bound);
    
        figure;
            semilogx(data.physical_duration{i},data.physical_threshold{i},'-s',Color='red');
            hold on 
            semilogx(data.physical_duration{i},data.physical_tau{i},'-x',Color='blue');
            t = ['$d/\delta=$ ',num2str(data.critical_value(i),'%6.5e')];
            title(t,'interpreter','latex')
            hold off
            figure_name = ['Critical Droplet Size of Re_tau = ',num2str(data.Reynolds_number),' ( Ratio = ',num2str(data.ratio(i),'%3.2f'),')','.pdf'];
            filename = fullfile(data_fold_path,figure_name);
            exportgraphics(gcf, filename, 'ContentType', 'vector');
            close(gcf);
    end
    
%%%% End of Calculating Critical Droplet Size %%%%

%%

%%%%%%%% Save Data %%%%%%%%%

result_data_name = ['Re_tau = ',num2str(data.Reynolds_number),'.mat'];
filename = fullfile(data_fold_path,result_data_name);
save(filename,"data");

%%%% End of Saving Data %%%%