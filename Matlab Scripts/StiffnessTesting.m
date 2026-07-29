%% Alginate Stiffness Calculation

clear all; close all; %clc

root = "/Users/averywilliamson/Library/CloudStorage/GoogleDrive-averyw@andrew.cmu.edu/Shared drives/CMU BORG - Hydrogel Printing and Biodegradable Robots/ALGINATE CHARACTERIZATION/MTS Experiments/Ring/";

inventory_file = "Inventory 2.xlsx";

inventory = readtable(root+inventory_file,"ReadVariableNames",true);

folder_b = root + "BUOYANCY/"; % read buoyancy data
files_b = {dir(folder_b).name};
csvfile_b = files_b{3};
data_b = readmatrix(folder_b + csvfile_b,"NumHeaderLine",8);
t_b = data_b(1:8703,1);
L_b = data_b(1:8703,2);
F_b = data_b(1:8703,3);

for sample = 21:(length(inventory.Folder))

    folder = root + inventory.Folder(sample) + "/";  % access sample folder
    files = {dir(folder).name};                 % read files in folder
    csvfile = files{3};                         % select csv
    if csvfile == ".DS_Store"
        csvfile = files{4};
    end
    sample_name = inventory.SampleID(sample);        % record sample ID
    sample_name = sample_name{1};               % isolate str


    %% Sample Parameters
    height = inventory.height_mm_(sample); % [mm]
   
    
    OD = inventory.finalOD_px_(sample) / inventory.Scale_pxpmm_(sample);
    ID = inventory.finalID_px_(sample) / inventory.Scale_pxpmm_(sample);
    
    A = (pi()*(OD/2)^2) - (pi()*(ID/2)^2);

    %% Read Data
    
    data = readmatrix(folder + csvfile,"NumHeaderLine",8);
    t = data(1:8703,1);
    L = data(1:8703,2);
    F = data(1:8703,3);
    F = F - F_b;
    
    figure("Position",[100,100,1800,300],"Color","w");
    subplot(1,5,1);
    plot(L,F);
    title('Raw Curves')
    xlabel("Displacement [mm]")
    ylabel('Force [N]')
    
    % find the starting points of the different cycles
    [pks,locs] = findpeaks(L);
    locs = [1;locs;length(L)];
    
    % separating out the last five cycles
    up_curves = cell(5,1);
    down_curves = cell(5,1);
    subplot(1,5,2); hold all
    minL = 0; % used to ensure the L range are the same for all curves
    maxL = 0;
    
    for i=1:5
        
        ind = (locs(end-5 + (i-1))-150) : (locs(end-5 + i)-150); % indices of the current cycle from L = 0 to L=0
        
    
        % data for current cycle
        t_ = t(ind) - t(ind(1)); % time from t0 of cycle
        L_ = L(ind);
        F_ = F(ind);
        
        % subtracting off offset in load cell and finding when the force goes above 0
        noise_floor_mean = mean(F_(t_ < 1)); % getting mean force in first 10 secs (no contact)
        noise_floor_std = std(F_(t_ < 1));
        F_ = F_ - noise_floor_mean; % resetting 0 point
    
        contact = find(F_ > 8*noise_floor_std); % where the sample is in contact
        peak = find(L_ == max(L_)); % location of max force (differentiating loading and unloading curves)
        
        % plotting all force data and the force data for the contact region
        plot(t_,F_,'--','Color',[1,1,1]*0.2);
        plot(t_(contact),F_(contact),'.r')
    
        % splitting out the loading and unloading curves
        up_curves{i} = {L_(contact(1):peak),F_(contact(1):peak)};
        down_curves{i} = {L_(peak:contact(end)),F_(peak:contact(end))};
    
        % updating the range of displacements where all curves have positive
        % force
        minL = max([minL, min(L_(contact(1):peak))]);
        maxL = max([maxL, min(L_(peak:contact(end)))]);
        
    end
    
    ylabel("Force [N]")
    xlabel("Time [s]")
    title({"Isolation Contact Regions","(Red = Contact Established)"})
    
    % Truncate curves to the same L range, plot 5 upcurves
    subplot(1,5,3); hold all
    for i=1:5
        ind = up_curves{i}{1}(:,1) > minL;
        up_curves{i}{1} = up_curves{i}{1}(ind,1);
        up_curves{i}{2} = up_curves{i}{2}(ind,1);
        
        plot(up_curves{i}{1}(:,1),up_curves{i}{2}(:,1))
    end
    ylabel("Force [N]")
    xlabel("Displacement [mm]")
    title("Loading Curves")
    
    
    %% Calculate Stress and Strain
    A0 = A; % original area
    L0_offset = minL; % how much to adjust the rest length by to account for plastic deformation
    
    L0 = height - L0_offset;
    lam0 = 1 - L0_offset / height;
    A = A0 / lam0; % adjusted area due to plastic deformation
    
    N = 250; % number of interpolation points
    L = linspace(1.4, 1.49,N); % using the same linear region for all curves
    
    P = zeros(5,N); % engineering stress [kPa] (1000*force [N] / Area [mm^2])
    for i=1:5
        P(i,:) = interp1(up_curves{i}{1}(:,1),up_curves{i}{2}(:,1),L,"linear","extrap") / A;
                 % (all L vals in linear, all F vals in linear, over
                 % linspace L)
    end
    % L = L - L0_offset + 0.001; % accounting for offset % 0.001 accounts for zero point not have a sign 
    S = (L - L0_offset)/ L0; % strain
    
    % Getting average and standard deviation
    meanP = mean(P,1);
    stdP = std(P,[],1);
    
    err_P = [meanP + stdP, fliplr(meanP - stdP)];
    % err_L = [lam, fliplr(lam)];
    
    subplot(1,5,4);hold all
    plot(S, P, "-");
    % linear_approx = polyfit(S,meanP,1);
    slope = (meanP(end)-meanP(1))/(S(end) - S(1));
    E = slope*S;
    offset = E(1) - meanP(1);
    E = E - offset;
    plot(S, E, 'Color', 'r', 'LineWidth', 3);
   
    % fill(err_L,err_P,'k',"EdgeColor","none","HandleVisibility","off","FaceAlpha",0.2);
    ylabel("Effective Stress [N/mm^2]")
    xlabel("Strain [mm/mm]")
    title({"Stress vs Strain", "(Red = Average Linear Approx)"})
    str = {"Eeff = " + slope + "MPa"};
    text(S(200), E(100), str)

    %% save files and values to inventory
    saveas(gcf,folder + "/figs/Stress-Strain Curves.svg",'svg'); % save svgs to drive folder
    str = "P" + (sample+1);
    writematrix(slope,root + inventory_file,'Sheet',1,'Range', str); % add E to inventory sheet
    
    

end



