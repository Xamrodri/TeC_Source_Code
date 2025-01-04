%%% Calculate monthly max snow depth runoff fluxes for analysis at the catcment scale

% Where to store runoff figures:
dir_fig_runoff = [dir_fig '\Snow_cover'];
if ~exist(dir_fig_runoff, 'dir'); mkdir(dir_fig_runoff); end 

% Compute hydrological years starts and ends
period_start = min(find(month(date_seas) == 10 & day(date_seas) == 1,1), find(month(date_seas) == 11 & day(date_seas) == 1,1));
period_end = find(month(date_seas) == 9,1,'last');

Years_no_seas = unique(year(date_seas(period_start:period_end)));

Hydro_year_label = strcat(string(num2str(Years_no_seas(1:end-1)-100*floor(Years_no_seas(1:end-1)/100),'%02d ')), '/', ...
    string(num2str(Years_no_seas(2:end)-100*floor(Years_no_seas(2:end)/100),'%02d ')));

Month_seas = month(date_seas(period_start:period_end));

hydro_month = [10,11,12,1,2,3,4,5,6,7,8,9];
hydro_month_labels = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'};

[SND_yearly_ym] = deal(NaN(numel(Years_no_seas)-1, 12,size(DTM,1), size(DTM,2)));

month_m = date_m-calmonths(1);

for yy = 1:length(Years_no_seas)-1
    for mm = 1:12
      if hydro_month(mm) > 9
        ind_start_d = find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy), hydro_month(mm),15);
      else
        ind_start_d = find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy)+1, hydro_month(mm),15);
      end 

    SND_ym_map = nanmax(snow_depth(:,:,ind_start:ind_end),[],3); SND_ym_map(MASK ~=1) = NaN;
    SND_yearly_ym(yy,mm,:,:) = SND_ym_map;
    end
end 

figure
imagesc()
