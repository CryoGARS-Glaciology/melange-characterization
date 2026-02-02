function [vel_seas] = create_seasonal_speed_profiles(vel_dir,inland_idx,inland_vel,ref_adjust,vterm_idx,term_idx,Tyrs,Tmos,vdtmin,vdtmax,vfilter,years,seasons,zyrs,zmos)
%%% Extract seasonal speed profiles from transect of ITS_LIVE velocities

%INPUTS:
%vel_dir = velocity data directory
%inland_idx = index or the point along the profile where the referencing begins
%inland_vel = starting reference for velocity pts (code set-up that it grabs data for all points less than this reference... points move inland)
%ref_adjust = determines how much you are shifting the points in the seasonal profiles when converting to flow-following reference frame
%vterm_idx = annual terminus reference point
%vdtmin = minimum time separation for image pairs used to calculate velocity (days)
%vdtmax = maximum time separation for image pairs used to calculate velocity (days)
%vfilter = flag to determine whether to calculate seasonal speed profiles using all years of data ('all') or only years with elevation data ('annual')
%years = years that you want in your speed dataset
%seasons = months corresponding to each season (recommend: seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};)
%zyrs = years with elevation data (create vector with all years if you aren't actually considering elevations as well)
%zmos = months in each year with elevation data (create vector with all months that is equal to the years vector if you aren't actually considering elevations as well)

%OUTPUTS:
%vel_seas = average seasonal speed profile for each year (m/yr)

%loop through the pts
vel_pts = dir(vel_dir);
vel_seas = NaN(max(inland_idx),4,length(years));
for i = 1:length(vel_pts)
    if contains(vel_pts(i).name,'velocity')
        pt_ref = str2num(vel_pts(i).name(end-5:end-4));

        if pt_ref <= inland_vel %setting this <= and making rel_ref = inland_idx - pt_ref+1 gives the inland pt on the glacier

            %read the file
            V = readtable([vel_dir,vel_pts(i).name]);

            %filter out all the velocities based on temporal resolution
            short_dts = find(V.days_dt>vdtmin & V.days_dt<vdtmax); %get rid of all velocities with coarse temporal resolution
            vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
            vmos = month(vel_dates); vyrs = year(vel_dates);
            vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;

            %calculate seasonal average speeds for each year at the point
            if contains(vfilter,'all')

                %grab velocities for all years
                for p = 1:length(years)
                    yr_idx = find(vyrs == years(p));
                    if ~isempty(yr_idx)

                        %isolate the velocities for that year
                        for k = 1:length(yr_idx)
                            v_temp(k,1) = vels(yr_idx(k));
                        end

                        %calculate seasonal average
                        mos_yr = vmos(yr_idx);
                        for k = 1:4
                            if ~isempty(zmos(zyrs==years(p)))
                                %use the DEMs from that year to come up with the position wrt the terminus
                                rel_ref = vterm_idx(p)-pt_ref+ref_adjust; %edit if monthly terminus positions are available
                                if rel_ref >=1
                                    vel_seas(rel_ref,k,p) = nanmean(v_temp(ismember(mos_yr,seasons(k,:))==1));
                                end
                            end
                        end

                        clear vels_yr mos_yr v_temp;
                    end
                    clear yr_idx
                end

            elseif contains(vfilter,'annual')
                %only use velocities from years with DEMs
                for p = 1:length(years)
                    if sum(ismember(zyrs,years(p))) > 0 %only extract velocities from years with DEMs
                        yr_idx = find(vyrs == years(p));
                        if ~isempty(yr_idx)

                            %isolate the velocities for that year
                            for k = 1:length(yr_idx)
                                v_temp(k,1) = vels(yr_idx(k));
                            end

                            %calculate seasonal average if there is a DEM
                            %for that year and within that season
                            mos_yr = vmos(yr_idx);
                            for k = 1:4
                                if ~isempty(zmos(zyrs==years(p) & ismember(zmos,seasons(k,:))))
                                    %use the DEMs from that season to come
                                    %up with the position wrt the terminus
                                    rel_ref = round(median([inland_idx(zyrs==years(p) & ismember(zmos,seasons(k,:))), term_idx(Tyrs==years(p) & ismember(Tmos,seasons(k,:)))]))-pt_ref+ref_adjust;
                                    if rel_ref >=1
                                        vel_seas(rel_ref,k,p) = nanmean(v_temp(ismember(mos_yr,seasons(k,:))==1));
                                    end
                                end
                            end

                            clear vels_yr mos_yr v_temp;
                        end
                        clear yr_idx
                    end
                end
            else
                error('Unkown velocity sampling strategy! vfilter must be ''annual'' or ''all''')
            end

            clear V short_dts vel_dates vel_dts vmos vyrs vels rel_ref;
        end
        clear pt_ref;
    end
end




end