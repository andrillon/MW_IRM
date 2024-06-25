function ROI_data_mars = extract_roi_Marsbar(roifilepath,des_path)

marsbar('on');

% MARSBAR
rois = maroi('load_cell', roifilepath); % make maroi ROI objects
des = mardo(des_path);  % make mardo design object
mY = get_marsy(rois{:}, des, 'mean'); % extract data into marsy data object
ROI_data_mars  = summary_data(mY);  % get summary time course(s)

end


