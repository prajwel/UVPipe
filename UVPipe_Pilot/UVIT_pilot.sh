#!/bin/bash


#combinescript=create_combined_events_lists.py
combinescript=background_removed__create_combined_events_lists.py
#combinescript=bright_source__create_combined_events_lists.py

precombining=pre_combining.log
postcombining=post_combining.log
extralog=extra.log

python check_CALDB_configuration.py LEVL1AS1UVT*tar_V*
if [ $? -eq 1 ]; then
    exit
fi

python groom.py > $extralog 2>&1
for d in uvt_*; do python cascade_sieve.py "$d"; done >> $extralog 2>&1

echo "Running modify_episode_events_lists.py" > $precombining
python modify_episode_events_lists.py >> $precombining 2>&1

echo "" >> $precombining
echo "Running combine_episodes_from_same_orbit.py" >> $precombining
python combine_episodes_from_same_orbit.py >> $precombining 2>&1

python L2_auto_report.py >> $extralog 2>&1

echo "" >> $precombining
echo "Running prune_L2.py" >> $precombining
python prune_L2.py >> $precombining 2>&1

echo "" >> $precombining
echo "Running create_coadded_VIS_images.py" >> $precombining
python create_coadded_VIS_images.py >> $precombining 2>&1

python to_check_UL2P_run.py >> $extralog 2>&1

echo "" >> $precombining
echo "Running update_coadded_VIS_image_headers.py" >> $precombining
python update_coadded_VIS_image_headers.py >> $precombining 2>&1

echo "Running $combinescript" > $postcombining
python $combinescript >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_exposure_maps.py" >> $postcombining
python create_combined_exposure_maps.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_products_in_detector_coordinates.py" >> $postcombining
python create_combined_products_in_detector_coordinates.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_VIS_images.py" >> $postcombining
python create_combined_VIS_images.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running UV_source_detection_for_AstrometryNet.py" >> $postcombining
python UV_source_detection_for_AstrometryNet.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running VIS_source_detection_for_AstrometryNet.py" >> $postcombining
python VIS_source_detection_for_AstrometryNet.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running modify_VIS_sources_to_UV_orientation.py" >> $postcombining
python modify_VIS_sources_to_UV_orientation.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running run_AstrometryNet_for_VIS.py" >> $postcombining
python run_AstrometryNet_for_VIS.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running run_AstrometryNet_for_UV.py" >> $postcombining
python run_AstrometryNet_for_UV.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running fixed_plate_scale_WCS_solution.py" >> $postcombining
python fixed_plate_scale_WCS_solution.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running apply_WCS_to_combined_products.py" >> $postcombining
python apply_WCS_to_combined_products.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_RA_Dec_combined_exposure_maps.py" >> $postcombining
python create_RA_Dec_combined_exposure_maps.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_products_in_RA_Dec_coordinates.py" >> $postcombining
python create_combined_products_in_RA_Dec_coordinates.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_RA_Dec_combined_VIS_images.py" >> $postcombining
python create_RA_Dec_combined_VIS_images.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running find_bad_combined_events_lists_using_curvit.py" >> $postcombining
python find_bad_combined_events_lists_using_curvit.py >> $postcombining 2>&1

cat $precombining $postcombining > combining.log

rm -rf output
