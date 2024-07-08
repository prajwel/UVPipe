#!/bin/bash


combinescript=create_combined_events_lists.py
#combinescript=background_removed__create_combined_events_lists.py
#combinescript=bright_source__create_combined_events_lists.py

precombining=pre_combining.log
postcombining=post_combining.log
extralog=extra.log

./check_CALDB_configuration.py LEVL1AS1UVT*tar_V*
if [ $? -eq 1 ]; then
    exit
fi

./groom.py > $extralog 2>&1
for d in uvt_*; do ./cascade_sieve.py "$d"; done >> $extralog 2>&1

echo "Running modify_episode_events_lists.py" > $precombining
./modify_episode_events_lists.py >> $precombining 2>&1

echo "" >> $precombining
echo "Running combine_episodes_from_same_orbit.py" >> $precombining
./combine_episodes_from_same_orbit.py >> $precombining 2>&1

./L2_auto_report.py >> $extralog 2>&1

echo "" >> $precombining
echo "Running prune_L2.py" >> $precombining
./prune_L2.py >> $precombining 2>&1

echo "" >> $precombining
echo "Running create_coadded_VIS_images.py" >> $precombining
./create_coadded_VIS_images.py >> $precombining 2>&1

./to_check_UL2P_run.py >> $extralog 2>&1

echo "" >> $precombining
echo "Running update_coadded_VIS_image_headers.py" >> $precombining
./update_coadded_VIS_image_headers.py >> $precombining 2>&1

echo "Running $combinescript" > $postcombining
./$combinescript >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_exposure_maps.py" >> $postcombining
./create_combined_exposure_maps.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_products_in_detector_coordinates.py" >> $postcombining
./create_combined_products_in_detector_coordinates.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_VIS_images.py" >> $postcombining
./create_combined_VIS_images.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running UV_source_detection_for_AstrometryNet.py" >> $postcombining
./UV_source_detection_for_AstrometryNet.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running VIS_source_detection_for_AstrometryNet.py" >> $postcombining
./VIS_source_detection_for_AstrometryNet.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running modify_VIS_sources_to_UV_orientation.py" >> $postcombining
./modify_VIS_sources_to_UV_orientation.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running run_AstrometryNet_for_VIS.py" >> $postcombining
./run_AstrometryNet_for_VIS.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running run_AstrometryNet_for_UV.py" >> $postcombining
./run_AstrometryNet_for_UV.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running fixed_plate_scale_WCS_solution.py" >> $postcombining
./fixed_plate_scale_WCS_solution.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running apply_WCS_to_combined_products.py" >> $postcombining
./apply_WCS_to_combined_products.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_RA_Dec_combined_exposure_maps.py" >> $postcombining
./create_RA_Dec_combined_exposure_maps.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_combined_products_in_RA_Dec_coordinates.py" >> $postcombining
./create_combined_products_in_RA_Dec_coordinates.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running create_RA_Dec_combined_VIS_images.py" >> $postcombining
./create_RA_Dec_combined_VIS_images.py >> $postcombining 2>&1

echo "" >> $postcombining
echo "Running find_bad_combined_events_lists_using_curvit.py" >> $postcombining
./find_bad_combined_events_lists_using_curvit.py >> $postcombining 2>&1

cat $precombining $postcombining > combining.log

rm -rf output
