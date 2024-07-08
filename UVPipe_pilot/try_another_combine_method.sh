#!/bin/bash

rm curve_orbitwise_* source_* *orbit_* *animation* combined_* RA_Dec_combined_* UV_oriented_combined_VIS_for_*

#combinescript=create_combined_events_lists.py
combinescript=background_removed__create_combined_events_lists.py
#combinescript=bright_source__create_combined_events_lists.py

precombining=pre_combining.log
postcombining=post_combining.log


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
