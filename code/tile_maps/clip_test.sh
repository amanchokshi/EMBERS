echo "-40dBm test"
p tile_maps.py --start_date=2019-09-12 --stop_date=2020-03-16 --out_dir=./../../outputs/tile_maps/40dBm/tile_maps_raw/ --rfe_clip=-40
p tile_maps_norm.py --out_dir=./../../outputs/tile_maps/40dBm/tile_maps_norm/ --map_dir=./../../outputs/tile_maps/40dBm/tile_maps_raw
p plot_tile_maps.py --out_dir=./../../outputs/tile_maps/40dBm --map_dir=./../../outputs/tile_maps/40dBm/tile_maps_norm
p compare_beams.py --out_dir=./../../outputs/tile_maps/40dBm/compare_beams/ --map_dir=../../outputs/tile_maps/40dBm/tile_maps_norm/

echo "-30dBm test"
p tile_maps.py --start_date=2019-09-12 --stop_date=2020-03-16 --out_dir=./../../outputs/tile_maps/30dBm/tile_maps_raw/ --rfe_clip=-30
p tile_maps_norm.py --out_dir=./../../outputs/tile_maps/30dBm/tile_maps_norm/ --map_dir=./../../outputs/tile_maps/20dBm/tile_maps_raw
p plot_tile_maps.py --out_dir=./../../outputs/tile_maps/30dBm --map_dir=./../../outputs/tile_maps/30dBm/tile_maps_norm
p compare_beams.py --out_dir=./../../outputs/tile_maps/30dBm/compare_beams/ --map_dir=../../outputs/tile_maps/30dBm/tile_maps_norm/

echo "-20dBm test"
p tile_maps.py --start_date=2019-09-12 --stop_date=2020-03-16 --out_dir=./../../outputs/tile_maps/20dBm/tile_maps_raw/ --rfe_clip=-20
p tile_maps_norm.py --out_dir=./../../outputs/tile_maps/20dBm/tile_maps_norm/ --map_dir=./../../outputs/tile_maps/20dBm/tile_maps_raw
p plot_tile_maps.py --out_dir=./../../outputs/tile_maps/20dBm --map_dir=./../../outputs/tile_maps/20dBm/tile_maps_norm
p compare_beams.py --out_dir=./../../outputs/tile_maps/20dBm/compare_beams/ --map_dir=../../outputs/tile_maps/20dBm/tile_maps_norm/


