mkdir ./../../../MWA_Beam/
git clone https://github.com/MWATelescope/mwa_pb.git ./../../../MWA_Beam/
wget http://cerberus.mwa128t.org/mwa_full_embedded_element_pattern.h5 -P ./../../../MWA_Beam/mwa_pb/data
cd ./../../../MWA_Beam/
python setup.py install
cd ./../mwa-satellites/code/tile_analysis

