# Create local conda env installation
conda env create python=3.7 -f env.yml --prefix ./beam-env
conda config --set env_prompt '({name})'

# Activate beam_env
conda activate ./beam-env

# install mwa primary beam repo and download FEE data
mkdir ./beam-env/lib/python3.7/site-packages/MWA_PB/
git clone https://github.com/MWATelescope/mwa_pb.git ./beam-env/lib/python3.7/site-packages/MWA_PB/
wget http://cerberus.mwa128t.org/mwa_full_embedded_element_pattern.h5 -P ./beam-env/lib/python3.7/site-packages/MWA_PB/mwa_pb/data
cd ./beam-env/lib/python3.7/site-packages/MWA_PB/
python setup.py install
conda deactivate
cd ../../../../../

