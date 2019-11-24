conda env create --name sat-env pyton=3.7 -f env.yml

echo export PYTHONPATH=$PYTHONPATH:$PWD/code >> ~/.bashrc
echo export PYTHONPATH=$PYTHONPATH:$PWD/code >> ~/.zshrc
