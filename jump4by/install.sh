#!/bin/bash 

# the execute script % 
JUMP2PATH=~/.jump2
if [ ! -d "$JUMP2PATH" ]; then
mkdir -p $JUMP2PATH/bin
cp bin/jump2 $JUMP2PATH/bin
mkdir -p $JUMP2PATH/env
cp bin/cluster_config.example $JUMP2PATH/env/cluster_config.sh
mkdir -p $JUMP2PATH/utils
cp bin/set_parellel_env.sh $JUMP2PATH/utils
fi

# uninstall the old version % 
echo y | pip uninstall jump2 
rm -rf build 

export PATH=~/.jump2/bin:$PATH
# install the laster version % 
python setup.py build 
python setup.py install --user 
