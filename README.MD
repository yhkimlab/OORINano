KAIST EE YHKLab (Prof. Yong-Hoon Kim)
OORINano in development: Dr. Joonho Park, Ms. SeungHeun Yu, 

###### Download and setting
### Download
git clone https://github.com/yhkimlab/OORINano oorinano

### setting
#set PYTHONPATH to the parent directory of oorinano or mv oorinano to any $PYTHONPATH such as .../site-packages; e.g.
export PYTHONPATH=$path_to_before_oorinano:$PYTHONPATH
#set PATH to oorinano/bin or mv oorinano/bin/* to your bin directory
export PATH=$path_to_oorinano/bin

#modify aux/env.py for pathway to simulator

###### Copy examples and refer to readme file in each example directory
### copy example directory to your_work_dir and run examples
cp -r Examples $work_dir

### or partly copy present working examples
cp -r Examples/catalysis  $your_work_dir
cp -r Examples/qtnegf     $your_work_dir

### Before run each example, check the vasp and siesta path in ~/utils/env.py




