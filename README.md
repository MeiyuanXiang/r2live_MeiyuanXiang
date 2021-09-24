# r2live_MeiyuanXiang
r2live相关论文、代码中文注释以及代码改动

# 参考
https://github.com/hku-mars/r2live 

# 环境
1. Ubuntu（测试了Ubuntu16.04.5、Ubuntu18.04）  
2. ROS (测试了kinetic、melodic)  
3. Ceres Solver（测试ceres-solver-1.14.0）  
4. livox_ros_driver  

# 编译
1. 下载源码 git clone https://github.com/MeiyuanXiang/r2live_MeiyuanXiang.git  
2. 将r2live_MeiyuanXiang\src下的r2live拷贝到ros工程空间src文件夹内，例如~/catkin_ws/src/  
3. cd ~/catkin_ws  
4. catkin_make  
5. source ~/catkin_ws/devel/setup.bash  

# Bag数据
https://drive.google.com/drive/folders/1LpoX6_05Zic-mRLOD38EO0w2ABXI1rrW 
harbor.bag、hku_main_building.bag和indoor_aggressive.bag  

# 运行
roslaunch r2live demo.launch  
rosbag play harbor.bag  
