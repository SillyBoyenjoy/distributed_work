# distributed_work
The distributed work includes master/client.

建议使用cplex12.7.0以上版本

# 主机程序：master/目录
  请在master目录下输入命令make run进行计算
  请在makedir/目录Makefile.tcpip文件中修改参数，包括lp文件、是否剪枝与从机ip端口等，部分参数说明：
  -iteration  定期从全局最优解从新启动改善全局最差解
  -cut        在gap达到3%后开始剪枝
  -polish     从某个解开始热启动
  parmipopt-transport-args修改连接从机的端口，ip地址与从机ip对应，端口号必须对应从机开启监听的端口

  lp文件请存放于data/目录中，并于Makefile.tcpip中修改，日志文件保存在log/目录中，求解结果保存在sol/目录中
  可随时中止程序运行并保存当前最优解，按下回车后输入quit即可

# 从机程序：client/目录
  请在client目录下输入命令make run打开监听，输入命令make clean关闭监听
  请在makedir/目录Makefile.tcpip文件中修改监听的端口号，请在9525~9528这个范围内选择端口
  client打开监听后自动在后台监听，待有主机发起计算后开始占用cpu与内存资源开始计算
