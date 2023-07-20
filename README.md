# NoahMP-Sensitivity-NCP

## 《华北平原农田生态系统的Noah-MP模型参数敏感性分析》实现代码

参考[grey-nearing/NoahMP-Sensitivity (github.com)](https://github.com/grey-nearing/NoahMP-Sensitivity)对模型代码及敏感性分析代码进行完善
修改处均用 `wangxi`注释

## 说明

1. `POINT_NOAHMP`文件夹内为模型代码，可以使用`makefile`进行编译或在Visual Studio内编译。
2. `YC`与`LC`文件夹内为禹城站与栾城站的敏感性分析代码，用编译好的模型替换`main.exe`后，运行`Sobol_Wrapper.m`即可。
3. `YC`与`LC`文件夹中的`restart_original.txt`为多轮预热后的重启文件，若需预热要给定土壤水初始场`soil_init.txt`。
4. `YC`与`LC`文件夹中的`obs.txt`并非两站点的真实观测，而是由模型默认参数运行后得到的结果。
   
### 如有问题可联系作者（1553638571@qq.com）