# sim2cpp

This repository provides a tutorial on how to convert a Simulink model into C++ code. 


## System Requirements

- Operating system: Windows/Linux/MacOS;

- Matlab (Simulink/Stateflow) version: >= 2020a. (Matlab license needed)

    - [Embedded Coder Toolbox](https://www.mathworks.com/products/embedded-coder.html)

- C++ compiler (e.g., Clang for MacOS)


## The Steps to Convert a Simulink Model to C++ Code

Here, I use Matlab (Simulink/Stateflow) 2024b and a toy Simulink model `Simulink2Code.slx` to demonstrate how to convert a Simulink model into C++ code. 
The `Simulink2Code` model implements a simple function: given two input signals `x` and `y`, it adds them and outputs a signal `result` after addition.

### 1. Install a C++ Compiler

- For MacOS users, use the command `Clang --version` or `g++ --version` to check if the C++ compiler is installed.
  If not, use the command `xcode-select --install` to install Xcode Command Line Tools.
- For Windows/Linux users, please check it yourself, as I don't have a Windows/Linux development environment.

### 2. Convert a Simulink Model to C++

Before converting the model `Simulink2Code.slx` to C++ code, it is neccessary to configure `Simulink2Code.slx`.

- Click `MODELING`/`Model Settings`;
  
    1. In the `Solver` Tab of `Configuration Parameters`, set the `Type` under the Sub-Tab `Solver selection` as `Fixed-step`.
       The users can set the `Fixed-step size (fundamental sample time)` under the Sub-Tab `Solver details` according to the needs of the simulation.
       Here, I set the `Fixed-step size (fundamental sample time)` as 1 (second);

    2. In the `Hardware Implementation` Tab of `Configuration Parameters`, set the `Device vendor` and `Device type` according to the user's own development environment.
       Here, I choose `Apple` and `ARM64`;
       
    3. In the `Code Generation` Tab of `Configuration Parameters`, set the `System target file` as `ert.tlc` and set the `Language` as `C++` under the Sub-Tab 'Target selection'.
       The users can select `Language standard` according to their needs. Here, I choose `C++14 (ISO)`;
       
    4. In the `Code Generation` Tab of `Configuration Parameters`, check the item `Generate code only` and `Package code and artifacts` under the Sub-Tab `Build process`.
       The option `Generate code only` means we intend to perform the compliation later; The option `Package code and artifacts` allows the generation of a zip file with a user-specified name,
       which contains all of the neccessary files to run the C++ code. The users can specify the name of this zip file under the Sub-Tab `Advanced parameters`.
       
    5. In the `Code Generation` Tab of `Configuration Parameters`, check the option `Create code generation report` and `Open report automatically` under the Sub-Tab `Report` of the Tab `Code Generation`.
       Once the C++ code is generated, the user can view the C++ code generation report directly.

    6. Click 'OK' at the bottom of `Configuration Parameters` to save the user's settings.
       
- Click `APPS` on the top of the navigation bar of the opened `Simulink2Code.slx`, go to `CODE GENERATION`/`Embedded Coder`;
  A new option `C++ CODE` will appear in the navigation bar of the opened `Simulink2Code.slx`; 
  
- The users can customize the Class Name and Namespace in the button `Code Interface`. The default `Class Name` is the name of this simulink model `Simulink2Code`.

- Click `Generate Code` button, choose `Generate code only`. 

- **Finally, a folder named `Simulink2Code_ert_rtw` and a zip file named `Simulink2Code` will be generated in the current directory. And a `Generation Report` window will pop up**.

### 3. Custom a Main Function

Now, navigate to the `Simulink2Code_ert_rtw` folder. `Simulink2Code.h` and `Simulink2Code.app` define the necessary header file and source file. 
We need to write a main function to run the C++ code-based simulink model. An example `main.app` has been provided in this toy demo.

### 4. Compile 

Run the following command to compile the C++ code. 
`g++ -o Simulink2Code main.cpp Simulink2Code.cpp -I. -std=c++14`

### 5. Execute

Run the following command to exeute the generated executable file.
`./Simulink2Code`

## An Instance of automatic_transmission

I also attached an automatic transmission model, a commonly-used simulink model in falsification community.
**I have made the necessary file arrangements and modifications.**
Please navigate to the folder `AT/automatic_transmission_ert_rtw` and run the following command to compile the C++ code.
```
g++ -o automatic_transmission main.cpp automatic_transmission.cpp automatic_transmission_data.cpp  -I. -std=c++14
./automatic_transmission
```

## An Instance of Adaptive Cruise Control System 

I attached an Adaptive Cruise Control (ACC) system using Model Predictive Control (MPC) from Simulink.
The `mpcACCsystem` folder contains official configuration files and the original reference model. 
`mpcACCsystem.slx` is used to generate C++ code, and the resulting C++ files are placed in the `mpcACCsystem_ert_rtw` directory.

**I have made the necessary file arrangements and modifications.**
Please navigate to the folder `AdaptiveCruiseControlExample/mpcACCsystem_ert_rtw` and run the following command to compile the C++ code.
```
g++ -o mpcACCsystem main.cpp mpcACCsystem.cpp mpcACCsystem_data.cpp rtGetNaN.cpp rt_nonfinite.cpp -I. -std=c++14
./mpcACCsystem
```
The output is printed in a time-stamped format showing key system signals at each simulation step, including a_lead, d_rel, v_rel, v_ego, and a_ego.
```
t = 0.00, a_lead = 0.000, d_rel = 40.375, v_rel = 4.994, v_ego = 20.006, a_ego = 0.852
t = 0.10, ...
```

In the MATLAB command window, you can also run the Simulink model directly using the command 
```
sim(mpcACCsystem)
```

The simulation results will be stored in the `logsout` variable in the workspace.


## Related Links

- [Simulink code generation to c++ with Visual Studio build] (https://www.youtube.com/watch?v=xfLyc7BhoCk)

- [Simulink Tutorial - 21 - Code Generation From Model] (https://www.youtube.com/watch?v=PITAD2Mduw4)
  
- [Simulink/S-Function生成Windows可执行的C++代码] (https://www.bilibili.com/video/BV1v3411A7xX/?share_source=copy_web&vd_source=20273743de6c2b7275d51e34b7ccf476)

- [Adaptive Cruise Control System Using Model Predictive Control from Simulink] (https://www.mathworks.com/help/mpc/ug/adaptive-cruise-control-using-model-predictive-controller.html)



