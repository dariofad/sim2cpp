# sim2cpp

This repository provides a tutorial on coverting a Simulink model to C++ code. 


## System Requirements

- Operating system: Windows/Linux/MacOS;

- Matlab (Simulink/Stateflow) version: >= 2020a. (Matlab license needed)

- [Embedded Coder Toolbox](https://www.mathworks.com/products/embedded-coder.html)

- C++ compiler (e.g., Clang for MacOS; GCC for Windows)


## Covert a Simulink Model to C++

Here, we will use a simple Simulink model `Simulink2Code.slx` to demostrate how to convert this model to C++ code.

## Custom Main Function

`Simulink2Code.h` and `Simulink2Code.app` define the necessary header file and source file. 
We need to write our own main function to run the C++ code-based simulink model. An example `main.app` has been provided in this demo.

## Compile 

Run the following command to compile the C++ code. 

`g++ -o Simulink2Code main.cpp Simulink2Code.cpp -I. -std=c++14`


## Execute

Run the following command to exeute the generated executable file.

`./Simulink2Code`


## Related Links

- [Simulink code generation to c++ with Visual Studio build] (https://www.youtube.com/watch?v=xfLyc7BhoCk)

- [Simulink Tutorial - 21 - Code Generation From Model] (https://www.youtube.com/watch?v=PITAD2Mduw4)
  
- [Simulink/S-Function生成Windows可执行的C++代码] (https://www.bilibili.com/video/BV1v3411A7xX/?share_source=copy_web&vd_source=20273743de6c2b7275d51e34b7ccf476)






