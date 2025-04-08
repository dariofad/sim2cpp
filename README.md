# sim2cpp

This repository provides a tutorial on coverting a Simulink model to C++ code. 


## System Requirements

- Operating system: Windows/Linux/MacOS;

- Matlab (Simulink/Stateflow) version: >= 2020a. (Matlab license needed)

- [Embedded Coder Toolbox](https://www.mathworks.com/products/embedded-coder.html)

- C++ compiler (e.g., Clang for MacOS)


## Usage

Here, I use Matlab 2024b and a simple Simulink model `Simulink2Code.slx` to demonstrate how to convert this model to C++ code. 
The model `Simulink2Code.slx` implements a simple function that adds two input signals `x` and `y`, then outputs a signal `result` after addition.

### Install a C++ Compiler

For MacOS users, use the command `Clang --version` or `g++ --version` to check if the C++ compiler is installed.
For Windows/Linux users, please check it yourself, as I don't have a PC with Windows/Linux installed.

### Convert a Simulink Model to C++

Before converting `Simulink2Code.slx` to C++ code, it is neccessary to configure the Simulink model `Simulink2Code.slx`.

- Before C++ code generation, click `MODELING`/`Model Settings`;
  
    1. In the `Solver` Tab of `Configuration Parameters`, set the `Type` under the Sub-Tab `Solver selection` as `Fixed-step`.
       The users can set the `Fixed-step size (fundamental sample time)` under the Sub-Tab `Solver details` according to the needs of the simulation.
       Here, I set the `Fixed-step size (fundamental sample time)` as 1 (second);
       
    2. In the `Code Generation` Tab `Configuration Parameters`, set the `System target file` as `ert.tlc` and set the `Language` as `C++` under the Sub-Tab 'Target selection'.
       The users can select `Language standard` according to their needs. Here, I choose `C++14 (ISO)`;
       
    3. In the `Code Generation` Tab `Configuration Parameters`, check the item `Generate code only` and `Package code and artifacts` under the Sub-Tab `Build process`.
       The option `Generate code only` means we intend to perform the compliation later; The option `Package code and artifacts` means we
       
    4. In the `Code Generation` Tab `Configuration Parameters`, set the `La`
       
- Click `APPS` on the top of the navigation bar of the opened `Simulink2Code.slx`, go to `CODE GENERATION`/`Embedded Coder`;
  
- A new option `C++ CODE` will appear in the navigation bar of the opened `Simulink2Code.slx`;
  
- A new option ``

### Custom Main Function

`Simulink2Code.h` and `Simulink2Code.app` define the necessary header file and source file. 
We need to write our own main function to run the C++ code-based simulink model. An example `main.app` has been provided in this demo.

### Compile 

Run the following command to compile the C++ code. 
`g++ -o Simulink2Code main.cpp Simulink2Code.cpp -I. -std=c++14`

### Execute

Run the following command to exeute the generated executable file.
`./Simulink2Code`

## Related Links

- [Simulink code generation to c++ with Visual Studio build] (https://www.youtube.com/watch?v=xfLyc7BhoCk)

- [Simulink Tutorial - 21 - Code Generation From Model] (https://www.youtube.com/watch?v=PITAD2Mduw4)
  
- [Simulink/S-Function生成Windows可执行的C++代码] (https://www.bilibili.com/video/BV1v3411A7xX/?share_source=copy_web&vd_source=20273743de6c2b7275d51e34b7ccf476)






