// #include <cstdio>
// #include "mpcACCsystem.h"
// 
// static mpcACCsystem model; // 实例化模型对象
// 
// int main(int argc, const char *argv[]) {
//   (void)(argc);
//   (void)(argv);
// 
//   model.initialize(); // 初始化模型
// 
//   double Ts = 0.1; // 假设采样周期为0.1s
//   int steps = 100;
// 
//   for (int i = 0; i < steps; ++i) {
//     double t = i * Ts;
// 
//     // 执行一步仿真
//     model.step();
// 
//     // 打印内部变量（假设你想看的是优化器输出 u 和正弦输入 a_lead）
//     // 这些变量都在 mpcACCsystem_B 结构体中
//     printf("t = %.2f, a_lead = %.3f, control_u = %.3f\n",
//        t,
//        model.getALead(),
//        model.getControlU());
//   }
// 
//   model.terminate(); // 清理资源
//   return 0;
// }

#include <iostream>
#include "mpcACCsystem.h"

int main() {
    mpcACCsystem model;
    model.initialize();

    const double Ts = 0.1;
    const int steps = 801;

    for (int i = 0; i < steps; ++i) {
        model.step();

        // 访问外部输出信号
        const auto& y = model.getExternalOutputs();

        printf("t = %.2f, a_lead = %.3f, d_rel = %.3f, v_rel = %.3f, v_ego = %.3f, a_ego = %.3f\n",
               i * Ts,
               y.a_lead,
               y.d_rel,
               y.v_rel,
               y.v_ego,
               y.a_ego);
    }

    model.terminate();
    return 0;
}