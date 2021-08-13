# Complementary Dynamics cpp

Public c++ code release for the SIGGRAPH Asia 2020 paper [Complementary Dynamics](https://www.dgp.toronto.edu/projects/complementary-dynamics/) by Jiayi Eris Zhang, Seungbae Bang, David I.W. Levin and Alec Jacobson.

# Dependency

[libigl](https://github.com/libigl)\
[Bartels](https://github.com/dilevin/Bartels)

# Setup

Before running the code, you need to setup the following.\
Compile static library `libbartels` with Bartels:

    cd Bartles
    mkdir build
    cd build
    cmake .. -Dbartels_USE_STATIC_LIBRARY=ON
    make -j8

Compile static library `libigl` with libigl:

    cd libigl
    mkdir build
    cd build
    cmake .. -DLIBIGL_USE_STATIC_LIBRARY=ON
    make -j8

# Run

Run examples by `./complementary_dynamics example_model` (with `example_model` being one example in provided in `examples` folder). For instance, try running following command

    ./complementary_dynamics sphere

all parameter for physical model is defined in `examples\example.json`. try tweaking parameter and see how it effects the result.
