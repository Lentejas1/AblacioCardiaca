cmake_minimum_required(VERSION 3.23)
project(AblacioCardiaca)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)

add_executable(AblacioCardiaca
        .idea/.gitignore
        .idea/AblacioCardiaca.iml
        .idea/modules.xml
        .idea/vcs.xml
        .idea/workspace.xml
        AblacioCardiaca.cpp
        GraficosAblacion.py)
