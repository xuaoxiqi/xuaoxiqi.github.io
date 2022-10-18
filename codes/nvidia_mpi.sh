#/bin/bash

MYAPP=/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/comm_libs/mpi

export PATH=${MYAPP}/bin:$HOME/.local/bin:$HOME/bin:$PATH
export LD_LIBRARY_PATH=${MYAPP}/lib:$LD_LIBRARY_PATH
export INCLUDE=${MYAPP}/include:$INCLUDE
export CPATH=${MYAPP}/include:$CPATH
export MANPATH=${MTAPP}/share/man:$MANPATH
