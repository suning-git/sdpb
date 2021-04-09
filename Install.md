This guide is for building SDPD.

* [Installation](#installation)

# Installation

1. sdpd is based on sdpb 2.4.0. Try to built sdpd 2.4.0 first and all its required libraries.

2. On the cluster of Perimeter Institute, I use following commands to build sdpd

        module load openmpi/gcc-9/64/4.0.5
        module load cmake/3.18.1
        module load boost/1.67.0
        module load python/3.7
        module load trilinos/13.0.0
        module load openblas
        
        git clone https://github.com/suning1985/sdpb.git
        cd sdpb
        git checkout sdpd_experimental
        
        ./waf configure --prefix /home/nsu2/packages/install/sdpd --elemental-dir=/cm/shared/apps/elemental/sdpb-fork/ --boost-dir $BOOST_ROOT --rapidjson-incdir=/cm/shared/apps/rapidjson/1.1.0/include
        ./waf install

    The commands are essentially the same from sdpb 2.4.0, except downloading. A binary named sdpd will appear in the folder specified by --prefix.

On the Caltech cluster, I use following commands to build sdpd

        module load cmake/3.10.2 gcc/7.3.0 openmpi/3.0.0 boost/1_68_0-gcc730
        
        git clone https://gitlab.com/bootstrapcollaboration/elemental.git
        mkdir -p elemental/build
        cd elemental/build
        export CXX=mpicxx
        export CC=mpicc
        cmake .. -DCMAKE_INSTALL_PREFIX=/central/groups/dssimmon/ning/packages/install
        make && make install
        
        cd ../..
        git clone https://github.com/Tencent/rapidjson.git

        git clone https://github.com/suning1985/sdpb.git
        cd sdpb
        git checkout sdpdd_experimental
        
        ./waf configure --prefix /central/groups/dssimmon/ning/sdpdd_experimental/install --elemental-dir=/central/groups/dssimmon/ning/packages/install --rapidjson-incdir=/central/groups/dssimmon/ning/packages/rapidjson/include
        ./waf install
